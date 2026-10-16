//! Expansion of typed matrix DAGs without evaluating held matrix operations.
use super::{BlockGrid, MatrixExpr, MatrixKind, MatrixOperand, SliceSpec};
use crate::expand::{expand_with_options, ExpandBudget, ExpandOptions};
use crate::expr::Expr;
use crate::matrix::{Matrix, MatrixError, MResult};
use rustc_hash::FxHashMap;

impl MatrixExpr {
    pub fn expand(&self) -> MResult<Self> {
        self.expand_with_options(&ExpandOptions::default(), &mut ExpandBudget::new(100_000))
    }

    /// Apply the scalar expansion hints and distribute ordered matrix products.
    /// Matrix powers and held unary operations retain their mathematical node.
    pub fn expand_with_options(&self, options: &ExpandOptions, budget: &mut ExpandBudget) -> MResult<Self> {
        Context { options, budget, matrices: FxHashMap::default(), scalars: FxHashMap::default() }.matrix(self)
    }
}

struct Context<'a> {
    options: &'a ExpandOptions,
    budget: &'a mut ExpandBudget,
    // Owners prevent allocation identities from being recycled during this
    // traversal. Semantic equality may identify different bound-variable names.
    matrices: FxHashMap<usize, (MatrixExpr, MatrixExpr)>,
    scalars: FxHashMap<usize, (Expr, Expr)>,
}

fn failure(error: impl std::fmt::Display) -> MatrixError { MatrixError::Value(error.to_string()) }

impl Context<'_> {
    fn scalar(&mut self, expression: &Expr) -> MResult<Expr> {
        if let Some((_, value)) = self.scalars.get(&expression.dag_id()) { return Ok(value.clone()); }
        let value = expand_with_options(expression, self.options, self.budget).map_err(failure)?;
        self.scalars.insert(expression.dag_id(), (expression.clone(), value.clone()));
        Ok(value)
    }

    fn scalar_child(&mut self, expression: &Expr) -> MResult<Expr> {
        if self.options.deep { self.scalar(expression) } else { Ok(expression.clone()) }
    }

    fn matrix_child(&mut self, expression: &MatrixExpr) -> MResult<MatrixExpr> {
        if self.options.deep { self.matrix(expression) } else { Ok(expression.clone()) }
    }

    fn operand(&mut self, operand: &MatrixOperand) -> MResult<MatrixOperand> {
        match operand {
            MatrixOperand::Scalar(expression) => self.scalar_child(expression).map(MatrixOperand::Scalar),
            MatrixOperand::Matrix(expression) => self.matrix_child(expression).map(MatrixOperand::Matrix),
        }
    }

    fn slice_spec(&mut self, spec: &SliceSpec) -> MResult<SliceSpec> {
        if !self.options.deep { return Ok(spec.clone()); }
        // Stop is derived from the affine range, not an independently
        // re-normalized raw Python bound. Preserve the selected progression.
        SliceSpec::parts(self.scalar(&spec.start)?, self.scalar(&spec.step)?, self.scalar(&spec.length)?)
    }

    fn product(&mut self, coefficient: Expr, factors: &[MatrixExpr]) -> MResult<MatrixExpr> {
        // Canonical MatMul combines adjacent explicit matrices. Reserve that
        // dense arithmetic too; it is not a single cheap abstract monomial.
        let mut pending: Vec<_> = factors.iter().rev().collect();
        let mut previous: Option<(usize, usize)> = None;
        while let Some(factor) = pending.pop() {
            match factor.kind() {
                MatrixKind::Mul(_, nested) => pending.extend(nested.iter().rev()),
                MatrixKind::Identity => {},
                MatrixKind::Explicit(matrix) => {
                    if let Some((rows, inner)) = previous {
                        self.budget.claim_product(rows, matrix.cols()).map_err(failure)?;
                        let entries = rows.checked_mul(matrix.cols()).ok_or_else(|| MatrixError::Value("Matrix expansion allocation limit exceeded".into()))?;
                        self.budget.claim_terms(entries).map_err(failure)?;
                        self.budget.claim_product(entries, inner).map_err(failure)?;
                        previous = Some((rows, matrix.cols()));
                    } else { previous = Some((matrix.rows(), matrix.cols())); }
                }
                _ => previous = None,
            }
        }
        let mut result = MatrixExpr::mul(coefficient, factors, true)?;
        // Multiplying distributed terms can introduce new scalar products,
        // e.g. (x+y)*A times (x-y)*B. Expand their combined coefficient too.
        if self.options.deep {
            if let MatrixKind::Mul(coefficient, factors) = result.kind() {
                let expanded = self.scalar(coefficient)?;
                if !expanded.same_presentation(coefficient) {
                    result = MatrixExpr::mul(expanded, factors, true)?;
                }
            }
        }
        Ok(result)
    }

    fn distributed_product(&mut self, coefficient: Expr, factors: Vec<MatrixExpr>) -> MResult<MatrixExpr> {
        if !self.options.mul || !factors.iter().any(|factor| matches!(factor.kind(), MatrixKind::Add(_))) {
            return self.product(coefficient, &factors);
        }
        // Check every intermediate Cartesian product before reserving space.
        let mut count = 1usize;
        for (position, factor) in factors.iter().enumerate() {
            let choices = match factor.kind() { MatrixKind::Add(terms) => terms.len(), _ => 1 };
            self.budget.claim_product(count, choices).map_err(failure)?;
            count = count.checked_mul(choices).ok_or_else(|| MatrixError::Value("Matrix expansion term limit exceeded".into()))?;
            self.budget.claim_terms(count).map_err(failure)?;
            if choices > 1 { self.budget.claim_product(count, position + 1).map_err(failure)?; }
        }
        self.budget.claim_product(count, factors.len()).map_err(failure)?;
        let mut terms = vec![Vec::new()];
        for factor in factors {
            if let MatrixKind::Add(choices) = factor.kind() {
                let mut next = Vec::with_capacity(terms.len() * choices.len());
                for prefix in terms {
                    for choice in choices {
                        let mut chain = prefix.clone();
                        chain.push(choice.clone());
                        next.push(chain);
                    }
                }
                terms = next;
            } else {
                // Appending a single choice does not copy an ever longer
                // prefix for every factor in a long ordered product.
                for prefix in &mut terms { prefix.push(factor.clone()); }
            }
        }
        let products = terms.iter().map(|factors| self.product(coefficient.clone(), factors)).collect::<MResult<Vec<_>>>()?;
        MatrixExpr::add(&products, true)
    }

    fn matrix(&mut self, expression: &MatrixExpr) -> MResult<MatrixExpr> {
        if let Some((_, value)) = self.matrices.get(&expression.dag_id()) { return Ok(value.clone()); }
        if let Some(value) = self.budget.matrix_result(expression, self.options) { return Ok(value); }
        self.budget.claim_terms(1).map_err(failure)?;
        self.budget.claim_product(1, 1).map_err(failure)?;
        let result = match expression.kind() {
            MatrixKind::Symbol(name) => {
                let rows = self.scalar_child(expression.rows())?;
                let cols = self.scalar_child(expression.cols())?;
                if rows.same_presentation(expression.rows()) && cols.same_presentation(expression.cols()) { expression.clone() }
                else { MatrixExpr::symbol(name, rows, cols)? }
            }
            MatrixKind::Identity | MatrixKind::Zero | MatrixKind::One => {
                let rows = self.scalar_child(expression.rows())?;
                let cols = self.scalar_child(expression.cols())?;
                if rows.same_presentation(expression.rows()) && cols.same_presentation(expression.cols()) { expression.clone() }
                else { MatrixExpr::new(expression.kind().clone(), rows, cols) }
            }
            MatrixKind::Explicit(matrix) => {
                // A concrete matrix expands its entries, even with deep=False,
                // matching MatrixBase.expand. A nested matrix is visited only
                // when its abstract parent's deep flag permits it.
                self.budget.claim_terms(matrix.data().len()).map_err(failure)?;
                let entries = matrix.data().iter().map(|entry| self.scalar(entry)).collect::<MResult<Vec<_>>>()?;
                if entries.iter().zip(matrix.data()).all(|(a,b)| a.same_presentation(b)) { expression.clone() }
                else { MatrixExpr::explicit(Matrix::new(matrix.rows(), matrix.cols(), entries)?) }
            }
            MatrixKind::Add(terms) => {
                self.budget.claim_terms(terms.len()).map_err(failure)?;
                let terms = terms.iter().map(|term| self.matrix_child(term)).collect::<MResult<Vec<_>>>()?;
                MatrixExpr::add(&terms, true)?
            }
            MatrixKind::Mul(coefficient, factors) => {
                let coefficient = self.scalar_child(coefficient)?;
                let factors = factors.iter().map(|factor| self.matrix_child(factor)).collect::<MResult<Vec<_>>>()?;
                self.distributed_product(coefficient, factors)?
            }
            MatrixKind::Pow(base, exponent) => {
                let mapped = self.matrix_child(base)?;
                let exponent_mapped = self.scalar_child(exponent)?;
                if mapped.dag_id() == base.dag_id() && exponent_mapped.same_presentation(exponent) { expression.clone() }
                else { MatrixExpr::new(MatrixKind::Pow(mapped.clone(), exponent_mapped), mapped.rows().clone(), mapped.cols().clone()) }
            }
            MatrixKind::Transpose(base) | MatrixKind::Adjoint(base) | MatrixKind::Inverse(base) | MatrixKind::Conjugate(base) => {
                let mapped = self.matrix_child(base)?;
                if mapped.dag_id() == base.dag_id() { expression.clone() }
                else {
                    let (kind, rows, cols) = match expression.kind() {
                        MatrixKind::Transpose(_) => (MatrixKind::Transpose(mapped.clone()), mapped.cols().clone(), mapped.rows().clone()),
                        MatrixKind::Adjoint(_) => (MatrixKind::Adjoint(mapped.clone()), mapped.cols().clone(), mapped.rows().clone()),
                        MatrixKind::Inverse(_) => (MatrixKind::Inverse(mapped.clone()), mapped.rows().clone(), mapped.cols().clone()),
                        _ => (MatrixKind::Conjugate(mapped.clone()), mapped.rows().clone(), mapped.cols().clone()),
                    };
                    MatrixExpr::new(kind, rows, cols)
                }
            }
            MatrixKind::HadamardProduct(factors) | MatrixKind::Kronecker(factors) => {
                let mapped = factors.iter().map(|factor| self.matrix_child(factor)).collect::<MResult<Vec<_>>>()?;
                if MatrixExpr::same_children(&mapped, factors) { expression.clone() }
                else if matches!(expression.kind(), MatrixKind::HadamardProduct(_)) { MatrixExpr::hadamard_product(&mapped, false)? }
                else { MatrixExpr::kronecker(&mapped, false)? }
            }
            MatrixKind::HadamardPower(base, exponent) => {
                let b = self.operand(base)?;
                let e = self.operand(exponent)?;
                if b.same_presentation(base) && e.same_presentation(exponent) { expression.clone() }
                else { MatrixExpr::hadamard_power(b, e, false)? }
            }
            MatrixKind::BlockDiag(blocks) => {
                let mapped = blocks.iter().map(|block| self.matrix_child(block)).collect::<MResult<Vec<_>>>()?;
                if MatrixExpr::same_children(&mapped, blocks) { expression.clone() }
                else { MatrixExpr::block_diag(&mapped)? }
            }
            MatrixKind::Block(grid) => {
                let mapped = grid.data().iter().map(|block| self.matrix_child(block)).collect::<MResult<Vec<_>>>()?;
                if MatrixExpr::same_children(&mapped, grid.data()) { expression.clone() }
                else { MatrixExpr::block(BlockGrid::new(grid.shape().0, grid.shape().1, mapped)?)? }
            }
            MatrixKind::Diag(base) | MatrixKind::DiagonalOf(base) | MatrixKind::Diagonal(base) => {
                let mapped = self.matrix_child(base)?;
                if mapped.dag_id() == base.dag_id() { expression.clone() }
                else { match expression.kind() { MatrixKind::Diag(_) => MatrixExpr::diag(&mapped)?, MatrixKind::DiagonalOf(_) => MatrixExpr::diagonal_of(&mapped)?, _ => MatrixExpr::diagonal(&mapped)? } }
            }
            MatrixKind::Slice(base, rows, cols) => {
                let mapped = self.matrix_child(base)?;
                let r = self.slice_spec(rows)?;
                let c = self.slice_spec(cols)?;
                if mapped.dag_id() == base.dag_id() && r.same_presentation(rows) && c.same_presentation(cols) { expression.clone() }
                else { mapped.slice_normalized(r, c)? }
            }
        };
        self.budget.remember_matrix(expression, self.options, &result);
        self.matrices.insert(expression.dag_id(), (expression.clone(), result.clone()));
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn symbol(name: &str) -> MatrixExpr { MatrixExpr::symbol(name, Expr::int(2), Expr::int(2)).unwrap() }
    fn expand(e: &MatrixExpr, options: ExpandOptions, limit: usize) -> MResult<MatrixExpr> {
        e.expand_with_options(&options, &mut ExpandBudget::new(limit))
    }
    #[test]
    fn ordered_distribution_and_scalar_coefficients_share_options() {
        let (a,b,c) = (symbol("A"), symbol("B"), symbol("C"));
        let coefficient = crate::parser::parse_str("(x+y)**2").unwrap();
        let sum = MatrixExpr::add(&[a.clone(),b.clone()],true).unwrap();
        let input = MatrixExpr::mul(coefficient.clone(), &[sum.clone(), c.clone()], true).unwrap();
        let expanded_coefficient = crate::parser::parse_str("x**2+2*x*y+y**2").unwrap();
        let expected = MatrixExpr::add(&[
            MatrixExpr::mul(expanded_coefficient.clone(), &[a,c.clone()], true).unwrap(),
            MatrixExpr::mul(expanded_coefficient.clone(), &[b,c], true).unwrap(),
        ], true).unwrap();
        assert_eq!(input.expand().unwrap(),expected);
        let no_mul = expand(&input,ExpandOptions { mul:false, ..Default::default() },100).unwrap();
        assert!(matches!(no_mul.kind(),MatrixKind::Mul(c,_) if c==&expanded_coefficient));
        let shallow = expand(&input,ExpandOptions { deep:false, ..Default::default() },100).unwrap();
        assert!(matches!(shallow.kind(),MatrixKind::Add(terms) if terms.iter().all(|t|matches!(t.kind(),MatrixKind::Mul(c,_) if c==&coefficient))));
    }
    #[test]
    fn held_wrappers_and_matrix_powers_remain_held() {
        let (a,b,c)=(symbol("A"),symbol("B"),symbol("C"));
        let parent=MatrixExpr::mul(Expr::one(),&[MatrixExpr::add(&[a,b],true).unwrap(),c],true).unwrap();
        for held in [parent.transpose(false).unwrap(),parent.adjoint(false).unwrap(),parent.inverse(false).unwrap(),parent.pow(Expr::int(2),false).unwrap()] {
            let expanded=held.expand().unwrap();
            assert_eq!(expanded.kind_name(),held.kind_name());
            assert!(matches!(expanded.children()[0].kind(),MatrixKind::Add(_)));
            let shallow=expand(&held,ExpandOptions { deep:false,..Default::default() },100).unwrap();
            assert_eq!(shallow,held);
        }
        let held=MatrixExpr::hadamard_product(&[parent.clone(),parent.clone()],false).unwrap();
        assert!(matches!(held.expand().unwrap().kind(),MatrixKind::HadamardProduct(f) if f.len()==2));
    }
    #[test]
    fn checked_cartesian_products_reject_limits_without_expanding_powers() {
        let (a,b)=(symbol("A"),symbol("B"));
        let factors=(1..13).map(|i|MatrixExpr::add(&[a.clone(),MatrixExpr::mul(Expr::int(i),&[b.clone()],true).unwrap()],true).unwrap()).collect::<Vec<_>>();
        let expression=MatrixExpr::mul(Expr::one(),&factors,false).unwrap();
        assert!(expand(&expression,ExpandOptions::default(),100).is_err());
        assert!(expand(&expression,ExpandOptions { mul:false,..Default::default() },100).is_ok());
        let power=MatrixExpr::add(&[a,b],true).unwrap().pow(Expr::int(1_000_000_000),true).unwrap();
        assert_eq!(expand(&power,ExpandOptions::default(),100).unwrap(),power);
    }
}
