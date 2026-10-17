//! Certified, condensed singular value decomposition over exact scalars.
use super::functions::budget::{raw_children, Budget};
use super::{MResult, Matrix, MatrixError};
use crate::assumptions::{exact_queries, is_negative, is_zero};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::func;
use crate::number::Number;
use crate::ops::neg;
use rustc_hash::{FxHashMap, FxHashSet};

fn unsupported(message: &str) -> MatrixError {
    MatrixError::NotImplemented(format!(
        "Cannot certify an exact singular value decomposition: {message}"
    ))
}

/// All zero/rank decisions are structural. In particular this implementation
/// does not use the generic QR/nullspace numeric pivot heuristics.
struct Context {
    budget: Budget,
    // Each root is proved positive by a nonzero Hermitian norm or a nonzero
    // eigenvalue of a Gram matrix. Conjugation can use that local proof even
    // when the global assumptions engine cannot establish the sign.
    positive: FxHashMap<Expr, Expr>,
}
impl Context {
    fn new() -> Self {
        Self {
            budget: Budget::new(),
            positive: FxHashMap::default(),
        }
    }
    fn clean(&mut self, value: &Matrix) -> MResult<Matrix> {
        self.budget.matrix(value)?;
        self.budget.allocate(value.data().len().saturating_mul(2))?;
        self.budget
            .work((value.data().len() as u64).saturating_mul(self.positive.len() as u64))?;
        let value = if self.positive.is_empty() {
            value.clone()
        } else {
            value.applyfunc(|entry| crate::subs::xreplace(entry, &self.positive))
        };
        self.budget.cancel(&value)
    }
    fn mul(&mut self, a: &Matrix, b: &Matrix) -> MResult<Matrix> {
        self.budget.product_preflight(a, b)?;
        let value = self.budget.mul(a, b)?;
        self.clean(&value)
    }
    fn scale(&mut self, a: &Matrix, scalar: &Expr) -> MResult<Matrix> {
        let value = self.budget.scale(a, scalar)?;
        self.clean(&value)
    }
    fn sub(&mut self, a: &Matrix, b: &Matrix) -> MResult<Matrix> {
        let minus = self.budget.scale(b, &Expr::minus_one())?;
        let difference = self.budget.add(a, &minus)?;
        self.clean(&difference)
    }
    fn adjoint(&mut self, value: &Matrix) -> MResult<Matrix> {
        self.budget.matrix(value)?;
        self.budget.allocate(value.data().len())?;
        let adjoint = value.adjoint();
        self.clean(&adjoint)
    }
    fn certify_equal(&mut self, a: &Matrix, b: &Matrix, message: &str) -> MResult<()> {
        if a.shape() != b.shape() || !self.sub(a, b)?.is_zero_matrix() {
            Err(unsupported(message))
        } else {
            Ok(())
        }
    }
    fn root_of_positive(&mut self, value: &Expr) -> MResult<Expr> {
        if is_zero(value) == Some(true) || is_negative(value) == Some(true) {
            return Err(unsupported(
                "a purported positive norm or eigenvalue vanished",
            ));
        }
        let root = self.budget.scalar_power(value, &Expr::rational(1, 2))?;
        let conjugate = func(Func::Conjugate, vec![root.clone()]);
        if conjugate != root {
            self.positive.insert(conjugate, root.clone());
        }
        Ok(root)
    }
    fn inverse(&mut self, value: &Expr) -> MResult<Expr> {
        let inverse = self.budget.scalar_power(value, &Expr::minus_one())?;
        let conjugate = func(Func::Conjugate, vec![value.clone()]);
        if conjugate == *value || self.positive.get(&conjugate) == Some(value) {
            // pow(sqrt(lambda), -1) may flatten to lambda**(-1/2).
            // Its conjugate then no longer contains conjugate(sqrt(lambda)).
            // Carry the same proved reality through the reciprocal explicitly.
            let conjugate_inverse = func(Func::Conjugate, vec![inverse.clone()]);
            if conjugate_inverse != inverse {
                self.positive.insert(conjugate_inverse, inverse.clone());
            }
        }
        Ok(inverse)
    }
    fn columns(&mut self, rows: usize, columns: &[Matrix]) -> MResult<Matrix> {
        self.budget.allocate(
            rows.checked_mul(columns.len())
                .ok_or_else(|| unsupported("matrix size overflow"))?,
        )?;
        let result = if columns.is_empty() {
            Matrix::zeros(rows, 0)
        } else {
            Matrix::hstack(&columns.iter().collect::<Vec<_>>())?
        };
        self.budget.matrix(&result)?;
        Ok(result)
    }

    fn zero(&mut self, value: &Expr) -> MResult<Option<bool>> {
        if let Some(zero) = is_zero(value) {
            return Ok(Some(zero));
        }
        // A nonzero exact real or imaginary component proves a complex
        // pivot nonzero without approximating its modulus.
        let real = func(Func::Re, vec![value.clone()]);
        let imaginary = func(Func::Im, vec![value.clone()]);
        self.budget.expression(&real)?;
        self.budget.expression(&imaginary)?;
        let r = is_zero(&real);
        let i = is_zero(&imaginary);
        if r == Some(false) || i == Some(false) {
            Ok(Some(false))
        } else if r == Some(true) && i == Some(true) {
            Ok(Some(true))
        } else {
            Ok(None)
        }
    }

    /// Partial column pivoting may defer an undecidable entry and use a later
    /// proven pivot. Rank is accepted only once every remaining row is zero.
    fn rref(&mut self, input: &Matrix) -> MResult<(Matrix, Vec<usize>)> {
        let mut matrix = self.clean(input)?;
        let (rows, cols) = matrix.shape();
        let mut pivots = Vec::new();
        for col in 0..cols {
            let row = pivots.len();
            if row == rows {
                break;
            }
            let mut pivot = None;
            for i in row..rows {
                if self.zero(matrix.get(i, col))? == Some(false) {
                    pivot = Some(i);
                    break;
                }
            }
            let Some(pivot) = pivot else {
                continue;
            };
            if pivot != row {
                for j in 0..cols {
                    let entry = matrix.get(row, j).clone();
                    matrix.set(row, j, matrix.get(pivot, j).clone());
                    matrix.set(pivot, j, entry);
                }
            }
            let inverse = self.inverse(matrix.get(row, col))?;
            let normalized = self.scale(&matrix.row(row), &inverse)?;
            for j in 0..cols {
                matrix.set(row, j, normalized.get(0, j).clone());
            }
            for i in 0..rows {
                if i == row || matrix.get(i, col).is_zero() {
                    continue;
                }
                let scaled = self.scale(&normalized, matrix.get(i, col))?;
                let difference = self.sub(&matrix.row(i), &scaled)?;
                for j in 0..cols {
                    matrix.set(i, j, difference.get(0, j).clone());
                }
            }
            pivots.push(col);
        }
        for i in pivots.len()..rows {
            for j in 0..cols {
                if self.zero(matrix.get(i, j))? != Some(true) {
                    return Err(unsupported(
                        "the rank depends on an undecidable zero condition",
                    ));
                }
            }
        }
        for (row, &pivot) in pivots.iter().enumerate() {
            for i in 0..rows {
                if matrix.get(i, pivot) != &if i == row { Expr::one() } else { Expr::zero() } {
                    return Err(unsupported("row reduction could not be certified"));
                }
            }
        }
        Ok((matrix, pivots))
    }
    fn kernel(&mut self, matrix: &Matrix) -> MResult<Vec<Matrix>> {
        let (reduced, pivots) = self.rref(matrix)?;
        let mut columns = Vec::new();
        for free in (0..matrix.cols()).filter(|col| !pivots.contains(col)) {
            self.budget.allocate(matrix.cols())?;
            let mut vector = vec![Expr::zero(); matrix.cols()];
            vector[free] = Expr::one();
            for (row, &pivot) in pivots.iter().enumerate() {
                vector[pivot] = neg(reduced.get(row, free));
            }
            columns.push(Matrix::column(vector));
        }
        Ok(columns)
    }
}

fn validate(input: &Matrix) -> MResult<()> {
    let mut pending: Vec<_> = input.data().to_vec();
    let mut seen = FxHashSet::default();
    let mut floating = false;
    while let Some(entry) = pending.pop() {
        if !seen.insert(entry.dag_id()) {
            continue;
        }
        match entry.kind() {
            Kind::Sym(symbol) if symbol.facts().get("finite") == Some(false) => {
                return Err(MatrixError::Value(
                    "Singular value decomposition requires finite entries".into(),
                ));
            }
            Kind::Func(Func::Undefined(head),_) if !head.internal&&head.facts().get("finite")==Some(false)=>{
                return Err(MatrixError::Value("Singular value decomposition requires finite entries".into()));
            }
            Kind::Const(
                Constant::Infinity
                | Constant::NegInfinity
                | Constant::ComplexInfinity
                | Constant::NaN,
            ) => {
                return Err(MatrixError::Value(
                    "Singular value decomposition requires finite entries".into(),
                ));
            }
            Kind::Num(Number::Float(value)) if !value.is_finite() => {
                return Err(MatrixError::Value(
                    "Singular value decomposition requires finite entries".into(),
                ));
            }
            Kind::Num(n) => floating |= n.is_float(),
            Kind::Add(c, terms) => {
                floating |= c.is_float() || terms.iter().any(|(_, c)| c.is_float());
                if matches!(c,Number::Float(v) if !v.is_finite())
                    || terms
                        .iter()
                        .any(|(_, c)| matches!(c,Number::Float(v) if !v.is_finite()))
                {
                    return Err(MatrixError::Value(
                        "Singular value decomposition requires finite entries".into(),
                    ));
                }
            }
            Kind::Mul(c, _) => {
                floating |= c.is_float();
                if matches!(c,Number::Float(v) if !v.is_finite()) {
                    return Err(MatrixError::Value(
                        "Singular value decomposition requires finite entries".into(),
                    ));
                }
            }
            _ => {}
        }
        pending.extend(raw_children(&entry));
    }
    if floating {
        return Err(unsupported(
            "floating-point input is not an exact scalar domain",
        ));
    }
    if input
        .data()
        .iter()
        .any(|entry| crate::decompose::commutativity(entry) != Some(true))
    {
        return Err(unsupported("entries must be commutative scalars"));
    }
    Ok(())
}

impl Matrix {
    /// Condensed exact SVD: `A = U*S*V.H`, where the columns of U and V
    /// are orthonormal and S contains only positive singular values.
    /// Unresolved spectra, rank conditions, and certificates are explicit errors.
    pub fn singular_value_decomposition(&self) -> MResult<(Matrix, Matrix, Matrix)> {
        exact_queries(|| self.svd_exact())
    }
    fn svd_exact(&self) -> MResult<(Matrix, Matrix, Matrix)> {
        let mut context = Context::new();
        context.budget.allocate(self.data().len())?;
        context.budget.matrix(self)?;
        validate(self)?;
        let (m, n) = self.shape();
        if self.is_empty() || self.is_zero_matrix() {
            return Ok((
                Matrix::zeros(m, 0),
                Matrix::zeros(0, 0),
                Matrix::zeros(n, 0),
            ));
        }
        let input = context.clean(self)?;
        let (_, pivots) = context.rref(&input)?;
        let rank = pivots.len();
        if rank == 0 {
            return Ok((
                Matrix::zeros(m, 0),
                Matrix::zeros(0, 0),
                Matrix::zeros(n, 0),
            ));
        }
        let adjoint = context.adjoint(&input)?;
        let gram = if m >= n {
            context.mul(&adjoint, &input)?
        } else {
            context.mul(&input, &adjoint)?
        };
        let dimension = gram.rows();
        let height = context.budget.matrix(&gram)?;
        context.budget.interpolation(dimension, height)?;
        context.budget.spectral_preflight(&gram)?;
        let mut eigenvalues = gram
            .eigenvals()
            .map_err(|error| unsupported(&format!("unresolved spectrum ({error})")))?;
        eigenvalues.sort_by(|a, b| crate::sortkey::cmp_default(&a.0, &b.0));
        if eigenvalues.iter().map(|(_, mult)| mult).sum::<usize>() != dimension {
            return Err(unsupported("the exact spectrum is incomplete"));
        }
        let zero_count: usize = eigenvalues
            .iter()
            .filter(|(value, _)| is_zero(value) == Some(true))
            .map(|(_, mult)| mult)
            .sum();
        if zero_count != dimension - rank {
            return Err(unsupported(
                "the spectrum cannot be split at zero with the certified rank",
            ));
        }
        let identity = context.budget.identity(dimension)?;
        let mut vectors = Vec::new();
        let mut singular = Vec::new();
        for (lambda, multiplicity) in eigenvalues {
            context.budget.expression(&lambda)?;
            if is_zero(&lambda) == Some(true) {
                continue;
            }
            // The Gram matrix is positive semidefinite. Its certified rank,
            // complete spectrum and exact zero multiplicity prove that each
            // remaining eigenvalue is strictly positive, without approximating it.
            let scaled_identity = context.scale(&identity, &lambda)?;
            let shifted = context.sub(&gram, &scaled_identity)?;
            let basis = context.kernel(&shifted)?;
            if basis.len() != multiplicity {
                return Err(unsupported("an eigenspace is incomplete"));
            }
            let sigma = context.root_of_positive(&lambda)?;
            let mut orthogonal: Vec<Matrix> = Vec::new();
            for vector in basis {
                let image = context.mul(&gram, &vector)?;
                let expected = context.scale(&vector, &lambda)?;
                context.certify_equal(
                    &image,
                    &expected,
                    "an eigenvector does not satisfy G*q=lambda*q",
                )?;
                let mut vector = vector;
                for previous in &orthogonal {
                    let h = context.adjoint(previous)?;
                    let projection = context.mul(&h, &vector)?;
                    let component = context.scale(previous, projection.get(0, 0))?;
                    vector = context.sub(&vector, &component)?;
                }
                let h = context.adjoint(&vector)?;
                let norm = context.mul(&h, &vector)?;
                // The kernel basis has independent unit free coordinates;
                // Gram-Schmidt preserves independence, hence this norm > 0.
                let root = context.root_of_positive(norm.get(0, 0))?;
                let inverse = context.inverse(&root)?;
                let normalized = context.scale(&vector, &inverse)?;
                singular.push(sigma.clone());
                vectors.push(normalized.clone());
                orthogonal.push(normalized);
            }
        }
        let known = context.columns(dimension, &vectors)?;
        let s = context.budget.diagonal(&singular)?;
        let inverse = singular
            .iter()
            .map(|s| context.inverse(s))
            .collect::<MResult<Vec<_>>>()?;
        let inverse = context.budget.diagonal(&inverse)?;
        let (u, v) = if m >= n {
            let left = context.mul(&input, &known)?;
            (context.mul(&left, &inverse)?, known)
        } else {
            let right = context.mul(&adjoint, &known)?;
            (known, context.mul(&right, &inverse)?)
        };
        let expected = context.budget.identity(rank)?;
        let uh = context.adjoint(&u)?;
        let vh = context.adjoint(&v)?;
        let uhu = context.mul(&uh, &u)?;
        let vhv = context.mul(&vh, &v)?;
        context.certify_equal(&uhu, &expected, "U.H*U is not the identity")?;
        context.certify_equal(&vhv, &expected, "V.H*V is not the identity")?;
        let us = context.mul(&u, &s)?;
        let reconstructed = context.mul(&us, &vh)?;
        context.certify_equal(&reconstructed, &input, "U*S*V.H does not reconstruct A")?;
        Ok((u, s, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn m(rows: &[&[i64]]) -> Matrix {
        Matrix::from_rows(
            rows.iter()
                .map(|row| row.iter().map(|n| Expr::int(*n)).collect())
                .collect(),
        )
        .unwrap()
    }
    fn certificate(a: &Matrix) {
        let (u, s, v) = a.singular_value_decomposition().unwrap();
        let mut c = Context::new();
        for i in 0..s.rows() {
            let root = s.get(i, i);
            let conj = func(Func::Conjugate, vec![root.clone()]);
            if conj != *root {
                c.positive.insert(conj, root.clone());
            }
        }
        let uh = c.adjoint(&u).unwrap();
        let vh = c.adjoint(&v).unwrap();
        let uhu = c.mul(&uh, &u).unwrap();
        let vhv = c.mul(&vh, &v).unwrap();
        c.certify_equal(&uhu, &Matrix::eye(s.rows()), "left")
            .unwrap();
        c.certify_equal(&vhv, &Matrix::eye(s.rows()), "right")
            .unwrap();
        let us = c.mul(&u, &s).unwrap();
        let output = c.mul(&us, &vh).unwrap();
        c.certify_equal(&output, a, "reconstruct").unwrap();
    }
    #[test]
    fn certified_real_complex_rectangular_and_repeated() {
        for a in [
            m(&[&[1, 2], &[2, 4]]),
            m(&[&[1, 2], &[3, 4]]),
            m(&[&[1, 2, 3], &[2, 4, 6]]),
            m(&[&[1, 0], &[0, 1], &[1, 1]]),
            m(&[&[0, 0, 0], &[0, 0, 0]]),
        ] {
            certificate(&a);
        }
        certificate(&Matrix::new(1, 1, vec![crate::ops::add2(&Expr::one(), &Expr::i())]).unwrap());
        let i = Expr::i();
        certificate(
            &Matrix::from_rows(vec![vec![Expr::one(), i.clone()], vec![i, Expr::one()]]).unwrap(),
        );
        let projector = Matrix::eye(3)
            .sub(&Matrix::ones(3, 3).scale(&Expr::rational(1, 3)))
            .unwrap();
        certificate(&projector);
    }
    #[test]
    fn proved_parameters_and_complex_rectangles() {
        let p = Expr::symbol_with(
            "p",
            crate::Assumptions::from_kwargs([("positive", true)]).unwrap(),
        );
        certificate(&m(&[&[1, 0, 2], &[0, 3, 0]]).scale(&p));
        let r = Expr::symbol_with(
            "r",
            crate::Assumptions::from_kwargs([("nonzero", true)]).unwrap(),
        );
        certificate(&Matrix::diag(&[r, Expr::zero()]));
        let z = Expr::symbol("z");
        let row = Matrix::from_rows(vec![vec![z.clone(), Expr::one()]]).unwrap();
        row.singular_value_decomposition().unwrap();
        row.transpose().singular_value_decomposition().unwrap();
        let complex = Matrix::from_rows(vec![
            vec![Expr::one(), Expr::i(), Expr::int(2)],
            vec![Expr::i(), Expr::one(), Expr::zero()],
        ])
        .unwrap();
        certificate(&complex);
        certificate(&complex.transpose());
        let v = Matrix::column(vec![Expr::one(), Expr::i(), Expr::one()]);
        let projector = Matrix::eye(3)
            .sub(&v.matmul(&v.adjoint()).unwrap().scale(&Expr::rational(1, 3)))
            .unwrap();
        certificate(&projector);
    }
    #[test]
    fn declared_nonfinite_symbols_are_rejected_before_spectrum() {
        for flags in [
            vec![("finite", false)],
            vec![("infinite", true)],
            vec![("extended_real", true), ("finite", false)],
        ] {
            let facts = crate::assumption_facts::FactSet::from_kwargs(flags).unwrap();
            let symbol = crate::Symbol::with_facts("z", facts).unwrap();
            let z = Expr::from_kind(Kind::Sym(symbol));
            let input = Matrix::new(1, 2, vec![z, Expr::one()]).unwrap();
            assert!(matches!(
                input.singular_value_decomposition(),
                Err(MatrixError::Value(_))
            ));
        }
        for flags in [vec![], vec![("extended_real", true)]] {
            let facts = crate::assumption_facts::FactSet::from_kwargs(flags).unwrap();
            let symbol = crate::Symbol::with_facts("z", facts).unwrap();
            let z = Expr::from_kind(Kind::Sym(symbol));
            assert!(Matrix::new(1, 2, vec![z, Expr::one()])
                .unwrap()
                .singular_value_decomposition()
                .is_ok());
        }
    }
    #[test]
    fn empties_rank_conditions_and_held_growth() {
        for (m, n) in [(0, 0), (0, 3), (4, 0)] {
            let (u, s, v) = Matrix::zeros(m, n).singular_value_decomposition().unwrap();
            assert_eq!(u.shape(), (m, 0));
            assert_eq!(s.shape(), (0, 0));
            assert_eq!(v.shape(), (n, 0));
        }
        assert!(Matrix::new(1, 1, vec![Expr::symbol("x")])
            .unwrap()
            .singular_value_decomposition()
            .is_err());
        let held = Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2), Expr::int(1_000_000_000)));
        assert!(Matrix::new(1, 1, vec![held])
            .unwrap()
            .singular_value_decomposition()
            .is_err());
        assert!(Matrix::new(1, 1, vec![Expr::num(Number::float(1.0))])
            .unwrap()
            .singular_value_decomposition()
            .is_err());
    }
}
