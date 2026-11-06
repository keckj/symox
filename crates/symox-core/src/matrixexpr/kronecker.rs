//! `combine_kronecker` (sympy/matrices/expressions/kronecker.py).
//!
//! SymPy applies `exhaust(bottom_up(exhaust(condition(haskron, typed({MatAdd:
//! kronecker_mat_add, MatMul: kronecker_mat_mul, MatPow:
//! kronecker_mat_pow})))))` and ends with `doit()`. Two of its rules are
//! not identities and are deliberately narrowed here:
//!
//! * `kronecker_mat_add` rewrites `KroneckerProduct(A, B) +
//!   KroneckerProduct(C, D)` as `KroneckerProduct(A + C, B + D)`, which is
//!   false in general (the cross terms `A⊗D + C⊗B` appear); sums are left
//!   unchanged.
//! * `kronecker_mat_pow` distributes any exponent over the factors; the
//!   mixed-product property only proves it for integer exponents (a principal
//!   power of a Kronecker product is not the product of principal powers, e.g.
//!   square roots of `-1` factors), so only integer exponents are distributed.
use super::{MatrixExpr, MatrixKind};
use crate::expr::Expr;
use crate::matrix::MResult;

/// SymPy's `exhaust` has no bound; a non-terminating rewrite would be a bug,
/// reported as an error instead of looping.
const MAX_PASSES: usize = 256;

impl MatrixExpr {
    fn has_kronecker(&self) -> bool {
        let mut found = false;
        self.visit_nodes(&mut |node| if matches!(node.kind(), MatrixKind::Kronecker(_)) { found = true; });
        found
    }

    /// `KroneckerProduct.has_matching_shape`.
    fn kronecker_matching_shape(a: &[Self], b: &[Self], outer_a: &Self, outer_b: &Self) -> bool {
        outer_a.cols() == outer_b.rows() && a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.cols() == y.rows())
    }

    /// One application of the typed rule, `None` when it does not apply.
    fn kronecker_rule(&self) -> MResult<Option<Self>> {
        match self.kind() {
            // kronecker_mat_mul: adjacent products are merged with
            // `A._kronecker_mul(B)`, the factor-wise products when the shapes
            // match and `A*B` otherwise; the result is `factor*MatMul(...)`.
            MatrixKind::Mul(coefficient, factors) => {
                let mut matrices = factors.clone();
                let mut changed = false;
                let mut i = 0;
                while i + 1 < matrices.len() {
                    let (a, b) = (&matrices[i], &matrices[i + 1]);
                    if let (MatrixKind::Kronecker(x), MatrixKind::Kronecker(y)) = (a.kind(), b.kind()) {
                        let merged = if Self::kronecker_matching_shape(x, y, a, b) {
                            let products = x.iter().zip(y).map(|(p, q)| Self::mul(Expr::one(), &[p.clone(), q.clone()], true)).collect::<MResult<Vec<_>>>()?;
                            Self::kronecker(&products, false)?
                        } else {
                            Self::mul(Expr::one(), &[a.clone(), b.clone()], true)?
                        };
                        matrices[i] = merged;
                        matrices.remove(i + 1);
                        changed = true;
                    } else {
                        i += 1;
                    }
                }
                if !changed { return Ok(None); }
                let product = if matrices.len() == 1 && coefficient.is_one() { matrices.remove(0) } else { Self::mul(coefficient.clone(), &matrices, false)? };
                Ok(Some(product))
            }
            // kronecker_mat_pow for square factors, integer exponents only.
            MatrixKind::Pow(base, exponent) => {
                let MatrixKind::Kronecker(factors) = base.kind() else { return Ok(None) };
                if crate::assumptions::is_integer(exponent) != Some(true) { return Ok(None); }
                if factors.iter().any(|f| f.rows() != f.cols()) { return Ok(None); }
                let powers = factors.iter().map(|f| f.pow(exponent.clone(), false)).collect::<MResult<Vec<_>>>()?;
                Ok(Some(Self::kronecker(&powers, false)?))
            }
            _ => Ok(None),
        }
    }

    /// `exhaust(condition(haskron, rule))` at one node.
    fn kronecker_exhaust(&self) -> MResult<Self> {
        let mut node = self.clone();
        for _ in 0..MAX_PASSES {
            if !node.has_kronecker() { return Ok(node); }
            match node.kronecker_rule()? {
                Some(next) if next != node => node = next,
                _ => return Ok(node),
            }
        }
        Err(crate::matrix::MatrixError::Value("combine_kronecker did not reach a fixed point".into()))
    }

    /// `bottom_up(...)`: the children first, then the node.
    fn kronecker_bottom_up(&self) -> MResult<Self> {
        let old = self.children();
        let children = old.iter().map(Self::kronecker_bottom_up).collect::<MResult<Vec<_>>>()?;
        let rebuilt = if Self::same_children(&children, &old) { self.clone() } else {
            let mut index = 0;
            self.map_children(&mut |_| { let child = children[index].clone(); index += 1; child }, &mut |e| e.clone())?
        };
        rebuilt.kronecker_exhaust()
    }

    /// `combine_kronecker(expr)`: the rewrite to a fixed point, then `doit()`.
    pub fn combine_kronecker(&self) -> MResult<Self> {
        let mut expr = self.clone();
        for _ in 0..MAX_PASSES {
            let next = expr.kronecker_bottom_up()?;
            if next == expr { return expr.doit(true); }
            expr = next;
        }
        Err(crate::matrix::MatrixError::Value("combine_kronecker did not reach a fixed point".into()))
    }
}
