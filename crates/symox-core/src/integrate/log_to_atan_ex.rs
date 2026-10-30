//! SymPy's `log_to_atan(f, g)` when `g` is a nonzero `EX` constant.
//!
//! The reference returns `2*atan(p.as_expr())` with `p, q = f.div(g)`: the
//! remainder is zero and each coefficient of `p` is the `EX` field quotient
//! of a coefficient of `f` by the constant, `Expression(a/c).simplify()`
//! (`cancel` then a shallow `expand`). The division is done coefficient by
//! coefficient, which is the reference quotient without the remainder pass.

use crate::expr::Expr;
use crate::ops;
use crate::poly::expression_domain::EX;
use crate::poly::convert::PResult;

/// `2*atan(p)` for `f` given by its coefficients (highest degree first) in
/// `x`, divided by the constant `c`.
pub fn constant_divisor(coefficients: &[Expr], c: &Expr, x: &Expr) -> PResult<Expr> {
    let divisor = EX::new(c.clone());
    let degree = coefficients.len().saturating_sub(1);
    let mut terms = Vec::with_capacity(coefficients.len());
    for (index, a) in coefficients.iter().enumerate() {
        let q = EX::new(a.clone()).div(&divisor)?;
        let power = (degree - index) as i64;
        terms.push(ops::mul(&[q.as_expr().clone(), ops::pow(x, &Expr::int(power))]));
    }
    Ok(ops::mul(&[Expr::int(2), crate::functions::atan(&ops::add(&terms))]))
}
