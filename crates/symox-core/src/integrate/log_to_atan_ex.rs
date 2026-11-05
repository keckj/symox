//! SymPy's `log_to_atan(f, g)` over `EX`.
//!
//! The reference converts `f + I*g` logarithms to arctangents by the
//! Euclidean recursion
//!
//! ```text
//! if f.degree() < g.degree(): f, g = -g, f
//! p, q = f.div(g)
//! if q.is_zero: return 2*atan(p.as_expr())
//! s, t, h = g.gcdex(-f)
//! u = (f*s + g*t).quo(h)
//! return 2*atan(u.as_expr()) + log_to_atan(s, t)
//! ```
//!
//! on `EX` polynomials, whose arithmetic is [`ExpressionPoly`]'s.
//!
//! When `g` is a nonzero constant the reference remainder is zero and each
//! coefficient of `p` is the `EX` field quotient of a coefficient of `f` by
//! the constant, `Expression(a/c).simplify()` (`cancel` then a shallow
//! `expand`). That quotient is computed coefficient by coefficient
//! ([`constant_divisor`]): the native `EX` normal form does not always
//! reduce `(a/c)*c - a` to zero on nested radicals, so the remainder pass of
//! the general division is skipped there. Fed the reference operands, this
//! reproduces the reference result.
//!
//! Exception (kept from the previous Python driver, now computed here): the
//! native `EX` normal form does not reach SymPy's `cancel` form of nested
//! Cardano radicals, and on such inputs the field quotient can be much
//! longer than SymPy's. The constant-divisor step therefore keeps the
//! shorter (printed length, ties to the plain quotient) of the field quotient
//! and the plain quotient `a/c` of each coefficient; both are exact.

use crate::expr::Expr;
use crate::ops;
use crate::poly::convert::{PResult, PolyError};
use crate::poly::expression_domain::{ExpressionPoly, EX};

/// `2*atan(p)` for `f` given by its coefficients (highest degree first) in
/// `x`, divided by the constant `c`.
pub fn constant_divisor(coefficients: &[Expr], c: &Expr, x: &Expr) -> PResult<Expr> {
    let field = ops::mul(&[Expr::int(2), crate::functions::atan(&constant_quotient(coefficients, c, x)?)]);
    let plain = ops::mul(&[Expr::int(2), crate::functions::atan(&plain_quotient(coefficients, c, x))]);
    if crate::printer::to_str(&field).len() < crate::printer::to_str(&plain).len() {
        return Ok(field);
    }
    Ok(plain)
}

fn constant_quotient(coefficients: &[Expr], c: &Expr, x: &Expr) -> PResult<Expr> {
    let divisor = EX::new(c.clone());
    let degree = coefficients.len().saturating_sub(1);
    let mut terms = Vec::with_capacity(coefficients.len());
    for (index, a) in coefficients.iter().enumerate() {
        let q = EX::new(a.clone()).div(&divisor)?;
        let power = (degree - index) as i64;
        terms.push(ops::mul(&[q.as_expr().clone(), ops::pow(x, &Expr::int(power))]));
    }
    Ok(ops::add(&terms))
}

/// `sum(a/c * x**k)` with the automatic evaluation of `a/c` only.
fn plain_quotient(coefficients: &[Expr], c: &Expr, x: &Expr) -> Expr {
    let inverse = ops::pow(c, &Expr::int(-1));
    let degree = coefficients.len().saturating_sub(1);
    let terms: Vec<Expr> = coefficients
        .iter()
        .enumerate()
        .map(|(index, a)| {
            let power = (degree - index) as i64;
            ops::mul(&[ops::mul(&[a.clone(), inverse.clone()]), ops::pow(x, &Expr::int(power))])
        })
        .collect();
    ops::add(&terms)
}

/// `log_to_atan(f, g)` for univariate `EX` polynomials of the same
/// generator, `g` nonzero. `check` is called once per Euclidean step.
pub fn log_to_atan<E: From<PolyError>>(
    f: &ExpressionPoly,
    g: &ExpressionPoly,
    check: &mut impl FnMut() -> Result<(), E>,
) -> Result<Expr, E> {
    if f.gens.len() != 1 || f.gens != g.gens {
        return Err(PolyError::Unsupported("univariate EX polynomials of one generator expected".into()).into());
    }
    let x = f.gens[0].clone();
    let (mut f, mut g) = (f.clone(), g.clone());
    let mut terms = Vec::new();
    loop {
        check()?;
        if f.degree(0) < g.degree(0) {
            (f, g) = (g.neg(), f);
        }
        if g.is_zero() {
            // `f.div(0)` raises in the reference.
            return Err(PolyError::Unsupported("log_to_atan with a zero divisor".into()).into());
        }
        if g.degree(0) == 0 {
            terms.push(constant_divisor(&f.all_coeffs(), &g.lc(), &x)?);
            break;
        }
        let (p, q) = f.div(&g)?;
        if q.is_zero() {
            terms.push(ops::mul(&[Expr::int(2), crate::functions::atan(&p.as_expr())]));
            break;
        }
        let (s, t, h) = g.gcdex(&f.neg())?;
        let (u, _) = f.mul(&s)?.add(&g.mul(&t)?)?.div(&h)?;
        terms.push(ops::mul(&[Expr::int(2), crate::functions::atan(&u.as_expr())]));
        (f, g) = (s, t);
    }
    Ok(ops::add(&terms))
}
