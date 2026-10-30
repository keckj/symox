//! `Order` of a leading term whose coefficient is a rational function of
//! `log(x)` (`O(log(x) + 1/log(x))` is `O((log(x)**2 + 1)/log(x))`).
//!
//! `Order.__new__` factors a sum first and then takes its leading term,
//! which keeps the `log(x)` factors (they are constant for
//! `as_leading_term`) and drops only the factors free of `x`. The leading
//! term `x**a*c(log(x))` is therefore `x**a` times the factored `c`, with
//! the constant factors removed, when `c` is not a single power of
//! `log(x)`; for `c = k*log(x)**b` it is `x**a*log(x)**b`, also for a
//! negative `b` (`O(x + 1/log(x))` is `O(1/log(x))`). A sum that stays a
//! sum once factored is left to the dominant-term path.

use super::engine::{degree_in, depends};
use crate::expr::{Expr, Kind};
use crate::ops::{mul, mul2, mul_factor, pow};
use crate::subs::xreplace;

/// `Some(b)` when `c = k*l**b` with `k` free of `l`.
fn log_power(c: &Expr, l: &Expr) -> Option<i64> {
    if c == l {
        return Some(1);
    }
    match c.kind() {
        Kind::Pow(b, p) if b == l => p.as_i64(),
        Kind::Mul(_, facs) => {
            let mut found = None;
            for (b, p) in facs {
                let f = mul_factor(b, p);
                if !depends(&f, l) {
                    continue;
                }
                if found.is_some() {
                    return None;
                }
                found = Some(log_power(&f, l)?);
            }
            found
        }
        _ => None,
    }
}

/// The order term at `x = 0+` of `e0` when its leading coefficient depends
/// on `log(x)` as described in the module notes; `None` otherwise.
pub fn log_coefficient_order(e0: &Expr, x: &Expr) -> Option<Expr> {
    if e0.may_have_order() || !depends(e0, x) {
        return None;
    }
    let lt = super::leading_term(e0, x, &Expr::zero(), false).ok()?;
    let (c, l) = (&lt.coeff_logt, &lt.logt);
    if !depends(c, l) || depends(c, x) {
        return None;
    }
    // a sum that stays a sum once factored keeps its dominant term
    // (`Add.extract_leading_order`: `O(log(x) + 1)` is `O(log(x))`)
    let whole = crate::poly::api::factor(&crate::poly::api::together(e0)).ok()?;
    if matches!(whole.kind(), Kind::Add(..)) {
        return None;
    }
    let xa = pow(x, &Expr::num(lt.exp.clone()));
    let logx = crate::functions::log(x);
    let back = |e: &Expr| xreplace(e, &std::iter::once((l.clone(), logx.clone())).collect());
    if let Some(b) = log_power(c, l) {
        // positive powers are already built by the monomial path
        if b >= 0 {
            return None;
        }
        return Some(mul2(&xa, &pow(&logx, &Expr::int(b))));
    }
    let t = crate::poly::api::together(c);
    let (num, den) = crate::simplify::numer_denom(&t);
    degree_in(&num, l)?;
    degree_in(&den, l)?;
    let f = crate::poly::api::factor(&t).ok()?;
    let keep: Vec<Expr> = match f.kind() {
        Kind::Mul(_, facs) => facs.iter().map(|(b, p)| mul_factor(b, p)).filter(|g| depends(g, l)).collect(),
        _ => vec![f.clone()],
    };
    Some(mul2(&xa, &back(&mul(&keep))))
}
