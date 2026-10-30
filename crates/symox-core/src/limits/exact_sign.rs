//! The sign of a constant that numerical evaluation cannot certify because
//! it is an exact zero hidden behind logarithms of perfect powers:
//! `log(4)/log(2) - 2`, the exponent difference of `4**x/2**(2*x)` in the
//! Gruntz leading term. The reference `sign(e, x)` applies `logcombine`
//! and its assumptions; here the logarithms of rational numbers are
//! expanded (`log(4) = 2*log(2)`, an identity of the real logarithm of a
//! positive rational), the result is put over a common denominator and its
//! sign is decided again. Nothing is decided without proof: the zero must
//! come out of the canonical arithmetic and a nonzero value is decided by
//! the certified numerical sign.

use crate::expr::{Expr, Func, Kind};
use crate::series::coeff::C;
use crate::series::engine::sign_of;

fn has_number_log(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Log, a) if a.len() == 1 && a[0].is_number() => true,
        _ => !e.is_atom() && e.args().iter().any(has_number_log),
    }
}

/// Whether `e` is an exact zero once the logarithms of rationals are
/// expanded; checked before the numerical sign, which cannot certify a zero
/// and would only give up after raising its precision.
pub fn vanishes_after_log_expansion(e: &Expr) -> bool {
    if !has_number_log(e) {
        return false;
    }
    let ex = crate::simp::expand_log(e, false);
    ex != *e && (ex.is_zero() || crate::simplify::together(&ex).is_zero())
}

/// `Some(sign)` once the logarithms of rationals of `e` (free of the limit
/// variable) are expanded, `None` when it stays undecided.
pub fn sign_after_log_expansion(e: &Expr) -> Option<i32> {
    if !has_number_log(e) {
        return None;
    }
    let ex = crate::simp::expand_log(e, false);
    if ex == *e {
        return None;
    }
    for cand in [ex.clone(), crate::simplify::together(&ex)] {
        if cand.is_zero() {
            return Some(0);
        }
        if let Some(s) = sign_of(&C::from_expr(cand.clone())) {
            return Some(s);
        }
        if crate::assumptions::is_positive(&cand) == Some(true) {
            return Some(1);
        }
        if crate::assumptions::is_negative(&cand) == Some(true) {
            return Some(-1);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::{div, sub};

    #[test]
    fn perfect_power_logs_cancel() {
        let l = |n: i64| crate::functions::log(&Expr::int(n));
        let e = sub(&div(&l(4), &l(2)), &Expr::int(2));
        assert_eq!(sign_after_log_expansion(&e), Some(0));
        assert!(vanishes_after_log_expansion(&e));
        let e = sub(&div(&l(8), &l(2)), &Expr::int(2));
        assert_eq!(sign_after_log_expansion(&e), Some(1));
    }
}
