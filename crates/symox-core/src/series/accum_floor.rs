//! The error of `Expr.series` for `floor`/`ceiling` of an argument that
//! diverges at the expansion point (`floor(x).series(x, oo)`).
//!
//! The reference maps the point to `0+`, where `floor._eval_nseries` and
//! `ceiling._eval_nseries` of an infinite `arg0` return the expansion of
//! the argument plus `AccumBounds(-1, 0)` (`AccumBounds(0, 1)`). `series`
//! then combines that bound with terms in `x`, and `AccumBounds` raises
//! `ValueError("Only real AccumulationBounds are supported")` unless the
//! variable is known to be real of one sign (a positive symbol gives an
//! `AccumBounds` value, not reproduced). It is raised whenever such a node
//! is reached from the root through sums and products
//! (`floor(1/x) + sin(x)`, `floor(x)*sin(1/x)` at `oo`); an integer power
//! of it raises `NotImplementedError` in `Pow._eval_nseries`
//! (`1/floor(x)` at `oo`). Under any other function the reference result
//! differs (`exp(floor(x))` is returned unexpanded) and nothing is
//! decided here.

use super::engine::depends;
use crate::assumptions::{is_negative, is_positive};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::limits::{limit, Dir};
use crate::ops::{add, neg, pow};
use crate::subs::subs1;

/// The exception the reference raises.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AccumError {
    /// `ValueError("Only real AccumulationBounds are supported")`
    Value,
    /// `NotImplementedError` of a power of the bounded node
    NotImplemented,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Walk {
    /// no step node of a diverging argument
    Free,
    Raise(AccumError),
    /// a node the rule does not decide
    Unknown,
}

fn diverging_step(e: &Expr, x: &Expr) -> Option<bool> {
    let Kind::Func(Func::Floor | Func::Ceiling, a) = e.kind() else { return None };
    if a.len() != 1 || !depends(&a[0], x) {
        return None;
    }
    let l = limit(&a[0], x, &Expr::zero(), Dir::Plus).ok();
    Some(matches!(l, Some(v) if v.is_constant(Constant::Infinity) || v.is_constant(Constant::NegInfinity)))
}

fn has_step(e: &Expr, x: &Expr) -> bool {
    depends(e, x) && (matches!(e.kind(), Kind::Func(Func::Floor | Func::Ceiling, _)) || e.args().iter().any(|a| has_step(a, x)))
}

fn walk(e: &Expr, x: &Expr) -> Walk {
    if !has_step(e, x) {
        return Walk::Free;
    }
    if let Some(div) = diverging_step(e, x) {
        return if div { Walk::Raise(AccumError::Value) } else { Walk::Unknown };
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) => {
            let mut out = Walk::Free;
            for a in e.args() {
                match walk(&a, x) {
                    Walk::Free => {}
                    Walk::Unknown => return Walk::Unknown,
                    Walk::Raise(k) => {
                        // the first factor in the reference order is not
                        // tracked: mixed kinds are not decided
                        if matches!(out, Walk::Raise(j) if j != k) {
                            return Walk::Unknown;
                        }
                        out = Walk::Raise(k);
                    }
                }
            }
            out
        }
        Kind::Pow(b, k) => {
            if depends(k, x) || k.as_i64().is_none() {
                return Walk::Unknown;
            }
            match walk(b, x) {
                Walk::Raise(_) => Walk::Raise(AccumError::NotImplemented),
                w => w,
            }
        }
        _ => Walk::Unknown,
    }
}

/// The reference exception of `e.series(x, x0, dir=dir)` when it comes from
/// a bounded `floor`/`ceiling` node; `None` when not decided.
pub fn accum_series_error(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> Option<AccumError> {
    if !x.is_symbol() || is_positive(x) == Some(true) || is_negative(x) == Some(true) {
        return None;
    }
    let u = Expr::dummy("x", crate::symbol::Assumptions::NONE);
    let v = if x0.is_constant(Constant::Infinity) {
        pow(&u, &Expr::minus_one())
    } else if x0.is_constant(Constant::NegInfinity) {
        neg(&pow(&u, &Expr::minus_one()))
    } else if depends(x0, x) || x0.is_constant(Constant::ComplexInfinity) {
        return None;
    } else if minus {
        add(&[x0.clone(), neg(&u)])
    } else {
        add(&[x0.clone(), u.clone()])
    };
    let e0 = subs1(e, x, &v);
    match walk(&e0, &u) {
        Walk::Raise(k) => Some(k),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::functions::func;

    #[test]
    fn floor_at_infinity_is_a_value_error() {
        let x = Expr::symbol("x");
        let e = func(Func::Floor, vec![x.clone()]);
        assert_eq!(accum_series_error(&e, &x, &Expr::oo(), false), Some(AccumError::Value));
        let inv = pow(&e, &Expr::minus_one());
        assert_eq!(accum_series_error(&inv, &x, &Expr::oo(), false), Some(AccumError::NotImplemented));
        let ex = func(Func::Exp, vec![e]);
        assert_eq!(accum_series_error(&ex, &x, &Expr::oo(), false), None);
    }
}
