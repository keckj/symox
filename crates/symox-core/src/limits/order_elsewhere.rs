//! An order term `O(g, (x, p))` in a limit at another point `x0 != p`.
//!
//! The reference `Limit.doit` takes the leading term of the expression at
//! `x0`, and `Order._eval_as_leading_term` of `O(g)` at its own finite
//! point `p` where `g -> 0` is `0` there: the order term drops out
//! (`limit(1 + O(x**2), x, 1)` is 1, also at `oo`, at `-1` and at a
//! symbolic point). An order term at an infinite point, or one whose `g`
//! does not vanish at `p`, keeps the limit unevaluated (the reference
//! returns it unevaluated or an `O(1)` quirk).

use super::{limit, Dir};
use crate::expr::{Constant, Expr};

/// `true` when `O(g, (x, p))` drops out of a limit at `x0`.
pub(super) fn drops_out(g: &Expr, vars: &[(Expr, Expr)], x: &Expr, x0: &Expr) -> bool {
    let [(v, p)] = vars else { return false };
    if v != x || p == x0 {
        return false;
    }
    let infinite = |e: &Expr| {
        e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity) || e.is_constant(Constant::ComplexInfinity)
    };
    if infinite(p) || crate::traverse::has(p, x) {
        return false;
    }
    matches!(limit(g, x, p, Dir::Plus), Ok(l) if l.is_zero())
}
