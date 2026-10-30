//! The moving point of `Order._eval_subs` for a trigonometric replacement
//! (`O(x**2).subs(x, cos(x))` is `O((x - pi/2)**2, (x, pi/2))`).
//!
//! The reference solves `old = new(d)` with `solveset`, takes the first
//! element of the solution set and its limit at `old -> point`. For
//! `new = a*F(var) + b` with `F` one of `cos`, `sin`, `tan` the solution
//! set is a union of image sets over the integers whose first element is
//! the principal inverse at `n = 0`: `acos(y)`, `asin(y)`, `atan(y)` with
//! `y = (old - b)/a`. From `oo` or `-oo`, `acos` and `asin` move the
//! point to `+-oo*I`; there the reference's `Order(e, (x, +-oo*I))` of an
//! expression in which `x` only occurs as `F(x)` is `O(1)`, the leading
//! term of `F(d + oo*I)` in `d` being free of `d`. Other replacements and
//! other moved points are left to the host solver.

use crate::expr::{Expr, Func, Kind};
use crate::ops::{div, sub};
use crate::series::engine::depends;
use crate::subs::subs1;

fn trig_node(e: &Expr, var: &Expr) -> Option<Expr> {
    if let Kind::Func(Func::Cos | Func::Sin | Func::Tan, a) = e.kind() {
        if a.len() == 1 && a[0] == *var {
            return Some(e.clone());
        }
    }
    e.args().iter().find_map(|a| trig_node(a, var))
}

/// `limit(F^-1((old - b)/a), old, point)`, or `None` when `new` is not of
/// that form or the limit is not decided.
pub fn moving_point(old: &Expr, new: &Expr, var: &Expr, point: &Expr) -> Option<Expr> {
    let infinite = |e: &Expr| {
        e.is_constant(crate::expr::Constant::Infinity)
            || e.is_constant(crate::expr::Constant::NegInfinity)
            || e.is_constant(crate::expr::Constant::ComplexInfinity)
    };
    if !old.is_symbol() || point.is_constant(crate::expr::Constant::ComplexInfinity) {
        return None;
    }
    let node = trig_node(new, var)?;
    let t = Expr::dummy("t", crate::symbol::Assumptions::NONE);
    let r = crate::subs::xreplace(new, &std::iter::once((node.clone(), t.clone())).collect());
    if depends(&r, var) {
        return None;
    }
    let a = crate::diff::diff(&r, &t);
    if a.is_zero() || depends(&a, &t) || depends(&a, old) {
        return None;
    }
    let b = subs1(&r, &t, &Expr::zero());
    if depends(&b, old) {
        return None;
    }
    let y = div(&sub(old, &b), &a);
    let Kind::Func(f, _) = node.kind() else { return None };
    let inv = match f {
        Func::Cos => Func::Acos,
        Func::Sin => Func::Asin,
        _ => Func::Atan,
    };
    let d = crate::functions::func(inv, vec![y]);
    let p = crate::limits::limit(&d, old, point, crate::limits::Dir::Plus).ok()?;
    if crate::assumptions::is_finite(&p) == Some(true) || (infinite(point) && imaginary_infinity(&p)) {
        return Some(p);
    }
    None
}

/// `oo*I` or `-oo*I`.
pub fn imaginary_infinity(p: &Expr) -> bool {
    let oo = Expr::oo();
    let i = Expr::i();
    *p == crate::ops::mul(&[oo.clone(), i.clone()]) || *p == crate::ops::mul(&[Expr::int(-1), oo, i])
}

/// Whether every occurrence of `var` in `e` is a `cos(var)`, `sin(var)` or
/// `tan(var)` node.
fn only_trig_of(e: &Expr, var: &Expr) -> bool {
    if e == var {
        return false;
    }
    if let Kind::Func(Func::Cos | Func::Sin | Func::Tan, a) = e.kind() {
        if a.len() == 1 && a[0] == *var {
            return true;
        }
    }
    e.args().iter().all(|a| only_trig_of(a, var))
}

/// The order built at an imaginary infinite point: `O(1, (var, p))` for an
/// expression of `F(var)` only (`O(x**2, (x, oo)).subs(x, cos(x))` is
/// `O(1, (x, oo*I))`), `None` otherwise.
pub fn order_at_imaginary_infinity(e: &Expr, var: &Expr, p: &Expr) -> Option<Expr> {
    (imaginary_infinity(p) && depends(e, var) && only_trig_of(e, var))
        .then(|| crate::series::order::order(&Expr::one(), &[(var.clone(), p.clone())]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cosine_moves_to_half_pi() {
        let x = Expr::symbol("x");
        let p = moving_point(&x, &crate::functions::func(Func::Cos, vec![x.clone()]), &x, &Expr::zero());
        assert_eq!(p, Some(div(&Expr::pi(), &Expr::int(2))));
    }

    #[test]
    fn cosine_from_infinity_moves_to_imaginary_infinity() {
        let x = Expr::symbol("x");
        let cos = crate::functions::func(Func::Cos, vec![x.clone()]);
        let p = moving_point(&x, &cos, &x, &Expr::oo()).unwrap();
        assert_eq!(p.to_string(), "oo*I");
        let e = crate::ops::pow(&cos, &Expr::int(2));
        assert_eq!(order_at_imaginary_infinity(&e, &x, &p).unwrap().to_string(), "O(1, (x, oo*I))");
        assert!(order_at_imaginary_infinity(&crate::ops::mul(&[x.clone(), cos]), &x, &p).is_none());
    }
}
