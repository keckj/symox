//! Antiderivatives that need special functions (`sin(x)/x -> Si(x)`, ...).

use super::*;
use crate::ops::{div, mul2, pow, sqrt};

/// `a` with `u == a*x` (`a` free of `x`, non zero).
fn linear_through_origin(u: &Expr, x: &Expr) -> Option<Expr> {
    if u == x {
        return Some(Expr::one());
    }
    let a = div(u, x);
    if crate::traverse::has(&a, x) || a.is_zero() {
        return None;
    }
    Some(a)
}

/// Integrate `f` with respect to `x` when it is one of the classic
/// non-elementary integrands.
pub fn table(f: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Mul(c, facs) = f.kind() else {
        // 1/log(a x)
        if let Kind::Pow(b, e) = f.kind() {
            if e.is_minus_one() {
                if let Kind::Func(Func::Log, la) = b.kind() {
                    let a = linear_through_origin(&la[0], x)?;
                    return Some(div(&super::li(&la[0]), &a));
                }
            }
        }
        return exp_square(f, x);
    };
    if !c.is_one() || facs.len() != 2 {
        return None;
    }
    // g(u) * x^-1
    let (g, xp) = if facs[1].0 == *x && facs[1].1.is_minus_one() {
        (crate::ops::mul_factor(&facs[0].0, &facs[0].1), true)
    } else if facs[0].0 == *x && facs[0].1.is_minus_one() {
        (crate::ops::mul_factor(&facs[1].0, &facs[1].1), true)
    } else {
        (Expr::zero(), false)
    };
    if !xp {
        return None;
    }
    let Kind::Func(h, a) = g.kind() else { return None };
    let u = &a[0];
    linear_through_origin(u, x)?;
    let head = match h {
        Func::Sin => SpecialFn::Si,
        Func::Cos => SpecialFn::Ci,
        Func::Sinh => SpecialFn::Shi,
        Func::Cosh => SpecialFn::Chi,
        Func::Exp => SpecialFn::Ei,
        _ => return None,
    };
    Some(sf(head, vec![u.clone()]))
}

/// `exp(c x^2)` with `c > 0` -> `sqrt(pi) erfi(sqrt(c) x) / (2 sqrt(c))`.
fn exp_square(f: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Func(Func::Exp, a) = f.kind() else { return None };
    let c = div(&a[0], &x.powi(2));
    if crate::traverse::has(&c, x) {
        return None;
    }
    if crate::evalf::eval_f64(&c).is_some_and(|v| v > 0.0) {
        let sc = sqrt(&c);
        return Some(div(
            &mul2(&sqrt(&Expr::pi()), &super::erfi(&mul2(&sc, x))),
            &mul2(&Expr::int(2), &sc),
        ));
    }
    let _ = pow;
    None
}
