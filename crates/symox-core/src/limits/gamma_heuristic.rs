//! `gamma` and `factorial` of an argument tending to `-oo`, where the
//! Gruntz algorithm raises `PoleError` and the reference `heuristics`
//! evaluates the function at the limit of its argument.
//!
//! `Limit.doit` rewrites `factorial(a)` as `gamma(a + 1)`, Gruntz fails on
//! the poles of `gamma` at the negative integers and `heuristics` returns
//! `gamma(limit(a))`, which stays the unevaluated `gamma(-oo)`
//! (`limit(factorial(x), x, -oo)`). The sum and the product rules of
//! `heuristics` give `c0 + c1*gamma(-oo)**k` for `c0 + c1*F(a)**k` with
//! `c0`, `c1` free of the variable and an integer `k != 0`
//! (`1/gamma(-oo)`, `1 + gamma(-oo)`). Other combinations (`gamma(x)/x`,
//! which the reference folds to 0 through `0*gamma(-oo)`) stay
//! unevaluated.

use super::gruntz::Gruntz;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{add, mul2, pow};
use crate::series::engine::depends;
use crate::special::SpecialFn;

/// `F(a)` with `F` in `gamma`, `factorial` and `a -> -oo`.
fn gamma_at_neg_infinity(f: &Expr, t: &Expr) -> bool {
    let a = match f.kind() {
        Kind::Func(Func::Gamma, a) | Kind::Func(Func::Special(SpecialFn::Factorial), a) if a.len() == 1 => &a[0],
        _ => return false,
    };
    let real = !crate::traverse::has(a, &Expr::i())
        && crate::traverse::free_symbols(a).iter().all(|v| v == t || crate::assumptions::is_real(v) == Some(true));
    if !real || !depends(a, t) {
        return false;
    }
    let mut g = Gruntz::new(t);
    matches!(g.limitinf(a), Ok(l) if l.is_constant(Constant::NegInfinity))
}

/// `(c1, k)` for `c1*F(a)**k`.
fn scaled_power(e: &Expr, t: &Expr) -> Option<(Expr, Expr)> {
    let (c, b, k) = match e.kind() {
        Kind::Mul(c, facs) => {
            let mut dep = facs.iter().filter(|(b, p)| depends(b, t) || depends(p, t));
            let (b, k) = dep.next()?;
            if dep.next().is_some() {
                return None;
            }
            let rest: Vec<Expr> = facs
                .iter()
                .filter(|(b2, _)| b2 != b)
                .map(|(b, p)| crate::ops::mul_factor(b, p))
                .chain(std::iter::once(Expr::num(c.clone())))
                .collect();
            (crate::ops::mul(&rest), b.clone(), k.clone())
        }
        Kind::Pow(b, k) => (Expr::one(), b.clone(), k.clone()),
        _ => (Expr::one(), e.clone(), Expr::one()),
    };
    if k.as_i64().is_none_or(|v| v == 0) || !gamma_at_neg_infinity(&b, t) {
        return None;
    }
    Some((c, k))
}

/// The `heuristics` value of `e` (variable `t -> oo`), see the module notes.
pub(super) fn gamma_value(e: &Expr, t: &Expr) -> Option<Expr> {
    let terms = if e.is_add() { e.args() } else { vec![e.clone()] };
    let (dep, free): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(|a| depends(a, t));
    let [d] = dep.as_slice() else { return None };
    let (c, k) = scaled_power(d, t)?;
    let g = crate::functions::func(Func::Gamma, vec![Expr::neg_oo()]);
    let mut parts = free;
    parts.push(mul2(&c, &pow(&g, &k)));
    Some(add(&parts))
}
