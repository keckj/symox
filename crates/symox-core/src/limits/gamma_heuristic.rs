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

/// `loggamma(a)` with `a -> -oo` along the reals.
fn loggamma_at_neg_infinity(f: &Expr, t: &Expr) -> bool {
    let Kind::Func(Func::LogGamma, a) = f.kind() else { return false };
    let a = &a[0];
    let real = !crate::traverse::has(a, &Expr::i())
        && crate::traverse::free_symbols(a).iter().all(|v| v == t || crate::assumptions::is_real(v) == Some(true));
    if !real || !depends(a, t) {
        return false;
    }
    let mut g = Gruntz::new(t);
    matches!(g.limitinf(a), Ok(l) if l.is_constant(Constant::NegInfinity))
}

/// The reference `heuristics` on `loggamma(a)`, `a -> -oo`: Gruntz raises
/// `PoleError` (the asymptotic expansion of `loggamma` holds for a
/// positive argument only), and `heuristics` evaluates each argument of the
/// products, sums and powers at its limit, with `loggamma(-oo) = zoo`
/// (`loggamma.eval`), then rebuilds the expression: `x*loggamma(x)` and
/// `loggamma(x) + 3` give `zoo`, `1/loggamma(x)` gives 0, while
/// `loggamma(x)/x` (`0*zoo`), `loggamma(x) + x` (`zoo - oo`) and
/// `loggamma(x) + loggamma(x + 1)` give `nan` and stay unevaluated. Only
/// products, sums and powers with a constant exponent around the
/// `loggamma` nodes are evaluated; other functions of them decline.
pub(super) fn loggamma_value(e: &Expr, t: &Expr) -> Option<Expr> {
    fn has_pole(e: &Expr, t: &Expr) -> bool {
        loggamma_at_neg_infinity(e, t) || e.args().iter().any(|a| has_pole(a, t))
    }
    fn value(e: &Expr, t: &Expr) -> Option<Expr> {
        if !depends(e, t) {
            return Some(e.clone());
        }
        if loggamma_at_neg_infinity(e, t) {
            return Some(Expr::zoo());
        }
        if !has_pole(e, t) {
            // `limit(a)` of an argument free of the pole; an infinite value
            // inside a larger expression (`oo*sign(y)`) is the reference's
            // `l.has(oo) and l.is_finite is None` exit.
            let mut g = Gruntz::new(t);
            let l = g.limitinf(e).ok()?;
            if depends(&l, t) || l.is_constant(Constant::NaN) {
                return None;
            }
            let atomic_inf = l.is_constant(Constant::Infinity)
                || l.is_constant(Constant::NegInfinity)
                || l.is_constant(Constant::ComplexInfinity);
            if !atomic_inf && has_inf(&l) {
                return None;
            }
            return Some(l);
        }
        match e.kind() {
            Kind::Add(..) => {
                let v: Option<Vec<Expr>> = e.args().iter().map(|a| value(a, t)).collect();
                Some(add(&v?))
            }
            Kind::Mul(..) => {
                let v: Option<Vec<Expr>> = e.args().iter().map(|a| value(a, t)).collect();
                Some(crate::ops::mul(&v?))
            }
            Kind::Pow(b, k) if !depends(k, t) => Some(pow(&value(b, t)?, k)),
            _ => None,
        }
    }
    fn has_inf(e: &Expr) -> bool {
        matches!(e.kind(), Kind::Const(c) if c.is_infinite()) || e.args().iter().any(has_inf)
    }
    if !has_pole(e, t) {
        return None;
    }
    let r = value(e, t)?;
    (!depends(&r, t) && !r.is_constant(Constant::NaN)).then_some(r)
}
