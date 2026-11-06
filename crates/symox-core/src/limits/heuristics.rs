//! The reference `heuristics(e, z, z0, dir)` of `sympy.series.limits`:
//! the limit of an `Add`, `Mul`, `Pow` or function from the limits of its
//! arguments (`sin(1/x) + atan(x)` at `0+` is `AccumBounds(-1, 1)`), or of
//! `e(1/z)` at `0+` for `z0 = oo`. `None` is the reference `None` (a held
//! or undefined argument limit, or a `nan` that the rewritings do not
//! remove).

use super::{limit, Dir};
use crate::expr::{Constant, Expr, Kind};
use crate::ops::pow;

fn has_infinity(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(c) => *c == Constant::Infinity,
        Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(has_infinity),
    }
}

fn is_nan(e: &Expr) -> bool {
    e.is_constant(Constant::NaN)
}

/// The reference `limit(a, z, z0, dir)` that `heuristics` calls: the
/// native limit, then the accumulation bounds of an oscillating
/// expression (`sin(1/x)` at `0+` is `AccumBounds(-1, 1)`).
fn full_limit(a: &Expr, z: &Expr, z0: &Expr, dir: Dir) -> Option<Expr> {
    if let Ok(l) = limit(a, z, z0, dir) {
        return Some(l);
    }
    let d = match dir {
        Dir::Plus => "+",
        Dir::Minus => "-",
        Dir::Both => "+-",
    };
    super::accum::oscillating_limit(a, z, z0, d).ok().flatten()
}

/// `e.func(*r)` when an argument limit is an accumulation bounds: the
/// bounds arithmetic of `Add`, `Mul` and `Pow`, and a function applied to
/// the bounds (`AccumBounds` evaluation); an arithmetic error is `None`.
fn rebuild_bounds(e: &Expr, r: &[Expr]) -> Option<Expr> {
    use super::accum;
    match e.kind() {
        Kind::Add(..) => accum::add_ctor(r).ok(),
        Kind::Mul(..) => accum::mul_ctor(r).ok(),
        Kind::Pow(..) => accum::op_pow(&r[0], &r[1]).ok(),
        Kind::Func(f, _) if r.len() == 1 => accum::apply(f.clone(), &r[0]).ok(),
        _ => None,
    }
}

/// `heuristics(e, z, z0, dir)`.
pub fn heuristics(e: &Expr, z: &Expr, z0: &Expr, dir: Dir) -> Option<Expr> {
    if z0.is_constant(Constant::Infinity) {
        // limit(e.subs(z, 1/z), z, 0, '+'); a held limit gives None
        let sub = crate::subs::subs1(e, z, &pow(z, &Expr::minus_one()));
        return limit(&sub, z, &Expr::zero(), Dir::Plus).ok();
    }
    if !matches!(e.kind(), Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Func(..)) {
        return None;
    }
    let args = e.args();
    let mut r = Vec::with_capacity(args.len());
    for a in &args {
        let l = full_limit(a, z, z0, dir)?;
        if has_infinity(&l) && crate::assumptions::is_finite(&l).is_none() {
            if e.is_add() {
                // factor_terms, then together, then factor, until a Mul
                let mut m = crate::simptools::exprtools::factor_terms(e, false, false, false, true);
                if !m.is_mul() {
                    m = crate::simplify::together(&m);
                }
                if !m.is_mul() {
                    m = crate::poly::api::factor(e).ok()?;
                }
                if m.is_mul() {
                    return heuristics(&m, z, z0, dir);
                }
            }
            return None;
        }
        if is_nan(&l) {
            return None;
        }
        r.push(l);
    }
    if r.is_empty() {
        return None;
    }
    let mut rv = if r.iter().any(super::accum::is_bounds) { rebuild_bounds(e, &r)? } else { crate::simptools::rebuild(e, &r) };
    if is_nan(&rv) && e.is_mul() && r.iter().any(super::accum::is_bounds) {
        // the bounded factors times the limit of the simplified others
        let (r2, e2): (Vec<_>, Vec<_>) = r.iter().zip(&args).partition(|(l, _)| super::accum::is_bounds(l));
        let r2: Vec<Expr> = r2.into_iter().map(|(l, _)| l.clone()).collect();
        let e2: Vec<Expr> = e2.into_iter().map(|(_, a)| a.clone()).collect();
        if !e2.is_empty() {
            let e3 = crate::simplify::simplify(&crate::ops::mul(&e2));
            let l = full_limit(&e3, z, z0, dir)?;
            let mut v = l;
            for b in &r2 {
                v = super::accum::op_mul(&v, b).ok()?;
            }
            rv = v;
        }
    }
    if is_nan(&rv) {
        let rat_e = crate::simptools::ratsimp::ratsimp(e);
        if is_nan(&rat_e) || rat_e == *e {
            return None;
        }
        return limit(&rat_e, z, z0, dir).ok();
    }
    Some(rv)
}
