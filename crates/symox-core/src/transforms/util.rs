//! Small helpers shared by the transform engines.

use crate::expr::{Expr, Func, Kind};
use crate::functions::exp;
use crate::ops::{add, div, mul, neg, sub};
use crate::special::SpecialFn;

/// Recursion budget of the rule engines (rules may call each other).
pub const MAX_DEPTH: u32 = 14;

/// Does `e` depend on `x` (structurally)?
pub fn has(e: &Expr, x: &Expr) -> bool {
    crate::binders::depends_on(e, x)
}

/// Does `e` contain the sub-expression `p`?
pub fn contains(e: &Expr, p: &Expr) -> bool {
    crate::traverse::has(e, p)
}

pub fn terms(e: &Expr) -> Vec<Expr> {
    if e.is_add() { e.args() } else { vec![e.clone()] }
}

pub fn factors(e: &Expr) -> Vec<Expr> {
    if e.is_mul() { e.args() } else { vec![e.clone()] }
}

/// `e == k*f` with `k` free of `x` (multiplicative split).
pub fn split_mul(e: &Expr, x: &Expr) -> (Expr, Expr) {
    let mut k = Vec::new();
    let mut f = Vec::new();
    for a in factors(e) {
        if has(&a, x) { f.push(a) } else { k.push(a) }
    }
    (mul(&k), mul(&f))
}

/// Dense coefficients `[c0, c1, ...]` of a polynomial in `x`.
pub fn poly(e: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    let mut c = crate::integrate::poly_coeffs(e, x)?;
    while c.len() > 1 && c.last().is_some_and(|z| z.is_zero()) {
        c.pop();
    }
    Some(c)
}

/// `u == a*x + b` with `a != 0` free of `x`: returns `(a, b)`.
pub fn linear(u: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    if !has(u, x) {
        return None;
    }
    let c = poly(u, x)?;
    if c.len() != 2 || c[1].is_zero() {
        return None;
    }
    Some((c[1].clone(), c[0].clone()))
}

/// `u == a*x` exactly (no offset).
pub fn scaled(u: &Expr, x: &Expr) -> Option<Expr> {
    let (a, b) = linear(u, x)?;
    if b.is_zero() { Some(a) } else { None }
}

/// `e == x**n` with `n` free of `x`, `x` itself giving `n = 1`.
pub fn power_of(e: &Expr, x: &Expr) -> Option<Expr> {
    if e == x {
        return Some(Expr::one());
    }
    match e.kind() {
        Kind::Pow(b, n) if b == x && !has(n, x) => Some(n.clone()),
        _ => None,
    }
}

/// Split `f` into `x**n * rest` where `rest` has no plain power of `x`.
pub fn split_power(f: &Expr, x: &Expr) -> (Expr, Expr) {
    let mut n = Vec::new();
    let mut rest = Vec::new();
    for a in factors(f) {
        match power_of(&a, x) {
            Some(p) => n.push(p),
            None => rest.push(a),
        }
    }
    (add(&n), mul(&rest))
}

pub fn truth(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Bool(b) => Some(*b),
        _ => None,
    }
}

pub fn is_true(e: &Expr) -> bool {
    truth(e) == Some(true)
}

pub fn tru() -> Expr {
    Expr::boolean(true)
}

pub fn pos(e: &Expr) -> bool {
    crate::assumptions::is_positive(e) == Some(true)
}
pub fn negv(e: &Expr) -> bool {
    crate::assumptions::is_negative(e) == Some(true)
}
pub fn nonneg(e: &Expr) -> bool {
    crate::assumptions::is_nonnegative(e) == Some(true)
}
pub fn nonpos(e: &Expr) -> bool {
    crate::assumptions::is_nonpositive(e) == Some(true)
}
pub fn real(e: &Expr) -> bool {
    crate::assumptions::is_real(e) == Some(true)
}

/// Is `e` a positive real number (or `re(e) > 0`)? Used for conditions
/// such as `re(a) > 0`.
pub fn re_pos(e: &Expr) -> bool {
    pos(e) || pos(&crate::functions::re(e))
}
pub fn re_nonneg(e: &Expr) -> bool {
    nonneg(e) || nonneg(&crate::functions::re(e))
}

/// `re(n) > c` decided or as a relational condition.
pub fn re_gt(n: &Expr, c: &Expr) -> Expr {
    let r = crate::relational::gt(&crate::functions::re(n), c);
    if truth(&r).is_some() {
        return r;
    }
    let r2 = crate::relational::gt(n, c);
    if truth(&r2).is_some() { r2 } else { r }
}

pub fn special(s: SpecialFn, args: Vec<Expr>) -> Expr {
    crate::special::apply(s, args)
}
pub fn besselj(n: &Expr, x: &Expr) -> Expr {
    special(SpecialFn::BesselJ, vec![n.clone(), x.clone()])
}
pub fn besselk(n: &Expr, x: &Expr) -> Expr {
    special(SpecialFn::BesselK, vec![n.clone(), x.clone()])
}
pub fn uppergamma(a: &Expr, x: &Expr) -> Expr {
    special(SpecialFn::UpperGamma, vec![a.clone(), x.clone()])
}
pub fn lowergamma(a: &Expr, x: &Expr) -> Expr {
    special(SpecialFn::LowerGamma, vec![a.clone(), x.clone()])
}
pub fn ei(x: &Expr) -> Expr {
    special(SpecialFn::Ei, vec![x.clone()])
}
pub fn digamma(x: &Expr) -> Expr {
    special(SpecialFn::Digamma, vec![x.clone()])
}

pub fn int(i: i64) -> Expr {
    Expr::int(i)
}
pub fn rat(n: i64, d: i64) -> Expr {
    Expr::rational(n, d)
}

/// Rewrite `sin, cos, sinh, cosh` (and `tan`-free products) in terms of
/// exponentials.
pub fn rewrite_exp(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| rewrite_exp(c));
    let i = Expr::i();
    match e.kind() {
        Kind::Func(Func::Sin, a) => {
            let u = mul(&[i.clone(), a[0].clone()]);
            div(&sub(&exp(&u), &exp(&neg(&u))), &mul(&[int(2), i]))
        }
        Kind::Func(Func::Cos, a) => {
            let u = mul(&[i, a[0].clone()]);
            div(&add(&[exp(&u), exp(&neg(&u))]), &int(2))
        }
        Kind::Func(Func::Sinh, a) => div(&sub(&exp(&a[0]), &exp(&neg(&a[0]))), &int(2)),
        Kind::Func(Func::Cosh, a) => div(&add(&[exp(&a[0]), exp(&neg(&a[0]))]), &int(2)),
        _ => e,
    }
}

/// Does `e` contain one of the trigonometric/hyperbolic heads?
pub fn has_trig(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |n| {
        if matches!(n.kind(), Kind::Func(Func::Sin | Func::Cos | Func::Sinh | Func::Cosh, _)) {
            found = true;
        }
    });
    found
}

pub fn has_head(e: &Expr, pred: &dyn Fn(&Kind) -> bool) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |n| {
        if pred(n.kind()) {
            found = true;
        }
    });
    found
}

/// Write `term` as `p * exp(c1*x + c0)` with `p` free of `x`;
/// returns `(p*exp(c0), c1)`.
pub fn as_exp_term(term: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    let mut expo = Vec::new();
    let mut rest = Vec::new();
    for f in factors(term) {
        match f.kind() {
            Kind::Func(Func::Exp, a) => expo.push(a[0].clone()),
            Kind::Pow(b, n) if matches!(b.kind(), Kind::Func(Func::Exp, _)) && !has(n, x) => {
                expo.push(mul(&[n.clone(), b.args()[0].clone()]))
            }
            Kind::Pow(b, n) if b.is_constant(crate::expr::Constant::E) => expo.push(n.clone()),
            _ => rest.push(f),
        }
    }
    let p = mul(&rest);
    if has(&p, x) {
        return None;
    }
    let u = crate::expand::expand(&add(&expo));
    let c = poly(&u, x)?;
    match c.len() {
        1 => Some((mul(&[p, exp(&c[0])]), Expr::zero())),
        2 => Some((mul(&[p, exp(&c[0])]), c[1].clone())),
        _ => None,
    }
}

/// Number of nodes of `e` (a crude complexity measure).
pub fn size(e: &Expr) -> usize {
    let mut n = 0;
    crate::traverse::preorder(e, &mut |_| n += 1);
    n
}
