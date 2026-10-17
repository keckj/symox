//! Inverse Laplace transform.
//!
//! Rational functions go through partial fractions (first and second
//! order denominators are inverted in closed form), then a table of simple
//! rules (`(s+a)**-c`, `1/(s**2+a**2)**2`, `1/(s*(s+a)**b)`), expansion,
//! and the program rules: `exp(-a*s)` time shift, `s**n` time
//! differentiation and `Derivative(F, s)` frequency differentiation.

use super::util::*;
use crate::expr::{Expr, Func, Kind};
use crate::functions::{cos, cosh, dirac_delta, exp, gamma, heaviside, sin, sinh};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};
use crate::special::SpecialFn;

type R = (Expr, Expr);

/// Inverse Laplace transform of `F` with respect to `s`: `(f(t), cond)`.
/// `plane` is only recorded in held results (`None` prints as `_None`).
pub fn inverse_laplace_transform(f: &Expr, s: &Expr, t: &Expr, plane: Option<&Expr>, simplify: bool) -> R {
    let plane = plane.cloned().unwrap_or_else(super::none_plane);
    let (r, c) = ilt(f, s, t, &plane, true, 0);
    let r = if simplify { crate::simplify::simplify(&r) } else { r };
    (r, c)
}

fn held_ilt(f: &Expr, s: &Expr, t: &Expr, plane: &Expr) -> R {
    (super::held(SpecialFn::InverseLaplaceTransform, vec![f.clone(), s.clone(), t.clone(), plane.clone()]), tru())
}

fn is_rational(e: &Expr, s: &Expr) -> bool {
    if !has(e, s) {
        return true;
    }
    match e.kind() {
        Kind::Sym(_) => true,
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().all(|a| is_rational(a, s)),
        Kind::Pow(b, n) => n.is_integer_number() && is_rational(b, s),
        _ => false,
    }
}

fn has_exp(e: &Expr) -> bool {
    has_head(e, &|k| matches!(k, Kind::Func(Func::Exp, _)))
}

pub(crate) fn ilt(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, dorational: bool, depth: u32) -> R {
    if depth > MAX_DEPTH {
        return held_ilt(f, s, t, plane);
    }
    let mut out = Vec::new();
    let mut conds = Vec::new();
    for term in terms(f) {
        let term = if has_exp(&term) {
            // move time shifts into exp(-a*s) numerator factors
            let ms = neg(s);
            let flipped = crate::subs::subs1(&term, s, &ms);
            let tog = crate::poly::api::together(&flipped);
            crate::subs::subs1(&tog, s, &ms)
        } else {
            term
        };
        let (k, g) = split_mul(&term, s);
        let r = (if dorational && is_rational(&term, s) { rational(&g, s, t, plane, depth) } else { None })
            .or_else(|| simple(&g, s, t))
            .or_else(|| expand_rule(&g, s, t, plane, depth))
            .or_else(|| prog(&g, s, t, plane, depth))
            .unwrap_or_else(|| held_ilt(&g, s, t, plane));
        out.push(mul(&[k, r.0]));
        conds.push(r.1);
    }
    let cond = if conds.len() == 1 { conds.pop().unwrap() } else { crate::logic::and(conds) };
    (add(&out), cond)
}

fn rational(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, depth: u32) -> Option<R> {
    let ap = crate::poly::api::apart(f, Some(s)).unwrap_or_else(|_| f.clone());
    let ht = heaviside(t);
    let mut out = Vec::new();
    let mut conds = Vec::new();
    for term in terms(&ap) {
        let (n, d) = crate::poly::api::as_numer_denom(&term);
        let (Some(dc), Some(nc)) = (poly(&d, s), poly(&n, s)) else { return None };
        let lead = dc.last().unwrap().clone();
        let dc: Vec<Expr> = dc.iter().map(|c| div(c, &lead)).collect();
        let nc: Vec<Expr> = nc.iter().map(|c| div(c, &lead)).collect();
        match dc.len() {
            1 => {
                for (i, c) in nc.iter().enumerate() {
                    out.push(mul(&[c.clone(), dirac_delta_n(t, i)]));
                }
            }
            2 if nc.len() == 1 => {
                out.push(mul(&[ht.clone(), nc[0].clone(), exp(&neg(&mul(&[dc[0].clone(), t.clone()])))]));
            }
            3 if nc.len() <= 2 => {
                let a = div(&dc[1], &int(2));
                let b0 = crate::expand::expand(&sub(&dc[0], &pow(&a, &int(2))));
                let b = crate::poly::api::factor(&b0).unwrap_or(b0);
                let (l, m) = if nc.len() == 1 { (Expr::zero(), nc[0].clone()) } else { (nc[1].clone(), nc[0].clone()) };
                let ea = exp(&neg(&mul(&[a.clone(), t.clone()])));
                let r = if b.is_zero() {
                    mul(&[
                        add(&[mul(&[m, t.clone()]), mul(&[l, sub(&Expr::one(), &mul(&[a.clone(), t.clone()]))])]),
                        ea,
                    ])
                } else {
                    let hyp = negv(&b);
                    let bb = if hyp { neg(&b) } else { b };
                    let bs = crate::simplify::simplify(&sqrt(&bb));
                    let bst = mul(&[bs.clone(), t.clone()]);
                    let (c, sn) = if hyp { (cosh(&bst), sinh(&bst)) } else { (cos(&bst), sin(&bst)) };
                    let ml = sub(&m, &mul(&[a.clone(), l.clone()]));
                    add(&[mul(&[l, ea.clone(), c]), mul(&[div(&ml, &bs), ea, sn])])
                };
                out.push(mul(&[ht.clone(), r]));
            }
            _ => {
                let (r, c) = ilt(&term, s, t, plane, false, depth + 1);
                out.push(r);
                conds.push(c);
            }
        }
    }
    let cond = if conds.is_empty() { tru() } else { crate::logic::and(conds) };
    Some((add(&out), cond))
}

fn dirac_delta_n(t: &Expr, n: usize) -> Expr {
    if n == 0 {
        dirac_delta(t)
    } else {
        crate::functions::func(Func::DiracDelta, vec![t.clone(), int(n as i64)])
    }
}

fn simple(f: &Expr, s: &Expr, t: &Expr) -> Option<R> {
    let ht = heaviside(t);
    if f.is_one() {
        return Some((dirac_delta(t), tru()));
    }
    let done = |e: Expr| Some((mul(&[ht.clone(), e]), tru()));
    // (alpha*s + a)**(-c)
    if let Kind::Pow(b, e) = f.kind() {
        if !has(e, s) {
            if let Some((al, a)) = linear(b, s) {
                let c = neg(e);
                if al.is_one() {
                    return done(div(&mul(&[pow(t, &sub(&c, &Expr::one())), exp(&neg(&mul(&[a, t.clone()])))]), &gamma(&c)));
                }
                if pos(&al) {
                    let sh = div(&a, &al);
                    let r = div(&mul(&[pow(t, &sub(&c, &Expr::one())), exp(&neg(&mul(&[sh, t.clone()])))]), &gamma(&c));
                    return done(mul(&[pow(&al, e), r]));
                }
            }
            // 1/(s**2 + a**2)**2
            if e == &int(-2) {
                if let Some(pc) = poly(b, s) {
                    if pc.len() == 3 && pc[1].is_zero() && pc[2].is_one() {
                        let a = sqrt(&pc[0]);
                        let at = mul(&[a.clone(), t.clone()]);
                        let r = div(&sub(&sin(&at), &mul(&[at.clone(), cos(&at)])), &mul(&[int(2), pow(&a, &int(3))]));
                        return done(r);
                    }
                }
            }
        }
        return None;
    }
    // 1/(s*(s+a)**b)
    if f.is_mul() {
        let fs = factors(f);
        if fs.len() == 2 {
            for (i, j) in [(0usize, 1usize), (1, 0)] {
                if power_of(&fs[i], s) != Some(int(-1)) {
                    continue;
                }
                if let Kind::Pow(base, e) = fs[j].kind() {
                    if has(e, s) {
                        continue;
                    }
                    if let Some((al, a)) = linear(base, s) {
                        if al.is_one() && !a.is_zero() {
                            let b = neg(e);
                            let r = div(
                                &lowergamma(&b, &mul(&[a.clone(), t.clone()])),
                                &mul(&[pow(&a, &b), gamma(&b)]),
                            );
                            return done(r);
                        }
                    }
                }
            }
        }
    }
    None
}

fn expand_rule(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, depth: u32) -> Option<R> {
    if f.is_add() {
        return None;
    }
    let r = crate::expand::expand(f);
    if r.is_add() {
        return Some(ilt(&r, s, t, plane, true, depth + 1));
    }
    if is_rational(f, s) {
        if let Ok(r) = crate::poly::api::apart(f, Some(s)) {
            if r.is_add() {
                return Some(ilt(&r, s, t, plane, true, depth + 1));
            }
        }
    }
    None
}

fn prog(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, depth: u32) -> Option<R> {
    time_shift(f, s, t, plane, depth)
        .or_else(|| freq_shift(f, s, t, plane))
        .or_else(|| time_diff(f, s, t, plane, depth))
        .or_else(|| freq_diff(f, s, t, plane, depth))
}

fn time_shift(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, depth: u32) -> Option<R> {
    if !has(f, s) {
        return Some((mul(&[f.clone(), dirac_delta(t)]), tru()));
    }
    if !has_exp(f) {
        return None;
    }
    let fs = factors(f);
    for (i, e) in fs.iter().enumerate() {
        let Kind::Func(Func::Exp, a) = e.kind() else { continue };
        let Some((a1, a0)) = linear(&a[0], s) else { continue };
        let g = mul(&fs.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, x)| x.clone()).collect::<Vec<_>>());
        if !negv(&a1) {
            return Some(held_ilt(f, s, t, plane));
        }
        let shifted = add(&[t.clone(), a1]);
        let k = exp(&a0);
        if fs.len() == 1 {
            return Some((mul(&[k, dirac_delta(&shifted)]), tru()));
        }
        let (r, c) = ilt(&g, s, &shifted, plane, true, depth + 1);
        return Some((mul(&[k, r]), c));
    }
    None
}

fn freq_shift(f: &Expr, s: &Expr, t: &Expr, plane: &Expr) -> Option<R> {
    let Kind::Func(head @ Func::Undefined(_), args) = f.kind() else { return None };
    if args.len() != 1 {
        return None;
    }
    let (a1, c) = linear(&args[0], s)?;
    if !a1.is_one() || c.is_zero() {
        return None;
    }
    // F(s + c) <-> exp(-c*t) f(t)
    let base = crate::functions::func(head.clone(), vec![s.clone()]);
    let (h, _) = held_ilt(&base, s, t, plane);
    Some((mul(&[exp(&neg(&mul(&[c, t.clone()]))), h]), tru()))
}

fn time_diff(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, depth: u32) -> Option<R> {
    let (n, g) = split_power(f, s);
    let k = n.as_i64().filter(|k| *k > 0)?;
    let (r, c) = ilt(&g, s, t, plane, true, depth + 1);
    let r = crate::subs::subs1(&r, &heaviside(t), &Expr::one());
    let d = crate::diff::diff_n(&r, t, k as u32);
    let held = has_head(&r, &|k| matches!(k, Kind::Func(Func::Special(SpecialFn::InverseLaplaceTransform), _)));
    Some((if held { d } else { mul(&[heaviside(t), d]) }, c))
}

fn freq_diff(f: &Expr, s: &Expr, t: &Expr, plane: &Expr, depth: u32) -> Option<R> {
    let Kind::Derivative(g, vars) = f.kind() else { return None };
    if vars.len() != 1 || &vars[0].0 != s {
        return None;
    }
    let (r, c) = ilt(g, s, t, plane, false, depth + 1);
    Some((mul(&[pow(&neg(t), &int(vars[0].1 as i64)), r]), c))
}
