//! Forward Laplace transform `F(s) = int_0^oo f(t) exp(-s*t) dt`.
//!
//! Mirrors SymPy's rule engine: linearity, a table of simple rules
//! (Bateman, Tables of Integral Transforms, ch. 4), then the program
//! rules (time shift by Heaviside, delta masking, time scaling, frequency
//! shift, trigonometric pole pairing, time differentiation, frequency
//! differentiation), then expansion. Unknown terms stay held as
//! `LaplaceTransform(f, t, s)` with plane `-oo`.

use super::util::*;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{
    abs, asinh, atan, dirac_delta, erfc, exp, gamma, heaviside, im, log, max, re,
};
use crate::ops::{add, div, inv, mul, neg, pow, sqrt, sub};
use crate::special::SpecialFn;

type R = (Expr, Expr, Expr);

/// Laplace transform of `f` with respect to `t`: `(F(s), plane, cond)`.
pub fn laplace_transform(f: &Expr, t: &Expr, s: &Expr, simplify: bool) -> R {
    let (r, p, c) = lt(f, t, s, 0);
    let r = if simplify { crate::simplify::simplify(&r) } else { r };
    (r, p, c)
}

fn held_lt(f: &Expr, t: &Expr, s: &Expr) -> R {
    (super::held(SpecialFn::LaplaceTransform, vec![f.clone(), t.clone(), s.clone()]), Expr::neg_oo(), tru())
}

fn ok(f: Expr, plane: Expr) -> Option<R> {
    Some((f, plane, tru()))
}

pub(crate) fn lt(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> R {
    if depth > MAX_DEPTH {
        return held_lt(f, t, s);
    }
    let mut pieces: Vec<(Expr, Expr)> = Vec::new();
    for ff in terms(f) {
        let (k, ft) = split_mul(&ff, t);
        if matches!(ft.kind(), Kind::Piecewise(_)) && !contains(&ft, &dirac_delta(t)) {
            if let Some(h) = piecewise_to_heaviside(&ft, t) {
                for term in terms(&crate::expand::expand(&h)) {
                    let (k1, f1) = split_mul(&term, t);
                    pieces.push((mul(&[k.clone(), k1]), f1));
                }
                continue;
            }
        }
        pieces.push((k, ft));
    }
    let mut out = Vec::new();
    let mut planes = Vec::new();
    let mut conds = Vec::new();
    for (k, ft) in pieces {
        let h0 = heaviside(t);
        let ft = if contains(&ft, &h0) && !contains(&ft, &dirac_delta(t)) {
            crate::subs::subs1(&ft, &h0, &Expr::one())
        } else {
            ft
        };
        let r = simple(&ft, t, s)
            .or_else(|| prog(&ft, t, s, depth))
            .or_else(|| expand_rule(&ft, t, s, depth))
            .unwrap_or_else(|| held_lt(&ft, t, s));
        out.push(mul(&[k, r.0]));
        planes.push(r.1);
        conds.push(r.2);
    }
    let plane = if planes.len() == 1 { planes.pop().unwrap() } else { max(planes) };
    let cond = if conds.len() == 1 { conds.pop().unwrap() } else { crate::logic::and(conds) };
    (add(&out), plane, cond)
}

/// Convert a Piecewise in `t` into a sum with Heaviside factors (valid in
/// the Laplace sense: single points do not matter).
pub(crate) fn piecewise_to_heaviside(f: &Expr, t: &Expr) -> Option<Expr> {
    let Kind::Piecewise(pieces) = f.kind() else { return None };
    let mut prev: Vec<Expr> = Vec::new();
    let mut out = Vec::new();
    for (e, c) in pieces {
        let hc = cond_heaviside(c, t)?;
        let mut fac = vec![e.clone(), hc.clone()];
        for p in &prev {
            fac.push(sub(&Expr::one(), p));
        }
        out.push(mul(&fac));
        if is_true(c) {
            break;
        }
        prev.push(hc);
    }
    Some(add(&out))
}

fn cond_heaviside(c: &Expr, t: &Expr) -> Option<Expr> {
    use crate::expr::RelOp;
    match c.kind() {
        Kind::Bool(true) => Some(Expr::one()),
        Kind::Bool(false) => Some(Expr::zero()),
        Kind::Rel(op, l, r) => {
            let d = sub(l, r);
            if !has(&d, t) {
                return None;
            }
            match op {
                RelOp::Gt | RelOp::Ge => Some(heaviside(&d)),
                RelOp::Lt | RelOp::Le => Some(heaviside(&neg(&d))),
                _ => None,
            }
        }
        Kind::Func(Func::Logic(crate::logic::LogicOp::And), args) => {
            let mut v = Vec::new();
            for a in args {
                v.push(cond_heaviside(a, t)?);
            }
            Some(mul(&v))
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------
// simple rules
// ---------------------------------------------------------------------

fn s2pw(s: &Expr, w: &Expr) -> Expr {
    add(&[pow(s, &int(2)), pow(w, &int(2))])
}
fn s2mw(s: &Expr, w: &Expr) -> Expr {
    sub(&pow(s, &int(2)), &pow(w, &int(2)))
}
fn absim(w: &Expr) -> Expr {
    abs(&im(w))
}
fn absre(w: &Expr) -> Expr {
    abs(&re(w))
}

fn simple(f: &Expr, t: &Expr, s: &Expr) -> Option<R> {
    if !has(f, t) {
        return ok(div(f, s), Expr::zero());
    }
    if f == t {
        return ok(inv(&pow(s, &int(2))), Expr::zero());
    }
    match f.kind() {
        Kind::Pow(b, n) if !has(n, t) => simple_pow(b, n, t, s),
        Kind::Func(head, args) => simple_func(head, args, t, s),
        Kind::Mul(..) => simple_mul(f, t, s),
        _ => None,
    }
}

/// `t**n` rule; open conditions are returned (SymPy's integration
/// fallback gives `re(n) > -1`).
fn power_rule(n: &Expr, s: &Expr) -> Option<R> {
    let c = re_gt(n, &int(-1));
    if truth(&c) == Some(false) {
        return None;
    }
    let np1 = add(&[n.clone(), Expr::one()]);
    Some((div(&gamma(&np1), &pow(s, &np1)), Expr::zero(), c))
}

fn simple_pow(b: &Expr, n: &Expr, t: &Expr, s: &Expr) -> Option<R> {
    if b == t {
        return power_rule(n, s);
    }
    let z = Expr::zero();
    if let Some((a, c)) = linear(b, t) {
        if c.is_zero() {
            if pos(&a) {
                let (r, p, cc) = power_rule(n, s)?;
                return Some((mul(&[pow(&a, n), r]), p, cc));
            }
            return None;
        }
        // (a*t + c)**n
        let ba = div(&c, &a);
        let bas = mul(&[ba.clone(), s.clone()]);
        if !(pos(&a) && pos(&c)) && !(pos(&ba) && n == &int(-1)) {
            return None;
        }
        if n == &int(-1) {
            return ok(neg(&div(&mul(&[exp(&bas), ei(&neg(&bas))]), &a)), z);
        }
        if n == &rat(-1, 2) {
            let pi = Expr::pi();
            return ok(
                div(&mul(&[sqrt(&div(&mul(&[a.clone(), pi]), s)), exp(&bas), erfc(&sqrt(&bas))]), &a),
                z,
            );
        }
        if n == &rat(-3, 2) {
            let pi = Expr::pi();
            let r = sub(
                &mul(&[int(2), pow(&c, &rat(-1, 2))]),
                &mul(&[int(2), sqrt(&div(&mul(&[pi, s.clone()]), &a)), exp(&bas), erfc(&sqrt(&bas))]),
            );
            return ok(div(&r, &a), z);
        }
        if n.as_i64().is_some_and(|k| k > 0) {
            return None; // expansion gives the polynomial table
        }
        if truth(&crate::relational::gt(n, &int(-1))) != Some(true) {
            return None;
        }
        let np1 = add(&[n.clone(), Expr::one()]);
        return ok(
            div(&mul(&[pow(&a, n), exp(&bas), uppergamma(&np1, &bas)]), &pow(s, &np1)),
            z,
        );
    }
    if n == &int(2) {
        if let Kind::Func(h, args) = b.kind() {
            if args.len() == 1 {
                if let Some(w) = scaled(&args[0], t) {
                    let w2 = pow(&w, &int(2));
                    match h {
                        Func::Cos => {
                            return ok(
                                div(&add(&[pow(s, &int(2)), mul(&[int(2), w2.clone()])]), &mul(&[add(&[pow(s, &int(2)), mul(&[int(4), w2])]), s.clone()])),
                                mul(&[int(2), absim(&w)]),
                            )
                        }
                        Func::Sinh => {
                            return ok(
                                div(&mul(&[int(2), w2.clone()]), &sub(&pow(s, &int(3)), &mul(&[int(4), w2, s.clone()]))),
                                mul(&[int(2), absre(&w)]),
                            )
                        }
                        Func::Cosh => {
                            return ok(
                                div(&sub(&pow(s, &int(2)), &mul(&[int(2), w2.clone()])), &sub(&pow(s, &int(3)), &mul(&[int(4), w2, s.clone()]))),
                                mul(&[int(2), absre(&w)]),
                            )
                        }
                        Func::Log if pos(&w) => {
                            let l = log(&div(&mul(&[exp(&Expr::constant(Constant::EulerGamma)), s.clone()]), &w));
                            let pi2 = pow(&Expr::pi(), &int(2));
                            return ok(div(&add(&[pow(&l, &int(2)), div(&pi2, &int(6))]), s), z);
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    None
}

fn simple_func(head: &Func, args: &[Expr], t: &Expr, s: &Expr) -> Option<R> {
    let z = Expr::zero();
    match head {
        Func::DiracDelta => {
            let (a, c) = linear(&args[0], t)?;
            let order = match args.get(1) {
                None => 0,
                Some(k) => k.as_i64().filter(|k| *k >= 0)?,
            };
            let b = neg(&c);
            if (pos(&a) && nonneg(&b)) || (negv(&a) && nonpos(&b)) {
                let n = int(order);
                let r = div(
                    &mul(&[pow(s, &n), exp(&neg(&div(&mul(&[s.clone(), b]), &a)))]),
                    &mul(&[pow(&a, &n), abs(&a)]),
                );
                return ok(r, Expr::neg_oo());
            }
            if (negv(&a) && nonneg(&b)) || (pos(&a) && nonpos(&b)) {
                return ok(z, Expr::neg_oo());
            }
            None
        }
        Func::Heaviside => {
            let (a, c) = linear(&args[0], t)?;
            let b = neg(&c);
            let e = exp(&neg(&div(&mul(&[s.clone(), b.clone()]), &a)));
            if pos(&a) && pos(&b) {
                return ok(div(&e, s), z);
            }
            if negv(&a) && negv(&b) {
                return ok(div(&sub(&Expr::one(), &e), s), z);
            }
            if pos(&a) && nonpos(&b) {
                return ok(inv(s), z);
            }
            if negv(&a) && pos(&b) {
                return ok(z.clone(), z);
            }
            None
        }
        Func::Exp => {
            let c = poly(&args[0], t);
            if let Some(c) = c {
                if c.len() == 2 {
                    return ok(div(&exp(&c[0]), &sub(s, &c[1])), re(&c[1]));
                }
                if c.len() == 3 {
                    // exp(-a*t**2 + b*t + c0) = exp(c0) * exp(b*t) * exp(-a*t**2)
                    let a = neg(&c[2]);
                    if !re_pos(&a) {
                        return None;
                    }
                    let sh = sub(s, &c[1]);
                    let four_a = mul(&[int(4), a.clone()]);
                    let g = mul(&[
                        sqrt(&div(&Expr::pi(), &four_a)),
                        exp(&div(&pow(&sh, &int(2)), &four_a)),
                        erfc(&div(&sh, &sqrt(&four_a))),
                    ]);
                    return ok(mul(&[exp(&c[0]), g]), z);
                }
                return None;
            }
            // exp(-a/t)
            let at = crate::expand::expand(&mul(&[args[0].clone(), t.clone()]));
            if !has(&at, t) {
                let a = neg(&at);
                if re_nonneg(&a) {
                    let as_ = mul(&[a.clone(), s.clone()]);
                    return ok(mul(&[int(2), sqrt(&div(&a, s)), besselk(&Expr::one(), &mul(&[int(2), sqrt(&as_)]))]), z);
                }
            }
            None
        }
        Func::Log => {
            let (a, c) = linear(&args[0], t)?;
            if c.is_zero() {
                if !pos(&a) {
                    return None;
                }
                let l = log(&div(&mul(&[exp(&Expr::constant(Constant::EulerGamma)), s.clone()]), &a));
                return ok(neg(&div(&l, s)), z);
            }
            // log(a*t + c) = log(c) + log(1 + a*t/c)
            if !(pos(&a) && pos(&c)) {
                return None;
            }
            let q = div(&mul(&[s.clone(), c.clone()]), &a);
            let r = div(&sub(&log(&c), &mul(&[exp(&q), ei(&neg(&q))])), s);
            ok(r, z)
        }
        Func::Sin => {
            let w = scaled(&args[0], t)?;
            ok(div(&w, &s2pw(s, &w)), absim(&w))
        }
        Func::Cos => {
            let w = scaled(&args[0], t)?;
            ok(div(s, &s2pw(s, &w)), absim(&w))
        }
        Func::Sinh => {
            let w = scaled(&args[0], t)?;
            ok(div(&w, &s2mw(s, &w)), absre(&w))
        }
        Func::Cosh => {
            let w = scaled(&args[0], t)?;
            ok(div(s, &s2mw(s, &w)), absre(&w))
        }
        Func::Abs => {
            if let Kind::Func(Func::Sin, a) = args[0].kind() {
                let w = scaled(&a[0], t)?;
                if pos(&w) {
                    let r = mul(&[
                        div(&w, &s2pw(s, &w)),
                        crate::functions::coth(&div(&mul(&[Expr::pi(), s.clone()]), &mul(&[int(2), w.clone()]))),
                    ]);
                    return ok(r, z);
                }
            }
            None
        }
        Func::Erf => {
            let a = scaled(&args[0], t)?;
            if !pos(&a) {
                return None;
            }
            let two_a = mul(&[int(2), a]);
            ok(div(&mul(&[exp(&pow(&div(s, &two_a), &int(2))), erfc(&div(s, &two_a))]), s), z)
        }
        Func::Special(sf) => {
            let (nu, arg) = (&args[0], args.get(1)?);
            let a = scaled(arg, t)?;
            if has(nu, t) {
                return None;
            }
            match sf {
                SpecialFn::BesselJ if is_true(&re_gt(nu, &int(-1))) => {
                    let r = sqrt(&s2pw(s, &a));
                    ok(div(&pow(&a, nu), &mul(&[r.clone(), pow(&add(&[s.clone(), r]), nu)])), absim(&a))
                }
                SpecialFn::BesselI if is_true(&re_gt(nu, &int(-1))) => {
                    let r = sqrt(&s2mw(s, &a));
                    ok(div(&pow(&a, nu), &mul(&[r.clone(), pow(&add(&[s.clone(), r]), nu)])), absre(&a))
                }
                SpecialFn::BesselY if nu.is_zero() => {
                    let r = mul(&[rat(-2, 1), inv(&Expr::pi()), asinh(&div(s, &a)), inv(&sqrt(&s2pw(s, &a)))]);
                    ok(r, absim(&a))
                }
                SpecialFn::BesselK if nu.is_zero() => {
                    let q = sqrt(&s2mw(s, &a));
                    ok(div(&log(&div(&add(&[s.clone(), q.clone()]), &a)), &q), neg(&re(&a)))
                }
                _ => None,
            }
        }
        _ => None,
    }
}

fn simple_mul(f: &Expr, t: &Expr, s: &Expr) -> Option<R> {
    let z = Expr::zero();
    let (n, g) = split_power(f, t);
    if n.is_zero() {
        // products of two trigonometric functions
        let fs = factors(&g);
        if fs.len() == 2 {
            let pick = |e: &Expr| -> Option<(Func, Expr)> {
                if let Kind::Func(h @ (Func::Sin | Func::Cos), a) = e.kind() {
                    return Some((h.clone(), scaled(&a[0], t)?));
                }
                None
            };
            let (p, q) = (pick(&fs[0])?, pick(&fs[1])?);
            let (p, q) = if p.0 == Func::Sin && q.0 == Func::Cos { (q, p) } else { (p, q) };
            let (a, b) = (p.1, q.1);
            let s2 = pow(s, &int(2));
            let d = mul(&[
                inv(&add(&[s2.clone(), pow(&add(&[a.clone(), b.clone()]), &int(2))])),
                inv(&add(&[s2.clone(), pow(&sub(&a, &b), &int(2))])),
            ]);
            let plane = add(&[absim(&a), absim(&b)]);
            let (a2, b2) = (pow(&a, &int(2)), pow(&b, &int(2)));
            let num = match (&p.0, &q.0) {
                (Func::Sin, Func::Sin) => mul(&[int(2), a.clone(), b.clone(), s.clone()]),
                (Func::Cos, Func::Sin) => mul(&[b.clone(), add(&[s2.clone(), neg(&a2), b2])]),
                (Func::Cos, Func::Cos) => mul(&[s.clone(), add(&[s2, a2, b2])]),
                _ => return None,
            };
            return ok(mul(&[num, d]), plane);
        }
        return None;
    }
    if g.is_mul() {
        return None;
    }
    let one = Expr::one();
    let np1 = add(&[n.clone(), one.clone()]);
    let sq = |e: &Expr| pow(e, &int(2));
    // t**n * (t + b)**-1
    if let Kind::Pow(base, e) = g.kind() {
        if e == &int(-1) && n == rat(1, 2) {
            if let Some((a, b)) = linear(base, t) {
                if a.is_one() && pos(&b) {
                    let bs = mul(&[b.clone(), s.clone()]);
                    let r = sub(&sqrt(&div(&Expr::pi(), s)), &mul(&[Expr::pi(), sqrt(&b), exp(&bs), erfc(&sqrt(&bs))]));
                    return ok(r, z);
                }
            }
        }
        if e == &int(2) && (n == int(-1) || n == int(-2)) {
            if let Kind::Func(Func::Sin, a) = base.kind() {
                let w = scaled(&a[0], t)?;
                let l = log(&add(&[one.clone(), div(&mul(&[int(4), sq(&w)]), &sq(s))]));
                let plane = mul(&[int(2), absim(&w)]);
                if n == int(-1) {
                    return ok(div(&l, &int(4)), plane);
                }
                let r = sub(&mul(&[w.clone(), atan(&div(&mul(&[int(2), w]), s))]), &div(&mul(&[s.clone(), l]), &int(4)));
                return ok(r, plane);
            }
        }
        return None;
    }
    let Kind::Func(head, args) = g.kind() else { return None };
    match head {
        Func::Exp => {
            if let Some(c) = poly(&args[0], t) {
                if c.len() == 2 {
                    let sa = sub(s, &c[1]);
                    if n.is_one() {
                        return ok(div(&exp(&c[0]), &sq(&sa)), re(&c[1]));
                    }
                    if !is_true(&re_gt(&n, &int(-1))) {
                        return None;
                    }
                    return ok(div(&mul(&[gamma(&np1), exp(&c[0])]), &pow(&sa, &np1)), re(&c[1]));
                }
                if c.len() == 3 && c[1].is_zero() && n.is_one() {
                    // t*exp(-a*t**2)
                    let a = neg(&c[2]);
                    if !re_pos(&a) {
                        return None;
                    }
                    let four_a = mul(&[int(4), a.clone()]);
                    let g = mul(&[
                        sqrt(&div(&Expr::pi(), &four_a)),
                        exp(&div(&sq(s), &four_a)),
                        erfc(&div(s, &sqrt(&four_a))),
                    ]);
                    let r = sub(&inv(&mul(&[int(2), a.clone()])), &mul(&[div(s, &mul(&[int(2), a])), g]));
                    return ok(mul(&[exp(&c[0]), r]), z);
                }
                return None;
            }
            let at = crate::expand::expand(&mul(&[args[0].clone(), t.clone()]));
            if has(&at, t) {
                return None;
            }
            let a = neg(&at);
            let sas = mul(&[int(2), sqrt(&mul(&[a.clone(), s.clone()]))]);
            if n == rat(1, 2) && re_nonneg(&a) {
                let r = mul(&[rat(1, 2), sqrt(&div(&Expr::pi(), &pow(s, &int(3)))), add(&[one.clone(), sas.clone()]), exp(&neg(&sas))]);
                return ok(r, z);
            }
            if n == rat(-1, 2) && re_nonneg(&a) {
                return ok(mul(&[sqrt(&div(&Expr::pi(), s)), exp(&neg(&sas))]), z);
            }
            if n == rat(-3, 2) && re_pos(&a) {
                return ok(mul(&[sqrt(&div(&Expr::pi(), &a)), exp(&neg(&sas))]), z);
            }
            if re_pos(&a) {
                let r = mul(&[int(2), pow(&div(&a, s), &div(&np1, &int(2))), besselk(&np1, &sas)]);
                return ok(r, z);
            }
            None
        }
        Func::Log if &args[0] == t => {
            if n == rat(-1, 2) {
                let r = neg(&mul(&[sqrt(&div(&Expr::pi(), s)), log(&mul(&[int(4), s.clone(), exp(&Expr::constant(Constant::EulerGamma))]))]));
                return ok(r, z);
            }
            if !is_true(&re_gt(&n, &int(-1))) {
                return None;
            }
            let r = mul(&[gamma(&np1), pow(s, &neg(&np1)), sub(&digamma(&np1), &log(s))]);
            ok(r, z)
        }
        Func::Sin if n == int(-1) => {
            let w = scaled(&args[0], t)?;
            ok(atan(&div(&w, s)), absim(&w))
        }
        Func::Sinh | Func::Cosh => {
            let a = scaled(&args[0], t)?;
            let sinh = matches!(head, Func::Sinh);
            if sinh && n == int(-1) {
                return ok(div(&log(&div(&add(&[s.clone(), a.clone()]), &sub(s, &a))), &int(2)), absre(&a));
            }
            let lim = if sinh { int(-2) } else { int(-1) };
            if truth(&crate::relational::gt(&n, &lim)) != Some(true) {
                return None;
            }
            let m = neg(&np1);
            let (p, q) = (pow(&sub(s, &a), &m), pow(&add(&[s.clone(), a.clone()]), &m));
            let inner = if sinh { sub(&p, &q) } else { add(&[p, q]) };
            ok(mul(&[div(&gamma(&np1), &int(2)), inner]), abs(&a))
        }
        Func::Special(sf @ (SpecialFn::BesselJ | SpecialFn::BesselI)) => {
            let (nu, a) = (&args[0], scaled(&args[1], t)?);
            let plus = matches!(sf, SpecialFn::BesselJ);
            let q = if plus { s2pw(s, &a) } else { s2mw(s, &a) };
            let plane = if plus { absim(&a) } else { absre(&a) };
            let sqpi = sqrt(&Expr::pi());
            if crate::expand::expand(&sub(&n, nu)).is_zero() && is_true(&re_gt(nu, &rat(-1, 2))) {
                let h = add(&[nu.clone(), rat(1, 2)]);
                let r = mul(&[pow(&int(2), nu), inv(&sqpi), gamma(&h), pow(&a, nu), pow(&q, &neg(&h))]);
                return ok(r, plane);
            }
            if crate::expand::expand(&sub(&n, &add(&[nu.clone(), one.clone()]))).is_zero() && is_true(&re_gt(nu, &int(-1))) {
                let h = add(&[nu.clone(), rat(3, 2)]);
                let r = mul(&[pow(&int(2), &add(&[nu.clone(), one])), inv(&sqpi), gamma(&h), pow(&a, nu), s.clone(), pow(&q, &neg(&h))]);
                return ok(r, plane);
            }
            None
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------
// program rules
// ---------------------------------------------------------------------

fn prog(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    rule_heaviside(f, t, s, depth)
        .or_else(|| rule_delta(f, t, s))
        .or_else(|| rule_timescale(f, t, s, depth))
        .or_else(|| rule_exp(f, t, s, depth))
        .or_else(|| rule_trig(f, t, s, depth))
        .or_else(|| rule_diff(f, t, s, depth))
        .or_else(|| rule_sdiff(f, t, s, depth))
}

/// Remove the factor at `i` from `fs` and multiply the rest.
fn others(fs: &[Expr], i: usize) -> Expr {
    mul(&fs.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, e)| e.clone()).collect::<Vec<_>>())
}

fn rule_heaviside(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    let fs = factors(f);
    for (i, h) in fs.iter().enumerate() {
        let Kind::Func(Func::Heaviside, a) = h.kind() else { continue };
        let Some((c, d)) = linear(&a[0], t) else { continue };
        let g = others(&fs, i);
        if pos(&c) {
            let sh = neg(&div(&d, &c));
            if pos(&sh) {
                let g2 = crate::subs::subs1(&g, t, &add(&[t.clone(), sh.clone()]));
                let (r, p, cc) = lt(&g2, t, s, depth + 1);
                return Some((mul(&[exp(&neg(&mul(&[sh, s.clone()]))), r]), p, cc));
            }
            if negv(&sh) || sh.is_zero() {
                return Some(lt(&g, t, s, depth + 1));
            }
        } else if negv(&c) {
            let sh = neg(&div(&d, &c));
            if pos(&sh) {
                let w = mul(&[sub(&Expr::one(), &heaviside(&sub(t, &sh))), g]);
                return Some(lt(&crate::expand::expand(&w), t, s, depth + 1));
            }
            if negv(&sh) {
                return ok(Expr::zero(), Expr::zero());
            }
        }
    }
    None
}

fn rule_delta(f: &Expr, t: &Expr, s: &Expr) -> Option<R> {
    let fs = factors(f);
    for (i, h) in fs.iter().enumerate() {
        let Kind::Func(Func::DiracDelta, a) = h.kind() else { continue };
        let z = others(&fs, i);
        if has_head(&z, &|k| matches!(k, Kind::Func(Func::DiracDelta, _))) {
            continue;
        }
        let order = match a.get(1) {
            None => 0u32,
            Some(k) => match k.as_i64() {
                Some(k) if k >= 0 => k as u32,
                _ => continue,
            },
        };
        let h = mul(&[z.clone(), exp(&neg(&mul(&[s.clone(), t.clone()])))]);
        let eval_at = |x: &Expr| -> Option<Expr> {
            let d = if order == 0 { h.clone() } else { crate::diff::diff_n(&h, t, order) };
            let v = crate::subs::subs1(&d, t, x);
            if v.has_nonfinite() {
                return None;
            }
            Some(v)
        };
        if let Some((b, c)) = linear(&a[0], t) {
            let loc = neg(&div(&c, &b));
            if real(&loc) && nonneg(&loc) {
                let v = eval_at(&loc)?;
                let sign = if order % 2 == 1 { int(-1) } else { int(1) };
                let r = div(&mul(&[sign, v]), &mul(&[abs(&b), pow(&b, &int(order as i64))]));
                return ok(r, Expr::neg_oo());
            }
            if negv(&loc) || crate::assumptions::is_real(&loc) == Some(false) {
                return ok(Expr::zero(), Expr::neg_oo());
            }
            continue;
        }
        if order != 0 {
            continue;
        }
        if let Some(c) = poly(&a[0], t) {
            if c.len() < 3 {
                continue;
            }
            let Ok(rts) = crate::poly::roots::roots(&a[0], std::slice::from_ref(t)) else { continue };
            if rts.is_empty() || rts.iter().any(|(_, m)| *m != 1) || rts.iter().map(|_| 1).sum::<usize>() + 1 != c.len() {
                continue;
            }
            let slope = crate::diff::diff(&a[0], t);
            let mut out = Vec::new();
            for (x, _) in &rts {
                let rr = real(x);
                if crate::assumptions::is_real(x) == Some(false) || negv(x) {
                    continue;
                }
                if !(rr && nonneg(x)) {
                    return None;
                }
                let v = eval_at(x)?;
                out.push(div(&v, &abs(&crate::subs::subs1(&slope, t, x))));
            }
            return ok(add(&out), Expr::neg_oo());
        }
    }
    None
}

fn rule_timescale(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    let Kind::Func(head, args) = f.kind() else { return None };
    if args.len() != 1 || matches!(head, Func::Special(_) | Func::Logic(_) | Func::Set(_)) {
        return None;
    }
    let a = scaled(&args[0], t)?;
    if !pos(&a) || a.is_one() {
        return None;
    }
    let g = crate::functions::func(head.clone(), vec![t.clone()]);
    let (r, p, c) = lt(&g, t, &div(s, &a), depth + 1);
    Some((div(&r, &a), mul(&[a, p]), c))
}

fn rule_exp(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    let fs = factors(f);
    for (i, h) in fs.iter().enumerate() {
        let Kind::Func(Func::Exp, a) = h.kind() else { continue };
        let Some((c1, c0)) = linear(&a[0], t) else { continue };
        let z = others(&fs, i);
        let (r, p, cc) = lt(&z, t, &sub(s, &c1), depth + 1);
        return Some((mul(&[exp(&c0), r]), add(&[p, re(&c1)]), cc));
    }
    None
}

fn rule_trig(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    if !has_trig(f) {
        return None;
    }
    let mut trigs = Vec::new();
    let mut other = Vec::new();
    for fac in factors(f) {
        if has_trig(&fac) || has_head(&fac, &|k| matches!(k, Kind::Func(Func::Exp, _))) {
            trigs.push(fac);
        } else {
            other.push(fac);
        }
    }
    let x1 = crate::expand::expand(&rewrite_exp(&mul(&trigs)));
    let mut xm: Vec<(Expr, Expr)> = Vec::new();
    for term in terms(&x1) {
        let (k, a) = if !has(&term, t) { (term.clone(), Expr::zero()) } else { as_exp_term(&term, t)? };
        let a = crate::expand::expand(&a);
        if let Some(e) = xm.iter_mut().find(|(_, b)| b == &a) {
            e.0 = add(&[e.0.clone(), k]);
        } else {
            xm.push((k, a));
        }
    }
    xm.retain(|(k, _)| !k.is_zero());
    let g = mul(&other);
    if !has(&g, t) {
        let (r, p) = pair_poles(xm, s);
        return Some((mul(&[g, r]), p, tru()));
    }
    let (gr, gp, gc) = lt(&g, t, s, depth + 1);
    let mut out = Vec::new();
    let mut planes = Vec::new();
    for (k, a) in xm {
        out.push(mul(&[k, crate::subs::subs1(&gr, s, &sub(s, &a))]));
        planes.push(add(&[gp.clone(), re(&a)]));
    }
    Some((add(&out), max(planes), gc))
}

/// Sum `k_i/(s - a_i)` combining complex-conjugate and symmetric poles
/// into real quadratic denominators; returns the sum and `Max(re(a_i))`.
fn pair_poles(mut xm: Vec<(Expr, Expr)>, s: &Expr) -> (Expr, Expr) {
    let mut results = Vec::new();
    let mut planes = Vec::new();
    let ex = crate::expand::expand;
    while let Some((k0, a0)) = xm.pop() {
        let conj = sub(&re(&a0), &mul(&[Expr::i(), im(&a0)]));
        let nonreal = !ex(&sub(&conj, &a0)).is_zero();
        let j = xm
            .iter()
            .position(|(_, a)| nonreal && ex(&sub(a, &conj)).is_zero())
            .or_else(|| xm.iter().position(|(_, a)| !a0.is_zero() && ex(&add(&[a.clone(), a0.clone()])).is_zero()));
        planes.push(re(&a0));
        match j {
            Some(j) => {
                let (k1, a1) = xm.remove(j);
                planes.push(re(&a1));
                let c1 = ex(&add(&[k0.clone(), k1.clone()]));
                let c0 = ex(&neg(&add(&[mul(&[k0, a1.clone()]), mul(&[k1, a0.clone()])])));
                let num = add(&[mul(&[c1, s.clone()]), c0]);
                let den = ex(&mul(&[sub(s, &a0), sub(s, &a1)]));
                results.push(div(&num, &den));
            }
            None => results.push(div(&k0, &sub(s, &a0))),
        }
    }
    let plane = if planes.len() == 1 { planes.pop().unwrap() } else { max(planes) };
    (add(&results), plane)
}

fn rule_diff(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    let Kind::Derivative(g, vars) = f.kind() else { return None };
    if vars.len() != 1 || &vars[0].0 != t {
        return None;
    }
    let Kind::Func(Func::Undefined(_), gargs) = g.kind() else { return None };
    if gargs.iter().filter(|a| has(a, t)).count() != 1 {
        return None;
    }
    let n = vars[0].1;
    let mut d = Vec::new();
    for k in 0..n {
        let y = if k == 0 {
            crate::subs::subs1(g, t, &Expr::zero())
        } else {
            crate::subs::subs1(&crate::diff::derivative_unevaluated(g, vec![(t.clone(), k)]), t, &Expr::zero())
        };
        d.push(mul(&[pow(s, &int((n - k - 1) as i64)), y]));
    }
    let (r, p, c) = lt(g, t, s, depth + 1);
    Some((sub(&mul(&[pow(s, &int(n as i64)), r]), &add(&d)), p, c))
}

fn has_held_lt(e: &Expr) -> bool {
    has_head(e, &|k| matches!(k, Kind::Func(Func::Special(SpecialFn::LaplaceTransform), _)))
}

fn rule_sdiff(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    if !f.is_mul() {
        return None;
    }
    let mut pfac = Vec::new();
    let mut ofac = Vec::new();
    for fac in factors(f) {
        if poly(&fac, t).is_some() {
            pfac.push(fac);
        } else {
            ofac.push(fac);
        }
    }
    if pfac.is_empty() {
        return None;
    }
    let pc = poly(&mul(&pfac), t)?;
    if pc.len() < 2 {
        return None;
    }
    let (r0, p, c) = lt(&mul(&ofac), t, s, depth + 1);
    let held = has_held_lt(&r0);
    let mut out = vec![mul(&[pc[0].clone(), r0.clone()])];
    let mut cur = r0.clone();
    for (k, ck) in pc.iter().enumerate().skip(1) {
        let dk = if held {
            mul(&[pow(&int(-1), &int(k as i64)), crate::diff::derivative_unevaluated(&r0, vec![(s.clone(), k as u32)])])
        } else {
            cur = neg(&crate::diff::diff(&cur, s));
            cur.clone()
        };
        out.push(mul(&[ck.clone(), dk]));
    }
    Some((add(&out), p, c))
}

fn expand_rule(f: &Expr, t: &Expr, s: &Expr, depth: u32) -> Option<R> {
    let r = crate::expand::expand(f);
    if r.is_add() || &r != f {
        return Some(lt(&r, t, s, depth + 1));
    }
    let r = crate::expand::expand(&crate::simp::expand_trig(f));
    if r.is_add() {
        return Some(lt(&r, t, s, depth + 1));
    }
    None
}
