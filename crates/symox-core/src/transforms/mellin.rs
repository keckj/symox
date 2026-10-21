//! Mellin transform `F(s) = int_0^oo x**(s-1) f(x) dx` with its fundamental
//! strip `a < Re(s) < b`, and the inverse transform.
//!
//! The forward engine decomposes `f = C * x**n * g(b*x**p)` where `g` is
//! one of the classical Mellin pairs (exponential, Heaviside cut offs,
//! binomials `(1 + u)**(-rho)`, `log(1 + u)`, `sin`, `cos`, Bessel `J` and
//! `K`, `erfc`) and applies the shift and scaling theorems. The inverse
//! engine recognises the same families on the `s` side (gamma factors,
//! reflection products `Gamma(z)*Gamma(rho - z)`, `pi/sin(pi*s)`, simple
//! poles `1/(s + n)` and scaling factors `c**(-s)`).

use super::util::*;
use crate::expr::{Expr, Func, Kind};
use crate::functions::{exp, gamma, heaviside, sin};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};
use crate::special::SpecialFn;

type R = (Expr, Expr, Expr, Expr);

/// Mellin transform of `f`: `Some((F, (a, b), cond))`, `None` when no rule
/// applies.
pub fn mellin_transform(f: &Expr, x: &Expr, s: &Expr, simplify: bool) -> Option<(Expr, (Expr, Expr), Expr)> {
    let (r, a, b, c) = mt(f, x, s, 0)?;
    if r.has_nonfinite() {
        return None;
    }
    let r = gamma_reflect(&r);
    let r = if simplify { super::fourier::finish(&r) } else { r };
    Some((r, (a, b), c))
}

fn oo() -> Expr {
    Expr::oo()
}

fn strip_max(a: &Expr, b: &Expr) -> Expr {
    crate::functions::max(vec![a.clone(), b.clone()])
}
fn strip_min(a: &Expr, b: &Expr) -> Expr {
    crate::functions::min(vec![a.clone(), b.clone()])
}

fn mt(f: &Expr, x: &Expr, s: &Expr, depth: u32) -> Option<R> {
    if depth > MAX_DEPTH || !has(f, x) {
        return None;
    }
    if f.is_add() {
        let mut vals = Vec::new();
        let mut a = Expr::neg_oo();
        let mut b = oo();
        let mut conds = Vec::new();
        for t in terms(f) {
            let (v, a1, b1, c) = mt(&t, x, s, depth + 1)?;
            vals.push(v);
            a = strip_max(&a, &a1);
            b = strip_min(&b, &b1);
            conds.push(c);
        }
        if nonneg(&sub(&a, &b)) {
            return None;
        }
        return Some((add(&vals), a, b, crate::logic::and(conds)));
    }
    let (c, g) = split_mul(f, x);
    if !c.is_one() {
        let (v, a, b, cd) = mt(&g, x, s, depth + 1)?;
        return Some((mul(&[c, v]), a, b, cd));
    }
    let (n, g) = split_power(&g, x);
    if !n.is_zero() {
        if !has(&g, x) {
            return None;
        }
        // M[x**n g](s) = G(s + n)
        let sn = add(&[s.clone(), n.clone()]);
        let (v, a, b, cd) = base(&g, x, &sn)?;
        // Re(s + n) in (a, b): the strip moves by Re(n)
        let rn = crate::functions::re(&n);
        let shift = |e: &Expr| if *e == Expr::oo() || *e == Expr::neg_oo() { e.clone() } else { sub(e, &rn) };
        return Some((v, shift(&a), shift(&b), cd));
    }
    if let Some(r) = base(&g, x, s) {
        return Some(r);
    }
    // products with a polynomial factor, e.g. (x + 1)*exp(-x)
    if g.is_mul() && factors(&g).iter().any(|a| a.is_add() && poly(a, x).is_some()) {
        let e = crate::expand::expand(&g);
        if e.is_add() {
            return mt(&e, x, s, depth + 1);
        }
    }
    None
}

/// `u == b * x**p` with `p != 0` and `b` free of `x`: returns `(b, p)`.
fn monomial(u: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    let (b, r) = split_mul(u, x);
    let p = power_of(&r, x)?;
    if p.is_zero() || has(&p, x) {
        return None;
    }
    Some((b, p))
}

/// `u == alpha + b*x**p`: returns `(alpha, b, p)`.
fn affine_monomial(u: &Expr, x: &Expr) -> Option<(Expr, Expr, Expr)> {
    let mut alpha = Vec::new();
    let mut m = None;
    for t in terms(u) {
        if has(&t, x) {
            if m.is_some() {
                return None;
            }
            m = Some(monomial(&t, x)?);
        } else {
            alpha.push(t);
        }
    }
    let (b, p) = m?;
    Some((add(&alpha), b, p))
}

/// Transform of `g0(b*x**p)` from the table entry `(G, a0, b0, cond)` of
/// `g0`: `(1/|p|) * b**(-s/p) * G(s/p)`.
fn scale(table: &dyn Fn(&Expr) -> Option<R>, b: &Expr, p: &Expr, s: &Expr) -> Option<R> {
    let sp = div(s, p);
    let (g, a0, b0, cd) = table(&sp)?;
    let (a, bb) = if pos(p) {
        (mul(&[p.clone(), a0]), mul(&[p.clone(), b0]))
    } else if negv(p) {
        (mul(&[p.clone(), b0]), mul(&[p.clone(), a0]))
    } else {
        return None;
    };
    let v = mul(&[div(&Expr::one(), &abs_num(p)), pow(b, &neg(&sp)), g]);
    Some((v, a, bb, cd))
}

fn abs_num(p: &Expr) -> Expr {
    if negv(p) { neg(p) } else { p.clone() }
}

/// The table of Mellin pairs for a single factor `g` (with general argument).
fn base(g: &Expr, x: &Expr, s: &Expr) -> Option<R> {
    let pi = Expr::pi();
    match g.kind() {
        Kind::Func(Func::Exp, a) => {
            // exp(-b*x**p)
            let (b, p) = monomial(&neg(&a[0]), x)?;
            let rb = crate::functions::re(&b);
            if negv(&rb) || rb.is_zero() {
                return None;
            }
            let cond = if re_pos(&b) { tru() } else { crate::relational::gt(&rb, &Expr::zero()) };
            let (v, a, bb, _) = scale(&|z| Some((gamma(z), Expr::zero(), oo(), tru())), &b, &p, s)?;
            Some((v, a, bb, cond))
        }
        Kind::Func(Func::Heaviside, a) => {
            let (al, b, p) = affine_monomial(&a[0], x)?;
            if pos(&al) && negv(&b) {
                // H(1 - u), u = -b/alpha x**p
                let u = neg(&div(&b, &al));
                scale(&|z| Some((div(&Expr::one(), z), Expr::zero(), oo(), tru())), &u, &p, s)
            } else if negv(&al) && pos(&b) {
                let u = neg(&div(&b, &al));
                scale(&|z| Some((neg(&div(&Expr::one(), z)), Expr::neg_oo(), Expr::zero(), tru())), &u, &p, s)
            } else {
                None
            }
        }
        Kind::Pow(bs, e) if !has(e, x) => {
            let (al, b, p) = affine_monomial(bs, x)?;
            if !pos(&al) || !pos(&b) {
                return None;
            }
            let rho = neg(e);
            let rr = crate::functions::re(&rho);
            if !pos(&rr) {
                return None;
            }
            let u = div(&b, &al);
            let (v, a, bb, cd) = scale(
                &|z| Some((div(&mul(&[gamma(z), gamma(&sub(&rho, z))]), &gamma(&rho)), Expr::zero(), rr.clone(), tru())),
                &u,
                &p,
                s,
            )?;
            Some((mul(&[pow(&al, e), v]), a, bb, cd))
        }
        Kind::Func(Func::Log, a) => {
            let (al, b, p) = affine_monomial(&a[0], x)?;
            if !al.is_one() || !pos(&b) {
                return None;
            }
            scale(
                &|z| Some((div(&pi, &mul(&[z.clone(), sin(&mul(&[pi.clone(), z.clone()]))])), int(-1), Expr::zero(), tru())),
                &b,
                &p,
                s,
            )
        }
        Kind::Func(Func::Sin, a) => {
            let (b, p) = monomial(&a[0], x)?;
            if !pos(&b) {
                return None;
            }
            scale(
                &|z| {
                    Some((
                        div(
                            &mul(&[pow(&int(2), &sub(z, &Expr::one())), sqrt(&Expr::pi()), gamma(&div(&add(&[z.clone(), Expr::one()]), &int(2)))]),
                            &gamma(&sub(&Expr::one(), &div(z, &int(2)))),
                        ),
                        int(-1),
                        Expr::one(),
                        tru(),
                    ))
                },
                &b,
                &p,
                s,
            )
        }
        Kind::Func(Func::Cos, a) => {
            let (b, p) = monomial(&a[0], x)?;
            if !pos(&b) {
                return None;
            }
            scale(
                &|z| {
                    Some((
                        div(
                            &mul(&[pow(&int(2), &sub(z, &Expr::one())), sqrt(&Expr::pi()), gamma(&div(z, &int(2)))]),
                            &gamma(&div(&sub(&Expr::one(), z), &int(2))),
                        ),
                        Expr::zero(),
                        Expr::one(),
                        tru(),
                    ))
                },
                &b,
                &p,
                s,
            )
        }
        Kind::Func(Func::Erfc, a) => {
            let (b, p) = monomial(&a[0], x)?;
            if !pos(&b) {
                return None;
            }
            scale(
                &|z| Some((div(&gamma(&div(&add(&[z.clone(), Expr::one()]), &int(2))), &mul(&[z.clone(), sqrt(&Expr::pi())])), Expr::zero(), oo(), tru())),
                &b,
                &p,
                s,
            )
        }
        Kind::Func(Func::Special(SpecialFn::BesselJ), a) if !has(&a[0], x) => {
            let nu = a[0].clone();
            let (b, p) = monomial(&a[1], x)?;
            if !pos(&b) {
                return None;
            }
            scale(
                &|z| {
                    Some((
                        div(
                            &mul(&[pow(&int(2), &sub(z, &Expr::one())), gamma(&div(&add(&[nu.clone(), z.clone()]), &int(2)))]),
                            &gamma(&add(&[div(&sub(&nu, z), &int(2)), Expr::one()])),
                        ),
                        neg(&nu),
                        rat(3, 2),
                        tru(),
                    ))
                },
                &b,
                &p,
                s,
            )
        }
        Kind::Func(Func::Special(SpecialFn::BesselK), a) if !has(&a[0], x) => {
            let nu = a[0].clone();
            let (b, p) = monomial(&a[1], x)?;
            if !pos(&b) {
                return None;
            }
            let anu = if nonneg(&nu) { nu.clone() } else if nonpos(&nu) { neg(&nu) } else { crate::functions::abs(&nu) };
            scale(
                &|z| {
                    Some((
                        mul(&[
                            pow(&int(2), &sub(z, &int(2))),
                            gamma(&div(&sub(z, &nu), &int(2))),
                            gamma(&div(&add(&[z.clone(), nu.clone()]), &int(2))),
                        ]),
                        anu.clone(),
                        oo(),
                        tru(),
                    ))
                },
                &b,
                &p,
                s,
            )
        }
        _ => None,
    }
}

/// `Gamma(z)*Gamma(m - z)` with an integer `m >= 1` becomes
/// `pi/sin(pi*z) * prod_{j=1}^{m-1} (j - z)` (also in denominators).
/// Legendre duplication: `Gamma(u)/Gamma(u/2 + 1/2)` or `Gamma(u)/Gamma(u/2)`
/// in a product is rewritten with `Gamma(u) = 2**(u-1) Gamma(u/2)
/// Gamma(u/2 + 1/2)/sqrt(pi)`, so that the denominator cancels.
pub(crate) fn gamma_duplicate(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| gamma_duplicate(c));
    if !e.is_mul() {
        return e;
    }
    let fs = e.args();
    let num: Vec<Expr> = fs
        .iter()
        .filter_map(|f| match f.kind() {
            Kind::Func(Func::Gamma, a) => Some(a[0].clone()),
            _ => None,
        })
        .collect();
    let den: Vec<Expr> = fs
        .iter()
        .filter_map(|f| match f.kind() {
            Kind::Pow(b, n) if n.as_i64() == Some(-1) => match b.kind() {
                Kind::Func(Func::Gamma, a) => Some(a[0].clone()),
                _ => None,
            },
            _ => None,
        })
        .collect();
    for u in &num {
        for v in &den {
            let d = crate::expand::expand(&sub(&mul(&[int(2), v.clone()]), u));
            if d.is_zero() || d.is_one() {
                let h = div(u, &int(2));
                let rep = div(
                    &mul(&[pow(&int(2), &sub(u, &Expr::one())), gamma(&h), gamma(&add(&[h.clone(), rat(1, 2)]))]),
                    &sqrt(&Expr::pi()),
                );
                let g = gamma(u);
                let mut out = Vec::new();
                let mut done = false;
                for f in &fs {
                    if !done && *f == g {
                        out.push(rep.clone());
                        done = true;
                    } else {
                        out.push(f.clone());
                    }
                }
                return mul(&out);
            }
        }
    }
    e
}

pub(crate) fn gamma_reflect(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| gamma_reflect(c));
    if !e.is_mul() {
        return e;
    }
    let mut fs: Vec<Option<Expr>> = e.args().into_iter().map(Some).collect();
    let arg_exp = |f: &Expr| -> Option<(Expr, i64)> {
        match f.kind() {
            Kind::Func(Func::Gamma, a) => Some((a[0].clone(), 1)),
            Kind::Pow(b, n) => match (b.kind(), n.as_i64()) {
                (Kind::Func(Func::Gamma, a), Some(m @ (1 | -1))) => Some((a[0].clone(), m)),
                _ => None,
            },
            _ => None,
        }
    };
    let mut out = Vec::new();
    let mut changed = false;
    for j in 0..fs.len() {
        let Some(fj) = fs[j].clone() else { continue };
        let Some((zj, ej)) = arg_exp(&fj) else { continue };
        for l in j + 1..fs.len() {
            let Some(fl) = fs[l].clone() else { continue };
            let Some((zl, el)) = arg_exp(&fl) else { continue };
            if el != ej {
                continue;
            }
            let m = crate::expand::expand(&add(&[zj.clone(), zl.clone()]));
            let Some(m) = m.as_i64() else { continue };
            if m < 1 || m > 8 {
                continue;
            }
            // choose z with a positive leading coefficient when possible
            let negative = |z: &Expr| terms(&crate::expand::expand(z)).iter().any(|t| !matches!(t.kind(), Kind::Num(_)) && crate::ops::could_extract_minus_sign(t));
            let z = if negative(&zj) { zl.clone() } else { zj.clone() };
            let mut parts = vec![div(&Expr::pi(), &sin(&mul(&[Expr::pi(), z.clone()])))];
            for i in 1..m {
                parts.push(sub(&int(i), &z));
            }
            let r = mul(&parts);
            out.push(if ej == 1 { r } else { div(&Expr::one(), &r) });
            fs[j] = None;
            fs[l] = None;
            changed = true;
            break;
        }
    }
    if !changed {
        return e;
    }
    out.extend(fs.into_iter().flatten());
    mul(&out)
}

// ---------------------------------------------------------------------------
// inverse transform
// ---------------------------------------------------------------------------

/// Inverse Mellin transform of `F(s)` on the strip `(a, b)`:
/// `Some((f, cond))`, `None` when no rule applies.
pub fn inverse_mellin_transform(f: &Expr, s: &Expr, x: &Expr, strip: (&Expr, &Expr), simplify: bool) -> Option<(Expr, Expr)> {
    let r = imt(f, s, x, strip.0, strip.1, 0)?;
    if r.has_nonfinite() {
        return None;
    }
    let r = if simplify { super::fourier::finish(&r) } else { r };
    Some((r, tru()))
}

fn imt(f: &Expr, s: &Expr, x: &Expr, a: &Expr, b: &Expr, depth: u32) -> Option<Expr> {
    if depth > MAX_DEPTH || !has(f, s) {
        return None;
    }
    if f.is_add() {
        let mut vals = Vec::new();
        for t in terms(f) {
            vals.push(imt(&t, s, x, a, b, depth + 1)?);
        }
        return Some(add(&vals));
    }
    let (c, g) = split_mul(f, s);
    if !c.is_one() {
        return Some(mul(&[c, imt(&g, s, x, a, b, depth + 1)?]));
    }
    // scaling factors c**(gamma*s): F(s) * c**(-s) <-> f(c*x)
    let mut sc = Vec::new();
    let mut rest = Vec::new();
    for fct in factors(&g) {
        match fct.kind() {
            Kind::Pow(bs, e) if !has(bs, s) => {
                if let Some((ga, de)) = linear(e, s) {
                    sc.push(pow(bs, &neg(&ga)));
                    rest.push(pow(bs, &de));
                    continue;
                }
                rest.push(fct.clone())
            }
            Kind::Func(Func::Exp, ar) => {
                if let Some((ga, de)) = linear(&ar[0], s) {
                    sc.push(exp(&neg(&ga)));
                    rest.push(exp(&de));
                    continue;
                }
                rest.push(fct.clone())
            }
            _ => rest.push(fct.clone()),
        }
    }
    if !sc.is_empty() {
        let cs = mul(&sc);
        let h = imt(&mul(&rest), s, x, a, b, depth + 1)?;
        return Some(crate::subs::subs1(&h, x, &mul(&[cs, x.clone()])));
    }
    let g = absorb(&reflect_sin(&g, s), s);
    let fs = factors(&g);
    let mut gam = Vec::new();
    let mut other = Vec::new();
    for f in &fs {
        match f.kind() {
            Kind::Func(Func::Gamma, ar) => gam.push(ar[0].clone()),
            _ => other.push(f.clone()),
        }
    }
    let k = mul(&other);
    if has(&k, s) {
        // simple pole 1/(alpha*s + beta)
        if gam.is_empty() {
            if let Kind::Pow(bs, e) = k.kind() {
                if e.as_i64() == Some(-1) {
                    let (al, be) = linear(bs, s)?;
                    let n = div(&be, &al);
                    let pole = neg(&n);
                    let xn = pow(x, &n);
                    if nonneg(&sub(a, &pole)) {
                        return Some(div(&mul(&[xn, heaviside(&sub(&Expr::one(), x))]), &al));
                    }
                    if nonpos(&sub(b, &pole)) {
                        return Some(neg(&div(&mul(&[xn, heaviside(&sub(x, &Expr::one()))]), &al)));
                    }
                }
            }
        }
        return None;
    }
    match gam.len() {
        1 => {
            // Gamma(alpha*s + beta) <-> x**(beta/alpha) exp(-x**(1/alpha)) / |alpha|
            let (al, be) = linear(&gam[0], s)?;
            let n = div(&be, &al);
            let pole = neg(&n);
            if pos(&al) {
                if negv(&sub(a, &pole)) {
                    return None;
                }
            } else if negv(&al) {
                if pos(&sub(b, &pole)) {
                    return None;
                }
            } else {
                return None;
            }
            let p = div(&Expr::one(), &al);
            Some(mul(&[k, div(&Expr::one(), &abs_num(&al)), pow(x, &n), exp(&neg(&pow(x, &p)))]))
        }
        2 => {
            // Gamma(alpha*s + beta) * Gamma(gamma - alpha*s)
            let (a1, b1) = linear(&gam[0], s)?;
            let (a2, b2) = linear(&gam[1], s)?;
            if !crate::expand::expand(&add(&[a1.clone(), a2.clone()])).is_zero() {
                return None;
            }
            let (al, be, ga) = if pos(&a1) { (a1, b1, b2) } else if pos(&a2) { (a2, b2, b1) } else { return None };
            let n = div(&be, &al);
            let rho = add(&[be, ga]);
            if !pos(&crate::functions::re(&rho)) {
                return None;
            }
            let p = div(&Expr::one(), &al);
            let lo = neg(&n);
            let hi = sub(&mul(&[p.clone(), rho.clone()]), &n);
            if negv(&sub(a, &lo)) || pos(&sub(b, &hi)) {
                return None;
            }
            Some(mul(&[k, p.clone(), gamma(&rho), pow(x, &n), pow(&add(&[Expr::one(), pow(x, &p)]), &neg(&rho))]))
        }
        _ => None,
    }
}

/// `1/sin(pi*z)` becomes `Gamma(z)*Gamma(1 - z)/pi`.
fn reflect_sin(g: &Expr, s: &Expr) -> Expr {
    let mut out = Vec::new();
    for f in factors(g) {
        if let Kind::Pow(b, e) = f.kind() {
            if let (Kind::Func(Func::Sin, a), Some(-1)) = (b.kind(), e.as_i64()) {
                let z = crate::expand::expand(&div(&a[0], &Expr::pi()));
                if linear(&z, s).is_some() {
                    out.push(gamma(&z));
                    out.push(gamma(&sub(&Expr::one(), &z)));
                    out.push(div(&Expr::one(), &Expr::pi()));
                    continue;
                }
            }
        }
        out.push(f);
    }
    mul(&out)
}

/// Absorb linear factors into gamma functions: `z * Gamma(z) = Gamma(z + 1)`.
fn absorb(g: &Expr, s: &Expr) -> Expr {
    let mut fs = factors(g);
    loop {
        let mut done = true;
        'outer: for j in 0..fs.len() {
            if linear(&fs[j], s).is_none() {
                continue;
            }
            for l in 0..fs.len() {
                if let Kind::Func(Func::Gamma, a) = fs[l].kind() {
                    let z = a[0].clone();
                    let d = crate::expand::expand(&sub(&fs[j], &z));
                    if d.is_zero() {
                        fs[l] = gamma(&add(&[z, Expr::one()]));
                        fs.remove(j);
                        done = false;
                        break 'outer;
                    }
                    let d = crate::expand::expand(&add(&[fs[j].clone(), z.clone()]));
                    if d.is_zero() {
                        // -z * Gamma(z) = -Gamma(z + 1)
                        fs[l] = neg(&gamma(&add(&[z, Expr::one()])));
                        fs.remove(j);
                        done = false;
                        break 'outer;
                    }
                }
            }
        }
        if done {
            break;
        }
    }
    mul(&fs)
}

