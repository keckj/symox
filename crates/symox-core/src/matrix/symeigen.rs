//! Eigenvalues and eigenvectors of matrices over rational function fields
//! `K = QQ(x, y, ...)`, reproducing SymPy's `DomainMatrix` based
//! `eigenvects` and its `roots()` for symbolic quadratics.
//!
//! Only "easy" spectra are supported: characteristic polynomials that split
//! into linear factors found among simple candidates (diagonal entries,
//! zero, square roots of discriminants) and at most one irreducible
//! quadratic factor.

use std::cmp::Ordering;

use super::ff::{berkowitz, Ring};
use super::int::{Int, Rat};
use super::mpoly::{MPoly, Mono};
use super::ratfunc::{as_numer_denom, cancel, cancel_polys, PolyCtx};
use super::{MResult, Matrix, MatrixError};
use crate::expr::{Expr, Kind};
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};

// ----------------------------------------------------------------------
// the field K
// ----------------------------------------------------------------------

/// A reduced rational function `n / d` (positive leading coefficient of
/// `d`).
#[derive(Clone, Debug, PartialEq)]
pub struct KF {
    pub n: MPoly,
    pub d: MPoly,
}

impl KF {
    pub fn new(n: MPoly, d: MPoly) -> KF {
        let (p, q) = cancel_polys(&n, &d);
        KF { n: p, d: q }
    }

    pub fn from_poly(p: MPoly) -> KF {
        let nv = p.nvars;
        KF { n: p, d: MPoly::one(nv) }
    }

    pub fn zero(nv: usize) -> KF {
        KF { n: MPoly::zero(nv), d: MPoly::one(nv) }
    }

    pub fn one(nv: usize) -> KF {
        KF { n: MPoly::one(nv), d: MPoly::one(nv) }
    }

    pub fn from_expr(ctx: &PolyCtx, e: &Expr) -> Option<KF> {
        let (a, b) = as_numer_denom(e);
        let qa = ctx.to_qpoly(&a)?;
        let qb = ctx.to_qpoly(&b)?;
        if qb.p.is_zero() {
            return None;
        }
        Some(KF::new(qa.p.mul_int(&qb.den), qb.p.mul_int(&qa.den)))
    }

    pub fn to_expr(&self, ctx: &PolyCtx) -> Expr {
        let ne = ctx.to_expr(&self.n);
        if self.d.is_one() {
            ne
        } else {
            div(&ne, &ctx.to_expr(&self.d))
        }
    }

    pub fn is_zero(&self) -> bool {
        self.n.is_zero()
    }

    pub fn add(&self, o: &KF) -> KF {
        if self.d == o.d {
            return KF::new(self.n.add(&o.n), self.d.clone());
        }
        KF::new(self.n.mul(&o.d).add(&o.n.mul(&self.d)), self.d.mul(&o.d))
    }

    pub fn sub(&self, o: &KF) -> KF {
        self.add(&o.neg())
    }

    pub fn neg(&self) -> KF {
        KF { n: self.n.neg(), d: self.d.clone() }
    }

    pub fn mul(&self, o: &KF) -> KF {
        KF::new(self.n.mul(&o.n), self.d.mul(&o.d))
    }

    pub fn inv(&self) -> KF {
        KF::new(self.d.clone(), self.n.clone())
    }

    pub fn div(&self, o: &KF) -> KF {
        KF::new(self.n.mul(&o.d), self.d.mul(&o.n))
    }
}

impl Ring for KF {
    fn zero_like(&self) -> KF {
        KF::zero(self.n.nvars)
    }
    fn one_like(&self) -> KF {
        KF::one(self.n.nvars)
    }
    fn is_zero(&self) -> bool {
        KF::is_zero(self)
    }
    fn is_one(&self) -> bool {
        self.n.is_one() && self.d.is_one()
    }
    fn add(&self, o: &KF) -> KF {
        KF::add(self, o)
    }
    fn sub(&self, o: &KF) -> KF {
        KF::sub(self, o)
    }
    fn mul(&self, o: &KF) -> KF {
        KF::mul(self, o)
    }
    fn neg(&self) -> KF {
        KF::neg(self)
    }
    fn exquo(&self, o: &KF) -> KF {
        KF::div(self, o)
    }
}

/// Square root of a polynomial with integer coefficients, if it is a
/// perfect square.
pub fn mpoly_sqrt(p: &MPoly) -> Option<MPoly> {
    let n = p.nvars;
    if p.is_zero() {
        return Some(p.clone());
    }
    let (lm, lc) = p.terms.first()?;
    if lm.iter().any(|e| e % 2 == 1) || lc.is_negative() {
        return None;
    }
    let c = lc.exact_sqrt()?;
    let m: Mono = lm.iter().map(|e| e / 2).collect();
    let mut s = MPoly::monomial(n, m, c);
    for _ in 0..p.terms.len() + 2 {
        let r = p.sub(&s.mul(&s));
        if r.is_zero() {
            return Some(s);
        }
        let (rm, rc) = &r.terms[0];
        let (sm, sc) = &s.terms[0];
        let mut nm = Mono::with_capacity(n);
        for (a, b) in rm.iter().zip(sm.iter()) {
            if a < b {
                return None;
            }
            nm.push(a - b);
        }
        let two_sc = sc.mul(&Int::from(2));
        if !rc.is_divisible_by(&two_sc) {
            return None;
        }
        s = s.add(&MPoly::monomial(n, nm, rc.divexact(&two_sc)));
    }
    None
}

fn sqrt_k(x: &KF) -> Option<KF> {
    // sqrt(n/d) = sqrt(n*d)/d
    let s = mpoly_sqrt(&x.n.mul(&x.d))?;
    Some(KF::new(s, x.d.clone()))
}

// ----------------------------------------------------------------------
// light factorisation (SymPy's factor() for the cases we need)
// ----------------------------------------------------------------------

/// `coeff * prod(f_i ** k_i)`; `factors` are primitive polynomials with
/// positive leading coefficient (as expressions).
#[derive(Clone, Debug)]
pub struct Factored {
    pub coeff: Rat,
    pub factors: Vec<(Expr, i64)>,
}

impl Factored {
    /// The (canonical) expression.
    pub fn to_expr(&self) -> Expr {
        let mut v = vec![Expr::num(self.coeff.to_number())];
        for (f, k) in &self.factors {
            v.push(pow(f, &Expr::int(*k)));
        }
        mul(&v)
    }

    /// SymPy's `Mul.make_args(factor(...))` (the coefficient multiplying a
    /// single sum stays separate, like SymPy's unevaluated `_keep_coeff`).
    fn make_args(&self) -> Vec<Expr> {
        let e = self.to_expr();
        if self.factors.len() == 1 && self.factors[0].1 == 1 && !self.coeff.is_one() && !self.coeff.neg().is_one() {
            return vec![Expr::num(self.coeff.to_number()), self.factors[0].0.clone()];
        }
        match e.kind() {
            Kind::Mul(..) => e.args(),
            _ => vec![e],
        }
    }
}

/// Univariate rational root splitting of a primitive integer polynomial in
/// a single variable `v`: returns linear factors and the remaining part.
fn split_univariate(p: &MPoly, v: usize) -> (Vec<(MPoly, i64)>, MPoly) {
    let nv = p.nvars;
    let d = p.degree(v) as usize;
    let mut coeffs: Vec<Rat> = vec![Rat::zero(); d + 1];
    for (m, c) in &p.terms {
        coeffs[d - m[v] as usize] = Rat::from_int(c.clone());
    }
    let mut out = Vec::new();
    let roots = match super::eigen::rational_roots_pub(&coeffs) {
        Ok(r) => r,
        Err(_) => return (vec![], p.clone()),
    };
    let mut rest = p.clone();
    for r in roots {
        // factor q*v - p with multiplicity
        let mut mono = Mono::from_elem(0, nv);
        mono[v] = 1;
        let lin = MPoly::from_terms(nv, vec![(mono, r.den.clone()), (Mono::from_elem(0, nv), r.num.neg())]);
        let mut k = 0;
        while let Some(q) = rest.divexact(&lin) {
            rest = q;
            k += 1;
        }
        if k > 0 {
            out.push((lin, k));
        }
    }
    (out, rest)
}

/// Factor an integer polynomial into (sign*content, factors).
fn factor_poly(p: &MPoly, ctx: &PolyCtx) -> (Int, Vec<(Expr, i64)>) {
    let nv = p.nvars;
    if p.is_zero() {
        return (Int::ZERO, vec![]);
    }
    let mut c = p.content();
    if p.lc().is_negative() {
        c = c.neg();
    }
    let mut prim = p.divexact_int(&c);
    let mut out: Vec<(Expr, i64)> = Vec::new();
    // monomial content
    if !prim.is_constant() {
        let mut m = prim.terms[0].0.clone();
        for (t, _) in &prim.terms[1..] {
            for (a, b) in m.iter_mut().zip(t.iter()) {
                *a = (*a).min(*b);
            }
        }
        for (v, &e) in m.iter().enumerate() {
            if e > 0 {
                out.push((ctx.gens[v].clone(), e as i64));
            }
        }
        if m.iter().any(|&e| e > 0) {
            let terms = prim
                .terms
                .iter()
                .map(|(t, k)| (t.iter().zip(m.iter()).map(|(a, b)| a - b).collect(), k.clone()))
                .collect();
            prim = MPoly { nvars: nv, terms };
        }
    }
    if prim.is_constant() {
        return (c, out);
    }
    // perfect powers
    let mut k = 1i64;
    while let Some(s) = mpoly_sqrt(&prim) {
        if s.is_constant() {
            break;
        }
        prim = s;
        k *= 2;
    }
    let used: Vec<usize> = (0..nv).filter(|&v| prim.uses_var(v)).collect();
    if used.len() == 1 {
        let (lins, rest) = split_univariate(&prim, used[0]);
        for (l, m) in lins {
            out.push((ctx.to_expr(&l), m * k));
        }
        if !rest.is_constant() {
            out.push((ctx.to_expr(&rest.normalize_sign()), k));
        } else if rest.lc().is_negative() {
            c = c.neg();
        }
    } else {
        out.push((ctx.to_expr(&prim), k));
    }
    (c, out)
}

/// SymPy's `factor()` restricted to content, monomial factors, perfect
/// squares and univariate rational roots.
pub fn factor_light(e: &Expr) -> Factored {
    if let Some(n) = e.as_number() {
        if let Some(r) = Rat::from_number(n) {
            return Factored { coeff: r, factors: vec![] };
        }
    }
    let (n, d) = as_numer_denom(e);
    let Some((ctx, ps)) = PolyCtx::from_exprs(&[n, d], false) else {
        return Factored { coeff: Rat::one(), factors: vec![(e.clone(), 1)] };
    };
    let (cn, mut fn_) = factor_poly(&ps[0].p, &ctx);
    let (cd, fd) = factor_poly(&ps[1].p, &ctx);
    // value = (cn * fn / qn.den) / (cd * fd / qd.den)
    let coeff = Rat::new(cn.mul(&ps[1].den), cd.mul(&ps[0].den));
    for (f, k) in fd {
        fn_.push((f, -k));
    }
    Factored { coeff, factors: fn_ }
}

/// SymPy's `_sqrt` in `roots_quadratic`, applied to a factored radicand.
fn sqrt_factored(f: &Factored) -> Expr {
    let args = f.make_args();
    let mut co = Vec::new();
    let mut other = Vec::new();
    for a in &args {
        match a.kind() {
            Kind::Pow(b, x) if x.as_i64().is_some_and(|k| k % 2 == 0) => {
                co.push(pow(b, &Expr::int(x.as_i64().unwrap() / 2)));
            }
            _ => other.push(a.clone()),
        }
    }
    if co.is_empty() {
        return sqrt_keep_coeff(f);
    }
    mul2(&mul(&co), &sqrt(&mul(&other)))
}

/// `sqrt` of SymPy's (possibly unevaluated) `coeff * sum`: the positive
/// numeric part is taken out, a negative sign goes into the sum.
fn sqrt_keep_coeff(f: &Factored) -> Expr {
    if f.factors.len() == 1 && f.factors[0].1 == 1 && !f.coeff.is_one() {
        let c = &f.coeff;
        let ac = if c.is_negative() { c.neg() } else { c.clone() };
        let inner = if c.is_negative() { neg(&f.factors[0].0) } else { f.factors[0].0.clone() };
        return mul2(&sqrt(&Expr::num(ac.to_number())), &sqrt(&inner));
    }
    sqrt(&f.to_expr())
}

/// `factor_terms` on the radicands of square roots: take the rational
/// content out (`sqrt(4*y + 4) -> 2*sqrt(y + 1)`).
fn factor_terms_radicals(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Pow(b, x) if *x == Expr::rational(1, 2) && b.is_add() => {
            let Kind::Add(c, terms) = b.kind() else { unreachable!() };
            let mut g: Option<crate::number::Number> = if c.is_zero() { None } else { Some(c.abs()) };
            for (_, k) in terms {
                if !k.is_rational() {
                    return e.clone();
                }
                g = Some(match g {
                    None => k.abs(),
                    Some(x) => x.gcd(k),
                });
            }
            let g = g.unwrap_or(crate::number::Number::ONE);
            if g.is_one() || g.is_zero() {
                return e.clone();
            }
            let inner = crate::expand::expand(&mul2(&Expr::num(g.inv().unwrap()), b));
            mul2(&sqrt(&Expr::num(g)), &sqrt(&inner))
        }
        _ if e.is_atom() => e.clone(),
        _ => crate::traverse::map_children(e, &mut |c| factor_terms_radicals(c)),
    }
}

// ----------------------------------------------------------------------
// roots of polynomials with symbolic coefficients
// ----------------------------------------------------------------------

/// SymPy's `roots_quadratic` for a composite (polynomial) domain.
fn roots_quadratic_sym(a: &Expr, b: &Expr, c: &Expr) -> Vec<Expr> {
    let ia = pow(a, &Expr::minus_one());
    if c.is_zero() {
        let r1 = factor_light(&neg(&mul2(b, &ia))).to_expr();
        return vec![Expr::zero(), r1];
    }
    if b.is_zero() {
        let r = factor_light(&neg(&mul2(c, &ia)));
        let rr = sqrt_factored(&r);
        return vec![neg(&rr), rr];
    }
    let d = factor_light(&sub(&pow(b, &Expr::int(2)), &mul(&[Expr::int(4), a.clone(), c.clone()])));
    let aa = mul2(&Expr::int(2), a);
    let iaa = pow(&aa, &Expr::minus_one());
    let bb = factor_light(&neg(&mul2(b, &iaa))).to_expr();
    let dd = factor_terms_radicals(&mul2(&sqrt_factored(&d), &iaa));
    let (mut r0, mut r1) = (sub(&bb, &dd), add(&[bb, dd]));
    if a.as_number().is_some_and(|n| n.is_negative()) {
        std::mem::swap(&mut r0, &mut r1);
    }
    vec![r0, r1]
}

/// Dense recursive representation (SymPy's `dmp`) used as a sort key.
#[derive(Clone, Debug, PartialEq)]
enum DKey {
    I(Int),
    L(Vec<DKey>),
}

fn dkey_cmp(a: &DKey, b: &DKey) -> Ordering {
    match (a, b) {
        (DKey::I(x), DKey::I(y)) => x.cmp(y),
        (DKey::L(x), DKey::L(y)) => {
            for (p, q) in x.iter().zip(y) {
                let c = dkey_cmp(p, q);
                if c != Ordering::Equal {
                    return c;
                }
            }
            x.len().cmp(&y.len())
        }
        (DKey::I(_), DKey::L(_)) => Ordering::Less,
        (DKey::L(_), DKey::I(_)) => Ordering::Greater,
    }
}

fn dmp_zero(u: usize) -> DKey {
    let mut r = DKey::L(vec![]);
    for _ in 0..u {
        r = DKey::L(vec![r]);
    }
    r
}

fn dmp(p: &MPoly, k: usize) -> DKey {
    let n = p.nvars;
    if k == n {
        return DKey::I(p.as_constant().unwrap_or(Int::ZERO));
    }
    if p.is_zero() {
        return dmp_zero(n - 1 - k);
    }
    let cs = p.coeffs_in(k);
    DKey::L(
        cs.iter()
            .rev()
            .map(|c| if k + 1 == n || !c.is_zero() { dmp(c, k + 1) } else { dmp_zero(n - 2 - k) })
            .collect(),
    )
}

/// Order of the linear factors `lambda - r` in SymPy's `factor_list`.
pub fn linear_factor_cmp(a: &Expr, b: &Expr) -> Ordering {
    let lam = Expr::symbol("_x0");
    let fa = sub(&lam, a);
    let fb = sub(&lam, b);
    let Some((_, ps)) = PolyCtx::from_exprs(&[fa, fb], false) else {
        return Ordering::Equal;
    };
    let norm = |q: &super::ratfunc::QPoly| {
        let c = q.p.content();
        let p = if c.is_zero() || c.is_one() { q.p.clone() } else { q.p.divexact_int(&c) };
        p.normalize_sign()
    };
    let (ka, kb) = (dmp(&norm(&ps[0]), 0), dmp(&norm(&ps[1]), 0));
    let la = match &ka {
        DKey::L(v) => v.len(),
        _ => 0,
    };
    let lb = match &kb {
        DKey::L(v) => v.len(),
        _ => 0,
    };
    la.cmp(&lb).then_with(|| dkey_cmp(&ka, &kb))
}

fn merge(rs: Vec<(Expr, usize)>) -> Vec<(Expr, usize)> {
    let mut out: Vec<(Expr, usize)> = Vec::new();
    for (r, k) in rs {
        if let Some(x) = out.iter_mut().find(|(e, _)| *e == r) {
            x.1 += k;
        } else {
            out.push((r, k));
        }
    }
    out
}

/// SymPy's `roots()` for a polynomial (coefficients highest degree first)
/// with symbolic coefficients of degree <= 2 after removing zero roots.
pub fn roots_symbolic(c: &[Expr]) -> MResult<Vec<(Expr, usize)>> {
    let mut c: Vec<Expr> = c.to_vec();
    while c.len() > 1 && c[0].is_zero() {
        c.remove(0);
    }
    let mut zeros = 0;
    while c.len() > 1 && c.last().unwrap().is_zero() {
        c.pop();
        zeros += 1;
    }
    let mut out: Vec<(Expr, usize)> = Vec::new();
    let n = c.len() - 1;
    if n > 0 {
        // preprocess_roots: clear denominators, primitive part
        let pc = preprocess(&c);
        let c = pc.as_deref().unwrap_or(&c);
        match n {
            1 => out.push((factor_light(&neg(&mul2(&c[1], &pow(&c[0], &Expr::minus_one())))).to_expr(), 1)),
            2 => {
                let (a, b, cc) = (&c[0], &c[1], &c[2]);
                if b.is_zero() {
                    for r in roots_quadratic_sym(a, b, cc) {
                        out.push((r, 1));
                    }
                } else {
                    let disc = cancel(&sub(&pow(b, &Expr::int(2)), &mul(&[Expr::int(4), a.clone(), cc.clone()])));
                    let iaa = pow(&mul2(&Expr::int(2), a), &Expr::minus_one());
                    let sq = if disc.is_zero() { Some(Expr::zero()) } else { expr_sqrt(&disc) };
                    match sq {
                        Some(s) if !s.is_zero() => {
                            // two linear factors (factor_list order)
                            let bb = neg(&mul2(b, &iaa));
                            let r0 = cancel(&sub(&bb, &mul2(&s, &iaa)));
                            let r1 = cancel(&add(&[bb, mul2(&s, &iaa)]));
                            let mut rs = vec![r0, r1];
                            rs.sort_by(linear_factor_cmp);
                            for r in rs {
                                out.push((r, 1));
                            }
                        }
                        _ => {
                            for r in roots_quadratic_sym(a, b, cc) {
                                out.push((r, 1));
                            }
                        }
                    }
                }
            }
            _ => {
                return Err(MatrixError::NotImplemented(
                    "roots of symbolic polynomials of degree > 2 are not implemented".into(),
                ))
            }
        }
    }
    if zeros > 0 {
        out.push((Expr::zero(), zeros));
    }
    Ok(merge(out))
}

/// Clear denominators and remove the polynomial content of the
/// coefficients (SymPy's `preprocess_roots`).
fn preprocess(c: &[Expr]) -> Option<Vec<Expr>> {
    let mut parts = Vec::with_capacity(2 * c.len());
    for e in c {
        let (n, d) = as_numer_denom(e);
        parts.push(n);
        parts.push(d);
    }
    let (ctx, ps) = PolyCtx::from_exprs(&parts, false)?;
    let nv = ctx.nvars();
    let ks: Vec<KF> = ps
        .chunks(2)
        .map(|nd| KF::new(nd[0].p.mul_int(&nd[1].den), nd[1].p.mul_int(&nd[0].den)))
        .collect();
    // common denominator
    let mut l = MPoly::one(nv);
    for k in &ks {
        if l.divexact(&k.d).is_none() {
            let g = super::mpoly::gcd(&l, &k.d);
            l = l.mul(&k.d.divexact(&g).unwrap());
        }
    }
    let polys: Vec<MPoly> = ks.iter().map(|k| k.n.mul(&l.divexact(&k.d).unwrap())).collect();
    let mut g = MPoly::zero(nv);
    for p in &polys {
        g = super::mpoly::gcd(&g, p);
    }
    if g.is_zero() {
        return None;
    }
    Some(polys.iter().map(|p| ctx.to_expr(&p.divexact(&g).unwrap())).collect())
}

/// Square root of a rational function expression when it is a square.
fn expr_sqrt(e: &Expr) -> Option<Expr> {
    if let Some(n) = e.as_number() {
        let r = Rat::from_number(n)?;
        if r.is_negative() {
            return None;
        }
        let a = r.num.exact_sqrt()?;
        let b = r.den.exact_sqrt()?;
        return Some(Expr::num(Rat::new(a, b).to_number()));
    }
    let (n, d) = as_numer_denom(e);
    let (ctx, ps) = PolyCtx::from_exprs(&[n, d], false)?;
    let k = KF::new(ps[0].p.mul_int(&ps[1].den), ps[1].p.mul_int(&ps[0].den));
    let s = sqrt_k(&k)?;
    Some(s.to_expr(&ctx))
}

// ----------------------------------------------------------------------
// eigenvectors over K
// ----------------------------------------------------------------------

enum KFactor {
    Linear(KF, usize),
    /// monic `l^2 + p l + q` (irreducible) and its primitive polynomial
    /// form (expressions, highest degree first) for `roots()`
    Quad(KF, KF, usize),
}

fn horner(p: &[KF], x: &KF) -> KF {
    let mut acc = p[0].zero_like();
    for c in p {
        acc = acc.mul(x).add(c);
    }
    acc
}

fn deflate(p: &[KF], r: &KF) -> Vec<KF> {
    let mut out = Vec::with_capacity(p.len() - 1);
    let mut acc = p[0].zero_like();
    for c in &p[..p.len() - 1] {
        acc = acc.mul(r).add(c);
        out.push(acc.clone());
    }
    out
}

/// Split a monic polynomial over K (highest degree first).
fn factor_k(p: &[KF], candidates: &[KF]) -> Option<Vec<KFactor>> {
    let nv = p[0].n.nvars;
    let mut p: Vec<KF> = p.to_vec();
    let mut out: Vec<KFactor> = Vec::new();
    let push_lin = |r: KF, out: &mut Vec<KFactor>| {
        for f in out.iter_mut() {
            if let KFactor::Linear(x, k) = f {
                if *x == r {
                    *k += 1;
                    return;
                }
            }
        }
        out.push(KFactor::Linear(r, 1));
    };
    while p.len() > 1 && p.last().unwrap().is_zero() {
        p.pop();
        push_lin(KF::zero(nv), &mut out);
    }
    while p.len() > 3 {
        let mut found = false;
        for c in candidates {
            if horner(&p, c).is_zero() {
                p = deflate(&p, c);
                push_lin(c.clone(), &mut out);
                found = true;
                break;
            }
        }
        if !found {
            return None;
        }
    }
    match p.len() {
        1 => {}
        2 => push_lin(p[1].neg().div(&p[0]), &mut out),
        3 => {
            let (a, b, c) = (&p[0], &p[1], &p[2]);
            let disc = b.mul(b).sub(&KF::from_poly(MPoly::constant(nv, Int::from(4))).mul(a).mul(c));
            let two_a = a.mul(&KF::from_poly(MPoly::constant(nv, Int::from(2))));
            if disc.is_zero() {
                let r = b.neg().div(&two_a);
                push_lin(r.clone(), &mut out);
                push_lin(r, &mut out);
            } else if let Some(s) = sqrt_k(&disc) {
                push_lin(b.neg().sub(&s).div(&two_a), &mut out);
                push_lin(b.neg().add(&s).div(&two_a), &mut out);
            } else {
                out.push(KFactor::Quad(b.div(a), c.div(a), 1));
            }
        }
        _ => unreachable!(),
    }
    Some(out)
}

/// Nullspace over K (standard basis from the reduced row echelon form).
fn nullspace_k(rows: usize, cols: usize, mut a: Vec<KF>) -> Vec<Vec<KF>> {
    let nv = a[0].n.nvars;
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p) = (r..rows).find(|&i| !a[i * cols + c].is_zero()) else {
            continue;
        };
        for k in 0..cols {
            a.swap(p * cols + k, r * cols + k);
        }
        let inv = a[r * cols + c].inv();
        for k in 0..cols {
            a[r * cols + k] = a[r * cols + k].mul(&inv);
        }
        for i in 0..rows {
            if i == r || a[i * cols + c].is_zero() {
                continue;
            }
            let f = a[i * cols + c].clone();
            for k in 0..cols {
                let v = a[i * cols + k].sub(&f.mul(&a[r * cols + k]));
                a[i * cols + k] = v;
            }
        }
        pivots.push(c);
        r += 1;
    }
    let mut basis = Vec::new();
    for free in (0..cols).filter(|c| !pivots.contains(c)) {
        let mut v = vec![KF::zero(nv); cols];
        v[free] = KF::one(nv);
        for (pr, &pc) in pivots.iter().enumerate() {
            v[pc] = v[pc].sub(&a[pr * cols + free]);
        }
        basis.push(v);
    }
    basis
}

/// Elements `a + b*l` of `K[l]/(l^2 + p l + q)`.
#[derive(Clone, Debug)]
struct QK {
    a: KF,
    b: KF,
}

struct ExtK {
    p: KF,
    q: KF,
}

impl ExtK {
    fn mul(&self, x: &QK, y: &QK) -> QK {
        let bd = x.b.mul(&y.b);
        QK { a: x.a.mul(&y.a).sub(&bd.mul(&self.q)), b: x.a.mul(&y.b).add(&x.b.mul(&y.a)).sub(&bd.mul(&self.p)) }
    }
    fn sub(&self, x: &QK, y: &QK) -> QK {
        QK { a: x.a.sub(&y.a), b: x.b.sub(&y.b) }
    }
    fn inv(&self, x: &QK) -> QK {
        let n = x.a.mul(&x.a).sub(&x.a.mul(&x.b).mul(&self.p)).add(&x.b.mul(&x.b).mul(&self.q));
        let ni = n.inv();
        QK { a: x.a.sub(&x.b.mul(&self.p)).mul(&ni), b: x.b.neg().mul(&ni) }
    }
    fn is_zero(x: &QK) -> bool {
        x.a.is_zero() && x.b.is_zero()
    }
}

fn nullspace_ext(f: &ExtK, rows: usize, cols: usize, mut a: Vec<QK>) -> Vec<Vec<QK>> {
    let nv = f.p.n.nvars;
    let zero = QK { a: KF::zero(nv), b: KF::zero(nv) };
    let one = QK { a: KF::one(nv), b: KF::zero(nv) };
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p) = (r..rows).find(|&i| !ExtK::is_zero(&a[i * cols + c])) else {
            continue;
        };
        for k in 0..cols {
            a.swap(p * cols + k, r * cols + k);
        }
        let inv = f.inv(&a[r * cols + c]);
        for k in 0..cols {
            a[r * cols + k] = f.mul(&a[r * cols + k], &inv);
        }
        for i in 0..rows {
            if i == r || ExtK::is_zero(&a[i * cols + c]) {
                continue;
            }
            let x = a[i * cols + c].clone();
            for k in 0..cols {
                let v = f.sub(&a[i * cols + k], &f.mul(&x, &a[r * cols + k]));
                a[i * cols + k] = v;
            }
        }
        pivots.push(c);
        r += 1;
    }
    let mut basis = Vec::new();
    for free in (0..cols).filter(|c| !pivots.contains(c)) {
        let mut v = vec![zero.clone(); cols];
        v[free] = one.clone();
        for (pr, &pc) in pivots.iter().enumerate() {
            v[pc] = f.sub(&v[pc], &a[pr * cols + free]);
        }
        basis.push(v);
    }
    basis
}

/// SymPy's `_eigenvects_DOM` for a matrix over `ZZ(gens)`/`QQ(gens)`.
/// `None` when the spectrum is not "easy".
pub fn eigenvects_field(m: &Matrix, ctx: &PolyCtx) -> Option<MResult<Vec<super::eigen::EigenVect>>> {
    let n = m.rows();
    let nv = ctx.nvars();
    let a: Vec<KF> = m.data().iter().map(|e| KF::from_expr(ctx, e)).collect::<Option<_>>()?;
    let one = KF::one(nv);
    let cp = berkowitz(&a, n, &one);
    let mut cands: Vec<KF> = Vec::new();
    for i in 0..n {
        let d = a[i * n + i].clone();
        if !cands.contains(&d) {
            cands.push(d);
        }
    }
    let factors = factor_k(&cp, &cands)?;
    let mut out: Vec<super::eigen::EigenVect> = Vec::new();
    for f in factors {
        match f {
            KFactor::Linear(r, k) => {
                let mut sh = a.clone();
                for i in 0..n {
                    sh[i * n + i] = sh[i * n + i].sub(&r);
                }
                let basis = nullspace_k(n, n, sh);
                let vs = basis
                    .iter()
                    .map(|v| Matrix::column(v.iter().map(|x| x.to_expr(ctx)).collect()))
                    .collect();
                out.push((r.to_expr(ctx), k, vs));
            }
            KFactor::Quad(p, q, k) => {
                let ext = ExtK { p: p.clone(), q: q.clone() };
                let mut sh: Vec<QK> = a.iter().map(|x| QK { a: x.clone(), b: KF::zero(nv) }).collect();
                for i in 0..n {
                    sh[i * n + i].b = sh[i * n + i].b.sub(&KF::one(nv));
                }
                let basis = nullspace_ext(&ext, n, n, sh);
                // roots of the minimal polynomial l^2 + p l + q
                let coeffs = vec![Expr::one(), p.to_expr(ctx), q.to_expr(ctx)];
                let rts = match roots_symbolic(&coeffs) {
                    Ok(r) => r,
                    Err(e) => return Some(Err(e)),
                };
                for (root, _) in rts {
                    let vs = basis
                        .iter()
                        .map(|v| {
                            Matrix::column(
                                v.iter()
                                    .map(|x| add(&[mul2(&x.b.to_expr(ctx), &root), x.a.to_expr(ctx)]))
                                    .collect(),
                            )
                        })
                        .collect();
                    out.push((root, k, vs));
                }
            }
        }
    }
    out.sort_by(|x, y| crate::sortkey::cmp_default(&x.0, &y.0));
    Some(Ok(out))
}
