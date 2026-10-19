//! Rational function helpers for the matrix algorithms: conversion between
//! expressions and [`MPoly`] (with SymPy's choice and order of generators),
//! `cancel`, SymPy's `dotprodsimp` intermediate simplification and a light
//! version of the `is_zero` assumption query.

use std::cmp::Ordering;

use rustc_hash::FxHashMap;

use super::int::{Int, Rat};
use super::mpoly::{gcd, MPoly, Mono};
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, add_from_terms, could_extract_minus_sign, div, mul, mul_factor, neg, pow};
use crate::printer::to_str;

// ----------------------------------------------------------------------
// generators
// ----------------------------------------------------------------------

/// SymPy's `_sort_gens` key.
fn gen_key(g: &Expr) -> (u32, String, u64) {
    let s = to_str(g);
    let digits = s.bytes().rev().take_while(|b| b.is_ascii_digit()).count();
    let (name, idx) = s.split_at(s.len() - digits);
    let index: u64 = if idx.is_empty() { 0 } else { idx.parse().unwrap_or(u64::MAX) };
    let order = if name.len() == 1 {
        let c = name.as_bytes()[0];
        match c {
            b'a'..=b'o' => 301 + (c - b'a') as u32,
            b'p'..=b'w' => 216 + (c - b'p') as u32,
            b'x'..=b'z' => 124 + (c - b'x') as u32,
            _ => 1000,
        }
    } else {
        1000
    };
    (order, name.to_string(), index)
}

pub fn cmp_gens(a: &Expr, b: &Expr) -> Ordering {
    gen_key(a).cmp(&gen_key(b))
}

/// SymPy's `decompose_power`: `factor = base**exp` with an integer `exp`.
fn decompose_power(f: &Expr, opaque_funcs: bool) -> (Expr, i64) {
    let (b, e) = match f.kind() {
        Kind::Pow(b, e) => (b.clone(), e.clone()),
        Kind::Func(Func::Exp, args) if !opaque_funcs => (Expr::e(), args[0].clone()),
        _ => return (f.clone(), 1),
    };
    let rebuild = |b: &Expr, e: &Expr| -> Expr {
        if b.is_constant(Constant::E) {
            crate::functions::exp(e)
        } else {
            pow(b, e)
        }
    };
    if let Some(n) = e.as_number() {
        if let Number::Float(_) = n {
            return (f.clone(), 1);
        }
        let (p, q) = n.numer_denom().unwrap();
        let p = match Int::from_big(p).to_i64() {
            Some(v) => v,
            None => return (f.clone(), 1),
        };
        let base = if q == num_bigint::BigInt::from(1) {
            b
        } else {
            rebuild(&b, &Expr::num(Number::rational(num_bigint::BigInt::from(1), q).unwrap()))
        };
        return (base, p);
    }
    // symbolic exponent: pull out the rational coefficient
    let (c, tail) = crate::ops::as_coeff_mul(&e);
    if !c.is_rational() || c.is_one() {
        return (f.clone(), 1);
    }
    if c.is_minus_one() {
        return (rebuild(&b, &tail), -1);
    }
    let (p, q) = c.numer_denom().unwrap();
    let p = match Int::from_big(p).to_i64() {
        Some(v) => v,
        None => return (f.clone(), 1),
    };
    let t = mul(&[Expr::num(Number::rational(num_bigint::BigInt::from(1), q).unwrap()), tail]);
    (rebuild(&b, &t), p)
}

/// A rational-coefficient polynomial `p / den` over the context generators.
#[derive(Clone, Debug)]
pub struct QPoly {
    pub p: MPoly,
    pub den: Int,
}

/// Generators of a set of polynomial expressions.
#[derive(Clone, Debug)]
pub struct PolyCtx {
    pub gens: Vec<Expr>,
    index: FxHashMap<Expr, usize>,
    opaque_funcs: bool,
}

type RawTerm = (Rat, Vec<(Expr, i64)>);

impl PolyCtx {
    /// Collect generators from the (expanded) expressions. Returns `None`
    /// when a coefficient is not rational (floats).
    pub fn from_exprs(exprs: &[Expr], opaque_funcs: bool) -> Option<(PolyCtx, Vec<QPoly>)> {
        let mut raws: Vec<Vec<RawTerm>> = Vec::with_capacity(exprs.len());
        let mut gens: Vec<Expr> = Vec::new();
        let mut seen: FxHashMap<Expr, ()> = FxHashMap::default();
        for e in exprs {
            let ex = expand(e);
            let terms = raw_terms(&ex, opaque_funcs)?;
            for (_, fs) in &terms {
                for (g, _) in fs {
                    if seen.insert(g.clone(), ()).is_none() {
                        gens.push(g.clone());
                    }
                }
            }
            raws.push(terms);
        }
        gens.sort_by(cmp_gens);
        let index: FxHashMap<Expr, usize> = gens.iter().enumerate().map(|(i, g)| (g.clone(), i)).collect();
        let ctx = PolyCtx { gens, index, opaque_funcs };
        let polys = raws.iter().map(|t| ctx.from_raw(t)).collect::<Option<Vec<_>>>()?;
        Some((ctx, polys))
    }

    pub fn nvars(&self) -> usize {
        self.gens.len()
    }

    fn from_raw(&self, terms: &[RawTerm]) -> Option<QPoly> {
        let n = self.gens.len();
        let mut den = Int::ONE;
        for (c, _) in terms {
            den = den.lcm(&c.den);
        }
        let mut out: Vec<(Mono, Int)> = Vec::with_capacity(terms.len());
        for (c, fs) in terms {
            let mut m = Mono::from_elem(0, n);
            for (g, e) in fs {
                let i = *self.index.get(g)?;
                if *e < 0 {
                    return None;
                }
                m[i] += *e as u32;
            }
            let k = c.num.mul(&den.divexact(&c.den));
            out.push((m, k));
        }
        Some(QPoly { p: MPoly::from_terms(n, out), den })
    }

    /// Convert another expression with the same generators.
    pub fn to_qpoly(&self, e: &Expr) -> Option<QPoly> {
        let ex = expand(e);
        let terms = raw_terms(&ex, self.opaque_funcs)?;
        self.from_raw(&terms)
    }

    /// Expression of an integer polynomial.
    pub fn to_expr(&self, p: &MPoly) -> Expr {
        self.to_expr_scaled(p, &Int::ONE)
    }

    /// Expression of `p / den`.
    pub fn to_expr_scaled(&self, p: &MPoly, den: &Int) -> Expr {
        if self.gens.iter().all(Expr::is_symbol) {
            // Distinct exponent vectors over distinct symbols give distinct
            // coefficient-free monomials: the sum can skip collection.
            let mut constant = Number::ZERO;
            let mut terms: Vec<(Expr, Number)> = Vec::with_capacity(p.terms.len());
            for (m, c) in &p.terms {
                let coeff = if den.is_one() { c.to_number() } else { Rat::new(c.clone(), den.clone()).to_number() };
                let fs: Vec<Expr> = m
                    .iter()
                    .enumerate()
                    .filter(|(_, e)| **e > 0)
                    .map(|(i, e)| if *e == 1 { self.gens[i].clone() } else { pow(&self.gens[i], &Expr::int(*e as i64)) })
                    .collect();
                if fs.is_empty() {
                    constant = coeff;
                } else {
                    terms.push((mul(&fs), coeff));
                }
            }
            return add_from_terms(constant, terms);
        }
        let mut terms: Vec<Expr> = Vec::with_capacity(p.terms.len());
        for (m, c) in &p.terms {
            let coeff = if den.is_one() {
                c.to_number()
            } else {
                Rat::new(c.clone(), den.clone()).to_number()
            };
            let mut fs: Vec<Expr> = Vec::with_capacity(m.len() + 1);
            fs.push(Expr::num(coeff));
            for (i, e) in m.iter().enumerate() {
                if *e > 0 {
                    fs.push(pow(&self.gens[i], &Expr::int(*e as i64)));
                }
            }
            terms.push(mul(&fs));
        }
        add(&terms)
    }
}

fn raw_terms(e: &Expr, opaque_funcs: bool) -> Option<Vec<RawTerm>> {
    let mut out = Vec::new();
    match e.kind() {
        Kind::Add(c, terms) => {
            for (t, k) in terms {
                let mut r = raw_term(t, opaque_funcs)?;
                r.0 = r.0.mul(&Rat::from_number(k)?);
                out.push(r);
            }
            if !c.is_zero() {
                out.push((Rat::from_number(c)?, Vec::new()));
            }
        }
        _ => {
            if !e.is_zero() {
                out.push(raw_term(e, opaque_funcs)?);
            }
        }
    }
    Some(out)
}

fn raw_term(t: &Expr, opaque_funcs: bool) -> Option<RawTerm> {
    let mut coeff = Rat::one();
    let mut fs: Vec<(Expr, i64)> = Vec::new();
    let push = |f: &Expr, fs: &mut Vec<(Expr, i64)>| {
        let (mut b, mut e) = decompose_power(f, opaque_funcs);
        if e < 0 {
            e = -e;
            b = pow(&b, &Expr::minus_one());
        }
        if let Some(x) = fs.iter_mut().find(|(g, _)| *g == b) {
            x.1 += e;
        } else {
            fs.push((b, e));
        }
    };
    match t.kind() {
        Kind::Num(n) => {
            coeff = Rat::from_number(n)?;
        }
        Kind::Mul(c, facs) => {
            coeff = Rat::from_number(c)?;
            for (b, x) in facs {
                push(&mul_factor(b, x), &mut fs);
            }
        }
        _ => push(t, &mut fs),
    }
    // floats hidden inside generators are fine; float coefficients are not
    Some((coeff, fs))
}

// ----------------------------------------------------------------------
// cancel
// ----------------------------------------------------------------------

/// Reduce `num/den` (integer polynomials) the way SymPy's
/// `PolyElement.cancel` does: divide by the gcd and make the leading
/// coefficient of the denominator positive.
pub fn cancel_polys(num: &MPoly, den: &MPoly) -> (MPoly, MPoly) {
    if num.is_zero() {
        return (num.clone(), MPoly::one(num.nvars));
    }
    let g = gcd(num, den);
    let (mut p, mut q) = if g.is_one() {
        (num.clone(), den.clone())
    } else {
        (num.divexact(&g).expect("gcd"), den.divexact(&g).expect("gcd"))
    };
    if q.lc().is_negative() {
        p = p.neg();
        q = q.neg();
    }
    (p, q)
}

/// `(p/cp) / (q/cq)` in lowest terms, as integer polynomials `(P, Q)`.
pub fn cancel_qpolys(f: &QPoly, g: &QPoly) -> (MPoly, MPoly) {
    let (p, q) = cancel_polys(&f.p, &g.p);
    // value = (p * g.den) / (q * f.den)
    let cg = g.den.gcd(&f.den);
    let cp = g.den.divexact(&cg);
    let cq = f.den.divexact(&cg);
    (p.mul_int(&cp), q.mul_int(&cq))
}

/// SymPy's `cancel` for rational functions with rational coefficients.
pub fn cancel(e: &Expr) -> Expr {
    cancel_opts(e, false)
}

/// SymPy's `as_numer_denom` (recursive through nested fractions).
pub fn as_numer_denom(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Num(n) => (Expr::num(n.numer()), Expr::num(n.denom())),
        Kind::Add(..) => {
            let parts: Vec<(Expr, Expr)> = e.args().iter().map(as_numer_denom).collect();
            if parts.iter().all(|(_, d)| d.is_one()) {
                return (e.clone(), Expr::one());
            }
            let mut dens: Vec<Expr> = Vec::new();
            for (_, d) in &parts {
                if !d.is_one() && !dens.contains(d) {
                    dens.push(d.clone());
                }
            }
            let den = mul(&dens);
            let nums: Vec<Expr> = parts
                .iter()
                .map(|(n, d)| {
                    let mut f = vec![n.clone()];
                    for dd in &dens {
                        if dd != d {
                            f.push(dd.clone());
                        }
                    }
                    mul(&f)
                })
                .collect();
            (add(&nums), den)
        }
        Kind::Mul(..) => {
            let mut num = Vec::new();
            let mut den = Vec::new();
            for f in e.args() {
                let (n, d) = as_numer_denom(&f);
                num.push(n);
                den.push(d);
            }
            (mul(&num), mul(&den))
        }
        Kind::Pow(b, x) => {
            if let Some(k) = x.as_number() {
                if k.is_integer() {
                    let (n, d) = as_numer_denom(b);
                    return if k.is_negative() {
                        let m = Expr::num(k.neg());
                        (pow(&d, &m), pow(&n, &m))
                    } else {
                        (pow(&n, x), pow(&d, x))
                    };
                }
                if k.is_negative() {
                    return (Expr::one(), pow(b, &Expr::num(k.neg())));
                }
                return (e.clone(), Expr::one());
            }
            if let Kind::Mul(c, _) = x.kind() {
                if c.is_negative() {
                    return (Expr::one(), pow(b, &neg(x)));
                }
            }
            (e.clone(), Expr::one())
        }
        Kind::Func(Func::Exp, args) => {
            if let Kind::Mul(c, _) = args[0].kind() {
                if c.is_negative() {
                    return (Expr::one(), crate::functions::exp(&neg(&args[0])));
                }
            }
            (e.clone(), Expr::one())
        }
        _ => (e.clone(), Expr::one()),
    }
}

pub fn cancel_opts(e: &Expr, opaque_funcs: bool) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => {}
        _ => return e.clone(),
    }
    let (n, d) = as_numer_denom(e);
    let Some((ctx, ps)) = PolyCtx::from_exprs(&[n, d], opaque_funcs) else {
        return e.clone();
    };
    if ctx.nvars() == 0 {
        return expand(e);
    }
    let (p, q) = cancel_qpolys(&ps[0], &ps[1]);
    let pe = ctx.to_expr(&p);
    if q.is_one() {
        return pe;
    }
    div(&pe, &ctx.to_expr(&q))
}

// ----------------------------------------------------------------------
// dotprodsimp
// ----------------------------------------------------------------------

/// SymPy's `fraction(a)` (non-exact) for a product: numerator and
/// denominator, without looking inside sums.
fn fraction(a: &Expr) -> (Expr, Expr) {
    let mut num = Vec::new();
    let mut den = Vec::new();
    for t in a.args() {
        match t.kind() {
            Kind::Num(n) if n.is_rational() && !n.is_integer() => {
                if !n.numer().is_one() {
                    num.push(Expr::num(n.numer()));
                }
                den.push(Expr::num(n.denom()));
            }
            Kind::Pow(b, x) => match x.as_number() {
                Some(n) if n.is_negative() => {
                    den.push(if n.is_minus_one() { b.clone() } else { pow(b, &Expr::num(n.neg())) })
                }
                Some(_) => num.push(t.clone()),
                None => match x.kind() {
                    Kind::Mul(c, _) if c.is_negative() => den.push(pow(b, &neg(x))),
                    _ => num.push(t.clone()),
                },
            },
            Kind::Func(Func::Exp, args) => match args[0].kind() {
                Kind::Mul(c, _) if c.is_negative() => den.push(crate::functions::exp(&neg(&args[0]))),
                _ => num.push(t.clone()),
            },
            _ => num.push(t.clone()),
        }
    }
    (mul(&num), mul(&den))
}

/// SymPy's `count_ops_alg` (from `dotprodsimp`).
pub fn count_ops_alg(e: &Expr) -> (usize, bool) {
    let mut ops = 0usize;
    let mut ratfunc = false;
    let mut stack = vec![e.clone()];
    while let Some(a) = stack.pop() {
        match a.kind() {
            Kind::Num(n) if n.is_rational() => {
                if !n.is_one() {
                    ops += n.is_negative() as usize + (!n.is_integer()) as usize;
                }
            }
            Kind::Mul(c, _) => {
                let mut a = a.clone();
                if could_extract_minus_sign(&a) {
                    ops += 1;
                    let _ = c;
                    a = neg(&a);
                }
                let (n, d) = fraction(&a);
                if n.is_integer_number() {
                    ops += 1 + n.as_number().unwrap().is_negative() as usize;
                    stack.push(d);
                } else if !d.is_one() {
                    if !d.is_integer_number() {
                        stack.push(d);
                        ratfunc = true;
                    }
                    ops += 1;
                    stack.push(n);
                } else {
                    let args = a.args();
                    ops += args.len().saturating_sub(1);
                    stack.extend(args);
                }
            }
            Kind::Add(..) => {
                let args = a.args();
                let la = args.len();
                let mut negs = 0;
                for ai in args {
                    if could_extract_minus_sign(&ai) {
                        negs += 1;
                        stack.push(neg(&ai));
                    } else {
                        stack.push(ai);
                    }
                }
                ops += la - (negs != la) as usize;
            }
            Kind::Pow(b, x) => {
                ops += 1;
                stack.push(b.clone());
                if !ratfunc {
                    ratfunc = match x.as_number() {
                        Some(n) => n.is_negative(),
                        None => true,
                    };
                }
            }
            _ => {}
        }
    }
    (ops, ratfunc)
}

/// SymPy's `dotprodsimp`: expand, then cancel rational functions when the
/// expression is large enough.
pub fn dotprodsimp(e: &Expr) -> Expr {
    let r = dotprodsimp_withsimp(e).0;
    if trace_dps() {
        eprintln!("DPS {} -> {}", e, r);
    }
    r
}

fn trace_dps() -> bool {
    static T: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *T.get_or_init(|| std::env::var_os("SYMOX_DPS_TRACE").is_some())
}

/// `dotprodsimp(e, withsimp=True)`: also returns whether the expression
/// may be simplified again later.
pub fn dotprodsimp_withsimp(e: &Expr) -> (Expr, bool) {
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => {}
        _ => return (e.clone(), false),
    }
    let mut ex = expand(e);
    let mut simplified = ex != *e;
    let (ops, ratfunc) = count_ops_alg(&ex);
    if ops >= 6 {
        if ratfunc {
            let c = cancel_opts(&ex, true);
            if c != ex {
                ex = c;
                simplified = true;
            }
        }
        return (ex, simplified);
    }
    if ops == 5 {
        if let Kind::Add(c, terms) = ex.kind() {
            if c.is_zero() && terms.len() == 2 {
                let a0 = crate::ops::add_term(&terms[0].0, &terms[0].1);
                let a1 = crate::ops::add_term(&terms[1].0, &terms[1].1);
                if last_arg_is_inverse(&a0) && last_arg_is_inverse(&a1) {
                    let t = sympy_together(&ex);
                    if count_ops_alg(&t).0 < ops {
                        return (t, true);
                    }
                    return (ex, simplified);
                }
            }
        }
    }
    (ex, true)
}

/// Approximation of SymPy's `Basic.compare` class order.
fn class_rank(e: &Expr) -> u8 {
    match e.kind() {
        Kind::Num(_) => 0,
        Kind::Const(_) => 1,
        Kind::Sym(_) => 2,
        Kind::Pow(..) => 3,
        Kind::Mul(..) => 4,
        Kind::Add(..) => 5,
        _ => 6,
    }
}

fn basic_cmp(a: &Expr, b: &Expr) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let r = class_rank(a).cmp(&class_rank(b));
    if r != Ordering::Equal {
        return r;
    }
    match (a.kind(), b.kind()) {
        (Kind::Sym(x), Kind::Sym(y)) => x.cmp(y),
        (Kind::Num(x), Kind::Num(y)) => x.cmp_value(y).unwrap_or(Ordering::Equal),
        (Kind::Pow(b1, e1), Kind::Pow(b2, e2)) => basic_cmp(b1, b2).then_with(|| basic_cmp(e1, e2)),
        _ => {
            let (x, y) = (a.args(), b.args());
            x.len().cmp(&y.len()).then_with(|| {
                for (p, q) in x.iter().zip(&y) {
                    let c = basic_cmp(p, q);
                    if c != Ordering::Equal {
                        return c;
                    }
                }
                Ordering::Equal
            })
        }
    }
}

/// Is the last argument of the product (in SymPy's internal order) a
/// power with exponent -1?
fn last_arg_is_inverse(t: &Expr) -> bool {
    let Kind::Mul(_, facs) = t.kind() else {
        return false;
    };
    let args: Vec<Expr> = facs.iter().map(|(b, x)| mul_factor(b, x)).collect();
    let last = args.iter().max_by(|a, b| basic_cmp(a, b));
    matches!(last.map(|l| l.kind()), Some(Kind::Pow(_, x)) if x.is_minus_one())
}

/// `Add.primitive()`: (positive rational content, primitive part).
fn add_primitive(e: &Expr) -> (Number, Expr) {
    let Kind::Add(c, terms) = e.kind() else {
        return (Number::ONE, e.clone());
    };
    let mut g: Option<Number> = if c.is_zero() { None } else { Some(c.abs()) };
    for (_, k) in terms {
        if !k.is_rational() {
            return (Number::ONE, e.clone());
        }
        g = Some(match g {
            None => k.abs(),
            Some(x) => x.gcd(k),
        });
    }
    let g = g.unwrap_or(Number::ONE);
    if g.is_one() || g.is_zero() {
        return (Number::ONE, e.clone());
    }
    let inv = Expr::num(g.inv().unwrap());
    (g, crate::expand::expand(&mul(&[inv, e.clone()])))
}

/// SymPy's `Term`: `coeff * numer / denom` with factor maps.
struct STerm {
    coeff: Number,
    numer: Vec<(Expr, i64)>,
    denom: Vec<(Expr, i64)>,
}

fn fmap_add(m: &mut Vec<(Expr, i64)>, b: Expr, e: i64) {
    if let Some(x) = m.iter_mut().find(|(g, _)| *g == b) {
        x.1 += e;
    } else {
        m.push((b, e));
    }
    m.retain(|(_, e)| *e != 0);
}

fn fmap_gcd(a: &[(Expr, i64)], b: &[(Expr, i64)]) -> Vec<(Expr, i64)> {
    a.iter()
        .filter_map(|(g, e)| b.iter().find(|(h, _)| h == g).map(|(_, f)| (g.clone(), (*e).min(*f))))
        .collect()
}

fn fmap_lcm(a: &[(Expr, i64)], b: &[(Expr, i64)]) -> Vec<(Expr, i64)> {
    let mut out: Vec<(Expr, i64)> = a.to_vec();
    for (g, e) in b {
        match out.iter_mut().find(|(h, _)| h == g) {
            Some(x) => x.1 = x.1.max(*e),
            None => out.push((g.clone(), *e)),
        }
    }
    out
}

fn fmap_expr(m: &[(Expr, i64)]) -> Expr {
    mul(&m.iter().map(|(b, e)| pow(b, &Expr::int(*e))).collect::<Vec<_>>())
}

impl STerm {
    fn new(t: &Expr) -> STerm {
        let (mut coeff, rest) = crate::ops::as_coeff_mul(t);
        let mut numer = Vec::new();
        let mut denom = Vec::new();
        let factors: Vec<Expr> = match rest.kind() {
            Kind::Mul(_, facs) => facs.iter().map(|(b, x)| mul_factor(b, x)).collect(),
            _ if rest.is_one() => vec![],
            _ => vec![rest.clone()],
        };
        for f in factors {
            let (mut base, e) = decompose_power(&f, false);
            if base.is_add() {
                let (cont, prim) = add_primitive(&base);
                if !cont.is_one() {
                    coeff = coeff.mul(&cont.pow_i64(e).unwrap());
                }
                base = prim;
            }
            if e > 0 {
                fmap_add(&mut numer, base, e);
            } else {
                fmap_add(&mut denom, base, -e);
            }
        }
        STerm { coeff, numer, denom }
    }

    /// `self / o` with the common factors of numerator and denominator
    /// removed (`Factors.normal`).
    fn quo(&self, o: &STerm) -> STerm {
        let coeff = self.coeff.div(&o.coeff).unwrap_or(Number::ONE);
        let mut numer = self.numer.clone();
        let mut denom = self.denom.clone();
        for (b, e) in &o.denom {
            fmap_add(&mut numer, b.clone(), *e);
        }
        for (b, e) in &o.numer {
            fmap_add(&mut denom, b.clone(), *e);
        }
        // normal(): cancel common bases
        let mut n2 = Vec::new();
        for (b, e) in numer {
            if let Some(d) = denom.iter_mut().find(|(g, _)| *g == b) {
                let m = e.min(d.1);
                d.1 -= m;
                if e - m != 0 {
                    n2.push((b, e - m));
                }
            } else {
                n2.push((b, e));
            }
        }
        denom.retain(|(_, e)| *e != 0);
        STerm { coeff, numer: n2, denom }
    }
}

/// SymPy's `gcd_terms` for a list of terms.
fn gcd_terms(terms: &[Expr]) -> Expr {
    let ts: Vec<STerm> = terms.iter().filter(|t| !t.is_zero()).map(STerm::new).collect();
    if ts.is_empty() {
        return Expr::zero();
    }
    let (cont_c, cont_n, cont_d, numer, denom) = if ts.len() == 1 {
        let t = &ts[0];
        (t.coeff.clone(), t.numer.clone(), t.denom.clone(), Expr::one(), Expr::one())
    } else {
        let mut cc = ts[0].coeff.abs();
        let mut cn = ts[0].numer.clone();
        let mut cd = ts[0].denom.clone();
        for t in &ts[1..] {
            cc = cc.gcd(&t.coeff);
            cn = fmap_gcd(&cn, &t.numer);
            cd = fmap_gcd(&cd, &t.denom);
        }
        let cont = STerm { coeff: cc.clone(), numer: cn.clone(), denom: cd.clone() };
        let qs: Vec<STerm> = ts.iter().map(|t| t.quo(&cont)).collect();
        let mut den = qs[0].denom.clone();
        for t in &qs[1..] {
            den = fmap_lcm(&den, &t.denom);
        }
        let mut numers = Vec::new();
        for t in &qs {
            // numer * (den / t.denom)
            let mut f = t.numer.clone();
            for (b, e) in &den {
                let have = t.denom.iter().find(|(g, _)| g == b).map(|(_, x)| *x).unwrap_or(0);
                if e - have != 0 {
                    fmap_add(&mut f, b.clone(), e - have);
                }
            }
            numers.push(mul(&[Expr::num(t.coeff.clone()), fmap_expr(&f)]));
        }
        (cc, cn, cd, add(&numers), fmap_expr(&den))
    };
    let (pc, numer) = if numer.is_add() { add_primitive(&numer) } else { (Number::ONE, numer) };
    let coeff = cont_c.mul(&pc);
    let factors = mul(&[fmap_expr(&cont_n), pow(&fmap_expr(&cont_d), &Expr::minus_one())]);
    let body = mul(&[factors, numer, pow(&denom, &Expr::minus_one())]);
    mul(&[Expr::num(coeff), body])
}

/// SymPy's `together` (combination by `gcd_terms`).
pub fn sympy_together(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Add(..) => {
            let parts: Vec<Expr> = e.args().iter().map(sympy_together).collect();
            gcd_terms(&parts)
        }
        Kind::Pow(b, x) => pow(&sympy_together(b), x),
        Kind::Mul(..) => mul(&e.args().iter().map(sympy_together).collect::<Vec<_>>()),
        _ => e.clone(),
    }
}

// ----------------------------------------------------------------------
// zero testing
// ----------------------------------------------------------------------

/// Approximation of SymPy's `expr.is_zero` (three-valued).
pub fn is_zero_assume(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_zero()),
        Kind::Const(Constant::NaN) => None,
        Kind::Const(_) => Some(false),
        Kind::Sym(_) => crate::assumptions::is_zero(e),
        Kind::Bool(_) => None,
        Kind::Func(Func::Exp, _) => Some(false),
        Kind::Pow(b, x) => {
            if is_zero_assume(b)==Some(true) && crate::assumptions::is_positive(x)==Some(true) {
                return Some(true);
            }
            if is_zero_assume(b) == Some(false) {
                if x.is_number() || b.is_number() || b.is_constant(Constant::E) {
                    return Some(false);
                }
            }
            if e.mask() == 0 {
                return numeric_zero_exact(e);
            }
            None
        }
        Kind::Add(constant,terms)=>{
            if e.mask()==0 {return numeric_zero_exact(e);}
            if matches!(constant,Number::Float(value) if !value.is_finite()) {return None;}
            let mut nonzero=usize::from(!constant.is_zero());
            for (term,coefficient) in terms {
                if coefficient.is_zero() {continue;}
                if matches!(coefficient,Number::Float(value) if !value.is_finite()) {return None;}
                match is_zero_assume(term) {
                    Some(true)=>{},Some(false)=>nonzero+=1,None=>return None,
                }
            }
            if nonzero==0 {Some(true)} else if nonzero==1 {Some(false)} else {None}
        }
        Kind::Mul(c, facs) => {
            if e.has_nonfinite() {return None;}
            if c.is_zero() {
                return Some(true);
            }
            let mut all_nonzero = true;
            for (b, x) in facs {
                match is_zero_assume(&mul_factor(b, x)) {
                    Some(false) => {}
                    Some(true) => return Some(true),
                    None => all_nonzero = false,
                }
            }
            if all_nonzero {
                Some(false)
            } else {
                None
            }
        }
        _ => {
            if e.mask() == 0 {
                numeric_zero_exact(e)
            } else {
                None
            }
        }
    }
}

fn numeric_zero_exact(e: &Expr) -> Option<bool> {
    if !crate::functions::is_numeric_value(e) {
        return None;
    }
    // A machine residual is not a proof that an exact pivot is nonzero.
    // Scaling a trig identity by a large integer amplifies rounding error
    // beyond any fixed threshold. Try the existing exact pivot simplifier
    // and retain None whenever no proof decides the query.
    if let Some(answer) = crate::assumptions::is_zero_exact(e) {
        return Some(answer);
    }
    if e.float_info().max_precision.is_some() { return None; }
    let simplified = super::linalg::pivot_simp(e);
    if simplified.is_zero() { return Some(true); }
    crate::assumptions::is_zero_exact(&simplified)
}

/// SymPy's `_is_zero_after_expand_mul`.
pub fn is_zero_expand(e: &Expr) -> bool {
    if e.is_zero() {
        return true;
    }
    if e.is_atom() {
        return false;
    }
    expand(e).is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(s: &str) -> Expr {
        parse_str(s).unwrap()
    }

    #[test]
    fn cancel_basic() {
        assert_eq!(cancel(&p("(x**2-1)/(x-1)")).to_string(), "x + 1");
        assert_eq!(cancel(&p("(2*x+2)/(4*x)")).to_string(), "(x + 1)/(2*x)");
        assert_eq!(cancel(&p("x/(1 - x)")).to_string(), "-x/(x - 1)");
        assert_eq!(cancel(&p("(x/2+1)/x")).to_string(), "(x + 2)/(2*x)");
        assert_eq!(cancel(&p("1/(y - x)")).to_string(), "-1/(x - y)");
    }

    #[test]
    fn gens_order() {
        let mut g = vec![p("a"), p("y"), p("x"), p("p"), p("x10"), p("x2"), p("sin(x)")];
        g.sort_by(cmp_gens);
        let s: Vec<String> = g.iter().map(|e| e.to_string()).collect();
        assert_eq!(s, ["x", "x2", "x10", "y", "p", "a", "sin(x)"]);
    }

    #[test]
    fn ops_count() {
        assert_eq!(count_ops_alg(&p("1 + 1/x")).0, 2);
        assert!(count_ops_alg(&p("1 + 1/x")).1);
    }

    #[test]
    fn scaled_exact_trig_identity_is_a_zero_pivot() {
        for scale in [0,20,100,500] {
            let value=p(&format!("10**{scale}*(sin(3)**2+cos(3)**2-1)"));
            assert_eq!(is_zero_assume(&value),Some(true));
            let matrix=crate::matrix::Matrix::new(1,1,vec![value]).unwrap();
            assert_eq!(matrix.rank(),0);
            assert!(matrix.inv().is_err());
            assert!(matrix.inverse_ge().is_err());
            assert!(matrix.inverse_lu().is_err());
        }
    }
}
