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
use crate::ops::{add, could_extract_minus_sign, div, mul, mul_factor, neg, pow};
use crate::printer::to_str;
use crate::simplify::numer_denom;

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

pub fn cancel_opts(e: &Expr, opaque_funcs: bool) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => {}
        _ => return e.clone(),
    }
    let (n, d) = numer_denom(e);
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

/// `(numerator, denominator)` of SymPy's `fraction` for a product.
fn fraction(a: &Expr) -> (Expr, Expr) {
    numer_denom(a)
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
    dotprodsimp_withsimp(e).0
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
                let has_inv = |t: &Expr| match t.kind() {
                    Kind::Mul(_, facs) => facs.iter().any(|(_, x)| x.is_minus_one()),
                    _ => false,
                };
                let a0 = crate::ops::add_term(&terms[0].0, &terms[0].1);
                let a1 = crate::ops::add_term(&terms[1].0, &terms[1].1);
                if has_inv(&a0) && has_inv(&a1) {
                    let t = crate::simplify::together(&ex);
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

// ----------------------------------------------------------------------
// zero testing
// ----------------------------------------------------------------------

/// Approximation of SymPy's `expr.is_zero` (three-valued).
pub fn is_zero_assume(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_zero()),
        Kind::Const(Constant::NaN) => None,
        Kind::Const(_) => Some(false),
        Kind::Sym(_) => None,
        Kind::Bool(_) => None,
        Kind::Func(Func::Exp, _) => Some(false),
        Kind::Pow(b, x) => {
            if is_zero_assume(b) == Some(false) {
                if x.is_number() || b.is_number() || b.is_constant(Constant::E) {
                    return Some(false);
                }
            }
            if e.mask() == 0 {
                return numeric_nonzero(e);
            }
            None
        }
        Kind::Mul(c, facs) => {
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
                numeric_nonzero(e)
            } else {
                None
            }
        }
    }
}

fn numeric_nonzero(e: &Expr) -> Option<bool> {
    if !crate::functions::is_numeric_value(e) {
        return None;
    }
    let z = crate::evalf::eval_c64(e)?;
    if !z.re.is_finite() || !z.im.is_finite() {
        return None;
    }
    if z.norm() > 1e-12 {
        Some(false)
    } else {
        None
    }
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
}
