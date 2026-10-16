//! Conversion between expressions and polynomials.
//!
//! Generators are detected like SymPy's `Poly(expr)`: after expansion,
//! every factor of every term that is not a rational number is decomposed
//! as `base**k` (`x**(3/2)` -> `(sqrt(x), 3)`, `exp(2*x)` -> `(exp(x), 2)`,
//! `1/x` -> `(1/x, 1)`) and the bases become the generators, sorted with
//! SymPy's generator order (`x, y, z` first, then `p..w`, then `a..o`,
//! then everything else by name).

use std::cmp::Ordering;
use std::fmt;

use num_bigint::BigInt;
use num_traits::One;
use rustc_hash::FxHashMap;

use super::int::{Int, Rat};
use super::mpoly::MPoly;
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::exp as fexp;
use crate::number::Number;
use crate::ops::{add, as_coeff_mul, mul, mul_factor, pow};
use crate::printer::to_str;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PolyError {
    /// The expression cannot be made into a polynomial (no generators...).
    Failed(String),
    /// A generator appears in a coefficient position.
    Generator(String),
    /// Operation not supported for these arguments.
    Unsupported(String),
}

impl fmt::Display for PolyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolyError::Failed(s) | PolyError::Generator(s) | PolyError::Unsupported(s) => write!(f, "{s}"),
        }
    }
}

impl std::error::Error for PolyError {}

pub type PResult<T> = Result<T, PolyError>;

// ----------------------------------------------------------------------
// generator order
// ----------------------------------------------------------------------

fn gens_order(name: &str) -> u32 {
    let mut ch = name.chars();
    let (Some(c), None) = (ch.next(), ch.next()) else {
        return 1000;
    };
    match c {
        'a'..='o' => 301 + (c as u32 - 'a' as u32),
        'p'..='w' => 216 + (c as u32 - 'p' as u32),
        'x' => 124,
        'y' => 125,
        'z' => 126,
        _ => 1000,
    }
}

/// SymPy's `_sort_gens` key: `(order, name, index)`.
pub fn gen_key(g: &Expr) -> (u32, String, u64) {
    let s = to_str(g);
    let digits = s.bytes().rev().take_while(|b| b.is_ascii_digit()).count();
    let (name, idx) = s.split_at(s.len() - digits);
    let index = if idx.is_empty() { 0 } else { idx.parse::<u64>().unwrap_or(u64::MAX) };
    (gens_order(name), name.to_string(), index)
}

pub fn sort_gens(gens: &mut [Expr]) {
    gens.sort_by_cached_key(gen_key);
}

// ----------------------------------------------------------------------
// decomposition of terms
// ----------------------------------------------------------------------

/// `(base, k)` with `factor = base**k` (SymPy's `decompose_power`), with
/// negative `k` turned into `(1/base, -k)`.
pub fn decompose_power(b: &Expr, e: &Expr) -> (Expr, i64) {
    let (base, k) = decompose_raw(b, e);
    if k < 0 {
        (pow(&base, &Expr::minus_one()), -k)
    } else {
        (base, k)
    }
}

fn decompose_raw(b: &Expr, e: &Expr) -> (Expr, i64) {
    match e.kind() {
        Kind::Num(Number::Int(k)) => {
            // a power of exp(...) is canonicalised into exp(k*...)
            (b.clone(), *k)
        }
        Kind::Num(Number::Rational(r)) => {
            let q = r.denom().clone();
            let p: Option<i64> = num_traits::ToPrimitive::to_i64(r.numer());
            match (p, num_traits::ToPrimitive::to_i64(&q)) {
                (Some(p), Some(q)) => (pow(b, &Expr::rational(1, q)), p),
                _ => (mul_factor(b, e), 1),
            }
        }
        Kind::Num(_) => (mul_factor(b, e), 1),
        _ => {
            let (c, tail) = as_coeff_mul(e);
            if !c.is_rational() || c.is_one() {
                return (mul_factor(b, e), 1);
            }
            if c.is_minus_one() {
                return (pow(b, &tail), -1);
            }
            let (p, q) = c.numer_denom().unwrap();
            match (num_traits::ToPrimitive::to_i64(&p), num_traits::ToPrimitive::to_i64(&q)) {
                (Some(p), Some(q)) => {
                    let t = mul(&[Expr::rational(1, q), tail]);
                    (pow(b, &t), p)
                }
                _ => (mul_factor(b, e), 1),
            }
        }
    }
}

/// Decomposition of a single (non numeric) factor `f` as `base**k`, `k`
/// possibly negative.
pub fn decompose_factor_signed(f: &Expr) -> (Expr, i64) {
    match f.kind() {
        Kind::Pow(b, e) => decompose_raw(b, e),
        Kind::Func(Func::Exp, args) => {
            let a = &args[0];
            if a.is_number() {
                if let Some(k) = a.as_i64() {
                    return (Expr::e(), k);
                }
                return (f.clone(), 1);
            }
            let (c, tail) = as_coeff_mul(a);
            if !c.is_rational() || c.is_one() {
                return (f.clone(), 1);
            }
            let (p, q) = c.numer_denom().unwrap();
            match (num_traits::ToPrimitive::to_i64(&p), num_traits::ToPrimitive::to_i64(&q)) {
                (Some(p), Some(q)) => (fexp(&mul(&[Expr::rational(1, q), tail])), p),
                _ => (f.clone(), 1),
            }
        }
        _ => (f.clone(), 1),
    }
}

/// Decomposition of a single (non numeric) factor `f` as `base**k` with
/// `k > 0` (negative powers give the base `1/base`).
pub fn decompose_factor(f: &Expr) -> (Expr, i64) {
    let (b, k) = decompose_factor_signed(f);
    if k < 0 {
        (pow(&b, &Expr::minus_one()), -k)
    } else {
        (b, k)
    }
}

/// A term: numeric coefficient and multiplicative factors.
pub struct RawTerm {
    pub coeff: Number,
    pub factors: Vec<Expr>,
}

/// Terms of an (already expanded) expression.
pub fn raw_terms(e: &Expr) -> Vec<RawTerm> {
    fn split(t: &Expr, k: &Number) -> RawTerm {
        match t.kind() {
            Kind::Num(n) => RawTerm { coeff: k.mul(n), factors: vec![] },
            Kind::Mul(c, facs) => RawTerm {
                coeff: k.mul(c),
                factors: facs.iter().map(|(b, x)| mul_factor(b, x)).collect(),
            },
            _ => RawTerm { coeff: k.clone(), factors: vec![t.clone()] },
        }
    }
    match e.kind() {
        Kind::Add(c, terms) => {
            let mut v: Vec<RawTerm> = terms.iter().map(|(t, k)| split(t, k)).collect();
            if !c.is_zero() {
                v.push(RawTerm { coeff: c.clone(), factors: vec![] });
            }
            v
        }
        _ => vec![split(e, &Number::ONE)],
    }
}

/// Exact rational value of a float.
pub fn float_to_rat(f: f64) -> Rat {
    if f == 0.0 || !f.is_finite() {
        return Rat::zero();
    }
    let bits = f.to_bits();
    let sign = if (bits >> 63) != 0 { -1 } else { 1 };
    let exp = ((bits >> 52) & 0x7ff) as i64;
    let mant = bits & ((1u64 << 52) - 1);
    let (m, e) = if exp == 0 { (mant, -1074) } else { (mant | (1u64 << 52), exp - 1075) };
    let m = BigInt::from(m) * sign;
    if e >= 0 {
        Rat::from(Int::from_big(m << (e as usize)))
    } else {
        Rat::new(Int::from_big(m), Int::from_big(BigInt::one() << ((-e) as usize)))
    }
}

pub fn number_to_rat(n: &Number) -> Rat {
    match n {
        Number::Float(f) => float_to_rat(*f),
        Number::MpFloat(value, _) => {
            let exact = value.to_ratio_exact();
            Rat::new(Int::from_big(exact.numer().clone()), Int::from_big(exact.denom().clone()))
        }
        _ => Rat::from_number(n).unwrap(),
    }
}

fn check_number(e: &Expr) -> PResult<()> {
    match e.kind() {
        Kind::Const(c) if c.is_infinite() || *c == Constant::NaN => {
            Err(PolyError::Failed(format!("cannot build a polynomial from {}", to_str(e))))
        }
        Kind::Rel(..) | Kind::Bool(_) | Kind::Piecewise(_) => {
            Err(PolyError::Failed(format!("cannot build a polynomial from {}", to_str(e))))
        }
        _ => Ok(()),
    }
}

/// Result of a conversion: polynomials over a common ring.
#[derive(Clone, Debug)]
pub struct Converted {
    pub polys: Vec<MPoly<Rat>>,
    /// All variables: the generators first, then the coefficient
    /// ("domain") generators.
    pub gens: Vec<Expr>,
    /// Number of true generators.
    pub ngens: usize,
    /// Some coefficient was a float.
    pub float: bool,
}

fn factor_needs_expand(b: &Expr, x: &Expr) -> bool {
    if b.is_add() {
        return true;
    }
    if x.is_add() {
        return true;
    }
    match b.kind() {
        Kind::Func(Func::Exp, a) => a[0].is_add(),
        Kind::Mul(..) => true,
        _ => false,
    }
}

/// Is `expand(e)` possibly different from `e` (as far as the polynomial
/// structure is concerned)?
pub fn needs_expand(e: &Expr) -> bool {
    let term = |t: &Expr| -> bool {
        match t.kind() {
            Kind::Mul(_, facs) => facs.iter().any(|(b, x)| factor_needs_expand(b, x)),
            Kind::Pow(b, x) => factor_needs_expand(b, x),
            Kind::Func(Func::Exp, a) => a[0].is_add(),
            _ => false,
        }
    };
    match e.kind() {
        Kind::Add(_, terms) => terms.iter().any(|(t, _)| term(t)),
        _ => term(e),
    }
}

fn expand_if_needed(e: &Expr) -> Expr {
    if needs_expand(e) {
        expand(e)
    } else {
        e.clone()
    }
}

/// Convert expressions with automatic generator detection.
pub fn from_exprs(exprs: &[Expr]) -> PResult<Converted> {
    let mut float = false;
    let mut index: FxHashMap<Expr, usize> = FxHashMap::default();
    let mut gens: Vec<Expr> = Vec::new();
    let mut all_terms: Vec<Vec<(Rat, Vec<(usize, i64)>)>> = Vec::new();
    for e in exprs {
        check_number(e)?;
        let ex = expand_if_needed(e);
        let mut terms = Vec::new();
        for t in raw_terms(&ex) {
            if t.coeff.is_float() {
                float = true;
            }
            let mut mono = Vec::new();
            for f in &t.factors {
                check_number(f)?;
                let (b, k) = decompose_factor(f);
                let id = *index.entry(b.clone()).or_insert_with(|| {
                    gens.push(b.clone());
                    gens.len() - 1
                });
                mono.push((id, k));
            }
            terms.push((number_to_rat(&t.coeff), mono));
        }
        all_terms.push(terms);
    }
    let mut order: Vec<usize> = (0..gens.len()).collect();
    let keys: Vec<_> = gens.iter().map(gen_key).collect();
    order.sort_by(|&a, &b| keys[a].cmp(&keys[b]).then_with(|| gens[a].cmp_struct(&gens[b])));
    let mut pos = vec![0usize; gens.len()];
    for (i, &g) in order.iter().enumerate() {
        pos[g] = i;
    }
    let n = gens.len();
    let sorted: Vec<Expr> = order.iter().map(|&i| gens[i].clone()).collect();
    let polys = all_terms
        .into_iter()
        .map(|terms| {
            let ts = terms
                .into_iter()
                .map(|(c, mono)| {
                    let mut e = vec![0u32; n];
                    for (g, k) in mono {
                        e[pos[g]] += k as u32;
                    }
                    (e, c)
                })
                .collect();
            MPoly::from_terms(n, ts)
        })
        .collect();
    Ok(Converted { polys, gens: sorted, ngens: n, float })
}

/// Convert expressions with respect to the given generators. Factors free
/// of the generators go to the coefficients, which are themselves
/// polynomials in automatically detected coefficient generators.
pub fn from_exprs_gens(exprs: &[Expr], gens: &[Expr]) -> PResult<Converted> {
    let k = gens.len();
    let mut index: FxHashMap<Expr, usize> = FxHashMap::default();
    for (i, g) in gens.iter().enumerate() {
        index.insert(g.clone(), i);
    }
    let mut float = false;
    // per expression: list of (gen monomial, coefficient expression factors, numeric coeff)
    let mut dom_index: FxHashMap<Expr, usize> = FxHashMap::default();
    let mut dom_gens: Vec<Expr> = Vec::new();
    #[allow(clippy::type_complexity)]
    let mut all: Vec<Vec<(Vec<u32>, Rat, Vec<(usize, i64)>)>> = Vec::new();
    for e in exprs {
        check_number(e)?;
        let ex = expand_if_needed(e);
        let mut terms = Vec::new();
        for t in raw_terms(&ex) {
            if t.coeff.is_float() {
                float = true;
            }
            let mut mono = vec![0u32; k];
            let mut dmono = Vec::new();
            for f in &t.factors {
                check_number(f)?;
                let (b, kk) = decompose_factor(f);
                if let Some(&i) = index.get(&b) {
                    mono[i] += kk as u32;
                    continue;
                }
                if gens.iter().any(|g| crate::traverse::has(f, g)) {
                    return Err(PolyError::Generator(format!(
                        "{} contains an element of the set of generators.",
                        to_str(f)
                    )));
                }
                let id = *dom_index.entry(b.clone()).or_insert_with(|| {
                    dom_gens.push(b.clone());
                    dom_gens.len() - 1
                });
                dmono.push((id, kk));
            }
            terms.push((mono, number_to_rat(&t.coeff), dmono));
        }
        all.push(terms);
    }
    let mut order: Vec<usize> = (0..dom_gens.len()).collect();
    let keys: Vec<_> = dom_gens.iter().map(gen_key).collect();
    order.sort_by(|&a, &b| keys[a].cmp(&keys[b]).then_with(|| dom_gens[a].cmp_struct(&dom_gens[b])));
    let mut pos = vec![0usize; dom_gens.len()];
    for (i, &g) in order.iter().enumerate() {
        pos[g] = i;
    }
    let n = k + dom_gens.len();
    let mut allgens: Vec<Expr> = gens.to_vec();
    allgens.extend(order.iter().map(|&i| dom_gens[i].clone()));
    let polys = all
        .into_iter()
        .map(|terms| {
            let ts = terms
                .into_iter()
                .map(|(mono, c, dm)| {
                    let mut e = mono;
                    e.resize(n, 0);
                    for (g, kk) in dm {
                        e[k + pos[g]] += kk as u32;
                    }
                    (e, c)
                })
                .collect();
            MPoly::from_terms(n, ts)
        })
        .collect();
    Ok(Converted { polys, gens: allgens, ngens: k, float })
}

pub fn rat_to_expr(c: &Rat) -> Expr {
    Expr::num(c.to_number())
}

pub fn rat_to_expr_float(c: &Rat, float: bool) -> Expr {
    if float {
        Expr::float(c.to_f64())
    } else {
        rat_to_expr(c)
    }
}

/// Expression of a polynomial in the given generators.
pub fn to_expr(p: &MPoly<Rat>, gens: &[Expr]) -> Expr {
    to_expr_f(p, gens, false)
}

pub fn to_expr_f(p: &MPoly<Rat>, gens: &[Expr], float: bool) -> Expr {
    let mut terms = Vec::with_capacity(p.len());
    for (e, c) in p.iter() {
        let mut f = Vec::with_capacity(gens.len() + 1);
        f.push(rat_to_expr_float(c, float));
        for (v, &k) in e.iter().enumerate() {
            if k != 0 {
                f.push(pow(&gens[v], &Expr::int(k as i64)));
            }
        }
        terms.push(mul(&f));
    }
    add(&terms)
}

pub fn to_expr_int(p: &MPoly<Int>, gens: &[Expr]) -> Expr {
    to_expr(&p.to_rat(), gens)
}

/// Re-embed polynomials over `from` generators into `to` generators
/// (every generator of `from` must be in `to`).
pub fn reembed(p: &MPoly<Rat>, from: &[Expr], to: &[Expr]) -> MPoly<Rat> {
    let perm: Vec<usize> = from.iter().map(|g| to.iter().position(|h| h == g).expect("generator")).collect();
    p.permute(&perm, to.len())
}

/// Compare two generator lists (for unification) with SymPy's order.
pub fn cmp_gens(a: &Expr, b: &Expr) -> Ordering {
    gen_key(a).cmp(&gen_key(b)).then_with(|| a.cmp_struct(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn gens_of(s: &str) -> Vec<String> {
        let c = from_exprs(&[parse_str(s).unwrap()]).unwrap();
        c.gens.iter().map(to_str).collect()
    }

    #[test]
    fn generators() {
        assert_eq!(gens_of("x + y + z + t + a + b + c"), vec!["x", "y", "z", "t", "a", "b", "c"]);
        assert_eq!(gens_of("pi*x + sin(y)"), vec!["x", "pi", "sin(y)"]);
        assert_eq!(gens_of("1/x + x"), vec!["x", "1/x"]);
        assert_eq!(gens_of("x**(3/2) + x"), vec!["x", "sqrt(x)"]);
        assert_eq!(gens_of("exp(2*x) + 1"), vec!["exp(x)"]);
        assert_eq!(gens_of("x1 + x10 + x2"), vec!["x1", "x2", "x10"]);
        let c = from_exprs(&[parse_str("(x + 1)**2").unwrap()]).unwrap();
        assert_eq!(to_str(&to_expr(&c.polys[0], &c.gens)), "x**2 + 2*x + 1");
        let c = from_exprs_gens(&[parse_str("x*y + 1").unwrap()], &[Expr::symbol("x")]).unwrap();
        assert_eq!(c.ngens, 1);
        assert_eq!(c.gens.len(), 2);
        assert!(from_exprs_gens(&[parse_str("sin(x)*x").unwrap()], &[Expr::symbol("x")]).is_err());
        assert_eq!(float_to_rat(0.5), Rat::new(Int::S(1), Int::S(2)));
    }
}
