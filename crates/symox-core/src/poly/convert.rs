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
use num_traits::{One, ToPrimitive, Zero};
use rustc_hash::{FxHashMap, FxHashSet};

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
    /// A polynomial or residue-class element has no multiplicative inverse.
    NotInvertible(String),
    /// Division does not stay in the selected coefficient ring.
    ExactQuotientFailed(String),
    /// A reference exception raised by the dense `EX` algorithms, named by
    /// its class (`PolynomialDivisionFailed`, `ZeroDivisionError`, ...).
    Reference { kind: &'static str, message: String },
}

impl fmt::Display for PolyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolyError::Failed(s) | PolyError::Generator(s) | PolyError::Unsupported(s) | PolyError::NotInvertible(s) | PolyError::ExactQuotientFailed(s) => write!(f, "{s}"),
            PolyError::Reference { message, .. } => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for PolyError {}

pub type PResult<T> = Result<T, PolyError>;

/// Monomial storage is u32. A representation limit must never truncate an
/// exponent or wrap a sum of exponents into a different polynomial.
pub(crate) fn add_exponent(slot: &mut u32, exponent: i64) -> PResult<()> {
    let limit = || PolyError::Unsupported("Polynomial exponent storage limit exceeded (32-bit monomials)".into());
    let exponent = u32::try_from(exponent).map_err(|_| limit())?;
    *slot = slot.checked_add(exponent).ok_or_else(limit)?;
    Ok(())
}

/// Keep the ordinary u32 path cheap, but combine exact coefficients before
/// enforcing the storage bound when a held source contains oversized terms.
fn source_polynomial(nvars: usize, terms: Vec<(Rat, Vec<(usize, i64)>)>) -> PResult<MPoly<Rat>> {
    let narrow = terms.iter().map(|(c, powers)| {
        let mut exponents = vec![0u32; nvars];
        if !c.is_zero() {
            for &(i, k) in powers { add_exponent(&mut exponents[i], k)?; }
        }
        Ok((exponents, c.clone()))
    }).collect::<PResult<Vec<_>>>();
    if let Ok(narrow) = narrow { return Ok(MPoly::from_terms(nvars, narrow)); }
    let mut wide: FxHashMap<Vec<BigInt>, Rat> = FxHashMap::default();
    for (c, powers) in terms {
        if c.is_zero() { continue; }
        let mut exponents = vec![BigInt::zero(); nvars];
        for (i, k) in powers { exponents[i] += k; }
        wide.entry(exponents).and_modify(|a| *a = a.add(&c)).or_insert(c);
    }
    let narrow = wide.into_iter().filter(|(_, c)| !c.is_zero()).map(|(exponents, c)| {
        let exponents = exponents.into_iter().map(|e| e.to_u32().ok_or_else(||
            PolyError::Unsupported("Polynomial exponent storage limit exceeded (32-bit monomials)".into())
        )).collect::<PResult<Vec<_>>>()?;
        Ok((exponents, c))
    }).collect::<PResult<Vec<_>>>()?;
    Ok(MPoly::from_terms(nvars, narrow))
}

// ----------------------------------------------------------------------
// generator order
// ----------------------------------------------------------------------

/// Conventional rank of a generator's base name (before a numeric suffix).
pub fn gens_order(name: &str) -> u32 {
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

/// The generators of `_parallel_dict_from_expr_no_gens` with
/// `greedy=False`: after expansion every multiplicative factor that is a
/// number (`sqrt(2)`, `pi`, `CRootOf(...)`) stays in the coefficients and
/// every other factor contributes its `decompose_power` base, so
/// `x - sqrt(2)*t` has the generators `t, x`.
pub fn nongreedy_generators(expressions: &[Expr]) -> Vec<Expr> {
    let mut gens: Vec<Expr> = Vec::new();
    for e in expressions {
        let e = expand(e);
        let terms = if e.is_add() { e.args() } else { vec![e] };
        for t in terms {
            let factors = if t.is_mul() { t.args() } else { vec![t] };
            for f in factors {
                if crate::functions::is_numeric_value(&f) {
                    continue;
                }
                let base = decompose_factor(&f).0;
                if !gens.contains(&base) {
                    gens.push(base);
                }
            }
        }
    }
    sort_gens(&mut gens);
    gens
}

// ----------------------------------------------------------------------
// decomposition of terms
// ----------------------------------------------------------------------

/// `(base, k)` with `factor = base**k` (SymPy's `decompose_power`), with
/// negative `k` turned into `(1/base, -k)`.
pub fn decompose_power(b: &Expr, e: &Expr) -> (Expr, i64) {
    let (base, k) = decompose_raw(b, e);
    if k < 0 {
        match k.checked_neg() {
            Some(k) => (pow(&base, &Expr::minus_one()), k),
            None => (mul_factor(b,e),1),
        }
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
        match k.checked_neg() {
            Some(k) => (pow(&b, &Expr::minus_one()), k),
            None => (f.clone(),1),
        }
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
            // As in the reference `_dict_from_expr`, a number factor of an
            // unevaluated product is a coefficient: Mul(-1, 2, x,
            // evaluate=False) is -2*x, not a polynomial in the generator 2.
            Kind::Mul(c, facs) => {
                let mut coeff = k.mul(c);
                let mut factors = Vec::with_capacity(facs.len());
                for (b, x) in facs {
                    match b.kind() {
                        Kind::Num(n) if x.is_one() && !n.is_nan() => coeff = coeff.mul(n),
                        _ => factors.push(mul_factor(b, x)),
                    }
                }
                RawTerm { coeff, factors }
            }
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

/// SymPy's `RR.to_rational(a)` (53-bit context, `limit=True`): the exact
/// value when its denominator is at most `int(1/tol)` with
/// `tol = 100*2**-52`, otherwise the closest fraction with a bounded
/// denominator (mpmath's continued-fraction `limit_denominator`).
pub fn rr_limit_rational(r: &Rat) -> Rat {
    let max_denom = BigInt::from(45_035_996_273_704u64);
    let (p, q) = (r.n.to_big(), r.d.to_big());
    if q <= max_denom {
        return r.clone();
    }
    let (mut p0, mut q0, mut p1, mut q1) = (BigInt::from(0), BigInt::from(1), BigInt::from(1), BigInt::from(0));
    let (mut n, mut d) = (p.clone(), q.clone());
    loop {
        let a = num_integer::Integer::div_floor(&n, &d);
        let q2 = &q0 + &a * &q1;
        if q2 > max_denom {
            break;
        }
        let p2 = &p0 + &a * &p1;
        p0 = std::mem::replace(&mut p1, p2);
        q0 = std::mem::replace(&mut q1, q2);
        let rem = &n - &a * &d;
        n = std::mem::replace(&mut d, rem);
    }
    let k = num_integer::Integer::div_floor(&(&max_denom - &q0), &q1);
    let rat = |a: BigInt, b: BigInt| Rat::new(Int::from_big(a), Int::from_big(b));
    let bound1 = rat(&p0 + &k * &p1, &q0 + &k * &q1);
    let bound2 = rat(p1, q1);
    if bound2.is_zero() || bound1.is_zero() {
        return r.clone();
    }
    if bound2.sub(r).abs().sub(&bound1.sub(r).abs()).is_negative() || bound2.sub(r).abs() == bound1.sub(r).abs() {
        bound2
    } else {
        bound1
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

/// SymPy's `_dict_from_expr` rejects noncommutative expressions before any
/// conversion: polynomial arithmetic would commute their factors.
pub(crate) fn check_commutative(e: &Expr) -> PResult<()> {
    if crate::decompose::commutes(e) {
        Ok(())
    } else {
        Err(PolyError::Reference {
            kind: "PolynomialError",
            message: "non-commutative expressions are not supported".into(),
        })
    }
}

pub(crate) fn check_number(e: &Expr) -> PResult<()> {
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
    /// The largest binary precision of the float coefficients (the
    /// reference `RealField(prec)` of `construct_domain`), 0 without floats.
    pub prec: u32,
}

fn factor_needs_expand(b: &Expr, x: &Expr, seen: &mut FxHashSet<usize>) -> bool {
    if b.is_add() {
        return true;
    }
    if x.is_add() {
        return true;
    }
    match b.kind() {
        Kind::Mul(..) => true,
        _ => function_needs_expand(b, seen) || needs_expand_in(x, seen),
    }
}

/// The reference expands deeply before building a polynomial: the
/// arguments of a function are expanded, and `exp` of a sum splits, so
/// `exp(x*(x + 1))` gives the generators `exp(x)` and `exp(x**2)`. `seen`
/// holds the function nodes already found not to need it (shared
/// subexpressions are visited once).
fn function_needs_expand(e: &Expr, seen: &mut FxHashSet<usize>) -> bool {
    let Kind::Func(head, args) = e.kind() else { return false };
    if !seen.insert(e.dag_id()) {
        return false;
    }
    match head {
        Func::Exp => args[0].is_add() || needs_expand_in(&args[0], seen),
        _ => args.iter().any(|a| needs_expand_in(a, seen)),
    }
}

/// Is `expand(e)` possibly different from `e` (as far as the polynomial
/// structure is concerned)?
pub fn needs_expand(e: &Expr) -> bool {
    needs_expand_in(e, &mut FxHashSet::default())
}

fn needs_expand_in(e: &Expr, seen: &mut FxHashSet<usize>) -> bool {
    let mut term = |t: &Expr| -> bool {
        match t.kind() {
            Kind::Mul(_, facs) => facs.iter().any(|(b, x)| factor_needs_expand(b, x, seen)),
            Kind::Pow(b, x) => factor_needs_expand(b, x, seen),
            Kind::Func(..) => function_needs_expand(t, seen),
            _ => false,
        }
    };
    match e.kind() {
        Kind::Add(_, terms) => terms.iter().any(|(t, _)| term(t)),
        _ => term(e),
    }
}

pub(crate) fn expand_if_needed(e: &Expr) -> Expr {
    if needs_expand(e) {
        expand(e)
    } else {
        e.clone()
    }
}

/// The factor of a term that the reference names when it contains a
/// generator without being a power of one: the first such factor in the
/// order of the arguments of the term.
pub(crate) fn generator_clash<'a>(factors: &'a [Expr], gens: &[Expr]) -> Option<&'a Expr> {
    factors
        .iter()
        .filter(|f| !gens.contains(&decompose_factor(f).0) && gens.iter().any(|g| crate::traverse::has_structural(f, g)))
        .min_by(|a, b| crate::traverse::args_cmp(a, b))
}

/// Convert expressions with automatic generator detection.
pub fn from_exprs(exprs: &[Expr]) -> PResult<Converted> {
    let mut float = false;
    let mut prec = 0u32;
    let mut index: FxHashMap<Expr, usize> = FxHashMap::default();
    let mut gens: Vec<Expr> = Vec::new();
    let mut all_terms: Vec<Vec<(Rat, Vec<(usize, i64)>)>> = Vec::new();
    for e in exprs {
        check_commutative(e)?;
        check_number(e)?;
        let ex = expand_if_needed(e);
        let mut terms = Vec::new();
        for t in raw_terms(&ex) {
            if t.coeff.is_float() {
                float = true;
                prec = prec.max(t.coeff.float_precision().unwrap_or(53));
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
    let polys = all_terms.into_iter().map(|terms| {
        source_polynomial(n, terms.into_iter().map(|(c, mono)|
            (c, mono.into_iter().map(|(g, k)| (pos[g], k)).collect())
        ).collect())
    }).collect::<PResult<Vec<_>>>()?;
    Ok(Converted { polys, gens: sorted, ngens: n, float, prec })
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
    let mut prec = 0u32;
    // per expression: list of (gen monomial, coefficient expression factors, numeric coeff)
    let mut dom_index: FxHashMap<Expr, usize> = FxHashMap::default();
    let mut dom_gens: Vec<Expr> = Vec::new();
    #[allow(clippy::type_complexity)]
    let mut all: Vec<Vec<(Vec<(usize, i64)>, Rat, Vec<(usize, i64)>)>> = Vec::new();
    for e in exprs {
        check_commutative(e)?;
        check_number(e)?;
        let ex = expand_if_needed(e);
        let mut terms = Vec::new();
        for t in raw_terms(&ex) {
            if t.coeff.is_float() {
                float = true;
                prec = prec.max(t.coeff.float_precision().unwrap_or(53));
            }
            let mut mono = Vec::new();
            let mut dmono = Vec::new();
            for f in &t.factors {
                check_number(f)?;
                let (b, kk) = decompose_factor(f);
                if let Some(&i) = index.get(&b) {
                    mono.push((i, kk));
                    continue;
                }
                if gens.iter().any(|g| crate::traverse::has_structural(f, g)) {
                    return Err(PolyError::Generator(format!(
                        "{} contains an element of the set of generators.",
                        to_str(generator_clash(&t.factors, gens).unwrap_or(f))
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
    let polys = all.into_iter().map(|terms| {
        source_polynomial(n, terms.into_iter().map(|(mut mono, c, dm)| {
            mono.extend(dm.into_iter().map(|(g, kk)| (k + pos[g], kk)));
            (c, mono)
        }).collect())
    }).collect::<PResult<Vec<_>>>()?;
    Ok(Converted { polys, gens: allgens, ngens: k, float, prec })
}

pub fn rat_to_expr(c: &Rat) -> Expr {
    Expr::num(c.to_number())
}

pub fn rat_to_expr_float(c: &Rat, float: bool) -> Expr {
    rat_to_expr_prec(c, if float { 53 } else { 0 })
}

/// A coefficient as a `Float` of `prec` bits (the reference
/// `RealField(prec).to_sympy`), or exactly when `prec` is 0.
pub fn rat_to_expr_prec(c: &Rat, prec: u32) -> Expr {
    if prec > 0 && prec != 53 {
        let value = crate::bigfloat::BigFloat::from_ratio(&c.n.to_big(), &c.d.to_big(), prec);
        return Expr::num(crate::number::Number::mp_float(value, prec));
    }
    if prec > 0 {
        let f = c.to_f64();
        if f.is_finite() && (f != 0.0 || c.is_zero()) {
            return Expr::float(f);
        }
        // outside the double range: a 53-bit float with a wider exponent,
        // as the reference RR domain keeps it
        let value = crate::bigfloat::BigFloat::from_ratio(&c.n.to_big(), &c.d.to_big(), 53);
        Expr::num(crate::number::Number::mp_float(value, 53))
    } else {
        rat_to_expr(c)
    }
}

/// Expression of a polynomial in the given generators.
pub fn to_expr(p: &MPoly<Rat>, gens: &[Expr]) -> Expr {
    to_expr_f(p, gens, false)
}

pub fn to_expr_f(p: &MPoly<Rat>, gens: &[Expr], float: bool) -> Expr {
    to_expr_prec(p, gens, if float { 53 } else { 0 })
}

/// Expression of a polynomial whose coefficients are `Float`s of `prec`
/// bits (exact rationals when `prec` is 0).
pub fn to_expr_prec(p: &MPoly<Rat>, gens: &[Expr], prec: u32) -> Expr {
    let mut terms = Vec::with_capacity(p.len());
    for (e, c) in p.iter() {
        let mut f = Vec::with_capacity(gens.len() + 1);
        f.push(rat_to_expr_prec(c, prec));
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

    #[test]
    fn exponent_storage_never_truncates_or_wraps() {
        let x=Expr::symbol("x");let y=Expr::symbol("y");
        let too_large=pow(&x,&Expr::int(1i64<<32));
        assert!(matches!(from_exprs(&[too_large.clone()]),Err(PolyError::Unsupported(_))));
        assert!(matches!(from_exprs_gens(&[too_large.clone()],std::slice::from_ref(&x)),Err(PolyError::Unsupported(_))));
        assert!(matches!(super::super::expression_domain::dict_from_expr(&too_large,std::slice::from_ref(&x)),Err(PolyError::Unsupported(_))));
        let legal=pow(&x,&Expr::int(u32::MAX as i64));
        let converted=from_exprs_gens(&[legal],std::slice::from_ref(&x)).unwrap();
        assert_eq!(converted.polys[0].iter().next().unwrap().0[0],u32::MAX);
        let held=Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,vec![(x.clone(),Expr::int(u32::MAX as i64)),(x.clone(),Expr::one())]));
        assert!(matches!(from_exprs_gens(&[held.clone()],std::slice::from_ref(&x)),Err(PolyError::Unsupported(_))));
        assert!(matches!(super::super::expression_domain::dict_from_expr(&held,std::slice::from_ref(&x)),Err(PolyError::Unsupported(_))));
        let coefficient=mul(&[x.clone(),pow(&y,&Expr::int(1i64<<32))]);
        assert!(matches!(from_exprs_gens(&[coefficient],std::slice::from_ref(&x)),Err(PolyError::Unsupported(_))));
    }

    #[test]
    fn held_generator_factors_accumulate_their_whole_degree() {
        let x=Expr::symbol("x");
        let held=Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,vec![(x.clone(),Expr::one()),(x.clone(),Expr::int(2))]));
        let dictionary=super::super::expression_domain::dict_from_expr(&held,std::slice::from_ref(&x)).unwrap();
        assert_eq!(dictionary,vec![(vec![3],Expr::one())]);
        let extreme=pow(&x,&Expr::int(i64::MIN));
        assert_eq!(decompose_factor(&extreme),(extreme.clone(),1));
    }

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
