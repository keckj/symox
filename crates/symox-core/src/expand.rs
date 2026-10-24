//! Expansion of products and powers of sums.
//!
//! The heavy lifting is done on a sparse polynomial representation over
//! "generators" (symbols, functions, radicals, ...) with exponent vectors
//! stored inline, which makes `expand((x+y+z+1)**20)` orders of magnitude
//! faster than a naive term-by-term expression multiplication.

use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use crate::expr::{Expr, Func, Kind};
use crate::functions::exp;
use crate::number::Number;
use crate::ops::{add, add_from_terms, mul, mul_factor, pow};

mod options;
mod memo;
pub use options::{expand_with_options, ExpandBudget, ExpandError, ExpandOptions};

type Mono = SmallVec<[i64; 8]>;

/// Default `expand()` (hints `mul`, `multinomial`, `power_exp`).
pub fn expand(e: &Expr) -> Expr {
    arithmetic_expand(e, &ExpandOptions::default())
}

fn arithmetic_expand(e: &Expr, options: &ExpandOptions) -> Expr {
    arithmetic_expand_with(e, options, &mut |_| None)
}
fn arithmetic_expand_with(
    e: &Expr,
    options: &ExpandOptions,
    opaque: &mut dyn FnMut(&Expr) -> Option<Expr>,
) -> Expr {
    let mut cache = FxHashMap::default();
    exf(e, &mut cache, options, opaque)
}

/// `ex` for a subexpression that is not flattened into its parent (a power
/// base, a function argument, the root): its terms are final, so the
/// multiplication hint now applies to each of them as a whole product.
fn exf(
    e: &Expr,
    cache: &mut FxHashMap<usize, (Expr, Expr)>,
    options: &ExpandOptions,
    opaque: &mut dyn FnMut(&Expr) -> Option<Expr>,
) -> Expr {
    let r = ex(e, cache, options, opaque);
    if options.deep && options.mul { fractional_sum_denominators(r) } else { r }
}

fn ex(
    e: &Expr,
    cache: &mut FxHashMap<usize, (Expr, Expr)>,
    options: &ExpandOptions,
    opaque: &mut dyn FnMut(&Expr) -> Option<Expr>,
) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    if let Some((_, v)) = cache.get(&e.dag_id()) {
        return v.clone();
    }
    if matches!(
        e.kind(),
        Kind::Func(
            Func::MatrixProxy(_)
                | Func::MatrixElement(_)
                | Func::MatrixTrace(_)
                | Func::MatrixDeterminant(_)
                | Func::MatrixDotProduct(_, _),
            _
        )
    ) {
        if let Some(value) = opaque(e) {
            cache.insert(e.dag_id(), (e.clone(), value.clone()));
            return value;
        }
    }
    let r = match e.kind() {
        Kind::Add(c, terms) => {
            let mut v: Vec<Expr> = Vec::with_capacity(terms.len() + 1);
            let mut unchanged: Vec<(Expr, Number)> = Vec::new();
            for (t, k) in terms {
                let nt = if options.deep {
                    ex(t, cache, options, opaque)
                } else {
                    t.clone()
                };
                // A held numeric multiple of a sum, like the reference
                // `(2 - I)*(1/5)` for `1/(2 + I)`, is a product that the
                // multiplication hint distributes.
                let held = options.deep && options.mul && nt.is_add() && !k.is_one();
                if nt.same_presentation(t) && !held {
                    unchanged.push((t.clone(), k.clone()));
                } else if held && !crate::decompose::commutes(&nt) {
                    // a number times a noncommutative sum stays a product
                    // when built, so its terms are distributed one by one
                    for a in nt.args() {
                        v.push(mul(&[Expr::num(k.clone()), a]));
                    }
                } else {
                    v.push(mul(&[Expr::num(k.clone()), nt]));
                }
            }
            if v.is_empty() {
                e.clone()
            } else {
                v.push(add_from_terms(c.clone(), unchanged));
                add(&v)
            }
        }
        Kind::Mul(c, facs) => {
            // `c*f` for a number `c` distributes over a sum `f` when it is
            // built, so the terms of `f` are not final yet
            let flattened = !c.is_one() && facs.len() == 1;
            let factors: Vec<Expr> = e
                .args()
                .iter()
                .map(|f| {
                    if options.deep && flattened {
                        ex(f, cache, options, opaque)
                    } else if options.deep {
                        exf(f, cache, options, opaque)
                    } else {
                        f.clone()
                    }
                })
                .collect();
            let m = mul(&factors);
            let mut r = if !options.mul {
                m
            } else if m.is_mul() && m.has_noncommutative_host() && factors.iter().any(|f| f.is_add()) {
                // a number times a noncommutative sum stays a product when
                // built (Mul.flatten), so it is distributed here
                expand_product(&factors)
            } else if m.is_mul() {
                expand_mul_node(&m)
            } else if m.is_add() && !factors.iter().any(|f| f.is_add()) {
                m
            } else if options.mul && factors.iter().any(|f| f.is_add()) {
                expand_product(&factors)
            } else {
                m
            };
            // rebuilding the product may merge expanded sums into new
            // powers ((x**2 + 2*x + 1)*(x + 1)**2): SymPy's fixed point loop
            // expands them again
            for _ in 0..3 {
                if !options.deep
                    || !options.mul
                    || !options.multinomial
                    || r.same_presentation(e)
                    || !has_unexpanded_power(&r)
                {
                    break;
                }
                r = ex(&r, cache, options, opaque);
            }
            r
        }
        Kind::Pow(b, x) => {
            let b2 = if options.deep {
                exf(b, cache, options, opaque)
            } else {
                b.clone()
            };
            let x2 = if options.deep {
                exf(x, cache, options, opaque)
            } else {
                x.clone()
            };
            // SymPy's multinomial hint leaves base**-1 alone, and no other
            // hint applies to it: with an unchanged sum as base the node is
            // returned as is (a held 1/(0.99 + 0.2*I) stays held).
            if x2.is_minus_one() && b2.is_add() && b2 == *b {
                e.clone()
            } else {
                expand_pow(&b2, &x2, cache, options, opaque)
            }
        }
        Kind::Func(Func::Exp, args) => {
            let a = if options.deep {
                exf(&args[0], cache, options, opaque)
            } else {
                args[0].clone()
            };
            match a.kind() {
                Kind::Add(..) if options.power_exp && crate::decompose::commutes(&a) => {
                    let parts: Vec<Expr> = a.args().iter().map(exp).collect();
                    mul(&parts)
                }
                _ => exp(&a),
            }
        }
        // AddWithLimits._eval_expand_basic: an Integral or a Sum of a sum
        // is the sum of the Integrals (Sums) of its terms, when the summand
        // commutes and the limits are not known to be infinite (or force).
        // A held log(x, b, evaluate=False) is rebuilt, which evaluates it.
        Kind::Func(Func::Undefined(head), args) if args.len() == 2 && crate::functions::is_held_log_base(head) => {
            exf(&crate::functions::log_base(&args[0], &args[1]), cache, options, opaque)
        }
        Kind::Func(Func::Integral | Func::Sum, _) if options.basic => {
            let rebuilt = if options.deep {
                crate::traverse::map_children(e, &mut |c| exf(c, cache, options, opaque))
            } else {
                e.clone()
            };
            split_add_with_limits(&rebuilt, options.force).unwrap_or(rebuilt)
        }
        _ if options.deep => {
            crate::traverse::map_children(e, &mut |c| exf(c, cache, options, opaque))
        }
        _ => e.clone(),
    };
    cache.insert(e.dag_id(), (e.clone(), r.clone()));
    r
}

/// SymPy's `ExprWithLimits.has_finite_limits`: False when a bound is
/// infinite, None when a bound may be infinite and the variable may be too.
fn has_finite_limits(limits: &[Expr]) -> Option<bool> {
    let infinite = |e: &Expr| crate::assumptions::fact(e, "infinite");
    let mut unknown = false;
    for limit in limits {
        let parts = limit.args();
        if parts.len() == 3 {
            if parts[1..].iter().any(|b| infinite(b) == Some(true)) {
                return Some(false);
            }
            if parts[1..].iter().any(|b| infinite(b).is_none()) && infinite(&parts[0]).is_none() {
                unknown = true;
            }
        } else if parts.first().is_some_and(|v| infinite(v).is_none()) {
            unknown = true;
        }
    }
    if unknown { None } else { Some(true) }
}

fn split_add_with_limits(e: &Expr, force: bool) -> Option<Expr> {
    let Kind::Func(head @ (Func::Integral | Func::Sum), args) = e.kind() else { return None };
    let summand = args.first()?;
    if !summand.is_add()
        || !(force || (crate::decompose::commutes(summand) && has_finite_limits(&args[1..]) != Some(false)))
    {
        return None;
    }
    let parts: Vec<Expr> = summand
        .args()
        .iter()
        .map(|term| {
            let mut pieces = Vec::with_capacity(args.len());
            pieces.push(term.clone());
            pieces.extend(args[1..].iter().cloned());
            crate::functions::func(head.clone(), pieces)
        })
        .collect();
    Some(add(&parts))
}

/// The reference expands every term `k*t` of a sum as a whole product: a
/// fractional coefficient joins the denominator of `t` (`fraction`), and a
/// number times a sum is distributed, so `1/(29791*(a + b))` becomes
/// `1/(29791*a + 29791*b)`.
fn fractional_sum_denominators(r: Expr) -> Expr {
    if let Kind::Mul(k, _) = r.kind() {
        return if k.is_rational() && !k.is_integer() && sum_denominator(&r) && has_multi_denominator(&r) {
            expand_mul_node(&r)
        } else {
            r
        };
    }
    let Kind::Add(_, terms) = r.kind() else { return r };
    if !terms.iter().any(|(t, k)| k.is_rational() && !k.is_integer() && sum_denominator(t)) {
        return r;
    }
    let parts: Vec<Expr> = r
        .args()
        .iter()
        .map(|t| if t.is_mul() && has_multi_denominator(t) { expand_mul_node(t) } else { t.clone() })
        .collect();
    add(&parts)
}

fn expand_pow(
    b: &Expr,
    x: &Expr,
    cache: &mut FxHashMap<usize, (Expr, Expr)>,
    options: &ExpandOptions,
    opaque: &mut dyn FnMut(&Expr) -> Option<Expr>,
) -> Expr {
    // x**(a+b) -> x**a * x**b, only when the base is known to be non-zero
    // (SymPy's power_exp hint without assumptions)
    let nonzero_base = match b.kind() {
        Kind::Num(n) => !n.is_zero(),
        Kind::Const(c) => !c.is_infinite() && *c != crate::expr::Constant::NaN,
        _ => false,
    };
    if let (Kind::Add(..), true) = (x.kind(), nonzero_base && options.power_exp) {
        let parts: Vec<Expr> = x.args().iter().map(|t| pow(b, t)).collect();
        let m = mul(&parts);
        return if m.is_mul() {
            let factors: Vec<Expr> = m
                .args()
                .iter()
                .map(|f| ex(f, cache, options, opaque))
                .collect();
            if options.mul && factors.iter().any(|f| f.is_add()) {
                expand_product(&factors)
            } else {
                mul(&factors)
            }
        } else if let Kind::Pow(b3, x3) = m.kind() {
            if b3.is_add() && x3.is_integer_number() {
                expand_pow(b3, x3, cache, options, opaque)
            } else {
                m
            }
        } else {
            m
        };
    }
    if let (Kind::Add(..), Some(n), true) = (b.kind(), x.as_number(), options.multinomial) {
        if let Some(k) = n.to_i64() {
            if k > 0 {
                return u32::try_from(k)
                    .map(|k| power_poly(b, k))
                    .unwrap_or_else(|_| pow(b, x));
            } else if k < 0 {
                return u32::try_from(k.unsigned_abs())
                    .map(|k| pow(&power_poly(b, k), &Expr::minus_one()))
                    .unwrap_or_else(|_| pow(b, x));
            }
        } else if n.is_rational() && !n.is_integer() {
            // (x+1)**(5/2) -> (x+1)**2 * sqrt(x+1), expanded
            let k = n.floor();
            if let Some(k) = k.to_i64() {
                if k >= 1 {
                    let frac = Expr::num(n.sub(&Number::int(k)));
                    let radical = pow(b, &frac);
                    let Ok(k) = u32::try_from(k) else {
                        return pow(b, x);
                    };
                    let factors = [power_poly(b, k), radical];
                    // Multinomial expansion includes this distribution even
                    // when the independent multiplication hint is disabled.
                    return expand_product(&factors);
                }
            }
        }
    }
    // SymPy's Pow._eval_expand_power_base of a product that does not
    // commute: an integer power repeats the ordered factors.
    if let (Kind::Mul(..), Some(k)) = (b.kind(), x.as_number().filter(|n| n.is_integer()).and_then(Number::to_i64)) {
        if options.power_base && b.has_noncommutative_host() && crate::decompose::commutativity(b) != Some(true) {
            let (commutative, ordered) = crate::decompose::args_cnc(b, false);
            let copies = usize::try_from(k.unsigned_abs()).ok().filter(|n| n.saturating_mul(ordered.len()) <= 1 << 16);
            if let Some(copies) = copies {
                let unit: Vec<Expr> = if k > 0 { ordered } else { ordered.iter().rev().map(|f| pow(f, &Expr::minus_one())).collect() };
                let mut factors: Vec<Expr> = (0..copies).flat_map(|_| unit.iter().cloned()).collect();
                if !commutative.is_empty() { factors.insert(0, pow(&mul(&commutative), x)); }
                return if options.mul && factors.iter().any(|f| f.is_add()) { expand_product(&factors) } else { mul(&factors) };
            }
        }
    }
    let p = pow(b, x);
    // (a*b)**2 may become a product of expandable powers
    if let Kind::Mul(..) = p.kind() {
        if !options.deep
            || !options.multinomial
            || !matches!(b.kind(), Kind::Mul(..))
            || !x.is_integer_number()
        {
            return p;
        }
        let factors: Vec<Expr> = p
            .args()
            .iter()
            .map(|f| ex(f, cache, options, opaque))
            .collect();
        if options.mul && factors.iter().any(|f| f.is_add()) {
            return expand_product(&factors);
        }
        return mul(&factors);
    }
    p
}

// ----------------------------------------------------------------------
// sparse polynomial machinery
// ----------------------------------------------------------------------

struct Gens {
    list: Vec<Expr>,
    index: FxHashMap<Expr, usize>,
}

impl Gens {
    fn new() -> Self {
        Gens {
            list: Vec::new(),
            index: FxHashMap::default(),
        }
    }
    fn id(&mut self, g: Expr) -> usize {
        if let Some(i) = self.index.get(&g) {
            return *i;
        }
        let i = self.list.len();
        self.index.insert(g.clone(), i);
        self.list.push(g);
        i
    }
}

/// A term as (coefficient, [(gen id, exponent)]).
type RawTerm = (Number, SmallVec<[(usize, i32); 8]>);

fn gen_of(b: &Expr, x: &Expr, gens: &mut Gens) -> (usize, i32) {
    if let Some(n) = x.as_number() {
        if let Some(k) = n.to_i64() {
            if let Ok(k) = i32::try_from(k) {
                return (gens.id(b.clone()), k);
            }
        }
        if let Number::Rational(r) = n {
            if let (Ok(p), true) = (
                i32::try_from(num_traits::ToPrimitive::to_i64(r.numer()).unwrap_or(i64::MAX)),
                r.denom() < &num_bigint::BigInt::from(1 << 20),
            ) {
                let q = num_traits::ToPrimitive::to_i64(r.denom()).unwrap();
                let g = Expr::from_kind(Kind::Pow(b.clone(), Expr::rational(1, q)));
                return (gens.id(g), p);
            }
        }
    }
    (gens.id(mul_factor(b, x)), 1)
}

fn raw_terms(e: &Expr, gens: &mut Gens) -> Vec<RawTerm> {
    let one_term = |t: &Expr, k: &Number, gens: &mut Gens| -> RawTerm {
        let mut m: SmallVec<[(usize, i32); 8]> = SmallVec::new();
        let mut c = k.clone();
        match t.kind() {
            Kind::Num(n) => c = c.mul(n),
            Kind::Mul(mc, facs) => {
                c = c.mul(mc);
                for (b, x) in facs {
                    m.push(gen_of(b, x, gens));
                }
            }
            Kind::Pow(b, x) => m.push(gen_of(b, x, gens)),
            _ => m.push((gens.id(t.clone()), 1)),
        }
        (c, m)
    };
    match e.kind() {
        Kind::Add(c, terms) => {
            let mut v: Vec<RawTerm> = terms.iter().map(|(t, k)| one_term(t, k, gens)).collect();
            if !c.is_zero() {
                v.push((c.clone(), SmallVec::new()));
            }
            v
        }
        _ => vec![one_term(e, &Number::ONE, gens)],
    }
}

type Poly = FxHashMap<Mono, Number>;

fn to_poly(raw: Vec<RawTerm>, n: usize) -> Poly {
    let mut p: Poly = FxHashMap::with_capacity_and_hasher(raw.len(), Default::default());
    for (c, m) in raw {
        let mut mono: Mono = SmallVec::from_elem(0, n);
        for (g, k) in m {
            mono[g] += i64::from(k);
        }
        match p.get_mut(&mono) {
            Some(v) => *v = v.add(&c),
            None => {
                p.insert(mono, c);
            }
        }
    }
    p.retain(|_, c| !c.is_zero());
    p
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    let (a, b) = if a.len() < b.len() { (b, a) } else { (a, b) };
    let mut out: Poly = FxHashMap::with_capacity_and_hasher(
        a.len().saturating_mul(b.len()).min(1 << 22),
        Default::default(),
    );
    let mut mono: Mono = SmallVec::new();
    for (mb, cb) in b {
        for (ma, ca) in a {
            mono.clear();
            mono.extend(ma.iter().zip(mb.iter()).map(|(x, y)| x + y));
            let c = ca.mul(cb);
            match out.get_mut(&mono) {
                Some(v) => *v = v.add(&c),
                None => {
                    out.insert(mono.clone(), c);
                }
            }
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

fn poly_pow(p: &Poly, n: u32) -> Poly {
    if n == 1 {
        return p.clone();
    }
    // repeated multiplication by the (small) base beats squaring for
    // sparse inputs
    let mut r = p.clone();
    for _ in 1..n {
        r = poly_mul(&r, p);
    }
    r
}

/// Generators that can be assembled into a `Mul` without re-canonicalising.
fn is_simple_gen(g: &Expr) -> bool {
    match g.kind() {
        Kind::Sym(_) | Kind::Add(..) => true,
        Kind::Func(f, _) => !matches!(f, Func::Exp),
        _ => false,
    }
}

fn from_poly(p: Poly, gens: &Gens) -> Expr {
    let n = gens.list.len();
    if gens.list.iter().all(is_simple_gen) {
        // gens in canonical order: monomial factors come out sorted
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&i, &j| gens.list[i].cmp(&gens.list[j]));
        let mut coeff = Number::ZERO;
        let mut extra: Vec<Expr> = Vec::new();
        let mut terms: Vec<(Expr, Number)> = Vec::with_capacity(p.len());
        for (m, c) in p {
            let mut facs: Vec<(Expr, Expr)> = Vec::new();
            for &i in &order {
                if m[i] != 0 {
                    facs.push((gens.list[i].clone(), Expr::int(m[i] as i64)));
                }
            }
            let t = match facs.len() {
                0 => {
                    coeff = coeff.add(&c);
                    continue;
                }
                1 if facs[0].0.is_add() && facs[0].1.is_one() => {
                    // a bare sum must be flattened into the result
                    extra.push(mul(&[Expr::num(c), facs[0].0.clone()]));
                    continue;
                }
                1 => mul_factor(&facs[0].0, &facs[0].1),
                _ => Expr::from_kind(Kind::Mul(Number::ONE, facs)),
            };
            terms.push((t, c));
        }
        let s = add_from_terms(coeff, terms);
        if extra.is_empty() {
            return s;
        }
        extra.push(s);
        return add(&extra);
    }
    let mut v: Vec<Expr> = Vec::with_capacity(p.len());
    for (m, c) in p {
        let mut f: Vec<Expr> = Vec::with_capacity(n + 1);
        f.push(Expr::num(c));
        for (i, k) in m.iter().enumerate() {
            if *k != 0 {
                f.push(pow(&gens.list[i], &Expr::int(*k as i64)));
            }
        }
        v.push(mul(&f));
    }
    // products of radicals may produce new sums (e.g. I*I = -1): flatten
    add(&v)
}

/// SymPy's `fraction(expr)` (non-exact) for a product: numerator and
/// denominator factors.
fn fraction_parts(m: &Expr) -> (Vec<Expr>, Vec<Expr>) {
    let mut num = Vec::new();
    let mut den = Vec::new();
    for t in m.args() {
        let (b, x) = match t.kind() {
            Kind::Pow(b, x) => (b.clone(), x.clone()),
            Kind::Func(Func::Exp, a) => (crate::expr::Expr::e(), a[0].clone()),
            Kind::Num(n) if n.is_rational() && !n.is_integer() => {
                let p = n.numer();
                if !p.is_one() {
                    num.push(Expr::num(p));
                }
                den.push(Expr::num(n.denom()));
                continue;
            }
            _ => {
                num.push(t);
                continue;
            }
        };
        let rebuild = |b: &Expr, x: &Expr| {
            if b.is_constant(crate::expr::Constant::E) {
                exp(x)
            } else {
                pow(b, x)
            }
        };
        match x.as_number() {
            Some(n) if n.is_negative() => {
                den.push(if n.is_minus_one() {
                    b
                } else {
                    rebuild(&b, &Expr::num(n.neg()))
                });
            }
            Some(_) => num.push(t),
            None => match x.kind() {
                Kind::Mul(c, _) if c.is_negative() => {
                    den.push(rebuild(&b, &crate::ops::neg(&x)));
                }
                _ => num.push(t),
            },
        }
    }
    (num, den)
}

/// SymPy's `Mul._eval_expand_mul`: a product of several denominators is
/// expanded into a single one, then sums in the numerator are distributed.
fn expand_mul_node(m: &Expr) -> Expr {
    let (num, den) = fraction_parts(m);
    let dm = mul(&den);
    let expr = if dm.is_mul() {
        let n = mul(&num);
        let n2 = if n.is_mul() { expand_mul_node(&n) } else { n };
        let d2 = expand_mul_node(&dm);
        crate::ops::div(&n2, &d2)
    } else if den.len() >= 2 {
        // the reference rebuilds `n/d` whatever `d` became: a number
        // distributed over a sum denominator, `1/(2*(x/2 + 1))`, leaves one
        crate::ops::div(&mul(&num), &dm)
    } else {
        m.clone()
    };
    if !expr.is_mul() {
        return expr;
    }
    let factors = expr.args();
    if factors.iter().any(|f| f.is_add()) {
        let r = expand_product(&factors);
        // SymPy repeats the expansion until a fixed point: products of
        // several denominators created by the distribution are combined
        if let Kind::Add(..) = r.kind() {
            if r.args().iter().any(has_multi_denominator) {
                let parts: Vec<Expr> = r
                    .args()
                    .iter()
                    .map(|t| {
                        if has_multi_denominator(t) {
                            expand_mul_node(t)
                        } else {
                            t.clone()
                        }
                    })
                    .collect();
                return add(&parts);
            }
        }
        r
    } else {
        expr
    }
}

/// Has the term a sum raised to a negative power among its factors?
fn sum_denominator(t: &Expr) -> bool {
    let neg = |b: &Expr, x: &Expr| b.is_add() && x.as_number().is_some_and(|n| n.is_negative());
    match t.kind() {
        Kind::Pow(b, x) => neg(b, x),
        Kind::Mul(_, facs) => facs.iter().any(|(b, x)| neg(b, x)),
        _ => false,
    }
}

/// Integer powers of sums left in a product or a sum of products.
fn has_unexpanded_power(e: &Expr) -> bool {
    let bad = |b: &Expr, x: &Expr| b.is_add() && x.as_i64().is_some_and(|k| k.unsigned_abs() >= 2);
    match e.kind() {
        Kind::Pow(b, x) => bad(b, x),
        Kind::Mul(_, facs) => facs.iter().any(|(b, x)| bad(b, x)),
        Kind::Add(_, terms) => terms.iter().any(|(t, _)| match t.kind() {
            Kind::Mul(_, facs) => facs.iter().any(|(b, x)| bad(b, x)),
            Kind::Pow(b, x) => bad(b, x),
            _ => false,
        }),
        _ => false,
    }
}

/// Does the product `t` have a denominator made of several factors?
fn has_multi_denominator(t: &Expr) -> bool {
    if !t.is_mul() {
        return false;
    }
    let (_, den) = fraction_parts(t);
    den.len() >= 2
}

/// Expand a product of already expanded factors.
pub fn expand_product(factors: &[Expr]) -> Expr {
    if factors.iter().any(Expr::has_noncommutative_host)
        && factors.iter().any(|factor| crate::decompose::commutativity(factor) != Some(true)) {
        return noncommutative_product(factors);
    }
    if factors.iter().any(|f| f.is_add() && f.float_info().max_precision.is_some()) {
        return float_product(factors);
    }
    let mut gens = Gens::new();
    let raws: Vec<Vec<RawTerm>> = factors.iter().map(|f| raw_terms(f, &mut gens)).collect();
    let n = gens.list.len();
    let mut acc: Option<Poly> = None;
    for r in raws {
        let p = to_poly(r, n);
        acc = Some(match acc {
            None => p,
            Some(a) => poly_mul(&a, &p),
        });
    }
    let r = from_poly(acc.unwrap_or_default(), &gens);
    if gens.list.iter().all(is_simple_gen) {
        r
    } else {
        // radical products may need a second pass ((1+sqrt(2))*(1-sqrt(2)))
        renormalize(&r)
    }
}

/// A product with Float sums, distributed in the order of the reference
/// `Mul._expandsums`: the sums are split in halves, each half is expanded
/// and the pairwise products are added at once. Float coefficients are
/// rounded at every addition, so a different association changes the last
/// bits of the result (and the roots of ill-conditioned polynomials).
fn float_product(factors: &[Expr]) -> Expr {
    fn sums_terms(sums: &[Expr]) -> Vec<Expr> {
        if sums.len() == 1 {
            return crate::traverse::ordered_args(&sums[0]);
        }
        let half = sums.len() / 2;
        let (left, right) = (sums_terms(&sums[..half]), sums_terms(&sums[half..]));
        let terms: Vec<Expr> = left.iter().flat_map(|a| right.iter().map(move |b| crate::ops::mul2(a, b))).collect();
        let added = add(&terms);
        if added.is_add() { crate::traverse::ordered_args(&added) } else { vec![added] }
    }
    // The reference iterates the arguments of the rebuilt product, in their
    // canonical (printed) order, not the internal storage order.
    let ordered = mul(factors);
    if !ordered.is_mul() {
        return ordered;
    }
    let (sums, plain): (Vec<Expr>, Vec<Expr>) = crate::traverse::ordered_args(&ordered).into_iter().partition(Expr::is_add);
    let plain = mul(&plain);
    let terms: Vec<Expr> = sums_terms(&sums)
        .iter()
        .map(|term| {
            let t = crate::ops::mul2(&plain, term);
            if t.is_mul() && t.args().iter().any(Expr::is_add) { expand_mul_node(&t) } else { t }
        })
        .collect();
    add(&terms)
}

/// Distribute a product whose factors may not commute, keeping the order
/// of the factors in every term (SymPy's `Mul._eval_expand_mul`).
fn noncommutative_product(factors: &[Expr]) -> Expr {
    let mut terms = vec![Expr::one()];
    for factor in factors {
        let summands = if factor.is_add() { factor.args() } else { vec![factor.clone()] };
        terms = terms.iter().flat_map(|term| summands.iter().map(move |summand| crate::ops::mul2(term, summand))).collect();
    }
    add(&terms)
}

/// Expand `b**n` for a sum `b` and positive integer `n`.
pub fn power_poly(b: &Expr, n: u32) -> Expr {
    if b.has_noncommutative_host() && crate::decompose::commutativity(b) != Some(true) {
        return noncommutative_product(&vec![b.clone(); n as usize]);
    }
    if n >= 2 && b.is_add() && b.float_info().max_precision.is_some() && !b.args().iter().any(crate::series::order::is_order) {
        // The reference `basic_from_dict(multinomial_coefficients(m, n), *terms)`:
        // each term is built from powers of the Float summands, so their
        // rounding follows the reference instead of repeated multiplication.
        let p = crate::traverse::ordered_args(b);
        let terms: Vec<Expr> = crate::poly::inexact_ground::multinomial_coefficients(p.len(), n as usize)
            .into_iter()
            .map(|(ks, c)| {
                let mut f = vec![Expr::num(Number::from_bigint(c))];
                f.extend(p.iter().zip(ks).filter(|(_, k)| *k > 0).map(|(g, k)| pow(g, &Expr::int(k as i64))));
                mul(&f)
            })
            .collect();
        return add(&terms);
    }
    let mut gens = Gens::new();
    let raw = raw_terms(b, &mut gens);
    let p = to_poly(raw, gens.list.len());
    let r = from_poly(poly_pow(&p, n), &gens);
    if gens.list.iter().all(is_simple_gen) {
        r
    } else {
        renormalize(&r)
    }
}

fn renormalize(e: &Expr) -> Expr {
    // terms whose product produced a sum (e.g. sqrt(2)*(..)) are distributed
    if let Kind::Add(..) = e.kind() {
        let parts: Vec<Expr> = e
            .args()
            .iter()
            .map(|t| {
                if let Kind::Mul(_, facs) = t.kind() {
                    if facs.iter().any(|(b, x)| b.is_add() && x.is_one()) {
                        return expand(t);
                    }
                }
                t.clone()
            })
            .collect();
        return add(&parts);
    }
    e.clone()
}
