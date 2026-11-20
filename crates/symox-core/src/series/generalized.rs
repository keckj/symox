//! Sparse series with exact real exponents beyond finite ramification.
//!
//! The ordinary engine remains the first choice. This conservative fallback
//! uses certified exponent comparisons and tracks absolute error valuations
//! through products and analytic powers. It only expands along x -> 0+;
//! branches and undecidable exponent orders never use floating point guesses.

use std::cmp::Ordering;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive};

use crate::expr::{Expr, Func, Kind};
use crate::functions as fx;
use crate::number::Number;
use crate::ops::{add, div, mul2, neg, pow, sub};
use super::engine::depends;

const TERMS: usize = 256;
const DEPTH: usize = 80;

/// Exact sign of a rational linear combination of real rational logarithms.
/// Clear coefficient denominators, then use monotonicity of log to compare
/// the resulting positive rational product with one.
fn log_sum_sign(e: &Expr) -> Option<i32> {
    let terms = if e.is_add() { e.args() } else { vec![e.clone()] };
    let mut logs = Vec::new();
    let mut denominator = BigInt::one();
    for term in terms {
        let (coefficient, head) = match term.kind() {
            Kind::Func(Func::Log, _) => (Number::ONE, term.clone()),
            Kind::Mul(c, fs) if fs.len() == 1 && fs[0].1.is_one() => (c.clone(), fs[0].0.clone()),
            _ => return None,
        };
        let Kind::Func(Func::Log, args) = head.kind() else { return None };
        let [argument] = args.as_slice() else { return None };
        let base = argument.as_number()?.to_ratio()?;
        if !base.is_positive() { return None; }
        let (p, q) = coefficient.numer_denom()?;
        denominator = denominator.lcm(&q);
        if denominator.bits() > 32 { return None; }
        logs.push((argument.clone(), p, q));
    }
    let mut product = Expr::one();
    let mut size = 0u64;
    for (base, p, q) in logs {
        let k = (p * (&denominator / q)).to_i64()?;
        if k.unsigned_abs() > 2048 { return None; }
        let (bn, bd) = base.as_number()?.numer_denom()?;
        size = size.checked_add((bn.bits() + bd.bits()).checked_mul(k.unsigned_abs())?)?;
        if size > 65536 { return None; }
        product = mul2(&product, &pow(&base, &Expr::int(k)));
    }
    Some(match product.as_number()?.cmp_value(&Number::ONE)? {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    })
}


// Cached numeric logarithms depend only on exact rational atoms and the
// enclosure precision. The cache is bounded and never stores assumption facts.
thread_local! {
    static LOG_INTERVALS: std::cell::RefCell<rustc_hash::FxHashMap<(Number, u32), crate::poly::enclosure::Iv>> =
        std::cell::RefCell::new(rustc_hash::FxHashMap::default());
}

/// For 0 <= t <= 1/3, log((1+t)/(1-t)) is the atanh series.
/// After m terms its positive remainder is bounded by
/// 2*t**(2m+1)/((2m+1)*(1-t*t)) <= 3*9**(-m)/(4*(2m+1))
/// < 2**(-3m). All included terms and operations round outward.
fn atanh_log_interval(t: &crate::poly::int::Rat, precision: u32) -> crate::poly::enclosure::Iv {
    use crate::poly::enclosure::Iv;
    use crate::poly::int::Rat;
    let work = precision + 32;
    let count = ((precision + 34) / 3) as usize;
    let t = Iv::from_rat(t, work);
    let square = t.mul(&t, work);
    let mut power = t.clone();
    let mut sum = Iv::zero();
    for j in 0..count {
        let coefficient = Rat::from_number(&Number::rational_i64(1, (2*j + 1) as i64).unwrap()).unwrap();
        sum = sum.add(&power.mul_rat(&coefficient, work), work);
        power = power.mul(&square, work);
    }
    let remainder = Iv { lo: BigInt::from(0), hi: BigInt::one(), e: -3 * count as i64 };
    sum.mul_rat(&Rat::from_i64(2), work).add(&remainder, work)
}

/// Range-reduce a positive rational r = 2**k * u, 1 <= u < 2.
/// Both log(u) and log(2) use the certified positive atanh tail above.
fn rational_log_interval(number: &Number, precision: u32) -> Option<crate::poly::enclosure::Iv> {
    use crate::poly::int::Rat;
    let key = (number.clone(), precision);
    if let Some(answer) = LOG_INTERVALS.with(|cache| cache.borrow().get(&key).cloned()) { return Some(answer); }
    let (mut numerator, mut denominator) = number.numer_denom()?;
    if !numerator.is_positive() || numerator.bits() + denominator.bits() > 65536 { return None; }
    let mut k = numerator.bits() as i64 - denominator.bits() as i64;
    if k.unsigned_abs() > 4096 { return None; }
    if k >= 0 { denominator <<= k as usize; } else { numerator <<= (-k) as usize; }
    if numerator < denominator { numerator <<= 1usize; k -= 1; }
    if k.unsigned_abs() > 4096 { return None; }
    let t = Number::rational(&numerator - &denominator, &numerator + &denominator)?;
    let mut answer = atanh_log_interval(&Rat::from_number(&t)?, precision);
    if k != 0 {
        let third = Rat::from_number(&Number::rational_i64(1, 3)?)?;
        let log_two = atanh_log_interval(&third, precision);
        answer = answer.add(&log_two.mul_rat(&Rat::from_i64(k), precision + 32), precision + 32);
    }
    LOG_INTERVALS.with(|cache| {
        let mut cache = cache.borrow_mut();
        if cache.len() >= 64 { cache.clear(); }
        cache.insert(key, answer.clone());
    });
    Some(answer)
}

fn integer_interval_power(base: &crate::poly::enclosure::Iv, k: i64, precision: u32) -> Option<crate::poly::enclosure::Iv> {
    use crate::poly::enclosure::Iv;
    if k.unsigned_abs() > 2048 { return None; }
    let mut base = if k < 0 { base.inv(precision)? } else { base.clone() };
    let mut power = k.unsigned_abs();
    let mut result = Iv::int(BigInt::one());
    while power != 0 {
        if power & 1 != 0 { result = result.mul(&base, precision); }
        power >>= 1;
        if power != 0 { base = base.mul(&base, precision); }
    }
    Some(result)
}

/// Real expressions generated by exact rationals, positive rational logs,
/// algebraic constants, and field arithmetic; no numerical midpoint tests.
fn logarithmic_interval(e: &Expr, precision: u32) -> Option<crate::poly::enclosure::Iv> {
    use crate::poly::enclosure::Iv;
    use crate::poly::int::Rat;
    fn visit(e: &Expr, precision: u32, depth: usize, cache: &mut rustc_hash::FxHashMap<Expr, Iv>) -> Option<Iv> {
        if depth > 64 || cache.len() > 1024 { return None; }
        if let Some(value) = cache.get(e) { return Some(value.clone()); }
        let algebraic = || {
            let enclosure = crate::poly::enclosure::enclose(e, precision).ok()?;
            enclosure.im.is_exact_zero().then_some(enclosure.re)
        };
        let result = match e.kind() {
            Kind::Num(n) => Iv::from_rat(&Rat::from_number(n)?, precision),
            Kind::Add(c, terms) => {
                let mut result = Iv::from_rat(&Rat::from_number(c)?, precision);
                for (term, coefficient) in terms {
                    result = result.add(&visit(term, precision, depth + 1, cache)?.mul_rat(&Rat::from_number(coefficient)?, precision), precision);
                }
                result
            }
            Kind::Mul(c, factors) => {
                let mut result = Iv::from_rat(&Rat::from_number(c)?, precision);
                for (base, exponent) in factors {
                    let factor = if let Some(k) = exponent.as_i64() {
                        integer_interval_power(&visit(base, precision, depth + 1, cache)?, k, precision)?
                    } else {
                        visit(&pow(base, exponent), precision, depth + 1, cache)?
                    };
                    result = result.mul(&factor, precision);
                }
                result
            }
            Kind::Pow(base, exponent) => {
                if let Some(k) = exponent.as_i64() {
                    integer_interval_power(&visit(base, precision, depth + 1, cache)?, k, precision)?
                } else { algebraic()? }
            }
            Kind::Func(Func::Log, args) if args.len() == 1 => rational_log_interval(args[0].as_number()?, precision)?,
            _ => algebraic()?,
        };
        cache.insert(e.clone(), result.clone());
        Some(result)
    }
    visit(e, precision, 0, &mut rustc_hash::FxHashMap::default())
}

fn sign(e: &Expr, depth: usize) -> Option<i32> {
    if depth > 12 { return None; }
    let e = crate::expand::expand(e);
    if e.is_zero() { return Some(0); }
    if let Some(n) = e.as_number().filter(|n| n.is_rational()) {
        return Some(if n.is_negative() { -1 } else { 1 });
    }
    // The exact-query bank excludes numerical sampling and tolerance tests.
    let structural = crate::assumptions::exact_queries(|| {
        if crate::assumptions::is_positive(&e) == Some(true) { Some(1) }
        else if crate::assumptions::is_negative(&e) == Some(true) { Some(-1) }
        else if crate::assumptions::is_zero(&e) == Some(true) { Some(0) }
        else { None }
    });
    if structural.is_some() { return structural; }
    if let Some(s) = log_sum_sign(&e) { return Some(s); }
    let (numerator, denominator) = crate::simplify::numer_denom(&e);
    if !denominator.is_one() && numerator != e {
        let d = sign(&denominator, depth + 1)?;
        if d == 0 { return None; }
        return sign(&numerator, depth + 1).map(|n| n * d);
    }
    if let Kind::Mul(c, factors) = e.kind() {
        let mut answer = if c.is_negative() { -1 } else { 1 };
        for (base, exponent) in factors {
            let k = exponent.as_i64()?;
            let s = sign(base, depth + 1)?;
            if s == 0 && k < 0 { return None; }
            answer *= if k % 2 == 0 { i32::from(s != 0) } else { s };
        }
        return Some(answer);
    }
    if !crate::traverse::free_symbols(&e).is_empty() { return None; }
    for precision in [64, 128, 256] {
        if let Some(bounds) = logarithmic_interval(&e, precision) {
            if bounds.positive() { return Some(1); }
            if bounds.negative() { return Some(-1); }
            if bounds.is_exact_zero() { return Some(0); }
        }
    }
    None
}

fn compare(a: &Expr, b: &Expr) -> Option<Ordering> {
    Some(match sign(&sub(a, b), 0)? {
        -1 => Ordering::Less,
        0 => Ordering::Equal,
        _ => Ordering::Greater,
    })
}

fn real_exponent(e: &Expr, x: &Expr) -> bool {
    !depends(e, x) && !e.has_nonfinite() && crate::assumptions::exact_queries(|| {
        crate::assumptions::is_real(e) == Some(true)
            && crate::assumptions::is_finite(e) == Some(true)
    })
}

fn polynomial(e: &Expr, x: &Expr) -> bool {
    if !depends(e, x) || e == x { return true; }
    let power = |base: &Expr, exponent: &Expr| {
        exponent.as_i64().is_some_and(|k| k >= 0) && polynomial(base, x)
    };
    match e.kind() {
        Kind::Add(_, terms) => terms.iter().all(|(term, _)| polynomial(term, x)),
        Kind::Mul(_, factors) => factors.iter().all(|(base, exponent)| power(base, exponent)),
        Kind::Pow(base, exponent) => power(base, exponent),
        _ => false,
    }
}

#[derive(Clone)]
struct Sparse {
    /// Ordered (exponent, coefficient) pairs; no zero coefficients.
    terms: Vec<(Expr, Expr)>,
    /// None is exact. Otherwise the error is O(x**precision).
    precision: Option<Expr>,
}

impl Sparse {
    fn zero() -> Self { Self { terms: Vec::new(), precision: None } }
    fn monomial(exponent: Expr, coefficient: Expr) -> Self {
        if coefficient.is_zero() { return Self::zero(); }
        Self { terms: vec![(exponent, coefficient)], precision: None }
    }
    fn constant(coefficient: Expr) -> Self { Self::monomial(Expr::zero(), coefficient) }

    fn insert(&mut self, exponent: Expr, coefficient: Expr) -> Option<()> {
        if coefficient.is_zero() { return Some(()); }
        for i in 0..self.terms.len() {
            match compare(&exponent, &self.terms[i].0)? {
                Ordering::Less => {
                    if self.terms.len() >= TERMS { return None; }
                    self.terms.insert(i, (exponent, coefficient));
                    return Some(());
                }
                Ordering::Equal => {
                    self.terms[i].1 = crate::expand::expand(&add(&[self.terms[i].1.clone(), coefficient]));
                    if self.terms[i].1.is_zero() { self.terms.remove(i); }
                    return Some(());
                }
                Ordering::Greater => {}
            }
        }
        if self.terms.len() >= TERMS { return None; }
        self.terms.push((exponent, coefficient));
        Some(())
    }

    fn restrict_precision(&mut self, other: Option<Expr>) -> Option<()> {
        if let Some(p) = other {
            if match &self.precision { None => true, Some(q) => compare(&p, q)? == Ordering::Less } {
                self.precision = Some(p);
            }
        }
        Some(())
    }

    fn valuation(&self) -> Option<Expr> {
        match (self.terms.first(), &self.precision) {
            (Some((v, _)), Some(p)) => Some(if compare(v, p)? == Ordering::Less { v.clone() } else { p.clone() }),
            (Some((v, _)), None) => Some(v.clone()),
            (None, p) => p.clone(),
        }
    }

    fn known_lead(&self) -> Option<(Expr, Expr)> {
        let (valuation, coefficient) = self.terms.first()?.clone();
        if let Some(p) = &self.precision {
            if compare(&valuation, p)? != Ordering::Less { return None; }
        }
        Some((valuation, coefficient))
    }

    fn clip(mut self, cutoff: &Expr, x: &Expr) -> Option<Self> {
        let mut keep = Vec::new();
        let mut dropped = false;
        for (exponent, coefficient) in self.terms {
            if compare(&exponent, cutoff)? == Ordering::Less {
                keep.push((exponent, coefficient));
            } else {
                // Logarithmic coefficients need a separate error scale.
                if depends(&coefficient, x) { return None; }
                dropped = true;
            }
        }
        self.terms = keep;
        if dropped { self.restrict_precision(Some(cutoff.clone()))?; }
        Some(self)
    }

    fn add(&self, other: &Self) -> Option<Self> {
        let mut result = self.clone();
        result.restrict_precision(other.precision.clone())?;
        for (exponent, coefficient) in &other.terms { result.insert(exponent.clone(), coefficient.clone())?; }
        Some(result)
    }

    fn scale(&self, coefficient: &Expr, shift: &Expr, x: &Expr) -> Option<Self> {
        if coefficient.is_zero() { return Some(Self::zero()); }
        if self.precision.is_some() && depends(coefficient, x) { return None; }
        Some(Self {
            terms: self.terms.iter().map(|(a, c)| (add(&[a.clone(), shift.clone()]), mul2(c, coefficient))).collect(),
            precision: self.precision.as_ref().map(|p| add(&[p.clone(), shift.clone()])),
        })
    }

    fn mul(&self, other: &Self, cutoff: &Expr, x: &Expr) -> Option<Self> {
        if self.terms.is_empty() && self.precision.is_none() || other.terms.is_empty() && other.precision.is_none() {
            return Some(Self::zero());
        }
        // Multiplying an error by an unbounded log coefficient cannot retain
        // its bare power Order. Ordinary logarithmic series keep their engine.
        if self.precision.is_some() && other.terms.iter().any(|(_, c)| depends(c, x))
            || other.precision.is_some() && self.terms.iter().any(|(_, c)| depends(c, x)) { return None; }
        let mut result = Self::zero();
        if let Some(p) = &self.precision { result.restrict_precision(Some(add(&[p.clone(), other.valuation()?])))?; }
        if let Some(p) = &other.precision { result.restrict_precision(Some(add(&[p.clone(), self.valuation()?])))?; }
        for (a, c) in &self.terms {
            for (b, d) in &other.terms { result.insert(add(&[a.clone(), b.clone()]), mul2(c, d))?; }
        }
        result.clip(cutoff, x)
    }
}

struct Expansion<'a> { x: &'a Expr, cutoff: Expr }

impl Expansion<'_> {
    fn analytic_power(&self, base: Sparse, exponent: &Expr) -> Option<Sparse> {
        let (valuation, coefficient) = base.known_lead()?;
        if depends(&coefficient, self.x) || crate::assumptions::exact_queries(|| crate::assumptions::is_zero(&coefficient)) != Some(false) {
            return None;
        }
        let shift = mul2(&valuation, exponent);
        if !real_exponent(&shift, self.x) { return None; }
        // Integer powers preserve every constant branch. Other powers need a
        // positive real leading coefficient on this positive real path.
        if !exponent.as_number().is_some_and(|n| n.is_integer()) && sign(&coefficient, 0)? != 1 { return None; }
        let prefactor = pow(&coefficient, exponent);
        let normalized = base.scale(&div(&Expr::one(), &coefficient), &neg(&valuation), self.x)?;
        let mut tail = normalized.add(&Sparse::constant(Expr::minus_one()))?;
        let relative = sub(&self.cutoff, &shift);
        if tail.terms.is_empty() && tail.precision.is_none() { return Some(Sparse::monomial(shift, prefactor)); }
        if sign(&tail.valuation()?, 0)? != 1 { return None; }
        // Expand (1 + tail)**p using its convergent local binomial series.
        // Each coefficient is exact; the first omitted power bounds the error.
        tail = tail.clip(&relative, self.x)?;
        let mut result = Sparse::constant(Expr::one());
        let mut power = Sparse::constant(Expr::one());
        let mut coefficient = Expr::one();
        for k in 1..=TERMS {
            power = power.mul(&tail, &relative, self.x)?;
            coefficient = div(&mul2(&coefficient, &sub(exponent, &Expr::int(k as i64 - 1))), &Expr::int(k as i64));
            if coefficient.is_zero() { return result.scale(&prefactor, &shift, self.x); }
            if power.terms.is_empty() {
                result.restrict_precision(power.precision.clone())?;
                return result.scale(&prefactor, &shift, self.x);
            }
            result = result.add(&power.scale(&coefficient, &Expr::zero(), self.x)?)?;
        }
        None
    }

    fn logarithm(&self, base: Sparse) -> Option<Sparse> {
        let (valuation, coefficient) = base.known_lead()?;
        if depends(&coefficient, self.x) || sign(&coefficient, 0)? != 1 { return None; }
        let main = add(&[fx::log(&coefficient), mul2(&valuation, &fx::log(self.x))]);
        let normalized = base.scale(&div(&Expr::one(), &coefficient), &neg(&valuation), self.x)?;
        let tail = normalized.add(&Sparse::constant(Expr::minus_one()))?;
        let mut result = Sparse::constant(main);
        if tail.terms.is_empty() && tail.precision.is_none() { return Some(result); }
        if sign(&tail.valuation()?, 0)? != 1 { return None; }
        let mut power = Sparse::constant(Expr::one());
        for k in 1..=TERMS {
            power = power.mul(&tail, &self.cutoff, self.x)?;
            if power.terms.is_empty() {
                result.restrict_precision(power.precision.clone())?;
                return Some(result);
            }
            let coefficient = Expr::rational(if k % 2 == 0 { -1 } else { 1 }, k as i64);
            result = result.add(&power.scale(&coefficient, &Expr::zero(), self.x)?)?;
        }
        None
    }

    fn expand(&self, e: &Expr, depth: usize) -> Option<Sparse> {
        if depth > DEPTH { return None; }
        if !depends(e, self.x) { return Some(Sparse::constant(e.clone())); }
        if e == self.x { return Some(Sparse::monomial(Expr::one(), Expr::one())); }
        let result = match e.kind() {
            Kind::Add(c, terms) => {
                let mut result = Sparse::constant(Expr::num(c.clone()));
                for (term, coefficient) in terms {
                    let expanded = self.expand(term, depth + 1)?.scale(&Expr::num(coefficient.clone()), &Expr::zero(), self.x)?;
                    result = result.add(&expanded)?;
                }
                result
            }
            Kind::Mul(c, factors) => {
                let mut result = Sparse::constant(Expr::num(c.clone()));
                for (base, exponent) in factors {
                    let expanded = self.power(base, exponent, depth + 1)?;
                    result = result.mul(&expanded, &self.cutoff, self.x)?;
                }
                // Mul._eval_nseries keeps O(x**n) when rebuilding a regular
                // nonpolynomial product changes its structure, even if the
                // rebuilt finite sum is mathematically exact. Add propagates
                // this marker; unchanged monomials and polynomials stay exact.
                if result.precision.is_none() && !result.terms.is_empty()
                    && !polynomial(e, self.x)
                    && result.terms.iter().all(|(_, c)| !depends(c, self.x))
                {
                    let rebuilt = add(&result.terms.iter().map(|(a, c)| mul2(c, &pow(self.x, a))).collect::<Vec<_>>());
                    if rebuilt != *e { result.precision = Some(self.cutoff.clone()); }
                }
                result
            }
            Kind::Pow(base, exponent) => self.power(base, exponent, depth + 1)?,
            Kind::Func(Func::Log, args) if args.len() == 1 => self.logarithm(self.expand(&args[0], depth + 1)?)?,
            _ => return None,
        };
        result.clip(&self.cutoff, self.x)
    }

    fn power(&self, base: &Expr, exponent: &Expr, depth: usize) -> Option<Sparse> {
        // A factor independent of the expansion variable remains its exact
        // coefficient, with the same domain exclusions as in the input.
        // For example y**-1 needs no new assertion that y is nonzero.
        if !depends(base, self.x) && !depends(exponent, self.x) {
            return Some(Sparse::constant(pow(base, exponent)));
        }
        if !real_exponent(exponent, self.x) { return None; }
        if base == self.x { return Some(Sparse::monomial(exponent.clone(), Expr::one())); }
        let expanded = self.expand(base, depth + 1)?;
        if exponent.is_one() { return Some(expanded); }
        self.analytic_power(expanded, exponent)
    }
}

/// A failed certificate leaves the original series error/fallback untouched.
pub(super) fn series(e: &Expr, x: &Expr, n: i64) -> Option<Expr> {
    fn has_float(e: &Expr) -> bool {
        e.as_number().is_some_and(Number::is_float) || e.args().iter().any(has_float)
    }
    if !x.is_symbol() || n < 1 || n > 128 || e.has_host() || e.has_nonfinite()
        || has_float(e) || crate::decompose::commutativity(e) != Some(true) || super::contains_order(e) { return None; }
    let target = Expr::int(n);
    // Pole factors can lower precision. Retry the entire finite computation
    // until the tracked absolute error reaches the requested order.
    for extra in [4, 8, 16, 32, 64] {
        let context = Expansion { x, cutoff: Expr::int(n + extra) };
        let Some(result) = context.expand(e, 0) else { continue };
        if let Some(p) = &result.precision {
            if compare(p, &target)? == Ordering::Less { continue; }
        }
        let result = result.clip(&target, x)?;
        let mut terms: Vec<Expr> = result.terms.into_iter().map(|(a, c)| mul2(&c, &pow(x, &a))).collect();
        if result.precision.is_some() { terms.push(super::order(&pow(x, &target), &[(x.clone(), Expr::zero())])); }
        return Some(add(&terms));
    }
    None
}
