//! Certified real-root analysis over QQ/ZZ. No floating-point root guesses.
//! Sparse storage and checked arithmetic bound work before dense allocations
//! or unbounded coefficient growth; isolating intervals are rational.
use super::fraction::Polynomial;
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::polyobj::{Domain, Poly};
use super::pseudo::{Budget, PseudoCoeff};
use std::collections::VecDeque;

type P = MPoly<Rat>;
const MAX_DEGREE: i64 = 4096;
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RootError {
    Domain(String),
    Invalid(String),
    Refinement(String),
    Limit(String),
}
pub type RootResult<T> = Result<T, RootError>;
impl From<super::convert::PolyError> for RootError {
    fn from(error: super::convert::PolyError) -> Self {
        let text = error.to_string();
        if text.contains("resource limit") {
            Self::Limit(text)
        } else {
            Self::Invalid(text)
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Bound {
    NegativeInfinity,
    Rational(Rat),
    PositiveInfinity,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RootInterval {
    pub lower: Rat,
    pub upper: Rat,
    pub multiplicity: u32,
}
#[derive(Clone, Debug, Default)]
pub struct IsolationOptions {
    pub epsilon: Option<Rat>,
    pub lower: Option<Rat>,
    pub upper: Option<Rat>,
    pub fast: bool,
}
#[derive(Clone, Debug, Default)]
pub struct RefinementOptions {
    pub epsilon: Option<Rat>,
    pub steps: Option<u64>,
    pub fast: bool,
    pub check_squarefree: bool,
}
struct Counter {
    polynomial: P,
    chain: Vec<P>,
}
struct Engine {
    budget: Budget,
    retained: u64,
}
impl Engine {
    fn new() -> Self {
        Self {
            budget: Budget::new(),
            retained: 0,
        }
    }
    fn charge(&mut self, a: &Rat, b: &Rat) -> RootResult<()> {
        self.budget
            .charge(1 + (a.n.bits() + a.d.bits() + b.n.bits() + b.d.bits()) / 64)?;
        Ok(())
    }
    fn add(&mut self, a: &Rat, b: &Rat) -> RootResult<Rat> {
        self.charge(a, b)?;
        Ok(a.checked_add(b)?)
    }
    fn sub(&mut self, a: &Rat, b: &Rat) -> RootResult<Rat> {
        self.add(a, &b.neg())
    }
    fn mul(&mut self, a: &Rat, b: &Rat) -> RootResult<Rat> {
        self.charge(a, b)?;
        Ok(a.checked_mul(b)?)
    }
    fn div(&mut self, a: &Rat, b: &Rat) -> RootResult<Rat> {
        self.charge(a, b)?;
        Ok(a.checked_div(b)?)
    }
    fn pow(&mut self, a: &Rat, mut n: u32) -> RootResult<Rat> {
        let mut out = Rat::one();
        let mut base = a.clone();
        while n != 0 {
            if n & 1 != 0 {
                out = self.mul(&out, &base)?;
            }
            n >>= 1;
            if n != 0 {
                base = self.mul(&base, &base)?;
            }
        }
        Ok(out)
    }
    fn eval(&mut self, p: &P, x: &Rat) -> RootResult<Rat> {
        let mut out = Rat::zero();
        let mut previous = p.degree(0).max(0) as u32;
        for (e, c) in p.iter() {
            let power = self.pow(x, previous - e[0])?;
            out = self.mul(&out, &power)?;
            out = self.add(&out, c)?;
            previous = e[0];
        }
        let power = self.pow(x, previous)?;
        self.mul(&out, &power)
    }
    fn derivative(&mut self, p: &P) -> RootResult<P> {
        let mut terms = Vec::with_capacity(p.len());
        for (e, c) in p.iter() {
            if e[0] != 0 {
                let value = self.mul(c, &Rat::from_i64(e[0] as i64))?;
                terms.push((vec![e[0] - 1], value));
            }
        }
        let result = P::from_terms(1, terms);
        self.budget.check(&result)?;
        Ok(result)
    }
    fn remainder(&mut self, a: &P, b: &P) -> RootResult<P> {
        if b.is_zero() {
            return Err(RootError::Invalid(
                "division by zero in real-root analysis".into(),
            ));
        }
        let mut r = a.clone();
        while !r.is_zero() && r.degree(0) >= b.degree(0) {
            self.budget.charge(1)?;
            let degree = r.degree(0);
            let coefficient = self.div(&r.lc(), &b.lc())?;
            let term = P::monomial(&[(degree - b.degree(0)) as u32], coefficient);
            let product = self.budget.mul(&term, b)?;
            r = self.budget.merge(&r, &product, true)?;
            if !r.is_zero() && r.degree(0) >= degree {
                return Err(RootError::Invalid(
                    "real-root division failed to lower its degree".into(),
                ));
            }
        }
        Ok(r)
    }
    fn quotient(&mut self, a: &P, b: &P) -> RootResult<P> {
        Ok(self.budget.exact_div(a, b)?)
    }
    fn monic(&mut self, a: &P) -> RootResult<P> {
        if a.is_zero() {
            return Ok(a.clone());
        }
        self.quotient(a, &P::constant(a.lc(), 1))
    }
    fn gcd(&mut self, a: &P, b: &P) -> RootResult<P> {
        let (mut a, mut b) = (a.clone(), b.clone());
        while !b.is_zero() {
            let r = self.remainder(&a, &b)?;
            a = b;
            b = self.monic(&r)?;
        }
        self.monic(&a)
    }
    fn counter(&mut self, f: &P) -> RootResult<Counter> {
        let derivative = self.derivative(f)?;
        let gcd = self.gcd(f, &derivative)?;
        let squarefree = if f.is_zero() {
            f.clone()
        } else {
            let q = self.quotient(f, &gcd)?;
            self.monic(&q)?
        };
        let derivative = self.derivative(&squarefree)?;
        let mut chain = vec![squarefree.clone()];
        self.budget.retain(&squarefree, &mut self.retained)?;
        if !derivative.is_zero() {
            self.budget.retain(&derivative, &mut self.retained)?;
            chain.push(derivative);
        }
        while chain.len() > 1 {
            let n = chain.len();
            let remainder = self.remainder(&chain[n - 2], &chain[n - 1])?.neg();
            if remainder.is_zero() {
                break;
            }
            self.budget.retain(&remainder, &mut self.retained)?;
            chain.push(remainder);
        }
        Ok(Counter {
            polynomial: squarefree,
            chain,
        })
    }
    fn variations(&mut self, counter: &Counter, point: &Bound) -> RootResult<i64> {
        let (mut changes, mut previous) = (0, 0);
        for polynomial in &counter.chain {
            let sign = match point {
                Bound::Rational(x) => {
                    let value = self.eval(polynomial, x)?;
                    if value.is_zero() {
                        0
                    } else if value.is_negative() {
                        -1
                    } else {
                        1
                    }
                }
                Bound::NegativeInfinity | Bound::PositiveInfinity => {
                    if polynomial.is_zero() {
                        0
                    } else {
                        let sign = if polynomial.lc().is_negative() { -1 } else { 1 };
                        if matches!(point, Bound::NegativeInfinity) && polynomial.degree(0) % 2 == 1
                        {
                            -sign
                        } else {
                            sign
                        }
                    }
                }
            };
            if sign != 0 {
                if previous != 0 && previous != sign {
                    changes += 1;
                }
                previous = sign;
            }
        }
        Ok(changes)
    }
    fn count(&mut self, counter: &Counter, lower: &Bound, upper: &Bound) -> RootResult<usize> {
        if lower > upper {
            return Err(RootError::Invalid(
                "real-root bounds must be in increasing order".into(),
            ));
        }
        if counter.polynomial.degree(0) <= 0 {
            return Ok(0);
        }
        let mut count = self.variations(counter, lower)? - self.variations(counter, upper)?;
        if let Bound::Rational(x) = lower {
            if self.eval(&counter.polynomial, x)?.is_zero() {
                count += 1;
            }
        }
        usize::try_from(count).map_err(|_| RootError::Invalid("negative Sturm root count".into()))
    }
    fn open_count(&mut self, counter: &Counter, a: &Rat, b: &Rat) -> RootResult<usize> {
        let mut count = self.variations(counter, &Bound::Rational(a.clone()))?
            - self.variations(counter, &Bound::Rational(b.clone()))?;
        if self.eval(&counter.polynomial, b)?.is_zero() {
            count -= 1;
        }
        usize::try_from(count)
            .map_err(|_| RootError::Invalid("negative open Sturm root count".into()))
    }
    fn midpoint(&mut self, a: &Rat, b: &Rat) -> RootResult<Rat> {
        let sum = self.add(a, b)?;
        self.div(&sum, &Rat::from_i64(2))
    }
    // The least-denominator rational in the closed interval, obtained by a
    // continued-fraction descent. It detects rational singleton roots without
    // integer factorization or enumerating divisors of enormous coefficients.
    fn simplest(&mut self, lower: &Rat, upper: &Rat) -> RootResult<Rat> {
        let (mut a, mut b) = (lower.clone(), upper.clone());
        let mut prefix = Vec::new();
        let mut out = loop {
            self.budget.charge(1)?;
            let floor = Rat::from(a.floor());
            if a == floor {
                break a;
            }
            if floor != Rat::from(b.floor()) {
                break self.add(&floor, &Rat::one())?;
            }
            prefix.push(floor.clone());
            let low = self.sub(&a, &floor)?;
            let high = self.sub(&b, &floor)?;
            a = self.div(&Rat::one(), &high)?;
            b = self.div(&Rat::one(), &low)?;
        };
        for floor in prefix.into_iter().rev() {
            let inverse = self.div(&Rat::one(), &out)?;
            out = self.add(&floor, &inverse)?;
        }
        Ok(out)
    }
    fn radius(&mut self, f: &P) -> RootResult<Rat> {
        let mut maximum = Rat::zero();
        let leading = f.lc();
        for (e, c) in f.iter() {
            if e[0] as i64 != f.degree(0) {
                let ratio = self.div(&c.abs(), &leading.abs())?;
                if ratio > maximum {
                    maximum = ratio;
                }
            }
        }
        self.add(&Rat::from(maximum.floor()), &Rat::from_i64(2))
    }
    fn rational_threshold(&mut self, f: &P) -> RootResult<Rat> {
        let mut denominator = Int::one();
        for coefficient in f.coeffs() {
            self.budget.charge(1 + coefficient.d.bits() / 64)?;
            let gcd = denominator.gcd(&coefficient.d);
            denominator = denominator.checked_div(&gcd)?.checked_mul(&coefficient.d)?;
        }
        let leading = f.lc();
        let bound = denominator
            .checked_div(&leading.d)?
            .checked_mul(&leading.n)?
            .abs();
        let denominator = bound.checked_mul(&bound)?.checked_mul(&Int::from(2))?;
        Ok(Rat::new(Int::one(), denominator))
    }
    fn multiplicity_factors(&mut self, f: &P) -> RootResult<Vec<(Counter, u32)>> {
        if f.degree(0) <= 0 {
            return Ok(Vec::new());
        }
        let derivative = self.derivative(f)?;
        let mut common = self.gcd(f, &derivative)?;
        let mut remaining = self.quotient(f, &common)?;
        remaining = self.monic(&remaining)?;
        let mut factors = Vec::new();
        let mut multiplicity = 1u32;
        while remaining.degree(0) > 0 {
            let overlap = self.gcd(&remaining, &common)?;
            let factor = self.quotient(&remaining, &overlap)?;
            if factor.degree(0) > 0 {
                factors.push((self.counter(&factor)?, multiplicity));
            }
            remaining = overlap;
            common = self.quotient(&common, &remaining)?;
            multiplicity += 1;
            self.budget.charge(1)?;
        }
        Ok(factors)
    }
    fn split(&mut self, counter: &Counter, a: Rat, b: Rat) -> RootResult<(Rat, Rat)> {
        let middle = self.midpoint(&a, &b)?;
        if self.eval(&counter.polynomial, &middle)?.is_zero() {
            return Ok((middle.clone(), middle));
        }
        if self.open_count(counter, &a, &middle)? == 1 {
            Ok((a, middle))
        } else {
            Ok((middle, b))
        }
    }
}
fn validate(p: &Polynomial) -> RootResult<Poly> {
    if p.field.is_some()
        || p.algebraic.is_some()
        || p.modulus.is_some()
        || p.domain == Domain::RR
        || p.ngens != 1
        || p.base.gens.len() != 1
    {
        return Err(RootError::Domain(
            "certified real-root analysis requires a univariate polynomial over QQ or ZZ".into(),
        ));
    }
    if p.rep.degree(0) > MAX_DEGREE {
        return Err(RootError::Limit(format!(
            "real-root resource limit: degree exceeds {MAX_DEGREE}"
        )));
    }
    Budget::new().check(&p.rep)?;
    Ok(p.base.clone())
}
fn epsilon(epsilon: &Option<Rat>) -> RootResult<()> {
    if epsilon.as_ref().is_some_and(|e| e <= &Rat::zero()) {
        return Err(RootError::Invalid("eps must be a positive rational".into()));
    }
    if let Some(e) = epsilon {
        e.validate_size()?;
    }
    Ok(())
}

pub fn sturm(p: &Polynomial, auto: bool) -> RootResult<Vec<Polynomial>> {
    let mut template = validate(p)?;
    if template.domain == Domain::ZZ && !auto {
        return Err(RootError::Domain(
            "cannot compute a Sturm sequence over ZZ with auto=False".into(),
        ));
    }
    template.domain = Domain::QQ;
    let mut engine = Engine::new();
    let counter = engine.counter(&template.rep)?;
    Ok(counter
        .chain
        .into_iter()
        .map(|r| template.with_rep(r).into())
        .collect())
}
pub fn count_roots(p: &Polynomial, lower: Bound, upper: Bound) -> RootResult<usize> {
    for bound in [&lower, &upper] {
        if let Bound::Rational(value) = bound {
            value.validate_size()?;
        }
    }
    let template = validate(p)?;
    let mut engine = Engine::new();
    let counter = engine.counter(&template.rep)?;
    engine.count(&counter, &lower, &upper)
}
pub fn intervals(p: &Polynomial, options: &IsolationOptions) -> RootResult<Vec<RootInterval>> {
    let template = validate(p)?;
    epsilon(&options.epsilon)?;
    for value in [&options.lower, &options.upper].into_iter().flatten() {
        value.validate_size()?;
    }
    if let (Some(a), Some(b)) = (&options.lower, &options.upper) {
        if a > b {
            return Err(RootError::Invalid(
                "real-root bounds must be in increasing order".into(),
            ));
        }
    }
    let mut engine = Engine::new();
    let counter = engine.counter(&template.rep)?;
    if counter.polynomial.degree(0) <= 0 {
        return Ok(Vec::new());
    }
    let factors = engine.multiplicity_factors(&template.rep)?;
    let threshold = engine.rational_threshold(&template.rep)?;
    let radius = engine.radius(&counter.polynomial)?;
    let a = options
        .lower
        .as_ref()
        .map_or_else(|| radius.neg(), |v| v.max(&radius.neg()).clone());
    let b = options
        .upper
        .as_ref()
        .map_or_else(|| radius.clone(), |v| v.min(&radius).clone());
    if a > b {
        return Ok(Vec::new());
    }
    let expected = engine.count(
        &counter,
        &Bound::Rational(a.clone()),
        &Bound::Rational(b.clone()),
    )?;
    let mut roots = Vec::<(Rat, Rat)>::new();
    if engine.eval(&counter.polynomial, &a)?.is_zero() {
        roots.push((a.clone(), a.clone()));
    }
    if b != a && engine.eval(&counter.polynomial, &b)?.is_zero() {
        roots.push((b.clone(), b.clone()));
    }
    let mut queue = VecDeque::new();
    if a < b {
        let count = engine.open_count(&counter, &a, &b)?;
        if count != 0 {
            queue.push_back((a, b, count));
        }
    }
    while let Some((a, b, count)) = queue.pop_front() {
        engine.budget.charge(1)?;
        if count == 1 {
            let candidate = engine.simplest(&a, &b)?;
            if candidate > a
                && candidate < b
                && engine.eval(&counter.polynomial, &candidate)?.is_zero()
            {
                roots.push((candidate.clone(), candidate));
                continue;
            }
            let width = engine.sub(&b, &a)?;
            let ready = width < threshold && options.epsilon.as_ref().is_none_or(|e| &width < e);
            if ready
                && !engine.eval(&counter.polynomial, &a)?.is_zero()
                && !engine.eval(&counter.polynomial, &b)?.is_zero()
            {
                roots.push((a, b));
                continue;
            }
        }
        let midpoint = engine.midpoint(&a, &b)?;
        let at_middle = engine.eval(&counter.polynomial, &midpoint)?.is_zero();
        if at_middle {
            roots.push((midpoint.clone(), midpoint.clone()));
        }
        let left = engine.open_count(&counter, &a, &midpoint)?;
        let right = engine.open_count(&counter, &midpoint, &b)?;
        if left + right + usize::from(at_middle) != count {
            return Err(RootError::Invalid(
                "Sturm subdivision did not preserve its root count".into(),
            ));
        }
        if left != 0 {
            queue.push_back((a, midpoint.clone(), left));
        }
        if right != 0 {
            queue.push_back((midpoint, b, right));
        }
        // fast=True explores newly isolated children immediately; normal mode
        // uses breadth-first subdivision. Both have the same exact certificate.
        if options.fast && queue.len() > 1 {
            if let Some(last) = queue.pop_back() {
                queue.push_front(last);
            }
        }
    }
    roots.sort();
    roots.dedup();
    if roots.len() != expected {
        return Err(RootError::Invalid(
            "incomplete real-root isolation certificate".into(),
        ));
    }
    roots
        .into_iter()
        .map(|(lower, upper)| {
            let mut multiplicity = None;
            for (factor, power) in &factors {
                if engine.count(
                    factor,
                    &Bound::Rational(lower.clone()),
                    &Bound::Rational(upper.clone()),
                )? == 1
                {
                    if multiplicity.replace(*power).is_some() {
                        return Err(RootError::Invalid(
                            "overlapping multiplicity factors".into(),
                        ));
                    }
                }
            }
            Ok(RootInterval {
                lower,
                upper,
                multiplicity: multiplicity.ok_or_else(|| {
                    RootError::Invalid("missing root multiplicity certificate".into())
                })?,
            })
        })
        .collect()
}
pub fn refine_root(
    p: &Polynomial,
    lower: Rat,
    upper: Rat,
    options: &RefinementOptions,
) -> RootResult<(Rat, Rat)> {
    lower.validate_size()?;
    upper.validate_size()?;
    let template = validate(p)?;
    epsilon(&options.epsilon)?;
    if lower > upper {
        return Err(RootError::Invalid(
            "real-root bounds must be in increasing order".into(),
        ));
    }
    let mut engine = Engine::new();
    let counter = engine.counter(&template.rep)?;
    if options.check_squarefree && counter.polynomial.degree(0) != template.rep.degree(0) {
        return Err(RootError::Refinement(
            "polynomial must be square-free".into(),
        ));
    }
    if engine.count(
        &counter,
        &Bound::Rational(lower.clone()),
        &Bound::Rational(upper.clone()),
    )? != 1
    {
        return Err(RootError::Refinement(
            "there must be exactly one real root in the interval".into(),
        ));
    }
    if engine.eval(&counter.polynomial, &lower)?.is_zero() {
        return Ok((lower.clone(), lower));
    }
    if engine.eval(&counter.polynomial, &upper)?.is_zero() {
        return Ok((upper.clone(), upper));
    }
    let (mut a, mut b) = (lower, upper);
    let mut steps = options.steps.unwrap_or(if options.epsilon.is_some() {
        u64::MAX
    } else {
        1
    });
    while a < b && steps != 0 {
        let width = engine.sub(&b, &a)?;
        if options.epsilon.as_ref().is_some_and(|e| &width < e) {
            break;
        }
        (a, b) = engine.split(&counter, a, b)?;
        if options.fast && a < b {
            (a, b) = engine.split(&counter, a, b)?;
        }
        steps -= 1;
    }
    Ok((a, b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{parser::parse_str, Expr};
    fn p(source: &str) -> Polynomial {
        Poly::from_expr(&parse_str(source).unwrap(), &[Expr::symbol("x")])
            .unwrap()
            .into()
    }
    fn r(n: i64, d: i64) -> Rat {
        Rat::new(Int::from(n), Int::from(d))
    }
    fn b(n: i64) -> Bound {
        Bound::Rational(Rat::from_i64(n))
    }
    #[test]
    fn exact_sturm_sequence_and_closed_endpoint_counts() {
        let f = p("(x-1)^3*(x+2)^2");
        let sequence = sturm(&f, true).unwrap();
        assert_eq!(
            sequence.iter().map(|p| p.rep.clone()).collect::<Vec<_>>(),
            vec![
                p("x^2+x-2").rep.clone(),
                p("2*x+1").rep.clone(),
                p("9/4").rep.clone()
            ]
        );
        assert!(sequence.iter().all(|p| p.domain == Domain::QQ));
        assert!(matches!(sturm(&f, false), Err(RootError::Domain(_))));
        assert_eq!(
            count_roots(&f, Bound::NegativeInfinity, Bound::PositiveInfinity).unwrap(),
            2
        );
        assert_eq!(count_roots(&f, b(-2), b(1)).unwrap(), 2);
        assert_eq!(count_roots(&f, b(1), b(1)).unwrap(), 1);
        assert_eq!(count_roots(&f, b(-1), b(0)).unwrap(), 0);
        assert_eq!(count_roots(&p("x^5-x-1"), b(1), b(2)).unwrap(), 1);
        assert_eq!(
            count_roots(
                &p("x^4+1"),
                Bound::NegativeInfinity,
                Bound::PositiveInfinity
            )
            .unwrap(),
            0
        );
        assert!(sturm(&p("0"), true).unwrap()[0].is_zero());
        assert!(sturm(&p("7"), true).unwrap()[0].rep.is_one());
    }
    #[test]
    fn isolation_certifies_completeness_multiplicity_and_precisions() {
        let f = p("(x-1)^3*(x+2)^2*(x^2-2)");
        for fast in [false, true] {
            let roots = intervals(
                &f,
                &IsolationOptions {
                    epsilon: Some(r(1, 10000)),
                    fast,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(roots.len(), 4);
            assert_eq!(
                roots
                    .iter()
                    .find(|v| v.lower == r(-2, 1) && v.upper == r(-2, 1))
                    .unwrap()
                    .multiplicity,
                2
            );
            assert_eq!(
                roots
                    .iter()
                    .find(|v| v.lower == r(1, 1) && v.upper == r(1, 1))
                    .unwrap()
                    .multiplicity,
                3
            );
            assert!(roots.windows(2).all(|pair| pair[0].upper < pair[1].lower));
            for root in roots {
                assert_eq!(
                    count_roots(
                        &f,
                        Bound::Rational(root.lower.clone()),
                        Bound::Rational(root.upper.clone())
                    )
                    .unwrap(),
                    1
                );
                assert!(root.upper.sub(&root.lower) < r(1, 10000));
            }
        }
        let roots = intervals(
            &p("(3*x-1)*(3*x+1)*(1000000*x-1)*(1000000*x-2)"),
            &Default::default(),
        )
        .unwrap();
        assert_eq!(
            roots.iter().map(|v| v.lower.clone()).collect::<Vec<_>>(),
            vec![r(-1, 3), r(1, 1000000), r(1, 500000), r(1, 3)]
        );
        assert!(roots
            .iter()
            .all(|v| v.lower == v.upper && v.multiplicity == 1));
        let clipped = intervals(
            &f,
            &IsolationOptions {
                lower: Some(r(-2, 1)),
                upper: Some(r(1, 1)),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(clipped.len(), 3);
        assert_eq!(
            intervals(
                &f,
                &IsolationOptions {
                    lower: Some(r(1, 1)),
                    upper: Some(r(1, 1)),
                    ..Default::default()
                }
            )
            .unwrap()
            .len(),
            1
        );
    }
    #[test]
    fn refinement_keeps_the_same_unique_root_and_bounds_steps() {
        let f = p("x^2-2");
        assert_eq!(
            refine_root(
                &f,
                r(1, 1),
                r(2, 1),
                &RefinementOptions {
                    steps: Some(0),
                    ..Default::default()
                }
            )
            .unwrap(),
            (r(1, 1), r(2, 1))
        );
        assert_eq!(
            refine_root(&f, r(1, 1), r(2, 1), &Default::default()).unwrap(),
            (r(1, 1), r(3, 2))
        );
        let tiny = Rat::new(
            Int::one(),
            Int::from_big(num_bigint::BigInt::from(1u32) << 200usize),
        );
        for fast in [false, true] {
            let (a, b) = refine_root(
                &f,
                r(1, 1),
                r(2, 1),
                &RefinementOptions {
                    epsilon: Some(tiny.clone()),
                    fast,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(a >= r(1, 1) && b <= r(2, 1) && b.sub(&a) < tiny);
            assert_eq!(
                count_roots(&f, Bound::Rational(a), Bound::Rational(b)).unwrap(),
                1
            );
        }
        let limited = refine_root(
            &f,
            r(1, 1),
            r(2, 1),
            &RefinementOptions {
                epsilon: Some(tiny),
                steps: Some(1),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(limited, (r(1, 1), r(3, 2)));
        assert_eq!(
            refine_root(&p("3*x-1"), r(1, 3), r(1, 3), &Default::default()).unwrap(),
            (r(1, 3), r(1, 3))
        );
    }
    #[test]
    fn invalid_inputs_and_resource_bounds_fail_without_guessing() {
        assert!(matches!(
            count_roots(&p("x"), b(1), b(0)),
            Err(RootError::Invalid(_))
        ));
        assert!(matches!(
            refine_root(&p("x^2-2"), r(1, 1), r(1, 1), &Default::default()),
            Err(RootError::Refinement(_))
        ));
        assert!(matches!(
            refine_root(&p("x^2-2"), r(-2, 1), r(2, 1), &Default::default()),
            Err(RootError::Refinement(_))
        ));
        assert!(matches!(
            intervals(
                &p("x"),
                &IsolationOptions {
                    epsilon: Some(Rat::zero()),
                    ..Default::default()
                }
            ),
            Err(RootError::Invalid(_))
        ));
        assert!(matches!(
            refine_root(
                &p("(x-1)^2"),
                r(0, 1),
                r(2, 1),
                &RefinementOptions {
                    check_squarefree: true,
                    ..Default::default()
                }
            ),
            Err(RootError::Refinement(_))
        ));
        assert!(matches!(
            count_roots(&p("x^4097+1"), b(-2), b(2)),
            Err(RootError::Limit(_))
        ));
        let finite: Polynomial = p("x^2+1").base.set_modulus(5).into();
        assert!(matches!(sturm(&finite, true), Err(RootError::Domain(_))));
    }
}
