//! Real-branch classification using outward rational enclosures. Approximate
//! floating point comparisons at -1/e must not drop a real Lambert branch.
use crate::expr::{Constant, Expr, Func, Kind};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, ToPrimitive, Zero};

fn rational(n: &crate::number::Number) -> Option<Q> {
    if let Some(q) = n.to_ratio() {
        return Some(q);
    }
    if matches!(n, crate::number::Number::Float(f) if !f.is_finite()) {
        return None;
    }
    let value = n.to_bigfloat(n.float_precision()?);
    if value.exp.unsigned_abs() > 1_000_000 {
        return None;
    }
    let exp = value.exp.small()?;
    Some(if exp >= 0 {
        Q::from_integer(value.man << exp as usize)
    } else {
        Q::new(
            value.man,
            BigInt::one() << exp.unsigned_abs() as usize,
        )
    })
}

#[derive(Clone)]
struct Interval {
    lo: Q,
    hi: Q,
}
impl Interval {
    fn point(q: Q) -> Self {
        Self {
            lo: q.clone(),
            hi: q,
        }
    }
    fn add(&self, b: &Self) -> Self {
        Self {
            lo: &self.lo + &b.lo,
            hi: &self.hi + &b.hi,
        }
    }
    fn neg(&self) -> Self {
        Self {
            lo: -&self.hi,
            hi: -&self.lo,
        }
    }
    fn mul(&self, b: &Self) -> Self {
        let products = [
            &self.lo * &b.lo,
            &self.lo * &b.hi,
            &self.hi * &b.lo,
            &self.hi * &b.hi,
        ];
        Self {
            lo: products.iter().min().unwrap().clone(),
            hi: products.iter().max().unwrap().clone(),
        }
    }
    fn inv(&self) -> Option<Self> {
        if self.lo <= Q::zero() && self.hi >= Q::zero() {
            return None;
        }
        Some(Self {
            lo: self.hi.recip(),
            hi: self.lo.recip(),
        })
    }
    fn round(&self, bits: u32) -> Self {
        let scale = BigInt::one() << bits;
        let lower = (self.lo.numer() * &scale).div_floor(self.lo.denom());
        let upper = (self.hi.numer() * &scale).div_ceil(self.hi.denom());
        Self {
            lo: Q::new(lower, scale.clone()),
            hi: Q::new(upper, scale),
        }
    }
    fn powi(&self, exponent: i64, bits: u32) -> Option<Self> {
        if exponent.unsigned_abs() > 4096 {
            return None;
        }
        let mut n = exponent.unsigned_abs();
        let mut base = if exponent < 0 {
            self.inv()?
        } else {
            self.clone()
        };
        let mut result = Self::point(Q::one());
        while n > 0 {
            if n & 1 != 0 {
                result = result.mul(&base).round(bits);
            }
            n >>= 1;
            if n > 0 {
                base = base.mul(&base).round(bits);
            }
        }
        Some(result)
    }
    fn sqrt(&self, bits: u32) -> Option<Self> {
        self.nth_root(2, bits)
    }
    fn nth_root(&self, degree: u32, bits: u32) -> Option<Self> {
        if self.lo < Q::zero() {
            return None;
        }
        if degree == 0 || degree > 4096 || bits.checked_mul(degree)? > 1_000_000 {
            return None;
        }
        let scale = BigInt::one() << bits;
        let lower = ((self.lo.numer() << (degree * bits)) / self.lo.denom()).nth_root(degree);
        let upper = ((self.hi.numer() << (degree * bits)) / self.hi.denom()).nth_root(degree) + 1;
        Some(Self {
            lo: Q::new(lower, scale.clone()),
            hi: Q::new(upper, scale),
        })
    }
}

fn exp(q: &Q, bits: u32) -> Option<Interval> {
    if q.is_negative() {
        return exp(&(-q), bits)?.inv().map(|r| r.round(bits));
    }
    if *q > Q::from_integer(1_000_000.into()) {
        return None;
    }
    let mut a = q.clone();
    let mut squarings = 0;
    while a > Q::one() {
        a /= BigInt::from(2);
        squarings += 1;
    }
    let working = bits + 32 + 2 * squarings;
    let epsilon = Q::new(BigInt::one(), BigInt::one() << working);
    let mut sum = Q::one();
    let mut term = Q::one();
    for n in 1..=(working + 32) {
        term = term * &a / BigInt::from(n);
        sum += &term;
        let next = &term * &a / BigInt::from(n + 1);
        if next <= epsilon {
            let mut result = Interval {
                lo: sum.clone(),
                hi: sum + next * BigInt::from(2),
            }
            .round(working);
            for _ in 0..squarings {
                result = result.mul(&result).round(working);
            }
            return Some(result.round(bits));
        }
    }
    None
}

/// Outward enclosure of atanh(t) for a rational |t| <= 1/3 from the series
/// sum t^(2n+1)/(2n+1); after the last term the remainder is at most
/// |t|^(2n+1)/((2n+1)(1 - t^2)). Partial sums are rounded outward.
fn atanh_small(t: &Q, bits: u32) -> Option<Interval> {
    let working = bits + 32;
    let square = t * t;
    let complement = Q::one() - &square;
    let threshold = Q::new(BigInt::one(), BigInt::one() << (bits + 8));
    let mut power = Interval::point(t.clone());
    let mut sum = Interval::point(Q::zero());
    for n in 0u32..1_000_000 {
        let divisor = Interval::point(Q::new(BigInt::one(), BigInt::from(2 * n + 1)));
        sum = sum.add(&power.mul(&divisor)).round(working);
        power = power.mul(&Interval::point(square.clone())).round(working);
        let magnitude = power.lo.abs().max(power.hi.abs());
        let tail = magnitude / (Q::from_integer(BigInt::from(2 * n + 3)) * &complement);
        if tail < threshold {
            return Some(Interval { lo: &sum.lo - &tail, hi: &sum.hi + &tail }.round(bits));
        }
    }
    None
}

/// Outward enclosure of log(q) for a positive rational, as k*log(2) plus
/// 2*atanh((m - 1)/(m + 1)) with q = 2^k*m and m in [2/3, 4/3].
fn log_point(q: &Q, bits: u32) -> Option<Interval> {
    if *q <= Q::zero() {
        return None;
    }
    let working = bits + 16;
    let shift = q.numer().bits() as i64 - q.denom().bits() as i64;
    let mut k = shift;
    let mut m = if shift >= 0 {
        q / Q::from_integer(BigInt::one() << shift as usize)
    } else {
        q * Q::from_integer(BigInt::one() << shift.unsigned_abs() as usize)
    };
    let (upper, lower) = (Q::new(4.into(), 3.into()), Q::new(2.into(), 3.into()));
    while m > upper {
        m /= Q::from_integer(2.into());
        k += 1;
    }
    while m < lower {
        m *= Q::from_integer(2.into());
        k -= 1;
    }
    let two = Interval::point(Q::from_integer(2.into()));
    let t = (&m - Q::one()) / (&m + Q::one());
    let mut result = atanh_small(&t, working)?.mul(&two);
    if k != 0 {
        let log2 = atanh_small(&Q::new(BigInt::one(), 3.into()), working + 64)?.mul(&two);
        result = result.add(&log2.mul(&Interval::point(Q::from_integer(k.into()))));
    }
    Some(result.round(bits))
}

/// Binary splitting of the exact arctan hypergeometric sum. Return P,Q,T
/// with product of term ratios P/Q and partial sum T/Q. Integer products
/// are balanced instead of performing a large division for every term.
fn atan_split(denominator: u32, start: u32, end: u32) -> (BigInt, BigInt, BigInt) {
    if end - start == 1 {
        if start == 0 {
            return (BigInt::one(), denominator.into(), BigInt::one());
        }
        let p = -BigInt::from(2 * start - 1);
        let q = BigInt::from(2 * start + 1) * denominator * denominator;
        return (p.clone(), q, p);
    }
    let middle = start + (end - start) / 2;
    let (pa, qa, ta) = atan_split(denominator, start, middle);
    let (pb, qb, tb) = atan_split(denominator, middle, end);
    let t = ta * &qb + &pa * tb;
    (pa * pb, qa * qb, t)
}

/// Alternating-series remainder bounds, followed by one outward rounding.
/// With L=floor(log2(n)), n^(2N+1) >= 2^((2N+1)L), so the selected N
/// makes the first omitted term smaller than 2^(-bits-16).
fn atan_reciprocal(denominator: u32, bits: u32) -> Interval {
    debug_assert!(denominator >= 2);
    let lower_log = denominator.ilog2();
    let terms = (bits + 16).div_ceil(2 * lower_log) + 1;
    let (_, q, t) = atan_split(denominator, 0, terms);
    let sum = Q::new(t, q);
    let tail = Q::new(BigInt::one(),
        BigInt::from(2 * terms + 1) * BigInt::from(denominator).pow(2 * terms + 1));
    let value = if terms % 2 == 0 {
        Interval { lo: sum.clone(), hi: sum + tail }
    } else {
        Interval { lo: &sum - tail, hi: sum }
    };
    value.round(bits)
}

fn pi(bits: u32) -> Interval {
    // Machin's exact identity: pi = 16 atan(1/5) - 4 atan(1/239).
    let a = atan_reciprocal(5, bits + 8);
    let b = atan_reciprocal(239, bits + 8);
    a.mul(&Interval::point(Q::from_integer(16.into())))
        .add(&b.mul(&Interval::point(Q::from_integer((-4).into())))).round(bits)
}

fn enclosure(e: &Expr, bits: u32, depth: u32) -> Option<Interval> {
    if depth > 64 {
        return None;
    }
    let recurse = |e: &Expr| enclosure(e, bits + 8, depth + 1);
    Some(match e.kind() {
        Kind::Num(n) => Interval::point(rational(n)?),
        Kind::Const(Constant::E) => exp(&Q::one(), bits)?,
        Kind::Const(Constant::Pi) => pi(bits),
        Kind::Add(c, terms) => {
            let mut result = Interval::point(rational(c)?);
            for (term, coefficient) in terms {
                result = result
                    .add(&recurse(term)?.mul(&Interval::point(rational(coefficient)?)))
                    .round(bits);
            }
            result
        }
        Kind::Mul(c, factors) => {
            let mut result = Interval::point(rational(c)?);
            for (base, exponent) in factors {
                result = result
                    .mul(&recurse(&crate::ops::mul_factor(base, exponent))?)
                    .round(bits);
            }
            result
        }
        Kind::Pow(base, exponent) => {
            let base = recurse(base)?;
            if let Some(k) = exponent.as_i64() {
                base.powi(k, bits)?
            } else if *exponent == Expr::rational(1, 2) {
                base.sqrt(bits)?
            } else if *exponent == Expr::rational(-1, 2) {
                base.sqrt(bits)?.inv()?
            } else {
                let rational = exponent.as_number()?.to_ratio()?;
                let power = rational.numer().to_i64()?;
                base.nth_root(rational.denom().to_u32()?, bits)?
                    .powi(power, bits)?
            }
        }
        Kind::Func(Func::Exp, args) => {
            let argument = recurse(&args[0])?;
            Interval {
                lo: exp(&argument.lo, bits)?.lo,
                hi: exp(&argument.hi, bits)?.hi,
            }
        }
        Kind::Func(Func::Log, args) if args.len() == 1 => {
            let argument = recurse(&args[0])?;
            if argument.lo <= Q::zero() {
                return None;
            }
            Interval {
                lo: log_point(&argument.lo, bits)?.lo,
                hi: log_point(&argument.hi, bits)?.hi,
            }
        }
        _ => return None,
    })
}

/// Return an integer only when both outward rational endpoints prove it.
/// mode -1/0/1 means floor/truncate/ceiling respectively. Failure to isolate
/// an integer cell leaves the caller's expression undecided.
pub(crate) fn integer_round(e: &Expr, mode: i8) -> Option<BigInt> {
    for bits in [32,64,128,256,512,1024,2048,4096,8192,16384,32768] {
        let value = enclosure(e, bits, 0)?;
        let round = |q: &Q| match mode {
            -1 => q.numer().div_floor(q.denom()),
            1 => q.numer().div_ceil(q.denom()),
            _ => q.to_integer(),
        };
        let lower = round(&value.lo);
        if lower == round(&value.hi) { return Some(lower); }
    }
    None
}

/// Sign certified by rational endpoints; unsupported expressions stay unknown.
pub(crate) fn real_sign(e:&Expr)->Option<i8>{
    for bits in [32,64,128,256,512,1024,2048,4096,8192,16384,32768]{
        let value=enclosure(e,bits,0)?;
        if value.lo>Q::zero(){return Some(1);}
        if value.hi<Q::zero(){return Some(-1);}
        if value.lo.is_zero()&&value.hi.is_zero(){return Some(0);}
    }
    None
}

/// W_0 is real on [-1/e,infinity); W_-1 on [-1/e,0). None is an
/// unproved domain, never evidence for selecting or rejecting a real branch.
pub(crate) fn real(z: &Expr, branch: i64) -> Option<bool> {
    let declared = crate::assumptions::exact_queries(|| crate::assumptions::is_real(z));
    if declared == Some(false) {
        return Some(false);
    }
    if !matches!(branch, 0 | -1) {
        return declared.filter(|r| *r).map(|_| false);
    }
    let nonnegative = crate::assumptions::exact_queries(|| crate::assumptions::is_nonnegative(z));
    if nonnegative == Some(true) {
        return Some(branch == 0);
    }
    if crate::ops::add2(z, &crate::functions::exp(&Expr::minus_one())).is_zero() {
        return Some(true);
    }
    for bits in [32, 64, 128, 256, 512, 1024, 2048, 4096] {
        let argument = enclosure(z, bits, 0)?;
        let boundary = exp(&Q::one(), bits)?.inv()?.neg();
        if argument.hi < boundary.lo {
            return Some(false);
        }
        if branch == -1 && argument.lo >= Q::zero() {
            return Some(false);
        }
        if argument.lo >= boundary.hi && (branch == 0 || argument.hi < Q::zero()) {
            return Some(true);
        }
    }
    None
}

/// Halley iteration entirely in native multiprecision arithmetic. Branch
/// selection and initial estimates do not round through a machine float.
pub(crate) fn mp(args: &[Expr], precision: u32) -> Option<crate::evalf_mp::MpComplex> {
    use crate::bigfloat::{pi, BigFloat as B};
    use crate::evalf_mp::{ev, MpComplex as C};
    let branch = args.get(1).map_or(Some(0), Expr::as_i64)?;
    let mut working = precision.checked_add(64)?;
    let mut z = ev(&args[0], working)?;
    if z.re.is_zero() && z.im.is_zero() {
        return (branch == 0).then(|| C::real(B::zero()));
    }
    let one = C::real(B::one());
    let two = C::real(B::from_i64(2));
    let mut near_branchpoint = None;
    if (branch == 0 || (branch == -1 && z.im.sign() >= 0)) && z.re.sign() < 0 {
        let q_expr = crate::ops::mul2(
            &Expr::int(2),
            &crate::ops::add2(&Expr::one(), &crate::ops::mul2(&Expr::e(), &args[0])),
        );
        let q = ev(&q_expr, working)?;
        if q.abs(working).is_zero() || q.abs(working).mag() < 1 {
            if q.re.is_zero() && q.im.is_zero() {
                return Some(C::real(B::from_i64(-1)));
            }
            let extra = q.abs(working).mag().saturating_neg().max(0) as u32;
            if extra > 4096 {
                return None;
            }
            working += extra;
            z = ev(&args[0], working)?;
            let delta = ev(&crate::ops::sqrt(&q_expr), working)?;
            near_branchpoint = Some(
                C::real(B::from_i64(-1))
                    .add(&if branch == -1 { delta.neg() } else { delta }, working),
            );
        }
    }
    let mut w = if let Some(seed) = near_branchpoint {
        seed
    } else if branch == 0 && z.re.sign() > 0 {
        z.add(&one, working).ln(working)?
    } else if branch == 0 && z.abs(working).mag() < -2 {
        z.clone()
    } else if branch == 0 && !z.is_real() && z.abs(working).mag() < 0 {
        z.mul(&one.sub(&z, working), working)
    } else if branch == 0 && z.is_real() && real(&args[0], 0) == Some(true) {
        C::real(B::from_f64(-0.5)) // exactly -1/2; no approximate seed digits
    } else if branch == -1 && z.is_real() && real(&args[0], -1) == Some(true) {
        let l = z.re.neg().ln(working)?;
        C::real(l.sub(&l.neg().ln(working)?, working))
    } else {
        let shifted = z.ln(working)?.add(
            &C {
                re: B::zero(),
                im: pi(working).mul(&B::from_i64(branch), working).ldexp(1),
            },
            working,
        );
        shifted.sub(&shifted.ln(working)?, working)
    };
    for _ in 0..128 {
        let exponential = w.exp(working)?;
        let residual = w.mul(&exponential, working).sub(&z, working);
        let wp1 = w.add(&one, working);
        let correction = w
            .add(&two, working)
            .mul(&residual, working)
            .div(&wp1.scale(&B::from_i64(2), working), working)?;
        let denominator = exponential.mul(&wp1, working).sub(&correction, working);
        let step = residual.div(&denominator, working)?;
        w = w.sub(&step, working);
        let threshold = w.abs(working).mag().max(0) - precision as i128 - 24;
        if step.abs(working).mag() < threshold {
            let residual = w
                .mul(&w.exp(working)?, working)
                .sub(&z, working)
                .abs(working);
            if residual.mag() < z.abs(working).mag().max(0) - precision as i128 - 16 {
                return Some(w.round_to(precision));
            }
        }
    }
    None
}

#[cfg(test)]
mod enclosure_tests {
    use super::*;

    #[test]
    fn certified_sign_resolves_tiny_exact_and_transcendental_differences(){
        let tiny=Expr::num(crate::number::Number::rational(BigInt::one(),BigInt::from(10u32).pow(100)).unwrap());
        assert_eq!(real_sign(&tiny),Some(1));
        assert_eq!(real_sign(&crate::ops::sub(&crate::functions::exp(&tiny),&Expr::one())),Some(1));
        assert_eq!(real_sign(&crate::ops::sub(&crate::functions::exp(&crate::ops::neg(&tiny)),&Expr::one())),Some(-1));
        assert_eq!(real_sign(&Expr::zero()),Some(0));
        assert_eq!(real_sign(&Expr::symbol("x")),None);
        assert_eq!(real_sign(&Expr::i()),None);
    }

    fn direct_sum(n: u32, terms: u32) -> Q {
        (0..terms).map(|k| Q::new(
            if k % 2 == 0 { BigInt::one() } else { -BigInt::one() },
            BigInt::from(2 * k + 1) * BigInt::from(n).pow(2 * k + 1),
        )).sum()
    }

    #[test]
    fn binary_split_matches_exact_series() {
        for n in [2, 5, 239] {
            for terms in 1..24 {
                let (_, q, t) = atan_split(n, 0, terms);
                assert_eq!(Q::new(t, q), direct_sum(n, terms));
            }
        }
    }

    #[test]
    fn outward_atan_contains_independent_alternating_bounds() {
        for n in [2, 5, 239] {
            for bits in [16, 32, 64, 128] {
                let result = atan_reciprocal(n, bits);
                let a = direct_sum(n, 100);
                let b = direct_sum(n, 101);
                assert!(result.lo <= a && a <= result.hi);
                assert!(result.lo <= b && b <= result.hi);
                assert!(result.hi - result.lo <= Q::new(2.into(), BigInt::one() << bits));
            }
        }
    }

    #[test]
    fn log_enclosure_contains_the_value_and_decides_the_lambert_domain() {
        for (q, value) in [(Q::from_integer(4.into()), 4f64.ln()), (Q::new(1.into(), 27.into()), (1.0f64 / 27.0).ln()),
                           (Q::new(5.into(), 7.into()), (5.0f64 / 7.0).ln()), (Q::from_integer(1000.into()), 1000f64.ln())] {
            let enclosure = log_point(&q, 80).unwrap();
            assert!(enclosure.lo.to_f64().unwrap() <= value + 1e-15 && enclosure.hi.to_f64().unwrap() >= value - 1e-15);
            assert!((&enclosure.hi - &enclosure.lo) < Q::new(BigInt::one(), BigInt::one() << 70));
        }
        let minus_log4 = crate::ops::neg(&crate::functions::log(&Expr::int(4)));
        assert_eq!(real(&minus_log4, 0), Some(false));
        assert_eq!(real(&minus_log4, -1), Some(false));
        let quarter_log2 = crate::ops::div(&crate::functions::log(&Expr::int(2)), &Expr::int(-4));
        assert_eq!(real(&quarter_log2, -1), Some(true));
    }

    #[test]
    fn pi_interval_agrees_with_independent_chudnovsky_value() {
        for bits in [32, 128, 512, 2048] {
            let interval = pi(bits);
            let value = crate::bigfloat::pi(bits + 64);
            let dyadic = value.to_ratio_exact();
            assert!(interval.lo <= dyadic && dyadic <= interval.hi);
            assert!(interval.hi - interval.lo <= Q::new(2.into(), BigInt::one() << bits));
        }
    }
}
