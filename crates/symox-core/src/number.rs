//! Exact and floating point numbers.
//!
//! `Number` is the numeric tower used by expressions: machine integers with a
//! transparent promotion to big integers, exact rationals and IEEE doubles.
//! All constructors normalise their result so that structural equality is
//! value equality for exact numbers (an integer that fits in an `i64` is
//! always stored as `Int`, a rational always has a denominator > 1).

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

#[derive(Clone, Debug)]
pub enum Number {
    /// Integer fitting in a machine word.
    Int(i64),
    /// Integer that does not fit in an `i64`.
    BigInt(BigInt),
    /// Reduced rational with denominator > 1.
    Rational(BigRational),
    /// Double precision float (never `-0.0`).
    Float(f64),
}

use Number::*;

impl Number {
    pub const ZERO: Number = Int(0);
    pub const ONE: Number = Int(1);
    pub const MINUS_ONE: Number = Int(-1);

    #[inline]
    pub fn int(i: i64) -> Number {
        Int(i)
    }

    pub fn from_bigint(b: BigInt) -> Number {
        match b.to_i64() {
            Some(i) => Int(i),
            None => BigInt(b),
        }
    }

    pub fn float(f: f64) -> Number {
        // normalise -0.0 so that hashing/equality are consistent
        Float(if f == 0.0 { 0.0 } else { f })
    }

    /// Build `n/d` in lowest terms. Returns `None` if `d == 0`.
    pub fn rational(n: BigInt, d: BigInt) -> Option<Number> {
        if d.is_zero() {
            return None;
        }
        let r = BigRational::new(n, d);
        Some(Self::from_ratio(r))
    }

    pub fn rational_i64(n: i64, d: i64) -> Option<Number> {
        if d == 0 {
            return None;
        }
        let g = n.gcd(&d);
        let (mut n, mut d) = (n / g, d / g);
        if d < 0 {
            // i64::MIN corner case handled by the big path
            if n == i64::MIN || d == i64::MIN {
                return Self::rational(BigInt::from(n), BigInt::from(d));
            }
            n = -n;
            d = -d;
        }
        if d == 1 {
            Some(Int(n))
        } else {
            Some(Rational(BigRational::new_raw(BigInt::from(n), BigInt::from(d))))
        }
    }

    pub fn from_ratio(r: BigRational) -> Number {
        if r.denom().is_one() {
            Self::from_bigint(r.numer().clone())
        } else {
            Rational(r)
        }
    }

    // ------------------------------------------------------------------
    // predicates
    // ------------------------------------------------------------------

    #[inline]
    pub fn is_zero(&self) -> bool {
        match self {
            Int(i) => *i == 0,
            Float(f) => *f == 0.0,
            _ => false,
        }
    }

    #[inline]
    pub fn is_one(&self) -> bool {
        matches!(self, Int(1))
    }

    #[inline]
    pub fn is_minus_one(&self) -> bool {
        matches!(self, Int(-1))
    }

    #[inline]
    pub fn is_integer(&self) -> bool {
        matches!(self, Int(_) | BigInt(_))
    }

    #[inline]
    pub fn is_rational(&self) -> bool {
        !matches!(self, Float(_))
    }

    #[inline]
    pub fn is_float(&self) -> bool {
        matches!(self, Float(_))
    }

    pub fn is_negative(&self) -> bool {
        match self {
            Int(i) => *i < 0,
            BigInt(b) => b.is_negative(),
            Rational(r) => r.is_negative(),
            Float(f) => *f < 0.0,
        }
    }

    pub fn is_positive(&self) -> bool {
        match self {
            Int(i) => *i > 0,
            BigInt(b) => b.is_positive(),
            Rational(r) => r.is_positive(),
            Float(f) => *f > 0.0,
        }
    }

    pub fn is_nan(&self) -> bool {
        matches!(self, Float(f) if f.is_nan())
    }

    pub fn to_i64(&self) -> Option<i64> {
        match self {
            Int(i) => Some(*i),
            _ => None,
        }
    }

    pub fn to_f64(&self) -> f64 {
        match self {
            Int(i) => *i as f64,
            BigInt(b) => b.to_f64().unwrap_or(f64::NAN),
            Rational(r) => ratio_to_f64(r),
            Float(f) => *f,
        }
    }

    pub fn to_bigint(&self) -> Option<BigInt> {
        match self {
            Int(i) => Some(BigInt::from(*i)),
            BigInt(b) => Some(b.clone()),
            _ => None,
        }
    }

    pub fn to_ratio(&self) -> Option<BigRational> {
        match self {
            Int(i) => Some(BigRational::from_integer(BigInt::from(*i))),
            BigInt(b) => Some(BigRational::from_integer(b.clone())),
            Rational(r) => Some(r.clone()),
            Float(_) => None,
        }
    }

    /// Numerator and denominator of an exact number.
    pub fn numer_denom(&self) -> Option<(BigInt, BigInt)> {
        match self {
            Int(i) => Some((BigInt::from(*i), BigInt::one())),
            BigInt(b) => Some((b.clone(), BigInt::one())),
            Rational(r) => Some((r.numer().clone(), r.denom().clone())),
            Float(_) => None,
        }
    }

    pub fn numer(&self) -> Number {
        match self {
            Rational(r) => Self::from_bigint(r.numer().clone()),
            _ => self.clone(),
        }
    }

    pub fn denom(&self) -> Number {
        match self {
            Rational(r) => Self::from_bigint(r.denom().clone()),
            _ => Int(1),
        }
    }

    // ------------------------------------------------------------------
    // arithmetic
    // ------------------------------------------------------------------

    pub fn neg(&self) -> Number {
        match self {
            Int(i) => match i.checked_neg() {
                Some(v) => Int(v),
                None => BigInt(-BigInt::from(*i)),
            },
            BigInt(b) => Self::from_bigint(-b),
            Rational(r) => Rational(-r),
            Float(f) => Self::float(-f),
        }
    }

    pub fn abs(&self) -> Number {
        if self.is_negative() {
            self.neg()
        } else {
            self.clone()
        }
    }

    pub fn add(&self, o: &Number) -> Number {
        match (self, o) {
            (Int(a), Int(b)) => match a.checked_add(*b) {
                Some(v) => Int(v),
                None => BigInt(BigInt::from(*a) + BigInt::from(*b)),
            },
            (Float(a), _) => Self::float(a + o.to_f64()),
            (_, Float(b)) => Self::float(self.to_f64() + b),
            (Rational(_), _) | (_, Rational(_)) => {
                Self::from_ratio(self.to_ratio().unwrap() + o.to_ratio().unwrap())
            }
            _ => Self::from_bigint(self.to_bigint().unwrap() + o.to_bigint().unwrap()),
        }
    }

    pub fn sub(&self, o: &Number) -> Number {
        match (self, o) {
            (Int(a), Int(b)) => match a.checked_sub(*b) {
                Some(v) => Int(v),
                None => BigInt(BigInt::from(*a) - BigInt::from(*b)),
            },
            _ => self.add(&o.neg()),
        }
    }

    pub fn mul(&self, o: &Number) -> Number {
        match (self, o) {
            (Int(a), Int(b)) => match a.checked_mul(*b) {
                Some(v) => Int(v),
                None => BigInt(BigInt::from(*a) * BigInt::from(*b)),
            },
            (Float(a), _) => Self::float(a * o.to_f64()),
            (_, Float(b)) => Self::float(self.to_f64() * b),
            (Rational(_), _) | (_, Rational(_)) => {
                Self::from_ratio(self.to_ratio().unwrap() * o.to_ratio().unwrap())
            }
            _ => Self::from_bigint(self.to_bigint().unwrap() * o.to_bigint().unwrap()),
        }
    }

    /// Division. Returns `None` when dividing an exact number by exact zero.
    pub fn div(&self, o: &Number) -> Option<Number> {
        match (self, o) {
            (Float(a), _) => Some(Self::float(a / o.to_f64())),
            (_, Float(b)) => Some(Self::float(self.to_f64() / b)),
            _ if o.is_zero() => None,
            (Int(a), Int(b)) => Self::rational_i64(*a, *b),
            _ => Some(Self::from_ratio(
                self.to_ratio().unwrap() / o.to_ratio().unwrap(),
            )),
        }
    }

    pub fn inv(&self) -> Option<Number> {
        Int(1).div(self)
    }

    /// `self ** n` for an integer exponent.
    /// Returns `None` for `0 ** negative` with an exact zero.
    pub fn pow_i64(&self, n: i64) -> Option<Number> {
        if let Float(f) = self {
            return Some(Self::float(f.powf(n as f64)));
        }
        if n < 0 {
            if self.is_zero() {
                return None;
            }
            let p = self.pow_u64(n.unsigned_abs());
            return p.inv();
        }
        Some(self.pow_u64(n as u64))
    }

    fn pow_u64(&self, n: u64) -> Number {
        match self {
            Int(a) => {
                if let Ok(n32) = u32::try_from(n) {
                    if let Some(v) = a.checked_pow(n32) {
                        return Int(v);
                    }
                }
                Self::from_bigint(big_pow(&BigInt::from(*a), n))
            }
            BigInt(b) => Self::from_bigint(big_pow(b, n)),
            Rational(r) => {
                let num = big_pow(r.numer(), n);
                let den = big_pow(r.denom(), n);
                Self::from_ratio(BigRational::new_raw(num, den))
            }
            Float(f) => Self::float(f.powf(n as f64)),
        }
    }

    /// floor of an exact number
    pub fn floor(&self) -> Number {
        match self {
            Rational(r) => Self::from_bigint(r.floor().to_integer()),
            Float(f) => Self::float(f.floor()),
            _ => self.clone(),
        }
    }

    /// Compare numerical values (`None` only for NaN).
    pub fn cmp_value(&self, o: &Number) -> Option<Ordering> {
        match (self, o) {
            (Int(a), Int(b)) => Some(a.cmp(b)),
            (Float(_), _) | (_, Float(_)) => self.to_f64().partial_cmp(&o.to_f64()),
            _ => Some(self.to_ratio().unwrap().cmp(&o.to_ratio().unwrap())),
        }
    }

    pub fn gcd(&self, o: &Number) -> Number {
        match (self, o) {
            (Int(a), Int(b)) => {
                let g = a.gcd(b);
                if g < 0 {
                    // only possible for i64::MIN
                    Self::from_bigint(BigInt::from(a.unsigned_abs().gcd(&b.unsigned_abs())))
                } else {
                    Int(g)
                }
            }
            _ => match (self.to_bigint(), o.to_bigint()) {
                (Some(a), Some(b)) => Self::from_bigint(a.gcd(&b)),
                _ => {
                    // gcd of rationals: gcd(numerators)/lcm(denominators)
                    match (self.numer_denom(), o.numer_denom()) {
                        (Some((an, ad)), Some((bn, bd))) => {
                            Self::rational(an.gcd(&bn), ad.lcm(&bd)).unwrap()
                        }
                        _ => Int(1),
                    }
                }
            },
        }
    }

    /// Structural ordering used for canonical sorting (not numerical).
    pub(crate) fn cmp_struct(&self, o: &Number) -> Ordering {
        fn rank(n: &Number) -> u8 {
            match n {
                Int(_) => 0,
                BigInt(_) => 1,
                Rational(_) => 2,
                Float(_) => 3,
            }
        }
        match (self, o) {
            (Int(a), Int(b)) => a.cmp(b),
            (BigInt(a), BigInt(b)) => a.cmp(b),
            (Rational(a), Rational(b)) => a.cmp(b),
            (Float(a), Float(b)) => a.to_bits().cmp(&b.to_bits()),
            _ => rank(self).cmp(&rank(o)),
        }
    }
}

fn big_pow(b: &BigInt, n: u64) -> BigInt {
    num_traits::pow::Pow::pow(b, n)
}

fn ratio_to_f64(r: &BigRational) -> f64 {
    match (r.numer().to_f64(), r.denom().to_f64()) {
        (Some(n), Some(d)) if n.is_finite() && d.is_finite() => n / d,
        _ => {
            // scale down huge values
            let nb = r.numer().bits() as i64;
            let db = r.denom().bits() as i64;
            let shift = (nb - db) - 60;
            let (n, d) = if shift > 0 {
                (r.numer().clone(), r.denom() << (shift as usize))
            } else {
                (r.numer() << ((-shift) as usize), r.denom().clone())
            };
            let q = (n / d).to_f64().unwrap_or(f64::NAN);
            q * 2f64.powi(shift as i32)
        }
    }
}

impl PartialEq for Number {
    fn eq(&self, o: &Number) -> bool {
        match (self, o) {
            (Int(a), Int(b)) => a == b,
            (BigInt(a), BigInt(b)) => a == b,
            (Rational(a), Rational(b)) => a == b,
            (Float(a), Float(b)) => a.to_bits() == b.to_bits(),
            _ => false,
        }
    }
}
impl Eq for Number {}

impl Hash for Number {
    fn hash<H: Hasher>(&self, h: &mut H) {
        match self {
            Int(i) => {
                0u8.hash(h);
                i.hash(h)
            }
            BigInt(b) => {
                1u8.hash(h);
                b.hash(h)
            }
            Rational(r) => {
                2u8.hash(h);
                r.numer().hash(h);
                r.denom().hash(h)
            }
            Float(f) => {
                3u8.hash(h);
                f.to_bits().hash(h)
            }
        }
    }
}

impl fmt::Display for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Int(i) => write!(f, "{i}"),
            BigInt(b) => write!(f, "{b}"),
            Rational(r) => write!(f, "{}/{}", r.numer(), r.denom()),
            Float(x) => write!(f, "{}", format_float(*x)),
        }
    }
}

/// Format a float with 15 significant digits like SymPy's default printer.
pub fn format_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return "0.0".into();
    }
    let a = x.abs();
    if (1e-4..1e16).contains(&a) {
        let digits = 15i32 - 1 - a.log10().floor() as i32;
        let digits = digits.max(0) as usize;
        let mut s = format!("{x:.digits$}");
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.push('0');
            }
        } else {
            s.push_str(".0");
        }
        s
    } else {
        let s = format!("{x:.14e}");
        // mantissa cleanup: 1.50000000000000e-5 -> 1.5e-5
        let (m, e) = s.split_once('e').unwrap();
        let mut m = m.to_string();
        if m.contains('.') {
            while m.ends_with('0') {
                m.pop();
            }
            if m.ends_with('.') {
                m.push('0');
            }
        }
        let e: i32 = e.parse().unwrap();
        format!("{m}e{}{}", if e >= 0 { "+" } else { "-" }, e.abs())
    }
}

/// Format a float with exactly 15 significant digits (SymPy's full
/// precision output for a standalone `Float`).
pub fn format_float_full(x: f64) -> String {
    if !x.is_finite() || x == 0.0 {
        return if x == 0.0 { "0.0".into() } else { format_float(x) };
    }
    let a = x.abs();
    let e10 = a.log10().floor() as i32;
    // rounding may bump the exponent (9.99999999999999999 -> 10.0)
    let rounded: f64 = format!("{a:.14e}").parse().unwrap_or(a);
    let e10 = if rounded >= 10f64.powi(e10 + 1) { e10 + 1 } else { e10 };
    if (-4..15).contains(&e10) {
        let prec = (14 - e10) as usize;
        format!("{x:.prec$}")
    } else {
        let s = format!("{x:.14e}");
        let (m, e) = s.split_once('e').unwrap();
        let e: i32 = e.parse().unwrap();
        format!("{m}e{}{}", if e >= 0 { "+" } else { "-" }, e.abs())
    }
}

// ----------------------------------------------------------------------
// integer helpers used for exact radicals
// ----------------------------------------------------------------------

/// Exact `k`-th root of a non-negative integer if it exists.
pub fn exact_root(n: &BigInt, k: u32) -> Option<BigInt> {
    if n.sign() == Sign::Minus {
        return None;
    }
    let r = n.nth_root(k);
    if num_traits::pow::Pow::pow(&r, k) == *n {
        Some(r)
    } else {
        None
    }
}

/// Partial factorisation by trial division up to `limit`.
/// Returns the list of `(prime, multiplicity)` and the unfactored cofactor.
pub fn trial_factor(n: &BigInt, limit: u64) -> (Vec<(BigInt, u32)>, BigInt) {
    let mut n = n.abs();
    let mut out = Vec::new();
    if n.is_zero() {
        return (out, n);
    }
    // fast path in u64
    if let Some(mut m) = n.to_u64() {
        let push = |p: u64, m: &mut u64, out: &mut Vec<(BigInt, u32)>| {
            let mut k = 0;
            while *m % p == 0 {
                *m /= p;
                k += 1;
            }
            if k > 0 {
                out.push((BigInt::from(p), k));
            }
        };
        push(2, &mut m, &mut out);
        push(3, &mut m, &mut out);
        let mut p = 5u64;
        while p <= limit && p.saturating_mul(p) <= m {
            push(p, &mut m, &mut out);
            push(p + 2, &mut m, &mut out);
            p += 6;
        }
        if m > 1 && (m as u128) < (limit as u128 + 1) * (limit as u128 + 1) {
            out.push((BigInt::from(m), 1));
            m = 1;
        }
        return (out, BigInt::from(m));
    }
    let mut p = 2u64;
    while p <= limit {
        let bp = BigInt::from(p);
        if &bp * &bp > n {
            break;
        }
        let mut k = 0;
        loop {
            let (q, r) = n.div_rem(&bp);
            if r.is_zero() {
                n = q;
                k += 1;
            } else {
                break;
            }
        }
        if k > 0 {
            out.push((bp, k));
        }
        p += if p == 2 { 1 } else { 2 };
    }
    if !n.is_one() {
        let lim = BigInt::from(limit);
        if n <= &lim * &lim {
            out.push((n, 1));
            n = BigInt::one();
        }
    }
    (out, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arith() {
        let a = Number::rational_i64(1, 2).unwrap();
        let b = Number::rational_i64(1, 3).unwrap();
        assert_eq!(a.add(&b), Number::rational_i64(5, 6).unwrap());
        assert_eq!(a.mul(&Int(2)), Int(1));
        assert_eq!(Int(i64::MAX).add(&Int(1)).sub(&Int(1)), Int(i64::MAX));
        assert_eq!(Int(2).pow_i64(-2).unwrap(), Number::rational_i64(1, 4).unwrap());
        assert!(Int(0).pow_i64(-1).is_none());
    }

    #[test]
    fn factor() {
        let (f, c) = trial_factor(&BigInt::from(360), 1000);
        assert_eq!(c, BigInt::one());
        assert_eq!(
            f,
            vec![(BigInt::from(2), 3), (BigInt::from(3), 2), (BigInt::from(5), 1)]
        );
    }

    #[test]
    fn float_fmt() {
        assert_eq!(format_float(1.5), "1.5");
        assert_eq!(format_float(2.0), "2.0");
        assert_eq!(format_float(1.0 / 3.0), "0.333333333333333");
        assert_eq!(format_float(std::f64::consts::PI), "3.14159265358979");
        assert_eq!(format_float_full(2.0), "2.00000000000000");
        assert_eq!(format_float_full(0.5), "0.500000000000000");
        assert_eq!(format_float_full(1e-5), "1.00000000000000e-5");
        assert_eq!(format_float_full(1e20), "1.00000000000000e+20");
        assert_eq!(format_float_full(123.456), "123.456000000000");
        assert_eq!(format_float_full(0.0001), "0.000100000000000000");
        assert_eq!(format_float_full(1e15), "1.00000000000000e+15");
    }
}
