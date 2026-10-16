//! Integers with a machine word fast path (`Int`) and exact rationals built
//! on top of them (`Rat`), the coefficient rings `ZZ` and `QQ` of the
//! polynomial module.

use std::cmp::Ordering;
use std::fmt;

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive};

use crate::number::Number;

/// An integer: `S` for values fitting in an `i64`, `B` otherwise
/// (invariant: a `B` never fits in an `i64`).
#[derive(Clone, PartialEq, Eq, Hash)]
pub enum Int {
    S(i64),
    B(Box<BigInt>),
}

impl fmt::Debug for Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Int::S(i) => write!(f, "{i}"),
            Int::B(b) => write!(f, "{b}"),
        }
    }
}

impl fmt::Display for Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl From<i64> for Int {
    #[inline]
    fn from(i: i64) -> Int {
        Int::S(i)
    }
}

impl From<BigInt> for Int {
    #[inline]
    fn from(b: BigInt) -> Int {
        Int::from_big(b)
    }
}

impl Int {
    #[inline]
    pub fn zero() -> Int {
        Int::S(0)
    }
    #[inline]
    pub fn one() -> Int {
        Int::S(1)
    }

    pub fn from_big(b: BigInt) -> Int {
        match b.to_i64() {
            Some(i) => Int::S(i),
            None => Int::B(Box::new(b)),
        }
    }

    pub fn from_u64(u: u64) -> Int {
        if u <= i64::MAX as u64 {
            Int::S(u as i64)
        } else {
            Int::B(Box::new(BigInt::from(u)))
        }
    }

    pub fn from_i128(v: i128) -> Int {
        if let Ok(i) = i64::try_from(v) {
            Int::S(i)
        } else {
            Int::B(Box::new(BigInt::from(v)))
        }
    }

    pub fn to_big(&self) -> BigInt {
        match self {
            Int::S(i) => BigInt::from(*i),
            Int::B(b) => (**b).clone(),
        }
    }

    #[inline]
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Int::S(i) => Some(*i),
            Int::B(_) => None,
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        matches!(self, Int::S(0))
    }
    #[inline]
    pub fn is_one(&self) -> bool {
        matches!(self, Int::S(1))
    }
    #[inline]
    pub fn is_minus_one(&self) -> bool {
        matches!(self, Int::S(-1))
    }

    #[inline]
    pub fn signum(&self) -> i32 {
        match self {
            Int::S(i) => i.signum() as i32,
            Int::B(b) => match b.sign() {
                Sign::Minus => -1,
                Sign::NoSign => 0,
                Sign::Plus => 1,
            },
        }
    }

    #[inline]
    pub fn is_negative(&self) -> bool {
        self.signum() < 0
    }
    #[inline]
    pub fn is_positive(&self) -> bool {
        self.signum() > 0
    }

    pub fn neg(&self) -> Int {
        match self {
            Int::S(i) => match i.checked_neg() {
                Some(v) => Int::S(v),
                None => Int::B(Box::new(-BigInt::from(*i))),
            },
            Int::B(b) => Int::from_big(-(**b).clone()),
        }
    }

    pub fn abs(&self) -> Int {
        if self.is_negative() {
            self.neg()
        } else {
            self.clone()
        }
    }

    #[inline]
    pub fn add(&self, o: &Int) -> Int {
        if let (Int::S(a), Int::S(b)) = (self, o) {
            if let Some(v) = a.checked_add(*b) {
                return Int::S(v);
            }
        }
        Int::from_big(self.to_big() + o.to_big())
    }

    #[inline]
    pub fn sub(&self, o: &Int) -> Int {
        if let (Int::S(a), Int::S(b)) = (self, o) {
            if let Some(v) = a.checked_sub(*b) {
                return Int::S(v);
            }
        }
        Int::from_big(self.to_big() - o.to_big())
    }

    #[inline]
    pub fn mul(&self, o: &Int) -> Int {
        match (self, o) {
            (Int::S(a), Int::S(b)) => match a.checked_mul(*b) {
                Some(v) => Int::S(v),
                None => Int::from_i128(*a as i128 * *b as i128),
            },
            (Int::S(0), _) | (_, Int::S(0)) => Int::S(0),
            (Int::S(a), Int::B(b)) | (Int::B(b), Int::S(a)) => Int::from_big(&**b * *a),
            (Int::B(a), Int::B(b)) => Int::from_big(&**a * &**b),
        }
    }

    #[inline]
    pub fn mul_i64(&self, k: i64) -> Int {
        self.mul(&Int::S(k))
    }

    /// `self += o`
    #[inline]
    pub fn add_assign(&mut self, o: &Int) {
        if let (Int::S(a), Int::S(b)) = (&*self, o) {
            if let Some(v) = a.checked_add(*b) {
                *self = Int::S(v);
                return;
            }
        }
        *self = self.add(o);
    }

    /// `self += a * b`
    #[inline]
    pub fn add_mul(&mut self, a: &Int, b: &Int) {
        if let (Int::S(s), Int::S(x), Int::S(y)) = (&*self, a, b) {
            let v = *s as i128 + *x as i128 * *y as i128;
            *self = Int::from_i128(v);
            return;
        }
        *self = self.add(&a.mul(b));
    }

    /// `self -= a * b`
    #[inline]
    pub fn sub_mul(&mut self, a: &Int, b: &Int) {
        if let (Int::S(s), Int::S(x), Int::S(y)) = (&*self, a, b) {
            let v = *s as i128 - *x as i128 * *y as i128;
            *self = Int::from_i128(v);
            return;
        }
        *self = self.sub(&a.mul(b));
    }

    /// Truncated division and remainder (like Rust's `/` and `%`).
    pub fn div_rem_trunc(&self, o: &Int) -> (Int, Int) {
        if let (Int::S(a), Int::S(b)) = (self, o) {
            if let (Some(q), Some(r)) = (a.checked_div(*b), a.checked_rem(*b)) {
                return (Int::S(q), Int::S(r));
            }
        }
        let (q, r) = self.to_big().div_rem(&o.to_big());
        (Int::from_big(q), Int::from_big(r))
    }

    /// Floor division and the non-negative (for `o > 0`) remainder.
    pub fn div_mod_floor(&self, o: &Int) -> (Int, Int) {
        if let (Int::S(a), Int::S(b)) = (self, o) {
            if *b != 0 && !(*a == i64::MIN && *b == -1) {
                return (Int::S(a.div_floor(b)), Int::S(a.mod_floor(b)));
            }
        }
        let (q, r) = self.to_big().div_mod_floor(&o.to_big());
        (Int::from_big(q), Int::from_big(r))
    }

    /// Remainder in `[0, |m|)`.
    pub fn modulo(&self, m: &Int) -> Int {
        match (self, m) {
            (Int::S(a), Int::S(b)) if *b != 0 => Int::S(a.rem_euclid(*b)),
            (Int::S(a), Int::B(b)) if b.is_positive() => {
                if *a >= 0 {
                    self.clone()
                } else {
                    Int::from_big(&**b + *a)
                }
            }
            (Int::B(a), Int::S(b)) if *b > 0 => {
                let r = (**a).mod_floor(&BigInt::from(*b));
                Int::S(r.to_i64().unwrap())
            }
            (Int::B(a), Int::B(b)) if b.is_positive() => Int::from_big((**a).mod_floor(b)),
            _ => {
                let mb = m.to_big().abs();
                Int::from_big(self.to_big().mod_floor(&mb))
            }
        }
    }

    /// Symmetric remainder in `(-m/2, m/2]` for `m > 0`.
    pub fn symmetric_mod(&self, m: &Int) -> Int {
        let r = self.modulo(m);
        // r > m/2  <=>  2r > m
        if r.mul_i64(2).cmp(m) == Ordering::Greater {
            r.sub(m)
        } else {
            r
        }
    }

    /// Exact quotient, `None` if `o` does not divide `self` (or `o == 0`).
    pub fn div_exact(&self, o: &Int) -> Option<Int> {
        if o.is_zero() {
            return None;
        }
        if o.is_one() {
            return Some(self.clone());
        }
        let (q, r) = self.div_rem_trunc(o);
        if r.is_zero() {
            Some(q)
        } else {
            None
        }
    }

    /// Exact quotient; the caller guarantees divisibility.
    pub fn div_exact_unchecked(&self, o: &Int) -> Int {
        if o.is_one() {
            return self.clone();
        }
        self.div_rem_trunc(o).0
    }

    pub fn gcd(&self, o: &Int) -> Int {
        if let (Int::S(a), Int::S(b)) = (self, o) {
            let g = a.unsigned_abs().gcd(&b.unsigned_abs());
            return Int::from_u64(g);
        }
        if self.is_zero() {
            return o.abs();
        }
        if o.is_zero() {
            return self.abs();
        }
        if self.is_one() || o.is_one() {
            return Int::one();
        }
        Int::from_big(self.to_big().gcd(&o.to_big()))
    }

    pub fn lcm(&self, o: &Int) -> Int {
        if self.is_zero() || o.is_zero() {
            return Int::zero();
        }
        let g = self.gcd(o);
        self.div_exact_unchecked(&g).mul(o).abs()
    }

    /// Extended gcd: `(g, s, t)` with `s*a + t*b = g`, `g >= 0`.
    pub fn gcdex(&self, o: &Int) -> (Int, Int, Int) {
        let (mut r0, mut r1) = (self.clone(), o.clone());
        let (mut s0, mut s1) = (Int::one(), Int::zero());
        let (mut t0, mut t1) = (Int::zero(), Int::one());
        while !r1.is_zero() {
            let (q, r) = r0.div_rem_trunc(&r1);
            r0 = std::mem::replace(&mut r1, r);
            let s2 = s0.sub(&q.mul(&s1));
            s0 = std::mem::replace(&mut s1, s2);
            let t2 = t0.sub(&q.mul(&t1));
            t0 = std::mem::replace(&mut t1, t2);
        }
        if r0.is_negative() {
            (r0.neg(), s0.neg(), t0.neg())
        } else {
            (r0, s0, t0)
        }
    }

    pub fn pow(&self, n: u32) -> Int {
        if let Int::S(a) = self {
            if let Some(v) = a.checked_pow(n) {
                return Int::S(v);
            }
        }
        Int::from_big(num_traits::pow::Pow::pow(&self.to_big(), n))
    }

    /// Number of bits of `|self|`.
    pub fn bits(&self) -> u64 {
        match self {
            Int::S(i) => 64 - i.unsigned_abs().leading_zeros() as u64,
            Int::B(b) => b.bits(),
        }
    }

    pub fn to_f64(&self) -> f64 {
        match self {
            Int::S(i) => *i as f64,
            Int::B(b) => b.to_f64().unwrap_or(if b.is_negative() { f64::NEG_INFINITY } else { f64::INFINITY }),
        }
    }

    /// `log2(|self|)` approximately (for bounds); `-inf` for zero.
    pub fn log2(&self) -> f64 {
        match self {
            Int::S(i) => (i.unsigned_abs() as f64).log2(),
            Int::B(b) => {
                let bits = b.bits();
                let shift = bits.saturating_sub(60);
                let top = (b.abs() >> shift).to_f64().unwrap();
                top.log2() + shift as f64
            }
        }
    }

    /// Residue modulo a small modulus `m > 0`, in `[0, m)`.
    #[inline]
    pub fn rem_u64(&self, m: u64) -> u64 {
        match self {
            Int::S(i) => {
                let r = (*i as i128).rem_euclid(m as i128);
                r as u64
            }
            Int::B(b) => {
                let r = (**b).mod_floor(&BigInt::from(m));
                r.to_u64().unwrap()
            }
        }
    }

    /// Integer square root of a non-negative integer (floor).
    pub fn isqrt(&self) -> Int {
        match self {
            Int::S(i) if *i >= 0 => {
                let mut r = (*i as f64).sqrt() as i64;
                while r > 0 && r.checked_mul(r).is_none_or(|v| v > *i) {
                    r -= 1;
                }
                while (r + 1).checked_mul(r + 1).is_some_and(|v| v <= *i) {
                    r += 1;
                }
                Int::S(r)
            }
            _ => Int::from_big(self.to_big().sqrt()),
        }
    }

    /// Exact `k`-th root if it exists.
    pub fn exact_root(&self, k: u32) -> Option<Int> {
        if k == 1 {
            return Some(self.clone());
        }
        if self.is_negative() {
            if k % 2 == 0 {
                return None;
            }
            return self.neg().exact_root(k).map(|r| r.neg());
        }
        let b = self.to_big();
        let r = b.nth_root(k);
        if num_traits::pow::Pow::pow(&r, k) == b {
            Some(Int::from_big(r))
        } else {
            None
        }
    }

    pub fn to_number(&self) -> Number {
        match self {
            Int::S(i) => Number::Int(*i),
            Int::B(b) => Number::BigInt((**b).clone()),
        }
    }

    pub fn from_number(n: &Number) -> Option<Int> {
        match n {
            Number::Int(i) => Some(Int::S(*i)),
            Number::BigInt(b) => Some(Int::from_big(b.clone())),
            _ => None,
        }
    }
}

impl Ord for Int {
    fn cmp(&self, o: &Int) -> Ordering {
        match (self, o) {
            (Int::S(a), Int::S(b)) => a.cmp(b),
            (Int::B(a), Int::B(b)) => a.cmp(b),
            (Int::S(_), Int::B(b)) => {
                if b.is_negative() {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            }
            (Int::B(a), Int::S(_)) => {
                if a.is_negative() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
        }
    }
}

impl PartialOrd for Int {
    fn partial_cmp(&self, o: &Int) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

// ----------------------------------------------------------------------
// rationals
// ----------------------------------------------------------------------

/// A reduced fraction `n/d` with `d > 0`.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Rat {
    pub n: Int,
    pub d: Int,
}

impl fmt::Debug for Rat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.d.is_one() {
            write!(f, "{:?}", self.n)
        } else {
            write!(f, "{:?}/{:?}", self.n, self.d)
        }
    }
}

impl From<Int> for Rat {
    fn from(n: Int) -> Rat {
        Rat { n, d: Int::one() }
    }
}

impl Rat {
    pub fn zero() -> Rat {
        Rat { n: Int::zero(), d: Int::one() }
    }
    pub fn one() -> Rat {
        Rat { n: Int::one(), d: Int::one() }
    }
    pub fn from_i64(i: i64) -> Rat {
        Rat { n: Int::S(i), d: Int::one() }
    }

    /// `n/d` reduced; panics if `d == 0`.
    pub fn new(n: Int, d: Int) -> Rat {
        assert!(!d.is_zero(), "zero denominator");
        if d.is_one() {
            return Rat { n, d };
        }
        let g = n.gcd(&d);
        let (mut n, mut d) = if g.is_one() {
            (n, d)
        } else {
            (n.div_exact_unchecked(&g), d.div_exact_unchecked(&g))
        };
        if d.is_negative() {
            n = n.neg();
            d = d.neg();
        }
        Rat { n, d }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
    #[inline]
    pub fn is_one(&self) -> bool {
        self.n.is_one() && self.d.is_one()
    }
    #[inline]
    pub fn is_integer(&self) -> bool {
        self.d.is_one()
    }
    #[inline]
    pub fn is_negative(&self) -> bool {
        self.n.is_negative()
    }

    pub fn neg(&self) -> Rat {
        Rat { n: self.n.neg(), d: self.d.clone() }
    }

    pub fn abs(&self) -> Rat {
        Rat { n: self.n.abs(), d: self.d.clone() }
    }

    pub fn add(&self, o: &Rat) -> Rat {
        if self.d.is_one() && o.d.is_one() {
            return Rat { n: self.n.add(&o.n), d: Int::one() };
        }
        if self.d == o.d {
            return Rat::new(self.n.add(&o.n), self.d.clone());
        }
        let g = self.d.gcd(&o.d);
        if g.is_one() {
            // already reduced: gcd(n1 d2 + n2 d1, d1 d2) = 1
            return Rat {
                n: self.n.mul(&o.d).add(&o.n.mul(&self.d)),
                d: self.d.mul(&o.d),
            };
        }
        let d1 = self.d.div_exact_unchecked(&g);
        let d2 = o.d.div_exact_unchecked(&g);
        let n = self.n.mul(&d2).add(&o.n.mul(&d1));
        Rat::new(n, d1.mul(&o.d))
    }

    pub fn sub(&self, o: &Rat) -> Rat {
        self.add(&o.neg())
    }

    pub fn mul(&self, o: &Rat) -> Rat {
        if self.d.is_one() && o.d.is_one() {
            return Rat { n: self.n.mul(&o.n), d: Int::one() };
        }
        if self.is_zero() || o.is_zero() {
            return Rat::zero();
        }
        let g1 = self.n.gcd(&o.d);
        let g2 = o.n.gcd(&self.d);
        let n = self.n.div_exact_unchecked(&g1).mul(&o.n.div_exact_unchecked(&g2));
        let d = self.d.div_exact_unchecked(&g2).mul(&o.d.div_exact_unchecked(&g1));
        Rat { n, d }
    }

    pub fn mul_int(&self, k: &Int) -> Rat {
        if self.d.is_one() {
            return Rat { n: self.n.mul(k), d: Int::one() };
        }
        Rat::new(self.n.mul(k), self.d.clone())
    }

    pub fn inv(&self) -> Rat {
        assert!(!self.is_zero(), "division by zero");
        if self.n.is_negative() {
            Rat { n: self.d.neg(), d: self.n.neg() }
        } else {
            Rat { n: self.d.clone(), d: self.n.clone() }
        }
    }

    pub fn div(&self, o: &Rat) -> Rat {
        self.mul(&o.inv())
    }

    pub fn pow(&self, k: u32) -> Rat {
        Rat { n: self.n.pow(k), d: self.d.pow(k) }
    }

    pub fn to_f64(&self) -> f64 {
        if self.d.is_one() {
            return self.n.to_f64();
        }
        let (a, b) = (self.n.to_f64(), self.d.to_f64());
        if a.is_finite() && b.is_finite() {
            a / b
        } else {
            self.to_number().to_f64()
        }
    }

    pub fn to_number(&self) -> Number {
        if self.d.is_one() {
            self.n.to_number()
        } else {
            Number::rational(self.n.to_big(), self.d.to_big()).unwrap()
        }
    }

    pub fn from_number(n: &Number) -> Option<Rat> {
        match n {
            Number::Int(i) => Some(Rat::from_i64(*i)),
            Number::BigInt(b) => Some(Rat::from(Int::from_big(b.clone()))),
            Number::Rational(r) => Some(Rat {
                n: Int::from_big(r.numer().clone()),
                d: Int::from_big(r.denom().clone()),
            }),
            Number::Float(_) => None,
        }
    }

    pub fn floor(&self) -> Int {
        self.n.div_mod_floor(&self.d).0
    }
}

impl Ord for Rat {
    fn cmp(&self, o: &Rat) -> Ordering {
        if self.d == o.d {
            return self.n.cmp(&o.n);
        }
        self.n.mul(&o.d).cmp(&o.n.mul(&self.d))
    }
}

impl PartialOrd for Rat {
    fn partial_cmp(&self, o: &Rat) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

// ----------------------------------------------------------------------
// the coefficient abstraction
// ----------------------------------------------------------------------

/// A commutative coefficient ring with exact division where possible.
pub trait Coeff: Clone + PartialEq + Eq + std::hash::Hash + fmt::Debug + Send + Sync + 'static {
    fn zero() -> Self;
    fn one() -> Self;
    fn from_i64(i: i64) -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    fn is_negative(&self) -> bool;
    /// `self / o` when exact (always for fields with `o != 0`).
    fn div_exact(&self, o: &Self) -> Option<Self>;
    fn add_assign(&mut self, o: &Self) {
        *self = self.add(o);
    }
    fn add_mul(&mut self, a: &Self, b: &Self) {
        *self = self.add(&a.mul(b));
    }
    fn to_rat(&self) -> Rat;
}

/// A field.
pub trait Field: Coeff {
    fn inv(&self) -> Self;
    fn div(&self, o: &Self) -> Self {
        self.mul(&o.inv())
    }
}

impl Coeff for Int {
    #[inline]
    fn zero() -> Self {
        Int::S(0)
    }
    #[inline]
    fn one() -> Self {
        Int::S(1)
    }
    #[inline]
    fn from_i64(i: i64) -> Self {
        Int::S(i)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        Int::is_zero(self)
    }
    #[inline]
    fn is_one(&self) -> bool {
        Int::is_one(self)
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        Int::add(self, o)
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        Int::sub(self, o)
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Int::mul(self, o)
    }
    #[inline]
    fn neg(&self) -> Self {
        Int::neg(self)
    }
    #[inline]
    fn is_negative(&self) -> bool {
        Int::is_negative(self)
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        Int::div_exact(self, o)
    }
    #[inline]
    fn add_assign(&mut self, o: &Self) {
        Int::add_assign(self, o)
    }
    #[inline]
    fn add_mul(&mut self, a: &Self, b: &Self) {
        Int::add_mul(self, a, b)
    }
    fn to_rat(&self) -> Rat {
        Rat::from(self.clone())
    }
}

impl Coeff for Rat {
    fn zero() -> Self {
        Rat::zero()
    }
    fn one() -> Self {
        Rat::one()
    }
    fn from_i64(i: i64) -> Self {
        Rat::from_i64(i)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        Rat::is_zero(self)
    }
    #[inline]
    fn is_one(&self) -> bool {
        Rat::is_one(self)
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        Rat::add(self, o)
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        Rat::sub(self, o)
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Rat::mul(self, o)
    }
    #[inline]
    fn neg(&self) -> Self {
        Rat::neg(self)
    }
    fn is_negative(&self) -> bool {
        Rat::is_negative(self)
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        if o.is_zero() {
            None
        } else {
            Some(Rat::div(self, o))
        }
    }
    fn to_rat(&self) -> Rat {
        self.clone()
    }
}

impl Field for Rat {
    fn inv(&self) -> Self {
        Rat::inv(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn int_basic() {
        let a = Int::S(i64::MAX);
        let b = a.add(&Int::one());
        assert!(matches!(b, Int::B(_)));
        assert_eq!(b.sub(&Int::one()), a);
        assert_eq!(Int::S(12).gcd(&Int::S(-18)), Int::S(6));
        assert_eq!(Int::S(-7).symmetric_mod(&Int::S(5)), Int::S(-2));
        assert_eq!(Int::S(7).symmetric_mod(&Int::S(5)), Int::S(2));
        assert_eq!(Int::S(-7).modulo(&Int::S(5)), Int::S(3));
        let big = Int::S(3).pow(50);
        assert_eq!(big.div_exact(&Int::S(3).pow(48)), Some(Int::S(9)));
        assert_eq!(Int::S(17).isqrt(), Int::S(4));
        let (g, s, t) = Int::S(240).gcdex(&Int::S(46));
        assert_eq!(g, Int::S(2));
        assert_eq!(s.mul(&Int::S(240)).add(&t.mul(&Int::S(46))), Int::S(2));
    }

    #[test]
    fn rat_basic() {
        let a = Rat::new(Int::S(1), Int::S(2));
        let b = Rat::new(Int::S(1), Int::S(3));
        assert_eq!(a.add(&b), Rat::new(Int::S(5), Int::S(6)));
        assert_eq!(a.mul(&b), Rat::new(Int::S(1), Int::S(6)));
        assert_eq!(a.sub(&a), Rat::zero());
        assert_eq!(Rat::new(Int::S(2), Int::S(-4)), Rat::new(Int::S(-1), Int::S(2)));
        assert!(a > b);
    }
}
