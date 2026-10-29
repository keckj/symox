//! A compact integer: machine word with transparent promotion to a big
//! integer. Used by the exact (fraction-free) matrix algorithms and the
//! internal polynomial arithmetic, where most coefficients stay small.

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::number::Number;

#[derive(Clone, Debug)]
pub enum Int {
    S(i64),
    /// Never fits in an `i64`.
    B(BigInt),
}

use Int::{B, S};

#[inline]
fn norm(b: BigInt) -> Int {
    match b.to_i64() {
        Some(v) => S(v),
        None => B(b),
    }
}

impl Int {
    pub const ZERO: Int = S(0);
    pub const ONE: Int = S(1);

    #[inline]
    pub fn from_i64(v: i64) -> Int {
        S(v)
    }

    pub fn from_big(b: BigInt) -> Int {
        norm(b)
    }

    pub fn to_big(&self) -> BigInt {
        match self {
            S(v) => BigInt::from(*v),
            B(b) => b.clone(),
        }
    }

    pub fn to_number(&self) -> Number {
        match self {
            S(v) => Number::Int(*v),
            B(b) => Number::from_bigint(b.clone()),
        }
    }

    /// The integer value of an exact integer `Number`.
    pub fn from_number(n: &Number) -> Option<Int> {
        match n {
            Number::Int(v) => Some(S(*v)),
            Number::BigInt(b) => Some(norm(b.clone())),
            _ => None,
        }
    }

    pub fn to_f64(&self) -> f64 {
        match self {
            S(v) => *v as f64,
            B(b) => b.to_f64().unwrap_or(f64::NAN),
        }
    }

    pub fn to_i64(&self) -> Option<i64> {
        match self {
            S(v) => Some(*v),
            B(_) => None,
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        matches!(self, S(0))
    }

    #[inline]
    pub fn is_one(&self) -> bool {
        matches!(self, S(1))
    }

    #[inline]
    pub fn is_minus_one(&self) -> bool {
        matches!(self, S(-1))
    }

    #[inline]
    pub fn is_negative(&self) -> bool {
        match self {
            S(v) => *v < 0,
            B(b) => b.is_negative(),
        }
    }

    pub fn signum(&self) -> i32 {
        match self {
            S(v) => v.signum() as i32,
            B(b) => {
                if b.is_negative() {
                    -1
                } else {
                    1
                }
            }
        }
    }

    pub fn neg(&self) -> Int {
        match self {
            S(v) => match v.checked_neg() {
                Some(r) => S(r),
                None => B(-BigInt::from(*v)),
            },
            B(b) => norm(-b),
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
        match (self, o) {
            (S(a), S(b)) => match a.checked_add(*b) {
                Some(r) => S(r),
                None => B(BigInt::from(*a) + BigInt::from(*b)),
            },
            (S(a), B(b)) | (B(b), S(a)) => norm(b + BigInt::from(*a)),
            (B(a), B(b)) => norm(a + b),
        }
    }

    #[inline]
    pub fn sub(&self, o: &Int) -> Int {
        match (self, o) {
            (S(a), S(b)) => match a.checked_sub(*b) {
                Some(r) => S(r),
                None => B(BigInt::from(*a) - BigInt::from(*b)),
            },
            (S(a), B(b)) => norm(BigInt::from(*a) - b),
            (B(a), S(b)) => norm(a - BigInt::from(*b)),
            (B(a), B(b)) => norm(a - b),
        }
    }

    #[inline]
    pub fn mul(&self, o: &Int) -> Int {
        match (self, o) {
            (S(a), S(b)) => match a.checked_mul(*b) {
                Some(r) => S(r),
                None => B(BigInt::from(*a) * BigInt::from(*b)),
            },
            (S(0), _) | (_, S(0)) => S(0),
            (S(a), B(b)) | (B(b), S(a)) => norm(b * BigInt::from(*a)),
            (B(a), B(b)) => norm(a * b),
        }
    }

    /// `self * a - o * b`, the elimination kernel.
    #[inline]
    pub fn mul_sub(&self, a: &Int, o: &Int, b: &Int) -> Int {
        if let (S(x), S(y), S(z), S(w)) = (self, a, o, b) {
            let r = (*x as i128) * (*y as i128) - (*z as i128) * (*w as i128);
            if let Ok(v) = i64::try_from(r) {
                return S(v);
            }
            return B(BigInt::from(r));
        }
        self.mul(a).sub(&o.mul(b))
    }

    /// Exact division (the result is unspecified when `o` does not divide
    /// `self`).
    #[inline]
    pub fn divexact(&self, o: &Int) -> Int {
        match (self, o) {
            (S(a), S(b)) => match a.checked_div(*b) {
                Some(r) => S(r),
                None => norm(BigInt::from(*a) / BigInt::from(*b)),
            },
            (B(a), S(b)) => norm(a / BigInt::from(*b)),
            (S(a), B(b)) => norm(BigInt::from(*a) / b),
            (B(a), B(b)) => norm(a / b),
        }
    }

    /// Exact division of `self * a - o * b` by `d` (Bareiss step).
    #[inline]
    pub fn mul_sub_div(&self, a: &Int, o: &Int, b: &Int, d: &Int) -> Int {
        if let (S(x), S(y), S(z), S(w), S(dd)) = (self, a, o, b, d) {
            let r = (*x as i128) * (*y as i128) - (*z as i128) * (*w as i128);
            let q = r / (*dd as i128);
            if let Ok(v) = i64::try_from(q) {
                return S(v);
            }
            return B(BigInt::from(q));
        }
        let r = self.mul_sub(a, o, b);
        if d.is_one() {
            r
        } else {
            r.divexact(d)
        }
    }

    /// Floor division (for non-negative operands this is exact division
    /// rounded down).
    pub fn divexact_floor(&self, o: &Int) -> Int {
        let (q, r) = self.div_rem(o);
        if !r.is_zero() && (r.is_negative() != o.is_negative()) {
            q.sub(&Int::ONE)
        } else {
            q
        }
    }

    /// Euclidean division with remainder of the same sign as the divisor
    /// is not needed: truncated division.
    pub fn div_rem(&self, o: &Int) -> (Int, Int) {
        match (self, o) {
            (S(a), S(b)) if !(*a == i64::MIN && *b == -1) => (S(a / b), S(a % b)),
            _ => {
                let (q, r) = self.to_big().div_rem(&o.to_big());
                (norm(q), norm(r))
            }
        }
    }

    pub fn is_divisible_by(&self, o: &Int) -> bool {
        if o.is_zero() {
            return self.is_zero();
        }
        self.div_rem(o).1.is_zero()
    }

    /// Non-negative gcd.
    pub fn gcd(&self, o: &Int) -> Int {
        match (self, o) {
            (S(a), S(b)) => {
                let g = a.unsigned_abs().gcd(&b.unsigned_abs());
                if g <= i64::MAX as u64 {
                    S(g as i64)
                } else {
                    B(BigInt::from(g))
                }
            }
            (S(0), x) | (x, S(0)) => x.abs(),
            (S(a), B(b)) | (B(b), S(a)) => {
                let r = (b % BigInt::from(*a)).to_i64().unwrap();
                S(a.unsigned_abs().gcd(&r.unsigned_abs()) as i64).abs_fix()
            }
            (B(a), B(b)) => norm(lehmer_gcd(a.magnitude().clone(), b.magnitude().clone()).into()),
        }
    }

    fn abs_fix(self) -> Int {
        match self {
            S(v) if v < 0 => S(v).neg(),
            x => x,
        }
    }

    pub fn lcm(&self, o: &Int) -> Int {
        if self.is_zero() || o.is_zero() {
            return S(0);
        }
        let g = self.gcd(o);
        self.divexact(&g).mul(o).abs()
    }

    pub fn pow(&self, n: u32) -> Int {
        let mut r = S(1);
        let mut b = self.clone();
        let mut n = n;
        while n > 0 {
            if n & 1 == 1 {
                r = r.mul(&b);
            }
            n >>= 1;
            if n > 0 {
                b = b.mul(&b);
            }
        }
        r
    }

    /// Integer square root when `self` is a perfect square.
    pub fn exact_sqrt(&self) -> Option<Int> {
        if self.is_negative() {
            return None;
        }
        let b = self.to_big();
        let r = b.sqrt();
        if &r * &r == b {
            Some(norm(r))
        } else {
            None
        }
    }

    pub fn bits(&self) -> u64 {
        match self {
            S(v) => 64 - v.unsigned_abs().leading_zeros() as u64,
            B(b) => b.bits(),
        }
    }
}

// ----------------------------------------------------------------------
// Lehmer gcd
// ----------------------------------------------------------------------

/// Lehmer's gcd (much faster than the binary gcd for multi-word numbers);
/// the shared implementation is `crate::bigint_gcd`.
pub fn lehmer_gcd(a: num_bigint::BigUint, b: num_bigint::BigUint) -> num_bigint::BigUint {
    crate::bigint_gcd::gcd_biguint(&a, &b)
}

impl From<i64> for Int {
    fn from(v: i64) -> Int {
        S(v)
    }
}

impl PartialEq for Int {
    #[inline]
    fn eq(&self, o: &Int) -> bool {
        match (self, o) {
            (S(a), S(b)) => a == b,
            (B(a), B(b)) => a == b,
            _ => false,
        }
    }
}
impl Eq for Int {}

impl Hash for Int {
    fn hash<H: Hasher>(&self, h: &mut H) {
        match self {
            S(v) => v.hash(h),
            B(b) => b.hash(h),
        }
    }
}

impl Ord for Int {
    fn cmp(&self, o: &Int) -> Ordering {
        match (self, o) {
            (S(a), S(b)) => a.cmp(b),
            (S(_), B(b)) => {
                if b.is_negative() {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            }
            (B(a), S(_)) => {
                if a.is_negative() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
            (B(a), B(b)) => a.cmp(b),
        }
    }
}
impl PartialOrd for Int {
    fn partial_cmp(&self, o: &Int) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

impl fmt::Display for Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            S(v) => write!(f, "{v}"),
            B(b) => write!(f, "{b}"),
        }
    }
}

impl Zero for Int {
    fn zero() -> Int {
        S(0)
    }
    fn is_zero(&self) -> bool {
        Int::is_zero(self)
    }
}

impl One for Int {
    fn one() -> Int {
        S(1)
    }
}

impl std::ops::Add for Int {
    type Output = Int;
    fn add(self, o: Int) -> Int {
        Int::add(&self, &o)
    }
}

impl std::ops::Mul for Int {
    type Output = Int;
    fn mul(self, o: Int) -> Int {
        Int::mul(&self, &o)
    }
}

/// Exact rational built on [`Int`]: `num / den`, `den > 0`, reduced.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Rat {
    pub num: Int,
    pub den: Int,
}

impl Rat {
    pub fn zero() -> Rat {
        Rat { num: S(0), den: S(1) }
    }
    pub fn one() -> Rat {
        Rat { num: S(1), den: S(1) }
    }
    pub fn from_int(i: Int) -> Rat {
        Rat { num: i, den: S(1) }
    }
    pub fn new(num: Int, den: Int) -> Rat {
        assert!(!den.is_zero(), "zero denominator");
        let g = num.gcd(&den);
        let (mut n, mut d) = if g.is_one() || g.is_zero() {
            (num, den)
        } else {
            (num.divexact(&g), den.divexact(&g))
        };
        if d.is_negative() {
            n = n.neg();
            d = d.neg();
        }
        if n.is_zero() {
            d = S(1);
        }
        Rat { num: n, den: d }
    }
    pub fn from_number(n: &Number) -> Option<Rat> {
        match n {
            Number::Int(v) => Some(Rat::from_int(S(*v))),
            Number::BigInt(b) => Some(Rat::from_int(norm(b.clone()))),
            Number::Rational(r) => Some(Rat {
                num: norm(r.numer().clone()),
                den: norm(r.denom().clone()),
            }),
            Number::Float(_) | Number::MpFloat(..) => None,
        }
    }
    pub fn to_number(&self) -> Number {
        if self.den.is_one() {
            self.num.to_number()
        } else {
            // already reduced with a positive denominator
            Number::from_ratio(num_rational::BigRational::new_raw(self.num.to_big(), self.den.to_big()))
        }
    }
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
    pub fn is_one(&self) -> bool {
        self.num.is_one() && self.den.is_one()
    }
    pub fn is_integer(&self) -> bool {
        self.den.is_one()
    }
    pub fn is_negative(&self) -> bool {
        self.num.is_negative()
    }
    pub fn neg(&self) -> Rat {
        Rat { num: self.num.neg(), den: self.den.clone() }
    }
    pub fn add(&self, o: &Rat) -> Rat {
        if self.den.is_one() && o.den.is_one() {
            return Rat::from_int(self.num.add(&o.num));
        }
        if self.den == o.den {
            return Rat::new(self.num.add(&o.num), self.den.clone());
        }
        Rat::new(
            self.num.mul(&o.den).add(&o.num.mul(&self.den)),
            self.den.mul(&o.den),
        )
    }
    pub fn sub(&self, o: &Rat) -> Rat {
        self.add(&o.neg())
    }
    pub fn mul(&self, o: &Rat) -> Rat {
        if self.den.is_one() && o.den.is_one() {
            return Rat::from_int(self.num.mul(&o.num));
        }
        Rat::new(self.num.mul(&o.num), self.den.mul(&o.den))
    }
    pub fn inv(&self) -> Rat {
        Rat::new(self.den.clone(), self.num.clone())
    }
    pub fn div(&self, o: &Rat) -> Rat {
        Rat::new(self.num.mul(&o.den), self.den.mul(&o.num))
    }
    pub fn cmp_value(&self, o: &Rat) -> Ordering {
        self.num.mul(&o.den).cmp(&o.num.mul(&self.den))
    }
    pub fn to_f64(&self) -> f64 {
        self.to_number().to_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overflow_promotes() {
        let a = Int::from(i64::MAX);
        let b = a.add(&Int::ONE);
        assert!(matches!(b, B(_)));
        assert_eq!(b.sub(&Int::ONE), a);
        let c = a.mul(&a).divexact(&a);
        assert_eq!(c, a);
        assert_eq!(Int::from(12).gcd(&Int::from(-18)), Int::from(6));
        assert_eq!(Int::from(i64::MIN).gcd(&Int::from(0)).to_big(), BigInt::from(i64::MIN).abs());
    }

    #[test]
    fn rationals() {
        let r = Rat::new(Int::from(6), Int::from(-4));
        assert_eq!(r.num, Int::from(-3));
        assert_eq!(r.den, Int::from(2));
        assert!(r.add(&Rat::new(Int::from(3), Int::from(2))).is_zero());
    }
}

#[cfg(test)]
mod gcd_tests {
    use super::*;

    #[test]
    fn lehmer_matches_binary() {
        let mut seed = 3u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seed
        };
        for _ in 0..200 {
            let mk = |k: usize, next: &mut dyn FnMut() -> u64| {
                let mut v = num_bigint::BigUint::from(1u32);
                for _ in 0..k {
                    v = v * num_bigint::BigUint::from(next()) + num_bigint::BigUint::from(next() % 1000);
                }
                v
            };
            let g = mk(1 + (next() % 4) as usize, &mut next);
            let a = mk(1 + (next() % 8) as usize, &mut next) * &g;
            let b = mk(1 + (next() % 8) as usize, &mut next) * &g;
            assert_eq!(lehmer_gcd(a.clone(), b.clone()), a.gcd(&b));
        }
    }
}
