//! Binary exponents of unbounded size, as mpmath's Python-integer exponents.
//!
//! Exponents up to `SMALL_LIMIT` in magnitude are stored inline as i128,
//! so ordinary values never allocate for their exponent and the common
//! arithmetic is two machine operations. Larger exponents are big integers.
//! The representation is canonical: a value is `Small` exactly when it fits
//! the inline range, so derived equality and hashing are by value.
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive};

/// Inline exponents satisfy |e| <= SMALL_LIMIT, so the sum or difference of
/// two of them, or of one and any i64 offset, never overflows i128.
pub const SMALL_LIMIT: i128 = 1 << 120;

#[derive(Clone, PartialEq, Eq, Hash)]
pub enum Exponent {
    Small(i128),
    Big(Box<BigInt>),
}

impl fmt::Debug for Exponent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Exponent::Small(v) => write!(f, "{v}"),
            Exponent::Big(v) => write!(f, "{v}"),
        }
    }
}

impl fmt::Display for Exponent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl Default for Exponent {
    fn default() -> Self {
        Exponent::Small(0)
    }
}

impl Exponent {
    pub const ZERO: Exponent = Exponent::Small(0);

    #[inline]
    pub fn from_i128(v: i128) -> Exponent {
        if v.unsigned_abs() <= SMALL_LIMIT as u128 { Exponent::Small(v) } else { Exponent::Big(Box::new(BigInt::from(v))) }
    }

    pub fn from_bigint(v: BigInt) -> Exponent {
        match v.to_i128() {
            Some(s) if s.unsigned_abs() <= SMALL_LIMIT as u128 => Exponent::Small(s),
            _ => Exponent::Big(Box::new(v)),
        }
    }

    /// The exponent as i128 when it is in the inline range.
    #[inline]
    pub fn small(&self) -> Option<i128> {
        match self {
            Exponent::Small(v) => Some(*v),
            Exponent::Big(_) => None,
        }
    }

    pub fn to_bigint(&self) -> BigInt {
        match self {
            Exponent::Small(v) => BigInt::from(*v),
            Exponent::Big(v) => (**v).clone(),
        }
    }

    /// i128 value clamped to +-2**124. Callers that only size guard bits or
    /// test a range keep their former behaviour for wide exponents.
    #[inline]
    pub fn saturating_i128(&self) -> i128 {
        match self {
            Exponent::Small(v) => *v,
            Exponent::Big(v) => if v.is_negative() { -(1 << 124) } else { 1 << 124 },
        }
    }

    /// |e| as u128, u128::MAX for exponents outside the inline range.
    #[inline]
    pub fn unsigned_abs(&self) -> u128 {
        match self {
            Exponent::Small(v) => v.unsigned_abs(),
            Exponent::Big(_) => u128::MAX,
        }
    }

    pub fn to_usize(&self) -> Option<usize> {
        self.small().and_then(|v| usize::try_from(v).ok())
    }

    pub fn to_i64(&self) -> Option<i64> {
        self.small().and_then(|v| i64::try_from(v).ok())
    }

    #[inline]
    pub fn is_negative(&self) -> bool {
        match self {
            Exponent::Small(v) => *v < 0,
            Exponent::Big(v) => v.is_negative(),
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        matches!(self, Exponent::Small(0))
    }

    pub fn is_odd(&self) -> bool {
        match self {
            Exponent::Small(v) => v & 1 != 0,
            Exponent::Big(v) => v.is_odd(),
        }
    }

    /// Number of bits of |e|.
    pub fn bits(&self) -> u64 {
        match self {
            Exponent::Small(v) => u64::from(128 - v.unsigned_abs().leading_zeros()),
            Exponent::Big(v) => v.bits(),
        }
    }

    pub fn abs(&self) -> Exponent {
        if self.is_negative() { -self } else { self.clone() }
    }

    /// Floor division by two.
    pub fn half_floor(&self) -> Exponent {
        match self {
            Exponent::Small(v) => Exponent::Small(v.div_euclid(2)),
            Exponent::Big(v) => Exponent::from_bigint(v.div_floor(&BigInt::from(2))),
        }
    }

    pub fn rem_euclid_u64(&self, modulus: u64) -> u64 {
        match self {
            Exponent::Small(v) => v.rem_euclid(modulus as i128) as u64,
            Exponent::Big(v) => v.mod_floor(&BigInt::from(modulus)).to_u64().unwrap(),
        }
    }

    pub fn min_ref<'a>(&'a self, other: &'a Exponent) -> &'a Exponent {
        if self <= other { self } else { other }
    }

    pub fn max_ref<'a>(&'a self, other: &'a Exponent) -> &'a Exponent {
        if self >= other { self } else { other }
    }

    #[inline]
    fn add_ref(&self, other: &Exponent) -> Exponent {
        if let (Exponent::Small(a), Exponent::Small(b)) = (self, other) {
            return Exponent::from_i128(a + b);
        }
        Exponent::from_bigint(self.to_bigint() + other.to_bigint())
    }

    #[inline]
    fn add_i128(&self, other: i128) -> Exponent {
        if let Exponent::Small(a) = self {
            if let Some(s) = a.checked_add(other) { return Exponent::from_i128(s); }
        }
        Exponent::from_bigint(self.to_bigint() + BigInt::from(other))
    }
}

impl From<i128> for Exponent {
    #[inline]
    fn from(v: i128) -> Self { Exponent::from_i128(v) }
}
impl From<i64> for Exponent {
    #[inline]
    fn from(v: i64) -> Self { Exponent::Small(v as i128) }
}
impl From<i32> for Exponent {
    #[inline]
    fn from(v: i32) -> Self { Exponent::Small(v as i128) }
}
impl From<u32> for Exponent {
    #[inline]
    fn from(v: u32) -> Self { Exponent::Small(v as i128) }
}
impl From<u64> for Exponent {
    #[inline]
    fn from(v: u64) -> Self { Exponent::Small(v as i128) }
}
impl From<BigInt> for Exponent {
    fn from(v: BigInt) -> Self { Exponent::from_bigint(v) }
}
impl From<&BigInt> for Exponent {
    fn from(v: &BigInt) -> Self { Exponent::from_bigint(v.clone()) }
}

impl Ord for Exponent {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Exponent::Small(a), Exponent::Small(b)) => a.cmp(b),
            (Exponent::Big(a), Exponent::Big(b)) => a.cmp(b),
            // Canonical form: a big exponent lies outside the inline range.
            (Exponent::Big(a), Exponent::Small(_)) => if a.is_negative() { Ordering::Less } else { Ordering::Greater },
            (Exponent::Small(_), Exponent::Big(b)) => if b.is_negative() { Ordering::Greater } else { Ordering::Less },
        }
    }
}
impl PartialOrd for Exponent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}
impl PartialEq<i128> for Exponent {
    #[inline]
    fn eq(&self, other: &i128) -> bool { matches!(self, Exponent::Small(v) if v == other) }
}
impl PartialOrd<i128> for Exponent {
    #[inline]
    fn partial_cmp(&self, other: &i128) -> Option<Ordering> {
        Some(match self {
            Exponent::Small(v) => v.cmp(other),
            Exponent::Big(v) => if v.is_negative() { Ordering::Less } else { Ordering::Greater },
        })
    }
}

macro_rules! binop {
    ($trait:ident, $method:ident, $body:expr) => {
        impl $trait<&Exponent> for &Exponent {
            type Output = Exponent;
            #[inline]
            fn $method(self, other: &Exponent) -> Exponent { let f: fn(&Exponent, &Exponent) -> Exponent = $body; f(self, other) }
        }
        impl $trait<Exponent> for Exponent {
            type Output = Exponent;
            #[inline]
            fn $method(self, other: Exponent) -> Exponent { (&self).$method(&other) }
        }
        impl $trait<&Exponent> for Exponent {
            type Output = Exponent;
            #[inline]
            fn $method(self, other: &Exponent) -> Exponent { (&self).$method(other) }
        }
        impl $trait<Exponent> for &Exponent {
            type Output = Exponent;
            #[inline]
            fn $method(self, other: Exponent) -> Exponent { self.$method(&other) }
        }
    };
}
binop!(Add, add, |a, b| a.add_ref(b));
binop!(Sub, sub, |a, b| a.add_ref(&-b));

macro_rules! intop {
    ($t:ty) => {
        impl Add<$t> for &Exponent {
            type Output = Exponent;
            #[inline]
            fn add(self, other: $t) -> Exponent { self.add_i128(other as i128) }
        }
        impl Add<$t> for Exponent {
            type Output = Exponent;
            #[inline]
            fn add(self, other: $t) -> Exponent { self.add_i128(other as i128) }
        }
        impl Sub<$t> for &Exponent {
            type Output = Exponent;
            #[inline]
            fn sub(self, other: $t) -> Exponent {
                match (other as i128).checked_neg() {
                    Some(n) => self.add_i128(n),
                    None => Exponent::from_bigint(self.to_bigint() - BigInt::from(other)),
                }
            }
        }
        impl Sub<$t> for Exponent {
            type Output = Exponent;
            #[inline]
            fn sub(self, other: $t) -> Exponent { (&self) - other }
        }
        impl AddAssign<$t> for Exponent {
            #[inline]
            fn add_assign(&mut self, other: $t) { *self = &*self + other; }
        }
        impl SubAssign<$t> for Exponent {
            #[inline]
            fn sub_assign(&mut self, other: $t) { *self = &*self - other; }
        }
        impl Mul<$t> for &Exponent {
            type Output = Exponent;
            #[inline]
            fn mul(self, other: $t) -> Exponent {
                if let Exponent::Small(a) = self {
                    if let Some(p) = a.checked_mul(other as i128) { return Exponent::from_i128(p); }
                }
                Exponent::from_bigint(self.to_bigint() * BigInt::from(other))
            }
        }
        impl Mul<$t> for Exponent {
            type Output = Exponent;
            #[inline]
            fn mul(self, other: $t) -> Exponent { (&self) * other }
        }
    };
}
impl AddAssign<&Exponent> for Exponent {
    #[inline]
    fn add_assign(&mut self, other: &Exponent) { *self = &*self + other; }
}
impl AddAssign<Exponent> for Exponent {
    #[inline]
    fn add_assign(&mut self, other: Exponent) { *self = &*self + &other; }
}
impl SubAssign<&Exponent> for Exponent {
    #[inline]
    fn sub_assign(&mut self, other: &Exponent) { *self = &*self - other; }
}
impl SubAssign<Exponent> for Exponent {
    #[inline]
    fn sub_assign(&mut self, other: Exponent) { *self = &*self - &other; }
}

intop!(i128);
intop!(i64);
intop!(i32);
intop!(u32);
intop!(u64);
intop!(usize);

impl Neg for &Exponent {
    type Output = Exponent;
    #[inline]
    fn neg(self) -> Exponent {
        match self {
            Exponent::Small(v) => Exponent::Small(-v),
            Exponent::Big(v) => Exponent::from_bigint(-(**v).clone()),
        }
    }
}
impl Neg for Exponent {
    type Output = Exponent;
    #[inline]
    fn neg(self) -> Exponent { -&self }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_arithmetic_crosses_the_inline_range() {
        let limit = Exponent::from(SMALL_LIMIT);
        assert!(matches!(limit, Exponent::Small(_)));
        let past = &limit + 1i64;
        assert!(matches!(past, Exponent::Big(_)));
        assert_eq!(&past - 1i64, limit);
        assert!(matches!(&past - 1i64, Exponent::Small(_)));
        assert!(past > SMALL_LIMIT);
        assert!(-&past < -SMALL_LIMIT);
        assert!(past > limit && -&past < -&limit);
        assert_eq!(&past + &-&past, Exponent::ZERO);
        assert_eq!(Exponent::from(i128::MAX) + i128::MAX, Exponent::from(BigInt::from(i128::MAX) * 2));
        assert_eq!(Exponent::from(BigInt::from(-7)).half_floor(), Exponent::from(-4i64));
        assert_eq!((&past + &past).half_floor(), past);
    }
}
