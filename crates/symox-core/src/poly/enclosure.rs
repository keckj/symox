//! Certified enclosures of algebraic numbers.
//!
//! An [`Iv`] is a real interval `[lo*2^e, hi*2^e]` with integer endpoints
//! and a shared binary exponent; every operation rounds outward, so the
//! exact value of an expression always lies in the computed rectangle.
//! A [`CBox`] is a complex rectangle. An exact zero imaginary part stays
//! exact through real operations, which decides the branch of real radicals.
//!
//! Enclosures cover rationals, `I`, the golden ratio and the tribonacci
//! constant, sums, products, integer and rational powers (principal
//! branch), `exp(I*pi*r)`, `cos`, `sin` and `tan` of rational multiples of
//! `pi`, and indexed polynomial roots. Transcendental constants use the
//! correctly rounded binary floating point values with a margin far larger
//! than their rounding error.

use std::cell::RefCell;
use std::collections::HashMap;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use super::int::{Int, Rat};
use super::rootisolation as iso;
use crate::bigfloat::{self, BigFloat};
use crate::expr::{Constant, Expr, Func, Kind};

/// Why an enclosure could not be computed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EncError {
    /// The working precision cannot separate the value from a branch cut or
    /// from zero; a larger precision may succeed.
    Precision,
    /// The expression is not built from supported algebraic operations.
    NotAlgebraic(Expr),
}

fn shr_floor(m: &BigInt, k: u64) -> BigInt {
    if m.is_negative() {
        -((-m - 1u32) >> k) - 1
    } else {
        m >> k
    }
}
fn shr_ceil(m: &BigInt, k: u64) -> BigInt {
    -shr_floor(&-m, k)
}
/// `m * 2^k`, rounded down when `k < 0`.
/// A float exponent as a dyadic interval exponent. Isolation works on
/// values whose binary exponents fit in 64 bits.
fn exp64(e: &crate::bigfloat::Exponent) -> i64 {
    e.to_i64().expect("binary exponent outside the root enclosure range")
}

fn shift_floor(m: &BigInt, k: i64) -> BigInt {
    if k >= 0 {
        m << (k as u64)
    } else {
        shr_floor(m, (-k) as u64)
    }
}
fn shift_ceil(m: &BigInt, k: i64) -> BigInt {
    if k >= 0 {
        m << (k as u64)
    } else {
        shr_ceil(m, (-k) as u64)
    }
}
fn div_ceil(a: &BigInt, b: &BigInt) -> BigInt {
    -((-a).div_floor(b))
}

/// A real interval `[lo*2^e, hi*2^e]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Iv {
    pub lo: BigInt,
    pub hi: BigInt,
    pub e: i64,
}

impl Iv {
    pub fn zero() -> Iv {
        Iv { lo: BigInt::zero(), hi: BigInt::zero(), e: 0 }
    }
    pub fn int(n: BigInt) -> Iv {
        Iv { lo: n.clone(), hi: n, e: 0 }
    }
    pub fn point(m: BigInt, e: i64) -> Iv {
        Iv { lo: m.clone(), hi: m, e }
    }
    pub fn from_bigfloat(x: &BigFloat) -> Iv {
        Iv::point(x.man.clone(), exp64(&x.exp))
    }
    /// `[c - 2^k, c + 2^k]`.
    pub fn ball(c: &BigFloat, k: i64) -> Iv {
        let ce = exp64(&c.exp);
        let e = ce.min(k);
        let m = shift_floor(&c.man, ce - e);
        let r = BigInt::one() << ((k - e) as u64);
        Iv { lo: &m - &r, hi: &m + &r, e }
    }
    pub fn from_rat(r: &Rat, prec: u32) -> Iv {
        let n = r.n.to_big();
        let d = r.d.to_big();
        if d.is_one() {
            return Iv::int(n);
        }
        let tz = d.trailing_zeros().unwrap_or(0);
        if d == (BigInt::one() << tz) {
            return Iv::point(n, -(tz as i64));
        }
        let k = prec as u64 + d.bits() + 2;
        let s = &n << k;
        Iv { lo: s.div_floor(&d), hi: div_ceil(&s, &d), e: -(k as i64) }
    }
    pub fn from_int(n: &Int) -> Iv {
        Iv::int(n.to_big())
    }
    pub fn is_exact_zero(&self) -> bool {
        self.lo.is_zero() && self.hi.is_zero()
    }
    pub fn is_point(&self) -> bool {
        self.lo == self.hi
    }
    pub fn contains_zero(&self) -> bool {
        !self.lo.is_positive() && !self.hi.is_negative()
    }
    pub fn positive(&self) -> bool {
        self.lo.is_positive()
    }
    pub fn negative(&self) -> bool {
        self.hi.is_negative()
    }
    pub fn round(mut self, prec: u32) -> Iv {
        let bits = self.lo.bits().max(self.hi.bits());
        if bits > prec as u64 + 2 {
            let k = bits - prec as u64 - 2;
            self.lo = shr_floor(&self.lo, k);
            self.hi = shr_ceil(&self.hi, k);
            self.e += k as i64;
        }
        self
    }
    fn aligned(&self, o: &Iv) -> (BigInt, BigInt, BigInt, BigInt, i64) {
        let e = self.e.min(o.e);
        (
            &self.lo << ((self.e - e) as u64),
            &self.hi << ((self.e - e) as u64),
            &o.lo << ((o.e - e) as u64),
            &o.hi << ((o.e - e) as u64),
            e,
        )
    }
    pub fn add(&self, o: &Iv, prec: u32) -> Iv {
        if self.is_exact_zero() {
            return o.clone();
        }
        if o.is_exact_zero() {
            return self.clone();
        }
        let (a, b, c, d, e) = self.aligned(o);
        Iv { lo: a + c, hi: b + d, e }.round(prec)
    }
    pub fn neg(&self) -> Iv {
        Iv { lo: -&self.hi, hi: -&self.lo, e: self.e }
    }
    pub fn sub(&self, o: &Iv, prec: u32) -> Iv {
        self.add(&o.neg(), prec)
    }
    pub fn mul(&self, o: &Iv, prec: u32) -> Iv {
        if self.is_exact_zero() || o.is_exact_zero() {
            return Iv::zero();
        }
        let p = [&self.lo * &o.lo, &self.lo * &o.hi, &self.hi * &o.lo, &self.hi * &o.hi];
        let lo = p.iter().min().unwrap().clone();
        let hi = p.iter().max().unwrap().clone();
        Iv { lo, hi, e: self.e + o.e }.round(prec)
    }
    pub fn sqr(&self, prec: u32) -> Iv {
        let (a, b) = (&self.lo * &self.lo, &self.hi * &self.hi);
        let (lo, hi) = if self.lo.is_negative() && self.hi.is_positive() {
            (BigInt::zero(), a.max(b))
        } else if a < b {
            (a, b)
        } else {
            (b, a)
        };
        Iv { lo, hi, e: 2 * self.e }.round(prec)
    }
    pub fn mul_rat(&self, r: &Rat, prec: u32) -> Iv {
        self.mul(&Iv::from_rat(r, prec), prec)
    }
    pub fn inv(&self, prec: u32) -> Option<Iv> {
        if self.contains_zero() {
            return None;
        }
        if self.is_point() && self.lo.magnitude().is_one() {
            return Some(Iv::point(self.lo.clone(), -self.e));
        }
        let k = prec as u64 + self.hi.bits().max(self.lo.bits()) + 2;
        let n = BigInt::one() << k;
        Some(
            Iv { lo: n.div_floor(&self.hi), hi: div_ceil(&n, &self.lo), e: -self.e - k as i64 }
                .round(prec),
        )
    }
    /// The smallest interval containing both.
    pub fn hull(&self, o: &Iv) -> Iv {
        let (a, b, c, d, e) = self.aligned(o);
        Iv { lo: a.min(c), hi: b.max(d), e }
    }
    /// Upper bound of `|x|` on the interval, as `(m, e)` with value `m*2^e`.
    pub fn mag(&self) -> (BigInt, i64) {
        (BigInt::from(self.lo.magnitude().max(self.hi.magnitude()).clone()), self.e)
    }
    pub fn lo_rat(&self) -> Rat {
        dyadic_rat(&self.lo, self.e)
    }
    pub fn hi_rat(&self) -> Rat {
        dyadic_rat(&self.hi, self.e)
    }
    pub fn mid(&self, prec: u32) -> BigFloat {
        BigFloat::round(&self.lo + &self.hi, i128::from(self.e - 1), prec)
    }
    /// Whether the width is at most `2^k`.
    pub fn width_below(&self, k: i64) -> bool {
        let w = &self.hi - &self.lo;
        w.is_zero() || (w.bits() as i64 + self.e) <= k
    }
    /// The binary exponent of the magnitude, `None` for an interval that
    /// contains zero.
    pub fn magnitude_exponent(&self) -> Option<i64> {
        if self.contains_zero() {
            None
        } else {
            Some(self.lo.magnitude().min(self.hi.magnitude()).bits() as i64 + self.e)
        }
    }
    pub fn within(&self, lo: &Rat, hi: &Rat) -> bool {
        &self.lo_rat() >= lo && &self.hi_rat() <= hi
    }
}

fn dyadic_rat(m: &BigInt, e: i64) -> Rat {
    if e >= 0 {
        Rat::new(Int::from_big(m << (e as u64)), Int::one())
    } else {
        Rat::new(Int::from_big(m.clone()), Int::from_big(BigInt::one() << ((-e) as u64)))
    }
}

/// A complex rectangle.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CBox {
    pub re: Iv,
    pub im: Iv,
}

impl CBox {
    pub fn real(re: Iv) -> CBox {
        CBox { re, im: Iv::zero() }
    }
    pub fn zero() -> CBox {
        CBox::real(Iv::zero())
    }
    pub fn one() -> CBox {
        CBox::real(Iv::int(BigInt::one()))
    }
    pub fn i() -> CBox {
        CBox { re: Iv::zero(), im: Iv::int(BigInt::one()) }
    }
    pub fn from_rat(r: &Rat, prec: u32) -> CBox {
        CBox::real(Iv::from_rat(r, prec))
    }
    pub fn is_real(&self) -> bool {
        self.im.is_exact_zero()
    }
    pub fn is_exact_zero(&self) -> bool {
        self.re.is_exact_zero() && self.im.is_exact_zero()
    }
    pub fn contains_zero(&self) -> bool {
        self.re.contains_zero() && self.im.contains_zero()
    }
    pub fn add(&self, o: &CBox, prec: u32) -> CBox {
        CBox { re: self.re.add(&o.re, prec), im: self.im.add(&o.im, prec) }
    }
    pub fn sub(&self, o: &CBox, prec: u32) -> CBox {
        CBox { re: self.re.sub(&o.re, prec), im: self.im.sub(&o.im, prec) }
    }
    pub fn neg(&self) -> CBox {
        CBox { re: self.re.neg(), im: self.im.neg() }
    }
    pub fn mul(&self, o: &CBox, prec: u32) -> CBox {
        let re = self.re.mul(&o.re, prec).sub(&self.im.mul(&o.im, prec), prec);
        let im = self.re.mul(&o.im, prec).add(&self.im.mul(&o.re, prec), prec);
        CBox { re, im }
    }
    pub fn sqr(&self, prec: u32) -> CBox {
        let re = self.re.sqr(prec).sub(&self.im.sqr(prec), prec);
        let p = self.re.mul(&self.im, prec);
        CBox { re, im: p.add(&p, prec) }
    }
    pub fn scale(&self, r: &Iv, prec: u32) -> CBox {
        CBox { re: self.re.mul(r, prec), im: self.im.mul(r, prec) }
    }
    /// Enclosure of `|z|^2`.
    pub fn norm2(&self, prec: u32) -> Iv {
        self.re.sqr(prec).add(&self.im.sqr(prec), prec)
    }
    pub fn inv(&self, prec: u32) -> Option<CBox> {
        if self.is_real() {
            return Some(CBox::real(self.re.inv(prec)?));
        }
        if self.re.is_exact_zero() {
            return Some(CBox { re: Iv::zero(), im: self.im.inv(prec)?.neg() });
        }
        let d = self.norm2(prec).inv(prec)?;
        Some(CBox { re: self.re.mul(&d, prec), im: self.im.neg().mul(&d, prec) })
    }
    pub fn powi(&self, n: i64, prec: u32) -> Option<CBox> {
        let mut base = if n < 0 { self.inv(prec)? } else { self.clone() };
        let mut k = n.unsigned_abs();
        let mut result = CBox::one();
        while k != 0 {
            if k & 1 == 1 {
                result = result.mul(&base, prec);
            }
            k >>= 1;
            if k != 0 {
                base = base.sqr(prec);
            }
        }
        Some(result)
    }
    pub fn round(self, prec: u32) -> CBox {
        CBox { re: self.re.round(prec), im: self.im.round(prec) }
    }
    pub fn hull(&self, o: &CBox) -> CBox {
        CBox { re: self.re.hull(&o.re), im: self.im.hull(&o.im) }
    }
    /// Whether the widths are at most `2^k`.
    pub fn width_below(&self, k: i64) -> bool {
        self.re.width_below(k) && self.im.width_below(k)
    }
    /// A rectangle that certainly excludes zero.
    pub fn excludes_zero(&self) -> bool {
        !self.contains_zero()
    }
    /// Whether the two rectangles certainly do not intersect.
    pub fn disjoint(&self, o: &CBox) -> bool {
        self.re.hi_rat() < o.re.lo_rat()
            || o.re.hi_rat() < self.re.lo_rat()
            || self.im.hi_rat() < o.im.lo_rat()
            || o.im.hi_rat() < self.im.lo_rat()
    }
    pub fn mid(&self, prec: u32) -> (BigFloat, BigFloat) {
        (self.re.mid(prec), self.im.mid(prec))
    }
}

/// Horner evaluation of a polynomial with interval coefficients, lowest
/// degree first.
pub fn horner(coefficients: &[CBox], z: &CBox, prec: u32) -> CBox {
    let mut value = CBox::zero();
    for c in coefficients.iter().rev() {
        value = value.mul(z, prec).add(c, prec);
    }
    value
}

/// Horner evaluation of a rational polynomial, lowest degree first.
pub fn horner_rat(coefficients: &[Rat], z: &CBox, prec: u32) -> CBox {
    let mut value = CBox::zero();
    for c in coefficients.iter().rev() {
        value = value.mul(z, prec).add(&CBox::from_rat(c, prec), prec);
    }
    value
}

// ---------------------------------------------------------------------------
// Transcendental constants with margins.

fn margin_exponent(prec: u32) -> i64 {
    -(prec as i64) - 4
}

/// Enclosures of `cos(pi*r)` and `sin(pi*r)`.
pub fn cos_sin_pi(r: &Rat, prec: u32) -> (Iv, Iv) {
    // Reduce r modulo 2 exactly.
    let two = Rat::from_i64(2);
    let k = r.div(&two).floor();
    let r = r.sub(&two.mul(&Rat::new(k, Int::one())));
    let quarter = |n: i64| Rat::new(Int::from(n), Int::from(2));
    let one = || Iv::int(BigInt::one());
    let m1 = || Iv::int(-BigInt::one());
    if r.is_zero() {
        return (one(), Iv::zero());
    }
    if r == quarter(1) {
        return (Iv::zero(), one());
    }
    if r == Rat::one() {
        return (m1(), Iv::zero());
    }
    if r == quarter(3) {
        return (Iv::zero(), m1());
    }
    let wp = prec + 40 + r.n.bits() as u32;
    let p = bigfloat::pi(wp);
    let x = BigFloat::from_int(&r.n.to_big())
        .mul(&p, wp)
        .div(&BigFloat::from_int(&r.d.to_big()), wp)
        .expect("nonzero denominator");
    let (s, c) = x.sin_cos(wp);
    let m = margin_exponent(prec);
    (Iv::ball(&c, m).round(prec), Iv::ball(&s, m).round(prec))
}

// ---------------------------------------------------------------------------
// Radicals.

/// Enclosure of `x^(1/q)` for a positive interval `x`.
fn real_root(x: &Iv, q: u32, prec: u32) -> Iv {
    if x.is_point() {
        if let Some(e) = (x.e % q as i64 == 0).then_some(x.e / q as i64) {
            let r = x.lo.nth_root(q);
            if r.pow(q) == x.lo {
                return Iv::point(r, e);
            }
        }
    }
    let target = q as i64 * (prec as i64 + 8);
    let s0 = target - x.hi.bits() as i64;
    let s = s0 + (x.e - s0).rem_euclid(q as i64);
    let lo = shift_floor(&x.lo, s);
    let hi = shift_ceil(&x.hi, s);
    let e = (x.e - s) / q as i64;
    let a = if lo.is_positive() { lo.nth_root(q) } else { BigInt::zero() };
    let mut b = hi.nth_root(q);
    if b.pow(q) < hi {
        b += 1u32;
    }
    Iv { lo: a, hi: b, e }.round(prec)
}

/// The principal value of `z^(1/q)`.
pub fn principal_root(z: &CBox, q: u32, prec: u32) -> Result<CBox, EncError> {
    if q == 1 || z.is_exact_zero() {
        return Ok(z.clone());
    }
    if z.is_real() {
        if z.re.positive() {
            return Ok(CBox::real(real_root(&z.re, q, prec)));
        }
        if z.re.negative() {
            let r = real_root(&z.re.neg(), q, prec);
            if q == 2 {
                return Ok(CBox { re: Iv::zero(), im: r });
            }
            let (c, s) = cos_sin_pi(&Rat::new(Int::one(), Int::from(q as i64)), prec);
            return Ok(CBox { re: r.mul(&c, prec), im: r.mul(&s, prec) });
        }
        return Err(EncError::Precision);
    }
    if !(z.re.positive() || !z.im.contains_zero()) {
        // The rectangle meets the negative real axis, the branch cut.
        return Err(EncError::Precision);
    }
    if z.re.is_exact_zero() && q == 2 {
        // sqrt(i*y) = sqrt(|y|/2)*(1 + sign(y)*i).
        let half = z.im.mul(&Iv::point(BigInt::one(), -1), prec);
        let positive = z.im.positive();
        let r = real_root(&if positive { half } else { half.neg() }, 2, prec);
        return Ok(CBox { re: r.clone(), im: if positive { r } else { r.neg() } });
    }
    // Approximate the principal root at the center, then certify a disc
    // around it: some root of w^q - z0 lies within |w^q - z0|/|w|^(q-1) of
    // w, and the disc must lie in the open sector |arg| < pi/q, which holds
    // exactly one root for every z0 off the branch cut.
    let wp = prec + 40;
    let (a, b) = z.mid(wp);
    let modulus = a.mul(&a, wp).add(&b.mul(&b, wp), wp).sqrt(wp).ok_or(EncError::Precision)?;
    if modulus.is_zero() {
        return Err(EncError::Precision);
    }
    let angle = BigFloat::atan2(&b, &a, wp).div_i64(q as i64, wp);
    let rho = modulus.ln(wp).ok_or(EncError::Precision)?.div_i64(q as i64, wp).exp(wp);
    let (s, c) = angle.sin_cos(wp);
    let w = CBox {
        re: Iv::from_bigfloat(&rho.mul(&c, wp)),
        im: Iv::from_bigfloat(&rho.mul(&s, wp)),
    };
    let wp2 = wp + 16;
    let defect = w.powi(q as i64, wp2).ok_or(EncError::Precision)?.sub(z, wp2).norm2(wp2);
    let base = w.norm2(wp2);
    let denominator = base.powi_real(q as i64 - 1, wp2);
    let ratio = defect.mul(&denominator.inv(wp2).ok_or(EncError::Precision)?, wp2);
    let radius = sqrt_upper(&ratio);
    let disc = CBox {
        re: w.re.add(&radius_iv(&radius), wp2),
        im: w.im.add(&radius_iv(&radius), wp2),
    };
    let (cos_q, sin_q) = cos_sin_pi(&Rat::new(Int::one(), Int::from(q as i64)), wp2);
    // Im(u*e^{i*phi}) > 0 and Im(u*e^{-i*phi}) < 0 with phi = pi/q.
    let upper = disc.re.mul(&sin_q, wp2).add(&disc.im.mul(&cos_q, wp2), wp2);
    let lower = disc.im.mul(&cos_q, wp2).sub(&disc.re.mul(&sin_q, wp2), wp2);
    if !(upper.positive() && lower.negative()) {
        return Err(EncError::Precision);
    }
    Ok(disc.round(prec))
}

impl Iv {
    fn powi_real(&self, n: i64, prec: u32) -> Iv {
        let mut result = Iv::int(BigInt::one());
        for _ in 0..n {
            result = result.mul(self, prec);
        }
        result
    }
}

/// An upper bound of the square root of the upper end of a nonnegative
/// interval, as `(m, e)` with value `m*2^e`.
fn sqrt_upper(x: &Iv) -> (BigInt, i64) {
    let mut hi = x.hi.clone();
    let mut e = x.e;
    if hi.is_negative() {
        return (BigInt::zero(), 0);
    }
    if e % 2 != 0 {
        hi <<= 1u32;
        e -= 1;
    }
    // Keep enough bits for a tight root.
    let extra = 64u64;
    hi <<= 2 * extra;
    e -= 2 * extra as i64;
    let mut r = hi.sqrt();
    if &r * &r < hi {
        r += 1u32;
    }
    (r, e / 2)
}
fn radius_iv(r: &(BigInt, i64)) -> Iv {
    Iv { lo: -&r.0, hi: r.0.clone(), e: r.1 }
}

// ---------------------------------------------------------------------------
// Indexed roots of integer polynomials.

fn sign_at_rat(f: &[BigInt], r: &Rat) -> i32 {
    let (p, q) = (r.n.to_big(), r.d.to_big());
    // sum c_i p^(n-i) q^i
    let mut acc = BigInt::zero();
    let mut qpow = BigInt::one();
    let n = f.len() - 1;
    let mut terms = Vec::with_capacity(n + 1);
    for _ in 0..=n {
        terms.push(qpow.clone());
        qpow *= &q;
    }
    for (i, c) in f.iter().enumerate() {
        acc = acc * &p + c * &terms[i];
    }
    match acc.sign() {
        num_bigint::Sign::Minus => -1,
        num_bigint::Sign::NoSign => 0,
        num_bigint::Sign::Plus => 1,
    }
}

#[derive(Clone)]
enum RootState {
    Real { f: Vec<BigInt>, lo: Rat, hi: Rat, exact: bool },
    Complex { f: Vec<Rat>, isolating: iso::ComplexBox },
}

thread_local! {
    static ROOTS: RefCell<HashMap<(Vec<Int>, usize), (RootState, Option<(u32, CBox)>)>> = RefCell::new(HashMap::new());
}

fn root_state(f: &[Int], index: usize) -> Option<RootState> {
    let (reals, g) = iso::isolate_real_roots_sqf(f, None, None, None, false).ok()?;
    if index < reals.len() {
        let (s, t) = reals[index].clone();
        let g: Vec<BigInt> = g.iter().map(Int::to_big).collect();
        let exact = s == t;
        return Some(RootState::Real { f: g, lo: s, hi: t, exact });
    }
    let rats: Vec<Rat> = f.iter().map(|c| Rat::new(c.clone(), Int::one())).collect();
    let boxes = iso::isolate_complex_roots_sqf(&rats, None, None, None).ok()?;
    let b = boxes.into_iter().nth(index - reals.len())?;
    Some(RootState::Complex { f: rats, isolating: b })
}

/// Enclosure of the real root of `f` (leading first) in the isolating
/// interval, to absolute width `2^-prec` relative to its magnitude.
fn refine_real(f: &[BigInt], lo: &mut Rat, hi: &mut Rat, prec: u32) -> Iv {
    let mut slo = sign_at_rat(f, lo);
    if slo == 0 {
        *hi = lo.clone();
        return Iv::from_rat(lo, prec);
    }
    let shi = sign_at_rat(f, hi);
    if shi == 0 {
        *lo = hi.clone();
        return Iv::from_rat(hi, prec);
    }
    let scale = {
        let a = lo.abs();
        let b = hi.abs();
        let m = if a < b { b } else { a };
        if m < Rat::one() {
            Rat::one()
        } else {
            m
        }
    };
    let tolerance = scale.div(&Rat::new(Int::from_big(BigInt::one() << (prec as u64 + 2)), Int::one()));
    // Newton steps from the midpoint, accepted only with a sign change.
    let n = f.len() - 1;
    let df: Vec<BigInt> = f[..n].iter().enumerate().map(|(i, c)| c * BigInt::from((n - i) as u64)).collect();
    let fq: Vec<Rat> = f.iter().map(|c| Rat::new(Int::from_big(c.clone()), Int::one())).collect();
    let dfq: Vec<Rat> = df.iter().map(|c| Rat::new(Int::from_big(c.clone()), Int::one())).collect();
    let mut guard = 0;
    while hi.sub(lo) > tolerance {
        guard += 1;
        if guard > 100_000 {
            break;
        }
        // Try a Newton step at working precision once the interval is small.
        let width = hi.sub(lo);
        if width < scale.div(&Rat::from_i64(1 << 20)) && !df.is_empty() {
            let wp = prec + 32;
            let mid = lo.add(hi).div(&Rat::from_i64(2));
            let mut x = BigFloat::from_ratio(&mid.n.to_big(), &mid.d.to_big(), wp);
            let mut ok = true;
            for _ in 0..64 {
                let (p, d) = bf_horner(&fq, &dfq, &x, wp);
                if d.is_zero() {
                    ok = false;
                    break;
                }
                let step = p.div(&d, wp).unwrap();
                x = x.sub(&step, wp);
                if step.is_zero() || step.mag() < x.mag().max(0) - wp as i128 + 4 {
                    break;
                }
            }
            if ok {
                let k = -(prec as i64) - 2 + scale_exponent(&scale);
                let xe = exp64(&x.exp);
                let e = xe.min(k);
                let m = shift_floor(&x.man, xe - e);
                let r = BigInt::one() << ((k - e) as u64);
                let a = dyadic_rat(&(&m - &r), e);
                let b = dyadic_rat(&(&m + &r), e);
                if &a >= lo && &b <= hi {
                    let sa = sign_at_rat(f, &a);
                    let sb = sign_at_rat(f, &b);
                    if sa == 0 {
                        *lo = a.clone();
                        *hi = a;
                        break;
                    }
                    if sb == 0 {
                        *lo = b.clone();
                        *hi = b;
                        break;
                    }
                    if sa != sb {
                        *lo = a;
                        *hi = b;
                        break;
                    }
                }
            }
        }
        let mid = lo.add(hi).div(&Rat::from_i64(2));
        let sm = sign_at_rat(f, &mid);
        if sm == 0 {
            *lo = mid.clone();
            *hi = mid;
            break;
        }
        if sm == slo {
            *lo = mid;
            slo = sm;
        } else {
            *hi = mid;
        }
    }
    let a = Iv::from_rat(lo, prec + 8);
    let b = Iv::from_rat(hi, prec + 8);
    Iv { lo: a.lo.clone(), hi: a.lo.clone(), e: a.e }.hull(&Iv { lo: b.hi.clone(), hi: b.hi, e: b.e }).round(prec)
}

fn scale_exponent(r: &Rat) -> i64 {
    r.n.bits() as i64 - r.d.bits() as i64 + 1
}

fn bf_horner(f: &[Rat], df: &[Rat], x: &BigFloat, wp: u32) -> (BigFloat, BigFloat) {
    let conv = |c: &Rat| BigFloat::from_ratio(&c.n.to_big(), &c.d.to_big(), wp);
    let mut p = BigFloat::zero();
    for c in f {
        p = p.mul(x, wp).add(&conv(c), wp);
    }
    let mut d = BigFloat::zero();
    for c in df {
        d = d.mul(x, wp).add(&conv(c), wp);
    }
    (p, d)
}

/// Complex arithmetic on pairs of binary floats for Newton iteration.
fn c_horner(f: &[Rat], z: &(BigFloat, BigFloat), wp: u32) -> ((BigFloat, BigFloat), (BigFloat, BigFloat)) {
    let conv = |c: &Rat| BigFloat::from_ratio(&c.n.to_big(), &c.d.to_big(), wp);
    let mul = |a: &(BigFloat, BigFloat), b: &(BigFloat, BigFloat)| {
        (
            a.0.mul(&b.0, wp).sub(&a.1.mul(&b.1, wp), wp),
            a.0.mul(&b.1, wp).add(&a.1.mul(&b.0, wp), wp),
        )
    };
    let mut p = (BigFloat::zero(), BigFloat::zero());
    let mut d = (BigFloat::zero(), BigFloat::zero());
    for c in f {
        d = mul(&d, z);
        d = (d.0.add(&p.0, wp), d.1.add(&p.1, wp));
        p = mul(&p, z);
        p = (p.0.add(&conv(c), wp), p.1);
    }
    (p, d)
}

/// Whether the closed rectangle spanned by `disc` and the isolating box
/// holds a single root. The disc holds a root and the box holds the
/// indexed one, so they are then the same root. This accepts a root on an
/// edge of its isolating box, where no disc fits inside the box.
fn single_root_in_hull(f: &[Rat], disc: &CBox, (ax, bx, ay, by): (&Rat, &Rat, &Rat, &Rat)) -> bool {
    let min = |a: Rat, b: &Rat| if &a < b { a } else { b.clone() };
    let max = |a: Rat, b: &Rat| if &a > b { a } else { b.clone() };
    let inf = (min(disc.re.lo_rat(), ax), min(disc.im.lo_rat(), ay));
    let sup = (max(disc.re.hi_rat(), bx), max(disc.im.hi_rat(), by));
    iso::count_complex_roots(f, Some(inf), Some(sup), &iso::Exclude::None).is_ok_and(|n| n == 1)
}

fn refine_complex(f: &[Rat], isolating: &mut iso::ComplexBox, prec: u32) -> Option<CBox> {
    let n = f.len() - 1;
    let lowest_first: Vec<Rat> = f.iter().rev().cloned().collect();
    let derivative: Vec<Rat> = lowest_first
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c.mul(&Rat::from_i64(i as i64)))
        .collect();
    for _attempt in 0..64 {
        let wp = prec + 32;
        let (ax, bx, ay, by) = (isolating.ax(), isolating.bx(), isolating.ay(), isolating.by());
        let two = Rat::from_i64(2);
        let cx = ax.add(&bx).div(&two);
        let cy = ay.add(&by).div(&two);
        let mut z = (
            BigFloat::from_ratio(&cx.n.to_big(), &cx.d.to_big(), wp),
            BigFloat::from_ratio(&cy.n.to_big(), &cy.d.to_big(), wp),
        );
        let mut converged = false;
        for _ in 0..200 {
            let (p, d) = c_horner(f, &z, wp);
            let dn = d.0.mul(&d.0, wp).add(&d.1.mul(&d.1, wp), wp);
            if dn.is_zero() {
                break;
            }
            // step = p/d = p*conj(d)/|d|^2
            let sr = p.0.mul(&d.0, wp).add(&p.1.mul(&d.1, wp), wp).div(&dn, wp).unwrap();
            let si = p.1.mul(&d.0, wp).sub(&p.0.mul(&d.1, wp), wp).div(&dn, wp).unwrap();
            z = (z.0.sub(&sr, wp), z.1.sub(&si, wp));
            let size = z.0.mag().max(z.1.mag()).max(0);
            let step = if sr.is_zero() { i128::MIN } else { sr.mag() }
                .max(if si.is_zero() { i128::MIN } else { si.mag() });
            if step < size - wp as i128 + 8 {
                converged = true;
                break;
            }
        }
        if converged {
            let w = CBox { re: Iv::from_bigfloat(&z.0), im: Iv::from_bigfloat(&z.1) };
            let wp2 = wp + 16;
            let value = horner_rat(&lowest_first, &w, wp2).norm2(wp2);
            let slope = horner_rat(&derivative, &w, wp2).norm2(wp2);
            if let Some(inverse) = slope.inv(wp2) {
                let ratio = value.mul(&inverse, wp2).mul(&Iv::int(BigInt::from((n * n) as u64)), wp2);
                let radius = sqrt_upper(&ratio);
                let disc = CBox {
                    re: w.re.add(&radius_iv(&radius), wp2),
                    im: w.im.add(&radius_iv(&radius), wp2),
                };
                if (disc.re.within(&ax, &bx) && disc.im.within(&ay, &by))
                    || single_root_in_hull(f, &disc, (&ax, &bx, &ay, &by))
                {
                    let size = disc.re.mag().0.bits() as i64 + disc.re.mag().1;
                    let size = size.max(disc.im.mag().0.bits() as i64 + disc.im.mag().1).max(0);
                    if disc.width_below(size - prec as i64) {
                        return Some(disc.round(prec + 8));
                    }
                }
            }
        }
        for _ in 0..4 {
            *isolating = isolating.refine().ok()?;
        }
    }
    None
}

/// Enclosure of the indexed root of the integer polynomial `f` (leading
/// coefficient first).
pub fn root_enclosure(f: &[Int], index: usize, prec: u32) -> Option<CBox> {
    let key = (f.to_vec(), index);
    let cached = ROOTS.with(|r| r.borrow().get(&key).cloned());
    let (mut state, best) = match cached {
        Some(v) => v,
        None => (root_state(f, index)?, None),
    };
    if let Some((p, b)) = &best {
        if *p >= prec {
            return Some(b.clone());
        }
    }
    let enclosure = match &mut state {
        RootState::Real { f, lo, hi, exact } => {
            if *exact {
                CBox::real(Iv::from_rat(lo, prec))
            } else {
                CBox::real(refine_real(f, lo, hi, prec))
            }
        }
        RootState::Complex { f, isolating } => refine_complex(f, isolating, prec)?,
    };
    ROOTS.with(|r| r.borrow_mut().insert(key, (state, Some((prec, enclosure.clone())))));
    Some(enclosure)
}

// ---------------------------------------------------------------------------
// Expressions.

fn rational(e: &Expr) -> Option<Rat> {
    e.as_number().and_then(Rat::from_number)
}

/// The rational `r` with `e == r*pi`, if any.
pub fn pi_multiple(e: &Expr) -> Option<Rat> {
    if let Kind::Const(Constant::Pi) = e.kind() {
        return Some(Rat::one());
    }
    if let Kind::Mul(c, factors) = e.kind() {
        if let [(b, x)] = factors.as_slice() {
            if matches!(b.kind(), Kind::Const(Constant::Pi)) && x.is_one() {
                return Rat::from_number(c);
            }
        }
    }
    None
}

/// The rational `r` with `e == I*pi*r`, if any.
pub fn i_pi_multiple(e: &Expr) -> Option<Rat> {
    let Kind::Mul(c, factors) = e.kind() else {
        return None;
    };
    let c = Rat::from_number(c)?;
    let mut pi = false;
    let mut i = false;
    for (b, x) in factors {
        if !x.is_one() {
            return None;
        }
        match b.kind() {
            Kind::Const(Constant::Pi) if !pi => pi = true,
            Kind::Const(Constant::I) if !i => i = true,
            _ => return None,
        }
    }
    (pi && i).then_some(c)
}

/// Enclosure of an algebraic numeric expression at `prec` bits.
pub fn enclose(e: &Expr, prec: u32) -> Result<CBox, EncError> {
    let not = || EncError::NotAlgebraic(e.clone());
    match e.kind() {
        Kind::Num(n) => Rat::from_number(n).map(|r| CBox::from_rat(&r, prec)).ok_or_else(not),
        Kind::Const(Constant::I) => Ok(CBox::i()),
        Kind::Const(Constant::GoldenRatio) => {
            let s = principal_root(&CBox::from_rat(&Rat::from_i64(5), prec), 2, prec)?;
            Ok(s.add(&CBox::one(), prec).scale(&Iv::point(BigInt::one(), -1), prec))
        }
        Kind::Const(Constant::TribonacciConstant) => {
            let f: Vec<Int> = [1, -1, -1, -1].iter().map(|&c| Int::from(c)).collect();
            root_enclosure(&f, 0, prec).ok_or(EncError::Precision)
        }
        Kind::Add(c, terms) => {
            let mut s = CBox::from_rat(&Rat::from_number(c).ok_or_else(not)?, prec);
            for (t, k) in terms {
                let k = Rat::from_number(k).ok_or_else(not)?;
                let v = enclose(t, prec)?;
                s = s.add(&v.scale(&Iv::from_rat(&k, prec), prec), prec);
            }
            Ok(s)
        }
        Kind::Mul(c, factors) => {
            let mut s = CBox::from_rat(&Rat::from_number(c).ok_or_else(not)?, prec);
            for (b, x) in factors {
                s = s.mul(&enclose_pow(b, x, e, prec)?, prec);
            }
            Ok(s)
        }
        Kind::Pow(b, x) => enclose_pow(b, x, e, prec),
        Kind::Func(Func::Exp, args) if args.len() == 1 => {
            let r = i_pi_multiple(&args[0]).ok_or_else(not)?;
            let (c, s) = cos_sin_pi(&r, prec);
            Ok(CBox { re: c, im: s })
        }
        Kind::Func(f @ (Func::Cos | Func::Sin | Func::Tan), args) if args.len() == 1 => {
            let r = pi_multiple(&args[0]).ok_or_else(not)?;
            let (c, s) = cos_sin_pi(&r, prec);
            Ok(CBox::real(match f {
                Func::Cos => c,
                Func::Sin => s,
                _ => s.mul(&c.inv(prec).ok_or(EncError::Precision)?, prec),
            }))
        }
        Kind::Func(Func::Undefined(_), _) if super::rootof::is_rootof(e) => {
            let (expr, index) = super::rootof::parts(e).ok_or_else(not)?;
            let (_, f) = super::rootof::integer_poly(&expr).ok_or_else(not)?;
            root_enclosure(&f, index, prec).ok_or(EncError::Precision)
        }
        _ => Err(not()),
    }
}

fn enclose_pow(b: &Expr, x: &Expr, whole: &Expr, prec: u32) -> Result<CBox, EncError> {
    let r = rational(x).ok_or_else(|| EncError::NotAlgebraic(whole.clone()))?;
    let base = enclose(b, prec)?;
    let q = r.d.as_i64().filter(|&q| q > 0 && q <= u32::MAX as i64).ok_or(EncError::Precision)? as u32;
    let p = r.n.as_i64().ok_or(EncError::Precision)?;
    let root = principal_root(&base, q, prec)?;
    root.powi(p, prec).ok_or(EncError::Precision)
}

/// Enclosure of `e` whose rectangle is at most `2^-bits` wide relative to
/// a unit scale, trying increasing working precisions.
pub fn enclose_to(e: &Expr, bits: u32, max_prec: u32) -> Result<CBox, EncError> {
    let mut prec = bits + 16;
    loop {
        match enclose(e, prec) {
            Ok(b) if b.width_below(-(bits as i64)) => return Ok(b),
            Ok(_) | Err(EncError::Precision) => {}
            Err(other) => return Err(other),
        }
        if prec >= max_prec {
            return Err(EncError::Precision);
        }
        prec = (prec * 2).min(max_prec);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn expr(s: &str) -> Expr {
        crate::parser::parse_str(s).unwrap()
    }
    fn contains(b: &CBox, re: f64, im: f64) -> bool {
        let tol = 1e-12;
        b.re.lo_rat().to_f64() <= re + tol
            && re - tol <= b.re.hi_rat().to_f64()
            && b.im.lo_rat().to_f64() <= im + tol
            && im - tol <= b.im.hi_rat().to_f64()
    }

    #[test]
    fn radicals_follow_the_principal_branch() {
        for (s, re, im) in [
            ("2^(1/4)", 2f64.powf(0.25), 0.0),
            ("(-2)^(1/3)", 2f64.powf(1.0 / 3.0) * 0.5, 2f64.powf(1.0 / 3.0) * 3f64.sqrt() / 2.0),
            ("sqrt(2)+sqrt(3)", 2f64.sqrt() + 3f64.sqrt(), 0.0),
            ("(1+I)^(1/3)", 1.0842150814913512, 0.2905145555072514),
            ("(-1-I)^(1/2)", 0.45508986056222733, -1.09868411346781),
        ] {
            let b = enclose(&expr(s), 128).unwrap();
            assert!(contains(&b, re, im), "{s}: {b:?}");
            assert!(b.width_below(-100), "{s}");
        }
        assert!(enclose(&expr("2^(1/3)"), 64).unwrap().is_real());
        assert!(enclose(&expr("sqrt(-3)"), 64).unwrap().re.is_exact_zero());
    }

    #[test]
    fn indexed_roots_are_enclosed() {
        let f: Vec<Int> = [1, 0, 0, 0, -1, -1].iter().map(|&c| Int::from(c)).collect();
        let r = root_enclosure(&f, 0, 200).unwrap();
        assert!(contains(&r, 1.1673039782614187, 0.0));
        assert!(r.width_below(-190));
        let c = root_enclosure(&f, 1, 200).unwrap();
        assert!(contains(&c, -0.7648844336005847, -0.35247154603172626));
        let d = root_enclosure(&f, 2, 300).unwrap();
        assert!(contains(&d, -0.7648844336005847, 0.35247154603172626));
        assert!(d.width_below(-290));
    }

    #[test]
    fn trigonometric_values_of_rational_angles() {
        let b = enclose(&expr("cos(pi/7)"), 100).unwrap();
        assert!(contains(&b, (std::f64::consts::PI / 7.0).cos(), 0.0));
        let b = enclose(&expr("exp(2*I*pi/5)"), 100).unwrap();
        let t = 2.0 * std::f64::consts::PI / 5.0;
        assert!(contains(&b, t.cos(), t.sin()));
    }
}
