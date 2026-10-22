//! Arbitrary precision binary floating point (the `mpmath` replacement).
//!
//! A [`BigFloat`] is `man * 2**exp` with an arbitrary precision integer
//! mantissa. Every operation takes the target precision in bits and rounds
//! the result to nearest (ties to even). Elementary functions use argument
//! reduction followed by Taylor series evaluated with a few guard bits;
//! `pi` uses the Chudnovsky series with binary splitting.

use std::cmp::Ordering;

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

mod decimal;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BigFloat {
    /// Signed mantissa; zero is represented by `man == 0, exp == 0`.
    pub man: BigInt,
    pub exp: i128,
}

/// Convert decimal digits to bits of precision (with no guard bits).
pub fn dps_to_prec(dps: u32) -> u32 {
    (((dps as f64 + 1.0) * std::f64::consts::LOG2_10).round() as u32).max(2)
}

/// Convert bits of precision to decimal digits (SymPy's `prec_to_dps`).
pub fn prec_to_dps(prec: u32) -> u32 {
    ((prec as f64 / std::f64::consts::LOG2_10 - 1.0).round().max(1.0)) as u32
}

/// Largest binary exponent a materialized value may carry. mpmath has no
/// bound; this one keeps exponent arithmetic exact in i128 while covering
/// every magnitude whose exponent is itself representable in memory.
pub const MAX_EXPONENT: i128 = 1 << 100;

/// Exponents above this use checked interval arithmetic instead of the
/// plain kernels.
pub const WIDE_EXPONENT: i128 = (i64::MAX / 8) as i128;

fn bitlen(n: &BigInt) -> i128 {
    n.bits() as i128
}

impl BigFloat {
    pub fn zero() -> BigFloat {
        BigFloat { man: BigInt::zero(), exp: 0 }
    }

    pub fn one() -> BigFloat {
        BigFloat { man: BigInt::one(), exp: 0 }
    }

    pub fn is_zero(&self) -> bool {
        self.man.is_zero()
    }

    pub fn is_negative(&self) -> bool {
        self.man.is_negative()
    }

    pub fn sign(&self) -> i32 {
        match self.man.sign() {
            Sign::Minus => -1,
            Sign::NoSign => 0,
            Sign::Plus => 1,
        }
    }

    /// Normalise: strip trailing zero bits of the mantissa.
    fn norm(man: BigInt, exp: i128) -> BigFloat {
        if man.is_zero() {
            return BigFloat::zero();
        }
        let tz = man.trailing_zeros().unwrap_or(0);
        if tz > 0 {
            BigFloat { man: man >> tz, exp: exp + tz as i128 }
        } else {
            BigFloat { man, exp }
        }
    }

    /// Round `man * 2**exp` to `prec` bits (nearest, ties to even).
    pub fn round(man: BigInt, exp: i128, prec: u32) -> BigFloat {
        let bl = bitlen(&man);
        let excess = bl - prec as i128;
        if excess <= 0 {
            return BigFloat::norm(man, exp);
        }
        let neg = man.is_negative();
        let mut m = man.abs();
        let shift = excess as usize;
        let half = BigInt::one() << (shift - 1);
        let low = &m & ((BigInt::one() << shift) - 1u32);
        m >>= shift;
        match low.cmp(&half) {
            Ordering::Greater => m += 1u32,
            Ordering::Equal => {
                if m.is_odd() {
                    m += 1u32
                }
            }
            Ordering::Less => {}
        }
        if neg {
            m = -m;
        }
        BigFloat::norm(m, exp + excess)
    }

    pub fn from_int(n: &BigInt) -> BigFloat {
        BigFloat::norm(n.clone(), 0)
    }

    pub fn from_i64(n: i64) -> BigFloat {
        BigFloat::norm(BigInt::from(n), 0)
    }

    pub fn from_f64(x: f64) -> BigFloat {
        if x == 0.0 || !x.is_finite() {
            return BigFloat::zero();
        }
        let bits = x.to_bits();
        let sign = if bits >> 63 == 0 { 1i64 } else { -1 };
        let e = ((bits >> 52) & 0x7ff) as i64;
        let frac = bits & ((1u64 << 52) - 1);
        let (m, e) = if e == 0 {
            (frac as i64, -1074)
        } else {
            ((frac | (1u64 << 52)) as i64, e - 1075)
        };
        BigFloat::norm(BigInt::from(sign * m), e as i128)
    }

    /// `n / d` rounded to `prec` bits.
    pub fn from_ratio(n: &BigInt, d: &BigInt, prec: u32) -> BigFloat {
        Self::from_int(n).div(&Self::from_int(d), prec).expect("nonzero denominator")
    }

    /// Parse a decimal without an intermediate machine float.
    pub fn from_decimal(s: &str, prec: u32) -> Result<BigFloat, String> {
        decimal::from_decimal(s,prec)
    }

    pub fn to_ratio_exact(&self) -> num_rational::BigRational {
        if self.exp >= 0 {
            num_rational::BigRational::from_integer(&self.man << self.exp as usize)
        } else {
            num_rational::BigRational::new(self.man.clone(), BigInt::one() << (-self.exp) as usize)
        }
    }

    pub fn to_f64(&self) -> f64 {
        if self.man.is_zero() {
            return 0.0;
        }
        let sign = if self.is_negative() { 1u64 << 63 } else { 0 };
        // Subnormals have a fixed quantum of 2**-1074. Round directly to
        // that grid to avoid both powi underflow and double rounding.
        if self.mag() < -1021 {
            let shift = self.exp + 1074;
            let mut mantissa = self.man.abs();
            if shift >= 0 {
                mantissa <<= shift as usize;
            } else {
                let discarded = shift.unsigned_abs();
                if discarded > u128::from(mantissa.bits()) + 1 { return f64::from_bits(sign); }
                let half = BigInt::one() << (discarded as usize - 1);
                let remainder = &mantissa & ((BigInt::one() << discarded as usize) - 1);
                mantissa >>= discarded as usize;
                if remainder > half || (remainder == half && mantissa.is_odd()) { mantissa += 1; }
            }
            return f64::from_bits(sign | mantissa.to_u64().unwrap_or(0));
        }
        let r = BigFloat::round(self.man.clone(), self.exp, 53);
        let exponent = r.mag() - 1;
        if exponent > 1023 { return f64::from_bits(sign | (0x7ff << 52)); }
        let mantissa = r.man.abs().to_u64().unwrap() << (53 - r.man.bits());
        f64::from_bits(sign | (((exponent + 1023) as u64) << 52) | (mantissa & ((1u64 << 52) - 1)))
    }

    pub fn neg(&self) -> BigFloat {
        BigFloat { man: -&self.man, exp: self.exp }
    }

    pub fn abs(&self) -> BigFloat {
        BigFloat { man: self.man.abs(), exp: self.exp }
    }

    /// Magnitude: `floor(log2 |x|) + 1` (0 for zero).
    pub fn mag(&self) -> i128 {
        if self.is_zero() {
            i128::MIN / 4
        } else {
            self.exp.saturating_add(bitlen(&self.man))
        }
    }

    /// Round a stored literal, returning None when the exponent cannot fit.
    pub fn checked_round_to(&self,precision:u32)->Option<BigFloat>{decimal::checked_round(self,precision)}
    pub fn checked_mul(&self,other:&BigFloat,precision:u32)->Option<BigFloat>{decimal::checked_product(self,other,precision)}
    pub fn checked_add(&self,other:&BigFloat,precision:u32)->Option<BigFloat>{decimal::checked_sum(self,other,precision)}

    pub fn add(&self, o: &BigFloat, prec: u32) -> BigFloat {
        if self.is_zero() {
            return BigFloat::round(o.man.clone(), o.exp, prec);
        }
        if o.is_zero() {
            return BigFloat::round(self.man.clone(), self.exp, prec);
        }
        // if one operand is far below the other, keep a sticky contribution
        let (a, b) = if self.mag() >= o.mag() { (self, o) } else { (o, self) };
        let gap = a.mag() - b.mag();
        if gap > prec as i128 + 4 {
            let shift = (prec as i128 + 4 - bitlen(&a.man)).max(0);
            let m = (&a.man << (shift as usize + 2)) + if b.is_negative() { -1 } else { 1 };
            return BigFloat::round(m, a.exp - shift as i128 - 2, prec);
        }
        let e = a.exp.min(b.exp);
        let m = (&a.man << ((a.exp - e) as usize)) + (&b.man << ((b.exp - e) as usize));
        BigFloat::round(m, e, prec)
    }

    pub fn sub(&self, o: &BigFloat, prec: u32) -> BigFloat {
        self.add(&o.neg(), prec)
    }

    pub fn mul(&self, o: &BigFloat, prec: u32) -> BigFloat {
        BigFloat::round(&self.man * &o.man, self.exp + o.exp, prec)
    }

    pub fn mul_i64(&self, k: i64, prec: u32) -> BigFloat {
        BigFloat::round(&self.man * BigInt::from(k), self.exp, prec)
    }

    /// Multiply by `2**k` (exact).
    pub fn ldexp<K: Into<i128>>(&self, k: K) -> BigFloat {
        let k = k.into();
        if self.is_zero() {
            return self.clone();
        }
        BigFloat { man: self.man.clone(), exp: self.exp + k }
    }

    pub fn div(&self, o: &BigFloat, prec: u32) -> Option<BigFloat> {
        if self.exp.unsigned_abs()>WIDE_EXPONENT as u128||o.exp.unsigned_abs()>WIDE_EXPONENT as u128{
            return decimal::checked_quotient(self,o,prec);
        }
        if o.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some(BigFloat::zero());
        }
        let shift = (prec as i128 + 3 + bitlen(&o.man) - bitlen(&self.man)).max(0);
        let num = &self.man << (shift as usize);
        let (q, r) = num.div_rem(&o.man);
        let sticky = if r.is_zero() { BigInt::zero() } else { BigInt::one() };
        let q = if q.is_negative() || (q.is_zero() && (num.is_negative() != o.man.is_negative())) {
            (q << 1usize) - sticky
        } else {
            (q << 1usize) + sticky
        };
        Some(BigFloat::round(q, self.exp - o.exp - shift - 1, prec))
    }

    pub fn div_i64(&self, k: i64, prec: u32) -> BigFloat {
        self.div(&BigFloat::from_i64(k), prec).unwrap()
    }

    pub fn cmp(&self, o: &BigFloat) -> Ordering {
        let sign=self.sign().cmp(&o.sign());
        if sign!=Ordering::Equal { return sign; }
        if self.is_zero() { return Ordering::Equal; }
        let magnitude=(i128::from(self.exp)+i128::from(self.man.bits()))
            .cmp(&(i128::from(o.exp)+i128::from(o.man.bits())));
        if magnitude!=Ordering::Equal {
            return if self.is_negative() { magnitude.reverse() } else { magnitude };
        }
        if self.exp==o.exp { return self.man.cmp(&o.man); }
        // Equal leading bit positions bound alignment by the mantissa sizes,
        // even when the stored exponents differ by billions elsewhere.
        let low=self.exp.min(o.exp);
        let left=&self.man << (i128::from(self.exp)-i128::from(low)) as usize;
        let right=&o.man << (i128::from(o.exp)-i128::from(low)) as usize;
        left.cmp(&right)
    }

    /// Square root (None for negative input).
    pub fn sqrt(&self, prec: u32) -> Option<BigFloat> {
        if self.is_negative() {
            return None;
        }
        if self.is_zero() {
            return Some(BigFloat::zero());
        }
        // make the exponent even and the mantissa large enough
        let mut shift = 2 * (prec as i128 + 4) - bitlen(&self.man);
        if shift < 0 {
            shift = 0;
        }
        if (self.exp - shift) % 2 != 0 {
            shift += 1;
        }
        let m = &self.man << (shift as usize);
        let r = m.sqrt();
        let sticky = if &r * &r == m { BigInt::zero() } else { BigInt::one() };
        Some(BigFloat::round((r << 1usize) + sticky, (self.exp - shift) / 2 - 1, prec))
    }

    /// Truncate towards zero to an integer.
    pub fn to_int_trunc(&self) -> BigInt {
        if self.exp >= 0 {
            &self.man << (self.exp as usize)
        } else {
            let s = (-self.exp) as usize;
            if self.man.is_negative() {
                -((-&self.man) >> s)
            } else {
                &self.man >> s
            }
        }
    }

    /// Round to the nearest integer.
    pub fn to_int_round(&self) -> BigInt {
        if self.exp >= 0 {
            return self.to_int_trunc();
        }
        let half = BigFloat { man: BigInt::one(), exp: -1 };
        let p = (bitlen(&self.man) + 4) as u32;
        let v = if self.is_negative() { self.sub(&half, p) } else { self.add(&half, p) };
        v.to_int_trunc()
    }

    /// Decimal representation with `dps` significant digits, SymPy style
    /// (fixed notation for moderate exponents, `e` notation otherwise).
    /// mpmath's `to_str(x, 0)` as SymPy prints a Float of fewer than 5 bits:
    /// no digit, only the decimal exponent rounded to the nearest power.
    pub fn to_string_no_digits(&self) -> String {
        if self.is_zero() {return "0.e+0".into();}
        let (neg, digits, e10) = self.decimal_digits(3);
        let e10=if matches!(digits.as_bytes().first(),Some(b'5'..=b'9')){e10+1}else{e10};
        format!("{}0.e{}{}",if neg{"-"}else{""},if e10<0{"-"}else{"+"},e10.unsigned_abs())
    }

    pub fn to_string_dps(&self, dps: u32, strip_zeros: bool) -> String {
        let dps=dps.max(1);
        if self.is_zero() {
            return if strip_zeros { "0.0".into() } else { "0.".to_string() + &"0".repeat(dps.saturating_sub(1) as usize) };
        }
        let (neg, digits, e10) = self.decimal_digits(dps);
        let mut digits = digits;
        let mut s = String::new();
        if neg {
            s.push('-');
        }
        // digits d1 d2 ... with value 0.d1d2... * 10**(e10+1) i.e. d1.d2... * 10**e10
        let minimum=(-(dps as i128/3)).min(-5);
        let fixed=minimum<e10 && e10<dps as i128;
        if strip_zeros {
            while digits.len() > 1 && digits.ends_with('0') {
                digits.pop();
            }
        }
        if fixed {
            if e10 >= 0 {
                let int_len = (e10 + 1) as usize;
                if digits.len() <= int_len {
                    s.push_str(&digits);
                    s.push_str(&"0".repeat(int_len - digits.len()));
                    s.push_str(".0");
                    if !strip_zeros {
                        s.pop();
                        s.pop();
                        s.push('.');
                    }
                } else {
                    s.push_str(&digits[..int_len]);
                    s.push('.');
                    s.push_str(&digits[int_len..]);
                }
            } else {
                s.push_str("0.");
                s.push_str(&"0".repeat((-e10 - 1) as usize));
                s.push_str(&digits);
            }
        } else {
            s.push_str(&digits[..1]);
            s.push('.');
            if digits.len() > 1 {
                s.push_str(&digits[1..]);
            } else {
                s.push('0');
            }
            s.push_str(&format!("e{}{}", if e10 >= 0 { "+" } else { "-" }, e10.abs()));
        }
        s
    }

    /// `(negative, dps digits, decimal exponent of the first digit)`.
    fn decimal_digits(&self, dps: u32) -> (bool, String, i128) {
        decimal::digits(self,dps)
    }
}

// ----------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------

/// Chudnovsky binary splitting: returns (P, Q, T).
fn chud_bs(a: u64, b: u64) -> (BigInt, BigInt, BigInt) {
    const C3_24: u64 = 640320u64 * 640320 * 640320 / 24;
    if b - a == 1 {
        if a == 0 {
            return (BigInt::one(), BigInt::one(), BigInt::from(13591409u64));
        }
        let a_ = BigInt::from(a);
        let p: BigInt =
            (BigInt::from(6 * a) - 5) * (BigInt::from(2 * a) - 1) * (BigInt::from(6 * a) - 1);
        let q: BigInt = &a_ * &a_ * &a_ * BigInt::from(C3_24);
        let mut t: BigInt = &p * (BigInt::from(13591409u64) + BigInt::from(545140134u64) * &a_);
        if a % 2 == 1 {
            t = -t;
        }
        return (p, q, t);
    }
    let m = (a + b) / 2;
    let (p1, q1, t1) = chud_bs(a, m);
    let (p2, q2, t2) = chud_bs(m, b);
    (&p1 * &p2, &q1 * &q2, &q2 * t1 + &p1 * t2)
}

/// pi to `prec` bits.
pub fn pi(prec: u32) -> BigFloat {
    let wp = prec + 20;
    let terms = (wp as f64 / 47.11).ceil() as u64 + 2;
    let (_, q, t) = chud_bs(0, terms);
    // pi = 426880 * sqrt(10005) * Q / T
    let s = BigFloat::from_i64(10005).sqrt(wp).unwrap();
    let num = BigFloat::from_int(&q).mul_i64(426880, wp).mul(&s, wp);
    num.div(&BigFloat::from_int(&t), wp).unwrap().round_to(prec)
}

/// ln(2) to `prec` bits, using 2*atanh(1/3).
///
/// The ratio 1/9 needs about three times fewer terms than the former 1/2
/// series. Fixed-point summands are exact floors, so the omitted fractions
/// and geometric tail give an explicit enclosure. Retry if its endpoints
/// do not yet round to the same target value.
pub fn ln2(prec: u32) -> BigFloat {
    let mut wp = prec + 32 + (u32::BITS-prec.max(1).leading_zeros());
    loop {
        let mut power = (BigInt::one() << wp as usize) / 3u32;
        let mut sum = BigInt::zero();
        let mut terms = 0u64;
        while !power.is_zero() {
            sum += &power / (2*terms+1);
            power /= 9u32;
            terms += 1;
        }
        let lower = sum * 2u32;
        // Each included summand loses <1 unit. At the first omitted power,
        // the remaining unscaled geometric tail is <9/8 units.
        let upper = &lower + BigInt::from(2*terms+3);
        let rounded = BigFloat::round(lower,-(wp as i128),prec);
        if rounded == BigFloat::round(upper,-(wp as i128),prec) { return rounded; }
        wp += 32;
    }
}

impl BigFloat {
    pub fn round_to(&self, prec: u32) -> BigFloat {
        BigFloat::round(self.man.clone(), self.exp, prec)
    }

    /// e**x. Panics if the binary exponent exceeds `MAX_EXPONENT`; callers
    /// evaluating untrusted magnitudes should use `checked_exp`.
    pub fn exp(&self, prec: u32) -> BigFloat {
        self.checked_exp(prec).expect("exponential exponent exceeds the BigFloat exponent range")
    }

    pub fn checked_exp(&self, prec: u32) -> Option<BigFloat> {
        if self.is_zero() {
            return Some(BigFloat::one());
        }
        let mag = self.mag();
        if mag > 102 { return None; }
        let wp = prec + 20 + mag.max(0) as u32;
        // x = k*ln2 + r, |r| <= ln2/2
        let l2 = ln2(wp + mag.max(0) as u32);
        let k = self.div(&l2, wp).unwrap().to_int_round();
        let ki = k.to_i128()?;
        let kf = BigFloat::from_int(&k);
        let r = self.sub(&kf.mul(&l2, wp + bitlen(&k) as u32), wp);
        // further reduce: r / 2^s
        let s = ((wp as f64).sqrt() as i64).max(1);
        let r = r.ldexp(-s);
        let wp2 = wp + s as u32;
        // Taylor series in fixed point
        let sum = taylor_exp(&r, wp2);
        let mut y = sum;
        for _ in 0..s {
            y = y.mul(&y, wp2);
        }
        let exponent = y.exp.checked_add(ki)?;
        if exponent.unsigned_abs() > MAX_EXPONENT as u128 { return None; }
        Some(BigFloat { man: y.man, exp: exponent }.round_to(prec))
    }

    /// Natural logarithm (None for x <= 0).
    pub fn ln(&self, prec: u32) -> Option<BigFloat> {
        if self.sign() <= 0 {
            return None;
        }
        if *self == BigFloat::one() {
            return Some(BigFloat::zero());
        }
        // Near one the reduction below cancels as many bits as ln(x) is
        // small; the exactly rounded x - 1 keeps them.
        let d = self.sub(&BigFloat::one(), prec + 30);
        if d.mag() < -4 {
            return Some(d.ln_1p_small(prec));
        }
        self.ln_core(prec)
    }

    /// ln(1 + x) for x > -1, with relative accuracy for tiny x.
    pub fn ln_1p(&self, prec: u32) -> Option<BigFloat> {
        if self.is_zero() {
            return Some(BigFloat::zero());
        }
        if self.mag() < -4 {
            return Some(self.ln_1p_small(prec));
        }
        let guard = (bitlen(&self.man) as u32).saturating_add(prec + 40);
        self.add(&BigFloat::one(), guard).ln(prec)
    }

    fn ln_1p_small(&self, prec: u32) -> BigFloat {
        let small = (-self.mag()) as u32;
        if (small as f64) <= (prec as f64).sqrt() + 8.0 {
            // 1 + x is exact at this width; ln_core loses at most `small` bits.
            let wp = prec + small + 40;
            let x = self.add(&BigFloat::one(), wp + small);
            return x.ln_core(prec + small + 10).unwrap().round_to(prec);
        }
        // ln(1 + x) = 2 atanh(x / (2 + x)), with |x / (2 + x)| < 2**-small.
        let wp = prec + 20;
        let z = self.div(&self.add(&BigFloat::from_i64(2), wp), wp).unwrap();
        let at = if (small as i128) > prec as i128 / 2 + 8 {
            let z3 = z.mul(&z, wp).mul(&z, wp).div_i64(3, wp);
            z.add(&z3, wp)
        } else {
            atanh_series(&z, wp + small + 10)
        };
        at.ldexp(1).round_to(prec)
    }

    fn ln_core(&self, prec: u32) -> Option<BigFloat> {
        let wp = prec + 30;
        // x = m * 2^e with m in [0.5, 1): ln x = ln m + e ln 2
        let e = self.mag();
        let m = self.ldexp(-e);
        // ln m via the atanh series: ln m = 2 atanh((m-1)/(m+1)), after
        // reducing m close to 1 with square roots
        let mut m = m;
        let mut roots = 0i64;
        let target = BigFloat { man: BigInt::one(), exp: -((wp as f64).sqrt() as i128).max(8) };
        loop {
            let d = m.sub(&BigFloat::one(), wp).abs();
            if d.cmp(&target) != Ordering::Greater || roots > 400 {
                break;
            }
            m = m.sqrt(wp).unwrap();
            roots += 1;
        }
        let z = m.sub(&BigFloat::one(), wp).div(&m.add(&BigFloat::one(), wp), wp).unwrap();
        // The fixed-point series has an absolute error; z is below one by
        // about sqrt(wp) bits after the reduction.
        let guard = if z.is_zero() { 0 } else { (-z.mag()).max(0) as u32 };
        let at = atanh_series(&z, wp + guard);
        let lnm = at.ldexp(1 + roots);
        let res = lnm.add(&ln2(wp).mul(&BigFloat::norm(BigInt::from(e), 0), wp), wp);
        Some(res.round_to(prec))
    }

    /// (sin x, cos x)
    pub fn sin_cos(&self, prec: u32) -> (BigFloat, BigFloat) {
        if self.is_zero() {
            return (BigFloat::zero(), BigFloat::one());
        }
        // |sin(x)/x-1| <= x^2/6 and |cos(x)-1| <= x^2/2. In this
        // range these errors are far below the requested relative precision.
        // Avoid an absolute fixed-point grid rounding a tiny argument to zero.
        if self.mag() < -(prec as i128+20) {
            return (self.round_to(prec),BigFloat::one());
        }
        let mag = self.mag().max(0) as u32;
        // Reduce modulo pi/2. Near a multiple of pi/2 the reduced argument
        // cancels; repeat with as many more bits of pi as were lost.
        let mut extra = 0u32;
        let (wp, k, r) = loop {
            let wp = prec + 30 + mag + extra;
            let p = pi(wp + mag);
            let half_pi = p.ldexp(-1);
            let k = self.div(&half_pi, wp).unwrap().to_int_round();
            let r = self.sub(&BigFloat::from_int(&k).mul(&half_pi, wp + bitlen(&k) as u32), wp);
            let lost = if r.is_zero() { None } else { Some((-r.mag()).max(0) as u32) };
            match lost {
                Some(lost) if lost <= extra + 10 => break (wp, k, r),
                Some(lost) if extra < 1 << 22 => extra = lost + 20,
                None if extra < 1 << 22 => extra = extra * 2 + 64,
                _ => break (wp, k, r),
            }
        };
        // r in [-pi/4, pi/4]; halve s times, use cos double angle
        let s = ((wp as f64).sqrt() as i64 / 2).max(1);
        let rs = r.ldexp(-s);
        let relative_guard=if rs.is_zero(){0}else{(-rs.mag()).max(0) as u32};
        let wp2 = wp + 2 * s as u32 + relative_guard;
        let (mut sn, mut cs) = taylor_sin_cos(&rs, wp2);
        for _ in 0..s {
            // sin 2a = 2 sin a cos a, cos 2a = 1 - 2 sin^2 a
            let s2 = sn.mul(&cs, wp2).ldexp(1);
            let c2 = BigFloat::one().sub(&sn.mul(&sn, wp2).ldexp(1), wp2);
            sn = s2;
            cs = c2;
        }
        let q = k.mod_floor(&BigInt::from(4)).to_i64().unwrap();
        let (s_, c_) = match q {
            0 => (sn, cs),
            1 => (cs, sn.neg()),
            2 => (sn.neg(), cs.neg()),
            _ => (cs.neg(), sn),
        };
        (s_.round_to(prec), c_.round_to(prec))
    }

    pub fn sin(&self, prec: u32) -> BigFloat {
        self.sin_cos(prec).0
    }

    pub fn cos(&self, prec: u32) -> BigFloat {
        self.sin_cos(prec).1
    }

    pub fn tan(&self, prec: u32) -> Option<BigFloat> {
        let (s, c) = self.sin_cos(prec + 10);
        s.div(&c, prec)
    }

    /// atan(x)
    pub fn atan(&self, prec: u32) -> BigFloat {
        if self.is_zero() {
            return BigFloat::zero();
        }
        // atan(x) = x - x**3/3 + ...: the next term is below 2**-2prec.
        if self.mag() < -(prec as i128 / 2 + 8) {
            let wp = prec + 10;
            let x3 = self.mul(self, wp).mul(self, wp).div_i64(3, wp);
            return self.sub(&x3, prec);
        }
        let wp = prec + 30 + (-self.mag()).max(0) as u32;
        let neg = self.is_negative();
        let mut x = self.abs();
        let mut add_half_pi = false;
        if x.cmp(&BigFloat::one()) == Ordering::Greater {
            // atan(x) = pi/2 - atan(1/x)
            x = BigFloat::one().div(&x, wp).unwrap();
            add_half_pi = true;
        }
        // reduce: atan(x) = 2 atan(x / (1 + sqrt(1 + x^2)))
        let mut doublings = 0i64;
        let target = BigFloat { man: BigInt::one(), exp: -((wp as f64).sqrt() as i128).max(8) };
        while x.cmp(&target) == Ordering::Greater && doublings < 200 {
            let d = BigFloat::one().add(&x.mul(&x, wp), wp).sqrt(wp).unwrap();
            x = x.div(&BigFloat::one().add(&d, wp), wp).unwrap();
            doublings += 1;
        }
        let guard = if x.is_zero() { 0 } else { (-x.mag()).max(0) as u32 };
        let mut r = atan_series(&x, wp + guard).ldexp(doublings);
        if add_half_pi {
            r = pi(wp).ldexp(-1).sub(&r, wp);
        }
        if neg {
            r = r.neg();
        }
        r.round_to(prec)
    }

    /// atan2(y, x)
    pub fn atan2(y: &BigFloat, x: &BigFloat, prec: u32) -> BigFloat {
        let wp = prec + 10;
        if x.is_zero() {
            if y.is_zero() {
                return BigFloat::zero();
            }
            let hp = pi(prec).ldexp(-1);
            return if y.is_negative() { hp.neg() } else { hp };
        }
        let a = y.div(x, wp).unwrap().atan(wp);
        if !x.is_negative() {
            return a.round_to(prec);
        }
        let p = pi(wp);
        if y.is_negative() {
            a.sub(&p, prec)
        } else {
            a.add(&p, prec)
        }
    }

    /// 1 - x**2 as (1 - x)(1 + x), exact up to the final rounding.
    fn one_minus_square(&self, prec: u32) -> BigFloat {
        let one = BigFloat::one();
        one.sub(self, prec + 4).mul(&one.add(self, prec + 4), prec)
    }

    pub fn asin(&self, prec: u32) -> Option<BigFloat> {
        let wp = prec + 20;
        let c = self.one_minus_square(wp);
        if c.is_negative() {
            return None;
        }
        if c.is_zero() {
            let hp = pi(prec).ldexp(-1);
            return Some(if self.is_negative() { hp.neg() } else { hp });
        }
        Some(BigFloat::atan2(self, &c.sqrt(wp).unwrap(), prec))
    }

    pub fn acos(&self, prec: u32) -> Option<BigFloat> {
        let wp = prec + 20;
        let c = self.one_minus_square(wp);
        if c.is_negative() {
            return None;
        }
        // atan2 keeps the relative accuracy of acos(x) near x = 1, where
        // pi/2 - asin(x) cancels.
        Some(BigFloat::atan2(&c.sqrt(wp).unwrap(), self, prec))
    }

    /// x**y for x > 0 (general real powers).
    pub fn powf(&self, y: &BigFloat, prec: u32) -> Option<BigFloat> {
        if self.is_zero() {
            return if y.sign() > 0 { Some(BigFloat::zero()) } else { None };
        }
        // An exponent above 2**(2**20) only has a representable power for
        // a base within an unrepresentably small distance of one.
        if y.mag() > 1 << 20 { return None; }
        let wp = prec + 20 + y.mag().max(0) as u32;
        let l = self.ln(wp)?;
        l.mul(y, wp).checked_exp(prec)
    }

    /// x**n for an integer n.
    pub fn powi(&self, n: i64, prec: u32) -> Option<BigFloat> {
        decimal::checked_integer_power(self,n,prec)
    }
}

/// Fixed point helper: value of `x` as an integer scaled by 2^wp.
fn to_fixed(x: &BigFloat, wp: u32) -> BigInt {
    let s = x.exp + wp as i128;
    if s >= 0 {
        &x.man << (s as usize)
    } else {
        let sh = (-s) as usize;
        if x.man.is_negative() {
            -((-&x.man) >> sh)
        } else {
            &x.man >> sh
        }
    }
}

fn from_fixed(v: BigInt, wp: u32) -> BigFloat {
    BigFloat::norm(v, -(wp as i128))
}

fn taylor_exp(r: &BigFloat, wp: u32) -> BigFloat {
    let x = to_fixed(r, wp);
    let one = BigInt::one() << (wp as usize);
    let mut sum = one.clone();
    let mut term = one;
    let mut k = 1u64;
    loop {
        term = (&term * &x >> (wp as usize)) / BigInt::from(k);
        if term.is_zero() {
            break;
        }
        sum += &term;
        k += 1;
    }
    from_fixed(sum, wp)
}

fn taylor_sin_cos(r: &BigFloat, wp: u32) -> (BigFloat, BigFloat) {
    let x = to_fixed(r, wp);
    let x2 = (&x * &x) >> (wp as usize);
    let one = BigInt::one() << (wp as usize);
    // sin
    let mut s = x.clone();
    let mut term = x;
    let mut k = 1u64;
    loop {
        term = -((&term * &x2) >> (wp as usize)) / BigInt::from((2 * k) * (2 * k + 1));
        if term.is_zero() {
            break;
        }
        s += &term;
        k += 1;
    }
    let mut c = one.clone();
    let mut term = one;
    let mut k = 1u64;
    loop {
        term = -((&term * &x2) >> (wp as usize)) / BigInt::from((2 * k - 1) * (2 * k));
        if term.is_zero() {
            break;
        }
        c += &term;
        k += 1;
    }
    (from_fixed(s, wp), from_fixed(c, wp))
}

fn atanh_series(z: &BigFloat, wp: u32) -> BigFloat {
    // atanh z = z + z^3/3 + z^5/5 + ...
    let x = to_fixed(z, wp);
    let x2 = (&x * &x) >> (wp as usize);
    let mut sum = x.clone();
    let mut p = x;
    let mut k = 3u64;
    loop {
        p = (&p * &x2) >> (wp as usize);
        let t = &p / BigInt::from(k);
        if t.is_zero() {
            break;
        }
        sum += t;
        k += 2;
    }
    from_fixed(sum, wp)
}

fn atan_series(z: &BigFloat, wp: u32) -> BigFloat {
    let x = to_fixed(z, wp);
    let x2 = (&x * &x) >> (wp as usize);
    let mut sum = x.clone();
    let mut p = x;
    let mut k = 3u64;
    let mut sign = -1i32;
    loop {
        p = (&p * &x2) >> (wp as usize);
        let t = &p / BigInt::from(k);
        if t.is_zero() {
            break;
        }
        if sign < 0 {
            sum -= t;
        } else {
            sum += t;
        }
        sign = -sign;
        k += 2;
    }
    from_fixed(sum, wp)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(x: &BigFloat, dps: u32) -> String {
        x.to_string_dps(dps, false)
    }

    #[test]
    fn comparisons_do_not_allocate_exponent_sized_integers() {
        let tiny=BigFloat{man:BigInt::one(),exp:-1_000_000_000};
        let huge=BigFloat{man:BigInt::one(),exp:1_000_000_000};
        assert_eq!(tiny.cmp(&BigFloat::zero()),Ordering::Greater);
        assert_eq!(tiny.cmp(&BigFloat::one()),Ordering::Less);
        assert_eq!(huge.cmp(&BigFloat::one()),Ordering::Greater);
        let same=BigFloat{man:BigInt::from(2),exp:999_999_999};
        assert_eq!(huge.cmp(&same),Ordering::Equal);
        let negative=BigFloat{man:BigInt::from(-1),exp:-1_000_000_000};
        assert_eq!(negative.cmp(&BigFloat::zero()),Ordering::Less);
        assert_eq!(negative.cmp(&BigFloat::from_i64(-1)),Ordering::Greater);
        let number=crate::number::Number::MpFloat(tiny,100);
        assert_eq!(number.cmp_value(&crate::number::Number::ZERO),Some(Ordering::Greater));
        assert_eq!(number.cmp_value(&crate::number::Number::ONE),Some(Ordering::Less));
    }

    #[test]
    fn constants() {
        let p = pi(dps_to_prec(50));
        assert_eq!(s(&p, 50), "3.1415926535897932384626433832795028841971693993751");
        let e = BigFloat::one().exp(dps_to_prec(40));
        assert_eq!(s(&e, 40), "2.718281828459045235360287471352662497757");
        assert_eq!(s(&ln2(dps_to_prec(30)), 30), "0.693147180559945309417232121458");
    }

    #[test]
    fn functions() {
        let p = dps_to_prec(30);
        let two = BigFloat::from_i64(2);
        assert_eq!(s(&two.sqrt(p).unwrap(), 30), "1.41421356237309504880168872421");
        assert_eq!(s(&BigFloat::from_i64(10).ln(p).unwrap(), 30), "2.30258509299404568401799145468");
        assert_eq!(s(&BigFloat::one().sin(p), 30), "0.841470984807896506652502321630");
        assert_eq!(s(&BigFloat::one().cos(p), 30), "0.540302305868139717400936607443");
        assert_eq!(s(&BigFloat::one().atan(p), 30), "0.785398163397448309615660845820");
        assert_eq!(s(&BigFloat::from_i64(100).sin(p), 30), "-0.506365641109758793656557610460");
        assert!((BigFloat::from_f64(0.1).exp(53).to_f64() - 0.1f64.exp()).abs() < 1e-16);
    }

    #[test]
    fn printing() {
        assert_eq!(BigFloat::from_f64(2.0).to_string_dps(15, false), "2.00000000000000");
        assert_eq!(BigFloat::from_f64(1e-7).to_string_dps(15, false), "1.00000000000000e-7");
        assert_eq!(BigFloat::from_f64(123.5).to_string_dps(5, false), "123.50");
    }

    #[test]
    fn machine_float_conversion_preserves_subnormals_and_rounds_ties_even() {
        let mut patterns = vec![1, 2, (1u64 << 52) - 1, 1u64 << 52, 0x7fef_ffff_ffff_ffff];
        let mut seed = 0x123456789abcdefu64;
        for _ in 0..1000 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            if (seed >> 52) & 0x7ff != 0x7ff { patterns.push(seed); }
        }
        for bits in patterns {
            let value = f64::from_bits(bits);
            assert_eq!(BigFloat::from_f64(value).to_f64().to_bits(), value.to_bits());
        }
        assert_eq!(BigFloat { man: BigInt::from(1), exp: -1075 }.to_f64().to_bits(), 0);
        assert_eq!(BigFloat { man: BigInt::from(3), exp: -1075 }.to_f64().to_bits(), 2);
        assert_eq!(BigFloat { man: BigInt::from(-3), exp: -1075 }.to_f64().to_bits(), (1u64 << 63) | 2);
    }

    #[test]
    fn unrepresentable_exponentials_and_powers_are_checked() {
        let huge = BigFloat::from_decimal("1e100", 200).unwrap();
        assert!(huge.checked_exp(200).is_none());
        assert!(huge.neg().checked_exp(200).is_none());
        assert!(BigFloat::from_i64(10).powf(&huge, 200).is_none());
        // 10**(2**63-1) has a binary exponent beyond i64, as in mpmath.
        let power = BigFloat::from_i64(10).powi(i64::MAX, 200).unwrap();
        assert!(power.mag() > i128::from(i64::MAX));
        let large = BigFloat::from_decimal("1e1000000000000000000", 200).unwrap();
        assert!(large.powi(i64::MAX, 200).is_none());
        assert_eq!(BigFloat::zero().powi(2, 200), Some(BigFloat::zero()));
        let ordinary = BigFloat::from_i64(1000).checked_exp(200).unwrap();
        assert!(ordinary.mag() > 1400);
    }

    #[test]
    fn negligible_addends_do_not_corrupt_short_mantissas() {
        let small = BigFloat::one().ldexp(-1000);
        for precision in [53, 169, 336] {
            assert_eq!(BigFloat::one().add(&small, precision), BigFloat::one());
            assert_eq!(BigFloat::one().sub(&small, precision), BigFloat::one());
            assert_eq!(BigFloat::from_i64(-1).add(&small, precision), BigFloat::from_i64(-1));
            assert_eq!(BigFloat::from_i64(-1).sub(&small, precision), BigFloat::from_i64(-1));
        }
    }
}
