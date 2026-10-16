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

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BigFloat {
    /// Signed mantissa; zero is represented by `man == 0, exp == 0`.
    pub man: BigInt,
    pub exp: i64,
}

/// Convert decimal digits to bits of precision (with no guard bits).
pub fn dps_to_prec(dps: u32) -> u32 {
    (((dps as f64 + 1.0) * std::f64::consts::LOG2_10).round() as u32).max(2)
}

/// Convert bits of precision to decimal digits (SymPy's `prec_to_dps`).
pub fn prec_to_dps(prec: u32) -> u32 {
    ((prec as f64 / std::f64::consts::LOG2_10 - 1.0).round().max(1.0)) as u32
}

fn bitlen(n: &BigInt) -> i64 {
    n.bits() as i64
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
    fn norm(man: BigInt, exp: i64) -> BigFloat {
        if man.is_zero() {
            return BigFloat::zero();
        }
        let tz = man.trailing_zeros().unwrap_or(0);
        if tz > 0 {
            BigFloat { man: man >> tz, exp: exp + tz as i64 }
        } else {
            BigFloat { man, exp }
        }
    }

    /// Round `man * 2**exp` to `prec` bits (nearest, ties to even).
    pub fn round(man: BigInt, exp: i64, prec: u32) -> BigFloat {
        let bl = bitlen(&man);
        let excess = bl - prec as i64;
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
        BigFloat::norm(BigInt::from(sign * m), e)
    }

    /// `n / d` rounded to `prec` bits.
    pub fn from_ratio(n: &BigInt, d: &BigInt, prec: u32) -> BigFloat {
        Self::from_int(n).div(&Self::from_int(d), prec).expect("nonzero denominator")
    }

    /// Parse a decimal without an intermediate machine float.
    pub fn from_decimal(s: &str, prec: u32) -> Result<BigFloat, String> {
        let s = s.trim();
        let (mantissa, exponent) = match s.find(['e', 'E']) {
            Some(i) => (&s[..i], s[i + 1..].parse::<i64>().map_err(|_| "invalid decimal exponent")?),
            None => (s, 0),
        };
        let (negative, mantissa) = match mantissa.as_bytes().first() {
            Some(b'-') => (true, &mantissa[1..]),
            Some(b'+') => (false, &mantissa[1..]),
            _ => (false, mantissa),
        };
        let mut digits = String::new();
        let mut fractional = None;
        for c in mantissa.chars() {
            if c == '.' && fractional.is_none() { fractional = Some(0i64); }
            else if c.is_ascii_digit() {
                digits.push(c);
                if let Some(ref mut count) = fractional { *count += 1; }
            } else { return Err("invalid decimal literal".into()); }
        }
        if digits.is_empty() { return Err("invalid decimal literal".into()); }
        let power = exponent.checked_sub(fractional.unwrap_or(0)).ok_or("decimal exponent overflow")?;
        if power.unsigned_abs() > 1_000_000 { return Err("decimal exponent exceeds supported magnitude (1000000)".into()); }
        let mut n: BigInt = digits.parse().map_err(|_| "invalid decimal literal")?;
        if negative { n = -n; }
        let ten = num_traits::pow::Pow::pow(BigInt::from(10), power.unsigned_abs());
        Ok(if power >= 0 { BigFloat::from_int(&(n * ten)).round_to(prec) }
            else { BigFloat::from_ratio(&n, &ten, prec) })
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
                if discarded > mantissa.bits() + 1 { return f64::from_bits(sign); }
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
    pub fn mag(&self) -> i64 {
        if self.is_zero() {
            i64::MIN / 4
        } else {
            self.exp + bitlen(&self.man)
        }
    }

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
        if gap > prec as i64 + 4 {
            let shift = (prec as i64 + 4 - bitlen(&a.man)).max(0);
            let m = (&a.man << (shift as usize + 2)) + if b.is_negative() { -1 } else { 1 };
            return BigFloat::round(m, a.exp - shift as i64 - 2, prec);
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
    pub fn ldexp(&self, k: i64) -> BigFloat {
        if self.is_zero() {
            return self.clone();
        }
        BigFloat { man: self.man.clone(), exp: self.exp + k }
    }

    pub fn div(&self, o: &BigFloat, prec: u32) -> Option<BigFloat> {
        if o.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some(BigFloat::zero());
        }
        let shift = (prec as i64 + 3 + bitlen(&o.man) - bitlen(&self.man)).max(0);
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
        let d = self.sub(o, (bitlen(&self.man).max(bitlen(&o.man)) + 8 + (self.exp - o.exp).unsigned_abs() as i64) as u32);
        d.sign().cmp(&0)
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
        let mut shift = 2 * (prec as i64 + 4) - bitlen(&self.man);
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
    pub fn to_string_dps(&self, dps: u32, strip_zeros: bool) -> String {
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
        let fixed = (-5..15).contains(&e10) || (dps > 15 && e10 >= -5 && e10 < dps as i64);
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
    fn decimal_digits(&self, dps: u32) -> (bool, String, i64) {
        let neg = self.is_negative();
        let a = self.abs();
        // estimate decimal exponent
        let mut e10 = ((a.mag() - 1) as f64 * std::f64::consts::LOG10_2).floor() as i64;
        loop {
            // scaled = a * 10**(dps - 1 - e10), rounded to integer
            let k = dps as i64 - 1 - e10;
            let prec = (a.mag().unsigned_abs() as f64 + (k.unsigned_abs() as f64) * 3.33 + 64.0) as u32 + dps * 4;
            let scaled = if k >= 0 {
                a.mul(&BigFloat::from_int(&num_traits::pow::Pow::pow(BigInt::from(10), k as u64)), prec)
            } else {
                a.div(&BigFloat::from_int(&num_traits::pow::Pow::pow(BigInt::from(10), (-k) as u64)), prec)
                    .unwrap()
            };
            let n = scaled.to_int_round();
            let s = n.to_string();
            if s.len() as u32 > dps {
                e10 += 1;
                continue;
            }
            if (s.len() as u32) < dps {
                e10 -= 1;
                continue;
            }
            return (neg, s, e10);
        }
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

/// ln(2) to `prec` bits: ln 2 = sum 1/(k 2^k).
pub fn ln2(prec: u32) -> BigFloat {
    let wp = prec + 20;
    // fixed point: ln2 = sum_{k>=1} 2^(wp-k)/k
    let one = BigInt::one() << (wp as usize);
    let mut s = BigInt::zero();
    let mut k = 1u64;
    loop {
        let term = (&one >> (k as usize)) / BigInt::from(k);
        if term.is_zero() {
            break;
        }
        s += term;
        k += 1;
    }
    BigFloat::round(s, -(wp as i64), prec)
}

impl BigFloat {
    pub fn round_to(&self, prec: u32) -> BigFloat {
        BigFloat::round(self.man.clone(), self.exp, prec)
    }

    /// e**x. Panics if the binary exponent cannot fit in i64; callers
    /// evaluating untrusted magnitudes should use `checked_exp`.
    pub fn exp(&self, prec: u32) -> BigFloat {
        self.checked_exp(prec).expect("exponential exponent exceeds the BigFloat i64 range")
    }

    pub fn checked_exp(&self, prec: u32) -> Option<BigFloat> {
        if self.is_zero() {
            return Some(BigFloat::one());
        }
        let mag = self.mag();
        if mag > 64 { return None; }
        let wp = prec + 20 + mag.max(0) as u32;
        // x = k*ln2 + r, |r| <= ln2/2
        let l2 = ln2(wp + mag.max(0) as u32);
        let k = self.div(&l2, wp).unwrap().to_int_round();
        let ki = k.to_i64()?;
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
        let wp = prec + 30;
        // x = m * 2^e with m in [0.5, 1): ln x = ln m + e ln 2
        let e = self.mag();
        let m = self.ldexp(-e);
        // ln m via the atanh series: ln m = 2 atanh((m-1)/(m+1)), after
        // reducing m close to 1 with square roots
        let mut m = m;
        let mut roots = 0i64;
        let target = BigFloat { man: BigInt::one(), exp: -((wp as f64).sqrt() as i64).max(8) };
        loop {
            let d = m.sub(&BigFloat::one(), wp).abs();
            if d.cmp(&target) != Ordering::Greater || roots > 400 {
                break;
            }
            m = m.sqrt(wp).unwrap();
            roots += 1;
        }
        let z = m.sub(&BigFloat::one(), wp).div(&m.add(&BigFloat::one(), wp), wp).unwrap();
        let at = atanh_series(&z, wp);
        let lnm = at.ldexp(1 + roots);
        let res = lnm.add(&ln2(wp).mul_i64(e, wp), wp);
        Some(res.round_to(prec))
    }

    /// (sin x, cos x)
    pub fn sin_cos(&self, prec: u32) -> (BigFloat, BigFloat) {
        if self.is_zero() {
            return (BigFloat::zero(), BigFloat::one());
        }
        let mag = self.mag().max(0) as u32;
        let wp = prec + 30 + mag;
        // reduce modulo pi/2
        let p = pi(wp + mag);
        let half_pi = p.ldexp(-1);
        let k = self.div(&half_pi, wp).unwrap().to_int_round();
        let r = self.sub(&BigFloat::from_int(&k).mul(&half_pi, wp + bitlen(&k) as u32), wp);
        // r in [-pi/4, pi/4]; halve s times, use cos double angle
        let s = ((wp as f64).sqrt() as i64 / 2).max(1);
        let rs = r.ldexp(-s);
        let wp2 = wp + 2 * s as u32;
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
        let wp = prec + 30;
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
        let target = BigFloat { man: BigInt::one(), exp: -((wp as f64).sqrt() as i64).max(8) };
        while x.cmp(&target) == Ordering::Greater && doublings < 200 {
            let d = BigFloat::one().add(&x.mul(&x, wp), wp).sqrt(wp).unwrap();
            x = x.div(&BigFloat::one().add(&d, wp), wp).unwrap();
            doublings += 1;
        }
        let mut r = atan_series(&x, wp).ldexp(doublings);
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

    pub fn asin(&self, prec: u32) -> Option<BigFloat> {
        let wp = prec + 20;
        let one = BigFloat::one();
        let c = one.sub(&self.mul(self, wp), wp);
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
        let a = self.asin(wp)?;
        Some(pi(wp).ldexp(-1).sub(&a, prec))
    }

    /// x**y for x > 0 (general real powers).
    pub fn powf(&self, y: &BigFloat, prec: u32) -> Option<BigFloat> {
        if self.is_zero() {
            return if y.sign() > 0 { Some(BigFloat::zero()) } else { None };
        }
        let wp = prec + 20 + y.mag().max(0) as u32;
        let l = self.ln(wp)?;
        l.mul(y, wp).checked_exp(prec)
    }

    /// x**n for an integer n.
    pub fn powi(&self, n: i64, prec: u32) -> Option<BigFloat> {
        let wp = prec + 10 + 64 - (n.unsigned_abs().leading_zeros()) ;
        let mut base = self.clone();
        let mut e = n.unsigned_abs();
        let mut acc = BigFloat::one();
        while e > 0 {
            if e & 1 == 1 {
                acc.exp.checked_add(base.exp)?;
                acc = acc.mul(&base, wp);
                if !acc.is_zero() && acc.mag().unsigned_abs() > 3_321_930 { return None; }
            }
            e >>= 1;
            if e > 0 {
                base.exp.checked_add(base.exp)?;
                base = base.mul(&base, wp);
                if !base.is_zero() && base.mag().unsigned_abs() > 3_321_930 { return None; }
            }
        }
        if n < 0 {
            BigFloat::one().div(&acc, prec)
        } else {
            Some(acc.round_to(prec))
        }
    }
}

/// Fixed point helper: value of `x` as an integer scaled by 2^wp.
fn to_fixed(x: &BigFloat, wp: u32) -> BigInt {
    let s = x.exp + wp as i64;
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
    BigFloat::norm(v, -(wp as i64))
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
        assert!(BigFloat::from_i64(10).powi(i64::MAX, 200).is_none());
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
