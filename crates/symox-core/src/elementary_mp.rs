//! Multiprecision elementary functions.
//!
//! Each component of a value carries the requested relative accuracy for
//! an exactly known argument: small values come from cancellation-free
//! formulas (log1p forms, products such as (1 - x)(1 + x), the Hull,
//! Fairgrieve and Tang algorithm for the complex inverse sine and cosine).
//! Branch cuts and the values on them follow mpmath, through which SymPy
//! evaluates these functions. `evaluate` then supplies the argument with
//! as many bits as the condition of the function at that point demands.

use std::cmp::Ordering;

use crate::bigfloat::{pi, BigFloat};
use crate::evalf_mp::MpComplex;
use crate::expr::Expr;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Elementary {
    Exp,
    Log,
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Sech,
    Csch,
    Asin,
    Acos,
    Atan,
    Acot,
    Asec,
    Acsc,
    Asinh,
    Acosh,
    Atanh,
    Acoth,
    Asech,
    Acsch,
}

use Elementary as F;

impl Elementary {
    pub(crate) fn of_func(f: &crate::expr::Func) -> Option<Elementary> {
        use crate::expr::Func;
        Some(match f {
            Func::Exp => F::Exp,
            Func::Log => F::Log,
            Func::Sin => F::Sin,
            Func::Cos => F::Cos,
            Func::Tan => F::Tan,
            Func::Cot => F::Cot,
            Func::Sec => F::Sec,
            Func::Csc => F::Csc,
            Func::Sinh => F::Sinh,
            Func::Cosh => F::Cosh,
            Func::Tanh => F::Tanh,
            Func::Coth => F::Coth,
            Func::Asin => F::Asin,
            Func::Acos => F::Acos,
            Func::Atan => F::Atan,
            Func::Acot => F::Acot,
            Func::Asinh => F::Asinh,
            Func::Acosh => F::Acosh,
            Func::Atanh => F::Atanh,
            _ => return None,
        })
    }

    pub(crate) fn of_special(s: crate::special::SpecialFn) -> Option<Elementary> {
        use crate::special::SpecialFn as S;
        Some(match s {
            S::Sech => F::Sech,
            S::Csch => F::Csch,
            S::Asec => F::Asec,
            S::Acsc => F::Acsc,
            S::Asech => F::Asech,
            S::Acsch => F::Acsch,
            S::Acoth => F::Acoth,
            _ => return None,
        })
    }
}

/// Bits of argument and working precision beyond which periodic reduction
/// is not attempted (the value then stays symbolic).
const REDUCTION_LIMIT: i128 = 1 << 30;
/// Working precision bound for one elementary value.
const PRECISION_LIMIT: i128 = 1 << 24;

fn one() -> BigFloat {
    BigFloat::one()
}

fn zero() -> BigFloat {
    BigFloat::zero()
}

fn half_pi(p: u32) -> BigFloat {
    pi(p + 2).ldexp(-1).round_to(p)
}

fn cx(re: BigFloat, im: BigFloat) -> MpComplex {
    MpComplex { re, im }
}

fn real(re: BigFloat) -> MpComplex {
    MpComplex::real(re)
}

fn signed(x: BigFloat, negative: bool) -> BigFloat {
    if negative { x.neg() } else { x }
}

fn tiny(x: &BigFloat, p: u32) -> bool {
    !x.is_zero() && x.mag() < -(p as i128 / 2 + 8)
}

/// Bits by which |x| is below one.
fn smallness(x: &BigFloat) -> u32 {
    if x.is_zero() { 0 } else { (-x.mag()).clamp(0, PRECISION_LIMIT) as u32 }
}

/// x + (a/b) x**3 for |x| below 2**(-p/2-8), where the next term is below
/// 2**(-2p) relative to x.
fn cubic(x: &BigFloat, a: i64, b: i64, p: u32) -> BigFloat {
    let wp = p + 10;
    let x3 = x.mul(x, wp).mul(x, wp).mul_i64(a, wp).div_i64(b, wp);
    x.add(&x3, p)
}

fn inv(x: &BigFloat, p: u32) -> Option<BigFloat> {
    one().div(x, p)
}

fn exceeds(x: &BigFloat, limit: i128) -> bool {
    !x.is_zero() && x.mag() > limit
}

// ---------------------------------------------------------------- real

fn sinh_cosh(x: &BigFloat, p: u32) -> Option<(BigFloat, BigFloat)> {
    if x.is_zero() {
        return Some((zero(), one()));
    }
    if tiny(x, p) {
        return Some((cubic(x, 1, 6, p), one()));
    }
    let wp = p + 20 + smallness(x);
    let e = x.checked_exp(wp)?;
    let r = inv(&e, wp)?;
    Some((e.sub(&r, wp).ldexp(-1).round_to(p), e.add(&r, wp).ldexp(-1).round_to(p)))
}

fn tanh_real(x: &BigFloat, p: u32) -> Option<BigFloat> {
    if x.is_zero() {
        return Some(zero());
    }
    if tiny(x, p) {
        return Some(cubic(x, -1, 3, p));
    }
    // 1 - |tanh x| < 2 exp(-2|x|) is below half an ulp of one.
    if x.mag() > (64 - p.leading_zeros()) as i128 + 2 {
        return Some(signed(one(), x.is_negative()));
    }
    let wp = p + 20 + smallness(x);
    let e = x.ldexp(1).checked_exp(wp)?;
    e.sub(&one(), wp).div(&e.add(&one(), wp), p)
}

fn asinh_real(x: &BigFloat, p: u32) -> Option<BigFloat> {
    if x.is_zero() {
        return Some(zero());
    }
    if x.is_negative() {
        return asinh_real(&x.neg(), p).map(|v| v.neg());
    }
    if tiny(x, p) {
        return Some(cubic(x, -1, 6, p));
    }
    // ln(1 + u), u = x + x**2/(1 + sqrt(1 + x**2)).
    let wp = p + 20;
    let x2 = x.mul(x, wp);
    let s = x2.add(&one(), wp).sqrt(wp)?;
    let u = x.add(&x2.div(&s.add(&one(), wp), wp)?, wp);
    u.ln_1p(p)
}

/// acosh(w) from w - 1 >= 0 and w + 1, both known with relative accuracy.
fn acosh_parts(dm1: &BigFloat, dp1: &BigFloat, p: u32) -> Option<BigFloat> {
    if dm1.is_zero() {
        return Some(zero());
    }
    let wp = p + 20;
    let u = dm1.add(&dm1.mul(dp1, wp).sqrt(wp)?, wp);
    u.ln_1p(p)
}

fn acosh_ge1(x: &BigFloat, p: u32) -> Option<BigFloat> {
    let wp = p + 20;
    acosh_parts(&x.sub(&one(), wp), &x.add(&one(), wp), p)
}

/// 0.5 ln((1 + t)/(1 - t)) written as 0.5 ln1p(2t/(1 - t)) from t and
/// 1 - t, for |t| < 1.
fn atanh_parts(t: &BigFloat, one_minus_t: &BigFloat, p: u32) -> Option<BigFloat> {
    if tiny(t, p) {
        return Some(cubic(t, 1, 3, p));
    }
    let wp = p + 20 + smallness(t);
    Some(t.ldexp(1).div(one_minus_t, wp)?.ln_1p(p)?.ldexp(-1))
}

/// The value of an elementary function at a real argument; out of the real
/// domain the complex value on the branch cut as mpmath returns it.
fn real_value(f: Elementary, x: &BigFloat, p: u32) -> Option<MpComplex> {
    let wp = p + 20;
    let unit = x.abs().cmp(&one());
    let negative = x.is_negative();
    Some(match f {
        F::Exp => real(x.checked_exp(p)?),
        F::Log if x.is_zero() => return None,
        F::Log if negative => cx(x.neg().ln(p)?, pi(p)),
        F::Log => real(x.ln(p)?),
        F::Sin | F::Cos | F::Tan | F::Cot | F::Sec | F::Csc => {
            if exceeds(x, REDUCTION_LIMIT) {
                return None;
            }
            let (s, c) = x.sin_cos(wp);
            real(match f {
                F::Sin => s.round_to(p),
                F::Cos => c.round_to(p),
                F::Tan => s.div(&c, p)?,
                F::Cot => c.div(&s, p)?,
                F::Sec => inv(&c, p)?,
                _ => inv(&s, p)?,
            })
        }
        F::Sinh | F::Cosh | F::Sech | F::Csch => {
            let (s, c) = sinh_cosh(x, wp)?;
            real(match f {
                F::Sinh => s.round_to(p),
                F::Cosh => c.round_to(p),
                F::Sech => inv(&c, p)?,
                _ => inv(&s, p)?,
            })
        }
        F::Tanh => real(tanh_real(x, p)?),
        F::Coth => real(inv(&tanh_real(x, wp)?, p)?),
        F::Asin | F::Acos if unit == Ordering::Greater => {
            let a = acosh_ge1(&x.abs(), p)?;
            match (f, negative) {
                (F::Asin, false) => cx(half_pi(p), a.neg()),
                (F::Asin, true) => cx(half_pi(p).neg(), a),
                (_, false) => cx(zero(), a),
                (_, true) => cx(pi(p), a.neg()),
            }
        }
        F::Asin => real(x.asin(p)?),
        F::Acos => real(x.acos(p)?),
        F::Atan => real(x.atan(p)),
        F::Acot if x.is_zero() => real(half_pi(p)),
        F::Acot => real(inv(x, wp)?.atan(p)),
        F::Asinh => real(asinh_real(x, p)?),
        F::Acosh => match x.cmp(&one()) {
            Ordering::Less if unit == Ordering::Less || unit == Ordering::Equal => cx(zero(), x.acos(p)?),
            Ordering::Less => cx(acosh_ge1(&x.neg(), p)?, pi(p)),
            _ => real(acosh_ge1(x, p)?),
        },
        F::Atanh | F::Acoth if unit == Ordering::Equal || x.is_zero() && f == F::Acoth => return None,
        F::Atanh if unit == Ordering::Less => real(atanh_parts(x, &one().sub(x, wp), p)?),
        F::Acoth if unit == Ordering::Greater => {
            // atanh(1/x) = 0.5 ln1p(2/(|x| - 1)) with the sign of x.
            let a = x.abs();
            let v = BigFloat::from_i64(2).div(&a.sub(&one(), wp), wp)?.ln_1p(p)?.ldexp(-1);
            real(signed(v, negative))
        }
        F::Atanh | F::Acoth => {
            // Beyond the unit interval: atanh(t), |t| > 1, has imaginary
            // part -pi/2 for t > 1 and pi/2 for t < -1 (mpmath).
            let a = x.abs();
            let re = if f == F::Atanh {
                BigFloat::from_i64(2).div(&a.sub(&one(), wp), wp)?.ln_1p(p)?.ldexp(-1)
            } else {
                // t = 1/x: 2/(|t| - 1) = 2|x|/(1 - |x|).
                a.ldexp(1).div(&one().sub(&a, wp), wp)?.ln_1p(p)?.ldexp(-1)
            };
            cx(signed(re, negative), signed(half_pi(p), !negative))
        }
        F::Asec | F::Acsc | F::Asech | F::Acsch if x.is_zero() => return None,
        F::Asec | F::Acsc if unit != Ordering::Less => {
            // With s = sqrt(x**2 - 1): asec x = atan(s) or pi - atan(s),
            // acsc x = sign(x) atan(1/s).
            let a = x.abs();
            let s = a.sub(&one(), wp).mul(&a.add(&one(), wp), wp).sqrt(wp)?;
            if f == F::Asec {
                let t = s.atan(wp);
                real(if negative { pi(wp).sub(&t, p) } else { t.round_to(p) })
            } else if s.is_zero() {
                real(signed(half_pi(p), negative))
            } else {
                real(signed(inv(&s, wp)?.atan(p), negative))
            }
        }
        F::Asec | F::Acsc => {
            // acos(w), asin(w) with |w| = 1/|x| > 1.
            let a = x.abs();
            let ach = acosh_parts(&one().sub(&a, wp).div(&a, wp)?, &one().add(&a, wp).div(&a, wp)?, p)?;
            match (f, negative) {
                (F::Asec, false) => cx(zero(), ach),
                (F::Asec, true) => cx(pi(p), ach.neg()),
                (_, false) => cx(half_pi(p), ach.neg()),
                (_, true) => cx(half_pi(p).neg(), ach),
            }
        }
        F::Asech => {
            // acosh(w), w = 1/x.
            let a = x.abs();
            if unit != Ordering::Greater {
                let ach = acosh_parts(&one().sub(&a, wp).div(&a, wp)?, &one().add(&a, wp).div(&a, wp)?, p)?;
                if negative { cx(ach, pi(p)) } else { real(ach) }
            } else {
                // |w| < 1: acosh(w) = i acos(w), acos(w) = atan2(sqrt(x**2 - 1), sign x).
                let s = a.sub(&one(), wp).mul(&a.add(&one(), wp), wp).sqrt(wp)?;
                let t = s.atan(wp);
                cx(zero(), if negative { pi(wp).sub(&t, p) } else { t.round_to(p) })
            }
        }
        F::Acsch => real(asinh_real(&inv(x, wp)?, p)?),
    })
}

// ------------------------------------------------------------- complex

fn complex_inv(z: &MpComplex, p: u32) -> Option<MpComplex> {
    z.inv(p)
}

fn square(x: &BigFloat, p: u32) -> BigFloat {
    x.mul(x, p)
}

pub(crate) fn log(z: &MpComplex, p: u32) -> Option<MpComplex> {
    if z.is_real() {
        return real_value(F::Log, &z.re, p);
    }
    let wp = p + 20;
    let size = z.re.mag().max(z.im.mag());
    let re = if (-1..=2).contains(&size) {
        // 0.5 ln1p(|z|**2 - 1), with (u - 1)(u + 1) + v**2 for the larger
        // component u: exact near |z| = 1.
        let (u, v) = if z.re.abs().cmp(&z.im.abs()) == Ordering::Less { (&z.im, &z.re) } else { (&z.re, &z.im) };
        let s = u.sub(&one(), wp).mul(&u.add(&one(), wp), wp).add(&square(v, wp), wp);
        if s.cmp(&one().neg()) != Ordering::Greater {
            return None;
        }
        s.ln_1p(p)?.ldexp(-1)
    } else {
        square(&z.re, wp).add(&square(&z.im, wp), wp).ln(p)?.ldexp(-1)
    };
    Some(cx(re, BigFloat::atan2(&z.im, &z.re, p)))
}

/// asin(z) and acos(z) for a non-real z (Hull, Fairgrieve and Tang).
fn asin_acos(z: &MpComplex, p: u32) -> Option<(MpComplex, MpComplex)> {
    let wp = p + 30;
    let x = z.re.abs();
    let y = z.im.abs();
    let half = |v: BigFloat| v.ldexp(-1);
    let xp1 = x.add(&one(), wp);
    let xm1 = x.sub(&one(), wp);
    let y2 = square(&y, wp);
    let r = square(&xp1, wp).add(&y2, wp).sqrt(wp)?;
    let s = square(&xm1, wp).add(&y2, wp).sqrt(wp)?;
    let a = half(r.add(&s, wp));
    let b = x.div(&a, wp)?;
    let apx = a.add(&x, wp);
    let r_xp1 = r.add(&xp1, wp);
    let below = x.cmp(&one()) != Ordering::Greater;
    let crossover = BigFloat::from_ratio(&6417.into(), &10000.into(), 64);
    let (asin_re, acos_re) = if b.cmp(&crossover) != Ordering::Greater {
        (b.asin(wp)?, b.acos(wp)?)
    } else {
        let d = if below {
            // s - (x - 1) = s + (1 - x)
            half(apx.mul(&y2.div(&r_xp1, wp)?.add(&s.sub(&xm1, wp), wp), wp)).sqrt(wp)?
        } else {
            let t = apx.div(&r_xp1, wp)?.add(&apx.div(&s.add(&xm1, wp), wp)?, wp);
            y.mul(&half(t).sqrt(wp)?, wp)
        };
        (x.div(&d, wp)?.atan(wp), d.div(&x, wp)?.atan(wp))
    };
    let im = if a.cmp(&BigFloat::from_ratio(&3.into(), &2.into(), 64)) != Ordering::Greater {
        let am1 = if x.cmp(&one()) == Ordering::Less {
            half(y2.div(&r_xp1, wp)?.add(&y2.div(&s.sub(&xm1, wp), wp)?, wp))
        } else {
            half(y2.div(&r_xp1, wp)?.add(&s.add(&xm1, wp), wp))
        };
        am1.add(&am1.mul(&a.add(&one(), wp), wp).sqrt(wp)?, wp).ln_1p(p)?
    } else {
        a.add(&a.sub(&one(), wp).mul(&a.add(&one(), wp), wp).sqrt(wp)?, wp).ln(p)?
    };
    let (xneg, yneg) = (z.re.is_negative(), z.im.is_negative());
    let asin = cx(signed(asin_re, xneg).round_to(p), signed(im.clone(), yneg));
    let acos_re = if xneg { pi(wp).sub(&acos_re, p) } else { acos_re.round_to(p) };
    let acos = cx(acos_re, signed(im, !yneg));
    Some((asin, acos))
}

fn asin_any(z: &MpComplex, p: u32) -> Option<MpComplex> {
    if z.is_real() { real_value(F::Asin, &z.re, p) } else { Some(asin_acos(z, p)?.0) }
}

fn acos_any(z: &MpComplex, p: u32) -> Option<MpComplex> {
    if z.is_real() { real_value(F::Acos, &z.re, p) } else { Some(asin_acos(z, p)?.1) }
}

fn acosh_any(z: &MpComplex, p: u32) -> Option<MpComplex> {
    if z.is_real() {
        return real_value(F::Acosh, &z.re, p);
    }
    // acosh z = i acos z for Im z > 0 and -i acos z for Im z < 0.
    let c = asin_acos(z, p)?.1;
    Some(if z.im.is_negative() { cx(c.im, c.re.neg()) } else { cx(c.im.neg(), c.re) })
}

fn atan_any(z: &MpComplex, p: u32) -> Option<MpComplex> {
    if z.is_real() {
        return real_value(F::Atan, &z.re, p);
    }
    let wp = p + 20;
    let (x, y) = (&z.re, &z.im);
    let one_minus_y = one().sub(y, wp);
    let x2 = square(x, wp);
    let q = x2.add(&square(&one_minus_y, wp), wp);
    if q.is_zero() {
        return None;
    }
    // Im atan z = ln(|1 + i z|**2 / |1 - i z|**2)/4, as ln1p near one.
    let u = y.ldexp(2).div(&q, wp)?;
    let im = if u.mag() < 0 {
        u.ln_1p(p)?.ldexp(-2)
    } else {
        x2.add(&square(&one().add(y, wp), wp), wp).div(&q, wp)?.ln(p)?.ldexp(-2)
    };
    let re = if x.is_zero() {
        match y.abs().cmp(&one()) {
            Ordering::Less => zero(),
            Ordering::Equal => return None,
            // On the cut: the side of the logarithmic formula.
            Ordering::Greater => signed(half_pi(p), y.is_negative()),
        }
    } else {
        let d = one_minus_y.mul(&one().add(y, wp), wp).sub(&x2, wp);
        BigFloat::atan2(&x.ldexp(1), &d, wp).ldexp(-1).round_to(p)
    };
    Some(cx(re, im))
}

/// i z
fn times_i(z: &MpComplex) -> MpComplex {
    cx(z.im.neg(), z.re.clone())
}

/// -i z
fn times_minus_i(z: &MpComplex) -> MpComplex {
    cx(z.im.clone(), z.re.neg())
}

fn complex_value(f: Elementary, z: &MpComplex, p: u32) -> Option<MpComplex> {
    let wp = p + 20;
    let (x, y) = (&z.re, &z.im);
    Some(match f {
        F::Exp => {
            if exceeds(y, REDUCTION_LIMIT) {
                return None;
            }
            z.exp(p)?
        }
        F::Log => log(z, p)?,
        F::Sin | F::Cos | F::Tan | F::Cot | F::Sec | F::Csc => {
            if exceeds(x, REDUCTION_LIMIT) {
                return None;
            }
            let (sx, cx_) = x.sin_cos(wp);
            let (sy, cy) = sinh_cosh(y, wp)?;
            let m = |a: &BigFloat, b: &BigFloat| a.mul(b, wp);
            match f {
                F::Sin => cx(m(&sx, &cy).round_to(p), m(&cx_, &sy).round_to(p)),
                F::Cos => cx(m(&cx_, &cy).round_to(p), m(&sx, &sy).neg().round_to(p)),
                F::Sec => complex_inv(&cx(m(&cx_, &cy), m(&sx, &sy).neg()), p)?,
                F::Csc => complex_inv(&cx(m(&sx, &cy), m(&cx_, &sy)), p)?,
                F::Tan => {
                    let d = square(&cx_, wp).add(&square(&sy, wp), wp);
                    cx(m(&sx, &cx_).div(&d, p)?, m(&sy, &cy).div(&d, p)?)
                }
                _ => {
                    let d = square(&sx, wp).add(&square(&sy, wp), wp);
                    cx(m(&sx, &cx_).div(&d, p)?, m(&sy, &cy).neg().div(&d, p)?)
                }
            }
        }
        F::Sinh | F::Cosh | F::Tanh | F::Coth | F::Sech | F::Csch => {
            if exceeds(y, REDUCTION_LIMIT) {
                return None;
            }
            let (sx, cx_) = sinh_cosh(x, wp)?;
            let (sy, cy) = y.sin_cos(wp);
            let m = |a: &BigFloat, b: &BigFloat| a.mul(b, wp);
            match f {
                F::Sinh => cx(m(&sx, &cy).round_to(p), m(&cx_, &sy).round_to(p)),
                F::Cosh => cx(m(&cx_, &cy).round_to(p), m(&sx, &sy).round_to(p)),
                F::Csch => complex_inv(&cx(m(&sx, &cy), m(&cx_, &sy)), p)?,
                F::Sech => complex_inv(&cx(m(&cx_, &cy), m(&sx, &sy)), p)?,
                F::Tanh => {
                    let d = square(&sx, wp).add(&square(&cy, wp), wp);
                    cx(m(&sx, &cx_).div(&d, p)?, m(&sy, &cy).div(&d, p)?)
                }
                _ => {
                    let d = square(&sx, wp).add(&square(&sy, wp), wp);
                    cx(m(&sx, &cx_).div(&d, p)?, m(&sy, &cy).neg().div(&d, p)?)
                }
            }
        }
        F::Asin => asin_acos(z, p)?.0,
        F::Acos => asin_acos(z, p)?.1,
        F::Atan => atan_any(z, p)?,
        F::Acot => atan_any(&complex_inv(z, wp)?, p)?,
        // asinh z = -i asin(i z), atanh z = -i atan(i z)
        F::Asinh => times_minus_i(&asin_any(&times_i(z), p)?),
        F::Atanh => times_minus_i(&atan_any(&times_i(z), p)?),
        F::Acosh => acosh_any(z, p)?,
        F::Asec => acos_any(&complex_inv(z, wp)?, p)?,
        F::Acsc => asin_any(&complex_inv(z, wp)?, p)?,
        F::Asech => acosh_any(&complex_inv(z, wp)?, p)?,
        F::Acsch => {
            let w = complex_inv(z, wp)?;
            times_minus_i(&asin_any(&times_i(&w), p)?)
        }
        F::Acoth => {
            let w = complex_inv(z, wp)?;
            times_minus_i(&atan_any(&times_i(&w), p)?)
        }
    })
}

/// The value of `f` at an exactly known `z`, each component accurate to
/// about `p` bits.
pub(crate) fn value(f: Elementary, z: &MpComplex, p: u32) -> Option<MpComplex> {
    if z.is_real() { real_value(f, &z.re, p) } else { complex_value(f, z, p) }
}

// ----------------------------------------------------------- condition

fn magnitude(z: &MpComplex) -> Option<i128> {
    [&z.re, &z.im].into_iter().filter(|v| !v.is_zero()).map(BigFloat::mag).max()
}

fn product(a: MpComplex, b: MpComplex) -> Option<i128> {
    magnitude(&a.mul(&b, 64))
}

const LP: u32 = 64;

/// (1 - z)(1 + z)
fn one_minus_square(z: &MpComplex) -> Option<i128> {
    product(cx(one().sub(&z.re, LP), z.im.neg()), cx(one().add(&z.re, LP), z.im.clone()))
}

/// (z - 1)(z + 1)
fn square_minus_one(z: &MpComplex) -> Option<i128> {
    one_minus_square(z)
}

/// 1 + z**2 = (1 + i z)(1 - i z)
fn one_plus_square(z: &MpComplex) -> Option<i128> {
    product(cx(one().sub(&z.im, LP), z.re.clone()), cx(one().add(&z.im, LP), z.re.neg()))
}

/// Approximate bits by which the relative error of the value exceeds the
/// relative error of the argument: log2 |z f'(z)/f(z)|, bounded from above.
enum Condition {
    Bits(i128),
    /// The computed value is zero.
    Zero,
    /// The argument is exactly a singular point of the derivative.
    Singular,
}

fn condition(f: Elementary, z: &MpComplex, v: &MpComplex) -> Condition {
    let mz = magnitude(z).unwrap_or(-(1 << 40));
    let Some(mv) = magnitude(v) else { return Condition::Zero };
    match condition_bits(f, mz, mv, z) {
        Some(bits) => Condition::Bits(bits),
        None => Condition::Singular,
    }
}

fn condition_bits(f: Elementary, mz: i128, mv: i128, z: &MpComplex) -> Option<i128> {
    Some(match f {
        F::Exp => mz.max(0),
        F::Log => -mv,
        // |cot z|, |tan z| <= 1 + 1/|v| (similarly for the hyperbolic case)
        F::Sin | F::Cos | F::Sinh | F::Cosh => mz + (1 - mv).max(1),
        // f'/f = 1/v + v
        F::Tan | F::Cot | F::Tanh | F::Coth => mz + mv.abs() + 1,
        // |tan z| <= 1 + |sec z|
        F::Sec | F::Csc | F::Sech | F::Csch => mz + (mv + 1).max(1),
        F::Asin | F::Acos => mz - one_minus_square(z)? / 2 - mv,
        F::Atan | F::Acot => mz - one_plus_square(z)? - mv,
        F::Asinh => mz - one_plus_square(z)? / 2 - mv,
        F::Acosh => mz - square_minus_one(z)? / 2 - mv,
        F::Atanh | F::Acoth => mz - one_minus_square(z)? - mv,
        F::Asec | F::Acsc => -square_minus_one(z)? / 2 - mv,
        F::Asech => -one_minus_square(z)? / 2 - mv,
        F::Acsch => -one_plus_square(z)? / 2 - mv,
    } + 2)
}

/// `f(arg)` to `p` bits: the argument is evaluated again with more bits
/// while the condition of `f` at the argument, or a value that cancels to
/// zero, shows that the working bits do not determine the result.
pub(crate) fn evaluate(f: Elementary, arg: &Expr, p: u32) -> Option<MpComplex> {
    let target = p + 10;
    let zero_limit = (4 * p + 512) as i128;
    let mut q = target;
    for _ in 0..12 {
        let z = crate::evalf_mp::ev(arg, q)?;
        // No value at a rounded argument that sits on a singularity (an
        // argument within 2^-q of 1 for atanh): the exact argument may lie
        // off it, so retry with more bits up to the cancellation bound.
        let Some(v) = value(f, &z, q) else {
            if q as i128 >= zero_limit {
                return None;
            }
            q = (2 * q + 32).min(zero_limit as u32);
            continue;
        };
        let need = match condition(f, &z, &v) {
            Condition::Bits(bits) => target as i128 + bits.max(0),
            // A finite value at a branch point of an exact argument.
            Condition::Singular => return Some(v.round_to(p)),
            Condition::Zero => {
                // An exactly cancelling value: retry with more bits up to a
                // bound, then keep it (SymPy returns such values as zero).
                if q as i128 >= zero_limit {
                    return Some(v);
                }
                (2 * q as i128 + 32).min(zero_limit)
            }
        };
        // log(z) off the real axis: the real part ln|z| is far below |log z|
        // near |z| = 1 and gets its own relative accuracy, as mpmath's
        // componentwise log does for an exact z. A real part that keeps
        // cancelling (|z| = 1 exactly) is bounded like an exact zero.
        let need = match (f, magnitude(&v)) {
            (F::Log, Some(mv)) if !v.im.is_zero() && !v.re.is_zero() && v.re.mag() < mv => {
                let real = target as i128 + (mv - v.re.mag()) + 2;
                if real > need && real > q as i128 && q as i128 >= zero_limit {
                    return Some(v.round_to(p));
                }
                need.max(real.min(zero_limit))
            }
            _ => need,
        };
        if need <= q as i128 {
            return Some(v.round_to(p));
        }
        if need > PRECISION_LIMIT {
            return None;
        }
        q = need as u32 + 10;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bf(s: &str) -> BigFloat {
        BigFloat::from_decimal(s, 200).unwrap()
    }

    fn show(v: &MpComplex) -> String {
        format!("{} {}", v.re.to_string_dps(15, false), v.im.to_string_dps(15, false))
    }

    #[test]
    fn tiny_arguments_keep_relative_accuracy() {
        let t = bf("1e-20");
        for f in [F::Sin, F::Sinh, F::Tan, F::Tanh, F::Asin, F::Atan, F::Asinh, F::Atanh] {
            let v = value(f, &real(t.clone()), 60).unwrap();
            assert_eq!(v.re.to_string_dps(12, false), "1.00000000000e-20", "{f:?}");
        }
        let l = value(F::Log, &real(one().add(&t, 200)), 60).unwrap();
        assert_eq!(l.re.to_string_dps(12, false), "1.00000000000e-20");
    }

    #[test]
    fn values_on_branch_cuts_follow_mpmath() {
        let x = real(bf("2.5"));
        assert_eq!(show(&value(F::Asin, &x, 60).unwrap()), "1.57079632679490 -1.56679923697241");
        let acos = value(F::Acos, &x, 60).unwrap();
        assert!(acos.re.is_zero());
        assert_eq!(acos.im.to_string_dps(15, false), "1.56679923697241");
        assert_eq!(show(&value(F::Atanh, &x, 60).unwrap()), "0.423648930193602 -1.57079632679490");
        assert_eq!(show(&value(F::Acosh, &real(bf("-2.5")), 60).unwrap()), "1.56679923697241 3.14159265358979");
    }
}
