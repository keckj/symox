//! Double-double arithmetic (~106 bit significand) used where double
//! precision series suffer from cancellation.

use std::ops::{Add, Div, Mul, Neg, Sub};
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct DD {
    pub hi: f64,
    pub lo: f64,
}

#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let e = (a - (s - bb)) + (b - bb);
    (s, e)
}

#[inline]
fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

#[inline]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

impl DD {
    pub const ZERO: DD = DD { hi: 0.0, lo: 0.0 };
    pub const ONE: DD = DD { hi: 1.0, lo: 0.0 };

    #[inline]
    pub fn new(x: f64) -> DD {
        DD { hi: x, lo: 0.0 }
    }

    #[inline]
    pub fn from_parts(hi: f64, lo: f64) -> DD {
        let (h, l) = quick_two_sum(hi, lo);
        DD { hi: h, lo: l }
    }

    #[inline]
    pub fn to_f64(self) -> f64 {
        self.hi + self.lo
    }

    #[inline]
    pub fn add_f64(self, b: f64) -> DD {
        let (s, e) = two_sum(self.hi, b);
        let e = e + self.lo;
        let (h, l) = quick_two_sum(s, e);
        DD { hi: h, lo: l }
    }

    #[inline]
    pub fn mul_f64(self, b: f64) -> DD {
        let (p, e) = two_prod(self.hi, b);
        let e = e + self.lo * b;
        let (h, l) = quick_two_sum(p, e);
        DD { hi: h, lo: l }
    }

    pub fn div_f64(self, b: f64) -> DD {
        let q1 = self.hi / b;
        let r = self - DD::new(b) * DD::new(q1);
        let q2 = r.hi / b;
        let r = r - DD::new(b) * DD::new(q2);
        let q3 = r.hi / b;
        DD::from_parts(q1, q2).add_f64(q3)
    }

    pub fn recip(self) -> DD {
        DD::ONE / self
    }

    pub fn abs(self) -> DD {
        if self.hi < 0.0 {
            -self
        } else {
            self
        }
    }

    pub fn sqr(self) -> DD {
        self * self
    }

    pub fn sqrt(self) -> DD {
        if self.hi <= 0.0 {
            return DD::new(self.hi.sqrt());
        }
        let x = 1.0 / self.hi.sqrt();
        let ax = self.hi * x;
        let d = self - DD::new(ax).sqr();
        DD::new(ax).add_f64(d.hi * (x * 0.5))
    }

    pub fn powi(self, n: i64) -> DD {
        let mut base = self;
        let mut k = n.unsigned_abs();
        let mut acc = DD::ONE;
        while k > 0 {
            if k & 1 == 1 {
                acc = acc * base;
            }
            k >>= 1;
            if k > 0 {
                base = base.sqr();
            }
        }
        if n < 0 {
            acc.recip()
        } else {
            acc
        }
    }

    /// `e^x` to double-double accuracy.
    pub fn exp(self) -> DD {
        if self.hi > 709.8 {
            return DD::new(f64::INFINITY);
        }
        if self.hi < -745.2 {
            return DD::ZERO;
        }
        let c = consts();
        let k = (self.hi / std::f64::consts::LN_2).round();
        let r = self - c.ln2 * DD::new(k);
        // r / 2^9, Taylor, then square 9 times
        let r = r.mul_f64(1.0 / 512.0);
        let mut term = r;
        let mut sum = r;
        let mut i = 2.0;
        loop {
            term = (term * r).div_f64(i);
            sum = sum + term;
            if term.hi.abs() < 1e-34 {
                break;
            }
            i += 1.0;
            if i > 30.0 {
                break;
            }
        }
        // (1 + s)^2 - 1 = s (2 + s)
        for _ in 0..9 {
            sum = sum * sum.add_f64(2.0);
        }
        let v = sum.add_f64(1.0);
        DD { hi: v.hi * 2f64.powi(k as i32), lo: v.lo * 2f64.powi(k as i32) }
    }

    /// Natural logarithm (Newton correction of the double result).
    pub fn ln(self) -> DD {
        if self.hi <= 0.0 {
            return DD::new(f64::NAN);
        }
        let l0 = DD::new(self.hi.ln());
        // l = l0 + x e^{-l0} - 1
        let t = self * (-l0).exp();
        l0 + t.add_f64(-1.0)
    }
}

impl Add for DD {
    type Output = DD;
    #[inline]
    fn add(self, b: DD) -> DD {
        let (s, e) = two_sum(self.hi, b.hi);
        let (t, f) = two_sum(self.lo, b.lo);
        let e = e + t;
        let (s, e) = quick_two_sum(s, e);
        let e = e + f;
        let (h, l) = quick_two_sum(s, e);
        DD { hi: h, lo: l }
    }
}

impl Sub for DD {
    type Output = DD;
    #[inline]
    fn sub(self, b: DD) -> DD {
        self + (-b)
    }
}

impl Neg for DD {
    type Output = DD;
    #[inline]
    fn neg(self) -> DD {
        DD { hi: -self.hi, lo: -self.lo }
    }
}

impl Mul for DD {
    type Output = DD;
    #[inline]
    fn mul(self, b: DD) -> DD {
        let (p, e) = two_prod(self.hi, b.hi);
        let e = e + (self.hi * b.lo + self.lo * b.hi);
        let (h, l) = quick_two_sum(p, e);
        DD { hi: h, lo: l }
    }
}

impl Div for DD {
    type Output = DD;
    fn div(self, b: DD) -> DD {
        let q1 = self.hi / b.hi;
        let r = self - b.mul_f64(q1);
        let q2 = r.hi / b.hi;
        let r = r - b.mul_f64(q2);
        let q3 = r.hi / b.hi;
        DD::from_parts(q1, q2).add_f64(q3)
    }
}

/// Constants in double-double precision, computed once with the
/// arbitrary precision code.
pub struct Consts {
    pub pi: DD,
    pub ln2: DD,
    pub euler: DD,
    /// `ln(2 pi) / 2`
    pub half_ln_2pi: DD,
    pub sqrt_pi: DD,
}

fn bf_to_dd(x: &crate::bigfloat::BigFloat) -> DD {
    let hi = x.to_f64();
    let rest = x.sub(&crate::bigfloat::BigFloat::from_f64(hi), 140);
    DD { hi, lo: rest.to_f64() }
}

pub fn consts() -> &'static Consts {
    static C: OnceLock<Consts> = OnceLock::new();
    C.get_or_init(|| {
        let p = 140;
        let pi = crate::bigfloat::pi(p);
        let ln2 = crate::bigfloat::ln2(p);
        let eg = crate::evalf_mp::euler_gamma(p);
        let two_pi = pi.ldexp(1);
        let hl = two_pi.ln(p).unwrap().ldexp(-1);
        let sp = pi.sqrt(p).unwrap();
        Consts {
            pi: bf_to_dd(&pi),
            ln2: bf_to_dd(&ln2),
            euler: bf_to_dd(&eg),
            half_ln_2pi: bf_to_dd(&hl),
            sqrt_pi: bf_to_dd(&sp),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let x = DD::new(2.0).sqrt();
        let y = x * x - DD::new(2.0);
        assert!(y.hi.abs() < 1e-30);
        let e = DD::ONE.exp();
        assert!((e.hi - std::f64::consts::E).abs() < 1e-16);
        let l = e.ln() - DD::ONE;
        assert!(l.hi.abs() < 1e-30, "{l:?}");
        let c = consts();
        assert!((c.pi.hi - std::f64::consts::PI).abs() == 0.0);
        assert!(c.pi.lo.abs() > 0.0 && c.pi.lo.abs() < 1e-15);
    }
}
