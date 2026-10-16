//! Resultants and discriminants (subresultant PRS over `ZZ[others]`).

use super::gcd::{content_in, prem_mz};
use super::int::Int;
use super::mpoly::MPoly;

/// Resultant of `a` and `b` with respect to `x_v` (Cohen, Algorithm
/// 3.3.7); the result is free of `x_v`.
pub fn resultant_mz(a: &MPoly<Int>, b: &MPoly<Int>, v: usize) -> MPoly<Int> {
    let n = a.nvars;
    if a.is_zero() || b.is_zero() {
        return MPoly::zero(n);
    }
    let (mut a, mut b) = (a.clone(), b.clone());
    let mut s = 1i64;
    if a.degree(v) < b.degree(v) {
        std::mem::swap(&mut a, &mut b);
        if a.degree(v) % 2 == 1 && b.degree(v) % 2 == 1 {
            s = -1;
        }
    }
    let (da0, db0) = (a.degree(v) as u32, b.degree(v) as u32);
    if db0 == 0 {
        // res(a, c) = c^deg(a)
        return b.pow(da0).scale(&Int::S(s));
    }
    let ca = content_in(&a, v);
    let cb = content_in(&b, v);
    let a1 = a.div_exact(&ca).unwrap();
    let b1 = b.div_exact(&cb).unwrap();
    let t = ca.pow(db0).mul(&cb.pow(da0));
    let (mut a, mut b) = (a1, b1);
    let mut g = MPoly::<Int>::one(n);
    let mut h = MPoly::<Int>::one(n);
    loop {
        let (da, db) = (a.degree(v), b.degree(v));
        let delta = (da - db) as u32;
        if da % 2 == 1 && db % 2 == 1 {
            s = -s;
        }
        let r = prem_mz(&a, &b, v);
        a = b;
        let den = g.mul(&h.pow(delta));
        b = r.div_exact(&den).expect("subresultant division");
        g = a.lc_in(v);
        // h = h^(1-delta) g^delta
        h = if delta == 0 {
            h
        } else {
            g.pow(delta).div_exact(&h.pow(delta - 1)).expect("subresultant h")
        };
        if b.is_zero() {
            return MPoly::zero(n);
        }
        if b.degree(v) == 0 {
            break;
        }
    }
    let da = a.degree(v) as u32;
    let lb = b.lc_in(v);
    let hh = if da == 0 {
        h.mul(&MPoly::one(n))
    } else {
        lb.pow(da).div_exact(&h.pow(da - 1)).expect("subresultant final")
    };
    hh.mul(&t).scale(&Int::S(s))
}

/// Discriminant with respect to `x_v`.
pub fn discriminant_mz(f: &MPoly<Int>, v: usize) -> MPoly<Int> {
    let n = f.nvars;
    let d = f.degree(v);
    if d <= 0 {
        return MPoly::zero(n);
    }
    let d = d as u32;
    let r = resultant_mz(f, &f.derivative(v), v);
    let lc = f.lc_in(v);
    let q = r.div_exact(&lc).expect("discriminant division");
    if (d * (d - 1) / 2) % 2 == 1 {
        q.neg()
    } else {
        q
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn res() {
        let x = MPoly::<Int>::var(0, 2);
        let y = MPoly::<Int>::var(1, 2);
        let one = MPoly::<Int>::one(2);
        // res_x(x^2 - 1, x - 2) = 3
        let r = resultant_mz(&x.pow(2).sub(&one), &x.sub(&one.scale(&Int::S(2))), 0);
        assert_eq!(r, one.scale(&Int::S(3)));
        // res_x(x^2 + y, x y - 1) = y^3 + 1
        let r = resultant_mz(&x.pow(2).add(&y), &x.mul(&y).sub(&one), 0);
        assert_eq!(r, y.pow(3).add(&one));
        // disc(a x^2 + b x + c) with a = 1: b^2 - 4c  (y = b, use constant c = 3)
        let f = x.pow(2).add(&y.mul(&x)).add(&one.scale(&Int::S(3)));
        assert_eq!(discriminant_mz(&f, 0), y.pow(2).sub(&one.scale(&Int::S(12))));
    }
}
