//! Exact gcd of big integers and reduced rational arithmetic.
//!
//! num-bigint computes `gcd` with Stein's binary algorithm, which needs one
//! full-length subtraction and shift per bit of the operands: quadratic in
//! the bit length with a large constant, and pathological when one operand
//! is much smaller than the other (`gcd(a, 2**k)` with a huge odd `a`
//! subtracts 1 from `a` bit by bit). This module computes the same unique
//! non-negative gcd with:
//!
//! - the removal of the common power of two first (`gcd(a, 2**k * m)` with
//!   `m = 1` needs only the trailing zeros of `a`);
//! - one Euclidean division when the operands have very different sizes;
//! - Lehmer's algorithm (Knuth, TAOCP vol. 2, 4.5.2, Algorithm L) on 64-bit
//!   limbs otherwise: the quotient sequence is computed on the leading 63
//!   bits and applied as a unimodular 2x2 matrix, one linear pass per about
//!   30 bits of reduction. Every matrix is checked to give non-negative
//!   results that fit; any other case takes an exact Euclidean step.
//!
//! The rational helpers produce the same reduced `BigRational` (positive
//! denominator, `0/1` for zero) as `num_rational::Ratio`, with the gcds of
//! Henrici's formulas (the ones of Python's `fractions`) instead of a full
//! reduction of the unreduced result.

use std::cmp::Ordering;

use num_bigint::{BigInt, BigUint, Sign};
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    let k = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    loop {
        b >>= b.trailing_zeros();
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        b -= a;
        if b == 0 {
            return a << k;
        }
    }
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        if a <= u64::MAX as u128 && b <= u64::MAX as u128 {
            return gcd_u64(a as u64, b as u64) as u128;
        }
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn normalize(v: &mut Vec<u64>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn bit_length(v: &[u64]) -> u64 {
    match v.last() {
        None => 0,
        Some(top) => 64 * (v.len() as u64 - 1) + (64 - top.leading_zeros() as u64),
    }
}

/// The 64 bits of `v` starting at bit `shift`.
fn bits_at(v: &[u64], shift: u64) -> u64 {
    let w = (shift / 64) as usize;
    let r = shift % 64;
    let lo = v.get(w).copied().unwrap_or(0);
    if r == 0 {
        return lo;
    }
    let hi = v.get(w + 1).copied().unwrap_or(0);
    (lo >> r) | (hi << (64 - r))
}

fn cmp_limbs(a: &[u64], b: &[u64]) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| a.iter().rev().cmp(b.iter().rev()))
}

fn to_u128(v: &[u64]) -> u128 {
    v.first().copied().unwrap_or(0) as u128 | (v.get(1).copied().unwrap_or(0) as u128) << 64
}

fn from_limbs(v: &[u64]) -> BigUint {
    let mut digits = Vec::with_capacity(2 * v.len());
    for &limb in v {
        digits.push(limb as u32);
        digits.push((limb >> 32) as u32);
    }
    BigUint::new(digits)
}

/// `out = a*u + b*v` for cofactors of opposite signs (or zero) with
/// `|a|, |b| < 2**63`. Returns false (leaving `out` unspecified) when the
/// result is negative; it never overflows `i128`.
fn lin_comb(u: &[u64], v: &[u64], a: i64, b: i64, out: &mut Vec<u64>) -> bool {
    out.clear();
    let (a, b) = (a as i128, b as i128);
    let mut carry: i128 = 0;
    for (i, &ui) in u.iter().enumerate() {
        let vi = v.get(i).copied().unwrap_or(0);
        let acc = a * ui as i128 + b * vi as i128 + carry;
        out.push(acc as u64);
        carry = acc >> 64;
    }
    if carry < 0 {
        return false;
    }
    if carry > 0 {
        out.push(carry as u64);
        if carry >> 64 != 0 {
            return false;
        }
    }
    normalize(out);
    true
}

/// One exact Euclidean step: `(u, v) <- (v, u mod v)`, `v != 0`.
fn division_step(u: &mut Vec<u64>, v: &mut Vec<u64>) {
    let r = if v.len() == 1 {
        let d = v[0] as u128;
        let mut rem: u128 = 0;
        for &limb in u.iter().rev() {
            rem = ((rem << 64) | limb as u128) % d;
        }
        if rem == 0 { Vec::new() } else { vec![rem as u64] }
    } else {
        let mut r = (from_limbs(u) % from_limbs(v)).to_u64_digits();
        normalize(&mut r);
        r
    };
    std::mem::swap(u, v);
    *v = r;
}

/// Lehmer gcd of `u >= v`, both normalized limb vectors.
fn lehmer(mut u: Vec<u64>, mut v: Vec<u64>) -> BigUint {
    let mut t = Vec::with_capacity(u.len() + 1);
    let mut w = Vec::with_capacity(u.len() + 1);
    loop {
        if v.is_empty() {
            return from_limbs(&u);
        }
        if u.len() <= 2 {
            return BigUint::from(gcd_u128(to_u128(&u), to_u128(&v)));
        }
        let shift = bit_length(&u) - 63;
        let mut x = bits_at(&u, shift) as i128;
        let mut y = bits_at(&v, shift) as i128;
        let (mut a, mut b, mut c, mut d) = (1i128, 0i128, 0i128, 1i128);
        loop {
            let (d1, d2) = (y + c, y + d);
            let (n1, n2) = (x + a, x + b);
            if d1 <= 0 || d2 <= 0 || n1 < 0 || n2 < 0 {
                break;
            }
            let q = n1 / d1;
            if q != n2 / d2 {
                break;
            }
            let nc = a - q * c;
            a = c;
            c = nc;
            let nd = b - q * d;
            b = d;
            d = nd;
            let ny = x - q * y;
            x = y;
            y = ny;
        }
        let fits = |z: i128| z.unsigned_abs() < 1u128 << 63;
        let alternating = |p: i128, q: i128| p == 0 || q == 0 || (p < 0) != (q < 0);
        if b == 0
            || !(fits(a) && fits(b) && fits(c) && fits(d))
            || !alternating(a, b)
            || !alternating(c, d)
            || !lin_comb(&u, &v, a as i64, b as i64, &mut t)
            || !lin_comb(&u, &v, c as i64, d as i64, &mut w)
        {
            division_step(&mut u, &mut v);
            continue;
        }
        // The matrix [[a, b], [c, d]] is unimodular, so gcd(t, w) = gcd(u, v).
        std::mem::swap(&mut u, &mut t);
        std::mem::swap(&mut v, &mut w);
        if cmp_limbs(&u, &v) == Ordering::Less {
            std::mem::swap(&mut u, &mut v);
        }
    }
}

/// Non-negative gcd of two natural numbers.
pub fn gcd_biguint(a: &BigUint, b: &BigUint) -> BigUint {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }
    if let (Some(x), Some(y)) = (a.to_u128(), b.to_u128()) {
        return BigUint::from(gcd_u128(x, y));
    }
    if a.is_one() || b.is_one() {
        return BigUint::one();
    }
    // gcd(2**i * a', 2**j * b') = 2**min(i, j) * gcd(a', b') for odd a', b'.
    let (za, zb) = (a.trailing_zeros().unwrap_or(0), b.trailing_zeros().unwrap_or(0));
    let k = za.min(zb);
    let a = if za > 0 { a >> za } else { a.clone() };
    let b = if zb > 0 { b >> zb } else { b.clone() };
    if a.is_one() || b.is_one() {
        return BigUint::one() << k;
    }
    let (mut u, mut v) = if a >= b { (a, b) } else { (b, a) };
    // Very unbalanced operands: one exact division first.
    if u.bits() > v.bits() + 64 {
        let r = &u % &v;
        if r.is_zero() {
            return v << k;
        }
        u = v;
        v = r;
    }
    let mut ul = u.to_u64_digits();
    let mut vl = v.to_u64_digits();
    normalize(&mut ul);
    normalize(&mut vl);
    let g = lehmer(ul, vl);
    if k > 0 { g << k } else { g }
}

/// Non-negative gcd of two integers (`gcd(0, 0) = 0`).
pub fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    if let (Some(x), Some(y)) = (a.to_i64(), b.to_i64()) {
        return BigInt::from(gcd_u64(x.unsigned_abs(), y.unsigned_abs()));
    }
    BigInt::from_biguint(Sign::Plus, gcd_biguint(a.magnitude(), b.magnitude()))
}

/// Non-negative lcm of two integers (`0` if either is zero).
pub fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    if a.is_zero() || b.is_zero() {
        return BigInt::zero();
    }
    let g = gcd(a, b);
    let q = if g.is_one() { a.abs() } else { a.abs() / &g };
    q * b.abs()
}

/// `a / b` for a `b` known to divide `a`.
fn div_exact(a: &BigInt, b: &BigInt) -> BigInt {
    if b.is_one() { a.clone() } else { a / b }
}

/// `n/d` in lowest terms with a positive denominator; panics if `d == 0`
/// (like `Ratio::new`).
pub fn ratio(n: BigInt, d: BigInt) -> BigRational {
    assert!(!d.is_zero(), "denominator == 0");
    if n.is_zero() {
        return BigRational::zero();
    }
    let g = gcd(&n, &d);
    let (mut n, mut d) = if g.is_one() { (n, d) } else { (n / &g, d / &g) };
    if d.is_negative() {
        n = -n;
        d = -d;
    }
    BigRational::new_raw(n, d)
}

/// Sum of two reduced rationals, reduced.
pub fn ratio_add(x: &BigRational, y: &BigRational) -> BigRational {
    let (na, da) = (x.numer(), x.denom());
    let (nb, db) = (y.numer(), y.denom());
    if da.is_one() && db.is_one() {
        return BigRational::from_integer(na + nb);
    }
    if na.is_zero() {
        return y.clone();
    }
    if nb.is_zero() {
        return x.clone();
    }
    let g = gcd(da, db);
    if g.is_one() {
        // gcd(na*db + nb*da, da*db) = 1 for coprime reduced denominators.
        return BigRational::new_raw(na * db + nb * da, da * db);
    }
    let s = div_exact(da, &g);
    let t = na * div_exact(db, &g) + nb * &s;
    if t.is_zero() {
        return BigRational::zero();
    }
    let g2 = gcd(&t, &g);
    if g2.is_one() {
        return BigRational::new_raw(t, s * db);
    }
    BigRational::new_raw(div_exact(&t, &g2), s * div_exact(db, &g2))
}

/// Difference of two reduced rationals, reduced.
pub fn ratio_sub(x: &BigRational, y: &BigRational) -> BigRational {
    ratio_add(x, &-y)
}

/// Product of two reduced rationals, reduced.
pub fn ratio_mul(x: &BigRational, y: &BigRational) -> BigRational {
    let (na, da) = (x.numer(), x.denom());
    let (nb, db) = (y.numer(), y.denom());
    if na.is_zero() || nb.is_zero() {
        return BigRational::zero();
    }
    if da.is_one() && db.is_one() {
        return BigRational::from_integer(na * nb);
    }
    let g1 = gcd(na, db);
    let g2 = gcd(nb, da);
    let n = div_exact(na, &g1) * div_exact(nb, &g2);
    let d = div_exact(da, &g2) * div_exact(db, &g1);
    BigRational::new_raw(n, d)
}

/// Quotient of two reduced rationals, reduced; panics if `y == 0`.
pub fn ratio_div(x: &BigRational, y: &BigRational) -> BigRational {
    assert!(!y.numer().is_zero(), "division by zero");
    let inv = if y.numer().is_negative() {
        BigRational::new_raw(-y.denom(), -y.numer())
    } else {
        BigRational::new_raw(y.denom().clone(), y.numer().clone())
    };
    ratio_mul(x, &inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_integer::Integer;

    /// Deterministic pseudo-random big integers (xorshift).
    struct Gen(u64);
    impl Gen {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn big(&mut self, limbs: usize) -> BigUint {
            let mut digits = Vec::new();
            for _ in 0..limbs {
                let x = self.next();
                digits.push(x as u32);
                digits.push((x >> 32) as u32);
            }
            BigUint::new(digits)
        }
    }

    #[test]
    fn matches_binary_gcd() {
        let mut g = Gen(0x9e37_79b9_7f4a_7c15);
        for case in 0..400 {
            let (la, lb, lc) = (1 + case % 23, 1 + (case * 7) % 19, case % 11);
            let common = g.big(lc) | BigUint::one();
            let mut a = g.big(la) * &common;
            let mut b = g.big(lb) * &common;
            if case % 5 == 0 {
                a <<= case % 130;
                b <<= (case * 3) % 70;
            }
            if case % 17 == 0 {
                b = BigUint::one() << (case * 13);
            }
            if case % 19 == 0 {
                b = a.clone() * 3u32;
            }
            assert_eq!(gcd_biguint(&a, &b), a.gcd(&b), "case {case}");
            assert_eq!(gcd_biguint(&b, &a), a.gcd(&b), "case {case}");
            let (sa, sb) = (BigInt::from(a.clone()), -BigInt::from(b.clone()));
            assert_eq!(gcd(&sa, &sb), sa.gcd(&sb));
            assert_eq!(lcm(&sa, &sb), sa.lcm(&sb));
        }
        assert_eq!(gcd_biguint(&BigUint::zero(), &BigUint::from(6u32)), BigUint::from(6u32));
        assert_eq!(gcd(&BigInt::zero(), &BigInt::zero()), BigInt::zero());
        assert_eq!(gcd(&BigInt::from(i64::MIN), &BigInt::from(i64::MIN)), BigInt::from(i64::MIN).abs());
    }

    #[test]
    fn large_quotients_and_fibonacci_operands() {
        // Consecutive Fibonacci numbers (all quotients 1) and operands with a
        // quotient wider than a limb (forces Euclidean steps).
        let (mut f0, mut f1) = (BigUint::zero(), BigUint::one());
        for _ in 0..3000 {
            let f2 = &f0 + &f1;
            f0 = f1;
            f1 = f2;
        }
        assert_eq!(gcd_biguint(&f1, &f0), BigUint::one());
        let mut g = Gen(12345);
        let v = g.big(9);
        let u = &v * (BigUint::one() << 70u32) + g.big(3);
        assert_eq!(gcd_biguint(&u, &v), u.gcd(&v));
        let p = g.big(40) | BigUint::one();
        assert_eq!(gcd_biguint(&(&p * g.big(30)), &(&p * g.big(31))) % &p, BigUint::zero());
    }

    #[test]
    fn rational_operations_match_ratio() {
        let mut g = Gen(77);
        let mut values = Vec::new();
        for case in 0..60usize {
            let n = BigInt::from(g.big(1 + case % 4)) - BigInt::from(g.big(1 + case % 3));
            let d = BigInt::from(g.big(1 + case % 5)) * BigInt::from(1 + case as u64 % 6);
            let d = if case % 4 == 0 { -d } else { d };
            let d = if d.is_zero() { BigInt::one() } else { d };
            let n = if case % 9 == 0 { BigInt::zero() } else { n };
            assert_eq!(ratio(n.clone(), d.clone()), BigRational::new(n.clone(), d.clone()));
            values.push(BigRational::new(n, d));
        }
        values.push(BigRational::from_integer(BigInt::from(7)));
        values.push(BigRational::new(BigInt::from(5), BigInt::from(12)));
        values.push(BigRational::new(BigInt::from(-1), BigInt::from(12)));
        for x in &values {
            for y in &values {
                assert_eq!(ratio_add(x, y), x + y);
                assert_eq!(ratio_sub(x, y), x - y);
                assert_eq!(ratio_mul(x, y), x * y);
                if !y.is_zero() {
                    assert_eq!(ratio_div(x, y), x / y);
                }
            }
        }
    }
}
