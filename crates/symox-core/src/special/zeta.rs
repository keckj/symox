//! Automatic evaluation of the zeta function family.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;

use super::integers as ints;
use super::*;
use crate::functions::log;
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    match s {
        SpecialFn::Zeta => eval_zeta(&a[0], a.get(1)),
        SpecialFn::DirichletEta => eval_eta(&a[0]),
        SpecialFn::Polylog => eval_polylog(&a[0], &a[1]),
        _ => None,
    }
}

/// `zeta(2n)` as a rational multiple of `pi**(2n)`.
pub(crate) fn zeta_even(n2: i64) -> Option<Expr> {
    // zeta(2n) = (-1)^(n+1) B_{2n} (2 pi)^{2n} / (2 (2n)!)
    if n2 <= 0 || n2 % 2 != 0 || n2 > 4000 {
        return None;
    }
    let b = ints::bernoulli_number(n2 as usize);
    let two_pow = BigInt::from(1) << (n2 as usize);
    let mut c = b * BigRational::from_integer(two_pow)
        / BigRational::from_integer(ints::factorial(n2 as u64) * 2);
    if (n2 / 2) % 2 == 0 {
        c = -c;
    }
    Some(mul2(&Expr::num(Number::from_ratio(c.abs())), &Expr::pi().powi(n2)))
}

fn eval_zeta(s: &Expr, a: Option<&Expr>) -> Option<Expr> {
    if let Some(a) = a {
        if a.is_one() {
            return Some(super::zeta(s));
        }
    }
    if s.is_one() {
        return Some(Expr::zoo());
    }
    if is_inf(s) {
        return Some(Expr::one());
    }
    if let Some(a) = a {
        if is_inf(a) {
            return Some(Expr::zero());
        }
    }
    let sint = as_small(s);
    let one = Expr::one();
    let aa = a.cloned().unwrap_or_else(|| one.clone());
    if let Some(k) = sint {
        if k <= 0 {
            // bernoulli(1-s, a) / (s - 1)
            if 1 - k > 2000 {
                return None;
            }
            let b = super::bernoulli_poly(&Expr::int(1 - k), &aa);
            return Some(div(&b, &Expr::int(k - 1)));
        }
    }
    if a.is_none() {
        if let Some(k) = sint {
            if k % 2 == 0 {
                return zeta_even(k);
            }
        }
        return None;
    }
    let a = a.unwrap();
    if let (Some(k), Some(m)) = (sint, as_small(a)) {
        if m > 0 {
            if m > 100_000 {
                return None;
            }
            let h = super::harmonic2(&Expr::int(m - 1), s);
            let _ = k;
            return Some(sub(&super::zeta(s), &h));
        }
    }
    if let Some(m) = as_small(a) {
        if m <= 0 && (is_integer(s) == Some(false) || is_nonpositive(s) == Some(false)) {
            return Some(Expr::nan());
        }
    }
    None
}

fn eval_eta(s: &Expr) -> Option<Expr> {
    if s.is_one() {
        return Some(log(&Expr::int(2)));
    }
    let z = super::zeta(s);
    if has_head(&z, SpecialFn::Zeta) {
        return None;
    }
    Some(mul2(&sub(&Expr::one(), &pow(&Expr::int(2), &sub(&Expr::one(), s))), &z))
}

fn eval_polylog(s: &Expr, z: &Expr) -> Option<Expr> {
    if z.is_number() || z.is_constant(Constant::I) || crate::functions::is_numeric_value(z) {
        if z.is_one() {
            return Some(super::zeta(s));
        }
        if z.is_minus_one() {
            return Some(neg(&super::dirichlet_eta(s)));
        }
        if z.is_zero() {
            return Some(Expr::zero());
        }
        if as_small(s) == Some(2) {
            if let Some(v) = dilog_table(z) {
                return Some(v);
            }
        }
    }
    // SymPy's `z.equals(1)` also accepts a float one.
    if z.as_number().is_some_and(|n| n.is_float() && n.sub(&crate::number::Number::int(1)).is_zero()) {
        return Some(super::zeta(s));
    }
    if s.is_zero() {
        return Some(div(z, &sub(&Expr::one(), z)));
    }
    if s.is_minus_one() {
        return Some(div(z, &sub(&Expr::one(), z).powi(2)));
    }
    None
}

fn dilog_table(z: &Expr) -> Option<Expr> {
    let pi2 = Expr::pi().powi(2);
    let log2 = log(&Expr::int(2));
    let s5 = sqrt(&Expr::int(5));
    let half = rat(1, 2);
    if *z == half {
        return Some(sub(&div(&pi2, &Expr::int(12)), &div(&log2.powi(2), &Expr::int(2))));
    }
    if *z == Expr::int(2) {
        return Some(sub(&div(&pi2, &Expr::int(4)), &mul(&[Expr::i(), Expr::pi(), log2])));
    }
    let g = mul2(&half, &sub(&s5, &Expr::one())); // (sqrt(5) - 1)/2
    let gp = mul2(&half, &add(&[s5.clone(), Expr::one()])); // (sqrt(5) + 1)/2
    if *z == neg(&g) {
        return Some(add(&[neg(&div(&pi2, &Expr::int(15))), div(&log(&g).powi(2), &Expr::int(2))]));
    }
    if *z == neg(&gp) {
        return Some(sub(&neg(&div(&pi2, &Expr::int(10))), &log(&gp).powi(2)));
    }
    if *z == mul2(&half, &sub(&Expr::int(3), &s5)) {
        return Some(sub(&div(&pi2, &Expr::int(15)), &log(&g).powi(2)));
    }
    if *z == g {
        return Some(sub(&div(&pi2, &Expr::int(10)), &log(&g).powi(2)));
    }
    None
}
