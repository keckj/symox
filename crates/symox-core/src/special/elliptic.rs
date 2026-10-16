//! Complete and incomplete elliptic integrals, Lambert W.

use super::*;
use crate::functions::{exp, gamma, log};
use crate::ops::{add, div, mul2, neg, pow, sqrt, sub};

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        EllipticK => eval_k(&a[0]),
        EllipticE => {
            if a.len() == 1 {
                eval_e1(&a[0])
            } else {
                eval_e2(&a[0], &a[1])
            }
        }
        EllipticF => eval_f(&a[0], &a[1]),
        EllipticPi => {
            if a.len() == 2 {
                eval_pi2(&a[0], &a[1])
            } else {
                eval_pi3(&a[0], &a[1], &a[2])
            }
        }
        LambertW => eval_lambertw(a),
        _ => None,
    }
}

fn any_inf(m: &Expr) -> bool {
    match m.kind() {
        Kind::Const(c) => c.is_infinite(),
        _ => extract_i(m).is_some_and(|t| is_inf(&t) || is_neg_inf(&t)),
    }
}

fn eval_k(m: &Expr) -> Option<Expr> {
    if m.is_zero() {
        return Some(div(&Expr::pi(), &Expr::int(2)));
    }
    if *m == rat(1, 2) {
        return Some(div(
            &mul2(&Expr::int(8), &pow(&Expr::pi(), &rat(3, 2))),
            &gamma(&rat(-1, 4)).powi(2),
        ));
    }
    if m.is_one() {
        return Some(Expr::zoo());
    }
    if m.is_minus_one() {
        return Some(div(
            &gamma(&rat(1, 4)).powi(2),
            &mul2(&Expr::int(4), &sqrt(&mul2(&Expr::int(2), &Expr::pi()))),
        ));
    }
    if any_inf(m) {
        return Some(Expr::zero());
    }
    None
}

fn eval_e1(m: &Expr) -> Option<Expr> {
    if m.is_zero() {
        return Some(div(&Expr::pi(), &Expr::int(2)));
    }
    if m.is_one() {
        return Some(Expr::one());
    }
    if is_inf(m) {
        return Some(mul2(&Expr::i(), &Expr::oo()));
    }
    if is_neg_inf(m) {
        return Some(Expr::oo());
    }
    if m.is_constant(Constant::ComplexInfinity) {
        return Some(Expr::zoo());
    }
    None
}

/// `2 z / pi` if it is an integer.
fn half_periods(z: &Expr) -> Option<Expr> {
    let k = div(&mul2(&Expr::int(2), z), &Expr::pi());
    if is_integer(&k) == Some(true) {
        Some(k)
    } else {
        None
    }
}

fn eval_e2(z: &Expr, m: &Expr) -> Option<Expr> {
    if z.is_zero() {
        return Some(Expr::zero());
    }
    if m.is_zero() {
        return Some(z.clone());
    }
    if let Some(k) = half_periods(z) {
        return Some(mul2(&k, &super::elliptic_e(m)));
    }
    if is_inf(m) || is_neg_inf(m) {
        return Some(Expr::zoo());
    }
    if could_extract_minus_sign(z) {
        return Some(neg(&super::elliptic_e2(&neg(z), m)));
    }
    None
}

fn eval_f(z: &Expr, m: &Expr) -> Option<Expr> {
    if z.is_zero() {
        return Some(Expr::zero());
    }
    if m.is_zero() {
        return Some(z.clone());
    }
    if let Some(k) = half_periods(z) {
        return Some(mul2(&k, &super::elliptic_k(m)));
    }
    if is_inf(m) || is_neg_inf(m) {
        return Some(Expr::zero());
    }
    if could_extract_minus_sign(z) {
        return Some(neg(&super::elliptic_f(&neg(z), m)));
    }
    None
}

fn eval_pi2(n: &Expr, m: &Expr) -> Option<Expr> {
    if n.is_zero() {
        return Some(super::elliptic_k(m));
    }
    if n.is_one() {
        return Some(Expr::zoo());
    }
    if m.is_zero() {
        return Some(div(&Expr::pi(), &mul2(&Expr::int(2), &sqrt(&sub(&Expr::one(), n)))));
    }
    if m.is_one() {
        return Some(div(&Expr::neg_oo(), &crate::functions::sign(&sub(n, &Expr::one()))));
    }
    if n == m {
        return Some(div(&super::elliptic_e(n), &sub(&Expr::one(), n)));
    }
    if is_inf(n) || is_neg_inf(n) || is_inf(m) || is_neg_inf(m) {
        return Some(Expr::zero());
    }
    None
}

fn eval_pi3(n: &Expr, z: &Expr, m: &Expr) -> Option<Expr> {
    if n.is_zero() {
        return Some(super::elliptic_f(z, m));
    }
    if let Some(k) = half_periods(z) {
        return Some(mul2(&k, &super::elliptic_pi(n, m)));
    }
    if is_inf(n) || is_neg_inf(n) || is_inf(m) || is_neg_inf(m) {
        return Some(Expr::zero());
    }
    if could_extract_minus_sign(z) {
        return Some(neg(&super::elliptic_pi3(n, &neg(z), m)));
    }
    None
}

fn eval_lambertw(a: &[Expr]) -> Option<Expr> {
    let x = &a[0];
    let k = a.get(1);
    if let Some(k) = k {
        if k.is_zero() {
            return Some(super::lambertw(x));
        }
    }
    let k0 = k.is_none();
    let e = Expr::e();
    if k0 {
        if x.is_zero() {
            return Some(Expr::zero());
        }
        if x.is_constant(Constant::E) {
            return Some(Expr::one());
        }
        if *x == neg(&exp(&Expr::minus_one())) {
            return Some(Expr::minus_one());
        }
        if *x == neg(&div(&log(&Expr::int(2)), &Expr::int(2))) {
            return Some(neg(&log(&Expr::int(2))));
        }
        // W(c log c) = log c for rational c >= 1
        if let Some(c) = match_w_log_w(x) {
            return Some(log(&c));
        }
        if *x == neg(&div(&Expr::pi(), &Expr::int(2))) {
            return Some(mul2(&Expr::i(), &div(&Expr::pi(), &Expr::int(2))));
        }
        if *x == exp(&add(&[Expr::one(), e.clone()])) {
            return Some(e);
        }
        if is_inf(x) {
            return Some(Expr::oo());
        }
        // W(e^(1/n)/n) = 1/n
        if let Some(n) = match_exp_inv(x) {
            return Some(pow(&n, &Expr::minus_one()));
        }
        return None;
    }
    let k = k.unwrap();
    if is_zero(k) == Some(false) && x.is_zero() {
        return Some(Expr::neg_oo());
    }
    if k.is_minus_one() {
        if *x == neg(&div(&Expr::pi(), &Expr::int(2))) {
            return Some(neg(&mul2(&Expr::i(), &div(&Expr::pi(), &Expr::int(2)))));
        }
        if *x == neg(&exp(&Expr::minus_one())) {
            return Some(Expr::minus_one());
        }
        if *x == mul2(&Expr::int(-2), &exp(&Expr::int(-2))) {
            return Some(Expr::int(-2));
        }
    }
    None
}

/// `x == c*log(c)` with `c` a rational number `>= 1` (or `c = 1/e`).
fn match_w_log_w(x: &Expr) -> Option<Expr> {
    if let Kind::Mul(c, f) = x.kind() {
        if f.len() == 1 && f[0].1.is_one() {
            if let Kind::Func(Func::Log, a) = f[0].0.kind() {
                if let Kind::Num(n) = a[0].kind() {
                    if n == c && !n.is_float() && n.cmp_value(&Number::ONE) != Some(std::cmp::Ordering::Less) {
                        return Some(a[0].clone());
                    }
                }
            }
        }
    }
    None
}

/// `x == e^(1/n)/n` for a positive integer `n`.
fn match_exp_inv(x: &Expr) -> Option<Expr> {
    for n in 2..=12 {
        let ne = Expr::int(n);
        let cand = div(&exp(&pow(&ne, &Expr::minus_one())), &ne);
        if *x == cand {
            return Some(ne);
        }
    }
    None
}
