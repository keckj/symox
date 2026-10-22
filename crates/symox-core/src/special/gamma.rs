//! Automatic evaluation of the gamma function family.

use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive};

use super::integers as ints;
use super::*;
use crate::functions::{cos, cot, erf, erfc, exp, gamma, log, loggamma, sin};
use crate::ops::sqrt;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};

fn numq(r: num_rational::BigRational) -> Expr {
    Expr::num(Number::from_ratio(r))
}

fn num(n: BigInt) -> Expr {
    Expr::num(Number::from_bigint(n))
}

fn euler_gamma() -> Expr {
    Expr::constant(Constant::EulerGamma)
}

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        Digamma => Some(super::polygamma(&Expr::zero(), &a[0])),
        Trigamma => Some(super::polygamma(&Expr::one(), &a[0])),
        Polygamma => eval_polygamma(&a[0], &a[1]),
        Beta => {
            if a.len() == 1 {
                return Some(super::beta(&a[0], &a[0]));
            }
            let (x, y) = (&a[0], &a[1]);
            if is_int_number(x) && is_int_number(y) {
                return Some(div(&mul2(&gamma(x), &gamma(y)), &gamma(&add(&[x.clone(), y.clone()]))));
            }
            None
        }
        LowerGamma => eval_lowergamma(&a[0], &a[1]),
        UpperGamma => eval_uppergamma(&a[0], &a[1]),
        _ => None,
    }
}

/// `gamma(x)` special values (SymPy's `gamma.eval`).
pub(super) fn eval_gamma(a: &Expr) -> Option<Expr> {
    match a.kind() {
        Kind::Num(n) if n.is_integer() && !n.is_float() => {
            if !n.is_positive() {
                return Some(Expr::zoo());
            }
            let k = n.to_i64()?;
            if k > 100_000 {
                return None;
            }
            Some(num(ints::factorial(k as u64 - 1)))
        }
        Kind::Num(n) if n.is_rational() && !n.is_float() && n.denom() == Number::int(2) => {
            // gamma(k + 1/2)
            let k = n.floor().to_i64()?;
            if k.abs() > 100_000 {
                return None;
            }
            let sp = sqrt(&Expr::pi());
            let (num_, den) = if k >= 0 {
                // (2k-1)!! / 2^k
                let k = k as u64;
                (if k == 0 { BigInt::from(1) } else { ints::factorial2(2 * k - 1) }, BigInt::from(1) << (k as usize))
            } else {
                // gamma(1/2 - m) = (-2)^m / (2m-1)!! sqrt(pi)
                let m = (-k) as u64;
                let mut p = BigInt::from(1) << (m as usize);
                if m % 2 == 1 {
                    p = -p;
                }
                (p, ints::factorial2(2 * m - 1))
            };
            Some(mul2(&Expr::num(Number::rational(num_, den).unwrap()), &sp))
        }
        Kind::Const(Constant::Infinity) => Some(Expr::oo()),
        _ => {
            if !a.is_number() && is_integer(a) == Some(true) && is_nonpositive(a) == Some(true) {
                return Some(Expr::zoo());
            }
            None
        }
    }
}

/// `loggamma(x)` special values.
pub(super) fn eval_loggamma(z: &Expr) -> Option<Expr> {
    if is_integer(z) == Some(true) {
        if is_nonpositive(z) == Some(true) {
            return Some(Expr::oo());
        }
        if is_positive(z) == Some(true) {
            return Some(log(&gamma(z)));
        }
    } else if let Some(n) = as_exact(z) {
        // half integers p/2 with p > 0
        if n.denom() == Number::int(2) && n.is_positive() {
            let p = n.numer();
            let pe = Expr::num(p.clone());
            let v = mul(&[
                sqrt(&Expr::pi()),
                pow(&Expr::int(2), &sub(&Expr::one(), &pe)),
                gamma(&pe),
                pow(&gamma(&div(&add(&[pe.clone(), Expr::one()]), &Expr::int(2))), &Expr::minus_one()),
            ]);
            return Some(log(&v));
        }
    } else if !z.is_number() && crate::assumptions::is_rational(z) == Some(true) {
        // SymPy: p, q = z.as_numer_denom(); half-integral values when
        // p.is_positive and q == 2.
        let (p, q) = crate::simptools::as_numer_denom(z);
        if q == Expr::int(2) && is_positive(&p) == Some(true) {
            let v = mul(&[
                sqrt(&Expr::pi()),
                pow(&Expr::int(2), &sub(&Expr::one(), &p)),
                gamma(&p),
                pow(&gamma(&div(&add(&[p.clone(), Expr::one()]), &Expr::int(2))), &Expr::minus_one()),
            ]);
            return Some(log(&v));
        }
    }
    if is_inf(z) {
        return Some(Expr::oo());
    }
    if is_neg_inf(z) || z.is_constant(Constant::ComplexInfinity) {
        return Some(Expr::zoo());
    }
    // abs(z) is oo
    if extract_i(z).is_some_and(|t| is_inf(&t) || is_neg_inf(&t)) {
        return Some(Expr::zoo());
    }
    None
}

fn eval_polygamma(n: &Expr, z: &Expr) -> Option<Expr> {
    let nn = as_small(n);
    if nn == Some(-1) {
        return Some(sub(&loggamma(z), &div(&log(&mul2(&Expr::int(2), &Expr::pi())), &Expr::int(2))));
    }
    let nn = match nn {
        Some(k) if k >= 0 => k,
        _ => return None,
    };
    if nn > 2000 {
        return None;
    }
    if is_inf(z) {
        return Some(if nn == 0 { Expr::oo() } else { Expr::zero() });
    }
    // polygamma(0, -oo) and polygamma(0, +-I*oo)
    if nn == 0 && (is_neg_inf(z) || extract_i(z).is_some_and(|t| is_inf(&t) || is_neg_inf(&t))) {
        return Some(Expr::oo());
    }
    let zq = as_exact(z)?;
    if zq.is_integer() && !zq.is_positive() {
        return Some(Expr::zoo());
    }
    let sign_fact = |nn: i64| -> Expr {
        // (-1)^(n+1) n!
        let f = num(ints::factorial(nn as u64));
        if nn % 2 == 0 {
            neg(&f)
        } else {
            f
        }
    };
    if zq.is_one() {
        if nn == 0 {
            return Some(neg(&euler_gamma()));
        }
        return Some(mul2(&sign_fact(nn), &super::zeta(&Expr::int(nn + 1))));
    }
    if *z == rat(1, 2) && nn > 0 {
        let c = (BigInt::from(1) << ((nn + 1) as usize)) - 1;
        return Some(mul(&[sign_fact(nn), num(c), super::zeta(&Expr::int(nn + 1))]));
    }
    if zq.is_integer() {
        let m = zq.to_i64()?;
        if m > 100_000 {
            return None;
        }
        if nn == 0 {
            return Some(sub(&numq(ints::harmonic_number(m as u64 - 1, 1)), &euler_gamma()));
        }
        let h = numq(ints::harmonic_number(m as u64 - 1, (nn + 1) as u32));
        return Some(mul2(&sign_fact(nn), &sub(&super::zeta(&Expr::int(nn + 1)), &h)));
    }
    if nn == 0 {
        let (p, q) = zq.numer_denom()?;
        let q = q.to_i64()?;
        if [2, 3, 4, 6].contains(&q) {
            return digamma_rational(&p, q);
        }
    }
    None
}

/// Gauss' digamma theorem for `p/q` (small `q`).
fn digamma_rational(p: &BigInt, q: i64) -> Option<Expr> {
    let qb = BigInt::from(q);
    // z = n + p0/q with 0 < p0 < q
    let mut p0 = p % &qb;
    if p0.is_negative() {
        p0 += &qb;
    }
    let n = ((p - &p0) / &qb).to_i64()?;
    if n.abs() > 10_000 {
        return None;
    }
    let p0 = p0.to_i64()?;
    let pi = Expr::pi();
    let mut terms = vec![
        neg(&euler_gamma()),
        neg(&div(&mul2(&pi, &cot(&mul2(&rat(p0, q), &pi))), &Expr::int(2))),
        neg(&log(&Expr::int(q))),
    ];
    for k in 1..q {
        let c = cos(&mul2(&rat(2 * k * p0, q), &pi));
        if c.is_zero() {
            continue;
        }
        let s = mul2(&Expr::int(2), &sin(&mul2(&rat(k, q), &pi)));
        terms.push(mul2(&c, &log(&s)));
    }
    let part1 = add(&terms);
    let z0 = rat(p0, q);
    if n > 0 {
        let v: Vec<Expr> = (0..n).map(|k| pow(&add(&[z0.clone(), Expr::int(k)]), &Expr::minus_one())).collect();
        Some(add(&[part1, add(&v)]))
    } else if n < 0 {
        let v: Vec<Expr> = (0..(-n))
            .map(|k| pow(&sub(&z0, &Expr::int(1 + k)), &Expr::minus_one()))
            .collect();
        Some(sub(&part1, &add(&v)))
    } else {
        Some(part1)
    }
}

/// The branch information of `lowergamma.eval` / `uppergamma.eval`:
/// `z = nx*exp_polar(2*I*pi*n)`.
fn incomplete_gamma_branch(lower: bool, a: &Expr, z: &Expr) -> Option<Expr> {
    let build = |x: &Expr| if lower { super::lowergamma(a, x) } else { super::uppergamma(a, x) };
    let (nx, n) = super::hyperexpand::extract_branch_factor(z, false);
    let integer = crate::assumptions::fact(a, "integer") == Some(true);
    let two_pi_i_n = mul(&[Expr::int(2), Expr::pi(), Expr::i(), n.clone()]);
    if integer && crate::assumptions::fact(a, "positive") == Some(true) {
        let nx = super::polar::unpolarify(z, false);
        return (nx != *z).then(|| build(&nx));
    }
    if integer && crate::assumptions::fact(a, "nonpositive") == Some(true) {
        if n.is_zero() {
            return None;
        }
        let c = div(&mul2(&two_pi_i_n, &pow(&Expr::minus_one(), &neg(a))), &super::factorial(&neg(a)));
        return Some(add(&[if lower { c } else { neg(&c) }, build(&nx)]));
    }
    if n.is_zero() {
        return None;
    }
    let e = exp(&mul2(&two_pi_i_n, a));
    Some(if lower {
        mul2(&e, &build(&nx))
    } else {
        add(&[mul2(&gamma(a), &sub(&Expr::one(), &e)), mul2(&e, &build(&nx))])
    })
}

fn eval_lowergamma(a: &Expr, x: &Expr) -> Option<Expr> {
    if x.is_zero() {
        return Some(Expr::zero());
    }
    if let Some(v) = incomplete_gamma_branch(true, a, x) {
        return Some(v);
    }
    let aq = as_exact(a)?;
    if aq.is_one() {
        return Some(sub(&Expr::one(), &exp(&neg(x))));
    }
    if *a == rat(1, 2) {
        return Some(mul2(&sqrt(&Expr::pi()), &erf(&sqrt(x))));
    }
    let a2 = aq.mul(&Number::int(2));
    if !a2.is_integer() {
        return None;
    }
    let b = sub(a, &Expr::one());
    let is_int = aq.is_integer();
    if as_exact(&b)?.is_positive() {
        if is_int {
            let m = aq.to_i64()?;
            if m > 1000 {
                return None;
            }
            let fb = num(ints::factorial(m as u64 - 1));
            let s: Vec<Expr> = (0..m)
                .map(|k| div(&x.powi(k), &num(ints::factorial(k as u64))))
                .collect();
            return Some(sub(&fb, &mul(&[exp(&neg(x)), fb.clone(), add(&s)])));
        }
        // half integer > 1
        let m = aq.add(&q_half()).to_i64()?; // a + 1/2
        if m > 1000 {
            return None;
        }
        let s: Vec<Expr> = (1..m)
            .map(|k| div(&pow(x, &rat(2 * k - 1, 2)), &gamma(&rat(2 * k + 1, 2))))
            .collect();
        let inner = sub(
            &div(&super::lowergamma(&rat(1, 2), x), &sqrt(&Expr::pi())),
            &mul2(&exp(&neg(x)), &add(&s)),
        );
        return Some(mul2(&gamma(a), &inner));
    }
    if !is_int {
        // a = 1/2 - m, m >= 1
        let m = q_half().sub(aq).to_i64()?;
        if m > 1000 {
            return None;
        }
        let t1 = div(
            &mul(&[pow(&Expr::minus_one(), &Expr::int(m)), Expr::pi(), erf(&sqrt(x))]),
            &gamma(&sub(&Expr::one(), a)),
        );
        let s: Vec<Expr> = (1..=m)
            .map(|k| {
                div(
                    &mul2(&pow(x, &add(&[Expr::int(k - 1), a.clone()])), &gamma(a)),
                    &gamma(&add(&[a.clone(), Expr::int(k)])),
                )
            })
            .collect();
        return Some(add(&[t1, mul2(&exp(&neg(x)), &add(&s))]));
    }
    None
}

fn q_half() -> Number {
    Number::rational_i64(1, 2).unwrap()
}

fn eval_uppergamma(a: &Expr, z: &Expr) -> Option<Expr> {
    if is_inf(z) {
        return Some(Expr::zero());
    }
    let special_z = matches!(z.kind(), Kind::Const(Constant::NaN | Constant::Infinity)) || z.is_zero();
    if !special_z {
        if let Some(v) = incomplete_gamma_branch(false, a, z) {
            return Some(v);
        }
    }
    if let Some(aq) = as_exact(a) {
        if aq.is_zero() && is_positive(z) == Some(true) {
            return Some(neg(&super::ei(&neg(z))));
        }
        if aq.is_one() {
            return Some(exp(&neg(z)));
        }
        if *a == rat(1, 2) {
            return Some(mul2(&sqrt(&Expr::pi()), &erfc(&sqrt(z))));
        }
        let a2 = aq.mul(&Number::int(2));
        if a2.is_integer() {
            let b = aq.sub(&Number::ONE);
            let is_int = aq.is_integer();
            if b.is_positive() {
                if is_int {
                    let m = aq.to_i64()?;
                    if m > 1000 {
                        return None;
                    }
                    let fb = num(ints::factorial(m as u64 - 1));
                    let s: Vec<Expr> = (0..m)
                        .map(|k| div(&z.powi(k), &num(ints::factorial(k as u64))))
                        .collect();
                    return Some(mul(&[exp(&neg(z)), fb, add(&s)]));
                }
                // half integer > 1
                let m = aq.sub(&q_half()).to_i64()?; // a - 1/2
                if m > 1000 {
                    return None;
                }
                let s: Vec<Expr> = (0..m)
                    .map(|k| {
                        div(
                            &mul2(&gamma(&rat(-1 - 2 * k, 2)), &pow(&neg(z), &Expr::int(k))),
                            &gamma(&sub(&Expr::one(), a)),
                        )
                    })
                    .collect();
                let sgn = pow(&Expr::minus_one(), &sub(a, &rat(3, 2)));
                return Some(add(&[
                    mul2(&gamma(a), &erfc(&sqrt(z))),
                    mul(&[sgn, exp(&neg(z)), sqrt(z), add(&s)]),
                ]));
            }
            if b.is_integer() {
                // a <= 0 integer: expint(-b, z) z^(b+1)
                let bb = Expr::num(b.clone());
                return Some(mul2(
                    &super::expint(&neg(&bb), z),
                    &pow(z, &add(&[bb, Expr::one()])),
                ));
            }
            if !is_int {
                let m = q_half().sub(aq).to_i64()?;
                if m > 1000 {
                    return None;
                }
                let t1 = div(
                    &mul(&[pow(&Expr::minus_one(), &Expr::int(m)), Expr::pi(), erfc(&sqrt(z))]),
                    &gamma(&sub(&Expr::one(), a)),
                );
                let s: Vec<Expr> = (0..m)
                    .map(|k| div(&mul2(&z.powi(k), &gamma(a)), &gamma(&add(&[a.clone(), Expr::int(k + 1)]))))
                    .collect();
                return Some(sub(&t1, &mul(&[pow(z, a), exp(&neg(z)), add(&s)])));
            }
        }
    }
    if is_zero(a) == Some(true) && is_positive(z) == Some(true) {
        return Some(neg(&super::ei(&neg(z))));
    }
    if z.is_zero() {
        if let Some(aq) = as_exact(a) {
            if aq.is_positive() {
                return Some(gamma(a));
            }
        } else if is_positive(a) == Some(true) {
            return Some(gamma(a));
        }
    }
    None
}
