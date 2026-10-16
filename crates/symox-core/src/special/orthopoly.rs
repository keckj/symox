//! Classical orthogonal polynomials.
//!
//! For a given integer degree the polynomial is returned explicitly (as
//! SymPy does, collected in powers of `x` with expanded coefficients);
//! otherwise special values and symmetries are applied.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

use super::combinatorial::poly_q;
use super::integers as ints;
use super::*;
use crate::expand::expand;
use crate::functions::{cos, exp, gamma};
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};

const MAX_DEGREE: i64 = 2000;

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        Legendre => eval_legendre(&a[0], &a[1]),
        AssocLegendre => eval_assoc_legendre(&a[0], &a[1], &a[2]),
        ChebyshevT => eval_chebyshevt(&a[0], &a[1]),
        ChebyshevU => eval_chebyshevu(&a[0], &a[1]),
        Hermite => eval_hermite(&a[0], &a[1]),
        Laguerre => eval_laguerre(&a[0], &a[1]),
        AssocLaguerre => eval_assoc_laguerre(&a[0], &a[1], &a[2]),
        Jacobi => eval_jacobi(&a[0], &a[1], &a[2], &a[3]),
        Gegenbauer => eval_gegenbauer(&a[0], &a[1], &a[2]),
        Ynm => eval_ynm(&a[0], &a[1], &a[2], &a[3]),
        _ => None,
    }
}

fn degree(n: &Expr) -> Option<i64> {
    as_small(n).filter(|k| k.abs() <= MAX_DEGREE)
}

fn qi(n: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

fn numq(r: BigRational) -> Expr {
    Expr::num(Number::from_ratio(r))
}

/// Rational coefficients of `P_n`.
pub(crate) fn legendre_coeffs(n: usize) -> Vec<BigRational> {
    // (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
    let mut p0 = vec![qi(1)];
    if n == 0 {
        return p0;
    }
    let mut p1 = vec![qi(0), qi(1)];
    for k in 1..n {
        let mut next = vec![BigRational::zero(); k + 2];
        let kk = k as i64;
        for (i, c) in p1.iter().enumerate() {
            next[i + 1] += c * qi(2 * kk + 1);
        }
        for (i, c) in p0.iter().enumerate() {
            next[i] -= c * qi(kk);
        }
        for c in next.iter_mut() {
            *c /= qi(kk + 1);
        }
        p0 = p1;
        p1 = next;
    }
    p1
}

/// Integer coefficients from a three-term recurrence
/// `p_{k+1} = alpha x p_k - beta(k) p_{k-1}`.
fn recurrence(n: usize, p0: Vec<BigInt>, p1: Vec<BigInt>, alpha: i64, beta: impl Fn(i64) -> i64) -> Vec<BigInt> {
    if n == 0 {
        return p0;
    }
    let (mut a, mut b) = (p0, p1);
    for k in 1..n {
        let mut next = vec![BigInt::zero(); b.len() + 1];
        for (i, c) in b.iter().enumerate() {
            next[i + 1] += c * alpha;
        }
        let bk = beta(k as i64);
        for (i, c) in a.iter().enumerate() {
            next[i] -= c * bk;
        }
        a = b;
        b = next;
    }
    b
}

fn zq(v: Vec<BigInt>) -> Vec<BigRational> {
    v.into_iter().map(BigRational::from_integer).collect()
}

fn chebyshevt_coeffs(n: usize) -> Vec<BigRational> {
    zq(recurrence(n, vec![BigInt::one()], vec![BigInt::zero(), BigInt::one()], 2, |_| 1))
}

fn chebyshevu_coeffs(n: usize) -> Vec<BigRational> {
    zq(recurrence(n, vec![BigInt::one()], vec![BigInt::zero(), BigInt::from(2)], 2, |_| 1))
}

fn hermite_coeffs(n: usize) -> Vec<BigRational> {
    zq(recurrence(n, vec![BigInt::one()], vec![BigInt::zero(), BigInt::from(2)], 2, |k| 2 * k))
}

fn laguerre_coeffs(n: usize) -> Vec<BigRational> {
    // L_n = sum_k C(n,k) (-1)^k x^k / k!
    let mut c = Vec::with_capacity(n + 1);
    let mut binom = BigInt::one();
    for k in 0..=n {
        let mut v = BigRational::new(binom.clone(), ints::factorial(k as u64));
        if k % 2 == 1 {
            v = -v;
        }
        c.push(v);
        binom = binom * BigInt::from((n - k) as u64) / BigInt::from(k as u64 + 1);
    }
    c
}

/// Polynomial in `x` with expression coefficients (constant first); each
/// coefficient is expanded and the result collected like SymPy's
/// `Poly.as_expr()`.
fn poly_e(c: &[Expr], x: &Expr) -> Expr {
    let terms: Vec<Expr> = c
        .iter()
        .enumerate()
        .map(|(i, k)| (i, expand(k)))
        .filter(|(_, k)| !k.is_zero())
        .map(|(i, k)| mul2(&k, &x.powi(i as i64)))
        .collect();
    add(&terms)
}

/// Rising factorial `(a)_k` as an (unexpanded) product.
fn poch(a: &Expr, k: i64) -> Expr {
    let v: Vec<Expr> = (0..k).map(|i| add(&[a.clone(), Expr::int(i)])).collect();
    mul(&v)
}

fn eval_legendre(n: &Expr, x: &Expr) -> Option<Expr> {
    if let Some(k) = degree(n) {
        if k < 0 {
            return Some(super::legendre(&Expr::int(-k - 1), x));
        }
        return Some(poly_q(&legendre_coeffs(k as usize), x));
    }
    if n.is_number() {
        return None;
    }
    if could_extract_minus_sign(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &super::legendre(n, &neg(x))));
    }
    if could_extract_minus_sign(n) && !could_extract_minus_sign(&sub(&neg(n), &Expr::one())) {
        return Some(super::legendre(&sub(&neg(n), &Expr::one()), x));
    }
    if x.is_zero() {
        let half = rat(1, 2);
        return Some(div(
            &sqrt(&Expr::pi()),
            &mul2(&gamma(&sub(&half, &div(n, &Expr::int(2)))), &gamma(&add(&[Expr::one(), div(n, &Expr::int(2))]))),
        ));
    }
    if x.is_one() {
        return Some(Expr::one());
    }
    if is_inf(x) {
        return Some(Expr::oo());
    }
    None
}

fn eval_assoc_legendre(n: &Expr, m: &Expr, x: &Expr) -> Option<Expr> {
    if could_extract_minus_sign(m) {
        let mm = neg(m);
        let f = div(
            &super::factorial(&add(&[mm.clone(), n.clone()])),
            &super::factorial(&sub(n, &mm)),
        );
        return Some(mul(&[pow(&Expr::minus_one(), &mm), f, super::assoc_legendre(n, &mm, x)]));
    }
    if m.is_zero() {
        return Some(super::legendre(n, x));
    }
    if x.is_zero() {
        let two = Expr::int(2);
        return Some(div(
            &mul2(&pow(&two, m), &sqrt(&Expr::pi())),
            &mul2(
                &gamma(&div(&sub(&sub(&Expr::one(), m), n), &two)),
                &gamma(&sub(&Expr::one(), &div(&sub(m, n), &two))),
            ),
        ));
    }
    let (nn, mm) = (degree(n)?, degree(m)?);
    if nn < 0 || mm > nn {
        return None;
    }
    // (-1)^m (1 - x^2)^(m/2) d^m/dx^m P_n(x)
    let mut c = legendre_coeffs(nn as usize);
    for _ in 0..mm {
        c = c.iter().enumerate().skip(1).map(|(i, v)| v * qi(i as i64)).collect();
    }
    if c.is_empty() {
        return Some(Expr::zero());
    }
    let p = poly_q(&c, x);
    let w = pow(&sub(&Expr::one(), &x.powi(2)), &rat(mm, 2));
    Some(mul(&[pow(&Expr::minus_one(), &Expr::int(mm)), w, p]))
}

fn eval_chebyshevt(n: &Expr, x: &Expr) -> Option<Expr> {
    if let Some(k) = degree(n) {
        return Some(poly_q(&chebyshevt_coeffs(k.unsigned_abs() as usize), x));
    }
    if n.is_number() {
        return None;
    }
    if could_extract_minus_sign(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &super::chebyshevt(n, &neg(x))));
    }
    if could_extract_minus_sign(n) {
        return Some(super::chebyshevt(&neg(n), x));
    }
    if x.is_zero() {
        return Some(cos(&mul(&[rat(1, 2), Expr::pi(), n.clone()])));
    }
    if x.is_one() {
        return Some(Expr::one());
    }
    if is_inf(x) {
        return Some(Expr::oo());
    }
    None
}

fn eval_chebyshevu(n: &Expr, x: &Expr) -> Option<Expr> {
    if let Some(k) = degree(n) {
        if k < 0 {
            if k == -1 {
                return Some(Expr::zero());
            }
            return Some(neg(&poly_q(&chebyshevu_coeffs((-k - 2) as usize), x)));
        }
        return Some(poly_q(&chebyshevu_coeffs(k as usize), x));
    }
    if n.is_number() {
        return None;
    }
    if could_extract_minus_sign(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &super::chebyshevu(n, &neg(x))));
    }
    if could_extract_minus_sign(n) {
        if n.is_minus_one() {
            return Some(Expr::zero());
        }
        let m2 = sub(&neg(n), &Expr::int(2));
        if !could_extract_minus_sign(&m2) {
            return Some(neg(&super::chebyshevu(&m2, x)));
        }
    }
    if x.is_zero() {
        return Some(cos(&mul(&[rat(1, 2), Expr::pi(), n.clone()])));
    }
    if x.is_one() {
        return Some(add(&[Expr::one(), n.clone()]));
    }
    if is_inf(x) {
        return Some(Expr::oo());
    }
    None
}

fn eval_hermite(n: &Expr, x: &Expr) -> Option<Expr> {
    if let Some(k) = degree(n) {
        if k < 0 {
            return None;
        }
        return Some(poly_q(&hermite_coeffs(k as usize), x));
    }
    if n.is_number() {
        return None;
    }
    if could_extract_minus_sign(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &super::hermite(n, &neg(x))));
    }
    if x.is_zero() {
        return Some(div(
            &mul2(&pow(&Expr::int(2), n), &sqrt(&Expr::pi())),
            &gamma(&div(&sub(&Expr::one(), n), &Expr::int(2))),
        ));
    }
    if is_inf(x) {
        return Some(Expr::oo());
    }
    None
}

fn eval_laguerre(n: &Expr, x: &Expr) -> Option<Expr> {
    if let Some(k) = degree(n) {
        if k < 0 {
            return Some(mul2(&exp(x), &super::laguerre(&Expr::int(-k - 1), &neg(x))));
        }
        return Some(poly_q(&laguerre_coeffs(k as usize), x));
    }
    if n.is_number() {
        return None;
    }
    if x.is_zero() {
        return Some(Expr::one());
    }
    if is_neg_inf(x) {
        return Some(Expr::oo());
    }
    if is_inf(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &Expr::oo()));
    }
    None
}

fn eval_assoc_laguerre(n: &Expr, alpha: &Expr, x: &Expr) -> Option<Expr> {
    if alpha.is_zero() {
        return Some(super::laguerre(n, x));
    }
    if let Some(k) = degree(n) {
        if k < 0 {
            return None;
        }
        // L_n^a(x) = sum_k (-1)^k (a+k+1)_{n-k} / ((n-k)! k!) x^k
        let c: Vec<Expr> = (0..=k)
            .map(|j| {
                let r = BigRational::new(
                    if j % 2 == 0 { BigInt::one() } else { -BigInt::one() },
                    ints::factorial((k - j) as u64) * ints::factorial(j as u64),
                );
                mul2(&numq(r), &poch(&add(&[alpha.clone(), Expr::int(j + 1)]), k - j))
            })
            .collect();
        return Some(poly_e(&c, x));
    }
    if n.is_number() {
        return None;
    }
    if x.is_zero() {
        return Some(super::binomial(&add(&[n.clone(), alpha.clone()]), alpha));
    }
    None
}

fn eval_gegenbauer(n: &Expr, a: &Expr, x: &Expr) -> Option<Expr> {
    if is_integer(n) == Some(true) && is_negative(n) == Some(true) {
        return Some(Expr::zero());
    }
    if *a == rat(1, 2) {
        return Some(super::legendre(n, x));
    }
    if a.is_one() {
        return Some(super::chebyshevu(n, x));
    }
    if a.is_minus_one() {
        return Some(Expr::zero());
    }
    if let Some(k) = degree(n) {
        if k < 0 {
            return Some(Expr::zero());
        }
        // C_n^a(x) = sum_j (-1)^j (a)_{n-j} / (j! (n-2j)!) (2x)^{n-2j}
        let mut c = vec![Expr::zero(); k as usize + 1];
        for j in 0..=(k / 2) {
            let p = (k - 2 * j) as usize;
            let mut r = BigRational::new(
                BigInt::one() << p,
                ints::factorial(j as u64) * ints::factorial(p as u64),
            );
            if j % 2 == 1 {
                r = -r;
            }
            c[p] = mul2(&numq(r), &poch(a, k - j));
        }
        return Some(poly_e(&c, x));
    }
    if n.is_number() {
        return None;
    }
    if x.is_minus_one() {
        if let Some(av) = as_exact(a) {
            if av.cmp_value(&Number::rational_i64(1, 2).unwrap()) == Some(std::cmp::Ordering::Greater) {
                return Some(Expr::zoo());
            }
        }
        return None;
    }
    if could_extract_minus_sign(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &super::gegenbauer(n, a, &neg(x))));
    }
    if x.is_zero() {
        let two = Expr::int(2);
        return Some(div(
            &mul(&[pow(&two, n), sqrt(&Expr::pi()), gamma(&add(&[a.clone(), div(n, &two)]))]),
            &mul(&[
                gamma(&div(&sub(&Expr::one(), n), &two)),
                gamma(&add(&[n.clone(), Expr::one()])),
                gamma(a),
            ]),
        ));
    }
    if x.is_one() {
        let two_a = mul2(&Expr::int(2), a);
        return Some(div(
            &gamma(&add(&[two_a.clone(), n.clone()])),
            &mul2(&gamma(&two_a), &gamma(&add(&[n.clone(), Expr::one()]))),
        ));
    }
    if is_inf(x) {
        return None;
    }
    None
}

fn eval_jacobi(n: &Expr, a: &Expr, b: &Expr, x: &Expr) -> Option<Expr> {
    let half = rat(1, 2);
    if a == b {
        if *a == rat(-1, 2) {
            return Some(mul(&[
                super::rf(&half, n),
                pow(&super::factorial(n), &Expr::minus_one()),
                super::chebyshevt(n, x),
            ]));
        }
        if a.is_zero() {
            return Some(super::legendre(n, x));
        }
        if *a == half {
            return Some(mul(&[
                super::rf(&rat(3, 2), n),
                pow(&super::factorial(&add(&[n.clone(), Expr::one()])), &Expr::minus_one()),
                super::chebyshevu(n, x),
            ]));
        }
        return Some(mul(&[
            super::rf(&add(&[a.clone(), Expr::one()]), n),
            pow(&super::rf(&add(&[mul2(&Expr::int(2), a), Expr::one()]), n), &Expr::minus_one()),
            super::gegenbauer(n, &add(&[a.clone(), half]), x),
        ]));
    }
    if *b == neg(a) {
        let one = Expr::one();
        return Some(mul(&[
            gamma(&add(&[n.clone(), a.clone(), one.clone()])),
            pow(&gamma(&add(&[n.clone(), one.clone()])), &Expr::minus_one()),
            pow(&add(&[one.clone(), x.clone()]), &div(a, &Expr::int(2))),
            pow(&sub(&one, x), &neg(&div(a, &Expr::int(2)))),
            super::assoc_legendre(n, &neg(a), x),
        ]));
    }
    if let Some(k) = degree(n) {
        if k < 0 {
            return None;
        }
        // P = 1/n! sum_m C(n,m) (a+m+1)_{n-m} (a+b+n+1)_m ((x-1)/2)^m
        let k = k as usize;
        let mut c = vec![Expr::zero(); k + 1];
        let mut acc: Vec<Vec<Expr>> = vec![Vec::new(); k + 1];
        let apb = add(&[a.clone(), b.clone(), Expr::int(k as i64 + 1)]);
        let mut binom = BigInt::one();
        for m in 0..=k {
            let coef = mul(&[
                numq(BigRational::new(binom.clone(), ints::factorial(k as u64))),
                poch(&add(&[a.clone(), Expr::int(m as i64 + 1)]), (k - m) as i64),
                poch(&apb, m as i64),
            ]);
            let coef = expand(&coef);
            // ((x-1)/2)^m = 2^-m sum_j C(m,j) x^j (-1)^(m-j)
            let mut bj = BigInt::one();
            for j in 0..=m {
                let mut r = BigRational::new(bj.clone(), BigInt::one() << m);
                if (m - j) % 2 == 1 {
                    r = -r;
                }
                acc[j].push(mul2(&numq(r), &coef));
                bj = bj * BigInt::from((m - j) as u64) / BigInt::from(j as u64 + 1);
            }
            binom = binom * BigInt::from((k - m) as u64) / BigInt::from(m as u64 + 1);
        }
        for (j, v) in acc.into_iter().enumerate() {
            c[j] = add(&v);
        }
        return Some(poly_e(&c, x));
    }
    if n.is_number() {
        return None;
    }
    if could_extract_minus_sign(x) {
        return Some(mul2(&pow(&Expr::minus_one(), n), &super::jacobi(n, b, a, &neg(x))));
    }
    if x.is_zero() {
        let one = Expr::one();
        return Some(mul(&[
            pow(&Expr::int(2), &neg(n)),
            gamma(&add(&[a.clone(), n.clone(), one.clone()])),
            pow(&mul2(&gamma(&add(&[a.clone(), one.clone()])), &super::factorial(n)), &Expr::minus_one()),
            super::hyper(&[sub(&neg(b), n), neg(n)], &[add(&[a.clone(), one])], &Expr::minus_one()),
        ]));
    }
    if x.is_one() {
        return Some(div(&super::rf(&add(&[a.clone(), Expr::one()]), n), &super::factorial(n)));
    }
    None
}

fn eval_ynm(n: &Expr, m: &Expr, theta: &Expr, phi: &Expr) -> Option<Expr> {
    let i = Expr::i();
    if could_extract_minus_sign(m) {
        let mm = neg(m);
        return Some(mul(&[
            pow(&Expr::minus_one(), &mm),
            exp(&mul(&[Expr::int(-2), i, mm.clone(), phi.clone()])),
            super::ynm(n, &mm, theta, phi),
        ]));
    }
    if could_extract_minus_sign(theta) {
        return Some(super::ynm(n, m, &neg(theta), phi));
    }
    if could_extract_minus_sign(phi) {
        let p = neg(phi);
        return Some(mul2(
            &exp(&mul(&[Expr::int(-2), i, m.clone(), p.clone()])),
            &super::ynm(n, m, theta, &p),
        ));
    }
    None
}
