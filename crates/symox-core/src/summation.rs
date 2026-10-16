//! Symbolic summation: `summation(f, (k, a, b))`.
//!
//! Handles constants, polynomials (Faulhaber's formula via Bernoulli
//! numbers), geometric terms, polynomial times geometric terms, zeta values
//! and small numeric ranges.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::diff::diff;
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{abs, undefined};
use crate::integrate::poly_coeffs;
use crate::number::Number;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sub};
use crate::relational::{eq, lt, piecewise};
use crate::subs::subs1;
use crate::binders::depends_on as has;

/// Bernoulli numbers `B_0 .. B_n` (with `B_1 = -1/2`).
pub fn bernoulli_numbers(n: usize) -> Vec<BigRational> {
    // Akiyama–Tanigawa gives B_n^+ (B_1 = +1/2)
    let mut a: Vec<BigRational> = Vec::with_capacity(n + 1);
    let mut out = Vec::with_capacity(n + 1);
    for m in 0..=n {
        a.push(BigRational::new(BigInt::one(), BigInt::from(m as u64 + 1)));
        for j in (1..=m).rev() {
            let t = (&a[j - 1] - &a[j]) * BigRational::from_integer(BigInt::from(j as u64));
            a[j - 1] = t;
        }
        out.push(a[0].clone());
    }
    if n >= 1 {
        out[1] = -out[1].clone();
    }
    out
}

/// Bernoulli number `B_n` as an expression.
pub fn bernoulli(n: usize) -> Expr {
    Expr::num(Number::from_ratio(bernoulli_numbers(n).pop().unwrap()))
}

fn binom(n: u64, k: u64) -> BigInt {
    let mut r = BigInt::one();
    for i in 0..k {
        r = r * BigInt::from(n - i) / BigInt::from(i + 1);
    }
    r
}

/// `F_p(n) = sum_{k=1}^{n} k^p` as a polynomial in `n`.
fn faulhaber(p: usize, n: &Expr) -> Expr {
    if p == 0 {
        return n.clone();
    }
    let b = bernoulli_numbers(p);
    let mut terms = Vec::with_capacity(p + 1);
    for (j, bj) in b.iter().enumerate() {
        // B_j^+ (B_1 = +1/2)
        let bj = if j == 1 { -bj.clone() } else { bj.clone() };
        if bj.is_zero() {
            continue;
        }
        let c = BigRational::from_integer(binom(p as u64 + 1, j as u64)) * bj
            / BigRational::from_integer(BigInt::from(p as u64 + 1));
        terms.push(mul2(&Expr::num(Number::from_ratio(c)), &n.powi((p + 1 - j) as i64)));
    }
    add(&terms)
}

/// `Sum(f, (k, a, b))` placeholder.
pub fn sum_unevaluated(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Expr {
    crate::binders::held(Func::Sum,f,vec![vec![k.clone(),a.clone(),b.clone()]]).expect("valid summation limits")
}

/// `sum_{k=a}^{b} f`.
pub fn summation(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Expr {
    try_summation(f, k, a, b).unwrap_or_else(|| sum_unevaluated(f, k, a, b))
}

pub fn try_summation(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    if matches!(a.kind(),Kind::Const(c) if c.is_infinite()) || b.is_constant(Constant::NegInfinity) {return None;}
    if b.is_constant(Constant::Infinity) {
        return infinite_sum(f, k, a);
    }
    if !has(f, k) {
        return Some(mul2(f, &add(&[b.clone(), neg(a), Expr::one()])));
    }
    // Karr convention: a backwards sum is the negative sum over its interior.
    // Decide only an exact integer difference, never approximate endpoint order.
    if let Some(distance)=sub(b,a).as_i64() {
        if distance == -1 {return Some(Expr::zero());}
        if distance < -1 {
            return try_summation(f,k,&add2(b,&Expr::one()),&sub(a,&Expr::one())).map(|value|neg(&value));
        }
        if distance <= 1999 && poly_coeffs(f,k).is_none() {
            return Some(add(&(0..=distance).map(|offset|subs1(f,k,&add2(a,&Expr::int(offset)))).collect::<Vec<_>>()));
        }
    }
    // constant factors (kept outside like SymPy)
    if let Kind::Mul(..) = f.kind() {
        let (c, rest): (Vec<Expr>, Vec<Expr>) = f.args().into_iter().partition(|x| !has(x, k));
        if !c.is_empty() {
            return Some(mul2(&mul(&c), &try_summation(&mul(&rest), k, a, b)?));
        }
    }
    // polynomial in k
    if let Some(c) = poly_coeffs(f, k) {
        let am1 = sub(a, &Expr::one());
        let mut terms = Vec::new();
        for (p, cp) in c.iter().enumerate() {
            if cp.is_zero() {
                continue;
            }
            let s = sub(&faulhaber(p, b), &faulhaber(p, &am1));
            terms.push(mul2(cp, &s));
        }
        return Some(expand(&add(&terms)));
    }
    // linearity
    if let Kind::Add(..) = f.kind() {
        let mut v = Vec::new();
        for t in f.args() {
            v.push(try_summation(&t, k, a, b)?);
        }
        return Some(add(&v));
    }
    // c * q**k * poly(k)
    if let Some((poly, q)) = poly_times_geometric(f, k) {
        return geometric_sum(&poly, &q, k, a, b);
    }
    None
}

/// Split `f = P(k) * q**k` (q free of k).
fn poly_times_geometric(f: &Expr, k: &Expr) -> Option<(Expr, Expr)> {
    let factors = match f.kind() {
        Kind::Mul(..) => f.args(),
        _ => vec![f.clone()],
    };
    let mut q: Option<Expr> = None;
    let mut rest = Vec::new();
    for fac in factors {
        match fac.kind() {
            Kind::Pow(base, e) if !has(base, k) && has(e, k) => {
                // base**(alpha*k + beta)
                let alpha = diff(e, k);
                if has(&alpha, k) {
                    return None;
                }
                let beta = expand(&sub(e, &mul2(&alpha, k)));
                if has(&beta, k) {
                    return None;
                }
                let qq = pow(base, &alpha);
                q = Some(match q {
                    None => qq,
                    Some(p) => mul2(&p, &qq),
                });
                rest.push(pow(base, &beta));
            }
            Kind::Func(Func::Exp, args) if has(&args[0], k) => {
                let alpha = diff(&args[0], k);
                if has(&alpha, k) {
                    return None;
                }
                let beta = expand(&sub(&args[0], &mul2(&alpha, k)));
                let qq = crate::functions::exp(&alpha);
                q = Some(match q {
                    None => qq,
                    Some(p) => mul2(&p, &qq),
                });
                rest.push(crate::functions::exp(&beta));
            }
            _ => rest.push(fac),
        }
    }
    let poly = mul(&rest);
    poly_coeffs(&poly, k)?;
    Some((poly, q?))
}

/// sum_{k=a}^{b} P(k) q^k via (q d/dq)^p of the geometric sum.
fn geometric_sum(poly: &Expr, q: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    let c = poly_coeffs(poly, k)?;
    let qs = Expr::symbol("_q");
    let g = div(
        &sub(&pow(&qs, a), &pow(&qs, &add2(b, &Expr::one()))),
        &sub(&Expr::one(), &qs),
    );
    let mut terms = Vec::new();
    let mut d = g.clone();
    for (p, cp) in c.iter().enumerate() {
        if p > 0 {
            d = mul2(&qs, &diff(&d, &qs));
        }
        if !cp.is_zero() {
            terms.push(mul2(cp, &d));
        }
    }
    let s = add(&terms);
    let numeric_q = crate::functions::is_numeric_value(q);
    let r = subs1(&s, &qs, q);
    let r = if numeric_q && c.len() > 1 { crate::simplify::together(&expand(&r)) } else { r };
    if numeric_q {
        if q.is_one() {
            return None;
        }
        return Some(if c.len() > 1 { tidy(&r) } else { r });
    }
    // symbolic ratio: the q = 1 case separately, like SymPy
    let at_one = try_summation(poly, k, a, b)?;
    Some(piecewise(vec![(at_one, eq(q, &Expr::one())), (r, Expr::boolean(true))]))
}

fn tidy(e: &Expr) -> Expr {
    let ex = expand(e);
    if crate::traverse::node_count(&ex) <= crate::traverse::node_count(e) {
        ex
    } else {
        e.clone()
    }
}

fn infinite_sum(f: &Expr, k: &Expr, a: &Expr) -> Option<Expr> {
    if !has(f, k) {
        return if f.is_zero() { Some(Expr::zero()) } else { None };
    }
    // geometric with |q| < 1
    if let Some((poly, q)) = poly_times_geometric(f, k) {
        let c = poly_coeffs(&poly, k)?;
        let qs = Expr::symbol("_q");
        let mut g = div(&pow(&qs, a), &sub(&Expr::one(), &qs));
        let mut terms = Vec::new();
        for (p, cp) in c.iter().enumerate() {
            if p > 0 {
                g = mul2(&qs, &diff(&g, &qs));
            }
            if !cp.is_zero() {
                terms.push(mul2(cp, &g));
            }
        }
        let s = subs1(&add(&terms), &qs, &q);
        if crate::functions::is_numeric_value(&q) {
            let qv = crate::evalf::eval_c64(&q)?;
            if qv.norm() < 1.0 {
                return Some(tidy(&crate::simplify::together(&s)));
            }
            return Some(Expr::oo());
        }
        return Some(piecewise(vec![
            (s, lt(&abs(&q), &Expr::one())),
            (sum_unevaluated(f, k, a, &Expr::oo()), Expr::boolean(true)),
        ]));
    }
    // zeta(p) for 1/k^p
    if let Kind::Pow(base, e) = f.kind() {
        if base == k && a.is_one() {
            if let Some(n) = e.as_i64() {
                if n < -1 && n % 2 == 0 {
                    let p = (-n) as usize;
                    // zeta(2m) = (-1)^(m+1) B_2m (2 pi)^2m / (2 (2m)!)
                    let b = bernoulli(p);
                    let mut fact = BigInt::one();
                    for i in 2..=p as u64 {
                        fact *= i;
                    }
                    let sign = if (p / 2) % 2 == 1 { 1 } else { -1 };
                    let c = mul(&[
                        Expr::int(sign),
                        b,
                        pow(&Expr::int(2), &Expr::int(p as i64 - 1)),
                        Expr::num(Number::from_ratio(BigRational::new(BigInt::one(), fact))),
                    ]);
                    return Some(mul2(&c, &Expr::pi().powi(p as i64)));
                }
                if n < -1 {
                    return Some(undefined("zeta", vec![Expr::int(-n)]));
                }
                if n >= -1 {
                    return Some(Expr::oo());
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn s(f: &str, a: &str, b: &str) -> String {
        let k = Expr::symbol("k");
        summation(&parse_str(f).unwrap(), &k, &parse_str(a).unwrap(), &parse_str(b).unwrap()).to_string()
    }

    #[test]
    fn sums() {
        assert_eq!(s("k", "1", "n"), "n**2/2 + n/2");
        assert_eq!(s("k**2", "1", "n"), "n**3/3 + n**2/2 + n/6");
        assert_eq!(s("k**3", "0", "n"), "n**4/4 + n**3/2 + n**2/4");
        assert_eq!(s("2*k+1", "0", "n-1"), "n**2");
        assert_eq!(s("k", "1", "10"), "55");
        assert_eq!(s("k**4", "1", "n"), "n**5/5 + n**4/2 + n**3/3 - n/30");
        assert_eq!(s("x", "1", "n"), "n*x");
        assert_eq!(s("1/k**2", "1", "oo"), "pi**2/6");
        assert_eq!(s("(1/2)**k", "0", "oo"), "2");
        assert_eq!(s("2**k", "0", "n"), "2**(n + 1) - 1");
        assert_eq!(s("k**2", "m", "n"), "-m**3/3 + m**2/2 - m/6 + n**3/3 + n**2/2 + n/6");
    }
}
