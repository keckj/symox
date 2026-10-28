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
use crate::functions::abs;
use crate::integrate::poly_coeffs;
use crate::number::Number;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sub};
use crate::relational::{eq, lt, piecewise};
use crate::subs::subs1;
use crate::binders::depends_on as has;

/// Tangent numbers `T_0 .. T_half` (Brent-Harvey): only products of
/// integers by small factors, no rational arithmetic.
fn tangent_numbers(half: usize) -> Vec<BigInt> {
    let mut t: Vec<BigInt> = vec![BigInt::zero(); half + 1];
    if half >= 1 {
        t[1] = BigInt::one();
    }
    for k in 2..=half {
        t[k] = &t[k - 1] * (k as u64 - 1);
    }
    for k in 2..=half {
        for j in k..=half {
            let carry = &t[j - 1] * ((j - k) as u64);
            t[j] *= (j - k + 2) as u64;
            t[j] += carry;
        }
    }
    t
}

/// `B_m` (with `B_1 = -1/2`) from the tangent numbers:
/// `B_2k = (-1)^(k-1) 2k T_k / (4^k (4^k - 1))`.
fn bernoulli_from_tangent(m: usize, t: &[BigInt]) -> BigRational {
    match m {
        0 => BigRational::one(),
        1 => BigRational::new(-BigInt::one(), BigInt::from(2)),
        _ if m % 2 == 1 => BigRational::zero(),
        _ => {
            let k = m / 2;
            let four_k = BigInt::one() << (2 * k);
            let r = BigRational::new(&t[k] * (2 * k as u64), &four_k * (&four_k - BigInt::one()));
            if k % 2 == 0 { -r } else { r }
        }
    }
}

/// `B_n` for a large even `n` like mpmath's `bernfrac`: the denominator
/// is known exactly (von Staudt-Clausen) and the numerator is the rounding
/// of `D * 2 n! zeta(n) / (2 pi)^n`, evaluated with enough bits to be exact.
fn bernoulli_zeta(n: usize) -> Option<BigRational> {
    use crate::bigfloat::{pi, BigFloat};
    let is_prime = |p: usize| p >= 2 && (2..).take_while(|d| d * d <= p).all(|d| p % d != 0);
    let mut den = BigInt::one();
    for d in (1..=n).filter(|d| n % d == 0) {
        if is_prime(d + 1) {
            den *= (d + 1) as u64;
        }
    }
    let nf = n as f64;
    let log2_fact = crate::special::integers::log2_factorial(nf);
    let bits = log2_fact + 1.0 - nf*(2.0*std::f64::consts::PI).log2() + den.bits() as f64;
    let bitlen = |v: u64| 64 - v.leading_zeros();
    let fact = crate::ntheory::factorial(n as u64);
    let mut guard = 64 + 2*bitlen(n as u64);
    for _ in 0..4 {
        let prec = bits.max(0.0) as u32 + guard;
        // zeta(n) = sum k^-n, truncated once k^-n is below the working precision
        let terms = (2f64.powf((prec as f64 + 8.0)/nf)).ceil() as u64 + 1;
        // k^-n ~ 2^-(n log2 k) only needs the bits that reach 2^-wp; each
        // correctly rounded term then errs by at most 2^-(wp + guard bits)
        let wp = prec + bitlen(terms) + 8;
        let mut zeta = BigFloat::one();
        for k in 2..=terms {
            let drop = (nf*(k as f64).log2()).floor();
            let pk = (wp as f64 + 8.0 - drop).max(32.0) as u32;
            let term = BigFloat::from_i64(k as i64).powi(-(n as i64), pk)?;
            zeta = zeta.add(&term, wp);
        }
        let two_pi_n = pi(wp).mul_i64(2, wp).powi(n as i64, wp)?;
        let value = BigFloat::from_int(&(&fact*&den*2u32)).mul(&zeta, wp).div(&two_pi_n, wp)?;
        let num = value.to_int_round();
        // certificate: the error analysis keeps the value within 2^-guard of
        // the integer numerator; retry with more bits if it is not within 1/8
        let gap = value.sub(&BigFloat::from_int(&num), wp).abs();
        if gap.cmp(&BigFloat::from_parts(BigInt::one(), -3i64)) == std::cmp::Ordering::Less {
            return Some(BigRational::new(if n % 4 == 0 { -num } else { num }, den));
        }
        guard *= 2;
    }
    None
}

/// Bernoulli numbers `B_0 .. B_n` (with `B_1 = -1/2`).
pub fn bernoulli_numbers(n: usize) -> Vec<BigRational> {
    let t = tangent_numbers(n / 2);
    (0..=n).map(|m| bernoulli_from_tangent(m, &t)).collect()
}

/// Bernoulli number `B_n` (with `B_1 = -1/2`).
pub fn bernoulli_rational(n: usize) -> BigRational {
    if n > 1 && n % 2 == 1 {
        BigRational::zero()
    } else if n >= 120 {
        bernoulli_zeta(n).unwrap_or_else(|| bernoulli_from_tangent(n, &tangent_numbers(n / 2)))
    } else {
        bernoulli_from_tangent(n, &tangent_numbers(n / 2))
    }
}

/// Bernoulli number `B_n` as an expression.
pub fn bernoulli(n: usize) -> Expr {
    Expr::num(Number::from_ratio(bernoulli_rational(n)))
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

/// SymPy's `piecewise._canonical_coeff`: `-2*x + 1 < 0` becomes `x > 1/2`.
pub(crate) fn canonical_coeff(r: &Expr) -> Expr {
    use crate::decompose::{as_coeff_add_head, as_coeff_mul_head};
    let r = crate::logic::canonical(r);
    let Kind::Rel(op, lhs, rhs) = r.kind() else { return r };
    if matches!(rhs.kind(), Kind::Bool(_)) {
        return r.clone();
    }
    let (b, l) = as_coeff_add_head(lhs, true);
    let (m, l) = as_coeff_mul_head(&l, true);
    let rhs = div(&sub(rhs, &b), &m);
    let op = if m.as_number().is_some_and(Number::is_negative) { crate::logic::reversed_op(*op) } else { *op };
    crate::relational::rel(op, &l, &rhs)
}

/// Is `t` a product with a KroneckerDelta factor that depends on `k`?
fn has_delta_in(t: &Expr, k: &Expr) -> bool {
    use crate::special::SpecialFn;
    let factors = match t.kind() {
        Kind::Mul(..) => t.args(),
        Kind::Func(Func::Sum, args) => return has_delta_in(&args[0], k),
        _ => vec![t.clone()],
    };
    factors.iter().any(|x| matches!(x.kind(), Kind::Func(Func::Special(SpecialFn::KroneckerDelta), _)) && has(x, k))
}

/// Sum of `f*KroneckerDelta(u, v)` where `u - v` is linear in `k` (SymPy's
/// `deltasummation`): the delta selects the single index solving `u == v`.
fn delta_sum(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    use crate::special::SpecialFn;
    let is_delta = |x: &Expr| {
        matches!(x.kind(), Kind::Func(Func::Special(SpecialFn::KroneckerDelta), args) if (args.len() == 2 || args.len() == 3) && has(x, k))
    };
    let factors = match f.kind() {
        Kind::Mul(..) => f.args(),
        _ => vec![f.clone()],
    };
    let delta = factors.iter().find(|x| is_delta(x))?.clone();
    let args = delta.args();
    let d = sub(&args[0], &args[1]);
    let c = diff(&d, k);
    if has(&c, k) || c.is_zero() {
        return None;
    }
    let r = expand(&sub(&d, &mul2(&c, k)));
    if has(&r, k) {
        return None;
    }
    let value = div(&neg(&r), &c);
    let mut rest = factors;
    let position = rest.iter().position(|x| *x == delta)?;
    rest.remove(position);
    let term = subs1(&mul(&rest), k, &value);
    // deltasummation: a delta range inside the summation limits needs no Piecewise.
    if let Some(range) = args.get(2).and_then(|r| crate::special::tuple_items(r)) {
        if range.len() == 2
            && crate::assumptions::is_nonpositive(&sub(a, &range[0])) == Some(true)
            && crate::assumptions::is_nonnegative(&sub(b, &range[1])) == Some(true)
        {
            return Some(term);
        }
    }
    // Interval(a, b).as_relational(value), as SymPy builds it
    let right = if b.is_constant(Constant::Infinity) { lt(&value, b) } else { crate::relational::le(&value, b) };
    let cond = crate::logic::and(vec![canonical_coeff(&crate::relational::le(a, &value)), canonical_coeff(&right)]);
    Some(piecewise(vec![(term, cond), (Expr::zero(), Expr::boolean(true))]))
}

pub fn try_summation(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    if matches!(a.kind(),Kind::Const(c) if c.is_infinite()) || b.is_constant(Constant::NegInfinity) {return None;}
    if let Some(value) = delta_sum(f, k, a, b) {
        return Some(value);
    }
    // A delta in k inside a held inner sum whose limits do not depend on k:
    // sum over k first (the order of finite sums is exchangeable).
    if let Kind::Func(Func::Sum, args) = f.kind() {
        if args.len() == 2 {
            if let Some(limit) = crate::binders::limits(args).pop().filter(|l| l.len() == 3) {
                if !has(&limit[1], k) && !has(&limit[2], k) && limit[0] != *k && has_delta_in(&args[0], k) {
                    let inner = summation(&args[0], k, a, b);
                    if has(&inner, k) { return None; }
                    return Some(summation(&inner, &limit[0], &limit[1], &limit[2]));
                }
            }
        }
    }
    // deltasummation of an Add: each term on its own, a term without a delta
    // stays a Sum (piecewise_fold(Add(*[deltasummation(h) for h in g.args]))).
    if let Kind::Add(..) = f.kind() {
        let terms = f.args();
        if terms.iter().any(|t| has_delta_in(t, k)) {
            return Some(add(&terms.iter().map(|t| summation(t, k, a, b)).collect::<Vec<_>>()));
        }
    }
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
    // k**(-m) from an integer start a >= 1: harmonic numbers
    if let Kind::Pow(base, exponent) = f.kind() {
        let m = exponent.as_i64().filter(|m| *m < 0);
        let start = a.as_i64().filter(|a| *a >= 1);
        if let (true, Some(m), Some(_)) = (base == k, m, start) {
            let am1 = sub(a, &Expr::one());
            let h = |x: &Expr| if m == -1 { crate::special::harmonic(x) } else { crate::special::harmonic2(x, &Expr::int(-m)) };
            return Some(sub(&h(b), &h(&am1)));
        }
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
                    return Some(crate::special::sf(crate::special::SpecialFn::Zeta, vec![Expr::int(-n)]));
                }
                if n >= -1 {
                    return Some(Expr::oo());
                }
            }
        }
    }
    rational_infinite(f, k, a)
}

/// `Σ_{k=a}^oo R(k)` for a convergent rational `R` whose partial fractions
/// have rational poles at integer distances: the simple poles telescope to
/// a finite sum and the higher ones reduce to `zeta(m)` minus a head.
fn rational_infinite(f: &Expr, k: &Expr, a: &Expr) -> Option<Expr> {
    let start = a.as_i64()?;
    let (num, den) = crate::simplify::numer_denom(&crate::simplify::together(f));
    let (pn, pd) = (poly_coeffs(&num, k)?, poly_coeffs(&den, k)?);
    if pd.len() < pn.len() + 2 {
        return None;
    }
    let parts = crate::poly::api::apart(f, Some(k)).ok()?;
    let terms = match parts.kind() {
        Kind::Add(..) => parts.args(),
        _ => vec![parts.clone()],
    };
    // (coefficient, pole shift r, order m) for c/(k + r)^m
    let mut poles: Vec<(Expr, BigRational, usize)> = Vec::new();
    for t in terms {
        let (n, d) = crate::simplify::numer_denom(&t);
        let c = poly_coeffs(&d, k)?;
        let m = c.len() - 1;
        if m == 0 || has(&n, k) {
            return None;
        }
        let lead = &c[m];
        let r = div(&c[m - 1], &mul2(&Expr::int(m as i64), lead));
        let r = r.as_number()?.to_ratio()?;
        let shifted = pow(&add2(k, &Expr::num(Number::from_ratio(r.clone()))), &Expr::int(m as i64));
        if !expand(&sub(&d, &mul2(lead, &shifted))).is_zero() {
            return None;
        }
        // a pole inside the summation range makes the sum meaningless
        let at = BigRational::from_integer(BigInt::from(start)) + &r;
        if at.is_integer() && at <= BigRational::zero() {
            return None;
        }
        poles.push((div(&n, lead), r, m));
    }
    let mut out = Vec::new();
    let simple: Vec<&(Expr, BigRational, usize)> = poles.iter().filter(|p| p.2 == 1).collect();
    if let Some(r0) = simple.iter().map(|p| p.1.clone()).min() {
        for (c, r, _) in &simple {
            let shift = r - &r0;
            if !shift.is_integer() {
                return None;
            }
            let shift = i64::try_from(shift.to_integer()).ok()?;
            for j in 0..shift {
                let at = add(&[a.clone(), Expr::num(Number::from_ratio(r0.clone())), Expr::int(j)]);
                out.push(neg(&div(c, &at)));
            }
        }
    }
    for (c, r, m) in poles.iter().filter(|p| p.2 > 1) {
        let n0 = BigRational::from_integer(BigInt::from(start)) + r;
        if !n0.is_integer() {
            return None;
        }
        let n0 = i64::try_from(n0.to_integer()).ok().filter(|n| *n >= 1)?;
        let zeta = infinite_sum(&pow(k, &Expr::int(-(*m as i64))), k, &Expr::one())?;
        let head: Vec<Expr> = (1..n0).map(|j| pow(&Expr::int(j), &Expr::int(-(*m as i64)))).collect();
        out.push(mul2(c, &sub(&zeta, &add(&head))));
    }
    Some(tidy(&add(&out)))
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
