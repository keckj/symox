//! Affine exponential and logarithmic equations, solved through Lambert W.
//! solve() follows SymPy's finite principal-branch convention; it includes W_-1
//! when its real domain is established, rather than claiming all complex branches.
use crate::expr::{Expr, Func, Kind};
use crate::functions::exp;
use crate::ops::{add, div, mul, mul2, neg, sub};
use crate::traverse::has;

fn affine(e: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    let c = super::poly_coeffs(&crate::expand::expand(e), x)?;
    if c.len() > 2 {
        return None;
    }
    Some((
        c.get(1).cloned().unwrap_or_else(Expr::zero),
        c.first().cloned().unwrap_or_else(Expr::zero),
    ))
}

fn factors(e: &Expr) -> Vec<Expr> {
    if e.is_mul() {
        e.args()
    } else {
        vec![e.clone()]
    }
}

/// Coefficient, affine base, positive integral power, exponential argument.
fn exponential_term(e: &Expr, x: &Expr) -> Option<(Expr, Expr, u32, Expr)> {
    let mut exponential = None;
    let mut rest = Vec::new();
    for factor in factors(e) {
        if let Some(argument) = exponential_argument(&factor, x) {
            if exponential.is_some() {
                return None;
            }
            exponential = Some(argument);
            continue;
        }
        rest.push(factor);
    }
    let (independent, dependent): (Vec<_>, Vec<_>) =
        rest.into_iter().partition(|factor| !has(factor, x));
    let coefficient = mul(&independent);
    let product = mul(&dependent);
    let (base, power) = match product.kind() {
        Kind::Pow(base, power) if power.as_i64().is_some_and(|p| p > 0 && p <= 4096) => {
            (base.clone(), power.as_i64()? as u32)
        }
        _ => (product, 1),
    };
    Some((coefficient, base, power, exponential?))
}

/// The argument `u` of a factor `exp(u)`, `b**u = exp(u*log(b))` with `b`
/// free of the unknown (the reference `_solve_lambert` rewrites such powers
/// with `exp`), or an integral power of one of them.
fn exponential_argument(factor: &Expr, x: &Expr) -> Option<Expr> {
    match factor.kind() {
        Kind::Func(Func::Exp, args) if has(&args[0], x) => Some(args[0].clone()),
        Kind::Pow(base, power) if !has(base, x) && has(power, x) => {
            Some(mul2(power, &crate::functions::log(base)))
        }
        Kind::Pow(base, power) if has(base, x) && power.as_i64().is_some() => {
            exponential_argument(base, x).map(|argument| mul2(power, &argument))
        }
        _ => None,
    }
}

fn branches(argument: &Expr) -> Option<Vec<Expr>> {
    let mut values = vec![crate::special::lambertw(argument)];
    if argument.is_zero() {
        return Some(values);
    }
    match crate::special::lambert::real(argument, -1) {
        Some(true) => values.push(crate::special::lambertw_k(argument, &Expr::minus_one())),
        None if crate::functions::is_numeric_value(argument) => return None,
        _ => {}
    }
    Some(values)
}

/// `k*(a*x+b)**n` with `k` free of the unknown, `n` a positive integer
/// greater than one and `a*x+b` affine (a single linear term is left to the
/// affine exponential form).
fn monomial(e: &Expr, x: &Expr) -> Option<(Expr, Expr, u32)> {
    if e.is_add() {
        return None;
    }
    let (independent, dependent): (Vec<_>, Vec<_>) =
        factors(e).into_iter().partition(|factor| !has(factor, x));
    let product = mul(&dependent);
    let Kind::Pow(base, power) = product.kind() else { return None };
    let n = power.as_i64().filter(|p| *p > 1 && *p <= 4096)?;
    let (a, _) = affine(base, x)?;
    if a.is_zero() {
        return None;
    }
    Some((mul(&independent), base.clone(), n as u32))
}

pub(crate) fn solve(f: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    let terms = if f.is_add() {
        f.args()
    } else {
        vec![f.clone()]
    };
    // (a*x+b)^n exp(c*x+d) = R. Taking every n-th root preserves all
    // polynomial branches before the finite Lambert-branch convention.
    for (i, term) in terms.iter().enumerate() {
        let Some((coefficient, base, power, exponent)) = exponential_term(term, x) else {
            continue;
        };
        let rest = add(&terms
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, t)| t.clone())
            .collect::<Vec<_>>());
        let (c, d) = affine(&exponent, x)?;
        if c.is_zero() {
            continue;
        }
        // k*(a*x+b)**n + C*exp(c*x+d) = 0 is k*(a*x+b)**n*exp(-c*x-d) + C = 0.
        let (coefficient, base, power, c, d, rest) = match monomial(&rest, x) {
            Some((k, monomial_base, n)) if !has(&base, x) && !has(&coefficient, x) => {
                let amplitude = mul2(&coefficient, &crate::ops::pow(&base, &Expr::int(power as i64)));
                (k, monomial_base, n, neg(&c), neg(&d), amplitude)
            }
            _ => (coefficient, base, power, c, d, rest),
        };
        if let Some((a, b)) = affine(&base, x) {
            if !a.is_zero() && !has(&rest, x) {
                let n = Expr::int(power as i64);
                let q = mul2(
                    &div(&neg(&rest), &coefficient),
                    &exp(&sub(&div(&mul2(&c, &b), &a), &d)),
                );
                let radicals = super::binomial_roots(&Expr::one(), &neg(&q), power);
                let mut result = Vec::new();
                for radical in radicals {
                    let argument = div(&mul2(&c, &radical), &mul2(&a, &n));
                    for w in branches(&argument)? {
                        result.push(sub(&div(&mul2(&n, &w), &c), &div(&b, &a)));
                    }
                }
                return Some(result);
            }
        }
        // a*x+b + C exp(c*x+d) = 0.
        if !has(&base, x) {
            let (a, b) = affine(&rest, x)?;
            if a.is_zero() {
                continue;
            }
            let amplitude = mul2(
                &coefficient,
                &crate::ops::pow(&base, &Expr::int(power as i64)),
            );
            let argument = mul2(
                &div(&mul2(&c, &amplitude), &a),
                &exp(&sub(&d, &div(&mul2(&c, &b), &a))),
            );
            return Some(
                branches(&argument)?
                    .into_iter()
                    .map(|w| sub(&neg(&div(&w, &c)), &div(&b, &a)))
                    .collect(),
            );
        }
    }
    // A log(b*x+c) + d*x+e = 0.
    for (i, term) in terms.iter().enumerate() {
        let mut logarithm = None;
        let mut coefficients = Vec::new();
        for factor in factors(term) {
            if let Kind::Func(Func::Log, args) = factor.kind() {
                if has(&args[0], x) && logarithm.is_none() {
                    logarithm = Some(args[0].clone());
                    continue;
                }
            }
            coefficients.push(factor);
        }
        let Some(argument) = logarithm else {
            continue;
        };
        let rest = add(&terms
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, t)| t.clone())
            .collect::<Vec<_>>());
        let a = mul(&coefficients);
        if has(&a, x) {
            if let Some(result) = power_logarithm(&coefficients, &argument, &rest, x) {
                return Some(result);
            }
            continue;
        }
        if a.is_zero() {
            continue;
        }
        let (b, c) = affine(&argument, x)?;
        let (d, e) = affine(&rest, x)?;
        if b.is_zero() || d.is_zero() {
            continue;
        }
        let ab = mul2(&a, &b);
        let z = mul2(
            &div(&d, &ab),
            &exp(&sub(&div(&mul2(&c, &d), &ab), &div(&e, &a))),
        );
        return Some(
            branches(&z)?
                .into_iter()
                .map(|w| sub(&div(&mul2(&a, &w), &d), &div(&c, &b)))
                .collect(),
        );
    }
    None
}

/// k*(p*x + q)**n*log(c*x + d) + rest = 0 with p*x + q proportional to
/// c*x + d and rest free of x. With L = c*x + d, L**n*log(L) = S, so
/// n*log(L)*exp(n*log(L)) = n*S and log(L) = W(n*S)/n, the form the
/// reference _solve_lambert returns (x**2*log(x) - 1 gives exp(W(2)/2)).
fn power_logarithm(coefficients: &[Expr], argument: &Expr, rest: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    if has(rest, x) {
        return None;
    }
    let (independent, dependent): (Vec<_>, Vec<_>) =
        coefficients.iter().cloned().partition(|factor| !has(factor, x));
    let k = mul(&independent);
    let product = mul(&dependent);
    let (base, power) = match product.kind() {
        Kind::Pow(base, power) if power.as_i64().is_some_and(|p| p > 0 && p <= 4096) => {
            (base.clone(), power.as_i64()?)
        }
        _ => (product, 1),
    };
    let (p, q) = affine(&base, x)?;
    let (c, d) = affine(argument, x)?;
    if p.is_zero() || c.is_zero() || k.is_zero() {
        return None;
    }
    if !crate::expand::expand(&sub(&mul2(&p, &d), &mul2(&q, &c))).is_zero() {
        return None;
    }
    let n = Expr::int(power);
    let ratio = div(&p, &c);
    let s = div(&neg(rest), &mul2(&k, &crate::ops::pow(&ratio, &n)));
    Some(
        branches(&mul2(&n, &s))?
            .into_iter()
            .map(|w| div(&sub(&exp(&div(&w, &n)), &d), &c))
            .collect(),
    )
}

/// The reference replaces a Lambert solution by nsimplify(s) when that
/// satisfies the equation and has no more operations. For a real numeric
/// solution nsimplify finds a rational, then (mpmath identify) forms such as
/// exp(r) and sqrt(r): exp(W(log(4))) is 2 and exp(W(3*exp(3))) is exp(3).
/// A candidate is taken only when it agrees with the solution to 200 bits.
pub(crate) fn simpler_rational(s: &Expr) -> Expr {
    if !crate::functions::is_numeric_value(s) || s.is_rational_number() {
        return s.clone();
    }
    let Some(value) = crate::evalf::eval_f64(s) else { return s.clone() };
    let rational = |v: f64| crate::simp::rational_approx(v, 10_000, 1e-14).map(|(p, q)| Expr::rational(p, q));
    let mut candidates = Vec::new();
    candidates.extend(rational(value));
    if value > 0.0 {
        candidates.extend(rational(value.ln()).filter(|r| !r.is_zero()).map(|r| exp(&r)));
        candidates.extend(rational(value * value).map(|r| crate::ops::sqrt(&r)));
    }
    let Some(exact) = crate::evalf_mp::eval_mp(s, 256) else { return s.clone() };
    let limit = value.abs().max(1.0).log2().ceil() as i128 - 200;
    let close = |part: &crate::bigfloat::BigFloat| part.is_zero() || part.mag() < limit;
    for candidate in candidates {
        if candidate == *s {
            continue;
        }
        let Some(approximation) = crate::evalf_mp::eval_mp(&candidate, 256) else { continue };
        let difference = exact.sub(&approximation, 256);
        if close(&difference.re) && close(&difference.im) {
            return candidate;
        }
    }
    s.clone()
}
