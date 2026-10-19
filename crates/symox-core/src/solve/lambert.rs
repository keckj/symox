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
        if let Kind::Func(Func::Exp, args) = factor.kind() {
            if has(&args[0], x) {
                if exponential.is_some() {
                    return None;
                }
                exponential = Some(args[0].clone());
                continue;
            }
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
        let a = mul(&coefficients);
        if has(&a, x) || a.is_zero() {
            continue;
        }
        let (b, c) = affine(&argument, x)?;
        let rest = add(&terms
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, t)| t.clone())
            .collect::<Vec<_>>());
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
