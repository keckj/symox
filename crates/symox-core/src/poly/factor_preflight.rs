//! A DAG-aware estimate before expression-to-polynomial expansion. The later
//! checked factor engine still validates actual arithmetic and storage. This
//! first boundary prevents held powers from allocating their expansion first.
use super::convert::{PResult, PolyError};
use super::factor_budget::Limits;
use crate::expr::{Expr, Kind};
use crate::number::Number;
use rustc_hash::FxHashMap;

#[derive(Clone)]
struct Size {
    degree: u64,
    terms: u64,
    bits: u64,
    atoms: Vec<Expr>,
}
fn limit(s: &str) -> PolyError {
    PolyError::Unsupported(format!(
        "factorization resource limit before expansion: {s}"
    ))
}
fn number_bits(n: &Number) -> u64 {
    match n {
        Number::Int(n) => 64 - n.unsigned_abs().leading_zeros() as u64,
        Number::BigInt(n) => n.bits(),
        Number::Rational(n) => n.numer().bits().max(n.denom().bits()),
        Number::Float(_) => 1075,
        Number::MpFloat(n, _) => n.man.bits().saturating_add(u64::try_from(n.exp.unsigned_abs()).unwrap_or(u64::MAX)),
    }
}
fn children(e: &Expr) -> Vec<Expr> {
    // Do not call Expr::args(): Mul::args constructs powers and could evaluate
    // a held giant numeric power before this boundary has inspected it.
    match e.kind() {
        Kind::Add(_, xs) => xs.iter().map(|(e, _)| e.clone()).collect(),
        Kind::Mul(_, xs) => xs
            .iter()
            .flat_map(|(a, b)| [a.clone(), b.clone()])
            .collect(),
        Kind::Pow(a, b) => vec![a.clone(), b.clone()],
        Kind::Func(_, xs) => xs.clone(),
        Kind::Derivative(e, xs) => std::iter::once(e.clone())
            .chain(xs.iter().map(|(e, _)| e.clone()))
            .collect(),
        Kind::Rel(_, a, b) => vec![a.clone(), b.clone()],
        Kind::Piecewise(xs) => xs
            .iter()
            .flat_map(|(a, b)| [a.clone(), b.clone()])
            .collect(),
        _ => vec![],
    }
}
fn atoms(a: &mut Vec<Expr>, b: &[Expr]) -> PResult<()> {
    for e in b {
        if !a.contains(e) {
            a.push(e.clone());
            if a.len() > 8 {
                return Err(limit("more than eight polynomial/coefficient atoms"));
            }
        }
    }
    Ok(())
}
fn monomials(degree: u64, n: usize, cap: u64) -> u64 {
    let mut result = 1u128;
    for j in 1..=n as u64 {
        result = result.saturating_mul((degree + j) as u128) / j as u128;
        if result > cap as u128 {
            return cap.saturating_add(1);
        }
    }
    result as u64
}
fn power(base: &Size, b: &Expr, e: &Expr, limits: Limits) -> PResult<Size> {
    let Some(value) = e.as_number().filter(|e| e.is_integer()) else {
        // A nonintegral or symbolic power is a coefficient/generator atom,
        // not an instruction to expand to that exponent.
        return Ok(Size {
            degree: u64::from(!base.atoms.is_empty()),
            terms: 1,
            bits: base.bits,
            atoms: if base.atoms.is_empty() {
                vec![]
            } else {
                vec![crate::expr::Expr::from_unevaluated_kind(Kind::Pow(
                    b.clone(),
                    e.clone(),
                ))]
            },
        });
    };
    if b.is_one() || b.is_minus_one() || b.is_zero() {
        return Ok(Size {
            degree: 0,
            terms: 1,
            bits: 1,
            atoms: vec![],
        });
    }
    let n = value
        .to_i64()
        .map(i64::unsigned_abs)
        .ok_or_else(|| limit("power exponent"))?;
    let degree = base.degree.saturating_mul(n);
    let bits = base
        .bits
        .saturating_add(64 - base.terms.max(1).leading_zeros() as u64)
        .saturating_mul(n);
    if degree > limits.degree as u64 {
        return Err(limit("power degree"));
    }
    if bits > limits.coefficient_bits {
        return Err(limit("power coefficient size"));
    }
    let direct = if n > u32::MAX as u64 {
        u64::MAX
    } else {
        base.terms.saturating_pow(n as u32)
    };
    Ok(Size {
        degree,
        terms: direct.min(monomials(degree, base.atoms.len(), limits.terms as u64)),
        bits,
        atoms: base.atoms.clone(),
    })
}
pub fn check(e: &Expr, limits: Limits) -> PResult<()> {
    let mut memo = FxHashMap::<usize, Size>::default();
    let mut stack = vec![(e.clone(), false)];
    while let Some((e, ready)) = stack.pop() {
        if memo.contains_key(&e.dag_id()) {
            continue;
        }
        super::factor_budget::charge(1);
        if !ready {
            if memo.len().saturating_add(stack.len()) > limits.terms {
                return Err(limit("expression DAG size"));
            }
            stack.push((e.clone(), true));
            for child in children(&e) {
                if !memo.contains_key(&child.dag_id()) {
                    stack.push((child, false));
                }
            }
            continue;
        }
        let mut result = match e.kind() {
            Kind::Num(n) => Size {
                degree: 0,
                terms: 1,
                bits: number_bits(n),
                atoms: vec![],
            },
            Kind::Const(_) | Kind::Bool(_) => Size {
                degree: 0,
                terms: 1,
                bits: 1,
                atoms: vec![],
            },
            Kind::Add(c, xs) => {
                let mut size = Size {
                    degree: 0,
                    terms: u64::from(!c.is_zero()),
                    bits: number_bits(c),
                    atoms: vec![],
                };
                for (child, c) in xs {
                    let s = &memo[&child.dag_id()];
                    atoms(&mut size.atoms, &s.atoms)?;
                    size.degree = size.degree.max(s.degree);
                    size.terms = size.terms.saturating_add(s.terms);
                    size.bits = size.bits.max(s.bits.saturating_add(number_bits(c)));
                }
                size.bits = size
                    .bits
                    .saturating_add(64 - size.terms.max(1).leading_zeros() as u64);
                size
            }
            Kind::Mul(c, xs) => {
                let mut size = Size {
                    degree: 0,
                    terms: 1,
                    bits: number_bits(c),
                    atoms: vec![],
                };
                for (b, e) in xs {
                    let s = power(&memo[&b.dag_id()], b, e, limits)?;
                    atoms(&mut size.atoms, &s.atoms)?;
                    size.degree = size.degree.saturating_add(s.degree);
                    size.terms = size.terms.saturating_mul(s.terms);
                    size.bits = size
                        .bits
                        .saturating_add(s.bits)
                        .saturating_add(64 - s.terms.max(1).leading_zeros() as u64);
                }
                size
            }
            Kind::Pow(b, e) => power(&memo[&b.dag_id()], b, e, limits)?,
            _ => Size {
                degree: 1,
                terms: 1,
                bits: 1,
                atoms: vec![e.clone()],
            },
        };
        result.terms = result.terms.min(monomials(
            result.degree,
            result.atoms.len(),
            limits.terms as u64,
        ));
        if result.degree > limits.degree as u64 {
            return Err(limit("polynomial degree"));
        }
        if result.terms > limits.terms as u64 {
            return Err(limit("expanded terms"));
        }
        if result.bits > limits.coefficient_bits {
            return Err(limit("coefficient size"));
        }
        if result
            .terms
            .saturating_mul(256 + result.bits + 32 * result.atoms.len() as u64)
            > limits.storage_bits
        {
            return Err(limit("expanded storage"));
        }
        memo.insert(e.dag_id(), result);
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn held_powers_are_checked_without_constructing_expanded_arguments() {
        let x = Expr::symbol("x");
        for b in [Expr::int(2), crate::ops::add2(&x, &Expr::one())] {
            let e = Expr::from_unevaluated_kind(Kind::Pow(b, Expr::int(1_000_000_000)));
            assert!(check(&e, Limits::default())
                .unwrap_err()
                .to_string()
                .contains("before expansion"));
        }
        let e = Expr::from_unevaluated_kind(Kind::Pow(
            crate::ops::add2(&x, &Expr::one()),
            Expr::int(100),
        ));
        assert!(check(&e, Limits::default()).is_ok());
        let mut e = x;
        for _ in 0..40 {
            e = Expr::from_unevaluated_kind(Kind::Add(
                Number::ZERO,
                vec![(e.clone(), Number::ONE), (e, Number::ONE)],
            ));
        }
        assert!(check(&e, Limits::default()).is_ok());
    }
}
