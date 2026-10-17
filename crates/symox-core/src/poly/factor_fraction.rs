//! Bounded structural rational splitting for fractional factor lists.
//!
//! No polynomial gcd or rational-function constructor is used here. In
//! particular (x²-1)/(x-1) retains its denominator. Sums use a structural LCM,
//! so equal literal denominator bases share powers, while differently written
//! (even algebraically equal) bases remain distinct until factorization.
use super::convert::{PResult, PolyError};
use super::factor_budget::{self, Limits};
use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, mul, pow};
use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub struct Parts {
    pub numerator: Expr,
    pub denominator: Expr,
}
fn limit(reason: &str) -> PolyError {
    PolyError::Unsupported(format!(
        "factorization resource limit: rational split {reason}"
    ))
}
fn check(e: &Expr, limits: Limits) -> PResult<()> {
    super::factor_preflight::check(e, limits)
}
fn product(xs: &[Expr], limits: Limits) -> PResult<Expr> {
    let raw = Expr::from_unevaluated_kind(Kind::Mul(
        Number::ONE,
        xs.iter().cloned().map(|x| (x, Expr::one())).collect(),
    ));
    check(&raw, limits)?;
    let result = mul(xs);
    check(&result, limits)?;
    Ok(result)
}
fn power(base: &Expr, exponent: u64, limits: Limits) -> PResult<Expr> {
    if exponent > limits.degree as u64 {
        return Err(limit("power exponent"));
    }
    let exponent = Expr::int(exponent as i64);
    check(
        &Expr::from_unevaluated_kind(Kind::Pow(base.clone(), exponent.clone())),
        limits,
    )?;
    let result = pow(base, &exponent);
    check(&result, limits)?;
    Ok(result)
}
fn numeric(n: &Number) -> Parts {
    Parts {
        numerator: Expr::num(n.numer()),
        denominator: Expr::num(n.denom()),
    }
}
fn number_bits(n: &Number) -> u64 {
    match n {
        Number::Int(_) | Number::Float(_) => 64,
        Number::BigInt(n) => n.bits(),
        Number::Rational(n) => n.numer().bits() + n.denom().bits(),
        Number::MpFloat(n, _) => n.man.bits().saturating_add(n.exp.unsigned_abs()),
    }
}
// Memoization prevents exponential traversal, but it must not retain an
// unbounded collection of distinct intermediate expressions either. Count
// shared nodes once, including the stored canonical operand vectors.
fn retain(
    e: &Expr,
    seen: &mut FxHashMap<usize, Expr>,
    storage: &mut u64,
    limits: Limits,
) -> PResult<()> {
    let mut stack = vec![e.clone()];
    while let Some(e) = stack.pop() {
        if seen.contains_key(&e.dag_id()) {
            continue;
        }
        // Hold a reference as long as its identity is accounted, preventing
        // allocator address reuse from being mistaken for an existing node.
        seen.insert(e.dag_id(), e.clone());
        factor_budget::charge(1);
        let mut bits = 256;
        match e.kind() {
            Kind::Num(n) => bits += number_bits(n),
            Kind::Add(c, xs) => {
                bits += number_bits(c);
                for (b, c) in xs {
                    bits += 128 + number_bits(c);
                    stack.push(b.clone());
                }
            }
            Kind::Mul(c, xs) => {
                bits += number_bits(c);
                for (b, x) in xs {
                    bits += 128;
                    stack.push(b.clone());
                    stack.push(x.clone());
                }
            }
            Kind::Pow(b, x) => {
                bits += 128;
                stack.push(b.clone());
                stack.push(x.clone());
            }
            Kind::Func(_, xs) => {
                bits += xs.len() as u64 * 64;
                stack.extend(xs.iter().cloned());
            }
            _ => {}
        }
        *storage = storage.saturating_add(bits);
        if *storage > limits.storage_bits || seen.len() > limits.terms {
            return Err(limit("retained DAG storage"));
        }
    }
    Ok(())
}
fn multiply(a: &Parts, b: &Parts, limits: Limits) -> PResult<Parts> {
    Ok(Parts {
        numerator: product(&[a.numerator.clone(), b.numerator.clone()], limits)?,
        denominator: product(&[a.denominator.clone(), b.denominator.clone()], limits)?,
    })
}

#[derive(Clone)]
struct Denominator {
    coefficient: Number,
    factors: Vec<(Expr, u64)>,
}
impl Denominator {
    fn new(e: &Expr) -> PResult<Self> {
        let (coefficient, factors) = match e.kind() {
            Kind::Num(n) => (n.clone(), vec![]),
            Kind::Mul(c, xs) => (c.clone(), xs.clone()),
            Kind::Pow(b, x) => (Number::ONE, vec![(b.clone(), x.clone())]),
            _ => (Number::ONE, vec![(e.clone(), Expr::one())]),
        };
        if coefficient.is_zero() {
            return Err(PolyError::Failed(
                "zero denominator in fractional factorization".into(),
            ));
        }
        let factors = factors
            .into_iter()
            .map(|(b, e)| {
                let exponent = e.as_i64().filter(|e| *e > 0).ok_or_else(|| {
                    PolyError::Unsupported("nonintegral structural denominator power".into())
                })?;
                Ok((b, exponent as u64))
            })
            .collect::<PResult<Vec<_>>>()?;
        Ok(Self {
            coefficient,
            factors,
        })
    }
    fn exponent(&self, b: &Expr) -> u64 {
        self.factors
            .iter()
            .find(|(x, _)| x == b)
            .map_or(0, |(_, k)| *k)
    }
    fn lcm(&mut self, other: &Self, limits: Limits) -> PResult<()> {
        factor_budget::charge(
            (self.factors.len() as u64)
                .saturating_mul(other.factors.len() as u64)
                .saturating_add(1),
        );
        // Denominators created by numeric() have integer coefficients.
        let gcd = self.coefficient.gcd(&other.coefficient).abs();
        let quotient = other
            .coefficient
            .div(&gcd)
            .ok_or_else(|| PolyError::Failed("invalid numeric denominator".into()))?;
        let c = product(
            &[Expr::num(self.coefficient.clone()), Expr::num(quotient)],
            limits,
        )?;
        self.coefficient = c
            .as_number()
            .ok_or_else(|| PolyError::Failed("invalid numeric denominator coefficient".into()))?
            .abs();
        for (b, k) in &other.factors {
            if let Some((_, n)) = self.factors.iter_mut().find(|(x, _)| x == b) {
                *n = (*n).max(*k);
            } else {
                self.factors.push((b.clone(), *k));
            }
        }
        Ok(())
    }
    fn quotient(&self, other: &Self, limits: Limits) -> PResult<Expr> {
        factor_budget::charge(
            (self.factors.len() as u64)
                .saturating_mul(other.factors.len() as u64)
                .saturating_add(1),
        );
        let c = self
            .coefficient
            .div(&other.coefficient)
            .ok_or_else(|| PolyError::Failed("invalid numeric denominator quotient".into()))?;
        let mut factors = vec![Expr::num(c)];
        for (b, k) in &self.factors {
            let exponent = k
                .checked_sub(other.exponent(b))
                .ok_or_else(|| PolyError::Failed("invalid denominator LCM".into()))?;
            if exponent > 0 {
                factors.push(power(b, exponent, limits)?);
            }
        }
        product(&factors, limits)
    }
    fn expression(&self, limits: Limits) -> PResult<Expr> {
        let mut factors = vec![Expr::num(self.coefficient.clone())];
        for (b, k) in &self.factors {
            factors.push(power(b, *k, limits)?);
        }
        product(&factors, limits)
    }
}
fn sum(
    parts: &[Parts],
    protected: &[Expr],
    limits: Limits,
    retained: &mut FxHashMap<usize, Expr>,
    storage: &mut u64,
) -> PResult<Parts> {
    if parts.is_empty() {
        return Ok(numeric(&Number::ZERO));
    }
    // Reserve the temporary operand arrays and cloned denominator tables
    // before allocation. They can be large even when their Expr nodes share
    // almost all of a DAG. The evolving LCM and vector capacity growth are
    // included conservatively in the four-fold factor allowance.
    let factor_count = parts
        .iter()
        .map(|p| match p.denominator.kind() {
            Kind::Mul(_, xs) => xs.len() as u64,
            _ => 1,
        })
        .fold(0u64, u64::saturating_add);
    let temporary_bits = factor_count
        .saturating_mul(512)
        .saturating_add((parts.len() as u64).saturating_mul(512));
    if temporary_bits > limits.storage_bits.saturating_sub(*storage) {
        return Err(limit("denominator table storage"));
    }
    let limits = Limits {
        storage_bits: limits.storage_bits - temporary_bits,
        ..limits
    };
    factor_budget::charge(factor_count);
    let denominators = parts
        .iter()
        .map(|p| Denominator::new(&p.denominator))
        .collect::<PResult<Vec<_>>>()?;
    let mut common = denominators[0].clone();
    for d in &denominators[1..] {
        common.lcm(d, limits)?;
    }
    let mut terms = Vec::with_capacity(parts.len());
    for (p, d) in parts.iter().zip(&denominators) {
        let step_limits = Limits {
            storage_bits: limits.storage_bits.saturating_sub(*storage) / 4,
            ..limits
        };
        let term = product(
            &[p.numerator.clone(), common.quotient(d, step_limits)?],
            step_limits,
        )?;
        retain(&term, retained, storage, limits)?;
        terms.push(term);
    }
    let raw = Expr::from_unevaluated_kind(Kind::Add(
        Number::ZERO,
        terms.iter().cloned().map(|e| (e, Number::ONE)).collect(),
    ));
    let operation_limits = Limits {
        storage_bits: limits.storage_bits.saturating_sub(*storage) / 4,
        ..limits
    };
    check(&raw, operation_limits)?;
    // This is shallow structural content extraction, not polynomial gcd or
    // cancellation against the denominator. A formal non-symbol generator
    // must not be reinterpreted as a rational factor by the legacy helper.
    let numerator = if protected.iter().any(|g| !g.is_symbol()) {
        add(&terms)
    } else {
        // The shallow legacy helper clones factor vectors and uses linear
        // searches for structural content. Charge their quadratic worst case
        // before entering it, including the numeric coefficient size.
        let mut factors = 0u64;
        let mut coefficient_bits = 1u64;
        for term in &terms {
            let (c, n) = match term.kind() {
                Kind::Mul(c, xs) => (Some(c), xs.len() as u64),
                Kind::Num(c) => (Some(c), 1),
                _ => (None, 1),
            };
            factors = factors.saturating_add(n + 1);
            if let Some(c) = c {
                coefficient_bits = coefficient_bits.max(number_bits(c));
            }
        }
        factor_budget::charge(
            factors
                .saturating_mul(factors)
                .saturating_mul(coefficient_bits / 64 + 1),
        );
        super::api::gcd_terms(&terms)
    };
    check(&numerator, limits)?;
    Ok(Parts {
        numerator,
        denominator: common.expression(limits)?,
    })
}
fn exponent(e: &Expr, limits: Limits) -> PResult<Option<i64>> {
    if !e.is_integer_number() {
        return Ok(None);
    }
    let n = e
        .as_i64()
        .filter(|n| n.unsigned_abs() <= limits.degree as u64)
        .ok_or_else(|| limit("power exponent"))?;
    Ok(Some(n))
}
fn raised(p: &Parts, n: i64, limits: Limits) -> PResult<Parts> {
    let (numerator, denominator) = if n < 0 {
        (&p.denominator, &p.numerator)
    } else {
        (&p.numerator, &p.denominator)
    };
    if denominator.is_zero() {
        return Err(PolyError::Failed(
            "zero denominator in fractional factorization".into(),
        ));
    }
    Ok(Parts {
        numerator: power(numerator, n.unsigned_abs(), limits)?,
        denominator: power(denominator, n.unsigned_abs(), limits)?,
    })
}

/// Generic expression factorization extracts numerical algebraic units before
/// imposing a coefficient domain. Pi, E and nonnumeric functions stay atoms.
pub fn scalar_content(e: &Expr, limits: Limits) -> PResult<(Expr, Expr)> {
    let scalar = |e: &Expr| match e.kind() {
        Kind::Num(_) | Kind::Const(crate::expr::Constant::I) => true,
        Kind::Add(_, xs) => xs
            .iter()
            .all(|(b, _)| b.is_constant(crate::expr::Constant::I)),
        Kind::Pow(b, x) => b.is_number() && x.is_number(),
        _ => false,
    };
    if scalar(e) {
        return Ok((e.clone(), Expr::one()));
    }
    let Kind::Mul(c, xs) = e.kind() else {
        return Ok((Expr::one(), e.clone()));
    };
    let (mut units, mut factors) = (vec![Expr::num(c.clone())], vec![]);
    for (b, x) in xs {
        let factor = if x.is_one() {
            b.clone()
        } else {
            let raw = Expr::from_unevaluated_kind(Kind::Pow(b.clone(), x.clone()));
            check(&raw, limits)?;
            pow(b, x)
        };
        if scalar(&factor) {
            units.push(factor);
        } else {
            factors.push(factor);
        }
    }
    let result = (product(&units, limits)?, product(&factors, limits)?);
    let rebuilt = product(&[result.0.clone(), result.1.clone()], limits)?;
    // Held Mul nodes may distribute a numerical coefficient only when rebuilt.
    // Both expansions have already passed the expression resource preflight.
    if rebuilt != *e && crate::expand::expand(&rebuilt) != crate::expand::expand(e) {
        return Err(PolyError::Failed(
            "scalar extraction failed reconstruction certificate".into(),
        ));
    }
    Ok(result)
}

/// Generator order of source components, before a generic square-free list
/// groups them. The decomposition itself remains global and certified.
pub fn component_generators(e: &Expr, limits: Limits) -> PResult<Vec<Expr>> {
    check(e, limits)?;
    let mut components = if matches!(e.kind(), Kind::Mul(..)) {
        e.args()
    } else {
        vec![e.clone()]
    };
    crate::printing::order::sort_factors(&mut components);
    let mut generators = vec![];
    for e in components {
        if e.is_number() {
            continue;
        }
        let p = super::polyobj::Poly::from_expr(&e, &[])?;
        for g in p.true_gens() {
            if !generators.contains(g) {
                generators.push(g.clone());
            }
        }
    }
    Ok(generators)
}

/// Split only rational arithmetic. Unknown functions/powers remain polynomial
/// atoms; explicit generators (including 1/x) are opaque formal generators.
pub fn split(expression: &Expr, protected: &[Expr], limits: Limits) -> PResult<Parts> {
    factor_budget::run(limits, || {
        check(expression, limits)?;
        let mut memo = FxHashMap::<usize, Parts>::default();
        let (mut retained, mut storage) = (FxHashMap::default(), 0u64);
        let mut stack = vec![(expression.clone(), false, 0usize)];
        while let Some((e, ready, depth)) = stack.pop() {
            if memo.contains_key(&e.dag_id()) {
                continue;
            }
            factor_budget::charge(1);
            if depth > 256 || memo.len().saturating_add(stack.len()) > limits.terms {
                return Err(limit("DAG depth or size"));
            }
            let atomic = protected.contains(&e);
            let mut children = vec![];
            if !atomic {
                match e.kind() {
                    Kind::Add(_, xs) => children.extend(xs.iter().map(|(e, _)| e.clone())),
                    Kind::Mul(_, xs) => {
                        for (b, x) in xs {
                            if exponent(x, limits)?.is_some() {
                                children.push(b.clone());
                            }
                        }
                    }
                    Kind::Pow(b, x) if exponent(x, limits)?.is_some() => children.push(b.clone()),
                    _ => {}
                }
            }
            if !ready && !children.is_empty() {
                stack.push((e.clone(), true, depth));
                for child in children {
                    if !memo.contains_key(&child.dag_id()) {
                        stack.push((child, false, depth + 1));
                    }
                }
                continue;
            }
            // Reserve headroom for numerator, denominator and their temporary
            // canonical operand vectors *before* the arithmetic allocations.
            let work_limits = Limits {
                storage_bits: limits.storage_bits.saturating_sub(storage) / 4,
                ..limits
            };
            let atom = || Parts {
                numerator: e.clone(),
                denominator: Expr::one(),
            };
            let result = if atomic {
                atom()
            } else {
                match e.kind() {
                    Kind::Num(n) => numeric(n),
                    Kind::Add(c, xs) => {
                        if (xs.len() as u64 + 1).saturating_mul(128) > work_limits.storage_bits {
                            return Err(limit("sum operand storage"));
                        }
                        let mut parts = Vec::with_capacity(xs.len() + 1);
                        if !c.is_zero() {
                            parts.push(numeric(c));
                        }
                        for (e, c) in xs {
                            let l = Limits {
                                storage_bits: limits.storage_bits.saturating_sub(storage) / 4,
                                ..limits
                            };
                            let part = multiply(&numeric(c), &memo[&e.dag_id()], l)?;
                            retain(&part.numerator, &mut retained, &mut storage, limits)?;
                            retain(&part.denominator, &mut retained, &mut storage, limits)?;
                            parts.push(part);
                        }
                        sum(&parts, protected, limits, &mut retained, &mut storage)?
                    }
                    Kind::Mul(c, xs) => {
                        let mut result = numeric(c);
                        for (b, x) in xs {
                            let p = if let Some(n) = exponent(x, limits)? {
                                raised(&memo[&b.dag_id()], n, work_limits)?
                            } else {
                                Parts {
                                    numerator: pow(b, x),
                                    denominator: Expr::one(),
                                }
                            };
                            result = multiply(&result, &p, work_limits)?;
                        }
                        result
                    }
                    Kind::Pow(b, x) => {
                        if let Some(n) = exponent(x, limits)? {
                            raised(&memo[&b.dag_id()], n, work_limits)?
                        } else {
                            atom()
                        }
                    }
                    Kind::Func(Func::Exp, args)
                        if args.len() == 1 && crate::ops::could_extract_minus_sign(&args[0]) =>
                    {
                        Parts {
                            numerator: Expr::one(),
                            denominator: crate::functions::exp(&crate::ops::neg(&args[0])),
                        }
                    }
                    _ => atom(),
                }
            };
            if result.denominator.is_zero() {
                return Err(PolyError::Failed(
                    "zero denominator in fractional factorization".into(),
                ));
            }
            check(&result.numerator, limits)?;
            check(&result.denominator, limits)?;
            retain(&result.numerator, &mut retained, &mut storage, limits)?;
            retain(&result.denominator, &mut retained, &mut storage, limits)?;
            memo.insert(e.dag_id(), result);
        }
        memo.remove(&expression.dag_id())
            .ok_or_else(|| PolyError::Failed("missing rational split root".into()))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    fn p(s: &str) -> Expr {
        parse_str(s).unwrap()
    }
    #[test]
    fn structural_lcm_preserves_uncancelled_polynomial_factors() {
        for (e, n, d) in [
            ("(x^2-1)/(x-1)", "x^2-1", "x-1"),
            ("1/x+1/x^2", "x+1", "x^2"),
            ("1/(x+1)+1/(x^2+2*x+1)", "x^2+3*x+2", "(x+1)*(x^2+2*x+1)"),
        ] {
            let r = split(&p(e), &[], Limits::default()).unwrap();
            assert_eq!(
                crate::expand::expand(&r.numerator),
                crate::expand::expand(&p(n))
            );
            assert_eq!(
                crate::expand::expand(&r.denominator),
                crate::expand::expand(&p(d))
            );
        }
    }
    #[test]
    fn shared_budget_and_formal_atoms() {
        assert!(split(
            &p("(x+1)^30"),
            &[],
            Limits {
                work: 4,
                ..Default::default()
            }
        )
        .is_err());
        let atom = p("1/x");
        let r = split(&atom, std::slice::from_ref(&atom), Limits::default()).unwrap();
        assert_eq!(r.numerator, atom);
        assert!(r.denominator.is_one());
        let held = Expr::from_unevaluated_kind(Kind::Pow(p("x+1"), Expr::int(1_000_000_000)));
        assert!(split(&held, &[], Limits::default()).is_err());
        let mut deep = p("x");
        for _ in 0..260 {
            deep = Expr::from_unevaluated_kind(Kind::Add(Number::ONE, vec![(deep, Number::ONE)]));
        }
        assert!(split(&deep, &[], Limits::default())
            .unwrap_err()
            .to_string()
            .contains("DAG depth"));
    }
    #[test]
    fn shared_dag_is_split_once_without_unfolding() {
        let mut e = p("1/x");
        for _ in 0..40 {
            e = Expr::from_unevaluated_kind(Kind::Add(
                Number::ZERO,
                vec![(e.clone(), Number::ONE), (e, Number::ONE)],
            ));
        }
        let r = split(&e, &[], Limits::default()).unwrap();
        assert_eq!(r.numerator, Expr::int(1i64 << 40));
        assert_eq!(r.denominator, p("x"));
    }
    #[test]
    fn numerical_unit_extraction_handles_held_distribution() {
        let e =
            Expr::from_unevaluated_kind(Kind::Mul(Number::Int(2), vec![(p("x+1"), Expr::one())]));
        let (c, f) = scalar_content(&e, Limits::default()).unwrap();
        assert_eq!(c, Expr::int(2));
        assert_eq!(f, p("x+1"));
    }
    #[test]
    fn denominator_tables_and_lcm_work_share_the_outer_budget() {
        let input = p("1/(x+1)+1/(x+2)+1/(x+3)+1/(x+4)");
        assert!(split(
            &input,
            &[],
            Limits {
                storage_bits: 10_000,
                ..Default::default()
            }
        )
        .is_err());
        // A permitted input still reconstructs with a single structural LCM.
        let parts = split(&input, &[], Limits::default()).unwrap();
        assert_eq!(
            crate::expand::expand(&parts.denominator),
            crate::expand::expand(&p("(x+1)*(x+2)*(x+3)*(x+4)"))
        );
    }
}
