//! Distributional primitives for one simple real affine Dirac delta.
//!
//! For an analytic multiplier g near the real root r of a*x+b,
//! g(x)*delta(a*x+b) = g(r)*delta(a*x+b).  Consequently a primitive is
//! g(r)*Heaviside(a*x+b)/a, including negative a. Other cases retain the
//! integrator's original fallback; in particular, products of distributions
//! and singular/discontinuous multipliers are not assigned point values.
use crate::{Expr, Func, Kind};
use crate::assumptions as query;
use crate::functions::heaviside;
use crate::ops::{div, mul, neg};
use crate::subs::subs1;

thread_local! {
    static EXACT_RULE_ALLOWED: std::cell::Cell<bool> = const { std::cell::Cell::new(true) };
}

/// Scope the new exact rule without suppressing any established integration
/// method. Inexact input cannot become accidentally eligible after u-substitution.
pub(super) fn with_input<R>(expression: &Expr, calculation: impl FnOnce() -> R) -> R {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) { EXACT_RULE_ALLOWED.with(|flag| flag.set(self.0)); }
    }
    let previous = EXACT_RULE_ALLOWED.with(|flag|
        flag.replace(flag.get() && crate::evalf_mp::max_precision(expression).is_none()));
    let _restore = Restore(previous);
    calculation()
}

pub(super) fn antiderivative(expression: &Expr, variable: &Expr) -> Option<Expr> {
    if !EXACT_RULE_ALLOWED.with(|flag| flag.get()) { return None; }
    match expression.kind() {
        Kind::Func(Func::DiracDelta, _) => {},
        Kind::Mul(_, factors) if factors.iter().any(|(base, exponent)|
            exponent.is_one() && matches!(base.kind(), Kind::Func(Func::DiracDelta, _))) => {},
        _ => return None,
    }
    let factors = if expression.is_mul() { expression.args() } else { vec![expression.clone()] };
    let mut delta = None;
    let mut multiplier = Vec::new();
    for factor in factors {
        if let Kind::Func(Func::DiracDelta, args) = factor.kind() {
            if delta.is_some() || args.is_empty() || args.len() > 2 { return None; }
            let order = match args.get(1).map(|k| k.as_number().cloned()) {
                None => 0,
                Some(Some(crate::number::Number::Int(k))) if k >= 0 => k,
                Some(_) => return None,
            };
            delta = Some((args[0].clone(), order));
        } else {
            multiplier.push(factor);
        }
    }
    let (argument, order) = delta?;
    if order > 0 {
        return derivative_antiderivative(&argument, order, &multiplier, variable);
    }
    let (slope, offset) = super::linear(&argument, variable)?;
    if !query::exact_queries(|| query::is_real(&slope) == Some(true)
        && query::is_zero(&slope) == Some(false)) { return None; }
    let root = div(&neg(&offset), &slope);
    if query::exact_queries(|| query::is_real(&root)) != Some(true) { return None; }
    let multiplier = mul(&multiplier);
    if !crate::matrix::functions::scalar_substitution_preflight(&multiplier, variable, &root) {
        return None;
    }
    if !analytic_at(&multiplier, variable, &root, 0) { return None; }
    let value = subs1(&multiplier, variable, &root);
    if super::has_bad(&value) { return None; }
    Some(div(&mul(&[value, heaviside(&argument)]), &slope))
}

/// `g(x)*DiracDelta(a*x + b, n)` with `n > 0`, by repeated integration by
/// parts as in the reference `deltaintegrate`: with `p` the root and
/// `c_j = binomial(n, j)*(-1)**j*g^(j)(p)`, the primitive's largest term is
/// `c_n*Heaviside(x - p)` when `c_n != 0`, otherwise
/// `c_j*DiracDelta(x - p, n - j - 1)` for the largest `j` with `c_j != 0`.
/// The scale `1/(a**n*|a|)` of `DiracDelta(a*x + b, n)` is kept, and a lone
/// delta keeps the reference form `DiracDelta(a*x + b, n - 1)/a`.
fn derivative_antiderivative(argument: &Expr, order: i64, multiplier: &[Expr], variable: &Expr) -> Option<Expr> {
    if order > 64 { return None; }
    let (slope, offset) = super::linear(argument, variable)?;
    if !query::exact_queries(|| query::is_real(&slope) == Some(true)
        && query::is_zero(&slope) == Some(false)) { return None; }
    let root = div(&neg(&offset), &slope);
    if query::exact_queries(|| query::is_real(&root)) != Some(true) { return None; }
    let delta = |k: i64, at: &Expr| crate::functions::func(Func::DiracDelta, if k == 0 { vec![at.clone()] } else { vec![at.clone(), Expr::int(k)] });
    if multiplier.is_empty() {
        return Some(div(&delta(order - 1, argument), &slope));
    }
    let multiplier = mul(multiplier);
    if !crate::matrix::functions::scalar_substitution_preflight(&multiplier, variable, &root) {
        return None;
    }
    if !analytic_at(&multiplier, variable, &root, 0) { return None; }
    let scale = div(&Expr::one(), &mul(&[crate::ops::pow(&slope, &Expr::int(order)), crate::functions::abs(&slope)]));
    let shifted = crate::ops::sub(variable, &root);
    let mut derivatives = vec![multiplier];
    for _ in 0..order {
        let last = derivatives.last()?.clone();
        derivatives.push(crate::diff::diff(&last, variable));
    }
    for j in (0..=order).rev() {
        let value = subs1(&derivatives[j as usize], variable, &root);
        if super::has_bad(&value) { return None; }
        let coefficient = mul(&[
            Expr::int(binomial(order, j)),
            Expr::int(if j % 2 == 0 { 1 } else { -1 }),
            value,
        ]);
        // as in the reference, only a proven zero moves to a lower term
        if query::exact_queries(|| query::is_zero(&coefficient)) == Some(true) {
            continue;
        }
        let term = if j == order { heaviside(&shifted) } else { delta(order - j - 1, &shifted) };
        return Some(mul(&[scale, coefficient, term]));
    }
    Some(Expr::zero())
}

fn binomial(n: i64, k: i64) -> i64 {
    let mut c: i128 = 1;
    for i in 0..k as i128 {
        c = c * (n as i128 - i) / (i + 1);
    }
    c as i64
}

fn analytic_at(expression: &Expr, variable: &Expr, root: &Expr, depth: usize) -> bool {
    if depth > 128 { return false; }
    if matches!(expression.kind(), Kind::Bool(_) | Kind::Rel(..)
        | Kind::Func(Func::Logic(_) | Func::Set(_) | Func::Predicate(_) | Func::AppliedPredicate, _)) {
        return false;
    }
    if super::free(expression, variable) {
        return !super::has_bad(expression)
            && query::exact_queries(|| query::is_finite(expression)) != Some(false);
    }
    if expression == variable { return true; }
    match expression.kind() {
        Kind::Add(..) | Kind::Mul(..) => expression.args().iter().all(|arg|
            analytic_at(arg, variable, root, depth+1)),
        Kind::Func(Func::Exp | Func::Sin | Func::Cos | Func::Sinh | Func::Cosh, args) =>
            args.iter().all(|arg| analytic_at(arg, variable, root, depth+1)),
        Kind::Pow(base, exponent) if super::free(exponent, variable) => {
            if !analytic_at(base, variable, root, depth+1) { return false; }
            if exponent.as_number().is_some_and(|n| n.is_integer() && !n.is_negative()) {
                return true;
            }
            let value = subs1(base, variable, root);
            query::exact_queries(|| query::is_positive(&value) == Some(true)
                || query::is_integer(exponent) == Some(true) && query::is_zero(&value) == Some(false))
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::functions::{dirac_delta, exp};
    use crate::ops::{add2, mul2, sub};

    #[test]
    fn simple_support_substitutes_only_the_analytic_weight() {
        let x = Expr::symbol("x");
        let argument = sub(&mul2(&Expr::int(-2), &x), &Expr::int(-2));
        let input = mul2(&exp(&x), &dirac_delta(&argument));
        let expected = div(&mul2(&exp(&Expr::one()), &heaviside(&argument)), &Expr::int(-2));
        assert_eq!(with_input(&input, || antiderivative(&input, &x)), Some(expected));
    }

    #[test]
    fn singular_and_discontinuous_weights_decline_the_new_rule() {
        let x = Expr::symbol("x");
        let delta = dirac_delta(&x);
        for weight in [div(&Expr::one(), &x), heaviside(&x), delta.clone()] {
            let input = mul2(&weight, &delta);
            assert!(with_input(&input, || antiderivative(&input, &x)).is_none());
        }
    }

    #[test]
    fn inexact_scope_is_restored_after_nested_calculation() {
        let x = Expr::symbol("x");
        let exact = dirac_delta(&add2(&x, &Expr::one()));
        let inexact = mul2(&Expr::float(1.5), &exact);
        with_input(&inexact, || {
            assert!(with_input(&exact, || antiderivative(&exact, &x)).is_none());
        });
        assert!(with_input(&exact, || antiderivative(&exact, &x)).is_some());
    }
}
