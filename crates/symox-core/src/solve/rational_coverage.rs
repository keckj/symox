//! Exact coverage certificates for already returned rational roots.
//!
//! This does not solve an equation or filter candidates. It proves that a
//! nonconstant QQ polynomial has exactly the supplied distinct rational roots,
//! so solve's later high-order-root classification need not solve it again.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::expr::{Expr, Kind};
use crate::poly::{convert, dense, Rat};

// This is an optional postprocessing optimization, not a solving limit.
// Declining its dense proof keeps the existing complete solver route intact.
const MAX_DENSE_CERTIFICATE_DEGREE: i64 = 4096;

/// Conservative monomial-exponent bounds before any polynomial expansion.
/// Besides bounding optional work, checked arithmetic prevents the converter's
/// u32 monomial representation from truncating an enormous source exponent.
fn qq_span(expression: &Expr, variable: &Expr,
    memo: &mut FxHashMap<usize, (Expr, Option<(u32, u32)>)>) -> Option<(u32, u32)> {
    if let Some((_, result)) = memo.get(&expression.dag_id()) { return *result; }
    let result = (|| {
        if expression == variable { return Some((1, 1)); }
        match expression.kind() {
            Kind::Num(number) if number.is_rational() => Some((0, 0)),
            Kind::Add(constant, terms) if constant.is_rational() => {
                let mut span = if constant.is_zero() { None } else { Some((0, 0)) };
                for (term, coefficient) in terms {
                    if !coefficient.is_rational() { return None; }
                    let (low, high) = qq_span(term, variable, memo)?;
                    span = Some(match span {
                        Some((a, b)) => (low.min(a), high.max(b)),
                        None => (low, high),
                    });
                }
                Some(span.unwrap_or((0, 0)))
            }
            Kind::Mul(coefficient, factors) if coefficient.is_rational() => {
                let (mut low, mut high) = (0u32, 0u32);
                for (base, exponent) in factors {
                    let power = u32::try_from(exponent.as_i64()?).ok()?;
                    if power != 0 {
                        let (a, b) = qq_span(base, variable, memo)?;
                        low = low.checked_add(a.checked_mul(power)?)?;
                        high = high.checked_add(b.checked_mul(power)?)?;
                    }
                }
                Some((low, high))
            }
            Kind::Pow(base, exponent) => {
                let power = u32::try_from(exponent.as_i64()?).ok()?;
                if power == 0 { return Some((0, 0)); }
                let (low, high) = qq_span(base, variable, memo)?;
                Some((low.checked_mul(power)?, high.checked_mul(power)?))
            }
            _ => None,
        }
    })();
    memo.insert(expression.dag_id(), (expression.clone(), result));
    result
}

/// Whether `candidates` contain every root, and only roots, of a nonconstant
/// polynomial over QQ. False means no certificate, never absence of solutions.
/// Root checks, assumptions, options and return-shape handling remain with the
/// caller. In particular an assumption-filtered proper subset fails this proof.
pub fn complete_rational_roots(equation: &Expr, variable: &Expr, candidates: &[Expr]) -> bool {
    let Kind::Sym(symbol) = variable.kind() else { return false };
    // Host subclasses and adapted scalar atoms may customize solve_linear's
    // derivative/substitution protocol. This certificate only targets ordinary
    // commutative Symbols; decline before consulting any host assumption hook.
    if symbol.host().is_some() || symbol.scalar_identity().is_some()
        || symbol.is_dummy() || symbol.is_wild()
        || symbol.facts().get("commutative") != Some(true) || candidates.is_empty() {
        return false;
    }
    let mut seen = FxHashSet::default();
    let mut roots = Vec::new();
    for candidate in candidates {
        let Some(number) = candidate.as_number().filter(|number| number.is_rational()) else { return false };
        if seen.insert(candidate.clone()) {
            let Some(rational) = Rat::from_number(number) else { return false };
            roots.push(rational);
        }
    }
    let source = super::as_expr(equation);
    if !crate::decompose::commutes(&source) { return false; }
    let mut pending = vec![source];
    let mut factors_seen = FxHashSet::default();
    let mut spans = FxHashMap::default();
    let mut used = vec![false; roots.len()];
    let mut nonconstant = false;
    while let Some(factor) = pending.pop() {
        if !factors_seen.insert(factor.clone()) { continue; }
        match factor.kind() {
            Kind::Mul(..) => {
                pending.extend(factor.args());
                continue;
            }
            Kind::Pow(base, exponent) if exponent.as_number().is_some_and(|number|
                number.is_integer() && number.is_positive()) => {
                // An exact positive integer power repeats every multiplicity
                // of its base, without changing the root set. Keep it factored
                // so (x-r)**1000000 does not require a million-term expansion.
                pending.push(base.clone());
                continue;
            }
            _ => {}
        }
        let Some((valuation_bound, degree_bound)) = qq_span(&factor, variable, &mut spans) else { return false };
        if i64::from(degree_bound - valuation_bound) > MAX_DENSE_CERTIFICATE_DEGREE { return false; }
        let Ok(mut converted) = convert::from_exprs_gens(std::slice::from_ref(&factor), std::slice::from_ref(variable)) else { return false };
        if converted.float || converted.gens != [variable.clone()] { return false; }
        let Some(mut polynomial) = converted.polys.pop() else { return false };
        if polynomial.is_zero() { return false; }
        if polynomial.is_constant() { continue; }
        nonconstant = true;

        // Removing x's valuation is all its repeated exact divisions at once.
        // It keeps sparse monomials sparse even when their degree is enormous.
        let valuation = polynomial.iter().map(|(exponent, _)| exponent[0]).min().unwrap_or(0);
        if valuation != 0 {
            let Some(index) = roots.iter().position(Rat::is_zero) else { return false };
            used[index] = true;
            polynomial = polynomial.div_monomial(&[valuation]);
        }
        if polynomial.is_constant() { continue; }
        if roots.iter().all(Rat::is_zero) { return false; }
        if polynomial.degree(0) > MAX_DENSE_CERTIFICATE_DEGREE { return false; }

        let mut remaining = dense::from_mpoly(&polynomial, 0);
        for (index, root) in roots.iter().enumerate() {
            if root.is_zero() { continue; }
            let divisor = [root.neg(), Rat::one()];
            while remaining.len() > 1 {
                let (quotient, remainder) = dense::divrem(&remaining, &divisor);
                if !remainder.is_empty() { break; }
                used[index] = true;
                remaining = quotient;
            }
        }
        // Every linear factor was removed with its full multiplicity. A
        // residual positive degree could hide other roots, rational or not.
        if remaining.len() != 1 || remaining[0].is_zero() { return false; }
    }
    nonconstant && used.into_iter().all(|present| present)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(source: &str) -> Expr { parse_str(source).unwrap() }
    fn certify(source: &str, roots: &[&str]) -> bool {
        complete_rational_roots(&p(source), &Expr::symbol("x"), &roots.iter().map(|root| p(root)).collect::<Vec<_>>())
    }

    #[test]
    fn repeated_exact_divisions_cover_all_multiplicities() {
        for source in ["(2*x-1)^3*(3*x+2)^2*(x-4)^4", "-7*(x-1/2)*(x+2/3)*(x-4)/3"] {
            assert!(certify(source, &["1/2", "-2/3", "4"]));
        }
        assert!(certify("x^5*(x-1)^3", &["0", "1"]));
        assert!(certify("(x-1)*(x-2)", &["1", "1", "2"]));
        assert!(certify("(x-1)^1000000", &["1"]));
        assert!(certify("x^1000000", &["0"]));
        assert!(certify("x^100000001-x^100000000", &["0", "1"]));
        assert!(certify("x^3-6*x^2+11*x-6", &["1", "2", "3"]));
        assert!(certify("x^3-3*x^2+3*x-1", &["1"]));
        let expanded = crate::expand::expand(&p("(2*x-1)^3*(3*x+2)^2"));
        assert!(complete_rational_roots(&expanded, &Expr::symbol("x"), &[p("1/2"), p("-2/3")]));
    }

    #[test]
    fn incomplete_extraneous_inexact_and_nonqq_inputs_do_not_certify() {
        for (source, roots) in [
            ("(x-1)*(x-2)", vec!["1"]),
            ("(x-1)*(x-2)", vec!["1", "2", "3"]),
            ("(x-1)*(x^2-2)", vec!["1"]),
            ("x^2-2", vec!["sqrt(2)", "-sqrt(2)"]),
            ("a*(x-1)", vec!["1"]),
            ("sin(1)*(x-1)", vec!["1"]),
            ("1.0*x-1", vec!["1"]),
            ("x-1", vec!["1.0"]),
            ("(x^2-1)/(x-1)", vec!["-1"]),
            ("x^3", vec!["1"]),
            ("x^4097-1", vec!["1"]),
            ("x^4294967297-1", vec!["1"]),
            ("0", vec!["1"]),
            ("1", vec!["1"]),
        ] {
            assert!(!certify(source, &roots), "{source} with {roots:?}");
        }
        assert!(!certify("x-1", &[]));
    }

    #[test]
    fn special_symbol_identities_and_host_protocols_decline_without_callbacks() {
        use crate::symbol::{Assumptions, Symbol};
        use crate::assumption_facts::FactSet;
        use std::sync::Arc;

        #[derive(Debug)]
        struct Hooked;
        impl crate::host::HostClass for Hooked {
            fn identity(&self) -> u64 { 0x51515 }
            fn class_name(&self) -> &str { "CoverageHookedSymbol" }
            fn as_any(&self) -> &dyn std::any::Any { self }
            fn fact(&self, _value: &Expr, _name: &str) -> Option<Option<bool>> {
                panic!("the certificate must decline before consulting a host hook")
            }
        }
        let noncommutative = Symbol::with_facts("x", FactSet::from_kwargs([("commutative", false)]).unwrap()).unwrap();
        for symbol in [
            Symbol::new("x").with_host(Arc::new(Hooked)),
            Symbol::new("x").with_scalar_identity(Expr::int(7)),
            noncommutative,
            Symbol::dummy("x", Assumptions::NONE),
            Symbol::wild("x", vec![]),
        ] {
            let variable = Expr::from_kind(Kind::Sym(symbol));
            assert!(!complete_rational_roots(&variable, &variable, &[Expr::zero()]));
        }
    }
}
