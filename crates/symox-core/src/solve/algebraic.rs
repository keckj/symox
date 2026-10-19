//! Complete exact roots over QQ when radical formulae do not suffice.
use crate::assumption_facts::FactSet;
use crate::expr::{Expr, Func, Kind, UndefinedFunction};
use crate::functions::func;
use crate::ops::div;
use crate::poly::{dense, fraction::Polynomial, int::Rat, polyobj::Poly, real_roots};

/// Factoring separates repeated factors before indices are assigned. This
/// returns every distinct root; multiplicities belong to roots(), not solve().
pub(super) fn roots(f: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    if crate::evalf_mp::max_precision(f).is_some() {
        return None;
    }
    let p = Poly::from_expr(f, std::slice::from_ref(x)).ok()?;
    if p.all_coeffs().ok()?.iter().any(|c| !c.is_rational_number()) {
        return None;
    }
    let (_, factors) = p.factor_list();
    let mut roots = Vec::new();
    for (factor, _) in factors {
        let expression = factor.as_expr();
        let coefficients = super::poly_coeffs(&expression, x)?;
        if let Some(radicals) = super::poly_roots(&coefficients)
            .or_else(|| super::complete_roots(&expression, x, &coefficients))
        {
            roots.extend(radicals);
        } else {
            roots.extend(indexed_roots(&factor)?);
        }
    }
    Some(roots)
}

fn indexed_roots(p: &Poly) -> Option<Vec<Expr>> {
    let coefficients = p.all_coeffs().ok()?;
    let leading = coefficients.first()?;
    let monic = crate::special::tuple(coefficients.iter().map(|c| div(c, leading)).collect());
    let polynomial: Polynomial = p.clone().into();
    let intervals = real_roots::intervals(&polynomial, &Default::default()).ok()?;
    let mut roots = Vec::new();
    for interval in &intervals {
        if interval.lower == interval.upper {
            roots.push(Expr::num(interval.lower.to_number()));
            continue;
        }
        let mut facts = vec![("algebraic", true), ("real", true), ("rational", false)];
        if !interval.lower.is_negative() {
            facts.push(("positive", true));
        }
        if interval.upper.is_negative() || interval.upper.is_zero() {
            facts.push(("negative", true));
        }
        roots.push(node(monic.clone(), roots.len(), facts)?);
    }
    while roots.len() < coefficients.len() - 1 {
        roots.push(node(
            monic.clone(),
            roots.len(),
            vec![("algebraic", true), ("real", false), ("rational", false)],
        )?);
    }
    Some(roots)
}

fn node(coefficients: Expr, index: usize, facts: Vec<(&str, bool)>) -> Option<Expr> {
    let head =
        UndefinedFunction::with_facts("ComplexRootOf", FactSet::from_kwargs(facts).ok()?, None)
            .ok()?;
    Some(func(
        Func::Undefined(head),
        vec![coefficients, Expr::int(index as i64)],
    ))
}

/// Prove a polynomial residual by reduction modulo the root's defining
/// polynomial. A nonconstant remainder is undecidable for a reducible defining
/// polynomial; it must never be interpreted as proof that this root fails.
pub(super) fn polynomial_zero(f: &Expr, x: &Expr, root: &Expr) -> Option<bool> {
    let Kind::Func(Func::Undefined(head), args) = root.kind() else {
        return None;
    };
    if head.name.as_ref() != "ComplexRootOf" || args.len() != 2 {
        return None;
    }
    let rational = |e: Expr| e.as_number().and_then(Rat::from_number);
    let coefficients: Vec<Rat> = super::poly_coeffs(f, x)?
        .into_iter()
        .map(rational)
        .collect::<Option<_>>()?;
    let defining: Vec<Rat> = args[0]
        .args()
        .into_iter()
        .rev()
        .map(rational)
        .collect::<Option<_>>()?;
    if defining.len() < 2 || defining.last()?.is_zero() {
        return None;
    }
    let (_, remainder) = dense::divrem(&coefficients, &defining);
    if remainder.is_empty() {
        Some(true)
    } else if remainder.len() == 1 {
        Some(false)
    } else {
        None
    }
}
