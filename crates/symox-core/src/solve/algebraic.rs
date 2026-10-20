//! Complete exact roots over QQ when radical formulae do not suffice.
use crate::expr::Expr;
use crate::poly::{dense, int::{Int, Rat}, polyobj::Poly, rootof};

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

/// The reference's `Poly.all_roots` of an irreducible factor whose roots
/// are not given by radicals: root objects of the primitive integer
/// polynomial in the generator `x`, scaled by its integer basis.
fn indexed_roots(p: &Poly) -> Option<Vec<Expr>> {
    let coefficients: Vec<Rat> = p
        .all_coeffs()
        .ok()?
        .iter()
        .map(|c| c.as_number().and_then(Rat::from_number))
        .collect::<Option<_>>()?;
    let mut integers = crate::poly::rootisolation::clear_denoms(&coefficients);
    let mut content = Int::zero();
    for c in &integers {
        content = content.gcd(c);
    }
    if content.is_zero() {
        return None;
    }
    if integers[0].is_negative() {
        content = content.neg();
    }
    for c in &mut integers {
        *c = c.div_exact_unchecked(&content);
    }
    rootof::all_roots_irreducible(&integers)
}

/// Prove a polynomial residual by reduction modulo the root's defining
/// polynomial. A nonconstant remainder is undecidable for a reducible defining
/// polynomial; it must never be interpreted as proof that this root fails.
pub(super) fn polynomial_zero(f: &Expr, x: &Expr, root: &Expr) -> Option<bool> {
    if !rootof::is_rootof(root) {
        return None;
    }
    let args = root.args();
    let rational = |e: Expr| e.as_number().and_then(Rat::from_number);
    let coefficients: Vec<Rat> = super::poly_coeffs(f, x)?
        .into_iter()
        .map(rational)
        .collect::<Option<_>>()?;
    let (_, defining) = rootof::integer_poly(&args[0])?;
    let mut defining: Vec<Rat> = defining.into_iter().map(|c| Rat::new(c, Int::one())).collect();
    defining.reverse();
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
