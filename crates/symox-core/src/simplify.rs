//! Simplification routines.
//!
//! This first version contains the building blocks other modules rely on:
//! `numer_denom`, `radsimp` (denominator rationalisation), `together`
//! and a light `simplify`. The polynomial based routines (`cancel`,
//! `factor`, ...) live in [`crate::poly`].

use crate::expr::{Expr, Func, Kind};
use rustc_hash::FxHashMap;
use crate::expand::expand;
use crate::ops::{add, div, mul, mul_factor, neg, pow};

/// Split `e` into `(numerator, denominator)` without combining sums
/// (SymPy's `as_numer_denom` for products and powers; sums are combined
/// over a common denominator).
pub fn numer_denom(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Num(n) => (Expr::num(n.numer()), Expr::num(n.denom())),
        Kind::Pow(b, x) => {
            if let Some(n) = x.as_number() {
                if n.is_negative() {
                    return (Expr::one(), pow(b, &Expr::num(n.neg())));
                }
            }
            if let Kind::Mul(c, _) = x.kind() {
                if c.is_negative() {
                    return (Expr::one(), pow(b, &neg(x)));
                }
            }
            (e.clone(), Expr::one())
        }
        Kind::Mul(c, facs) => {
            let mut num = vec![Expr::num(c.numer())];
            let mut den = vec![Expr::num(c.denom())];
            for (b, x) in facs {
                let (n, d) = numer_denom(&mul_factor(b, x));
                num.push(n);
                den.push(d);
            }
            (mul(&num), mul(&den))
        }
        Kind::Add(..) => {
            let parts: Vec<(Expr, Expr)> = e.args().iter().map(numer_denom).collect();
            if parts.iter().all(|(_, d)| d.is_one()) {
                return (e.clone(), Expr::one());
            }
            // common denominator: product of distinct denominators
            let mut dens: Vec<Expr> = Vec::new();
            for (_, d) in &parts {
                if !d.is_one() && !dens.contains(d) {
                    dens.push(d.clone());
                }
            }
            let den = mul(&dens);
            let mut nums = Vec::with_capacity(parts.len());
            for (n, d) in &parts {
                let mut f = vec![n.clone()];
                for dd in &dens {
                    if dd != d {
                        f.push(dd.clone());
                    }
                }
                nums.push(mul(&f));
            }
            (add(&nums), den)
        }
        _ => (e.clone(), Expr::one()),
    }
}

/// Combine a sum of fractions over a common denominator.
pub fn together(e: &Expr) -> Expr {
    let (n, d) = numer_denom(e);
    if d.is_one() {
        return n;
    }
    div(&n, &d)
}

fn has_sqrt(e: &Expr) -> bool {
    match e.kind() {
        Kind::Pow(b, x) => {
            (b.is_number() && *x == Expr::rational(1, 2)) || has_sqrt(b)
        }
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().any(has_sqrt),
        _ => false,
    }
}

/// Rationalise denominators of the form `a + b*sqrt(c)`.
pub fn radsimp(e: &Expr) -> Expr {
    let (n, d) = numer_denom(e);
    if d.is_one() || !has_sqrt(&d) {
        return e.clone();
    }
    let mut num = n;
    let mut den = d;
    for _ in 0..4 {
        if !has_sqrt(&den) {
            break;
        }
        let Kind::Add(c, terms) = den.kind() else {
            break;
        };
        // conjugate: flip the sign of the radical terms
        let mut conj_terms = Vec::new();
        let mut flipped = false;
        for (t, k) in terms {
            if !flipped && has_sqrt(t) {
                conj_terms.push(mul(&[Expr::num(k.neg()), t.clone()]));
                flipped = true;
            } else {
                conj_terms.push(mul(&[Expr::num(k.clone()), t.clone()]));
            }
        }
        conj_terms.push(Expr::num(c.clone()));
        let conj = add(&conj_terms);
        num = expand(&mul(&[num, conj.clone()]));
        den = expand(&mul(&[den, conj]));
    }
    if has_sqrt(&den) {
        return e.clone();
    }
    expand(&div(&num, &den))
}

/// SymPy's simplify candidates from the polynomial module: `cancel(e)`
/// and `factor_terms(cancel(e))` (common content and monomials pulled
/// out, a bare `number*(sum)` product being distributed again).
fn poly_candidates(e: &Expr) -> Vec<Expr> {
    let Ok(c) = crate::poly::api::cancel(e) else {
        return vec![];
    };
    let ft = crate::poly::api::together(&c);
    let ft = match ft.kind() {
        Kind::Mul(k, facs) if !k.is_one() && facs.len() == 1 && facs[0].0.is_add() && facs[0].1.is_one() => {
            mul(&[Expr::num(k.clone()), facs[0].0.clone()])
        }
        _ => ft,
    };
    vec![c, ft]
}

/// Release evaluation barriers before algebraic simplification. Expansion
/// alone deliberately keeps UnevaluatedExpr, matching SymPy's behavior.
fn release_unevaluated(e: &Expr, cache: &mut FxHashMap<Expr, Expr>) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    if let Some(cached) = cache.get(e) {
        return cached.clone();
    }
    let result = match e.kind() {
        Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") && args.len() == 1 => {
            release_unevaluated(&args[0], cache)
        }
        _ => crate::traverse::map_children(e, &mut |child| release_unevaluated(child, cache)),
    };
    cache.insert(e.clone(), result.clone());
    result
}

/// A light general purpose simplifier (to be extended).
pub fn simplify(e: &Expr) -> Expr {
    let evaluated = release_unevaluated(e, &mut FxHashMap::default());
    let e = &evaluated;
    let t = crate::simp::trigsimp(e);
    // Prefer the cancelled rational form when equally small, including
    // its deterministic positive leading denominator coefficient.
    let mut candidates = poly_candidates(e);
    candidates.extend([
        e.clone(),
        expand(e),
        together(e),
        crate::simp::powsimp(e),
        t.clone(),
        crate::simp::trigsimp(&expand(e)),
    ]);
    let selected = candidates
        .into_iter()
        .min_by_key(|c| crate::traverse::node_count(c))
        .unwrap();
    // Sign normalization can add a numeric -1 node while keeping the same
    // arithmetic cost. Apply it after choosing candidates, so that size
    // comparisons cannot restore a negative-looking denominator.
    let (numerator, denominator) = numer_denom(&selected);
    if denominator.is_add() {
        // Factor arrays use the arena's fast internal order. Polynomial
        // conversion supplies a deterministic order of symbols and leading
        // monomial, including names with Unicode subscripts.
        let negative = crate::poly::convert::from_exprs(std::slice::from_ref(&denominator))
            .ok()
            .map(|c| c.polys[0].lc().is_negative())
            .unwrap_or_else(|| crate::ops::could_extract_minus_sign(&denominator));
        if negative {
            return div(&neg(&numerator), &neg(&denominator));
        }
    }
    selected
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn barycentric_ratio_denominator_sign() {
        let expr = crate::parser::parse_str("e11/(e11*e20-e10*e21)").unwrap();
        let expected = crate::parser::parse_str("-e11/(e10*e21-e11*e20)").unwrap();
        assert_eq!(simplify(&expr), expected);
        let a = Expr::symbol("e₀₀");
        let b = Expr::symbol("e₀₁");
        let c = Expr::symbol("e₂₀");
        let d = Expr::symbol("e₂₁");
        let determinant = crate::ops::sub(&mul(&[a, d]), &mul(&[b.clone(), c]));
        let ratio = div(&neg(&b), &neg(&determinant));
        assert_eq!(simplify(&ratio), div(&b, &determinant));
    }

    #[test]
    fn simplify_releases_unevaluated_expr() {
        let x = Expr::symbol("x");
        let inner = add(&[x.clone(), Expr::one()]);
        let held = crate::functions::undefined("UnevaluatedExpr", vec![inner.clone()]);
        assert_eq!(expand(&held), held);
        assert_eq!(simplify(&crate::ops::sub(&held, &inner)), Expr::zero());
        assert_eq!(crate::diff::diff(&held, &x), crate::diff::derivative_unevaluated(&held, vec![(x, 1)]));
    }
}
