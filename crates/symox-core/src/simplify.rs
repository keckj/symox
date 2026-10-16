//! Simplification routines.
//!
//! This first version contains the building blocks other modules rely on:
//! `numer_denom`, `radsimp` (denominator rationalisation), `together`
//! and a light `simplify`. The polynomial based routines (`cancel`,
//! `factor`, ...) live in [`crate::poly`].

use crate::expr::{Expr, Kind};
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

/// A light general purpose simplifier (to be extended).
pub fn simplify(e: &Expr) -> Expr {
    let t = crate::simp::trigsimp(e);
    let mut candidates = vec![
        e.clone(),
        expand(e),
        together(e),
        crate::simp::powsimp(e),
        t.clone(),
        crate::simp::trigsimp(&expand(e)),
    ];
    candidates.extend(poly_candidates(e));
    candidates
        .into_iter()
        .min_by_key(|c| crate::traverse::node_count(c))
        .unwrap()
}
