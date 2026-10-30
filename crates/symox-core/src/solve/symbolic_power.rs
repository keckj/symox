//! The single-generator branch of SymPy 1.14's `_solve` for a generator
//! `x**e` whose exponent `e` is symbolic and free of `x`: the polynomial in
//! `g = x**e` is solved for `t = g` and every root `r` (the negative and
//! complex ones included) is inverted with `_vsolve(g - t, x) = [t**(1/e)]`.
//! `x**(2*y) - 3` is `t**2 - 3` with `t = x**y`, which gives
//! `[3**(1/(2*y)), (-sqrt(3))**(1/y)]`.

use crate::expr::{Expr, Kind};
use crate::ops::{div, pow};
use crate::subs::subs1;
use crate::traverse::{free_symbols, has};

/// None when `f` is not a polynomial of degree >= 2 in one generator
/// `x**e` with a symbolic exponent (the earlier paths apply then).
pub(super) fn solve(f: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    let mut candidate = false;
    crate::traverse::preorder(f, &mut |p| {
        if let Kind::Pow(b, e) = p.kind() {
            candidate |= b == x && !has(e, x) && !free_symbols(e).is_empty();
        }
    });
    if !candidate {
        return None;
    }
    let p = crate::poly::polyobj::Poly::from_expr(f, &[]).ok()?;
    let gens: Vec<&Expr> = p.gens.iter().filter(|g| has(g, x)).collect();
    let [g] = gens.as_slice() else { return None };
    let Kind::Pow(base, e) = g.kind() else { return None };
    if base != x || has(e, x) || free_symbols(e).is_empty() {
        return None;
    }
    let t = Expr::dummy("t", crate::Assumptions::NONE);
    let ft = subs1(f, g, &t);
    if has(&ft, x) {
        return None;
    }
    let pt = crate::poly::polyobj::Poly::from_expr(&ft, std::slice::from_ref(&t)).ok()?;
    if pt.degree(0).unwrap_or(0) < 2 {
        return None;
    }
    let roots = super::solve_rec(&ft, &t, depth + 1)?;
    let inv = div(&Expr::one(), e);
    let mut out: Vec<Expr> = Vec::new();
    for r in roots {
        let v = pow(&r, &inv);
        if !out.contains(&v) {
            out.push(v);
        }
    }
    Some(out)
}
