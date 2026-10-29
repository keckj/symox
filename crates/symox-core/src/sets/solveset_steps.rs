//! Result checks of the reference `solveset` (SymPy 1.14,
//! `sympy/solvers/solveset.py` and `sympy/solvers/inequalities.py`):
//!
//! - `_solveset(..., _check=True)`: the finite solutions of the top-level
//!   call are kept only when `domain_check` of the original equation
//!   accepts them (no subexpression becomes infinite). For `Eq(lhs, rhs)`
//!   both sides are checked, so a denominator that cancels in `lhs - rhs`
//!   still excludes its zeros.
//! - `f.expand().is_zero`: an equation that vanishes identically is solved
//!   by the whole domain.
//! - `solve_univariate_inequality`: the singularities of the denominators
//!   of the original relational (`denoms(expr, gen)`, both sides) are
//!   excluded from the solution.

use super::*;
use crate::subs::subs1;
use crate::traverse::has;

fn is_rootof(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Undefined(head), args) if crate::poly::rootof::is_rootof_head(head, args.len()))
}

fn proved_infinite(e: &Expr) -> bool {
    crate::assumptions::is_finite(e) == Some(false)
}

/// The reference `_domain_check`.
fn subexpressions_finite(f: &Expr, x: &Expr, p: &Expr) -> bool {
    if f.is_atom() && crate::assumptions::is_finite(f) == Some(true) {
        return true;
    }
    // a relational has no `is_infinite`: only its sides are checked
    if !matches!(f.kind(), Kind::Rel(..)) && proved_infinite(&if has(f, x) { subs1(f, x, p) } else { f.clone() }) {
        return false;
    }
    if let Kind::Piecewise(pieces) = f.kind() {
        for (expr, condition) in pieces {
            match subs1(condition, x, p).kind() {
                Kind::Bool(false) => continue,
                Kind::Bool(true) => return subexpressions_finite(expr, x, p),
                // the reference gives an undecided case the benefit of the
                // doubt
                _ => return true,
            }
        }
        return false;
    }
    f.args().iter().all(|a| subexpressions_finite(a, x, p))
}

/// SymPy's `domain_check(f, symbol, p)`.
pub(super) fn domain_check(f: &Expr, x: &Expr, p: &Expr) -> bool {
    if proved_infinite(p) {
        return false;
    }
    subexpressions_finite(f, x, p)
}

/// The symbol-dependent core of the original equation (`fx` of the
/// reference `_check`): the dependent terms of a sum, then the dependent
/// factors of a product. A relational is checked as it is.
fn dependent_core(f: &Expr, x: &Expr) -> Expr {
    if matches!(f.kind(), Kind::Rel(..)) {
        return f.clone();
    }
    let f = match f.kind() {
        Kind::Add(..) => crate::ops::add(&f.args().into_iter().filter(|t| has(t, x)).collect::<Vec<_>>()),
        _ => f.clone(),
    };
    match f.kind() {
        Kind::Mul(..) => crate::ops::mul(&f.args().into_iter().filter(|t| has(t, x)).collect::<Vec<_>>()),
        _ => f,
    }
}

fn check_finite_parts(result: &Expr, fx: &Expr, x: &Expr) -> Expr {
    if let Some(elements) = finite_elements(result) {
        let kept: Vec<Expr> = elements.iter().filter(|s| is_rootof(s) || domain_check(fx, x, s)).cloned().collect();
        if kept.len() == elements.len() {
            return result.clone();
        }
        return finite_set(kept);
    }
    match as_set(result) {
        Some((SetOp::Union, parts)) => {
            let checked: Vec<Expr> = parts.iter().map(|s| check_finite_parts(s, fx, x)).collect();
            if checked.as_slice() == parts { result.clone() } else { union(checked) }
        }
        Some((SetOp::Intersection, parts)) if parts.first().is_some_and(|s| finite_elements(s).is_some()) => {
            let first = check_finite_parts(&parts[0], fx, x);
            if first == parts[0] { return result.clone(); }
            let mut all = vec![first];
            all.extend(parts[1..].iter().cloned());
            intersection(all)
        }
        Some((SetOp::Complement, parts)) if parts.len() == 2 && finite_elements(&parts[0]).is_some() => {
            let first = check_finite_parts(&parts[0], fx, x);
            if first == parts[0] { result.clone() } else { complement(&first, &parts[1]) }
        }
        _ => result.clone(),
    }
}

/// The `_check` step of the top-level `solveset` call: finite solutions
/// where a subexpression of the original equation `f` is infinite are
/// removed. Undecided candidates are kept, as in the reference; a
/// ConditionSet is returned unchanged.
pub(super) fn check_solutions(f: &Expr, x: &Expr, result: Expr) -> Expr {
    if is_op(&result, SetOp::ConditionSet) || !has(f, x) {
        return result;
    }
    let fx = dependent_core(f, x);
    check_finite_parts(&result, &fx, x)
}

/// `f.expand().is_zero` at the start of the reference `_solveset`.
pub(super) fn vanishes_identically(f: &Expr) -> bool {
    !matches!(f.kind(), Kind::Rel(..) | Kind::Bool(_)) && crate::expand::expand(f).is_zero()
}

/// The denominators of every subexpression of `e` that depend on `x`
/// (the reference `denoms(expr, gen)`), split into factors.
fn denominators(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Rel(..) | Kind::Bool(_) | Kind::Func(Func::Logic(_) | Func::Set(_), _)) || !has(s, x) {
            return;
        }
        let (_, d) = crate::simptools::radsimp::fraction(s, false);
        if d.is_one() {
            return;
        }
        for factor in crate::simptools::mul_args(&d) {
            if has(&factor, x) && !out.contains(&factor) {
                out.push(factor);
            }
        }
    });
}

fn periodic_in(e: &Expr, x: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc, a) = s.kind() {
            if has(&a[0], x) {
                found = true;
            }
        }
    });
    found
}

/// The real singularities of the denominators of both sides of the
/// original relational, which the reference `solve_univariate_inequality`
/// excludes from the solution (`x + 1/x > -2 + 1/x` excludes `x = 0`).
/// `roots` solves `d == 0` over the reals. A denominator whose real zeros
/// are not a finite set of numbers makes the inequality unsolved (`Err`),
/// as the reference `solvify` raises there; periodic denominators are left
/// to the periodic solver.
pub(super) fn relational_singularities(
    ineq: &Expr,
    x: &Expr,
    roots: &dyn Fn(&Expr) -> Option<Expr>,
) -> Result<Vec<Expr>, ()> {
    let Kind::Rel(_, lhs, rhs) = ineq.kind() else { return Ok(Vec::new()) };
    let mut dens = Vec::new();
    denominators(lhs, x, &mut dens);
    denominators(rhs, x, &mut dens);
    let mut points: Vec<Expr> = Vec::new();
    for d in dens {
        if periodic_in(&d, x) {
            continue;
        }
        let zeros = roots(&d).ok_or(())?;
        if is_empty_set(&zeros) {
            continue;
        }
        let elements = finite_elements(&zeros).ok_or(())?;
        if elements.iter().any(|v| !crate::logic::is_number_like(v)) {
            return Err(());
        }
        for v in elements {
            if !points.contains(v) {
                points.push(v.clone());
            }
        }
    }
    Ok(points)
}
