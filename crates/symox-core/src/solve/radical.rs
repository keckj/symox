//! The branches of SymPy 1.14's `_solve` (`sympy/solvers/solvers.py`) that
//! an algebraic equation with radicals of the unknown takes, in the
//! reference order, before the `_tsolve` steps:
//!
//! 1. the generators of `Poly(f_num)` that contain the unknown. With one
//!    generator `g` (`sqrt(y)` for `y**(3/2) - x`), `roots` of the
//!    polynomial in `g` and the inversion `_vsolve(g - t, symbol)` of every
//!    root, which the reference makes on its `tsolve` pass and thus through
//!    `unrad(g - t)`: `y**(3/2) = x` has the three solutions
//!    `x**(2/3)` and `(x**(1/3)*w)**2` for both complex cube roots `w` of 1.
//!    If the inversion fails, the branch below applies.
//! 2. several generators (`y` and `sqrt(y**3)` for
//!    `-2*y - (C1 + x)*sqrt(y**3)`): the fallback of `_solve` through
//!    `unrad(f_num, symbol)` and the solution of the radical-free equation,
//!    with its change of variable when there is one
//!    (`inv = _vsolve(ieq, symbol)[0]`, `inv.subs(isym, xi)` for every root
//!    `xi` of `_solve(eq, isym)`).
//!
//! Admission: every subexpression that contains the unknown is a sum, a
//! product, the unknown itself or a power with a Rational exponent, and one
//! exponent is not an integer. The equations mixing radicals with functions
//! of the unknown keep the other native paths (see the parity contract).
//!
//! `None` means "not taken here": when `unrad` raises or returns None, or
//! the radical-free equation cannot be solved, the caller goes on to the
//! `_tsolve` steps exactly like the reference; it never means "no solution".
//! The candidates are checked by the caller (`solve_checked`), as the
//! reference `_solve` checks them after its fallback.
//!
//! The reference leaves the unrad solutions (and the solutions of the
//! generator branch when there are more than two) unsimplified unless the
//! `simplify` flag is given; they are recorded for the Python wrapper
//! (`take_unsimplified`), which simplifies the other solutions.

use std::cell::{Cell, RefCell};

use crate::expr::{Expr, Kind};
use crate::ops::{add, mul, pow, sub};
use crate::poly::polyobj::Poly;
use crate::subs::subs1;
use crate::traverse::has;

use super::unrad::{self, UFail, UnradHost};

thread_local! {
    static UNSIMPLIFIED: RefCell<Vec<Expr>> = const { RefCell::new(Vec::new()) };
    static UNRAD_DISABLED: Cell<bool> = const { Cell::new(false) };
    static BRANCH_OFF: Cell<bool> = const { Cell::new(false) };
}

/// The solutions recorded since the last call: the ones the reference
/// leaves unsimplified (unrad results, more than two generator roots).
pub fn take_unsimplified() -> Vec<Expr> {
    UNSIMPLIFIED.with(|u| std::mem::take(&mut *u.borrow_mut()))
}

/// Run `f` with the reference `_unrad=False` flag: the outermost `_solve`
/// (and the factors of a product) skip the unrad fallback; the nested
/// calls keep it (the reference pops the flag).
pub fn without_unrad<T>(disabled: bool, f: impl FnOnce() -> T) -> T {
    let previous = UNRAD_DISABLED.with(|c| c.replace(disabled));
    let r = f();
    UNRAD_DISABLED.with(|c| c.set(previous));
    r
}

/// Run `f` without these branches: the native `solveset` takes the
/// candidates of `solve` as a finite set, whereas the reference
/// `_solveset` reaches `_solve_radical` for these equations, which keeps an
/// undecided candidate in a ConditionSet (`radical_step`).
pub fn without_radical_branch<T>(f: impl FnOnce() -> T) -> T {
    let previous = BRANCH_OFF.with(|c| c.replace(true));
    let r = f();
    BRANCH_OFF.with(|c| c.set(previous));
    r
}

pub(super) fn record(values: &[Expr]) {
    UNSIMPLIFIED.with(|u| {
        let mut u = u.borrow_mut();
        for v in values {
            if !u.contains(v) {
                u.push(v.clone());
            }
        }
    });
}

/// Every subexpression containing `x` is `x`, a sum, a product or a power
/// with a Rational exponent; `radical` is set when an exponent is not an
/// integer.
fn algebraic(e: &Expr, x: &Expr, radical: &mut bool) -> bool {
    if e == x || !has(e, x) {
        return true;
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().all(|a| algebraic(a, x, radical)),
        Kind::Pow(b, n) => {
            let Some(n) = n.as_number() else { return false };
            if !n.is_rational() {
                return false;
            }
            if !n.is_integer() {
                *radical = true;
            }
            algebraic(b, x, radical)
        }
        _ => false,
    }
}

/// `_vsolve(e, s, check=False, simplify=False)` for the native `unrad`.
struct Host {
    depth: u32,
}

impl UnradHost for Host {
    type Error = ();
    fn vsolve(&mut self, eq: &Expr, x: &Expr) -> Result<Vec<Expr>, UFail<()>> {
        nested_checked(eq, x, self.depth + 1)
            .ok_or_else(|| UFail::NotImplemented(String::new()))
    }
    fn take(&mut self, _e: &Expr) -> Option<Result<bool, UFail<()>>> {
        None
    }
}

/// The two branches above for the numerator `num` of the equation.
pub(super) fn solve(num: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    if BRANCH_OFF.with(Cell::get) {
        return None;
    }
    let mut radical = false;
    if !algebraic(num, x, &mut radical) || !radical || x.as_sym().is_none() {
        return None;
    }
    let poly = Poly::from_expr(num, &[]).ok()?;
    let gens: Vec<Expr> = poly.true_gens().iter().filter(|g| has(g, x)).cloned().collect();
    if let [g] = gens.as_slice() {
        if g != x {
            if let Some(found) = single_generator(num, g, x, depth) {
                return Some(found);
            }
        }
    }
    // several generators (the radicals of an algebraic equation always
    // have different bases or a denominator q > 1 in their exponents), or
    // a generator that could not be inverted
    if depth == 0 && UNRAD_DISABLED.with(Cell::get) {
        return None;
    }
    unrad_branch(num, x, depth)
}

/// The unrad fallback of `_solve` for radicals of the unknown mixed with
/// other functions of it (`sqrt(x)*exp(x) - 1`, `cosh(x)/x - sqrt(2)`):
/// the reference reaches it when the generator branch found nothing, and
/// before `_tsolve`. Taken inside the native `_tsolve` scope only
/// (`native_tsolve::scope`); None when `unrad` declines or the
/// radical-free equation is not solved, so that the `_tsolve` steps follow.
pub(super) fn solve_mixed(num: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    if !super::native_tsolve::enabled() || BRANCH_OFF.with(Cell::get) || x.as_sym().is_none() {
        return None;
    }
    if depth == 0 && UNRAD_DISABLED.with(Cell::get) {
        return None;
    }
    let mut radical = false;
    // the algebraic equations take `solve` above
    if algebraic(num, x, &mut radical) || !super::native_tsolve::has_radical(num, x) {
        return None;
    }
    unrad_branch(num, x, depth)
}

/// A nested `_vsolve(eq, symbol, check=False)`: the native `_solve`, then,
/// inside the native `_tsolve` scope, its `_tsolve` steps.
fn nested_checked(eq: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    super::solve_checked_at(eq, x, false, depth).or_else(|| {
        let options = super::native_tsolve::SolveOptions { check: false, ..Default::default() };
        super::native_tsolve::fallback(eq, x, options)
    })
}

/// `roots(Poly(f_num, g, composite=True))` and the inversion of `g`.
fn single_generator(num: &Expr, g: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    let poly = Poly::from_expr(num, std::slice::from_ref(g)).ok()?;
    if poly.ngens != 1 {
        return None;
    }
    let t = Expr::dummy("t", crate::Assumptions::NONE);
    let terms: Vec<Expr> = poly.terms().into_iter()
        .map(|(m, c)| mul(&[c, pow(&t, &Expr::int(i64::from(m[0])))]))
        .collect();
    let ft = add(&terms);
    let roots = super::solve_rec(&ft, &t, depth + 1)?;
    // `_vsolve(g - t, symbol)` on the reference `tsolve` pass: unrad, then
    // the `_tsolve` inversion.
    let s = Expr::dummy("t", crate::Assumptions::NONE);
    let relation = sub(g, &s);
    let inverse = unrad_branch(&relation, x, depth + 1)
        .or_else(|| super::invert(&relation, x, depth + 1))?;
    let mut out: Vec<Expr> = Vec::new();
    for i in &inverse {
        for r in &roots {
            let v = subs1(i, &s, r);
            if !out.contains(&v) {
                out.push(v);
            }
        }
    }
    if out.len() > 2 {
        record(&out);
    }
    Some(out)
}

/// Whether `unrad(numer(together(eq)), x)` raises NotImplementedError: the
/// only case where the established resultant elimination of the Python
/// wrapper (`_unrad_resultant`) can still reduce an equation that the
/// native `_solve` left unsolved.
pub(super) fn unrad_raises_ni(eq: &Expr, x: &Expr) -> bool {
    let (num, _) = crate::simplify::numer_denom(&crate::simplify::together(eq));
    let mut host = Host { depth: 8 };
    matches!(
        unrad::unrad(&mut host, &num, std::slice::from_ref(x), false, false),
        Err(UFail::NotImplemented(_))
    )
}

/// The unrad fallback of `_solve`; None when `unrad` raises or gives
/// nothing, or when the radical-free equation is not solved (the reference
/// then tries `_tsolve`).
fn unrad_branch(num: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    if depth > 8 {
        return None;
    }
    let mut host = Host { depth };
    let (eq, cov) = match unrad::unrad(&mut host, num, std::slice::from_ref(x), false, false) {
        Ok(Some(u)) => u,
        // `unrad` returned None, or raised: NotImplementedError and
        // ValueError take the reference to `_tsolve`; the other errors
        // (PolynomialError, IndexError) are declined to the same steps.
        _ => return None,
    };
    let values = if let [isym, ieq] = cov.as_slice() {
        // `inv = _vsolve(ieq, symbol, **flags)[0]`
        let inverse = nested_solve(ieq, x, depth)?;
        let inverse = inverse.first()?.clone();
        let mut out: Vec<Expr> = Vec::new();
        for xi in nested_solve(&eq, isym, depth)? {
            let v = subs1(&inverse, isym, &xi);
            if !out.contains(&v) {
                out.push(v);
            }
        }
        out
    } else {
        nested_solve(&eq, x, depth)?
    };
    record(&values);
    Some(values)
}

/// A nested `_vsolve(eq, symbol, **flags)` of the unrad fallback, with the
/// default `simplify` of the reference `_solve`: every factor of a product
/// is solved by its own `_solve` call, whose solutions are simplified unless
/// there are more than two of them (the high-order roots of the polynomial
/// branch). The candidates are checked by the outer call.
fn nested_solve(eq: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    let factors: Vec<Expr> = if eq.is_mul() {
        eq.args().into_iter().filter(|f| has(f, x)).collect()
    } else {
        vec![eq.clone()]
    };
    let mut out: Vec<Expr> = Vec::new();
    for factor in factors {
        let mut found = nested_checked(&factor, x, depth + 1)?;
        if found.len() <= 2 {
            found = found.iter().map(crate::simplify::simplify).collect();
        }
        for v in found {
            if !out.contains(&v) {
                out.push(v);
            }
        }
    }
    Some(out)
}
