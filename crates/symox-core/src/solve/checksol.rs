//! SymPy 1.14 `checksol(f, {symbol: sol}, **flags)` for one expression and
//! one symbol, with the reference flags `numerical`, `force`, `simplify`
//! and `minimal` (all read with the reference defaults by the callers).
//!
//! The attempts follow the reference loop:
//! 0. substitution; an illegal value (`nan`, `zoo`, `oo`, `-oo`) is False;
//! 1. a value that still has free symbols must be constant in `symbol`;
//!    the content is removed and the numerator expanded;
//! 2. (not with `minimal`) the candidate (with `simplify`) and the
//!    substituted value are simplified; with `force` the value is posified
//!    and expanded;
//! 3. a value without non-integer powers or functions cannot be zero; else
//!    the assumption system decides `is_zero` (on the unposified value).
//!
//! After each of attempts 0-2 a Rational value is decided exactly and a
//! numerical value (with `numerical`) by the reference absolute test
//! `abs(val.n(18).n(12, chop=True)) < 1e-9`.
//!
//! Where the reference would need `is_constant` with the symbol still
//! present, or a numerical value does not evaluate to a number, None
//! (undecided) is returned; the reference may return False there. A False
//! is only returned on a ground the reference also returns False on.

use crate::expr::{Expr, Func, Kind};
use crate::subs::subs1;
use crate::traverse::{free_symbols, has};

/// The reference `checksol` flags.
#[derive(Clone, Copy, Debug)]
pub struct CheckFlags {
    pub numerical: bool,
    pub force: bool,
    pub simplify: bool,
    pub minimal: bool,
}

impl Default for CheckFlags {
    fn default() -> Self {
        CheckFlags { numerical: true, force: true, simplify: true, minimal: false }
    }
}

fn has_illegal(e: &Expr) -> bool {
    has(e, &Expr::nan()) || has(e, &Expr::zoo()) || has(e, &Expr::oo()) || has(e, &Expr::neg_oo())
}

fn is_zero(e: &Expr) -> Option<bool> {
    crate::assumptions::settle_nonreal(e, crate::assumptions::is_zero(e))
}

fn is_number(e: &Expr) -> bool {
    crate::symbolic_number::is_symbolic_number(e)
}

fn mexpand(e: &Expr) -> Expr {
    crate::solve::tsolve::mexpand(e, true)
}

/// `p.is_Function` of the reference traversal (an applied function class,
/// not an unevaluated calculus object or a container).
fn is_function(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(f, _) if !matches!(f, Func::GenericExpr | Func::SymbolicDerivative | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _)) && !matches!(f, Func::Undefined(head) if head.is_internal("Tuple")))
}

fn subs_pairs(e: &Expr, reps: &[(Expr, Expr)]) -> Expr {
    if reps.is_empty() {
        e.clone()
    } else {
        crate::subs::subs(e, reps)
    }
}

/// `abs(val.n(18).n(12, chop=True)) < 1e-9`, None when the value does not
/// evaluate to a number.
pub fn numeric_small(val: &Expr) -> Option<bool> {
    let p18 = crate::bigfloat::dps_to_prec(18);
    let p12 = crate::bigfloat::dps_to_prec(12);
    let cutoff = -(p12 as i128);
    // The value is evaluated at 18 digits under the 12-digit chop cutoff
    // directly, so that an exact-zero residual is not chased through the
    // whole cancellation budget before being chopped. A chopped part is
    // proven below the cutoff (far below 1e-9) and an unchopped one is
    // resolved to 18 digits, so the decision is the one of the two-step
    // evaluation, which is kept when this one gives no number.
    let direct = crate::evalf_mp::evalf_chop(val, p18, cutoff);
    let v = if free_symbols(&direct).is_empty() && !has_illegal(&direct) && crate::evalf::eval_c64(&direct).is_some() {
        crate::evalf_mp::evalf_chop(&direct, p12, cutoff)
    } else {
        let v = crate::evalf_mp::evalf(val, p18);
        crate::evalf_mp::evalf_chop(&v, p12, cutoff)
    };
    if !free_symbols(&v).is_empty() || has_illegal(&v) {
        return None;
    }
    let c = crate::evalf::eval_c64(&v)?;
    let m = c.norm();
    if m.is_nan() {
        return None;
    }
    Some(m < 1e-9)
}

/// The end of one attempt: a Rational is decided exactly, a number
/// numerically (with `numerical`).
fn decide(val: &Expr, numerical: bool) -> Option<Option<bool>> {
    if val.is_rational_number() {
        return Some(Some(val.is_zero()));
    }
    if numerical && is_number(val) {
        return Some(numeric_small(val));
    }
    None
}

/// SymPy 1.14 `checksol(f, sym, sol, **flags)` for one Expr `f`.
pub fn checksol(f: &Expr, sym: &Expr, sol: &Expr, flags: CheckFlags) -> Option<bool> {
    if is_number(f) {
        return is_zero(f);
    }
    if f.is_zero() {
        return Some(true);
    }
    if has_illegal(sol) {
        return Some(false);
    }
    // attempt 0
    let mut val = subs1(f, sym, sol);
    if has_illegal(&val) {
        return Some(false);
    }
    if let Some(r) = decide(&val, flags.numerical) {
        return r;
    }
    // attempt 1
    if !is_number(&val) {
        if free_symbols(&val).contains(sym) {
            // `is_constant(sym)` with the symbol still present is not
            // decided here
            return None;
        }
        if let Ok((_, prim)) = crate::decompose::as_content_primitive(&val, false, true) {
            val = prim;
        }
        val = mexpand(&crate::simptools::as_numer_denom(&val).0);
    }
    if let Some(r) = decide(&val, flags.numerical) {
        return r;
    }
    // attempt 2
    if flags.minimal {
        return None;
    }
    let sol = if flags.simplify { crate::simplify::simplify(sol) } else { sol.clone() };
    val = crate::simplify::simplify(&subs1(f, sym, &sol));
    let mut reps: Vec<(Expr, Expr)> = Vec::new();
    if flags.force {
        let (v, r) = crate::simptools::misc::posify(&val);
        reps = r;
        val = v;
        let ex = mexpand(&val);
        if is_number(&ex) {
            val = ex;
        }
    }
    if let Some(r) = decide(&val, flags.numerical) {
        return r;
    }
    // attempt 3
    let mut saw = false;
    crate::traverse::preorder(&crate::simp::expand_mul_deep(&val, true), &mut |p| {
        if saw {
            return;
        }
        match p.kind() {
            Kind::Pow(_, x) if !x.is_integer_number() => saw = true,
            _ if is_function(p) => saw = true,
            _ => {}
        }
    });
    if !saw {
        return Some(false);
    }
    if flags.force {
        val = subs_pairs(&val, &reps);
    }
    let nz = is_zero(&val).map(|z| !z)?;
    let mut lambert = false;
    crate::traverse::preorder(&val, &mut |e| {
        if matches!(e.kind(), Kind::Func(Func::Special(s), _) if s.name() == "LambertW") {
            lambert = true;
        }
    });
    if is_number(&val) && lambert {
        return None;
    }
    Some(!nz)
}

/// The check of the reference `_solve` (`checksol(f_num, {x: s},
/// simplify=False) is not False`) for a candidate whose substituted value
/// still has free symbols (a parametric candidate; the numerical ones are
/// verified by the caller): true when `checksol` proves it invalid.
/// `x**(y/2) + x**y - 6` drops `3**(2/y)` as in the reference.
pub fn rejects_parametric(f_num: &Expr, x: &Expr, s: &Expr) -> bool {
    if free_symbols(&subs1(f_num, x, s)).is_empty() {
        return false;
    }
    checksol(f_num, x, s, CheckFlags { simplify: false, ..Default::default() }) == Some(false)
}
