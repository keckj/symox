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

/// The outcome of the reference `checksol` on one equation.
#[derive(Clone, Debug, PartialEq)]
pub enum Verdict {
    /// True, False or None, returned without a warning.
    Done(Option<bool>),
    /// None at the end of the attempt loop: the reference warns ("could not
    /// verify solution") when called with `warn=True`.
    Unverified,
    /// The substituted value is not an Expr (a Relational other than `Eq`,
    /// a Boolean or a set): the reference raises `AttributeError: '<type>'
    /// object has no attribute 'is_constant'` on it.
    NoIsConstant(Expr),
    /// Not decided here the way the reference would decide it (a number
    /// that does not evaluate numerically); the caller falls back.
    Decline(&'static str),
}

impl Verdict {
    /// The internal (solve) reading: everything undecided is None.
    pub fn result(&self) -> Option<bool> {
        match self {
            Verdict::Done(r) => *r,
            _ => None,
        }
    }
}

/// The end of one attempt: a Rational is decided exactly, a number
/// numerically (with `numerical`).
fn decide(val: &Expr, numerical: bool) -> Option<Verdict> {
    if val.is_rational_number() {
        return Some(Verdict::Done(Some(val.is_zero())));
    }
    if numerical && is_number(val) {
        return Some(match numeric_small(val) {
            Some(r) => Verdict::Done(Some(r)),
            None => Verdict::Decline("numerical value not evaluated"),
        });
    }
    None
}

/// Not an `Expr` in the reference (no `is_constant`, `is_number` False).
fn is_non_expr(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Rel(..) | Kind::Bool(_))
        || matches!(e.kind(), Kind::Func(Func::Logic(_) | Func::Set(_) | Func::Predicate(_) | Func::AppliedPredicate, _))
}

/// `Mul.make_args(denom(p))` of every Expr node `p` (SymPy's `denoms(eq)`
/// without symbols).
fn all_denoms(f: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(f, &mut |n| {
        if is_non_expr(n) {
            return;
        }
        let d = crate::simptools::radsimp::fraction(n, false).1;
        if !d.is_one() {
            for m in if d.is_mul() { crate::traverse::ordered_args(&d) } else { vec![d] } {
                if !out.contains(&m) {
                    out.push(m);
                }
            }
        }
    });
    out
}

/// SymPy 1.14 `Expr.is_constant(*wrt, simplify=simplify)`. The random
/// real and complex evaluations of the reference are left out: they can
/// only return False, and the derivative test that follows them returns
/// False for every expression they would catch.
pub fn is_constant_ref(e: &Expr, wrt: &[Expr], simplify: bool) -> Option<bool> {
    is_constant_with(e, wrt, simplify, &mut native_simplify).unwrap_or(None)
}

/// The simplifier of the native callers: `crate::simplify::simplify`.
pub fn native_simplify(e: &Expr) -> Option<Expr> {
    Some(crate::simplify::simplify(e))
}

/// A simplifier: the reference `simplify(e)`, or None when it failed (the
/// whole check is then declined).
pub type Simplifier<'a> = dyn FnMut(&Expr) -> Option<Expr> + 'a;

/// `is_constant_ref` with the simplifier `simp`; Err when it failed.
pub fn is_constant_with(e: &Expr, wrt: &[Expr], simplify: bool, simp: &mut Simplifier<'_>) -> Result<Option<bool>, ()> {
    if is_number(e) {
        return Ok(Some(true));
    }
    let free = free_symbols(e);
    if free.is_empty() {
        return Ok(Some(true));
    }
    let mut keys: Vec<Expr> = Vec::new();
    for w in wrt {
        if !keys.contains(w) {
            keys.push(w.clone());
        }
    }
    if !keys.is_empty() && !keys.iter().any(|w| free.contains(w)) {
        return Ok(Some(true));
    }
    let wrt: Vec<Expr> = if keys.is_empty() { free.iter().cloned().collect() } else { keys };
    let expr = if simplify { simp(e).ok_or(())? } else { e.clone() };
    if is_zero(&expr) == Some(true) {
        return Ok(Some(true));
    }
    let same = wrt.len() == free.len() && wrt.iter().all(|w| free.contains(w));
    if same {
        let at = |v: i64| {
            let pairs: Vec<(Expr, Expr)> = free.iter().map(|s| (s.clone(), Expr::int(v))).collect();
            crate::subs::simultaneous(&expr, &pairs)
        };
        let a = at(0);
        if !a.is_constant(crate::expr::Constant::NaN) {
            let b = at(1);
            if !b.is_constant(crate::expr::Constant::NaN) && crate::simptools::methods::equals(&b, &a) == Some(false) {
                return Ok(Some(false));
            }
        }
    }
    for w in &wrt {
        let mut deriv = crate::diff::diff(&expr, w);
        if simplify {
            deriv = simp(&deriv).ok_or(())?;
        }
        if !deriv.is_zero() {
            return Ok(Some(false));
        }
    }
    // fuzzy_not(fuzzy_or(den.is_zero for den in denoms(self)))
    let mut any_none = false;
    for d in all_denoms(e) {
        match is_zero(&d) {
            Some(true) => return Ok(Some(false)),
            None => any_none = true,
            Some(false) => {}
        }
    }
    Ok(if any_none { None } else { Some(true) })
}

/// `f.subs(sol)` for a dict `sol`: SymPy sorts an unordered mapping so that
/// larger keys come first (`-_nodes`, then `default_sort_key`) and applies
/// the pairs in turn.
pub fn subs_mapping(f: &Expr, sol: &[(Expr, Expr)]) -> Expr {
    if sol.len() == 1 {
        return subs1(f, &sol[0].0, &sol[0].1);
    }
    let mut pairs = sol.to_vec();
    pairs.sort_by_cached_key(|(old, _)| (std::cmp::Reverse(crate::traverse::node_count(old)), crate::sortkey::sort_key(old)));
    subs_pairs(f, &pairs)
}

/// SymPy 1.14 `checksol(f, sol, **flags)` for one equation `f` and a
/// solution mapping `sol` (in the caller's dict order). With `simplify`,
/// the reference replaces the values of `sol` by their simplified form
/// when it reaches attempt 2 (`sol[k] = simplify(sol[k])`); this is done in
/// place here and `simplified` is set.
pub fn checksol_mapping(f: &Expr, sol: &mut [(Expr, Expr)], flags: CheckFlags, simplified: &mut bool) -> Verdict {
    checksol_mapping_with(f, sol, flags, simplified, &mut native_simplify)
}

/// `checksol_mapping` with the simplifier `simp` (the public binding passes
/// the public `simplify`, as the reference calls it).
pub fn checksol_mapping_with(
    f: &Expr, sol: &mut [(Expr, Expr)], flags: CheckFlags, simplified: &mut bool, simp: &mut Simplifier<'_>,
) -> Verdict {
    const FAILED: Verdict = Verdict::Decline("simplify failed");
    if !is_non_expr(f) && is_number(f) {
        return Verdict::Done(is_zero(f));
    }
    let mut f = f.clone();
    if let Kind::Rel(op @ (crate::expr::RelOp::Eq | crate::expr::RelOp::Ne), l, r) = f.kind() {
        let (op, mut l, mut r) = (*op, l.clone(), r.clone());
        if matches!(r.kind(), Kind::Bool(_)) {
            std::mem::swap(&mut l, &mut r);
            f = Expr::from_kind(Kind::Rel(op, l.clone(), r.clone()));
        }
        if matches!(l.kind(), Kind::Bool(_)) {
            f = subs_mapping(&f, sol);
            // SymPy's `is_Boolean`: only BooleanAtom and BooleanFunction
            // (a Relational, a Symbol or Contains are not)
            if !matches!(f.kind(), Kind::Bool(_) | Kind::Func(Func::Logic(_), _)) {
                return Verdict::Done(None);
            }
        } else if op == crate::expr::RelOp::Eq {
            f = crate::ops::sub(&l, &r);
        }
    }
    if let Kind::Bool(b) = f.kind() {
        return Verdict::Done(Some(*b));
    }
    if !matches!(f.kind(), Kind::Rel(..)) && f.is_zero() {
        return Verdict::Done(Some(true));
    }
    if sol.iter().any(|(_, v)| has_illegal(v)) {
        return Verdict::Done(Some(false));
    }
    let keys: Vec<Expr> = sol.iter().map(|(k, _)| k.clone()).collect();
    // attempt 0
    let mut val = subs_mapping(&f, sol);
    if has_illegal(&val) {
        return Verdict::Done(Some(false));
    }
    if is_non_expr(&val) {
        return Verdict::NoIsConstant(val);
    }
    if let Some(r) = decide(&val, flags.numerical) {
        return r;
    }
    // attempt 1
    if !is_number(&val) {
        match is_constant_with(&val, &keys, !flags.minimal, simp) {
            Ok(Some(true)) => {}
            Ok(_) => return Verdict::Done(Some(false)),
            Err(()) => return FAILED,
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
        return Verdict::Done(None);
    }
    if flags.simplify {
        let mut new = Vec::with_capacity(sol.len());
        for (_, v) in sol.iter() {
            match simp(v) {
                Some(s) => new.push(s),
                None => return FAILED,
            }
        }
        for ((_, v), s) in sol.iter_mut().zip(new) {
            *v = s;
        }
        *simplified = true;
    }
    val = match simp(&subs_mapping(&f, sol)) {
        Some(v) => v,
        None => return FAILED,
    };
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
            Kind::Piecewise(_) => saw = true,
            _ if is_function(p) => saw = true,
            _ => {}
        }
    });
    if !saw {
        return Verdict::Done(Some(false));
    }
    if flags.force {
        val = subs_pairs(&val, &reps);
    }
    let Some(nz) = is_zero(&val).map(|z| !z) else {
        return Verdict::Unverified;
    };
    let mut lambert = false;
    crate::traverse::preorder(&val, &mut |e| {
        if matches!(e.kind(), Kind::Func(Func::Special(s), _) if s.name() == "LambertW") {
            lambert = true;
        }
    });
    if is_number(&val) && lambert {
        return Verdict::Done(None);
    }
    Verdict::Done(Some(!nz))
}

/// The result of the reference `checksol` on a list of equations.
#[derive(Clone, Debug)]
pub struct SystemOutcome {
    /// `fuzzy_and` of the per-equation results, or the first error or
    /// decline.
    pub verdict: Verdict,
    /// One entry per equation that ended undecided at the end of the
    /// attempt loop (a warning with `warn=True`): the solution values at
    /// that moment when they had been simplified, None when they were
    /// still the caller's.
    pub warnings: Vec<Option<Vec<Expr>>>,
    /// Whether the values of the mapping were replaced by their
    /// simplified form (they are in the `sol` argument).
    pub simplified: bool,
}

/// SymPy 1.14 `checksol(fs, sol, **flags)` for an iterable of equations:
/// `fuzzy_and(checksol(fi, sol, **flags) for fi in fs)`, which stops at the
/// first False. The mapping is shared, so a simplification made for one
/// equation is seen by the next ones.
pub fn checksol_system(fs: &[Expr], sol: &mut [(Expr, Expr)], flags: CheckFlags) -> SystemOutcome {
    checksol_system_with(fs, sol, flags, &mut native_simplify)
}

/// `checksol_system` with the simplifier `simp`.
pub fn checksol_system_with(fs: &[Expr], sol: &mut [(Expr, Expr)], flags: CheckFlags, simp: &mut Simplifier<'_>) -> SystemOutcome {
    let mut simplified = false;
    let mut warnings = Vec::new();
    let mut rv = Some(true);
    for f in fs {
        let v = checksol_mapping_with(f, sol, flags, &mut simplified, simp);
        let r = match v {
            Verdict::Done(r) => r,
            Verdict::Unverified => {
                warnings.push(simplified.then(|| sol.iter().map(|(_, v)| v.clone()).collect()));
                None
            }
            other => return SystemOutcome { verdict: other, warnings, simplified },
        };
        if r == Some(false) {
            return SystemOutcome { verdict: Verdict::Done(Some(false)), warnings, simplified };
        }
        if rv == Some(true) {
            rv = r;
        }
    }
    SystemOutcome { verdict: Verdict::Done(rv), warnings, simplified }
}

/// SymPy 1.14 `checksol(f, sym, sol, **flags)` for one Expr `f`.
pub fn checksol(f: &Expr, sym: &Expr, sol: &Expr, flags: CheckFlags) -> Option<bool> {
    let mut pairs = [(sym.clone(), sol.clone())];
    let mut simplified = false;
    checksol_mapping(f, &mut pairs, flags, &mut simplified).result()
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
