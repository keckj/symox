//! The `_tsolve` fallback of SymPy 1.14's `_solve` with native recursive
//! calls: the `_tsolve` steps of [`super::tsolve`] run with a host whose
//! `solve`, `nsimplify`, `equals` and `count_ops` are native, so that the
//! public `solve` no longer goes back to Python for them.
//!
//! The binding runs the native solver inside [`scope`]. Within it,
//!
//! * the native `_solve` takes the unrad fallback for radicals mixed with
//!   other functions of the unknown (`radical::solve_mixed`), between the
//!   generator branch and the `_tsolve` steps, as the reference does;
//! * the nested `_vsolve` calls of that fallback and of the `_tsolve` steps
//!   end with the `_tsolve` steps ([`fallback`]);
//! * the top-level call ends with [`fallback`] too.
//!
//! Every computation the native host cannot reproduce exactly *declines*:
//! a nested equation that the Python wrapper of `solve` would rewrite first
//! (Floats, hyperbolic functions, `Abs`, `Piecewise`, held derivatives,
//! integrals, sums or products, undefined functions), a nested equation
//! that the native chain cannot solve while the Python fallbacks could
//! still apply (a radical of the unknown: `_unrad_resultant`), `nsimplify`
//! of a number that is neither an integer nor a small rational
//! (`mpmath.identify`), an undecided `equals`, and every error other than
//! the PolynomialError the reference catches. A decline is recorded for
//! the whole call; the binding then solves again without these branches
//! (the previous native solver) and leaves the rest to the Python path, so
//! a decline never turns into an empty or partial answer.

use std::cell::{Cell, RefCell};

use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::traverse::{free_symbols, has, preorder};

use super::tsolve::{self as ts, Fail, Flags, Host, Options};

/// Computations that the native host takes from the embedding layer:
/// `mpmath.identify` for `nsimplify` of a number (`nsimplify`) and the
/// `minimal_polynomial` tail of `Expr.equals` (`equals_zero`, reached only
/// when the native `equals` is undecided). A hook returning None declines.
pub struct Hooks {
    pub nsimplify: Box<dyn Fn(&Expr) -> Option<Expr>>,
    pub equals_zero: Box<dyn Fn(&Expr) -> Option<Option<bool>>>,
}

thread_local! {
    static HOOKS: RefCell<Option<std::rc::Rc<Hooks>>> = const { RefCell::new(None) };
    static ENABLED: Cell<bool> = const { Cell::new(false) };
    static DECLINED: Cell<bool> = const { Cell::new(false) };
    /// The reference `tsolve_saw` flag: equations whose `_tsolve` steps
    /// are running (an equation met again is not handled again).
    static SAW: RefCell<Vec<(Expr, Expr)>> = const { RefCell::new(Vec::new()) };
}

/// Nesting bound of the `_tsolve` steps (the reference relies on
/// `tsolve_saw` alone; deeper chains decline to the Python path).
const MAX_SAW: usize = 12;

/// The flags of a `solve` call that the native fallback reads.
#[derive(Clone, Copy, Debug)]
pub struct SolveOptions {
    pub check: bool,
    /// `None`: the flag was not given.
    pub simplify: Option<bool>,
    pub force: bool,
    pub bivariate: bool,
}

impl Default for SolveOptions {
    fn default() -> Self {
        SolveOptions { check: true, simplify: None, force: true, bivariate: true }
    }
}

/// Run `f` with (`enabled`) or without the native fallbacks; returns its
/// value and whether a computation declined to the Python path.
pub fn scope<T>(enabled: bool, f: impl FnOnce() -> T) -> (T, bool) {
    let previous = ENABLED.with(|c| c.replace(enabled));
    let declined = DECLINED.with(|c| c.replace(false));
    let saw = SAW.with(|s| std::mem::take(&mut *s.borrow_mut()));
    let r = f();
    let now = DECLINED.with(|c| c.replace(declined));
    SAW.with(|s| *s.borrow_mut() = saw);
    ENABLED.with(|c| c.set(previous));
    (r, now)
}

/// Run `f` with the `hooks` of the embedding layer.
pub fn with_hooks<T>(hooks: Hooks, f: impl FnOnce() -> T) -> T {
    let previous = HOOKS.with(|h| h.replace(Some(std::rc::Rc::new(hooks))));
    let r = f();
    HOOKS.with(|h| *h.borrow_mut() = previous);
    r
}

fn hooks() -> Option<std::rc::Rc<Hooks>> {
    HOOKS.with(|h| h.borrow().clone())
}

pub(super) fn enabled() -> bool {
    ENABLED.with(Cell::get)
}

fn decline_why(why: &str, e: &Expr) {
    if std::env::var_os("SYMOX_TSOLVE_TRACE").is_some() {
        eprintln!("native tsolve declines ({why}): {e}");
    }
    decline();
}

fn decline() {
    DECLINED.with(|c| c.set(true));
}

fn declined() -> bool {
    DECLINED.with(Cell::get)
}

/// A Python-only computation was needed.
pub struct Decline;

type R<T> = Result<T, Fail<Decline>>;

fn is_hyperbolic(f: &Func) -> bool {
    matches!(f, Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth) || matches!(f.name(), "sech" | "csch")
}

/// Nodes that the Python layer of `solve` handles before or after the
/// native solver: `Piecewise` (`_solve_piecewise`), held derivatives,
/// integrals, sums, products and substitutions (`_masked_solution`),
/// undefined functions (`inverse()` of host classes), `Abs`
/// (`_solve_rewrite_abs`), hyperbolic functions and Floats (the rewrites
/// of the public `solve`).
fn python_layer(e: &Expr, x: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if found {
            return;
        }
        found = match s.kind() {
            Kind::Piecewise(_) | Kind::Derivative(..) => true,
            Kind::Num(Number::Float(_)) => true,
            Kind::Func(f, _) => {
                matches!(f, Func::Subs | Func::SymbolicDerivative | Func::Integral | Func::Sum | Func::Product
                    | Func::Undefined(_) | Func::Abs | Func::Limit)
                    || (is_hyperbolic(f) && has(s, x))
            }
            _ => false,
        };
    });
    found || e.float_info().max_precision.is_some()
}

/// A power containing the unknown whose exponent leads with a Rational of
/// denominator `q > 1` (`_Q(pow) != 1` in the `_take` of the reference
/// `unrad`: `sqrt(x)`, `x**(y/2)`): the equations that `unrad` (and the
/// Python `_unrad_solution`) can take.
pub(super) fn has_radical(e: &Expr, x: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if let Kind::Pow(_, n) = s.kind() {
            let c = match n.as_number() {
                Some(c) => c.clone(),
                None => crate::ops::as_coeff_mul(n).0,
            };
            if c.is_rational() && !c.is_integer() && has(s, x) {
                found = true;
            }
        }
    });
    found
}

/// The number of distinct radicals of the unknown (`_radical_generators`).
fn radical_generators(e: &Expr, x: &Expr) -> usize {
    let mut seen: Vec<Expr> = Vec::new();
    preorder(e, &mut |s| {
        if let Kind::Pow(b, n) = s.kind() {
            if n.as_number().is_some_and(|n| n.is_rational() && !n.is_integer()) && has(b, x) && !seen.contains(s) {
                seen.push(s.clone());
            }
        }
    });
    seen.len()
}

struct NativeHost {
    opts: SolveOptions,
}

impl Host for NativeHost {
    type Error = Decline;

    fn solve(&mut self, eq: &Expr, x: &Expr, flags: Flags) -> R<Vec<Expr>> {
        let opts = match flags {
            Flags::Plain => SolveOptions::default(),
            Flags::Same => self.opts,
            Flags::NoBivariate => SolveOptions { bivariate: false, ..self.opts },
            Flags::Popped => SolveOptions { bivariate: true, ..self.opts },
        };
        nested_solve(eq, x, opts)
    }

    /// SymPy's `nsimplify(e)`: with free symbols, only Floats change
    /// (`_real_to_rational`); a number goes through `pslq` at 15 digits,
    /// which keeps a small rational, and then `mpmath.identify`, which is
    /// not native (the `nsimplify` hook).
    fn nsimplify(&mut self, e: &Expr) -> R<Expr> {
        if e.float_info().max_precision.is_some() {
            decline_why("nsimplify float", e);
            return Err(Fail::Host(Decline));
        }
        if !free_symbols(e).is_empty() {
            return Ok(e.clone());
        }
        if let Some(n) = e.as_number() {
            if n.is_integer() {
                return Ok(e.clone());
            }
            if n.is_rational() {
                let small = |v: Number| v.to_i64().is_some_and(|v| v.abs() <= 1000);
                if small(n.numer()) && small(n.denom()) {
                    return Ok(e.clone());
                }
            }
        }
        match hooks().and_then(|h| (h.nsimplify)(e)) {
            Some(v) => Ok(v),
            None => {
                decline_why("nsimplify hook", e);
                Err(Fail::Host(Decline))
            }
        }
    }

    fn equals_zero(&mut self, e: &Expr) -> R<Option<bool>> {
        if let Some(v) = crate::simptools::methods::equals(e, &Expr::zero()) {
            return Ok(Some(v));
        }
        match hooks().and_then(|h| (h.equals_zero)(e)) {
            Some(v) => Ok(v),
            None => {
                decline_why("equals hook", e);
                Err(Fail::Host(Decline))
            }
        }
    }

    fn count_ops(&mut self, e: &Expr) -> usize {
        crate::simptools::methods::count_ops(e)
    }

    /// The inverse of a host function class is a Python method.
    fn inverse(&mut self, _f: &Expr, _rhs: &Expr) -> R<Option<Expr>> {
        decline_why("inverse", _f);
        Err(Fail::Host(Decline))
    }
}

/// A nested `solve(eq, x, **flags)` of the `_tsolve` steps: the native
/// solver, then [`fallback`], then the simplification of the public
/// `solve` (`_simplified_solutions`).
fn nested_solve(eq: &Expr, x: &Expr, opts: SolveOptions) -> R<Vec<Expr>> {
    let eq = super::as_expr(eq);
    if python_layer(&eq, x) {
        decline_why("python layer (nested)", &eq);
        return Err(Fail::Host(Decline));
    }
    let saved = super::take_unsimplified();
    let found = super::without_assumption_filter(|| super::without_unrad(false, || {
        if !opts.bivariate {
            // `solve(..., bivariate=False)` is not taken by the native
            // binding: the reference flag only reaches `_tsolve`.
            return fallback(&eq, x, opts);
        }
        super::solve_atom(&eq, x, opts.check)
            .or_else(|| super::solve_with_check(&eq, x, opts.check))
            .or_else(|| fallback(&eq, x, opts))
    }));
    let unsimplified = super::take_unsimplified();
    super::radical::record(&saved);
    if declined() {
        return Err(Fail::Host(Decline));
    }
    match found {
        Some(values) => Ok(simplified(&eq, x, values, &unsimplified, opts.simplify)),
        // SymPy raises here; the Python wrapper still has the resultant
        // elimination when `unrad` raises NotImplementedError.
        None if has_radical(&eq, x) && super::radical::unrad_raises_ni(&eq, x) => {
            decline_why("nested radical unsolved", &eq);
            Err(Fail::Host(Decline))
        }
        None => Err(Fail::NotImplemented(None, format!("could not solve {eq} for {x}"))),
    }
}

/// `_simplified_solutions` of the public `solve` for a nested call.
fn simplified(eq: &Expr, x: &Expr, values: Vec<Expr>, unsimplified: &[Expr], flag: Option<bool>) -> Vec<Expr> {
    if flag == Some(false) {
        return values;
    }
    if flag.is_none() && radical_generators(eq, x) > 1 {
        return values;
    }
    // more than two roots of one generator stay as they are
    // (`_high_order_roots`), as the unrad solutions
    let (num, _) = crate::simplify::numer_denom(&crate::simplify::together(eq));
    let high_order = flag.is_none() && !num.is_mul() && values.len() > 2;
    let mut out: Vec<Expr> = Vec::new();
    for v in values {
        let keep = high_order || (flag.is_none() && unsimplified.contains(&v));
        let v = if keep { v } else { crate::simplify::simplify(&v) };
        if !out.contains(&v) {
            out.push(v);
        }
    }
    out.sort_by_cached_key(crate::sortkey::sort_key);
    out
}

/// The `_tsolve` steps for `eq = 0` in the Symbol `x` with the native host:
/// the solutions, or None when no step applies or when the computation
/// declined (recorded for [`scope`]).
pub fn fallback(eq: &Expr, x: &Expr, opts: SolveOptions) -> Option<Vec<Expr>> {
    if !enabled() || x.as_sym().is_none() {
        return None;
    }
    let f = super::as_expr(eq);
    if !has(&f, x) {
        return None;
    }
    if python_layer(&f, x) {
        decline_why("python layer", &f);
        return None;
    }
    let key = (f.clone(), x.clone());
    let depth = SAW.with(|s| s.borrow().len());
    if SAW.with(|s| s.borrow().contains(&key)) {
        return None;
    }
    if depth >= MAX_SAW {
        decline_why("depth", &f);
        return None;
    }
    SAW.with(|s| s.borrow_mut().push(key));
    let mut host = NativeHost { opts };
    let options = Options { force: opts.force, bivariate: opts.bivariate, check: opts.check, numerical: true };
    let r = ts::tsolve_step(&mut host, &f, x, options);
    SAW.with(|s| s.borrow_mut().pop());
    match r {
        Ok(v) => v,
        // the reference `_solve` catches PolynomialError around `_tsolve`
        Err(Fail::Poly(_)) => None,
        // NotImplementedError leaves the Python fallbacks of `solve`
        // (masking, host inverses, `_unrad_solution`) to try
        Err(Fail::NotImplemented(..)) if !has_radical(&f, x) => None,
        Err(e) => {
            let why = match &e {
                Fail::NotImplemented(_, m) => format!("NotImplementedError with a radical: {m}"),
                Fail::Value(_, m) => format!("ValueError: {m}"),
                Fail::Type(m) => format!("TypeError: {m}"),
                Fail::Poly(_) => "PolynomialError".to_string(),
                Fail::Host(_) => "host".to_string(),
            };
            decline_why(&why, &f);
            None
        }
    }
}
