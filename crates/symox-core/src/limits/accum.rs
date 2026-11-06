//! Accumulation bounds and the limits whose value is one.
//!
//! SymPy 1.14 (sympy/calculus/accumulationbounds.py) represents the closed
//! real interval of the values an oscillating expression keeps taking as
//! `AccumBounds(min, max)`, an expression class with its own arithmetic.
//! The class itself is a host class (its Python facade carries the
//! printing, the assumptions and the `isinstance` identity); its nodes are
//! applications of the registered host head `AccumulationBounds(min, max)`.
//! Everything computed with them lives here:
//!
//! * the class methods `__new__`, `__add__`, `__sub__`, `__rsub__`,
//!   `__neg__`, `__mul__`, `__truediv__`, `__rtruediv__`, `__pow__`,
//!   `__rpow__`, `__abs__`, `__contains__`, `intersection`, `union`
//!   and the `_eval_is_ge` dispatch, line by line;
//! * the coefficient folding of `Add.flatten` and `Mul.flatten` for terms
//!   holding bounds;
//! * the `eval` branches of the elementary functions that return bounds
//!   (`sin(oo)` is `AccumBounds(-1, 1)`, `exp(AccumBounds(a, b))` is
//!   `AccumBounds(exp(a), exp(b))`, ...);
//! * the oscillating limit: the reference reaches
//!   `limit(sin(x), x, oo) = AccumBounds(-1, 1)` through the leading terms
//!   of sums, products, powers and functions combined with the bounds
//!   arithmetic. When the native engine declines a limit, the expression is
//!   evaluated term by term: subexpressions without an oscillating function
//!   go to the native engine, an oscillating function of an argument tending
//!   to an infinity becomes the bounds of its `eval`, and the parts combine
//!   with the bounds arithmetic. Any part left undecided, any nan and any
//!   `oo - oo` keeps the limit unevaluated.
//!
//! Python operators are emulated by [`op_add`] and its siblings: a bounds
//! operand dispatches to the bounds method (its reflected method when it is
//! on the right, as Python calls the subclass first), other operands use
//! the native arithmetic followed by the coefficient folding.

use std::sync::RwLock;

use crate::expr::{Constant, Expr, Func, Kind, RelOp, UndefinedFunction};
use crate::number::Number;
use crate::ops;
use crate::special::SpecialFn;

use super::{Dir, LimitError};

/// The Python exception a bounds computation raises.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BoundsError {
    Value(String),
    Type(String),
    NotImplemented(String),
    Overflow(String),
    Attribute(String),
    /// An exception of another class raised by a host callback; the host
    /// keeps the exception object.
    Host,
    /// The walker's internal decline (`_Fail`); never leaves the walker.
    Fail,
}

pub type BRes<T> = Result<T, BoundsError>;

fn type_error<T>(msg: impl Into<String>) -> BRes<T> { Err(BoundsError::Type(msg.into())) }
fn value_error<T>(msg: impl Into<String>) -> BRes<T> { Err(BoundsError::Value(msg.into())) }

/// Host services the bounds arithmetic needs from the Python facade.
#[derive(Clone, Copy)]
pub struct Hooks {
    /// SymPy's `Min(*args)` (`false`) or `Max(*args)` (`true`) for arguments
    /// the native lattice functions do not take: nested `Min`/`Max` nodes
    /// and host instances, which SymPy's `MinMaxBase` algorithm collapses.
    pub lattice: fn(bool, &[Expr]) -> BRes<Expr>,
    /// `type(e)` as Python prints it, for SymPy's comparison errors.
    pub describe: fn(&Expr) -> String,
}

static HEAD: RwLock<Option<UndefinedFunction>> = RwLock::new(None);
static HOOKS: RwLock<Option<Hooks>> = RwLock::new(None);

/// Register the host head of `AccumulationBounds(min, max)`.
pub fn register(head: UndefinedFunction) {
    *HEAD.write().unwrap_or_else(|p| p.into_inner()) = Some(head);
}

/// Register the host services.
pub fn set_hooks(hooks: Hooks) {
    *HOOKS.write().unwrap_or_else(|p| p.into_inner()) = Some(hooks);
}

fn hooks() -> Option<Hooks> { *HOOKS.read().unwrap_or_else(|p| p.into_inner()) }

/// Whether the bounds class is registered.
pub fn registered() -> bool { HEAD.read().unwrap_or_else(|p| p.into_inner()).is_some() }

const CLASS: &str = "AccumulationBounds";

/// Whether `e` is an `AccumulationBounds` node.
pub fn is_bounds(e: &Expr) -> bool {
    e.is_host() && crate::host::class_name(e) == Some(CLASS)
}

/// Whether `e` holds an `AccumulationBounds` node.
pub fn contains_bounds(e: &Expr) -> bool {
    if !e.has_host() { return false; }
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found |= is_bounds(s));
    found
}

/// A new bounds node `AccumulationBounds(lo, hi)` (no validation: see
/// [`ab_new`]).
fn node(lo: Expr, hi: Expr) -> BRes<Expr> {
    let guard = HEAD.read().unwrap_or_else(|p| p.into_inner());
    let Some(head) = guard.as_ref() else {
        return Err(BoundsError::NotImplemented("AccumulationBounds is not registered".into()));
    };
    let head = match head.host().and_then(|host| host.applied()) {
        Some(handle) => head.clone().with_host(handle),
        None => head.clone(),
    };
    Ok(Expr::from_kind(Kind::Func(Func::Undefined(head), vec![lo, hi])))
}

fn ends(ab: &Expr) -> (Expr, Expr) {
    let args = ab.args();
    (args[0].clone(), args[1].clone())
}

// --- Python attribute emulation ---------------------------------------------

fn oo() -> Expr { Expr::oo() }
fn noo() -> Expr { Expr::neg_oo() }
fn is_oo(e: &Expr) -> bool { e.is_constant(Constant::Infinity) }
fn is_noo(e: &Expr) -> bool { e.is_constant(Constant::NegInfinity) }
fn is_nan(e: &Expr) -> bool { e.is_constant(Constant::NaN) }
fn is_zoo(e: &Expr) -> bool { e.is_constant(Constant::ComplexInfinity) }
fn is_one(e: &Expr) -> bool { *e == Expr::one() }

fn fact(e: &Expr, name: &str) -> Option<bool> { crate::assumptions::fact(e, name) }
fn yes(e: &Expr, name: &str) -> bool { fact(e, name) == Some(true) }
fn settled(e: &Expr, v: Option<bool>) -> Option<bool> { crate::assumptions::settle_nonreal(e, v) }
fn zero(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_zero(e)) }
fn positive(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_positive(e)) }
fn negative(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_negative(e)) }
fn nonnegative(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_nonnegative(e)) }
fn integer(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_integer(e)) }
fn rational(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_rational(e)) }
fn real(e: &Expr) -> Option<bool> { settled(e, crate::assumptions::is_real(e)) }
fn finite(e: &Expr) -> Option<bool> { crate::assumptions::is_finite(e) }

/// `e.is_number`.
fn is_number(e: &Expr) -> bool { crate::symbolic_number::is_symbolic_number(e) }
/// `e.is_Number`.
fn is_number_class(e: &Expr) -> bool {
    e.is_number() || matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN))
}
/// `e.is_comparable`.
fn is_comparable(e: &Expr) -> bool { crate::integrate::piecewise::is_comparable(e) }

/// SymPy's `.args`, as the binding reports them.
pub fn python_args(e: &Expr) -> Vec<Expr> {
    if let Kind::Func(Func::Heaviside, args) = e.kind() {
        if args.len() == 1 { return vec![args[0].clone(), Expr::rational(1, 2)]; }
    }
    crate::traverse::ordered_args(e)
}

/// `bool(e)` of a Python expression.
fn truthy(e: &Expr) -> BRes<bool> {
    match e.kind() {
        Kind::Bool(b) => Ok(*b),
        Kind::Num(n) => Ok(!n.is_zero()),
        Kind::Rel(..) => type_error(format!("cannot determine truth value of Relational: {}", crate::printer::to_str(e))),
        _ => Ok(true),
    }
}

fn make_args_add(e: &Expr) -> Vec<Expr> { if e.is_add() { python_args(e) } else { vec![e.clone()] } }
fn make_args_mul(e: &Expr) -> Vec<Expr> { if e.is_mul() { python_args(e) } else { vec![e.clone()] } }

// --- constructors and operators ---------------------------------------------

/// SymPy's coefficient folding of `Add.flatten` (`add`) or `Mul.flatten`
/// for arguments holding bounds: Numbers and bounds fold into the leading
/// coefficient through the bounds arithmetic. `None` when no argument
/// holds bounds at the top, or the coefficient is not a Number or bounds.
pub fn fold(add: bool, args: &[Expr]) -> BRes<Option<Expr>> {
    let nested = |e: &Expr| if add { e.is_add() } else { e.is_mul() };
    let found = args.iter().any(|arg| {
        is_bounds(arg)
            || (nested(arg) && python_args(arg).iter().any(is_bounds))
            || (add && arg.is_mul() && python_args(arg).first().is_some_and(is_bounds))
    });
    if !found { return Ok(None); }
    let mut flat = Vec::new();
    for arg in args {
        if nested(arg) && crate::decompose::commutativity(arg) == Some(true) {
            flat.extend(python_args(arg));
        } else {
            flat.push(arg.clone());
        }
    }
    let mut coeff = if add { Expr::zero() } else { Expr::one() };
    let mut rest = Vec::new();
    for mut arg in flat {
        if add && arg.is_mul() && python_args(&arg).first().is_some_and(is_bounds) {
            // A Number times bounds distributed over a sum.
            if let Some(folded) = fold(false, &python_args(&arg))? {
                if truthy(&folded)? { arg = folded; }
            }
        }
        if is_bounds(&arg) {
            coeff = if add { ab_add(&arg, &coeff)? } else { ab_mul(&arg, &coeff)? };
        } else if is_number_class(&arg) {
            coeff = if add { op_add(&coeff, &arg)? } else { op_mul(&coeff, &arg)? };
        } else {
            rest.push(arg);
        }
        if !(is_bounds(&coeff) || is_number_class(&coeff)) { return Ok(None); }
    }
    if is_nan(&coeff) || rest.is_empty() { return Ok(Some(coeff)); }
    let mut all = vec![coeff];
    all.extend(rest);
    Ok(Some(if add { ops::add(&all) } else { ops::mul(&all) }))
}

/// The native result of an operator, folded as `Add.flatten` would.
fn native_fold(e: Expr) -> BRes<Expr> {
    if crate::traverse::has_unfolded_bounds(&e) {
        if let Some(folded) = fold(e.is_add(), &python_args(&e))? { return Ok(folded); }
    }
    Ok(e)
}

/// `a + b`.
pub fn op_add(a: &Expr, b: &Expr) -> BRes<Expr> {
    if is_bounds(a) { return ab_add(a, b); }
    if is_bounds(b) { return ab_add(b, a); }
    native_fold(ops::add(&[a.clone(), b.clone()]))
}

/// `a - b`.
pub fn op_sub(a: &Expr, b: &Expr) -> BRes<Expr> {
    if is_bounds(a) { return ab_sub(a, b); }
    if is_bounds(b) { return ab_rsub(b, a); }
    native_fold(ops::sub(a, b))
}

/// `a*b`.
pub fn op_mul(a: &Expr, b: &Expr) -> BRes<Expr> {
    if is_bounds(a) { return ab_mul(a, b); }
    if is_bounds(b) { return ab_mul(b, a); }
    native_fold(ops::mul(&[a.clone(), b.clone()]))
}

/// `a/b`.
pub fn op_div(a: &Expr, b: &Expr) -> BRes<Expr> {
    if is_bounds(a) { return ab_div(a, b); }
    if is_bounds(b) { return ab_rdiv(b, a); }
    native_fold(ops::div(a, b))
}

/// `a**b`.
pub fn op_pow(a: &Expr, b: &Expr) -> BRes<Expr> {
    if is_bounds(a) { return ab_pow(a, b); }
    if is_bounds(b) { return ab_rpow(b, a); }
    native_fold(ops::pow(a, b))
}

/// `-a`.
pub fn op_neg(a: &Expr) -> BRes<Expr> {
    if is_bounds(a) { return ab_neg(a); }
    Ok(ops::neg(a))
}

/// `Add(*v)`.
pub fn add_ctor(v: &[Expr]) -> BRes<Expr> {
    Ok(match fold(true, v)? { Some(e) => e, None => ops::add(v) })
}

/// `Mul(*v)`.
pub fn mul_ctor(v: &[Expr]) -> BRes<Expr> {
    Ok(match fold(false, v)? { Some(e) => e, None => ops::mul(v) })
}

/// `Add(*v, evaluate=False)` (`add`) or `Mul(*v, evaluate=False)`.
pub fn held(add: bool, v: Vec<Expr>) -> Expr {
    if v.len() < 2 { return if add { ops::add(&v) } else { ops::mul(&v) }; }
    if !v.iter().any(|a| a.has_held_arithmetic() || a.has_matrix()) {
        let canonical = if add { ops::add(&v) } else { ops::mul(&v) };
        if python_args(&canonical) == v { return canonical; }
    }
    if add {
        Expr::from_unevaluated_kind(Kind::Add(Number::Int(0), v.into_iter().map(|a| (a, Number::Int(1))).collect()))
    } else {
        Expr::from_unevaluated_kind(Kind::Mul(Number::Int(1), v.into_iter().map(|a| (a, Expr::one())).collect()))
    }
}

/// `Pow(a, b, evaluate=False)`.
fn held_pow(a: &Expr, b: &Expr) -> Expr { Expr::from_unevaluated_kind(Kind::Pow(a.clone(), b.clone())) }

/// Whether the native `Min`/`Max` take `arg` (SymPy's `_native_argument`
/// of the lattice classes).
fn native_lattice_argument(arg: &Expr) -> bool {
    if matches!(arg.kind(), Kind::Func(Func::Min | Func::Max, _)) || arg.has_host() { return false; }
    if is_number_class(arg) { return !is_nan(arg); }
    fact(arg, "extended_real") != Some(false) && !(is_number(arg) && !is_comparable(arg))
}

/// `Min(*args)` (`max` false) or `Max(*args)`.
pub fn lattice(max: bool, args: &[Expr]) -> BRes<Expr> {
    if !args.iter().all(native_lattice_argument) {
        if let Some(hooks) = hooks() { return (hooks.lattice)(max, args); }
    }
    let name = if max { "Max" } else { "Min" };
    for value in args {
        let numeric_nonreal = is_number(value) && real(value) == Some(false) && !is_oo(value) && !is_noo(value);
        if is_nan(value) || numeric_nonreal || fact(value, "extended_real") == Some(false) {
            return value_error(format!("{name} requires extended-real arguments"));
        }
    }
    Ok(crate::functions::func(if max { Func::Max } else { Func::Min }, args.to_vec()))
}

// --- comparisons ------------------------------------------------------------

/// The evaluated relation `a op b` of Python's ordering operators.
pub fn rel_value(op: RelOp, a: &Expr, b: &Expr) -> BRes<Expr> {
    for side in [a, b] {
        if fact(side, "extended_real") == Some(false) {
            return type_error(format!("Invalid comparison of non-real {}", crate::printer::to_str(side)));
        }
        if is_nan(side) { return type_error("Invalid NaN comparison"); }
    }
    let native = crate::relational::rel(op, a, b);
    if !(a.has_host() || b.has_host()) || !matches!(native.kind(), Kind::Rel(..)) {
        return Ok(native);
    }
    if !(contains_bounds(a) || contains_bounds(b)) {
        // Other host classes keep their own protocol in the facade.
        return Ok(native);
    }
    let decided = match op {
        RelOp::Ge => protocol_is_ge(a, b)?,
        RelOp::Le => protocol_is_ge(b, a)?,
        RelOp::Gt => protocol_is_ge(b, a)?.map(|v| !v),
        RelOp::Lt => protocol_is_ge(a, b)?.map(|v| !v),
        _ => None,
    };
    Ok(decided.map_or(native, Expr::boolean))
}

/// `bool(a op b)`.
pub fn decide(op: RelOp, a: &Expr, b: &Expr) -> BRes<bool> { truthy(&rel_value(op, a, b)?) }

/// SymPy's comparison protocol `is_ge` for operands with bounds (the
/// numeric `_n2` step never applies: bounds are not numbers).
fn protocol_is_ge(l: &Expr, r: &Expr) -> BRes<Option<bool>> {
    if let Some(v) = dispatch_is_ge(l, r)? { return Ok(Some(v)); }
    if yes(l, "extended_real") && yes(r, "extended_real") {
        if (yes(l, "infinite") && yes(l, "extended_positive")) || (yes(r, "infinite") && yes(r, "extended_negative")) {
            return Ok(Some(true));
        }
        let diff = op_sub(l, r)?;
        if !is_nan(&diff) {
            if let Some(v) = fact(&diff, "extended_nonnegative") { return Ok(Some(v)); }
        }
    }
    Ok(None)
}

/// `sympy.core.relational.is_ge`.
pub fn is_ge(l: &Expr, r: &Expr) -> BRes<Option<bool>> {
    if let Some(v) = dispatch_is_ge(l, r)? { return Ok(Some(v)); }
    if is_nan(l) || is_nan(r) { return Ok(None); }
    if !yes(l, "extended_real") || !yes(r, "extended_real") { return Ok(None); }
    let result = rel_value(RelOp::Ge, l, r)?;
    match result.kind() {
        Kind::Bool(v) => Ok(Some(*v)),
        _ => Ok(fact(&op_sub(l, r)?, "extended_nonnegative")),
    }
}

fn is_lt(l: &Expr, r: &Expr) -> BRes<Option<bool>> { Ok(is_ge(l, r)?.map(|v| !v)) }

fn describe(e: &Expr) -> String {
    match hooks() {
        Some(hooks) => (hooks.describe)(e),
        None => format!("<class '{}'>", e.class_name()),
    }
}

/// The `_eval_is_ge` registrations of the bounds.
pub fn dispatch_is_ge(l: &Expr, r: &Expr) -> BRes<Option<bool>> {
    match (is_bounds(l), is_bounds(r)) {
        (true, true) => {
            let ((lmin, lmax), (rmin, rmax)) = (ends(l), ends(r));
            if is_ge(&lmin, &rmax)? == Some(true) { return Ok(Some(true)); }
            if is_lt(&lmax, &rmin)? == Some(true) { return Ok(Some(false)); }
            Ok(None)
        }
        (true, false) => {
            if !yes(r, "extended_real") {
                return type_error(format!("Invalid comparison of {} {}", describe(r), crate::printer::to_str(r)));
            }
            if is_comparable(r) {
                let (lmin, lmax) = ends(l);
                if is_ge(&lmin, r)? == Some(true) { return Ok(Some(true)); }
                if is_lt(&lmax, r)? == Some(true) { return Ok(Some(false)); }
            }
            Ok(None)
        }
        (false, true) => {
            if !yes(l, "extended_real") {
                return type_error(format!("Invalid comparison of {} {}", describe(l), crate::printer::to_str(l)));
            }
            if is_comparable(l) {
                let (rmin, rmax) = ends(r);
                // is_le(rhs.max, lhs) and is_gt(rhs.min, lhs)
                if is_ge(l, &rmax)? == Some(true) { return Ok(Some(true)); }
                if is_ge(l, &rmin)?.map(|v| !v) == Some(true) { return Ok(Some(false)); }
            }
            Ok(None)
        }
        (false, false) => Ok(None),
    }
}

/// Python's builtin `min(items)`.
fn py_min(items: &[Expr]) -> BRes<Expr> {
    let mut best = items[0].clone();
    for item in &items[1..] {
        if decide(RelOp::Lt, item, &best)? { best = item.clone(); }
    }
    Ok(best)
}

/// Python's builtin `max(items)`.
fn py_max(items: &[Expr]) -> BRes<Expr> {
    let mut best = items[0].clone();
    for item in &items[1..] {
        if decide(RelOp::Gt, item, &best)? { best = item.clone(); }
    }
    Ok(best)
}

// --- the bounds class -------------------------------------------------------

/// `e.is_extended_real`, completed for the products, sums and powers whose
/// realness the kernel leaves unknown but the reference derives.
fn extended_real(e: &Expr) -> Option<bool> {
    let value = fact(e, "extended_real");
    let args = python_args(e);
    if value.is_some() || args.is_empty() { return value; }
    if e.is_mul() || e.is_add() {
        return if args.iter().all(|a| extended_real(a) == Some(true)) { Some(true) } else { None };
    }
    if e.is_pow() {
        let (base, exp) = (&args[0], &args[1]);
        if yes(base, "extended_nonnegative") && yes(exp, "extended_real") && finite(base) == Some(true) && finite(exp) == Some(true) {
            return Some(true);
        }
    }
    None
}

/// `AccumBounds(min, max)`.
pub fn ab_new(min: &Expr, max: &Expr) -> BRes<Expr> {
    if extended_real(min) != Some(true) || extended_real(max) != Some(true) {
        return value_error("Only real AccumulationBounds are supported");
    }
    if max == min { return Ok(max.clone()); }
    let bad = if is_number(max) && is_number(min) {
        is_comparable(max) && is_comparable(min) && decide(RelOp::Lt, max, min)?
    } else {
        yes(&op_sub(max, min)?, "extended_negative")
    };
    if bad { return value_error("Lower limit should be smaller than upper limit"); }
    node(min.clone(), max.clone())
}

fn ab(min: Expr, max: Expr) -> BRes<Expr> { ab_new(&min, &max) }
fn whole() -> BRes<Expr> { ab(noo(), oo()) }

/// `AccumBounds.__add__` (and `__radd__`).
pub fn ab_add(s: &Expr, o: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    if is_bounds(o) {
        let (omin, omax) = ends(o);
        return ab(add_ctor(&[min, omin])?, add_ctor(&[max, omax])?);
    }
    if (is_oo(o) && is_noo(&min)) || (is_noo(o) && is_oo(&max)) { return whole(); }
    if yes(o, "extended_real") {
        if is_noo(&min) && is_oo(&max) { return whole(); }
        if is_noo(&min) { return ab(noo(), op_add(&max, o)?); }
        if is_oo(&max) { return ab(op_add(&min, o)?, oo()); }
        return ab(add_ctor(&[min, o.clone()])?, add_ctor(&[max, o.clone()])?);
    }
    Ok(held(true, vec![s.clone(), o.clone()]))
}

/// `AccumBounds.__neg__`.
pub fn ab_neg(s: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    ab(op_neg(&max)?, op_neg(&min)?)
}

/// `AccumBounds.__sub__`.
pub fn ab_sub(s: &Expr, o: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    if is_bounds(o) {
        let (omin, omax) = ends(o);
        return ab(add_ctor(&[min, op_neg(&omax)?])?, add_ctor(&[max, op_neg(&omin)?])?);
    }
    if (is_noo(o) && is_noo(&min)) || (is_oo(o) && is_oo(&max)) { return whole(); }
    if yes(o, "extended_real") {
        if is_noo(&min) && is_oo(&max) { return whole(); }
        if is_noo(&min) { return ab(noo(), op_sub(&max, o)?); }
        if is_oo(&max) { return ab(op_sub(&min, o)?, oo()); }
        let no = op_neg(o)?;
        return ab(add_ctor(&[min, no.clone()])?, add_ctor(&[max, no])?);
    }
    add_ctor(&[s.clone(), op_neg(o)?])
}

/// `AccumBounds.__rsub__`: `o - s`.
pub fn ab_rsub(s: &Expr, o: &Expr) -> BRes<Expr> { op_add(&ab_neg(s)?, o) }

fn is_whole(s: &Expr) -> bool {
    let (min, max) = ends(s);
    is_noo(&min) && is_oo(&max)
}

/// `AccumBounds.__mul__` (and `__rmul__`).
pub fn ab_mul(s: &Expr, o: &Expr) -> BRes<Expr> {
    if is_whole(s) { return Ok(s.clone()); }
    let (min, max) = ends(s);
    if is_bounds(o) {
        if is_whole(o) { return Ok(o.clone()); }
        let mut v: Vec<Expr> = Vec::new();
        for a in [&min, &max] {
            let vi = op_mul(o, a)?;
            let parts = python_args(&vi);
            for p in if parts.is_empty() { vec![vi] } else { parts } {
                if !v.contains(&p) { v.push(p); }
            }
        }
        return ab(lattice(false, &v)?, lattice(true, &v)?);
    }
    if is_oo(o) {
        if zero(&min) == Some(true) { return ab(Expr::zero(), oo()); }
        if zero(&max) == Some(true) { return ab(noo(), Expr::zero()); }
    }
    if is_noo(o) {
        if zero(&min) == Some(true) { return ab(noo(), Expr::zero()); }
        if zero(&max) == Some(true) { return ab(Expr::zero(), oo()); }
    }
    if yes(o, "extended_real") {
        if zero(o) == Some(true) {
            if is_oo(&max) { return ab(Expr::zero(), oo()); }
            if is_noo(&min) { return ab(noo(), Expr::zero()); }
            return Ok(Expr::zero());
        }
        if yes(o, "extended_positive") {
            return ab(mul_ctor(&[min, o.clone()])?, mul_ctor(&[max, o.clone()])?);
        } else if yes(o, "extended_negative") {
            return ab(mul_ctor(&[max, o.clone()])?, mul_ctor(&[min, o.clone()])?);
        }
    }
    if crate::series::order::is_order(o) { return Ok(o.clone()); }
    Ok(held(false, vec![s.clone(), o.clone()]))
}

/// `Mul(*factors)` as `Mul.flatten` folds numbers and bounds.
fn flat_mul(factors: &[Expr]) -> BRes<Expr> {
    let mut coeff = Expr::one();
    let mut rest = Vec::new();
    for factor in factors {
        for term in make_args_mul(factor) {
            if is_bounds(&term) {
                coeff = ab_mul(&term, &coeff)?;
            } else if is_number_class(&term) {
                coeff = op_mul(&coeff, &term)?;
            } else {
                rest.push(term);
            }
        }
    }
    if rest.is_empty() { return Ok(coeff); }
    let product = mul_ctor(&rest)?;
    if !is_bounds(&coeff) { return op_mul(&coeff, &product); }
    let mut v = vec![coeff];
    v.extend(make_args_mul(&product));
    Ok(held(false, v))
}

fn inv(e: &Expr) -> BRes<Expr> { op_div(&Expr::one(), e) }

/// `AccumBounds.__truediv__`.
pub fn ab_div(s: &Expr, o: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    let ext = |e: &Expr, n: &str| yes(e, n);
    if is_bounds(o) {
        let (omin, omax) = ends(o);
        if positive(&omin) == Some(true) || negative(&omax) == Some(true) {
            return op_mul(s, &ab(inv(&omax)?, inv(&omin)?)?);
        }
        if ext(&min, "extended_nonpositive") && ext(&max, "extended_nonnegative")
            && ext(&omin, "extended_nonpositive") && ext(&omax, "extended_nonnegative") {
            if zero(&min) == Some(true) && zero(&omin) == Some(true) { return ab(Expr::zero(), oo()); }
            if zero(&max) == Some(true) && zero(&omin) == Some(true) { return ab(noo(), Expr::zero()); }
            return whole();
        }
        if ext(&max, "extended_negative") {
            if ext(&omin, "extended_negative") {
                if zero(&omax) == Some(true) { return ab(op_div(&max, &omin)?, oo()); }
                if ext(&omax, "extended_positive") { return whole(); }
            }
            if zero(&omin) == Some(true) && ext(&omax, "extended_positive") {
                return ab(noo(), op_div(&max, &omax)?);
            }
        }
        if ext(&min, "extended_positive") {
            if ext(&omin, "extended_negative") {
                if zero(&omax) == Some(true) { return ab(noo(), op_div(&min, &omin)?); }
                if ext(&omax, "extended_positive") { return whole(); }
            }
            if zero(&omin) == Some(true) && ext(&omax, "extended_positive") {
                return ab(op_div(&min, &omax)?, oo());
            }
        }
    } else if ext(o, "extended_real") {
        if is_oo(o) || is_noo(o) {
            if *s == whole()? { return whole(); }
            if is_oo(&max) { return ab(lattice(false, &[Expr::zero(), o.clone()])?, lattice(true, &[Expr::zero(), o.clone()])?); }
            if is_noo(&min) {
                let no = op_neg(o)?;
                return ab(lattice(false, &[Expr::zero(), no.clone()])?, lattice(true, &[Expr::zero(), no])?);
            }
        }
        if ext(o, "extended_positive") { return ab(op_div(&min, o)?, op_div(&max, o)?); }
        if ext(o, "extended_negative") { return ab(op_div(&max, o)?, op_div(&min, o)?); }
    }
    let r = inv(o)?;
    if is_zoo(&r) { return Ok(held(false, vec![s.clone(), r])); }
    flat_mul(&[s.clone(), r])
}

/// `AccumBounds.__rtruediv__`: `o/s`.
pub fn ab_rdiv(s: &Expr, o: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    if yes(o, "extended_real") {
        if zero(o) == Some(true) { return Ok(Expr::zero()); }
        if yes(&min, "extended_nonpositive") && yes(&max, "extended_nonnegative") {
            if zero(&min) == Some(true) {
                if yes(o, "extended_positive") { return ab(mul_ctor(&[o.clone(), inv(&max)?])?, oo()); }
                if yes(o, "extended_negative") { return ab(noo(), mul_ctor(&[o.clone(), inv(&max)?])?); }
            }
            if zero(&max) == Some(true) {
                if yes(o, "extended_positive") { return ab(noo(), mul_ctor(&[o.clone(), inv(&min)?])?); }
                if yes(o, "extended_negative") { return ab(mul_ctor(&[o.clone(), inv(&min)?])?, oo()); }
            }
            return whole();
        }
        let (a, b) = (op_div(o, &min)?, op_div(o, &max)?);
        return ab(lattice(false, &[a.clone(), b.clone()])?, lattice(true, &[a, b])?);
    }
    Ok(held(false, vec![o.clone(), inv(s)?]))
}

/// `AccumBounds.__pow__`.
pub fn ab_pow(s: &Expr, o: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    let one = Expr::one();
    let minus_one = Expr::minus_one();
    if is_oo(o) {
        if yes(&min, "extended_nonnegative") {
            if decide(RelOp::Lt, &max, &one)? { return Ok(Expr::zero()); }
            if decide(RelOp::Gt, &min, &one)? { return Ok(oo()); }
            return ab(Expr::zero(), oo());
        } else if yes(&max, "extended_negative") {
            if decide(RelOp::Gt, &min, &minus_one)? { return Ok(Expr::zero()); }
            if decide(RelOp::Lt, &max, &minus_one)? { return Ok(Expr::zoo()); }
            return Ok(Expr::nan());
        } else {
            if decide(RelOp::Gt, &min, &minus_one)? {
                if decide(RelOp::Lt, &max, &one)? { return Ok(Expr::zero()); }
                return ab(Expr::zero(), oo());
            }
            return whole();
        }
    }
    if is_noo(o) { return op_pow(&inv(s)?, &oo()); }
    if nonnegative(&op_sub(&max, &min)?) == Some(true) && nonnegative(&min) == Some(true) && nonnegative(o) == Some(true) {
        return ab(op_pow(&min, o)?, op_pow(&max, o)?);
    }
    if zero(o) == Some(true) { return Ok(Expr::one()); }
    if o.as_number().is_some_and(|n| n.is_integer()) || integer(o) == Some(true) {
        if yes(&min, "extended_positive") {
            let (a, b) = (op_pow(&min, o)?, op_pow(&max, o)?);
            return ab(lattice(false, &[a.clone(), b.clone()])?, lattice(true, &[a, b])?);
        } else if yes(&max, "extended_negative") {
            let (a, b) = (op_pow(&max, o)?, op_pow(&min, o)?);
            return ab(lattice(false, &[a.clone(), b.clone()])?, lattice(true, &[a, b])?);
        }
        let parity = crate::special::modulo(o, &Expr::int(2));
        if parity == Expr::zero() {
            if yes(o, "extended_negative") {
                if zero(&min) == Some(true) { return ab(op_pow(&max, o)?, oo()); }
                if zero(&max) == Some(true) { return ab(op_pow(&min, o)?, oo()); }
                return op_pow(&inv(s)?, &op_neg(o)?);
            }
            let (a, b) = (op_pow(&min, o)?, op_pow(&max, o)?);
            return ab(Expr::zero(), lattice(true, &[a, b])?);
        } else if parity == Expr::one() {
            if yes(o, "extended_negative") {
                if zero(&min) == Some(true) { return ab(op_pow(&max, o)?, oo()); }
                if zero(&max) == Some(true) { return ab(noo(), op_pow(&min, o)?); }
                return op_pow(&inv(s)?, &op_neg(o)?);
            }
            return ab(op_pow(&min, o)?, op_pow(&max, o)?);
        }
    }
    // non-integer exponent: 0**neg or neg**frac yields complex
    if (is_number(o) || rational(o) == Some(true))
        && (yes(&min, "extended_nonnegative") || (yes(o, "extended_nonnegative") && yes(&min, "extended_nonnegative"))) {
        let (num, den) = crate::simptools::as_numer_denom(o);
        if num == Expr::one() {
            let r = inv(&den)?;
            return ab(op_pow(&min, &r)?, op_pow(&max, &r)?);
        } else if den != Expr::one() {
            return op_pow(&op_pow(s, &num)?, &inv(&den)?);
        }
    }
    if is_bounds(o) {
        let (omin, _) = ends(o);
        if yes(&min, "extended_positive") || (yes(&min, "extended_nonnegative") && yes(&omin, "extended_nonnegative")) {
            let (oa, ob) = ends(o);
            let p = [op_pow(s, &oa)?, op_pow(s, &ob)?];
            if !p.iter().any(|i| i.is_pow()) {
                let mut a = Vec::new();
                for i in &p {
                    let args = python_args(i);
                    if args.is_empty() { a.push(i.clone()); } else { a.extend(args); }
                }
                // try: AccumBounds(min(a), max(a)) except TypeError: pass
                let attempt = (|| -> BRes<Expr> { ab_new(&py_min(&a)?, &py_max(&a)?) })();
                match attempt {
                    Err(BoundsError::Type(_)) => {}
                    other => return other,
                }
            }
        }
    }
    Ok(held_pow(s, o))
}

/// `AccumBounds.__rpow__`: `o**s`.
pub fn ab_rpow(s: &Expr, o: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    if real(o) == Some(true) && yes(o, "extended_nonnegative") && yes(&op_sub(&max, &min)?, "extended_positive") {
        if is_one(o) { return Ok(Expr::one()); }
        if yes(o, "extended_positive") {
            let (mut a, mut b) = (op_pow(o, &min)?, op_pow(o, &max)?);
            // min(a, b) != a
            let smaller = if decide(RelOp::Lt, &b, &a)? { b.clone() } else { a.clone() };
            if smaller != a { std::mem::swap(&mut a, &mut b); }
            return ab_new(&a, &b);
        }
        if zero(o) == Some(true) {
            if zero(&min) == Some(true) { return ab(Expr::zero(), Expr::one()); }
            if yes(&min, "extended_positive") { return Ok(Expr::zero()); }
        }
    }
    Ok(held_pow(o, s))
}

/// `AccumBounds.__abs__`.
pub fn ab_abs(s: &Expr) -> BRes<Expr> {
    let (min, max) = ends(s);
    if yes(&max, "extended_negative") { return ab_neg(s); }
    if yes(&min, "extended_negative") {
        let a = crate::functions::func(Func::Abs, vec![min]);
        return ab(Expr::zero(), lattice(true, &[a, max])?);
    }
    Ok(s.clone())
}

/// `other in AccumBounds(...)`.
pub fn ab_contains(s: &Expr, o: &Expr) -> BRes<bool> {
    let (min, max) = ends(s);
    if is_oo(o) || is_noo(o) { return Ok(is_noo(&min) || is_oo(&max)); }
    let le = rel_value(RelOp::Le, &min, o)?;
    let ge = rel_value(RelOp::Ge, &max, o)?;
    // And(min <= other, max >= other) must evaluate to a truth value.
    let falsity = |e: &Expr| matches!(e.kind(), Kind::Bool(false));
    let truth = |e: &Expr| matches!(e.kind(), Kind::Bool(true));
    if falsity(&le) || falsity(&ge) { return Ok(false); }
    if truth(&le) && truth(&ge) { return Ok(true); }
    type_error("input failed to evaluate")
}

/// The elements of `items` (a FiniteSet, in its order) inside `s`:
/// `s.intersection(FiniteSet(*items))`.
pub fn ab_members(s: &Expr, items: &[Expr]) -> BRes<Vec<Expr>> {
    let mut members = Vec::new();
    for item in items {
        if ab_contains(s, item)? && !members.contains(item) { members.push(item.clone()); }
    }
    Ok(members)
}

/// `AccumBounds.intersection` of two bounds; `Ok(None)` is SymPy's
/// fall-through (the method returns None); the empty set is
/// `Ok(Some(None))`.
pub fn ab_intersection(s: &Expr, o: &Expr) -> BRes<Option<Option<Expr>>> {
    let ((smin, smax), (omin, omax)) = (ends(s), ends(o));
    if decide(RelOp::Lt, &smax, &omin)? || decide(RelOp::Gt, &smin, &omax)? { return Ok(Some(None)); }
    if decide(RelOp::Le, &smin, &omin)? {
        if decide(RelOp::Le, &smax, &omax)? { return Ok(Some(Some(ab_new(&omin, &smax)?))); }
        if decide(RelOp::Gt, &smax, &omax)? { return Ok(Some(Some(o.clone()))); }
    }
    if decide(RelOp::Le, &omin, &smin)? {
        if decide(RelOp::Lt, &omax, &smax)? { return Ok(Some(Some(ab_new(&smin, &omax)?))); }
        if decide(RelOp::Gt, &omax, &smax)? { return Ok(Some(Some(s.clone()))); }
    }
    Ok(None)
}

/// `AccumBounds.union` of two bounds; `None` is SymPy's fall-through.
pub fn ab_union(s: &Expr, o: &Expr) -> BRes<Option<Expr>> {
    let ((smin, smax), (omin, omax)) = (ends(s), ends(o));
    if decide(RelOp::Le, &smin, &omin)? && decide(RelOp::Ge, &smax, &omin)? {
        return Ok(Some(ab_new(&smin, &lattice(true, &[smax, omax])?)?));
    }
    if decide(RelOp::Le, &omin, &smin)? && decide(RelOp::Ge, &omax, &smin)? {
        return Ok(Some(ab_new(&omin, &lattice(true, &[smax, omax])?)?));
    }
    Ok(None)
}

// --- function evaluation ----------------------------------------------------

/// The heads whose `eval` returns accumulation bounds.
pub const HEADS: [&str; 18] = ["sin", "cos", "tan", "cot", "sec", "csc", "exp", "log", "frac",
    "atanh", "asech", "atan", "sinh", "cosh", "tanh", "coth", "sech", "csch"];
const MUL_HEADS: [&str; 7] = ["exp", "sinh", "cosh", "tanh", "coth", "sech", "csch"];
/// The `_make` gate: `exp`, `cos`, `sec` and the hyperbolic functions of a
/// product reach the bounds branches.
const ACCUM_MUL_HEADS: [&str; 9] = ["exp", "cos", "sec", "sinh", "cosh", "tanh", "coth", "sech", "csch"];
const OSCILLATING: [&str; 7] = ["sin", "cos", "tan", "cot", "sec", "csc", "frac"];

fn pi() -> Expr { Expr::pi() }
fn pi_times(n: i64, d: i64) -> Expr { ops::mul(&[pi(), Expr::rational(n, d)]) }

/// The one-argument checks of the `_apply` binding (`py_apply`).
fn validate(f: &Func, arg: &Expr) -> BRes<()> {
    let v = [arg.clone()];
    if matches!(f, Func::Floor | Func::Ceiling) {
        if let Some(Number::MpFloat(value, _)) = arg.as_number() {
            if !value.integer_part_fits() { return Err(BoundsError::Overflow("too many digits in integer".into())); }
        }
    }
    if let Func::Special(s) = f {
        if !s.nargs().contains(&1) {
            return type_error(format!("{} accepts argument counts {:?}, got 1", f.name(), s.nargs()));
        }
        if let Err((type_err, msg)) = crate::special::numtheory::check_args(*s, &v) {
            return if type_err { type_error(msg) } else { value_error(msg) };
        }
        if let Some(msg) = crate::special::orthopoly_eval_error(*s, &v) { return value_error(msg); }
    }
    if let Some(n) = f.arity() {
        if n != 1 { return type_error(format!("{} takes exactly {n} argument(s)", f.name())); }
    }
    if matches!(f, Func::DiracDelta | Func::Heaviside) {
        let imaginary = crate::functions::func(Func::Im, vec![arg.clone()]);
        if crate::functions::imaginary_part_is_zero(&imaginary) == Some(false) {
            return value_error(format!("Function defined only for Real Values. Complex part: {imaginary}  found in {arg} ."));
        }
    }
    Ok(())
}

/// `f(arg)` through the public function head: the bounds branches of its
/// `eval`, then the native evaluation.
pub fn apply(f: Func, arg: &Expr) -> BRes<Expr> {
    let name = f.name().to_string();
    if HEADS.contains(&name.as_str())
        && (is_bounds(arg) || finite(arg) == Some(false) || (arg.is_mul() && ACCUM_MUL_HEADS.contains(&name.as_str()))) {
        if let Some(value) = evaluate(&name, arg)? { return Ok(value); }
    }
    validate(&f, arg)?;
    Ok(crate::functions::func(f, vec![arg.clone()]))
}

fn call(name: &str, arg: &Expr) -> BRes<Expr> {
    let f = Func::from_name(name).ok_or_else(|| BoundsError::NotImplemented(format!("unknown function {name}")))?;
    apply(f, arg)
}

/// `arg.is_extended_real and arg.is_finite is False`, completed for
/// `oo*sign(s)` with a nonzero real s.
fn real_infinity(arg: &Expr) -> bool {
    if yes(arg, "extended_real") && finite(arg) == Some(false) { return true; }
    if arg.is_mul() {
        let (coeff, rest) = crate::decompose::as_coeff_mul_head(arg, false);
        return (is_oo(&coeff) || is_noo(&coeff)) && fact(&rest, "extended_real") == Some(true) && zero(&rest) == Some(false);
    }
    false
}

/// Whether `s` meets one of `points` (`s.intersection(FiniteSet(...))` is
/// not empty); a bounds constructor returning a non-bounds value has no
/// `intersection` method.
fn meets(min: &Expr, max: &Expr, points: &[Expr]) -> BRes<bool> {
    let s = ab_new(min, max)?;
    if !is_bounds(&s) {
        return Err(BoundsError::Attribute(format!("'{}' object has no attribute 'intersection'", s.class_name())));
    }
    Ok(!ab_members(&s, points)?.is_empty())
}

fn eval_sin(arg: &Expr) -> BRes<Option<Expr>> {
    if is_oo(arg) || is_noo(arg) { return ab(Expr::minus_one(), Expr::one()).map(Some); }
    if !is_bounds(arg) { return Ok(None); }
    let (mut min, mut max) = ends(arg);
    let d = call("floor", &op_div(&min, &op_mul(&Expr::int(2), &pi())?)?)?;
    let shift = op_mul(&op_mul(&d, &Expr::int(2))?, &pi())?;
    if !is_noo(&min) { min = op_sub(&min, &shift)?; }
    if !is_oo(&max) { max = op_sub(&max, &shift)?; }
    let peaks = [op_div(&pi(), &Expr::int(2))?, pi_times(5, 2)];
    let first = meets(&min, &max, &peaks)?;
    if first && meets(&min, &max, &[pi_times(3, 2), pi_times(7, 2)])? {
        return ab(Expr::minus_one(), Expr::one()).map(Some);
    }
    let ends = || -> BRes<[Expr; 2]> { Ok([call("sin", &min)?, call("sin", &max)?]) };
    if first { return ab(lattice(false, &ends()?)?, Expr::one()).map(Some); }
    if meets(&min, &max, &[pi_times(3, 2), pi_times(4, 1)])? {
        return ab(Expr::minus_one(), lattice(true, &ends()?)?).map(Some);
    }
    let v = ends()?;
    ab(lattice(false, &v)?, lattice(true, &v)?).map(Some)
}

fn eval_cos(arg: &Expr) -> BRes<Option<Expr>> {
    if is_oo(arg) || is_noo(arg) { return ab(Expr::minus_one(), Expr::one()).map(Some); }
    if is_bounds(arg) { return call("sin", &op_add(arg, &op_div(&pi(), &Expr::int(2))?)?).map(Some); }
    if real_infinity(arg) { return ab(Expr::minus_one(), Expr::one()).map(Some); }
    Ok(None)
}

fn eval_tan(arg: &Expr) -> BRes<Option<Expr>> {
    if is_oo(arg) || is_noo(arg) { return whole().map(Some); }
    if !is_bounds(arg) { return Ok(None); }
    let (mut min, mut max) = ends(arg);
    let d = call("floor", &op_div(&min, &pi())?)?;
    let shift = op_mul(&d, &pi())?;
    if !is_noo(&min) { min = op_sub(&min, &shift)?; }
    if !is_oo(&max) { max = op_sub(&max, &shift)?; }
    if meets(&min, &max, &[op_div(&pi(), &Expr::int(2))?, pi_times(3, 2)])? { return whole().map(Some); }
    ab(call("tan", &min)?, call("tan", &max)?).map(Some)
}

fn eval_cot(arg: &Expr) -> BRes<Option<Expr>> {
    if is_oo(arg) || is_noo(arg) { return whole().map(Some); }
    if is_bounds(arg) { return op_neg(&call("tan", &op_add(arg, &op_div(&pi(), &Expr::int(2))?)?)?).map(Some); }
    Ok(None)
}

fn eval_exp(arg: &Expr) -> BRes<Option<Expr>> {
    if is_bounds(arg) {
        let (min, max) = ends(arg);
        return ab(call("exp", &min)?, call("exp", &max)?).map(Some);
    }
    if arg.is_mul() {
        // exp.eval: an infinite coefficient of a number, exp(I*oo) is nan.
        let (coeff, mut terms) = crate::decompose::as_coeff_mul_head(arg, false);
        if (is_oo(&coeff) || is_noo(&coeff)) && is_number(&terms) {
            if is_noo(&coeff) { terms = op_neg(&terms)?; }
            let re = call("re", &terms)?;
            if zero(&re) == Some(true) && terms != Expr::zero() { return Ok(Some(Expr::nan())); }
            if positive(&re) == Some(true) && call("im", &terms)? != Expr::zero() { return Ok(Some(Expr::zoo())); }
            if negative(&re) == Some(true) { return Ok(Some(Expr::zero())); }
        }
    }
    Ok(None)
}

fn eval_log(arg: &Expr) -> BRes<Option<Expr>> {
    if !is_bounds(arg) { return Ok(None); }
    let (min, max) = ends(arg);
    if positive(&min) == Some(true) { return ab(call("log", &min)?, call("log", &max)?).map(Some); }
    if zero(&min) == Some(true) { return ab(noo(), call("log", &max)?).map(Some); }
    // the reference `log.eval`: any other lower end (negative or of unknown
    // sign) is `nan`.
    Ok(Some(Expr::nan()))
}

fn eval_hyperbolic(name: &str, arg: &Expr) -> BRes<Option<Expr>> {
    if !arg.is_mul() { return Ok(None); }
    let i = Expr::i();
    let c = op_div(arg, &i)?;
    if !(is_oo(&c) || is_noo(&c)) { return Ok(None); }
    Ok(Some(match name {
        "sinh" => op_mul(&i, &call("sin", &c)?)?,
        "cosh" => call("cos", &c)?,
        "tanh" => op_mul(&i, &call("tan", &c)?)?,
        "coth" => op_mul(&op_neg(&i)?, &call("cot", &c)?)?,
        "sech" => inv(&call("cos", &c)?)?,
        // an odd reciprocal: csch(-I*oo) is -csch(I*oo), and 1/(I*sin(oo))
        // keeps I as a factor of the bounds.
        _ if is_oo(&c) => held(false, vec![i, whole()?]),
        _ => held(false, vec![Expr::minus_one(), i, whole()?]),
    }))
}

/// The reference value of `name(arg)` when it is (or involves)
/// accumulation bounds, else `None` (`_accum_eval.evaluate`).
pub fn evaluate(name: &str, arg: &Expr) -> BRes<Option<Expr>> {
    if !HEADS.contains(&name) { return Ok(None); }
    // Only infinite arguments and the bounds themselves reach the branches.
    if !is_bounds(arg) && !yes(arg, "infinite") && finite(arg) != Some(false) {
        let args = python_args(arg);
        if !arg.is_mul() || args.is_empty() || !yes(&args[0], "infinite") { return Ok(None); }
        if !MUL_HEADS.contains(&name) && !(matches!(name, "cos" | "sec") && real_infinity(arg)) { return Ok(None); }
    }
    match name {
        "sin" => eval_sin(arg),
        "cos" => eval_cos(arg),
        "tan" => eval_tan(arg),
        "cot" => eval_cot(arg),
        "sec" => match eval_cos(arg)? { Some(t) => inv(&t).map(Some), None => Ok(None) },
        "csc" => match eval_sin(arg)? { Some(t) => inv(&t).map(Some), None => Ok(None) },
        "exp" => eval_exp(arg),
        "log" => eval_log(arg),
        "frac" => if is_oo(arg) || is_noo(arg) { ab(Expr::zero(), Expr::one()).map(Some) } else { Ok(None) },
        "atanh" | "asech" => {
            if !is_zoo(arg) { return Ok(None); }
            let half = op_div(&pi(), &Expr::int(2))?;
            op_mul(&Expr::i(), &ab(op_neg(&half)?, half)?).map(Some)
        }
        "atan" => {
            if !is_zoo(arg) { return Ok(None); }
            let half = op_div(&pi(), &Expr::int(2))?;
            ab(op_neg(&half)?, half).map(Some)
        }
        _ => eval_hyperbolic(name, arg),
    }
}

// --- the oscillating limit --------------------------------------------------

fn head_name(e: &Expr) -> Option<&str> {
    match e.kind() {
        Kind::Func(f, _) => Some(f.name()),
        _ => None,
    }
}

fn is_oscillating_head(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc | Func::Special(SpecialFn::Frac), _))
}

/// Whether an oscillating function of `z` occurs in `e`.
pub fn oscillates(e: &Expr, z: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if !found && is_oscillating_head(s) && crate::traverse::has(s, z) { found = true; }
    });
    found
}

fn any_node(e: &Expr, p: &dyn Fn(&Expr) -> bool) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| { if !found && p(s) { found = true; } });
    found
}

fn has_nan(e: &Expr) -> bool { any_node(e, &is_nan) }

fn bounded_value(value: Expr) -> BRes<Expr> {
    if has_nan(&value) || any_node(&value, &|s| matches!(s.kind(), Kind::Func(Func::Limit, _))) {
        return Err(BoundsError::Fail);
    }
    Ok(value)
}

/// The native limit as the `_limit` binding computes it, recording its
/// declines for the audits.
fn native_limit(e: &Expr, z: &Expr, z0: &Expr, dir: Dir) -> Option<Expr> {
    let held = super::held_limit(e, z, z0, dir);
    super::validate_held(&held).ok()?;
    let args = held.args();
    let d = Dir::parse(args[3].as_sym()?.name())?;
    match super::limit(&args[0], &args[1], &args[2], d) {
        Ok(r) => Some(r),
        Err(LimitError::Unknown(s)) => {
            super::trace::record(&s, &held);
            None
        }
        Err(LimitError::Differ(..)) => None,
    }
}

fn native(e: &Expr, z: &Expr, z0: &Expr, dir: Dir) -> BRes<Expr> {
    match native_limit(e, z, z0, dir) {
        Some(r) => bounded_value(r),
        None => Err(BoundsError::Fail),
    }
}

fn is_number_like(v: &Expr) -> bool { is_bounds(v) || (is_number(v) && yes(v, "extended_real")) }

/// Numbers and bounds combine with the bounds arithmetic, as the reference
/// flatten does for numeric terms; symbolic parts stay in the constructor.
fn combine(values: Vec<Expr>, add: bool) -> BRes<Expr> {
    let op = |a: &Expr, b: &Expr| if add { op_add(a, b) } else { op_mul(a, b) };
    let rebuild = |v: &[Expr]| if add { add_ctor(v) } else { mul_ctor(v) };
    let (numeric, symbolic): (Vec<Expr>, Vec<Expr>) = values.into_iter().partition(is_number_like);
    if numeric.is_empty() { return rebuild(&symbolic); }
    let mut combined = numeric[0].clone();
    for v in &numeric[1..] { combined = op(&combined, v)?; }
    if symbolic.is_empty() { return Ok(combined); }
    if is_number_like(&combined) && !is_bounds(&combined) {
        return op(&combined, &rebuild(&symbolic)?);
    }
    let mut all = vec![combined];
    all.extend(symbolic);
    rebuild(&all)
}

/// `value` holds an oscillating function left unevaluated at an infinite
/// argument, such as `sin(oo*sign(s))`.
fn held_oscillation(value: &Expr) -> bool {
    any_node(value, &|node| {
        head_name(node).is_some_and(|n| OSCILLATING.contains(&n))
            && python_args(node).first().is_some_and(|a| crate::traverse::has(a, &oo()) || crate::traverse::has(a, &noo()))
    })
}

/// `oo*sign(s)` up to a sign: the reference heuristics value of the
/// argument `s*x` at an infinity.
fn real_infinity_times_sign(arg: &Expr) -> bool {
    crate::sortkey::ordered_factors(arg).iter().all(|factor| {
        is_oo(factor) || is_noo(factor) || *factor == Expr::minus_one() || matches!(factor.kind(), Kind::Func(Func::Sign, _))
    })
}

/// Whether the walker evaluates applications of `f` (the public heads).
fn walkable(f: &Func) -> bool {
    !matches!(f, Func::Undefined(_) | Func::Order | Func::Logic(_) | Func::Set(_) | Func::Predicate(_) | Func::AppliedPredicate
        | Func::GenericExpr | Func::ExpPolar)
}

fn walk(e: &Expr, z: &Expr, z0: &Expr, dir: Dir, used: &mut bool) -> BRes<Expr> {
    if !crate::traverse::has(e, z) { return Ok(e.clone()); }
    if e == z { return Ok(z0.clone()); }
    if !oscillates(e, z) { return native(e, z, z0, dir); }
    if e.is_add() {
        let mut values = Vec::new();
        for a in python_args(e) { values.push(walk(&a, z, z0, dir, used)?); }
        if values.iter().any(is_oo) && values.iter().any(is_noo) { return Err(BoundsError::Fail); }
        return bounded_value(combine(values, true)?);
    }
    if e.is_mul() {
        let mut values = Vec::new();
        for a in python_args(e) { values.push(walk(&a, z, z0, dir, used)?); }
        return bounded_value(combine(values, false)?);
    }
    if let Kind::Pow(eb, ee) = e.kind() {
        let base = walk(eb, z, z0, dir, used)?;
        let exp = walk(ee, z, z0, dir, used)?;
        if is_bounds(&base) && exp == Expr::zero() && crate::traverse::has(ee, z) {
            // The reference leading term keeps the vanishing exponent:
            // AccumBounds(a, b)**t with t -> 0+ has the bounds of a**t and
            // b**t, where a zero endpoint stays 0 (oo for t -> 0-) and an
            // infinite endpoint stays oo (0 for t -> 0-), so frac(x)**(1/x)
            // at oo is AccumBounds(0, 1), not 1.
            let side = native(&op_div(&Expr::one(), ee)?, z, z0, dir)?;
            let (bmin, bmax) = ends(&base);
            if !(is_oo(&side) || is_noo(&side)) || nonnegative(&bmin) != Some(true) { return Err(BoundsError::Fail); }
            let pos = is_oo(&side);
            let endpoint = |v: &Expr| if *v == Expr::zero() {
                if pos { Expr::zero() } else { oo() }
            } else if is_oo(v) {
                if pos { oo() } else { Expr::zero() }
            } else { Expr::one() };
            let points = [endpoint(&bmin), endpoint(&bmax)];
            return bounded_value(ab_new(&py_min(&points)?, &py_max(&points)?)?);
        }
        if is_bounds(&exp) || (is_bounds(&base) && rational(&exp) != Some(true) && nonnegative(&ends(&base).0) != Some(true)) {
            // A base that may be negative under an irrational, infinite or
            // oscillating exponent is not real; the reference keeps such
            // limits unevaluated (issues 14276, 18442).
            return Err(BoundsError::Fail);
        }
        return bounded_value(op_pow(&base, &exp)?);
    }
    if let Kind::Func(f, _) = e.kind() {
        let args = python_args(e);
        let function = !matches!(f, Func::GenericExpr | Func::SymbolicDerivative | Func::Subs | Func::Limit | Func::Integral
            | Func::Sum | Func::Product | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_)
            | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _))
            && !matches!(f, Func::Undefined(head) if head.is_internal("Tuple"));
        if function && args.len() == 1 {
            let name = f.name().to_string();
            let arg = walk(&args[0], z, z0, dir, used)?;
            if !walkable(f) { return Err(BoundsError::Fail); }
            let oscillating = OSCILLATING.contains(&name.as_str());
            if oscillating && !is_bounds(&arg) && fact(&arg, "infinite") != Some(true) && finite(&arg) != Some(false)
                && !real_infinity(&arg) {
                // A finite limit of the argument: the native engine knows
                // better (continuity at the point), unless it was bounds.
                return native(e, z, z0, dir);
            }
            // exp(sin(oo*sign(s))): the reference leaves it unevaluated
            if held_oscillation(&arg) { return Err(BoundsError::Fail); }
            let value = apply(f.clone(), &arg)?;
            if contains_bounds(&value) { *used = true; }
            if oscillating && held_oscillation(&value) {
                // sin(oo*sign(s)) stands for AccumBounds(-1, 1) once the
                // sign of s is known; an argument that may be finite is
                // undecided.
                if finite(&arg) == Some(false) && real_infinity_times_sign(&arg) {
                    *used = true;
                    return Ok(value);
                }
                return Err(BoundsError::Fail);
            }
            return bounded_value(value);
        }
    }
    Err(BoundsError::Fail)
}

fn checked(value: Expr) -> Option<Expr> {
    if has_nan(&value) || is_zoo(&value) { None } else { Some(value) }
}

/// The limit of an expression affine in bounded oscillating parts
/// ([`super::bounded_affine`]), combined with the bounds arithmetic as the
/// reference leading terms combine through `Mul` and `Add` flatten.
pub fn bounded_affine_value(e: &Expr, z: &Expr, z0: &Expr, dir: &str) -> BRes<Option<Expr>> {
    let minus = match dir { "+" => false, "-" => true, _ => return Ok(None) };
    let Some((c0, coeffs)) = super::bounded_affine(e, z, z0, minus) else { return Ok(None) };
    if coeffs.is_empty() { return Ok(checked(c0)); }
    let (c, lo, hi) = &coeffs[0];
    let mut value = mul_ctor(&[c.clone(), ab_new(lo, hi)?])?;
    for (c, lo, hi) in &coeffs[1..] {
        value = op_add(&value, &mul_ctor(&[c.clone(), ab_new(lo, hi)?])?)?;
    }
    for term in make_args_add(&c0) { value = op_add(&value, &term)?; }
    Ok(checked(value))
}

/// The limit of `e` with accumulation bounds (`_accum_limit.oscillating_limit`),
/// or `None`. Errors are the ones SymPy raises from the bounds arithmetic.
pub fn oscillating_limit(e: &Expr, z: &Expr, z0: &Expr, dir: &str) -> BRes<Option<Expr>> {
    if !matches!(dir, "+" | "-" | "+-") || !z.is_symbol() || is_zoo(z0) || is_nan(z0) { return Ok(None); }
    if !registered() { return Ok(None); }
    if !oscillates(e, z) { return bounded_affine_value(e, z, z0, dir); }
    // The native decomposition (affine in bounded sin/cos values, their
    // integer powers and rounding remainders) is tried first; the walker
    // only takes the forms it declines.
    if let Some(value) = bounded_affine_value(e, z, z0, dir)? { return Ok(Some(value)); }
    if dir == "+-" {
        if yes(z0, "infinite") { return Ok(None); }
        let right = oscillating_limit(e, z, z0, "+")?;
        let left = oscillating_limit(e, z, z0, "-")?;
        return Ok(match (right, left) { (Some(r), Some(l)) if r == l => Some(r), _ => None });
    }
    let d = if dir == "+" { Dir::Plus } else { Dir::Minus };
    let mut used = false;
    let value = match walk(e, z, z0, d, &mut used) {
        Ok(value) => value,
        Err(BoundsError::Fail | BoundsError::Type(_) | BoundsError::Value(_) | BoundsError::NotImplemented(_)) => {
            return bounded_affine_value(e, z, z0, dir);
        }
        Err(err) => return Err(err),
    };
    if !used { return Ok(None); }
    Ok(checked(value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::host::{HostClass, HostFormat};
    use std::sync::Arc;

    #[derive(Debug)]
    struct Bounds;

    impl HostClass for Bounds {
        fn identity(&self) -> u64 { 0xAB }
        fn class_name(&self) -> &str { CLASS }
        fn print(&self, _: &Expr, _: HostFormat) -> Option<String> { None }
        fn applied(&self) -> Option<Arc<dyn HostClass>> { Some(Arc::new(Bounds)) }
        fn expression(&self) -> bool { true }
        fn fact(&self, _: &Expr, name: &str) -> Option<Option<bool>> {
            Some(match name {
                "extended_real" | "commutative" => Some(true),
                "complex" | "imaginary" => Some(false),
                _ => None,
            })
        }
        fn as_any(&self) -> &dyn std::any::Any { self }
    }

    fn setup() {
        register(UndefinedFunction::new(CLASS, crate::symbol::Assumptions::NONE, None).with_host(Arc::new(Bounds)));
    }

    fn bounds(lo: Expr, hi: Expr) -> Expr { ab_new(&lo, &hi).unwrap() }

    #[test]
    fn sine_at_infinity() {
        setup();
        let x = Expr::symbol("x");
        let sin = crate::functions::func(Func::Sin, vec![x.clone()]);
        let value = oscillating_limit(&sin, &x, &Expr::oo(), "-").unwrap().unwrap();
        assert!(is_bounds(&value));
        assert_eq!(ends(&value), (Expr::minus_one(), Expr::one()));
        let inner = crate::functions::func(Func::Sin, vec![ops::div(&Expr::one(), &x)]);
        let value = oscillating_limit(&inner, &x, &Expr::zero(), "+").unwrap().unwrap();
        assert_eq!(ends(&value), (Expr::minus_one(), Expr::one()));
    }

    #[test]
    fn arithmetic() {
        setup();
        let b = bounds(Expr::int(1), Expr::int(2));
        assert_eq!(ends(&op_add(&b, &Expr::int(1)).unwrap()), (Expr::int(2), Expr::int(3)));
        assert_eq!(ends(&op_mul(&Expr::int(-2), &b).unwrap()), (Expr::int(-4), Expr::int(-2)));
        let q = op_div(&Expr::one(), &b).unwrap();
        assert_eq!(ends(&q), (Expr::rational(1, 2), Expr::one()));
        let sq = op_pow(&bounds(Expr::int(-1), Expr::int(1)), &Expr::int(2)).unwrap();
        assert_eq!(ends(&sq), (Expr::zero(), Expr::one()));
        assert_eq!(ab_new(&Expr::int(2), &Expr::int(1)), value_error("Lower limit should be smaller than upper limit"));
        assert!(ab_contains(&b, &Expr::rational(3, 2)).unwrap());
        assert!(!ab_contains(&b, &Expr::int(3)).unwrap());
        assert!(decide(RelOp::Lt, &b, &Expr::int(3)).unwrap());
        assert!(decide(RelOp::Lt, &b, &Expr::rational(3, 2)).is_err());
    }
}
