//! Python bindings of the accumulation bounds (`symox::limits::accum`):
//! the class methods of `AccumulationBounds`, the bounds branches of the
//! elementary functions and the oscillating limit.

use std::sync::atomic::{AtomicBool, Ordering};

use pyo3::exceptions::{PyAttributeError, PyNotImplementedError, PyOverflowError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use symox::expr::{Expr, RelOp};
use symox::limits::accum::{self as ab, BRes, BoundsError, Hooks};

use crate::host::PyHostHead;
use crate::{to_expr, PyExpr};

fn w(e: Expr) -> PyExpr { PyExpr { e } }

/// Whether `accumulationbounds._mark_created` ran: the first bounds a
/// native computation returns enable the coefficient folding of the
/// Python constructors, as the reference constructor does.
static MARKED: AtomicBool = AtomicBool::new(false);

fn mark(py: Python<'_>, e: &Expr) -> PyResult<()> {
    if MARKED.load(Ordering::Relaxed) || !ab::contains_bounds(e) { return Ok(()); }
    py.import("symox.calculus.accumulationbounds")?.getattr("_mark_created")?.call0()?;
    MARKED.store(true, Ordering::Relaxed);
    Ok(())
}

fn error(py: Python<'_>, err: BoundsError) -> PyErr {
    match err {
        BoundsError::Value(m) => PyValueError::new_err(m),
        BoundsError::Type(m) => PyTypeError::new_err(m),
        BoundsError::NotImplemented(m) => PyNotImplementedError::new_err(m),
        BoundsError::Overflow(m) => PyOverflowError::new_err(m),
        BoundsError::Attribute(m) => PyAttributeError::new_err(m),
        BoundsError::Host => match crate::host::raise_pending() {
            Err(e) => e,
            Ok(()) => PyValueError::new_err("host callback failed"),
        },
        BoundsError::Fail => {
            let _ = py;
            PyNotImplementedError::new_err("undecided accumulation bounds")
        }
    }
}

/// The Python result of a bounds computation.
fn finish<T>(py: Python<'_>, value: BRes<T>, exprs: impl Fn(&T) -> Vec<&Expr>) -> PyResult<T> {
    let value = match value {
        Ok(v) => v,
        Err(err) => return Err(error(py, err)),
    };
    crate::host::raise_pending()?;
    for e in exprs(&value) { mark(py, e)?; }
    Ok(value)
}

fn expr_result(py: Python<'_>, value: BRes<Expr>) -> PyResult<PyExpr> {
    finish(py, value, |e| vec![e]).map(w)
}

fn opt_result(py: Python<'_>, value: BRes<Option<Expr>>) -> PyResult<Option<PyExpr>> {
    finish(py, value, |e| e.iter().collect()).map(|e| e.map(w))
}

fn bounds_error(py: Python<'_>, err: PyErr) -> BoundsError {
    let message = || err.value(py).str().map(|s| s.to_string()).unwrap_or_default();
    if err.is_instance_of::<PyValueError>(py) { return BoundsError::Value(message()); }
    if err.is_instance_of::<PyTypeError>(py) { return BoundsError::Type(message()); }
    if err.is_instance_of::<PyNotImplementedError>(py) { return BoundsError::NotImplemented(message()); }
    if err.is_instance_of::<PyOverflowError>(py) { return BoundsError::Overflow(message()); }
    if err.is_instance_of::<PyAttributeError>(py) { return BoundsError::Attribute(message()); }
    crate::host::report(py, err, None);
    BoundsError::Host
}

/// SymPy's `Min`/`Max` lattice algorithm of the facade for the arguments
/// the native lattice functions do not take.
fn lattice(max: bool, args: &[Expr]) -> BRes<Expr> {
    Python::with_gil(|py| {
        let run = || -> PyResult<Expr> {
            let class = py.import("symox")?.getattr(if max { "Max" } else { "Min" })?;
            let objects = args.iter().map(|a| crate::expr_object(py, a.clone())).collect::<PyResult<Vec<_>>>()?;
            let value = class.call1(PyTuple::new(py, objects)?)?;
            to_expr(&value)
        };
        run().map_err(|err| bounds_error(py, err))
    })
}

fn describe(e: &Expr) -> String {
    Python::with_gil(|py| {
        crate::expr_object(py, e.clone())
            .and_then(|o| o.bind(py).get_type().repr().map(|r| r.to_string()))
            .unwrap_or_else(|_| format!("<class '{}'>", e.class_name()))
    })
}

/// Register the `AccumulationBounds` host head.
#[pyfunction]
fn _register_accum_bounds(head: &Bound<'_, PyHostHead>) {
    ab::register(head.get().head().clone());
    ab::set_hooks(Hooks { lattice, describe });
}

/// `AccumBounds(min, max)` (validated; equal ends give the end).
#[pyfunction]
fn _accum_new(py: Python<'_>, min: &Bound<'_, PyAny>, max: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let (min, max) = (to_expr(min)?, to_expr(max)?);
    expr_result(py, ab::ab_new(&min, &max))
}

/// The bounds operator `op` with `this` the bounds operand: `add`, `sub`,
/// `rsub`, `mul`, `truediv`, `rtruediv`, `pow`, `rpow`.
#[pyfunction]
fn _accum_op(py: Python<'_>, op: &str, this: &Bound<'_, PyAny>, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let (s, o) = (to_expr(this)?, to_expr(other)?);
    if !ab::is_bounds(&s) { return Err(PyTypeError::new_err("the operand is not AccumulationBounds")); }
    let value = match op {
        "add" => ab::ab_add(&s, &o),
        "sub" => ab::ab_sub(&s, &o),
        "rsub" => ab::ab_rsub(&s, &o),
        "mul" => ab::ab_mul(&s, &o),
        "truediv" => ab::ab_div(&s, &o),
        "rtruediv" => ab::ab_rdiv(&s, &o),
        "pow" => ab::ab_pow(&s, &o),
        "rpow" => ab::ab_rpow(&s, &o),
        _ => return Err(PyValueError::new_err(format!("unknown bounds operator {op}"))),
    };
    expr_result(py, value)
}

/// `-AccumBounds(...)` (`neg`) or `abs(AccumBounds(...))` (`abs`).
#[pyfunction]
fn _accum_unary(py: Python<'_>, op: &str, this: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let s = to_expr(this)?;
    if !ab::is_bounds(&s) { return Err(PyTypeError::new_err("the operand is not AccumulationBounds")); }
    let value = match op {
        "neg" => ab::ab_neg(&s),
        "abs" => ab::ab_abs(&s),
        _ => return Err(PyValueError::new_err(format!("unknown bounds operator {op}"))),
    };
    expr_result(py, value)
}

/// `other in AccumBounds(...)`.
#[pyfunction]
fn _accum_contains(py: Python<'_>, this: &Bound<'_, PyAny>, other: &Bound<'_, PyAny>) -> PyResult<bool> {
    let (s, o) = (to_expr(this)?, to_expr(other)?);
    finish(py, ab::ab_contains(&s, &o), |_| vec![])
}

/// The elements of a FiniteSet inside the bounds, in order.
#[pyfunction]
fn _accum_members(py: Python<'_>, this: &Bound<'_, PyAny>, items: Vec<Bound<'_, PyAny>>) -> PyResult<Vec<PyExpr>> {
    let s = to_expr(this)?;
    let items = items.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    finish(py, ab::ab_members(&s, &items), |_| vec![]).map(|v| v.into_iter().map(w).collect())
}

/// The intersection of two bounds: `(True, bounds)`, `(True, None)` for the
/// empty set, `(False, None)` when the reference method falls through.
#[pyfunction]
fn _accum_intersection(py: Python<'_>, this: &Bound<'_, PyAny>, other: &Bound<'_, PyAny>) -> PyResult<(bool, Option<PyExpr>)> {
    let (s, o) = (to_expr(this)?, to_expr(other)?);
    let value = finish(py, ab::ab_intersection(&s, &o), |v| v.iter().flatten().collect())?;
    Ok(match value {
        Some(v) => (true, v.map(w)),
        None => (false, None),
    })
}

/// The union of two bounds; `None` when the reference method falls through.
#[pyfunction]
fn _accum_union(py: Python<'_>, this: &Bound<'_, PyAny>, other: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    let (s, o) = (to_expr(this)?, to_expr(other)?);
    opt_result(py, ab::ab_union(&s, &o))
}

/// The `_eval_is_ge` registrations of the bounds.
#[pyfunction]
fn _accum_is_ge(py: Python<'_>, lhs: &Bound<'_, PyAny>, rhs: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    let (l, r) = (to_expr(lhs)?, to_expr(rhs)?);
    finish(py, ab::dispatch_is_ge(&l, &r), |_| vec![])
}

/// The coefficient folding of `Add` (`add`) or `Mul` for arguments holding
/// bounds; `None` when it does not apply.
#[pyfunction]
fn _accum_fold(py: Python<'_>, op: &str, args: Vec<Bound<'_, PyAny>>) -> PyResult<Option<PyExpr>> {
    let args = args.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    opt_result(py, ab::fold(op == "Add", &args))
}

/// The bounds value of `name(arg)` (the bounds branches of `eval`), or
/// `None`.
#[pyfunction]
fn _accum_evaluate(py: Python<'_>, name: &str, arg: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    let arg = to_expr(arg)?;
    opt_result(py, ab::evaluate(name, &arg))
}

/// The limit of `e` with accumulation bounds, or `None` when undecided.
#[pyfunction]
fn _oscillating_limit(py: Python<'_>, e: &Bound<'_, PyAny>, z: &Bound<'_, PyAny>, z0: &Bound<'_, PyAny>, dir: &str) -> PyResult<Option<PyExpr>> {
    let (e, z, z0) = (to_expr(e)?, to_expr(z)?, to_expr(z0)?);
    opt_result(py, ab::oscillating_limit(&e, &z, &z0, dir))
}

/// Whether `e` oscillates in `z` (an oscillating function of `z` occurs).
#[pyfunction]
fn _accum_oscillates(e: &Bound<'_, PyAny>, z: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(ab::oscillates(&to_expr(e)?, &to_expr(z)?))
}

/// `bool(a op b)` with the bounds comparison protocol.
#[pyfunction]
fn _accum_decide(py: Python<'_>, op: &str, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<bool> {
    let (a, b) = (to_expr(a)?, to_expr(b)?);
    let op = match op {
        "<" => RelOp::Lt,
        "<=" => RelOp::Le,
        ">" => RelOp::Gt,
        ">=" => RelOp::Ge,
        _ => return Err(PyValueError::new_err(format!("unknown comparison {op}"))),
    };
    finish(py, ab::decide(op, &a, &b), |_| vec![])
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    for f in [
        wrap_pyfunction!(_register_accum_bounds, m)?,
        wrap_pyfunction!(_accum_new, m)?,
        wrap_pyfunction!(_accum_op, m)?,
        wrap_pyfunction!(_accum_unary, m)?,
        wrap_pyfunction!(_accum_contains, m)?,
        wrap_pyfunction!(_accum_members, m)?,
        wrap_pyfunction!(_accum_intersection, m)?,
        wrap_pyfunction!(_accum_union, m)?,
        wrap_pyfunction!(_accum_is_ge, m)?,
        wrap_pyfunction!(_accum_fold, m)?,
        wrap_pyfunction!(_accum_evaluate, m)?,
        wrap_pyfunction!(_oscillating_limit, m)?,
        wrap_pyfunction!(_accum_oscillates, m)?,
        wrap_pyfunction!(_accum_decide, m)?,
    ] {
        m.add_function(f)?;
    }
    Ok(())
}
