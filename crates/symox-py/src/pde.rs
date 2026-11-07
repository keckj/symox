//! Bindings of the native port of `sympy/solvers/pde.py`
//! (`symox::ode::pde`).
//!
//! The public callables `simplify`, `expand`, `integrate`, `solve`, `doit`
//! and `dsolve` arrive in a dict (the `lie_group` hooks); an error they raise
//! propagates unchanged. `None` declines: the caller keeps the reference
//! body.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use symox::ode::lie_group::{Exc, LieError};
use symox::ode::pde::{self, Checked, Classified, Pdsolved, Sol};
use symox::Expr;

use crate::ode_lie::{settle, LieHooks};
use crate::{to_expr, wrap as w};

fn exprs(v: &[Bound<'_, PyAny>]) -> Option<Vec<Expr>> {
    v.iter().map(|e| to_expr(e).ok()).collect()
}

fn opt_expr(v: Option<&Bound<'_, PyAny>>) -> Result<Option<Expr>, ()> {
    match v {
        None => Ok(None),
        Some(v) if v.is_none() => Ok(None),
        Some(v) => to_expr(v).map(Some).map_err(|_| ()),
    }
}

/// The `classify_pde` result: the dict (`dict=True`) or the tuple of hints.
fn classified(py: Python<'_>, c: Classified, dict: bool) -> PyResult<PyObject> {
    let hints = PyTuple::new(py, c.hints.iter().copied())?;
    if c.tuple_only || !dict {
        return Ok(hints.into_any().unbind());
    }
    let out = PyDict::new(py);
    out.set_item("order", c.order)?;
    if c.order == 0 {
        out.set_item("default", py.None())?;
        return Ok(out.into_any().unbind());
    }
    if let Some(m) = &c.matched {
        let r = PyDict::new(py);
        for (k, v) in &m.entries {
            r.set_item(w(k.clone()), w(v.clone()))?;
        }
        r.set_item("b", w(m.b.clone()))?;
        r.set_item("c", w(m.c.clone()))?;
        r.set_item("d", w(m.d.clone()))?;
        if let Some(e) = &m.e {
            r.set_item("e", w(e.clone()))?;
        }
        for h in &c.hints {
            out.set_item(*h, &r)?;
        }
    }
    out.set_item("default", c.hints.first().copied())?;
    out.set_item("ordered_hints", hints)?;
    Ok(out.into_any().unbind())
}

/// `classify_pde(eq, func, dict, prep=prep)`; `None` declines.
#[pyfunction]
#[pyo3(signature = (eq, func, dict, prep, hooks))]
fn pde_classify(py: Python<'_>, eq: &Bound<'_, PyAny>, func: Option<&Bound<'_, PyAny>>, dict: bool, prep: bool, hooks: &Bound<'_, PyDict>) -> PyResult<PyObject> {
    let (Ok(eq), Ok(func)) = (to_expr(eq), opt_expr(func)) else { return Ok(py.None()) };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = pde::classify_pde(&eq, func.as_ref(), prep, &mut hk);
    match settle(&hk, r)? {
        Some(c) => classified(py, c, dict),
        None => Ok(py.None()),
    }
}

fn nie_instance(py: Python<'_>, e: LieError<PyErr>) -> PyObject {
    match e {
        LieError::Hook(e) => e.value(py).clone().into_any().unbind(),
        LieError::Raise(_, msg) => PyNotImplementedError::new_err(msg).value(py).clone().into_any().unbind(),
        _ => PyNotImplementedError::new_err("").value(py).clone().into_any().unbind(),
    }
}

/// `pdsolve(eq, func, hint, solvefun, prep=prep)`: `('one', solution)`,
/// `('all', dict)` or `('unclassified', message)` (the classification
/// `NotImplementedError`); `None` declines. `solvefun` is a sample
/// application of the arbitrary function, or `None` for the default one.
#[pyfunction]
#[pyo3(signature = (eq, func, hint, prep, solvefun, hooks))]
fn pde_pdsolve(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    func: Option<&Bound<'_, PyAny>>,
    hint: &str,
    prep: bool,
    solvefun: Option<&Bound<'_, PyAny>>,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let (Ok(eq), Ok(func), Ok(solvefun)) = (to_expr(eq), opt_expr(func), opt_expr(solvefun)) else { return Ok(py.None()) };
    let solvefun = solvefun.unwrap_or_else(|| pde::default_solvefun(&eq, func.as_ref().unwrap_or(&eq)));
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = pde::pdsolve(&eq, func.as_ref(), hint, prep, &solvefun, &mut hk);
    let Some(v) = settle(&hk, r)? else { return Ok(py.None()) };
    let pair = |tag: &str, v: PyObject| -> PyResult<PyObject> { Ok(PyTuple::new(py, [tag.into_pyobject(py)?.into_any().unbind(), v])?.into_any().unbind()) };
    match v {
        Pdsolved::One(e) => pair("one", w(e).into_pyobject(py)?.into_any().unbind()),
        Pdsolved::Unclassified(m) => pair("unclassified", m.into_pyobject(py)?.into_any().unbind()),
        Pdsolved::All { order, default, entries } => {
            let out = PyDict::new(py);
            out.set_item("order", order)?;
            out.set_item("default", default)?;
            let mut failed = Vec::new();
            for (h, r) in entries {
                match r {
                    Ok(e) => out.set_item(h, w(e))?,
                    Err(e) => failed.push((h, nie_instance(py, e))),
                }
            }
            for (h, e) in failed {
                out.set_item(h, e)?;
            }
            pair("all", out.into_any().unbind())
        }
    }
}

/// `pde_<hint>(eq, func, order, match, solvefun)` with the matched values
/// `b`, `c`, `d`, `e`; `None` declines.
#[pyfunction]
#[pyo3(signature = (hint, eq, func, b, c, d, e, solvefun, hooks))]
#[allow(clippy::too_many_arguments)]
fn pde_solve_hint(
    py: Python<'_>,
    hint: &str,
    eq: &Bound<'_, PyAny>,
    func: &Bound<'_, PyAny>,
    b: &Bound<'_, PyAny>,
    c: &Bound<'_, PyAny>,
    d: &Bound<'_, PyAny>,
    e: Option<&Bound<'_, PyAny>>,
    solvefun: &Bound<'_, PyAny>,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let Some(v) = exprs(&[eq.clone(), func.clone(), b.clone(), c.clone(), d.clone(), solvefun.clone()]) else { return Ok(py.None()) };
    let Ok(e) = opt_expr(e) else { return Ok(py.None()) };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = pde::solve_hint(hint, &v[0], &v[1], &v[2], &v[3], &v[4], e.as_ref(), &v[5], &mut hk);
    match settle(&hk, r)? {
        Some(e) => Ok(w(e).into_pyobject(py)?.into_any().unbind()),
        None => Ok(py.None()),
    }
}

fn checked(py: Python<'_>, c: Checked) -> PyResult<PyObject> {
    match c {
        Checked::One(ok, e) => Ok(PyTuple::new(py, [ok.into_pyobject(py)?.to_owned().into_any().unbind(), w(e).into_pyobject(py)?.into_any().unbind()])?.into_any().unbind()),
        Checked::List(v) => {
            let out = PyList::empty(py);
            for c in v {
                out.append(checked(py, c)?)?;
            }
            Ok(out.into_any().unbind())
        }
    }
}

/// `checkpdesol(pde, sol, func, solve_for_func)`: `sol` is one expression,
/// or a list when `seq` (the result is then a list, which the caller turns
/// into the type of the given sequence); `None` declines.
#[pyfunction]
#[pyo3(signature = (pde, sol, seq, func, solve_for_func, hooks))]
fn pde_checkpdesol(
    py: Python<'_>,
    pde: &Bound<'_, PyAny>,
    sol: &Bound<'_, PyAny>,
    seq: bool,
    func: Option<&Bound<'_, PyAny>>,
    solve_for_func: bool,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let (Ok(pde_e), Ok(func)) = (to_expr(pde), opt_expr(func)) else { return Ok(py.None()) };
    let sol = if seq {
        let Ok(list) = sol.downcast::<PyList>() else { return Ok(py.None()) };
        let items: Vec<_> = list.iter().collect();
        let Some(v) = exprs(&items) else { return Ok(py.None()) };
        Sol::Seq(v)
    } else {
        let Ok(e) = to_expr(sol) else { return Ok(py.None()) };
        Sol::One(e)
    };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = pde::checkpdesol(&pde_e, &sol, func.as_ref(), solve_for_func, &mut hk);
    let r = match r {
        Err(LieError::Raise(Exc::NotImplemented, msg)) if msg.starts_with(pde::UNABLE_PREFIX) && !hk.foreign => {
            let text = py.import("symox.utilities.misc")?.getattr("filldedent")?.call1((msg,))?;
            return Err(PyNotImplementedError::new_err(text.unbind()));
        }
        r => r,
    };
    match settle(&hk, r)? {
        Some(c) => checked(py, c),
        None => Ok(py.None()),
    }
}

/// `pde_separate(eq, fun, sep, strategy)`: `(result,)` with `result` the
/// list `[lhs, rhs]` or `None`; `None` declines.
#[pyfunction]
#[pyo3(signature = (eq, fun, sep, strategy, hooks))]
fn pde_separate(py: Python<'_>, eq: &Bound<'_, PyAny>, fun: &Bound<'_, PyAny>, sep: Vec<Bound<'_, PyAny>>, strategy: &str, hooks: &Bound<'_, PyDict>) -> PyResult<PyObject> {
    let (Ok(eq), Ok(fun), Some(sep)) = (to_expr(eq), to_expr(fun), exprs(&sep)) else { return Ok(py.None()) };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = pde::pde_separate(&eq, &fun, &sep, strategy, &mut hk);
    match settle(&hk, r)? {
        Some(Some([l, r])) => Ok(PyTuple::new(py, [PyList::new(py, [w(l), w(r)])?.into_any().unbind()])?.into_any().unbind()),
        Some(None) => Ok(PyTuple::new(py, [py.None()])?.into_any().unbind()),
        None => Ok(py.None()),
    }
}

/// `ode_order(expr, func)` of `deutils`.
#[pyfunction]
fn pde_ode_order(e: &Bound<'_, PyAny>, func: &Bound<'_, PyAny>) -> PyResult<u32> {
    Ok(pde::ode_order(&to_expr(e)?, &to_expr(func)?))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pde_classify, m)?)?;
    m.add_function(wrap_pyfunction!(pde_pdsolve, m)?)?;
    m.add_function(wrap_pyfunction!(pde_solve_hint, m)?)?;
    m.add_function(wrap_pyfunction!(pde_checkpdesol, m)?)?;
    m.add_function(wrap_pyfunction!(pde_separate, m)?)?;
    m.add_function(wrap_pyfunction!(pde_ode_order, m)?)?;
    Ok(())
}
