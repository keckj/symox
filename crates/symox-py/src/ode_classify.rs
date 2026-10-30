//! Python bindings for the native first-order `classify_ode` matches.

use num_bigint::BigInt;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::ode::classify::{self, Hooks, Outcome};
use symox::Expr;

use crate::{to_expr, wrap as w};

fn outcome<'py, T>(py: Python<'py>, o: &Outcome<T>, parts: impl Fn(&T) -> Vec<Expr>) -> PyResult<PyObject> {
    Ok(match o {
        Outcome::Match(v) => PyTuple::new(py, parts(v).into_iter().map(w))?.into_any().unbind(),
        Outcome::NoMatch => false.into_pyobject(py)?.to_owned().into_any().unbind(),
        Outcome::Decline => py.None(),
    })
}

/// First-order matches of `eq_expanded` in `fx = f(x)`, `df = f(x).diff(x)`.
/// `cancel`, `expand` and `count_ops` are the public callables the
/// reference matcher uses. Each entry of the returned dict is a tuple of the
/// wild values (a match), `False` (no match) or `None` (not decided).
#[pyfunction]
fn ode_first_order_matches(
    py: Python<'_>,
    eq_expanded: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    df: &Bound<'_, PyAny>,
    cancel: &Bound<'_, PyAny>,
    expand: &Bound<'_, PyAny>,
    count_ops: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let eq = to_expr(eq_expanded)?;
    let (fx, x, df) = (to_expr(fx)?, to_expr(x)?, to_expr(df)?);
    let call = |f: &Bound<'_, PyAny>, e: &Expr| -> Option<Expr> { f.call1((w(e.clone()),)).ok().and_then(|v| to_expr(&v).ok()) };
    let mut cancel_hook = |e: &Expr| call(cancel, e);
    let mut expand_hook = |e: &Expr| call(expand, e);
    // A failing count makes every guard reject, which declines.
    let mut count_hook = |e: &Expr| -> BigInt {
        count_ops.call1((w(e.clone()),)).and_then(|v| v.extract::<BigInt>()).unwrap_or_else(|_| BigInt::from(u64::MAX))
    };
    let mut hooks = Hooks { cancel: &mut cancel_hook, expand: &mut expand_hook, count_ops: &mut count_hook };
    let r = classify::first_order(&eq, &fx, &x, &df, &mut hooks);
    let out = PyDict::new(py);
    out.set_item("de", outcome(py, &r.de, |(d, e)| vec![d.clone(), e.clone()])?)?;
    out.set_item("linear", outcome(py, &r.linear, |(p, q)| vec![p.clone(), q.clone()])?)?;
    out.set_item("bernoulli", outcome(py, &r.bernoulli, |(p, q, n)| [p.clone(), q.clone()].into_iter().chain(n.clone()).collect())?)?;
    out.set_item("riccati_special", outcome(py, &r.riccati_special, |_| Vec::new())?)?;
    Ok(out.into_any().unbind())
}

/// `d + e*df` on an already collected expression (`d` excluding the entries
/// of `d_exclude`), with the same result convention.
#[pyfunction]
fn ode_split_de(py: Python<'_>, collected: &Bound<'_, PyAny>, df: &Bound<'_, PyAny>, d_exclude: Vec<Bound<'_, PyAny>>, count_ops: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let c = to_expr(collected)?;
    let df = to_expr(df)?;
    let exclude = d_exclude.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    let mut none = |_: &Expr| -> Option<Expr> { None };
    let mut none2 = |_: &Expr| -> Option<Expr> { None };
    let mut count_hook = |e: &Expr| -> BigInt {
        count_ops.call1((w(e.clone()),)).and_then(|v| v.extract::<BigInt>()).unwrap_or_else(|_| BigInt::from(u64::MAX))
    };
    let mut hooks = Hooks { cancel: &mut none, expand: &mut none2, count_ops: &mut count_hook };
    let r = classify::split_de(&c, &df, &exclude, &mut hooks);
    outcome(py, &r, |(d, e)| vec![d.clone(), e.clone()])
}

/// The ordered `classify_ode` hint list from the matched solver hints in
/// matching order: the `_Integral` entries are added and the result follows
/// `allhints`.
#[pyfunction]
fn ode_ordered_hints(py: Python<'_>, matched: Vec<String>) -> PyResult<PyObject> {
    let registered = classify::registered_hints(&matched);
    let ordered = classify::ordered_hints(&registered);
    Ok(PyList::new(py, ordered)?.into_any().unbind())
}

/// The reference `match_riccati(eq, f, x)`: `(b0, b1, b2)`, `False` or
/// `None` (not decided).
#[pyfunction]
fn ode_match_riccati(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    df: &Bound<'_, PyAny>,
    cancel: &Bound<'_, PyAny>,
    expand: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let eq = to_expr(eq)?;
    let (fx, x, df) = (to_expr(fx)?, to_expr(x)?, to_expr(df)?);
    let call = |f: &Bound<'_, PyAny>, e: &Expr| -> Option<Expr> { f.call1((w(e.clone()),)).ok().and_then(|v| to_expr(&v).ok()) };
    let mut cancel_hook = |e: &Expr| call(cancel, e);
    let mut expand_hook = |e: &Expr| call(expand, e);
    let mut count_hook = |_: &Expr| -> BigInt { BigInt::from(u64::MAX) };
    let mut hooks = Hooks { cancel: &mut cancel_hook, expand: &mut expand_hook, count_ops: &mut count_hook };
    let r = classify::riccati(&eq, &fx, &x, &df, &mut hooks);
    outcome(py, &r, |(a, b, c)| vec![a.clone(), b.clone(), c.clone()])
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ode_first_order_matches, m)?)?;
    m.add_function(wrap_pyfunction!(ode_split_de, m)?)?;
    m.add_function(wrap_pyfunction!(ode_ordered_hints, m)?)?;
    m.add_function(wrap_pyfunction!(ode_match_riccati, m)?)?;
    Ok(())
}
