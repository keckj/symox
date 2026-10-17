//! Bindings for `sympy.solvers.recurr`: `rsolve` and the Abramov/Petkovsek
//! solvers `rsolve_poly`, `rsolve_ratio` and `rsolve_hyper`, plus
//! `difference_delta`.
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use symox::recurr::{self, RsolveError, Solution};
use symox::Expr;
use crate::{to_expr, wrap};

fn coefficients(coeffs: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
    coeffs.try_iter()?.map(|c| to_expr(&c?)).collect()
}

fn wants_symbols(hints: Option<&Bound<'_, PyDict>>) -> PyResult<bool> {
    match hints.map(|h| h.get_item("symbols")).transpose()?.flatten() {
        Some(v) => v.is_truthy(),
        None => Ok(false),
    }
}

fn solution(py: Python<'_>, s: Option<Solution>, symbols: bool) -> PyResult<PyObject> {
    Ok(match s {
        None => py.None(),
        Some((e, c)) if symbols => {
            let c = PyList::new(py, c.into_iter().map(wrap))?;
            PyTuple::new(py, [wrap(e).into_pyobject(py)?.into_any(), c.into_any()])?.into_any().unbind()
        }
        Some((e, _)) => wrap(e).into_pyobject(py)?.into_any().unbind(),
    })
}

#[pyfunction]
#[pyo3(signature = (coeffs, f, n, shift=0, **hints))]
fn rsolve_poly(py: Python<'_>, coeffs: &Bound<'_, PyAny>, f: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>, shift: usize, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (coeffs, f, n) = (coefficients(coeffs)?, to_expr(f)?, to_expr(n)?);
    solution(py, recurr::rsolve_poly(&coeffs, &f, &n, shift), wants_symbols(hints)?)
}

#[pyfunction]
#[pyo3(signature = (coeffs, f, n, **hints))]
fn rsolve_ratio(py: Python<'_>, coeffs: &Bound<'_, PyAny>, f: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (coeffs, f, n) = (coefficients(coeffs)?, to_expr(f)?, to_expr(n)?);
    solution(py, recurr::rsolve_ratio(&coeffs, &f, &n), wants_symbols(hints)?)
}

#[pyfunction]
#[pyo3(signature = (coeffs, f, n, **hints))]
fn rsolve_hyper(py: Python<'_>, coeffs: &Bound<'_, PyAny>, f: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (coeffs, f, n) = (coefficients(coeffs)?, to_expr(f)?, to_expr(n)?);
    solution(py, recurr::rsolve_hyper(&coeffs, &f, &n), wants_symbols(hints)?)
}

/// `difference_delta(expr, n=None, step=1)`: `expr(n + step) - expr(n)`.
#[pyfunction]
#[pyo3(signature = (expr, n=None, step=None))]
fn difference_delta(py: Python<'_>, expr: &Bound<'_, PyAny>, n: Option<&Bound<'_, PyAny>>, step: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
    let n = n.filter(|n| !n.is_none()).map(to_expr).transpose()?;
    let step = step.map(to_expr).transpose()?.unwrap_or_else(Expr::one);
    let e = recurr::difference_delta(&to_expr(expr)?, n.as_ref(), &step).map_err(PyValueError::new_err)?;
    Ok(wrap(e).into_pyobject(py)?.into_any().unbind())
}

/// `hypersimp(f, k)`: the simplified term ratio `f(k + 1)/f(k)` of a
/// hypergeometric term, `None` otherwise.
#[pyfunction]
fn hypersimp(py: Python<'_>, f: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    match recurr::hypersimp(&to_expr(f)?, &to_expr(k)?) {
        Some(r) => Ok(wrap(symox::simplify::simplify(&r)).into_pyobject(py)?.into_any().unbind()),
        None => Ok(py.None()),
    }
}

/// `hypersimilar(f, g, k)`: whether `f/g` is a rational function of `k`.
#[pyfunction]
fn hypersimilar(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(recurr::hypersimilar(&to_expr(f)?, &to_expr(g)?, &to_expr(k)?))
}

/// `rsolve(f, y(n), init=None)`; `init` is a dict keyed by indices or
/// terms `y(i)`, or a list of the values `y(0), y(1), ...`.
#[pyfunction]
#[pyo3(signature = (f, y, init=None))]
fn rsolve(py: Python<'_>, f: &Bound<'_, PyAny>, y: &Bound<'_, PyAny>, init: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
    let (f, y) = (to_expr(f)?, to_expr(y)?);
    let init: Option<Vec<(Expr, Expr)>> = match init {
        None => None,
        Some(i) if i.is_none() => None,
        Some(i) => Some(if let Ok(d) = i.downcast::<PyDict>() {
            d.iter().map(|(k, v)| Ok((to_expr(&k)?, to_expr(&v)?))).collect::<PyResult<_>>()?
        } else if let Ok(l) = i.downcast::<PyList>() {
            l.iter().enumerate().map(|(k, v)| Ok((Expr::int(k as i64), to_expr(&v)?))).collect::<PyResult<_>>()?
        } else {
            // any other mapping, through `init.items()` as SymPy does
            let mut pairs = Vec::new();
            for pair in i.call_method0("items")?.try_iter()? {
                let (k, v): (Bound<'_, PyAny>, Bound<'_, PyAny>) = pair?.extract()?;
                pairs.push((to_expr(&k)?, to_expr(&v)?));
            }
            pairs
        }),
    };
    match recurr::rsolve(&f, &y, init.as_deref()) {
        Ok(Some(e)) => Ok(wrap(e).into_pyobject(py)?.into_any().unbind()),
        Ok(None) => Ok(py.None()),
        Err(RsolveError::Value(m)) => Err(PyValueError::new_err(m)),
        Err(RsolveError::Type(m)) => Err(PyTypeError::new_err(m)),
    }
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(rsolve, module)?)?;
    module.add_function(wrap_pyfunction!(rsolve_poly, module)?)?;
    module.add_function(wrap_pyfunction!(rsolve_ratio, module)?)?;
    module.add_function(wrap_pyfunction!(rsolve_hyper, module)?)?;
    module.add_function(wrap_pyfunction!(difference_delta, module)?)?;
    module.add_function(wrap_pyfunction!(hypersimp, module)?)?;
    module.add_function(wrap_pyfunction!(hypersimilar, module)?)?;
    Ok(())
}
