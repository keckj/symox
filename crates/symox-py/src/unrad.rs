//! Python binding of the native `unrad` (`symox::solve::unrad`). The
//! `_vsolve` calls and a caller's `_take` call back into Python; everything
//! else is native.

use pyo3::exceptions::{PyAttributeError, PyIndexError, PyNotImplementedError};
use pyo3::prelude::*;

use symox::solve::unrad::{self as ur, UFail, UnradHost};
use symox::Expr;

use crate::{to_expr, wrap};

struct PyHost<'py> {
    py: Python<'py>,
    solve: Bound<'py, PyAny>,
    take: Option<Bound<'py, PyAny>>,
}

pub(crate) type F = UFail<PyErr>;

pub(crate) fn host_err(py: Python<'_>, e: PyErr) -> F {
    if e.is_instance_of::<PyNotImplementedError>(py) {
        UFail::NotImplemented(String::new())
    } else {
        UFail::Host(e)
    }
}

pub(crate) fn to_py_err(py: Python<'_>, e: F) -> PyErr {
    match e {
        UFail::NotImplemented(m) => PyNotImplementedError::new_err(m),
        UFail::Poly(m) => crate::poly::reference_error(py, "PolynomialError", &m),
        UFail::Index(m) => PyIndexError::new_err(m),
        UFail::Attribute(m) => PyAttributeError::new_err(m),
        UFail::Host(e) => e,
    }
}

impl UnradHost for PyHost<'_> {
    type Error = PyErr;
    fn vsolve(&mut self, eq: &Expr, x: &Expr) -> Result<Vec<Expr>, F> {
        let py = self.py;
        let kw = pyo3::types::PyDict::new(py);
        kw.set_item("check", false).map_err(UFail::Host)?;
        kw.set_item("simplify", false).map_err(UFail::Host)?;
        let r = self.solve.call((wrap(eq.clone()), wrap(x.clone())), Some(&kw)).map_err(|e| host_err(py, e))?;
        let mut out = Vec::new();
        for item in r.try_iter().map_err(UFail::Host)? {
            out.push(to_expr(&item.map_err(UFail::Host)?).map_err(UFail::Host)?);
        }
        Ok(out)
    }
    fn take(&mut self, e: &Expr) -> Option<Result<bool, F>> {
        let py = self.py;
        let t = self.take.as_ref()?;
        Some(t.call1((wrap(e.clone()),)).and_then(|r| r.is_truthy()).map_err(|e| host_err(py, e)))
    }
}

/// `unrad(eq, *syms, _reverse=..., _take=...)` for an Expr `eq`; returns
/// None or `(eq, cov)`.
#[pyfunction]
#[pyo3(signature = (eq, syms, reverse=false, take=None, solve=None))]
fn _unrad(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    syms: Vec<Bound<'_, PyAny>>,
    reverse: bool,
    take: Option<Bound<'_, PyAny>>,
    solve: Option<Bound<'_, PyAny>>,
) -> PyResult<PyObject> {
    let e = to_expr(eq)?;
    let syms: Vec<Expr> = syms.iter().map(to_expr).collect::<PyResult<_>>()?;
    let solve = match solve {
        Some(s) => s,
        None => py.import("symox")?.getattr("solve")?,
    };
    let custom = take.is_some();
    let mut h = PyHost { py, solve, take };
    match ur::unrad(&mut h, &e, &syms, reverse, custom) {
        Ok(None) => Ok(py.None()),
        Ok(Some((eq, cov))) => {
            let cov: Vec<_> = cov.into_iter().map(wrap).collect();
            Ok((wrap(eq), cov).into_pyobject(py)?.into_any().unbind())
        }
        Err(e) => Err(to_py_err(py, e)),
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_unrad, m)?)?;
    Ok(())
}
