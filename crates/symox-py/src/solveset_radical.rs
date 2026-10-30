//! Python binding of the native radical and negative-base steps of
//! `solveset` (`symox::sets::solveset_radical`). The recursive `solveset`,
//! `solve` (for `unrad`'s `_vsolve`), `imageset` and `intersect` calls go
//! back to Python callbacks that take and return core expressions;
//! everything else is native.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use symox::sets::{self as st, StepHost};
use symox::solve::unrad::{UFail, UnradHost};
use symox::Expr;

use crate::unrad::{host_err, to_py_err, F};
use crate::{to_expr, wrap};

struct Host<'py> {
    py: Python<'py>,
    solveset: Bound<'py, PyAny>,
    imageset: Bound<'py, PyAny>,
    intersect: Bound<'py, PyAny>,
    solve: Bound<'py, PyAny>,
    fallback: Option<Bound<'py, PyAny>>,
}

fn back(r: PyResult<Bound<'_, PyAny>>) -> Result<Expr, F> {
    to_expr(&r.map_err(UFail::Host)?).map_err(UFail::Host)
}

impl UnradHost for Host<'_> {
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
    fn take(&mut self, _e: &Expr) -> Option<Result<bool, F>> {
        None
    }
}

impl StepHost for Host<'_> {
    fn solveset(&mut self, f: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, F> {
        back(self.solveset.call1((wrap(f.clone()), wrap(x.clone()), wrap(domain.clone()))))
    }
    fn imageset(&mut self, lam: &Expr, base: &Expr) -> Result<Expr, F> {
        back(self.imageset.call1((wrap(lam.clone()), wrap(base.clone()))))
    }
    fn intersect(&mut self, a: &Expr, b: &Expr) -> Result<Expr, F> {
        back(self.intersect.call1((wrap(a.clone()), wrap(b.clone()))))
    }
    fn is_value_error(&self, e: &PyErr) -> bool {
        e.is_instance_of::<PyValueError>(self.py)
    }
    fn unrad_fallback(&mut self, f: &Expr, x: &Expr) -> Result<Option<(Expr, Vec<Expr>)>, F> {
        let py = self.py;
        let Some(cb) = self.fallback.as_ref() else { return Ok(None) };
        let r = cb.call1((wrap(f.clone()), wrap(x.clone()))).map_err(|e| host_err(py, e))?;
        if r.is_none() {
            return Ok(None);
        }
        let eq = to_expr(&r.get_item(0).map_err(UFail::Host)?).map_err(UFail::Host)?;
        let mut cov = Vec::new();
        for item in r.get_item(1).map_err(UFail::Host)?.try_iter().map_err(UFail::Host)? {
            cov.push(to_expr(&item.map_err(UFail::Host)?).map_err(UFail::Host)?);
        }
        Ok(Some((eq, cov)))
    }
}

type Step = fn(&mut Host<'_>, &Expr, &Expr, &Expr) -> Result<Option<Expr>, F>;

#[allow(clippy::too_many_arguments)]
fn run<'py>(
    py: Python<'py>,
    step: Step,
    f: &Bound<'py, PyAny>,
    symbol: &Bound<'py, PyAny>,
    domain: &Bound<'py, PyAny>,
    solveset: Bound<'py, PyAny>,
    imageset: Bound<'py, PyAny>,
    intersect: Bound<'py, PyAny>,
    solve: Bound<'py, PyAny>,
    fallback: Option<Bound<'py, PyAny>>,
) -> PyResult<Option<crate::PyExpr>> {
    let (f, x, d) = (to_expr(f)?, to_expr(symbol)?, to_expr(domain)?);
    let mut h = Host { py, solveset, imageset, intersect, solve, fallback };
    match step(&mut h, &f, &x, &d) {
        Ok(r) => Ok(r.map(wrap)),
        Err(e) => Err(to_py_err(py, e)),
    }
}

/// `radical_step(f, symbol, domain, ...)`: the reference `_solve_radical`
/// result, or None to keep the native one. The callbacks take and return
/// core expressions; `unrad_fallback(f, x)` is the elimination used where
/// the reference `unrad` raises NotImplementedError.
#[pyfunction]
#[pyo3(signature = (f, symbol, domain, solveset, imageset, intersect, solve, unrad_fallback=None))]
#[allow(clippy::too_many_arguments)]
fn _solveset_radical_step<'py>(
    py: Python<'py>,
    f: &Bound<'py, PyAny>,
    symbol: &Bound<'py, PyAny>,
    domain: &Bound<'py, PyAny>,
    solveset: Bound<'py, PyAny>,
    imageset: Bound<'py, PyAny>,
    intersect: Bound<'py, PyAny>,
    solve: Bound<'py, PyAny>,
    unrad_fallback: Option<Bound<'py, PyAny>>,
) -> PyResult<Option<crate::PyExpr>> {
    run(py, |h, f, x, d| st::radical_step(h, f, x, d), f, symbol, domain, solveset, imageset, intersect, solve, unrad_fallback)
}

/// `negative_base_step(f, symbol, domain, ...)`, or None.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn _solveset_negative_base_step<'py>(
    py: Python<'py>,
    f: &Bound<'py, PyAny>,
    symbol: &Bound<'py, PyAny>,
    domain: &Bound<'py, PyAny>,
    solveset: Bound<'py, PyAny>,
    imageset: Bound<'py, PyAny>,
    intersect: Bound<'py, PyAny>,
    solve: Bound<'py, PyAny>,
) -> PyResult<Option<crate::PyExpr>> {
    run(py, |h, f, x, d| st::negative_base_step(h, f, x, d), f, symbol, domain, solveset, imageset, intersect, solve, None)
}

/// `negative_base_guard(f, symbol, domain, result)`.
#[pyfunction]
fn _solveset_negative_base_guard(f: &Bound<'_, PyAny>, symbol: &Bound<'_, PyAny>, domain: &Bound<'_, PyAny>, result: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    Ok(wrap(st::negative_base_guard(&to_expr(f)?, &to_expr(symbol)?, &to_expr(domain)?, &to_expr(result)?)))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_solveset_radical_step, m)?)?;
    m.add_function(wrap_pyfunction!(_solveset_negative_base_step, m)?)?;
    m.add_function(wrap_pyfunction!(_solveset_negative_base_guard, m)?)
}
