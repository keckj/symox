//! Python bindings for targeted simplification routines.

use pyo3::prelude::*;

use symox::simp;

use crate::{to_expr, wrap, PyExpr};

#[pyfunction]
fn trigsimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simp::trigsimp(&to_expr(e)?)))
}

#[pyfunction]
fn expand_trig(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simp::expand_trig(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(signature = (e, deep=true, combine="all", force=false))]
fn powsimp(e: &Bound<'_, PyAny>, deep: bool, combine: &str, force: bool) -> PyResult<PyExpr> {
    let _ = (deep, combine, force);
    Ok(wrap(simp::powsimp(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(signature = (e, deep=true, force=false))]
fn expand_log(e: &Bound<'_, PyAny>, deep: bool, force: bool) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(wrap(simp::expand_log(&to_expr(e)?, force)))
}

#[pyfunction]
#[pyo3(signature = (e, force=false))]
fn logcombine(e: &Bound<'_, PyAny>, force: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::logcombine(&to_expr(e)?, force)))
}

#[pyfunction]
#[pyo3(signature = (e, deep=true, force=false))]
fn expand_power_base(e: &Bound<'_, PyAny>, deep: bool, force: bool) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(wrap(simp::expand_power_base(&to_expr(e)?, force)))
}

#[pyfunction]
#[pyo3(signature = (e, deep=true))]
fn expand_power_exp(e: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(wrap(simp::expand_power_exp(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(signature = (e, deep=true))]
fn expand_mul(e: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(wrap(simp::expand_mul(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(signature = (e, constants=None, tolerance=None, full=false, rational=None))]
fn nsimplify(
    e: &Bound<'_, PyAny>,
    constants: Option<&Bound<'_, PyAny>>,
    tolerance: Option<f64>,
    full: bool,
    rational: Option<bool>,
) -> PyResult<PyExpr> {
    let _ = (tolerance, full, rational);
    let mut cs = Vec::new();
    if let Some(c) = constants {
        for x in c.try_iter()? {
            cs.push(to_expr(&x?)?);
        }
    }
    Ok(wrap(simp::nsimplify(&to_expr(e)?, &cs)))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(trigsimp, m)?)?;
    m.add_function(wrap_pyfunction!(expand_trig, m)?)?;
    m.add_function(wrap_pyfunction!(powsimp, m)?)?;
    m.add_function(wrap_pyfunction!(expand_log, m)?)?;
    m.add_function(wrap_pyfunction!(logcombine, m)?)?;
    m.add_function(wrap_pyfunction!(expand_power_base, m)?)?;
    m.add_function(wrap_pyfunction!(expand_power_exp, m)?)?;
    m.add_function(wrap_pyfunction!(expand_mul, m)?)?;
    m.add_function(wrap_pyfunction!(nsimplify, m)?)?;
    Ok(())
}
