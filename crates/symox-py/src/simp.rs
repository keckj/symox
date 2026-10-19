//! Python bindings for targeted simplification routines.

use pyo3::prelude::*;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};

use symox::simp;

use crate::{to_expr, wrap, PyExpr};

#[pyfunction]
fn trigsimp(expr: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simp::trigsimp(&to_expr(expr)?)))
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true))]
fn expand_trig(expr: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::expand_trig_deep(&to_expr(expr)?, deep)))
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true, combine="all", force=false))]
fn powsimp(expr: &Bound<'_, PyAny>, deep: bool, combine: &str, force: bool) -> PyResult<PyExpr> {
    match combine {
        "all" | "exp" => {},
        "base" => return Err(PyNotImplementedError::new_err("powsimp combine='base' is not implemented")),
        _ => return Err(PyValueError::new_err("combine must be 'all', 'exp' or 'base'")),
    }
    if force { return Err(PyNotImplementedError::new_err("powsimp force=True is not implemented")); }
    Ok(wrap(simp::powsimp_deep(&to_expr(expr)?, deep)))
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true, force=false))]
fn expand_log(expr: &Bound<'_, PyAny>, deep: bool, force: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::expand_log_deep(&to_expr(expr)?, force, deep)))
}

#[pyfunction]
#[pyo3(signature = (expr, force=false))]
fn logcombine(expr: &Bound<'_, PyAny>, force: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::logcombine(&to_expr(expr)?, force)))
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true, force=false))]
fn expand_power_base(expr: &Bound<'_, PyAny>, deep: bool, force: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::expand_power_base_deep(&to_expr(expr)?, force, deep)))
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true))]
fn expand_power_exp(expr: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::expand_power_exp_deep(&to_expr(expr)?, deep)))
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true))]
fn expand_mul(expr: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    Ok(wrap(simp::expand_mul_deep(&to_expr(expr)?, deep)))
}

#[pyfunction]
#[pyo3(signature = (expr, constants=None, tolerance=None, full=false, rational=None))]
fn nsimplify(
    expr: &Bound<'_, PyAny>,
    constants: Option<&Bound<'_, PyAny>>,
    tolerance: Option<f64>,
    full: bool,
    rational: Option<bool>,
) -> PyResult<PyExpr> {
    if let Some(tolerance) = tolerance {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(PyValueError::new_err("tolerance must be finite and positive"));
        }
        return Err(PyNotImplementedError::new_err("nsimplify tolerance is not implemented"));
    }
    if full { return Err(PyNotImplementedError::new_err("nsimplify full=True is not implemented")); }
    if rational.is_some() { return Err(PyNotImplementedError::new_err("nsimplify rational option is not implemented")); }
    let mut cs = Vec::new();
    if let Some(c) = constants {
        for x in c.try_iter()? {
            cs.push(to_expr(&x?)?);
        }
    }
    Ok(wrap(simp::nsimplify(&to_expr(expr)?, &cs)))
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
