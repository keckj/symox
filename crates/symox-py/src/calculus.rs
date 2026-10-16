//! Local complex-analysis bindings.
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use crate::{to_expr, wrap, PyExpr};

#[pyfunction]
fn residue(expr: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, x0: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let (expression, variable, point) = (to_expr(expr)?, to_expr(x)?, to_expr(x0)?);
    if !variable.is_symbol() || symox::traverse::has(&point, &variable) {
        return Err(PyValueError::new_err("residue requires a symbol and an independent point"));
    }
    symox::calculus::residue(&expression, &variable, &point).map(wrap).map_err(PyNotImplementedError::new_err)
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(residue, module)?)?;
    Ok(())
}
