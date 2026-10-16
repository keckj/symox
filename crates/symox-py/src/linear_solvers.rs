//! Binding helpers; public argument normalization lives in _linear.py.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;
use symox::linear::{self, LinearError};
use crate::{to_expr, wrap, PyExpr};

pyo3::create_exception!(_core, NonlinearError, PyValueError);

fn error(value: LinearError) -> PyErr {
    match value {
        LinearError::Nonlinear(e) => NonlinearError::new_err(format!("nonlinear term: {e}")),
        LinearError::Invalid(message) => PyValueError::new_err(message),
        LinearError::Matrix(error) => PyValueError::new_err(error.to_string()),
    }
}

fn expressions(value: &Bound<'_, PyAny>) -> PyResult<Vec<symox::Expr>> {
    value.try_iter()?.map(|item| to_expr(&item?)).collect()
}

#[pyfunction]
fn _linear_eq_coefficients<'py>(py: Python<'py>, equations: &Bound<'_, PyAny>, symbols: &Bound<'_, PyAny>) -> PyResult<(Bound<'py, PyList>, Bound<'py, PyList>)> {
    let (a, b) = linear::linear_eq_to_matrix(&expressions(equations)?, &expressions(symbols)?).map_err(error)?;
    Ok((PyList::new(py, a.data().iter().cloned().map(wrap))?, PyList::new(py, b.data().iter().cloned().map(wrap))?))
}

#[pyfunction]
fn _linsolve_matrix(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, symbols: &Bound<'_, PyAny>) -> PyResult<Option<Vec<PyExpr>>> {
    let a = a.downcast::<crate::matrices::PyMatrix>().map_err(|_| PyTypeError::new_err("Expected a Matrix"))?;
    let b = b.downcast::<crate::matrices::PyMatrix>().map_err(|_| PyTypeError::new_err("Expected a Matrix"))?;
    Ok(linear::linsolve(&a.borrow().m, &b.borrow().m, &expressions(symbols)?).map_err(error)?.map(|values| values.into_iter().map(wrap).collect()))
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("NonlinearError", module.py().get_type::<NonlinearError>())?;
    module.add_function(wrap_pyfunction!(_linear_eq_coefficients, module)?)?;
    module.add_function(wrap_pyfunction!(_linsolve_matrix, module)?)?;
    Ok(())
}
