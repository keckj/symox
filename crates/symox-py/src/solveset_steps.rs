//! `_solveset_modular`: the native `_solve_modular` step of `solveset`.
use pyo3::exceptions::{PyAttributeError, PyTypeError};
use pyo3::prelude::*;
use symox::sets::{self as st, ModularError};

/// The reference result of a modular equation, or None when `f` is not
/// modular in `symbol`. Raises the reference TypeError / NotInvertible.
#[pyfunction]
fn _solveset_modular(f: &Bound<'_, PyAny>, symbol: &Bound<'_, PyAny>, domain: &Bound<'_, PyAny>) -> PyResult<Option<crate::PyExpr>> {
    let (f, x, domain) = (crate::to_expr(f)?, crate::to_expr(symbol)?, crate::to_expr(domain)?);
    match st::solveset_modular(&f, &x, &domain) {
        Ok(r) => Ok(r.map(crate::wrap)),
        Err(ModularError::Type(message)) => Err(PyTypeError::new_err(message)),
        Err(ModularError::NotInvertible(message)) => Err(crate::poly::an_rational::NotInvertible::new_err(message)),
        Err(ModularError::Attribute(message)) => Err(PyAttributeError::new_err(message)),
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_solveset_modular, m)?)
}
