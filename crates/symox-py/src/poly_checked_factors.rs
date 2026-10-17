//! Checked typed factor views; no expression-printer round trips.
use super::{division, err, w, PolyValue, PyPoly};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use symox::poly::{checked_factors as native, factor_budget::Limits};
pub(super) fn ground(p: &PolyValue, py: Python<'_>) -> PyResult<PyObject> {
    let result = PyDict::new(py);
    for (root, multiplicity) in native::ground_roots(p, Limits::default()).map_err(err)? {
        result.set_item(w(root), multiplicity)?;
    }
    Ok(result.into_any().unbind())
}
pub(super) fn include(p: &PolyValue, squarefree: bool, all: bool) -> PyResult<Vec<(PyPoly, u32)>> {
    Ok(native::include(p, squarefree, all, Limits::default())
        .map_err(err)?
        .into_iter()
        .map(|(p, k)| (PyPoly { p }, k))
        .collect())
}
#[pyfunction]
#[pyo3(signature=(f,*gens,**kwargs))]
fn ground_roots<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyObject> {
    if kwargs.is_some_and(|kw| kw.contains("polys").unwrap_or(false)) {
        return Err(division::FlagError::new_err(
            "'polys' flag is not allowed in ground_roots",
        ));
    }
    let (p, _) = division::prepare_unary(f, gens, kwargs)?;
    ground(&p, f.py())
}
pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ground_roots, m)?)?;
    Ok(())
}
