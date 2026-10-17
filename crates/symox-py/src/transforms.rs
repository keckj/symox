//! Integral transform bindings (`symox::transforms`).
use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use symox::special::SpecialFn;
use crate::{to_expr, wrap, PyExpr};

fn head(name: &str) -> PyResult<SpecialFn> {
    SpecialFn::from_name(name)
        .filter(|s| s.is_transform())
        .ok_or_else(|| PyValueError::new_err(format!("unknown integral transform {name}")))
}

/// `doit` of the held transform `name(*args)`: the result followed by the
/// auxiliary items (plane/strip and condition). Raises
/// `NotImplementedError` (SymPy's `IntegralTransformError`) when the
/// transform cannot be computed and `needeval` is set.
#[pyfunction]
#[pyo3(signature = (name, args, simplify=None, needeval=false))]
fn _transform_doit(name: &str, args: Vec<Bound<'_, PyAny>>, simplify: Option<bool>, needeval: bool) -> PyResult<Vec<PyExpr>> {
    let s = head(name)?;
    let args = args.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    if Some(args.len()) != s.arity() {
        return Err(PyTypeError::new_err(format!(
            "{name} takes exactly {} arguments ({} given)",
            s.arity().unwrap_or(0),
            args.len()
        )));
    }
    symox::transforms::transform_doit(s, &args, simplify, needeval)
        .map(|v| v.into_iter().map(wrap).collect())
        .map_err(PyNotImplementedError::new_err)
}

/// Held transform node `name(*args)`.
#[pyfunction]
fn _transform_held(name: &str, args: Vec<Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let s = head(name)?;
    let args = args.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    if Some(args.len()) != s.arity() {
        return Err(PyTypeError::new_err(format!(
            "{name} takes exactly {} arguments ({} given)",
            s.arity().unwrap_or(0),
            args.len()
        )));
    }
    Ok(wrap(symox::transforms::held(s, args)))
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(_transform_doit, module)?)?;
    module.add_function(wrap_pyfunction!(_transform_held, module)?)?;
    Ok(())
}

/// Argument `i` of a held transform node (`function`, `function_variable`,
/// `transform_variable`).
pub(crate) fn argument(e: &symox::Expr, i: usize) -> Option<PyExpr> {
    symox::transforms::as_held(e).and_then(|(_, args)| args.get(i).cloned()).map(wrap)
}

pub(crate) fn not_a_transform(name: &str) -> PyErr {
    pyo3::exceptions::PyAttributeError::new_err(format!("only integral transforms have {name}"))
}

/// Attributes implemented by the Python glue (`symox._transforms`).
pub(crate) fn python_attribute(py: Python<'_>, e: &PyExpr, name: &str) -> PyResult<PyObject> {
    if symox::transforms::as_held(&e.e).is_none() {
        return Err(not_a_transform(name));
    }
    Ok(py.import("symox._transforms")?.getattr("_attribute")?.call1((e.clone(), name))?.unbind())
}
