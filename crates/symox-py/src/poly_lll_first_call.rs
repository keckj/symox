//! Native evidence for exact first-call LLL admission; no classifier warmup.
#[path="lll_cold_getters.rs"]
mod getters;
#[path="poly_lll_proof_support.rs"]
mod support;
#[path="lll_cold_reference.rs"]
mod reference;
#[path="lll_cold_graph.rs"]
mod graph;
use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyDict, PyTuple};

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    support::register(module)?;
    reference::register(module)?;
    graph::register(module)?;
    let py = module.py();
    let api = getters::Getters::load();
    let helper = PyCFunction::new_closure(py, Some(c"_lll_safe_function_state"), None,
        move |args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>|
              -> PyResult<Option<(PyObject, PyObject, PyObject, PyObject, PyObject)>> {
            if kwargs.is_some() || args.len() != 1 {
                return Err(pyo3::exceptions::PyTypeError::new_err("expected one function"));
            }
            let Some(api) = api else { return Ok(None); };
            let object = args.get_item(0)?;
            if !api.is_function(&object) { return Ok(None); }
            Ok(Some((api.code(&object)?.unbind(), api.globals(&object)?.unbind(),
                api.defaults(&object)?.unbind(), api.kwdefaults(&object)?.unbind(),
                api.closure(&object)?.unbind())))
        })?;
    module.add("_lll_safe_function_state", helper)
}
