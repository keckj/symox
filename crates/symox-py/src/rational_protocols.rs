//! Read-only authentication of Python function snapshots used by rational tools.
//! Only exact tuples and canonical Python functions reach the identity loop.
//! The caller retains its original protocol loop when dependencies are replaced.
use pyo3::exceptions::{PyAttributeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyDict, PyTuple, PyType};

fn matches(
    py: Python<'_>, states: &Bound<'_, PyAny>, function_type: &Bound<'_, PyAny>,
    dependencies: &Bound<'_, PyAny>, expected: &[Py<PyAny>], canonical_function: &Py<PyType>,
) -> PyResult<Option<bool>> {
    py.check_signals()?;
    if !function_type.is(canonical_function.bind(py)) { return Ok(None); }
    let Ok(dependencies) = dependencies.downcast_exact::<PyTuple>() else { return Ok(None); };
    if dependencies.len() != expected.len() { return Ok(None); }
    for (dependency, original) in dependencies.iter().zip(expected) {
        if !dependency.is(original.bind(py)) { return Ok(None); }
    }
    let Ok(states) = states.downcast_exact::<PyTuple>() else { return Ok(Some(false)); };
    for state in states.iter() {
        py.check_signals()?;
        let Ok(state) = state.downcast_exact::<PyTuple>() else { return Ok(Some(false)); };
        if state.len() != 6 { return Ok(Some(false)); }
        let function = state.get_item(0)?;
        if !function.get_type().is(canonical_function.bind(py)) { return Ok(Some(false)); }
        let Ok(contents) = state.get_item(3)?.downcast_into_exact::<PyTuple>() else { return Ok(Some(false)); };
        for (attribute, index) in [(pyo3::intern!(py,"__code__"),1),(pyo3::intern!(py,"__closure__"),2),(pyo3::intern!(py,"__defaults__"),4),(pyo3::intern!(py,"__kwdefaults__"),5)] {
            if !function.getattr(attribute)?.is(&state.get_item(index)?) { return Ok(Some(false)); }
        }
        let closure = state.get_item(2)?;
        if closure.is_none() {
            if !contents.is_empty() { return Ok(Some(false)); }
            continue;
        }
        // Identity with function.__closure__ above proves these are actual cells,
        // so neither arbitrary iterators nor user attribute hooks run here.
        let closure = closure.downcast_exact::<PyTuple>()?;
        if closure.len() != contents.len() { return Ok(Some(false)); }
        for (cell, value) in closure.iter().zip(contents.iter()) {
            py.check_signals()?;
            if !cell.getattr(pyo3::intern!(py,"cell_contents"))?.is(&value) { return Ok(Some(false)); }
        }
    }
    py.check_signals()?;
    Ok(Some(true))
}

/// Native builtin ownership is read with C getters, so a foreign function
/// cannot authenticate itself through its name/module attributes or hooks.
fn builtin_in_namespace(function: &Bound<'_, PyAny>, namespace: &Bound<'_, PyDict>) -> bool {
    if function.downcast_exact::<PyCFunction>().is_err() { return false; }
    unsafe {
        let owner = pyo3::ffi::PyCFunction_GetSelf(function.as_ptr());
        !owner.is_null() && pyo3::ffi::PyModule_CheckExact(owner) != 0
            && pyo3::ffi::PyModule_GetDict(owner) == namespace.as_ptr()
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let builtins = unsafe { Bound::from_borrowed_ptr(py, pyo3::ffi::PyEval_GetBuiltins()) };
    let builtins = builtins.downcast::<PyDict>()?;
    let canonical_function = py.eval(c"lambda: None", None, None)?.get_type().unbind();
    let mut expected = Vec::new();
    for name in ["type", "all", "len", "zip", "tuple"] {
        expected.push(builtins.get_item(name)?.map_or_else(|| py.None(), Bound::unbind));
    }
    // A user function installed before importing the extension is not a stock
    // builtin. Refuse the optimization instead of capturing it as trusted.
    let enabled = expected[0].bind(py).is(py.get_type::<PyType>())
        && expected[4].bind(py).is(py.get_type::<PyTuple>())
        && expected[3].as_ptr() == std::ptr::addr_of_mut!(pyo3::ffi::PyZip_Type).cast()
        && [1usize,2].into_iter().all(|i| builtin_in_namespace(expected[i].bind(py), builtins))
        && expected[1].bind(py).getattr("__name__")?.extract::<String>()? == "all"
        && expected[2].bind(py).getattr("__name__")?.extract::<String>()? == "len";
    let helper = PyCFunction::new_closure(py, Some(c"_atan_states_match"), None,
        move |args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>| -> PyResult<Option<bool>> {
            if kwargs.is_some() || args.len() != 3 { return Err(PyTypeError::new_err("expected states, function type and dependencies")); }
            if !enabled { return Ok(None); }
            let py = args.py();
            match matches(py, &args.get_item(0)?, &args.get_item(1)?, &args.get_item(2)?, &expected, &canonical_function) {
                Err(error) if error.is_instance_of::<PyAttributeError>(py)
                    || error.is_instance_of::<PyValueError>(py) || error.is_instance_of::<PyTypeError>(py) => Ok(Some(false)),
                result => result,
            }
        })?;
    m.add("_atan_states_match", helper)
}
