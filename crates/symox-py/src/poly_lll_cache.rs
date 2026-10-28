//! Read-only rejection of a known-cold ABC cache, before expensive admission.
//! This never classifies a coefficient, mutates a cache, or admits native LLL.
use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyDict, PyModule, PySet, PyString, PyTuple, PyType};

#[pyclass(frozen)]
struct CacheView {
    #[pyo3(get)]
    _abc_impl: Py<PyAny>,
}

fn native_dump(function: &Bound<'_, PyAny>, owner: &Bound<'_, PyAny>) -> PyResult<bool> {
    let Ok(owner) = owner.downcast_exact::<PyModule>() else { return Ok(false); };
    if function.downcast_exact::<PyCFunction>().is_err()
        || unsafe { pyo3::ffi::PyCFunction_GetSelf(function.as_ptr()) } != owner.as_ptr()
    { return Ok(false); }
    // Validate all keys before the sole name lookup, without encoding unrelated
    // metadata or dispatching hashing/equality on a foreign key.
    let namespace = owner.dict();
    for (i, (key, _)) in namespace.iter().enumerate() {
        if i % 64 == 0 { function.py().check_signals()?; }
        if key.downcast_exact::<PyString>().is_err() { return Ok(false); }
    }
    let Some(name) = namespace.get_item("__name__")? else { return Ok(false); };
    let Ok(name) = name.downcast_exact::<PyString>() else { return Ok(false); };
    // Exact unicode equality works on codepoints, including unpaired surrogates.
    Ok(name.as_any().eq("_abc")? && function.getattr("__name__")?.extract::<String>()? == "_get_dump")
}

/// Bootstrap captures are evidence only if the captured objects really are
/// immutable builtin primitives. Their live identities are checked by Python.
#[pyfunction]
fn _lll_bootstrap_stock(lookups: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<bool> {
    let Ok(lookups) = lookups.downcast_exact::<PyTuple>() else { return Ok(false); };
    let namespace = unsafe { Bound::from_borrowed_ptr(py, pyo3::ffi::PyEval_GetBuiltins()) };
    let Ok(namespace) = namespace.downcast_exact::<PyDict>() else { return Ok(false); };
    let factories = ["type", "object", "dict", "set", "frozenset", "bool", "float", "complex",
        "str", "bytes", "bytearray", "list", "tuple", "int", "range", "enumerate", "super", "zip"];
    for (i, state) in lookups.iter().enumerate() {
        if i % 64 == 0 { py.check_signals()?; }
        let Ok(state) = state.downcast_exact::<PyTuple>() else { return Ok(false); };
        if state.len() != 3 { return Ok(false); }
        let name = state.get_item(1)?;
        let Ok(name) = name.downcast_exact::<PyString>() else { return Ok(false); };
        let name = name.extract::<String>()?;
        let value = state.get_item(2)?;
        if factories.contains(&name.as_str()) {
            let Ok(class) = value.downcast_exact::<PyType>() else { return Ok(false); };
            if unsafe { pyo3::ffi::PyType_GetFlags(class.as_type_ptr()) }
                & pyo3::ffi::Py_TPFLAGS_HEAPTYPE != 0
            { return Ok(false); }
            if class.getattr("__name__")?.extract::<String>()? != name
                || class.getattr("__module__")?.extract::<String>()? != "builtins"
            { return Ok(false); }
        } else {
            if value.downcast_exact::<PyCFunction>().is_err() { return Ok(false); }
            let owner = unsafe { pyo3::ffi::PyCFunction_GetSelf(value.as_ptr()) };
            if owner.is_null() || unsafe { pyo3::ffi::PyModule_CheckExact(owner) } == 0
                || unsafe { pyo3::ffi::PyModule_GetDict(owner) } != namespace.as_ptr()
                || value.getattr("__name__")?.extract::<String>()? != name
            { return Ok(false); }
        }
    }
    Ok(true)
}

#[pyfunction]
fn _lll_cache_is_cold(
    basic: &Bound<'_, PyAny>, element: &Bound<'_, PyAny>,
    dump: &Bound<'_, PyAny>, owner: &Bound<'_, PyAny>, py: Python<'_>,
) -> PyResult<Option<bool>> {
    py.check_signals()?;
    if basic.downcast::<PyType>().is_err() || !native_dump(dump, owner)? {
        return Ok(None);
    }
    // Invoke type's immutable native dictionary descriptor directly. Neither
    // Basic's metaclass lookup hooks nor an overridden __dict__ are consulted.
    let descriptor = py.get_type::<PyType>().getattr("__dict__")?.get_item("__dict__")?;
    let dictionary = descriptor.call_method1("__get__", (basic,))?;
    // Validate every key before a lookup: a str subclass could otherwise
    // supply hashing/equality hooks, even when _abc_impl was encountered first.
    // No key encoding is needed, so unrelated surrogate names remain harmless.
    for (i, key) in dictionary.try_iter()?.enumerate() {
        if i % 64 == 0 { py.check_signals()?; }
        if key?.downcast_exact::<PyString>().is_err() { return Ok(None); }
    }
    let data = match dictionary.get_item("_abc_impl") {
        Ok(data) => data.unbind(),
        Err(error) if error.is_instance_of::<pyo3::exceptions::PyKeyError>(py) => return Ok(None),
        Err(error) => return Err(error),
    };
    // _abc._get_dump reads only _abc_impl. A native holder exposes the actual
    // data without running user attribute hooks on Basic or its metaclass.
    // CPython versions rejecting a holder simply retain the original guard.
    let view = Bound::new(py, CacheView { _abc_impl: data })?;
    let state = match dump.call1((view,)) {
        Ok(state) => state,
        Err(error) if error.is_instance_of::<pyo3::exceptions::PyTypeError>(py) => return Ok(None),
        Err(error) => return Err(error),
    };
    let Ok(state) = state.downcast_exact::<PyTuple>() else { return Ok(None); };
    if state.len() != 4 { return Ok(None); }
    let negative = state.get_item(2)?;
    let Ok(negative) = negative.downcast_exact::<PySet>() else { return Ok(None); };
    for (i, reference) in negative.iter().enumerate() {
        if i % 64 == 0 { py.check_signals()?; }
        // The stable C getter cannot invoke a weakref subclass's __call__.
        // Its borrowed result remains valid while reference and the GIL live.
        #[allow(deprecated)] // abi3-py39 must also work before PyWeakref_GetRef.
        let referent = unsafe { pyo3::ffi::PyWeakref_GetObject(reference.as_ptr()) };
        if referent.is_null() { return Err(PyErr::fetch(py)); }
        if referent == element.as_ptr() { return Ok(Some(false)); }
    }
    py.check_signals()?;
    Ok(Some(true))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(_lll_bootstrap_stock, module)?)?;
    module.add_function(wrap_pyfunction!(_lll_cache_is_cold, module)?)
}
