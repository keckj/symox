//! Metadata-only building blocks. No classification and no verdict cache.
use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyDict, PyModule, PyString, PyTuple, PyType};

#[pyfunction]
pub(super) fn _lll_module_registry(py:Python<'_>)->PyResult<Option<PyObject>> {
    let modules=unsafe {Bound::<PyAny>::from_borrowed_ptr(py,pyo3::ffi::PyImport_GetModuleDict())};
    let Ok(dictionary)=modules.downcast_exact::<PyDict>() else {return Ok(None);};
    if dictionary.keys().iter().any(|key|key.downcast_exact::<PyString>().is_err()){return Ok(None);}
    Ok(Some(modules.unbind()))
}

#[pyfunction]
pub(super) fn _lll_loaded_modules(names:&Bound<'_,PyAny>, py:Python<'_>)->PyResult<Option<PyObject>> {
    let Ok(names)=names.downcast_exact::<PyTuple>() else {return Ok(None);};
    let modules=unsafe { Bound::<PyAny>::from_borrowed_ptr(py,pyo3::ffi::PyImport_GetModuleDict()) };
    let Ok(modules)=modules.downcast_exact::<PyDict>() else {return Ok(None);};
    // No dictionary lookup precedes validation of every possible colliding key.
    if modules.keys().iter().any(|key|key.downcast_exact::<PyString>().is_err()) {return Ok(None);}
    let mut selected=Vec::with_capacity(names.len());
    for name in names.iter() {
        if name.downcast_exact::<PyString>().is_err(){return Ok(None);}
        let Some(value)=modules.get_item(name)? else {return Ok(None);};
        let Ok(module)=value.downcast_exact::<PyModule>() else {return Ok(None);};
        let namespace=module.dict();
        if namespace.keys().iter().any(|key|key.downcast_exact::<PyString>().is_err()){return Ok(None);}
        selected.push((module.clone().unbind(),namespace.unbind()));
    }
    Ok(Some(PyTuple::new(py,selected)?.into_any().unbind()))
}

#[pyfunction]
pub(super) fn _lll_weakref_target(reference:&Bound<'_,PyAny>)->PyResult<PyObject> {
    #[allow(deprecated)]
    let target=unsafe {pyo3::ffi::PyWeakref_GetObject(reference.as_ptr())};
    if target.is_null(){return Err(PyErr::fetch(reference.py()));}
    Ok(unsafe {Bound::<PyAny>::from_borrowed_ptr(reference.py(),target)}.unbind())
}

#[pyfunction]
pub(super) fn _lll_native_owner(value:&Bound<'_,PyAny>,name:&Bound<'_,PyAny>,owner:&Bound<'_,PyAny>)->PyResult<bool> {
    if name.downcast_exact::<PyString>().is_err(){return Ok(false);}
    let Ok(owner)=owner.downcast_exact::<PyModule>() else {return Ok(false);};
    let namespace=owner.dict();
    if namespace.keys().iter().any(|key|key.downcast_exact::<PyString>().is_err()){return Ok(false);}
    let Some(module_name)=namespace.get_item("__name__")? else {return Ok(false);};
    if module_name.downcast_exact::<PyString>().is_err(){return Ok(false);}
    if let Ok(class)=value.downcast::<PyType>() {
        if !class.get_type().is(&value.py().get_type::<PyType>())
            || unsafe{pyo3::ffi::PyType_GetFlags(class.as_type_ptr())}&pyo3::ffi::Py_TPFLAGS_HEAPTYPE!=0 {
            return Ok(false);
        }
        return Ok(class.getattr("__name__")?.eq(name)? && class.getattr("__module__")?.eq(module_name)?);
    }
    if value.downcast_exact::<PyCFunction>().is_err()
        || unsafe{pyo3::ffi::PyCFunction_GetSelf(value.as_ptr())}!=owner.as_ptr(){return Ok(false);}
    value.getattr("__name__")?.eq(name)
}

fn authenticated_dump(dump: &Bound<'_,PyAny>) -> PyResult<bool> {
    if dump.downcast_exact::<PyCFunction>().is_err() { return Ok(false); }
    let owner = unsafe { pyo3::ffi::PyCFunction_GetSelf(dump.as_ptr()) };
    if owner.is_null() || unsafe { pyo3::ffi::PyModule_CheckExact(owner) } == 0 {
        return Ok(false);
    }
    let owner = unsafe { Bound::<PyAny>::from_borrowed_ptr(dump.py(), owner) };
    let namespace = owner.downcast_exact::<PyModule>()?.dict();
    for (key,_) in namespace.iter() {
        if key.downcast_exact::<PyString>().is_err() { return Ok(false); }
    }
    let Some(name) = namespace.get_item("__name__")? else { return Ok(false); };
    if name.downcast_exact::<PyString>().is_err() || !name.eq("_abc")? { return Ok(false); }
    let name = dump.getattr("__name__")?;
    Ok(name.downcast_exact::<PyString>().is_ok() && name.eq("_get_dump")?)
}

pub(super) fn class_state<'py>(value: &Bound<'py, PyAny>) -> PyResult<Option<(Bound<'py, PyAny>, Bound<'py, PyAny>, Bound<'py, PyAny>, Bound<'py, PyType>)>> {
    let py = value.py();
    if value.downcast::<PyType>().is_err() { return Ok(None); }
    let descriptors = py.get_type::<PyType>().getattr("__dict__")?;
    let namespace = descriptors.get_item("__dict__")?.call_method1("__get__", (value,))?;
    for key in namespace.try_iter()? {
        if key?.downcast_exact::<PyString>().is_err() { return Ok(None); }
    }
    let mro = descriptors.get_item("__mro__")?.call_method1("__get__", (value,))?;
    let bases = descriptors.get_item("__bases__")?.call_method1("__get__", (value,))?;
    Ok(Some((namespace, mro, bases, value.get_type())))
}

#[pyfunction]
fn _lll_class_state(value: &Bound<'_, PyAny>) -> PyResult<Option<(PyObject,PyObject,PyObject,PyObject)>> {
    Ok(class_state(value)?.map(|(a,b,c,d)|(a.unbind(),b.unbind(),c.unbind(),d.into_any().unbind())))
}

#[pyfunction]
pub(super) fn _lll_type_subclasses(value: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<Option<PyObject>> {
    if value.downcast::<PyType>().is_err() { return Ok(None); }
    let native = py.get_type::<PyType>().getattr("__dict__")?.get_item("__subclasses__")?;
    Ok(Some(native.call1((value,))?.unbind()))
}

#[pyfunction]
pub(super) fn _lll_instance_state(value: &Bound<'_, PyAny>) -> PyResult<Option<(PyObject,bool)>> {
    // Public stable getter reads the actual instance dictionary rather than a
    // Python __dict__ property or user __getattribute__ implementation.
    let Some(api)=super::getters::Getters::load() else { return Ok(None); };
    let Some(dictionary)=api.dictionary(value)? else { return Ok(None); };
    let class = value.get_type();
    let generic = unsafe { pyo3::ffi::PyType_GetSlot(class.as_type_ptr(), pyo3::ffi::Py_tp_getattro) }
        == pyo3::ffi::PyObject_GenericGetAttr as *const () as *mut std::ffi::c_void;
    Ok(Some((dictionary.unbind(),generic)))
}

#[pyclass(frozen)]
struct ABCView { #[pyo3(get)] _abc_impl: PyObject }

#[pyfunction]
pub(super) fn _lll_abc_state(value: &Bound<'_, PyAny>, dump: &Bound<'_,PyAny>, py: Python<'_>) -> PyResult<Option<PyObject>> {
    if !authenticated_dump(dump)? { return Ok(None); }
    let Some((namespace,_,_,_)) = class_state(value)? else { return Ok(None); };
    let data = match namespace.get_item("_abc_impl") {
        Ok(data) => data,
        Err(error) if error.is_instance_of::<pyo3::exceptions::PyKeyError>(py) => return Ok(None),
        Err(error) => return Err(error),
    };
    // Only the authenticated native getter sees this private holder.
    let view = Bound::new(py, ABCView { _abc_impl: data.unbind() })?;
    let state = match dump.call1((view,)) {
        Ok(state) => state,
        Err(error) if error.is_instance_of::<pyo3::exceptions::PyTypeError>(py) => return Ok(None),
        Err(error) => return Err(error),
    };
    if state.downcast_exact::<PyTuple>().is_err() { return Ok(None); }
    Ok(Some(state.unbind()))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(_lll_loaded_modules,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_module_registry,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_weakref_target,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_native_owner,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_class_state,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_type_subclasses,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_instance_state,module)?)?;
    module.add_function(wrap_pyfunction!(_lll_abc_state,module)?)
}
