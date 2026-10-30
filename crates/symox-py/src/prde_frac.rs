//! Wave 30 (prdefrac): helpers for the parametric Risch DE admission.
//!
//! `_str_keyed(d)` is the namespace check of the splitfactor/PolyMatrix
//! admissions ("an exact dict whose keys are all exact str"), done without
//! a Python-level loop. It reads keys only: no hashing, no equality, no
//! Python callbacks.
//!
//! `_function_states_unchanged(records)` is `_polymatrix_native._unchanged`
//! over a tuple of `(function, code, defaults, kwdefaults, cell_contents)`
//! records, reading the live fields through the public CPython getters of
//! `poly_function_metadata` (no `object.__getattr__` audit event). It returns
//! None when those getters are unavailable (the caller then runs the Python
//! check), and False for any malformed record or empty cell.
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

// The same public-getter module as the DDM admissions (resolved once, no
// cached function data); included here so the shared file stays untouched.
#[path = "poly_function_metadata.rs"]
#[allow(dead_code)]
mod function_metadata;

#[pyfunction]
fn _str_keyed(value: &Bound<'_, PyAny>) -> bool {
    let Ok(d) = value.downcast_exact::<PyDict>() else { return false };
    // Borrowed references only, under the GIL; nothing here runs Python code,
    // so the dictionary cannot change during the scan.
    let mut pos: pyo3::ffi::Py_ssize_t = 0;
    let mut key: *mut pyo3::ffi::PyObject = std::ptr::null_mut();
    let mut val: *mut pyo3::ffi::PyObject = std::ptr::null_mut();
    unsafe {
        while pyo3::ffi::PyDict_Next(d.as_ptr(), &mut pos, &mut key, &mut val) != 0 {
            if pyo3::ffi::PyUnicode_CheckExact(key) == 0 {
                return false;
            }
        }
    }
    true
}

fn record_unchanged(api: &function_metadata::Getters, record: &Bound<'_, PyAny>) -> PyResult<bool> {
    let Ok(r) = record.downcast_exact::<PyTuple>() else { return Ok(false) };
    if r.len() != 5 {
        return Ok(false);
    }
    let function = r.get_item(0)?;
    if !api.is_function(&function) {
        return Ok(false);
    }
    let py = record.py();
    let (code, defaults, keywords, contents) = (r.get_item(1)?, r.get_item(2)?, r.get_item(3)?, r.get_item(4)?);
    if !api.field(&function, 0)?.is(&code)
        || !api.field(&function, 2)?.is(&defaults)
        || !api.field(&function, 3)?.is(&keywords)
    {
        return Ok(false);
    }
    if !defaults.is_none() && !defaults.is_exact_instance_of::<PyTuple>()
        || !keywords.is_none() && !keywords.is_exact_instance_of::<PyDict>()
    {
        return Ok(false);
    }
    let Ok(contents) = contents.downcast_exact::<PyTuple>() else { return Ok(false) };
    let cells = api.field(&function, 4)?;
    if cells.is_none() {
        return Ok(contents.is_empty());
    }
    let Ok(cells) = cells.downcast_exact::<PyTuple>() else { return Ok(false) };
    if cells.len() != contents.len() {
        return Ok(false);
    }
    let name = pyo3::intern!(py, "cell_contents");
    for (cell, value) in cells.iter().zip(contents.iter()) {
        match cell.getattr(name) {
            Ok(current) if current.is(&value) => {}
            _ => return Ok(false),
        }
    }
    Ok(true)
}

#[pyfunction]
fn _function_states_unchanged(records: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    let Some(api) = function_metadata::api() else { return Ok(None) };
    let Ok(records) = records.downcast_exact::<PyTuple>() else { return Ok(Some(false)) };
    for record in records.iter() {
        if !record_unchanged(api, &record)? {
            return Ok(Some(false));
        }
    }
    Ok(Some(true))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_str_keyed, m)?)?;
    m.add_function(wrap_pyfunction!(_function_states_unchanged, m)?)
}
