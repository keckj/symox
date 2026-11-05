//! Native loops of the closed-QQ splitfactor admission
//! (`python/symox/integrals/_risch_splitfactor_native.reduce`).
//!
//! The Python admission authenticates, on every call, the namespaces, the
//! function states and the global lookups recorded by `register_root` (about
//! 50 namespaces, 400 function states and 3 000 lookups). These two functions
//! run the same checks, in the same order, with the same identity
//! comparisons and the same attribute reads, without a Python-level loop:
//!
//! - `_risch_splitfactor_scopes(modules, scopes, builtins, builtin_aliases)`:
//!   every module is an exact module whose `__dict__` is the recorded one,
//!   every scope is an exact dict with exact str keys, every recorded
//!   builtin is unchanged. Returns a bool.
//! - `_risch_splitfactor_states(states, lookups, builtins, rational)`: the
//!   function-state loop and `_states_match` of `reduce`, then the lookup
//!   loop. Returns None when the admission must decline, otherwise the tuple
//!   of lookup indices whose current value is not the recorded object but
//!   lies in the rationaltools namespace (the caller decides those with the
//!   Python `_equivalent`, exactly as before). It raises the exceptions the
//!   Python loops raise (other than the `ValueError` of an empty cell, which
//!   declines). Without the public CPython function getters it returns None
//!   (the Python reader `_ddm_function_metadata` then answers None too).
//!
//! Nothing is cached: every field is read anew on each call.
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyModule, PyTuple};

#[path = "poly_function_metadata.rs"]
#[allow(dead_code)]
mod function_metadata;

fn plain(value: &Bound<'_, PyAny>) -> bool {
    let Ok(d) = value.downcast_exact::<PyDict>() else { return false };
    let mut pos: pyo3::ffi::Py_ssize_t = 0;
    let mut key: *mut pyo3::ffi::PyObject = std::ptr::null_mut();
    let mut val: *mut pyo3::ffi::PyObject = std::ptr::null_mut();
    // Borrowed references under the GIL; no Python code runs in the scan.
    unsafe {
        while pyo3::ffi::PyDict_Next(d.as_ptr(), &mut pos, &mut key, &mut val) != 0 {
            if pyo3::ffi::PyUnicode_CheckExact(key) == 0 {
                return false;
            }
        }
    }
    true
}

/// `d.get(key)` of an exact dict (None when absent), as the Python method.
fn get<'py>(d: &Bound<'py, PyDict>, key: &Bound<'py, PyAny>) -> PyResult<Option<Bound<'py, PyAny>>> {
    d.get_item(key)
}

#[pyfunction]
fn _risch_splitfactor_scopes(modules: &Bound<'_, PyAny>, scopes: &Bound<'_, PyAny>, builtins: &Bound<'_, PyAny>,
    builtin_aliases: &Bound<'_, PyAny>) -> PyResult<bool> {
    let (Ok(modules), Ok(scopes), Ok(builtins), Ok(aliases)) = (modules.downcast_exact::<PyTuple>(),
        scopes.downcast_exact::<PyTuple>(), builtins.downcast_exact::<PyDict>(), builtin_aliases.downcast_exact::<PyTuple>())
    else { return Ok(false) };
    for item in modules.iter() {
        let Ok(item) = item.downcast_exact::<PyTuple>() else { return Ok(false) };
        if item.len() != 2 {
            return Ok(false);
        }
        let module = item.get_item(0)?;
        let Ok(module) = module.downcast_exact::<PyModule>() else { return Ok(false) };
        if !module.dict().is(&item.get_item(1)?) {
            return Ok(false);
        }
    }
    for scope in scopes.iter() {
        if !plain(&scope) {
            return Ok(false);
        }
    }
    for item in aliases.iter() {
        let Ok(item) = item.downcast_exact::<PyTuple>() else { return Ok(false) };
        if item.len() != 2 {
            return Ok(false);
        }
        match get(builtins, &item.get_item(0)?)? {
            Some(v) if v.is(&item.get_item(1)?) => {}
            _ => return Ok(false),
        }
    }
    Ok(true)
}

/// `plain(d)`, remembering the dictionaries already scanned in this call
/// (no Python code runs during the call, so they cannot change meanwhile).
fn plain_memo(value: &Bound<'_, PyAny>, seen: &mut Vec<usize>) -> bool {
    let key = value.as_ptr() as usize;
    if seen.contains(&key) {
        return true;
    }
    let ok = plain(value);
    if ok {
        seen.push(key);
    }
    ok
}

/// One state of `reduce`'s first loop and of `_states_match`.
fn state_matches(api: &function_metadata::Getters, state: &Bound<'_, PyAny>, builtins: &Bound<'_, PyDict>,
    first: bool, seen: &mut Vec<usize>) -> PyResult<bool> {
    let Ok(s) = state.downcast_exact::<PyTuple>() else { return Ok(false) };
    if s.len() != 7 {
        return Ok(false);
    }
    let py = state.py();
    let function = s.get_item(0)?;
    if !api.is_function(&function) {
        return Ok(false);
    }
    let fields = (0..5).map(|i| api.field(&function, i)).collect::<PyResult<Vec<_>>>()?;
    for i in 0..4 {
        if !fields[i].is(&s.get_item(i + 1)?) {
            return Ok(false);
        }
    }
    if !fields[2].is_none() && !fields[2].is_exact_instance_of::<PyTuple>()
        || !fields[3].is_none() && !fields[3].is_exact_instance_of::<PyDict>()
    {
        return Ok(false);
    }
    if first {
        return Ok(fields[4].is_none() || fields[4].is_exact_instance_of::<PyTuple>());
    }
    if !function.getattr(pyo3::intern!(py, "__builtins__"))?.is(builtins) || !plain_memo(&fields[1], seen) {
        return Ok(false);
    }
    let current = function.getattr(pyo3::intern!(py, "__dict__"))?;
    let Ok(current) = current.downcast_exact::<PyDict>() else { return Ok(false) };
    let Ok(attributes) = s.get_item(6)?.downcast_into_exact::<PyTuple>() else { return Ok(false) };
    if current.len() != attributes.len() {
        return Ok(false);
    }
    for ((key, value), pair) in current.iter().zip(attributes.iter()) {
        let Ok(pair) = pair.downcast_exact::<PyTuple>() else { return Ok(false) };
        if pair.len() != 2 || !key.is(&pair.get_item(0)?) || !value.is(&pair.get_item(1)?) {
            return Ok(false);
        }
    }
    let Ok(contents) = s.get_item(5)?.downcast_into_exact::<PyTuple>() else { return Ok(false) };
    let closure = if fields[4].is_none() { PyTuple::empty(py) } else {
        match fields[4].downcast_exact::<PyTuple>() { Ok(t) => t.clone(), Err(_) => return Ok(false) }
    };
    if closure.len() != contents.len() {
        return Ok(false);
    }
    let name = pyo3::intern!(py, "cell_contents");
    for (cell, value) in closure.iter().zip(contents.iter()) {
        match cell.getattr(name) {
            Ok(v) if v.is(&value) => {}
            Ok(_) => return Ok(false),
            Err(e) if e.is_instance_of::<pyo3::exceptions::PyValueError>(py) => return Ok(false),
            Err(e) => return Err(e),
        }
    }
    Ok(true)
}

#[pyfunction]
fn _risch_splitfactor_states(py: Python<'_>, states: &Bound<'_, PyAny>, lookups: &Bound<'_, PyAny>, builtins: &Bound<'_, PyAny>,
    rational: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    let Some(api) = function_metadata::api() else { return Ok(None) };
    let (Ok(states), Ok(lookups), Ok(builtins)) = (states.downcast_exact::<PyTuple>(),
        lookups.downcast_exact::<PyTuple>(), builtins.downcast_exact::<PyDict>()) else { return Ok(None) };
    let mut seen = Vec::new();
    for first in [true, false] {
        for state in states.iter() {
            if !state_matches(api, &state, builtins, first, &mut seen)? {
                return Ok(None);
            }
        }
    }
    let mut pending = Vec::new();
    for (i, item) in lookups.iter().enumerate() {
        let Ok(item) = item.downcast_exact::<PyTuple>() else { return Ok(None) };
        if item.len() != 4 {
            return Ok(None);
        }
        let scope = item.get_item(0)?;
        let Ok(dict) = scope.downcast_exact::<PyDict>() else { return Ok(None) };
        let name = item.get_item(1)?;
        let current = get(dict, &name)?;
        let present = item.get_item(2)?;
        // `(name in scope) is not present`: present is the bool recorded by
        // register_root.
        let is_in = current.is_some();
        if !present.is(pyo3::types::PyBool::new(py, is_in).as_any()) {
            return Ok(None);
        }
        let fallback = get(builtins, &name)?;
        let current = current.or(fallback).unwrap_or_else(|| py.None().into_bound(py));
        if !current.is(&item.get_item(3)?) {
            if !scope.is(rational) {
                return Ok(None);
            }
            pending.push(i);
        }
    }
    Ok(Some(PyTuple::new(py, pending)?.into_any().unbind()))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let _ = function_metadata::api();
    m.add_function(wrap_pyfunction!(_risch_splitfactor_scopes, m)?)?;
    m.add_function(wrap_pyfunction!(_risch_splitfactor_states, m)?)
}
