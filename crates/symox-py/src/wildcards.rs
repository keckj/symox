//! Owned Python predicate tuples with weak interning of live symbol identities.
use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use std::collections::HashMap;
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc, Mutex, OnceLock, Weak,
};
use symox::expr::Expr;
use symox::symbol::{WildcardError, WildcardProperties};

type Registry = HashMap<Option<isize>, Vec<Weak<Properties>>>;
fn registry() -> &'static Mutex<Registry> {
    static VALUE: OnceLock<Mutex<Registry>> = OnceLock::new();
    VALUE.get_or_init(|| Mutex::new(HashMap::new()))
}
struct Properties {
    id: u64,
    hash: Option<isize>,
    values: Py<PyTuple>,
}
impl std::fmt::Debug for Properties {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("PythonWildcardProperties")
            .field(&self.id)
            .finish()
    }
}
impl Drop for Properties {
    fn drop(&mut self) {
        let mut entries = registry().lock().unwrap_or_else(|p| p.into_inner());
        if let Some(bucket) = entries.get_mut(&self.hash) {
            bucket.retain(|entry| entry.strong_count() != 0);
            if bucket.is_empty() {
                entries.remove(&self.hash);
            }
        }
    }
}
impl WildcardProperties for Properties {
    fn identity(&self) -> u64 {
        self.id
    }
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
    fn evaluate(&self, candidate: &Expr) -> Result<bool, WildcardError> {
        Python::with_gil(|py| -> PyResult<bool> {
            let candidate = crate::expr_object(py, candidate.clone())?;
            for predicate in self.values.bind(py).iter() {
                if !predicate.call1((&candidate,))?.is_truthy()? {
                    return Ok(false);
                }
            }
            // A successful reference binding inserts Wild into a Python dict.
            // Preserve its hash error even when the native caller consumes
            // bindings directly (replace, find or collect).
            self.values.bind(py).hash()?;
            Ok(true)
        })
        .map_err(|e| Box::new(e) as WildcardError)
    }
}

pub(crate) fn intern(values: &Bound<'_, PyTuple>) -> PyResult<Option<Arc<dyn WildcardProperties>>> {
    if values.is_empty() {
        return Ok(None);
    }
    let py = values.py();
    let hash = match values.hash() {
        Ok(hash) => Some(hash),
        Err(error) if error.is_instance_of::<PyTypeError>(py) => None,
        Err(error) => return Err(error),
    };
    // Python equality may reenter the constructor, so it runs without the lock.
    let candidates: Vec<_> = {
        let entries = registry().lock().unwrap_or_else(|p| p.into_inner());
        entries
            .get(&hash)
            .into_iter()
            .flatten()
            .filter_map(Weak::upgrade)
            .collect()
    };
    for candidate in candidates {
        if values.eq(candidate.values.bind(py))? {
            return Ok(Some(candidate));
        }
    }
    static NEXT: AtomicU64 = AtomicU64::new(1);
    let value = Arc::new(Properties {
        id: NEXT.fetch_add(1, Ordering::Relaxed),
        hash,
        values: values.clone().unbind(),
    });
    registry()
        .lock()
        .unwrap_or_else(|p| p.into_inner())
        .entry(hash)
        .or_default()
        .push(Arc::downgrade(&value));
    Ok(Some(value))
}

pub(crate) fn properties<'py>(
    py: Python<'py>,
    symbol: &symox::symbol::Symbol,
) -> PyResult<Bound<'py, PyTuple>> {
    match symbol.wild_properties() {
        None => Ok(PyTuple::empty(py)),
        Some(value) => value
            .as_any()
            .downcast_ref::<Properties>()
            .map(|value| value.values.bind(py).clone())
            .ok_or_else(|| {
                PyNotImplementedError::new_err(
                    "native wildcard predicates have no Python representation",
                )
            }),
    }
}

pub(crate) fn to_python_error(error: WildcardError) -> PyErr {
    if error.is::<symox::matches::ProtocolTypeError>() {
        return PyTypeError::new_err(error.to_string());
    }
    match error.downcast::<PyErr>() {
        Ok(error) => *error,
        Err(error) => PyValueError::new_err(error.to_string()),
    }
}

/// Mirror Python hash errors for predicate tuples, including nested Wild atoms.
pub(crate) fn check_hashable(py: Python<'_>, expression: &Expr) -> PyResult<()> {
    let mut pending = vec![expression.clone()];
    let mut seen = rustc_hash::FxHashMap::default();
    while let Some(node) = pending.pop() {
        // Visited nodes stay alive (see `Expr::dag_id`).
        if !node.has_wild_properties() || seen.insert(node.dag_id(), node.clone()).is_some() {
            continue;
        }
        if let Some(symbol) = node.as_sym() {
            properties(py, symbol)?.hash()?;
            pending.extend(symbol.exclude().iter().cloned());
        } else {
            pending.extend(node.args());
        }
    }
    Ok(())
}
