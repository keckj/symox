//! Python bindings for `dsolve`.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::Expr;

use crate::{to_expr, wrap};

#[pyfunction]
#[pyo3(signature = (eq, func=None, ics=None, **_hints))]
fn dsolve(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    func: Option<&Bound<'_, PyAny>>,
    ics: Option<&Bound<'_, PyDict>>,
    _hints: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyObject> {
    if let Some(hints) = _hints {
        for (key, value) in hints.iter() {
            let key: String = key.extract()?;
            let supported = match key.as_str() {
                "hint" => value.extract::<String>()? == "default",
                "simplify" => value.is_truthy()?,
                _ => false,
            };
            if !supported {
                return Err(PyNotImplementedError::new_err(format!("dsolve option '{key}' with this value is not implemented")));
            }
        }
    }
    let e = to_expr(eq)?;
    let f = match func {
        Some(f) => to_expr(f)?,
        None => find_function(&e).ok_or_else(|| PyNotImplementedError::new_err("cannot determine the unknown function"))?,
    };
    let mut conds = Vec::new();
    if let Some(d) = ics {
        for (k, v) in d.iter() {
            conds.push((to_expr(&k)?, to_expr(&v)?));
        }
    }
    let sols = symox::ode::dsolve(&e, &f, &conds)
        .ok_or_else(|| PyNotImplementedError::new_err(format!("dsolve: cannot solve {e}")))?;
    if sols.len() == 1 {
        return Ok(wrap(sols[0].clone()).into_pyobject(py)?.into_any().unbind());
    }
    Ok(PyList::new(py, sols.into_iter().map(wrap))?.into_any().unbind())
}

fn find_function(e: &Expr) -> Option<Expr> {
    let mut pending=vec![e.clone()];
    let mut seen=rustc_hash::FxHashMap::default();
    while let Some(node)=pending.pop() {
        if seen.insert(node.dag_id(),node.clone()).is_some() {continue;}
        if seen.len()>100000 {return None;}
        if let symox::Kind::Func(symox::Func::Undefined(head),args)=node.kind() {
            if !head.internal && args.len()==1 && args[0].is_symbol() {
                return Some(node);
            }
        }
        pending.extend(symox::inspection::scalar_children(&node).into_iter().rev());
    }
    None
}

/// `checkodesol(ode, sol)` -> `(True, 0)` when the solution checks.
#[pyfunction]
#[pyo3(signature = (ode, sol, func=None))]
fn checkodesol<'py>(py: Python<'py>, ode: &Bound<'py, PyAny>, sol: &Bound<'py, PyAny>, func: Option<&Bound<'py, PyAny>>) -> PyResult<Bound<'py, PyTuple>> {
    let e = to_expr(ode)?;
    let s = to_expr(sol)?;
    let f = match func {
        Some(f) => to_expr(f)?,
        None => find_function(&e).ok_or_else(|| PyNotImplementedError::new_err("cannot determine the unknown function"))?,
    };
    let y = match s.kind() {
        symox::Kind::Rel(_, _, r) => r.clone(),
        _ => s.clone(),
    };
    let ok = symox::ode::checkodesol_checked(&e, &f, &y).ok_or_else(|| PyNotImplementedError::new_err(format!("checkodesol: ODE order above {} exceeds the native resource limit", symox::ode::CHECKODESOL_MAX_ORDER)))?;
    PyTuple::new(py, [ok.into_pyobject(py)?.to_owned().into_any(), wrap(if ok { Expr::zero() } else { y }).into_pyobject(py)?.into_any()])
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dsolve, m)?)?;
    m.add_function(wrap_pyfunction!(checkodesol, m)?)?;
    Ok(())
}
