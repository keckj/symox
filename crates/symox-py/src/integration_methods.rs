//! The Risch, manual and heurisch methods of `integrate`, provided to the
//! native integrator by the Python ports in `symox.integrals`.

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use symox::expr::{Expr, Func, Kind};
use symox::integrate::host_methods::{self, HostMethod, HostOutcome};
use symox::integrate::Conds;

use crate::{to_expr, wrap, PyExprNode};

fn conds_name(conds: Conds) -> &'static str {
    match conds {
        Conds::Piecewise => "piecewise",
        Conds::Separate => "separate",
        Conds::None => "none",
    }
}

/// The host callback of the native integrator. The ports catch what the
/// reference catches (NotImplementedError from `risch_integrate`,
/// PolynomialError from `heurisch`, ValueError and PolynomialError from
/// `manualintegrate`). When a `True` option selects the method, the port is
/// the whole reference computation and any other exception is raised when
/// the computation returns to Python, as the reference raises it out of
/// `integrate` (`integrate(a**(t - s)/(a**t + 1), t, heurisch=True)` raises
/// the TypeError of `factor_list`). In the default chain, where the native
/// rules precede the ported stages, an ordinary exception only leaves the
/// stage to the next method; an interrupt (or another BaseException) is
/// always raised.
fn integration_host(method: HostMethod, f: &Expr, x: &Expr, conds: Conds) -> HostOutcome {
    if symox::host::failure_pending() {
        return HostOutcome::Unavailable;
    }
    Python::with_gil(|py| {
        let run = || -> PyResult<HostOutcome> {
            let name = match method {
                HostMethod::Risch => "risch",
                HostMethod::Manual => "manual",
                HostMethod::Term => "term",
                HostMethod::ManualOnly => "manual_only",
                HostMethod::HeurischOnly => "heurisch_only",
            };
            let module = py.import("symox.integrals._methods")?;
            let result = module.getattr("_native_host")?
                .call1((name, wrap(f.clone()), wrap(x.clone()), conds_name(conds)))?;
            if result.is_none() {
                return Ok(HostOutcome::Unavailable);
            }
            if let Ok(pair) = result.downcast::<PyTuple>() {
                return Ok(HostOutcome::NonElementary(to_expr(&pair.get_item(0)?)?, to_expr(&pair.get_item(1)?)?));
            }
            Ok(HostOutcome::Value(to_expr(&result)?))
        };
        match run() {
            Ok(outcome) => outcome,
            Err(error) => {
                let selected = match method {
                    HostMethod::ManualOnly | HostMethod::HeurischOnly => true,
                    HostMethod::Risch => host_methods::methods().risch == Some(true),
                    HostMethod::Manual | HostMethod::Term => false,
                };
                if selected || !error.is_instance_of::<pyo3::exceptions::PyException>(py) {
                    crate::host::report(py, error, None);
                }
                HostOutcome::Unavailable
            }
        }
    })
}

pub(crate) fn register() {
    host_methods::register_integration_host(integration_host);
}

/// The `(risch, manual, heurisch)` options in force for the current
/// integration (None: default).
#[pyfunction]
pub(crate) fn _integration_methods() -> (Option<bool>, Option<bool>, Option<bool>) {
    let m = host_methods::methods();
    (m.risch, m.manual, m.heurisch)
}

/// Record a held integral node as a `NonElementaryIntegral`.
#[pyfunction]
pub(crate) fn _mark_nonelementary(e: &Bound<'_, PyAny>) -> PyResult<()> {
    host_methods::mark_nonelementary(&to_expr(e)?);
    Ok(())
}

/// Whether `e` is a node recorded as a `NonElementaryIntegral`.
#[pyfunction]
pub(crate) fn _is_nonelementary(e: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(host_methods::is_nonelementary(&to_expr(e)?))
}

/// The `NonElementaryIntegral` object of a recorded node.
pub(crate) fn nonelementary_object(py: Python<'_>, e: &Expr) -> PyResult<Option<PyObject>> {
    if !matches!(e.kind(), Kind::Func(Func::Integral, _)) || !host_methods::is_nonelementary(e) {
        return Ok(None);
    }
    let class = py.import("symox.integrals.risch")?.getattr("NonElementaryIntegral")?;
    let value = Bound::new(py, PyExprNode { e: e.clone() })?;
    Ok(Some(class.call_method1("_from_native", (value,))?.unbind()))
}
