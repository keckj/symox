//! Python bindings for the native characteristic-equation solutions of
//! linear ODEs (`symox::ode::constcoeff`). Each function returns `None`
//! when the kernel declines, and the caller then runs the reference body.

use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};

use symox::ode::constcoeff::{self, CollectTerm};
use symox::ode::{undetermined, variation};
use symox::Expr;

use crate::{to_expr, wrap as w};

fn coefficient_list(coeffs: &Bound<'_, PyAny>) -> PyResult<Vec<(u32, Expr)>> {
    let mut out = Vec::new();
    for item in coeffs.try_iter()? {
        let (i, c): (u32, Bound<'_, PyAny>) = item?.extract()?;
        out.push((i, to_expr(&c)?));
    }
    Ok(out)
}

fn expr_list(values: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
    values.try_iter()?.map(|v| to_expr(&v?)).collect()
}

fn terms_to_py<'py>(py: Python<'py>, terms: &[CollectTerm]) -> PyResult<Bound<'py, PyList>> {
    let out = PyList::empty(py);
    for (i, a, b) in terms {
        let i = i.into_pyobject(py)?.into_any();
        let a = w(a.clone()).into_pyobject(py)?.into_any();
        let b = w(b.clone()).into_pyobject(py)?.into_any();
        out.append(PyTuple::new(py, [i, a, b])?)?;
    }
    Ok(out)
}

/// `_get_const_characteristic_eq_sols`: `coeffs` are the `(i, r[i])` pairs
/// of the integer keys `i >= 0`. Returns `(gensols, collectterms)` or `None`.
#[pyfunction]
fn ode_const_char_sols(py: Python<'_>, coeffs: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, order: usize) -> PyResult<PyObject> {
    let coeffs = coefficient_list(coeffs)?;
    let x = to_expr(x)?;
    Ok(match constcoeff::const_characteristic_sols(&coeffs, &x, order) {
        Some(s) => {
            let gensols = PyList::new(py, s.gensols.into_iter().map(w))?;
            PyTuple::new(py, [gensols.into_any(), terms_to_py(py, &s.collectterms)?.into_any()])?.into_any().unbind()
        }
        None => py.None(),
    })
}

/// `_get_simplified_sol` on the right-hand side `rhs`: `collectterms` is a
/// sequence of `(i, reroot, imroot)`. Returns the simplified rhs or `None`.
#[pyfunction]
fn ode_const_simplified_sol(py: Python<'_>, rhs: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, collectterms: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let rhs = to_expr(rhs)?;
    let x = to_expr(x)?;
    let mut terms = Vec::new();
    for item in collectterms.try_iter()? {
        let (i, a, b): (u32, Bound<'_, PyAny>, Bound<'_, PyAny>) = item?.extract()?;
        terms.push((i, to_expr(&a)?, to_expr(&b)?));
    }
    Ok(match constcoeff::simplified_solution(&rhs, &x, &terms) {
        Some(e) => w(e).into_pyobject(py)?.into_any().unbind(),
        None => py.None(),
    })
}

/// `_get_euler_characteristic_eq_sols`: `coeffs` are the `(i, match[i])`
/// pairs of the integer keys `i >= 0`, `constants` the `2*degree` numbered
/// constants. Returns `(gsol_rhs, gensols)` or `None`.
#[pyfunction]
fn ode_euler_char_sols(py: Python<'_>, coeffs: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, constants: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let coeffs = coefficient_list(coeffs)?;
    let x = to_expr(x)?;
    let constants = expr_list(constants)?;
    Ok(match constcoeff::euler_characteristic_sols(&coeffs, &x, &constants) {
        Some((rhs, gensols)) => {
            let gensols = PyList::new(py, gensols.into_iter().map(w))?;
            PyTuple::new(py, [w(rhs).into_pyobject(py)?.into_any(), gensols.into_any()])?.into_any().unbind()
        }
        None => py.None(),
    })
}

/// `_undetermined_coefficients_match(forcing, x, f(x), eq_homogeneous)['trialset']`
/// (`homogeneous` is `None` for `eq_homogeneous is S.Zero`). Returns the
/// list of trial terms or `None`.
#[pyfunction]
#[pyo3(signature = (forcing, x, fx, homogeneous=None))]
fn ode_uc_trialset(py: Python<'_>, forcing: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, fx: &Bound<'_, PyAny>, homogeneous: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
    let forcing = to_expr(forcing)?;
    let x = to_expr(x)?;
    let fx = to_expr(fx)?;
    let homogeneous = homogeneous.map(to_expr).transpose()?;
    Ok(match undetermined::trial_set(&forcing, &x, &fx, homogeneous.as_ref()) {
        Some(terms) => PyList::new(py, terms.into_iter().map(w))?.into_any().unbind(),
        None => py.None(),
    })
}

/// `_solve_undetermined_coefficients`: the right-hand side
/// `gsol.rhs + psol` or `None`.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn ode_uc_solve(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    order: usize,
    gensols_len: usize,
    gsol_rhs: &Bound<'_, PyAny>,
    trialset: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let eq = to_expr(eq)?;
    let fx = to_expr(fx)?;
    let x = to_expr(x)?;
    let gsol_rhs = to_expr(gsol_rhs)?;
    let trialset = expr_list(trialset)?;
    Ok(match undetermined::solve_undetermined(&eq, &fx, &x, order, gensols_len, &gsol_rhs, &trialset) {
        Some(e) => w(e).into_pyobject(py)?.into_any().unbind(),
        None => py.None(),
    })
}

enum Callback {
    Py(PyErr),
    Decline,
}

fn callback<'py>(f: &Bound<'py, PyAny>) -> impl FnMut(&Expr) -> Result<Expr, Callback> + 'py {
    let f = f.clone();
    move |e: &Expr| {
        let out = f.call1((w(e.clone()),)).map_err(Callback::Py)?;
        to_expr(&out).map_err(|_| Callback::Decline)
    }
}

/// `_solve_variation_of_parameters`: the right-hand side
/// `homogen_sol.rhs + psol` or `None`. `simplify_wronskian` and
/// `simplify_psol` are the reference's simplifications, called only when
/// `simplify_flag` is set.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn ode_vop_solve(
    py: Python<'_>,
    roots: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    forcing: &Bound<'_, PyAny>,
    lead: &Bound<'_, PyAny>,
    order: usize,
    homogen_rhs: &Bound<'_, PyAny>,
    simplify_flag: bool,
    simplify_wronskian: &Bound<'_, PyAny>,
    simplify_psol: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let roots = expr_list(roots)?;
    let x = to_expr(x)?;
    let forcing = to_expr(forcing)?;
    let lead = to_expr(lead)?;
    let homogen_rhs = to_expr(homogen_rhs)?;
    let mut sw = callback(simplify_wronskian);
    let mut sp = callback(simplify_psol);
    match variation::variation_of_parameters(&roots, &x, &forcing, &lead, order, &homogen_rhs, simplify_flag, &mut sw, &mut sp) {
        Ok(Some(e)) => Ok(w(e).into_pyobject(py)?.into_any().unbind()),
        Ok(None) | Err(Callback::Decline) => Ok(py.None()),
        Err(Callback::Py(err)) => Err(err),
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ode_uc_trialset, m)?)?;
    m.add_function(wrap_pyfunction!(ode_uc_solve, m)?)?;
    m.add_function(wrap_pyfunction!(ode_vop_solve, m)?)?;
    m.add_function(wrap_pyfunction!(ode_const_char_sols, m)?)?;
    m.add_function(wrap_pyfunction!(ode_const_simplified_sol, m)?)?;
    m.add_function(wrap_pyfunction!(ode_euler_char_sols, m)?)?;
    Ok(())
}
