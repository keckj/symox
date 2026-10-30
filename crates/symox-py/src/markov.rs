//! Finite Markov chain adapters: the kernels live in symox::stats::markov.
//! Every function returns None when the kernel declines; the caller then
//! keeps its general path.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use symox::matrix::Matrix;
use symox::stats::markov::{self, Classes, Fundamental, PathQuery};
use symox::Expr;

use crate::matrices::{wrapm_immutable, PyMatrix};
use crate::{to_expr, wrap, PyExpr};

const SINGULAR: &str = "The fundamental matrix doesn't exist.";

fn classes(py: Python<'_>, t: &Matrix) -> PyResult<Option<Classes>> {
    let mut interrupted = None;
    let result = {
        let mut check = || match py.check_signals() {
            Ok(()) => true,
            Err(e) => {
                interrupted = Some(e);
                false
            }
        };
        markov::communication_classes(t, &mut check)
    };
    match interrupted {
        Some(e) => Err(e),
        None => Ok(result),
    }
}

#[pyfunction]
fn _markov_communication_classes(py: Python<'_>, t: PyRef<'_, PyMatrix>) -> PyResult<Option<Vec<(Vec<usize>, bool, u64)>>> {
    Ok(classes(py, &t.m)?.map(|c| {
        c.classes.into_iter().zip(c.recurrent).zip(c.periods).map(|((s, r), p)| (s, r, p)).collect()
    }))
}

#[pyfunction]
fn _markov_decompose(py: Python<'_>, t: PyRef<'_, PyMatrix>) -> PyResult<Option<(Vec<usize>, PyMatrix, PyMatrix, PyMatrix)>> {
    let Some(c) = classes(py, &t.m)? else { return Ok(None) };
    let (states, a, b, cc) = markov::decompose(&t.m, &c);
    Ok(Some((states, wrapm_immutable(a), wrapm_immutable(b), wrapm_immutable(cc))))
}

#[pyfunction]
fn _markov_canonical_form(py: Python<'_>, t: PyRef<'_, PyMatrix>) -> PyResult<Option<(Vec<usize>, PyMatrix)>> {
    if t.m.rows() == 0 {
        return Ok(None);
    }
    let Some(c) = classes(py, &t.m)? else { return Ok(None) };
    let (states, a, b, cc) = markov::decompose(&t.m, &c);
    Ok(Some((states, wrapm_immutable(markov::canonical(&a, &b, &cc)))))
}

#[pyfunction]
fn _markov_stationary(t: PyRef<'_, PyMatrix>) -> Option<PyMatrix> {
    markov::stationary_rational(&t.m).map(wrapm_immutable)
}

#[pyfunction]
fn _markov_ctmc_limiting(g: PyRef<'_, PyMatrix>) -> Option<PyMatrix> {
    markov::limiting_rational(&g.m).map(wrapm_immutable)
}

fn fundamental_result(r: Option<Fundamental>) -> PyResult<Option<PyMatrix>> {
    match r {
        None => Ok(None),
        Some(Fundamental::Matrix(m)) => Ok(Some(wrapm_immutable(m))),
        Some(Fundamental::Singular) => Err(PyValueError::new_err(SINGULAR)),
    }
}

#[pyfunction]
fn _markov_fundamental_matrix(py: Python<'_>, t: PyRef<'_, PyMatrix>) -> PyResult<Option<PyMatrix>> {
    if t.m.rows() == 0 {
        return Ok(None);
    }
    let Some(c) = classes(py, &t.m)? else { return Ok(None) };
    fundamental_result(markov::fundamental_rational(&t.m, &c))
}

#[pyfunction]
fn _markov_absorbing_probabilities(py: Python<'_>, t: PyRef<'_, PyMatrix>) -> PyResult<Option<PyMatrix>> {
    if t.m.rows() == 0 {
        return Ok(None);
    }
    let Some(c) = classes(py, &t.m)? else { return Ok(None) };
    fundamental_result(markov::absorbing_rational(&t.m, &c))
}

/// The path sums of a P (`expectation=False`) or E query; `(numerator,
/// denominator)` or None when an indicator is undecided.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn _markov_path_sums(
    py: Python<'_>,
    kernels: Vec<PyRef<'_, PyMatrix>>,
    initial: usize,
    coded: Vec<Bound<'_, PyAny>>,
    labels: Vec<Bound<'_, PyAny>>,
    atoms: Vec<(Bound<'_, PyAny>, usize)>,
    gate: Option<Bound<'_, PyAny>>,
    event: Bound<'_, PyAny>,
    expectation: bool,
) -> PyResult<Option<(PyExpr, PyExpr)>> {
    let kernels: Vec<Matrix> = kernels.iter().map(|k| k.m.clone()).collect();
    if kernels.iter().any(|k| !markov::admitted(k)) {
        return Ok(None);
    }
    let convert = |v: &[Bound<'_, PyAny>]| v.iter().map(to_expr).collect::<PyResult<Vec<Expr>>>();
    let coded = convert(&coded)?;
    let labels = convert(&labels)?;
    let atoms = atoms.iter().map(|(a, s)| Ok((to_expr(a)?, *s))).collect::<PyResult<Vec<_>>>()?;
    let gate = gate.as_ref().map(to_expr).transpose()?;
    let event = to_expr(&event)?;
    let host = |e: &Expr| e.has_host() || e.has_matrix();
    if coded.iter().chain(&labels).chain(atoms.iter().map(|(a, _)| a)).chain(gate.iter()).chain(std::iter::once(&event)).any(host) {
        return Ok(None);
    }
    let q = PathQuery { kernels: &kernels, initial, coded: &coded, labels: &labels, atoms: &atoms, gate: gate.as_ref(), event: &event, expectation };
    let mut subs = |e: &Expr, pairs: &[(Expr, Expr)]| crate::checked_subs(py, e, pairs, true);
    let mut check = || py.check_signals();
    Ok(markov::enumerate_paths(&q, &mut subs, &mut check)?.map(|(n, d)| (wrap(n), wrap(d))))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_markov_communication_classes, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_decompose, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_canonical_form, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_stationary, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_ctmc_limiting, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_fundamental_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_absorbing_probabilities, m)?)?;
    m.add_function(wrap_pyfunction!(_markov_path_sums, m)?)?;
    Ok(())
}
