//! `_piecewise_fold` and `_fold_piecewise_function`: the native public
//! `piecewise_fold` and the `piecewise_fold` step of the `Sum`, `Product`
//! and `Integral` constructors. `None` means the kernel declined and the
//! Python code runs.
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use symox::simptools::piecewise_fold_public::{self as pf, FoldHost};
use symox::{Expr, Kind};

enum Failure {
    Python(PyErr),
    Decline,
}

struct Host<'py> {
    py: Python<'py>,
    facts: Option<symox::assumption_facts::FactSet>,
}

impl<'py> Host<'py> {
    fn new(py: Python<'py>) -> Self {
        Host { py, facts: None }
    }
}

impl FoldHost for Host<'_> {
    type Error = Failure;

    fn rebuild(&mut self, expr: &Expr, args: Vec<Expr>) -> Result<Expr, Failure> {
        let run = || -> PyResult<Option<Expr>> {
            let func = crate::expr_object(self.py, expr.clone())?.into_bound(self.py).getattr("func")?;
            let args = args.into_iter().map(|a| crate::expr_object(self.py, a)).collect::<PyResult<Vec<_>>>()?;
            let rebuilt = func.call1(PyTuple::new(self.py, args)?)?;
            if !rebuilt.is_instance_of::<crate::PyExprNode>() {
                return Ok(None);
            }
            Ok(crate::to_expr(&rebuilt).ok())
        };
        match run() {
            Ok(Some(e)) => Ok(e),
            Ok(None) => Err(Failure::Decline),
            Err(e) => Err(Failure::Python(e)),
        }
    }

    fn ite_condition(&mut self, condition: &Expr) -> Result<Expr, Failure> {
        let run = || -> PyResult<Option<Expr>> {
            let simplify_logic = crate::loaded_symox(self.py)?.getattr("simplify_logic")?;
            let nnf = crate::expr_object(self.py, condition.clone())?.into_bound(self.py).call_method0("to_nnf")?;
            let kwargs = pyo3::types::PyDict::new(self.py);
            kwargs.set_item("form", "cnf")?;
            let c = simplify_logic.call((nnf,), Some(&kwargs))?;
            if !c.is_instance_of::<crate::PyExprNode>() {
                return Ok(None);
            }
            Ok(crate::to_expr(&c).ok())
        };
        match run() {
            Ok(Some(e)) => Ok(e),
            Ok(None) => Err(Failure::Decline),
            Err(e) => Err(Failure::Python(e)),
        }
    }

    fn decline(&mut self) -> Failure {
        Failure::Decline
    }

    fn dummy(&mut self) -> Expr {
        // `Dummy()`: the facts of a plain Symbol, the shared counter
        let facts = *self.facts.get_or_insert_with(|| symox::Symbol::new("").given_facts());
        match symox::Symbol::anonymous_dummy_with_facts(facts) {
            Ok(s) => Expr::from_kind(Kind::Sym(s)),
            Err(_) => unreachable!("the facts of a plain symbol are consistent"),
        }
    }
}

/// The input object itself when it holds no Piecewise (the Python code
/// returns it unchanged), the folded expression, or None when declined.
fn finish(py: Python<'_>, result: Result<Expr, Failure>) -> PyResult<PyObject> {
    match result {
        Ok(e) => crate::expr_object(py, e),
        Err(Failure::Decline) => Ok(py.None()),
        Err(Failure::Python(e)) => Err(e),
    }
}

/// `piecewise_fold(expr, evaluate)` with `evaluate` None or a truth value.
#[pyfunction]
#[pyo3(signature = (expr, evaluate))]
fn _piecewise_fold(py: Python<'_>, expr: &Bound<'_, PyAny>, evaluate: Option<bool>) -> PyResult<PyObject> {
    let native = crate::to_expr(expr)?;
    if !pf::admitted(&native) {
        return Ok(py.None());
    }
    if !pf::contains_piecewise(&native) {
        return Ok(expr.clone().unbind());
    }
    finish(py, pf::piecewise_fold(&mut Host::new(py), &native, evaluate))
}

/// `_fold_piecewise_function(function, limits)` with the limit variables
/// (the free symbols of `function` when there are no limits).
#[pyfunction]
fn _fold_piecewise_function(py: Python<'_>, function: &Bound<'_, PyAny>, variables: Vec<Bound<'_, PyAny>>) -> PyResult<PyObject> {
    let native = crate::to_expr(function)?;
    if !pf::admitted(&native) {
        return Ok(py.None());
    }
    if !pf::contains_piecewise(&native) {
        return Ok(function.clone().unbind());
    }
    let variables = variables.iter().map(crate::to_expr).collect::<PyResult<Vec<_>>>()?;
    finish(py, pf::fold_limits_function(&mut Host::new(py), &native, &variables))
}

/// `expr.has(Piecewise)` of a native expression, None when it holds host
/// objects or matrices (the Python query decides).
#[pyfunction]
fn _has_piecewise(expr: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    let expr = crate::to_expr(expr)?;
    if expr.has_host() || expr.has_matrix() {
        return Ok(None);
    }
    Ok(Some(pf::contains_piecewise(&expr)))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_piecewise_fold, m)?)?;
    m.add_function(wrap_pyfunction!(_fold_piecewise_function, m)?)?;
    m.add_function(wrap_pyfunction!(_has_piecewise, m)?)
}
