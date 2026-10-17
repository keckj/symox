//! Thin bindings for the native assumptions engine.
use crate::{replacement_expr, to_expr, wrap, PyExpr};
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::PyTuple,
};
use symox::{
    ask::{self, AskError, AskOptions, Pred},
    expr::Expr,
};

pub(crate) fn error(e: AskError) -> PyErr {
    match e {
        AskError::Value(m) => PyValueError::new_err(m),
        AskError::Type(m) => PyTypeError::new_err(m),
    }
}

#[pyfunction]
#[pyo3(signature = (name, undefined=false))]
fn _predicate(name: &str, undefined: bool) -> PyResult<PyExpr> {
    let p = if undefined {
        Pred::undefined(name)
    } else {
        Pred::from_name(name).ok_or_else(|| PyValueError::new_err("unknown predicate"))?
    };
    Ok(wrap(ask::predicate(p)))
}
#[pyfunction]
fn _predicate_names() -> Vec<&'static str> {
    Pred::all().map(Pred::name).collect()
}

#[pyfunction]
fn _applied_predicate(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    ask::applied_from_args(
        args.iter()
            .map(|arg| replacement_expr(&arg))
            .collect::<PyResult<_>>()?,
    )
    .map(wrap)
    .map_err(error)
}

fn context(value: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
    value.try_iter()?.map(|arg| to_expr(&arg?)).collect()
}
#[pyfunction]
fn _ask(
    proposition: &Bound<'_, PyAny>,
    assumptions: &Bound<'_, PyAny>,
    ctx: &Bound<'_, PyAny>,
) -> PyResult<Option<bool>> {
    ask::ask_with(
        &to_expr(proposition)?,
        &to_expr(assumptions)?,
        &AskOptions {
            context: context(ctx)?,
            external: None,
        },
    )
    .map_err(error)
}
#[pyfunction]
fn _refine(
    expr: &Bound<'_, PyAny>,
    assumptions: &Bound<'_, PyAny>,
    ctx: &Bound<'_, PyAny>,
) -> PyResult<PyExpr> {
    ask::refine_with(
        &to_expr(expr)?,
        &to_expr(assumptions)?,
        &AskOptions {
            context: context(ctx)?,
            external: None,
        },
    )
    .map(wrap)
    .map_err(error)
}
pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_predicate, m)?)?;
    m.add_function(wrap_pyfunction!(_predicate_names, m)?)?;
    m.add_function(wrap_pyfunction!(_applied_predicate, m)?)?;
    m.add_function(wrap_pyfunction!(_ask, m)?)?;
    m.add_function(wrap_pyfunction!(_refine, m)?)?;
    Ok(())
}
