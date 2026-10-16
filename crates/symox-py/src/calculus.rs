//! Local complex-analysis bindings.
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use crate::{to_expr, wrap, PyExpr};

#[pyfunction]
fn residue(expr: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, x0: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let (expression, variable, point) = (to_expr(expr)?, to_expr(x)?, to_expr(x0)?);
    if !variable.is_symbol() || symox::traverse::has(&point, &variable) {
        return Err(PyValueError::new_err("residue requires a symbol and an independent point"));
    }
    symox::calculus::residue(&expression, &variable, &point).map(wrap).map_err(PyNotImplementedError::new_err)
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(residue, module)?)?;
    module.add_function(wrap_pyfunction!(_calculus, module)?)?;
    Ok(())
}

#[pyfunction]
pub(crate) fn _calculus(name:&str,body:&Bound<'_,PyAny>,limits:&Bound<'_,pyo3::types::PyTuple>)->PyResult<PyExpr> {
    use symox::{Func,Kind};
    use pyo3::types::{PyList,PyTuple};
    let kind=match name {"Integral"=>Func::Integral,"Sum"=>Func::Sum,_=>return Err(PyValueError::new_err("unknown calculus constructor"))};
    let mut body=to_expr(body)?;
    let mut parsed=Vec::new();
    for limit in limits.iter() {
        let entries=if limit.downcast::<PyTuple>().is_ok()||limit.downcast::<PyList>().is_ok() {
            let raw=limit.try_iter()?.collect::<PyResult<Vec<_>>>()?;
            if raw.iter().any(|value|value.is_none()) {
                if kind!=Func::Integral || raw.len()!=3 {return Err(PyNotImplementedError::new_err("None endpoints require an Integral (x,a,b) limit"));}
                let variable=to_expr(&raw[0])?;
                match (raw[1].is_none(),raw[2].is_none()) {
                    (true,true)=>vec![variable],
                    (true,false)=>vec![variable,to_expr(&raw[2])?],
                    (false,true)=>{body=symox::ops::neg(&body);vec![variable,to_expr(&raw[1])?]},
                    _=>unreachable!(),
                }
            } else {raw.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?}
        } else {
            let value=to_expr(&limit)?;
            match value.kind() {Kind::Func(Func::Undefined(head),args) if head.is_internal("Tuple")=>args.clone(),_=>vec![value]}
        };
        parsed.push(entries);
    }
    if parsed.is_empty() && kind==Func::Integral {
        let variables=symox::traverse::free_symbols(&body);
        if variables.len()!=1 {return Err(PyValueError::new_err("specify integration variables explicitly"));}
        parsed.push(variables);
    }
    symox::binders::held(kind,&body,parsed).map(wrap).map_err(PyValueError::new_err)
}
