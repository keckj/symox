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
    module.add_function(wrap_pyfunction!(product, module)?)?;
    module.add_function(wrap_pyfunction!(_is_convergent, module)?)?;
    module.add_function(wrap_pyfunction!(_is_absolutely_convergent, module)?)?;
    module.add_function(wrap_pyfunction!(_limits_index, module)?)?;
    module.add_function(wrap_pyfunction!(_reverse_order, module)?)?;
    module.add_function(wrap_pyfunction!(_change_index, module)?)?;
    module.add_function(wrap_pyfunction!(_reorder_limit, module)?)?;
    module.add_function(wrap_pyfunction!(_limits_property, module)?)?;
    Ok(())
}

#[pyfunction]
pub(crate) fn _calculus(name:&str,body:&Bound<'_,PyAny>,limits:&Bound<'_,pyo3::types::PyTuple>)->PyResult<PyExpr> {
    use symox::{Func,Kind};
    use pyo3::types::{PyList,PyTuple};
    let kind=match name {"Integral"=>Func::Integral,"Sum"=>Func::Sum,"Product"=>Func::Product,_=>return Err(PyValueError::new_err("unknown calculus constructor"))};
    let mut body=to_expr(body)?;
    // SymPy's ExprWithLimits returns nan before looking at the limits
    if body.is_constant(symox::expr::Constant::NaN) {return Ok(wrap(body));}
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

/// `product(f, *limits)`: `Product(f, *limits).doit(deep=False)`.
#[pyfunction]
#[pyo3(signature = (f, *limits))]
fn product(f:&Bound<'_,PyAny>,limits:&Bound<'_,pyo3::types::PyTuple>)->PyResult<PyExpr> {
    let held=_calculus("Product",f,limits)?;
    match held.e.kind() {
        symox::Kind::Func(symox::Func::Product,args)=>Ok(wrap(symox::products::doit(args,false))),
        _=>Ok(held),
    }
}

/// The term and single limit of a Sum or Product node.
fn one_limit(e:&symox::Expr,heads:&[symox::Func],attribute:&str)->PyResult<(symox::Func,symox::Expr,Vec<symox::Expr>)> {
    use symox::Kind;
    match e.kind() {
        Kind::Func(f,args) if heads.contains(f) && args.len()==2=>Ok((f.clone(),args[0].clone(),args[1].args())),
        Kind::Func(f,_) if heads.contains(f)=>Err(PyNotImplementedError::new_err("convergence of multiple limits is not implemented")),
        _=>Err(pyo3::exceptions::PyAttributeError::new_err(format!("'Expr' object has no attribute '{attribute}'"))),
    }
}

/// `Sum.is_convergent()` and `Product.is_convergent()` for one-limit nodes.
#[pyfunction]
fn _is_convergent(e:&Bound<'_,PyAny>)->PyResult<PyExpr> {
    use symox::Func;
    let (f,term,limit)=one_limit(&to_expr(e)?,&[Func::Sum,Func::Product],"is_convergent")?;
    let check=if f==Func::Sum {symox::convergence::sum_is_convergent} else {symox::convergence::is_convergent};
    check(&term,&limit[0],&limit[1],&limit[2]).map(|b|wrap(symox::Expr::boolean(b))).map_err(PyNotImplementedError::new_err)
}

/// `Sum.is_absolutely_convergent()`.
#[pyfunction]
fn _is_absolutely_convergent(e:&Bound<'_,PyAny>)->PyResult<PyExpr> {
    let (_,term,limit)=one_limit(&to_expr(e)?,&[symox::Func::Sum],"is_absolutely_convergent")?;
    symox::convergence::sum_is_absolutely_convergent(&term,&limit[0],&limit[1],&limit[2])
        .map(|b|wrap(symox::Expr::boolean(b))).map_err(PyNotImplementedError::new_err)
}

fn limits_error(e:symox::intlimits::LimitsError)->PyErr {
    match e {
        symox::intlimits::LimitsError::Value(m)=>PyValueError::new_err(m),
        symox::intlimits::LimitsError::Reorder(m)=>PyNotImplementedError::new_err(m),
    }
}

/// `expr.index(x)` of a Sum or Product.
#[pyfunction]
fn _limits_index(e:&Bound<'_,PyAny>,x:&Bound<'_,PyAny>)->PyResult<usize> {
    symox::intlimits::index(&to_expr(e)?,&to_expr(x)?).map_err(limits_error)
}

/// `expr.reverse_order(*indices)` with the indices already resolved.
#[pyfunction]
fn _reverse_order(e:&Bound<'_,PyAny>,indices:Vec<usize>)->PyResult<PyExpr> {
    symox::intlimits::reverse_order(&to_expr(e)?,&indices).map(wrap).map_err(limits_error)
}

/// `expr.change_index(var, trafo, newvar=None)`.
#[pyfunction]
#[pyo3(signature = (e, var, trafo, newvar=None))]
fn _change_index(e:&Bound<'_,PyAny>,var:&Bound<'_,PyAny>,trafo:&Bound<'_,PyAny>,newvar:Option<&Bound<'_,PyAny>>)->PyResult<PyExpr> {
    let newvar=newvar.filter(|v|!v.is_none()).map(to_expr).transpose()?;
    symox::intlimits::change_index(&to_expr(e)?,&to_expr(var)?,&to_expr(trafo)?,newvar.as_ref()).map(wrap).map_err(limits_error)
}

/// `expr.reorder_limit(x, y)`; a failed interchange raises NotImplementedError.
#[pyfunction]
fn _reorder_limit(e:&Bound<'_,PyAny>,x:usize,y:usize)->PyResult<PyExpr> {
    symox::intlimits::reorder_limit(&to_expr(e)?,x,y).map(wrap).map_err(limits_error)
}

/// `has_finite_limits`, `has_reversed_limits` or `has_empty_sequence`.
#[pyfunction]
fn _limits_property(e:&Bound<'_,PyAny>,name:&str)->PyResult<Option<bool>> {
    let e=to_expr(e)?;
    let value=match name {
        "has_finite_limits"=>symox::intlimits::has_finite_limits(&e),
        "has_reversed_limits"=>symox::intlimits::has_reversed_limits(&e),
        "has_empty_sequence"=>symox::intlimits::has_empty_sequence(&e),
        _=>return Err(PyValueError::new_err("unknown limits property")),
    };
    value.map_err(|_|pyo3::exceptions::PyAttributeError::new_err(format!("'Expr' object has no attribute '{name}'")))
}
