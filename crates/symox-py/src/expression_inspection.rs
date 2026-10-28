//! Scalar inspection adapters. Symbol assumptions are handled in lib.rs.
use pyo3::prelude::*;
use pyo3::types::{PySet, PyTuple, PyType};
use pyo3::exceptions::PyNotImplementedError;

use symox::expr::{Expr, Func, Kind};
use crate::{Head, expr_object, to_expr};

pub(crate) fn classify(e:&Expr, syms:&Bound<'_,PyTuple>, kind:symox::inspection::PolynomialKind)->PyResult<Option<bool>>{
    if e.has_matrix(){
        return Err(PyNotImplementedError::new_err("Polynomial inspection of typed matrix parents is not implemented"));
    }
    let generators=if syms.is_empty(){symox::traverse::free_symbols(e)}else{
        syms.iter().map(|value|to_expr(&value)).collect::<PyResult<Vec<_>>>()?
    };
    if generators.iter().any(Expr::has_matrix){
        return Err(PyNotImplementedError::new_err("Polynomial inspection requires scalar generators"));
    }
    Ok(symox::inspection::classify(e,&generators,kind))
}

enum Filter {
    Class(String), Function(Func), Functions, Python(PyObject),
    // A builtin head class: its native nodes, and instances of user classes
    // derived from it (SymPy's atoms selects subclass instances).
    Head(Func, PyObject),
}
fn class_matches(e:&Expr,name:&str)->bool{
    match name{
        "Symbol"=>e.is_symbol(),
        "Integer"=>e.is_integer_number(),
        "Rational"=>e.is_rational_number(),
        "Number"=>matches!(e.kind(),Kind::Num(_))||matches!(e.class_name(),"Infinity"|"NegativeInfinity"|"NaN"),
        "NumberSymbol"=>matches!(e.kind(),Kind::Const(c) if c.is_number_symbol()),
        _=>e.class_name()==name,
    }
}
fn function_node(e:&Expr)->bool{
    match e.kind(){
        Kind::Func(Func::Undefined(head),_)=>!head.is_internal("Tuple")&&!head.is_internal("UnevaluatedExpr")&&!head.is_internal("Order"),
        Kind::Func(Func::SymbolicDerivative|Func::Set(_)|Func::Logic(_)|Func::Max|Func::Min|Func::Subs|Func::Limit|Func::Integral|Func::Sum|Func::Product|Func::MatrixElement(_)|Func::MatrixProxy(_)|Func::MatrixTrace(_)|Func::MatrixDeterminant(_)|Func::MatrixDotProduct(..),_)=>false,
        Kind::Func(_,_)=>true,
        Kind::Piecewise(_)=>true,
        _=>false,
    }
}
fn filter(py:Python<'_>,value:&Bound<'_,PyAny>)->PyResult<Filter>{
    if let Ok(atom)=value.extract::<PyRef<'_, crate::PyExprNode>>(){
        return Ok(match atom.e.kind(){
            Kind::Func(function,_) if function_node(&atom.e)=>Filter::Function(function.clone()),
            _=>Filter::Class(atom.e.class_name().into()),
        });
    }
    if value.is(&py.import("symox")?.getattr("Function")?){return Ok(Filter::Functions);}
    if let Some(head)=crate::function_head(value)?{
        return Ok(Filter::Function(Func::Undefined(head)));
    }
    if let Ok(head)=value.extract::<PyRef<'_,Head>>(){
        return Ok(match &head.kind{
            crate::HeadKind::Func(function)=>Filter::Function(function.clone()),
            crate::HeadKind::Atom(atom)=>Filter::Class(atom.class_name().into()),
            crate::HeadKind::Add=>Filter::Class("Add".into()),
            crate::HeadKind::Mul=>Filter::Class("Mul".into()),
            crate::HeadKind::Pow=>Filter::Class("Pow".into()),
            crate::HeadKind::Derivative=>Filter::Class("Derivative".into()),
            crate::HeadKind::Piecewise=>Filter::Class("Piecewise".into()),
            crate::HeadKind::Rel(relation)=>Filter::Class(Expr::from_unevaluated_kind(Kind::Rel(*relation,Expr::zero(),Expr::zero())).class_name().into()),
        });
    }
    if value.is_callable(){
        if let Ok(name)=value.getattr("__name__").and_then(|name|name.extract::<String>()){
            // The public Expr class is a category checked by isinstance, not
            // the generic expression head of the same name.
            // Derivative has machine and expression-count representations.
            // Its public class/instance protocol must select both, rather
            // than interpreting its name as a callable Func head. A user
            // callable with the same name also retains its Python protocol.
            let category=name=="Derivative"||(name=="Expr"&&value.is_instance_of::<PyType>());
            if let Some(function)=Func::from_name(&name).filter(|_|!category){
                if value.is_instance_of::<PyType>(){return Ok(Filter::Head(function,value.clone().unbind()));}
                return Ok(Filter::Function(function));
            }
            if crate::defining_module(&value).starts_with("sympy."){
                if matches!(name.as_str(),"Number"|"NumberSymbol"|"Integer"|"Rational"|"Float"|"Symbol"|"Dummy"|"Wild"|"Zero"|"One"|"NegativeOne"|"Half"|"Add"|"Mul"|"Pow"|"Tuple"){
                    return Ok(Filter::Class(name));
                }
                if name=="Function"{return Ok(Filter::Functions);}
            }
        }
    }
    let class=if value.is_instance_of::<PyType>(){value.clone()}else{value.get_type().into_any()};
    Ok(Filter::Python(class.unbind()))
}

pub(crate) fn atoms(py:Python<'_>, e:&Expr, types:&Bound<'_,PyTuple>)->PyResult<PyObject>{
    if e.has_matrix(){
        return Err(PyNotImplementedError::new_err("Atoms of typed matrix parents are not implemented"));
    }
    let filters=types.iter().map(|value|filter(py,&value)).collect::<PyResult<Vec<_>>>()?;
    let result=PySet::empty(py)?;
    for(node,children)in symox::inspection::scalar_nodes(e){
        let mut include=filters.is_empty()&&children.is_empty();
        let mut object=None;
        for filter in &filters{
            include=match filter{
                Filter::Class(name)=>class_matches(&node,name),
                Filter::Function(function)=>matches!(node.kind(),Kind::Func(actual,_) if actual==function),
                Filter::Functions=>function_node(&node),
                Filter::Head(function,class)=>if node.is_host(){
                    let value=object.get_or_insert(expr_object(py,node.clone())?);
                    value.bind(py).is_instance(class.bind(py))?
                }else{matches!(node.kind(),Kind::Func(actual,_) if actual==function)},
                Filter::Python(class)=>{
                    let value=object.get_or_insert(expr_object(py,node.clone())?);
                    value.bind(py).is_instance(class.bind(py))?
                },
            };
            if include{break;}
        }
        if include{result.add(match object{Some(value)=>value,None=>expr_object(py,node)?})?;}
    }
    Ok(result.into_any().unbind())
}


/// Public scalar count_ops; None preserves the Python host/matrix protocol.
#[pyfunction]
#[pyo3(signature=(value,visual,kind=None))]
pub(crate) fn _native_count_ops(py:Python<'_>,value:&Bound<'_,PyAny>,visual:bool,kind:Option<&str>)->PyResult<Option<PyObject>> {
    use symox::inspection::operation_counts::{count,CountError};
    let mut expression=to_expr(value)?;
    if let Some(kind)=kind {
        expression=symox::inspection::operation_counts::presentation(kind,expression.args())
            .ok_or_else(||pyo3::exceptions::PyValueError::new_err("unsupported operation-count presentation"))?;
    }
    let counts=count(&expression).map_err(|error|match error {
        CountError::Cycle=>pyo3::exceptions::PyValueError::new_err(error.message()),
        _=>PyNotImplementedError::new_err(error.message()),
    })?;
    counts.map(|counts|if visual {
        // String case conversion is an API convention. Use this Python's
        // Unicode version, including undefined function names from user code.
        let terms=counts.labels.iter().map(|(label,number)| {
            let label:String=pyo3::types::PyString::new(py,label).call_method0("upper")?.extract()?;
            Ok(symox::ops::mul2(&Expr::num(symox::number::Number::from_bigint(number.clone())),&Expr::symbol(&label)))
        }).collect::<PyResult<Vec<_>>>()?;
        expr_object(py,symox::ops::add(&terms))
    } else {
        Ok(counts.total().into_pyobject(py)?.into_any().unbind())
    }).transpose()
}
