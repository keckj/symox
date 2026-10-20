//! Recognize bounded exact monomials before ordinary scalar Add dispatch.
use pyo3::prelude::*;
use pyo3::types::{PyBytes,PyTuple};
use symox::{Expr,Kind};
use symox::number::Number;

fn coefficient(n:&Number)->bool{
    match n{
        Number::Int(_)=>true,
        Number::BigInt(n)=>n.bits()<=256,
        Number::Rational(n)=>n.numer().bits()<=256&&n.denom().bits()<=256,
        _=>false,
    }
}
fn symbol(e:&Expr)->bool{
    e.as_sym().is_some_and(|s|!s.is_wild()&&!s.is_dummy()&&s.name().len()<=256)
}
fn exponent(e:&Expr)->bool{
    match e.as_number(){
        Some(Number::Int(n))=>*n>=0,
        Some(Number::BigInt(n))=>n.sign()!=num_bigint::Sign::Minus&&n.bits()<=256,
        _=>false,
    }
}
fn monomial(e:&Expr)->bool{
    if !e.is_canonical()||e.has_nonfinite()||e.has_wild_properties(){return false;}
    match e.kind(){
        Kind::Num(n)=>coefficient(n),Kind::Sym(_)=>symbol(e),
        Kind::Pow(b,p)=>symbol(b)&&exponent(p),
        Kind::Mul(c,factors)=>coefficient(c)&&factors.len()<=8
            &&factors.iter().all(|(b,p)|symbol(b)&&exponent(p)),
        _=>false,
    }
}

#[pyfunction]
pub(crate) fn _exact_expr_args(args:&Bound<'_,PyTuple>)->bool{
    let kind=args.py().get_type::<crate::PyExprNode>();
    args.iter().all(|value|value.get_type().is(&kind))
}

/// An injective identity key, with no Python hashing or conversion hooks.
/// Owners are retained by the bounded Python cache until this key is evicted.
#[pyfunction]
pub(crate) fn _exact_expr_key<'py>(args:&Bound<'py,PyTuple>)->Option<Bound<'py,PyBytes>>{
    if !(2..=1024).contains(&args.len()){return None;}
    let kind=args.py().get_type::<crate::PyExprNode>();
    let mut key=Vec::with_capacity(args.len()*std::mem::size_of::<usize>());
    for value in args.iter(){
        if !value.get_type().is(&kind){return None;}
        key.extend_from_slice(&(value.as_ptr() as usize).to_ne_bytes());
    }
    Some(PyBytes::new(args.py(),&key))
}

#[pyfunction]
pub(crate) fn _monomial_add(args:&Bound<'_,PyTuple>)->PyResult<Option<crate::PyExpr>>{
    if !(2..=1024).contains(&args.len()){return Ok(None);}
    let mut values=Vec::with_capacity(args.len());
    let kind=args.py().get_type::<crate::PyExprNode>();
    for value in args.iter(){
        if !value.get_type().is(&kind){return Ok(None);}
        let value=value.extract::<PyRef<'_,crate::PyExprNode>>()?;
        if !monomial(&value.e){return Ok(None);}
        values.push(value.e.clone());
    }
    Ok(Some(crate::wrap(symox::ops::add(&values))))
}
