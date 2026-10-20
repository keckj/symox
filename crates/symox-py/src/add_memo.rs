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

/// Plain node objects: exactly the base node type or the nominal SymPy
/// class native wrapping gives that node, with no host node inside. User
/// subclasses, host trees and other wrappers keep full constructor dispatch.
/// A nominal class is only ever allocated for its own nodes, so each class
/// needs a single nominal lookup per argument tuple.
struct PlainNodes<'py>{
    base:Bound<'py,pyo3::types::PyType>,
    verified:Vec<*mut pyo3::ffi::PyTypeObject>,
}

impl<'py> PlainNodes<'py>{
    fn new(py:Python<'py>)->Self{
        PlainNodes{base:py.get_type::<crate::PyExprNode>(),verified:Vec::new()}
    }

    fn node<'a>(&mut self,value:&'a Bound<'py,PyAny>)->Option<&'a Expr>{
        let pointer=unsafe{pyo3::ffi::Py_TYPE(value.as_ptr())};
        if pointer==self.base.as_type_ptr()||self.verified.contains(&pointer){
            // SAFETY: the base node type and the verified nominal classes
            // are PyExprNode types.
            let e=&unsafe{value.downcast_unchecked::<crate::PyExprNode>()}.get().e;
            return (!e.has_host()).then_some(e);
        }
        let node=value.downcast::<crate::PyExprNode>().ok()?;
        let e=&node.get().e;
        if e.has_host(){return None;}
        let class=value.get_type();
        if crate::node_type(value.py(),e).is_some_and(|nominal|class.is(&nominal)){
            if self.verified.len()<8{self.verified.push(pointer);}
            Some(e)
        } else {None}
    }
}

#[pyfunction]
pub(crate) fn _exact_expr_args(args:&Bound<'_,PyTuple>)->bool{
    let mut plain=PlainNodes::new(args.py());
    args.iter().all(|value|plain.node(&value).is_some())
}

/// An injective identity key, with no Python hashing or conversion hooks.
/// Owners are retained by the bounded Python cache until this key is evicted.
#[pyfunction]
pub(crate) fn _exact_expr_key<'py>(args:&Bound<'py,PyTuple>)->Option<Bound<'py,PyBytes>>{
    if !(2..=1024).contains(&args.len()){return None;}
    let mut key=Vec::with_capacity(args.len()*std::mem::size_of::<usize>());
    let mut plain=PlainNodes::new(args.py());
    for value in args.iter(){
        plain.node(&value)?;
        key.extend_from_slice(&(value.as_ptr() as usize).to_ne_bytes());
    }
    Some(PyBytes::new(args.py(),&key))
}

#[pyfunction]
pub(crate) fn _monomial_add(args:&Bound<'_,PyTuple>)->PyResult<Option<crate::PyExpr>>{
    if !(2..=1024).contains(&args.len()){return Ok(None);}
    let mut values=Vec::with_capacity(args.len());
    let mut plain=PlainNodes::new(args.py());
    for value in args.iter(){
        let Some(value)=plain.node(&value) else {return Ok(None);};
        if !monomial(value){return Ok(None);}
        values.push(value.clone());
    }
    Ok(Some(crate::wrap(symox::ops::add(&values))))
}
