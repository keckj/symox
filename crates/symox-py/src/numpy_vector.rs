//! Bounded NumPy code generation for elementary float64 array kernels.
use pyo3::prelude::*;
use pyo3::types::PyList;
use symox::{Expr,Kind};
use symox::expr::{Func,Constant};
use symox::number::Number;

fn number(n:&Number)->bool{
    match n{
        Number::Int(_)=>true,
        Number::BigInt(n)=>n.bits()<=53,
        Number::Rational(n)=>n.numer().bits()<=256&&n.denom().bits()<=256,
        Number::Float(n)=>n.is_finite(),
        Number::MpFloat(n,p)=>*p<=256&&n.to_f64().is_finite(),
    }
}

pub(crate) fn supports(e:&Expr)->bool{
    fn walk(e:&Expr,memo:&mut rustc_hash::FxHashMap<usize,bool>)->bool{
        if let Some(value)=memo.get(&e.dag_id()){return *value;}
        let value=match e.kind(){
            Kind::Num(n)=>number(n),Kind::Sym(_)=>true,
            Kind::Const(Constant::Pi|Constant::E)=>true,
            Kind::Add(c,terms)=>number(c)&&terms.iter().all(|(term,c)|number(c)&&walk(term,memo)),
            Kind::Mul(c,factors)=>number(c)&&factors.iter().all(|(b,p)|walk(b,memo)&&walk(p,memo)),
            Kind::Pow(b,p)=>walk(b,memo)&&walk(p,memo),
            Kind::Func(Func::Sin|Func::Cos|Func::Exp|Func::Log,args)=>args.iter().all(|e|walk(e,memo)),
            _=>false,
        };
        memo.insert(e.dag_id(),value);value
    }
    !e.has_held_arithmetic()&&!e.has_matrix()&&!e.has_wild_properties()
        &&walk(e,&mut rustc_hash::FxHashMap::default())
        &&symox::traverse::node_count(e)<=4096
}

pub(crate) fn dispatcher(py:Python<'_>,arguments:&[Expr],expression:&Expr)->PyResult<Py<PyAny>>{
    let arguments=PyList::new(py,arguments.iter().cloned().map(crate::wrap))?;
    Ok(py.import("symox._numpy_vector")?.getattr("compile")?
        .call1((arguments,crate::wrap(expression.clone())))?.unbind())
}
