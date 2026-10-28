//! Canonical exact field (QQ, GF(p), QQ<a>) triangular solve, retaining
//! actual scratch/output publication.
use pyo3::prelude::*;
use pyo3::types::{PyInt,PyList,PyTuple,PyDict,PyString};
use symox::poly::{lu_solve,domains::ExactDomain};
use super::{sequence,entries,rationals,scalar,Scalar,element_object,field_domain,field_values,field_object,with_field_scope,FieldKind,FieldPayload};
use crate::poly::domains::PyDomainElement;

/// Read actual caller lookups without evaluating late globals too early.
/// In particular missing/replaced reversed/error names must not precede the
/// established forward substitution and consistency failure. Real frame
/// dictionaries have immutable C getters; no user module/attribute hooks run.
#[pyfunction]
fn _ddm_ilu_solve_dependencies(py:Python<'_>)->PyResult<Option<PyObject>> {
    let frame=unsafe {pyo3::ffi::PyEval_GetFrame()};
    if frame.is_null(){return Ok(None);}
    let frame=unsafe {Bound::from_borrowed_ptr(py,frame.cast())};
    let globals=frame.getattr(pyo3::intern!(py,"f_globals"))?;
    let builtins=frame.getattr(pyo3::intern!(py,"f_builtins"))?;
    let (Ok(globals),Ok(builtins))=(globals.downcast_exact::<PyDict>(),builtins.downcast_exact::<PyDict>()) else{return Ok(None)};
    // A non-string dictionary key could execute user equality during lookup.
    if globals.iter().chain(builtins.iter()).any(|(key,_)|!key.is_exact_instance_of::<PyString>()){return Ok(None);}
    let mut values=Vec::new();
    for name in ["range","reversed","DMNonInvertibleMatrixError"] {
        let value=if let Some(value)=globals.get_item(name)?{Some(value)}else{builtins.get_item(name)?};
        values.push(value.map_or_else(||py.None(),Bound::unbind));
    }
    Ok(Some(PyTuple::new(py,values)?.into_any().unbind()))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)] // Explicit real prefix objects, never cached mutable proof state.
fn _ddm_ilu_solve(x:&Bound<'_,PyAny>,l:&Bound<'_,PyAny>,u:&Bound<'_,PyAny>,
    b:&Bound<'_,PyAny>,y:&Bound<'_,PyAny>,dimensions:&Bound<'_,PyAny>,
    stock:&Bound<'_,PyAny>,exception:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<bool> {
    py.check_signals()?;
    let Ok(dimensions)=dimensions.downcast_exact::<PyTuple>() else{return Ok(false)};
    if dimensions.len()!=3{return Ok(false);}
    let mut sizes=Vec::new();
    for v in dimensions.iter(){
        if !v.is_exact_instance_of::<PyInt>(){return Ok(false);}
        let Ok(v)=v.extract::<usize>() else{return Ok(false)};sizes.push(v);
    }
    let (m,n,o)=(sizes[0],sizes[1],sizes[2]);
    if m<n{return Ok(false);}
    let stock=if stock.is_none(){None}else{Some(stock)};
    let Some(xrows)=sequence(x,stock,true) else{return Ok(false)};
    let Some(yrows)=sequence(y,None,true) else{return Ok(false)};
    let Some(lrows)=sequence(l,stock,false) else{return Ok(false)};
    let Some(urows)=sequence(u,stock,false) else{return Ok(false)};
    let Some(brows)=sequence(b,stock,false) else{return Ok(false)};
    let mut destinations=rustc_hash::FxHashSet::default();
    for row in xrows.iter().chain(&yrows){
        if !row.is_exact_instance_of::<PyList>()||!destinations.insert(row.as_ptr() as usize){return Ok(false);}
    }
    if lrows.iter().chain(&urows).chain(&brows).any(|row|destinations.contains(&(row.as_ptr() as usize))){return Ok(false);}
    let Some(xcells)=entries(&xrows,true) else{return Ok(false)};
    let Some(ycells)=entries(&yrows,true) else{return Ok(false)};
    let Some(lcells)=entries(&lrows,false) else{return Ok(false)};
    let Some(ucells)=entries(&urows,false) else{return Ok(false)};
    let Some(bcells)=entries(&brows,false) else{return Ok(false)};
    if xcells.len()!=n||ycells.len()!=m||lcells.len()!=m||ucells.len()!=m||bcells.len()!=m
        ||xcells.iter().chain(&ycells).chain(&bcells).any(|row|row.len()!=o)
        ||lcells.iter().any(|row|row.len()!=m)||ucells.iter().any(|row|row.len()!=n){return Ok(false);}
    // User destructors in overwritten arbitrary placeholders can mutate source
    // rows/protocols. None and canonical immutable QQ cells cannot do that.
    // Frozen native DomainElement payloads have no user destructor either.
    if xcells.iter().flatten().any(|v|!v.is_none()&&!v.is_exact_instance_of::<PyDomainElement>())
        ||ycells.iter().flatten().any(|v|!v.is_none()){return Ok(false);}
    let xr=xrows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
    let yr=yrows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
    if let Some(first)=ucells.first().and_then(|row|row.first()).filter(|v|scalar(v)!=Some(Scalar::QQ)) {
        // QQ_I division uses the Python Gaussian wrapper: not admitted here.
        let Some((domain,kind))=field_domain(first,false) else{return Ok(false)};
        return match kind {
            FieldKind::Finite=>field_solve::<symox::poly::f5b::Fp>(&lcells,&ucells,&bcells,o,&xr,&yr,&domain,exception,py),
            FieldKind::Algebraic=>field_solve::<symox::poly::algebraic::AN>(&lcells,&ucells,&bcells,o,&xr,&yr,&domain,exception,py),
        };
    }
    if xcells.iter().flatten().any(|v|!v.is_none()&&scalar(v)!=Some(Scalar::QQ)){return Ok(false);}
    let Some(l)=rationals(&lcells,py)? else{return Ok(false)};
    let Some(u)=rationals(&ucells,py)? else{return Ok(false)};
    let Some(b)=rationals(&bcells,py)? else{return Ok(false)};
    match lu_solve::solve(&l,&u,&b,o,&mut ||py.check_signals(),&mut |event| {
        match event {
            lu_solve::Event::Y(i,k,value,from_b)=>{
                if from_b {yr[i].set_item(k,&bcells[i][k])}
                else {yr[i].set_item(k,element_object(py,value.clone(),ExactDomain::QQ)?)}
            },
            lu_solve::Event::X(i,k,value)=>xr[i].set_item(k,element_object(py,value.clone(),ExactDomain::QQ)?),
        }
    }) {
        Ok(done)=>Ok(done),
        Err(lu_solve::Failure::Interrupted(error))=>Err(error),
        Err(lu_solve::Failure::NonInvertible)=>Err(PyErr::from_value(exception.call0()?)),
    }
}

#[allow(clippy::too_many_arguments)]
fn field_solve<C:FieldPayload>(lcells:&[Vec<Bound<'_,PyAny>>],ucells:&[Vec<Bound<'_,PyAny>>],bcells:&[Vec<Bound<'_,PyAny>>],
    o:usize,xr:&[Bound<'_,PyList>],yr:&[Bound<'_,PyList>],domain:&symox::poly::domains::ExactDomain,
    exception:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<bool> {
    let Some(l)=field_values::<C>(lcells,domain,py)? else{return Ok(false)};
    let Some(u)=field_values::<C>(ucells,domain,py)? else{return Ok(false)};
    let Some(b)=field_values::<C>(bcells,domain,py)? else{return Ok(false)};
    let outcome=with_field_scope(domain,||Ok(lu_solve::solve(&l,&u,&b,o,&mut ||py.check_signals(),&mut |event| {
        match event {
            lu_solve::Event::Y(i,k,value,from_b)=>{
                if from_b {yr[i].set_item(k,&bcells[i][k])}
                else {yr[i].set_item(k,field_object(py,value.clone(),domain)?)}
            },
            lu_solve::Event::X(i,k,value)=>xr[i].set_item(k,field_object(py,value.clone(),domain)?),
        }
    })))?;
    match outcome {
        None=>Ok(false),
        Some(Ok(done))=>Ok(done),
        Some(Err(lu_solve::Failure::Interrupted(error)))=>Err(error),
        Some(Err(lu_solve::Failure::NonInvertible))=>Err(PyErr::from_value(exception.call0()?)),
    }
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_ddm_ilu_solve,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_ilu_solve_dependencies,module)?)?;
    Ok(())
}
