//! Smith normal decomposition over stock ZZ, and the integer Hermite kernels
//! of the number-field module helpers. Admission lives in Python.
use pyo3::prelude::*;
use pyo3::types::{PyDict,PyInt,PyList,PyModule,PyString,PyTuple};
use symox::poly::{smith,hnf,hnf_modular,domains::ExactDomain,Int,Rat};
use super::{entries,integers,sequence,element_object,integer_object,Scalar};

fn shape(value:&Bound<'_,PyAny>)->PyResult<Option<(usize,usize)>> {
    let Ok(value)=value.downcast_exact::<PyTuple>() else{return Ok(None)};
    if value.len()!=2{return Ok(None);}
    let (m,n)=(value.get_item(0)?,value.get_item(1)?);
    if !m.is_exact_instance_of::<PyInt>()||!n.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let (Ok(m),Ok(n))=(m.extract::<usize>(),n.extract::<usize>()) else{return Ok(None)};
    Ok(Some((m,n)))
}

/// Exact list of exact lists, distinct rows, with `kind` scalars and the given shape.
fn matrix(a:&Bound<'_,PyAny>,m:usize,n:usize,kind:Scalar)->PyResult<Option<Vec<Vec<Int>>>> {
    let py=a.py();
    let Some(rows)=sequence(a,None,true) else{return Ok(None)};
    let mut seen=rustc_hash::FxHashSet::default();
    if rows.iter().any(|r|!r.is_exact_instance_of::<PyList>()||!seen.insert(r.as_ptr() as usize)){return Ok(None);}
    let Some(cells)=entries(&rows,true) else{return Ok(None)};
    if cells.len()!=m||cells.iter().any(|r|r.len()!=n){return Ok(None);}
    integers(&cells,kind,py)
}

fn zz(py:Python<'_>,value:Int)->PyResult<PyObject> {element_object(py,Rat::from(value),ExactDomain::ZZ)}

fn rows_object(py:Python<'_>,rows:Vec<Vec<Int>>,element:fn(Python<'_>,Int)->PyResult<PyObject>)->PyResult<PyObject> {
    let out=PyList::empty(py);
    for row in rows {
        py.check_signals()?;
        let line=PyList::empty(py);
        for value in row {line.append(element(py,value)?)?;}
        out.append(line)?;
    }
    Ok(out.into_any().unbind())
}

/// `invariants` (full false) or `(invariants, s, t)`; None declines.
#[pyfunction]
fn _ddm_smith(a:&Bound<'_,PyAny>,dimensions:&Bound<'_,PyAny>,full:bool,py:Python<'_>)->PyResult<Option<PyObject>> {
    py.check_signals()?;
    let Some((m,n))=shape(dimensions)? else{return Ok(None)};
    let Some(values)=matrix(a,m,n,Scalar::ZZ)? else{return Ok(None)};
    let Some((invs,transforms))=smith::decompose(values,m,n,full,&mut ||py.check_signals())? else{return Ok(None)};
    let invs=invs.into_iter().map(|v|zz(py,v)).collect::<PyResult<Vec<_>>>()?;
    let invs=PyTuple::new(py,invs)?.into_any().unbind();
    match transforms {
        Some((s,t))=>Ok(Some(PyTuple::new(py,[invs,rows_object(py,s,zz)?,rows_object(py,t,zz)?])?.into_any().unbind())),
        None=>Ok(Some(invs)),
    }
}

fn int(py:Python<'_>,value:Int)->PyResult<PyObject> {integer_object(py,&value)}

/// Number-field ordinary HNF on exact Python ints: the cropped rows, or None.
#[pyfunction]
fn _nf_hnf(a:&Bound<'_,PyAny>,dimensions:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<Option<PyObject>> {
    py.check_signals()?;
    let Some((m,n))=shape(dimensions)? else{return Ok(None)};
    let Some(mut values)=matrix(a,m,n,Scalar::Integer)? else{return Ok(None)};
    let Some(k)=hnf::ordinary(&mut values,n,&mut ||py.check_signals(),&mut |_|Ok(()))? else{return Ok(None)};
    let rows=values.into_iter().map(|row|row[k..].to_vec()).collect();
    Ok(Some((rows_object(py,rows,int)?,n-k).into_pyobject(py)?.into_any().unbind()))
}

/// Number-field modulo-D HNF on exact Python ints: the m x m rows, or None.
#[pyfunction]
fn _nf_hnf_modular(a:&Bound<'_,PyAny>,dimensions:&Bound<'_,PyAny>,modulus:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<Option<PyObject>> {
    py.check_signals()?;
    let Some((m,n))=shape(dimensions)? else{return Ok(None)};
    if !modulus.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let Some(mut values)=matrix(a,m,n,Scalar::Integer)? else{return Ok(None)};
    let Some(d)=integers(&[vec![modulus.clone()]],Scalar::Integer,py)? else{return Ok(None)};
    let d=d[0][0].clone();
    let Some(w)=hnf_modular::modular(&mut values,n,&d,&mut ||py.check_signals(),&mut |_|Ok(()))? else{return Ok(None)};
    Ok(Some(rows_object(py,w,int)?))
}

/// Identity proof of a captured protocol: every `(dict, name, expected)` entry
/// still binds `expected` (`missing` for an absent name), and every
/// `(function, metadata)` still has the same five public metadata fields.
/// Reads only: exact-dict lookups of exact-str names and the public getters.
#[pyfunction]
fn _snapshot_unchanged(entries:&Bound<'_,PyAny>,functions:&Bound<'_,PyAny>,missing:&Bound<'_,PyAny>)->PyResult<bool> {
    let Some(api)=super::function_metadata::api() else{return Ok(false)};
    let (Ok(entries),Ok(functions))=(entries.downcast_exact::<PyTuple>(),functions.downcast_exact::<PyTuple>()) else{return Ok(false)};
    for entry in entries.iter() {
        let Ok(entry)=entry.downcast_exact::<PyTuple>() else{return Ok(false)};
        if entry.len()!=3{return Ok(false);}
        let container=entry.get_item(0)?;let name=entry.get_item(1)?;let expected=entry.get_item(2)?;
        let Ok(dict)=container.downcast_exact::<PyDict>() else{return Ok(false)};
        if !name.is_exact_instance_of::<PyString>(){return Ok(false);}
        let same=match dict.get_item(&name)? {Some(value)=>value.is(&expected),None=>expected.is(missing)};
        if !same{return Ok(false);}
    }
    for entry in functions.iter() {
        let Ok(entry)=entry.downcast_exact::<PyTuple>() else{return Ok(false)};
        if entry.len()!=2{return Ok(false);}
        let function=entry.get_item(0)?;let metadata=entry.get_item(1)?;
        let Ok(metadata)=metadata.downcast_exact::<PyTuple>() else{return Ok(false)};
        if metadata.len()!=5||!api.is_function(&function){return Ok(false);}
        for index in 0..5 {if !api.field(&function,index)?.is(&metadata.get_item(index)?){return Ok(false);}}
    }
    Ok(true)
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_snapshot_unchanged,module)?)?;
    module.add_function(wrap_pyfunction!(_ddm_smith,module)?)?;
    module.add_function(wrap_pyfunction!(_nf_hnf,module)?)?;
    module.add_function(wrap_pyfunction!(_nf_hnf_modular,module)?)?;
    Ok(())
}
