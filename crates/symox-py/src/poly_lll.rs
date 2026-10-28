//! Exact LLL publication after the public constructor/conversion prefix.
#[path="poly_lll_cache.rs"]
mod cache;
#[path="poly_lll_first_call.rs"]
mod first_call;
use pyo3::prelude::*;
use pyo3::types::{PyInt,PyList};
use symox::poly::{Int,Rat,domains::ExactDomain,lll};
use super::{scalar,integers,element_object,Scalar};
type PublishedRows=(Vec<PyObject>,Option<Vec<PyObject>>);

fn rows<'py>(value:&Bound<'py,PyAny>,stock:&Bound<'py,PyAny>)->Option<Vec<Vec<Bound<'py,PyAny>>>> {
    if !value.get_type().is(stock){return None;}
    value.downcast::<PyList>().ok()?.iter().map(|row|{
        if !row.is_exact_instance_of::<PyList>(){return None;}
        Some(row.downcast::<PyList>().ok()?.iter().collect())
    }).collect()
}

fn publish(py:Python<'_>,values:Vec<Vec<Int>>,sources:&[Option<usize>],
    original:&[Bound<'_,PyAny>])->PyResult<Vec<PyObject>> {
    values.into_iter().zip(sources).map(|(row,source)|{
        py.check_signals()?;
        if let Some(source)=source{return Ok(original[*source].clone().unbind());}
        let out=PyList::empty(py);
        for (i,value) in row.into_iter().enumerate(){
            if i%64==0{py.check_signals()?;}
            out.append(element_object(py,Rat::from(value),ExactDomain::ZZ)?)?;
        }
        Ok(out.into_any().unbind())
    }).collect()
}

/// None declines noncanonical prefixes. Input matrices remain unchanged until
/// both result row vectors have been constructed successfully by the caller.
#[pyfunction]
#[allow(clippy::too_many_arguments)] // Prefix objects and protocol classes remain explicit, never global.
fn _ddm_lll_native(y:&Bound<'_,PyAny>,columns:&Bound<'_,PyAny>,delta:&Bound<'_,PyAny>,
    transform:&Bound<'_,PyAny>,first:&Bound<'_,PyAny>,stock:&Bound<'_,PyAny>,
    rank_error:&Bound<'_,PyAny>,py:Python<'_>)->PyResult<Option<PublishedRows>> {
    py.check_signals()?;
    if !columns.is_exact_instance_of::<PyInt>()||scalar(delta)!=Some(Scalar::QQ){return Ok(None);}
    let Ok(columns)=columns.extract::<usize>() else{return Ok(None)};
    let Some(input)=rows(y,stock) else{return Ok(None)};
    if input.is_empty()||input.iter().any(|r|r.len()!=columns){return Ok(None);}
    let Some(values)=integers(&input,Scalar::ZZ,py)? else{return Ok(None)};
    // A temporary conversion hook may restore itself before selection. Never
    // replay or discard its actual result when it differs from this input.
    if !first.is_exact_instance_of::<PyList>(){return Ok(None);}
    let converted=first.downcast::<PyList>()?;
    if converted.len()!=columns{return Ok(None);}
    for (i,value) in converted.iter().enumerate(){
        if i%64==0{py.check_signals()?;}
        if scalar(&value)!=Some(Scalar::QQ){return Ok(None);}
        let value=value.downcast::<crate::poly::domains::PyDomainElement>()?.get();
        let symox::poly::domains::ElementValue::Rational(value)=&value.value.value else{unreachable!()};
        if !value.d.is_one()||value.n!=values[0][i]{return Ok(None);}
    }
    let t_original=if transform.is_none(){None}else{
        let Some(t)=rows(transform,stock) else{return Ok(None)};
        let Some(t_values)=integers(&t,Scalar::ZZ,py)? else{return Ok(None)};
        if t_values.len()!=values.len()||t_values.iter().any(|row|row.len()!=values.len()){return Ok(None);}
        for (i,row) in t_values.iter().enumerate(){for (j,value) in row.iter().enumerate(){
            if value!=&Int::from(i64::from(i==j)){return Ok(None);}
        }}
        Some(transform.downcast::<PyList>()?.iter().collect::<Vec<_>>())
    };
    let delta=delta.downcast::<crate::poly::domains::PyDomainElement>()?.get();
    let symox::poly::domains::ElementValue::Rational(delta)=&delta.value.value else{unreachable!()};
    let result=match lll::checked(&values,columns,delta,t_original.is_some(),&mut ||py.check_signals()) {
        Ok(result)=>result,
        Err(lll::Error::Interrupted(error))=>return Err(error),
        Err(lll::Error::Rank)=>return Err(PyErr::from_value(rank_error.call1(("input matrix contains linearly dependent rows",))?)),
        Err(lll::Error::Delta|lll::Error::Shape)=>return Ok(None),
    };
    let original=y.downcast::<PyList>()?.iter().collect::<Vec<_>>();
    let basis=publish(py,result.basis,&result.row_sources,&original)?;
    let transform=match (result.transform,t_original){
        (Some(t),Some(original))=>Some(publish(py,t,&result.row_sources,&original)?),
        _=>None,
    };
    Ok(Some((basis,transform)))
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    cache::register(module)?;
    first_call::register(module)?;
    module.add_function(wrap_pyfunction!(_ddm_lll_native,module)?)?;Ok(())
}
