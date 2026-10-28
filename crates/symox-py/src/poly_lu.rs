//! Exact field LU (QQ, GF(p), QQ<a>) representation checks and in-place publication.
use pyo3::prelude::*;
use pyo3::types::{PyList,PyTuple,PyInt};
use symox::poly::{domains::ExactDomain,lu};
use super::{sequence,entries,rationals,element_object,scalar,Scalar,field_domain,field_values,field_object,with_field_scope,FieldKind,FieldPayload};

#[pyfunction]
#[pyo3(signature=(a,stock=None,dimensions=None))]
fn _ddm_ilu<'py>(a:&Bound<'py,PyAny>,stock:Option<&Bound<'py,PyAny>>,dimensions:Option<&Bound<'py,PyAny>>,py:Python<'py>)
    ->PyResult<Option<Vec<(usize,usize)>>> {
    py.check_signals()?;
    let Some(rows)=sequence(a,stock,true) else{return Ok(None)};
    let mut seen=rustc_hash::FxHashSet::default();
    for row in &rows {
        if !row.is_exact_instance_of::<PyList>()||!seen.insert(row.as_ptr() as usize){return Ok(None);}
    }
    let Some(cells)=entries(&rows,true) else{return Ok(None)};
    if let Some(first)=cells.first(){if cells.iter().any(|row|row.len()!=first.len()){return Ok(None);}}
    // The real Python length prefix may use user callbacks. Its actual
    // dimensions must agree with these canonical rows before any mutation;
    // refusal resumes below that prefix without replaying those callbacks.
    if let Some(dimensions)=dimensions.filter(|value|!value.is_none()) {
        let Ok(dimensions)=dimensions.downcast_exact::<PyTuple>() else{return Ok(None)};
        if dimensions.len()!=2{return Ok(None);}
        let m=dimensions.get_item(0)?;let n=dimensions.get_item(1)?;
        if !m.is_exact_instance_of::<PyInt>()||!n.is_exact_instance_of::<PyInt>(){return Ok(None);}
        let (Ok(m),Ok(n))=(m.extract::<usize>(),n.extract::<usize>()) else{return Ok(None)};
        if m!=cells.len()||n!=cells.first().map_or(0,Vec::len){return Ok(None);}
    }
    if let Some(first)=cells.first().and_then(|row|row.first()).filter(|v|scalar(v)!=Some(Scalar::QQ)) {
        // QQ_I division routes through the Python Gaussian wrapper, whose
        // conversion callbacks are not authenticated here: Python keeps it.
        let Some((domain,kind))=field_domain(first,false) else{return Ok(None)};
        let outer=a.downcast::<PyList>()?;
        let output=rows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
        return match kind {
            FieldKind::Finite=>field_lu::<symox::poly::f5b::Fp>(&cells,outer,output,&domain,py),
            FieldKind::Algebraic=>field_lu::<symox::poly::algebraic::AN>(&cells,outer,output,&domain,py),
        };
    }
    let Some(mut values)=rationals(&cells,py)? else{return Ok(None)};
    let outer=a.downcast::<PyList>()?;
    let mut output=rows.into_iter().map(|row|row.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
    lu::in_place(&mut values,&mut ||py.check_signals(),&mut |event| {
        match event {
            lu::Event::Swap(i,j)=>{
                outer.set_item(i,&output[j])?;
                outer.set_item(j,&output[i])?;
                output.swap(i,j);
                Ok(())
            },
            lu::Event::Set(i,j,value)=>output[i].set_item(j,element_object(py,value.clone(),ExactDomain::QQ)?),
        }
    })
}

fn field_lu<C:FieldPayload>(cells:&[Vec<Bound<'_,PyAny>>],outer:&Bound<'_,PyList>,mut output:Vec<Bound<'_,PyList>>,
    domain:&symox::poly::domains::ExactDomain,py:Python<'_>)->PyResult<Option<Vec<(usize,usize)>>> {
    let Some(mut values)=field_values::<C>(cells,domain,py)? else{return Ok(None)};
    with_field_scope(domain,||lu::in_place(&mut values,&mut ||py.check_signals(),&mut |event| {
        match event {
            lu::Event::Swap(i,j)=>{
                outer.set_item(i,&output[j])?;outer.set_item(j,&output[i])?;output.swap(i,j);Ok(())
            },
            lu::Event::Set(i,j,value)=>output[i].set_item(j,field_object(py,value.clone(),domain)?),
        }
    })).map(Option::flatten)
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_ddm_ilu,module)?)?;
    Ok(())
}
