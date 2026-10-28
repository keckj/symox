//! QQ LU representation checks and exact in-place publication.
use pyo3::prelude::*;
use pyo3::types::{PyList,PyTuple,PyInt};
use symox::poly::{domains::ExactDomain,lu};
use super::{sequence,entries,rationals,element_object};

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

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_ddm_ilu,module)?)?;
    Ok(())
}
