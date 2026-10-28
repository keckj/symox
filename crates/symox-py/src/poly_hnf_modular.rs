//! Exact modulo-D HNF, preserving working cells, W insertion and modulus objects.
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyInt, PyList, PyModule, PyString, PyTuple};
use pyo3::exceptions::PyTypeError;
use symox::poly::{hnf_modular, domains::ExactDomain, Int, Rat};
use super::{entries, integers, scalar, sequence, element_object, Scalar};

fn dictionary<'py>(v:&Bound<'py,PyAny>)->Option<Bound<'py,PyDict>> {
    let d=v.downcast_exact::<PyDict>().ok()?;
    if d.iter().any(|(k,_)|!k.is_exact_instance_of::<PyString>()){return None;}
    Some(d.clone())
}
fn lookup<'py>(g:&Bound<'py,PyDict>,b:&Bound<'py,PyDict>,name:&str)->PyResult<PyObject> {
    let v=if let Some(v)=g.get_item(name)?{Some(v)}else{b.get_item(name)?};
    Ok(v.map_or_else(||g.py().None(),Bound::unbind))
}
fn function_dicts<'py>(f:&Bound<'py,PyAny>,reader:&Bound<'py,PyAny>)->PyResult<Option<(Bound<'py,PyDict>,Bound<'py,PyDict>)>> {
    let data=reader.call1((f,))?;
    let Ok(data)=data.downcast_exact::<PyTuple>() else{return Ok(None)};
    if data.len()!=5{return Ok(None);}
    let globals=data.get_item(1)?;let Some(g)=dictionary(&globals) else{return Ok(None)};
    // The verified optional service above accepts exact Python functions only.
    // This public readonly member does not audit on the verified runtimes.
    let builtins=f.getattr(pyo3::intern!(f.py(),"__builtins__"))?;
    let Some(b)=dictionary(&builtins) else{return Ok(None)};Ok(Some((g,b)))
}
fn dependencies(py:Python<'_>,print:&Bound<'_,PyAny>,metadata:Option<&PyObject>)->PyResult<Option<PyObject>> {
    let Some(metadata)=metadata else{return Ok(None)};
    let frame=unsafe{pyo3::ffi::PyEval_GetFrame()};if frame.is_null(){return Ok(None);}
    let frame=unsafe{Bound::from_borrowed_ptr(py,frame.cast())};
    let g=frame.getattr(pyo3::intern!(py,"f_globals"))?;let b=frame.getattr(pyo3::intern!(py,"f_builtins"))?;
    let (Some(g),Some(b))=(dictionary(&g),dictionary(&b)) else{return Ok(None)};
    let mut out=Vec::new();
    for name in ["range","len","ZZ","_gcdex","add_columns","symmetric_residue"]{out.push(lookup(&g,&b,name)?);}
    let reader=metadata.bind(py);
    let Some((gg,gb))=function_dicts(out[3].bind(py),reader)? else{return Ok(None)};
    let Some((ag,ab))=function_dicts(out[4].bind(py),reader)? else{return Ok(None)};
    let Some((eg,eb))=function_dicts(print,reader)? else{return Ok(None)};
    out.push(lookup(&gg,&gb,"ZZ")?);out.push(lookup(&ag,&ab,"range")?);out.push(lookup(&ag,&ab,"len")?);
    out.push(lookup(&eg,&eb,"type")?);out.push(lookup(&eg,&eb,"int")?);
    Ok(Some(PyTuple::new(py,out)?.into_any().unbind()))
}

#[pyfunction]
fn _ddm_hnf_modular(a:&Bound<'_,PyAny>,w:&Bound<'_,PyAny>,modulus:&Bound<'_,PyAny>,shape:&Bound<'_,PyAny>,defaultdict:&Bound<'_,PyAny>)->PyResult<bool> {
    let py=a.py();py.check_signals()?;
    let Ok(shape)=shape.downcast_exact::<PyTuple>() else{return Ok(false)};
    if shape.len()!=2{return Ok(false);}
    let (m,n)=(shape.get_item(0)?,shape.get_item(1)?);
    if !m.is_exact_instance_of::<PyInt>()||!n.is_exact_instance_of::<PyInt>(){return Ok(false);}
    let (Ok(m),Ok(n))=(m.extract::<usize>(),n.extract::<usize>()) else{return Ok(false)};
    if n<m||!a.is_exact_instance_of::<PyList>()||!w.get_type().is(defaultdict){return Ok(false);}
    let Ok(output)=w.downcast::<PyDict>() else{return Ok(false)};
    if !output.is_empty()||!w.getattr(pyo3::intern!(py,"default_factory"))?.is(py.get_type::<PyDict>()){return Ok(false);}
    let Some(rows)=sequence(a,None,true) else{return Ok(false)};
    let mut seen=rustc_hash::FxHashSet::default();
    if rows.iter().any(|r|!r.is_exact_instance_of::<PyList>()||!seen.insert(r.as_ptr() as usize)){return Ok(false);}
    let Some(cells)=entries(&rows,true) else{return Ok(false)};
    if cells.len()!=m||cells.iter().any(|r|r.len()!=n)||scalar(modulus)!=Some(Scalar::ZZ){return Ok(false);}
    let Some(mut values)=integers(&cells,Scalar::ZZ,py)? else{return Ok(false)};
    let Some(mut d)=integers(&[vec![modulus.clone()]],Scalar::ZZ,py)? else{return Ok(false)};
    let d=d.remove(0).remove(0);if d<=Int::zero(){return Ok(false);}
    let rows=rows.iter().map(|r|r.downcast::<PyList>().unwrap().clone()).collect::<Vec<_>>();
    let mut current_modulus=modulus.clone().unbind();
    // Original range(m-1,-1,-1) yields one key per column, shared by every
    // row dictionary. Original range(m) supplies separate outer keys. Values
    // above CPython's small-int cache must remain distinct across these sets.
    let row_keys=(0..m).map(|i|i.into_pyobject(py).unwrap().into_any().unbind()).collect::<Vec<_>>();
    let col_keys=(0..m).map(|i|i.into_pyobject(py).unwrap().into_any().unbind()).collect::<Vec<_>>();
    let mut output_rows:Vec<Option<Bound<'_,PyDict>>>=vec![None;m];
    let mut write=|i:usize,j:usize,value:&Bound<'_,PyAny>|->PyResult<()> {
        let row=if let Some(row)=&output_rows[i]{row.clone()}else{
            let row=PyDict::new(py);output.set_item(row_keys[i].bind(py),&row)?;output_rows[i]=Some(row.clone());row
        };
        row.set_item(col_keys[j].bind(py),value)
    };
    let result=hnf_modular::modular(&mut values,n,&d,&mut ||py.check_signals(),&mut |event|{
        match event {
            hnf_modular::Event::Work(i,j,v)=>rows[i].set_item(j,element_object(py,Rat::from(v.clone()),ExactDomain::ZZ)?),
            hnf_modular::Event::WorkModulus(i,j)=>rows[i].set_item(j,current_modulus.bind(py)),
            hnf_modular::Event::Output(i,j,v)=>write(i,j,element_object(py,Rat::from(v.clone()),ExactDomain::ZZ)?.bind(py)),
            hnf_modular::Event::OutputModulus(i,j)=>write(i,j,current_modulus.bind(py)),
            hnf_modular::Event::ModulusUpdate(v)=>{current_modulus=element_object(py,Rat::from(v.clone()),ExactDomain::ZZ)?;Ok(())},
        }
    })?;Ok(result.is_some())
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_ddm_hnf_modular,module)?)?;
    let metadata=if let Some(d)=dictionary(module.dict().as_any()){d.get_item("_ddm_function_metadata")?.map(Bound::unbind)}else{None};
    module.add("_ddm_hnf_modular_dependencies",pyo3::types::PyCFunction::new_closure(module.py(),Some(c"_ddm_hnf_modular_dependencies"),None,
        move |args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>|->PyResult<Option<PyObject>> {
            if kwargs.is_some()||args.len()!=1{return Err(PyTypeError::new_err("expected comparison helper"));}
            dependencies(args.py(),&args.get_item(0)?,metadata.as_ref())
        })?)?;Ok(())
}
