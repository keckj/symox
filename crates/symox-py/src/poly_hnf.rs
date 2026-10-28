//! Ordinary ZZ HNF publication and noncalling authentication of helper lookups.
use pyo3::prelude::*;
use pyo3::exceptions::PyTypeError;
use pyo3::types::{PyCFunction,PyDict,PyInt,PyList,PyModule,PyString,PyTuple};
use symox::poly::{hnf,domains::ExactDomain};
use super::{sequence,entries,integers,Scalar,element_object,integer_object};

fn exact_dictionary<'py>(value:&Bound<'py,PyAny>)->Option<Bound<'py,PyDict>> {
    let dict=value.downcast_exact::<PyDict>().ok()?;
    if dict.iter().any(|(key,_)|!key.is_exact_instance_of::<PyString>()){return None;}
    Some(dict.clone())
}

fn lookup<'py>(globals:&Bound<'py,PyDict>,builtins:&Bound<'py,PyDict>,name:&str)->PyResult<PyObject> {
    let value=if let Some(v)=globals.get_item(name)?{Some(v)}else{builtins.get_item(name)?};
    Ok(value.map_or_else(||globals.py().None(),Bound::unbind))
}

fn function_globals<'py>(function:&Bound<'py,PyAny>,kind:&Bound<'py,PyAny>)
    ->PyResult<Option<(Bound<'py,PyDict>,Bound<'py,PyDict>)>> {
    if !function.get_type().is(kind){return Ok(None);}
    let globals=function.getattr(pyo3::intern!(function.py(),"__globals__"))?;
    let Some(globals)=exact_dictionary(&globals) else{return Ok(None)};
    let Some(builtins)=globals.get_item("__builtins__")? else{return Ok(None)};
    let builtins=if let Some(b)=exact_dictionary(&builtins){b}
        else if let Ok(module)=builtins.downcast_exact::<PyModule>(){module.dict()}
        else{return Ok(None)};
    if builtins.iter().any(|(key,_)|!key.is_exact_instance_of::<PyString>()){return Ok(None);}
    Ok(Some((globals,builtins)))
}

fn dependencies(py:Python<'_>,print_function:&Bound<'_,PyAny>,function_type:&Bound<'_,PyAny>)
    ->PyResult<Option<PyObject>> {
    let frame=unsafe {pyo3::ffi::PyEval_GetFrame()};
    if frame.is_null(){return Ok(None);}
    let frame=unsafe {Bound::from_borrowed_ptr(py,frame.cast())};
    let globals=frame.getattr(pyo3::intern!(py,"f_globals"))?;
    let builtins=frame.getattr(pyo3::intern!(py,"f_builtins"))?;
    let (Some(g),Some(b))=(exact_dictionary(&globals),exact_dictionary(&builtins)) else{return Ok(None)};
    let names=["range","len","ZZ","_gcdex","add_columns"];
    let mut out=Vec::new();for name in names{out.push(lookup(&g,&b,name)?);}
    let Some((gg,gb))=function_globals(out[3].bind(py),function_type)? else{return Ok(None)};
    let Some((ag,ab))=function_globals(out[4].bind(py),function_type)? else{return Ok(None)};
    let Some((eg,eb))=function_globals(print_function,function_type)? else{return Ok(None)};
    out.push(lookup(&gg,&gb,"ZZ")?);
    out.push(lookup(&ag,&ab,"range")?);out.push(lookup(&ag,&ab,"len")?);
    out.push(lookup(&eg,&eb,"type")?);out.push(lookup(&eg,&eb,"int")?);
    Ok(Some(PyTuple::new(py,out)?.into_any().unbind()))
}

#[pyfunction]
#[pyo3(signature=(a,dimensions,stock=None,record_operations=false))]
fn _ddm_hnf(a:&Bound<'_,PyAny>,dimensions:&Bound<'_,PyAny>,stock:Option<&Bound<'_,PyAny>>,
    record_operations:bool,py:Python<'_>)->PyResult<Option<PyObject>> {
    py.check_signals()?;
    let Ok(dimensions)=dimensions.downcast_exact::<PyTuple>() else{return Ok(None)};
    if dimensions.len()!=2{return Ok(None);}
    let m=dimensions.get_item(0)?;let n=dimensions.get_item(1)?;
    if !m.is_exact_instance_of::<PyInt>()||!n.is_exact_instance_of::<PyInt>(){return Ok(None);}
    let (Ok(m),Ok(n))=(m.extract::<usize>(),n.extract::<usize>()) else{return Ok(None)};
    let Some(rows)=sequence(a,stock,true) else{return Ok(None)};
    let mut seen=rustc_hash::FxHashSet::default();
    for row in &rows{if !row.is_exact_instance_of::<PyList>()||!seen.insert(row.as_ptr() as usize){return Ok(None);}}
    let Some(cells)=entries(&rows,true) else{return Ok(None)};
    if cells.len()!=m||cells.iter().any(|row|row.len()!=n){return Ok(None);}
    let Some(mut values)=integers(&cells,Scalar::ZZ,py)? else{return Ok(None)};
    let rows=rows.into_iter().map(|r|r.downcast_into::<PyList>().unwrap()).collect::<Vec<_>>();
    let operations=PyList::empty(py);
    let k=hnf::ordinary(&mut values,n,&mut ||py.check_signals(),&mut |event| {
        match event {
            hnf::Event::Set(i,j,value)=>rows[i].set_item(j,element_object(py,symox::poly::Rat::from(value.clone()),ExactDomain::ZZ)?),
            hnf::Event::Operation(i,j,a,b,c,d)=>{
                if record_operations {
                    let record=PyTuple::new(py,[i.into_pyobject(py)?.into_any().unbind(),j.into_pyobject(py)?.into_any().unbind(),
                        integer_object(py,a)?,integer_object(py,b)?,integer_object(py,c)?,integer_object(py,d)?])?;
                    operations.append(record)?;
                }Ok(())
            },
        }
    })?;
    let Some(k)=k else{return Ok(None)};
    if record_operations {Ok(Some((k,operations).into_pyobject(py)?.into_any().unbind()))}
    else{Ok(Some(k.into_pyobject(py)?.into_any().unbind()))}
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_ddm_hnf,module)?)?;
    let py=module.py();let function_type=py.eval(c"lambda: None",None,None)?.get_type().unbind();
    module.add("_ddm_hnf_dependencies",PyCFunction::new_closure(py,Some(c"_ddm_hnf_dependencies"),None,
        move |args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>|->PyResult<Option<PyObject>> {
            if kwargs.is_some()||args.len()!=1{return Err(PyTypeError::new_err("expected comparison helper"));}
            dependencies(args.py(),&args.get_item(0)?,function_type.bind(args.py()))
        })?)?;
    Ok(())
}
