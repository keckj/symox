//! Immutable reference DATA only. No runtime compiler or serializer is called.
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyBytes, PyCFunction, PyInt, PyString, PyTuple};
pub(super) enum CodeValue {
    None, Bool(bool), Int(i64), String(&'static str), Bytes(&'static [u8]),
    Tuple(&'static [CodeValue]), Code(&'static [(&'static str, CodeValue)]),
}
#[path="lll_cold_code_specs.rs"]
mod specs;

// All operations below are on exact immutable builtin values. In particular,
// no Python function.__code__ read, serializer, compiler, user equality or
// string encoding is needed. co_filename/line positions are intentionally not
// semantic fields; a caller must separately decline active tracing/profiling.
fn matches(value: &Bound<'_, PyAny>, expected: &CodeValue, api: &super::getters::Getters) -> PyResult<bool> {
    match expected {
        CodeValue::None => Ok(value.is_none()),
        CodeValue::Bool(expected) => Ok(value.downcast_exact::<PyBool>().is_ok()
            && value.is_truthy()? == *expected),
        CodeValue::Int(expected) => Ok(value.downcast_exact::<PyInt>().is_ok()
            && value.eq(*expected)?),
        CodeValue::String(expected) => Ok(value.downcast_exact::<PyString>().is_ok()
            && value.eq(*expected)?),
        CodeValue::Bytes(expected) => Ok(value.downcast_exact::<PyBytes>()
            .is_ok_and(|value|value.as_bytes() == *expected)),
        CodeValue::Tuple(expected) => {
            let Ok(value) = value.downcast_exact::<PyTuple>() else { return Ok(false); };
            if value.len() != expected.len() { return Ok(false); }
            for (value, expected) in value.iter().zip(*expected) {
                if !matches(&value, expected, api)? { return Ok(false); }
            }
            Ok(true)
        },
        CodeValue::Code(fields) => {
            if !api.is_code(value) { return Ok(false); }
            for (name, expected) in *fields {
                if !matches(&value.getattr(*name)?, expected, api)? { return Ok(false); }
            }
            Ok(true)
        },
    }
}

pub(super) fn function_code_matches(function:&Bound<'_,PyAny>,label:&str,api:&super::getters::Getters)->PyResult<bool> {
    if !api.is_function(function){return Ok(false);}
    let Some((_,expected))=specs::REFERENCES.iter().find(|(name,_)|*name==label) else {return Ok(false);};
    matches(&api.code(function)?,expected,api)
}

#[pyfunction]
fn _lll_code_matches(function:&Bound<'_,PyAny>, label:&Bound<'_,PyAny>)->PyResult<Option<bool>> {
    let Some(api)=super::getters::Getters::load() else { return Ok(None); };
    if !api.is_function(function) || label.downcast_exact::<PyString>().is_err() { return Ok(Some(false)); }
    // A consumed label is compared as Unicode, never encoded to UTF-8.
    for (name, expected) in specs::REFERENCES {
        if label.eq(*name)? { return Ok(Some(matches(&api.code(function)?, expected, &api)?)); }
    }
    Ok(Some(false))
}
fn value(py: Python<'_>, item: &CodeValue) -> PyResult<PyObject> {
    let (tag,payload): (&str,PyObject)=match item {
        CodeValue::None=>("none",py.None()),
        CodeValue::Bool(v)=>("bool",v.into_pyobject(py)?.to_owned().into_any().unbind()),
        CodeValue::Int(v)=>("int",v.into_pyobject(py)?.into_any().unbind()),
        CodeValue::String(v)=>("str",v.into_pyobject(py)?.into_any().unbind()),
        CodeValue::Bytes(v)=>("bytes",PyBytes::new(py,v).into_any().unbind()),
        CodeValue::Tuple(values)=>{
            let values=values.iter().map(|v|value(py,v)).collect::<PyResult<Vec<_>>>()?;
            ("tuple",PyTuple::new(py,values)?.into_any().unbind())
        },
        CodeValue::Code(fields)=>{
            let fields=fields.iter().map(|(name,v)|Ok((*name,value(py,v)?).into_pyobject(py)?.into_any().unbind())).collect::<PyResult<Vec<PyObject>>>()?;
            ("code",PyTuple::new(py,fields)?.into_any().unbind())
        },
    };
    Ok((tag,payload).into_pyobject(py)?.into_any().unbind())
}
pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    module.add_function(wrap_pyfunction!(_lll_code_matches,module)?)?;
    let py=module.py();
    let refs=specs::REFERENCES.iter().map(|(name,v)|Ok((*name,value(py,v)?).into_pyobject(py)?.into_any().unbind())).collect::<PyResult<Vec<PyObject>>>()?;
    let reference=PyTuple::new(py,refs)?.unbind();
    let getter=PyCFunction::new_closure(py,Some(c"_lll_canonical_code_data"),None,
        move |args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,pyo3::types::PyDict>>|->PyResult<PyObject>{
            if !args.is_empty() || kwargs.is_some(){return Err(pyo3::exceptions::PyTypeError::new_err("expected no arguments"));}
            Ok(reference.clone_ref(args.py()).into_any())
        })?;
    module.add("_lll_canonical_code_data",getter)
}
