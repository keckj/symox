//! Optional public CPython function getters, without audit callbacks/layouts.
//! Only CPython 3.10.12 and 3.10.20/Linux x86_64 are verified for this route.
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use std::sync::OnceLock;
type Getter=unsafe extern "C" fn(*mut pyo3::ffi::PyObject)->*mut pyo3::ffi::PyObject;

#[derive(Clone,Copy)]
pub(super) struct Getters { functions:[Getter;5],function_type:usize }
impl Getters {
    #[cfg(all(target_os="linux",target_arch="x86_64"))]
    fn load()->Option<Self> {
        use std::ffi::{c_char,c_int,c_void,CStr};
        // Public C API: docs.python.org/3.10/c-api/function.html. Resolve it
        // optionally: these symbols are not hard dependencies of the abi3 DSO.
        let version=unsafe {CStr::from_ptr(pyo3::ffi::Py_GetVersion())}.to_bytes();
        if !(version.starts_with(b"3.10.12 ")||version.starts_with(b"3.10.20 "))
            ||version.windows(4).any(|v|v==b"PyPy"){return None;}
        #[link(name="dl")]
        unsafe extern "C" {
            fn dlopen(name:*const c_char,flags:c_int)->*mut c_void;
            fn dlsym(handle:*mut c_void,name:*const c_char)->*mut c_void;
            fn dlclose(handle:*mut c_void)->c_int;
        }
        unsafe {
            let handle=dlopen(std::ptr::null(),1); // RTLD_LAZY, process image.
            if handle.is_null(){return None;}
            let functions=[c"PyFunction_GetCode",c"PyFunction_GetGlobals",c"PyFunction_GetDefaults",
                c"PyFunction_GetKwDefaults",c"PyFunction_GetClosure"].map(|name|dlsym(handle,name.as_ptr()));
            let function_type=dlsym(handle,c"PyFunction_Type".as_ptr());
            let _=dlclose(handle); // Interpreter image stays loaded for its life.
            if function_type.is_null()||functions.iter().any(|p|p.is_null()){return None;}
            Some(Self{functions:functions.map(|p|std::mem::transmute::<*mut c_void,Getter>(p)),function_type:function_type as usize})
        }
    }
    #[cfg(not(all(target_os="linux",target_arch="x86_64")))]
    fn load()->Option<Self>{None}

    pub(super) fn is_function(&self,function:&Bound<'_,PyAny>)->bool {
        function.get_type().as_ptr() as usize==self.function_type
    }

    /// Shared native consumers use the same public API identities and reread
    /// each live field. No function metadata or admission result is cached.
    pub(super) fn field<'py>(&self,function:&Bound<'py,PyAny>,index:usize)->PyResult<Bound<'py,PyAny>> {
        if !self.is_function(function)||index>=self.functions.len(){
            return Err(pyo3::exceptions::PyTypeError::new_err("expected a native function metadata field"));
        }
        let py=function.py();let ptr=unsafe {(self.functions[index])(function.as_ptr())};
        if ptr.is_null(){
            if PyErr::occurred(py){Err(PyErr::fetch(py))}else{Ok(py.None().into_bound(py))}
        }else{Ok(unsafe {Bound::from_borrowed_ptr(py,ptr)})}
    }

    fn read(&self,function:&Bound<'_,PyAny>)->PyResult<Option<PyObject>> {
        if function.get_type().as_ptr() as usize!=self.function_type{return Ok(None);}
        let py=function.py();let mut values=Vec::with_capacity(5);
        // Read every field anew; only public borrowed-reference getters are
        // used. No cached function data, struct layout or Python attribute read.
        for getter in self.functions {
            let ptr=unsafe {getter(function.as_ptr())};
            let value=if ptr.is_null(){
                if PyErr::occurred(py){return Err(PyErr::fetch(py));}py.None()
            }else{unsafe {Bound::from_borrowed_ptr(py,ptr)}.unbind()};
            values.push(value);
        }
        Ok(Some(PyTuple::new(py,values)?.into_any().unbind()))
    }
}

// Only the immutable public API addresses/type identity are retained. Each
// function's current metadata is read anew and is never memoized.
static GETTERS:OnceLock<Option<Getters>>=OnceLock::new();

pub(super) fn api()->Option<&'static Getters> {
    GETTERS.get_or_init(Getters::load).as_ref()
}

#[pyfunction]
fn _ddm_function_metadata(function:&Bound<'_,PyAny>)->PyResult<Option<PyObject>> {
    match GETTERS.get_or_init(Getters::load){Some(getters)=>getters.read(function),None=>Ok(None)}
}

pub(super) fn register(module:&Bound<'_,PyModule>)->PyResult<()> {
    let _=GETTERS.get_or_init(Getters::load);
    module.add_function(wrap_pyfunction!(_ddm_function_metadata,module)?)?;
    Ok(())
}
