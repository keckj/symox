//! Optional, documented CPython function getters; never a private object layout.
//! This experimental path is limited to the two verified CPython 3.10/Linux
//! x86_64 environments. Missing APIs leave the complete historical path intact.
use pyo3::prelude::*;
use std::ffi::{c_char, c_int, c_void};
type DictGetter = unsafe extern "C" fn(*mut pyo3::ffi::PyObject, *mut c_void) -> *mut pyo3::ffi::PyObject;

#[derive(Clone, Copy)]
pub(super) struct Getters {
    functions: &'static crate::poly::ddm::function_metadata::Getters,
    code_type: usize,
    instance_dict: DictGetter,
}

impl Getters {
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    pub(super) fn load() -> Option<Self> {
        // Reuse the accepted shared function metadata API and its verified
        // CPython3.10.12/.20 whitelist. Only these two additional documented
        // symbols are optional dependencies of the cold graph verifier.
        let functions=crate::poly::ddm::function_metadata::api()?;
        #[link(name = "dl")]
        unsafe extern "C" {
            fn dlopen(name: *const c_char, flags: c_int) -> *mut c_void;
            fn dlsym(handle: *mut c_void, name: *const c_char) -> *mut c_void;
            fn dlclose(handle: *mut c_void) -> c_int;
        }
        unsafe {
            let handle = dlopen(std::ptr::null(), 1); // RTLD_LAZY, process image only.
            if handle.is_null() { return None; }
            let symbols=[c"PyObject_GenericGetDict",c"PyCode_Type"]
                .map(|name|dlsym(handle,name.as_ptr()));
            let _=dlclose(handle);
            if symbols.iter().any(|symbol|symbol.is_null()){return None;}
            Some(Self{
                functions,
                code_type:symbols[1] as usize,
                instance_dict:std::mem::transmute::<*mut c_void,DictGetter>(symbols[0]),
            })
        }
    }

    #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
    pub(super) fn load() -> Option<Self> { None }

    pub(super) fn is_function(&self, value: &Bound<'_,PyAny>) -> bool {
        self.functions.is_function(value)
    }

    pub(super) fn is_code(&self, value: &Bound<'_,PyAny>) -> bool {
        value.get_type().as_ptr() as usize == self.code_type
    }

    pub(super) fn dictionary<'py>(&self, value: &Bound<'py,PyAny>) -> PyResult<Option<Bound<'py,PyAny>>> {
        // Public API (new reference), available on the verified 3.10 runtime.
        let object=unsafe { (self.instance_dict)(value.as_ptr(),std::ptr::null_mut()) };
        if !object.is_null() { return Ok(Some(unsafe { Bound::from_owned_ptr(value.py(),object) })); }
        let error=PyErr::fetch(value.py());
        if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(value.py())
            || error.is_instance_of::<pyo3::exceptions::PyTypeError>(value.py()) { return Ok(None); }
        Err(error)
    }

    pub(super) fn code<'py>(&self,f:&Bound<'py,PyAny>)->PyResult<Bound<'py,PyAny>>{self.functions.field(f,0)}
    pub(super) fn globals<'py>(&self,f:&Bound<'py,PyAny>)->PyResult<Bound<'py,PyAny>>{self.functions.field(f,1)}
    pub(super) fn defaults<'py>(&self,f:&Bound<'py,PyAny>)->PyResult<Bound<'py,PyAny>>{self.functions.field(f,2)}
    pub(super) fn kwdefaults<'py>(&self,f:&Bound<'py,PyAny>)->PyResult<Bound<'py,PyAny>>{self.functions.field(f,3)}
    pub(super) fn closure<'py>(&self,f:&Bound<'py,PyAny>)->PyResult<Bound<'py,PyAny>>{self.functions.field(f,4)}
}
