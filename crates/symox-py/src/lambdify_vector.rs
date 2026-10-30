//! Float64 NumPy array calls of a compiled function, on the array buffers.
//!
//! A call whose arguments are C-contiguous native float64 `numpy.ndarray`s
//! (of one common shape for the NumPy kernels) and Python `float`/`int`
//! scalars reads the input buffers in place and writes the result straight
//! into a new array, in Rust, with the GIL released on large arrays:
//!
//! - the elementary kernels of `lambdify(..., 'numpy')` evaluate the printed
//!   NumPy text with NumPy's own operations (`symox::lambdify_vec`), bit for
//!   bit, which the Python ufunc evaluation computed before;
//! - the other kernels run the compiled program (`Compiled::call_many`) as
//!   before, without the copies of the arguments and of the result.
//!
//! Any other call (other dtypes, layouts, subclasses, broadcasting, 0-d
//! arrays, a text the native evaluation declines) returns `None` and keeps
//! the established path, with its errors.
use std::sync::Mutex;

use numpy::{PyArrayDescrMethods, PyArrayDyn, PyArrayMethods, PyUntypedArray, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use pyo3::types::{PyFloat, PyInt, PyList, PyTuple};
use symox::lambdify_vec::{self, NumpyName, NumpyText, Plan};

use crate::PyCompiled;

/// The NumPy text of a `numpy_vectors` kernel and its plans per pattern of
/// scalar arguments.
#[derive(Default)]
pub(crate) struct VectorCache {
    text: GILOnceCell<Option<(NumpyText, Vec<NumpyName>)>>,
    plans: Mutex<Vec<(u64, Option<std::sync::Arc<Plan>>)>>,
    // the NumPy object evaluation lacks a function: complex arrays go to
    // SymPy's generated function
    complex_reference: std::sync::atomic::AtomicBool,
}

/// Arrays at least this long are evaluated without the GIL.
const RELEASE: usize = 4096;

/// Whether NumPy's float64 `sin`, `cos`, `exp`, `log` and `power` loops are
/// the C library functions here (NumPy uses its own SIMD kernels on
/// AVX-512 processors): checked once on probe values, bit for bit.
fn numpy_libm(py: Python<'_>) -> bool {
    static LIBM: GILOnceCell<bool> = GILOnceCell::new();
    *LIBM.get_or_init(py, || probe(py).unwrap_or(false))
}

fn probe(py: Python<'_>) -> PyResult<bool> {
    let np = py.import("numpy")?;
    let core = np.getattr("_core").or_else(|_| np.getattr("core"))?;
    let features = core.getattr("_multiarray_umath")?.getattr("__cpu_features__")?;
    for name in ["AVX512F", "AVX512_SKX"] {
        if features.get_item(name).and_then(|v| v.is_truthy()).unwrap_or(true) {
            return Ok(false);
        }
    }
    let mut xs: Vec<f64> = vec![0.0, -0.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, f64::MIN_POSITIVE,
        5e-324, -5e-324, 1e-300, 1e300, -1e300, 1e22, 710.0, -745.5, 1e-8, 0.5, 1.0, -1.0, 2.0];
    for i in 0..2048 {
        let t = i as f64;
        xs.push(t * 0.0137 - 14.0);
        xs.push((t * 0.731).sin() * 10f64.powf((t % 40.0) - 20.0));
    }
    let ys: Vec<f64> = xs.iter().enumerate().map(|(i, x)| match i % 4 {
        0 => 2.5,
        1 => -0.75 + (i as f64) * 1e-3,
        2 => 3.0,
        _ => x.abs().sqrt(),
    }).collect();
    let a = numpy::PyArray1::from_slice(py, &xs);
    let b = numpy::PyArray1::from_slice(py, &ys);
    let kwargs = pyo3::types::PyDict::new(py);
    kwargs.set_item("all", "ignore")?;
    let errstate = np.call_method("errstate", (), Some(&kwargs))?;
    errstate.call_method0("__enter__")?;
    let check = || -> PyResult<bool> {
        let unary: [(&str, unsafe extern "C" fn(f64) -> f64); 4] =
            [("sin", libm::sin), ("cos", libm::cos), ("exp", libm::exp), ("log", libm::log)];
        for (name, f) in unary {
            let r = np.getattr(name)?.call1((&a,))?;
            let r = r.downcast::<numpy::PyArray1<f64>>()?.readonly();
            if !r.as_slice()?.iter().zip(&xs).all(|(r, x)| same(*r, unsafe { f(*x) })) {
                return Ok(false);
            }
        }
        let r = np.getattr("power")?.call1((&a, &b))?;
        let r = r.downcast::<numpy::PyArray1<f64>>()?.readonly();
        Ok(r.as_slice()?.iter().zip(xs.iter().zip(&ys)).all(|(r, (x, y))| same(*r, unsafe { libm::pow(*x, *y) })))
    };
    let result = check();
    errstate.call_method1("__exit__", (py.None(), py.None(), py.None()))?;
    result
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

mod libm {
    extern "C" {
        pub fn sin(x: f64) -> f64;
        pub fn cos(x: f64) -> f64;
        pub fn exp(x: f64) -> f64;
        pub fn log(x: f64) -> f64;
        pub fn pow(x: f64, y: f64) -> f64;
    }
}

/// The bindings of the `numpy.<name>` attributes of the text: the original
/// ufunc of that name, or a Python float.
fn bindings(py: Python<'_>, names: &[String]) -> PyResult<Vec<NumpyName>> {
    let np = py.import("numpy")?;
    let ufunc = np.getattr("ufunc")?;
    names.iter().map(|name| {
        let Ok(value) = np.getattr(name.as_str()) else { return Ok(NumpyName::Other) };
        if let Ok(f) = value.downcast_exact::<PyFloat>() {
            return Ok(NumpyName::Float(f.value()));
        }
        if value.get_type().is(&ufunc) && value.getattr("__name__")?.extract::<String>().is_ok_and(|n| &n == name) {
            return Ok(NumpyName::Ufunc);
        }
        Ok(NumpyName::Other)
    }).collect()
}

/// `numpy.ndarray`, once NumPy is imported (an array argument implies it).
fn ndarray_type(py: Python<'_>) -> PyResult<Option<&Py<pyo3::types::PyType>>> {
    static NDARRAY: GILOnceCell<Py<pyo3::types::PyType>> = GILOnceCell::new();
    if let Some(t) = NDARRAY.get(py) {
        return Ok(Some(t));
    }
    if !py.import("sys")?.getattr("modules")?.contains("numpy")? {
        return Ok(None);
    }
    Ok(Some(NDARRAY.get_or_init(py, || py.get_type::<PyUntypedArray>().unbind())))
}

enum Arg<'py> {
    Scalar(f64),
    Array(Bound<'py, PyArrayDyn<f64>>),
}

/// The native evaluation of the call, or `None` to keep the established path.
pub(crate) fn call(py: Python<'_>, f: &PyCompiled, args: &Bound<'_, PyTuple>) -> PyResult<Option<PyObject>> {
    if args.len() > 64 || args.is_empty() || args.len() != f.c.n_args() {
        return Ok(None);
    }
    // scalar calls keep their path at no cost (and NumPy is not imported)
    if args.iter().all(|a| a.is_exact_instance_of::<PyFloat>() || a.is_exact_instance_of::<PyInt>()) {
        return Ok(None);
    }
    let Some(ndarray) = ndarray_type(py)? else { return Ok(None) };
    let ndarray = ndarray.bind(py);
    if !args.iter().any(|a| a.get_type().is(ndarray)) {
        return Ok(None);
    }
    let mut values = Vec::with_capacity(args.len());
    let mut first: Option<usize> = None;
    let (mut complex, mut deferred) = (false, false);
    for (i, a) in args.iter().enumerate() {
        if let Ok(v) = a.downcast_exact::<PyFloat>() {
            values.push(Arg::Scalar(v.value()));
        } else if a.is_exact_instance_of::<PyInt>() {
            // Python ints of an integer polynomial keep its exact evaluation
            if f.exact.get(py).is_some() {
                return Ok(None);
            }
            // an int beyond the float range keeps its OverflowError
            let Ok(v) = a.extract::<f64>() else { return Ok(None) };
            values.push(Arg::Scalar(v));
        } else if a.get_type().is(&ndarray) {
            let dtype = unsafe { a.downcast_unchecked::<PyUntypedArray>() }.dtype();
            if dtype.kind() == b'c' {
                // complex arrays keep NumPy's object evaluation (as before)
                complex = true;
                values.push(Arg::Scalar(0.0));
                continue;
            }
            // native float64 only (the dtype check of `downcast`, cheaper)
            if dtype.num() != numpy::npyffi::NPY_TYPES::NPY_DOUBLE as std::os::raw::c_int
                || dtype.is_native_byteorder() != Some(true) {
                return Ok(None);
            }
            let array = unsafe { a.downcast_unchecked::<PyArrayDyn<f64>>() };
            if array.ndim() == 0 || !array.is_c_contiguous() {
                deferred = true;
            }
            first.get_or_insert(i);
            values.push(Arg::Array(array.clone()));
        } else {
            return Ok(None);
        }
    }
    if complex {
        if f.vector.complex_reference.load(std::sync::atomic::Ordering::Relaxed) {
            if let Some(g) = f.reference_function(py)? {
                return Ok(Some(g.call(args, None)?.unbind()));
            }
        }
        return match f.object_evaluator(py)?.bind(py).call(args, None) {
            Ok(result) => Ok(Some(result.unbind())),
            // a function NumPy lacks (`numpy.gamma`): SymPy's generated
            // function (SciPy's `gamma`) evaluates the complex arrays
            Err(error) if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(py) => {
                match f.reference_function(py)? {
                    Some(g) => {
                        let result = g.call(args, None).map_err(|_| error)?;
                        f.vector.complex_reference.store(true, std::sync::atomic::Ordering::Relaxed);
                        Ok(Some(result.unbind()))
                    }
                    None => Err(error),
                }
            }
            Err(error) => Err(error),
        };
    }
    if deferred {
        return Ok(None);
    }
    let Some(first) = first else { return Ok(None) };
    let Arg::Array(lead) = &values[first] else { unreachable!() };
    let shape = lead.shape().to_vec();
    let n = lead.len();
    for v in &values {
        if let Arg::Array(a) = v {
            // the NumPy kernels broadcast shapes; the compiled program
            // requires equal sizes (its error stays on the established path)
            let same = if f.numpy_vectors { a.shape() == shape.as_slice() } else { a.len() == n };
            if !same {
                return Ok(None);
            }
        }
    }
    let scalars: Vec<f64> = values.iter().map(|v| if let Arg::Scalar(s) = v { *s } else { 0.0 }).collect();
    // the input buffers are read in place (they are not written)
    let slices: Vec<&[f64]> = values.iter().map(|v| match v {
        Arg::Array(a) => unsafe { a.as_slice() }.unwrap_or(&[]),
        Arg::Scalar(_) => &[],
    }).collect();
    if f.numpy_vectors {
        let mask = values.iter().enumerate().fold(0u64, |m, (i, v)| if matches!(v, Arg::Scalar(_)) { m | 1 << i } else { m });
        return match plan(py, f, mask)? {
            Some(plan) => {
                let out = unsafe { PyArrayDyn::<f64>::new(py, shape.as_slice(), false) };
                let dst = unsafe { out.as_slice_mut() }.expect("new C-contiguous array");
                if n >= RELEASE {
                    py.allow_threads(|| plan.eval(&scalars, &slices, dst));
                } else {
                    plan.eval(&scalars, &slices, dst);
                }
                Ok(Some(out.into_any().unbind()))
            }
            None => Ok(None),
        };
    }
    // the compiled program: scalars broadcast against the arrays
    let broadcast: Vec<Vec<f64>> = values.iter().map(|v| match v {
        Arg::Scalar(s) => vec![*s; n],
        Arg::Array(_) => Vec::new(),
    }).collect();
    let inputs: Vec<&[f64]> = slices.iter().zip(&broadcast).map(|(s, b)| if b.is_empty() { *s } else { b.as_slice() }).collect();
    let outs: Vec<Bound<'_, PyArrayDyn<f64>>> = (0..f.c.n_outputs())
        .map(|_| unsafe { PyArrayDyn::<f64>::new(py, shape.as_slice(), false) }).collect();
    {
        let mut dsts: Vec<&mut [f64]> = outs.iter().map(|o| unsafe { o.as_slice_mut() }.expect("new C-contiguous array")).collect();
        if n >= RELEASE {
            py.allow_threads(|| f.c.call_many(&inputs, &mut dsts));
        } else {
            f.c.call_many(&inputs, &mut dsts);
        }
    }
    if f.scalar_out {
        Ok(Some(outs.into_iter().next().expect("one output").into_any().unbind()))
    } else {
        Ok(Some(PyList::new(py, outs)?.into_any().unbind()))
    }
}

impl PyCompiled {
    /// The NumPy evaluation of the expressions on Python objects (complex
    /// and object arrays, NumPy subclasses, external scalars).
    pub(crate) fn object_evaluator<'py>(&'py self, py: Python<'py>) -> PyResult<&'py Py<PyAny>> {
        self.object_evaluator.get_or_try_init(py, || {
            let arguments = PyList::new(py, self.arguments.iter().cloned().map(crate::wrap))?;
            let expressions = PyList::new(py, self.exprs.iter().cloned().map(crate::wrap))?;
            let expression = if self.scalar_out { expressions.get_item(0)? } else { expressions.into_any() };
            Ok::<_, PyErr>(py.import("symox.utilities.lambdify")?
                .getattr("_compile_object")?.call1((arguments, expression))?.unbind())
        })
    }
}

/// The plan of a `numpy_vectors` kernel for the scalar arguments of `mask`.
fn plan(py: Python<'_>, f: &PyCompiled, mask: u64) -> PyResult<Option<std::sync::Arc<Plan>>> {
    if let Some((_, p)) = f.vector.plans.lock().unwrap().iter().find(|(m, _)| *m == mask) {
        return Ok(p.clone());
    }
    let text = f.vector.text.get_or_try_init(py, || -> PyResult<_> {
        let arguments = PyList::new(py, f.arguments.iter().cloned().map(crate::wrap))?;
        let source = py.import("symox._numpy_vector")?.getattr("source")?
            .call1((arguments, crate::wrap(f.exprs[0].clone())))?;
        let Ok(source) = source.extract::<String>() else { return Ok(None) };
        let Some(text) = lambdify_vec::parse_lambda(&source) else { return Ok(None) };
        let bindings = bindings(py, text.names())?;
        Ok(Some((text, bindings)))
    })?;
    let plan = match text {
        Some((text, bindings)) => {
            let scalar: Vec<bool> = (0..f.c.n_args()).map(|i| mask >> i & 1 == 1).collect();
            text.plan(&scalar, bindings, numpy_libm(py)).map(std::sync::Arc::new)
        }
        None => None,
    };
    f.vector.plans.lock().unwrap().push((mask, plan.clone()));
    Ok(plan)
}
