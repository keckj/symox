//! Python subclasses of native atoms, function heads and expressions.
//!
//! A user subclass of `Symbol`, `Dummy`, `Function`, `Expr`, `Atom` or
//! `Basic` is a real subtype of `_core.Expr`. Native nodes carry an opaque
//! handle to its class, so every expression rebuilt by native code returns
//! instances of that class again.
//!
//! Handles keep a weak reference to the Python instance of their node: the
//! instance owns the native node, and a strong back reference would form a
//! cycle the collector cannot see. Instead every Python expression object
//! containing host nodes holds strong, collector-visible references to their
//! live instances (`_symox_hosts`), as SymPy expressions hold their
//! arguments. Native containers that are not Python expression objects
//! (matrices, tuples, applications of plain undefined functions) do not hold
//! them: when an instance has been collected while its node lives on there
//! only, a fresh instance of the same class is created without the collected
//! attributes.
use std::cmp::Ordering as Order;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use pyo3::impl_::pyclass_init::PyObjectInit;
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use pyo3::types::{PyDict, PyString, PyTuple, PyType, PyWeakrefMethods, PyWeakrefReference};
use symox::expr::{Expr, Func, Kind, UndefinedFunction, WeakExpr};
use symox::host::{HostClass, HostFormat};

use crate::{to_expr, PyExprNode};

const STR: u32 = 1;
const LATEX: u32 = 2;
const RENAMES: u32 = 4;
const EVAL: u32 = 8;
const FDIFF: u32 = 16;
const EXPR: u32 = 32;
const ATOM: u32 = 64;
const DERIV: u32 = 128;
const HASHABLE: u32 = 256;

static NEXT_CLASS: AtomicU64 = AtomicU64::new(1);
static NEXT_HANDLE: AtomicU64 = AtomicU64::new(1);

/// The Python `HostTree` type: plain expression objects owning host nodes.
static HOST_TREE: GILOnceCell<Py<PyType>> = GILOnceCell::new();

// Instances returned by Python constructors during a native operation (an
// evaluated `subs` result) are only weakly bound to their nodes until a
// Python object holds them. The most recent ones are kept alive so that their
// instance attributes survive until the operation's result is wrapped.
const RECENT_CAPACITY: usize = 1024;
static RECENT: Mutex<std::collections::VecDeque<Py<PyAny>>> = Mutex::new(std::collections::VecDeque::new());

fn keep_recent(value: &Bound<'_, PyAny>) {
    let evicted = match RECENT.lock() {
        Ok(mut recent) => {
            recent.push_back(value.clone().unbind());
            if recent.len() > RECENT_CAPACITY { recent.pop_front() } else { None }
        }
        Err(_) => None,
    };
    // Dropped outside the lock: finalizers may create expressions.
    drop(evicted);
}

#[derive(Debug)]
struct ClassInfo {
    id: u64,
    name: String,
    class: Py<PyType>,
    flags: u32,
    class_key: Option<(u8, u32)>,
}

/// Registration of one Python class extending a native atom or head.
#[pyclass(name = "HostClass", module = "symox._core", frozen)]
pub(crate) struct PyHostClass {
    info: Arc<ClassInfo>,
}

#[pymethods]
impl PyHostClass {
    #[getter]
    fn flags(&self) -> u32 { self.info.flags }
}

/// A native function head defined by a Python class.
#[pyclass(name = "HostHead", module = "symox._core", frozen)]
pub(crate) struct PyHostHead {
    head: UndefinedFunction,
}

/// The handle of one host atom, shared by every node of that atom.
#[derive(Debug)]
struct AtomHost {
    class: Arc<ClassInfo>,
    serial: u64,
    node: OnceLock<WeakExpr>,
    instance: Mutex<Option<Py<PyWeakrefReference>>>,
}

/// The node of an application and its live Python instance.
#[derive(Debug)]
struct Application {
    node: WeakExpr,
    instance: Py<PyWeakrefReference>,
}

/// The handle of a host head. The head registered by the class has no
/// application; every application built from Python gets its own handle,
/// bound to its node and instance.
#[derive(Debug)]
struct HeadHost {
    class: Arc<ClassInfo>,
    serial: u64,
    application: Mutex<Option<Application>>,
}

/// The registered Python `HostTree` type.
#[pyfunction]
fn _host_register(py: Python<'_>, tree: &Bound<'_, PyType>) -> PyResult<()> {
    if !tree.is_subclass_of::<PyExprNode>()? {
        return Err(pyo3::exceptions::PyTypeError::new_err("HostTree must extend the native Expr type"));
    }
    let _ = HOST_TREE.set(py, tree.clone().unbind());
    Ok(())
}

/// A new instance of `class` owning `e`, created without calling the class
/// constructor (which may evaluate or reset attributes). The instance holds
/// the live instances of the host nodes below it.
fn instantiate<'py>(py: Python<'py>, class: &Bound<'py, PyType>, e: Expr) -> PyResult<Bound<'py, PyAny>> {
    if !class.is_subclass_of::<PyExprNode>()? {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "{} is not a subtype of the native Expr type", class.name()?)));
    }
    let mut hosts = Vec::new();
    if e.has_host() {
        let mut seen = Vec::new();
        for arg in e.args() { collect(py, &arg, &mut hosts, &mut seen)?; }
    }
    // SAFETY: `class` is a subtype of the native Expr type, whose layout
    // the initializer fills; the new reference is owned by the result.
    let value = unsafe {
        let object = PyClassInitializer::from(PyExprNode { e }).into_new_object(py, class.as_type_ptr())?;
        Bound::from_owned_ptr(py, object)
    };
    if !hosts.is_empty() { hold(py, &value, hosts); }
    Ok(value)
}

/// Store `hosts` in the `_symox_hosts` slot (or instance dictionary) of
/// `value`, bypassing class `__setattr__` overrides. Classes without either
/// storage do not hold their host nodes.
fn hold(py: Python<'_>, value: &Bound<'_, PyAny>, hosts: Vec<PyObject>) {
    let stored = PyTuple::new(py, hosts).and_then(|hosts| {
        let name = PyString::intern(py, "_symox_hosts");
        // SAFETY: all pointers are valid, borrowed references.
        let status = unsafe { pyo3::ffi::PyObject_GenericSetAttr(value.as_ptr(), name.as_ptr(), hosts.as_ptr()) };
        if status < 0 { Err(PyErr::fetch(py)) } else { Ok(()) }
    });
    drop(stored);
}

/// The live instances of the outermost host nodes occurring in `e`; host
/// atoms are revived when needed.
fn collect(py: Python<'_>, e: &Expr, hosts: &mut Vec<PyObject>, seen: &mut Vec<usize>) -> PyResult<()> {
    if !e.has_host() { return Ok(()); }
    let mut push = |value: Bound<'_, PyAny>, hosts: &mut Vec<PyObject>| {
        let key = value.as_ptr() as usize;
        if !seen.contains(&key) { seen.push(key); hosts.push(value.unbind()); }
    };
    match e.kind() {
        Kind::Sym(symbol) => {
            if let Some(host) = symbol.host().and_then(|host| host.as_any().downcast_ref::<AtomHost>()) {
                push(host.object(py, e)?, hosts);
            }
            return Ok(());
        }
        Kind::Func(Func::Undefined(head), _) => {
            if let Some(value) = head.host().and_then(|host| host.as_any().downcast_ref::<HeadHost>())
                .and_then(|host| host.live(py, e)) {
                push(value, hosts);
                return Ok(());
            }
        }
        _ => {}
    }
    for arg in e.args() { collect(py, &arg, hosts, seen)?; }
    Ok(())
}

/// The plain Python expression object of `e`, which contains host nodes.
pub(crate) fn tree_object(py: Python<'_>, e: Expr) -> PyResult<PyObject> {
    match HOST_TREE.get(py) {
        Some(tree) => Ok(instantiate(py, tree.bind(py), e)?.unbind()),
        None => Ok(Bound::new(py, PyExprNode { e })?.into_any().unbind()),
    }
}

fn report(py: Python<'_>, error: PyErr, context: Option<&Bound<'_, PyAny>>) {
    // Native rebuilds cannot propagate exceptions. Report them like
    // exceptions raised in destructors and keep the native result.
    error.write_unraisable(py, context);
}

/// Whether the current thread holds the GIL. Native workers comparing host
/// nodes must not wait for it while the calling thread holds it.
fn holds_gil() -> bool {
    extern "C" { fn PyGILState_Check() -> std::os::raw::c_int; }
    // SAFETY: PyGILState_Check (stable ABI) may be called from any thread.
    unsafe { PyGILState_Check() == 1 }
}

impl AtomHost {
    fn live<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyAny>> {
        self.instance.lock().ok()?.as_ref().and_then(|reference| reference.bind(py).upgrade())
    }

    fn object<'py>(&self, py: Python<'py>, e: &Expr) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = self.live(py) { return Ok(value); }
        let value = instantiate(py, self.class.class.bind(py), e.clone())?;
        self.attach(&value)?;
        Ok(value)
    }

    fn attach(&self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let reference = PyWeakrefReference::new(value)?.unbind();
        let previous = self.instance.lock().ok().and_then(|mut slot| slot.replace(reference));
        drop(previous);
        Ok(())
    }

    /// The instance of this atom, live or revived from its node.
    fn instance<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyAny>> {
        if let Some(value) = self.live(py) { return Some(value); }
        let node = self.node.get()?.upgrade()?;
        self.object(py, &node).ok()
    }
}

impl HeadHost {
    fn applied(&self) -> HeadHost {
        HeadHost { class: self.class.clone(), serial: NEXT_HANDLE.fetch_add(1, Ordering::Relaxed), application: Mutex::new(None) }
    }

    /// The live instance bound to node `e`.
    fn live<'py>(&self, py: Python<'py>, e: &Expr) -> Option<Bound<'py, PyAny>> {
        let application = self.application.lock().ok()?;
        let application = application.as_ref().filter(|application| application.node.is(e))?;
        application.instance.bind(py).upgrade()
    }

    /// Bind `value` as the instance of node `e`, unless this handle is
    /// bound to another live node.
    fn bind(&self, e: &Expr, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let reference = PyWeakrefReference::new(value)?.unbind();
        let previous = match self.application.lock() {
            Ok(mut slot) => {
                let free = slot.as_ref().is_none_or(|application| application.node.is(e) || application.node.upgrade().is_none());
                if free { slot.replace(Application { node: e.downgrade(), instance: reference }) } else { None }
            }
            Err(_) => None,
        };
        drop(previous);
        Ok(())
    }

    fn object<'py>(&self, py: Python<'py>, e: &Expr) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = self.live(py, e) { return Ok(value); }
        let value = instantiate(py, self.class.class.bind(py), e.clone())?;
        self.bind(e, &value)?;
        Ok(value)
    }

    /// The instance of the node bound to this handle, live or revived.
    fn instance<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyAny>> {
        let node = self.application.lock().ok()?.as_ref()?.node.upgrade()?;
        self.object(py, &node).ok()
    }
}

/// `_hashable_content()` of the instances of two handles of one class, or
/// `None` when either is unavailable.
fn contents<'py>(py: Python<'py>, a: Option<Bound<'py, PyAny>>, b: Option<Bound<'py, PyAny>>) -> Option<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
    let content = |value: Bound<'py, PyAny>| value.call_method0("_hashable_content").map_err(|error| drop(error)).ok();
    let _ = py;
    Some((content(a?)?, content(b?)?))
}

/// SymPy equality of two instances of a class defining `_hashable_content`.
fn same_content<'py>(py: Python<'py>, a: Option<Bound<'py, PyAny>>, b: Option<Bound<'py, PyAny>>) -> bool {
    contents(py, a, b).is_some_and(|(a, b)| a.eq(&b).unwrap_or(false))
}

/// A total order of instances that are not equal: their contents when
/// comparable, otherwise their creation order.
fn order_content<'py>(py: Python<'py>, a: Option<Bound<'py, PyAny>>, b: Option<Bound<'py, PyAny>>, serials: (u64, u64)) -> Order {
    if let Some((a, b)) = contents(py, a, b) {
        if a.eq(&b).unwrap_or(false) { return Order::Equal; }
        if let Ok(less) = a.lt(&b) { return if less { Order::Less } else { Order::Greater }; }
    }
    serials.0.cmp(&serials.1)
}

fn print_object(py: Python<'_>, object: &Bound<'_, PyAny>, format: HostFormat) -> Option<String> {
    let result = py.import("symox._host").and_then(|module| module.getattr("_print")?.call1((object, format.name())));
    match result {
        Ok(text) if text.is_none() => None,
        Ok(text) => match text.extract::<String>() {
            Ok(text) => Some(text),
            Err(error) => { report(py, error, Some(object)); None }
        },
        Err(error) => { report(py, error, Some(object)); None }
    }
}

fn wants(flags: u32, format: HostFormat) -> bool {
    match format {
        HostFormat::Str => flags & STR != 0,
        HostFormat::Latex => flags & LATEX != 0,
        HostFormat::Repr => false,
    }
}

impl HostClass for AtomHost {
    fn identity(&self) -> u64 { self.class.id }
    fn class_name(&self) -> &str { &self.class.name }
    fn current_name(&self, atom: &Expr) -> Option<String> {
        if self.class.flags & RENAMES == 0 { return None; }
        let Kind::Sym(symbol) = atom.kind() else { return None };
        Python::with_gil(|py| {
            let value = self.live(py)?;
            match value.getattr("name").and_then(|name| name.extract::<String>()) {
                Ok(name) if name != symbol.name() => Some(name),
                Ok(_) => None,
                Err(error) => { report(py, error, Some(&value)); None }
            }
        })
    }
    fn print(&self, expr: &Expr, format: HostFormat) -> Option<String> {
        if !wants(self.class.flags, format) { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => print_object(py, &value, format),
            Err(error) => { report(py, error, None); None }
        })
    }
    fn same(&self, other: &dyn HostClass) -> bool {
        if self.class.flags & HASHABLE == 0 { return true; }
        let Some(other) = other.as_any().downcast_ref::<AtomHost>() else { return false };
        holds_gil() && Python::with_gil(|py| same_content(py, self.instance(py), other.instance(py)))
    }
    fn order(&self, other: &dyn HostClass) -> Order {
        if self.class.flags & HASHABLE == 0 { return Order::Equal; }
        let Some(other) = other.as_any().downcast_ref::<AtomHost>() else { return Order::Equal };
        let serials = (self.serial, other.serial);
        if !holds_gil() { return serials.0.cmp(&serials.1); }
        Python::with_gil(|py| order_content(py, self.instance(py), other.instance(py), serials))
    }
    fn as_any(&self) -> &dyn std::any::Any { self }
}

impl HostClass for HeadHost {
    fn identity(&self) -> u64 { self.class.id }
    fn class_name(&self) -> &str { &self.class.name }
    fn print(&self, expr: &Expr, format: HostFormat) -> Option<String> {
        if !wants(self.class.flags, format) { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => print_object(py, &value, format),
            Err(error) => { report(py, error, None); None }
        })
    }
    fn evaluate(&self, args: &[Expr]) -> Option<Expr> {
        if self.class.flags & EVAL == 0 { return None; }
        Python::with_gil(|py| {
            let class = self.class.class.bind(py);
            let result = args.iter().map(|arg| crate::scalar_object(py, arg.clone())).collect::<PyResult<Vec<_>>>()
                .and_then(|args| PyTuple::new(py, args))
                .and_then(|args| class.call1(args))
                .and_then(|value| { keep_recent(&value); to_expr(&value) });
            match result {
                Ok(value) => Some(value),
                Err(error) => { report(py, error, Some(class.as_any())); None }
            }
        })
    }
    fn fdiff(&self, expr: &Expr, index: usize) -> Option<Expr> {
        if self.class.flags & FDIFF == 0 { return None; }
        Python::with_gil(|py| {
            let value = match self.object(py, expr) {
                Ok(value) => value,
                Err(error) => { report(py, error, None); return None; }
            };
            match value.call_method1("fdiff", (index,)).and_then(|partial| to_expr(&partial)) {
                Ok(partial) => Some(partial),
                Err(error) => {
                    // ArgumentIndexError means the class has no closed form
                    // for this argument, so the derivative stays held.
                    let unknown = py.import("symox.core.function")
                        .and_then(|module| module.getattr("ArgumentIndexError"))
                        .is_ok_and(|kind| error.matches(py, kind).unwrap_or(false));
                    if !unknown { report(py, error, Some(&value)); }
                    None
                }
            }
        })
    }
    fn same(&self, other: &dyn HostClass) -> bool {
        if self.class.flags & HASHABLE == 0 { return true; }
        let Some(other) = other.as_any().downcast_ref::<HeadHost>() else { return false };
        holds_gil() && Python::with_gil(|py| same_content(py, self.instance(py), other.instance(py)))
    }
    fn order(&self, other: &dyn HostClass) -> Order {
        if self.class.flags & HASHABLE == 0 { return Order::Equal; }
        let Some(other) = other.as_any().downcast_ref::<HeadHost>() else { return Order::Equal };
        let serials = (self.serial, other.serial);
        if !holds_gil() { return serials.0.cmp(&serials.1); }
        Python::with_gil(|py| order_content(py, self.instance(py), other.instance(py), serials))
    }
    fn expression(&self) -> bool { self.class.flags & EXPR != 0 }
    fn atomic(&self) -> bool { self.class.flags & ATOM != 0 }
    fn class_key(&self) -> Option<(u8, u32)> { self.class.class_key }
    fn derivative(&self, expr: &Expr, x: &Expr) -> Option<Expr> {
        if self.class.flags & DERIV == 0 { return None; }
        Python::with_gil(|py| {
            let value = match self.object(py, expr) {
                Ok(value) => value,
                Err(error) => { report(py, error, None); return None; }
            };
            let result = crate::scalar_object(py, x.clone())
                .and_then(|x| value.call_method1("_eval_derivative", (x,)))
                .and_then(|partial| if partial.is_none() { Ok(None) } else { to_expr(&partial).map(Some) });
            match result {
                Ok(partial) => partial,
                Err(error) => { report(py, error, Some(&value)); None }
            }
        })
    }
    fn as_any(&self) -> &dyn std::any::Any { self }
}

/// The Python object of a host atom or application, when `e` has a host.
pub(crate) fn object(py: Python<'_>, e: &Expr) -> PyResult<Option<PyObject>> {
    match e.kind() {
        Kind::Sym(symbol) => match symbol.host().and_then(|host| host.as_any().downcast_ref::<AtomHost>()) {
            Some(host) => Ok(Some(host.object(py, e)?.unbind())),
            None => Ok(None),
        },
        Kind::Func(Func::Undefined(head), _) => match head.host().and_then(|host| host.as_any().downcast_ref::<HeadHost>()) {
            Some(host) => Ok(Some(host.object(py, e)?.unbind())),
            None => Ok(None),
        },
        _ => Ok(None),
    }
}

/// The Python class owning a host atom or application.
pub(crate) fn class_of<'py>(py: Python<'py>, e: &Expr) -> Option<Bound<'py, PyType>> {
    let info = match e.kind() {
        Kind::Sym(symbol) => &symbol.host()?.as_any().downcast_ref::<AtomHost>()?.class,
        Kind::Func(Func::Undefined(head), _) => &head.host()?.as_any().downcast_ref::<HeadHost>()?.class,
        _ => return None,
    };
    Some(info.class.bind(py).clone())
}

/// The head registered by host class `cls` itself (not by a base class).
pub(crate) fn class_head(cls: &Bound<'_, PyType>) -> PyResult<Option<UndefinedFunction>> {
    let Some(head) = cls.getattr(pyo3::intern!(cls.py(), "__dict__"))?.get_item("_symox_head").ok() else { return Ok(None) };
    let Ok(head) = head.downcast::<PyHostHead>() else { return Ok(None) };
    let head = head.get().head.clone();
    let owned = head.host().and_then(|host| host.as_any().downcast_ref::<HeadHost>())
        .is_some_and(|host| host.class.class.bind(cls.py()).is(cls));
    Ok(owned.then_some(head))
}

/// A new unevaluated application of a host head, bound to its instance.
pub(crate) fn apply<'py>(py: Python<'py>, head: &UndefinedFunction, args: &Bound<'py, PyTuple>) -> PyResult<Bound<'py, PyAny>> {
    let args = args.iter().map(|arg| crate::scalar_expr(&arg)).collect::<PyResult<Vec<_>>>()?;
    if !head.accepts(args.len()) {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "{} accepts argument counts {:?}, got {}", head.name, head.nargs.as_ref().unwrap(), args.len())));
    }
    let Some(host) = head.host().and_then(|host| host.as_any().downcast_ref::<HeadHost>()) else {
        return Err(pyo3::exceptions::PyTypeError::new_err("not a host head"));
    };
    let class = host.class.class.bind(py).clone();
    let head = head.clone().with_host(Arc::new(host.applied()));
    let e = Expr::from_kind(Kind::Func(Func::Undefined(head), args));
    let value = instantiate(py, &class, e.clone())?;
    // Bind the handle actually present in the node.
    if let Kind::Func(Func::Undefined(head), _) = e.kind() {
        if let Some(host) = head.host().and_then(|host| host.as_any().downcast_ref::<HeadHost>()) { host.bind(&e, &value)?; }
    }
    Ok(value)
}

/// Register `cls` as a native host class with the given hook flags and,
/// for expression classes, the first two fields of its `class_key()`.
#[pyfunction]
#[pyo3(signature = (cls, flags, class_key=None))]
fn _host_class(cls: &Bound<'_, PyType>, flags: u32, class_key: Option<(u8, u32)>) -> PyResult<PyHostClass> {
    let name = cls.getattr("__name__")?.extract::<String>()?;
    let id = NEXT_CLASS.fetch_add(1, Ordering::Relaxed);
    Ok(PyHostClass { info: Arc::new(ClassInfo { id, name, class: cls.clone().unbind(), flags, class_key }) })
}

/// A new instance of the host class owning the native atom `atom`.
#[pyfunction]
fn _host_atom(py: Python<'_>, host: &PyHostClass, atom: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let Kind::Sym(symbol) = to_expr(atom)?.kind().clone() else {
        return Err(pyo3::exceptions::PyTypeError::new_err("host atoms must be symbols"));
    };
    let handle = Arc::new(AtomHost {
        class: host.info.clone(), serial: NEXT_HANDLE.fetch_add(1, Ordering::Relaxed),
        node: OnceLock::new(), instance: Mutex::new(None),
    });
    let e = Expr::from_kind(Kind::Sym(symbol.with_host(handle.clone())));
    let _ = handle.node.set(e.downgrade());
    let value = instantiate(py, host.info.class.bind(py), e)?;
    handle.attach(&value)?;
    Ok(value.unbind())
}

/// The native head defined by a host function or expression class.
#[pyfunction]
#[pyo3(signature = (host, name, assumptions, nargs=None))]
fn _host_head(host: &PyHostClass, name: &str, assumptions: &Bound<'_, PyDict>, nargs: Option<Vec<usize>>) -> PyResult<PyHostHead> {
    let handle = HeadHost { class: host.info.clone(), serial: 0, application: Mutex::new(None) };
    let head = crate::function_spec(name, assumptions, nargs)?.with_host(Arc::new(handle));
    Ok(PyHostHead { head })
}

/// The unevaluated application of a host head; the class evaluates first.
#[pyfunction]
fn _host_apply(py: Python<'_>, head: &PyHostHead, args: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    Ok(apply(py, &head.head, args)?.unbind())
}

/// Whether a host atom or application occurs in `expr`, and whether `expr`
/// itself is one.
#[pyfunction]
fn _host_contains(expr: &Bound<'_, PyAny>) -> PyResult<(bool, bool)> {
    let e = to_expr(expr)?;
    Ok((e.has_host(), e.is_host()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyHostClass>()?;
    module.add_class::<PyHostHead>()?;
    module.add_function(wrap_pyfunction!(_host_register, module)?)?;
    module.add_function(wrap_pyfunction!(_host_class, module)?)?;
    module.add_function(wrap_pyfunction!(_host_atom, module)?)?;
    module.add_function(wrap_pyfunction!(_host_head, module)?)?;
    module.add_function(wrap_pyfunction!(_host_apply, module)?)?;
    module.add_function(wrap_pyfunction!(_host_contains, module)?)?;
    Ok(())
}
