//! Python subclasses of native atoms, function heads and expressions.
//!
//! A user subclass of `Symbol`, `Dummy`, `Function`, `Expr`, `Atom` or
//! `Basic` is a real subtype of `_core.Expr`. Native nodes carry an opaque
//! handle to its class, so every expression rebuilt by native code returns
//! instances of that class again.
//!
//! Each node owns the Python instance bound to it: its handle holds a strong
//! reference, so the instance (identity and attributes) lives exactly as
//! long as some expression, Python object or native container references
//! the node. The instance itself owns its node, so this forms a reference
//! cycle that the collector sees through `traverse`: a Python expression
//! object reports the instances of the host nodes that are reachable only
//! through it, and `clear` releases them.
use std::cmp::Ordering as Order;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use pyo3::impl_::pyclass_init::PyObjectInit;
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use pyo3::types::{PyDict, PyTuple, PyType};
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
const SUBS: u32 = 512;
const ASSUME: u32 = 1024;
const EVALF: u32 = 2048;
const FUNC: u32 = 4096;
const NONCOMMUTATIVE: u32 = 8192;
const UNCOMMUTED: u32 = 16384;
const COMMUTATIVE_HOOK: u32 = 32768;

fn may_not_commute(flags: u32) -> bool { flags & (NONCOMMUTATIVE | UNCOMMUTED | COMMUTATIVE_HOOK) != 0 }

static NEXT_CLASS: AtomicU64 = AtomicU64::new(1);
static NEXT_HANDLE: AtomicU64 = AtomicU64::new(1);

/// The Python `HostTree` type: plain expression objects owning host nodes.
static HOST_TREE: GILOnceCell<Py<PyType>> = GILOnceCell::new();

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

/// A Python instance bound to a node it owns. `contents` is the address of
/// the instance's native payload, which identifies it during collection.
#[derive(Debug)]
struct Binding {
    node: WeakExpr,
    instance: Py<PyAny>,
    contents: usize,
}

/// The handle of one host atom, shared by every node of that atom: they all
/// denote one Python instance.
#[derive(Debug)]
struct AtomHost {
    class: Arc<ClassInfo>,
    serial: u64,
    bound: Mutex<Option<Binding>>,
}

/// The handle of a host head. The head registered by the class binds no
/// instance; every application gets its own handle, bound to the instance
/// of its node.
#[derive(Debug)]
struct HeadHost {
    class: Arc<ClassInfo>,
    serial: u64,
    bound: Mutex<Option<Binding>>,
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
/// constructor (which may evaluate or reset attributes).
fn instantiate<'py>(py: Python<'py>, class: &Bound<'py, PyType>, e: Expr) -> PyResult<Bound<'py, PyAny>> {
    if !class.is_subclass_of::<PyExprNode>()? {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "{} is not a subtype of the native Expr type", class.name()?)));
    }
    // SAFETY: `class` is a subtype of the native Expr type, whose layout
    // the initializer fills; the new reference is owned by the result.
    unsafe {
        let object = PyClassInitializer::from(PyExprNode { e }).into_new_object(py, class.as_type_ptr())?;
        Ok(Bound::from_owned_ptr(py, object))
    }
}

/// The instance bound in `slot`, when its node is `e` (any node for atoms).
fn bound_instance<'py>(py: Python<'py>, slot: &Mutex<Option<Binding>>, e: Option<&Expr>) -> Option<Bound<'py, PyAny>> {
    let guard = slot.lock().ok()?;
    let binding = guard.as_ref().filter(|binding| e.is_none_or(|e| binding.node.is(e)))?;
    Some(binding.instance.bind(py).clone())
}

/// Bind `value`, a new instance owning node `e`, in `slot` unless an
/// instance is bound already. Returns the instance of `e`: the one bound
/// before when it denotes `e` (always for atoms), else `value`.
fn bind<'py>(slot: &Mutex<Option<Binding>>, e: &Expr, value: Bound<'py, PyAny>, atom: bool) -> Bound<'py, PyAny> {
    let py = value.py();
    let contents = match value.downcast::<PyExprNode>() {
        Ok(native) => native.get() as *const PyExprNode as usize,
        Err(_) => return value,
    };
    // No Python object is created or released while the slot is locked.
    let previous = {
        let Ok(mut guard) = slot.lock() else { return value };
        match guard.as_ref() {
            Some(binding) if atom || binding.node.is(e) => Some(binding.instance.bind(py).clone()),
            Some(_) => None,
            None => {
                *guard = Some(Binding { node: e.downgrade(), instance: value.clone().unbind(), contents });
                None
            }
        }
    };
    previous.unwrap_or(value)
}

/// The binding slot of a host handle of this module.
fn slot(handle: &Arc<dyn HostClass>) -> Option<&Mutex<Option<Binding>>> {
    let handle = handle.as_any();
    handle.downcast_ref::<AtomHost>().map(|host| &host.bound)
        .or_else(|| handle.downcast_ref::<HeadHost>().map(|host| &host.bound))
}

/// The handles whose bound instances are reachable only through the Python
/// expression object with payload address `me`, which owns `root`.
///
/// A node belongs to the region of that object when every strong reference
/// to it comes from the object or from nodes of the region, except for the
/// reference held by its own bound instance (when that is another object).
/// A handle belongs to it when every reference to the handle comes from
/// nodes of the region. References this walk cannot attribute (native
/// containers, caches, other objects) keep instances out of the region, so
/// the collector may defer, never break, a live binding.
fn exclusive(me: usize, root: &Expr) -> Vec<&Arc<dyn HostClass>> {
    let extra = |e: &Expr| -> usize {
        e.host_handle().and_then(slot).and_then(|slot| slot.try_lock().ok())
            .is_some_and(|guard| guard.as_ref().is_some_and(|binding| binding.node.is(e) && binding.contents != me))
            as usize
    };
    let mut handles: Vec<&Arc<dyn HostClass>> = Vec::new();
    let mut holders: rustc_hash::FxHashMap<usize, usize> = Default::default();
    let mut owners: rustc_hash::FxHashMap<usize, usize> = Default::default();
    let mut stack = Vec::new();
    if root.strong_count() == 1 + extra(root) { stack.push(root); }
    while let Some(node) = stack.pop() {
        if let Some(handle) = node.host_handle() {
            let held = holders.entry(Arc::as_ptr(handle) as *const () as usize).or_insert(0);
            *held += 1;
            if *held == Arc::strong_count(handle) { handles.push(handle); }
        }
        node.for_each_owned(&mut |child| {
            if !child.has_host() { return; }
            let owned = owners.entry(child.dag_id()).or_insert(0);
            *owned += 1;
            if *owned + extra(child) == child.strong_count() { stack.push(child); }
        });
    }
    handles
}

/// Report to the collector the host instances reachable only through the
/// Python expression object with payload address `me`.
pub(crate) fn traverse(me: usize, root: &Expr, visit: &pyo3::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
    for handle in exclusive(me, root) {
        let Some(guard) = slot(handle).and_then(|slot| slot.try_lock().ok()) else { continue };
        if let Some(binding) = guard.as_ref() { visit.call(&binding.instance)?; }
    }
    Ok(())
}

/// Release the instances reported by `traverse`: the collector found them
/// unreachable together with this object.
pub(crate) fn clear(me: usize, root: &Expr) {
    let released = exclusive(me, root).into_iter()
        .filter_map(|handle| slot(handle)?.try_lock().ok()?.take())
        .collect::<Vec<_>>();
    drop(released);
}

/// The object of `e`, which contains host nodes, as an instance of the
/// nominal node class `class` (a HostTree subtype).
pub(crate) fn instance_of(py: Python<'_>, class: &Bound<'_, PyType>, e: Expr) -> PyResult<PyObject> {
    Ok(instantiate(py, class, e)?.unbind())
}

/// Whether `class` is the registered plain `HostTree` type.
pub(crate) fn is_tree_type(py: Python<'_>, class: &Bound<'_, PyType>) -> bool {
    HOST_TREE.get(py).is_some_and(|tree| class.is(tree.bind(py)))
}

/// The plain Python expression object of `e`, which contains host nodes.
pub(crate) fn tree_object(py: Python<'_>, e: Expr) -> PyResult<PyObject> {
    match HOST_TREE.get(py) {
        Some(tree) => Ok(instantiate(py, tree.bind(py), e)?.unbind()),
        None => Ok(Bound::new(py, PyExprNode { e })?.into_any().unbind()),
    }
}

thread_local! {
    /// The first exception raised by a host callback of this thread that
    /// was not yet raised to the caller.
    static PENDING: std::cell::RefCell<Option<PyErr>> = const { std::cell::RefCell::new(None) };
}

/// Record an exception raised by a host callback. Native algorithms cannot
/// unwind through their callers, so they continue with the native fallback
/// and the exception is raised when the result returns to Python (see
/// `raise_pending`), as the reference raises it out of the computation.
fn report(_py: Python<'_>, error: PyErr, _context: Option<&Bound<'_, PyAny>>) {
    let extra = PENDING.with(|pending| {
        let mut pending = pending.borrow_mut();
        if pending.is_none() { *pending = Some(error); None } else { Some(error) }
    });
    symox::host::note_failure();
    drop(extra);
}

/// Raise the exception of a host callback that failed during the native
/// computation now returning to Python.
#[inline]
pub(crate) fn raise_pending() -> PyResult<()> {
    if !symox::host::failure_pending() { return Ok(()); }
    symox::host::clear_failure();
    match PENDING.with(|pending| pending.borrow_mut().take()) {
        Some(error) => Err(error),
        None => Ok(()),
    }
}

/// `value`, unless a host callback failed while computing it.
#[inline]
pub(crate) fn checked<T>(value: T) -> PyResult<T> {
    raise_pending()?;
    Ok(value)
}

/// Whether the current thread holds the GIL. Native workers comparing host
/// nodes must not wait for it while the calling thread holds it.
fn holds_gil() -> bool {
    extern "C" { fn PyGILState_Check() -> std::os::raw::c_int; }
    // SAFETY: PyGILState_Check (stable ABI) may be called from any thread.
    unsafe { PyGILState_Check() == 1 }
}

impl AtomHost {
    /// The instance of this atom, when one is bound.
    fn instance<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyAny>> {
        bound_instance(py, &self.bound, None)
    }

    /// The instance of this atom, created for node `e` when none is bound.
    fn object<'py>(&self, py: Python<'py>, e: &Expr) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = self.instance(py) { return Ok(value); }
        let value = instantiate(py, self.class.class.bind(py), e.clone())?;
        Ok(bind(&self.bound, e, value, true))
    }
}

impl HeadHost {
    fn fresh(&self) -> HeadHost {
        HeadHost { class: self.class.clone(), serial: NEXT_HANDLE.fetch_add(1, Ordering::Relaxed), bound: Mutex::new(None) }
    }

    /// The instance bound to this handle.
    fn instance<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyAny>> {
        bound_instance(py, &self.bound, None)
    }

    /// The instance of node `e`. A node sharing the handle of another bound
    /// node gets a new, unbound instance.
    fn object<'py>(&self, py: Python<'py>, e: &Expr) -> PyResult<Bound<'py, PyAny>> {
        if let Some(value) = bound_instance(py, &self.bound, Some(e)) { return Ok(value); }
        let value = instantiate(py, self.class.class.bind(py), e.clone())?;
        Ok(bind(&self.bound, e, value, false))
    }
}

/// `_hashable_content()` of the instances of two handles of one class, or
/// `None` when either is unavailable or raised (the error is reported).
fn contents<'py>(py: Python<'py>, a: Option<Bound<'py, PyAny>>, b: Option<Bound<'py, PyAny>>) -> Option<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
    let content = |value: Bound<'py, PyAny>| value.call_method0("_hashable_content")
        .map_err(|error| report(py, error, Some(&value))).ok();
    Some((content(a?)?, content(b?)?))
}

/// SymPy equality of two instances of a class defining `_hashable_content`.
fn same_content<'py>(py: Python<'py>, a: Option<Bound<'py, PyAny>>, b: Option<Bound<'py, PyAny>>) -> bool {
    contents(py, a, b).is_some_and(|(a, b)| a.eq(&b).unwrap_or_else(|error| { report(py, error, None); false }))
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

/// An optional expression returned by a host hook (`None` declines).
fn optional_expr(value: Bound<'_, PyAny>) -> PyResult<Option<Expr>> {
    if value.is_none() { Ok(None) } else { to_expr(&value).map(Some) }
}

/// SymPy `_eval_subs(old, new)` of the instance `value`.
fn call_subs(py: Python<'_>, value: &Bound<'_, PyAny>, old: &Expr, new: &Expr) -> Option<Expr> {
    let result = crate::scalar_object(py, old.clone())
        .and_then(|old| Ok((old, crate::scalar_object(py, new.clone())?)))
        .and_then(|(old, new)| value.call_method1(pyo3::intern!(py, "_eval_subs"), (old, new)))
        .and_then(optional_expr);
    match result {
        Ok(value) => value,
        Err(error) => { report(py, error, Some(value)); None }
    }
}

/// The answer of the `_eval_is_<fact>` handlers of the instance `value`
/// (SymPy `_ask`), or `None` to keep the native facts.
fn call_fact(py: Python<'_>, value: &Bound<'_, PyAny>, name: &str) -> Option<Option<bool>> {
    let result = py.import("symox._host").and_then(|module| module.getattr("ask")?.call1((value, name)))
        .and_then(|answer| answer.extract::<Option<bool>>());
    match result {
        Ok(Some(answer)) => Some(Some(answer)),
        Ok(None) => None,
        Err(error) => { report(py, error, Some(value)); None }
    }
}

/// SymPy `_eval_evalf(prec)` of the instance `value` of node `e`.
fn call_evalf(py: Python<'_>, value: &Bound<'_, PyAny>, e: &Expr, precision: u32) -> Option<Expr> {
    match value.call_method1(pyo3::intern!(py, "_eval_evalf"), (precision,)).and_then(optional_expr) {
        // A hook returning its own node cannot be evaluated further.
        Ok(result) => result.filter(|result| result != e),
        Err(error) => { report(py, error, Some(value)); None }
    }
}

impl HostClass for AtomHost {
    fn identity(&self) -> u64 { self.class.id }
    fn class_name(&self) -> &str { &self.class.name }
    fn current_name(&self, atom: &Expr) -> Option<String> {
        if self.class.flags & RENAMES == 0 { return None; }
        let Kind::Sym(symbol) = atom.kind() else { return None };
        Python::with_gil(|py| {
            let value = self.instance(py)?;
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
    fn substitutes(&self) -> bool { self.class.flags & SUBS != 0 }
    fn subs(&self, expr: &Expr, old: &Expr, new: &Expr) -> Option<Expr> {
        if self.class.flags & SUBS == 0 { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => call_subs(py, &value, old, new),
            Err(error) => { report(py, error, None); None }
        })
    }
    fn fact(&self, expr: &Expr, name: &str) -> Option<Option<bool>> {
        if self.class.flags & ASSUME == 0 { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => call_fact(py, &value, name),
            Err(error) => { report(py, error, None); None }
        })
    }
    fn may_not_commute(&self) -> bool { may_not_commute(self.class.flags) }
    fn commutative(&self, expr: &Expr) -> Option<Option<bool>> {
        let flags = self.class.flags;
        if flags & NONCOMMUTATIVE != 0 { return Some(Some(false)); }
        if flags & UNCOMMUTED != 0 { return Some(None); }
        if flags & COMMUTATIVE_HOOK == 0 { return None; }
        Some(self.fact(expr, "commutative").flatten())
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
        if self.class.flags & (EVAL | FUNC) == 0 { return None; }
        Python::with_gil(|py| {
            // A rebuild calls `self.func(*args)` as SymPy's xreplace and
            // subs do: an overridden `func` of the instance rebuilt, else
            // the class, whose constructor evaluates.
            let class = self.class.class.bind(py).clone().into_any();
            let callee = match self.instance(py) {
                Some(instance) if self.class.flags & FUNC != 0 => match instance.getattr(pyo3::intern!(py, "func")) {
                    Ok(callee) => callee,
                    Err(error) => { report(py, error, Some(&instance)); return None; }
                },
                _ if self.class.flags & EVAL == 0 => return None,
                _ => class,
            };
            let result = args.iter().map(|arg| crate::scalar_object(py, arg.clone())).collect::<PyResult<Vec<_>>>()
                .and_then(|args| PyTuple::new(py, args))
                .and_then(|args| callee.call1(args))
                .and_then(|value| to_expr(&value));
            match result {
                Ok(value) => Some(value),
                Err(error) => { report(py, error, Some(&callee)); None }
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
    fn applied(&self) -> Option<Arc<dyn HostClass>> { Some(Arc::new(self.fresh())) }
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
    fn differentiates(&self) -> bool { self.class.flags & DERIV != 0 }
    fn substitutes(&self) -> bool { self.class.flags & SUBS != 0 }
    fn subs(&self, expr: &Expr, old: &Expr, new: &Expr) -> Option<Expr> {
        if self.class.flags & SUBS == 0 { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => call_subs(py, &value, old, new),
            Err(error) => { report(py, error, None); None }
        })
    }
    fn fact(&self, expr: &Expr, name: &str) -> Option<Option<bool>> {
        if self.class.flags & ASSUME == 0 { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => call_fact(py, &value, name),
            Err(error) => { report(py, error, None); None }
        })
    }
    fn may_not_commute(&self) -> bool { may_not_commute(self.class.flags) }
    fn commutative(&self, expr: &Expr) -> Option<Option<bool>> {
        let flags = self.class.flags;
        if flags & NONCOMMUTATIVE != 0 { return Some(Some(false)); }
        if flags & UNCOMMUTED != 0 { return Some(None); }
        if flags & COMMUTATIVE_HOOK == 0 { return None; }
        Some(self.fact(expr, "commutative").flatten())
    }
    fn evaluates_numerically(&self) -> bool { self.class.flags & EVALF != 0 }
    fn evalf(&self, expr: &Expr, precision: u32) -> Option<Expr> {
        if self.class.flags & EVALF == 0 { return None; }
        Python::with_gil(|py| match self.object(py, expr) {
            Ok(value) => call_evalf(py, &value, expr, precision),
            Err(error) => { report(py, error, None); None }
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
    let handle = Arc::new(host.fresh());
    let e = Expr::from_kind(Kind::Func(Func::Undefined(head.clone().with_host(handle.clone())), args));
    let value = instantiate(py, &class, e.clone())?;
    Ok(bind(&handle.bound, &e, value, false))
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
        class: host.info.clone(), serial: NEXT_HANDLE.fetch_add(1, Ordering::Relaxed), bound: Mutex::new(None),
    });
    let e = Expr::from_kind(Kind::Sym(symbol.with_host(handle.clone())));
    let value = instantiate(py, host.info.class.bind(py), e.clone())?;
    Ok(bind(&handle.bound, &e, value, true).unbind())
}

/// The native head defined by a host function or expression class.
#[pyfunction]
#[pyo3(signature = (host, name, assumptions, nargs=None))]
fn _host_head(host: &PyHostClass, name: &str, assumptions: &Bound<'_, PyDict>, nargs: Option<Vec<usize>>) -> PyResult<PyHostHead> {
    let handle = HeadHost { class: host.info.clone(), serial: 0, bound: Mutex::new(None) };
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
