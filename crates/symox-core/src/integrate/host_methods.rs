//! Integration methods provided by the host language, and the proven
//! non-elementary integrals they produce.
//!
//! The Python bindings install ports of the reference `risch_integrate`,
//! `manualintegrate`, `heurisch` and the term-loop stages of
//! `Integral._eval_integral`. The default chain uses the Risch algorithm
//! and then the reference term loop when the native rules give no
//! antiderivative; the `risch`, `manual` and `heurisch` options of
//! `integrate` select or exclude these methods.
//!
//! A non-elementary integral is recorded by node identity: only the nodes
//! built from a completed Risch proof (or by the explicit
//! `NonElementaryIntegral` constructor) are reported as such, and an equal
//! `Integral` built elsewhere is not.

use std::cell::Cell;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

use rustc_hash::FxHashMap;

use super::Conds;
use crate::expr::{Expr, Func, Kind, WeakExpr};

/// A method of the reference integrator implemented by the host.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HostMethod {
    /// `risch_integrate(f, x, separate_integral=True)`.
    Risch,
    /// The `manualintegrate` stage of the default term loop.
    Manual,
    /// The stages of the default term loop before the Meijer G method, on
    /// one term (`ratint` to `heurisch`).
    Term,
    /// `_eval_integral` with `manual=True` (the other methods off).
    ManualOnly,
    /// `_eval_integral` with `heurisch=True` (the other methods off).
    HeurischOnly,
}

/// The answer of a host method.
pub enum HostOutcome {
    /// The method does not apply (NotImplementedError, or no result).
    Unavailable,
    /// An antiderivative (for `manualintegrate`, possibly with unevaluated
    /// integrals in it).
    Value(Expr),
    /// `risch_integrate(separate_integral=True)` proved the second integrand
    /// non-elementary: the integral is `first + Integral(second, x)`.
    NonElementary(Expr, Expr),
}

pub type IntegrationHost = fn(HostMethod, &Expr, &Expr, Conds) -> HostOutcome;

static HOST: OnceLock<IntegrationHost> = OnceLock::new();

/// Install the host implementation of the Risch, manual and heurisch
/// methods.
pub fn register_integration_host(host: IntegrationHost) {
    let _ = HOST.set(host);
}

pub(crate) fn host() -> Option<IntegrationHost> {
    HOST.get().copied()
}

/// The `risch`, `manual` and `heurisch` options in force (None: default).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Methods {
    pub risch: Option<bool>,
    pub manual: Option<bool>,
    pub heurisch: Option<bool>,
}

thread_local! {
    static METHODS: Cell<Methods> = const { Cell::new(Methods { risch: None, manual: None, heurisch: None }) };
    static CONDS: Cell<Conds> = const { Cell::new(Conds::Piecewise) };
}

pub fn methods() -> Methods {
    METHODS.with(|m| m.get())
}

pub(crate) fn conds() -> Conds {
    CONDS.with(|c| c.get())
}

struct Restore<T: Copy + 'static>(&'static std::thread::LocalKey<Cell<T>>, T);
impl<T: Copy + 'static> Drop for Restore<T> {
    fn drop(&mut self) {
        let value = self.1;
        self.0.with(|cell| cell.set(value));
    }
}

fn scoped<T: Copy + 'static, R>(key: &'static std::thread::LocalKey<Cell<T>>, value: T, body: impl FnOnce() -> R) -> R {
    let previous = key.with(|cell| cell.replace(value));
    let _restore = Restore(key, previous);
    body()
}

/// Run `body` with the given method options (as `integrate(..., risch=...,
/// manual=..., heurisch=...)` and `conds`).
pub fn with_methods<R>(methods: Methods, conds: Conds, body: impl FnOnce() -> R) -> R {
    scoped(&METHODS, methods, || scoped(&CONDS, conds, body))
}

/// Run `body` as `i.doit(risch=False)`: the Risch method excluded and the
/// other options at their defaults.
pub(crate) fn doit_without_risch<R>(body: impl FnOnce() -> R) -> R {
    with_methods(Methods { risch: Some(false), manual: None, heurisch: None }, Conds::Piecewise, body)
}

struct Registry {
    nodes: FxHashMap<usize, WeakExpr>,
    live_after_purge: usize,
}

static REGISTRY: Mutex<Option<Registry>> = Mutex::new(None);
static ANY_MARKED: AtomicBool = AtomicBool::new(false);

/// Record the integral node `e` as a proven (or explicitly constructed)
/// non-elementary integral.
pub fn mark_nonelementary(e: &Expr) {
    if !matches!(e.kind(), Kind::Func(Func::Integral, _)) {
        return;
    }
    let Ok(mut guard) = REGISTRY.lock() else { return };
    let registry = guard.get_or_insert_with(|| Registry { nodes: FxHashMap::default(), live_after_purge: 0 });
    if registry.nodes.len() > 2 * registry.live_after_purge + 64 {
        registry.nodes.retain(|_, weak| weak.upgrade().is_some());
        registry.live_after_purge = registry.nodes.len();
    }
    registry.nodes.insert(e.dag_id(), e.downgrade());
    ANY_MARKED.store(true, Ordering::Release);
}

/// Whether `e` is a node recorded by [`mark_nonelementary`].
pub fn is_nonelementary(e: &Expr) -> bool {
    if !ANY_MARKED.load(Ordering::Acquire) || !matches!(e.kind(), Kind::Func(Func::Integral, _)) {
        return false;
    }
    let Ok(guard) = REGISTRY.lock() else { return false };
    guard.as_ref().and_then(|r| r.nodes.get(&e.dag_id())).is_some_and(|weak| weak.is(e))
}

/// The method selected by a `True` option, if any.
pub(crate) fn selected() -> Option<HostMethod> {
    let m = methods();
    if m.risch == Some(true) {
        Some(HostMethod::Risch)
    } else if m.manual == Some(true) {
        Some(HostMethod::ManualOnly)
    } else if m.heurisch == Some(true) {
        Some(HostMethod::HeurischOnly)
    } else {
        None
    }
}

/// `xreplace(e, m)` keeping the non-elementary record of the integrals in
/// `e` on their images.
pub(crate) fn xreplace_marked(e: &Expr, m: &FxHashMap<Expr, Expr>) -> Expr {
    fn marked_below(e: &Expr) -> bool {
        let mut found = false;
        crate::traverse::preorder(e, &mut |s| found |= is_nonelementary(s));
        found
    }
    fn go(e: &Expr, m: &FxHashMap<Expr, Expr>) -> Expr {
        if let Some(r) = m.get(e) {
            return r.clone();
        }
        if is_nonelementary(e) {
            let r = crate::subs::xreplace(e, m);
            mark_nonelementary(&r);
            return r;
        }
        if !marked_below(e) {
            return crate::subs::xreplace(e, m);
        }
        crate::simptools::map_args(e, &mut |a| go(a, m))
    }
    if !ANY_MARKED.load(Ordering::Acquire) {
        return crate::subs::xreplace(e, m);
    }
    go(e, m)
}

/// A held `Integral(f, x)` recorded as non-elementary.
pub(crate) fn nonelementary_integral(f: &Expr, x: &Expr) -> Expr {
    let e = super::integral_unevaluated(f, x);
    mark_nonelementary(&e);
    e
}
