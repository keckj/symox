//! Opaque handles to host-language classes carried by atoms and function heads.
//!
//! A host (for instance a Python subclass of `Symbol` or `Function`) keeps its
//! own identity and behaviour. Native identity includes `identity()`, so two
//! atoms with the same name but different host classes are distinct, as SymPy
//! compares `type(self)` before `_hashable_content`. Every callback is
//! optional: returning `None` keeps the native behaviour.

use std::fmt;

use crate::expr::Expr;

/// Printer families that a host class may customize.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum HostFormat {
    Str,
    Repr,
    Latex,
    /// A one-line unicode pretty form of an atom or a nullary application.
    Pretty,
    /// A one-line ASCII pretty form of an atom or a nullary application.
    PrettyAscii,
}

impl HostFormat {
    pub fn name(self) -> &'static str {
        match self {
            HostFormat::Str => "str",
            HostFormat::Repr => "repr",
            HostFormat::Latex => "latex",
            HostFormat::Pretty => "pretty",
            HostFormat::PrettyAscii => "pretty_ascii",
        }
    }
}

pub trait HostClass: Send + Sync + fmt::Debug {
    /// Stable identity of the host class. Equal identities denote one class.
    fn identity(&self) -> u64;
    /// Host class name, used to order atoms deterministically.
    fn class_name(&self) -> &str;
    /// The current host-side name of an atom, when it was renamed.
    fn current_name(&self, _atom: &Expr) -> Option<String> { None }
    /// A complete rendering of `expr` (an atom or an application of this head).
    fn print(&self, _expr: &Expr, _format: HostFormat) -> Option<String> { None }
    /// Host evaluation of an application of this head (SymPy `eval`).
    /// `None` keeps the application unevaluated.
    fn evaluate(&self, _args: &[Expr]) -> Option<Expr> { None }
    /// Partial derivative of an application with respect to argument
    /// `index`, counted from 1 (SymPy `fdiff`).
    fn fdiff(&self, _expr: &Expr, _index: usize) -> Option<Expr> { None }
    /// Whether two handles of one class denote equal objects, for classes
    /// whose identity includes host state (SymPy `_hashable_content`).
    /// Only called when both identities are equal.
    fn same(&self, _other: &dyn HostClass) -> bool { true }
    /// The applications of this head are plain expressions (a SymPy `Expr`
    /// or `Basic` subclass that is not a `Function`): no chain rule, and
    /// `\operatorname{name}\left(args\right)` in LaTeX.
    fn expression(&self) -> bool { false }
    /// The applications of this head are atoms (SymPy `Atom` subclasses).
    fn atomic(&self) -> bool { false }
    /// The atoms of this class are numbers (SymPy `Number` and
    /// `NumberSymbol` subclasses): `is_number` is True.
    fn number(&self) -> bool { false }
    /// The first two fields of the class `class_key()`, when they differ
    /// from those of a function head.
    fn class_key(&self) -> Option<(u8, u32)> { None }
    /// `_eval_derivative` of an application of an expression class, or of
    /// a function class that `differentiates`; `None` holds the derivative.
    fn derivative(&self, _expr: &Expr, _x: &Expr) -> Option<Expr> { None }
    /// A total order of handles of one class that are not `same`, so that
    /// native ordering stays consistent with equality.
    fn order(&self, _other: &dyn HostClass) -> std::cmp::Ordering { std::cmp::Ordering::Equal }
    /// A fresh handle of this class for a new application of this head
    /// built natively, so that each application binds its own host object.
    fn applied(&self) -> Option<std::sync::Arc<dyn HostClass>> { None }
    /// Whether `subs` must consult `subs` on this node (SymPy `_eval_subs`),
    /// including where the replaced pattern does not occur.
    fn substitutes(&self) -> bool { false }
    /// SymPy `_eval_subs(old, new)` of this atom or application; `None`
    /// falls back to native substitution of its arguments.
    fn subs(&self, _expr: &Expr, _old: &Expr, _new: &Expr) -> Option<Expr> { None }
    /// The host's answer to the assumption `name` of this atom or
    /// application (SymPy `_ask` over its `_eval_is_<fact>` handlers):
    /// `Some(answer)` when the host decides the query, `None` to keep the
    /// native facts of its class.
    fn fact(&self, _expr: &Expr, _name: &str) -> Option<Option<bool>> { None }
    /// Whether applications of a function head differentiate through
    /// `derivative` (an overridden SymPy `_eval_derivative`) instead of the
    /// chain rule over `fdiff`.
    fn differentiates(&self) -> bool { false }
    /// Whether applications of a function head are elements of a named
    /// symbolic array (SymPy `ArrayElement`, `_diff_wrt`): the first argument
    /// names the array, the others are the indices.
    fn array_element(&self) -> bool { false }
    /// Whether numerical evaluation consults `evalf` (SymPy `_eval_evalf`).
    fn evaluates_numerically(&self) -> bool { false }
    /// SymPy `_eval_evalf(prec)` of an application, at binary precision
    /// `precision`: a numerical value, or `None` to keep it unevaluated.
    fn evalf(&self, _expr: &Expr, _precision: u32) -> Option<Expr> { None }
    /// Whether instances of this class may fail to commute: the class
    /// declares `is_commutative` False or None, has no commutative default
    /// (SymPy `Expr`, `Atom` and `Basic` subclasses), or decides it with an
    /// `_eval_is_commutative` handler.
    fn may_not_commute(&self) -> bool { false }
    /// SymPy `is_commutative` of this atom or application when its class
    /// decides it, `None` when it follows from the native facts and, for a
    /// function, from the arguments.
    fn commutative(&self, _expr: &Expr) -> Option<Option<bool>> { None }
    fn as_any(&self) -> &dyn std::any::Any;
}

/// Equality of optional host handles: one class, and equal host state.
pub(crate) fn same(a: Option<&std::sync::Arc<dyn HostClass>>, b: Option<&std::sync::Arc<dyn HostClass>>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => a.identity() == b.identity()
            && (std::ptr::addr_eq(std::sync::Arc::as_ptr(a), std::sync::Arc::as_ptr(b)) || a.same(&**b)),
        _ => false,
    }
}

/// Order of optional host handles with equal class keys (see `HostClass::order`).
pub(crate) fn order(a: Option<&std::sync::Arc<dyn HostClass>>, b: Option<&std::sync::Arc<dyn HostClass>>) -> std::cmp::Ordering {
    match (a, b) {
        (Some(a), Some(b)) if a.identity() == b.identity()
            && !std::ptr::addr_eq(std::sync::Arc::as_ptr(a), std::sync::Arc::as_ptr(b)) => a.order(&**b),
        _ => std::cmp::Ordering::Equal,
    }
}

/// The head handle of an application of a host expression class.
pub fn expression_head(e: &Expr) -> Option<&std::sync::Arc<dyn HostClass>> {
    match e.kind() {
        crate::expr::Kind::Func(crate::expr::Func::Undefined(head), _) => head.host().filter(|host| host.expression()),
        _ => None,
    }
}

/// The host rendering of an atom or application owned by a host class.
pub fn print(e: &Expr, format: HostFormat) -> Option<String> {
    if !e.has_host() { return None; }
    match e.kind() {
        crate::expr::Kind::Sym(symbol) => symbol.host()?.print(e, format),
        crate::expr::Kind::Func(crate::expr::Func::Undefined(head), _) => head.host()?.print(e, format),
        _ => None,
    }
}

/// The current display name of an atom, which a host object may rename.
pub fn atom_name<'a>(e: &Expr, symbol: &'a crate::symbol::Symbol) -> std::borrow::Cow<'a, str> {
    match symbol.host().and_then(|host| host.current_name(e)) {
        Some(name) => std::borrow::Cow::Owned(name),
        None => std::borrow::Cow::Borrowed(symbol.name()),
    }
}

/// The host class name of an atom or application, used by `srepr`.
pub fn class_name(e: &Expr) -> Option<&str> {
    match e.kind() {
        crate::expr::Kind::Sym(symbol) => symbol.host().map(|host| host.class_name()),
        crate::expr::Kind::Func(crate::expr::Func::Undefined(head), _) => head.host().map(|host| host.class_name()),
        _ => None,
    }
}

thread_local! {
    static FAILED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Record that a host callback raised. Results derived while a failure is
/// pending are returned to the host, which reports the failure, but are not
/// cached by native nodes.
pub fn note_failure() { FAILED.with(|failed| failed.set(true)); }
/// Whether a host failure is pending on this thread.
pub fn failure_pending() -> bool { FAILED.with(|failed| failed.get()) }
/// Clear the pending host failure of this thread.
pub fn clear_failure() { FAILED.with(|failed| failed.set(false)); }

/// The host's answer to assumption `name` of a host atom or application.
#[inline]
pub fn fact(e: &Expr, name: &str) -> Option<Option<bool>> {
    if !e.is_host() { return None; }
    e.host_handle()?.fact(e, name)
}
