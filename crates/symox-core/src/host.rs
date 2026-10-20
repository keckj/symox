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
}

impl HostFormat {
    pub fn name(self) -> &'static str {
        match self {
            HostFormat::Str => "str",
            HostFormat::Repr => "repr",
            HostFormat::Latex => "latex",
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
    /// The first two fields of the class `class_key()`, when they differ
    /// from those of a function head.
    fn class_key(&self) -> Option<(u8, u32)> { None }
    /// `_eval_derivative` of an application of an expression class, or
    /// `None` to hold the derivative.
    fn derivative(&self, _expr: &Expr, _x: &Expr) -> Option<Expr> { None }
    /// A total order of handles of one class that are not `same`, so that
    /// native ordering stays consistent with equality.
    fn order(&self, _other: &dyn HostClass) -> std::cmp::Ordering { std::cmp::Ordering::Equal }
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
