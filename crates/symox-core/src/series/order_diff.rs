//! `Order._eval_derivative`: `O(f, *args).diff(x)` is `O(f.diff(x), *args)`
//! (`O(x**2).diff(x)` is `O(x)`), for an order term that holds `x` in its
//! expression or among its variables; otherwise the derivative is `0`.

use crate::expr::{Expr, Func, Kind};

/// The derivative of an order term, `None` for any other expression.
pub fn order_derivative(e: &Expr, x: &Expr) -> Option<Expr> {
    if !matches!(e.kind(), Kind::Func(Func::Order, _)) || !x.is_symbol() {
        return None;
    }
    let (inner, vars) = crate::series::order::order_parts(e)?;
    let free = crate::traverse::free_symbols(&inner);
    if !free.contains(x) && !vars.iter().any(|(v, _)| v == x) {
        return Some(Expr::zero());
    }
    Some(crate::series::order::order(&crate::diff::diff(&inner, x), &vars))
}
