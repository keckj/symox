//! The variable/order plan of `Derivative.__new__` for a matrix
//! expression, an explicit matrix or an array derived by scalar Symbols
//! (the symbolic-order derivatives of `tensor/array/array_derivatives.py`):
//! the merged `(Symbol, count)` pairs of the arguments, the zero exit and
//! the canonical order (`Derivative._sort_variable_count`). The caller
//! derives the integer orders and builds the held derivative.

use crate::{Expr, Kind};

/// One argument of `diff(*variables)`.
pub enum OrderItem {
    /// An int or Integer following its variable (`diff(x, 2)`).
    Order(num_bigint::BigInt),
    /// A variable with its count (1 for a bare variable).
    Variable(Expr, Expr),
}

/// `_parse`: `Ok(None)` when an argument is not a scalar Symbol with an
/// admissible order (another path applies), `Err` for a negative order.
pub fn parse_orders(items: Vec<OrderItem>) -> Result<Option<Vec<(Expr, Expr)>>, &'static str> {
    let mut pairs: Vec<(Expr, Expr)> = Vec::new();
    let one = Expr::one();
    for item in items {
        let (v, count) = match item {
            OrderItem::Order(n) => {
                let Some(last) = pairs.last_mut() else { return Ok(None) };
                if last.1 != one || n < num_bigint::BigInt::from(0) {
                    return Ok(None);
                }
                if n == num_bigint::BigInt::from(0) {
                    pairs.pop();
                } else {
                    last.1 = Expr::num(crate::number::Number::from_bigint(n));
                }
                continue;
            }
            OrderItem::Variable(v, count) => (v, count),
        };
        if !matches!(v.kind(), Kind::Sym(_)) || v.has_matrix() {
            return Ok(None);
        }
        if let Kind::Num(n) = count.kind() {
            if !n.is_integer() {
                return Ok(None);
            }
        }
        if count.has_matrix() || count.has_host() {
            return Ok(None);
        }
        if crate::assumptions::fact(&count, "negative") == Some(true) {
            return Err("order of differentiation must be nonnegative");
        }
        match pairs.last_mut() {
            Some(last) if last.0 == v => {
                let total = crate::ops::add(&[last.1.clone(), count]);
                if total.is_zero() {
                    pairs.pop();
                } else {
                    last.1 = total;
                }
            }
            _ if !count.is_zero() => pairs.push((v, count)),
            _ => {}
        }
    }
    Ok(Some(pairs))
}

/// Whether a count is an Integer (`_symbolic_orders`).
pub fn is_integer_order(count: &Expr) -> bool {
    matches!(count.kind(), Kind::Num(n) if n.is_integer())
}

/// The plan of `_evaluate(expr, pairs, cls)`.
pub enum OrderPlan {
    /// No variables: the expression itself.
    Unchanged,
    /// A positive order of a variable `expr` does not contain.
    Zero,
    /// The canonical pairs, and the index of the first symbolic order
    /// (`ordered.len()` when every order is an Integer): the integer
    /// prefix is derived, the rest held.
    Derive(Vec<(Expr, Expr)>, usize),
}

/// The zero exit with the free symbols `free` of the expression, then the
/// pairs sorted by `default_sort_key` of their variables with equal
/// neighbours merged.
pub fn plan_orders(pairs: Vec<(Expr, Expr)>, free: &[Expr]) -> OrderPlan {
    if pairs.is_empty() {
        return OrderPlan::Unchanged;
    }
    if pairs.iter().any(|(v, c)| crate::assumptions::fact(c, "positive") == Some(true) && !free.contains(v)) {
        return OrderPlan::Zero;
    }
    let mut sorted = pairs;
    sorted.sort_by(|a, b| crate::sortkey::cmp_default(&a.0, &b.0));
    let mut ordered: Vec<(Expr, Expr)> = Vec::with_capacity(sorted.len());
    for (v, c) in sorted {
        match ordered.last_mut() {
            Some(last) if last.0 == v => last.1 = crate::ops::add(&[last.1.clone(), c]),
            _ => ordered.push((v, c)),
        }
    }
    let split = ordered.iter().position(|(_, c)| !is_integer_order(c)).unwrap_or(ordered.len());
    OrderPlan::Derive(ordered, split)
}
