//! Relationals (`Eq`, `Lt`, ...) and `Piecewise`.

use std::cmp::Ordering;

use crate::expr::{Expr, Kind, RelOp};
use crate::functions::{contains_i, is_numeric_value};
use crate::ops::sub;

/// Build `lhs op rhs`, evaluating to a boolean when decidable.
pub fn rel(op: RelOp, a: &Expr, b: &Expr) -> Expr {
    if let Some(v) = decide(op, a, b) {
        return Expr::boolean(v);
    }
    Expr::from_kind(Kind::Rel(op, a.clone(), b.clone()))
}

fn decide(op: RelOp, a: &Expr, b: &Expr) -> Option<bool> {
    if a == b {
        return Some(matches!(op, RelOp::Eq | RelOp::Le | RelOp::Ge));
    }
    let d = sub(a, b);
    if let Kind::Num(n) = d.kind() {
        let c = n.cmp_value(&crate::number::Number::ZERO)?;
        return Some(cmp_holds(op, c));
    }
    if is_numeric_value(&d) && !contains_i(&d) {
        let v = crate::evalf::eval_f64(&d)?;
        if v.abs() > 1e-10 {
            return Some(cmp_holds(op, v.partial_cmp(&0.0)?));
        }
    }
    if matches!(op, RelOp::Eq | RelOp::Ne) && is_numeric_value(&d) {
        if let Some(z) = crate::evalf::eval_c64(&d) {
            if z.norm() > 1e-10 {
                return Some(op == RelOp::Ne);
            }
        }
    }
    None
}

fn cmp_holds(op: RelOp, c: Ordering) -> bool {
    match op {
        RelOp::Eq => c == Ordering::Equal,
        RelOp::Ne => c != Ordering::Equal,
        RelOp::Lt => c == Ordering::Less,
        RelOp::Le => c != Ordering::Greater,
        RelOp::Gt => c == Ordering::Greater,
        RelOp::Ge => c != Ordering::Less,
    }
}

pub fn eq(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Eq, a, b)
}
pub fn ne(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Ne, a, b)
}
pub fn lt(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Lt, a, b)
}
pub fn le(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Le, a, b)
}
pub fn gt(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Gt, a, b)
}
pub fn ge(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Ge, a, b)
}

/// `Piecewise((e1, c1), (e2, c2), ...)`.
pub fn piecewise(pieces: Vec<(Expr, Expr)>) -> Expr {
    let mut out: Vec<(Expr, Expr)> = Vec::with_capacity(pieces.len());
    for (e, c) in pieces {
        match c.kind() {
            Kind::Bool(false) => continue,
            Kind::Bool(true) => {
                if out.is_empty() {
                    return e;
                }
                out.push((e, c));
                break;
            }
            _ => out.push((e, c)),
        }
    }
    if out.is_empty() {
        return Expr::nan();
    }
    Expr::from_kind(Kind::Piecewise(out))
}
