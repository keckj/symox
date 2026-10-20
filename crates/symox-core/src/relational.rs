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
    use crate::expr::Constant;
    if a.is_constant(Constant::NaN) || b.is_constant(Constant::NaN) {
        return match op { RelOp::Eq=>Some(false),RelOp::Ne=>Some(true),_=>None };
    }
    if a == b {
        return Some(matches!(op, RelOp::Eq | RelOp::Le | RelOp::Ge));
    }
    if matches!(op, RelOp::Eq | RelOp::Ne) {
        // Boolean atoms are distinct from arithmetic 0/1. Symbols retain the
        // reference's undecided Boolean/arithmetic role; numerical coercion
        // must not turn a mixed-domain equality into a subtraction of booleans.
        let boolean = |e: &Expr| matches!(e.kind(), Kind::Bool(_) | Kind::Rel(..)
            | Kind::Func(crate::Func::Logic(_) | crate::Func::AppliedPredicate
                | crate::Func::Set(crate::sets::SetOp::Contains), _));
        if matches!((a.kind(), b.kind()), (Kind::Bool(_), Kind::Bool(_)))
            || (!a.is_symbol() && !b.is_symbol() && boolean(a) != boolean(b)) {
            return Some(op == RelOp::Ne);
        }
    }
    if let (Some(a), Some(b)) = (a.as_number(), b.as_number()) {
        return Some(cmp_holds(op, a.cmp_value(b)?));
    }
    let infinity=|value:&Expr| {
        if value.is_constant(Constant::Infinity) { Some(1) }
        else if value.is_constant(Constant::NegInfinity) { Some(-1) }
        else { None }
    };
    match (infinity(a),infinity(b)) {
        (Some(x),Some(y))=>return Some(cmp_holds(op,x.cmp(&y))),
        (Some(sign),None)=>return infinite_decision(op,sign,b),
        (None,Some(sign))=>{
            let swapped=match op {RelOp::Lt=>RelOp::Gt,RelOp::Le=>RelOp::Ge,RelOp::Gt=>RelOp::Lt,RelOp::Ge=>RelOp::Le,other=>other};
            return infinite_decision(swapped,sign,a);
        }
        (None,None)=>{},
    }
    let d = sub(a, b);
    if !d.is_number() {
        use crate::assumptions as asm;
        // A sign fact that fails, or the complementary one, decides an
        // ordering only for an extended-real difference: a non-real one is
        // neither nonnegative nor negative (SymPy's ``is_ge`` requires
        // extended-real sides).
        let real = || asm::fact(&d, "extended_real") == Some(true);
        let sign = |holds: Option<bool>, opposite: &dyn Fn() -> Option<bool>| match holds {
            Some(true) => Some(true),
            _ if !real() => None,
            Some(false) => Some(false),
            None => opposite().map(|b| !b),
        };
        let r = match op {
            RelOp::Gt => sign(asm::is_positive(&d), &|| asm::is_nonpositive(&d)),
            RelOp::Ge => sign(asm::is_nonnegative(&d), &|| asm::is_negative(&d)),
            RelOp::Lt => sign(asm::is_negative(&d), &|| asm::is_nonnegative(&d)),
            RelOp::Le => sign(asm::is_nonpositive(&d), &|| asm::is_positive(&d)),
            RelOp::Eq => asm::is_zero(&d),
            RelOp::Ne => asm::is_zero(&d).map(|b| !b),
        };
        if r.is_some() {
            return r;
        }
    }
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

fn infinite_decision(op:RelOp, sign:i32, other:&Expr)->Option<bool> {
    use crate::assumptions as assumptions;
    let real=assumptions::is_real(other)==Some(true);
    let extended=assumptions::fact(other,"extended_real");
    let finite=real || assumptions::is_finite(other)==Some(true);
    if matches!(op,RelOp::Eq|RelOp::Ne) {
        return (finite || extended==Some(false)).then_some(op==RelOp::Ne);
    }
    if real || extended==Some(true) {
        if finite {return Some(cmp_holds(op,if sign>0{Ordering::Greater}else{Ordering::Less}));}
        // An extended-real symbol can itself be an infinity.
        return match (sign>0,op) {
            (true,RelOp::Ge)|(false,RelOp::Le)=>Some(true),
            (true,RelOp::Lt)|(false,RelOp::Gt)=>Some(false),
            _=>None,
        };
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

#[cfg(test)]
mod boolean_equality_tests {
    use super::*;

    #[test]
    fn boolean_atoms_do_not_equal_arithmetic_integers() {
        for boolean in [false, true] {
            for number in [-1, 0, 1, 2] {
                let boolean = Expr::boolean(boolean);
                let number = Expr::int(number);
                assert_eq!(eq(&boolean, &number), Expr::boolean(false));
                assert_eq!(eq(&number, &boolean), Expr::boolean(false));
                assert_eq!(ne(&boolean, &number), Expr::boolean(true));
            }
        }
        assert_eq!(eq(&Expr::boolean(false), &Expr::boolean(true)), Expr::boolean(false));
    }

    #[test]
    fn symbols_can_still_represent_boolean_variables() {
        let x = Expr::symbol("x");
        assert!(matches!(eq(&x, &Expr::boolean(true)).kind(), Kind::Rel(RelOp::Eq, ..)));
    }
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
        let c = crate::logic::boolean_piecewise_to_ite(&c).unwrap_or(c);
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
