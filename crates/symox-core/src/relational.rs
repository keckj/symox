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
    let equality = matches!(op, RelOp::Eq | RelOp::Ne);
    // SymPy's is_ge decides an ordering of symbolic sides only when both
    // are extended real (`y < y` stays unevaluated for a complex `y`).
    let extended_real = |e: &Expr| crate::assumptions::fact(e, "extended_real") == Some(true);
    if a == b && (equality || extended_real(a)) {
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
    // Otherwise `_n2` compares two comparable numbers numerically.
    if !equality && !(extended_real(a) && extended_real(b)) && !(is_numeric_value(a) && is_numeric_value(b)) {
        return None;
    }
    let d = sub(a, b);
    // is_eq: `z is False and dif.is_commutative` (issue 10728); a
    // noncommutative difference is nonzero by the assumption rules, and its
    // ratio test is undecided too, so Eq(A*B, B*A) stays unevaluated.
    if matches!(op, RelOp::Eq | RelOp::Ne) && d.has_noncommutative_host()
        && crate::decompose::commutativity(&d) != Some(true) {
        return None;
    }
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
        // SymPy's is_gt/is_ge and the zero test of is_eq read the extended
        // sign of the difference: Abs(x) + pi is extended positive even when
        // its finiteness is unknown.
        let ext = |name: &str| asm::fact(&d, name) == Some(true);
        let r = match op {
            RelOp::Gt if ext("extended_positive") => Some(true),
            RelOp::Gt if ext("extended_nonpositive") => Some(false),
            RelOp::Ge if ext("extended_nonnegative") => Some(true),
            RelOp::Ge if ext("extended_negative") => Some(false),
            RelOp::Lt if ext("extended_negative") => Some(true),
            RelOp::Lt if ext("extended_nonnegative") => Some(false),
            RelOp::Le if ext("extended_nonpositive") => Some(true),
            RelOp::Le if ext("extended_positive") => Some(false),
            RelOp::Eq | RelOp::Ne if ext("extended_positive") || ext("extended_negative") => Some(op == RelOp::Ne),
            _ => None,
        };
        if r.is_some() {
            return r;
        }
        // SymPy's is_eq then reads the ratio: a finite nonzero numerator
        // over a finite denominator is not zero (-1/s**2 - 1 for real s).
        if matches!(op, RelOp::Eq | RelOp::Ne) {
            let (n, den) = crate::simplify::numer_denom(&d);
            let infinite = |e: &Expr| asm::is_finite(e).map(|f| !f);
            let rv = if asm::is_zero(&n) == Some(true) {
                asm::is_zero(&den).map(|z| !z)
            } else if asm::is_finite(&n) == Some(true) {
                if infinite(&den) == Some(true) {
                    Some(true)
                } else if asm::is_zero(&n) == Some(false) {
                    infinite(&den).or_else(|| {
                        // The condition that makes the denominator infinite
                        // decides inequality when it does not make the
                        // original sides equal.
                        let (l, r) = clear_coefficients(&den, &Expr::oo())?;
                        let pair = [(l, r)];
                        let (la, lb) = (crate::subs::subs(a, &pair), crate::subs::subs(b, &pair));
                        if &la == a && &lb == b {
                            return None;
                        }
                        match decide(RelOp::Eq, &la, &lb) {
                            Some(false) => Some(false),
                            _ => None,
                        }
                    })
                } else {
                    None
                }
            } else {
                None
            };
            // rv is SymPy's "difference is zero" verdict.
            if let Some(z) = rv {
                return Some((op == RelOp::Eq) == z);
            }
        }
    }
    // A numeric difference of symbolic sides orders them only when both are
    // extended real (SymPy's is_ge): s < s + 1 stays unevaluated for a
    // complex s.
    if !matches!(op, RelOp::Eq | RelOp::Ne) && !(a.is_number() && b.is_number()) {
        let real = |e: &Expr| crate::assumptions::fact(e, "extended_real") == Some(true);
        if !(real(a) && real(b)) {
            return None;
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

/// SymPy's ``clear_coefficients(expr, rhs)`` for an expression with free
/// symbols: strip rational content and additive constants into ``rhs``.
fn clear_coefficients(expr: &Expr, rhs: &Expr) -> Option<(Expr, Expr)> {
    use crate::ops::{div, neg};
    if crate::traverse::free_symbols(expr).is_empty() {
        return None;
    }
    let (mut expr, mut rhs) = (expr.clone(), rhs.clone());
    let mut was: Option<Expr> = None;
    while !expr.is_zero() && was.as_ref() != Some(&expr) {
        was = Some(expr.clone());
        let (m, e) = crate::decompose::as_content_primitive(&expr, false, true).ok()?;
        rhs = div(&rhs, &Expr::num(m));
        let (c, e) = crate::decompose::as_coeff_add_head(&e, true);
        rhs = sub(&rhs, &c);
        expr = e;
    }
    let expr = crate::simptools::misc::signsimp(&expr, false);
    if crate::decompose::could_extract_minus_sign(&expr) {
        return Some((neg(&expr), neg(&rhs)));
    }
    Some((expr, rhs))
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

/// `Piecewise((e1, c1), (e2, c2), ...)`: `Piecewise.__new__` followed by
/// `Piecewise.eval` (`_piecewise_collapse_arguments`).
pub fn piecewise(pieces: Vec<(Expr, Expr)>) -> Expr {
    piecewise_depth(pieces, 0)
}

fn piecewise_depth(pieces: Vec<(Expr, Expr)>, depth: usize) -> Expr {
    let mut out: Vec<(Expr, Expr)> = Vec::with_capacity(pieces.len());
    for (e, c) in pieces {
        // `ExprCondPair`: a condition containing a Piecewise is folded and
        // the folded Piecewise of booleans becomes an ITE
        let c = if !matches!(c.kind(), Kind::Piecewise(_)) && condition_has_piecewise(&c, &mut rustc_hash::FxHashSet::default()) {
            crate::integrate::piecewise::piecewise_fold(&c, Some(true))
        } else {
            c
        };
        let c = crate::logic::boolean_piecewise_to_ite(&c).unwrap_or(c);
        match c.kind() {
            Kind::Bool(false) => continue,
            Kind::Bool(true) => {
                out.push((e, c));
                break;
            }
            _ => out.push((e, c)),
        }
    }
    if out.is_empty() {
        return Expr::nan();
    }
    if out.len() == 1 && matches!(out[0].1.kind(), Kind::Bool(true)) {
        return out.pop().unwrap().0;
    }
    let collapsed = collapse_arguments(&out);
    if collapsed.is_empty() {
        return Expr::nan();
    }
    if collapsed != out && depth < 16 {
        return piecewise_depth(collapsed, depth + 1);
    }
    Expr::from_kind(Kind::Piecewise(out))
}

fn condition_has_piecewise(e: &Expr, seen: &mut rustc_hash::FxHashSet<usize>) -> bool {
    if !seen.insert(e.dag_id()) {
        return false;
    }
    matches!(e.kind(), Kind::Piecewise(_)) || e.args().iter().any(|a| condition_has_piecewise(a, seen))
}

/// `cond.replace(lambda _: _.is_Relational, _canonical_coeff)`.
fn canonical_conditions(c: &Expr) -> Expr {
    use crate::expr::Func;
    match c.kind() {
        Kind::Rel(..) => crate::summation::canonical_coeff(c),
        Kind::Func(Func::Logic(op), args) => {
            let new: Vec<Expr> = args.iter().map(canonical_conditions).collect();
            if new == *args {
                c.clone()
            } else {
                crate::logic::eval(*op, new)
            }
        }
        _ => c.clone(),
    }
}

/// SymPy's `_piecewise_collapse_arguments`.
fn collapse_arguments(args: &[(Expr, Expr)]) -> Vec<(Expr, Expr)> {
    use crate::expr::Func;
    use crate::logic::LogicOp;
    let is_true = |c: &Expr| matches!(c.kind(), Kind::Bool(true));
    let mut newargs: Vec<(Expr, Expr)> = Vec::with_capacity(args.len());
    let mut current: rustc_hash::FxHashSet<Expr> = rustc_hash::FxHashSet::default();
    for (expr, cond) in args {
        let mut cond = canonical_conditions(cond);
        let mut expr = expr.clone();
        if let Kind::Piecewise(inner) = expr.kind() {
            let inner = inner.clone();
            let mut unmatching: Vec<(Expr, Expr)> = Vec::new();
            for (e, c) in inner {
                if current.contains(&c) {
                    continue;
                }
                if c == cond {
                    if !is_true(&c) {
                        expr = if unmatching.is_empty() {
                            e
                        } else {
                            unmatching.push((e, c));
                            piecewise(unmatching)
                        };
                    }
                    break;
                }
                unmatching.push((e, c));
            }
        }
        let and_args: Vec<Expr> = match cond.kind() {
            Kind::Func(Func::Logic(LogicOp::And), a) => a.clone(),
            _ => Vec::new(),
        };
        if current.contains(&cond) || and_args.iter().any(|a| current.contains(a)) {
            continue;
        }
        if !and_args.is_empty() {
            let mut nonredundant = Vec::with_capacity(and_args.len());
            let mut falsified = false;
            for c in and_args {
                if let Kind::Rel(op, a, b) = c.kind() {
                    if current.contains(&crate::logic::canonical(&crate::logic::negated_rel(&c))) {
                        continue;
                    }
                    let weak = match op {
                        RelOp::Lt => Some(RelOp::Le),
                        RelOp::Gt => Some(RelOp::Ge),
                        _ => None,
                    };
                    if let Some(w) = weak {
                        if current.contains(&rel(w, a, b)) {
                            falsified = true;
                            break;
                        }
                    }
                }
                nonredundant.push(c);
            }
            cond = if falsified { Expr::boolean(false) } else { crate::logic::and(nonredundant) };
        } else if matches!(cond.kind(), Kind::Rel(..))
            && current.contains(&crate::logic::canonical(&crate::logic::negated_rel(&cond)))
        {
            cond = Expr::boolean(true);
        }
        current.insert(cond.clone());
        if let Some(last) = newargs.last_mut() {
            if last.0 == expr {
                let mut orcond = crate::logic::or2(&cond, &last.1);
                if matches!(orcond.kind(), Kind::Func(Func::Logic(LogicOp::And | LogicOp::Or), _)) {
                    orcond = crate::logic::distribute_and_over_or(&orcond);
                }
                last.1 = orcond;
                continue;
            } else if last.1 == cond {
                continue;
            }
        }
        newargs.push((expr, cond));
    }
    newargs
}
