//! Polar lifting and argument reduction on the logarithm surface.
use crate::{Expr, Func, Kind};
use crate::expr::Constant;
use crate::functions::{self as f, func};
use crate::ops::{add, div, mul, mul2, neg, sub};
use super::SpecialFn;

pub(crate) fn is_polar(x: &Expr) -> bool {
    match x.kind() {
        Kind::Func(Func::ExpPolar, _) => true,
        Kind::Func(Func::Special(SpecialFn::PolarLift | SpecialFn::PrincipalBranch), _) => true,
        Kind::Sym(s) => s.facts().get("polar") == Some(true),
        Kind::Mul(c, factors) => c.is_positive() && factors.iter().any(|(b, _)| is_polar(b))
            && factors.iter().all(|(b, _)| is_polar(b) || crate::assumptions::is_positive(b) == Some(true)),
        Kind::Pow(b, _) => is_polar(b),
        _ => false,
    }
}
fn exp_polar(x: Expr) -> Expr { func(Func::ExpPolar, vec![x]) }
fn lift(x: &Expr) -> Expr { func(Func::Special(SpecialFn::PolarLift), vec![x.clone()]) }
fn angle(x: &Expr, p: &Expr) -> Expr { func(Func::Special(SpecialFn::PeriodicArgument), vec![x.clone(), p.clone()]) }

/// A resolved unwrapped angle; an ordinary arg(z) remains legitimate here.
fn unwrapped(x: &Expr) -> Option<Expr> {
    match x.kind() {
        Kind::Func(Func::ExpPolar, a) => Some(f::im(&a[0])),
        Kind::Func(Func::Special(SpecialFn::PolarLift), a) => Some(f::arg(&a[0])),
        Kind::Func(Func::Special(SpecialFn::PrincipalBranch), a) => Some(angle(&a[0], &a[1])),
        Kind::Mul(_, _) => {
            let angles: Option<Vec<_>> = x.args().iter().map(unwrapped).collect();
            Some(add(&angles?))
        }
        Kind::Pow(b, e) if is_polar(b) => Some(add(&[
            mul2(&f::re(e), &unwrapped(b)?), mul2(&f::im(e), &f::log(&f::abs(b))),
        ])),
        _ if !is_polar(x) => Some(f::arg(x)),
        _ => None,
    }
}
fn unresolved_angle(x: &Expr) -> bool {
    matches!(x.kind(), Kind::Func(Func::Arg | Func::Atan | Func::Atan2, _))
        || matches!(x.kind(), Kind::Func(Func::Special(SpecialFn::PeriodicArgument), _))
        || x.args().iter().any(unresolved_angle)
}
fn reduce(u: &Expr, period: &Expr) -> Option<Expr> {
    if period.is_constant(Constant::Infinity) { return Some(u.clone()); }
    let turns = func(Func::Ceiling, vec![sub(&div(u, period), &Expr::rational(1,2))]);
    if turns.is_integer_number() { Some(sub(u, &mul2(&turns, period))) } else { None }
}
pub(crate) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    let x = &a[0];
    match s {
        SpecialFn::PolarLift => {
            if is_polar(x) { return Some(x.clone()); }
            if x.is_zero() { return None; }
            let argument = f::arg(x);
            if !unresolved_angle(&argument) && crate::traverse::free_symbols(x).is_empty() {
                let half = div(&Expr::pi(), &Expr::int(2));
                if argument.is_zero() || argument == Expr::pi() || argument == half || argument == neg(&half) {
                    return Some(mul2(&f::abs(x), &exp_polar(mul2(&Expr::i(), &argument))));
                }
            }
            if crate::assumptions::is_positive(x) == Some(true) {
                return Some(mul2(x, &exp_polar(Expr::zero())));
            }
            if x.is_mul() {
                let (known, rest): (Vec<_>, Vec<_>) = x.args().into_iter().partition(|v|
                    is_polar(v) || crate::assumptions::is_positive(v) == Some(true));
                if !known.is_empty() {
                    let mut factors = known;
                    if !rest.is_empty() { factors.push(lift(&mul(&rest))); }
                    else if !factors.iter().any(is_polar) { factors.push(exp_polar(Expr::zero())); }
                    return Some(mul(&factors));
                }
            }
            None
        }
        SpecialFn::PeriodicArgument => {
            let p = &a[1];
            if crate::assumptions::is_positive(p) != Some(true) && !p.is_constant(Constant::Infinity) { return None; }
            if let Kind::Func(Func::Special(SpecialFn::PrincipalBranch), args) = x.kind() {
                if p.is_constant(Constant::Infinity) { return Some(angle(&args[0], &args[1])); }
            }
            if x.is_mul() {
                let rest: Vec<_> = x.args().into_iter().filter(|v| crate::assumptions::is_positive(v) != Some(true)).collect();
                let reduced = mul(&rest);
                if reduced != *x { return Some(angle(&reduced, p)); }
            }
            let u = unwrapped(x)?;
            if unresolved_angle(&u) { return None; }
            reduce(&u, p)
        }
        SpecialFn::PrincipalBranch => {
            let p = &a[1];
            if p.is_constant(Constant::Infinity) { return Some(x.clone()); }
            if crate::assumptions::is_positive(p) != Some(true) { return None; }
            if let Kind::Func(Func::Special(SpecialFn::PolarLift), args) = x.kind() {
                return Some(func(Func::Special(s), vec![args[0].clone(), p.clone()]));
            }
            if x.is_zero() { return None; }
            let u = unwrapped(x)?;
            if unresolved_angle(&u) { return None; }
            let reduced = reduce(&u, p)?;
            // A magnitude on the complex plane and an explicit polar phase
            // avoid dropping the winding or introducing principal powers.
            Some(mul2(&f::abs(x), &exp_polar(mul2(&Expr::i(), &reduced))))
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn floats_keep_the_polar_head() {
        let ep = exp_polar(Expr::float(0.25));
        assert!(ep.is_func(&Func::ExpPolar));
        assert_eq!(crate::evalf::eval_f64(&ep), Some(0.25_f64.exp()));
    }
    #[test]
    fn winding_and_boundary() {
        let ep = exp_polar(mul(&[Expr::int(5), Expr::i(), Expr::pi()]));
        assert_eq!(angle(&ep, &Expr::oo()), mul2(&Expr::int(5), &Expr::pi()));
        assert_eq!(angle(&ep, &mul2(&Expr::int(2), &Expr::pi())), Expr::pi());
        let minus = exp_polar(neg(&mul2(&Expr::i(), &Expr::pi())));
        assert_eq!(angle(&minus, &mul2(&Expr::int(2), &Expr::pi())), Expr::pi());
    }
    #[test]
    fn lift_and_log_preserve_sheets() {
        let value = lift(&Expr::int(-4));
        assert!(is_polar(&value));
        assert_eq!(f::abs(&value), Expr::int(4));
        let z = Expr::symbol("z");
        let ep = exp_polar(z.clone());
        assert_eq!(f::log(&ep), z.clone());
        assert_eq!(crate::ops::pow(&ep, &Expr::rational(1,2)), exp_polar(div(&z, &Expr::int(2))));
    }
}
