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
/// SymPy gives up on an unbranched argument still holding a periodic
/// argument or an arctangent; an ordinary arg(z) is a legitimate value.
fn unresolved_angle(x: &Expr) -> bool {
    matches!(x.kind(), Kind::Func(Func::Atan | Func::Atan2, _))
        || matches!(x.kind(), Kind::Func(Func::Special(SpecialFn::PeriodicArgument), _))
        || x.args().iter().any(unresolved_angle)
}
fn has_head(x: &Expr, s: SpecialFn) -> bool {
    matches!(x.kind(), Kind::Func(Func::Special(g), _) if *g == s) || x.args().iter().any(|a| has_head(a, s))
}
fn reduce(u: &Expr, period: &Expr) -> Option<Expr> {
    if period.is_constant(Constant::Infinity) { return Some(u.clone()); }
    let turns = func(Func::Ceiling, vec![sub(&div(u, period), &Expr::rational(1,2))]);
    if !has_ceiling(&turns) { Some(sub(u, &mul2(&turns, period))) } else { None }
}
fn has_ceiling(x: &Expr) -> bool {
    matches!(x.kind(), Kind::Func(Func::Ceiling, _)) || x.args().iter().any(has_ceiling)
}
fn positive(x: &Expr) -> bool { crate::assumptions::is_positive(x) == Some(true) }
/// `periodic_argument(z, period)` evaluated like the reference.
fn periodic_argument(x: &Expr, p: &Expr) -> Option<Expr> {
    if !positive(p) && !p.is_constant(Constant::Infinity) { return None; }
    let infinite = p.is_constant(Constant::Infinity);
    if let Kind::Func(Func::Special(SpecialFn::PrincipalBranch), args) = x.kind() {
        if infinite { return Some(angle(&args[0], &args[1])); }
    }
    if let Kind::Func(Func::Special(SpecialFn::PolarLift), args) = x.kind() {
        let full = infinite || crate::assumptions::is_nonnegative(&sub(p, &mul2(&Expr::int(2), &Expr::pi()))) == Some(true);
        if full { return Some(angle(&args[0], p)); }
    }
    if x.is_mul() {
        let all = x.args();
        let rest: Vec<_> = all.iter().filter(|v| !positive(v)).cloned().collect();
        if rest.len() != all.len() { return Some(angle(&mul(&rest), p)); }
    }
    let u = unwrapped(x)?;
    if unresolved_angle(&u) { return None; }
    reduce(&u, p)
}
fn unbranched(x: &Expr) -> Expr { angle(x, &Expr::oo()) }
/// `principal_branch(z, period)` evaluated like the reference.
fn principal_branch(x: &Expr, p: &Expr) -> Option<Expr> {
    let pb = |a: &Expr| func(Func::Special(SpecialFn::PrincipalBranch), vec![a.clone(), p.clone()]);
    if let Kind::Func(Func::Special(SpecialFn::PolarLift), args) = x.kind() {
        return Some(pb(&args[0]));
    }
    if p.is_constant(Constant::Infinity) { return Some(x.clone()); }
    if !positive(p) { return None; }
    let ub = unbranched(x);
    let barg = angle(x, p);
    if ub != barg && !has_head(&ub, SpecialFn::PeriodicArgument) && !has_head(&barg, SpecialFn::PeriodicArgument) {
        // Lift every non-symbol factor; plain symbols stay on the plane.
        fn strip(e: &Expr) -> Expr {
            if let Kind::Func(Func::Special(SpecialFn::PolarLift), a) = e.kind() {
                let inner = strip(&a[0]);
                return if inner.is_symbol() { inner } else { lift(&inner) };
            }
            crate::traverse::map_children(e, &mut |c| strip(c))
        }
        let pl = strip(&lift(x));
        if !has_head(&pl, SpecialFn::PolarLift) {
            let ub = unbranched(&pl);
            let mut res = if ub != barg { mul2(&exp_polar(mul2(&Expr::i(), &sub(&barg, &ub))), &pl) } else { pl };
            if !is_polar(&res) && !has_exp_polar(&res) { res = mul2(&res, &exp_polar(Expr::zero())); }
            return Some(res);
        }
    }
    let free = crate::traverse::free_symbols(x);
    let (mut c, mut m) = (Vec::new(), Vec::new());
    if free.is_empty() { c.push(x.clone()); } else {
        let factors = if x.is_mul() { x.args() } else { vec![x.clone()] };
        for y in factors {
            if crate::traverse::free_symbols(&y).is_empty() || positive(&y) { c.push(y); } else { m.push(y); }
        }
    }
    let c = mul(&c);
    let argument = angle(&c, p);
    if has_head(&argument, SpecialFn::PeriodicArgument) { return None; }
    let number = crate::traverse::free_symbols(&argument).is_empty();
    let rest = mul(&m);
    if number && (unbranched(&c) != argument || (argument.is_zero() && !m.is_empty() && !c.is_one())) {
        if argument.is_zero() { return Some(mul2(&f::abs(&c), &pb(&rest))); }
        return Some(mul2(&pb(&mul2(&exp_polar(mul2(&Expr::i(), &argument)), &rest)), &f::abs(&c)));
    }
    if number && m.is_empty() {
        let half = div(p, &Expr::int(2));
        let inside = crate::assumptions::is_negative(&sub(&f::abs(&argument), &half)) == Some(true);
        if inside || argument == half { return Some(mul2(&exp_polar(mul2(&argument, &Expr::i())), &f::abs(&c))); }
    }
    None
}
fn has_exp_polar(x: &Expr) -> bool {
    x.is_func(&Func::ExpPolar) || x.args().iter().any(has_exp_polar)
}
/// `x.has(exp_polar, polar_lift)`.
pub(crate) fn has_polar_func(x: &Expr) -> bool {
    matches!(x.kind(), Kind::Func(Func::ExpPolar | Func::Special(SpecialFn::PolarLift), _))
        || x.args().iter().any(has_polar_func)
}
pub(crate) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    let x = &a[0];
    match s {
        SpecialFn::PolarLift => {
            if is_polar(x) { return Some(x.clone()); }
            if x.is_zero() { return None; }
            let argument = f::arg(x);
            if !unresolved_angle(&argument) && !matches!(argument.kind(), Kind::Func(Func::Arg, _))
                && crate::traverse::free_symbols(x).is_empty() {
                let half = div(&Expr::pi(), &Expr::int(2));
                if argument.is_zero() || argument == Expr::pi() || argument == half || argument == neg(&half) {
                    return Some(mul2(&f::abs(x), &exp_polar(mul2(&Expr::i(), &argument))));
                }
            }
            if positive(x) {
                return Some(mul2(x, &exp_polar(Expr::zero())));
            }
            if x.is_mul() {
                let (known, rest): (Vec<_>, Vec<_>) = x.args().into_iter().partition(|v| is_polar(v) || positive(v));
                if !known.is_empty() {
                    let mut factors = known;
                    if !rest.is_empty() { factors.push(lift(&mul(&rest))); }
                    else if !factors.iter().any(is_polar) { factors.push(exp_polar(Expr::zero())); }
                    return Some(mul(&factors));
                }
            }
            None
        }
        SpecialFn::PeriodicArgument => periodic_argument(x, &a[1]),
        SpecialFn::PrincipalBranch => principal_branch(x, &a[1]),
        _ => None,
    }
}

/// Functions the reference marks `unbranched`: they act on the projection
/// to the complex plane, so polar arguments may be unpolarified inside.
fn unbranched_head(f: &Func) -> bool {
    use SpecialFn as S;
    matches!(f, Func::Re | Func::Im | Func::Abs | Func::Exp | Func::ExpPolar
        | Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth
        | Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc
        | Func::Erf | Func::Erfc | Func::Gamma)
        || matches!(f, Func::Special(S::Sech | S::Csch | S::AiryAi | S::AiryBi | S::AiryAiPrime
            | S::AiryBiPrime | S::Erfi | S::FresnelS | S::FresnelC | S::Beta))
}

/// `polarify(eq, subs=False, lift)`: numbers become polar numbers under the
/// standard argument; with `lift`, symbols and sums are lifted as well.
pub fn polarify(e: &Expr, lift_all: bool) -> Expr { polarify_rec(e, lift_all, false) }
fn polarify_rec(e: &Expr, lift_all: bool, pause: bool) -> Expr {
    if is_polar(e) { return e.clone(); }
    if !pause && is_number_expr(e) { return lift(e); }
    match e.kind() {
        Kind::Sym(_) if !pause && lift_all => lift(e),
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => e.clone(),
        Kind::Add(..) => {
            let r = add(&e.args().iter().map(|a| polarify_rec(a, lift_all, true)).collect::<Vec<_>>());
            if lift_all { lift(&r) } else { r }
        }
        // the `hyper` and `meijerg` constructors unpolarify the parameters
        // they are rebuilt with
        Kind::Func(Func::Special(SpecialFn::Hyper | SpecialFn::MeijerG), args) => {
            let last = args.len().saturating_sub(1);
            let new: Vec<Expr> = args
                .iter()
                .enumerate()
                .map(|(i, a)| if i == last { polarify_rec(a, lift_all, false) } else { unpolarify_params(&polarify_rec(a, lift_all, false)) })
                .collect();
            crate::simptools::rebuild(e, &new)
        }
        Kind::Func(..) => crate::simptools::map_args(e, &mut |a| polarify_rec(a, lift_all, false)),
        _ => crate::simptools::map_args(e, &mut |a| polarify_rec(a, lift_all, pause)),
    }
}
/// The parameter tuples (possibly nested) of `hyper`/`meijerg` with each
/// parameter unpolarified.
fn unpolarify_params(t: &Expr) -> Expr {
    match super::tuple_items(t) {
        Some(items) => {
            let new: Vec<Expr> = items.iter().map(unpolarify_params).collect();
            crate::simptools::rebuild(t, &new)
        }
        None => unpolarify(t, false),
    }
}
fn is_number_expr(e: &Expr) -> bool {
    // A Tuple is not an Expr, so it is never `is_number`.
    !matches!(e.kind(), Kind::Bool(_) | Kind::Rel(..)) && !crate::logic::is_boolean(e) && super::tuple_items(e).is_none()
        && crate::traverse::free_symbols(e).is_empty() && crate::functions::is_numeric_value(e)
}

/// `unpolarify(eq, exponents_only)`: project polar numbers to the plane
/// wherever the value is unchanged by doing so.
pub fn unpolarify(e: &Expr, exponents_only: bool) -> Expr {
    let mut eq = e.clone();
    loop {
        let res = unpolarify_rec(&eq, exponents_only, exponents_only);
        if res == eq { break; }
        eq = res;
    }
    // exp_polar(0) -> 1 and polar_lift(0) -> 0 are always correct.
    fn finish(e: &Expr) -> Expr {
        match e.kind() {
            Kind::Func(Func::ExpPolar, a) if a[0].is_zero() => Expr::one(),
            Kind::Func(Func::Special(SpecialFn::PolarLift), a) if a[0].is_zero() => Expr::zero(),
            _ => crate::simptools::map_args(e, &mut |c| finish(c)),
        }
    }
    finish(&eq)
}
fn unpolarify_rec(e: &Expr, exponents_only: bool, pause: bool) -> Expr {
    if e.is_atom() { return e.clone(); }
    if !pause {
        match e.kind() {
            Kind::Func(Func::ExpPolar, a) => return f::exp(&unpolarify_rec(&a[0], exponents_only, false)),
            Kind::Func(Func::Special(SpecialFn::PrincipalBranch), a) if a[1] == mul2(&Expr::int(2), &Expr::pi()) =>
                return unpolarify_rec(&a[0], exponents_only, false),
            Kind::Add(..) | Kind::Mul(..) | Kind::Bool(_) => {
                return crate::simptools::map_args(e, &mut |c| unpolarify_rec(c, exponents_only, false));
            }
            Kind::Rel(op, l, r) => {
                let equality = matches!(op, crate::expr::RelOp::Eq | crate::expr::RelOp::Ne);
                if !equality || l.is_zero() || r.is_zero() {
                    return crate::relational::rel(*op, &unpolarify_rec(l, exponents_only, false), &unpolarify_rec(r, exponents_only, false));
                }
            }
            Kind::Func(Func::Logic(_), _) => {
                return crate::simptools::map_args(e, &mut |c| unpolarify_rec(c, exponents_only, false));
            }
            Kind::Func(Func::Special(SpecialFn::PolarLift), a) => return unpolarify_rec(&a[0], exponents_only, false),
            _ => {}
        }
    }
    match e.kind() {
        Kind::Pow(b, x) => {
            let expo = unpolarify_rec(x, exponents_only, false);
            let integer = crate::assumptions::is_integer(&expo) == Some(true);
            let base = unpolarify_rec(b, exponents_only, !(integer && !pause));
            crate::ops::pow(&base, &expo)
        }
        Kind::Func(h, _) if unbranched_head(h) =>
            crate::simptools::map_args(e, &mut |c| unpolarify_rec(c, exponents_only, exponents_only)),
        _ => crate::simptools::map_args(e, &mut |c| unpolarify_rec(c, exponents_only, true)),
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
