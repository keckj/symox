//! `AccumBounds` values of `series`/`nseries` for `frac`, `floor` and
//! `ceiling` of an argument that is infinite at the expansion point.
//!
//! The reference maps the point to `0+` and, for an infinite
//! `arg0 = arg.subs(x, 0)`, `frac._eval_nseries` returns
//! `AccumBounds(0, 1) + O(x**n)` (`O(1)` for `n <= 0`), and
//! `floor._eval_nseries` (`ceiling`) returns the expansion `s` of the
//! argument plus `AccumBounds(-1, 0)` (`AccumBounds(0, 1)`), `O(1)` in
//! place of the bounds for `n <= 0`. The bounds cannot be built here (a
//! host node); the value is returned as `(bounds, rest)` and the Python
//! front end adds `AccumBounds(*bounds)` to `rest`, whose host arithmetic
//! is the reference one (`AccumBounds(-1, 0) + 1/p` with a positive `p` is
//! `AccumBounds(-1 + 1/p, 1/p)`).
//!
//! Admitted natively:
//! * sums `c0 + sum(c_i*frac(a_i)) + r` with real numbers `c_i`, each
//!   `a_i.subs(x, 0)` infinite after the point is moved to 0, and `r`
//!   free of step functions with an infinite argument and expanded by the
//!   native series (`frac(1/x) + x` is `AccumBounds(0, 1) + x + O(x**6)`);
//! * `floor(a)`/`ceiling(a)` alone, `a` a Laurent polynomial in the moved
//!   variable, for `series` of a variable known positive or negative (any
//!   other variable is the `ValueError` of `accum_floor`) and for
//!   `nseries` at `0+`.
//!
//! * a sum holding `c*frac(a)` whose `a.subs(x, 0)` is `nan`
//!   (`frac(x**2 + x)` at `oo`, `frac(exp(1/x))` at 0): `nan`, the value of
//!   the reference, once the rest of the sum has a native expansion.
//!
//! Everything else declines (the previous `PoleError`): products
//! of a bounded node by a variable, and `floor` of an argument whose
//! expansion is not exact (the reference keeps an unevaluated sum).

use super::engine::depends;
use crate::assumptions::{is_negative, is_positive};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{add, neg, pow, sub};
use crate::special::SpecialFn;
use crate::subs::{subs1, xreplace};

/// `(bounds, rest)`: the value is `AccumBounds(lo, hi) + rest`, or `rest`
/// alone when `bounds` is `None`.
pub type AccumValue = (Option<(Expr, Expr)>, Expr);

fn infinite(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity) || e.is_constant(Constant::ComplexInfinity)
}

/// The step functions of the module notes.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Step {
    Floor,
    Ceiling,
    Frac,
}

fn step_of(e: &Expr) -> Option<(Step, Expr)> {
    let (s, a) = match e.kind() {
        Kind::Func(Func::Floor, a) => (Step::Floor, a),
        Kind::Func(Func::Ceiling, a) => (Step::Ceiling, a),
        Kind::Func(Func::Special(SpecialFn::Frac), a) => (Step::Frac, a),
        _ => return None,
    };
    (a.len() == 1).then(|| (s, a[0].clone()))
}

fn is_step(e: &Expr) -> bool {
    step_of(e).is_some()
}

fn has_step(e: &Expr, u: &Expr) -> bool {
    depends(e, u) && (is_step(e) || e.args().iter().any(|a| has_step(a, u)))
}

/// `F(a)` of a step function whose argument is infinite at `u = 0`.
fn infinite_step(e: &Expr, u: &Expr) -> Option<(Step, Expr)> {
    let (f, a) = step_of(e)?;
    if !depends(&a, u) {
        return None;
    }
    infinite(&subs1(&a, u, &Expr::zero())).then_some((f, a))
}

/// `(c, F(a))` for `c*F(a)` with a real number `c`.
fn scaled_frac(t: &Expr, u: &Expr) -> Option<Expr> {
    let (c, f) = match t.kind() {
        Kind::Mul(c, facs) if facs.len() == 1 && facs[0].1.is_one() => (Expr::num(c.clone()), facs[0].0.clone()),
        _ => (Expr::one(), t.clone()),
    };
    match infinite_step(&f, u)? {
        (Step::Frac, _) if c.as_number().is_some() => Some(c),
        _ => None,
    }
}

/// `c*frac(a)` with a real number `c` and `a.subs(x, 0)` equal to `nan`
/// (`frac(x**2 + x)` at `oo`, `frac(exp(1/x))` at 0): `frac._eval_nseries`
/// then expands `arg - nan`, and the value is `nan`.
fn nan_frac(t: &Expr, u: &Expr) -> bool {
    let (c, f) = match t.kind() {
        Kind::Mul(c, facs) if facs.len() == 1 && facs[0].1.is_one() => (Some(c.clone()), facs[0].0.clone()),
        _ => (None, t.clone()),
    };
    if c.is_some_and(|c| c.is_zero()) {
        return false;
    }
    match step_of(&f) {
        Some((Step::Frac, a)) => depends(&a, u) && subs1(&a, u, &Expr::zero()).is_constant(Constant::NaN),
        _ => false,
    }
}

/// The order term `O((x - x0)**n)` of an expansion at `x0` (`O(1)` for
/// `n <= 0`).
fn order_at(x: &Expr, x0: &Expr, n: i64) -> Expr {
    let k = n.max(0);
    let base = if infinite(x0) { pow(x, &Expr::minus_one()) } else { sub(x, x0) };
    crate::series::order::order(&pow(&base, &Expr::int(k)), &[(x.clone(), x0.clone())])
}

/// The value of `e.series(x, x0, n, dir)` (`nseries` at `0+` when
/// `nseries`) described in the module notes; `None` when not decided.
pub fn accum_series_value(e: &Expr, x: &Expr, x0: &Expr, n: i64, minus: bool, nseries: bool) -> Option<AccumValue> {
    if !x.is_symbol() || depends(x0, x) || x0.is_constant(Constant::ComplexInfinity) || x0.is_constant(Constant::NaN) {
        return None;
    }
    let u = Expr::dummy("u", crate::symbol::Assumptions::POSITIVE);
    let (v, back) = if x0.is_constant(Constant::Infinity) {
        (pow(&u, &Expr::minus_one()), pow(x, &Expr::minus_one()))
    } else if x0.is_constant(Constant::NegInfinity) {
        (neg(&pow(&u, &Expr::minus_one())), neg(&pow(x, &Expr::minus_one())))
    } else if minus {
        (sub(x0, &u), sub(x0, x))
    } else {
        (add(&[x0.clone(), u.clone()]), sub(x, x0))
    };
    let e0 = subs1(e, x, &v);
    if !has_step(&e0, &u) {
        return None;
    }
    // floor / ceiling alone
    if let Some((f @ (Step::Floor | Step::Ceiling), a)) = infinite_step(&e0, &u) {
        if !nseries && is_positive(x) != Some(true) && is_negative(x) != Some(true) {
            return None;
        }
        // an exact (Laurent polynomial) expansion of the argument
        let s = crate::series::series(&a, &u, &Expr::zero(), n.max(1), false).ok()?;
        if s.may_have_order() || crate::series::engine::degree_in(&crate::expand::expand(&crate::ops::mul2(&s, &pow(&u, &Expr::int(64)))), &u).is_none() {
            return None;
        }
        let back_map = |e: &Expr| xreplace(e, &std::iter::once((u.clone(), back.clone())).collect());
        let s = back_map(&s);
        if n <= 0 {
            return Some((None, add(&[s, order_at(x, x0, n)])));
        }
        let b = if f == Step::Floor { (Expr::minus_one(), Expr::zero()) } else { (Expr::zero(), Expr::one()) };
        return Some((Some(b), s));
    }
    // c0 + sum(c_i*frac(a_i)) + r
    let terms = if e0.is_add() { e0.args() } else { vec![e0.clone()] };
    let (mut lo, mut hi) = (Expr::zero(), Expr::zero());
    let mut rest = Vec::new();
    let mut found = false;
    let mut nan = false;
    for t in terms {
        if !has_step(&t, &u) {
            rest.push(t);
            continue;
        }
        if nan_frac(&t, &u) {
            nan = true;
            continue;
        }
        let c = scaled_frac(&t, &u)?;
        found = true;
        if c.as_number().is_some_and(|k| k.is_negative()) {
            lo = add(&[lo, c]);
        } else {
            hi = add(&[hi, c]);
        }
    }
    if !found && !nan {
        return None;
    }
    // the rest is expanded in the original variable
    let back_map = |e: &Expr| xreplace(e, &std::iter::once((u.clone(), back.clone())).collect());
    let r = back_map(&add(&rest));
    let rs = crate::series::series(&r, x, x0, n, minus).ok()?;
    if nan {
        return Some((None, Expr::nan()));
    }
    let rs = add(&[rs, order_at(x, x0, n)]);
    if n <= 0 {
        return Some((None, rs));
    }
    Some((Some((lo, hi)), rs))
}
