//! Native second-order matches of `classify_ode`.
//!
//! Ports of the SymPy 1.14 matching steps of `2nd_hypergeometric`
//! (`match_2nd_hypergeometric`), `2nd_linear_airy`, `2nd_linear_bessel` (its
//! rejection of a constant leading coefficient) and of the order-two power
//! series hints (`2nd_power_series_ordinary`, `2nd_power_series_regular`).
//!
//! The wild matching uses the same matcher as the public `Expr.match`
//! (`crate::pattern::pattern_match`), the collection the same kernel as the
//! public `collect`; `cancel` is the public callable supplied by the binding.
//! `Decline` leaves the decision to the reference matcher.

use crate::diff::derivative_unevaluated;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{add, div, mul2, pow, sub};
use crate::pattern::pattern_match;
use crate::poly::collect::{self, Collection, Options};
use crate::traverse::{free_symbols, has_structural as has, preorder};

use super::classify::Outcome;
use super::linear::{linear_coefficients, unmodelled};

/// `cancel`; `None` declines.
pub type Cancel<'a> = &'a mut dyn FnMut(&Expr) -> Option<Expr>;

/// Special functions, whose substitution can raise in the public `subs`.
pub(super) fn has_special(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Func(Func::Special(_), _)) {
            found = true;
        }
    });
    found
}

/// The reference `e != 0` is false only on an exact zero; `None` for a
/// float, whose comparison the port does not model.
pub(super) fn nonzero(e: &Expr) -> Option<bool> {
    if has_float(e) {
        return None;
    }
    Some(!e.is_zero())
}

pub(super) fn has_float(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Num(n) if n.is_float()) {
            found = true;
        }
    });
    found
}

pub(super) fn d(fx: &Expr, x: &Expr, n: u32) -> Expr {
    derivative_unevaluated(fx, vec![(x.clone(), n)])
}

pub(super) fn is_polynomial(e: &Expr) -> bool {
    crate::inspection::classify(e, &free_symbols(e), crate::inspection::PolynomialKind::Polynomial) == Some(true)
}

/// `collect(eq, [f'', f', f]).match(a3*f'' + b3*f' + c3*f)` followed by the
/// reference's polynomial check: `[a3, b3, c3]`. A coefficient that is not
/// polynomial (the reference then matches the expanded numerator again)
/// declines.
pub fn deq_match(eq: &Expr, fx: &Expr, x: &Expr) -> Outcome<[Expr; 3]> {
    if unmodelled(eq) || eq.has_matrix() {
        return Outcome::Decline;
    }
    let (d2, d1) = (d(fx, x, 2), d(fx, x, 1));
    let keys = [d2.clone(), d1.clone(), fx.clone()];
    let exclude = vec![fx.clone(), d1.clone(), d2.clone()];
    let collected = match collect::run::<()>(eq, &keys, Options { exact: false, discover: false, distribute_order_term: true }, true, None) {
        Ok(Collection::Expression(v)) => v,
        _ => return Outcome::Decline,
    };
    let names = ["a3", "b3", "c3"];
    let wilds: Vec<Expr> = names.iter().map(|n| Expr::wild(n, exclude.clone())).collect();
    let pattern = add(&[mul2(&wilds[0], &d2), mul2(&wilds[1], &d1), mul2(&wilds[2], fx)]);
    let m = match pattern_match(&collected, &pattern) {
        Ok(Some(m)) => m,
        Ok(None) => return Outcome::NoMatch,
        Err(_) => return Outcome::Decline,
    };
    let mut values = Vec::with_capacity(3);
    for w in &wilds {
        match m.get(w) {
            Some(v) => values.push(v.clone()),
            None => return Outcome::Decline,
        }
    }
    if !values.iter().all(is_polynomial) {
        return Outcome::Decline;
    }
    Outcome::Match([values[0].clone(), values[1].clone(), values[2].clone()])
}

/// `match_2nd_hypergeometric(eq, f(x))`: `(A, B)` with
/// `A = cancel(b3/a3)`, `B = cancel(c3/a3)`; `NoMatch` for its `[]`.
pub fn hypergeometric_ab(eq: &Expr, fx: &Expr, x: &Expr, cancel: Cancel<'_>) -> Outcome<(Expr, Expr)> {
    let [a3, b3, c3] = match deq_match(eq, fx, x) {
        Outcome::Match(v) => v,
        Outcome::NoMatch => return Outcome::NoMatch,
        Outcome::Decline => return Outcome::Decline,
    };
    match nonzero(&a3) {
        Some(true) => {}
        Some(false) => return Outcome::NoMatch,
        None => return Outcome::Decline,
    }
    match (cancel(&div(&b3, &a3)), cancel(&div(&c3, &a3))) {
        (Some(a), Some(b)) => Outcome::Match((a, b)),
        _ => Outcome::Decline,
    }
}

/// The `2nd_linear_airy` match on `eq_high_order_free`: `(b, m)` of
/// `f'' - (b + m*x)*f`.
pub fn airy(eq: &Expr, fx: &Expr, x: &Expr, order: u32, cancel: Cancel<'_>) -> Outcome<(Expr, Expr)> {
    if order != 2 {
        return Outcome::NoMatch;
    }
    let m = match linear_coefficients(eq, fx, x, order) {
        Outcome::Match(m) => m,
        Outcome::NoMatch => return Outcome::NoMatch,
        Outcome::Decline => return Outcome::Decline,
    };
    let (c0, c1, c2) = (&m[1], &m[2], &m[3]);
    match nonzero(c2) {
        Some(true) => {}
        Some(false) => return Outcome::NoMatch,
        None => return Outcome::Decline,
    }
    match crate::assumptions::is_zero(c1) {
        Some(true) => {}
        _ => return Outcome::NoMatch,
    }
    let Some(q) = cancel(&div(c0, c2)) else {
        return Outcome::Decline;
    };
    let df = d(fx, x, 1);
    let exclude = vec![x.clone(), fx.clone(), df];
    let a4 = Expr::wild("a4", exclude.clone());
    let b4 = Expr::wild("b4", exclude);
    let pattern = add(&[a4.clone(), mul2(&b4, x)]);
    match pattern_match(&q, &pattern) {
        Ok(Some(r)) => match (r.get(&a4), r.get(&b4)) {
            (Some(b), Some(slope)) => match nonzero(slope) {
                Some(true) => Outcome::Match((b.clone(), slope.clone())),
                Some(false) => Outcome::NoMatch,
                None => Outcome::Decline,
            },
            _ => Outcome::Decline,
        },
        Ok(None) => Outcome::NoMatch,
        Err(_) => Outcome::Decline,
    }
}

/// The `2nd_linear_bessel` match on `eq_high_order_free`, decided natively
/// only when it is rejected before the factorization: no linear second-order
/// form, a zero leading coefficient, or a leading coefficient free of `x`
/// (the reference `factor(a3).match(a4*(x - b)**b4)` then binds `b4 = 0` or
/// fails, both rejections). Never `Match`.
pub fn bessel_rejected(eq: &Expr, fx: &Expr, x: &Expr) -> Outcome<()> {
    match deq_match(eq, fx, x) {
        Outcome::NoMatch => Outcome::NoMatch,
        Outcome::Decline => Outcome::Decline,
        Outcome::Match([a3, _, _]) => {
            if has_float(&a3) {
                Outcome::Decline
            } else if a3.is_zero() || !has(&a3, x) {
                Outcome::NoMatch
            } else {
                Outcome::Decline
            }
        }
    }
}

/// The order-two power series hints at `x0 = point`.
#[derive(Clone, Debug)]
pub enum PowerSeries2 {
    /// `2nd_power_series_ordinary` with the wild values `[a3, b3, c3]`.
    Ordinary([Expr; 3]),
    /// `2nd_power_series_regular` with `p`, `q`.
    Regular(Expr, Expr),
    /// Neither hint.
    Neither,
}

fn singular(e: &Expr) -> bool {
    [Constant::Infinity, Constant::NaN, Constant::ComplexInfinity, Constant::NegInfinity]
        .into_iter()
        .any(|c| has(e, &Expr::from_kind(Kind::Const(c))))
}

/// The order-two branch of `classify_ode` on `reduced_eq`.
pub fn power_series2(reduced_eq: &Expr, fx: &Expr, x: &Expr, point: &Expr, cancel: Cancel<'_>) -> Option<PowerSeries2> {
    let [a3, b3, c3] = match deq_match(reduced_eq, fx, x) {
        Outcome::Match(v) => v,
        Outcome::NoMatch => return Some(PowerSeries2::Neither),
        Outcome::Decline => return None,
    };
    match nonzero(&a3) {
        Some(true) => {}
        Some(false) => return Some(PowerSeries2::Neither),
        None => return None,
    }
    if has_special(&a3) || has_special(&b3) || has_special(&c3) || has_special(point) || unmodelled(point) {
        return None;
    }
    let p = cancel(&div(&b3, &a3))?;
    let q = cancel(&div(&c3, &a3))?;
    let at = |e: &Expr| crate::subs::subs1(e, x, point);
    if !singular(&at(&p)) && !singular(&at(&q)) {
        return Some(PowerSeries2::Ordinary([a3, b3, c3]));
    }
    let shift = sub(x, point);
    let p = cancel(&mul2(&shift, &p))?;
    if singular(&at(&p)) {
        return Some(PowerSeries2::Neither);
    }
    let q = cancel(&mul2(&pow(&shift, &Expr::int(2)), &q))?;
    if singular(&at(&q)) {
        return Some(PowerSeries2::Neither);
    }
    Some(PowerSeries2::Regular(p, q))
}
