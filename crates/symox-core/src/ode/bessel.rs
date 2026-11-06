//! Native `2nd_linear_bessel` solver.
//!
//! Port of SymPy 1.14 `SecondLinearBessel._matches` and
//! `SecondLinearBessel._get_general_solution`
//! (`sympy/solvers/ode/single.py`): the collected second-order match, its
//! retry on the expanded numerator, the shift of the singular point, the
//! normalization `a3 -> x**2`, the wild matches of the Bessel form and the
//! construction of `x**((1 - c4)/2)*(C1*besselj(..) + C2*bessely(..))`.
//!
//! The collection is the public `collect` kernel. The wild matches are the
//! native port of the reference `Expr.match` ([`crate::matches`]); an input
//! it declines is not modelled. `factor`, `simplify`, `cancel` and `expand`
//! are the public callables of [`Hooks`]; their errors propagate.
//! `Outcome::Decline` leaves the decision to the reference body.

use crate::expr::{Expr, Kind, RelOp};
use crate::ops::{add, div, mul2, pow, sqrt, sub};
use crate::poly::collect::{self, Collection, Options};
use crate::relational::rel;
use crate::special::{besselj, bessely};
use crate::subs::subs1;

use super::algebraic::Hooks;
use super::classify::Outcome;
use super::linear::unmodelled;
use super::second::{d, has_float, has_special, is_polynomial, nonzero};

/// The reference `self.rn`.
#[derive(Clone, Debug, PartialEq)]
pub struct Bessel {
    pub n: Expr,
    pub a4: Expr,
    pub c4: Expr,
    pub d4: Expr,
    pub b4: Expr,
}

/// `None` for no match, `Err(())` for a match the port does not model.
type Matched<T> = Result<Option<T>, ()>;

/// `collect(eq, [f'', f', f]).match(a3*f'' + b3*f' + c3*f)`.
fn deq<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, hooks: &mut H) -> Result<Matched<[Expr; 3]>, H::Error> {
    let (d2, d1) = (d(fx, x, 2), d(fx, x, 1));
    let keys = [d2.clone(), d1.clone(), fx.clone()];
    let exclude = vec![fx.clone(), d1.clone(), d2.clone()];
    let collected = match collect::run::<()>(eq, &keys, Options { exact: false, discover: false, distribute_order_term: true }, true, None) {
        Ok(Collection::Expression(v)) => v,
        _ => return Ok(Err(())),
    };
    let wilds: Vec<Expr> = ["a3", "b3", "c3"].iter().map(|n| Expr::wild(n, exclude.clone())).collect();
    let pattern = add(&[mul2(&wilds[0], &d2), mul2(&wilds[1], &d1), mul2(&wilds[2], fx)]);
    Ok(matched(&collected, &pattern, &[&wilds[0], &wilds[1], &wilds[2]], hooks)?.map(|v| v.map(|v| [v[0].clone(), v[1].clone(), v[2].clone()])))
}

/// `e.match(pattern)` (the native port of the reference matcher) read at
/// `wilds`; a declined match or a missing wild is not modelled.
fn matched<H: Hooks>(e: &Expr, pattern: &Expr, wilds: &[&Expr], _hooks: &mut H) -> Result<Matched<Vec<Expr>>, H::Error> {
    Ok(match crate::matches::reference_match(e, pattern, false) {
        Ok(Some(m)) => wilds.iter().map(|w| m.get(*w).cloned().ok_or(())).collect::<Result<Vec<_>, _>>().map(Some),
        Ok(None) => Ok(None),
        Err(_) => Err(()),
    })
}

/// The reference `if point:`: false only for a numeric zero.
fn truthy(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => !n.is_zero(),
        _ => true,
    }
}

macro_rules! decline {
    ($e:expr) => {
        match $e {
            Ok(v) => v,
            Err(()) => return Ok(Outcome::Decline),
        }
    };
}

/// `SecondLinearBessel._matches` on `eq_high_order_free` of order `order`.
pub fn bessel_match<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, order: u32, hooks: &mut H) -> Result<Outcome<Bessel>, H::Error> {
    if unmodelled(eq) || eq.has_matrix() || eq.has_host() || !x.is_symbol() {
        return Ok(Outcome::Decline);
    }
    let df = d(fx, x, 1);
    let mut r = decline!(deq(eq, fx, x, hooks)?);
    if order == 2 {
        if let Some(v) = &r {
            if !v.iter().all(is_polynomial) {
                let (n, _) = crate::simptools::as_numer_denom(eq);
                let eq = hooks.expand(&n)?;
                if unmodelled(&eq) || eq.has_host() {
                    return Ok(Outcome::Decline);
                }
                r = decline!(deq(&eq, fx, x, hooks)?);
            }
        }
    }
    let Some([a3, b3, c3]) = r else {
        return Ok(Outcome::NoMatch);
    };
    match nonzero(&a3) {
        Some(true) => {}
        Some(false) => return Ok(Outcome::NoMatch),
        None => return Ok(Outcome::Decline),
    }
    let ex = |extra: bool| -> Vec<Expr> {
        let mut v = vec![fx.clone(), df.clone()];
        if extra {
            v.insert(0, x.clone());
        }
        v
    };
    let a = Expr::wild("a", ex(false));
    let b = Expr::wild("b", ex(true));
    let a4 = Expr::wild("a4", ex(true));
    let b4 = Expr::wild("b4", ex(true));
    let c4 = Expr::wild("c4", ex(true));
    let d4 = Expr::wild("d4", ex(true));

    // leading coeff of f(x).diff(x, 2)
    let fa3 = hooks.factor(&a3)?;
    let coeff = decline!(matched(&fa3, &mul2(&a4, &pow(&sub(x, &b), &b4)), &[&a4, &b, &b4], hooks)?);
    let Some(coeff) = coeff else {
        return Ok(Outcome::NoMatch);
    };
    let (ca4, point, cb4) = (&coeff[0], &coeff[1], &coeff[2]);
    // coeff[b4] == 0 means constant coefficient
    match nonzero(cb4) {
        Some(true) => {}
        Some(false) => return Ok(Outcome::NoMatch),
        None => return Ok(Outcome::Decline),
    }
    let (mut a3, mut b3, mut c3) = (a3, b3, c3);
    if has_float(point) {
        return Ok(Outcome::Decline);
    }
    if truthy(point) {
        if [&a3, &b3, &c3].iter().any(|e| has_special(e)) || has_special(point) || unmodelled(point) {
            return Ok(Outcome::Decline);
        }
        let shifted = add(&[x.clone(), point.clone()]);
        a3 = hooks.simplify(&subs1(&a3, x, &shifted))?;
        b3 = hooks.simplify(&subs1(&b3, x, &shifted))?;
        c3 = hooks.simplify(&subs1(&c3, x, &shifted))?;
    }
    // making a3 in the form of x**2
    let scale = mul2(ca4, &pow(x, &add(&[Expr::int(-2), cb4.clone()])));
    let _a3 = hooks.cancel(&div(&a3, &scale))?;
    let b3 = hooks.cancel(&div(&b3, &scale))?;
    let c3 = hooks.cancel(&div(&c3, &scale))?;
    // checking if b3 is of form c*(x-b)
    let fb3 = hooks.factor(&b3)?;
    let Some(coeff1) = decline!(matched(&fb3, &mul2(&a4, x), &[&a4], hooks)?) else {
        return Ok(Outcome::NoMatch);
    };
    // c3 maybe of very complex form so I am simply checking (a - b) form
    let ec3 = hooks.expand(&c3)?;
    let Some(coeff2_) = decline!(matched(&ec3, &sub(&a, &b), &[&a, &b], hooks)?) else {
        return Ok(Outcome::NoMatch);
    };
    // matching with standard form for c3
    let pattern = mul2(&pow(&c4, &Expr::int(2)), &pow(x, &mul2(&Expr::int(2), &a4)));
    let fa = hooks.factor(&coeff2_[0])?;
    let Some(coeff2) = decline!(matched(&fa, &pattern, &[&c4, &a4], hooks)?) else {
        return Ok(Outcome::NoMatch);
    };
    let n = match nonzero(&coeff2_[1]) {
        None => return Ok(Outcome::Decline),
        Some(false) => Expr::zero(),
        Some(true) => match decline!(matched(&hooks.factor(&coeff2_[1])?, &pow(&d4, &Expr::int(2)), &[&d4], hooks)?) {
            Some(v) => v[0].clone(),
            // `None[d4]` raises TypeError in the reference.
            None => return Ok(Outcome::Decline),
        },
    };
    Ok(Outcome::Match(Bessel { n, a4: coeff2[0].clone(), d4: coeff2[1].clone(), c4: coeff1[0].clone(), b4: point.clone() }))
}

/// `SecondLinearBessel._get_general_solution` with the numbered constants
/// `c1`, `c2`: `Eq(f(x), ...)`. A `c4` that is not a rational number
/// declines (the reference `Rational(1 - c4, 2)` rejects a symbol and
/// rationalizes a float).
pub fn bessel_solution(m: &Bessel, fx: &Expr, x: &Expr, c1: &Expr, c2: &Expr) -> Option<Expr> {
    let Kind::Num(c4) = m.c4.kind() else {
        return None;
    };
    if !c4.is_rational() || [&m.n, &m.a4, &m.d4, &m.b4].iter().any(|e| has_special(e) || unmodelled(e) || e.has_host()) {
        return None;
    }
    let quarter = Expr::rational(1, 4);
    let n = sqrt(&add(&[pow(&m.n, &Expr::int(2)), mul2(&quarter, &pow(&sub(&m.c4, &Expr::one()), &Expr::int(2)))]));
    let exponent = div(&sub(&Expr::one(), &m.c4), &Expr::int(2));
    let order = div(&n, &m.d4);
    let arg = div(&mul2(&m.a4, &pow(x, &m.d4)), &m.d4);
    let body = mul2(&pow(x, &exponent), &add(&[mul2(c1, &besselj(&order, &arg)), mul2(c2, &bessely(&order, &arg))]));
    let rhs = subs1(&body, x, &sub(x, &m.b4));
    Some(rel(RelOp::Eq, fx, &rhs))
}

