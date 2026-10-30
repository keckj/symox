//! Native `nth_order_reducible` and `Riccati_special_minus2` matches of
//! `classify_ode` (SymPy 1.14 `NthOrderReducible._matches` and the
//! `RiccatiSpecial` pattern match).
//!
//! The substitution is the kernel of the public `subs`, the wild match the
//! matcher of the public `Expr.match`. `Decline` leaves the decision to the
//! reference matcher.

use crate::diff::derivative_unevaluated;
use crate::expr::{Expr, Kind};
use crate::ops::{add, mul, mul2, pow};
use crate::pattern::pattern_match;
use crate::traverse::{has_structural as has, preorder};

use super::classify::Outcome;
use super::linear::unmodelled;

/// `NthOrderReducible._matches` on `eq_preprocessed`: the smallest order of
/// the derivatives of `fx` in `x` alone when there are at least two such
/// orders and `fx` appears only through derivatives of at least that order.
pub fn nth_order_reducible(eq: &Expr, fx: &Expr, x: &Expr) -> Outcome<u32> {
    if unmodelled(eq) || eq.has_matrix() {
        return Outcome::Decline;
    }
    let mut ords: Vec<u32> = Vec::new();
    preorder(eq, &mut |s| {
        if let Kind::Derivative(e, vars) = s.kind() {
            if e == fx && vars.len() == 1 && vars[0].0 == *x && !ords.contains(&vars[0].1) {
                ords.push(vars[0].1);
            }
        }
    });
    if ords.len() < 2 {
        return Outcome::NoMatch;
    }
    let smallest = *ords.iter().min().expect("two orders");
    let dummy = Expr::dummy("D", Default::default());
    let replaced = crate::subs::subs1(eq, &derivative_unevaluated(fx, vec![(x.clone(), smallest)]), &dummy);
    if has(&replaced, fx) {
        Outcome::NoMatch
    } else {
        Outcome::Match(smallest)
    }
}

/// `collected.match(a*df + b*f**2 + c*f/x + d/x**2)` with the wilds of
/// `RiccatiSpecial` (`a`, `b` excluding `x`, `f`, `df` and `0`; `c`, `d`
/// excluding `x`, `f`, `df`), on `collected = eq_expanded.collect([df, f],
/// func=cancel)`: `[a, b, c, d]`.
pub fn riccati_special(collected: &Expr, fx: &Expr, x: &Expr, df: &Expr) -> Outcome<[Expr; 4]> {
    if unmodelled(collected) || collected.has_matrix() {
        return Outcome::Decline;
    }
    let base = vec![x.clone(), fx.clone(), df.clone()];
    let mut with_zero = base.clone();
    with_zero.push(Expr::zero());
    let a = Expr::wild("a", with_zero.clone());
    let b = Expr::wild("b", with_zero);
    let c = Expr::wild("c", base.clone());
    let d = Expr::wild("d", base);
    let pattern = add(&[
        mul2(&a, df),
        mul2(&b, &pow(fx, &Expr::int(2))),
        mul(&[c.clone(), fx.clone(), pow(x, &Expr::int(-1))]),
        mul2(&d, &pow(x, &Expr::int(-2))),
    ]);
    match pattern_match(collected, &pattern) {
        Ok(Some(m)) => match (m.get(&a), m.get(&b), m.get(&c), m.get(&d)) {
            (Some(a), Some(b), Some(c), Some(d)) => Outcome::Match([a.clone(), b.clone(), c.clone(), d.clone()]),
            _ => Outcome::Decline,
        },
        Ok(None) => Outcome::NoMatch,
        Err(_) => Outcome::Decline,
    }
}
