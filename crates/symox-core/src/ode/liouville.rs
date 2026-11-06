//! Native `Liouville` and `2nd_nonlinear_autonomous_conserved` matches of
//! `classify_ode` (SymPy 1.14 `SinglePatternODESolver._matches` at order
//! two).
//!
//! The collection is the kernel of the public `collect` with the public
//! `cancel` as `func`, the wild match the matcher of the public
//! `Expr.match`, and `expand` the public callable. The solvers' `_verify`
//! steps stay with the caller. `Decline` leaves the decision to the
//! reference matcher.

use crate::diff::derivative_unevaluated;
use crate::expr::Expr;
use crate::ops::{add, div, mul2, pow};
use crate::matches::reference_match;

use super::classify::{blocked, coeff, collected, Hooks, Outcome};
use super::linear::unmodelled;

/// The order-two pattern matches of one equation.
#[derive(Clone, Debug)]
pub struct SecondOrderPatterns {
    /// `(d, e, k)` of `d*f'' + e*f'**2 + k*f'`.
    pub liouville: Outcome<[Expr; 3]>,
    /// `fy` of `f'' + fy`.
    pub conserved: Outcome<Expr>,
}

fn matched<const N: usize>(e: &Expr, pattern: &Expr, wilds: &[Expr; N]) -> Outcome<[Expr; N]> {
    match reference_match(e, pattern, false) {
        Ok(Some(m)) => {
            let mut out = Vec::with_capacity(N);
            for w in wilds {
                match m.get(w) {
                    Some(v) => out.push(v.clone()),
                    None => return Outcome::Decline,
                }
            }
            match out.try_into() {
                Ok(v) => Outcome::Match(v),
                Err(_) => Outcome::Decline,
            }
        }
        Ok(None) => Outcome::NoMatch,
        Err(_) => Outcome::Decline,
    }
}

/// Both matches on `eq_expanded` (the expanded prepared equation) of order
/// two in `fx = f(x)`. As in the reference, the collection is on
/// `[f', f]`; the `Liouville` pattern has a wild coefficient on `f''` and
/// matches the collected equation, the conserved pattern has coefficient
/// one and matches `collect(expand(eq/eq.coeff(f'')))`.
pub fn second_order_patterns(eq_expanded: &Expr, fx: &Expr, x: &Expr, hooks: &mut Hooks<'_>) -> SecondOrderPatterns {
    let declined = SecondOrderPatterns { liouville: Outcome::Decline, conserved: Outcome::Decline };
    if unmodelled(eq_expanded) || eq_expanded.has_matrix() || blocked(eq_expanded) {
        return declined;
    }
    let d1 = derivative_unevaluated(fx, vec![(x.clone(), 1)]);
    let d2 = derivative_unevaluated(fx, vec![(x.clone(), 2)]);
    let liouville = match collected(eq_expanded, &d1, fx, hooks) {
        Some(c) => {
            let d = Expr::wild("d", vec![d1.clone(), d2.clone()]);
            let e = Expr::wild("e", vec![d1.clone()]);
            let k = Expr::wild("k", vec![d1.clone()]);
            let pattern = add(&[mul2(&d, &d2), mul2(&e, &pow(&d1, &Expr::int(2))), mul2(&k, &d1)]);
            matched(&c, &pattern, &[d, e, k])
        }
        None => Outcome::Decline,
    };
    let lead = coeff(eq_expanded, &d2);
    let conserved = if lead.is_zero() || blocked(&lead) {
        Outcome::Decline
    } else {
        match (hooks.expand)(&div(eq_expanded, &lead)).and_then(|e| collected(&e, &d1, fx, hooks)) {
            Some(c) => {
                let fy = Expr::wild("fy", vec![Expr::zero(), d1.clone(), d2.clone()]);
                let pattern = add(&[d2.clone(), fy.clone()]);
                match matched(&c, &pattern, &[fy]) {
                    Outcome::Match([v]) => Outcome::Match(v),
                    Outcome::NoMatch => Outcome::NoMatch,
                    Outcome::Decline => Outcome::Decline,
                }
            }
            None => Outcome::Decline,
        }
    };
    SecondOrderPatterns { liouville, conserved }
}
