//! Dominant terms of a one-variable `Order` whose argument has no power
//! leading term (`O(exp(1/x) + 1)`, `O(2**x + factorial(x), (x, oo))`).
//!
//! `Order.__new__` moves the point to 0, factors a sum and, while it stays a
//! sum, keeps the terms of `Add.extract_leading_order`: the order of every
//! term, with a term dropped when the order of a kept one contains it.
//! `Order.contains` decides `O(f) ⊇ O(g)` by `limit(f/g, x, 0+) != 0`
//! (Gruntz only, `heuristics=False`). A term that is a power with an
//! exponent depending on `x` and no power leading term is written
//! `exp(k*log(b))` (the `PoleError` branch of `Order.__new__`), which gives
//! `O(exp(x*log(3)), (x, oo))` for `O(3**x, (x, oo))`.
//!
//! Admitted natively: sums whose factored form stays a sum and where every
//! containment needed is decided by the native limit. An undecided
//! containment keeps the previous result (the whole sum), and is a gap,
//! never a proof that the terms are incomparable.

use std::cell::Cell;

use super::engine::depends;
use crate::expr::{Expr, Func, Kind};
use crate::ops::{add, div, mul2};

thread_local! {
    static ACTIVE: Cell<bool> = const { Cell::new(false) };
}

/// Runs `f` unless a dominant-term computation is already running on this
/// thread (the limits computed for the containments build orders too).
fn guarded<T>(f: impl FnOnce() -> Option<T>) -> Option<T> {
    if ACTIVE.with(|a| a.replace(true)) {
        return None;
    }
    let r = f();
    ACTIVE.with(|a| a.set(false));
    r
}

/// `exp(k*log(b))` for a power `b**k` whose exponent depends on `x`.
pub(super) fn pole_power(e0: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Pow(b, k) = e0.kind() else { return None };
    if !depends(k, x) {
        return None;
    }
    let arg = crate::simp::expand_log(&mul2(k, &crate::functions::log(b)), false);
    Some(crate::functions::func(Func::Exp, vec![arg]))
}

/// `Some(true)` when `O(f)` contains `O(g)` at `x = 0+`, `Some(false)` when
/// it does not, `None` when the native limit declines.
fn contains_at_zero(f: &Expr, g: &Expr, x: &Expr) -> Option<bool> {
    if f == g {
        return Some(true);
    }
    let ratio = div(f, g);
    let l = crate::limits::limit(&ratio, x, &Expr::zero(), crate::limits::Dir::Plus).ok()?;
    if depends(&l, x) {
        return None;
    }
    Some(!l.is_zero())
}

/// The dominant terms of the sum `e0` at `x = 0+` (see the module notes);
/// `None` keeps the previous result.
pub(super) fn dominant_terms(e0: &Expr, x: &Expr) -> Option<Expr> {
    if !e0.is_add() || e0.may_have_order() {
        return None;
    }
    guarded(|| {
        let f = crate::poly::api::factor(e0).ok()?;
        if !f.is_add() {
            return None;
        }
        let zero = [(x.clone(), Expr::zero())];
        let mut kept: Vec<Expr> = Vec::new();
        for t in f.args() {
            let o = super::order::order(&t, &zero);
            let (oe, _) = super::order::order_parts(&o)?;
            if oe.is_zero() {
                continue;
            }
            let mut drop = false;
            for k in &kept {
                if *k != oe && contains_at_zero(k, &oe, x)? {
                    drop = true;
                    break;
                }
            }
            if drop {
                continue;
            }
            let mut next = vec![oe.clone()];
            for k in &kept {
                if *k != oe && contains_at_zero(&oe, k, x)? {
                    continue;
                }
                next.push(k.clone());
            }
            kept = next;
        }
        let r = add(&kept);
        (r != *e0).then_some(r)
    })
}
