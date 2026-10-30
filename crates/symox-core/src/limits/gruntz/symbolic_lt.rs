//! Leading terms with symbolic (irrational) powers of the mrv variable,
//! where the integer-exponent series engine does not apply (`w**log(3)`).
//!
//! The reference reaches these through `as_leading_term` of each node:
//!
//! * `log(a)` with `a -> 1` is `log(1 + (a - 1)) ~ a - 1`
//!   (`log._eval_as_leading_term`), so `log(1 + w**p)/w ~ w**(p - 1)`;
//! * `F(g)` with `g -> 0` and `F` analytic with `F(0) = 0` is
//!   `F^(m)(0)/m! * g**m` (`sin(w**p) ~ w**p`);
//! * a sum whose leading terms cancel (`exp(w**p) - 1`) is expanded:
//!   each analytic `F(g)` with `g -> 0` is replaced by its Taylor
//!   polynomial of order `k` plus `rho*g**(k + 1)`, `rho` a fresh symbol
//!   standing for the remainder ratio (which tends to the next Taylor
//!   coefficient, a finite number). The leading term of the result is
//!   accepted only when it is free of `rho`: the remainder terms are then
//!   of strictly higher order. Otherwise the next order is tried.
//!
//! Every step is exact; an undecided sign or a remaining cancellation is
//! reported as unknown, never guessed.

use super::super::{unknown, LRes};
use super::Gruntz;
use crate::expr::{Expr, Func, Kind};
use crate::ops::{add, mul2, pow, sub};
use crate::series::engine::depends;
use crate::subs::{subs1, xreplace};
use rustc_hash::FxHashMap;

/// Analytic at 0, with the Taylor series of the engine.
fn analytic_at_zero(f: &Func) -> bool {
    matches!(
        f,
        Func::Exp
            | Func::Sin
            | Func::Cos
            | Func::Tan
            | Func::Sinh
            | Func::Cosh
            | Func::Tanh
            | Func::Asin
            | Func::Atan
            | Func::Asinh
            | Func::Atanh
            | Func::Erf
    )
}

const MAX_TAYLOR: i64 = 4;

thread_local! {
    /// Nesting of the Taylor replacement (each level replaces nodes, the
    /// bound keeps the recursion finite).
    static TAYLOR_DEPTH: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

impl Gruntz {
    /// `log(a)` with `lt(a) = 1`: the leading term of `a - 1`, which must
    /// vanish.
    pub(super) fn log_near_one(&mut self, a: &Expr, w: &Expr, logw: &Expr) -> LRes<(Expr, Expr)> {
        let r = crate::expand::expand(&sub(a, &Expr::one()));
        if r.is_zero() || !depends(&r, w) {
            return unknown("cancellation in log leading term");
        }
        let (c, e) = self.lt_fallback(&r, w, logw)?;
        if self.sign(&e)? != 1 {
            return unknown("cancellation in log leading term");
        }
        Ok((c, e))
    }

    /// `F(g)` with `lt(g) = c1*w**e1`, `e1 > 0`, `F` analytic at 0 and
    /// `F(0) = 0`: `a_m*c1**m*w**(m*e1)` from the first nonzero Taylor
    /// coefficient `a_m`.
    pub(super) fn func_at_zero(&mut self, f: &Func, c1: &Expr, e1: &Expr) -> LRes<(Expr, Expr)> {
        if !analytic_at_zero(f) {
            return unknown("zero of a function in symbolic leading term");
        }
        let s = Expr::dummy("s", crate::symbol::Assumptions::NONE);
        let fs = crate::functions::func(f.clone(), vec![s.clone()]);
        let lt = match crate::series::leading_term(&fs, &s, &Expr::zero(), false) {
            Ok(lt) => lt,
            Err(_) => return unknown("zero of a function in symbolic leading term"),
        };
        let m = lt.exp.clone();
        if !m.is_integer() || !m.is_positive() || depends(&lt.coeff, &s) || lt.coeff.is_zero() {
            return unknown("zero of a function in symbolic leading term");
        }
        let me = Expr::num(m);
        Ok((mul2(&lt.coeff, &pow(c1, &me)), mul2(&me, e1)))
    }

    /// A sum whose leading terms cancel: Taylor replacement of its
    /// analytic nodes at a vanishing argument (see the module notes).
    pub(super) fn cancelled_sum(&mut self, f: &Expr, w: &Expr, logw: &Expr) -> LRes<(Expr, Expr)> {
        if TAYLOR_DEPTH.with(|d| d.get()) >= 2 {
            return unknown("cancellation in symbolic leading term");
        }
        // one shared variable for the outer functions F(s) of every node
        let s = Expr::dummy("s", crate::symbol::Assumptions::NONE);
        let mut nodes: Vec<(Expr, Expr, Expr)> = vec![];
        self.vanishing_nodes(f, w, logw, &s, &mut nodes);
        if nodes.is_empty() {
            return unknown("cancellation in symbolic leading term");
        }
        TAYLOR_DEPTH.with(|d| d.set(d.get() + 1));
        let r = self.taylor_orders(f, w, logw, &s, &nodes);
        TAYLOR_DEPTH.with(|d| d.set(d.get() - 1));
        r
    }

    fn taylor_orders(&mut self, f: &Expr, w: &Expr, logw: &Expr, s: &Expr, nodes: &[(Expr, Expr, Expr)]) -> LRes<(Expr, Expr)> {
        for k in 1..=MAX_TAYLOR {
            let rho = Expr::dummy("rho", crate::symbol::Assumptions::NONE);
            let mut m = FxHashMap::default();
            let mut ok = true;
            for (node, outer, g) in nodes {
                // outer(s) is F(s) (or log(1 + s) for a log near one)
                let t = match crate::series::series(outer, s, &Expr::zero(), k + 1, false) {
                    Ok(t) => crate::series::order::remove_o(&t),
                    Err(_) => {
                        ok = false;
                        break;
                    }
                };
                let poly = subs1(&t, s, g);
                let rem = mul2(&rho, &pow(g, &Expr::int(k + 1)));
                m.insert(node.clone(), add(&[poly, rem]));
            }
            if !ok {
                return unknown("cancellation in symbolic leading term");
            }
            let g = crate::expand::expand(&xreplace(f, &m));
            if g.is_zero() || !depends(&g, w) {
                continue;
            }
            match self.lt_fallback(&g, w, logw) {
                Ok((c, e)) if !depends(&c, &rho) && !depends(&c, s) && !depends(&e, s) && !c.is_zero() => return Ok((c, e)),
                Ok(_) => continue,
                Err(_) => continue,
            }
        }
        unknown("cancellation in symbolic leading term")
    }

    /// Outermost nodes `F(g)` (`F` analytic at 0, `g -> 0`) and `log(a)`
    /// (`a -> 1`), with `F(s)` and the vanishing argument.
    fn vanishing_nodes(&mut self, e: &Expr, w: &Expr, logw: &Expr, s: &Expr, out: &mut Vec<(Expr, Expr, Expr)>) {
        if !depends(e, w) || out.iter().any(|(n, _, _)| n == e) {
            return;
        }
        if let Kind::Func(f, a) = e.kind() {
            if a.len() == 1 && depends(&a[0], w) {
                if analytic_at_zero(f) {
                    if let Ok((_, e1)) = self.lt_fallback(&a[0], w, logw) {
                        if self.sign(&e1).ok() == Some(1) {
                            out.push((e.clone(), crate::functions::func(f.clone(), vec![s.clone()]), a[0].clone()));
                            return;
                        }
                    }
                } else if matches!(f, Func::Log) {
                    if let Ok((c1, e1)) = self.lt_fallback(&a[0], w, logw) {
                        if c1 == Expr::one() && e1.is_zero() {
                            let h = crate::expand::expand(&sub(&a[0], &Expr::one()));
                            let log1p = crate::functions::func(Func::Log, vec![add(&[Expr::one(), s.clone()])]);
                            out.push((e.clone(), log1p, h));
                            return;
                        }
                    }
                }
            }
        }
        for a in e.args() {
            self.vanishing_nodes(&a, w, logw, s, out);
        }
    }
}
