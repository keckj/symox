//! `_inverse_laplace_irrational`: the irrational algebraic functions of
//! Bateman (1954), section 5.3.
//!
//! Every factor of `F(s)` must match `(a*s**m + b)**n` with `a, b, m, n`
//! free of `s`; the factors are split into zeros (`n > 0`) and poles
//! (`n < 0`), sorted by `(n, b != 0, b)`, and the first structurally
//! matching rule decides. Its inner condition (a positive shift) failing
//! leaves the function to the next stage, as the reference's `elif` chain
//! does.
//!
//! Where the reference formula is wrong the rule is corrected, never
//! extended: rules 5.3.14 and 5.3.15 divide the shift by the leading
//! coefficient of the cubed pole and scale the whole value (the reference
//! uses `b` and applies the factor to one term only, which is right only for
//! a unit coefficient), rule 5.3.14 requires the zero to be `sqrt(a*s)`,
//! rule 5.3.2 requires the square-root pole to be `sqrt(a*s)`, and the
//! closed forms of rules 5.3.6, 5.3.10 and 5.3.16 are the correct ones (the
//! reference values fail `L[f](s) == F(s)`). The values agree with the
//! reference wherever its value is correct.

use super::util::*;
use crate::expr::{Expr, Func, Kind};
use crate::functions::{exp, func, heaviside};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};

type R = (Expr, Expr);

#[derive(Clone, Debug)]
struct Term {
    a: Expr,
    m: Expr,
    b: Expr,
    n: Expr,
}

/// `base == a*s**m + b` with `a, m, b` free of `s`.
fn match_binomial(base: &Expr, s: &Expr) -> Option<(Expr, Expr, Expr)> {
    let mut constant = Vec::new();
    let mut dependent = Vec::new();
    for term in terms(base) {
        if has(&term, s) { dependent.push(term) } else { constant.push(term) }
    }
    if dependent.len() != 1 {
        return None;
    }
    let (a, g) = split_mul(&dependent[0], s);
    let m = power_of(&g, s)?;
    if m.is_zero() {
        return None;
    }
    Some((a, m, add(&constant)))
}

/// `f == (a*s**m + b)**n`, read through `as_base_exp` like the reference's
/// `Pow.matches`.
fn match_factor(f: &Expr, s: &Expr) -> Option<Term> {
    let (base, n) = match f.kind() {
        Kind::Pow(b, e) => (b.clone(), e.clone()),
        _ => (f.clone(), Expr::one()),
    };
    if has(&n, s) {
        return None;
    }
    let (a, m, b) = match_binomial(&base, s)?;
    Some(Term { a, m, b, n })
}

/// `x < y` decided, `None` where the reference's comparison of relationals
/// would raise.
fn less(x: &Expr, y: &Expr) -> Option<bool> {
    let d = sub(x, y);
    if let Some(n) = d.as_number() {
        return Some(n.is_negative());
    }
    match crate::assumptions::is_negative(&d) {
        Some(true) => Some(true),
        Some(false) if crate::assumptions::is_real(&d) == Some(true) => Some(false),
        _ => None,
    }
}

/// `key(x) < key(y)` for the tuple key `(n, b != 0, b)`, deciding only
/// the comparisons the reference makes.
fn key_less(x: &Term, y: &Term) -> Option<bool> {
    if x.n != y.n {
        return less(&x.n, &y.n);
    }
    let (bx, by) = (!x.b.is_zero(), !y.b.is_zero());
    if bx != by {
        return Some(!bx);
    }
    if x.b == y.b {
        return Some(false);
    }
    less(&x.b, &y.b)
}

/// Stable insertion sort by `key_less`; `None` on an undecidable comparison.
fn sorted(mut v: Vec<Term>) -> Option<Vec<Term>> {
    for i in 1..v.len() {
        let mut j = i;
        while j > 0 && key_less(&v[j], &v[j - 1])? {
            v.swap(j, j - 1);
            j -= 1;
        }
    }
    Some(v)
}

/// `x0*x1*...` evaluated left to right like the reference's operators, so
/// that a number times a sum is distributed exactly where it is there.
fn ml(xs: &[Expr]) -> Expr {
    xs.iter().skip(1).fold(xs[0].clone(), |acc, x| mul(&[acc, x.clone()]))
}

fn erf(x: &Expr) -> Expr {
    func(Func::Erf, vec![x.clone()])
}
fn erfc(x: &Expr) -> Expr {
    func(Func::Erfc, vec![x.clone()])
}

pub(crate) fn irrational(f: &Expr, s: &Expr, t: &Expr) -> Option<R> {
    let mut constants = Expr::one();
    let mut zeros = Vec::new();
    let mut poles = Vec::new();
    for fct in crate::sortkey::ordered_factors(f) {
        let term = match_factor(&fct, s)?;
        if term.a.is_zero() {
            constants = mul(&[constants, fct]);
        } else if pos(&term.n) {
            zeros.push(term);
        } else if negv(&term.n) {
            poles.push(term);
        } else {
            return None;
        }
    }
    let poles = sorted(poles)?;
    let zeros = sorted(zeros)?;
    let c = constants;
    let half = rat(1, 2);
    let one = Expr::one();
    let pi = Expr::pi();
    let st = sqrt(t);
    let sqrt_pi = sqrt(&pi);
    let is = |e: &Expr, v: &Expr| e == v;
    let q = |n: i64| int(n);
    let sq = |e: &Expr| pow(e, &q(2));
    // exp(a**2*t)*erfc(a*sqrt(t))
    let ee = |a: &Expr| mul(&[exp(&mul(&[sq(a), t.clone()])), erfc(&mul(&[a.clone(), st.clone()]))]);
    let p = &poles;
    let z = &zeros;
    let mut result: Option<Expr> = None;

    if p.len() == 1 && z.is_empty() {
        let p0 = &p[0];
        if is(&p0.n, &q(-1)) && is(&p0.m, &half) {
            let a_ = div(&p0.b, &p0.a);
            let k_ = mul(&[div(&one, &p0.a), c.clone()]);
            if pos(&a_) {
                result = Some(sub(&div(&div(&k_, &sqrt_pi), &st), &mul(&[k_.clone(), a_.clone(), ee(&a_)])));
            }
        } else if is(&p0.n, &q(-2)) && is(&p0.m, &half) {
            let a_sq = div(&p0.b, &p0.a);
            let a_ = sq(&a_sq);
            let k_ = mul(&[div(&one, &sq(&p0.a)), c.clone()]);
            if pos(&a_sq) {
                let at = mul(&[a_.clone(), t.clone()]);
                // the reference has 1 - ... + (1 - 2*a*t)*exp(a*t)*(erf(...) - 1)
                result = Some(mul(&[k_, add(&[
                    neg(&mul(&[div(&q(2), &sqrt_pi), sqrt(&a_), st.clone()])),
                    mul(&[add(&[one.clone(), mul(&[q(2), at.clone()])]), exp(&at), erfc(&mul(&[sqrt(&a_), st.clone()]))]),
                ])]));
            }
        } else if is(&p0.n, &q(-3)) && is(&p0.m, &half) {
            let a_ = div(&p0.b, &p0.a);
            let k_ = mul(&[div(&one, &pow(&p0.a, &q(3))), c.clone()]);
            if pos(&a_) {
                let a2t = mul(&[sq(&a_), t.clone()]);
                result = Some(mul(&[k_, sub(
                    &mul(&[div(&q(2), &sqrt_pi), add(&[a2t.clone(), one.clone()]), st.clone()]),
                    &mul(&[a_.clone(), t.clone(), add(&[mul(&[q(2), a2t.clone()]), q(3)]), ee(&a_)]),
                )]));
            }
        } else if is(&p0.n, &q(-4)) && is(&p0.m, &half) {
            let a_ = div(&p0.b, &p0.a);
            let k_ = div(&mul(&[div(&one, &pow(&p0.a, &q(4))), c.clone()]), &q(3));
            if pos(&a_) {
                let poly = add(&[mul(&[q(4), pow(&a_, &q(4)), sq(t)]), mul(&[q(12), sq(&a_), t.clone()]), q(3)]);
                result = Some(mul(&[k_, sub(
                    &mul(&[t.clone(), poly, ee(&a_)]),
                    // the reference has a**3*t**(5/2) for a*t**(3/2)
                    &mul(&[div(&q(2), &sqrt_pi), a_.clone(), pow(t, &rat(3, 2)), add(&[mul(&[q(2), sq(&a_), t.clone()]), q(5)])]),
                )]));
            }
        } else if is(&p0.n, &neg(&half)) && is(&p0.m, &q(2)) {
            let a_ = sqrt(&div(&p0.b, &p0.a));
            let k_ = mul(&[div(&one, &sqrt(&p0.a)), c.clone()]);
            result = Some(mul(&[k_, besselj(&Expr::zero(), &mul(&[a_, t.clone()]))]));
        }
    } else if p.len() == 1 && z.len() == 1 {
        let (p0, z0) = (&p[0], &z[0]);
        if is(&p0.n, &q(-3)) && is(&p0.m, &half) && is(&z0.n, &half) && is(&z0.m, &one) && z0.b.is_zero() {
            // sqrt(az*s)/(a0*sqrt(s) + b0)**3
            let a_ = div(&p0.b, &p0.a);
            let k_ = mul(&[div(&sqrt(&z0.a), &pow(&p0.a, &q(3))), c.clone()]);
            let a2t = mul(&[sq(&a_), t.clone()]);
            result = Some(mul(&[k_, sub(
                &mul(&[add(&[mul(&[q(2), pow(&a_, &q(4)), sq(t)]), mul(&[q(5), a2t.clone()]), one.clone()]), ee(&a_)]),
                &mul(&[div(&q(2), &sqrt_pi), a_.clone(), add(&[a2t, q(2)]), st.clone()]),
            )]));
        } else if is(&p0.n, &q(-1)) && is(&p0.m, &one) && is(&z0.n, &half) && is(&z0.m, &one) {
            let a_ = div(&z0.b, &z0.a);
            let b_ = div(&p0.b, &p0.a);
            let k_ = mul(&[div(&sqrt(&z0.a), &p0.a), c.clone()]);
            let d = sub(&a_, &b_);
            result = Some(mul(&[k_, add(&[
                div(&div(&exp(&mul(&[neg(&a_), t.clone()])), &st), &sqrt_pi),
                mul(&[sqrt(&d), exp(&mul(&[neg(&b_), t.clone()])), erf(&mul(&[sqrt(&d), st.clone()]))]),
            ])]));
        }
    } else if p.len() == 2 && z.is_empty() {
        let (p0, p1) = (&p[0], &p[1]);
        if is(&p0.n, &q(-1)) && is(&p0.m, &one) && is(&p1.n, &neg(&half)) && is(&p1.m, &one) && p1.b.is_zero() {
            // 5.3.1
            let a_ = neg(&div(&p0.b, &p0.a));
            let k_ = mul(&[div(&div(&one, &sqrt(&p1.a)), &p0.a), c.clone()]);
            if pos(&a_) {
                result = Some(mul(&[div(&k_, &sqrt(&a_)), exp(&mul(&[a_.clone(), t.clone()])), erf(&mul(&[sqrt(&a_), st.clone()]))]));
            }
        } else if is(&p0.n, &q(-1)) && is(&p0.m, &one) && p0.b.is_zero() && is(&p1.n, &q(-1)) && is(&p1.m, &half) {
            // 5.3.5
            let a_ = div(&p1.b, &p1.a);
            let k_ = mul(&[div(&div(&div(&one, &p0.a), &p1.a), &a_), c.clone()]);
            if pos(&a_) {
                result = Some(mul(&[k_, sub(&one, &ee(&a_))]));
            }
        } else if is(&p0.n, &q(-1)) && is(&p0.m, &half) && is(&p1.n, &neg(&half)) && is(&p1.m, &one) && p1.b.is_zero() {
            // 5.3.7
            let a_ = div(&p0.b, &p0.a);
            let k_ = mul(&[div(&one, &mul(&[p0.a.clone(), sqrt(&p1.a)])), c.clone()]);
            if pos(&a_) {
                result = Some(mul(&[k_, ee(&a_)]));
            }
        } else if is(&p0.n, &rat(-3, 2)) && is(&p0.m, &one) && p0.b.is_zero() && is(&p1.n, &q(-1)) && is(&p1.m, &half) {
            // 5.3.8
            let a_ = div(&p1.b, &p1.a);
            let k_ = mul(&[div(&div(&one, &mul(&[pow(&p0.a, &rat(3, 2)), p1.a.clone()])), &sq(&a_)), c.clone()]);
            if pos(&a_) {
                result = Some(mul(&[k_, add(&[
                    mul(&[div(&q(2), &sqrt_pi), a_.clone(), st.clone()]),
                    ee(&a_),
                    neg(&one),
                ])]));
            }
        } else if is(&p0.n, &q(-2)) && is(&p0.m, &half) && is(&p1.n, &q(-1)) && is(&p1.m, &one) && p1.b.is_zero() {
            // 5.3.11
            let a_sq = div(&p0.b, &p0.a);
            let a_ = sq(&a_sq);
            let k_ = mul(&[div(&div(&one, &sq(&p0.a)), &p1.a), c.clone()]);
            if pos(&a_sq) {
                let inv = div(&one, &a_);
                result = Some(mul(&[k_, add(&[
                    inv.clone(),
                    mul(&[sub(&mul(&[q(2), t.clone()]), &inv), exp(&mul(&[a_.clone(), t.clone()])), erfc(&mul(&[sqrt(&a_), st.clone()]))]),
                    neg(&mul(&[div(&div(&q(2), &sqrt_pi), &sqrt(&a_)), st.clone()])),
                ])]));
            }
        } else if is(&p0.n, &q(-2)) && is(&p0.m, &half) && is(&p1.n, &neg(&half)) && is(&p1.m, &one) && p1.b.is_zero() {
            // 5.3.12
            let a_ = div(&p0.b, &p0.a);
            let k_ = mul(&[div(&div(&one, &sq(&p0.a)), &sqrt(&p1.a)), c.clone()]);
            if pos(&a_) {
                result = Some(mul(&[k_, sub(
                    &mul(&[div(&q(2), &sqrt_pi), st.clone()]),
                    &ml(&[q(2), a_.clone(), t.clone(), ee(&a_)]),
                )]));
            }
        } else if is(&p0.n, &q(-3)) && is(&p0.m, &half) && is(&p1.n, &neg(&half)) && is(&p1.m, &one) && p1.b.is_zero() {
            // 5.3.15: 1/(sqrt(a1*s)*(a0*sqrt(s) + b0)**3)
            let a_ = div(&p0.b, &p0.a);
            let k_ = div(&div(&c, &sqrt(&p1.a)), &pow(&p0.a, &q(3)));
            result = Some(mul(&[k_, sub(
                &mul(&[add(&[mul(&[q(2), sq(&a_), t.clone()]), one.clone()]), t.clone(), ee(&a_)]),
                &mul(&[div(&q(2), &sqrt_pi), a_.clone(), pow(t, &rat(3, 2))]),
            )]));
        } else if is(&p0.n, &q(-1)) && is(&p0.m, &one) && is(&p1.n, &neg(&half)) && is(&p1.m, &one) {
            // 5.3.23
            let a_ = div(&p0.b, &p0.a);
            let b_ = div(&p1.b, &p1.a);
            let k_ = div(&div(&c, &sqrt(&p1.a)), &p0.a);
            let d = sub(&b_, &a_);
            result = Some(mul(&[k_, div(&one, &sqrt(&d)), exp(&mul(&[neg(&a_), t.clone()])), erf(&mul(&[sqrt(&d), st.clone()]))]));
        }
    } else if p.len() == 2 && z.len() == 1 {
        let (p0, p1, z0) = (&p[0], &p[1], &z[0]);
        if is(&p0.n, &q(-1)) && is(&p0.m, &one) && is(&p1.n, &q(-1)) && is(&p1.m, &half)
            && is(&z0.n, &half) && is(&z0.m, &one) && z0.b.is_zero() {
            // 5.3.6
            let a_sq = div(&p1.b, &p1.a);
            let a_ = sq(&a_sq);
            let b_ = neg(&div(&p0.b, &p0.a));
            let k_ = mul(&[div(&div(&div(&sqrt(&z0.a), &p0.a), &p1.a), &sub(&a_, &b_)), c.clone()]);
            if pos(&a_sq) && pos(&b_) {
                result = Some(mul(&[k_, add(&[
                    mul(&[a_.clone(), exp(&mul(&[a_.clone(), t.clone()])), erfc(&mul(&[sqrt(&a_), st.clone()]))]),
                    // erf, where the reference has erfc
                    mul(&[sqrt(&a_), sqrt(&b_), exp(&mul(&[b_.clone(), t.clone()])), erf(&mul(&[sqrt(&b_), st.clone()]))]),
                    neg(&mul(&[b_.clone(), exp(&mul(&[b_.clone(), t.clone()]))])),
                ])]));
            }
        } else if is(&p0.n, &q(-1)) && is(&p0.m, &one) && p0.b.is_zero() && is(&p1.n, &q(-1)) && is(&p1.m, &half)
            && is(&z0.n, &one) && is(&z0.m, &half) {
            // 5.3.17
            let a_num = div(&z0.b, &z0.a);
            let a_ = div(&p1.b, &p1.a);
            if add(&[a_.clone(), a_num]).is_zero() {
                let k_ = mul(&[div(&div(&z0.a, &p0.a), &p1.a), c.clone()]);
                result = Some(mul(&[k_, sub(&mul(&[q(2), ee(&a_)]), &one)]));
            }
        } else if is(&p1.n, &q(-1)) && is(&p1.m, &one) && p1.b.is_zero() && is(&p0.n, &q(-2)) && is(&p0.m, &half)
            && is(&z0.n, &q(2)) && is(&z0.m, &half) {
            // 5.3.18
            let a_num = div(&z0.b, &z0.a);
            let a_ = div(&p0.b, &p0.a);
            if add(&[a_.clone(), a_num]).is_zero() {
                let k_ = mul(&[div(&div(&sq(&z0.a), &p1.a), &sq(&p0.a)), c.clone()]);
                result = Some(mul(&[k_, add(&[
                    one.clone(),
                    mul(&[q(8), sq(&a_), t.clone(), ee(&a_)]),
                    neg(&mul(&[div(&q(8), &sqrt_pi), a_.clone(), st.clone()])),
                ])]));
            }
        } else if is(&p1.n, &q(-1)) && is(&p1.m, &one) && p1.b.is_zero() && is(&p0.n, &q(-3)) && is(&p0.m, &half)
            && is(&z0.n, &q(3)) && is(&z0.m, &half) {
            // 5.3.19
            let a_num = div(&z0.b, &z0.a);
            let a_ = div(&p0.b, &p0.a);
            if add(&[a_.clone(), a_num]).is_zero() {
                let k_ = mul(&[div(&div(&pow(&z0.a, &q(3)), &p1.a), &pow(&p0.a, &q(3))), c.clone()]);
                let a2t = mul(&[sq(&a_), t.clone()]);
                result = Some(mul(&[k_, add(&[
                    ml(&[q(2), add(&[mul(&[q(8), pow(&a_, &q(4)), sq(t)]), mul(&[q(8), a2t.clone()]), one.clone()]), ee(&a_)]),
                    neg(&mul(&[div(&q(8), &sqrt_pi), a_.clone(), st.clone(), add(&[mul(&[q(2), a2t]), one.clone()])])),
                    neg(&one),
                ])]));
            }
        }
    } else if p.len() == 3 && z.is_empty() {
        let (p0, p1, p2) = (&p[0], &p[1], &p[2]);
        if is(&p0.n, &q(-1)) && p0.b.is_zero() && is(&p0.m, &one) && is(&p1.n, &q(-1)) && is(&p1.m, &one)
            && is(&p2.n, &neg(&half)) && is(&p2.m, &one) && p2.b.is_zero() {
            // 5.3.2
            let a_ = neg(&div(&p1.b, &p1.a));
            let k_ = mul(&[div(&div(&div(&one, &p0.a), &p1.a), &sqrt(&p2.a)), c.clone()]);
            if pos(&a_) {
                result = Some(mul(&[k_, sub(
                    &mul(&[pow(&a_, &rat(-3, 2)), exp(&mul(&[a_.clone(), t.clone()])), erf(&mul(&[sqrt(&a_), st.clone()]))]),
                    &mul(&[div(&div(&q(2), &a_), &sqrt_pi), st.clone()]),
                )]));
            }
        } else if is(&p0.n, &q(-1)) && is(&p0.m, &one) && is(&p1.n, &q(-1)) && is(&p1.m, &half)
            && is(&p2.n, &neg(&half)) && is(&p2.m, &one) && p2.b.is_zero() {
            // 5.3.9
            let a_sq = div(&p1.b, &p1.a);
            let a_ = sq(&a_sq);
            let b_ = neg(&div(&p0.b, &p0.a));
            let k_ = div(
                &div(&div(&div(&one, &p0.a), &p1.a), &sqrt(&p2.a)),
                &mul(&[sqrt(&b_), sub(&a_, &b_)]),
            );
            let k_ = mul(&[k_, c.clone()]);
            if pos(&a_sq) && pos(&b_) {
                result = Some(mul(&[k_, add(&[
                    mul(&[sqrt(&b_), exp(&mul(&[a_.clone(), t.clone()])), erfc(&mul(&[sqrt(&a_), st.clone()]))]),
                    mul(&[sqrt(&a_), exp(&mul(&[b_.clone(), t.clone()])), erf(&mul(&[sqrt(&b_), st.clone()]))]),
                    neg(&mul(&[sqrt(&b_), exp(&mul(&[b_.clone(), t.clone()]))])),
                ])]));
            }
        }
    }
    result.map(|r| (mul(&[heaviside(t), r]), tru()))
}
