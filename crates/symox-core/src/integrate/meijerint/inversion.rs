//! Laplace inversion by Meijer G-functions: the reference's
//! `meijerint_inversion` with `_rewrite_inversion`, `_int_inversion`,
//! `_inflate_fox_h` and `_check_antecedents_inversion` ([L] 5.7-5.10).
//!
//! `F(s)` is written as `fac*s**a*G(eta*s**b)`, the Bromwich integral of
//! each G-function is again a (Fox H, inflated to Meijer G) function of
//! `t`, and the result is returned under the convergence conditions of the
//! inversion integral. Exponential factors `exp(-c*s)` become time shifts.
//! A condition that simplifies to False makes the method give up.

use super::conds::my_unpolarify;
use super::{flip_g, get_coeff_exp, inflate_g, is_false, powdenest_polar, rewrite1, G};
use crate::expr::Constant;
use crate::functions as fx;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::relational::{eq, ge, gt, le, lt, ne};
use crate::special::SpecialFn;
use crate::{Expr, Func, Kind};

fn has_x(e: &Expr, x: &Expr) -> bool {
    crate::traverse::has(e, x)
}

fn b(v: bool) -> Expr {
    Expr::boolean(v)
}

fn and(v: Vec<Expr>) -> Expr {
    crate::logic::and(v)
}

fn or(v: Vec<Expr>) -> Expr {
    crate::logic::or(v)
}

/// `bool(e)` of a comparison; `None` where the reference raises TypeError.
fn truth(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Bool(v) => Some(*v),
        _ => None,
    }
}

/// `_is_analytic`: no `Heaviside` or `Abs` of the variable.
fn is_analytic(f: &Expr, x: &Expr) -> bool {
    let mut ok = true;
    crate::traverse::preorder(f, &mut |e| {
        if let Kind::Func(Func::Heaviside | Func::Abs, _) = e.kind() {
            if has_x(e, x) {
                ok = false;
            }
        }
    });
    ok
}

/// `exp(a)` in either native representation.
fn exp_arg(e: &Expr) -> Option<Expr> {
    match e.kind() {
        Kind::Func(Func::Exp, a) => Some(a[0].clone()),
        Kind::Pow(base, a) if base.is_constant(Constant::E) => Some(a.clone()),
        _ => None,
    }
}

/// `_rewrite_inversion`: absorb `po == x**s` into `g`.
fn rewrite_inversion(fac: &Expr, po: &Expr, g: &G, x: &Expr) -> Option<(Expr, G)> {
    let (_, s) = get_coeff_exp(po, x)?;
    let (a, bb) = get_coeff_exp(&g.z, x)?;
    let by = div(&s, &bb);
    let tr = |l: &[Expr]| l.iter().map(|t| add(&[t.clone(), by.clone()])).collect::<Vec<_>>();
    let c = powdenest_polar(&div(fac, &pow(&a, &by)));
    Some((c, G::new(tr(&g.an), tr(&g.ao), tr(&g.bm), tr(&g.bo), g.z.clone())?))
}

/// `_inflate_fox_h`: the function with integrand `d/Gamma(a*s)` as `C*G`.
fn inflate_fox_h(g: &G, a: &Expr) -> Option<(Expr, G)> {
    let n = a.as_number().filter(|n| n.is_rational())?;
    if n.is_negative() {
        return inflate_fox_h(&flip_g(g)?, &neg(a));
    }
    let (p, q) = n.numer_denom()?;
    let p: i64 = num_traits::ToPrimitive::to_i64(&p)?;
    let q: i64 = num_traits::ToPrimitive::to_i64(&q)?;
    let (d, g) = inflate_g(g, q)?;
    let pe = Expr::int(p);
    let d = div(
        &d,
        &mul(&[
            pow(&mul(&[Expr::int(2), Expr::pi()]), &Expr::rational(1 - p, 2)),
            pow(&pe, &Expr::rational(-1, 2)),
        ]),
    );
    let z = div(&g.z, &pow(&pe, &pe));
    let mut bo = g.bo.clone();
    bo.extend((0..p).map(|k| Expr::rational(k + 1, p)));
    Some((d, G::new(g.an.clone(), g.ao.clone(), g.bm.clone(), bo, z)?))
}

/// `_int_inversion`: the Laplace inversion integral of `g`.
fn int_inversion(g: &G, x: &Expr, t: &Expr) -> Option<Expr> {
    let (bb, a) = get_coeff_exp(&g.z, x)?;
    let h = G::new(g.an.clone(), g.ao.clone(), g.bm.clone(), g.bo.clone(), div(&bb, &pow(t, &a)))?;
    let (c, h) = inflate_fox_h(&h, &neg(&a))?;
    Some(mul(&[div(&c, t), h.to_expr()]))
}

/// `statement_half`: convergence of `z**a*exp(b*z**c)` on one side.
fn statement_half(a: &Expr, bb: &Expr, c: &Expr, z: &Expr, x: &Expr, plus: bool) -> Option<Expr> {
    let (coeff, exponent) = get_coeff_exp(z, x)?;
    let a = mul(&[a.clone(), exponent.clone()]);
    let bb = mul(&[bb.clone(), pow(&coeff, c)]);
    let c = mul(&[c.clone(), exponent]);
    let rc = fx::re(&c);
    let half = mul(&[Expr::i(), rc.clone(), Expr::pi(), Expr::rational(1, 2)]);
    let w = mul(&[bb.clone(), fx::exp(&if plus { half } else { neg(&half) })]);
    let zero = Expr::zero();
    let m1 = Expr::minus_one();
    Some(or(vec![
        and(vec![or(vec![eq(&bb, &zero), le(&rc, &zero)]), le(&fx::re(&a), &m1)]),
        and(vec![ne(&bb, &zero), eq(&fx::im(&c), &zero), gt(&rc, &zero), lt(&fx::re(&w), &zero)]),
        and(vec![
            ne(&bb, &zero),
            eq(&fx::im(&c), &zero),
            gt(&rc, &zero),
            le(&fx::re(&w), &zero),
            le(&fx::re(&a), &m1),
        ]),
    ]))
}

/// `statement`: convergence of `z**a*exp(b*z**c)` on both sides.
fn statement(a: &Expr, bb: &Expr, c: &Expr, z: &Expr, x: &Expr) -> Option<Expr> {
    Some(and(vec![statement_half(a, bb, c, z, x, true)?, statement_half(a, bb, c, z, x, false)?]))
}

/// `_check_antecedents_inversion`; `None` where the reference raises.
fn check_antecedents_inversion(g: &G, x: &Expr) -> Option<Expr> {
    let z = &g.z;
    let (_, e) = get_coeff_exp(z, x)?;
    if truth(&lt(&e, &Expr::zero()))? {
        return check_antecedents_inversion(&flip_g(g)?, x);
    }
    let (m, n, p, q) = (g.bm.len() as i64, g.an.len() as i64, g.p() as i64, g.q() as i64);
    let tau = m + n - p;
    let nu = q - m - n;
    let rho = Expr::rational(tau - nu, 2);
    let sigma = q - p;
    let epsilon = match sigma {
        1 => Expr::rational(1, 2),
        s if s > 1 => Expr::one(),
        _ => Expr::nan(),
    };
    let delta = g.delta();
    // is the computation valid for these parameters?
    let valid = truth(&ge(&delta, &div(&e, &Expr::int(2))))? || (p >= 1 && p >= q);
    if !valid {
        return Some(b(false));
    }
    // condition A
    for a in &g.an {
        for bm in &g.bm {
            if crate::assumptions::fact(&sub(a, bm), "integer") == Some(true) && truth(&gt(a, bm))? {
                return Some(b(false));
            }
        }
    }
    let e_of = |z: &Expr| -> Option<Expr> {
        let mut v = Vec::new();
        for a in &g.an {
            v.push(statement(&sub(a, &Expr::one()), &Expr::zero(), &Expr::zero(), z, x)?);
        }
        Some(and(v))
    };
    if p >= q {
        return e_of(z);
    }
    let pi = Expr::pi();
    let ipi = |k: i64| fx::exp(&mul(&[Expr::i(), pi.clone(), Expr::int(k)]));
    let theta = div(
        &add(&[Expr::rational(1 - sigma, 2), add(&g.bq()), neg(&add(&g.ap()))]),
        &Expr::int(sigma),
    );
    let ms = Expr::int(-sigma);
    let inv_sigma = Expr::rational(1, sigma);
    let h = |z: &Expr| statement(&theta, &ms, &inv_sigma, z, x);
    let hp = |z: &Expr| statement_half(&theta, &ms, &inv_sigma, z, x, true);
    let hm = |z: &Expr| statement_half(&theta, &ms, &inv_sigma, z, x, false);
    let int = Expr::int;
    let half_pi = div(&pi, &int(2));
    let zero = Expr::zero();
    let times_pi = |c: Expr| mul(&[c, pi.clone()]);
    let zm = |k: i64| mul(&[z.clone(), ipi(k)]);
    let mut conds = Vec::new();
    // Theorem 1
    conds.push(and(vec![
        b(1 <= n),
        b(1 <= m),
        ge(&sub(&times_pi(rho.clone()), &delta), &half_pi),
        gt(&delta, &zero),
        e_of(&zm(nu + 1))?,
    ]));
    // Theorem 2, statements (2) and (3)
    conds.push(and(vec![
        b(p + 1 <= m),
        b(m + 1 <= q),
        gt(&delta, &zero),
        lt(&delta, &half_pi),
        b(n == 0),
        ge(&sub(&times_pi(int(m - p + 1)), &delta), &half_pi),
        hp(&zm(q - m))?,
        hm(&zm(m - q))?,
    ]));
    // Theorem 2, statement (5)
    conds.push(and(vec![
        b(m == q),
        b(n == 0),
        gt(&delta, &zero),
        ge(&sub(&times_pi(add(&[int(sigma), epsilon.clone()])), &delta), &half_pi),
        h(z)?,
    ]));
    // Theorem 3, statements (6) and (7)
    conds.push(and(vec![
        or(vec![
            and(vec![b(p <= q - 2), b(1 <= tau), le(&int(tau), &Expr::rational(sigma, 2))]),
            and(vec![b(p + 1 <= m + n), le(&int(m + n), &Expr::rational(p + q, 2))]),
        ]),
        gt(&delta, &zero),
        lt(&delta, &half_pi),
        ge(&sub(&times_pi(int(tau + 1)), &delta), &half_pi),
        hp(&zm(nu))?,
        hm(&zm(-nu))?,
    ]));
    // Theorem 4, statements (10) and (11)
    conds.push(and(vec![
        b(1 <= m),
        gt(&rho, &zero),
        gt(&delta, &zero),
        lt(&add(&[delta.clone(), times_pi(rho.clone())]), &half_pi),
        ge(&sub(&times_pi(add(&[int(tau), epsilon])), &delta), &half_pi),
        hp(&zm(nu))?,
        hm(&zm(-nu))?,
    ]));
    // the trivial case
    conds.push(b(m == 0));
    Some(or(conds))
}

/// `meijerint_inversion(f, x, t)`: the inverse Laplace transform of `f(x)`
/// as `Piecewise((res, cond), (InverseLaplaceTransform(f, x, t, None),
/// True))`, a bare (possibly conditional) delta function, or `None`.
pub(crate) fn meijerint_inversion(f: &Expr, x: &Expr, t_: &Expr) -> Option<Expr> {
    let f_ = f.clone();
    let facts = crate::assumption_facts::FactSet::from_kwargs([("polar", true)]).ok()?;
    let t = Expr::from_kind(Kind::Sym(crate::symbol::Symbol::dummy_with_facts("t", facts, None).ok()?));
    let mut f = crate::subs::subs1(f, t_, &t);
    if !is_analytic(&f, x) {
        return None;
    }
    // exponential factors are time shifts
    let mut shift = Expr::zero();
    let args = if f.is_mul() {
        Some(crate::simptools::mul_args(&f))
    } else if exp_arg(&f).is_some() {
        Some(vec![f.clone()])
    } else {
        None
    };
    if let Some(mut args) = args {
        let mut newargs = Vec::new();
        let mut exponentials = Vec::new();
        while let Some(arg) = args.pop() {
            if let Some(ea) = exp_arg(&arg) {
                let arg2 = crate::expand::expand(&arg);
                if arg2.is_mul() {
                    args.extend(crate::simptools::mul_args(&arg2));
                    continue;
                }
                match get_coeff_exp(&ea, x) {
                    Some((a, e)) if e.is_one() => exponentials.push(a),
                    _ => newargs.push(arg),
                }
            } else if let Kind::Pow(base, ex) = arg.kind() {
                let arg2 = crate::expand::expand(&arg);
                if arg2.is_mul() {
                    args.extend(crate::simptools::mul_args(&arg2));
                    continue;
                }
                if !has_x(base, x) {
                    if let Some((a, e)) = get_coeff_exp(ex, x) {
                        if e.is_one() {
                            // the reference also keeps the factor in `f`,
                            // which would count the shift twice
                            exponentials.push(mul(&[a, fx::log(base)]));
                            continue;
                        }
                    }
                }
                newargs.push(arg);
            } else {
                newargs.push(arg);
            }
        }
        shift = add(&exponentials);
        f = mul(&newargs);
    }
    let shifted_t = add(&[t.clone(), shift.clone()]);
    if !has_x(&f, x) {
        let cond = eq(&fx::im(&shift), &Expr::zero());
        if is_false(&cond) {
            return None;
        }
        let res = mul(&[f, fx::dirac_delta(&shifted_t)]);
        return Some(crate::relational::piecewise(vec![(crate::subs::subs1(&res, &t, t_), cond)]));
    }
    let (fac, po, terms, mut cond) = rewrite1(&f, x)?;
    let mut res = Vec::new();
    for (c, s, g) in terms {
        let (c, g) = rewrite_inversion(&mul(&[fac.clone(), c]), &mul(&[po.clone(), pow(x, &s)]), &g, x)?;
        res.push(mul(&[c, int_inversion(&g, x, &t)?]));
        cond = and(vec![cond, check_antecedents_inversion(&g, x)?]);
        if is_false(&cond) {
            break;
        }
    }
    let cond = my_unpolarify(&cond);
    if is_false(&cond) {
        return None;
    }
    let mut res = my_unpolarify(&crate::special::hyperexpand(&add(&res)));
    let mut has_heaviside = false;
    crate::traverse::preorder(&res, &mut |e| has_heaviside |= matches!(e.kind(), Kind::Func(Func::Heaviside, _)));
    if !has_heaviside {
        res = mul(&[res, fx::heaviside(&t)]);
    }
    let res = crate::subs::subs1(&res, &t, &shifted_t);
    let cond = if matches!(cond.kind(), Kind::Bool(_)) { cond } else { crate::subs::subs1(&cond, &t, &shifted_t) };
    let held = crate::special::raw(
        SpecialFn::InverseLaplaceTransform,
        vec![crate::subs::subs1(&f_, &t, t_), x.clone(), t_.clone(), crate::transforms::none_plane()],
    );
    Some(crate::relational::piecewise(vec![(crate::subs::subs1(&res, &t, t_), cond), (held, b(true))]))
}
