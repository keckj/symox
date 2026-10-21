//! The `HyperRep_*` representatives of hypergeometric functions and their
//! `nonrep` / `nonrepsmall` rewrites (SymPy's `sympy.functions.special.hyper`).
//!
//! A representative carries its argument on the Riemann surface of the
//! logarithm; the rewrites extract the branch factor `exp_polar(I*pi*k)` and
//! return the elementary form valid on that sheet. `nonrep` also covers the
//! exterior of the unit disc with a Piecewise, `nonrepsmall` gives the form
//! valid near the origin only.

use crate::functions as f;
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};
use crate::special::polar::unpolarify;
use crate::special::SpecialFn;
use crate::{Expr, Func, Kind};

use super::{half, int, ipi};

/// Whether `s` is one of the HyperRep classes.
pub(crate) fn is_hyperrep(s: SpecialFn) -> bool {
    use SpecialFn::*;
    matches!(s, HyperRepPower1 | HyperRepPower2 | HyperRepLog1 | HyperRepLog2 | HyperRepAtanh | HyperRepAsin1
        | HyperRepAsin2 | HyperRepSqrts1 | HyperRepSqrts2 | HyperRepCosasin | HyperRepSinasin)
}

/// `HyperRep.eval`: the unpolarified arguments when they differ.
pub(crate) fn hyperrep_eval(args: &[Expr]) -> Option<Vec<Expr>> {
    let n = args.len();
    let new: Vec<Expr> = args.iter().enumerate().map(|(i, a)| if i + 1 < n { unpolarify(a, false) } else { a.clone() }).collect();
    (new != args).then_some(new)
}

/// `Expr.extract_branch_factor(allow_half)`: split `exp_polar(I*pi*2*n)`
/// (or an odd multiple of `I*pi` when `allow_half`) off a polar number.
pub(crate) fn extract_branch_factor(e: &Expr, allow_half: bool) -> (Expr, Expr) {
    let mut n = Expr::zero();
    let mut res = Expr::one();
    let mut exps = Vec::new();
    for arg in crate::simptools::mul_args(e) {
        if let Kind::Func(Func::ExpPolar, a) = arg.kind() {
            exps.push(a[0].clone());
        } else {
            res = crate::ops::mul2(&res, &arg);
        }
    }
    let mut piimult = Expr::zero();
    let mut extras: Vec<Expr> = Vec::new();
    let ipi = ipi();
    while let Some(ex) = exps.pop() {
        if ex.is_add() {
            exps.extend(crate::simptools::add_args(&ex));
            continue;
        }
        if ex.is_mul() {
            if let Some(c) = crate::simptools::as_coefficient(&ex, &ipi) {
                piimult = add(&[piimult, c]);
                continue;
            }
        }
        extras.push(ex);
    }
    let (coeff, tail) = if crate::simptools::is_num(&piimult) {
        (piimult, Vec::new())
    } else {
        let syms: Vec<Expr> = crate::traverse::free_symbols(&piimult).into_iter().collect();
        crate::decompose::as_coeff_add(&piimult, &syms)
    };
    let branchfact = mul(&[f::ceiling(&sub(&div(&coeff, &int(2)), &half())), int(2)]);
    n = add(&[n, div(&branchfact, &int(2))]);
    let mut c = sub(&coeff, &branchfact);
    if allow_half {
        if let Some(nc) = crate::simptools::extract_additively(&c, &Expr::one()) {
            n = add(&[n, half()]);
            c = nc;
        }
    }
    let mut inner = vec![c];
    inner.extend(tail);
    let newexp = add(&[mul(&[ipi, add(&inner)]), add(&extras)]);
    if !newexp.is_zero() {
        res = crate::ops::mul2(&res, &f::func(Func::ExpPolar, vec![newexp]));
    }
    (res, n)
}

fn is_odd(n: &Expr) -> bool { crate::assumptions::fact(n, "odd") == Some(true) }
fn is_even(n: &Expr) -> bool { crate::assumptions::fact(n, "even") == Some(true) }
fn is_integer(a: &Expr) -> bool { crate::assumptions::is_integer(a) == Some(true) }

fn one() -> Expr { Expr::one() }
fn i() -> Expr { Expr::i() }
fn pi() -> Expr { Expr::pi() }
fn p1(x: &Expr) -> Expr { add(&[one(), x.clone()]) }
fn m1(x: &Expr) -> Expr { sub(&one(), x) }

/// The four branch forms of one representative. `p` holds the parameters
/// before the argument.
fn small(s: SpecialFn, p: &[Expr], x: &Expr) -> Expr {
    use SpecialFn::*;
    match s {
        HyperRepPower1 => pow(&m1(x), &p[0]),
        HyperRepPower2 => {
            let a = &p[0];
            mul(&[pow(&int(2), &sub(&mul(&[int(2), a.clone()]), &one())),
                pow(&p1(&sqrt(&m1(x))), &sub(&one(), &mul(&[int(2), a.clone()])))])
        }
        HyperRepLog1 => f::log(&m1(x)),
        HyperRepAtanh => div(&f::atanh(&sqrt(x)), &sqrt(x)),
        HyperRepAsin1 => div(&f::asin(&sqrt(x)), &sqrt(x)),
        HyperRepAsin2 => div(&small(HyperRepAsin1, &[], x), &small(HyperRepPower1, &[half()], x)),
        HyperRepSqrts1 => {
            let a2 = mul(&[int(2), p[0].clone()]);
            div(&add(&[pow(&m1(&sqrt(x)), &a2), pow(&p1(&sqrt(x)), &a2)]), &int(2))
        }
        HyperRepSqrts2 => {
            let a2 = mul(&[int(2), p[0].clone()]);
            div(&mul(&[sqrt(x), sub(&pow(&m1(&sqrt(x)), &a2), &pow(&p1(&sqrt(x)), &a2))]), &int(2))
        }
        HyperRepLog2 => f::log(&add(&[half(), div(&sqrt(&m1(x)), &int(2))])),
        HyperRepCosasin => f::cos(&mul(&[int(2), p[0].clone(), f::asin(&sqrt(x))])),
        HyperRepSinasin => mul(&[div(&sqrt(x), &sqrt(&m1(x))), f::sin(&mul(&[int(2), p[0].clone(), f::asin(&sqrt(x))]))]),
        _ => unreachable!(),
    }
}

fn small_minus(s: SpecialFn, p: &[Expr], x: &Expr) -> Expr {
    use SpecialFn::*;
    match s {
        HyperRepPower1 => pow(&p1(x), &p[0]),
        HyperRepPower2 => {
            let a = &p[0];
            mul(&[pow(&int(2), &sub(&mul(&[int(2), a.clone()]), &one())),
                pow(&p1(&sqrt(&p1(x))), &sub(&one(), &mul(&[int(2), a.clone()])))])
        }
        HyperRepLog1 => f::log(&p1(x)),
        HyperRepAtanh => div(&f::atan(&sqrt(x)), &sqrt(x)),
        HyperRepAsin1 => div(&f::asinh(&sqrt(x)), &sqrt(x)),
        HyperRepAsin2 => div(&small_minus(HyperRepAsin1, &[], x), &small_minus(HyperRepPower1, &[half()], x)),
        HyperRepSqrts1 => mul(&[pow(&p1(x), &p[0]), f::cos(&mul(&[int(2), p[0].clone(), f::atan(&sqrt(x))]))]),
        HyperRepSqrts2 => mul(&[sqrt(x), pow(&p1(x), &p[0]), f::sin(&mul(&[int(2), p[0].clone(), f::atan(&sqrt(x))]))]),
        HyperRepLog2 => f::log(&add(&[half(), div(&sqrt(&p1(x)), &int(2))])),
        HyperRepCosasin => f::cosh(&mul(&[int(2), p[0].clone(), f::asinh(&sqrt(x))])),
        HyperRepSinasin => mul(&[neg(&div(&sqrt(x), &sqrt(&p1(x)))), f::sinh(&mul(&[int(2), p[0].clone(), f::asinh(&sqrt(x))]))]),
        _ => unreachable!(),
    }
}

fn big(s: SpecialFn, p: &[Expr], x: &Expr, n: &Expr) -> Expr {
    use SpecialFn::*;
    let two = int(2);
    match s {
        HyperRepPower1 => {
            let a = &p[0];
            if is_integer(a) {
                return small(s, p, x);
            }
            mul(&[pow(&sub(x, &one()), a), f::exp(&mul(&[sub(&mul(&[two, n.clone()]), &one()), pi(), i(), a.clone()]))])
        }
        HyperRepPower2 => {
            let a = &p[0];
            let (sgn, n) = if is_odd(n) { (one(), sub(n, &one())) } else { (Expr::minus_one(), n.clone()) };
            mul(&[pow(&two, &sub(&mul(&[two.clone(), a.clone()]), &one())),
                pow(&add(&[one(), mul(&[sgn, i(), sqrt(&sub(x, &one()))])]), &sub(&one(), &mul(&[two.clone(), a.clone()]))),
                f::exp(&mul(&[int(-2), n, pi(), i(), a.clone()]))])
        }
        HyperRepLog1 => add(&[f::log(&sub(x, &one())), mul(&[sub(&mul(&[two, n.clone()]), &one()), pi(), i()])]),
        HyperRepAtanh => {
            let q = sqrt(x);
            let h = div(&mul(&[i(), pi()]), &two);
            if is_even(n) {
                div(&add(&[f::func(Func::Special(SpecialFn::Acoth), vec![q.clone()]), h]), &q)
            } else {
                div(&sub(&f::func(Func::Special(SpecialFn::Acoth), vec![q.clone()]), &h), &q)
            }
        }
        HyperRepAsin1 => {
            let q = sqrt(x);
            mul(&[pow(&Expr::minus_one(), n), add(&[div(&mul(&[sub(&half(), n), pi()]), &q), div(&mul(&[i(), f::acosh(&q)]), &q)])])
        }
        HyperRepAsin2 => div(&big(HyperRepAsin1, &[], x, n), &big(HyperRepPower1, &[half()], x, n)),
        HyperRepSqrts1 => {
            let a = &p[0];
            let q = sqrt(x);
            let e = |k: Expr| f::exp(&mul(&[two.clone(), pi(), i(), k, a.clone()]));
            let a2 = mul(&[two.clone(), a.clone()]);
            if is_even(n) {
                div(&add(&[mul(&[pow(&add(&[q.clone(), one()]), &a2), e(n.clone())]),
                    mul(&[pow(&sub(&q, &one()), &a2), e(sub(n, &one()))])]), &two)
            } else {
                let n = sub(n, &one());
                div(&add(&[mul(&[pow(&sub(&q, &one()), &a2), e(add(&[n.clone(), one()]))]),
                    mul(&[pow(&add(&[q.clone(), one()]), &a2), e(n)])]), &two)
            }
        }
        HyperRepSqrts2 => {
            let a = &p[0];
            let q = sqrt(x);
            let e = |k: Expr| f::exp(&mul(&[two.clone(), pi(), i(), a.clone(), k]));
            let a2 = mul(&[two.clone(), a.clone()]);
            if is_even(n) {
                mul(&[div(&q, &two), sub(&mul(&[pow(&sub(&q, &one()), &a2), e(sub(n, &one()))]),
                    &mul(&[pow(&add(&[q.clone(), one()]), &a2), e(n.clone())]))])
            } else {
                let n = sub(n, &one());
                mul(&[div(&q, &two), sub(&mul(&[pow(&sub(&q, &one()), &a2), e(add(&[n.clone(), one()]))]),
                    &mul(&[pow(&add(&[q.clone(), one()]), &a2), e(n)]))])
            }
        }
        HyperRepLog2 => {
            let q = sqrt(x);
            let head = add(&[mul(&[sub(n, &half()), pi(), i()]), f::log(&div(&q, &two))]);
            let t = mul(&[i(), f::asin(&div(&one(), &q))]);
            if is_even(n) { add(&[head, t]) } else { sub(&head, &t) }
        }
        HyperRepCosasin => {
            let a = &p[0];
            f::cosh(&add(&[mul(&[two.clone(), a.clone(), f::acosh(&sqrt(x))]),
                mul(&[a.clone(), pi(), i(), sub(&mul(&[two, n.clone()]), &one())])]))
        }
        HyperRepSinasin => {
            let a = &p[0];
            mul(&[neg(&div(&one(), &sqrt(&sub(&one(), &div(&one(), x))))),
                f::sinh(&add(&[mul(&[two.clone(), a.clone(), f::acosh(&sqrt(x))]),
                    mul(&[a.clone(), pi(), i(), sub(&mul(&[two, n.clone()]), &one())])]))])
        }
        _ => unreachable!(),
    }
}

fn big_minus(s: SpecialFn, p: &[Expr], x: &Expr, n: &Expr) -> Expr {
    use SpecialFn::*;
    let two = int(2);
    match s {
        HyperRepPower1 => {
            let a = &p[0];
            if is_integer(a) {
                return small_minus(s, p, x);
            }
            mul(&[pow(&p1(x), a), f::exp(&mul(&[two, n.clone(), pi(), i(), a.clone()]))])
        }
        HyperRepPower2 => {
            let a = &p[0];
            let sgn = if is_odd(n) { Expr::minus_one() } else { one() };
            mul(&[sgn.clone(), pow(&two, &sub(&mul(&[two.clone(), a.clone()]), &one())),
                pow(&add(&[sqrt(&p1(x)), sgn]), &sub(&one(), &mul(&[two.clone(), a.clone()]))),
                f::exp(&mul(&[int(-2), pi(), i(), a.clone(), n.clone()]))])
        }
        HyperRepLog1 => add(&[f::log(&p1(x)), mul(&[two, n.clone(), pi(), i()])]),
        HyperRepAtanh => {
            let q = sqrt(x);
            if is_even(n) { div(&f::atan(&q), &q) } else { div(&sub(&f::atan(&q), &pi()), &q) }
        }
        HyperRepAsin1 => {
            let q = sqrt(x);
            mul(&[pow(&Expr::minus_one(), n), add(&[div(&f::asinh(&q), &q), div(&mul(&[n.clone(), pi(), i()]), &q)])])
        }
        HyperRepAsin2 => div(&big_minus(HyperRepAsin1, &[], x, n), &big_minus(HyperRepPower1, &[half()], x, n)),
        HyperRepSqrts1 => {
            let a = &p[0];
            let q = sqrt(x);
            let head = mul(&[pow(&p1(x), a), f::exp(&mul(&[two.clone(), pi(), i(), n.clone(), a.clone()]))]);
            let arg = mul(&[two.clone(), a.clone(), f::atan(&q)]);
            if is_even(n) {
                mul(&[head, f::cos(&arg)])
            } else {
                mul(&[head, f::cos(&sub(&arg, &mul(&[two, pi(), a.clone()])))])
            }
        }
        HyperRepSqrts2 => {
            let a = &p[0];
            let q = sqrt(x);
            let head = mul(&[pow(&p1(x), a), f::exp(&mul(&[two.clone(), pi(), i(), n.clone(), a.clone()])), q.clone()]);
            let arg = mul(&[two.clone(), a.clone(), f::atan(&q)]);
            if is_even(n) {
                mul(&[head, f::sin(&arg)])
            } else {
                mul(&[head, f::sin(&sub(&arg, &mul(&[two, pi(), a.clone()])))])
            }
        }
        HyperRepLog2 => {
            let h = mul(&[pi(), i(), n.clone()]);
            if is_even(n) {
                add(&[h, f::log(&add(&[half(), div(&sqrt(&p1(x)), &two)]))])
            } else {
                add(&[h, f::log(&sub(&div(&sqrt(&p1(x)), &two), &half()))])
            }
        }
        HyperRepCosasin => {
            let a = &p[0];
            f::cosh(&add(&[mul(&[two.clone(), a.clone(), f::asinh(&sqrt(x))]), mul(&[two, a.clone(), pi(), i(), n.clone()])]))
        }
        HyperRepSinasin => {
            let a = &p[0];
            mul(&[neg(&div(&one(), &sqrt(&add(&[one(), div(&one(), x)])))),
                f::sinh(&add(&[mul(&[two.clone(), a.clone(), f::asinh(&sqrt(x))]), mul(&[two, a.clone(), pi(), i(), n.clone()])]))])
        }
        _ => unreachable!(),
    }
}

/// `_eval_rewrite_as_nonrep` (`small_only == false`) or
/// `_eval_rewrite_as_nonrepsmall` of one representative.
fn rewrite_one(s: SpecialFn, args: &[Expr], small_only: bool) -> Expr {
    let (params, z) = args.split_at(args.len() - 1);
    let (x, mut n) = extract_branch_factor(&z[0], true);
    let minus = !n.is_integer_number();
    if small_only {
        return if minus { small_minus(s, params, &x) } else { small(s, params, &x) };
    }
    if minus {
        n = sub(&n, &half());
    }
    let (sm, bg) = if minus {
        (small_minus(s, params, &x), big_minus(s, params, &x, &n))
    } else {
        (small(s, params, &x), big(s, params, &x, &n))
    };
    if bg == sm {
        return sm;
    }
    crate::relational::piecewise(vec![
        (bg, crate::relational::gt(&f::abs(&x), &one())),
        (sm, Expr::boolean(true)),
    ])
}

/// `expr.rewrite('nonrep')` / `expr.rewrite('nonrepsmall')`: rebuild
/// bottom-up, replacing every representative.
pub(crate) fn rewrite(e: &Expr, small_only: bool) -> Expr {
    if !has_hyperrep(e) {
        return e.clone();
    }
    crate::simptools::bottom_up(e, &mut |x| match x.kind() {
        Kind::Func(Func::Special(s), a) if is_hyperrep(*s) => rewrite_one(*s, a, small_only),
        _ => x.clone(),
    })
}

fn has_hyperrep(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Special(s), _) if is_hyperrep(*s) => true,
        _ => !e.is_atom() && e.args().iter().any(has_hyperrep),
    }
}

