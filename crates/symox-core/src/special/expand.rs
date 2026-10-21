//! `expand_func`, `besselsimp` and `hyperexpand`.

use super::*;
use crate::functions::{cos, cosh, erf, exp, gamma, log, sin, sinh};
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};
use crate::traverse::map_children;
use rustc_hash::FxHashMap;

fn rebuild(e: &Expr, f: &dyn Fn(&Expr) -> Option<Expr>) -> Expr {
    // Meijer expansions retain the original function as their domain fallback.
    // Expanding that fallback again would nest the same piecewise expression.
    if let Kind::Piecewise(branches)=e.kind() {
        let is_fallback=|value:&Expr,condition:&Expr|*condition==Expr::boolean(true)
            && matches!(value.kind(),Kind::Func(Func::Special(SpecialFn::MeijerG),_));
        if branches.iter().any(|(value,condition)|is_fallback(value,condition)) {
            let mapped=branches.iter().map(|(value,condition)|{
                (if is_fallback(value,condition){value.clone()}else{rebuild(value,f)},rebuild(condition,f))
            }).collect();
            return crate::relational::piecewise(mapped);
        }
    }
    let r = if e.is_atom() { e.clone() } else { map_children(e, &mut |c| rebuild(c, f)) };
    match f(&r) {
        Some(v) => v,
        None => r,
    }
}

/// Split `x + n` with integer `n` (returns `(x, n)`).
fn split_int_shift(e: &Expr) -> Option<(Expr, i64)> {
    if let Kind::Add(c, terms) = e.kind() {
        if let Some(n) = c.to_i64() {
            if c.is_integer() && !c.is_float() && n != 0 {
                return Some((crate::ops::add_from_terms(Number::ZERO, terms.clone()), n));
            }
        }
    }
    None
}

/// SymPy's `expand_func`: rewrite special functions in terms of simpler
/// ones where a closed form is known.
pub fn expand_func(e: &Expr) -> Expr {
    rebuild(e, &expand_one)
}

pub fn expand_func_deep(e: &Expr, deep: bool) -> Expr {
    if deep { expand_func(e) } else { expand_one(e).unwrap_or_else(|| e.clone()) }
}

fn expand_one(e: &Expr) -> Option<Expr> {
    use SpecialFn::*;
    if let Kind::Const(c) = e.kind() {
        // GoldenRatio / TribonacciConstant._eval_expand_func
        return match c {
            crate::expr::Constant::GoldenRatio => Some(div(&add(&[Expr::one(), sqrt(&Expr::int(5))]), &Expr::int(2))),
            crate::expr::Constant::TribonacciConstant => {
                let s33 = mul2(&Expr::int(3), &sqrt(&Expr::int(33)));
                let third = Expr::rational(1, 3);
                let a = pow(&sub(&Expr::int(19), &s33), &third);
                let b = pow(&add(&[Expr::int(19), s33]), &third);
                Some(div(&add(&[Expr::one(), a, b]), &Expr::int(3)))
            }
            _ => None,
        };
    }
    let Kind::Func(f, a) = e.kind() else { return None };
    match f {
        Func::Gamma => {
            let (x, n) = split_int_shift(&a[0])?;
            if n.abs() > 1000 {
                return None;
            }
            if n > 0 {
                let mut v: Vec<Expr> = (0..n).map(|k| add(&[x.clone(), Expr::int(k)])).collect();
                v.push(gamma(&x));
                Some(mul(&v))
            } else {
                let v: Vec<Expr> = (1..=(-n)).map(|k| sub(&x, &Expr::int(k))).collect();
                Some(div(&gamma(&x), &mul(&v)))
            }
        }
        Func::Special(s) => match s {
            BesselJ | BesselY | BesselI | BesselK if a.len() == 2 => {
                let v = bessel_expand(*s, &a[0], &a[1], &mut FxHashMap::default(), 0);
                (v != *e).then_some(v)
            }
            Binomial => {
                let k = as_small(&a[1])?;
                if !(0..=1000).contains(&k) {
                    return None;
                }
                let n = &a[0];
                let v: Vec<Expr> = (0..k).map(|i| sub(n, &Expr::int(i))).collect();
                Some(div(&mul(&v), &super::factorial(&Expr::int(k))))
            }
            Polygamma => {
                let m = as_small(&a[0])?;
                if m < 0 {
                    return None;
                }
                let (x, n) = split_int_shift(&a[1])?;
                if n <= 0 || n > 1000 {
                    return None;
                }
                // psi^(m)(x+n) = psi^(m)(x) + (-1)^m m! sum 1/(x+k)^(m+1)
                let mf = super::factorial(&Expr::int(m));
                let c = if m % 2 == 0 { mf } else { neg(&mf) };
                let v: Vec<Expr> = (0..n)
                    .map(|k| pow(&add(&[x.clone(), Expr::int(k)]), &Expr::int(-(m + 1))))
                    .collect();
                Some(add(&[super::polygamma(&a[0], &x), mul2(&c, &add(&v))]))
            }
            Harmonic if a.len() == 1 => {
                let (x, n) = split_int_shift(&a[0])?;
                if n <= 0 || n > 1000 {
                    return None;
                }
                let v: Vec<Expr> = (1..=n).map(|k| pow(&add(&[x.clone(), Expr::int(k)]), &Expr::minus_one())).collect();
                Some(add(&[super::harmonic(&x), add(&v)]))
            }
            Beta if a.len() == 2 => Some(div(
                &mul2(&gamma(&a[0]), &gamma(&a[1])),
                &gamma(&add(&[a[0].clone(), a[1].clone()])),
            )),
            Erf2 => Some(sub(&erf(&a[1]), &erf(&a[0]))),
            Polylog if a.len() == 2 => {
                let r = super::hyperexpand::polylog_expand(&a[0], &a[1]);
                (r != super::polylog(&a[0], &a[1])).then_some(r)
            }
            Lerchphi if a.len() == 3 => {
                let r = super::hyperexpand::lerchphi_expand(&a[0], &a[1], &a[2]);
                (r != super::lerchphi(&a[0], &a[1], &a[2])).then_some(r)
            }
            Expint => {
                let n = as_small(&a[0])?;
                if !(2..=200).contains(&n) {
                    return None;
                }
                let z = &a[1];
                let nf = super::factorial(&Expr::int(n - 1));
                let t1 = div(&mul2(&pow(&neg(z), &Expr::int(n - 1)), &super::expint(&Expr::one(), z)), &nf);
                let s: Vec<Expr> = (0..=(n - 2))
                    .map(|k| mul2(&super::factorial(&Expr::int(n - k - 2)), &pow(&neg(z), &Expr::int(k))))
                    .collect();
                Some(add(&[t1, div(&mul2(&exp(&neg(z)), &add(&s)), &nf)]))
            }
            Jn | Yn => spherical_closed(*s, &a[0], &a[1]),
            Hn1 | Hn2 => {
                let j=super::jn(&a[0],&a[1]);let y=super::yn(&a[0],&a[1]);
                Some(expand_func(&add(&[j,mul(&[Expr::i(),int(if *s==Hn1{1}else{-1}),y])])) )
            },
            Ynm => {
                let (n, m, th, ph) = (&a[0], &a[1], &a[2], &a[3]);
                let pre = sqrt(&mul(&[
                    add(&[mul2(&Expr::int(2), n), Expr::one()]),
                    pow(&mul2(&Expr::int(4), &Expr::pi()), &Expr::minus_one()),
                    super::factorial(&sub(n, m)),
                    pow(&super::factorial(&add(&[n.clone(), m.clone()])), &Expr::minus_one()),
                ]));
                let p = match (as_small(n), as_small(m)) {
                    (Some(nn), Some(mm)) if nn >= 0 && mm >= 0 && mm <= nn && nn <= 200 => {
                        // (-1)^m sin^m(theta) P_n^(m)(cos theta)
                        let c = super::orthopoly::legendre_deriv_coeffs(nn as usize, mm as usize);
                        let poly = super::combinatorial::poly_q(&c, &cos(th));
                        mul(&[Expr::int(if mm % 2 == 0 { 1 } else { -1 }), sin(th).powi(mm), poly])
                    }
                    _ => super::assoc_legendre(n, m, &cos(th)),
                };
                Some(mul(&[pre, exp(&mul(&[Expr::i(), m.clone(), ph.clone()])), p]))
            }
            Hyper => hyperexpand_one(a),
            MeijerG => super::advanced::meijer_expand(a),
            _ => None,
        },
        _ => None,
    }
}

/// Coefficient functions `f_n(z)` of the spherical Bessel functions.
fn fn_poly(n: i64, z: &Expr) -> Expr {
    // f_0 = 1/z, f_1 = 1/z^2, f_{n-1} + f_{n+1} = (2n+1)/z f_n
    let zi = pow(z, &Expr::minus_one());
    if n >= 0 {
        let (mut a, mut b) = (zi.clone(), z.powi(-2));
        if n == 0 {
            return a;
        }
        for k in 1..n {
            let c = crate::expand::expand(&sub(&mul(&[Expr::int(2 * k + 1), zi.clone(), b.clone()]), &a));
            a = b;
            b = c;
        }
        b
    } else {
        // downward: f_{n-1} = (2n+1)/z f_n - f_{n+1}
        let (mut b, mut a) = (zi.clone(), z.powi(-2)); // b = f_0, a = f_1
        let mut k = 0;
        while k > n {
            let c = crate::expand::expand(&sub(&mul(&[Expr::int(2 * k + 1), zi.clone(), b.clone()]), &a));
            a = b;
            b = c;
            k -= 1;
        }
        b
    }
}

fn spherical_closed(s: SpecialFn, n: &Expr, z: &Expr) -> Option<Expr> {
    let k = as_small(n)?;
    if !(0..=100).contains(&k) {
        return None;
    }
    let sgn = |p: i64| if p % 2 == 0 { Expr::one() } else { Expr::minus_one() };
    Some(if s == SpecialFn::Jn {
        // j_n = f_n sin z + (-1)^(n+1) f_{-n-1} cos z
        add(&[mul2(&fn_poly(k, z), &sin(z)), mul(&[sgn(k + 1), fn_poly(-k - 1, z), cos(z)])])
    } else {
        // y_n = (-1)^(n+1) f_{-n-1} sin z - f_n cos z
        sub(&mul(&[sgn(k + 1), fn_poly(-k - 1, z), sin(z)]), &mul2(&fn_poly(k, z), &cos(z)))
    })
}

/// A `BesselBase` node: the Bessel, Hankel and spherical Bessel functions.
pub(crate) fn is_bessel_base(e: &Expr) -> bool {
    use SpecialFn::*;
    matches!(e.kind(), Kind::Func(Func::Special(BesselJ | BesselY | BesselI | BesselK | Hankel1 | Hankel2 | Jn | Yn), _))
}

/// Whether `e` contains a `BesselBase` function.
pub(crate) fn has_bessel(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |x| {
        if !found && is_bessel_base(x) {
            found = true;
        }
    });
    found
}

/// The `(_a, _b)` constants of the Bessel recurrences.
fn bessel_ab(s: SpecialFn) -> (i64, i64) {
    match s {
        SpecialFn::BesselI => (-1, 1),
        SpecialFn::BesselK => (1, -1),
        _ => (1, 1),
    }
}

/// `BesselBase._eval_expand_func`: the recurrence down to the orders of
/// `(-1, 1]` (up from below `-1`), memoised by order.
fn bessel_expand(s: SpecialFn, nu: &Expr, z: &Expr, memo: &mut FxHashMap<Expr, Expr>, depth: usize) -> Expr {
    if let Some(v) = memo.get(nu) {
        return v.clone();
    }
    let this = super::apply(s, vec![nu.clone(), z.clone()]);
    if depth > 1000 || crate::assumptions::is_real(nu) != Some(true) {
        return this;
    }
    let (a, b) = bessel_ab(s);
    let (a, b) = (Expr::int(a), Expr::int(b));
    let one = Expr::one();
    let sub_expand = |k: Expr, memo: &mut FxHashMap<Expr, Expr>| {
        let f = super::apply(s, vec![k.clone(), z.clone()]);
        // the reference expands the evaluated node, which may no longer be
        // this function of `z`
        match f.kind() {
            Kind::Func(Func::Special(h), args) if *h == s && args[1] == *z => bessel_expand(s, &args[0], z, memo, depth + 1),
            _ => f,
        }
    };
    let r = if crate::assumptions::is_positive(&sub(nu, &one)) == Some(true) {
        let m2 = sub_expand(sub(nu, &Expr::int(2)), memo);
        let m1 = sub_expand(sub(nu, &one), memo);
        add(&[
            mul(&[neg(&a), b.clone(), m2]),
            div(&mul(&[Expr::int(2), a.clone(), sub(nu, &one), m1]), z),
        ])
    } else if crate::assumptions::is_negative(&add(&[nu.clone(), one.clone()])) == Some(true) {
        let p1 = sub_expand(add(&[nu.clone(), one.clone()]), memo);
        let p2 = sub_expand(add(&[nu.clone(), Expr::int(2)]), memo);
        sub(&div(&mul(&[Expr::int(2), b.clone(), add(&[nu.clone(), one.clone()]), p1]), z), &mul(&[a, b, p2]))
    } else {
        this
    };
    memo.insert(nu.clone(), r.clone());
    r
}

/// `Mul.make_args(z)` meets `factors`.
fn has_factor(z: &Expr, factors: &[Expr]) -> bool {
    crate::simptools::mul_args(z).iter().any(|f| factors.contains(f))
}

fn bessel_parts(e: &Expr, heads: &[SpecialFn]) -> Option<(SpecialFn, Expr, Expr)> {
    match e.kind() {
        Kind::Func(Func::Special(s), a) if heads.contains(s) && a.len() == 2 => Some((*s, a[0].clone(), a[1].clone())),
        _ => None,
    }
}

fn ipi(c: Expr) -> Expr {
    mul(&[c, Expr::i(), Expr::pi()])
}

fn exp_polar(c: Expr) -> Expr {
    crate::functions::func(Func::ExpPolar, vec![ipi(c)])
}

/// `_bessel_simp_recursion` for one head: `f(a)` is replaced by
/// `2*(a + 1)/x*f(a + 1) - f(a + 2)` while both shifted terms occur.
fn bessel_recursion(head: SpecialFn, mut expr: Expr) -> Expr {
    for _ in 0..1000 {
        let mut found: Vec<Expr> = Vec::new();
        crate::traverse::preorder(&expr, &mut |x| {
            if bessel_parts(x, &[head]).is_some() && !found.contains(x) {
                found.push(x.clone());
            }
        });
        // sorted by `re(order)`; an undecidable comparison ends the pass
        let keys: Vec<Expr> = found.iter().map(|b| crate::functions::re(&b.args()[0])).collect();
        let mut idx: Vec<usize> = (0..found.len()).collect();
        let mut failed = false;
        // insertion sort with the reference's `<`
        for i in 1..idx.len() {
            let mut j = i;
            while j > 0 {
                let lt = crate::relational::lt(&keys[idx[j]], &keys[idx[j - 1]]);
                match lt.kind() {
                    Kind::Bool(true) => {
                        idx.swap(j, j - 1);
                        j -= 1;
                    }
                    Kind::Bool(false) => break,
                    _ => {
                        failed = true;
                        break;
                    }
                }
            }
            if failed {
                break;
            }
        }
        if failed {
            return expr;
        }
        let mut changed = false;
        for &i in &idx {
            let ba = &found[i];
            let (a, x) = (ba.args()[0].clone(), ba.args()[1].clone());
            let a1 = add(&[a.clone(), Expr::one()]);
            let bap1 = super::apply(head, vec![a1.clone(), x.clone()]);
            let bap2 = super::apply(head, vec![add(&[a.clone(), Expr::int(2)]), x.clone()]);
            if crate::traverse::has(&expr, &bap1) && crate::traverse::has(&expr, &bap2) {
                let rep = sub(&div(&mul(&[Expr::int(2), a1, bap1]), &x), &bap2);
                expr = crate::subs::subs1(&expr, ba, &rep);
                changed = true;
                break;
            }
        }
        if !changed {
            return expr;
        }
    }
    expr
}

/// SymPy's `besselsimp`: factors `I` and `-1` (and their polar forms) are
/// taken out of the arguments of `besselj` and `besseli`, half integer
/// orders are written with elementary functions, integer orders above one
/// are reduced by the recurrence, the `besselj`/`bessely` recursion is
/// applied and a changed result is factored.
pub fn besselsimp(e: &Expr) -> Expr {
    use SpecialFn::{BesselI, BesselJ, BesselK, BesselY};
    let orig = e.clone();
    let half = Expr::rational(1, 2);
    let ifactors = [Expr::i(), exp_polar(half.clone()), exp_polar(neg(&half))];
    let minusfactors = [Expr::minus_one(), exp_polar(Expr::one())];
    // besselj(nu, z).rewrite(besseli) and besseli(nu, z).rewrite(besselj)
    let to_other = |head: SpecialFn| {
        let ifactors = ifactors.clone();
        move |x: &Expr| -> Option<Expr> {
            let (_, nu, z) = bessel_parts(x, &[head])?;
            if !has_factor(&z, &ifactors) {
                return None;
            }
            let (other, sgn) = if head == BesselJ { (BesselI, 1) } else { (BesselJ, -1) };
            let lift = super::apply(SpecialFn::PolarLift, vec![mul2(&Expr::int(-sgn), &Expr::i())]);
            Some(mul2(
                &crate::functions::exp(&ipi(div(&mul2(&Expr::int(sgn), &nu), &Expr::int(2)))),
                &super::apply(other, vec![nu, mul2(&lift, &z)]),
            ))
        }
    };
    let tominus = |head: SpecialFn| {
        let minusfactors = minusfactors.clone();
        move |x: &Expr| -> Option<Expr> {
            let (_, nu, z) = bessel_parts(x, &[head])?;
            if !has_factor(&z, &minusfactors) {
                return None;
            }
            Some(mul2(
                &crate::functions::exp(&ipi(nu.clone())),
                &super::apply(head, vec![nu, mul2(&exp_polar(Expr::minus_one()), &z)]),
            ))
        }
    };
    let mut expr = rebuild(e, &to_other(BesselJ));
    expr = rebuild(&expr, &to_other(BesselI));
    expr = rebuild(&expr, &tominus(BesselJ));
    expr = rebuild(&expr, &tominus(BesselI));
    for head in [BesselJ, BesselY, BesselI, BesselK] {
        expr = rebuild(&expr, &|x| {
            let (s, nu, z) = bessel_parts(x, &[head])?;
            let n = as_exact(&nu)?;
            if n.denom() == Number::int(2) {
                return half_order(s, n, &z);
            }
            if n.is_integer() && n.to_i64().is_some_and(|k| k > 1) {
                let v = bessel_expand(s, &nu, &z, &mut FxHashMap::default(), 0);
                return Some(crate::expand::expand(&v));
            }
            None
        });
    }
    for head in [BesselJ, BesselY] {
        let present = bessel_head_present(&expr, head);
        if present {
            expr = bessel_recursion(head, expr);
        }
    }
    if expr != orig {
        expr = crate::poly::api::factor(&expr).unwrap_or(expr);
    }
    expr
}

fn bessel_head_present(e: &Expr, head: SpecialFn) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |x| {
        if !found && bessel_parts(x, &[head]).is_some() {
            found = true;
        }
    });
    found
}

fn half_order(s: SpecialFn, nu: &Number, z: &Expr) -> Option<Expr> {
    use SpecialFn::*;
    let two = Expr::int(2);
    let pre = mul(&[sqrt(&two), pow(&Expr::pi(), &rat(-1, 2)), pow(z, &rat(-1, 2))]);
    let half = Number::rational_i64(1, 2).unwrap();
    let mhalf = half.neg();
    if *nu == half || *nu == mhalf {
        let up = *nu == half;
        return Some(match s {
            BesselJ => mul2(&pre, &if up { sin(z) } else { cos(z) }),
            BesselY => {
                if up {
                    neg(&mul2(&pre, &cos(z)))
                } else {
                    mul2(&pre, &sin(z))
                }
            }
            BesselI => mul2(&pre, &if up { sinh(z) } else { cosh(z) }),
            _ => mul(&[sqrt(&Expr::pi()), pow(&two, &rat(-1, 2)), pow(z, &rat(-1, 2)), exp(&neg(z))]),
        });
    }
    // general half integer: via spherical functions, written SymPy style as
    // sqrt(2) * P(z) / (sqrt(pi) * z**(nu)) with P free of fractional powers
    let n = nu.sub(&half).to_i64()?;
    if n.abs() > 50 {
        return None;
    }
    let pre2 = mul(&[sqrt(&two), pow(&Expr::pi(), &rat(-1, 2)), sqrt(z)]);
    let normalized = |sign: i64, sph: Expr, k: i64| -> Expr {
        // sqrt(2z/pi) * sign * sph  =  sqrt(2)/sqrt(pi) * z^(-k-1/2) * (sign z^(k+1) sph)
        let p = crate::expand::expand(&mul(&[Expr::int(sign), z.powi(k + 1), sph]));
        let first_neg = crate::printer::ordered_terms(&p)
            .first()
            .is_some_and(|t| crate::ops::as_coeff_mul(t).0.is_negative());
        let (p, sg) = if first_neg { (neg(&p), -1) } else { (p, 1) };
        mul(&[Expr::int(sg), sqrt(&two), pow(&Expr::pi(), &rat(-1, 2)), pow(z, &rat(-2 * k - 1, 2)), p])
    };
    match s {
        BesselJ | BesselY => {
            // J_{n+1/2} = sqrt(2z/pi) j_n(z), Y_{n+1/2} = sqrt(2z/pi) y_n(z),
            // J_{-m-1/2} = (-1)^(m+1) Y_{m+1/2}, Y_{-m-1/2} = (-1)^m J_{m+1/2}
            let (head, k, sign) = if n >= 0 {
                (if s == BesselJ { Jn } else { Yn }, n, 1)
            } else {
                let m = -n - 1;
                if s == BesselJ {
                    (Yn, m, if m % 2 == 0 { -1 } else { 1 })
                } else {
                    (Jn, m, if m % 2 == 0 { 1 } else { -1 })
                }
            };
            let sph = spherical_closed(head, &Expr::int(k), z)?;
            let _ = &pre2;
            Some(normalized(sign, sph, k))
        }
        BesselI => {
            // I_{n+1/2}(z) = sqrt(2z/pi) i_n(z), i_n(z) = (-i)^n j_n(iz); I_{-n-1/2} via y_n
            let (head, k) = if n >= 0 { (Jn, n) } else { (Yn, -n - 1) };
            let iz = mul2(&Expr::i(), z);
            let sph = spherical_closed(head, &Expr::int(k), &iz)?;
            let f = if n >= 0 {
                pow(&neg(&Expr::i()), &Expr::int(k))
            } else {
                // I_{-k-1/2}(z) = sqrt(2z/pi) (-i)^(k+1) y_k(iz)
                pow(&neg(&Expr::i()), &Expr::int(k + 1))
            };
            Some(normalized(1, crate::expand::expand(&mul2(&f, &sph)), k))
        }
        BesselK => {
            // K_{n+1/2}(z) = sqrt(pi/(2z)) e^-z sum_k (n+k)!/(k!(n-k)!(2z)^k), K_{-nu} = K_nu
            let k = if n >= 0 { n } else { -n - 1 };
            let v: Vec<Expr> = (0..=k)
                .map(|j| {
                    let c = super::integers::factorial((k + j) as u64)
                        / (super::integers::factorial(j as u64) * super::integers::factorial((k - j) as u64));
                    let c = Expr::num(Number::rational(c, num_bigint::BigInt::from(1) << (j as usize)).unwrap());
                    mul2(&c, &z.powi(k - j))
                })
                .collect();
            Some(mul(&[
                sqrt(&Expr::pi()),
                pow(&two, &rat(-1, 2)),
                pow(z, &rat(-2 * k - 1, 2)),
                exp(&neg(z)),
                add(&v),
            ]))
        }
        _ => None,
    }
}

/// SymPy's `hyperexpand(e)` with the default keywords; an expression whose
/// expansion raises in the reference is returned unchanged.
pub fn hyperexpand(e: &Expr) -> Expr {
    super::hyperexpand::hyperexpand(e)
}

#[allow(dead_code)]
fn hyperexpand_one(a: &[Expr]) -> Option<Expr> {
    let ap = tuple_items(&a[0])?;
    let bq = tuple_items(&a[1])?;
    let z = &a[2];
    let one = Expr::one();
    match (ap.len(), bq.len()) {
        (0, 0) => Some(exp(z)),
        (1, 0) => Some(pow(&sub(&one, z), &neg(&ap[0]))),
        (1, 1) if ap[0].is_one() && bq[0] == Expr::int(2) => {
            Some(sub(&div(&exp(z), z), &pow(z, &Expr::minus_one())))
        }
        (2, 1) if ap[0].is_one() && ap[1].is_one() && bq[0] == Expr::int(2) => {
            Some(neg(&div(&log(&sub(&one, z)), z)))
        }
        (2, 1) if ap[0] == rat(1, 2) && ap[1] == rat(1, 2) && bq[0] == rat(3, 2) => {
            let s = sqrt(z);
            Some(div(&crate::functions::asin(&s), &s))
        }
        _ => None,
    }
}
