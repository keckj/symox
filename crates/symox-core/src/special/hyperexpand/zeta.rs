//! `lerchphi._eval_expand_func` and `polylog._eval_expand_func`, which the
//! Lerch phi origins of `hyperexpand` use for their basis.

use num_traits::{One, Signed, ToPrimitive};

use crate::functions::{self as f, func};
use crate::number::Number;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::special::{hurwitz_zeta, lerchphi, polylog};
use crate::{Expr, Func, Kind};

fn rat_of(e: &Expr) -> Option<num_rational::BigRational> {
    match e.kind() {
        Kind::Num(n) if n.is_rational() => n.to_ratio(),
        _ => None,
    }
}

fn small_int(e: &Expr) -> Option<i64> {
    match e.kind() {
        Kind::Num(Number::Int(i)) => Some(*i),
        _ => None,
    }
}

fn exp_polar(x: Expr) -> Expr { func(Func::ExpPolar, vec![x]) }

/// `polylog(s, z)._eval_expand_func()`.
pub(crate) fn polylog_expand(s: &Expr, z: &Expr) -> Expr {
    if s.is_one() {
        return neg(&f::log(&sub(&Expr::one(), z)));
    }
    if let Some(n) = small_int(s) {
        if n <= 0 {
            let u = Expr::dummy("u", Default::default());
            let mut start = div(&u, &sub(&Expr::one(), &u));
            for _ in 0..-n {
                start = mul(&[u.clone(), crate::diff::diff(&start, &u)]);
            }
            let e = crate::simptools::expand_mul(&start);
            return crate::subs::subs(&e, &[(u, z.clone())]);
        }
    }
    polylog(s, z)
}

/// `lerchphi(z, s, a)._eval_expand_func()`.
pub(crate) fn lerchphi_expand(z: &Expr, s: &Expr, a: &Expr) -> Expr {
    if z.is_one() {
        return hurwitz_zeta(s, a);
    }
    if let Some(n) = small_int(s) {
        if n <= 0 {
            let t = Expr::dummy("t", Default::default());
            let p = crate::poly::polyobj::Poly::from_exprs(&[pow(&add(&[t.clone(), a.clone()]), &Expr::int(-n))], std::slice::from_ref(&t));
            if let Ok(Ok(mut coeffs)) = p.map(|p| p[0].all_coeffs()) {
                coeffs.reverse();
                let mut start = div(&Expr::one(), &sub(&Expr::one(), &t));
                let mut res = Expr::zero();
                for c in coeffs {
                    res = add(&[res, mul(&[c, start.clone()])]);
                    start = mul(&[t.clone(), crate::diff::diff(&start, &t)]);
                }
                return crate::subs::subs(&res, &[(t, z.clone())]);
            }
        }
    }
    if let Some(ar) = rat_of(a) {
        let mut ar = ar;
        let mut addend = Expr::zero();
        let mut mult = Expr::one();
        let one = num_rational::BigRational::one();
        if ar > one {
            let mut n = ar.floor().to_integer();
            if num_rational::BigRational::from_integer(n.clone()) == ar {
                n -= 1;
            }
            let Some(nn) = n.to_i64().filter(|v| *v <= 100_000) else { return lerchphi(z, s, a) };
            ar -= num_rational::BigRational::from_integer(n);
            let a2 = Expr::num(Number::from_ratio(ar.clone()));
            mult = pow(z, &Expr::int(-nn));
            addend = add(&(0..nn)
                .map(|k| neg(&div(&pow(z, &Expr::int(k - nn)), &pow(&add(&[a2.clone(), Expr::int(k)]), s))))
                .collect::<Vec<_>>());
        } else if !ar.is_positive() {
            let n: num_bigint::BigInt = (-ar.clone()).floor().to_integer() + num_bigint::BigInt::one();
            let Some(nn) = n.to_i64().filter(|v| *v <= 100_000) else { return lerchphi(z, s, a) };
            ar += num_rational::BigRational::from_integer(n);
            let a2 = Expr::num(Number::from_ratio(ar.clone()));
            mult = pow(z, &Expr::int(nn));
            addend = add(&(0..nn)
                .map(|k| div(&pow(z, &Expr::int(nn - 1 - k)), &pow(&sub(&a2, &Expr::int(k + 1)), s)))
                .collect::<Vec<_>>());
        }
        let (m, n) = (ar.numer().clone(), ar.denom().clone());
        let Some(nq) = n.to_i64().filter(|v| *v <= 10_000) else { return lerchphi(z, s, a) };
        let m = Expr::num(Number::from_bigint(m));
        let nq_e = Expr::int(nq);
        let zet = exp_polar(div(&mul(&[Expr::int(2), Expr::pi(), Expr::i()]), &nq_e));
        let root = pow(z, &div(&Expr::one(), &nq_e));
        let up_zet = crate::special::polar::unpolarify(&zet, false);
        let mut addargs = Vec::new();
        for k in 0..nq {
            let ke = Expr::int(k);
            let mut p = polylog(s, &mul(&[pow(&zet, &ke), root.clone()]));
            if let Kind::Func(Func::Special(crate::special::SpecialFn::Polylog), pa) = p.kind() {
                p = polylog_expand(&pa[0].clone(), &pa[1].clone());
            }
            addargs.push(div(&p, &pow(&mul(&[pow(&up_zet, &ke), root.clone()]), &m)));
        }
        return add(&[addend, mul(&[mult, pow(&nq_e, &sub(s, &Expr::one())), add(&addargs)])]);
    }
    let pq = |p: i64, q: i64| Some((p, q));
    let exp_case = match z.kind() {
        Kind::Func(Func::Exp, ea) => rat_of(&div(&ea[0], &mul(&[Expr::pi(), Expr::i()]))).and_then(|_| {
            let r = rat_of(&div(&ea[0], &mul(&[Expr::int(2), Expr::pi(), Expr::i()])))?;
            Some((r.numer().to_i64()?, r.denom().to_i64()?))
        }),
        _ if z.is_minus_one() => pq(1, 2),
        _ if *z == Expr::i() => pq(1, 4),
        _ if *z == neg(&Expr::i()) => pq(-1, 4),
        _ => None,
    };
    if let Some((p, q)) = exp_case {
        if q.abs() <= 10_000 {
            let qe = Expr::int(q);
            let terms: Vec<Expr> = (0..q)
                .map(|k| {
                    mul(&[
                        f::exp(&div(&mul(&[Expr::int(2), Expr::pi(), Expr::i(), Expr::int(k), Expr::int(p)]), &qe)),
                        div(&Expr::one(), &pow(&qe, s)),
                        hurwitz_zeta(s, &div(&add(&[Expr::int(k), a.clone()]), &qe)),
                    ])
                })
                .collect();
            return add(&terms);
        }
    }
    lerchphi(z, s, a)
}
