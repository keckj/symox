//! Partial derivatives of the special functions.

use super::*;
use crate::functions::{cos, cosh, cot, exp, gamma, log, sin, sinh};
use crate::ops::sqrt;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};

fn one() -> Expr {
    Expr::one()
}

fn pg0(x: &Expr) -> Expr {
    super::polygamma(&Expr::zero(), x)
}

fn p1(x: &Expr) -> Expr {
    add(&[x.clone(), Expr::one()])
}

fn m1(x: &Expr) -> Expr {
    sub(x, &Expr::one())
}

/// `d f(args) / d args[i]`, `None` when SymPy leaves it unevaluated.
pub fn partial(s: SpecialFn, a: &[Expr], i: usize) -> Option<Expr> {
    use SpecialFn::*;
    if !s.nargs().contains(&a.len()) || i >= a.len() {
        return None;
    }
    let sp = sqrt(&Expr::pi());
    let e = || sf(s, a.to_vec());
    Some(match (s, i) {
        (Factorial, 0) => mul2(&gamma(&p1(&a[0])), &pg0(&p1(&a[0]))),
        (Binomial, 0) => {
            let (n, k) = (&a[0], &a[1]);
            mul2(&sub(&pg0(&p1(n)), &pg0(&p1(&sub(n, k)))), &e())
        }
        (Binomial, 1) => {
            let (n, k) = (&a[0], &a[1]);
            mul2(&sub(&pg0(&p1(&sub(n, k))), &pg0(&p1(k))), &e())
        }
        (Harmonic, 0) => {
            if a.len() == 1 {
                super::hurwitz_zeta(&Expr::int(2), &p1(&a[0]))
            } else {
                let m = &a[1];
                mul2(m, &super::hurwitz_zeta(&p1(m), &p1(&a[0])))
            }
        }
        (Catalan, 0) => {
            let x = &a[0];
            mul2(
                &add(&[pg0(&add(&[x.clone(), rat(1, 2)])), neg(&pg0(&add(&[x.clone(), int(2)]))), log(&int(4))]),
                &e(),
            )
        }
        (Polygamma, 1) => super::polygamma(&p1(&a[0]), &a[1]),
        (Digamma, 0) => super::polygamma(&one(), &a[0]),
        (Trigamma, 0) => super::polygamma(&int(2), &a[0]),
        (Beta, _) if a.len() == 2 => {
            let (x, y) = (&a[0], &a[1]);
            let xy = add(&[x.clone(), y.clone()]);
            let v = if i == 0 { x } else { y };
            mul2(&sub(&pg0(v), &pg0(&xy)), &e())
        }
        (LowerGamma, 1) => mul2(&pow(&a[1], &m1(&a[0])), &exp(&neg(&a[1]))),
        (UpperGamma, 1) => neg(&mul2(&pow(&a[1], &m1(&a[0])), &exp(&neg(&a[1])))),
        (Zeta, 1) => neg(&mul2(&a[0], &super::hurwitz_zeta(&p1(&a[0]), &a[1]))),
        (Polylog, 1) => div(&super::polylog(&m1(&a[0]), &a[1]), &a[1]),
        (Lerchphi, 0) => {
            let (z, s_, aa) = (&a[0], &a[1], &a[2]);
            div(&sub(&super::lerchphi(z, &m1(s_), aa), &mul2(aa, &e())), z)
        }
        (Lerchphi, 2) => neg(&mul2(&a[1], &super::lerchphi(&a[0], &p1(&a[1]), &a[2]))),
        (Erfi, 0) => div(&mul2(&int(2), &exp(&a[0].powi(2))), &sp),
        (Erf2, 0) => neg(&div(&mul2(&int(2), &exp(&neg(&a[0].powi(2)))), &sp)),
        (Erf2, 1) => div(&mul2(&int(2), &exp(&neg(&a[1].powi(2)))), &sp),
        (Erfinv, 0) => mul(&[sp, rat(1, 2), exp(&e().powi(2))]),
        (Erfcinv, 0) => mul(&[sp, rat(-1, 2), exp(&e().powi(2))]),
        (FresnelS, 0) => sin(&mul(&[Expr::pi(), a[0].powi(2), rat(1, 2)])),
        (FresnelC, 0) => cos(&mul(&[Expr::pi(), a[0].powi(2), rat(1, 2)])),
        (Ei, 0) => div(&exp(&a[0]), &a[0]),
        (Li | LiOffset, 0) => pow(&log(&a[0]), &Expr::minus_one()),
        (Si, 0) => div(&sin(&a[0]), &a[0]),
        (Ci, 0) => div(&cos(&a[0]), &a[0]),
        (Shi, 0) => div(&sinh(&a[0]), &a[0]),
        (Chi, 0) => div(&cosh(&a[0]), &a[0]),
        (Expint, 1) => neg(&super::expint(&m1(&a[0]), &a[1])),
        (Expint, 0) => {
            // SymPy: -z**(nu - 1)*meijerg([], [1, 1], [0, 0, 1 - nu], [], z)
            let (nu, z) = (&a[0], &a[1]);
            let t = super::tuple;
            let g = sf(
                MeijerG,
                vec![
                    t(vec![t(vec![]), t(vec![one(), one()])]),
                    t(vec![t(vec![Expr::zero(), Expr::zero(), sub(&one(), nu)]), t(vec![])]),
                    z.clone(),
                ],
            );
            neg(&mul2(&pow(z, &m1(nu)), &g))
        }
        (E1, 0) => neg(&super::expint(&Expr::zero(), &a[0])),
        (BesselJ | BesselY | BesselI | BesselK | Hankel1 | Hankel2, 1) => {
            let (nu, z) = (&a[0], &a[1]);
            let lo = sf(s, vec![m1(nu), z.clone()]);
            let hi = sf(s, vec![p1(nu), z.clone()]);
            let (ca, cb) = match s {
                BesselI => (-1, 1),
                BesselK => (1, -1),
                _ => (1, 1),
            };
            // b/2 f(nu-1) - a/2 f(nu+1)
            sub(&mul2(&rat(cb, 2), &lo), &mul2(&rat(ca, 2), &hi))
        }
        (Jn | Yn, 1) => {
            // f_{n-1}(z) - (n+1)/z f_n(z)
            let (n, z) = (&a[0], &a[1]);
            sub(&sf(s, vec![m1(n), z.clone()]), &mul2(&div(&p1(n), z), &e()))
        }
        (AiryAi, 0) => super::airyaiprime(&a[0]),
        (AiryBi, 0) => super::airybiprime(&a[0]),
        (AiryAiPrime, 0) => mul2(&a[0], &super::airyai(&a[0])),
        (AiryBiPrime, 0) => mul2(&a[0], &super::airybi(&a[0])),
        (EllipticK, 0) => {
            let m = &a[0];
            let omm = sub(&one(), m);
            div(
                &sub(&super::elliptic_e(m), &mul2(&omm, &super::elliptic_k(m))),
                &mul(&[int(2), m.clone(), omm]),
            )
        }
        (EllipticE, _) if a.len() == 1 => {
            let m = &a[0];
            div(&sub(&super::elliptic_e(m), &super::elliptic_k(m)), &mul2(&int(2), m))
        }
        (EllipticE, 0) => sqrt(&sub(&one(), &mul2(&a[1], &sin(&a[0]).powi(2)))),
        (EllipticE, 1) => {
            let (z, m) = (&a[0], &a[1]);
            div(&sub(&super::elliptic_e2(z, m), &super::elliptic_f(z, m)), &mul2(&int(2), m))
        }
        (EllipticF, 0) => pow(&sub(&one(), &mul2(&a[1], &sin(&a[0]).powi(2))), &rat(-1, 2)),
        (EllipticF, 1) => {
            let (z, m) = (&a[0], &a[1]);
            let omm = sub(&one(), m);
            let d = sqrt(&sub(&one(), &mul2(m, &sin(z).powi(2))));
            add(&[
                div(&super::elliptic_e2(z, m), &mul(&[int(2), m.clone(), omm.clone()])),
                neg(&div(&super::elliptic_f(z, m), &mul2(&int(2), m))),
                neg(&div(&sin(&mul2(&int(2), z)), &mul(&[int(4), omm, d]))),
            ])
        }
        (EllipticPi, 0) if a.len() == 2 => {
            let (n, m) = (&a[0], &a[1]);
            div(
                &add(&[
                    super::elliptic_e(m),
                    div(&mul2(&sub(m, n), &super::elliptic_k(m)), n),
                    div(&mul2(&sub(&n.powi(2), m), &e()), n),
                ]),
                &mul(&[int(2), sub(m, n), m1(n)]),
            )
        }
        (EllipticPi, 1) if a.len() == 2 => {
            let (n, m) = (&a[0], &a[1]);
            div(
                &add(&[div(&super::elliptic_e(m), &m1(m)), e()]),
                &mul2(&int(2), &sub(n, m)),
            )
        }
        (EllipticPi, 1) => {
            let (n, z, m) = (&a[0], &a[1], &a[2]);
            pow(
                &mul2(
                    &sub(&one(), &mul2(n, &sin(z).powi(2))),
                    &sqrt(&sub(&one(), &mul2(m, &sin(z).powi(2)))),
                ),
                &Expr::minus_one(),
            )
        }
        (LambertW, 0) => {
            let w = e();
            div(&w, &mul2(&a[0], &add(&[one(), w.clone()])))
        }
        (Hyper, 2) => {
            let ap = tuple_items(&a[0])?;
            let bq = tuple_items(&a[1])?;
            let fac = div(&mul(ap), &mul(bq));
            let nap: Vec<Expr> = ap.iter().map(p1).collect();
            let nbq: Vec<Expr> = bq.iter().map(p1).collect();
            mul2(&fac, &super::hyper(&nap, &nbq, &a[2]))
        }
        (Legendre, 1) => {
            let (n, x) = (&a[0], &a[1]);
            div(
                &mul2(n, &sub(&mul2(x, &e()), &super::legendre(&m1(n), x))),
                &m1(&x.powi(2)),
            )
        }
        (AssocLegendre, 2) => {
            let (n, m, x) = (&a[0], &a[1], &a[2]);
            div(
                &sub(
                    &mul(&[n.clone(), x.clone(), e()]),
                    &mul2(&add(&[m.clone(), n.clone()]), &super::assoc_legendre(&m1(n), m, x)),
                ),
                &m1(&x.powi(2)),
            )
        }
        (ChebyshevT, 1) => mul2(&a[0], &super::chebyshevu(&m1(&a[0]), &a[1])),
        (ChebyshevU, 1) => {
            let (n, x) = (&a[0], &a[1]);
            div(
                &sub(&mul2(&p1(n), &super::chebyshevt(&p1(n), x)), &mul2(x, &e())),
                &m1(&x.powi(2)),
            )
        }
        (Hermite, 1) => mul(&[int(2), a[0].clone(), super::hermite(&m1(&a[0]), &a[1])]),
        (Laguerre, 1) => neg(&super::assoc_laguerre(&m1(&a[0]), &one(), &a[1])),
        (AssocLaguerre, 2) => neg(&super::assoc_laguerre(&m1(&a[0]), &p1(&a[1]), &a[2])),
        (Jacobi, 3) => {
            let (n, al, be, x) = (&a[0], &a[1], &a[2], &a[3]);
            let c = div(&add(&[al.clone(), be.clone(), n.clone(), one()]), &int(2));
            mul2(&c, &super::jacobi(&m1(n), &p1(al), &p1(be), x))
        }
        (Gegenbauer, 2) => {
            let (n, al, x) = (&a[0], &a[1], &a[2]);
            mul(&[int(2), al.clone(), super::gegenbauer(&m1(n), &p1(al), x)])
        }
        (Ynm, 2) => {
            let (n, m, th, ph) = (&a[0], &a[1], &a[2], &a[3]);
            add(&[
                mul(&[m.clone(), cot(th), e()]),
                mul(&[
                    sqrt(&mul2(&sub(n, m), &add(&[m.clone(), n.clone(), one()]))),
                    exp(&neg(&mul2(&Expr::i(), ph))),
                    super::ynm(n, &p1(m), th, ph),
                ]),
            ])
        }
        (Ynm, 3) => mul(&[Expr::i(), a[1].clone(), e()]),
        (Sech, 0) => neg(&mul2(&crate::functions::tanh(&a[0]), &e())),
        (Csch, 0) => neg(&mul2(&crate::functions::coth(&a[0]), &e())),
        (Asec | Acsc, 0) => {
            let x = &a[0];
            let d = div(&one(), &mul2(&x.powi(2), &sqrt(&sub(&one(), &x.powi(-2)))));
            if s == Asec { d } else { neg(&d) }
        }
        (Asech, 0) => {
            let x = &a[0];
            neg(&div(&one(), &mul2(x, &sqrt(&sub(&one(), &x.powi(2))))))
        }
        (Acsch, 0) => {
            let x = &a[0];
            neg(&div(&one(), &mul2(&x.powi(2), &sqrt(&add(&[one(), x.powi(-2)])))))
        }
        (Acoth, 0) => div(&one(), &sub(&one(), &a[0].powi(2))),
        (Sinc, 0) => {
            let x = &a[0];
            sub(&div(&cos(x), x), &div(&sin(x), &x.powi(2)))
        }
        _ => return super::advanced::partial(s,a,i),
    })
}
