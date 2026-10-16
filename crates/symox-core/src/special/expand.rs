//! `expand_func`, `besselsimp` and `hyperexpand`.

use super::*;
use crate::functions::{cos, cosh, erf, exp, gamma, log, sin, sinh};
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};
use crate::traverse::map_children;

fn rebuild(e: &Expr, f: &dyn Fn(&Expr) -> Option<Expr>) -> Expr {
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

fn expand_one(e: &Expr) -> Option<Expr> {
    use SpecialFn::*;
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
            Polylog => {
                let s = as_small(&a[0])?;
                let z = &a[1];
                if s == 1 {
                    return Some(neg(&log(&sub(&Expr::one(), z))));
                }
                if s <= 0 && s > -200 {
                    let n = (-s) as u64;
                    let u = div(z, &sub(&Expr::one(), z));
                    let v: Vec<Expr> = (1..=n + 1)
                        .map(|k| {
                            let c = super::integers::factorial(k - 1) * super::integers::stirling2(n + 1, k);
                            mul2(&Expr::num(Number::from_bigint(c)), &u.powi(k as i64))
                        })
                        .collect();
                    return Some(add(&v));
                }
                None
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

/// Rewrite Bessel functions of half integer order in elementary terms.
pub fn besselsimp(e: &Expr) -> Expr {
    rebuild(e, &|x| {
        let Kind::Func(Func::Special(s), a) = x.kind() else { return None };
        if !matches!(s, SpecialFn::BesselJ | SpecialFn::BesselY | SpecialFn::BesselI | SpecialFn::BesselK) {
            return None;
        }
        let nu = as_exact(&a[0])?;
        if nu.denom() != Number::int(2) {
            return None;
        }
        half_order(*s, nu, &a[1])
    })
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

/// Closed forms of a few generalised hypergeometric functions.
pub fn hyperexpand(e: &Expr) -> Expr {
    rebuild(e, &|x| match x.kind() {
        Kind::Func(Func::Special(SpecialFn::Hyper), a) => hyperexpand_one(a),
        _ => None,
    })
}

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
