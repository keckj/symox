//! Automatic evaluation of Bessel type functions.

use super::*;
use crate::functions::gamma;
use crate::ops::{div, mul, mul2, neg, pow};

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        BesselJ => eval_ji(&a[0], &a[1], true),
        BesselI => eval_ji(&a[0], &a[1], false),
        BesselY => eval_y(&a[0], &a[1]),
        BesselK => eval_k(&a[0], &a[1]),
        Jn => {
            let (nu, z) = (&a[0], &a[1]);
            if z.is_zero() {
                if nu.is_zero() {
                    return Some(Expr::one());
                }
                if is_integer(nu) == Some(true) {
                    if is_positive(nu) == Some(true) {
                        return Some(Expr::zero());
                    }
                    return Some(Expr::zoo());
                }
            }
            if is_inf(z) || is_neg_inf(z) {
                return Some(Expr::zero());
            }
            None
        }
        AiryAi | AiryBi | AiryAiPrime | AiryBiPrime => eval_airy(s, &a[0]),
        _ => None,
    }
}

fn re_positive(nu: &Expr) -> Option<bool> {
    match nu.kind() {
        Kind::Num(n) => Some(n.is_positive()),
        _ => {
            if crate::assumptions::is_real(nu) == Some(true) {
                is_positive(nu)
            } else {
                None
            }
        }
    }
}

fn re_negative(nu: &Expr) -> Option<bool> {
    match nu.kind() {
        Kind::Num(n) => Some(n.is_negative()),
        _ => {
            if crate::assumptions::is_real(nu) == Some(true) {
                is_negative(nu)
            } else {
                None
            }
        }
    }
}

/// Bessel J (`j = true`) and I.
fn eval_ji(nu: &Expr, z: &Expr, j: bool) -> Option<Expr> {
    let head = if j { SpecialFn::BesselJ } else { SpecialFn::BesselI };
    if z.is_zero() {
        if nu.is_zero() {
            return Some(Expr::one());
        }
        if (is_integer(nu) == Some(true) && is_zero(nu) == Some(false)) || re_positive(nu) == Some(true) {
            return Some(Expr::zero());
        }
        if re_negative(nu) == Some(true) && is_integer(nu) != Some(true) {
            return Some(Expr::zoo());
        }
        return None;
    }
    if j {
        if is_inf(z) || is_neg_inf(z) {
            return Some(Expr::zero());
        }
    } else {
        if is_inf(z) {
            return Some(Expr::oo());
        }
        if is_neg_inf(z) {
            return Some(mul2(&pow(&Expr::minus_one(), nu), &Expr::oo()));
        }
    }
    if could_extract_minus_sign(z) {
        let mz = neg(z);
        return Some(mul(&[pow(z, nu), pow(&mz, &neg(nu)), sf(head, vec![nu.clone(), mz])]));
    }
    if is_integer(nu) == Some(true) {
        if could_extract_minus_sign(nu) {
            let mn = neg(nu);
            return Some(if j {
                mul2(&pow(&Expr::minus_one(), &mn), &sf(head, vec![mn, z.clone()]))
            } else {
                sf(head, vec![mn, z.clone()])
            });
        }
        if let Some(t) = extract_i(z) {
            if !t.is_zero() {
                return Some(if j {
                    mul2(&pow(&Expr::i(), nu), &super::besseli(nu, &t))
                } else {
                    mul2(&pow(&Expr::i(), &neg(nu)), &super::besselj(nu, &neg(&t)))
                });
            }
        }
    }
    // Branch handling.
    if is_integer(nu) == Some(true) {
        let newz = super::polar::unpolarify(z, false);
        if newz != *z {
            return Some(sf(head, vec![nu.clone(), newz]));
        }
    } else {
        let (newz, n) = super::hyperexpand::extract_branch_factor(z, false);
        if !n.is_zero() {
            let f = crate::functions::exp(&mul(&[Expr::int(2), n, Expr::pi(), nu.clone(), Expr::i()]));
            return Some(mul2(&f, &sf(head, vec![nu.clone(), newz])));
        }
    }
    let nnu = super::polar::unpolarify(nu, false);
    if nnu != *nu {
        return Some(sf(head, vec![nnu, z.clone()]));
    }
    None
}

fn eval_y(nu: &Expr, z: &Expr) -> Option<Expr> {
    if z.is_zero() {
        if nu.is_zero() {
            return Some(Expr::neg_oo());
        }
        if nu.is_number() {
            return Some(Expr::zoo());
        }
        if crate::assumptions::is_real(nu) == Some(true) && is_zero(nu) == Some(false) {
            return Some(Expr::zoo());
        }
        return None;
    }
    if is_inf(z) || is_neg_inf(z) {
        return Some(Expr::zero());
    }
    // bessely(nu, +-I*oo) = exp(+-I*pi*(nu + 1)/2)*oo
    if let Some(t) = extract_i(z) {
        let sign = if is_inf(&t) { Some(1) } else if is_neg_inf(&t) { Some(-1) } else { None };
        if let Some(sign) = sign {
            let phase = mul(&[
                Expr::int(sign),
                Expr::i(),
                Expr::pi(),
                crate::ops::add(&[nu.clone(), Expr::one()]),
                rat(1, 2),
            ]);
            return Some(mul2(&crate::functions::exp(&phase), &Expr::oo()));
        }
    }
    if is_integer(nu) == Some(true) && could_extract_minus_sign(nu) {
        let mn = neg(nu);
        return Some(mul2(&pow(&Expr::minus_one(), &mn), &super::bessely(&mn, z)));
    }
    None
}

fn eval_k(nu: &Expr, z: &Expr) -> Option<Expr> {
    if z.is_zero() {
        if nu.is_zero() {
            return Some(Expr::oo());
        }
        if nu.is_number() {
            return Some(Expr::zoo());
        }
        if crate::assumptions::is_real(nu) == Some(true) && is_zero(nu) == Some(false) {
            return Some(Expr::zoo());
        }
        return None;
    }
    if is_inf(z) || extract_i(z).is_some_and(|t| is_inf(&t) || is_neg_inf(&t)) {
        return Some(Expr::zero());
    }
    if is_integer(nu) == Some(true) && could_extract_minus_sign(nu) {
        return Some(super::besselk(&neg(nu), z));
    }
    None
}

fn eval_airy(s: SpecialFn, z: &Expr) -> Option<Expr> {
    use SpecialFn::*;
    let three = Expr::int(3);
    if z.is_zero() {
        return Some(match s {
            AiryAi => div(&Expr::one(), &mul2(&pow(&three, &rat(2, 3)), &gamma(&rat(2, 3)))),
            AiryBi => div(&Expr::one(), &mul2(&pow(&three, &rat(1, 6)), &gamma(&rat(2, 3)))),
            AiryAiPrime => neg(&div(&Expr::one(), &mul2(&pow(&three, &rat(1, 3)), &gamma(&rat(1, 3))))),
            _ => div(&pow(&three, &rat(1, 6)), &gamma(&rat(1, 3))),
        });
    }
    if is_inf(z) {
        return Some(match s {
            AiryAi | AiryAiPrime => Expr::zero(),
            _ => Expr::oo(),
        });
    }
    if is_neg_inf(z) && matches!(s, AiryAi | AiryBi | AiryBiPrime) {
        return Some(Expr::zero());
    }
    None
}
