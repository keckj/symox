//! Double precision (real and complex) evaluation of the special functions.

use num_complex::Complex64 as C64;

use super::{tuple_items, SpecialFn};
use crate::expr::{Expr, Func, Kind};

pub mod bessel;
pub mod dd;
pub mod erf;
pub mod gamma;
pub mod misc;

pub(crate) const PI: f64 = std::f64::consts::PI;

#[inline]
pub(crate) fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

#[inline]
pub(crate) fn cr(re: f64) -> C64 {
    C64::new(re, 0.0)
}

#[inline]
pub(crate) fn is_real_c(z: C64) -> bool {
    z.im == 0.0
}

/// `sin(pi x)` with exact argument reduction.
pub fn sinpi(x: f64) -> f64 {
    if x == x.floor() {
        return 0.0;
    }
    let r = x - 2.0 * (x / 2.0).floor(); // [0, 2)
    let (r, sign) = if r >= 1.0 { (r - 1.0, -1.0) } else { (r, 1.0) };
    let r = if r > 0.5 { 1.0 - r } else { r };
    sign * (PI * r).sin()
}

/// `cos(pi x)` with exact argument reduction.
pub fn cospi(x: f64) -> f64 {
    let r = (x - 2.0 * (x / 2.0).floor()).abs(); // [0, 2)
    let r = if r > 1.0 { 2.0 - r } else { r }; // [0, 1]
    if r == 0.5 {
        return 0.0;
    }
    if r > 0.5 {
        -(PI * (1.0 - r)).cos()
    } else {
        (PI * r).cos()
    }
}

/// Real double result (for `lambdify`) of a one-argument function.
pub fn real1(s: SpecialFn, x: f64) -> f64 {
    apply(s, &[cr(x)]).filter(|z| z.im == 0.0).map(|z| z.re).unwrap_or(f64::NAN)
}

/// Real double result of a two-argument function.
pub fn real2(s: SpecialFn, x: f64, y: f64) -> f64 {
    apply(s, &[cr(x), cr(y)]).filter(|z| z.im == 0.0).map(|z| z.re).unwrap_or(f64::NAN)
}

/// Real double result of an n-argument function.
pub fn realn(s: SpecialFn, a: &[f64]) -> f64 {
    let v: Vec<C64> = a.iter().map(|x| cr(*x)).collect();
    apply(s, &v).filter(|z| z.im == 0.0).map(|z| z.re).unwrap_or(f64::NAN)
}

/// Evaluate `s(args)` with numeric expression arguments to a complex double.
pub fn eval_c64(s: SpecialFn, args: &[Expr]) -> Option<C64> {
    if s == SpecialFn::Hyper {
        let ap: Option<Vec<C64>> = tuple_items(&args[0])?.iter().map(crate::evalf::eval_c64).collect();
        let bq: Option<Vec<C64>> = tuple_items(&args[1])?.iter().map(crate::evalf::eval_c64).collect();
        let z = crate::evalf::eval_c64(&args[2])?;
        return misc::hyper_c(&ap?, &bq?, z);
    }
    if s == SpecialFn::Bell {
        return None;
    }
    let a: Option<Vec<C64>> = args.iter().map(crate::evalf::eval_c64).collect();
    apply(s, &a?)
}

/// Numerical evaluation as an expression (`None` if not possible or not
/// finite).
pub fn eval_exprs(s: SpecialFn, args: &[Expr]) -> Option<Expr> {
    let z = eval_c64(s, args)?;
    if !z.re.is_finite() || !z.im.is_finite() {
        return None;
    }
    Some(crate::evalf::c64_to_expr(z))
}

fn int_arg(z: C64) -> Option<i64> {
    if z.im == 0.0 && z.re == z.re.floor() && z.re.abs() < 1e15 {
        Some(z.re as i64)
    } else {
        None
    }
}

/// Core of the numeric dispatch.
pub fn apply(s: SpecialFn, a: &[C64]) -> Option<C64> {
    use SpecialFn::*;
    let x = *a.first()?;
    let r = match s {
        Factorial => gamma::factorial_c(x)?,
        Factorial2 => {
            // x!! = 2^{x/2} (2/pi)^{(1 - cos(pi x))/4} Gamma(x/2 + 1)
            let n = int_arg(x)?;
            if n < -1 {
                return None;
            }
            let mut r = 1.0;
            let mut k = n;
            while k > 1 {
                r *= k as f64;
                k -= 2;
            }
            cr(r)
        }
        Subfactorial => {
            let n = int_arg(x)?;
            if n < 0 {
                return None;
            }
            let v = super::integers::subfactorial(n as u64);
            cr(num_traits::ToPrimitive::to_f64(&v)?)
        }
        Binomial => gamma::binomial_c(x, a[1])?,
        RisingFactorial => gamma::rf_c(x, a[1])?,
        FallingFactorial => gamma::ff_c(x, a[1])?,
        Harmonic => gamma::harmonic_c(x, a.get(1).copied().unwrap_or(cr(1.0)))?,
        Bernoulli => {
            if let Some(n) = int_arg(x) {
                if n >= 0 && a.len() == 1 {
                    let b = super::integers::bernoulli_number(n as usize);
                    return Some(cr(ratio_f64(&b)));
                }
            }
            gamma::bernoulli_c(x, a.get(1).copied())?
        }
        Euler => {
            if let Some(n) = int_arg(x) {
                if n >= 0 && a.len() == 1 {
                    let e = super::integers::euler_number(n as usize);
                    return Some(cr(num_traits::ToPrimitive::to_f64(&e)?));
                }
            }
            gamma::euler_c(x, a.get(1).copied())?
        }
        Catalan => gamma::catalan_c(x)?,
        Fibonacci => {
            if a.len() == 2 {
                return None;
            }
            gamma::fibonacci_c(x)
        }
        Lucas => gamma::lucas_c(x),
        Bell | Partition => {
            let n = int_arg(x)?;
            if n < 0 || a.len() != 1 {
                return None;
            }
            let v = if s == Bell {
                super::integers::bell_number(n as u64)
            } else {
                super::integers::partition_number(n as u64)
            };
            cr(num_traits::ToPrimitive::to_f64(&v)?)
        }
        Polygamma => gamma::polygamma_c(x, a[1])?,
        Digamma => gamma::digamma_c(x)?,
        Trigamma => gamma::polygamma_c(cr(1.0), x)?,
        Beta => gamma::beta_c(x, *a.get(1).unwrap_or(&x))?,
        LowerGamma | UpperGamma => {
            let (aa, z) = (x, a[1]);
            if !(is_real_c(aa) && is_real_c(z)) {
                return None;
            }
            cr(if s == LowerGamma {
                gamma::lowergamma_real(aa.re, z.re)?
            } else {
                gamma::uppergamma_real(aa.re, z.re)?
            })
        }
        Zeta => {
            if a.len() == 2 {
                gamma::hurwitz_zeta(x, a[1])?
            } else {
                gamma::zeta_c(x)?
            }
        }
        DirichletEta => gamma::eta_c(x)?,
        Polylog => gamma::polylog_c(x, a[1])?,
        Lerchphi => gamma::lerchphi_c(x, a[1], a[2])?,
        Erfi => erf::erfi_c(x),
        Erf2 => erf::erf_c(a[1]) - erf::erf_c(x),
        Erfinv => cr(erf::erfinv_real(real(x)?)?),
        Erfcinv => cr(erf::erfcinv_real(real(x)?)?),
        FresnelS => erf::fresnel_c(x)?.0,
        FresnelC => erf::fresnel_c(x)?.1,
        Ei => erf::ei_c(x)?,
        Li => erf::li_c(x)?,
        LiOffset => erf::li_c(x)? - erf::li_c(cr(2.0))?,
        Si => erf::trig_integral(0, x)?,
        Ci => erf::trig_integral(1, x)?,
        Shi => erf::trig_integral(2, x)?,
        Chi => erf::trig_integral(3, x)?,
        Expint => erf::expint_c(x, a[1])?,
        E1 => erf::e1_c(x)?,
        BesselJ => bessel::besselj_c(x, a[1])?,
        BesselY => bessel::bessely_c(x, a[1])?,
        BesselI => bessel::besseli_c(x, a[1])?,
        BesselK => bessel::besselk_cc(x, a[1])?,
        Hankel1 => bessel::hankel_c(x, a[1], true)?,
        Hankel2 => bessel::hankel_c(x, a[1], false)?,
        Jn => bessel::spherical_c(x, a[1], true)?,
        Yn => bessel::spherical_c(x, a[1], false)?,
        AiryAi => bessel::airy_c(0, x)?,
        AiryBi => bessel::airy_c(1, x)?,
        AiryAiPrime => bessel::airy_c(2, x)?,
        AiryBiPrime => bessel::airy_c(3, x)?,
        EllipticK => misc::ellipk_c(x)?,
        EllipticE => {
            if a.len() == 1 {
                misc::ellipe_c(x)?
            } else {
                misc::ellipe2_c(x, a[1])?
            }
        }
        EllipticF => misc::ellipf_c(x, a[1])?,
        EllipticPi => {
            if a.len() == 2 {
                misc::ellippi_c(x, None, a[1])?
            } else {
                misc::ellippi_c(x, Some(a[1]), a[2])?
            }
        }
        LambertW => {
            let k = match a.get(1) {
                Some(k) => int_arg(*k)?,
                None => 0,
            };
            misc::lambertw_c(x, k)?
        }
        Hyper => return None,
        Legendre => misc::legendre_c(x, a[1])?,
        AssocLegendre => misc::assoc_legendre_c(x, a[1], a[2])?,
        ChebyshevT => misc::chebyshev_c(true, x, a[1])?,
        ChebyshevU => misc::chebyshev_c(false, x, a[1])?,
        Hermite => misc::hermite_c(x, a[1])?,
        Laguerre => {
            int_arg(x)?;
            misc::laguerre_c(x, cr(0.0), a[1])?
        }
        AssocLaguerre => misc::laguerre_c(x, a[1], a[2])?,
        Jacobi => misc::jacobi_c(x, a[1], a[2], a[3])?,
        Gegenbauer => misc::gegenbauer_c(x, a[1], a[2])?,
        Ynm => misc::ynm_c(x, a[1], a[2], a[3])?,
    };
    Some(r)
}

fn real(z: C64) -> Option<f64> {
    if z.im == 0.0 {
        Some(z.re)
    } else {
        None
    }
}

fn ratio_f64(r: &num_rational::BigRational) -> f64 {
    use num_traits::ToPrimitive;
    let n = r.numer().to_f64().unwrap_or(f64::NAN);
    let d = r.denom().to_f64().unwrap_or(f64::NAN);
    if n.is_finite() && d.is_finite() {
        n / d
    } else {
        // big values: scale
        let nb = r.numer().bits() as i64;
        let db = r.denom().bits() as i64;
        let sh = (nb - db) - 60;
        let num = if sh > 0 { r.numer() >> (sh as usize) } else { r.numer() << ((-sh) as usize) };
        (num.to_f64().unwrap() / r.denom().to_f64().unwrap_or(f64::INFINITY)) * 2f64.powi(sh as i32)
    }
}

/// Is `e` a `Func::Special` node? (helper for the evaluators)
pub fn special_head(e: &Expr) -> Option<SpecialFn> {
    match e.kind() {
        Kind::Func(Func::Special(s), _) => Some(*s),
        _ => None,
    }
}
