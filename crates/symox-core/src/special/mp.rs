//! Arbitrary precision evaluation through exact gamma/elementary identities.
//!
//! Unsupported families return `None`; never pad a double precision result
//! with digits that were not computed.

use crate::evalf_mp::{ev, MpComplex};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

use crate::expr::{Expr, Kind};
use crate::number::Number;
use crate::functions::{cos, gamma};
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sqrt, sub};

use super::SpecialFn;

/// Evaluate supported heads to `prec` bits without changing their symbolic form.
/// Gamma ratios require real arguments; elementary identities may be complex.
pub fn eval(s: SpecialFn, args: &[Expr], prec: u32) -> Option<MpComplex> {
    if !s.nargs().contains(&args.len()) {
        return None;
    }
    if s == SpecialFn::LambertW {
        return super::lambert::mp(args, prec);
    }
    // SymPy hands mpmath the exact binary value of a Float argument. The
    // identities below are built from exact rationals, so that a factor such
    // as gamma(0.75) is not evaluated at the Float's own precision.
    let args: Vec<Expr> = args.iter().map(exact_floats).collect();
    let args = &args[..];
    let x = &args[0];
    let one = Expr::one();
    let plus_one = |a: &Expr| add2(a, &one);
    let expr = match s {
        SpecialFn::Factorial => gamma(&plus_one(x)),
        SpecialFn::Binomial => div(
            &gamma(&plus_one(x)),
            &mul2(&gamma(&plus_one(&args[1])), &gamma(&plus_one(&sub(x, &args[1])))),
        ),
        SpecialFn::RisingFactorial => div(&gamma(&add2(x, &args[1])), &gamma(x)),
        SpecialFn::FallingFactorial => div(&gamma(&plus_one(x)), &gamma(&plus_one(&sub(x, &args[1])))),
        SpecialFn::Beta => {
            let y = args.get(1).unwrap_or(x);
            div(&mul2(&gamma(x), &gamma(y)), &gamma(&add2(x, y)))
        }
        SpecialFn::Catalan => div(
            &mul2(&pow(&Expr::int(4), x), &gamma(&add2(x, &Expr::rational(1, 2)))),
            &mul2(&sqrt(&Expr::pi()), &gamma(&add2(x, &Expr::int(2)))),
        ),
        SpecialFn::Fibonacci if args.len() == 1 => {
            let phi = div(&add2(&one, &sqrt(&Expr::int(5))), &Expr::int(2));
            div(
                &sub(&pow(&phi, x), &mul2(&cos(&mul2(&Expr::pi(), x)), &pow(&phi, &neg(x)))),
                &sqrt(&Expr::int(5)),
            )
        }
        SpecialFn::Lucas => {
            let phi = div(&add2(&one, &sqrt(&Expr::int(5))), &Expr::int(2));
            add2(&pow(&phi, x), &mul2(&cos(&mul2(&Expr::pi(), x)), &pow(&phi, &neg(x))))
        }
        SpecialFn::Mod | SpecialFn::Rem | SpecialFn::Sech | SpecialFn::Csch | SpecialFn::Asec | SpecialFn::Acsc
        | SpecialFn::Asech | SpecialFn::Acsch | SpecialFn::Acoth | SpecialFn::Sinc | SpecialFn::Frac
        | SpecialFn::KroneckerDelta => super::elementary::mp_rewrite(s, args, prec)?,
        _ => return None,
    };
    ev(&expr, prec)
}

/// `e` with each Float leaf replaced by the exact rational value of its
/// binary representation (moderate magnitudes only).
fn exact_floats(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Num(n) if n.is_float() => {
            let value = n.to_bigfloat(n.float_precision().unwrap_or(53));
            if value.is_zero() {
                return Expr::zero();
            }
            if value.mag().abs() > 4096 || value.exp.abs() > 8192 {
                return e.clone();
            }
            let man = value.man.clone();
            let ratio = if value.exp >= 0 {
                BigRational::from_integer(man << (value.exp as usize))
            } else {
                BigRational::new(man, BigInt::one() << ((-value.exp) as usize))
            };
            Expr::num(Number::from_ratio(ratio))
        }
        Kind::Add(..) => add(&e.args().iter().map(exact_floats).collect::<Vec<_>>()),
        Kind::Mul(..) => mul(&e.args().iter().map(exact_floats).collect::<Vec<_>>()),
        Kind::Pow(base, exponent) => pow(&exact_floats(base), &exact_floats(exponent)),
        _ => e.clone(),
    }
}
