//! Arbitrary precision evaluation through exact gamma/elementary identities.
//!
//! Unsupported families return `None`; never pad a double precision result
//! with digits that were not computed.

use crate::evalf_mp::{ev, MpComplex};
use crate::expr::Expr;
use crate::functions::{cos, gamma};
use crate::ops::{add2, div, mul2, neg, pow, sqrt, sub};

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
