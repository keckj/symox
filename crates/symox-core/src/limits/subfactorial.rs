//! A certified, uniformly bounded remainder for derangements.
//!
//! On the principal continuation `!a = Gamma(a + 1, -1)/e`, for real
//! `a > -1`,
//!
//! `!a - Gamma(a + 1)/e = exp(I*pi*a)/e * integral_0^1 u^a exp(u) du`.
//!
//! Consequently its modulus is at most `1/(a + 1)`. This also follows
//! directly from the alternating exponential series for integral `a`.
//! An affine occurrence may therefore be replaced by `Gamma(a + 1)/e`
//! only after proving that its coefficient divided by `a + 1` tends to
//! zero. Keeping this certificate is essential: multiplying the remainder
//! by `a` exposes its nonvanishing oscillation.

use super::{affine_parts, gruntz::Gruntz};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions as fx;
use crate::ops::{add, div, mul2};
use crate::series::engine::depends;
use crate::special::SpecialFn;

pub(super) fn value(e: &Expr, t: &Expr) -> Option<Expr> {
    fn collect(e: &Expr, t: &Expr, out: &mut Vec<Expr>) {
        if !depends(e, t) || out.contains(e) {
            return;
        }
        if matches!(e.kind(), Kind::Func(Func::Special(SpecialFn::Subfactorial), _)) {
            out.push(e.clone());
            return;
        }
        for a in e.args() {
            collect(&a, t, out);
        }
    }

    let mut nodes = Vec::new();
    collect(e, t, &mut nodes);
    if nodes.is_empty() {
        return None;
    }
    // This also excludes products, powers, and singular functions of the
    // remainders; continuity alone would not justify discarding an error
    // in such expressions.
    let (constant, coefficients) = affine_parts(e, &nodes)?;
    let mut g = Gruntz::new(t);
    let positive = Expr::dummy("subfactorial_path", crate::symbol::Assumptions::POSITIVE);
    let mut terms = vec![constant];
    for (node, coefficient) in nodes.iter().zip(coefficients) {
        let Kind::Func(_, args) = node.kind() else { return None };
        let [a] = args.as_slice() else { return None };
        // A syntactically real formula need not be real on this path:
        // log(-t), for example, contains no explicit I. Ask assumptions
        // after making the positive limit variable explicit.
        if crate::assumptions::is_real(&crate::subs::subs1(a, t, &positive)) != Some(true)
            || !g.limitinf(a).ok()?.is_constant(Constant::Infinity)
        {
            return None;
        }
        let a1 = add(&[a.clone(), Expr::one()]);
        if !g.limitinf(&div(&coefficient, &a1)).ok()?.is_zero() {
            return None;
        }
        let main = mul2(&fx::gamma(&a1), &fx::exp(&Expr::minus_one()));
        terms.push(mul2(&coefficient, &main));
    }
    g.limitinf(&add(&terms)).ok()
}
