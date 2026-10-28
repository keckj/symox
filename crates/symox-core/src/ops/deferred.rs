//! Exact held numeric coefficients produced by large-order recipes.
//! Ordinary pow/num_pow are unchanged. The budget selects a representation
//! only for these new recipes; explicit doit/int use the complete evaluator.
use crate::{expr::Kind, number::Number, Expr};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

pub(crate) const RECIPE_BITS: u64 = 1 << 20;

fn eligible(base: &Expr, exponent: &Expr) -> bool {
    exponent.as_number().is_some_and(Number::is_integer)
        && !base.has_host()
        && !base.has_nonfinite()
        && !base.has_matrix()
        && base.float_info().max_precision.is_none()
        && crate::functions::is_numeric_value(base)
        && crate::assumptions::exact_queries(|| crate::assumptions::is_finite(base)) == Some(true)
        && (!exponent.as_number().is_some_and(Number::is_negative)
            || crate::assumptions::exact_queries(|| crate::assumptions::is_zero(base))
                == Some(false))
}
/// Restore an existing recipe without eager arithmetic. No host/Float power
/// can acquire this metadata through the pickle protocol.
pub fn restore(base: &Expr, exponent: &Expr) -> Option<Expr> {
    if !eligible(base, exponent) {
        return None;
    }
    if exponent.is_zero()
        || exponent.is_one()
        || base.is_zero()
        || base.is_one()
        || base.is_minus_one()
    {
        return Some(super::pow(base, exponent));
    }
    Some(Expr::deferred_numeric_power(base.clone(), exponent.clone()))
}
/// Preserve a recipe when a transformation changes its operands. Trivial,
/// inexact, nonfinite, symbolic and host cases retain the ordinary evaluator.
pub fn rebuild(base: &Expr, exponent: &Expr) -> Expr {
    let Some(n) = exponent.as_number().and_then(Number::to_bigint) else {
        return super::pow(base, exponent);
    };
    if n.is_zero()
        || n.is_one()
        || base.is_zero()
        || base.is_one()
        || base.is_minus_one()
        || !eligible(base, exponent)
    {
        return super::pow(base, exponent);
    }
    let large = if let Some((a, b)) = base.as_number().and_then(Number::numer_denom) {
        n.abs() * BigInt::from(a.bits().max(b.bits()).max(1)) > BigInt::from(RECIPE_BITS)
    } else {
        n.abs() > BigInt::from(64u8)
    };
    if large {
        Expr::deferred_numeric_power(base.clone(), exponent.clone())
    } else {
        super::pow(base, exponent)
    }
}
/// The recipe already proved that holding is appropriate. Source constants
/// and variable polynomial powers use their existing representation instead.
pub(crate) fn hold(base: &Expr, exponent: &Expr) -> Expr {
    restore(base, exponent)
        .unwrap_or_else(|| Expr::from_unevaluated_kind(Kind::Pow(base.clone(), exponent.clone())))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn held_recipe_has_no_effect_on_direct_numeric_power_or_explicit_doit() {
        let b = Expr::int(2);
        let n = Expr::int(12);
        let held = restore(&b, &n).unwrap();
        assert!(held.is_deferred_numeric_power());
        assert_eq!(super::super::pow(&b, &n), Expr::int(4096));
        assert_eq!(crate::functions::doit(&held, true), Expr::int(4096));
        assert_eq!(crate::functions::doit(&held, false), Expr::int(4096));
        let x = Expr::symbol("x");
        let product = super::super::mul2(&held, &x);
        assert!(product.has_deferred_numeric_power());
        assert!(crate::subs::subs1(&product, &x, &Expr::one()).is_deferred_numeric_power());
        assert_eq!(
            crate::subs::subs1(&product, &x, &Expr::zero()),
            Expr::zero()
        );
        assert_eq!(super::super::div(&held, &held), Expr::one());
    }
    #[test]
    fn reconstruction_restores_exact_trivial_and_float_paths() {
        let n = Expr::int(1_000_000);
        assert_eq!(rebuild(&Expr::zero(), &n), Expr::zero());
        assert_eq!(rebuild(&Expr::one(), &n), Expr::one());
        assert_eq!(rebuild(&Expr::minus_one(), &n), Expr::one());
        assert_eq!(rebuild(&Expr::int(2), &Expr::zero()), Expr::one());
        assert_eq!(rebuild(&Expr::int(2), &Expr::one()), Expr::int(2));
        assert_eq!(rebuild(&Expr::float(2.0), &Expr::int(3)), Expr::float(8.0));
        assert!(restore(&Expr::symbol("x"), &n).is_none());
        assert!(restore(&Expr::float(2.0), &n).is_none());
    }
}
