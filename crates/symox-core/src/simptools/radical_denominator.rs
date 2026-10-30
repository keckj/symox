//! The denominator step of SymPy 1.14's `simplify`.
//!
//! After its chain of `shorter` candidates and the hollow-Mul pass, SymPy's
//! `simplify` rationalizes a denominator that is a sum:
//!
//! ```python
//! numer, denom = expr.as_numer_denom()
//! if denom.is_Add:
//!     n, d = fraction(radsimp(1/denom, symbolic=False, max_terms=1))
//!     if n is not S.One:
//!         expr = (numer*n).expand()/d
//! ```
//!
//! `radsimp` itself first replaces a numeric sum in a denominator that is
//! exactly a rational number (`prepare_numeric_denominators`). A sum of
//! square roots is decided natively; any other exact form is passed to the
//! caller's `fallback`, which may decline it.

use crate::expr::Expr;
use crate::ops;

use super::radsimp::{fraction, prepare_numeric_denominators, radsimp};

/// SymPy's `simplify` step on `expr`. `check` is polled between nodes and
/// `fallback` recognizes a numeric denominator that is not a sum of square
/// roots (`Ok(None)` keeps it). An error of either stops the step.
///
/// The reference tests `n is not S.One`. Its own pipeline has already
/// negated a denominator that the native one may not have, so a denominator
/// proven to be a Number (e.g. `3 - (-1)**(2/3) + (-1)**(1/3) == 4`) is also
/// taken.
pub fn rationalize_sum_denominator<E>(
    expr: &Expr,
    check: &mut impl FnMut() -> Result<(), E>,
    fallback: &mut impl FnMut(&Expr) -> Result<Option<Expr>, E>,
) -> Result<Expr, E> {
    let (numer, denom) = super::as_numer_denom(expr);
    if !denom.is_add() {
        return Ok(expr.clone());
    }
    let inverse = prepare_numeric_denominators(&ops::inv(&denom), check, fallback)?;
    let (n, d) = fraction(&radsimp(&inverse, false, 1), false);
    if !n.is_one() || (super::is_number_atom(&d) && d != denom) {
        let product = ops::mul(&[numer, n]);
        let options = crate::expand::ExpandOptions::default();
        let mut budget = crate::expand::ExpandBudget::unlimited();
        // SymPy's `expand` has no size limit; the unlimited budget never
        // fails, but an error keeps the expression unchanged.
        let Ok(expanded) = crate::expand::expand_with_options(&product, &options, &mut budget) else {
            return Ok(expr.clone());
        };
        return Ok(ops::div(&expanded, &d));
    }
    Ok(expr.clone())
}

/// The step without a host: a numeric denominator outside the sums of
/// square roots declines the step (the expression is returned unchanged),
/// since its exact recognition needs `nsimplify`.
pub fn rationalize_sum_denominator_native(expr: &Expr) -> Expr {
    let mut check = || Ok::<(), ()>(());
    let mut decline = |_: &Expr| Err::<Option<Expr>, ()>(());
    rationalize_sum_denominator(expr, &mut check, &mut decline).unwrap_or_else(|_| expr.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn step(s: &str) -> String {
        crate::printer::to_str(&rationalize_sum_denominator_native(&parse_str(s).unwrap()))
    }

    #[test]
    fn numeric_radical_denominators_follow_sympy() {
        // SymPy 1.14 simplify values.
        assert_eq!(step("(2 + sqrt(5))/(50*(1 + sqrt(5)))"), "sqrt(5)/200 + 3/200");
        assert_eq!(step("1/(1 + sqrt(2))"), "-1 + sqrt(2)");
        assert_eq!(step("(1 + sqrt(3))/(2 - sqrt(3))"), "5 + 3*sqrt(3)");
    }

    #[test]
    fn other_denominators_are_kept() {
        for s in ["x/(sqrt(2) + sqrt(3))", "1/(x + sqrt(2))", "1/(x + 1)", "x + 1", "sqrt(2)/3"] {
            let e = parse_str(s).unwrap();
            assert_eq!(rationalize_sum_denominator_native(&e), e, "{s}");
        }
    }
}
