//! The repeat loop of SymPy's `expand`.
//!
//! After applying its hints once, SymPy's `Expr.expand` repeats the
//! `multinomial` and `mul` hints until the expression no longer changes.
//! One pass can create new work for the next one: distributing
//! `(x + y)**(2/3)*(1 + (x + y)**(2/3))` builds `(x + y)**(4/3)`, which the
//! multinomial hint turns into `x*(x + y)**(1/3) + y*(x + y)**(1/3)`.
//! This module detects such nodes so the expander repeats only when needed.

use rustc_hash::FxHashSet;

use crate::expr::{Expr, Kind};

use super::ExpandOptions;

/// Bound on the extra passes; each one only acts on the nodes built by the
/// previous one, and a pass without change stops the loop.
pub(super) const MAX_PASSES: usize = 8;

/// Whether a further `mul`/`multinomial` pass may change `e`: a sum raised
/// to a rational power above 1, or a product holding a commutative sum.
pub(super) fn needs_pass(e: &Expr, options: &ExpandOptions) -> bool {
    let mut seen = FxHashSet::default();
    scan(e, options, &mut seen)
}

fn sum_power(base: &Expr, exponent: &Expr, options: &ExpandOptions, factor: bool) -> bool {
    if !base.is_add() {
        return false;
    }
    let Some((p, q)) = exponent.as_number().and_then(|n| n.numer_denom()) else {
        return false;
    };
    let unit = p == q;
    (options.multinomial && p > q)
        || (factor && unit && options.mul && crate::decompose::commutes(base))
}

fn scan(e: &Expr, options: &ExpandOptions, seen: &mut FxHashSet<usize>) -> bool {
    if e.is_atom() || !seen.insert(e.dag_id()) {
        return false;
    }
    let found = match e.kind() {
        Kind::Mul(_, factors) => factors.iter().any(|(b, x)| sum_power(b, x, options, true)),
        Kind::Pow(b, x) => sum_power(b, x, options, false),
        _ => false,
    };
    if found {
        return true;
    }
    // Without `deep` only the root is expanded, as in the reference.
    options.deep && e.args().iter().any(|a| scan(a, options, seen))
}

#[cfg(test)]
mod tests {
    use crate::parser::parse_str;
    use crate::printer::to_str;

    fn expanded(s: &str) -> String {
        to_str(&crate::expand::expand(&parse_str(s).unwrap()))
    }

    #[test]
    fn products_of_fractional_sum_powers_are_expanded_again() {
        // SymPy 1.14 expand values.
        assert_eq!(expanded("(x + y)**(2/3)*(1 + (x + y)**(2/3))"), "x*(x + y)**(1/3) + y*(x + y)**(1/3) + (x + y)**(2/3)");
        assert_eq!(
            expanded("(3*sqrt(93) + 31)**(2/3)*(1 + (3*sqrt(93) + 31)**(2/3))"),
            "(3*sqrt(93) + 31)**(2/3) + 3*sqrt(93)*(3*sqrt(93) + 31)**(1/3) + 31*(3*sqrt(93) + 31)**(1/3)"
        );
        assert_eq!(expanded("(x + y)*((x + y)**(1/3) + 1)"), "x*(x + y)**(1/3) + x + y*(x + y)**(1/3) + y");
    }
}
