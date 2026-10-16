//! Local complex analysis through the native Laurent-series engine.

use crate::expr::Expr;
use crate::ops::add;
use crate::traverse::has;

/// Coefficient of `(x-point)^-1` at an isolated pole (including higher poles).
/// Branch points and essential singularities unsupported by the series engine
/// return an error instead of a guessed coefficient.
pub fn residue(expression: &Expr, variable: &Expr, point: &Expr) -> Result<Expr, String> {
    if !variable.is_symbol() || has(point, variable) {
        return Err("residue requires a symbol and a point independent of it".into());
    }
    let (terms, _) = crate::series::series_coeffs(expression, variable, point, 1, false)
        .map_err(|e| format!("residue: {e}"))?;
    let mut coefficients = Vec::new();
    for (coefficient, exponent) in terms {
        if has(&coefficient, variable) || !exponent.is_integer() {
            return Err("residue requires a Laurent expansion at an isolated singularity".into());
        }
        if exponent.is_minus_one() {
            coefficients.push(coefficient);
        }
    }
    Ok(add(&coefficients))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    #[test]
    fn isolated_poles_and_regular_points() {
        let x = Expr::symbol("x");
        for (input, point, answer) in [
            ("exp(x)/x**5", "0", "1/24"),
            ("1/(x*(x+1))", "0", "1"),
            ("1/(x*(x+1))", "-1", "-1"),
            ("sin(x)/x", "0", "0"),
            ("1/sin(x)", "0", "1"),
        ] {
            let result = residue(&parse_str(input).unwrap(), &x, &parse_str(point).unwrap()).unwrap();
            assert_eq!(result, parse_str(answer).unwrap(), "{input}");
        }
        assert!(residue(&parse_str("sqrt(x)/x").unwrap(), &x, &Expr::zero()).is_err());
    }
}
