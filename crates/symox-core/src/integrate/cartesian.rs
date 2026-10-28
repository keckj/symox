//! Exact Cartesian numeric coefficients for heurisch's conjugate-log pairs.
//!
//! Principal square roots stay compact in the general root solver. Heurisch
//! needs their algebraic real and imaginary parts locally, before collecting
//! conjugate polynomial factors into real logarithms and inverse tangents.

use rustc_hash::FxHashMap;

use crate::complex_parts::{as_real_imag, PartsError};
use crate::expr::{Constant, Expr, Kind};
use crate::ops::{add, div, mul, mul2, pow, sqrt, sub};

type Parts = (Expr, Expr);

#[derive(Default)]
struct Splitter {
    // Retain the original nodes along with their physical keys: sharing and
    // temporary expression allocation cannot invalidate the memo entries.
    memo: FxHashMap<usize, (Expr, Option<Parts>)>,
}

impl Splitter {
    fn parts(&mut self, coefficient: &Expr) -> Result<Option<Parts>, PartsError> {
        if let Some((_, result)) = self.memo.get(&coefficient.dag_id()) {
            return Ok(result.clone());
        }
        let result = match coefficient.kind() {
            Kind::Add(..) => {
                let mut real = Vec::new();
                let mut imaginary = Vec::new();
                for term in coefficient.args() {
                    let Some((r, i)) = self.parts(&term)? else { return Ok(None) };
                    real.push(r);
                    imaginary.push(i);
                }
                Some((add(&real), add(&imaginary)))
            }
            Kind::Mul(..) => {
                let (mut real, mut imaginary) = (Expr::one(), Expr::zero());
                for factor in coefficient.args() {
                    let Some((r, i)) = self.parts(&factor)? else { return Ok(None) };
                    let next_real = sub(&mul2(&real, &r), &mul2(&imaginary, &i));
                    imaginary = add(&[mul2(&real, &i), mul2(&imaginary, &r)]);
                    real = next_real;
                }
                Some((real, imaginary))
            }
            Kind::Pow(base, exponent) if exponent == &Expr::rational(1, 2)
                && coefficient.float_info().max_precision.is_none() => {
                if let Some((real, imaginary)) = self.parts(base)? {
                    if crate::assumptions::is_zero(&imaginary) == Some(false) {
                        // For a+I*b with b != 0, r=sqrt(a*a+b*b)>|a|.
                        // sqrt((r+a)/2)>0, and b/(2*root_real) chooses the
                        // sign of the principal root without a floating test.
                        let radius = sqrt(&crate::expand::expand(&add(&[
                            pow(&real, &Expr::int(2)), pow(&imaginary, &Expr::int(2)),
                        ])));
                        let root_real = sqrt(&div(&add(&[radius, real]), &Expr::int(2)));
                        Some((root_real.clone(), div(&imaginary, &mul2(&Expr::int(2), &root_real))))
                    } else {
                        as_real_imag(coefficient, true, None)?
                    }
                } else {
                    None
                }
            }
            _ => as_real_imag(coefficient, true, None)?,
        };
        self.memo.insert(coefficient.dag_id(), (coefficient.clone(), result.clone()));
        Ok(result)
    }
}

/// Split only closed numeric polynomial coefficients into their exact real
/// and imaginary components. Unsupported polynomial forms, parameters, and
/// closed undefined functions keep their original expressions.
pub fn numeric_coefficients(polynomial: &Expr, generators: &[Expr]) -> Result<Expr, PartsError> {
    numeric_coefficients_with_status(polynomial, generators).map(|(value, _)| value)
}

/// Also report whether any coefficient was successfully split. This is a
/// semantic processing flag, not expression inequality: z+I is processed
/// even though its Cartesian form is already identical to the input.
pub fn numeric_coefficients_with_status(polynomial: &Expr, generators: &[Expr]) -> Result<(Expr, bool), PartsError> {
    crate::assumptions::exact_queries(|| {
        if generators.is_empty() || crate::poly::convert::check_commutative(polynomial).is_err() {
            return Ok((polynomial.clone(), false));
        }
        let Ok(coefficients) = crate::poly::expression_domain::dict_from_expr(polynomial, generators) else {
            return Ok((polynomial.clone(), false));
        };
        let mut splitter = Splitter::default();
        let mut terms = Vec::new();
        let mut changed = false;
        for (monomial, coefficient) in coefficients {
            let power = mul(&generators.iter().zip(monomial).map(|(generator, exponent)|
                pow(generator, &Expr::int(i64::from(exponent)))).collect::<Vec<_>>());
            if crate::symbolic_number::is_symbolic_number(&coefficient)
                && crate::traverse::free_symbols(&coefficient).is_empty()
                && (crate::traverse::has(&coefficient, &Expr::i())
                    || crate::assumptions::is_real(&coefficient) == Some(false)) {
                if let Some((real, imaginary)) = splitter.parts(&coefficient)? {
                    if !crate::traverse::has(&real, &Expr::i()) && !crate::traverse::has(&imaginary, &Expr::i()) {
                        terms.push(mul2(&real, &power));
                        terms.push(mul(&[Expr::constant(Constant::I), imaginary, power]));
                        changed = true;
                        continue;
                    }
                }
            }
            terms.push(mul2(&coefficient, &power));
        }
        Ok((if changed { add(&terms) } else { polynomial.clone() }, changed))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(value: &str) -> Expr { parse_str(value).unwrap() }

    #[test]
    fn principal_square_roots_expose_exact_conjugate_components() {
        let z = Expr::symbol("z");
        for (source, expected) in [
            ("z + sqrt(1/2 + sqrt(3)*I/2)", "z + sqrt(3)/2 + I/2"),
            ("z + sqrt(1/2 - sqrt(3)*I/2)", "z + sqrt(3)/2 - I/2"),
            ("z + sqrt(-3 + 4*I)", "z + 1 + 2*I"),
            ("z + sqrt(-3 - 4*I)", "z + 1 - 2*I"),
        ] {
            assert_eq!(numeric_coefficients(&p(source), std::slice::from_ref(&z)).unwrap(), p(expected));
        }
    }

    #[test]
    fn unknown_coefficients_and_nonpolynomial_forms_keep_their_presentation() {
        let z = Expr::symbol("z");
        for source in ["z + f(1) + I", "z + sqrt(a + I)", "1/z + sqrt(1 + I)", "sin(z) + I"] {
            let original = p(source);
            assert_eq!(numeric_coefficients(&original, std::slice::from_ref(&z)).unwrap(), original);
        }
    }

    #[test]
    fn processing_status_distinguishes_equal_cartesian_forms_from_untouched_inputs() {
        let z = Expr::symbol("z");
        for (source, processed) in [("z + 1", false), ("z + f(1) + I", false), ("z + I", true)] {
            let original = p(source);
            assert_eq!(numeric_coefficients_with_status(&original, std::slice::from_ref(&z)).unwrap(), (original, processed));
        }
    }
}
