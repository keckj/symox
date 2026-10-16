//! Strict linear coefficient extraction and parametric exact linear solutions.

use rustc_hash::FxHashMap;

use crate::expr::{Expr, Func, Kind, RelOp};
use crate::matrix::{Matrix, MatrixError};
use crate::ops::{add, add2, mul2, neg, sub};

#[derive(Clone, Debug)]
pub enum LinearError {
    Nonlinear(Expr),
    Invalid(String),
    Matrix(MatrixError),
}

#[derive(Clone)]
struct Form {
    coefficients: Vec<Expr>,
    constant: Expr,
}

fn depends(e: &Expr, variables: &FxHashMap<Expr, usize>, cache: &mut FxHashMap<Expr, bool>) -> bool {
    if let Some(result) = cache.get(e) {
        return *result;
    }
    let result = variables.contains_key(e) || e.args().iter().any(|a| depends(a, variables, cache));
    cache.insert(e.clone(), result);
    result
}

fn form(
    e: &Expr,
    variables: &FxHashMap<Expr, usize>,
    dependency: &mut FxHashMap<Expr, bool>,
    cache: &mut FxHashMap<Expr, Form>,
) -> Result<Form, LinearError> {
    if let Some(result) = cache.get(e) {
        return Ok(result.clone());
    }
    let n = variables.len();
    let mut result = Form { coefficients: vec![Expr::zero(); n], constant: Expr::zero() };
    if let Some(&index) = variables.get(e) {
        result.coefficients[index] = Expr::one();
    } else if !depends(e, variables, dependency) {
        result.constant = e.clone();
    } else {
        match e.kind() {
            Kind::Add(c, terms) => {
                result.constant = Expr::num(c.clone());
                for (term, coefficient) in terms {
                    let part = form(term, variables, dependency, cache)?;
                    let coefficient = Expr::num(coefficient.clone());
                    for (slot, value) in result.coefficients.iter_mut().zip(part.coefficients) {
                        *slot = add2(slot, &mul2(&coefficient, &value));
                    }
                    result.constant = add2(&result.constant, &mul2(&coefficient, &part.constant));
                }
            }
            Kind::Mul(c, factors) => {
                let mut scalar = Expr::num(c.clone());
                let mut dependent = None;
                for (base, exponent) in factors {
                    if depends(base, variables, dependency) || depends(exponent, variables, dependency) {
                        if !exponent.is_one() || dependent.is_some() {
                            return Err(LinearError::Nonlinear(e.clone()));
                        }
                        dependent = Some(form(base, variables, dependency, cache)?);
                    } else {
                        scalar = mul2(&scalar, &crate::ops::pow(base, exponent));
                    }
                }
                if let Some(part) = dependent {
                    result.coefficients = part.coefficients.iter().map(|x| mul2(&scalar, x)).collect();
                    result.constant = mul2(&scalar, &part.constant);
                } else {
                    result.constant = scalar;
                }
            }
            Kind::Pow(base, exponent) if exponent.is_one() => {
                result = form(base, variables, dependency, cache)?;
            }
            Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") && args.len() == 1 => {
                result = form(&args[0], variables, dependency, cache)?;
            }
            _ => return Err(LinearError::Nonlinear(e.clone())),
        }
    }
    cache.insert(e.clone(), result.clone());
    Ok(result)
}

/// Extract A and b from equations A*x=b without cancelling nonlinear factors.
pub fn linear_eq_to_matrix(equations: &[Expr], symbols: &[Expr]) -> Result<(Matrix, Matrix), LinearError> {
    if symbols.is_empty() {
        return Err(LinearError::Invalid("Symbols must be given".into()));
    }
    let mut variables = FxHashMap::default();
    for (index, symbol) in symbols.iter().enumerate() {
        if variables.insert(symbol.clone(), index).is_some() {
            return Err(LinearError::Invalid("Symbols must be unique".into()));
        }
    }
    let mut coefficients = Vec::with_capacity(equations.len() * symbols.len());
    let mut constants = Vec::with_capacity(equations.len());
    let mut dependency = FxHashMap::default();
    let mut cache = FxHashMap::default();
    for equation in equations {
        let expression = match equation.kind() {
            Kind::Rel(RelOp::Eq, lhs, rhs) => sub(lhs, rhs),
            Kind::Bool(_) | Kind::Rel(..) => return Err(LinearError::Invalid("Expected expressions or equality equations".into())),
            _ => equation.clone(),
        };
        let result = form(&expression, &variables, &mut dependency, &mut cache)?;
        coefficients.extend(result.coefficients);
        constants.push(neg(&result.constant));
    }
    Ok((Matrix::new_unchecked(equations.len(), symbols.len(), coefficients),
        Matrix::new_unchecked(equations.len(), 1, constants)))
}

/// Solve A*x=b, retaining supplied variables in free columns.
pub fn linsolve(a: &Matrix, b: &Matrix, symbols: &[Expr]) -> Result<Option<Vec<Expr>>, LinearError> {
    if b.cols() != 1 || a.rows() != b.rows() || symbols.len() != a.cols() {
        return Err(LinearError::Invalid("Linear system matrix dimensions do not match the symbols".into()));
    }
    let n = a.cols();
    let augmented = Matrix::hstack(&[a, b]).map_err(LinearError::Matrix)?;
    let (reduced, pivots) = augmented.rref();
    if pivots.iter().any(|&column| column >= n) {
        return Ok(None);
    }
    let free: Vec<usize> = (0..n).filter(|column| !pivots.contains(column)).collect();
    let mut result = symbols.to_vec();
    for (row, &column) in pivots.iter().enumerate() {
        let terms: Vec<Expr> = free.iter().map(|&j| mul2(reduced.get(row, j), &symbols[j])).collect();
        // Preserve the affine sum rather than collecting its terms into one
        // fraction; callers can inspect the coefficient of each free unknown.
        result[column] = sub(reduced.get(row, n), &add(&terms));
    }
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    #[test]
    fn coefficients_and_parametric_solution() {
        let symbols = vec![Expr::symbol("x"), Expr::symbol("y")];
        let expression = parse_str("a*(x+2*y)-3").unwrap();
        let (a, b) = linear_eq_to_matrix(&[expression], &symbols).unwrap();
        assert_eq!(a.get(0, 0).to_string(), "a");
        assert_eq!(a.get(0, 1).to_string(), "2*a");
        let solution = linsolve(&a, &b, &symbols).unwrap().unwrap();
        assert_eq!(solution[1], symbols[1]);
        assert_eq!(crate::expand::expand(&sub(&mul2(&a.get(0, 0), &solution[0]), &sub(b.get(0, 0), &mul2(a.get(0, 1), &solution[1])))), Expr::zero());
    }

    #[test]
    fn nonlinear_cancellation_is_not_silently_performed() {
        let x = Expr::symbol("x");
        for source in ["x**2", "1/x", "sin(x)", "(x**2-1)/(x-1)"] {
            assert!(matches!(linear_eq_to_matrix(&[parse_str(source).unwrap()], &[x.clone()]), Err(LinearError::Nonlinear(_))));
        }
    }
}
