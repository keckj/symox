//! Rodrigues' identity for numeric real 2x2 and 3x3 skew matrices.
//!
//! Float entries are lifted to their exact dyadic values before products,
//! avoiding underflow in the squared angle. Only the final result is rounded
//! to the input precision by the native multiprecision evaluator.
use super::budget::Budget;
use crate::expr::Expr;
use crate::functions::{cos, sin};
use crate::matrix::{MResult, Matrix, MatrixError};
use crate::number::Number;
use crate::ops::{add, div, mul2, pow, sqrt};

pub(super) fn evaluate(matrix: &Matrix) -> MResult<Option<Matrix>> {
    if !matches!(matrix.shape(), (2, 2) | (3, 3)) {
        return Ok(None);
    }
    if !matrix
        .data()
        .iter()
        .all(|entry| entry.as_number().is_some())
        || !matrix.data().iter().any(|entry| {
            entry
                .as_number()
                .and_then(Number::float_precision)
                .is_some()
        })
    {
        return Ok(None);
    }
    // Reject other matrices before dyadic allocation or new resource checks;
    // their existing exact/spectral paths keep their accepted domains.
    for i in 0..matrix.rows() {
        if !matrix.get(i, i).as_number().unwrap().is_zero() {
            return Ok(None);
        }
        for j in i + 1..matrix.cols() {
            if matrix
                .get(i, j)
                .as_number()
                .unwrap()
                .cmp_value(&matrix.get(j, i).as_number().unwrap().neg())
                != Some(std::cmp::Ordering::Equal)
            {
                return Ok(None);
            }
        }
    }
    let mut budget = Budget::new();
    budget.matrix(matrix)?;
    let mut precision = None;
    let mut entries = Vec::with_capacity(matrix.data().len());
    for entry in matrix.data() {
        let Some(number) = entry.as_number() else {
            return Ok(None);
        };
        let exact = match number {
            Number::Float(value) if value.is_finite() => {
                precision = Some(precision.unwrap_or(0).max(53));
                crate::bigfloat::BigFloat::from_f64(*value).to_ratio_exact()
            }
            Number::MpFloat(value, bits) => {
                // Check before allocating a dyadic denominator or numerator.
                if value.exp.unsigned_abs().saturating_add(u128::from(value.man.bits())) > 1_000_000 {
                    return Err(MatrixError::NotImplemented(
                        "Exact matrix function resource limit exceeded".into(),
                    ));
                }
                precision = Some(precision.unwrap_or(0).max(*bits));
                value.to_ratio_exact()
            }
            _ => match number.to_ratio() {
                Some(ratio) => ratio,
                None => return Ok(None),
            },
        };
        entries.push(Expr::num(Number::from_ratio(exact)));
    }
    let n = matrix.rows();
    let exact = Matrix::new(n, n, entries)?;
    budget.matrix(&exact)?;
    for i in 0..n {
        if !exact.get(i, i).is_zero() {
            return Ok(None);
        }
        for j in i + 1..n {
            if !add(&[exact.get(i, j).clone(), exact.get(j, i).clone()]).is_zero() {
                return Ok(None);
            }
        }
    }
    let mut squares = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            squares.push(budget.scalar_power(exact.get(i, j), &Expr::int(2))?);
        }
    }
    let theta = sqrt(&add(&squares));
    budget.expression(&theta)?;
    let sinc = |argument: Expr| {
        if argument.is_zero() {
            Expr::one()
        } else {
            div(&sin(&argument), &argument)
        }
    };
    let identity = budget.identity(n)?;
    let linear = budget.scale(&exact, &sinc(theta.clone()))?;
    let result = if n == 2 {
        // A^2 = -theta^2 I.
        let constant = budget.scale(&identity, &cos(&theta))?;
        budget.add(&constant, &linear)?
    } else {
        // A^3 = -theta^2 A. The half-angle coefficient avoids subtracting
        // nearly equal cosines. The zero angle is handled before division.
        let half_sinc = sinc(mul2(&Expr::rational(1, 2), &theta));
        let coefficient = mul2(&Expr::rational(1, 2), &pow(&half_sinc, &Expr::int(2)));
        let square = budget.mul(&exact, &exact)?;
        let quadratic = budget.scale(&square, &coefficient)?;
        let first = budget.add(&identity, &linear)?;
        budget.add(&first, &quadratic)?
    };
    let Some(precision) = precision else {
        return Ok(Some(result));
    };
    // Large rotations need extra absolute accuracy in sqrt(sum(a_ij^2))
    // before trigonometric argument reduction. Retain guard bits also for
    // the quadratic products and cancellation in the diagonal entries.
    let magnitude = exact
        .data()
        .iter()
        .filter_map(Expr::as_number)
        .filter(|number| !number.is_zero())
        .map(|number| u64::try_from(number.to_bigfloat(64).mag().unsigned_abs()).unwrap_or(u64::MAX))
        .max()
        .unwrap_or(0);
    let working = 2 * u64::from(precision)
        .saturating_add(2 * magnitude)
        .saturating_add(64);
    if working > 1_000_000 {
        return Err(MatrixError::NotImplemented(
            "Exact matrix function resource limit exceeded".into(),
        ));
    }
    budget.work(working)?;
    let mut output = Vec::with_capacity(n * n);
    for entry in result.data() {
        let mut working = working;
        let rounded = loop {
            if working.saturating_mul(2) > 1_000_000 {
                return Err(MatrixError::NotImplemented(
                    "Exact matrix function resource limit exceeded".into(),
                ));
            }
            budget.work(3 * working)?;
            let evaluate = |bits: u64| {
                crate::evalf_mp::eval_mp(entry, bits as u32).ok_or_else(|| {
                    MatrixError::NotImplemented(
                        "Cannot evaluate the skew matrix exponential at the requested precision"
                            .into(),
                    )
                })
            };
            let value = evaluate(working)?.round_to(precision);
            let refined = evaluate(2 * working)?.round_to(precision);
            if value.re == refined.re
                && value.im == refined.im
                && (entry.is_zero() || !refined.re.is_zero() || !refined.im.is_zero())
            {
                break refined;
            }
            working *= 2;
        };
        output.push(crate::evalf_mp::to_expr(rounded, precision));
    }
    Ok(Some(Matrix::new(n, n, output)?))
}
