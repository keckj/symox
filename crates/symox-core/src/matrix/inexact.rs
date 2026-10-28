//! Dense Gauss-Jordan inversion over finite RR and CC scalar domains.
//! Values remain binary floats at the selected domain precision throughout.

use super::{Matrix, MatrixError};
use crate::bigfloat::BigFloat;
use crate::complex_parts::literal_parts;
use crate::evalf_mp::MpComplex;
use crate::expr::Expr;

fn scalar_parts(entry: &Expr) -> Option<(Expr, Expr)> {
    let parts = literal_parts(entry).or_else(|| crate::complex_parts::as_real_imag(entry, true, None).ok().flatten())?;
    for part in [&parts.0, &parts.1] {
        let number = part.as_number()?;
        if matches!(number, crate::number::Number::Float(value) if !value.is_finite()) { return None; }
    }
    Some(parts)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InexactDomain {
    pub complex: bool,
    pub precision: u32,
}

/// The same scalar domain admitted by construct_domain: rational real and
/// imaginary components, at least one Float, and no nonfinite components.
pub fn domain(matrix: &Matrix) -> Option<InexactDomain> {
    let mut precision = None;
    let mut complex = false;
    for entry in matrix.data() {
        let (re, im) = scalar_parts(entry)?;
        complex |= !im.is_zero();
        for part in [re, im] {
            let number = part.as_number()?;
            if let Some(p) = number.float_precision() {
                precision = Some(precision.unwrap_or(0).max(p));
            }
        }
    }
    Some(InexactDomain { complex, precision: precision? })
}

#[derive(Debug)]
pub enum InexactError<E> {
    Matrix(MatrixError),
    Interrupted(E),
}

fn is_zero(value: &MpComplex) -> bool { value.re.is_zero() && value.im.is_zero() }

/// mpmath's abs(mpc): exact squares, sum at p+4, then sqrt at p.
/// Comparing these rounded magnitudes also retains the first row on ties.
fn magnitude(value: &MpComplex, precision: u32) -> BigFloat {
    if value.im.is_zero() { return value.re.abs().round_to(precision); }
    if value.re.is_zero() { return value.im.abs().round_to(precision); }
    let square = |v: &BigFloat| v.mul(v, (2*v.man.bits() + 2) as u32);
    square(&value.re).add(&square(&value.im), precision + 4).sqrt(precision).expect("positive sum of squares")
}

/// ddm_irref([M | I]), including its partial pivot search and singularity
/// detection in the augmented identity columns. The callback is checked
/// during input conversion and every elimination row, preserving interrupts.
pub fn inverse<E>(matrix: &Matrix, selected: InexactDomain, mut check: impl FnMut() -> Result<(), E>) -> Result<Matrix, InexactError<E>> {
    let failure = |e| InexactError::Matrix(e);
    if matrix.rows() != matrix.cols() {
        return Err(failure(MatrixError::NonSquare("A Matrix must be square to invert.".into())));
    }
    if selected.precision == 0 || selected.precision > u32::MAX - 16 {
        return Err(failure(MatrixError::Value("invalid inexact matrix precision".into())));
    }
    let n = matrix.rows();
    let width = 2*n;
    let precision = selected.precision;
    // Float(..., dps) and RealField/ComplexField.to_sympy both pass through
    // prec_to_dps -> dps_to_prec. This presentation precision may differ
    // slightly from p (e.g. 139 -> 140, 211 -> 213). Arithmetic remains at p.
    let presentation = crate::bigfloat::dps_to_prec(crate::bigfloat::prec_to_dps(precision));
    let zero = MpComplex::real(BigFloat::zero());
    let one = MpComplex::real(BigFloat::one());
    let mut entries = vec![zero.clone(); n*width];
    for row in 0..n {
        check().map_err(InexactError::Interrupted)?;
        for col in 0..n {
            let (re, im) = scalar_parts(matrix.get(row, col)).ok_or_else(|| failure(MatrixError::NotImplemented("inexact inversion requires finite rational or Float scalar components".into())))?;
            if !selected.complex && !im.is_zero() {
                return Err(failure(MatrixError::Value("complex entry in a real inexact matrix".into())));
            }
            entries[row*width + col] = MpComplex {
                re: re.as_number().unwrap().to_bigfloat(presentation).round_to(precision),
                im: im.as_number().unwrap().to_bigfloat(presentation).round_to(precision),
            };
        }
        entries[row*width + n + row] = one.clone();
    }
    let mut row = 0;
    for col in 0..width {
        if row == n { break; }
        check().map_err(InexactError::Interrupted)?;
        let mut pivot = row;
        let mut largest = magnitude(&entries[row*width + col], precision);
        for candidate in row + 1..n {
            let size = magnitude(&entries[candidate*width + col], precision);
            if size.cmp(&largest).is_gt() { pivot = candidate; largest = size; }
        }
        if is_zero(&entries[pivot*width + col]) { continue; }
        if col >= n { return Err(failure(super::non_invertible())); }
        if row != pivot {
            for j in 0..width { entries.swap(row*width + j, pivot*width + j); }
        }
        // mpc**-1 rounds its exact squared norm at p+10, then each
        // division at p. MpComplex::div follows that same convention.
        let reciprocal = one.div(&entries[row*width + col], precision).expect("nonzero pivot");
        for j in col..width {
            entries[row*width + j] = entries[row*width + j].mul(&reciprocal, precision);
        }
        for other in 0..n {
            check().map_err(InexactError::Interrupted)?;
            if other == row || is_zero(&entries[other*width + col]) { continue; }
            let factor = entries[other*width + col].clone();
            entries[other*width + col] = zero.clone();
            for j in col + 1..width {
                entries[other*width + j] = entries[other*width + j].sub(&factor.mul(&entries[row*width + j], precision), precision);
            }
        }
        row += 1;
    }
    if row != n { return Err(failure(super::non_invertible())); }
    let result = (0..n).flat_map(|r| {
        let entries = &entries;
        (0..n).map(move |c| crate::evalf_mp::to_expr(entries[r*width + n + c].clone(), presentation))
    }).collect();
    Ok(Matrix::new_unchecked(n, n, result))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pivoted_real_and_complex_inverses_keep_float_coefficients() {
        for complex in [false, true] {
            let precision = 140;
            let value = |s| crate::parser::parse_str_with_precision(s, precision).unwrap();
            let entries = if complex { ["0.0", "2.0*I", "4.0", "0.0"] } else { ["0.0", "2.0", "4.0", "0.0"] };
            let matrix = Matrix::new(2,2,entries.map(value).to_vec()).unwrap();
            assert_eq!(domain(&matrix), Some(InexactDomain { complex, precision }));
            let result = inverse(&matrix, domain(&matrix).unwrap(), || Ok::<_, ()>(())).unwrap();
            assert_eq!(result.get(0,0), &Expr::zero());
            assert_eq!(result.get(1,1), &Expr::zero());
            assert_eq!(result.get(0,1), &value("0.25"));
            assert_eq!(result.get(1,0), &value(if complex { "-0.5*I" } else { "0.5" }));
        }
    }

    #[test]
    fn singularity_and_interruptions_remain_distinct() {
        let matrix = Matrix::new(2,2,[1,2,2,4].map(Expr::int).to_vec()).unwrap();
        let selected = InexactDomain { complex: false, precision: 140 };
        assert!(matches!(inverse(&matrix, selected, || Ok::<_, ()>(())), Err(InexactError::Matrix(MatrixError::NonInvertible(_)))));
        let mut calls = 0;
        let result = inverse(&matrix, selected, || { calls += 1; if calls == 4 { Err("interrupted") } else { Ok(()) } });
        assert!(matches!(result, Err(InexactError::Interrupted("interrupted"))));
        assert_eq!(calls, 4);
    }

    #[test]
    fn domain_proof_requires_a_float_and_finite_rational_parts() {
        for (source, expected) in [
            ("1/3", None),
            ("1+I", None),
            ("1.25", Some(InexactDomain { complex: false, precision: 140 })),
            ("1.25+I/3", Some(InexactDomain { complex: true, precision: 140 })),
            ("oo", None),
            ("nan", None),
            ("x+1.25", None),
            ("sqrt(2)+1.25*I", None),
        ] {
            let entry = crate::parser::parse_str_with_precision(source, 140).unwrap();
            assert_eq!(domain(&Matrix::new(1,1,vec![entry]).unwrap()), expected, "{source}");
        }
        for precision in 1..=3 {
            let value = Expr::num(crate::number::Number::mp_float(BigFloat::from_i64(2), precision));
            let matrix = Matrix::new(1,1,vec![value]).unwrap();
            let result = inverse(&matrix, domain(&matrix).unwrap(), || Ok::<_, ()>(())).unwrap();
            assert_eq!(result.get(0,0).as_number().unwrap().float_precision(), Some(7));
            assert_eq!(result.get(0,0).as_number().unwrap().to_bigfloat(10), BigFloat::one().ldexp(-1));
        }
    }
}
