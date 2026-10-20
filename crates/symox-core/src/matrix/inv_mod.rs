//! Certified inversion over Z/mZ, including composite and even moduli.
//!
//! Bezout row combinations create a unit pivot when no entry of a column is
//! individually a unit. All matrix entries are kept reduced modulo m.
use super::{MResult, Matrix, MatrixError};
use crate::bigfloat::BigFloat;
use crate::{Expr, Number};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

const MAX_STORAGE_BITS: u64 = 128_000_000;
const MAX_WORK: u64 = 100_000_000;
fn limited() -> MatrixError {
    MatrixError::NotImplemented("Exact modular matrix inverse resource limit exceeded".into())
}
fn domain_error() -> MatrixError {
    MatrixError::Value("inv_mod: matrix entries must be integers".into())
}
fn noninvertible() -> MatrixError {
    MatrixError::NonInvertible("Matrix is not invertible modulo the specified modulus".into())
}
fn words(bits: u64) -> u64 {
    bits.div_ceil(64).max(1)
}

struct Budget {
    work: u64,
    storage: u64,
    reserved: u64,
}
impl Budget {
    fn new() -> Self {
        Self {
            work: MAX_WORK,
            storage: MAX_STORAGE_BITS,
            reserved: 0,
        }
    }
    fn claim(&mut self, cost: u64) -> MResult<()> {
        self.work = self.work.checked_sub(cost).ok_or_else(limited)?;
        Ok(())
    }
    fn shape(&mut self, rows: usize, cols: usize, modulus: &BigInt) -> MResult<()> {
        if rows != cols {
            return Err(super::non_square());
        }
        if modulus <= &BigInt::zero() {
            return Err(MatrixError::Value(
                "modulus must be a positive integer".into(),
            ));
        }
        if rows == 0 {
            return Ok(());
        }
        let n = u64::try_from(rows).map_err(|_| limited())?;
        let cells = n.checked_mul(n).ok_or_else(limited)?;
        // Five copies bound original residues, [A|I], inverse, and output.
        // Modulo one only the shared-zero output is needed; certificates are
        // tautological in that ring and no cubic elimination is performed.
        self.reserved = if modulus.is_one() {
            cells.saturating_mul(64).saturating_add(512)
        } else {
            cells
                .saturating_mul(5)
                .saturating_mul(modulus.bits().saturating_add(256))
                .saturating_add(modulus.bits().saturating_mul(16))
        };
        if self.reserved > self.storage {
            return Err(limited());
        }
        let estimate = if modulus.is_one() {
            cells
        } else {
            n.saturating_pow(3)
                .saturating_mul(24)
                .saturating_mul(words(modulus.bits()).saturating_pow(2))
        };
        if estimate > self.work {
            return Err(limited());
        }
        Ok(())
    }
    fn scratch(&self, bits: u64) -> MResult<()> {
        if self.reserved.saturating_add(bits.saturating_mul(4)) > self.storage {
            Err(limited())
        } else {
            Ok(())
        }
    }
    fn reduce(&mut self, value: &BigInt, modulus: &BigInt) -> MResult<BigInt> {
        if modulus.is_one() {
            return Ok(BigInt::zero());
        }
        self.scratch(value.bits())?;
        self.claim(words(value.bits()).saturating_mul(words(modulus.bits())))?;
        if !value.is_negative() && value < modulus {
            Ok(value.clone())
        } else {
            Ok(value.mod_floor(modulus))
        }
    }
    fn mul(&mut self, a: &BigInt, b: &BigInt, modulus: &BigInt) -> MResult<BigInt> {
        let bits = a.bits().saturating_add(b.bits()).saturating_add(1);
        self.scratch(bits)?;
        self.claim(words(a.bits()).saturating_mul(words(b.bits())))?;
        let product = a * b;
        self.reduce(&product, modulus)
    }
    fn add(&mut self, a: &BigInt, b: &BigInt, modulus: &BigInt) -> MResult<BigInt> {
        self.claim(words(a.bits().max(b.bits())))?;
        let value = a + b;
        self.reduce(&value, modulus)
    }
    fn subtract(&mut self, a: &BigInt, b: &BigInt, modulus: &BigInt) -> MResult<BigInt> {
        self.claim(words(a.bits().max(b.bits())))?;
        let value = a - b;
        self.reduce(&value, modulus)
    }
    fn linear(
        &mut self,
        a: &BigInt,
        x: &BigInt,
        b: &BigInt,
        y: &BigInt,
        modulus: &BigInt,
    ) -> MResult<BigInt> {
        let left = self.mul(a, x, modulus)?;
        let right = self.mul(b, y, modulus)?;
        self.add(&left, &right, modulus)
    }
    /// Exact extended Euclid, with accounting before every division/product.
    fn bezout(&mut self, a: &BigInt, b: &BigInt) -> MResult<(BigInt, BigInt, BigInt)> {
        self.scratch(a.bits().max(b.bits()).saturating_mul(2).saturating_add(4))?;
        let (mut old_r, mut r) = (a.clone(), b.clone());
        let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());
        let (mut old_t, mut t) = (BigInt::zero(), BigInt::one());
        while !r.is_zero() {
            self.claim(words(old_r.bits()).saturating_mul(words(r.bits())))?;
            let q = &old_r / &r;
            self.claim(
                words(q.bits()).saturating_mul(
                    words(r.bits())
                        .saturating_add(words(s.bits()))
                        .saturating_add(words(t.bits())),
                ),
            )?;
            let next_r = &old_r - &q * &r;
            let next_s = &old_s - &q * &s;
            let next_t = &old_t - &q * &t;
            (old_r, r) = (r, next_r);
            (old_s, s) = (s, next_s);
            (old_t, t) = (t, next_t);
        }
        Ok((old_r, old_s, old_t))
    }
    fn unit_inverse(&mut self, value: &BigInt, modulus: &BigInt) -> MResult<Option<BigInt>> {
        if value.is_one() {
            return Ok(Some(BigInt::one()));
        }
        if value.is_zero() {
            return Ok(None);
        }
        let (g, coefficient, _) = self.bezout(value, modulus)?;
        if g.is_one() {
            Ok(Some(self.reduce(&coefficient, modulus)?))
        } else {
            Ok(None)
        }
    }
    fn power_of_two(&mut self, mut exponent: u64, modulus: &BigInt) -> MResult<BigInt> {
        let mut result = BigInt::one();
        let mut base = BigInt::from(2).mod_floor(modulus);
        while exponent > 0 {
            if exponent & 1 != 0 {
                result = self.mul(&result, &base, modulus)?;
            }
            exponent >>= 1;
            if exponent > 0 {
                base = self.mul(&base, &base, modulus)?;
            }
        }
        Ok(result)
    }
    fn float_residue(&mut self, number: &BigFloat, modulus: &BigInt) -> MResult<BigInt> {
        if number.man.is_zero() {
            return Ok(BigInt::zero());
        }
        if number.exp < 0 {
            self.scratch(number.man.bits())?;
            self.claim(words(number.man.bits()))?;
            let shift = u64::try_from(number.exp.unsigned_abs()).unwrap_or(u64::MAX);
            if number.man.trailing_zeros().unwrap_or(0) < shift {
                return Err(domain_error());
            }
            if modulus.is_one() {
                return Ok(BigInt::zero());
            }
            let shift = usize::try_from(shift).map_err(|_| limited())?;
            self.reduce(&(&number.man >> shift), modulus)
        } else {
            if modulus.is_one() {
                return Ok(BigInt::zero());
            }
            let mantissa = self.reduce(&number.man, modulus)?;
            let power = self.power_of_two(number.exp as u64, modulus)?;
            self.mul(&mantissa, &power, modulus)
        }
    }
    fn entry(&mut self, entry: &Expr, modulus: &BigInt) -> MResult<BigInt> {
        self.claim(1)?;
        // Held arithmetic is intentionally not coerced, matching SymPy's ZZ
        // conversion. Inspecting the tag is sufficient even for giant DAGs.
        match entry.as_number() {
            Some(Number::Int(value)) => self.reduce(&BigInt::from(*value), modulus),
            Some(Number::BigInt(value)) => self.reduce(value, modulus),
            Some(Number::Float(value)) if value.is_finite() && value.fract() == 0.0 => {
                self.float_residue(&BigFloat::from_f64(*value), modulus)
            }
            Some(Number::MpFloat(value, _)) => self.float_residue(value, modulus),
            _ => Err(domain_error()),
        }
    }
}

/// Shared by the dense call and the sparse wrapper before densification.
pub fn preflight(rows: usize, cols: usize, modulus: &BigInt) -> MResult<()> {
    Budget::new().shape(rows, cols, modulus)
}

impl Matrix {
    /// Exact inverse modulo any positive integer, without factoring it.
    /// Integer-valued floating entries are reduced from their dyadic form;
    /// their (possibly enormous) integer values are never materialized.
    pub fn inv_mod(&self, modulus: &BigInt) -> MResult<Matrix> {
        let n = self.rows();
        let mut budget = Budget::new();
        budget.shape(n, self.cols(), modulus)?;
        if n == 0 || modulus.is_one() {
            // Validate the ZZ domain without retaining an unnecessary dense
            // array of BigInt zeros in the trivial ring.
            for entry in self.data() {
                budget.entry(entry, modulus)?;
            }
            return Ok(Matrix::zeros(n, n));
        }
        let original = self
            .data()
            .iter()
            .map(|entry| budget.entry(entry, modulus))
            .collect::<MResult<Vec<_>>>()?;
        let width = n.checked_mul(2).ok_or_else(limited)?;
        let mut augmented = vec![BigInt::zero(); n.checked_mul(width).ok_or_else(limited)?];
        for i in 0..n {
            augmented[i * width..i * width + n].clone_from_slice(&original[i * n..(i + 1) * n]);
            augmented[i * width + n + i] = BigInt::one();
        }
        for col in 0..n {
            let pivot = (col..n)
                .find(|&i| !augmented[i * width + col].is_zero())
                .ok_or_else(noninvertible)?;
            if pivot != col {
                for j in 0..width {
                    augmented.swap(col * width + j, pivot * width + j);
                }
            }
            let mut inverse = budget.unit_inverse(&augmented[col * width + col], modulus)?;
            for row in col + 1..n {
                if inverse.is_some() {
                    break;
                }
                let a = augmented[col * width + col].clone();
                let b = augmented[row * width + col].clone();
                if b.is_zero() {
                    continue;
                }
                let (g, s, t) = budget.bezout(&a, &b)?;
                budget.claim(
                    words(a.bits().max(b.bits()))
                        .saturating_mul(words(g.bits()))
                        .saturating_mul(2),
                )?;
                let lower_left = -(&b / &g);
                let lower_right = &a / &g;
                for j in 0..width {
                    let x = augmented[col * width + j].clone();
                    let y = augmented[row * width + j].clone();
                    let upper = budget.linear(&s, &x, &t, &y, modulus)?;
                    let lower = budget.linear(&lower_left, &x, &lower_right, &y, modulus)?;
                    augmented[col * width + j] = upper;
                    augmented[row * width + j] = lower;
                }
                inverse = budget.unit_inverse(&augmented[col * width + col], modulus)?;
            }
            let inverse = inverse.ok_or_else(noninvertible)?;
            for j in 0..width {
                augmented[col * width + j] =
                    budget.mul(&augmented[col * width + j], &inverse, modulus)?;
            }
            for row in 0..n {
                if row == col {
                    continue;
                }
                let factor = augmented[row * width + col].clone();
                if factor.is_zero() {
                    continue;
                }
                for j in 0..width {
                    let product = budget.mul(&factor, &augmented[col * width + j], modulus)?;
                    augmented[row * width + j] =
                        budget.subtract(&augmented[row * width + j], &product, modulus)?;
                }
            }
        }
        let inverse: Vec<_> = (0..n)
            .flat_map(|i| augmented[i * width + n..(i + 1) * width].iter().cloned())
            .collect();
        for (a, b) in [(&original, &inverse), (&inverse, &original)] {
            for i in 0..n {
                for j in 0..n {
                    let mut value = BigInt::zero();
                    for k in 0..n {
                        if !a[i * n + k].is_zero() && !b[k * n + j].is_zero() {
                            let product = budget.mul(&a[i * n + k], &b[k * n + j], modulus)?;
                            value = budget.add(&value, &product, modulus)?;
                        }
                    }
                    let expected = if i == j {
                        BigInt::one()
                    } else {
                        BigInt::zero()
                    };
                    if value != expected {
                        return Err(MatrixError::NotImplemented(
                            "Unable to certify the modular matrix inverse".into(),
                        ));
                    }
                }
            }
        }
        Matrix::new(
            n,
            n,
            inverse
                .into_iter()
                .map(|value| Expr::num(Number::from_bigint(value)))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn matrix(data: &[&[i64]]) -> Matrix {
        Matrix::from_rows(
            data.iter()
                .map(|row| row.iter().map(|n| Expr::int(*n)).collect())
                .collect(),
        )
        .unwrap()
    }
    #[test]
    fn composite_bezout_and_small_exhaustive() {
        for modulus in [4, 6, 9] {
            let m = BigInt::from(modulus);
            for a in 0..4 {
                for b in 0..4 {
                    for c in 0..4 {
                        for d in 0..4 {
                            let input = matrix(&[&[a, b], &[c, d]]);
                            let det = BigInt::from(a * d - b * c);
                            let result = input.inv_mod(&m);
                            if det.gcd(&m).is_one() {
                                assert!(result.is_ok(), "{input}: {result:?}");
                            } else {
                                assert!(matches!(result, Err(MatrixError::NonInvertible(_))));
                            }
                        }
                    }
                }
            }
        }
        assert!(matrix(&[&[6, 1, 0], &[10, 0, 1], &[15, 1, 1]])
            .inv_mod(&BigInt::from(30))
            .is_ok());
        assert_eq!(
            matrix(&[&[2, 3], &[3, 2]])
                .inv_mod(&BigInt::from(6))
                .unwrap(),
            matrix(&[&[2, 3], &[3, 2]])
        );
    }
    #[test]
    fn compact_float_and_trivial_ring() {
        let value = Expr::num(Number::MpFloat(
            BigFloat {
                man: BigInt::one(),
                exp: 1_000_000_000,
            },
            80,
        ));
        let input = Matrix::new(2, 2, vec![value, Expr::one(), Expr::one(), Expr::zero()]).unwrap();
        let inverse = input.inv_mod(&BigInt::from(127)).unwrap();
        let residue = BigInt::from(2).modpow(&BigInt::from(1_000_000_000), &BigInt::from(127));
        assert_eq!(
            inverse,
            matrix(&[&[0, 1], &[1, 0]])
                .add(&Matrix::diag(&[
                    Expr::zero(),
                    Expr::num(Number::from_bigint(
                        (-residue).mod_floor(&BigInt::from(127))
                    ))
                ]))
                .unwrap()
        );
        let extreme = Expr::num(Number::MpFloat(
            BigFloat {
                man: BigInt::one(),
                exp: i64::MAX.into(),
            },
            80,
        ));
        assert_eq!(
            Matrix::new(1, 1, vec![extreme])
                .unwrap()
                .inv_mod(&BigInt::from(127))
                .unwrap(),
            Matrix::eye(1)
        );
        let integer_float = Expr::num(Number::MpFloat(
            BigFloat {
                man: BigInt::from(-4),
                exp: -1,
            },
            80,
        ));
        assert_eq!(
            Matrix::new(1, 1, vec![integer_float])
                .unwrap()
                .inv_mod(&BigInt::from(5))
                .unwrap(),
            matrix(&[&[2]])
        );
        assert_eq!(input.inv_mod(&BigInt::one()).unwrap(), Matrix::zeros(2, 2));
        assert_eq!(
            Matrix::zeros(0, 0).inv_mod(&BigInt::from(6)).unwrap(),
            Matrix::zeros(0, 0)
        );
    }
    #[test]
    fn domain_and_resource_preflight_without_materializing() {
        let held = Expr::from_unevaluated_kind(crate::expr::Kind::Pow(
            Expr::int(2),
            Expr::int(1_000_000_000),
        ));
        assert!(matches!(
            Matrix::new(1, 1, vec![held])
                .unwrap()
                .inv_mod(&BigInt::from(5)),
            Err(MatrixError::Value(_))
        ));
        let mut budget = Budget::new();
        budget.storage = 100;
        assert!(budget.shape(2, 2, &BigInt::from(6)).is_err());
        let mut budget = Budget::new();
        budget.work = 10;
        assert!(budget.shape(2, 2, &BigInt::from(6)).is_err());
        let mut budget = Budget::new();
        budget.work = 1;
        assert!(budget
            .mul(&BigInt::from(123), &BigInt::from(456), &BigInt::from(1009))
            .is_err());
        let mut budget = Budget::new();
        budget.work = 1;
        let wide_float = BigFloat {
            man: BigInt::one() << 128,
            exp: -128,
        };
        assert!(budget.float_residue(&wide_float, &BigInt::one()).is_err());
        assert!(preflight(1_000_000, 1_000_000, &BigInt::from(6)).is_err());
        assert!(Matrix::new(
            1,
            1,
            vec![Expr::num(Number::MpFloat(
                BigFloat {
                    man: BigInt::one(),
                    exp: -1
                },
                80
            ))]
        )
        .unwrap()
        .inv_mod(&BigInt::from(5))
        .is_err());
    }
}
