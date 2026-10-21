//! Exact deterministic row-basis LLL, with an integer inverse certificate.
//!
//! Reduction order and nearest-integer ties follow SymPy's DDM algorithm:
//! floor(mu + 1/2), descending size reductions, and adjacent Lovasz swaps.
use super::Matrix;
use crate::bigfloat::BigFloat;
use crate::{Expr, Number};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LllError {
    Delta,
    Parameter,
    Shape,
    Domain,
    Rank,
    Resource,
    Certificate,
}
impl fmt::Display for LllError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Delta => "delta must lie in range (0.25, 1)",
            Self::Parameter => "LLL delta must be a finite rational number",
            Self::Shape => "input matrix must have shape (m, n) with m <= n",
            Self::Domain => {
                "input matrix entries must belong to ZZ (exact integers or integral floats)"
            }
            Self::Rank => "input matrix contains linearly dependent rows",
            Self::Resource => "Exact matrix LLL resource limit exceeded",
            Self::Certificate => "Unable to certify exact LLL reduction and unimodularity",
        })
    }
}
impl std::error::Error for LllError {}
pub type LllResult<T> = Result<T, LllError>;

const MAX_BITS: u64 = 1_000_000;
const MAX_STORAGE: u64 = 128_000_000;
const MAX_WORK: u64 = 100_000_000;
fn words(bits: u64) -> u64 {
    bits.div_ceil(64).max(1)
}
fn qbits(q: &Q) -> u64 {
    q.numer().bits().saturating_add(q.denom().bits())
}

// One budget covers conversion, reduction, the independent Gram-Schmidt
// certificate, and both integer-inverse certificates. All raw products and
// normalization/GCD scratch are bounded *before* entering BigInt/BigRational.
struct Budget {
    work: u64,
    slots: u64,
    height: u64,
    output_storage: u64,
}
impl Budget {
    fn new() -> Self {
        Self {
            work: MAX_WORK,
            slots: 32,
            height: 64,
            output_storage: 0,
        }
    }
    fn claim(&mut self, work: u64) -> LllResult<()> {
        self.work = self.work.checked_sub(work).ok_or(LllError::Resource)?;
        Ok(())
    }
    fn scratch(&self, bits: u64) -> LllResult<()> {
        if bits > MAX_BITS
            || self
                .slots
                .saturating_mul(self.height.saturating_add(256))
                .saturating_add(self.output_storage)
                .saturating_add(bits.saturating_mul(16))
                > MAX_STORAGE
        {
            Err(LllError::Resource)
        } else {
            Ok(())
        }
    }
    fn record(&mut self, bits: u64) -> LllResult<()> {
        self.height = self.height.max(bits);
        self.scratch(0)
    }
    fn shape(&mut self, m: usize, n: usize) -> LllResult<()> {
        if m > n {
            return Err(LllError::Shape);
        }
        if m == 0 {
            return Ok(());
        }
        let m = m as u64;
        let n = n as u64;
        self.slots = m
            .saturating_mul(n)
            .saturating_mul(6)
            .saturating_add(m.saturating_mul(m).saturating_mul(6))
            .saturating_add(m.saturating_mul(6))
            .saturating_add(32);
        // Reserve actual native Expr/Arc headers before constructing either
        // output; BigInt payloads are covered by the scalar-cell reserve.
        let header = (std::mem::size_of::<crate::expr::Node>() + 4 * std::mem::size_of::<usize>())
            as u64
            * 8;
        self.output_storage = m.saturating_mul(n.saturating_add(m)).saturating_mul(header);
        self.scratch(0)?;
        self.claim(self.slots)?;
        let work = m
            .saturating_mul(m)
            .saturating_mul(n)
            .saturating_add(m.saturating_pow(3).saturating_mul(3));
        if work > self.work {
            return Err(LllError::Resource);
        }
        Ok(())
    }
    fn clone_int(&mut self, value: &BigInt) -> LllResult<BigInt> {
        self.scratch(value.bits())?;
        self.record(value.bits())?;
        self.claim(words(value.bits()))?;
        Ok(value.clone())
    }
    fn clone_q(&mut self, value: &Q) -> LllResult<Q> {
        self.scratch(qbits(value))?;
        self.claim(words(qbits(value)))?;
        Ok(value.clone())
    }
    fn clone_rows(&mut self, rows: &[Vec<BigInt>]) -> LllResult<Vec<Vec<BigInt>>> {
        rows.iter()
            .map(|row| row.iter().map(|v| self.clone_int(v)).collect())
            .collect()
    }
    fn qabs(&mut self, value: &Q) -> LllResult<Q> {
        self.scratch(qbits(value))?;
        self.claim(words(qbits(value)))?;
        Ok(value.abs())
    }
    fn add(&mut self, a: &BigInt, b: &BigInt) -> LllResult<BigInt> {
        self.scratch(a.bits().max(b.bits()).saturating_add(1))?;
        self.claim(words(a.bits().max(b.bits())))?;
        let value = a + b;
        self.record(value.bits())?;
        Ok(value)
    }
    fn sub(&mut self, a: &BigInt, b: &BigInt) -> LllResult<BigInt> {
        self.scratch(a.bits().max(b.bits()).saturating_add(1))?;
        self.claim(words(a.bits().max(b.bits())))?;
        let value = a - b;
        self.record(value.bits())?;
        Ok(value)
    }
    fn mul(&mut self, a: &BigInt, b: &BigInt) -> LllResult<BigInt> {
        if a.is_zero() || b.is_zero() {
            self.claim(1)?;
            return Ok(BigInt::zero());
        }
        self.scratch(a.bits().saturating_add(b.bits()))?;
        self.claim(words(a.bits()).saturating_mul(words(b.bits())))?;
        let value = a * b;
        self.record(value.bits())?;
        Ok(value)
    }
    fn quotient(&mut self, a: &BigInt, b: &BigInt) -> LllResult<BigInt> {
        self.scratch(a.bits().saturating_add(1))?;
        self.claim(
            words(a.bits())
                .saturating_mul(words(b.bits()))
                .saturating_mul(64),
        )?;
        let value = a.div_floor(b);
        self.record(value.bits())?;
        Ok(value)
    }
    fn rational(&mut self, n: BigInt, d: BigInt) -> LllResult<Q> {
        if d.is_zero() {
            return Err(LllError::Rank);
        }
        let bits = n.bits().saturating_add(d.bits());
        self.scratch(bits)?;
        self.claim(words(bits).saturating_pow(2).saturating_mul(64))?;
        let value = Q::new(n, d);
        self.record(qbits(&value))?;
        Ok(value)
    }
    fn integer(&mut self, a: &BigInt) -> LllResult<Q> {
        Ok(Q::from_integer(self.clone_int(a)?))
    }
    fn qadd(&mut self, a: &Q, b: &Q) -> LllResult<Q> {
        let left = self.mul(a.numer(), b.denom())?;
        let right = self.mul(b.numer(), a.denom())?;
        let n = self.add(&left, &right)?;
        let d = self.mul(a.denom(), b.denom())?;
        self.rational(n, d)
    }
    fn qsub(&mut self, a: &Q, b: &Q) -> LllResult<Q> {
        let left = self.mul(a.numer(), b.denom())?;
        let right = self.mul(b.numer(), a.denom())?;
        let n = self.sub(&left, &right)?;
        let d = self.mul(a.denom(), b.denom())?;
        self.rational(n, d)
    }
    fn qmul(&mut self, a: &Q, b: &Q) -> LllResult<Q> {
        let n = self.mul(a.numer(), b.numer())?;
        let d = self.mul(a.denom(), b.denom())?;
        self.rational(n, d)
    }
    fn qdiv(&mut self, a: &Q, b: &Q) -> LllResult<Q> {
        if b.is_zero() {
            return Err(LllError::Rank);
        }
        let n = self.mul(a.numer(), b.denom())?;
        let d = self.mul(a.denom(), b.numer())?;
        self.rational(n, d)
    }
    fn cmp(&mut self, a: &Q, b: &Q) -> LllResult<Ordering> {
        let left = self.mul(a.numer(), b.denom())?;
        let right = self.mul(b.numer(), a.denom())?;
        self.claim(words(left.bits().max(right.bits())))?;
        Ok(left.cmp(&right))
    }
    fn small(&mut self, a: &Q) -> LllResult<bool> {
        self.scratch(a.numer().bits())?;
        self.claim(words(a.numer().bits()))?;
        let twice = self.mul(&a.numer().abs(), &BigInt::from(2))?;
        Ok(twice <= *a.denom())
    }
    fn closest(&mut self, a: &Q) -> LllResult<BigInt> {
        let shifted = self.qadd(a, &Q::new(BigInt::one(), BigInt::from(2)))?;
        self.quotient(shifted.numer(), shifted.denom())
    }
    fn float_integer(&mut self, value: &BigFloat) -> LllResult<BigInt> {
        self.scratch(value.man.bits())?;
        self.claim(words(value.man.bits()))?;
        if value.man.is_zero() {
            return Ok(BigInt::zero());
        }
        if value.exp >= 0 {
            let bits = value.man.bits().saturating_add(u64::try_from(value.exp).unwrap_or(u64::MAX));
            self.scratch(bits)?;
            self.record(bits)?;
            self.claim(words(bits))?;
            let shift = usize::try_from(value.exp).map_err(|_| LllError::Resource)?;
            Ok(&value.man << shift)
        } else {
            let shift = u64::try_from(value.exp.unsigned_abs()).unwrap_or(u64::MAX);
            if value.man.trailing_zeros().unwrap_or(0) < shift {
                return Err(LllError::Domain);
            }
            let shift = usize::try_from(shift).map_err(|_| LllError::Resource)?;
            let result = &value.man >> shift;
            self.record(result.bits())?;
            Ok(result)
        }
    }
    fn entry(&mut self, entry: &Expr) -> LllResult<BigInt> {
        self.claim(1)?;
        match entry.as_number() {
            Some(Number::Int(n)) => Ok(BigInt::from(*n)),
            Some(Number::BigInt(n)) => self.clone_int(n),
            Some(Number::Float(f)) if f.is_finite() && f.fract() == 0.0 => {
                self.float_integer(&BigFloat::from_f64(*f))
            }
            Some(Number::MpFloat(f, _)) => self.float_integer(f),
            _ => Err(LllError::Domain),
        }
    }
    fn delta(&mut self, delta: &Expr) -> LllResult<Q> {
        let q = match delta.as_number() {
            Some(Number::Int(n)) => Q::from_integer(BigInt::from(*n)),
            Some(Number::BigInt(n)) => self.integer(n)?,
            Some(Number::Rational(q)) => {
                self.scratch(qbits(q))?;
                self.record(qbits(q))?;
                self.claim(words(qbits(q)))?;
                q.clone()
            }
            Some(Number::Float(f)) if f.is_finite() => self.float_delta(&BigFloat::from_f64(*f))?,
            Some(Number::MpFloat(f, _)) => self.float_delta(f)?,
            _ => return Err(LllError::Parameter),
        };
        if self.cmp(&q, &Q::new(BigInt::one(), BigInt::from(4)))? != Ordering::Greater
            || self.cmp(&q, &Q::one())? != Ordering::Less
        {
            return Err(LllError::Delta);
        }
        Ok(q)
    }
    fn float_delta(&mut self, f: &BigFloat) -> LllResult<Q> {
        // Only values in (1/4,1) can be valid; reject compact extremes before
        // materializing a dyadic denominator or shifting the mantissa.
        let bits = f.man.bits();
        self.scratch(bits)?;
        self.claim(words(bits))?;
        let scale = f.exp + i128::from(bits);
        if !f.man.is_positive() || !(-1..=0).contains(&scale) {
            return Err(LllError::Delta);
        }
        let shift = u64::try_from(f.exp.unsigned_abs()).unwrap_or(u64::MAX);
        self.scratch(bits.saturating_add(shift))?;
        self.claim(words(shift))?;
        let d = BigInt::one() << usize::try_from(shift).map_err(|_| LllError::Resource)?;
        let n = self.clone_int(&f.man)?;
        let q = self.rational(n, d)?;
        self.limit_denominator(&q)
    }
    fn limit_denominator(&mut self, q: &Q) -> LllResult<Q> {
        // SymPy QQ.from_sympy(Float) uses RR's default 53-bit max denominator.
        let max = BigInt::from(45_035_996_273_704u64);
        if q.denom() <= &max {
            return self.clone_q(q);
        }
        let (mut p0, mut q0, mut p1, mut q1) =
            (BigInt::zero(), BigInt::one(), BigInt::one(), BigInt::zero());
        let (mut n, mut d) = (self.clone_int(q.numer())?, self.clone_int(q.denom())?);
        loop {
            let a = self.quotient(&n, &d)?;
            let product = self.mul(&a, &q1)?;
            let q2 = self.add(&q0, &product)?;
            if q2 > max {
                break;
            }
            let product = self.mul(&a, &p1)?;
            let p2 = self.add(&p0, &product)?;
            (p0, q0, p1, q1) = (p1, q1, p2, q2);
            let product = self.mul(&a, &d)?;
            let remainder = self.sub(&n, &product)?;
            (n, d) = (d, remainder);
        }
        let difference = self.sub(&max, &q0)?;
        let k = self.quotient(&difference, &q1)?;
        let product = self.mul(&k, &p1)?;
        let p = self.add(&p0, &product)?;
        let product = self.mul(&k, &q1)?;
        let d = self.add(&q0, &product)?;
        let bound1 = self.rational(p, d)?;
        let bound2 = self.rational(p1, q1)?;
        if bound1.is_zero() || bound2.is_zero() {
            return self.clone_q(q);
        }
        let difference1 = self.qsub(&bound1, q)?;
        let difference2 = self.qsub(&bound2, q)?;
        let difference1 = self.qabs(&difference1)?;
        let difference2 = self.qabs(&difference2)?;
        if self.cmp(&difference2, &difference1)? != Ordering::Greater {
            Ok(bound2)
        } else {
            Ok(bound1)
        }
    }
}

fn gram_schmidt(
    y: &[Vec<BigInt>],
    n: usize,
    budget: &mut Budget,
) -> LllResult<(Vec<Vec<Q>>, Vec<Q>)> {
    let m = y.len();
    let mut star = vec![vec![Q::zero(); n]; m];
    let mut mu = vec![vec![Q::zero(); m]; m];
    let mut g = vec![Q::zero(); m];
    for i in 0..m {
        for j in 0..n {
            star[i][j] = budget.integer(&y[i][j])?;
        }
        for j in 0..i {
            let mut dot = Q::zero();
            for z in 0..n {
                let integer = budget.integer(&y[i][z])?;
                let term = budget.qmul(&integer, &star[j][z])?;
                dot = budget.qadd(&dot, &term)?;
            }
            mu[i][j] = budget.qdiv(&dot, &g[j])?;
            for z in 0..n {
                let term = budget.qmul(&mu[i][j], &star[j][z])?;
                star[i][z] = budget.qsub(&star[i][z], &term)?;
            }
        }
        for z in 0..n {
            let term = budget.qmul(&star[i][z], &star[i][z])?;
            g[i] = budget.qadd(&g[i], &term)?;
        }
        // This additionally rejects SymPy's accidental acceptance of one
        // all-zero row: a row lattice basis must be linearly independent.
        if g[i].is_zero() {
            return Err(LllError::Rank);
        }
    }
    Ok((mu, g))
}
fn lovasz(mu: &[Vec<Q>], g: &[Q], delta: &Q, k: usize, b: &mut Budget) -> LllResult<bool> {
    let square = b.qmul(&mu[k][k - 1], &mu[k][k - 1])?;
    let difference = b.qsub(delta, &square)?;
    let rhs = b.qmul(&difference, &g[k - 1])?;
    Ok(b.cmp(&g[k], &rhs)? != Ordering::Less)
}
struct Basis {
    y: Vec<Vec<BigInt>>,
    t: Vec<Vec<BigInt>>,
    inverse: Vec<Vec<BigInt>>,
    mu: Vec<Vec<Q>>,
    g: Vec<Q>,
}
impl Basis {
    fn reduce(&mut self, k: usize, j: usize, b: &mut Budget) -> LllResult<()> {
        let r = b.closest(&self.mu[k][j])?;
        for z in 0..self.y[k].len() {
            let product = b.mul(&r, &self.y[j][z])?;
            self.y[k][z] = b.sub(&self.y[k][z], &product)?;
        }
        let rational = b.integer(&r)?;
        for z in 0..j {
            let product = b.qmul(&rational, &self.mu[j][z])?;
            self.mu[k][z] = b.qsub(&self.mu[k][z], &product)?;
        }
        self.mu[k][j] = b.qsub(&self.mu[k][j], &rational)?;
        for z in 0..self.t.len() {
            let product = b.mul(&r, &self.t[j][z])?;
            self.t[k][z] = b.sub(&self.t[k][z], &product)?;
            // T <- E*T, inverse <- inverse*E^-1 (column j += r column k).
            let product = b.mul(&r, &self.inverse[z][k])?;
            self.inverse[z][j] = b.add(&self.inverse[z][j], &product)?;
        }
        Ok(())
    }
    fn swap(&mut self, k: usize, b: &mut Budget) -> LllResult<()> {
        let nu = b.clone_q(&self.mu[k][k - 1])?;
        let square = b.qmul(&nu, &nu)?;
        let product = b.qmul(&square, &self.g[k - 1])?;
        let alpha = b.qadd(&self.g[k], &product)?;
        let beta = b.qdiv(&self.g[k - 1], &alpha)?;
        self.mu[k][k - 1] = b.qmul(&nu, &beta)?;
        self.g[k] = b.qmul(&self.g[k], &beta)?;
        self.g[k - 1] = alpha;
        self.y.swap(k, k - 1);
        self.t.swap(k, k - 1);
        for row in &mut self.inverse {
            row.swap(k, k - 1);
        }
        b.claim((k - 1) as u64)?;
        let (lower, upper) = self.mu.split_at_mut(k);
        for j in 0..k - 1 {
            std::mem::swap(&mut lower[k - 1][j], &mut upper[0][j]);
        }
        for i in k + 1..self.mu.len() {
            let xi = b.clone_q(&self.mu[i][k])?;
            let product = b.qmul(&nu, &xi)?;
            self.mu[i][k] = b.qsub(&self.mu[i][k - 1], &product)?;
            let product = b.qmul(&self.mu[k][k - 1], &self.mu[i][k])?;
            self.mu[i][k - 1] = b.qadd(&product, &xi)?;
        }
        Ok(())
    }
    fn certify(
        &self,
        original: &[Vec<BigInt>],
        n: usize,
        delta: &Q,
        b: &mut Budget,
    ) -> LllResult<()> {
        let m = self.y.len();
        let (mu, g) = gram_schmidt(&self.y, n, b)?;
        for i in 0..m {
            for j in 0..i {
                if !b.small(&mu[i][j])? {
                    return Err(LllError::Certificate);
                }
            }
            if i > 0 && !lovasz(&mu, &g, delta, i, b)? {
                return Err(LllError::Certificate);
            }
        }
        for i in 0..m {
            for j in 0..n {
                let mut sum = BigInt::zero();
                for (k, row) in original.iter().enumerate() {
                    let product = b.mul(&self.t[i][k], &row[j])?;
                    sum = b.add(&sum, &product)?;
                }
                if sum != self.y[i][j] {
                    return Err(LllError::Certificate);
                }
            }
        }
        // An integer two-sided inverse proves det(T)=+/-1, beyond T*A=B.
        for (left, right) in [(&self.t, &self.inverse), (&self.inverse, &self.t)] {
            for (i, row) in left.iter().enumerate() {
                for j in 0..m {
                    let mut sum = BigInt::zero();
                    for k in 0..m {
                        let product = b.mul(&row[k], &right[k][j])?;
                        sum = b.add(&sum, &product)?;
                    }
                    let expected = if i == j {
                        BigInt::one()
                    } else {
                        BigInt::zero()
                    };
                    if sum != expected {
                        return Err(LllError::Certificate);
                    }
                }
            }
        }
        Ok(())
    }
}

fn output(
    rows: usize,
    cols: usize,
    data: Vec<Vec<BigInt>>,
    budget: &mut Budget,
) -> LllResult<Matrix> {
    let cells = rows.checked_mul(cols).ok_or(LllError::Resource)?;
    let mut entries = Vec::with_capacity(cells);
    for value in data.into_iter().flatten() {
        budget.claim(words(value.bits()).saturating_mul(4))?;
        entries.push(Expr::num(Number::from_bigint(value)));
    }
    Matrix::new(rows, cols, entries).map_err(|_| LllError::Certificate)
}

/// Sparse wrappers use this before any dense allocation or entry conversion.
pub fn preflight(rows: usize, cols: usize, delta: &Expr) -> LllResult<()> {
    let mut b = Budget::new();
    b.delta(delta)?;
    b.shape(rows, cols)
}
impl Matrix {
    pub fn lll(&self, delta: &Expr) -> LllResult<Matrix> {
        Ok(self.lll_transform(delta)?.0)
    }
    pub fn lll_transform(&self, delta: &Expr) -> LllResult<(Matrix, Matrix)> {
        let (m, n) = self.shape();
        let mut b = Budget::new();
        let delta = b.delta(delta)?;
        b.shape(m, n)?;
        if m == 0 {
            return Ok((self.clone(), Matrix::zeros(0, 0)));
        }
        let mut original = Vec::with_capacity(m);
        for row in self.data().chunks(n) {
            original.push(
                row.iter()
                    .map(|e| b.entry(e))
                    .collect::<LllResult<Vec<_>>>()?,
            );
        }
        if m == 1 {
            if original[0].iter().all(BigInt::is_zero) {
                return Err(LllError::Rank);
            }
            // There are no size/Lovasz inequalities for a singleton basis.
            // The identity transform and its inverse certify the lattice
            // without squaring an otherwise harmless very large coefficient.
            let basis = output(1, n, original, &mut b)?;
            return Ok((basis, Matrix::eye(1)));
        }
        let (mu, g) = gram_schmidt(&original, n, &mut b)?;
        let mut t = vec![vec![BigInt::zero(); m]; m];
        for (i, row) in t.iter_mut().enumerate() {
            row[i] = BigInt::one();
        }
        let mut basis = Basis {
            y: b.clone_rows(&original)?,
            inverse: b.clone_rows(&t)?,
            t,
            mu,
            g,
        };
        let mut k = 1;
        while k < m {
            b.claim(1)?;
            if !b.small(&basis.mu[k][k - 1])? {
                basis.reduce(k, k - 1, &mut b)?;
            }
            if lovasz(&basis.mu, &basis.g, &delta, k, &mut b)? {
                for j in (0..k - 1).rev() {
                    if !b.small(&basis.mu[k][j])? {
                        basis.reduce(k, j, &mut b)?;
                    }
                }
                k += 1;
            } else {
                basis.swap(k, &mut b)?;
                k = (k - 1).max(1);
            }
        }
        basis.certify(&original, n, &delta, &mut b)?;
        Ok((
            output(m, n, basis.y, &mut b)?,
            output(m, m, basis.t, &mut b)?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn matrix(values: &[&[i64]]) -> Matrix {
        Matrix::from_rows(
            values
                .iter()
                .map(|r| r.iter().map(|n| Expr::int(*n)).collect())
                .collect(),
        )
        .unwrap()
    }
    #[test]
    fn classic_basis_ties_and_exact_integer_inverse_certificate() {
        let delta = Expr::rational(3, 4);
        let a = matrix(&[
            &[1, 0, 0, 0, -20160],
            &[0, 1, 0, 0, 33768],
            &[0, 0, 1, 0, 39578],
            &[0, 0, 0, 1, 47757],
        ]);
        let expected = matrix(&[
            &[10, -3, -2, 8, -4],
            &[3, -9, 8, 1, -11],
            &[-3, 13, -9, -3, -9],
            &[-12, -7, -11, 9, -1],
        ]);
        let (basis, t) = a.lll_transform(&delta).unwrap();
        assert_eq!(basis, expected);
        assert_eq!(t.matmul(&a).unwrap(), basis);
        assert!(t.det().unwrap() == Expr::one() || t.det().unwrap() == Expr::minus_one());
        assert_eq!(a.lll(&delta).unwrap(), basis);
        let (basis, t) = matrix(&[&[2, 0], &[3, 2]]).lll_transform(&delta).unwrap();
        assert_eq!(basis, matrix(&[&[2, 0], &[-1, 2]]));
        assert_eq!(t, matrix(&[&[1, 0], &[-2, 1]]));
        let (basis, t) = matrix(&[&[2, 0], &[-3, 2]]).lll_transform(&delta).unwrap();
        assert_eq!(basis, matrix(&[&[2, 0], &[-1, 2]]));
        assert_eq!(t, matrix(&[&[1, 0], &[1, 1]]));
    }
    #[test]
    fn small_exact_bases_domains_shapes_and_empty() {
        let delta = Expr::rational(3, 4);
        for seed in 0..160i64 {
            let numbers: Vec<_> = (0..6i64)
                .map(|j| ((seed * (j + 3) + j * j + 11) % 13) - 6)
                .collect();
            let a = Matrix::new(2, 3, numbers.iter().map(|v| Expr::int(*v)).collect()).unwrap();
            let independent = (0..3).any(|i| {
                (i + 1..3).any(|j| numbers[i] * numbers[3 + j] != numbers[j] * numbers[3 + i])
            });
            let result = a.lll_transform(&delta);
            if independent {
                let (b, t) = result.unwrap();
                assert_eq!(t.matmul(&a).unwrap(), b);
            } else {
                assert_eq!(result, Err(LllError::Rank));
            }
        }
        assert_eq!(
            Matrix::zeros(0, 7).lll_transform(&delta).unwrap(),
            (Matrix::zeros(0, 7), Matrix::zeros(0, 0))
        );
        assert_eq!(Matrix::zeros(1, 3).lll(&delta), Err(LllError::Rank));
        assert_eq!(Matrix::zeros(3, 0).lll(&delta), Err(LllError::Shape));
        assert_eq!(
            Matrix::eye(1).lll(&Expr::rational(1, 4)),
            Err(LllError::Delta)
        );
        assert_eq!(Matrix::eye(1).lll(&Expr::one()), Err(LllError::Delta));
        assert_eq!(matrix(&[&[1, 2], &[2, 4]]).lll(&delta), Err(LllError::Rank));
        let input = Matrix::new(1, 2, vec![Expr::num(Number::float(2.0)), Expr::int(-3)]).unwrap();
        assert_eq!(input.lll(&delta).unwrap(), matrix(&[&[2, -3]]));
        for entry in [
            Expr::rational(1, 2),
            Expr::symbol("x"),
            Expr::i(),
            Expr::num(Number::float(1.5)),
        ] {
            assert_eq!(
                Matrix::new(1, 1, vec![entry]).unwrap().lll(&delta),
                Err(LllError::Domain)
            );
        }
    }
    #[test]
    fn unimodularity_is_verified_beyond_reconstruction_and_reducedness() {
        let original = vec![
            vec![BigInt::one(), BigInt::zero()],
            vec![BigInt::zero(), BigInt::one()],
        ];
        let nonunit = vec![
            vec![BigInt::from(2), BigInt::zero()],
            vec![BigInt::zero(), BigInt::from(2)],
        ];
        let basis = Basis {
            y: nonunit.clone(),
            t: nonunit,
            inverse: original.clone(),
            mu: vec![],
            g: vec![],
        };
        let delta = Q::new(BigInt::from(3), BigInt::from(4));
        // 2*I is reduced and T*I=2*I holds, but this is a proper sublattice.
        assert_eq!(
            basis.certify(&original, 2, &delta, &mut Budget::new()),
            Err(LllError::Certificate)
        );
    }

    #[test]
    fn delta_float_coercion_matches_rational_field_and_extremes_are_checked() {
        let mut b = Budget::new();
        for (value, n, d) in [(0.75, 3, 4), (0.3, 3, 10), (0.7, 7, 10)] {
            assert_eq!(
                b.delta(&Expr::num(Number::float(value))).unwrap(),
                Q::new(BigInt::from(n), BigInt::from(d))
            );
        }
        for exp in [1_000_000_000, -1_000_000_000, i64::MAX, i64::MIN] {
            let value = Expr::num(Number::MpFloat(
                BigFloat {
                    man: BigInt::one(),
                    exp: exp.into(),
                },
                80,
            ));
            assert_eq!(b.delta(&value), Err(LllError::Delta));
        }
        let a = matrix(&[&[2, 0], &[3, 2]]);
        assert_eq!(
            a.lll(&Expr::num(Number::float(0.3))).unwrap(),
            a.lll(&Expr::rational(3, 10)).unwrap()
        );
    }
    #[test]
    fn budgets_precede_allocations_products_gcd_and_latent_evaluation() {
        let delta = Expr::rational(3, 4);
        assert_eq!(
            preflight(usize::MAX, usize::MAX, &delta),
            Err(LllError::Resource)
        );
        assert_eq!(preflight(0, usize::MAX, &delta), Ok(()));
        let held =
            Expr::from_unevaluated_kind(crate::Kind::Pow(Expr::int(2), Expr::int(1_000_000_000)));
        assert_eq!(
            Matrix::new(1, 1, vec![held.clone()]).unwrap().lll(&delta),
            Err(LllError::Domain)
        );
        assert_eq!(Matrix::eye(1).lll(&held), Err(LllError::Parameter));
        let huge = Expr::num(Number::MpFloat(
            BigFloat {
                man: BigInt::one(),
                exp: 1_000_000_000,
            },
            80,
        ));
        assert_eq!(
            Matrix::new(1, 1, vec![huge]).unwrap().lll(&delta),
            Err(LllError::Resource)
        );
        let mut header_budget = Budget::new();
        header_budget.shape(1, 1).unwrap();
        assert!(
            header_budget.output_storage >= 2 * std::mem::size_of::<crate::expr::Node>() as u64 * 8
        );
        let mut b = Budget::new();
        b.work = 0;
        assert_eq!(
            b.mul(&BigInt::from(2), &BigInt::from(3)),
            Err(LllError::Resource)
        );
        assert_eq!(
            b.rational(BigInt::from(7), BigInt::from(11)),
            Err(LllError::Resource)
        );
        assert_eq!(b.clone_q(&Q::one()), Err(LllError::Resource));
        assert_eq!(
            b.clone_rows(&[vec![BigInt::one()]]),
            Err(LllError::Resource)
        );
        assert_eq!(
            b.float_integer(&BigFloat {
                man: BigInt::one() << 128,
                exp: -128
            }),
            Err(LllError::Resource)
        );
        let mut b = Budget::new();
        b.height = MAX_BITS;
        assert_eq!(
            b.mul(&(BigInt::one() << 600_000), &(BigInt::one() << 600_000)),
            Err(LllError::Resource)
        );
    }
}
