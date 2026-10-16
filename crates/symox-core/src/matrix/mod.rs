//! Dense symbolic matrices with SymPy compatible semantics.
//!
//! [`Matrix`] stores its entries row-major as canonical [`Expr`]s. Linear
//! algebra follows SymPy's algorithms (and therefore its output forms):
//! exact rational matrices are handled by fraction-free integer
//! elimination, polynomial matrices by fraction-free elimination over an
//! internal multivariate polynomial ring, everything else by the generic
//! expression based algorithms with SymPy's `dotprodsimp` intermediate
//! simplification.
//!
//! ```
//! use symox::matrix::Matrix;
//! use symox::prelude::*;
//!
//! let m = Matrix::from_rows(vec![
//!     vec![Expr::int(1), Expr::int(2)],
//!     vec![Expr::int(3), Expr::int(4)],
//! ]).unwrap();
//! assert_eq!(m.det().unwrap().to_string(), "-2");
//! assert_eq!(m.inv().unwrap().to_string(), "Matrix([[-2, 1], [3/2, -1/2]])");
//! ```

pub mod domain;
pub mod eigen;
pub mod ff;
pub mod int;
pub mod linalg;
pub mod mpoly;
pub mod ratfunc;

use std::fmt;

use crate::diff::diff;
use crate::evalf::evalf;
use crate::expand::expand;
use crate::expr::{Expr, Func};
use crate::functions::func;
use crate::ops::{add, mul2, neg, pow};
use crate::printer::to_str;
use crate::subs::subs;

pub use linalg::DetMethod;

/// Errors raised by matrix operations (mirroring SymPy's exceptions).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MatrixError {
    /// `ShapeError`
    Shape(String),
    /// `NonSquareMatrixError`
    NonSquare(String),
    /// `NonInvertibleMatrixError`
    NonInvertible(String),
    /// `NotImplementedError`
    NotImplemented(String),
    /// `ValueError`
    Value(String),
    /// `IndexError`
    Index(String),
    /// `TypeError`
    Type(String),
    /// `MatrixError`
    Other(String),
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            MatrixError::Shape(s)
            | MatrixError::NonSquare(s)
            | MatrixError::NonInvertible(s)
            | MatrixError::NotImplemented(s)
            | MatrixError::Value(s)
            | MatrixError::Index(s)
            | MatrixError::Type(s)
            | MatrixError::Other(s) => s,
        };
        f.write_str(s)
    }
}

impl std::error::Error for MatrixError {}

pub type MResult<T> = Result<T, MatrixError>;

pub(crate) fn non_square() -> MatrixError {
    MatrixError::NonSquare(String::new())
}

pub(crate) fn non_invertible() -> MatrixError {
    MatrixError::NonInvertible("Matrix det == 0; not invertible.".into())
}

/// A dense matrix of expressions (row-major).
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<Expr>,
}

impl Matrix {
    // ------------------------------------------------------------------
    // construction
    // ------------------------------------------------------------------

    /// Matrix from row-major data (`data.len() == rows * cols`).
    pub fn new(rows: usize, cols: usize, data: Vec<Expr>) -> MResult<Matrix> {
        if data.len() != rows * cols {
            return Err(MatrixError::Value(format!(
                "List length should be equal to rows*columns: {} != {}",
                data.len(),
                rows * cols
            )));
        }
        Ok(Matrix { rows, cols, data })
    }

    pub(crate) fn new_unchecked(rows: usize, cols: usize, data: Vec<Expr>) -> Matrix {
        debug_assert_eq!(data.len(), rows * cols);
        Matrix { rows, cols, data }
    }

    /// Matrix from a list of rows (all of the same length).
    pub fn from_rows(rows: Vec<Vec<Expr>>) -> MResult<Matrix> {
        let r = rows.len();
        let c = rows.first().map(|v| v.len()).unwrap_or(0);
        let mut data = Vec::with_capacity(r * c);
        for row in rows {
            if row.len() != c {
                return Err(MatrixError::Value("mismat row length".into()));
            }
            data.extend(row);
        }
        Ok(Matrix { rows: r, cols: c, data })
    }

    /// Column vector.
    pub fn column(v: Vec<Expr>) -> Matrix {
        let n = v.len();
        Matrix { rows: n, cols: 1, data: v }
    }

    /// Row vector.
    pub fn row_vector(v: Vec<Expr>) -> Matrix {
        let n = v.len();
        Matrix { rows: 1, cols: n, data: v }
    }

    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> Expr) -> Matrix {
        let mut data = Vec::with_capacity(rows * cols);
        for i in 0..rows {
            for j in 0..cols {
                data.push(f(i, j));
            }
        }
        Matrix { rows, cols, data }
    }

    pub fn zeros(rows: usize, cols: usize) -> Matrix {
        Matrix { rows, cols, data: vec![Expr::zero(); rows * cols] }
    }

    pub fn ones(rows: usize, cols: usize) -> Matrix {
        Matrix { rows, cols, data: vec![Expr::one(); rows * cols] }
    }

    pub fn eye(n: usize) -> Matrix {
        Matrix::from_fn(n, n, |i, j| if i == j { Expr::one() } else { Expr::zero() })
    }

    /// Diagonal matrix with the given diagonal entries.
    pub fn diag(entries: &[Expr]) -> Matrix {
        let n = entries.len();
        Matrix::from_fn(n, n, |i, j| if i == j { entries[i].clone() } else { Expr::zero() })
    }

    /// Block diagonal matrix.
    pub fn diag_blocks(blocks: &[Matrix]) -> Matrix {
        let r: usize = blocks.iter().map(|b| b.rows).sum();
        let c: usize = blocks.iter().map(|b| b.cols).sum();
        let mut m = Matrix::zeros(r, c);
        let (mut r0, mut c0) = (0, 0);
        for b in blocks {
            for i in 0..b.rows {
                for j in 0..b.cols {
                    m.data[(r0 + i) * c + c0 + j] = b.get(i, j).clone();
                }
            }
            r0 += b.rows;
            c0 += b.cols;
        }
        m
    }

    // ------------------------------------------------------------------
    // access
    // ------------------------------------------------------------------

    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }

    pub fn is_vector(&self) -> bool {
        self.rows == 1 || self.cols == 1
    }

    /// Row-major entries.
    pub fn data(&self) -> &[Expr] {
        &self.data
    }

    pub fn into_data(self) -> Vec<Expr> {
        self.data
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &Expr {
        &self.data[i * self.cols + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: Expr) {
        self.data[i * self.cols + j] = v;
    }

    pub fn row(&self, i: usize) -> Matrix {
        Matrix::row_vector(self.data[i * self.cols..(i + 1) * self.cols].to_vec())
    }

    pub fn col(&self, j: usize) -> Matrix {
        Matrix::column((0..self.rows).map(|i| self.get(i, j).clone()).collect())
    }

    pub fn to_rows(&self) -> Vec<Vec<Expr>> {
        (0..self.rows)
            .map(|i| self.data[i * self.cols..(i + 1) * self.cols].to_vec())
            .collect()
    }

    /// Sub-matrix with the given rows and columns (in that order).
    pub fn extract(&self, rows: &[usize], cols: &[usize]) -> MResult<Matrix> {
        if rows.iter().any(|&r| r >= self.rows) || cols.iter().any(|&c| c >= self.cols) {
            return Err(MatrixError::Index("index out of range".into()));
        }
        let mut data = Vec::with_capacity(rows.len() * cols.len());
        for &r in rows {
            for &c in cols {
                data.push(self.get(r, c).clone());
            }
        }
        Ok(Matrix { rows: rows.len(), cols: cols.len(), data })
    }

    /// `self[r0:r1, c0:c1]`.
    pub fn slice(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Matrix {
        let rows: Vec<usize> = (r0..r1.min(self.rows)).collect();
        let cols: Vec<usize> = (c0..c1.min(self.cols)).collect();
        self.extract(&rows, &cols).unwrap()
    }

    pub fn reshape(&self, rows: usize, cols: usize) -> MResult<Matrix> {
        if rows * cols != self.data.len() {
            return Err(MatrixError::Value("Invalid reshape parameters".into()));
        }
        Ok(Matrix { rows, cols, data: self.data.clone() })
    }

    // ------------------------------------------------------------------
    // structural operations
    // ------------------------------------------------------------------

    pub fn transpose(&self) -> Matrix {
        Matrix::from_fn(self.cols, self.rows, |i, j| self.get(j, i).clone())
    }

    /// Entry-wise complex conjugate.
    pub fn conjugate(&self) -> Matrix {
        self.applyfunc(|e| func(Func::Conjugate, vec![e.clone()]))
    }

    /// Conjugate transpose (`H`).
    pub fn adjoint(&self) -> Matrix {
        self.transpose().conjugate()
    }

    pub fn trace(&self) -> MResult<Expr> {
        if !self.is_square() {
            return Err(non_square());
        }
        Ok(add(&(0..self.rows).map(|i| self.get(i, i).clone()).collect::<Vec<_>>()))
    }

    pub fn diagonal(&self) -> Matrix {
        let n = self.rows.min(self.cols);
        Matrix::row_vector((0..n).map(|i| self.get(i, i).clone()).collect())
    }

    /// Horizontal concatenation.
    pub fn hstack(ms: &[&Matrix]) -> MResult<Matrix> {
        if ms.is_empty() {
            return Ok(Matrix::zeros(0, 0));
        }
        let r = ms[0].rows;
        if ms.iter().any(|m| m.rows != r) {
            return Err(MatrixError::Shape(
                "`self` and `other` must have the same number of rows.".into(),
            ));
        }
        let c: usize = ms.iter().map(|m| m.cols).sum();
        let mut data = Vec::with_capacity(r * c);
        for i in 0..r {
            for m in ms {
                data.extend_from_slice(&m.data[i * m.cols..(i + 1) * m.cols]);
            }
        }
        Ok(Matrix { rows: r, cols: c, data })
    }

    /// Vertical concatenation.
    pub fn vstack(ms: &[&Matrix]) -> MResult<Matrix> {
        if ms.is_empty() {
            return Ok(Matrix::zeros(0, 0));
        }
        let c = ms[0].cols;
        if ms.iter().any(|m| m.cols != c) {
            return Err(MatrixError::Shape(
                "`self` and `other` must have the same number of columns.".into(),
            ));
        }
        let r: usize = ms.iter().map(|m| m.rows).sum();
        let mut data = Vec::with_capacity(r * c);
        for m in ms {
            data.extend_from_slice(&m.data);
        }
        Ok(Matrix { rows: r, cols: c, data })
    }

    pub fn row_join(&self, o: &Matrix) -> MResult<Matrix> {
        if self.cols == 0 && self.rows != o.rows {
            return Ok(o.clone());
        }
        Matrix::hstack(&[self, o])
    }

    pub fn col_join(&self, o: &Matrix) -> MResult<Matrix> {
        if self.rows == 0 && self.cols != o.cols {
            return Ok(o.clone());
        }
        Matrix::vstack(&[self, o])
    }

    /// Insert the rows of `o` before row `pos` (SymPy's `row_insert`).
    pub fn row_insert(&self, pos: isize, o: &Matrix) -> MResult<Matrix> {
        if self.data.is_empty() {
            return Ok(o.clone());
        }
        let pos = clamp_insert(pos, self.rows);
        if self.cols != o.cols {
            return Err(MatrixError::Shape(
                "The matrices have incompatible number of columns".into(),
            ));
        }
        let mut data = Vec::with_capacity(self.data.len() + o.data.len());
        data.extend_from_slice(&self.data[..pos * self.cols]);
        data.extend_from_slice(&o.data);
        data.extend_from_slice(&self.data[pos * self.cols..]);
        Ok(Matrix { rows: self.rows + o.rows, cols: self.cols, data })
    }

    /// Insert the columns of `o` before column `pos`.
    pub fn col_insert(&self, pos: isize, o: &Matrix) -> MResult<Matrix> {
        if self.data.is_empty() {
            return Ok(o.clone());
        }
        let pos = clamp_insert(pos, self.cols);
        if self.rows != o.rows {
            return Err(MatrixError::Shape(
                "The matrices have incompatible number of rows".into(),
            ));
        }
        let c = self.cols + o.cols;
        let mut data = Vec::with_capacity(self.rows * c);
        for i in 0..self.rows {
            let r = &self.data[i * self.cols..(i + 1) * self.cols];
            data.extend_from_slice(&r[..pos]);
            data.extend_from_slice(&o.data[i * o.cols..(i + 1) * o.cols]);
            data.extend_from_slice(&r[pos..]);
        }
        Ok(Matrix { rows: self.rows, cols: c, data })
    }

    pub fn row_del(&mut self, i: usize) -> MResult<()> {
        if i >= self.rows {
            return Err(MatrixError::Index("Index out of range".into()));
        }
        self.data.drain(i * self.cols..(i + 1) * self.cols);
        self.rows -= 1;
        Ok(())
    }

    pub fn col_del(&mut self, j: usize) -> MResult<()> {
        if j >= self.cols {
            return Err(MatrixError::Index("Index out of range".into()));
        }
        let mut data = Vec::with_capacity(self.rows * (self.cols - 1));
        for i in 0..self.rows {
            for k in 0..self.cols {
                if k != j {
                    data.push(self.get(i, k).clone());
                }
            }
        }
        self.data = data;
        self.cols -= 1;
        Ok(())
    }

    pub fn row_swap(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }
        for k in 0..self.cols {
            self.data.swap(i * self.cols + k, j * self.cols + k);
        }
    }

    pub fn col_swap(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }
        for k in 0..self.rows {
            self.data.swap(k * self.cols + i, k * self.cols + j);
        }
    }

    /// Apply a list of row swaps (SymPy's `permute_rows` with a swap list).
    pub fn permute_rows_swaps(&self, swaps: &[(usize, usize)]) -> Matrix {
        let mut m = self.clone();
        for &(a, b) in swaps {
            m.row_swap(a, b);
        }
        m
    }

    // ------------------------------------------------------------------
    // predicates
    // ------------------------------------------------------------------

    pub fn is_zero_matrix(&self) -> bool {
        self.data.iter().all(|e| e.is_zero())
    }

    pub fn is_symmetric(&self) -> bool {
        self.is_square()
            && (0..self.rows).all(|i| (0..i).all(|j| self.get(i, j) == self.get(j, i)))
    }

    pub fn is_diagonal(&self) -> bool {
        (0..self.rows).all(|i| (0..self.cols).all(|j| i == j || self.get(i, j).is_zero()))
    }

    pub fn is_upper(&self) -> bool {
        (0..self.rows).all(|i| (0..i.min(self.cols)).all(|j| self.get(i, j).is_zero()))
    }

    pub fn is_lower(&self) -> bool {
        (0..self.rows).all(|i| (i + 1..self.cols).all(|j| self.get(i, j).is_zero()))
    }

    pub fn free_symbols(&self) -> Vec<Expr> {
        let mut out: Vec<Expr> = Vec::new();
        for e in &self.data {
            for s in crate::traverse::free_symbols(e) {
                if !out.contains(&s) {
                    out.push(s);
                }
            }
        }
        out.sort_by(|a, b| a.cmp_struct(b));
        out
    }

    // ------------------------------------------------------------------
    // arithmetic
    // ------------------------------------------------------------------

    pub fn add(&self, o: &Matrix) -> MResult<Matrix> {
        if self.shape() != o.shape() {
            return Err(MatrixError::Shape(format!(
                "Matrix size mismatch: {:?} + {:?}.",
                self.shape(),
                o.shape()
            )));
        }
        Ok(Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().zip(&o.data).map(|(a, b)| crate::ops::add2(a, b)).collect(),
        })
    }

    pub fn sub(&self, o: &Matrix) -> MResult<Matrix> {
        if self.shape() != o.shape() {
            return Err(MatrixError::Shape(format!(
                "Matrix size mismatch: {:?} + {:?}.",
                self.shape(),
                o.shape()
            )));
        }
        Ok(Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().zip(&o.data).map(|(a, b)| crate::ops::sub(a, b)).collect(),
        })
    }

    pub fn neg(&self) -> Matrix {
        self.applyfunc(neg)
    }

    /// `self * k` (entry * scalar).
    pub fn scale(&self, k: &Expr) -> Matrix {
        self.applyfunc(|e| mul2(e, k))
    }

    /// `k * self` (scalar * entry).
    pub fn scale_left(&self, k: &Expr) -> Matrix {
        self.applyfunc(|e| mul2(k, e))
    }

    /// Matrix product (entries are sums of products, not expanded).
    pub fn matmul(&self, o: &Matrix) -> MResult<Matrix> {
        if self.cols != o.rows {
            return Err(MatrixError::Shape(format!(
                "Matrix size mismatch: {:?} * {:?}.",
                self.shape(),
                o.shape()
            )));
        }
        if let (Some(a), Some(b)) = (linalg::RatMat::from_matrix(self), linalg::RatMat::from_matrix(o)) {
            return Ok(a.matmul(&b).to_matrix());
        }
        let n = self.cols;
        let mut data = Vec::with_capacity(self.rows * o.cols);
        let mut terms: Vec<Expr> = Vec::with_capacity(n);
        for i in 0..self.rows {
            for j in 0..o.cols {
                terms.clear();
                for k in 0..n {
                    let a = self.get(i, k);
                    let b = o.get(k, j);
                    if a.is_zero() || b.is_zero() {
                        continue;
                    }
                    terms.push(mul2(a, b));
                }
                data.push(add(&terms));
            }
        }
        Ok(Matrix { rows: self.rows, cols: o.cols, data })
    }

    /// Entry-wise (Hadamard) product.
    pub fn multiply_elementwise(&self, o: &Matrix) -> MResult<Matrix> {
        if self.shape() != o.shape() {
            return Err(MatrixError::Shape("Matrix size mismatch".into()));
        }
        Ok(Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().zip(&o.data).map(|(a, b)| mul2(a, b)).collect(),
        })
    }

    /// Integer power (negative powers go through the inverse).
    pub fn powi(&self, n: i64) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare("Power of a non-square matrix".into()));
        }
        if n == 0 {
            return Ok(Matrix::eye(self.rows));
        }
        if n == 1 {
            return Ok(self.clone());
        }
        if self.is_diagonal() {
            let e = Expr::int(n);
            return Ok(Matrix::from_fn(self.rows, self.cols, |i, j| {
                if i == j {
                    pow(self.get(i, i), &e)
                } else {
                    Expr::zero()
                }
            }));
        }
        if self.rows == 1 {
            return Ok(Matrix::new_unchecked(1, 1, vec![pow(&self.data[0], &Expr::int(n))]));
        }
        let (base, k) = if n < 0 { (self.inv()?, n.unsigned_abs()) } else { (self.clone(), n as u64) };
        if let Some(r) = linalg::RatMat::from_matrix(&base) {
            return Ok(r.pow(k).to_matrix());
        }
        let mut prevsimp = vec![true; base.data.len()];
        Ok(pow_rec(&base, k, &mut prevsimp))
    }

    // ------------------------------------------------------------------
    // element-wise operations
    // ------------------------------------------------------------------

    pub fn applyfunc(&self, f: impl Fn(&Expr) -> Expr) -> Matrix {
        Matrix { rows: self.rows, cols: self.cols, data: self.data.iter().map(f).collect() }
    }

    pub fn try_applyfunc<E>(&self, mut f: impl FnMut(&Expr) -> Result<Expr, E>) -> Result<Matrix, E> {
        let mut data = Vec::with_capacity(self.data.len());
        for e in &self.data {
            data.push(f(e)?);
        }
        Ok(Matrix { rows: self.rows, cols: self.cols, data })
    }

    pub fn subs(&self, pairs: &[(Expr, Expr)]) -> Matrix {
        self.applyfunc(|e| subs(e, pairs))
    }

    pub fn xreplace(&self, map: &rustc_hash::FxHashMap<Expr, Expr>) -> Matrix {
        self.applyfunc(|e| crate::subs::xreplace(e, map))
    }

    pub fn diff(&self, x: &Expr) -> Matrix {
        self.applyfunc(|e| diff(e, x))
    }

    pub fn diff_many(&self, vars: &[(Expr, u32)]) -> Matrix {
        self.applyfunc(|e| crate::diff::diff_many(e, vars))
    }

    pub fn expand(&self) -> Matrix {
        self.applyfunc(expand)
    }

    pub fn evalf(&self) -> Matrix {
        self.applyfunc(evalf)
    }

    pub fn simplify(&self) -> Matrix {
        self.applyfunc(crate::simplify::simplify)
    }

    /// Cancel common factors of every entry (rational functions).
    pub fn cancel(&self) -> Matrix {
        self.applyfunc(ratfunc::cancel)
    }

    // ------------------------------------------------------------------
    // vectors
    // ------------------------------------------------------------------

    /// Dot product of two vectors of the same length (no conjugation).
    pub fn dot(&self, o: &Matrix) -> MResult<Expr> {
        if !self.is_vector() || !o.is_vector() {
            // SymPy: matrix-matrix "dot" is a matrix product flattened
            return Err(MatrixError::Shape(
                "Dimensions incorrect for dot product".into(),
            ));
        }
        if self.len() != o.len() {
            return Err(MatrixError::Shape(format!(
                "Dimensions incorrect for dot product: {:?}, {:?}",
                self.shape(),
                o.shape()
            )));
        }
        let terms: Vec<Expr> = self.data.iter().zip(&o.data).map(|(a, b)| mul2(a, b)).collect();
        Ok(add(&terms))
    }

    /// Cross product of 3-vectors (shape of `self`).
    pub fn cross(&self, o: &Matrix) -> MResult<Matrix> {
        if self.len() != 3 || o.len() != 3 {
            return Err(MatrixError::Shape(format!(
                "Dimensions incorrect for cross product: {:?} x {:?}",
                self.shape(),
                o.shape()
            )));
        }
        let (a, b) = (&self.data, &o.data);
        let c = |i: usize, j: usize| crate::ops::sub(&mul2(&a[i], &b[j]), &mul2(&a[j], &b[i]));
        Ok(Matrix { rows: self.rows, cols: self.cols, data: vec![c(1, 2), c(2, 0), c(0, 1)] })
    }

    /// Frobenius norm (2-norm for vectors), as SymPy's `norm()`.
    pub fn norm(&self) -> Expr {
        let terms: Vec<Expr> = self
            .data
            .iter()
            .map(|e| pow(&func(Func::Abs, vec![e.clone()]), &Expr::int(2)))
            .collect();
        crate::ops::sqrt(&add(&terms))
    }

    // ------------------------------------------------------------------
    // calculus
    // ------------------------------------------------------------------

    /// Jacobian of the vector `self` with respect to the variables `vars`.
    pub fn jacobian(&self, vars: &[Expr]) -> MResult<Matrix> {
        if !self.is_vector() {
            return Err(MatrixError::Type(
                "``self`` must be a row or a column matrix".into(),
            ));
        }
        let fs = &self.data;
        Ok(Matrix::from_fn(fs.len(), vars.len(), |i, j| diff(&fs[i], &vars[j])))
    }

    // ------------------------------------------------------------------
    // printing
    // ------------------------------------------------------------------

    /// SymPy's `str(M)`: `Matrix([[1, 2], [3, 4]])`.
    pub fn to_str(&self) -> String {
        if self.rows == 0 || self.cols == 0 {
            return format!("Matrix({}, {}, [])", self.rows, self.cols);
        }
        let mut s = String::from("Matrix([");
        for i in 0..self.rows {
            if i > 0 {
                s.push_str(", ");
            }
            s.push('[');
            for j in 0..self.cols {
                if j > 0 {
                    s.push_str(", ");
                }
                s.push_str(&to_str(self.get(i, j)));
            }
            s.push(']');
        }
        s.push_str("])");
        s
    }

    /// SymPy's `repr(M)` (multi-line, right aligned columns).
    pub fn to_repr(&self) -> String {
        if self.rows == 0 || self.cols == 0 {
            return format!("Matrix({}, {}, [])", self.rows, self.cols);
        }
        let strs: Vec<String> = self.data.iter().map(to_str).collect();
        let mut width = vec![0usize; self.cols];
        for i in 0..self.rows {
            for j in 0..self.cols {
                width[j] = width[j].max(strs[i * self.cols + j].chars().count());
            }
        }
        let mut lines = Vec::with_capacity(self.rows);
        for i in 0..self.rows {
            let mut l = String::from("[");
            for j in 0..self.cols {
                if j > 0 {
                    l.push_str(", ");
                }
                let s = &strs[i * self.cols + j];
                for _ in s.chars().count()..width[j] {
                    l.push(' ');
                }
                l.push_str(s);
            }
            l.push(']');
            lines.push(l);
        }
        if self.rows == 1 {
            format!("Matrix([{}])", lines[0])
        } else {
            format!("Matrix([\n{}])", lines.join(",\n"))
        }
    }
}

fn clamp_insert(pos: isize, n: usize) -> usize {
    let n = n as isize;
    let p = if pos < 0 { (pos + n).max(0) } else { pos.min(n) };
    p as usize
}

/// SymPy's `_eval_pow_by_recursion_dotprodsimp`.
fn pow_rec(m: &Matrix, n: u64, prevsimp: &mut [bool]) -> Matrix {
    if n == 1 {
        return m.clone();
    }
    let (a, b) = if n % 2 == 1 {
        (m.clone(), pow_rec(m, n - 1, prevsimp))
    } else {
        let h = pow_rec(m, n / 2, prevsimp);
        (h.clone(), h)
    };
    let p = a.matmul(&b).unwrap();
    let data = p
        .data
        .iter()
        .enumerate()
        .map(|(i, e)| {
            if prevsimp[i] {
                let (s, flag) = ratfunc::dotprodsimp_withsimp(e);
                prevsimp[i] = flag;
                s
            } else {
                e.clone()
            }
        })
        .collect();
    Matrix { rows: p.rows, cols: p.cols, data }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_str())
    }
}

// ----------------------------------------------------------------------
// free functions
// ----------------------------------------------------------------------

/// Jacobian of `exprs` with respect to `vars`.
pub fn jacobian(exprs: &[Expr], vars: &[Expr]) -> Matrix {
    Matrix::from_fn(exprs.len(), vars.len(), |i, j| diff(&exprs[i], &vars[j]))
}

/// Hessian matrix of `f` with respect to `vars`.
pub fn hessian(f: &Expr, vars: &[Expr]) -> Matrix {
    let n = vars.len();
    let first: Vec<Expr> = vars.iter().map(|v| diff(f, v)).collect();
    let mut m = Matrix::zeros(n, n);
    for i in 0..n {
        for j in i..n {
            let d = diff(&first[i], &vars[j]);
            m.set(i, j, d.clone());
            m.set(j, i, d);
        }
    }
    m
}

/// Wronskian of `fns` with respect to `x` (determinant of derivatives).
pub fn wronskian(fns: &[Expr], x: &Expr) -> MResult<Expr> {
    let n = fns.len();
    if n == 0 {
        return Ok(Expr::one());
    }
    let mut rows: Vec<Vec<Expr>> = Vec::with_capacity(n);
    // W[i, j] = d^j f_i / dx^j
    for f in fns {
        let mut row = Vec::with_capacity(n);
        let mut d = f.clone();
        for j in 0..n {
            if j > 0 {
                d = diff(&d, x);
            }
            row.push(d.clone());
        }
        rows.push(row);
    }
    Matrix::from_rows(rows)?.det()
}

#[cfg(test)]
mod tests;
