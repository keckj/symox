//! Further decompositions: Cholesky, LDL, QR (Gram-Schmidt),
//! diagonalisation, and norms (ports of SymPy's dense algorithms).

use super::ratfunc::is_zero_assume;
use super::{non_square, MResult, Matrix, MatrixError};
use crate::expand::expand;
use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::ops::{add, div, mul2, pow, sqrt, sub};

/// Conjugate, evaluated for numeric expressions (`conjugate(sqrt(2))` is
/// `sqrt(2)`, `I` becomes `-I`).
fn conj(e: &Expr) -> Expr {
    if e.mask() == 0 && crate::functions::is_numeric_value(e) {
        if let Some(z) = crate::evalf::eval_c64(e) {
            let mut map = rustc_hash::FxHashMap::default();
            map.insert(Expr::i(), crate::ops::neg(&Expr::i()));
            let c = crate::subs::xreplace(e, &map);
            if let Some(w) = crate::evalf::eval_c64(&c) {
                if (w - z.conj()).norm() <= 1e-12 * (1.0 + z.norm()) {
                    return c;
                }
            }
        }
    }
    func(Func::Conjugate, vec![e.clone()])
}

fn abs(e: &Expr) -> Expr {
    func(Func::Abs, vec![e.clone()])
}

/// `e.is_positive is False` for the cases we can decide.
fn is_not_positive(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => !n.is_positive(),
        _ => {
            if e.mask() == 0 && crate::functions::is_numeric_value(e) {
                if let Some(z) = crate::evalf::eval_c64(e) {
                    return z.im.abs() < 1e-12 && z.re <= 0.0;
                }
            }
            false
        }
    }
}

impl Matrix {
    /// `M == M.H` (decidable cases only; `None` when undecided).
    pub fn is_hermitian(&self) -> Option<bool> {
        if !self.is_square() {
            return Some(false);
        }
        let n = self.rows();
        let mut undecided = false;
        for i in 0..n {
            for j in i..n {
                let d = expand(&sub(self.get(i, j), &conj(self.get(j, i))));
                match is_zero_assume(&d) {
                    Some(true) => {}
                    Some(false) => return Some(false),
                    None => undecided = true,
                }
            }
        }
        if undecided {
            None
        } else {
            Some(true)
        }
    }

    fn check_hermitian(&self, hermitian: bool) -> MResult<()> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare("Matrix must be square.".into()));
        }
        if hermitian {
            if self.is_hermitian() != Some(true) {
                return Err(MatrixError::Value("Matrix must be Hermitian.".into()));
            }
        } else if !self.is_symmetric() {
            return Err(MatrixError::Value("Matrix must be symmetric.".into()));
        }
        Ok(())
    }

    /// Cholesky decomposition `M = L * L.H`.
    pub fn cholesky(&self, hermitian: bool) -> MResult<Matrix> {
        self.check_hermitian(hermitian)?;
        let n = self.rows();
        let mut l = Matrix::zeros(n, n);
        let c = |e: &Expr| if hermitian { conj(e) } else { e.clone() };
        for i in 0..n {
            for j in 0..i {
                let s: Vec<Expr> = (0..j).map(|k| mul2(l.get(i, k), &c(l.get(j, k)))).collect();
                let v = mul2(&pow(l.get(j, j), &Expr::minus_one()), &sub(self.get(i, j), &add(&s)));
                l.set(i, j, v);
            }
            let s: Vec<Expr> = (0..i).map(|k| mul2(l.get(i, k), &c(l.get(i, k)))).collect();
            let lii2 = sub(self.get(i, i), &add(&s));
            if hermitian && is_not_positive(&lii2) {
                return Err(MatrixError::NonPositiveDefinite("Matrix must be positive-definite".into()));
            }
            l.set(i, i, sqrt(&lii2));
        }
        Ok(l)
    }

    /// LDL decomposition `M = L * D * L.H`.
    pub fn ldl_decomposition(&self, hermitian: bool) -> MResult<(Matrix, Matrix)> {
        self.check_hermitian(hermitian)?;
        let n = self.rows();
        let mut d = Matrix::zeros(n, n);
        let mut l = Matrix::eye(n);
        let c = |e: &Expr| if hermitian { conj(e) } else { e.clone() };
        for i in 0..n {
            for j in 0..i {
                let s: Vec<Expr> = (0..j)
                    .map(|k| crate::ops::mul(&[l.get(i, k).clone(), c(l.get(j, k)), d.get(k, k).clone()]))
                    .collect();
                let v = mul2(&pow(d.get(j, j), &Expr::minus_one()), &sub(self.get(i, j), &add(&s)));
                l.set(i, j, v);
            }
            let s: Vec<Expr> = (0..i)
                .map(|k| crate::ops::mul(&[l.get(i, k).clone(), c(l.get(i, k)), d.get(k, k).clone()]))
                .collect();
            let dii = sub(self.get(i, i), &add(&s));
            if hermitian && is_not_positive(&dii) {
                return Err(MatrixError::NonPositiveDefinite("Matrix must be positive-definite".into()));
            }
            d.set(i, i, dii);
        }
        Ok((l, d))
    }

    /// QR decomposition by Gram-Schmidt (SymPy's `QRdecomposition`).
    pub fn qr_decomposition(&self) -> MResult<(Matrix, Matrix)> {
        let (rows, cols) = self.shape();
        let mut q = self.clone();
        let mut r = Matrix::zeros(cols, cols);
        let mut ranked = Vec::new();
        let col = |m: &Matrix, j: usize| m.col(j);
        let hdot = |u: &Matrix, v: &Matrix| -> Expr {
            let terms: Vec<Expr> = u.data().iter().zip(v.data()).map(|(a, b)| mul2(&conj(a), b)).collect();
            add(&terms)
        };
        for j in 0..cols {
            for i in 0..j {
                let qi = col(&q, i);
                if qi.is_zero_matrix() {
                    continue;
                }
                let qj = col(&q, j);
                let rij = expand(&div(&hdot(&qi, &qj), &hdot(&qi, &qi)));
                r.set(i, j, rij.clone());
                for k in 0..rows {
                    let v = sub(q.get(k, j), &mul2(q.get(k, i), &rij));
                    q.set(k, j, v);
                }
            }
            for k in 0..rows {
                let v = expand(q.get(k, j));
                q.set(k, j, v);
            }
            if !col(&q, j).data().iter().all(|e| is_zero_assume(e) == Some(true)) {
                ranked.push(j);
                r.set(j, j, Expr::one());
            }
        }
        let all_rows: Vec<usize> = (0..rows).collect();
        let all_cols: Vec<usize> = (0..cols).collect();
        let mut q = q.extract(&all_rows, &ranked)?;
        let mut r = r.extract(&ranked, &all_cols)?;
        for i in 0..q.cols() {
            let norm = q.col(i).norm();
            let inv = pow(&norm, &Expr::minus_one());
            for k in 0..q.rows() {
                let v = mul2(q.get(k, i), &inv);
                q.set(k, i, v);
            }
            for k in 0..r.cols() {
                let v = mul2(r.get(i, k), &norm);
                r.set(i, k, v);
            }
        }
        Ok((q, r))
    }

    /// `(P, D)` with `M = P * D * P**-1`.
    pub fn diagonalize(&self) -> MResult<(Matrix, Matrix)> {
        self.diagonalize_with_options(false, false, false)
    }

    pub fn diagonalize_with_options(&self, reals_only: bool, sort: bool, normalize: bool) -> MResult<(Matrix, Matrix)> {
        if !self.is_square() {
            return Err(non_square());
        }
        let mut ev = self.eigenvects()?;
        if sort {
            ev.sort_by(|a, b| crate::sortkey::cmp_default(&a.0, &b.0));
        }
        let mut cols: Vec<Matrix> = Vec::new();
        let mut diag: Vec<Expr> = Vec::new();
        for (val, mult, basis) in &ev {
            if basis.len() != *mult || (reals_only && crate::assumptions::is_real(val) != Some(true)) {
                return Err(MatrixError::Other("Matrix is not diagonalizable".into()));
            }
            for v in basis {
                diag.push(val.clone());
                cols.push(if normalize { v.scale(&pow(&v.norm(), &Expr::minus_one())) } else { v.clone() });
            }
        }
        let refs: Vec<&Matrix> = cols.iter().collect();
        Ok((Matrix::hstack(&refs)?, Matrix::diag(&diag)))
    }

    /// Unsupported eigenproblems must propagate their error rather than
    /// incorrectly claiming the matrix is not diagonalizable.
    pub fn is_diagonalizable(&self, reals_only: bool) -> MResult<bool> {
        if !self.is_square() {
            return Ok(false);
        }
        let ev = self.eigenvects()?;
        Ok(ev.iter().all(|(v, m, b)| b.len() == *m
            && (!reals_only || crate::assumptions::is_real(v) == Some(true))))
    }

    /// Column-stacked vector.
    pub fn vec(&self) -> Matrix {
        self.transpose().reshape(self.len(), 1).unwrap()
    }

    /// Non-zero entries (row-major).
    pub fn values(&self) -> Vec<Expr> {
        self.data().iter().filter(|e| !e.is_zero()).cloned().collect()
    }

    /// SymPy's `norm(ord)`: `None`/2 (vectors), 1, `oo`, `-oo` and the
    /// Frobenius norm of matrices.
    pub fn norm_ord(&self, ord: Option<&str>) -> MResult<Expr> {
        let vals = {
            let v = self.values();
            if v.is_empty() {
                vec![Expr::zero()]
            } else {
                v
            }
        };
        let is_vec = self.rows() == 1 || self.cols() == 1;
        let two = Expr::int(2);
        let maxf = |v: Vec<Expr>| func(Func::Max, v);
        let minf = |v: Vec<Expr>| func(Func::Min, v);
        if is_vec {
            return Ok(match ord {
                None | Some("2") => sqrt(&add(&vals.iter().map(|e| pow(&abs(e), &two)).collect::<Vec<_>>())),
                Some("1") => add(&vals.iter().map(abs).collect::<Vec<_>>()),
                Some("oo") => maxf(vals.iter().map(abs).collect()),
                Some("-oo") => minf(vals.iter().map(abs).collect()),
                Some(o) => {
                    let p = crate::parser::parse_str(o).map_err(|e| MatrixError::Value(e.to_string()))?;
                    pow(&add(&vals.iter().map(|e| pow(&abs(e), &p)).collect::<Vec<_>>()), &pow(&p, &Expr::minus_one()))
                }
            });
        }
        match ord {
            Some("1") => {
                let m = self.applyfunc(abs);
                Ok(maxf((0..m.cols()).map(|j| add(m.col(j).data())).collect()))
            }
            Some("oo") => {
                let m = self.applyfunc(abs);
                Ok(maxf((0..m.rows()).map(|i| add(m.row(i).data())).collect()))
            }
            None | Some("f") | Some("fro") | Some("frobenius") | Some("vector") => self.vec().norm_ord(Some("2")),
            _ => Err(MatrixError::NotImplemented("Matrix Norms under development".into())),
        }
    }
}
