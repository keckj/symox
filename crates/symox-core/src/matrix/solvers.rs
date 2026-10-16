//! Exact rank factorizations, generalized inverses and structured solves.
use super::ratfunc::{dotprodsimp, is_zero_assume};
use super::{non_invertible, non_square, MResult, Matrix, MatrixError};
use crate::expr::Expr;
use crate::ops::{add, div, mul2, sub};
use std::sync::atomic::{AtomicU64, Ordering};

impl Matrix {
    fn check_rhs(&self, rhs: &Matrix) -> MResult<()> {
        if self.rows() != rhs.rows() {
            return Err(MatrixError::Shape(
                "Matrices must have the same number of rows.".into(),
            ));
        }
        Ok(())
    }

    /// Full rank factorization `A = C F`, with the original pivot columns in C.
    pub fn rank_decomposition(&self) -> MResult<(Matrix, Matrix)> {
        let (rref, pivots) = self.rref();
        let c = self.extract(&(0..self.rows()).collect::<Vec<_>>(), &pivots)?;
        let f = rref.slice(0, pivots.len(), 0, self.cols());
        Ok((c, f))
    }

    /// Moore-Penrose inverse by a full rank decomposition, over exact expressions.
    pub fn pinv(&self) -> MResult<Matrix> {
        let (c, f) = self.rank_decomposition()?;
        if c.cols() == 0 {
            return Ok(Matrix::zeros(self.cols(), self.rows()));
        }
        let ch = c.adjoint();
        let fh = f.adjoint();
        // F+ C+ = F.H (F F.H)^-1 (C.H C)^-1 C.H.
        fh.matmul(&f.matmul(&fh)?.inv()?)?
            .matmul(&ch.matmul(&c)?.inv()?)?
            .matmul(&ch)
            .map(|m| m.cancel())
    }

    /// All least-squares solutions, with arbitrary null-space components.
    pub fn pinv_solve(&self, rhs: &Matrix, arbitrary: Option<&Matrix>) -> MResult<Matrix> {
        self.check_rhs(rhs)?;
        if arbitrary.is_some_and(|w| w.shape() != (self.cols(), rhs.cols())) {
            return Err(MatrixError::Shape(
                "Arbitrary matrix must have shape (A.cols, rhs.cols).".into(),
            ));
        }
        let p = self.pinv()?;
        let particular = p.matmul(rhs)?;
        let projector = Matrix::eye(self.cols()).sub(&p.matmul(self)?)?.cancel();
        if projector.is_zero_matrix() {
            return Ok(particular);
        }
        let generated;
        let w = match arbitrary {
            Some(w) => w,
            None => {
                // Expressions currently identify symbols by name; fresh suffixes give
                // independent arbitrary matrices on successive calls as Dummy does.
                static NEXT: AtomicU64 = AtomicU64::new(0);
                let used = [self.free_symbols(), rhs.free_symbols()].concat();
                loop {
                    let id = NEXT.fetch_add(1, Ordering::Relaxed);
                    let candidate = Matrix::from_fn(self.cols(), rhs.cols(), |i, j| {
                        Expr::symbol(&format!("_w{i}_{j}_{id}"))
                    });
                    if !candidate.data().iter().any(|e| used.contains(e)) {
                        generated = candidate;
                        break;
                    }
                }
                &generated
            }
        };
        particular.add(&projector.matmul(w)?)
    }

    pub fn lower_triangular_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        self.triangular_solve(rhs, true)
    }

    pub fn upper_triangular_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        self.triangular_solve(rhs, false)
    }

    fn triangular_solve(&self, rhs: &Matrix, lower: bool) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(non_square());
        }
        self.check_rhs(rhs)?;
        if lower && !self.is_lower() {
            return Err(MatrixError::Value(
                "Matrix must be lower triangular.".into(),
            ));
        }
        if !lower && !self.is_upper() {
            return Err(MatrixError::Type("Matrix is not upper triangular.".into()));
        }
        let n = self.rows();
        let mut result = Matrix::zeros(n, rhs.cols());
        for j in 0..rhs.cols() {
            for pos in 0..n {
                let i = if lower { pos } else { n - 1 - pos };
                if is_zero_assume(self.get(i, i)) == Some(true) {
                    return Err(if lower {
                        MatrixError::Type("Matrix must be non-singular.".into())
                    } else {
                        MatrixError::Value("Matrix must be non-singular.".into())
                    });
                }
                let indices = if lower { 0..i } else { i + 1..n };
                let terms: Vec<_> = indices
                    .map(|k| mul2(self.get(i, k), result.get(k, j)))
                    .collect();
                result.set(
                    i,
                    j,
                    dotprodsimp(&div(&sub(rhs.get(i, j), &add(&terms)), self.get(i, i))),
                );
            }
        }
        Ok(result)
    }

    pub fn diagonal_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        if !self.is_diagonal() {
            return Err(MatrixError::Type("Matrix should be diagonal".into()));
        }
        if self.rows() != rhs.rows() {
            return Err(MatrixError::Type("Size mismatch".into()));
        }
        if !self.is_square() {
            return Err(non_square());
        }
        Ok(Matrix::from_fn(rhs.rows(), rhs.cols(), |i, j| {
            div(rhs.get(i, j), self.get(i, i))
        }))
    }

    /// Least squares for full-column-rank systems, using the Hermitian QR factor.
    pub fn qr_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        self.check_rhs(rhs)?;
        if self.rows() < self.cols() {
            return Err(MatrixError::NotImplemented(
                "Underdetermined systems require pinv_solve or gauss_jordan_solve.".into(),
            ));
        }
        // Keep the orthogonal columns unnormalized. This is the same QR
        // solve with their squared norms carried separately, and avoids
        // introducing square roots and conjugates of square roots into
        // exact symbolic systems.
        let n = self.cols();
        let mut q = self.clone();
        let mut r = Matrix::eye(n);
        let mut norms = Vec::with_capacity(n);
        for j in 0..n {
            for i in 0..j {
                let projection = q.col(i).adjoint().matmul(&q.col(j))?;
                let coefficient = dotprodsimp(&div(projection.get(0, 0), &norms[i]));
                r.set(i, j, coefficient.clone());
                for k in 0..self.rows() {
                    q.set(
                        k,
                        j,
                        dotprodsimp(&sub(q.get(k, j), &mul2(q.get(k, i), &coefficient))),
                    );
                }
            }
            let column = q.col(j);
            let square = dotprodsimp(column.adjoint().matmul(&column)?.get(0, 0));
            if is_zero_assume(&square) == Some(true) {
                return Err(non_invertible());
            }
            norms.push(square);
        }
        let projected = q.adjoint().matmul(rhs)?;
        let y = Matrix::from_fn(n, rhs.cols(), |i, j| {
            dotprodsimp(&div(projected.get(i, j), &norms[i]))
        });
        r.upper_triangular_solve(&y).map(|m| m.cancel())
    }

    /// Structured solve, reformulating general matrices as Hermitian normal equations.
    fn factor_solve(&self, rhs: &Matrix, ldl: bool, reform: bool) -> MResult<Matrix> {
        self.check_rhs(rhs)?;
        if self.rows() < self.cols() {
            return Err(MatrixError::NotImplemented(
                "Underdetermined systems require pinv_solve or gauss_jordan_solve.".into(),
            ));
        }
        if self.cols() == 0 {
            return Ok(Matrix::zeros(0, rhs.cols()));
        }
        if !self.is_square() || (!self.is_symmetric() && self.is_hermitian() != Some(true)) {
            if reform {
                return Err(MatrixError::NotImplemented(
                    "Unable to establish Hermitian normal equations for these symbolic entries."
                        .into(),
                ));
            }
            let h = self.adjoint();
            let normal = h.matmul(self)?.cancel();
            if normal.rank() < self.cols() {
                return Err(non_invertible());
            }
            return normal.factor_solve(&h.matmul(rhs)?, ldl, true);
        }
        if self.rank() < self.cols() {
            return Err(non_invertible());
        }
        let hermitian = !self.is_symmetric();
        let decomposition = if ldl {
            self.ldl_decomposition(hermitian).map(|(l, d)| (l, Some(d)))
        } else {
            self.cholesky(hermitian).map(|l| (l, None))
        }
        .and_then(|(l, d)| {
            let diagonal = d.as_ref().unwrap_or(&l);
            let mut valid =
                (0..diagonal.rows()).all(|i| is_zero_assume(diagonal.get(i, i)) != Some(true));
            for entry in l.data().iter().chain(d.iter().flat_map(|m| m.data())) {
                crate::traverse::preorder(entry, &mut |e| {
                    valid &= !e.is_constant(crate::expr::Constant::NaN)
                        && !e.is_constant(crate::expr::Constant::ComplexInfinity);
                });
            }
            if valid {
                Ok((l, d))
            } else {
                Err(non_invertible())
            }
        });
        let (l, d) = match decomposition {
            Ok(result) => result,
            // A Hermitian indefinite matrix need not admit this unpivoted
            // factorization. Its Gram matrix does when it has full column rank.
            Err(_) if !reform => {
                let h = self.adjoint();
                return h
                    .matmul(self)?
                    .cancel()
                    .factor_solve(&h.matmul(rhs)?, ldl, true);
            }
            Err(e) => return Err(e),
        };
        if !ldl
            && hermitian
            && (0..l.rows()).any(|i| crate::assumptions::is_real(l.get(i, i)) != Some(true))
        {
            // Square-root-free Cholesky avoids unprovable conjugates of
            // symbolic diagonal roots while retaining the same exact solve.
            return self.factor_solve(rhs, true, reform);
        }
        let y = l.lower_triangular_solve(rhs)?;
        let z = if let Some(d) = d {
            d.diagonal_solve(&y)?
        } else {
            y
        };
        let upper = if hermitian {
            l.adjoint()
        } else {
            l.transpose()
        };
        upper.upper_triangular_solve(&z).map(|m| m.cancel())
    }

    pub fn cholesky_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        self.factor_solve(rhs, false, false)
    }

    pub fn ldl_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        self.factor_solve(rhs, true, false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(rows: &[&[i64]]) -> Matrix {
        Matrix::from_rows(
            rows.iter()
                .map(|r| r.iter().map(|v| Expr::int(*v)).collect())
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn rank_and_penrose_identities() {
        let a = m(&[&[1, 2, 3], &[2, 4, 6]]);
        let (c, f) = a.rank_decomposition().unwrap();
        assert_eq!(c.shape(), (2, 1));
        assert_eq!(c.matmul(&f).unwrap(), a);
        let p = a.pinv().unwrap();
        assert_eq!(a.matmul(&p).unwrap().matmul(&a).unwrap(), a);
        assert_eq!(p.matmul(&a).unwrap().matmul(&p).unwrap(), p);
        let ap = a.matmul(&p).unwrap();
        let pa = p.matmul(&a).unwrap();
        assert_eq!(ap.adjoint(), ap);
        assert_eq!(pa.adjoint(), pa);
    }

    #[test]
    fn structured_solves() {
        let a = m(&[&[1, 0], &[1, 1], &[1, 2]]);
        let b = m(&[&[1, 2], &[2, 4], &[2, 3]]);
        for x in [a.qr_solve(&b), a.cholesky_solve(&b), a.ldl_solve(&b)] {
            let residual = a.matmul(&x.unwrap()).unwrap().sub(&b).unwrap();
            assert!(a
                .adjoint()
                .matmul(&residual)
                .unwrap()
                .cancel()
                .is_zero_matrix());
        }
        let zero = Matrix::zeros(0, 3);
        assert_eq!(zero.pinv().unwrap().shape(), (3, 0));
        assert_eq!(zero.rank_decomposition().unwrap().0.shape(), (0, 0));
    }

    #[test]
    fn complex_hermitian_inverse_and_least_squares() {
        let a = Matrix::from_rows(vec![
            vec![Expr::one(), Expr::i()],
            vec![Expr::i(), Expr::minus_one()],
        ])
        .unwrap();
        let p = a.pinv().unwrap();
        assert!(a
            .matmul(&p)
            .unwrap()
            .matmul(&a)
            .unwrap()
            .sub(&a)
            .unwrap()
            .cancel()
            .is_zero_matrix());
        let ap = a.matmul(&p).unwrap();
        let pa = p.matmul(&a).unwrap();
        assert!(ap.adjoint().sub(&ap).unwrap().cancel().is_zero_matrix());
        assert!(pa.adjoint().sub(&pa).unwrap().cancel().is_zero_matrix());
        let full = Matrix::from_rows(vec![
            vec![Expr::one(), Expr::i()],
            vec![Expr::zero(), Expr::int(2)],
            vec![Expr::int(2), Expr::one()],
        ])
        .unwrap();
        let b = m(&[&[1], &[2], &[3]]);
        for (method, solution) in [
            ("QR", full.qr_solve(&b)),
            ("CH", full.cholesky_solve(&b)),
            ("LDL", full.ldl_solve(&b)),
        ] {
            let residual = full.matmul(&solution.unwrap()).unwrap().sub(&b).unwrap();
            let normal = full.adjoint().matmul(&residual).unwrap().cancel();
            assert!(normal.is_zero_matrix(), "{method}: {normal}");
        }
    }
}
