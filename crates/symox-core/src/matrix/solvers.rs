//! Exact rank factorizations, generalized inverses and structured solves.
use super::ratfunc::{dotprodsimp, is_zero_assume};
use super::{non_invertible, non_square, MResult, Matrix, MatrixError};
use crate::expr::Expr;
use crate::ops::{add, div, mul2, sub};

impl Matrix {
    /// A defined rational specialization of full column rank proves generic
    /// full rank for a rational-function matrix. Failed probes prove nothing:
    /// callers retain the original exact rank computation in that case.
    fn full_rank_specialization(&self) -> bool {
        use crate::expr::Kind;
        let mut budget = super::functions::budget::Budget::new();
        if budget.allocate(self.data().len()).is_err() {
            return false;
        }
        let mut pending = self
            .data()
            .iter()
            .cloned()
            .map(|entry| (entry, false))
            .collect::<Vec<_>>();
        let mut seen = rustc_hash::FxHashSet::default();
        let mut depths = rustc_hash::FxHashMap::default();
        while let Some((entry, expanded)) = pending.pop() {
            let children = super::functions::budget::raw_children(&entry);
            if expanded {
                let depth = 1 + children
                    .iter()
                    .map(|child| depths[&child.dag_id()])
                    .max()
                    .unwrap_or(0usize);
                // The original rank algorithm remains the fallback for deep
                // DAGs; do not enter recursive probe preflights on those DAGs.
                if depth > 128 {
                    return false;
                }
                depths.insert(entry.dag_id(), depth);
                continue;
            }
            if !seen.insert(entry.dag_id()) {
                continue;
            }
            if seen.len() > 100_000 {
                return false;
            }
            let rational = match entry.kind() {
                Kind::Num(number) => number.is_rational(),
                Kind::Sym(symbol) => {
                    symbol.facts().get("commutative") != Some(false)
                        && symbol.facts().get("finite") != Some(false)
                }
                Kind::Add(constant, terms) => {
                    constant.is_rational()
                        && terms
                            .iter()
                            .all(|(_, coefficient)| coefficient.is_rational())
                }
                Kind::Mul(coefficient, factors) => {
                    coefficient.is_rational()
                        && factors.iter().all(|(_, power)| power.is_integer_number())
                }
                Kind::Pow(_, power) => power.is_integer_number(),
                _ => false,
            };
            if !rational {
                return false;
            }
            pending.push((entry, true));
            pending.extend(children.into_iter().map(|child| (child, false)));
        }
        if budget.matrix(self).is_err() {
            return false;
        }
        let symbols = self.free_symbols();
        for attempt in 0..3u64 {
            let mut seed = 0x9e3779b97f4a7c15u64.wrapping_add(attempt);
            let mut replacements = rustc_hash::FxHashMap::default();
            for symbol in &symbols {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let value =
                    if crate::assumptions::exact_queries(|| crate::assumptions::is_zero(symbol))
                        == Some(true)
                    {
                        0
                    } else {
                        1 + ((seed >> 32) % 13) as i64
                    };
                replacements.insert(symbol.clone(), Expr::int(value));
            }
            // Probe only when substitution cannot force enormous held powers.
            // A local probe limit selects the original rank path, never an
            // unsupported result or a claim of rank deficiency.
            for entry in self.data() {
                if budget
                    .substitutions_preflight(entry, &replacements)
                    .is_err()
                {
                    return false;
                }
            }
            let specialized = self.applyfunc(|entry| crate::subs::xreplace(entry, &replacements));
            let Ok(height) = budget.matrix(&specialized) else {
                return false;
            };
            if budget.interpolation(self.cols(), height).is_err() {
                return false;
            }
            if specialized.data().iter().all(Expr::is_rational_number)
                && specialized.rank() == self.cols()
            {
                return true;
            }
        }
        false
    }

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
                // SymPy: symbols('w:{rows}_:{cols}', cls=Dummy), printed `_w0_0`.
                generated = Matrix::from_fn(self.cols(), rhs.cols(), |i, j| {
                    Expr::dummy(&format!("w{i}_{j}"), Default::default())
                });
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
        self.triangular_solve_form(rhs, lower, true)
    }

    fn triangular_solve_form(&self, rhs: &Matrix, lower: bool, simplify: bool) -> MResult<Matrix> {
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
                let value = div(&sub(rhs.get(i, j), &add(&terms)), self.get(i, i));
                result.set(i, j, if simplify { dotprodsimp(&value) } else { value });
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
        // Avoid an expanded symbolic RREF only after an exact full-rank
        // specialization certifies that its generic rank cannot be deficient.
        // Inconclusive probes retain the established rank/error path.
        let compact = self.data().iter().any(|value| !value.is_rational_number());
        if !self.is_square() || (!self.is_symmetric() && self.is_hermitian() != Some(true)) {
            if reform {
                return Err(MatrixError::NotImplemented(
                    "Unable to establish Hermitian normal equations for these symbolic entries."
                        .into(),
                ));
            }
            let h = self.adjoint();
            let normal = h.matmul(self)?;
            let normal = if compact { normal } else { normal.cancel() };
            if !compact && normal.rank() < self.cols() {
                return Err(non_invertible());
            }
            return normal.factor_solve(&h.matmul(rhs)?, ldl, true);
        }
        if (!compact || !self.full_rank_specialization()) && self.rank() < self.cols() {
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
            let mut pending = l
                .data()
                .iter()
                .chain(d.iter().flat_map(|m| m.data()))
                .cloned()
                .collect::<Vec<_>>();
            let mut visited = rustc_hash::FxHashSet::default();
            while let Some(entry) = pending.pop() {
                if !visited.insert(entry.dag_id()) {
                    continue;
                }
                valid &= !entry.is_constant(crate::expr::Constant::NaN)
                    && !entry.is_constant(crate::expr::Constant::ComplexInfinity);
                pending.extend(entry.args());
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
                let normal = h.matmul(self)?;
                return (if compact { normal } else { normal.cancel() }).factor_solve(
                    &h.matmul(rhs)?,
                    ldl,
                    true,
                );
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
        let y = l.triangular_solve_form(rhs, true, !compact)?;
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
        let solved = upper.triangular_solve_form(&z, false, !compact)?;
        Ok(if compact { solved } else { solved.cancel() })
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
