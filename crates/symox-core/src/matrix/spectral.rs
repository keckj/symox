//! Exact spectral quantities and Jordan chains, sharing the native eigensolver.
use super::ratfunc::is_zero_assume;
use super::{non_square, MResult, Matrix, MatrixError};
use crate::expr::{Expr, Func};
use crate::functions::func;
use crate::ops::{div, pow, sqrt};

impl Matrix {
    fn small_gram(&self) -> MResult<Matrix> {
        let h = self.adjoint();
        if self.rows() >= self.cols() {
            h.matmul(self)
        } else {
            self.matmul(&h)
        }
    }

    /// Singular values with multiplicities, padded to `cols` as in SymPy.
    pub fn singular_values(&self) -> MResult<Vec<Expr>> {
        let mut values = Vec::new();
        for (lambda, multiplicity) in self.small_gram()?.cancel().eigenvals()? {
            values.extend(std::iter::repeat_n(sqrt(&lambda), multiplicity));
        }
        values.resize(self.cols(), Expr::zero());
        values.sort_by(|a, b| crate::sortkey::cmp_default(b, a));
        Ok(values)
    }

    pub fn condition_number(&self) -> MResult<Expr> {
        if self.is_empty() {
            return Ok(Expr::zero());
        }
        let values = self.singular_values()?;
        Ok(div(
            &func(Func::Max, values.clone()),
            &func(Func::Min, values),
        ))
    }

    /// Spectral projector formula for the Moore-Penrose inverse.
    /// Normalization of eigenvectors is deliberately unnecessary: a Hermitian
    /// matrix's eigenspaces are orthogonal, and its projectors are invariant
    /// under the choice of basis inside any repeated eigenspace.
    pub fn pinv_ed(&self) -> MResult<Matrix> {
        if self.is_empty() || self.is_zero_matrix() {
            return Ok(Matrix::zeros(self.cols(), self.rows()));
        }
        let gram = self.small_gram()?.cancel();
        let (p, d) = gram.diagonalize_with_options(false, true, false)?;
        let inverse = Matrix::diag(
            &(0..d.rows())
                .map(|i| {
                    let value = d.get(i, i);
                    if is_zero_assume(value) == Some(true) {
                        Expr::zero()
                    } else {
                        pow(value, &Expr::minus_one())
                    }
                })
                .collect::<Vec<_>>(),
        );
        let gram_pinv = p.matmul(&inverse)?.matmul(&p.inv()?)?.cancel();
        let h = self.adjoint();
        if self.rows() >= self.cols() {
            gram_pinv.matmul(&h)
        } else {
            h.matmul(&gram_pinv)
        }
        .map(|m| m.cancel())
    }

    /// Jordan form with an exact change of basis. Blocks are ordered by
    /// eigenvalue and then decreasing size, as in SymPy. For each block size
    /// k, choose heads in ker(N^k)/(ker(N^(k-1)) + N ker(N^(k+1))). This
    /// quotient prevents dependent chains when several block sizes coexist.
    pub fn jordan_decomposition(&self) -> MResult<(Matrix, Matrix, Vec<Matrix>)> {
        if !self.is_square() {
            return Err(non_square());
        }
        let n = self.rows();
        if n == 0 {
            return Ok((Matrix::eye(0), Matrix::eye(0), Vec::new()));
        }
        let mut eigenvalues = self.eigenvals()?;
        eigenvalues.sort_by(|a, b| crate::sortkey::cmp_default(&a.0, &b.0));
        let mut columns = Vec::new();
        let mut blocks = Vec::new();
        for (value, multiplicity) in eigenvalues {
            let nil = self.sub(&Matrix::eye(n).scale(&value))?.cancel();
            let mut kernels: Vec<Vec<Matrix>> = vec![Vec::new()];
            let mut power = Matrix::eye(n);
            for _ in 1..=multiplicity {
                power = power.matmul(&nil)?.cancel();
                kernels.push(power.nullspace());
                if kernels.last().unwrap().len() == multiplicity {
                    break;
                }
            }
            if kernels.last().unwrap().len() != multiplicity {
                return Err(MatrixError::NotImplemented(format!(
                    "Cannot determine the generalized eigenspace for {value}"
                )));
            }
            let max_size = kernels.len() - 1;
            kernels.push(kernels[max_size].clone());
            for size in (1..=max_size).rev() {
                let count = 2 * kernels[size].len();
                let subtract = kernels[size - 1].len() + kernels[size + 1].len();
                if count < subtract {
                    return Err(MatrixError::NotImplemented(
                        "Inconsistent symbolic Jordan ranks".into(),
                    ));
                }
                let count = count - subtract;
                if count == 0 {
                    continue;
                }
                let mut excluded = kernels[size - 1].clone();
                for v in &kernels[size + 1] {
                    excluded.push(nil.matmul(v)?.cancel());
                }
                let mut space = column_stack(n, &excluded)?;
                let mut rank = space.rank();
                let mut found = 0;
                for head in &kernels[size] {
                    let extended = space.row_join(head)?;
                    let new_rank = extended.rank();
                    if new_rank == rank {
                        continue;
                    }
                    let mut chain = vec![head.clone()];
                    for _ in 1..size {
                        chain.push(nil.matmul(chain.last().unwrap())?.cancel());
                    }
                    columns.extend(chain.into_iter().rev());
                    blocks.push(Matrix::from_fn(size, size, |i, j| {
                        if i == j {
                            value.clone()
                        } else if j == i + 1 {
                            Expr::one()
                        } else {
                            Expr::zero()
                        }
                    }));
                    found += 1;
                    space = extended;
                    rank = new_rank;
                    if found == count {
                        break;
                    }
                }
                if found != count {
                    return Err(MatrixError::NotImplemented(
                        "Unable to choose independent symbolic Jordan chains".into(),
                    ));
                }
            }
        }
        let p = column_stack(n, &columns)?;
        let j = Matrix::diag_blocks(&blocks);
        if p.shape() != (n, n)
            || p.rank() != n
            || !self
                .matmul(&p)?
                .sub(&p.matmul(&j)?)?
                .cancel()
                .is_zero_matrix()
        {
            return Err(MatrixError::NotImplemented(
                "Unable to verify symbolic Jordan chains exactly".into(),
            ));
        }
        Ok((p, j, blocks))
    }
}

fn column_stack(rows: usize, columns: &[Matrix]) -> MResult<Matrix> {
    if columns.is_empty() {
        return Ok(Matrix::zeros(rows, 0));
    }
    Matrix::hstack(&columns.iter().collect::<Vec<_>>())
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
    fn spectral_and_ed() {
        let a = m(&[&[1, 2, 3], &[2, 4, 6]]);
        assert_eq!(a.pinv_ed().unwrap(), a.pinv().unwrap());
        let d = m(&[&[3, 0], &[0, 4]]);
        assert_eq!(
            d.singular_values().unwrap(),
            vec![Expr::int(4), Expr::int(3)]
        );
        assert_eq!(
            d.condition_number().unwrap(),
            crate::ops::div(&Expr::int(4), &Expr::int(3))
        );
    }
    #[test]
    fn jordan_mixed_chain_lengths() {
        let block = |size, value| {
            Matrix::from_fn(size, size, |i, j| {
                if i == j {
                    Expr::int(value)
                } else if j == i + 1 {
                    Expr::one()
                } else {
                    Expr::zero()
                }
            })
        };
        let j = Matrix::diag_blocks(&[block(3, 2), block(2, 2), block(1, 2), block(2, -1)]);
        let p = Matrix::from_fn(8, 8, |i, j| {
            if i == j {
                Expr::one()
            } else if j == i + 1 {
                Expr::int((i + 1) as i64)
            } else {
                Expr::zero()
            }
        });
        let a = p.matmul(&j).unwrap().matmul(&p.inv().unwrap()).unwrap();
        let (q, k, blocks) = a.jordan_decomposition().unwrap();
        assert_eq!(
            blocks.iter().map(|b| b.rows()).collect::<Vec<_>>(),
            vec![2, 3, 2, 1]
        );
        assert_eq!(a.matmul(&q).unwrap(), q.matmul(&k).unwrap());
        assert_eq!(q.rank(), 8);
    }
}
