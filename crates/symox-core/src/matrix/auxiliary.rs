//! SymPy `MatrixBase` utilities that do not belong to a larger algorithm
//! family: triangular parts, rotations, Berkowitz vectors, permanents,
//! fraction-free LU, Householder reductions and small classic constructors.
//!
//! The algorithms follow SymPy's step order so that symbolic results have the
//! same structure (see LICENSE-SYMPY). Exponential or cubic kernels check a
//! conservative size bound and raise `NotImplemented` instead of running away.
use super::ratfunc::is_zero_assume;
use super::{MResult, Matrix, MatrixError};
use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::number::Number;
use crate::ops::{add, div, mul2, neg, pow, sub};
use num_bigint::BigInt;

/// Largest square size accepted by the exponential Ryser permanent.
pub const MAX_PERMANENT: usize = 20;

fn non_square(message: &str) -> MatrixError {
    MatrixError::NonSquare(message.into())
}

fn resource(name: &str) -> MatrixError {
    MatrixError::NotImplemented(format!(
        "{name}: matrix size exceeds the native resource limit"
    ))
}

fn compact(e: &Expr) -> Expr {
    let cancelled = crate::poly::api::cancel(e).unwrap_or_else(|_| e.clone());
    crate::simplify::radsimp(&cancelled)
}

fn check_shape(name: &str, rows: usize, cols: usize) -> MResult<()> {
    if rows
        .checked_mul(cols)
        .is_none_or(|n| n > super::structural::MAX_ITEMS)
    {
        Err(resource(name))
    } else {
        Ok(())
    }
}

/// Structural truthiness of SymPy entries (`bool(expr)`): only zero is false.
fn truthy(e: &Expr) -> bool {
    !e.is_zero()
}

/// `LeviCivita(*args)` for integer arguments (SymPy's `eval_levicivita`):
/// zero on repeated indices, otherwise prod(a_j - a_i) / prod(i!).
fn levi_civita(args: &[i128]) -> Expr {
    let mut numerator: i128 = 1;
    let mut denominator: i128 = 1;
    for i in 0..args.len() {
        for j in i + 1..args.len() {
            numerator *= args[j] - args[i];
        }
        denominator *= (1..=i as i128).product::<i128>();
    }
    if numerator == 0 {
        return Expr::zero();
    }
    div(
        &Expr::num(Number::from_bigint(BigInt::from(numerator))),
        &Expr::num(Number::from_bigint(BigInt::from(denominator))),
    )
}

impl Matrix {
    /// Entries on or above (`upper`) / below the `k`-th diagonal.
    pub fn triangular_part(&self, k: i64, upper: bool) -> Matrix {
        Matrix::from_fn(self.rows(), self.cols(), |i, j| {
            let offset = j as i128 - i as i128;
            let keep = if upper {
                offset >= k as i128
            } else {
                offset <= k as i128
            };
            if keep {
                self.get(i as usize, j as usize).clone()
            } else {
                Expr::zero()
            }
        })
    }

    /// Rotate by `k` quarter turns clockwise (SymPy's `rot90`).
    pub fn rot90(&self, k: i64) -> Matrix {
        let (r, c) = self.shape();
        match k.rem_euclid(4) {
            0 => self.clone(),
            1 => Matrix::from_fn(c, r, |i, j| self.get(r - 1 - j, i).clone()),
            2 => Matrix::from_fn(r, c, |i, j| self.get(r - 1 - i, c - 1 - j).clone()),
            _ => Matrix::from_fn(c, r, |i, j| self.get(j, c - 1 - i).clone()),
        }
    }

    /// Half-vectorization of a square matrix, column by column below the
    /// diagonal (inclusive when `diagonal`). The caller checks symmetry.
    pub fn vech(&self, diagonal: bool) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(non_square("Matrix must be square"));
        }
        let c = self.cols();
        let mut v = Vec::new();
        for j in 0..c {
            for i in (if diagonal { j } else { j + 1 })..c {
                v.push(self.get(i, j).clone());
            }
        }
        Ok(Matrix::column(v))
    }

    /// Skew-symmetric cross-product matrix of a 3-vector.
    pub fn hat(&self) -> MResult<Matrix> {
        if self.shape() != (3, 1) {
            return Err(MatrixError::Shape(format!(
                "Dimensions incorrect, expected (3, 1), got ({}, {})",
                self.rows(),
                self.cols()
            )));
        }
        let (x, y, z) = (self.get(0, 0), self.get(1, 0), self.get(2, 0));
        Matrix::new(
            3,
            3,
            vec![
                Expr::zero(),
                neg(z),
                y.clone(),
                z.clone(),
                Expr::zero(),
                neg(x),
                neg(y),
                x.clone(),
                Expr::zero(),
            ],
        )
    }

    /// Inverse of `hat`; requires a 3x3 skew-symmetric matrix.
    pub fn vee(&self) -> MResult<Matrix> {
        if self.shape() != (3, 3) {
            return Err(MatrixError::Shape(format!(
                "Dimensions incorrect, expected (3, 3), got ({}, {})",
                self.rows(),
                self.cols()
            )));
        }
        if self.is_anti_symmetric_with(&mut |e| Ok(crate::simplify::simplify(e)))? != Some(true) {
            return Err(MatrixError::Value("Matrix is not skew-symmetric".into()));
        }
        Matrix::new(
            3,
            1,
            vec![
                self.get(2, 1).clone(),
                self.get(0, 2).clone(),
                self.get(1, 0).clone(),
            ],
        )
    }

    /// Three-valued anti-symmetry test over the nonzero entries, with
    /// `simplify` applied to every `M[i, j] + M[j, i]`.
    pub fn is_anti_symmetric_with(
        &self,
        simplify: &mut dyn FnMut(&Expr) -> MResult<Expr>,
    ) -> MResult<Option<bool>> {
        if !self.is_square() {
            return Ok(Some(false));
        }
        let mut result = Some(true);
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                if !truthy(self.get(i, j)) {
                    continue;
                }
                let s = simplify(&add(&[self.get(i, j).clone(), self.get(j, i).clone()]))?;
                match crate::assumptions::is_zero(&s) {
                    Some(true) => {}
                    Some(false) => return Ok(Some(false)),
                    None => result = None,
                }
            }
        }
        Ok(result)
    }

    /// Dual of an antisymmetric matrix (SymPy's `dual`, LeviCivita based).
    pub fn dual(&self) -> MResult<Matrix> {
        // SymPy sizes the result by the row count and reads the leading
        // block without a shape check.
        let n = self.rows();
        let mut work = Matrix::zeros(n, n);
        if self.is_square() && self.is_symmetric() {
            return Ok(work);
        }
        if n >= 2 && self.cols() < n {
            return Err(MatrixError::Index(format!(
                "Index out of range: a[{}]",
                self.cols().max(1)
            )));
        }
        if n > 512 {
            return Err(resource("dual"));
        }
        for i in 1..n {
            for j in 1..n {
                let terms: Vec<Expr> = (1..n)
                    .map(|k| {
                        mul2(
                            &levi_civita(&[i as i128, j as i128, 0, k as i128]),
                            self.get(0, k),
                        )
                    })
                    .collect();
                let acum = add(&terms);
                work.set(i, j, acum.clone());
                work.set(j, i, neg(&acum));
            }
        }
        for l in 1..n {
            let mut terms = Vec::new();
            for a in 1..n {
                for b in 1..n {
                    terms.push(mul2(
                        &levi_civita(&[0, l as i128, a as i128, b as i128]),
                        self.get(a, b),
                    ));
                }
            }
            let acum = div(&add(&terms), &Expr::int(2));
            work.set(0, l, neg(&acum));
            work.set(l, 0, acum);
        }
        Ok(work)
    }

    /// Permanent by Ryser's formula, in SymPy's subset order, simplified.
    pub fn permanent(&self) -> MResult<Expr> {
        let mut m = self.clone();
        if m.rows() > m.cols() {
            m = m.transpose();
        }
        let (rows, n) = m.shape();
        if n > MAX_PERMANENT {
            return Err(resource("per"));
        }
        let binomial = |a: usize, b: usize| -> Expr {
            if b > a {
                return Expr::zero();
            }
            let mut value = BigInt::from(1);
            for t in 0..b {
                value = value * (a - t) / (t + 1);
            }
            Expr::num(Number::from_bigint(value))
        };
        let mut terms = Vec::new();
        for size in 1..=rows {
            for subset in combinations(n, size) {
                let mut product = Vec::with_capacity(rows + 1);
                for i in 0..rows {
                    product.push(add(&subset
                        .iter()
                        .map(|&j| m.get(i, j).clone())
                        .collect::<Vec<_>>()));
                }
                let sign = if size % 2 == 0 {
                    Expr::one()
                } else {
                    Expr::minus_one()
                };
                product.push(sign);
                product.push(binomial(n - size, rows - size));
                terms.push(crate::ops::mul(&product));
            }
        }
        let mut perm = add(&terms);
        if rows % 2 == 1 {
            perm = neg(&perm);
        }
        Ok(crate::simplify::simplify(&perm))
    }

    /// SymPy's `berkowitz()` vectors: characteristic polynomial coefficient
    /// lists of the leading principal submatrices, starting with `(1,)`.
    pub fn berkowitz_vectors(&self) -> MResult<Vec<Vec<Expr>>> {
        let mut berk = vec![vec![Expr::one()]];
        if self.data().is_empty() {
            return Ok(berk);
        }
        if !self.is_square() {
            return Err(non_square("Matrix must be square"));
        }
        let n_total = self.rows();
        if n_total > 512 {
            return Err(resource("berkowitz"));
        }
        let mut a = self.clone();
        let mut transforms: Vec<Matrix> = vec![Matrix::zeros(0, 0); n_total.saturating_sub(1)];
        for n in (2..=n_total).rev() {
            let k = n - 1;
            let r = a.slice(k, k + 1, 0, k).neg();
            let c = a.slice(0, k, k, k + 1);
            let diag = neg(a.get(k, k));
            a = a.slice(0, k, 0, k);
            let mut items = vec![c];
            for i in 0..n - 2 {
                let next = a.matmul(&items[i])?;
                items.push(next);
            }
            let mut values = vec![Expr::one(), diag];
            for b in &items {
                values.push(r.matmul(b)?.get(0, 0).clone());
            }
            let mut t = Matrix::zeros(n + 1, n);
            for i in 0..n {
                for (offset, value) in values[..n - i + 1].iter().enumerate() {
                    t.set(i + offset, i, value.clone());
                }
            }
            transforms[k - 1] = t;
        }
        let mut polys = vec![Matrix::column(vec![Expr::one(), neg(a.get(0, 0))])];
        for (i, t) in transforms.iter().enumerate() {
            let next = t.matmul(&polys[i])?;
            polys.push(next);
        }
        berk.extend(polys.into_iter().map(|p| p.data().to_vec()));
        Ok(berk)
    }

    /// Signed constant coefficients of the Berkowitz vectors.
    pub fn berkowitz_minors(&self) -> MResult<Vec<Expr>> {
        Ok(self
            .berkowitz_vectors()?
            .into_iter()
            .enumerate()
            .map(|(i, p)| {
                let last = p.last().cloned().unwrap_or_else(Expr::zero);
                if i % 2 == 0 {
                    last
                } else {
                    neg(&last)
                }
            })
            .collect())
    }

    /// Fraction-free LU decomposition `P*A = L*D**-1*U` (SymPy's
    /// `LUdecompositionFF`), returning `(P, L, DD, U)`.
    pub fn lu_decomposition_ff(&self) -> MResult<(Matrix, Matrix, Matrix, Matrix)> {
        let (n, m) = self.shape();
        if n > 2048 || m > 2048 {
            return Err(resource("LUdecompositionFF"));
        }
        let mut u = self.clone();
        let mut l = Matrix::eye(n);
        let mut p = Matrix::eye(n);
        let mut dd = Matrix::zeros(n, n);
        let mut oldpivot = Expr::one();
        for k in 0..n.saturating_sub(1) {
            if k >= m {
                return Err(MatrixError::Index(format!("Index out of range: a[{k}]")));
            }
            if u.get(k, k).is_zero() {
                let kpivot = (k + 1..n)
                    .find(|&r| truthy(u.get(r, k)))
                    .ok_or_else(|| MatrixError::Value("Matrix is not full rank".into()))?;
                for j in k..m {
                    let (a, b) = (u.get(k, j).clone(), u.get(kpivot, j).clone());
                    u.set(k, j, b);
                    u.set(kpivot, j, a);
                }
                for j in 0..k {
                    let (a, b) = (l.get(k, j).clone(), l.get(kpivot, j).clone());
                    l.set(k, j, b);
                    l.set(kpivot, j, a);
                }
                for j in 0..n {
                    let (a, b) = (p.get(k, j).clone(), p.get(kpivot, j).clone());
                    p.set(k, j, b);
                    p.set(kpivot, j, a);
                }
            }
            let ukk = u.get(k, k).clone();
            l.set(k, k, ukk.clone());
            dd.set(k, k, mul2(&oldpivot, &ukk));
            for i in k + 1..n {
                let uik = u.get(i, k).clone();
                l.set(i, k, uik.clone());
                for j in k + 1..m {
                    let v = div(
                        &sub(&mul2(&ukk, u.get(i, j)), &mul2(u.get(k, j), &uik)),
                        &oldpivot,
                    );
                    u.set(i, j, v);
                }
                u.set(i, k, Expr::zero());
            }
            oldpivot = ukk;
        }
        if n > 0 {
            dd.set(n - 1, n - 1, oldpivot);
        }
        Ok((p, l, dd, u))
    }

    /// Combined `L\U` storage and row swaps (SymPy's `LUdecomposition_Simple`).
    pub fn lu_decomposition_simple(
        &self,
        rankcheck: bool,
    ) -> MResult<(Matrix, Vec<(usize, usize)>)> {
        if self.data().is_empty() {
            return Ok((Matrix::zeros(self.rows(), self.cols()), vec![]));
        }
        super::linalg::lu_simple(self, rankcheck)
    }

    /// Companion matrix of a monic polynomial given by all its coefficients,
    /// leading coefficient first.
    pub fn companion(coeffs: &[Expr]) -> MResult<Matrix> {
        let size = coeffs.len().saturating_sub(1);
        if size == 0 || !coeffs[0].is_one() {
            return Err(MatrixError::Value(
                "Companion matrices require a monic polynomial of positive degree".into(),
            ));
        }
        check_shape("companion", size, size)?;
        Ok(Matrix::from_fn(size, size, |i, j| {
            if j == size - 1 {
                neg(&coeffs[coeffs.len() - 1 - i])
            } else if i == j + 1 {
                Expr::one()
            } else {
                Expr::zero()
            }
        }))
    }

    /// Jordan block with `eigenvalue` on the diagonal and ones on the upper
    /// (or `lower`) band.
    pub fn jordan_block(size: usize, eigenvalue: &Expr, lower: bool) -> MResult<Matrix> {
        check_shape("jordan_block", size, size)?;
        Ok(Matrix::from_fn(size, size, |i, j| {
            if i == j {
                eigenvalue.clone()
            } else if (!lower && i + 1 == j) || (lower && j + 1 == i) {
                Expr::one()
            } else {
                Expr::zero()
            }
        }))
    }

    /// Wilkinson's tridiagonal eigenvalue test matrices `(W-, W+)`.
    pub fn wilkinson(n: usize) -> MResult<(Matrix, Matrix)> {
        let size = n
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| resource("wilkinson"))?;
        check_shape("wilkinson", size, size)?;
        let band = |i: usize, j: usize| i + 1 == j || j + 1 == i;
        let make = |absolute: bool| {
            Matrix::from_fn(size, size, |i, j| {
                if i == j {
                    let d = i as i64 - n as i64;
                    Expr::int(if absolute { d.abs() } else { d })
                } else if band(i, j) {
                    Expr::one()
                } else {
                    Expr::zero()
                }
            })
        };
        Ok((make(false), make(true)))
    }

    /// Whether the characteristic polynomial is `x**n`.
    pub fn is_nilpotent(&self) -> MResult<bool> {
        if self.data().is_empty() {
            return Ok(true);
        }
        if !self.is_square() {
            return Err(non_square("Nilpotency is valid only for square matrices"));
        }
        let cp = self.charpoly("_x")?;
        Ok(cp.coeffs.iter().skip(1).all(|c| c.is_zero()))
    }

    /// Whether any entry contains a symbol (SymPy's `has(Symbol)`).
    pub fn is_symbolic(&self) -> bool {
        let mut seen = rustc_hash::FxHashSet::default();
        let mut pending = self.data().to_vec();
        while let Some(e) = pending.pop() {
            if !seen.insert(e.clone()) {
                continue;
            }
            if matches!(e.kind(), Kind::Sym(_)) {
                return true;
            }
            pending.extend(e.args());
        }
        false
    }

    /// Leading diagonal blocks (SymPy's `get_diag_blocks`).
    pub fn diagonal_blocks(&self) -> Vec<Matrix> {
        let mut blocks = Vec::new();
        let mut m = self.clone();
        loop {
            let (r, c) = m.shape();
            let mut split = None;
            for i in 1..=r {
                let (right, bottom) = if i == 1 {
                    (
                        m.slice(0, 1.min(r), i.min(c), c),
                        m.slice(i, r, 0, 1.min(c)),
                    )
                } else {
                    (m.slice(0, i, i.min(c), c), m.slice(i, r, 0, i.min(c)))
                };
                if right.data().iter().any(truthy) || bottom.data().iter().any(truthy) {
                    continue;
                }
                split = Some(i);
                break;
            }
            let Some(i) = split else { break };
            let block = m.slice(0, i, 0, i.min(c));
            let done = block.shape() == m.shape();
            blocks.push(block);
            if done {
                break;
            }
            m = m.slice(i, r, i.min(c), c);
        }
        blocks
    }

    /// Connected components of the undirected graph of nonzero entries.
    pub fn connected_components(&self) -> MResult<Vec<Vec<usize>>> {
        if !self.is_square() {
            return Err(non_square(""));
        }
        let n = self.rows();
        let mut gmap = vec![Vec::new(); n];
        for i in 0..n {
            for j in 0..n {
                if truthy(self.get(i, j)) {
                    gmap[i].push(j);
                    gmap[j].push(i);
                }
            }
        }
        Ok(super::linalg::tarjan(&gmap))
    }

    /// Recursive 2x2 block inverse for size > 20, LU inverse otherwise.
    pub fn inverse_block(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(non_square("A Matrix must be square to invert."));
        }
        let i = self.rows();
        if i <= 20 {
            return self.inverse_lu();
        }
        let h = i / 2;
        let a = self.slice(0, h, 0, h);
        let b = self.slice(0, h, h, i);
        let c = self.slice(h, i, 0, h);
        let d = self.slice(h, i, h, i);
        let fallback = |e: MatrixError| match e {
            MatrixError::NonInvertible(_) => self.inverse_lu(),
            e => Err(e),
        };
        let d_inv = match d.inverse_block() {
            Ok(v) => v,
            Err(e) => return fallback(e),
        };
        let b_d_i = b.matmul(&d_inv)?;
        let bdc = b_d_i.matmul(&c)?;
        let a_n = match a.sub(&bdc)?.inverse_block() {
            Ok(v) => v,
            Err(e) => return fallback(e),
        };
        let b_n = a_n.neg().matmul(&b_d_i)?;
        let dc = d_inv.matmul(&c)?;
        let c_n = dc.neg().matmul(&a_n)?;
        let d_n = d_inv.add(&dc.matmul(&b_n.neg())?)?;
        let top = Matrix::hstack(&[&a_n, &b_n])?;
        let bottom = Matrix::hstack(&[&c_n, &d_n])?;
        Matrix::vstack(&[&top, &bottom])
    }

    /// `(H, P)` with `H` upper Hessenberg and `A = P*H*P.H` (Householder).
    pub fn upper_hessenberg_decomposition(&self) -> MResult<(Matrix, Matrix)> {
        if !self.is_square() {
            return Err(non_square("Matrix must be square."));
        }
        let n = self.cols();
        if n > 256 {
            return Err(resource("upper_hessenberg_decomposition"));
        }
        let mut p = Matrix::eye(n);
        let mut h = self.clone();
        let two = Expr::int(2);
        for j in 0..n.saturating_sub(2) {
            let mut u = h.slice(j + 1, n, j, j + 1);
            if u.slice(1, u.rows(), 0, 1)
                .data()
                .iter()
                .all(|e| is_zero_assume(e) == Some(true))
            {
                continue;
            }
            let s = func(Func::Sign, vec![u.get(0, 0).clone()]);
            let shifted = if !s.is_zero() {
                add(&[u.get(0, 0).clone(), mul2(&s, &u.norm())])
            } else {
                add(&[u.get(0, 0).clone(), u.norm()])
            };
            u.set(0, 0, shifted);
            // Use the unnormalized reflector: dividing by ||u|| twice
            // introduces nested square roots and much larger expression DAGs.
            let vh = u.adjoint();
            let den = crate::expand::expand(&add(&u
                .data()
                .iter()
                .zip(vh.data())
                .map(|(x, y)| mul2(x, y))
                .collect::<Vec<_>>()));
            let two_v = u.scale_left(&div(&two, &den));
            let lower = h.slice(j + 1, n, 0, n);
            let lower = lower.sub(&two_v.matmul(&vh.matmul(&lower)?)?)?;
            h.set_block(j + 1, 0, &lower);
            let right = h.slice(0, n, j + 1, n);
            let right = right.sub(&right.matmul(&two_v)?.matmul(&vh)?)?;
            h.set_block(0, j + 1, &right);
            let right = p.slice(0, n, j + 1, n);
            let right = right.sub(&right.matmul(&two_v)?.matmul(&vh)?)?;
            p.set_block(0, j + 1, &right);
            for i in j + 2..n {
                h.set(i, j, Expr::zero());
            }
        }
        Ok((h, p))
    }

    /// Householder bidiagonalization `A = U*B*V` (upper) as in SymPy;
    /// `transforms` additionally accumulates `U` and `V`.
    pub fn householder_bidiagonal(&self, transforms: bool) -> MResult<(Matrix, Matrix, Matrix)> {
        let (m, n) = self.shape();
        if m > 128 || n > 128 {
            return Err(resource("bidiagonalize"));
        }
        let mut a = self.clone();
        let mut u = Matrix::eye(m);
        let mut v = Matrix::eye(n);
        for i in 0..m.min(n) {
            let (hv, bet) = householder_vector(&a.slice(i, m, i, i + 1))?;
            let hh = Matrix::eye(m - i)
                .sub(&hv.scale_left(&bet).matmul(&hv.adjoint())?)?
                .applyfunc(compact);
            let block = hh.matmul(&a.slice(i, m, i, n))?.applyfunc(compact);
            a.set_block(i, i, &block);
            for row in i + 1..m {
                a.set(row, i, Expr::zero());
            }
            if transforms {
                let mut temp = Matrix::eye(m);
                temp.set_block(i, i, &hh);
                u = u.matmul(&temp)?.applyfunc(compact);
            }
            if i + 1 + 2 <= n {
                let (hv, bet) = householder_vector(&a.slice(i, i + 1, i + 1, n).adjoint())?;
                let hh = Matrix::eye(n - i - 1)
                    .sub(&hv.scale_left(&bet).matmul(&hv.adjoint())?)?
                    .applyfunc(compact);
                let block = a.slice(i, m, i + 1, n).matmul(&hh)?.applyfunc(compact);
                a.set_block(i, i + 1, &block);
                for col in i + 2..n {
                    a.set(i, col, Expr::zero());
                }
                if transforms {
                    let mut temp = Matrix::eye(n);
                    temp.set_block(i + 1, i + 1, &hh);
                    v = temp.matmul(&v)?.applyfunc(compact);
                }
            }
        }
        Ok((u, a, v))
    }

    /// Reduced row echelon forms of `[A | I | rhs]`, split as SymPy's
    /// `rref_rhs` does (including its whole-matrix result for empty `rhs`).
    pub fn rref_rhs(&self, rhs: &Matrix) -> MResult<(Matrix, Matrix)> {
        let eye = Matrix::eye(self.rows());
        let (r, _) = Matrix::hstack(&[self, &eye, rhs])?.rref();
        let left = r.slice(0, r.rows(), 0, self.cols());
        let right = if rhs.cols() == 0 {
            r.clone()
        } else {
            r.slice(0, r.rows(), r.cols() - rhs.cols(), r.cols())
        };
        Ok((left, right))
    }

    /// Cramer's rule, one determinant per unknown and right-hand side.
    pub fn cramer_solve(
        &self,
        rhs: &Matrix,
        det: &mut dyn FnMut(&Matrix) -> MResult<Expr>,
    ) -> MResult<Matrix> {
        if rhs.rows() != self.rows() {
            return Err(MatrixError::Shape("Matrix size mismatch.".into()));
        }
        let det_m = det(self)?;
        let mut x = Matrix::zeros(rhs.rows(), rhs.cols());
        for sol in 0..rhs.cols() {
            for col in 0..rhs.rows() {
                let replaced = Matrix::from_fn(self.rows(), self.cols(), |i, j| {
                    if j == col {
                        rhs.get(i, sol).clone()
                    } else {
                        self.get(i, j).clone()
                    }
                });
                x.set(col, sol, div(&det(&replaced)?, &det_m));
            }
        }
        Ok(x)
    }

    /// The vector divided by its norm, or zeros for a zero norm.
    pub fn normalized_with(
        &self,
        is_zero: &mut dyn FnMut(&Expr) -> MResult<bool>,
    ) -> MResult<Matrix> {
        if self.rows() != 1 && self.cols() != 1 {
            return Err(MatrixError::Shape(
                "A Matrix must be a vector to normalize.".into(),
            ));
        }
        let norm = self.norm();
        if is_zero(&norm)? {
            Ok(Matrix::zeros(self.rows(), self.cols()))
        } else {
            Ok(self.applyfunc(|e| div(e, &norm)))
        }
    }

    /// Projection of `self` onto `v`: `v * (self.dot(v) / v.dot(v))`.
    pub fn project(&self, v: &Matrix) -> MResult<Matrix> {
        let k = div(&self.dot(v)?, &v.dot(v)?);
        Ok(v.scale(&k))
    }

    /// Gram-Schmidt on vectors (SymPy's `orthogonalize`); dependent vectors
    /// are dropped unless `rankcheck` turns them into an error.
    pub fn orthogonalize(
        vecs: &[Matrix],
        normalize: bool,
        rankcheck: bool,
    ) -> MResult<Vec<Matrix>> {
        if vecs.is_empty() {
            return Ok(vec![]);
        }
        let row_vectors = vecs.iter().all(|v| v.rows() == 1);
        let columns: Vec<Matrix> = vecs.iter().map(|v| v.vec()).collect();
        let refs: Vec<&Matrix> = columns.iter().collect();
        let (q, _) = Matrix::hstack(&refs)?.qr_decomposition_optional(normalize)?;
        if rankcheck && q.cols() < vecs.len() {
            return Err(MatrixError::Value(
                "GramSchmidt: vector set not linearly independent".into(),
            ));
        }
        Ok((0..q.cols())
            .map(|i| {
                if row_vectors {
                    q.col(i).transpose()
                } else {
                    q.col(i)
                }
            })
            .collect())
    }

    /// Rows built from `ntop` matrices that fill the shape irregularly.
    pub fn irregular(ntop: usize, matrices: &[Matrix]) -> MResult<Matrix> {
        if ntop > matrices.len() {
            return Err(MatrixError::Index("pop from empty list".into()));
        }
        if ntop == 0 && matrices.iter().any(|m| m.rows() != 0) {
            return Err(MatrixError::Value(
                "A nonempty irregular matrix requires at least one top block".into(),
            ));
        }
        let mut queue: std::collections::VecDeque<usize> = (0..matrices.len()).collect();
        let mut remaining: Vec<usize> = matrices.iter().map(|m| m.rows()).collect();
        let mut active = Vec::with_capacity(ntop);
        for _ in 0..ntop {
            active.push(
                queue
                    .pop_front()
                    .ok_or_else(|| MatrixError::Index("pop from empty list".into()))?,
            );
        }
        let cols = active.iter().try_fold(0usize, |n, &i| {
            n.checked_add(matrices[i].cols())
                .ok_or_else(|| resource("irregular"))
        })?;
        let mut rows: Vec<Vec<Expr>> = Vec::new();
        if remaining.iter().any(|&r| r != 0) {
            check_shape("irregular", 1, cols)?;
        }
        while remaining.iter().any(|&r| r != 0) {
            let mut row = Vec::with_capacity(cols);
            for slot in 0..active.len() {
                let j = active[slot];
                let b = &matrices[j];
                if remaining[j] > 0 {
                    let r = b.rows() - remaining[j];
                    row.extend((0..b.cols()).map(|c| b.get(r, c).clone()));
                    remaining[j] -= 1;
                } else {
                    return Err(MatrixError::Index(
                        "Irregular blocks do not fill the matrix".into(),
                    ));
                }
                if remaining[j] == 0 {
                    if let Some(next) = queue.pop_front() {
                        active[slot] = next;
                    }
                }
            }
            if row.len() != cols {
                return Err(MatrixError::Value(
                    "Matrices provided do not appear to fill the space completely.".into(),
                ));
            }
            check_shape(
                "irregular",
                rows.len()
                    .checked_add(1)
                    .ok_or_else(|| resource("irregular"))?,
                cols,
            )?;
            rows.push(row);
        }
        if rows.is_empty() {
            return Ok(Matrix::zeros(0, 0));
        }
        Matrix::from_rows(rows)
    }

    /// Overwrite the block starting at `(r0, c0)` with `block`.
    pub fn set_block(&mut self, r0: usize, c0: usize, block: &Matrix) {
        for i in 0..block.rows() {
            for j in 0..block.cols() {
                self.set(r0 + i, c0 + j, block.get(i, j).clone());
            }
        }
    }
}

/// Householder vector and coefficient of a column (SymPy's
/// `_householder_vector`).
fn householder_vector(x: &Matrix) -> MResult<(Matrix, Expr)> {
    if x.rows() == 0 || x.cols() != 1 {
        return Err(MatrixError::Shape(
            "Householder input must be a nonempty column vector".into(),
        ));
    }
    let square = |v: &Matrix| {
        crate::expand::expand(&add(&v
            .data()
            .iter()
            .map(|e| mul2(&func(Func::Conjugate, vec![e.clone()]), e))
            .collect::<Vec<_>>()))
    };
    let mut v = x.clone();
    let norm_squared = square(x);
    if is_zero_assume(&norm_squared) == Some(true)
        || x.slice(1, x.rows(), 0, 1)
            .data()
            .iter()
            .all(|e| is_zero_assume(e) == Some(true))
    {
        v.set(0, 0, Expr::one());
        return Ok((v, Expr::zero()));
    }
    let first = x.get(0, 0);
    let phase = if is_zero_assume(first) == Some(true) {
        Expr::one()
    } else {
        func(Func::Sign, vec![first.clone()])
    };
    // Add the phase-aligned norm. Subtracting it can create an exact zero
    // pivot which is still held symbolically, followed by division by zero.
    let norm = crate::ops::sqrt(&norm_squared);
    v.set(0, 0, add(&[first.clone(), mul2(&phase, &norm)]));
    let beta = div(&Expr::int(2), &square(&v));
    Ok((v, beta))
}

/// Determinant of the Casoratian matrix `[seq_j(n + i)]` (or at `i` when
/// `zero`).
pub fn casoratian(seqs: &[Expr], n: &Expr, zero: bool) -> MResult<Expr> {
    let k = seqs.len();
    let m = Matrix::from_fn(k, k, |i, j| {
        let at = if zero {
            Expr::int(i as i64)
        } else {
            add(&[n.clone(), Expr::int(i as i64)])
        };
        crate::subs::subs(&seqs[j], &[(n.clone(), at)])
    });
    m.det_method(super::DetMethod::Bareiss)
}

/// Three-valued outcome of SymPy's `Expr.equals`.
#[derive(Clone, Debug, PartialEq)]
pub enum Equality {
    True,
    False,
    /// Undecided; carries the simplified difference (`failing_expression`).
    Unknown(Expr),
}

fn has_add_or_mod(e: &Expr) -> bool {
    match e.kind() {
        Kind::Add(..) => true,
        Kind::Func(f, _) if f.name() == "Mod" => true,
        _ => e.args().iter().any(has_add_or_mod),
    }
}

/// Conservative port of `Expr.equals`: simplification to zero proves
/// equality; a difference without sums, a nonzero constant, or two exact
/// evaluations that differ prove inequality; anything else is undecided.
pub fn expr_equals(a: &Expr, b: &Expr) -> Equality {
    if a == b {
        return Equality::True;
    }
    difference_is_zero(&crate::simplify::simplify(&sub(a, b)), 0)
}

fn difference_is_zero(diff: &Expr, depth: usize) -> Equality {
    if diff.is_zero() {
        return Equality::True;
    }
    if !has_add_or_mod(diff) {
        return Equality::False;
    }
    if let Kind::Mul(_, factors) = diff.kind() {
        if factors.len() > 1 && depth < 8 {
            let mut any_zero = false;
            let mut decided = true;
            for (base, exp) in factors {
                let factor = pow(base, exp);
                match difference_is_zero(&crate::simplify::simplify(&factor), depth + 1) {
                    Equality::True => any_zero = true,
                    Equality::False => {}
                    Equality::Unknown(_) => decided = false,
                }
            }
            if decided {
                return if any_zero {
                    Equality::True
                } else {
                    Equality::False
                };
            }
        }
    }
    let free = crate::traverse::free_symbols(diff);
    if free.is_empty() {
        return match crate::assumptions::is_zero(diff) {
            Some(true) => Equality::True,
            Some(false) => Equality::False,
            None => Equality::Unknown(diff.clone()),
        };
    }
    let evaluate = |value: &dyn Fn(usize) -> Expr| -> Option<Expr> {
        let pairs: Vec<(Expr, Expr)> = free
            .iter()
            .enumerate()
            .map(|(k, s)| (s.clone(), value(k)))
            .collect();
        let v = crate::simplify::simplify(&crate::subs::subs(diff, &pairs));
        (crate::traverse::free_symbols(&v).is_empty() && !v.has_nonfinite()).then_some(v)
    };
    let points = [
        evaluate(&|_| Expr::zero()),
        evaluate(&|k| Expr::int(k as i64 + 2)),
        evaluate(&|k| Expr::rational(1, k as i64 + 3)),
    ];
    let values: Vec<&Expr> = points.iter().flatten().collect();
    for i in 0..values.len() {
        for j in i + 1..values.len() {
            if crate::assumptions::is_zero(&sub(values[i], values[j])) == Some(false) {
                return Equality::False;
            }
        }
    }
    Equality::Unknown(diff.clone())
}

impl Matrix {
    /// Matrix `equals`: `False` on any unequal entry, else the first
    /// undecided entry result, else `True`.
    pub fn equals(&self, other: &Matrix) -> Equality {
        if self.shape() != other.shape() {
            return Equality::False;
        }
        let mut result = Equality::True;
        for (a, b) in self.data().iter().zip(other.data()) {
            match expr_equals(a, b) {
                Equality::False => return Equality::False,
                Equality::True => {}
                unknown => {
                    if result == Equality::True {
                        result = unknown;
                    }
                }
            }
        }
        result
    }
}

/// All `size`-subsets of `0..n` in lexicographic (itertools) order.
fn combinations(n: usize, size: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if size > n {
        return out;
    }
    let mut idx: Vec<usize> = (0..size).collect();
    loop {
        out.push(idx.clone());
        let mut i = size;
        loop {
            if i == 0 {
                return out;
            }
            i -= 1;
            if idx[i] != i + n - size {
                break;
            }
            if i == 0 {
                return out;
            }
        }
        idx[i] += 1;
        for j in i + 1..size {
            idx[j] = idx[j - 1] + 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(rows: usize, cols: usize, v: &[i64]) -> Matrix {
        Matrix::new(rows, cols, v.iter().map(|&x| Expr::int(x)).collect()).unwrap()
    }

    #[test]
    fn combinations_order() {
        assert_eq!(combinations(3, 2), vec![vec![0, 1], vec![0, 2], vec![1, 2]]);
        assert_eq!(combinations(2, 0), vec![Vec::<usize>::new()]);
        assert!(combinations(1, 2).is_empty());
    }

    #[test]
    fn permanent_berkowitz_and_rotations() {
        let a = m(3, 3, &[1, 2, 0, 3, 4, 5, 0, 6, 7]);
        assert_eq!(a.permanent().unwrap(), Expr::int(100));
        let berk = a.berkowitz_vectors().unwrap();
        assert_eq!(
            berk.last().unwrap(),
            &vec![Expr::int(1), Expr::int(-12), Expr::int(3), Expr::int(44)]
        );
        assert_eq!(
            a.berkowitz_minors().unwrap().last().unwrap(),
            &Expr::int(-44)
        );
        assert_eq!(a.rot90(4), a);
        assert_eq!(a.rot90(1).rot90(3), a);
        assert_eq!(a.rot90(1), m(3, 3, &[0, 3, 1, 6, 4, 2, 7, 5, 0]));
    }

    #[test]
    fn fraction_free_lu_reconstructs() {
        let a = m(3, 3, &[0, 2, 1, 3, 4, 5, 1, 6, 7]);
        let (p, l, dd, u) = a.lu_decomposition_ff().unwrap();
        let lhs = p.matmul(&a).unwrap();
        let rhs = l.matmul(&dd.inv().unwrap()).unwrap().matmul(&u).unwrap();
        assert_eq!(
            lhs.sub(&rhs)
                .unwrap()
                .data()
                .iter()
                .filter(|e| !e.is_zero())
                .count(),
            0
        );
        assert!(m(2, 2, &[0, 0, 0, 1]).lu_decomposition_ff().is_err());
    }

    #[test]
    fn levi_civita_matches_sympy() {
        assert_eq!(levi_civita(&[0, 1, 2, 3]), Expr::one());
        assert_eq!(levi_civita(&[1, 0, 2, 3]), Expr::minus_one());
        assert_eq!(levi_civita(&[0, 1, 2, 5]), Expr::int(10));
        assert_eq!(levi_civita(&[0, 0, 2, 3]), Expr::zero());
    }
}
