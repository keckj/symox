//! Linear algebra on [`Matrix`]: determinants, inverses, row reduction,
//! subspaces, LU, linear solvers and characteristic polynomials.
//!
//! Every operation dispatches on the ground domain of the entries:
//! exact rationals and polynomial matrices use fraction-free algorithms
//! (whose results coincide with SymPy's canonical outputs), the rest uses
//! faithful ports of SymPy's generic expression algorithms.

use super::domain::{matrix_domain, Domain, PolyMat};
use super::ff::{bareiss_det, berkowitz, rref_den, Ring};
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::ratfunc::{
    cancel, cancel_polys, dotprodsimp, is_zero_assume, is_zero_expand, PolyCtx,
};
use super::{non_invertible, non_square, MResult, Matrix, MatrixError};
use crate::expand::expand;
use crate::expr::{Expr, Kind};
use crate::number::Number;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};

pub use super::domain::RatMat;

/// Smallest dimension for which the multimodular kernels are used.
const MODULAR_MIN: usize = 6;

/// Determinant algorithms (SymPy's `method=` keyword).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DetMethod {
    Bareiss,
    Berkowitz,
    LU,
    DomainGE,
}

impl DetMethod {
    pub fn from_name(s: &str) -> MResult<DetMethod> {
        Ok(match s.to_lowercase().as_str() {
            "bareiss" | "bareis" => DetMethod::Bareiss,
            "berkowitz" => DetMethod::Berkowitz,
            "lu" | "det_lu" => DetMethod::LU,
            "domain-ge" => DetMethod::DomainGE,
            "bird" | "laplace" => DetMethod::Berkowitz,
            _ => {
                return Err(MatrixError::Value(format!(
                    "Determinant method '{s}' unrecognized"
                )))
            }
        })
    }
}

// ----------------------------------------------------------------------
// exact rational kernels
// ----------------------------------------------------------------------

impl RatMat {
    /// Exact determinant.
    pub fn det(&self) -> Rat {
        let n = self.rows;
        let (a, scales) = self.to_int_rows();
        let d = bareiss_det(a, n, &Int::ONE);
        let mut s = Int::ONE;
        for l in &scales {
            s = s.mul(l);
        }
        Rat::new(d, s)
    }

    /// Reduced row echelon form and pivot columns.
    pub fn rref(&self) -> (RatMat, Vec<usize>) {
        let (mut a, _) = self.to_int_rows();
        if self.rows.min(self.cols) >= MODULAR_MIN {
            if let Some((x, d, pivots)) = super::modular::rref_int(&a, self.rows, self.cols) {
                let r = pivots.len();
                let mut data = vec![Rat::zero(); self.rows * self.cols];
                for i in 0..r {
                    for j in 0..self.cols {
                        let v = &x[i * self.cols + j];
                        if !v.is_zero() {
                            data[i * self.cols + j] = Rat::new(v.clone(), d.clone());
                        }
                    }
                }
                return (RatMat { rows: self.rows, cols: self.cols, data }, pivots);
            }
        }
        let (den, pivots) = rref_den(&mut a, self.rows, self.cols, &Int::ONE);
        let data = a.into_iter().map(|x| if x.is_zero() { Rat::zero() } else { Rat::new(x, den.clone()) }).collect();
        (RatMat { rows: self.rows, cols: self.cols, data }, pivots)
    }

    pub fn rank(&self) -> usize {
        let (mut a, _) = self.to_int_rows();
        rref_den(&mut a, self.rows, self.cols, &Int::ONE).1.len()
    }

    /// Exact inverse (`None` when singular).
    pub fn inv(&self) -> Option<RatMat> {
        let n = self.rows;
        let (a, scales) = self.to_int_rows();
        if n >= MODULAR_MIN {
            if let Some(res) = super::modular::inv_int(&a, n) {
                let (adj, det) = res?;
                let data = adj
                    .iter()
                    .enumerate()
                    .map(|(k, x)| {
                        if x.is_zero() {
                            Rat::zero()
                        } else {
                            Rat::new(x.mul(&scales[k % n]), det.clone())
                        }
                    })
                    .collect();
                return Some(RatMat { rows: n, cols: n, data });
            }
        }
        self.inv_ff_scaled(a, scales)
    }

    /// Fraction-free inverse (always available).
    pub fn inv_ff(&self) -> Option<RatMat> {
        let (a, scales) = self.to_int_rows();
        self.inv_ff_scaled(a, scales)
    }

    fn inv_ff_scaled(&self, a: Vec<Int>, scales: Vec<Int>) -> Option<RatMat> {
        let n = self.rows;
        let w = 2 * n;
        let mut aug = Vec::with_capacity(n * w);
        for i in 0..n {
            aug.extend_from_slice(&a[i * n..(i + 1) * n]);
            for j in 0..n {
                aug.push(if i == j { Int::ONE } else { Int::ZERO });
            }
        }
        let (den, pivots) = rref_den(&mut aug, n, w, &Int::ONE);
        if pivots.len() < n || pivots[n - 1] != n - 1 {
            return None;
        }
        // A'^{-1} = X / den,  A^{-1} = A'^{-1} diag(scales)
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                let x = &aug[i * w + n + j];
                data.push(if x.is_zero() {
                    Rat::zero()
                } else {
                    Rat::new(x.mul(&scales[j]), den.clone())
                });
            }
        }
        Some(RatMat { rows: n, cols: n, data })
    }

    /// Fraction-free reduced row echelon form (always available).
    pub fn rref_ff(&self) -> (RatMat, Vec<usize>) {
        let (mut a, _) = self.to_int_rows();
        let (den, pivots) = rref_den(&mut a, self.rows, self.cols, &Int::ONE);
        let data = a.into_iter().map(|x| if x.is_zero() { Rat::zero() } else { Rat::new(x, den.clone()) }).collect();
        (RatMat { rows: self.rows, cols: self.cols, data }, pivots)
    }

    /// Characteristic polynomial coefficients `[1, c1, ..., cn]`.
    pub fn charpoly(&self) -> Vec<Rat> {
        let n = self.rows;
        // clear all denominators with one common scale d: A = B/d,
        // det(tI - B/d) = d^-n det(d t I - B) => c_k = c'_k / d^k
        let mut d = Int::ONE;
        for r in &self.data {
            d = d.lcm(&r.den);
        }
        let b: Vec<Int> = self.data.iter().map(|r| r.num.mul(&d.divexact(&r.den))).collect();
        let cp = berkowitz(&b, n, &Int::ONE);
        let mut out = Vec::with_capacity(n + 1);
        let mut dk = Int::ONE;
        for c in cp {
            out.push(Rat::new(c, dk.clone()));
            dk = dk.mul(&d);
        }
        out
    }
}

// ----------------------------------------------------------------------
// polynomial kernels
// ----------------------------------------------------------------------

fn poly_det(pm: &PolyMat, scales: &[Int]) -> Expr {
    let n = pm.rows;
    let one = MPoly::one(pm.ctx.nvars());
    let d = bareiss_det(pm.data.clone(), n, &one);
    let mut s = Int::ONE;
    for l in scales {
        s = s.mul(l);
    }
    pm.ctx.to_expr_scaled(&d, &s)
}

// ----------------------------------------------------------------------
// SymPy's pivot search
// ----------------------------------------------------------------------

/// The "simplify" used in the pivot search (enough to decide zero-ness of
/// rational functions).
fn pivot_simp(e: &Expr) -> Expr {
    let x = expand(e);
    if x.is_atom() {
        return x;
    }
    let c = cancel(&x);
    if c.is_number() {
        return c;
    }
    // Pythagorean identities (what SymPy's simplify finds for rotation
    // like matrices); only used when they decide the value is a number
    let t = expand(&trig_squares(&c));
    if t.is_number() {
        return t;
    }
    c
}

/// Light stand-in for SymPy's `simplify` on generic (EX) coefficients:
/// expansion, cancellation and the Pythagorean identities, keeping the
/// smallest result.
pub(crate) fn ex_simplify(e: &Expr) -> Expr {
    let c0 = cancel(&expand(e));
    if c0.is_atom() {
        return c0;
    }
    let mut best = c0.clone();
    let mut size = crate::traverse::node_count(&c0);
    for cand in [expand(&trig_squares(&c0)), expand(&trig_squares_rev(&c0))] {
        let n = crate::traverse::node_count(&cand);
        if n < size {
            size = n;
            best = cand;
        }
    }
    best
}

/// Rewrite `cos(u)**(2k)` as `(1 - sin(u)**2)**k`.
fn trig_squares_rev(e: &Expr) -> Expr {
    use crate::expr::Func;
    if e.is_atom() {
        return e.clone();
    }
    if let Kind::Pow(b, x) = e.kind() {
        if let (Kind::Func(Func::Cos, args), Some(n)) = (b.kind(), x.as_i64()) {
            if n >= 2 && n % 2 == 0 {
                let s = crate::functions::sin(&args[0]);
                return pow(&sub(&Expr::one(), &pow(&s, &Expr::int(2))), &Expr::int(n / 2));
            }
        }
    }
    crate::traverse::map_children(e, &mut |c| trig_squares_rev(c))
}

/// Rewrite `sin(u)**(2k)` as `(1 - cos(u)**2)**k` and `cosh(u)**(2k)` as
/// `(1 + sinh(u)**2)**k`.
fn trig_squares(e: &Expr) -> Expr {
    use crate::expr::Func;
    if e.is_atom() {
        return e.clone();
    }
    if let Kind::Pow(b, x) = e.kind() {
        if let (Kind::Func(f, args), Some(n)) = (b.kind(), x.as_i64()) {
            if n >= 2 && n % 2 == 0 {
                let two = Expr::int(2);
                let k = Expr::int(n / 2);
                match f {
                    Func::Sin => {
                        let c = crate::functions::cos(&args[0]);
                        return pow(&sub(&Expr::one(), &pow(&c, &two)), &k);
                    }
                    Func::Cosh => {
                        let s = crate::functions::sinh(&args[0]);
                        return pow(&add(&[Expr::one(), pow(&s, &two)]), &k);
                    }
                    _ => {}
                }
            }
        }
    }
    crate::traverse::map_children(e, &mut |c| trig_squares(c))
}

struct Pivot {
    offset: Option<usize>,
    value: Option<Expr>,
    newly: Vec<(usize, Expr)>,
}

fn is_float_or_int(e: &Expr) -> Option<bool> {
    match e.as_number() {
        Some(Number::Float(_)) => Some(true),
        Some(Number::Int(_)) | Some(Number::BigInt(_)) => Some(false),
        _ => None,
    }
}

/// SymPy's `_find_reasonable_pivot` with `iszerofunc`.
fn find_pivot(col: &[Expr], iszero: &dyn Fn(&Expr) -> Option<bool>, simp: bool) -> Pivot {
    // numerical column: partial pivoting
    let kinds: Vec<Option<bool>> = col.iter().map(is_float_or_int).collect();
    if !col.is_empty() && kinds.iter().all(|k| k.is_some()) && kinds.iter().any(|k| *k == Some(true)) {
        let abs: Vec<f64> = col.iter().map(|e| e.as_number().unwrap().to_f64().abs()).collect();
        let mut best = 0;
        for (i, v) in abs.iter().enumerate() {
            if *v > abs[best] {
                best = i;
            }
        }
        if abs[best] == 0.0 {
            let newly = col
                .iter()
                .enumerate()
                .filter(|(_, e)| !e.is_zero())
                .map(|(i, _)| (i, Expr::zero()))
                .collect();
            return Pivot { offset: None, value: None, newly };
        }
        return Pivot { offset: Some(best), value: Some(col[best].clone()), newly: vec![] };
    }
    let mut newly = Vec::new();
    let mut possible: Vec<Option<bool>> = Vec::with_capacity(col.len());
    for (i, x) in col.iter().enumerate() {
        let z = iszero(x);
        if z == Some(false) {
            return Pivot { offset: Some(i), value: Some(x.clone()), newly };
        }
        possible.push(z);
    }
    if possible.iter().all(|z| *z == Some(true)) {
        return Pivot { offset: None, value: None, newly };
    }
    if simp {
        for (i, x) in col.iter().enumerate() {
            if possible[i].is_some() {
                continue;
            }
            let s = pivot_simp(x);
            let z = iszero(&s);
            if z.is_some() {
                newly.push((i, s.clone()));
            }
            if z == Some(false) {
                return Pivot { offset: Some(i), value: Some(s), newly };
            }
            possible[i] = z;
        }
        if possible.iter().all(|z| *z == Some(true)) {
            return Pivot { offset: None, value: None, newly };
        }
    }
    let i = possible.iter().position(|z| z.is_none()).unwrap();
    Pivot { offset: Some(i), value: Some(col[i].clone()), newly }
}

fn iszero_assume(e: &Expr) -> Option<bool> {
    is_zero_assume(e)
}

fn iszero_expand(e: &Expr) -> Option<bool> {
    Some(is_zero_expand(e))
}

// ----------------------------------------------------------------------
// generic algorithms (faithful ports of SymPy's Matrix routines)
// ----------------------------------------------------------------------

/// SymPy's recursive Bareiss with `dotprodsimp` on every entry.
fn bareiss_expr(m: &Matrix) -> Expr {
    fn rec(mat: Matrix, cumm: Expr) -> Expr {
        let n = mat.rows();
        if n == 0 {
            return Expr::one();
        }
        if n == 1 {
            return mat.get(0, 0).clone();
        }
        let col: Vec<Expr> = (0..n).map(|i| mat.get(i, 0).clone()).collect();
        let p = find_pivot(&col, &iszero_expand, false);
        let (Some(pos), Some(pv)) = (p.offset, p.value) else {
            return Expr::zero();
        };
        let rows: Vec<usize> = (0..n).filter(|&i| i != pos).collect();
        let mut data = Vec::with_capacity((n - 1) * (n - 1));
        let inv_cumm = if cumm.is_one() { None } else { Some(pow(&cumm, &Expr::minus_one())) };
        for &i in &rows {
            for j in 0..n - 1 {
                let r = sub(
                    &mul2(&pv, mat.get(i, j + 1)),
                    &mul2(mat.get(pos, j + 1), mat.get(i, 0)),
                );
                let r = match &inv_cumm {
                    Some(ic) => mul2(&r, ic),
                    None => r,
                };
                data.push(dotprodsimp(&r));
            }
        }
        let sub_m = Matrix::new_unchecked(n - 1, n - 1, data);
        let d = rec(sub_m, pv);
        if pos % 2 == 1 {
            mul2(&Expr::minus_one(), &d)
        } else {
            d
        }
    }
    rec(m.clone(), Expr::one())
}

/// SymPy's `_berkowitz_vector` (plain arithmetic, no simplification).
fn berkowitz_vector_expr(m: &Matrix) -> Vec<Expr> {
    let n = m.rows();
    if n == 0 {
        return vec![Expr::one()];
    }
    if n == 1 {
        return vec![Expr::one(), neg(m.get(0, 0))];
    }
    let a = m.get(0, 0).clone();
    let r = m.slice(0, 1, 1, n);
    let c = m.slice(1, n, 0, 1);
    let sub_a = m.slice(1, n, 1, n);
    let mut diags = vec![c];
    for i in 0..n - 2 {
        let d = sub_a.matmul_plain(&diags[i]);
        diags.push(d);
    }
    let neg_r = r.neg();
    let mut ds: Vec<Expr> = vec![Expr::one(), neg(&a)];
    for d in &diags {
        ds.push(neg_r.matmul_plain(d).get(0, 0).clone());
    }
    let sub_v = berkowitz_vector_expr(&sub_a);
    // toeplitz (n+1) x n times sub_v (n)
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let mut terms = Vec::new();
        for j in 0..n {
            if j > i {
                continue;
            }
            let t = &ds[i - j];
            if t.is_zero() || sub_v[j].is_zero() {
                continue;
            }
            terms.push(mul2(t, &sub_v[j]));
        }
        out.push(add(&terms));
    }
    out
}

impl Matrix {
    /// Matrix product with expression arithmetic only.
    pub(crate) fn matmul_plain(&self, o: &Matrix) -> Matrix {
        let n = self.cols();
        let skip_zeros = !self.has_nonfinite() && !o.has_nonfinite();
        Matrix::from_fn(self.rows(), o.cols(), |i, j| {
            let terms: Vec<Expr> = (0..n)
                .filter(|&k| !skip_zeros || (!self.get(i, k).is_zero() && !o.get(k, j).is_zero()))
                .map(|k| if skip_zeros { mul2(self.get(i, k), o.get(k, j)) }
                    else { super::nonfinite_product(self.get(i, k), o.get(k, j)) })
                .collect();
            add(&terms)
        })
    }
}

/// SymPy's `_row_reduce_list`. Returns `(matrix, pivot columns, swaps)`.
#[allow(clippy::too_many_arguments)]
fn row_reduce(
    m: &Matrix,
    normalize_last: bool,
    normalize: bool,
    zero_above: bool,
    simp: bool,
) -> (Matrix, Vec<usize>, Vec<(usize, usize)>) {
    let (rows, cols) = m.shape();
    let mut mat: Vec<Expr> = m.data().to_vec();
    let mut pivot_cols = Vec::new();
    let mut swaps = Vec::new();
    let (mut piv_row, mut piv_col) = (0usize, 0usize);
    while piv_col < cols && piv_row < rows {
        let col: Vec<Expr> = (piv_row..rows).map(|i| mat[i * cols + piv_col].clone()).collect();
        let p = find_pivot(&col, &iszero_assume, simp);
        for (off, v) in p.newly {
            mat[(off + piv_row) * cols + piv_col] = v;
        }
        let Some(off) = p.offset else {
            piv_col += 1;
            continue;
        };
        let mut pivot_val = p.value.unwrap();
        pivot_cols.push(piv_col);
        if off != 0 {
            for k in 0..cols {
                mat.swap(piv_row * cols + k, (off + piv_row) * cols + k);
            }
            swaps.push((piv_row, off + piv_row));
        }
        if !normalize_last {
            let (i, j) = (piv_row, piv_col);
            mat[i * cols + j] = Expr::one();
            for p in i * cols + j + 1..(i + 1) * cols {
                mat[p] = dotprodsimp(&div(&mat[p], &pivot_val));
            }
            pivot_val = Expr::one();
        }
        for row in 0..rows {
            if row == piv_row {
                continue;
            }
            if !zero_above && row < piv_row {
                continue;
            }
            let val = mat[row * cols + piv_col].clone();
            if iszero_assume(&val) == Some(true) {
                continue;
            }
            for k in 0..cols {
                let a = &mat[row * cols + k];
                let b = &mat[piv_row * cols + k];
                let v = sub(&mul2(&pivot_val, a), &mul2(&val, b));
                mat[row * cols + k] = dotprodsimp(&v);
            }
        }
        piv_row += 1;
    }
    if normalize_last && normalize {
        for (pi, &pj) in pivot_cols.iter().enumerate() {
            let pv = mat[pi * cols + pj].clone();
            mat[pi * cols + pj] = Expr::one();
            for p in pi * cols + pj + 1..(pi + 1) * cols {
                mat[p] = dotprodsimp(&div(&mat[p], &pv));
            }
        }
    }
    (Matrix::new_unchecked(rows, cols, mat), pivot_cols, swaps)
}

/// SymPy's `_LUdecomposition_Simple`: combined LU matrix and row swaps.
fn lu_simple(m: &Matrix, rankcheck: bool) -> MResult<(Matrix, Vec<(usize, usize)>)> {
    let (rows, cols) = m.shape();
    if rows == 0 || cols == 0 {
        return Ok((Matrix::zeros(rows, cols), vec![]));
    }
    let rank_err = || {
        MatrixError::Value(
            "Rank of matrix is strictly less than number of rows or columns. Pass keyword argument rankcheck=False to compute the LU decomposition of this matrix.".into(),
        )
    };
    let mut lu = m.clone();
    let mut swaps = Vec::new();
    let mut pivot_col = 0usize;
    for pivot_row in 0..rows - 1 {
        let mut found: Option<(usize, Expr)> = None;
        while pivot_col != cols {
            // naive pivot search
            let mut indet: Option<(usize, Expr)> = None;
            for r in pivot_row..rows {
                let v = lu.get(r, pivot_col);
                match is_zero_assume(v) {
                    Some(false) => {
                        found = Some((r - pivot_row, v.clone()));
                        break;
                    }
                    None if indet.is_none() => indet = Some((r - pivot_row, v.clone())),
                    _ => {}
                }
            }
            if found.is_none() {
                found = indet;
            }
            if found.is_some() {
                break;
            }
            pivot_col += 1;
        }
        if rankcheck && pivot_col != pivot_row {
            return Err(rank_err());
        }
        let Some((off, _)) = found else {
            return Ok((lu, swaps));
        };
        let cand = pivot_row + off;
        if pivot_row != cand {
            swaps.push((pivot_row, cand));
            for k in 0..pivot_row {
                let a = lu.get(pivot_row, k).clone();
                let b = lu.get(cand, k).clone();
                lu.set(pivot_row, k, b);
                lu.set(cand, k, a);
            }
            for k in pivot_col..cols {
                let a = lu.get(pivot_row, k).clone();
                let b = lu.get(cand, k).clone();
                lu.set(pivot_row, k, b);
                lu.set(cand, k, a);
            }
        }
        let start_col = pivot_col + 1;
        for row in pivot_row + 1..rows {
            let l = div(lu.get(row, pivot_col), lu.get(pivot_row, pivot_col));
            lu.set(row, pivot_row, l.clone());
            for c in start_col..cols {
                let v = sub(lu.get(row, c), &mul2(&l, lu.get(pivot_row, c)));
                lu.set(row, c, v);
            }
        }
        if pivot_row != pivot_col {
            for row in pivot_row + 1..rows {
                lu.set(row, pivot_col, Expr::zero());
            }
        }
        pivot_col += 1;
        if pivot_col == cols {
            return Ok((lu, swaps));
        }
    }
    if rankcheck {
        let k = rows.min(cols) - 1;
        if is_zero_assume(lu.get(k, k)) == Some(true) {
            return Err(rank_err());
        }
    }
    Ok((lu, swaps))
}

fn swaps_parity(s: &[(usize, usize)]) -> bool {
    s.len() % 2 == 1
}

/// Strongly connected components (SymPy's Tarjan order).
pub fn scc(m: &Matrix) -> Vec<Vec<usize>> {
    let n = m.rows();
    let gmap: Vec<Vec<usize>> = (0..n)
        .map(|i| (0..n).filter(|&j| !m.get(i, j).is_zero()).collect())
        .collect();
    let mut lowlink = vec![usize::MAX; n];
    let mut index = vec![usize::MAX; n];
    let mut on_stack = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut comps = Vec::new();
    let mut callstack: Vec<(usize, usize)> = Vec::new();
    let mut counter = 0usize;
    for v in 0..n {
        if index[v] != usize::MAX {
            continue;
        }
        // start(v)
        index[v] = stack.len();
        lowlink[v] = stack.len();
        let _ = counter;
        counter += 1;
        stack.push(v);
        on_stack[v] = true;
        callstack.push((v, 0));
        while let Some(&(v1, it)) = callstack.last() {
            if it < gmap[v1].len() {
                callstack.last_mut().unwrap().1 += 1;
                let v2 = gmap[v1][it];
                if index[v2] == usize::MAX {
                    index[v2] = stack.len();
                    lowlink[v2] = stack.len();
                    stack.push(v2);
                    on_stack[v2] = true;
                    callstack.push((v2, 0));
                } else if on_stack[v2] {
                    lowlink[v1] = lowlink[v1].min(index[v2]);
                }
            } else {
                // finish(v1)
                if lowlink[v1] == index[v1] {
                    let mut comp = vec![];
                    loop {
                        let w = stack.pop().unwrap();
                        on_stack[w] = false;
                        comp.push(w);
                        if w == v1 {
                            break;
                        }
                    }
                    comp.reverse();
                    comps.push(comp);
                }
                callstack.pop();
                if let Some(&(p, _)) = callstack.last() {
                    lowlink[p] = lowlink[p].min(lowlink[v1]);
                }
            }
        }
    }
    comps
}

/// Determinant of a 2x2 / 3x3 matrix by the explicit formula.
fn det_small(m: &Matrix) -> Expr {
    let g = |i: usize, j: usize| m.get(i, j);
    let e = if m.rows() == 2 {
        sub(&mul2(g(0, 0), g(1, 1)), &mul2(g(0, 1), g(1, 0)))
    } else {
        let t = |a: (usize, usize), b: (usize, usize), c: (usize, usize)| {
            mul(&[g(a.0, a.1).clone(), g(b.0, b.1).clone(), g(c.0, c.1).clone()])
        };
        add(&[
            t((0, 0), (1, 1), (2, 2)),
            t((0, 1), (1, 2), (2, 0)),
            t((0, 2), (1, 0), (2, 1)),
            neg(&t((0, 2), (1, 1), (2, 0))),
            neg(&t((0, 0), (1, 2), (2, 1))),
            neg(&t((0, 1), (1, 0), (2, 2))),
        ])
    };
    dotprodsimp(&e)
}

/// `uniquely_named_symbol(name, [M])` with the `_` prefix modifier.
pub(crate) fn unique_symbol(name: &str, ms: &[&Matrix]) -> Expr {
    let mut names: Vec<String> = Vec::new();
    for m in ms {
        for e in m.data() {
            crate::traverse::preorder(e, &mut |s| match s.kind() {
                Kind::Sym(n) => names.push(n.to_string()),
                Kind::Func(crate::expr::Func::Undefined(n), _) => names.push(n.to_string()),
                _ => {}
            });
        }
    }
    let mut x = name.to_string();
    while names.iter().any(|n| *n == x) {
        x = format!("_{x}");
    }
    Expr::symbol(&x)
}

/// Characteristic polynomial: coefficients (highest degree first) and the
/// SymPy domain name of the `PurePoly`.
#[derive(Clone, Debug)]
pub struct CharPoly {
    pub var: Expr,
    pub coeffs: Vec<Expr>,
    pub domain: String,
}

impl CharPoly {
    /// `PurePoly.as_expr()`: `sum(c_k * var**k)` (coefficients kept as
    /// they are, like SymPy's conversion from the ground domain).
    pub fn as_expr(&self) -> Expr {
        let n = self.coeffs.len() - 1;
        let terms: Vec<Expr> = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| mul2(c, &pow(&self.var, &Expr::int((n - i) as i64))))
            .collect();
        add(&terms)
    }
}

impl Matrix {
    // ------------------------------------------------------------------
    // determinant
    // ------------------------------------------------------------------

    /// Determinant (Bareiss, SymPy's default).
    pub fn det(&self) -> MResult<Expr> {
        self.det_method(DetMethod::Bareiss)
    }

    pub fn det_method(&self, method: DetMethod) -> MResult<Expr> {
        if !self.is_square() {
            return Err(non_square());
        }
        let n = self.rows();
        match n {
            0 => return Ok(Expr::one()),
            1 => return Ok(self.get(0, 0).clone()),
            2 | 3 => return Ok(det_small(self)),
            _ => {}
        }
        if method == DetMethod::DomainGE || RatMat::from_matrix(self).is_some() {
            if let Some(r) = RatMat::from_matrix(self) {
                // the value is unique; SCC splitting is irrelevant
                return Ok(Expr::num(r.det().to_number()));
            }
        }
        let comps = scc(self);
        let dom = if method == DetMethod::Bareiss { Some(matrix_domain(self)) } else { None };
        let mut dets = Vec::with_capacity(comps.len());
        for b in &comps {
            let blk = self.extract(b, b)?;
            let d = match method {
                DetMethod::Bareiss | DetMethod::DomainGE => {
                    let mut r = None;
                    if let Some(Domain::Poly { ctx, .. }) = &dom {
                        if blk.rows() > 3 {
                            if let Some((pm, scales)) = PolyMat::from_matrix(&blk, ctx) {
                                r = Some(poly_det(&pm, &scales));
                            }
                        }
                    }
                    match r {
                        Some(r) => r,
                        None => bareiss_expr(&blk),
                    }
                }
                DetMethod::Berkowitz => {
                    let v = berkowitz_vector_expr(&blk);
                    let last = v.last().unwrap().clone();
                    if (v.len() - 1) % 2 == 1 {
                        mul2(&Expr::minus_one(), &last)
                    } else {
                        last
                    }
                }
                DetMethod::LU => {
                    let (lu, swaps) = lu_simple(&blk, false)?;
                    let k = lu.rows() - 1;
                    if is_zero_assume(lu.get(k, k)) == Some(true) {
                        Expr::zero()
                    } else {
                        let mut d = if swaps_parity(&swaps) { Expr::minus_one() } else { Expr::one() };
                        for i in 0..lu.rows() {
                            d = mul2(&d, lu.get(i, i));
                        }
                        d
                    }
                }
            };
            dets.push(d);
        }
        Ok(mul(&dets))
    }

    // ------------------------------------------------------------------
    // inverse
    // ------------------------------------------------------------------

    /// Inverse (SymPy's default method selection).
    pub fn inv(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare("A Matrix must be square to invert.".into()));
        }
        let n = self.rows();
        if n == 0 {
            return Ok(self.clone());
        }
        if let Some(r) = RatMat::from_matrix(self) {
            return r.inv().map(|x| x.to_matrix()).ok_or_else(non_invertible);
        }
        match matrix_domain(self) {
            Domain::Poly { ctx, .. } => {
                if let Some((pm, scales)) = PolyMat::from_matrix(self, &ctx) {
                    let scales: Vec<MPoly> =
                        scales.iter().map(|s| MPoly::constant(ctx.nvars(), s.clone())).collect();
                    return poly_inverse(&pm, &scales);
                }
                self.inverse_ge()
            }
            Domain::Frac { ctx, .. } => {
                if let Some((pm, scales)) = PolyMat::from_frac_matrix(self, &ctx) {
                    return poly_inverse(&pm, &scales);
                }
                self.inverse_ge()
            }
            _ => self.inverse_ge(),
        }
    }

    /// Inverse by Gauss-Jordan elimination of `[M | I]` (SymPy's `GE`).
    pub fn inverse_ge(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare("A Matrix must be square to invert.".into()));
        }
        let n = self.rows();
        let big = Matrix::hstack(&[self, &Matrix::eye(n)])?;
        let (red, _) = big.rref_generic(true);
        if (0..n).any(|j| is_zero_assume(red.get(j, j)) == Some(true)) {
            return Err(non_invertible());
        }
        Ok(red.slice(0, n, n, 2 * n))
    }

    /// Inverse via the LU decomposition (SymPy's `LU`).
    pub fn inverse_lu(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare("A Matrix must be square to invert.".into()));
        }
        self.lu_solve(&Matrix::eye(self.rows()))
    }

    /// Inverse via the adjugate (SymPy's `ADJ`).
    pub fn inverse_adj(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare("A Matrix must be square to invert.".into()));
        }
        let d = self.det_method(DetMethod::Berkowitz)?;
        if is_zero_assume(&d) == Some(true) {
            return Err(non_invertible());
        }
        Ok(self.adjugate()?.applyfunc(|e| div(e, &d)))
    }

    /// Adjugate (classical adjoint), from Berkowitz minors.
    pub fn adjugate(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(non_square());
        }
        let n = self.rows();
        if let Some(r) = RatMat::from_matrix(self) {
            if n > 0 {
                // adj(A) = det(A) * A^{-1} when invertible; otherwise minors
                if let Some(inv) = r.inv() {
                    let d = r.det();
                    let data = inv.data.iter().map(|x| Expr::num(x.mul(&d).to_number())).collect();
                    return Ok(Matrix::new_unchecked(n, n, data));
                }
            }
        }
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                // adj[i][j] = cofactor(j, i)
                data.push(self.cofactor(j, i)?);
            }
        }
        Ok(Matrix::new_unchecked(n, n, data))
    }

    pub fn minor_submatrix(&self, i: usize, j: usize) -> Matrix {
        let rows: Vec<usize> = (0..self.rows()).filter(|&r| r != i).collect();
        let cols: Vec<usize> = (0..self.cols()).filter(|&c| c != j).collect();
        self.extract(&rows, &cols).unwrap()
    }

    pub fn minor(&self, i: usize, j: usize) -> MResult<Expr> {
        if !self.is_square() {
            return Err(non_square());
        }
        self.minor_submatrix(i, j).det_method(DetMethod::Berkowitz)
    }

    pub fn cofactor(&self, i: usize, j: usize) -> MResult<Expr> {
        let m = self.minor(i, j)?;
        Ok(if (i + j) % 2 == 1 { mul2(&Expr::minus_one(), &m) } else { m })
    }

    pub fn cofactor_matrix(&self) -> MResult<Matrix> {
        if !self.is_square() {
            return Err(non_square());
        }
        let n = self.rows();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                data.push(self.cofactor(i, j)?);
            }
        }
        Ok(Matrix::new_unchecked(n, n, data))
    }

    // ------------------------------------------------------------------
    // row reduction
    // ------------------------------------------------------------------

    /// Reduced row echelon form and the pivot columns.
    pub fn rref(&self) -> (Matrix, Vec<usize>) {
        if let Some(r) = RatMat::from_matrix(self) {
            let (m, p) = r.rref();
            return (m.to_matrix(), p);
        }
        self.rref_generic(true)
    }

    fn rref_generic(&self, simp: bool) -> (Matrix, Vec<usize>) {
        if let Some(r) = RatMat::from_matrix(self) {
            let (m, p) = r.rref();
            return (m.to_matrix(), p);
        }
        let (m, p, _) = row_reduce(self, true, true, true, simp);
        (m, p)
    }

    /// Echelon form (not reduced, not normalised) and pivots.
    pub fn echelon_form(&self) -> (Matrix, Vec<usize>) {
        let (m, p, _) = row_reduce(self, true, false, false, true);
        (m, p)
    }

    pub fn rank(&self) -> usize {
        let (rows, cols) = self.shape();
        if rows == 0 || cols == 0 {
            return 0;
        }
        if let Some(r) = RatMat::from_matrix(self) {
            return r.rank();
        }
        let zeros: Vec<Option<bool>> = self.data().iter().map(is_zero_assume).collect();
        if (rows <= 1 || cols <= 1) && zeros.contains(&Some(false)) {
            return 1;
        }
        if rows == 2 && cols == 2 {
            if !zeros.contains(&Some(false)) && !zeros.contains(&None) {
                return 0;
            }
            if let Ok(d) = self.det() {
                let dz = is_zero_assume(&d);
                if dz == Some(true) && zeros.contains(&Some(false)) {
                    return 1;
                }
                if dz == Some(false) {
                    return 2;
                }
            }
        }
        // permute complicated columns to the right
        let mut cx: Vec<(usize, usize)> = (0..cols)
            .map(|j| ((0..rows).filter(|&i| is_zero_assume(self.get(i, j)).is_none()).count(), j))
            .collect();
        cx.sort();
        let perm: Vec<usize> = cx.into_iter().map(|(_, j)| j).collect();
        let all_rows: Vec<usize> = (0..rows).collect();
        let pm = self.extract(&all_rows, &perm).unwrap();
        let (_, pivots, _) = row_reduce(&pm, true, false, false, true);
        pivots.len()
    }

    /// Basis of the null space (column vectors).
    pub fn nullspace(&self) -> Vec<Matrix> {
        let (reduced, pivots) = self.rref();
        let cols = self.cols();
        let mut basis = Vec::new();
        for free in (0..cols).filter(|c| !pivots.contains(c)) {
            let mut v = vec![Expr::zero(); cols];
            v[free] = Expr::one();
            for (pr, &pc) in pivots.iter().enumerate() {
                v[pc] = sub(&v[pc], reduced.get(pr, free));
            }
            basis.push(Matrix::column(v));
        }
        basis
    }

    /// Basis of the column space (columns of `self` at the pivots).
    pub fn columnspace(&self) -> Vec<Matrix> {
        let pivots = if let Some(r) = RatMat::from_matrix(self) { r.rref().1 } else { self.echelon_form().1 };
        pivots.iter().map(|&j| self.col(j)).collect()
    }

    /// Basis of the row space (non-zero rows of the echelon form).
    pub fn rowspace(&self) -> Vec<Matrix> {
        let (red, pivots) = self.echelon_form();
        (0..pivots.len()).map(|i| red.row(i)).collect()
    }

    // ------------------------------------------------------------------
    // LU and linear systems
    // ------------------------------------------------------------------

    /// `(L, U, row swaps)` with `P*M = L*U`.
    pub fn lu_decomposition(&self) -> MResult<(Matrix, Matrix, Vec<(usize, usize)>)> {
        self.lu_decomposition_with_rankcheck(false)
    }

    pub fn lu_decomposition_with_rankcheck(&self, rankcheck: bool) -> MResult<(Matrix, Matrix, Vec<(usize, usize)>)> {
        let (combined, p) = lu_simple(self, rankcheck)?;
        let (r, c) = combined.shape();
        let l = Matrix::from_fn(r, r, |i, j| {
            if i < j {
                Expr::zero()
            } else if i == j {
                Expr::one()
            } else if j < c {
                combined.get(i, j).clone()
            } else {
                Expr::zero()
            }
        });
        let u = Matrix::from_fn(r, c, |i, j| if i > j { Expr::zero() } else { combined.get(i, j).clone() });
        Ok((l, u, p))
    }

    /// Solve `self * x = rhs` by LU decomposition (SymPy's `LUsolve`).
    pub fn lu_solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        if rhs.rows() != self.rows() {
            return Err(MatrixError::Shape(
                "``M`` and ``rhs`` must have the same number of rows.".into(),
            ));
        }
        let (m, n) = self.shape();
        if m < n {
            return Err(MatrixError::NotImplemented("Underdetermined systems not supported.".into()));
        }
        let (a, perm) = lu_simple(self, true).map_err(|_| non_invertible())?;
        let mut b = rhs.permute_rows_swaps(&perm);
        let bc = b.cols();
        for i in 0..m {
            for j in 0..i.min(n) {
                let scale = a.get(i, j).clone();
                for k in 0..bc {
                    let v = sub(b.get(i, k), &mul2(&scale, b.get(j, k)));
                    b.set(i, k, v);
                }
            }
        }
        if m > n {
            for i in n..m {
                for j in 0..bc {
                    if is_zero_assume(b.get(i, j)) != Some(true) {
                        return Err(MatrixError::Value("The system is inconsistent.".into()));
                    }
                }
            }
            b = b.slice(0, n, 0, bc);
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                let scale = a.get(i, j).clone();
                for k in 0..bc {
                    let v = sub(b.get(i, k), &mul2(&scale, b.get(j, k)));
                    b.set(i, k, v);
                }
            }
            let inv = pow(a.get(i, i), &Expr::minus_one());
            for k in 0..bc {
                let v = mul2(&inv, b.get(i, k));
                b.set(i, k, v);
            }
        }
        Ok(b)
    }

    /// Gauss-Jordan solution `(solution, parameters)`; the solution may
    /// contain the free parameters `tau0, tau1, ...`.
    pub fn gauss_jordan_solve(&self, b: &Matrix) -> MResult<(Matrix, Matrix)> {
        let aug = Matrix::hstack(&[self, b])?;
        let bcols = b.cols();
        let (row, col) = (self.rows(), self.cols());
        let _ = row;
        let (a_full, pivots) = aug.rref_generic(true);
        let a = a_full.slice(0, a_full.rows(), 0, col);
        let v = a_full.slice(0, a_full.rows(), col, col + bcols);
        let pivots: Vec<usize> = pivots.into_iter().filter(|&p| p < col).collect();
        let rank = pivots.len();
        let free: Vec<usize> = (0..col).filter(|c| !pivots.contains(c)).collect();
        let tail = v.slice(rank, v.rows(), 0, bcols);
        if !tail.data().iter().all(|e| is_zero_assume(e) == Some(true)) {
            return Err(MatrixError::Value("Linear system has no solution".into()));
        }
        // free parameters
        let base = {
            let mut names: Vec<String> = Vec::new();
            for e in aug.data() {
                for s in crate::traverse::free_symbols(e) {
                    names.push(s.to_string().trim_end_matches(|c: char| c.is_ascii_digit()).to_string());
                }
            }
            let mut x = "tau".to_string();
            while names.contains(&x) {
                x = format!("_{x}");
            }
            x
        };
        let ntau = (col - rank) * bcols;
        let taus: Vec<Expr> = (0..ntau).map(|k| Expr::symbol(&format!("{base}{k}"))).collect();
        let tau = Matrix::new_unchecked(col - rank, bcols, taus);
        let vv = a.extract(&(0..rank).collect::<Vec<_>>(), &free)?;
        let vt = v.slice(0, rank, 0, bcols);
        let top = vt.sub(&vv.matmul_plain(&tau))?;
        let free_sol = Matrix::vstack(&[&top, &tau])?;
        let permutation: Vec<usize> = pivots.iter().chain(free.iter()).cloned().collect();
        let mut sol = Matrix::zeros(col, bcols);
        for k in 0..col {
            for j in 0..bcols {
                sol.set(permutation[k], j, free_sol.get(k, j).clone());
            }
        }
        Ok((sol, tau))
    }

    /// Solve `self * x = rhs` for a unique solution (SymPy's `solve`).
    pub fn solve(&self, rhs: &Matrix) -> MResult<Matrix> {
        match self.gauss_jordan_solve(rhs) {
            Ok((sol, tau)) => {
                if !tau.is_empty() {
                    return Err(MatrixError::NonInvertible(
                        "Matrix det == 0; not invertible. Try ``M.gauss_jordan_solve(rhs)`` to obtain a parametric solution.".into(),
                    ));
                }
                Ok(sol)
            }
            Err(MatrixError::Value(_)) => Err(non_invertible()),
            Err(e) => Err(e),
        }
    }

    // ------------------------------------------------------------------
    // characteristic polynomial
    // ------------------------------------------------------------------

    /// Characteristic polynomial `det(lambda*I - M)` in the symbol named
    /// `name` (renamed with a `_` prefix if it clashes).
    pub fn charpoly(&self, name: &str) -> MResult<CharPoly> {
        if !self.is_square() {
            return Err(non_square());
        }
        let var = unique_symbol(name, &[self]);
        let n = self.rows();
        if let Some(r) = RatMat::from_matrix(self) {
            let cs = r.charpoly();
            let qq = cs.iter().any(|c| !c.is_integer()) || !r.is_integer();
            return Ok(CharPoly {
                var,
                coeffs: cs.iter().map(|c| Expr::num(c.to_number())).collect(),
                domain: if qq { "QQ".into() } else { "ZZ".into() },
            });
        }
        let dom = matrix_domain(self);
        match &dom {
            Domain::Poly { ctx, .. } => {
                if let Some(cs) = poly_charpoly(self, ctx) {
                    return Ok(CharPoly { var, coeffs: cs, domain: dom.name() });
                }
            }
            Domain::Frac { ctx, .. } => {
                if let Some(cs) = frac_charpoly(self, ctx) {
                    return Ok(CharPoly { var, coeffs: cs, domain: dom.name() });
                }
            }
            _ => {}
        }
        // generic: Berkowitz over expressions, simplified coefficients
        let v = berkowitz_vector_expr(self);
        let coeffs: Vec<Expr> = v.iter().map(ex_simplify).collect();
        let _ = n;
        let domain = super::domain::construct_domain(&coeffs).name();
        Ok(CharPoly { var, coeffs, domain })
    }
}

fn poly_charpoly(m: &Matrix, ctx: &PolyCtx) -> Option<Vec<Expr>> {
    let n = m.rows();
    // common denominator d for all entries: A = B/d
    let mut qs = Vec::with_capacity(m.len());
    let mut d = Int::ONE;
    for e in m.data() {
        let q = ctx.to_qpoly(e)?;
        d = d.lcm(&q.den);
        qs.push(q);
    }
    let b: Vec<MPoly> = qs.iter().map(|q| q.p.mul_int(&d.divexact(&q.den))).collect();
    let one = MPoly::one(ctx.nvars());
    let cp = berkowitz(&b, n, &one);
    let mut dk = Int::ONE;
    let mut out = Vec::with_capacity(n + 1);
    for c in cp {
        out.push(ctx.to_expr_scaled(&c, &dk));
        dk = dk.mul(&d);
    }
    Some(out)
}

fn frac_charpoly(m: &Matrix, ctx: &PolyCtx) -> Option<Vec<Expr>> {
    let n = m.rows();
    let nv = ctx.nvars();
    // A = B / D with a common polynomial denominator D
    let mut nums = Vec::with_capacity(m.len());
    let mut dens = Vec::with_capacity(m.len());
    for e in m.data() {
        let (a, b) = super::ratfunc::as_numer_denom(e);
        let pa = ctx.to_qpoly(&a)?;
        let pb = ctx.to_qpoly(&b)?;
        nums.push(pa.p.mul_int(&pb.den));
        dens.push(pb.p.mul_int(&pa.den));
    }
    let mut dd = MPoly::one(nv);
    for den in &dens {
        if dd.divexact(den).is_none() {
            let g = super::mpoly::gcd(&dd, den);
            dd = dd.mul(&den.divexact(&g).unwrap());
        }
    }
    let b: Vec<MPoly> = nums.iter().zip(&dens).map(|(x, y)| x.mul(&dd.divexact(y).unwrap())).collect();
    let one = MPoly::one(nv);
    let cp = berkowitz(&b, n, &one);
    let mut dk = MPoly::one(nv);
    let mut out = Vec::with_capacity(n + 1);
    for c in cp {
        let (p, q) = cancel_polys(&c, &dk);
        let pe = ctx.to_expr(&p);
        out.push(if q.is_one() { pe } else { div(&pe, &ctx.to_expr(&q)) });
        dk = dk.mul(&dd);
    }
    Some(out)
}

/// Inverse of a polynomial matrix `A = diag(scales)^{-1} * pm` with every
/// entry reduced to lowest terms.
fn poly_inverse(pm: &PolyMat, scales: &[MPoly]) -> MResult<Matrix> {
    let n = pm.rows;
    let nv = pm.ctx.nvars();
    let one = MPoly::one(nv);
    let zero = MPoly::zero(nv);
    // (adjugate, determinant): cofactors for small matrices, fraction-free
    // Gauss-Jordan otherwise
    let (adj, den): (Vec<MPoly>, MPoly) = if n <= 4 {
        let det = bareiss_det(pm.data.clone(), n, &one);
        if det.is_zero() {
            return Err(non_invertible());
        }
        let mut adj = vec![zero.clone(); n * n];
        for i in 0..n {
            for j in 0..n {
                // adj[j][i] = (-1)^(i+j) * minor(i, j)
                let mut sub = Vec::with_capacity((n - 1) * (n - 1));
                for r in (0..n).filter(|&r| r != i) {
                    for c in (0..n).filter(|&c| c != j) {
                        sub.push(pm.get(r, c).clone());
                    }
                }
                let m = bareiss_det(sub, n - 1, &one);
                adj[j * n + i] = if (i + j) % 2 == 1 { m.neg() } else { m };
            }
        }
        (adj, det)
    } else {
        let w = 2 * n;
        let mut aug = Vec::with_capacity(n * w);
        for i in 0..n {
            for j in 0..n {
                aug.push(pm.get(i, j).clone());
            }
            for j in 0..n {
                aug.push(if i == j { one.clone() } else { zero.clone() });
            }
        }
        let (den, pivots) = rref_den(&mut aug, n, w, &one);
        if pivots.len() < n || pivots[n - 1] != n - 1 {
            return Err(non_invertible());
        }
        let mut adj = Vec::with_capacity(n * n);
        for i in 0..n {
            adj.extend_from_slice(&aug[i * w + n..(i + 1) * w]);
        }
        (adj, den)
    };
    let mut data = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let x = &adj[i * n + j];
            if x.is_zero() {
                data.push(Expr::zero());
                continue;
            }
            let num = x.mul(&scales[j]);
            let (p, q) = cancel_polys(&num, &den);
            let pe = pm.ctx.to_expr(&p);
            data.push(if q.is_one() { pe } else { div(&pe, &pm.ctx.to_expr(&q)) });
        }
    }
    Ok(Matrix::new_unchecked(n, n, data))
}

#[allow(dead_code)]
fn ring_check<R: Ring>(_: &R) {}
