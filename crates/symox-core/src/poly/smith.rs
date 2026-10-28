//! Exact Smith normal decomposition over a principal ideal domain, in the
//! reference operation order of
//! `sympy/polys/matrices/normalforms.py::_smith_normal_decomp` (1.14).
//!
//! Every row/column operation, pivot choice, unit normalisation, recursive
//! block product and divisibility repair is the reference one, so the
//! invariants and the transformation matrices `s`, `t` (with `s*m*t` the
//! Smith form) are identical, not merely equivalent. The ring operations are
//! those of [`Pid`]: typed integers for ZZ, and the native element arithmetic
//! (`DomainCoeff`) for the other admitted domains.
use super::Int;

/// The domain operations read by the reference body: element `+`, `*`,
/// unary `-`, `Domain.div`, `Domain.gcdex`, `Domain.gcd`, `Domain.exquo` and
/// the pivot normalisation (`canonical_unit`, or `1/m[0][0]` over a field).
pub trait Pid: Clone {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn add(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    /// `self += a*b` (the matrix product accumulation).
    fn add_mul(&mut self, a: &Self, b: &Self) { *self = self.add(&a.mul(b)); }
    /// `Domain.div(self, o)`: quotient and remainder.
    fn div(&self, o: &Self) -> (Self, Self);
    /// `Domain.gcdex(self, o)` as `(x, y, g)`.
    fn gcdex(&self, o: &Self) -> (Self, Self, Self);
    /// `Domain.gcd(self, o)`.
    fn gcd(&self, o: &Self) -> Self;
    /// `Domain.exquo(self, o)`.
    fn exquo(&self, o: &Self) -> Self;
    /// The factor `c` the reference multiplies a non-zero pivot by, or None
    /// when `c == domain.one` (no multiplication).
    fn unit_factor(&self) -> Option<Self>;
    /// The reference pivot search tests the row/column *index* against the
    /// domain zero (`ind[0] != zero`): over GF(p) an index that is a multiple
    /// of p compares equal to zero, the row swap is skipped and the column
    /// search runs instead.
    fn index_is_zero(i: usize) -> bool { i == 0 }
    /// True when the reference test `m[0][0] != 0` (a Python int) holds for
    /// a zero pivot and the field normalisation `1 / m[0][0]` then raises:
    /// the kernel declines so the Python body raises the reference error.
    fn zero_pivot_raises() -> bool { false }
}

pub type Matrix<C = Int> = Vec<Vec<C>>;

/// `(invariants, Some((s, t)))` when `full`, `(invariants, None)` otherwise.
pub type Decomposition<C = Int> = (Vec<C>, Option<(Matrix<C>, Matrix<C>)>);

/// The reference `ZZ.gcdex` (python ground types, `sympy.external.ntheory.gcdext`
/// reordered as `(x, y, g)`), computed on absolute values with sign correction.
pub fn gcdex(a: &Int, b: &Int) -> (Int, Int, Int) {
    if a.is_zero() || b.is_zero() {
        let g = if a.is_zero() { b.abs() } else { a.abs() };
        if g.is_zero() {
            return (Int::zero(), Int::zero(), Int::zero());
        }
        return (a.div_mod_floor(&g).0, b.div_mod_floor(&g).0, g);
    }
    let x_sign = a.is_negative();
    let y_sign = b.is_negative();
    let (mut a, mut b) = (a.abs(), b.abs());
    let (mut x, mut r) = (Int::one(), Int::zero());
    let (mut y, mut s) = (Int::zero(), Int::one());
    while !b.is_zero() {
        let (q, c) = a.div_mod_floor(&b);
        a = std::mem::replace(&mut b, c);
        let next = x.sub(&q.mul(&r));
        x = std::mem::replace(&mut r, next);
        let next = y.sub(&q.mul(&s));
        y = std::mem::replace(&mut s, next);
    }
    (if x_sign { x.neg() } else { x }, if y_sign { y.neg() } else { y }, a)
}

impl Pid for Int {
    fn zero() -> Self { Int::zero() }
    fn one() -> Self { Int::one() }
    fn is_zero(&self) -> bool { Int::is_zero(self) }
    fn add(&self, o: &Self) -> Self { Int::add(self, o) }
    fn mul(&self, o: &Self) -> Self { Int::mul(self, o) }
    fn neg(&self) -> Self { Int::neg(self) }
    fn add_mul(&mut self, a: &Self, b: &Self) { Int::add_mul(self, a, b) }
    fn div(&self, o: &Self) -> (Self, Self) { self.div_mod_floor(o) }
    fn gcdex(&self, o: &Self) -> (Self, Self, Self) { gcdex(self, o) }
    fn gcd(&self, o: &Self) -> Self { gcdex(self, o).2 }
    fn exquo(&self, o: &Self) -> Self { self.div_exact_unchecked(o) }
    // ZZ.canonical_unit: -1 for a negative pivot.
    fn unit_factor(&self) -> Option<Self> { self.is_negative().then(|| Int::from(-1i64)) }
}

fn eye<C: Pid>(n: usize) -> Matrix<C> {
    (0..n).map(|i| (0..n).map(|j| if i == j { C::one() } else { C::zero() }).collect()).collect()
}

// replace m[i, :] by a*m[i, :] + b*m[j, :] and m[j, :] by c*m[i, :] + d*m[j, :]
fn add_rows<C: Pid>(m: &mut Matrix<C>, i: usize, j: usize, a: &C, b: &C, c: &C, d: &C) {
    for k in 0..m[0].len() {
        let e = m[i][k].clone();
        m[i][k] = a.mul(&e).add(&b.mul(&m[j][k]));
        m[j][k] = c.mul(&e).add(&d.mul(&m[j][k]));
    }
}

// replace m[:, i] by a*m[:, i] + b*m[:, j] and m[:, j] by c*m[:, i] + d*m[:, j]
fn add_columns<C: Pid>(m: &mut Matrix<C>, i: usize, j: usize, a: &C, b: &C, c: &C, d: &C) {
    for row in m.iter_mut() {
        let e = row[i].clone();
        row[i] = a.mul(&e).add(&b.mul(&row[j]));
        row[j] = c.mul(&e).add(&d.mul(&row[j]));
    }
}

fn product<C: Pid>(a: &Matrix<C>, b: &Matrix<C>) -> Matrix<C> {
    let inner = b.len();
    let cols = if inner == 0 { 0 } else { b[0].len() };
    a.iter().map(|row| (0..cols).map(|j| {
        let mut sum = C::zero();
        for (l, value) in row.iter().enumerate().take(inner) { sum.add_mul(value, &b[l][j]); }
        sum
    }).collect()).collect()
}

/// Block diagonal `diag(1, small)` of order `n`.
fn bordered<C: Pid>(small: Matrix<C>, n: usize) -> Matrix<C> {
    let mut out = Vec::with_capacity(n);
    out.push((0..n).map(|j| if j == 0 { C::one() } else { C::zero() }).collect());
    for row in small {
        let mut line = Vec::with_capacity(n);
        line.push(C::zero());
        line.extend(row);
        out.push(line);
    }
    out
}

/// Smith decomposition of the `rows x cols` matrix `m` (consumed: the
/// reference works on a private copy). `None` declines a malformed shape or
/// an unreachable zero divisor, before any result is produced.
pub fn decompose<C: Pid, E>(m: Matrix<C>, rows: usize, cols: usize, full: bool,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Option<Decomposition<C>>, E> {
    check()?;
    if m.len() != rows || m.iter().any(|row| row.len() != cols) {
        return Ok(None);
    }
    recurse(m, rows, cols, full, check)
}

fn recurse<C: Pid, E>(mut m: Matrix<C>, rows: usize, cols: usize, full: bool,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Option<Decomposition<C>>, E> {
    check()?;
    if rows == 0 || cols == 0 {
        return Ok(Some((Vec::new(), full.then(|| (eye(rows), eye(cols))))));
    }
    let (mut s, mut t) = if full { (eye(rows), eye(cols)) } else { (Vec::new(), Vec::new()) };
    let zero = C::zero();
    let one = C::one();

    // Permute the rows and columns until m[0,0] is non-zero if possible.
    match (0..rows).find(|&i| !m[i][0].is_zero()) {
        Some(i) if !C::index_is_zero(i) => {
            m.swap(0, i);
            if full { s.swap(0, i); }
        }
        _ => {
            if let Some(j) = (0..cols).find(|&j| !m[0][j].is_zero()) {
                if !C::index_is_zero(j) {
                    for row in m.iter_mut() { row.swap(0, j); }
                    if full { for row in t.iter_mut() { row.swap(0, j); } }
                }
            }
        }
    }

    // Make the first row and column except m[0,0] zero.
    while (1..cols).any(|i| !m[0][i].is_zero()) || (1..rows).any(|i| !m[i][0].is_zero()) {
        check()?;
        let mut pivot = m[0][0].clone();
        for j in 1..rows {
            if m[j][0].is_zero() { continue; }
            if pivot.is_zero() { return Ok(None); }
            let (d, r) = m[j][0].div(&pivot);
            if r.is_zero() {
                let d = d.neg();
                add_rows(&mut m, 0, j, &one, &zero, &d, &one);
                if full { add_rows(&mut s, 0, j, &one, &zero, &d, &one); }
            } else {
                let (a, b, g) = pivot.gcdex(&m[j][0]);
                let d_0 = m[j][0].exquo(&g);
                let d_j = pivot.exquo(&g).neg();
                add_rows(&mut m, 0, j, &a, &b, &d_0, &d_j);
                if full { add_rows(&mut s, 0, j, &a, &b, &d_0, &d_j); }
                pivot = g;
            }
        }
        let mut pivot = m[0][0].clone();
        for j in 1..cols {
            if m[0][j].is_zero() { continue; }
            if pivot.is_zero() { return Ok(None); }
            let (d, r) = m[0][j].div(&pivot);
            if r.is_zero() {
                let d = d.neg();
                add_columns(&mut m, 0, j, &one, &zero, &d, &one);
                if full { add_columns(&mut t, 0, j, &one, &zero, &d, &one); }
            } else {
                let (a, b, g) = pivot.gcdex(&m[0][j]);
                let d_0 = m[0][j].exquo(&g);
                let d_j = pivot.exquo(&g).neg();
                add_columns(&mut m, 0, j, &a, &b, &d_0, &d_j);
                if full { add_columns(&mut t, 0, j, &a, &b, &d_0, &d_j); }
                pivot = g;
            }
        }
    }

    // c = canonical_unit(m[0][0]) (1/m[0][0] over a field); applied unless c == 1.
    if m[0][0].is_zero() && C::zero_pivot_raises() { return Ok(None); }
    if !m[0][0].is_zero() {
        if let Some(c) = m[0][0].unit_factor() {
            m[0][0] = m[0][0].mul(&c);
            if full { s[0] = s[0].iter().map(|e| e.mul(&c)).collect(); }
        }
    }

    let mut invs = Vec::new();
    if rows != 1 && cols != 1 {
        let lower: Matrix<C> = m[1..].iter().map(|row| row[1..].to_vec()).collect();
        let Some((small, transforms)) = recurse(lower, rows - 1, cols - 1, full, check)? else {
            return Ok(None);
        };
        invs = small;
        if let Some((s_small, t_small)) = transforms {
            s = product(&bordered(s_small, rows), &s);
            t = product(&t, &bordered(t_small, cols));
        }
    }

    let pivot = m[0][0].clone();
    let result = if !pivot.is_zero() {
        let mut result = Vec::with_capacity(invs.len() + 1);
        result.push(pivot);
        result.extend(invs);
        // In case m[0] doesn't divide the invariants of the rest of the matrix.
        for i in 0..result.len() - 1 {
            check()?;
            let (a, b) = (result[i].clone(), result[i + 1].clone());
            if a.is_zero() { return Ok(None); }
            if b.is_zero() || b.div(&a).1.is_zero() { break; }
            let (x, y, d) = if full { a.gcdex(&b) } else { (C::zero(), C::zero(), a.gcd(&b)) };
            let alpha = a.div(&d).0;
            if full {
                let beta = b.div(&d).0;
                add_rows(&mut s, i, i + 1, &one, &zero, &x, &one);
                add_columns(&mut t, i, i + 1, &one, &y, &zero, &one);
                add_rows(&mut s, i, i + 1, &one, &alpha.neg(), &zero, &one);
                add_columns(&mut t, i, i + 1, &one, &zero, &beta.neg(), &one);
                add_rows(&mut s, i, i + 1, &zero, &one, &one.neg(), &zero);
            }
            result[i + 1] = b.mul(&alpha);
            result[i] = d;
        }
        result
    } else {
        if full {
            if rows > 1 { s.rotate_left(1); }
            if cols > 1 { for row in t.iter_mut() { row.rotate_left(1); } }
        }
        invs.push(pivot);
        invs
    };
    Ok(Some((result, full.then_some((s, t)))))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ints(rows: &[&[i64]]) -> Matrix {
        rows.iter().map(|row| row.iter().copied().map(Int::from).collect()).collect()
    }
    fn run(rows: &[&[i64]], r: usize, c: usize) -> Decomposition {
        decompose(ints(rows), r, c, true, &mut || Ok::<_, ()>(())).unwrap().unwrap()
    }
    fn certificate(a: &Matrix, (invs, st): &Decomposition) {
        let (s, t) = st.as_ref().unwrap();
        let p = product(&product(s, a), t);
        for (i, row) in p.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let expected = if i == j && i < invs.len() { invs[i].clone() } else { Int::zero() };
                assert!(*value == expected);
            }
        }
    }

    #[test]
    fn reference_gcdex_signs() {
        let cases = [(4, 6, (-1, 1, 2)), (-4, 6, (1, 1, 2)), (6, -4, (1, 1, 2)), (0, -5, (0, -1, 5)),
            (-5, 0, (-1, 0, 5)), (0, 0, (0, 0, 0)), (3, 9, (1, 0, 3)), (-12, 18, (1, 1, 6))];
        for (a, b, (x, y, g)) in cases {
            let r = gcdex(&Int::from(a), &Int::from(b));
            assert!(r == (Int::from(x), Int::from(y), Int::from(g)));
        }
    }

    #[test]
    fn reference_docstring_decomposition() {
        let a = ints(&[&[12, 6, 4], &[3, 9, 6], &[2, 16, 14]]);
        let d = run(&[&[12, 6, 4], &[3, 9, 6], &[2, 16, 14]], 3, 3);
        assert!(d.0 == ints(&[&[1, 10, 30]])[0]);
        let (s, t) = d.1.clone().unwrap();
        assert!(s == ints(&[&[0, 1, -1], &[1, -4, 0], &[0, -2, 3]]));
        assert!(t == ints(&[&[1, 1, 10], &[0, -1, -2], &[0, 1, 3]]));
        certificate(&a, &d);
    }

    #[test]
    fn zero_rows_columns_and_empty_shapes() {
        for (rows, r, c) in [(vec![vec![0i64, 0], vec![0, 0]], 2, 2), (vec![vec![0, 0, 3]], 1, 3),
            (vec![vec![0], vec![-4]], 2, 1), (vec![vec![0, 0], vec![0, 6], vec![4, 0]], 3, 2)] {
            let a: Matrix = rows.iter().map(|row| row.iter().copied().map(Int::from).collect()).collect();
            let d = decompose(a.clone(), r, c, true, &mut || Ok::<_, ()>(())).unwrap().unwrap();
            certificate(&a, &d);
        }
        let d = decompose(Vec::<Vec<Int>>::new(), 0, 3, true, &mut || Ok::<_, ()>(())).unwrap().unwrap();
        assert!(d.0.is_empty());
        assert!(d.1.unwrap().1 == ints(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]));
        assert!(decompose(ints(&[&[1], &[2, 3]]), 2, 2, false, &mut || Ok::<_, ()>(())).unwrap().is_none());
    }
}
