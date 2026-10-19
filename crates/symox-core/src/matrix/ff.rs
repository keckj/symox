//! Fraction-free algorithms over an integral domain with exact division:
//! Bareiss determinant, fraction-free Gauss-Jordan (`rref_den`) and the
//! division-free Berkowitz characteristic polynomial.

use super::int::Int;
use super::mpoly::MPoly;

/// Integral domain with exact division.
pub trait Ring: Clone + PartialEq {
    fn zero_like(&self) -> Self;
    fn one_like(&self) -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    /// Exact division (`o` divides `self`).
    fn exquo(&self, o: &Self) -> Self;
    /// `(self*a - o*b) / d`.
    /// `sum a_i * b_i` over the pairs.
    fn sum_products(pairs: &[(&Self, &Self)], zero: &Self) -> Self {
        let mut acc = zero.clone();
        for (x, y) in pairs {
            acc = acc.add(&x.mul(y));
        }
        acc
    }
    fn mul_sub_div(&self, a: &Self, o: &Self, b: &Self, d: &Self) -> Self {
        let r = self.mul(a).sub(&o.mul(b));
        if d.is_one() {
            r
        } else {
            r.exquo(d)
        }
    }
}

impl Ring for Int {
    fn zero_like(&self) -> Int {
        Int::ZERO
    }
    fn one_like(&self) -> Int {
        Int::ONE
    }
    fn is_zero(&self) -> bool {
        Int::is_zero(self)
    }
    fn is_one(&self) -> bool {
        Int::is_one(self)
    }
    fn add(&self, o: &Int) -> Int {
        Int::add(self, o)
    }
    fn sub(&self, o: &Int) -> Int {
        Int::sub(self, o)
    }
    fn mul(&self, o: &Int) -> Int {
        Int::mul(self, o)
    }
    fn neg(&self) -> Int {
        Int::neg(self)
    }
    fn exquo(&self, o: &Int) -> Int {
        self.divexact(o)
    }
    fn mul_sub_div(&self, a: &Int, o: &Int, b: &Int, d: &Int) -> Int {
        Int::mul_sub_div(self, a, o, b, d)
    }
}

impl Ring for MPoly {
    fn zero_like(&self) -> MPoly {
        MPoly::zero(self.nvars)
    }
    fn one_like(&self) -> MPoly {
        MPoly::one(self.nvars)
    }
    fn is_zero(&self) -> bool {
        MPoly::is_zero(self)
    }
    fn is_one(&self) -> bool {
        MPoly::is_one(self)
    }
    fn add(&self, o: &MPoly) -> MPoly {
        MPoly::add(self, o)
    }
    fn sub(&self, o: &MPoly) -> MPoly {
        MPoly::sub(self, o)
    }
    fn mul(&self, o: &MPoly) -> MPoly {
        MPoly::mul(self, o)
    }
    fn neg(&self) -> MPoly {
        MPoly::neg(self)
    }
    fn exquo(&self, o: &MPoly) -> MPoly {
        self.divexact(o).expect("inexact division in fraction-free elimination")
    }
    fn sum_products(pairs: &[(&MPoly, &MPoly)], zero: &MPoly) -> MPoly {
        MPoly::sum_products(pairs, zero.nvars)
    }
}

/// Bareiss determinant of the `n x n` row-major matrix `a`, following
/// SymPy's recursion (pivot = first non-zero entry of the first column,
/// the pivot row is removed keeping the order of the other rows).
pub fn bareiss_det<R: Ring>(mut a: Vec<R>, n: usize, one: &R) -> R {
    if n == 0 {
        return one.clone();
    }
    let mut size = n;
    let mut cumm = one.clone();
    let mut negate = false;
    loop {
        if size == 1 {
            let d = a[0].clone();
            return if negate { d.neg() } else { d };
        }
        let piv = match (0..size).find(|&i| !a[i * size].is_zero()) {
            Some(p) => p,
            None => return one.zero_like(),
        };
        if piv % 2 == 1 {
            negate = !negate;
        }
        let pv = a[piv * size].clone();
        let m = size - 1;
        let mut b: Vec<R> = Vec::with_capacity(m * m);
        for i in 0..size {
            if i == piv {
                continue;
            }
            let ai0 = &a[i * size];
            for j in 0..m {
                let v = pv.mul_sub_div(&a[i * size + j + 1], &a[piv * size + j + 1], ai0, &cumm);
                b.push(v);
            }
        }
        a = b;
        size = m;
        cumm = pv;
    }
}

/// Fraction-free Gauss-Jordan elimination. Returns `(den, pivots)` with
/// `a` overwritten so that `a / den` is the reduced row echelon form.
pub fn rref_den<R: Ring>(a: &mut [R], rows: usize, cols: usize, one: &R) -> (R, Vec<usize>) {
    let mut pivots = Vec::new();
    let mut prev = one.clone();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p) = (r..rows).find(|&i| !a[i * cols + c].is_zero()) else {
            continue;
        };
        if p != r {
            for k in 0..cols {
                a.swap(p * cols + k, r * cols + k);
            }
        }
        let pv = a[r * cols + c].clone();
        for i in 0..rows {
            if i == r {
                continue;
            }
            let f = a[i * cols + c].clone();
            if f.is_zero() {
                // entries still need the scaling by pv / prev
                if !(pv == prev) {
                    for k in 0..cols {
                        let v = &a[i * cols + k];
                        if !v.is_zero() {
                            a[i * cols + k] = v.mul(&pv).exquo(&prev);
                        }
                    }
                }
                continue;
            }
            for k in 0..cols {
                let v = a[i * cols + k].mul_sub_div(&pv, &a[r * cols + k], &f, &prev);
                a[i * cols + k] = v;
            }
        }
        prev = pv;
        pivots.push(c);
        r += 1;
    }
    // rows of the pivot block are scaled by `prev`: entries a/den with
    // den = prev for all pivot rows (the pivots all equal `prev`)
    (prev, pivots)
}

/// Berkowitz algorithm: coefficients `[1, c_1, ..., c_n]` of
/// `det(t*I - A) = t^n + c_1 t^{n-1} + ... + c_n`.
pub fn berkowitz<R: Ring>(a: &[R], n: usize, one: &R) -> Vec<R> {
    let zero = one.zero_like();
    if n == 0 {
        return vec![one.clone()];
    }
    // vector of the characteristic polynomial of the leading k x k block
    let mut v: Vec<R> = vec![one.clone(), a[0].neg()];
    for k in 1..n {
        // A_{k+1} = [[A_k, C], [R, a_kk]] with R row k (cols 0..k), C col k (rows 0..k)
        let r: Vec<R> = (0..k).map(|j| a[k * n + j].clone()).collect();
        let c: Vec<R> = (0..k).map(|i| a[i * n + k].clone()).collect();
        let akk = a[k * n + k].clone();
        // Toeplitz column: t_0 = 1, t_1 = -a_kk, t_{i+2} = -R A^i C
        let mut t: Vec<R> = Vec::with_capacity(k + 2);
        t.push(one.clone());
        t.push(akk.neg());
        let mut x = c.clone(); // A^i C
        for i in 0..k {
            let s = dot(&r, &x, &zero);
            t.push(s.neg());
            if i + 1 < k {
                // x = A_k x
                let mut y = Vec::with_capacity(k);
                for row in 0..k {
                    let pairs: Vec<(&R, &R)> = (0..k)
                        .map(|col| (&a[row * n + col], &x[col]))
                        .filter(|(e, xc)| !e.is_zero() && !xc.is_zero())
                        .collect();
                    y.push(R::sum_products(&pairs, &zero));
                }
                x = y;
            }
        }
        // new vector = T * v, T lower triangular Toeplitz (k+2) x (k+1)
        let mut nv: Vec<R> = Vec::with_capacity(k + 2);
        for i in 0..k + 2 {
            let pairs: Vec<(&R, &R)> = (0..=i.min(k))
                .filter(|&j| i - j < t.len())
                .map(|j| (&t[i - j], &v[j]))
                .filter(|(tt, vj)| !tt.is_zero() && !vj.is_zero())
                .collect();
            nv.push(R::sum_products(&pairs, &zero));
        }
        v = nv;
    }
    v
}

fn dot<R: Ring>(a: &[R], b: &[R], zero: &R) -> R {
    let pairs: Vec<(&R, &R)> = a.iter().zip(b).filter(|(x, y)| !x.is_zero() && !y.is_zero()).collect();
    R::sum_products(&pairs, zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ints(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&x| Int::from(x)).collect()
    }

    #[test]
    fn det_int() {
        let a = ints(&[2, 0, 1, 1, 3, 2, 1, 1, 2]);
        assert_eq!(bareiss_det(a, 3, &Int::ONE), Int::from(6));
        let b = ints(&[0, 1, 1, 0]);
        assert_eq!(bareiss_det(b, 2, &Int::ONE), Int::from(-1));
        let c = ints(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(bareiss_det(c, 3, &Int::ONE), Int::from(0));
    }

    #[test]
    fn rref_int() {
        let mut a = ints(&[1, 2, 3, 4, 5, 6, 7, 8, 10]);
        let (den, piv) = rref_den(&mut a, 3, 3, &Int::ONE);
        assert_eq!(piv, vec![0, 1, 2]);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[i * 3 + j], if i == j { den.clone() } else { Int::ZERO });
            }
        }
        let mut b = ints(&[1, 2, 3, 2, 4, 6, 1, 0, 1]);
        let (den, piv) = rref_den(&mut b, 3, 3, &Int::ONE);
        assert_eq!(piv, vec![0, 1]);
        // rref = [[1,0,1],[0,1,1],[0,0,0]]
        assert_eq!(b[2].divexact(&den), Int::ONE);
        assert_eq!(b[5].divexact(&den), Int::ONE);
        assert!(b[6..].iter().all(|x| x.is_zero()));
    }

    #[test]
    fn charpoly_int() {
        let a = ints(&[1, 2, 3, 4]);
        let cp = berkowitz(&a, 2, &Int::ONE);
        assert_eq!(cp, ints(&[1, -5, -2]));
        let b = ints(&[2, 0, 0, 0, 3, 4, 0, 4, 9]);
        let cp = berkowitz(&b, 3, &Int::ONE);
        // (t-2)(t^2-12t+11) = t^3 -14t^2 +35t -22
        assert_eq!(cp, ints(&[1, -14, 35, -22]));
    }
}
