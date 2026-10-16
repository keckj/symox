//! Multimodular (CRT) algorithms for integer matrices: reduced row echelon
//! form and inverse. Results are exact and deterministic: the number of
//! primes is chosen from Hadamard's bound, unlucky primes are detected,
//! and the rank is certified (otherwise the callers fall back to
//! fraction-free elimination).

use std::sync::OnceLock;

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{ToPrimitive, Zero};

use super::int::Int;

// ----------------------------------------------------------------------
// primes and Montgomery arithmetic
// ----------------------------------------------------------------------

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % p == 0 {
            return n == p;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// The `k`-th largest prime below `2^62` (generated lazily, cached).
pub fn prime(k: usize) -> u64 {
    use std::sync::Mutex;
    static PRIMES: OnceLock<Mutex<Vec<u64>>> = OnceLock::new();
    let cell = PRIMES.get_or_init(|| Mutex::new(Vec::new()));
    let mut ps = cell.lock().unwrap();
    while ps.len() <= k {
        let mut n = match ps.last() {
            Some(&p) => p - 2,
            None => (1u64 << 62) - 1,
        };
        while !is_prime(n) {
            n -= 2;
        }
        ps.push(n);
    }
    ps[k]
}

/// Montgomery arithmetic modulo an odd `p < 2^63`.
#[derive(Clone, Copy, Debug)]
pub struct Mont {
    pub p: u64,
    pinv: u64,
    r2: u64,
}

impl Mont {
    pub fn new(p: u64) -> Mont {
        // pinv = -p^{-1} mod 2^64 (Newton iteration)
        let mut inv: u64 = 1;
        for _ in 0..7 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let pinv = inv.wrapping_neg();
        let r = ((1u128 << 64) % p as u128) as u64;
        let r2 = mulmod(r, r, p);
        Mont { p, pinv, r2 }
    }

    #[inline(always)]
    pub fn redc(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.pinv);
        let u = ((t + m as u128 * self.p as u128) >> 64) as u64;
        if u >= self.p {
            u - self.p
        } else {
            u
        }
    }

    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.redc(a as u128 * b as u128)
    }

    #[inline(always)]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }

    #[inline(always)]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }

    pub fn to_mont(&self, a: u64) -> u64 {
        self.mul(a % self.p, self.r2)
    }

    pub fn from_mont(&self, a: u64) -> u64 {
        self.redc(a as u128)
    }

    /// Inverse of a Montgomery-form element.
    pub fn inv(&self, a: u64) -> u64 {
        let x = self.from_mont(a);
        let i = powmod(x, self.p - 2, self.p);
        self.to_mont(i)
    }

    pub fn one(&self) -> u64 {
        self.to_mont(1)
    }

    /// Reduce an integer and convert to Montgomery form.
    pub fn reduce(&self, x: &Int) -> u64 {
        let r = match x {
            Int::S(v) => v.rem_euclid(self.p as i64) as u64,
            Int::B(b) => {
                let m = BigInt::from(self.p);
                let r = ((b % &m) + &m) % &m;
                r.to_u64().unwrap()
            }
        };
        self.to_mont(r)
    }
}

// ----------------------------------------------------------------------
// elimination modulo p
// ----------------------------------------------------------------------

/// In-place reduced row echelon form modulo `p` (Montgomery form).
/// Returns `(pivot columns, original row index of every row after the
/// row swaps)`.
fn rref_mod(a: &mut [u64], rows: usize, cols: usize, mp: &Mont) -> (Vec<usize>, Vec<usize>) {
    let mut order: Vec<usize> = (0..rows).collect();
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(pr) = (r..rows).find(|&i| a[i * cols + c] != 0) else {
            continue;
        };
        if pr != r {
            for k in 0..cols {
                a.swap(pr * cols + k, r * cols + k);
            }
            order.swap(pr, r);
        }
        let inv = mp.inv(a[r * cols + c]);
        for k in c..cols {
            a[r * cols + k] = mp.mul(a[r * cols + k], inv);
        }
        let (before, rest) = a.split_at_mut(r * cols);
        let (prow, after) = rest.split_at_mut(cols);
        let prow: &[u64] = prow;
        let elim = |row: &mut [u64]| {
            let f = row[c];
            if f == 0 {
                return;
            }
            for k in c..cols {
                let v = prow[k];
                if v != 0 {
                    row[k] = mp.sub(row[k], mp.mul(f, v));
                }
            }
        };
        for row in before.chunks_mut(cols) {
            elim(row);
        }
        for row in after.chunks_mut(cols) {
            elim(row);
        }
        pivots.push(c);
        r += 1;
    }
    (pivots, order)
}

/// Determinant of the leading `n x n` (row-major, Montgomery form) matrix
/// and its inverse appended as columns (in-place Gauss-Jordan on
/// `[A | I]`). Returns `None` when singular modulo `p`.
fn inv_mod(a: &[u64], n: usize, mp: &Mont) -> Option<(u64, Vec<u64>)> {
    let w = 2 * n;
    let one = mp.one();
    let mut m = vec![0u64; n * w];
    for i in 0..n {
        m[i * w..i * w + n].copy_from_slice(&a[i * n..(i + 1) * n]);
        m[i * w + n + i] = one;
    }
    let mut det = one;
    for c in 0..n {
        let pr = (c..n).find(|&i| m[i * w + c] != 0)?;
        if pr != c {
            for k in 0..w {
                m.swap(pr * w + k, c * w + k);
            }
            det = mp.sub(0, det);
        }
        let pv = m[c * w + c];
        det = mp.mul(det, pv);
        let inv = mp.inv(pv);
        for k in c..w {
            m[c * w + k] = mp.mul(m[c * w + k], inv);
        }
        let (before, rest) = m.split_at_mut(c * w);
        let (prow, after) = rest.split_at_mut(w);
        let prow: &[u64] = prow;
        let elim = |row: &mut [u64]| {
            let f = row[c];
            if f == 0 {
                return;
            }
            for k in c..w {
                let v = prow[k];
                if v != 0 {
                    row[k] = mp.sub(row[k], mp.mul(f, v));
                }
            }
        };
        for row in before.chunks_mut(w) {
            elim(row);
        }
        for row in after.chunks_mut(w) {
            elim(row);
        }
    }
    let mut inv = Vec::with_capacity(n * n);
    for i in 0..n {
        inv.extend_from_slice(&m[i * w + n..(i + 1) * w]);
    }
    Some((det, inv))
}

// ----------------------------------------------------------------------
// CRT
// ----------------------------------------------------------------------

/// Incremental Chinese remaindering of a vector of residues.
struct Crt {
    modulus: BigUint,
    vals: Vec<BigUint>,
}

impl Crt {
    fn new(n: usize) -> Crt {
        Crt { modulus: BigUint::from(1u32), vals: vec![BigUint::zero(); n] }
    }

    fn bits(&self) -> u64 {
        self.modulus.bits()
    }

    /// Add residues `r` (plain form) modulo the prime `p`.
    fn add(&mut self, r: &[u64], p: u64) {
        let mmod = (&self.modulus % p).to_u64().unwrap();
        let minv = powmod(mmod, p - 2, p);
        for (x, &ri) in self.vals.iter_mut().zip(r) {
            let xm = (&*x % p).to_u64().unwrap();
            let diff = if ri >= xm { ri - xm } else { ri + p - xm };
            let k = mulmod(diff, minv, p);
            if k != 0 {
                *x += &self.modulus * k;
            }
        }
        self.modulus *= p;
    }

    /// Symmetric representatives.
    fn finish(self) -> Vec<Int> {
        let half = &self.modulus >> 1;
        let m = BigInt::from_biguint(Sign::Plus, self.modulus.clone());
        self.vals
            .into_iter()
            .map(|x| {
                if x > half {
                    Int::from_big(BigInt::from_biguint(Sign::Plus, x) - &m)
                } else {
                    Int::from_big(BigInt::from_biguint(Sign::Plus, x))
                }
            })
            .collect()
    }
}

/// `log2` of the Euclidean norm of a row (upper bound).
fn row_norm_bits(row: &[Int]) -> f64 {
    let mut s = 0f64;
    let mut big = false;
    let mut maxb = 0u64;
    for x in row {
        let b = x.bits();
        maxb = maxb.max(b);
        if b > 500 {
            big = true;
        } else {
            let f = x.to_f64();
            s += f * f;
        }
    }
    if big {
        maxb as f64 + 0.5 * (row.len() as f64).log2() + 1.0
    } else if s == 0.0 {
        0.0
    } else {
        0.5 * s.log2() * (1.0 + 1e-12) + 1e-9
    }
}

// ----------------------------------------------------------------------
// public kernels
// ----------------------------------------------------------------------

/// Reduced row echelon form of an integer matrix: `(X, D, pivots)` with
/// `rref[:r] = X / D` (`X` is `r x cols`). `None` when the modular
/// computation could not be certified (use fraction-free elimination).
pub fn rref_int(a: &[Int], rows: usize, cols: usize) -> Option<(Vec<Int>, Int, Vec<usize>)> {
    if rows == 0 || cols == 0 {
        return Some((vec![], Int::ONE, vec![]));
    }
    let mp0 = Mont::new(prime(0));
    let mut w: Vec<u64> = a.iter().map(|x| mp0.reduce(x)).collect();
    let (pivots, order) = rref_mod(&mut w, rows, cols, &mp0);
    let r = pivots.len();
    if r == 0 {
        // certify: the matrix must be exactly zero
        return if a.iter().all(|x| x.is_zero()) { Some((vec![], Int::ONE, vec![])) } else { None };
    }
    let srows: Vec<usize> = order[..r].to_vec();
    let free: Vec<usize> = (0..cols).filter(|c| !pivots.contains(c)).collect();
    // X has identity (times D) in the pivot columns: only the free columns
    // and D need reconstruction
    let mut x_free: Vec<Int> = Vec::new();
    let mut d = Int::ONE;
    if !free.is_empty() {
        let bits: f64 = srows.iter().map(|&i| row_norm_bits(&a[i * cols..(i + 1) * cols])).sum();
        let need = bits.ceil() as u64 + 2;
        let sub: Vec<Int> = srows.iter().flat_map(|&i| a[i * cols..(i + 1) * cols].iter().cloned()).collect();
        let nvals = r * free.len() + 1;
        let mut crt = Crt::new(nvals);
        let mut k = 0usize;
        let mut misses = 0;
        while crt.bits() <= need {
            let mp = Mont::new(prime(k));
            k += 1;
            let mut m: Vec<u64> = sub.iter().map(|x| mp.reduce(x)).collect();
            let (pv, _) = rref_mod(&mut m, r, cols, &mp);
            if pv != pivots {
                misses += 1;
                if misses > 8 {
                    return None;
                }
                continue;
            }
            // D mod p = det(A[S, P]) mod p
            let mut dm: Vec<u64> = Vec::with_capacity(r * r);
            for i in 0..r {
                for &c in &pivots {
                    dm.push(mp.reduce(&sub[i * cols + c]));
                }
            }
            let Some((dp, _)) = det_mod(&mut dm, r, &mp) else {
                misses += 1;
                if misses > 8 {
                    return None;
                }
                continue;
            };
            let mut res = Vec::with_capacity(nvals);
            for i in 0..r {
                for &c in &free {
                    res.push(mp.from_mont(mp.mul(m[i * cols + c], dp)));
                }
            }
            res.push(mp.from_mont(dp));
            crt.add(&res, mp.p);
        }
        let vals = crt.finish();
        d = vals[nvals - 1].clone();
        x_free = vals[..nvals - 1].to_vec();
        if d.is_zero() {
            return None;
        }
    }
    // assemble X
    let mut x = vec![Int::ZERO; r * cols];
    for (i, &pc) in pivots.iter().enumerate() {
        x[i * cols + pc] = d.clone();
    }
    for i in 0..r {
        for (fj, &c) in free.iter().enumerate() {
            x[i * cols + c] = x_free[i * free.len() + fj].clone();
        }
    }
    // certify the rank: every other row of A must be a combination of the
    // rows of X / D:  D * A[i, :] == sum_k A[i, P_k] * X[k, :]
    if r < rows {
        for i in (0..rows).filter(|i| !srows.contains(i)) {
            let row = &a[i * cols..(i + 1) * cols];
            for j in &free {
                let mut s = Int::ZERO;
                for (kk, &pc) in pivots.iter().enumerate() {
                    let f = &row[pc];
                    if !f.is_zero() {
                        s = s.add(&f.mul(&x[kk * cols + j]));
                    }
                }
                if s != d.mul(&row[*j]) {
                    return None;
                }
            }
            // pivot columns hold trivially: sum_k A[i,P_k] * D*delta = D*A[i,P_j]
        }
    }
    Some((x, d, pivots))
}

fn det_mod(m: &mut [u64], n: usize, mp: &Mont) -> Option<(u64, ())> {
    let mut det = mp.one();
    for c in 0..n {
        let pr = (c..n).find(|&i| m[i * n + c] != 0)?;
        if pr != c {
            for k in 0..n {
                m.swap(pr * n + k, c * n + k);
            }
            det = mp.sub(0, det);
        }
        let pv = m[c * n + c];
        det = mp.mul(det, pv);
        let inv = mp.inv(pv);
        for i in c + 1..n {
            let f = mp.mul(m[i * n + c], inv);
            if f == 0 {
                continue;
            }
            for k in c..n {
                let v = m[c * n + k];
                if v != 0 {
                    m[i * n + k] = mp.sub(m[i * n + k], mp.mul(f, v));
                }
            }
        }
    }
    Some((det, ()))
}

/// Exact inverse of an integer matrix as `(adj, det)` with
/// `A^{-1} = adj / det`. `Some(None)` when singular, `None` when the
/// modular computation gave up.
pub fn inv_int(a: &[Int], n: usize) -> Option<Option<(Vec<Int>, Int)>> {
    let bits: f64 = (0..n).map(|i| row_norm_bits(&a[i * n..(i + 1) * n])).sum();
    if (0..n).any(|i| a[i * n..(i + 1) * n].iter().all(|x| x.is_zero())) {
        return Some(None);
    }
    let need = bits.ceil() as u64 + 2;
    let mut crt = Crt::new(n * n + 1);
    let mut k = 0usize;
    let mut misses = 0;
    while crt.bits() <= need {
        let mp = Mont::new(prime(k));
        k += 1;
        let m: Vec<u64> = a.iter().map(|x| mp.reduce(x)).collect();
        let Some((det, inv)) = inv_mod(&m, n, &mp) else {
            misses += 1;
            if misses >= 3 && crt.bits() <= 1 {
                // probably singular: decide exactly
                let d = super::ff::bareiss_det(a.to_vec(), n, &Int::ONE);
                if d.is_zero() {
                    return Some(None);
                }
            }
            if misses > 16 {
                return None;
            }
            continue;
        };
        let mut res: Vec<u64> = inv.iter().map(|&v| mp.from_mont(mp.mul(v, det))).collect();
        res.push(mp.from_mont(det));
        crt.add(&res, mp.p);
    }
    let mut vals = crt.finish();
    let det = vals.pop().unwrap();
    if det.is_zero() {
        return Some(None);
    }
    Some(Some((vals, det)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn montgomery() {
        let p = prime(0);
        let mp = Mont::new(p);
        let a = mp.to_mont(123456789);
        let b = mp.to_mont(987654321);
        assert_eq!(mp.from_mont(mp.mul(a, b)), mulmod(123456789, 987654321, p));
        assert_eq!(mp.from_mont(mp.mul(a, mp.inv(a))), 1);
        assert!(is_prime(p) && p < (1 << 62));
    }

    #[test]
    fn rref_and_inverse() {
        let v: Vec<Int> = [1, 2, 3, 4, 5, 6, 7, 8, 9].iter().map(|&x| Int::from(x)).collect();
        let (x, d, piv) = rref_int(&v, 3, 3).unwrap();
        assert_eq!(piv, vec![0, 1]);
        // rref row 0 = [1, 0, -1]
        assert_eq!(Rat_(&x[2], &d), (-1, 1));
        let w: Vec<Int> = [2, 1, 1, 3].iter().map(|&x| Int::from(x)).collect();
        let (adj, det) = inv_int(&w, 2).unwrap().unwrap();
        assert_eq!(det, Int::from(5));
        assert_eq!(adj, [3, -1, -1, 2].iter().map(|&x| Int::from(x)).collect::<Vec<_>>());
        assert!(inv_int(&v, 3).unwrap().is_none());
    }

    #[allow(non_snake_case)]
    fn Rat_(a: &Int, b: &Int) -> (i64, i64) {
        let r = super::super::int::Rat::new(a.clone(), b.clone());
        (r.num.to_i64().unwrap(), r.den.to_i64().unwrap())
    }
}
