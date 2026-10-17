//! Integer kernels of `sympy.discrete`: the number theoretic transform,
//! the Walsh-Hadamard and Moebius (zeta) transforms and linear recurrences
//! with integer data. Symbolic sequences stay on the Python side.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use super::isprime;
use super::residue::primitive_root;

type R<T> = Result<T, String>;

fn bitrev(i: usize, bits: u32) -> usize {
    if bits == 0 {
        return 0;
    }
    i.reverse_bits() >> (usize::BITS - bits)
}

fn pad_pow2<T: Clone>(a: &mut Vec<T>, zero: T) {
    let n = a.len();
    if n >= 2 && n & (n - 1) != 0 {
        a.resize(n.next_power_of_two(), zero);
    }
}

/// SymPy's `ntt(seq, prime)` / `intt(seq, prime)`.
pub fn ntt(seq: &[BigInt], p: &BigInt, inverse: bool) -> R<Vec<BigInt>> {
    if !isprime(p) {
        return Err("Expected prime modulus for Number Theoretic Transform".into());
    }
    let mut a: Vec<BigInt> = seq.iter().map(|x| x.mod_floor(p)).collect();
    let n0 = a.len();
    if n0 < 1 {
        return Ok(a);
    }
    let n = n0.next_power_of_two();
    let b = n.trailing_zeros();
    let pm1: BigInt = p - 1;
    if !(&pm1 % n).is_zero() {
        return Err("Expected prime modulus of the form (m*2**k + 1)".into());
    }
    a.resize(n, BigInt::zero());
    for i in 1..n {
        let j = bitrev(i, b);
        if i < j {
            a.swap(i, j);
        }
    }
    let pr = primitive_root(p, true)?.ok_or("no primitive root")?;
    let mut rt = pr.modpow(&(&pm1 / n), p);
    if inverse {
        rt = rt.modpow(&(p - 2), p);
    }
    let mut w = vec![BigInt::one(); (n / 2).max(1)];
    for i in 1..n / 2 {
        w[i] = (&w[i - 1] * &rt) % p;
    }
    let mut h = 2;
    while h <= n {
        let (hf, ut) = (h / 2, n / h);
        for i in (0..n).step_by(h) {
            for j in 0..hf {
                let u = a[i + j].clone();
                let v = &a[i + j + hf] * &w[ut * j];
                a[i + j] = (&u + &v).mod_floor(p);
                a[i + j + hf] = (u - v).mod_floor(p);
            }
        }
        h *= 2;
    }
    if inverse {
        let rv = BigInt::from(n).modpow(&(p - 2), p);
        for x in a.iter_mut() {
            *x = (&*x * &rv) % p;
        }
    }
    Ok(a)
}

/// Unnormalised Walsh-Hadamard transform of an integer sequence (padded to
/// a power of two). The inverse transform is this result divided by its
/// length.
pub fn fwht(seq: &[BigInt]) -> Vec<BigInt> {
    let mut a = seq.to_vec();
    if a.len() < 2 {
        return a;
    }
    pad_pow2(&mut a, BigInt::zero());
    let n = a.len();
    let mut h = 2;
    while h <= n {
        let hf = h / 2;
        for i in (0..n).step_by(h) {
            for j in 0..hf {
                let u = a[i + j].clone();
                let v = a[i + j + hf].clone();
                a[i + j] = &u + &v;
                a[i + j + hf] = u - v;
            }
        }
        h *= 2;
    }
    a
}

/// Zeta (`sgn = 1`) or Moebius (`sgn = -1`) transform over the subset or
/// superset lattice of an integer sequence.
pub fn mobius(seq: &[BigInt], sgn: i32, subset: bool) -> Vec<BigInt> {
    let mut a = seq.to_vec();
    if a.len() < 2 {
        return a;
    }
    pad_pow2(&mut a, BigInt::zero());
    let n = a.len();
    let mut i = 1;
    while i < n {
        for j in 0..n {
            if ((j & i) != 0) == subset {
                let t = a[j ^ i].clone();
                if sgn > 0 {
                    a[j] += t;
                } else {
                    a[j] -= t;
                }
            }
        }
        i *= 2;
    }
    a
}

/// SymPy's `linrec_coeffs(c, n)` with integer coefficients.
pub fn linrec_coeffs(c: &[BigInt], n: &BigInt) -> Vec<BigInt> {
    let k = c.len();
    if k == 0 {
        return vec![];
    }
    // Exponent digits from the top until the remaining prefix drops below k.
    let mut prefixes = vec![n.clone()];
    while prefixes.last().unwrap() >= &BigInt::from(k) {
        let last = prefixes.last().unwrap().clone();
        prefixes.push(last >> 1);
    }
    let base = prefixes.pop().unwrap();
    let base: usize = (&base).try_into().unwrap_or(0);
    let mut u = vec![BigInt::zero(); k];
    if base < k {
        u[base] = BigInt::one();
    }
    while let Some(m) = prefixes.pop() {
        let offset = if m.is_odd() { 1 } else { 0 };
        let mut w = vec![BigInt::zero(); 2 * u.len() - 1 + offset];
        for (i, p) in u.iter().enumerate() {
            if p.is_zero() {
                continue;
            }
            for (j, q) in u.iter().enumerate() {
                w[offset + i + j] += p * q;
            }
        }
        for j in (k..w.len()).rev() {
            let wj = w[j].clone();
            if wj.is_zero() {
                continue;
            }
            for (i, ci) in c.iter().enumerate() {
                w[j - i - 1] += &wj * ci;
            }
        }
        w.truncate(k);
        u = w;
    }
    u
}

/// SymPy's `linrec(coeffs, init, n)` with integer data (`init` already
/// padded to the order, `n >= 0`).
pub fn linrec(c: &[BigInt], init: &[BigInt], n: &BigInt) -> BigInt {
    let k = c.len();
    if n < &BigInt::from(k) {
        let i: usize = n.try_into().unwrap_or(0);
        return init[i].clone();
    }
    linrec_coeffs(c, n).iter().zip(init).map(|(u, v)| u * v).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(xs: &[i64]) -> Vec<BigInt> {
        xs.iter().map(|&x| BigInt::from(x)).collect()
    }

    #[test]
    fn transforms() {
        let p = BigInt::from(3 * 2i64.pow(8) + 1);
        let t = ntt(&v(&[1, 2, 3, 4]), &p, false).unwrap();
        assert_eq!(t, v(&[10, 643, 767, 122]));
        assert_eq!(ntt(&t, &p, true).unwrap(), v(&[1, 2, 3, 4]));
        assert_eq!(fwht(&v(&[4, 2, 2, 0, 0, 2, -2, 0])), v(&[8, 0, 8, 0, 8, 8, 0, 0]));
        assert_eq!(mobius(&v(&[1, 2, 3]), 1, true), v(&[1, 3, 4, 6]));
        assert_eq!(mobius(&v(&[1, 2, 3]), 1, false), v(&[6, 2, 3, 0]));
        assert_eq!(linrec(&v(&[1, 1]), &v(&[0, 1]), &BigInt::from(10)), BigInt::from(55));
        assert_eq!(linrec(&v(&[1, 1]), &v(&[0, 1]), &BigInt::from(1)), BigInt::from(1));
        assert_eq!(linrec(&v(&[2, -1, 3]), &v(&[1, 0, 2]), &BigInt::from(7)), BigInt::from(123));
    }
}
