//! Exact integer and rational sequences used by the special functions.

use std::sync::{Mutex, OnceLock};

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// `n!`
pub fn factorial(n: u64) -> BigInt {
    crate::ntheory::factorial(n)
}

/// Double factorial `n!!` for `n >= 0`.
pub fn factorial2(n: u64) -> BigInt {
    fn prod(a: u64, b: u64) -> BigInt {
        // product of a, a+2, ..., b (same parity)
        if a > b {
            return BigInt::one();
        }
        let cnt = (b - a) / 2 + 1;
        if cnt <= 16 {
            let mut r = BigInt::one();
            let mut i = a;
            while i <= b {
                r *= i;
                i += 2;
            }
            return r;
        }
        let m = a + 2 * (cnt / 2);
        prod(a, m - 2) * prod(m, b)
    }
    if n < 2 {
        return BigInt::one();
    }
    prod(if n % 2 == 0 { 2 } else { 1 }, n)
}

/// Subfactorial (number of derangements) `!n`.
pub fn subfactorial(n: u64) -> BigInt {
    let (mut a, mut b) = (BigInt::one(), BigInt::zero()); // !0, !1
    if n == 0 {
        return a;
    }
    for k in 2..=n {
        let c = (&a + &b) * BigInt::from(k - 1);
        a = b;
        b = c;
    }
    b
}

/// Binomial coefficient for integers (SymPy semantics for negative `n`).
pub fn binomial(n: &BigInt, k: &BigInt) -> BigInt {
    crate::ntheory::binomial(n, k)
}

fn bernoulli_cache() -> &'static Mutex<Vec<BigRational>> {
    static C: OnceLock<Mutex<Vec<BigRational>>> = OnceLock::new();
    C.get_or_init(|| Mutex::new(Vec::new()))
}

/// Bernoulli number `B_n` with SymPy >= 1.12 convention `B_1 = +1/2`.
pub fn bernoulli_number(n: usize) -> BigRational {
    if n == 1 {
        return BigRational::new(BigInt::one(), BigInt::from(2));
    }
    if n > 1 && n % 2 == 1 {
        return BigRational::zero();
    }
    let mut c = bernoulli_cache().lock().unwrap();
    if c.len() <= n {
        // recompute with headroom (Akiyama-Tanigawa on even indices is
        // not simpler, use the full table)
        let m = (n + 1).max(2 * c.len()).max(32);
        let mut b = crate::summation::bernoulli_numbers(m);
        if b.len() > 1 {
            b[1] = BigRational::new(BigInt::one(), BigInt::from(2));
        }
        *c = b;
    }
    c[n].clone()
}

/// Euler (secant) number `E_n` (`E_0 = 1, E_2 = -1, E_4 = 5`, odd ones 0).
pub fn euler_number(n: usize) -> BigInt {
    if n % 2 == 1 {
        return BigInt::zero();
    }
    // Seidel's boustrophedon gives the zigzag numbers A_m; E_{2k} = (-1)^k A_{2k}
    let mut row = vec![BigInt::one()];
    for m in 1..=n {
        let mut next = Vec::with_capacity(m + 1);
        next.push(BigInt::zero());
        // alternate direction
        for j in 0..m {
            let v = &next[j] + &row[m - 1 - j];
            next.push(v);
        }
        row = next;
    }
    let a = row.last().unwrap().clone();
    if (n / 2) % 2 == 1 {
        -a
    } else {
        a
    }
}

/// Catalan number.
pub fn catalan_number(n: u64) -> BigInt {
    let b = binomial(&BigInt::from(2 * n), &BigInt::from(n));
    b / BigInt::from(n + 1)
}

/// `(F_n, F_{n+1})` by fast doubling.
fn fib_pair(n: u64) -> (BigInt, BigInt) {
    if n == 0 {
        return (BigInt::zero(), BigInt::one());
    }
    let (a, b) = fib_pair(n / 2);
    let c = &a * (&b * 2 - &a);
    let d = &a * &a + &b * &b;
    if n % 2 == 0 {
        (c, d)
    } else {
        let e = &c + &d;
        (d, e)
    }
}

/// Fibonacci number for any integer index.
pub fn fibonacci_number(n: i64) -> BigInt {
    let f = fib_pair(n.unsigned_abs()).0;
    if n < 0 && n % 2 == 0 {
        -f
    } else {
        f
    }
}

/// Lucas number for any integer index.
pub fn lucas_number(n: i64) -> BigInt {
    let m = n.unsigned_abs();
    let (a, b) = fib_pair(m);
    // L_m = 2 F_{m+1} - F_m
    let l: BigInt = b * 2 - a;
    if n < 0 && m % 2 == 1 {
        -l
    } else {
        l
    }
}

/// Bell number `B_n`.
pub fn bell_number(n: u64) -> BigInt {
    // Bell triangle
    let mut row = vec![BigInt::one()];
    for _ in 0..n {
        let mut next = Vec::with_capacity(row.len() + 1);
        next.push(row.last().unwrap().clone());
        for v in &row {
            let x = next.last().unwrap() + v;
            next.push(x);
        }
        row = next;
    }
    row[0].clone()
}

fn partition_cache() -> &'static Mutex<Vec<BigInt>> {
    static C: OnceLock<Mutex<Vec<BigInt>>> = OnceLock::new();
    C.get_or_init(|| Mutex::new(vec![BigInt::one()]))
}

/// Number of integer partitions `p(n)` (Euler's pentagonal recurrence).
pub fn partition_number(n: u64) -> BigInt {
    let n = n as usize;
    let mut p = partition_cache().lock().unwrap();
    while p.len() <= n {
        let m = p.len();
        let mut s = BigInt::zero();
        let mut k = 1usize;
        loop {
            let g1 = k * (3 * k - 1) / 2;
            if g1 > m {
                break;
            }
            let pos = k % 2 == 1;
            if pos {
                s += &p[m - g1];
            } else {
                s -= &p[m - g1];
            }
            let g2 = k * (3 * k + 1) / 2;
            if g2 <= m {
                if pos {
                    s += &p[m - g2];
                } else {
                    s -= &p[m - g2];
                }
            }
            k += 1;
        }
        p.push(s);
    }
    p[n].clone()
}

/// Unsigned Stirling numbers of the first kind `c(n, k)`.
pub fn stirling1(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let k = k as usize;
    let mut row = vec![BigInt::zero(); k + 1];
    row[0] = BigInt::one();
    for i in 0..n {
        for j in (0..=k.min(i as usize + 1)).rev() {
            let mut v = &row[j] * BigInt::from(i);
            if j > 0 {
                v += &row[j - 1];
            }
            row[j] = v;
        }
    }
    row[k].clone()
}

/// Stirling numbers of the second kind `S(n, k)`.
pub fn stirling2(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    if n == 0 {
        return BigInt::one();
    }
    // explicit formula S(n,k) = 1/k! sum (-1)^j C(k,j) (k-j)^n
    let mut s = BigInt::zero();
    let mut c = BigInt::one();
    for j in 0..=k {
        let t = &c * num_traits::pow(BigInt::from(k - j), n as usize);
        if j % 2 == 0 {
            s += t;
        } else {
            s -= t;
        }
        c = c * BigInt::from(k - j) / BigInt::from(j + 1);
    }
    s / factorial(k)
}

/// SymPy's `stirling(n, k, d=None, kind=2, signed=False)`.
pub fn stirling(n: u64, k: u64, d: Option<u64>, kind: u32, signed: bool) -> BigInt {
    if let Some(d) = d {
        // d-reduced second kind numbers
        if n < d || k < d {
            return if n == k { BigInt::one() } else { BigInt::zero() };
        }
        return stirling2(n - d + 1, k - d + 1);
    }
    if kind == 1 {
        let s = stirling1(n, k);
        if signed && (n - k.min(n)) % 2 == 1 {
            -s
        } else {
            s
        }
    } else {
        stirling2(n, k)
    }
}

/// `multinomial_coefficients(m, n)`: map from exponent tuples `(k_1..k_m)`
/// with `sum k_i = n` to `n!/(k_1!...k_m!)`.
pub fn multinomial_coefficients(m: usize, n: u64) -> Vec<(Vec<u64>, BigInt)> {
    let mut out = Vec::new();
    if m == 0 {
        if n == 0 {
            out.push((vec![], BigInt::one()));
        }
        return out;
    }
    let nf = factorial(n);
    let mut cur = vec![0u64; m];
    fn rec(i: usize, left: u64, cur: &mut Vec<u64>, nf: &BigInt, out: &mut Vec<(Vec<u64>, BigInt)>) {
        let m = cur.len();
        if i == m - 1 {
            cur[i] = left;
            let mut d = BigInt::one();
            for k in cur.iter() {
                d *= factorial(*k);
            }
            out.push((cur.clone(), nf / d));
            return;
        }
        for k in (0..=left).rev() {
            cur[i] = k;
            rec(i + 1, left - k, cur, nf, out);
        }
    }
    rec(0, n, &mut cur, &nf, &mut out);
    out
}

/// Generalised harmonic number `H_n^{(m)} = sum_{k=1}^n k^{-m}` for `m >= 0`.
pub fn harmonic_number(n: u64, m: u32) -> BigRational {
    let mut s = BigRational::zero();
    for k in 1..=n {
        let d = num_traits::pow(BigInt::from(k), m as usize);
        s += BigRational::new(BigInt::one(), d);
    }
    s
}

/// Coefficients (constant first) of the Bernoulli polynomial `B_n(x)`.
pub fn bernoulli_poly_coeffs(n: usize) -> Vec<BigRational> {
    // B_n(x) = sum_k C(n,k) b_k x^(n-k) with b_1 = -1/2
    let mut c = vec![BigRational::zero(); n + 1];
    let mut binom = BigInt::one();
    for k in 0..=n {
        let mut b = bernoulli_number(k);
        if k == 1 {
            b = -b;
        }
        c[n - k] = BigRational::from_integer(binom.clone()) * b;
        binom = binom * BigInt::from(n - k) / BigInt::from(k + 1);
    }
    c
}

/// Coefficients of the Euler polynomial `E_n(x)`.
pub fn euler_poly_coeffs(n: usize) -> Vec<BigRational> {
    // E_n(x) = 2/(n+1) (B_{n+1}(x) - 2^{n+1} B_{n+1}(x/2))
    let b = bernoulli_poly_coeffs(n + 1);
    let two = BigInt::from(2);
    let mut c = vec![BigRational::zero(); n + 2];
    for (j, bj) in b.iter().enumerate() {
        // 2^{n+1} B(x/2): coefficient of x^j is bj 2^{n+1-j}
        let p = num_traits::pow(two.clone(), n + 1 - j);
        c[j] = bj - bj * BigRational::from_integer(p);
    }
    let f = BigRational::new(two, BigInt::from(n as u64 + 1));
    let mut out: Vec<BigRational> = c.into_iter().map(|x| x * &f).collect();
    out.truncate(n + 1);
    out
}

/// Coefficients of the Fibonacci polynomial `F_n(x)` (`F_1 = 1, F_2 = x`).
pub fn fibonacci_poly_coeffs(n: u64) -> Vec<BigInt> {
    if n == 0 {
        return vec![BigInt::zero()];
    }
    let mut a: Vec<BigInt> = vec![BigInt::zero()]; // F_0
    let mut b: Vec<BigInt> = vec![BigInt::one()]; // F_1
    for _ in 1..n {
        // F_{k+1} = x F_k + F_{k-1}
        let mut c = vec![BigInt::zero(); b.len() + 1];
        for (i, v) in b.iter().enumerate() {
            c[i + 1] += v;
        }
        for (i, v) in a.iter().enumerate() {
            c[i] += v;
        }
        a = b;
        b = c;
    }
    b
}

/// Coefficients of the Bell (Touchard) polynomial `B_n(x) = sum S(n,k) x^k`.
pub fn bell_poly_coeffs(n: u64) -> Vec<BigInt> {
    (0..=n).map(|k| stirling2(n, k)).collect()
}

/// `p! / q!`-free helper: exact rational `n` as `(numer, denom)` i64 pair.
pub fn ratio_i64(r: &BigRational) -> Option<(i64, i64)> {
    Some((r.numer().to_i64()?, r.denom().to_i64()?))
}

/// `gcd`-reduced rational from integers.
pub fn q(n: i64, d: i64) -> BigRational {
    BigRational::new(BigInt::from(n), BigInt::from(d))
}

/// `true` if `n` is even.
pub fn is_even(n: &BigInt) -> bool {
    n.is_even()
}

/// absolute value
pub fn babs(n: &BigInt) -> BigInt {
    n.abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequences() {
        assert_eq!(factorial2(7), BigInt::from(105));
        assert_eq!(factorial2(8), BigInt::from(384));
        assert_eq!(subfactorial(5), BigInt::from(44));
        assert_eq!(euler_number(10), BigInt::from(-50521));
        assert_eq!(euler_number(4), BigInt::from(5));
        assert_eq!(catalan_number(30), BigInt::from(3814986502092304i64));
        assert_eq!(fibonacci_number(100).to_string(), "354224848179261915075");
        assert_eq!(fibonacci_number(-5), BigInt::from(5));
        assert_eq!(fibonacci_number(-4), BigInt::from(-3));
        assert_eq!(lucas_number(10), BigInt::from(123));
        assert_eq!(lucas_number(-3), BigInt::from(-4));
        assert_eq!(bell_number(20), BigInt::from(51724158235372i64));
        assert_eq!(partition_number(100), BigInt::from(190569292));
        assert_eq!(stirling(5, 2, None, 2, false), BigInt::from(15));
        assert_eq!(stirling(5, 2, None, 1, false), BigInt::from(50));
        assert_eq!(stirling(5, 2, None, 1, true), BigInt::from(-50));
        assert_eq!(stirling(10, 3, None, 2, false), BigInt::from(9330));
        assert_eq!(stirling(5, 2, Some(2), 2, false), BigInt::from(1));
        assert_eq!(bernoulli_number(50).to_string(), "495057205241079648212477525/66");
        assert_eq!(multinomial_coefficients(3, 2).len(), 6);
    }
}
