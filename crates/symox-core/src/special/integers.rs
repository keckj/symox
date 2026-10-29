//! Exact integer and rational sequences used by the special functions.

use std::sync::{Mutex, OnceLock};

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// Largest exact integer, in bits (32 MiB), that the automatic evaluation
/// of the combinatorial sequences materialises. Above it a call stays
/// unevaluated: the reference (pure Python integers) already needs minutes
/// at a quarter of this size (`factorial(4*10**6)`, 8.2e7 bits, 40 s) and
/// does not finish `factorial(10**7)` (2.2e8 bits) in 150 s.
pub const MAX_EXACT_BITS: f64 = 268_435_456.0;

/// `true` when an exact result of about `bits` bits is materialised.
pub fn fits_exact(bits: f64) -> bool {
    bits.is_finite() && bits <= MAX_EXACT_BITS
}

/// Upper bound of `log2(n!)` (Stirling series truncated after `1/(12n)`).
pub fn log2_factorial(n: f64) -> f64 {
    if n < 2.0 {
        return 0.0;
    }
    (n * n.ln() - n + 0.5 * (2.0 * std::f64::consts::PI * n).ln() + 1.0 / (12.0 * n)) / std::f64::consts::LN_2
}

/// Upper bound of `log2(binomial(n, k))` for `0 <= k <= n`.
pub fn log2_binomial(n: f64, k: f64) -> f64 {
    let m = k.min(n - k);
    if m <= 0.0 {
        return 0.0;
    }
    if n < 1e12 {
        // Stirling difference; the slack covers the truncation errors
        log2_factorial(n) - log2_factorial(m) - log2_factorial(n - m) + 2.0
    } else {
        // binomial(n, m) <= (e n / m)^m
        m * ((n / m).log2() + std::f64::consts::LOG2_E)
    }
}

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
    // !k = k*!(k-1) + (-1)^k: compose the affine maps x -> k*x + (-1)^k
    // over [a, b] by binary splitting, as (prod k, offset).
    fn affine(a: u64, b: u64) -> (BigInt, BigInt) {
        if b - a < 16 {
            let (mut m, mut c) = (BigInt::one(), BigInt::zero());
            for k in a..=b {
                m *= k;
                c *= k;
                if k % 2 == 0 {
                    c += 1u32;
                } else {
                    c -= 1u32;
                }
            }
            return (m, c);
        }
        let mid = a + (b - a) / 2;
        let (ml, cl) = affine(a, mid);
        let (mr, cr) = affine(mid + 1, b);
        let c = &mr * cl + cr;
        (mr * ml, c)
    }
    if n == 0 {
        return BigInt::one();
    }
    let (m, c) = affine(1, n);
    m + c
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
    if n >= 120 && n > 2 * c.len() {
        // isolated large index: the zeta formula beats rebuilding the table
        return crate::summation::bernoulli_rational(n);
    }
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

/// `E_n` for a large even `n = 2m` from the Dirichlet beta function, the
/// analogue of the zeta method of the Bernoulli numbers:
/// `E_2m = (-1)^m 4^(m+1) (2m)! beta(2m+1) / pi^(2m+1)` is an integer, the
/// rounding of a value computed with enough bits; `None` if the value is
/// not certified within 1/8 of an integer after the retries.
fn euler_beta(n: usize) -> Option<BigInt> {
    use crate::bigfloat::{pi, BigFloat};
    let m = n / 2;
    let nf = n as f64;
    let s = n as i64 + 1;
    let bits = 2.0 * (m as f64 + 1.0) + log2_factorial(nf) - (nf + 1.0) * std::f64::consts::PI.log2();
    let bitlen = |v: u64| 64 - v.leading_zeros();
    let fact = factorial(n as u64);
    let mut guard = 64 + 2 * bitlen(n as u64);
    for _ in 0..4 {
        let prec = bits.max(0.0) as u32 + guard;
        // beta(s) = sum_k (-1)^k (2k+1)^-s, truncated below 2^-(prec+8)
        let last = (2f64.powf((prec as f64 + 8.0) / (s as f64))).ceil() as u64 + 2;
        let wp = prec + bitlen(last) + 8;
        let mut beta = BigFloat::one();
        let mut j = 3u64;
        let mut negative = true;
        while j <= last {
            let drop = (s as f64 * (j as f64).log2()).floor();
            let pj = (wp as f64 + 8.0 - drop).max(32.0) as u32;
            let term = BigFloat::from_i64(j as i64).powi(-s, pj)?;
            beta = if negative { beta.sub(&term, wp) } else { beta.add(&term, wp) };
            negative = !negative;
            j += 2;
        }
        let pi_s = pi(wp).powi(s, wp)?;
        let scale = BigFloat::from_int(&(&fact << (2 * (m + 1))));
        let value = scale.mul(&beta, wp).div(&pi_s, wp)?;
        let num = value.to_int_round();
        let gap = value.sub(&BigFloat::from_int(&num), wp).abs();
        if gap.cmp(&BigFloat::from_parts(BigInt::one(), -3i64)) == std::cmp::Ordering::Less {
            return Some(if m % 2 == 1 { -num } else { num });
        }
        guard *= 2;
    }
    None
}

/// Euler (secant) number `E_n` (`E_0 = 1, E_2 = -1, E_4 = 5`, odd ones 0).
pub fn euler_number(n: usize) -> BigInt {
    if n % 2 == 1 {
        return BigInt::zero();
    }
    if n >= 120 {
        if let Some(e) = euler_beta(n) {
            return e;
        }
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
    let b = crate::ntheory::binomial_u64(2 * n, n);
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

/// Tribonacci number `T_n` (`T_0 = 0, T_1 = T_2 = 1`): `x^n` is reduced
/// modulo `x^3 - x^2 - x - 1` by binary powering, so `T_n = a + b` for
/// `x^n = a x^2 + b x + c`.
pub fn tribonacci_number(n: u64) -> BigInt {
    let (mut a, mut b, mut c) = (BigInt::zero(), BigInt::zero(), BigInt::one());
    for bit in (0..64 - n.leading_zeros()).rev() {
        // square, with x^3 = x^2 + x + 1 and x^4 = 2x^2 + 2x + 1
        let aa = &a * &a;
        let ab2 = (&a * &b) << 1usize;
        let ac2 = (&a * &c) << 1usize;
        let bb = &b * &b;
        let bc2 = (&b * &c) << 1usize;
        let cc = &c * &c;
        let aa2 = &aa << 1usize;
        let na = &aa2 + &ab2 + &bb + &ac2;
        let nb = &aa2 + &ab2 + &bc2;
        let nc = aa + ab2 + cc;
        (a, b, c) = (na, nb, nc);
        if (n >> bit) & 1 == 1 {
            // multiply by x
            let na = &a + &b;
            let nb = &a + &c;
            (a, b, c) = (na, nb, a);
        }
    }
    a + b
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
        s = crate::bigint_gcd::ratio_add(&s, &BigRational::new_raw(BigInt::one(), d));
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
        let seidel = |n: usize| {
            let mut row = vec![BigInt::one()];
            for m in 1..=n {
                let mut next = vec![BigInt::zero()];
                for j in 0..m {
                    let v = &next[j] + &row[m - 1 - j];
                    next.push(v);
                }
                row = next;
            }
            let a = row.last().unwrap().clone();
            if (n / 2) % 2 == 1 { -a } else { a }
        };
        for n in [120usize, 122, 124, 250, 400, 402] {
            assert_eq!(euler_beta(n), Some(seidel(n)), "E_{n}");
        }
        assert_eq!(catalan_number(30), BigInt::from(3814986502092304i64));
        assert_eq!(fibonacci_number(100).to_string(), "354224848179261915075");
        assert_eq!(fibonacci_number(-5), BigInt::from(5));
        assert_eq!(fibonacci_number(-4), BigInt::from(-3));
        assert_eq!(lucas_number(10), BigInt::from(123));
        assert_eq!(lucas_number(-3), BigInt::from(-4));
        let mut t = (BigInt::zero(), BigInt::one(), BigInt::one());
        for k in 0..200u64 {
            assert_eq!(tribonacci_number(k), t.0);
            t = (t.1.clone(), t.2.clone(), &t.0 + &t.1 + &t.2);
        }
        let (mut a, mut b) = (BigInt::one(), BigInt::zero());
        assert_eq!(subfactorial(0), BigInt::one());
        for k in 2..300u64 {
            let c = (&a + &b) * BigInt::from(k - 1);
            a = b;
            b = c;
            assert_eq!(subfactorial(k), b);
        }
        let mut c = BigInt::one();
        for n in 0..300u64 {
            assert_eq!(catalan_number(n), c);
            c = c * BigInt::from(2 * (2 * n + 1)) / BigInt::from(n + 2);
        }
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
