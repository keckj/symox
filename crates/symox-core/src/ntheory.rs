//! Number theory: primality, factorisation and arithmetic functions.

use std::collections::BTreeMap;

pub mod bbp;
pub mod congruence;
pub mod discrete;
pub mod ecm;
pub mod factor;
pub mod generate;
pub mod misc;
pub mod mt;
pub mod primetest;
pub mod qs;
pub mod residue;
mod memo;

thread_local! {
    static FACTOR_MEMO: std::cell::RefCell<memo::Memo<BigInt,BTreeMap<BigInt,u32>>> = std::cell::RefCell::new(memo::Memo::new());
    static PRIMEPI_MEMO: std::cell::RefCell<memo::Memo<u64,u64>> = std::cell::RefCell::new(memo::Memo::new());
}

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

// ----------------------------------------------------------------------
// machine word arithmetic
// ----------------------------------------------------------------------

#[inline]
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64 % m;
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

/// Deterministic Miller–Rabin for 64-bit integers.
pub fn is_prime_u64(n: u64) -> bool {
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
    'outer: for a in [2u64, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        let a = a % n;
        if a == 0 {
            continue;
        }
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Pollard–Brent rho: a non-trivial factor of a composite `n`.
fn pollard_brent_u64(n: u64) -> u64 {
    if n % 2 == 0 {
        return 2;
    }
    let mut c = 1u64;
    loop {
        let f = |x: u64| (mulmod(x, x, n) + c) % n;
        let (mut y, mut r, mut q) = (2u64, 1u64, 1u64);
        let (mut g, mut x, mut ys);
        let m = 128;
        loop {
            x = y;
            for _ in 0..r {
                y = f(y);
            }
            let mut k = 0;
            loop {
                ys = y;
                for _ in 0..m.min(r - k) {
                    y = f(y);
                    q = mulmod(q, x.abs_diff(y), n);
                }
                g = gcd_u64(q, n);
                k += m;
                if k >= r || g > 1 {
                    break;
                }
            }
            r *= 2;
            if g > 1 {
                break;
            }
        }
        if g == n {
            loop {
                ys = f(ys);
                g = gcd_u64(x.abs_diff(ys), n);
                if g > 1 {
                    break;
                }
            }
        }
        if g != n {
            return g;
        }
        c += 1;
    }
}

fn factor_u64(n: u64, out: &mut BTreeMap<u64, u32>) {
    if n == 1 {
        return;
    }
    let mut n = n;
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] {
        while n % p == 0 {
            *out.entry(p).or_insert(0) += 1;
            n /= p;
        }
    }
    if n == 1 {
        return;
    }
    let mut stack = vec![n];
    while let Some(m) = stack.pop() {
        if m == 1 {
            continue;
        }
        if is_prime_u64(m) {
            *out.entry(m).or_insert(0) += 1;
            continue;
        }
        let d = pollard_brent_u64(m);
        stack.push(d);
        stack.push(m / d);
    }
}

// ----------------------------------------------------------------------
// big integers
// ----------------------------------------------------------------------

fn big_powmod(a: &BigInt, e: &BigInt, m: &BigInt) -> BigInt {
    a.modpow(e, m)
}

/// Strong probable prime test to base `a`.
fn sprp(n: &BigInt, a: &BigInt) -> bool {
    let one = BigInt::one();
    let nm1 = n - &one;
    let s = nm1.trailing_zeros().unwrap_or(0);
    let d = &nm1 >> s;
    let mut x = big_powmod(a, &d, n);
    if x == one || x == nm1 {
        return true;
    }
    for _ in 1..s {
        x = (&x * &x) % n;
        if x == nm1 {
            return true;
        }
    }
    false
}

pub(crate) fn jacobi_big(a: &BigInt, n: &BigInt) -> i32 {
    // n odd positive
    let mut a = a.mod_floor(n);
    let mut n = n.clone();
    let mut t = 1;
    let three = BigInt::from(3);
    let five = BigInt::from(5);
    let eight = BigInt::from(8);
    while !a.is_zero() {
        while a.is_even() {
            a >>= 1;
            let r = n.mod_floor(&eight);
            if r == three || r == five {
                t = -t;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a.mod_floor(&BigInt::from(4)) == three && n.mod_floor(&BigInt::from(4)) == three {
            t = -t;
        }
        a = a.mod_floor(&n);
    }
    if n.is_one() {
        t
    } else {
        0
    }
}

/// Strong Lucas probable prime test (Selfridge parameters).
fn strong_lucas(n: &BigInt) -> bool {
    // find D with jacobi(D, n) = -1
    let mut d = BigInt::from(5);
    loop {
        let j = jacobi_big(&d, n);
        if j == -1 {
            break;
        }
        if j == 0 && d.abs() != *n {
            return false;
        }
        d = if d.is_positive() { -(d + 2u32) } else { -(d - 2u32) };
        if d.abs() > BigInt::from(1000) {
            // n is likely a perfect square
            if n.sqrt().pow(2) == *n {
                return false;
            }
        }
    }
    let p: BigInt = BigInt::one();
    let q: BigInt = (BigInt::one() - &d) / 4u32;
    let np1: BigInt = n + 1;
    let s = np1.trailing_zeros().unwrap_or(0);
    let dd: BigInt = &np1 >> s;
    // compute U_d, V_d, Q^d by binary method
    let (mut u, mut v, mut qk): (BigInt, BigInt, BigInt) = (BigInt::one(), p.clone(), q.clone());
    let inv2: BigInt = (n + 1u32) / 2u32;
    let bits = dd.bits();
    for i in (0..bits - 1).rev() {
        u = (&u * &v).mod_floor(n);
        v = (&v * &v - &qk * 2u32).mod_floor(n);
        qk = (&qk * &qk).mod_floor(n);
        if dd.bit(i) {
            let u2 = ((&p * &u + &v) * &inv2).mod_floor(n);
            let v2 = ((&d * &u + &p * &v) * &inv2).mod_floor(n);
            u = u2;
            v = v2;
            qk = (&qk * &q).mod_floor(n);
        }
    }
    if u.is_zero() || v.is_zero() {
        return true;
    }
    for _ in 1..s {
        v = (&v * &v - &qk * 2u32).mod_floor(n);
        if v.is_zero() {
            return true;
        }
        qk = (&qk * &qk).mod_floor(n);
    }
    false
}

/// Primality test: deterministic below 2^64, Baillie–PSW above (no known
/// counterexample).
pub fn isprime(n: &BigInt) -> bool {
    if n.sign() != Sign::Plus {
        return false;
    }
    if let Some(m) = n.to_u64() {
        return is_prime_u64(m);
    }
    for p in [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71] {
        if (n % p).is_zero() {
            return false;
        }
    }
    sprp(n, &BigInt::from(2)) && strong_lucas(n)
}

/// Pollard–Brent rho on a composite `n`; gives up (`None`) once the cycle
/// length passes `max_r` (`None` = never).
fn pollard_brent_big(n: &BigInt, max_r: Option<u64>) -> Option<BigInt> {
    let mut c = BigInt::one();
    loop {
        let f = |x: &BigInt| (x * x + &c).mod_floor(n);
        let (mut y, mut r, mut q) = (BigInt::from(2), 1u64, BigInt::one());
        let (mut g, mut x, mut ys);
        let m = 128u64;
        loop {
            x = y.clone();
            for _ in 0..r {
                y = f(&y);
            }
            let mut k = 0;
            loop {
                ys = y.clone();
                for _ in 0..m.min(r - k) {
                    y = f(&y);
                    q = (&q * (&x - &y).abs()).mod_floor(n);
                }
                g = q.gcd(n);
                k += m;
                if k >= r || !g.is_one() {
                    break;
                }
            }
            r *= 2;
            if !g.is_one() {
                break;
            }
            if max_r.is_some_and(|mr| r > mr) {
                return None;
            }
        }
        if g == *n {
            loop {
                ys = f(&ys);
                g = (&x - &ys).abs().gcd(n);
                if !g.is_one() {
                    break;
                }
            }
        }
        if g != *n {
            return Some(g);
        }
        c += 1;
    }
}

/// A non-trivial factor of a composite, non-power `n`: a short rho run
/// for small factors, then ECM with growing bounds (SymPy's schedule).
fn split_big(n: &BigInt) -> BigInt {
    if let Some(d) = pollard_brent_big(n, Some(1 << 16)) {
        return d;
    }
    let seed = BigInt::from(1234);
    let mut b1 = 10_000u64;
    while b1 <= 50_000_000 {
        if let Ok(Some(d)) = ecm::one_factor(n, b1, 100 * b1, 50, &seed) {
            return d;
        }
        b1 *= 5;
    }
    pollard_brent_big(n, None).unwrap()
}

/// Prime factorisation `{p: e}` of `|n|` (`n != 0`).
pub fn factorint(n: &BigInt) -> BTreeMap<BigInt, u32> {
    let key = n.abs();
    if let Some(value) = FACTOR_MEMO.with(|memo|memo.borrow().get(&key)) { return value; }
    let result = factorint_uncached(&key);
    if key.bits() <= 4096 { FACTOR_MEMO.with(|memo|memo.borrow_mut().insert(key,result.clone())); }
    result
}

fn factorint_uncached(n: &BigInt) -> BTreeMap<BigInt,u32> {
    let mut out: BTreeMap<BigInt, u32> = BTreeMap::new();
    let n = n.abs();
    if n.is_zero() || n.is_one() {
        return out;
    }
    if let Some(m) = n.to_u64() {
        let mut o = BTreeMap::new();
        factor_u64(m, &mut o);
        return o.into_iter().map(|(p, e)| (BigInt::from(p), e)).collect();
    }
    let mut n = n;
    // trial division by small primes
    let mut p = 2u64;
    while p < 10_000 {
        let bp = BigInt::from(p);
        while (&n % &bp).is_zero() {
            n /= &bp;
            *out.entry(bp.clone()).or_insert(0) += 1;
        }
        p += if p == 2 { 1 } else { 2 };
    }
    let mut stack = vec![n];
    while let Some(m) = stack.pop() {
        if m.is_one() {
            continue;
        }
        if let Some(u) = m.to_u64() {
            let mut o = BTreeMap::new();
            factor_u64(u, &mut o);
            for (q, e) in o {
                *out.entry(BigInt::from(q)).or_insert(0) += e;
            }
            continue;
        }
        if isprime(&m) {
            *out.entry(m).or_insert(0) += 1;
            continue;
        }
        // perfect powers make rho slow
        let mut split = false;
        for k in 2..=(m.bits() as u32).min(64) {
            let r = m.nth_root(k);
            if r.pow(k) == m {
                for _ in 0..k {
                    stack.push(r.clone());
                }
                split = true;
                break;
            }
        }
        if split {
            continue;
        }
        let d = split_big(&m);
        stack.push(&m / &d);
        stack.push(d);
    }
    out
}

/// Primes in `[a, b)` (segmented sieve).
pub fn primerange(a: u64, b: u64) -> Vec<u64> {
    if b <= 2 || b <= a {
        return vec![];
    }
    let lim = (b as f64).sqrt() as u64 + 1;
    let mut small = vec![true; (lim + 1) as usize];
    let mut base = Vec::new();
    for i in 2..=lim {
        if small[i as usize] {
            base.push(i);
            let mut j = i * i;
            while j <= lim {
                small[j as usize] = false;
                j += i;
            }
        }
    }
    let a = a.max(2);
    let mut out = Vec::new();
    let seg = 1u64 << 18;
    let mut lo = a;
    while lo < b {
        let hi = (lo + seg).min(b);
        let mut mark = vec![true; (hi - lo) as usize];
        for &p in &base {
            if p * p >= hi {
                break;
            }
            let start = ((lo + p - 1) / p * p).max(p * p);
            let mut j = start;
            while j < hi {
                mark[(j - lo) as usize] = false;
                j += p;
            }
        }
        for (i, m) in mark.iter().enumerate() {
            if *m {
                out.push(lo + i as u64);
            }
        }
        lo = hi;
    }
    out
}

/// Smallest prime > n.
pub fn nextprime(n: &BigInt) -> BigInt {
    let mut m: BigInt = if n < &BigInt::from(2) { BigInt::from(2) } else { n + 1 };
    if m == BigInt::from(2) {
        return m;
    }
    if m.is_even() {
        m += 1;
    }
    while !isprime(&m) {
        m += 2;
    }
    m
}

/// Largest prime < n (None if n <= 2).
pub fn prevprime(n: &BigInt) -> Option<BigInt> {
    if n <= &BigInt::from(2) {
        return None;
    }
    if *n == BigInt::from(3) {
        return Some(BigInt::from(2));
    }
    let mut m: BigInt = n - 1;
    if m.is_even() {
        m -= 1;
    }
    while !isprime(&m) {
        m -= 2;
    }
    Some(m)
}

/// The n-th prime (1-based).
pub fn prime(n: u64) -> u64 {
    assert!(n >= 1);
    if n < 6 {
        return [2, 3, 5, 7, 11][(n - 1) as usize];
    }
    let nf = n as f64;
    let ub = (nf * (nf.ln() + nf.ln().ln())) as u64 + 10;
    primerange(2, ub)[(n - 1) as usize]
}

/// Number of primes <= n (Lucy Hedgehog's O(n^(3/4)) algorithm).
pub fn primepi(n: u64) -> u64 {
    if let Some(value) = PRIMEPI_MEMO.with(|memo|memo.borrow().get(&n)) { return value; }
    let value = primepi_uncached(n);
    PRIMEPI_MEMO.with(|memo|memo.borrow_mut().insert(n,value));
    value
}

fn primepi_uncached(n: u64) -> u64 {
    if n < 2 {
        return 0;
    }
    let r = (n as f64).sqrt() as u64;
    let r = if (r + 1) * (r + 1) <= n { r + 1 } else if r * r > n { r - 1 } else { r };
    // values n // i for i in 1..=r, then the small values r'..=1
    let mut v: Vec<u64> = (1..=r).map(|i| n / i).collect();
    let last = *v.last().unwrap();
    v.extend((1..last).rev());
    let idx = |x: u64| -> usize {
        if x <= last.saturating_sub(1).max(r) && x < last {
            // small values are stored at the end, in decreasing order
            v.len() - x as usize
        } else {
            (n / x - 1) as usize
        }
    };
    let mut s: Vec<u64> = v.iter().map(|&x| x - 1).collect();
    for p in 2..=r {
        let ip = idx(p);
        if s[ip] > s[idx(p - 1)] {
            let sp = s[idx(p - 1)];
            let p2 = p * p;
            for k in 0..v.len() {
                let x = v[k];
                if x < p2 {
                    break;
                }
                let j = idx(x / p);
                s[k] -= s[j] - sp;
            }
        }
    }
    s[0]
}

/// Euler's totient.
pub fn totient(n: &BigInt) -> BigInt {
    let mut r = n.abs();
    for p in factorint(n).keys() {
        r = r / p * (p - 1);
    }
    r
}

/// Sorted positive divisors.
pub fn divisors(n: &BigInt) -> Vec<BigInt> {
    let mut ds = vec![BigInt::one()];
    for (p, e) in factorint(n) {
        let mut next = Vec::with_capacity(ds.len() * (e as usize + 1));
        for d in &ds {
            let mut pk = BigInt::one();
            for _ in 0..=e {
                next.push(d * &pk);
                pk *= &p;
            }
        }
        ds = next;
    }
    ds.sort();
    ds
}

/// Number of divisors.
pub fn divisor_count(n: &BigInt) -> BigInt {
    factorint(n).values().fold(BigInt::one(), |acc, e| acc * (e + 1))
}

/// Sum of the k-th powers of the divisors.
pub fn divisor_sigma(n: &BigInt, k: u32) -> BigInt {
    let mut r = BigInt::one();
    for (p, e) in factorint(n) {
        let pk = p.pow(k);
        let mut s = BigInt::zero();
        let mut t = BigInt::one();
        for _ in 0..=e {
            s += &t;
            t *= &pk;
        }
        r *= s;
    }
    r
}

/// Modular inverse of `a` mod `m`.
pub fn mod_inverse(a: &BigInt, m: &BigInt) -> Option<BigInt> {
    let e = a.extended_gcd(m);
    if !e.gcd.is_one() && e.gcd != -BigInt::one() {
        return None;
    }
    Some(e.x.mod_floor(m))
}

/// Jacobi symbol (a/n) for odd positive n.
pub fn jacobi_symbol(a: &BigInt, n: &BigInt) -> Option<i32> {
    if n.is_even() || !n.is_positive() {
        return None;
    }
    Some(jacobi_big(a, n))
}

/// Legendre symbol (a/p) for an odd prime p.
pub fn legendre_symbol(a: &BigInt, p: &BigInt) -> Option<i32> {
    if !isprime(p) || *p == BigInt::from(2) {
        return None;
    }
    Some(jacobi_big(a, p))
}

/// Möbius function.
pub fn mobius(n: &BigInt) -> i32 {
    let f = factorint(n);
    if f.values().any(|e| *e > 1) {
        0
    } else if f.len() % 2 == 0 {
        1
    } else {
        -1
    }
}

/// `(root, exact)` with root = floor(n^(1/k)).
pub fn integer_nthroot(n: &BigInt, k: u32) -> (BigInt, bool) {
    let r = n.nth_root(k);
    let exact = r.pow(k) == *n;
    (r, exact)
}

/// `n = b**e` with the largest possible `e`, if `n` is a perfect power.
pub fn perfect_power(n: &BigInt) -> Option<(BigInt, u32)> {
    if n <= &BigInt::one() {
        return None;
    }
    let f = factorint(n);
    let g = f.values().fold(0u32, |g, e| g.gcd(e));
    if g < 2 {
        return None;
    }
    let b = f.iter().fold(BigInt::one(), |acc, (p, e)| acc * p.pow(e / g));
    Some((b, g))
}

/// Multiplicative order of `a` modulo `n`.
pub fn n_order(a: &BigInt, n: &BigInt) -> Option<BigInt> {
    if !a.gcd(n).is_one() {
        return None;
    }
    let phi = totient(n);
    let mut order = phi.clone();
    for (p, _) in factorint(&phi) {
        while (&order % &p).is_zero() && a.modpow(&(&order / &p), n).is_one() {
            order /= &p;
        }
    }
    Some(order)
}

/// Smallest primitive root modulo `p` (prime).
pub fn primitive_root(p: &BigInt) -> Option<BigInt> {
    if !isprime(p) {
        return None;
    }
    if *p == BigInt::from(2) {
        return Some(BigInt::one());
    }
    let phi: BigInt = p - 1;
    let fs: Vec<BigInt> = factorint(&phi).into_keys().collect();
    let mut g = BigInt::from(2);
    loop {
        if fs.iter().all(|q| !g.modpow(&(&phi / q), p).is_one()) {
            return Some(g);
        }
        g += 1;
    }
}

/// n!
pub fn factorial(n: u64) -> BigInt {
    // binary splitting product
    fn prod(a: u64, b: u64) -> BigInt {
        if b <= a + 8 {
            let mut r = BigInt::one();
            for i in a..=b {
                r *= i;
            }
            return r;
        }
        let m = (a + b) / 2;
        prod(a, m) * prod(m + 1, b)
    }
    if n < 2 {
        BigInt::one()
    } else {
        prod(2, n)
    }
}

/// Binomial coefficient C(n, k) for integers (n may be negative).
pub fn binomial(n: &BigInt, k: &BigInt) -> BigInt {
    if k.is_negative() {
        return BigInt::zero();
    }
    if n.is_negative() {
        // C(n, k) = (-1)^k C(k - n - 1, k)
        let r = binomial(&(k - n - 1), k);
        return if k.is_odd() { -r } else { r };
    }
    if k > n {
        return BigInt::zero();
    }
    let k = if (n - k) < *k { n - k } else { k.clone() };
    let mut r = BigInt::one();
    if let Some(kk) = k.to_u64() {
        for i in 0..kk {
            r = r * (n - i) / (i + 1);
        }
    } else {
        // An unrepresentable loop bound must never be mistaken for zero
        // iterations (which returned the false exact answer 1). Symbolic
        // callers retain such enormous outputs unevaluated before entering
        // this low-level exact integer routine.
        let mut i = BigInt::zero();
        while i < k {
            r = r * (n - &i) / (&i + 1);
            i += 1;
        }
    }
    r
}

/// Continued fraction expansion of a rational p/q.
pub fn continued_fraction(p: &BigInt, q: &BigInt) -> Vec<BigInt> {
    let (mut p, mut q) = (p.clone(), q.clone());
    let mut out = Vec::new();
    while !q.is_zero() {
        let (a, r) = p.div_mod_floor(&q);
        out.push(a);
        p = q;
        q = r;
    }
    out
}

/// Chinese remainder theorem: x ≡ r_i mod m_i (moduli need not be coprime).
pub fn crt(rs: &[BigInt], ms: &[BigInt]) -> Option<(BigInt, BigInt)> {
    let mut x = BigInt::zero();
    let mut m = BigInt::one();
    for (r, mi) in rs.iter().zip(ms) {
        let g = m.gcd(mi);
        let diff = r - &x;
        if !(&diff % &g).is_zero() {
            return None;
        }
        let l = &m / &g * mi;
        let inv = mod_inverse(&(&m / &g), &(mi / &g)).unwrap_or_else(BigInt::zero);
        let t = ((&diff / &g) * inv).mod_floor(&(mi / &g));
        x = (&x + &m * t).mod_floor(&l);
        m = l;
    }
    Some((x, m))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(n: &str) -> BigInt {
        n.parse().unwrap()
    }

    #[test]
    fn primes() {
        assert!(is_prime_u64(2) && is_prime_u64(97) && !is_prime_u64(561));
        assert!(is_prime_u64(18446744073709551557));
        assert!(isprime(&b("170141183460469231731687303715884105727"))); // 2^127-1
        assert!(!isprime(&b("170141183460469231731687303715884105729")));
        assert_eq!(primerange(1, 30), vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
        assert_eq!(prime(1000), 7919);
        assert_eq!(primepi(1_000_000), 78498);
        assert_eq!(primepi(10_000_000_000), 455052511);
        for n in [0u64, 1, 2, 3, 10, 11, 100, 101, 1000, 7919, 65536] {
            assert_eq!(primepi(n), primerange(2, n + 1).len() as u64, "n={n}");
        }
        assert_eq!(nextprime(&b("100")), b("101"));
    }

    #[test]
    fn factoring() {
        let f = factorint(&b("600851475143"));
        assert_eq!(f.into_iter().collect::<Vec<_>>(), vec![(b("71"), 1), (b("839"), 1), (b("1471"), 1), (b("6857"), 1)]);
        let n = b("1000000016000000063"); // 1000000007 * 1000000009
        let f = factorint(&n);
        assert_eq!(f.len(), 2);
        let n = b("18446744073709551617"); // 2^64 + 1
        let f = factorint(&n);
        assert_eq!(f.keys().cloned().collect::<Vec<_>>(), vec![b("274177"), b("67280421310721")]);
        assert_eq!(totient(&b("36")), b("12"));
        assert_eq!(divisors(&b("12")), vec![b("1"), b("2"), b("3"), b("4"), b("6"), b("12")]);
        assert_eq!(mod_inverse(&b("3"), &b("11")), Some(b("4")));
        assert_eq!(factorial(20), b("2432902008176640000"));
        assert_eq!(binomial(&b("10"), &b("3")), b("120"));
        assert_eq!(primitive_root(&b("23")), Some(b("5")));
        assert_eq!(n_order(&b("2"), &b("7")), Some(b("3")));
        assert_eq!(crt(&[b("2"), b("3"), b("2")], &[b("3"), b("5"), b("7")]), Some((b("23"), b("105"))));
    }
}
