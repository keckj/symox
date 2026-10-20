//! Factorisation helpers ported from SymPy's `ntheory.factor_`:
//! multiplicities, Pollard's rho and p-1 methods (with SymPy's seeds),
//! divisor enumeration orders, perfect powers and the abundance family.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::mt::PyRandom;
use super::{factorint, isprime, primerange};

/// Exponents of the known Mersenne primes `2**p - 1`.
pub const MERSENNE_PRIME_EXPONENTS: [u64; 52] = [
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
    9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
    1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
    32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933, 136279841,
];

/// SymPy's `factorint` dictionary as an ordered list: `{0: 1}` for zero,
/// `{}` for one, prime factors in increasing order then `-1` for negative
/// numbers.
pub fn factor_list(n: &BigInt) -> Vec<(BigInt, u32)> {
    if n.is_zero() {
        return vec![(BigInt::zero(), 1)];
    }
    let mut v: Vec<(BigInt, u32)> = factorint(n).into_iter().collect();
    if n.is_negative() {
        v.push((-BigInt::one(), 1));
    }
    v
}

/// `(n / f**m, m)` with `m` maximal (SymPy's `remove`); `None` for `f < 2`.
pub fn remove(n: &BigInt, f: &BigInt) -> Option<(BigInt, u64)> {
    if f < &BigInt::from(2) {
        return None;
    }
    if n.is_zero() {
        return Some((BigInt::zero(), 0));
    }
    let mut n = n.clone();
    let mut m = 0u64;
    loop {
        let (q, r) = n.div_rem(f);
        if !r.is_zero() {
            return Some((n, m));
        }
        n = q;
        m += 1;
    }
}

/// Number of trailing zero bits of `|n|` (0 for 0).
pub fn trailing(n: &BigInt) -> u64 {
    n.trailing_zeros().unwrap_or(0)
}

/// Multiplicity of `p` in `n` for integers (`Err` carries SymPy's message).
pub fn multiplicity(p: &BigInt, n: &BigInt) -> Result<u64, String> {
    if n.is_zero() {
        return Err(format!("no such integer exists: multiplicity of {n} is not-defined"));
    }
    remove(n, p).map(|r| r.1).ok_or_else(|| "factor must be > 1".to_string())
}

/// Multiplicity of `p` in `n!`.
pub fn multiplicity_in_factorial(p: &BigInt, n: &BigInt) -> Result<BigInt, String> {
    if !p.is_positive() {
        return Err(format!("expecting positive integer got {p}"));
    }
    if n.is_negative() {
        return Err(format!("expecting non-negative integer got {n}"));
    }
    // f[v] = max prime with exponent v
    let mut f: Vec<(u32, BigInt)> = Vec::new();
    for (k, v) in factorint(p) {
        match f.iter_mut().find(|(e, _)| *e == v) {
            Some(slot) => {
                if k > slot.1 {
                    slot.1 = k
                }
            }
            None => f.push((v, k)),
        }
    }
    let mut best: Option<BigInt> = None;
    for (v, k) in f {
        // (n + k - digit_sum(n, k)) // (k - 1) // v
        let mut s = BigInt::zero();
        let mut m = n.clone();
        while !m.is_zero() {
            let (q, r) = m.div_rem(&k);
            s += r;
            m = q;
        }
        let val = ((n - s).div_floor(&(&k - 1))).div_floor(&BigInt::from(v));
        best = Some(match best {
            Some(b) if b <= val => b,
            _ => val,
        });
    }
    // p == 1: SymPy takes min() of an empty sequence
    best.ok_or_else(|| "min() arg is an empty sequence".to_string())
}

fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    a.gcd(b)
}

/// Pollard's rho method with SymPy's defaults and pseudo-random retries.
///
/// `f` is an optional user iteration function; `Err` values from it are
/// propagated unchanged.
#[allow(clippy::too_many_arguments)]
pub fn pollard_rho<E>(
    n: &BigInt,
    s: &BigInt,
    a: &BigInt,
    retries: i64,
    seed: &BigInt,
    max_steps: Option<&BigInt>,
    mut f: Option<&mut dyn FnMut(&BigInt) -> Result<BigInt, E>>,
) -> Result<Option<BigInt>, E> {
    let mut rng = PyRandom::new(&(seed + retries));
    let mut v = s.clone();
    let mut a = a.clone();
    let max_steps = max_steps.filter(|m| !m.is_zero());
    for _ in 0..(retries + 1).max(0) {
        let mut u = v.clone();
        let mut j = BigInt::zero();
        loop {
            if let Some(ms) = max_steps {
                if &j > ms {
                    break;
                }
            }
            j += 1;
            let step = |x: &BigInt, f: &mut Option<&mut dyn FnMut(&BigInt) -> Result<BigInt, E>>| -> Result<BigInt, E> {
                match f {
                    Some(g) => g(x),
                    None => Ok((x.modpow(&BigInt::from(2), n) + &a).mod_floor(n)),
                }
            };
            u = step(&u, &mut f)?;
            let v1 = step(&v, &mut f)?;
            v = step(&v1, &mut f)?;
            let g = gcd(&(&u - &v), n);
            if g.is_one() {
                continue;
            }
            if &g == n {
                break;
            }
            return Ok(Some(g));
        }
        v = rng.randint(&BigInt::zero(), &(n - 1));
        a = rng.randint(&BigInt::one(), &(n - 3));
        f = None;
    }
    Ok(None)
}

/// Pollard's p-1 method with SymPy's defaults and pseudo-random retries.
pub fn pollard_pm1(n: &BigInt, b: u64, a: &BigInt, retries: i64, seed: &BigInt) -> Option<BigInt> {
    let mut rng = PyRandom::new(&(seed + b));
    let mut a = a.clone();
    let primes = primerange(2, b + 1);
    let lb = (b as f64).ln();
    for _ in 0..(retries + 1).max(0) {
        let mut am = a.clone();
        for &p in &primes {
            let e = (lb / (p as f64).ln()) as u32;
            am = am.modpow(&BigInt::from(p).pow(e), n);
        }
        let g = gcd(&(&am - 1), n);
        if g > BigInt::one() && &g < n {
            return Some(g);
        }
        a = rng.randint(&BigInt::from(2), &(n - 2));
    }
    None
}

/// Prime factors of `|n|` in increasing order.
pub fn primefactors(n: &BigInt) -> Vec<BigInt> {
    if n.abs() <= BigInt::one() {
        return vec![];
    }
    factorint(n).into_keys().collect()
}

/// Divisors of `|n|` in SymPy's generator order (`proper` excludes `|n|`).
pub fn divisors_unsorted(n: &BigInt, proper: bool) -> Vec<BigInt> {
    let n = n.abs();
    if n <= BigInt::one() {
        return if !proper && n.is_one() { vec![BigInt::one()] } else { vec![] };
    }
    let fs: Vec<(BigInt, u32)> = factorint(&n).into_iter().collect();
    let mut l = vec![BigInt::one()];
    for (p, e) in fs.iter().rev() {
        let mut pows = vec![BigInt::one()];
        for _ in 0..*e {
            let x = pows.last().unwrap() * p;
            pows.push(x);
        }
        let mut next = Vec::with_capacity(l.len() * pows.len());
        for q in &l {
            for pw in &pows {
                next.push(pw * q);
            }
        }
        l = next;
    }
    if proper {
        l.retain(|d| d != &n);
    }
    l
}

/// SymPy's `_divisor_sigma` evaluated on the raw `factorint` dictionary
/// (including the `-1` and `0` keys).
pub fn divisor_sigma_raw(n: &BigInt, k: u32) -> BigInt {
    let fs = factor_list(n);
    if k == 0 {
        return fs.iter().fold(BigInt::one(), |acc, (_, e)| acc * (e + 1));
    }
    let mut r = BigInt::one();
    for (p, e) in fs {
        let num: BigInt = p.pow(k * (e + 1)) - 1;
        let den: BigInt = p.pow(k) - 1;
        r *= num.div_floor(&den);
    }
    r
}

/// `divisor_sigma(n) - 2*n`.
pub fn abundance(n: &BigInt) -> BigInt {
    divisor_sigma_raw(n, 1) - n * 2
}

/// Lucas-Lehmer test for `2**p - 1` (`p` an odd prime).
fn lucas_lehmer(p: u64) -> bool {
    let m: BigInt = (BigInt::one() << p) - 1;
    let mut s = BigInt::from(4);
    for _ in 0..p - 2 {
        s = (&s * &s - 2u32).mod_floor(&m);
    }
    s.is_zero()
}

/// Is `n` a Mersenne prime?
pub fn is_mersenne_prime(n: &BigInt) -> bool {
    if n < &BigInt::one() {
        return false;
    }
    let n1: BigInt = n + 1;
    if !(n & &n1).is_zero() {
        return false;
    }
    let p = n.bits();
    if MERSENNE_PRIME_EXPONENTS.contains(&p) {
        return true;
    }
    if p < 65_000_000 || !isprime(&BigInt::from(p)) {
        return false;
    }
    lucas_lehmer(p)
}

/// Is `n` a perfect number?  `Err` is SymPy's (theoretical) odd perfect
/// number report.
pub fn is_perfect(n: &BigInt) -> Result<bool, String> {
    if n < &BigInt::one() {
        return Ok(false);
    }
    if n.is_even() {
        let m = (n.bits() + 1) >> 1;
        let v: BigInt = (BigInt::one() << (m - 1)) * ((BigInt::one() << m) - 1);
        if &v != n {
            return Ok(false);
        }
        return Ok(MERSENNE_PRIME_EXPONENTS.contains(&m) || is_mersenne_prime(&((BigInt::one() << m) - 1)));
    }
    if n.bits() < 6644 && *n < BigInt::from(10).pow(2000) {
        return Ok(false);
    }
    if (n % 105u32).is_zero() {
        return Ok(false);
    }
    if [(12u32, 1u32), (468, 117), (324, 81)].iter().all(|(m, r)| n % m != BigInt::from(*r)) {
        return Ok(false);
    }
    if abundance(n).is_zero() {
        return Err(format!("odd perfect number found: {n}"));
    }
    Ok(false)
}

/// `is_carmichael`: `Some(b)` is a bool result, `None` is SymPy's `0`
/// (returned for even `n >= 561`).
pub fn is_carmichael(n: &BigInt) -> Option<bool> {
    if n < &BigInt::from(561) {
        return Some(false);
    }
    if n.is_even() {
        return None;
    }
    if isprime(n) {
        return Some(false);
    }
    let n1: BigInt = n - 1;
    Some(factorint(n).iter().all(|(p, e)| *e == 1 && (&n1 % (p - 1u32)).is_zero()))
}

/// Exact integer `k`-th root (`floor`) of `n >= 0` and exactness.
pub fn iroot(n: &BigInt, k: u32) -> (BigInt, bool) {
    let r = n.nth_root(k);
    let exact = r.pow(k) == *n;
    (r, exact)
}

fn smallest_prime_factor(e: u64) -> u64 {
    let mut p = 2;
    while p * p <= e {
        if e % p == 0 {
            return p;
        }
        p += 1;
    }
    e
}

/// Largest power decomposition `n = b**e` (`e >= 2` maximal) for `n > 0`.
pub(crate) fn perfect_power_big(n: &BigInt) -> Option<(BigInt, u64)> {
    thread_local! {
        static MEMO: std::cell::RefCell<super::memo::Memo<BigInt,Option<(BigInt,u64)>>> = std::cell::RefCell::new(super::memo::Memo::new());
    }
    if let Some(value) = MEMO.with(|memo|memo.borrow().get(n)) { return value; }
    let value = perfect_power_big_uncached(n);
    if n.bits() <= 4096 { MEMO.with(|memo|memo.borrow_mut().insert(n.clone(),value.clone())); }
    value
}

fn perfect_power_big_uncached(n: &BigInt) -> Option<(BigInt,u64)> {
    if n <= &BigInt::from(3) {
        return None;
    }
    let mut b = n.clone();
    let mut e = 1u64;
    let mut p = 2u32;
    if is_mersenne_nonpower(&b) { return None; }
    while b.bits() as u32 >= p {
        let (r, exact) = iroot(&b, p);
        if exact && r > BigInt::one() {
            b = r;
            e *= p as u64;
            if is_mersenne_nonpower(&b) { break; }
            continue;
        }
        p += 1;
        while !super::is_prime_u64(p as u64) {
            p += 1;
        }
    }
    if e > 1 {
        Some((b, e))
    } else {
        None
    }
}

/// A positive 2**k-1 (k>=2) is never a nontrivial perfect power. Even
/// powers of odd numbers are 1 modulo 8, rather than 7 (k=2 gives 3).
/// For an odd prime q, a**q+1 has the odd factor
/// a**(q-1)-a**(q-2)+...-a+1 > 1, so cannot be a power of two.
/// The caller supplies n>1. This is an exact bit test, not a probable-prime
/// shortcut; composite Mersenne numbers satisfy the same argument.
fn is_mersenne_nonpower(n:&BigInt)->bool{
    n.magnitude().trailing_ones()==n.bits()
}

/// SymPy's `perfect_power(n, candidates, big)` for integers.  Returns
/// `None` where SymPy returns `False`.
pub fn perfect_power(n: &BigInt, candidates: Option<&[BigInt]>, big: bool) -> Option<(BigInt, u64)> {
    if n.is_negative() {
        let m = -n;
        match candidates {
            None => {
                let (b, e) = perfect_power(&m, None, true)?;
                let e2 = e & e.wrapping_neg();
                let (b, e) = (b.pow(e2 as u32), e / e2);
                if e <= 1 {
                    return None;
                }
                if big || isprime(&BigInt::from(e)) {
                    return Some((-b, e));
                }
                let p = smallest_odd_prime_factor(e);
                Some((-(b.pow((e / p) as u32)), p))
            }
            Some(c) => {
                let odd: Vec<BigInt> = c.iter().filter(|i| i.is_odd()).cloned().collect();
                if odd.is_empty() {
                    return None;
                }
                let (b, e) = perfect_power(&m, Some(&odd), big)?;
                Some((-b, e))
            }
        }
    } else {
        match candidates {
            None => {
                let (b, e) = perfect_power_big(n)?;
                if big {
                    return Some((b, e));
                }
                let p = smallest_prime_factor(e);
                Some((b.pow((e / p) as u32), p))
            }
            Some(c) => {
                if n <= &BigInt::from(3) {
                    return None;
                }
                // floor(log2(n)) + 2; a root of exponent bits(n) would be < 2
                let max_possible = n.bits() as i64 + 1;
                let not_square = matches!((n % 10u32).to_u32().unwrap_or(0), 2 | 3 | 7 | 8);
                let min_possible = 2 + not_square as i64;
                let mut cs: Vec<i64> = c
                    .iter()
                    .filter_map(|i| i.to_i64())
                    .filter(|i| min_possible <= *i && *i < max_possible)
                    .collect();
                cs.sort();
                cs.dedup();
                if n.is_even() {
                    let e = trailing(n) as i64;
                    cs.retain(|i| e % i == 0);
                }
                if big {
                    cs.reverse();
                }
                for e in cs {
                    let (r, ok) = iroot(n, e as u32);
                    if ok {
                        return Some((r, e as u64));
                    }
                }
                None
            }
        }
    }
}

fn smallest_odd_prime_factor(e: u64) -> u64 {
    let mut p = 3;
    while p <= e {
        if e % p == 0 && super::is_prime_u64(p) {
            return p;
        }
        p += 2;
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mersenne_bases_keep_complete_maximal_power_decompositions(){
        for bits in 2..=128 {
            let base=(BigInt::one()<<bits)-1u32;
            assert_eq!(perfect_power_big_uncached(&base),None);
            for exponent in [2,3,6,12,25] {
                assert_eq!(perfect_power_big_uncached(&base.pow(exponent)),Some((base.clone(),exponent as u64)));
            }
        }
    }

    fn b(n: i64) -> BigInt {
        BigInt::from(n)
    }

    #[test]
    fn rho_and_pm1() {
        let r = pollard_rho::<()>(&b(8051), &b(2), &b(1), 5, &b(1234), None, None).unwrap();
        assert_eq!(r, Some(b(97)));
        let n: BigInt = "21477639576571".parse().unwrap();
        assert_eq!(pollard_pm1(&n, 10, &b(2), 0, &b(1234)), None);
        assert_eq!(pollard_pm1(&n, 1787, &b(2), 0, &b(1234)), Some(b(4410317)));
    }

    #[test]
    fn divisors_and_powers() {
        assert_eq!(divisors_unsorted(&b(12), false), vec![b(1), b(2), b(4), b(3), b(6), b(12)]);
        assert_eq!(perfect_power(&b(-8), None, true), Some((b(-2), 3)));
        assert_eq!(perfect_power(&b(64), None, false), Some((b(8), 2)));
        assert_eq!(perfect_power(&b(64), None, true), Some((b(2), 6)));
        assert_eq!(perfect_power(&b(12), None, true), None);
        assert_eq!(abundance(&b(12)), b(4));
        assert_eq!(is_carmichael(&b(561)), Some(true));
        assert!(is_perfect(&b(8128)).unwrap());
        assert_eq!(multiplicity_in_factorial(&b(2), &b(10)).unwrap(), b(8));
    }
}
