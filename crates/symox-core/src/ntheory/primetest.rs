//! Probable-prime tests ported from SymPy's `ntheory.primetest` (Fermat,
//! Euler, Euler-Jacobi, Miller-Rabin, Lucas, strong and extra strong
//! Lucas, Proth), with SymPy's argument checks and error messages.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use super::jacobi_big;

type R<T> = Result<T, String>;

fn b(n: i64) -> BigInt {
    BigInt::from(n)
}

fn bit_scan1(n: &BigInt) -> u64 {
    n.trailing_zeros().unwrap_or(0)
}

/// Exact perfect square test for `n >= 0`.
pub fn is_square(n: &BigInt) -> bool {
    if n.is_negative() {
        return false;
    }
    let r = n.sqrt();
    &r * &r == *n
}

fn jacobi(a: &BigInt, n: &BigInt) -> i32 {
    jacobi_big(a, n)
}

fn halve(x: BigInt, n: &BigInt) -> BigInt {
    let x = x.mod_floor(n);
    if x.is_odd() {
        (x + n) >> 1
    } else {
        x >> 1
    }
}

/// The modular Lucas sequence `(U_k, V_k, Q**k)` modulo odd `n`.
pub fn lucas_sequence(n: &BigInt, p: &BigInt, q: &BigInt, k: &BigInt) -> (BigInt, BigInt, BigInt) {
    if k.is_zero() {
        return (BigInt::zero(), b(2), BigInt::one());
    }
    let d: BigInt = p * p - q * 4;
    let mut u = BigInt::one();
    let mut v = p.clone();
    let mut qk = q.mod_floor(n);
    let bits = k.bits();
    for i in (0..bits - 1).rev() {
        u = (&u * &v).mod_floor(n);
        v = (&v * &v - &qk * BigInt::from(2)).mod_floor(n);
        qk = (&qk * &qk).mod_floor(n);
        if k.bit(i) {
            let u2 = halve(&u * p + &v, n);
            let v2 = halve(&v * p + &u * &d, n);
            u = u2;
            v = v2;
            qk = (&qk * q).mod_floor(n);
        }
    }
    (u.mod_floor(n), v.mod_floor(n), qk.mod_floor(n))
}

/// Shared prologue of gmpy-style `is_*_prp(n, a)`: an error, an early
/// answer (`Ok(Err(bool))`) or the reduced base (`Ok(Ok(a))`).
fn prp_prep(name: &str, n: &BigInt, a: &BigInt) -> R<Result<BigInt, bool>> {
    if a < &b(2) {
        return Err(format!("{name}() requires 'a' greater than or equal to 2"));
    }
    if n < &b(1) {
        return Err(format!("{name}() requires 'n' be greater than 0"));
    }
    if n.is_one() {
        return Ok(Err(false));
    }
    if n.is_even() {
        return Ok(Err(n == &b(2)));
    }
    let a = a.mod_floor(n);
    if !a.gcd(n).is_one() {
        return Err(format!("{name}() requires gcd(n,a) == 1"));
    }
    Ok(Ok(a))
}

/// SymPy's `is_fermat_pseudoprime(n, a)`.
pub fn is_fermat_pseudoprime(n: &BigInt, a: &BigInt) -> R<bool> {
    if a.is_one() {
        return Ok(n == &b(2) || n.is_odd());
    }
    Ok(match prp_prep("is_fermat_prp", n, a)? {
        Err(r) => r,
        Ok(a) => a.modpow(&(n - 1), n).is_one(),
    })
}

/// SymPy's `is_euler_pseudoprime(n, a)`.
pub fn is_euler_pseudoprime(n: &BigInt, a: &BigInt) -> R<bool> {
    if a < &b(1) {
        return Err("a should be an integer greater than 0".into());
    }
    if n < &b(1) {
        return Err("n should be an integer greater than 0".into());
    }
    if n.is_one() {
        return Ok(false);
    }
    if a.is_one() {
        return Ok(n == &b(2) || n.is_odd());
    }
    if n.is_even() {
        return Ok(n == &b(2));
    }
    if !n.gcd(a).is_one() {
        return Err("The two numbers should be relatively prime".into());
    }
    let r = a.modpow(&((n - 1) / 2), n);
    Ok(r.is_one() || r == n - 1)
}

/// SymPy's `is_euler_jacobi_pseudoprime(n, a)`.
pub fn is_euler_jacobi_pseudoprime(n: &BigInt, a: &BigInt) -> R<bool> {
    if a.is_one() {
        return Ok(n == &b(2) || n.is_odd());
    }
    Ok(match prp_prep("is_euler_prp", n, a)? {
        Err(r) => r,
        Ok(a) => a.modpow(&(n >> 1), n) == b(jacobi(&a, n) as i64).mod_floor(n),
    })
}

fn strong_test(n: &BigInt, base: &BigInt, s: u64, t: &BigInt) -> bool {
    let nm1: BigInt = n - 1;
    let mut x = base.modpow(t, n);
    if x.is_one() || x == nm1 {
        return true;
    }
    for _ in 1..s {
        x = (&x * &x) % n;
        if x == nm1 {
            return true;
        }
        if x.is_one() {
            return false;
        }
    }
    false
}

/// SymPy's `mr(n, bases)` Miller-Rabin test.
pub fn mr(n: &BigInt, bases: &[BigInt]) -> bool {
    if n < &b(2) || (n > &b(2) && n.is_even()) {
        return false;
    }
    let nm1: BigInt = n - 1;
    let s = bit_scan1(&nm1);
    let t: BigInt = n >> s;
    for base in bases {
        let base = if base >= n { base % n } else { base.clone() };
        if base >= b(2) && !strong_test(n, &base, s, &t) {
            return false;
        }
    }
    true
}

/// Selfridge parameter search shared by the (strong) Lucas tests:
/// `Ok(Some(D))` with `jacobi(D, n) == -1`, `Ok(None)` when n is composite.
fn selfridge_d(n: &BigInt, name: &str) -> R<Option<BigInt>> {
    let mut d = 5i64;
    while d < 1_000_000 {
        let dd = if d & 2 != 0 { -d } else { d };
        let big = b(dd);
        let j = jacobi(&big, n);
        if j == -1 {
            return Ok(Some(big));
        }
        if j == 0 && !(big.mod_floor(n)).is_zero() {
            return Ok(None);
        }
        if d == 13 && is_square(n) {
            return Ok(None);
        }
        d += 2;
    }
    Err(format!("appropriate value for D cannot be found in {name}()"))
}

/// SymPy's `is_lucas_prp(n)` (Selfridge parameters).
pub fn is_lucas_prp(n: &BigInt) -> R<bool> {
    if n < &b(2) {
        return Ok(false);
    }
    if n.is_even() {
        return Ok(n == &b(2));
    }
    let Some(d) = selfridge_d(n, "is_selfridge_prp")? else {
        return Ok(false);
    };
    let q: BigInt = (BigInt::one() - &d) / 4;
    Ok(lucas_sequence(n, &BigInt::one(), &q, &(n + 1)).0.is_zero())
}

/// SymPy's `is_strong_lucas_prp(n)` (Selfridge parameters).
pub fn is_strong_lucas_prp(n: &BigInt) -> R<bool> {
    if n < &b(2) {
        return Ok(false);
    }
    if n.is_even() {
        return Ok(n == &b(2));
    }
    let Some(d) = selfridge_d(n, "is_strong_selfridge_prp")? else {
        return Ok(false);
    };
    let q: BigInt = (BigInt::one() - &d) / 4;
    let np1: BigInt = n + 1;
    let s = bit_scan1(&np1);
    let (u, mut v, mut qk) = lucas_sequence(n, &BigInt::one(), &q, &(&np1 >> s));
    if u.is_zero() || v.is_zero() {
        return Ok(true);
    }
    for _ in 1..s {
        v = (&v * &v - &qk * BigInt::from(2)).mod_floor(n);
        if v.is_zero() {
            return Ok(true);
        }
        qk = qk.modpow(&b(2), n);
    }
    Ok(false)
}

/// SymPy's `is_extra_strong_lucas_prp(n)`.
pub fn is_extra_strong_lucas_prp(n: &BigInt) -> bool {
    if n == &b(2) {
        return true;
    }
    if n < &b(2) || n.is_even() || is_square(n) {
        return false;
    }
    let mut p = 3i64;
    let pp = loop {
        let d = b(p * p - 4);
        let j = jacobi(&d, n);
        if j == -1 {
            break b(p);
        }
        if j == 0 && !d.mod_floor(n).is_zero() {
            return false;
        }
        p += 1;
    };
    let np1: BigInt = n + 1;
    let s = bit_scan1(&np1);
    let k: BigInt = &np1 >> s;
    let (u, mut v, _) = lucas_sequence(n, &pp, &BigInt::one(), &k);
    if u.is_zero() && (v == b(2) || v == n - 2) {
        return true;
    }
    for _ in 1..s {
        if v.is_zero() {
            return true;
        }
        v = (&v * &v - BigInt::from(2)).mod_floor(n);
    }
    false
}

/// SymPy's `proth_test(n)`.
pub fn proth_test(n: &BigInt) -> R<bool> {
    let err = || Err("n is not Proth number".to_string());
    if n < &b(3) {
        return err();
    }
    let nm1: BigInt = n - 1;
    let m = bit_scan1(&nm1);
    let k: BigInt = n >> m;
    if m < k.bits() {
        return err();
    }
    if (n % 3u32).is_zero() {
        return Ok(n == &b(3));
    }
    if !(&k % 3u32).is_zero() {
        return Ok(b(3).modpow(&(n >> 1), n) == nm1);
    }
    if is_square(n) {
        return Ok(false);
    }
    let mut a = b(5);
    while &a < n {
        let j = jacobi(&a, n);
        if j == -1 {
            return Ok(a.modpow(&(n >> 1), n) == nm1);
        }
        if j == 0 {
            return Ok(false);
        }
        a += 1;
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pseudoprimes() {
        assert_eq!(is_fermat_pseudoprime(&b(341), &b(2)), Ok(true));
        assert_eq!(is_euler_pseudoprime(&b(341), &b(2)), Ok(true));
        assert_eq!(is_euler_jacobi_pseudoprime(&b(561), &b(2)), Ok(true));
        assert!(mr(&b(1373653), &[b(2)]));
        assert!(mr(&b(1373653), &[b(2), b(3)]));
        assert!(!mr(&b(1373653), &[b(5)]));
        assert_eq!(is_lucas_prp(&b(323)), Ok(true));
        assert_eq!(is_strong_lucas_prp(&b(5459)), Ok(true));
        assert!(is_extra_strong_lucas_prp(&b(989)));
        assert_eq!(proth_test(&b(97)), Ok(true));
        assert!(proth_test(&b(7)).is_err());
    }
}
