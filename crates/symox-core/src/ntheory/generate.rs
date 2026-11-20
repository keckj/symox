//! Prime generation helpers ported from SymPy's `ntheory.generate`.

use num_bigint::BigInt;
use num_traits::One;

use super::{prime, primepi, primerange};

/// `nextprime(n, ith)` (SymPy semantics, `ith >= 1`).
pub fn nextprime_ith(n: &BigInt, ith: u64) -> BigInt {
    let mut i = ith;
    let mut n = n.clone();
    if n < BigInt::from(2) {
        n = BigInt::from(2);
        i -= 1;
        if i == 0 {
            return n;
        }
    }
    loop {
        n = super::nextprime(&n);
        i -= 1;
        if i == 0 {
            return n;
        }
    }
}

/// `primorial(n, nth)`: product of the first `n` primes, or of the
/// primes `<= n`.
pub fn primorial(n: u64, nth: bool) -> BigInt {
    if nth {
        if n == 0 {
            return BigInt::one();
        }
        let last = prime(n);
        super::product_u64(&primerange(2, last + 1))
    } else {
        super::product_u64(&primerange(2, n.saturating_add(1)))
    }
}

/// Number of composites `<= n`.
pub fn compositepi(n: u64) -> u64 {
    if n < 4 {
        return 0;
    }
    n - primepi(n) - 1
}

/// The `nth` composite number (`nth >= 1`).
pub fn composite(nth: u64) -> u64 {
    // smallest a with a - pi(a) - 1 >= nth; that a is composite
    let (mut lo, mut hi) = (4u64, 2 * nth + 4);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if mid - primepi(mid) - 1 >= nth {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Euler totients of `lo..hi`.
pub fn totient_range(lo: u64, hi: u64) -> Vec<u64> {
    if hi <= lo {
        return vec![];
    }
    let mut t: Vec<u64> = (0..hi).collect();
    for i in 2..hi as usize {
        if t[i] == i as u64 {
            let mut j = i;
            while j < hi as usize {
                t[j] -= t[j] / i as u64;
                j += i;
            }
        }
    }
    t.drain(..lo as usize);
    t
}

/// Mobius function on `lo..hi`.
pub fn mobius_range(lo: u64, hi: u64) -> Vec<i64> {
    if hi <= lo {
        return vec![];
    }
    let h = hi as usize;
    let mut m = vec![0i64; h];
    if h > 1 {
        m[1] = 1;
    }
    for i in 1..h {
        let mi = m[i];
        if mi != 0 {
            let mut j = 2 * i;
            while j < h {
                m[j] -= mi;
                j += i;
            }
        }
    }
    m.drain(..lo as usize);
    m
}

/// Largest prime `< n` (`n >= 3`).
pub fn prevprime(n: &BigInt) -> Option<BigInt> {
    super::prevprime(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generate() {
        assert_eq!(nextprime_ith(&BigInt::from(2), 2), BigInt::from(5));
        assert_eq!(nextprime_ith(&BigInt::from(-5), 1), BigInt::from(2));
        assert_eq!(nextprime_ith(&BigInt::from(1), 3), BigInt::from(5));
        assert_eq!(primorial(4, true), BigInt::from(210));
        assert_eq!(primorial(10, false), BigInt::from(210));
        assert_eq!(composite(1), 4);
        assert_eq!(composite(36), 52);
        assert_eq!(compositepi(25), 15);
        assert_eq!(totient_range(7, 12), vec![6, 4, 6, 4, 10]);
        assert_eq!(mobius_range(7, 18), vec![-1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1]);
    }
}
