//! Mersenne Twister (MT19937) reproducing CPython's `random.Random`.
//!
//! SymPy's `pollard_rho` and `pollard_pm1` reseed a `random.Random` with a
//! fixed seed and draw `randint` values between retries; reproducing the
//! exact stream makes the native ports return the same factors.

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Signed, Zero};

const N: usize = 624;
const M: usize = 397;

/// CPython-compatible MT19937 generator.
pub struct PyRandom {
    mt: [u32; N],
    index: usize,
}

impl PyRandom {
    fn init_genrand(s: u32) -> PyRandom {
        let mut mt = [0u32; N];
        mt[0] = s;
        for i in 1..N {
            mt[i] = 1812433253u32.wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_add(i as u32);
        }
        PyRandom { mt, index: N }
    }

    fn init_by_array(key: &[u32]) -> PyRandom {
        let mut r = PyRandom::init_genrand(19650218);
        let mt = &mut r.mt;
        let (mut i, mut j) = (1usize, 0usize);
        let mut k = N.max(key.len());
        while k > 0 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
            k -= 1;
        }
        k = N - 1;
        while k > 0 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1566083941)).wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
            k -= 1;
        }
        mt[0] = 0x8000_0000;
        r
    }

    /// `random.Random(seed)` for an integer seed (the absolute value is used).
    pub fn new(seed: &BigInt) -> PyRandom {
        let mut key = seed.abs().to_u32_digits().1;
        if key.is_empty() {
            key.push(0);
        }
        PyRandom::init_by_array(&key)
    }

    fn genrand_u32(&mut self) -> u32 {
        if self.index >= N {
            let mt = &mut self.mt;
            for kk in 0..N {
                let y = (mt[kk] & 0x8000_0000) | (mt[(kk + 1) % N] & 0x7fff_ffff);
                let mut v = mt[(kk + M) % N] ^ (y >> 1);
                if y & 1 != 0 {
                    v ^= 0x9908_b0df;
                }
                mt[kk] = v;
            }
            self.index = 0;
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// `getrandbits(k)`.
    pub fn getrandbits(&mut self, k: u64) -> BigUint {
        if k == 0 {
            return BigUint::zero();
        }
        if k <= 32 {
            return BigUint::from(self.genrand_u32() >> (32 - k));
        }
        let words = ((k - 1) / 32 + 1) as usize;
        let mut v = Vec::with_capacity(words);
        let mut left = k as i64;
        for _ in 0..words {
            let mut r = self.genrand_u32();
            if left < 32 {
                r >>= 32 - left;
            }
            v.push(r);
            left -= 32;
        }
        BigUint::new(v)
    }

    /// `_randbelow(n)` (Python >= 3.10 algorithm, `n > 0`).
    pub fn randbelow(&mut self, n: &BigUint) -> BigUint {
        let k = n.bits();
        loop {
            let r = self.getrandbits(k);
            if &r < n {
                return r;
            }
        }
    }

    /// `randint(a, b)` (inclusive bounds, `a <= b`).
    pub fn randint(&mut self, a: &BigInt, b: &BigInt) -> BigInt {
        let width: BigInt = b - a + 1;
        if width.sign() != Sign::Plus {
            return a.clone();
        }
        let w = width.magnitude().clone();
        if w.is_one() {
            return a.clone();
        }
        a + BigInt::from(self.randbelow(&w))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_cpython() {
        let mut r = PyRandom::new(&BigInt::from(1234));
        let v: Vec<u32> = (0..3).map(|_| r.genrand_u32()).collect();
        assert_eq!(v, vec![4150886329, 3342196574, 1892932127]);
        let mut r = PyRandom::new(&BigInt::from(1234));
        assert_eq!(r.randint(&BigInt::from(0), &BigInt::from(1_000_000)), BigInt::from(815965));
        let big: BigInt = "1000000000000000000000000000000".parse().unwrap();
        assert_eq!(r.randint(&BigInt::from(5), &big).to_string(), "79821698239120684918437567012");
        assert_eq!(r.getrandbits(70).to_string(), "1124191341926742866032");
        assert_eq!(PyRandom::new(&BigInt::from(0)).genrand_u32(), 3626764237);
        assert_eq!(PyRandom::new(&(BigInt::from(1u64 << 40) + 5)).genrand_u32(), 2166296868);
    }
}
