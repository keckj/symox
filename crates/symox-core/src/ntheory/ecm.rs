//! Lenstra's elliptic curve factorization with Montgomery curves, ported
//! from SymPy's `ntheory.ecm` (same curve choice from the seeded Mersenne
//! twister, same stage 1 / stage 2 bounds).

use num_bigint::BigInt;
use num_integer::{Integer, Roots};
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::factor::perfect_power_big;
use super::mt::PyRandom;
use super::{isprime, mod_inverse, primerange};

/// A projective point `(X : Z)` on `b*y**2 = x**3 + a*x**2 + x` with
/// `a24 = (a + 2) / 4`, modulo `m`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point {
    pub x: BigInt,
    pub z: BigInt,
    pub a24: BigInt,
    pub m: BigInt,
}

impl Point {
    pub fn new(x: BigInt, z: BigInt, a24: BigInt, m: BigInt) -> Point {
        Point { x, z, a24, m }
    }

    /// Point equality as in SymPy (`X1 * Z2 == X2 * Z1` modulo `m`).
    pub fn same(&self, o: &Point) -> bool {
        self.a24 == o.a24
            && self.m == o.m
            && (&self.x * &o.z).mod_floor(&self.m) == (&o.x * &self.z).mod_floor(&self.m)
    }

    /// Differential addition `self + q` knowing `diff = self - q`.
    pub fn add(&self, q: &Point, diff: &Point) -> Point {
        let u = (&self.x - &self.z) * (&q.x + &q.z);
        let v = (&self.x + &self.z) * (&q.x - &q.z);
        let (add, sub) = (&u + &v, &u - &v);
        let x = (&diff.z * &add * &add).mod_floor(&self.m);
        let z = (&diff.x * &sub * &sub).mod_floor(&self.m);
        Point::new(x, z, self.a24.clone(), self.m.clone())
    }

    pub fn double(&self) -> Point {
        let s: BigInt = &self.x + &self.z;
        let d: BigInt = &self.x - &self.z;
        let u = (&s * &s).mod_floor(&self.m);
        let v = (&d * &d).mod_floor(&self.m);
        let diff: BigInt = &u - &v;
        let x = (&u * &v).mod_floor(&self.m);
        let z = (&diff * (&v + &self.a24 * &diff)).mod_floor(&self.m);
        Point::new(x, z, self.a24.clone(), self.m.clone())
    }

    /// `k * self` by the Montgomery ladder (`k >= 1`).
    pub fn mont_ladder(&self, k: &BigInt) -> Point {
        let mut q = self.clone();
        let mut r = self.double();
        let bits = k.bits();
        for i in (0..bits.saturating_sub(1)).rev() {
            if k.bit(i) {
                q = r.add(&q, self);
                r = r.double();
            } else {
                r = q.add(&r, self);
                q = q.double();
            }
        }
        q
    }
}

/// Python 3.10's message for `randint(6, n - 1)` on an empty range.
fn empty_range(n: &BigInt) -> String {
    format!("empty range for randrange() (6, {}, {})", n, n - 6)
}

pub(crate) fn one_factor(n: &BigInt, b1: u64, b2: u64, max_curve: u64, seed: &BigInt) -> Result<Option<BigInt>, String> {
    let mut rng = PyRandom::new(seed);
    if n - 1 < BigInt::from(6) {
        return Err(empty_range(n));
    }
    let d = std::cmp::min(b2.sqrt(), (b1 / 2).saturating_sub(1)) as usize;
    if d < 2 {
        return Err("list index out of range".into());
    }
    let mut k = BigInt::one();
    for p in primerange(2, b1 + 1) {
        let e = ((b1 as f64).ln() / (p as f64).ln()) as u32;
        k *= BigInt::from(p).pow(e);
    }
    let dd = d as u64;
    let mut deltas_list: Vec<Vec<usize>> = vec![];
    let mut r = b1 + 2 * dd;
    while r < b2 + 2 * dd {
        let mut ds: Vec<usize> = primerange(r - 2 * dd, r + 2 * dd).into_iter().map(|q| (q.abs_diff(r) >> 1) as usize).collect();
        ds.sort();
        ds.dedup();
        deltas_list.push(ds);
        r += 4 * dd;
    }
    let lo = BigInt::from(6);
    let hi: BigInt = n - 1;
    for _ in 0..max_curve {
        let sigma = rng.randint(&lo, &hi);
        let u = (&sigma * &sigma - BigInt::from(5)).mod_floor(n);
        let v = (&sigma * BigInt::from(4)).mod_floor(n);
        let u3 = u.modpow(&BigInt::from(3), n);
        let den: BigInt = (BigInt::from(16) * &u3 * &v).mod_floor(n);
        let a24 = match mod_inverse(&den, n) {
            Some(inv) => ((&v - &u).modpow(&BigInt::from(3), n) * (BigInt::from(3) * &u + &v) * inv).mod_floor(n),
            None => {
                let g = (BigInt::from(2) * &u3 * &v).gcd(n);
                if &g == n {
                    continue;
                }
                return Ok(Some(g));
            }
        };
        let q = Point::new(u3, v.modpow(&BigInt::from(3), n), a24, n.clone());
        let q = q.mont_ladder(&k);
        let g = q.z.gcd(n);
        if !g.is_one() && &g != n {
            return Ok(Some(g));
        } else if &g == n {
            continue;
        }
        let q2 = q.double();
        let mut s: Vec<Point> = Vec::with_capacity(d);
        s.push(q.clone());
        s.push(q2.add(&q, &q));
        for i in 2..d {
            let next = s[i - 1].add(&q2, &s[i - 2]);
            s.push(next);
        }
        let beta: Vec<BigInt> = s.iter().map(|p| (&p.x * &p.z).mod_floor(n)).collect();
        let mut g = BigInt::one();
        let w = q.mont_ladder(&BigInt::from(4 * dd));
        let mut t = q.mont_ladder(&BigInt::from(b1 - 2 * dd));
        let mut rr = q.mont_ladder(&BigInt::from(b1 + 2 * dd));
        for deltas in &deltas_list {
            let alpha = (&rr.x * &rr.z).mod_floor(n);
            for &delta in deltas {
                let sd = &s[delta];
                let f = (&rr.x - &sd.x) * (&rr.z + &sd.z) - &alpha + &beta[delta];
                g = (g * f).mod_floor(n);
            }
            let next = rr.add(&w, &t);
            t = std::mem::replace(&mut rr, next);
        }
        let g = n.gcd(&g);
        if !g.is_one() && &g != n {
            return Ok(Some(g));
        }
    }
    Ok(None)
}

/// SymPy's `ecm(n, B1, B2, max_curve, seed)`: the set of prime factors of
/// `n` (sorted), with SymPy's error messages.
pub fn ecm(n: &BigInt, b1: &BigInt, b2: &BigInt, max_curve: u64, seed: &BigInt) -> Result<Vec<BigInt>, String> {
    if b1.is_odd() || b2.is_odd() {
        return Err("both bounds must be even".into());
    }
    if n.is_zero() {
        return Err("n should be nonzero".into());
    }
    let (b1, b2) = match (b1.to_u64(), b2.to_u64()) {
        (Some(x), Some(y)) => (x, y),
        _ => return Err("bounds too large".into()),
    };
    const TF_LIMIT: u64 = 100_000;
    let mut factors: Vec<BigInt> = vec![];
    let mut n = n.clone();
    for p in primerange(2, TF_LIMIT) {
        if (&n % p).is_zero() {
            factors.push(BigInt::from(p));
            while (&n % p).is_zero() {
                n /= p;
            }
        }
    }
    let mut queue: Vec<BigInt> = vec![];
    fn check(m: BigInt, factors: &mut Vec<BigInt>, queue: &mut Vec<BigInt>) {
        if isprime(&m) {
            factors.push(m);
            return;
        }
        if m.is_positive() {
            if let Some((b, _)) = perfect_power_big(&m) {
                return check(b, factors, queue);
            }
        }
        queue.push(m);
    }
    check(n, &mut factors, &mut queue);
    while let Some(n) = queue.pop() {
        let Some(f) = one_factor(&n, b1, b2, max_curve, seed)? else {
            return Err("Increase the bounds".into());
        };
        let cof = &n / &f;
        check(f, &mut factors, &mut queue);
        check(cof, &mut factors, &mut queue);
    }
    factors.sort();
    factors.dedup();
    Ok(factors)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(s: &str) -> BigInt {
        s.parse().unwrap()
    }

    #[test]
    fn factors() {
        let r = ecm(&b("398883434337287"), &b("10000"), &b("100000"), 200, &b("1234")).unwrap();
        assert_eq!(r, vec![b("4009823"), b("99476569")]);
        assert_eq!(ecm(&b("1"), &b("10000"), &b("100000"), 200, &b("1234")).unwrap_err(), "empty range for randrange() (6, 1, -5)");
        let p = Point::new(b("11"), b("16"), b("7"), b("29"));
        let q = p.double();
        assert_eq!((q.x.clone(), q.z.clone()), (b("13"), b("10")));
        assert_eq!(p.mont_ladder(&b("3")), q.add(&p, &p));
    }
}
