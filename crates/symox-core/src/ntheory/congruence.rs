//! Polynomial congruences and binomials modulo composites, ported from
//! SymPy's `residue_ntheory.polynomial_congruence` (`gf_csolve`) and
//! `residue_ntheory.binomial_mod`.

use std::collections::BTreeSet;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::residue::linear_congruence;
use super::{crt, factorint, mod_inverse};

fn b(n: i64) -> BigInt {
    BigInt::from(n)
}

/// Horner evaluation of `f` (coefficients from the leading one) at `x`.
fn eval(f: &[BigInt], x: &BigInt) -> BigInt {
    f.iter().fold(BigInt::zero(), |acc, c| acc * x + c)
}

fn eval_mod(f: &[BigInt], x: &BigInt, p: &BigInt) -> BigInt {
    f.iter().fold(BigInt::zero(), |acc, c| (acc * x + c).mod_floor(p))
}

fn derivative(f: &[BigInt]) -> Vec<BigInt> {
    let d = f.len().saturating_sub(1);
    f.iter().take(d).enumerate().map(|(i, c)| c * (d - i)).collect()
}

// ---- dense polynomials over GF(p), leading coefficient first ----------

fn strip(mut f: Vec<BigInt>) -> Vec<BigInt> {
    let k = f.iter().take_while(|c| c.is_zero()).count();
    f.drain(..k);
    f
}

fn gf_rem(f: &[BigInt], g: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    let mut r: Vec<BigInt> = f.to_vec();
    if g.is_empty() {
        return r;
    }
    let inv = mod_inverse(&g[0], p).expect("monic divisor");
    while r.len() >= g.len() {
        let c = (&r[0] * &inv).mod_floor(p);
        if !c.is_zero() {
            for (i, gi) in g.iter().enumerate() {
                r[i] = (&r[i] - &c * gi).mod_floor(p);
            }
        }
        r.remove(0);
    }
    strip(r)
}

fn gf_quo(f: &[BigInt], g: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    let mut r: Vec<BigInt> = f.to_vec();
    let inv = mod_inverse(&g[0], p).expect("monic divisor");
    let mut q = vec![];
    while r.len() >= g.len() {
        let c = (&r[0] * &inv).mod_floor(p);
        for (i, gi) in g.iter().enumerate() {
            r[i] = (&r[i] - &c * gi).mod_floor(p);
        }
        r.remove(0);
        q.push(c);
    }
    strip(q)
}

fn gf_mul(f: &[BigInt], g: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    if f.is_empty() || g.is_empty() {
        return vec![];
    }
    let mut r = vec![BigInt::zero(); f.len() + g.len() - 1];
    for (i, a) in f.iter().enumerate() {
        for (j, c) in g.iter().enumerate() {
            r[i + j] += a * c;
        }
    }
    strip(r.into_iter().map(|c| c.mod_floor(p)).collect())
}

fn gf_monic(f: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    if f.is_empty() {
        return f;
    }
    let inv = mod_inverse(&f[0], p).expect("nonzero lead");
    f.into_iter().map(|c| (c * &inv).mod_floor(p)).collect()
}

fn gf_gcd(f: &[BigInt], g: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    let (mut a, mut c) = (f.to_vec(), g.to_vec());
    while !c.is_empty() {
        let r = gf_rem(&a, &c, p);
        a = c;
        c = r;
    }
    gf_monic(a, p)
}

/// `base**e mod (f, p)`.
fn gf_pow_mod(base: &[BigInt], e: &BigInt, f: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    let mut r = vec![BigInt::one()];
    let mut x = gf_rem(base, f, p);
    for i in 0..e.bits() {
        if e.bit(i) {
            r = gf_rem(&gf_mul(&r, &x, p), f, p);
        }
        x = gf_rem(&gf_mul(&x, &x, p), f, p);
    }
    r
}

fn gf_sub_one(mut f: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    if f.is_empty() {
        return vec![(p - 1u32).mod_floor(p)];
    }
    let last = f.len() - 1;
    f[last] = (&f[last] - 1u32).mod_floor(p);
    strip(f)
}

/// Distinct roots in `[0, p)` of `f` modulo the prime `p` (sorted).
fn roots_mod_prime(f: &[BigInt], p: &BigInt) -> Result<Vec<BigInt>, String> {
    let g = strip(f.iter().map(|c| c.mod_floor(p)).collect());
    if g.is_empty() {
        return match p.to_u64().filter(|&n| n <= 1 << 24) {
            Some(n) => Ok((0..n).map(BigInt::from).collect()),
            None => Err(format!("every residue modulo {p} is a root")),
        };
    }
    let small = p.to_u64().filter(|&n| n <= 4096 || (g.len() as u64).pow(2) >= n);
    if let Some(n) = small {
        return Ok((0..n).map(BigInt::from).filter(|x| eval_mod(&g, x, p).is_zero()).collect());
    }
    let mut roots = BTreeSet::new();
    let mut g = g;
    if g.last().is_some_and(|c| c.is_zero()) {
        roots.insert(BigInt::zero());
        while g.last().is_some_and(|c| c.is_zero()) {
            g.pop();
        }
    }
    let x = vec![BigInt::one(), BigInt::zero()];
    let h = gf_sub_one(gf_pow_mod(&x, &(p - 1u32), &g, p), p);
    let mut todo = vec![gf_gcd(&g, &h, p)];
    let mut a = BigInt::one();
    while let Some(f) = todo.pop() {
        match f.len() {
            0 | 1 => continue,
            2 => {
                roots.insert((-&f[1] * mod_inverse(&f[0], p).unwrap()).mod_floor(p));
                continue;
            }
            _ => {}
        }
        loop {
            let s = gf_sub_one(gf_pow_mod(&[BigInt::one(), a.clone()], &((p - 1u32) >> 1), &f, p), p);
            a += 1;
            let d = gf_gcd(&f, &s, p);
            if 1 < d.len() && d.len() < f.len() {
                todo.push(gf_quo(&f, &d, p));
                todo.push(d);
                break;
            }
        }
    }
    Ok(roots.into_iter().collect())
}

/// SymPy's `csolve_prime`: roots of `f` modulo `p**e`, lifted from the roots
/// modulo `p` with the same linear congruences.
fn roots_mod_prime_power(f: &[BigInt], p: &BigInt, e: u32) -> Result<Vec<BigInt>, String> {
    let x1 = roots_mod_prime(f, p)?;
    if e == 1 {
        return Ok(x1);
    }
    let df = derivative(f);
    let mut out = vec![];
    let mut stack: Vec<(BigInt, u32)> = x1.into_iter().map(|x| (x, 1)).collect();
    while let Some((x, s)) = stack.pop() {
        if s == e {
            out.push(x);
            continue;
        }
        let ps = p.pow(s);
        let alpha = eval(&df, &x);
        let beta = -(eval(f, &x).div_floor(&ps));
        for v in linear_congruence(&alpha, &beta, p) {
            stack.push((&x + v * &ps, s + 1));
        }
    }
    out.sort();
    Ok(out)
}

/// SymPy's `gf_csolve(f, n)`: every root of `f` modulo `n` (sorted).
pub fn csolve(f: &[BigInt], n: &BigInt) -> Result<Vec<BigInt>, String> {
    let fac = factorint(n);
    let mut sols: Vec<(BigInt, BigInt)> = vec![(BigInt::zero(), BigInt::one())];
    for (p, e) in &fac {
        let q = p.pow(*e);
        let rs = roots_mod_prime_power(f, p, *e)?;
        let mut next = Vec::with_capacity(sols.len() * rs.len());
        for (x, m) in &sols {
            for r in &rs {
                next.push(crt(&[x.clone(), r.clone()], &[m.clone(), q.clone()]).unwrap());
            }
        }
        sols = next;
    }
    let mut v: Vec<BigInt> = sols.into_iter().map(|(x, _)| x).collect();
    v.sort();
    Ok(v)
}

// ---- binomial coefficients modulo k ------------------------------------

/// `C(n, m) mod p` by Lucas' theorem.
fn binomial_mod_prime(n: &BigInt, m: &BigInt, p: &BigInt) -> BigInt {
    let (mut n, mut m) = (n.clone(), m.clone());
    let mut r = BigInt::one();
    while !m.is_zero() {
        let (nq, ni) = n.div_mod_floor(p);
        let (mq, mi) = m.div_mod_floor(p);
        if mi > ni {
            return BigInt::zero();
        }
        let k = std::cmp::min(mi.clone(), &ni - &mi);
        let (mut num, mut den) = (BigInt::one(), BigInt::one());
        let mut i = BigInt::zero();
        while i < k {
            num = (num * (&ni - &i)).mod_floor(p);
            i += 1;
            den = (den * &i).mod_floor(p);
        }
        r = (r * num * mod_inverse(&den, p).unwrap()).mod_floor(p);
        n = nq;
        m = mq;
    }
    r
}

/// Product of the integers in `[1, x]` coprime to `p`, modulo `q = p**e`.
struct CoprimeFactorial {
    p: BigInt,
    q: BigInt,
    table: Option<Vec<u64>>,
    full: BigInt,
}

impl CoprimeFactorial {
    fn new(p: &BigInt, q: &BigInt) -> Self {
        let mut s = CoprimeFactorial { p: p.clone(), q: q.clone(), table: None, full: BigInt::zero() };
        if let (Some(qq), Some(pp)) = (q.to_u64(), p.to_u64()) {
            if qq <= 1 << 22 {
                let mut t = Vec::with_capacity(qq as usize + 1);
                let mut acc = 1u128;
                t.push(1u64);
                for i in 1..=qq {
                    if i % pp != 0 {
                        acc = acc * i as u128 % qq as u128;
                    }
                    t.push(acc as u64);
                }
                s.table = Some(t);
            }
        }
        s.full = s.partial(q);
        s
    }

    fn partial(&self, x: &BigInt) -> BigInt {
        if let Some(t) = &self.table {
            return BigInt::from(t[x.to_usize().unwrap()]);
        }
        let mut acc = BigInt::one();
        let mut i = BigInt::one();
        while &i <= x {
            if !(&i % &self.p).is_zero() {
                acc = (acc * &i) % &self.q;
            }
            i += 1;
        }
        acc
    }

    fn eval(&self, x: &BigInt) -> BigInt {
        let (d, r) = x.div_mod_floor(&self.q);
        (self.full.modpow(&d, &self.q) * self.partial(&r)).mod_floor(&self.q)
    }
}

/// `C(n, m) mod p**e`.
fn binomial_mod_prime_power(n: &BigInt, m: &BigInt, p: &BigInt, e: u32) -> BigInt {
    if e == 1 {
        return binomial_mod_prime(n, m, p);
    }
    let q = p.pow(e);
    let r0: BigInt = n - m;
    // Kummer: the p-adic valuation is the number of carries of m + (n - m).
    let (mut a, mut c, mut carry, mut v) = (m.clone(), r0.clone(), 0u32, 0u64);
    while !(a.is_zero() && c.is_zero() && carry == 0) {
        let (aq, ad) = a.div_mod_floor(p);
        let (cq, cd) = c.div_mod_floor(p);
        let s: BigInt = ad + cd + carry;
        carry = if &s >= p { 1 } else { 0 };
        v += carry as u64;
        a = aq;
        c = cq;
    }
    if v >= e as u64 {
        return BigInt::zero();
    }
    let f = CoprimeFactorial::new(p, &q);
    let (mut nn, mut mm, mut rr) = (n.clone(), m.clone(), r0);
    let (mut num, mut den) = (BigInt::one(), BigInt::one());
    while !nn.is_zero() {
        num = (num * f.eval(&nn)).mod_floor(&q);
        den = (den * f.eval(&mm) * f.eval(&rr)).mod_floor(&q);
        nn /= p;
        mm /= p;
        rr /= p;
    }
    (num * mod_inverse(&den, &q).unwrap() * p.pow(v as u32)).mod_floor(&q)
}

/// SymPy's `binomial_mod(n, m, k)`: `binomial(n, m) % k`.
pub fn binomial_mod(n: &BigInt, m: &BigInt, k: &BigInt) -> Result<BigInt, String> {
    if k < &b(1) {
        return Err("k is required to be positive".into());
    }
    if n.is_negative() || m.is_negative() || m > n {
        return Ok(BigInt::zero());
    }
    let mut x = BigInt::zero();
    let mut md = BigInt::one();
    for (p, e) in factorint(k) {
        let r = binomial_mod_prime_power(n, m, &p, e);
        let (nx, nm) = crt(&[x, r], &[md, p.pow(e)]).unwrap();
        x = nx;
        md = nm;
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binom(n: u64, k: u64) -> BigInt {
        (0..k).fold(BigInt::one(), |acc, i| acc * (n - i) / (i + 1))
    }

    #[test]
    fn binomials() {
        for k in [1u64, 2, 4, 8, 9, 12, 27, 49, 100, 360, 1024, 7 * 7 * 7 * 11] {
            for n in 0..40u64 {
                for m in 0..=n {
                    let want = binom(n, m) % k;
                    assert_eq!(binomial_mod(&b(n as i64), &b(m as i64), &b(k as i64)).unwrap(), want, "{n} {m} {k}");
                }
            }
        }
    }

    #[test]
    fn congruences() {
        let f = [b(1), b(0), b(0), b(-4)];
        let want: Vec<BigInt> = (0..60).filter(|x| (x * x * x - 4) % 60 == 0).map(b).collect();
        assert_eq!(csolve(&f, &b(60)).unwrap(), want);
        let p = b(1_000_003);
        // (x - 5)(x - 7)(x^2 + 1) has roots 5, 7 and sqrt(-1) when p = 1 mod 4
        let f = [b(1), b(-12), b(36), b(-12), b(35)];
        let r = csolve(&f, &p).unwrap();
        assert!(r.contains(&b(5)) && r.contains(&b(7)));
        for x in &r {
            assert!(eval_mod(&f, x, &p).is_zero());
        }
    }
}
