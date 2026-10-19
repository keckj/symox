//! Polynomial congruences and binomials modulo composites, ported from
//! SymPy's `residue_ntheory.polynomial_congruence` (`gf_csolve`) and
//! `residue_ntheory.binomial_mod` (BSD-3-Clause, see LICENSE-SYMPY).

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

fn signed_modpow(base: &BigInt, exponent: &BigInt, modulus: &BigInt) -> BigInt {
    if exponent.is_negative() {
        mod_inverse(base, modulus).expect("unit factorial is invertible")
            .modpow(&(-exponent), modulus)
    } else {
        base.modpow(exponent, modulus)
    }
}

/// Integer-valued Lagrange weights for nodes 0,1,...,degree or their squares.
fn interpolation_weight(u: &BigInt, j: u32, degree: u32, squared: bool) -> BigInt {
    let (mut numerator, mut denominator) = (u.clone(), BigInt::from(j));
    let (u_node, j_node) = if squared {
        (u * u, BigInt::from(j) * j)
    } else {
        (u.clone(), BigInt::from(j))
    };
    for i in 1..=degree {
        if i == j { continue; }
        let node = if squared { BigInt::from(i) * i } else { BigInt::from(i) };
        numerator *= &u_node - &node;
        denominator *= &j_node - node;
    }
    debug_assert!((&numerator % &denominator).is_zero());
    numerator / denominator
}

/// Granville's prime-power unit-factorial interpolation, following the pinned
/// reference algorithm. Work is polynomial in the exponent and linear in p,
/// rather than linear in p**e. Extra p-adic precision handles exceptional
/// Bernoulli denominators and the sign correction for p=2.
struct GranvilleFactorial {
    p: BigInt,
    q: BigInt,
    exponent: u32,
    half: u32,
    working_modulus: BigInt,
    block_factorials: Vec<BigInt>,
    factorials: Vec<BigInt>,
}

impl GranvilleFactorial {
    fn new(p: &BigInt, q: &BigInt, exponent: u32) -> Self {
        let mut half = exponent / 2;
        let exceptional = (half == 1 && *p == b(2))
            || BigInt::from(2 * half + 1) == *p
            || BigInt::from(2 * half + 1) == p * p;
        let working_exponent = if exceptional {
            if exponent % 2 == 1 { half += 1; }
            2 * half
        } else { 2 * half + 1 };
        let working_modulus = p.pow(working_exponent);
        let mut block_factorials = Vec::new();
        let mut product = BigInt::one();
        for j in 1..=half {
            let mut value: BigInt = BigInt::from(j - 1) * p + 1;
            let end = BigInt::from(j) * p;
            while value < end {
                product = (product * &value).mod_floor(&working_modulus);
                value += 1;
            }
            block_factorials.push(product.clone());
        }
        Self { p: p.clone(), q: q.clone(), exponent, half, working_modulus,
            block_factorials, factorials: vec![BigInt::one()] }
    }

    fn factorial(&mut self, n: &BigInt) -> BigInt {
        if let Some(index) = n.to_usize() {
            while self.factorials.len() <= index {
                let next = self.factorials.len();
                self.factorials.push((&self.factorials[next - 1] * next).mod_floor(&self.q));
            }
            return self.factorials[index].clone();
        }
        // No narrowing conversion or false result for enormous integer inputs.
        let (mut product, mut value) = (BigInt::one(), BigInt::one());
        while &value <= n {
            product = (product * &value).mod_floor(&self.q);
            value += 1;
        }
        product
    }

    fn full_blocks(&self, u: &BigInt) -> BigInt {
        let mut product = BigInt::one();
        let mut sign: BigInt = u / 2;
        for (index, factorial) in self.block_factorials.iter().enumerate() {
            let j = index as u32 + 1;
            let exponent = interpolation_weight(u, j, self.half, true);
            product = (product * signed_modpow(factorial, &exponent, &self.working_modulus))
                .mod_floor(&self.working_modulus);
            if self.p == b(2) { sign += &exponent * (j / 2); }
        }
        if self.p == b(2) && sign.is_odd() { product = -product; }
        product.mod_floor(&self.q)
    }

    fn partial_block(&self, u: &BigInt, v: &BigInt, factorial: &BigInt) -> BigInt {
        if v.is_zero() { return BigInt::one(); }
        let inverse = mod_inverse(factorial, &self.q).expect("v < p");
        let mut product = BigInt::one();
        for j in 1..self.exponent {
            let mut coefficient = inverse.clone();
            let start = BigInt::from(j) * &self.p;
            let mut value = BigInt::one();
            while &value <= v {
                coefficient = (coefficient * (&start + &value)).mod_floor(&self.q);
                value += 1;
            }
            let exponent = interpolation_weight(u, j, self.exponent - 1, false);
            product = (product * signed_modpow(&coefficient, &exponent, &self.q)).mod_floor(&self.q);
        }
        product
    }

    fn eval(&mut self, n: &BigInt) -> BigInt {
        let (u, v) = n.div_mod_floor(&self.p);
        let factorial = self.factorial(&v);
        if u.is_zero() { return factorial; }
        (&factorial * self.full_blocks(&u) * self.partial_block(&u, &v, &factorial)).mod_floor(&self.q)
    }
}

/// Product of the integers in `[1, x]` coprime to `p`, modulo `q = p**e`.
struct CoprimeFactorial {
    q: BigInt,
    table: Option<Vec<u64>>,
    granville: Option<GranvilleFactorial>,
    full: BigInt,
}

impl CoprimeFactorial {
    fn new(p: &BigInt, q: &BigInt, exponent: u32) -> Self {
        let mut s = CoprimeFactorial { q: q.clone(), table: None, granville: None,
            full: if *p == b(2) && exponent >= 3 { BigInt::one() } else { q - 1 } };
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
        if s.table.is_none() { s.granville = Some(GranvilleFactorial::new(p, q, exponent)); }
        s
    }

    fn partial(&mut self, x: &BigInt) -> BigInt {
        if let Some(t) = &self.table {
            return BigInt::from(t[x.to_usize().unwrap()]);
        }
        self.granville.as_mut().expect("table or interpolation").eval(x)
    }

    fn eval(&mut self, x: &BigInt) -> BigInt {
        let (d, r) = x.div_mod_floor(&self.q);
        let sign = if self.full.is_one() || d.is_even() { BigInt::one() } else { self.full.clone() };
        (sign * self.partial(&r)).mod_floor(&self.q)
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
    let mut f = CoprimeFactorial::new(p, &q, e);
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

/// Exact n! modulo a positive modulus, grouping primes with equal valuations.
/// This follows the reference factorial._facmod strategy rather than creating
/// the factorial, whose bit length can dwarf the modulus and final result.
pub(crate) fn factorial_mod(n: &BigInt, modulus: &BigInt) -> Option<BigInt> {
    let n = n.to_u64()?;
    let root = num_integer::Roots::sqrt(&n);
    let mut grouped = std::collections::BTreeMap::<u64, BigInt>::new();
    let mut result = BigInt::one();
    for prime in super::primerange(2, n.checked_add(1)?) {
        let (mut valuation, mut quotient) = (0u64, n / prime);
        while quotient != 0 { valuation += quotient; quotient /= prime; }
        if valuation < root {
            let product = grouped.entry(valuation).or_insert_with(BigInt::one);
            *product = (&*product * prime).mod_floor(modulus);
        } else {
            result = (result * BigInt::from(prime).modpow(&valuation.into(), modulus)).mod_floor(modulus);
        }
    }
    for (exponent, base) in grouped {
        if base.is_zero() { return Some(BigInt::zero()); }
        result = (result * base.modpow(&exponent.into(), modulus)).mod_floor(modulus);
    }
    Some(result.mod_floor(modulus))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binom(n: u64, k: u64) -> BigInt {
        (0..k).fold(BigInt::one(), |acc, i| acc * (n - i) / (i + 1))
    }

    #[test]
    fn granville_unit_factorials_match_independent_products() {
        for prime in [2, 3, 5, 7, 11, 19] {
            let p = b(prime);
            for exponent in 2..=10 {
                let q = p.pow(exponent);
                let mut factorial = GranvilleFactorial::new(&p, &q, exponent);
                let mut exact = BigInt::one();
                for n in 0..=100 {
                    if n % prime != 0 { exact = (exact * n).mod_floor(&q); }
                    assert_eq!(factorial.eval(&b(n)), exact, "p={prime} e={exponent} n={n}");
                }
            }
        }
    }

    #[test]
    fn huge_binomials_mod_prime_powers_without_materialization() {
        for (n, k, modulus, expected) in [
            ("1000000000000000000", "1000000000000", "10000600009", "3744312326"),
            ("1000000000000000000000000000000", "1000000000000000", "4294967296", "4194304000"),
            ("100000000000000000099", "1000000001", "3486784401", "39051072"),
            ("100000000000000000999", "10000000000", "244140625", "9765625"),
        ] {
            let parse = |s: &str| s.parse::<BigInt>().unwrap();
            assert_eq!(binomial_mod(&parse(n), &parse(k), &parse(modulus)).unwrap(), parse(expected));
        }
    }

    #[test]
    fn modular_factorials_match_exact_small_factorials() {
        let mut exact = BigInt::one();
        for n in 0..100 {
            if n > 0 { exact *= n; }
            for modulus in [1, 2, 3, 4, 7, 24, 97, 100003] {
                assert_eq!(factorial_mod(&b(n), &b(modulus)).unwrap(), exact.mod_floor(&b(modulus)));
            }
        }
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
