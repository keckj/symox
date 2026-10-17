//! Self-initializing quadratic sieve, ported from SymPy's `ntheory.qs`
//! (same factor base, seeded polynomial choice, sieve thresholds, partial
//! relations and GF(2) elimination, hence the same factors in the same
//! order).

use std::collections::HashMap;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::factor::perfect_power_big;
use super::mt::PyRandom;
use super::residue::sqrt_mod_iter;
use super::{isprime, mod_inverse, primerange};

/// `math.log` of an integer as CPython computes it.
fn ln(n: &BigInt) -> f64 {
    match n.to_f64() {
        Some(f) if f.is_finite() => f.ln(),
        _ => {
            let bits = n.bits();
            let shift = bits - 53;
            let m = (n >> shift).to_f64().unwrap();
            m.ln() + (shift as f64) * std::f64::consts::LN_2
        }
    }
}

struct Fb {
    prime: u64,
    tmem_p: u64,
    log_p: i64,
    soln1: Option<u64>,
    soln2: u64,
    b_ainv: Vec<u64>,
}

#[derive(Clone)]
struct Bits(Vec<u64>);

impl Bits {
    fn new(n: usize) -> Bits {
        Bits(vec![0; n / 64 + 1])
    }
    fn get(&self, i: usize) -> bool {
        self.0[i / 64] >> (i % 64) & 1 == 1
    }
    fn flip(&mut self, i: usize) {
        self.0[i / 64] ^= 1 << (i % 64);
    }
    fn xor(&mut self, o: &Bits) {
        for (a, b) in self.0.iter_mut().zip(&o.0) {
            *a ^= *b;
        }
    }
    fn intersects(&self, o: &Bits) -> bool {
        self.0.iter().zip(&o.0).any(|(a, b)| a & b != 0)
    }
}

fn pmod(a: &BigInt, p: u64) -> u64 {
    a.mod_floor(&BigInt::from(p)).to_u64().unwrap()
}

fn inv_mod(a: u64, p: u64) -> u64 {
    mod_inverse(&BigInt::from(a), &BigInt::from(p)).unwrap().to_u64().unwrap()
}

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}

struct Relation {
    u: BigInt,
    v: BigInt,
    vec: Bits,
}

fn check_smoothness(num: &BigInt, fb: &[Fb], width: usize) -> (Bits, BigInt) {
    let mut vec = Bits::new(width);
    let mut num = num.clone();
    if num.is_negative() {
        num = -num;
        vec.flip(0);
    }
    for (i, f) in fb.iter().enumerate() {
        let p = BigInt::from(f.prime);
        if !(&num % &p).is_zero() {
            continue;
        }
        let mut e = 1;
        num /= &p;
        while (&num % &p).is_zero() {
            e += 1;
            num /= &p;
        }
        if e % 2 == 1 {
            vec.flip(i + 1);
        }
    }
    (vec, num)
}

/// SymPy's `qs_factor(N, prime_bound, M, ERROR_TERM, seed)` as an
/// insertion-ordered list of `(factor, multiplicity)`.
pub fn qs_factor(n: &BigInt, prime_bound: u64, m: u64, error_term: f64, seed: &BigInt) -> Result<Vec<(BigInt, u32)>, String> {
    if n < &BigInt::from(2) {
        return Err("N should be greater than 1".into());
    }
    let mut factors: Vec<(BigInt, u32)> = vec![];
    let mut n = n.clone();
    if n.is_even() {
        let mut e = 1;
        n /= 2;
        while n.is_even() && !n.is_zero() {
            n /= 2;
            e += 1;
        }
        factors.push((BigInt::from(2), e));
    }
    if isprime(&n) {
        factors.push((n, 1));
        return Ok(factors);
    }
    if let Some((b, e)) = perfect_power_big(&n) {
        factors.push((b, e as u32));
        return Ok(factors);
    }
    if n <= BigInt::one() {
        return Err("N should not be a power of 2".into());
    }
    let mut n_copy = n.clone();
    let mut rng = PyRandom::new(seed);

    // factor base
    let mut fb: Vec<Fb> = vec![];
    let (mut idx_1000, mut idx_5000): (Option<usize>, Option<usize>) = (None, None);
    for p in primerange(1, prime_bound) {
        let bp = BigInt::from(p);
        if !n.modpow(&BigInt::from((p - 1) / 2), &bp).is_one() {
            continue;
        }
        if p > 1000 && idx_1000.is_none() {
            idx_1000 = Some(fb.len().wrapping_sub(1));
        }
        if p > 5000 && idx_5000.is_none() {
            idx_5000 = Some(fb.len().wrapping_sub(1));
        }
        let t = sqrt_mod_iter(&n, &bp).into_iter().min().unwrap().to_u64().unwrap();
        let log_p = ((p as f64).ln() * 1024.0).round() as i64;
        fb.push(Fb { prime: p, tmem_p: t, log_p, soln1: None, soln2: 0, b_ainv: vec![] });
    }
    if fb.is_empty() {
        return Err("list index out of range".into());
    }
    let width = fb.len() + 1;
    let threshold = fb.len() * 105 / 100;
    let mm = m as i64;
    let approx_val = ln(&(&n * 2)) / 2.0 - (m as f64).ln();
    let start = idx_1000.filter(|&i| i != 0).unwrap_or(0);
    let end = idx_5000.filter(|&i| i != 0).unwrap_or(fb.len() - 1);
    let accumulated = ((m as f64).ln() + ln(&n) / 2.0 - error_term) * 1024.0;
    let partial_bound = BigInt::from(128 * fb.last().unwrap().prime);

    let mut smooth: Vec<Relation> = vec![];
    let mut partial: HashMap<BigInt, Relation> = HashMap::new();
    let nb = n.clone();
    'outer: loop {
        // choose `a` as a product of factor base primes close to sqrt(2N)/M
        let (mut best_a, mut best_q, mut best_ratio): (BigInt, Vec<usize>, Option<f64>) = (BigInt::one(), vec![], None);
        for _ in 0..50 {
            let mut a = BigInt::one();
            let mut q: Vec<usize> = vec![];
            while ln(&a) < approx_val {
                let mut rand_p = 0usize;
                while rand_p == 0 || q.contains(&rand_p) {
                    rand_p = rng.randint(&BigInt::from(start), &BigInt::from(end)).to_usize().unwrap();
                }
                a *= fb[rand_p].prime;
                q.push(rand_p);
            }
            let ratio = (ln(&a) - approx_val).exp();
            if best_ratio.is_none_or(|r| (ratio - 1.0).abs() < (r - 1.0).abs()) {
                best_q = q;
                best_a = a;
                best_ratio = Some(ratio);
            }
        }
        let a = best_a;
        let mut bl: Vec<BigInt> = vec![];
        for &val in &best_q {
            let ql = fb[val].prime;
            let aq = &a / ql;
            let mut gamma = mulmod(fb[val].tmem_p, inv_mod(pmod(&aq, ql), ql), ql);
            if 2 * gamma > ql {
                gamma = ql - gamma;
            }
            bl.push(aq * gamma);
        }
        let mut bcur: BigInt = bl.iter().sum();
        for f in fb.iter_mut() {
            let p = f.prime;
            if pmod(&a, p) == 0 {
                f.soln1 = None;
                continue;
            }
            let a_inv = inv_mod(pmod(&a, p), p);
            f.b_ainv = bl.iter().map(|be| mulmod(pmod(&(be * 2), p), a_inv, p)).collect();
            let bm = pmod(&bcur, p);
            f.soln1 = Some(mulmod(a_inv, (f.tmem_p + p - bm) % p, p));
            f.soln2 = mulmod(a_inv, (2 * p - f.tmem_p - bm) % p, p);
        }
        // SymPy computes `range(1, 2**(len(B) - 1))`, a TypeError for `B == []`
        let npolys: u64 = if bl.is_empty() { 1 } else { 1 << (bl.len() - 1) };
        for i in 0..npolys {
            if i > 0 {
                let v = i.trailing_zeros() as usize;
                let neg_pow: i64 = 2 * ((i >> (v + 1)) % 2) as i64 - 1;
                bcur = &bcur + 2 * neg_pow * &bl[v];
                for f in fb.iter_mut() {
                    let Some(s1) = f.soln1 else { continue };
                    let p = f.prime as i128;
                    let d = neg_pow as i128 * f.b_ainv[v] as i128;
                    f.soln1 = Some((s1 as i128 - d).rem_euclid(p) as u64);
                    f.soln2 = (f.soln2 as i128 - d).rem_euclid(p) as u64;
                }
            }
            let a2 = &a * &a;
            let ab = 2 * &a * &bcur;
            let b2 = &bcur * &bcur - &nb;
            // sieve
            let mut sieve = vec![0i64; 2 * m as usize + 1];
            for f in &fb {
                let Some(s1) = f.soln1 else { continue };
                let p = f.prime as usize;
                let mut idx = ((m + s1) % f.prime) as usize;
                while idx < 2 * m as usize {
                    sieve[idx] += f.log_p;
                    idx += p;
                }
                if f.prime == 2 {
                    continue;
                }
                let mut idx = ((m + f.soln2) % f.prime) as usize;
                while idx < 2 * m as usize {
                    sieve[idx] += f.log_p;
                    idx += p;
                }
            }
            // trial division
            let mut proper: Vec<BigInt> = vec![];
            for (k, &val) in sieve.iter().enumerate() {
                if (val as f64) < accumulated {
                    continue;
                }
                let x = BigInt::from(k as i64 - mm);
                let v = (&a2 * &x + &ab) * &x + &b2;
                let (mut vec, num) = check_smoothness(&v, &fb, width);
                if num.is_one() {
                    smooth.push(Relation { u: &a * &x + &bcur, v, vec });
                } else if num < partial_bound && isprime(&num) {
                    if (&nb % &num).is_zero() {
                        if !proper.contains(&num) {
                            proper.push(num);
                        }
                        continue;
                    }
                    let mut u = &a * &x + &bcur;
                    if let Some(prev) = partial.remove(&num) {
                        u = (u * &prev.u * mod_inverse(&num, &nb).unwrap()).mod_floor(&nb);
                        let v = (v * &prev.v).div_floor(&(&num * &num));
                        vec.xor(&prev.vec);
                        smooth.push(Relation { u, v, vec });
                    } else {
                        partial.insert(num, Relation { u, v, vec });
                    }
                }
            }
            proper.sort();
            for p in proper {
                if !(&n_copy % &p).is_zero() {
                    continue;
                }
                let mut e = 1;
                n_copy /= &p;
                while (&n_copy % &p).is_zero() {
                    n_copy /= &p;
                    e += 1;
                }
                factors.push((p, e));
            }
            if threshold <= smooth.len() {
                break 'outer;
            }
        }
        if bl.is_empty() {
            return Err("'float' object cannot be interpreted as an integer".into());
        }
    }

    // GF(2) elimination exactly as SymPy's `_find_factor`
    let rows = smooth.len();
    let mut matrix: Vec<Bits> = smooth.iter().map(|r| r.vec.clone()).collect();
    let mut mark = vec![false; rows];
    for pos in 0..width {
        for i in 0..rows {
            if matrix[i].get(pos) {
                let mut add_col = matrix[i].clone();
                add_col.flip(pos);
                let mut unit = Bits::new(width);
                unit.flip(pos);
                matrix[i] = unit;
                mark[i] = true;
                for j in i + 1..rows {
                    if matrix[j].get(pos) {
                        matrix[j].xor(&add_col);
                    }
                }
                break;
            }
        }
    }
    for i in 0..rows {
        if mark[i] {
            continue;
        }
        let (mut u, mut v) = (smooth[i].u.clone(), smooth[i].v.clone());
        for j in 0..rows {
            if mark[j] && matrix[i].intersects(&matrix[j]) {
                u *= &smooth[j].u;
                v *= &smooth[j].v;
            }
        }
        if v.is_negative() {
            return Err("isqrt() argument must be nonnegative".into());
        }
        let v = v.sqrt();
        let g = (&u - &v).gcd(&nb);
        if g > BigInt::one() && g < nb && (&n_copy % &g).is_zero() {
            let mut e = 1;
            n_copy /= &g;
            while (&n_copy % &g).is_zero() {
                n_copy /= &g;
                e += 1;
            }
            factors.push((g, e));
            if n_copy.is_one() || isprime(&n_copy) {
                break;
            }
        }
    }
    if !n_copy.is_one() {
        factors.push((n_copy, 1));
    }
    Ok(factors)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn semiprime() {
        let n: BigInt = "25645121643901801".parse().unwrap();
        let f = qs_factor(&n, 2000, 10000, 25.0, &BigInt::from(1234)).unwrap();
        let mut ps: Vec<String> = f.iter().map(|(p, _)| p.to_string()).collect();
        ps.sort();
        assert_eq!(ps, vec!["4753701529", "5394769"]);
    }
}
