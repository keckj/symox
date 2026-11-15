//! Dense univariate polynomials over `GF(p)` for word sized primes
//! (`p < 2^32`, so that products fit in a `u64`), and their
//! factorisation (distinct degree + Cantor–Zassenhaus equal degree).

use super::int::Int;

pub type GfPoly = Vec<u64>;

#[inline]
pub fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    (a * b) % p
}

pub fn powmod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64 % p;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, p);
        }
        a = mulmod(a, a, p);
        e >>= 1;
    }
    r
}

pub fn invmod(a: u64, p: u64) -> u64 {
    // extended Euclid on signed values
    let (mut r0, mut r1) = (a as i64, p as i64);
    let (mut s0, mut s1) = (1i64, 0i64);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }
    debug_assert_eq!(r0, 1, "not invertible");
    s0.rem_euclid(p as i64) as u64
}

pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    // deterministic Miller-Rabin for 64 bit
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    let mpow = |mut a: u128, mut e: u64| -> u128 {
        let m = n as u128;
        let mut r = 1u128;
        a %= m;
        while e > 0 {
            if e & 1 == 1 {
                r = r * a % m;
            }
            a = a * a % m;
            e >>= 1;
        }
        r
    };
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = mpow(a as u128, d);
        if x == 1 || x == (n - 1) as u128 {
            continue;
        }
        for _ in 0..s - 1 {
            x = x * x % n as u128;
            if x == (n - 1) as u128 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

pub fn next_prime(mut n: u64) -> u64 {
    n += 1;
    while !is_prime(n) {
        n += 1;
    }
    n
}

pub fn trim(p: &mut GfPoly) {
    while p.last() == Some(&0) {
        p.pop();
    }
}

pub fn from_int(f: &[Int], p: u64) -> GfPoly {
    let mut r: GfPoly = f.iter().map(|c| c.rem_u64(p)).collect();
    trim(&mut r);
    r
}

/// Symmetric lift to integers.
pub fn to_int_sym(f: &[u64], p: u64) -> Vec<Int> {
    f.iter()
        .map(|&c| if c > p / 2 { Int::S(c as i64 - p as i64) } else { Int::S(c as i64) })
        .collect()
}

pub fn add(a: &[u64], b: &[u64], p: u64) -> GfPoly {
    let n = a.len().max(b.len());
    let mut r = Vec::with_capacity(n);
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0) + b.get(i).copied().unwrap_or(0);
        r.push(if x >= p { x - p } else { x });
    }
    trim(&mut r);
    r
}

pub fn sub(a: &[u64], b: &[u64], p: u64) -> GfPoly {
    let n = a.len().max(b.len());
    let mut r = Vec::with_capacity(n);
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0);
        let y = b.get(i).copied().unwrap_or(0);
        r.push(if x >= y { x - y } else { x + p - y });
    }
    trim(&mut r);
    r
}

pub fn scale(a: &[u64], k: u64, p: u64) -> GfPoly {
    let mut r: GfPoly = a.iter().map(|&c| mulmod(c, k, p)).collect();
    trim(&mut r);
    r
}

pub fn mul(a: &[u64], b: &[u64], p: u64) -> GfPoly {
    super::factor_budget::dense(a.len().saturating_add(b.len()).saturating_sub(1));
    super::factor_budget::charge(a.len().saturating_mul(b.len()) as u64);
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let n = a.len() + b.len() - 1;
    // accumulate in u128 to postpone reductions
    let mut acc = vec![0u128; n];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            acc[i + j] += (x * y) as u128;
        }
    }
    let mut r: GfPoly = acc.into_iter().map(|v| (v % p as u128) as u64).collect();
    trim(&mut r);
    r
}

pub fn monic(a: &[u64], p: u64) -> GfPoly {
    match a.last() {
        None => vec![],
        Some(&l) if l == 1 => a.to_vec(),
        Some(&l) => scale(a, invmod(l, p), p),
    }
}

pub fn divrem(a: &[u64], b: &[u64], p: u64) -> (GfPoly, GfPoly) {
    super::factor_budget::dense(a.len());
    super::factor_budget::charge(super::factor_budget::division_work(a.len(), b.len()));
    assert!(!b.is_empty());
    if a.len() < b.len() {
        return (vec![], a.to_vec());
    }
    let db = b.len() - 1;
    let inv = invmod(b[db], p);
    let mut q = vec![0u64; a.len() - db];
    // lazy reduction: subtract c*b_j as the addition of c*(p - b_j), only
    // the leading coefficient of each step is reduced
    let nb: Vec<u64> = b[..db].iter().map(|&x| if x == 0 { 0 } else { p - x }).collect();
    if p < (1 << 16) {
        let mut r: Vec<u64> = a.to_vec();
        for k in (0..q.len()).rev() {
            let c = mulmod(r[k + db] % p, inv, p);
            if c != 0 {
                for (j, &bj) in nb.iter().enumerate() {
                    r[k + j] += c * bj;
                }
            }
            q[k] = c;
        }
        let mut rr: GfPoly = r[..db].iter().map(|&x| x % p).collect();
        trim(&mut rr);
        trim(&mut q);
        return (q, rr);
    }
    let mut r: Vec<u128> = a.iter().map(|&x| x as u128).collect();
    let pp = p as u128;
    for k in (0..q.len()).rev() {
        let c = mulmod((r[k + db] % pp) as u64, inv, p);
        if c != 0 {
            for (j, &bj) in nb.iter().enumerate() {
                r[k + j] += (c * bj) as u128;
            }
        }
        q[k] = c;
    }
    let mut rr: GfPoly = r[..db].iter().map(|&x| (x % pp) as u64).collect();
    trim(&mut rr);
    trim(&mut q);
    (q, rr)
}

pub fn rem(a: &[u64], b: &[u64], p: u64) -> GfPoly {
    if a.len() < b.len() {
        return a.to_vec();
    }
    divrem(a, b, p).1
}

pub fn gcd(a: &[u64], b: &[u64], p: u64) -> GfPoly {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    while !b.is_empty() {
        let r = rem(&a, &b, p);
        a = b;
        b = r;
    }
    monic(&a, p)
}

/// `(g, s, t)` with `s a + t b = g` monic.
pub fn gcdex(a: &[u64], b: &[u64], p: u64) -> (GfPoly, GfPoly, GfPoly) {
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    let (mut s0, mut s1): (GfPoly, GfPoly) = (vec![1], vec![]);
    let (mut t0, mut t1): (GfPoly, GfPoly) = (vec![], vec![1]);
    while !r1.is_empty() {
        let (q, r) = divrem(&r0, &r1, p);
        r0 = std::mem::replace(&mut r1, r);
        let s2 = sub(&s0, &mul(&q, &s1, p), p);
        s0 = std::mem::replace(&mut s1, s2);
        let t2 = sub(&t0, &mul(&q, &t1, p), p);
        t0 = std::mem::replace(&mut t1, t2);
    }
    if r0.is_empty() {
        return (r0, s0, t0);
    }
    let l = invmod(*r0.last().unwrap(), p);
    (scale(&r0, l, p), scale(&s0, l, p), scale(&t0, l, p))
}

pub fn derivative(a: &[u64], p: u64) -> GfPoly {
    let mut r: GfPoly = a
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, &c)| mulmod(c, i as u64 % p, p))
        .collect();
    trim(&mut r);
    r
}

/// `a^e mod f`.
pub fn powmod_poly(a: &[u64], mut e: u128, f: &[u64], p: u64) -> GfPoly {
    let mut base = rem(a, f, p);
    let mut r: GfPoly = vec![1];
    while e > 0 {
        if e & 1 == 1 {
            r = rem(&mul(&r, &base, p), f, p);
        }
        e >>= 1;
        if e > 0 {
            base = rem(&mul(&base, &base, p), f, p);
        }
    }
    r
}

/// Is `f` square-free modulo `p` (and of the same degree)?
pub fn is_squarefree(f: &[u64], p: u64) -> bool {
    let d = derivative(f, p);
    if d.is_empty() {
        return false;
    }
    gcd(f, &d, p).len() == 1
}

/// Distinct degree factorisation of a monic square-free `f`:
/// `[(g, d)]` with `g` the product of the irreducible factors of degree `d`.
pub fn ddf(f: &[u64], p: u64) -> Vec<(GfPoly, usize)> {
    let mut out = Vec::new();
    let mut f = f.to_vec();
    let x: GfPoly = vec![0, 1];
    let mut h = x.clone();
    let mut d = 0usize;
    let n0 = f.len() - 1;
    // Frobenius matrix rows x^(i p) mod f (built by repeated shifts, cheap
    // for small p); h^p = Σ h_i x^(i p)
    let frob: Option<Vec<GfPoly>> = if p < (1 << 16) && (p as usize) < 4 * n0 * n0 {
        Some(frobenius_rows(&f, p))
    } else {
        None
    };
    while f.len() > 1 {
        d += 1;
        if 2 * d > f.len() - 1 {
            break;
        }
        h = match &frob {
            Some(q) => {
                let mut acc = vec![0u64; n0];
                for (i, &hi) in h.iter().enumerate() {
                    if hi == 0 {
                        continue;
                    }
                    for (k, &qk) in q[i].iter().enumerate() {
                        acc[k] += hi * qk;
                        if acc[k] >= (1u64 << 62) {
                            acc[k] %= p;
                        }
                    }
                }
                let mut r: GfPoly = acc.into_iter().map(|v| v % p).collect();
                trim(&mut r);
                rem(&r, &f, p)
            }
            None => powmod_poly(&h, p as u128, &f, p),
        };
        let g = gcd(&f, &sub(&h, &x, p), p);
        if g.len() > 1 {
            f = divrem(&f, &g, p).0;
            h = rem(&h, &f, p);
            out.push((g, d));
        }
    }
    if f.len() > 1 {
        let df = f.len() - 1;
        out.push((f, df));
    }
    out
}

/// `x^(i p) mod f` for `i < deg f`.
fn frobenius_rows(f: &[u64], p: u64) -> Vec<GfPoly> {
    let n = f.len() - 1;
    super::factor_budget::charge((n as u64).saturating_mul(n as u64).saturating_mul(p));
    let lc_inv = invmod(f[n], p);
    let fm: Vec<u64> = f.iter().map(|&c| mulmod(c, lc_inv, p)).collect();
    // multiply a reduced polynomial (length n) by x modulo f
    let shift = |v: &mut Vec<u64>| {
        let top = v[n - 1];
        for k in (1..n).rev() {
            v[k] = v[k - 1];
        }
        v[0] = 0;
        if top != 0 {
            for k in 0..n {
                let t = mulmod(top, fm[k], p);
                v[k] = if v[k] >= t { v[k] - t } else { v[k] + p - t };
            }
        }
    };
    let mut rows = Vec::with_capacity(n);
    let mut cur = vec![0u64; n];
    cur[0] = 1;
    rows.push(cur.clone());
    for _ in 1..n {
        for _ in 0..p {
            shift(&mut cur);
        }
        rows.push(cur.clone());
    }
    rows
}

/// Small deterministic pseudo random generator (xorshift).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1)
    }
    pub fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
}

/// Equal degree factorisation (Cantor–Zassenhaus, `p` odd): split the
/// monic product `f` of irreducibles of degree `d`.
pub fn edf(f: &[u64], d: usize, p: u64, rng: &mut Rng) -> Vec<GfPoly> {
    let n = f.len() - 1;
    if n == d {
        return vec![f.to_vec()];
    }
    if p == 2 {
        return edf_char2(f, d, rng);
    }
    loop {
        super::factor_budget::event("finite-field splitting attempts");
        let mut a: GfPoly = (0..n).map(|_| rng.next() % p).collect();
        trim(&mut a);
        if a.len() < 2 {
            continue;
        }
        let mut g = gcd(f, &a, p);
        if g.len() == 1 {
            // a^((p^d - 1)/2) = (a^(1 + p + ... + p^(d-1)))^((p-1)/2)
            let mut t = a.clone();
            let mut s = a.clone();
            for _ in 1..d {
                t = powmod_poly(&t, p as u128, f, p);
                s = rem(&mul(&s, &t, p), f, p);
            }
            let b = powmod_poly(&s, ((p - 1) / 2) as u128, f, p);
            g = gcd(f, &sub(&b, &[1], p), p);
        }
        if g.len() > 1 && g.len() < f.len() {
            let h = divrem(f, &g, p).0;
            let mut r = edf(&g, d, p, rng);
            r.extend(edf(&monic(&h, p), d, p, rng));
            return r;
        }
    }
}

fn edf_char2(f: &[u64], d: usize, rng: &mut Rng) -> Vec<GfPoly> {
    let n = f.len() - 1;
    if n == d {
        return vec![f.to_vec()];
    }
    loop {
        super::factor_budget::event("finite-field splitting attempts");
        let mut a: GfPoly = (0..n).map(|_| rng.next() % 2).collect();
        trim(&mut a);
        if a.len() < 2 {
            continue;
        }
        // trace map a + a^2 + ... + a^(2^(d-1))
        let mut t = a.clone();
        let mut s = a.clone();
        for _ in 1..d {
            t = rem(&mul(&t, &t, 2), f, 2);
            s = add(&s, &t, 2);
        }
        let g = gcd(f, &s, 2);
        if g.len() > 1 && g.len() < f.len() {
            let h = divrem(f, &g, 2).0;
            let mut r = edf_char2(&g, d, rng);
            r.extend(edf_char2(&h, d, rng));
            return r;
        }
    }
}

/// Factor a monic square-free polynomial into monic irreducibles.
pub fn factor_sqf(f: &[u64], p: u64) -> Vec<GfPoly> {
    let mut rng = Rng::new(0x5eed ^ p);
    let mut out = Vec::new();
    for (g, d) in ddf(f, p) {
        out.extend(edf(&g, d, p, &mut rng));
    }
    out.sort();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_mod_p() {
        // x^4 - 1 mod 5 = (x-1)(x-2)(x-3)(x-4)
        let f = vec![4u64, 0, 0, 0, 1];
        let fs = factor_sqf(&f, 5);
        assert_eq!(fs.len(), 4);
        let mut prod = vec![1u64];
        for g in &fs {
            prod = mul(&prod, g, 5);
        }
        assert_eq!(prod, f);
        // x^2 + 1 irreducible mod 3
        assert_eq!(factor_sqf(&[1, 0, 1], 3).len(), 1);
        // mod 2
        let f2 = mul(&[1, 1, 1], &[1, 1, 0, 1], 2);
        assert_eq!(factor_sqf(&f2, 2).len(), 2);
        assert!(is_prime(1_000_000_007));
        assert!(!is_prime(1_000_000_008));
    }
}

/// Square-free decomposition over `GF(p)`: `(lc, [(monic factor, k)])`.
pub fn sqf_list(f: &[u64], p: u64) -> (u64, Vec<(GfPoly, u32)>) {
    let mut f = f.to_vec();
    trim(&mut f);
    if f.is_empty() {
        return (0, vec![]);
    }
    let lc = *f.last().unwrap();
    f = monic(&f, p);
    let mut factors: Vec<(GfPoly, u32)> = Vec::new();
    if f.len() <= 1 {
        return (lc, factors);
    }
    let mut n = 1u32;
    loop {
        let mut sqf = false;
        let df = derivative(&f, p);
        if !df.is_empty() {
            let mut g = gcd(&f, &df, p);
            let mut h = divrem(&f, &g, p).0;
            let mut i = 1u32;
            while h.len() > 1 {
                let gg = gcd(&g, &h, p);
                let hh = divrem(&h, &gg, p).0;
                if hh.len() > 1 {
                    factors.push((monic(&hh, p), i * n));
                }
                g = divrem(&g, &gg, p).0;
                h = gg;
                i += 1;
            }
            if g.len() <= 1 {
                sqf = true;
            } else {
                f = g;
            }
        }
        if sqf {
            break;
        }
        // f is a p-th power: f(x) = g(x^p) = g(x)^p
        let r = p as usize;
        let d = (f.len() - 1) / r;
        f = (0..=d).map(|i| f[i * r]).collect();
        n *= p as u32;
        if f.len() <= 1 {
            break;
        }
    }
    (lc, factors)
}

/// Full factorisation over `GF(p)`: `(lc, [(monic irreducible, k)])`.
pub fn factor_list(f: &[u64], p: u64) -> (u64, Vec<(GfPoly, u32)>) {
    let (lc, parts) = sqf_list(f, p);
    let mut out = Vec::new();
    for (g, k) in parts {
        for h in factor_sqf(&g, p) {
            out.push((h, k));
        }
    }
    (lc, out)
}
