//! Univariate factorisation over `ZZ`: modular factorisation, multifactor
//! quadratic Hensel lifting and Zassenhaus recombination (with a constant
//! term test), plus the cyclotomic shortcut for `x^n ± 1`.

use rustc_hash::FxHashMap;

use super::dense;
use super::gf;
use super::int::Int;
use super::sqf::sqf_list_zx;

// ----------------------------------------------------------------------
// polynomials modulo an integer m
// ----------------------------------------------------------------------

fn reduce(p: &[Int], m: &Int) -> Vec<Int> {
    let mut r: Vec<Int> = p.iter().map(|c| c.modulo(m)).collect();
    dense::trim(&mut r);
    r
}

fn mul_m(a: &[Int], b: &[Int], m: &Int) -> Vec<Int> {
    reduce(&dense::mul_int(a, b), m)
}

fn add_m(a: &[Int], b: &[Int], m: &Int) -> Vec<Int> {
    reduce(&dense::add(a, b), m)
}

fn sub_m(a: &[Int], b: &[Int], m: &Int) -> Vec<Int> {
    reduce(&dense::sub(a, b), m)
}

/// Division by a monic polynomial modulo `m`.
fn divrem_monic_m(a: &[Int], b: &[Int], m: &Int) -> (Vec<Int>, Vec<Int>) {
    if a.len() < b.len() {
        return (vec![], a.to_vec());
    }
    let db = b.len() - 1;
    let mut r = a.to_vec();
    let mut q = vec![Int::zero(); a.len() - db];
    for k in (0..q.len()).rev() {
        let c = r[k + db].modulo(m);
        if !c.is_zero() {
            for (j, bj) in b.iter().enumerate() {
                r[k + j].sub_mul(&c, bj);
                r[k + j] = r[k + j].modulo(m);
            }
        }
        q[k] = c;
    }
    r.truncate(db);
    (reduce(&q, m), reduce(&r, m))
}

fn inv_mod(a: &Int, m: &Int) -> Int {
    let (g, s, _) = a.modulo(m).gcdex(m);
    debug_assert!(g.is_one());
    s.modulo(m)
}

/// One quadratic Hensel step (von zur Gathen–Gerhard 15.10): from
/// `f = g h mod m`, `s g + t h = 1 mod m` (`h` monic) to the same
/// relations modulo `m2` (any multiple of `m` dividing `m^2`).
fn hensel_step(
    f: &[Int],
    g: &[Int],
    h: &[Int],
    s: &[Int],
    t: &[Int],
    m2: &Int,
) -> (Vec<Int>, Vec<Int>, Vec<Int>, Vec<Int>) {
    let e = sub_m(f, &mul_m(g, h, m2), m2);
    let (q, r) = divrem_monic_m(&mul_m(s, &e, m2), h, m2);
    let g2 = add_m(&add_m(g, &mul_m(t, &e, m2), m2), &mul_m(&q, g, m2), m2);
    let h2 = add_m(h, &r, m2);
    let b = sub_m(&add_m(&mul_m(s, &g2, m2), &mul_m(t, &h2, m2), m2), &[Int::one()], m2);
    let (c, d) = divrem_monic_m(&mul_m(s, &b, m2), &h2, m2);
    let s2 = sub_m(s, &d, m2);
    let t2 = sub_m(&sub_m(t, &mul_m(t, &b, m2), m2), &mul_m(&c, &g2, m2), m2);
    (g2, h2, s2, t2)
}

fn gf_to_int(p: &[u64]) -> Vec<Int> {
    p.iter().map(|&c| Int::from_u64(c)).collect()
}

/// Lift the monic modular factors `fs` of `f` (`f = lc(f) Π fs mod p`) to
/// monic factors modulo `pl = p^l`.
fn hensel_lift(f: &[Int], fs: &[Vec<u64>], p: u64, pl: &Int) -> Vec<Vec<Int>> {
    if fs.len() == 1 {
        let l = inv_mod(f.last().unwrap(), pl);
        return vec![reduce(&dense::scale(f, &l), pl)];
    }
    let k = fs.len() / 2;
    let pm = Int::from_u64(p);
    let mut g0: Vec<u64> = vec![f.last().unwrap().rem_u64(p)];
    for a in &fs[..k] {
        g0 = gf::mul(&g0, a, p);
    }
    let mut h0: Vec<u64> = vec![1];
    for a in &fs[k..] {
        h0 = gf::mul(&h0, a, p);
    }
    let (_, s0, t0) = gf::gcdex(&g0, &h0, p);
    let (mut g, mut h, mut s, mut t) = (gf_to_int(&g0), gf_to_int(&h0), gf_to_int(&s0), gf_to_int(&t0));
    // exponent sequence 1 -> ... -> l with e_{i+1} <= 2 e_i
    let l = {
        let mut l = 0u32;
        let mut m = Int::one();
        while &m < pl {
            m = m.mul(&pm);
            l += 1;
        }
        l
    };
    let mut seq = vec![l];
    while *seq.last().unwrap() > 1 {
        let e = *seq.last().unwrap();
        seq.push(e.div_ceil(2));
    }
    seq.reverse();
    for &e in &seq[1..] {
        let m2 = pm.pow(e);
        let r = hensel_step(f, &g, &h, &s, &t, &m2);
        g = r.0;
        h = r.1;
        s = r.2;
        t = r.3;
    }
    let g = reduce(&g, pl);
    let h = reduce(&h, pl);
    let mut out = hensel_lift(&g, &fs[..k], p, pl);
    out.extend(hensel_lift(&h, &fs[k..], p, pl));
    out
}

// ----------------------------------------------------------------------
// Montgomery arithmetic modulo an odd M < 2^126 (fast Hensel lifting)
// ----------------------------------------------------------------------

#[inline]
fn mul_wide(a: u128, b: u128) -> (u128, u128) {
    let (a0, a1) = (a as u64 as u128, a >> 64);
    let (b0, b1) = (b as u64 as u128, b >> 64);
    let p00 = a0 * b0;
    let p01 = a0 * b1;
    let p10 = a1 * b0;
    let p11 = a1 * b1;
    let mid = (p00 >> 64) + (p01 as u64 as u128) + (p10 as u64 as u128);
    let lo = (p00 as u64 as u128) | (mid << 64);
    let hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
    (hi, lo)
}

/// Arithmetic in `Z/M` with Montgomery representatives (`R = 2^128`).
struct Mont {
    m: u128,
    minv: u128, // -M^{-1} mod 2^128
    r2: u128,   // R^2 mod M
    mi: Int,
}

impl Mont {
    fn new(m: &Int) -> Option<Mont> {
        let mb = m.to_big();
        if mb.bits() > 125 || !mb.bit(0) {
            return None;
        }
        let mu: u128 = num_traits::ToPrimitive::to_u128(&mb)?;
        let mut inv: u128 = mu;
        for _ in 0..7 {
            inv = inv.wrapping_mul(2u128.wrapping_sub(mu.wrapping_mul(inv)));
        }
        let r2b = (num_bigint::BigInt::from(1) << 256usize) % &mb;
        let r2 = num_traits::ToPrimitive::to_u128(&r2b)?;
        Some(Mont { m: mu, minv: inv.wrapping_neg(), r2, mi: m.clone() })
    }

    #[inline]
    fn redc(&self, hi: u128, lo: u128) -> u128 {
        let u = lo.wrapping_mul(self.minv);
        let (uh, ul) = mul_wide(u, self.m);
        let (_, carry) = lo.overflowing_add(ul);
        let mut t = hi + uh + carry as u128;
        if t >= self.m {
            t -= self.m;
        }
        t
    }

    #[inline]
    fn mul(&self, a: u128, b: u128) -> u128 {
        let (h, l) = mul_wide(a, b);
        self.redc(h, l)
    }

    #[inline]
    fn add(&self, a: u128, b: u128) -> u128 {
        let s = a + b;
        if s >= self.m {
            s - self.m
        } else {
            s
        }
    }

    #[inline]
    fn sub(&self, a: u128, b: u128) -> u128 {
        if a >= b {
            a - b
        } else {
            a + self.m - b
        }
    }

    fn to_mont(&self, x: &Int) -> u128 {
        let r = x.modulo(&self.mi);
        let v = num_traits::ToPrimitive::to_u128(&r.to_big()).unwrap();
        self.mul(v, self.r2)
    }

    fn one(&self) -> u128 {
        self.mul(1, self.r2)
    }

    fn trim(p: &mut Vec<u128>) {
        while p.last() == Some(&0) {
            p.pop();
        }
    }

    fn pmul(&self, a: &[u128], b: &[u128]) -> Vec<u128> {
        if a.is_empty() || b.is_empty() {
            return vec![];
        }
        let mut r = vec![0u128; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                r[i + j] = self.add(r[i + j], self.mul(x, y));
            }
        }
        Self::trim(&mut r);
        r
    }

    fn padd(&self, a: &[u128], b: &[u128]) -> Vec<u128> {
        let n = a.len().max(b.len());
        let mut r: Vec<u128> = (0..n)
            .map(|i| self.add(a.get(i).copied().unwrap_or(0), b.get(i).copied().unwrap_or(0)))
            .collect();
        Self::trim(&mut r);
        r
    }

    fn psub(&self, a: &[u128], b: &[u128]) -> Vec<u128> {
        let n = a.len().max(b.len());
        let mut r: Vec<u128> = (0..n)
            .map(|i| self.sub(a.get(i).copied().unwrap_or(0), b.get(i).copied().unwrap_or(0)))
            .collect();
        Self::trim(&mut r);
        r
    }

    /// Division by a monic polynomial.
    fn pdivrem(&self, a: &[u128], b: &[u128]) -> (Vec<u128>, Vec<u128>) {
        if a.len() < b.len() {
            return (vec![], a.to_vec());
        }
        let db = b.len() - 1;
        let mut r = a.to_vec();
        let mut q = vec![0u128; a.len() - db];
        for k in (0..q.len()).rev() {
            let c = r[k + db];
            if c != 0 {
                for (j, &bj) in b.iter().enumerate() {
                    r[k + j] = self.sub(r[k + j], self.mul(c, bj));
                }
            }
            q[k] = c;
        }
        r.truncate(db);
        Self::trim(&mut r);
        Self::trim(&mut q);
        (q, r)
    }

    fn step(&self, f: &[u128], g: &[u128], h: &[u128], s: &[u128], t: &[u128]) -> (Vec<u128>, Vec<u128>, Vec<u128>, Vec<u128>) {
        let e = self.psub(f, &self.pmul(g, h));
        let (q, r) = self.pdivrem(&self.pmul(s, &e), h);
        let g2 = self.padd(&self.padd(g, &self.pmul(t, &e)), &self.pmul(&q, g));
        let h2 = self.padd(h, &r);
        let one = [self.one()];
        let b = self.psub(&self.padd(&self.pmul(s, &g2), &self.pmul(t, &h2)), &one);
        let (c, d) = self.pdivrem(&self.pmul(s, &b), &h2);
        let s2 = self.psub(s, &d);
        let t2 = self.psub(&self.psub(t, &self.pmul(t, &b)), &self.pmul(&c, &g2));
        (g2, h2, s2, t2)
    }

}

/// Plain residues (`u128` in `[0, m)`) to Montgomery form and back.
fn to_m(ctx: &Mont, v: &[u128]) -> Vec<u128> {
    let mut r: Vec<u128> = v.iter().map(|&c| ctx.mul(c % ctx.m, ctx.r2)).collect();
    Mont::trim(&mut r);
    r
}

fn from_m(ctx: &Mont, v: &[u128]) -> Vec<u128> {
    let mut r: Vec<u128> = v.iter().map(|&c| ctx.redc(0, c)).collect();
    Mont::trim(&mut r);
    r
}

/// Montgomery version of [`hensel_lift`]: `ctxs` are the successive
/// moduli `p^{e_i}` (the last one is `p^l`), `f` holds residues modulo
/// `p^l`.
fn hensel_lift_mont(ctxs: &[Mont], f: &[u128], fs: &[Vec<u64>], p: u64) -> Vec<Vec<u128>> {
    let top = ctxs.last().unwrap();
    if fs.len() == 1 {
        let lc = Int::from_big(num_bigint::BigInt::from(*f.last().unwrap()));
        let l = inv_mod(&lc, &top.mi);
        let lm = top.to_mont(&l);
        let fm = to_m(top, f);
        let v: Vec<u128> = fm.iter().map(|&c| top.mul(c, lm)).collect();
        return vec![from_m(top, &v)];
    }
    let k = fs.len() / 2;
    let lcp = (*f.last().unwrap() % p as u128) as u64;
    let mut g0: Vec<u64> = vec![lcp];
    for a in &fs[..k] {
        g0 = gf::mul(&g0, a, p);
    }
    let mut h0: Vec<u64> = vec![1];
    for a in &fs[k..] {
        h0 = gf::mul(&h0, a, p);
    }
    let (_, s0, t0) = gf::gcdex(&g0, &h0, p);
    let cv = |a: &[u64]| -> Vec<u128> { a.iter().map(|&c| c as u128).collect() };
    let (mut g, mut h, mut s, mut t) = (cv(&g0), cv(&h0), cv(&s0), cv(&t0));
    for ctx in ctxs {
        let fm = to_m(ctx, f);
        let r = ctx.step(&fm, &to_m(ctx, &g), &to_m(ctx, &h), &to_m(ctx, &s), &to_m(ctx, &t));
        g = from_m(ctx, &r.0);
        h = from_m(ctx, &r.1);
        s = from_m(ctx, &r.2);
        t = from_m(ctx, &r.3);
    }
    let mut out = hensel_lift_mont(ctxs, &g, &fs[..k], p);
    out.extend(hensel_lift_mont(ctxs, &h, &fs[k..], p));
    out
}

/// Hensel lifting with the fastest available arithmetic.
fn hensel_lift_any(f: &[Int], fs: &[Vec<u64>], p: u64, pl: &Int, l: u32) -> Vec<Vec<Int>> {
    if pl.bits() <= 125 && p % 2 == 1 {
        let mut seq = vec![l.max(1)];
        while *seq.last().unwrap() > 1 {
            let e = *seq.last().unwrap();
            seq.push(e.div_ceil(2));
        }
        seq.reverse();
        let pm = Int::from_u64(p);
        let ctxs: Option<Vec<Mont>> = seq[1..].iter().map(|&e| Mont::new(&pm.pow(e))).collect();
        if let Some(ctxs) = ctxs.filter(|c| !c.is_empty()) {
            let fm: Vec<u128> = f
                .iter()
                .map(|c| num_traits::ToPrimitive::to_u128(&c.modulo(pl).to_big()).unwrap())
                .collect();
            let lifted = hensel_lift_mont(&ctxs, &fm, fs, p);
            return lifted
                .into_iter()
                .map(|v| {
                    let mut r: Vec<Int> = v.into_iter().map(|c| Int::from_big(num_bigint::BigInt::from(c))).collect();
                    dense::trim(&mut r);
                    r
                })
                .collect();
        }
    }
    hensel_lift(f, fs, p, pl)
}

// ----------------------------------------------------------------------
// cyclotomic polynomials
// ----------------------------------------------------------------------

fn divisors(n: u64) -> Vec<u64> {
    let mut d = Vec::new();
    let mut i = 1;
    while i * i <= n {
        if n % i == 0 {
            d.push(i);
            if i * i != n {
                d.push(n / i);
            }
        }
        i += 1;
    }
    d.sort();
    d
}

/// The `n`-th cyclotomic polynomial.
pub fn cyclotomic(n: u64, memo: &mut FxHashMap<u64, Vec<Int>>) -> Vec<Int> {
    if let Some(p) = memo.get(&n) {
        return p.clone();
    }
    let mut f = vec![Int::zero(); n as usize + 1];
    f[0] = Int::S(-1);
    f[n as usize] = Int::one();
    for d in divisors(n) {
        if d < n {
            let c = cyclotomic(d, memo);
            f = dense::div_exact_int(&f, &c).expect("cyclotomic division");
        }
    }
    memo.insert(n, f.clone());
    f
}

/// Factors of `x^n - 1` or `x^n + 1` (`None` for other polynomials).
pub fn cyclotomic_factor(f: &[Int]) -> Option<Vec<Vec<Int>>> {
    let n = f.len() as u64 - 1;
    if n < 1 || !f.last().unwrap().is_one() {
        return None;
    }
    if f[1..f.len() - 1].iter().any(|c| !c.is_zero()) {
        return None;
    }
    let mut memo = FxHashMap::default();
    if f[0].is_minus_one() {
        Some(divisors(n).into_iter().map(|d| cyclotomic(d, &mut memo)).collect())
    } else if f[0].is_one() {
        let ds = divisors(n);
        let mut out = Vec::new();
        for d in divisors(2 * n) {
            if !ds.contains(&d) {
                out.push(cyclotomic(d, &mut memo));
            }
        }
        Some(out)
    } else {
        None
    }
}

// ----------------------------------------------------------------------
// Zassenhaus
// ----------------------------------------------------------------------

fn subset_sums(degrees: &[usize], n: usize) -> Vec<bool> {
    let mut s = vec![false; n + 1];
    s[0] = true;
    for &d in degrees {
        for k in (d..=n).rev() {
            if s[k - d] {
                s[k] = true;
            }
        }
    }
    s
}

/// Irreducible factors of a primitive square-free polynomial with
/// positive leading coefficient and degree >= 1.
pub fn zassenhaus(f: &[Int]) -> Vec<Vec<Int>> {
    let n = f.len() - 1;
    if n <= 1 {
        return vec![f.to_vec()];
    }
    if f[0].is_zero() {
        let mut out = vec![vec![Int::zero(), Int::one()]];
        out.extend(zassenhaus(&f[1..]));
        return out;
    }
    let lc = f.last().unwrap().clone();
    // choose a prime: square-free reduction, few factors
    let mut best: Option<(u64, Vec<(Vec<u64>, usize)>, usize)> = None;
    let mut possible = vec![true; n + 1];
    let mut tried = 0;
    let mut p = 2u64;
    while tried < 7 {
        p = gf::next_prime(p);
        if lc.rem_u64(p) == 0 {
            continue;
        }
        let fp = gf::from_int(f, p);
        if fp.len() != f.len() || !gf::is_squarefree(&fp, p) {
            continue;
        }
        tried += 1;
        let fm = gf::monic(&fp, p);
        let dd = gf::ddf(&fm, p);
        let mut degs = Vec::new();
        for (g, d) in &dd {
            for _ in 0..(g.len() - 1) / d {
                degs.push(*d);
            }
        }
        let count = degs.len();
        let ss = subset_sums(&degs, n);
        for k in 0..=n {
            possible[k] &= ss[k];
        }
        if count == 1 || (1..n).all(|k| !possible[k]) {
            return vec![f.to_vec()];
        }
        if best.as_ref().is_none_or(|b| count < b.2) {
            best = Some((p, dd, count));
        }
        if p > 1000 && tried >= 3 {
            break;
        }
    }
    let (p, dd, _) = best.unwrap();
    let mut rng = gf::Rng::new(0xfac7 ^ p);
    let mut fs: Vec<Vec<u64>> = Vec::new();
    for (g, d) in dd {
        fs.extend(gf::edf(&g, d, p, &mut rng));
    }
    // coefficient bound for the factors (times the leading coefficient)
    let norm2: f64 = f.iter().map(|c| c.to_f64() * c.to_f64()).sum::<f64>().sqrt();
    let norm2_log = if norm2.is_finite() {
        norm2.log2()
    } else {
        f.iter().map(|c| c.log2()).fold(0.0, f64::max) + ((n + 1) as f64).log2() / 2.0
    };
    let bits = 1.0 + lc.log2() + n as f64 + norm2_log + 2.0;
    let l = (bits / (p as f64).log2()).ceil().max(1.0) as u32;
    let pl = Int::from_u64(p).pow(l);
    let lifted = hensel_lift_any(f, &fs, p, &pl, l);
    recombine(f, lifted, &pl)
}

fn sym(p: &[Int], m: &Int) -> Vec<Int> {
    let mut r: Vec<Int> = p.iter().map(|c| c.symmetric_mod(m)).collect();
    dense::trim(&mut r);
    r
}

fn next_comb(idx: &mut [usize], n: usize) -> bool {
    let s = idx.len();
    for k in (0..s).rev() {
        if idx[k] != k + n - s {
            idx[k] += 1;
            for j in k + 1..s {
                idx[j] = idx[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

fn recombine(f: &[Int], lifted: Vec<Vec<Int>>, pl: &Int) -> Vec<Vec<Int>> {
    let mut factors: Vec<Vec<Int>> = Vec::new();
    let mut fcur = f.to_vec();
    let mut t: Vec<usize> = (0..lifted.len()).collect();
    let mut s = 1usize;
    while 2 * s <= t.len() {
        let mut found: Option<(Vec<usize>, Vec<Int>, Vec<Int>)> = None;
        let b = fcur.last().unwrap().clone();
        let tc = b.mul(&fcur[0]);
        // iterate over the s-subsets of t
        let mut idx: Vec<usize> = (0..s).collect();
        loop {
            // constant term test
            let mut c = b.clone();
            for &i in &idx {
                c = c.mul(&lifted[t[i]][0]).modulo(pl);
            }
            let c = c.symmetric_mod(pl);
            if !c.is_zero() && tc.div_exact(&c).is_some() {
                let mut g = vec![b.clone()];
                for &i in &idx {
                    g = mul_m(&g, &lifted[t[i]], pl);
                }
                let g = dense::primitive(&sym(&g, pl)).1;
                if let Some(q) = dense::div_exact_int(&fcur, &g) {
                    found = Some((idx.clone(), g, q));
                    break;
                }
            }
            if !next_comb(&mut idx, t.len()) {
                break;
            }
        }
        match found {
            Some((sel, g, q)) => {
                factors.push(g);
                fcur = q;
                let sel_set: Vec<usize> = sel.iter().map(|&i| t[i]).collect();
                t.retain(|i| !sel_set.contains(i));
            }
            None => s += 1,
        }
    }
    if fcur.len() > 1 {
        let (_, pp) = dense::primitive(&fcur);
        factors.push(pp);
    }
    factors
}

/// Full factorisation of a univariate integer polynomial:
/// `(content, [(factor, multiplicity)])`, factors primitive with positive
/// leading coefficients (unsorted).
pub fn factor_zx(f: &[Int]) -> (Int, Vec<(Vec<Int>, u32)>) {
    if f.is_empty() {
        return (Int::zero(), vec![]);
    }
    let (c, mut g) = dense::primitive(f);
    if g.len() == 1 {
        return (c.mul(&g[0]), vec![]);
    }
    let mut out: Vec<(Vec<Int>, u32)> = Vec::new();
    // powers of x
    let j = g.iter().take_while(|x| x.is_zero()).count();
    if j > 0 {
        g.drain(..j);
        out.push((vec![Int::zero(), Int::one()], j as u32));
    }
    if g.len() == 1 {
        return (c, out);
    }
    if g.len() == 2 {
        out.push((g, 1));
        return (c, out);
    }
    if let Some(fs) = cyclotomic_factor(&g) {
        for h in fs {
            out.push((h, 1));
        }
        return (c, out);
    }
    let (_, parts) = sqf_list_zx(&g);
    for (a, k) in parts {
        let irr = if let Some(fs) = cyclotomic_factor(&a) { fs } else { zassenhaus(&a) };
        for h in irr {
            out.push((h, k));
        }
    }
    (c, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zi(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&x| Int::S(x)).collect()
    }

    fn prod(fs: &[(Vec<Int>, u32)]) -> Vec<Int> {
        let mut r = vec![Int::one()];
        for (f, k) in fs {
            for _ in 0..*k {
                r = dense::mul_int(&r, f);
            }
        }
        r
    }

    #[test]
    fn small() {
        // x^4 - 1
        let (c, fs) = factor_zx(&zi(&[-1, 0, 0, 0, 1]));
        assert_eq!(c, Int::one());
        assert_eq!(fs.len(), 3);
        // 6x^2 + 5x + 1 = (2x+1)(3x+1)
        let (c, fs) = factor_zx(&zi(&[1, 5, 6]));
        assert_eq!(c, Int::one());
        assert_eq!(fs.len(), 2);
        assert_eq!(prod(&fs), zi(&[1, 5, 6]));
        // x^2 + 1 irreducible
        let (_, fs) = factor_zx(&zi(&[1, 0, 1]));
        assert_eq!(fs.len(), 1);
        // (x^2-2)(x^2-3) needs recombination modulo most primes
        let f = dense::mul_int(&zi(&[-2, 0, 1]), &zi(&[-3, 0, 1]));
        let (_, fs) = factor_zx(&f);
        assert_eq!(fs.len(), 2);
        assert_eq!(prod(&fs), f);
    }

    #[test]
    fn bigger() {
        // (x+1)^3 (x^3 + 2x + 7)^2 (5x - 3)
        let f = prod(&[(zi(&[1, 1]), 3), (zi(&[7, 2, 0, 1]), 2), (zi(&[-3, 5]), 1)]);
        let (c, fs) = factor_zx(&dense::scale(&f, &Int::S(-6)));
        assert_eq!(c, Int::S(-6));
        assert_eq!(prod(&fs), f);
        assert_eq!(fs.len(), 3);
        // Wilkinson-like product of 12 linear factors
        let mut fs2 = Vec::new();
        for i in 1..=12 {
            fs2.push((zi(&[-i, 1]), 1));
        }
        let w = prod(&fs2);
        let (_, fs) = factor_zx(&w);
        assert_eq!(fs.len(), 12);
        let x20 = {
            let mut v = vec![Int::zero(); 21];
            v[0] = Int::S(-1);
            v[20] = Int::one();
            v
        };
        let (_, fs) = factor_zx(&x20);
        assert_eq!(fs.len(), 6);
    }
}

#[cfg(test)]
mod mont_tests {
    use super::*;

    #[test]
    fn mont_matches_int_lifting() {
        let mut f = vec![Int::one()];
        for i in 1..=8i64 {
            f = dense::mul_int(&f, &[Int::S(-i), Int::one()]);
        }
        let p = 11u64;
        let fp = gf::monic(&gf::from_int(&f, p), p);
        let fs = gf::factor_sqf(&fp, p);
        let l = 12u32;
        let pl = Int::from_u64(p).pow(l);
        let a = hensel_lift(&f, &fs, p, &pl);
        let b = hensel_lift_any(&f, &fs, p, &pl, l);
        assert_eq!(a, b);
        let ctx = Mont::new(&pl).unwrap();
        let x = Int::S(123456789);
        let y = Int::S(987654321);
        let r = ctx.from_mont(ctx.mul(ctx.to_mont(&x), ctx.to_mont(&y)));
        assert_eq!(r, x.mul(&y).modulo(&pl));
    }
}
