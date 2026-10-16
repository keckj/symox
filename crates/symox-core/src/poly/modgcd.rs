//! Brown's modular gcd for multivariate integer polynomials.
//!
//! The gcd is computed modulo word sized primes (`PGCD`: evaluation of the
//! last variable at points of `GF(p)`, recursive gcd, Newton interpolation),
//! the images are combined with the Chinese remainder theorem, and the
//! candidate is verified by exact division over `ZZ` (Geddes, Czapor,
//! Labahn, Algorithms 7.1/7.2).

use std::cell::Cell;
use std::fmt;

use super::gf;
use super::int::{Coeff, Field, Int, Rat};
use super::mpoly::MPoly;

thread_local! {
    static MODULUS: Cell<u64> = const { Cell::new(2147483647) };
}

#[inline]
fn modulus() -> u64 {
    MODULUS.with(|m| m.get())
}

/// An element of `GF(p)` for the current thread's modulus.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Zp(pub u64);

impl fmt::Debug for Zp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Coeff for Zp {
    fn zero() -> Self {
        Zp(0)
    }
    fn one() -> Self {
        Zp(1)
    }
    fn from_i64(i: i64) -> Self {
        Zp(i.rem_euclid(modulus() as i64) as u64)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    #[inline]
    fn is_one(&self) -> bool {
        self.0 == 1
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        let p = modulus();
        let s = self.0 + o.0;
        Zp(if s >= p { s - p } else { s })
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        let p = modulus();
        Zp(if self.0 >= o.0 { self.0 - o.0 } else { self.0 + p - o.0 })
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Zp((self.0 * o.0) % modulus())
    }
    #[inline]
    fn neg(&self) -> Self {
        if self.0 == 0 {
            *self
        } else {
            Zp(modulus() - self.0)
        }
    }
    fn is_negative(&self) -> bool {
        false
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        if o.0 == 0 {
            None
        } else {
            Some(self.mul(&o.inv()))
        }
    }
    fn to_rat(&self) -> Rat {
        Rat::from(Int::from_u64(self.0))
    }
}

impl Field for Zp {
    fn inv(&self) -> Self {
        Zp(gf::invmod(self.0, modulus()))
    }
}

fn reduce_poly(f: &MPoly<Int>, p: u64) -> MPoly<Zp> {
    f.map_coeffs(|c| Zp(c.rem_u64(p)))
}

/// Univariate gcd (in variable `v`) of polynomials over `GF(p)`, monic.
fn uni_gcd(a: &MPoly<Zp>, b: &MPoly<Zp>, v: usize) -> MPoly<Zp> {
    let p = modulus();
    let n = a.nvars;
    let da: Vec<u64> = super::dense::from_mpoly(a, v).iter().map(|c| c.0).collect();
    let db: Vec<u64> = super::dense::from_mpoly(b, v).iter().map(|c| c.0).collect();
    let g = gf::gcd(&da, &db, p);
    let gz: Vec<Zp> = g.iter().map(|&c| Zp(c)).collect();
    super::dense::to_mpoly(&gz, v, n)
}

/// Content of `f` seen as a polynomial in `vars` with coefficients in
/// `GF(p)[x_y]`: gcd of the coefficient polynomials (univariate in `y`).
fn content_y(f: &MPoly<Zp>, y: usize) -> MPoly<Zp> {
    let n = f.nvars;
    // group terms by the exponents of the variables other than y
    let mut index: rustc_hash::FxHashMap<Vec<u32>, usize> = Default::default();
    let mut groups: Vec<Vec<(Vec<u32>, Zp)>> = Vec::new();
    for (e, c) in f.iter() {
        let mut key = e.to_vec();
        key[y] = 0;
        let mut ye = vec![0u32; n];
        ye[y] = e[y];
        let id = *index.entry(key).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[id].push((ye, *c));
    }
    // small groups first: their gcd is cheaper and often 1
    groups.sort_by_key(|g| g.len());
    let mut g: Option<MPoly<Zp>> = None;
    for ts in groups {
        let c = MPoly::from_terms(n, ts);
        g = Some(match g {
            None => c.monic(),
            Some(h) => uni_gcd(&h, &c, y),
        });
        if g.as_ref().unwrap().is_one() {
            break;
        }
    }
    g.unwrap_or_else(|| MPoly::one(n))
}

/// Leading coefficient of `f` as a polynomial in `vars \ {y}` (lex order),
/// i.e. a univariate polynomial in `y`.
fn lc_y(f: &MPoly<Zp>, y: usize) -> MPoly<Zp> {
    let n = f.nvars;
    let mut key: Option<Vec<u32>> = None;
    let mut ts = Vec::new();
    for (e, c) in f.iter() {
        let mut k = e.to_vec();
        k[y] = 0;
        match &key {
            None => {
                key = Some(k.clone());
            }
            Some(kk) => {
                if *kk != k {
                    // terms are sorted with y possibly in between: compare
                    // lex on the other variables
                    if k > *kk {
                        key = Some(k.clone());
                        ts.clear();
                    } else {
                        continue;
                    }
                }
            }
        }
        let mut ye = vec![0u32; n];
        ye[y] = e[y];
        ts.push((ye, *c));
    }
    MPoly::from_terms(n, ts)
}

/// Leading monomial ignoring variable `y`.
fn lead_mono_wo(f: &MPoly<Zp>, y: usize) -> Vec<u32> {
    let mut best: Option<Vec<u32>> = None;
    for (e, _) in f.iter() {
        let mut k = e.to_vec();
        k[y] = 0;
        if best.as_ref().is_none_or(|b| k > *b) {
            best = Some(k);
        }
    }
    best.unwrap_or_default()
}

fn eval_y(f: &MPoly<Zp>, y: usize, a: Zp) -> MPoly<Zp> {
    f.eval(y, &a)
}

/// Divide `f` by the univariate (in `y`) polynomial `c`.
fn div_y(f: &MPoly<Zp>, c: &MPoly<Zp>) -> MPoly<Zp> {
    if c.is_one() {
        return f.clone();
    }
    f.div_exact(c).expect("content divides")
}

/// `PGCD`: gcd over `GF(p)`, normalised to be monic in lex order.
fn pgcd(a: &MPoly<Zp>, b: &MPoly<Zp>, depth: usize) -> Option<MPoly<Zp>> {
    let n = a.nvars;
    if a.is_zero() {
        return Some(b.monic());
    }
    if b.is_zero() {
        return Some(a.monic());
    }
    let da = a.degrees();
    let db = b.degrees();
    let vars: Vec<usize> = (0..n).filter(|&v| da[v] > 0 || db[v] > 0).collect();
    if vars.is_empty() {
        return Some(MPoly::one(n));
    }
    if vars.len() == 1 {
        return Some(uni_gcd(a, b, vars[0]));
    }
    if depth > 64 {
        return None;
    }
    let y = *vars.last().unwrap();
    // contents with respect to the other variables (univariate in y)
    let ca = content_y(a, y);
    let cb = content_y(b, y);
    let a = div_y(a, &ca);
    let b = div_y(b, &cb);
    let c = uni_gcd(&ca, &cb, y);
    let g = uni_gcd(&lc_y(&a, y), &lc_y(&b, y), y);
    let gd: Vec<Zp> = super::dense::from_mpoly(&g, y);
    let p = modulus();
    let bound = a.degree(y).min(b.degree(y)) as usize + g.degree(y).max(0) as usize + 1;
    let mut q: Option<MPoly<Zp>> = None;
    let mut qmono: Vec<u32> = Vec::new();
    let mut pts: Vec<Zp> = Vec::new();
    // Newton basis product Π (y - a_i)
    let mut basis: Vec<Zp> = vec![Zp(1)];
    let mut alpha = 0u64;
    let mut tries = 0usize;
    loop {
        let mut stable = false;
        alpha += 1;
        tries += 1;
        if alpha >= p || tries > 4 * bound + 64 {
            return None;
        }
        let av = Zp(alpha);
        let ga = super::dense::eval(&gd, &av);
        if ga.is_zero() {
            continue;
        }
        let ea = eval_y(&a, y, av);
        let eb = eval_y(&b, y, av);
        let Some(ci) = pgcd(&ea, &eb, depth + 1) else { continue };
        let ci = ci.scale(&ga);
        let mono = lead_mono_wo(&ci, y);
        if ci.is_constant() {
            // gcd is the content part only
            return Some(c.monic());
        }
        match &q {
            Some(_) if mono > qmono => continue, // unlucky point
            Some(_) if mono == qmono => {
                // Newton interpolation: q += (ci - q(alpha)) * basis / basis(alpha)
                let qv = eval_y(q.as_ref().unwrap(), y, av);
                let diff = ci.sub(&qv);
                if diff.is_zero() {
                    // early termination: the interpolant did not change
                    stable = true;
                } else {
                    let bv = super::dense::eval(&basis, &av);
                    let k = bv.inv();
                    let bp = super::dense::to_mpoly(&basis, y, n);
                    let qq = q.take().unwrap().add(&diff.mul(&bp).scale(&k));
                    q = Some(qq);
                }
            }
            _ => {
                // first point or a lower degree image: restart
                q = Some(ci);
                qmono = mono;
                pts.clear();
                basis = vec![Zp(1)];
            }
        }
        pts.push(av);
        basis = super::dense::mul(&basis, &[av.neg(), Zp(1)]);
        if pts.len() >= bound || stable {
            let cand = q.as_ref().unwrap();
            let cc = content_y(cand, y);
            let pp = div_y(cand, &cc);
            // inner levels rely on the degree bound; the final candidate
            // is verified over ZZ by the caller
            if depth > 0 || (a.div_exact(&pp).is_some() && b.div_exact(&pp).is_some()) {
                return Some(pp.mul(&c).monic());
            }
            if pts.len() > 2 * bound + 8 {
                return None;
            }
            if pts.len() < bound {
                // early termination failed the division test: go on
                continue;
            }
        }
    }
}

/// Modular gcd of two primitive integer polynomials (several variables).
/// `None` if the computation did not converge (the caller falls back).
pub fn modular_gcd_mz(a: &MPoly<Int>, b: &MPoly<Int>) -> Option<MPoly<Int>> {
    let n = a.nvars;
    let g = a.lc().gcd(&b.lc());
    let mut p: u64 = (1u64 << 31) - 1;
    let mut acc: Option<MPoly<Int>> = None; // CRT accumulator (coefficients in [0, m))
    let mut accmono: Vec<u32> = Vec::new();
    let mut m = Int::one();
    let old = modulus();
    let mut result = None;
    for _ in 0..400 {
        p = prev_prime(p);
        if g.rem_u64(p) == 0 || a.lc().rem_u64(p) == 0 || b.lc().rem_u64(p) == 0 {
            continue;
        }
        MODULUS.with(|c| c.set(p));
        let ap = reduce_poly(a, p);
        let bp = reduce_poly(b, p);
        let Some(gp) = pgcd(&ap, &bp, 0) else { continue };
        if gp.is_constant() {
            result = Some(MPoly::one(n));
            break;
        }
        // normalise the leading coefficient to g mod p
        let gp = gp.scale(&Zp(g.rem_u64(p)));
        let mono = gp.lm().to_vec();
        let gpi: MPoly<Int> = gp.map_coeffs(|c| Int::from_u64(c.0));
        match &acc {
            Some(_) if mono > accmono => continue,
            Some(cur) if mono == accmono => {
                // CRT: x = cur mod m, x = gp mod p
                let inv = Zp(m.rem_u64(p)).inv();
                let mp = m.mul(&Int::from_u64(p));
                let mut changed = false;
                let mut terms: Vec<(Vec<u32>, Int)> = Vec::new();
                // union of monomials
                let mut i = 0;
                let mut j = 0;
                while i < cur.len() || j < gpi.len() {
                    let (e, xc, yc) = if j >= gpi.len() || (i < cur.len() && cur.exp(i) > gpi.exp(j)) {
                        let r = (cur.exp(i).to_vec(), cur.coeff(i).clone(), 0u64);
                        i += 1;
                        r
                    } else if i >= cur.len() || gpi.exp(j) > cur.exp(i) {
                        let r = (gpi.exp(j).to_vec(), Int::zero(), gpi.coeff(j).as_i64().unwrap() as u64);
                        j += 1;
                        r
                    } else {
                        let r = (cur.exp(i).to_vec(), cur.coeff(i).clone(), gpi.coeff(j).as_i64().unwrap() as u64);
                        i += 1;
                        j += 1;
                        r
                    };
                    let xr = xc.rem_u64(p);
                    let delta = ((yc + p - xr) % p) * inv.0 % p;
                    let v = xc.add(&m.mul(&Int::from_u64(delta)));
                    if !changed && v.symmetric_mod(&mp) != xc.symmetric_mod(&m) {
                        changed = true;
                    }
                    terms.push((e, v));
                }
                m = mp;
                acc = Some(MPoly::from_terms(n, terms));
                if changed {
                    continue;
                }
            }
            _ => {
                acc = Some(gpi);
                accmono = mono;
                m = Int::from_u64(p);
                continue;
            }
        }
        // stable: test the candidate
        let cand = acc.as_ref().unwrap().map_coeffs(|c| c.symmetric_mod(&m));
        let (_, pp) = cand.primitive();
        if a.div_exact(&pp).is_some() && b.div_exact(&pp).is_some() {
            result = Some(pp);
            break;
        }
    }
    MODULUS.with(|c| c.set(old));
    result
}

fn prev_prime(mut p: u64) -> u64 {
    loop {
        p -= 1;
        if gf::is_prime(p) {
            return p;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modgcd() {
        let x = MPoly::<Int>::var(0, 3);
        let y = MPoly::<Int>::var(1, 3);
        let z = MPoly::<Int>::var(2, 3);
        let one = MPoly::<Int>::one(3);
        let g = x.add(&y).add(&z).add(&one).pow(3);
        let a = g.mul(&x.sub(&y).add(&z.scale(&Int::S(2))).pow(2));
        let b = g.mul(&x.mul(&y).sub(&z).add(&one.scale(&Int::S(2))));
        let h = modular_gcd_mz(&a, &b).unwrap();
        assert_eq!(h, g);
        // coprime
        let h = modular_gcd_mz(&x.add(&one), &y.add(&one)).unwrap();
        assert!(h.is_one());
        // big coefficients
        let big = Int::S(10).pow(30);
        let g2 = x.scale(&big).add(&y.mul(&z)).add(&one);
        let a2 = g2.mul(&x.add(&z));
        let b2 = g2.mul(&y.sub(&one));
        assert_eq!(modular_gcd_mz(&a2, &b2).unwrap(), g2);
    }
}
