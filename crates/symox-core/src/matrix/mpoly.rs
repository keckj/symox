//! A small sparse multivariate polynomial ring `Z[x_1, ..., x_n]` used by
//! the matrix algorithms (fraction-free elimination needs exact division,
//! `cancel` needs gcds).
//!
//! Terms are kept sorted in decreasing lexicographic order of their
//! exponent vectors (variable 0 is the most significant), without zero
//! coefficients. This module is self-contained on purpose; a full
//! polynomial module can replace it later.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use super::int::Int;

pub type Mono = SmallVec<[u32; 6]>;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct MPoly {
    pub nvars: usize,
    pub terms: Vec<(Mono, Int)>,
}

#[derive(PartialEq, Eq)]
struct LexMono(Mono);

impl Ord for LexMono {
    fn cmp(&self, o: &Self) -> Ordering {
        lex(&self.0, &o.0)
    }
}

impl PartialOrd for LexMono {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

#[inline]
fn lex(a: &Mono, b: &Mono) -> Ordering {
    a.as_slice().cmp(b.as_slice())
}

#[inline]
fn mono_mul(a: &Mono, b: &Mono) -> Mono {
    a.iter().zip(b.iter()).map(|(x, y)| crate::poly::mpoly::exponent_add(*x,*y)).collect()
}

#[inline]
fn mono_div(a: &Mono, b: &Mono) -> Option<Mono> {
    let mut r = Mono::with_capacity(a.len());
    for (x, y) in a.iter().zip(b.iter()) {
        if x < y {
            return None;
        }
        r.push(x - y);
    }
    Some(r)
}

impl MPoly {
    pub fn zero(nvars: usize) -> MPoly {
        MPoly { nvars, terms: Vec::new() }
    }

    pub fn constant(nvars: usize, c: Int) -> MPoly {
        if c.is_zero() {
            return MPoly::zero(nvars);
        }
        MPoly { nvars, terms: vec![(Mono::from_elem(0, nvars), c)] }
    }

    pub fn one(nvars: usize) -> MPoly {
        MPoly::constant(nvars, Int::ONE)
    }

    pub fn var(nvars: usize, i: usize) -> MPoly {
        let mut m = Mono::from_elem(0, nvars);
        m[i] = 1;
        MPoly { nvars, terms: vec![(m, Int::ONE)] }
    }

    pub fn monomial(nvars: usize, m: Mono, c: Int) -> MPoly {
        if c.is_zero() {
            return MPoly::zero(nvars);
        }
        MPoly { nvars, terms: vec![(m, c)] }
    }

    /// Build from unsorted terms (duplicates are combined).
    pub fn from_terms(nvars: usize, terms: Vec<(Mono, Int)>) -> MPoly {
        let mut map: FxHashMap<Mono, Int> = FxHashMap::default();
        for (m, c) in terms {
            match map.get_mut(&m) {
                Some(v) => *v = v.add(&c),
                None => {
                    map.insert(m, c);
                }
            }
        }
        let mut t: Vec<(Mono, Int)> = map.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        t.sort_unstable_by(|a, b| lex(&b.0, &a.0));
        MPoly { nvars, terms: t }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn is_one(&self) -> bool {
        self.terms.len() == 1 && self.terms[0].1.is_one() && self.terms[0].0.iter().all(|e| *e == 0)
    }

    pub fn is_constant(&self) -> bool {
        self.terms.is_empty() || (self.terms.len() == 1 && self.terms[0].0.iter().all(|e| *e == 0))
    }

    pub fn as_constant(&self) -> Option<Int> {
        if self.terms.is_empty() {
            return Some(Int::ZERO);
        }
        if self.is_constant() {
            return Some(self.terms[0].1.clone());
        }
        None
    }

    pub fn is_monomial(&self) -> bool {
        self.terms.len() == 1
    }

    /// Leading coefficient (lex).
    pub fn lc(&self) -> Int {
        self.terms.first().map(|t| t.1.clone()).unwrap_or(Int::ZERO)
    }

    pub fn neg(&self) -> MPoly {
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(m, c)| (m.clone(), c.neg())).collect(),
        }
    }

    pub fn add(&self, o: &MPoly) -> MPoly {
        self.merge(o, false)
    }

    pub fn sub(&self, o: &MPoly) -> MPoly {
        self.merge(o, true)
    }

    fn merge(&self, o: &MPoly, negate: bool) -> MPoly {
        if o.terms.is_empty() {
            return self.clone();
        }
        if self.terms.is_empty() {
            return if negate { o.neg() } else { o.clone() };
        }
        let mut out = Vec::with_capacity(self.terms.len() + o.terms.len());
        let (mut i, mut j) = (0, 0);
        let (a, b) = (&self.terms, &o.terms);
        while i < a.len() && j < b.len() {
            match lex(&a[i].0, &b[j].0) {
                Ordering::Greater => {
                    out.push(a[i].clone());
                    i += 1;
                }
                Ordering::Less => {
                    let c = if negate { b[j].1.neg() } else { b[j].1.clone() };
                    out.push((b[j].0.clone(), c));
                    j += 1;
                }
                Ordering::Equal => {
                    let c = if negate { a[i].1.sub(&b[j].1) } else { a[i].1.add(&b[j].1) };
                    if !c.is_zero() {
                        out.push((a[i].0.clone(), c));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend(a[i..].iter().cloned());
        for t in &b[j..] {
            let c = if negate { t.1.neg() } else { t.1.clone() };
            out.push((t.0.clone(), c));
        }
        MPoly { nvars: self.nvars, terms: out }
    }

    pub fn mul_int(&self, c: &Int) -> MPoly {
        if c.is_zero() {
            return MPoly::zero(self.nvars);
        }
        if c.is_one() {
            return self.clone();
        }
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(m, x)| (m.clone(), x.mul(c))).collect(),
        }
    }

    pub fn divexact_int(&self, c: &Int) -> MPoly {
        if c.is_one() {
            return self.clone();
        }
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(m, x)| (m.clone(), x.divexact(c))).collect(),
        }
    }

    pub fn mul_term(&self, m: &Mono, c: &Int) -> MPoly {
        if c.is_zero() {
            return MPoly::zero(self.nvars);
        }
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(x, k)| (mono_mul(x, m), k.mul(c))).collect(),
        }
    }

    pub fn mul(&self, o: &MPoly) -> MPoly {
        if self.is_zero() || o.is_zero() {
            return MPoly::zero(self.nvars);
        }
        if !self.product_fits(o) {
            return Self::wide_sum_products(&[(self,o)],self.nvars);
        }
        if self.terms.len() == 1 {
            return o.mul_term(&self.terms[0].0, &self.terms[0].1);
        }
        if o.terms.len() == 1 {
            return self.mul_term(&o.terms[0].0, &o.terms[0].1);
        }
        if self.terms.len() * o.terms.len() <= 256 {
            // small products: sort and merge instead of hashing
            let mut prods: Vec<(Mono, Int)> = Vec::with_capacity(self.terms.len() * o.terms.len());
            for (ma, ca) in &self.terms {
                for (mb, cb) in &o.terms {
                    prods.push((mono_mul(ma, mb), ca.mul(cb)));
                }
            }
            prods.sort_unstable_by(|a, b| lex(&b.0, &a.0));
            let mut out: Vec<(Mono, Int)> = Vec::with_capacity(prods.len());
            for (m, c) in prods {
                match out.last_mut() {
                    Some(last) if last.0 == m => last.1 = last.1.add(&c),
                    _ => out.push((m, c)),
                }
            }
            out.retain(|(_, c)| !c.is_zero());
            return MPoly { nvars: self.nvars, terms: out };
        }
        let mut map: FxHashMap<Mono, Int> =
            FxHashMap::with_capacity_and_hasher(self.terms.len() * o.terms.len(), Default::default());
        for (ma, ca) in &self.terms {
            for (mb, cb) in &o.terms {
                let m = mono_mul(ma, mb);
                let p = ca.mul(cb);
                match map.get_mut(&m) {
                    Some(v) => *v = v.add(&p),
                    None => {
                        map.insert(m, p);
                    }
                }
            }
        }
        let mut t: Vec<(Mono, Int)> = map.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        t.sort_unstable_by(|a, b| lex(&b.0, &a.0));
        MPoly { nvars: self.nvars, terms: t }
    }

    /// `sum a_i * b_i`, accumulated in one table instead of a chain of merges.
    pub fn sum_products(pairs: &[(&MPoly, &MPoly)], nvars: usize) -> MPoly {
        if pairs.iter().any(|(a,b)|!a.product_fits(b)) {
            return Self::wide_sum_products(pairs,nvars);
        }
        let size: usize = pairs.iter().map(|(a, b)| a.terms.len() * b.terms.len()).sum();
        if pairs.len() < 2 || size <= 256 {
            let mut acc = MPoly::zero(nvars);
            for (a, b) in pairs {
                acc = acc.add(&a.mul(b));
            }
            return acc;
        }
        let mut map: FxHashMap<Mono, Int> = FxHashMap::with_capacity_and_hasher(size, Default::default());
        for (a, b) in pairs {
            for (ma, ca) in &a.terms {
                for (mb, cb) in &b.terms {
                    let m = mono_mul(ma, mb);
                    let p = ca.mul(cb);
                    match map.get_mut(&m) {
                        Some(v) => *v = v.add(&p),
                        None => {
                            map.insert(m, p);
                        }
                    }
                }
            }
        }
        let mut t: Vec<(Mono, Int)> = map.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        t.sort_unstable_by(|a, b| lex(&b.0, &a.0));
        MPoly { nvars, terms: t }
    }

    fn product_fits(&self, other:&Self)->bool {
        self.is_zero() || other.is_zero() || (0..self.nvars).all(|v|self.degree(v).checked_add(other.degree(v)).is_some())
    }

    /// Matrix dot products can cancel high terms across several products.
    /// Keep all those products in one wide sparse accumulator before checking.
    fn wide_sum_products(pairs:&[(&MPoly,&MPoly)], nvars:usize)->MPoly {
        let mut accumulated:FxHashMap<Vec<u64>,Int>=FxHashMap::default();
        for (a,b) in pairs {
            for (ma,ca) in &a.terms {
                for (mb,cb) in &b.terms {
                    let exponent=ma.iter().zip(mb).map(|(a,b)|u64::from(*a)+u64::from(*b)).collect();
                    let product=ca.mul(cb);
                    if product.is_zero() { continue; }
                    let coefficient=accumulated.entry(exponent).or_insert(Int::ZERO);
                    *coefficient=coefficient.add(&product);
                }
            }
        }
        let terms=accumulated.into_iter().filter(|(_,c)|!c.is_zero()).map(|(m,c)| {
            let exponent=m.into_iter().map(|e|u32::try_from(e).unwrap_or_else(|_|crate::poly::mpoly::exponent_overflow())).collect();
            (exponent,c)
        }).collect();
        MPoly::from_terms(nvars,terms)
    }

    pub fn pow(&self, n: u32) -> MPoly {
        let mut r = MPoly::one(self.nvars);
        let mut b = self.clone();
        let mut n = n;
        while n > 0 {
            if n & 1 == 1 {
                r = r.mul(&b);
            }
            n >>= 1;
            if n > 0 {
                b = b.mul(&b);
            }
        }
        r
    }

    /// `self * a - o * b`.
    pub fn mul_sub(&self, a: &MPoly, o: &MPoly, b: &MPoly) -> MPoly {
        if !self.product_fits(a) || !o.product_fits(b) {
            return Self::wide_sum_products(&[(self,a),(&o.neg(),b)],self.nvars);
        }
        self.mul(a).sub(&o.mul(b))
    }

    /// Exact division; `None` when `o` does not divide `self`.
    pub fn divexact(&self, o: &MPoly) -> Option<MPoly> {
        if o.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some(MPoly::zero(self.nvars));
        }
        if let Some(c) = o.as_constant() {
            if c.is_one() {
                return Some(self.clone());
            }
            let mut terms = Vec::with_capacity(self.terms.len());
            for (m, x) in &self.terms {
                let (q, r) = x.div_rem(&c);
                if !r.is_zero() {
                    return None;
                }
                terms.push((m.clone(), q));
            }
            return Some(MPoly { nvars: self.nvars, terms });
        }
        if o.terms.len() == 1 {
            let (om, oc) = &o.terms[0];
            let mut terms = Vec::with_capacity(self.terms.len());
            for (m, x) in &self.terms {
                let (q, r) = x.div_rem(oc);
                if !r.is_zero() {
                    return None;
                }
                terms.push((mono_div(m, om)?, q));
            }
            return Some(MPoly { nvars: self.nvars, terms });
        }
        let (lm, lc) = (&o.terms[0].0, &o.terms[0].1);
        // quick degree checks
        for v in 0..self.nvars {
            if self.degree(v) < o.degree(v) {
                return None;
            }
        }
        // heap based division: the remainder lives in a hash map, the
        // candidate leading monomials in a max-heap
        let mut rem: FxHashMap<Mono, Int> =
            FxHashMap::with_capacity_and_hasher(self.terms.len() * 2, Default::default());
        let mut heap: BinaryHeap<LexMono> = BinaryHeap::with_capacity(self.terms.len() * 2);
        for (m, c) in &self.terms {
            rem.insert(m.clone(), c.clone());
            heap.push(LexMono(m.clone()));
        }
        let mut q: Vec<(Mono, Int)> = Vec::new();
        while let Some(LexMono(m)) = heap.pop() {
            let Some(c) = rem.remove(&m) else { continue };
            if c.is_zero() {
                continue;
            }
            let qm = mono_div(&m, lm)?;
            let (qc, r) = c.div_rem(lc);
            if !r.is_zero() {
                return None;
            }
            for (gm, gc) in &o.terms[1..] {
                let key = mono_mul(&qm, gm);
                let v = qc.mul(gc);
                match rem.get_mut(&key) {
                    Some(x) => *x = x.sub(&v),
                    None => {
                        rem.insert(key.clone(), v.neg());
                        heap.push(LexMono(key));
                    }
                }
            }
            q.push((qm, qc));
        }
        Some(MPoly { nvars: self.nvars, terms: q })
    }

    /// Degree in variable `v` (-1 for zero is represented as 0).
    pub fn degree(&self, v: usize) -> u32 {
        self.terms.iter().map(|(m, _)| m[v]).max().unwrap_or(0)
    }

    pub fn total_degree(&self) -> u64 {
        self.terms.iter().map(|(m, _)| m.iter().map(|&k| u64::from(k)).sum::<u64>()).max().unwrap_or(0)
    }

    pub fn uses_var(&self, v: usize) -> bool {
        self.terms.iter().any(|(m, _)| m[v] > 0)
    }

    /// gcd of the integer coefficients (non-negative).
    pub fn content(&self) -> Int {
        let mut g = Int::ZERO;
        for (_, c) in &self.terms {
            g = g.gcd(c);
            if g.is_one() {
                break;
            }
        }
        g
    }

    /// Component-wise minimum exponent vector.
    fn min_mono(&self) -> Mono {
        let mut m = self.terms[0].0.clone();
        for (t, _) in &self.terms[1..] {
            for (a, b) in m.iter_mut().zip(t.iter()) {
                if *b < *a {
                    *a = *b;
                }
            }
        }
        m
    }

    fn div_mono(&self, m: &Mono) -> MPoly {
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(x, c)| (mono_div(x, m).unwrap(), c.clone())).collect(),
        }
    }

    /// Split into coefficients with respect to variable `v`:
    /// `self = Σ c_k v^k`, returned dense, index = k.
    pub fn coeffs_in(&self, v: usize) -> Vec<MPoly> {
        let d = self.degree(v) as usize;
        let mut buckets: Vec<Vec<(Mono, Int)>> = vec![Vec::new(); d + 1];
        for (m, c) in &self.terms {
            let k = m[v] as usize;
            let mut mm = m.clone();
            mm[v] = 0;
            buckets[k].push((mm, c.clone()));
        }
        // terms stay sorted inside each bucket since removing the same
        // variable preserves the relative lex order... except when v is
        // not the most significant variable; re-sort to be safe.
        buckets
            .into_iter()
            .map(|mut t| {
                t.sort_unstable_by(|a, b| lex(&b.0, &a.0));
                MPoly { nvars: self.nvars, terms: t }
            })
            .collect()
    }

    fn from_coeffs_in(nvars: usize, v: usize, cs: &[MPoly]) -> MPoly {
        let mut terms = Vec::new();
        for (k, c) in cs.iter().enumerate() {
            for (m, x) in &c.terms {
                let mut mm = m.clone();
                let k=u32::try_from(k).unwrap_or_else(|_|crate::poly::mpoly::exponent_overflow());
                mm[v] = crate::poly::mpoly::exponent_add(mm[v],k);
                terms.push((mm, x.clone()));
            }
        }
        terms.sort_unstable_by(|a, b| lex(&b.0, &a.0));
        MPoly { nvars, terms }
    }

    /// Make the leading coefficient positive.
    pub fn normalize_sign(self) -> MPoly {
        if self.lc().is_negative() {
            self.neg()
        } else {
            self
        }
    }

    /// Evaluate variable `v` at integer `a`.
    pub fn eval_var(&self, v: usize, a: &Int) -> MPoly {
        let mut terms = Vec::with_capacity(self.terms.len());
        for (m, c) in &self.terms {
            let mut mm = m.clone();
            let e = mm[v];
            mm[v] = 0;
            terms.push((mm, c.mul(&a.pow(e))));
        }
        MPoly::from_terms(self.nvars, terms)
    }
}

// ----------------------------------------------------------------------
// gcd
// ----------------------------------------------------------------------

/// Greatest common divisor with positive leading coefficient.
pub fn gcd(f: &MPoly, g: &MPoly) -> MPoly {
    if f.is_zero() {
        return g.clone().normalize_sign();
    }
    if g.is_zero() {
        return f.clone().normalize_sign();
    }
    let n = f.nvars;
    let cf = f.content();
    let cg = g.content();
    let c = cf.gcd(&cg);
    if f.is_constant() || g.is_constant() {
        return MPoly::constant(n, c);
    }
    let f = f.divexact_int(&cf);
    let g = g.divexact_int(&cg);
    // common monomial factor
    let mf = f.min_mono();
    let mg = g.min_mono();
    let m: Mono = mf.iter().zip(mg.iter()).map(|(a, b)| *a.min(b)).collect();
    let f = f.div_mono(&mf);
    let g = g.div_mono(&mg);
    let h = gcd_prim(&f, &g);
    h.mul_term(&m, &c).normalize_sign()
}

// ----------------------------------------------------------------------
// heuristic gcd (Char, Geddes, Gonnet), as SymPy's dmp_zz_heu_gcd
// ----------------------------------------------------------------------

fn max_norm(f: &MPoly) -> Int {
    let mut m = Int::ZERO;
    for (_, c) in &f.terms {
        let a = c.abs();
        if a > m {
            m = a;
        }
    }
    m
}

fn isqrt(x: &Int) -> Int {
    Int::from_big(x.to_big().sqrt())
}

/// Recover a polynomial in `v` from its value at `v = x` (symmetric
/// base-`x` digits of every coefficient).
fn interpolate(h: &MPoly, x: &Int, v: usize) -> MPoly {
    let half = x.divexact(&Int::from(2));
    let mut terms = Vec::new();
    for (m, c) in &h.terms {
        let mut c = c.clone();
        let mut k = 0u32;
        while !c.is_zero() {
            let (_, mut d) = c.div_rem(x);
            if d.is_negative() {
                d = d.add(x);
            }
            if d > half {
                d = d.sub(x);
            }
            if !d.is_zero() {
                let mut mm = m.clone();
                mm[v] = crate::poly::mpoly::exponent_add(mm[v],k);
                terms.push((mm, d.clone()));
            }
            c = c.sub(&d).divexact(x);
            k += 1;
        }
    }
    MPoly::from_terms(h.nvars, terms)
}

fn primitive_part(f: &MPoly) -> MPoly {
    let c = f.content();
    if c.is_zero() || c.is_one() {
        f.clone()
    } else {
        f.divexact_int(&c)
    }
}

fn heu_gcd(f: &MPoly, g: &MPoly) -> Option<MPoly> {
    let vars: Vec<usize> = (0..f.nvars).filter(|&v| f.uses_var(v) || g.uses_var(v)).collect();
    heu_rec(f, g, &vars, 0)
}

fn heu_rec(f: &MPoly, g: &MPoly, vars: &[usize], depth: usize) -> Option<MPoly> {
    let n = f.nvars;
    if f.is_zero() || g.is_zero() {
        return None;
    }
    if f.is_constant() || g.is_constant() || vars.is_empty() {
        return Some(MPoly::constant(n, f.content().gcd(&g.content())));
    }
    if depth > 12 {
        return None;
    }
    let c = f.content().gcd(&g.content());
    let (f, g) = if c.is_one() { (f.clone(), g.clone()) } else { (f.divexact_int(&c), g.divexact_int(&c)) };
    let (fnorm, gnorm) = (max_norm(&f), max_norm(&g));
    let b = fnorm.clone().min(gnorm.clone()).mul(&Int::from(2)).add(&Int::from(29));
    let lcf = f.lc().abs();
    let lcg = g.lc().abs();
    let t = fnorm.divexact_floor(&lcf).min(gnorm.divexact_floor(&lcg)).mul(&Int::from(2)).add(&Int::from(4));
    let s99 = isqrt(&b).mul(&Int::from(99));
    let mut x = b.min(s99).max(t);
    let v0 = vars[0];
    for _ in 0..6 {
        let ff = f.eval_var(v0, &x);
        let gg = g.eval_var(v0, &x);
        if !ff.is_zero() && !gg.is_zero() {
            if let Some(h) = heu_rec(&ff, &gg, &vars[1..], depth + 1) {
                let hh = primitive_part(&interpolate(&h, &x, v0));
                if !hh.is_zero() {
                    if let (Some(_), Some(_)) = (f.divexact(&hh), g.divexact(&hh)) {
                        return Some(hh.mul_int(&c));
                    }
                }
                // cofactor based candidates
                if let Some(cff) = ff.divexact(&h) {
                    let cf = interpolate(&cff, &x, v0);
                    if !cf.is_zero() {
                        if let Some(hh) = f.divexact(&cf) {
                            if g.divexact(&hh).is_some() {
                                return Some(hh.mul_int(&c));
                            }
                        }
                    }
                }
                if let Some(cfg) = gg.divexact(&h) {
                    let cg = interpolate(&cfg, &x, v0);
                    if !cg.is_zero() {
                        if let Some(hh) = g.divexact(&cg) {
                            if f.divexact(&hh).is_some() {
                                return Some(hh.mul_int(&c));
                            }
                        }
                    }
                }
            }
        }
        // next evaluation point: 73794*x*sqrt(sqrt(x))//27011
        x = x.mul(&Int::from(73794)).mul(&isqrt(&isqrt(&x))).divexact_floor(&Int::from(27011));
    }
    None
}

/// gcd of two polynomials with unit content and no monomial content
/// (the result is only defined up to sign).
fn gcd_prim(f: &MPoly, g: &MPoly) -> MPoly {
    let n = f.nvars;
    if f.is_constant() || g.is_constant() {
        return MPoly::one(n);
    }
    if f == g {
        return f.clone();
    }
    if f.terms.len() == 1 || g.terms.len() == 1 {
        // a monomial has no non-monomial factors and monomial content
        // was removed
        return MPoly::one(n);
    }
    // variables used by only one of the two polynomials
    for v in 0..n {
        let (uf, ug) = (f.uses_var(v), g.uses_var(v));
        if uf && !ug {
            let c = content_in(f, v);
            return gcd(&c, g);
        }
        if ug && !uf {
            let c = content_in(g, v);
            return gcd(f, &c);
        }
    }
    if let Some(h) = heu_gcd(f, g) {
        return h;
    }
    // try trivial divisibility
    if f.terms.len() >= g.terms.len() {
        if f.divexact(g).is_some() {
            return g.clone();
        }
    } else if g.divexact(f).is_some() {
        return f.clone();
    }
    // main variable: the used variable of smallest degree
    let mut best: Option<(u32, usize)> = None;
    for v in 0..n {
        if f.uses_var(v) {
            let d = f.degree(v).max(g.degree(v));
            if best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, v));
            }
        }
    }
    let v = match best {
        Some((_, v)) => v,
        None => return MPoly::one(n),
    };
    let fc = f.coeffs_in(v);
    let gc = g.coeffs_in(v);
    let cont_f = gcd_list(&fc, n);
    let cont_g = gcd_list(&gc, n);
    let cc = gcd(&cont_f, &cont_g);
    let pf: Vec<MPoly> = fc.iter().map(|c| c.divexact(&cont_f).unwrap()).collect();
    let pg: Vec<MPoly> = gc.iter().map(|c| c.divexact(&cont_g).unwrap()).collect();
    let h = subresultant_gcd(pf, pg, n);
    let hp = MPoly::from_coeffs_in(n, v, &h);
    cc.mul(&hp)
}

fn gcd_list(cs: &[MPoly], n: usize) -> MPoly {
    let mut g = MPoly::zero(n);
    for c in cs {
        if c.is_zero() {
            continue;
        }
        g = gcd(&g, c);
        if g.is_one() {
            break;
        }
    }
    g
}

/// Content of `f` with respect to variable `v` (gcd of the coefficients).
pub fn content_in(f: &MPoly, v: usize) -> MPoly {
    gcd_list(&f.coeffs_in(v), f.nvars)
}

// univariate (in the main variable) helpers over the coefficient ring

fn trim(a: &mut Vec<MPoly>) {
    while a.last().is_some_and(|c| c.is_zero()) {
        a.pop();
    }
}

fn prem(a: &[MPoly], b: &[MPoly], n: usize) -> Vec<MPoly> {
    let db = b.len() - 1;
    let lcb = &b[db];
    let mut r: Vec<MPoly> = a.to_vec();
    trim(&mut r);
    if r.len() < b.len() {
        return r;
    }
    let mut e = (r.len() - b.len() + 1) as u32;
    while !r.is_empty() && r.len() > db {
        let dr = r.len() - 1;
        let t = r[dr].clone();
        let k = dr - db;
        // r = lcb*r - t*x^k*b
        for c in r.iter_mut() {
            *c = c.mul(lcb);
        }
        for (i, bc) in b.iter().enumerate() {
            let s = t.mul(bc);
            r[i + k] = r[i + k].sub(&s);
        }
        trim(&mut r);
        e -= 1;
    }
    if e > 0 {
        let f = lcb.pow(e);
        for c in r.iter_mut() {
            *c = c.mul(&f);
        }
    }
    let _ = n;
    r
}

fn subresultant_gcd(a: Vec<MPoly>, b: Vec<MPoly>, n: usize) -> Vec<MPoly> {
    let (mut a, mut b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    trim(&mut a);
    trim(&mut b);
    let mut g = MPoly::one(n);
    let mut h = MPoly::one(n);
    loop {
        let delta = (a.len() - b.len()) as u32;
        let r = prem(&a, &b, n);
        if r.is_empty() {
            break;
        }
        if r.len() == 1 {
            // constant remainder: coprime
            return vec![MPoly::one(n)];
        }
        a = b;
        let d = g.mul(&h.pow(delta));
        b = r.iter().map(|c| c.divexact(&d).expect("subresultant division")).collect();
        g = a.last().unwrap().clone();
        h = if delta == 0 {
            h
        } else if delta == 1 {
            g.clone()
        } else {
            g.pow(delta).divexact(&h.pow(delta - 1)).expect("subresultant h")
        };
    }
    // primitive part of b
    let c = gcd_list(&b, n);
    b.iter().map(|x| x.divexact(&c).unwrap()).collect()
}

/// `(gcd, f/gcd, g/gcd)`.
pub fn cofactors(f: &MPoly, g: &MPoly) -> (MPoly, MPoly, MPoly) {
    let h = gcd(f, g);
    if h.is_zero() {
        return (h, MPoly::zero(f.nvars), MPoly::zero(f.nvars));
    }
    let p = f.divexact(&h).expect("gcd divides f");
    let q = g.divexact(&h).expect("gcd divides g");
    (h, p, q)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(n: usize, terms: &[(&[u32], i64)]) -> MPoly {
        MPoly::from_terms(
            n,
            terms.iter().map(|(m, c)| (Mono::from_slice(m), Int::from(*c))).collect(),
        )
    }

    #[test]
    fn mul_div() {
        let x = MPoly::var(2, 0);
        let y = MPoly::var(2, 1);
        let a = x.add(&y).add(&MPoly::one(2));
        let b = x.sub(&y);
        let ab = a.mul(&b);
        assert_eq!(ab.divexact(&a).unwrap(), b);
        assert_eq!(ab.divexact(&b).unwrap(), a);
        assert!(ab.divexact(&x).is_none());
    }

    #[test]
    fn gcds() {
        let x = MPoly::var(3, 0);
        let y = MPoly::var(3, 1);
        let z = MPoly::var(3, 2);
        let a = x.add(&y).mul(&x.sub(&z)); // (x+y)(x-z)
        let b = x.add(&y).mul(&y.add(&z)).mul_int(&Int::from(2)); // 2(x+y)(y+z)
        let g = gcd(&a, &b);
        assert_eq!(g, x.add(&y));
        let c = p(1, &[(&[2], 1), (&[0], -1)]); // x^2-1
        let d = p(1, &[(&[2], 1), (&[1], 2), (&[0], 1)]); // x^2+2x+1
        assert_eq!(gcd(&c, &d), p(1, &[(&[1], 1), (&[0], 1)]));
        let e = p(1, &[(&[2], 6), (&[0], -6)]);
        let f = p(1, &[(&[1], 4), (&[0], 4)]);
        assert_eq!(gcd(&e, &f), p(1, &[(&[1], 2), (&[0], 2)]));
        // coprime
        let u = x.mul(&y).add(&MPoly::one(3));
        let v = x.add(&z);
        assert!(gcd(&u, &v).is_one());
        // higher degree
        let w = x.pow(3).sub(&y.pow(3)); // (x-y)(x^2+xy+y^2)
        let t = x.pow(2).sub(&y.pow(2)); // (x-y)(x+y)
        assert_eq!(gcd(&w, &t), x.sub(&y));
    }
}

#[cfg(test)]
mod exponent_storage_tests {
    use super::*;
    #[test]
    fn matrix_exponent_storage_preserves_cancelled_dot_products() {
        let p=MPoly::monomial(1,Mono::from_slice(&[1<<31]),Int::ONE);
        let q=p.neg();
        let result=crate::poly::mpoly::checked_operation(||Ok(MPoly::sum_products(&[(&p,&p),(&p,&q)],1))).unwrap();
        assert!(result.is_zero());
        assert!(crate::poly::mpoly::checked_operation(||Ok(p.mul(&p))).is_err());
        assert!(crate::poly::mpoly::checked_operation(||Ok(p.mul_sub(&p,&p,&p))).unwrap().is_zero());
    }
}

#[cfg(test)]
mod wide_total_degree_tests {
    use super::*;
    #[test]
    fn wide_total_degree_matrix_metadata_is_exact() {
        let p=MPoly::monomial(2,Mono::from_slice(&[1<<31,1<<31]),Int::ONE);
        assert_eq!(p.total_degree(),1u64<<32);
        assert_eq!(MPoly::zero(2).total_degree(),0);
    }
}
