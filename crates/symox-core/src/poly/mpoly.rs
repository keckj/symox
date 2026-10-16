//! Sparse multivariate polynomials.
//!
//! Terms are stored in two flat arrays (exponent vectors with a fixed
//! stride, and coefficients), sorted in decreasing lexicographic order of
//! the exponent vectors, without zero coefficients. Multiplication and
//! exact division pack exponent vectors into a single `u64` whenever the
//! degrees allow it (variable 0 in the most significant bits, so that
//! integer comparison of packed keys is the lex order).

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fmt;

use rustc_hash::FxHashMap;

use super::int::{Coeff, Field, Int, Rat};

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct MPoly<C> {
    pub nvars: usize,
    exps: Vec<u32>,
    coeffs: Vec<C>,
}

/// Packing of exponent vectors into `u64` keys.
#[derive(Clone, Debug)]
pub(crate) struct Packer {
    shifts: Vec<u32>,
    masks: Vec<u64>,
}

impl Packer {
    /// A packer able to represent exponents up to `maxdeg[v]` per variable.
    pub(crate) fn new(maxdeg: &[u32]) -> Option<Packer> {
        let bits: Vec<u32> = maxdeg.iter().map(|&d| 32 - d.leading_zeros()).collect();
        let total: u32 = bits.iter().sum();
        if total > 64 {
            return None;
        }
        let n = maxdeg.len();
        let mut shifts = vec![0u32; n];
        let mut masks = vec![0u64; n];
        let mut acc = 0u32;
        for v in (0..n).rev() {
            shifts[v] = acc;
            masks[v] = if bits[v] == 0 { 0 } else { (u64::MAX >> (64 - bits[v])) << acc };
            acc += bits[v];
        }
        Some(Packer { shifts, masks })
    }

    #[inline]
    pub(crate) fn pack(&self, e: &[u32]) -> u64 {
        let mut k = 0u64;
        for (v, &x) in e.iter().enumerate() {
            if x != 0 {
                k |= (x as u64) << self.shifts[v];
            }
        }
        k
    }

    #[inline]
    pub(crate) fn unpack(&self, k: u64, out: &mut [u32]) {
        for v in 0..out.len() {
            out[v] = ((k & self.masks[v]) >> self.shifts[v]) as u32;
        }
    }

    /// Does the monomial `a` divide `b` (both packed)?
    #[inline]
    pub(crate) fn divides(&self, a: u64, b: u64) -> bool {
        self.masks.iter().all(|&m| (a & m) <= (b & m))
    }
}

impl<C: Coeff> MPoly<C> {
    pub fn zero(nvars: usize) -> Self {
        MPoly { nvars, exps: Vec::new(), coeffs: Vec::new() }
    }

    pub fn constant(c: C, nvars: usize) -> Self {
        if c.is_zero() {
            return Self::zero(nvars);
        }
        MPoly { nvars, exps: vec![0; nvars], coeffs: vec![c] }
    }

    pub fn one(nvars: usize) -> Self {
        Self::constant(C::one(), nvars)
    }

    /// The variable `x_i`.
    pub fn var(i: usize, nvars: usize) -> Self {
        let mut e = vec![0; nvars];
        e[i] = 1;
        MPoly { nvars, exps: e, coeffs: vec![C::one()] }
    }

    pub fn monomial(exps: &[u32], c: C) -> Self {
        if c.is_zero() {
            return Self::zero(exps.len());
        }
        MPoly { nvars: exps.len(), exps: exps.to_vec(), coeffs: vec![c] }
    }

    /// Build from terms in arbitrary order (duplicates are combined).
    pub fn from_terms(nvars: usize, mut terms: Vec<(Vec<u32>, C)>) -> Self {
        terms.sort_by(|a, b| b.0.cmp(&a.0));
        let mut exps = Vec::with_capacity(terms.len() * nvars);
        let mut coeffs: Vec<C> = Vec::with_capacity(terms.len());
        let mut last: Option<Vec<u32>> = None;
        for (e, c) in terms {
            debug_assert_eq!(e.len(), nvars);
            if last.as_ref() == Some(&e) {
                coeffs.last_mut().unwrap().add_assign(&c);
            } else {
                exps.extend_from_slice(&e);
                coeffs.push(c);
                last = Some(e);
            }
        }
        let mut p = MPoly { nvars, exps, coeffs };
        p.remove_zeros();
        p
    }

    /// Build from terms already sorted (decreasing lex) and distinct, with
    /// non-zero coefficients.
    pub fn from_sorted(nvars: usize, exps: Vec<u32>, coeffs: Vec<C>) -> Self {
        debug_assert_eq!(exps.len(), nvars * coeffs.len());
        MPoly { nvars, exps, coeffs }
    }

    fn remove_zeros(&mut self) {
        if self.coeffs.iter().all(|c| !c.is_zero()) {
            return;
        }
        let n = self.nvars;
        let mut w = 0;
        for r in 0..self.coeffs.len() {
            if !self.coeffs[r].is_zero() {
                if w != r {
                    self.coeffs.swap(w, r);
                    for v in 0..n {
                        self.exps[w * n + v] = self.exps[r * n + v];
                    }
                }
                w += 1;
            }
        }
        self.coeffs.truncate(w);
        self.exps.truncate(w * n);
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    #[inline]
    pub fn exp(&self, i: usize) -> &[u32] {
        &self.exps[i * self.nvars..(i + 1) * self.nvars]
    }

    #[inline]
    pub fn coeff(&self, i: usize) -> &C {
        &self.coeffs[i]
    }

    pub fn coeffs(&self) -> &[C] {
        &self.coeffs
    }

    pub fn coeffs_mut(&mut self) -> &mut [C] {
        &mut self.coeffs
    }

    pub fn iter(&self) -> impl Iterator<Item = (&[u32], &C)> + '_ {
        (0..self.len()).map(move |i| (self.exp(i), &self.coeffs[i]))
    }

    /// Leading coefficient in lex order (zero for the zero polynomial).
    pub fn lc(&self) -> C {
        self.coeffs.first().cloned().unwrap_or_else(C::zero)
    }

    pub fn lm(&self) -> &[u32] {
        self.exp(0)
    }

    pub fn is_constant(&self) -> bool {
        self.is_zero() || (self.len() == 1 && self.exp(0).iter().all(|&e| e == 0))
    }

    pub fn constant_value(&self) -> Option<C> {
        if self.is_zero() {
            Some(C::zero())
        } else if self.is_constant() {
            Some(self.coeffs[0].clone())
        } else {
            None
        }
    }

    /// The constant term.
    pub fn constant_term(&self) -> C {
        match self.len() {
            0 => C::zero(),
            l if self.exp(l - 1).iter().all(|&e| e == 0) => self.coeffs[l - 1].clone(),
            _ => C::zero(),
        }
    }

    pub fn is_one(&self) -> bool {
        self.is_constant() && !self.is_zero() && self.coeffs[0].is_one()
    }

    pub fn is_monomial(&self) -> bool {
        self.len() == 1
    }

    /// Degree in variable `v`, `-1` for the zero polynomial.
    pub fn degree(&self, v: usize) -> i64 {
        if self.is_zero() {
            return -1;
        }
        let n = self.nvars;
        if v == 0 {
            return self.exps[0] as i64;
        }
        (0..self.len()).map(|i| self.exps[i * n + v]).max().unwrap_or(0) as i64
    }

    /// Maximal exponent of each variable.
    pub fn degrees(&self) -> Vec<u32> {
        let n = self.nvars;
        let mut d = vec![0u32; n];
        for i in 0..self.len() {
            for v in 0..n {
                d[v] = d[v].max(self.exps[i * n + v]);
            }
        }
        d
    }

    pub fn total_degree(&self) -> i64 {
        if self.is_zero() {
            return -1;
        }
        self.iter().map(|(e, _)| e.iter().map(|&x| x as i64).sum::<i64>()).max().unwrap()
    }

    /// Variables actually occurring.
    pub fn used_vars(&self) -> Vec<usize> {
        let d = self.degrees();
        (0..self.nvars).filter(|&v| d[v] > 0).collect()
    }

    pub fn neg(&self) -> Self {
        MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self.coeffs.iter().map(|c| c.neg()).collect(),
        }
    }

    pub fn scale(&self, k: &C) -> Self {
        if k.is_zero() {
            return Self::zero(self.nvars);
        }
        if k.is_one() {
            return self.clone();
        }
        let mut p = MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self.coeffs.iter().map(|c| c.mul(k)).collect(),
        };
        p.remove_zeros();
        p
    }

    /// Divide every coefficient exactly by `k` (`None` if not exact).
    pub fn div_coeff(&self, k: &C) -> Option<Self> {
        if k.is_one() {
            return Some(self.clone());
        }
        let mut coeffs = Vec::with_capacity(self.len());
        for c in &self.coeffs {
            coeffs.push(c.div_exact(k)?);
        }
        Some(MPoly { nvars: self.nvars, exps: self.exps.clone(), coeffs })
    }

    fn merge(&self, o: &Self, negate: bool) -> Self {
        let n = self.nvars;
        debug_assert_eq!(n, o.nvars);
        let (la, lb) = (self.len(), o.len());
        let mut exps = Vec::with_capacity((la + lb) * n);
        let mut coeffs = Vec::with_capacity(la + lb);
        let (mut i, mut j) = (0, 0);
        while i < la && j < lb {
            match self.exp(i).cmp(o.exp(j)) {
                Ordering::Greater => {
                    exps.extend_from_slice(self.exp(i));
                    coeffs.push(self.coeffs[i].clone());
                    i += 1;
                }
                Ordering::Less => {
                    exps.extend_from_slice(o.exp(j));
                    coeffs.push(if negate { o.coeffs[j].neg() } else { o.coeffs[j].clone() });
                    j += 1;
                }
                Ordering::Equal => {
                    let c = if negate {
                        self.coeffs[i].sub(&o.coeffs[j])
                    } else {
                        self.coeffs[i].add(&o.coeffs[j])
                    };
                    if !c.is_zero() {
                        exps.extend_from_slice(self.exp(i));
                        coeffs.push(c);
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        while i < la {
            exps.extend_from_slice(self.exp(i));
            coeffs.push(self.coeffs[i].clone());
            i += 1;
        }
        while j < lb {
            exps.extend_from_slice(o.exp(j));
            coeffs.push(if negate { o.coeffs[j].neg() } else { o.coeffs[j].clone() });
            j += 1;
        }
        MPoly { nvars: n, exps, coeffs }
    }

    pub fn add(&self, o: &Self) -> Self {
        if self.is_zero() {
            return o.clone();
        }
        if o.is_zero() {
            return self.clone();
        }
        self.merge(o, false)
    }

    pub fn sub(&self, o: &Self) -> Self {
        if o.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return o.neg();
        }
        self.merge(o, true)
    }

    /// Multiply by the monomial `c * x^e`.
    pub fn mul_term(&self, e: &[u32], c: &C) -> Self {
        if c.is_zero() || self.is_zero() {
            return Self::zero(self.nvars);
        }
        let n = self.nvars;
        let mut exps = self.exps.clone();
        for i in 0..self.len() {
            for v in 0..n {
                exps[i * n + v] += e[v];
            }
        }
        let mut p = MPoly {
            nvars: n,
            exps,
            coeffs: self.coeffs.iter().map(|x| x.mul(c)).collect(),
        };
        p.remove_zeros();
        p
    }

    pub fn mul(&self, o: &Self) -> Self {
        let n = self.nvars;
        if self.is_zero() || o.is_zero() {
            return Self::zero(n);
        }
        if o.len() == 1 {
            return self.mul_term(o.exp(0), &o.coeffs[0]);
        }
        if self.len() == 1 {
            return o.mul_term(self.exp(0), &self.coeffs[0]);
        }
        let da = self.degrees();
        let db = o.degrees();
        let dsum: Vec<u32> = da.iter().zip(&db).map(|(a, b)| a + b).collect();
        if let Some(pk) = Packer::new(&dsum) {
            let ka: Vec<u64> = (0..self.len()).map(|i| pk.pack(self.exp(i))).collect();
            let kb: Vec<u64> = (0..o.len()).map(|i| pk.pack(o.exp(i))).collect();
            let maxkey = pk.pack(&dsum);
            let pairs = self.len().saturating_mul(o.len());
            let mut out: Vec<(u64, C)>;
            if maxkey < (1 << 22) && (maxkey as usize) <= pairs.saturating_mul(8) {
                // dense accumulation
                let mut acc: Vec<C> = vec![C::zero(); maxkey as usize + 1];
                for (i, ca) in self.coeffs.iter().enumerate() {
                    let k1 = ka[i];
                    for (j, cb) in o.coeffs.iter().enumerate() {
                        acc[(k1 + kb[j]) as usize].add_mul(ca, cb);
                    }
                }
                out = Vec::new();
                for (k, c) in acc.into_iter().enumerate().rev() {
                    if !c.is_zero() {
                        out.push((k as u64, c));
                    }
                }
            } else {
                let mut map: FxHashMap<u64, C> =
                    FxHashMap::with_capacity_and_hasher(pairs.min(1 << 20), Default::default());
                for (i, ca) in self.coeffs.iter().enumerate() {
                    let k1 = ka[i];
                    for (j, cb) in o.coeffs.iter().enumerate() {
                        map.entry(k1 + kb[j]).or_insert_with(C::zero).add_mul(ca, cb);
                    }
                }
                out = map.into_iter().filter(|(_, c)| !c.is_zero()).collect();
                out.sort_unstable_by(|a, b| b.0.cmp(&a.0));
            }
            return Self::from_packed(&pk, n, out);
        }
        // generic path
        let mut map: FxHashMap<Vec<u32>, C> = FxHashMap::default();
        let mut buf = vec![0u32; n];
        for i in 0..self.len() {
            for j in 0..o.len() {
                for v in 0..n {
                    buf[v] = self.exp(i)[v] + o.exp(j)[v];
                }
                map.entry(buf.clone()).or_insert_with(C::zero).add_mul(&self.coeffs[i], &o.coeffs[j]);
            }
        }
        Self::from_terms(n, map.into_iter().collect())
    }

    pub(crate) fn from_packed(pk: &Packer, n: usize, terms: Vec<(u64, C)>) -> Self {
        let mut exps = vec![0u32; terms.len() * n];
        let mut coeffs = Vec::with_capacity(terms.len());
        for (i, (k, c)) in terms.into_iter().enumerate() {
            pk.unpack(k, &mut exps[i * n..(i + 1) * n]);
            coeffs.push(c);
        }
        MPoly { nvars: n, exps, coeffs }
    }

    pub fn square(&self) -> Self {
        self.mul(self)
    }

    pub fn pow(&self, k: u32) -> Self {
        if k == 0 {
            return Self::one(self.nvars);
        }
        if self.len() == 1 {
            let e: Vec<u32> = self.exp(0).iter().map(|&x| x * k).collect();
            let mut c = C::one();
            for _ in 0..k {
                c = c.mul(&self.coeffs[0]);
            }
            return Self::monomial(&e, c);
        }
        // repeated multiplication by the base is faster for sparse bases
        let mut r = self.clone();
        for _ in 1..k {
            r = r.mul(self);
        }
        r
    }

    /// Exact division by `g` (`None` if `g` does not divide `self`).
    pub fn div_exact(&self, g: &Self) -> Option<Self> {
        let n = self.nvars;
        if g.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some(self.clone());
        }
        if g.len() == 1 {
            // monomial divisor
            let ge = g.exp(0);
            let gc = &g.coeffs[0];
            let mut exps = Vec::with_capacity(self.exps.len());
            let mut coeffs = Vec::with_capacity(self.len());
            for (e, c) in self.iter() {
                for v in 0..n {
                    if e[v] < ge[v] {
                        return None;
                    }
                    exps.push(e[v] - ge[v]);
                }
                coeffs.push(c.div_exact(gc)?);
            }
            return Some(MPoly { nvars: n, exps, coeffs });
        }
        let df = self.degrees();
        let dg = g.degrees();
        if df.iter().zip(&dg).any(|(a, b)| a < b) {
            return None;
        }
        let pk = Packer::new(&df)?;
        Some(self.div_heap(g, &pk, &df, &dg)?)
    }

    /// Johnson's heap division, exact case.
    fn div_heap(&self, g: &Self, pk: &Packer, df: &[u32], dg: &[u32]) -> Option<Self> {
        let n = self.nvars;
        let fk: Vec<u64> = (0..self.len()).map(|i| pk.pack(self.exp(i))).collect();
        let gk: Vec<u64> = (0..g.len()).map(|i| pk.pack(g.exp(i))).collect();
        let glc = &g.coeffs[0];
        let g0 = gk[0];
        let mut qk: Vec<u64> = Vec::new();
        let mut qc: Vec<C> = Vec::new();
        // heap of (key, i (quotient index), j (divisor index))
        let mut heap: BinaryHeap<(u64, u32, u32)> = BinaryHeap::new();
        let mut fi = 0usize;
        let mut buf = vec![0u32; n];
        // bound on quotient exponents
        let qmax: Vec<u32> = df.iter().zip(dg).map(|(a, b)| a - b).collect();
        loop {
            let hk = heap.peek().map(|t| t.0);
            let fkey = if fi < fk.len() { Some(fk[fi]) } else { None };
            let m = match (fkey, hk) {
                (None, None) => break,
                (Some(a), None) => a,
                (None, Some(b)) => b,
                (Some(a), Some(b)) => a.max(b),
            };
            let mut c = C::zero();
            if fkey == Some(m) {
                c = self.coeffs[fi].clone();
                fi += 1;
            }
            while let Some(&(k, i, j)) = heap.peek() {
                if k != m {
                    break;
                }
                heap.pop();
                let prod = qc[i as usize].mul(&g.coeffs[j as usize]);
                c = c.sub(&prod);
                let j2 = j as usize + 1;
                if j2 < gk.len() {
                    heap.push((qk[i as usize] + gk[j2], i, j2 as u32));
                }
            }
            if c.is_zero() {
                continue;
            }
            // new quotient term m / lm(g)
            if !pk.divides(g0, m) {
                return None;
            }
            let qkey = m - g0;
            pk.unpack(qkey, &mut buf);
            if buf.iter().zip(&qmax).any(|(a, b)| a > b) {
                return None;
            }
            let cq = c.div_exact(glc)?;
            let idx = qk.len() as u32;
            qk.push(qkey);
            qc.push(cq);
            if gk.len() > 1 {
                heap.push((qkey + gk[1], idx, 1));
            }
        }
        let terms: Vec<(u64, C)> = qk.into_iter().zip(qc).collect();
        Some(Self::from_packed(pk, n, terms))
    }

    /// Map coefficients into another ring.
    pub fn map_coeffs<D: Coeff>(&self, f: impl Fn(&C) -> D) -> MPoly<D> {
        let mut p = MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self.coeffs.iter().map(f).collect(),
        };
        p.remove_zeros();
        p
    }

    /// Partial derivative with respect to `v`.
    pub fn derivative(&self, v: usize) -> Self {
        let n = self.nvars;
        let mut exps = Vec::with_capacity(self.exps.len());
        let mut coeffs = Vec::with_capacity(self.len());
        for (e, c) in self.iter() {
            if e[v] == 0 {
                continue;
            }
            let k = C::from_i64(e[v] as i64);
            let nc = c.mul(&k);
            if nc.is_zero() {
                continue;
            }
            let start = exps.len();
            exps.extend_from_slice(e);
            exps[start + v] -= 1;
            coeffs.push(nc);
        }
        // order is preserved: lowering the same variable by one keeps lex
        // order among terms with e[v] > 0
        let _ = n;
        MPoly { nvars: self.nvars, exps, coeffs }
    }

    /// Substitute `x_v = a` (the variable stays in the ring, with degree 0).
    pub fn eval(&self, v: usize, a: &C) -> Self {
        let n = self.nvars;
        if !a.is_zero() && (v + 1..n).all(|w| self.iter().all(|(e, _)| e[w] == 0)) {
            return self.eval_last(v, a);
        }
        if a.is_zero() {
            let mut exps = Vec::new();
            let mut coeffs = Vec::new();
            for (e, c) in self.iter() {
                if e[v] == 0 {
                    exps.extend_from_slice(e);
                    coeffs.push(c.clone());
                }
            }
            return MPoly { nvars: n, exps, coeffs };
        }
        let mut powers: Vec<C> = vec![C::one()];
        let mut map: FxHashMap<Vec<u32>, C> = FxHashMap::default();
        for (e, c) in self.iter() {
            let k = e[v] as usize;
            while powers.len() <= k {
                let p = powers.last().unwrap().mul(a);
                powers.push(p);
            }
            let mut ne = e.to_vec();
            ne[v] = 0;
            map.entry(ne).or_insert_with(C::zero).add_mul(c, &powers[k]);
        }
        Self::from_terms(n, map.into_iter().collect())
    }

    /// `eval` when no variable after `v` occurs: terms differing only in
    /// `x_v` are contiguous, so the evaluation is a single pass.
    fn eval_last(&self, v: usize, a: &C) -> Self {
        let n = self.nvars;
        let mut powers: Vec<C> = vec![C::one()];
        let mut exps: Vec<u32> = Vec::with_capacity(self.exps.len());
        let mut coeffs: Vec<C> = Vec::with_capacity(self.len());
        let mut i = 0;
        while i < self.len() {
            let e = self.exp(i);
            let mut acc = C::zero();
            let mut j = i;
            while j < self.len() && self.exp(j)[..v] == e[..v] {
                let k = self.exp(j)[v] as usize;
                while powers.len() <= k {
                    let p = powers.last().unwrap().mul(a);
                    powers.push(p);
                }
                acc.add_mul(&self.coeffs[j], &powers[k]);
                j += 1;
            }
            if !acc.is_zero() {
                exps.extend_from_slice(&e[..v]);
                exps.extend(std::iter::repeat_n(0, n - v));
                coeffs.push(acc);
            }
            i = j;
        }
        MPoly { nvars: n, exps, coeffs }
    }

    /// Evaluate all variables.
    pub fn eval_all(&self, vals: &[C]) -> C {
        let mut s = C::zero();
        for (e, c) in self.iter() {
            let mut t = c.clone();
            for (v, &k) in e.iter().enumerate() {
                for _ in 0..k {
                    t = t.mul(&vals[v]);
                }
            }
            s.add_assign(&t);
        }
        s
    }

    /// Split into coefficients with respect to `v`: `[(k, c_k)]` with
    /// `self = Σ c_k x_v^k`, decreasing `k`, `c_k` free of `x_v`.
    pub fn coeffs_in(&self, v: usize) -> Vec<(u32, Self)> {
        let n = self.nvars;
        if v == 0 {
            // already grouped
            let mut out: Vec<(u32, Self)> = Vec::new();
            for (e, c) in self.iter() {
                let k = e[0];
                if out.last().map(|x| x.0) != Some(k) {
                    out.push((k, Self::zero(n)));
                }
                let p = &mut out.last_mut().unwrap().1;
                p.exps.push(0);
                p.exps.extend_from_slice(&e[1..]);
                p.coeffs.push(c.clone());
            }
            return out;
        }
        let mut groups: FxHashMap<u32, Vec<(Vec<u32>, C)>> = FxHashMap::default();
        for (e, c) in self.iter() {
            let mut ne = e.to_vec();
            ne[v] = 0;
            groups.entry(e[v]).or_default().push((ne, c.clone()));
        }
        let mut out: Vec<(u32, Self)> = groups
            .into_iter()
            .map(|(k, ts)| {
                // terms are already in lex order once x_v is removed
                let mut exps = Vec::with_capacity(ts.len() * n);
                let mut coeffs = Vec::with_capacity(ts.len());
                for (e, c) in ts {
                    exps.extend_from_slice(&e);
                    coeffs.push(c);
                }
                (k, MPoly { nvars: n, exps, coeffs })
            })
            .collect();
        out.sort_by(|a, b| b.0.cmp(&a.0));
        out
    }

    /// Inverse of [`coeffs_in`].
    pub fn from_coeffs_in(v: usize, n: usize, cs: &[(u32, Self)]) -> Self {
        let mut acc = Self::zero(n);
        for (k, c) in cs {
            if c.is_zero() {
                continue;
            }
            let mut e = vec![0u32; n];
            e[v] = *k;
            acc = acc.add(&c.mul_term(&e, &C::one()));
        }
        acc
    }

    /// Leading coefficient with respect to `x_v`.
    pub fn lc_in(&self, v: usize) -> Self {
        let n = self.nvars;
        if self.is_zero() {
            return Self::zero(n);
        }
        let d = self.degree(v) as u32;
        let mut exps = Vec::new();
        let mut coeffs = Vec::new();
        for (e, c) in self.iter() {
            if e[v] == d {
                let s = exps.len();
                exps.extend_from_slice(e);
                exps[s + v] = 0;
                coeffs.push(c.clone());
            }
        }
        MPoly { nvars: n, exps, coeffs }
    }

    /// Coefficient of `x_v^k` (free of `x_v`).
    pub fn coeff_of(&self, v: usize, k: u32) -> Self {
        let n = self.nvars;
        let mut exps = Vec::new();
        let mut coeffs = Vec::new();
        for (e, c) in self.iter() {
            if e[v] == k {
                let s = exps.len();
                exps.extend_from_slice(e);
                exps[s + v] = 0;
                coeffs.push(c.clone());
            }
        }
        MPoly { nvars: n, exps, coeffs }
    }

    /// Minimal exponent of each variable (the monomial content).
    pub fn min_exps(&self) -> Vec<u32> {
        let n = self.nvars;
        if self.is_zero() {
            return vec![0; n];
        }
        let mut m = self.exp(0).to_vec();
        for i in 1..self.len() {
            for v in 0..n {
                m[v] = m[v].min(self.exps[i * n + v]);
            }
        }
        m
    }

    /// Divide by the monomial `x^e` (caller guarantees divisibility).
    pub fn div_monomial(&self, e: &[u32]) -> Self {
        let n = self.nvars;
        let mut p = self.clone();
        for i in 0..p.len() {
            for v in 0..n {
                p.exps[i * n + v] -= e[v];
            }
        }
        p
    }

    /// Re-embed into a ring with `new_n` variables, variable `v` going to
    /// `perm[v]`.
    pub fn permute(&self, perm: &[usize], new_n: usize) -> Self {
        let mut terms = Vec::with_capacity(self.len());
        for (e, c) in self.iter() {
            let mut ne = vec![0u32; new_n];
            for (v, &k) in e.iter().enumerate() {
                if k != 0 {
                    ne[perm[v]] += k;
                }
            }
            terms.push((ne, c.clone()));
        }
        Self::from_terms(new_n, terms)
    }

    /// Substitute `x_v -> x_v + a` (Taylor shift in one variable).
    pub fn shift(&self, v: usize, a: &C) -> Self {
        if a.is_zero() {
            return self.clone();
        }
        let n = self.nvars;
        let cs = self.coeffs_in(v);
        // Horner in (x_v + a)
        let mut lin = Self::var(v, n);
        lin = lin.add(&Self::constant(a.clone(), n));
        let mut acc = Self::zero(n);
        let mut idx = 0;
        let top = cs.first().map(|x| x.0).unwrap_or(0);
        for k in (0..=top).rev() {
            acc = acc.mul(&lin);
            if idx < cs.len() && cs[idx].0 == k {
                acc = acc.add(&cs[idx].1);
                idx += 1;
            }
        }
        acc
    }

    /// Coefficient of the monomial `e` (zero if absent).
    pub fn coeff_at(&self, e: &[u32]) -> C {
        // binary search in decreasing order
        let (mut lo, mut hi) = (0usize, self.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            match self.exp(mid).cmp(e) {
                Ordering::Equal => return self.coeffs[mid].clone(),
                Ordering::Greater => lo = mid + 1,
                Ordering::Less => hi = mid,
            }
        }
        C::zero()
    }

    /// Terms of total degree `<= d` in the variables `vars`.
    pub fn truncate_degree(&self, vars: &[usize], d: u32) -> Self {
        let mut exps = Vec::new();
        let mut coeffs = Vec::new();
        for (e, c) in self.iter() {
            let s: u32 = vars.iter().map(|&v| e[v]).sum();
            if s <= d {
                exps.extend_from_slice(e);
                coeffs.push(c.clone());
            }
        }
        MPoly { nvars: self.nvars, exps, coeffs }
    }
}

impl<C: Field> MPoly<C> {
    /// Make the leading coefficient one.
    pub fn monic(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }
        let l = self.coeffs[0].inv();
        self.scale(&l)
    }
}

impl MPoly<Int> {
    /// Gcd of the coefficients (non-negative).
    pub fn content(&self) -> Int {
        let mut g = Int::zero();
        for c in &self.coeffs {
            g = g.gcd(c);
            if g.is_one() {
                break;
            }
        }
        g
    }

    /// `(content, primitive part)`, the content carrying the sign of the
    /// leading coefficient (so the primitive part has a positive one).
    pub fn primitive(&self) -> (Int, Self) {
        if self.is_zero() {
            return (Int::zero(), self.clone());
        }
        let mut g = self.content();
        if self.coeffs[0].is_negative() {
            g = g.neg();
        }
        if g.is_one() {
            return (g, self.clone());
        }
        let p = MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self.coeffs.iter().map(|c| c.div_exact_unchecked(&g)).collect(),
        };
        (g, p)
    }

    pub fn max_norm(&self) -> Int {
        self.coeffs.iter().map(|c| c.abs()).max().unwrap_or_else(Int::zero)
    }

    pub fn to_rat(&self) -> MPoly<Rat> {
        MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self.coeffs.iter().map(|c| Rat::from(c.clone())).collect(),
        }
    }
}

impl MPoly<Rat> {
    /// `(d, d*self)` with `d` the lcm of the denominators.
    pub fn clear_denoms(&self) -> (Int, MPoly<Int>) {
        let mut d = Int::one();
        for c in &self.coeffs {
            if !c.d.is_one() {
                d = d.lcm(&c.d);
            }
        }
        let p = MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self
                .coeffs
                .iter()
                .map(|c| if d.is_one() { c.n.clone() } else { c.n.mul(&d.div_exact_unchecked(&c.d)) })
                .collect(),
        };
        (d, p)
    }

    pub fn is_integral(&self) -> bool {
        self.coeffs.iter().all(|c| c.d.is_one())
    }

    /// Integer polynomial if all coefficients are integers.
    pub fn to_int(&self) -> Option<MPoly<Int>> {
        if !self.is_integral() {
            return None;
        }
        Some(MPoly {
            nvars: self.nvars,
            exps: self.exps.clone(),
            coeffs: self.coeffs.iter().map(|c| c.n.clone()).collect(),
        })
    }
}

impl<C: Coeff> fmt::Debug for MPoly<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        for (i, (e, c)) in self.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{c:?}")?;
            for (v, &k) in e.iter().enumerate() {
                if k == 1 {
                    write!(f, "*x{v}")?;
                } else if k > 1 {
                    write!(f, "*x{v}^{k}")?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(n: usize, t: &[(&[u32], i64)]) -> MPoly<Int> {
        MPoly::from_terms(n, t.iter().map(|(e, c)| (e.to_vec(), Int::S(*c))).collect())
    }

    #[test]
    fn arith() {
        // (x + y)^2
        let a = p(2, &[(&[1, 0], 1), (&[0, 1], 1)]);
        let s = a.mul(&a);
        assert_eq!(s, p(2, &[(&[2, 0], 1), (&[1, 1], 2), (&[0, 2], 1)]));
        assert_eq!(s.div_exact(&a), Some(a.clone()));
        let b = p(2, &[(&[1, 0], 1), (&[0, 0], 1)]);
        assert_eq!(s.div_exact(&b), None);
        assert_eq!(s.sub(&s), MPoly::zero(2));
        let c = a.pow(5);
        assert_eq!(c.len(), 6);
        assert_eq!(c.div_exact(&a.pow(3)), Some(a.pow(2)));
        assert_eq!(s.derivative(0), p(2, &[(&[1, 0], 2), (&[0, 1], 2)]));
        assert_eq!(s.eval(1, &Int::S(1)), p(2, &[(&[2, 0], 1), (&[1, 0], 2), (&[0, 0], 1)]));
        let cs = s.coeffs_in(1);
        assert_eq!(cs.len(), 3);
        assert_eq!(MPoly::from_coeffs_in(1, 2, &cs), s);
        assert_eq!(b.shift(0, &Int::S(-1)), p(2, &[(&[1, 0], 1)]));
    }
}
