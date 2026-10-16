//! Truncated Laurent series in one variable `t` with integer exponents.
//!
//! Puiseux series in `x` are represented through `x = t^D` for a
//! ramification index `D` chosen by the expansion engine; this module only
//! knows about integer exponents.
//!
//! A [`Ser`] is `Σ c[i] t^(val+i) + O(t^prec)`, or an exact (terminating)
//! sum when `exact` is set. All elementary function kernels are `O(n²)`
//! recurrences derived from first order differential equations (J.C.P.
//! Miller's power recurrence, `E' = g'E` for `exp`, etc.); nothing is
//! computed by repeated differentiation.

use super::coeff::C;
use crate::number::Number;

/// "Infinite" precision of exact series.
pub const INF: i64 = i64::MAX / 4;

/// Series expansion failures.
#[derive(Clone, Debug, PartialEq)]
pub enum SErr {
    /// A finer ramification is needed: restart with `D * k`.
    NeedD(i64),
    /// Expansion point is a pole/essential singularity of a function.
    Pole(String),
    /// Expression not supported by the expansion engine.
    Unsupported(String),
    /// Could not reach the requested precision (e.g. undecidable
    /// cancellation of the leading term).
    Precision,
}

impl std::fmt::Display for SErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SErr::NeedD(k) => write!(f, "needs ramification {k}"),
            SErr::Pole(s) => write!(f, "Cannot expand {s} around the expansion point"),
            SErr::Unsupported(s) => write!(f, "series expansion of {s} is not supported"),
            SErr::Precision => write!(f, "could not determine the leading term"),
        }
    }
}

pub type SRes<T> = Result<T, SErr>;

#[derive(Clone, Debug)]
pub struct Ser {
    /// exponent of `c[0]`
    pub val: i64,
    pub c: Vec<C>,
    /// absolute precision (`O(t^prec)`), [`INF`] when exact
    pub prec: i64,
    pub exact: bool,
}

impl Ser {
    pub fn zero_exact() -> Ser {
        Ser { val: 0, c: vec![], prec: INF, exact: true }
    }

    /// `O(t^prec)`
    pub fn zero_inexact(prec: i64) -> Ser {
        Ser { val: prec, c: vec![], prec, exact: false }
    }

    pub fn constant(c: C) -> Ser {
        Ser::monomial(c, 0)
    }

    pub fn monomial(c: C, k: i64) -> Ser {
        if c.is_zero() {
            return Ser::zero_exact();
        }
        Ser { val: k, c: vec![c], prec: INF, exact: true }
    }

    pub fn from_coeffs(val: i64, c: Vec<C>, prec: i64, exact: bool) -> Ser {
        Ser { val, c, prec: if exact { INF } else { prec }, exact }.norm()
    }

    /// Is this the zero series (exact zero or `O(t^p)` without terms)?
    pub fn is_zero(&self) -> bool {
        self.c.is_empty()
    }

    pub fn is_exact_zero(&self) -> bool {
        self.exact && self.c.is_empty()
    }

    /// Known number of terms after the leading one (relative precision).
    pub fn rel(&self) -> i64 {
        if self.exact {
            INF
        } else {
            self.prec - self.val
        }
    }

    /// Leading coefficient (`None` for zero series).
    pub fn lead(&self) -> Option<&C> {
        self.c.first()
    }

    /// Coefficient of `t^k` (zero outside the stored range).
    pub fn at(&self, k: i64) -> C {
        if k < self.val {
            return C::zero();
        }
        let i = (k - self.val) as usize;
        self.c.get(i).cloned().unwrap_or_else(C::zero)
    }

    /// Strip leading zeros, trim trailing zeros of exact series and
    /// enforce `val + len <= prec`.
    pub fn norm(mut self) -> Ser {
        let lead = self.c.iter().position(|c| !c.is_zero());
        match lead {
            None => {
                self.c.clear();
                if self.exact {
                    self.val = 0;
                    self.prec = INF;
                } else {
                    self.val = self.prec;
                }
                return self;
            }
            Some(k) if k > 0 => {
                self.c.drain(..k);
                self.val += k as i64;
            }
            _ => {}
        }
        if self.exact {
            while self.c.last().is_some_and(|c| c.is_zero()) {
                self.c.pop();
            }
            self.prec = INF;
        } else {
            if self.val >= self.prec {
                self.c.clear();
                self.val = self.prec;
                return self;
            }
            let n = (self.prec - self.val) as usize;
            self.c.truncate(n);
        }
        self
    }

    /// Truncate to absolute precision `n`.
    pub fn trunc(&self, n: i64) -> Ser {
        if self.exact {
            if self.c.is_empty() || self.val + self.c.len() as i64 <= n {
                return self.clone();
            }
        } else if self.prec <= n {
            return self.clone();
        }
        if self.val >= n {
            return Ser::zero_inexact(n);
        }
        let k = ((n - self.val) as usize).min(self.c.len());
        Ser { val: self.val, c: self.c[..k].to_vec(), prec: n, exact: false }.norm()
    }

    pub fn neg(&self) -> Ser {
        Ser { val: self.val, c: self.c.iter().map(C::neg).collect(), prec: self.prec, exact: self.exact }
    }

    pub fn scale(&self, k: &C) -> Ser {
        if k.is_zero() {
            return if self.exact { Ser::zero_exact() } else { Ser::zero_inexact(self.prec) };
        }
        if k.is_one() {
            return self.clone();
        }
        Ser {
            val: self.val,
            c: self.c.iter().map(|c| c.mul(k)).collect(),
            prec: self.prec,
            exact: self.exact,
        }
        .norm()
    }

    /// Multiply by `t^s`.
    pub fn shift(&self, s: i64) -> Ser {
        let mut r = self.clone();
        if r.c.is_empty() && r.exact {
            return r;
        }
        r.val += s;
        if !r.exact {
            r.prec += s;
        }
        r
    }

    pub fn add(&self, o: &Ser) -> Ser {
        if self.is_exact_zero() {
            return o.clone();
        }
        if o.is_exact_zero() {
            return self.clone();
        }
        let exact = self.exact && o.exact;
        let prec = self.prec.min(o.prec);
        let lo = match (self.c.is_empty(), o.c.is_empty()) {
            (true, true) => return Ser::zero_inexact(prec),
            (true, false) => o.val,
            (false, true) => self.val,
            (false, false) => self.val.min(o.val),
        };
        let hi = if exact {
            (self.val + self.c.len() as i64).max(o.val + o.c.len() as i64)
        } else {
            let h1 = if self.exact { self.val + self.c.len() as i64 } else { self.prec };
            let h2 = if o.exact { o.val + o.c.len() as i64 } else { o.prec };
            h1.max(h2).min(prec)
        };
        if hi <= lo {
            return if exact { Ser::zero_exact() } else { Ser::zero_inexact(prec) };
        }
        let mut c = Vec::with_capacity((hi - lo) as usize);
        for k in lo..hi {
            let a = self.get(k);
            let b = o.get(k);
            c.push(match (a, b) {
                (Some(a), Some(b)) => a.add(b),
                (Some(a), None) => a.clone(),
                (None, Some(b)) => b.clone(),
                (None, None) => C::zero(),
            });
        }
        Ser { val: lo, c, prec: if exact { INF } else { prec }, exact }.norm()
    }

    fn get(&self, k: i64) -> Option<&C> {
        if k < self.val {
            return None;
        }
        self.c.get((k - self.val) as usize)
    }

    pub fn sub(&self, o: &Ser) -> Ser {
        self.add(&o.neg())
    }

    pub fn mul(&self, o: &Ser) -> Ser {
        if self.is_exact_zero() || o.is_exact_zero() {
            return Ser::zero_exact();
        }
        let exact = self.exact && o.exact;
        let val = self.val + o.val;
        let prec = if exact {
            INF
        } else {
            (self.val.saturating_add(o.prec)).min(o.val.saturating_add(self.prec))
        };
        if self.c.is_empty() || o.c.is_empty() {
            return Ser::zero_inexact(prec);
        }
        let (la, lb) = (self.c.len(), o.c.len());
        let len = if exact {
            la + lb - 1
        } else {
            let l = prec - val;
            if l <= 0 {
                return Ser::zero_inexact(prec);
            }
            (l as usize).min(la + lb - 1)
        };
        let nz_a: Vec<usize> = (0..la).filter(|&i| !self.c[i].is_zero()).collect();
        let mut c = vec![C::zero(); len];
        // fast path: all numbers
        let all_num = self.c.iter().all(|c| matches!(c, C::N(_))) && o.c.iter().all(|c| matches!(c, C::N(_)));
        if all_num {
            let bn: Vec<&Number> = o.c.iter().map(|c| c.as_number().unwrap()).collect();
            for (k, ck) in c.iter_mut().enumerate() {
                let mut acc = Number::ZERO;
                for &i in &nz_a {
                    if i > k {
                        break;
                    }
                    let j = k - i;
                    if j >= lb {
                        continue;
                    }
                    if bn[j].is_zero() {
                        continue;
                    }
                    acc = super::coeff::nadd(&acc, &super::coeff::nmul(self.c[i].as_number().unwrap(), bn[j]));
                }
                *ck = C::N(acc);
            }
        } else {
            for (k, ck) in c.iter_mut().enumerate() {
                let mut acc = C::zero();
                for &i in &nz_a {
                    if i > k {
                        break;
                    }
                    let j = k - i;
                    if j >= lb || o.c[j].is_zero() {
                        continue;
                    }
                    acc = acc.add(&self.c[i].mul(&o.c[j]));
                }
                *ck = acc;
            }
        }
        Ser { val, c, prec, exact }.norm()
    }

    /// Inverse to absolute precision (at most) `n`. Requires a non-zero
    /// leading coefficient.
    pub fn inv(&self, n: i64) -> SRes<Ser> {
        let a0 = self.lead().ok_or(SErr::Precision)?;
        let inv0 = a0.inv().ok_or(SErr::Precision)?;
        let val = -self.val;
        if self.exact && self.c.len() == 1 {
            return Ok(Ser::monomial(inv0, val));
        }
        let rel = self.rel().min(n - val);
        if rel <= 0 {
            return Ok(Ser::zero_inexact(n.min(val + self.rel())));
        }
        let len = rel as usize;
        let mut b: Vec<C> = Vec::with_capacity(len);
        b.push(inv0.clone());
        let ninv0 = inv0.neg();
        let nz: Vec<usize> = (1..self.c.len()).filter(|&i| !self.c[i].is_zero()).collect();
        for k in 1..len {
            let mut acc = C::zero();
            for &i in &nz {
                if i > k {
                    break;
                }
                acc = acc.add(&self.c[i].mul(&b[k - i]));
            }
            b.push(acc.mul(&ninv0));
        }
        Ok(Ser { val, c: b, prec: val + rel, exact: false }.norm())
    }

    /// d/dt
    pub fn deriv(&self) -> Ser {
        let c: Vec<C> = self
            .c
            .iter()
            .enumerate()
            .map(|(i, c)| c.mul_int(self.val + i as i64))
            .collect();
        if self.c.is_empty() {
            return if self.exact { Ser::zero_exact() } else { Ser::zero_inexact(self.prec - 1) };
        }
        Ser { val: self.val - 1, c, prec: if self.exact { INF } else { self.prec - 1 }, exact: self.exact }.norm()
    }

    /// Antiderivative with zero constant term. Fails on a `1/t` term.
    pub fn integ(&self) -> SRes<Ser> {
        let mut c = Vec::with_capacity(self.c.len());
        for (i, ci) in self.c.iter().enumerate() {
            let k = self.val + i as i64;
            if k == -1 {
                if !ci.is_zero() {
                    return Err(SErr::Unsupported("integration producing a logarithm".into()));
                }
                c.push(C::zero());
            } else {
                c.push(ci.div_int(k + 1));
            }
        }
        if self.c.is_empty() {
            return Ok(if self.exact { Ser::zero_exact() } else { Ser::zero_inexact(self.prec + 1) });
        }
        Ok(Ser { val: self.val + 1, c, prec: if self.exact { INF } else { self.prec + 1 }, exact: self.exact }
            .norm())
    }

    /// Dense coefficients of `t^0 .. t^(len-1)` (requires `val >= 0`).
    pub fn dense0(&self, len: usize) -> Vec<C> {
        (0..len as i64).map(|k| self.at(k)).collect()
    }
}

// ----------------------------------------------------------------------
// kernels on dense coefficient vectors starting at t^0
// ----------------------------------------------------------------------

fn nonzero_from1(g: &[C]) -> Vec<usize> {
    (1..g.len()).filter(|&i| !g[i].is_zero()).collect()
}

/// `exp(g)` for `g[0] == 0` (the constant term of `g` is ignored).
pub fn exp_kernel(g: &[C], len: usize) -> Vec<C> {
    let nz = nonzero_from1(g);
    // ig[i] = i * g[i]
    let ig: Vec<C> = (0..g.len()).map(|i| g[i].mul_int(i as i64)).collect();
    let mut e: Vec<C> = Vec::with_capacity(len);
    e.push(C::one());
    for k in 1..len {
        let mut acc = C::zero();
        for &i in &nz {
            if i > k {
                break;
            }
            acc = acc.add(&ig[i].mul(&e[k - i]));
        }
        e.push(acc.div_int(k as i64));
    }
    e
}

/// `log(a)` for `a[0] == 1`.
pub fn log_kernel(a: &[C], len: usize) -> Vec<C> {
    let nz = nonzero_from1(a);
    let mut l: Vec<C> = Vec::with_capacity(len);
    l.push(C::zero());
    for k in 1..len {
        let mut acc = C::zero();
        for &i in &nz {
            if i >= k {
                break;
            }
            if l[k - i].is_zero() {
                continue;
            }
            acc = acc.add(&a[i].mul(&l[k - i]).mul_int((k - i) as i64));
        }
        let ak = if k < a.len() { a[k].clone() } else { C::zero() };
        l.push(ak.sub(&acc.div_int(k as i64)));
    }
    l
}

/// `a ** r` for `a[0] == 1` (J.C.P. Miller's recurrence).
pub fn pow_kernel(a: &[C], r: &C, len: usize) -> Vec<C> {
    let nz = nonzero_from1(a);
    let mut p: Vec<C> = Vec::with_capacity(len);
    p.push(C::one());
    match r {
        C::N(rn) => {
            let r1 = rn.add(&Number::ONE);
            for k in 1..len {
                let mut acc = C::zero();
                for &i in &nz {
                    if i > k {
                        break;
                    }
                    // (r+1) i - k
                    let f = r1.mul(&Number::int(i as i64)).sub(&Number::int(k as i64));
                    if f.is_zero() || p[k - i].is_zero() {
                        continue;
                    }
                    acc = acc.add(&a[i].mul(&p[k - i]).mul_num(&f));
                }
                p.push(acc.div_int(k as i64));
            }
        }
        C::E(_) => {
            let r1 = r.add(&C::one());
            for k in 1..len {
                let mut acc = C::zero();
                for &i in &nz {
                    if i > k {
                        break;
                    }
                    let f = r1.mul_int(i as i64).sub(&C::int(k as i64));
                    acc = acc.add(&a[i].mul(&p[k - i]).mul(&f));
                }
                p.push(acc.div_int(k as i64));
            }
        }
    }
    p
}

/// `(sin(g), cos(g))` (or `(sinh, cosh)`) for `g[0] == 0`.
pub fn sincos_kernel(g: &[C], len: usize, hyperbolic: bool) -> (Vec<C>, Vec<C>) {
    let nz = nonzero_from1(g);
    let ig: Vec<C> = (0..g.len()).map(|i| g[i].mul_int(i as i64)).collect();
    let mut s: Vec<C> = Vec::with_capacity(len);
    let mut c: Vec<C> = Vec::with_capacity(len);
    s.push(C::zero());
    c.push(C::one());
    for k in 1..len {
        let mut as_ = C::zero();
        let mut ac = C::zero();
        for &i in &nz {
            if i > k {
                break;
            }
            if !c[k - i].is_zero() {
                as_ = as_.add(&ig[i].mul(&c[k - i]));
            }
            if !s[k - i].is_zero() {
                ac = ac.add(&ig[i].mul(&s[k - i]));
            }
        }
        s.push(as_.div_int(k as i64));
        let ck = ac.div_int(k as i64);
        c.push(if hyperbolic { ck } else { ck.neg() });
    }
    (s, c)
}

/// Solution of `T' = alpha (1 + sigma T²) g'` with `T(0) = t0`
/// (`tan`: alpha = sigma = 1, `tanh`: 1, -1, `cot`: -1, 1, `coth`: 1, -1).
pub fn tan_kernel(g: &[C], len: usize, t0: &C, alpha: i64, sigma: i64) -> Vec<C> {
    let nz = nonzero_from1(g);
    let ig: Vec<C> = (0..g.len()).map(|i| g[i].mul_int(i as i64 * alpha)).collect();
    let mut t: Vec<C> = Vec::with_capacity(len);
    let mut q: Vec<C> = Vec::with_capacity(len);
    t.push(t0.clone());
    for k in 1..len {
        // q_{k-1}
        let j = k - 1;
        let mut sq = C::zero();
        for l in 0..=j {
            if t[l].is_zero() || t[j - l].is_zero() {
                continue;
            }
            sq = sq.add(&t[l].mul(&t[j - l]));
        }
        let mut qj = sq.mul_int(sigma);
        if j == 0 {
            qj = qj.add(&C::one());
        }
        q.push(qj);
        let mut acc = C::zero();
        for &i in &nz {
            if i > k {
                break;
            }
            if q[k - i].is_zero() {
                continue;
            }
            acc = acc.add(&ig[i].mul(&q[k - i]));
        }
        t.push(acc.div_int(k as i64));
    }
    t
}
