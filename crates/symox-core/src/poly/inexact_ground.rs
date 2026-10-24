//! Polynomials over the inexact composite domains `RR[y...]`, `CC[y...]`,
//! `RR(y...)` and `CC(y...)` at any working precision.
//!
//! Like the composite Gaussian grounds, these polynomials are stored as `EX`
//! polynomials whose coefficients are kept in the normal form of the
//! reference ground element (`Domain.to_sympy`). Every algorithm converts the
//! coefficients into ground elements that reproduce the reference
//! `PolyElement`/`FracElement` objects over `RealField(prec)` and
//! `ComplexField(prec)`: an insertion-ordered sparse dictionary of `mpf`/`mpc`
//! scalars rounded at the domain precision with mpmath's nearest rounding.
//! The ring operations follow `sympy.polys.rings` and `sympy.polys.fields`
//! term by term, ground gcds go through the exact domain obtained by the
//! reference `to_rational` (with its bounded denominator), and the dense
//! algorithms (`dup_*`) are the reference ones over the composite domain, so
//! rounding happens at the same places as in the reference.
use super::algebraic::{AlgebraicScope, AN};
use super::convert::{PResult, PolyError};
use super::expression_domain::ExpressionPoly;
use super::gaussian as gi;
use super::gaussian_ground::GaussianGround;
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use crate::bigfloat::{dps_to_prec, prec_to_dps, BigFloat as BF};
use crate::expr::{Expr, Kind};
use crate::number::Number;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use rustc_hash::FxHashMap;
use std::cell::RefCell;

/// The inexact ground of a composite domain: `RealField(precision)` or
/// `ComplexField(precision)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Inexact {
    pub precision: u32,
    pub complex: bool,
}

fn reference(kind: &'static str, message: impl Into<String>) -> PolyError {
    PolyError::Reference { kind, message: message.into() }
}
fn zero_division() -> PolyError {
    reference("ZeroDivisionError", "polynomial division")
}

// ----------------------------------------------------------------------
// scalars: mpf / mpc at the domain precision
// ----------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct C {
    re: BF,
    im: BF,
}

/// The exact product of two binary floats.
fn xmul(a: &BF, b: &BF) -> BF {
    BF::round(&a.man * &b.man, &a.exp + &b.exp, u32::MAX)
}

fn bitlen_u64(n: u64) -> u64 {
    64 - n.leading_zeros() as u64
}

/// mpmath `mpf_pow_int` for a nonnegative exponent (nearest rounding).
fn pow_real(s: &BF, n: u64, prec: u32) -> BF {
    if n == 0 {
        return BF::one();
    }
    if n == 1 {
        return s.round_to(prec);
    }
    if s.is_zero() {
        return BF::zero();
    }
    if n == 2 {
        return BF::round(&s.man * &s.man, &s.exp + &s.exp, prec);
    }
    let neg = s.is_negative() && n % 2 == 1;
    let sign = |m: BigInt| if neg { -m } else { m };
    let mut man = s.man.abs();
    let mut exp = s.exp.clone();
    if man.is_one() {
        return BF::round(sign(BigInt::one()), &exp * n, u32::MAX);
    }
    let bc = man.bits();
    if bc.saturating_mul(n) < 1000 {
        let m = num_traits::pow(man, n as usize);
        return BF::round(sign(m), &exp * n, prec);
    }
    let workprec = prec as u64 + 4 * bitlen_u64(n) + 4;
    let (mut pm, mut pe) = (BigInt::one(), crate::bigfloat::Exponent::ZERO);
    let mut n = n;
    loop {
        if n & 1 == 1 {
            pm = &pm * &man;
            pe += &exp;
            let pbc = pm.bits();
            if pbc > workprec {
                let sh = pbc - workprec;
                pm >>= sh;
                pe += sh;
            }
            n -= 1;
            if n == 0 {
                break;
            }
        }
        man = &man * &man;
        exp = &exp + &exp;
        let bc = man.bits();
        if bc > workprec {
            let sh = bc - workprec;
            man >>= sh;
            exp += sh;
        }
        n /= 2;
    }
    BF::round(sign(pm), pe, prec)
}

impl Inexact {
    fn p(&self) -> u32 {
        self.precision
    }
    /// The precision of the `Float` objects returned by `to_sympy`.
    fn float_prec(&self) -> u32 {
        dps_to_prec(prec_to_dps(self.precision))
    }
    pub fn display(&self) -> &'static str {
        if self.complex {
            "CC"
        } else {
            "RR"
        }
    }
    fn zero(&self) -> C {
        C { re: BF::zero(), im: BF::zero() }
    }
    fn one(&self) -> C {
        C { re: BF::one(), im: BF::zero() }
    }
    fn is_zero(&self, a: &C) -> bool {
        a.re.is_zero() && a.im.is_zero()
    }
    fn add(&self, a: &C, b: &C) -> C {
        C { re: a.re.add(&b.re, self.p()), im: a.im.add(&b.im, self.p()) }
    }
    fn sub(&self, a: &C, b: &C) -> C {
        C { re: a.re.sub(&b.re, self.p()), im: a.im.sub(&b.im, self.p()) }
    }
    fn neg(&self, a: &C) -> C {
        C { re: a.re.neg(), im: a.im.neg() }
    }
    fn mul(&self, a: &C, b: &C) -> C {
        if !self.complex {
            return C { re: a.re.mul(&b.re, self.p()), im: BF::zero() };
        }
        let (p, q, r, s) = (xmul(&a.re, &b.re), xmul(&a.im, &b.im), xmul(&a.re, &b.im), xmul(&a.im, &b.re));
        C { re: p.sub(&q, self.p()), im: r.add(&s, self.p()) }
    }
    fn div(&self, a: &C, b: &C) -> PResult<C> {
        let zd = || reference("ZeroDivisionError", "");
        if !self.complex {
            return Ok(C { re: a.re.div(&b.re, self.p()).ok_or_else(zd)?, im: BF::zero() });
        }
        let wp = self.p() + 10;
        let mag = xmul(&b.re, &b.re).add(&xmul(&b.im, &b.im), wp);
        let t = xmul(&a.re, &b.re).add(&xmul(&a.im, &b.im), wp);
        let u = xmul(&a.im, &b.re).sub(&xmul(&a.re, &b.im), wp);
        Ok(C { re: t.div(&mag, self.p()).ok_or_else(zd)?, im: u.div(&mag, self.p()).ok_or_else(zd)? })
    }
    /// `mpf_mul_int`/`mpc_mul_int`.
    fn mul_int(&self, a: &C, k: &BigInt) -> C {
        let m = |x: &BF| BF::round(&x.man * k, x.exp.clone(), self.p());
        C { re: m(&a.re), im: m(&a.im) }
    }
    fn from_int(&self, n: &BigInt) -> C {
        C { re: BF::round(n.clone(), 0, self.p()), im: BF::zero() }
    }
    fn from_rat(&self, r: &Rat) -> C {
        self.from_parts(r, &Rat::zero())
    }
    /// `RR.from_QQ` and `CC.from_QQ_I`: `dtype(numerator)/denominator`
    /// for each part.
    fn from_parts(&self, x: &Rat, y: &Rat) -> C {
        let part = |r: &Rat| -> BF {
            if r.is_zero() {
                return BF::zero();
            }
            let n = BF::round(r.n.to_big(), 0, self.p());
            n.div(&BF::from_int(&r.d.to_big()), self.p()).expect("nonzero denominator")
        };
        let re = part(x);
        let im = part(y);
        if self.complex && !y.is_zero() {
            // `dtype(x)/q + dtype(0, y)/s`: an exact sum of disjoint parts
            return self.add(&C { re, im: BF::zero() }, &C { re: BF::zero(), im });
        }
        C { re, im }
    }
    /// `mpc_pow_int`/`mpf_pow_int` for a nonnegative exponent.
    fn pow(&self, a: &C, n: u64) -> PResult<C> {
        let p = self.p();
        if a.im.is_zero() {
            return Ok(C { re: pow_real(&a.re, n, p), im: BF::zero() });
        }
        if a.re.is_zero() {
            let v = pow_real(&a.im, n, p);
            return Ok(match n % 4 {
                0 => C { re: v, im: BF::zero() },
                1 => C { re: BF::zero(), im: v },
                2 => C { re: v.neg(), im: BF::zero() },
                _ => C { re: BF::zero(), im: v.neg() },
            });
        }
        if n == 0 {
            return Ok(self.one());
        }
        if n == 1 {
            return Ok(C { re: a.re.round_to(p), im: a.im.round_to(p) });
        }
        if n == 2 {
            let re = xmul(&a.re, &a.re).sub(&xmul(&a.im, &a.im), p);
            let im = a.re.mul(&a.im, p).ldexp(1);
            return Ok(C { re, im });
        }
        let (mut am, mut ae) = (a.re.man.clone(), a.re.exp.clone());
        let (mut bm, mut be) = (a.im.man.clone(), a.im.exp.clone());
        let de = &ae - &be;
        let size = (n as u128).saturating_mul(de.unsigned_abs() + am.bits().max(bm.bits()) as u128);
        if size >= 10000 {
            return Err(PolyError::Unsupported("large complex powers through mpc_exp(mpc_log) are not implemented".into()));
        }
        if de > 0 {
            am <<= de.to_usize().unwrap();
            ae = be.clone();
        } else {
            bm <<= (-&de).to_usize().unwrap();
            be = ae.clone();
        }
        let (mut wre, mut wim) = (BigInt::one(), BigInt::zero());
        let (mut x, mut y) = (am, bm);
        let mut k = n;
        while k > 0 {
            if k & 1 == 1 {
                let r = &wre * &x - &wim * &y;
                let i = &wim * &x + &wre * &y;
                wre = r;
                wim = i;
                k -= 1;
            }
            let x2 = &x * &x - &y * &y;
            let y2 = BigInt::from(2) * &x * &y;
            x = x2;
            y = y2;
            k /= 2;
        }
        Ok(C { re: BF::round(wre, &ae * n, p), im: BF::round(wim, &be * n, p) })
    }
    fn is_negative(&self, a: &C) -> bool {
        !self.complex && a.re.is_negative()
    }
    fn is_nonnegative(&self, a: &C) -> bool {
        !self.complex && !a.re.is_negative()
    }
    /// `Domain.canonical_unit`: `-1` for a negative real number.
    fn canonical_unit(&self, a: &C) -> C {
        if self.is_negative(a) {
            self.neg(&self.one())
        } else {
            self.one()
        }
    }
    fn max_denom(&self, prec: u32) -> BigInt {
        ((BigInt::one() << prec as usize) / 200u32).max(BigInt::from(99))
    }
    /// `RealField.to_rational(a, limit=True)` with the denominator bound of
    /// a context of precision `prec`.
    fn to_rational(&self, x: &BF, prec: u32) -> Rat {
        let r = x.to_ratio_exact();
        let (p, q) = (r.numer().clone(), r.denom().clone());
        let max_denom = self.max_denom(prec);
        let rat = |a: BigInt, b: BigInt| Rat::new(Int::from_big(a), Int::from_big(b));
        if q <= max_denom {
            return rat(p, q);
        }
        let (mut p0, mut q0, mut p1, mut q1) = (BigInt::zero(), BigInt::one(), BigInt::one(), BigInt::zero());
        let (mut n, mut d) = (p.clone(), q.clone());
        loop {
            let a = num_integer::Integer::div_floor(&n, &d);
            let q2 = &q0 + &a * &q1;
            if q2 > max_denom {
                break;
            }
            let p2 = &p0 + &a * &p1;
            p0 = std::mem::replace(&mut p1, p2);
            q0 = std::mem::replace(&mut q1, q2);
            let rem = &n - &a * &d;
            n = std::mem::replace(&mut d, rem);
        }
        let k = num_integer::Integer::div_floor(&(&max_denom - &q0), &q1);
        let number = rat(p, q);
        let bound1 = rat(&p0 + &k * &p1, &q0 + &k * &q1);
        let bound2 = rat(p1, q1);
        if bound2.is_zero() || bound1.is_zero() {
            return number;
        }
        let d2 = bound2.sub(&number).abs();
        let d1 = bound1.sub(&number).abs();
        if !d1.sub(&d2).is_negative() {
            bound2
        } else {
            bound1
        }
    }
    /// `QQ.convert`/`QQ_I.convert` of a scalar (the parts of a complex
    /// number are rationalized as plain `mpf` objects, with the default
    /// 53-bit context bound).
    fn exact_parts(&self, a: &C) -> (Rat, Rat) {
        if self.complex {
            (self.to_rational(&a.re, 53), self.to_rational(&a.im, 53))
        } else {
            (self.to_rational(&a.re, self.p()), Rat::zero())
        }
    }
    /// `to_sympy`: `Float(re, dps) + I*Float(im, dps)`.
    fn to_expr(&self, a: &C) -> Expr {
        let fp = self.float_prec();
        let f = |x: &BF| Number::mp_float(x.round_to(fp), fp);
        if !self.complex {
            return Expr::num(f(&a.re));
        }
        let re = if a.re.is_zero() { None } else { Some(Expr::num(f(&a.re))) };
        let im = if a.im.is_zero() { None } else { Some(Expr::from_kind(Kind::Mul(f(&a.im), vec![(Expr::i(), Expr::one())]))) };
        match (re, im) {
            (Some(r), Some(i)) => crate::ops::add(&[r, i]),
            (Some(r), None) => r,
            (None, Some(i)) => i,
            (None, None) => Expr::num(f(&BF::zero())),
        }
    }
    /// `from_sympy` of a number: `evalf(n=dps)` then the context rounding.
    fn from_number(&self, e: &Expr) -> Option<C> {
        if !crate::traverse::free_symbols(e).is_empty() {
            return None;
        }
        let v = crate::evalf_mp::eval_mp_default(e, self.float_prec())?;
        let c = C { re: v.re.round_to(self.p()), im: v.im.round_to(self.p()) };
        if !self.complex && !c.im.is_zero() {
            return None;
        }
        Some(c)
    }
}

// ----------------------------------------------------------------------
// ground polynomials: insertion-ordered sparse dictionaries
// ----------------------------------------------------------------------

#[derive(Clone, Debug, Default)]
pub struct R {
    slots: Vec<Option<(Vec<u32>, C)>>,
    index: FxHashMap<Vec<u32>, usize>,
}

impl PartialEq for R {
    fn eq(&self, o: &Self) -> bool {
        self.len() == o.len() && self.iter().all(|(k, v)| o.get(k) == Some(v))
    }
}

impl R {
    fn len(&self) -> usize {
        self.index.len()
    }
    fn is_empty(&self) -> bool {
        self.index.is_empty()
    }
    fn get(&self, k: &[u32]) -> Option<&C> {
        self.index.get(k).and_then(|&i| self.slots[i].as_ref().map(|s| &s.1))
    }
    fn set(&mut self, k: Vec<u32>, c: C) {
        if let Some(&i) = self.index.get(&k) {
            self.slots[i].as_mut().unwrap().1 = c;
        } else {
            self.index.insert(k.clone(), self.slots.len());
            self.slots.push(Some((k, c)));
        }
    }
    fn del(&mut self, k: &[u32]) {
        if let Some(i) = self.index.remove(k) {
            self.slots[i] = None;
        }
    }
    fn iter(&self) -> impl Iterator<Item = (&Vec<u32>, &C)> + '_ {
        self.slots.iter().filter_map(|s| s.as_ref().map(|(k, c)| (k, c)))
    }
    /// `p.copy()`, dropping tombstones.
    fn copy(&self) -> R {
        let mut r = R::default();
        for (k, c) in self.iter() {
            r.set(k.clone(), c.clone());
        }
        r
    }
    fn from_iter(terms: impl IntoIterator<Item = (Vec<u32>, C)>) -> R {
        let mut r = R::default();
        for (k, c) in terms {
            r.set(k, c);
        }
        r
    }
    fn leading_expv(&self) -> Option<Vec<u32>> {
        self.iter().map(|(k, _)| k).max().cloned()
    }
}

/// The polynomial ring `K[y...]` over an inexact scalar domain.
#[derive(Clone, Debug)]
struct Ring {
    k: Inexact,
    m: usize,
    symbols: Vec<Expr>,
}

fn mono_add(a: &[u32], b: &[u32]) -> Vec<u32> {
    a.iter().zip(b).map(|(x, y)| x + y).collect()
}

/// `binomial_coefficients(n)` in the reference dictionary order.
fn binomial_coefficients(n: usize) -> Vec<(Vec<usize>, BigInt)> {
    let mut d: Vec<(Vec<usize>, BigInt)> = vec![(vec![0, n], BigInt::one()), (vec![n, 0], BigInt::one())];
    let mut a = BigInt::one();
    for k in 1..=n / 2 {
        a = a * BigInt::from(n - k + 1) / BigInt::from(k);
        for key in [vec![k, n - k], vec![n - k, k]] {
            if let Some(slot) = d.iter_mut().find(|(x, _)| *x == key) {
                slot.1 = a.clone();
            } else {
                d.push((key, a.clone()));
            }
        }
    }
    d
}

/// The reference `multinomial_coefficients(m, n)` in dictionary order.
pub(crate) fn multinomial_coefficients(m: usize, n: usize) -> Vec<(Vec<usize>, BigInt)> {
    if m == 0 {
        return if n == 0 { vec![(vec![], BigInt::one())] } else { vec![] };
    }
    if m == 2 {
        return binomial_coefficients(n);
    }
    if m >= 2 * n && n > 1 {
        // `multinomial_coefficients_iterator`
        let mc = multinomial_coefficients(n, n);
        let lookup = |t: &[usize]| -> BigInt {
            let mut key: Vec<usize> = t.iter().copied().filter(|&x| x != 0).collect();
            key.resize(n, 0);
            mc.iter().find(|(k, _)| *k == key).map(|(_, v)| v.clone()).unwrap_or_else(BigInt::zero)
        };
        let mut out = Vec::new();
        let mut t = vec![0usize; m];
        t[0] = n;
        out.push((t.clone(), lookup(&t)));
        let mut j = if n > 0 { 0 } else { m };
        while j + 1 < m {
            let tj = t[j];
            if j > 0 {
                t[j] = 0;
                t[0] = tj;
            }
            if tj > 1 {
                t[j + 1] += 1;
                j = 0;
            } else {
                j += 1;
                t[j] += 1;
            }
            t[0] -= 1;
            out.push((t.clone(), lookup(&t)));
        }
        return out;
    }
    let mut t = vec![0usize; m];
    t[0] = n;
    let mut r: Vec<(Vec<usize>, BigInt)> = vec![(t.clone(), BigInt::one())];
    let mut idx: FxHashMap<Vec<usize>, usize> = FxHashMap::default();
    idx.insert(t.clone(), 0);
    let mut j = if n > 0 { 0 } else { m };
    while j + 1 < m {
        let tj = t[j];
        if j > 0 {
            t[j] = 0;
            t[0] = tj;
        }
        let (start, mut v);
        if tj > 1 {
            t[j + 1] += 1;
            j = 0;
            start = 1;
            v = BigInt::zero();
        } else {
            j += 1;
            start = j + 1;
            v = r[idx[&t]].1.clone();
            t[j] += 1;
        }
        for k in start..m {
            if t[k] > 0 {
                t[k] -= 1;
                v += &r[idx[&t]].1;
                t[k] += 1;
            }
        }
        t[0] -= 1;
        let value = v * BigInt::from(tj) / BigInt::from(n - t[0]);
        if let Some(&i) = idx.get(&t) {
            r[i].1 = value;
        } else {
            idx.insert(t.clone(), r.len());
            r.push((t.clone(), value));
        }
    }
    r
}

impl Ring {
    fn zm(&self) -> Vec<u32> {
        vec![0; self.m]
    }
    fn zero(&self) -> R {
        R::default()
    }
    fn one(&self) -> R {
        self.ground(self.k.one())
    }
    fn ground(&self, c: C) -> R {
        let mut r = R::default();
        if !self.k.is_zero(&c) {
            r.set(self.zm(), c);
        }
        r
    }
    fn gen(&self, i: usize) -> R {
        let mut e = self.zm();
        e[i] = 1;
        R::from_iter([(e, self.k.one())])
    }
    fn is_one(&self, a: &R) -> bool {
        a.len() == 1 && a.get(&self.zm()) == Some(&self.k.one())
    }
    fn lc(&self, a: &R) -> C {
        a.leading_expv().and_then(|e| a.get(&e).cloned()).unwrap_or_else(|| self.k.zero())
    }
    fn add(&self, a: &R, b: &R) -> R {
        let mut p = a.copy();
        for (k, v) in b.iter() {
            let c = match p.get(k) {
                Some(x) => self.k.add(x, v),
                None => v.clone(),
            };
            if self.k.is_zero(&c) {
                p.del(k);
            } else {
                p.set(k.clone(), c);
            }
        }
        p
    }
    fn sub(&self, a: &R, b: &R) -> R {
        let mut p = a.copy();
        for (k, v) in b.iter() {
            let c = match p.get(k) {
                Some(x) => self.k.sub(x, v),
                None => self.k.neg(v),
            };
            if self.k.is_zero(&c) {
                p.del(k);
            } else {
                p.set(k.clone(), c);
            }
        }
        p
    }
    fn neg(&self, a: &R) -> R {
        R::from_iter(a.iter().map(|(k, v)| (k.clone(), self.k.neg(v))))
    }
    fn strip(&self, p: R) -> R {
        if p.iter().any(|(_, v)| self.k.is_zero(v)) {
            R::from_iter(p.iter().filter(|(_, v)| !self.k.is_zero(v)).map(|(k, v)| (k.clone(), v.clone())))
        } else {
            p
        }
    }
    fn acc(&self, p: &mut R, e: Vec<u32>, v: C) {
        let c = match p.get(&e) {
            Some(x) => self.k.add(x, &v),
            None => v,
        };
        p.set(e, c);
    }
    fn mul(&self, a: &R, b: &R) -> R {
        if a.is_empty() || b.is_empty() {
            return self.zero();
        }
        let mut p = R::default();
        for (e1, v1) in a.iter() {
            for (e2, v2) in b.iter() {
                self.acc(&mut p, mono_add(e1, e2), self.k.mul(v1, v2));
            }
        }
        self.strip(p)
    }
    fn square(&self, a: &R) -> PResult<R> {
        let keys: Vec<&Vec<u32>> = a.iter().map(|(k, _)| k).collect();
        let mut p = R::default();
        for i in 0..keys.len() {
            let k1 = keys[i];
            let pk = a.get(k1).unwrap();
            for k2 in keys.iter().take(i) {
                self.acc(&mut p, mono_add(k1, k2), self.k.mul(pk, a.get(k2).unwrap()));
            }
        }
        let mut p = R::from_iter(p.iter().map(|(k, v)| (k.clone(), self.k.add(v, v))));
        for (k, v) in a.iter() {
            self.acc(&mut p, mono_add(k, k), self.k.pow(v, 2)?);
        }
        Ok(self.strip(p))
    }
    fn pow(&self, a: &R, n: u64) -> PResult<R> {
        if n == 0 {
            return Ok(self.one());
        }
        if a.len() == 1 {
            let (monom, coeff) = a.iter().next().unwrap();
            let e: Vec<u32> = monom.iter().map(|x| x * n as u32).collect();
            let c = if *coeff == self.k.one() { coeff.clone() } else { self.k.pow(coeff, n)? };
            return Ok(R::from_iter([(e, c)]));
        }
        match n {
            1 => Ok(a.copy()),
            2 => self.square(a),
            3 => Ok(self.mul(a, &self.square(a)?)),
            _ if a.len() <= 5 => self.pow_multinomial(a, n),
            _ => {
                let mut p = self.one();
                let mut c = a.clone();
                let mut n = n;
                loop {
                    if n & 1 == 1 {
                        p = self.mul(&p, &c);
                        n -= 1;
                        if n == 0 {
                            break;
                        }
                    }
                    c = self.square(&c)?;
                    n /= 2;
                }
                Ok(p)
            }
        }
    }
    fn pow_multinomial(&self, a: &R, n: u64) -> PResult<R> {
        let terms: Vec<(&Vec<u32>, &C)> = a.iter().collect();
        let mut poly = R::default();
        for (multinomial, mc) in multinomial_coefficients(terms.len(), n as usize) {
            let mut monom = self.zm();
            let mut coeff: Option<C> = None;
            for (&exp, (m, c)) in multinomial.iter().zip(&terms) {
                if exp > 0 {
                    for (x, y) in monom.iter_mut().zip(m.iter()) {
                        *x += y * exp as u32;
                    }
                    let t = self.k.pow(c, exp as u64)?;
                    coeff = Some(match coeff {
                        None => self.k.mul_int(&t, &mc),
                        Some(v) => self.k.mul(&v, &t),
                    });
                }
            }
            let coeff = coeff.unwrap_or_else(|| self.k.from_int(&mc));
            let c = match poly.get(&monom) {
                Some(x) => self.k.add(x, &coeff),
                None => coeff,
            };
            if !self.k.is_zero(&c) {
                poly.set(monom, c);
            } else {
                poly.del(&monom);
            }
        }
        Ok(poly)
    }
    fn mul_ground(&self, a: &R, x: &C) -> R {
        if self.k.is_zero(x) {
            return self.zero();
        }
        self.strip(R::from_iter(a.iter().map(|(k, v)| (k.clone(), self.k.mul(v, x)))))
    }
    fn quo_ground(&self, a: &R, x: &C) -> PResult<R> {
        if self.k.is_zero(x) {
            return Err(zero_division());
        }
        if a.is_empty() || *x == self.k.one() {
            return Ok(a.clone());
        }
        let terms = a.iter().map(|(k, v)| Ok((k.clone(), self.k.div(v, x)?))).collect::<PResult<Vec<_>>>()?;
        Ok(self.strip(R::from_iter(terms)))
    }
    fn monic(&self, a: &R) -> PResult<R> {
        if a.is_empty() {
            return Ok(a.clone());
        }
        self.quo_ground(a, &self.lc(a))
    }
    fn term_div(&self, am: &[u32], bm: &[u32]) -> Option<Vec<u32>> {
        if bm.iter().all(|&x| x == 0) {
            return Some(am.to_vec());
        }
        am.iter().zip(bm).map(|(a, b)| a.checked_sub(*b)).collect()
    }
    /// `PolyElement.div` by a single divisor.
    fn div(&self, a: &R, g: &R) -> PResult<(R, R)> {
        if g.is_empty() {
            return Err(zero_division());
        }
        if a.is_empty() {
            return Ok((self.zero(), self.zero()));
        }
        let gl = g.leading_expv().unwrap();
        let gc = g.get(&gl).unwrap().clone();
        let mut p = a.copy();
        let (mut q, mut r) = (R::default(), R::default());
        let mut steps = 0usize;
        while !p.is_empty() {
            steps += 1;
            if steps > 100_000 {
                return Err(PolyError::Unsupported(
                    "inexact polynomial division does not terminate (the reference loops forever)".into(),
                ));
            }
            let expv = p.leading_expv().unwrap();
            let lc = p.get(&expv).unwrap().clone();
            match self.term_div(&expv, &gl) {
                Some(m) => {
                    let c = self.k.div(&lc, &gc)?;
                    iadd_monom(&self.k, &mut q, m.clone(), c.clone());
                    let mc = self.k.neg(&c);
                    for (k, v) in g.iter() {
                        let ka = mono_add(k, &m);
                        let coeff = match p.get(&ka) {
                            Some(x) => self.k.add(x, &self.k.mul(v, &mc)),
                            None => self.k.mul(v, &mc),
                        };
                        if self.k.is_zero(&coeff) {
                            p.del(&ka);
                        } else {
                            p.set(ka, coeff);
                        }
                    }
                }
                None => {
                    iadd_monom(&self.k, &mut r, expv.clone(), lc);
                    p.del(&expv);
                }
            }
        }
        Ok((q, r))
    }
    fn quo(&self, a: &R, g: &R) -> PResult<R> {
        Ok(self.div(a, g)?.0)
    }
    fn rem(&self, a: &R, g: &R) -> PResult<R> {
        Ok(self.div(a, g)?.1)
    }
    fn exquo(&self, a: &R, g: &R) -> PResult<R> {
        let (q, r) = self.div(a, g)?;
        if r.is_empty() {
            Ok(q)
        } else {
            Err(PolyError::ExactQuotientFailed(format!("{} does not divide {}", self.to_expr(g, &self.symbols), self.to_expr(a, &self.symbols))))
        }
    }
    fn is_nonnegative(&self, a: &R) -> bool {
        self.k.is_nonnegative(&self.lc(a))
    }
    fn is_negative(&self, a: &R) -> bool {
        self.k.is_negative(&self.lc(a))
    }
    fn gcd_zero(&self, g: &R) -> (R, R, R) {
        if self.is_nonnegative(g) {
            (g.clone(), self.zero(), self.one())
        } else {
            (self.neg(g), self.zero(), self.neg(&self.one()))
        }
    }
    fn gcd_monom(&self, f: &R, g: &R) -> PResult<(R, R, R)> {
        let (mf, cf) = f.iter().next().map(|(k, v)| (k.clone(), v.clone())).unwrap();
        let mut mgcd = mf.clone();
        for (mg, _) in g.iter() {
            for (x, y) in mgcd.iter_mut().zip(mg) {
                *x = (*x).min(*y);
            }
        }
        let one = self.k.one();
        let sub = |a: &[u32]| -> Vec<u32> { a.iter().zip(&mgcd).map(|(x, y)| x - y).collect() };
        let h = R::from_iter([(mgcd.clone(), one.clone())]);
        let cff = R::from_iter([(sub(&mf), self.k.div(&cf, &one)?)]);
        let cfg = R::from_iter(g.iter().map(|(mg, cg)| Ok((sub(mg), self.k.div(cg, &one)?))).collect::<PResult<Vec<_>>>()?);
        Ok((h, cff, cfg))
    }
    /// `PolyElement.cofactors`.
    fn cofactors(&self, f: &R, g: &R) -> PResult<(R, R, R)> {
        if f.is_empty() && g.is_empty() {
            return Ok((self.zero(), self.zero(), self.zero()));
        }
        if f.is_empty() {
            return Ok(self.gcd_zero(g));
        }
        if g.is_empty() {
            let (h, cfg, cff) = self.gcd_zero(f);
            return Ok((h, cff, cfg));
        }
        if f.len() == 1 {
            return self.gcd_monom(f, g);
        }
        if g.len() == 1 {
            let (h, cfg, cff) = self.gcd_monom(g, f)?;
            return Ok((h, cff, cfg));
        }
        self.exact_gcd(f, g)
    }
    /// `dmp_inner_gcd` over `RR`/`CC`: through `QQ`/`QQ_I` and back.
    fn exact_gcd(&self, f: &R, g: &R) -> PResult<(R, R, R)> {
        let m = self.m;
        if !self.k.complex {
            let conv = |p: &R| -> MPoly<Rat> {
                MPoly::from_terms(m, p.iter().map(|(e, c)| (e.clone(), self.k.exact_parts(c).0)).filter(|(_, c)| !c.is_zero()).collect())
            };
            let (h, a, b) = super::gcd::gcd_mq(&conv(f), &conv(g));
            let back = |p: &MPoly<Rat>| -> R {
                let mut terms: Vec<(Vec<u32>, C)> = p.iter().map(|(e, c)| (e.to_vec(), self.k.from_rat(c))).collect();
                terms.reverse();
                self.strip(R::from_iter(terms))
            };
            return Ok((back(&h), back(&a), back(&b)));
        }
        let _s = AlgebraicScope::new(gi::field());
        let conv = |p: &R| -> MPoly<AN> {
            MPoly::from_terms(
                m,
                p.iter()
                    .map(|(e, c)| {
                        let (x, y) = self.k.exact_parts(c);
                        (e.clone(), gi::make(x, y))
                    })
                    .filter(|(_, c)| !c.coefficients.iter().all(Rat::is_zero))
                    .collect(),
            )
        };
        let (fe, ge) = (conv(f), conv(g));
        let (h, a, b) = if fe.is_zero() || ge.is_zero() {
            return Err(PolyError::Unsupported("degenerate complex gcd".into()));
        } else {
            let h = super::algebraic_multi::gcd(&fe, &ge)?;
            let a = super::algebraic_multi::exact(&fe, &h)?;
            let b = super::algebraic_multi::exact(&ge, &h)?;
            (h, a, b)
        };
        let back = |p: &MPoly<AN>| -> R {
            let mut terms: Vec<(Vec<u32>, C)> = p
                .iter()
                .map(|(e, c)| {
                    let (x, y) = gi::parts(c);
                    (e.to_vec(), self.k.from_parts(&x, &y))
                })
                .collect();
            terms.reverse();
            self.strip(R::from_iter(terms))
        };
        Ok((back(&h), back(&a), back(&b)))
    }
    fn gcd(&self, f: &R, g: &R) -> PResult<R> {
        Ok(self.cofactors(f, g)?.0)
    }
    /// `PolyElement.lcm` over a field ground.
    fn lcm(&self, f: &R, g: &R) -> PResult<R> {
        let h = self.quo(&self.mul(f, g), &self.gcd(f, g)?)?;
        self.monic(&h)
    }
    /// `PolyElement.cancel` (the ground has no associated ring).
    fn cancel(&self, f: &R, g: &R) -> PResult<(R, R)> {
        if f.is_empty() {
            return Ok((f.clone(), self.one()));
        }
        let (_, p, q) = self.cofactors(f, g)?;
        if self.is_negative(&q) {
            return Ok((self.neg(&p), self.neg(&q)));
        }
        Ok((p, q))
    }
    fn to_expr(&self, a: &R, symbols: &[Expr]) -> Expr {
        let terms: Vec<Expr> = a
            .iter()
            .map(|(e, c)| {
                let mut f = vec![self.k.to_expr(c)];
                for (s, &k) in symbols.iter().zip(e) {
                    if k > 0 {
                        f.push(crate::ops::pow(s, &Expr::int(k as i64)));
                    }
                }
                crate::ops::mul(&f)
            })
            .collect();
        crate::ops::add(&terms)
    }
}

fn iadd_monom(k: &Inexact, p: &mut R, e: Vec<u32>, c: C) {
    match p.get(&e).cloned() {
        None => p.set(e, c),
        Some(x) => {
            let s = k.add(&x, &c);
            if k.is_zero(&s) {
                p.del(&e);
            } else {
                p.set(e, s);
            }
        }
    }
}

// ----------------------------------------------------------------------
// ground elements: PolyElement (ring) or FracElement (field)
// ----------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
pub struct E {
    n: R,
    d: R,
}

/// The composite ground domain `K[y...]` or `K(y...)`.
#[derive(Clone, Debug)]
pub struct Dom {
    r: Ring,
    field: bool,
    symbols: Vec<Expr>,
    ground: GaussianGround,
}

type Key = (Expr, Inexact, bool, Vec<Expr>);
thread_local! {
    static CACHE: RefCell<FxHashMap<Key, E>> = RefCell::new(FxHashMap::default());
}

impl Dom {
    pub fn new(g: &GaussianGround) -> Option<Dom> {
        let k = g.inexact?;
        Some(Dom { r: Ring { k, m: g.symbols.len(), symbols: g.symbols.clone() }, field: !g.ring, symbols: g.symbols.clone(), ground: g.clone() })
    }
    fn k(&self) -> &Inexact {
        &self.r.k
    }
    fn display(&self) -> String {
        self.ground.display()
    }
    fn as_ring(&self) -> Dom {
        Dom { field: false, ground: self.ground.to_ring(), ..self.clone() }
    }
    fn as_field(&self) -> Dom {
        Dom { field: true, ground: self.ground.to_field(), ..self.clone() }
    }
    fn raw(&self, n: R) -> E {
        E { n, d: self.r.one() }
    }
    fn zero(&self) -> E {
        self.raw(self.r.zero())
    }
    fn one(&self) -> E {
        self.raw(self.r.one())
    }
    fn from_c(&self, c: C) -> E {
        self.raw(self.r.ground(c))
    }
    fn from_int(&self, n: i64) -> E {
        self.from_c(self.k().from_int(&BigInt::from(n)))
    }
    fn is_zero(&self, a: &E) -> bool {
        a.n.is_empty()
    }
    fn is_one(&self, a: &E) -> bool {
        self.r.is_one(&a.n) && self.r.is_one(&a.d)
    }
    /// `FracElement.new`: `raw_new(*numer.cancel(denom))`.
    fn frac(&self, n: R, d: R) -> PResult<E> {
        let (n, d) = self.r.cancel(&n, &d)?;
        Ok(E { n, d })
    }
    fn add(&self, a: &E, b: &E) -> PResult<E> {
        if !self.field {
            return Ok(self.raw(self.r.add(&a.n, &b.n)));
        }
        if self.is_zero(b) {
            return Ok(a.clone());
        }
        if self.is_zero(a) {
            return Ok(b.clone());
        }
        if a.d == b.d {
            return self.frac(self.r.add(&a.n, &b.n), a.d.clone());
        }
        self.frac(self.r.add(&self.r.mul(&a.n, &b.d), &self.r.mul(&a.d, &b.n)), self.r.mul(&a.d, &b.d))
    }
    fn sub(&self, a: &E, b: &E) -> PResult<E> {
        if !self.field {
            return Ok(self.raw(self.r.sub(&a.n, &b.n)));
        }
        if self.is_zero(b) {
            return Ok(a.clone());
        }
        if self.is_zero(a) {
            return Ok(self.neg(b));
        }
        if a.d == b.d {
            return self.frac(self.r.sub(&a.n, &b.n), a.d.clone());
        }
        self.frac(self.r.sub(&self.r.mul(&a.n, &b.d), &self.r.mul(&a.d, &b.n)), self.r.mul(&a.d, &b.d))
    }
    fn neg(&self, a: &E) -> E {
        E { n: self.r.neg(&a.n), d: a.d.clone() }
    }
    fn mul(&self, a: &E, b: &E) -> PResult<E> {
        if !self.field {
            return Ok(self.raw(self.r.mul(&a.n, &b.n)));
        }
        if self.is_zero(a) || self.is_zero(b) {
            return Ok(self.zero());
        }
        self.frac(self.r.mul(&a.n, &b.n), self.r.mul(&a.d, &b.d))
    }
    fn pow(&self, a: &E, n: i64) -> PResult<E> {
        if !self.field {
            if n < 0 {
                return Err(reference("ValueError", "exponent must be a non-negative integer"));
            }
            return Ok(self.raw(self.r.pow(&a.n, n as u64)?));
        }
        if n >= 0 {
            return Ok(E { n: self.r.pow(&a.n, n as u64)?, d: self.r.pow(&a.d, n as u64)? });
        }
        if self.is_zero(a) {
            return Err(reference("ZeroDivisionError", ""));
        }
        let k = n.unsigned_abs();
        Ok(E { n: self.r.pow(&a.d, k)?, d: self.r.pow(&a.n, k)? })
    }
    fn sqr(&self, a: &E) -> PResult<E> {
        self.pow(a, 2)
    }
    fn div_field(&self, a: &E, b: &E) -> PResult<E> {
        if self.is_zero(b) {
            return Err(reference("ZeroDivisionError", ""));
        }
        self.frac(self.r.mul(&a.n, &b.d), self.r.mul(&a.d, &b.n))
    }
    /// `K.quo`: `a // b` in the ring, `a / b` in the field.
    fn quo(&self, a: &E, b: &E) -> PResult<E> {
        if self.field {
            return self.div_field(a, b);
        }
        Ok(self.raw(self.r.quo(&a.n, &b.n)?))
    }
    fn exquo(&self, a: &E, b: &E) -> PResult<E> {
        if self.field {
            return self.div_field(a, b);
        }
        Ok(self.raw(self.r.exquo(&a.n, &b.n)?))
    }
    fn rem(&self, a: &E, b: &E) -> PResult<E> {
        if self.field {
            return Ok(self.zero());
        }
        Ok(self.raw(self.r.rem(&a.n, &b.n)?))
    }
    fn gcd(&self, a: &E, b: &E) -> PResult<E> {
        if !self.field {
            return Ok(self.raw(self.r.gcd(&a.n, &b.n)?));
        }
        let p = self.r.gcd(&a.n, &b.n)?;
        let q = self.r.lcm(&a.d, &b.d)?;
        self.frac(p, q)
    }
    fn lcm(&self, a: &E, b: &E) -> PResult<E> {
        if !self.field {
            return Ok(self.raw(self.r.lcm(&a.n, &b.n)?));
        }
        let p = self.r.lcm(&a.n, &b.n)?;
        let q = self.r.gcd(&a.d, &b.d)?;
        self.frac(p, q)
    }
    /// `Domain.cofactors`.
    fn cofactors(&self, a: &E, b: &E) -> PResult<(E, E, E)> {
        let g = self.gcd(a, b)?;
        Ok((g.clone(), self.quo(a, &g)?, self.quo(b, &g)?))
    }
    fn canonical_unit(&self, a: &E) -> E {
        self.from_c(self.k().canonical_unit(&self.r.lc(&a.n)))
    }
    fn is_negative(&self, a: &E) -> bool {
        self.r.is_negative(&a.n)
    }
    fn is_nonnegative(&self, a: &E) -> bool {
        self.r.is_nonnegative(&a.n)
    }

    // ---- expressions ----

    fn coercion(&self, e: &Expr) -> PolyError {
        GaussianGround::coercion_error(&self.ground, e)
    }
    pub fn to_expr(&self, a: &E) -> Expr {
        let e = if self.field {
            crate::ops::div(&self.r.to_expr(&a.n, &self.symbols), &self.r.to_expr(&a.d, &self.symbols))
        } else {
            self.r.to_expr(&a.n, &self.symbols)
        };
        // A fraction with a constant denominator prints like several distinct
        // representations (`0.5/0.5` and `1.0/1`), so only exact round trips
        // are remembered.
        if self.field && a.d.iter().all(|(k, _)| k.iter().all(|&x| x == 0)) {
            return e;
        }
        CACHE.with(|c| {
            let mut c = c.borrow_mut();
            if c.len() > 50_000 {
                c.clear();
            }
            c.insert((e.clone(), *self.k(), self.field, self.symbols.clone()), a.clone());
        });
        e
    }
    /// `from_sympy`: the reference `_rebuild_expr`.
    pub fn from_expr(&self, e: &Expr) -> PResult<E> {
        if let Some(v) = CACHE.with(|c| c.borrow().get(&(e.clone(), *self.k(), self.field, self.symbols.clone())).cloned()) {
            return Ok(v);
        }
        self.rebuild(e)
    }
    fn rebuild(&self, e: &Expr) -> PResult<E> {
        if let Some(i) = self.symbols.iter().position(|s| s == e) {
            return Ok(self.raw(self.r.gen(i)));
        }
        match e.kind() {
            Kind::Add(..) => {
                let args = e.args();
                let mut acc = self.rebuild(&args[0])?;
                for a in &args[1..] {
                    acc = self.add(&acc, &self.rebuild(a)?)?;
                }
                return Ok(acc);
            }
            Kind::Mul(..) => {
                let args = e.args();
                let mut acc = self.rebuild(&args[0])?;
                for a in &args[1..] {
                    acc = self.mul(&acc, &self.rebuild(a)?)?;
                }
                return Ok(acc);
            }
            Kind::Pow(b, x) => {
                if let Some(n) = x.as_number().filter(|n| n.is_integer()).and_then(Number::to_i64) {
                    if (!self.field && n > 1) || (self.field && n != 1) {
                        return self.pow(&self.rebuild(b)?, n);
                    }
                }
            }
            _ => {}
        }
        match self.k().from_number(e) {
            Some(c) => Ok(self.from_c(c)),
            None => Err(self.coercion(e)),
        }
    }
    /// `from_FractionField` into the ring: the numerator divided by the
    /// denominator when the division is exact.
    fn retract(&self, a: &E) -> Option<E> {
        let (q, r) = self.r.div(&a.n, &a.d).ok()?;
        r.is_empty().then(|| self.raw(q))
    }
}

// ----------------------------------------------------------------------
// dense univariate polynomials over the composite domain
// ----------------------------------------------------------------------

type U = Vec<E>;

fn degree(f: &[E]) -> i64 {
    f.len() as i64 - 1
}

impl Dom {
    fn strip(&self, mut f: U) -> U {
        let k = f.iter().take_while(|c| self.is_zero(c)).count();
        f.drain(..k);
        f
    }
    fn lc_of(&self, f: &[E]) -> E {
        f.first().cloned().unwrap_or_else(|| self.zero())
    }
    fn tc_of(&self, f: &[E]) -> E {
        f.last().cloned().unwrap_or_else(|| self.zero())
    }
    fn d_neg(&self, f: &[E]) -> U {
        f.iter().map(|c| self.neg(c)).collect()
    }
    fn d_add(&self, f: &[E], g: &[E]) -> PResult<U> {
        if f.is_empty() {
            return Ok(g.to_vec());
        }
        if g.is_empty() {
            return Ok(f.to_vec());
        }
        let (df, dg) = (degree(f), degree(g));
        if df == dg {
            return Ok(self.strip(f.iter().zip(g).map(|(a, b)| self.add(a, b)).collect::<PResult<_>>()?));
        }
        let k = (df - dg).unsigned_abs() as usize;
        let (mut h, f, g) = if df > dg { (f[..k].to_vec(), &f[k..], g) } else { (g[..k].to_vec(), f, &g[k..]) };
        for (a, b) in f.iter().zip(g) {
            h.push(self.add(a, b)?);
        }
        Ok(h)
    }
    fn d_sub(&self, f: &[E], g: &[E]) -> PResult<U> {
        if f.is_empty() {
            return Ok(self.d_neg(g));
        }
        if g.is_empty() {
            return Ok(f.to_vec());
        }
        let (df, dg) = (degree(f), degree(g));
        if df == dg {
            return Ok(self.strip(f.iter().zip(g).map(|(a, b)| self.sub(a, b)).collect::<PResult<_>>()?));
        }
        let k = (df - dg).unsigned_abs() as usize;
        let (mut h, f, g) = if df > dg { (f[..k].to_vec(), &f[k..], g) } else { (self.d_neg(&g[..k]), f, &g[k..]) };
        for (a, b) in f.iter().zip(g) {
            h.push(self.sub(a, b)?);
        }
        Ok(h)
    }
    fn d_mul(&self, f: &[E], g: &[E]) -> PResult<U> {
        if f == g {
            return self.d_sqr(f);
        }
        if f.is_empty() || g.is_empty() {
            return Ok(Vec::new());
        }
        let (df, dg) = (f.len() - 1, g.len() - 1);
        let mut h = Vec::with_capacity(df + dg + 1);
        for i in 0..=df + dg {
            let mut coeff = self.zero();
            for j in i.saturating_sub(dg)..=df.min(i) {
                coeff = self.add(&coeff, &self.mul(&f[j], &g[i - j])?)?;
            }
            h.push(coeff);
        }
        Ok(self.strip(h))
    }
    fn d_sqr(&self, f: &[E]) -> PResult<U> {
        if f.is_empty() {
            return Ok(Vec::new());
        }
        let df = f.len() - 1;
        let mut h = Vec::new();
        for i in 0..=2 * df {
            let mut c = self.zero();
            let jmin = i.saturating_sub(df);
            let jmax0 = i.min(df);
            let n = jmax0 - jmin + 1;
            let jmax = (jmin + n / 2) as i64 - 1;
            let mut j = jmin as i64;
            while j <= jmax {
                c = self.add(&c, &self.mul(&f[j as usize], &f[i - j as usize])?)?;
                j += 1;
            }
            c = self.add(&c, &c)?;
            if n & 1 == 1 {
                let elem = &f[(jmax + 1) as usize];
                c = self.add(&c, &self.sqr(elem)?)?;
            }
            h.push(c);
        }
        Ok(self.strip(h))
    }
    fn d_pow(&self, f: &[E], n: u32) -> PResult<U> {
        if n == 0 {
            return Ok(vec![self.one()]);
        }
        if n == 1 || f.is_empty() || (f.len() == 1 && self.is_one(&f[0])) {
            return Ok(f.to_vec());
        }
        let mut g = vec![self.one()];
        let mut f = f.to_vec();
        let mut n = n;
        loop {
            let m = n;
            n /= 2;
            if m % 2 == 1 {
                g = self.d_mul(&g, &f)?;
                if n == 0 {
                    break;
                }
            }
            f = self.d_sqr(&f)?;
        }
        Ok(g)
    }
    fn add_term(&self, f: &[E], c: &E, i: usize) -> PResult<U> {
        if self.is_zero(c) {
            return Ok(f.to_vec());
        }
        let n = f.len();
        if n > 0 && i == n - 1 {
            let mut h = vec![self.add(&f[0], c)?];
            h.extend_from_slice(&f[1..]);
            return Ok(self.strip(h));
        }
        if i >= n {
            let mut h = vec![c.clone()];
            h.extend(std::iter::repeat_n(self.zero(), i - n));
            h.extend_from_slice(f);
            return Ok(h);
        }
        let m = n - i - 1;
        let mut h = f.to_vec();
        h[m] = self.add(&f[m], c)?;
        Ok(h)
    }
    fn mul_term(&self, f: &[E], c: &E, i: usize) -> PResult<U> {
        if self.is_zero(c) || f.is_empty() {
            return Ok(Vec::new());
        }
        let mut h: U = f.iter().map(|x| self.mul(x, c)).collect::<PResult<_>>()?;
        h.extend(std::iter::repeat_n(self.zero(), i));
        Ok(h)
    }
    fn mul_ground(&self, f: &[E], c: &E) -> PResult<U> {
        if self.is_zero(c) || f.is_empty() {
            return Ok(Vec::new());
        }
        f.iter().map(|x| self.mul(x, c)).collect()
    }
    fn quo_ground(&self, f: &[E], c: &E) -> PResult<U> {
        if self.is_zero(c) {
            return Err(zero_division());
        }
        f.iter().map(|x| self.quo(x, c)).collect()
    }
    fn exquo_ground(&self, f: &[E], c: &E) -> PResult<U> {
        if self.is_zero(c) {
            return Err(zero_division());
        }
        f.iter().map(|x| self.exquo(x, c)).collect()
    }
    fn d_monic(&self, f: &[E]) -> PResult<U> {
        if f.is_empty() {
            return Ok(Vec::new());
        }
        let lc = &f[0];
        if self.is_one(lc) {
            return Ok(f.to_vec());
        }
        self.exquo_ground(f, lc)
    }
    fn d_content(&self, f: &[E]) -> PResult<E> {
        let mut cont = self.zero();
        for c in f {
            cont = self.gcd(&cont, c)?;
            if self.is_one(&cont) {
                break;
            }
        }
        Ok(cont)
    }
    fn d_primitive(&self, f: &[E]) -> PResult<(E, U)> {
        if f.is_empty() {
            return Ok((self.zero(), Vec::new()));
        }
        let cont = self.d_content(f)?;
        if self.is_one(&cont) {
            return Ok((cont, f.to_vec()));
        }
        let q = self.quo_ground(f, &cont)?;
        Ok((cont, q))
    }
    fn division_failed(&self, f: &[E], g: &[E]) -> PolyError {
        reference("PolynomialDivisionFailed", format!("{:?}\u{1f}{:?}", f.len(), g.len()))
    }
    fn rr_div(&self, f: &[E], g: &[E]) -> PResult<(U, U)> {
        let (df, dg) = (degree(f), degree(g));
        if g.is_empty() {
            return Err(zero_division());
        }
        let (mut q, mut r, mut dr) = (Vec::new(), f.to_vec(), df);
        if df < dg {
            return Ok((q, r));
        }
        let lc_g = g[0].clone();
        loop {
            let lc_r = self.lc_of(&r);
            if !self.is_zero(&self.rem(&lc_r, &lc_g)?) {
                break;
            }
            let c = self.exquo(&lc_r, &lc_g)?;
            let j = (dr - dg) as usize;
            q = self.add_term(&q, &c, j)?;
            let h = self.mul_term(g, &c, j)?;
            r = self.d_sub(&r, &h)?;
            let pdr = dr;
            dr = degree(&r);
            if dr < dg {
                break;
            } else if dr >= pdr {
                return Err(self.division_failed(f, g));
            }
        }
        Ok((q, r))
    }
    fn ff_div(&self, f: &[E], g: &[E]) -> PResult<(U, U)> {
        let (df, dg) = (degree(f), degree(g));
        if g.is_empty() {
            return Err(zero_division());
        }
        let (mut q, mut r, mut dr) = (Vec::new(), f.to_vec(), df);
        if df < dg {
            return Ok((q, r));
        }
        let lc_g = g[0].clone();
        loop {
            let lc_r = self.lc_of(&r);
            let c = self.exquo(&lc_r, &lc_g)?;
            let j = (dr - dg) as usize;
            q = self.add_term(&q, &c, j)?;
            let h = self.mul_term(g, &c, j)?;
            r = self.d_sub(&r, &h)?;
            let pdr = dr;
            dr = degree(&r);
            if dr < dg {
                break;
            } else if dr == pdr {
                // remove the leading term created by rounding error
                r = self.strip(r[1..].to_vec());
                dr = degree(&r);
                if dr < dg {
                    break;
                }
            } else if dr > pdr {
                return Err(self.division_failed(f, g));
            }
        }
        Ok((q, r))
    }
    fn d_div(&self, f: &[E], g: &[E]) -> PResult<(U, U)> {
        if self.field {
            self.ff_div(f, g)
        } else {
            self.rr_div(f, g)
        }
    }
    fn d_quo(&self, f: &[E], g: &[E]) -> PResult<U> {
        Ok(self.d_div(f, g)?.0)
    }
    fn d_exquo(&self, f: &[E], g: &[E]) -> PResult<U> {
        let (q, r) = self.d_div(f, g)?;
        if r.is_empty() {
            Ok(q)
        } else {
            Err(PolyError::ExactQuotientFailed(format!("{:?} does not divide {:?}", g.len(), f.len())))
        }
    }
    fn d_pdiv(&self, f: &[E], g: &[E]) -> PResult<(U, U)> {
        let (df, dg) = (degree(f), degree(g));
        if g.is_empty() {
            return Err(zero_division());
        }
        let (mut q, mut r, mut dr) = (Vec::new(), f.to_vec(), df);
        if df < dg {
            return Ok((q, r));
        }
        let mut n = df - dg + 1;
        let lc_g = g[0].clone();
        loop {
            let lc_r = self.lc_of(&r);
            let j = (dr - dg) as usize;
            n -= 1;
            let qq = self.mul_ground(&q, &lc_g)?;
            q = self.add_term(&qq, &lc_r, j)?;
            let rr = self.mul_ground(&r, &lc_g)?;
            let gg = self.mul_term(g, &lc_r, j)?;
            r = self.d_sub(&rr, &gg)?;
            let pdr = dr;
            dr = degree(&r);
            if dr < dg {
                break;
            } else if dr >= pdr {
                return Err(self.division_failed(f, g));
            }
        }
        let c = self.pow(&lc_g, n)?;
        Ok((self.mul_ground(&q, &c)?, self.mul_ground(&r, &c)?))
    }
    fn d_prem(&self, f: &[E], g: &[E]) -> PResult<U> {
        let (df, dg) = (degree(f), degree(g));
        if g.is_empty() {
            return Err(zero_division());
        }
        let (mut r, mut dr) = (f.to_vec(), df);
        if df < dg {
            return Ok(r);
        }
        let mut n = df - dg + 1;
        let lc_g = g[0].clone();
        loop {
            let lc_r = self.lc_of(&r);
            let j = (dr - dg) as usize;
            n -= 1;
            let rr = self.mul_ground(&r, &lc_g)?;
            let gg = self.mul_term(g, &lc_r, j)?;
            r = self.d_sub(&rr, &gg)?;
            let pdr = dr;
            dr = degree(&r);
            if dr < dg {
                break;
            } else if dr >= pdr {
                return Err(self.division_failed(f, g));
            }
        }
        self.mul_ground(&r, &self.pow(&lc_g, n)?)
    }
    fn inner_subresultants(&self, f: &[E], g: &[E]) -> PResult<(Vec<U>, Vec<E>)> {
        let (mut f, mut g) = (f.to_vec(), g.to_vec());
        let (mut n, mut m) = (degree(&f), degree(&g));
        if n < m {
            std::mem::swap(&mut f, &mut g);
            std::mem::swap(&mut n, &mut m);
        }
        if f.is_empty() {
            return Ok((Vec::new(), Vec::new()));
        }
        if g.is_empty() {
            return Ok((vec![f], vec![self.one()]));
        }
        let mut rs = vec![f.clone(), g.clone()];
        let mut d = n - m;
        let b = self.pow(&self.neg(&self.one()), d + 1)?;
        let mut h = self.d_prem(&f, &g)?;
        h = self.mul_ground(&h, &b)?;
        let mut lc = self.lc_of(&g);
        let mut c = self.pow(&lc, d)?;
        let mut s = vec![self.one(), c.clone()];
        c = self.neg(&c);
        while !h.is_empty() {
            let k = degree(&h);
            rs.push(h.clone());
            f = std::mem::replace(&mut g, h);
            d = m - k;
            m = k;
            let b = self.mul(&self.neg(&lc), &self.pow(&c, d)?)?;
            h = self.d_prem(&f, &g)?;
            h = self.quo_ground(&h, &b)?;
            lc = self.lc_of(&g);
            if d > 1 {
                let q = self.pow(&c, d - 1)?;
                c = self.quo(&self.pow(&self.neg(&lc), d)?, &q)?;
            } else {
                c = self.neg(&lc);
            }
            s.push(self.neg(&c));
        }
        Ok((rs, s))
    }
    fn prs_resultant(&self, f: &[E], g: &[E]) -> PResult<(E, Vec<U>)> {
        if f.is_empty() || g.is_empty() {
            return Ok((self.zero(), Vec::new()));
        }
        let (rs, s) = self.inner_subresultants(f, g)?;
        if degree(rs.last().unwrap()) > 0 {
            return Ok((self.zero(), rs));
        }
        Ok((s.last().unwrap().clone(), rs))
    }
    fn d_diff(&self, f: &[E], m: usize) -> PResult<U> {
        if m == 0 {
            return Ok(f.to_vec());
        }
        let mut n = degree(f);
        if n < m as i64 {
            return Ok(Vec::new());
        }
        let mut deriv = Vec::new();
        for coeff in &f[..f.len() - m] {
            let mut k = n;
            for i in (n - m as i64 + 1..n).rev() {
                k *= i;
            }
            deriv.push(self.mul(&self.from_int(k), coeff)?);
            n -= 1;
        }
        Ok(self.strip(deriv))
    }
    fn d_integrate(&self, f: &[E], m: usize) -> PResult<U> {
        if m == 0 || f.is_empty() {
            return Ok(f.to_vec());
        }
        let mut g: U = vec![self.zero(); m];
        for (i, c) in f.iter().rev().enumerate() {
            let mut n = (i + 1) as i64;
            for j in 1..m {
                n *= (i + j + 1) as i64;
            }
            g.insert(0, self.exquo(c, &self.from_int(n))?);
        }
        Ok(g)
    }
    fn d_eval(&self, f: &[E], a: &E) -> PResult<E> {
        if self.is_zero(a) {
            return Ok(self.tc_of(f));
        }
        let mut result = self.zero();
        for c in f {
            result = self.mul(&result, a)?;
            result = self.add(&result, c)?;
        }
        Ok(result)
    }
    fn d_compose(&self, f: &[E], g: &[E]) -> PResult<U> {
        if g.len() <= 1 {
            return Ok(self.strip(vec![self.d_eval(f, &self.lc_of(g))?]));
        }
        if f.is_empty() {
            return Ok(Vec::new());
        }
        let mut h = vec![f[0].clone()];
        for c in &f[1..] {
            h = self.d_mul(&h, g)?;
            h = self.add_term(&h, c, 0)?;
        }
        Ok(h)
    }
    fn trivial_gcd(&self, f: &[E], g: &[E]) -> PResult<Option<(U, U, U)>> {
        if f.is_empty() && g.is_empty() {
            return Ok(Some((Vec::new(), Vec::new(), Vec::new())));
        }
        if self.field {
            if f.is_empty() {
                return Ok(Some((self.d_monic(g)?, Vec::new(), vec![self.lc_of(g)])));
            }
            if g.is_empty() {
                return Ok(Some((self.d_monic(f)?, vec![self.lc_of(f)], Vec::new())));
            }
            return Ok(None);
        }
        if f.is_empty() {
            return Ok(Some(if self.is_nonnegative(&self.lc_of(g)) {
                (g.to_vec(), Vec::new(), vec![self.one()])
            } else {
                (self.d_neg(g), Vec::new(), vec![self.neg(&self.one())])
            }));
        }
        if g.is_empty() {
            return Ok(Some(if self.is_nonnegative(&self.lc_of(f)) {
                (f.to_vec(), vec![self.one()], Vec::new())
            } else {
                (self.d_neg(f), vec![self.neg(&self.one())], Vec::new())
            }));
        }
        Ok(None)
    }
    /// `dup_inner_gcd` over the composite domain (`dup_rr_prs_gcd` or
    /// `dup_ff_prs_gcd`).
    fn inner_gcd(&self, f: &[E], g: &[E]) -> PResult<(U, U, U)> {
        if let Some(r) = self.trivial_gcd(f, g)? {
            return Ok(r);
        }
        if self.field {
            let (rs, _) = self.inner_subresultants(f, g)?;
            let h = self.d_monic(rs.last().unwrap())?;
            let cff = self.d_quo(f, &h)?;
            let cfg = self.d_quo(g, &h)?;
            return Ok((h, cff, cfg));
        }
        let (fc, ff) = self.d_primitive(f)?;
        let (gc, gg) = self.d_primitive(g)?;
        let mut c = self.gcd(&fc, &gc)?;
        let (rs, _) = self.inner_subresultants(&ff, &gg)?;
        let (_, h) = self.d_primitive(rs.last().unwrap())?;
        c = self.mul(&c, &self.canonical_unit(&self.lc_of(&h)))?;
        let h = self.mul_ground(&h, &c)?;
        let cff = self.d_quo(f, &h)?;
        let cfg = self.d_quo(g, &h)?;
        Ok((h, cff, cfg))
    }
    fn d_gcd(&self, f: &[E], g: &[E]) -> PResult<U> {
        Ok(self.inner_gcd(f, g)?.0)
    }
    fn d_lcm(&self, f: &[E], g: &[E]) -> PResult<U> {
        if self.field {
            let h = self.d_quo(&self.d_mul(f, g)?, &self.d_gcd(f, g)?)?;
            return self.d_monic(&h);
        }
        if f.is_empty() || g.is_empty() {
            return Ok(Vec::new());
        }
        let (fc, ff) = self.d_primitive(f)?;
        let (gc, gg) = self.d_primitive(g)?;
        let c = self.lcm(&fc, &gc)?;
        let h = self.d_quo(&self.d_mul(&ff, &gg)?, &self.d_gcd(&ff, &gg)?)?;
        let u = self.canonical_unit(&self.lc_of(&h));
        self.mul_ground(&h, &self.mul(&c, &u)?)
    }
    /// `dup_clear_denoms(f, K0, K1=ring, convert=True)`.
    fn clear_denoms_ring(&self, f: &[E]) -> PResult<(R, Vec<R>)> {
        let ring = self.as_ring();
        let mut common = self.r.one();
        for c in f {
            common = self.r.lcm(&common, &c.d)?;
        }
        if self.r.is_one(&common) {
            return Ok((common, f.iter().map(|c| c.n.clone()).collect()));
        }
        let out = f.iter().map(|c| Ok(self.r.mul(&c.n, &self.r.quo(&common, &c.d)?))).collect::<PResult<Vec<_>>>()?;
        let _ = ring;
        Ok((common, out))
    }
    fn half_gcdex(&self, f: &[E], g: &[E]) -> PResult<(U, U)> {
        if !self.field {
            return Err(reference("DomainError", format!("Cannot compute half extended GCD over {}", self.display())));
        }
        let (mut a, mut b): (U, U) = (vec![self.one()], Vec::new());
        let (mut f, mut g) = (f.to_vec(), g.to_vec());
        while !g.is_empty() {
            let (q, r) = self.d_div(&f, &g)?;
            f = std::mem::replace(&mut g, r);
            let nb = self.d_sub(&a, &self.d_mul(&q, &b)?)?;
            a = std::mem::replace(&mut b, nb);
        }
        let a = self.quo_ground(&a, &self.lc_of(&f))?;
        let f = self.d_monic(&f)?;
        Ok((a, f))
    }
    fn gcdex(&self, f: &[E], g: &[E]) -> PResult<(U, U, U)> {
        let (s, h) = self.half_gcdex(f, g)?;
        let ff = self.d_sub(&h, &self.d_mul(&s, f)?)?;
        let t = self.d_quo(&ff, g)?;
        Ok((s, t, h))
    }
    fn invert(&self, f: &[E], g: &[E]) -> PResult<U> {
        let (s, h) = self.half_gcdex(f, g)?;
        if h.len() == 1 && self.is_one(&h[0]) {
            Ok(self.d_div(&s, g)?.1)
        } else {
            Err(PolyError::NotInvertible("zero divisor".into()))
        }
    }
    fn sqf_list(&self, f: &[E], all: bool) -> PResult<(E, Vec<(U, usize)>)> {
        let (coeff, f) = if self.field {
            (self.lc_of(f), self.d_monic(f)?)
        } else {
            let (mut coeff, mut f) = self.d_primitive(f)?;
            if self.is_negative(&self.lc_of(&f)) {
                f = self.d_neg(&f);
                coeff = self.neg(&coeff);
            }
            (coeff, f)
        };
        if degree(&f) <= 0 {
            return Ok((coeff, Vec::new()));
        }
        let mut result = Vec::new();
        let mut i = 1;
        let h = self.d_diff(&f, 1)?;
        let (_, mut p, mut q) = self.inner_gcd(&f, &h)?;
        loop {
            let d = self.d_diff(&p, 1)?;
            let h = self.d_sub(&q, &d)?;
            if h.is_empty() {
                result.push((p, i));
                break;
            }
            let (g, np, nq) = self.inner_gcd(&p, &h)?;
            p = np;
            q = nq;
            if all || degree(&g) > 0 {
                result.push((g, i));
            }
            i += 1;
        }
        let deg: i64 = result.iter().map(|(p, k)| *k as i64 * degree(p)).sum();
        if deg != degree(&f) {
            return Err(reference("AssertionError", ""));
        }
        Ok((coeff, result))
    }
    fn sqf_part(&self, f: &[E]) -> PResult<U> {
        if f.is_empty() {
            return Ok(Vec::new());
        }
        let f = if self.is_negative(&self.lc_of(f)) { self.d_neg(f) } else { f.to_vec() };
        let gcd = self.d_gcd(&f, &self.d_diff(&f, 1)?)?;
        let sqf = self.d_quo(&f, &gcd)?;
        if self.field {
            self.d_monic(&sqf)
        } else {
            Ok(self.d_primitive(&sqf)?.1)
        }
    }

    // ---- conversion to and from EX polynomials ----

    fn dense(&self, p: &ExpressionPoly) -> PResult<U> {
        let terms = p.terms();
        let Some(top) = terms.first().map(|(m, _)| m[0] as usize) else { return Ok(Vec::new()) };
        let mut f = vec![self.zero(); top + 1];
        for (m, c) in terms {
            f[top - m[0] as usize] = self.from_expr(&c)?;
        }
        Ok(self.strip(f))
    }
    fn poly(&self, f: &[E], gens: &[Expr]) -> PResult<ExpressionPoly> {
        let n = f.len();
        let terms = f
            .iter()
            .enumerate()
            .filter(|(_, c)| !self.is_zero(c))
            .map(|(i, c)| (vec![(n - 1 - i) as u32], self.to_expr(c)))
            .collect();
        let mut p = ExpressionPoly::from_terms(terms, gens)?;
        p.ground = Some(self.ground.clone());
        Ok(p)
    }
}

// ----------------------------------------------------------------------
// public entry points
// ----------------------------------------------------------------------

fn univariate(p: &ExpressionPoly) -> PResult<()> {
    if p.gens.len() != 1 {
        return Err(PolyError::Unsupported(
            "multivariate polynomials over inexact composite domains are only supported for ring arithmetic".into(),
        ));
    }
    Ok(())
}

fn dom(p: &ExpressionPoly) -> Dom {
    Dom::new(p.ground.as_ref().expect("inexact ground")).expect("inexact ground")
}

/// Whether `p` carries an inexact composite ground.
pub fn is_inexact(p: &ExpressionPoly) -> bool {
    p.ground.as_ref().is_some_and(|g| g.inexact.is_some())
}

/// Convert every coefficient of `p` into the ground of `g`.
pub fn convert_poly(g: &GaussianGround, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let d = Dom::new(g).expect("inexact ground");
    let mut terms = Vec::new();
    for (m, c) in p.terms() {
        let e = d.from_expr(&c)?;
        if !d.is_zero(&e) {
            terms.push((m, d.to_expr(&e)));
        }
    }
    let mut out = ExpressionPoly::from_terms(terms, &p.gens)?;
    out.ground = Some(g.clone());
    Ok(out)
}

/// The ground element of `g` for `e`, as its normal form.
pub fn element_expr(g: &GaussianGround, e: &Expr) -> PResult<Expr> {
    let d = Dom::new(g).expect("inexact ground");
    Ok(d.to_expr(&d.from_expr(e)?))
}

/// The numerator and denominator of an element of `RR(y...)`, `CC(y...)`
/// (the reference `FracElement.numer` and `.denom`) as ring expressions.
pub fn element_parts(g: &GaussianGround, e: &Expr) -> PResult<(Expr, Expr)> {
    let d = Dom::new(g).expect("inexact ground");
    let a = d.from_expr(e)?;
    Ok((d.r.to_expr(&a.n, &d.symbols), d.r.to_expr(&a.d, &d.symbols)))
}

/// Domain element operations of `RR[y]`, `CC(y)`, ...: `op` is one of
/// `+ - * / quo rem exquo gcd lcm`; `/` is the exact quotient in a ring.
pub fn element_op(g: &GaussianGround, a: &Expr, b: &Expr, op: &str) -> PResult<Expr> {
    let d = Dom::new(g).expect("inexact ground");
    let (x, y) = (d.from_expr(a)?, d.from_expr(b)?);
    let r = match op {
        "+" => d.add(&x, &y)?,
        "-" => d.sub(&x, &y)?,
        "*" => d.mul(&x, &y)?,
        "/" | "exquo" => d.exquo(&x, &y)?,
        "quo" => d.quo(&x, &y)?,
        "rem" => d.rem(&x, &y)?,
        "gcd" => d.gcd(&x, &y)?,
        "lcm" => d.lcm(&x, &y)?,
        _ => return Err(PolyError::Unsupported(format!("unknown coefficient operation {op}"))),
    };
    Ok(d.to_expr(&r))
}

/// Sparse arithmetic in several generators: `dmp_*` over the composite
/// domain, through the recursive dense representation.
mod multi {
    use super::*;

    #[derive(Clone, Debug, PartialEq)]
    pub enum D {
        L(Vec<E>),
        N(Vec<D>),
    }

    pub fn is_zero(f: &D) -> bool {
        match f {
            D::L(v) => v.is_empty(),
            D::N(v) => v.is_empty(),
        }
    }
    fn zero_like(f: &D) -> D {
        match f {
            D::L(_) => D::L(Vec::new()),
            D::N(_) => D::N(Vec::new()),
        }
    }
    fn strip(v: Vec<D>) -> Vec<D> {
        let k = v.iter().take_while(|c| is_zero(c)).count();
        v.into_iter().skip(k).collect()
    }
    pub fn add(k: &Dom, f: &D, g: &D) -> PResult<D> {
        Ok(match (f, g) {
            (D::L(a), D::L(b)) => D::L(k.d_add(a, b)?),
            (D::N(a), D::N(b)) => {
                if a.is_empty() {
                    return Ok(g.clone());
                }
                if b.is_empty() {
                    return Ok(f.clone());
                }
                let (da, db) = (a.len(), b.len());
                if da == db {
                    D::N(strip(a.iter().zip(b).map(|(x, y)| add(k, x, y)).collect::<PResult<_>>()?))
                } else {
                    let kk = da.abs_diff(db);
                    let (mut h, x, y) = if da > db { (a[..kk].to_vec(), &a[kk..], &b[..]) } else { (b[..kk].to_vec(), &a[..], &b[kk..]) };
                    for (p, q) in x.iter().zip(y) {
                        h.push(add(k, p, q)?);
                    }
                    D::N(h)
                }
            }
            _ => unreachable!(),
        })
    }
    pub fn neg(k: &Dom, f: &D) -> D {
        match f {
            D::L(a) => D::L(k.d_neg(a)),
            D::N(a) => D::N(a.iter().map(|x| neg(k, x)).collect()),
        }
    }
    pub fn sub(k: &Dom, f: &D, g: &D) -> PResult<D> {
        Ok(match (f, g) {
            (D::L(a), D::L(b)) => D::L(k.d_sub(a, b)?),
            (D::N(a), D::N(b)) => {
                if a.is_empty() {
                    return Ok(neg(k, g));
                }
                if b.is_empty() {
                    return Ok(f.clone());
                }
                let (da, db) = (a.len(), b.len());
                if da == db {
                    D::N(strip(a.iter().zip(b).map(|(x, y)| sub(k, x, y)).collect::<PResult<_>>()?))
                } else {
                    let kk = da.abs_diff(db);
                    let (mut h, x, y) =
                        if da > db { (a[..kk].to_vec(), &a[kk..], &b[..]) } else { (b[..kk].iter().map(|x| neg(k, x)).collect(), &a[..], &b[kk..]) };
                    for (p, q) in x.iter().zip(y) {
                        h.push(sub(k, p, q)?);
                    }
                    D::N(h)
                }
            }
            _ => unreachable!(),
        })
    }
    pub fn mul(k: &Dom, f: &D, g: &D) -> PResult<D> {
        match (f, g) {
            (D::L(a), D::L(b)) => Ok(D::L(k.d_mul(a, b)?)),
            (D::N(a), D::N(b)) => {
                if f == g {
                    return sqr(k, f);
                }
                if a.is_empty() || b.is_empty() {
                    return Ok(zero_like(f));
                }
                let (df, dg) = (a.len() - 1, b.len() - 1);
                let z = zero_like(&a[0]);
                let mut h = Vec::new();
                for i in 0..=df + dg {
                    let mut coeff = z.clone();
                    for j in i.saturating_sub(dg)..=df.min(i) {
                        coeff = add(k, &coeff, &mul(k, &a[j], &b[i - j])?)?;
                    }
                    h.push(coeff);
                }
                Ok(D::N(strip(h)))
            }
            _ => unreachable!(),
        }
    }
    pub fn sqr(k: &Dom, f: &D) -> PResult<D> {
        match f {
            D::L(a) => Ok(D::L(k.d_sqr(a)?)),
            D::N(a) => {
                if a.is_empty() {
                    return Ok(f.clone());
                }
                let df = a.len() - 1;
                let z = zero_like(&a[0]);
                let mut h = Vec::new();
                for i in 0..=2 * df {
                    let mut c = z.clone();
                    let jmin = i.saturating_sub(df);
                    let jmax0 = i.min(df);
                    let n = jmax0 - jmin + 1;
                    let jmax = (jmin + n / 2) as i64 - 1;
                    let mut j = jmin as i64;
                    while j <= jmax {
                        c = add(k, &c, &mul(k, &a[j as usize], &a[i - j as usize])?)?;
                        j += 1;
                    }
                    c = mul_ground(k, &c, &k.from_int(2))?;
                    if n & 1 == 1 {
                        c = add(k, &c, &sqr(k, &a[(jmax + 1) as usize])?)?;
                    }
                    h.push(c);
                }
                Ok(D::N(strip(h)))
            }
        }
    }
    fn mul_ground(k: &Dom, f: &D, c: &E) -> PResult<D> {
        Ok(match f {
            D::L(a) => D::L(k.mul_ground(a, c)?),
            D::N(a) => D::N(a.iter().map(|x| mul_ground(k, x, c)).collect::<PResult<_>>()?),
        })
    }
    pub fn pow(k: &Dom, f: &D, n: u32) -> PResult<D> {
        let one = |f: &D| -> D {
            fn build(f: &D, k: &Dom) -> D {
                match f {
                    D::L(_) => D::L(vec![k.one()]),
                    D::N(v) => D::N(vec![build(v.first().unwrap_or(&D::L(Vec::new())), k)]),
                }
            }
            build(f, k)
        };
        if n == 0 {
            return Ok(one_of_depth(k, f));
        }
        let _ = one;
        if n == 1 || is_zero(f) || *f == one_of_depth(k, f) {
            return Ok(f.clone());
        }
        let mut g = one_of_depth(k, f);
        let mut f = f.clone();
        let mut n = n;
        loop {
            let m = n;
            n /= 2;
            if m % 2 == 1 {
                g = mul(k, &g, &f)?;
                if n == 0 {
                    break;
                }
            }
            f = sqr(k, &f)?;
        }
        Ok(g)
    }
    fn depth(f: &D) -> usize {
        match f {
            D::L(_) => 0,
            D::N(v) => 1 + v.first().map(depth).unwrap_or(0),
        }
    }
    fn one_at(k: &Dom, u: usize) -> D {
        if u == 0 {
            D::L(vec![k.one()])
        } else {
            D::N(vec![one_at(k, u - 1)])
        }
    }
    fn one_of_depth(k: &Dom, f: &D) -> D {
        one_at(k, depth(f))
    }
    pub fn from_poly(k: &Dom, p: &ExpressionPoly) -> PResult<D> {
        let u = p.gens.len() - 1;
        let terms: Vec<(Vec<u32>, E)> = p.terms().into_iter().map(|(m, c)| Ok((m, k.from_expr(&c)?))).collect::<PResult<_>>()?;
        Ok(build(k, &terms, u, 0))
    }
    fn build(k: &Dom, terms: &[(Vec<u32>, E)], u: usize, level: usize) -> D {
        let top = terms.iter().map(|(m, _)| m[level]).max();
        let Some(top) = top else {
            return if level == u { D::L(Vec::new()) } else { D::N(Vec::new()) };
        };
        if level == u {
            let mut v = vec![k.zero(); top as usize + 1];
            for (m, c) in terms {
                v[(top - m[level]) as usize] = c.clone();
            }
            return D::L(v);
        }
        let mut out = Vec::new();
        for d in (0..=top).rev() {
            let sub: Vec<(Vec<u32>, E)> = terms.iter().filter(|(m, _)| m[level] == d).cloned().collect();
            out.push(build(k, &sub, u, level + 1));
        }
        D::N(strip(out))
    }
    pub fn to_terms(f: &D, prefix: &mut Vec<u32>, out: &mut Vec<(Vec<u32>, E)>) {
        match f {
            D::L(v) => {
                let n = v.len();
                for (i, c) in v.iter().enumerate() {
                    if !c.n.is_empty() {
                        let mut m = prefix.clone();
                        m.push((n - 1 - i) as u32);
                        out.push((m, c.clone()));
                    }
                }
            }
            D::N(v) => {
                let n = v.len();
                for (i, c) in v.iter().enumerate() {
                    prefix.push((n - 1 - i) as u32);
                    to_terms(c, prefix, out);
                    prefix.pop();
                }
            }
        }
    }
    pub fn to_poly(k: &Dom, f: &D, gens: &[Expr]) -> PResult<ExpressionPoly> {
        let mut out = Vec::new();
        to_terms(f, &mut Vec::new(), &mut out);
        let terms = out.into_iter().map(|(m, c)| (m, k.to_expr(&c))).collect();
        let mut p = ExpressionPoly::from_terms(terms, gens)?;
        p.ground = Some(k.ground.clone());
        Ok(p)
    }
}

fn binary(
    f: &ExpressionPoly,
    g: &ExpressionPoly,
    uni: impl Fn(&Dom, &[E], &[E]) -> PResult<U>,
    many: impl Fn(&Dom, &multi::D, &multi::D) -> PResult<multi::D>,
) -> PResult<ExpressionPoly> {
    let k = dom(f);
    if f.gens.len() == 1 {
        return k.poly(&uni(&k, &k.dense(f)?, &k.dense(g)?)?, &f.gens);
    }
    multi::to_poly(&k, &many(&k, &multi::from_poly(&k, f)?, &multi::from_poly(&k, g)?)?, &f.gens)
}

pub fn add(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    binary(f, g, |k, a, b| k.d_add(a, b), |k, a, b| multi::add(k, a, b))
}
pub fn sub(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    binary(f, g, |k, a, b| k.d_sub(a, b), |k, a, b| multi::sub(k, a, b))
}
pub fn mul(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    binary(f, g, |k, a, b| k.d_mul(a, b), |k, a, b| multi::mul(k, a, b))
}
pub fn neg(f: &ExpressionPoly) -> PResult<ExpressionPoly> {
    binary(f, f, |k, a, _| Ok(k.d_neg(a)), |k, a, _| Ok(multi::neg(k, a)))
}
pub fn pow(f: &ExpressionPoly, n: u32) -> PResult<ExpressionPoly> {
    binary(f, f, |k, a, _| k.d_pow(a, n), |k, a, _| multi::pow(k, a, n))
}

fn pair<T>(f: &ExpressionPoly, g: &ExpressionPoly, op: impl Fn(&Dom, &[E], &[E]) -> PResult<T>) -> PResult<(Dom, T)> {
    univariate(f)?;
    let k = dom(f);
    let r = op(&k, &k.dense(f)?, &k.dense(g)?)?;
    Ok((k, r))
}

/// `Poly.div` with `auto`: over a ring the division runs in the fraction
/// field and the results are retracted when every coefficient lies in the
/// ring.
pub fn div(f: &ExpressionPoly, g: &ExpressionPoly, auto: bool) -> PResult<(ExpressionPoly, ExpressionPoly)> {
    univariate(f)?;
    let k = dom(f);
    let (a, b) = (k.dense(f)?, k.dense(g)?);
    if k.field || !auto {
        let (q, r) = k.d_div(&a, &b)?;
        return Ok((k.poly(&q, &f.gens)?, k.poly(&r, &f.gens)?));
    }
    let kf = k.as_field();
    let lift = |v: &[E]| -> U { v.to_vec() };
    let (q, r) = kf.d_div(&lift(&a), &lift(&b))?;
    let retract = |v: &[E]| -> Option<U> { v.iter().map(|c| kf.retract(c)).collect() };
    match (retract(&q), retract(&r)) {
        (Some(q), Some(r)) => Ok((k.poly(&q, &f.gens)?, k.poly(&r, &f.gens)?)),
        _ => Ok((kf.poly(&q, &f.gens)?, kf.poly(&r, &f.gens)?)),
    }
}
/// `Poly.rem`/`Poly.quo` with `auto`: each result is retracted on its own.
pub fn div_part(f: &ExpressionPoly, g: &ExpressionPoly, auto: bool, quotient: bool) -> PResult<ExpressionPoly> {
    univariate(f)?;
    let k = dom(f);
    let (a, b) = (k.dense(f)?, k.dense(g)?);
    let pick = |(q, r): (U, U)| if quotient { q } else { r };
    if k.field || !auto {
        return k.poly(&pick(k.d_div(&a, &b)?), &f.gens);
    }
    let kf = k.as_field();
    let v = pick(kf.d_div(&a, &b)?);
    match v.iter().map(|c| kf.retract(c)).collect::<Option<U>>() {
        Some(v) => k.poly(&v, &f.gens),
        None => kf.poly(&v, &f.gens),
    }
}
pub fn exquo(f: &ExpressionPoly, g: &ExpressionPoly, auto: bool) -> PResult<ExpressionPoly> {
    univariate(f)?;
    let k = dom(f);
    let (a, b) = (k.dense(f)?, k.dense(g)?);
    if k.field || !auto {
        return k.poly(&k.d_exquo(&a, &b)?, &f.gens);
    }
    let kf = k.as_field();
    let v = kf.d_exquo(&a, &b)?;
    match v.iter().map(|c| kf.retract(c)).collect::<Option<U>>() {
        Some(v) => k.poly(&v, &f.gens),
        None => kf.poly(&v, &f.gens),
    }
}
pub fn pdiv(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(ExpressionPoly, ExpressionPoly)> {
    let (k, (q, r)) = pair(f, g, |k, a, b| k.d_pdiv(a, b))?;
    Ok((k.poly(&q, &f.gens)?, k.poly(&r, &f.gens)?))
}
pub fn prem(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, r) = pair(f, g, |k, a, b| k.d_prem(a, b))?;
    k.poly(&r, &f.gens)
}
pub fn pexquo(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, (q, r)) = pair(f, g, |k, a, b| k.d_pdiv(a, b))?;
    if !r.is_empty() {
        return Err(PolyError::ExactQuotientFailed(format!("{} does not divide {}", g.as_expr(), f.as_expr())));
    }
    k.poly(&q, &f.gens)
}
pub fn gcd(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, h) = pair(f, g, |k, a, b| k.d_gcd(a, b))?;
    k.poly(&h, &f.gens)
}
pub fn lcm(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, h) = pair(f, g, |k, a, b| k.d_lcm(a, b))?;
    k.poly(&h, &f.gens)
}
pub fn cofactors(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(ExpressionPoly, ExpressionPoly, ExpressionPoly)> {
    let (k, (h, a, b)) = pair(f, g, |k, a, b| k.inner_gcd(a, b))?;
    Ok((k.poly(&h, &f.gens)?, k.poly(&a, &f.gens)?, k.poly(&b, &f.gens)?))
}
/// `dmp_cancel(f, g, include=True)` over `RR` (53 bits) for the expression
/// `cancel` of float input, whose coefficients are given as the exact values
/// of their doubles: `RR` has no associated ring, so no denominators are
/// cleared; the gcd is computed over `QQ` after `to_rational(limit=True)`
/// (`_dmp_inner_gcd` of an inexact domain), the cofactors of the monic gcd
/// come back through `RR.from_QQ` and only the signs of the leading
/// coefficients are normalized. `x/(0.3*x + 0.1)` keeps its coefficients,
/// where clearing the binary denominators gave `2**55`-sized integers.
fn rr_inner_gcd(f: &MPoly<Rat>, g: &MPoly<Rat>, precision: u32) -> (MPoly<Rat>, MPoly<Rat>) {
    let k = Inexact { precision, complex: false };
    let n = f.nvars;
    let to_qq = |p: &MPoly<Rat>| -> MPoly<Rat> {
        MPoly::from_terms(
            n,
            p.iter()
                .map(|(e, c)| {
                    let x = BF::from_ratio(&c.n.to_big(), &c.d.to_big(), precision);
                    (e.to_vec(), k.to_rational(&x, precision))
                })
                .filter(|(_, c)| !c.is_zero())
                .collect(),
        )
    };
    let back = |p: &MPoly<Rat>| -> MPoly<Rat> {
        MPoly::from_terms(
            n,
            p.iter()
                .map(|(e, c)| {
                    let r = k.from_rat(c).re.to_ratio_exact();
                    (e.to_vec(), Rat::new(Int::from_big(r.numer().clone()), Int::from_big(r.denom().clone())))
                })
                .filter(|(_, c)| !c.is_zero())
                .collect(),
        )
    };
    let (_, p, q) = super::gcd::gcd_mq(&to_qq(f), &to_qq(g));
    (back(&p), back(&q))
}

/// The cofactors of `dmp_inner_gcd` over `RR`: the reference computes them
/// over `QQ` after `RR.to_rational` of every coefficient, and rounds back.
fn rr_cofactors(f: &MPoly<Rat>, g: &MPoly<Rat>) -> (MPoly<Rat>, MPoly<Rat>) {
    rr_inner_gcd(f, g, 53)
}

/// `Poly.cancel` over `RR` (`dmp_cancel` over a field without associated
/// ring): `(negated, p, q)`, with `cp = -1` when `negated`, `cq = 1`.
pub fn rr_poly_cancel(f: &MPoly<Rat>, g: &MPoly<Rat>) -> (bool, MPoly<Rat>, MPoly<Rat>) {
    let (p, q) = rr_cofactors(f, g);
    let (pn, qn) = (p.lc().is_negative(), q.lc().is_negative());
    if pn && qn {
        (false, p.neg(), q.neg())
    } else if pn {
        (true, p.neg(), q)
    } else if qn {
        (true, p, q.neg())
    } else {
        (false, p, q)
    }
}

pub fn rr_cancel(f: &MPoly<Rat>, g: &MPoly<Rat>) -> (MPoly<Rat>, MPoly<Rat>) {
    rr_cancel_prec(f, g, 53)
}

/// `rr_cancel` over `RealField(precision)`.
pub fn rr_cancel_prec(f: &MPoly<Rat>, g: &MPoly<Rat>, precision: u32) -> (MPoly<Rat>, MPoly<Rat>) {
    let (mut p, mut q) = rr_inner_gcd(f, g, precision);
    // `p_neg`/`q_neg`: a negative denominator is negated with the
    // numerator; a negative numerator alone is negated and multiplied by
    // `cp = -1` again (include=True), so it is kept.
    if q.lc().is_negative() {
        q = q.neg();
        p = p.neg();
    }
    (p, q)
}

/// `RealField(precision).to_rational(x)` (with its bounded denominator) of
/// the coefficient `r` rounded at `precision`.
pub fn rr_to_rational(r: &Rat, precision: u32) -> Rat {
    let k = Inexact { precision, complex: false };
    let x = BF::from_ratio(&r.n.to_big(), &r.d.to_big(), precision);
    k.to_rational(&x, precision)
}

/// `a/b` in `RealField(precision)` arithmetic, as an exact rational.
pub fn rr_quotient(a: &Rat, b: &Rat, precision: u32) -> Rat {
    let k = Inexact { precision, complex: false };
    let q = k.div(&k.from_rat(a), &k.from_rat(b)).expect("nonzero divisor").re.to_ratio_exact();
    Rat::new(Int::from_big(q.numer().clone()), Int::from_big(q.denom().clone()))
}

/// `Poly.cancel(g, include=False)`: `(cp, cq, p, q)` with `cp/cq` the
/// coefficient ratio in the associated ring.
pub fn cancel(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(Expr, Expr, ExpressionPoly, ExpressionPoly)> {
    univariate(f)?;
    let k = dom(f);
    let (a, b) = (k.dense(f)?, k.dense(g)?);
    if !k.field {
        let (mut cp, cq) = (k.one(), k.one());
        let (_, mut p, mut q) = k.inner_gcd(&a, &b)?;
        let (pn, qn) = (k.is_negative(&k.lc_of(&p)), k.is_negative(&k.lc_of(&q)));
        if pn && qn {
            p = k.d_neg(&p);
            q = k.d_neg(&q);
        } else if pn {
            cp = k.neg(&cp);
            p = k.d_neg(&p);
        } else if qn {
            cp = k.neg(&cp);
            q = k.d_neg(&q);
        }
        return Ok((k.to_expr(&cp), k.to_expr(&cq), k.poly(&p, &f.gens)?, k.poly(&q, &f.gens)?));
    }
    let kr = k.as_ring();
    let (cq0, fa) = k.clear_denoms_ring(&a)?;
    let (cp0, gb) = k.clear_denoms_ring(&b)?;
    let fa: U = fa.into_iter().map(|n| kr.raw(n)).collect();
    let gb: U = gb.into_iter().map(|n| kr.raw(n)).collect();
    let (_, p, q) = kr.inner_gcd(&fa, &gb)?;
    let (_, mut cp, cq) = kr.cofactors(&kr.raw(cp0), &kr.raw(cq0))?;
    let mut p: U = p.into_iter().map(|c| k.raw(c.n)).collect();
    let mut q: U = q.into_iter().map(|c| k.raw(c.n)).collect();
    let (pn, qn) = (k.is_negative(&k.lc_of(&p)), k.is_negative(&k.lc_of(&q)));
    if pn && qn {
        p = k.d_neg(&p);
        q = k.d_neg(&q);
    } else if pn {
        cp = kr.neg(&cp);
        p = k.d_neg(&p);
    } else if qn {
        cp = kr.neg(&cp);
        q = k.d_neg(&q);
    }
    Ok((kr.to_expr(&cp), kr.to_expr(&cq), k.poly(&p, &f.gens)?, k.poly(&q, &f.gens)?))
}
pub fn half_gcdex(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(ExpressionPoly, ExpressionPoly)> {
    let (k, (s, h)) = pair(f, g, |k, a, b| k.half_gcdex(a, b))?;
    Ok((k.poly(&s, &f.gens)?, k.poly(&h, &f.gens)?))
}
pub fn gcdex(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(ExpressionPoly, ExpressionPoly, ExpressionPoly)> {
    let (k, (s, t, h)) = pair(f, g, |k, a, b| k.gcdex(a, b))?;
    Ok((k.poly(&s, &f.gens)?, k.poly(&t, &f.gens)?, k.poly(&h, &f.gens)?))
}
pub fn invert(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, s) = pair(f, g, |k, a, b| k.invert(a, b))?;
    k.poly(&s, &f.gens)
}
/// The fraction-field version of `f` (`Poly.to_field`).
pub fn to_field(f: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let g = f.ground.as_ref().expect("inexact ground").to_field();
    convert_poly(&g, f)
}
fn unary<T>(f: &ExpressionPoly, op: impl Fn(&Dom, &[E]) -> PResult<T>) -> PResult<(Dom, T)> {
    univariate(f)?;
    let k = dom(f);
    let r = op(&k, &k.dense(f)?)?;
    Ok((k, r))
}
pub fn diff(f: &ExpressionPoly, m: usize) -> PResult<ExpressionPoly> {
    let (k, r) = unary(f, |k, a| k.d_diff(a, m))?;
    k.poly(&r, &f.gens)
}
/// `Poly.integrate` with `auto`: over a ring the integral lives in the
/// fraction field.
pub fn integrate(f: &ExpressionPoly, m: usize, auto: bool) -> PResult<ExpressionPoly> {
    univariate(f)?;
    let k = dom(f);
    let k = if auto && !k.field { k.as_field() } else { k };
    let a = k.dense(f)?;
    k.poly(&k.d_integrate(&a, m)?, &f.gens)
}
pub fn monic(f: &ExpressionPoly) -> PResult<ExpressionPoly> {
    univariate(f)?;
    let k = dom(f).as_field();
    let a = k.dense(f)?;
    k.poly(&k.d_monic(&a)?, &f.gens)
}
pub fn content(f: &ExpressionPoly) -> PResult<Expr> {
    let (k, c) = unary(f, |k, a| k.d_content(a))?;
    Ok(k.to_expr(&c))
}
pub fn primitive(f: &ExpressionPoly) -> PResult<(Expr, ExpressionPoly)> {
    let (k, (c, p)) = unary(f, |k, a| k.d_primitive(a))?;
    Ok((k.to_expr(&c), k.poly(&p, &f.gens)?))
}
pub fn mul_ground(f: &ExpressionPoly, c: &Expr) -> PResult<ExpressionPoly> {
    let k = dom(f);
    let c = k.from_expr(c)?;
    if f.gens.len() != 1 {
        let terms = f.terms().into_iter().map(|(m, x)| Ok((m, k.to_expr(&k.mul(&k.from_expr(&x)?, &c)?)))).collect::<PResult<Vec<_>>>()?;
        let mut p = ExpressionPoly::from_terms(terms.into_iter().filter(|(_, x)| !x.is_zero()).collect(), &f.gens)?;
        p.ground = f.ground.clone();
        return Ok(p);
    }
    k.poly(&k.mul_ground(&k.dense(f)?, &c)?, &f.gens)
}
pub fn quo_ground(f: &ExpressionPoly, c: &Expr) -> PResult<ExpressionPoly> {
    let (k, r) = unary(f, |k, a| {
        let c = k.from_expr(c)?;
        k.quo_ground(a, &c)
    })?;
    k.poly(&r, &f.gens)
}
pub fn exquo_ground(f: &ExpressionPoly, c: &Expr) -> PResult<ExpressionPoly> {
    let (k, r) = unary(f, |k, a| {
        let c = k.from_expr(c)?;
        k.exquo_ground(a, &c)
    })?;
    k.poly(&r, &f.gens)
}
/// `Poly.clear_denoms()`: `(coeff, f)` with `f` still over the field.
pub fn clear_denoms(f: &ExpressionPoly) -> PResult<(Expr, ExpressionPoly)> {
    univariate(f)?;
    let k = dom(f);
    let a = k.dense(f)?;
    if !k.field {
        return Ok((Expr::one(), f.clone()));
    }
    let kr = k.as_ring();
    let (common, out) = k.clear_denoms_ring(&a)?;
    let out: U = out.into_iter().map(|n| k.raw(n)).collect();
    Ok((kr.to_expr(&kr.raw(common)), k.poly(&out, &f.gens)?))
}
pub fn compose(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, r) = pair(f, g, |k, a, b| k.d_compose(a, b))?;
    k.poly(&r, &f.gens)
}
pub fn sqf_list(f: &ExpressionPoly, all: bool) -> PResult<(Expr, Vec<(ExpressionPoly, usize)>)> {
    let (k, (c, fs)) = unary(f, |k, a| k.sqf_list(a, all))?;
    let fs = fs.into_iter().map(|(p, i)| Ok((k.poly(&p, &f.gens)?, i))).collect::<PResult<Vec<_>>>()?;
    Ok((k.to_expr(&c), fs))
}
pub fn sqf_list_include(f: &ExpressionPoly, all: bool) -> PResult<Vec<(ExpressionPoly, usize)>> {
    let (k, (c, fs)) = unary(f, |k, a| k.sqf_list(a, all))?;
    let mut out = Vec::new();
    if let Some((p, 1)) = fs.first() {
        out.push((k.poly(&k.mul_ground(p, &c)?, &f.gens)?, 1));
        for (p, i) in &fs[1..] {
            out.push((k.poly(p, &f.gens)?, *i));
        }
    } else {
        out.push((k.poly(&k.strip(vec![c]), &f.gens)?, 1));
        for (p, i) in &fs {
            out.push((k.poly(p, &f.gens)?, *i));
        }
    }
    Ok(out)
}
pub fn sqf_part(f: &ExpressionPoly) -> PResult<ExpressionPoly> {
    let (k, r) = unary(f, |k, a| k.sqf_part(a))?;
    k.poly(&r, &f.gens)
}
pub fn sqf_p(f: &ExpressionPoly) -> PResult<bool> {
    let (_, r) = unary(f, |k, a| {
        if a.is_empty() {
            return Ok(true);
        }
        Ok(degree(&k.d_gcd(a, &k.d_diff(a, 1)?)?) == 0)
    })?;
    Ok(r)
}
pub fn subresultants(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<Vec<ExpressionPoly>> {
    let (k, (rs, _)) = pair(f, g, |k, a, b| k.inner_subresultants(a, b))?;
    rs.iter().map(|p| k.poly(p, &f.gens)).collect()
}
pub fn resultant(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<Expr> {
    let (k, (r, _)) = pair(f, g, |k, a, b| k.prs_resultant(a, b))?;
    Ok(k.to_expr(&r))
}
pub fn resultant_prs(f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(Expr, Vec<ExpressionPoly>)> {
    let (k, (r, rs)) = pair(f, g, |k, a, b| k.prs_resultant(a, b))?;
    Ok((k.to_expr(&r), rs.iter().map(|p| k.poly(p, &f.gens)).collect::<PResult<_>>()?))
}
pub fn discriminant(f: &ExpressionPoly) -> PResult<Expr> {
    let (k, r) = unary(f, |k, a| {
        let d = degree(a);
        if d <= 0 {
            return Ok(k.zero());
        }
        let s = if ((d * (d - 1)) / 2) % 2 == 0 { 1 } else { -1 };
        let c = k.lc_of(a);
        let (r, _) = k.prs_resultant(a, &k.d_diff(a, 1)?)?;
        k.quo(&r, &k.mul(&c, &k.from_int(s))?)
    })?;
    Ok(k.to_expr(&r))
}
/// `Poly.eval(x, a)` for a univariate polynomial.
pub fn eval(f: &ExpressionPoly, a: &Expr) -> PResult<Expr> {
    let (k, r) = unary(f, |k, p| {
        let a = k.from_expr(a)?;
        k.d_eval(p, &a)
    })?;
    Ok(k.to_expr(&r))
}
/// `dup_terms_gcd`-free `LC`/`TC` helpers are plain coefficient reads; the
/// reference raises for the norms over a composite domain.
pub fn abs_error() -> PolyError {
    reference("TypeError", "bad operand type for abs(): 'PolyElement'")
}
pub fn factor_error(g: &GaussianGround) -> PolyError {
    if g.ring {
        reference("TypeError", "unsupported operand type(s) for *: 'PolyElement' and 'PolyElement'")
    } else {
        reference("TypeError", "bad operand type for abs(): 'FracElement'")
    }
}

// ----------------------------------------------------------------------
// inference
// ----------------------------------------------------------------------

/// The reference `_construct_composite` for float coefficients: `RR`/`CC`
/// (at the largest float precision) as a polynomial ring or, with a
/// nonconstant denominator or `field`, a fraction field over the other
/// coefficient generators. `None` when the reference selects another
/// domain.
pub fn infer_coefficients(coefficients: &[Expr], field: bool) -> Option<GaussianGround> {
    let i = Expr::i();
    let (mut numers, mut denoms) = (Vec::new(), Vec::new());
    for c in coefficients {
        let (n, d) = super::api::as_numer_denom(c);
        numers.push(n);
        denoms.push(d);
    }
    let mut parts = numers.clone();
    parts.extend(denoms.iter().cloned());
    let converted = super::convert::from_exprs(&parts).ok()?;
    if !converted.float {
        return None;
    }
    let complex = converted.gens.contains(&i);
    let symbols: Vec<Expr> = converted.gens.iter().filter(|g| **g != i).cloned().collect();
    if symbols.is_empty() {
        return None;
    }
    for g in &symbols {
        if crate::traverse::free_symbols(g).is_empty() && super::gaussian_ground::is_algebraic_number(g) {
            return None;
        }
    }
    let mut seen = Vec::new();
    for g in &symbols {
        for s in crate::traverse::free_symbols(g) {
            if seen.contains(&s) {
                return None;
            }
            seen.push(s);
        }
    }
    // the largest `Float` precision among the coefficients
    let mut precision = 0u32;
    for p in &parts {
        crate::traverse::preorder(p, &mut |e: &Expr| {
            if let Some(n) = e.as_number().and_then(|n| n.float_precision()) {
                precision = precision.max(n);
            }
        });
    }
    let precision = if precision == 0 { 53 } else { precision };
    let k = coefficients.len();
    let symbol_index: Vec<usize> = converted.gens.iter().enumerate().filter(|(_, g)| **g != i).map(|(j, _)| j).collect();
    let fractions = field || converted.polys[k..].iter().any(|d| d.iter().any(|(e, _)| symbol_index.iter().any(|&j| e[j] > 0)));
    Some(GaussianGround { integral: false, ring: !fractions, symbols, inexact: Some(Inexact { precision, complex }) })
}

/// `infer_coefficients` for the coefficients of `e` in `gens`.
pub fn infer(e: &Expr, gens: &[Expr], field: bool) -> Option<GaussianGround> {
    let coefficients: Vec<Expr> = if gens.is_empty() {
        vec![e.clone()]
    } else {
        super::expression_domain::dict_from_expr(e, gens).ok()?.into_iter().map(|(_, c)| c).collect()
    };
    infer_coefficients(&coefficients, field)
}

// Univariate polynomials over the plain field `RR` (53 bits): the
// reference `dup_half_gcdex`, `dup_gcdex` and `dup_invert` in floating
// point, with `dup_ff_div` dropping the leading term created by rounding.
// Polynomials are exact rational images of their `mpf` coefficients.
type RV = Vec<C>;

fn rr_strip(k: Inexact, mut f: RV) -> RV {
    let n = f.iter().take_while(|c| k.is_zero(c)).count();
    f.drain(..n);
    f
}
fn rr_dense(k: Inexact, f: &MPoly<Rat>) -> RV {
    let n = f.iter().map(|(e, _)| e[0] as usize).max().unwrap_or(0);
    let mut v = vec![k.zero(); if f.is_zero() { 0 } else { n + 1 }];
    for (e, c) in f.iter() {
        v[n - e[0] as usize] = k.from_rat(c);
    }
    rr_strip(k, v)
}
fn rr_sparse(k: Inexact, f: &RV) -> MPoly<Rat> {
    let n = f.len();
    MPoly::from_terms(
        1,
        f.iter()
            .enumerate()
            .filter(|(_, c)| !k.is_zero(c))
            .map(|(i, c)| {
                let r = c.re.to_ratio_exact();
                (vec![(n - 1 - i) as u32], Rat::new(Int::from_big(r.numer().clone()), Int::from_big(r.denom().clone())))
            })
            .collect(),
    )
}
/// `dup_sub`.
fn rr_sub(k: Inexact, f: &RV, g: &RV) -> RV {
    let (df, dg) = (f.len(), g.len());
    let mut h = Vec::with_capacity(df.max(dg));
    if df > dg {
        h.extend(f[..df - dg].iter().cloned());
        h.extend(f[df - dg..].iter().zip(g).map(|(a, b)| k.sub(a, b)));
    } else {
        h.extend(g[..dg - df].iter().map(|b| k.neg(b)));
        h.extend(f.iter().zip(&g[dg - df..]).map(|(a, b)| k.sub(a, b)));
    }
    rr_strip(k, h)
}
/// `dup_mul` (the schoolbook product of the reference below degree 100).
fn rr_mul(k: Inexact, f: &RV, g: &RV) -> RV {
    if f.is_empty() || g.is_empty() {
        return vec![];
    }
    let (df, dg) = (f.len() - 1, g.len() - 1);
    let mut h = Vec::with_capacity(df + dg + 1);
    for i in 0..=df + dg {
        let mut c = k.zero();
        for j in i.saturating_sub(dg)..=df.min(i) {
            c = k.add(&c, &k.mul(&f[j], &g[i - j]));
        }
        h.push(c);
    }
    rr_strip(k, h)
}
/// `dup_ff_div` over an inexact field.
fn rr_div(k: Inexact, f: &RV, g: &RV) -> PResult<(RV, RV)> {
    if g.is_empty() {
        return Err(zero_division());
    }
    let dg = g.len() as i64 - 1;
    let mut dr = f.len() as i64 - 1;
    let (mut q, mut r): (RV, RV) = (vec![], f.clone());
    if dr < dg {
        return Ok((q, r));
    }
    let lc_g = g[0].clone();
    loop {
        let c = k.div(&r[0], &lc_g)?;
        let j = (dr - dg) as usize;
        if q.is_empty() {
            q = vec![k.zero(); j + 1];
            q[0] = c.clone();
        } else {
            let m = q.len() - j - 1;
            q[m] = k.add(&q[m], &c);
        }
        let mut h: RV = g.iter().map(|a| k.mul(a, &c)).collect();
        h.extend(std::iter::repeat(k.zero()).take(j));
        r = rr_sub(k, &r, &h);
        let prev = dr;
        dr = r.len() as i64 - 1;
        if dr < dg {
            break;
        } else if dr == prev {
            r = rr_strip(k, r[1..].to_vec());
            dr = r.len() as i64 - 1;
            if dr < dg {
                break;
            }
        } else if dr > prev {
            return Err(reference("PolynomialDivisionFailed", "polynomial division failed"));
        }
    }
    Ok((q, r))
}
fn rr_monic(k: Inexact, f: &RV) -> PResult<RV> {
    if f.is_empty() || k.is_zero(&k.sub(&f[0], &k.one())) {
        return Ok(f.clone());
    }
    f.iter().map(|c| k.div(c, &f[0])).collect()
}
fn rr_half_gcdex_dense(k: Inexact, f: &RV, g: &RV) -> PResult<(RV, RV)> {
    let (mut f, mut g) = (f.clone(), g.clone());
    let (mut a, mut b): (RV, RV) = (vec![k.one()], vec![]);
    while !g.is_empty() {
        let (q, r) = rr_div(k, &f, &g)?;
        f = std::mem::replace(&mut g, r);
        let nb = rr_sub(k, &a, &rr_mul(k, &q, &b));
        a = std::mem::replace(&mut b, nb);
    }
    let lc = f.first().cloned().unwrap_or_else(|| k.zero());
    let a = a.iter().map(|c| k.div(c, &lc)).collect::<PResult<RV>>()?;
    Ok((a, rr_monic(k, &f)?))
}
/// `dup_half_gcdex` over `RR`.
pub fn rr_half_gcdex(f: &MPoly<Rat>, g: &MPoly<Rat>) -> PResult<(MPoly<Rat>, MPoly<Rat>)> {
    rr_half_gcdex_prec(f, g, 53)
}
/// `dup_half_gcdex` over `RealField(precision)`.
pub fn rr_half_gcdex_prec(f: &MPoly<Rat>, g: &MPoly<Rat>, precision: u32) -> PResult<(MPoly<Rat>, MPoly<Rat>)> {
    let k = Inexact { precision, complex: false };
    let (s, h) = rr_half_gcdex_dense(k, &rr_dense(k, f), &rr_dense(k, g))?;
    Ok((rr_sparse(k, &s), rr_sparse(k, &h)))
}
/// `dup_gcdex` over `RR`.
pub fn rr_gcdex(f: &MPoly<Rat>, g: &MPoly<Rat>) -> PResult<(MPoly<Rat>, MPoly<Rat>, MPoly<Rat>)> {
    rr_gcdex_prec(f, g, 53)
}
/// `dup_gcdex` over `RealField(precision)`.
pub fn rr_gcdex_prec(f: &MPoly<Rat>, g: &MPoly<Rat>, precision: u32) -> PResult<(MPoly<Rat>, MPoly<Rat>, MPoly<Rat>)> {
    let k = Inexact { precision, complex: false };
    let (a, b) = (rr_dense(k, f), rr_dense(k, g));
    let (s, h) = rr_half_gcdex_dense(k, &a, &b)?;
    let big_f = rr_sub(k, &h, &rr_mul(k, &s, &a));
    // The reference takes `dup_quo` here: the rounding remainder is dropped.
    let (t, _) = rr_div(k, &big_f, &b)?;
    Ok((rr_sparse(k, &s), rr_sparse(k, &t), rr_sparse(k, &h)))
}
/// `dup_invert` over `RR`.
pub fn rr_invert(f: &MPoly<Rat>, g: &MPoly<Rat>) -> PResult<MPoly<Rat>> {
    rr_invert_prec(f, g, 53)
}
/// `dup_invert` over `RealField(precision)`.
pub fn rr_invert_prec(f: &MPoly<Rat>, g: &MPoly<Rat>, precision: u32) -> PResult<MPoly<Rat>> {
    let k = Inexact { precision, complex: false };
    let b = rr_dense(k, g);
    let (s, h) = rr_half_gcdex_dense(k, &rr_dense(k, f), &b)?;
    if h.len() == 1 && k.is_zero(&k.sub(&h[0], &k.one())) {
        Ok(rr_sparse(k, &rr_div(k, &s, &b)?.1))
    } else {
        Err(PolyError::NotInvertible("zero divisor".into()))
    }
}
