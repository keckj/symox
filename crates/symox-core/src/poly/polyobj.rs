//! The high level polynomial object (SymPy's `Poly`): a sparse polynomial
//! over `QQ` together with its generators and its ground domain.
//!
//! Coefficient domains like `ZZ[y]` are represented by extra variables
//! ("domain generators") after the true generators.

use std::cmp::Ordering;

use super::convert::{self, from_exprs, from_exprs_gens, to_expr_f, PResult, PolyError};
use super::factor_multi::factor_list_mz;
use super::gcd::gcd_mz;
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::sqf::sqf_list_mz;
use crate::expr::{Expr, Kind};
use crate::printer::to_str;

#[derive(Clone, Copy, PartialEq, Eq, Debug, PartialOrd, Ord)]
pub enum Domain {
    ZZ,
    QQ,
    RR,
}

#[derive(Clone, Debug)]
pub struct Poly {
    pub rep: MPoly<Rat>,
    pub gens: Vec<Expr>,
    pub ngens: usize,
    pub domain: Domain,
    /// `Some(p)` for polynomials over `GF(p)` (symmetric representatives).
    pub modulus: Option<u64>,
}

/// Symmetric representatives modulo `p` of the coefficients.
pub fn reduce_mod(rep: &MPoly<Rat>, p: u64) -> MPoly<Rat> {
    rep.map_coeffs(|c| {
        let n = c.n.rem_u64(p);
        let d = c.d.rem_u64(p);
        let v = if d == 1 { n } else { super::gf::mulmod(n, super::gf::invmod(d, p), p) };
        let s = if v > p / 2 { v as i64 - p as i64 } else { v as i64 };
        Rat::from_i64(s)
    })
}

fn natural_domain(p: &MPoly<Rat>, float: bool) -> Domain {
    if float {
        Domain::RR
    } else if p.is_integral() {
        Domain::ZZ
    } else {
        Domain::QQ
    }
}

/// SymPy's `_unify_gens`.
pub fn unify_gens(f: &[Expr], g: &[Expr]) -> Vec<Expr> {
    if f == g {
        return f.to_vec();
    }
    let mut f: Vec<Expr> = f.to_vec();
    let mut g: Vec<Expr> = g.to_vec();
    let common: Vec<Expr> = f.iter().filter(|x| g.contains(x)).cloned().collect();
    let mut k = 0;
    for gi in g.iter_mut() {
        if common.contains(gi) {
            *gi = common[k].clone();
            k += 1;
        }
    }
    let mut gens = Vec::new();
    for c in &common {
        let i = f.iter().position(|x| x == c).unwrap();
        gens.extend(f.drain(..i));
        f.remove(0);
        let i = g.iter().position(|x| x == c).unwrap();
        gens.extend(g.drain(..i));
        g.remove(0);
        gens.push(c.clone());
    }
    gens.extend(f);
    gens.extend(g);
    gens
}

impl Poly {
    pub fn new(rep: MPoly<Rat>, gens: Vec<Expr>, ngens: usize, domain: Domain) -> Poly {
        debug_assert_eq!(rep.nvars, gens.len());
        Poly { rep, gens, ngens, domain, modulus: None }
    }

    /// `Poly(expr, *gens)`.
    pub fn from_expr(e: &Expr, gens: &[Expr]) -> PResult<Poly> {
        let c = if gens.is_empty() { from_exprs(std::slice::from_ref(e))? } else { from_exprs_gens(std::slice::from_ref(e), gens)? };
        if c.ngens == 0 {
            return Err(PolyError::Failed(format!(
                "can't initialize from '{}' without generators",
                to_str(e)
            )));
        }
        let rep = c.polys.into_iter().next().unwrap();
        let domain = natural_domain(&rep, c.float);
        Ok(Poly { rep, gens: c.gens, ngens: c.ngens, domain, modulus: None })
    }

    /// Several expressions over a common ring (`parallel_poly_from_expr`).
    pub fn from_exprs(es: &[Expr], gens: &[Expr]) -> PResult<Vec<Poly>> {
        let c = if gens.is_empty() { from_exprs(es)? } else { from_exprs_gens(es, gens)? };
        if c.ngens == 0 {
            return Err(PolyError::Failed("can't initialize without generators".into()));
        }
        let mut dom = Domain::ZZ;
        for p in &c.polys {
            dom = dom.max(natural_domain(p, c.float));
        }
        Ok(c.polys
            .into_iter()
            .map(|rep| Poly { rep, gens: c.gens.clone(), ngens: c.ngens, domain: dom, modulus: None })
            .collect())
    }

    /// The same polynomial over `GF(p)`.
    pub fn set_modulus(&self, p: u64) -> Poly {
        Poly {
            rep: reduce_mod(&self.rep, p),
            gens: self.gens.clone(),
            ngens: self.ngens,
            domain: Domain::ZZ,
            modulus: Some(p),
        }
    }

    fn require_gf(&self) -> PResult<u64> {
        let p = self.modulus.unwrap();
        if self.gens.len() != 1 {
            return Err(PolyError::Unsupported("multivariate polynomials over finite fields".into()));
        }
        Ok(p)
    }

    /// Dense coefficients in `[0, p)` (univariate).
    pub fn to_gf(&self) -> Vec<u64> {
        let p = self.modulus.unwrap();
        let d = super::dense::from_mpoly(&self.rep, 0);
        let mut v: Vec<u64> = d.iter().map(|c| c.n.rem_u64(p)).collect();
        super::gf::trim(&mut v);
        v
    }

    pub fn from_gf(&self, v: &[u64]) -> Poly {
        let p = self.modulus.unwrap();
        let ints = super::gf::to_int_sym(v, p);
        let rats = super::dense::to_rat(&ints);
        let mut rr = rats;
        super::dense::trim(&mut rr);
        self.with_rep(super::dense::to_mpoly(&rr, 0, self.gens.len()))
    }

    pub fn with_rep(&self, rep: MPoly<Rat>) -> Poly {
        let rep = match self.modulus {
            Some(p) => reduce_mod(&rep, p),
            None => rep,
        };
        Poly { rep, gens: self.gens.clone(), ngens: self.ngens, domain: self.domain, modulus: self.modulus }
    }

    pub fn true_gens(&self) -> &[Expr] {
        &self.gens[..self.ngens]
    }

    pub fn domain_gens(&self) -> &[Expr] {
        &self.gens[self.ngens..]
    }

    pub fn is_float(&self) -> bool {
        self.domain == Domain::RR
    }

    pub fn as_expr(&self) -> Expr {
        to_expr_f(&self.rep, &self.gens, self.is_float())
    }

    pub fn is_zero(&self) -> bool {
        self.rep.is_zero()
    }

    pub fn domain_str(&self) -> String {
        let base = match self.domain {
            Domain::ZZ => "ZZ",
            Domain::QQ => "QQ",
            Domain::RR => "RR",
        };
        if self.ngens == self.gens.len() {
            base.to_string()
        } else {
            let ds: Vec<String> = self.domain_gens().iter().map(to_str).collect();
            format!("{base}[{}]", ds.join(","))
        }
    }

    /// Degree in the `i`-th generator (`None` for the zero polynomial,
    /// i.e. `-oo`).
    pub fn degree(&self, i: usize) -> Option<u32> {
        if self.rep.is_zero() {
            None
        } else {
            Some(self.rep.degree(i) as u32)
        }
    }

    pub fn gen_index(&self, g: &Expr) -> Option<usize> {
        self.true_gens().iter().position(|x| x == g)
    }

    pub fn total_degree(&self) -> Option<u32> {
        if self.rep.is_zero() {
            return None;
        }
        Some(
            self.rep
                .iter()
                .map(|(e, _)| e[..self.ngens].iter().sum::<u32>())
                .max()
                .unwrap(),
        )
    }

    /// Terms `(monomial in the true generators, coefficient)` in
    /// decreasing lex order.
    pub fn terms(&self) -> Vec<(Vec<u32>, Expr)> {
        let k = self.ngens;
        let n = self.gens.len();
        let mut out: Vec<(Vec<u32>, Vec<(Vec<u32>, Rat)>)> = Vec::new();
        for (e, c) in self.rep.iter() {
            let m = e[..k].to_vec();
            let mut de = vec![0u32; n];
            de[k..].copy_from_slice(&e[k..]);
            if out.last().map(|x| &x.0) != Some(&m) {
                out.push((m, Vec::new()));
            }
            out.last_mut().unwrap().1.push((de, c.clone()));
        }
        out.into_iter()
            .map(|(m, ts)| {
                let p = MPoly::from_terms(n, ts);
                (m, to_expr_f(&p, &self.gens, self.is_float()))
            })
            .collect()
    }

    pub fn coeffs(&self) -> Vec<Expr> {
        let t = self.terms();
        if t.is_empty() {
            return vec![Expr::zero()];
        }
        t.into_iter().map(|x| x.1).collect()
    }

    pub fn monoms(&self) -> Vec<Vec<u32>> {
        let t = self.terms();
        if t.is_empty() {
            return vec![vec![0; self.ngens]];
        }
        t.into_iter().map(|x| x.0).collect()
    }

    /// Dense coefficients, highest degree first (univariate only).
    pub fn all_coeffs(&self) -> PResult<Vec<Expr>> {
        if self.ngens != 1 {
            return Err(PolyError::Unsupported("multivariate polynomials not supported".into()));
        }
        let t = self.terms();
        if t.is_empty() {
            return Ok(vec![Expr::zero()]);
        }
        let d = t[0].0[0] as usize;
        let mut v = vec![Expr::zero(); d + 1];
        for (m, c) in t {
            v[d - m[0] as usize] = c;
        }
        Ok(v)
    }

    /// Leading coefficient.
    pub fn lc(&self) -> Expr {
        self.terms().into_iter().next().map(|x| x.1).unwrap_or_else(Expr::zero)
    }

    /// Trailing coefficient.
    pub fn tc(&self) -> Expr {
        match self.terms().into_iter().last() {
            Some((m, c)) if m.iter().all(|&e| e == 0) => c,
            _ => Expr::zero(),
        }
    }

    /// Coefficient of a monomial in the true generators.
    pub fn coeff_monomial(&self, m: &[u32]) -> Expr {
        for (mm, c) in self.terms() {
            if mm == m {
                return c;
            }
        }
        Expr::zero()
    }

    /// Bring two polynomials to common generators and domain.
    pub fn unify(&self, o: &Poly) -> (Poly, Poly) {
        if self.gens == o.gens && self.ngens == o.ngens {
            let d = self.domain.max(o.domain);
            let m = self.modulus.or(o.modulus);
            let mut a = self.clone();
            let mut b = o.clone();
            a.domain = d;
            b.domain = d;
            if m.is_some() && (a.modulus.is_none() || b.modulus.is_none()) {
                a = a.set_modulus(m.unwrap());
                b = b.set_modulus(m.unwrap());
            }
            return (a, b);
        }
        let tg = unify_gens(self.true_gens(), o.true_gens());
        let mut dg: Vec<Expr> = Vec::new();
        for g in self.domain_gens().iter().chain(o.domain_gens()) {
            if !tg.contains(g) && !dg.contains(g) {
                dg.push(g.clone());
            }
        }
        convert::sort_gens(&mut dg);
        let ngens = tg.len();
        let mut all = tg;
        all.extend(dg);
        let a = convert::reembed(&self.rep, &self.gens, &all);
        let b = convert::reembed(&o.rep, &o.gens, &all);
        let d = self.domain.max(o.domain);
        (
            Poly { rep: a, gens: all.clone(), ngens, domain: d, modulus: self.modulus.or(o.modulus) },
            Poly { rep: b, gens: all, ngens, domain: d, modulus: self.modulus.or(o.modulus) },
        )
    }

    pub fn add(&self, o: &Poly) -> Poly {
        let (a, b) = self.unify(o);
        a.with_rep(a.rep.add(&b.rep))
    }

    pub fn sub(&self, o: &Poly) -> Poly {
        let (a, b) = self.unify(o);
        a.with_rep(a.rep.sub(&b.rep))
    }

    pub fn mul(&self, o: &Poly) -> Poly {
        let (a, b) = self.unify(o);
        a.with_rep(a.rep.mul(&b.rep))
    }

    pub fn neg(&self) -> Poly {
        self.with_rep(self.rep.neg())
    }

    pub fn pow(&self, k: u32) -> Poly {
        self.with_rep(self.rep.pow(k))
    }

    pub fn mul_ground(&self, c: &Rat) -> Poly {
        let mut p = self.with_rep(self.rep.scale(c));
        if !c.is_integer() && p.domain == Domain::ZZ {
            p.domain = Domain::QQ;
        }
        p
    }

    /// Derivative with respect to the `i`-th generator.
    pub fn diff(&self, i: usize) -> Poly {
        self.with_rep(self.rep.derivative(i))
    }

    pub fn has_domain_gens(&self) -> bool {
        self.ngens != self.gens.len()
    }

    /// Polynomial division with remainder (recursive in the generators,
    /// over the field of fractions of the domain). As in SymPy, results
    /// over `ZZ` are retracted to `ZZ` when their coefficients allow it.
    pub fn div(&self, o: &Poly) -> PResult<(Poly, Poly)> {
        let (a, b) = self.unify(o);
        if b.rep.is_zero() {
            return Err(PolyError::Unsupported("polynomial division by zero".into()));
        }
        if a.modulus.is_some() {
            let p = a.require_gf()?;
            let (q, r) = super::gf::divrem(&a.to_gf(), &b.to_gf(), p);
            return Ok((a.from_gf(&q), a.from_gf(&r)));
        }
        let (q, r) = div_rec(&a.rep, &b.rep, 0);
        if a.has_domain_gens() && !r.is_zero() && r.degree(0) >= b.rep.degree(0) {
            // the quotient over the fraction field of the coefficient ring
            // would need rational function coefficients
            return Err(PolyError::Unsupported(format!(
                "division over {} with non-polynomial quotient is not supported",
                a.domain_str()
            )));
        }
        let mut dom = a.domain;
        if dom == Domain::ZZ && !(q.is_integral() && r.is_integral()) {
            dom = Domain::QQ;
        }
        let mk = |rep: MPoly<Rat>| Poly { rep, gens: a.gens.clone(), ngens: a.ngens, domain: dom, modulus: None };
        Ok((mk(q), mk(r)))
    }

    pub fn rem(&self, o: &Poly) -> PResult<Poly> {
        if self.modulus.is_some() || o.modulus.is_some() {
            return Ok(self.div(o)?.1);
        }
        let (_, mut r) = self.div(o)?;
        if self.domain.max(o.domain) == Domain::ZZ && r.rep.is_integral() {
            r.domain = Domain::ZZ;
        }
        Ok(r)
    }

    pub fn quo(&self, o: &Poly) -> PResult<Poly> {
        if self.modulus.is_some() || o.modulus.is_some() {
            return Ok(self.div(o)?.0);
        }
        let (mut q, _) = self.div(o)?;
        if self.domain.max(o.domain) == Domain::ZZ && q.rep.is_integral() {
            q.domain = Domain::ZZ;
        }
        Ok(q)
    }

    /// Pseudo-division with respect to the first generator:
    /// `lc(g)**(deg f - deg g + 1) * f = q*g + r`.
    pub fn pdiv(&self, o: &Poly) -> PResult<(Poly, Poly)> {
        let (a, b) = self.unify(o);
        if b.rep.is_zero() {
            return Err(PolyError::Unsupported("polynomial division by zero".into()));
        }
        let n = a.gens.len();
        let f = &a.rep;
        let g = &b.rep;
        let db = g.degree(0);
        let lb = g.lc_in(0);
        let mut q = MPoly::zero(n);
        let mut r = f.clone();
        let mut steps = (f.degree(0) - db + 1).max(0) as u32;
        while !r.is_zero() && r.degree(0) >= db {
            let j = (r.degree(0) - db) as u32;
            let lr = r.lc_in(0);
            let mut e = vec![0u32; n];
            e[0] = j;
            let t = lr.mul_term(&e, &Rat::one());
            q = q.mul(&lb).add(&t);
            r = r.mul(&lb).sub(&t.mul(g));
            steps = steps.saturating_sub(1);
        }
        if steps > 0 {
            let c = lb.pow(steps);
            q = q.mul(&c);
            r = r.mul(&c);
        }
        Ok((a.with_rep(q), a.with_rep(r)))
    }

    /// Exact quotient (error if not exact).
    pub fn exquo(&self, o: &Poly) -> PResult<Poly> {
        let (a, b) = self.unify(o);
        match a.rep.div_exact(&b.rep) {
            Some(q) if b.domain != Domain::ZZ || q.is_integral() => Ok(a.with_rep(q)),
            _ => Err(PolyError::Unsupported(format!(
                "{} does not divide {}",
                to_str(&o.as_expr()),
                to_str(&self.as_expr())
            ))),
        }
    }

    /// Gcd (positive leading coefficient over `ZZ`, monic over `QQ`).
    pub fn gcd(&self, o: &Poly) -> Poly {
        let (a, b) = self.unify(o);
        if let Some(p) = a.modulus {
            return a.from_gf(&super::gf::gcd(&a.to_gf(), &b.to_gf(), p));
        }
        a.with_rep(gcd_rep(&a.rep, &b.rep, a.domain))
    }

    pub fn lcm(&self, o: &Poly) -> Poly {
        let (a, b) = self.unify(o);
        if let Some(p) = a.modulus {
            let (x, y) = (a.to_gf(), b.to_gf());
            if x.is_empty() || y.is_empty() {
                return a.from_gf(&[]);
            }
            let g = super::gf::gcd(&x, &y, p);
            let l = super::gf::divrem(&super::gf::mul(&x, &y, p), &g, p).0;
            return a.from_gf(&super::gf::monic(&l, p));
        }
        a.with_rep(lcm_rep(&a.rep, &b.rep, a.domain))
    }

    /// `(cofactor gcd, f/gcd, g/gcd)`.
    pub fn cofactors(&self, o: &Poly) -> (Poly, Poly, Poly) {
        let (a, b) = self.unify(o);
        let h = gcd_rep(&a.rep, &b.rep, a.domain);
        if h.is_zero() {
            return (a.with_rep(h), a.with_rep(MPoly::zero(a.gens.len())), a.with_rep(MPoly::zero(a.gens.len())));
        }
        let ca = a.rep.div_exact(&h).unwrap();
        let cb = b.rep.div_exact(&h).unwrap();
        (a.with_rep(h), a.with_rep(ca), a.with_rep(cb))
    }

    /// Integer content and primitive part.
    pub fn primitive(&self) -> (Rat, Poly) {
        if self.rep.is_zero() {
            return (Rat::zero(), self.clone());
        }
        let (d, z) = self.rep.clear_denoms();
        let c = z.content();
        let pp = z.div_coeff(&c).unwrap();
        (Rat::new(c, d), self.with_rep(pp.to_rat()))
    }

    /// `(coefficient, [(factor, multiplicity)])` sorted like SymPy.
    pub fn factor_list(&self) -> (Rat, Vec<(Poly, u32)>) {
        if let Some(p) = self.modulus {
            return self.gf_list(p, false);
        }
        let (d, z) = self.rep.clear_denoms();
        let (c, fs) = factor_list_mz(&z);
        let coeff = Rat::new(c, d);
        let mut out: Vec<(Poly, u32)> = fs.into_iter().map(|(f, k)| (self.with_rep(f.to_rat()), k)).collect();
        sort_factors_poly(&mut out, self.ngens);
        (coeff, out)
    }

    /// Factorisation / square-free decomposition over `GF(p)` (factors
    /// sorted by their coefficient lists in `[0, p)`, like SymPy).
    fn gf_list(&self, p: u64, sqf: bool) -> (Rat, Vec<(Poly, u32)>) {
        let f = self.to_gf();
        let (lc, mut fs) = if sqf { super::gf::sqf_list(&f, p) } else { super::gf::factor_list(&f, p) };
        fs.sort_by(|(a, ka), (b, kb)| {
            a.len().cmp(&b.len()).then(ka.cmp(kb)).then_with(|| a.iter().rev().cmp(b.iter().rev()))
        });
        let lcs = if lc > p / 2 { lc as i64 - p as i64 } else { lc as i64 };
        (Rat::from_i64(lcs), fs.into_iter().map(|(g, k)| (self.from_gf(&g), k)).collect())
    }

    /// Square-free decomposition.
    pub fn sqf_list(&self) -> (Rat, Vec<(Poly, u32)>) {
        if let Some(p) = self.modulus {
            return self.gf_list(p, true);
        }
        let (d, z) = self.rep.clear_denoms();
        let (c, mut parts) = sqf_list_mz(&z);
        // powers of the variables are part of the decomposition in SymPy
        // (they are handled by the general algorithm)
        parts.retain(|(p, _)| !p.is_constant());
        let coeff = Rat::new(c, d);
        let mut out: Vec<(Poly, u32)> = parts.into_iter().map(|(f, k)| (self.with_rep(f.to_rat()), k)).collect();
        sort_factors_poly(&mut out, self.ngens);
        (coeff, out)
    }

    /// `f(g)` substituting `g` for the first generator.
    pub fn compose(&self, g: &Poly) -> Poly {
        let (a, b) = self.unify(g);
        let n = a.gens.len();
        let cs = a.rep.coeffs_in(0);
        let mut acc = MPoly::zero(n);
        let top = cs.first().map(|x| x.0).unwrap_or(0);
        let mut idx = 0;
        for k in (0..=top).rev() {
            acc = acc.mul(&b.rep);
            if idx < cs.len() && cs[idx].0 == k {
                acc = acc.add(&cs[idx].1);
                idx += 1;
            }
        }
        a.with_rep(acc)
    }

    /// Functional decomposition (univariate).
    pub fn decompose(&self) -> PResult<Vec<Poly>> {
        if self.gens.len() != 1 {
            return Err(PolyError::Unsupported("univariate polynomial expected".into()));
        }
        let d = super::dense::from_mpoly(&self.rep, 0);
        if d.len() <= 2 {
            return Ok(vec![self.clone()]);
        }
        let fs = super::roots::decompose(&d);
        Ok(fs
            .into_iter()
            .map(|c| {
                let mut p = self.with_rep(super::dense::to_mpoly(&c, 0, 1));
                if p.domain == Domain::ZZ && !p.rep.is_integral() {
                    p.domain = Domain::QQ;
                }
                p
            })
            .collect())
    }

    /// `f(x + a)` in the first generator.
    pub fn shift(&self, a: &Rat) -> Poly {
        let mut p = self.with_rep(self.rep.shift(0, a));
        if !a.is_integer() && p.domain == Domain::ZZ {
            p.domain = Domain::QQ;
        }
        p
    }

    /// Antiderivative with respect to the `i`-th generator.
    pub fn integrate(&self, i: usize) -> Poly {
        let n = self.gens.len();
        let mut terms = Vec::new();
        for (e, c) in self.rep.iter() {
            let mut ne = e.to_vec();
            ne[i] += 1;
            let k = ne[i] as i64;
            terms.push((ne, c.div(&Rat::from_i64(k))));
        }
        let mut p = self.with_rep(MPoly::from_terms(n, terms));
        if p.domain == Domain::ZZ {
            p.domain = Domain::QQ;
        }
        p
    }

    /// Product of the distinct square-free factors.
    pub fn sqf_part(&self) -> Poly {
        let (_, fs) = self.sqf_list();
        let mut acc = MPoly::one(self.gens.len());
        for (f, _) in fs {
            acc = acc.mul(&f.rep);
        }
        self.with_rep(acc)
    }

    pub fn max_norm(&self) -> Rat {
        self.rep.coeffs().iter().map(|c| c.abs()).max().unwrap_or_else(Rat::zero)
    }

    pub fn l1_norm(&self) -> Rat {
        self.rep.coeffs().iter().fold(Rat::zero(), |a, c| a.add(&c.abs()))
    }

    /// `(d, d*f)` with integer coefficients.
    pub fn clear_denoms(&self) -> (Int, Poly) {
        let (d, z) = self.rep.clear_denoms();
        (d, self.with_rep(z.to_rat()))
    }

    /// Total degrees (in the generators) of the terms.
    fn term_degrees(&self) -> Vec<u32> {
        self.terms().iter().map(|(m, _)| m.iter().sum()).collect()
    }

    pub fn is_linear(&self) -> bool {
        self.term_degrees().iter().all(|&d| d <= 1)
    }

    pub fn is_quadratic(&self) -> bool {
        self.term_degrees().iter().all(|&d| d <= 2)
    }

    pub fn is_homogeneous(&self) -> bool {
        let d = self.term_degrees();
        d.windows(2).all(|w| w[0] == w[1])
    }

    /// SymPy's `str(Poly)`.
    pub fn to_sympy_string(&self) -> String {
        let atom = |g: &Expr| matches!(g.kind(), Kind::Sym(_) | Kind::Const(_));
        let gens_s: Vec<String> = self
            .true_gens()
            .iter()
            .map(|g| if atom(g) { to_str(g) } else { format!("({})", to_str(g)) })
            .collect();
        let mut terms: Vec<String> = Vec::new();
        let tm = self.terms();
        let tm = if tm.is_empty() { vec![(vec![0; self.ngens], Expr::zero())] } else { tm };
        for (m, c) in tm {
            let mut sm = Vec::new();
            for (i, &e) in m.iter().enumerate() {
                if e == 1 {
                    sm.push(gens_s[i].clone());
                } else if e > 1 {
                    sm.push(format!("{}**{}", gens_s[i], e));
                }
            }
            let s_monom = sm.join("*");
            let s_coeff;
            if c.is_add() {
                s_coeff = if s_monom.is_empty() { to_str(&c) } else { format!("({})", to_str(&c)) };
            } else {
                if !s_monom.is_empty() {
                    if c.is_one() {
                        terms.push("+".into());
                        terms.push(s_monom);
                        continue;
                    }
                    if c.is_minus_one() {
                        terms.push("-".into());
                        terms.push(s_monom);
                        continue;
                    }
                }
                s_coeff = match c.kind() {
                    Kind::Num(crate::number::Number::Float(f)) => crate::number::format_float(*f),
                    _ => to_str(&c),
                };
            }
            let s_term = if s_monom.is_empty() { s_coeff } else { format!("{s_coeff}*{s_monom}") };
            if let Some(rest) = s_term.strip_prefix('-') {
                terms.push("-".into());
                terms.push(rest.to_string());
            } else {
                terms.push("+".into());
                terms.push(s_term);
            }
        }
        if terms[0] == "-" || terms[0] == "+" {
            let m = terms.remove(0);
            if m == "-" {
                terms[0] = format!("-{}", terms[0]);
            }
        }
        let gl: Vec<String> = self.true_gens().iter().map(to_str).collect();
        if let Some(p) = self.modulus {
            return format!("Poly({}, {}, modulus={p})", terms.join(" "), gl.join(", "));
        }
        format!(
            "Poly({}, {}, domain='{}')",
            terms.join(" "),
            gl.join(", "),
            self.domain_str()
        )
    }

    /// Value at `x_i = a` (rational), as a polynomial in the remaining
    /// generators (or a ground element).
    pub fn eval_at(&self, i: usize, a: &Rat) -> Poly {
        let r = self.rep.eval(i, a);
        let mut gens = self.gens.clone();
        gens.remove(i);
        let perm: Vec<usize> = (0..self.gens.len()).map(|v| if v < i { v } else if v > i { v - 1 } else { 0 }).collect();
        let rep = r.permute(&perm, gens.len());
        let mut dom = self.domain;
        if !a.is_integer() && dom == Domain::ZZ {
            dom = Domain::QQ;
        }
        Poly { rep, gens, ngens: self.ngens - 1, domain: dom, modulus: self.modulus }
    }

    /// Structural equality (same generators and representation).
    pub fn equals(&self, o: &Poly) -> bool {
        self.gens == o.gens && self.ngens == o.ngens && self.rep == o.rep
    }
}

/// Recursive division over a field (`dmp_ff_div`): quotient terms are
/// produced while the leading coefficients divide.
pub fn div_rec(f: &MPoly<Rat>, g: &MPoly<Rat>, v: usize) -> (MPoly<Rat>, MPoly<Rat>) {
    let n = f.nvars;
    if v >= n || (f.is_constant() && g.is_constant()) {
        let gc = g.constant_value().unwrap_or_else(Rat::zero);
        if let Some(fc) = f.constant_value() {
            if !gc.is_zero() {
                return (MPoly::constant(fc.div(&gc), n), MPoly::zero(n));
            }
        }
        if v >= n {
            return (MPoly::zero(n), f.clone());
        }
    }
    let dg = g.degree(v);
    let mut q = MPoly::zero(n);
    let mut r = f.clone();
    if r.is_zero() || r.degree(v) < dg {
        return (q, r);
    }
    let lcg = g.lc_in(v);
    loop {
        let dr = r.degree(v);
        let lcr = r.lc_in(v);
        let (c, rr) = div_rec(&lcr, &lcg, v + 1);
        if !rr.is_zero() || c.is_zero() {
            break;
        }
        let mut e = vec![0u32; n];
        e[v] = (dr - dg) as u32;
        let t = c.mul_term(&e, &Rat::one());
        q = q.add(&t);
        r = r.sub(&t.mul(g));
        if r.is_zero() || r.degree(v) < dg {
            break;
        }
        if r.degree(v) >= dr {
            break;
        }
    }
    (q, r)
}

pub fn gcd_rep(a: &MPoly<Rat>, b: &MPoly<Rat>, dom: Domain) -> MPoly<Rat> {
    let n = a.nvars;
    if a.is_zero() && b.is_zero() {
        return MPoly::zero(n);
    }
    let (_, za) = a.clear_denoms();
    let (_, zb) = b.clear_denoms();
    let (h, _, _) = gcd_mz(&za, &zb);
    let h = h.to_rat();
    if dom == Domain::ZZ && a.is_integral() && b.is_integral() {
        h
    } else {
        h.monic()
    }
}

pub fn lcm_rep(a: &MPoly<Rat>, b: &MPoly<Rat>, dom: Domain) -> MPoly<Rat> {
    let n = a.nvars;
    if a.is_zero() || b.is_zero() {
        return MPoly::zero(n);
    }
    if dom == Domain::ZZ && a.is_integral() && b.is_integral() {
        let za = a.to_int().unwrap();
        let zb = b.to_int().unwrap();
        let (ca, pa) = za.primitive();
        let (cb, pb) = zb.primitive();
        let c = ca.lcm(&cb);
        let (_, ha, _) = gcd_mz(&pa, &pb);
        let l = ha.mul(&pb).scale(&c);
        let l = if l.lc().is_negative() { l.neg() } else { l };
        return l.to_rat();
    }
    let (h, _, _) = gcd_mz(&a.clear_denoms().1, &b.clear_denoms().1);
    let l = a.mul(b).div_exact(&h.to_rat()).unwrap();
    l.monic()
}

/// Compare dense recursive representations like Python compares nested
/// lists (`v` is the current level).
pub fn rep_cmp(a: &MPoly<Rat>, b: &MPoly<Rat>, v: usize) -> Ordering {
    let n = a.nvars;
    if v >= n {
        let x = a.constant_value().unwrap_or_else(Rat::zero);
        let y = b.constant_value().unwrap_or_else(Rat::zero);
        return x.cmp(&y);
    }
    let la = a.degree(v) + 1;
    let lb = b.degree(v) + 1;
    let m = la.min(lb);
    for i in 0..m {
        let ca = a.coeff_of(v, (la - 1 - i) as u32);
        let cb = b.coeff_of(v, (lb - 1 - i) as u32);
        let c = rep_cmp(&ca, &cb, v + 1);
        if c != Ordering::Equal {
            return c;
        }
    }
    la.cmp(&lb)
}

/// SymPy's `_sort_factors` (multiple=True): by `(len(rep), k, rep)`.
pub fn sort_factors_poly(fs: &mut [(Poly, u32)], _ngens: usize) {
    fs.sort_by(|(a, ka), (b, kb)| {
        let la = a.rep.degree(0) + 1;
        let lb = b.rep.degree(0) + 1;
        la.cmp(&lb).then(ka.cmp(kb)).then_with(|| rep_cmp(&a.rep, &b.rep, 0))
    });
}

/// SymPy's `_sorted_factors` for `factor_list`/`sqf_list` results.
pub fn sorted_factors(fs: &mut [(Poly, u32)], sqf: bool) {
    fs.sort_by(|(a, ka), (b, kb)| {
        let la = a.rep.degree(0) + 1;
        let lb = b.rep.degree(0) + 1;
        let main = if sqf {
            ka.cmp(kb).then(la.cmp(&lb)).then(a.gens.len().cmp(&b.gens.len()))
        } else {
            la.cmp(&lb).then(a.gens.len().cmp(&b.gens.len())).then(ka.cmp(kb))
        };
        main.then_with(|| a.domain_str().cmp(&b.domain_str()))
            .then_with(|| rep_cmp(&a.rep, &b.rep, 0))
    });
}

pub fn int_poly(p: &Poly) -> Option<MPoly<Int>> {
    p.rep.to_int()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn poly(s: &str) -> Poly {
        Poly::from_expr(&parse_str(s).unwrap(), &[]).unwrap()
    }

    #[test]
    fn printing() {
        assert_eq!(poly("x**2 + 1").to_sympy_string(), "Poly(x**2 + 1, x, domain='ZZ')");
        assert_eq!(poly("x/2 + 1").to_sympy_string(), "Poly(1/2*x + 1, x, domain='QQ')");
        assert_eq!(poly("sqrt(2)*x").to_sympy_string(), "Poly(x*(sqrt(2)), x, sqrt(2), domain='ZZ')");
        assert_eq!(poly("pi*x + sin(y)").to_sympy_string(), "Poly(x*pi + (sin(y)), x, pi, sin(y), domain='ZZ')");
        assert_eq!(poly("1/x + x").to_sympy_string(), "Poly(x + (1/x), x, 1/x, domain='ZZ')");
        assert_eq!(poly("x**2*y - 3*x + y/3").to_sympy_string(), "Poly(x**2*y - 3*x + 1/3*y, x, y, domain='QQ')");
        let p = Poly::from_expr(&parse_str("x*y + 1").unwrap(), &[Expr::symbol("x")]).unwrap();
        assert_eq!(p.to_sympy_string(), "Poly(y*x + 1, x, domain='ZZ[y]')");
        let p = Poly::from_expr(&parse_str("x**2 + 2*x*y + 1").unwrap(), &[Expr::symbol("y")]).unwrap();
        assert_eq!(p.to_sympy_string(), "Poly(2*x*y + x**2 + 1, y, domain='ZZ[x]')");
        assert_eq!(poly("2*x - 1.5").to_sympy_string(), "Poly(2.0*x - 1.5, x, domain='RR')");
        let a = poly("x");
        let b = Poly::from_expr(&parse_str("y").unwrap(), &[]).unwrap();
        assert_eq!(a.add(&b).to_sympy_string(), "Poly(x + y, x, y, domain='ZZ')");
    }

    #[test]
    fn division() {
        let (q, r) = poly("x**2 + 1").div(&poly("2*x")).unwrap();
        assert_eq!(q.to_sympy_string(), "Poly(1/2*x, x, domain='QQ')");
        assert_eq!(r.to_sympy_string(), "Poly(1, x, domain='QQ')");
    }

    #[test]
    fn factor_list_order() {
        let (c, fs) = poly("x**4*y - x**2*y**3").factor_list();
        assert_eq!(c, Rat::one());
        let s: Vec<String> = fs.iter().map(|(p, k)| format!("{}^{k}", to_str(&p.as_expr()))).collect();
        assert_eq!(s, vec!["y^1", "x - y^1", "x + y^1", "x^2"]);
    }
}
