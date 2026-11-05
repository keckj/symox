//! `is_log_deriv_k_t_radical_in_field` and `parametric_log_deriv_heu` of
//! SymPy `prde.py` on a differential tower `K = Q(params)(x, t_1, ..., t_k)`
//! (exp, primitive and base levels), over the rational function field of
//! every symbol involved ([`RF`]).
//!
//! The stock bodies are reproduced step by step:
//! - `residue_reduce` (subresultant PRS of `d` and `a - z*Dd` over
//!   `K[t][z]`, Yun square-free decomposition of the resultant over
//!   `K[z]`, the splitting `S_i = gcd(p_i, kappa(p_i))`, the `invert`
//!   normalization, which gives `h(rho) = S(t, rho)/lc_t(S)(rho)`);
//! - the rational root check, `p = f - sum(rho*Dh/h)`, `p in k[t]` and the
//!   degree bound;
//! - the exp, primitive and base branches (recursively, one level down) and
//!   the stock exceptions of the tan, other and unknown cases;
//! - `parametric_log_deriv_heu` (`solve` for one unknown, `splitfactor`).
//!
//! Every Dummy the stock bodies create is recorded in order in
//! [`Tower::trace`] (`z`, one unnamed Dummy per `real_roots` call, `i`,
//! `c1`), so the caller can create them as the stock body does.
//! A step that is not reproduced (possible raise of the stock body, inputs
//! outside the field) returns [`Stop::Decline`]; the caller then runs the
//! Python body unchanged.
use crate::integrate::prde_const::derivation;
use crate::poly::dense;
use crate::poly::gcd::{gcd_mz, gcd_mz_only};
use crate::poly::int::{Coeff, Field, Int, Rat};
use crate::poly::mpoly::MPoly;
use crate::poly::ratfunc::RF;

/// Why a kernel stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Stop {
    /// Not reproduced: run the stock body.
    Decline,
    /// The stock body raises `NotImplementedError(message)`.
    NotImplemented(String),
    /// The stock body raises `ValueError(message)`.
    Value(String),
}

pub type Out<T> = Result<T, Stop>;

/// A univariate polynomial in the variable of one level, coefficients in
/// the field below (free of that variable).
pub type UP = Vec<RF>;

/// One level of the tower: its variable index in the field, `D(t)` and
/// `DE.cases[level]`.
#[derive(Clone, Debug)]
pub struct Level {
    pub var: usize,
    pub d: RF,
    pub case: String,
}

pub struct Tower<'a> {
    levels: Vec<Level>,
    d_up: Vec<UP>,
    n: usize,
    /// Names of the Dummies created by the stock body, in order (`""` for
    /// `Dummy()`).
    pub trace: Vec<&'static str>,
    check: &'a mut dyn FnMut() -> bool,
}

fn konst(r: &Rat, n: usize) -> RF {
    RF::new(MPoly::constant(r.n.clone(), n), MPoly::constant(r.d.clone(), n))
}

fn int_rf(i: &Int, n: usize) -> RF {
    RF::from_poly(MPoly::constant(i.clone(), n))
}

/// The rational number of a constant `RF`.
fn rational(f: &RF) -> Option<Rat> {
    Some(Rat::new(f.n.constant_value()?, f.d.constant_value()?))
}

fn poly_up(p: &MPoly<Int>, v: usize) -> UP {
    let cs = p.coeffs_in(v);
    let Some(top) = cs.first().map(|c| c.0 as usize) else { return Vec::new() };
    let mut out = vec![RF::zero(); top + 1];
    for (k, c) in cs {
        out[k as usize] = RF::from_poly(c);
    }
    dense::trim(&mut out);
    out
}

/// `f` as a polynomial in `x_v` (`None` if its denominator involves `x_v`).
pub fn up_of(f: &RF, v: usize) -> Option<UP> {
    if f.d.degree(v) > 0 {
        return None;
    }
    let inv = RF::from_poly(f.d.clone()).inv();
    Some(poly_up(&f.n, v).iter().map(|c| c.mul(&inv)).collect())
}

/// Reduced numerator and denominator of `f` as polynomials in `x_v`.
fn split(f: &RF, v: usize) -> (UP, UP) {
    (poly_up(&f.n, v), poly_up(&f.d, v))
}

fn up_rf(u: &[RF], v: usize, n: usize) -> RF {
    // Horner in x_v.
    let x = RF::from_poly(MPoly::var(v, n));
    let mut acc = RF::zero();
    for c in u.iter().rev() {
        acc = acc.mul(&x).add(c);
    }
    acc
}

fn rf_pow(f: &RF, e: &Int) -> Option<RF> {
    let k = e.as_i64()?;
    let base = if k < 0 { f.inv() } else { f.clone() };
    let k = u32::try_from(k.unsigned_abs()).ok()?;
    Some(RF { n: base.n.pow(k), d: base.d.pow(k) })
}

fn deg(p: &[RF]) -> isize {
    dense::deg(p)
}

fn monic(p: &[RF]) -> UP {
    if p.is_empty() {
        return Vec::new();
    }
    dense::monic(p)
}

fn quo(a: &[RF], b: &[RF]) -> UP {
    dense::divrem(a, b).0
}

fn lcm_up(a: &[RF], b: &[RF]) -> UP {
    let g = dense::gcd_field(a, b);
    monic(&dense::mul(&quo(a, &g), b))
}

fn is_rational_poly(p: &[RF]) -> Option<Vec<Rat>> {
    p.iter().map(rational).collect()
}

/// Yun square-free decomposition over the coefficient field: monic
/// factors of positive degree with their multiplicities.
fn sqf_field(f: &[RF]) -> Vec<(UP, usize)> {
    let mut out = Vec::new();
    if deg(f) < 1 {
        return out;
    }
    let f = monic(f);
    let df = dense::derivative(&f);
    let a0 = dense::gcd_field(&f, &df);
    let mut b = quo(&f, &a0);
    let c = quo(&df, &a0);
    let mut d = dense::sub(&c, &dense::derivative(&b));
    let mut i = 1;
    while deg(&b) > 0 {
        let a = dense::gcd_field(&b, &d);
        let nb = quo(&b, &a);
        let nc = quo(&d, &a);
        d = dense::sub(&nc, &dense::derivative(&nb));
        if deg(&a) > 0 {
            out.push((monic(&a), i));
        }
        b = nb;
        i += 1;
    }
    out
}

/// One `(s, h)` pair of `H` in `residue_reduce`.
struct Entry {
    s: Vec<Rat>,
    /// `None`: `h = d`; otherwise the subresultant `S(t, z)`.
    sub: Option<MPoly<RF>>,
    /// `Poly(s, z).monic()` may be over `QQ(...)` with symbols: the stock
    /// `real_roots` then raises `NotImplementedError` for a non-integral
    /// `s` (see `symbolic_first`).
    symbolic: bool,
}

/// `solve([p[i] - c1*q[i] ...], c1)`: the unique solution, `None` when the
/// system has none or does not determine `c1`.
fn solve_ratio(eqs: &[(RF, RF)]) -> Option<RF> {
    let c = eqs.iter().find(|(_, q)| !q.is_zero()).map(|(p, q)| p.div(q))?;
    eqs.iter().all(|(p, q)| p.sub(&c.mul(q)).is_zero()).then_some(c)
}

fn mpoly_lcm(a: &MPoly<Int>, b: &MPoly<Int>) -> MPoly<Int> {
    if a.is_constant() && a.constant_value().is_some_and(|c| c.abs().is_one()) {
        return b.clone();
    }
    let (_, _, bg) = gcd_mz(a, b);
    a.mul(&bg)
}

impl<'a> Tower<'a> {
    /// `levels[0..=top]` (bottom first, level 0 is the base `x` with
    /// `D = 1`); `None` when some `D(t_k)` is not a polynomial in `t_k`.
    pub fn new(levels: Vec<Level>, n: usize, check: &'a mut dyn FnMut() -> bool) -> Option<Self> {
        let mut d_up = Vec::with_capacity(levels.len());
        for l in &levels {
            d_up.push(up_of(&l.d, l.var)?);
        }
        Some(Tower { levels, d_up, n, trace: Vec::new(), check })
    }

    fn tick(&mut self) -> Out<()> {
        if (self.check)() { Ok(()) } else { Err(Stop::Decline) }
    }

    fn active(&self, k: usize) -> Vec<(usize, RF)> {
        self.levels[..k].iter().map(|l| (l.var, l.d.clone())).collect()
    }

    /// `kappa`: the derivation of the levels below `k` on `K_{k-1}`.
    fn kappa(&self, c: &RF, k: usize) -> RF {
        derivation(c, &self.active(k))
    }

    /// The total derivation of `p in K_{k-1}[t_k]`.
    fn deriv(&self, p: &[RF], k: usize) -> UP {
        let below = self.active(k);
        let mut out: UP = p.iter().map(|c| derivation(c, &below)).collect();
        dense::trim(&mut out);
        dense::add(&out, &dense::mul(&dense::derivative(p), &self.d_up[k]))
    }

    fn rf(&self, p: &[RF], k: usize) -> RF {
        up_rf(p, self.levels[k].var, self.n)
    }

    fn lower_vars(&self, k: usize) -> Vec<usize> {
        self.levels[..k].iter().map(|l| l.var).collect()
    }

    /// The stock `PolynomialError` of `gcd(a, s**j, *kkinv)` in the
    /// subresultant normalization: the denominators of `lc_t(S)` (as a
    /// polynomial in `z`) have a multi-term lcm involving a lower monomial.
    fn lc_denominator_raises(&self, lc: &[RF], k: usize) -> bool {
        let mut l = MPoly::one(self.n);
        for c in lc {
            if !c.is_zero() {
                l = mpoly_lcm(&l, &c.d);
            }
        }
        if l.len() <= 1 {
            return false;
        }
        let used = l.used_vars();
        self.lower_vars(k).iter().any(|v| used.contains(v))
    }

    /// Whether the constant first entry `(c, 1)` of `r.sqf_list_include()`
    /// gives a non-unit `S = gcd(c, kappa(c))` (an extra `(1, h)` pair of
    /// `H`), following the domain of `Poly(r, z)`.
    fn constant_entry(&self, r: &[RF], k: usize) -> bool {
        let coeffs: Vec<&RF> = r.iter().filter(|c| !c.is_zero()).collect();
        if coeffs.iter().any(|c| !c.d.is_constant()) {
            return false;
        }
        let mut ring_qq = coeffs.iter().any(|c| !c.d.is_one());
        if coeffs.iter().all(|c| c.n.is_constant()) {
            if ring_qq {
                return false;
            }
            let mut g = Int::zero();
            for c in &coeffs {
                g = g.gcd(&c.n.constant_value().unwrap());
            }
            return !g.abs().is_one();
        }
        let mut c = MPoly::zero(self.n);
        for v in &coeffs {
            c = if c.is_zero() { v.n.clone() } else { gcd_mz_only(&c, &v.n) };
        }
        let kc = self.kappa(&RF::from_poly(c.clone()), k);
        let g = if kc.is_zero() {
            c
        } else {
            if !kc.d.is_constant() {
                return false;
            }
            ring_qq |= !kc.d.is_one();
            gcd_mz_only(&c, &kc.n)
        };
        if ring_qq {
            !g.is_constant()
        } else {
            !(g.is_constant() && g.constant_value().is_some_and(|v| v.abs().is_one()))
        }
    }

    /// The content kept by `S = gcd(p, kappa(p))` of `splitfactor_sqf` for
    /// the first factor `p = content(r)*m` of `r.sqf_list_include()`
    /// (multiplicity 1, `m` primitive), with whether the ring of that `gcd`
    /// has rational coefficients (`QQ[...]`, where the integer part of the
    /// content is a unit). `None` when `Poly(r, z)` or `kappa(p)` has a
    /// symbolic denominator: the `gcd` is then a field gcd, without content.
    fn first_content(&self, r: &[RF], m: &[RF], k: usize) -> Option<(MPoly<Int>, bool)> {
        let coeffs: Vec<&RF> = r.iter().filter(|c| !c.is_zero()).collect();
        if coeffs.iter().any(|c| !c.d.is_constant()) {
            return None;
        }
        let mut ring_qq = coeffs.iter().any(|c| !c.d.is_one());
        let mut c = MPoly::zero(self.n);
        for v in &coeffs {
            c = if c.is_zero() { v.n.clone() } else { gcd_mz_only(&c, &v.n) };
        }
        // The primitive numerator of `m`.
        let mut l = MPoly::one(self.n);
        for v in m.iter().filter(|v| !v.is_zero()) {
            l = mpoly_lcm(&l, &v.d);
        }
        let lf = RF::from_poly(l);
        let nums: Vec<MPoly<Int>> = m.iter().filter(|v| !v.is_zero()).map(|v| v.mul(&lf).n).collect();
        let mut g = MPoly::zero(self.n);
        for v in &nums {
            g = if g.is_zero() { v.clone() } else { gcd_mz_only(&g, v) };
        }
        let cf = RF::from_poly(c.clone());
        let mut kc = MPoly::zero(self.n);
        for v in &nums {
            let (_, m0, _) = gcd_mz(v, &g);
            let dv = self.kappa(&cf.mul(&RF::from_poly(m0)), k);
            if dv.is_zero() {
                continue;
            }
            if !dv.d.is_constant() {
                return None;
            }
            ring_qq |= !dv.d.is_one();
            kc = if kc.is_zero() { dv.n } else { gcd_mz_only(&kc, &dv.n) };
        }
        Some((if kc.is_zero() { c } else { gcd_mz_only(&c, &kc) }, ring_qq))
    }

    /// Whether that kept content is not constant. Then `Poly(S, z).monic()`
    /// is over `QQ(...)`, where `real_roots` raises
    /// `NotImplementedError("sorted roots not supported over QQ")` for a
    /// non-integral monic `s`: such an entry is left to the stock body.
    fn symbolic_first(&self, r: &[RF], m: &[RF], k: usize) -> bool {
        self.first_content(r, m, k).is_some_and(|(kept, _)| !kept.is_constant())
    }

    /// Whether a constant `S` of the first factor is not one (`S.is_one`
    /// false): the stock loop of `residue_reduce` then still processes it.
    fn first_constant_nonunit(&self, r: &[RF], m: &[RF], k: usize) -> bool {
        match self.first_content(r, m, k) {
            None => false,
            Some((kept, true)) => !kept.is_constant(),
            Some((kept, false)) => !(kept.is_constant() && kept.constant_value().is_some_and(|v| v.abs().is_one())),
        }
    }

    /// `residue_reduce(a, d, DE, z)` at level `k` for a reduced `a/d`,
    /// `d` monic: the pairs of `H` and `b`.
    fn residue_reduce(&mut self, a: &[RF], d: &[RF], k: usize) -> Out<(Vec<Entry>, bool)> {
        self.tick()?;
        if a.is_empty() || deg(d) == 0 {
            return Ok((Vec::new(), true));
        }
        let a = dense::rem(a, d);
        let dd = self.deriv(d, k);
        let lift = |p: &[RF], ze: u32| -> Vec<(Vec<u32>, RF)> {
            p.iter().enumerate().filter(|(_, c)| !c.is_zero()).map(|(i, c)| (vec![i as u32, ze], c.clone())).collect()
        };
        let big_d = MPoly::from_terms(2, lift(d, 0));
        let mut q = lift(&a, 0);
        q.extend(lift(&dense::neg(&dd), 1));
        let big_q = MPoly::from_terms(2, q);
        let prs = if deg(&dd) <= deg(d) {
            crate::poly::prs::subresultants(&big_d, &big_q, 0)
        } else {
            crate::poly::prs::subresultants(&big_q, &big_d, 0)
        }
        .map_err(|_| Stop::Decline)?;
        if prs.resultant.degree(0) > 0 {
            return Err(Stop::Decline);
        }
        let r = dense::from_mpoly(&prs.resultant, 1);
        if r.is_empty() {
            return Err(Stop::Decline);
        }
        self.tick()?;
        let factors = sqf_field(&r);
        let kappas: Vec<UP> = factors.iter().map(|(m, _)| {
            let mut v: UP = m.iter().map(|c| self.kappa(c, k)).collect();
            dense::trim(&mut v);
            v
        }).collect();
        let b = kappas.iter().all(|v| v.is_empty());
        let nd = deg(d) as usize;
        let last = |i: usize| prs.sequence.iter().rev().find(|p| p.degree(0) == i as i64).cloned();
        let mut out = Vec::new();
        // The constant entry: `Sp` starts with `(c, 1)` when the smallest
        // multiplicity exceeds 1 (or `r` is constant in `z`).
        let first = factors.first().map(|f| f.1);
        let prepended = first != Some(1);
        let constant = prepended && self.constant_entry(&r, k);
        if constant {
            if nd == 1 {
                if b {
                    out.push(Entry { s: vec![Rat::one()], sub: None, symbolic: false });
                }
            } else if let Some(sub) = last(1) {
                let lc = dense::from_mpoly(&sub.lc_in(0), 1);
                if self.lc_denominator_raises(&lc, k) {
                    return Err(Stop::Decline);
                }
                if b {
                    out.push(Entry { s: vec![Rat::one()], sub: Some(sub), symbolic: false });
                }
            }
        }
        for (j, ((m, i), km)) in factors.iter().zip(&kappas).enumerate() {
            let i = *i;
            let symbolic = j == 0 && i == 1 && self.symbolic_first(&r, m, k);
            let s = if km.is_empty() { m.clone() } else { monic(&dense::gcd_field(m, km)) };
            if deg(&s) < 1 {
                // A constant S_i is one unless it is the content kept for
                // the first factor; a non-unit one still reaches the stock
                // loop (its kkinv gcd), without a root (b is false: kappa
                // of the factor is not zero).
                if i != nd && j == 0 && i == 1 && !prepended && self.first_constant_nonunit(&r, m, k) {
                    if let Some(sub) = last(i) {
                        let lc = dense::from_mpoly(&sub.lc_in(0), 1);
                        if self.lc_denominator_raises(&lc, k) {
                            return Err(Stop::Decline);
                        }
                    }
                }
                continue;
            }
            let Some(sq) = is_rational_poly(&s) else { return Err(Stop::Decline) };
            if i == nd {
                out.push(Entry { s: sq, sub: None, symbolic });
                continue;
            }
            let Some(sub) = last(i) else { continue };
            let lc = dense::from_mpoly(&sub.lc_in(0), 1);
            if self.lc_denominator_raises(&lc, k) || deg(&dense::gcd_field(&lc, &s)) != 0 {
                return Err(Stop::Decline);
            }
            out.push(Entry { s: sq, sub: Some(sub), symbolic });
        }
        Ok((out, b))
    }

    /// `is_log_deriv_k_t_radical_in_field(fa, fd, DE, case)` at level `k`.
    /// The `u` of the answer is the stock value (cancelled at `k >= 1`).
    fn is_log_deriv(&mut self, fa: &[RF], fd: &[RF], k: usize, case: &str) -> Out<Option<(Int, RF)>> {
        self.tick()?;
        if fd.is_empty() {
            return Err(Stop::Decline);
        }
        // fa, fd = fa.cancel(fd, include=True); d monic.
        let g = dense::gcd_field(fa, fd);
        let (a, d) = (quo(fa, &g), quo(fd, &g));
        let l = dense::lc(&d).inv();
        let a = dense::scale(&a, &l);
        let d = dense::scale(&d, &l);
        self.trace.push("z");
        let (h, b) = self.residue_reduce(&a, &d, k)?;
        if !b {
            return Ok(None);
        }
        let mut roots = Vec::with_capacity(h.len());
        let mut ok = true;
        for e in &h {
            self.trace.push("");
            if e.symbolic && e.s.iter().any(|c| !c.d.is_one()) {
                return Err(Stop::Decline);
            }
            match crate::integrate::prde_log::rational_split(&e.s) {
                Some(r) => roots.push(r),
                None => ok = false,
            }
        }
        if !ok {
            return Ok(None);
        }
        self.trace.push("i");
        let n = self.n;
        let mut terms: Vec<(UP, Rat)> = Vec::new();
        for (e, rs) in h.iter().zip(&roots) {
            for rho in rs {
                let hr = match &e.sub {
                    None => d.clone(),
                    Some(sub) => monic(&dense::from_mpoly(&sub.eval(1, &konst(rho, n)), 0)),
                };
                if hr.is_empty() {
                    return Err(Stop::Decline);
                }
                terms.push((hr, rho.clone()));
            }
        }
        let var = self.levels[k].var;
        let mut p = self.rf(&a, k).div(&self.rf(&d, k));
        for (hr, rho) in &terms {
            let dh = self.rf(&self.deriv(hr, k), k);
            p = p.sub(&konst(rho, n).mul(&dh).div(&self.rf(hr, k)));
        }
        if p.d.degree(var) > 0 {
            return Ok(None);
        }
        let pdeg = p.n.degree(var).max(0) as isize;
        if !p.is_zero() && pdeg >= deg(&self.d_up[k]).max(1) {
            return Ok(None);
        }
        let level_case = self.levels[k].case.clone();
        let case = if case == "auto" { level_case.as_str() } else { case };
        let (num, u) = match case {
            "exp" => {
                let dk = &self.d_up[k];
                if level_case != "exp" || k == 0 || deg(dk) != 1 || !dk[0].is_zero() {
                    return Err(Stop::Decline);
                }
                let eta = dk[1].clone();
                let below = self.levels[k - 1].var;
                let (pa, pd) = split(&p, below);
                let (wa, wd) = split(&eta, below);
                let Some((nn, e, v)) = self.heu(&pa, &pd, &wa, &wd, k - 1)? else { return Ok(None) };
                let t = RF::from_poly(MPoly::var(var, n));
                let Some(te) = rf_pow(&t, &e) else { return Err(Stop::Decline) };
                (nn, v.mul(&te))
            }
            "primitive" => {
                // Any level above the base (the stock branch only lowers the
                // level); `p` free of `t` (always so unless `deg D >= 2`).
                if k == 0 || p.n.degree(var) > 0 {
                    return Err(Stop::Decline);
                }
                let below = self.levels[k - 1].var;
                let (pa, pd) = split(&p, below);
                let Some(r) = self.is_log_deriv(&pa, &pd, k - 1, "auto")? else { return Ok(None) };
                r
            }
            "base" => {
                let sqf = deg(&d) < 1 || deg(&dense::gcd_field(&d, &dense::derivative(&d))) == 0;
                if !sqf || deg(&a) >= deg(&d) {
                    return Ok(None);
                }
                if k != 0 {
                    // Any level: the stock `u` is the uncancelled product
                    // `Mul(*[Pow(h, rho*n)])` of Polys, `1` without terms.
                    if terms.is_empty() {
                        return Ok(Some((Int::one(), RF::one())));
                    }
                    return Err(Stop::Decline);
                }
                let mut nn = Int::one();
                for (_, rho) in &terms {
                    nn = nn.lcm(&rho.d).abs();
                }
                let mut u = RF::one();
                for (hr, rho) in &terms {
                    let e = rho.mul(&Rat::from(nn.clone()));
                    let Some(f) = rf_pow(&self.rf(hr, k), &e.n) else { return Err(Stop::Decline) };
                    u = u.mul(&f);
                }
                return Ok(Some((nn, u)));
            }
            "tan" => {
                return Err(Stop::NotImplemented(
                    "The hypertangent case is not yet implemented for is_log_deriv_k_t_radical_in_field()".into(),
                ))
            }
            "other_linear" | "other_nonlinear" => {
                return Err(Stop::Value(format!("The {} case is not supported in this function.", case)))
            }
            other => {
                return Err(Stop::Value(format!(
                    "case must be one of {{'primitive', 'exp', 'tan', 'base', 'auto'}}, not {}",
                    other
                )))
            }
        };
        let mut cd = num.clone();
        for (_, rho) in &terms {
            cd = cd.lcm(&rho.d).abs();
        }
        let Some(m) = cd.div_exact(&num) else {
            return Err(Stop::Value("Inexact division".into()));
        };
        let Some(mut out) = rf_pow(&u, &m) else { return Err(Stop::Decline) };
        for (hr, rho) in &terms {
            let e = rho.mul(&Rat::from(cd.clone()));
            let Some(f) = rf_pow(&self.rf(hr, k), &e.n) else { return Err(Stop::Decline) };
            out = out.mul(&f);
        }
        Ok(Some((cd, out)))
    }

    /// The finishing step of `parametric_log_deriv_heu` for `c1 = c`.
    #[allow(clippy::too_many_arguments)]
    fn heu_finish(&mut self, c: Option<RF>, fa: &[RF], fd: &[RF], wa: &[RF], wd: &[RF], k: usize)
        -> Out<Option<(Int, Int, RF)>> {
        let Some(c) = c.as_ref().and_then(rational) else { return Ok(None) };
        let n = self.n;
        let (mm, nn) = (c.n.clone(), c.d.clone());
        let nfa = dense::sub(
            &dense::scale(&dense::mul(fa, wd), &int_rf(&nn, n)),
            &dense::scale(&dense::mul(wa, fd), &int_rf(&mm, n)),
        );
        let nfd = dense::mul(fd, wd);
        let Some((q, v)) = self.is_log_deriv(&nfa, &nfd, k, "auto")? else { return Ok(None) };
        if v.is_zero() || q.is_zero() {
            return Ok(None);
        }
        Ok(Some((q.mul(&nn), q.mul(&mm), v)))
    }

    /// `splitfactor(p, DE)` at level `k` (`p` of positive degree), up to
    /// constant factors. A constant piece goes through the stock `gcd` over
    /// the generators `1/t_i` (`kinv`); the plain-symbol rational functions
    /// of the tower never hold a factor that `Poly(..., *kinv)` rejects, and
    /// the constant factors of `ln`, `ls` do not change `z` up to a unit,
    /// hence neither the remainders nor the solution `c1` of the caller.
    fn splitfactor(&mut self, p: &[RF], k: usize) -> Out<(UP, UP)> {
        self.tick()?;
        let one = vec![RF::one()];
        if p.is_empty() {
            return Ok((Vec::new(), one));
        }
        if deg(p) == 0 {
            return Ok((p.to_vec(), one));
        }
        let dp = self.deriv(p, k);
        if dp.is_empty() {
            return Ok((p.to_vec(), one));
        }
        let h = dense::gcd_field(p, &dp);
        let g = dense::gcd_field(p, &dense::derivative(p));
        let (s, rem) = dense::divrem(&h, &g);
        if !rem.is_empty() {
            return Err(Stop::Decline);
        }
        if deg(&s) == 0 {
            return Ok((p.to_vec(), one));
        }
        let (qn, qs) = self.splitfactor(&quo(p, &s), k)?;
        Ok((qn, dense::mul(&qs, &s)))
    }

    /// `parametric_log_deriv_heu(fa, fd, wa, wd, DE)` at level `k`.
    fn heu(&mut self, fa: &[RF], fd: &[RF], wa: &[RF], wd: &[RF], k: usize) -> Out<Option<(Int, Int, RF)>> {
        self.tick()?;
        self.trace.push("c1");
        if fd.is_empty() || wd.is_empty() {
            return Err(Stop::Decline);
        }
        let p = quo(fa, fd);
        let q = quo(wa, wd);
        let b = (deg(&self.d_up[k]) - 1).max(0);
        let c = deg(&p).max(deg(&q));
        if deg(&q) > b {
            let eqs: Vec<(RF, RF)> = ((b + 1)..=c)
                .map(|i| (p.get(i as usize).cloned().unwrap_or_else(RF::zero), q.get(i as usize).cloned().unwrap_or_else(RF::zero)))
                .collect();
            let s = solve_ratio(&eqs);
            return self.heu_finish(s, fa, fd, wa, wd, k);
        }
        if deg(&p) > b {
            return Ok(None);
        }
        let l = lcm_up(&monic(fd), &monic(wd));
        if deg(&l) < 1 {
            // `z = ls*gcd(ln, ln')` of a constant `l` is free of `t`.
            return Ok(None);
        }
        let (ln, ls) = self.splitfactor(&l, k)?;
        let z = dense::mul(&ls, &dense::gcd_field(&ln, &dense::derivative(&ln)));
        if deg(&z) < 1 {
            return Ok(None);
        }
        let r1 = dense::rem(&dense::mul(fa, &quo(&l, fd)), &z);
        let r2 = dense::rem(&dense::mul(wa, &quo(&l, wd)), &z);
        let eqs: Vec<(RF, RF)> = (0..deg(&z) as usize)
            .map(|i| (r1.get(i).cloned().unwrap_or_else(RF::zero), r2.get(i).cloned().unwrap_or_else(RF::zero)))
            .collect();
        let s = solve_ratio(&eqs);
        self.heu_finish(s, fa, fd, wa, wd, k)
    }

    fn top(&self) -> usize {
        self.levels.len() - 1
    }

    /// `is_log_deriv_k_t_radical_in_field(fa, fd, DE, case)` at the top
    /// level (`k >= 1`); `fa`, `fd` are polynomials in its variable.
    pub fn log_deriv_top(&mut self, fa: &[RF], fd: &[RF], case: &str) -> Out<Option<(Int, RF)>> {
        let k = self.top();
        if k == 0 {
            return Err(Stop::Decline);
        }
        self.is_log_deriv(fa, fd, k, case)
    }

    /// `parametric_log_deriv_heu(fa, fd, wa, wd, DE)` at the top level
    /// (`k >= 1`).
    pub fn heu_top(&mut self, fa: &[RF], fd: &[RF], wa: &[RF], wd: &[RF]) -> Out<Option<(Int, Int, RF)>> {
        let k = self.top();
        if k == 0 {
            return Err(Stop::Decline);
        }
        self.heu(fa, fd, wa, wd, k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::ratfunc::Scope;

    // Variables: 0 = x, 1 = t0.
    fn var(i: usize) -> RF {
        RF::from_poly(MPoly::var(i, 2))
    }
    fn c(v: i64) -> RF {
        RF::from_poly(MPoly::constant(Int::from(v), 2))
    }
    fn exp_tower(check: &mut dyn FnMut() -> bool) -> Tower<'_> {
        // D = [1, t0]: t0 = exp(x).
        let levels = vec![
            Level { var: 0, d: c(1), case: "base".into() },
            Level { var: 1, d: var(1), case: "exp".into() },
        ];
        Tower::new(levels, 2, check).unwrap()
    }

    fn log_tower(check: &mut dyn FnMut() -> bool) -> Tower<'_> {
        // D = [1, 1/x]: t0 = log(x).
        let levels = vec![
            Level { var: 0, d: c(1), case: "base".into() },
            Level { var: 1, d: var(0).inv(), case: "primitive".into() },
        ];
        Tower::new(levels, 2, check).unwrap()
    }
    fn q(n: i64, d: i64) -> RF {
        c(n).mul(&c(d).inv())
    }

    // The expected values are those of SymPy 1.14 on the same calls.
    #[test]
    fn exp_level_heuristic_gives_up() {
        let _s = Scope::new(2);
        // f = 1 = D(t0)/t0 and f = t0/(t0 + 1): the heuristic of
        // `parametric_log_deriv` returns None (as in the reference).
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        assert_eq!(tw.log_deriv_top(&[c(1)], &[c(1)], "auto").unwrap(), None);
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        assert_eq!(tw.log_deriv_top(&[c(0), c(1)], &[c(1), c(1)], "auto").unwrap(), None);
        // f = 2/t0: no solution, Dummies z, one real_roots Dummy, i.
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        assert_eq!(tw.log_deriv_top(&[c(2)], &[c(0), c(1)], "auto").unwrap(), None);
        assert_eq!(tw.trace, vec!["z", "", "i"]);
    }

    #[test]
    fn log_level_answers() {
        let _s = Scope::new(2);
        // f = -1/(2*x): (2, 1/x).
        let mut ok = || true;
        let mut tw = log_tower(&mut ok);
        let (n, u) = tw.log_deriv_top(&[c(-1)], &[c(2).mul(&var(0))], "auto").unwrap().unwrap();
        assert_eq!(n, Int::from(2));
        assert_eq!(u, var(0).inv());
        // f = 3/(x*(2*t0 - 1)) = (3/2)*D(2*t0 - 1)/(2*t0 - 1):
        // (2, t0**3 - 3*t0**2/2 + 3*t0/4 - 1/8).
        let mut ok = || true;
        let mut tw = log_tower(&mut ok);
        let (n, u) = tw.log_deriv_top(&[c(3)], &[c(0).sub(&var(0)), c(2).mul(&var(0))], "auto").unwrap().unwrap();
        assert_eq!(n, Int::from(2));
        let t = var(1);
        let want = t.mul(&t).mul(&t).sub(&q(3, 2).mul(&t).mul(&t)).add(&q(3, 4).mul(&t)).sub(&q(1, 8));
        assert_eq!(u, want);
    }

    fn tower<'a>(d: RF, case: &str, check: &'a mut dyn FnMut() -> bool) -> Tower<'a> {
        let levels = vec![
            Level { var: 0, d: c(1), case: "base".into() },
            Level { var: 1, d, case: case.into() },
        ];
        Tower::new(levels, 2, check).unwrap()
    }

    #[test]
    fn base_case_above_the_base_level() {
        let _s = Scope::new(2);
        // case='base' on an exp level, f = 0: no residue term, (1, 1).
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        assert_eq!(tw.log_deriv_top(&[], &[c(1)], "base").unwrap(), Some((Int::one(), RF::one())));
        // f = 1/t0 on a log level (D t0 = 1/x), case='base': the residue
        // x is not a constant, SymPy gives None.
        let mut ok = || true;
        let mut tw = log_tower(&mut ok);
        assert_eq!(tw.log_deriv_top(&[c(1)], &[c(0), c(1)], "base"), Ok(None));
        // f = 1/(x*t0): SymPy gives (1, t0), the uncancelled product
        // `Mul(*[Pow(h, rho*n)])`, left to the stock body.
        let mut ok = || true;
        let mut tw = log_tower(&mut ok);
        assert_eq!(tw.log_deriv_top(&[c(1)], &[c(0), var(0)], "base"), Err(Stop::Decline));
    }

    #[test]
    fn heuristic_constant_lcm() {
        let _s = Scope::new(2);
        // fd = -2*x, wd = 3 free of t0: z = ls*gcd(ln, ln') has no t0.
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        assert_eq!(tw.heu_top(&[], &[c(-2).mul(&var(0))], &[], &[c(3)]).unwrap(), None);
        assert_eq!(tw.trace, vec!["c1"]);
    }

    #[test]
    fn kept_content_of_the_residue_gcd() {
        let _s = Scope::new(2);
        // D = x*t0 (t0 = exp(x**2/2)): r = -3x/2*z + 3x/4 keeps no content in
        // gcd(r, kappa(r)), the answer is computed: n = 6.
        let x = var(0);
        let t = var(1);
        let mut ok = || true;
        let mut tw = tower(x.mul(&t), "exp", &mut ok);
        let fa = [c(-3), c(6).mul(&x).mul(&x).add(&c(2))];
        let fd = [c(-18).mul(&x), c(12).mul(&x)];
        let (n, _) = tw.log_deriv_top(&fa, &fd, "auto").unwrap().unwrap();
        assert_eq!(n, Int::from(6));
        // D = 1 + t0**2 (tan): r = -x**2*z + 3x**2/2 keeps x, the stock
        // real_roots raises over QQ[x]: left to the stock body.
        let mut ok = || true;
        let mut tw = tower(c(1).add(&t.mul(&t)), "tan", &mut ok);
        let fa = [c(0), c(0), c(-3)];
        let fd = [c(2).mul(&x), c(-2)];
        assert_eq!(tw.log_deriv_top(&fa, &fd, "auto"), Err(Stop::Decline));
    }

    #[test]
    fn tan_and_unknown_cases_raise() {
        let _s = Scope::new(2);
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        assert!(matches!(tw.log_deriv_top(&[c(0)], &[c(1)], "tan"), Err(Stop::NotImplemented(_))));
        let mut ok = || true;
        let mut tw = exp_tower(&mut ok);
        let e = tw.log_deriv_top(&[c(0)], &[c(1)], "bogus").unwrap_err();
        assert_eq!(e, Stop::Value("case must be one of {'primitive', 'exp', 'tan', 'base', 'auto'}, not bogus".into()));
    }
}
