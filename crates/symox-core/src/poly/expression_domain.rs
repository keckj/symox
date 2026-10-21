//! The reference expression domain `EX` and its dense polynomials.
//!
//! An element wraps an arbitrary expression. Construction never simplifies;
//! arithmetic follows `ExpressionDomain.Expression` literally: sums and
//! differences short-circuit on a structural zero, products on a structural
//! zero or two numbers, and every other result is normalized by `cancel`
//! followed by a shallow `expand(mul=True)` with all other hints off. Zero
//! detection is `not ex.is_zero`, so an undecided coefficient counts as
//! nonzero, exactly as in the reference.
//!
//! Polynomials use the reference recursive dense representation (`dmp`) and
//! the algorithms below are line-by-line ports of `densebasic`,
//! `densearith`, `densetools`, `euclidtools` and `sqfreetools` for a field
//! without an associated ring: division is `dmp_ff_div`, gcds are
//! subresultant PRS gcds, contents are trivial. The reference errors
//! (`PolynomialDivisionFailed`, `ExactQuotientFailed`, `NotInvertible`,
//! `MultivariatePolynomialError`, `ZeroDivisionError`) are reported through
//! [`PolyError::Reference`].
use super::convert::{PResult, PolyError};
use std::collections::BTreeMap;
use crate::expand::{expand_with_options, ExpandBudget, ExpandOptions};
use crate::expr::{Expr, Func, Kind};

/// An element of the reference `EX` domain.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct EX(pub Expr);

fn reference(kind: &'static str, message: impl Into<String>) -> PolyError {
    PolyError::Reference { kind, message: message.into() }
}

/// `Expression.simplify`: `ex.cancel().expand(**eflags)`.
fn simplify(e: Expr) -> PResult<EX> {
    let cancelled = super::api::cancel(&e)?;
    let options = ExpandOptions {
        deep: false,
        basic: false,
        mul: true,
        multinomial: false,
        power_base: false,
        power_exp: false,
        log: false,
        trig: false,
        func: false,
        force: false,
    };
    let expanded = expand_with_options(&cancelled, &options, &mut ExpandBudget::default())
        .map_err(|error| PolyError::Unsupported(error.0))?;
    Ok(EX(expanded))
}

impl EX {
    pub fn new(e: Expr) -> Self {
        EX(e)
    }
    pub fn zero() -> Self {
        EX(Expr::zero())
    }
    pub fn one() -> Self {
        EX(Expr::one())
    }
    pub fn int(n: i64) -> Self {
        EX(Expr::int(n))
    }
    pub fn as_expr(&self) -> &Expr {
        &self.0
    }
    /// `bool(EX(a))`: `not a.is_zero`; an undecided zero test is truthy.
    pub fn truthy(&self) -> bool {
        crate::assumptions::is_zero(&self.0) != Some(true)
    }
    /// `a == EX.zero`: structural comparison with the integer zero.
    fn is_structural_zero(&self) -> bool {
        self.0 == Expr::zero()
    }
    /// `K.is_one(a)`: structural comparison with the integer one.
    pub fn is_one(&self) -> bool {
        self.0 == Expr::one()
    }
    pub fn add(&self, g: &Self) -> PResult<Self> {
        if g.is_structural_zero() {
            return Ok(self.clone());
        }
        if self.is_structural_zero() {
            return Ok(g.clone());
        }
        simplify(crate::ops::add(&[self.0.clone(), g.0.clone()]))
    }
    pub fn sub(&self, g: &Self) -> PResult<Self> {
        if g.is_structural_zero() {
            return Ok(self.clone());
        }
        if self.is_structural_zero() {
            return Ok(g.neg());
        }
        simplify(crate::ops::sub(&self.0, &g.0))
    }
    pub fn mul(&self, g: &Self) -> PResult<Self> {
        if self.is_structural_zero() || g.is_structural_zero() {
            return Ok(EX::zero());
        }
        if self.0.as_number().is_some() && g.0.as_number().is_some() {
            return Ok(EX(crate::ops::mul(&[self.0.clone(), g.0.clone()])));
        }
        simplify(crate::ops::mul(&[self.0.clone(), g.0.clone()]))
    }
    pub fn neg(&self) -> Self {
        EX(crate::ops::neg(&self.0))
    }
    pub fn pow(&self, n: i64) -> PResult<Self> {
        simplify(crate::ops::pow(&self.0, &Expr::int(n)))
    }
    pub fn pow_expr(&self, n: &Expr) -> PResult<Self> {
        simplify(crate::ops::pow(&self.0, n))
    }
    /// `K.quo`, `K.exquo` and `/` of the field `EX`.
    pub fn div(&self, g: &Self) -> PResult<Self> {
        simplify(crate::ops::div(&self.0, &g.0))
    }
    pub fn abs(&self) -> Self {
        EX(crate::functions::func(Func::Abs, vec![self.0.clone()]))
    }
    /// `EX.is_negative`: `could_extract_minus_sign`.
    pub fn is_negative(&self) -> bool {
        crate::decompose::could_extract_minus_sign(&self.0)
    }
    /// `EX.is_positive`, `is_nonnegative` and `is_nonpositive` read the
    /// numeric factor of `as_coeff_mul`.
    pub fn is_positive(&self) -> bool {
        crate::ops::as_coeff_mul(&self.0).0.is_positive()
    }
    pub fn is_nonnegative(&self) -> bool {
        let c = crate::ops::as_coeff_mul(&self.0).0;
        c.is_positive() || c.is_zero()
    }
    pub fn numer(&self) -> Self {
        EX(super::api::as_numer_denom(&self.0).0)
    }
    pub fn denom(&self) -> Self {
        EX(super::api::as_numer_denom(&self.0).1)
    }
    /// `Expression.lcm`: the expression-level `lcm`.
    pub fn lcm(&self, g: &Self) -> PResult<Self> {
        Ok(EX(super::api::lcm(&self.0, &g.0, &[])?))
    }
    fn key(&self) -> crate::sortkey::Key {
        crate::sortkey::sort_key(&self.0)
    }
}

impl std::fmt::Display for EX {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "EX({})", crate::printer::to_str(&self.0))
    }
}

/// A recursive dense polynomial: `L` at level 0, `N` above.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Dmp {
    L(Vec<EX>),
    N(Vec<Dmp>),
}

impl Dmp {
    fn l(&self) -> &Vec<EX> {
        match self {
            Dmp::L(v) => v,
            Dmp::N(_) => unreachable!("dense level mismatch"),
        }
    }
    fn n(&self) -> &Vec<Dmp> {
        match self {
            Dmp::N(v) => v,
            Dmp::L(_) => unreachable!("dense level mismatch"),
        }
    }
    fn len(&self) -> usize {
        match self {
            Dmp::L(v) => v.len(),
            Dmp::N(v) => v.len(),
        }
    }
}

impl std::fmt::Display for Dmp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Dmp::L(v) => write!(f, "[{}]", v.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")),
            Dmp::N(v) => write!(f, "[{}]", v.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")),
        }
    }
}

fn division_failed(f: &Dmp, g: &Dmp) -> PolyError {
    reference("PolynomialDivisionFailed", format!("{f}\u{1f}{g}"))
}
fn zero_division() -> PolyError {
    reference("ZeroDivisionError", "polynomial division")
}

// ----------------------------------------------------------------------
// densebasic
// ----------------------------------------------------------------------

pub fn dmp_zero(u: usize) -> Dmp {
    let mut r = Dmp::L(Vec::new());
    for _ in 0..u {
        r = Dmp::N(vec![r]);
    }
    r
}
pub fn dmp_zero_p(f: &Dmp, u: usize) -> bool {
    let mut f = f;
    for _ in 0..u {
        let v = f.n();
        if v.len() != 1 {
            return false;
        }
        f = &v[0];
    }
    f.l().is_empty()
}
pub fn dmp_ground(c: EX, u: usize) -> Dmp {
    if !c.truthy() {
        return dmp_zero(u);
    }
    let mut r = Dmp::L(vec![c]);
    for _ in 0..u {
        r = Dmp::N(vec![r]);
    }
    r
}
pub fn dmp_one(u: usize) -> Dmp {
    dmp_ground(EX::one(), u)
}
fn dmp_ground_p(f: &Dmp, c: Option<&EX>, u: usize) -> bool {
    if let Some(c) = c {
        if !c.truthy() {
            return dmp_zero_p(f, u);
        }
    }
    let mut f = f;
    for _ in 0..u {
        let v = f.n();
        if v.len() != 1 {
            return false;
        }
        f = &v[0];
    }
    match c {
        None => f.l().len() <= 1,
        Some(c) => f.l().len() == 1 && &f.l()[0] == c,
    }
}
pub fn dmp_one_p(f: &Dmp, u: usize) -> bool {
    dmp_ground_p(f, Some(&EX::one()), u)
}
/// `dmp_degree`, with `-1` for the reference `-oo`.
pub fn dmp_degree(f: &Dmp, u: usize) -> i64 {
    if dmp_zero_p(f, u) {
        -1
    } else {
        f.len() as i64 - 1
    }
}
fn dup_degree(f: &[EX]) -> i64 {
    f.len() as i64 - 1
}
fn dup_strip(f: Vec<EX>) -> Vec<EX> {
    if f.is_empty() || f[0].truthy() {
        return f;
    }
    let i = f.iter().position(EX::truthy).unwrap_or(f.len());
    f[i..].to_vec()
}
fn dmp_strip(f: Vec<Dmp>, u: usize) -> Dmp {
    let whole = Dmp::N(f);
    if dmp_zero_p(&whole, u) {
        return whole;
    }
    let Dmp::N(f) = whole else { unreachable!() };
    let i = f.iter().position(|c| !dmp_zero_p(c, u - 1)).unwrap_or(f.len());
    if i == f.len() {
        dmp_zero(u)
    } else {
        Dmp::N(f[i..].to_vec())
    }
}
fn strip(f: Dmp, u: usize) -> Dmp {
    match f {
        Dmp::L(v) => Dmp::L(dup_strip(v)),
        Dmp::N(v) => dmp_strip(v, u),
    }
}
fn dup_lc(f: &[EX]) -> EX {
    f.first().cloned().unwrap_or_else(EX::zero)
}
fn dup_tc(f: &[EX]) -> EX {
    f.last().cloned().unwrap_or_else(EX::zero)
}
pub fn dmp_ground_lc(f: &Dmp, u: usize) -> EX {
    let mut f = f;
    for _ in 0..u {
        f = &f.n()[0];
    }
    dup_lc(f.l())
}
pub fn dmp_ground_tc(f: &Dmp, u: usize) -> EX {
    let mut f = f;
    for _ in 0..u {
        f = f.n().last().expect("nonempty level");
    }
    dup_tc(f.l())
}
/// `dmp_LC` above level zero.
fn dmp_lc(f: &Dmp) -> Dmp {
    f.n()[0].clone()
}
fn dmp_zeros(n: usize, u: usize) -> Vec<Dmp> {
    (0..n).map(|_| dmp_zero(u)).collect()
}
pub fn dmp_degree_in(f: &Dmp, j: usize, u: usize) -> i64 {
    fn rec(g: &Dmp, v: usize, i: usize, j: usize) -> i64 {
        if i == j {
            return dmp_degree(g, v);
        }
        g.n().iter().map(|c| rec(c, v - 1, i + 1, j)).max().unwrap_or(-1)
    }
    if j == 0 {
        return dmp_degree(f, u);
    }
    rec(f, u, 0, j)
}
pub fn dmp_degree_list(f: &Dmp, u: usize) -> Vec<i64> {
    fn rec(g: &Dmp, v: usize, i: usize, degs: &mut [i64]) {
        degs[i] = degs[i].max(dmp_degree(g, v));
        if v > 0 {
            for c in g.n() {
                rec(c, v - 1, i + 1, degs);
            }
        }
    }
    let mut degs = vec![-1; u + 1];
    rec(f, u, 0, &mut degs);
    degs
}

/// `dmp_to_dict`: monomials with truthy coefficients only.
pub fn dmp_to_dict(f: &Dmp, u: usize) -> Vec<(Vec<u32>, EX)> {
    let mut out = Vec::new();
    fn rec(f: &Dmp, u: usize, prefix: &mut Vec<u32>, out: &mut Vec<(Vec<u32>, EX)>) {
        if u == 0 {
            let v = f.l();
            let n = v.len();
            for k in 0..n {
                let c = &v[n - 1 - k];
                if c.truthy() {
                    let mut m = prefix.clone();
                    m.push(k as u32);
                    out.push((m, c.clone()));
                }
            }
            return;
        }
        let n = dmp_degree(f, u);
        for k in 0..=n.max(-1) {
            if n < 0 {
                break;
            }
            prefix.push(k as u32);
            rec(&f.n()[(n - k) as usize], u - 1, prefix, out);
            prefix.pop();
        }
    }
    rec(f, u, &mut Vec::new(), &mut out);
    out
}

/// `dmp_from_dict`: absent monomials are zero; interior coefficients keep
/// their value even when falsy, as in the reference.
pub fn dmp_from_dict(f: &[(Vec<u32>, EX)], u: usize) -> Dmp {
    if u == 0 {
        if f.is_empty() {
            return Dmp::L(Vec::new());
        }
        let n = f.iter().map(|(m, _)| m[0]).max().unwrap();
        let mut h = Vec::with_capacity(n as usize + 1);
        for k in (0..=n).rev() {
            h.push(f.iter().rev().find(|(m, _)| m[0] == k).map(|(_, c)| c.clone()).unwrap_or_else(EX::zero));
        }
        return Dmp::L(dup_strip(h));
    }
    if f.is_empty() {
        return dmp_zero(u);
    }
    let n = f.iter().map(|(m, _)| m[0]).max().unwrap();
    let mut h = Vec::with_capacity(n as usize + 1);
    for k in (0..=n).rev() {
        let mut tail: Vec<(Vec<u32>, EX)> = Vec::new();
        for (m, c) in f.iter().filter(|(m, _)| m[0] == k) {
            let t = m[1..].to_vec();
            if let Some(slot) = tail.iter_mut().find(|(s, _)| *s == t) {
                slot.1 = c.clone();
            } else {
                tail.push((t, c.clone()));
            }
        }
        if tail.is_empty() {
            h.push(dmp_zero(u - 1));
        } else {
            h.push(dmp_from_dict(&tail, u - 1));
        }
    }
    dmp_strip(h, u)
}

fn igcd(a: u32, b: u32) -> u32 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

pub fn dmp_multi_deflate(polys: &[Dmp], u: usize) -> (Vec<u32>, Vec<Dmp>) {
    if u == 0 {
        let mut big = 0;
        for p in polys {
            let v = p.l();
            if dup_degree(v) <= 0 {
                return (vec![1], polys.to_vec());
            }
            let mut g = 0;
            for i in 0..v.len() {
                if !v[v.len() - 1 - i].truthy() {
                    continue;
                }
                g = igcd(g, i as u32);
                if g == 1 {
                    return (vec![1], polys.to_vec());
                }
            }
            big = igcd(big, g);
        }
        let out = polys.iter().map(|p| Dmp::L(p.l().iter().step_by(big as usize).cloned().collect())).collect();
        return (vec![big], out);
    }
    let mut dicts = Vec::new();
    let mut b = vec![0u32; u + 1];
    for p in polys {
        let f = dmp_to_dict(p, u);
        if !dmp_zero_p(p, u) {
            for (m, _) in &f {
                for (i, e) in m.iter().enumerate() {
                    b[i] = igcd(b[i], *e);
                }
            }
        }
        dicts.push(f);
    }
    for x in b.iter_mut() {
        if *x == 0 {
            *x = 1;
        }
    }
    if b.iter().all(|&x| x == 1) {
        return (b, polys.to_vec());
    }
    let out = dicts
        .into_iter()
        .map(|f| {
            let h: Vec<_> = f.into_iter().map(|(m, c)| (m.iter().zip(&b).map(|(a, d)| a / d).collect(), c)).collect();
            dmp_from_dict(&h, u)
        })
        .collect();
    (b, out)
}

fn dup_inflate(f: &[EX], m: u32) -> Vec<EX> {
    if m == 1 || f.is_empty() {
        return f.to_vec();
    }
    let mut result = vec![f[0].clone()];
    for c in &f[1..] {
        for _ in 0..m - 1 {
            result.push(EX::zero());
        }
        result.push(c.clone());
    }
    result
}
pub fn dmp_inflate(f: &Dmp, m: &[u32], u: usize) -> Dmp {
    fn rec(g: &Dmp, m: &[u32], v: usize, i: usize) -> Dmp {
        if v == 0 {
            return Dmp::L(dup_inflate(g.l(), m[i]));
        }
        let g: Vec<Dmp> = g.n().iter().map(|c| rec(c, m, v - 1, i + 1)).collect();
        let mut result = vec![g[0].clone()];
        for c in &g[1..] {
            for _ in 1..m[i] {
                result.push(dmp_zero(v - 1));
            }
            result.push(c.clone());
        }
        Dmp::N(result)
    }
    if u == 0 {
        return Dmp::L(dup_inflate(f.l(), m[0]));
    }
    if m.iter().all(|&x| x == 1) {
        return f.clone();
    }
    rec(f, m, u, 0)
}

/// `dmp_terms_gcd`.
pub fn dmp_terms_gcd(f: &Dmp, u: usize) -> (Vec<u32>, Dmp) {
    if dmp_ground_tc(f, u).truthy() || dmp_zero_p(f, u) {
        return (vec![0; u + 1], f.clone());
    }
    let dict = dmp_to_dict(f, u);
    let mut g = dict[0].0.clone();
    for (m, _) in &dict[1..] {
        for (a, b) in g.iter_mut().zip(m) {
            *a = (*a).min(*b);
        }
    }
    if g.iter().all(|&x| x == 0) {
        return (g, f.clone());
    }
    let h: Vec<_> = dict.into_iter().map(|(m, c)| (m.iter().zip(&g).map(|(a, b)| a - b).collect(), c)).collect();
    (g.clone(), dmp_from_dict(&h, u))
}

// ----------------------------------------------------------------------
// densearith
// ----------------------------------------------------------------------

fn dup_add_term(f: &[EX], c: &EX, i: usize) -> PResult<Vec<EX>> {
    if !c.truthy() {
        return Ok(f.to_vec());
    }
    let n = f.len();
    if n > 0 && i == n - 1 {
        let mut h = vec![f[0].add(c)?];
        h.extend_from_slice(&f[1..]);
        return Ok(dup_strip(h));
    }
    if i >= n {
        let mut h = vec![c.clone()];
        h.extend((0..i - n).map(|_| EX::zero()));
        h.extend_from_slice(f);
        return Ok(h);
    }
    let m = n - i - 1;
    let mut h = f[..m].to_vec();
    h.push(f[m].add(c)?);
    h.extend_from_slice(&f[m + 1..]);
    Ok(h)
}
fn dmp_add_term(f: &Dmp, c: &Dmp, i: usize, u: usize) -> PResult<Dmp> {
    if u == 0 {
        let c = c.l().first().cloned().unwrap_or_else(EX::zero);
        return Ok(Dmp::L(dup_add_term(f.l(), &c, i)?));
    }
    let v = u - 1;
    if dmp_zero_p(c, v) {
        return Ok(f.clone());
    }
    let fv = f.n();
    let n = fv.len();
    if i == n - 1 {
        let mut h = vec![dmp_add(&fv[0], c, v)?];
        h.extend_from_slice(&fv[1..]);
        return Ok(dmp_strip(h, u));
    }
    if i >= n {
        let mut h = vec![c.clone()];
        h.extend(dmp_zeros(i - n, v));
        h.extend_from_slice(fv);
        return Ok(Dmp::N(h));
    }
    let m = n - i - 1;
    let mut h = fv[..m].to_vec();
    h.push(dmp_add(&fv[m], c, v)?);
    h.extend_from_slice(&fv[m + 1..]);
    Ok(Dmp::N(h))
}
fn dup_mul_term(f: &[EX], c: &EX, i: usize) -> PResult<Vec<EX>> {
    if !c.truthy() || f.is_empty() {
        return Ok(Vec::new());
    }
    let mut h = f.iter().map(|cf| cf.mul(c)).collect::<PResult<Vec<_>>>()?;
    h.extend((0..i).map(|_| EX::zero()));
    Ok(h)
}
/// `dmp_mul_term` with a coefficient one level down.
fn dmp_mul_term(f: &Dmp, c: &Dmp, i: usize, u: usize) -> PResult<Dmp> {
    if u == 0 {
        let c = dup_lc(c.l());
        return Ok(Dmp::L(dup_mul_term(f.l(), &c, i)?));
    }
    let v = u - 1;
    if dmp_zero_p(f, u) {
        return Ok(f.clone());
    }
    if dmp_zero_p(c, v) {
        return Ok(dmp_zero(u));
    }
    let mut h = f.n().iter().map(|cf| dmp_mul(cf, c, v)).collect::<PResult<Vec<_>>>()?;
    h.extend(dmp_zeros(i, v));
    Ok(Dmp::N(h))
}
fn dup_mul_ground(f: &[EX], c: &EX) -> PResult<Vec<EX>> {
    if !c.truthy() || f.is_empty() {
        return Ok(Vec::new());
    }
    f.iter().map(|cf| cf.mul(c)).collect()
}
pub fn dmp_mul_ground(f: &Dmp, c: &EX, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_mul_ground(f.l(), c)?));
    }
    Ok(Dmp::N(f.n().iter().map(|cf| dmp_mul_ground(cf, c, u - 1)).collect::<PResult<_>>()?))
}
fn dup_quo_ground(f: &[EX], c: &EX) -> PResult<Vec<EX>> {
    if !c.truthy() {
        return Err(zero_division());
    }
    if f.is_empty() {
        return Ok(Vec::new());
    }
    f.iter().map(|cf| cf.div(c)).collect()
}
pub fn dmp_quo_ground(f: &Dmp, c: &EX, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_quo_ground(f.l(), c)?));
    }
    Ok(Dmp::N(f.n().iter().map(|cf| dmp_quo_ground(cf, c, u - 1)).collect::<PResult<_>>()?))
}
/// `dmp_exquo_ground`: for a field, the same quotients as `quo_ground`.
pub fn dmp_exquo_ground(f: &Dmp, c: &EX, u: usize) -> PResult<Dmp> {
    dmp_quo_ground(f, c, u)
}
fn dup_neg(f: &[EX]) -> Vec<EX> {
    f.iter().map(EX::neg).collect()
}
pub fn dmp_neg(f: &Dmp, u: usize) -> Dmp {
    match f {
        Dmp::L(v) => Dmp::L(dup_neg(v)),
        Dmp::N(v) => Dmp::N(v.iter().map(|c| dmp_neg(c, u - 1)).collect()),
    }
}
fn dup_add(f: &[EX], g: &[EX]) -> PResult<Vec<EX>> {
    if f.is_empty() {
        return Ok(g.to_vec());
    }
    if g.is_empty() {
        return Ok(f.to_vec());
    }
    let (df, dg) = (f.len(), g.len());
    if df == dg {
        return Ok(dup_strip(f.iter().zip(g).map(|(a, b)| a.add(b)).collect::<PResult<_>>()?));
    }
    let k = df.abs_diff(dg);
    let (mut h, f, g) = if df > dg { (f[..k].to_vec(), &f[k..], g) } else { (g[..k].to_vec(), f, &g[k..]) };
    for (a, b) in f.iter().zip(g) {
        h.push(a.add(b)?);
    }
    Ok(h)
}
pub fn dmp_add(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_add(f.l(), g.l())?));
    }
    let df = dmp_degree(f, u);
    if df < 0 {
        return Ok(g.clone());
    }
    let dg = dmp_degree(g, u);
    if dg < 0 {
        return Ok(f.clone());
    }
    let v = u - 1;
    let (fv, gv) = (f.n(), g.n());
    if df == dg {
        return Ok(dmp_strip(fv.iter().zip(gv).map(|(a, b)| dmp_add(a, b, v)).collect::<PResult<_>>()?, u));
    }
    let k = (df - dg).unsigned_abs() as usize;
    let (mut h, f, g) = if df > dg { (fv[..k].to_vec(), &fv[k..], &gv[..]) } else { (gv[..k].to_vec(), &fv[..], &gv[k..]) };
    for (a, b) in f.iter().zip(g) {
        h.push(dmp_add(a, b, v)?);
    }
    Ok(Dmp::N(h))
}
fn dup_sub(f: &[EX], g: &[EX]) -> PResult<Vec<EX>> {
    if f.is_empty() {
        return Ok(dup_neg(g));
    }
    if g.is_empty() {
        return Ok(f.to_vec());
    }
    let (df, dg) = (f.len(), g.len());
    if df == dg {
        return Ok(dup_strip(f.iter().zip(g).map(|(a, b)| a.sub(b)).collect::<PResult<_>>()?));
    }
    let k = df.abs_diff(dg);
    let (mut h, f, g) = if df > dg { (f[..k].to_vec(), &f[k..], g) } else { (dup_neg(&g[..k]), f, &g[k..]) };
    for (a, b) in f.iter().zip(g) {
        h.push(a.sub(b)?);
    }
    Ok(h)
}
pub fn dmp_sub(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_sub(f.l(), g.l())?));
    }
    let df = dmp_degree(f, u);
    if df < 0 {
        return Ok(dmp_neg(g, u));
    }
    let dg = dmp_degree(g, u);
    if dg < 0 {
        return Ok(f.clone());
    }
    let v = u - 1;
    let (fv, gv) = (f.n(), g.n());
    if df == dg {
        return Ok(dmp_strip(fv.iter().zip(gv).map(|(a, b)| dmp_sub(a, b, v)).collect::<PResult<_>>()?, u));
    }
    let k = (df - dg).unsigned_abs() as usize;
    let (mut h, f, g) = if df > dg {
        (fv[..k].to_vec(), &fv[k..], &gv[..])
    } else {
        (gv[..k].iter().map(|c| dmp_neg(c, v)).collect(), &fv[..], &gv[k..])
    };
    for (a, b) in f.iter().zip(g) {
        h.push(dmp_sub(a, b, v)?);
    }
    Ok(Dmp::N(h))
}
fn dup_slice(f: &[EX], m: usize, n: usize) -> Vec<EX> {
    let k = f.len();
    let big_m = k.saturating_sub(m);
    let big_n = k.saturating_sub(n);
    let mut f: Vec<EX> = if big_n < big_m { f[big_n..big_m].to_vec() } else { Vec::new() };
    while !f.is_empty() && f[0] == EX::zero() {
        f.remove(0);
    }
    if f.is_empty() {
        return f;
    }
    f.extend((0..m).map(|_| EX::zero()));
    f
}
fn dup_rshift(f: &[EX], n: usize) -> Vec<EX> {
    f[..f.len().saturating_sub(n)].to_vec()
}
fn dup_lshift(f: &[EX], n: usize) -> Vec<EX> {
    if f.is_empty() {
        return Vec::new();
    }
    let mut h = f.to_vec();
    h.extend((0..n).map(|_| EX::zero()));
    h
}
fn dup_mul(f: &[EX], g: &[EX]) -> PResult<Vec<EX>> {
    if f == g {
        return dup_sqr(f);
    }
    if f.is_empty() || g.is_empty() {
        return Ok(Vec::new());
    }
    let (df, dg) = (f.len() - 1, g.len() - 1);
    let n = df.max(dg) + 1;
    if n < 100 {
        let mut h = Vec::with_capacity(df + dg + 1);
        for i in 0..=df + dg {
            let mut coeff = EX::zero();
            for j in i.saturating_sub(dg)..=df.min(i) {
                coeff = coeff.add(&f[j].mul(&g[i - j])?)?;
            }
            h.push(coeff);
        }
        return Ok(dup_strip(h));
    }
    let n2 = n / 2;
    let (fl, gl) = (dup_slice(f, 0, n2), dup_slice(g, 0, n2));
    let fh = dup_rshift(&dup_slice(f, n2, n), n2);
    let gh = dup_rshift(&dup_slice(g, n2, n), n2);
    let (lo, hi) = (dup_mul(&fl, &gl)?, dup_mul(&fh, &gh)?);
    let mid = dup_mul(&dup_add(&fl, &fh)?, &dup_add(&gl, &gh)?)?;
    let mid = dup_sub(&mid, &dup_add(&lo, &hi)?)?;
    dup_add(&dup_add(&lo, &dup_lshift(&mid, n2))?, &dup_lshift(&hi, 2 * n2))
}
pub fn dmp_mul(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_mul(f.l(), g.l())?));
    }
    if f == g {
        return dmp_sqr(f, u);
    }
    let df = dmp_degree(f, u);
    if df < 0 {
        return Ok(f.clone());
    }
    let dg = dmp_degree(g, u);
    if dg < 0 {
        return Ok(g.clone());
    }
    let (df, dg) = (df as usize, dg as usize);
    let v = u - 1;
    let (fv, gv) = (f.n(), g.n());
    let mut h = Vec::with_capacity(df + dg + 1);
    for i in 0..=df + dg {
        let mut coeff = dmp_zero(v);
        for j in i.saturating_sub(dg)..=df.min(i) {
            coeff = dmp_add(&coeff, &dmp_mul(&fv[j], &gv[i - j], v)?, v)?;
        }
        h.push(coeff);
    }
    Ok(dmp_strip(h, u))
}
fn dup_sqr(f: &[EX]) -> PResult<Vec<EX>> {
    let df = f.len() as i64 - 1;
    let mut h = Vec::new();
    for i in 0..=(2 * df).max(-1) {
        if df < 0 {
            break;
        }
        let mut c = EX::zero();
        let jmin = (i - df).max(0);
        let jmax = i.min(df);
        let n = jmax - jmin + 1;
        let jmax = jmin + n / 2 - 1;
        for j in jmin..=jmax {
            c = c.add(&f[j as usize].mul(&f[(i - j) as usize])?)?;
        }
        c = c.add(&c)?;
        if n & 1 == 1 {
            let elem = &f[(jmax + 1) as usize];
            c = c.add(&elem.pow(2)?)?;
        }
        h.push(c);
    }
    Ok(dup_strip(h))
}
pub fn dmp_sqr(f: &Dmp, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_sqr(f.l())?));
    }
    let df = dmp_degree(f, u);
    if df < 0 {
        return Ok(f.clone());
    }
    let v = u - 1;
    let fv = f.n();
    let mut h = Vec::new();
    for i in 0..=2 * df {
        let mut c = dmp_zero(v);
        let jmin = (i - df).max(0);
        let jmax = i.min(df);
        let n = jmax - jmin + 1;
        let jmax = jmin + n / 2 - 1;
        for j in jmin..=jmax {
            c = dmp_add(&c, &dmp_mul(&fv[j as usize], &fv[(i - j) as usize], v)?, v)?;
        }
        c = dmp_mul_ground(&c, &EX::int(2), v)?;
        if n & 1 == 1 {
            let elem = dmp_sqr(&fv[(jmax + 1) as usize], v)?;
            c = dmp_add(&c, &elem, v)?;
        }
        h.push(c);
    }
    Ok(dmp_strip(h, u))
}
pub fn dmp_pow(f: &Dmp, n: u32, u: usize) -> PResult<Dmp> {
    if n == 0 {
        return Ok(dmp_one(u));
    }
    if n == 1 || dmp_zero_p(f, u) || (if u == 0 { f.l().len() == 1 && f.l()[0] == EX::one() } else { dmp_one_p(f, u) }) {
        return Ok(f.clone());
    }
    let mut g = dmp_one(u);
    let mut f = f.clone();
    let mut n = n;
    loop {
        let m = n;
        n /= 2;
        if m % 2 == 1 {
            g = dmp_mul(&g, &f, u)?;
            if n == 0 {
                break;
            }
        }
        f = dmp_sqr(&f, u)?;
    }
    Ok(g)
}

/// `dmp_pdiv` / `dmp_prem`; the quotient is skipped for `prem`.
fn pseudo(f: &Dmp, g: &Dmp, u: usize, quotient: bool) -> PResult<(Dmp, Dmp)> {
    let df = dmp_degree(f, u);
    let dg = dmp_degree(g, u);
    if dg < 0 {
        return Err(zero_division());
    }
    let mut q = dmp_zero(u);
    let mut r = f.clone();
    let mut dr = df;
    if df < dg {
        return Ok((q, r));
    }
    let mut big_n = df - dg + 1;
    if u == 0 {
        let lc_g = dup_lc(g.l());
        loop {
            let lc_r = dup_lc(r.l());
            let j = (dr - dg) as usize;
            big_n -= 1;
            if quotient {
                let big_q = dup_mul_ground(q.l(), &lc_g)?;
                q = Dmp::L(dup_add_term(&big_q, &lc_r, j)?);
            }
            let big_r = dup_mul_ground(r.l(), &lc_g)?;
            let big_g = dup_mul_term(g.l(), &lc_r, j)?;
            r = Dmp::L(dup_sub(&big_r, &big_g)?);
            let previous = dr;
            dr = dmp_degree(&r, 0);
            if dr < dg {
                break;
            }
            if dr >= previous {
                return Err(division_failed(f, g));
            }
        }
        let c = lc_g.pow(big_n)?;
        if quotient {
            q = Dmp::L(dup_mul_ground(q.l(), &c)?);
        }
        r = Dmp::L(dup_mul_ground(r.l(), &c)?);
        return Ok((q, r));
    }
    let lc_g = dmp_lc(g);
    loop {
        let lc_r = dmp_lc(&r);
        let j = (dr - dg) as usize;
        big_n -= 1;
        if quotient {
            let big_q = dmp_mul_term(&q, &lc_g, 0, u)?;
            q = dmp_add_term(&big_q, &lc_r, j, u)?;
        }
        let big_r = dmp_mul_term(&r, &lc_g, 0, u)?;
        let big_g = dmp_mul_term(g, &lc_r, j, u)?;
        r = dmp_sub(&big_r, &big_g, u)?;
        let previous = dr;
        dr = dmp_degree(&r, u);
        if dr < dg {
            break;
        }
        if dr >= previous {
            return Err(division_failed(f, g));
        }
    }
    let c = dmp_pow(&lc_g, big_n as u32, u - 1)?;
    if quotient {
        q = dmp_mul_term(&q, &c, 0, u)?;
    }
    r = dmp_mul_term(&r, &c, 0, u)?;
    Ok((q, r))
}
pub fn dmp_pdiv(f: &Dmp, g: &Dmp, u: usize) -> PResult<(Dmp, Dmp)> {
    if u == 0 && g.l().is_empty() {
        return Err(zero_division());
    }
    pseudo(f, g, u, true)
}
pub fn dmp_prem(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    if u == 0 && g.l().is_empty() {
        return Err(zero_division());
    }
    Ok(pseudo(f, g, u, false)?.1)
}

/// `dmp_ff_div`.
pub fn dmp_div(f: &Dmp, g: &Dmp, u: usize) -> PResult<(Dmp, Dmp)> {
    let df = dmp_degree(f, u);
    let dg = dmp_degree(g, u);
    if u == 0 {
        if g.l().is_empty() {
            return Err(zero_division());
        }
        let mut q: Vec<EX> = Vec::new();
        let mut r = f.l().clone();
        let mut dr = df;
        if df < dg {
            return Ok((Dmp::L(q), Dmp::L(r)));
        }
        let lc_g = dup_lc(g.l());
        loop {
            let lc_r = dup_lc(&r);
            let c = lc_r.div(&lc_g)?;
            let j = (dr - dg) as usize;
            q = dup_add_term(&q, &c, j)?;
            let h = dup_mul_term(g.l(), &c, j)?;
            r = dup_sub(&r, &h)?;
            let previous = dr;
            dr = dup_degree(&r);
            if dr < dg {
                break;
            }
            if dr >= previous {
                return Err(division_failed(f, g));
            }
        }
        return Ok((Dmp::L(q), Dmp::L(r)));
    }
    if dg < 0 {
        return Err(zero_division());
    }
    let mut q = dmp_zero(u);
    let mut r = f.clone();
    let mut dr = df;
    if df < dg {
        return Ok((q, r));
    }
    let lc_g = dmp_lc(g);
    let v = u - 1;
    loop {
        let lc_r = dmp_lc(&r);
        let (c, rest) = dmp_div(&lc_r, &lc_g, v)?;
        if !dmp_zero_p(&rest, v) {
            break;
        }
        let j = (dr - dg) as usize;
        q = dmp_add_term(&q, &c, j, u)?;
        let h = dmp_mul_term(g, &c, j, u)?;
        r = dmp_sub(&r, &h, u)?;
        let previous = dr;
        dr = dmp_degree(&r, u);
        if dr < dg {
            break;
        }
        if dr >= previous {
            return Err(division_failed(f, g));
        }
    }
    Ok((q, r))
}
pub fn dmp_rem(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    Ok(dmp_div(f, g, u)?.1)
}
pub fn dmp_quo(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    Ok(dmp_div(f, g, u)?.0)
}
pub fn dmp_exquo(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    let (q, r) = dmp_div(f, g, u)?;
    if dmp_zero_p(&r, u) {
        Ok(q)
    } else {
        Err(reference("ExactQuotientFailed", format!("{g} does not divide {f}")))
    }
}
fn max_by_key(values: Vec<EX>) -> EX {
    let mut best: Option<(crate::sortkey::Key, EX)> = None;
    for v in values {
        let k = v.key();
        // Python's max keeps the first of equal items: replace only when
        // the current maximum compares strictly below the candidate.
        if best.as_ref().is_none_or(|(bk, _)| *bk < k) {
            best = Some((k, v));
        }
    }
    best.map(|(_, v)| v).unwrap_or_else(EX::zero)
}
pub fn dmp_max_norm(f: &Dmp, u: usize) -> EX {
    if u == 0 {
        let v = f.l();
        if v.is_empty() {
            return EX::zero();
        }
        return max_by_key(v.iter().map(EX::abs).collect());
    }
    max_by_key(f.n().iter().map(|c| dmp_max_norm(c, u - 1)).collect())
}
/// Python's `sum` starting from the integer zero: the first addition is a
/// reflected `0 + EX`, which always simplifies.
fn python_sum(values: Vec<EX>) -> PResult<EX> {
    let mut it = values.into_iter();
    let Some(first) = it.next() else { return Ok(EX::zero()) };
    let mut total = simplify(first.0)?;
    for v in it {
        total = total.add(&v)?;
    }
    Ok(total)
}
pub fn dmp_l1_norm(f: &Dmp, u: usize) -> PResult<EX> {
    if u == 0 {
        let v = f.l();
        if v.is_empty() {
            return Ok(EX::zero());
        }
        return python_sum(v.iter().map(EX::abs).collect());
    }
    python_sum(f.n().iter().map(|c| dmp_l1_norm(c, u - 1)).collect::<PResult<_>>()?)
}

// ----------------------------------------------------------------------
// densetools
// ----------------------------------------------------------------------

fn falling(n: i64, m: usize) -> i64 {
    let mut k = n;
    let mut i = n - 1;
    while i > n - m as i64 {
        k *= i;
        i -= 1;
    }
    k
}
fn dup_diff(f: &[EX], m: usize) -> PResult<Vec<EX>> {
    if m == 0 {
        return Ok(f.to_vec());
    }
    let mut n = dup_degree(f);
    if n < m as i64 {
        return Ok(Vec::new());
    }
    let mut deriv = Vec::new();
    for coeff in &f[..f.len() - m] {
        let k = if m == 1 { n } else { falling(n, m) };
        deriv.push(EX::int(k).mul(coeff)?);
        n -= 1;
    }
    Ok(dup_strip(deriv))
}
pub fn dmp_diff(f: &Dmp, m: usize, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_diff(f.l(), m)?));
    }
    if m == 0 {
        return Ok(f.clone());
    }
    let mut n = dmp_degree(f, u);
    if n < m as i64 {
        return Ok(dmp_zero(u));
    }
    let v = u - 1;
    let fv = f.n();
    let mut deriv = Vec::new();
    for coeff in &fv[..fv.len() - m] {
        let k = if m == 1 { n } else { falling(n, m) };
        deriv.push(dmp_mul_ground(coeff, &EX::int(k), v)?);
        n -= 1;
    }
    Ok(dmp_strip(deriv, u))
}
pub fn dmp_diff_in(f: &Dmp, m: usize, j: usize, u: usize) -> PResult<Dmp> {
    fn rec(g: &Dmp, m: usize, v: usize, i: usize, j: usize) -> PResult<Dmp> {
        if i == j {
            return dmp_diff(g, m, v);
        }
        let h = g.n().iter().map(|c| rec(c, m, v - 1, i + 1, j)).collect::<PResult<Vec<_>>>()?;
        Ok(dmp_strip(h, v))
    }
    rec(f, m, u, 0, j)
}
fn dup_integrate(f: &[EX], m: usize) -> PResult<Vec<EX>> {
    if m == 0 || f.is_empty() {
        return Ok(f.to_vec());
    }
    let mut g: Vec<EX> = (0..m).map(|_| EX::zero()).collect();
    for (i, c) in f.iter().rev().enumerate() {
        let mut n = i as i64 + 1;
        for j in 1..m as i64 {
            n *= i as i64 + j + 1;
        }
        g.insert(0, c.div(&EX::int(n))?);
    }
    Ok(g)
}
pub fn dmp_integrate(f: &Dmp, m: usize, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(dup_integrate(f.l(), m)?));
    }
    if m == 0 || dmp_zero_p(f, u) {
        return Ok(f.clone());
    }
    let v = u - 1;
    let mut g = dmp_zeros(m, v);
    for (i, c) in f.n().iter().rev().enumerate() {
        let mut n = i as i64 + 1;
        for j in 1..m as i64 {
            n *= i as i64 + j + 1;
        }
        g.insert(0, dmp_quo_ground(c, &EX::int(n), v)?);
    }
    Ok(Dmp::N(g))
}
pub fn dmp_integrate_in(f: &Dmp, m: usize, j: usize, u: usize) -> PResult<Dmp> {
    fn rec(g: &Dmp, m: usize, v: usize, i: usize, j: usize) -> PResult<Dmp> {
        if i == j {
            return dmp_integrate(g, m, v);
        }
        let h = g.n().iter().map(|c| rec(c, m, v - 1, i + 1, j)).collect::<PResult<Vec<_>>>()?;
        Ok(dmp_strip(h, v))
    }
    rec(f, m, u, 0, j)
}
fn dup_eval(f: &[EX], a: &EX) -> PResult<EX> {
    if !a.truthy() {
        return Ok(dup_tc(f));
    }
    let mut result = EX::zero();
    for c in f {
        result = result.mul(a)?;
        result = result.add(c)?;
    }
    Ok(result)
}
/// `dmp_eval` at the main variable; the result is one level down (a
/// ground element is returned as a level-zero constant list).
fn dmp_eval(f: &Dmp, a: &EX, u: usize) -> PResult<Dmp> {
    if u == 0 {
        return Ok(Dmp::L(vec![dup_eval(f.l(), a)?]));
    }
    if !a.truthy() {
        return Ok(f.n().last().cloned().unwrap_or_else(|| dmp_zero(u - 1)));
    }
    let v = u - 1;
    let fv = f.n();
    let mut result = fv[0].clone();
    for coeff in &fv[1..] {
        result = dmp_mul_ground(&result, a, v)?;
        result = dmp_add(&result, coeff, v)?;
    }
    Ok(result)
}
/// `dup_eval` result for a univariate polynomial.
pub fn dup_eval_ground(f: &Dmp, a: &EX) -> PResult<EX> {
    dup_eval(f.l(), a)
}
/// `dmp_eval_in` for `u > 0`.
pub fn dmp_eval_in(f: &Dmp, a: &EX, j: usize, u: usize) -> PResult<Dmp> {
    fn rec(g: &Dmp, a: &EX, v: usize, i: usize, j: usize) -> PResult<Dmp> {
        if i == j {
            return dmp_eval(g, a, v);
        }
        let h = g.n().iter().map(|c| rec(c, a, v - 1, i + 1, j)).collect::<PResult<Vec<_>>>()?;
        Ok(strip(rebuild(h, v - 1), v - 1))
    }
    rec(f, a, u, 0, j)
}
/// Children of an `N` level after a recursive evaluation removed one level
/// below: level-zero constants become ground lists.
fn rebuild(children: Vec<Dmp>, v: usize) -> Dmp {
    if v == 0 {
        // Each child is a level-zero singleton holding the value.
        Dmp::L(children.into_iter().map(|c| c.l().first().cloned().unwrap_or_else(EX::zero)).collect())
    } else {
        Dmp::N(children)
    }
}
pub fn dmp_ground_monic(f: &Dmp, u: usize) -> PResult<Dmp> {
    if dmp_zero_p(f, u) {
        return Ok(f.clone());
    }
    let lc = dmp_ground_lc(f, u);
    if lc.is_one() {
        return Ok(f.clone());
    }
    dmp_exquo_ground(f, &lc, u)
}
/// `dmp_ground_content` over `EX`: `K.gcd` is one.
pub fn dmp_ground_content(f: &Dmp, u: usize) -> EX {
    if dmp_zero_p(f, u) {
        EX::zero()
    } else {
        EX::one()
    }
}
pub fn dmp_ground_primitive(f: &Dmp, u: usize) -> PResult<(EX, Dmp)> {
    if dmp_zero_p(f, u) {
        return Ok((EX::zero(), f.clone()));
    }
    let cont = dmp_ground_content(f, u);
    if cont.is_one() {
        return Ok((cont, f.clone()));
    }
    Ok((cont.clone(), dmp_quo_ground(f, &cont, u)?))
}
pub fn dmp_compose(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    if u == 0 {
        let (fv, gv) = (f.l(), g.l());
        if gv.len() <= 1 {
            return Ok(Dmp::L(dup_strip(vec![dup_eval(fv, &dup_lc(gv))?])));
        }
        if fv.is_empty() {
            return Ok(Dmp::L(Vec::new()));
        }
        let mut h = vec![fv[0].clone()];
        for c in &fv[1..] {
            h = dup_mul(&h, gv)?;
            h = dup_add_term(&h, c, 0)?;
        }
        return Ok(Dmp::L(h));
    }
    if dmp_zero_p(f, u) {
        return Ok(f.clone());
    }
    let fv = f.n();
    let mut h = Dmp::N(vec![fv[0].clone()]);
    for c in &fv[1..] {
        h = dmp_mul(&h, g, u)?;
        h = dmp_add_term(&h, c, 0, u)?;
    }
    Ok(h)
}
pub fn dup_shift(f: &Dmp, a: &EX) -> PResult<Dmp> {
    let mut f = f.l().clone();
    let n = f.len() as i64 - 1;
    let mut i = n;
    while i > 0 {
        for j in 0..i as usize {
            let t = a.mul(&f[j])?;
            f[j + 1] = f[j + 1].add(&t)?;
        }
        i -= 1;
    }
    Ok(Dmp::L(f))
}
/// `dmp_clear_denoms` for `EX` (no associated ring: `K1 = EX`).
pub fn dmp_clear_denoms(f: &Dmp, u: usize) -> PResult<(EX, Dmp)> {
    fn rec(g: &Dmp, v: usize) -> PResult<EX> {
        let mut common = EX::one();
        if v == 0 {
            for c in g.l() {
                common = common.lcm(&c.denom())?;
            }
        } else {
            for c in g.n() {
                common = common.lcm(&rec(c, v - 1)?)?;
            }
        }
        Ok(common)
    }
    if u == 0 {
        let fv = f.l();
        let mut common = EX::one();
        for c in fv {
            common = common.lcm(&c.denom())?;
        }
        if common.is_one() {
            return Ok((common, f.clone()));
        }
        let h = fv.iter().map(|c| c.numer().mul(&common.div(&c.denom())?)).collect::<PResult<Vec<_>>>()?;
        return Ok((common, Dmp::L(dup_strip(h))));
    }
    let common = rec(f, u)?;
    let f = if common.is_one() { f.clone() } else { dmp_mul_ground(f, &common, u)? };
    Ok((common, f))
}

// ----------------------------------------------------------------------
// euclidtools
// ----------------------------------------------------------------------

fn dup_sub_mul(f: &[EX], g: &[EX], h: &[EX]) -> PResult<Vec<EX>> {
    dup_sub(f, &dup_mul(g, h)?)
}
fn dup_monic(f: &[EX]) -> PResult<Vec<EX>> {
    if f.is_empty() {
        return Ok(Vec::new());
    }
    let lc = dup_lc(f);
    if lc.is_one() {
        return Ok(f.to_vec());
    }
    dup_quo_ground(f, &lc)
}
pub fn dup_half_gcdex(f: &Dmp, g: &Dmp) -> PResult<(Dmp, Dmp)> {
    let mut f = f.clone();
    let mut g = g.clone();
    let mut a = vec![EX::one()];
    let mut b: Vec<EX> = Vec::new();
    while !g.l().is_empty() {
        let (q, r) = dmp_div(&f, &g, 0)?;
        f = g;
        g = r;
        let next = dup_sub_mul(&a, q.l(), &b)?;
        a = b;
        b = next;
    }
    let a = dup_quo_ground(&a, &dup_lc(f.l()))?;
    let f = dup_monic(f.l())?;
    Ok((Dmp::L(a), Dmp::L(f)))
}
pub fn dup_gcdex(f: &Dmp, g: &Dmp) -> PResult<(Dmp, Dmp, Dmp)> {
    let (s, h) = dup_half_gcdex(f, g)?;
    let big_f = dup_sub_mul(h.l(), s.l(), f.l())?;
    let t = dmp_quo(&Dmp::L(big_f), g, 0)?;
    Ok((s, t, h))
}
pub fn dup_invert(f: &Dmp, g: &Dmp) -> PResult<Dmp> {
    let (s, h) = dup_half_gcdex(f, g)?;
    if h.l().len() == 1 && h.l()[0] == EX::one() {
        dmp_rem(&s, g, 0)
    } else {
        Err(PolyError::NotInvertible("zero divisor".into()))
    }
}
/// `dmp_inner_subresultants`: the PRS and the subresultant coefficients.
pub fn dmp_inner_subresultants(f: &Dmp, g: &Dmp, u: usize) -> PResult<(Vec<Dmp>, Vec<Dmp>)> {
    let (mut f, mut g) = (f.clone(), g.clone());
    let mut n = dmp_degree(&f, u);
    let mut m = dmp_degree(&g, u);
    if n < m {
        std::mem::swap(&mut f, &mut g);
        std::mem::swap(&mut n, &mut m);
    }
    if dmp_zero_p(&f, u) {
        return Ok((Vec::new(), Vec::new()));
    }
    if u == 0 {
        // Ground coefficients are stored as level-zero singletons.
        if g.l().is_empty() {
            return Ok((vec![f], vec![Dmp::L(vec![EX::one()])]));
        }
        let mut r = vec![f.clone(), g.clone()];
        let mut d = n - m;
        let b = EX::int(-1).pow(d + 1)?;
        let mut h = dmp_prem(&f, &g, 0)?;
        h = Dmp::L(dup_mul_ground(h.l(), &b)?);
        let mut lc = dup_lc(g.l());
        let mut c = lc.pow(d)?;
        let mut s = vec![EX::one(), c.clone()];
        c = c.neg();
        while !h.l().is_empty() {
            let k = dmp_degree(&h, 0);
            r.push(h.clone());
            f = g;
            g = h;
            d = m - k;
            m = k;
            let b = lc.neg().mul(&c.pow(d)?)?;
            h = dmp_prem(&f, &g, 0)?;
            h = Dmp::L(dup_quo_ground(h.l(), &b)?);
            lc = dup_lc(g.l());
            if d > 1 {
                let q = c.pow(d - 1)?;
                c = lc.neg().pow(d)?.div(&q)?;
            } else {
                c = lc.neg();
            }
            s.push(c.neg());
        }
        return Ok((r, s.into_iter().map(|x| Dmp::L(vec![x])).collect()));
    }
    let v = u - 1;
    if dmp_zero_p(&g, u) {
        return Ok((vec![f], vec![dmp_ground(EX::one(), v)]));
    }
    let mut r = vec![f.clone(), g.clone()];
    let mut d = n - m;
    let b = dmp_pow(&dmp_ground(EX::int(-1), v), (d + 1) as u32, v)?;
    let mut h = dmp_prem(&f, &g, u)?;
    h = dmp_mul_term(&h, &b, 0, u)?;
    let mut lc = dmp_lc(&g);
    let mut c = dmp_pow(&lc, d as u32, v)?;
    let mut s = vec![dmp_ground(EX::one(), v), c.clone()];
    c = dmp_neg(&c, v);
    while !dmp_zero_p(&h, u) {
        let k = dmp_degree(&h, u);
        r.push(h.clone());
        f = g;
        g = h;
        d = m - k;
        m = k;
        let b = dmp_mul(&dmp_neg(&lc, v), &dmp_pow(&c, d as u32, v)?, v)?;
        h = dmp_prem(&f, &g, u)?;
        h = Dmp::N(h.n().iter().map(|ch| dmp_quo(ch, &b, v)).collect::<PResult<_>>()?);
        lc = dmp_lc(&g);
        if d > 1 {
            let p = dmp_pow(&dmp_neg(&lc, v), d as u32, v)?;
            let q = dmp_pow(&c, (d - 1) as u32, v)?;
            c = dmp_quo(&p, &q, v)?;
        } else {
            c = dmp_neg(&lc, v);
        }
        s.push(dmp_neg(&c, v));
    }
    Ok((r, s))
}
pub fn dmp_subresultants(f: &Dmp, g: &Dmp, u: usize) -> PResult<Vec<Dmp>> {
    Ok(dmp_inner_subresultants(f, g, u)?.0)
}
/// `dmp_prs_resultant`; at level zero the resultant is a singleton list.
pub fn dmp_prs_resultant(f: &Dmp, g: &Dmp, u: usize) -> PResult<(Dmp, Vec<Dmp>)> {
    if dmp_zero_p(f, u) || dmp_zero_p(g, u) {
        return Ok((if u == 0 { Dmp::L(vec![EX::zero()]) } else { dmp_zero(u - 1) }, Vec::new()));
    }
    let (r, s) = dmp_inner_subresultants(f, g, u)?;
    if dmp_degree(r.last().unwrap(), u) > 0 {
        return Ok((if u == 0 { Dmp::L(vec![EX::zero()]) } else { dmp_zero(u - 1) }, r));
    }
    Ok((s.last().unwrap().clone(), r))
}
pub fn dmp_resultant(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    Ok(dmp_prs_resultant(f, g, u)?.0)
}
/// `dmp_discriminant`; a level-zero result is a singleton list.
pub fn dmp_discriminant(f: &Dmp, u: usize) -> PResult<Dmp> {
    let d = dmp_degree(f, u);
    let s = if ((d * (d - 1)) / 2) % 2 == 0 { 1 } else { -1 };
    if u == 0 {
        if d <= 0 {
            return Ok(Dmp::L(vec![EX::zero()]));
        }
        let c = dup_lc(f.l());
        let r = dmp_resultant(f, &dmp_diff(f, 1, 0)?, 0)?;
        let r = dup_lc(r.l());
        return Ok(Dmp::L(vec![r.div(&c.mul(&EX::int(s))?)?]));
    }
    let v = u - 1;
    if d <= 0 {
        return Ok(dmp_zero(v));
    }
    let c = dmp_lc(f);
    let r = dmp_resultant(f, &dmp_diff(f, 1, u)?, u)?;
    let c = dmp_mul_ground(&c, &EX::int(s), v)?;
    dmp_quo(&r, &c, v)
}
fn dmp_content(f: &Dmp, u: usize) -> PResult<Dmp> {
    let v = u - 1;
    let fv = f.n();
    let mut cont = fv[0].clone();
    if dmp_zero_p(f, u) {
        return Ok(cont);
    }
    for c in &fv[1..] {
        cont = dmp_gcd(&cont, c, v)?;
        if dmp_one_p(&cont, v) {
            break;
        }
    }
    if dmp_ground_lc(&cont, v).is_negative() {
        Ok(dmp_neg(&cont, v))
    } else {
        Ok(cont)
    }
}
pub fn dmp_primitive(f: &Dmp, u: usize) -> PResult<(Dmp, Dmp)> {
    let cont = dmp_content(f, u)?;
    let v = u - 1;
    if dmp_zero_p(f, u) || dmp_one_p(&cont, v) {
        return Ok((cont, f.clone()));
    }
    Ok((cont.clone(), Dmp::N(f.n().iter().map(|c| dmp_quo(c, &cont, v)).collect::<PResult<_>>()?)))
}
type Triple = (Dmp, Dmp, Dmp);
fn dmp_ff_trivial_gcd(f: &Dmp, g: &Dmp, u: usize) -> PResult<Option<Triple>> {
    let zero_f = dmp_zero_p(f, u);
    let zero_g = dmp_zero_p(g, u);
    if u == 0 {
        if zero_f && zero_g {
            return Ok(Some((dmp_zero(0), dmp_zero(0), dmp_zero(0))));
        }
        if zero_f {
            return Ok(Some((Dmp::L(dup_monic(g.l())?), dmp_zero(0), Dmp::L(vec![dup_lc(g.l())]))));
        }
        if zero_g {
            return Ok(Some((Dmp::L(dup_monic(f.l())?), Dmp::L(vec![dup_lc(f.l())]), dmp_zero(0))));
        }
        return Ok(None);
    }
    if zero_f && zero_g {
        return Ok(Some((dmp_zero(u), dmp_zero(u), dmp_zero(u))));
    }
    if zero_f {
        return Ok(Some((dmp_ground_monic(g, u)?, dmp_zero(u), dmp_ground(dmp_ground_lc(g, u), u))));
    }
    if zero_g {
        return Ok(Some((dmp_ground_monic(f, u)?, dmp_ground(dmp_ground_lc(f, u), u), dmp_zero(u))));
    }
    dmp_simplify_gcd(f, g, u)
}
fn dmp_simplify_gcd(f: &Dmp, g: &Dmp, u: usize) -> PResult<Option<Triple>> {
    let df = dmp_degree(f, u);
    let dg = dmp_degree(g, u);
    if df > 0 && dg > 0 {
        return Ok(None);
    }
    let (big_f, big_g) = if df == 0 && dg == 0 {
        (dmp_lc(f), dmp_lc(g))
    } else if df == 0 {
        (dmp_lc(f), dmp_content(g, u)?)
    } else {
        (dmp_content(f, u)?, dmp_lc(g))
    };
    let v = u - 1;
    let h = dmp_gcd(&big_f, &big_g, v)?;
    let cff = Dmp::N(f.n().iter().map(|cf| dmp_quo(cf, &h, v)).collect::<PResult<_>>()?);
    let cfg = Dmp::N(g.n().iter().map(|cg| dmp_quo(cg, &h, v)).collect::<PResult<_>>()?);
    Ok(Some((Dmp::N(vec![h]), cff, cfg)))
}
fn dmp_ff_prs_gcd(f: &Dmp, g: &Dmp, u: usize) -> PResult<Triple> {
    if let Some(result) = dmp_ff_trivial_gcd(f, g, u)? {
        return Ok(result);
    }
    if u == 0 {
        let h = dmp_subresultants(f, g, 0)?.pop().expect("nonzero inputs");
        let h = Dmp::L(dup_monic(h.l())?);
        let cff = dmp_quo(f, &h, 0)?;
        let cfg = dmp_quo(g, &h, 0)?;
        return Ok((h, cff, cfg));
    }
    let (fc, big_f) = dmp_primitive(f, u)?;
    let (gc, big_g) = dmp_primitive(g, u)?;
    let h = dmp_subresultants(&big_f, &big_g, u)?.pop().expect("nonzero inputs");
    let (c, _, _) = dmp_ff_prs_gcd(&fc, &gc, u - 1)?;
    let (_, h) = dmp_primitive(&h, u)?;
    let h = dmp_mul_term(&h, &c, 0, u)?;
    let h = dmp_ground_monic(&h, u)?;
    let cff = dmp_quo(f, &h, u)?;
    let cfg = dmp_quo(g, &h, u)?;
    Ok((h, cff, cfg))
}
pub fn dmp_inner_gcd(f: &Dmp, g: &Dmp, u: usize) -> PResult<Triple> {
    if u == 0 {
        return dmp_ff_prs_gcd(f, g, 0);
    }
    let (j, polys) = dmp_multi_deflate(&[f.clone(), g.clone()], u);
    let (h, cff, cfg) = dmp_ff_prs_gcd(&polys[0], &polys[1], u)?;
    Ok((dmp_inflate(&h, &j, u), dmp_inflate(&cff, &j, u), dmp_inflate(&cfg, &j, u)))
}
pub fn dmp_gcd(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    Ok(dmp_inner_gcd(f, g, u)?.0)
}
pub fn dmp_lcm(f: &Dmp, g: &Dmp, u: usize) -> PResult<Dmp> {
    let h = dmp_quo(&dmp_mul(f, g, u)?, &dmp_gcd(f, g, u)?, u)?;
    dmp_ground_monic(&h, u)
}
/// `dmp_cancel` for a field without an associated ring.
pub fn dmp_cancel(f: &Dmp, g: &Dmp, u: usize) -> PResult<(EX, EX, Dmp, Dmp)> {
    let (mut cp, cq) = (EX::one(), EX::one());
    let (_, mut p, mut q) = dmp_inner_gcd(f, g, u)?;
    let p_neg = dmp_ground_lc(&p, u).is_negative();
    let q_neg = dmp_ground_lc(&q, u).is_negative();
    if p_neg && q_neg {
        p = dmp_neg(&p, u);
        q = dmp_neg(&q, u);
    } else if p_neg {
        cp = cp.neg();
        p = dmp_neg(&p, u);
    } else if q_neg {
        cp = cp.neg();
        q = dmp_neg(&q, u);
    }
    Ok((cp, cq, p, q))
}

// ----------------------------------------------------------------------
// sqfreetools
// ----------------------------------------------------------------------

pub fn dmp_sqf_p(f: &Dmp, u: usize) -> PResult<bool> {
    if dmp_zero_p(f, u) {
        return Ok(true);
    }
    if u == 0 {
        return Ok(dmp_degree(&dmp_gcd(f, &dmp_diff(f, 1, 0)?, 0)?, 0) == 0);
    }
    for i in 0..=u {
        let fp = dmp_diff_in(f, 1, i, u)?;
        if dmp_zero_p(&fp, u) {
            continue;
        }
        let gcd = dmp_gcd(f, &fp, u)?;
        if dmp_degree_in(&gcd, i, u) != 0 {
            return Ok(false);
        }
    }
    Ok(true)
}
pub fn dmp_sqf_part(f: &Dmp, u: usize) -> PResult<Dmp> {
    if dmp_zero_p(f, u) {
        return Ok(f.clone());
    }
    let f = if dmp_ground_lc(f, u).is_negative() { dmp_neg(f, u) } else { f.clone() };
    let gcd = if u == 0 {
        dmp_gcd(&f, &dmp_diff(&f, 1, 0)?, 0)?
    } else {
        let mut gcd = f.clone();
        for i in 0..=u {
            gcd = dmp_gcd(&gcd, &dmp_diff_in(&f, 1, i, u)?, u)?;
        }
        gcd
    };
    let sqf = dmp_quo(&f, &gcd, u)?;
    dmp_ground_monic(&sqf, u)
}
fn check_degrees(f: &Dmp, u: usize, result: &[(Dmp, usize)]) -> PResult<()> {
    let mut degs = vec![0i64; u + 1];
    for (fac, k) in result {
        for (d, e) in degs.iter_mut().zip(dmp_degree_list(fac, u)) {
            *d += *k as i64 * e;
        }
    }
    let expected = dmp_degree_list(f, u);
    let same = if u == 0 { degs[0] == expected[0] } else { degs == expected };
    if same {
        Ok(())
    } else {
        Err(reference("AssertionError", ""))
    }
}
pub fn dmp_sqf_list(f: &Dmp, u: usize, all: bool) -> PResult<(EX, Vec<(Dmp, usize)>)> {
    let f_orig = f.clone();
    let coeff = dmp_ground_lc(f, u);
    let f = dmp_ground_monic(f, u)?;
    if u == 0 {
        if dmp_degree(&f, 0) <= 0 {
            return Ok((coeff, Vec::new()));
        }
        let mut result = Vec::new();
        let mut i = 1;
        let h = dmp_diff(&f, 1, 0)?;
        let (_, mut p, mut q) = dmp_inner_gcd(&f, &h, 0)?;
        loop {
            let d = dmp_diff(&p, 1, 0)?;
            let h = dmp_sub(&q, &d, 0)?;
            if h.l().is_empty() {
                result.push((p, i));
                break;
            }
            let (g, np, nq) = dmp_inner_gcd(&p, &h, 0)?;
            p = np;
            q = nq;
            if all || dmp_degree(&g, 0) > 0 {
                result.push((g, i));
            }
            i += 1;
        }
        check_degrees(&f_orig, 0, &result)?;
        return Ok((coeff, result));
    }
    let deg = dmp_degree(&f, u);
    if deg < 0 {
        return Ok((coeff, Vec::new()));
    }
    let (content, f) = dmp_primitive(&f, u)?;
    let mut result: Vec<(usize, Dmp)> = Vec::new();
    if deg != 0 {
        let h = dmp_diff(&f, 1, u)?;
        let (_, mut p, mut q) = dmp_inner_gcd(&f, &h, u)?;
        let mut i = 1;
        loop {
            let d = dmp_diff(&p, 1, u)?;
            let h = dmp_sub(&q, &d, u)?;
            if dmp_zero_p(&h, u) {
                result.push((i, p));
                break;
            }
            let (g, np, nq) = dmp_inner_gcd(&p, &h, u)?;
            p = np;
            q = nq;
            if all || dmp_degree(&g, u) > 0 {
                result.push((i, g));
            }
            i += 1;
        }
    }
    let (coeff_content, result_content) = dmp_sqf_list(&content, u - 1, all)?;
    let coeff = coeff.mul(&coeff_content)?;
    for (fac, i) in result_content {
        let fac = Dmp::N(vec![fac]);
        if let Some(slot) = result.iter_mut().find(|(k, _)| *k == i) {
            slot.1 = dmp_mul(&slot.1, &fac, u)?;
        } else {
            result.push((i, fac));
        }
    }
    result.sort_by_key(|(k, _)| *k);
    let result: Vec<(Dmp, usize)> = result.into_iter().map(|(k, f)| (f, k)).collect();
    check_degrees(&f_orig, u, &result)?;
    Ok((coeff, result))
}
pub fn dmp_sqf_list_include(f: &Dmp, u: usize, all: bool) -> PResult<Vec<(Dmp, usize)>> {
    let (coeff, factors) = dmp_sqf_list(f, u, all)?;
    if let Some((first, 1)) = factors.first() {
        let g = dmp_mul_ground(first, &coeff, u)?;
        let mut out = vec![(g, 1)];
        out.extend(factors[1..].iter().cloned());
        return Ok(out);
    }
    let g = if u == 0 { Dmp::L(dup_strip(vec![coeff])) } else { dmp_ground(coeff, u) };
    let mut out = vec![(g, 1)];
    out.extend(factors);
    Ok(out)
}
pub fn dup_gff_list(f: &Dmp) -> PResult<Vec<(Dmp, usize)>> {
    if f.l().is_empty() {
        return Err(PolyError::Failed("greatest factorial factorization doesn't exist for a zero polynomial".into()));
    }
    let f = Dmp::L(dup_monic(f.l())?);
    if dmp_degree(&f, 0) == 0 {
        return Ok(Vec::new());
    }
    let mut g = dmp_gcd(&f, &dup_shift(&f, &EX::one())?, 0)?;
    let mut big_h = dup_gff_list(&g)?;
    for entry in big_h.iter_mut() {
        let (h, k) = entry.clone();
        g = dmp_mul(&g, &dup_shift(&h, &EX::int(-(k as i64)))?, 0)?;
        *entry = (h, k + 1);
    }
    let f = dmp_quo(&f, &g, 0)?;
    if dmp_degree(&f, 0) == 0 {
        return Ok(big_h);
    }
    let mut out = vec![(f, 1)];
    out.extend(big_h);
    Ok(out)
}
/// `dup_sturm`.
pub fn dup_sturm(f: &Dmp) -> PResult<Vec<Dmp>> {
    let f = dmp_sqf_part(f, 0)?;
    let mut sturm = vec![f.clone(), dmp_diff(&f, 1, 0)?];
    while !sturm.last().unwrap().l().is_empty() {
        let n = sturm.len();
        let s = dmp_rem(&sturm[n - 2], &sturm[n - 1], 0)?;
        sturm.push(dmp_neg(&s, 0));
    }
    sturm.pop();
    Ok(sturm)
}

// ----------------------------------------------------------------------
// Poly-level wrapper
// ----------------------------------------------------------------------

/// `Poly(..., domain=EX)`: a dense `EX` polynomial and its generators.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExpressionPoly {
    pub rep: Dmp,
    pub gens: Vec<Expr>,
    /// A composite Gaussian ground (`ZZ_I[y]`, `QQ_I(y)`, ...) whose
    /// coefficients this polynomial stores as expressions.
    pub ground: Option<super::gaussian_ground::GaussianGround>,
}

fn dup_from_raw_dict(g: &BTreeMap<usize, EX>) -> Vec<EX> {
    let Some(n) = g.keys().next_back().copied() else { return Vec::new() };
    dup_strip((0..=n).rev().map(|k| g.get(&k).cloned().unwrap_or_else(EX::zero)).collect())
}
fn dup_right_decompose(f: &[EX], s: usize) -> PResult<Vec<EX>> {
    let n = f.len() - 1;
    let lc = f[0].clone();
    let fd: BTreeMap<usize, EX> = f.iter().enumerate().filter(|(_, c)| c.truthy()).map(|(i, c)| (n - i, c.clone())).collect();
    let mut g: BTreeMap<usize, EX> = BTreeMap::new();
    g.insert(s, EX::one());
    let r = n / s;
    for i in 1..s {
        let mut coeff = EX::zero();
        for j in 0..i {
            let Some(fc) = fd.get(&(n + j - i)) else { continue };
            let Some(gc) = g.get(&(s - j)) else { continue };
            let k = EX::int(i as i64 - (r * j) as i64);
            coeff = coeff.add(&k.mul(fc)?.mul(gc)?)?;
        }
        let d = EX::int((i * r) as i64).mul(&lc)?;
        g.insert(s - i, coeff.div(&d)?);
    }
    Ok(dup_from_raw_dict(&g))
}
fn dup_left_decompose(f: &[EX], h: &[EX]) -> PResult<Option<Vec<EX>>> {
    let mut g: BTreeMap<usize, EX> = BTreeMap::new();
    let mut f = Dmp::L(f.to_vec());
    let h = Dmp::L(h.to_vec());
    let mut i = 0;
    while !f.l().is_empty() {
        let (q, r) = dmp_div(&f, &h, 0)?;
        if dmp_degree(&r, 0) > 0 {
            return Ok(None);
        }
        g.insert(i, r.l().first().cloned().unwrap_or_else(EX::zero));
        f = q;
        i += 1;
    }
    Ok(Some(dup_from_raw_dict(&g)))
}
fn dup_decompose_step(f: &[EX]) -> PResult<Option<(Vec<EX>, Vec<EX>)>> {
    let df = f.len().saturating_sub(1);
    for s in 2..df {
        if df % s != 0 {
            continue;
        }
        let h = dup_right_decompose(f, s)?;
        if let Some(g) = dup_left_decompose(f, &h)? {
            return Ok(Some((g, h)));
        }
    }
    Ok(None)
}
pub fn dup_decompose(f: &[EX]) -> PResult<Vec<Vec<EX>>> {
    let mut f = f.to_vec();
    let mut chain = Vec::new();
    while let Some((g, h)) = dup_decompose_step(&f)? {
        f = g;
        chain.insert(0, h);
    }
    chain.insert(0, f);
    Ok(chain)
}

/// `pure_complex`: a number `a + b*I` with numeric `a` and `b`.
fn pure_complex(e: &Expr) -> bool {
    let imaginary = |t: &Expr| t.is_constant(crate::expr::Constant::I);
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(_) => imaginary(e),
        Kind::Mul(_, factors) => factors.len() == 1 && imaginary(&factors[0].0) && factors[0].1.is_one(),
        Kind::Add(_, terms) => terms.len() == 1 && imaginary(&terms[0].0),
        _ => false,
    }
}

/// `parallel_poly_from_expr(exprs, *gens)` when the joint `construct_domain`
/// over every coefficient selects `EX`; `None` keeps the other constructors.
pub fn joint_selection(exprs: &[Expr], gens: &[Expr]) -> Option<Vec<ExpressionPoly>> {
    if gens.is_empty() {
        return None;
    }
    let polys: Vec<ExpressionPoly> = exprs.iter().map(|e| ExpressionPoly::from_expr(e, gens).ok()).collect::<Option<_>>()?;
    let coeffs: Vec<Expr> = polys.iter().flat_map(|p| p.terms().into_iter().map(|(_, c)| c)).collect();
    matches!(construct_selects_expression(&coeffs), Ok(true)).then_some(polys)
}

/// Whether the reference `construct_domain` (no extension, default
/// composite option, no field flag) selects `EX` for these coefficients.
///
/// `_construct_simple` accepts rational, float and pure complex numbers.
/// Otherwise `_construct_composite` reads the generators of the numerators
/// and denominators and gives up (selecting `EX`) when there are none, when
/// one is an algebraic number, or when two share a free symbol.
pub fn construct_selects_expression(coeffs: &[Expr]) -> PResult<bool> {
    if coeffs.iter().all(pure_complex) {
        return Ok(false);
    }
    let mut parts = Vec::with_capacity(2 * coeffs.len());
    for c in coeffs {
        let (n, d) = super::api::as_numer_denom(c);
        parts.push(n);
        parts.push(d);
    }
    let converted = super::convert::from_exprs(&parts)?;
    let gens: Vec<&Expr> = converted
        .gens
        .iter()
        .filter(|g| !g.is_constant(crate::expr::Constant::I))
        .collect();
    if gens.is_empty() {
        return Ok(true);
    }
    let mut seen: Vec<Expr> = Vec::new();
    for g in gens {
        if crate::functions::is_numeric_value(g) && crate::assumptions::fact(g, "algebraic") == Some(true) {
            return Ok(true);
        }
        let symbols = crate::traverse::free_symbols(g);
        if symbols.iter().any(|s| seen.contains(s)) {
            return Ok(true);
        }
        seen.extend(symbols);
    }
    Ok(false)
}

/// The generators of the fraction field the reference `_construct_composite`
/// selects for these coefficients when no domain is requested.
///
/// It is `Some` only when the composite construction applies (no `EX`
/// fallback, exact coefficients) and some coefficient denominator involves
/// a generator, which is when the reference builds `ZZ(...)` instead of the
/// polynomial ring `ZZ[...]`. Generators are read from the numerators and
/// denominators, so `1/pi` and `exp(-a)` both select the generator
/// `pi` or `exp(a)`.
pub fn composite_fraction_parameters(coeffs: &[Expr]) -> PResult<Option<Vec<Expr>>> {
    Ok(composite_parameters(coeffs)?.and_then(|(gens, fractions)| fractions.then_some(gens)))
}

/// The generators of the reference `_construct_composite` domain of these
/// coefficients and whether it is a fraction field; `None` when the
/// reference selects a simple or the expression domain.
pub fn composite_parameters(coeffs: &[Expr]) -> PResult<Option<(Vec<Expr>, bool)>> {
    if coeffs.iter().all(pure_complex) || construct_selects_expression(coeffs)? {
        return Ok(None);
    }
    let mut parts = Vec::with_capacity(2 * coeffs.len());
    let mut denominators = Vec::with_capacity(coeffs.len());
    for c in coeffs {
        let (n, d) = super::api::as_numer_denom(c);
        parts.push(n);
        denominators.push(d);
    }
    let k = parts.len();
    parts.extend(denominators);
    let converted = super::convert::from_exprs(&parts)?;
    if converted.float || converted.gens.is_empty() || converted.gens.iter().any(|g| g.is_constant(crate::expr::Constant::I)) {
        return Ok(None);
    }
    let fractions = converted.polys[k..].iter().any(|d| !d.is_constant());
    Ok(Some((converted.gens, fractions)))
}

/// `_dict_from_expr_if_gens` for one expression: coefficient sums are
/// plain `Add`s of `Mul(*coeff)` products, never simplified.
pub fn dict_from_expr(e: &Expr, gens: &[Expr]) -> PResult<Vec<(Vec<u32>, Expr)>> {
    super::convert::check_number(e)?;
    let ex = super::convert::expand_if_needed(e);
    let mut out: Vec<(Vec<u32>, Expr)> = Vec::new();
    for t in super::convert::raw_terms(&ex) {
        let mut monom = vec![0u32; gens.len()];
        let mut coeff = vec![Expr::num(t.coeff.clone())];
        for f in &t.factors {
            super::convert::check_number(f)?;
            let (b, k) = super::convert::decompose_factor(f);
            if let Some(i) = gens.iter().position(|g| *g == b) {
                monom[i] = k as u32;
                continue;
            }
            if gens.iter().any(|g| crate::traverse::has_structural(f, g)) {
                return Err(PolyError::Generator(format!(
                    "{} contains an element of the set of generators.",
                    crate::printer::to_str(super::convert::generator_clash(&t.factors, gens).unwrap_or(f))
                )));
            }
            coeff.push(f.clone());
        }
        let c = crate::ops::mul(&coeff);
        if let Some(slot) = out.iter_mut().find(|(m, _)| *m == monom) {
            slot.1 = crate::ops::add(&[slot.1.clone(), c]);
        } else {
            out.push((monom, c));
        }
    }
    Ok(out)
}

impl ExpressionPoly {
    pub fn u(&self) -> usize {
        self.gens.len() - 1
    }
    pub fn with_rep(&self, rep: Dmp) -> Self {
        Self { rep, gens: self.gens.clone(), ground: self.ground.clone() }
    }
    /// Coefficients become `EX` elements without simplification.
    pub fn from_terms(terms: Vec<(Vec<u32>, Expr)>, gens: &[Expr]) -> PResult<Self> {
        if gens.is_empty() {
            return Err(PolyError::Failed("can't initialize from 'dict' without generators".into()));
        }
        let dict: Vec<(Vec<u32>, EX)> = terms.into_iter().map(|(m, c)| (m, EX(c))).collect();
        Ok(Self { rep: dmp_from_dict(&dict, gens.len() - 1), gens: gens.to_vec(), ground: None })
    }
    /// `Poly.from_dict` after a monomial map: equal monomials are added in
    /// `EX`, which short-circuits on a structural zero.
    pub fn from_summed_terms(terms: Vec<(Vec<u32>, Expr)>, gens: &[Expr]) -> PResult<Self> {
        let mut acc: BTreeMap<Vec<u32>, EX> = BTreeMap::new();
        for (m, c) in terms {
            let c = EX(c);
            let sum = match acc.remove(&m) {
                Some(a) => a.add(&c)?,
                None => c,
            };
            acc.insert(m, sum);
        }
        Self::from_terms(acc.into_iter().map(|(m, c)| (m, c.0)).collect(), gens)
    }
    pub fn from_expr(e: &Expr, gens: &[Expr]) -> PResult<Self> {
        if gens.is_empty() {
            return Err(PolyError::Failed("can't initialize from 'dict' without generators".into()));
        }
        Self::from_terms(dict_from_expr(e, gens)?, gens)
    }
    pub fn zero(gens: &[Expr]) -> Self {
        Self { rep: dmp_zero(gens.len() - 1), gens: gens.to_vec(), ground: None }
    }
    /// `Poly.terms()`: lexicographically decreasing monomials.
    pub fn terms(&self) -> Vec<(Vec<u32>, Expr)> {
        let mut t: Vec<_> = dmp_to_dict(&self.rep, self.u()).into_iter().map(|(m, c)| (m, c.0)).collect();
        t.sort_by(|a, b| b.0.cmp(&a.0));
        t
    }
    /// `Poly.as_expr()`: `Add(*[Mul(coeff, *gens**monom)])` in dict order.
    pub fn as_expr(&self) -> Expr {
        let terms: Vec<Expr> = dmp_to_dict(&self.rep, self.u())
            .into_iter()
            .map(|(m, c)| {
                let mut f = vec![c.0];
                for (g, &k) in self.gens.iter().zip(&m) {
                    if k != 0 {
                        f.push(crate::ops::pow(g, &Expr::int(k as i64)));
                    }
                }
                crate::ops::mul(&f)
            })
            .collect();
        crate::ops::add(&terms)
    }
    pub fn is_zero(&self) -> bool {
        dmp_zero_p(&self.rep, self.u())
    }
    /// Reorder onto `gens` (a superset), through truthy terms only.
    pub fn reorder(&self, gens: &[Expr]) -> PResult<Self> {
        if self.gens == gens {
            return Ok(self.clone());
        }
        let positions: Vec<usize> = self
            .gens
            .iter()
            .map(|g| gens.iter().position(|h| h == g).ok_or_else(|| PolyError::Failed("generators do not unify".into())))
            .collect::<PResult<_>>()?;
        let dict: Vec<(Vec<u32>, EX)> = dmp_to_dict(&self.rep, self.u())
            .into_iter()
            .map(|(m, c)| {
                let mut n = vec![0; gens.len()];
                for (i, k) in positions.iter().zip(m) {
                    n[*i] = k;
                }
                (n, c)
            })
            .collect();
        Ok(Self { rep: dmp_from_dict(&dict, gens.len() - 1), gens: gens.to_vec(), ground: self.ground.clone() })
    }
    pub fn add(&self, o: &Self) -> PResult<Self> {
        Ok(self.with_rep(dmp_add(&self.rep, &o.rep, self.u())?))
    }
    pub fn sub(&self, o: &Self) -> PResult<Self> {
        Ok(self.with_rep(dmp_sub(&self.rep, &o.rep, self.u())?))
    }
    pub fn mul(&self, o: &Self) -> PResult<Self> {
        Ok(self.with_rep(dmp_mul(&self.rep, &o.rep, self.u())?))
    }
    pub fn neg(&self) -> Self {
        self.with_rep(dmp_neg(&self.rep, self.u()))
    }
    pub fn pow(&self, n: u32) -> PResult<Self> {
        Ok(self.with_rep(dmp_pow(&self.rep, n, self.u())?))
    }
    pub fn div(&self, o: &Self) -> PResult<(Self, Self)> {
        let (q, r) = dmp_div(&self.rep, &o.rep, self.u())?;
        Ok((self.with_rep(q), self.with_rep(r)))
    }
    /// `Poly.exquo`: a failure is re-raised with the `as_expr` forms.
    pub fn exquo(&self, o: &Self) -> PResult<Self> {
        match dmp_exquo(&self.rep, &o.rep, self.u()) {
            Ok(q) => Ok(self.with_rep(q)),
            Err(PolyError::Reference { kind: "ExactQuotientFailed", .. }) => Err(self.exquo_failed(o)),
            Err(e) => Err(e),
        }
    }
    fn exquo_failed(&self, o: &Self) -> PolyError {
        reference("ExactQuotientFailed", format!("{} does not divide {}", o.as_expr(), self.as_expr()))
    }
    pub fn pdiv(&self, o: &Self) -> PResult<(Self, Self)> {
        let (q, r) = dmp_pdiv(&self.rep, &o.rep, self.u())?;
        Ok((self.with_rep(q), self.with_rep(r)))
    }
    pub fn pexquo(&self, o: &Self) -> PResult<Self> {
        let (q, r) = self.pdiv(o)?;
        if r.is_zero() {
            Ok(q)
        } else {
            Err(self.exquo_failed(o))
        }
    }
    pub fn gcd(&self, o: &Self) -> PResult<Self> {
        Ok(self.with_rep(dmp_gcd(&self.rep, &o.rep, self.u())?))
    }
    pub fn lcm(&self, o: &Self) -> PResult<Self> {
        Ok(self.with_rep(dmp_lcm(&self.rep, &o.rep, self.u())?))
    }
    pub fn cofactors(&self, o: &Self) -> PResult<(Self, Self, Self)> {
        let (h, a, b) = dmp_inner_gcd(&self.rep, &o.rep, self.u())?;
        Ok((self.with_rep(h), self.with_rep(a), self.with_rep(b)))
    }
    pub fn cancel(&self, o: &Self) -> PResult<(Expr, Expr, Self, Self)> {
        let (cp, cq, p, q) = dmp_cancel(&self.rep, &o.rep, self.u())?;
        Ok((cp.0, cq.0, self.with_rep(p), self.with_rep(q)))
    }
    /// `PolyElement.cancel` in a sparse ring over `EX`, which has no
    /// associated ring: cofactors, then the canonical unit of `q`.
    pub fn ring_cancel(&self, o: &Self) -> PResult<(Self, Self)> {
        if self.is_zero() {
            let one = Self::from_terms(vec![(vec![0; self.gens.len()], Expr::one())], &self.gens)?;
            return Ok((self.clone(), one));
        }
        let (_, p, q) = self.cofactors(o)?;
        if dmp_ground_lc(&q.rep, q.u()).is_negative() {
            return Ok((p.neg(), q.neg()));
        }
        Ok((p, q))
    }
    fn univariate(&self) -> PResult<()> {
        if self.gens.len() == 1 {
            Ok(())
        } else {
            Err(reference("ValueError", "univariate polynomial expected"))
        }
    }
    pub fn half_gcdex(&self, o: &Self) -> PResult<(Self, Self)> {
        self.univariate()?;
        let (s, h) = dup_half_gcdex(&self.rep, &o.rep)?;
        Ok((self.with_rep(s), self.with_rep(h)))
    }
    pub fn gcdex(&self, o: &Self) -> PResult<(Self, Self, Self)> {
        self.univariate()?;
        let (s, t, h) = dup_gcdex(&self.rep, &o.rep)?;
        Ok((self.with_rep(s), self.with_rep(t), self.with_rep(h)))
    }
    pub fn invert(&self, o: &Self) -> PResult<Self> {
        self.univariate()?;
        Ok(self.with_rep(dup_invert(&self.rep, &o.rep)?))
    }
    pub fn diff(&self, m: usize, j: usize) -> PResult<Self> {
        Ok(self.with_rep(dmp_diff_in(&self.rep, m, j, self.u())?))
    }
    pub fn integrate(&self, m: usize, j: usize) -> PResult<Self> {
        Ok(self.with_rep(dmp_integrate_in(&self.rep, m, j, self.u())?))
    }
    pub fn monic(&self) -> PResult<Self> {
        Ok(self.with_rep(dmp_ground_monic(&self.rep, self.u())?))
    }
    pub fn lc(&self) -> Expr {
        dmp_ground_lc(&self.rep, self.u()).0
    }
    pub fn tc(&self) -> Expr {
        dmp_ground_tc(&self.rep, self.u()).0
    }
    pub fn content(&self) -> Expr {
        dmp_ground_content(&self.rep, self.u()).0
    }
    pub fn primitive(&self) -> PResult<(Expr, Self)> {
        let (c, p) = dmp_ground_primitive(&self.rep, self.u())?;
        Ok((c.0, self.with_rep(p)))
    }
    pub fn max_norm(&self) -> Expr {
        dmp_max_norm(&self.rep, self.u()).0
    }
    pub fn l1_norm(&self) -> PResult<Expr> {
        Ok(dmp_l1_norm(&self.rep, self.u())?.0)
    }
    pub fn mul_ground(&self, c: &Expr) -> PResult<Self> {
        Ok(self.with_rep(dmp_mul_ground(&self.rep, &EX(c.clone()), self.u())?))
    }
    pub fn quo_ground(&self, c: &Expr) -> PResult<Self> {
        Ok(self.with_rep(dmp_quo_ground(&self.rep, &EX(c.clone()), self.u())?))
    }
    pub fn exquo_ground(&self, c: &Expr) -> PResult<Self> {
        Ok(self.with_rep(dmp_exquo_ground(&self.rep, &EX(c.clone()), self.u())?))
    }
    /// `dmp_abs`: `Expression.__abs__` on every coefficient.
    pub fn abs(&self) -> Self {
        fn walk(f: &Dmp) -> Dmp {
            match f {
                Dmp::L(v) => Dmp::L(v.iter().map(EX::abs).collect()),
                Dmp::N(v) => Dmp::N(v.iter().map(walk).collect()),
            }
        }
        self.with_rep(walk(&self.rep))
    }
    /// `resultant(includePRS=True)`: the resultant and the subresultant
    /// polynomial remainder sequence.
    pub fn resultant_prs(&self, o: &Self) -> PResult<(Result<Expr, Self>, Vec<Self>)> {
        let (r, prs) = dmp_prs_resultant(&self.rep, &o.rep, self.u())?;
        Ok((self.drop_main(r)?, prs.into_iter().map(|p| self.with_rep(p)).collect()))
    }
    /// `dup_decompose`: a functional decomposition `f = f_1(f_2(...))`.
    pub fn decompose(&self) -> PResult<Vec<Self>> {
        if self.u() != 0 {
            return Err(reference("ValueError", "univariate polynomial expected"));
        }
        Ok(dup_decompose(self.rep.l())?.into_iter().map(|f| self.with_rep(Dmp::L(f))).collect())
    }
    pub fn clear_denoms(&self) -> PResult<(Expr, Self)> {
        let (c, p) = dmp_clear_denoms(&self.rep, self.u())?;
        Ok((c.0, self.with_rep(p)))
    }
    pub fn compose(&self, o: &Self) -> PResult<Self> {
        Ok(self.with_rep(dmp_compose(&self.rep, &o.rep, self.u())?))
    }
    pub fn shift(&self, a: &Expr) -> PResult<Self> {
        self.univariate()?;
        Ok(self.with_rep(dup_shift(&self.rep, &EX(a.clone()))?))
    }
    pub fn sqf_p(&self) -> PResult<bool> {
        dmp_sqf_p(&self.rep, self.u())
    }
    pub fn sqf_part(&self) -> PResult<Self> {
        Ok(self.with_rep(dmp_sqf_part(&self.rep, self.u())?))
    }
    pub fn sqf_list(&self, all: bool) -> PResult<(Expr, Vec<(Self, usize)>)> {
        let (c, f) = dmp_sqf_list(&self.rep, self.u(), all)?;
        Ok((c.0, f.into_iter().map(|(p, k)| (self.with_rep(p), k)).collect()))
    }
    pub fn sqf_list_include(&self, all: bool) -> PResult<Vec<(Self, usize)>> {
        Ok(dmp_sqf_list_include(&self.rep, self.u(), all)?.into_iter().map(|(p, k)| (self.with_rep(p), k)).collect())
    }
    pub fn gff_list(&self) -> PResult<Vec<(Self, usize)>> {
        self.univariate()?;
        Ok(dup_gff_list(&self.rep)?.into_iter().map(|(p, k)| (self.with_rep(p), k)).collect())
    }
    pub fn sturm(&self) -> PResult<Vec<Self>> {
        self.univariate()?;
        Ok(dup_sturm(&self.rep)?.into_iter().map(|p| self.with_rep(p)).collect())
    }
    pub fn terms_gcd(&self) -> (Vec<u32>, Self) {
        let (j, p) = dmp_terms_gcd(&self.rep, self.u());
        (j, self.with_rep(p))
    }
    pub fn subresultants(&self, o: &Self) -> PResult<Vec<Self>> {
        Ok(dmp_subresultants(&self.rep, &o.rep, self.u())?.into_iter().map(|p| self.with_rep(p)).collect())
    }
    /// The resultant: a ground expression for a univariate input, a
    /// polynomial in the remaining generators otherwise.
    pub fn resultant(&self, o: &Self) -> PResult<Result<Expr, Self>> {
        let r = dmp_resultant(&self.rep, &o.rep, self.u())?;
        self.drop_main(r)
    }
    pub fn discriminant(&self) -> PResult<Result<Expr, Self>> {
        let r = dmp_discriminant(&self.rep, self.u())?;
        self.drop_main(r)
    }
    /// An elimination result back over all generators, with the eliminated
    /// main generator at exponent zero.
    pub fn retain_main(&self, r: Result<Expr, Self>) -> PResult<Self> {
        let n = self.gens.len();
        let terms = match r {
            Ok(e) => vec![(vec![0; n], e)],
            Err(p) => p.terms().into_iter().map(|(m, c)| (std::iter::once(0).chain(m).collect(), c)).collect(),
        };
        Self::from_terms(terms, &self.gens)
    }
    /// Remove the first generator of a polynomial that does not involve it.
    pub fn eliminate_first(&self) -> PResult<Self> {
        if self.gens.len() < 2 || self.degree(0) > 0 {
            return Err(PolyError::Failed("elimination result still contains its first generator".into()));
        }
        Self::from_terms(self.terms().into_iter().map(|(m, c)| (m[1..].to_vec(), c)).collect(), &self.gens[1..])
    }
    fn drop_main(&self, r: Dmp) -> PResult<Result<Expr, Self>> {
        if self.u() == 0 {
            return Ok(Ok(r.l().first().cloned().unwrap_or_else(EX::zero).0));
        }
        Ok(Err(Self { rep: r, gens: self.gens[1..].to_vec(), ground: self.ground.clone() }))
    }
    /// `Poly.eval(gens[j], a)`: a ground expression or a polynomial in the
    /// remaining generators.
    pub fn eval(&self, j: usize, a: &Expr) -> PResult<Result<Expr, Self>> {
        let a = EX(a.clone());
        if self.u() == 0 {
            return Ok(Ok(dup_eval_ground(&self.rep, &a)?.0));
        }
        let rep = dmp_eval_in(&self.rep, &a, j, self.u())?;
        let mut gens = self.gens.clone();
        gens.remove(j);
        Ok(Err(Self { rep, gens, ground: self.ground.clone() }))
    }
    pub fn degree(&self, j: usize) -> i64 {
        dmp_degree_in(&self.rep, j, self.u())
    }
    pub fn degree_list(&self) -> Vec<i64> {
        dmp_degree_list(&self.rep, self.u())
    }
    /// The dense coefficient list of a univariate polynomial.
    pub fn all_coeffs(&self) -> Vec<Expr> {
        match &self.rep {
            Dmp::L(v) if v.is_empty() => vec![Expr::zero()],
            Dmp::L(v) => v.iter().map(|c| c.0.clone()).collect(),
            Dmp::N(_) => Vec::new(),
        }
    }
    /// The legacy metadata view: coefficients over the automatically
    /// detected coefficient generators, or opaque atoms when that fails.
    pub fn legacy(&self) -> super::polyobj::Poly {
        let terms = self.terms();
        if let Ok(p) = super::polyobj::Poly::from_coefficient_terms(terms.clone(), &self.gens) {
            return p;
        }
        let mut atoms: Vec<Expr> = Vec::new();
        let mut output = Vec::new();
        for (m, c) in &terms {
            let index = atoms.iter().position(|a| a == c).unwrap_or_else(|| {
                atoms.push(c.clone());
                atoms.len() - 1
            });
            let mut e = m.clone();
            e.extend((0..atoms.len()).map(|_| 0));
            output.push((e, index));
        }
        let n = self.gens.len() + atoms.len();
        let rep = super::mpoly::MPoly::from_terms(
            n,
            output
                .into_iter()
                .map(|(mut e, index)| {
                    e.resize(n, 0);
                    e[self.gens.len() + index] = 1;
                    (e, super::int::Rat::from_i64(1))
                })
                .collect(),
        );
        let mut gens = self.gens.clone();
        gens.extend(atoms);
        super::polyobj::Poly::new(rep, gens, self.gens.len(), super::polyobj::Domain::ZZ)
    }
}

/// Payload carried by an `EX` coefficient failure through the infallible
/// [`Coeff`](super::int::Coeff) interface of the sparse algorithms.
struct CoefficientFailure(PolyError);

fn unwind<T>(r: PResult<T>) -> T {
    match r {
        Ok(v) => v,
        // `resume_unwind` skips the panic hook: the failure is reported by
        // [`guarded`], not printed.
        Err(e) => std::panic::resume_unwind(Box::new(CoefficientFailure(e))),
    }
}

/// Run a sparse algorithm over `EX` coefficients, returning the first
/// coefficient failure as an error. Other panics propagate unchanged.
pub fn guarded<T>(f: impl FnOnce() -> T) -> PResult<T> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(v) => Ok(v),
        Err(payload) => match payload.downcast::<CoefficientFailure>() {
            Ok(failure) => Err(failure.0),
            Err(other) => std::panic::resume_unwind(other),
        },
    }
}

/// `EX` as the coefficient field of the sparse ring algorithms
/// (`PolyElement` over `EX`): zero tests are `bool(EX(a))` and every
/// operation is the reference `Expression` arithmetic. Use inside
/// [`guarded`].
impl super::int::Coeff for EX {
    fn zero() -> Self {
        EX::zero()
    }
    fn one() -> Self {
        EX::one()
    }
    fn from_i64(i: i64) -> Self {
        EX::int(i)
    }
    fn is_zero(&self) -> bool {
        !self.truthy()
    }
    fn is_one(&self) -> bool {
        EX::is_one(self)
    }
    fn add(&self, o: &Self) -> Self {
        unwind(EX::add(self, o))
    }
    fn sub(&self, o: &Self) -> Self {
        unwind(EX::sub(self, o))
    }
    fn mul(&self, o: &Self) -> Self {
        unwind(EX::mul(self, o))
    }
    fn neg(&self) -> Self {
        EX::neg(self)
    }
    fn is_negative(&self) -> bool {
        EX::is_negative(self)
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        Some(unwind(EX::div(self, o)))
    }
    fn to_rat(&self) -> super::int::Rat {
        match self.0.as_number().and_then(super::int::Rat::from_number) {
            Some(r) => r,
            None => unwind(Err(PolyError::Failed(format!("expected a rational coefficient, got {}", self.0)))),
        }
    }
}

impl super::int::Field for EX {
    fn inv(&self) -> Self {
        unwind(EX::one().div(self))
    }
    fn div(&self, o: &Self) -> Self {
        unwind(EX::div(self, o))
    }
}

impl ExpressionPoly {
    /// The sparse `PolyElement` of this polynomial.
    pub fn to_sparse(&self) -> super::mpoly::MPoly<EX> {
        let n = self.gens.len();
        super::mpoly::MPoly::from_terms(n, self.terms().into_iter().map(|(m, c)| (m, EX(c))).collect())
    }
    /// The dense polynomial of a sparse `PolyElement` in these generators.
    pub fn from_sparse(p: &super::mpoly::MPoly<EX>, gens: &[Expr]) -> PResult<Self> {
        let terms = (0..p.len()).map(|i| (p.exp(i).to_vec(), p.coeff(i).0.clone())).collect();
        Self::from_terms(terms, gens)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    fn p(s: &str, gens: &[&str]) -> ExpressionPoly {
        let gens: Vec<Expr> = gens.iter().map(|g| Expr::symbol(g)).collect();
        ExpressionPoly::from_expr(&parse_str(s).unwrap(), &gens).unwrap()
    }
    #[test]
    fn composite_parameters_follow_construct_domain() {
        let e = |s: &str| parse_str(s).unwrap();
        let names = |v: Vec<Expr>| v.iter().map(crate::printer::to_str).collect::<Vec<_>>();
        let (g, fractions) = composite_parameters(&[e("1/pi"), e("2*pi+1")]).unwrap().unwrap();
        assert_eq!((names(g), fractions), (vec!["pi".to_string()], true));
        let (g, fractions) = composite_parameters(&[e("exp(-a)"), e("exp(a)")]).unwrap().unwrap();
        assert_eq!((names(g), fractions), (vec!["exp(a)".to_string()], true));
        let (g, fractions) = composite_parameters(&[e("b"), e("a/b")]).unwrap().unwrap();
        assert_eq!((names(g), fractions), (vec!["a".to_string(), "b".to_string()], true));
        assert!(composite_parameters(&[e("sqrt(a)"), e("a")]).unwrap().is_none());
        assert!(composite_parameters(&[e("1/2"), e("I")]).unwrap().is_none());
    }

    #[test]
    fn reference_construct_domain_selection() {
        let e = |s: &str| parse_str(s).unwrap();
        assert!(construct_selects_expression(&[e("1"), e("2*sqrt(2)"), e("2")]).unwrap());
        assert!(construct_selects_expression(&[e("1.5*sqrt(2)")]).unwrap());
        assert!(construct_selects_expression(&[e("sin(a)"), e("cos(a)")]).unwrap());
        assert!(construct_selects_expression(&[e("sin(a)"), e("1/a")]).unwrap());
        assert!(!construct_selects_expression(&[e("sqrt(a)"), e("1")]).unwrap());
        assert!(!construct_selects_expression(&[e("pi"), e("1/2")]).unwrap());
        assert!(!construct_selects_expression(&[e("2+3*I"), e("1/2")]).unwrap());
        assert!(!construct_selects_expression(&[e("sin(a)"), e("b")]).unwrap());
    }
    #[test]
    fn radical_gcd_and_square_free_reduction() {
        let f = p("x**2+2*sqrt(2)*x+2", &["x"]);
        let g = p("x+sqrt(2)", &["x"]);
        assert_eq!(f.gcd(&g).unwrap().as_expr(), parse_str("x+sqrt(2)").unwrap());
        assert!(!f.sqf_p().unwrap());
        assert_eq!(f.sqf_part().unwrap().as_expr(), parse_str("x+sqrt(2)").unwrap());
        let (c, factors) = f.sqf_list(false).unwrap();
        assert_eq!(c, Expr::one());
        assert_eq!(factors.len(), 1);
        assert_eq!(factors[0].1, 2);
        assert_eq!(f.max_norm(), parse_str("2*sqrt(2)").unwrap());
    }
    #[test]
    fn division_by_a_monic_linear_factor() {
        let f = p("x**2", &["x"]);
        let g = p("x+1", &["x"]);
        let (q, r) = f.div(&g).unwrap();
        assert_eq!(q.as_expr(), parse_str("x-1").unwrap());
        assert_eq!(r.as_expr(), Expr::one());
    }
    #[test]
    fn multivariate_gcd_keeps_content() {
        let f = p("x*y+sqrt(2)*y", &["x", "y"]);
        let g = p("x**2*y+2*sqrt(2)*x*y+2*y", &["x", "y"]);
        assert_eq!(f.gcd(&g).unwrap().as_expr(), parse_str("x*y+sqrt(2)*y").unwrap());
    }
}
