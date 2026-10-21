//! Multivariate gcd, square-free decomposition and Trager factorization over
//! algebraic number fields.
//!
//! These follow the reference `dmp_ff_prs_gcd`, `dmp_sqf_list`,
//! `dmp_sqf_part`, `dmp_sqf_norm` and `dmp_ext_factor`: every result is
//! exact, gcds are made monic with respect to the lex-leading coefficient,
//! and factor lists use the reference `_sort_factors` order on the dense
//! recursive representation.
use super::algebraic::{NumberField, AN};
use super::convert::{PResult, PolyError};
use super::int::{Coeff, Field, Int, Rat};
use super::mpoly::MPoly;
use std::cmp::Ordering;
use std::sync::Arc;

/// Exact quotient in a polynomial ring over a field.
pub(crate) fn exact<C: Field>(a: &MPoly<C>, b: &MPoly<C>) -> PResult<MPoly<C>> {
    if let Some(q) = a.div_exact(b) {
        return Ok(q);
    }
    let (q, r) = super::fraction::recursive_division(a, b, 0, &|_| true);
    if r.is_zero() {
        Ok(q)
    } else {
        Err(PolyError::Failed(
            "inexact multivariate division over a field".into(),
        ))
    }
}

/// Pseudo-remainder of `f` by `g` in the variable `v`.
fn prem<C: Field>(f: &MPoly<C>, g: &MPoly<C>, v: usize) -> MPoly<C> {
    let n = f.nvars;
    let dg = g.degree(v);
    let lg = g.lc_in(v);
    let mut r = f.clone();
    while !r.is_zero() && r.degree(v) >= dg {
        super::factor_budget::event("algebraic pseudo-remainder");
        let dr = r.degree(v);
        let lr = r.lc_in(v);
        let mut e = vec![0; n];
        e[v] = (dr - dg) as u32;
        r = r.mul(&lg).sub(&g.mul(&lr).mul_term(&e, &C::one()));
    }
    r
}

/// Content and primitive part with respect to `v`, the content being the
/// monic gcd of the coefficients in the variables after `v`.
fn primitive<C: Field>(p: &MPoly<C>, v: usize) -> PResult<(MPoly<C>, MPoly<C>)> {
    let n = p.nvars;
    let mut content = MPoly::zero(n);
    for (_, c) in p.coeffs_in(v) {
        content = gcd_from(&content, &c, v + 1)?;
        if content.is_constant() {
            break;
        }
    }
    if content.is_zero() {
        return Ok((content, p.clone()));
    }
    let primitive = exact(p, &content)?;
    Ok((content, primitive))
}

/// Monic gcd of polynomials in which the variables before `v` are absent.
fn gcd_from<C: Field>(a: &MPoly<C>, b: &MPoly<C>, v: usize) -> PResult<MPoly<C>> {
    let n = a.nvars;
    if a.is_zero() {
        return Ok(b.monic());
    }
    if b.is_zero() {
        return Ok(a.monic());
    }
    let Some(v) = (v..n).find(|&i| a.degree(i) > 0 || b.degree(i) > 0) else {
        return Ok(MPoly::one(n));
    };
    let (ca, pa) = primitive(a, v)?;
    let (cb, pb) = primitive(b, v)?;
    let content = gcd_from(&ca, &cb, v + 1)?;
    let (pa, pb) = (pa.monic(), pb.monic());
    let (mut f, mut g) = if pa.degree(v) >= pb.degree(v) {
        (pa, pb)
    } else {
        (pb, pa)
    };
    while !g.is_zero() {
        let r = prem(&f, &g, v);
        f = g;
        // The primitive part is defined up to a unit of the field; making it
        // ground monic keeps the rational heights of the sequence bounded.
        g = if r.is_zero() { r } else { primitive(&r, v)?.1.monic() };
    }
    let part = if f.degree(v) <= 0 {
        MPoly::one(n)
    } else {
        f
    };
    Ok(part.mul(&content).monic())
}

/// Monic gcd over a field (zero when both inputs vanish).
pub fn gcd<C: Field>(a: &MPoly<C>, b: &MPoly<C>) -> PResult<MPoly<C>> {
    gcd_from(a, b, 0)
}

/// Square-free part made monic, as the reference `dmp_sqf_part`.
pub fn sqf_part<C: Field>(f: &MPoly<C>) -> PResult<MPoly<C>> {
    if f.is_zero() {
        return Ok(f.clone());
    }
    let mut g = f.clone();
    for i in 0..f.nvars {
        g = gcd(&g, &f.derivative(i))?;
    }
    Ok(exact(f, &g)?.monic())
}

fn sqf_from<C: Field>(f: &MPoly<C>, v: usize) -> PResult<Vec<(MPoly<C>, u32)>> {
    let n = f.nvars;
    if v >= n || f.is_zero() {
        return Ok(Vec::new());
    }
    let f = f.monic();
    let (content, f) = primitive(&f, v)?;
    let mut result: std::collections::BTreeMap<u32, MPoly<C>> = Default::default();
    if f.degree(v) > 0 {
        let h = f.derivative(v);
        let g = gcd(&f, &h)?;
        let mut p = exact(&f, &g)?;
        let mut q = exact(&h, &g)?;
        let mut i = 1;
        loop {
            super::factor_budget::event("algebraic square-free decomposition");
            let h = q.sub(&p.derivative(v));
            if h.is_zero() {
                result.insert(i, p);
                break;
            }
            let g = gcd(&p, &h)?;
            p = exact(&p, &g)?;
            q = exact(&h, &g)?;
            if g.degree(v) > 0 {
                result.insert(i, g);
            }
            i += 1;
        }
    }
    for (factor, i) in sqf_from(&content, v + 1)? {
        let entry = result.entry(i).or_insert_with(|| MPoly::one(n));
        *entry = entry.mul(&factor);
    }
    Ok(result.into_iter().map(|(i, p)| (p, i)).collect())
}

/// Square-free decomposition `(lc, [(factor, multiplicity)])` with monic
/// factors, as the reference `dmp_sqf_list` over a field.
pub fn sqf_list<C: Field>(f: &MPoly<C>) -> PResult<(C, Vec<(MPoly<C>, u32)>)> {
    if f.is_zero() {
        return Ok((C::zero(), Vec::new()));
    }
    Ok((f.lc(), sqf_from(f, 0)?))
}

/// Field norm `Res_T(m(T), f(x, T))` of a polynomial over `field`, up to a
/// nonzero rational constant.
fn norm(f: &MPoly<AN>, field: &NumberField) -> MPoly<Int> {
    let n = f.nvars;
    let modulus = MPoly::from_terms(
        n + 1,
        field
            .minimal
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(j, c)| {
                let mut e = vec![0; n + 1];
                e[0] = j as u32;
                (e, c.clone())
            })
            .collect(),
    );
    let value = MPoly::from_terms(
        n + 1,
        f.iter()
            .flat_map(|(e, a)| {
                a.coefficients
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_zero())
                    .map(move |(j, c)| {
                        let mut x = vec![j as u32];
                        x.extend_from_slice(e);
                        (x, c.clone())
                    })
            })
            .collect(),
    );
    let (_, m) = modulus.clear_denoms();
    let (_, value) = value.clear_denoms();
    let r = super::resultant::resultant_mz(&m, &value, 0);
    MPoly::from_terms(n, r.iter().map(|(e, c)| (e[1..].to_vec(), c.clone())).collect())
}

fn squarefree_mz(r: &MPoly<Int>) -> bool {
    if r.is_zero() {
        return true;
    }
    (0..r.nvars).all(|i| {
        let d = r.derivative(i);
        d.is_zero() || super::gcd::gcd_mz(r, &d).0.degree(i) <= 0
    })
}

/// Python-list comparison key of the reference dense recursive form.
#[derive(Clone, PartialEq, Eq)]
enum Dense {
    Leaf(Vec<Rat>),
    Node(Vec<Dense>),
}
impl Ord for Dense {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Dense::Leaf(a), Dense::Leaf(b)) => a.cmp(b),
            (Dense::Node(a), Dense::Node(b)) => a.cmp(b),
            (Dense::Leaf(_), Dense::Node(_)) => Ordering::Less,
            (Dense::Node(_), Dense::Leaf(_)) => Ordering::Greater,
        }
    }
}
impl PartialOrd for Dense {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn element_key(a: &AN) -> Vec<Rat> {
    let mut c = a.coefficients.clone();
    while c.last().is_some_and(|x| x.is_zero()) {
        c.pop();
    }
    c.reverse();
    c
}

fn dense_key(p: &MPoly<AN>, v: usize) -> Dense {
    dense_key_by(p, v, &element_key)
}

fn dense_key_by(p: &MPoly<AN>, v: usize, key: &dyn Fn(&AN) -> Vec<Rat>) -> Dense {
    if v >= p.nvars {
        return Dense::Leaf(key(&p.constant_value().unwrap_or_else(AN::zero)));
    }
    if p.is_zero() {
        return Dense::Node(Vec::new());
    }
    let groups = p.coeffs_in(v);
    let top = groups[0].0;
    let mut out = Vec::with_capacity(top as usize + 1);
    let mut groups = groups.into_iter().peekable();
    for k in (0..=top).rev() {
        if groups.peek().is_some_and(|(d, _)| *d == k) {
            out.push(dense_key_by(&groups.next().unwrap().1, v + 1, key));
        } else {
            out.push(dense_key_by(&MPoly::zero(p.nvars), v + 1, key));
        }
    }
    Dense::Node(out)
}

/// Compare two polynomials as the reference compares their dense
/// recursive lists.
pub fn compare_reps(a: &MPoly<AN>, b: &MPoly<AN>) -> Ordering {
    dense_key(a, 0).cmp(&dense_key(b, 0))
}

/// [`compare_reps`] with the ground elements compared by `key`.
pub fn compare_reps_by(a: &MPoly<AN>, b: &MPoly<AN>, key: impl Fn(&AN) -> Vec<Rat>) -> Ordering {
    dense_key_by(a, 0, &key).cmp(&dense_key_by(b, 0, &key))
}

/// Sort factors as the reference `_sort_factors(multiple=True)`: by the
/// length of the dense form in the main variable, the multiplicity and the
/// dense form itself.
pub fn sort_factors(factors: &mut [(MPoly<AN>, u32)]) {
    factors.sort_by_cached_key(|(p, k)| ((p.degree(0).max(-1) + 1), *k, dense_key(p, 0)));
}

/// `sort_factors` with the ground elements compared by `key`.
pub fn sort_factors_by(factors: &mut [(MPoly<AN>, u32)], key: impl Fn(&AN) -> Vec<Rat>) {
    factors.sort_by_cached_key(|(p, k)| ((p.degree(0).max(-1) + 1), *k, dense_key_by(p, 0, &key)));
}

fn trial_division(
    mut f: MPoly<AN>,
    factors: Vec<MPoly<AN>>,
) -> PResult<Vec<(MPoly<AN>, u32)>> {
    let mut result = Vec::with_capacity(factors.len());
    for factor in factors {
        let mut k = 0;
        while let Ok(q) = exact(&f, &factor) {
            f = q;
            k += 1;
        }
        if k == 0 {
            return Err(PolyError::Failed("trial division failed".into()));
        }
        result.push((factor, k));
    }
    sort_factors(&mut result);
    Ok(result)
}

/// Shift sequence of the reference `_dmp_sqf_norm_shifts`, truncated after
/// `bound` uniform shifts.
fn shifts(f: &MPoly<AN>, bound: usize) -> impl Iterator<Item = Vec<i64>> {
    let n = f.nvars;
    let mut order: Vec<(i64, usize)> = (0..n).map(|i| (f.degree(i), i)).collect();
    order.sort();
    let singles: Vec<Vec<i64>> = order
        .into_iter()
        .filter(|(d, _)| *d > 0)
        .map(|(_, i)| {
            let mut s = vec![0; n];
            s[i] = 1;
            s
        })
        .collect();
    std::iter::once(vec![0; n])
        .chain(singles)
        .chain((1..=bound as i64).map(move |j| vec![j; n]))
}

fn shifted(f: &MPoly<AN>, s: &[i64], alpha: &AN, sign: i64) -> MPoly<AN> {
    let mut g = f.clone();
    for (v, &k) in s.iter().enumerate() {
        if k != 0 {
            g = g.shift(v, &alpha.mul(&AN::from_i64(sign * k)));
        }
    }
    g
}

/// Trager factorization of a nonconstant polynomial that is not divisible by
/// any variable, as the reference `dmp_ext_factor`.
fn ext_factor(f: &MPoly<AN>, field: &Arc<NumberField>) -> PResult<(AN, Vec<(MPoly<AN>, u32)>)> {
    let lc = f.lc();
    let f = f.monic();
    if f.is_constant() {
        return Ok((lc, Vec::new()));
    }
    let original = f.clone();
    let f = sqf_part(&f)?;
    let alpha = AN::alpha(field.clone());
    let total = (0..f.nvars).map(|i| f.degree(i).max(0) as usize).sum::<usize>();
    let roots = field.degree().checked_mul(total).ok_or_else(|| {
        PolyError::Unsupported("algebraic factorization degree overflow".into())
    })?;
    let bound = roots
        .checked_mul(roots.saturating_sub(1))
        .ok_or_else(|| PolyError::Unsupported("algebraic factorization degree overflow".into()))?
        / 2
        + 1;
    for s in shifts(&f, bound) {
        super::factor_budget::event("algebraic norm shifts");
        let mut g = shifted(&f, &s, &alpha, -1);
        let r = norm(&g, field);
        if !squarefree_mz(&r) {
            continue;
        }
        let (_, factors) = super::factor_multi::factor_list_mz(&r);
        let candidates = if factors.len() <= 1 {
            vec![f.clone()]
        } else {
            let mut out = Vec::with_capacity(factors.len());
            for (h, _) in factors {
                let h = h.map_coeffs(|c| AN::rational(Rat::from(c.clone()), alpha.field.clone()));
                let h = gcd(&h, &g)?;
                g = exact(&g, &h)?;
                out.push(shifted(&h, &s, &alpha, 1));
            }
            out
        };
        return Ok((lc, trial_division(original, candidates)?));
    }
    Err(PolyError::Unsupported(
        "could not separate conjugate roots for exact algebraic factorization".into(),
    ))
}

/// Irreducible factorization `(lc, [(monic factor, multiplicity)])` over the
/// number field, as the reference `dmp_factor_list` for algebraic domains.
pub fn factor_list(f: &MPoly<AN>, field: &Arc<NumberField>) -> PResult<(AN, Vec<(MPoly<AN>, u32)>)> {
    let n = f.nvars;
    if f.is_zero() {
        return Ok((AN::zero(), Vec::new()));
    }
    // dmp_terms_gcd: strip the monomial content first.
    let mut j = vec![u32::MAX; n];
    for (e, _) in f.iter() {
        for v in 0..n {
            j[v] = j[v].min(e[v]);
        }
    }
    let monomial = MPoly::monomial(&j, AN::one());
    let reduced = exact(f, &monomial)?;
    let (lc, mut factors) = ext_factor(&reduced, field)?;
    for (v, &k) in j.iter().enumerate().rev() {
        if k > 0 {
            factors.insert(0, (MPoly::var(v, n), k));
        }
    }
    sort_factors(&mut factors);
    Ok((lc, factors))
}

#[cfg(test)]
mod tests {
    use super::super::algebraic::{AlgebraicPoly, AlgebraicScope};
    use super::*;
    use crate::expr::Expr;

    fn parse(s: &str) -> Expr {
        crate::parser::parse_str(s).unwrap()
    }

    #[test]
    fn rational_multivariate_gcd_and_sqf() {
        let x = MPoly::<Rat>::var(0, 2);
        let y = MPoly::<Rat>::var(1, 2);
        let a = x.add(&y).mul(&x.sub(&y)).mul(&x);
        let b = x.add(&y).mul(&x.add(&y)).scale(&Rat::from_i64(3));
        assert_eq!(gcd(&a, &b).unwrap(), x.add(&y));
        let f = x.add(&y).pow(2).mul(&y.pow(3));
        let (_, parts) = sqf_list(&f).unwrap();
        assert_eq!(parts, vec![(x.add(&y), 2), (y, 3)]);
    }

    #[test]
    fn trager_factorization_of_bivariate_quartic() {
        let x = parse("x");
        let y = parse("y");
        let field = NumberField::from_generators(&[parse("2^(1/4)")]).unwrap();
        let _scope = AlgebraicScope::new(field.clone());
        let p = AlgebraicPoly::from_expr(&parse("x^4-2*y^4"), &[x, y], field.clone()).unwrap();
        let (c, factors) = factor_list(&p.rep, &field).unwrap();
        assert_eq!(factors.len(), 3);
        let product = factors
            .iter()
            .fold(MPoly::constant(c, 2), |a, (b, k)| a.mul(&b.pow(*k)));
        assert_eq!(product, p.rep);
        assert!(factors.iter().all(|(f, _)| f.lc().is_one()));
    }
}
