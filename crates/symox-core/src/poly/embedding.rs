//! Certified embeddings of algebraic numbers in an arbitrary number field.
//!
//! A field `QQ(theta)` is given by its generator expression and the monic
//! minimal polynomial of `theta`. An element `e` belongs to the field
//! exactly when its minimal polynomial `p` has a linear factor
//! `x - r(theta)` over the field with `r(theta) == e`. The linear factors
//! come from the exact norm factorization over the field; the one equal to
//! `e` is identified on rigorous complex enclosures:
//!
//! * when `p` splits in the field, its roots are exactly the candidates,
//!   so a value rectangle meeting a single candidate rectangle decides;
//! * otherwise the roots of `p` are isolated first, and `e` equals a
//!   candidate exactly when both lie in the same isolating rectangle. No
//!   common rectangle proves that `e` is outside the field.
//!
//! This is the reference `to_number_field` with the factorization method
//! of `field_isomorphism`, using a certified comparison instead of a
//! floating-point one. A precision that never separates the rectangles is
//! reported as undecided, never as absence.
//!
//! `primitive_element` follows the reference algorithm: factor the current
//! minimal polynomial over the next extension, keep the factor vanishing
//! at the current generator and shift it until its norm is squarefree.

use std::cell::RefCell;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use super::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN};
use super::convert::{PResult, PolyError};
use super::dense;
use super::enclosure::{self, CBox, EncError};
use super::int::{Coeff, Field, Int, Rat};
use super::minpoly::{self, MinpolyError};
use crate::expr::{Expr, Kind};

const START_PREC: u32 = 64;
const MAX_PREC: u32 = 1 << 14;

const ABSENT: &str = "coefficient is not in the selected algebraic number field";
const OUTSIDE: &str = "coefficient is outside the selected algebraic subfield";

fn absent() -> PolyError {
    PolyError::Unsupported(ABSENT.into())
}

/// Whether `e` reports a proven absence from a number field, as opposed to
/// an undecided or unsupported membership.
pub fn is_absence(e: &PolyError) -> bool {
    matches!(e, PolyError::Unsupported(m) | PolyError::Failed(m) if m == ABSENT || m == OUTSIDE)
}

fn undecided(e: &Expr) -> PolyError {
    PolyError::Unsupported(format!(
        "could not certify the algebraic field membership of {e}"
    ))
}

/// The minimal polynomial error as a polynomial error, keeping the
/// reference message.
pub fn minpoly_error(e: MinpolyError) -> PolyError {
    PolyError::Unsupported(e.to_string())
}

/// Monic minimal polynomial over `QQ`, lowest degree first.
pub fn minimal_rational(e: &Expr) -> Result<Vec<Rat>, MinpolyError> {
    let p = minpoly::minimal_polynomial_qq(e)?;
    Ok(dense::monic(&dense::to_rat(&p)))
}

fn rational(e: &Expr) -> Option<Rat> {
    e.as_number().and_then(Rat::from_number)
}

/// `a*b mod m` over `QQ`.
fn mul_mod(a: &[Rat], b: &[Rat], m: &[Rat]) -> Vec<Rat> {
    dense::rem(&dense::mul(a, b), m)
}

fn inverse_mod(a: &[Rat], m: &[Rat]) -> PResult<Vec<Rat>> {
    if a.is_empty() {
        return Err(PolyError::Failed(
            "zero algebraic coefficient denominator".into(),
        ));
    }
    let (g, s, _) = dense::gcdex_field(a, m);
    if g != vec![Rat::one()] {
        return Err(PolyError::Failed("reducible number field modulus".into()));
    }
    Ok(s)
}

fn power_mod(a: &[Rat], k: i64, m: &[Rat]) -> PResult<Vec<Rat>> {
    let mut base = if k < 0 { inverse_mod(a, m)? } else { a.to_vec() };
    let mut k = k.unsigned_abs();
    let mut result = vec![Rat::one()];
    while k != 0 {
        if k & 1 != 0 {
            result = mul_mod(&result, &base, m);
        }
        k >>= 1;
        if k != 0 {
            base = mul_mod(&base, &base, m);
        }
    }
    Ok(result)
}

/// Roots in `QQ[T]/(minimal)` of the irreducible rational polynomial `p`.
fn linear_roots(p: &[Rat], field: &Arc<NumberField>) -> PResult<Vec<Vec<Rat>>> {
    let _scope = AlgebraicScope::new(field.clone());
    let coefficients: Vec<AN> = p
        .iter()
        .map(|c| AN::rational(c.clone(), field.clone()))
        .collect();
    let poly = AlgebraicPoly {
        rep: dense::to_mpoly(&coefficients, 0, 1),
        gens: vec![Expr::symbol("_x")],
        number_field: field.clone(),
        gaussian: None,
    };
    let mut roots = Vec::new();
    for f in poly.factor_squarefree()? {
        let f = dense::from_mpoly(&f.rep, 0);
        if f.len() == 2 {
            roots.push(f[0].mul(&f[1].inv()).neg().coefficients);
        }
    }
    Ok(roots)
}

type Key = (Vec<Rat>, Expr, Expr);
thread_local! {
    static CACHE: RefCell<FxHashMap<Key, Option<Vec<Rat>>>> = RefCell::new(FxHashMap::default());
}

/// Coefficients of `e` in the power basis of `generator`, `None` when `e`
/// is proven not to belong to the field.
fn embed(generator: &Expr, minimal: &[Rat], e: &Expr) -> PResult<Option<Vec<Rat>>> {
    let key = (minimal.to_vec(), generator.clone(), e.clone());
    if let Some(v) = CACHE.with(|c| c.borrow().get(&key).cloned()) {
        return Ok(v);
    }
    let result = embed_uncached(generator, minimal, e)?;
    CACHE.with(|c| {
        let mut c = c.borrow_mut();
        if c.len() > 4096 {
            c.clear();
        }
        c.insert(key, result.clone());
    });
    Ok(result)
}

fn embed_uncached(generator: &Expr, minimal: &[Rat], e: &Expr) -> PResult<Option<Vec<Rat>>> {
    let p = minimal_rational(e).map_err(minpoly_error)?;
    if p.len() == 2 {
        let mut r = vec![p[0].neg()];
        dense::trim(&mut r);
        return Ok(Some(r));
    }
    let degree = minimal.len() - 1;
    if degree % (p.len() - 1) != 0 {
        return Ok(None);
    }
    let field = NumberField::bare(generator.clone(), minimal.to_vec());
    let roots = linear_roots(&p, &field)?;
    if roots.is_empty() {
        return Ok(None);
    }
    let split = roots.len() == p.len() - 1;
    let (_, integer) = dense::clear_denoms(&p);
    let leading_first: Vec<Int> = integer.iter().rev().cloned().collect();
    let mut prec = START_PREC;
    loop {
        if let Some(answer) = attempt(generator, e, &roots, split, &leading_first, prec)? {
            return Ok(answer);
        }
        if prec >= MAX_PREC {
            return Err(undecided(e));
        }
        prec *= 2;
    }
}

fn enclose(e: &Expr, prec: u32) -> PResult<Option<CBox>> {
    match enclosure::enclose(e, prec) {
        Ok(b) => Ok(Some(b)),
        Err(EncError::Precision) => Ok(None),
        Err(EncError::NotAlgebraic(_)) => Err(undecided(e)),
    }
}

fn pairwise_disjoint(boxes: &[CBox]) -> bool {
    boxes
        .iter()
        .enumerate()
        .all(|(i, a)| boxes[i + 1..].iter().all(|b| a.disjoint(b)))
}

/// One certified comparison at `prec`; `None` when it does not decide.
fn attempt(
    generator: &Expr,
    e: &Expr,
    roots: &[Vec<Rat>],
    split: bool,
    f: &[Int],
    prec: u32,
) -> PResult<Option<Option<Vec<Rat>>>> {
    let Some(value) = enclose(e, prec)? else {
        return Ok(None);
    };
    let Some(theta) = enclose(generator, prec)? else {
        return Ok(None);
    };
    let candidates: Vec<CBox> = roots
        .iter()
        .map(|r| enclosure::horner_rat(r, &theta, prec))
        .collect();
    if split {
        if !pairwise_disjoint(&candidates) {
            return Ok(None);
        }
        let hits: Vec<usize> = (0..candidates.len())
            .filter(|&i| !candidates[i].disjoint(&value))
            .collect();
        return Ok((hits.len() == 1).then(|| Some(roots[hits[0]].clone())));
    }
    let degree = f.len() - 1;
    let mut isolating = Vec::with_capacity(degree);
    for k in 0..degree {
        match enclosure::root_enclosure(f, k, prec) {
            Some(b) => isolating.push(b),
            None => return Ok(None),
        }
    }
    if !pairwise_disjoint(&isolating) {
        return Ok(None);
    }
    let locate = |b: &CBox| {
        let hits: Vec<usize> = (0..degree).filter(|&k| !isolating[k].disjoint(b)).collect();
        (hits.len() == 1).then(|| hits[0])
    };
    let Some(target) = locate(&value) else {
        return Ok(None);
    };
    let mut found = None;
    for (r, c) in roots.iter().zip(&candidates) {
        let Some(k) = locate(c) else {
            return Ok(None);
        };
        if k == target {
            found = Some(r.clone());
        }
    }
    Ok(Some(found))
}

/// Membership of the numeric expression `e` in `QQ(generator)`:
/// `Some(coefficients)` (lowest power first) when it belongs to the field,
/// `None` when it is proven not to, an error when it is undecided.
///
/// Sums, products and integer powers are combined exactly from their
/// parts; any other subexpression, or a combination whose parts are not
/// all in the field, is embedded as a whole.
pub fn member(generator: &Expr, minimal: &[Rat], e: &Expr) -> PResult<Option<Vec<Rat>>> {
    if !crate::traverse::free_symbols(e).is_empty() {
        return Err(absent());
    }
    let key = (minimal.to_vec(), generator.clone(), e.clone());
    if let Some(v) = CACHE.with(|c| c.borrow().get(&key).cloned()) {
        return Ok(v);
    }
    match structural(generator, minimal, e) {
        Ok(Some(v)) => Ok(Some(v)),
        Ok(None) | Err(_) => embed(generator, minimal, e),
    }
}

fn structural(generator: &Expr, m: &[Rat], e: &Expr) -> PResult<Option<Vec<Rat>>> {
    if let Some(r) = rational(e) {
        let mut v = vec![r];
        dense::trim(&mut v);
        return Ok(Some(v));
    }
    if e == generator {
        return Ok(Some(dense::rem(&[Rat::zero(), Rat::one()], m)));
    }
    match e.kind() {
        Kind::Add(c, terms) => {
            let mut acc = vec![Rat::from_number(c).ok_or_else(absent)?];
            for (t, k) in terms {
                let k = Rat::from_number(k).ok_or_else(absent)?;
                let Some(v) = member(generator, m, t)? else {
                    return Ok(None);
                };
                acc = dense::add(&acc, &dense::scale(&v, &k));
            }
            dense::trim(&mut acc);
            Ok(Some(acc))
        }
        Kind::Mul(c, factors) => {
            let mut acc = vec![Rat::from_number(c).ok_or_else(absent)?];
            for (b, x) in factors {
                let factor = if x.is_one() {
                    b.clone()
                } else {
                    crate::ops::pow(b, x)
                };
                let Some(v) = member(generator, m, &factor)? else {
                    return Ok(None);
                };
                acc = mul_mod(&acc, &v, m);
            }
            dense::trim(&mut acc);
            Ok(Some(acc))
        }
        Kind::Pow(b, x) if x.as_i64().is_some() => {
            let k = x.as_i64().unwrap();
            let Some(v) = member(generator, m, b)? else {
                return Ok(None);
            };
            power_mod(&v, k, m).map(Some)
        }
        _ => embed(generator, m, e),
    }
}

/// Coefficients of `e` in `QQ(generator)`; absence and undecided cases
/// are both errors, with distinct messages.
pub fn coefficient(generator: &Expr, minimal: &[Rat], e: &Expr) -> PResult<Vec<Rat>> {
    member(generator, minimal, e)?.ok_or_else(absent)
}

/// The factor of `factors` (coefficients in `QQ(ext)`) vanishing at
/// `value`, certified by exclusion of all the others.
fn vanishing_factor(factors: &[Vec<AN>], value: &Expr, ext: &Expr) -> PResult<usize> {
    if factors.len() == 1 {
        return Ok(0);
    }
    let mut prec = START_PREC;
    loop {
        if let (Some(v), Some(t)) = (enclose(value, prec)?, enclose(ext, prec)?) {
            let zero: Vec<usize> = factors
                .iter()
                .enumerate()
                .filter(|(_, f)| {
                    let c: Vec<CBox> = f
                        .iter()
                        .map(|a| enclosure::horner_rat(&a.coefficients, &t, prec))
                        .collect();
                    !enclosure::horner(&c, &v, prec).excludes_zero()
                })
                .map(|(i, _)| i)
                .collect();
            if zero.len() == 1 {
                return Ok(zero[0]);
            }
        }
        if prec >= MAX_PREC {
            return Err(PolyError::Unsupported(format!(
                "multiple candidates for the minimal polynomial of {value}"
            )));
        }
        prec *= 2;
    }
}

/// `f(x + a)` over the current algebraic scope.
fn shift(p: &[AN], a: &AN) -> Vec<AN> {
    let mut out = Vec::new();
    for c in p.iter().rev() {
        out = dense::mul(&out, &[a.clone(), AN::one()]);
        out = dense::add(&out, std::slice::from_ref(c));
    }
    out
}

/// The reference `primitive_element(extension)`: a generator
/// `sum(c_i*extension[i])`, its monic minimal polynomial over `QQ` (lowest
/// degree first) and the integer weights `c_i`.
pub fn primitive_element(extension: &[Expr]) -> PResult<(Expr, Vec<Rat>, Vec<i64>)> {
    primitive_element_with_reps(extension).map(|(g, m, w, _)| (g, m, w))
}

/// `primitive_element` with the representation of every extension element
/// in the power basis of the primitive element (lowest power first), as
/// the reference `ex=True` computes them: at each step the new extension
/// element is the root of the linear `gcd(f(theta, a), p(a))` over the new
/// field, and the previous representations are evaluated at the old
/// generator. The representations are recorded as certified embeddings.
pub fn primitive_element_with_reps(
    extension: &[Expr],
) -> PResult<(Expr, Vec<Rat>, Vec<i64>, Vec<Vec<Rat>>)> {
    let Some(first) = extension.first() else {
        return Err(PolyError::Unsupported(
            "an algebraic extension needs a generator".into(),
        ));
    };
    let mut generator = first.clone();
    let mut weights = vec![1];
    let mut g = minimal_rational(first).map_err(minpoly_error)?;
    let mut reps = vec![reduced(vec![Rat::zero(), Rat::one()], &g)];
    for ext in &extension[1..] {
        if let Some(r) = rational(ext) {
            weights.push(0);
            reps.push(reduced(vec![r], &g));
            continue;
        }
        let p = minimal_rational(ext).map_err(minpoly_error)?;
        let field = NumberField::bare(ext.clone(), p.clone());
        let (f, s, norm) = {
            let _scope = AlgebraicScope::new(field.clone());
            let lifted: Vec<AN> = g
                .iter()
                .map(|c| AN::rational(c.clone(), field.clone()))
                .collect();
            let poly = AlgebraicPoly {
                rep: dense::to_mpoly(&lifted, 0, 1),
                gens: vec![Expr::symbol("_x")],
                number_field: field.clone(),
                gaussian: None,
            };
            let factors: Vec<Vec<AN>> = poly
                .factor_squarefree()?
                .iter()
                .map(|f| dense::monic(&dense::from_mpoly(&f.rep, 0)))
                .collect();
            let mut f = factors[vanishing_factor(&factors, &generator, ext)?].clone();
            let alpha = AN::alpha(field.clone());
            let mut s = 0i64;
            let norm = loop {
                let norm = super::algebraic::polynomial_norm(&f, &field)?;
                if dense::gcd_field(&norm, &dense::derivative(&norm)).len() <= 1 {
                    break norm;
                }
                f = shift(&f, &alpha.neg());
                s += 1;
            };
            (f, s, norm)
        };
        if s != 0 {
            generator = crate::ops::add2(
                &generator,
                &crate::ops::mul2(&Expr::int(s), ext),
            );
        }
        weights.push(s);
        g = dense::monic(&norm);
        let erep = switched_root(&f, &p, &generator, &g)?;
        let mut old = reduced(vec![Rat::zero(), Rat::one()], &g);
        old = dense::sub(&old, &dense::scale(&erep, &Rat::from_i64(s)));
        dense::trim(&mut old);
        reps = reps
            .iter()
            .map(|r| {
                let mut acc: Vec<Rat> = Vec::new();
                for c in r.iter().rev() {
                    acc = dense::add(&mul_mod(&acc, &old, &g), std::slice::from_ref(c));
                }
                reduced(acc, &g)
            })
            .collect();
        reps.push(erep);
    }
    for (ext, rep) in extension.iter().zip(&reps) {
        record(&generator, &g, ext, rep.clone());
    }
    Ok((generator, g, weights, reps))
}

fn reduced(mut p: Vec<Rat>, m: &[Rat]) -> Vec<Rat> {
    if p.len() >= m.len() {
        p = dense::rem(&p, m);
    }
    dense::trim(&mut p);
    p
}

/// The extension element `a` in `K = QQ(theta)`, from `f(x, a)` over
/// `QQ(a)` with `f(theta, a) == 0`: the root of the linear
/// `gcd(f(theta, y), p(y))` over `K` (the reference `_switch_domain` and
/// `_linsolve`).
fn switched_root(f: &[AN], p: &[Rat], theta: &Expr, minimal: &[Rat]) -> PResult<Vec<Rat>> {
    let k = NumberField::bare(theta.clone(), minimal.to_vec());
    let _scope = AlgebraicScope::new(k.clone());
    let width = f.iter().map(|c| c.coefficients.len()).max().unwrap_or(0);
    let h: Vec<AN> = (0..width)
        .map(|j| {
            let column: Vec<Rat> = f
                .iter()
                .map(|c| c.coefficients.get(j).cloned().unwrap_or_else(Rat::zero))
                .collect();
            AN::new(column, k.clone())
        })
        .collect();
    let p: Vec<AN> = p.iter().map(|c| AN::rational(c.clone(), k.clone())).collect();
    let mut h = h;
    dense::trim(&mut h);
    let gcd = dense::monic(&dense::gcd_field(&h, &p));
    if gcd.len() != 2 {
        return Err(PolyError::Failed(
            "the primitive element does not determine the extension element".into(),
        ));
    }
    Ok(gcd[0].neg().coefficients)
}

/// Records the certified embedding of `e` in `QQ(generator)`.
fn record(generator: &Expr, minimal: &[Rat], e: &Expr, rep: Vec<Rat>) {
    let key = (minimal.to_vec(), generator.clone(), e.clone());
    CACHE.with(|c| {
        let mut c = c.borrow_mut();
        if c.len() > 4096 {
            c.clear();
        }
        c.insert(key, Some(rep));
    });
}

/// Why `minimal_polynomial_over` produced no polynomial.
#[derive(Clone, Debug)]
pub enum OverError {
    Minpoly(MinpolyError),
    Poly(PolyError),
}

/// Rational content of an expression, the reference `Expr.primitive()`
/// content: the gcd of the numerators of the term coefficients over the
/// lcm of their denominators.
fn rational_content(terms: &[Expr]) -> Rat {
    let mut n = Int::from(0);
    let mut d = Int::from(1);
    for t in terms {
        let (c, _) = crate::ops::as_coeff_mul(t);
        if let Some(c) = Rat::from_number(&c) {
            n = n.gcd(&c.n);
            d = d.lcm(&c.d);
        }
    }
    if n.is_zero() {
        Rat::one()
    } else {
        Rat::new(n.abs(), d)
    }
}

/// Minimal polynomial of `ex` over the number field `field` as an
/// expression in `x`: the irreducible factor over the field of its minimal
/// polynomial over `QQ` that vanishes at `ex`, with the rational content
/// removed and a positive leading coefficient like the reference.
pub fn minimal_polynomial_over(ex: &Expr, field: &Arc<NumberField>, x: &Expr) -> Result<Expr, OverError> {
    let p = minimal_rational(ex).map_err(OverError::Minpoly)?;
    let _scope = AlgebraicScope::new(field.clone());
    let lifted: Vec<AN> = p
        .iter()
        .map(|c| AN::rational(c.clone(), field.clone()))
        .collect();
    let poly = AlgebraicPoly {
        rep: dense::to_mpoly(&lifted, 0, 1),
        gens: vec![x.clone()],
        number_field: field.clone(),
        gaussian: None,
    };
    let factors: Vec<Vec<AN>> = if p.len() == 2 {
        vec![lifted]
    } else {
        poly.factor_squarefree()
            .map_err(OverError::Poly)?
            .iter()
            .map(|f| dense::monic(&dense::from_mpoly(&f.rep, 0)))
            .collect()
    };
    let f = &factors[vanishing_factor(&factors, ex, &field.generator).map_err(OverError::Poly)?];
    let mut terms = Vec::new();
    for (k, c) in f.iter().enumerate().rev() {
        if c.is_zero() {
            continue;
        }
        let power = crate::ops::powi(x, k as i64);
        let value = c.as_expr();
        let parts = match value.kind() {
            Kind::Add(..) => crate::traverse::ordered_args(&value),
            _ => vec![value],
        };
        terms.extend(parts.iter().map(|t| crate::ops::mul2(t, &power)));
    }
    let content = rational_content(&terms);
    let scale = Expr::num(content.inv().to_number());
    let mut result = crate::ops::add(&terms.iter().map(|t| crate::ops::mul2(&scale, t)).collect::<Vec<_>>());
    let leading = crate::ops::mul2(&scale, &f.last().map(AN::as_expr).unwrap_or_else(Expr::one));
    if crate::assumptions::is_negative(&leading) == Some(true) {
        result = crate::ops::neg(&result);
    }
    Ok(result)
}


#[cfg(test)]
mod tests {
    use super::*;
    fn expr(s: &str) -> Expr {
        crate::parser::parse_str(s).unwrap()
    }
    fn rats(v: &[i64]) -> Vec<Rat> {
        v.iter().map(|&c| Rat::from_i64(c)).collect()
    }

    #[test]
    fn quartic_root_field_membership() {
        let theta = expr("2^(1/4)");
        let m = rats(&[-2, 0, 0, 0, 1]);
        assert_eq!(coefficient(&theta, &m, &expr("sqrt(2)")).unwrap(), rats(&[0, 0, 1]));
        assert_eq!(
            coefficient(&theta, &m, &expr("1/(1+2^(1/4))")).unwrap(),
            rats(&[-1, 1, -1, 1])
        );
        assert_eq!(member(&theta, &m, &expr("sqrt(3)")).unwrap(), None);
        assert_eq!(member(&theta, &m, &expr("I")).unwrap(), None);
        assert_eq!(member(&theta, &m, &expr("I*2^(1/4)")).unwrap(), None);
        assert_eq!(
            coefficient(&theta, &m, &expr("sqrt(2+2*sqrt(2))^2")).unwrap(),
            rats(&[2, 0, 2])
        );
    }

    #[test]
    fn nested_radicals_in_quadratic_fields() {
        let theta = expr("sqrt(2)");
        let m = rats(&[-2, 0, 1]);
        assert_eq!(coefficient(&theta, &m, &expr("sqrt(3+2*sqrt(2))")).unwrap(), rats(&[1, 1]));
        assert_eq!(member(&theta, &m, &expr("sqrt(3)")).unwrap(), None);
    }

    #[test]
    fn reference_primitive_elements() {
        let (g, m, w) = primitive_element(&[expr("2^(1/4)"), expr("sqrt(3)")]).unwrap();
        assert_eq!(g, expr("2^(1/4)+sqrt(3)"));
        assert_eq!(m, rats(&[49, 0, -180, 0, 50, 0, -12, 0, 1]));
        assert_eq!(w, vec![1, 1]);
        let (_, m, w) = primitive_element(&[expr("sqrt(2)"), expr("2^(1/4)")]).unwrap();
        assert_eq!(m, rats(&[2, -8, -4, 0, 1]));
        assert_eq!(w, vec![1, 1]);
        let (_, m, _) = primitive_element(&[expr("2^(1/4)"), expr("I")]).unwrap();
        assert_eq!(m, rats(&[1, 0, 28, 0, 2, 0, 4, 0, 1]));
    }

    #[test]
    fn primitive_element_representations_match_the_reference() {
        let q = |n: i64, d: i64| Rat::new(Int::from(n), Int::from(d));
        let (_, _, _, reps) =
            primitive_element_with_reps(&[expr("2^(1/4)"), expr("sqrt(3)")]).unwrap();
        let first = [0, -2567, 0, 855, 0, -187, 0, 15].map(|n| q(n, 1064));
        let second = [0, 3631, 0, -855, 0, 187, 0, -15].map(|n| q(n, 1064));
        assert_eq!(reps, vec![first.to_vec(), second.to_vec()]);
        let r = super::super::rootof::complex_root(&expr("x^5-x-1"), 0).unwrap();
        let (g, m, _, reps) = primitive_element_with_reps(&[r.clone(), expr("sqrt(2)")]).unwrap();
        let sqrt2 = [
            -432925246, 3328959399, -898953585, -2179515340, 8200610, 788960456, -12198575,
            -199591540, 1203600, 19545040,
        ]
        .map(|n| q(n, 903370613));
        assert_eq!(reps[1], sqrt2.to_vec());
        // The recorded representations answer membership directly.
        assert_eq!(coefficient(&g, &m, &r).unwrap(), reps[0]);
    }

    #[test]
    fn minimal_polynomials_over_number_fields() {
        let x = Expr::symbol("x");
        let k = NumberField::from_generator(&expr("sqrt(2)")).unwrap();
        for (e, expected) in [
            ("2^(1/4)", "x^2 - sqrt(2)"),
            ("sqrt(2)+sqrt(3)", "x^2 - 2*sqrt(2)*x - 1"),
            ("sqrt(3)", "x^2 - 3"),
            ("sqrt(2)", "x - sqrt(2)"),
            ("1/2", "2*x - 1"),
        ] {
            assert_eq!(minimal_polynomial_over(&expr(e), &k, &x).unwrap(), expr(expected), "{e}");
        }
    }
}

