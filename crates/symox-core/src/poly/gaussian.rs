//! Polynomials over the Gaussian integers `ZZ_I` and rationals `QQ_I`.
//!
//! Elements are stored in the number field `QQ<I>` as `x + y*I`. The
//! Gaussian rationals share every field algorithm with `QQ<I>`; the ring
//! `ZZ_I` follows the reference `densetools`/`euclidtools`/`sqfreetools`
//! ring algorithms literally: contents are normalized gcds, gcds of
//! univariate polynomials are `c*u*pp` with `u` the canonical unit of the
//! primitive part's leading coefficient, multivariate gcds are normalized
//! by their ground leading coefficient, `lcm`, `sqf_list` and the
//! `dup_zz_i_factor` content extraction keep the reference units.
use super::algebraic::{AlgebraicScope, NumberField, AN};
use super::convert::{PResult, PolyError};
use super::int::{Coeff, Field, Int, Rat};
use super::mpoly::MPoly;
use crate::expr::Expr;
use std::sync::Arc;

thread_local! {
    static FIELD: Arc<NumberField> =
        NumberField::from_generator(&Expr::i()).expect("QQ<I> is a quadratic field");
}

/// The field `QQ<I>` carrying Gaussian coefficients.
pub fn field() -> Arc<NumberField> {
    FIELD.with(Arc::clone)
}

/// Whether `field` is the representation of `QQ<I>` used for Gaussian
/// coefficients.
pub fn is_gaussian_field(field: &NumberField) -> bool {
    field.minimal == [Rat::one(), Rat::zero(), Rat::one()] && field.generator == Expr::i()
}

pub fn parts(a: &AN) -> (Rat, Rat) {
    (
        a.coefficients.first().cloned().unwrap_or_else(Rat::zero),
        a.coefficients.get(1).cloned().unwrap_or_else(Rat::zero),
    )
}

pub fn make(x: Rat, y: Rat) -> AN {
    AN::new(vec![x, y], field())
}

pub fn is_integral(a: &AN) -> bool {
    a.coefficients.iter().all(Rat::is_integer)
}

pub fn poly_is_integral(p: &MPoly<AN>) -> bool {
    p.coeffs().iter().all(is_integral)
}

/// The reference coercion failure for a conversion to `ZZ_I`: the first
/// non-integral component, scanning terms from the leading one and the real
/// part before the imaginary part.
pub fn integrality_error(p: &MPoly<AN>) -> Option<String> {
    p.coeffs().iter().find_map(|a| {
        let (x, y) = parts(a);
        [x, y].into_iter().find(|q| !q.is_integer()).map(|q| format!(
            "expected an integer, got {}",
            crate::printer::to_str(&crate::expr::Expr::num(q.to_number()))
        ))
    })
}

/// The reference `GaussianElement.quadrant`.
pub fn quadrant(a: &AN) -> usize {
    let (x, y) = parts(a);
    if y.is_zero() {
        if x.is_negative() {
            2
        } else {
            0
        }
    } else if !y.is_negative() {
        if !x.is_zero() && !x.is_negative() {
            0
        } else {
            1
        }
    } else if x.is_negative() {
        2
    } else {
        3
    }
}

/// `units[-quadrant]` with `units = (1, I, -1, -I)`.
pub fn canonical_unit(a: &AN) -> AN {
    let one = Rat::one();
    let zero = Rat::zero();
    match quadrant(a) {
        0 => make(one, zero),
        1 => make(zero, one.neg()),
        2 => make(one.neg(), zero),
        _ => make(zero, one),
    }
}

pub fn normalize(a: &AN) -> AN {
    a.mul(&canonical_unit(a))
}

/// Rounded Euclidean division of Gaussian integers, the reference
/// `GaussianInteger.__divmod__`.
pub fn divmod(a: &AN, b: &AN) -> (AN, AN) {
    let (x, y) = parts(b);
    let (sx, sy) = parts(a);
    let a1 = sx.mul(&x).add(&sy.mul(&y));
    let b1 = sx.neg().mul(&y).add(&sy.mul(&x));
    let c = x.mul(&x).add(&y.mul(&y));
    let two = Rat::from_i64(2);
    let round = |v: &Rat| {
        let shifted = two.mul(v).add(&c).div(&two.mul(&c));
        Rat::from(shifted.floor())
    };
    let q = make(round(&a1), round(&b1));
    let r = a.sub(&q.mul(b));
    (q, r)
}

/// Normalized Euclidean gcd in `ZZ_I`.
pub fn gcd(a: &AN, b: &AN) -> AN {
    let (mut a, mut b) = (a.clone(), b.clone());
    while !b.is_zero() {
        let r = divmod(&a, &b).1;
        a = b;
        b = r;
    }
    normalize(&a)
}

/// `(a*b)//gcd(a, b)`, not normalized.
pub fn lcm(a: &AN, b: &AN) -> AN {
    let g = gcd(a, b);
    if g.is_zero() {
        return g;
    }
    divmod(&a.mul(b), &g).0
}

/// Positive integer denominator of a Gaussian rational.
pub fn denom(a: &AN) -> Int {
    let (x, y) = parts(a);
    x.d.lcm(&y.d)
}

/// Gaussian integer `a*denom(a)`.
pub fn numer(a: &AN) -> AN {
    a.mul(&AN::rational(Rat::from(denom(a)), a.field.clone()))
}

/// The reference `QQ_I.gcd`: `gcd(numer)/lcm(denom)`.
pub fn field_gcd(a: &AN, b: &AN) -> AN {
    let p = gcd(&numer(a), &numer(b));
    let q = denom(a).lcm(&denom(b));
    p.mul(&AN::rational(Rat::new(Int::one(), q), p.field.clone()))
}

/// Ground content: the normalized gcd of all coefficients, folded from
/// zero (`dmp_ground_content`); `field` selects the `QQ_I` gcd.
pub fn ground_content(p: &MPoly<AN>, ring: bool) -> AN {
    let mut content = AN::zero();
    for c in p.coeffs() {
        content = if ring {
            gcd(&content, c)
        } else {
            field_gcd(&content, c)
        };
        if content.is_one() {
            break;
        }
    }
    content
}

/// `dmp_ground_primitive`.
pub fn ground_primitive(p: &MPoly<AN>, ring: bool) -> (AN, MPoly<AN>) {
    if p.is_zero() {
        return (AN::zero(), p.clone());
    }
    let content = ground_content(p, ring);
    if content.is_one() {
        return (content, p.clone());
    }
    let inverse = content.inv();
    (content, p.scale(&inverse))
}

/// `dmp_clear_denoms` into `ZZ_I`: the positive integer lcm of the
/// coefficient denominators and the scaled polynomial.
pub fn clear_denoms(p: &MPoly<AN>) -> (Int, MPoly<AN>) {
    let _scope = AlgebraicScope::new(field());
    let mut common = Int::one();
    for c in p.coeffs() {
        common = common.lcm(&denom(c));
    }
    if common.is_one() {
        return (common, p.clone());
    }
    let scale = AN::rational(Rat::from(common.clone()), field());
    (common, p.scale(&scale))
}

fn exact(a: &MPoly<AN>, b: &MPoly<AN>) -> MPoly<AN> {
    super::algebraic_multi::exact(a, b).expect("exact division by a Gaussian divisor")
}

fn constant(c: AN, n: usize) -> MPoly<AN> {
    MPoly::constant(c, n)
}

/// A gcd of `f` and `g` in `ZZ_I[x...]`, defined up to a unit:
/// `gcd(contents)*pp(gcd over QQ_I)`.
fn true_gcd(f: &MPoly<AN>, g: &MPoly<AN>) -> PResult<MPoly<AN>> {
    let c = gcd(&ground_content(f, true), &ground_content(g, true));
    let h = super::algebraic_multi::gcd(f, g)?;
    let (_, h) = clear_denoms(&h);
    let (_, h) = ground_primitive(&h, true);
    Ok(h.scale(&c))
}

/// Main-variable content at level `v` (`dmp_content`): the gcd of the
/// coefficients in `x_v`, folded from the leading one.
fn content_at(p: &MPoly<AN>, v: usize) -> PResult<MPoly<AN>> {
    let groups = p.coeffs_in(v);
    let Some((_, first)) = groups.first() else {
        return Ok(MPoly::zero(p.nvars));
    };
    let mut content = first.clone();
    for (_, c) in &groups[1..] {
        content = inner_gcd_at(&content, c, v + 1)?.0;
        if content.is_one() {
            break;
        }
    }
    Ok(content)
}

/// `dmp_primitive` at level `v`.
fn primitive_at(p: &MPoly<AN>, v: usize) -> PResult<(MPoly<AN>, MPoly<AN>)> {
    let content = content_at(p, v)?;
    if p.is_zero() || content.is_one() {
        return Ok((content, p.clone()));
    }
    Ok((content.clone(), exact(p, &content)))
}

/// The reference `dmp_inner_gcd` over `ZZ_I` for polynomials in which the
/// variables before `v` are absent: `(h, f/h, g/h)`.
fn inner_gcd_at(
    f: &MPoly<AN>,
    g: &MPoly<AN>,
    v: usize,
) -> PResult<(MPoly<AN>, MPoly<AN>, MPoly<AN>)> {
    let n = f.nvars;
    let minus_one = constant(AN::one().neg(), n);
    match (f.is_zero(), g.is_zero()) {
        (true, true) => return Ok((MPoly::zero(n), MPoly::zero(n), MPoly::zero(n))),
        // `K.is_nonnegative` is false for every Gaussian element.
        (true, false) => return Ok((g.neg(), MPoly::zero(n), minus_one)),
        (false, true) => return Ok((f.neg(), minus_one, MPoly::zero(n))),
        _ => {}
    }
    let h = if v + 1 >= n {
        // dup_rr_prs_gcd
        let c = gcd(&ground_content(f, true), &ground_content(g, true));
        let h = super::algebraic_multi::gcd(f, g)?;
        let (_, h) = clear_denoms(&h);
        let (_, h) = ground_primitive(&h, true);
        let unit = c.mul(&canonical_unit(&h.lc()));
        h.scale(&unit)
    } else if f.is_one() || g.is_one() {
        return Ok((MPoly::one(n), f.clone(), g.clone()));
    } else {
        let df = f.degree(v);
        let dg = g.degree(v);
        if df > 0 && dg > 0 {
            // dmp_rr_prs_gcd
            let h = true_gcd(f, g)?;
            let unit = canonical_unit(&h.lc());
            if unit.is_one() {
                h
            } else {
                h.scale(&unit)
            }
        } else {
            // _dmp_simplify_gcd
            let (a, b) = if df == 0 && dg == 0 {
                (f.clone(), g.clone())
            } else if df == 0 {
                (f.clone(), content_at(g, v)?)
            } else {
                (content_at(f, v)?, g.clone())
            };
            inner_gcd_at(&a, &b, v + 1)?.0
        }
    };
    let cff = exact(f, &h);
    let cfg = exact(g, &h);
    Ok((h, cff, cfg))
}

/// `(gcd, f/gcd, g/gcd)` in `ZZ_I[x...]`, with the reference normalization.
pub fn cofactors(
    f: &MPoly<AN>,
    g: &MPoly<AN>,
) -> PResult<(MPoly<AN>, MPoly<AN>, MPoly<AN>)> {
    let _scope = AlgebraicScope::new(field());
    inner_gcd_at(f, g, 0)
}

pub fn ring_gcd(f: &MPoly<AN>, g: &MPoly<AN>) -> PResult<MPoly<AN>> {
    Ok(cofactors(f, g)?.0)
}

/// `dup_rr_lcm`/`dmp_rr_lcm`.
pub fn ring_lcm(f: &MPoly<AN>, g: &MPoly<AN>) -> PResult<MPoly<AN>> {
    let _scope = AlgebraicScope::new(field());
    let n = f.nvars;
    if f.is_zero() || g.is_zero() {
        return Ok(MPoly::zero(n));
    }
    let (fc, fp) = ground_primitive(f, true);
    let (gc, gp) = ground_primitive(g, true);
    let c = lcm(&fc, &gc);
    let h = exact(&fp.mul(&gp), &inner_gcd_at(&fp, &gp, 0)?.0);
    if n == 1 {
        let unit = canonical_unit(&h.lc());
        Ok(h.scale(&c.mul(&unit)))
    } else {
        Ok(h.scale(&c))
    }
}

/// Ring division `dmp_rr_div`: the quotient stops at the first leading
/// coefficient that the divisor's does not divide in `ZZ_I`.
pub fn ring_div(f: &MPoly<AN>, g: &MPoly<AN>) -> PResult<(MPoly<AN>, MPoly<AN>)> {
    if g.is_zero() {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let _scope = AlgebraicScope::new(field());
    Ok(super::fraction::recursive_division(f, g, 0, &is_integral))
}

/// `dmp_sqf_list` over `ZZ_I` in the variables from `v`.
fn sqf_list_at(f: &MPoly<AN>, v: usize) -> PResult<(AN, Vec<(MPoly<AN>, u32)>)> {
    let n = f.nvars;
    let (coeff, f) = ground_primitive(f, true);
    if v + 1 >= n {
        // dup_sqf_list
        if f.degree(v) <= 0 {
            return Ok((coeff, Vec::new()));
        }
        let mut result = Vec::new();
        let h = f.derivative(v);
        let (_, mut p, mut q) = inner_gcd_at(&f, &h, v)?;
        let mut i = 1;
        loop {
            super::factor_budget::event("Gaussian square-free decomposition");
            let d = p.derivative(v);
            let h = q.sub(&d);
            if h.is_zero() {
                result.push((p, i));
                break;
            }
            let (g, p1, q1) = inner_gcd_at(&p, &h, v)?;
            p = p1;
            q = q1;
            if g.degree(v) > 0 {
                result.push((g, i));
            }
            i += 1;
        }
        return Ok((coeff, result));
    }
    let degree = f.degree(v);
    if degree < 0 {
        return Ok((coeff, Vec::new()));
    }
    let (content, f) = primitive_at(&f, v)?;
    let mut result: std::collections::BTreeMap<u32, MPoly<AN>> = Default::default();
    if degree != 0 {
        let h = f.derivative(v);
        let (_, mut p, mut q) = inner_gcd_at(&f, &h, v)?;
        let mut i = 1;
        loop {
            super::factor_budget::event("Gaussian square-free decomposition");
            let d = p.derivative(v);
            let h = q.sub(&d);
            if h.is_zero() {
                result.insert(i, p);
                break;
            }
            let (g, p1, q1) = inner_gcd_at(&p, &h, v)?;
            p = p1;
            q = q1;
            if g.degree(v) > 0 {
                result.insert(i, g);
            }
            i += 1;
        }
    }
    let (content_coeff, content_factors) = sqf_list_at(&content, v + 1)?;
    let coeff = coeff.mul(&content_coeff);
    for (factor, i) in content_factors {
        let entry = result.entry(i).or_insert_with(|| MPoly::one(n));
        *entry = entry.mul(&factor);
    }
    Ok((coeff, result.into_iter().map(|(i, p)| (p, i)).collect()))
}

/// Square-free decomposition over `ZZ_I` (`dmp_sqf_list`).
pub fn ring_sqf_list(f: &MPoly<AN>) -> PResult<(AN, Vec<(MPoly<AN>, u32)>)> {
    let _scope = AlgebraicScope::new(field());
    if f.is_zero() {
        return Ok((AN::zero(), Vec::new()));
    }
    sqf_list_at(f, 0)
}

/// The reference `_sort_factors` with Gaussian elements ordered by
/// `[y, x]`.
pub fn sort_factors(factors: &mut [(MPoly<AN>, u32)]) {
    super::algebraic_multi::sort_factors_by(factors, element_key);
}

/// The reference `GaussianElement` order, `[y, x]`.
pub fn element_key(a: &AN) -> Vec<Rat> {
    let (x, y) = parts(a);
    vec![y, x]
}

/// Compare two Gaussian polynomials as the reference compares their dense
/// recursive lists.
pub fn compare_reps(a: &MPoly<AN>, b: &MPoly<AN>) -> std::cmp::Ordering {
    let _scope = AlgebraicScope::new(field());
    super::algebraic_multi::compare_reps_by(a, b, element_key)
}

/// `dmp_factor_list` over `ZZ_I` (`ring`) or `QQ_I`, from the irreducible
/// monic factorization over `QQ<I>` given by `field_factor`.
pub fn factor_list(
    f: &MPoly<AN>,
    ring: bool,
    field_factor: impl Fn(&MPoly<AN>) -> PResult<(AN, Vec<(MPoly<AN>, u32)>)>,
) -> PResult<(AN, Vec<(MPoly<AN>, u32)>)> {
    let _scope = AlgebraicScope::new(field());
    let n = f.nvars;
    if f.is_zero() {
        return Ok((AN::zero(), Vec::new()));
    }
    let j = f.min_exps();
    let f = f.div_monomial(&j);
    let (content, f) = ground_primitive(&f, ring);
    let (mut coeff, factors) = field_factor(&f)?;
    let mut out = Vec::with_capacity(factors.len() + n);
    for (factor, i) in factors {
        if !ring {
            out.push((factor, i));
            continue;
        }
        // dup_zz_i_factor
        let (common, numerator) = clear_denoms(&factor);
        let (c, primitive) = ground_primitive(&numerator, true);
        let common = AN::rational(Rat::from(common), field());
        coeff = coeff.mul(&c.pow(i as i64)).div(&common.pow(i as i64));
        out.push((primitive, i));
    }
    for (v, &k) in j.iter().enumerate().rev() {
        if k > 0 {
            out.insert(0, (MPoly::var(v, n), k));
        }
    }
    if ring && !is_integral(&coeff) {
        return Err(PolyError::Failed(
            "Gaussian integer factorization produced a nonintegral unit".into(),
        ));
    }
    sort_factors(&mut out);
    Ok((coeff.mul(&content), out))
}

/// The domain the reference constructs without options for expressions
/// whose only nonrational coefficient atom is `I`: `ZZ_I` (`integral`) or
/// `QQ_I`, with the generators `gens` or, when none are given, the inferred
/// generators other than `I`. Floating coefficients select `CC` instead and
/// any other coefficient atom a composite domain, so both give `None`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DefaultSelection {
    pub integral: bool,
    pub gens: Vec<Expr>,
}

/// Whether the atom `I` occurs in `e`, visiting each shared node once so that
/// deeply shared expression graphs stay linear. `args()` builds the terms of
/// a sum and the factors of a product afresh, so each visited node is kept
/// alive in `seen`: a freed node's address could otherwise be reused by a
/// later node, which would then be skipped (`I` missed in `12 - 2*I`).
fn contains_imaginary_unit(e: &Expr) -> bool {
    let i = Expr::i();
    if e.has_matrix() {
        return crate::traverse::has(e, &i);
    }
    let mut seen = rustc_hash::FxHashMap::default();
    let mut pending = vec![e.clone()];
    while let Some(node) = pending.pop() {
        if node == i {
            return true;
        }
        if seen.insert(node.dag_id(), node.clone()).is_some() || super::rootof::is_rootof(&node) {
            continue;
        }
        pending.extend(node.args());
    }
    false
}

pub fn default_selection(exprs: &[Expr], gens: &[Expr]) -> Option<DefaultSelection> {
    let i = Expr::i();
    if !exprs.iter().any(contains_imaginary_unit) {
        return None;
    }
    let gens = if gens.is_empty() {
        let inferred = super::convert::from_exprs(exprs).ok()?;
        if inferred.float || !inferred.gens.contains(&i) {
            return None;
        }
        let gens: Vec<Expr> = inferred.gens.into_iter().filter(|g| g != &i).collect();
        if gens.is_empty() {
            return None;
        }
        gens
    } else {
        gens.to_vec()
    };
    let converted = super::convert::from_exprs_gens(exprs, &gens).ok()?;
    if converted.float || converted.gens[converted.ngens..] != [i] {
        return None;
    }
    let integral = converted.polys.iter().all(MPoly::is_integral);
    Some(DefaultSelection { integral, gens })
}

/// The reference `dmp_cancel` over `ZZ_I` (`ring`) or `QQ_I`:
/// `(cp, cq, p, q)` with `f/g = (cp*p)/(cq*q)`.
pub fn cancel(
    f: &MPoly<AN>,
    g: &MPoly<AN>,
    ring: bool,
) -> PResult<(AN, AN, MPoly<AN>, MPoly<AN>)> {
    let _scope = AlgebraicScope::new(field());
    if ring {
        let (_, p, q) = inner_gcd_at(f, g, 0)?;
        return Ok((AN::one(), AN::one(), p, q));
    }
    let (cq, f) = clear_denoms(f);
    let (cp, g) = clear_denoms(g);
    let (_, p, q) = inner_gcd_at(&f, &g, 0)?;
    let cq = AN::rational(Rat::from(cq), field());
    let cp = AN::rational(Rat::from(cp), field());
    let h = gcd(&cp, &cq);
    let cp = divmod(&cp, &h).0;
    let cq = divmod(&cq, &h).0;
    Ok((cp, cq, p, q))
}

/// The reference `PolyElement.cancel` over `ZZ_I`/`QQ_I`, used by the
/// expression `cancel`: ring cofactors, the scalar denominators folded back
/// in, then both sides multiplied by the unit that makes the leading
/// coefficient of the denominator canonical.
pub fn expr_cancel(f: &MPoly<AN>, g: &MPoly<AN>) -> PResult<(MPoly<AN>, MPoly<AN>)> {
    if g.is_zero() {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let _scope = AlgebraicScope::new(field());
    if f.is_zero() {
        return Ok((f.clone(), MPoly::one(f.nvars)));
    }
    let (cp, cq, p, q) = cancel(f, g, false)?;
    let (p, q) = (p.scale(&cp), q.scale(&cq));
    let unit = canonical_unit(&q.lc());
    if unit.is_one() {
        Ok((p, q))
    } else {
        Ok((p.scale(&unit), q.scale(&unit)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_selection_finds_i_inside_an_unexpanded_coefficient() {
        // `args()` builds `9*t` and `-2*I` afresh; with the visited set keyed
        // by addresses only, a reused address hid `I` for some names.
        for name in ["a", "b", "h", "i", "l", "m", "o", "p", "t", "w", "y"] {
            let e = crate::parser::parse_str(&format!("{name}**2*(12 - 2*I) + 9*{name}")).unwrap();
            let t = Expr::symbol(name);
            assert!(contains_imaginary_unit(&e), "{name}");
            let selected = default_selection(std::slice::from_ref(&e), std::slice::from_ref(&t));
            assert_eq!(selected, Some(DefaultSelection { integral: true, gens: vec![t] }), "{name}");
        }
    }

    #[test]
    fn expression_cancel_normalizes_the_denominator() {
        let _scope = AlgebraicScope::new(field());
        let x = MPoly::<AN>::var(0, 1);
        let i = constant(g(0, 1), 1);
        let one = MPoly::one(1);
        // (x**2 + 1)/(I*x + 1) = -I*x + 1
        let (p, q) = expr_cancel(&x.mul(&x).add(&one), &i.mul(&x).add(&one)).unwrap();
        assert_eq!(q, one);
        assert_eq!(p, constant(g(0, -1), 1).mul(&x).add(&one));
        // the Poly-level dmp_cancel keeps (1, x + I, I)
        let (_, _, p, q) = cancel(&x.mul(&x).add(&one), &i.mul(&x).add(&one), true).unwrap();
        assert_eq!(p, x.add(&i));
        assert_eq!(q, i);
    }

    fn g(x: i64, y: i64) -> AN {
        make(Rat::from_i64(x), Rat::from_i64(y))
    }

    #[test]
    fn element_normalization_and_euclid() {
        let _scope = AlgebraicScope::new(field());
        assert_eq!(normalize(&g(0, 2)), g(2, 0));
        assert_eq!(normalize(&g(-1, -1)), g(1, 1));
        assert_eq!(normalize(&g(1, -1)), g(1, 1));
        assert_eq!(gcd(&g(2, 0), &g(0, 2)), g(2, 0));
        assert_eq!(gcd(&g(3, 1), &g(2, 0)), g(1, 1));
        let (q, r) = divmod(&g(7, 3), &g(2, -1));
        assert_eq!(q.mul(&g(2, -1)).add(&r), g(7, 3));
        let (rx, ry) = parts(&r);
        assert!(rx.mul(&rx).add(&ry.mul(&ry)) < Rat::from_i64(5));
    }

    #[test]
    fn ring_gcd_lcm_and_factor_units() {
        let _scope = AlgebraicScope::new(field());
        let x = MPoly::<AN>::var(0, 1);
        let two = constant(g(2, 0), 1);
        let i = constant(g(0, 1), 1);
        // 2*x**2 + 2 and 2*x - 2*I
        let f = x.mul(&x).add(&MPoly::one(1)).mul(&two);
        let h = x.sub(&i).mul(&two);
        assert_eq!(ring_gcd(&f, &h).unwrap(), h);
        assert_eq!(ring_lcm(&f, &h).unwrap(), f);
        // gcd(0, g) = -g as in the reference
        assert_eq!(ring_gcd(&MPoly::zero(1), &h).unwrap(), h.neg());
        let (q, r) = ring_div(&x.mul(&x).add(&MPoly::one(1)), &x.scale(&g(2, 0)).sub(&i)).unwrap();
        assert!(q.is_zero());
        assert_eq!(r, x.mul(&x).add(&MPoly::one(1)));
    }
}
