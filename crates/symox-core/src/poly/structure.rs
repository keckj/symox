//! Sparse structural polynomial operations. These map exponent vectors only:
//! coefficients and exact coefficient descriptors are carried unchanged.
use super::convert::{PResult, PolyError};
use super::fraction::Polynomial;
use super::int::Coeff;
use super::mpoly::MPoly;
use crate::Expr;

const MAX_SLOTS: usize = 1_000_000;
const MAX_STORAGE: u64 = 64_000_000;
pub const MAX_DENSE: usize = 100_000;
fn resource(reason: &str) -> PolyError {
    PolyError::Unsupported(format!("polynomial structure resource limit: {reason}"))
}
fn check<C: Coeff>(p: &MPoly<C>) -> PResult<()> {
    if p.nvars > MAX_SLOTS
        || p.len()
            .checked_mul(p.nvars.max(1))
            .is_none_or(|n| n > MAX_SLOTS)
    {
        return Err(resource("exponent slots"));
    }
    let bits = p
        .coeffs()
        .iter()
        .map(|c| c.budget_weight().saturating_add(32 * p.nvars as u64))
        .fold(0u64, u64::saturating_add);
    if bits > MAX_STORAGE {
        return Err(resource("coefficient storage"));
    }
    Ok(())
}
fn with_support<T>(p: &Polynomial, f: impl FnOnce(Vec<Vec<u32>>) -> T) -> PResult<T> {
    let support = if let Some(a) = &p.expression {
        expression_check(a)?;
        a.terms().into_iter().map(|(e, _)| e).collect()
    } else if let Some(a) = &p.algebraic {
        check(&a.rep)?;
        a.rep.iter().map(|(e, _)| e.to_vec()).collect()
    } else if let Some(a) = &p.field {
        check(&a.rep)?;
        a.rep.iter().map(|(e, _)| e.to_vec()).collect()
    } else {
        check(&p.base.rep)?;
        p.base
            .rep
            .iter()
            .map(|(e, _)| e[..p.ngens].to_vec())
            .collect()
    };
    Ok(f(support))
}
fn expression_check(p: &super::expression_domain::ExpressionPoly) -> PResult<()> {
    let n = p.terms().len();
    if n.checked_mul(p.gens.len().max(1)).is_none_or(|k| k > MAX_SLOTS) {
        return Err(resource("exponent slots"));
    }
    Ok(())
}
fn changed<C: Coeff>(
    p: &MPoly<C>,
    main: usize,
    indices: &[usize],
    map: &impl Fn(&[u32]) -> Vec<u32>,
    inverse: &impl Fn(&[u32]) -> Vec<u32>,
) -> PResult<MPoly<C>> {
    check(p)?;
    let nvars = indices.len() + p.nvars - main;
    if p.len()
        .checked_mul(nvars.max(1))
        .is_none_or(|n| n > MAX_SLOTS)
    {
        return Err(resource("output exponent slots"));
    }
    let mut output = Vec::with_capacity(p.len());
    for (e, c) in p.iter() {
        let transformed = map(&e[..main]);
        if transformed.len() != indices.len() || inverse(&transformed) != e[..main] {
            return Err(PolyError::Failed(
                "structural transform failed its exponent reconstruction certificate".into(),
            ));
        }
        let mut powers = transformed;
        powers.extend_from_slice(&e[main..]);
        output.push((powers, c.clone()));
    }
    // An inverse exists on this support, so no terms can merge or disappear.
    let result = MPoly::from_terms(nvars, output);
    if result.len() != p.len() {
        return Err(PolyError::Failed(
            "structural transform changed its sparse support cardinality".into(),
        ));
    }
    Ok(result)
}
fn transform(
    p: &Polynomial,
    indices: &[usize],
    map: impl Fn(&[u32]) -> Vec<u32>,
    inverse: impl Fn(&[u32]) -> Vec<u32>,
) -> PResult<Polynomial> {
    if indices.is_empty() || indices.iter().any(|&i| i >= p.ngens) {
        return Err(PolyError::Failed(
            "invalid structural generator selection".into(),
        ));
    }
    let gens: Vec<_> = indices.iter().map(|&i| p.true_gens()[i].clone()).collect();
    if let Some(a) = &p.expression {
        expression_check(a)?;
        let terms = a.terms();
        let mut output = Vec::with_capacity(terms.len());
        for (e, c) in terms {
            let t = map(&e);
            if t.len() != indices.len() || inverse(&t) != e {
                return Err(PolyError::Failed(
                    "structural transform failed its exponent reconstruction certificate".into(),
                ));
            }
            output.push((t, c));
        }
        let mut result=super::expression_domain::ExpressionPoly::from_terms(output, &gens)?;
        result.ground=a.ground.clone();
        return Ok(result.into());
    }
    if let Some(a) = &p.algebraic {
        let mut result = a.with_rep(changed(&a.rep, a.gens.len(), indices, &map, &inverse)?);
        result.gens = gens;
        return Ok(result.into());
    }
    if let Some(a) = &p.field {
        let mut result = a.with_rep(changed(&a.rep, a.gens.len(), indices, &map, &inverse)?);
        result.gens = gens;
        return Ok(result.into());
    }
    let mut result = p
        .base
        .with_rep(changed(&p.base.rep, p.ngens, indices, &map, &inverse)?);
    result.gens = gens;
    result.gens.extend_from_slice(p.domain_gens());
    result.ngens = indices.len();
    Ok(result.into())
}
fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
pub fn deflate(p: &Polynomial) -> PResult<(Vec<u32>, Polynomial)> {
    let divisor = with_support(p, |support| {
        let mut gcds = vec![0; p.ngens];
        for e in support {
            for (g, &k) in gcds.iter_mut().zip(&e) {
                *g = gcd(*g, k);
            }
        }
        for g in &mut gcds {
            *g = (*g).max(1);
        }
        gcds
    })?;
    let indices: Vec<_> = (0..p.ngens).collect();
    let result = transform(
        p,
        &indices,
        |e| e.iter().zip(&divisor).map(|(k, d)| k / d).collect(),
        |e| e.iter().zip(&divisor).map(|(k, d)| k * d).collect(),
    )?;
    Ok((divisor, result))
}
pub fn terms_gcd(p: &Polynomial) -> PResult<(Vec<u32>, Polynomial)> {
    let minimum = with_support(p, |support| {
        let mut mins = support.first().cloned().unwrap_or_else(|| vec![0; p.ngens]);
        for e in support {
            for (g, &k) in mins.iter_mut().zip(&e) {
                *g = (*g).min(k);
            }
        }
        mins
    })?;
    let indices: Vec<_> = (0..p.ngens).collect();
    let result = transform(
        p,
        &indices,
        |e| e.iter().zip(&minimum).map(|(k, d)| k - d).collect(),
        |e| e.iter().zip(&minimum).map(|(k, d)| k + d).collect(),
    )?;
    Ok((minimum, result))
}
pub fn exclude(p: &Polynomial) -> PResult<Polynomial> {
    let needed = with_support(p, |support| {
        let mut needed = vec![false; p.ngens];
        for e in support {
            for (b, k) in needed.iter_mut().zip(e) {
                *b |= k != 0;
            }
        }
        needed
    })?;
    let indices: Vec<_> = (0..p.ngens).filter(|&i| needed[i]).collect();
    // SymPy retains the original generator tuple for all ground polynomials,
    // including zero: a public Poly never becomes a generator-free scalar.
    if indices.is_empty() || indices.len() == p.ngens {
        return Ok(p.clone());
    }
    transform(
        p,
        &indices,
        |e| indices.iter().map(|&i| e[i]).collect(),
        |e| {
            let mut old = vec![0; p.ngens];
            for (&i, &k) in indices.iter().zip(e) {
                old[i] = k;
            }
            old
        },
    )
}
pub fn reorder(p: &Polynomial, indices: &[usize]) -> PResult<Polynomial> {
    if indices.len() != p.ngens
        || indices
            .iter()
            .enumerate()
            .any(|(j, &i)| i >= p.ngens || indices[..j].contains(&i))
    {
        return Err(PolyError::Failed(
            "generators list can differ only up to order of elements".into(),
        ));
    }
    transform(
        p,
        indices,
        |e| indices.iter().map(|&i| e[i]).collect(),
        |e| {
            let mut old = vec![0; p.ngens];
            for (&i, &k) in indices.iter().zip(e) {
                old[i] = k;
            }
            old
        },
    )
}
fn remapped<C: Coeff>(
    p: &MPoly<C>,
    main: usize,
    width: usize,
    map: &impl Fn(&[u32]) -> Option<Vec<u32>>,
) -> PResult<MPoly<C>> {
    check(p)?;
    let nvars = width + p.nvars - main;
    if p.len()
        .checked_mul(nvars.max(1))
        .is_none_or(|n| n > MAX_SLOTS)
    {
        return Err(resource("output exponent slots"));
    }
    let mut output = Vec::with_capacity(p.len());
    for (e, c) in p.iter() {
        if let Some(mut powers) = map(&e[..main]) {
            debug_assert_eq!(powers.len(), width);
            powers.extend_from_slice(&e[main..]);
            output.push((powers, c.clone()));
        }
    }
    Ok(MPoly::from_terms(nvars, output))
}
/// Rebuild `p` over the generators `gens`, mapping each exponent vector of
/// the true generators through `map` (`None` drops the term). Terms that
/// land on the same exponent are added in the coefficient domain.
fn remap(
    p: &Polynomial,
    gens: Vec<Expr>,
    map: impl Fn(&[u32]) -> Option<Vec<u32>>,
) -> PResult<Polynomial> {
    let width = gens.len();
    if width == 0 {
        return Err(PolyError::Failed("a polynomial needs at least one generator".into()));
    }
    if let Some(a) = &p.expression {
        expression_check(a)?;
        let terms = a.terms().into_iter().filter_map(|(e, c)| map(&e).map(|t| (t, c))).collect();
        let mut result = super::expression_domain::ExpressionPoly::from_summed_terms(terms, &gens)?;
        result.ground = a.ground.clone();
        return Ok(result.into());
    }
    if let Some(a) = &p.algebraic {
        let mut result = a.with_rep(remapped(&a.rep, a.gens.len(), width, &map)?);
        result.gens = gens;
        return Ok(result.into());
    }
    if let Some(a) = &p.field {
        let mut result = a.with_rep(remapped(&a.rep, a.gens.len(), width, &map)?);
        result.gens = gens;
        return Ok(result.into());
    }
    let mut result = p
        .base
        .with_rep(remapped(&p.base.rep, p.ngens, width, &map)?);
    result.gens = gens;
    result.gens.extend_from_slice(p.domain_gens());
    result.ngens = width;
    Ok(result.into())
}
/// SymPy's `Poly.homogenize(s)`: every term is completed to the total degree
/// with powers of `s`, which is appended unless it already is a generator.
pub fn homogenize(p: &Polynomial, s: &Expr) -> PResult<Polynomial> {
    if p.domain_gens().contains(s)
        || p.field.as_ref().is_some_and(|f| f.parameters.contains(s))
    {
        return Err(PolyError::Failed(
            "homogenizing symbol occurs in the coefficient domain".into(),
        ));
    }
    let degree = with_support(p, |support| {
        support
            .iter()
            .map(|e| e.iter().map(|&k| k as u64).sum::<u64>())
            .max()
            .unwrap_or(0)
    })?;
    let mut gens = p.true_gens().to_vec();
    let index = match gens.iter().position(|g| g == s) {
        Some(i) => i,
        None => {
            gens.push(s.clone());
            gens.len() - 1
        }
    };
    // The total degree may exceed one exponent's storage width. Only the
    // new coordinate must fit: x**N*y**N + x**N can homogenize to
    // x**N*y**N + x**N*z**N when 2*N exceeds u32.
    let fits = with_support(p, |support| support.iter().all(|e| {
        let missing = degree - e.iter().map(|&k| u64::from(k)).sum::<u64>();
        u64::from(e.get(index).copied().unwrap_or(0)) + missing <= u64::from(u32::MAX)
    }))?;
    if !fits {
        return Err(super::mpoly::exponent_limit_error());
    }
    let width = gens.len();
    // The support preflight proves each resulting coordinate fits in u32.
    remap(p, gens, |e| {
        let mut powers = e.to_vec();
        powers.resize(width, 0);
        let missing = degree - e.iter().map(|&k| k as u64).sum::<u64>();
        powers[index] = (powers[index] as u64 + missing) as u32;
        Some(powers)
    })
}
/// SymPy's `Poly.slice(gen, m, n)`: the terms whose degree in generator
/// `index` lies in `[m, n)`.
pub fn slice(p: &Polynomial, index: usize, m: u32, n: u32) -> PResult<Polynomial> {
    if index >= p.ngens {
        return Err(PolyError::Failed("invalid polynomial generator index".into()));
    }
    remap(p, p.true_gens().to_vec(), |e| {
        (m <= e[index] && e[index] < n).then(|| e.to_vec())
    })
}
/// SymPy's `Poly.ltrim(gen)`: drop the generators before `index`, which must
/// not occur in any term. `None` reports that the polynomial cannot be trimmed.
pub fn ltrim(p: &Polynomial, index: usize) -> PResult<Option<Polynomial>> {
    if index >= p.ngens {
        return Err(PolyError::Failed("invalid polynomial generator index".into()));
    }
    if index == 0 {
        return Ok(Some(p.clone()));
    }
    let free = with_support(p, |support| {
        support.iter().all(|e| e[..index].iter().all(|&k| k == 0))
    })?;
    if !free {
        return Ok(None);
    }
    remap(p, p.true_gens()[index..].to_vec(), |e| Some(e[index..].to_vec())).map(Some)
}
#[derive(Debug, PartialEq, Eq)]
pub enum HomogeneousOrder {
    Zero,
    Degree(u64),
    Inhomogeneous,
}
pub fn homogeneous_order(p: &Polynomial) -> PResult<HomogeneousOrder> {
    with_support(p, |support| {
        let mut order = None;
        for e in support {
            let degree = e.iter().map(|&k| k as u64).sum::<u64>();
            if order.is_some_and(|k| k != degree) {
                return HomogeneousOrder::Inhomogeneous;
            }
            order = Some(degree);
        }
        order.map_or(HomogeneousOrder::Zero, HomogeneousOrder::Degree)
    })
}
pub fn has_only_indices(p: &Polynomial, indices: &[usize]) -> PResult<bool> {
    if indices.iter().any(|&i| i >= p.ngens) {
        return Err(PolyError::Failed(
            "invalid polynomial generator index".into(),
        ));
    }
    with_support(p, |support| {
        support.iter().all(|e| {
            e.iter()
                .enumerate()
                .all(|(i, &k)| k == 0 || indices.contains(&i))
        })
    })
}
/// The three dense views share one guard and one linear sparse-to-dense pass.
pub fn dense_terms(p: &Polynomial) -> PResult<Vec<(Vec<u32>, Expr)>> {
    if p.ngens != 1 {
        return Err(PolyError::Failed(
            "multivariate polynomials not supported".into(),
        ));
    }
    with_support(p, |_| ())?;
    let degree = if let Some(f) = &p.expression {
        f.degree(0)
    } else if let Some(f) = &p.field {
        f.rep.degree(0)
    } else if let Some(f) = &p.algebraic {
        f.rep.degree(0)
    } else {
        p.base.rep.degree(0)
    }
    .max(0) as usize;
    let size = degree
        .checked_add(1)
        .filter(|&n| n <= MAX_DENSE)
        .ok_or_else(|| resource("dense output length"))?;
    let terms = p.terms();
    let mut output: Vec<_> = (0..size)
        .rev()
        .map(|k| (vec![k as u32], Expr::zero()))
        .collect();
    for (m, c) in terms {
        output[degree - m[0] as usize].1 = c;
    }
    Ok(output)
}

/// The same exponent vectors and coefficients over new generators, as
/// SymPy's `Poly.new(rep, *gens)`. The generators must be distinct and must
/// not occur in the coefficient domain. Generators that are not symbols are
/// supported over the plain `ZZ`/`QQ`/`RR` (and `K[params]`) backends, whose
/// storage does not depend on re-reading the generators from expressions.
pub fn rename_gens(p: &Polynomial, gens: &[Expr]) -> PResult<Polynomial> {
    if gens.len() != p.ngens {
        return Err(PolyError::Failed(format!(
            "expected {} generators, got {}",
            p.ngens,
            gens.len()
        )));
    }
    if gens.iter().enumerate().any(|(j, g)| gens[..j].contains(g)) {
        return Err(PolyError::Failed("generators must be distinct expressions".into()));
    }
    let parameters: &[Expr] = if let Some(f) = &p.field {
        &f.parameters
    } else if p.algebraic.is_some() || p.expression.is_some() {
        &[]
    } else {
        p.domain_gens()
    };
    if gens.iter().any(|g| parameters.contains(g)) {
        return Err(PolyError::Failed("a generator also occurs in the coefficient domain".into()));
    }
    // Generator descriptors label exponent-vector slots. Renaming does not
    // reinterpret coefficients or require those labels to be bare symbols.
    remap(p, gens.to_vec(), |e| Some(e.to_vec()))
}

/// SymPy's `Poly.as_expr()`: every term is `Mul(coefficient, monomial)`, so
/// composite coefficients of a `K[params]` domain stay grouped.
pub fn sympy_expr(p: &Polynomial) -> Expr {
    if p.field.is_some() || p.algebraic.is_some() || p.expression.is_some() || !p.base.has_domain_gens() {
        return p.as_expr();
    }
    let gens = p.true_gens().to_vec();
    crate::ops::add(
        &p.terms()
            .into_iter()
            .map(|(e, c)| {
                let mut factors = vec![c];
                for (g, &k) in gens.iter().zip(&e) {
                    if k != 0 {
                        factors.push(crate::ops::pow(g, &Expr::int(k as i64)));
                    }
                }
                crate::ops::mul(&factors)
            })
            .collect::<Vec<_>>(),
    )
}

/// SymPy's `dup_cyclotomic_p` on descending integer coefficients: Graeffe's
/// root squaring `F(x**2) = +-f(x) f(-x)` maps a cyclotomic polynomial to
/// itself, to its mirror, or to a square of a cyclotomic polynomial.
pub fn cyclotomic_p(f: &[super::int::Int], irreducible: bool) -> bool {
    use super::int::{Int, Rat};
    fn strip(f: &[Int]) -> Vec<Int> {
        f.iter().skip_while(|c| c.is_zero()).cloned().collect()
    }
    fn sqr(f: &[Int]) -> Vec<Int> {
        let asc: Vec<Int> = f.iter().rev().cloned().collect();
        let mut r = super::dense::mul(&asc, &asc);
        r.reverse();
        r
    }
    fn sub(f: &[Int], g: &[Int]) -> Vec<Int> {
        let n = f.len().max(g.len());
        let at = |p: &[Int], i: usize| if i + p.len() >= n { p[i + p.len() - n].clone() } else { Int::zero() };
        strip(&(0..n).map(|i| at(f, i).sub(&at(g, i))).collect::<Vec<_>>())
    }
    fn positive(f: Vec<Int>) -> Vec<Int> {
        if f.first().is_some_and(|c| c.is_negative()) { f.iter().map(|c| c.neg()).collect() } else { f }
    }
    let (Some(lc), Some(tc)) = (f.first(), f.last()) else { return false };
    if !lc.is_one() || !(tc.is_one() || tc.is_minus_one()) {
        return false;
    }
    let n = f.len() - 1;
    if !irreducible {
        let asc: Vec<Int> = f.iter().rev().cloned().collect();
        let (c, factors) = super::factor_multi::factor_list_mz(&super::dense::to_mpoly(&asc, 0, 1));
        if !c.is_one() || factors.len() != 1 || factors[0].1 != 1
            || super::dense::from_mpoly(&factors[0].0, 0) != asc
        {
            return false;
        }
    }
    let even: Vec<Int> = f.iter().skip(n % 2).step_by(2).cloned().collect();
    let odd: Vec<Int> = if n >= 1 { f.iter().skip((n - 1) % 2).step_by(2).cloned().collect() } else { vec![] };
    let g = sqr(&strip(&even));
    let mut h = sqr(&strip(&odd));
    if !h.is_empty() {
        h.push(Int::zero());
    }
    let big = positive(sub(&g, &h));
    if big == f {
        return true;
    }
    let mirror = positive(f.iter().enumerate().map(|(i, c)| if (n - i) % 2 == 1 { c.neg() } else { c.clone() }).collect());
    if big == mirror && cyclotomic_p(&mirror, false) {
        return true;
    }
    if big.is_empty() {
        return false;
    }
    // Primitive square-free part with a positive leading coefficient.
    let asc: Vec<Rat> = big.iter().rev().map(|c| Rat::from(c.clone())).collect();
    let gcd = super::dense::gcd_field(&asc, &super::dense::derivative(&asc));
    let part = super::dense::divrem(&asc, &gcd).0;
    let (_, ints) = super::dense::clear_denoms(&part);
    let mut root: Vec<Int> = super::dense::primitive(&ints).1;
    root.reverse();
    sqr(&root) == big && cyclotomic_p(&root, false)
}

#[cfg(test)]
mod tests {
    use super::super::polyobj::Poly;
    use super::*;
    use crate::parser::parse_str;
    fn p(s: &str) -> Polynomial {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        Poly::from_expr(&parse_str(s).unwrap(), &[x, y])
            .unwrap()
            .into()
    }
    #[test]
    fn exponent_maps_reconstruct_and_preserve_coefficients() {
        let input = p("3*x^6*y^4+5*x^2*y^2");
        let (d, f) = deflate(&input).unwrap();
        assert_eq!(d, vec![2, 2]);
        assert_eq!(f.as_expr(), p("3*x^3*y^2+5*x*y").as_expr());
        let (d, f) = terms_gcd(&input).unwrap();
        assert_eq!(d, vec![2, 2]);
        assert_eq!(f.as_expr(), p("3*x^4*y^2+5").as_expr());
        let swapped = reorder(&input, &[1, 0]).unwrap();
        assert_eq!(swapped.as_expr(), input.as_expr());
        assert_eq!(reorder(&swapped, &[1, 0]).unwrap().rep, input.rep);
        assert_eq!(
            homogeneous_order(&p("x^2+y^2")).unwrap(),
            HomogeneousOrder::Degree(2)
        );
        assert_eq!(homogeneous_order(&p("0")).unwrap(), HomogeneousOrder::Zero);
    }
    #[test]
    fn homogenize_slice_and_ltrim_follow_sympy() {
        let z = Expr::symbol("z");
        let y = Expr::symbol("y");
        let h = homogenize(&p("x^3+x*y+1"), &z).unwrap();
        assert_eq!(h.true_gens().len(), 3);
        let x = Expr::symbol("x");
        let expected: Polynomial = Poly::from_expr(&parse_str("x^3+x*y*z+z^3").unwrap(), &[x.clone(), y.clone(), z])
            .unwrap()
            .into();
        assert_eq!(h.rep, expected.rep);
        // An existing generator absorbs the missing degree and merges terms.
        let merged = homogenize(&p("x*y+x"), &y).unwrap();
        assert_eq!(merged.as_expr(), p("2*x*y").as_expr());
        let s = slice(&p("x^3*y+x^2+x+y"), 0, 1, 3).unwrap();
        assert_eq!(s.as_expr(), p("x^2+x").as_expr());
        assert!(ltrim(&p("x^2*y"), 1).unwrap().is_none());
        let t = ltrim(&p("y^2+1"), 1).unwrap().unwrap();
        assert_eq!(t.true_gens(), &[y]);
    }
    #[test]
    fn sparse_large_exponents_do_not_trigger_dense_allocations() {
        let x = Expr::symbol("x");
        let p: Polynomial = Poly::from_expr(&parse_str("x^1000000000+1").unwrap(), &[x])
            .unwrap()
            .into();
        assert!(dense_terms(&p)
            .unwrap_err()
            .to_string()
            .contains("dense output length"));
        let (d, f) = deflate(&p).unwrap();
        assert_eq!(d, vec![1_000_000_000]);
        assert_eq!(f.rep.degree(0), 1);
    }

    #[test]
    fn cyclotomic_recognition_matches_sympy() {
        let v = |c: &[i64]| c.iter().map(|&k| super::super::int::Int::from_i128(k as i128)).collect::<Vec<_>>();
        assert!(super::cyclotomic_p(&v(&[1, 1, 1]), false));
        assert!(super::cyclotomic_p(&v(&[1, 0, 0, 0, 1]), false));
        assert!(super::cyclotomic_p(&v(&[1, -1, 0, 1, -1, 1, 0, -1, 1]), false)); // Phi_15
        assert!(!super::cyclotomic_p(&v(&[1, 0, 1, 0, 1]), false)); // Phi_3 * Phi_6
        assert!(!super::cyclotomic_p(&v(&[1, 2, 1]), false));
        assert!(!super::cyclotomic_p(&v(&[2, 1]), false));
    }
}

#[cfg(test)]
mod wide_total_degree_tests {
    use super::*;
    use super::super::{polyobj::{Poly,Domain},int::Rat};
    #[test]
    fn wide_total_degree_homogenization_checks_coordinates_only() {
        let n=1u32<<31;
        let x=Expr::symbol("x");let y=Expr::symbol("y");let z=Expr::symbol("z");
        let p:Polynomial=Poly::new(MPoly::from_terms(2,vec![(vec![n,n],Rat::one()),(vec![n,0],Rat::one())]),vec![x,y.clone()],2,Domain::ZZ).into();
        let h=homogenize(&p,&z).unwrap();
        assert_eq!(h.monoms(),vec![vec![n,n,0],vec![n,0,n]]);
        assert_eq!(homogeneous_order(&h).unwrap(),HomogeneousOrder::Degree(1u64<<32));
        assert_eq!(homogenize(&p,&y).unwrap().monoms(),vec![vec![n,n]]);
        let q=p.base.with_rep(p.base.rep.add(&MPoly::one(2)));
        assert!(homogenize(&q.into(),&z).is_err());
    }
}
