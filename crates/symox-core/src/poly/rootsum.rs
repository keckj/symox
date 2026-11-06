//! The rational case of the reference `RootSum` constructor
//! (`sympy.polys.rootoftools.RootSum._rational_case`).
//!
//! For a rational function `expr` of `var` and a polynomial of degree `n`,
//! the sum of `expr` over the roots is rewritten in the coefficients of the
//! polynomial: the roots become the symbols `r0, ..., r(n-1)`, the summed
//! fraction is brought over a common denominator, every coefficient of the
//! numerator and of the denominator (as polynomials in the remaining
//! generators over `QQ[r0, ...]`) is reduced by the elementary symmetric
//! polynomials (`PolyElement.symmetrize`, the remainder being dropped as in
//! the reference) and Vieta's formulas replace the elementary symmetric
//! polynomials. The result is `factor(p/q)`.
use std::collections::{BTreeMap, HashMap};

use super::convert::from_exprs;
use super::int::Rat;
use crate::expr::Expr;
use crate::ops::{add, div, mul, pow};

type Sparse = BTreeMap<Vec<u32>, Rat>;

fn sparse_mul(a: &Sparse, b: &Sparse) -> Sparse {
    let mut out = Sparse::new();
    for (ma, ca) in a {
        for (mb, cb) in b {
            let m: Vec<u32> = ma.iter().zip(mb).map(|(x, y)| x + y).collect();
            let c = ca.mul(cb);
            let entry = out.entry(m).or_insert_with(Rat::zero);
            *entry = entry.add(&c);
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

fn sparse_sub_assign(a: &mut Sparse, b: &Sparse) {
    for (m, c) in b {
        let entry = a.entry(m.clone()).or_insert_with(Rat::zero);
        *entry = entry.sub(c);
    }
    a.retain(|_, c| !c.is_zero());
}

/// The elementary symmetric polynomial `e_k` in `n` variables.
fn elementary(k: usize, n: usize) -> Sparse {
    let mut out = Sparse::new();
    fn walk(start: usize, left: usize, n: usize, m: &mut Vec<u32>, out: &mut Sparse) {
        if left == 0 {
            out.insert(m.clone(), Rat::one());
            return;
        }
        for i in start..n {
            m[i] = 1;
            walk(i + 1, left - 1, n, m, out);
            m[i] = 0;
        }
    }
    walk(0, k, n, &mut vec![0; n], &mut out);
    out
}

/// `PolyElement.symmetrize` on a polynomial in `n` variables: the terms
/// `(coeff, exponents)` of the symmetric part in the elementary symmetric
/// polynomials, in the order of the reduction. The remainder is dropped.
fn symmetrize(f: &Sparse, n: usize, elementary_polys: &[Sparse], powers: &mut HashMap<(usize, u32), Sparse>) -> Vec<(Rat, Vec<u32>)> {
    let mut f = f.clone();
    let weights: Vec<u32> = (1..=n as u32).rev().collect();
    let mut symmetric = Vec::new();
    while !f.is_empty() {
        let mut best: Option<(u32, Vec<u32>, Rat)> = None;
        // `f.terms()` lists the terms in decreasing lex order.
        for (monom, coeff) in f.iter().rev() {
            if (0..n.saturating_sub(1)).all(|i| monom[i] >= monom[i + 1]) {
                let height = weights.iter().zip(monom).map(|(w, m)| w * m).max().unwrap_or(0);
                if best.as_ref().map_or(true, |b| height > b.0) {
                    best = Some((height, monom.clone(), coeff.clone()));
                }
            }
        }
        let Some((_, monom, coeff)) = best else { break };
        let exponents: Vec<u32> = (0..n).map(|i| monom[i] - if i + 1 < n { monom[i + 1] } else { 0 }).collect();
        let mut product = Sparse::new();
        product.insert(vec![0; n], coeff.clone());
        for (i, &k) in exponents.iter().enumerate() {
            if k > 0 {
                let power = powers.entry((i, k)).or_insert_with(|| {
                    let mut acc = Sparse::new();
                    acc.insert(vec![0; n], Rat::one());
                    for _ in 0..k {
                        acc = sparse_mul(&acc, &elementary_polys[i]);
                    }
                    acc
                }).clone();
                product = sparse_mul(&product, &power);
            }
        }
        sparse_sub_assign(&mut f, &product);
        symmetric.push((coeff, exponents));
    }
    symmetric
}

/// The coefficients over `QQ[roots]` of `p`, a polynomial with generators
/// found automatically (`Poly(p, domain=QQ[roots], expand=False)`), as
/// `(monomial in the other generators, coefficient)`. `Err(())` when the
/// conversion would leave the reference domain (floats, generators built
/// from a root), `Ok(None)` on `GeneratorsNeeded`.
#[allow(clippy::type_complexity)]
fn coefficients(p: &Expr, roots: &[Expr]) -> Result<Option<(Vec<Expr>, Vec<(Vec<u32>, Sparse)>)>, ()> {
    let n = roots.len();
    let converted = from_exprs(std::slice::from_ref(p)).map_err(|_| ())?;
    if converted.float || converted.ngens != converted.gens.len() {
        return Err(());
    }
    let mut root_index = vec![None; converted.gens.len()];
    let mut others = Vec::new();
    for (j, g) in converted.gens.iter().enumerate() {
        if let Some(i) = roots.iter().position(|r| r == g) {
            root_index[j] = Some(i);
        } else {
            if roots.iter().any(|r| crate::traverse::has(g, r)) {
                return Err(());
            }
            others.push(j);
        }
    }
    let poly = &converted.polys[0];
    let mut rows: Vec<(Vec<u32>, Sparse)> = Vec::new();
    for (e, c) in poly.iter() {
        let outer: Vec<u32> = others.iter().map(|&j| e[j]).collect();
        let mut inner = vec![0u32; n];
        for (j, idx) in root_index.iter().enumerate() {
            if let Some(i) = idx {
                inner[*i] = e[j];
            }
        }
        let row = match rows.iter_mut().find(|r| r.0 == outer) {
            Some(r) => r,
            None => {
                rows.push((outer, Sparse::new()));
                rows.last_mut().unwrap()
            }
        };
        let entry = row.1.entry(inner).or_insert_with(Rat::zero);
        *entry = entry.add(c);
    }
    if others.is_empty() {
        return Ok(None);
    }
    let gens = others.iter().map(|&j| converted.gens[j].clone()).collect();
    Ok(Some((gens, rows)))
}

/// The whole polynomial over `QQ[roots]` (the `GeneratorsNeeded` branch).
fn whole(p: &Expr, roots: &[Expr]) -> Option<Sparse> {
    let n = roots.len();
    match coefficients(p, roots) {
        Ok(None) => {}
        _ => return None,
    }
    let converted = from_exprs(std::slice::from_ref(p)).ok()?;
    let mut out = Sparse::new();
    for (e, c) in converted.polys[0].iter() {
        let mut inner = vec![0u32; n];
        for (j, g) in converted.gens.iter().enumerate() {
            let i = roots.iter().position(|r| r == g)?;
            inner[i] = e[j];
        }
        let entry = out.entry(inner).or_insert_with(Rat::zero);
        *entry = entry.add(c);
    }
    out.retain(|_, c| !c.is_zero());
    Some(out)
}

fn monomial(gens: &[Expr], m: &[u32]) -> Expr {
    let mut factors = Vec::new();
    for (g, &k) in gens.iter().zip(m) {
        if k > 0 {
            factors.push(pow(g, &Expr::int(k as i64)));
        }
    }
    mul(&factors)
}

/// The reference `RootSum._rational_case(poly, Lambda(var, expr))` for a
/// polynomial with Vieta values `viete[i] = (-1)**(i+1) * c(i+1)/lc`
/// (`viete(poly, roots)`). `expand` is the public `Expr.expand`. `None`
/// when a step would leave the native representation (floats, a generator
/// built from a root, a root missing from every coefficient); the caller
/// then runs the reference body.
pub fn rational_case(expr: &Expr, var: &Expr, viete: &[Expr], expand: &dyn Fn(&Expr) -> Option<Expr>) -> Option<Expr> {
    let n = viete.len();
    if n == 0 {
        return None;
    }
    let roots: Vec<Expr> = (0..n).map(|i| Expr::symbol(&format!("r{i}"))).collect();
    let terms: Vec<Expr> = roots.iter().map(|r| crate::subs::subs(expr, &[(var.clone(), r.clone())])).collect();
    let f = add(&terms);
    let (p, q) = crate::simptools::as_numer_denom(&super::api::together(&f));
    let p = expand(&p)?;
    let q = expand(&q)?;
    // (gens, rows) or a single coefficient.
    let split = |e: &Expr| -> Option<(Option<Vec<Expr>>, Vec<(Vec<u32>, Sparse)>)> {
        match coefficients(e, &roots) {
            Ok(Some((gens, rows))) => Some((Some(gens), rows)),
            Ok(None) => Some((None, vec![(Vec::new(), whole(e, &roots)?)])),
            Err(()) => None,
        }
    };
    let (p_gens, p_rows) = split(&p)?;
    let (q_gens, q_rows) = split(&q)?;
    // `sring` over every coefficient: each root must be a generator.
    let mut present = vec![false; n];
    for (_, c) in p_rows.iter().chain(&q_rows) {
        for m in c.keys() {
            for (i, &k) in m.iter().enumerate() {
                if k > 0 {
                    present[i] = true;
                }
            }
        }
    }
    if present.iter().any(|x| !x) {
        return None;
    }
    let elementary_polys: Vec<Sparse> = (1..=n).map(|k| elementary(k, n)).collect();
    let mut powers = HashMap::new();
    let symbols: Vec<Expr> = (1..=n).map(|i| Expr::symbol(&format!("s{i}"))).collect();
    let values: Vec<(Expr, Expr)> = symbols.iter().cloned().zip(viete.iter().cloned()).collect();
    let mut reduce = |c: &Sparse| -> Expr {
        let sym = symmetrize(c, n, &elementary_polys, &mut powers);
        let parts: Vec<Expr> = sym.iter().map(|(coeff, ex)| mul(&[Expr::num(coeff.to_number()), monomial(&symbols, ex)])).collect();
        crate::subs::subs(&add(&parts), &values)
    };
    let mut rebuild = |gens: &Option<Vec<Expr>>, rows: &[(Vec<u32>, Sparse)]| -> Expr {
        match gens {
            None => reduce(&rows[0].1),
            Some(gens) => {
                let parts: Vec<Expr> = rows.iter().map(|(m, c)| mul(&[reduce(c), monomial(gens, m)])).collect();
                add(&parts)
            }
        }
    };
    let p = rebuild(&p_gens, &p_rows);
    let q = rebuild(&q_gens, &q_rows);
    super::api::factor(&div(&p, &q)).ok()
}
