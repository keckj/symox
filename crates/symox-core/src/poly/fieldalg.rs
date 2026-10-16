//! Ordered sparse algorithms over an arbitrary exact coefficient field.

use super::groebner::MonomialOrder;
use super::convert::{PResult, PolyError};
use super::int::Field;
use super::mpoly::MPoly;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;

fn divides(a: &[u32], b: &[u32]) -> bool {
    a.iter().zip(b).all(|(x, y)| x <= y)
}
fn difference(a: &[u32], b: &[u32]) -> Vec<u32> {
    a.iter().zip(b).map(|(x, y)| x - y).collect()
}

pub fn reduce<C: Field>(
    input: &MPoly<C>,
    basis: &[MPoly<C>],
    order: MonomialOrder,
) -> (Vec<MPoly<C>>, MPoly<C>) {
    reduce_filtered(input, basis, order, |_| true)
}

pub fn reduce_filtered<C: Field>(
    input: &MPoly<C>,
    basis: &[MPoly<C>],
    order: MonomialOrder,
    allow: impl Fn(&C) -> bool,
) -> (Vec<MPoly<C>>, MPoly<C>) {
    let mut f = input.clone();
    let mut q = vec![MPoly::zero(input.nvars); basis.len()];
    let mut r = MPoly::zero(input.nvars);
    let leads: Vec<_> = basis
        .iter()
        .map(|p| {
            if p.is_zero() {
                None
            } else {
                Some(order.leading_index(p))
            }
        })
        .collect();
    while !f.is_zero() {
        let i = order.leading_index(&f);
        let e = f.exp(i).to_vec();
        let c = f.coeff(i).clone();
        if let Some((j, g)) = basis.iter().enumerate().find(|(j, g)| {
            leads[*j].is_some_and(|k| divides(g.exp(k), &e) && allow(&c.div(g.coeff(k))))
        }) {
            let k = leads[j].unwrap();
            let t = MPoly::monomial(&difference(&e, g.exp(k)), c.div(g.coeff(k)));
            q[j] = q[j].add(&t);
            f = f.sub(&t.mul(g));
        } else {
            let t = MPoly::monomial(&e, c);
            r = r.add(&t);
            f = f.sub(&t);
        }
    }
    (q, r)
}

pub fn monic<C: Field>(p: &MPoly<C>, order: MonomialOrder) -> MPoly<C> {
    if p.is_zero() {
        p.clone()
    } else {
        p.scale(&p.coeff(order.leading_index(p)).inv())
    }
}

pub fn buchberger<C: Field>(input: &[MPoly<C>], order: MonomialOrder) -> Vec<MPoly<C>> {
    let mut basis = Vec::new();
    for f in input {
        let r = reduce(f, &basis, order).1;
        if r.is_zero() {
            continue;
        }
        if r.is_constant() {
            return vec![MPoly::one(r.nvars)];
        }
        basis.push(monic(&r, order));
    }
    let mut pairs = VecDeque::new();
    for j in 0..basis.len() {
        for i in 0..j {
            pairs.push_back((i, j));
        }
    }
    while let Some((i, j)) = pairs.pop_front() {
        let (f, g) = (&basis[i], &basis[j]);
        let a = order.leading_index(f);
        let b = order.leading_index(g);
        if f.exp(a)
            .iter()
            .zip(g.exp(b))
            .all(|(&x, &y)| x == 0 || y == 0)
        {
            continue;
        }
        let lcm: Vec<_> = f
            .exp(a)
            .iter()
            .zip(g.exp(b))
            .map(|(&x, &y)| x.max(y))
            .collect();
        let s = f
            .mul_term(&difference(&lcm, f.exp(a)), &f.coeff(a).inv())
            .sub(&g.mul_term(&difference(&lcm, g.exp(b)), &g.coeff(b).inv()));
        let r = reduce(&s, &basis, order).1;
        if r.is_zero() {
            continue;
        }
        if r.is_constant() {
            return vec![MPoly::one(r.nvars)];
        }
        let k = basis.len();
        for i in 0..k {
            pairs.push_back((i, k));
        }
        basis.push(monic(&r, order));
    }
    let minimal: Vec<_> = basis
        .iter()
        .enumerate()
        .filter(|(i, f)| {
            let a = f.exp(order.leading_index(f));
            !basis.iter().enumerate().any(|(j, g)| {
                let b = g.exp(order.leading_index(g));
                j != *i && divides(b, a) && (a != b || j < *i)
            })
        })
        .map(|(_, f)| f.clone())
        .collect();
    let mut out: Vec<_> = minimal
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let others: Vec<_> = minimal
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, g)| g.clone())
                .collect();
            monic(&reduce(f, &others, order).1, order)
        })
        .collect();
    out.sort_by(|a, b| order.compare(b.exp(order.leading_index(b)), a.exp(order.leading_index(a))));
    out
}

pub fn zero_dimensional<C: Field>(basis: &[MPoly<C>], order: MonomialOrder, n: usize) -> bool {
    !basis.iter().any(MPoly::is_constant)
        && (0..n).all(|v| {
            basis.iter().any(|p| {
                let e = p.exp(order.leading_index(p));
                e[v] > 0 && e.iter().enumerate().all(|(j, &k)| j == v || k == 0)
            })
        })
}

/// Change monomial order through exact linear dependence in the quotient
/// algebra. The coefficient type supplies the field operations, including
/// rational-function and algebraic number fields.
pub fn fglm<C: Field>(
    basis: &[MPoly<C>],
    source: MonomialOrder,
    target: MonomialOrder,
    nvars: usize,
) -> PResult<Vec<MPoly<C>>> {
    if source == target { return Ok(basis.to_vec()); }
    if !zero_dimensional(basis, source, nvars) {
        return Err(PolyError::Unsupported("cannot convert Groebner bases of ideals with positive dimension".into()));
    }
    let leading: Vec<_> = basis.iter().map(|p| p.exp(source.leading_index(p)).to_vec()).collect();
    let zero = vec![0; nvars];
    let mut staircase = Vec::new();
    let mut queue = VecDeque::from([zero.clone()]);
    let mut seen = FxHashSet::default();
    seen.insert(zero.clone());
    while let Some(exponent) = queue.pop_front() {
        if leading.iter().any(|m| divides(m, &exponent)) { continue; }
        staircase.push(exponent.clone());
        for variable in 0..nvars {
            let mut child = exponent.clone();
            child[variable] += 1;
            if seen.insert(child.clone()) { queue.push_back(child); }
        }
    }
    staircase.sort_by(|a,b| source.compare(a,b));
    let dimension = staircase.len();
    let indices: FxHashMap<_,_> = staircase.iter().cloned().enumerate().map(|(i,e)|(e,i)).collect();
    // matrices[variable][column][row] are exact multiplication matrices.
    let mut matrices = Vec::with_capacity(nvars);
    for variable in 0..nvars {
        let mut columns = Vec::with_capacity(dimension);
        for exponent in &staircase {
            let mut child = exponent.clone();
            child[variable] += 1;
            let remainder = reduce(&MPoly::monomial(&child, C::one()), basis, source).1;
            let mut column = vec![C::zero(); dimension];
            for (m,c) in remainder.iter() { column[indices[m]] = c.clone(); }
            columns.push(column);
        }
        matrices.push(columns);
    }
    let multiply = |variable:usize, vector:&[C]| {
        let mut out = vec![C::zero(); dimension];
        for (column,c) in vector.iter().enumerate() {
            if c.is_zero() { continue; }
            for (row,value) in out.iter_mut().enumerate() {
                *value = value.add(&matrices[variable][column][row].mul(c));
            }
        }
        out
    };
    let mut one = vec![C::zero(); dimension];
    one[indices[&zero]] = C::one();
    let mut candidates = vec![(zero.clone(), one)];
    seen.clear();
    seen.insert(zero);
    let mut accepted: Vec<Vec<u32>> = Vec::new();
    let mut rows: Vec<(usize, Vec<C>, Vec<C>)> = Vec::new();
    let mut result: Vec<MPoly<C>> = Vec::new();
    let mut target_leading: Vec<Vec<u32>> = Vec::new();
    while !candidates.is_empty() {
        let next = (0..candidates.len()).min_by(|&i,&j| target.compare(&candidates[i].0,&candidates[j].0)).unwrap();
        let (exponent,vector) = candidates.swap_remove(next);
        if target_leading.iter().any(|m| divides(m,&exponent)) { continue; }
        let mut residual = vector.clone();
        let mut coefficients = vec![C::zero(); accepted.len()];
        for (pivot,row,witness) in &rows {
            let c = residual[*pivot].clone();
            if c.is_zero() { continue; }
            for (a,b) in residual.iter_mut().zip(row) { *a = a.sub(&c.mul(b)); }
            for (a,b) in coefficients.iter_mut().zip(witness) { *a = a.add(&c.mul(b)); }
        }
        if let Some(pivot) = residual.iter().position(|c| !c.is_zero()) {
            let inverse = residual[pivot].inv();
            let row: Vec<_> = residual.iter().map(|c| c.mul(&inverse)).collect();
            let mut witness: Vec<_> = coefficients.iter().map(|c| c.neg().mul(&inverse)).collect();
            witness.push(inverse);
            rows.push((pivot,row,witness));
            accepted.push(exponent.clone());
            for variable in 0..nvars {
                let mut child = exponent.clone();
                child[variable] += 1;
                if seen.insert(child.clone()) { candidates.push((child,multiply(variable,&vector))); }
            }
        } else {
            let mut terms = vec![(exponent.clone(),C::one())];
            terms.extend(accepted.iter().zip(coefficients).filter(|(_,c)| !c.is_zero()).map(|(m,c)|(m.clone(),c.neg())));
            result.push(MPoly::from_terms(nvars,terms));
            target_leading.push(exponent);
        }
    }
    result.sort_by(|a,b| target.compare(b.exp(target.leading_index(b)),a.exp(target.leading_index(a))));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::convert::from_exprs;
    use crate::poly::int::Rat;

    #[test]
    fn fglm_generic_fields_preserve_the_ideal_in_every_order() {
        let expressions:Vec<_> = ["x*y-1", "y**2-x"].iter().map(|s|crate::parser::parse_str(s).unwrap()).collect();
        let input = from_exprs(&expressions).unwrap().polys;
        for source in [MonomialOrder::Lex,MonomialOrder::Grlex,MonomialOrder::Grevlex] {
            let original = buchberger::<Rat>(&input,source);
            for target in [MonomialOrder::Lex,MonomialOrder::Grlex,MonomialOrder::Grevlex] {
                let converted = fglm(&original,source,target,2).unwrap();
                assert_eq!(converted,buchberger(&input,target));
                assert!(original.iter().all(|p|reduce(p,&converted,target).1.is_zero()));
                assert!(converted.iter().all(|p|reduce(p,&original,source).1.is_zero()));
            }
        }
    }
}
