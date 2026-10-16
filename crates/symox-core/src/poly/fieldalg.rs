//! Ordered sparse algorithms over an arbitrary exact coefficient field.

use super::groebner::MonomialOrder;
use super::int::Field;
use super::mpoly::MPoly;
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
