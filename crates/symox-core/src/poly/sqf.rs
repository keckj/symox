//! Square-free decomposition (Yun's algorithm) over `ZZ`.

use super::dense;
use super::gcd::{gcd_mz, gcd_zx, primitive_in};
use super::int::Int;
use super::mpoly::MPoly;

/// Yun's algorithm for a primitive univariate polynomial of positive
/// degree: `[(a_i, i)]` with `f = Π a_i^i`, `a_i` square-free, pairwise
/// coprime, of positive degree (positive leading coefficients).
pub fn yun_zx(f: &[Int]) -> Vec<(Vec<Int>, u32)> {
    let mut out = Vec::new();
    if f.len() <= 1 {
        return out;
    }
    let df = dense::derivative(f);
    let (_, mut c, d0) = gcd_zx(f, &df);
    // d = f'/b - c'
    let mut d = dense::sub(&d0, &dense::derivative(&c));
    let mut i = 1u32;
    while c.len() > 1 {
        let (a, c2, d2) = gcd_zx(&c, &d);
        if a.len() > 1 {
            out.push((a, i));
        }
        c = c2;
        d = dense::sub(&d2, &dense::derivative(&c));
        i += 1;
    }
    out
}

/// Square-free decomposition of a univariate integer polynomial:
/// `(content, [(a_i, i)])` (content with the sign of the leading
/// coefficient).
pub fn sqf_list_zx(f: &[Int]) -> (Int, Vec<(Vec<Int>, u32)>) {
    if f.is_empty() {
        return (Int::zero(), vec![]);
    }
    let (c, p) = dense::primitive(f);
    if p.len() == 1 {
        return (c.mul(&p[0]), vec![]);
    }
    (c, yun_zx(&p))
}

/// Yun's algorithm with respect to `x_v` for a polynomial primitive in
/// `x_v` (all square-free parts have positive degree in `x_v`).
fn yun_mz(f: &MPoly<Int>, v: usize) -> Vec<(MPoly<Int>, u32)> {
    let mut out = Vec::new();
    let df = f.derivative(v);
    let (_, mut c, d0) = gcd_mz(f, &df);
    let mut d = d0.sub(&c.derivative(v));
    let mut i = 1u32;
    while c.degree(v) > 0 {
        let (a, c2, d2) = gcd_mz(&c, &d);
        if a.degree(v) > 0 {
            out.push((a, i));
        }
        c = c2;
        d = d2.sub(&c.derivative(v));
        i += 1;
    }
    out
}

/// Square-free decomposition of a multivariate integer polynomial,
/// following the variable order (main variable first, then the content
/// recursively): parts of equal multiplicity are multiplied together.
pub fn sqf_list_mz(f: &MPoly<Int>) -> (Int, Vec<(MPoly<Int>, u32)>) {
    if f.is_zero() {
        return (Int::zero(), vec![]);
    }
    let (c, p) = f.primitive();
    if p.is_constant() {
        return (c.mul(&p.lc()), vec![]);
    }
    let mut parts: Vec<(MPoly<Int>, u32)> = Vec::new();
    sqf_rec(&p, 0, &mut parts);
    parts.sort_by_key(|x| x.1);
    (c, parts)
}

fn sqf_rec(f: &MPoly<Int>, v0: usize, parts: &mut Vec<(MPoly<Int>, u32)>) {
    let n = f.nvars;
    let Some(v) = (v0..n).find(|&v| f.degree(v) > 0) else {
        return;
    };
    let (cont, pp) = primitive_in(f, v);
    let add = |a: MPoly<Int>, k: u32, parts: &mut Vec<(MPoly<Int>, u32)>| {
        if let Some(e) = parts.iter_mut().find(|x| x.1 == k) {
            e.0 = e.0.mul(&a);
        } else {
            parts.push((a, k));
        }
    };
    for (a, k) in yun_mz(&pp, v) {
        let a = if a.lc().is_negative() { a.neg() } else { a };
        add(a, k, parts);
    }
    if !cont.is_constant() {
        let mut sub = Vec::new();
        sqf_rec(&cont, v + 1, &mut sub);
        for (a, k) in sub {
            add(a, k, parts);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zi(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&x| Int::S(x)).collect()
    }

    #[test]
    fn yun() {
        // (x-1)^2 (x+1) = x^3 - x^2 - x + 1
        let f = zi(&[1, -1, -1, 1]);
        let r = sqf_list_zx(&f);
        assert_eq!(r.1, vec![(zi(&[1, 1]), 1), (zi(&[-1, 1]), 2)]);
    }

    #[test]
    fn yun_multi() {
        let x = MPoly::<Int>::var(0, 2);
        let y = MPoly::<Int>::var(1, 2);
        let a = x.add(&y);
        let b = x.sub(&y);
        let f = a.pow(2).mul(&b);
        let (_, parts) = sqf_list_mz(&f);
        assert_eq!(parts, vec![(b, 1), (a, 2)]);
        let g = y.add(&MPoly::one(2)).pow(3).mul(&x);
        let (_, parts) = sqf_list_mz(&g);
        assert_eq!(parts.len(), 2);
    }
}
