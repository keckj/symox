//! Hankel transform `F_nu(k) = int_0^oo f(r) J_nu(k*r) r dr` and its inverse
//! (the transform is its own inverse with the variables exchanged).
//!
//! Table (Bateman, Tables of Integral Transforms, vol. II, ch. 8):
//!
//! * `r**m`: `2**(m+1) Gamma((nu+m+2)/2) / (Gamma((nu-m)/2) k**(m+2))`,
//! * `exp(-a*r)`, `exp(-a*r)/r`,
//! * `r**nu * exp(-a*r**2)`,
//! * `r**nu / (r**2 + a**2)**(mu+1)` (Bessel `K` closed forms),
//! * `r**nu * Heaviside(a - r)` (disc).

use super::util::*;
use crate::expr::{Expr, Func, Kind};
use crate::functions::{exp, gamma};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};

type V = (Expr, Expr);

/// Direct evaluation of a Hankel transform (`inverse` only changes the
/// roles of the variables): `Some((F, cond))` or `None`.
pub fn hankel_type_transform(inverse: bool, f: &Expr, r: &Expr, k: &Expr, nu: &Expr, simplify: bool) -> Option<V> {
    let _ = inverse;
    // the transform variable is positive (the kernel is J_nu(k*r), k > 0)
    let kp = match k.kind() {
        Kind::Sym(_) => Expr::dummy("k", crate::symbol::Assumptions::POSITIVE),
        _ => k.clone(),
    };
    let (v, c) = ht(f, r, &kp, nu, 0)?;
    if v.has_nonfinite() {
        return None;
    }
    let v = if simplify { super::fourier::finish(&v) } else { v };
    let back = |e: &Expr| crate::subs::subs1(e, &kp, k);
    Some((back(&v), back(&c)))
}

fn ht(f: &Expr, r: &Expr, k: &Expr, nu: &Expr, depth: u32) -> Option<V> {
    if depth > MAX_DEPTH {
        return None;
    }
    if f.is_add() {
        let mut vals = Vec::new();
        let mut conds = Vec::new();
        for t in terms(f) {
            let (v, c) = ht(&t, r, k, nu, depth + 1)?;
            vals.push(v);
            conds.push(c);
        }
        return Some((add(&vals), crate::logic::and(conds)));
    }
    let (c, g) = split_mul(f, r);
    if g.is_one() {
        // constant: the pure power r**0
        return power(&Expr::zero(), &c, k, nu);
    }
    let (m, g) = split_power(&g, r);
    if g.is_one() {
        return power(&m, &c, k, nu);
    }
    let two = int(2);
    match g.kind() {
        Kind::Func(Func::Exp, a) => {
            let e = neg(&a[0]);
            if let Some(al) = scaled(&e, r) {
                if !re_pos(&al) {
                    return None;
                }
                let q = sqrt(&add(&[pow(&al, &two), pow(k, &two)]));
                if m.is_zero() {
                    // k**nu (a + nu*q) / (q**3 (a + q)**nu)
                    let v = div(
                        &mul(&[pow(k, nu), add(&[al.clone(), mul(&[nu.clone(), q.clone()])])]),
                        &mul(&[pow(&q, &int(3)), pow(&add(&[al.clone(), q.clone()]), nu)]),
                    );
                    return Some((mul(&[c, v]), tru()));
                }
                if m == int(-1) {
                    let v = div(&pow(&div(&sub(&q, &al), k), nu), &q);
                    return Some((mul(&[c, v]), tru()));
                }
                return None;
            }
            // r**nu exp(-a r**2)
            let p = poly(&crate::expand::expand(&e), r)?;
            if p.len() != 3 || !p[0].is_zero() || !p[1].is_zero() {
                return None;
            }
            let al = p[2].clone();
            if !re_pos(&al) || crate::expand::expand(&sub(&m, nu)) != Expr::zero() {
                return None;
            }
            let v = div(
                &mul(&[pow(k, nu), exp(&neg(&div(&pow(k, &two), &mul(&[int(4), al.clone()]))))]),
                &pow(&mul(&[two, al]), &add(&[nu.clone(), Expr::one()])),
            );
            Some((mul(&[c, v]), tru()))
        }
        Kind::Pow(b, e) if !has(e, r) => {
            // r**nu / (r**2 + a**2)**(mu + 1)
            let p = poly(&crate::expand::expand(b), r)?;
            if p.len() != 3 || !p[1].is_zero() || !p[2].is_one() || !pos(&p[0]) {
                return None;
            }
            if crate::expand::expand(&sub(&m, nu)) != Expr::zero() {
                return None;
            }
            let a = sqrt(&p[0]);
            let mu = sub(&neg(e), &Expr::one());
            let v = div(
                &mul(&[pow(&a, &sub(nu, &mu)), pow(k, &mu), besselk_half(&sub(nu, &mu), &mul(&[a, k.clone()]))]),
                &mul(&[pow(&two, &mu), gamma(&add(&[mu.clone(), Expr::one()]))]),
            );
            Some((mul(&[c, v]), tru()))
        }
        Kind::Func(Func::Heaviside, a) => {
            // disc: r**nu * Heaviside(a - r)
            let (sl, a0) = linear(&a[0], r)?;
            if !negv(&sl) || crate::expand::expand(&sub(&m, nu)) != Expr::zero() {
                return None;
            }
            let rad = neg(&div(&a0, &sl));
            if !pos(&rad) {
                return None;
            }
            let v = div(
                &mul(&[pow(&rad, &add(&[nu.clone(), Expr::one()])), besselj(&add(&[nu.clone(), Expr::one()]), &mul(&[k.clone(), rad]))]),
                k,
            );
            Some((mul(&[c, v]), tru()))
        }
        _ => None,
    }
}

/// `int_0^oo r**(m+1) J_nu(k r) dr`.
fn power(m: &Expr, c: &Expr, k: &Expr, nu: &Expr) -> Option<V> {
    let mu = add(&[m.clone(), Expr::one()]);
    let lo = crate::relational::gt(&crate::functions::re(&mu), &sub(&neg(&crate::functions::re(nu)), &Expr::one()));
    let hi = crate::relational::lt(&crate::functions::re(&mu), &rat(1, 2));
    let mut conds = Vec::new();
    for cd in [lo, hi] {
        match truth(&cd) {
            Some(false) => return None,
            Some(true) => {}
            None => conds.push(cd),
        }
    }
    let two = int(2);
    let v = div(
        &mul(&[pow(&two, &mu), gamma(&div(&add(&[nu.clone(), mu.clone(), Expr::one()]), &two))]),
        &mul(&[gamma(&div(&add(&[sub(nu, &mu), Expr::one()]), &two)), pow(k, &add(&[mu, Expr::one()]))]),
    );
    Some((mul(&[c.clone(), v]), crate::logic::and(conds)))
}

/// `besselk(nu, z)`, expanded in closed form for half integer orders.
fn besselk_half(nu: &Expr, z: &Expr) -> Expr {
    let n2 = mul(&[int(2), nu.clone()]);
    if let Some(m) = n2.as_i64() {
        if m % 2 != 0 {
            let n = (m.abs() - 1) / 2;
            return k_half(n, z);
        }
    }
    besselk(nu, z)
}

/// `K_{n+1/2}(z) = sqrt(pi/(2z)) exp(-z) sum_j (n+j)!/(j!(n-j)!(2z)**j)`.
fn k_half(n: i64, z: &Expr) -> Expr {
    let mut sum = Vec::new();
    let fact = |v: i64| -> Expr { Expr::int((1..=v).product::<i64>().max(1)) };
    for j in 0..=n {
        let c = div(&fact(n + j), &mul(&[fact(j), fact(n - j)]));
        sum.push(div(&c, &pow(&mul(&[int(2), z.clone()]), &int(j))));
    }
    mul(&[sqrt(&div(&Expr::pi(), &mul(&[int(2), z.clone()]))), exp(&neg(z)), add(&sum)])
}
