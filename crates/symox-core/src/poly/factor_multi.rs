//! Multivariate factorisation over `ZZ`.
//!
//! Square-free parts are factored by evaluating all variables but a main
//! one at an integer point, factoring the univariate image and lifting the
//! image factors back with a multivariate Hensel construction (in the style
//! of Wang's EEZ / Geddes–Czapor–Labahn, the leading coefficient problem is
//! solved by imposing the full leading coefficient on every factor). The
//! lifting is done over `QQ` with exact arithmetic, which removes the need
//! for coefficient bounds; when the image has spurious factors the lift
//! fails and two-way splits of the image factors are tried.

use super::dense;
use super::factor_uni::{factor_zx, zassenhaus};
use super::gcd::{content_in, gcd_zx, primitive_in};
use super::gf::Rng;
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::sqf::sqf_list_mz;

/// Full factorisation: `(content, [(factor, multiplicity)])`, factors
/// primitive with positive leading coefficients (unsorted).
pub fn factor_list_mz(f: &MPoly<Int>) -> (Int, Vec<(MPoly<Int>, u32)>) {
    let n = f.nvars;
    if f.is_zero() {
        return (Int::zero(), vec![]);
    }
    let (c, mut p) = f.primitive();
    if p.is_constant() {
        return (c.mul(&p.lc()), vec![]);
    }
    let mut out = Vec::new();
    let m = p.min_exps();
    if m.iter().any(|&e| e > 0) {
        for (v, &e) in m.iter().enumerate() {
            if e > 0 {
                out.push((MPoly::var(v, n), e));
            }
        }
        p = p.div_monomial(&m);
    }
    if p.is_constant() {
        return (c, out);
    }
    let used = p.used_vars();
    if used.len() == 1 {
        let v = used[0];
        let (_, fs) = factor_zx(&dense::from_mpoly(&p, v));
        for (g, k) in fs {
            out.push((dense::to_mpoly(&g, v, n), k));
        }
        return (c, out);
    }
    let (_, parts) = sqf_list_mz(&p);
    for (a, k) in parts {
        for g in factor_sqf_mz(&a) {
            out.push((g, k));
        }
    }
    (c, out)
}

fn normalize(g: MPoly<Int>) -> MPoly<Int> {
    let (_, p) = g.primitive();
    p
}

/// Irreducible factors of a primitive square-free polynomial.
pub fn factor_sqf_mz(f: &MPoly<Int>) -> Vec<MPoly<Int>> {
    let n = f.nvars;
    let f = normalize(f.clone());
    if f.is_constant() {
        return vec![];
    }
    let used = f.used_vars();
    if used.len() == 1 {
        let v = used[0];
        return zassenhaus(&dense::from_mpoly(&f, v))
            .into_iter()
            .map(|g| dense::to_mpoly(&g, v, n))
            .collect();
    }
    let degs = f.degrees();
    let v = *used.iter().min_by_key(|&&v| (degs[v], v)).unwrap();
    let cont = content_in(&f, v);
    if !cont.is_constant() {
        let rest = f.div_exact(&cont).expect("content divides");
        let mut out = factor_sqf_mz(&cont);
        out.extend(factor_sqf_mz(&rest));
        return out;
    }
    if degs[v] == 1 {
        return vec![f];
    }
    eez(&f, v).into_iter().map(normalize).collect()
}

fn eval_point(f: &MPoly<Int>, ys: &[usize], a: &[Int]) -> MPoly<Int> {
    let mut g = f.clone();
    for (y, x) in ys.iter().zip(a) {
        g = g.eval(*y, x);
    }
    g
}

fn eez(f: &MPoly<Int>, v: usize) -> Vec<MPoly<Int>> {
    let ys: Vec<usize> = f.used_vars().into_iter().filter(|&w| w != v).collect();
    let lc = f.lc_in(v);
    let d = f.degree(v);
    let mut rng = Rng::new(0x3ee2 + f.len() as u64);
    let mut best: Option<(Vec<Int>, Vec<Vec<Int>>)> = None;
    let mut good = 0;
    // SymPy's dmp_zz_wang widens the evaluation range until it finds
    // valid points; a fixed number of attempts would make an input whose
    // leading coefficient vanishes on the first ranges look irreducible.
    // The range keeps growing, so a nonzero leading coefficient and
    // discriminant eventually have a non-root (inside a budget scope the
    // attempt limit raises instead).
    for attempt in 0u64.. {
        if good >= 3 || (good >= 1 && attempt >= 60) {
            break;
        }
        super::factor_budget::event("multivariate evaluation attempts");
        let a: Vec<Int> = if attempt == 0 {
            vec![Int::zero(); ys.len()]
        } else {
            let b = 1 + attempt / 4;
            ys.iter()
                .map(|_| Int::S((rng.next() % (2 * b + 1)) as i64 - b as i64))
                .collect()
        };
        if eval_point(&lc, &ys, &a).is_zero() {
            continue;
        }
        let u = dense::from_mpoly(&eval_point(f, &ys, &a), v);
        if dense::deg(&u) != d as isize {
            continue;
        }
        let du = dense::derivative(&u);
        if gcd_zx(&u, &du).0.len() > 1 {
            continue;
        }
        let (_, fs) = factor_zx(&u);
        let us: Vec<Vec<Int>> = fs.into_iter().map(|(g, _)| g).collect();
        if us.len() == 1 {
            return vec![f.clone()];
        }
        good += 1;
        if best.as_ref().is_none_or(|b| us.len() < b.1.len()) {
            best = Some((a, us));
        }
    }
    let Some((a, us)) = best else {
        super::factor_budget::require(false,"uncertified multivariate evaluation");
        return vec![f.clone()];
    };
    if let Some(gs) = lift(f, v, &ys, &a, &us) {
        return gs;
    }
    // spurious image factors: look for a two way split
    let r = us.len();
    for s in 1..=r / 2 {
        let mut idx: Vec<usize> = (0..s).collect();
        loop {
            super::factor_budget::event("multivariate recombination");
            if !(2 * s == r && idx[0] != 0) {
                let mut u1 = vec![Int::one()];
                let mut u2 = vec![Int::one()];
                for (i, u) in us.iter().enumerate() {
                    if idx.contains(&i) {
                        u1 = dense::mul_int(&u1, u);
                    } else {
                        u2 = dense::mul_int(&u2, u);
                    }
                }
                if let Some(gs) = lift(f, v, &ys, &a, &[u1, u2]) {
                    let mut out = Vec::new();
                    for g in gs {
                        out.extend(factor_sqf_mz(&g));
                    }
                    return out;
                }
            }
            if !next_comb(&mut idx, r) {
                break;
            }
        }
    }
    vec![f.clone()]
}

fn next_comb(idx: &mut [usize], n: usize) -> bool {
    let s = idx.len();
    for k in (0..s).rev() {
        if idx[k] != k + n - s {
            idx[k] += 1;
            for j in k + 1..s {
                idx[j] = idx[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

/// Solutions of univariate multi-term diophantine equations
/// `Σ σ_i Π_{j≠i} a_j = c` with `deg σ_i < deg a_i`.
struct UniDio {
    a: Vec<Vec<Rat>>,
    s: Vec<Vec<Rat>>,
}

impl UniDio {
    fn new(a: Vec<Vec<Rat>>) -> Option<UniDio> {
        let r = a.len();
        let mut s = Vec::with_capacity(r);
        for i in 0..r {
            let mut b = vec![Rat::one()];
            for (j, aj) in a.iter().enumerate() {
                if j != i {
                    b = dense::rem(&dense::mul(&b, aj), &a[i]);
                }
            }
            let (g, si, _) = dense::gcdex_field(&b, &a[i]);
            if g.len() != 1 {
                return None;
            }
            s.push(dense::rem(&si, &a[i]));
        }
        Some(UniDio { a, s })
    }

    fn solve(&self, c: &[Rat]) -> Vec<Vec<Rat>> {
        self.a
            .iter()
            .zip(&self.s)
            .map(|(ai, si)| dense::rem(&dense::mul(c, si), ai))
            .collect()
    }
}

fn product(us: &[MPoly<Rat>], n: usize) -> MPoly<Rat> {
    let mut p = MPoly::one(n);
    for u in us {
        p = p.mul(u);
    }
    p
}

fn cofactors(a: &[MPoly<Rat>], n: usize) -> Vec<MPoly<Rat>> {
    (0..a.len())
        .map(|i| {
            let mut p = MPoly::one(n);
            for (j, aj) in a.iter().enumerate() {
                if j != i {
                    p = p.mul(aj);
                }
            }
            p
        })
        .collect()
}

fn ypow(y: usize, k: u32, n: usize) -> Vec<u32> {
    let mut e = vec![0u32; n];
    e[y] = k;
    e
}

/// Multivariate multi-term diophantine solver (variables `ys` evaluated at
/// zero, solutions truncated at degree `d` in each of them).
fn mtd(
    a: &[MPoly<Rat>],
    c: &MPoly<Rat>,
    ys: &[usize],
    d: u32,
    v: usize,
    dio: &UniDio,
) -> Vec<MPoly<Rat>> {
    let n = c.nvars;
    if ys.is_empty() {
        let cd = dense::from_mpoly(c, v);
        return dio.solve(&cd).iter().map(|s| dense::to_mpoly(s, v, n)).collect();
    }
    let ym = *ys.last().unwrap();
    let rest = &ys[..ys.len() - 1];
    let anew: Vec<MPoly<Rat>> = a.iter().map(|x| x.eval(ym, &Rat::zero())).collect();
    let cnew = c.eval(ym, &Rat::zero());
    let mut sigma = mtd(&anew, &cnew, rest, d, v, dio);
    let b = cofactors(a, n);
    let mut e = c.clone();
    for (s, bi) in sigma.iter().zip(&b) {
        e = e.sub(&s.mul(bi));
    }
    for m in 1..=d {
        if e.is_zero() {
            break;
        }
        let cm = e.coeff_of(ym, m);
        if cm.is_zero() {
            continue;
        }
        let ds = mtd(&anew, &cm, rest, d, v, dio);
        let mono = ypow(ym, m, n);
        for (i, dsi) in ds.into_iter().enumerate() {
            let t = dsi.mul_term(&mono, &Rat::one());
            e = e.sub(&t.mul(&b[i]));
            sigma[i] = sigma[i].add(&t);
        }
    }
    sigma
}

/// Lift the univariate factors `us` of `f(x_v, a)` to a factorisation of
/// `f`; `None` if they do not correspond to a true factorisation.
fn lift(f: &MPoly<Int>, v: usize, ys: &[usize], a: &[Int], us: &[Vec<Int>]) -> Option<Vec<MPoly<Int>>> {
    super::factor_budget::event("multivariate Hensel lifting");
    let n = f.nvars;
    let r = us.len();
    let mut ff = f.clone();
    for (y, x) in ys.iter().zip(a) {
        ff = ff.shift(*y, x);
    }
    let l = ff.lc_in(v);
    let l0 = l.constant_term();
    if l0.is_zero() {
        return None;
    }
    let fp = if r > 1 { ff.mul(&l.pow(r as u32 - 1)) } else { ff.clone() };
    let fq = fp.to_rat();
    let lq = l.to_rat();
    let mut u: Vec<MPoly<Rat>> = us
        .iter()
        .map(|g| {
            let s = Rat::new(l0.clone(), g.last().unwrap().clone());
            let gq: Vec<Rat> = g.iter().map(|c| s.mul(&Rat::from(c.clone()))).collect();
            dense::to_mpoly(&gq, v, n)
        })
        .collect();
    let dio = UniDio::new(u.iter().map(|x| dense::from_mpoly(x, v)).collect())?;
    let maxdeg = ys.iter().map(|&y| fp.degree(y)).max().unwrap_or(0).max(0) as u32;
    for j in 0..ys.len() {
        let yj = ys[j];
        let mut fj = fq.clone();
        let mut lj = lq.clone();
        for &y in &ys[j + 1..] {
            fj = fj.eval(y, &Rat::zero());
            lj = lj.eval(y, &Rat::zero());
        }
        let u1 = u.clone();
        for ui in u.iter_mut() {
            let du = ui.degree(v) as u32;
            let e = ypow(v, du, n);
            let old = ui.lc_in(v);
            *ui = ui.sub(&old.mul_term(&e, &Rat::one())).add(&lj.mul_term(&e, &Rat::one()));
        }
        let mut e = fj.sub(&product(&u, n));
        let dj = fj.degree(yj).max(0) as u32;
        for k in 1..=dj {
            super::factor_budget::event("multivariate Hensel coefficients");
            if e.is_zero() {
                break;
            }
            let c = e.coeff_of(yj, k);
            if c.is_zero() {
                continue;
            }
            let ds = mtd(&u1, &c, &ys[..j], maxdeg, v, &dio);
            let mono = ypow(yj, k, n);
            for (i, dsi) in ds.into_iter().enumerate() {
                u[i] = u[i].add(&dsi.mul_term(&mono, &Rat::one()));
            }
            e = fj.sub(&product(&u, n));
        }
        if !e.is_zero() {
            return None;
        }
    }
    if product(&u, n) != fq {
        return None;
    }
    let mut gs = Vec::with_capacity(r);
    for ui in u {
        let (_, zi) = ui.clear_denoms();
        let (_, mut g) = primitive_in(&zi, v);
        for (y, x) in ys.iter().zip(a) {
            g = g.shift(*y, &x.neg());
        }
        gs.push(normalize(g));
    }
    // final check against the original polynomial
    let mut p = MPoly::one(n);
    for g in &gs {
        p = p.mul(g);
    }
    if p != *f && p.neg() != *f {
        return None;
    }
    Some(gs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(f: &MPoly<Int>, expected: usize) {
        let (c, fs) = factor_list_mz(f);
        let mut p = MPoly::constant(c, f.nvars);
        for (g, k) in &fs {
            p = p.mul(&g.pow(*k));
        }
        assert_eq!(&p, f);
        let total: usize = fs.len();
        assert_eq!(total, expected, "{fs:?}");
    }

    #[test]
    fn bivariate() {
        let x = MPoly::<Int>::var(0, 2);
        let y = MPoly::<Int>::var(1, 2);
        let one = MPoly::<Int>::one(2);
        // x^2 - y^2
        check(&x.pow(2).sub(&y.pow(2)), 2);
        // x^3 y - x y^3
        check(&x.pow(3).mul(&y).sub(&x.mul(&y.pow(3))), 4);
        // (x+y)^2 (x-y)
        check(&x.add(&y).pow(2).mul(&x.sub(&y)), 2);
        // (x y + 1)(x^2 + y^2 + 1)(x - y + 3)
        let f = x.mul(&y).add(&one).mul(&x.pow(2).add(&y.pow(2)).add(&one)).mul(&x.sub(&y).add(&one.scale(&Int::S(3))));
        check(&f, 3);
        // non-monic leading coefficients: (y x^2 + 1)(y^2 x + x + 2)
        let f = y.mul(&x.pow(2)).add(&one).mul(&y.pow(2).mul(&x).add(&x).add(&one.scale(&Int::S(2))));
        check(&f, 2);
        // irreducible
        check(&x.pow(2).sub(&y.pow(3)).sub(&one.scale(&Int::S(2))), 1);
    }

    #[test]
    fn trivariate() {
        let x = MPoly::<Int>::var(0, 3);
        let y = MPoly::<Int>::var(1, 3);
        let z = MPoly::<Int>::var(2, 3);
        let s = x.add(&y).add(&z);
        check(&s.pow(4), 1);
        let f = x.mul(&y).add(&z).mul(&x.add(&y.mul(&z))).mul(&x.sub(&z).add(&MPoly::one(3)));
        check(&f, 3);
        // x^2 - y^2 z^2 = (x - yz)(x + yz)
        check(&x.pow(2).sub(&y.pow(2).mul(&z.pow(2))), 2);
    }
}
