//! Greatest common divisors over `ZZ` and `QQ`.
//!
//! Univariate: heuristic gcd (evaluation at a large integer, gcd of the
//! integers, reconstruction of the polynomial from its ξ-adic expansion)
//! with a modular (CRT) fallback. Multivariate: the same heuristic applied
//! recursively, one variable at a time, with a recursive primitive PRS
//! fallback.

use super::dense;
use super::gf;
use super::int::{Int, Rat};
use super::mpoly::MPoly;

const HEU_GCD_MAX: usize = 6;

/// Next evaluation point of the heuristic gcd.
fn next_xi(x: &Int) -> Int {
    x.mul_i64(73794).mul(&x.isqrt().isqrt()).div_mod_floor(&Int::S(27011)).0
}

fn initial_xi(f_norm: &Int, g_norm: &Int, f_lc: &Int, g_lc: &Int) -> Int {
    let mn = if f_norm < g_norm { f_norm } else { g_norm };
    let b = mn.mul_i64(2).add(&Int::S(29));
    let s = b.isqrt().mul_i64(99);
    let a = if b < s { b } else { s };
    let qf = f_norm.div_mod_floor(&f_lc.abs()).0;
    let qg = g_norm.div_mod_floor(&g_lc.abs()).0;
    let q = if qf < qg { qf } else { qg };
    let c = q.mul_i64(2).add(&Int::S(2));
    if a > c {
        a
    } else {
        c
    }
}

/// ξ-adic reconstruction of an integer polynomial from `h = H(ξ)`.
fn interpolate_int(h: &Int, x: &Int) -> Vec<Int> {
    let mut out = Vec::new();
    let mut h = h.clone();
    while !h.is_zero() {
        let g = h.symmetric_mod(x);
        h = h.sub(&g).div_exact_unchecked(x);
        out.push(g);
    }
    out
}

fn neg_if_needed(p: Vec<Int>) -> Vec<Int> {
    if p.last().is_some_and(|c| c.is_negative()) {
        dense::neg(&p)
    } else {
        p
    }
}

/// Univariate gcd over `ZZ` with cofactors: `(h, f/h, g/h)`, `h` with a
/// positive leading coefficient (or zero).
pub fn gcd_zx(f: &[Int], g: &[Int]) -> (Vec<Int>, Vec<Int>, Vec<Int>) {
    if f.is_empty() && g.is_empty() {
        return (vec![], vec![], vec![]);
    }
    if f.is_empty() {
        let (c, _) = dense::primitive(g);
        let s = if c.is_negative() { Int::S(-1) } else { Int::one() };
        return (dense::scale(g, &s), vec![], vec![s]);
    }
    if g.is_empty() {
        let (c, _) = dense::primitive(f);
        let s = if c.is_negative() { Int::S(-1) } else { Int::one() };
        return (dense::scale(f, &s), vec![s], vec![]);
    }
    if f.len() == 1 || g.len() == 1 {
        let c = dense::content(f).gcd(&dense::content(g));
        let cf: Vec<Int> = f.iter().map(|x| x.div_exact_unchecked(&c)).collect();
        let cg: Vec<Int> = g.iter().map(|x| x.div_exact_unchecked(&c)).collect();
        return (vec![c], cf, cg);
    }
    if f == g {
        let (c, pp) = dense::primitive(f);
        let _ = pp;
        let s = if c.is_negative() { Int::S(-1) } else { Int::one() };
        return (dense::scale(f, &s), vec![s.clone()], vec![s]);
    }
    if let Some(r) = heugcd_zx(f, g) {
        return r;
    }
    modular_gcd_zx(f, g)
}

fn heugcd_zx(f: &[Int], g: &[Int]) -> Option<(Vec<Int>, Vec<Int>, Vec<Int>)> {
    // extract the common content
    let cont = dense::content(f).gcd(&dense::content(g));
    let f: Vec<Int> = f.iter().map(|c| c.div_exact_unchecked(&cont)).collect();
    let g: Vec<Int> = g.iter().map(|c| c.div_exact_unchecked(&cont)).collect();
    let f_norm = dense::max_norm(&f);
    let g_norm = dense::max_norm(&g);
    let mut x = initial_xi(&f_norm, &g_norm, f.last().unwrap(), g.last().unwrap());
    let fin = |h: Vec<Int>, cff: Vec<Int>, cfg: Vec<Int>| {
        (dense::scale(&h, &cont), cff, cfg)
    };
    let degree = f.len().max(g.len()) as u64;
    let norm_bits = f_norm.bits().max(g_norm.bits());
    for _ in 0..HEU_GCD_MAX {
        // Horner evaluation at `x` builds integers of about `degree * bits(x)`
        // bits; when that would exceed the active factorization allowance the
        // exact modular gcd is used instead (SymPy's heuristic has no such
        // allowance, so this keeps large-coefficient inputs factorable).
        let xb = x.bits();
        let eval_bits = xb.saturating_mul(degree).saturating_add(norm_bits);
        let eval_work = degree.saturating_mul(degree).saturating_mul(xb).saturating_mul(xb) / 64;
        if !super::factor_budget::affordable(eval_bits, eval_work) {
            return None;
        }
        let ff = dense::eval(&f, &x);
        let gg = dense::eval(&g, &x);
        if !ff.is_zero() && !gg.is_zero() {
            let h = ff.gcd(&gg);
            let cff = ff.div_exact_unchecked(&h);
            let cfg = gg.div_exact_unchecked(&h);
            let hp = dense::primitive(&interpolate_int(&h, &x)).1;
            if let Some(a) = dense::div_exact_int(&f, &hp) {
                if let Some(b) = dense::div_exact_int(&g, &hp) {
                    return Some(fin(hp, a, b));
                }
            }
            let cffp = neg_if_needed(interpolate_int(&cff, &x));
            if let Some(h2) = dense::div_exact_int(&f, &cffp) {
                if let Some(b) = dense::div_exact_int(&g, &h2) {
                    return Some(fin(h2, cffp, b));
                }
            }
            let cfgp = neg_if_needed(interpolate_int(&cfg, &x));
            if let Some(h2) = dense::div_exact_int(&g, &cfgp) {
                if let Some(a) = dense::div_exact_int(&f, &h2) {
                    return Some(fin(h2, a, cfgp));
                }
            }
        }
        x = next_xi(&x);
    }
    None
}

/// Modular gcd (Brown/Collins) for univariate integer polynomials.
fn modular_gcd_zx(f: &[Int], g: &[Int]) -> (Vec<Int>, Vec<Int>, Vec<Int>) {
    let (cf, pf) = dense::primitive(f);
    let (cg, pg) = dense::primitive(g);
    let c = cf.gcd(&cg);
    let lcg = pf.last().unwrap().gcd(pg.last().unwrap());
    let mut p: u64 = (1 << 31) - 1;
    let mut modulus = Int::one();
    let mut acc: Vec<Int> = vec![];
    let mut acc_deg: isize = isize::MAX;
    loop {
        p = prev_prime(p);
        if lcg.rem_u64(p) == 0 {
            continue;
        }
        let a = gf::from_int(&pf, p);
        let b = gf::from_int(&pg, p);
        let mut h = gf::gcd(&a, &b, p);
        let d = dense::deg(&h);
        if d == 0 {
            let s = if c.is_negative() { c.neg() } else { c.clone() };
            let s = s.abs();
            return (
                vec![s.clone()],
                f.iter().map(|x| x.div_exact_unchecked(&s)).collect(),
                g.iter().map(|x| x.div_exact_unchecked(&s)).collect(),
            );
        }
        h = gf::scale(&h, lcg.rem_u64(p), p);
        if d > acc_deg {
            continue;
        }
        if d < acc_deg {
            acc = h.iter().map(|&v| Int::from_u64(v)).collect();
            modulus = Int::from_u64(p);
            acc_deg = d;
        } else {
            // CRT combine
            let pm = Int::from_u64(p);
            let inv = gf::invmod(modulus.rem_u64(p), p);
            let mp = modulus.mul(&pm);
            let mut changed = false;
            for (i, ai) in acc.iter_mut().enumerate() {
                let hi = h.get(i).copied().unwrap_or(0);
                let ar = ai.rem_u64(p);
                let delta = gf::mulmod((hi + p - ar) % p, inv, p);
                let v = ai.add(&modulus.mul(&Int::from_u64(delta)));
                if !changed && v.symmetric_mod(&mp) != ai.symmetric_mod(&modulus) {
                    changed = true;
                }
                *ai = v;
            }
            modulus = mp;
            if changed {
                continue;
            }
        }
        let cand: Vec<Int> = acc.iter().map(|v| v.symmetric_mod(&modulus)).collect();
        let hp = dense::primitive(&cand).1;
        if let (Some(a), Some(b)) = (dense::div_exact_int(&pf, &hp), dense::div_exact_int(&pg, &hp)) {
            let h = dense::scale(&hp, &c.abs());
            let sa = cf.div_exact_unchecked(&c.abs());
            let sb = cg.div_exact_unchecked(&c.abs());
            return (h, dense::scale(&a, &sa), dense::scale(&b, &sb));
        }
    }
}

fn prev_prime(mut p: u64) -> u64 {
    loop {
        p -= 1;
        if gf::is_prime(p) {
            return p;
        }
    }
}

// ----------------------------------------------------------------------
// multivariate
// ----------------------------------------------------------------------

/// Sign normalisation: positive leading coefficient.
fn normalize_sign(h: MPoly<Int>, a: MPoly<Int>, b: MPoly<Int>) -> (MPoly<Int>, MPoly<Int>, MPoly<Int>) {
    if h.lc().is_negative() {
        (h.neg(), a.neg(), b.neg())
    } else {
        (h, a, b)
    }
}

fn used_union(f: &MPoly<Int>, g: &MPoly<Int>) -> Vec<usize> {
    let df = f.degrees();
    let dg = g.degrees();
    (0..f.nvars).filter(|&v| df[v] > 0 || dg[v] > 0).collect()
}

/// Gcd over `ZZ` of two multivariate polynomials, with cofactors.
pub fn gcd_mz(f: &MPoly<Int>, g: &MPoly<Int>) -> (MPoly<Int>, MPoly<Int>, MPoly<Int>) {
    let n = f.nvars;
    if f.is_zero() && g.is_zero() {
        return (MPoly::zero(n), MPoly::zero(n), MPoly::zero(n));
    }
    if f.is_zero() {
        let s = if g.lc().is_negative() { Int::S(-1) } else { Int::one() };
        return (g.scale(&s), MPoly::zero(n), MPoly::constant(s, n));
    }
    if g.is_zero() {
        let s = if f.lc().is_negative() { Int::S(-1) } else { Int::one() };
        return (f.scale(&s), MPoly::constant(s, n), MPoly::zero(n));
    }
    if f == g {
        let s = if f.lc().is_negative() { Int::S(-1) } else { Int::one() };
        return (f.scale(&s), MPoly::constant(s.clone(), n), MPoly::constant(s, n));
    }
    if f.is_constant() || g.is_constant() {
        let c = f.content().gcd(&g.content());
        let a = f.div_coeff(&c).unwrap();
        let b = g.div_coeff(&c).unwrap();
        return (MPoly::constant(c, n), a, b);
    }
    // monomial content
    let mf = f.min_exps();
    let mg = g.min_exps();
    if mf.iter().any(|&e| e > 0) || mg.iter().any(|&e| e > 0) {
        let m: Vec<u32> = mf.iter().zip(&mg).map(|(a, b)| *a.min(b)).collect();
        let f1 = f.div_monomial(&mf);
        let g1 = g.div_monomial(&mg);
        let (h, a, b) = gcd_mz(&f1, &g1);
        let fa: Vec<u32> = mf.iter().zip(&m).map(|(x, y)| x - y).collect();
        let gb: Vec<u32> = mg.iter().zip(&m).map(|(x, y)| x - y).collect();
        let one = Int::one();
        return (h.mul_term(&m, &one), a.mul_term(&fa, &one), b.mul_term(&gb, &one));
    }
    let vars = used_union(f, g);
    if vars.len() == 1 {
        let v = vars[0];
        let (h, a, b) = gcd_zx(&dense::from_mpoly(f, v), &dense::from_mpoly(g, v));
        return (dense::to_mpoly(&h, v, n), dense::to_mpoly(&a, v, n), dense::to_mpoly(&b, v, n));
    }
    // a variable occurring in only one of the polynomials: the gcd is the
    // gcd of the other polynomial with the coefficients in that variable
    let df = f.degrees();
    let dg = g.degrees();
    for &v in &vars {
        if (df[v] == 0) != (dg[v] == 0) {
            let (with, without, swapped) = if df[v] > 0 { (f, g, false) } else { (g, f, true) };
            let mut h = without.clone();
            for (_, c) in with.coeffs_in(v) {
                h = gcd_mz(&h, &c).0;
                if h.is_one() {
                    break;
                }
            }
            let h = if h.lc().is_negative() { h.neg() } else { h };
            let a = with.div_exact(&h).expect("gcd divides");
            let b = without.div_exact(&h).expect("gcd divides");
            return if swapped { (h, b, a) } else { (h, a, b) };
        }
    }
    if use_modular(f, g, &vars) {
        let cont = f.content().gcd(&g.content());
        let fp = f.div_coeff(&cont).unwrap();
        let gp = g.div_coeff(&cont).unwrap();
        if let Some(h) = super::modgcd::modular_gcd_mz(&fp, &gp) {
            let h = h.scale(&cont);
            let h = if h.lc().is_negative() { h.neg() } else { h };
            let a = f.div_exact(&h).expect("gcd divides f");
            let b = g.div_exact(&h).expect("gcd divides g");
            return (h, a, b);
        }
    }
    if let Some(r) = heugcd_mz(f, g) {
        return normalize_sign(r.0, r.1, r.2);
    }
    let h = prs_gcd(f, g, &vars);
    let h = if h.lc().is_negative() { h.neg() } else { h };
    let a = f.div_exact(&h).expect("gcd divides f");
    let b = g.div_exact(&h).expect("gcd divides g");
    (h, a, b)
}

/// Heuristic choice between the heuristic gcd (cheap for small inputs) and
/// the modular gcd (robust against coefficient growth).
fn use_modular(f: &MPoly<Int>, g: &MPoly<Int>, vars: &[usize]) -> bool {
    if std::env::var_os("SYMOX_HEUGCD").is_some() {
        return false;
    }
    if vars.len() < 2 {
        return false;
    }
    if std::env::var_os("SYMOX_MODGCD").is_some() {
        return true;
    }
    // predicted size (bits) of the integers at the bottom of the heuristic
    // gcd recursion: every evaluated variable multiplies it by (1 + degree)
    let norm = f.max_norm().max(g.max_norm());
    let mut bits = norm.bits() as f64 + 2.0;
    let df = f.degrees();
    let dg = g.degrees();
    for &v in &vars[1..] {
        bits *= 1.0 + df[v].max(dg[v]) as f64;
    }
    bits > 6000.0
}

/// Gcd only.
pub fn gcd_mz_only(f: &MPoly<Int>, g: &MPoly<Int>) -> MPoly<Int> {
    gcd_mz(f, g).0
}

/// ξ-adic reconstruction in variable `v` of a polynomial image.
fn interpolate_mpoly(h: &MPoly<Int>, x: &Int, v: usize) -> MPoly<Int> {
    let n = h.nvars;
    let mut out: Vec<(Vec<u32>, Int)> = Vec::new();
    let mut h = h.clone();
    let mut i = 0u32;
    while !h.is_zero() {
        let g = h.map_coeffs(|c| c.symmetric_mod(x));
        for (e, c) in g.iter() {
            let mut ne = e.to_vec();
            ne[v] = i;
            out.push((ne, c.clone()));
        }
        h = h.sub(&g).map_coeffs(|c| c.div_exact_unchecked(x));
        i += 1;
    }
    MPoly::from_terms(n, out)
}

fn prim_part(p: &MPoly<Int>) -> MPoly<Int> {
    p.primitive().1
}

fn heugcd_mz(f: &MPoly<Int>, g: &MPoly<Int>) -> Option<(MPoly<Int>, MPoly<Int>, MPoly<Int>)> {
    let cont = f.content().gcd(&g.content());
    let f = f.div_coeff(&cont).unwrap();
    let g = g.div_coeff(&cont).unwrap();
    let vars = used_union(&f, &g);
    // eliminate the variable of smallest degree last: evaluate the one of
    // highest degree (fewest recursive images to interpolate is not
    // essential; any choice is correct)
    let v = *vars.last().unwrap();
    let f_norm = f.max_norm();
    let g_norm = g.max_norm();
    let mut x = initial_xi(&f_norm, &g_norm, &f.lc(), &g.lc());
    let fin = |h: MPoly<Int>, a: MPoly<Int>, b: MPoly<Int>| (h.scale(&cont), a, b);
    for _ in 0..HEU_GCD_MAX {
        let ff = f.eval(v, &x);
        let gg = g.eval(v, &x);
        if !ff.is_zero() && !gg.is_zero() {
            let (h, cff, cfg) = gcd_mz(&ff, &gg);
            let hp = prim_part(&interpolate_mpoly(&h, &x, v));
            if !hp.is_zero() {
                if let Some(a) = f.div_exact(&hp) {
                    if let Some(b) = g.div_exact(&hp) {
                        return Some(fin(hp, a, b));
                    }
                }
            }
            let cffp = interpolate_mpoly(&cff, &x, v);
            if !cffp.is_zero() {
                if let Some(h2) = f.div_exact(&cffp) {
                    if let Some(b) = g.div_exact(&h2) {
                        return Some(fin(h2, cffp, b));
                    }
                }
            }
            let cfgp = interpolate_mpoly(&cfg, &x, v);
            if !cfgp.is_zero() {
                if let Some(h2) = g.div_exact(&cfgp) {
                    if let Some(a) = f.div_exact(&h2) {
                        return Some(fin(h2, a, cfgp));
                    }
                }
            }
        }
        x = next_xi(&x);
    }
    None
}

/// Pseudo-remainder of `a` by `b` with respect to `x_v`.
pub fn prem_mz(a: &MPoly<Int>, b: &MPoly<Int>, v: usize) -> MPoly<Int> {
    let n = a.nvars;
    let db = b.degree(v);
    let lb = b.lc_in(v);
    let mut r = a.clone();
    let mut steps = (a.degree(v) - db + 1).max(0) as u32;
    while !r.is_zero() && r.degree(v) >= db {
        let dr = r.degree(v);
        let lr = r.lc_in(v);
        let mut e = vec![0u32; n];
        e[v] = (dr - db) as u32;
        r = r.mul(&lb).sub(&lr.mul(b).mul_term(&e, &Int::one()));
        steps = steps.saturating_sub(1);
    }
    if steps > 0 {
        r = r.mul(&lb.pow(steps));
    }
    r
}

/// Content with respect to `x_v` (gcd of the coefficients in `x_v`).
pub fn content_in(f: &MPoly<Int>, v: usize) -> MPoly<Int> {
    let mut h: Option<MPoly<Int>> = None;
    for (_, c) in f.coeffs_in(v) {
        h = Some(match h {
            None => c,
            Some(h) => gcd_mz(&h, &c).0,
        });
        if h.as_ref().unwrap().is_one() {
            break;
        }
    }
    let h = h.unwrap_or_else(|| MPoly::zero(f.nvars));
    if h.lc().is_negative() {
        h.neg()
    } else {
        h
    }
}

/// Primitive part with respect to `x_v`.
pub fn primitive_in(f: &MPoly<Int>, v: usize) -> (MPoly<Int>, MPoly<Int>) {
    let c = content_in(f, v);
    if c.is_one() || c.is_zero() {
        return (c, f.clone());
    }
    let p = f.div_exact(&c).expect("content divides");
    (c, p)
}

fn prs_gcd(f: &MPoly<Int>, g: &MPoly<Int>, vars: &[usize]) -> MPoly<Int> {
    let v = vars[0];
    let (cf, pf) = primitive_in(f, v);
    let (cg, pg) = primitive_in(g, v);
    let c = gcd_mz(&cf, &cg).0;
    let (mut a, mut b) = if pf.degree(v) >= pg.degree(v) { (pf, pg) } else { (pg, pf) };
    while !b.is_zero() && b.degree(v) > 0 {
        let r = prem_mz(&a, &b, v);
        a = b;
        b = if r.is_zero() { r } else { primitive_in(&r, v).1 };
    }
    let h = if b.is_zero() { primitive_in(&a, v).1 } else { MPoly::one(f.nvars) };
    h.mul(&c)
}

// ----------------------------------------------------------------------
// over QQ
// ----------------------------------------------------------------------

/// Monic gcd over `QQ` with cofactors.
pub fn gcd_mq(f: &MPoly<Rat>, g: &MPoly<Rat>) -> (MPoly<Rat>, MPoly<Rat>, MPoly<Rat>) {
    let n = f.nvars;
    if f.is_zero() && g.is_zero() {
        return (MPoly::zero(n), MPoly::zero(n), MPoly::zero(n));
    }
    let (_, fi) = f.clear_denoms();
    let (_, gi) = g.clear_denoms();
    let (h, _, _) = gcd_mz(&fi, &gi);
    let h = h.to_rat().monic();
    let a = f.div_exact(&h).expect("gcd divides");
    let b = g.div_exact(&h).expect("gcd divides");
    (h, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zi(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&x| Int::S(x)).collect()
    }

    fn p(n: usize, t: &[(&[u32], i64)]) -> MPoly<Int> {
        MPoly::from_terms(n, t.iter().map(|(e, c)| (e.to_vec(), Int::S(*c))).collect())
    }

    #[test]
    fn univariate() {
        let f = dense::mul_int(&zi(&[1, 1]), &zi(&[-3, 0, 2]));
        let g = dense::mul_int(&zi(&[1, 1]), &zi(&[5, 7]));
        let (h, a, b) = gcd_zx(&f, &g);
        assert_eq!(h, zi(&[1, 1]));
        assert_eq!(a, zi(&[-3, 0, 2]));
        assert_eq!(b, zi(&[5, 7]));
        let (h, _, _) = modular_gcd_zx(&f, &g);
        assert_eq!(h, zi(&[1, 1]));
        let (h, _, _) = gcd_zx(&zi(&[2, 2]), &zi(&[4, 0, -4]));
        assert_eq!(h, zi(&[2, 2]));
    }

    #[test]
    fn multivariate() {
        let xy = p(2, &[(&[1, 0], 1), (&[0, 1], 1)]);
        let a = p(2, &[(&[2, 0], 1), (&[0, 1], 3), (&[0, 0], 1)]);
        let b = p(2, &[(&[1, 1], 2), (&[0, 0], -7)]);
        let f = xy.mul(&a);
        let g = xy.mul(&b);
        let (h, ca, cb) = gcd_mz(&f, &g);
        assert_eq!(h, xy);
        assert_eq!(ca, a);
        assert_eq!(cb, b);
        let h2 = prs_gcd(&f, &g, &[0, 1]);
        assert_eq!(h2, xy);
        // with monomial content
        let f = f.mul_term(&[2, 1], &Int::S(3));
        let g = g.mul_term(&[1, 3], &Int::S(6));
        let (h, _, _) = gcd_mz(&f, &g);
        assert_eq!(h, xy.mul_term(&[1, 1], &Int::S(3)));
    }
}
