//! Roots of univariate polynomials.
//!
//! [`roots`] follows SymPy's `roots()`: zero roots, linear and binomial
//! shortcuts, factorisation over `QQ`, functional decomposition, and the
//! closed forms for quadratics, cubics (Cardano), quartics (Euler /
//! Ferrari) and cyclotomic factors. [`nroots`] computes all complex roots
//! numerically using precision-preserving Aberth–Ehrlich iteration. The
//! machine-float helper remains available for internal numerical estimates.

use num_complex::Complex64 as C64;

use super::api::{cancel, factor};
use super::convert::{PResult, PolyError};
use super::dense;
use super::factor_uni::cyclotomic;
use super::int::{Field, Int, Rat};
use super::polyobj::Poly;
use crate::expand::expand;
use crate::expr::{Expr, Kind};
use crate::functions::{cos, sin};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
/// `incomplete` is the reference `solve(..., incomplete=True)` flag: a
/// polynomial whose `roots` are incomplete and whose `all_roots` is not
/// supported (symbolic coefficients) yields the partial `roots` result. It is
/// off outside a public `solve` call, so internal callers still see `None`.
pub struct RadicalOptions{pub cubics:bool,pub quartics:bool,pub incomplete:bool}
impl Default for RadicalOptions{fn default()->Self{Self{cubics:true,quartics:true,incomplete:false}}}
thread_local!{static RADICAL_OPTIONS:std::cell::Cell<RadicalOptions>=std::cell::Cell::new(RadicalOptions{cubics:true,quartics:true,incomplete:false});}

pub fn radical_options()->RadicalOptions{RADICAL_OPTIONS.with(std::cell::Cell::get)}

/// Lexical native computation only: options cannot leak across calls, threads
/// or unwinding, and callers leave the scope before invoking Python callbacks.
pub fn with_radical_options<T>(options:RadicalOptions,compute:impl FnOnce()->T)->T{
    struct Restore(RadicalOptions);
    impl Drop for Restore{fn drop(&mut self){RADICAL_OPTIONS.with(|slot|slot.set(self.0));}}
    let _restore=Restore(RADICAL_OPTIONS.with(|slot|slot.replace(options)));
    compute()
}

fn is_neg(e: &Expr) -> bool {
    match e.as_number() {
        Some(n) => n.is_negative(),
        None => match crate::evalf::eval_c64(e) {
            Some(z) if e.mask() == 0 => z.im == 0.0 && z.re < 0.0,
            _ => false,
        },
    }
}

fn is_pos(e: &Expr) -> bool {
    match e.as_number() {
        Some(n) => n.is_positive(),
        None => match crate::evalf::eval_c64(e) {
            Some(z) if e.mask() == 0 => z.im.abs() <= 1e-15 * z.re.abs() && z.re > 0.0,
            _ => false,
        },
    }
}

fn is_real_num(e: &Expr) -> bool {
    e.mask() == 0 && crate::evalf::eval_c64(e).is_some_and(|z| z.im == 0.0)
}

fn root_n(b: &Expr, n: i64) -> Expr {
    pow(b, &Expr::rational(1, n))
}

/// Coefficients (highest degree first) as expressions, plus whether the
/// domain is numerical (no symbolic coefficients).
struct UPoly {
    c: Vec<Expr>,
    numerical: bool,
}

impl UPoly {
    fn deg(&self) -> usize {
        self.c.len() - 1
    }
    fn nth(&self, k: usize) -> Expr {
        self.c[self.deg() - k].clone()
    }
    fn length(&self) -> usize {
        self.c.iter().filter(|x| !x.is_zero()).count()
    }
    fn monic(&self) -> Vec<Expr> {
        let l = self.c[0].clone();
        self.c.iter().map(|x| simp(&div(x, &l), self.numerical)).collect()
    }
}

fn simp(e: &Expr, numerical: bool) -> Expr {
    if numerical {
        e.clone()
    } else {
        factor(e).unwrap_or_else(|_| e.clone())
    }
}

fn upoly_from_rats(c: &[Rat]) -> UPoly {
    // c low -> high
    UPoly { c: c.iter().rev().map(|r| Expr::num(r.to_number())).collect(), numerical: true }
}

/// The reference `CRootOf._roots_trivial` with radicals: the linear,
/// quadratic and binomial roots of a numerical polynomial (coefficients
/// highest degree first), `None` for any other polynomial.
pub fn roots_trivial(c: &[Expr]) -> Option<Vec<Expr>> {
    let f = UPoly { c: c.to_vec(), numerical: true };
    if c.len() < 2 || c[0].is_zero() {
        return None;
    }
    match f.deg() {
        1 => Some(roots_linear(&f)),
        2 => Some(roots_quadratic(&f)),
        _ if f.length() == 2 && !f.nth(0).is_zero() => Some(roots_binomial(&f)),
        _ => None,
    }
}

// ----------------------------------------------------------------------
// closed forms
// ----------------------------------------------------------------------

fn roots_linear(f: &UPoly) -> Vec<Expr> {
    let r = neg(&div(&f.nth(0), &f.nth(1)));
    vec![simp(&r, f.numerical)]
}

fn sqrt_squares(d: &Expr) -> Expr {
    // pull even powers out of the square root
    let args: Vec<Expr> = match d.kind() {
        Kind::Mul(..) => d.args(),
        _ => vec![d.clone()],
    };
    let mut co = Vec::new();
    let mut other = Vec::new();
    for a in args {
        if let Kind::Pow(b, x) = a.kind() {
            if let Some(k) = x.as_i64() {
                if k % 2 == 0 {
                    co.push(pow(b, &Expr::int(k / 2)));
                    continue;
                }
            }
        }
        other.push(a);
    }
    if co.is_empty() {
        return sqrt(d);
    }
    mul(&[mul(&co), sqrt(&mul(&other))])
}

/// The reference `expand_2arg`: distribute numeric coefficients over sums,
/// bottom up, and nothing else.
fn expand_2arg(e: &Expr) -> Expr {
    crate::simptools::bottom_up(e, &mut |t: &Expr| {
        if t.is_mul() {
            let (c, r) = crate::ops::as_coeff_mul(t);
            if r.is_add() {
                let c = Expr::num(c);
                return add(&r.args().iter().map(|x| mul(&[c.clone(), x.clone()])).collect::<Vec<_>>());
            }
        }
        t.clone()
    })
}

fn roots_quadratic(f: &UPoly) -> Vec<Expr> {
    let (a, b, c) = (f.c[0].clone(), f.c[1].clone(), f.c[2].clone());
    let num = f.numerical;
    if c.is_zero() {
        let r0 = Expr::zero();
        let r1 = neg(&div(&b, &a));
        if !num {
            return vec![r0, simp(&r1, false)];
        }
        if is_neg(&r1) {
            return vec![r1, r0];
        }
        return vec![r0, r1];
    }
    if b.is_zero() {
        let r = simp(&neg(&div(&c, &a)), num);
        let rr = sqrt_squares(&r);
        return vec![neg(&rr), rr];
    }
    let mut d = sub(&pow(&b, &Expr::int(2)), &mul(&[Expr::int(4), a.clone(), c.clone()]));
    let aa = mul(&[Expr::int(2), a.clone()]);
    let mut bb = neg(&div(&b, &aa));
    if !num {
        d = simp(&d, false);
        bb = simp(&bb, false);
    }
    // The reference `factor_terms(_sqrt(d)/A)`.
    let dd = crate::simptools::exprtools::factor_terms(&div(&sqrt_squares(&d), &aa), false, false, false, true);
    let mut r0 = sub(&bb, &dd);
    let mut r1 = add(&[bb, dd]);
    if is_neg(&a) {
        std::mem::swap(&mut r0, &mut r1);
    } else if !num {
        r0 = expand_2arg(&r0);
        r1 = expand_2arg(&r1);
    }
    vec![r0, r1]
}

/// Exact first quadrant values `(cos, sin)` of `p*pi/q` for the angles
/// SymPy evaluates but the core does not.
fn trig_table(p: i64, q: i64) -> Option<(Expr, Expr)> {
    let r = |a: i64, b: i64| Expr::rational(a, b);
    let s5 = sqrt(&Expr::int(5));
    let s2 = sqrt(&Expr::int(2));
    let a = add(&[r(1, 4), mul(&[r(1, 4), s5.clone()])]); // 1/4 + sqrt(5)/4
    let b = add(&[r(-1, 4), mul(&[r(1, 4), s5.clone()])]); // -1/4 + sqrt(5)/4
    let c = sqrt(&add(&[r(5, 8), mul(&[r(-1, 8), s5.clone()])])); // sqrt(5/8 - sqrt(5)/8)
    let d = sqrt(&add(&[r(5, 8), mul(&[r(1, 8), s5])])); // sqrt(sqrt(5)/8 + 5/8)
    let e = sqrt(&add(&[r(1, 2), mul(&[r(1, 4), s2.clone()])])); // sqrt(sqrt(2)/4 + 1/2)
    let f = sqrt(&add(&[r(1, 2), mul(&[r(-1, 4), s2])])); // sqrt(1/2 - sqrt(2)/4)
    Some(match (p, q) {
        (1, 5) => (a, c),
        (2, 5) => (b, d),
        (1, 10) => (d, b),
        (3, 10) => (c, a),
        (1, 8) => (e, f),
        (3, 8) => (f, e),
        _ => return None,
    })
}

/// `(cos(k*pi/n), sin(k*pi/n))` with exact values where known.
fn cos_sin_pi(k: i64, n: i64) -> (Expr, Expr) {
    // reduce k/n modulo 2 into [0, 2)
    let g = gcd_i(k, n).max(1);
    let (mut k, n) = (k / g, n / g);
    let n = if n < 0 { k = -k; -n } else { n };
    k = k.rem_euclid(2 * n);
    let (mut sc, mut ss) = (1i64, 1i64);
    if k > n {
        k = 2 * n - k;
        ss = -1;
    }
    // now angle in [0, pi]
    if 2 * k > n {
        k = n - k;
        sc = -1;
    }
    let g = gcd_i(k, n).max(1);
    let (p, q) = (k / g, n / g);
    let (c, s) = match trig_table(p, q) {
        Some(v) => v,
        None => {
            let ang = mul(&[Expr::rational(p, q), Expr::pi()]);
            (cos(&ang), sin(&ang))
        }
    };
    (mul(&[Expr::int(sc), c]), mul(&[Expr::int(ss), s]))
}

/// `exp(2*pi*I*k/n)` expanded (exact cos/sin values where known).
fn unit_root(k: i64, n: i64) -> Expr {
    let (c, s) = cos_sin_pi(2 * k, n);
    add(&[c, mul(&[Expr::i(), s])])
}

/// `root(base, n)` expanded into real and imaginary parts for real bases.
fn principal_root(base: &Expr, n: i64) -> Expr {
    if is_neg(base) && is_real_num(base) {
        let r = root_n(&neg(base), n);
        let (c, s) = cos_sin_pi(1, n);
        return expand(&mul(&[r, add(&[c, mul(&[Expr::i(), s])])]));
    }
    let r = root_n(base, n);
    if !is_real_num(base) {
        // the reference's alpha.expand(complex=True)
        let options = crate::expand::ExpandOptions::default();
        return crate::complex_parts::expand_complex(&r, &options).unwrap_or(r);
    }
    r
}

fn roots_binomial(f: &UPoly) -> Vec<Expr> {
    let n = f.deg() as i64;
    let a = f.nth(n as usize);
    let b = f.nth(0);
    let base = neg(&cancel(&div(&b, &a)).unwrap_or_else(|_| div(&b, &a)));
    let alpha = if base.mask() == 0 { principal_root(&base, n) } else { root_n(&base, n) };
    let negb = is_neg(&base);
    let even = n % 2 == 0;
    let big = negb && even && is_pos(&add(&[base.clone(), Expr::one()]));
    let mut ks: Vec<i64> = Vec::new();
    let mut imax = n / 2;
    if even {
        ks.push(imax);
        imax -= 1;
    }
    if !negb {
        ks.push(0);
    }
    for i in (1..=imax).rev() {
        if negb {
            ks.extend([i, -i]);
        } else {
            ks.extend([-i, i]);
        }
    }
    if negb {
        ks.push(0);
        if big {
            // reversed pairs (no effect on the values, kept for parity)
        }
    }
    ks.into_iter()
        .map(|k| {
            let zeta = unit_root(k, n);
            expand(&mul(&[alpha.clone(), zeta]))
        })
        .collect()
}

fn gcd_i(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        gcd_i(b, a % b)
    }
}

fn totient(n: i64) -> i64 {
    (1..=n).filter(|&k| gcd_i(k, n) == 1).count() as i64
}

/// Index `n` if `f` (integer coefficients, low -> high) is `Φ_n`.
fn cyclotomic_index(f: &[Int]) -> Option<i64> {
    let d = f.len() as i64 - 1;
    if d < 1 || !f.last().unwrap().is_one() {
        return None;
    }
    let mut memo = Default::default();
    for n in 1..=(6 * d * d + 10).min(100000) {
        if totient(n) == d && cyclotomic(n as u64, &mut memo) == f {
            return Some(n);
        }
    }
    None
}

fn roots_cyclotomic(n: i64) -> Vec<Expr> {
    let h = n / 2;
    let mut ks: Vec<i64> = (1..=n).filter(|&i| gcd_i(i, n) == 1).collect();
    ks.sort_by_key(|&x| if x <= h { (x, -1) } else { ((x - n).abs(), 1) });
    ks.iter().rev().map(|&k| unit_root(k, n)).collect()
}

fn roots_cubic(f: &UPoly) -> Vec<Expr> {
    let m = f.monic();
    let (a, b, c) = (m[1].clone(), m[2].clone(), m[3].clone());
    if c.is_zero() {
        let rs = roots_list(&UPoly { c: vec![Expr::one(), a.clone(), b.clone()], numerical: f.numerical });
        if rs.len() == 2 {
            return vec![rs[0].clone(), Expr::zero(), rs[1].clone()];
        }
    }
    let third = Expr::rational(1, 3);
    let p = sub(&b, &div(&pow(&a, &Expr::int(2)), &Expr::int(3)));
    let q = add(&[
        c.clone(),
        neg(&div(&mul(&[a.clone(), b.clone()]), &Expr::int(3))),
        mul(&[Expr::rational(2, 27), pow(&a, &Expr::int(3))]),
    ]);
    let pon3 = div(&p, &Expr::int(3));
    let aon3 = div(&a, &Expr::int(3));
    let mut u1: Option<Expr> = None;
    if p.is_zero() {
        if q.is_zero() {
            return vec![neg(&aon3); 3];
        }
        u1 = Some(if is_pos(&q) { neg(&root_n(&q, 3)) } else { root_n(&neg(&q), 3) });
    } else if q.is_zero() {
        let rs = roots_list(&UPoly { c: vec![Expr::one(), Expr::zero(), p.clone()], numerical: f.numerical });
        if rs.len() == 2 {
            return [rs[0].clone(), Expr::zero(), rs[1].clone()].iter().map(|t| sub(t, &aon3)).collect();
        }
    } else if is_real_num(&q) && is_neg(&q) {
        let inner = add(&[neg(&div(&q, &Expr::int(2))), sqrt(&add(&[div(&pow(&q, &Expr::int(2)), &Expr::int(4)), pow(&pon3, &Expr::int(3))]))]);
        u1 = Some(neg(&root_n(&inner, 3)));
    }
    let coeff = div(&mul(&[Expr::i(), sqrt(&Expr::int(3))]), &Expr::int(2));
    let mh = Expr::rational(-1, 2);
    match u1 {
        None => {
            let us = [Expr::one(), add(&[mh.clone(), coeff.clone()]), sub(&mh, &coeff)];
            let (bb, cc, dd) = (a.clone(), b.clone(), c.clone());
            let d0 = sub(&pow(&bb, &Expr::int(2)), &mul(&[Expr::int(3), cc.clone()]));
            let d1 = add(&[
                mul(&[Expr::int(2), pow(&bb, &Expr::int(3))]),
                mul(&[Expr::int(-9), bb.clone(), cc.clone()]),
                mul(&[Expr::int(27), dd.clone()]),
            ]);
            let disc = sub(&pow(&d1, &Expr::int(2)), &mul(&[Expr::int(4), pow(&d0, &Expr::int(3))]));
            let cbig = pow(&div(&add(&[d1.clone(), sqrt(&disc)]), &Expr::int(2)), &third);
            us.iter()
                .map(|uk| {
                    let t = add(&[bb.clone(), mul(&[uk.clone(), cbig.clone()]), div(&d0, &mul(&[cbig.clone(), uk.clone()]))]);
                    neg(&div(&t, &Expr::int(3)))
                })
                .collect()
        }
        Some(u1) => {
            let u2 = mul(&[u1.clone(), add(&[mh.clone(), coeff.clone()])]);
            let u3 = mul(&[u1.clone(), sub(&mh, &coeff)]);
            if p.is_zero() {
                return vec![sub(&u1, &aon3), sub(&u2, &aon3), sub(&u3, &aon3)];
            }
            [u1, u2, u3]
                .iter()
                .map(|u| add(&[neg(u), div(&pon3, u), neg(&aon3)]))
                .collect()
        }
    }
}

fn is_rational_expr(e: &Expr) -> bool {
    e.as_number().is_some_and(|n| n.is_rational())
}

fn roots_quartic_euler(p: &Expr, q: &Expr, r: &Expr, a: &Expr) -> Option<Vec<Expr>> {
    // 64 R^3 + 32 p R^2 + (4 p^2 - 16 r) R - q^2 = 0, rational roots only
    if ![p, q, r].iter().all(|x| is_rational_expr(x)) {
        return None;
    }
    let pr = super::convert::number_to_rat(p.as_number()?);
    let qr = super::convert::number_to_rat(q.as_number()?);
    let rr = super::convert::number_to_rat(r.as_number()?);
    let c3 = Rat::from_i64(64);
    let c2 = pr.mul(&Rat::from_i64(32));
    let c1 = pr.mul(&pr).mul(&Rat::from_i64(4)).sub(&rr.mul(&Rat::from_i64(16)));
    let c0 = qr.mul(&qr).neg();
    let rs = rational_roots(&[c0, c1, c2, c3]);
    let rs: Vec<Rat> = rs.into_iter().filter(|x| !x.is_zero()).collect();
    let big = rs.into_iter().max()?;
    let rx = Expr::num(big.to_number());
    let c1 = sqrt(&rx);
    let bq = div(&mul(&[neg(q), c1.clone()]), &mul(&[Expr::int(4), rx.clone()]));
    let aq = sub(&neg(&rx), &div(p, &Expr::int(2)));
    let c2 = sqrt(&add(&[aq.clone(), bq.clone()]));
    let c3 = sqrt(&sub(&aq, &bq));
    Some(vec![
        sub(&sub(&c1, &c2), a),
        sub(&sub(&neg(&c1), &c3), a),
        sub(&add(&[neg(&c1), c3]), a),
        sub(&add(&[c1, c2]), a),
    ])
}

fn rational_roots(c: &[Rat]) -> Vec<Rat> {
    // roots of an integer-cleared polynomial via factorisation
    let (_, z) = dense::clear_denoms(c);
    let (_, fs) = super::factor_uni::factor_zx(&z);
    let mut out = Vec::new();
    for (g, _) in fs {
        if g.len() == 2 {
            out.push(Rat::new(g[0].neg(), g[1].clone()));
        }
    }
    out
}

fn roots_quartic(f: &UPoly) -> Vec<Expr> {
    let m = f.monic();
    let (a, b, c, d) = (m[1].clone(), m[2].clone(), m[3].clone(), m[4].clone());
    let num = f.numerical;
    if d.is_zero() {
        let mut out = vec![Expr::zero()];
        out.extend(roots_list(&UPoly { c: vec![Expr::one(), a, b, c], numerical: num }));
        return out;
    }
    if !a.is_zero() && expand(&sub(&pow(&div(&c, &a), &Expr::int(2)), &d)).is_zero() {
        let mm = div(&c, &a);
        let g = UPoly { c: vec![Expr::one(), a.clone(), sub(&b, &mul(&[Expr::int(2), mm.clone()]))], numerical: num };
        let z = roots_quadratic(&g);
        let h1 = UPoly { c: vec![Expr::one(), neg(&z[0]), mm.clone()], numerical: num };
        let h2 = UPoly { c: vec![Expr::one(), neg(&z[1]), mm.clone()], numerical: num };
        let mut r = roots_quadratic(&h1);
        r.extend(roots_quadratic(&h2));
        return r;
    }
    let a2 = pow(&a, &Expr::int(2));
    let e = sub(&b, &div(&mul(&[Expr::int(3), a2.clone()]), &Expr::int(8)));
    let ff = expand(&add(&[c.clone(), mul(&[a.clone(), sub(&div(&a2, &Expr::int(8)), &div(&b, &Expr::int(2)))])]));
    let aon4 = div(&a, &Expr::int(4));
    let g = expand(&sub(
        &d,
        &mul(&[aon4.clone(), add(&[mul(&[a.clone(), sub(&div(&mul(&[Expr::int(3), a2.clone()]), &Expr::int(64)), &div(&b, &Expr::int(4)))]), c.clone()])]),
    ));
    if ff.is_zero() {
        let ys = roots_list(&UPoly { c: vec![Expr::one(), e.clone(), g.clone()], numerical: num });
        if ys.len() == 2 {
            let y1 = sqrt(&ys[0]);
            let y2 = sqrt(&ys[1]);
            return [neg(&y1), neg(&y2), y1, y2].iter().map(|t| sub(t, &aon4)).collect();
        }
    }
    if g.is_zero() {
        let mut y = vec![Expr::zero()];
        y.extend(roots_list(&UPoly { c: vec![Expr::one(), Expr::zero(), e.clone(), ff.clone()], numerical: num }));
        return y.iter().map(|t| sub(t, &aon4)).collect();
    }
    if let Some(s) = roots_quartic_euler(&e, &ff, &g, &aon4) {
        return s;
    }
    let p = expand(&sub(&neg(&div(&pow(&e, &Expr::int(2)), &Expr::int(12))), &g));
    let q = add(&[
        neg(&div(&pow(&e, &Expr::int(3)), &Expr::int(108))),
        div(&mul(&[e.clone(), g.clone()]), &Expr::int(3)),
        neg(&div(&pow(&ff, &Expr::int(2)), &Expr::int(8))),
    ]);
    let th = Expr::rational(1, 3);
    let ans = |y: &Expr| -> Vec<Expr> {
        let w = sqrt(&add(&[e.clone(), mul(&[Expr::int(2), y.clone()])]));
        let arg1 = add(&[mul(&[Expr::int(3), e.clone()]), mul(&[Expr::int(2), y.clone()])]);
        let arg2 = div(&mul(&[Expr::int(2), ff.clone()]), &w);
        let mut out = Vec::new();
        for s in [-1i64, 1] {
            let root = sqrt(&neg(&add(&[arg1.clone(), mul(&[Expr::int(s), arg2.clone()])])));
            for t in [-1i64, 1] {
                let v = div(&sub(&mul(&[Expr::int(s), w.clone()]), &mul(&[Expr::int(t), root.clone()])), &Expr::int(2));
                out.push(sub(&v, &aon4));
            }
        }
        out
    };
    let y1 = sub(&mul(&[e.clone(), Expr::rational(-5, 6)]), &pow(&q, &th));
    if p.is_zero() {
        return ans(&y1);
    }
    let root = sqrt(&add(&[div(&pow(&q, &Expr::int(2)), &Expr::int(4)), div(&pow(&p, &Expr::int(3)), &Expr::int(27))]));
    let r = add(&[neg(&div(&q, &Expr::int(2))), root]);
    let u = pow(&r, &th);
    let y2 = add(&[mul(&[e.clone(), Expr::rational(-5, 6)]), u.clone(), neg(&div(&div(&p, &u), &Expr::int(3)))]);
    ans(&y2)
}

/// `roots(..., multiple=True)` for an expression-coefficient polynomial
/// (used recursively by the closed forms).
fn roots_list(f: &UPoly) -> Vec<Expr> {
    let mut out = Vec::new();
    for (r, k) in roots_of(f) {
        for _ in 0..k {
            out.push(r.clone());
        }
    }
    out
}

/// Functional decomposition `f = g1(g2(...(gk)))` over `QQ` (Kozen–Landau).
pub fn decompose(f: &[Rat]) -> Vec<Vec<Rat>> {
    decompose_field(f,&|_|true).expect("rational coefficients have characteristic zero")
}

/// Bound dense allocation before converting a sparse polynomial.
pub fn decompose_preflight(variables:usize,degree:i64)->PResult<()> {
    if variables!=1 {return Err(PolyError::Unsupported("univariate polynomial expected".into()));}
    if degree>4096 {return Err(PolyError::Unsupported("functional decomposition is limited to degree 4096".into()));}
    Ok(())
}

/// Kozen--Landau decomposition over an exact field, optionally restricted to
/// coefficients that belong to its declared coefficient ring. Every accepted
/// pair is certified by exact polynomial division, without numeric sampling.
pub fn decompose_field<C:Field>(f:&[C],allow:&impl Fn(&C)->bool)->PResult<Vec<Vec<C>>> {
    fn right<C:Field>(f: &[C], s: usize,allow:&impl Fn(&C)->bool) -> PResult<Option<Vec<C>>> {
        let n = f.len() - 1;
        let lc = f[n].clone();
        let r = n / s;
        let mut g = vec![C::zero(); s + 1];
        g[s] = C::one();
        for i in 1..s {
            let mut coeff = C::zero();
            for j in 0..i {
                let fc = &f[n + j - i];
                let gc = &g[s - j];
                if fc.is_zero() || gc.is_zero() {
                    continue;
                }
                let k = C::from_i64(i as i64 - (r * j) as i64);
                coeff = coeff.add(&k.mul(fc).mul(gc));
            }
            let denominator=C::from_i64((i*r) as i64).mul(&lc);
            if denominator.is_zero() {return Err(PolyError::NotInvertible("zero divisor in functional decomposition".into()));}
            let value=coeff.div(&denominator);
            if !allow(&value) {return Ok(None);}
            g[s-i]=value;
        }
        Ok(Some(g))
    }
    fn left<C:Field>(f: &[C], h: &[C],allow:&impl Fn(&C)->bool) -> Option<Vec<C>> {
        let mut g = Vec::new();
        let mut f = f.to_vec();
        while !f.is_empty() {
            let (q, r) = dense::divrem(&f, h);
            if r.len() > 1 {
                return None;
            }
            let coefficient=r.first().cloned().unwrap_or_else(C::zero);
            if !allow(&coefficient) {return None;}
            g.push(coefficient);
            f = q;
        }
        dense::trim(&mut g);
        Some(g)
    }
    fn step<C:Field>(f: &[C],allow:&impl Fn(&C)->bool) -> PResult<Option<(Vec<C>, Vec<C>)>> {
        if f.len()<3 {return Ok(None);}
        let df = f.len() - 1;
        for s in 2..df {
            if df % s != 0 {
                continue;
            }
            let Some(h) = right(f, s,allow)? else {continue};
            if let Some(g) = left(f, &h,allow) {
                return Ok(Some((g, h)));
            }
        }
        Ok(None)
    }
    let mut fs: Vec<Vec<C>> = Vec::new();
    let mut f = f.to_vec();
    while let Some((g, h)) = step(&f,allow)? {
        f = g;
        fs.insert(0, h);
    }
    fs.insert(0, f);
    Ok(fs)
}

fn rats_of(f: &UPoly) -> Option<Vec<Rat>> {
    f.c.iter()
        .rev()
        .map(|e| e.as_number().filter(|n| n.is_rational()).map(super::convert::number_to_rat))
        .collect()
}

fn try_heuristics(f: &UPoly) -> Vec<Expr> {
    let n = f.deg();
    if n == 0 {
        return vec![];
    }
    if f.length() == 1 {
        return vec![Expr::zero(); n];
    }
    if f.length() == 2 {
        if n == 1 {
            return roots_linear(f).iter().map(|r| cancel(r).unwrap_or_else(|_| r.clone())).collect();
        }
        return roots_binomial(f);
    }
    let mut result = Vec::new();
    let mut f = UPoly { c: f.c.clone(), numerical: f.numerical };
    if let Some(rs) = rats_of(&f) {
        for i in [-1i64, 1] {
            let v = dense::eval(&rs, &Rat::from_i64(i));
            if v.is_zero() {
                let (q, _) = dense::divrem(&rs, &[Rat::from_i64(-i), Rat::one()]);
                f = upoly_from_rats(&q);
                result.push(Expr::int(i));
                break;
            }
        }
    }
    let n = f.deg();
    if n == 1 {
        result.extend(roots_linear(&f).iter().map(|r| cancel(r).unwrap_or_else(|_| r.clone())));
    } else if n == 2 {
        result.extend(roots_quadratic(&f).iter().map(|r| cancel(r).unwrap_or_else(|_| r.clone())));
    } else if let Some(idx) = rats_of(&f).and_then(|rs| {
        let (_, z) = dense::clear_denoms(&rs);
        let (_, pz) = dense::primitive(&z);
        if rs.last().is_some_and(|l| l.is_one()) {
            cyclotomic_index(&pz)
        } else {
            None
        }
    }) {
        result.extend(roots_cyclotomic(idx));
    } else if n == 3 && radical_options().cubics {
        result.extend(roots_cubic(&f));
    } else if n == 4 && radical_options().quartics {
        // Enabling a quartic formula includes its internal cubic resolvent.
        result.extend(with_radical_options(RadicalOptions::default(),||roots_quartic(&f)));
    }
    result
}

fn try_decompose(f: &[Rat]) -> Vec<Expr> {
    let fs = decompose(f);
    let mut roots = try_heuristics(&upoly_from_rats(&fs[0]));
    for h in &fs[1..] {
        let prev = std::mem::take(&mut roots);
        for r in prev {
            // h - r
            let mut c: Vec<Expr> = h.iter().rev().map(|x| Expr::num(x.to_number())).collect();
            let last = c.len() - 1;
            c[last] = sub(&c[last], &r);
            let numerical = c.iter().all(|x| x.is_number());
            roots.extend(try_heuristics(&UPoly { c, numerical: numerical || r.mask() == 0 }));
        }
    }
    roots
}

/// The reference `_integer_basis` of an integer polynomial given by its
/// terms `(degree, coefficient)`, highest degree first.
fn integer_basis(terms: &[(u32, num_bigint::BigInt)]) -> Option<num_bigint::BigInt> {
    use num_integer::Integer;
    use num_traits::{One, Signed, Zero};
    let mut coeffs: Vec<num_bigint::BigInt> = terms.iter().map(|t| t.1.abs()).collect();
    let mut monoms: Vec<u32> = terms.iter().map(|t| t.0).collect();
    if coeffs[0] < coeffs[coeffs.len() - 1] {
        coeffs.reverse();
        let n = monoms[0];
        monoms = monoms.iter().rev().map(|m| n - m).collect();
    } else {
        return None;
    }
    monoms.pop();
    coeffs.pop();
    if monoms.len() == 1 {
        let (r, exact) = crate::ntheory::integer_nthroot(&coeffs[0], monoms[0]);
        return exact.then_some(r);
    }
    let g = coeffs.iter().fold(num_bigint::BigInt::zero(), |acc, c| acc.gcd(c));
    let mut divs = crate::ntheory::divisors(&g);
    divs.sort();
    divs.retain(|d| !d.is_one());
    'outer: for div in divs.into_iter().rev() {
        for (m, c) in monoms.iter().zip(&coeffs) {
            if !(c % num_traits::pow(div.clone(), *m as usize)).is_zero() {
                continue 'outer;
            }
        }
        return Some(div);
    }
    None
}

/// The reference `preprocess_roots` over a polynomial coefficient ring:
/// clear denominators, take the primitive part and, when every coefficient
/// is a term, rescale `x = gen**ratio * y` for each generator whose
/// exponents are proportional to the degrees.  Returns the root factor and
/// the preprocessed polynomial.
fn preprocess_symbolic(f: &UPoly) -> Option<(Expr, UPoly)> {
    if super::expression_domain::construct_selects_expression(&f.c).ok()? {
        return None;
    }
    let x = Expr::dummy("x", Default::default());
    let n = f.deg();
    let terms: Vec<Expr> = f.c.iter().enumerate().map(|(i, c)| mul(&[c.clone(), pow(&x, &Expr::int((n - i) as i64))])).collect();
    let (numer, _) = super::api::as_numer_denom(&add(&terms));
    let converted = super::convert::from_exprs(std::slice::from_ref(&numer)).ok()?;
    if converted.float || converted.gens.len() < 2 {
        return None;
    }
    let xi = converted.gens.iter().position(|g| *g == x)?;
    let others: Vec<usize> = (0..converted.gens.len()).filter(|&i| i != xi).collect();
    let poly = &converted.polys[0];
    // Group the terms by degree in `x`, highest first.
    let mut rows: Vec<(u32, Vec<(Vec<u32>, Rat)>)> = Vec::new();
    for (e, c) in poly.iter() {
        let rest: Vec<u32> = others.iter().map(|&i| e[i]).collect();
        match rows.iter_mut().find(|r| r.0 == e[xi]) {
            Some(r) => r.1.push((rest, c.clone())),
            None => rows.push((e[xi], vec![(rest, c.clone())])),
        }
    }
    rows.sort_by(|a, b| b.0.cmp(&a.0));
    // Integer content.
    let mut den = Int::one();
    for r in &rows {
        for (_, c) in &r.1 {
            den = den.lcm(&c.d);
        }
    }
    let mut content = Int::zero();
    for r in &rows {
        for (_, c) in &r.1 {
            content = content.gcd(&c.n.mul(&den.div_exact_unchecked(&c.d)));
        }
    }
    let scale = Rat::new(den, content);
    let mut rows: Vec<(u32, Vec<(Vec<u32>, Rat)>)> = rows.into_iter().map(|(k, t)| (k, t.into_iter().map(|(e, c)| (e, c.mul(&scale))).collect())).collect();
    let monomial = rows.iter().all(|r| r.1.len() == 1);
    // Monomial content in the coefficient generators.
    for j in 0..others.len() {
        let low = rows.iter().flat_map(|r| r.1.iter().map(move |t| t.0[j])).min().unwrap_or(0);
        for r in rows.iter_mut() {
            for t in r.1.iter_mut() {
                t.0[j] -= low;
            }
        }
    }
    let mut coeff = Expr::one();
    if monomial {
        let base: Vec<u32> = rows.iter().map(|r| r.0).collect();
        let mut eliminated = vec![false; others.len()];
        for j in 0..others.len() {
            let mut strip: Vec<u32> = rows.iter().map(|r| r.1[0].0[j]).collect();
            let reverse = strip[0] < strip[strip.len() - 1];
            if reverse {
                strip.reverse();
            }
            let mut ratio: Option<i64> = None;
            let mut complete = true;
            for (&a, &b) in base.iter().zip(&strip) {
                if a == 0 && b == 0 {
                    continue;
                }
                if a == 0 || b == 0 || b % a != 0 {
                    complete = false;
                    break;
                }
                let r = (b / a) as i64;
                if ratio.is_some_and(|q| q != r) {
                    complete = false;
                    break;
                }
                ratio = Some(r);
            }
            if complete {
                let Some(mut ratio) = ratio else { continue };
                if reverse {
                    ratio = -ratio;
                }
                eliminated[j] = true;
                coeff = mul(&[coeff, pow(&converted.gens[others[j]], &Expr::int(-ratio))]);
            }
        }
        for r in rows.iter_mut() {
            for t in r.1.iter_mut() {
                for j in 0..others.len() {
                    if eliminated[j] {
                        t.0[j] = 0;
                    }
                }
            }
        }
    }
    // The reference `_integer_basis` once the domain retracts to ZZ.
    if rows.iter().all(|r| r.1.iter().all(|t| t.0.iter().all(|&d| d == 0))) {
        let terms: Vec<(u32, num_bigint::BigInt)> = rows.iter().map(|r| (r.0, r.1.iter().fold(Int::zero(), |acc, t| acc.add(&t.1.n)).to_big())).filter(|t| !num_traits::Zero::is_zero(&t.1)).collect();
        if let Some(basis) = integer_basis(&terms) {
            let n = rows[0].0;
            for r in rows.iter_mut() {
                let power = num_traits::pow(basis.clone(), (n - r.0) as usize);
                for t in r.1.iter_mut() {
                    t.1 = Rat::new(Int::from_big(num_integer::Integer::div_floor(&t.1.n.to_big(), &power)), Int::one());
                }
            }
            coeff = mul(&[coeff, Expr::num(Rat::new(Int::from_big(basis), Int::one()).to_number())]);
        }
    }
    let mut c = vec![Expr::zero(); rows[0].0 as usize + 1];
    let top = rows[0].0 as usize;
    for (k, t) in &rows {
        let parts: Vec<Expr> = t
            .iter()
            .map(|(e, q)| {
                let mut factors = vec![Expr::num(q.to_number())];
                for (j, &d) in e.iter().enumerate() {
                    if d > 0 {
                        factors.push(pow(&converted.gens[others[j]], &Expr::int(d as i64)));
                    }
                }
                mul(&factors)
            })
            .collect();
        c[top - *k as usize] = add(&parts);
    }
    let numerical = c.iter().all(|t| t.is_rational_number());
    Some((coeff, UPoly { c, numerical }))
}

/// The reference `_try_decompose` over the rational-function field of the
/// coefficients.
fn try_decompose_symbolic(f: &UPoly) -> Vec<Expr> {
    let Some(parts) = decompose_symbolic(f) else { return try_heuristics(f) };
    let mut roots = try_heuristics(&parts[0]);
    for h in &parts[1..] {
        let prev = std::mem::take(&mut roots);
        for r in prev {
            let mut c = h.c.clone();
            let last = c.len() - 1;
            c[last] = sub(&c[last], &r);
            roots.extend(try_heuristics(&UPoly { c, numerical: false }));
        }
    }
    roots
}

fn decompose_symbolic(f: &UPoly) -> Option<Vec<UPoly>> {
    decompose_preflight(1, f.deg() as i64).ok()?;
    let parameters = super::convert::from_exprs(&f.c).ok()?.gens;
    if parameters.is_empty() {
        return None;
    }
    let x = Expr::dummy("x", Default::default());
    let n = f.deg();
    let terms: Vec<Expr> = f.c.iter().enumerate().map(|(i, c)| mul(&[c.clone(), pow(&x, &Expr::int((n - i) as i64))])).collect();
    let p = super::fraction::FractionPoly::from_expr(&add(&terms), std::slice::from_ref(&x), &parameters, "ZZ(..)").ok()?;
    let _scope = super::ratfunc::Scope::new(parameters.len());
    let coefficients = dense::from_mpoly(&p.rep, 0);
    let parts = decompose_field(&coefficients, &|_| true).ok()?;
    Some(
        parts
            .iter()
            .map(|part| UPoly {
                c: part.iter().rev().map(|c| super::fraction::coefficient_expr(c, &parameters)).collect(),
                numerical: false,
            })
            .collect(),
    )
}

/// The factors in the variable of a polynomial with symbolic coefficients
/// over a polynomial or rational-function domain, with the total number of
/// factors of the cleared numerator (the reference `factor_list` length).
fn symbolic_factors(f: &UPoly) -> Option<(usize, Vec<(UPoly, u32)>, bool)> {
    // Over EX the reference also factors `Poly(f.as_expr())`, with the
    // algebraic coefficients as generators; only an irreducible simple
    // factor keeps the EX heuristics.
    let expression = super::expression_domain::construct_selects_expression(&f.c).ok()?;
    let x = Expr::dummy("x", Default::default());
    let n = f.deg();
    let terms: Vec<Expr> = f.c.iter().enumerate().map(|(i, c)| mul(&[c.clone(), pow(&x, &Expr::int((n - i) as i64))])).collect();
    let (numer, _) = super::api::as_numer_denom(&add(&terms));
    let (_, factors) = super::api::factor_list(&numer).ok()?;
    let count = factors.len();
    let mut out = Vec::new();
    for (g, k) in factors {
        if !crate::traverse::has(&g, &x) {
            continue;
        }
        let k = u32::try_from(k.as_i64()?).ok()?;
        let c = Poly::from_expr(&g, std::slice::from_ref(&x)).ok()?.all_coeffs().ok()?;
        out.push((UPoly { c, numerical: false }, k));
    }
    Some((count, out, expression))
}

fn update(result: &mut Vec<(Expr, u32)>, zeros: &mut u32, r: Expr, k: u32) {
    if r.is_zero() {
        *zeros += k;
    }
    if let Some(x) = result.iter_mut().find(|x| x.0 == r) {
        x.1 += k;
    } else {
        result.push((r, k));
    }
}

/// Roots of a polynomial with expression coefficients (highest first).
fn roots_of(f: &UPoly) -> Vec<(Expr, u32)> {
    roots_of_with(f, true)
}

fn roots_of_with(f: &UPoly, preprocess: bool) -> Vec<(Expr, u32)> {
    let mut c = f.c.clone();
    while c.first().is_some_and(|x| x.is_zero()) {
        c.remove(0);
    }
    let mut k = 0u32;
    while c.len() > 1 && c.last().unwrap().is_zero() {
        c.pop();
        k += 1;
    }
    let mut zeros = k;
    let mut result: Vec<(Expr, u32)> = Vec::new();
    if c.len() <= 1 {
        if zeros > 0 {
            result.push((Expr::zero(), zeros));
        }
        return result;
    }
    let mut f = UPoly { c, numerical: f.numerical };
    // preprocess: primitive integer polynomial for numerical domains
    if let Some(rs) = rats_of(&f) {
        let (_, z) = dense::clear_denoms(&rs);
        let (_, pz) = dense::primitive(&z);
        f = upoly_from_rats(&dense::to_rat(&pz));
        let n = f.deg();
        if n == 1 {
            for r in roots_linear(&f) {
                update(&mut result, &mut zeros, r, 1);
            }
        } else if f.length() == 2 {
            let rs = if n == 2 { roots_quadratic(&f) } else { roots_binomial(&f) };
            for r in rs {
                update(&mut result, &mut zeros, r, 1);
            }
        } else {
            let (_, fs) = super::factor_uni::factor_zx(&pz);
            let mut fs = fs;
            fs.sort_by(|(a, ka), (b, kb)| {
                a.len().cmp(&b.len()).then(ka.cmp(kb)).then_with(|| {
                    for (x, y) in a.iter().rev().zip(b.iter().rev()) {
                        let c = x.cmp(y);
                        if c != std::cmp::Ordering::Equal {
                            return c;
                        }
                    }
                    std::cmp::Ordering::Equal
                })
            });
            if fs.len() == 1 && n == 2 {
                for r in roots_quadratic(&f) {
                    update(&mut result, &mut zeros, r, 1);
                }
            } else if fs.len() == 1 && fs[0].1 == 1 {
                for r in try_decompose(&dense::to_rat(&pz)) {
                    update(&mut result, &mut zeros, r, 1);
                }
            } else {
                for (g, kk) in fs {
                    let gq = upoly_from_rats(&dense::to_rat(&g));
                    for r in try_heuristics(&gq) {
                        update(&mut result, &mut zeros, r, kk);
                    }
                }
            }
        }
    } else if let Some((coeff, g)) = preprocess.then(|| preprocess_symbolic(&f)).flatten().filter(|(_, g)| g.c != f.c) {
        for (r, k) in roots_of_with(&g, false) {
            let r = if coeff.is_one() { r } else { mul(&[coeff.clone(), r]) };
            update(&mut result, &mut zeros, r, k);
        }
    } else {
        // symbolic coefficients
        let n = f.deg();
        let parts = if n == 1 {
            vec![(roots_linear(&f), 1)]
        } else if f.length() == 2 {
            if n == 2 {
                vec![(roots_quadratic(&f), 1)]
            } else {
                vec![(roots_binomial(&f), 1)]
            }
        } else {
            match symbolic_factors(&f) {
                // The reference factors `Poly(f.as_expr())` in every generator
                // and applies the heuristics to each factor in the variable.
                Some((count, factors, _)) if !(count == 1 && (n == 2 || factors.first().is_some_and(|g| g.1 == 1))) => factors
                    .iter()
                    .map(|(g, k)| (try_heuristics(g), *k))
                    .collect(),
                Some((1, _, _)) if n == 2 => vec![(roots_quadratic(&f), 1)],
                // Over EX the reference applies `to_rational_coeffs` (not
                // ported) and otherwise the heuristics, without decomposition.
                Some((1, _, true)) => vec![(try_heuristics(&f), 1)],
                Some((1, _, false)) => vec![(try_decompose_symbolic(&f), 1)],
                _ => vec![(try_heuristics(&f), 1)],
            }
        };
        for (rs, k) in parts {
            for r in rs {
                update(&mut result, &mut zeros, r, k);
            }
        }
    }
    if zeros > 0 && !result.iter().any(|x| x.0.is_zero()) {
        result.push((Expr::zero(), zeros));
    } else if let Some(x) = result.iter_mut().find(|x| x.0.is_zero()) {
        x.1 = zeros;
    }
    result
}

/// `roots(f, *gens)`: `[(root, multiplicity)]` in SymPy's dict order.
pub fn roots(e: &Expr, gens: &[Expr]) -> PResult<Vec<(Expr, u32)>> {
    // `I` is a ZZ_I coefficient, not a generator, by default.
    let selected = if gens.is_empty() {
        super::gaussian::default_selection(std::slice::from_ref(e), gens).map(|s| s.gens)
    } else {
        None
    };
    let gens = selected.as_deref().unwrap_or(gens);
    let p = match Poly::from_expr(e, gens) {
        Ok(p) => p,
        Err(PolyError::Failed(_)) => return Ok(vec![]),
        Err(x) => return Err(x),
    };
    roots_poly(&p)
}

pub fn roots_poly(p: &Poly) -> PResult<Vec<(Expr, u32)>> {
    if p.ngens != 1 {
        return Err(PolyError::Unsupported("multivariate polynomials are not supported".into()));
    }
    if p.is_float() && !p.has_domain_gens() {
        let rs = nroots_poly(p, 15)?;
        let mut out: Vec<(Expr, u32)> = Vec::new();
        for r in rs {
            if let Some(x) = out.iter_mut().find(|x| x.0 == r) {
                x.1 += 1;
            } else {
                out.push((r, 1));
            }
        }
        return Ok(out);
    }
    let c = p.all_coeffs()?;
    let numerical = !p.has_domain_gens();
    Ok(roots_of(&UPoly { c, numerical }))
}

// ----------------------------------------------------------------------
// numerical roots
// ----------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct DD {
    hi: f64,
    lo: f64,
}

impl DD {
    fn new(x: f64) -> DD {
        DD { hi: x, lo: 0.0 }
    }
    fn two_sum(a: f64, b: f64) -> DD {
        let s = a + b;
        let bb = s - a;
        let e = (a - (s - bb)) + (b - bb);
        DD { hi: s, lo: e }
    }
    fn add(self, o: DD) -> DD {
        let s = DD::two_sum(self.hi, o.hi);
        let t = DD::two_sum(self.lo, o.lo);
        let mut r = DD::two_sum(s.hi, s.lo + t.hi);
        r = DD::two_sum(r.hi, r.lo + t.lo);
        r
    }
    fn neg(self) -> DD {
        DD { hi: -self.hi, lo: -self.lo }
    }
    fn sub(self, o: DD) -> DD {
        self.add(o.neg())
    }
    fn mul(self, o: DD) -> DD {
        let p = self.hi * o.hi;
        let e = self.hi.mul_add(o.hi, -p);
        let e = e + (self.hi * o.lo + self.lo * o.hi);
        DD::two_sum(p, e)
    }
    fn div(self, o: DD) -> DD {
        let q1 = self.hi / o.hi;
        let r = self.sub(o.mul(DD::new(q1)));
        let q2 = r.hi / o.hi;
        let r = r.sub(o.mul(DD::new(q2)));
        let q3 = r.hi / o.hi;
        DD::two_sum(q1, q2).add(DD::new(q3))
    }
    fn from_rat(r: &Rat) -> DD {
        // n/d with both converted in double-double
        let n = DD::from_int(&r.n);
        let d = DD::from_int(&r.d);
        n.div(d)
    }
    fn from_int(i: &Int) -> DD {
        let hi = i.to_f64();
        if !hi.is_finite() {
            return DD::new(hi);
        }
        let hb = <num_bigint::BigInt as num_traits::FromPrimitive>::from_f64(hi).unwrap();
        let rest = i.to_big() - hb;
        DD { hi, lo: num_traits::ToPrimitive::to_f64(&rest).unwrap_or(0.0) }
    }
}

#[derive(Clone, Copy, Debug)]
struct CDD {
    re: DD,
    im: DD,
}

impl CDD {
    fn add(self, o: CDD) -> CDD {
        CDD { re: self.re.add(o.re), im: self.im.add(o.im) }
    }
    fn mul(self, o: CDD) -> CDD {
        CDD {
            re: self.re.mul(o.re).sub(self.im.mul(o.im)),
            im: self.re.mul(o.im).add(self.im.mul(o.re)),
        }
    }
    fn div(self, o: CDD) -> CDD {
        let den = o.re.mul(o.re).add(o.im.mul(o.im));
        let re = self.re.mul(o.re).add(self.im.mul(o.im)).div(den);
        let im = self.im.mul(o.re).sub(self.re.mul(o.im)).div(den);
        CDD { re, im }
    }
}

/// All complex roots of a polynomial with rational coefficients (low ->
/// high), sorted like SymPy's `nroots`.
pub fn nroots_rat(c: &[Rat]) -> Vec<C64> {
    let mut c = c.to_vec();
    dense::trim(&mut c);
    let mut zeros = 0;
    while c.len() > 1 && c[0].is_zero() {
        c.remove(0);
        zeros += 1;
    }
    let n = c.len().saturating_sub(1);
    let mut out: Vec<C64> = vec![C64::new(0.0, 0.0); zeros];
    if n == 0 {
        return sort_roots(out);
    }
    // normalise for the f64 iteration
    let lc = c[n].to_f64();
    let cf: Vec<f64> = c.iter().map(|x| x.to_f64() / lc).collect();
    let cd: Vec<CDD> = c
        .iter()
        .map(|x| CDD { re: DD::from_rat(&x.div(&c[n])), im: DD::new(0.0) })
        .collect();
    let mut z = aberth(&cf);
    // Newton polishing in double-double
    for zi in z.iter_mut() {
        let mut x = CDD { re: DD::new(zi.re), im: DD::new(zi.im) };
        for _ in 0..6 {
            let (p, dp) = horner_dd(&cd, x);
            if dp.re.hi == 0.0 && dp.im.hi == 0.0 {
                break;
            }
            let step = p.div(dp);
            x = CDD { re: x.re.sub(step.re), im: x.im.sub(step.im) };
            if step.re.hi.abs() + step.im.hi.abs() <= 1e-34 * (x.re.hi.abs() + x.im.hi.abs()) {
                break;
            }
        }
        *zi = C64::new(x.re.hi + x.re.lo, x.im.hi + x.im.lo);
    }
    // clean up: real roots of real polynomials
    let scale = z.iter().map(|w| w.norm()).fold(0.0, f64::max).max(1.0);
    for w in z.iter_mut() {
        if w.im.abs() <= 1e-13 * w.norm().max(1e-300) || w.im.abs() < 1e-16 * scale && w.im.abs() < 1e-12 * w.re.abs() {
            // check by refining as a real root
            let mut x = DD::new(w.re);
            let cr: Vec<DD> = cd.iter().map(|q| q.re).collect();
            for _ in 0..8 {
                let (p, dp) = horner_real_dd(&cr, x);
                if dp.hi == 0.0 {
                    break;
                }
                let st = p.div(dp);
                x = x.sub(st);
                if st.hi.abs() <= 1e-34 * x.hi.abs() {
                    break;
                }
            }
            *w = C64::new(x.hi + x.lo, 0.0);
        } else if w.re.abs() <= 1e-15 * w.im.abs() {
            w.re = 0.0;
        }
    }
    // conjugate symmetry
    let m = z.len();
    for i in 0..m {
        if z[i].im > 0.0 {
            let mut best = None;
            let mut bd = f64::INFINITY;
            for j in 0..m {
                if z[j].im < 0.0 {
                    let d = (z[j] - z[i].conj()).norm();
                    if d < bd {
                        bd = d;
                        best = Some(j);
                    }
                }
            }
            if let Some(j) = best {
                if bd <= 1e-10 * z[i].norm().max(1.0) {
                    let avg_re = 0.5 * (z[i].re + z[j].re);
                    let avg_im = 0.5 * (z[i].im - z[j].im);
                    z[i] = C64::new(avg_re, avg_im);
                    z[j] = C64::new(avg_re, -avg_im);
                }
            }
        }
    }
    out.extend(z);
    sort_roots(out)
}

fn sort_roots(mut v: Vec<C64>) -> Vec<C64> {
    v.sort_by(|a, b| {
        let ka = (a.im != 0.0) as u8;
        let kb = (b.im != 0.0) as u8;
        ka.cmp(&kb)
            .then(a.re.partial_cmp(&b.re).unwrap_or(std::cmp::Ordering::Equal))
            .then(a.im.abs().partial_cmp(&b.im.abs()).unwrap_or(std::cmp::Ordering::Equal))
            .then(a.im.signum().partial_cmp(&b.im.signum()).unwrap_or(std::cmp::Ordering::Equal))
    });
    v
}

fn horner_dd(c: &[CDD], x: CDD) -> (CDD, CDD) {
    let zero = CDD { re: DD::new(0.0), im: DD::new(0.0) };
    let mut p = zero;
    let mut dp = zero;
    for k in (0..c.len()).rev() {
        dp = dp.mul(x).add(p);
        p = p.mul(x).add(c[k]);
    }
    (p, dp)
}

fn horner_real_dd(c: &[DD], x: DD) -> (DD, DD) {
    let mut p = DD::new(0.0);
    let mut dp = DD::new(0.0);
    for k in (0..c.len()).rev() {
        dp = dp.mul(x).add(p);
        p = p.mul(x).add(c[k]);
    }
    (p, dp)
}

/// Aberth–Ehrlich simultaneous iteration (monic coefficients, low -> high).
fn aberth(c: &[f64]) -> Vec<C64> {
    aberth_steps(c, 500)
}

pub(super) fn aberth_steps(c: &[f64], maxsteps: usize) -> Vec<C64> {
    let n = c.len() - 1;
    if n == 1 {
        return vec![C64::new(-c[0], 0.0)];
    }
    // initial points on a circle of radius from the Cauchy bound
    let r = c[..n].iter().map(|x| x.abs()).fold(0.0, f64::max).powf(1.0 / n as f64).max(1e-3);
    let r = (1.0 + r).min(2.0 * r.max(1.0));
    let mut z: Vec<C64> = (0..n)
        .map(|k| {
            let t = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64) + 0.4;
            C64::from_polar(r * 0.9, t)
        })
        .collect();
    let eval = |x: C64| -> (C64, C64) {
        let mut p = C64::new(0.0, 0.0);
        let mut dp = C64::new(0.0, 0.0);
        for k in (0..c.len()).rev() {
            dp = dp * x + p;
            p = p * x + c[k];
        }
        (p, dp)
    };
    for _ in 0..maxsteps {
        let mut maxstep: f64 = 0.0;
        for i in 0..n {
            let (p, dp) = eval(z[i]);
            if p.norm() == 0.0 {
                continue;
            }
            let ratio = p / dp;
            let mut s = C64::new(0.0, 0.0);
            for j in 0..n {
                if j != i {
                    let d = z[i] - z[j];
                    if d.norm() > 0.0 {
                        s += 1.0 / d;
                    }
                }
            }
            let w = ratio / (C64::new(1.0, 0.0) - ratio * s);
            if w.re.is_finite() && w.im.is_finite() {
                z[i] -= w;
                maxstep = maxstep.max(w.norm() / z[i].norm().max(1e-300));
            }
        }
        if maxstep < 1e-15 {
            break;
        }
    }
    z
}

/// `nroots(f)` as expressions.
pub fn nroots(e: &Expr, gens: &[Expr], n: u32) -> PResult<Vec<Expr>> {
    nroots_options(e, gens, n, 50, true)
}

pub fn nroots_options(e: &Expr, gens: &[Expr], n: u32, maxsteps: u32, cleanup: bool) -> PResult<Vec<Expr>> {
    // `I` is a ZZ_I coefficient, not a generator, by default.
    // With floats (the reference CC) the generators are those without `I`.
    let selected = if gens.is_empty() && crate::traverse::has(e, &Expr::i()) {
        super::convert::from_exprs(std::slice::from_ref(e)).ok().and_then(|c| {
            let gens: Vec<Expr> = c.gens.iter().filter(|g| **g != Expr::i()).cloned().collect();
            (gens.len() + 1 == c.gens.len() && !gens.is_empty()).then_some(gens)
        })
    } else {
        None
    };
    let gens = selected.as_deref().unwrap_or(gens);
    let p = Poly::from_expr(e, gens)?;
    nroots_poly_options(&p, n, maxsteps, cleanup)
}

/// `nroots` over ZZ_I/QQ_I: the exact Gaussian square-free factors are
/// solved numerically with complex coefficients.
fn gaussian_nroots(p: &Poly, n: u32, maxsteps: u32, cleanup: bool) -> PResult<Vec<Expr>> {
    use super::algebraic::{AlgebraicPoly, AlgebraicScope, AN};
    use super::mpoly::MPoly;
    let field = super::gaussian::field();
    let _scope = AlgebraicScope::new(field.clone());
    // The exact coefficients in x of the representation in (x, I); float
    // inputs are their exact binary values, like the reference evaluation.
    let mut terms: std::collections::BTreeMap<u32, (Rat, Rat)> = std::collections::BTreeMap::new();
    for (e, c) in p.rep.iter() {
        let (x, y) = match e[1] % 4 {
            0 => (c.clone(), Rat::zero()),
            1 => (Rat::zero(), c.clone()),
            2 => (c.neg(), Rat::zero()),
            _ => (Rat::zero(), c.neg()),
        };
        let entry = terms.entry(e[0]).or_insert_with(|| (Rat::zero(), Rat::zero()));
        *entry = (entry.0.add(&x), entry.1.add(&y));
    }
    let terms: Vec<(Vec<u32>, AN)> = terms
        .into_iter()
        .map(|(k, (x, y))| (vec![k], super::gaussian::make(x, y)))
        .filter(|(_, c)| !super::int::Coeff::is_zero(c))
        .collect();
    let a = AlgebraicPoly {
        rep: MPoly::from_terms(1, terms),
        gens: p.gens[..1].to_vec(),
        number_field: field,
        gaussian: Some(false),
    };
    let (_, factors) = a.sqf_list();
    let factors: Vec<(Vec<(Rat, Rat)>, u32)> = factors
        .into_iter()
        .map(|(f, k)| {
            let degree = f.rep.degree(0) as usize;
            let mut c = vec![(Rat::zero(), Rat::zero()); degree + 1];
            for (e, coefficient) in f.rep.iter() {
                c[e[0] as usize] = super::gaussian::parts(coefficient);
            }
            (c, k)
        })
        .collect();
    super::nroots_mp::gaussian_roots(&factors, n, maxsteps, cleanup)
}

/// `Poly.nroots` over EX. The multiplicities come from the exact square-free
/// decomposition over EX; each square-free factor is then solved from
/// Gaussian rational approximations of its coefficients. The reference
/// evaluates the coefficients to `n` digits only; here the working precision
/// is raised until two successive evaluations give the same rounded roots.
/// `domain` names the coefficient domain in the reference error message.
pub fn expression_nroots(p: &super::expression_domain::ExpressionPoly, domain: &str, n: u32, maxsteps: u32, cleanup: bool) -> PResult<Vec<Expr>> {
    use crate::bigfloat::dps_to_prec;
    if p.gens.len() != 1 {
        return Err(PolyError::Unsupported("multivariate polynomial".into()));
    }
    if p.degree(0) <= 0 {
        return Ok(vec![]);
    }
    let numeric = || PolyError::Reference { kind: "DomainError", message: format!("Numerical domain expected, got {domain}") };
    for c in p.all_coeffs() {
        if !crate::traverse::free_symbols(&c).is_empty() {
            return Err(numeric());
        }
        if crate::evalf_mp::eval_mp(&c, 64).is_none() {
            return Err(PolyError::Unsupported(format!("cannot evaluate the coefficient {} numerically", crate::printer::to_str(&c))));
        }
    }
    let (_, factors) = p.sqf_list(false)?;
    let rat = |b: &crate::bigfloat::BigFloat| {
        let q = b.to_ratio_exact();
        Rat::new(super::int::Int::from_big(q.numer().clone()), super::int::Int::from_big(q.denom().clone()))
    };
    let approximate = |wp: u32| -> PResult<Vec<Expr>> {
        let mut parts = Vec::with_capacity(factors.len());
        for (f, k) in &factors {
            let mut c = Vec::new();
            for e in f.all_coeffs().iter().rev() {
                let z = crate::evalf_mp::eval_mp(e, wp).ok_or_else(numeric)?;
                c.push((rat(&z.re), rat(&z.im)));
            }
            parts.push((c, *k as u32));
        }
        super::nroots_mp::gaussian_roots(&parts, n, maxsteps, cleanup)
    };
    let mut wp = dps_to_prec(n) + 64;
    let mut previous = approximate(wp)?;
    for _ in 0..4 {
        wp = 2 * wp + 64;
        let next = approximate(wp)?;
        if next == previous {
            return Ok(next);
        }
        previous = next;
    }
    Ok(previous)
}

pub fn nroots_poly(p: &Poly, n: u32) -> PResult<Vec<Expr>> {
    nroots_poly_options(p, n, 50, true)
}

pub fn nroots_poly_options(p: &Poly, n: u32, maxsteps: u32, cleanup: bool) -> PResult<Vec<Expr>> {
    if p.ngens == 1 && p.gens[1..] == [Expr::i()] && p.modulus.is_none() {
        return gaussian_nroots(p, n, maxsteps, cleanup);
    }
    if p.ngens != 1 || p.has_domain_gens() {
        return Err(PolyError::Unsupported("can't compute numerical roots of a multivariate or symbolic polynomial".into()));
    }
    let c = dense::from_mpoly(&p.rep, 0);
    super::nroots_mp::rational_roots(&c, n, maxsteps, cleanup)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    use crate::printer::to_str;

    #[test]
    fn radical_controls_preserve_complete_solve_fallback_and_restore_context(){
        let disabled=RadicalOptions{cubics:false,quartics:false,incomplete:false};
        let x=Expr::symbol("x");let cubic=parse_str("x^3-x+1").unwrap();
        with_radical_options(disabled,||{
            assert!(roots(&cubic,&[x.clone()]).unwrap().is_empty());
            let complete=crate::solve::solve(&cubic,&x).unwrap();
            assert_eq!(complete.len(),3);
            assert!(complete.iter().all(|root|root.class_name()=="ComplexRootOf"));
            assert_eq!(roots(&parse_str("x^4+1").unwrap(),&[x.clone()]).unwrap().len(),4);
            with_radical_options(RadicalOptions::default(),||assert_eq!(roots(&cubic,&[x.clone()]).unwrap().len(),3));
            assert_eq!(radical_options(),disabled);
            assert_eq!(std::thread::spawn(radical_options).join().unwrap(),RadicalOptions::default());
        });
        assert_eq!(radical_options(),RadicalOptions::default());
        let panic=std::panic::catch_unwind(||with_radical_options(disabled,||panic!("scope test")));
        assert!(panic.is_err());
        assert_eq!(radical_options(),RadicalOptions::default());
    }

    fn symbolic(s: &str) -> Vec<String> {
        let mut v: Vec<String> = roots(&parse_str(s).unwrap(), &[Expr::symbol("x")]).unwrap().iter().map(|(e, k)| format!("{}: {k}", to_str(e))).collect();
        v.sort();
        v
    }

    #[test]
    fn symbolic_coefficients_follow_the_reference_preprocessing() {
        // preprocess_roots rescales x -> a*x before the binomial solver.
        assert_eq!(symbolic("a*x^4-16*a^5"), ["-2*I*a: 1", "-2*a: 1", "2*I*a: 1", "2*a: 1"]);
        // factor_list over ZZ(a) splits the quadratic into linear factors.
        assert_eq!(symbolic("x^2+(a+1/a)*x+1"), ["-1/a: 1", "-a: 1"]);
        // A single irreducible quartic is decomposed over ZZ(a).
        assert_eq!(symbolic("x^4-a^2*x^2-2*a^2-1"), [
            "-sqrt(a**2/2 + sqrt(a**4 + 8*a**2 + 4)/2): 1",
            "-sqrt(a**2/2 - sqrt(a**4 + 8*a**2 + 4)/2): 1",
            "sqrt(a**2/2 + sqrt(a**4 + 8*a**2 + 4)/2): 1",
            "sqrt(a**2/2 - sqrt(a**4 + 8*a**2 + 4)/2): 1",
        ]);
    }

    fn r(s: &str) -> String {
        let rs = roots(&parse_str(s).unwrap(), &[]).unwrap();
        let v: Vec<String> = rs.iter().map(|(e, k)| format!("{}: {k}", to_str(e))).collect();
        format!("{{{}}}", v.join(", "))
    }

    #[test]
    fn exact() {
        assert_eq!(r("x**2 - 2"), "{-sqrt(2): 1, sqrt(2): 1}");
        assert_eq!(r("(x - 1)**2*(x + 2)"), "{-2: 1, 1: 2}");
        assert_eq!(r("x**2 - x - 1"), "{1/2 - sqrt(5)/2: 1, 1/2 + sqrt(5)/2: 1}");
        assert_eq!(r("x**2 + x + 1"), "{-1/2 - sqrt(3)*I/2: 1, -1/2 + sqrt(3)*I/2: 1}");
        assert_eq!(r("x**3 - 1"), "{1: 1, -1/2 - sqrt(3)*I/2: 1, -1/2 + sqrt(3)*I/2: 1}");
        assert_eq!(r("x**4 - 5*x**2 + 4"), "{2: 1, 1: 1, -1: 1, -2: 1}");
    }

    #[test]
    fn numeric() {
        let v = nroots_rat(&[Rat::from_i64(-2), Rat::zero(), Rat::zero(), Rat::one()]);
        assert!((v[0].re - 2f64.cbrt()).abs() < 1e-15);
        assert_eq!(v[0].im, 0.0);
        assert!(v[1].im < 0.0 && v[2].im > 0.0);
    }
}
