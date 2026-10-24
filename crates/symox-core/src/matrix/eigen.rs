//! Eigenvalues and eigenvectors.
//!
//! Root finding is self-contained: rational roots (rational root theorem),
//! quadratic formula and small binomials, reproducing the order and form of
//! SymPy's `roots()` for these cases. Characteristic polynomials whose
//! irreducible factors have degree > 2 are reported as not implemented.

use super::int::{Int, Rat};
use super::linalg::RatMat;
use super::ratfunc::{cancel, is_zero_assume};
use super::{MResult, Matrix, MatrixError};
use crate::expand::expand;
use crate::expr::{Expr, Kind};
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};

// ----------------------------------------------------------------------
// univariate polynomials over Q (coefficients highest degree first)
// ----------------------------------------------------------------------

type UPoly = Vec<Rat>;

fn trim(p: &mut UPoly) {
    while p.len() > 1 && p[0].is_zero() {
        p.remove(0);
    }
}

fn deg(p: &UPoly) -> usize {
    p.len() - 1
}

fn eval(p: &UPoly, x: &Rat) -> Rat {
    let mut acc = Rat::zero();
    for c in p {
        acc = acc.mul(x).add(c);
    }
    acc
}

/// Divide by the monic linear factor `(x - r)` (exact).
fn deflate(p: &UPoly, r: &Rat) -> UPoly {
    let mut out = Vec::with_capacity(p.len() - 1);
    let mut acc = Rat::zero();
    for c in &p[..p.len() - 1] {
        acc = acc.mul(r).add(c);
        out.push(acc.clone());
    }
    out
}

fn poly_rem(a: &UPoly, b: &UPoly) -> UPoly {
    let mut r = a.clone();
    trim(&mut r);
    let db = deg(b);
    while r.len() > db && !(r.len() == 1 && r[0].is_zero()) {
        if r.len() - 1 < db {
            break;
        }
        let q = r[0].div(&b[0]);
        for (i, c) in b.iter().enumerate() {
            r[i] = r[i].sub(&q.mul(c));
        }
        r.remove(0);
        if r.is_empty() {
            r.push(Rat::zero());
            break;
        }
        trim(&mut r);
        if r.len() - 1 < db {
            break;
        }
    }
    if r.is_empty() {
        r.push(Rat::zero());
    }
    r
}

fn poly_quo(a: &UPoly, b: &UPoly) -> UPoly {
    let mut r = a.clone();
    let db = deg(b);
    let mut q = Vec::new();
    while r.len() > db {
        let c = r[0].div(&b[0]);
        for (i, x) in b.iter().enumerate() {
            r[i] = r[i].sub(&c.mul(x));
        }
        r.remove(0);
        q.push(c);
    }
    if q.is_empty() {
        q.push(Rat::zero());
    }
    q
}

fn is_zero_poly(p: &UPoly) -> bool {
    p.iter().all(|c| c.is_zero())
}

fn monic(p: &UPoly) -> UPoly {
    let l = p[0].clone();
    p.iter().map(|c| c.div(&l)).collect()
}

fn poly_gcd(a: &UPoly, b: &UPoly) -> UPoly {
    let (mut a, mut b) = (a.clone(), b.clone());
    trim(&mut a);
    trim(&mut b);
    while !is_zero_poly(&b) {
        let r = poly_rem(&a, &b);
        a = b;
        b = r;
        trim(&mut b);
    }
    monic(&a)
}

fn derivative(p: &UPoly) -> UPoly {
    let n = deg(p);
    if n == 0 {
        return vec![Rat::zero()];
    }
    (0..n).map(|i| p[i].mul(&Rat::from_int(Int::from((n - i) as i64)))).collect()
}

/// Yun's square-free decomposition: `[(factor, multiplicity)]` (monic).
fn squarefree(p: &UPoly) -> Vec<(UPoly, usize)> {
    let f = monic(p);
    if deg(&f) == 0 {
        return vec![];
    }
    let mut out = Vec::new();
    let fp = derivative(&f);
    let mut a = poly_gcd(&f, &fp);
    let mut b = poly_quo(&f, &a);
    let mut c = poly_quo(&fp, &a);
    let mut d: UPoly = {
        let bp = derivative(&b);
        let mut t: UPoly = c.clone();
        // d = c - b'
        let off = t.len().max(bp.len());
        let mut x = vec![Rat::zero(); off];
        for (i, v) in t.iter().rev().enumerate() {
            x[off - 1 - i] = x[off - 1 - i].add(v);
        }
        for (i, v) in bp.iter().rev().enumerate() {
            x[off - 1 - i] = x[off - 1 - i].sub(v);
        }
        trim(&mut x);
        t = x;
        t
    };
    let mut i = 1;
    loop {
        if deg(&b) == 0 {
            break;
        }
        a = poly_gcd(&b, &d);
        if deg(&a) > 0 {
            out.push((a.clone(), i));
        }
        b = poly_quo(&b, &a);
        c = poly_quo(&d, &a);
        let bp = derivative(&b);
        let off = c.len().max(bp.len());
        let mut x = vec![Rat::zero(); off];
        for (k, v) in c.iter().rev().enumerate() {
            x[off - 1 - k] = x[off - 1 - k].add(v);
        }
        for (k, v) in bp.iter().rev().enumerate() {
            x[off - 1 - k] = x[off - 1 - k].sub(v);
        }
        trim(&mut x);
        d = x;
        i += 1;
        if i > 64 {
            break;
        }
    }
    out
}

/// Primitive integer polynomial with positive leading coefficient.
fn primitive_int(p: &UPoly) -> Vec<Int> {
    let mut l = Int::ONE;
    for c in p {
        l = l.lcm(&c.den);
    }
    let mut v: Vec<Int> = p.iter().map(|c| c.num.mul(&l.divexact(&c.den))).collect();
    let mut g = Int::ZERO;
    for c in &v {
        g = g.gcd(c);
    }
    if !g.is_zero() && !g.is_one() {
        v = v.iter().map(|c| c.divexact(&g)).collect();
    }
    if v[0].is_negative() {
        v = v.iter().map(|c| c.neg()).collect();
    }
    v
}

fn divisors(n: &Int) -> Option<Vec<Int>> {
    let n = n.abs();
    let v = n.to_i64()?;
    if v > 1_000_000_000_000 {
        return None;
    }
    let v = v as u64;
    let mut small = Vec::new();
    let mut large = Vec::new();
    let mut d = 1u64;
    while d * d <= v {
        if v % d == 0 {
            small.push(d);
            if d * d != v {
                large.push(v / d);
            }
        }
        d += 1;
    }
    large.reverse();
    small.extend(large);
    Some(small.into_iter().map(|x| Int::from(x as i64)).collect())
}

/// Rational roots of a square-free polynomial.
fn rational_roots(p: &UPoly) -> MResult<Vec<Rat>> {
    let ip = primitive_int(p);
    let n = ip.len() - 1;
    if n == 0 {
        return Ok(vec![]);
    }
    let mut roots = Vec::new();
    let mut a0 = ip[n].clone();
    let mut k = n;
    while a0.is_zero() && k > 0 {
        roots.push(Rat::zero());
        k -= 1;
        a0 = ip[k].clone();
    }
    let too_big = || MatrixError::NotImplemented("polynomial coefficients too large for root finding".into());
    let (Some(ps), Some(qs)) = (divisors(&a0), divisors(&ip[0])) else {
        // Divisor enumeration is only a shortcut: SymPy factors the
        // characteristic polynomial whatever its coefficient size.
        return rational_roots_by_factoring(&ip, roots).ok_or_else(too_big);
    };
    let pq: UPoly = ip.iter().map(|c| Rat::from_int(c.clone())).collect();
    for q in &qs {
        for pp in &ps {
            for s in [pp.clone(), pp.neg()] {
                let r = Rat::new(s, q.clone());
                if roots.contains(&r) {
                    continue;
                }
                if eval(&pq, &r).is_zero() {
                    roots.push(r);
                }
            }
        }
    }
    Ok(roots)
}

/// Rational roots of the primitive integer polynomial `ip` (leading
/// coefficient first) from its linear factors over Q, appended to `roots`.
fn rational_roots_by_factoring(ip: &[Int], mut roots: Vec<Rat>) -> Option<Vec<Rat>> {
    let x = Expr::symbol("_lambda");
    let n = ip.len() - 1;
    let e = add(&ip.iter().enumerate().map(|(i, c)| mul2(&int_expr(c), &x.powi((n - i) as i64))).collect::<Vec<_>>());
    let (_, factors) = crate::poly::api::factor_list(&e).ok()?;
    for (f, _) in factors {
        let c = crate::integrate::poly_coeffs(&f, &x)?;
        if c.len() != 2 {
            continue;
        }
        let b = c[0].as_number().and_then(Rat::from_number)?;
        let a = c[1].as_number().and_then(Rat::from_number)?;
        let r = b.neg().div(&a);
        if !roots.contains(&r) {
            roots.push(r);
        }
    }
    Some(roots)
}

/// A factor of the characteristic polynomial over Q.
#[derive(Clone, Debug)]
struct Factor {
    /// primitive integer coefficients, positive leading coefficient
    coeffs: Vec<Int>,
    mult: usize,
}

/// SymPy's `factor_list` (over Q) restricted to factors of degree <= 2.
fn factor_list(p: &UPoly) -> MResult<Vec<Factor>> {
    let mut out: Vec<Factor> = Vec::new();
    for (sf, m) in squarefree(p) {
        let mut rest = sf.clone();
        for r in rational_roots(&sf)? {
            rest = deflate(&rest, &r);
            let lin = vec![r.den.clone(), r.num.neg()];
            out.push(Factor { coeffs: lin, mult: m });
        }
        match deg(&rest) {
            0 => {}
            2 | 3 => out.push(Factor { coeffs: primitive_int(&rest), mult: m }),
            _ => {
                // Split the rational-root-free part over Z (e.g. into two
                // quadratics); factor_roots handles what stays irreducible.
                let x = Expr::symbol("_lambda");
                let ip = primitive_int(&rest);
                let n = ip.len() - 1;
                let e = add(&ip.iter().enumerate().map(|(i, c)| mul2(&int_expr(c), &x.powi((n - i) as i64))).collect::<Vec<_>>());
                let unsupported = || MatrixError::NotImplemented(
                    "eigenvalues of characteristic polynomials with irreducible factors of degree > 4 are not implemented".into());
                let (_, factors) = crate::poly::api::factor_list(&e).map_err(|_| unsupported())?;
                for (f, k) in factors {
                    let c = crate::integrate::poly_coeffs(&f, &x).ok_or_else(unsupported)?;
                    let mut coeffs: Vec<Int> = c.iter().rev()
                        .map(|c| c.as_number().and_then(Rat::from_number).map(|r| r.num))
                        .collect::<Option<_>>().ok_or_else(unsupported)?;
                    if coeffs[0].is_negative() {
                        coeffs = coeffs.iter().map(|c| c.neg()).collect();
                    }
                    if coeffs.len() > 5 {
                        return Err(unsupported());
                    }
                    out.push(Factor { coeffs, mult: m * k.as_i64().ok_or_else(unsupported)? as usize });
                }
            }
        }
    }
    out.sort_by(|a, b| {
        a.coeffs
            .len()
            .cmp(&b.coeffs.len())
            .then(a.mult.cmp(&b.mult))
            .then_with(|| a.coeffs.cmp(&b.coeffs))
    });
    Ok(out)
}

fn rat_expr(r: &Rat) -> Expr {
    Expr::num(r.to_number())
}

fn int_expr(i: &Int) -> Expr {
    Expr::num(i.to_number())
}

/// SymPy's `roots_quadratic` for numeric coefficients.
fn roots_quadratic(a: &Expr, b: &Expr, c: &Expr) -> Vec<Expr> {
    let ia = pow(a, &Expr::minus_one());
    if c.is_zero() {
        let r1 = neg(&mul2(b, &ia));
        let (mut r0, mut r1) = (Expr::zero(), r1);
        if r1.as_number().is_some_and(|n| n.is_negative()) {
            std::mem::swap(&mut r0, &mut r1);
        }
        return vec![r0, r1];
    }
    if b.is_zero() {
        let r = neg(&mul2(c, &ia));
        let rr = sqrt(&r);
        return vec![neg(&rr), rr];
    }
    let d = sub(&pow(b, &Expr::int(2)), &mul(&[Expr::int(4), a.clone(), c.clone()]));
    let aa = mul2(&Expr::int(2), a);
    let bb = neg(&mul2(b, &pow(&aa, &Expr::minus_one())));
    let dd = mul2(&sqrt(&d), &pow(&aa, &Expr::minus_one()));
    let (mut r0, mut r1) = (sub(&bb, &dd), add(&[bb, dd]));
    if a.as_number().is_some_and(|n| n.is_negative()) {
        std::mem::swap(&mut r0, &mut r1);
    }
    vec![r0, r1]
}

/// `exp(2*pi*I*k/n)` expanded, for small `n`.
fn root_of_unity(k: i64, n: i64) -> Option<Expr> {
    let k = k.rem_euclid(n);
    let half = Expr::rational(1, 2);
    let s3 = mul2(&sqrt(&Expr::int(3)), &half);
    let i = Expr::i();
    Some(match (n, k) {
        (_, 0) => Expr::one(),
        (2, 1) => Expr::minus_one(),
        (3, 1) => add(&[neg(&half), mul2(&s3, &i)]),
        (3, 2) => add(&[neg(&half), neg(&mul2(&s3, &i))]),
        (4, 1) => i,
        (4, 2) => Expr::minus_one(),
        (4, 3) => neg(&i),
        _ => return None,
    })
}

/// SymPy's `roots_binomial` (`a*x**n + b`) for small `n`.
pub(crate) fn roots_binomial(a: &Expr, b: &Expr, n: usize) -> MResult<Vec<Expr>> {
    let base = neg(&cancel(&mul2(b, &pow(a, &Expr::minus_one()))));
    let ni = n as i64;
    let neg_base = base.as_number().is_some_and(|x| x.is_negative());
    let alpha = if n == 2 {
        sqrt(&base)
    } else if !neg_base {
        pow(&base, &Expr::rational(1, ni))
    } else {
        return Err(MatrixError::NotImplemented("roots of binomials with negative base".into()));
    };
    let even = n % 2 == 0;
    let mut ks: Vec<i64> = Vec::new();
    let mut imax = ni / 2;
    if even {
        ks.push(imax);
        imax -= 1;
    }
    if !neg_base {
        ks.push(0);
    }
    for i in (1..=imax).rev() {
        if neg_base {
            ks.extend([i, -i]);
        } else {
            ks.extend([-i, i]);
        }
    }
    if neg_base {
        ks.push(0);
    }
    let mut out = Vec::new();
    for k in ks {
        let z = root_of_unity(k, ni).ok_or_else(|| {
            MatrixError::NotImplemented("roots of high degree binomials".into())
        })?;
        out.push(expand(&mul2(&alpha, &z)));
    }
    Ok(out)
}

/// Roots of an irreducible-or-linear factor (SymPy's `_try_heuristics`).
fn factor_roots(f: &[Int]) -> MResult<Vec<Expr>> {
    let n = f.len() - 1;
    let nonzero = f.iter().filter(|c| !c.is_zero()).count();
    if nonzero == 2 {
        if n == 1 {
            return Ok(vec![Expr::num(Rat::new(f[1].neg(), f[0].clone()).to_number())]);
        }
        return roots_binomial(&int_expr(&f[0]), &int_expr(&f[n]), n);
    }
    if n == 2 {
        return Ok(roots_quadratic(&int_expr(&f[0]), &int_expr(&f[1]), &int_expr(&f[2]))
            .iter()
            .map(cancel)
            .collect());
    }
    if n == 3 {
        return roots_cubic(f);
    }
    // Quartics (Ferrari/Euler) and decomposable higher degrees go through the
    // polynomial `roots`, which is what SymPy's eigenvals calls as well.
    let x = Expr::symbol("_lambda");
    let poly = add(&f.iter().enumerate().map(|(i, c)| mul2(&int_expr(c), &x.powi((n - i) as i64))).collect::<Vec<_>>());
    if let Ok(rs) = crate::poly::roots::roots(&poly, &[x]) {
        let found: Vec<Expr> = rs.into_iter().flat_map(|(r, k)| std::iter::repeat_n(r, k as usize)).collect();
        if found.len() == n {
            return Ok(found);
        }
    }
    Err(MatrixError::NotImplemented("roots of factors of degree > 4".into()))
}

/// SymPy's `roots_cubic` (Cardano) for rational coefficients.
fn roots_cubic(f: &[Int]) -> MResult<Vec<Expr>> {
    let lc = int_expr(&f[0]);
    let il = pow(&lc, &Expr::minus_one());
    let (a, b, c) = (mul2(&int_expr(&f[1]), &il), mul2(&int_expr(&f[2]), &il), mul2(&int_expr(&f[3]), &il));
    roots_cubic_expr(&a, &b, &c, &|p: &[Expr]| {
        let rs: Option<Vec<Rat>> = p.iter().map(|e| e.as_number().and_then(Rat::from_number)).collect();
        Ok(multiple(&roots_rational(&rs.unwrap())?))
    })
}

/// SymPy's `roots_cubic` for the monic cubic `x**3 + a*x**2 + b*x + c`
/// (`quad` computes `roots([1, u, v], multiple=True)`).
pub(crate) fn roots_cubic_expr(
    a: &Expr,
    b: &Expr,
    c: &Expr,
    quad: &dyn Fn(&[Expr]) -> MResult<Vec<Expr>>,
) -> MResult<Vec<Expr>> {
    let third = Expr::rational(1, 3);
    let three = Expr::int(3);
    if c.is_zero() {
        let rs = quad(&[Expr::one(), a.clone(), b.clone()])?;
        return Ok(vec![rs[0].clone(), Expr::zero(), rs[1].clone()]);
    }
    let a2 = pow(a, &Expr::int(2));
    let p = sub(b, &div(&a2, &three));
    let q = add(&[
        sub(c, &div(&mul2(a, b), &three)),
        div(&mul2(&Expr::int(2), &pow(a, &three)), &Expr::int(27)),
    ]);
    let pon3 = div(&p, &three);
    let aon3 = div(a, &three);
    let positive = |e: &Expr| e.as_number().is_some_and(|n| n.is_positive());
    let negative = |e: &Expr| e.as_number().is_some_and(|n| n.is_negative());
    let mut u1: Option<Expr> = None;
    if p.is_zero() {
        if q.is_zero() {
            let x = neg(&aon3);
            return Ok(vec![x.clone(), x.clone(), x]);
        }
        u1 = Some(if positive(&q) { neg(&pow(&q, &third)) } else { pow(&neg(&q), &third) });
    } else if q.is_zero() {
        let ys = quad(&[Expr::one(), Expr::zero(), p.clone()])?;
        let ma = neg(&aon3);
        return Ok(vec![add(&[ys[0].clone(), ma.clone()]), ma.clone(), add(&[ys[1].clone(), ma])]);
    } else if negative(&q) {
        let inner = add(&[div(&pow(&q, &Expr::int(2)), &Expr::int(4)), pow(&pon3, &three)]);
        let base = add(&[neg(&div(&q, &Expr::int(2))), sqrt(&inner)]);
        u1 = Some(neg(&pow(&base, &third)));
    }
    let coeff = mul(&[Expr::i(), sqrt(&three), Expr::rational(1, 2)]);
    let mhalf = Expr::rational(-1, 2);
    let Some(u1) = u1 else {
        let u2 = add(&[mhalf.clone(), coeff.clone()]);
        let u3 = sub(&mhalf, &coeff);
        let d0 = sub(&a2, &mul2(&three, b));
        // Python's left to right evaluation: (9*b)*c distributes 9 first
        let d1 = add(&[
            sub(&mul2(&Expr::int(2), &pow(a, &three)), &mul2(&mul2(&Expr::int(9), a), b)),
            mul2(&Expr::int(27), c),
        ]);
        let disc = sub(&pow(&d1, &Expr::int(2)), &mul2(&Expr::int(4), &pow(&d0, &three)));
        let cc = pow(&div(&add(&[d1.clone(), sqrt(&disc)]), &Expr::int(2)), &third);
        let mut out = Vec::new();
        for uk in [Expr::one(), u2, u3] {
            // -(b + uk*C + D0/C/uk)/3
            let t = div(&div(&d0, &cc), &uk);
            let s = add(&[a.clone(), mul2(&uk, &cc), t]);
            out.push(div(&neg(&s), &three));
        }
        return Ok(out);
    };
    let u2 = mul2(&u1, &add(&[mhalf.clone(), coeff.clone()]));
    let u3 = mul2(&u1, &sub(&mhalf, &coeff));
    let ma = neg(&aon3);
    if p.is_zero() {
        return Ok(vec![add(&[u1, ma.clone()]), add(&[u2, ma.clone()]), add(&[u3, ma])]);
    }
    Ok([u1, u2, u3]
        .iter()
        .map(|u| add(&[neg(u), div(&pon3, u), ma.clone()]))
        .collect())
}

/// `roots(..., multiple=True)` from the multiset form.
pub(crate) fn multiple(r: &[(Expr, usize)]) -> Vec<Expr> {
    let mut v = Vec::new();
    for (e, k) in r {
        for _ in 0..*k {
            v.push(e.clone());
        }
    }
    v
}

/// SymPy's `roots(p)` (as an insertion ordered multiset) for a polynomial
/// with rational coefficients given highest degree first.
pub fn roots_rational(coeffs: &[Rat]) -> MResult<Vec<(Expr, usize)>> {
    let mut p: UPoly = coeffs.to_vec();
    trim(&mut p);
    let mut result: Vec<(Expr, usize)> = Vec::new();
    let push = |r: Expr, k: usize, result: &mut Vec<(Expr, usize)>| {
        if let Some(x) = result.iter_mut().find(|(e, _)| *e == r) {
            x.1 += k;
        } else {
            result.push((r, k));
        }
    };
    // terms_gcd
    let mut zeros = 0;
    while p.len() > 1 && p.last().unwrap().is_zero() {
        p.pop();
        zeros += 1;
    }
    if deg(&p) > 0 {
        let ip = primitive_int(&p);
        let n = ip.len() - 1;
        let nonzero = ip.iter().filter(|c| !c.is_zero()).count();
        if n == 1 {
            push(Expr::num(Rat::new(ip[1].neg(), ip[0].clone()).to_number()), 1, &mut result);
        } else if nonzero == 2 {
            let rs = if n == 2 {
                roots_quadratic(&int_expr(&ip[0]), &int_expr(&ip[1]), &int_expr(&ip[2]))
            } else {
                roots_binomial(&int_expr(&ip[0]), &int_expr(&ip[n]), n)?
            };
            for r in rs {
                push(r, 1, &mut result);
            }
        } else {
            let fl = factor_list(&p)?;
            if fl.len() == 1 && n == 2 {
                for r in roots_quadratic(&int_expr(&ip[0]), &int_expr(&ip[1]), &int_expr(&ip[2])) {
                    push(r, 1, &mut result);
                }
            } else if fl.len() == 1 && fl[0].mult == 1 {
                // irreducible of degree > 2 (_try_decompose)
                for r in factor_roots(&fl[0].coeffs)? {
                    push(r, 1, &mut result);
                }
            } else {
                for f in &fl {
                    for r in factor_roots(&f.coeffs)? {
                        push(r, f.mult, &mut result);
                    }
                }
            }
        }
    }
    if zeros > 0 {
        push(Expr::zero(), zeros, &mut result);
    }
    Ok(result)
}

pub use crate::sortkey::cmp_default;

// ----------------------------------------------------------------------
// algebraic number fields Q[t]/(minpoly)
// ----------------------------------------------------------------------

/// `Q[t]/(t^n + m_{n-1} t^{n-1} + ... + m_0)`; elements are coefficient
/// vectors (low degree first) of length `n`.
struct AlgF {
    m: Vec<Rat>,
}

type AE = Vec<Rat>;

fn ptrim(p: &mut Vec<Rat>) {
    while p.len() > 1 && p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
}

/// Polynomial division (low degree first).
fn pdivmod(a: &[Rat], b: &[Rat]) -> (Vec<Rat>, Vec<Rat>) {
    let mut r = a.to_vec();
    ptrim(&mut r);
    let mut b = b.to_vec();
    ptrim(&mut b);
    let db = b.len() - 1;
    if r.len() < b.len() {
        return (vec![Rat::zero()], r);
    }
    let mut q = vec![Rat::zero(); r.len() - db];
    while r.len() > db && !(r.len() == 1 && r[0].is_zero()) {
        let k = r.len() - 1 - db;
        let c = r.last().unwrap().div(&b[db]);
        for (i, bc) in b.iter().enumerate() {
            r[i + k] = r[i + k].sub(&c.mul(bc));
        }
        q[k] = c;
        r.pop();
        ptrim(&mut r);
        if r.len() <= db && r.len() < b.len() {
            break;
        }
    }
    if r.is_empty() {
        r.push(Rat::zero());
    }
    (q, r)
}

fn pmul(a: &[Rat], b: &[Rat]) -> Vec<Rat> {
    let mut out = vec![Rat::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] = out[i + j].add(&x.mul(y));
        }
    }
    out
}

fn psub(a: &[Rat], b: &[Rat]) -> Vec<Rat> {
    let n = a.len().max(b.len());
    let mut out = vec![Rat::zero(); n];
    for (i, x) in a.iter().enumerate() {
        out[i] = x.clone();
    }
    for (i, y) in b.iter().enumerate() {
        out[i] = out[i].sub(y);
    }
    ptrim(&mut out);
    out
}

impl AlgF {
    fn n(&self) -> usize {
        self.m.len()
    }
    fn zero(&self) -> AE {
        vec![Rat::zero(); self.n()]
    }
    fn one(&self) -> AE {
        let mut v = self.zero();
        v[0] = Rat::one();
        v
    }
    fn modulus(&self) -> Vec<Rat> {
        let mut v = self.m.clone();
        v.push(Rat::one());
        v
    }
    fn reduce(&self, mut p: Vec<Rat>) -> AE {
        let n = self.n();
        while p.len() > n {
            let t = p.pop().unwrap();
            let k = p.len() - n;
            if !t.is_zero() {
                for (j, mj) in self.m.iter().enumerate() {
                    p[k + j] = p[k + j].sub(&t.mul(mj));
                }
            }
        }
        p.resize(n, Rat::zero());
        p
    }
    fn sub(&self, x: &AE, y: &AE) -> AE {
        x.iter().zip(y).map(|(a, b)| a.sub(b)).collect()
    }
    fn mul(&self, x: &AE, y: &AE) -> AE {
        self.reduce(pmul(x, y))
    }
    fn is_zero(x: &AE) -> bool {
        x.iter().all(|c| c.is_zero())
    }
    fn inv(&self, x: &AE) -> AE {
        // extended Euclid: s*x = r (mod modulus)
        let (mut r0, mut r1) = (self.modulus(), {
            let mut v = x.clone();
            ptrim(&mut v);
            v
        });
        let (mut s0, mut s1) = (vec![Rat::zero()], vec![Rat::one()]);
        while !(r1.len() == 1 && r1[0].is_zero()) {
            let (q, r) = pdivmod(&r0, &r1);
            r0 = std::mem::replace(&mut r1, r);
            let ns = psub(&s0, &pmul(&q, &s1));
            s0 = std::mem::replace(&mut s1, ns);
        }
        let c = r0[0].inv();
        self.reduce(s0.iter().map(|v| v.mul(&c)).collect())
    }
}

/// Nullspace over the algebraic field (standard basis from the rref).
fn nullspace_alg(f: &AlgF, rows: usize, cols: usize, mut a: Vec<AE>) -> Vec<Vec<AE>> {
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p) = (r..rows).find(|&i| !AlgF::is_zero(&a[i * cols + c])) else {
            continue;
        };
        for k in 0..cols {
            a.swap(p * cols + k, r * cols + k);
        }
        let inv = f.inv(&a[r * cols + c]);
        for k in 0..cols {
            a[r * cols + k] = f.mul(&a[r * cols + k], &inv);
        }
        for i in 0..rows {
            if i == r {
                continue;
            }
            let x = a[i * cols + c].clone();
            if AlgF::is_zero(&x) {
                continue;
            }
            for k in 0..cols {
                let v = f.sub(&a[i * cols + k], &f.mul(&x, &a[r * cols + k]));
                a[i * cols + k] = v;
            }
        }
        pivots.push(c);
        r += 1;
    }
    let mut basis = Vec::new();
    for free in (0..cols).filter(|c| !pivots.contains(c)) {
        let mut v = vec![f.zero(); cols];
        v[free] = f.one();
        for (pr, &pc) in pivots.iter().enumerate() {
            v[pc] = f.sub(&v[pc], &a[pr * cols + free]);
        }
        basis.push(v);
    }
    basis
}

// ----------------------------------------------------------------------
// public API
// ----------------------------------------------------------------------

/// `(eigenvalue, algebraic multiplicity, eigenvectors)`.
pub type EigenVect = (Expr, usize, Vec<Matrix>);

impl Matrix {
    /// Eigenvalues with multiplicities, in SymPy's dict order.
    pub fn eigenvals(&self) -> MResult<Vec<(Expr, usize)>> {
        if self.is_empty() {
            return Ok(vec![]);
        }
        if !self.is_square() {
            return Err(MatrixError::NonSquare(format!("{} must be a square matrix.", self.to_str())));
        }
        let is_rat = RatMat::from_matrix(self).is_some();
        if !is_rat && self.data().iter().any(|e| contains_float(e)) {
            return Err(MatrixError::NotImplemented(
                "numerical eigenvalues of float matrices are not implemented".into(),
            ));
        }
        let mut all: Vec<(Expr, usize)> = Vec::new();
        let push = |v: Expr, k: usize, all: &mut Vec<(Expr, usize)>| {
            if let Some(x) = all.iter_mut().find(|(e, _)| *e == v) {
                x.1 += k;
            } else {
                all.push((v, k));
            }
        };
        for b in super::linalg::scc(self) {
            if is_rat && b.len() == 1 {
                push(self.get(b[0], b[0]).clone(), 1, &mut all);
                continue;
            }
            let blk = self.extract(&b, &b)?;
            let eigs = block_roots(&blk)?;
            let total: usize = eigs.iter().map(|(_, k)| k).sum();
            if total != blk.rows() {
                return Err(MatrixError::Other(
                    "It is not always possible to express the eigenvalues of a matrix of size 5x5 or higher in radicals.".into(),
                ));
            }
            for (v, k) in eigs {
                push(v, k, &mut all);
            }
        }
        Ok(all)
    }

    /// Eigenvalues, eigenvector bases, sorted like SymPy.
    pub fn eigenvects(&self) -> MResult<Vec<EigenVect>> {
        if !self.is_square() {
            return Err(MatrixError::NonSquare(format!("{} must be a square matrix.", self.to_str())));
        }
        if let Some(r) = RatMat::from_matrix(self) {
            return eigenvects_rational(self, &r);
        }
        match super::domain::matrix_domain(self) {
            super::domain::Domain::Poly { ctx, .. } | super::domain::Domain::Frac { ctx, .. } => {
                if let Some(r) = super::symeigen::eigenvects_field(self, &ctx) {
                    return r;
                }
            }
            _ => {}
        }
        let mut vals = self.eigenvals()?;
        vals.sort_by(|a, b| cmp_default(&a.0, &b.0));
        let n = self.rows();
        let mut out = Vec::new();
        for (val, mult) in vals {
            let m = self.sub(&Matrix::eye(n).scale(&val))?;
            let vs = m.nullspace();
            if vs.is_empty() {
                return Err(MatrixError::NotImplemented(format!(
                    "Can't evaluate eigenvector for eigenvalue {val}"
                )));
            }
            out.push((val, mult, vs));
        }
        Ok(out)
    }
}

fn contains_float(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_float(),
        _ => e.args().iter().any(contains_float),
    }
}

/// Roots of the characteristic polynomial of a block.
fn block_roots(blk: &Matrix) -> MResult<Vec<(Expr, usize)>> {
    let cp = blk.charpoly("lambda")?;
    // rational coefficients?
    let rats: Option<Vec<Rat>> = cp
        .coeffs
        .iter()
        .map(|c| c.as_number().and_then(Rat::from_number))
        .collect();
    if let Some(r) = rats {
        return roots_rational(&r);
    }
    symbolic_roots(&cp.coeffs)
}

fn symbolic_roots(c: &[Expr]) -> MResult<Vec<(Expr, usize)>> {
    super::symeigen::roots_symbolic(c)
}

/// Rational roots of a polynomial with rational coefficients.
pub fn rational_roots_pub(p: &[Rat]) -> MResult<Vec<Rat>> {
    let mut v: UPoly = p.to_vec();
    trim(&mut v);
    let mut out = Vec::new();
    for (sf, _) in squarefree(&v) {
        out.extend(rational_roots(&sf)?);
    }
    Ok(out)
}

fn eigenvects_rational(m: &Matrix, r: &RatMat) -> MResult<Vec<EigenVect>> {
    let n = m.rows();
    let cp = r.charpoly();
    let fl = factor_list(&cp)?;
    let mut out: Vec<EigenVect> = Vec::new();
    for f in &fl {
        if f.coeffs.len() == 2 {
            let ev = Rat::new(f.coeffs[1].neg(), f.coeffs[0].clone());
            let e = rat_expr(&ev);
            let mut sh = r.clone();
            for i in 0..n {
                sh.data[i * n + i] = sh.data[i * n + i].sub(&ev);
            }
            let vs = sh.to_matrix().nullspace();
            out.push((e, f.mult, vs));
        } else {
            // algebraic eigenvalue: minimal polynomial = the factor
            let d = f.coeffs.len() - 1;
            let a0 = Rat::from_int(f.coeffs[0].clone());
            let mono: Vec<Rat> = (0..d).map(|k| Rat::from_int(f.coeffs[d - k].clone()).div(&a0)).collect();
            let fld = AlgF { m: mono };
            let mut theta = fld.zero();
            theta[1] = Rat::one();
            let mut a: Vec<AE> = r
                .data
                .iter()
                .map(|x| {
                    let mut v = fld.zero();
                    v[0] = x.clone();
                    v
                })
                .collect();
            for i in 0..n {
                a[i * n + i] = fld.sub(&a[i * n + i], &theta);
            }
            let basis = nullspace_alg(&fld, n, n, a);
            // roots of the minimal polynomial, SymPy's roots() order
            let mp: Vec<Rat> = f.coeffs.iter().map(|c| Rat::from_int(c.clone())).collect();
            let rts = roots_rational(&mp)?;
            if rts.iter().map(|(_, k)| k).sum::<usize>() != d {
                return Err(MatrixError::NotImplemented("eigenvectors for CRootOf eigenvalues".into()));
            }
            for (root, _) in rts {
                let vs: Vec<Matrix> = basis
                    .iter()
                    .map(|v| {
                        Matrix::column(
                            v.iter()
                                .map(|x| {
                                    let terms: Vec<Expr> = x
                                        .iter()
                                        .enumerate()
                                        .filter(|(_, c)| !c.is_zero())
                                        .map(|(k, c)| {
                                            if k == 0 {
                                                rat_expr(c)
                                            } else {
                                                mul2(&rat_expr(c), &pow(&root, &Expr::int(k as i64)))
                                            }
                                        })
                                        .collect();
                                    add(&terms)
                                })
                                .collect(),
                        )
                    })
                    .collect();
                out.push((root, f.mult, vs));
            }
        }
    }
    out.sort_by(|a, b| cmp_default(&a.0, &b.0));
    let _ = is_zero_assume;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(v: &[i64]) -> Vec<Rat> {
        v.iter().map(|&x| Rat::from_int(Int::from(x))).collect()
    }

    fn show(v: &[(Expr, usize)]) -> String {
        v.iter().map(|(e, k)| format!("{e}: {k}")).collect::<Vec<_>>().join(", ")
    }

    #[test]
    fn roots_basic() {
        assert_eq!(show(&roots_rational(&r(&[1, -5, -2])).unwrap()), "5/2 - sqrt(33)/2: 1, 5/2 + sqrt(33)/2: 1");
        assert_eq!(show(&roots_rational(&r(&[1, -6, 11, -6])).unwrap()), "3: 1, 2: 1, 1: 1");
        assert_eq!(show(&roots_rational(&r(&[1, 0, -2])).unwrap()), "-sqrt(2): 1, sqrt(2): 1");
        assert_eq!(show(&roots_rational(&r(&[1, -2, 1])).unwrap()), "1: 2");
        assert_eq!(show(&roots_rational(&r(&[1, 0, 1, 0])).unwrap()), "-I: 1, I: 1, 0: 1");
        assert_eq!(show(&roots_rational(&r(&[1, -3, 0, 4])).unwrap()), "-1: 1, 2: 2");
    }
}
