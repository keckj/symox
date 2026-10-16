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
use crate::ops::{add, mul, mul2, neg, pow, sqrt, sub};

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
    let ps = divisors(&a0).ok_or_else(too_big)?;
    let qs = divisors(&ip[0]).ok_or_else(too_big)?;
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
            2 => out.push(Factor { coeffs: primitive_int(&rest), mult: m }),
            _ => {
                return Err(MatrixError::NotImplemented(
                    "eigenvalues of characteristic polynomials with irreducible factors of degree > 2 are not implemented".into(),
                ))
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
fn roots_binomial(a: &Expr, b: &Expr, n: usize) -> MResult<Vec<Expr>> {
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
    Err(MatrixError::NotImplemented("roots of factors of degree > 2".into()))
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
                // irreducible of degree > 2
                return Err(MatrixError::NotImplemented(
                    "roots of irreducible polynomials of degree > 2".into(),
                ));
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
// quadratic number fields Q[t]/(t^2 + p t + q)
// ----------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
struct QF {
    a: Rat,
    b: Rat,
}

#[derive(Clone, Debug)]
struct Field2 {
    p: Rat,
    q: Rat,
}

impl Field2 {
    fn zero(&self) -> QF {
        QF { a: Rat::zero(), b: Rat::zero() }
    }
    fn one(&self) -> QF {
        QF { a: Rat::one(), b: Rat::zero() }
    }
    fn sub(&self, x: &QF, y: &QF) -> QF {
        QF { a: x.a.sub(&y.a), b: x.b.sub(&y.b) }
    }
    fn mul(&self, x: &QF, y: &QF) -> QF {
        // (a + b t)(c + d t) = ac + (ad+bc) t + bd t^2,  t^2 = -p t - q
        let bd = x.b.mul(&y.b);
        QF {
            a: x.a.mul(&y.a).sub(&bd.mul(&self.q)),
            b: x.a.mul(&y.b).add(&x.b.mul(&y.a)).sub(&bd.mul(&self.p)),
        }
    }
    fn inv(&self, x: &QF) -> QF {
        // (a + b t)^-1 = (a - b p - b t) / (a^2 - a b p + b^2 q)
        let n = x.a.mul(&x.a).sub(&x.a.mul(&x.b).mul(&self.p)).add(&x.b.mul(&x.b).mul(&self.q));
        let ni = n.inv();
        QF { a: x.a.sub(&x.b.mul(&self.p)).mul(&ni), b: x.b.neg().mul(&ni) }
    }
    fn is_zero(&self, x: &QF) -> bool {
        x.a.is_zero() && x.b.is_zero()
    }
}

/// Nullspace of a matrix over the quadratic field (standard basis from
/// the reduced row echelon form).
fn nullspace_qf(f: &Field2, rows: usize, cols: usize, mut a: Vec<QF>) -> Vec<Vec<QF>> {
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p) = (r..rows).find(|&i| !f.is_zero(&a[i * cols + c])) else {
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
            if f.is_zero(&x) {
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
            // quadratic minimal polynomial t^2 + p t + q (monic)
            let a0 = Rat::from_int(f.coeffs[0].clone());
            let p = Rat::from_int(f.coeffs[1].clone()).div(&a0);
            let q = Rat::from_int(f.coeffs[2].clone()).div(&a0);
            let fld = Field2 { p: p.clone(), q: q.clone() };
            let theta = QF { a: Rat::zero(), b: Rat::one() };
            let mut a: Vec<QF> = r.data.iter().map(|x| QF { a: x.clone(), b: Rat::zero() }).collect();
            for i in 0..n {
                a[i * n + i] = fld.sub(&a[i * n + i], &theta);
            }
            let basis = nullspace_qf(&fld, n, n, a);
            // roots of the minimal polynomial, SymPy's roots() order
            let rts = roots_rational(&[Rat::one(), p, q])?;
            for (root, _) in rts {
                let vs: Vec<Matrix> = basis
                    .iter()
                    .map(|v| {
                        Matrix::column(
                            v.iter()
                                .map(|x| add(&[mul2(&rat_expr(&x.b), &root), rat_expr(&x.a)]))
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
