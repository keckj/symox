//! Parametric Risch differential equation (SymPy `prde.py`), base case.
//!
//! Every polynomial carries the ZZ/QQ domain that the reference `Poly`
//! methods would give it (`unify`, `quo`/`div` retraction, `gcd`/`lcm`
//! normalization, `cancel(include=True)`), so callers can publish the same
//! `Poly` objects. A kernel returns `None` whenever the reference would raise
//! or take a branch that is not ported; the Python body then runs unchanged.
use crate::poly::{dense, Int, Rat};

/// Degree of the zero polynomial (SymPy `-oo`).
pub const NEG_INF: i64 = -(1 << 40);
/// Largest admitted degree bound `n` (loops are linear in it).
const MAX_N: i64 = 4096;

/// Dense polynomial in the active variable with its SymPy domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct P {
    /// Coefficients, lowest degree first, no trailing zeros.
    pub c: Vec<Rat>,
    /// Domain QQ (otherwise ZZ, all coefficients integral).
    pub qq: bool,
}

fn integral(c: &[Rat]) -> bool {
    c.iter().all(|x| x.is_integer())
}

fn trimmed(mut c: Vec<Rat>) -> Vec<Rat> {
    dense::trim(&mut c);
    c
}

/// Exact division `a / b` over QQ (b != 0).
fn divrem(a: &[Rat], b: &[Rat]) -> (Vec<Rat>, Vec<Rat>) {
    dense::divrem(a, b)
}

fn monic(a: &[Rat]) -> Vec<Rat> {
    if a.is_empty() {
        return Vec::new();
    }
    let l = dense::lc(a).inv();
    a.iter().map(|x| x.mul(&l)).collect()
}

fn gcd_monic(a: &[Rat], b: &[Rat]) -> Vec<Rat> {
    if a.is_empty() {
        return monic(b);
    }
    if b.is_empty() {
        return monic(a);
    }
    monic(&dense::gcd_field(a, b))
}

/// Integer content (positive) of an integral coefficient vector.
fn content(c: &[Rat]) -> Int {
    let mut g = Int::zero();
    for x in c {
        g = g.gcd(&x.n);
    }
    g.abs()
}

/// Least common denominator of the coefficients.
fn denominator(c: &[Rat]) -> Int {
    let mut l = Int::one();
    for x in c {
        l = l.lcm(&x.d);
    }
    l.abs()
}

fn scale(c: &[Rat], k: &Rat) -> Vec<Rat> {
    if k.is_zero() {
        return Vec::new();
    }
    c.iter().map(|x| x.mul(k)).collect()
}

fn int_rat(i: &Int) -> Rat {
    Rat::from(i.clone())
}

/// Primitive integral polynomial with a positive leading coefficient that is
/// an associate of `c` over QQ.
fn primitive_positive(c: &[Rat]) -> Vec<Rat> {
    if c.is_empty() {
        return Vec::new();
    }
    let cleared = scale(c, &int_rat(&denominator(c)));
    let k = content(&cleared);
    let mut out: Vec<Rat> = cleared.iter().map(|x| Rat::new(x.n.div_exact(&k).unwrap_or_else(Int::zero), Int::one())).collect();
    if dense::lc(&out).is_negative() {
        out = out.iter().map(|x| x.neg()).collect();
    }
    out
}

/// GCD over ZZ (`dup_inner_gcd`): content gcd times primitive gcd, positive
/// leading coefficient; a single zero argument returns the other one with
/// its content kept.
fn zz_gcd(f: &[Rat], g: &[Rat]) -> Vec<Rat> {
    if f.is_empty() && g.is_empty() {
        return Vec::new();
    }
    let sign = |p: &[Rat]| -> Vec<Rat> {
        if dense::lc(p).is_negative() {
            p.iter().map(|x| x.neg()).collect()
        } else {
            p.to_vec()
        }
    };
    if f.is_empty() {
        return sign(g);
    }
    if g.is_empty() {
        return sign(f);
    }
    let c = content(f).gcd(&content(g)).abs();
    scale(&primitive_positive(&dense::gcd_field(f, g)), &int_rat(&c))
}

impl P {
    pub fn new(c: Vec<Rat>, qq: bool) -> P {
        P { c: trimmed(c), qq }
    }
    /// Domain inferred from the coefficients, as `expr.as_poly(t)`.
    pub fn inferred(c: Vec<Rat>) -> P {
        let c = trimmed(c);
        let qq = !integral(&c);
        P { c, qq }
    }
    pub fn zero(qq: bool) -> P {
        P { c: Vec::new(), qq }
    }
    pub fn one(qq: bool) -> P {
        P { c: vec![Rat::one()], qq }
    }
    /// `Poly(r*t**n, t)`.
    pub fn monomial(r: &Rat, n: usize) -> P {
        if r.is_zero() {
            return P::zero(false);
        }
        let mut c = vec![Rat::zero(); n + 1];
        c[n] = r.clone();
        P::inferred(c)
    }
    fn retract(c: Vec<Rat>, unified_qq: bool) -> P {
        let c = trimmed(c);
        let qq = unified_qq || !integral(&c);
        P { c, qq }
    }
    pub fn is_zero(&self) -> bool {
        self.c.is_empty()
    }
    pub fn is_ground(&self) -> bool {
        self.c.len() <= 1
    }
    pub fn deg(&self) -> i64 {
        if self.c.is_empty() {
            NEG_INF
        } else {
            self.c.len() as i64 - 1
        }
    }
    pub fn lc(&self) -> Rat {
        dense::lc(&self.c)
    }
    pub fn nth(&self, i: i64) -> Rat {
        if i < 0 {
            return Rat::zero();
        }
        self.c.get(i as usize).cloned().unwrap_or_else(Rat::zero)
    }
    pub fn add(&self, o: &P) -> P {
        P::new(dense::add(&self.c, &o.c), self.qq || o.qq)
    }
    pub fn sub(&self, o: &P) -> P {
        P::new(dense::sub(&self.c, &o.c), self.qq || o.qq)
    }
    pub fn mul(&self, o: &P) -> P {
        P::new(dense::mul(&self.c, &o.c), self.qq || o.qq)
    }
    pub fn neg(&self) -> P {
        P::new(dense::neg(&self.c), self.qq)
    }
    /// `Poly.diff(t)`: same domain.
    pub fn diff(&self) -> P {
        P::new(dense::derivative(&self.c), self.qq)
    }
    /// Base-case `derivation(p, DE)` with `D = [1]`: `(p').as_poly(t)`.
    pub fn derivation(&self) -> P {
        P::inferred(dense::derivative(&self.c))
    }
    pub fn to_field(&self) -> P {
        P { c: self.c.clone(), qq: true }
    }
    pub fn quo(&self, o: &P) -> Option<P> {
        if o.is_zero() {
            return None;
        }
        Some(P::retract(divrem(&self.c, &o.c).0, self.qq || o.qq))
    }
    pub fn div(&self, o: &P) -> Option<(P, P)> {
        if o.is_zero() {
            return None;
        }
        let (q, r) = divrem(&self.c, &o.c);
        let qq = self.qq || o.qq || !integral(&q) || !integral(&r);
        Some((P::new(q, qq), P::new(r, qq)))
    }
    /// `Poly.exquo`; `None` for an inexact quotient (ExactQuotientFailed).
    pub fn exquo(&self, o: &P) -> Option<P> {
        if o.is_zero() {
            return None;
        }
        let (q, r) = divrem(&self.c, &o.c);
        if !r.is_empty() {
            return None;
        }
        Some(P::retract(q, self.qq || o.qq))
    }
    pub fn exquo_ground(&self, k: &Rat) -> Option<P> {
        if k.is_zero() || !self.qq {
            return None;
        }
        let inv = k.inv();
        Some(P::new(scale(&self.c, &inv), true))
    }
    pub fn gcd(&self, o: &P) -> P {
        if self.qq || o.qq {
            P::new(gcd_monic(&self.c, &o.c), true)
        } else {
            P::new(zz_gcd(&self.c, &o.c), false)
        }
    }
    /// `None` when both are zero (the reference divides by a zero gcd).
    pub fn lcm(&self, o: &P) -> Option<P> {
        let qq = self.qq || o.qq;
        if self.is_zero() && o.is_zero() {
            return if qq { None } else { Some(P::zero(false)) };
        }
        if self.is_zero() || o.is_zero() {
            return Some(P::zero(qq));
        }
        if qq {
            let h = divrem(&dense::mul(&self.c, &o.c), &gcd_monic(&self.c, &o.c)).0;
            return Some(P::new(monic(&h), true));
        }
        let (fc, gc) = (content(&self.c), content(&o.c));
        let fp = scale(&self.c, &int_rat(&fc).inv());
        let gp = scale(&o.c, &int_rat(&gc).inv());
        let h = divrem(&dense::mul(&fp, &gp), &zz_gcd(&fp, &gp)).0;
        let mut c = int_rat(&fc.lcm(&gc).abs());
        if dense::lc(&h).is_negative() {
            c = c.neg();
        }
        Some(P::new(scale(&h, &c), false))
    }
    /// `self.cancel(o, include=True)` of the stock `Poly.cancel`: both
    /// results are in the domain inferred from the coefficients of the two
    /// inputs (`Poly.from_exprs` in `an_rational::poly_cancel`).
    pub fn cancel(&self, o: &P) -> Option<(P, P)> {
        if o.is_zero() {
            return None;
        }
        let qq = !integral(&self.c) || !integral(&o.c);
        let (mut cq, mut cp) = (Int::one(), Int::one());
        let (f, g) = if qq {
            cq = denominator(&self.c);
            cp = denominator(&o.c);
            (scale(&self.c, &int_rat(&cq)), scale(&o.c, &int_rat(&cp)))
        } else {
            (self.c.clone(), o.c.clone())
        };
        let (mut p, mut q) = if f.is_empty() {
            (Vec::new(), vec![if dense::lc(&g).is_negative() { Rat::one().neg() } else { Rat::one() }])
        } else {
            let h = zz_gcd(&f, &g);
            (divrem(&f, &h).0, divrem(&g, &h).0)
        };
        if qq {
            let k = cp.gcd(&cq).abs();
            cp = cp.div_exact(&k)?;
            cq = cq.div_exact(&k)?;
        }
        let (pn, qn) = (dense::lc(&p).is_negative(), dense::lc(&q).is_negative());
        if pn && qn {
            p = dense::neg(&p);
            q = dense::neg(&q);
        } else if pn {
            cp = cp.neg();
            p = dense::neg(&p);
        } else if qn {
            cp = cp.neg();
            q = dense::neg(&q);
        }
        Some((P::new(scale(&p, &int_rat(&cp)), qq), P::new(scale(&q, &int_rat(&cq)), qq)))
    }
    /// `Poly.half_gcdex` (over QQ): `(s, h)` with `s*self = h (mod o)`.
    pub fn half_gcdex(&self, o: &P) -> Option<(P, P)> {
        if self.is_zero() && o.is_zero() {
            return None;
        }
        let (mut f, mut g) = (self.c.clone(), o.c.clone());
        let (mut a, mut b) = (vec![Rat::one()], Vec::new());
        while !g.is_empty() {
            let (q, r) = divrem(&f, &g);
            f = std::mem::replace(&mut g, r);
            let next = dense::sub(&a, &dense::mul(&q, &b));
            a = std::mem::replace(&mut b, next);
        }
        let l = dense::lc(&f).inv();
        Some((P::new(scale(&a, &l), true), P::new(monic(&f), true)))
    }
}

/// `gcdex_diophantine(a, b, c)`.
pub fn gcdex_diophantine(a: &P, b: &P, c: &P) -> Option<(P, P)> {
    let (s, g) = a.half_gcdex(b)?;
    let mut s = s.mul(&c.exquo(&g)?);
    if !s.is_zero() && s.deg() >= b.deg() {
        s = s.div(b)?.1;
    }
    let t = c.sub(&s.mul(a)).exquo(b)?;
    Some((s, t))
}

/// Dense rational matrix (row major). Publication uses PolyMatrix QQ.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mat {
    pub rows: usize,
    pub cols: usize,
    pub e: Vec<Rat>,
}

impl Mat {
    pub fn zeros(rows: usize, cols: usize) -> Mat {
        Mat { rows, cols, e: vec![Rat::zero(); rows * cols] }
    }
    pub fn eye(n: usize) -> Mat {
        let mut m = Mat::zeros(n, n);
        for i in 0..n {
            m.e[i * n + i] = Rat::one();
        }
        m
    }
    pub fn at(&self, i: usize, j: usize) -> &Rat {
        &self.e[i * self.cols + j]
    }
    pub fn from_fn(rows: usize, cols: usize, f: impl Fn(usize, usize) -> Rat) -> Mat {
        let mut e = Vec::with_capacity(rows * cols);
        for i in 0..rows {
            for j in 0..cols {
                e.push(f(i, j));
            }
        }
        Mat { rows, cols, e }
    }
    pub fn row_join(&self, o: &Mat) -> Option<Mat> {
        if self.rows != o.rows {
            return None;
        }
        Some(Mat::from_fn(self.rows, self.cols + o.cols, |i, j| {
            if j < self.cols { self.at(i, j).clone() } else { o.at(i, j - self.cols).clone() }
        }))
    }
    pub fn col_join(&self, o: &Mat) -> Option<Mat> {
        if self.cols != o.cols {
            return None;
        }
        let mut e = self.e.clone();
        e.extend(o.e.iter().cloned());
        Some(Mat { rows: self.rows + o.rows, cols: self.cols, e })
    }
    pub fn neg(&self) -> Mat {
        Mat { rows: self.rows, cols: self.cols, e: self.e.iter().map(|x| x.neg()).collect() }
    }
    /// Reduced row echelon form (same shape) and pivot columns, with the
    /// first nonzero entry of each column as pivot (PolyMatrix order).
    pub fn rref(&self) -> (Mat, Vec<usize>) {
        let mut rows: Vec<Vec<Rat>> = (0..self.rows).map(|i| self.e[i * self.cols..(i + 1) * self.cols].to_vec()).collect();
        let mut pivots = Vec::new();
        let mut r = 0;
        for j in 0..self.cols {
            if r == self.rows {
                break;
            }
            let Some(p) = (r..self.rows).find(|&i| !rows[i][j].is_zero()) else { continue };
            rows.swap(r, p);
            let inv = rows[r][j].inv();
            rows[r] = rows[r].iter().map(|x| x.mul(&inv)).collect();
            for i in 0..self.rows {
                if i != r && !rows[i][j].is_zero() {
                    let c = rows[i][j].clone();
                    let pivot_row = rows[r].clone();
                    for (x, y) in rows[i].iter_mut().zip(pivot_row.iter()) {
                        *x = x.sub(&c.mul(y));
                    }
                }
            }
            pivots.push(j);
            r += 1;
        }
        (Mat { rows: self.rows, cols: self.cols, e: rows.into_iter().flatten().collect() }, pivots)
    }
    /// Canonical nullspace basis (one per free column, 1 at that column).
    pub fn nullspace(&self) -> Vec<Vec<Rat>> {
        let (r, pivots) = self.rref();
        let mut out = Vec::new();
        for j in 0..self.cols {
            if pivots.contains(&j) {
                continue;
            }
            let mut v = vec![Rat::zero(); self.cols];
            v[j] = Rat::one();
            for (i, &p) in pivots.iter().enumerate() {
                if !r.at(i, j).is_zero() {
                    v[p] = r.at(i, j).neg();
                }
            }
            out.push(v);
        }
        out
    }
}

/// `constant_system(A, u, DE)` with ground rational entries: the derivation
/// loop never applies because no entry involves a monomial.
pub fn constant_system(a: &Mat, u: &Mat) -> Option<(Mat, Mat)> {
    if a.rows * a.cols == 0 {
        return Some((a.clone(), u.clone()));
    }
    let (au, _) = a.row_join(u)?.rref();
    let c = au.cols;
    if c == 0 {
        return None;
    }
    let left = Mat::from_fn(au.rows, c - 1, |i, j| au.at(i, j).clone());
    let right = Mat::from_fn(au.rows, 1, |i, _| au.at(i, c - 1).clone());
    Some((left, right))
}

fn constant_system_a(m: &Mat) -> Option<Mat> {
    Some(constant_system(m, &Mat::zeros(m.rows, 1))?.0)
}

/// Base-case `splitfactor(p, DE)`: every square factor is normal.
pub fn splitfactor(p: &P) -> Option<(P, P)> {
    let one = P::one(p.qq);
    if p.is_zero() {
        return Some((p.clone(), one));
    }
    if p.is_ground() {
        let s = p.gcd(&p.derivation());
        let n = p.exquo(&s)?;
        return Some((n, s));
    }
    Some((p.clone(), one))
}

/// Base-case `weak_normalizer(a, d, DE)` when the simple part `d1` of the
/// normal denominator is a constant (the resultant has no `z`). Declines
/// otherwise.
pub fn weak_normalizer(a: &P, d: &P) -> Option<(P, P, P)> {
    let (dn, _) = splitfactor(d)?;
    let g = dn.gcd(&dn.diff());
    let sqf = dn.quo(&g)?;
    let d1 = sqf.quo(&sqf.gcd(&g))?;
    if !d1.is_ground() || d1.is_zero() {
        return None;
    }
    gcdex_diophantine(&d.quo(&d1)?, &d1, a)?;
    Some((P::one(false), a.clone(), d.clone()))
}

pub struct NormalDenom {
    pub a: P,
    pub ba: P,
    pub bd: P,
    pub g: Vec<(P, P)>,
    pub h: P,
}

pub fn prde_normal_denom(fa: &P, fd: &P, g: &[(P, P)]) -> Option<NormalDenom> {
    if g.is_empty() {
        return None;
    }
    let (dn, _) = splitfactor(fd)?;
    let mut gd = P::one(false);
    for (_, d) in g {
        gd = gd.lcm(d)?;
    }
    let (en, _) = splitfactor(&gd)?;
    let p = dn.gcd(&en);
    let h = en.gcd(&en.diff()).quo(&p.gcd(&p.diff()))?;
    let a = dn.mul(&h);
    let c = a.mul(&h);
    let ba = a.mul(fa).sub(&dn.mul(&h.derivation()).mul(fd));
    let (ba, bd) = ba.cancel(fd)?;
    let g = g.iter().map(|(x, d)| c.mul(x).cancel(d)).collect::<Option<Vec<_>>>()?;
    Some(NormalDenom { a, ba, bd, g, h })
}

/// Base/primitive `prde_special_denom`: `(a, ba.quo(bd), G, Poly(1, t))`.
pub fn prde_special_denom(ba: &P, bd: &P) -> Option<(P, P)> {
    Some((ba.quo(bd)?, P::one(false)))
}

pub fn prde_linear_constraints(g: &[(P, P)]) -> Option<(Vec<P>, Mat)> {
    let m = g.len();
    let mut d = g.first()?.1.clone();
    for (_, x) in &g[1..] {
        d = d.lcm(x)?;
    }
    let d = d.to_field();
    let q = g.iter().map(|(a, x)| a.mul(&d.quo(x)?).div(&d)).collect::<Option<Vec<_>>>()?;
    let m_ = if q.iter().all(|(_, r)| r.is_zero()) {
        Mat::zeros(0, m)
    } else {
        let n = q.iter().map(|(_, r)| r.deg()).max()?;
        Mat::from_fn((n + 1) as usize, m, |i, j| q[j].1.nth(i as i64))
    };
    Some((q.into_iter().map(|(a, _)| a).collect(), m_))
}

pub fn poly_linear_constraints(p: &[P], d: &P) -> Option<(Vec<P>, Mat)> {
    let m = p.len();
    if m == 0 {
        return None;
    }
    let qr = p.iter().map(|x| x.div(d)).collect::<Option<Vec<_>>>()?;
    let mat = if qr.iter().all(|(_, r)| r.is_zero()) {
        Mat::zeros(0, m)
    } else {
        let n = qr.iter().map(|(_, r)| r.deg()).max()?;
        Mat::from_fn((n + 1) as usize, m, |i, j| qr[j].1.nth(i as i64))
    };
    Some((qr.into_iter().map(|(q, _)| q).collect(), mat))
}

pub struct Spde {
    pub a: P,
    pub b: P,
    pub q: Vec<P>,
    pub r: Vec<P>,
    pub n: i64,
}

pub fn prde_spde(a: &P, b: &P, q: &[P], n: i64) -> Option<Spde> {
    let rz = q.iter().map(|qi| gcdex_diophantine(b, a, qi)).collect::<Option<Vec<_>>>()?;
    let bb = b.add(&a.derivation());
    let qq = rz.iter().map(|(r, z)| z.sub(&r.derivation())).collect();
    let r = rz.into_iter().map(|(r, _)| r).collect();
    Some(Spde { a: a.clone(), b: bb, q: qq, r, n: n - a.deg() })
}

fn checked_n(n: i64) -> Option<i64> {
    (n <= MAX_N).then_some(n)
}

fn max_degree(q: &[P]) -> i64 {
    q.iter().map(P::deg).max().unwrap_or(NEG_INF)
}

/// Shared tail of the b-large and b-small (deg b > 0) kernels.
fn no_cancel_tail(q: &[P], m: usize) -> Option<Mat> {
    let dc = if q.iter().all(P::is_zero) { -1 } else { max_degree(q) };
    let mm = Mat::from_fn((dc + 1) as usize, m, |i, j| q[j].nth(i as i64));
    let a = constant_system_a(&mm)?;
    let c = Mat::eye(m);
    a.row_join(&Mat::zeros(a.rows, m))?.col_join(&c.row_join(&c.neg())?)
}

/// `prde_no_cancel_b_large(b, Q, n, DE)`; `q` holds the final `Q` values.
pub fn prde_no_cancel_b_large(b: &P, q: &mut [P], n: i64) -> Option<(Vec<P>, Mat)> {
    let n = checked_n(n)?;
    let m = q.len();
    if b.is_zero() || m == 0 {
        return None;
    }
    let db = b.deg();
    let blc = b.lc();
    let mut h = vec![P::zero(false); m];
    for big_n in (0..=n).rev() {
        for i in 0..m {
            let si = q[i].nth(big_n + db).div(&blc);
            let sitn = P::monomial(&si, big_n as usize);
            h[i] = h[i].add(&sitn);
            q[i] = q[i].sub(&sitn.derivation()).sub(&b.mul(&sitn));
        }
    }
    Some((h, no_cancel_tail(q, m)?))
}

/// Base-case `prde_no_cancel_b_small(b, Q, n, DE)` (`DE.d == 1`).
pub fn prde_no_cancel_b_small(b: &P, q: &mut [P], n: i64) -> Option<(Vec<P>, Mat)> {
    let n = checked_n(n)?;
    let m = q.len();
    if m == 0 {
        return None;
    }
    let mut h = vec![P::zero(false); m];
    for big_n in (1..=n).rev() {
        for i in 0..m {
            let si = q[i].nth(big_n - 1).div(&Rat::from_i64(big_n));
            let sitn = P::monomial(&si, big_n as usize);
            h[i] = h[i].add(&sitn);
            q[i] = q[i].sub(&sitn.derivation()).sub(&b.mul(&sitn));
        }
    }
    if b.deg() > 0 {
        let (db, blc) = (b.deg(), b.lc());
        for i in 0..m {
            let si = P::monomial(&q[i].nth(db).div(&blc), 0);
            h[i] = h[i].add(&si);
            q[i] = q[i].sub(&si.derivation()).sub(&b.mul(&si));
        }
        return Some((h, no_cancel_tail(q, m)?));
    }
    let f = P::one(true);
    let mut brow: Vec<Rat> = q.iter().map(|x| x.nth(0)).collect();
    brow.push(Rat::zero());
    let bm = Mat { rows: 1, cols: m + 1, e: brow };
    let d = max_degree(q);
    let a = if d > 0 {
        constant_system_a(&Mat::from_fn(d as usize, m, |i, j| q[j].nth(i as i64 + 1)))?
    } else {
        Mat::zeros(0, m)
    };
    let r = 1;
    let a = a.row_join(&Mat::zeros(a.rows, r + m))?;
    let bm = bm.row_join(&Mat::zeros(1, m))?;
    let c = Mat::eye(m).row_join(&Mat::zeros(m, r))?.row_join(&Mat::eye(m).neg())?;
    let mut out = vec![f];
    out.extend(h);
    Some((out, a.col_join(&bm)?.col_join(&c)?))
}

/// `(Matrix([polys]) * v)[0]` for a QQ vector `v`: a QQ polynomial.
fn combine(polys: &[P], v: &[Rat]) -> P {
    let mut acc = Vec::new();
    for (p, c) in polys.iter().zip(v) {
        acc = dense::add(&acc, &scale(&p.c, c));
    }
    P::new(acc, true)
}

/// `-eye(m) | V | zeros(m, k)` stacked over `zeros(B.rows, m) | B`.
fn assemble(m: usize, v: &[Vec<Rat>], k: usize, b: &Mat) -> Option<Mat> {
    let mut a = Mat::eye(m).neg();
    for vj in v {
        a = a.row_join(&Mat { rows: m, cols: 1, e: vj.clone() })?;
    }
    a = a.row_join(&Mat::zeros(m, k))?;
    a.col_join(&Mat::zeros(b.rows, m).row_join(b)?)
}

/// `param_poly_rischDE(a, b, q, n, DE)` in the base case.
pub fn param_poly_rischde(a: &P, b: &P, q: &[P], n: i64, depth: usize) -> Option<(Vec<P>, Mat)> {
    let m = q.len();
    if m == 0 || a.is_zero() || depth > 64 {
        return None;
    }
    let n = checked_n(n)?;
    if n < 0 {
        if q.iter().all(P::is_zero) {
            return Some((Vec::new(), Mat::zeros(1, m)));
        }
        let big_n = max_degree(q);
        let mm = Mat::from_fn((big_n + 1) as usize, m, |i, j| q[j].nth(i as i64));
        return Some((Vec::new(), constant_system_a(&mm)?));
    }
    if a.is_ground() {
        let al = a.lc();
        let b = b.to_field().exquo_ground(&al)?;
        let mut q = q.iter().map(|x| x.to_field().exquo_ground(&al)).collect::<Option<Vec<_>>>()?;
        return if !b.is_zero() {
            prde_no_cancel_b_large(&b, &mut q, n)
        } else {
            prde_no_cancel_b_small(&b, &mut q, n)
        };
    }
    let (mut a, mut b, mut q, mut n) = (a.clone(), b.clone(), q.to_vec(), n);
    let mut alpha = P::one(a.qq);
    let mut beta = vec![P::zero(a.qq); m];
    let mut d;
    loop {
        let s = prde_spde(&a, &b, &q, n)?;
        (a, b, q, n) = (s.a, s.b, s.q, s.n);
        beta = beta.iter().zip(&s.r).map(|(bi, ri)| bi.add(&alpha.mul(ri))).collect();
        alpha = alpha.mul(&a);
        d = a.gcd(&b);
        if !d.is_ground() || n < 0 {
            break;
        }
    }
    let (qq, mm) = poly_linear_constraints(&q, &d)?;
    let am = constant_system_a(&mm)?;
    let v = am.nullspace();
    if v.is_empty() {
        return Some((Vec::new(), Mat::eye(m)));
    }
    let r: Vec<P> = v.iter().map(|vj| combine(&qq, vj)).collect();
    let f: Vec<P> = v.iter().map(|vj| combine(&beta, vj)).collect();
    let (g, bm) = param_poly_rischde(&a.quo(&d)?, &b.quo(&d)?, &r, n, depth + 1)?;
    let mut h = f;
    h.extend(g.iter().map(|gk| alpha.mul(gk)));
    Some((h, assemble(m, &v, g.len(), &bm)?))
}

/// Base-case `bound_degree(a, b, Q, DE, parametric=True)`.
pub fn bound_degree(a: &P, b: &P, q: &[P]) -> Option<i64> {
    if a.is_zero() || q.is_empty() {
        return None;
    }
    let (da, db, dc) = (a.deg(), b.deg(), max_degree(q));
    let alpha = b.lc().neg().div(&a.lc());
    let mut n = 0.max(dc - db.max(da - 1));
    if db == da - 1 && alpha.is_integer() {
        let al = alpha.n.as_i64()?;
        n = 0.max(al).max(dc - db);
    }
    Some(n)
}

pub enum ParamRisch {
    /// `([], eye(m))`.
    NoSolution(usize),
    Solution(Vec<(P, P)>, Mat),
}

/// `param_rischDE(fa, fd, G, DE)` in the base case.
pub fn param_rischde(fa: &P, fd: &P, g: &[(P, P)]) -> Option<ParamRisch> {
    let m = g.len();
    if m == 0 || fd.is_zero() {
        return None;
    }
    let (q, fa, fd) = weak_normalizer(fa, fd)?;
    let mut gamma = q.clone();
    let g = g.iter().map(|(ga, gd)| q.mul(ga).cancel(gd)).collect::<Option<Vec<_>>>()?;
    let nd = prde_normal_denom(&fa, &fd, &g)?;
    gamma = gamma.mul(&nd.h);
    let (bb, hs) = prde_special_denom(&nd.ba, &nd.bd)?;
    gamma = gamma.mul(&hs);
    let aa = nd.a;
    let gg = aa.gcd(&bb);
    let (a, b) = (aa.quo(&gg)?, bb.quo(&gg)?);
    let g = nd.g.iter().map(|(x, d)| x.cancel(&d.mul(&gg))).collect::<Option<Vec<_>>>()?;
    let (q, mm) = prde_linear_constraints(&g)?;
    let mm = constant_system_a(&mm)?;
    let v = mm.nullspace();
    if v.is_empty() {
        return Some(ParamRisch::NoSolution(m));
    }
    let r: Vec<P> = v.iter().map(|vj| combine(&q, vj)).collect();
    let n = bound_degree(&a, &b, &r)?;
    let (h, bm) = param_poly_rischde(&a, &b, &r, n, 0)?;
    let am = assemble(m, &v, h.len(), &bm)?;
    let w = am.nullspace();
    let k = h.len();
    if k == 0 {
        // No relation left (`W` empty): `M` is `0 x m` and `C` is `eye(m)`.
        if w.is_empty() {
            return Some(ParamRisch::NoSolution(m));
        }
        // Otherwise `wl[-0:]` keeps the whole vector; leave that to the reference.
        return None;
    }
    let mut items = Vec::with_capacity(w.len() * (m + k));
    for wl in &w {
        items.extend(wl[..m].iter().cloned());
        items.extend(wl[wl.len() - k..].iter().cloned());
    }
    let big_m = Mat { rows: w.len(), cols: m + k, e: items };
    let nn = big_m.nullspace();
    let c = Mat { rows: nn.len(), cols: if nn.is_empty() { 0 } else { m + k }, e: nn.into_iter().flatten().collect() };
    let h = h.iter().map(|hk| hk.cancel(&gamma)).collect::<Option<Vec<_>>>()?;
    Some(ParamRisch::Solution(h, c))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::int::Coeff;
    fn z(c: &[i64]) -> P {
        P::new(c.iter().map(|&x| Rat::from_i64(x)).collect(), false)
    }
    #[test]
    fn cancel_include_matches_reference() {
        let (p, q) = z(&[4, 2]).cancel(&z(&[0, -6])).unwrap();
        assert_eq!((p, q), (z(&[-2, -1]), z(&[0, 3])));
        let half = P::new(vec![Rat::one(), Rat::new(Int::one(), Int::from_i64(2))], true);
        let tq = P::new(vec![Rat::zero(), Rat::new(Int::from_i64(3), Int::from_i64(4))], true);
        let (p, q) = half.cancel(&tq).unwrap();
        assert_eq!(p, P::new(vec![Rat::from_i64(4), Rat::from_i64(2)], true));
        assert_eq!(q, P::new(vec![Rat::zero(), Rat::from_i64(3)], true));
        let one = P::new(vec![Rat::one()], true);
        assert_eq!(one.cancel(&z(&[1])).unwrap(), (z(&[1]), z(&[1])));
    }
    #[test]
    fn zz_gcd_and_lcm_normalize_like_reference() {
        assert_eq!(z(&[-6, -4]).gcd(&z(&[])), z(&[6, 4]));
        assert_eq!(z(&[-6, -4]).lcm(&z(&[9, 6])).unwrap(), z(&[18, 12]));
        assert_eq!(z(&[0, 4]).lcm(&z(&[])).unwrap(), z(&[]));
    }
    #[test]
    fn nullspace_is_canonical() {
        let m = Mat { rows: 1, cols: 3, e: vec![Rat::one(), Rat::from_i64(2), Rat::zero()] };
        let v = m.nullspace();
        assert_eq!(v, vec![vec![Rat::from_i64(-2), Rat::one(), Rat::zero()], vec![Rat::zero(), Rat::zero(), Rat::one()]]);
    }
}
