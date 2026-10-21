//! Differential operators in `z*d/dz` (SymPy's `Operator` family) as
//! polynomials in the formal variable `_x`, and the machinery that applies
//! them to coefficient vectors.

use crate::ops::{add, div, mul, neg, sub};
use crate::Expr;

use super::{is_integer_number, is_nonpositive_integer, HResult, HxError};

/// A polynomial in `_x`, coefficients from the constant term up. Every
/// coefficient is kept in the canonical form of SymPy's coefficient domain
/// (a cancelled rational function of the parameters).
pub(crate) type XPoly = Vec<Expr>;

fn zero() -> Expr { Expr::zero() }

fn trim(mut p: XPoly) -> XPoly {
    while p.len() > 1 && p.last().is_some_and(Expr::is_zero) {
        p.pop();
    }
    if p.is_empty() {
        p.push(zero());
    }
    p
}

/// The canonical domain element of a coefficient.
pub(crate) fn norm_coeff(c: &Expr) -> Expr {
    if c.is_number() && c.as_number().is_some() {
        return c.clone();
    }
    let e = crate::expand::expand(c);
    if e.as_number().is_some() {
        return e;
    }
    crate::poly::api::cancel(&e).unwrap_or(e)
}

fn normalize(p: XPoly) -> XPoly {
    trim(p.iter().map(norm_coeff).collect())
}

fn padd(a: &XPoly, b: &XPoly) -> XPoly {
    let n = a.len().max(b.len());
    (0..n)
        .map(|k| match (a.get(k), b.get(k)) {
            (Some(x), Some(y)) => add(&[x.clone(), y.clone()]),
            (Some(x), None) | (None, Some(x)) => x.clone(),
            (None, None) => zero(),
        })
        .collect()
}

fn pneg(a: &XPoly) -> XPoly { a.iter().map(neg).collect() }

fn psub(a: &XPoly, b: &XPoly) -> XPoly { padd(a, &pneg(b)) }

fn pmul(a: &XPoly, b: &XPoly) -> XPoly {
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] = add(&[out[i + j].clone(), mul(&[x.clone(), y.clone()])]);
        }
    }
    normalize(out)
}

fn pscale(a: &XPoly, c: &Expr) -> XPoly { a.iter().map(|x| mul(&[c.clone(), x.clone()])).collect() }

fn plus_const(a: &XPoly, c: &Expr) -> XPoly {
    let mut r = a.clone();
    r[0] = add(&[r[0].clone(), c.clone()]);
    r
}

/// `p(q(_x))` for `p` in an auxiliary variable.
fn compose(p: &XPoly, q: &XPoly) -> XPoly {
    let mut r: XPoly = vec![p[p.len() - 1].clone()];
    for k in (0..p.len() - 1).rev() {
        r = plus_const(&pmul(&r, q), &p[k]);
    }
    normalize(r)
}

/// `Poly(n.all_coeffs()[:-1], A)`: the polynomial without its constant
/// term, divided by `A`.
fn shift_down(p: &XPoly) -> XPoly {
    if p.len() <= 1 {
        return vec![zero()];
    }
    p[1..].to_vec()
}

fn lin(c0: Expr, c1: Expr) -> XPoly { normalize(vec![c0, c1]) }

/// An operator `sum c_k (z d/dz)^k`.
#[derive(Clone, Debug)]
pub(crate) struct Op {
    pub poly: XPoly,
}

impl Op {
    fn new(p: XPoly) -> Op { Op { poly: normalize(p) } }

    /// `MultOperator(p)` for a `p` free of `_x`.
    pub fn mult(c: &Expr) -> Op { Op::new(vec![c.clone()]) }

    pub fn shift_a(ai: &Expr) -> HResult<Op> {
        if ai.is_zero() {
            return Err(HxError::Value("Cannot increment zero upper index.".into()));
        }
        Ok(Op::new(vec![Expr::one(), div(&Expr::one(), ai)]))
    }

    pub fn shift_b(bi: &Expr) -> HResult<Op> {
        if bi.is_one() {
            return Err(HxError::Value("Cannot decrement unit lower index.".into()));
        }
        Ok(Op::new(vec![Expr::one(), div(&Expr::one(), &sub(bi, &Expr::one()))]))
    }

    pub fn unshift_a(ap: &[Expr], bq: &[Expr], i: usize, z: &Expr) -> HResult<Op> {
        let mut ap = ap.to_vec();
        let ai = sub(&ap.remove(i), &Expr::one());
        if ai.is_zero() {
            return Err(HxError::Value("Cannot decrement unit upper index.".into()));
        }
        let mut m: XPoly = vec![mul(&[z.clone(), ai.clone()])];
        for a in &ap {
            m = pmul(&m, &lin(a.clone(), Expr::one()));
        }
        let d = lin(neg(&ai), ai.clone());
        let mut n = d.clone();
        for b in bq {
            n = pmul(&n, &normalize(plus_const(&d, &sub(b, &Expr::one()))));
        }
        let b0 = norm_coeff(&neg(&n[0]));
        if b0.is_zero() {
            return Err(HxError::Value("Cannot decrement upper index: cancels with lower".into()));
        }
        let n = compose(&shift_down(&n), &lin(Expr::one(), div(&Expr::one(), &ai)));
        Ok(Op::new(pscale(&psub(&n, &m), &div(&Expr::one(), &b0))))
    }

    pub fn unshift_b(ap: &[Expr], bq: &[Expr], i: usize, z: &Expr) -> HResult<Op> {
        let mut bq = bq.to_vec();
        let bi = add(&[bq.remove(i), Expr::one()]);
        if bi.is_zero() {
            return Err(HxError::Value("Cannot increment -1 lower index.".into()));
        }
        let bm1 = sub(&bi, &Expr::one());
        let mut m: XPoly = lin(zero(), bm1.clone());
        for b in &bq {
            m = pmul(&m, &lin(sub(b, &Expr::one()), Expr::one()));
        }
        let d = lin(neg(&bm1), bm1.clone());
        let mut n: XPoly = vec![z.clone()];
        for a in ap {
            n = pmul(&n, &normalize(plus_const(&d, a)));
        }
        let b0 = norm_coeff(&n[0]);
        if b0.is_zero() {
            return Err(HxError::Value("Cannot increment index: cancels with upper".into()));
        }
        let n = compose(&shift_down(&n), &lin(Expr::one(), div(&Expr::one(), &bm1)));
        Ok(Op::new(pscale(&psub(&m, &n), &div(&Expr::one(), &b0))))
    }

    pub fn meijer_shift_a(bi: &Expr) -> Op { Op::new(vec![bi.clone(), Expr::minus_one()]) }
    pub fn meijer_shift_b(bi: &Expr) -> Op { Op::new(vec![sub(&Expr::one(), bi), Expr::one()]) }
    pub fn meijer_shift_c(bi: &Expr) -> Op { Op::new(vec![neg(bi), Expr::one()]) }
    pub fn meijer_shift_d(bi: &Expr) -> Op { Op::new(vec![sub(bi, &Expr::one()), Expr::minus_one()]) }

    /// `MeijerUnShiftA` (decrement upper b index `i`).
    pub fn meijer_unshift_a(an: &[Expr], ap: &[Expr], bm: &[Expr], bq: &[Expr], i: usize, z: &Expr) -> HResult<Op> {
        let mut bm = bm.to_vec();
        let bi = sub(&bm.remove(i), &Expr::one());
        let mut m: XPoly = vec![Expr::one()];
        for b in &bm {
            m = pmul(&m, &lin(b.clone(), Expr::minus_one()));
        }
        for b in bq {
            m = pmul(&m, &lin(neg(b), Expr::one()));
        }
        let d = lin(bi.clone(), Expr::minus_one());
        let mut n: XPoly = vec![z.clone()];
        for a in an {
            n = pmul(&n, &normalize(plus_const(&d, &sub(&Expr::one(), a))));
        }
        for a in ap {
            n = pmul(&n, &normalize(plus_const(&pneg(&d), &sub(a, &Expr::one()))));
        }
        let b0 = norm_coeff(&n[0]);
        if b0.is_zero() {
            return Err(HxError::Value("Cannot decrement upper b index (cancels)".into()));
        }
        let n = compose(&shift_down(&n), &lin(bi, Expr::minus_one()));
        Ok(Op::new(pscale(&psub(&m, &n), &div(&Expr::one(), &b0))))
    }

    /// `MeijerUnShiftB` (increment upper a index `i`).
    pub fn meijer_unshift_b(an: &[Expr], ap: &[Expr], bm: &[Expr], bq: &[Expr], i: usize, z: &Expr) -> HResult<Op> {
        let mut an = an.to_vec();
        let ai = add(&[an.remove(i), Expr::one()]);
        let mut m: XPoly = vec![z.clone()];
        for a in &an {
            m = pmul(&m, &lin(sub(&Expr::one(), a), Expr::one()));
        }
        for a in ap {
            m = pmul(&m, &lin(sub(a, &Expr::one()), Expr::minus_one()));
        }
        let d = lin(sub(&ai, &Expr::one()), Expr::one());
        let mut n: XPoly = vec![Expr::one()];
        for b in bm {
            n = pmul(&n, &normalize(plus_const(&pneg(&d), b)));
        }
        for b in bq {
            n = pmul(&n, &normalize(plus_const(&d, &neg(b))));
        }
        let b0 = norm_coeff(&n[0]);
        if b0.is_zero() {
            return Err(HxError::Value("Cannot increment upper a index (cancels)".into()));
        }
        let n = compose(&shift_down(&n), &lin(sub(&Expr::one(), &ai), Expr::one()));
        Ok(Op::new(pscale(&psub(&m, &n), &div(&Expr::one(), &b0))))
    }

    /// `MeijerUnShiftC` (decrement lower b index `i`).
    pub fn meijer_unshift_c(an: &[Expr], ap: &[Expr], bm: &[Expr], bq: &[Expr], i: usize, z: &Expr) -> HResult<Op> {
        let mut bq = bq.to_vec();
        let bi = sub(&bq.remove(i), &Expr::one());
        let mut m: XPoly = vec![Expr::one()];
        for b in bm {
            m = pmul(&m, &lin(b.clone(), Expr::minus_one()));
        }
        for b in &bq {
            m = pmul(&m, &lin(neg(b), Expr::one()));
        }
        let d = lin(bi.clone(), Expr::one());
        let mut n: XPoly = vec![z.clone()];
        for a in an {
            n = pmul(&n, &normalize(plus_const(&d, &sub(&Expr::one(), a))));
        }
        for a in ap {
            n = pmul(&n, &normalize(plus_const(&pneg(&d), &sub(a, &Expr::one()))));
        }
        let b0 = norm_coeff(&n[0]);
        if b0.is_zero() {
            return Err(HxError::Value("Cannot decrement lower b index (cancels)".into()));
        }
        let n = compose(&shift_down(&n), &lin(neg(&bi), Expr::one()));
        Ok(Op::new(pscale(&psub(&m, &n), &div(&Expr::one(), &b0))))
    }

    /// `MeijerUnShiftD` (increment lower a index `i`).
    pub fn meijer_unshift_d(an: &[Expr], ap: &[Expr], bm: &[Expr], bq: &[Expr], i: usize, z: &Expr) -> HResult<Op> {
        let mut ap = ap.to_vec();
        let ai = add(&[ap.remove(i), Expr::one()]);
        let mut m: XPoly = vec![z.clone()];
        for a in an {
            m = pmul(&m, &lin(sub(&Expr::one(), a), Expr::one()));
        }
        for a in &ap {
            m = pmul(&m, &lin(sub(a, &Expr::one()), Expr::minus_one()));
        }
        let d = lin(sub(&ai, &Expr::one()), Expr::minus_one());
        let mut n: XPoly = vec![Expr::one()];
        for b in bm {
            n = pmul(&n, &normalize(plus_const(&pneg(&d), b)));
        }
        for b in bq {
            n = pmul(&n, &normalize(plus_const(&d, &neg(b))));
        }
        let b0 = norm_coeff(&n[0]);
        if b0.is_zero() {
            return Err(HxError::Value("Cannot increment lower a index (cancels)".into()));
        }
        let n = compose(&shift_down(&n), &lin(sub(&ai, &Expr::one()), Expr::minus_one()));
        Ok(Op::new(pscale(&psub(&m, &n), &div(&Expr::one(), &b0))))
    }

    /// `ReduceOrder(ai, bj)`.
    pub fn reduce_order(ai: &Expr, bj: &Expr) -> Option<Op> {
        let n = sub(ai, bj);
        if !is_integer_number(&n) {
            return None;
        }
        let n = n.as_i64()?;
        if n < 0 {
            return None;
        }
        if is_nonpositive_integer(bj) {
            return None;
        }
        let mut p: XPoly = vec![Expr::one()];
        for k in 0..n {
            let c = add(&[bj.clone(), Expr::int(k)]);
            p = pmul(&p, &lin(Expr::one(), div(&Expr::one(), &c)));
        }
        Some(Op::new(p))
    }

    /// `ReduceOrder._meijer(b, a, sign)`.
    fn meijer(b: &Expr, a: &Expr, sign: i64) -> Option<Op> {
        let n = sub(b, a);
        if crate::assumptions::is_negative(&n) == Some(true) || !is_integer_number(&n) {
            return None;
        }
        let n = n.as_i64()?;
        let mut p: XPoly = vec![Expr::one()];
        for k in 0..n {
            p = pmul(&p, &lin(add(&[a.clone(), Expr::int(k)]), Expr::int(sign)));
        }
        Some(Op::new(p))
    }

    pub fn meijer_minus(b: &Expr, a: &Expr) -> Option<Op> { Op::meijer(b, a, -1) }

    pub fn meijer_plus(a: &Expr, b: &Expr) -> Option<Op> {
        Op::meijer(&sub(&Expr::one(), a), &sub(&Expr::one(), b), 1)
    }

    /// `Poly(op._poly.subs({z: zz, _x: -_x}), _x)` (the second Slater sum).
    pub fn reflected(&self, z: &Expr, zz: &Expr) -> Op {
        let pairs = [(z.clone(), zz.clone())];
        Op::new(
            self.poly
                .iter()
                .enumerate()
                .map(|(k, c)| {
                    let c = crate::subs::subs(c, &pairs);
                    if k % 2 == 1 { neg(&c) } else { c }
                })
                .collect(),
        )
    }
}

/// Objects operators act on: scalars and coefficient row vectors.
pub(crate) trait Lin: Clone {
    fn scale(&self, c: &Expr) -> Self;
    fn plus(&self, o: &Self) -> Self;
}

impl Lin for Expr {
    fn scale(&self, c: &Expr) -> Self { mul(&[c.clone(), self.clone()]) }
    fn plus(&self, o: &Self) -> Self { add(&[self.clone(), o.clone()]) }
}

impl Lin for Vec<Expr> {
    fn scale(&self, c: &Expr) -> Self { self.iter().map(|e| mul(&[c.clone(), e.clone()])).collect() }
    fn plus(&self, o: &Self) -> Self { self.iter().zip(o).map(|(a, b)| add(&[a.clone(), b.clone()])).collect() }
}

/// `Operator.apply(obj, op)`.
pub(crate) fn apply<T: Lin>(o: &Op, obj: &T, f: &mut dyn FnMut(&T) -> HResult<T>) -> HResult<T> {
    let coeffs = &o.poly;
    let mut diffs = vec![obj.clone()];
    for _ in 1..coeffs.len() {
        let next = f(diffs.last().unwrap())?;
        diffs.push(next);
    }
    let mut r = diffs[0].scale(&coeffs[0]);
    for (c, d) in coeffs[1..].iter().zip(&diffs[1..]) {
        r = r.plus(&d.scale(c));
    }
    Ok(r)
}

/// `apply_operators(obj, ops, op)`: the last operator acts first.
pub(crate) fn apply_operators<T: Lin>(obj: T, ops: &[Op], f: &mut dyn FnMut(&T) -> HResult<T>) -> HResult<T> {
    let mut res = obj;
    for o in ops.iter().rev() {
        res = apply(o, &res, f)?;
    }
    Ok(res)
}

/// `make_simp(z)`: cancel a rational function of `z` over the domain of
/// its coefficients (`poly(numer, z).cancel(poly(denom, z))`).
pub(crate) fn make_simp(e: &Expr, z: &Expr) -> HResult<Expr> {
    use crate::poly::polyobj::Poly;
    let (n, d) = crate::simptools::as_numer_denom(e);
    let n = crate::expand::expand(&n);
    let ps = Poly::from_exprs(&[n.clone(), d.clone()], std::slice::from_ref(z))
        .map_err(|err| HxError::Polynomial(format!("{err:?}")))?;
    let (f, g) = (&ps[0], &ps[1]);
    if f.is_zero() {
        return Ok(Expr::zero());
    }
    let (cq, fz) = f.rep.clear_denoms();
    let (cp, gz) = g.rep.clear_denoms();
    let (_, pp, qq) = crate::poly::gcd::gcd_mz(&fz, &gz);
    let gg = cp.gcd(&cq);
    let cp = cp.div_exact_unchecked(&gg);
    let cq = cq.div_exact_unchecked(&gg);
    let mut c = Expr::num(crate::number::Number::rational(cp.to_big(), cq.to_big()).expect("nonzero denominator"));
    let mut pp = pp.to_rat();
    let mut qq = qq.to_rat();
    let p_neg = pp.lc().is_negative();
    let q_neg = qq.lc().is_negative();
    if p_neg && q_neg {
        pp = pp.neg();
        qq = qq.neg();
    } else if p_neg {
        c = neg(&c);
        pp = pp.neg();
    } else if q_neg {
        c = neg(&c);
        qq = qq.neg();
    }
    let pe = z_expr(&f.with_rep(pp), z);
    let qe = z_expr(&g.with_rep(qq), z);
    Ok(div(&mul(&[c, pe]), &qe))
}

/// `Poly.as_expr()` for a polynomial in `z` over a coefficient ring: each
/// power of `z` multiplies its (expanded) coefficient.
fn z_expr(p: &crate::poly::polyobj::Poly, z: &Expr) -> Expr {
    let terms: Vec<Expr> = p
        .terms()
        .into_iter()
        .map(|(m, c)| mul(&[c, crate::ops::pow(z, &Expr::int(m[0] as i64))]))
        .collect();
    add(&terms)
}
