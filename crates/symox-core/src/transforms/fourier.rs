//! Fourier, sine and cosine transforms with SymPy's normalisations:
//!
//! * `fourier_transform`: `F(k) = int_{-oo}^{oo} f(x) exp(-2*pi*I*x*k) dx`,
//! * `inverse_fourier_transform`: same with `exp(2*pi*I*x*k)`,
//! * `sine_transform` / `cosine_transform` and their inverses:
//!   `sqrt(2/pi) * int_0^oo f(x) sin(x*k) dx` (resp. `cos`).
//!
//! The engine evaluates `J(w) = int_{-oo}^{oo} g(x) exp(I*w*x) dx` for the
//! classical families (delta distributions, Gaussians times polynomials and
//! trigonometric factors, one sided and two sided exponentials, finite
//! intervals, Lorentzians, the sinc function, frequency shifts) and the half
//! line integrals `int_0^oo x**(mu-1) exp(-c*x) exp(I*w*x) dx`.

use super::util::*;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{abs, exp, gamma, sign};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};
use crate::special::SpecialFn;
use crate::symbol::Assumptions;

type V = (Expr, Expr);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FourierKind {
    Fourier,
    InverseFourier,
    Sine,
    InverseSine,
    Cosine,
    InverseCosine,
}

impl FourierKind {
    pub fn from_head(s: SpecialFn) -> Option<Self> {
        use SpecialFn::*;
        Some(match s {
            FourierTransform => FourierKind::Fourier,
            InverseFourierTransform => FourierKind::InverseFourier,
            SineTransform => FourierKind::Sine,
            InverseSineTransform => FourierKind::InverseSine,
            CosineTransform => FourierKind::Cosine,
            InverseCosineTransform => FourierKind::InverseCosine,
            _ => return None,
        })
    }

    /// SymPy's constants `(a, b)` of the kernel.
    fn ab(self) -> (Expr, Expr) {
        match self {
            FourierKind::Fourier => (Expr::one(), mul(&[int(-2), Expr::pi()])),
            FourierKind::InverseFourier => (Expr::one(), mul(&[int(2), Expr::pi()])),
            _ => (div(&sqrt(&int(2)), &sqrt(&Expr::pi())), Expr::one()),
        }
    }
}

/// Direct evaluation of a Fourier type transform: `Some((F, cond))`, or
/// `None` when no rule applies (SymPy's `IntegralTransformError`).
pub fn fourier_type_transform(kind: FourierKind, f: &Expr, x: &Expr, k: &Expr, simplify: bool) -> Option<V> {
    let (a, b) = kind.ab();
    let w = mul(&[b, k.clone()]);
    let (v, c) = match kind {
        FourierKind::Fourier | FourierKind::InverseFourier => full(f, x, &w, 0)?,
        FourierKind::Sine | FourierKind::InverseSine => sincos(f, x, &w, true, 0)?,
        FourierKind::Cosine | FourierKind::InverseCosine => sincos(f, x, &w, false, 0)?,
    };
    let r = mul(&[a, v]);
    if r.has_nonfinite() {
        return None;
    }
    Some((if simplify { finish(&r) } else { r }, c))
}

/// Simplification of a transform result.
pub(crate) fn finish(e: &Expr) -> Expr {
    let simp = crate::simplify::simplify;
    let first = exp_to_trig(&simp(e));
    let expanded = exp_to_trig(&expand_exp_args(&crate::expand::expand(e)));
    let mut cands = vec![first.clone(), simp(&expanded)];
    let gam = super::mellin::gamma_reflect(&super::mellin::gamma_duplicate(&first));
    if gam != first {
        cands.push(simp(&gam));
    }
    let best = cands.into_iter().min_by_key(size).unwrap();
    if size(&best) <= size(e) + 2 { best } else { e.clone() }
}

/// Expand the arguments of the exponentials of `e`.
fn expand_exp_args(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| expand_exp_args(c));
    match e.kind() {
        Kind::Func(Func::Exp, a) => exp(&crate::expand::expand(&a[0])),
        _ => e,
    }
}

/// Combine `c*exp(I*u) + c*exp(-I*u)` into `2*c*cos(u)` and
/// `c*exp(I*u) - c*exp(-I*u)` into `2*I*c*sin(u)`, bottom up.
pub(crate) fn exp_to_trig(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| exp_to_trig(c));
    if !e.is_add() {
        return e;
    }
    let i = Expr::i();
    // (coefficient, exponent u where exp(I*u), rest)
    let mut items: Vec<Option<(Expr, Expr, Expr)>> = Vec::new();
    let mut others = Vec::new();
    for t in e.args() {
        let mut ex = None;
        let mut rest = Vec::new();
        for fct in factors(&t) {
            match fct.kind() {
                Kind::Func(Func::Exp, a) if ex.is_none() => {
                    let u = crate::expand::expand(&div(&a[0], &i));
                    if !crate::traverse::has(&u, &i) {
                        ex = Some(u);
                        continue;
                    }
                    rest.push(fct.clone())
                }
                _ => rest.push(fct.clone()),
            }
        }
        match ex {
            Some(u) => {
                let r = mul(&rest);
                let (c, r) = match r.kind() {
                    Kind::Mul(..) => {
                        let fs = r.args();
                        if matches!(fs[0].kind(), Kind::Num(_)) {
                            (fs[0].clone(), mul(&fs[1..]))
                        } else {
                            (Expr::one(), r.clone())
                        }
                    }
                    Kind::Num(_) => (r.clone(), Expr::one()),
                    _ => (Expr::one(), r.clone()),
                };
                items.push(Some((c, u, r)))
            }
            None => others.push(t),
        }
    }
    let mut out = others;
    let mut changed = false;
    for j in 0..items.len() {
        let Some((c, u, r)) = items[j].clone() else { continue };
        let mut paired = false;
        for l in j + 1..items.len() {
            let Some((c2, u2, r2)) = items[l].clone() else { continue };
            if r2 != r || crate::expand::expand(&add(&[u.clone(), u2.clone()])) != Expr::zero() {
                continue;
            }
            if c2 == c {
                out.push(mul(&[int(2), c.clone(), r.clone(), crate::functions::cos(&u)]));
            } else if crate::expand::expand(&add(&[c.clone(), c2.clone()])).is_zero() {
                out.push(mul(&[int(2), i.clone(), c.clone(), r.clone(), crate::functions::sin(&u)]));
            } else {
                continue;
            }
            items[l] = None;
            paired = true;
            changed = true;
            break;
        }
        if !paired {
            out.push(mul(&[c, r, exp(&mul(&[i.clone(), u]))]));
        }
    }
    if changed { add(&out) } else { e }
}

fn and_all(c: Vec<Expr>) -> Expr {
    crate::logic::and(c)
}

/// Positive dummy used for the half line integrals.
fn pos_dummy(name: &str) -> Expr {
    Expr::dummy(name, Assumptions::POSITIVE)
}

/// Evaluate `Heaviside` and `sign` of a sign definite argument (after the
/// substitution of a positive dummy).
fn settle(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| settle(c));
    match e.kind() {
        Kind::Func(Func::Heaviside, a) if pos(&a[0]) => Expr::one(),
        Kind::Func(Func::Heaviside, a) if negv(&a[0]) => Expr::zero(),
        Kind::Func(Func::Sign, a) if pos(&a[0]) => Expr::one(),
        Kind::Func(Func::Sign, a) if negv(&a[0]) => int(-1),
        _ => e,
    }
}

/// `int_{-oo}^{oo} f(x) exp(I*w*x) dx`.
pub(crate) fn full(f: &Expr, x: &Expr, w: &Expr, depth: u32) -> Option<V> {
    if depth > MAX_DEPTH || !has(f, x) {
        return None;
    }
    if f.is_add() {
        if let Some(r) = heaviside_sum(f, x, w) {
            return Some(r);
        }
        let mut vals = Vec::new();
        let mut conds = Vec::new();
        for t in terms(f) {
            let (v, c) = full(&t, x, w, depth + 1)?;
            vals.push(v);
            conds.push(c);
        }
        return Some((add(&vals), and_all(conds)));
    }
    let (c, g) = split_mul(f, x);
    if !c.is_one() {
        let (v, cd) = full(&g, x, w, depth + 1)?;
        return Some((mul(&[c, v]), cd));
    }
    if let Some(r) = delta(&g, x, w) {
        return Some(r);
    }
    if let Some(r) = gaussian(&g, x, w) {
        return Some(r);
    }
    if let Some(r) = heaviside_rule(&g, x, w) {
        return Some(r);
    }
    if let Some(r) = halves(&g, x, w) {
        return Some(r);
    }
    if let Some(r) = lorentz(&g, x, w) {
        return Some(r);
    }
    if let Some(r) = sinc(&g, x, w) {
        return Some(r);
    }
    if let Some(r) = shift(&g, x, w, depth) {
        return Some(r);
    }
    if g.is_mul() && factors(&g).iter().any(|a| a.is_add()) {
        let e = crate::expand::expand(&g);
        if e.is_add() {
            return full(&e, x, w, depth + 1);
        }
    }
    None
}

/// `DiracDelta(alpha*x + beta, n) * h(x)`.
fn delta(g: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let fs = factors(g);
    let idx = fs.iter().position(|a| matches!(a.kind(), Kind::Func(Func::DiracDelta, _)))?;
    let (u, n) = match fs[idx].kind() {
        Kind::Func(Func::DiracDelta, a) => (a[0].clone(), a.get(1).and_then(|n| n.as_i64()).unwrap_or(0)),
        _ => return None,
    };
    if n < 0 {
        return None;
    }
    let (al, be) = linear(&u, x)?;
    let mut rest: Vec<Expr> = fs.clone();
    rest.remove(idx);
    let h = mul(&rest);
    if has_head(&h, &|k| matches!(k, Kind::Func(Func::DiracDelta, _))) {
        return None;
    }
    let x0 = neg(&div(&be, &al));
    let kern = mul(&[h, exp(&mul(&[Expr::i(), w.clone(), x.clone()]))]);
    let d = crate::diff::diff_n(&kern, x, n as u32);
    let at = crate::subs::subs1(&d, x, &x0);
    let fac = div(&pow(&neg(&inv_e(&al)), &int(n)), &abs(&al));
    Some((mul(&[fac, at]), tru()))
}

fn inv_e(e: &Expr) -> Expr {
    div(&Expr::one(), e)
}

/// One term `C * x**n * exp(q(x))` with `n` a non negative integer and `q`
/// a polynomial: returns `(C, n, coefficients of q)`.
fn exp_poly_term(t: &Expr, x: &Expr) -> Option<(Expr, u32, Vec<Expr>)> {
    let mut expo = Vec::new();
    let mut n = 0u32;
    let mut rest = Vec::new();
    for f in factors(t) {
        match f.kind() {
            Kind::Func(Func::Exp, a) => expo.push(a[0].clone()),
            Kind::Pow(b, e) if b.is_constant(Constant::E) => expo.push(e.clone()),
            _ => {
                if let Some(p) = power_of(&f, x) {
                    let m = p.as_i64()?;
                    if m < 0 {
                        return None;
                    }
                    n += m as u32;
                } else if has(&f, x) {
                    return None;
                } else {
                    rest.push(f)
                }
            }
        }
    }
    let q = poly(&crate::expand::expand(&add(&expo)), x)?;
    Some((mul(&rest), n, q))
}

/// Gaussians times polynomials (and trigonometric factors).
fn gaussian(g: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let mut quad = false;
    crate::traverse::preorder(g, &mut |n| {
        if let Kind::Func(Func::Exp, a) = n.kind() {
            if poly(&crate::expand::expand(&a[0]), x).is_some_and(|c| c.len() == 3) {
                quad = true;
            }
        }
    });
    if !quad {
        return None;
    }
    let e = crate::expand::expand(&rewrite_exp(g));
    let beta = Expr::dummy("beta", Assumptions::NONE);
    let mut vals = Vec::new();
    let mut conds = Vec::new();
    for t in terms(&e) {
        let (cf, n, q) = exp_poly_term(&t, x)?;
        if q.len() != 3 {
            return None;
        }
        let c = neg(&q[2]);
        let rc = crate::functions::re(&c);
        if negv(&rc) || rc.is_zero() || nonpos(&c) {
            return None;
        }
        if !re_pos(&c) {
            conds.push(crate::relational::gt(&rc, &Expr::zero()));
        }
        let base = mul(&[div(&sqrt(&Expr::pi()), &sqrt(&c)), exp(&div(&pow(&beta, &int(2)), &mul(&[int(4), c.clone()])))]);
        let d = crate::diff::diff_n(&base, &beta, n);
        let b = add(&[q[1].clone(), mul(&[Expr::i(), w.clone()])]);
        vals.push(mul(&[cf, exp(&q[0]), crate::subs::subs1(&d, &beta, &b)]));
    }
    Some((add(&vals), and_all(conds)))
}

/// `int_0^oo h(y) exp(I*w*y) dy` for sums of `C * y**nu * exp(-lambda*y)`;
/// `regular` allows undamped integer powers (Abel regularisation, used for
/// finite intervals and balanced Heaviside sums).
fn half(h: &Expr, y: &Expr, w: &Expr, regular: bool) -> Option<V> {
    let e = crate::expand::expand(&rewrite_exp(h));
    let mut vals = Vec::new();
    let mut conds = Vec::new();
    for t in terms(&e) {
        if t.is_zero() {
            continue;
        }
        let ((nu, cf), c1) = as_exp_term_pow(&t, y)?;
        let lam = neg(&c1);
        let pp = sub(&lam, &mul(&[Expr::i(), w.clone()]));
        let rl = crate::functions::re(&lam);
        let oscillatory = lam.is_zero() || rl.is_zero();
        if negv(&rl) || (oscillatory && !regular) {
            return None;
        }
        if !oscillatory && !re_pos(&lam) {
            conds.push(crate::relational::gt(&rl, &Expr::zero()));
        }
        match nu.as_i64() {
            Some(m) if m < 0 => return None,
            Some(_) => {}
            None => {
                if oscillatory {
                    return None;
                }
                let c = re_gt(&nu, &int(-1));
                match truth(&c) {
                    Some(false) => return None,
                    Some(true) => {}
                    None => conds.push(c),
                }
            }
        }
        let n1 = add(&[nu, Expr::one()]);
        vals.push(mul(&[cf, gamma(&n1), pow(&pp, &neg(&n1))]));
    }
    Some((add(&vals), and_all(conds)))
}

/// `term == C * y**nu * exp(c1*y + c0)`: returns `((nu, C*exp(c0)), c1)`.
fn as_exp_term_pow(t: &Expr, y: &Expr) -> Option<((Expr, Expr), Expr)> {
    let (nu, rest) = split_power(t, y);
    let (cf, c1) = as_exp_term(&rest, y)?;
    Some(((nu, cf), c1))
}

/// Functions involving `Abs(x)`, `Heaviside(+-x)` or `sign(x)`: split the
/// real line into the two half lines.
fn halves(g: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let xs = x.clone();
    if !has_head(g, &|k| match k {
        Kind::Func(Func::Abs | Func::Heaviside | Func::Sign, a) => has(&a[0], &xs),
        _ => false,
    }) {
        return None;
    }
    let y = pos_dummy("y");
    let gp = settle(&crate::subs::subs1(g, x, &y));
    let gn = settle(&crate::subs::subs1(g, x, &neg(&y)));
    let ys = y.clone();
    let bad = |e: &Expr| {
        has_head(e, &|k| match k {
            Kind::Func(Func::Abs | Func::Heaviside | Func::Sign, a) => has(&a[0], &ys),
            _ => false,
        })
    };
    if bad(&gp) || bad(&gn) {
        return None;
    }
    let mut vals = Vec::new();
    let mut conds = Vec::new();
    if !gp.is_zero() {
        let (v, c) = half(&gp, &y, w, false)?;
        vals.push(v);
        conds.push(c);
    }
    if !gn.is_zero() {
        let (v, c) = half(&gn, &y, &neg(w), false)?;
        vals.push(v);
        conds.push(c);
    }
    Some((add(&vals), and_all(conds)))
}

/// `Heaviside(alpha*x + beta)` factors with a non zero offset: rays and
/// finite intervals.
fn heaviside_rule(g: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let fs = factors(g);
    let mut lo: Option<Expr> = None;
    let mut hi: Option<Expr> = None;
    let mut rest = Vec::new();
    let mut found = false;
    for f in &fs {
        if let Kind::Func(Func::Heaviside, a) = f.kind() {
            if a.len() == 1 {
                if let Some((al, be)) = linear(&a[0], x) {
                    if !be.is_zero() {
                        let x0 = neg(&div(&be, &al));
                        found = true;
                        if pos(&al) {
                            lo = Some(match lo {
                                None => x0,
                                Some(l) => pick(&l, &x0, true)?,
                            });
                            continue;
                        } else if negv(&al) {
                            hi = Some(match hi {
                                None => x0,
                                Some(h) => pick(&h, &x0, false)?,
                            });
                            continue;
                        }
                        return None;
                    }
                }
            }
        }
        rest.push(f.clone());
    }
    if !found {
        return None;
    }
    let h = mul(&rest);
    let xs = x.clone();
    if has_head(&h, &|k| match k {
        Kind::Func(Func::Abs | Func::Heaviside | Func::Sign, a) => has(&a[0], &xs),
        _ => false,
    }) {
        return None;
    }
    let y = pos_dummy("y");
    let iw = mul(&[Expr::i(), w.clone()]);
    match (lo, hi) {
        (Some(l), None) => {
            let (v, c) = half(&crate::subs::subs1(&h, x, &add(&[y.clone(), l.clone()])), &y, w, false)?;
            Some((mul(&[exp(&mul(&[iw, l])), v]), c))
        }
        (None, Some(u)) => {
            let (v, c) = half(&crate::subs::subs1(&h, x, &sub(&u, &y)), &y, &neg(w), false)?;
            Some((mul(&[exp(&mul(&[iw, u])), v]), c))
        }
        (Some(l), Some(u)) => {
            let d = sub(&u, &l);
            if nonpos(&d) {
                return Some((Expr::zero(), tru()));
            }
            if !pos(&d) {
                return None;
            }
            let (v1, _) = half(&crate::subs::subs1(&h, x, &add(&[y.clone(), l.clone()])), &y, w, true)?;
            let (v2, _) = half(&crate::subs::subs1(&h, x, &add(&[y.clone(), u.clone()])), &y, w, true)?;
            Some((sub(&mul(&[exp(&mul(&[iw.clone(), l])), v1]), &mul(&[exp(&mul(&[iw, u])), v2])), tru()))
        }
        (None, None) => None,
    }
}

/// Larger (`max`) or smaller of two decidably ordered endpoints.
fn pick(a: &Expr, b: &Expr, max: bool) -> Option<Expr> {
    let d = sub(a, b);
    if nonneg(&d) {
        Some(if max { a.clone() } else { b.clone() })
    } else if nonpos(&d) {
        Some(if max { b.clone() } else { a.clone() })
    } else {
        None
    }
}

/// Balanced sums `sum c_i * Heaviside(alpha_i*x + beta_i)` (rectangles).
fn heaviside_sum(f: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let mut right = Vec::new();
    let mut left = Vec::new();
    let mut vals = Vec::new();
    let iw = mul(&[Expr::i(), w.clone()]);
    for t in terms(f) {
        let (c, h) = split_mul(&t, x);
        let Kind::Func(Func::Heaviside, a) = h.kind() else { return None };
        let (al, be) = linear(&a[0], x)?;
        let x0 = neg(&div(&be, &al));
        let e = exp(&mul(&[iw.clone(), x0]));
        if pos(&al) {
            // int_{x0}^oo exp(I*w*x) dx = -exp(I*w*x0)/(I*w)
            vals.push(neg(&div(&mul(&[c.clone(), e]), &iw)));
            right.push(c);
        } else if negv(&al) {
            vals.push(div(&mul(&[c.clone(), e]), &iw));
            left.push(c);
        } else {
            return None;
        }
    }
    if !crate::expand::expand(&add(&right)).is_zero() || !crate::expand::expand(&add(&left)).is_zero() {
        return None;
    }
    Some((add(&vals), tru()))
}

/// `(c1*x + c0) / (alpha*x**2 + beta*x + gamma)` with complex roots.
fn lorentz(g: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let (n, d) = crate::simplify::numer_denom(g);
    if !has(&d, x) {
        return None;
    }
    let nc = poly(&crate::expand::expand(&n), x)?;
    let dc = poly(&crate::expand::expand(&d), x)?;
    if dc.len() != 3 || nc.len() > 2 {
        return None;
    }
    let (ga, be, al) = (dc[0].clone(), dc[1].clone(), dc[2].clone());
    let disc = sub(&pow(&be, &int(2)), &mul(&[int(4), al.clone(), ga.clone()]));
    if !negv(&disc) || !real(&al) || !real(&be) {
        return None;
    }
    let c0 = nc[0].clone();
    let c1 = nc.get(1).cloned().unwrap_or_else(Expr::zero);
    let x0 = neg(&div(&be, &mul(&[int(2), al.clone()])));
    let q = div(&neg(&disc), &mul(&[int(4), pow(&al, &int(2))]));
    let sq = sqrt(&q);
    let damp = exp(&neg(&mul(&[sq.clone(), abs(w)])));
    let pi = Expr::pi();
    let v = add(&[
        div(&mul(&[add(&[c0, mul(&[c1.clone(), x0.clone()])]), pi.clone()]), &sq),
        mul(&[c1, Expr::i(), pi, sign(w)]),
    ]);
    let r = mul(&[exp(&mul(&[Expr::i(), w.clone(), x0])), div(&v, &al), damp]);
    Some((r, tru()))
}

/// `w > a` (`greater`) or `w < a`, solved for the symbol when `w = c*k`
/// with a sign-definite constant `c`.
fn rel(w: &Expr, a: &Expr, greater: bool) -> Expr {
    let fs = factors(w);
    let syms: Vec<&Expr> = fs.iter().filter(|f| matches!(f.kind(), Kind::Sym(_))).collect();
    if syms.len() == 1 {
        let k = syms[0].clone();
        let c = div(w, &k);
        if !has(&c, &k) && (pos(&c) || negv(&c)) {
            let b = div(a, &c);
            let g = greater == pos(&c);
            return if g { crate::relational::gt(&k, &b) } else { crate::relational::lt(&k, &b) };
        }
    }
    if greater { crate::relational::gt(w, a) } else { crate::relational::lt(w, a) }
}

/// `sin(alpha*x)/x`.
fn sinc(g: &Expr, x: &Expr, w: &Expr) -> Option<V> {
    let fs = factors(g);
    if fs.len() != 2 {
        return None;
    }
    let mut al = None;
    let mut inv = false;
    for f in &fs {
        match f.kind() {
            Kind::Func(Func::Sin, a) => al = Some(scaled(&a[0], x)?),
            Kind::Pow(b, e) if b == x && e.as_i64() == Some(-1) => inv = true,
            _ => return None,
        }
    }
    let al = al?;
    if !inv || !real(&al) {
        return None;
    }
    let aa = abs(&al);
    let pi = Expr::pi();
    let r = crate::relational::piecewise(vec![
        (Expr::zero(), crate::logic::or(vec![rel(w, &aa, true), rel(w, &neg(&aa), false)])),
        (mul(&[sign(&al), pi]), tru()),
    ]);
    Some((r, tru()))
}

/// Frequency shift: `exp(I*beta*x) * h(x)` with real `beta` (trigonometric
/// factors are rewritten as exponentials).
fn shift(g: &Expr, x: &Expr, w: &Expr, depth: u32) -> Option<V> {
    let has_lin_exp = factors(g).iter().any(|f| match f.kind() {
        Kind::Func(Func::Exp, a) => linear_part(&a[0], x).is_some_and(|c| !c.is_zero()),
        _ => false,
    });
    if !has_trig(g) && !has_lin_exp {
        return None;
    }
    let e = crate::expand::expand(&rewrite_exp(g));
    let mut vals = Vec::new();
    let mut conds = Vec::new();
    for t in terms(&e) {
        let mut lin = Vec::new();
        let mut rest = Vec::new();
        for f in factors(&t) {
            match f.kind() {
                Kind::Func(Func::Exp, a) => {
                    let c = linear_part(&a[0], x)?;
                    lin.push(c.clone());
                    let r = crate::expand::expand(&sub(&a[0], &mul(&[c, x.clone()])));
                    rest.push(exp(&r));
                }
                _ => rest.push(f),
            }
        }
        let c = add(&lin);
        let beta = crate::expand::expand(&div(&c, &Expr::i()));
        if !beta.is_zero() && !real(&beta) {
            return None;
        }
        let h = mul(&rest);
        if has_trig(&h) {
            return None;
        }
        if beta.is_zero() && t == *g {
            return None;
        }
        let (v, cd) = full(&h, x, &add(&[w.clone(), beta]), depth + 1)?;
        vals.push(v);
        conds.push(cd);
    }
    Some((add(&vals), and_all(conds)))
}

/// Coefficient of `x` in the expanded exponent `u` (the rest may be any
/// function of `x` without a linear part).
fn linear_part(u: &Expr, x: &Expr) -> Option<Expr> {
    let u = crate::expand::expand(u);
    let mut c = Vec::new();
    for t in terms(&u) {
        let (k, r) = split_mul(&t, x);
        if r == *x {
            c.push(k)
        } else if !has(&t, x) {
            continue;
        } else if poly(&t, x).is_some_and(|p| p.len() == 2) {
            return None;
        }
    }
    Some(add(&c))
}

/// `sqrt(2/pi) * int_0^oo f(x) sin(w*x) dx` (or `cos`) without the
/// constant: the integral alone.
fn sincos(f: &Expr, x: &Expr, w: &Expr, sine: bool, depth: u32) -> Option<V> {
    if depth > MAX_DEPTH || !has(f, x) {
        return None;
    }
    if f.is_add() {
        let mut vals = Vec::new();
        let mut conds = Vec::new();
        for t in terms(f) {
            let (v, c) = sincos(&t, x, w, sine, depth + 1)?;
            vals.push(v);
            conds.push(c);
        }
        return Some((add(&vals), and_all(conds)));
    }
    let (c, g) = split_mul(f, x);
    if !c.is_one() {
        let (v, cd) = sincos(&g, x, w, sine, depth + 1)?;
        return Some((mul(&[c, v]), cd));
    }
    // parity: even functions for the cosine, odd ones for the sine
    let gm = crate::subs::subs1(&g, x, &neg(x));
    let par = if sine { add(&[gm.clone(), g.clone()]) } else { sub(&gm, &g) };
    if crate::expand::expand(&par).is_zero() {
        if let Some((v, cd)) = full(&g, x, w, depth + 1) {
            let d = if sine { mul(&[int(2), Expr::i()]) } else { int(2) };
            return Some((div(&v, &d), cd));
        }
    }
    let e = crate::expand::expand(&rewrite_exp(&g));
    let mut vals = Vec::new();
    let mut conds = Vec::new();
    for t in terms(&e) {
        let ((nu, cf), c1) = as_exp_term_pow(&t, x)?;
        let mu = add(&[nu, Expr::one()]);
        let lam = neg(&c1);
        let pi = Expr::pi();
        if lam.is_zero() {
            // pure powers
            let (lo, hi) = if sine { (int(-1), int(1)) } else { (int(0), int(1)) };
            for cd in [re_gt(&mu, &lo), crate::relational::lt(&crate::functions::re(&mu), &hi)] {
                match truth(&cd) {
                    Some(false) => return None,
                    Some(true) => {}
                    None => conds.push(cd),
                }
            }
            let v = if sine && mu.is_zero() {
                div(&pi, &int(2))
            } else {
                let tr = if sine {
                    crate::functions::sin(&div(&mul(&[pi, mu.clone()]), &int(2)))
                } else {
                    crate::functions::cos(&div(&mul(&[pi, mu.clone()]), &int(2)))
                };
                mul(&[gamma(&mu), tr, pow(w, &neg(&mu))])
            };
            vals.push(mul(&[cf, v]));
            continue;
        }
        let rl = crate::functions::re(&lam);
        if negv(&rl) || rl.is_zero() {
            return None;
        }
        if !re_pos(&lam) {
            conds.push(crate::relational::gt(&rl, &Expr::zero()));
        }
        let lo = if sine { int(-1) } else { int(0) };
        let cd = re_gt(&mu, &lo);
        match truth(&cd) {
            Some(false) => return None,
            Some(true) => {}
            None => conds.push(cd),
        }
        let iw = mul(&[Expr::i(), w.clone()]);
        let ip = mul(&[gamma(&mu), pow(&sub(&lam, &iw), &neg(&mu))]);
        let im = mul(&[gamma(&mu), pow(&add(&[lam.clone(), iw]), &neg(&mu))]);
        let v = if sine { div(&sub(&ip, &im), &mul(&[int(2), Expr::i()])) } else { div(&add(&[ip, im]), &int(2)) };
        vals.push(mul(&[cf, v]));
    }
    Some((add(&vals), and_all(conds)))
}

