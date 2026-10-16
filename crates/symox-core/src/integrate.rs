//! Symbolic integration (indefinite and definite).
//!
//! A rule based integrator in the spirit of SymPy's `manualintegrate`:
//! linearity, power rule, a table of elementary antiderivatives composed
//! with linear arguments, rational functions, trigonometric powers,
//! integration by parts (LIATE), cyclic exponential/trigonometric
//! products and u-substitution driven by the integrand's sub-expressions.

use rustc_hash::FxHashSet;

use crate::diff::diff;
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::*;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sqrt, sub};
use crate::relational::{ne, piecewise};
use crate::simplify::numer_denom;
use crate::subs::subs1;
use crate::traverse::has;

const MAX_DEPTH: u32 = 10;

/// Indefinite integral of `f` with respect to `x`. Returns an unevaluated
/// `Integral(f, x)` when no rule applies.
pub fn integrate(f: &Expr, x: &Expr) -> Expr {
    match try_integrate(f, x) {
        Some(r) => r,
        None => integral_unevaluated(f, x),
    }
}

/// Definite integral of `f` over `[a, b]` (via the antiderivative).
pub fn integrate_definite(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Expr {
    match try_integrate(f, x) {
        Some(big_f) => {
            let fb = subs1(&big_f, x, b);
            let fa = subs1(&big_f, x, a);
            if has_bad(&fb) || has_bad(&fa) {
                return crate::functions::undefined(
                    "Integral",
                    vec![f.clone(), crate::functions::undefined("Tuple", vec![x.clone(), a.clone(), b.clone()])],
                );
            }
            expand_if_numeric(&sub(&fb, &fa))
        }
        None => crate::functions::undefined(
            "Integral",
            vec![f.clone(), crate::functions::undefined("Tuple", vec![x.clone(), a.clone(), b.clone()])],
        ),
    }
}

fn has_bad(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::NaN | Constant::ComplexInfinity) => true,
        _ => e.args().iter().any(has_bad),
    }
}

fn expand_if_numeric(e: &Expr) -> Expr {
    if crate::functions::is_numeric_value(e) {
        let x = expand(e);
        if crate::traverse::node_count(&x) <= crate::traverse::node_count(e) {
            return x;
        }
    }
    e.clone()
}

/// `Integral(f, x)` placeholder.
pub fn integral_unevaluated(f: &Expr, x: &Expr) -> Expr {
    crate::functions::undefined("Integral", vec![f.clone(), x.clone()])
}

/// Antiderivative, or `None` if the integrator gives up.
pub fn try_integrate(f: &Expr, x: &Expr) -> Option<Expr> {
    let mut seen = FxHashSet::default();
    integ(f, x, 0, &mut seen)
}

#[inline]
fn free(e: &Expr, x: &Expr) -> bool {
    !has(e, x)
}

/// `u == a*x + b` with `a != 0`, `a` and `b` free of `x`.
fn linear(u: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    if free(u, x) {
        return None;
    }
    let a = diff(u, x);
    if !free(&a, x) || a.is_zero() {
        return None;
    }
    let b = expand(&sub(u, &mul2(&a, x)));
    if !free(&b, x) {
        return None;
    }
    Some((a, b))
}

/// Dense coefficients `[c0, c1, ...]` of `e` as a polynomial in `x`.
pub fn poly_coeffs(e: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    if free(e, x) {
        return Some(vec![e.clone()]);
    }
    let ex = expand(e);
    let terms: Vec<Expr> = match ex.kind() {
        Kind::Add(..) => ex.args(),
        _ => vec![ex.clone()],
    };
    let mut coeffs: Vec<Vec<Expr>> = Vec::new();
    for t in terms {
        let (deg, c) = monomial_degree(&t, x)?;
        if deg > 10_000 {
            return None;
        }
        if coeffs.len() <= deg {
            coeffs.resize(deg + 1, Vec::new());
        }
        coeffs[deg].push(c);
    }
    Some(coeffs.into_iter().map(|v| add(&v)).collect())
}

fn monomial_degree(t: &Expr, x: &Expr) -> Option<(usize, Expr)> {
    if free(t, x) {
        return Some((0, t.clone()));
    }
    if t == x {
        return Some((1, Expr::one()));
    }
    match t.kind() {
        Kind::Pow(b, e) if b == x => {
            let k = e.as_i64()?;
            if k < 0 {
                return None;
            }
            Some((k as usize, Expr::one()))
        }
        Kind::Mul(c, facs) => {
            let mut deg = None;
            let mut rest = vec![Expr::num(c.clone())];
            for (b, e) in facs {
                if b == x {
                    let k = e.as_i64()?;
                    if k < 0 || deg.is_some() {
                        return None;
                    }
                    deg = Some(k as usize);
                } else {
                    let f = crate::ops::mul_factor(b, e);
                    if !free(&f, x) {
                        return None;
                    }
                    rest.push(f);
                }
            }
            Some((deg.unwrap_or(0), mul(&rest)))
        }
        _ => None,
    }
}

fn from_coeffs(c: &[Expr], x: &Expr) -> Expr {
    let terms: Vec<Expr> = c
        .iter()
        .enumerate()
        .filter(|(_, k)| !k.is_zero())
        .map(|(i, k)| mul2(k, &x.powi(i as i64)))
        .collect();
    add(&terms)
}

/// Polynomial division with exact coefficients: (quotient, remainder).
fn poly_divmod(n: &[Expr], d: &[Expr]) -> Option<(Vec<Expr>, Vec<Expr>)> {
    let dn = d.len().checked_sub(1)?;
    let lead = d.last()?.clone();
    if lead.is_zero() {
        return None;
    }
    let mut r: Vec<Expr> = n.to_vec();
    if r.len() <= dn {
        return Some((vec![Expr::zero()], r));
    }
    let mut q = vec![Expr::zero(); r.len() - dn];
    for i in (dn..r.len()).rev() {
        let c = expand(&div(&r[i], &lead));
        if c.is_zero() {
            continue;
        }
        q[i - dn] = c.clone();
        for j in 0..=dn {
            r[i - dn + j] = expand(&sub(&r[i - dn + j], &mul2(&c, &d[j])));
        }
    }
    r.truncate(dn);
    if r.is_empty() {
        r.push(Expr::zero());
    }
    Some((q, r))
}

fn trim(mut c: Vec<Expr>) -> Vec<Expr> {
    while c.len() > 1 && c.last().is_some_and(|x| x.is_zero()) {
        c.pop();
    }
    c
}

/// Split a product into (factors free of x, the rest).
fn split_const(f: &Expr, x: &Expr) -> (Expr, Expr) {
    match f.kind() {
        Kind::Mul(..) => {
            let (mut c, mut r) = (Vec::new(), Vec::new());
            for a in f.args() {
                if free(&a, x) {
                    c.push(a)
                } else {
                    r.push(a)
                }
            }
            (mul(&c), mul(&r))
        }
        _ => (Expr::one(), f.clone()),
    }
}

fn integ(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    if depth > MAX_DEPTH {
        return None;
    }
    if free(f, x) {
        return Some(mul2(f, x));
    }
    if f == x {
        return Some(div(&x.powi(2), &Expr::int(2)));
    }
    // guard against substitution loops
    if !seen.insert(f.clone()) {
        return None;
    }
    let r = integ_inner(f, x, depth, seen);
    seen.remove(f);
    r
}

fn integ_inner(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    // polynomials: integrate the expanded form term by term
    if let Some(c) = poly_coeffs(f, x) {
        let terms: Vec<Expr> = c
            .iter()
            .enumerate()
            .filter(|(_, k)| !k.is_zero())
            .map(|(i, k)| mul(&[k.clone(), x.powi(i as i64 + 1), Expr::rational(1, i as i64 + 1)]))
            .collect();
        return Some(add(&terms));
    }
    match f.kind() {
        Kind::Add(..) => {
            let mut v = Vec::new();
            for t in f.args() {
                v.push(integ(&t, x, depth + 1, seen)?);
            }
            return Some(add(&v));
        }
        Kind::Mul(..) => {
            let (c, rest) = split_const(f, x);
            if !c.is_one() {
                let r = integ(&rest, x, depth + 1, seen)?;
                // c * Piecewise((F/a, Ne(a, 0)), ...) with c a multiple of a:
                // the a = 0 branch is a constant, drop it (SymPy does too)
                if let Kind::Piecewise(p) = r.kind() {
                    if let Kind::Rel(crate::expr::RelOp::Ne, a, z) = p[0].1.kind() {
                        if z.is_zero() && !has(&div(&c, a), a) {
                            return Some(mul2(&c, &p[0].0));
                        }
                    }
                }
                return Some(mul2(&c, &r));
            }
        }
        _ => {}
    }
    if let Some(r) = table(f, x) {
        return Some(r);
    }
    if let Some(r) = rational(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = trig_powers(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = cyclic(f, x) {
        return Some(r);
    }
    if let Some(r) = usub(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = parts(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = special(f, x) {
        return Some(r);
    }
    // rewrite tan/cot/sec/csc products into sin/cos
    let rw = rewrite_trig(f);
    if rw != *f {
        if let Some(r) = integ(&rw, x, depth + 1, seen) {
            return Some(r);
        }
    }
    let ex = expand(f);
    if ex != *f {
        return integ(&ex, x, depth + 1, seen);
    }
    None
}

fn rewrite_trig(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Func(Func::Tan, a) => div(&sin(&a[0]), &cos(&a[0])),
        Kind::Func(Func::Cot, a) => div(&cos(&a[0]), &sin(&a[0])),
        Kind::Func(Func::Sec, a) => div(&Expr::one(), &cos(&a[0])),
        Kind::Func(Func::Csc, a) => div(&Expr::one(), &sin(&a[0])),
        _ if e.is_atom() => e.clone(),
        _ => crate::traverse::map_children(e, &mut rewrite_trig),
    }
}

/// Elementary antiderivatives of `F(a*x + b)`.
fn table(f: &Expr, x: &Expr) -> Option<Expr> {
    // power rule: (a x + b)^n
    if let Kind::Pow(base, n) = f.kind() {
        if free(n, x) {
            if let Some((a, _)) = linear(base, x) {
                if n.is_minus_one() {
                    return Some(div(&log(base), &a));
                }
                if n.is_number() {
                    let n1 = add2(n, &Expr::one());
                    return Some(div(&pow(base, &n1), &mul2(&a, &n1)));
                }
                // symbolic exponent: Piecewise like SymPy
                let n1 = add2(n, &Expr::one());
                return Some(piecewise(vec![
                    (div(&pow(base, &n1), &mul2(&a, &n1)), ne(n, &Expr::minus_one())),
                    (div(&log(base), &a), Expr::boolean(true)),
                ]));
            }
        }
        // c^(a x + b)
        if free(base, x) {
            if let Some((a, _)) = linear(n, x) {
                let lb = log(base);
                let r = div(f, &mul2(&a, &lb));
                if lb.is_number() || crate::functions::is_numeric_value(&lb) {
                    return Some(r);
                }
                return Some(piecewise(vec![(r, ne(&lb, &Expr::zero())), (x.clone(), Expr::boolean(true))]));
            }
        }
        // 1/sqrt(c - d x^2), 1/sqrt(x^2 + c), sqrt(c - d x^2), ...
        if let Some(r) = sqrt_quadratic(base, n, x) {
            return Some(r);
        }
        if let Some(r) = inverse_quadratic(base, n, x) {
            return Some(r);
        }
        // sec^2, csc^2, 1/cos^2, 1/sin^2 (SymPy writes tan as sin/cos)
        if let (Kind::Func(g, a), Some(k)) = (base.kind(), n.as_i64()) {
            if let Some((da, _)) = linear(&a[0], x) {
                let u = &a[0];
                let r = match (g, k) {
                    (Func::Sec, 2) | (Func::Cos, -2) => Some(div(&sin(u), &cos(u))),
                    (Func::Csc, 2) | (Func::Sin, -2) => Some(neg(&div(&cos(u), &sin(u)))),
                    (Func::Tan, 2) => Some(sub(&div(&sin(u), &cos(u)), u)),
                    (Func::Cot, 2) => Some(sub(&neg(&div(&cos(u), &sin(u))), u)),
                    _ => None,
                };
                if let Some(r) = r {
                    return Some(div(&r, &da));
                }
            }
        }
    }
    let Kind::Func(func, args) = f.kind() else {
        return None;
    };
    if args.len() != 1 {
        return None;
    }
    let u = &args[0];
    let (a, _) = linear(u, x)?;
    let sym = numeric_or_piecewise_guard(&a);
    let res = match func {
        Func::Exp => f.clone(),
        Func::Sin => neg(&cos(u)),
        Func::Cos => sin(u),
        Func::Tan => neg(&log(&cos(u))),
        Func::Cot => log(&sin(u)),
        Func::Sec => log(&add2(&tan(u), &sec(u))),
        Func::Csc => neg(&log(&add2(&cot(u), &csc(u)))),
        Func::Sinh => cosh(u),
        Func::Cosh => sinh(u),
        Func::Tanh => sub(u, &log(&add2(&tanh(u), &Expr::one()))),
        Func::Log => sub(&mul2(u, &log(u)), u),
        Func::Atan => sub(&mul2(u, &atan(u)), &div(&log(&add2(&u.powi(2), &Expr::one())), &Expr::int(2))),
        Func::Asin => add2(&mul2(u, &asin(u)), &sqrt(&sub(&Expr::one(), &u.powi(2)))),
        Func::Acos => sub(&mul2(u, &acos(u)), &sqrt(&sub(&Expr::one(), &u.powi(2)))),
        Func::Asinh => sub(&mul2(u, &asinh(u)), &sqrt(&add2(&u.powi(2), &Expr::one()))),
        Func::Atanh => add2(&mul2(u, &atanh(u)), &div(&log(&sub(&u.powi(2), &Expr::one())), &Expr::int(2))),
        Func::Erf => add2(
            &mul2(u, &erf(u)),
            &div(&exp(&neg(&u.powi(2))), &sqrt(&Expr::pi())),
        ),
        _ => return None,
    };
    let r = div(&res, &a);
    if sym {
        // exp(a*x) with symbolic a: SymPy returns a Piecewise
        let at_zero = mul2(&subs1(f, &a, &Expr::zero()), x);
        return Some(piecewise(vec![(r, ne(&a, &Expr::zero())), (at_zero, Expr::boolean(true))]));
    }
    Some(r)
}

fn numeric_or_piecewise_guard(a: &Expr) -> bool {
    !crate::functions::is_numeric_value(a)
}

/// Integrals involving sqrt of a quadratic `c + d x^2` (no linear term).
fn sqrt_quadratic(base: &Expr, n: &Expr, x: &Expr) -> Option<Expr> {
    let half = Expr::rational(1, 2);
    let mhalf = Expr::rational(-1, 2);
    if *n != half && *n != mhalf {
        return None;
    }
    let c = poly_coeffs(base, x)?;
    if c.len() != 3 || !c[1].is_zero() {
        return None;
    }
    let (c0, d) = (&c[0], &c[2]);
    let (c0v, dv) = (crate::evalf::eval_f64(c0)?, crate::evalf::eval_f64(d)?);
    if *n == mhalf {
        if c0v > 0.0 && dv < 0.0 {
            // 1/sqrt(c - k x^2) = asin(sqrt(k/c) x)/sqrt(k)
            let k = neg(d);
            return Some(div(&asin(&mul2(&sqrt(&div(&k, c0)), x)), &sqrt(&k)));
        }
        if c0v > 0.0 && dv > 0.0 {
            return Some(div(&asinh(&mul2(&sqrt(&div(d, c0)), x)), &sqrt(d)));
        }
        if c0v < 0.0 && dv > 0.0 {
            return Some(div(&acosh(&mul2(&sqrt(&div(d, &neg(c0))), x)), &sqrt(d)));
        }
        return None;
    }
    // sqrt(c + d x^2)
    let s = sqrt(base);
    if c0v > 0.0 && dv < 0.0 {
        let k = neg(d);
        let t = mul2(&sqrt(&div(&k, c0)), x);
        return Some(add2(
            &div(&mul2(x, &s), &Expr::int(2)),
            &div(&mul2(c0, &asin(&t)), &mul2(&Expr::int(2), &sqrt(&k))),
        ));
    }
    if c0v > 0.0 && dv > 0.0 {
        let t = mul2(&sqrt(&div(d, c0)), x);
        return Some(add2(
            &div(&mul2(x, &s), &Expr::int(2)),
            &div(&mul2(c0, &asinh(&t)), &mul2(&Expr::int(2), &sqrt(d))),
        ));
    }
    None
}

/// 1/(c + d x^2) -> atan / atanh forms.
fn inverse_quadratic(base: &Expr, n: &Expr, x: &Expr) -> Option<Expr> {
    if !n.is_minus_one() {
        return None;
    }
    let c = poly_coeffs(base, x)?;
    if c.len() != 3 || !c[1].is_zero() {
        return None;
    }
    let (c0, d) = (&c[0], &c[2]);
    let (c0v, dv) = (crate::evalf::eval_f64(c0)?, crate::evalf::eval_f64(d)?);
    if c0v * dv > 0.0 {
        // atan(sqrt(d/c) x) / sqrt(c d)
        let t = mul2(&sqrt(&div(d, c0)), x);
        let r = div(&atan(&t), &sqrt(&mul2(c0, d)));
        return Some(if c0v < 0.0 { neg(&r) } else { r });
    }
    None
}

/// Rational functions with a denominator of degree <= 2, or a power of a
/// linear polynomial.
fn rational(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let (num, den) = numer_denom(f);
    if free(&den, x) {
        return None;
    }
    let nc = poly_coeffs(&num, x)?;
    // denominator a power of a linear polynomial: substitute u = a x + b
    if let Kind::Pow(b, k) = den.kind() {
        if let (Some((a, c)), Some(_)) = (linear(b, x), k.as_i64()) {
            let u = Expr::symbol("_u");
            let xu = div(&sub(&u, &c), &a);
            let numu = expand(&subs1(&num, x, &xu));
            let g = expand(&div(&numu, &pow(&u, k)));
            let gi = integ(&g, &u, depth + 1, seen)?;
            return Some(div(&subs1(&gi, &u, b), &a));
        }
    }
    let dc = trim(poly_coeffs(&den, x)?);
    let nc = trim(nc);
    if dc.len() < 2 {
        return None;
    }
    let (q, r) = poly_divmod(&nc, &dc)?;
    let mut parts = Vec::new();
    let qp = from_coeffs(&q, x);
    if !qp.is_zero() {
        parts.push(integ(&qp, x, depth + 1, seen)?);
    }
    let r = trim(r);
    if r.len() == 1 && r[0].is_zero() {
        return Some(add(&parts));
    }
    match dc.len() {
        2 => {
            // r0 / (d1 x + d0)
            parts.push(div(&mul2(&r[0], &log(&from_coeffs(&dc, x))), &dc[1]));
        }
        3 => {
            let (c, b, a) = (&dc[0], &dc[1], &dc[2]);
            let p = r.get(1).cloned().unwrap_or_else(Expr::zero);
            let qq = r[0].clone();
            let disc = expand(&sub(&b.powi(2), &mul(&[Expr::int(4), a.clone(), c.clone()])));
            let dv = crate::evalf::eval_f64(&disc)?;
            let d_expr = from_coeffs(&dc, x);
            if dv.abs() < 1e-300 {
                return None;
            }
            if dv < 0.0 {
                // (p x + q)/(a x^2 + b x + c)
                //  = p/(2a) * log(den) + (q - p b/(2a)) * 2/sqrt(-disc) * atan((2 a x + b)/sqrt(-disc))
                let s = sqrt(&neg(&disc));
                let two_a = mul2(&Expr::int(2), a);
                let log_part = mul2(&div(&p, &two_a), &log(&d_expr));
                let k = sub(&qq, &div(&mul2(&p, b), &two_a));
                let atan_part = mul(&[
                    k,
                    Expr::int(2),
                    pow(&s, &Expr::minus_one()),
                    atan(&div(&add2(&mul2(&two_a, x), b), &s)),
                ]);
                if !p.is_zero() {
                    parts.push(log_part);
                }
                if !atan_part.is_zero() {
                    parts.push(expand(&atan_part));
                }
            } else {
                // two real roots r1, r2: A/(x - r1) + B/(x - r2)
                let s = sqrt(&disc);
                let two_a = mul2(&Expr::int(2), a);
                let r1 = div(&sub(&neg(b), &s), &two_a);
                let r2 = div(&add2(&neg(b), &s), &two_a);
                let numer = from_coeffs(&r, x);
                // A = num(r1) / (a (r1 - r2)), B = num(r2) / (a (r2 - r1))
                let a1 = expand(&div(&subs1(&numer, x, &r1), &mul2(a, &sub(&r1, &r2))));
                let b1 = expand(&div(&subs1(&numer, x, &r2), &mul2(a, &sub(&r2, &r1))));
                parts.push(mul2(&a1, &log(&sub(x, &r1))));
                parts.push(mul2(&b1, &log(&sub(x, &r2))));
            }
        }
        _ => {
            // a polynomial denominator of higher degree: needs full
            // factorisation (provided by the poly module once available)
            return None;
        }
    }
    Some(add(&parts))
}

/// sin(u)^m cos(u)^n with u linear.
fn trig_powers(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let facs: Vec<(Expr, Expr)> = match f.kind() {
        Kind::Mul(c, facs) if c.is_one() => facs.clone(),
        Kind::Pow(b, e) => vec![(b.clone(), e.clone())],
        _ => return None,
    };
    let mut arg: Option<Expr> = None;
    let (mut m, mut n) = (0i64, 0i64);
    for (b, e) in &facs {
        let k = e.as_i64()?;
        match b.kind() {
            Kind::Func(Func::Sin, a) | Kind::Func(Func::Cos, a) => {
                if arg.as_ref().is_some_and(|u| *u != a[0]) {
                    return None;
                }
                arg = Some(a[0].clone());
                if b.is_func(&Func::Sin) {
                    m += k
                } else {
                    n += k
                }
            }
            _ => return None,
        }
    }
    let u = arg?;
    let (a, _) = linear(&u, x)?;
    if m < 0 || n < 0 {
        return None;
    }
    let s = sin(&u);
    let c = cos(&u);
    let t = Expr::symbol("_t");
    if n % 2 == 1 {
        let g = mul(&[pow(&sub(&Expr::one(), &t.powi(2)), &Expr::int((n - 1) / 2)), t.powi(m)]);
        let gi = integ(&expand(&g), &t, depth + 1, seen)?;
        return Some(div(&subs1(&gi, &t, &s), &a));
    }
    if m % 2 == 1 {
        // sin^m cos^n = -(1 - t^2)^((m-1)/2) t^n dt, t = cos u
        let g = mul(&[
            Expr::minus_one(),
            pow(&sub(&Expr::one(), &t.powi(2)), &Expr::int((m - 1) / 2)),
            t.powi(n),
        ]);
        let gi = integ(&expand(&g), &t, depth + 1, seen)?;
        return Some(div(&subs1(&gi, &t, &c), &a));
    }
    // both even: reduction formula on the larger power
    if m >= 2 {
        // ∫ sin^m cos^n = -sin^(m-1) cos^(n+1)/(m+n) + (m-1)/(m+n) ∫ sin^(m-2) cos^n
        let mn = Expr::int(m + n);
        let first = neg(&div(&mul(&[s.powi(m - 1), c.powi(n + 1)]), &mn));
        let rest_f = mul(&[s.powi(m - 2), c.powi(n)]);
        let rest = integ(&rest_f, x, depth + 1, seen)?;
        return Some(add2(&div(&first, &a), &mul2(&Expr::rational(m - 1, m + n), &rest)));
    }
    if n >= 2 {
        let mn = Expr::int(m + n);
        let first = div(&mul(&[s.powi(m + 1), c.powi(n - 1)]), &mn);
        let rest_f = mul(&[s.powi(m), c.powi(n - 2)]);
        let rest = integ(&rest_f, x, depth + 1, seen)?;
        return Some(add2(&div(&first, &a), &mul2(&Expr::rational(n - 1, m + n), &rest)));
    }
    None
}

/// exp(a x) * sin(b x) and exp(a x) * cos(b x).
fn cyclic(f: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Mul(c, facs) = f.kind() else { return None };
    if !c.is_one() || facs.len() != 2 {
        return None;
    }
    let mut e_arg = None;
    let mut t: Option<(bool, Expr)> = None;
    for (b, e) in facs {
        if !e.is_one() {
            return None;
        }
        match b.kind() {
            Kind::Func(Func::Exp, a) => e_arg = Some(a[0].clone()),
            Kind::Func(Func::Sin, a) => t = Some((true, a[0].clone())),
            Kind::Func(Func::Cos, a) => t = Some((false, a[0].clone())),
            _ => return None,
        }
    }
    let (ea, (is_sin, ta)) = (e_arg?, t?);
    let (a, _) = linear(&ea, x)?;
    let (b, _) = linear(&ta, x)?;
    let den = add2(&a.powi(2), &b.powi(2));
    let ex = exp(&ea);
    let (s, co) = (sin(&ta), cos(&ta));
    let inner = if is_sin {
        sub(&mul2(&a, &s), &mul2(&b, &co))
    } else {
        add2(&mul2(&a, &co), &mul2(&b, &s))
    };
    Some(expand(&div(&mul2(&ex, &inner), &den)))
}

/// Candidate substitutions: function arguments, power bases and exponents,
/// and the functions themselves.
fn candidates(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    if free(e, x) || e == x {
        return;
    }
    match e.kind() {
        Kind::Func(_, args) => {
            out.push(e.clone());
            for a in args {
                if a != x {
                    out.push(a.clone());
                }
                candidates(a, x, out);
            }
        }
        Kind::Pow(b, n) => {
            if b != x {
                out.push(b.clone());
            }
            if !free(n, x) {
                out.push(n.clone());
            }
            if n.is_rational_number() && !n.is_integer_number() {
                out.push(e.clone());
            }
            candidates(b, x, out);
            candidates(n, x, out);
        }
        _ => {
            for a in e.args() {
                candidates(&a, x, out);
            }
        }
    }
}

fn usub(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let mut cands = Vec::new();
    candidates(f, x, &mut cands);
    let mut uniq = FxHashSet::default();
    cands.retain(|c| uniq.insert(c.clone()));
    // prefer bigger sub-expressions first, then SymPy's canonical order
    cands.sort_by_cached_key(|c| (std::cmp::Reverse(crate::traverse::node_count(c)), crate::sortkey::sort_key(c)));
    let usym = Expr::symbol("_u");
    for u in cands {
        if linear(&u, x).is_some() && !matches!(u.kind(), Kind::Func(..)) {
            continue;
        }
        let du = diff(&u, x);
        if du.is_zero() {
            continue;
        }
        let g = div(f, &du);
        let gu = subs1(&g, &u, &usym);
        let gu = if free(&gu, x) {
            gu
        } else {
            // try harder: cancel common factors by expanding the ratio
            let alt = subs1(&expand(&g), &u, &usym);
            if !free(&alt, x) {
                continue;
            }
            alt
        };
        if let Some(r) = integ(&gu, &usym, depth + 1, seen) {
            return Some(subs1(&r, &usym, &u));
        }
    }
    None
}

fn is_poly_factor(e: &Expr, x: &Expr) -> bool {
    poly_coeffs(e, x).is_some() && !free(e, x)
}

/// Integration by parts, choosing `u` by the LIATE rule.
fn parts(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let factors: Vec<Expr> = match f.kind() {
        Kind::Mul(..) => f.args(),
        Kind::Func(..) => vec![f.clone()],
        Kind::Pow(b, n) if n.is_integer_number() && matches!(b.kind(), Kind::Func(Func::Log, _)) => {
            vec![f.clone()]
        }
        _ => return None,
    };
    // rank: Log, Inverse trig, Algebraic, Trig, Exp
    fn rank(e: &Expr, x: &Expr) -> u8 {
        match e.kind() {
            Kind::Func(Func::Log, _) => 0,
            Kind::Pow(b, _) if matches!(b.kind(), Kind::Func(Func::Log, _)) => 0,
            Kind::Func(Func::Atan | Func::Asin | Func::Acos | Func::Asinh | Func::Atanh | Func::Acot, _) => 1,
            _ if is_poly_factor(e, x) => 2,
            Kind::Func(Func::Sin | Func::Cos | Func::Sinh | Func::Cosh, _) => 3,
            Kind::Func(Func::Exp, _) => 4,
            Kind::Pow(b, _) if !free(&diff(e, x), x) && free(b, x) => 4,
            _ => 5,
        }
    }
    let mut idx: Vec<usize> = (0..factors.len()).collect();
    idx.sort_by_key(|&i| rank(&factors[i], x));
    for &i in idx.iter().take(2) {
        let u = &factors[i];
        if rank(u, x) > 2 {
            continue;
        }
        let dv: Vec<Expr> = factors.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, e)| e.clone()).collect();
        let dv = mul(&dv);
        if rank(u, x) == 2 && dv.is_one() {
            continue;
        }
        let v = integ(&dv, x, depth + 1, seen)?;
        let du = diff(u, x);
        let rest = expand(&mul2(&v, &du));
        let ri = integ(&rest, x, depth + 1, seen)?;
        let r = sub(&mul2(u, &v), &ri);
        // polynomial * exp(..): SymPy factors the exponential out
        if let Kind::Func(Func::Exp, _) = dv.kind() {
            if rank(u, x) == 2 {
                let q = expand(&div(&r, &dv));
                if poly_coeffs(&q, x).is_some() {
                    return Some(mul2(&q, &dv));
                }
            }
        }
        return Some(r);
    }
    None
}

/// Non-elementary special cases.
fn special(f: &Expr, x: &Expr) -> Option<Expr> {
    // exp(-c x^2) -> sqrt(pi) erf(sqrt(c) x) / (2 sqrt(c))
    if let Kind::Func(Func::Exp, a) = f.kind() {
        let c = poly_coeffs(&a[0], x)?;
        if c.len() == 3 && c[0].is_zero() && c[1].is_zero() {
            let k = neg(&c[2]);
            if crate::evalf::eval_f64(&k).is_some_and(|v| v > 0.0) {
                let sk = sqrt(&k);
                return Some(div(&mul2(&sqrt(&Expr::pi()), &erf(&mul2(&sk, x))), &mul2(&Expr::int(2), &sk)));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn check(f: &str) {
        let x = Expr::symbol("x");
        let fe = parse_str(f).unwrap();
        let r = try_integrate(&fe, &x).unwrap_or_else(|| panic!("failed to integrate {f}"));
        // verify numerically that d/dx r == f at a few points
        let d = diff(&r, &x);
        for v in [0.3, 0.55, 0.7] {
            let a = eval_c64(&subs(&d, &[(x.clone(), Expr::float(v))]));
            let b = eval_c64(&subs(&fe, &[(x.clone(), Expr::float(v))]));
            match (a, b) {
                (Some(a), Some(b)) => assert!((a - b).norm() < 1e-8 * (1.0 + b.norm()), "{f}: got {r}"),
                _ => panic!("{f}: cannot evaluate {r}"),
            }
        }
    }

    #[test]
    fn basic() {
        for f in [
            "x**3 + 2*x", "1/x", "exp(2*x)", "sin(3*x+1)", "x*exp(x)", "x**2*sin(x)", "log(x)",
            "x*log(x)", "atan(x)", "1/(1+x**2)", "1/sqrt(1-x**2)", "sin(x)**2", "cos(x)**3",
            "sin(x)*cos(x)", "exp(x)*sin(x)", "2*x*exp(x**2)", "x/(x**2+1)", "cos(x)*exp(sin(x))",
            "1/(x**2-1)", "1/(x**2+2*x+5)", "(x+1)/(x**2+x+1)", "tan(x)", "log(x)/x", "1/(x*log(x))",
            "sqrt(1-x**2)", "asin(x)", "x**3*log(x)", "log(x)**2", "1/(2*x+3)**3", "x/(x+1)**2",
            "exp(x)/(1+exp(x))", "sec(x)**2", "sin(x)**3*cos(x)**2", "exp(-x**2)", "x*atan(x)",
        ] {
            check(f);
        }
    }
}
