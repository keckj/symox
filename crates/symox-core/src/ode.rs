//! Ordinary differential equations: `dsolve`.
//!
//! Supported classes: linear equations with constant coefficients
//! (homogeneous part from the characteristic roots, particular solutions
//! by undetermined coefficients), Cauchy–Euler equations, first-order
//! linear equations (integrating factor), separable and Bernoulli
//! first-order equations. Initial conditions solve for the constants.

use rustc_hash::FxHashMap;

use crate::diff::diff;
use crate::expand::expand;
use crate::expr::{Expr, Func, Kind};
use crate::functions::*;
use crate::integrate::{poly_coeffs, try_integrate};
use crate::ops::{add, div, mul, mul2, neg, pow, sub};
use crate::relational::eq;
use crate::simplify::together;
use crate::solve::{as_expr, poly_roots, solve, solve_system, SystemSolution};
use crate::subs::{subs, subs1, xreplace};
use crate::traverse::{has, map_children};

fn constant(i: usize) -> Expr {
    Expr::symbol(&format!("C{i}"))
}

/// Replace derivatives of `f` by symbols `_y0, _y1, ...`; returns the
/// rewritten expression and the order.
fn to_symbols(e: &Expr, f: &Expr, x: &Expr, order: &mut usize) -> Expr {
    if e == f {
        return Expr::symbol("_y0");
    }
    if let Kind::Derivative(g, vars) = e.kind() {
        if g == f && vars.len() == 1 && vars[0].0 == *x {
            let k = vars[0].1 as usize;
            *order = (*order).max(k);
            return Expr::symbol(&format!("_y{k}"));
        }
    }
    if e.is_atom() {
        return e.clone();
    }
    map_children(e, &mut |c| to_symbols(c, f, x, order))
}

fn ysym(k: usize) -> Expr {
    Expr::symbol(&format!("_y{k}"))
}

/// Solve the ODE `ode = 0` (or an `Eq`) for the function `f = f(x)`.
/// Returns the list of solutions `Eq(f(x), ...)`.
pub fn dsolve(ode: &Expr, f: &Expr, ics: &[(Expr, Expr)]) -> Option<Vec<Expr>> {
    let Kind::Func(Func::Undefined(_), args) = f.kind() else { return None };
    if args.len() != 1 {
        return None;
    }
    let x = args[0].clone();
    let mut order = 0;
    let e = expand(&to_symbols(&as_expr(ode), f, &x, &mut order));
    let e = together(&e);
    let (e, _) = crate::simplify::numer_denom(&e);
    let e = expand(&e);
    if order == 0 {
        return None;
    }
    let sols = if let Some(lin) = linear_coeffs(&e, order, &x) {
        vec![solve_linear(&lin, order, &x)?]
    } else if order == 1 {
        solve_first_order_nonlinear(&e, &x)?
    } else {
        return None;
    };
    let mut out = Vec::new();
    for s in sols {
        let s = if ics.is_empty() { s } else { apply_ics(&s, f, &x, ics)? };
        out.push(eq(f, &tidy(&s, &x)));
    }
    Some(out)
}

/// Coefficients `a_k(x)` and right-hand side `g(x)` of a linear ODE
/// `sum a_k y^(k) = g`.
struct Linear {
    a: Vec<Expr>,
    g: Expr,
}

fn linear_coeffs(e: &Expr, order: usize, x: &Expr) -> Option<Linear> {
    let mut a = Vec::with_capacity(order + 1);
    let mut rest = e.clone();
    for k in 0..=order {
        let y = ysym(k);
        let c = poly_coeffs(e, &y)?;
        if c.len() > 2 {
            return None;
        }
        let ak = c.get(1).cloned().unwrap_or_else(Expr::zero);
        if (0..=order).any(|j| has(&ak, &ysym(j))) {
            return None;
        }
        a.push(ak.clone());
        rest = expand(&sub(&rest, &mul2(&ak, &y)));
    }
    if (0..=order).any(|j| has(&rest, &ysym(j))) {
        return None;
    }
    let _ = x;
    Some(Linear { a, g: neg(&rest) })
}

fn solve_linear(l: &Linear, order: usize, x: &Expr) -> Option<Expr> {
    let lead = l.a[order].clone();
    let consts = l.a.iter().all(|c| !has(c, x));
    if consts {
        let (hom, roots) = constant_coeff_homogeneous(&l.a, x)?;
        let part = if l.g.is_zero() { Expr::zero() } else { undetermined(&l.a, &roots, &l.g, x).or_else(|| variation_of_parameters(&l.a, &hom, &l.g, x))? };
        let hom_sum = group_by_exp(&hom.iter().enumerate().map(|(i, b)| mul2(&constant(i + 1), b)).collect::<Vec<_>>());
        return Some(add(&[hom_sum, part]));
    }
    // Cauchy–Euler: a_k = c_k x^k
    if l.g.is_zero() {
        if let Some(r) = cauchy_euler(&l.a, x) {
            return Some(r);
        }
    }
    if order == 1 {
        // y' + p y = q
        let p = together(&div(&l.a[0], &lead));
        let q = together(&div(&l.g, &lead));
        let ip = try_integrate(&p, x)?;
        let mu = crate::simp::powsimp(&exp(&ip));
        let mu = simplify_exp_log(&mu);
        let inner = if q.is_zero() { Expr::zero() } else { try_integrate(&expand(&mul2(&mu, &q)), x)? };
        let y = div(&add(&[constant(1), inner]), &mu);
        return Some(expand(&y));
    }
    None
}

/// Sum the terms, factoring out an exponential shared by several of them:
/// C1*exp(-x) + C2*x*exp(-x) -> (C1 + C2*x)*exp(-x).
fn group_by_exp(terms: &[Expr]) -> Expr {
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    for t in terms {
        let factors = match t.kind() {
            Kind::Mul(..) => t.args(),
            _ => vec![t.clone()],
        };
        let (e, rest): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(|f| f.is_func(&Func::Exp));
        let key = mul(&e);
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, v)) => v.push(mul(&rest)),
            None => groups.push((key, vec![mul(&rest)])),
        }
    }
    add(&groups.into_iter().map(|(k, v)| mul2(&add(&v), &k)).collect::<Vec<_>>())
}

/// exp(-log(x)) -> 1/x and friends.
fn simplify_exp_log(e: &Expr) -> Expr {
    let e = crate::simp::expand_log(&crate::simp::logcombine(e, true), false);
    match e.kind() {
        Kind::Func(Func::Exp, a) => {
            let l = crate::simp::logcombine(&a[0], true);
            exp(&l)
        }
        _ => e,
    }
}

/// Characteristic roots with multiplicities.
fn char_roots(a: &[Expr]) -> Option<Vec<(Expr, usize)>> {
    let roots = poly_roots(a)?;
    let r = Expr::symbol("_r");
    let p = add(&a.iter().enumerate().map(|(i, c)| mul2(c, &r.powi(i as i64))).collect::<Vec<_>>());
    let mut out = Vec::new();
    for root in roots {
        // multiplicity: number of derivatives vanishing at the root
        let mut m = 0;
        let mut d = p.clone();
        loop {
            let v = expand(&subs1(&d, &r, &root));
            let zero = v.is_zero() || crate::evalf::eval_c64(&v).is_some_and(|z| z.norm() < 1e-12);
            if !zero {
                break;
            }
            m += 1;
            d = diff(&d, &r);
        }
        out.push((root, m.max(1)));
    }
    Some(out)
}

fn constant_coeff_homogeneous(a: &[Expr], x: &Expr) -> Option<(Vec<Expr>, Vec<(Expr, usize)>)> {
    let mut roots = char_roots(a)?;
    roots.sort_by_cached_key(|(r, _)| crate::sortkey::sort_key(r));
    let mut basis = Vec::new();
    let mut done: Vec<Expr> = Vec::new();
    for (r, m) in &roots {
        if done.contains(r) {
            continue;
        }
        let numeric = crate::functions::is_numeric_value(r);
        let re_ = if numeric { expand(&re(r)) } else { r.clone() };
        let im_ = if numeric { expand(&im(r)) } else { Expr::zero() };
        let is_real = im_.is_zero();
        if !is_real {
            // pair a +- b I: use the one with b > 0
            let conj = expand(&sub(&re_, &mul2(&Expr::i(), &im_)));
            done.push(conj);
            done.push(r.clone());
            let b = if crate::evalf::eval_f64(&im_).is_some_and(|v| v < 0.0) { neg(&im_) } else { im_.clone() };
            for k in 0..*m {
                let xk = x.powi(k as i64);
                let ex = exp(&mul2(&re_, x));
                basis.push(mul(&[xk.clone(), ex.clone(), sin(&mul2(&b, x))]));
                basis.push(mul(&[xk, ex, cos(&mul2(&b, x))]));
            }
        } else {
            done.push(r.clone());
            for k in 0..*m {
                basis.push(mul2(&x.powi(k as i64), &exp(&mul2(r, x))));
            }
        }
    }
    Some((basis, roots))
}

/// Apply the operator sum a_k D^k to `y`.
fn apply_op(a: &[Expr], y: &Expr, x: &Expr) -> Expr {
    let mut acc = Vec::with_capacity(a.len());
    let mut d = y.clone();
    for (k, c) in a.iter().enumerate() {
        if k > 0 {
            d = diff(&d, x);
        }
        if !c.is_zero() {
            acc.push(mul2(c, &d));
        }
    }
    expand(&add(&acc))
}

/// Split a forcing term into `x^k * exp(alpha x) * (sin|cos)(beta x)` data.
fn forcing_terms(g: &Expr, x: &Expr) -> Option<Vec<(usize, Expr, Expr)>> {
    let terms = match g.kind() {
        Kind::Add(..) => g.args(),
        _ => vec![g.clone()],
    };
    let mut out: Vec<(usize, Expr, Expr)> = Vec::new();
    for t in terms {
        let factors = match t.kind() {
            Kind::Mul(..) => t.args(),
            _ => vec![t.clone()],
        };
        let (mut deg, mut alpha, mut beta) = (0usize, Expr::zero(), Expr::zero());
        for fct in factors {
            if !has(&fct, x) {
                continue;
            }
            if fct == *x {
                deg += 1;
                continue;
            }
            match fct.kind() {
                Kind::Pow(b, e) if b == x => deg += usize::try_from(e.as_i64()?).ok()?,
                Kind::Func(Func::Exp, a) => {
                    let k = diff(&a[0], x);
                    if has(&k, x) {
                        return None;
                    }
                    alpha = add(&[alpha, k]);
                }
                Kind::Func(Func::Sin | Func::Cos, a) => {
                    let k = diff(&a[0], x);
                    if has(&k, x) || !beta.is_zero() {
                        return None;
                    }
                    beta = k;
                }
                _ => return None,
            }
        }
        // merge with an existing family of the same (alpha, beta)
        match out.iter_mut().find(|(_, a, b)| *a == alpha && *b == beta) {
            Some(entry) => entry.0 = entry.0.max(deg),
            None => out.push((deg, alpha, beta)),
        }
    }
    Some(out)
}

fn undetermined(a: &[Expr], roots: &[(Expr, usize)], g: &Expr, x: &Expr) -> Option<Expr> {
    let fams = forcing_terms(g, x)?;
    let mut unknowns = Vec::new();
    let mut trial_terms = Vec::new();
    let mut n = 0;
    for (deg, alpha, beta) in fams {
        let target = expand(&add(&[alpha.clone(), mul2(&Expr::i(), &beta)]));
        let s = roots.iter().find(|(r, _)| expand(&sub(r, &target)).is_zero()).map(|(_, m)| *m).unwrap_or(0);
        let ex = exp(&mul2(&alpha, x));
        for k in 0..=deg {
            let xk = x.powi((k + s) as i64);
            let mut push = |f: Expr| {
                let c = Expr::symbol(&format!("_A{n}"));
                n += 1;
                unknowns.push(c.clone());
                trial_terms.push(mul(&[c, xk.clone(), ex.clone(), f]));
            };
            if beta.is_zero() {
                push(Expr::one());
            } else {
                push(sin(&mul2(&beta, x)));
                push(cos(&mul2(&beta, x)));
            }
        }
    }
    let trial = add(&trial_terms);
    let resid = expand(&sub(&apply_op(a, &trial, x), g));
    // collect coefficients of the independent x-dependent parts
    let mut groups: FxHashMap<Expr, Vec<Expr>> = FxHashMap::default();
    let terms = match resid.kind() {
        Kind::Add(..) => resid.args(),
        _ => vec![resid.clone()],
    };
    for t in terms {
        let factors = match t.kind() {
            Kind::Mul(..) => t.args(),
            _ => vec![t.clone()],
        };
        let (dep, ind): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(|f| has(f, x));
        groups.entry(mul(&dep)).or_default().push(mul(&ind));
    }
    let eqs: Vec<Expr> = groups.into_values().map(|v| add(&v)).collect();
    let sol = solve_system(&eqs, &unknowns)?;
    let SystemSolution::Map(m) = sol else { return None };
    let map: FxHashMap<Expr, Expr> = m.into_iter().collect();
    let mut yp = xreplace(&trial, &map);
    for u in &unknowns {
        yp = subs1(&yp, u, &Expr::zero());
    }
    Some(expand(&yp))
}

fn variation_of_parameters(a: &[Expr], hom: &[Expr], g: &Expr, x: &Expr) -> Option<Expr> {
    // second order only: y_p = -y1 ∫ y2 g / W + y2 ∫ y1 g / W  (normalised)
    if hom.len() != 2 {
        return None;
    }
    let lead = a.last()?;
    let gn = div(g, lead);
    let (y1, y2) = (&hom[0], &hom[1]);
    let w = expand(&sub(&mul2(y1, &diff(y2, x)), &mul2(y2, &diff(y1, x))));
    let w = crate::simp::trigsimp(&w);
    let i1 = try_integrate(&expand(&div(&mul2(y2, &gn), &w)), x)?;
    let i2 = try_integrate(&expand(&div(&mul2(y1, &gn), &w)), x)?;
    Some(expand(&add(&[neg(&mul2(y1, &i1)), mul2(y2, &i2)])))
}

fn cauchy_euler(a: &[Expr], x: &Expr) -> Option<Expr> {
    // a_k = c_k x^k
    let mut c = Vec::with_capacity(a.len());
    for (k, ak) in a.iter().enumerate() {
        let ck = expand(&div(ak, &x.powi(k as i64)));
        if has(&ck, x) {
            return None;
        }
        c.push(ck);
    }
    // characteristic polynomial sum c_k r(r-1)...(r-k+1)
    let r = Expr::symbol("_r");
    let mut p = Vec::new();
    for (k, ck) in c.iter().enumerate() {
        let mut fall = Expr::one();
        for j in 0..k {
            fall = mul2(&fall, &sub(&r, &Expr::int(j as i64)));
        }
        p.push(mul2(ck, &fall));
    }
    let p = expand(&add(&p));
    let coeffs = poly_coeffs(&p, &r)?;
    let roots = char_roots(&coeffs)?;
    let mut basis = Vec::new();
    let mut done: Vec<Expr> = Vec::new();
    for (rt, m) in roots {
        if done.contains(&rt) {
            continue;
        }
        done.push(rt.clone());
        let im_ = expand(&im(&rt));
        if !im_.is_zero() {
            let re_ = expand(&re(&rt));
            done.push(expand(&sub(&re_, &mul2(&Expr::i(), &im_))));
            let b = if crate::evalf::eval_f64(&im_).is_some_and(|v| v < 0.0) { neg(&im_) } else { im_ };
            for k in 0..m {
                let lk = log(x).powi(k as i64);
                basis.push(mul(&[lk.clone(), pow(x, &re_), sin(&mul2(&b, &log(x)))]));
                basis.push(mul(&[lk, pow(x, &re_), cos(&mul2(&b, &log(x)))]));
            }
        } else {
            for k in 0..m {
                basis.push(mul2(&log(x).powi(k as i64), &pow(x, &rt)));
            }
        }
    }
    Some(add(&basis.iter().enumerate().map(|(i, b)| mul2(&constant(i + 1), b)).collect::<Vec<_>>()))
}

/// Separable (and Bernoulli) first-order equations.
fn solve_first_order_nonlinear(e: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    let y = ysym(0);
    let yp = ysym(1);
    // y' = F(x, y)
    let fs = solve(e, &yp)?;
    if fs.len() != 1 {
        return None;
    }
    let big_f = expand(&fs[0]);
    // separable: F = g(x) * h(y)
    let factors = match big_f.kind() {
        Kind::Mul(..) => big_f.args(),
        _ => vec![big_f.clone()],
    };
    let (gx, hy): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(|f| !has(f, &y));
    let (gx, hy) = (mul(&gx), mul(&hy));
    if has(&hy, x) {
        // Bernoulli: y' = q y^n - p y
        return bernoulli(&big_f, x);
    }
    let lhs = try_integrate(&div(&Expr::one(), &hy), &y)?;
    let rhs = add(&[try_integrate(&gx, x)?, constant(1)]);
    let sols = solve(&sub(&lhs, &rhs), &y)?;
    if sols.is_empty() {
        return None;
    }
    Some(sols)
}

fn bernoulli(f: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    let y = ysym(0);
    let terms = match f.kind() {
        Kind::Add(..) => f.args(),
        _ => return None,
    };
    let (mut p, mut q, mut n) = (Expr::zero(), Expr::zero(), None);
    for t in terms {
        let c = div(&t, &y);
        if !has(&c, &y) {
            p = sub(&p, &c);
            continue;
        }
        // t = q(x) y^n
        let factors = t.args();
        let (dep, ind): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(|f| has(f, &y));
        let d = mul(&dep);
        let k = match d.kind() {
            Kind::Pow(b, e) if *b == y && !has(e, x) => e.clone(),
            _ => return None,
        };
        if n.is_some() {
            return None;
        }
        n = Some(k);
        q = mul(&ind);
    }
    let n = n?;
    // v = y^(1-n): v' + (1-n) p v = (1-n) q
    let one_n = sub(&Expr::one(), &n);
    let l = Linear { a: vec![mul2(&one_n, &p), Expr::one()], g: mul2(&one_n, &q) };
    let v = solve_linear(&l, 1, x)?;
    Some(vec![pow(&v, &div(&Expr::one(), &one_n))])
}

fn apply_ics(sol: &Expr, f: &Expr, x: &Expr, ics: &[(Expr, Expr)]) -> Option<Expr> {
    let mut eqs = Vec::new();
    for (cond, val) in ics {
        // f(x0) or Subs(Derivative(f(x), x, k), x, x0) given as f'(x0)
        let (k, x0) = match cond.kind() {
            Kind::Func(Func::Undefined(n), a) if &**n == "Subs" => {
                let Kind::Derivative(_, vars) = a[0].kind() else { return None };
                (vars[0].1 as usize, a[2].clone())
            }
            Kind::Func(Func::Undefined(_), a) => (0usize, a[0].clone()),
            _ => return None,
        };
        let mut d = sol.clone();
        for _ in 0..k {
            d = diff(&d, x);
        }
        eqs.push(sub(&subs1(&d, x, &x0), val));
    }
    let _ = f;
    let consts: Vec<Expr> = (1..=10).map(constant).filter(|c| has(sol, c)).collect();
    match solve_system(&eqs, &consts)? {
        SystemSolution::Map(m) => Some(expand(&subs(sol, &m))),
        _ => None,
    }
}

fn tidy(e: &Expr, x: &Expr) -> Expr {
    let _ = x;
    e.clone()
}

/// Check a candidate solution `y` of the ODE (residual is zero).
pub fn checkodesol(ode: &Expr, f: &Expr, y: &Expr) -> bool {
    let Kind::Func(_, args) = f.kind() else { return false };
    let x = &args[0];
    let mut order = 0;
    let e = to_symbols(&as_expr(ode), f, x, &mut order);
    let mut map = Vec::new();
    let mut d = y.clone();
    for k in 0..=order {
        map.push((ysym(k), d.clone()));
        d = diff(&d, x);
    }
    let r = expand(&subs(&e, &map));
    if r.is_zero() {
        return true;
    }
    let r = crate::simp::trigsimp(&r);
    if r.is_zero() {
        return true;
    }
    // numeric spot check with arbitrary constants
    let mut rs = r;
    for i in 1..=10 {
        rs = subs1(&rs, &constant(i), &Expr::rational(3 * i as i64 + 1, 7));
    }
    [0.31, 0.77, 1.3].iter().all(|v| {
        let val = subs1(&rs, x, &Expr::float(*v));
        crate::evalf::eval_c64(&val).is_some_and(|z| z.norm() < 1e-8)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn ds(s: &str) -> Vec<Expr> {
        let f = parse_str("f(x)").unwrap();
        let ode = parse_str(s).unwrap();
        let sols = dsolve(&ode, &f, &[]).unwrap_or_else(|| panic!("dsolve failed on {s}"));
        for sol in &sols {
            let Kind::Rel(_, _, y) = sol.kind() else { panic!() };
            assert!(checkodesol(&ode, &f, y), "{s}: {sol} does not check");
        }
        sols
    }

    #[test]
    fn odes() {
        assert_eq!(ds("Derivative(f(x), x) - f(x)")[0].to_string(), "Eq(f(x), C1*exp(x))");
        assert_eq!(ds("Derivative(f(x), x) + 2*f(x)")[0].to_string(), "Eq(f(x), C1*exp(-2*x))");
        assert_eq!(ds("Derivative(f(x), x) - x")[0].to_string(), "Eq(f(x), C1 + x**2/2)");
        for s in [
            "Derivative(f(x), x) - x*f(x)",
            "Derivative(f(x), x) + f(x) - x",
            "Derivative(f(x), (x, 2)) - f(x)",
            "Derivative(f(x), (x, 2)) + f(x)",
            "Derivative(f(x), (x, 2)) + 2*Derivative(f(x), x) + f(x)",
            "Derivative(f(x), (x, 2)) - 3*Derivative(f(x), x) + 2*f(x)",
            "Derivative(f(x), (x, 2)) + 4*f(x) - sin(x)",
            "Derivative(f(x), (x, 2)) - f(x) - exp(2*x)",
            "Derivative(f(x), x) - f(x)**2",
            "Derivative(f(x), (x, 3)) - f(x)",
            "Derivative(f(x), x) - f(x)/x",
            "x**2*Derivative(f(x), (x, 2)) + x*Derivative(f(x), x) - f(x)",
            "Derivative(f(x), (x, 2)) - 2*Derivative(f(x), x) + 5*f(x)",
            "Derivative(f(x), x) - sin(x)*f(x)",
            "Derivative(f(x), (x, 2)) - x",
            "Derivative(f(x), (x, 2)) - f(x) - exp(x)",
            "Derivative(f(x), x) + f(x) - f(x)**2",
        ] {
            ds(s);
        }
    }
}
