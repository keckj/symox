//! Method of undetermined coefficients for linear ODEs with constant
//! coefficients, as the reference (SymPy 1.14,
//! `sympy/solvers/ode/nonhomogeneous.py`) computes it:
//!
//! * `_undetermined_coefficients_match`: the trial set of a forcing term,
//!   each linearly independent derivative of each additive term, multiplied
//!   by the power of `x` that makes it independent of the homogeneous
//!   solutions (`is_homogeneous_solution` loop);
//! * `_solve_undetermined_coefficients`: the substitution of the trial
//!   function in the equation, the linear system of the coefficients of the
//!   independent terms (`separatevars` grouping) and the particular solution
//!   `trialfunc.subs(coeffvals)` added to the homogeneous solution.
//!
//! The native domain is the forcing terms that are sums of monomials
//! `q * x**k * exp(a*x) * T(b*x)` with `q`, `a`, `b` rational, `k` a
//! non-negative integer and `T` one of `sin`, `cos` (or absent), and the
//! equations `Σ c_n * Derivative(f(x), (x, n))` with rational `c_n`. On
//! that domain:
//!
//! * `powsimp(expr, combine='exp')` is the identity and `_test_term` is true;
//! * `_get_trial_set` of a monomial is the set
//!   `{x**j * exp(a*x) * T'(b*x) : j <= k, T' in {sin, cos}}` (only `T' = T`
//!   absent when `b = 0`). Its mutable default `exprs=set()` is only
//!   modified when an additive term expands to a sum, which is excluded;
//! * `is_homogeneous_solution(x**j * exp(a*x) * T(b*x))` holds exactly when
//!   `a + i*b` is a root of the characteristic polynomial of multiplicity
//!   larger than `j`, so the loop multiplies the set by `x**m`, `m` the
//!   multiplicity of `a + i*b`;
//! * the grouped coefficients form a linear system over QQ; when it has a
//!   unique solution, `solve` returns it and `trialfunc.subs(coeffvals)` is
//!   the canonical sum `Σ v_i * t_i`.
//!
//! Every function returns `None` outside that domain, when the system is
//! inconsistent or underdetermined (the reference raises or returns a
//! parametric solution) and when the homogeneous equation is zero; the
//! caller then runs the reference body unchanged.

use std::collections::BTreeMap;

use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::number::Number;
use crate::ops::{add, add2, mul, mul2, pow};

type Q = BigRational;

/// The trigonometric factor of a monomial.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Trig {
    None,
    Sin,
    Cos,
}

/// `x**k * exp(a*x) * T(b*x)`; `b` is zero exactly when `trig` is `None`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Mono {
    k: u32,
    a: Q,
    trig: Trig,
    b: Q,
}

impl Mono {
    fn one() -> Self {
        Mono { k: 0, a: Q::zero(), trig: Trig::None, b: Q::zero() }
    }
}

/// A linear combination of monomials with rational coefficients.
type Lin = BTreeMap<Mono, Q>;

fn accumulate(out: &mut Lin, m: Mono, c: Q) {
    if c.is_zero() {
        return;
    }
    let entry = out.entry(m).or_insert_with(Q::zero);
    *entry += c;
}

fn rational(e: &Expr) -> Option<Q> {
    match e.kind() {
        Kind::Num(n) if n.is_rational() => n.to_ratio(),
        _ => None,
    }
}

/// `c` of an argument `c*x` (`x` itself is `c = 1`).
fn linear_coefficient(arg: &Expr, x: &Expr) -> Option<Q> {
    if arg == x {
        return Some(Q::one());
    }
    match arg.kind() {
        Kind::Mul(c, factors) if c.is_rational() && factors.len() == 1 && &factors[0].0 == x && factors[0].1 == Expr::one() => c.to_ratio(),
        _ => None,
    }
}

fn factor_into(base: &Expr, exp: &Expr, x: &Expr, m: &mut Mono) -> Option<()> {
    if base == x {
        let k = exp.as_number().filter(|n| n.is_integer() && n.is_positive())?.to_i64()?;
        m.k = m.k.checked_add(u32::try_from(k).ok()?)?;
        return Some(());
    }
    if exp != &Expr::one() {
        return None;
    }
    match base.kind() {
        Kind::Pow(b, e) => factor_into(b, e, x, m),
        Kind::Func(Func::Exp, args) if args.len() == 1 && m.a.is_zero() => {
            let a = linear_coefficient(&args[0], x)?;
            if a.is_zero() {
                return None;
            }
            m.a = a;
            Some(())
        }
        Kind::Func(f @ (Func::Sin | Func::Cos), args) if args.len() == 1 && m.trig == Trig::None => {
            let b = linear_coefficient(&args[0], x)?;
            if !b.is_positive() {
                return None;
            }
            m.trig = if *f == Func::Sin { Trig::Sin } else { Trig::Cos };
            m.b = b;
            Some(())
        }
        _ => None,
    }
}

/// `(q, monomial)` of a term `q * x**k * exp(a*x) * T(b*x)`.
fn parse_monomial(term: &Expr, x: &Expr) -> Option<(Q, Mono)> {
    let mut m = Mono::one();
    let q = match term.kind() {
        Kind::Num(_) => return Some((rational(term)?, m)),
        Kind::Mul(c, factors) => {
            if !c.is_rational() {
                return None;
            }
            for (b, e) in factors {
                factor_into(b, e, x, &mut m)?;
            }
            c.to_ratio()?
        }
        _ => {
            factor_into(term, &Expr::one(), x, &mut m)?;
            Q::one()
        }
    };
    Some((q, m))
}

fn number(q: &Q) -> Expr {
    Expr::num(Number::from_ratio(q.clone()))
}

fn mono_expr(m: &Mono, x: &Expr) -> Expr {
    let mut factors = Vec::new();
    if m.k > 0 {
        factors.push(pow(x, &Expr::int(m.k as i64)));
    }
    if !m.a.is_zero() {
        factors.push(func(Func::Exp, vec![mul2(&number(&m.a), x)]));
    }
    match m.trig {
        Trig::None => {}
        Trig::Sin => factors.push(func(Func::Sin, vec![mul2(&number(&m.b), x)])),
        Trig::Cos => factors.push(func(Func::Cos, vec![mul2(&number(&m.b), x)])),
    }
    if factors.is_empty() {
        Expr::one()
    } else {
        mul(&factors)
    }
}

/// The canonical expression of `m`, checked to parse back to `m` with a unit
/// coefficient.
fn checked_mono_expr(m: &Mono, x: &Expr) -> Option<Expr> {
    let e = mono_expr(m, x);
    match parse_monomial(&e, x) {
        Some((q, back)) if q.is_one() && &back == m => Some(e),
        _ => None,
    }
}

/// `Add.make_args(e)` as monomials (the numeric part of an `Add` first).
fn additive_terms(e: &Expr, x: &Expr) -> Option<Vec<(Q, Mono)>> {
    let mut out = Vec::new();
    match e.kind() {
        Kind::Add(c, terms) => {
            if !c.is_zero() {
                out.push((c.to_ratio().filter(|_| c.is_rational())?, Mono::one()));
            }
            for (t, k) in terms {
                if !k.is_rational() {
                    return None;
                }
                let (q, m) = parse_monomial(t, x)?;
                out.push((q * k.to_ratio()?, m));
            }
        }
        _ => out.push(parse_monomial(e, x)?),
    }
    Some(out)
}

/// The order `n` of `Derivative(f(x), (x, n))` (`0` for `f(x)` itself).
fn derivative_order(e: &Expr, fx: &Expr, x: &Expr) -> Option<usize> {
    if e == fx {
        return Some(0);
    }
    match e.kind() {
        Kind::Derivative(inner, vars) if inner == fx && vars.len() == 1 && &vars[0].0 == x && vars[0].1 > 0 => Some(vars[0].1 as usize),
        _ => None,
    }
}

/// `c * D` with `c` rational and `D` an order of `f(x)`.
fn linear_term(t: &Expr, fx: &Expr, x: &Expr) -> Option<(Q, usize)> {
    if let Some(n) = derivative_order(t, fx, x) {
        return Some((Q::one(), n));
    }
    match t.kind() {
        Kind::Mul(c, factors) if c.is_rational() && factors.len() == 1 && factors[0].1 == Expr::one() => {
            let n = derivative_order(&factors[0].0, fx, x)?;
            Some((c.to_ratio()?, n))
        }
        _ => None,
    }
}

/// The operator `Σ p_n D^n` (by order) and the forcing monomials of `eq`.
fn parse_equation(eq: &Expr, fx: &Expr, x: &Expr) -> Option<(Vec<Q>, Lin)> {
    let mut ops: Vec<Q> = Vec::new();
    let mut forcing = Lin::new();
    let mut push = |t: &Expr, k: Q| -> Option<()> {
        if let Some((c, n)) = linear_term(t, fx, x) {
            if ops.len() <= n {
                ops.resize(n + 1, Q::zero());
            }
            ops[n] += c * k;
            return Some(());
        }
        let (q, m) = parse_monomial(t, x)?;
        accumulate(&mut forcing, m, q * k);
        Some(())
    };
    match eq.kind() {
        Kind::Add(c, terms) => {
            if !c.is_zero() {
                push(&Expr::num(c.clone()), Q::one())?;
            }
            for (t, k) in terms {
                if !k.is_rational() {
                    return None;
                }
                push(t, k.to_ratio()?)?;
            }
        }
        _ => push(eq, Q::one())?,
    }
    forcing.retain(|_, c| !c.is_zero());
    Some((ops, forcing))
}

fn derivative(p: &Lin) -> Lin {
    let mut out = Lin::new();
    for (m, c) in p {
        if m.k > 0 {
            let mut d = m.clone();
            d.k -= 1;
            accumulate(&mut out, d, c * Q::from_integer(m.k.into()));
        }
        if !m.a.is_zero() {
            accumulate(&mut out, m.clone(), c * &m.a);
        }
        match m.trig {
            Trig::None => {}
            Trig::Sin => {
                let mut d = m.clone();
                d.trig = Trig::Cos;
                accumulate(&mut out, d, c * &m.b);
            }
            Trig::Cos => {
                let mut d = m.clone();
                d.trig = Trig::Sin;
                accumulate(&mut out, d, -(c * &m.b));
            }
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

/// `L[m]` for `L = Σ ops[n] D^n`.
fn apply(ops: &[Q], m: &Mono) -> Lin {
    let mut out = Lin::new();
    let mut d: Lin = std::iter::once((m.clone(), Q::one())).collect();
    for (n, c) in ops.iter().enumerate() {
        if n > 0 {
            d = derivative(&d);
        }
        if c.is_zero() {
            continue;
        }
        for (mm, v) in &d {
            accumulate(&mut out, mm.clone(), c * v);
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

/// Gaussian rationals `re + i*im`.
#[derive(Clone)]
struct Gauss {
    re: Q,
    im: Q,
}

impl Gauss {
    fn mul(&self, o: &Gauss) -> Gauss {
        Gauss { re: &self.re * &o.re - &self.im * &o.im, im: &self.re * &o.im + &self.im * &o.re }
    }
}

/// The multiplicity of `a + i*b` as a root of `Σ ops[n] t**n`; `None` for
/// the zero polynomial.
fn multiplicity(ops: &[Q], a: &Q, b: &Q) -> Option<u32> {
    let mut p: Vec<Q> = ops.to_vec();
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    if p.is_empty() {
        return None;
    }
    let z = Gauss { re: a.clone(), im: b.clone() };
    let mut mult = 0u32;
    loop {
        if p.is_empty() {
            return None;
        }
        // Horner evaluation of p at z
        let mut acc = Gauss { re: Q::zero(), im: Q::zero() };
        for c in p.iter().rev() {
            acc = acc.mul(&z);
            acc.re += c;
        }
        if !(acc.re.is_zero() && acc.im.is_zero()) {
            return Some(mult);
        }
        mult += 1;
        // derivative of p
        p = p.iter().enumerate().skip(1).map(|(n, c)| c * Q::from_integer((n as i64).into())).collect();
    }
}

/// `_undetermined_coefficients_match(forcing, x, f(x), eq_homogeneous)['trialset']`
/// when `test` is true; `homogeneous` is `None` for `eq_homogeneous is S.Zero`
/// (no multiplication by powers of `x`). The elements are returned in no
/// particular order (the reference returns a set).
pub fn trial_set(forcing: &Expr, x: &Expr, fx: &Expr, homogeneous: Option<&Expr>) -> Option<Vec<Expr>> {
    if !x.is_symbol() {
        return None;
    }
    let ops = match homogeneous {
        Some(h) => {
            let (ops, rest) = parse_equation(h, fx, x)?;
            if !rest.is_empty() {
                return None;
            }
            Some(ops)
        }
        None => None,
    };
    let mut set: Vec<Mono> = Vec::new();
    for (_, m) in additive_terms(forcing, x)? {
        let shift = match &ops {
            Some(ops) => multiplicity(ops, &m.a, &m.b)?,
            None => 0,
        };
        let trigs: &[Trig] = match m.trig {
            Trig::None => &[Trig::None],
            _ => &[Trig::Sin, Trig::Cos],
        };
        for j in 0..=m.k {
            for t in trigs {
                let mono = Mono { k: j.checked_add(shift)?, a: m.a.clone(), trig: *t, b: m.b.clone() };
                if !set.contains(&mono) {
                    set.push(mono);
                }
            }
        }
    }
    set.iter().map(|m| checked_mono_expr(m, x)).collect()
}

/// Solves `A v = rhs` over QQ; `None` unless the solution exists and is
/// unique.
fn solve_unique(mut rows: Vec<Vec<Q>>, n: usize) -> Option<Vec<Q>> {
    // each row: n coefficients followed by the right-hand side
    let mut pivot_cols = Vec::new();
    let mut r = 0;
    for col in 0..n {
        let Some(p) = (r..rows.len()).find(|&i| !rows[i][col].is_zero()) else {
            return None;
        };
        rows.swap(r, p);
        let inv = Q::one() / &rows[r][col];
        for v in rows[r].iter_mut() {
            *v *= &inv;
        }
        let pivot = rows[r].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i != r && !row[col].is_zero() {
                let f = row[col].clone();
                for (v, pv) in row.iter_mut().zip(pivot.iter()) {
                    *v -= &f * pv;
                }
            }
        }
        pivot_cols.push(col);
        r += 1;
    }
    // remaining rows must be consistent (all zero, including rhs)
    if rows[r..].iter().any(|row| !row[n].is_zero()) {
        return None;
    }
    Some((0..n).map(|i| rows[i][n].clone()).collect())
}

/// `_solve_undetermined_coefficients(eq, f(x), order, match, trialset)`:
/// the right-hand side `gsol.rhs + psol` of the returned equation.
/// `gensols_len` is `len(match['list'])` and `gsol_rhs` is
/// `match['sol'].rhs`.
pub fn solve_undetermined(eq: &Expr, fx: &Expr, x: &Expr, order: usize, gensols_len: usize, gsol_rhs: &Expr, trialset: &[Expr]) -> Option<Expr> {
    if gensols_len != order || !x.is_symbol() || trialset.is_empty() {
        return None;
    }
    let (ops, forcing) = parse_equation(eq, fx, x)?;
    let mut monos = Vec::with_capacity(trialset.len());
    for t in trialset {
        let (q, m) = parse_monomial(t, x)?;
        if !q.is_one() || monos.contains(&m) || checked_mono_expr(&m, x).as_ref() != Some(t) {
            return None;
        }
        monos.push(m);
    }
    let images: Vec<Lin> = monos.iter().map(|m| apply(&ops, m)).collect();
    let mut keys: Vec<&Mono> = Vec::new();
    for k in images.iter().flat_map(|l| l.keys()).chain(forcing.keys()) {
        if !keys.contains(&k) {
            keys.push(k);
        }
    }
    let n = monos.len();
    let rows: Vec<Vec<Q>> = keys
        .iter()
        .map(|k| {
            let mut row: Vec<Q> = images.iter().map(|l| l.get(*k).cloned().unwrap_or_else(Q::zero)).collect();
            row.push(-forcing.get(*k).cloned().unwrap_or_else(Q::zero));
            row
        })
        .collect();
    let values = solve_unique(rows, n)?;
    let terms: Vec<Expr> = values.iter().zip(trialset).filter(|(v, _)| !v.is_zero()).map(|(v, t)| mul2(&number(v), t)).collect();
    let psol = add(&terms);
    Some(add2(gsol_rhs, &psol))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(s: &str) -> Expr {
        parse_str(s).unwrap()
    }

    #[test]
    fn monomials_round_trip() {
        let x = p("x");
        for s in ["x**2*exp(x)*sin(x)", "exp(-x/2)", "cos(3*x)", "x", "1", "x**3*cos(x/2)"] {
            let e = p(s);
            let (q, m) = parse_monomial(&e, &x).unwrap();
            assert!(q.is_one());
            assert_eq!(mono_expr(&m, &x), e, "{s}");
        }
        assert!(parse_monomial(&p("exp(x + 1)"), &x).is_none());
        assert!(parse_monomial(&p("sin(x)*cos(x)"), &x).is_none());
        assert!(parse_monomial(&p("log(x)"), &x).is_none());
        assert!(parse_monomial(&p("a*x"), &x).is_none());
    }

    #[test]
    fn multiplicities() {
        let q = |n: i64| Q::from_integer(n.into());
        // (t - 1)**2 (t**2 + 1) = t**4 - 2 t**3 + 2 t**2 - 2 t + 1
        let ops = vec![q(1), q(-2), q(2), q(-2), q(1)];
        assert_eq!(multiplicity(&ops, &q(1), &q(0)), Some(2));
        assert_eq!(multiplicity(&ops, &q(0), &q(1)), Some(1));
        assert_eq!(multiplicity(&ops, &q(2), &q(0)), Some(0));
        assert_eq!(multiplicity(&[q(0)], &q(0), &q(0)), None);
    }

    #[test]
    fn trial_sets_and_solutions() {
        let x = p("x");
        let fx = p("f(x)");
        // f'' + f = sin(x): resonance, trial set {x*sin(x), x*cos(x)}
        let h = p("Derivative(f(x), x, 2) + f(x)");
        let mut ts = trial_set(&p("sin(x)"), &x, &fx, Some(&h)).unwrap();
        ts.sort_by(crate::sortkey::cmp_default);
        let mut want = vec![p("x*sin(x)"), p("x*cos(x)")];
        want.sort_by(crate::sortkey::cmp_default);
        assert_eq!(ts, want);
        let eq = p("Derivative(f(x), x, 2) + f(x) - sin(x)");
        let rhs = solve_undetermined(&eq, &fx, &x, 2, 2, &p("C1*sin(x) + C2*cos(x)"), &ts).unwrap();
        assert_eq!(rhs, p("C1*sin(x) + C2*cos(x) - x*cos(x)/2"));
        // f'' - 3 f' + 2 f = x*exp(x)
        let h = p("Derivative(f(x), x, 2) - 3*Derivative(f(x), x) + 2*f(x)");
        let ts = trial_set(&p("x*exp(x)"), &x, &fx, Some(&h)).unwrap();
        assert_eq!(ts.len(), 2);
        let eq = p("Derivative(f(x), x, 2) - 3*Derivative(f(x), x) + 2*f(x) - x*exp(x)");
        let rhs = solve_undetermined(&eq, &fx, &x, 2, 2, &p("C1*exp(x) + C2*exp(2*x)"), &ts).unwrap();
        assert_eq!(rhs, p("C1*exp(x) + C2*exp(2*x) - x**2*exp(x)/2 - x*exp(x)"));
        // underdetermined trial set (no shift): declined
        let ts = trial_set(&p("exp(x)"), &x, &fx, None).unwrap();
        assert!(solve_undetermined(&eq, &fx, &x, 2, 2, &p("C1*exp(x) + C2*exp(2*x)"), &ts).is_none());
        // forcing outside the domain
        assert!(trial_set(&p("log(x)"), &x, &fx, None).is_none());
        assert!(trial_set(&p("x*(x + 1)"), &x, &fx, None).is_none());
    }
}
