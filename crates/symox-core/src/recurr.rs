//! Linear recurrences with polynomial coefficients, after SymPy's
//! `sympy/solvers/recurr.py`: `rsolve_poly` and `rsolve_ratio`
//! (Abramov's polynomial and rational solutions), `rsolve_hyper`
//! (Petkovsek's hypergeometric solutions) and the `rsolve` front end.
//!
//! Hypergeometric terms are recognised natively: a term is decomposed into
//! a rational function of `n`, exponentials `b**(c*n)` and gamma-like
//! factors with linear arguments, from which the ratio `f(n+1)/f(n)` and
//! the similarity test (`f/g` rational) are exact.

use num_bigint::BigInt;

use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::inspection::{classify, PolynomialKind};
use crate::matrix::Matrix;
use crate::number::Number;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sub};
use crate::special::SpecialFn;
use crate::subs::{subs, subs1, xreplace};
use crate::traverse::{free_symbols, has};

/// Errors raised by [`rsolve`], named after the Python exception classes.
#[derive(Debug, Clone)]
pub enum RsolveError {
    Value(String),
    Type(String),
}

/// `(solution, arbitrary constants)` as returned with `symbols=True`.
pub type Solution = (Expr, Vec<Expr>);

fn int(i: i64) -> Expr {
    Expr::int(i)
}

fn c_symbol(i: usize) -> Expr {
    Expr::symbol(&format!("C{i}"))
}

fn shifted(e: &Expr, n: &Expr, k: &Expr) -> Expr {
    subs1(e, n, &add2(n, k))
}

fn shift_by(e: &Expr, n: &Expr, k: i64) -> Expr {
    shifted(e, n, &int(k))
}

fn is_polynomial(e: &Expr, n: &Expr) -> bool {
    classify(e, std::slice::from_ref(n), PolynomialKind::Polynomial) == Some(true)
}

fn is_rational(e: &Expr, n: &Expr) -> bool {
    classify(e, std::slice::from_ref(n), PolynomialKind::Rational) == Some(true)
}

/// Dense ascending coefficients of a polynomial in `n`, without trailing
/// zeros (empty for the zero polynomial).
fn coefficients(e: &Expr, n: &Expr) -> Option<Vec<Expr>> {
    let mut c: Vec<Expr> = crate::integrate::poly_coeffs(&expand(e), n)?.iter().map(expand).collect();
    while c.last().is_some_and(|x| x.is_zero()) {
        c.pop();
    }
    Some(c)
}

fn degree(e: &Expr, n: &Expr) -> Option<i64> {
    coefficients(e, n).and_then(|c| if c.is_empty() { None } else { Some(c.len() as i64 - 1) })
}

fn binomial_int(n: i64, k: i64) -> Expr {
    if k < 0 || n < 0 || k > n {
        return Expr::zero();
    }
    Expr::num(Number::from_bigint(crate::special::integers::binomial(&BigInt::from(n), &BigInt::from(k))))
}

/// Keys of `roots(e, x)` in SymPy's order (empty when nothing is found).
fn root_keys(e: &Expr, x: &Expr) -> Vec<Expr> {
    if !has(e, x) {
        return vec![];
    }
    crate::poly::roots::roots(e, std::slice::from_ref(x)).map(|r| r.into_iter().map(|(z, _)| z).collect()).unwrap_or_default()
}

/// `roots(e, x, filter='Z', predicate=lambda r: r >= 0)`.
fn nonnegative_integer_roots(e: &Expr, x: &Expr) -> Vec<i64> {
    root_keys(e, x).iter().filter_map(|z| z.as_i64()).filter(|z| *z >= 0).collect()
}

fn quo(a: &Expr, b: &Expr, n: &Expr) -> Option<Expr> {
    crate::poly::api::quo_expr(a, b, std::slice::from_ref(n)).ok()
}

fn gcd(a: &Expr, b: &Expr, n: &Expr) -> Option<Expr> {
    crate::poly::api::gcd(a, b, std::slice::from_ref(n)).ok()
}

fn numer_denom(e: &Expr) -> (Expr, Expr) {
    crate::poly::api::as_numer_denom(e)
}

/// Outcome of SymPy's `solve(equations, *symbols)` on a linear system.
enum Linear {
    /// No nontrivial equation: SymPy returns `[]`.
    Trivial,
    Inconsistent,
    /// Pivot unknowns in terms of the free ones.
    Map(Vec<(Expr, Expr)>),
}

fn solve_linear(equations: &[Expr], unknowns: &[Expr]) -> Option<Linear> {
    let equations: Vec<Expr> = equations.iter().map(expand).filter(|e| !e.is_zero()).collect();
    if equations.is_empty() {
        return Some(Linear::Trivial);
    }
    if unknowns.is_empty() {
        return Some(Linear::Inconsistent);
    }
    let (a, b) = crate::linear::linear_eq_to_matrix(&equations, unknowns).ok()?;
    let Some(values) = crate::linear::linsolve(&a, &b, unknowns).ok()? else {
        return Some(Linear::Inconsistent);
    };
    let mut map = Vec::new();
    for (u, v) in unknowns.iter().zip(values) {
        if v != *u {
            // Domain arithmetic in SymPy returns cancelled values with
            // rationalised algebraic coefficients.
            let v = if free_symbols(&v).is_empty() { crate::simplify::simplify(&v) } else { crate::poly::api::cancel(&v).unwrap_or(v) };
            map.push((u.clone(), v));
        }
    }
    Some(Linear::Map(map))
}

// ----------------------------------------------------------------------
// hypergeometric terms
// ----------------------------------------------------------------------

/// `rational * prod(base**(slope*n)) * prod(gamma(arg)**multiplicity)`.
struct Term {
    rational: Expr,
    exponentials: Vec<(Expr, Expr)>,
    gammas: Vec<(Expr, i64)>,
}

/// `(slope, intercept)` of an expression of degree at most one in `n`.
fn linear(e: &Expr, n: &Expr) -> Option<(Expr, Expr)> {
    let c = coefficients(e, n)?;
    match c.len() {
        0 => Some((Expr::zero(), Expr::zero())),
        1 => Some((Expr::zero(), c[0].clone())),
        2 => Some((c[1].clone(), c[0].clone())),
        _ => None,
    }
}

/// Gamma arguments (with multiplicities) of a gamma-like function.
fn gamma_form(e: &Expr) -> Option<Vec<(Expr, i64)>> {
    let Kind::Func(f, a) = e.kind() else { return None };
    let one = Expr::one();
    Some(match f {
        Func::Gamma => vec![(a[0].clone(), 1)],
        Func::Special(SpecialFn::Factorial) => vec![(add2(&a[0], &one), 1)],
        Func::Special(SpecialFn::RisingFactorial) => vec![(add2(&a[0], &a[1]), 1), (a[0].clone(), -1)],
        Func::Special(SpecialFn::FallingFactorial) => vec![(add2(&a[0], &one), 1), (add(&[a[0].clone(), neg(&a[1]), one]), -1)],
        Func::Special(SpecialFn::Binomial) => vec![
            (add2(&a[0], &one), 1),
            (add2(&a[1], &one), -1),
            (add(&[a[0].clone(), neg(&a[1]), one]), -1),
        ],
        _ => return None,
    })
}

fn decompose_factor(term: &mut Term, b: &Expr, e: &Expr, n: &Expr) -> Option<()> {
    let (b, e) = match b.kind() {
        Kind::Func(Func::Exp, a) => (Expr::constant(Constant::E), mul2(&a[0], e)),
        _ => (b.clone(), e.clone()),
    };
    if !has(&b, n) && !has(&e, n) {
        term.rational = mul2(&term.rational, &pow(&b, &e));
        return Some(());
    }
    if !has(&b, n) {
        let (slope, intercept) = linear(&e, n)?;
        term.rational = mul2(&term.rational, &pow(&b, &intercept));
        term.exponentials.push((b, slope));
        return Some(());
    }
    let k = e.as_i64()?;
    if let Some(gammas) = gamma_form(&b) {
        for (a, m) in gammas {
            linear(&a, n)?;
            if has(&a, n) {
                term.gammas.push((a, m * k));
            } else {
                let g = crate::functions::gamma(&a);
                if g.is_constant(Constant::ComplexInfinity) || g.is_constant(Constant::NaN) || g.is_zero() {
                    return None;
                }
                term.rational = mul2(&term.rational, &pow(&g, &int(m * k)));
            }
        }
        return Some(());
    }
    if is_rational(&b, n) {
        term.rational = mul2(&term.rational, &pow(&b, &int(k)));
        return Some(());
    }
    None
}

fn decompose(t: &Expr, n: &Expr) -> Option<Term> {
    let mut term = Term { rational: Expr::one(), exponentials: vec![], gammas: vec![] };
    match t.kind() {
        Kind::Mul(c, factors) => {
            term.rational = Expr::num(c.clone());
            for (b, e) in factors {
                decompose_factor(&mut term, b, e, n)?;
            }
        }
        Kind::Pow(b, e) => decompose_factor(&mut term, b, e, n)?,
        _ => decompose_factor(&mut term, t, &Expr::one(), n)?,
    }
    Some(term)
}

/// `gamma(a + s)/gamma(a)` for an integer `s`.
fn gamma_shift(a: &Expr, s: i64) -> Expr {
    if s >= 0 {
        mul(&(0..s).map(|j| add2(a, &int(j))).collect::<Vec<_>>())
    } else {
        div(&Expr::one(), &mul(&(1..=-s).map(|j| sub(a, &int(j))).collect::<Vec<_>>()))
    }
}

fn term_ratio(t: &Term, n: &Expr) -> Option<Expr> {
    let mut factors = vec![div(&shift_by(&t.rational, n, 1), &t.rational)];
    for (b, s) in &t.exponentials {
        factors.push(pow(b, s));
    }
    for (a, m) in &t.gammas {
        let (slope, _) = linear(a, n)?;
        let s = slope.as_i64()?;
        factors.push(pow(&gamma_shift(a, s), &int(*m)));
    }
    let r = crate::poly::api::cancel(&mul(&factors)).ok()?;
    is_rational(&r, n).then_some(r)
}

/// `t1/t0` as a rational function of `n` when the terms are similar.
fn relative(t1: &Term, t0: &Term, n: &Expr) -> Option<Expr> {
    let mut groups: Vec<(Expr, Expr)> = Vec::new();
    for (b, s, sign) in t1.exponentials.iter().map(|(b, s)| (b, s, 1)).chain(t0.exponentials.iter().map(|(b, s)| (b, s, -1))) {
        let s = mul2(&int(sign), s);
        match groups.iter_mut().find(|g| g.0 == *b) {
            Some(g) => g.1 = expand(&add2(&g.1, &s)),
            None => groups.push((b.clone(), s)),
        }
    }
    let base = mul(&groups.iter().map(|(b, s)| pow(b, s)).collect::<Vec<_>>());
    if !expand(&base).is_one() {
        return None;
    }
    let mut classes: Vec<(Expr, Vec<(i64, i64)>)> = Vec::new();
    for (a, m) in t1.gammas.iter().map(|(a, m)| (a, *m)).chain(t0.gammas.iter().map(|(a, m)| (a, -*m))) {
        let mut placed = false;
        for (a0, members) in classes.iter_mut() {
            if let Some(offset) = expand(&sub(a, a0)).as_i64() {
                members.push((offset, m));
                placed = true;
                break;
            }
        }
        if !placed {
            classes.push((a.clone(), vec![(0, m)]));
        }
    }
    let mut factors = vec![div(&t1.rational, &t0.rational)];
    for (a0, members) in &classes {
        if members.iter().map(|(_, m)| m).sum::<i64>() != 0 {
            return None;
        }
        for (offset, m) in members {
            factors.push(pow(&gamma_shift(a0, *offset), &int(*m)));
        }
    }
    let r = mul(&factors);
    is_rational(&r, n).then_some(r)
}

/// `hypersimp(f, n)`: the rational function `f(n+1)/f(n)` in factored
/// form, or `None` when `f` is not recognised as hypergeometric.
pub fn hypersimp(f: &Expr, n: &Expr) -> Option<Expr> {
    let r = hyper_ratio(f, n)?;
    Some(crate::poly::api::factor(&r).unwrap_or(r))
}

fn hyper_ratio(f: &Expr, n: &Expr) -> Option<Expr> {
    if f.is_zero() {
        return None;
    }
    if let Some(t) = decompose(f, n) {
        return term_ratio(&t, n);
    }
    let ex = expand(f);
    let Kind::Add(..) = ex.kind() else { return None };
    let terms: Vec<Term> = ex.args().iter().map(|t| decompose(t, n)).collect::<Option<_>>()?;
    let rel = add(&terms.iter().map(|t| relative(t, &terms[0], n)).collect::<Option<Vec<_>>>()?);
    if rel.is_zero() {
        return None;
    }
    let ratio = mul(&[term_ratio(&terms[0], n)?, shift_by(&rel, n, 1), div(&Expr::one(), &rel)]);
    let r = crate::poly::api::cancel(&ratio).ok()?;
    is_rational(&r, n).then_some(r)
}

/// `f.is_hypergeometric(n)`.
pub fn is_hypergeometric(f: &Expr, n: &Expr) -> bool {
    hyper_ratio(f, n).is_some()
}

/// `hypersimilar(f, g, n)`: whether `f/g` is a rational function of `n`.
pub fn hypersimilar(f: &Expr, g: &Expr, n: &Expr) -> bool {
    match (decompose(f, n), decompose(g, n)) {
        (Some(a), Some(b)) => relative(&a, &b, n).is_some(),
        _ => false,
    }
}

// ----------------------------------------------------------------------
// Abramov: polynomial solutions
// ----------------------------------------------------------------------

/// `_delta(p, k)`: the `k`-th forward difference of `p` at `a`.
fn delta(p: &Expr, k: i64, a: i64, n: &Expr) -> Expr {
    let mut b = Number::int(1);
    let mut d = vec![subs1(p, n, &int(a + k))];
    for i in 1..=k {
        b = b.mul(&Number::rational_i64(i - k - 1, i).unwrap());
        d.push(mul2(&Expr::num(b.clone()), &subs1(p, n, &int(a + k - i))));
    }
    add(&d)
}

/// Rename the remaining constants `C` to the leading names of `all`.
fn renumber(result: Expr, remaining: Vec<Expr>, all: &[Expr]) -> Solution {
    if remaining.as_slice() == all {
        return (result, remaining);
    }
    let map = remaining.iter().cloned().zip(all.iter().cloned()).collect();
    let result = xreplace(&result, &map);
    (result, all[..remaining.len()].to_vec())
}

/// `rsolve_poly(coeffs, f, n, shift, symbols=True)`: polynomial solutions of
/// `sum(coeffs[i]*y(n+i)) = f` with constants numbered from `shift`.
pub fn rsolve_poly(coeffs: &[Expr], f: &Expr, n: &Expr, shift: usize) -> Option<Solution> {
    if !is_polynomial(f, n) || coeffs.is_empty() {
        return None;
    }
    let homogeneous = f.is_zero();
    let r = coeffs.len() - 1;
    let polys: Vec<Expr> = (0..=r)
        .map(|i| expand(&add(&(i..=r).map(|j| mul2(&coeffs[j], &binomial_int(j as i64, i as i64))).collect::<Vec<_>>())))
        .collect();
    let mut terms: Vec<(Expr, Option<i64>)> = Vec::with_capacity(r + 1);
    for p in &polys {
        let c = coefficients(p, n)?;
        terms.push(match c.last() {
            Some(lc) => (lc.clone(), Some(c.len() as i64 - 1)),
            None => (Expr::zero(), None),
        });
    }
    let mut d = terms[0].1;
    let mut b = terms[0].1;
    for (i, t) in terms.iter().enumerate().skip(1) {
        if t.1 > d {
            d = t.1;
        }
        let shifted = t.1.map(|e| e - i as i64);
        if shifted > b {
            b = shifted;
        }
    }
    let (d, b) = (d?, b?);
    let x = Expr::dummy("x", Default::default());
    let mut degree_poly = Vec::new();
    for (i, t) in terms.iter().enumerate() {
        if t.1.map(|e| e - i as i64) == Some(b) {
            let ff = mul(&(0..i as i64).map(|j| sub(&x, &int(j))).collect::<Vec<_>>());
            degree_poly.push(mul2(&t.0, &ff));
        }
    }
    let mut candidates: Vec<i64> = nonnegative_integer_roots(&expand(&add(&degree_poly)), &x).into_iter().max().into_iter().collect();
    if homogeneous {
        candidates.push(-b - 1);
    } else {
        candidates.push(degree(f, n).unwrap_or(0) - b);
        candidates.push(-b - 1);
    }
    let big_n = candidates.into_iter().max().unwrap();
    if big_n < 0 {
        return if homogeneous { Some((Expr::zero(), vec![])) } else { None };
    }
    if big_n as usize <= r {
        let all: Vec<Expr> = (0..=big_n as usize).map(|i| c_symbol(i + shift)).collect();
        let y = add(&all.iter().enumerate().map(|(i, c)| mul2(c, &pow(n, &int(i as i64)))).collect::<Vec<_>>());
        let e = add(&coeffs.iter().enumerate().map(|(i, c)| mul2(c, &shift_by(&y, n, i as i64))).collect::<Vec<_>>());
        // solve_undetermined_coeffs: every power of n must vanish
        let system = coefficients(&sub(&e, f), n)?;
        let present: Vec<Expr> = all.iter().filter(|c| system.iter().any(|s| has(s, c))).cloned().collect();
        let (result, remaining) = match solve_linear(&system, &present)? {
            Linear::Map(sol) => {
                let remaining = all.iter().filter(|c| !sol.iter().any(|(u, _)| u == *c)).cloned().collect();
                (subs(&y, &sol), remaining)
            }
            Linear::Trivial => (y, all.clone()),
            // SymPy keeps the bare ansatz here, which is no solution.
            Linear::Inconsistent => return None,
        };
        return Some(renumber(result, remaining, &all));
    }
    let big_n = big_n as usize;
    let a_ord = r;
    let u = (big_n as i64 + a_ord as i64 + b + 1) as usize;
    let a = nonnegative_integer_roots(&polys[r], n).into_iter().max().map_or(0, |m| m + 1);
    let mut alpha: Vec<Expr> = Vec::new();
    for i in -(a_ord as i64)..=d {
        let mut ii = vec![Expr::one(); d as usize + 1];
        for k in 1..=d as usize {
            ii[k] = div(&mul2(&ii[k - 1], &add2(&x, &int(i - k as i64 + 1))), &int(k as i64));
        }
        let mut s = Vec::new();
        for (j, p) in polys.iter().enumerate() {
            for (k, ik) in ii.iter().enumerate() {
                let bin = binomial_int(k as i64, i + j as i64);
                if bin.is_zero() {
                    continue;
                }
                s.push(mul(&[ik.clone(), bin, delta(p, k as i64, a, n)]));
            }
        }
        alpha.push(add(&s));
    }
    let alpha_at = |i: i64, at: i64| subs1(&alpha[(i + a_ord as i64) as usize], &x, &int(at));
    let mut v = vec![vec![Expr::zero(); a_ord]; u];
    for (i, row) in v.iter_mut().enumerate().take(a_ord.min(u)) {
        row[i] = Expr::one();
    }
    let mut g = vec![Expr::zero(); u];
    for i in a_ord..u {
        let mut acc = vec![Expr::zero(); a_ord];
        let mut gi = Expr::zero();
        for k in 1..=(a_ord as i64 + b) {
            if (i as i64) - k < 0 {
                break;
            }
            let bk = alpha_at(k - a_ord as i64, i as i64 - k);
            let row = (i as i64 - k) as usize;
            for (j, accj) in acc.iter_mut().enumerate() {
                *accj = add2(accj, &mul2(&bk, &v[row][j]));
            }
            if !homogeneous {
                gi = add2(&gi, &mul2(&bk, &g[row]));
            }
        }
        let denom = alpha_at(-(a_ord as i64), i as i64);
        for (j, accj) in acc.iter().enumerate() {
            v[i][j] = neg(&div(accj, &denom));
        }
        if !homogeneous {
            g[i] = div(&sub(&delta(f, i as i64 - a_ord as i64, a, n), &gi), &denom);
        }
    }
    let mut p = vec![Expr::one(); u];
    for i in 1..u {
        p[i] = expand(&div(&mul2(&p[i - 1], &sub(n, &int(a + i as i64 - 1))), &int(i as i64)));
    }
    let q: Vec<Expr> = (0..a_ord).map(|i| add(&(0..u).map(|k| expand(&mul2(&v[k][i], &p[k]))).collect::<Vec<_>>())).collect();
    let h = add(&(0..u).map(|k| expand(&mul2(&g[k], &p[k]))).collect::<Vec<_>>());
    let all: Vec<Expr> = (0..a_ord).map(|i| c_symbol(i + shift)).collect();
    let gfun = |i: i64| add(&all.iter().zip(&q).map(|(c, qq)| mul2(c, &delta(qq, i, a, n))).collect::<Vec<_>>());
    let equations: Vec<Expr> = (big_n + 1..u)
        .map(|i| if homogeneous { gfun(i as i64) } else { add2(&gfun(i as i64), &delta(&h, i as i64, a, n)) })
        .collect();
    let solutions = if equations.is_empty() {
        vec![]
    } else {
        match solve_linear(&equations, &all)? {
            Linear::Map(s) => s,
            Linear::Trivial | Linear::Inconsistent => {
                return if homogeneous { Some((Expr::zero(), vec![])) } else { None };
            }
        }
    };
    let mut result = vec![if homogeneous { Expr::zero() } else { h }];
    let mut remaining = Vec::new();
    for (c, qq) in all.iter().zip(&q) {
        let s = match solutions.iter().find(|(u, _)| u == c) {
            Some((_, value)) => mul2(value, qq),
            None => {
                remaining.push(c.clone());
                mul2(c, qq)
            }
        };
        result.push(expand(&s));
    }
    Some(renumber(add(&result), remaining, &all))
}

// ----------------------------------------------------------------------
// Abramov: rational solutions
// ----------------------------------------------------------------------

/// `rsolve_ratio(coeffs, f, n, symbols=True)`: rational solutions.
pub fn rsolve_ratio(coeffs: &[Expr], f: &Expr, n: &Expr) -> Option<Solution> {
    if !is_polynomial(f, n) || coeffs.is_empty() {
        return None;
    }
    let r = coeffs.len() - 1;
    let mut big_a = expand(&shift_by(&coeffs[r], n, -(r as i64)));
    let mut big_b = coeffs[0].clone();
    let h = Expr::dummy("h", Default::default());
    let mut res = crate::poly::api::resultant(&big_a, &shifted(&big_b, n, &h), &[n.clone(), h.clone()]).ok()?;
    if !is_polynomial(&res, &h) {
        let (p, q) = numer_denom(&res);
        res = quo(&p, &q, &h)?;
    }
    let Some(max_root) = nonnegative_integer_roots(&res, &h).into_iter().max() else {
        return rsolve_poly(coeffs, f, n, 0);
    };
    let mut c = Expr::one();
    for i in (0..=max_root).rev() {
        let d = gcd(&big_a, &shift_by(&big_b, n, i), n)?;
        big_a = quo(&big_a, &d, n)?;
        big_b = quo(&big_b, &shift_by(&d, n, -i), n)?;
        let factors: Vec<Expr> = (0..=i).map(|j| shift_by(&d, n, -j)).collect();
        c = mul2(&c, &mul(&factors));
    }
    let mut denoms: Vec<Expr> = (0..=r).map(|i| shift_by(&c, n, i as i64)).collect();
    let mut numers = vec![Expr::zero(); r + 1];
    for i in 0..=r {
        let g = gcd(&coeffs[i], &denoms[i], n)?;
        numers[i] = quo(&coeffs[i], &g, n)?;
        denoms[i] = quo(&denoms[i], &g, n)?;
    }
    for (i, numer) in numers.iter_mut().enumerate() {
        let others: Vec<Expr> = denoms.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, d)| d.clone()).collect();
        *numer = mul2(numer, &mul(&others));
    }
    let (result, syms) = rsolve_poly(&numers, &mul2(f, &mul(&denoms)), n, 0)?;
    Some((crate::simplify::simplify(&div(&result, &c)), syms))
}

// ----------------------------------------------------------------------
// Petkovsek: hypergeometric solutions
// ----------------------------------------------------------------------

fn integer_value(e: &Expr) -> Option<BigInt> {
    e.as_number().filter(|x| x.is_integer()).and_then(|x| x.to_bigint())
}

/// Casoratian `det([s_j(n + i)])`, using SymPy's explicit small formulas.
fn casoratian(seqs: &[Expr], n: &Expr) -> Expr {
    let k = seqs.len();
    let m = |i: usize, j: usize| shift_by(&seqs[j], n, i as i64);
    match k {
        1 => m(0, 0),
        2 => sub(&mul2(&m(0, 0), &m(1, 1)), &mul2(&m(0, 1), &m(1, 0))),
        3 => add(&[
            mul(&[m(0, 0), m(1, 1), m(2, 2)]),
            mul(&[m(0, 1), m(1, 2), m(2, 0)]),
            mul(&[m(0, 2), m(1, 0), m(2, 1)]),
            neg(&mul(&[m(0, 2), m(1, 1), m(2, 0)])),
            neg(&mul(&[m(0, 0), m(1, 2), m(2, 1)])),
            neg(&mul(&[m(0, 1), m(1, 0), m(2, 2)])),
        ]),
        _ => {
            let rows = (0..k).map(|i| (0..k).map(|j| m(i, j)).collect()).collect();
            Matrix::from_rows(rows).and_then(|a| a.det()).unwrap_or_else(|_| Expr::one())
        }
    }
}

fn has_factorial_like(e: &Expr) -> bool {
    if let Kind::Func(Func::Special(SpecialFn::Factorial | SpecialFn::FallingFactorial | SpecialFn::RisingFactorial), _) = e.kind() {
        return true;
    }
    e.args().iter().any(has_factorial_like)
}

/// Operation count in the spirit of SymPy's `count_ops`, used to choose
/// between equivalent forms like `simplify` does.
fn count_ops(e: &Expr) -> usize {
    match e.kind() {
        Kind::Add(..) => {
            let a = e.args();
            a.len() - 1 + a.iter().map(count_ops).sum::<usize>()
        }
        Kind::Mul(c, _) => {
            let (num, den) = numer_denom(e);
            if !den.is_one() {
                return 1 + count_ops(&num) + count_ops(&den);
            }
            let a = e.args();
            let neg = usize::from(c.is_negative());
            a.len() - 1 + neg + a.iter().map(count_ops).sum::<usize>()
        }
        Kind::Pow(b, x) => {
            if x.as_number().is_some_and(|x| x.is_negative()) {
                let (num, den) = numer_denom(e);
                return 1 + count_ops(&num) + count_ops(&den);
            }
            1 + count_ops(b) + count_ops(x)
        }
        Kind::Func(..) => 1 + e.args().iter().map(count_ops).sum::<usize>(),
        _ => 0,
    }
}

fn has_gamma_like(e: &Expr) -> bool {
    gamma_form(e).is_some() || e.args().iter().any(has_gamma_like)
}

/// Combine the gamma-like factors of a product whose arguments differ by
/// integers (`rf(2, n - 1)/factorial(n - 1)` is `n`), like `combsimp`.
fn combsimp_term(t: &Expr) -> Option<Expr> {
    let factors: Vec<(Expr, Expr)> = match t.kind() {
        Kind::Mul(c, f) => std::iter::once((Expr::num(c.clone()), Expr::one())).chain(f.iter().cloned()).collect(),
        Kind::Pow(b, e) => vec![(b.clone(), e.clone())],
        _ => vec![(t.clone(), Expr::one())],
    };
    let mut out = Vec::new();
    let mut gammas: Vec<(Expr, i64)> = Vec::new();
    for (b, e) in factors {
        match (gamma_form(&b), e.as_i64()) {
            (Some(form), Some(m)) => gammas.extend(form.into_iter().map(|(a, k)| (a, k * m))),
            _ => out.push(pow(&b, &e)),
        }
    }
    if gammas.is_empty() {
        return None;
    }
    let mut classes: Vec<(Expr, Vec<(i64, i64)>)> = Vec::new();
    for (a, m) in gammas {
        if free_symbols(&a).is_empty() {
            let g = crate::functions::gamma(&a);
            if g.is_constant(Constant::ComplexInfinity) || g.is_constant(Constant::NaN) || g.is_zero() {
                return None;
            }
            out.push(pow(&g, &int(m)));
            continue;
        }
        match classes.iter_mut().find_map(|(a0, members)| expand(&sub(&a, a0)).as_i64().map(|o| (members, o))) {
            Some((members, offset)) => members.push((offset, m)),
            None => classes.push((a, vec![(0, m)])),
        }
    }
    for (a0, members) in classes {
        let low = members.iter().map(|(o, _)| *o).min().unwrap();
        let base = expand(&add2(&a0, &int(low)));
        let total: i64 = members.iter().map(|(_, m)| m).sum();
        for (offset, m) in members {
            out.push(pow(&gamma_shift(&base, offset - low), &int(m)));
        }
        if total != 0 {
            let g = if crate::assumptions::is_integer(&base) == Some(true) {
                crate::special::factorial(&sub(&base, &Expr::one()))
            } else {
                crate::functions::gamma(&base)
            };
            out.push(pow(&g, &int(total)));
        }
    }
    Some(mul(&out))
}

fn combsimp(e: &Expr) -> Expr {
    if !has_gamma_like(e) {
        return e.clone();
    }
    let terms = if let Kind::Add(..) = e.kind() { e.args() } else { vec![e.clone()] };
    let candidate = add(&terms.iter().map(|t| combsimp_term(t).unwrap_or_else(|| t.clone())).collect::<Vec<_>>());
    if count_ops(&candidate) as f64 > 1.7 * count_ops(e) as f64 {
        e.clone()
    } else {
        candidate
    }
}

/// The sign normalisation `simplify` applies to a kernel term: integer
/// offsets of a `(-1)**(e + c)` factor move into the coefficient, and a bare
/// `-1` coefficient is absorbed back, so `(-1)**(n - 1)*n` is
/// `(-1)**(n + 1)*n` while `(-1)**(n - 3)*n/6` is `-(-1)**n*n/6`.
fn sign_merge(t: &Expr) -> Expr {
    let (mut coeff, factors) = match t.kind() {
        Kind::Mul(c, f) => (Expr::num(c.clone()), f.clone()),
        Kind::Pow(b, e) => (Expr::one(), vec![(b.clone(), e.clone())]),
        _ => return t.clone(),
    };
    let mut out = Vec::new();
    let mut sign_exp = None;
    for (b, e) in &factors {
        if b.as_i64() == Some(-1) && sign_exp.is_none() {
            let (offset, rest) = match e.kind() {
                Kind::Add(..) => {
                    let (nums, others): (Vec<Expr>, Vec<Expr>) = e.args().into_iter().partition(|a| a.as_number().is_some());
                    (add(&nums), add(&others))
                }
                _ => (Expr::zero(), e.clone()),
            };
            if let Some(o) = offset.as_i64() {
                if o % 2 != 0 {
                    coeff = crate::ops::neg(&coeff);
                }
                sign_exp = Some(rest);
                continue;
            }
        }
        out.push(pow(b, e));
    }
    let Some(mut x) = sign_exp else {
        return t.clone();
    };
    if coeff.as_i64() == Some(-1) {
        x = add2(&x, &Expr::one());
        coeff = Expr::one();
    }
    out.push(coeff);
    out.push(pow(&Expr::int(-1), &x));
    mul(&out)
}

/// `rsolve_hyper(coeffs, f, n, symbols=True)`: hypergeometric solutions.
pub fn rsolve_hyper(coeffs: &[Expr], f: &Expr, n: &Expr) -> Option<Solution> {
    if coeffs.is_empty() {
        return None;
    }
    let r = coeffs.len() - 1;
    let mut result = Expr::zero();
    if !f.is_zero() {
        let mut inhomogeneous: Vec<Expr> = if let Kind::Add(..) = f.kind() {
            let mut similar: Vec<(Expr, Vec<Expr>)> = Vec::new();
            let ex = expand(f);
            let parts = if let Kind::Add(..) = ex.kind() { ex.args() } else { vec![ex] };
            for g in parts {
                if !is_hypergeometric(&g, n) {
                    return None;
                }
                match similar.iter_mut().find(|(h, _)| hypersimilar(&g, h, n)) {
                    Some((_, acc)) => acc.push(g),
                    None => similar.push((g, vec![])),
                }
            }
            similar.into_iter().map(|(g, acc)| add2(&g, &add(&acc))).collect()
        } else if is_hypergeometric(f, n) {
            vec![f.clone()]
        } else {
            return None;
        };
        for i in 0..inhomogeneous.len() {
            let g = inhomogeneous[i].clone();
            let mut coeff = Expr::one();
            let mut polys = coeffs.to_vec();
            let mut denoms = vec![Expr::one(); r + 1];
            let s = hypersimp(&g, n)?;
            for j in 1..=r {
                coeff = mul2(&coeff, &shift_by(&s, n, j as i64 - 1));
                let (p, q) = numer_denom(&coeff);
                polys[j] = mul2(&polys[j], &p);
                denoms[j] = q;
            }
            for (j, poly) in polys.iter_mut().enumerate() {
                let others: Vec<Expr> = denoms.iter().enumerate().filter(|(i, _)| *i != j).map(|(_, d)| d.clone()).collect();
                *poly = mul2(poly, &mul(&others));
            }
            let rhs = mul(&denoms);
            let particular = match rsolve_ratio(&polys, &rhs, n) {
                Some((value, syms)) => {
                    let zero: Vec<(Expr, Expr)> = syms.iter().map(|s| (s.clone(), Expr::zero())).collect();
                    Some(if zero.is_empty() { value } else { subs(&value, &zero) })
                }
                None => rsolve_poly(&polys, &rhs, n, 0).map(|(value, _)| value),
            };
            match particular {
                Some(value) if !value.is_zero() => inhomogeneous[i] = mul2(&inhomogeneous[i], &value),
                _ => return None,
            }
            result = crate::simplify::simplify(&combsimp(&add(&inhomogeneous)));
        }
    }
    let z = Expr::dummy("Z", Default::default());
    let p = coeffs[0].clone();
    let q = shift_by(&coeffs[r], n, 1 - r as i64);
    let p_factors = root_keys(&p, n);
    let q_factors = root_keys(&q, n);
    let mut factors: Vec<(Expr, Expr)> = p_factors.iter().map(|p| (sub(n, p), Expr::one())).collect();
    factors.push((Expr::one(), Expr::one()));
    for p in &p_factors {
        for q in &q_factors {
            if let (Some(pi), Some(qi)) = (integer_value(p), integer_value(q)) {
                if pi <= qi {
                    continue;
                }
            }
            factors.push((sub(n, p), sub(n, q)));
        }
    }
    factors.extend(q_factors.iter().map(|q| (Expr::one(), sub(n, q))));
    let mut kernel: Vec<Expr> = Vec::new();
    for (big_a, big_b) in &factors {
        let d = mul2(big_a, &shift_by(big_b, n, r as i64 - 1));
        let mut polys = Vec::with_capacity(r + 1);
        let mut degrees = Vec::new();
        for (i, c) in coeffs.iter().enumerate() {
            let a = mul(&(0..i).map(|j| shift_by(big_a, n, j as i64)).collect::<Vec<_>>());
            let b = mul(&(i..r).map(|j| shift_by(big_b, n, j as i64)).collect::<Vec<_>>());
            let poly = quo(&mul(&[c.clone(), a, b]), &d, n)?;
            let cs = coefficients(&poly, n)?;
            if !cs.is_empty() {
                degrees.push(cs.len() - 1);
            }
            polys.push(cs);
        }
        let Some(&dd) = degrees.iter().max() else { return None };
        let mut zpoly = Vec::new();
        for (i, cs) in polys.iter().enumerate() {
            if let Some(c) = cs.get(dd).filter(|c| !c.is_zero()) {
                zpoly.push(mul2(c, &pow(&z, &int(i as i64))));
            }
        }
        let poly_exprs: Vec<Expr> = polys.iter().map(|cs| add(&cs.iter().enumerate().map(|(k, c)| mul2(c, &pow(n, &int(k as i64)))).collect::<Vec<_>>())).collect();
        for zr in root_keys(&add(&zpoly), &z) {
            if zr.is_zero() {
                continue;
            }
            let recurr: Vec<Expr> = poly_exprs.iter().enumerate().map(|(i, p)| mul2(p, &pow(&zr, &int(i as i64)))).collect();
            let sol: Vec<Expr> = if dd == 0 && !add(&recurr.iter().enumerate().skip(1).map(|(j, c)| mul2(c, &int(j as i64))).collect::<Vec<_>>()).is_zero() {
                vec![c_symbol(0)]
            } else {
                let (value, syms) = rsolve_poly(&recurr, &Expr::zero(), n, 0)?;
                let collected = crate::poly::api::collect(&value, &syms);
                syms.iter().map(|s| coefficient(&collected, s)).collect()
            };
            for c in sol {
                let ratio = crate::simplify::simplify(&div(&mul(&[zr.clone(), big_a.clone(), shift_by(&c, n, 1)]), &mul2(big_b, &c)));
                let mut n0 = Expr::zero();
                for root in root_keys(&numer_denom(&ratio).1, n) {
                    if has(&root, &Expr::constant(Constant::I)) {
                        return None;
                    }
                    let candidate = add2(&root, &Expr::one());
                    if let (Some(x), Some(y)) = (n0.as_number(), candidate.as_number()) {
                        if x.cmp_value(&y) == Some(std::cmp::Ordering::Less) {
                            n0 = candidate;
                        }
                    }
                }
                let mut k = crate::products::product(&ratio, vec![vec![n.clone(), n0, sub(n, &Expr::one())]]).ok()?;
                if has_factorial_like(&k) {
                    k = sign_merge(&combsimp(&k));
                }
                let mut seqs = kernel.clone();
                seqs.push(k.clone());
                if !casoratian(&seqs, n).is_zero() {
                    kernel.push(k);
                }
            }
        }
    }
    kernel.sort_by(crate::sortkey::cmp_default);
    let symbols: Vec<Expr> = (0..kernel.len()).map(c_symbol).collect();
    let mut terms = vec![result];
    for (c, k) in symbols.iter().zip(&kernel) {
        terms.push(mul2(c, k));
    }
    Some((add(&terms), symbols))
}

/// `e.coeff(s)` for an expression linear in `s` after `collect`.
fn coefficient(e: &Expr, s: &Expr) -> Expr {
    let terms = if let Kind::Add(..) = e.kind() { e.args() } else { vec![e.clone()] };
    let mut out = Vec::new();
    for t in terms {
        if t == *s {
            out.push(Expr::one());
            continue;
        }
        let Kind::Mul(..) = t.kind() else { continue };
        let factors = t.args();
        if let Some(i) = factors.iter().position(|f| f == s) {
            let rest: Vec<Expr> = factors.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, f)| f.clone()).collect();
            let rest = mul(&rest);
            if !has(&rest, s) {
                out.push(rest);
            }
        }
    }
    add(&out)
}

// ----------------------------------------------------------------------
// rsolve
// ----------------------------------------------------------------------

/// An initial condition key: an integer index or a term `y(i)`.
fn init_index(k: &Expr, y: &Expr) -> Result<i64, RsolveError> {
    if let Some(i) = k.as_number().filter(|x| x.is_integer()).and_then(|_| k.as_i64()) {
        return Ok(i);
    }
    if let (Kind::Func(f, a), Kind::Func(g, _)) = (k.kind(), y.kind()) {
        if f == g && a.len() == 1 {
            if let Some(i) = a[0].as_i64() {
                return Ok(i);
            }
        }
    }
    Err(RsolveError::Value(format!("Integer or term expected, got '{k}'")))
}

/// `rsolve(f, y(n), init)`: solve the linear recurrence `f = 0` (or an
/// equation) for `y(n)`. `init` holds `(index or y(i), value)` pairs.
/// `Ok(None)` mirrors SymPy returning `None`.
pub fn rsolve(f: &Expr, y: &Expr, init: Option<&[(Expr, Expr)]>) -> Result<Option<Expr>, RsolveError> {
    let f = match f.kind() {
        Kind::Rel(RelOp::Eq, lhs, rhs) => sub(lhs, rhs),
        _ => f.clone(),
    };
    let Kind::Func(head, yargs) = y.kind() else {
        return Err(RsolveError::Value(format!("'{y}' is not a function application")));
    };
    let n = yargs.first().cloned().ok_or_else(|| RsolveError::Value(format!("'{y}' is not a function application")))?;
    let fname = head.name().to_string();
    let is_y = |e: &Expr| matches!(e.kind(), Kind::Func(g, _) if g == head);
    let ex = expand(&f);
    let terms = if let Kind::Add(..) = ex.kind() { ex.args() } else { vec![ex] };
    let mut h_part: Vec<(i64, Vec<Expr>)> = Vec::new();
    let mut i_part = Vec::new();
    for g in terms {
        if let Some(x) = g.as_number().filter(|x| !x.is_rational()) {
            // SymPy's `Float.as_coeff_mul(y)` reports the float as a factor
            return Err(RsolveError::Value(format!("'{fname}({n} + k)' expected, got '{}'", Expr::num(x.abs()))));
        }
        let factors = if let Kind::Mul(..) = g.kind() { g.args() } else { vec![g.clone()] };
        let (dep, coeff): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(|x| has_head(x, head));
        let coeff = mul(&coeff);
        if dep.is_empty() {
            i_part.push(coeff);
            continue;
        }
        let expected = |h: &Expr| RsolveError::Value(format!("'{fname}({n} + k)' expected, got '{h}'"));
        if dep.len() > 1 {
            return Err(expected(&mul(&dep)));
        }
        let h = &dep[0];
        if !is_y(h) || h.args().len() != 1 {
            return Err(expected(h));
        }
        let k = expand(&sub(&h.args()[0], &n));
        if has(&k, &n) {
            return Err(expected(h));
        }
        let Some(k) = k.as_i64() else {
            if k.as_number().is_some() {
                return Err(expected(h));
            }
            return Err(RsolveError::Type("Cannot convert symbols to int".into()));
        };
        match h_part.iter_mut().find(|(i, _)| *i == k) {
            Some((_, v)) => v.push(coeff),
            None => h_part.push((k, vec![coeff])),
        }
    }
    let mut h_part: Vec<(i64, Expr)> = h_part.into_iter().map(|(k, v)| (k, crate::simplify::simplify(&add(&v)))).collect();
    let mut i_part = add(&i_part);
    if !i_part.is_zero() && !is_hypergeometric(&i_part, &n) {
        let all = matches!(i_part.kind(), Kind::Add(..)) && {
            let ex = expand(&i_part);
            let parts = if let Kind::Add(..) = ex.kind() { ex.args() } else { vec![ex] };
            parts.iter().all(|t| is_hypergeometric(t, &n))
        };
        if !all {
            return Err(RsolveError::Value(format!("The independent term should be a sum of hypergeometric functions, got '{i_part}'")));
        }
    }
    let ns = std::slice::from_ref(&n);
    let mut common = Expr::one();
    for (_, coeff) in &h_part {
        if !is_rational(coeff, &n) {
            return Err(RsolveError::Value(format!("Polynomial or rational function expected, got '{coeff}'")));
        }
        if !is_polynomial(coeff, &n) {
            common = crate::poly::api::lcm(&common, &numer_denom(coeff).1, ns).map_err(|e| RsolveError::Value(format!("{e:?}")))?;
        }
    }
    let (i_numer, i_denom) = numer_denom(&i_part);
    if is_polynomial(&i_denom, &n) {
        common = crate::poly::api::lcm(&common, &i_denom, ns).map_err(|e| RsolveError::Value(format!("{e:?}")))?;
    }
    if !common.is_one() {
        for (_, coeff) in h_part.iter_mut() {
            let (numer, denom) = numer_denom(coeff);
            *coeff = mul2(&numer, &quo(&common, &denom, &n).ok_or_else(|| RsolveError::Value("division failed".into()))?);
        }
        i_part = mul2(&i_numer, &quo(&common, &i_denom, &n).ok_or_else(|| RsolveError::Value("division failed".into()))?);
    }
    let k_min = h_part.iter().map(|(k, _)| *k).min().ok_or_else(|| RsolveError::Value("min() arg is an empty sequence".into()))?;
    if k_min < 0 {
        let k = -k_min;
        i_part = expand(&shift_by(&i_part, &n, k));
        for (i, coeff) in h_part.iter_mut() {
            *i += k;
            *coeff = expand(&shift_by(coeff, &n, k));
        }
    }
    let k_max = h_part.iter().map(|(k, _)| *k).max().unwrap();
    let coeffs: Vec<Expr> = (0..=k_max).map(|i| h_part.iter().find(|(k, _)| *k == i).map_or_else(Expr::zero, |(_, c)| c.clone())).collect();
    let Some((mut solution, symbols)) = rsolve_hyper(&coeffs, &neg(&i_part), &n) else {
        return Ok(None);
    };
    let init = init.filter(|i| !i.is_empty());
    if let (false, Some(init)) = (symbols.is_empty(), init) {
        let mut equations = Vec::new();
        for (k, v) in init {
            let i = init_index(k, y)?;
            let mut eq = sub(&subs1(&solution, &n, &int(i)), v);
            if has(&eq, &Expr::constant(Constant::NaN)) {
                let l = crate::limits::limit(&solution, &n, &int(i), crate::limits::Dir::Plus).map_err(|e| RsolveError::Value(format!("{e:?}")))?;
                eq = sub(&l, v);
            }
            equations.push(eq);
        }
        match solve_linear(&equations, &symbols) {
            Some(Linear::Map(sol)) if !sol.is_empty() => solution = subs(&solution, &sol),
            _ => return Ok(None),
        }
    }
    Ok(Some(solution))
}

/// `difference_delta(expr, n, step)`: the forward difference
/// `expr(n + step) - expr(n)`; `n` defaults to the only free symbol.
pub fn difference_delta(e: &Expr, n: Option<&Expr>, step: &Expr) -> Result<Expr, String> {
    let n = match n {
        Some(n) => n.clone(),
        None => {
            let free = free_symbols(e);
            match free.len() {
                0 => return Ok(Expr::zero()),
                1 => free.into_iter().next().unwrap(),
                _ => {
                    return Err(format!(
                        "Since there is more than one variable in the expression, a variable must be supplied to take the difference of {e}"
                    ))
                }
            }
        }
    };
    if !free_symbols(step).is_empty() || crate::assumptions::is_finite(step) == Some(false) {
        return Err("Step should be a finite number.".into());
    }
    Ok(delta_of(e, &n, step))
}

/// SymPy's `_eval_difference_delta` rules (sums term by term, the product
/// rule on the first factor, `b**e(n)` and the new terms of a `Sum`), with
/// `expr(n + step) - expr(n)` when no rule applies or a rule gives zero.
fn delta_of(e: &Expr, n: &Expr, step: &Expr) -> Expr {
    let shifted = |x: &Expr| crate::subs::subs1(x, n, &add2(n, step));
    let rule = match e.kind() {
        Kind::Add(..) => Some(add(&e.args().iter().map(|a| delta_of(a, n, step)).collect::<Vec<_>>())),
        Kind::Mul(..) => {
            let args = e.args();
            let rest = mul(&args[1..]);
            Some(add2(&mul2(&shifted(&args[0]), &delta_of(&rest, n, step)), &mul2(&delta_of(&args[0], n, step), &rest)))
        }
        Kind::Pow(b, x) if has(x, n) && !has(b, n) => Some(mul2(&sub(&pow(b, &sub(&shifted(x), x)), &Expr::one()), e)),
        Kind::Func(Func::Sum, args) => {
            let mut limits = crate::binders::limits(args);
            let last = limits.pop().filter(|l| l.len() == 3);
            last.and_then(|l| {
                let body = if limits.is_empty() { args[0].clone() } else { crate::binders::held(Func::Sum, &args[0], limits).ok()? };
                let tail = vec![vec![l[0].clone(), add2(&l[2], &Expr::one()), shifted(&l[2])]];
                Some(crate::functions::doit(&crate::binders::held(Func::Sum, &body, tail).ok()?, true))
            })
        }
        _ => None,
    };
    match rule {
        Some(r) if !r.is_zero() => r,
        _ => sub(&shifted(e), e),
    }
}

fn has_head(e: &Expr, head: &Func) -> bool {
    if let Kind::Func(g, _) = e.kind() {
        if g == head {
            return true;
        }
    }
    e.args().iter().any(|a| has_head(a, head))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(s: &str) -> Expr {
        parse_str(s).unwrap()
    }

    #[test]
    fn hypergeometric_ratios() {
        let n = Expr::symbol("n");
        assert_eq!(hypersimp(&p("2**n*n"), &n).unwrap().to_string(), "2*(n + 1)/n");
        assert_eq!(hypersimp(&p("factorial(n)"), &n).unwrap(), p("n + 1"));
        assert!(hypersimp(&p("sin(n)"), &n).is_none());
        assert!(hypersimp(&p("2**(n**2)"), &n).is_none());
        assert!(hypersimilar(&p("2**n*n"), &p("2**n"), &n));
        assert!(!hypersimilar(&p("2**n"), &p("3**n"), &n));
        assert!(hypersimilar(&p("4**n"), &p("2**(2*n)"), &n));
    }

    #[test]
    fn abramov_and_petkovsek() {
        let n = Expr::symbol("n");
        let (r, c) = rsolve_poly(&[p("-1"), p("1")], &n, &n, 0).unwrap();
        assert_eq!(r, p("C0 + n**2/2 - n/2"));
        assert_eq!(c, vec![p("C0")]);
        let (r, _) = rsolve_hyper(&[p("-2"), p("1")], &Expr::zero(), &n).unwrap();
        assert_eq!(r, p("2**n*C0"));
        assert!(rsolve_poly(&[p("n"), p("-n")], &Expr::one(), &n, 0).is_none());
    }
}
