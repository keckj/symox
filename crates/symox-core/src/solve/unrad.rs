//! SymPy 1.14's `unrad` (`sympy/solvers/solvers.py`): remove the radicals
//! of an equation with symbolic powers, returning `(eq, cov)` or None.
//!
//! The port follows the reference step by step: the preconditioning
//! (`powdenest(factor_terms(eq, radical=True, clear=True))`, numerator,
//! `_mexpand`), the generators of `eq.as_poly()` that `_take` selects,
//! the grouping of the terms by their radicals, the `ok` cases (one group,
//! one common base with the change of variable `covsym**lcm = base`, two
//! groups with or without a change of variable, three groups with two cube
//! roots, the power-of-two cases), the repetition guard and `_canonical`.
//! The NotImplementedError cases are the reference ones.
//!
//! The `_vsolve(..., check=False, simplify=False)` calls and a custom
//! `_take` go through an [`UnradHost`].

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::subs::{subs1, xreplace};
use crate::traverse::free_symbols;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use rustc_hash::FxHashMap;

/// Why `unrad` stopped: the reference NotImplementedError, a polynomial
/// error, an index or attribute error of the reference code, or a host
/// error (from the solve callback).
pub enum UFail<E> {
    NotImplemented(String),
    Poly(String),
    Index(String),
    Attribute(String),
    Host(E),
}

type R<T, E> = Result<T, UFail<E>>;

fn ni<T, E>(msg: &str) -> R<T, E> {
    Err(UFail::NotImplemented(msg.to_string()))
}

/// The calls back into the public API.
pub trait UnradHost {
    type Error;
    /// `_vsolve(eq, x, check=False, simplify=False)`. A NotImplementedError
    /// of the solver is returned as `UFail::NotImplemented`.
    fn vsolve(&mut self, eq: &Expr, x: &Expr) -> R<Vec<Expr>, Self::Error>;
    /// A custom `_take` passed by the caller; None for the reference one.
    fn take(&mut self, e: &Expr) -> Option<R<bool, Self::Error>>;
}

/// The result: the radical-free equation and the change of variable
/// (empty, or `[p, relation]`).
pub type Unrad = (Expr, Vec<Expr>);

// ---------------------------------------------------------------- helpers

fn is_sym_pow(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Pow(b, _) if !b.is_constant(Constant::E))
}

/// The reference `as_base_exp`.
fn as_base_exp(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Pow(b, x) => {
            if let Some(n) = b.as_number() {
                if let Some((p, q)) = n.numer_denom() {
                    if p.is_one() && !q.is_one() {
                        return (Expr::num(Number::from_bigint(q)), neg(x));
                    }
                }
            }
            (b.clone(), x.clone())
        }
        Kind::Func(crate::expr::Func::Exp, a) if a.len() == 1 => (Expr::e(), a[0].clone()),
        _ => (e.clone(), Expr::one()),
    }
}

/// `e.as_coeff_Mul(rational)`.
fn as_coeff_mul(e: &Expr, rational: bool) -> (Number, Expr) {
    let ok = |n: &Number| !rational || n.is_rational();
    match e.kind() {
        Kind::Num(n) if ok(n) => (n.clone(), Expr::one()),
        Kind::Mul(c, _) if !c.is_one() && ok(c) => (c.clone(), div(e, &Expr::num(c.clone()))),
        _ => (Number::ONE, e.clone()),
    }
}

fn make_args_add(e: &Expr) -> Vec<Expr> {
    if e.is_add() {
        crate::traverse::ordered_args(e)
    } else {
        vec![e.clone()]
    }
}

fn make_args_mul(e: &Expr) -> Vec<Expr> {
    if e.is_mul() {
        crate::traverse::ordered_args(e)
    } else {
        vec![e.clone()]
    }
}

fn args_of(e: &Expr) -> Vec<Expr> {
    crate::traverse::ordered_args(e)
}

fn intersects(e: &Expr, syms: &[Expr]) -> bool {
    free_symbols(e).iter().any(|s| syms.contains(s))
}

fn has_any(e: &Expr, pats: &[Expr]) -> bool {
    pats.iter().any(|p| crate::traverse::has(e, p))
}

/// Structural `has_free` of a generator.
fn has_free(e: &Expr, gens: &[Expr]) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |n| {
        if !found && gens.contains(n) {
            found = true;
        }
    });
    found
}

/// `_Q(pow)`: the denominator of the leading Rational of the exponent.
fn q_of(p: &Expr) -> BigInt {
    let (_, x) = as_base_exp(p);
    let (c, _) = as_coeff_mul(&x, false);
    match c.numer_denom() {
        Some((_, q)) if c.is_rational() => q,
        _ => BigInt::one(),
    }
}

fn mexpand(e: &Expr) -> Expr {
    crate::solve::tsolve::mexpand(e, true)
}

fn numer(e: &Expr) -> Expr {
    crate::simptools::as_numer_denom(e).0
}

/// SymPy's `is_number`.
fn is_number(e: &Expr) -> bool {
    crate::symbolic_number::is_symbolic_number(e)
}

/// `decompose_power(factor)` with the negative-exponent rule of the
/// polynomial constructor.
fn decompose_power(f: &Expr) -> (Expr, BigInt) {
    let (b, x) = as_base_exp(f);
    let (base, e) = if let Some(n) = x.as_number() {
        match n.numer_denom() {
            Some((p, q)) => {
                if q.is_one() {
                    (b, p)
                } else {
                    (pow(&b, &Expr::num(Number::from_ratio(num_rational::BigRational::new(BigInt::one(), q)))), p)
                }
            }
            None => (f.clone(), BigInt::one()),
        }
    } else {
        let (c, tail) = as_coeff_mul(&x, true);
        if c == Number::MINUS_ONE {
            (pow(&b, &tail), BigInt::from(-1))
        } else if !c.is_one() {
            let (p, q) = c.numer_denom().unwrap();
            let t = mul(&[Expr::num(Number::from_ratio(num_rational::BigRational::new(BigInt::one(), q))), tail]);
            (pow(&b, &t), p)
        } else {
            (f.clone(), BigInt::one())
        }
    };
    if e.is_negative() {
        (pow(&base, &Expr::int(-1)), -e)
    } else {
        (base, e)
    }
}

/// The reference `_not_a_coeff` numbers stay generators.
fn is_coeff(f: &Expr) -> bool {
    match f.kind() {
        Kind::Num(n) => !n.is_nan(),
        Kind::Const(Constant::I) => true,
        _ => false,
    }
}

/// `Poly(e)` expands its argument first.
fn poly_expand(e: &Expr) -> Expr {
    let e = crate::simptools::expand_hints(e, true, true, true, true, true, true);
    mexpand(&e)
}

/// The generators of `e.as_poly()` (as a set, in first-seen order).
fn auto_gens(e: &Expr) -> Vec<Expr> {
    let e = poly_expand(e);
    let mut gens: Vec<Expr> = Vec::new();
    for t in make_args_add(&e) {
        for f in make_args_mul(&t) {
            if is_coeff(&f) {
                continue;
            }
            let (b, _) = decompose_power(&f);
            if !gens.contains(&b) {
                gens.push(b);
            }
        }
    }
    gens
}

/// `e.as_poly(*gens)` as `(monomial -> coefficient)` terms, None when a
/// factor contains a generator without being a power of one.
fn poly_terms(e: &Expr, gens: &[Expr]) -> Option<Vec<(Vec<BigInt>, Expr)>> {
    let e = poly_expand(e);
    let mut out: Vec<(Vec<BigInt>, Vec<Expr>)> = Vec::new();
    for t in make_args_add(&e) {
        let mut coeff = Vec::new();
        let mut monom = vec![BigInt::zero(); gens.len()];
        for f in make_args_mul(&t) {
            if matches!(f.kind(), Kind::Num(_)) {
                coeff.push(f);
                continue;
            }
            let (b, k) = decompose_power(&f);
            if let Some(i) = gens.iter().position(|g| *g == b) {
                monom[i] = k;
            } else if !has_free(&f, gens) {
                coeff.push(f);
            } else {
                return None;
            }
        }
        let c = mul(&coeff);
        match out.iter_mut().find(|(m, _)| *m == monom) {
            Some((_, v)) => v.push(c),
            None => out.push((monom, vec![c])),
        }
    }
    Some(out.into_iter().map(|(m, v)| (m, add(&v))).filter(|(_, c)| !c.is_zero()).collect())
}

fn poly_as_expr(terms: &[(Vec<BigInt>, Expr)], gens: &[Expr]) -> Expr {
    let mut v = Vec::new();
    for (m, c) in terms {
        let mut f = vec![c.clone()];
        for (g, k) in gens.iter().zip(m) {
            if !k.is_zero() {
                f.push(pow(g, &Expr::num(Number::from_bigint(k.clone()))));
            }
        }
        v.push(mul(&f));
    }
    add(&v)
}

/// `Poly(e, p).degree(p)`.
fn degree_in<E>(e: &Expr, p: &Expr) -> R<BigInt, E> {
    match poly_terms(e, std::slice::from_ref(p)) {
        Some(t) => Ok(t.iter().map(|(m, _)| m[0].clone()).max().unwrap_or_else(|| BigInt::from(-1))),
        None => Err(UFail::Poly(format!("{} contains an element of the set of generators.", crate::printer::to_str(e)))),
    }
}

/// `_rads_bases_lcm(poly)` over the generators of a poly.
struct RBL {
    rads: Vec<Expr>,
    bases: Vec<Expr>,
    lcm: BigInt,
}

fn rads_bases_lcm(gens: &[Expr]) -> RBL {
    let mut lcm = BigInt::one();
    let mut rads = Vec::new();
    let mut bases: Vec<Expr> = Vec::new();
    for g in gens {
        let q = q_of(g);
        if !q.is_one() {
            if !rads.contains(g) {
                rads.push(g.clone());
            }
            lcm = lcm.lcm(&q);
            let b = as_base_exp(g).0;
            let b = match g.kind() {
                Kind::Pow(b, _) => b.clone(),
                _ => b,
            };
            if !bases.contains(&b) {
                bases.push(b);
            }
        }
    }
    RBL { rads, bases, lcm }
}

fn set_eq(a: &[Expr], b: &[Expr]) -> bool {
    a.len() == b.len() && a.iter().all(|x| b.contains(x))
}

fn ordered(mut v: Vec<Expr>) -> Vec<Expr> {
    crate::logic::ordered(&mut v);
    v
}

fn is_power_of_two(n: &BigInt) -> bool {
    n.is_positive() && (n & (n - BigInt::one())).is_zero()
}

fn big(n: &BigInt) -> Expr {
    Expr::num(Number::from_bigint(n.clone()))
}

fn covsym() -> Expr {
    Expr::dummy("p", crate::symbol::Assumptions::NONNEGATIVE)
}

fn sym_name(p: &Expr) -> String {
    match p.kind() {
        Kind::Sym(s) => s.name().to_string(),
        _ => "p".to_string(),
    }
}

/// `_canonical(eq, cov)`.
fn canonical(eq: &Expr, cov: &[Expr]) -> Unrad {
    let mut eq = eq.clone();
    let mut cov = cov.to_vec();
    if cov.len() == 2 {
        let p = cov[0].clone();
        let np = Expr::dummy(&sym_name(&p), crate::symbol::Assumptions::NONE);
        let mut rep = FxHashMap::default();
        rep.insert(p, np.clone());
        eq = xreplace(&eq, &rep);
        cov = vec![np, xreplace(&cov[1], &rep)];
    }
    eq = crate::simptools::exprtools::factor_terms(&mexpand(&numer(&eq)), false, true, false, true);
    if eq.is_mul() {
        let mut args = Vec::new();
        for f in args_of(&eq) {
            if is_number(&f) {
                continue;
            }
            if is_sym_pow(&f) {
                args.push(as_base_exp_pow(&f).0);
            } else {
                args.push(f);
            }
        }
        eq = mul(&args);
    }
    let mut margs = make_args_mul(&eq);
    let mut changed = false;
    for m in margs.iter_mut() {
        if crate::decompose::could_extract_minus_sign(m) {
            *m = neg(m);
            changed = true;
        }
    }
    if changed {
        eq = if margs.len() == 1 { margs.pop().unwrap() } else { crate::simptools::unevaluated_mul(&margs) };
    }
    (eq, cov)
}

/// `(base, exp)` of a Pow node (not the Rational-base rewrite).
fn as_base_exp_pow(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        _ => (e.clone(), Expr::one()),
    }
}

// ---------------------------------------------------------------- unrad

/// The flags that the reference passes down its recursion.
struct State {
    cov: Vec<Expr>,
    nwas: Option<usize>,
    rpt: usize,
    reverse: bool,
    custom_take: bool,
    /// The symbols of the outermost frame, which its `_take` closure reads.
    take_syms: Option<Vec<Expr>>,
}

fn take<H: UnradHost>(h: &mut H, st: &State, d: &Expr, syms: &[Expr]) -> R<bool, H::Error> {
    if st.custom_take {
        if let Some(r) = h.take(d) {
            return r;
        }
    }
    let syms = st.take_syms.as_deref().unwrap_or(syms);
    for p in make_args_mul(d) {
        if !is_sym_pow(&p) {
            continue;
        }
        if q_of(&p).is_one() {
            continue;
        }
        if intersects(&p, syms) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `cov[:] = ...` of `_cov(p, e)`.
fn set_cov<H: UnradHost>(h: &mut H, cov: &mut Vec<Expr>, p: &Expr, e: &Expr) -> R<(), H::Error> {
    if cov.len() == 2 {
        let (oldp, olde) = (cov[0].clone(), cov[1].clone());
        let d = degree_in(e, p)?;
        if d == BigInt::from(1) || d == BigInt::from(2) {
            let s = h.vsolve(e, p)?;
            let Some(s0) = s.first() else { return Err(UFail::Index("list index out of range".into())) };
            *cov = vec![p.clone(), subs1(&olde, &oldp, s0)];
        } else {
            return ni("");
        }
    } else {
        *cov = vec![p.clone(), e.clone()];
    }
    Ok(())
}

fn norm2(a: &Expr, b: &Expr) -> Expr {
    add(&[pow(a, &Expr::int(2)), pow(b, &Expr::int(2)), mul(&[Expr::int(2), a.clone(), b.clone()])])
}

/// The unradicalized form of `a*root(A, 3) + b*root(B, 3) + d = c`.
fn zz(a: &Expr, aa: &Expr, b: &Expr, bb: &Expr, c: &Expr, d: &Expr) -> Expr {
    let p = |e: &Expr, k: i64| pow(e, &Expr::int(k));
    let t = |k: i64, f: &[Expr]| {
        let mut v = vec![Expr::int(k)];
        v.extend_from_slice(f);
        mul(&v)
    };
    let (a3, a6, a9) = (p(a, 3), p(a, 6), p(a, 9));
    let (b3, b6, b9) = (p(b, 3), p(b, 6), p(b, 9));
    let terms = vec![
        t(1, &[p(aa, 3), a9.clone()]), t(3, &[p(aa, 2), bb.clone(), a6.clone(), b3.clone()]),
        t(-3, &[p(aa, 2), a6.clone(), p(c, 3)]), t(9, &[p(aa, 2), a6.clone(), p(c, 2), d.clone()]),
        t(-9, &[p(aa, 2), a6.clone(), c.clone(), p(d, 2)]), t(3, &[p(aa, 2), a6.clone(), p(d, 3)]),
        t(3, &[aa.clone(), p(bb, 2), a3.clone(), b6.clone()]), t(21, &[aa.clone(), bb.clone(), a3.clone(), b3.clone(), p(c, 3)]),
        t(-63, &[aa.clone(), bb.clone(), a3.clone(), b3.clone(), p(c, 2), d.clone()]),
        t(63, &[aa.clone(), bb.clone(), a3.clone(), b3.clone(), c.clone(), p(d, 2)]),
        t(-21, &[aa.clone(), bb.clone(), a3.clone(), b3.clone(), p(d, 3)]), t(3, &[aa.clone(), a3.clone(), p(c, 6)]),
        t(-18, &[aa.clone(), a3.clone(), p(c, 5), d.clone()]), t(45, &[aa.clone(), a3.clone(), p(c, 4), p(d, 2)]),
        t(-60, &[aa.clone(), a3.clone(), p(c, 3), p(d, 3)]), t(45, &[aa.clone(), a3.clone(), p(c, 2), p(d, 4)]),
        t(-18, &[aa.clone(), a3.clone(), c.clone(), p(d, 5)]), t(3, &[aa.clone(), a3.clone(), p(d, 6)]),
        t(1, &[p(bb, 3), b9.clone()]), t(-3, &[p(bb, 2), b6.clone(), p(c, 3)]),
        t(9, &[p(bb, 2), b6.clone(), p(c, 2), d.clone()]), t(-9, &[p(bb, 2), b6.clone(), c.clone(), p(d, 2)]),
        t(3, &[p(bb, 2), b6.clone(), p(d, 3)]), t(3, &[bb.clone(), b3.clone(), p(c, 6)]),
        t(-18, &[bb.clone(), b3.clone(), p(c, 5), d.clone()]), t(45, &[bb.clone(), b3.clone(), p(c, 4), p(d, 2)]),
        t(-60, &[bb.clone(), b3.clone(), p(c, 3), p(d, 3)]), t(45, &[bb.clone(), b3.clone(), p(c, 2), p(d, 4)]),
        t(-18, &[bb.clone(), b3.clone(), c.clone(), p(d, 5)]), t(3, &[bb.clone(), b3.clone(), p(d, 6)]),
        t(-1, &[p(c, 9)]), t(9, &[p(c, 8), d.clone()]), t(-36, &[p(c, 7), p(d, 2)]), t(84, &[p(c, 6), p(d, 3)]),
        t(-126, &[p(c, 5), p(d, 4)]), t(126, &[p(c, 4), p(d, 5)]), t(-84, &[p(c, 3), p(d, 6)]),
        t(36, &[p(c, 2), p(d, 7)]), t(-9, &[c.clone(), p(d, 8)]), t(1, &[p(d, 9)]),
    ];
    mul(&[sub(c, d), add(&terms)])
}

/// SymPy's `unrad(eq, *syms, **flags)` for an Expr `eq` (an Equality is
/// moved to one side by the caller).
pub fn unrad<H: UnradHost>(h: &mut H, eq: &Expr, syms: &[Expr], reverse: bool, custom_take: bool) -> R<Option<Unrad>, H::Error> {
    let mut st = State { cov: Vec::new(), nwas: None, rpt: 0, reverse, custom_take, take_syms: None };
    unrad_rec(h, eq, syms, &mut st)
}

/// `unrad(eq, *syms, cov=cov, n=n, rpt=rpt)`: the reference recursion flags
/// given by the caller. Also returns the final content of the `cov` list,
/// which the reference updates in place (`cov[:] = ...`).
#[allow(clippy::too_many_arguments)]
pub fn unrad_with_flags<H: UnradHost>(
    h: &mut H, eq: &Expr, syms: &[Expr], reverse: bool, custom_take: bool, cov: Vec<Expr>, nwas: Option<usize>, rpt: usize,
) -> R<(Option<Unrad>, Vec<Expr>), H::Error> {
    let mut st = State { cov, nwas, rpt, reverse, custom_take, take_syms: None };
    let r = unrad_rec(h, eq, syms, &mut st)?;
    Ok((r, st.cov))
}

fn unrad_rec<H: UnradHost>(h: &mut H, eq0: &Expr, syms_in: &[Expr], st: &mut State) -> R<Option<Unrad>, H::Error> {
    let cov_in = st.cov.clone();
    let nwas = st.nwas;
    let mut rpt = st.rpt;
    // preconditioning
    let mut eq = crate::simptools::powsimp::powdenest(&crate::simptools::exprtools::factor_terms(eq0, true, true, false, true), false, false);
    eq = numer(&eq);
    eq = mexpand(&eq);
    if is_number(&eq) {
        return Ok(None);
    }
    let mut syms: Vec<Expr> = Vec::new();
    for s in syms_in {
        if !syms.contains(s) {
            syms.push(s.clone());
        }
    }
    if syms.is_empty() {
        syms = free_symbols(&eq);
    }
    let all_gens = auto_gens(&eq);
    let mut gens = Vec::new();
    for g in &all_gens {
        if take(h, st, g, &syms)? {
            gens.push(g.clone());
        }
    }
    if gens.is_empty() {
        return Ok(None);
    }
    let Some(pterms) = poly_terms(&eq, &gens) else { return Ok(None) };
    // `g.exp` of the reference: a generator taken by a custom `_take` that
    // is not a power has no `exp`
    for g in &gens {
        let x = match g.kind() {
            Kind::Pow(_, x) => x.clone(),
            Kind::Func(Func::Exp, a) => a[0].clone(),
            _ => return Err(UFail::Attribute(format!("'{}' object has no attribute 'exp'", g.class_name()))),
        };
        if has_any(&x, &syms) {
            return Ok(None);
        }
    }
    let pgens = gens.clone();
    let RBL { rads, bases, lcm } = rads_bases_lcm(&pgens);
    let cs = covsym();
    let mut newsyms: Vec<Expr> = Vec::new();
    for r in &rads {
        for s in free_symbols(r) {
            if syms.contains(&s) && !newsyms.contains(&s) {
                newsyms.push(s);
            }
        }
    }
    if !set_eq(&newsyms, &syms) {
        syms = newsyms;
    }
    let mut syms_sorted = syms.clone();
    crate::logic::ordered(&mut syms_sorted);
    let syms = syms_sorted;
    // get terms together that have common generators
    let pexpr = poly_as_expr(&pterms, &pgens);
    let mut groups: Vec<(Vec<usize>, Vec<Expr>)> = vec![(Vec::new(), Vec::new())];
    for t in make_args_add(&pexpr) {
        let key = if take(h, st, &t, &syms)? {
            let tg = auto_gens(&t);
            let mut k: Vec<usize> = rads.iter().enumerate().filter(|(_, r)| tg.contains(r)).map(|(i, _)| i).collect();
            k.sort();
            k
        } else {
            Vec::new()
        };
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, v)) => v.push(t),
            None => groups.push((key, vec![t])),
        }
    }
    let others = add(&groups.remove(0).1);
    let mut rterms: Vec<Expr> = groups.into_iter().map(|(_, v)| add(&v)).collect();
    rterms = ordered(rterms);
    rterms.reverse();

    let mut ok = false;
    let depth = crate::simptools::sqrtdenest::sqrt_depth(&eq);
    let lcm_e = big(&lcm);
    let two = BigInt::from(2);
    if rterms.len() == 1 && !(rterms[0].is_add() && lcm > two) {
        eq = sub(&pow(&rterms[0], &lcm_e), &pow(&neg(&others), &lcm_e));
        ok = true;
    } else {
        if rterms.len() == 1 && rterms[0].is_add() {
            rterms = args_of(&rterms[0]);
        }
        if bases.len() == 1 {
            let b = bases[0].clone();
            let mut xs = if syms.len() > 1 { free_symbols(&b) } else { syms.clone() };
            crate::logic::ordered(&mut xs);
            let x = xs[0].clone();
            let attempt: R<Unrad, H::Error> = (|| {
                let inv = h.vsolve(&sub(&pow(&cs, &lcm_e), &b), &x)?;
                if inv.is_empty() {
                    return ni("");
                }
                let e = subs1(&subs1(&pexpr, &b, &pow(&cs, &lcm_e)), &x, &inv[0]);
                let mut cov = cov_in.clone();
                set_cov(h, &mut cov, &cs, &sub(&pow(&cs, &lcm_e), &b))?;
                st.cov = cov.clone();
                Ok(canonical(&e, &cov))
            })();
            match attempt {
                Ok(r) => return Ok(Some(r)),
                Err(UFail::NotImplemented(_)) => {}
                Err(e) => return Err(e),
            }
        }
        if rterms.len() == 2 {
            if others.is_zero() {
                eq = sub(&pow(&rterms[0], &lcm_e), &pow(&neg(&rterms[1]), &lcm_e));
                ok = true;
            } else if !is_power_of_two(&lcm) {
                let (mut r0, mut r1) = (rterms[0].clone(), rterms[1].clone());
                if st.reverse {
                    std::mem::swap(&mut r0, &mut r1);
                }
                let mut i0 = rads_bases_lcm(&auto_gens(&r0));
                let mut i1 = rads_bases_lcm(&auto_gens(&r1));
                'outer: for rev in 0..2 {
                    if rev == 1 {
                        std::mem::swap(&mut i0, &mut i1);
                        std::mem::swap(&mut r0, &mut r1);
                    }
                    let rads1 = mul(&i1.rads);
                    let lcm1 = big(&i1.lcm);
                    let t1 = pow(&rads1, &lcm1);
                    let c = sub(&pow(&cs, &lcm1), &t1);
                    for x in &syms {
                        let attempt: R<(), H::Error> = (|| {
                            let sol = h.vsolve(&c, x)?;
                            if sol.is_empty() {
                                return ni("");
                            }
                            let neweq = add(&[subs1(&r0, x, &sol[0]), div(&mul(&[cs.clone(), r1.clone()]), &rads1), others.clone()]);
                            let tmp = unrad(h, &neweq, std::slice::from_ref(&cs), false, false)?;
                            let mut cov = st.cov.clone();
                            match tmp {
                                Some((e, newcov)) => {
                                    eq = e;
                                    if newcov.len() == 2 {
                                        let (newp, newc) = (newcov[0].clone(), newcov[1].clone());
                                        let s = h.vsolve(&newc, &cs)?;
                                        let Some(s0) = s.first() else { return Err(UFail::Index("list index out of range".into())) };
                                        set_cov(h, &mut cov, &newp, &subs1(&c, &cs, s0))?;
                                    } else {
                                        set_cov(h, &mut cov, &cs, &c)?;
                                    }
                                }
                                None => {
                                    eq = neweq;
                                    set_cov(h, &mut cov, &cs, &c)?;
                                }
                            }
                            st.cov = cov;
                            Ok(())
                        })();
                        match attempt {
                            Ok(()) => {
                                ok = true;
                                break;
                            }
                            Err(UFail::NotImplemented(_)) => {
                                if rev == 1 {
                                    return ni("no successful change of variable found");
                                }
                            }
                            Err(e) => return Err(e),
                        }
                    }
                    if ok {
                        break 'outer;
                    }
                }
            }
        } else if rterms.len() == 3 {
            let mut info: Vec<RBL> = rterms.iter().map(|t| rads_bases_lcm(&auto_gens(t))).collect();
            let three = BigInt::from(3);
            if info[0].lcm != three {
                let i = info.remove(0);
                info.push(i);
                let r = rterms.remove(0);
                rterms.push(r);
            } else if info[1].lcm != three {
                let i = info.remove(1);
                info.push(i);
                let r = rterms.remove(1);
                rterms.push(r);
            }
            if info[0].lcm == three && info[1].lcm == three {
                if !set_eq(&info[1].bases, &info[2].bases) {
                    info.swap(0, 1);
                    rterms.swap(0, 1);
                }
                if set_eq(&info[1].bases, &info[2].bases) {
                    eq = add(&[pow(&rterms[0], &Expr::int(3)), pow(&add(&[rterms[1].clone(), rterms[2].clone(), others.clone()]), &Expr::int(3))]);
                    ok = true;
                } else if info[2].lcm < BigInt::from(5) {
                    let tt = |i: usize| -> R<(Expr, Expr), H::Error> {
                        let b = mul(&info[i].rads);
                        let q = crate::poly::api::cancel(&div(&rterms[i], &b)).map_err(|e| UFail::Poly(format!("{e:?}")))?;
                        Ok((q, mul(&info[i].bases)))
                    };
                    let (aa, a_) = tt(0)?;
                    let (bb, b_) = tt(1)?;
                    let cc = neg(&rterms[2]);
                    eq = zz(&aa, &a_, &bb, &b_, &cc, &others);
                    ok = true;
                }
            }
        }
        // handle power-of-2 cases
        if !ok && is_power_of_two(&lcm) && ((others.is_zero() && rterms.len() == 4) || rterms.len() < 4) {
            match rterms.len() {
                4 => {
                    eq = sub(&norm2(&rterms[0], &rterms[1]), &norm2(&rterms[2], &rterms[3]));
                    ok = true;
                }
                3 => {
                    eq = sub(&norm2(&rterms[1], &rterms[2]), &norm2(&rterms[0], &others));
                    ok = true;
                }
                2 => {
                    eq = sub(&pow(&rterms[0], &Expr::int(2)), &norm2(&rterms[1], &others));
                    ok = true;
                }
                _ => {}
            }
        }
    }
    let new_depth = if ok { crate::simptools::sqrtdenest::sqrt_depth(&eq) } else { depth };
    rpt += 1;
    if !ok || (nwas == Some(rterms.len()) && new_depth == depth && rpt > 3) {
        return ni("Cannot remove all radicals");
    }
    st.nwas = Some(rterms.len());
    st.rpt = rpt;
    if st.take_syms.is_none() && !st.custom_take {
        st.take_syms = Some(syms.clone());
    }
    let neq = unrad_rec(h, &eq, &syms, st)?;
    let mut cov = st.cov.clone();
    if let Some((e, c)) = neq {
        eq = e;
        cov = c;
    }
    Ok(Some(canonical(&eq, &cov)))
}
