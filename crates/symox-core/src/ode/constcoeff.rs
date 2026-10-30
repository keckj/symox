//! Characteristic-equation solutions of linear ODEs with constant
//! coefficients and of Cauchy-Euler equations, as the reference (SymPy 1.14,
//! `sympy/solvers/ode/nonhomogeneous.py`) builds them:
//!
//! * `_get_const_characteristic_eq_sols`: the characteristic polynomial in a
//!   fresh dummy, its roots with multiplicities in the order of
//!   `roots(chareq, multiple=True)` (`ordered`: node count, then
//!   `default_sort_key`), or the indexed roots `rootof(chareq, k)` when the
//!   radicals are incomplete, and the fundamental solutions with
//!   `exp`/`sin`/`cos` together with the `collectterms` list;
//! * `_get_simplified_sol`: the collection of the general solution on the
//!   `collectterms` patterns followed by `powsimp`;
//! * `_get_euler_characteristic_eq_sols`: the indicial polynomial, the roots
//!   `rootof(chareq, k)`, the general solution with the numbered constants
//!   and the fundamental solutions.
//!
//! Every function returns `None` whenever the inputs leave the subset on which
//! the native construction is known to reproduce the reference expression;
//! the caller then runs the reference body unchanged.

use std::cmp::Ordering;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::ops::{add, mul, pow};
use crate::poly::collect::{self, Collection, Options};
use crate::poly::int::Int;
use crate::poly::polyobj::Poly;
use crate::poly::roots as proots;
use crate::poly::rootof;
use crate::symbol::Assumptions;

/// `(i, reroot, imroot)` entries of the reference `collectterms`.
pub type CollectTerm = (u32, Expr, Expr);

/// The fundamental solutions and the `collectterms` of a characteristic
/// equation.
#[derive(Clone, Debug)]
pub struct CharSols {
    pub gensols: Vec<Expr>,
    pub collectterms: Vec<CollectTerm>,
}

fn m(a: &Expr, b: &Expr) -> Expr {
    mul(&[a.clone(), b.clone()])
}

fn f1(f: Func, a: Expr) -> Expr {
    func(f, vec![a])
}

fn ipow(x: &Expr, i: u32) -> Expr {
    pow(x, &Expr::int(i as i64))
}

fn has_atan2(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |a| {
        if matches!(a.kind(), Kind::Func(Func::Atan2, _)) {
            found = true;
        }
    });
    found
}

/// A coefficient on which `Poly(chareq, x)` is a polynomial over QQ or a
/// polynomial ring QQ[a, b, ...] in the free symbols: rational numbers and
/// symbols combined by sums, products and non-negative integer powers.
fn admissible_coefficient(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_rational(),
        Kind::Sym(s) => !s.is_wild(),
        Kind::Add(c, terms) => c.is_rational() && terms.iter().all(|(t, k)| k.is_rational() && admissible_coefficient(t)),
        Kind::Mul(c, factors) => {
            c.is_rational()
                && factors.iter().all(|(b, x)| {
                    matches!(x.as_number(), Some(n) if n.is_integer() && !n.is_negative()) && admissible_coefficient(b)
                })
        }
        Kind::Pow(b, x) => matches!(x.as_number(), Some(n) if n.is_integer() && !n.is_negative()) && admissible_coefficient(b),
        _ => false,
    }
}

/// `ordered(roots)` of the reference `roots(..., multiple=True)`, each root
/// repeated by its multiplicity.
fn ordered_multiple(rs: Vec<(Expr, u32)>) -> Vec<Expr> {
    let mut rs: Vec<(usize, Expr, u32)> = rs.into_iter().map(|(r, k)| (crate::simptools::nodes(&r), r, k)).collect();
    rs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| crate::sortkey::cmp_default(&a.1, &b.1)));
    let mut out = Vec::new();
    for (_, r, k) in rs {
        for _ in 0..k {
            out.push(r.clone());
        }
    }
    out
}

/// The rational coefficients, leading first, of a univariate polynomial
/// over QQ (no coefficient generators).
fn rational_coeffs(p: &Poly) -> Option<Vec<num_rational::BigRational>> {
    if p.ngens != 1 || p.has_domain_gens() || p.is_float() {
        return None;
    }
    p.all_coeffs().ok()?.iter().map(|c| c.as_number()?.to_ratio()).collect()
}

/// `clear_denoms` and `primitive` of a polynomial over QQ, leading first.
fn primitive_integer(c: &[num_rational::BigRational]) -> Option<Vec<BigInt>> {
    let mut lcm = BigInt::one();
    for q in c {
        lcm = lcm.lcm(q.denom());
    }
    let ints: Vec<BigInt> = c.iter().map(|q| q.numer() * (&lcm / q.denom())).collect();
    let mut g = BigInt::zero();
    for v in &ints {
        g = g.gcd(v);
    }
    if g.is_zero() {
        return None;
    }
    Some(ints.into_iter().map(|v| v / &g).collect())
}

/// `[rootof(chareq, k) for k in range(degree)]` for an irreducible chareq
/// over QQ of degree at least three whose roots are not given by the
/// reference's trivial (linear, quadratic, binomial) formulas.
fn indexed_roots_irreducible(p: &Poly, gen: &Expr) -> Option<Vec<Expr>> {
    let c = rational_coeffs(p)?;
    let degree = c.len().checked_sub(1)?;
    if degree < 3 {
        return None;
    }
    let exprs: Vec<Expr> = c.iter().map(|q| Expr::num(crate::number::Number::from_ratio(q.clone()))).collect();
    if proots::roots_trivial(&exprs).is_some() {
        return None;
    }
    let (_, factors) = p.factor_list();
    if factors.len() != 1 || factors[0].1 != 1 || factors[0].0.degree(0)? as usize != degree {
        return None;
    }
    let ints = primitive_integer(&c)?;
    if !ints[0].is_positive() {
        return None;
    }
    let ints: Vec<Int> = ints.into_iter().map(Int::from_big).collect();
    let roots = rootof::all_roots_irreducible_in(&ints, gen)?;
    (roots.len() == degree).then_some(roots)
}

/// The reference `_get_const_characteristic_eq_sols(r, func, order)`.
/// `coeffs` holds the integer keys `i >= 0` of `r` with their values.
pub fn const_characteristic_sols(coeffs: &[(u32, Expr)], x: &Expr, order: usize) -> Option<CharSols> {
    if coeffs.iter().any(|(_, c)| !admissible_coefficient(c)) {
        return None;
    }
    let symbol = Expr::dummy("x", Assumptions::NONE);
    let chareq = add(&coeffs.iter().map(|(i, c)| m(c, &ipow(&symbol, *i))).collect::<Vec<_>>());
    let p = Poly::from_expr(&chareq, std::slice::from_ref(&symbol)).ok()?;
    if p.ngens != 1 || p.is_float() || p.gens.first() != Some(&symbol) {
        return None;
    }
    let found = proots::with_radical_options(proots::RadicalOptions::default(), || proots::roots_poly(&p)).ok()?;
    let mut chareqroots = ordered_multiple(found);
    if chareqroots.len() != order {
        chareqroots = indexed_roots_irreducible(&p, &symbol)?;
        if chareqroots.len() != order {
            return None;
        }
    }
    let all = p.all_coeffs().ok()?;
    let chareq_is_complex = !all.iter().all(|c| crate::assumptions::is_real(c) == Some(true));

    // Counter(chareqroots): distinct roots in first-appearance order.
    let mut counts: Vec<(Expr, u32)> = Vec::new();
    for r in &chareqroots {
        match counts.iter_mut().find(|(e, _)| e == r) {
            Some(entry) => entry.1 += 1,
            None => counts.push((r.clone(), 1)),
        }
    }
    let mut collectterms: Vec<CollectTerm> = Vec::new();
    let mut gensols: Vec<Expr> = Vec::new();
    let mut conjugate_roots: Vec<Expr> = Vec::new();
    for root in &chareqroots {
        let Some(pos) = counts.iter().position(|(e, _)| e == root) else {
            continue;
        };
        let (_, multiplicity) = counts.remove(pos);
        for i in 0..multiplicity {
            let xi = ipow(x, i);
            if chareq_is_complex {
                gensols.push(m(&xi, &f1(Func::Exp, m(root, x))));
                collectterms.insert(0, (i, root.clone(), Expr::zero()));
                continue;
            }
            let reroot = f1(Func::Re, root.clone());
            let imroot = f1(Func::Im, root.clone());
            if has_atan2(&imroot) && has_atan2(&reroot) {
                gensols.push(m(&xi, &f1(Func::Exp, m(root, x))));
                collectterms.insert(0, (i, root.clone(), Expr::zero()));
            } else {
                if conjugate_roots.contains(root) {
                    collectterms.insert(0, (i, reroot, imroot));
                    continue;
                }
                if imroot.is_zero() {
                    gensols.push(m(&xi, &f1(Func::Exp, m(&reroot, x))));
                    collectterms.insert(0, (i, reroot, Expr::zero()));
                    continue;
                }
                conjugate_roots.push(f1(Func::Conjugate, root.clone()));
                let e = m(&xi, &f1(Func::Exp, m(&reroot, x)));
                gensols.push(m(&e, &f1(Func::Sin, m(&f1(Func::Abs, imroot.clone()), x))));
                gensols.push(m(&e, &f1(Func::Cos, m(&imroot, x))));
                collectterms.insert(0, (i, reroot, imroot));
            }
        }
    }
    Some(CharSols { gensols, collectterms })
}

fn cmp_term(a: &CollectTerm, b: &CollectTerm) -> Ordering {
    let key = |t: &CollectTerm| [Expr::int(t.0 as i64), t.1.clone(), t.2.clone()];
    let (ka, kb) = (key(a), key(b));
    for (u, v) in ka.iter().zip(kb.iter()) {
        let o = crate::sortkey::cmp_default(u, v);
        if o != Ordering::Equal {
            return o;
        }
    }
    Ordering::Equal
}

fn collect_on(e: &Expr, pattern: &Expr) -> Option<Expr> {
    // the reference `collect(sol, pattern)`: exact=False, which does not
    // discover compound generators (that is exact=None)
    let options = Options { exact: false, discover: false, distribute_order_term: true };
    match collect::run::<()>(e, std::slice::from_ref(pattern), options, true, None).ok()? {
        Collection::Expression(v) => Some(v),
        Collection::Dictionary(_) => None,
    }
}

/// The reference `_get_simplified_sol([Eq(f(x), rhs)], f(x), collectterms)`,
/// returning the simplified right-hand side.
pub fn simplified_solution(rhs: &Expr, x: &Expr, collectterms: &[CollectTerm]) -> Option<Expr> {
    let mut terms = collectterms.to_vec();
    // collectterms.sort(key=default_sort_key); collectterms.reverse()
    terms.sort_by(cmp_term);
    terms.reverse();
    let mut sol = crate::simp::expand_mul_deep(rhs, true);
    for (i, reroot, imroot) in &terms {
        let base = m(&ipow(x, *i), &f1(Func::Exp, m(reroot, x)));
        sol = collect_on(&sol, &m(&base, &f1(Func::Sin, m(&f1(Func::Abs, imroot.clone()), x))))?;
        sol = collect_on(&sol, &m(&base, &f1(Func::Cos, m(imroot, x))))?;
    }
    for (i, reroot, _) in &terms {
        sol = collect_on(&sol, &m(&ipow(x, *i), &f1(Func::Exp, m(reroot, x))))?;
    }
    Some(crate::simptools::powsimp::powsimp(&sol, &crate::simptools::powsimp::PowsimpOptions::default()))
}

/// The roots `[rootof(chareq, k) for k in range(degree)]` of a polynomial
/// over QQ, when the reference returns explicit numbers: the trivial
/// (linear, quadratic, binomial) formulas, or only rational roots, which
/// come in increasing order repeated by multiplicity.
fn rootof_explicit(p: &Poly) -> Option<Vec<Expr>> {
    let c = rational_coeffs(p)?;
    let degree = c.len().checked_sub(1)?;
    if degree == 0 {
        return None;
    }
    let exprs: Vec<Expr> = c.iter().map(|q| Expr::num(crate::number::Number::from_ratio(q.clone()))).collect();
    if let Some(rs) = proots::roots_trivial(&exprs) {
        return (rs.len() == degree).then_some(rs);
    }
    let (_, factors) = p.factor_list();
    let mut rs: Vec<(num_rational::BigRational, u32)> = Vec::new();
    for (f, k) in &factors {
        let fc = rational_coeffs(f)?;
        if fc.len() != 2 || fc[0].is_zero() {
            return None;
        }
        rs.push((-(&fc[1] / &fc[0]), *k));
    }
    rs.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out = Vec::new();
    for (r, k) in rs {
        let e = Expr::num(crate::number::Number::from_ratio(r));
        for _ in 0..k {
            out.push(e.clone());
        }
    }
    (out.len() == degree).then_some(out)
}

/// The reference `_get_euler_characteristic_eq_sols(eq, func, match_obj)`:
/// the right-hand side of the general solution and the fundamental
/// solutions. `coeffs` holds the integer keys `i >= 0` of `match_obj`;
/// `constants` the `2*degree` numbered constants in increasing order.
pub fn euler_characteristic_sols(coeffs: &[(u32, Expr)], x: &Expr, constants: &[Expr]) -> Option<(Expr, Vec<Expr>)> {
    let symbol = Expr::dummy("x", Assumptions::NONE);
    let xs = pow(x, &symbol);
    let xms = pow(x, &crate::ops::neg(&symbol));
    let mut parts = Vec::new();
    for (i, c) in coeffs {
        let mut d = xs.clone();
        for _ in 0..*i {
            d = crate::diff::diff(&d, x);
        }
        parts.push(crate::expand::expand(&mul(&[c.clone(), d, xms.clone()])));
    }
    let chareq = add(&parts);
    let p = Poly::from_expr(&chareq, std::slice::from_ref(&symbol)).ok()?;
    if p.ngens != 1 || p.gens.first() != Some(&symbol) {
        return None;
    }
    let degree = p.degree(0)? as usize;
    if constants.len() != 2 * degree {
        return None;
    }
    let chareqroots = rootof_explicit(&p)?;
    let mut counts: Vec<(Expr, u32)> = Vec::new();
    for r in &chareqroots {
        match counts.iter_mut().find(|(e, _)| e == r) {
            Some(entry) => entry.1 += 1,
            None => counts.push((r.clone(), 1)),
        }
    }
    let ln = f1(Func::Log, x.clone());
    let mut next = constants.iter();
    let mut terms: Vec<Expr> = Vec::new();
    let mut collectterms: Vec<CollectTerm> = Vec::new();
    for (root, multiplicity) in &counts {
        for i in 0..*multiplicity {
            let lni = ipow(&ln, i);
            if rootof::is_rootof(root) {
                return None;
            } else if crate::assumptions::is_real(root) == Some(true) {
                terms.push(m(&m(&lni, &pow(x, root)), next.next()?));
                collectterms.insert(0, (i, root.clone(), Expr::zero()));
            } else {
                let reroot = f1(Func::Re, root.clone());
                let imroot = f1(Func::Im, root.clone());
                let c1 = next.next()?;
                let c2 = next.next()?;
                let s = m(c1, &f1(Func::Sin, m(&f1(Func::Abs, imroot.clone()), &ln)));
                let co = m(c2, &f1(Func::Cos, m(&imroot, &ln)));
                terms.push(m(&m(&lni, &pow(x, &reroot)), &add(&[s, co])));
                collectterms.insert(0, (i, reroot, imroot));
            }
        }
    }
    let mut gsol = Expr::zero();
    for t in terms {
        gsol = add(&[gsol, t]);
    }
    let mut gensols: Vec<Expr> = Vec::new();
    for (i, reroot, imroot) in &collectterms {
        let base = m(&ipow(&ln, *i), &pow(x, reroot));
        if imroot.is_zero() {
            gensols.push(base);
        } else {
            let sin_form = m(&base, &f1(Func::Sin, m(&f1(Func::Abs, imroot.clone()), &ln)));
            if gensols.contains(&sin_form) {
                gensols.push(m(&base, &f1(Func::Cos, m(imroot, &ln))));
            } else {
                gensols.push(sin_form);
            }
        }
    }
    Some((gsol, gensols))
}
