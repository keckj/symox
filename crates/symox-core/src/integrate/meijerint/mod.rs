//! Definite integration by Meijer G-functions, following SymPy's
//! `meijerint_definite` (Saxena's theorems for one and two G-functions).
//!
//! The integrand is rewritten as `C*x**s*G(eta*x**b)` or as a product of two
//! G-functions with tabulated rewrites; the integral over (0, oo) is then a
//! single G-function (or a product of gammas) valid under the tabulated
//! convergence conditions, which are returned alongside the value. Nothing is
//! reported without its condition, and a condition that simplifies to False
//! makes the method give up instead of returning a value.
//!
//! Indefinite integrals follow `meijerint_indefinite` (`indefinite`).
//!
//! When the table has no entry for a factor, `_rewrite_single` computes
//! its Mellin transform with Meijer G integration and reads the inverse
//! transform as a G-function (`recursive`).

use crate::expr::Constant;
use crate::functions as f;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};
use crate::special::polar::{polarify, unpolarify};
use crate::special::{tuple, tuple_items, SpecialFn};
use crate::{Expr, Func, Kind};

mod conds;
mod indefinite;
mod inversion;
mod recursive;
mod table;

pub(crate) use indefinite::meijerint_indefinite;
pub use indefinite::meijerint_indefinite_checked;

thread_local! {
    /// The first error the reference would raise out of a rewrite that
    /// the port abandons instead (class name, message).
    static PENDING: std::cell::RefCell<Option<(&'static str, String)>> = const { std::cell::RefCell::new(None) };
}

/// Records the error the reference raises at this point; the first one
/// wins, as the reference stops there.
pub(crate) fn set_pending(class: &'static str, message: String) {
    PENDING.with(|p| {
        let mut p = p.borrow_mut();
        if p.is_none() {
            *p = Some((class, message));
        }
    });
}

pub(crate) fn has_pending() -> bool {
    PENDING.with(|p| p.borrow().is_some())
}

pub(crate) fn take_pending() -> Option<(&'static str, String)> {
    PENDING.with(|p| p.borrow_mut().take())
}

/// The `_CoeffExpValueError` message of `_get_coeff_exp(expr, x)`.
pub(crate) fn coeff_exp_error(expr: &Expr, x: &Expr) -> String {
    let e = expand_power_base(&crate::simptools::powsimp::powsimp(expr, &Default::default()));
    let (_, m) = crate::decompose::as_coeff_mul(&e, std::slice::from_ref(x), false);
    match m.as_slice() {
        [m] if matches!(m.kind(), Kind::Pow(..)) => "expr not of form a*x**b".into(),
        _ => format!("expr not of form a*x**b: {expr}"),
    }
}
pub(crate) use inversion::meijerint_inversion;


/// A Meijer G-function `G^{m,n}_{p,q}(an, ao; bm, bo | z)`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct G {
    pub an: Vec<Expr>,
    pub ao: Vec<Expr>,
    pub bm: Vec<Expr>,
    pub bo: Vec<Expr>,
    pub z: Expr,
}

fn has_infinity(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity))
        || (!e.is_atom() && e.args().iter().any(has_infinity))
}

impl G {
    /// The `meijerg` constructor: parameters are unpolarified; infinite
    /// parameters and an `a_i - b_j` that is a positive integer are rejected
    /// (SymPy raises ValueError; the caller abandons the rewrite).
    pub fn new(an: Vec<Expr>, ao: Vec<Expr>, bm: Vec<Expr>, bo: Vec<Expr>, z: Expr) -> Option<G> {
        Self::try_new(an, ao, bm, bo, z).ok()
    }
    /// `new` with the reference ValueError message on a rejection.
    pub fn try_new(an: Vec<Expr>, ao: Vec<Expr>, bm: Vec<Expr>, bo: Vec<Expr>, z: Expr) -> Result<G, String> {
        // meijerg.__new__ stores each parameter list `ordered`.
        let up = |v: Vec<Expr>| crate::simptools::ordered(v.iter().map(|p| unpolarify(p, false)).collect());
        let (an, ao, bm, bo) = (up(an), up(ao), up(bm), up(bo));
        if an.iter().chain(&ao).chain(&bm).chain(&bo).any(has_infinity) {
            return Err("G-function parameters must be finite".into());
        }
        for a in &an {
            for b in &bm {
                let d = sub(a, b);
                if d.is_integer_number() && d.as_number().is_some_and(|n| n.is_positive()) {
                    return Err("no parameter a1, ..., an may differ from any b1, ..., bm by a positive integer".into());
                }
            }
        }
        Ok(G { an, ao, bm, bo, z })
    }
    pub fn ap(&self) -> Vec<Expr> { self.an.iter().chain(&self.ao).cloned().collect() }
    pub fn bq(&self) -> Vec<Expr> { self.bm.iter().chain(&self.bo).cloned().collect() }
    pub fn p(&self) -> usize { self.an.len() + self.ao.len() }
    pub fn q(&self) -> usize { self.bm.len() + self.bo.len() }
    pub fn nu(&self) -> Expr { sub(&add(&self.bq()), &add(&self.ap())) }
    pub fn delta(&self) -> Expr {
        let twice = 2 * (self.bm.len() + self.an.len()) as i64 - (self.p() + self.q()) as i64;
        Expr::rational(twice, 2)
    }
    pub fn to_expr(&self) -> Expr {
        crate::special::raw(SpecialFn::MeijerG, vec![
            tuple(vec![tuple(self.an.clone()), tuple(self.ao.clone())]),
            tuple(vec![tuple(self.bm.clone()), tuple(self.bo.clone())]),
            self.z.clone(),
        ])
    }
    pub fn from_expr(e: &Expr) -> Option<G> {
        let Kind::Func(Func::Special(SpecialFn::MeijerG), a) = e.kind() else { return None };
        let aa = tuple_items(&a[0])?;
        let bb = tuple_items(&a[1])?;
        Some(G {
            an: tuple_items(&aa[0])?.to_vec(),
            ao: tuple_items(&aa[1])?.to_vec(),
            bm: tuple_items(&bb[0])?.to_vec(),
            bo: tuple_items(&bb[1])?.to_vec(),
            z: a[2].clone(),
        })
    }
    fn with_z(&self, z: Expr) -> Option<G> {
        G::new(self.an.clone(), self.ao.clone(), self.bm.clone(), self.bo.clone(), z)
    }
    /// `get_period`: the period of the function in `log(z)`.
    pub fn period(&self) -> Expr {
        fn compute(l: &[Expr]) -> Option<i64> {
            let mut lcm: i64 = 1;
            for (i, b) in l.iter().enumerate() {
                let n = b.as_number().filter(|n| n.is_rational())?;
                for c in &l[i + 1..] {
                    if sub(b, c).is_integer_number() {
                        return None;
                    }
                }
                let (_, d) = n.numer_denom()?;
                let d: i64 = num_traits::ToPrimitive::to_i64(&d)?;
                lcm = num_integer::lcm(lcm, d);
            }
            Some(lcm)
        }
        let beta = compute(&self.bm);
        let alpha = compute(&self.an);
        let two_pi = |k: i64| mul2(&Expr::int(2 * k), &Expr::pi());
        let (p, q) = (self.p(), self.q());
        let pick = if p == q {
            match (alpha, beta) {
                (Some(a), Some(b)) => Some(num_integer::lcm(a, b)),
                _ => None,
            }
        } else if p < q { beta } else { alpha };
        pick.map_or_else(Expr::oo, two_pi)
    }
}

fn shifted(l: &[Expr], by: &Expr) -> Vec<Expr> { l.iter().map(|a| add(&[a.clone(), by.clone()])).collect() }
fn complement(l: &[Expr]) -> Vec<Expr> { l.iter().map(|a| sub(&Expr::one(), a)).collect() }

fn has_x(e: &Expr, x: &Expr) -> bool { crate::traverse::has(e, x) }

/// Fuzzy `is_polar` of the reference (`None` for most expressions).
fn polar_fuzzy(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Sym(s) => if s.facts().get("polar") == Some(true) { Some(true) } else { None },
        Kind::Func(Func::ExpPolar, _) => Some(true),
        Kind::Func(Func::Special(SpecialFn::PolarLift | SpecialFn::PrincipalBranch), _) => Some(true),
        Kind::Mul(..) => {
            let args = crate::simptools::mul_args(e);
            if !args.iter().any(|a| polar_fuzzy(a) == Some(true)) { return Some(false); }
            Some(args.iter().all(|a| polar_fuzzy(a) == Some(true) || crate::assumptions::is_positive(a) == Some(true)))
        }
        Kind::Pow(b, _) => polar_fuzzy(b),
        _ => None,
    }
}

/// `expand_power_base(e)` (deep, without force).
pub(crate) fn expand_power_base(e: &Expr) -> Expr {
    let e = if e.is_atom() { e.clone() } else { crate::simptools::map_args(e, &mut |c| expand_power_base(c)) };
    let Kind::Pow(b, ex) = e.kind() else { return e };
    if !b.is_mul() { return e; }
    let cargs = crate::simptools::mul_args(b);
    let (mut other, mut nonneg, mut neg_) = (Vec::new(), Vec::new(), Vec::new());
    for c in cargs {
        if crate::assumptions::fact(&c, "extended_real") == Some(false) { other.push(c); continue; }
        match polar_fuzzy(&c) {
            Some(true) => nonneg.push(c),
            None => match crate::assumptions::fact(&c, "extended_nonnegative") {
                Some(true) => nonneg.push(c),
                Some(false) => neg_.push(c),
                None => other.push(c),
            },
            Some(false) => other.push(c),
        }
    }
    let cargs = if crate::assumptions::is_integer(ex) == Some(true) {
        let mut all = nonneg;
        all.extend(neg_);
        all.extend(other.drain(..));
        all
    } else {
        if neg_.len() > 1 {
            let mut o = Expr::one();
            if other.is_empty() && neg_[0].as_number().is_some() {
                o = mul2(&o, &neg_.remove(0));
            }
            if neg_.len() % 2 == 1 { o = neg(&o); }
            for n in &neg_ { nonneg.push(neg(n)); }
            if !o.is_one() { other.push(o); }
        } else if !neg_.is_empty() && !other.is_empty() {
            if neg_[0].as_number().is_some() && neg_[0] != Expr::minus_one() {
                other.push(Expr::minus_one());
                nonneg.push(neg(&neg_[0]));
            } else {
                other.extend(neg_);
            }
        } else {
            other.extend(neg_);
        }
        nonneg
    };
    let mut parts: Vec<Expr> = cargs.iter().map(|c| pow(c, ex)).collect();
    if !other.is_empty() { parts.push(pow(&mul(&other), ex)); }
    mul(&parts)
}

/// `_get_coeff_exp`: `expr = c*x**b`; `None` when it has another form.
pub(crate) fn get_coeff_exp(expr: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    let e = expand_power_base(&crate::simptools::powsimp::powsimp(expr, &Default::default()));
    let (c, m) = crate::decompose::as_coeff_mul(&e, std::slice::from_ref(x), false);
    match m.as_slice() {
        [] => Some((c, Expr::zero())),
        [m] => match m.kind() {
            Kind::Pow(b, ex) if b == x => Some((c, ex.clone())),
            _ if m == x => Some((c, Expr::one())),
            _ => None,
        },
        _ => None,
    }
}

/// `_split_mul`: `f = fac * po * g` with `fac` constant and `po = x**s`.
fn split_mul(f: &Expr, x: &Expr) -> (Expr, Expr, Expr) {
    let (mut fac, mut po, mut g) = (Vec::new(), Vec::new(), Vec::new());
    let f = expand_power_base(f);
    for a in crate::simptools::mul_args(&f) {
        if &a == x {
            po.push(x.clone());
        } else if !has_x(&a, x) {
            fac.push(a);
        } else {
            if let Kind::Pow(base, ex) = a.kind() {
                if !has_x(ex, x) {
                    let (mut c, mut t) = crate::decompose::as_coeff_mul(base, std::slice::from_ref(x), false);
                    if t.as_slice() != std::slice::from_ref(x) {
                        (c, t) = crate::decompose::as_coeff_mul(&crate::simp::expand_mul(base), std::slice::from_ref(x), false);
                    }
                    if t.as_slice() == std::slice::from_ref(x) {
                        po.push(pow(x, ex));
                        fac.push(unpolarify(&polarify(&pow(&c, ex), false), false));
                        continue;
                    }
                }
            }
            g.push(a);
        }
    }
    // The reference accumulates with `*=`, i.e. evaluated products.
    let prod = |v: &[Expr]| v.iter().fold(Expr::one(), |acc, t| mul2(&acc, t));
    (prod(&fac), prod(&po), prod(&g))
}

/// `_mul_args`: factors with integer powers repeated.
fn mul_args_rep(f: &Expr) -> Vec<Expr> {
    let mut gs = Vec::new();
    for g in crate::simptools::mul_args(f) {
        if let Kind::Pow(b, e) = g.kind() {
            if let Some(n) = e.as_i64().filter(|_| e.is_integer_number()) {
                let base = if n < 0 { pow(b, &Expr::minus_one()) } else { b.clone() };
                for _ in 0..n.unsigned_abs() { gs.push(base.clone()); }
                continue;
            }
        }
        gs.push(g);
    }
    gs
}

/// `_mul_as_two_parts`, in the order of `multiset_partitions(gs, 2)`.
fn mul_as_two_parts(f: &Expr) -> Option<Vec<(Expr, Expr)>> {
    let gs = mul_args_rep(f);
    if gs.len() < 2 { return None; }
    if gs.len() == 2 { return Some(vec![(gs[0].clone(), gs[1].clone())]); }
    let n = gs.len();
    if gs.iter().all(|g| g == &gs[0]) {
        let mut out = Vec::new();
        for k in 1..=n / 2 {
            out.push((pow(&gs[0], &Expr::int(k as i64)), pow(&gs[0], &Expr::int((n - k) as i64))));
        }
        return Some(out);
    }
    let sorted = crate::simptools::ordered(gs);
    let mut elements: Vec<Expr> = Vec::new();
    let mut mult: Vec<usize> = Vec::new();
    for g in sorted {
        if elements.last() == Some(&g) { *mult.last_mut().unwrap() += 1; } else { elements.push(g); mult.push(1); }
    }
    // Enumerate first-part multiplicity vectors A with A >= mult - A
    // (lexicographically), in decreasing lexicographic order.
    let mut vectors: Vec<Vec<usize>> = Vec::new();
    let mut cur = mult.clone();
    loop {
        let rest: Vec<usize> = mult.iter().zip(&cur).map(|(m, c)| m - c).collect();
        let total: usize = cur.iter().sum();
        if total != 0 && total != n && cur >= rest { vectors.push(cur.clone()); }
        // decrement like a mixed-radix counter
        let mut i = cur.len();
        loop {
            if i == 0 { break; }
            i -= 1;
            if cur[i] > 0 { cur[i] -= 1; for j in i + 1..cur.len() { cur[j] = mult[j]; } i = usize::MAX; break; }
        }
        if i != usize::MAX { break; }
    }
    let part = |v: &[usize]| {
        let mut fs = Vec::new();
        for (e, k) in elements.iter().zip(v) { for _ in 0..*k { fs.push(e.clone()); } }
        mul(&fs)
    };
    Some(vectors.iter().map(|a| {
        let b: Vec<usize> = mult.iter().zip(a).map(|(m, c)| m - c).collect();
        (part(a), part(&b))
    }).collect())
}

/// `_exponents`: exponents of `x` (other than zero) in `e`.
fn exponents(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    if e == x { if !out.contains(&Expr::one()) { out.push(Expr::one()); } return; }
    if let Kind::Pow(b, ex) = e.kind() {
        if b == x { if !out.contains(ex) { out.push(ex.clone()); } return; }
    }
    if e.is_atom() { return; }
    for a in e.args() { exponents(&a, x, out); }
}

fn is_sympy_function(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Logic(_), _) => false,
        Kind::Func(Func::Undefined(n), _) if n.is_internal("Tuple") => false,
        Kind::Func(..) => true,
        _ => false,
    }
}

/// `_functions`: the heads of the functions of `x` in `e`.
fn functions(e: &Expr, x: &Expr) -> usize {
    let mut heads: Vec<Func> = Vec::new();
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(h, _) = s.kind() {
            if is_sympy_function(s) && has_x(s, x) && !heads.contains(h) { heads.push(h.clone()); }
        }
    });
    heads.len()
}

/// `_find_splitting_points`: the roots of innermost linear arguments.
pub(crate) fn find_splitting_points(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    // SymPy: `if not isinstance(expr, Expr): return`; a Symbol is an Expr
    if matches!(e.kind(), Kind::Rel(..) | Kind::Bool(_))
        || (!matches!(e.kind(), Kind::Sym(_)) && crate::logic::is_boolean(e)) { return; }
    if let Kind::Func(Func::Undefined(n), _) = e.kind() { if n.is_internal("Tuple") { return; } }
    if let Some((p, q)) = linear(e, x) {
        if !p.is_zero() {
            let r = neg(&div(&q, &p));
            if !out.contains(&r) { out.push(r); }
            return;
        }
    }
    if e.is_atom() { return; }
    for a in e.args() { find_splitting_points(&a, x, out); }
}

/// `e.match(p*x + q)` with `p`, `q` free of `x`.
fn linear(e: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    let (q, dep) = crate::decompose::as_independent(e, std::slice::from_ref(x), Some(true));
    if dep.is_zero() { return Some((Expr::zero(), q)); }
    let mut ps = Vec::new();
    for t in crate::simptools::add_args(&dep) {
        let (c, rest) = crate::decompose::as_independent(&t, std::slice::from_ref(x), Some(false));
        if &rest != x { return None; }
        ps.push(c);
    }
    Some((add(&ps), q))
}

/// `_inflate_g`: `g = C*h` with `h` a G-function of argument `z**n`.
fn inflate_g(g: &G, n: i64) -> Option<(Expr, G)> {
    let inflate = |l: &[Expr]| {
        let mut out = Vec::new();
        for a in l { for i in 0..n { out.push(div(&add(&[a.clone(), Expr::int(i)]), &Expr::int(n))); } }
        out
    };
    let v = Expr::int(g.p() as i64 - g.q() as i64);
    let nn = Expr::int(n);
    let c = div(
        &pow(&nn, &add(&[Expr::one(), g.nu(), div(&v, &Expr::int(2))])),
        &pow(&mul2(&Expr::int(2), &Expr::pi()), &mul2(&Expr::int(n - 1), &g.delta())),
    );
    let z = mul2(&pow(&g.z, &nn), &pow(&nn, &mul2(&nn, &v)));
    Some((c, G::new(inflate(&g.an), inflate(&g.ao), inflate(&g.bm), inflate(&g.bo), z)?))
}

/// `_flip_g`: a G-function of the inverse argument.
fn flip_g(g: &G) -> Option<G> {
    G::new(complement(&g.bm), complement(&g.bo), complement(&g.an), complement(&g.ao), pow(&g.z, &Expr::minus_one()))
}

fn principal_branch(e: &Expr, period: &Expr) -> Expr {
    crate::special::apply(SpecialFn::PrincipalBranch, vec![e.clone(), period.clone()])
}

/// `_my_principal_branch`.
fn my_principal_branch(e: &Expr, period: &Expr, full_pb: bool) -> Expr {
    let res = principal_branch(e, period);
    if full_pb { return res; }
    fn strip(e: &Expr) -> Expr {
        if let Kind::Func(Func::Special(SpecialFn::PrincipalBranch), a) = e.kind() { return strip(&a[0]); }
        if e.is_atom() { return e.clone(); }
        crate::simptools::map_args(e, &mut |c| strip(c))
    }
    strip(&res)
}

/// `_rewrite_saxena_1`: `fac*po*g = C*G(a*x)`.
fn rewrite_saxena_1(fac: &Expr, po: &Expr, g: &G, x: &Expr) -> Option<(Expr, G)> {
    let (_, s) = get_coeff_exp(po, x)?;
    let (a, b) = get_coeff_exp(&g.z, x)?;
    let a = my_principal_branch(&a, &g.period(), false);
    let shift = sub(&div(&add(&[s, Expr::one()]), &b), &Expr::one());
    let c = div(fac, &mul2(&f::abs(&b), &pow(&a, &shift)));
    Some((c, G::new(shifted(&g.an, &shift), shifted(&g.ao, &shift), shifted(&g.bm, &shift), shifted(&g.bo, &shift), mul2(&a, x))?))
}

/// `_int0oo_1`: the integral of `G(eta*x)` over (0, oo).
fn int0oo_1(g: &G, x: &Expr) -> Option<Expr> {
    let (eta, _) = get_coeff_exp(&g.z, x)?;
    let mut res = div(&Expr::one(), &eta);
    for b in &g.bm { res = mul2(&res, &f::gamma(&add(&[b.clone(), Expr::one()]))); }
    for a in &g.an { res = mul2(&res, &f::gamma(&sub(&neg(a), &Expr::zero()))); }
    for b in &g.bo { res = div(&res, &f::gamma(&neg(b))); }
    for a in &g.ao { res = div(&res, &f::gamma(&add(&[a.clone(), Expr::one()]))); }
    Some(crate::simptools::gammasimp::gammasimp(&unpolarify(&res, false)))
}

/// `powdenest(e, polar=True)`: denest treating every symbol as polar.
pub(crate) fn powdenest_polar(e: &Expr) -> Expr {
    let p = polarify(e, false);
    let mut fwd = Vec::new();
    let mut back = Vec::new();
    for s in crate::traverse::free_symbols(&p) {
        let name = match s.kind() { Kind::Sym(sym) => sym.name().to_string(), _ => continue };
        let facts = match crate::assumption_facts::FactSet::from_kwargs([("polar", true)]) { Ok(f) => f, Err(_) => return e.clone() };
        let Ok(d) = crate::symbol::Symbol::dummy_with_facts(&name, facts, None) else { return e.clone() };
        let d = Expr::from_kind(Kind::Sym(d));
        fwd.push((s.clone(), d.clone()));
        back.push((d, s));
    }
    let p = crate::subs::subs(&p, &fwd);
    let inner = crate::simptools::powsimp::powdenest(&unpolarify(&p, true), false, false);
    unpolarify(&crate::subs::subs(&inner, &back), false)
}

/// `_rewrite_saxena`: bring `fac*po*g1*g2` to `C*G1(a1*x)*G2(a2*x**b)`
/// with `G1` absorbing the power.
fn rewrite_saxena(fac: &Expr, po: &Expr, g1: &G, g2: &G, x: &Expr, full_pb: bool) -> Option<Option<(Expr, G, G)>> {
    let pb = |g: &G| -> Option<G> {
        let (a, b) = get_coeff_exp(&g.z, x)?;
        g.with_z(mul2(&my_principal_branch(&a, &g.period(), full_pb), &pow(x, &b)))
    };
    let (_, s) = get_coeff_exp(po, x)?;
    let (_, mut b1) = get_coeff_exp(&g1.z, x)?;
    let (_, mut b2) = get_coeff_exp(&g2.z, x)?;
    let (mut g1, mut g2) = (g1.clone(), g2.clone());
    let negative = |b: &Expr| b.as_number().is_some_and(|n| n.is_negative())
        || crate::assumptions::is_negative(b) == Some(true);
    if negative(&b1) { b1 = neg(&b1); g1 = flip_g(&g1)?; }
    if negative(&b2) { b2 = neg(&b2); g2 = flip_g(&g2)?; }
    let rational = |b: &Expr| b.as_number().filter(|n| n.is_rational()).and_then(|n| n.numer_denom());
    let (Some((m1, n1)), Some((m2, n2))) = (rational(&b1), rational(&b2)) else { return Some(None) };
    use num_traits::ToPrimitive;
    let (m1, n1, m2, n2) = (m1.to_i64()?, n1.to_i64()?, m2.to_i64()?, n2.to_i64()?);
    if m1 == 0 || m2 == 0 { return None; }
    let tau = num_integer::lcm(m1 * n2, m2 * n1);
    let r1 = tau / (m1 * n2);
    let r2 = tau / (m2 * n1);
    let (c1, g1) = inflate_g(&g1, r1)?;
    let (c2, g2) = inflate_g(&g2, r2)?;
    let g1 = pb(&g1)?;
    let g2 = pb(&g2)?;
    let fac = mul(&[fac.clone(), c1, c2]);
    let (a1, b) = get_coeff_exp(&g1.z, x)?;
    let (a2, _) = get_coeff_exp(&g2.z, x)?;
    let ex = sub(&div(&add(&[s, Expr::one()]), &b), &Expr::one());
    let fac = div(&fac, &mul2(&f::abs(&b), &pow(&a1, &ex)));
    let g1 = G::new(shifted(&g1.an, &ex), shifted(&g1.ao, &ex), shifted(&g1.bm, &ex), shifted(&g1.bo, &ex), mul2(&a1, x))?;
    let g2 = g2.with_z(mul2(&a2, x))?;
    Some(Some((powdenest_polar(&fac), g1, g2)))
}

/// `_int0oo`: the integral of `G1(eta*x)*G2(omega*x)` over (0, oo).
fn int0oo(g1: &G, g2: &G, x: &Expr) -> Option<Expr> {
    let (eta, _) = get_coeff_exp(&g1.z, x)?;
    let (omega, _) = get_coeff_exp(&g2.z, x)?;
    let negs = |l: &[Expr]| l.iter().map(neg).collect::<Vec<_>>();
    let mut a1 = negs(&g1.bm); a1.extend(g2.an.iter().cloned());
    let mut a2 = g2.ao.clone(); a2.extend(negs(&g1.bo));
    let mut b1 = negs(&g1.an); b1.extend(g2.bm.iter().cloned());
    let mut b2 = g2.bo.clone(); b2.extend(negs(&g1.ao));
    let g = G::new(a1, a2, b1, b2, div(&omega, &eta))?;
    Some(div(&g.to_expr(), &eta))
}

type Terms = Vec<(Expr, Expr, G)>;

/// `_rewrite_single(f, x, recursive)`.
fn rewrite_single(f: &Expr, x: &Expr, recursive: bool) -> Option<(Terms, Expr)> {
    if let Some(r) = table::rewrite_single(f, x) {
        return Some(r);
    }
    if !recursive {
        return None;
    }
    recursive::rewrite(f, x)
}

/// `_rewrite1`: `f = fac*po*sum(C*x**s*G)`.
fn rewrite1(f: &Expr, x: &Expr, recursive: bool) -> Option<(Expr, Expr, Terms, Expr)> {
    let (fac, po, g) = split_mul(f, x);
    let (terms, cond) = rewrite_single(&g, x, recursive)?;
    Some((fac, po, terms, cond))
}

/// `_rewrite2`: `f = fac*po*(sum G)*(sum G)`.
fn rewrite2(f: &Expr, x: &Expr) -> Option<(Expr, Expr, Terms, Terms, Expr)> {
    let (fac, po, g) = split_mul(f, x);
    if mul_args_rep(&g).iter().any(|e| table::rewrite_single(e, x).is_none()) { return None; }
    let l = mul_as_two_parts(&g)?;
    let key = |p: &(Expr, Expr)| {
        let ex = |e: &Expr| { let mut v = Vec::new(); exponents(e, x, &mut v); v.len() };
        let sp = |e: &Expr| { let mut v = Vec::new(); find_splitting_points(e, x, &mut v); v.len() };
        (
            ex(&p.0).max(ex(&p.1)),
            functions(&p.0, x).max(functions(&p.1, x)),
            sp(&p.0).max(sp(&p.1)),
            1 + crate::simptools::nodes(&p.0) + crate::simptools::nodes(&p.1),
        )
    };
    let mut l: Vec<_> = l.into_iter().map(|p| (key(&p), p)).collect();
    l.sort_by(|(ka, a), (kb, b)| ka.cmp(kb)
        .then_with(|| crate::sortkey::cmp_default(&a.0, &b.0))
        .then_with(|| crate::sortkey::cmp_default(&a.1, &b.1)));
    for recursive in [false, true] {
        for (_, (f1, f2)) in &l {
            let Some(r1) = rewrite_single(f1, x, recursive) else { continue };
            let Some(r2) = rewrite_single(f2, x, recursive) else { continue };
            let cond = crate::logic::and(vec![r1.1, r2.1]);
            if !is_false(&cond) {
                return Some((fac, po, r1.0, r2.0, cond));
            }
        }
    }
    None
}

pub(crate) fn is_false(e: &Expr) -> bool { matches!(e.kind(), Kind::Bool(false)) }

fn hyperexpand(e: &Expr) -> Expr { crate::special::hyperexpand(e) }

/// `_meijerint_definite_4`: the integral of `f` over (0, oo).
fn definite_4(f: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    definite_4_opt(f, x, false)
}

/// `_meijerint_definite_4(f, x, only_double)`: `only_double` skips the
/// single G-function rewrite (the recursive Mellin rewrite uses it).
fn definite_4_opt(f: &Expr, x: &Expr, only_double: bool) -> Option<(Expr, Expr)> {
    if let Some((fac, po, g, mut cond)) = if only_double { None } else { rewrite1(f, x, false) } {
        let mut res = Expr::zero();
        let mut ok = true;
        for (c, s, gf) in &g {
            if c.is_zero() { continue; }
            let Some((c, gf)) = rewrite_saxena_1(&mul2(&fac, c), &mul2(&po, &pow(x, s)), gf, x) else { ok = false; break };
            let Some(v) = int0oo_1(&gf, x) else { ok = false; break };
            res = add(&[res, mul2(&c, &v)]);
            cond = crate::logic::and(vec![cond, conds::check_antecedents_1(&gf, x, false)?]);
            if is_false(&cond) { break; }
        }
        if ok {
            let cond = conds::my_unpolarify(&cond);
            if !is_false(&cond) {
                return Some((conds::my_unpolarify(&hyperexpand(&res)), cond));
            }
        }
    }
    let (fac, po, g1, g2, cond0) = rewrite2(f, x)?;
    for full_pb in [false, true] {
        let mut cond = cond0.clone();
        let mut res = Expr::zero();
        'outer: for (c1, s1, f1) in &g1 {
            for (c2, s2, f2) in &g2 {
                let r = rewrite_saxena(&mul(&[fac.clone(), c1.clone(), c2.clone()]), &mul2(&po, &pow(x, &add(&[s1.clone(), s2.clone()]))), f1, f2, x, full_pb)?;
                let Some((c, f1_, f2_)) = r else { return None };
                cond = crate::logic::and(vec![cond, conds::check_antecedents(&f1_, &f2_, x)?]);
                if is_false(&cond) { break 'outer; }
                res = add(&[res, mul2(&c, &int0oo(&f1_, &f2_, x)?)]);
            }
        }
        let cond = conds::my_unpolarify(&cond);
        if !is_false(&cond) {
            return Some((conds::my_unpolarify(&hyperexpand(&res)), cond));
        }
    }
    None
}

/// `_meijerint_definite_3`: also try the terms of a sum separately.
fn definite_3(f: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    if let Some(r) = definite_4(f, x) {
        if !is_false(&r.1) { return Some(r); }
    }
    if f.is_add() {
        let mut vals = Vec::new();
        let mut cs = Vec::new();
        for g in crate::simptools::add_args(f) {
            let (v, c) = definite_4(&g, x)?;
            vals.push(v);
            cs.push(c);
        }
        let c = crate::logic::and(cs);
        if !is_false(&c) { return Some((add(&vals), c)); }
    }
    None
}

fn has_head_in(e: &Expr, heads: &[Func]) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(h, _) = s.kind() { if heads.contains(h) { found = true; } }
    });
    found
}

const TRIG: &[Func] = &[Func::Sin, Func::Cos, Func::Tan, Func::Cot, Func::Sec, Func::Csc];
fn hyperbolic_heads() -> Vec<Func> {
    vec![Func::Sinh, Func::Cosh, Func::Tanh, Func::Coth, Func::Special(SpecialFn::Sech), Func::Special(SpecialFn::Csch)]
}

/// `_guess_expansion`.
fn guess_expansion(f: &Expr) -> Vec<Expr> {
    let mut res = vec![f.clone()];
    let push = |e: Expr, res: &mut Vec<Expr>| if !res.contains(&e) { res.push(e) };
    push(crate::simp::expand_mul(f), &mut res);
    push(crate::expand::expand(f), &mut res);
    let mut trig: Vec<Func> = TRIG.to_vec();
    trig.extend(hyperbolic_heads());
    if has_head_in(f, &trig) {
        push(crate::simp::expand_mul(&crate::simp::expand_trig(f)), &mut res);
    }
    if has_head_in(f, &[Func::Sin, Func::Cos]) {
        push(crate::simptools::fu::sincos_to_sum(f), &mut res);
    }
    res
}

fn positive_dummy(name: &str) -> Expr {
    Expr::dummy(name, crate::symbol::Assumptions::POSITIVE)
}

/// `_meijerint_definite_2`: the integral over (0, oo).
fn definite_2(f: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    let d = positive_dummy("x");
    let f = crate::subs::subs1(f, x, &d);
    if f.is_zero() { return Some((Expr::zero(), Expr::boolean(true))); }
    for g in guess_expansion(&f) {
        if let Some(r) = definite_3(&g, &d) { return Some(r); }
    }
    None
}

/// `_has(res, meijerg)`: after `piecewise_fold`, every piece contains a
/// Meijer G-function. A node outside a Piecewise folds into the products of
/// the pieces of its arguments, so all of them contain one exactly when one
/// argument has a G-function in every piece.
fn has_meijerg(e: &Expr) -> bool {
    match e.kind() {
        Kind::Piecewise(pieces) => pieces.iter().all(|(v, c)| has_meijerg(v) || crate::special::has_head(c, SpecialFn::MeijerG)),
        Kind::Func(Func::Special(SpecialFn::MeijerG), _) => true,
        _ => e.args().iter().any(has_meijerg),
    }
}

fn is_inf(e: &Expr) -> bool { e.is_constant(Constant::Infinity) }
fn is_neg_inf(e: &Expr) -> bool { e.is_constant(Constant::NegInfinity) }

/// `meijerint_definite(f, x, a, b)`: `(value, condition)` or `None` when
/// the method does not apply.
pub(crate) fn meijerint_definite(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Option<(Expr, Expr)> {
    let mut singular = false;
    crate::traverse::preorder(f, &mut |s| {
        if matches!(s.kind(), Kind::Func(Func::DiracDelta, _)) {
            singular = true;
        }
        if let Kind::Func(Func::Undefined(u), _) = s.kind() {
            if &**u == "SingularityFunction" {
                singular = true;
            }
        }
    });
    if singular {
        return None;
    }
    let (f0, x0, a0, b0) = (f.clone(), x.clone(), a.clone(), b.clone());
    let d = Expr::dummy("x", crate::symbol::Assumptions::default());
    let mut f = crate::subs::subs1(f, x, &d);
    let x = d;
    if a == b { return Some((Expr::zero(), Expr::boolean(true))); }
    let mut results: Vec<(Expr, Expr)> = Vec::new();
    if is_neg_inf(a) && !is_inf(b) {
        return meijerint_definite(&crate::subs::subs1(&f, &x, &neg(&x)), &x, &neg(b), &neg(a));
    } else if is_neg_inf(a) {
        let mut points = Vec::new();
        find_splitting_points(&f, &x, &mut points);
        points.sort_by(|p, q| crate::sortkey::cmp_default(q, p));
        points.push(Expr::zero());
        for c in points {
            if crate::assumptions::fact(&c, "extended_real") != Some(true) { continue; }
            let Some((r1, c1)) = definite_2(&crate::subs::subs1(&f, &x, &add(&[x.clone(), c.clone()])), &x) else { continue };
            let Some((r2, c2)) = definite_2(&crate::subs::subs1(&f, &x, &sub(&c, &x)), &x) else { continue };
            let cond = conds::condsimp(&crate::logic::and(vec![c1, c2]), true);
            if is_false(&cond) { continue; }
            return Some((add(&[r1, r2]), cond));
        }
    } else if is_inf(a) {
        let (v, c) = meijerint_definite(&f, &x, b, &Expr::oo())?;
        return Some((neg(&v), c));
    } else if a.is_zero() && is_inf(b) {
        if let Some(res) = definite_2(&f, &x) {
            if has_meijerg(&res.0) { results.push(res); } else { return Some(res); }
        }
    } else {
        let mut b = b.clone();
        if is_inf(&b) {
            let mut points = Vec::new();
            find_splitting_points(&f, &x, &mut points);
            for split in points {
                if crate::assumptions::is_nonnegative(&sub(a, &split)) == Some(true) {
                    let g = mul2(&crate::subs::subs1(&f, &x, &add(&[x.clone(), split.clone()])),
                        &f::heaviside(&sub(&add(&[x.clone(), split.clone()]), a)));
                    if let Some(res) = definite_2(&g, &x) {
                        if has_meijerg(&res.0) { results.push(res); } else { return Some(res); }
                    }
                }
            }
        }
        f = crate::subs::subs1(&f, &x, &add(&[x.clone(), a.clone()]));
        b = sub(&b, a);
        if !is_inf(&b) {
            let phi = f::exp(&mul2(&Expr::i(), &f::arg(&b)));
            b = f::abs(&b);
            f = crate::subs::subs1(&f, &x, &mul2(&phi, &x));
            f = mul(&[f, f::heaviside(&sub(&b, &x)), phi]);
        }
        if let Some(res) = definite_2(&f, &x) {
            if has_meijerg(&res.0) { results.push(res); } else { return Some(res); }
        }
    }
    let mut hyp = hyperbolic_heads();
    hyp.truncate(6);
    if has_head_in(&f0, &hyp) {
        let rewritten = crate::special::rewrite::rewrite(
            &f0,
            &["sinh", "cosh", "tanh", "coth", "sech", "csch"].map(String::from),
            "exp",
            true,
        );
        if let Some((v, c)) = meijerint_definite(&rewritten, &x0, &a0, &b0) {
            let mut exps = Vec::new();
            crate::traverse::preorder(&v, &mut |s| if s.is_func(&Func::Exp) && !exps.contains(s) { exps.push(s.clone()) });
            let v = crate::poly::api::collect(&crate::simptools::exprtools::factor_terms(&v, false, false, false, true), &exps);
            return Some((v, c));
        }
    }
    results.into_iter().min_by(|(r1, c1), (r2, c2)| {
        (crate::simptools::nodes(r1) + crate::simptools::nodes(c1))
            .cmp(&(crate::simptools::nodes(r2) + crate::simptools::nodes(c2)))
            .then_with(|| crate::sortkey::cmp_default(r1, r2))
            .then_with(|| crate::sortkey::cmp_default(c1, c2))
    })
}

/// Internal probe of the Meijer G pipeline, used to compare each stage with
/// the reference implementation.
pub(crate) fn debug(op: &str, args: &[Expr]) -> Option<Vec<Expr>> {
    if op.starts_with("hx") {
        return crate::special::hyperexpand::debug(op, args);
    }
    let tuple = |v: Vec<Expr>| crate::special::tuple(v);
    match (op, args) {
        ("rewrite_single", [f, x]) => {
            let (terms, cond) = table::rewrite_single(f, x)?;
            let t = terms.into_iter().map(|(fac, po, g)| tuple(vec![fac, po, g.to_expr()])).collect();
            Some(vec![tuple(t), cond])
        }
        ("rewrite_single_recursive", [f, x]) => {
            let (terms, cond) = rewrite_single(f, x, true)?;
            let t = terms.into_iter().map(|(fac, po, g)| tuple(vec![fac, po, g.to_expr()])).collect();
            Some(vec![tuple(t), cond])
        }
        ("rewrite1", [f, x]) => {
            let (fac, po, gs, cond) = rewrite1(f, x, true)?;
            let t = gs.into_iter().map(|(c, s, g)| tuple(vec![c, s, g.to_expr()])).collect();
            Some(vec![fac, po, tuple(t), cond])
        }
        ("rewrite2", [f, x]) => {
            let (fac, po, g1, g2, cond) = rewrite2(f, x)?;
            let t1 = g1.into_iter().map(|(c, s, g)| tuple(vec![c, s, g.to_expr()])).collect();
            let t2 = g2.into_iter().map(|(c, s, g)| tuple(vec![c, s, g.to_expr()])).collect();
            Some(vec![fac, po, tuple(t1), tuple(t2), cond])
        }
        ("check_antecedents_1", [g, x]) => Some(vec![conds::check_antecedents_1(&G::from_expr(g)?, x, false)?]),
        ("check_antecedents_1", [g, x, h]) => {
            let helper = matches!(h.kind(), Kind::Bool(true));
            Some(vec![conds::check_antecedents_1(&G::from_expr(g)?, x, helper)?])
        }
        ("check_antecedents", [g1, g2, x]) => {
            Some(vec![conds::check_antecedents(&G::from_expr(g1)?, &G::from_expr(g2)?, x)?])
        }
        ("inversion", [f, x, t]) => Some(vec![meijerint_inversion(f, x, t)?]),
        ("indefinite", [f, x]) => Some(vec![meijerint_indefinite(f, x)?]),
        ("piecewise_fold", [e]) => Some(vec![crate::integrate::piecewise::piecewise_fold(e, Some(true))]),
        ("piecewise", pieces) => Some(vec![crate::relational::piecewise(
            pieces.chunks(2).map(|p| (p[0].clone(), p[1].clone())).collect(),
        )]),
        ("powdenest_polar", [e]) => Some(vec![powdenest_polar(e)]),
        ("polarify", [e]) => Some(vec![polarify(e, false)]),
        ("condsimp", [c]) => Some(vec![conds::condsimp(c, true)]),
        ("eval_cond", [c]) => Some(vec![conds::eval_cond(c)]),
        ("int0oo_1", [g, x]) => Some(vec![int0oo_1(&G::from_expr(g)?, x)?]),
        ("int0oo", [g1, g2, x]) => Some(vec![int0oo(&G::from_expr(g1)?, &G::from_expr(g2)?, x)?]),
        ("definite_4", [f, x]) => {
            let (v, c) = definite_4(f, x)?;
            Some(vec![v, c])
        }
        ("split_points", [f, x]) => {
            let mut v = Vec::new();
            find_splitting_points(f, x, &mut v);
            Some(v)
        }
        ("split_try", [f, x, a, split]) => {
            let g = mul2(&crate::subs::subs1(f, x, &add(&[x.clone(), split.clone()])),
                &f::heaviside(&sub(&add(&[x.clone(), split.clone()]), a)));
            let nn = crate::assumptions::is_nonnegative(&sub(a, split));
            let (v, c) = definite_2(&g, x)?;
            Some(vec![g, Expr::boolean(nn == Some(true)), v, c])
        }
        ("definite", [f, x, a, b]) => {
            let (v, c) = meijerint_definite(f, x, a, b)?;
            Some(vec![v, c])
        }
        _ => None,
    }
}
