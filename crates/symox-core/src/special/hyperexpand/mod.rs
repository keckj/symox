//! A port of SymPy's `sympy.simplify.hyperexpand`: closed forms of
//! generalised hypergeometric functions and Meijer G functions.
//!
//! A function is reduced in order, matched against a table of known origins
//! (or a Lerch phi or hypergeometric basis is built), and the target is
//! reached from the origin by a plan of differential shift operators acting
//! on the coefficient vector of the origin's basis. Meijer G functions use a
//! direct formula when one matches and Slater's theorem otherwise.

mod formulae;
mod hyperrep;
mod ops;
mod plan;
mod table;
mod zeta;
pub(crate) use zeta::{lerchphi_expand, polylog_expand};

use std::cell::RefCell;

use crate::ops::{add, div, mul, sub};
use crate::special::polar::{polarify, unpolarify};
use crate::special::{tuple, tuple_items, SpecialFn};
use crate::{Expr, Func, Kind};

pub(crate) use hyperrep::{extract_branch_factor, hyperrep_eval, rewrite as hyperrep_rewrite};

use ops::{apply_operators, make_simp, Op};

/// The exceptions the reference raises out of `hyperexpand`.
#[derive(Clone, Debug)]
pub(crate) enum HxError {
    Value(String),
    NotImplemented(String),
    Type(String),
    Polynomial(String),
}

pub(crate) type HResult<T> = Result<T, HxError>;

impl std::fmt::Display for HxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HxError::Value(s) => write!(f, "ValueError: {s}"),
            HxError::NotImplemented(s) => write!(f, "NotImplementedError: {s}"),
            HxError::Type(s) => write!(f, "TypeError: {s}"),
            HxError::Polynomial(s) => write!(f, "PolynomialError: {s}"),
        }
    }
}

/// The `rewrite` keyword of `hyperexpand`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Rewrite {
    Default,
    None,
    Nonrep,
    Nonrepsmall,
}

pub(crate) fn half() -> Expr { Expr::rational(1, 2) }
pub(crate) fn int(i: i64) -> Expr { Expr::int(i) }
pub(crate) fn ipi() -> Expr { mul(&[Expr::i(), Expr::pi()]) }

/// `e.is_Integer`: an exact integer number.
pub(crate) fn is_integer_number(e: &Expr) -> bool { crate::special::is_int_number(e) }

/// `e.is_integer and e.is_nonpositive`.
pub(crate) fn is_nonpositive_integer(e: &Expr) -> bool {
    crate::special::is_integer(e) == Some(true) && crate::special::is_nonpositive(e) == Some(true)
}

/// Sign of the number `e` (None when `e` is not a comparable number).
pub(crate) fn num_cmp(e: &Expr) -> Option<std::cmp::Ordering> {
    crate::simptools::num_sign(e)
}

/// `bool(a < b)` for numeric differences, as Python's sort needs it.
fn lt_num(a: &Expr, b: &Expr) -> HResult<bool> {
    match num_cmp(&sub(a, b)) {
        Some(o) => Ok(o == std::cmp::Ordering::Less),
        None => Err(HxError::Type("cannot determine truth value of Relational".into())),
    }
}

/// Python's stable `list.sort()` on SymPy numbers.
pub(crate) fn sort_num(v: &mut [Expr]) -> HResult<()> {
    let mut err = None;
    v.sort_by(|a, b| match (lt_num(a, b), lt_num(b, a)) {
        (Ok(true), _) => std::cmp::Ordering::Less,
        (_, Ok(true)) => std::cmp::Ordering::Greater,
        (Ok(false), Ok(false)) => std::cmp::Ordering::Equal,
        (Err(e), _) | (_, Err(e)) => {
            err.get_or_insert(e);
            std::cmp::Ordering::Equal
        }
    });
    match err {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

/// `sorted(v, key=default_sort_key)`.
pub(crate) fn sort_default(v: &mut [Expr]) { v.sort_by(crate::sortkey::cmp_default); }

/// `_mod1(x)`: the class of `x` modulo the integers.
pub(crate) fn mod1(x: &Expr) -> Expr {
    if x.as_number().is_some() {
        return sub(x, &crate::functions::floor(x));
    }
    let (c, rest) = crate::ops::as_coeff_add(x);
    let c = Expr::num(c);
    add(&[sub(&c, &crate::functions::floor(&c)), rest])
}

/// `sift(params, key)`: buckets in order of first appearance.
pub(crate) type Buckets = Vec<(Expr, Vec<Expr>)>;

pub(crate) fn sift(params: &[Expr], key: &mut dyn FnMut(&Expr) -> Expr) -> Buckets {
    let mut out: Buckets = Vec::new();
    for p in params {
        let k = key(p);
        match out.iter_mut().find(|(kk, _)| *kk == k) {
            Some((_, v)) => v.push(p.clone()),
            None => out.push((k, vec![p.clone()])),
        }
    }
    out
}

pub(crate) fn sift_mod1(params: &[Expr]) -> Buckets { sift(params, &mut mod1) }

pub(crate) fn bucket<'a>(b: &'a Buckets, k: &Expr) -> Option<&'a Vec<Expr>> {
    b.iter().find(|(kk, _)| kk == k).map(|(_, v)| v)
}

/// `Hyper_Function(ap, bq)`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct HF {
    pub ap: Vec<Expr>,
    pub bq: Vec<Expr>,
}

pub(crate) type Invariants = (usize, Vec<(Expr, usize)>, Vec<(Expr, usize)>);

impl HF {
    pub fn new(ap: Vec<Expr>, bq: Vec<Expr>) -> HF {
        HF { ap: ap.iter().map(crate::expand::expand).collect(), bq: bq.iter().map(crate::expand::expand).collect() }
    }

    pub fn sizes(&self) -> (usize, usize) { (self.ap.len(), self.bq.len()) }

    pub fn gamma(&self) -> usize {
        self.ap
            .iter()
            .filter(|x| crate::special::is_integer(x) == Some(true) && crate::special::is_negative(x) == Some(true))
            .count()
    }

    pub fn xreplace(&self, repl: &[(Expr, Expr)]) -> HF {
        let map: rustc_hash::FxHashMap<Expr, Expr> = repl.iter().cloned().collect();
        HF::new(
            self.ap.iter().map(|e| crate::subs::xreplace(e, &map)).collect(),
            self.bq.iter().map(|e| crate::subs::xreplace(e, &map)).collect(),
        )
    }

    pub fn build_invariants(&self) -> Invariants {
        let tr = |b: Buckets| {
            let mut v: Vec<(Expr, usize)> = b.into_iter().filter(|(_, v)| !v.is_empty()).map(|(k, v)| (k, v.len())).collect();
            if !v.iter().any(|(k, _)| is_mod_node(k)) {
                v.sort_by(|a, b| crate::sortkey::cmp_default(&a.0, &b.0));
            }
            v
        };
        (self.gamma(), tr(sift_mod1(&self.ap)), tr(sift_mod1(&self.bq)))
    }

    /// `difficulty(func)`: None for -1.
    pub fn difficulty(&self, func: &HF) -> HResult<Option<Expr>> {
        if self.gamma() != func.gamma() {
            return Ok(None);
        }
        let (oab, obb, ab, bb) = (sift_mod1(&self.ap), sift_mod1(&self.bq), sift_mod1(&func.ap), sift_mod1(&func.bq));
        let mut diff = Expr::zero();
        for (b, ob) in [(&ab, &oab), (&bb, &obb)] {
            let mut keys: Vec<Expr> = b.iter().map(|(k, _)| k.clone()).collect();
            for (k, _) in ob.iter() {
                if !keys.contains(k) {
                    keys.push(k.clone());
                }
            }
            for m in &keys {
                let (Some(l1), Some(l2)) = (bucket(b, m), bucket(ob, m)) else { return Ok(None) };
                if l1.len() != l2.len() {
                    return Ok(None);
                }
                let (mut l1, mut l2) = (l1.clone(), l2.clone());
                sort_num(&mut l1)?;
                sort_num(&mut l2)?;
                for (i, j) in l1.iter().zip(&l2) {
                    diff = add(&[diff, crate::functions::abs(&sub(i, j))]);
                }
            }
        }
        Ok(Some(diff))
    }

    pub fn is_suitable_origin(&self) -> bool {
        for a in &self.ap {
            for b in &self.bq {
                let d = sub(a, b);
                if crate::special::is_integer(&d) == Some(true) && crate::special::is_negative(&d) == Some(false) {
                    return false;
                }
            }
        }
        if self.ap.iter().any(Expr::is_zero) {
            return false;
        }
        !self.bq.iter().any(is_nonpositive_integer)
    }
}

fn is_mod_node(e: &Expr) -> bool { matches!(e.kind(), Kind::Func(Func::Special(crate::special::SpecialFn::Mod), _)) }

/// `G_Function(an, ap, bm, bq)`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct GF {
    pub an: Vec<Expr>,
    pub ap: Vec<Expr>,
    pub bm: Vec<Expr>,
    pub bq: Vec<Expr>,
}

impl GF {
    pub fn new(an: Vec<Expr>, ap: Vec<Expr>, bm: Vec<Expr>, bq: Vec<Expr>) -> GF {
        let ex = |v: Vec<Expr>| v.iter().map(crate::expand::expand).collect();
        GF { an: ex(an), ap: ex(ap), bm: ex(bm), bq: ex(bq) }
    }

    pub fn signature(&self) -> (usize, usize, usize, usize) { (self.an.len(), self.ap.len(), self.bm.len(), self.bq.len()) }

    /// `meijerg(an, ap, bm, bq, z)`: each list ordered and unpolarified.
    pub fn call(&self, z: &Expr) -> Expr {
        let prep = |v: &[Expr]| {
            let mut v = v.to_vec();
            v.sort_by(crate::sortkey::cmp_default);
            tuple(v.iter().map(|x| crate::special::polar::unpolarify(x, false)).collect())
        };
        crate::special::raw(SpecialFn::MeijerG, vec![
            tuple(vec![prep(&self.an), prep(&self.ap)]),
            tuple(vec![prep(&self.bm), prep(&self.bq)]),
            z.clone(),
        ])
    }

    /// `compute_buckets()`: each bucket sorted by `x - x0`, descending for
    /// `an` and `bq`.
    pub fn compute_buckets(&self) -> HResult<[Buckets; 4]> {
        let mut out: [Buckets; 4] = Default::default();
        for (k, (lis, flip)) in [(&self.an, true), (&self.ap, false), (&self.bm, false), (&self.bq, true)].into_iter().enumerate() {
            let mut d = sift_mod1(lis);
            for (_, items) in d.iter_mut() {
                let x0 = items[0].clone();
                let mut keyed: Vec<(Expr, Expr)> = items.iter().map(|x| (sub(x, &x0), x.clone())).collect();
                let mut err = None;
                keyed.sort_by(|a, b| {
                    let (a, b) = if flip { (&b.0, &a.0) } else { (&a.0, &b.0) };
                    match num_cmp(&sub(a, b)) {
                        Some(o) => o,
                        None => {
                            err.get_or_insert(HxError::Type("cannot determine truth value of Relational".into()));
                            std::cmp::Ordering::Equal
                        }
                    }
                });
                if let Some(e) = err {
                    return Err(e);
                }
                *items = keyed.into_iter().map(|(_, x)| x).collect();
            }
            out[k] = d;
        }
        Ok(out)
    }
}

/// A formula: `closed_form = C*B`, `z d/dz B = M*B`.
#[derive(Clone, Debug)]
pub(crate) struct Formula {
    pub func: HF,
    pub z: Expr,
    pub symbols: Vec<Expr>,
    pub b: Vec<Expr>,
    pub c: Vec<Expr>,
    pub m: Vec<Vec<Expr>>,
}

pub(crate) fn subs_all(v: &[Expr], pairs: &[(Expr, Expr)]) -> Vec<Expr> { v.iter().map(|e| crate::subs::subs(e, pairs)).collect() }

pub(crate) fn subs_mat(m: &[Vec<Expr>], pairs: &[(Expr, Expr)]) -> Vec<Vec<Expr>> { m.iter().map(|r| subs_all(r, pairs)).collect() }

/// `e.has(nan, oo, -oo, zoo)`.
pub(crate) fn has_nonfinite(e: &Expr) -> bool {
    use crate::traverse::has;
    has(e, &Expr::nan()) || has(e, &Expr::oo()) || has(e, &Expr::neg_oo()) || has(e, &Expr::zoo())
}

thread_local! {
    static Z0: Expr = Expr::dummy("z0", Default::default());
}

/// The default `z0` dummy of `_hyperexpand`.
fn default_z0() -> Expr { Z0.with(|z| z.clone()) }

/// `make_derivative_operator(M, z)`.
fn derivative_operator<'a>(m: &'a [Vec<Expr>], z: &'a Expr) -> impl FnMut(&Vec<Expr>) -> HResult<Vec<Expr>> + 'a {
    move |c: &Vec<Expr>| {
        let n = c.len();
        let mut out = Vec::with_capacity(n);
        for j in 0..n {
            let mut terms = vec![mul(&[z.clone(), crate::diff::diff(&c[j], z)])];
            for (i, ci) in c.iter().enumerate() {
                terms.push(mul(&[ci.clone(), m[i][j].clone()]));
            }
            out.push(make_simp(&add(&terms), z)?);
        }
        Ok(out)
    }
}

fn zdiff(z: &Expr) -> impl FnMut(&Expr) -> HResult<Expr> + '_ { move |f: &Expr| Ok(mul(&[z.clone(), crate::diff::diff(f, z)])) }

fn apply_rewrite(e: &Expr, rw: Rewrite) -> Expr {
    match rw {
        Rewrite::Nonrep => hyperrep_rewrite(e, false),
        Rewrite::Nonrepsmall => hyperrep_rewrite(e, true),
        Rewrite::None | Rewrite::Default => e.clone(),
    }
}

/// `expr.replace(hyper, hyperexpand_special)`.
fn replace_hyper_special(e: &Expr) -> Expr {
    crate::simptools::bottom_up(e, &mut |x| match x.kind() {
        Kind::Func(Func::Special(SpecialFn::Hyper), a) => match (tuple_items(&a[0]), tuple_items(&a[1])) {
            (Some(ap), Some(bq)) => plan::hyperexpand_special(ap, bq, &a[2]),
            _ => x.clone(),
        },
        _ => x.clone(),
    })
}

fn has_hyper(e: &Expr) -> bool { crate::special::has_head(e, SpecialFn::Hyper) }

/// SymPy's `_hyperexpand(func, z, ops0, z0, premult, prem, rewrite)`.
pub(crate) fn hyperexpand_func(
    func: HF,
    z: &Expr,
    ops0: &[Op],
    z0: &Expr,
    premult: &Expr,
    prem: &Expr,
    rewrite: Rewrite,
) -> HResult<Expr> {
    if crate::assumptions::is_zero(z) == Some(true) {
        return Ok(Expr::one());
    }
    let z = polarify(z, false);
    let rewrite = if rewrite == Rewrite::Default { Rewrite::Nonrepsmall } else { rewrite };

    let carryout_plan = |f: &Formula, ops: &[Op]| -> HResult<Expr> {
        let zz = [(f.z.clone(), z0.clone())];
        let m = subs_mat(&f.m, &zz);
        let c = apply_operators(subs_all(&f.c, &zz), ops, &mut derivative_operator(&m, z0))?;
        let mp: Vec<Vec<Expr>> = m
            .iter()
            .enumerate()
            .map(|(i, r)| r.iter().enumerate().map(|(j, e)| if i == j { add(&[e.clone(), prem.clone()]) } else { e.clone() }).collect())
            .collect();
        let mut c = apply_operators(c, ops0, &mut derivative_operator(&mp, z0))?;
        if premult.is_one() {
            c = c.iter().map(|e| make_simp(e, z0)).collect::<HResult<_>>()?;
        }
        let b = subs_all(&f.b, &zz);
        let mut r = Expr::zero();
        for (ci, bi) in c.iter().zip(&b) {
            r = add(&[r, mul(&[ci.clone(), bi.clone()])]);
        }
        let r = mul(&[r, premult.clone()]);
        let res = crate::subs::subs(&r, &[(z0.clone(), z.clone())]);
        Ok(apply_rewrite(&res, rewrite))
    };

    let (func, mut ops) = plan::reduce_order(&func);
    if let Some(res) = plan::try_polynomial(&func, z0)? {
        let p = apply_operators(res, &ops, &mut zdiff(z0))?;
        let p = apply_operators(mul(&[p, premult.clone()]), ops0, &mut zdiff(z0))?;
        let p = crate::subs::subs(&crate::simplify::simplify(&p), &[(z0.clone(), z.clone())]);
        return Ok(unpolarify(&p, false));
    }
    let mut p = Expr::zero();
    let mut func = func;
    if let Some((f2, nops, p2)) = plan::try_shifted_sum(&func, z0)? {
        func = f2;
        ops.extend(nops);
        p = p2;
    }
    let p = apply_operators(p, &ops, &mut zdiff(z0))?;
    let p = apply_operators(mul(&[p, premult.clone()]), ops0, &mut zdiff(z0))?;
    let p = crate::subs::subs(&crate::simplify::simplify(&p), &[(z0.clone(), z.clone())]);

    let uz = unpolarify(&z, false);
    if (uz.is_one() || uz.is_minus_one()) && func.sizes() == (2, 1) {
        let f = plan::build_hypergeometric_formula(&func)?;
        let r = replace_hyper_special(&carryout_plan(&f, &ops)?);
        if !has_hyper(&r) {
            return Ok(add(&[r, p]));
        }
    }
    let formula = match formulae::lookup_origin(&func)? {
        Some(f) => f,
        None => match plan::try_lerchphi(&func)? {
            Some(f) => f,
            None => plan::build_hypergeometric_formula(&func)?,
        },
    };
    ops.extend(plan::devise_plan(&func, &formula.func, z0)?);
    let r = add(&[carryout_plan(&formula, &ops)?, p]);
    Ok(replace_hyper_special(&crate::integrate::meijerint::powdenest_polar(&r)))
}

/// A Slater condition: a Python bool or a relational.
fn is_true(c: &Expr) -> bool { matches!(c.kind(), Kind::Bool(true)) }
fn is_false(c: &Expr) -> bool { matches!(c.kind(), Kind::Bool(false)) }

fn polar_lift(e: &Expr) -> Expr { crate::special::apply(SpecialFn::PolarLift, vec![e.clone()]) }

fn count_head(e: &Expr, s: SpecialFn) -> usize {
    let own = usize::from(matches!(e.kind(), Kind::Func(Func::Special(g), _) if *g == s));
    if e.is_atom() {
        return own;
    }
    own + e.args().iter().map(|a| count_head(a, s)).sum::<usize>()
}

/// SymPy's `_meijergexpand(func, z0, allow_hyper, rewrite, place)`.
pub(crate) fn meijergexpand_func(func: GF, z0: &Expr, allow_hyper: bool, rewrite: Rewrite, place: Option<&Expr>) -> HResult<Expr> {
    let rewrite = if rewrite == Rewrite::Default { Rewrite::None } else { rewrite };
    let func0 = func.clone();
    let z = Expr::dummy("z", Default::default());
    let (func, mut ops) = plan::reduce_order_meijer(&func);
    if let Some(f) = formulae::meijer_lookup_origin(&func)? {
        ops.extend(plan::devise_plan_meijer(&f.func, &func, &z)?);
        let zz = [(f.z.clone(), z.clone())];
        let m = subs_mat(&f.m, &zz);
        let c = apply_operators(subs_all(&f.c, &zz), &ops, &mut derivative_operator(&m, &z))?;
        let c: Vec<Expr> = c.iter().map(|e| make_simp(e, &z)).collect::<HResult<_>>()?;
        let b = subs_all(&f.b, &zz);
        let mut r = Expr::zero();
        for (ci, bi) in c.iter().zip(&b) {
            r = add(&[r, mul(&[ci.clone(), bi.clone()])]);
        }
        let r = crate::subs::subs(&r, &[(z.clone(), z0.clone())]);
        return Ok(crate::integrate::meijerint::powdenest_polar(&r));
    }

    let t = Expr::dummy("t", Default::default());
    let (slater1, cond1) = plan::do_slater(&func.an, &func.bm, &func.ap, &func.bq, &z, z0, &ops, &t)?;
    let tr = |l: &[Expr]| l.iter().map(|x| sub(&Expr::one(), x)).collect::<Vec<_>>();
    let ops2: Vec<Op> = ops.iter().map(|o| o.reflected(&z, &div(&Expr::one(), &t))).collect();
    let (slater2, cond2) =
        plan::do_slater(&tr(&func.bm), &tr(&func.an), &tr(&func.bq), &tr(&func.ap), &t, &div(&Expr::one(), z0), &ops2, &t)?;
    let mut slater1 = crate::integrate::meijerint::powdenest_polar(&crate::subs::subs(&slater1, &[(z.clone(), z0.clone())]));
    let mut slater2 =
        crate::integrate::meijerint::powdenest_polar(&crate::subs::subs(&slater2, &[(t.clone(), div(&Expr::one(), z0))]));
    let (mut cond1, mut cond2) = (cond1, cond2);
    if !matches!(cond2.kind(), Kind::Bool(_)) {
        cond2 = crate::subs::subs(&cond2, &[(t.clone(), div(&Expr::one(), &z))]);
    }
    let p = func.an.len() + func.ap.len();
    let q = func.bm.len() + func.bq.len();
    // delta = m + n - (p + q)/2, compared through 2*delta.
    let delta2 = 2 * (func.bm.len() + func.an.len()) as i64 - (p + q) as i64;
    // ``(re(nu) < -1) is not False`` always holds: a SymPy relational is
    // never the Python ``False``.
    if delta2 > 0 || (delta2 == 0 && p == q && polar_lift(z0) == polar_lift(&Expr::one())) {
        if !is_false(&cond1) {
            cond1 = Expr::boolean(true);
        }
        if !is_false(&cond2) {
            cond2 = Expr::boolean(true);
        }
    }
    let pick = |c: &Expr| if is_true(c) { if rewrite == Rewrite::None { Rewrite::Nonrep } else { rewrite } } else if rewrite == Rewrite::None { Rewrite::Nonrepsmall } else { rewrite };
    slater1 = apply_rewrite(&slater1, pick(&cond1));
    slater2 = apply_rewrite(&slater2, pick(&cond2));
    if !is_false(&cond1) && !is_false(&cond2) {
        if let Some(pl) = place {
            if pl.is_zero() {
                cond2 = Expr::boolean(false);
            }
            if *pl == Expr::zoo() {
                cond1 = Expr::boolean(false);
            }
        }
    }
    if !matches!(cond1.kind(), Kind::Bool(_)) {
        cond1 = crate::subs::subs(&cond1, &[(z.clone(), z0.clone())]);
    }
    if !matches!(cond2.kind(), Kind::Bool(_)) {
        cond2 = crate::subs::subs(&cond2, &[(z.clone(), z0.clone())]);
    }
    let weight = |e: &Expr, c: &Expr| -> (u8, usize, usize) {
        let mut c0 = if is_true(c) { 0 } else if is_false(c) { 1 } else { 2 };
        if has_nonfinite(e) {
            c0 = 3;
        }
        (c0, count_head(e, SpecialFn::Hyper), crate::simptools::methods::count_ops(e))
    };
    let w1 = weight(&slater1, &cond1);
    let w2 = weight(&slater2, &cond2);
    let wmin = w1.min(w2);
    if wmin.0 == 0 && wmin.1 <= 1 {
        return Ok(if w1 < w2 { slater1 } else { slater2 });
    }
    let r = crate::relational::piecewise(vec![(slater1, cond1), (slater2, cond2), (func0.call(z0), Expr::boolean(true))]);
    if w1.0.max(w2.0) <= 1 && w1.1.max(w2.1) <= 1 {
        return Ok(r);
    }
    if !has_hyper(&r) || allow_hyper {
        return Ok(r);
    }
    Ok(func0.call(z0))
}

/// SymPy's `hyperexpand(f, allow_hyper, rewrite, place)`.
pub(crate) fn hyperexpand_with(f: &Expr, allow_hyper: bool, rewrite: Rewrite, place: Option<&Expr>) -> HResult<Expr> {
    let err: RefCell<Option<HxError>> = RefCell::new(None);
    let z0 = default_z0();
    let r = crate::simptools::bottom_up(f, &mut |x| {
        if err.borrow().is_some() {
            return x.clone();
        }
        if let Kind::Func(Func::Special(SpecialFn::Hyper), a) = x.kind() {
            if let (Some(ap), Some(bq)) = (tuple_items(&a[0]), tuple_items(&a[1])) {
                match hyperexpand_func(HF::new(ap.to_vec(), bq.to_vec()), &a[2], &[], &z0, &Expr::one(), &Expr::zero(), rewrite) {
                    Ok(r) => return r,
                    Err(e) => {
                        *err.borrow_mut() = Some(e);
                    }
                }
            }
        }
        x.clone()
    });
    if let Some(e) = err.into_inner() {
        return Err(e);
    }
    let err: RefCell<Option<HxError>> = RefCell::new(None);
    let r = crate::simptools::bottom_up(&r, &mut |x| {
        if err.borrow().is_some() {
            return x.clone();
        }
        if let Some(g) = meijer_parts(x) {
            match meijergexpand_func(g, &x.args()[2], allow_hyper, rewrite, place) {
                Ok(r) if !has_nonfinite(&r) => return r,
                Ok(_) => {}
                Err(e) => *err.borrow_mut() = Some(e),
            }
        }
        x.clone()
    });
    match err.into_inner() {
        Some(e) => Err(e),
        None => Ok(r),
    }
}

fn meijer_parts(x: &Expr) -> Option<GF> {
    let Kind::Func(Func::Special(SpecialFn::MeijerG), a) = x.kind() else { return None };
    let aa = tuple_items(&a[0])?;
    let bb = tuple_items(&a[1])?;
    Some(GF::new(
        tuple_items(&aa[0])?.to_vec(),
        tuple_items(&aa[1])?.to_vec(),
        tuple_items(&bb[0])?.to_vec(),
        tuple_items(&bb[1])?.to_vec(),
    ))
}

/// `hyperexpand(f)` with the default keywords; an expression the reference
/// cannot expand without raising is returned unchanged.
pub fn hyperexpand(f: &Expr) -> Expr { hyperexpand_with(f, false, Rewrite::Default, None).unwrap_or_else(|_| f.clone()) }

/// Public `hyperexpand(f, allow_hyper, rewrite, place)`. `rewrite` is the
/// keyword's string (`None` disables the rewrite); errors come back as
/// `(exception class name, message)`.
pub fn hyperexpand_options(
    f: &Expr,
    allow_hyper: bool,
    rewrite: Option<&str>,
    place: Option<&Expr>,
) -> Result<Expr, (&'static str, String)> {
    let rewrite = match rewrite {
        None => Rewrite::None,
        Some("nonrep") => Rewrite::Nonrep,
        Some("nonrepsmall") => Rewrite::Nonrepsmall,
        Some(_) => Rewrite::Default,
    };
    hyperexpand_with(f, allow_hyper, rewrite, place).map_err(|e| match e {
        HxError::Value(m) => ("ValueError", m),
        HxError::NotImplemented(m) => ("NotImplementedError", m),
        HxError::Type(m) => ("TypeError", m),
        HxError::Polynomial(m) => ("PolynomialError", m),
    })
}

fn rewrite_kw(e: &Expr) -> Rewrite {
    match e.kind() {
        Kind::Sym(s) => match s.name() {
            "nonrep" => Rewrite::Nonrep,
            "nonrepsmall" => Rewrite::Nonrepsmall,
            "none" => Rewrite::None,
            _ => Rewrite::Default,
        },
        _ => Rewrite::Default,
    }
}

fn place_kw(e: &Expr) -> Option<Expr> {
    match e.kind() {
        Kind::Sym(s) if s.name() == "none" => None,
        _ => Some(e.clone()),
    }
}

/// Internal probes compared with the reference implementation.
pub(crate) fn debug(op: &str, args: &[Expr]) -> Option<Vec<Expr>> {
    let wrap = |r: HResult<Expr>| match r {
        Ok(e) => vec![e],
        Err(e) => {
            let text = e.to_string();
            let (kind, msg) = text.split_once(':').unwrap_or((&text, ""));
            vec![Expr::symbol(&format!("EXC_{kind}")), Expr::symbol(msg.trim())]
        }
    };
    match (op, args) {
        ("hx", [e, allow, rw, place]) => {
            let place = place_kw(place);
            Some(wrap(hyperexpand_with(e, matches!(allow.kind(), Kind::Bool(true)), rewrite_kw(rw), place.as_ref())))
        }
        ("hx_hyper", [ap, bq, z, z0, premult, prem]) => {
            let f = HF::new(tuple_items(ap)?.to_vec(), tuple_items(bq)?.to_vec());
            Some(wrap(hyperexpand_func(f, z, &[], z0, premult, prem, Rewrite::None)))
        }
        ("hx_meijer", [an, ap, bm, bq, z0, allow, rw, place]) => {
            let g = GF::new(tuple_items(an)?.to_vec(), tuple_items(ap)?.to_vec(), tuple_items(bm)?.to_vec(), tuple_items(bq)?.to_vec());
            let place = place_kw(place);
            Some(wrap(meijergexpand_func(g, z0, matches!(allow.kind(), Kind::Bool(true)), rewrite_kw(rw), place.as_ref())))
        }
        ("hx_table", []) => Some(formulae::table_check()),
        _ => None,
    }
}
