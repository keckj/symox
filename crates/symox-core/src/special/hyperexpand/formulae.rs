//! The formula collections: SymPy's `FormulaCollection` and
//! `MeijerFormulaCollection`, their instantiation and origin lookup.

use std::rc::Rc;

use rustc_hash::FxHashMap;

use crate::ops::{add, sub};
use crate::Expr;

use super::table::{RawHyper, RawMeijer, HYPER, MEIJER};
use super::{bucket, half, has_nonfinite, mod1, num_cmp, sift, sift_mod1, subs_all, subs_mat, Formula, HResult, HxError, Invariants, GF, HF};

#[derive(Clone, Copy, Debug)]
enum Matcher {
    UpperGamma,
    Si3113,
}

/// `MeijerFormula`.
#[derive(Clone, Debug)]
pub(crate) struct MeijerFormula {
    pub func: GF,
    pub z: Expr,
    pub b: Vec<Expr>,
    pub c: Vec<Expr>,
    pub m: Vec<Vec<Expr>>,
    matcher: Option<Matcher>,
}

struct Tables {
    hyper: Vec<Formula>,
    /// Concrete formulae, keyed by sizes and invariants; a later formula
    /// replaces an earlier one with the same key.
    concrete: Vec<((usize, usize), Invariants, usize)>,
    symbolic: Vec<((usize, usize), Vec<usize>)>,
    meijer: Vec<MeijerFormula>,
    dummies: Dummies,
}

#[derive(Clone)]
struct Dummies {
    a: Expr,
    rho: Expr,
}

thread_local! {
    static TABLES: Rc<Tables> = Rc::new(build());
}

fn tables() -> Rc<Tables> { TABLES.with(Rc::clone) }

fn parse(s: &str, map: &FxHashMap<Expr, Expr>) -> Expr {
    let e = crate::parser::parse_str(s).unwrap_or_else(|err| panic!("hyperexpand table entry {s:?}: {err:?}"));
    crate::subs::xreplace(&e, map)
}

fn build() -> Tables {
    let dummy = |n: &str| Expr::dummy(n, Default::default());
    let (a, b, c, z, rho) = (dummy("a"), dummy("b"), dummy("c"), dummy("z"), dummy("rho"));
    let mut map: FxHashMap<Expr, Expr> = FxHashMap::default();
    for (name, d) in [("hxa", &a), ("hxb", &b), ("hxc", &c), ("hxz", &z), ("hxrho", &rho)] {
        map.insert(Expr::symbol(name), d.clone());
    }
    let pv = |v: &[&str]| v.iter().map(|s| parse(s, &map)).collect::<Vec<_>>();
    let pm = |m: &[&[&str]]| m.iter().map(|r| pv(r)).collect::<Vec<_>>();
    let hyper: Vec<Formula> = HYPER
        .iter()
        .map(|r: &RawHyper| Formula {
            func: HF::new(pv(r.ap), pv(r.bq)),
            z: z.clone(),
            symbols: pv(r.symbols),
            b: pv(r.b),
            c: pv(r.c),
            m: pm(r.m),
        })
        .collect();
    let mut concrete: Vec<((usize, usize), Invariants, usize)> = Vec::new();
    let mut symbolic: Vec<((usize, usize), Vec<usize>)> = Vec::new();
    for (i, f) in hyper.iter().enumerate() {
        let sizes = f.func.sizes();
        if !f.symbols.is_empty() {
            match symbolic.iter_mut().find(|(s, _)| *s == sizes) {
                Some((_, v)) => v.push(i),
                None => symbolic.push((sizes, vec![i])),
            }
        } else {
            let inv = f.func.build_invariants();
            match concrete.iter_mut().find(|(s, iv, _)| *s == sizes && *iv == inv) {
                Some(slot) => slot.2 = i,
                None => concrete.push((sizes, inv, i)),
            }
        }
    }
    let meijer = MEIJER
        .iter()
        .map(|r: &RawMeijer| MeijerFormula {
            func: GF::new(pv(r.an), pv(r.ap), pv(r.bm), pv(r.bq)),
            z: z.clone(),
            b: pv(r.b),
            c: pv(r.c),
            m: pm(r.m),
            matcher: Some(match r.matcher {
                "detect_uppergamma" => Matcher::UpperGamma,
                _ => Matcher::Si3113,
            }),
        })
        .collect();
    Tables { hyper, concrete, symbolic, meijer, dummies: Dummies { a, rho } }
}

/// Every closed form of the table (a parse check for the debug probe).
pub(crate) fn table_check() -> Vec<Expr> {
    let t = tables();
    let mut out: Vec<Expr> = t.hyper.iter().map(|f| crate::special::tuple(f.b.clone())).collect();
    out.extend(t.meijer.iter().map(|f| crate::special::tuple(f.b.clone())));
    out
}

/// `solve(e, n)` for an `e` linear in `n`, as `find_instantiations` uses it.
fn solve_linear(e: &Expr, n: &Expr) -> HResult<Expr> {
    let e = crate::expand::expand(e);
    let k = crate::expand::expand(&crate::diff::diff(&e, n));
    if k.is_zero() || crate::traverse::has(&k, n) {
        return Err(HxError::Value("not enough values to unpack".into()));
    }
    let c0 = crate::subs::subs(&e, &[(n.clone(), Expr::zero())]);
    Ok(crate::simplify::simplify(&crate::ops::neg(&crate::ops::div(&c0, &k))))
}

fn xreplace_pairs(e: &Expr, repl: &[(Expr, Expr)]) -> Expr {
    let map: FxHashMap<Expr, Expr> = repl.iter().cloned().collect();
    crate::subs::xreplace(e, &map)
}

fn floor_ceil(v: &Expr, ceil: bool) -> HResult<i64> {
    let r = if ceil { crate::functions::ceiling(v) } else { crate::functions::floor(v) };
    r.as_i64().ok_or_else(|| HxError::Type("cannot determine truth value of Relational".into()))
}

fn num_min_max(v: &[Expr]) -> HResult<(Expr, Expr)> {
    let mut lo = v[0].clone();
    let mut hi = v[0].clone();
    for x in &v[1..] {
        match (num_cmp(&sub(x, &lo)), num_cmp(&sub(x, &hi))) {
            (Some(a), Some(b)) => {
                if a == std::cmp::Ordering::Less {
                    lo = x.clone();
                }
                if b == std::cmp::Ordering::Greater {
                    hi = x.clone();
                }
            }
            _ => return Err(HxError::Type("cannot determine truth value of Relational".into())),
        }
    }
    Ok((lo, hi))
}

/// `Formula.find_instantiations(func)`.
fn find_instantiations(f: &Formula, func: &HF) -> HResult<Vec<Vec<(Expr, Expr)>>> {
    if func.ap.len() != f.func.ap.len() || func.bq.len() != f.func.bq.len() {
        return Err(HxError::Type("Cannot instantiate other number of parameters".into()));
    }
    let mut symbol_values: Vec<&Vec<Expr>> = Vec::new();
    for a in &f.symbols {
        if f.func.ap.contains(a) {
            symbol_values.push(&func.ap);
        } else if f.func.bq.contains(a) {
            symbol_values.push(&func.bq);
        } else {
            return Err(HxError::Value(format!("At least one of the parameters of the formula must be equal to {a}")));
        }
    }
    let base_repl = product(&symbol_values);
    let (abuckets, bbuckets) = (sift_mod1(&func.ap), sift_mod1(&func.bq));
    let mut critical: Vec<Vec<Expr>> = f.symbols.iter().map(|_| vec![Expr::zero()]).collect();
    let mut result = Vec::new();
    let n = Expr::dummy("n", Default::default());
    for values in base_repl {
        let repl: Vec<(Expr, Expr)> = f.symbols.iter().cloned().zip(values).collect();
        let symb_a = sift(&f.func.ap, &mut |x| mod1(&xreplace_pairs(x, &repl)));
        let symb_b = sift(&f.func.bq, &mut |x| mod1(&xreplace_pairs(x, &repl)));
        for (bk, obk) in [(&abuckets, &symb_a), (&bbuckets, &symb_b)] {
            let mut keys: Vec<Expr> = bk.iter().map(|(k, _)| k.clone()).collect();
            for (k, _) in obk.iter() {
                if !keys.contains(k) {
                    keys.push(k.clone());
                }
            }
            for m in &keys {
                let (Some(targets), Some(oexprs)) = (bucket(bk, m), bucket(obk, m)) else { break };
                if targets.len() != oexprs.len() {
                    break;
                }
                for (si, (a, a0)) in repl.iter().enumerate() {
                    if !crate::traverse::free_symbols(a0).is_empty() {
                        continue;
                    }
                    let repl0: Vec<(Expr, Expr)> =
                        repl.iter().map(|(s, v)| if s == a { (s.clone(), add(&[v.clone(), n.clone()])) } else { (s.clone(), v.clone()) }).collect();
                    for expr in oexprs.iter().filter(|e| crate::traverse::has(e, a)) {
                        for target in targets {
                            let n0 = solve_linear(&sub(&xreplace_pairs(expr, &repl0), target), &n)?;
                            if !crate::traverse::free_symbols(&n0).is_empty() {
                                return Err(HxError::Value("Value should not be true".into()));
                            }
                            critical[si].push(n0);
                        }
                    }
                }
            }
        }
        let mut per_symbol: Vec<Vec<Expr>> = Vec::new();
        for ((_, a0), vals) in repl.iter().zip(&critical) {
            let (lo, hi) = num_min_max(vals)?;
            let (lo, hi) = (floor_ceil(&lo, false)?, floor_ceil(&hi, true)?);
            per_symbol.push((lo..=hi).map(|k| add(&[a0.clone(), Expr::int(k)])).collect());
        }
        let refs: Vec<&Vec<Expr>> = per_symbol.iter().collect();
        for l in product(&refs) {
            result.push(f.symbols.iter().cloned().zip(l).collect());
        }
    }
    Ok(result)
}

/// `itertools.product(*lists)`.
fn product(lists: &[&Vec<Expr>]) -> Vec<Vec<Expr>> {
    let mut out: Vec<Vec<Expr>> = vec![Vec::new()];
    for l in lists {
        let mut next = Vec::with_capacity(out.len() * l.len());
        for prefix in &out {
            for x in l.iter() {
                let mut p = prefix.clone();
                p.push(x.clone());
                next.push(p);
            }
        }
        out = next;
    }
    out
}

/// `subs(dict)`: keys ordered by node count (descending), then by
/// `default_sort_key`.
fn ordered_pairs(repl: &[(Expr, Expr)]) -> Vec<(Expr, Expr)> {
    let mut v = repl.to_vec();
    v.sort_by(|a, b| {
        crate::traverse::node_count(&b.0).cmp(&crate::traverse::node_count(&a.0)).then_with(|| crate::sortkey::cmp_default(&a.0, &b.0))
    });
    v
}

/// `FormulaCollection.lookup_origin(func)`.
pub(crate) fn lookup_origin(func: &HF) -> HResult<Option<Formula>> {
    let t = tables();
    let inv = func.build_invariants();
    let sizes = func.sizes();
    if let Some((_, _, i)) = t.concrete.iter().find(|(s, iv, _)| *s == sizes && *iv == inv) {
        return Ok(Some(t.hyper[*i].clone()));
    }
    let Some((_, list)) = t.symbolic.iter().find(|(s, _)| *s == sizes) else { return Ok(None) };
    let mut possible: Vec<(Expr, Vec<(Expr, Expr)>, usize, HF)> = Vec::new();
    for &i in list {
        let f = &t.hyper[i];
        for repl in find_instantiations(f, func)? {
            let func2 = f.func.xreplace(&repl);
            if !func2.is_suitable_origin() {
                continue;
            }
            let Some(diff) = func2.difficulty(func)? else { continue };
            possible.push((diff, repl, i, func2));
        }
    }
    let mut err = None;
    possible.sort_by(|a, b| match num_cmp(&sub(&a.0, &b.0)) {
        Some(o) => o,
        None => {
            err.get_or_insert(HxError::Type("cannot determine truth value of Relational".into()));
            std::cmp::Ordering::Equal
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    for (_, repl, i, func2) in possible {
        let f = &t.hyper[i];
        let pairs = ordered_pairs(&repl);
        let f2 = Formula { func: func2, z: f.z.clone(), symbols: vec![], b: subs_all(&f.b, &pairs), c: subs_all(&f.c, &pairs), m: subs_mat(&f.m, &pairs) };
        if !(f2.b.iter().chain(&f2.c).chain(f2.m.iter().flatten()).any(has_nonfinite)) {
            return Ok(Some(f2));
        }
    }
    Ok(None)
}

fn mod1s(e: &Expr) -> Expr { mod1(&crate::simplify::simplify(e)) }

fn gt0(e: &Expr) -> HResult<bool> {
    match num_cmp(e) {
        Some(o) => Ok(o == std::cmp::Ordering::Greater),
        None => Err(HxError::Type("cannot determine truth value of Relational".into())),
    }
}

fn detect_uppergamma(func: &GF, d: &Dummies) -> HResult<Option<(Vec<(Expr, Expr)>, GF)>> {
    let x = func.an[0].clone();
    let (mut y, mut z) = (func.bm[0].clone(), func.bm[1].clone());
    let mut swapped = false;
    if mod1s(&sub(&x, &y)).is_zero() {
        swapped = true;
        std::mem::swap(&mut y, &mut z);
    }
    if !mod1s(&sub(&x, &z)).is_zero() || gt0(&sub(&x, &z))? {
        return Ok(None);
    }
    let l = if swapped { vec![x.clone(), y.clone()] } else { vec![y.clone(), x.clone()] };
    Ok(Some((vec![(d.rho.clone(), y.clone()), (d.a.clone(), sub(&x, &y))], GF::new(vec![x], vec![], l, vec![]))))
}

fn detect_3113(func: &GF, d: &Dummies) -> HResult<Option<(Vec<(Expr, Expr)>, GF)>> {
    let x = func.an[0].clone();
    let (u, v, w) = (func.bm[0].clone(), func.bm[1].clone(), func.bm[2].clone());
    let (sig, x1, x2, y);
    if mod1s(&sub(&u, &v)).is_zero() {
        if mod1s(&sub(&v, &w)).is_zero() {
            return Ok(None);
        }
        sig = [half(), half(), Expr::zero()];
        (x1, x2, y) = (u, v, w);
    } else if mod1s(&sub(&x, &u)).is_zero() {
        sig = [half(), Expr::zero(), half()];
        (x1, y, x2) = (u, v, w);
    } else {
        sig = [Expr::zero(), half(), half()];
        (y, x1, x2) = (u, v, w);
    }
    if !mod1s(&sub(&x, &x1)).is_zero()
        || !mod1s(&sub(&x, &x2)).is_zero()
        || mod1s(&sub(&x, &y)) != half()
        || gt0(&sub(&x, &x1))?
        || gt0(&sub(&x, &x2))?
    {
        return Ok(None);
    }
    let bm = sig.iter().map(|t| add(&[x.clone(), crate::ops::neg(&half()), t.clone()])).collect();
    Ok(Some((vec![(d.a.clone(), x.clone())], GF::new(vec![x], vec![], bm, vec![]))))
}

/// `MeijerFormulaCollection.lookup_origin(func)`.
pub(crate) fn meijer_lookup_origin(func: &GF) -> HResult<Option<MeijerFormula>> {
    let t = tables();
    for f in t.meijer.iter().filter(|f| f.func.signature() == func.signature()) {
        let res = match f.matcher {
            Some(Matcher::UpperGamma) => detect_uppergamma(func, &t.dummies)?,
            Some(Matcher::Si3113) => detect_3113(func, &t.dummies)?,
            None => None,
        };
        if let Some((repl, newfunc)) = res {
            let pairs = ordered_pairs(&repl);
            return Ok(Some(MeijerFormula {
                func: newfunc,
                z: f.z.clone(),
                b: subs_all(&f.b, &pairs),
                c: subs_all(&f.c, &pairs),
                m: subs_mat(&f.m, &pairs),
                matcher: None,
            }));
        }
    }
    Ok(None)
}
