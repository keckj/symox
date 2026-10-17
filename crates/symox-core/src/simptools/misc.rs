//! `signsimp`, `posify`, `separatevars`, `hypersimp` and `hypersimilar`
//! from `sympy.simplify.simplify`.
use super::{args, bottom_up, held_coeff_mul, map_args, mul_args, preorder};
use crate::expr::{Expr, Kind};
use crate::number::Number;
use crate::sortkey::cmp_default;
use crate::{assumptions, decompose, ops, traverse};

/// `expr.xreplace(rule)` over SymPy's public arguments: a node found in
/// the rule is replaced as a whole, other nodes are rebuilt with the
/// evaluating constructors.
pub fn xreplace(e: &Expr, rule: &[(Expr, Expr)]) -> Expr {
    if rule.is_empty() {
        return e.clone();
    }
    if let Some((_, v)) = rule.iter().find(|(k, _)| k == e) {
        return v.clone();
    }
    if e.is_atom() {
        return e.clone();
    }
    map_args(e, &mut |a| xreplace(a, rule))
}

/// All `Add` nodes of the public tree (`expr.atoms(Add)`), without duplicates.
fn add_atoms(e: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    preorder(e, &mut |a| {
        if a.is_add() && !out.contains(a) {
            out.push(a.clone());
        }
    });
    out
}

/// `sub_pre` and `sub_post` of `sympy.simplify.cse_opts`.
fn sub_pre_post(e: &Expr) -> Expr {
    let mut reps: Vec<(Expr, Expr)> = Vec::new();
    for a in add_atoms(e) {
        if decompose::could_extract_minus_sign(&a) {
            let na = ops::neg(&a);
            if na.is_mul() {
                continue;
            }
            reps.push((a, held_coeff_mul(&Number::Int(-1), &na)));
        }
    }
    let e = xreplace(e, &reps);
    // Second pass: the persisting Adds are marked `1*-1*(-a)`. A marker
    // absorbed by an evaluating parent behaves like `-1*(-a)`; a marker
    // left in place is turned back into `-(-a)` by `sub_post`.
    let mut atoms = add_atoms(&e);
    atoms.sort_by(cmp_default);
    let mut negs: Vec<(Expr, Expr)> = Vec::new();
    let mut markers: Vec<(Expr, Expr)> = Vec::new();
    for a in atoms {
        if let Some((_, r)) = reps.iter().find(|(k, _)| *k == a) {
            negs.push((a, r.clone()));
        } else if decompose::could_extract_minus_sign(&a) {
            let na = ops::neg(&a);
            let marker = held_coeff_mul(&Number::Int(-1), &na);
            markers.push((marker.clone(), ops::neg(&na)));
            negs.push((a, marker));
        }
    }
    let e = xreplace(&e, &negs);
    xreplace(&e, &markers)
}

/// `expr.replace(lambda x: x.is_Mul and --x != x, lambda x: --x)`.
fn double_negate_muls(e: &Expr) -> Expr {
    bottom_up(e, &mut |x| {
        if x.is_mul() {
            let y = ops::neg(&ops::neg(x));
            if y != *x {
                return y;
            }
        }
        x.clone()
    })
}

fn is_expr_like(e: &Expr) -> bool {
    match e.kind() {
        Kind::Bool(_) => false,
        Kind::Func(crate::expr::Func::Logic(_) | crate::expr::Func::Set(_), _) => false,
        _ => true,
    }
}

/// `signsimp(expr, evaluate)`.
pub fn signsimp(expr: &Expr, evaluate: bool) -> Expr {
    if !is_expr_like(expr) || expr.is_atom() {
        return expr.clone();
    }
    let e = double_negate_muls(expr);
    let e = sub_pre_post(&e);
    if !is_expr_like(&e) || e.is_atom() {
        return e;
    }
    if e.is_add() {
        let parts: Vec<Expr> = args(&e).iter().map(|a| signsimp(a, true)).collect();
        let rv = ops::add(&parts);
        if !evaluate && rv.is_add() && decompose::could_extract_minus_sign(&rv) {
            return held_coeff_mul(&Number::Int(-1), &ops::neg(&rv));
        }
        return rv;
    }
    if evaluate {
        return double_negate_muls(&e);
    }
    e
}

/// `posify(eq)`: replaces symbols whose sign is unknown by positive
/// dummies; returns the new expression and the `(dummy, symbol)` pairs.
pub fn posify(eq: &Expr) -> (Expr, Vec<(Expr, Expr)>) {
    let mut fwd: Vec<(Expr, Expr)> = Vec::new();
    let mut syms = traverse::free_symbols(eq);
    syms.sort_by(cmp_default);
    for s in syms {
        if assumptions::is_positive(&s) != None {
            continue;
        }
        let Some(sym) = s.as_sym() else { continue };
        let mut kw: Vec<(&str, bool)> = sym.facts().iter().filter(|(k, _)| *k != "positive").collect();
        kw.push(("positive", true));
        let Ok(facts) = crate::assumption_facts::FactSet::from_kwargs(kw) else { continue };
        let Ok(d) = crate::symbol::Symbol::dummy_with_facts(sym.name(), facts, None) else { continue };
        fwd.push((s.clone(), Expr::from_kind(Kind::Sym(d))));
    }
    let out = crate::subs::subs(eq, &fwd);
    (out, fwd.into_iter().map(|(s, d)| (d, s)).collect())
}

/// `posify` of a sequence: one set of dummies shared by all entries.
pub fn posify_seq(eqs: &[Expr]) -> (Vec<Expr>, Vec<(Expr, Expr)>) {
    let mut syms: Vec<Expr> = Vec::new();
    for e in eqs {
        for s in traverse::free_symbols(e) {
            if !syms.contains(&s) {
                syms.push(s);
            }
        }
    }
    let mut fwd: Vec<(Expr, Expr)> = Vec::new();
    for s in syms {
        let (_, r) = posify(&s);
        fwd.extend(r.into_iter().map(|(d, s)| (s, d)));
    }
    let out = eqs.iter().map(|e| crate::subs::subs(e, &fwd)).collect();
    (out, fwd.into_iter().map(|(s, d)| (d, s)).collect())
}

fn unposify(e: &Expr, reps: &[(Expr, Expr)]) -> Expr {
    if reps.is_empty() {
        e.clone()
    } else {
        crate::subs::subs(e, reps)
    }
}

fn factor_or_same(e: &Expr) -> Expr {
    crate::poly::api::factor(e).unwrap_or_else(|_| e.clone())
}

fn cset(e: &Expr) -> Vec<Expr> {
    let mut c = decompose::args_cnc(e, true).0;
    let mut out: Vec<Expr> = Vec::new();
    for x in c.drain(..) {
        if !out.contains(&x) {
            out.push(x);
        }
    }
    out
}

fn separatevars_inner(expr: &Expr, force: bool) -> Expr {
    if let Kind::Func(crate::expr::Func::Abs, a) = expr.kind() {
        let arg = &a[0];
        if arg.is_mul() && !super::is_num(arg) {
            return match separatevars_dict(arg, None, force) {
                Some(d) => ops::mul(&d.iter().map(|(_, v)| crate::functions::func(crate::expr::Func::Abs, vec![v.clone()])).collect::<Vec<_>>()),
                None => expr.clone(),
            };
        }
    }
    if traverse::free_symbols(expr).len() < 2 {
        return expr.clone();
    }
    if expr.is_mul() {
        // SymPy passes `force` positionally (as `symbols`), so the
        // recursion runs with force=False.
        return map_args(expr, &mut |a| separatevars_inner(a, false));
    }
    let mut expr = expr.clone();
    if let Kind::Pow(b, x) = expr.kind() {
        expr = ops::pow(&separatevars_inner(b, force), x);
    }
    let options = crate::expand::ExpandOptions { mul: false, multinomial: false, force, ..Default::default() };
    expr = crate::expand::expand_with_options(&expr, &options, &mut crate::expand::ExpandBudget::default()).unwrap_or(expr);
    let (pe, reps) = if force { posify(&expr) } else { (expr.clone(), vec![]) };
    expr = unposify(&factor_or_same(&pe), &reps);
    if !expr.is_add() {
        return expr;
    }
    let terms = args(&expr);
    let mut commonc = cset(&terms[0]);
    for t in &terms[1..] {
        let c = cset(t);
        commonc.retain(|x| c.contains(x));
    }
    let commonc = ops::mul(&commonc);
    let commonc = ops::as_coeff_mul(&commonc).1;
    let commonc_set = cset(&commonc);
    let parts: Vec<Expr> = terms
        .iter()
        .map(|a| {
            let (c, nc) = decompose::args_cnc(a, true);
            let mut cs: Vec<Expr> = Vec::new();
            for x in c {
                if !commonc_set.contains(&x) && !cs.contains(&x) {
                    cs.push(x);
                }
            }
            ops::mul(&[ops::mul(&cs), ops::mul(&nc)])
        })
        .collect();
    let mut nonsepar = ops::add(&parts);
    if traverse::free_symbols(&nonsepar).len() > 1 {
        let (pe, reps) = if force { posify(&nonsepar) } else { (nonsepar.clone(), vec![]) };
        let f = unposify(&factor_or_same(&pe), &reps);
        if !f.is_add() {
            nonsepar = f;
        }
    }
    ops::mul(&[commonc, nonsepar])
}

/// `separatevars(expr, force=force)`.
pub fn separatevars(expr: &Expr, force: bool) -> Expr {
    separatevars_inner(expr, force)
}

/// `separatevars(expr, symbols, dict=True, force=force)`. Keys are
/// `Some(symbol)` or `None` for `'coeff'`; `symbols=None` in SymPy is
/// `Some(&[])` with `none_symbols`. Returns `None` when the expression
/// does not separate.
pub fn separatevars_dict(expr: &Expr, symbols: Option<&[Expr]>, force: bool) -> Option<Vec<(Option<Expr>, Expr)>> {
    let expr = separatevars_inner(expr, force);
    separatevars_dict_of(&expr, symbols)
}

fn separatevars_dict_of(expr: &Expr, symbols: Option<&[Expr]>) -> Option<Vec<(Option<Expr>, Expr)>> {
    let symbols: Vec<Expr> = match symbols {
        Some(s) if !s.is_empty() => s.to_vec(),
        _ => {
            let mut s = traverse::free_symbols(expr);
            if s.is_empty() {
                return None;
            }
            s.sort_by(cmp_default);
            s
        }
    };
    let mut ret: Vec<(Option<Expr>, Vec<Expr>)> = symbols.iter().map(|s| (Some(s.clone()), vec![])).collect();
    ret.push((None, vec![]));
    for i in mul_args(expr) {
        let fs = traverse::free_symbols(&i);
        let inter: Vec<&Expr> = symbols.iter().filter(|s| fs.contains(s)).collect();
        if inter.len() > 1 {
            return None;
        }
        let key = inter.first().map(|s| (*s).clone());
        ret.iter_mut().find(|(k, _)| *k == key).unwrap().1.push(i);
    }
    Some(ret.into_iter().map(|(k, v)| (k, ops::mul(&v))).collect())
}

/// `separatevars(expr, symbols=None, dict=True)`: `{'coeff': expr}`.
pub fn separatevars_dict_none(expr: &Expr, force: bool) -> Vec<(Option<Expr>, Expr)> {
    vec![(None, separatevars_inner(expr, force))]
}

/// `hypersimilar(f, g, k)`.
pub fn hypersimilar(f: &Expr, g: &Expr, k: &Expr) -> Option<bool> {
    let h = super::gammasimp::rewrite_gamma(&ops::div(f, g), true);
    let h = crate::special::expand::expand_func(&h);
    let options = crate::expand::ExpandOptions { basic: false, ..Default::default() };
    let h = crate::expand::expand_with_options(&h, &options, &mut crate::expand::ExpandBudget::default()).unwrap_or(h);
    crate::inspection::classify(&h, std::slice::from_ref(k), crate::inspection::PolynomialKind::Rational)
}

/// `hypersimp(f, k)`: the simplified ratio `f(k+1)/f(k)` when it is a
/// rational function of `k`, else `None`.
pub fn hypersimp(f: &Expr, k: &Expr) -> Option<Expr> {
    let shifted = crate::subs::subs(f, &[(k.clone(), ops::add(&[k.clone(), Expr::one()]))]);
    let g = ops::div(&shifted, f);
    let g = super::gammasimp::rewrite_gamma(&g, true);
    // piecewise_fold(g).args[-1][0]: every Piecewise takes its last piece
    let g = bottom_up(&g, &mut |node| match node.kind() {
        Kind::Piecewise(pieces) => pieces.last().map(|(x, _)| x.clone()).unwrap_or_else(|| node.clone()),
        _ => node.clone(),
    });
    let g = crate::special::expand::expand_func(&g);
    let g = super::powsimp::powsimp(&g, &super::powsimp::PowsimpOptions { deep: true, combine: super::powsimp::Combine::Exp, force: false });
    if crate::inspection::classify(&g, std::slice::from_ref(k), crate::inspection::PolynomialKind::Rational) == Some(true) {
        Some(super::methods::simplify_rational(&g))
    } else {
        None
    }
}
