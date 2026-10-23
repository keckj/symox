//! `sympy.simplify.gammasimp` and `sympy.simplify.combsimp`.
use super::{args, as_numer_denom, bottom_up, map_args, mul_args, ordered};
use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::number::Number;
use crate::special::SpecialFn;
use crate::{assumptions, decompose, ops, traverse};

fn gamma(x: &Expr) -> Expr {
    func(Func::Gamma, vec![x.clone()])
}

fn special(s: SpecialFn, a: Vec<Expr>) -> Expr {
    func(Func::Special(s), a)
}

fn is_gamma(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Gamma, _))
}

fn gamma_arg(e: &Expr) -> Option<&Expr> {
    match e.kind() {
        Kind::Func(Func::Gamma, a) => Some(&a[0]),
        _ => None,
    }
}

fn has_gamma(e: &Expr) -> bool {
    count_gamma(e) > 0
}

fn count_gamma(e: &Expr) -> usize {
    let mut n = 0;
    walk(e, &mut |a| {
        if is_gamma(a) {
            n += 1;
        }
    });
    n
}

/// Pre-order traversal that also enters `Piecewise` pieces.
fn walk(e: &Expr, f: &mut dyn FnMut(&Expr)) {
    if let Kind::Piecewise(pieces) = e.kind() {
        f(e);
        for (x, c) in pieces {
            walk(x, f);
            walk(c, f);
        }
        return;
    }
    f(e);
    if matches!(e.kind(), Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Func(..) | Kind::Rel(..)) {
        for a in args(e) {
            walk(&a, f);
        }
    }
}

/// Bottom-up rebuild that also enters `Piecewise` pieces.
fn replace(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    if let Kind::Piecewise(pieces) = e.kind() {
        let new: Vec<(Expr, Expr)> = pieces.iter().map(|(x, c)| (replace(x, f), replace(c, f))).collect();
        return f(&crate::relational::piecewise(new));
    }
    if e.is_atom() {
        return f(e);
    }
    let rebuilt = map_args(e, &mut |c| replace(c, f));
    f(&rebuilt)
}

/// SymPy's `isinstance(e, Function)` for the heads symox models.
fn is_function(e: &Expr) -> bool {
    match e.kind() {
        Kind::Piecewise(_) => true,
        Kind::Func(f, _) => match f {
            Func::Logic(_) | Func::Set(_) | Func::Order | Func::Integral | Func::Sum | Func::Subs | Func::Limit | Func::Max | Func::Min => false,
            Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..) => false,
            Func::Undefined(n) => !(n.is_internal("UnevaluatedExpr") || n.is_internal("Tuple")),
            _ => true,
        },
        _ => false,
    }
}

fn function_atoms(e: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    walk(e, &mut |a| {
        if is_function(a) && !out.contains(a) {
            out.push(a.clone());
        }
    });
    out
}

/// `expr.rewrite(gamma, piecewise=piecewise)` for the combinatorial functions.
pub fn rewrite_gamma(e: &Expr, piecewise: bool) -> Expr {
    replace(e, &mut |node| rewrite_gamma_node(node, piecewise).unwrap_or_else(|| node.clone()))
}

fn rewrite_gamma_node(e: &Expr, piecewise: bool) -> Option<Expr> {
    let Kind::Func(Func::Special(s), a) = e.kind() else { return None };
    let one = Expr::one();
    Some(match s {
        SpecialFn::Factorial => gamma(&ops::add(&[a[0].clone(), one])),
        SpecialFn::Binomial => {
            let (n, k) = (&a[0], &a[1]);
            ops::div(
                &gamma(&ops::add(&[n.clone(), one.clone()])),
                &ops::mul(&[gamma(&ops::add(&[k.clone(), one.clone()])), gamma(&ops::add(&[ops::sub(n, k), one]))]),
            )
        }
        SpecialFn::RisingFactorial => {
            let (x, k) = (&a[0], &a[1]);
            let main = ops::div(&gamma(&ops::add(&[x.clone(), k.clone()])), &gamma(x));
            let other = ops::mul(&[
                ops::pow(&Expr::minus_one(), k),
                ops::div(&gamma(&ops::sub(&one, x)), &gamma(&ops::add(&[ops::neg(k), ops::neg(x), one.clone()]))),
            ]);
            if !piecewise {
                if assumptions::is_nonpositive(x) == Some(true) {
                    return Some(other);
                }
                return Some(main);
            }
            crate::relational::piecewise(vec![(main, crate::relational::gt(x, &Expr::zero())), (other, Expr::boolean(true))])
        }
        SpecialFn::FallingFactorial => {
            let (x, k) = (&a[0], &a[1]);
            let main = ops::div(&gamma(&ops::add(&[x.clone(), one.clone()])), &gamma(&ops::add(&[ops::sub(x, k), one])));
            let other = ops::mul(&[ops::pow(&Expr::minus_one(), k), ops::div(&gamma(&ops::sub(k, x)), &gamma(&ops::neg(x)))]);
            if !piecewise {
                if assumptions::is_negative(x) == Some(true) {
                    return Some(other);
                }
                return Some(main);
            }
            crate::relational::piecewise(vec![(main, crate::relational::ge(x, &Expr::zero())), (other, Expr::boolean(true))])
        }
        SpecialFn::Beta => {
            let x = a[0].clone();
            let y = a.get(1).cloned().unwrap_or_else(|| x.clone());
            ops::div(&ops::mul(&[gamma(&x), gamma(&y)]), &gamma(&ops::add(&[x, y])))
        }
        SpecialFn::Catalan => {
            let n = &a[0];
            let half = Expr::rational(1, 2);
            ops::div(
                &ops::mul(&[ops::pow(&Expr::int(4), n), gamma(&ops::add(&[n.clone(), half.clone()]))]),
                &ops::mul(&[gamma(&half), gamma(&ops::add(&[n.clone(), Expr::int(2)]))]),
            )
        }
        _ => return None,
    })
}

/// `expr.rewrite(factorial)` for `gamma`.
fn rewrite_factorial(e: &Expr) -> Expr {
    replace(e, &mut |node| match gamma_arg(node) {
        Some(z) => special(SpecialFn::Factorial, vec![ops::sub(z, &Expr::one())]),
        None => node.clone(),
    })
}

fn full_expand(e: &Expr) -> Expr {
    crate::expand::expand(e)
}

fn int_of(n: &Number) -> Option<i64> {
    if n.is_integer() { n.to_i64() } else { None }
}

/// `_rf(1, b)` evaluated as far as SymPy's `_rf.eval` goes, with the
/// remaining `_rf(1, X)` written as `gamma(X + 1)`.
fn rf_one(b: &Expr) -> Expr {
    let one = Expr::one();
    if let Some(n) = b.as_number().and_then(int_of) {
        return rf_int(&one, n);
    }
    if b.is_add() {
        let (c, rest) = ops::as_coeff_add(b);
        if let Some(c) = int_of(&c) {
            let head = gamma(&ops::add(&[rest.clone(), one.clone()]));
            if c > 0 {
                return ops::mul(&[head, rf_int(&ops::add(&[one, rest]), c)]);
            } else if c < 0 {
                return ops::div(&head, &rf_int(&ops::add(&[one, rest, Expr::int(c)]), -c));
            }
        }
    }
    gamma(&ops::add(&[b.clone(), one]))
}

fn rf_int(a: &Expr, n: i64) -> Expr {
    if n == 0 {
        return Expr::one();
    }
    if n > 0 {
        ops::mul(&(0..n).map(|i| ops::add(&[a.clone(), Expr::int(i)])).collect::<Vec<_>>())
    } else {
        ops::inv(&ops::mul(&(1..=-n).map(|i| ops::sub(a, &Expr::int(i))).collect::<Vec<_>>()))
    }
}

/// `gamma(n) -> _rf(1, (n - 1).expand())` followed by the `_rf -> gamma` rewrite.
fn normalize_gammas(e: &Expr) -> Expr {
    replace(e, &mut |node| match gamma_arg(node) {
        Some(n) => rf_one(&full_expand(&ops::sub(n, &Expr::one()))),
        None => node.clone(),
    })
}

/// `factor(expr)`. SymPy's `factor` runs `together` on the whole product
/// before it expands each factor base into a polynomial; expanding the
/// bases first changes the result once numeric-base powers combine
/// (`2**s*6**s` is `12**s`).
fn factor(e: &Expr) -> Expr {
    crate::poly::api::factor(e).unwrap_or_else(|_| e.clone())
}

/// `expand_func(gamma(n))` for a Rational `n`.
fn expand_gamma_rational(node: &Expr, n: &Number) -> Expr {
    if let Some((p, q)) = n.numer_denom() {
        if !num_traits::One::is_one(&q) && num_traits::Signed::abs(&p) > q {
            let whole = num_integer::Integer::div_floor(&p, &q);
            let frac = Number::from_ratio(num_rational::BigRational::new(&p - &whole * &q, q.clone()));
            if let Some(whole) = num_traits::ToPrimitive::to_i64(&whole) {
                let f = Expr::num(frac);
                return ops::mul(&[gamma(&f), rf_int(&f, whole)]);
            }
        }
    }
    crate::special::expand::expand_func(node)
}

fn gamma_factor(x: &Expr) -> bool {
    if is_gamma(x) {
        return true;
    }
    if x.is_add() || x.is_mul() {
        return args(x).iter().any(gamma_factor);
    }
    if let Kind::Pow(b, e) = x.kind() {
        if assumptions::is_integer(e) == Some(true) || assumptions::is_positive(b) == Some(true) {
            return gamma_factor(b);
        }
    }
    false
}

fn gamma_rat(x: &Expr) -> Expr {
    let was = count_gamma(x);
    let xx = normalize_gammas(x);
    if count_gamma(&xx) < was { xx } else { x.clone() }
}

fn is_integer_expr(e: &Expr) -> Option<i64> {
    e.as_number().and_then(int_of)
}

fn remove_first(v: &mut Vec<Expr>, x: &Expr) {
    if let Some(i) = v.iter().position(|y| y == x) {
        v.remove(i);
    }
}

fn explicate(p: &Expr) -> (Option<bool>, Vec<Expr>) {
    if p.is_one() {
        return (None, vec![]);
    }
    let (b, e) = decompose::as_base_exp(p);
    if let Some(n) = is_integer_expr(&e) {
        let n = n.max(0) as usize;
        if let Some(a) = gamma_arg(&b) {
            (Some(true), vec![a.clone(); n])
        } else {
            (Some(false), vec![b; n])
        }
    } else {
        (Some(false), vec![p.clone()])
    }
}

fn reflect(gammas: &mut Vec<Expr>, numer: &mut Vec<Expr>, denom: &mut Vec<Expr>) {
    let mut new = Vec::new();
    while let Some(g1) = gammas.pop() {
        if assumptions::is_integer(&g1) == Some(true) {
            new.push(g1);
            continue;
        }
        let mut hit = false;
        for i in 0..gammas.len() {
            let n = ops::add(&[g1.clone(), gammas[i].clone(), Expr::minus_one()]);
            let Some(n) = is_integer_expr(&n) else { continue };
            numer.push(Expr::pi());
            denom.push(func(Func::Sin, vec![ops::mul(&[Expr::pi(), g1.clone()])]));
            gammas.remove(i);
            if n > 0 {
                numer.extend((0..n).map(|k| ops::add(&[Expr::one(), ops::neg(&g1), Expr::int(k)])));
            } else if n < 0 {
                denom.extend((0..-n).map(|k| ops::sub(&ops::neg(&g1), &Expr::int(k))));
            }
            hit = true;
            break;
        }
        if !hit {
            new.push(g1);
        }
    }
    *gammas = new;
}

fn duplicate(ng: &mut Vec<Expr>, dg: &mut Vec<Expr>, no: &mut Vec<Expr>, d_o: &mut Vec<Expr>) {
    loop {
        let mut found = None;
        'outer: for x in ng.iter() {
            for y in dg.iter() {
                let n = ops::sub(x, &ops::mul(&[Expr::int(2), y.clone()]));
                if let Some(n) = is_integer_expr(&n) {
                    found = Some((x.clone(), y.clone(), n));
                    break 'outer;
                }
            }
        }
        let Some((x, y, n)) = found else { break };
        remove_first(ng, &x);
        remove_first(dg, &y);
        let y2 = ops::mul(&[Expr::int(2), y.clone()]);
        if n > 0 {
            no.extend((0..n).map(|k| ops::add(&[y2.clone(), Expr::int(k)])));
        } else if n < 0 {
            d_o.extend((0..-n).map(|k| ops::sub(&ops::sub(&y2, &Expr::one()), &Expr::int(k))));
        }
        ng.push(ops::add(&[y.clone(), Expr::rational(1, 2)]));
        no.push(ops::pow(&Expr::int(2), &ops::sub(&y2, &Expr::one())));
        d_o.push(ops::sqrt(&Expr::pi()));
    }
}

fn frac_mod1(n: &Number) -> Number {
    n.sub(&n.floor())
}

/// `_run(coeffs)` of the multiplication theorem.
fn run(coeffs: &mut Vec<Number>) -> Option<(i64, Number, Vec<Number>)> {
    let mut u: Vec<Number> = Vec::new();
    for c in coeffs.iter() {
        if !u.contains(c) {
            u.push(c.clone());
        }
    }
    for i in 0..u.len() {
        let dj: Vec<(Number, usize)> = ((i + 1)..u.len()).map(|j| (frac_mod1(&u[j].sub(&u[i])), j)).collect();
        for (one, _) in &dj {
            let Some((p, q)) = one.numer_denom() else { continue };
            if !(num_traits::One::is_one(&p) && !num_traits::One::is_one(&q)) {
                continue;
            }
            let Some(n) = num_traits::ToPrimitive::to_i64(&q) else { continue };
            let mut got = vec![i];
            let mut get: Vec<i64> = (1..n).collect();
            let mut complete = false;
            for (d, j) in &dj {
                let m = d.mul(&Number::Int(n));
                if let Some(m) = int_of(&m) {
                    if let Some(pos) = get.iter().position(|g| *g == m) {
                        get.remove(pos);
                        got.push(*j);
                        if get.is_empty() {
                            complete = true;
                            break;
                        }
                    }
                }
            }
            if !complete {
                continue;
            }
            let mut vals = Vec::new();
            for j in got {
                let c = u[j].clone();
                if let Some(pos) = coeffs.iter().position(|x| *x == c) {
                    coeffs.remove(pos);
                }
                vals.push(c);
            }
            let first = vals.remove(0);
            return Some((n, first, vals));
        }
    }
    None
}

fn mult_thm(gammas: &mut Vec<Expr>, numer: &mut Vec<Expr>) {
    let mut rats: Vec<(Expr, Vec<Number>)> = Vec::new();
    for g in gammas.iter() {
        let (c, resid) = ops::as_coeff_add(g);
        match rats.iter_mut().find(|(r, _)| *r == resid) {
            Some((_, v)) => v.push(c),
            None => rats.push((resid, vec![c])),
        }
    }
    rats.sort_by(|a, b| crate::sortkey::cmp_default(&a.0, &b.0));
    let mut out = Vec::new();
    for (resid, coeffs) in rats.iter_mut() {
        coeffs.sort_by(|a, b| a.cmp_value(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut new = Vec::new();
        while let Some((n, ui, other)) = run(coeffs) {
            for u in &other {
                let con = ops::add(&[resid.clone(), Expr::num(u.clone()), Expr::minus_one()]);
                let count = u.sub(&ui).to_f64() as i64;
                for k in 0..count.max(0) {
                    numer.push(ops::sub(&con, &Expr::int(k)));
                }
            }
            let con = ops::mul(&[Expr::int(n), ops::add(&[resid.clone(), Expr::num(ui.clone())])]);
            numer.push(ops::mul(&[
                ops::pow(&ops::mul(&[Expr::int(2), Expr::pi()]), &Expr::rational(n - 1, 2)),
                ops::pow(&Expr::int(n), &ops::sub(&Expr::rational(1, 2), &con)),
            ]));
            new.push(con);
        }
        out.extend(coeffs.iter().map(|c| ops::add(&[resid.clone(), Expr::num(c.clone())])));
        out.extend(new);
    }
    *gammas = out;
}

fn compute_st(e: &Expr) -> (Vec<Expr>, Vec<Expr>) {
    let free = traverse::free_symbols(e);
    let mut t: Vec<Expr> = Vec::new();
    walk(e, &mut |a| {
        let item = if is_function(a) {
            Some(a.clone())
        } else if let Kind::Pow(_, x) = a.kind() {
            Some(x.clone())
        } else {
            None
        };
        if let Some(item) = item {
            if !t.contains(&item) {
                t.push(item);
            }
        }
    });
    (free, t)
}

fn same_set(a: &[Expr], b: &[Expr]) -> bool {
    a.len() == b.len() && a.iter().all(|x| b.contains(x))
}

fn find_fuzzy(l: &[Expr], x: &Expr) -> Option<Expr> {
    if l.is_empty() {
        return None;
    }
    let (s1, t1) = compute_st(x);
    for y in l {
        let (s2, t2) = compute_st(y);
        let meet = s1.iter().any(|s| s2.contains(s));
        if !same_set(&t1, &t2) || (!meet && (!s1.is_empty() || !s2.is_empty())) {
            continue;
        }
        let q = ops::div(x, y);
        let q = crate::poly::api::cancel(&q).unwrap_or(q);
        let a = traverse::free_symbols(&q).len();
        let b = s1.len();
        let c = s2.len();
        if a == 0 && (b > 0 || c > 0) {
            return Some(y.clone());
        }
    }
    None
}

fn fuzzy(gammas: &mut Vec<Expr>, numer: &mut Vec<Expr>, denom: &mut Vec<Expr>) {
    let mut new = Vec::new();
    while let Some(mut g) = gammas.pop() {
        let mut cont = true;
        while cont {
            cont = false;
            if let Some(y) = find_fuzzy(numer, &g) {
                remove_first(numer, &y);
                if y != g {
                    numer.push(ops::div(&y, &g));
                }
                g = ops::add(&[g, Expr::one()]);
                cont = true;
            }
            let gm1 = ops::sub(&g, &Expr::one());
            if let Some(y) = find_fuzzy(denom, &gm1) {
                remove_first(denom, &y);
                if y != gm1 {
                    numer.push(ops::div(&gm1, &y));
                }
                g = gm1;
                cont = true;
            }
        }
        new.push(g);
    }
    *gammas = new;
}

fn rule_gamma(expr: &Expr, mut level: u32, as_comb: bool) -> Expr {
    if expr.is_atom() {
        return expr.clone();
    }
    let mut expr = expr.clone();
    if level == 0 {
        expr = if let Kind::Piecewise(pieces) = expr.kind() {
            crate::relational::piecewise(pieces.iter().map(|(x, c)| (rule_gamma(x, 1, as_comb), c.clone())).collect())
        } else {
            map_args(&expr, &mut |x| rule_gamma(x, 1, as_comb))
        };
        level = 1;
    }
    if !expr.is_mul() {
        return expr;
    }
    if level == 1 {
        let margs = args(&expr);
        let (c, nc): (Vec<Expr>, Vec<Expr>) = margs.into_iter().partition(|a| decompose::commutativity(a) != Some(false));
        if c.is_empty() {
            return expr;
        }
        if !nc.is_empty() {
            return ops::mul(&[rule_gamma(&ops::mul(&c), 2, as_comb), ops::mul(&nc)]);
        }
        level = 2;
    }
    if level == 2 {
        let (t, f): (Vec<Expr>, Vec<Expr>) = args(&expr).into_iter().partition(gamma_factor);
        let gamma_ind = ops::mul(&f);
        let d = ops::mul(&t);
        let (mut nd, mut dd) = as_numer_denom(&d);
        for ipass in 0..2 {
            let mut a = ordered(mul_args(&nd));
            for i in 0..a.len() {
                if a[i].is_add() {
                    let terms: Vec<Expr> = args(&a[i]).iter().map(|t| rule_gamma(&gamma_rat(&ops::div(t, &dd)), 3, as_comb)).collect();
                    let (ni, d2) = as_numer_denom(&ops::add(&terms));
                    a[i] = ni;
                    dd = d2;
                    if !has_gamma(&dd) {
                        break;
                    }
                }
            }
            nd = ops::mul(&a);
            if ipass == 0 && !gamma_factor(&nd) {
                break;
            }
            std::mem::swap(&mut nd, &mut dd);
        }
        expr = ops::div(&ops::mul(&[gamma_ind, nd.clone()]), &dd);
        if !(expr.is_mul() && (gamma_factor(&dd) || gamma_factor(&nd))) {
            return expr;
        }
        level = 3;
    }
    if level == 3 {
        loop {
            let was = expr.clone();
            expr = rule_gamma(&expr, 4, as_comb);
            if expr == was {
                return expr;
            }
        }
    }
    let mut ng = Vec::new();
    let mut dg = Vec::new();
    let mut no = Vec::new();
    let mut d_o = Vec::new();
    let mut newargs = ordered(args(&expr));
    while let Some(p) = newargs.pop() {
        let (n, d) = as_numer_denom(&p);
        match explicate(&n) {
            (Some(true), l) => ng.extend(l),
            (Some(false), l) => no.extend(l),
            _ => {}
        }
        match explicate(&d) {
            (Some(true), l) => dg.extend(l),
            (Some(false), l) => d_o.extend(l),
            _ => {}
        }
    }
    if !as_comb {
        reflect(&mut ng, &mut no, &mut d_o);
        reflect(&mut dg, &mut d_o, &mut no);
        duplicate(&mut ng, &mut dg, &mut no, &mut d_o);
        duplicate(&mut dg, &mut ng, &mut d_o, &mut no);
        mult_thm(&mut ng, &mut no);
        mult_thm(&mut dg, &mut d_o);
    }
    fuzzy(&mut ng, &mut no, &mut d_o);
    fuzzy(&mut dg, &mut d_o, &mut no);
    let num = ops::mul(&ng.iter().map(gamma).collect::<Vec<_>>());
    let den = ops::mul(&dg.iter().map(gamma).collect::<Vec<_>>());
    ops::div(&ops::mul(&[ops::div(&num, &den), ops::mul(&no)]), &ops::mul(&d_o))
}

fn gammasimp_core(expr: &Expr, as_comb: bool, level: u32) -> Expr {
    let expr = normalize_gammas(expr);
    let was = factor(&expr);
    let mut out = rule_gamma(&was, level, as_comb);
    if out != was {
        out = factor(&out);
    }
    replace(&out, &mut |node| match gamma_arg(node) {
        Some(n) if n.is_rational_number() => expand_gamma_rational(node, n.as_number().unwrap()),
        _ => node.clone(),
    })
}

/// `_gammasimp(expr, as_comb)`.
fn gammasimp_inner(expr: &Expr, as_comb: bool) -> Expr {
    if let Kind::Piecewise(pieces) = expr.kind() {
        // the arguments are ExprCondPairs, simplified one level down
        let new = pieces.iter().map(|(x, c)| (gammasimp_core(x, as_comb, 1), c.clone())).collect();
        return crate::relational::piecewise(new);
    }
    gammasimp_core(expr, as_comb, 0)
}

/// `gammasimp(expr)`.
pub fn gammasimp(e: &Expr) -> Expr {
    let expr = rewrite_gamma(e, true);
    let f = function_atoms(&expr);
    if !f.iter().any(is_gamma) {
        return expr;
    }
    let others: Vec<Expr> = ordered(f.into_iter().filter(|a| !is_gamma(a)).collect());
    if others.is_empty() {
        return gammasimp_inner(&expr, false);
    }
    let mut to_dummy = Vec::new();
    let mut back = Vec::new();
    for fi in &others {
        let d = Expr::dummy("Dummy", crate::symbol::Assumptions(0));
        let simp = if let Kind::Piecewise(_) = fi.kind() {
            gammasimp_inner(fi, false)
        } else {
            map_args(fi, &mut |a| gammasimp_inner(a, false))
        };
        to_dummy.push((fi.clone(), d.clone()));
        back.push((d, simp));
    }
    let d = super::misc::xreplace(&expr, &to_dummy);
    super::misc::xreplace(&gammasimp_inner(&d, false), &back)
}

/// `_gamma_as_comb(expr)`.
fn gamma_as_comb(expr: &Expr) -> Expr {
    let expr = rewrite_factorial(expr);
    bottom_up(&expr, &mut |rv| {
        if !rv.is_mul() {
            return rv.clone();
        }
        let mut rvd: Vec<(Expr, Expr)> = Vec::new();
        for t in args(rv) {
            let (b, e) = decompose::as_base_exp(&t);
            match rvd.iter_mut().find(|(k, _)| *k == b) {
                Some((_, v)) => *v = ops::add(&[v.clone(), e]),
                None => rvd.push((b, e)),
            }
        }
        let mut nd: [Vec<Expr>; 2] = [Vec::new(), Vec::new()];
        for (k, v) in rvd.iter_mut() {
            if let Kind::Func(Func::Special(SpecialFn::Factorial), a) = k.kind() {
                if let Some(n) = is_integer_expr(v) {
                    if n > 0 {
                        nd[0].extend(std::iter::repeat_n(a[0].clone(), n as usize));
                    } else {
                        nd[1].extend(std::iter::repeat_n(a[0].clone(), (-n) as usize));
                    }
                    *v = Expr::zero();
                }
            }
        }
        if nd[0].is_empty() || nd[1].is_empty() {
            return rv.clone();
        }
        let mut hit = false;
        for m in 0..2 {
            let mut i = 0;
            while i < nd[m].len() {
                let ai = nd[m][i].clone();
                let mut matched = false;
                for j in (i + 1)..nd[m].len() {
                    let aj = nd[m][j].clone();
                    let sum = ops::add(&[ai.clone(), aj.clone()]);
                    if nd[1 - m].contains(&sum) {
                        hit = true;
                        remove_first(&mut nd[1 - m], &sum);
                        nd[m].remove(j);
                        nd[m].remove(i);
                        let low = if super::methods::count_ops(&ai) < super::methods::count_ops(&aj) { ai.clone() } else { aj.clone() };
                        let key = special(SpecialFn::Binomial, vec![sum, low]);
                        let delta = Expr::int(if m == 0 { -1 } else { 1 });
                        match rvd.iter_mut().find(|(k, _)| *k == key) {
                            Some((_, v)) => *v = ops::add(&[v.clone(), delta]),
                            None => rvd.push((key, delta)),
                        }
                        matched = true;
                        break;
                    }
                }
                if !matched {
                    i += 1;
                }
            }
        }
        if !hit {
            return rv.clone();
        }
        let mut num: Vec<Expr> = rvd.iter().map(|(k, v)| ops::pow(k, v)).collect();
        num.extend(nd[0].iter().map(|k| special(SpecialFn::Factorial, vec![k.clone()])));
        let den: Vec<Expr> = nd[1].iter().map(|k| special(SpecialFn::Factorial, vec![k.clone()])).collect();
        ops::div(&ops::mul(&num), &ops::mul(&den))
    })
}

/// `combsimp(expr)`.
pub fn combsimp(e: &Expr) -> Expr {
    let expr = rewrite_gamma(e, false);
    let mut non_integer = false;
    walk(&expr, &mut |a| {
        if let Some(n) = gamma_arg(a) {
            if assumptions::is_integer(n) != Some(true) {
                non_integer = true;
            }
        }
    });
    if non_integer {
        return gammasimp(&expr);
    }
    let expr = gammasimp_inner(&expr, true);
    gamma_as_comb(&expr)
}
