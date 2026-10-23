//! Simplification routines.
//!
//! This first version contains the building blocks other modules rely on:
//! `numer_denom`, `radsimp` (denominator rationalisation), `together`
//! and a light `simplify`. The polynomial based routines (`cancel`,
//! `factor`, ...) live in [`crate::poly`].

use crate::expr::{Expr, Func, Kind};
use rustc_hash::FxHashMap;
use crate::expand::expand;
use crate::ops::{add, div, mul, mul_factor, neg, pow};

/// Split `e` into `(numerator, denominator)` without combining sums
/// (SymPy's `as_numer_denom` for products and powers; sums are combined
/// over a common denominator).
pub fn numer_denom(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Num(n) => (Expr::num(n.numer()), Expr::num(n.denom())),
        Kind::Pow(b, x) => {
            if let Some(n) = x.as_number() {
                if n.is_negative() {
                    return (Expr::one(), pow(b, &Expr::num(n.neg())));
                }
            }
            if let Kind::Mul(c, _) = x.kind() {
                if c.is_negative() {
                    return (Expr::one(), pow(b, &neg(x)));
                }
            }
            (e.clone(), Expr::one())
        }
        Kind::Mul(c, facs) => {
            let mut num = vec![Expr::num(c.numer())];
            let mut den = vec![Expr::num(c.denom())];
            for (b, x) in facs {
                let (n, d) = numer_denom(&mul_factor(b, x));
                num.push(n);
                den.push(d);
            }
            (mul(&num), mul(&den))
        }
        Kind::Add(..) => {
            let parts: Vec<(Expr, Expr)> = e.args().iter().map(numer_denom).collect();
            if parts.iter().all(|(_, d)| d.is_one()) {
                return (e.clone(), Expr::one());
            }
            // common denominator: product of distinct denominators
            let mut dens: Vec<Expr> = Vec::new();
            for (_, d) in &parts {
                if !d.is_one() && !dens.contains(d) {
                    dens.push(d.clone());
                }
            }
            let den = mul(&dens);
            let mut nums = Vec::with_capacity(parts.len());
            for (n, d) in &parts {
                let mut f = vec![n.clone()];
                for dd in &dens {
                    if dd != d {
                        f.push(dd.clone());
                    }
                }
                nums.push(mul(&f));
            }
            (add(&nums), den)
        }
        _ => (e.clone(), Expr::one()),
    }
}

/// Combine a sum of fractions over a common denominator.
pub fn together(e: &Expr) -> Expr {
    let (n, d) = numer_denom(e);
    if d.is_one() {
        return n;
    }
    div(&n, &d)
}

fn has_sqrt(e: &Expr) -> bool {
    match e.kind() {
        Kind::Pow(b, x) => {
            (b.is_number() && *x == Expr::rational(1, 2)) || has_sqrt(b)
        }
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().any(has_sqrt),
        _ => false,
    }
}

/// Rationalise denominators of the form `a + b*sqrt(c)`.
pub fn radsimp(e: &Expr) -> Expr {
    let (n, d) = numer_denom(e);
    if d.is_one() || !has_sqrt(&d) {
        return e.clone();
    }
    let mut num = n;
    let mut den = d;
    for _ in 0..4 {
        if !has_sqrt(&den) {
            break;
        }
        let Kind::Add(c, terms) = den.kind() else {
            break;
        };
        // conjugate: flip the sign of the radical terms
        let mut conj_terms = Vec::new();
        let mut flipped = false;
        for (t, k) in terms {
            if !flipped && has_sqrt(t) {
                conj_terms.push(mul(&[Expr::num(k.neg()), t.clone()]));
                flipped = true;
            } else {
                conj_terms.push(mul(&[Expr::num(k.clone()), t.clone()]));
            }
        }
        conj_terms.push(Expr::num(c.clone()));
        let conj = add(&conj_terms);
        num = expand(&mul(&[num, conj.clone()]));
        den = expand(&mul(&[den, conj]));
    }
    if has_sqrt(&den) {
        return e.clone();
    }
    expand(&div(&num, &den))
}

/// SymPy's simplify candidates from the polynomial module: `cancel(e)`
/// and `factor_terms(cancel(e))` (common content and monomials pulled
/// out, a bare `number*(sum)` product being distributed again).
fn poly_candidates(e: &Expr) -> Vec<Expr> {
    let Ok(c) = crate::poly::api::cancel(e) else {
        return vec![];
    };
    let ft = crate::poly::api::together(&c);
    let ft = match ft.kind() {
        Kind::Mul(k, facs) if !k.is_one() && facs.len() == 1 && facs[0].0.is_add() && facs[0].1.is_one() => {
            mul(&[Expr::num(k.clone()), facs[0].0.clone()])
        }
        _ => ft,
    };
    vec![c, ft]
}

/// Release evaluation barriers before algebraic simplification. Expansion
/// alone deliberately keeps UnevaluatedExpr, matching SymPy's behavior.
fn release_unevaluated(e: &Expr, cache: &mut FxHashMap<Expr, Expr>) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    if let Some(cached) = cache.get(e) {
        return cached.clone();
    }
    let result = match e.kind() {
        Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") && args.len() == 1 => {
            release_unevaluated(&args[0], cache)
        }
        _ => crate::traverse::map_children(e, &mut |child| release_unevaluated(child, cache)),
    };
    cache.insert(e.clone(), result.clone());
    result
}

/// `expr.has(CombinatorialFunction, gamma)`.
fn has_gamma_or_combinatorial(e: &Expr) -> bool {
    use crate::special::SpecialFn::*;
    let mut found = false;
    crate::traverse::preorder(e, &mut |x| {
        if !found {
            found = matches!(
                x.kind(),
                Kind::Func(Func::Gamma, _)
                    | Kind::Func(
                        Func::Special(Factorial | Factorial2 | Subfactorial | Binomial | RisingFactorial | FallingFactorial),
                        _
                    )
            );
        }
    });
    found
}

/// SymPy's `hollow_mul` transform: a two-argument product of a number and
/// a sum is rebuilt with `Mul(*args)`, which distributes the number.
fn distribute_hollow_mul(e: &Expr) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    let e = crate::traverse::map_children(e, &mut |child| distribute_hollow_mul(child));
    match e.kind() {
        Kind::Mul(k, facs) if !k.is_one() && facs.len() == 1 && facs[0].0.is_add() && facs[0].1.is_one() => {
            mul(&[Expr::num(k.clone()), facs[0].0.clone()])
        }
        _ => e,
    }
}

/// The expression the reference `simplify` holds when it applies
/// `factor_terms(expr, sign=False)` unconditionally: the shortest of
/// `together(expr, deep=True)`, the cancelled forms and the
/// content-primitive input. When the native selection factors to the same
/// expression, both searches are in the same state and the factored form is
/// kept although it has more operations, as
/// `2*sqrt(T)*sqrt(k)*sqrt(R1 + R2)` for `sqrt(4*R1*T*k + 4*R2*T*k)` with
/// positive symbols. `None` when a step fails.
fn reference_factored_state(e: &Expr) -> Option<Expr> {
    let factored = distribute_hollow_mul(&crate::simptools::exprtools::factor_terms(&reference_state(e)?, false, false, false, false));
    // The later shorter steps (cancel, factor_terms and
    // expand_power_exp(expand_mul)) undo a factorization that has more
    // operations than its expansion: n**3 - n stays expanded.
    let short = shorter(&[factored.clone(), crate::poly::api::reference_cancel(&factored).ok()?]);
    let short = shorter(&[
        short.clone(),
        crate::simptools::exprtools::factor_terms(&short, false, false, false, true),
        crate::simp::expand_power_exp(&crate::simp::expand_mul(&short)),
    ]);
    (short == factored).then_some(factored)
}

/// SymPy's `shorter`: the first of the choices with the fewest operations.
fn shorter(choices: &[Expr]) -> Expr {
    use crate::simptools::methods::count_ops;
    let mut best = choices[0].clone();
    let mut best_ops = count_ops(&best);
    for c in &choices[1..] {
        let ops = count_ops(c);
        if ops < best_ops {
            best = c.clone();
            best_ops = ops;
        }
    }
    best
}

fn has_log(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |x| found |= matches!(x.kind(), Kind::Func(Func::Log, _)));
    found
}

/// `expr.has(TrigonometricFunction, HyperbolicFunction)`.
fn has_trig_or_hyperbolic(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |x| {
        if !found {
            found = matches!(
                x.kind(),
                Kind::Func(
                    Func::Sin
                        | Func::Cos
                        | Func::Tan
                        | Func::Cot
                        | Func::Sec
                        | Func::Csc
                        | Func::Sinh
                        | Func::Cosh
                        | Func::Tanh
                        | Func::Coth,
                    _
                )
            );
        }
    });
    found
}

/// The candidate of SymPy's own sequence of rewrites for an expression of
/// sums, products and powers: the shortest state of the rational steps,
/// `factor_terms(sign=False)`, `trigsimp(deep=True)` (Fu's `futrig`),
/// `combsimp` for gamma and combinatorial functions, the `powsimp`, `cancel`, `factor_terms` and `expand_power_exp(expand_mul)`
/// shorter steps, `exptrigsimp` and the `hollow_mul` pass. The native search
/// can expand factored powers of a sum with trigonometric terms into
/// multiple-angle sums, where this sequence keeps them factored, as SymPy
/// does: `(x + cos(y))**k` stays a power. `None` when a step fails.
fn reference_sequence(e: &Expr) -> Option<Expr> {
    use crate::simptools::powsimp::{Combine, PowsimpOptions};
    let mut expr = crate::simptools::exprtools::factor_terms(&reference_state(e)?, false, false, false, false);
    if has_trig_or_hyperbolic(&expr) {
        expr = crate::simptools::trigsimp::futrig(&expr, true);
    }
    // `if expr.has(log): expr = shorter(expand_log(expr, deep=True),
    // logcombine(expr))`, which does not compare with `expr` itself
    // (log(4)/2 becomes log(2)).
    if has_log(&expr) {
        expr = shorter(&[crate::simp::expand_log(&expr, false), crate::simp::logcombine(&expr, false)]);
    }
    // `if expr.has(CombinatorialFunction, gamma): expr = combsimp(expr)`.
    if has_gamma_or_combinatorial(&expr) {
        expr = crate::simptools::gammasimp::combsimp(&expr);
    }
    let short = shorter(&[
        crate::simptools::powsimp::powsimp(&expr, &PowsimpOptions { deep: true, combine: Combine::Exp, force: false }),
        crate::simptools::powsimp::powsimp(&expr, &Default::default()),
        expr,
    ]);
    let short = shorter(&[short.clone(), crate::poly::api::reference_cancel(&short).ok()?]);
    let short = shorter(&[
        short.clone(),
        crate::simptools::exprtools::factor_terms(&short, false, false, false, true),
        crate::simp::expand_power_exp(&crate::simp::expand_mul(&short)),
    ]);
    let short = if has_exp_or_trig(&short) { crate::simptools::trigsimp::exptrigsimp(&short) } else { short };
    Some(distribute_hollow_mul(&short))
}

/// The shortest state of the reference `simplify` after its rational
/// steps (signsimp, `normal` bottom up, content, cancel and together),
/// before `factor_terms`.
fn reference_state(e: &Expr) -> Option<Expr> {
    // The cancellation of the bottom-up normal form runs under the factor
    // budget: algebraic gcds over many dependent numeric radicals
    // (12**(1/6), 12**(5/6), sqrt(3), sqrt(229) in a Ferrari root residual)
    // can blow up, and this state is only one candidate among several.
    crate::poly::factor_budget::run(crate::poly::factor_budget::Limits::default(), || {
        reference_state_unbudgeted(e).ok_or_else(|| crate::poly::convert::PolyError::Unsupported("no reference state".into()))
    })
    .ok()
}

fn reference_state_unbudgeted(e: &Expr) -> Option<Expr> {
    let start = crate::simptools::misc::signsimp(e, true);
    let start = crate::simptools::bottom_up(&start, &mut |w| crate::simptools::normal(w));
    let (content, primitive) =
        crate::decompose::as_content_primitive(&crate::simptools::powsimp::powsimp(&start, &Default::default()), false, false).ok()?;
    let start = mul(&[Expr::num(content), primitive]);
    let cancelled = crate::poly::api::reference_cancel(&start).ok()?;
    let expanded = crate::poly::api::reference_cancel(&crate::simptools::mexpand(&cancelled)).ok()?;
    let expr1 = shorter(&[cancelled, expanded]);
    let expr2 = shorter(&[crate::poly::api::together(&start), crate::poly::api::together(&expr1)]);
    Some(shorter(&[expr2, expr1, start]))
}

/// `expr.has(TrigonometricFunction, HyperbolicFunction, ExpBase, exp)`.
fn has_exp_or_trig(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |x| {
        if !found {
            found = matches!(
                x.kind(),
                Kind::Func(
                    Func::Exp
                        | Func::Sin
                        | Func::Cos
                        | Func::Tan
                        | Func::Cot
                        | Func::Sec
                        | Func::Csc
                        | Func::Sinh
                        | Func::Cosh
                        | Func::Tanh
                        | Func::Coth,
                    _
                )
            );
        }
    });
    found
}

/// A light general purpose simplifier (to be extended).
pub fn simplify(e: &Expr) -> Expr {
    let evaluated = release_unevaluated(e, &mut FxHashMap::default());
    if evaluated.is_atom() {
        return evaluated;
    }
    // `BesselBase._eval_simplify` is besselsimp; other expressions with
    // Bessel functions pass through besselsimp before the final choice.
    if crate::special::expand::is_bessel_base(&evaluated) {
        return crate::special::besselsimp(&evaluated);
    }
    let evaluated = if crate::special::expand::has_bessel(&evaluated) {
        crate::special::besselsimp(&evaluated)
    } else {
        evaluated
    };
    let e = &evaluated;
    if e.is_atom() {
        return e.clone();
    }
    // SymPy's simplify returns 0 at once for an expression known to be
    // zero; a numeric value decides it exactly (root objects included).
    if crate::functions::is_numeric_value(e) && crate::assumptions::is_zero(e) == Some(true) {
        return Expr::zero();
    }
    let mut trig_memo = crate::simp::TrigMemo::default();
    let t = crate::simp::trigsimp_cached(e, &mut trig_memo);
    // Zero is already the smallest possible result. In particular, do not
    // run polynomial conversion and Fu after proving a trig identity.
    if t.is_zero() {
        return t;
    }
    // Prefer the cancelled rational form when equally small, including
    // its deterministic positive leading denominator coefficient.
    let mut candidates = poly_candidates(e);
    let cancelled = candidates.first().cloned();
    if candidates.iter().any(Expr::is_zero) {
        return Expr::zero();
    }
    let expanded = expand(e);
    if expanded.is_zero() {
        return expanded;
    }
    let expanded_trig = if expanded == *e { t.clone() } else { crate::simp::trigsimp_cached(&expanded, &mut trig_memo) };
    if expanded_trig.is_zero() {
        return expanded_trig;
    }
    candidates.extend([
        e.clone(),
        expanded,
        together(e),
        crate::simp::powsimp(e),
        t.clone(),
        expanded_trig,
    ]);
    if has_gamma_or_combinatorial(e) {
        // SymPy's simplify passes gamma and combinatorial functions through
        // combsimp (the duplication and reflection formulas).
        candidates.push(crate::simptools::gammasimp::combsimp(e));
    }
    if crate::simp::has_trig(e) {
        // SymPy's simplify also tries the Fu-based trigsimp, which rewrites
        // tan/cot/sec/csc through sin and cos (tan(x)*cos(x) -> sin(x)).
        candidates.push(crate::simptools::trigsimp::futrig(e, true));
    }
    // SymPy's own sequence keeps the factored powers that the expansion
    // based candidates lose. It comes last, so a tie keeps the native form.
    let reference = reference_sequence(e);
    candidates.extend(reference.clone());
    let selected = candidates
        .into_iter()
        // Ties in size go to SymPy's count_ops measure, which counts an
        // n-factor product as n - 1 multiplications (3*exp(a)*exp(b) is
        // larger than 3*exp(a + b)).
        .min_by_key(|c| (crate::traverse::node_count(c), crate::simptools::methods::count_ops(c)))
        .unwrap();
    // SymPy's shorter() keeps the cancelled fraction on a count_ops tie: a
    // numerator sign pulled out (-(2*x + y)/d for (-2*x - y)/d) has fewer
    // nodes but is not shorter.
    let selected = match cancelled {
        Some(c) if c != selected && crate::simptools::methods::count_ops(&c) <= crate::simptools::methods::count_ops(&selected) => {
            let (n1, d1) = numer_denom(&selected);
            let (n0, d0) = numer_denom(&c);
            if d1 == d0 && d0.is_add() && expand(&n1) == expand(&n0) { c } else { selected }
        }
        _ => selected,
    };
    // The cancelled form splits the powers of a sum exponent
    // (C1*exp(-1)*exp(-x) for C1*exp(-x - 1)); SymPy's shorter() keeps the
    // input on a count_ops tie, so a selection that only splits the powers
    // of the input (powsimp gives the input back) must be shorter.
    let selected = if selected != *e
        && crate::simptools::methods::count_ops(e) <= crate::simptools::methods::count_ops(&selected)
        && crate::simp::powsimp(&selected) == *e
    {
        e.clone()
    } else {
        selected
    };
    // SymPy then applies factor_terms(expr, sign=False) (x - x*y is
    // x*(1 - y)); its final hollow_mul pass distributes a remaining
    // number*(sum) product again. SymPy applies it unconditionally, but the
    // native selection above is not SymPy's chain of rewrites, so the
    // result is kept only when it is not larger and keeps a proven real
    // value provable: factoring the radicals out of a Cardano conjugate pair
    // z**(7/3) + z**(-7/3) would hide that the root is real.
    let factored = distribute_hollow_mul(&crate::simptools::exprtools::factor_terms(&selected, false, false, false, false));
    let keeps_realness = || {
        !crate::functions::contains_i(&selected)
            || crate::assumptions::is_real(&selected) != Some(true)
            || crate::assumptions::is_real(&factored) == Some(true)
    };
    // On a count_ops tie, SymPy's own sequence decides: when it ends on the
    // selected form, its shorter steps undid the factorization
    // (-2*x*f(x) + x*exp(-x**2) stays a sum, which solve's form-dependent
    // factorable classification of dsolve relies on).
    let factored_ops = crate::simptools::methods::count_ops(&factored);
    let selected_ops = crate::simptools::methods::count_ops(&selected);
    let tie_keeps_reference = factored_ops == selected_ops && reference.as_ref() == Some(&selected);
    let selected = if ((factored_ops <= selected_ops && !tie_keeps_reference)
        || (factored != selected && reference_factored_state(e).as_ref() == Some(&factored)))
        && keeps_realness()
    {
        factored
    } else {
        selected
    };
    // SymPy applies powsimp after its rational and expansion steps, so the
    // powers in an expanded or cancelled form are combined as well.
    let mut selected = selected;
    // The light pass merges equal bases and exponentials; SymPy's full
    // powsimp also merges equal exponents (sqrt(2)*sqrt(a + b) is
    // sqrt(2*a + 2*b)).
    for combined in [
        crate::simp::powsimp(&selected),
        crate::simptools::powsimp::powsimp(&selected, &Default::default()),
    ] {
        // SymPy compares these with count_ops alone: a negated combined
        // power (-sqrt(2*a + 2*b)) has more nodes but fewer operations.
        let ops = |c: &Expr| (crate::simptools::methods::count_ops(c), crate::traverse::node_count(c));
        if combined != selected && ops(&combined) < ops(&selected) {
            selected = combined;
        }
    }
    // SymPy then applies exptrigsimp, without a size comparison, when the
    // shortest form has trigonometric, hyperbolic or exponential functions:
    // im(exp(2*I)/(1 - exp(4*I))) becomes 1/(2*sin(2)). Its final ratio test
    // against the input is applied by the public wrapper.
    if has_exp_or_trig(&selected) {
        selected = crate::simptools::trigsimp::exptrigsimp(&selected);
    }
    // SymPy reaches its sign through two steps. The input first passes
    // through signsimp, which negates a sum with more negative terms (else
    // the canonical order of the sum and its negation decides). Its cancel
    // candidate then makes the denominator's leading coefficient positive in
    // the polynomial generator order (x, y, z before other names). The
    // smaller of the two by count_ops wins, the signsimp form on a tie:
    // b/(b*x - a*y) keeps its sign, e11/(e11*e20 - e10*e21) and x/(1 - x)
    // flip, and (C1 + x)/(C1 - x) is kept although x leads C1.
    let (numerator, denominator) = numer_denom(&selected);
    if denominator.is_add() {
        let flipped = || div(&neg(&numerator), &neg(&denominator));
        let signsimp_form =
            if crate::decompose::could_extract_minus_sign(&denominator) { flipped() } else { selected.clone() };
        let cancel_negative = crate::poly::convert::from_exprs(std::slice::from_ref(&denominator))
            .ok()
            .map(|c| c.polys[0].lc().is_negative())
            .unwrap_or(false);
        let cancel_form = if cancel_negative { flipped() } else { selected.clone() };
        let ops = crate::simptools::methods::count_ops;
        return if ops(&cancel_form) < ops(&signsimp_form) { cancel_form } else { signsimp_form };
    }
    selected
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn barycentric_ratio_denominator_sign() {
        let expr = crate::parser::parse_str("e11/(e11*e20-e10*e21)").unwrap();
        let expected = crate::parser::parse_str("-e11/(e10*e21-e11*e20)").unwrap();
        assert_eq!(simplify(&expr), expected);
        let a = Expr::symbol("e₀₀");
        let b = Expr::symbol("e₀₁");
        let c = Expr::symbol("e₂₀");
        let d = Expr::symbol("e₂₁");
        let determinant = crate::ops::sub(&mul(&[a, d]), &mul(&[b.clone(), c]));
        let ratio = div(&neg(&b), &neg(&determinant));
        assert_eq!(simplify(&ratio), div(&b, &determinant));
    }

    #[test]
    fn simplify_releases_unevaluated_expr() {
        let x = Expr::symbol("x");
        let inner = add(&[x.clone(), Expr::one()]);
        let held = crate::functions::undefined("UnevaluatedExpr", vec![inner.clone()]);
        assert_eq!(expand(&held), held);
        assert_eq!(simplify(&crate::ops::sub(&held, &inner)), Expr::zero());
        assert_eq!(crate::diff::diff(&held, &x), crate::diff::derivative_unevaluated(&held, vec![(x, 1)]));
    }
}
