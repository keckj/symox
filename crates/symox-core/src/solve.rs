//! Equation solving: `solve` (univariate, systems), `nsolve`.

use rustc_hash::FxHashMap;

use crate::diff::diff;
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::functions::*;
use crate::integrate::poly_coeffs;
use crate::number::Number;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sqrt, sub};
use crate::simplify::{numer_denom, together};
use crate::subs::{subs, subs1};
use crate::traverse::{free_symbols, has};

mod algebraic;
mod recast;
mod rational_coverage;
pub use rational_coverage::complete_rational_roots;
mod trig;
pub use recast::{solve_atom,solve_atom_system};
pub(crate) mod lambert;
pub mod tsolve;
pub mod checksol;
mod symbolic_power;
pub mod unrad;

/// `lhs - rhs` for an equation (`Eq`) or the expression itself.
pub fn as_expr(eq: &Expr) -> Expr {
    match eq.kind() {
        Kind::Rel(RelOp::Eq, a, b) => sub(a, b),
        Kind::Bool(value) => if *value { Expr::zero() } else { Expr::one() },
        _ => eq.clone(),
    }
}

fn sort_solutions(v: &mut Vec<Expr>) {
    let mut seen: Vec<Expr> = Vec::new();
    v.retain(|e| {
        if seen.contains(e) {
            false
        } else {
            seen.push(e.clone());
            true
        }
    });
    v.sort_by_cached_key(crate::sortkey::sort_key);
}

/// Solve `eq = 0` for `x`. Returns `None` when no method applies.
pub fn solve(eq: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    solve_with_check(eq, x, true)
}

pub fn solve_with_check(eq: &Expr, x: &Expr, check: bool) -> Option<Vec<Expr>> {
    if let Kind::Bool(b) = eq.kind() {
        return if *b { None } else { Some(vec![]) };
    }
    // The reference `solve` expands trig angles once, on the whole equation.
    let f = trig::expand_angles(&as_expr(eq), x);
    let mut sols = solve_checked(&f, x, check)?;
    // The reference `solve` rewrites Eq(lhs, rhs) as the unevaluated
    // Add(lhs, -rhs) and rejects the solutions that make a denominator of
    // any of its subexpressions zero (`_simple_dens`), including the
    // denominators that cancel in lhs - rhs: Eq(x + 1/x, 1/x) has none.
    if check {
        if let Kind::Rel(RelOp::Eq, lhs, rhs) = eq.kind() {
            let mut dens = simple_dens(lhs, x);
            for d in simple_dens(rhs, x) {
                if !dens.contains(&d) { dens.push(d); }
            }
            if !dens.is_empty() {
                sols.retain(|s| !dens.iter().any(|d| denominator_vanishes(d, x, s)));
            }
        }
    }
    Some(sols)
}

fn solve_checked(f: &Expr, x: &Expr, check: bool) -> Option<Vec<Expr>> {
    let f = f.clone();
    // The reference `_solve` solves each factor of a product on its own (so
    // each factor keeps its own root forms), then rejects the solutions that
    // make a denominator of the product zero.
    let product = f.is_mul() && has(&f, x);
    let mut sols = if product {
        let mut sols = Vec::new();
        for factor in f.args() {
            if matches!(factor.kind(), Kind::Const(c) if c.is_infinite()) {
                return Some(vec![]);
            }
            if has(&factor, x) {
                sols.extend(solve_checked(&factor, x, check)?);
            }
        }
        sols
    } else {
        solve_rec(&f, x, 0)?
    };
    // Infinities reached by log(0), for example, are limiting values rather
    // than finite solutions. This restriction also applies with check=False.
    // A candidate such as `2*atan(I) = oo*I` (from the tan half-angle
    // rewrite) is not a finite solution either.
    sols.retain(|s| !matches!(s.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN) && !has_infinity(s));
    // An equation with an infinite coefficient (1 - oo*x) has no finite
    // solution where it evaluates to nan: the reference `solve_linear`
    // and `checksol` reject the candidate `x = 0` (`1 - oo*0 = nan`).
    if !product && has_infinity(&f) {
        sols.retain(|s| verify(&f, x, s));
    }
    if check && product {
        // Each factor's solutions were checked against that factor. Like the
        // reference, the union is then only filtered by the denominators of
        // every subexpression (`_simple_dens`), not checked against the whole
        // product: `x*log(x)` keeps `0`, where `log` is undefined.
        let dens = simple_dens(&f, x);
        sols.retain(|s| !dens.iter().any(|d| denominator_vanishes(d, x, s)) && matches_assumptions(x, s));
    } else if check {
        // drop solutions where the original expression is undefined
        let (_, den) = numer_denom(&together(&f));
        sols.retain(|s| {
            let d = expand(&subs1(&den, x, s));
            !(d.is_zero() || algebraic::polynomial_zero(&den, x, s) == Some(true)
                || matches!(d.kind(), Kind::Const(Constant::ComplexInfinity | Constant::NaN)))
                && matches_assumptions(x, s)
        });
        // Numerical substitution only verifies the original equation; it is
        // never used to infer a candidate's symbolic assumption properties.
        // The reference `_solve` sets `check = False` for the roots of a
        // numerator polynomial in the unknown itself (only the denominators
        // are checked): those roots are exact, and evaluating the residual of
        // large radical roots chases an exact zero through every precision.
        let (numerator, _) = numer_denom(&together(&f));
        if poly_coeffs(&numerator, x).is_none() {
            sols.retain(|s| verify(&f, x, s) && !checksol::rejects_parametric(&numerator, x, s));
        }
    }
    sort_solutions(&mut sols);
    Some(sols)
}

/// The reference `_simple_dens(f, symbols)`: the denominator factors of
/// every subexpression of `f` that depend on `x`, with a power of numeric
/// exponent replaced by its base (and `foo**0` dropped).
fn simple_dens(f: &Expr, x: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(f, &mut |s| {
        if matches!(s.kind(), Kind::Rel(..) | Kind::Bool(_) | Kind::Func(Func::Logic(_), _)) {
            return;
        }
        let (_, d) = crate::simptools::radsimp::fraction(s, false);
        if d.is_one() {
            return;
        }
        for factor in crate::simptools::mul_args(&d) {
            if !has(&factor, x) {
                continue;
            }
            let factor = match factor.kind() {
                Kind::Pow(base, exponent) if exponent.is_number() => {
                    if exponent.is_zero() { continue; }
                    base.clone()
                }
                _ => factor,
            };
            if !out.contains(&factor) {
                out.push(factor);
            }
        }
    });
    out
}

/// `checksol(den, {x: s})` is true: the denominator provably vanishes at
/// the candidate. An undecided value keeps the candidate, as in the
/// reference.
fn denominator_vanishes(den: &Expr, x: &Expr, s: &Expr) -> bool {
    let d = expand(&subs1(den, x, s));
    if d.is_zero() || algebraic::polynomial_zero(den, x, s) == Some(true) {
        return true;
    }
    if has(&d, x) || !crate::functions::is_numeric_value(&d) {
        return false;
    }
    if let Some(zero) = crate::assumptions::exact_queries(|| crate::assumptions::is_zero(&d)) {
        return zero;
    }
    crate::evalf_mp::eval_mp(&d, 64).is_some_and(|v| v.abs(64).cmp(&crate::bigfloat::BigFloat::one().ldexp(-30)).is_lt())
}

/// Compare a logarithm's algebraic argument with one without evaluating
/// the logarithm. For log(Abs(z)), compare |z|^2 instead: both comparisons
/// have the same sign. A rectangle containing zero proves nothing.
fn logarithm_argument_shift(argument: &Expr) -> Option<crate::poly::enclosure::CBox> {
    use crate::poly::enclosure::{enclose, CBox};
    if !free_symbols(argument).is_empty() {
        return None;
    }
    if let Kind::Func(Func::Abs, args) = argument.kind() {
        let modulus_squared = enclose(&args[0], 64).ok()?.norm2(64);
        if !modulus_squared.positive() { return None; }
        Some(CBox::real(modulus_squared).sub(&CBox::one(), 64))
    } else {
        if !enclose(argument, 64).ok()?.excludes_zero() { return None; }
        enclose(&sub(argument, &Expr::one()), 64).ok()
    }
}

/// Additional exact zero decisions needed by candidate-assumption checks.
/// In particular Im(-I*log(z)) = -log(Abs(z)) is nonzero whenever a
/// certified algebraic enclosure separates |z|^2 from one. Products retain
/// an undecided factor instead of declaring it nonzero from a sample.
fn certified_candidate_zero(value: &Expr) -> Option<bool> {
    use crate::poly::enclosure::enclose;
    let rectangle = match value.kind() {
        Kind::Func(Func::Log, args) if args.len() == 1 => logarithm_argument_shift(&args[0])?,
        Kind::Mul(coefficient, factors) if !coefficient.is_zero() => {
            for (base, exponent) in factors {
                // A nonzero finite number to a rational power is nonzero.
                // This only supplies a negative zero verdict, so a zero base
                // or an unsupported exponent remains undecided here.
                let ordinary = crate::assumptions::is_zero(base) == Some(false)
                    && crate::assumptions::is_finite(base) == Some(true);
                if !exponent.is_rational_number()
                    || !(ordinary || certified_candidate_zero(base) == Some(false)) {
                    return None;
                }
            }
            return Some(false);
        }
        _ if free_symbols(value).is_empty() => enclose(value, 64).ok()?,
        _ => return None,
    };
    if rectangle.is_exact_zero() { Some(true) }
    else if rectangle.excludes_zero() { Some(false) }
    else { None }
}

/// The sign of a real logarithm follows the argument's position relative
/// to one. This extends the algebraic sign check without approximating a
/// transcendental expression or rejecting symbolic parameter signs.
fn certified_candidate_sign(value: &Expr) -> Option<i32> {
    use crate::poly::enclosure::enclose;
    let interval = match value.kind() {
        Kind::Func(Func::Log, args) if args.len() == 1 => {
            let positive = if let Kind::Func(Func::Abs, inner) = args[0].kind() {
                enclose(&inner[0], 64).ok()?.norm2(64).positive()
            } else {
                let argument = enclose(&args[0], 64).ok()?;
                argument.im.is_exact_zero() && argument.re.positive()
            };
            if !positive { return None; }
            let shifted = logarithm_argument_shift(&args[0])?;
            if !shifted.im.is_exact_zero() { return None; }
            shifted.re
        }
        Kind::Mul(coefficient, factors) if coefficient.is_rational() => {
            let mut sign = if coefficient.is_positive() { 1 } else if coefficient.is_negative() { -1 } else { 0 };
            for (base, exponent) in factors {
                if !exponent.is_one() { return None; }
                sign *= certified_candidate_sign(base)?;
            }
            return Some(sign);
        }
        _ if free_symbols(value).is_empty() => {
            let rectangle = enclose(value, 64).ok()?;
            if !rectangle.im.is_exact_zero() { return None; }
            rectangle.re
        }
        _ => return None,
    };
    if interval.positive() { Some(1) }
    else if interval.negative() { Some(-1) }
    else if interval.is_exact_zero() { Some(0) }
    else { None }
}

fn matches_assumptions(variable: &Expr, value: &Expr) -> bool {
    use crate::assumptions as query;
    let Some(symbol) = variable.as_sym() else { return true; };
    let assumptions = symbol.facts();
    if assumptions.len()==1 { return true; } // default commutativity only
    query::exact_queries(|| {
        let algebraic=query::quadratic_facts(value);
        let enclosure=std::cell::OnceCell::new();
        let rectangle=||enclosure.get_or_init(|| {
            if free_symbols(value).is_empty() { crate::poly::enclosure::enclose(value,64).ok() }
            else { None }
        }).as_ref();
        let integer = if value.as_number().is_some_and(Number::is_float) { None } else { query::is_integer(value).or_else(||algebraic.and_then(|f|f.integer)) };
        let rational = if value.as_number().is_some_and(Number::is_float) { None } else { query::is_rational(value).or_else(||algebraic.and_then(|f|f.rational)) };
        let real = query::is_real(value).or_else(||algebraic.and_then(|f|f.real)).or_else(|| {
            // A compact nested principal radical need not have an explicit
            // Cartesian form. Decide its reality only from an outward-rounded
            // algebraic enclosure: an imaginary interval excluding zero is a
            // proof of non-reality, whereas containing zero is inconclusive.
            // This keeps real-domain filtering without expanding atan2/arg
            // expressions or rejecting tiny imaginary parts in machine precision.
            let enclosure = rectangle()?;
            if enclosure.im.is_exact_zero() { Some(true) }
            else if !enclosure.im.contains_zero() { Some(false) }
            else { None }
        }).or_else(|| {
            let imaginary_part = im(value);
            if matches!(imaginary_part.kind(), Kind::Func(Func::Im, _)) { None }
            else { query::is_zero(&imaginary_part).or_else(|| certified_candidate_zero(&imaginary_part)) }
        });
        let zero = query::is_zero(value).or_else(||algebraic.and_then(|f|f.zero));
        let nonzero=match real {Some(false)=>Some(false),Some(true)=>zero.map(|z|!z),None=>zero.filter(|z|*z).map(|_|false)};
        let imaginary = algebraic.and_then(|f|f.imaginary).or_else(||if zero == Some(true) || real == Some(true) { Some(false) }
            else {
                let real_part = re(value);
                if matches!(real_part.kind(), Kind::Func(Func::Re, _)) { None }
                else { query::is_zero(&real_part).and_then(|is_zero| if !is_zero { Some(false) } else { real.map(|real| !real) }) }
            });
        let prime = query::fact(value,"prime");
        let finite = query::is_finite(value);
        let complex = finite.filter(|finite| !finite).or_else(|| {
            if real == Some(true) || imaginary == Some(true) { Some(true) } else { None }
        });
        let irrational = match (real, rational) {
            (Some(false), _) | (_, Some(true)) => Some(false),
            (Some(true), Some(false)) => Some(true),
            _ => None,
        };
        let certified_sign=|| {
            if real != Some(true) { return None; }
            let interval=&rectangle()?.re;
            if interval.positive() { Some(1) }
            else if interval.negative() { Some(-1) }
            else if interval.is_exact_zero() { Some(0) }
            else { None }
        };
        let positive=query::is_positive(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s>0))
            .or_else(||certified_sign().or_else(||certified_candidate_sign(value)).map(|s|s>0));
        let negative=query::is_negative(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s<0))
            .or_else(||certified_sign().or_else(||certified_candidate_sign(value)).map(|s|s<0));
        let nonnegative=query::is_nonnegative(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s>=0))
            .or_else(||certified_sign().or_else(||certified_candidate_sign(value)).map(|s|s>=0));
        let nonpositive=query::is_nonpositive(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s<=0))
            .or_else(||certified_sign().or_else(||certified_candidate_sign(value)).map(|s|s<=0));
        assumptions.iter().all(|(name,expected)|{
            let known=match name {
                "real"=>real,"integer"=>integer,"rational"=>rational,"irrational"=>irrational,
                "zero"=>zero,"nonzero"=>nonzero,"positive"=>positive,"negative"=>negative,
                "nonnegative"=>nonnegative,"nonpositive"=>nonpositive,"prime"=>prime,
                "finite"=>finite,"complex"=>complex,"imaginary"=>imaginary,
                _=>query::fact(value,name),
            };
            known.is_none_or(|known|known==expected)
        })
    })
}

fn verify(f: &Expr, x: &Expr, s: &Expr) -> bool {
    if let Some(zero) = algebraic::polynomial_zero(f, x, s) {
        return zero;
    }
    verify_equation(f, &[(x.clone(), s.clone())])
}

/// Check exact identities first. A numerical rejection of an exact candidate
/// requires a stable residual at two increasing multiprecision evaluations,
/// well above equation-scaled rounding error. Undecidable cases are retained.
/// Approximate data uses its declared precision and an equation-relative scale.
fn verify_equation(f: &Expr, assignments: &[(Expr, Expr)]) -> bool {
    let v = subs(f, assignments);
    if v.is_zero() { return true; }
    if matches!(v.kind(), Kind::Const(Constant::ComplexInfinity | Constant::NaN)) { return false; }
    // the reference `checksol` rejects a value holding nan, zoo or +-oo
    if tsolve::has_illegal(&v) { return false; }
    let precision = assignments.iter().filter_map(|(_, value)| crate::evalf_mp::max_precision(value))
        .chain(crate::evalf_mp::max_precision(f)).max();
    let Some(precision) = precision else {
        if let Some(zero) = crate::assumptions::exact_queries(|| crate::assumptions::is_zero(&v)) {
            return zero;
        }
        if !crate::functions::is_numeric_value(&v) { return true; }
        let mut previous = None;
        for bits in [160, 320] {
            let scale = residual_scale(f, assignments, bits);
            // An exactly cancelling component (the real part of
            // `log(x) + log(x - 1)` at `1/2 - sqrt(5)/2`, whose residual is
            // `2*I*pi`) is resolved by chopping it below the tolerance, as
            // the reference `checksol` evaluates `val.n(18).n(12, chop=True)`.
            let cutoff = scale.mag().saturating_add(24 - bits as i128);
            // A sum cancelling more than `bits` bits of its largest term is
            // already within the tolerance below, so a larger cancellation
            // budget only chases an exact zero (a residual at exact radical
            // roots) through thousands of bits; the chopped evaluation then
            // decides it.
            let Some(value) = crate::evalf_mp::eval_mp_cancelling(&v, bits, bits)
                .or_else(|| crate::evalf_mp::eval_mp_chopped(&v, bits, cutoff)) else { return true; };
            let magnitude = value.abs(bits);
            if magnitude.cmp(&scale.ldexp(24 - bits as i64)).is_le() { return true; }
            if let Some(prior) = previous {
                let change = value.sub(&prior, bits).abs(bits);
                if change.cmp(&magnitude.ldexp(-24)).is_le() { return false; }
            }
            previous = Some(value);
        }
        return true;
    };
    let working = precision + 64;
    let Some(residual) = crate::evalf_mp::eval_mp(&v, working) else { return true; };
    let scale = residual_scale(f, assignments, working);
    residual.abs(working).cmp(&scale.ldexp(8 - precision as i64)).is_le()
}

fn residual_scale(f: &Expr, assignments: &[(Expr, Expr)], working: u32) -> crate::bigfloat::BigFloat {
    let mut scale = crate::bigfloat::BigFloat::one();
    for term in if f.is_add() { f.args() } else { vec![f.clone()] } {
        if let Some(value) = crate::evalf_mp::eval_mp(&subs(&term, assignments), working) {
            let magnitude = value.abs(working);
            if magnitude.cmp(&scale).is_gt() { scale = magnitude; }
        }
    }
    scale
}

fn solve_rec(f: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    if depth > 8 {
        return None;
    }
    if !has(f, x) {
        return if f.is_zero() { None } else { Some(vec![]) };
    }
    let (num, _) = numer_denom(&together(f));
    let num = if has(&num, x) { num } else { return Some(vec![]) };
    // numer_denom multiplies distinct denominators (log(a)/6 - 1/2 gives
    // 2*log(a) - 6 over 12); SymPy's together uses their lcm. A spurious
    // integer content would let logcombine raise the arguments to a power
    // and add roots of other branches, so it is removed.
    let num = without_integer_content(&num);
    // polynomial in x
    if let Some(c) = poly_coeffs(&num, x) {
        // The reference solves polynomials through `roots` over the composite
        // coefficient domain, which factors symbolic coefficients first, and
        // through `Poly.all_roots` (root objects) when `roots` is incomplete.
        let degree = c.iter().rposition(|t| !t.is_zero()).unwrap_or(0);
        let exact = crate::evalf_mp::max_precision(&num).is_none();
        if exact && degree >= 3 {
            if let Some(r) = complete_roots(&num, x, &c).or_else(|| algebraic::all_roots(&num, x)) {
                return Some(r);
            }
        } else if c.iter().any(|t| !t.is_rational_number()) && (c.len() > 3 || splits_quadratic(&num, x)) {
            if let Some(r) = complete_roots(&num, x, &c) {
                return Some(r);
            }
        }
        // An incomplete `roots` result with symbolic coefficients is
        // reported by the reference as a partial solution list; the closed
        // forms below are tried before the solver declares it unsupported.
        return poly_roots(&c).or_else(|| complete_roots(&num, x, &c))
            .or_else(|| algebraic::all_roots(&num, x))
            .or_else(|| partial_roots(&num, x));
    }
    // product = 0 -> union of the factors' solutions
    if let Kind::Mul(_, facs) = num.kind() {
        let mut out = Vec::new();
        for (b, e) in facs {
            if !has(b, x) {
                continue;
            }
            if e.as_number().is_some_and(|n| n.is_negative()) {
                continue;
            }
            // exp(u) never vanishes. SymPy's as_numer_denom moves exp(-u)
            // to the denominator, which the numerator here keeps.
            if matches!(b.kind(), Kind::Func(Func::Exp, _)) {
                continue;
            }
            out.extend(solve_rec(b, x, depth + 1)?);
        }
        return Some(out);
    }
    // several trig generators and no other function: rewrite to tan
    {
        let mut gens = Vec::new();
        generators(&num, x, &mut gens);
        if let Some(newf) = trig::tan_rewrite(&num, x, &gens) {
            return solve_rec(&newf, x, depth + 1).map(|r| r.iter().map(tabulated_atan).collect());
        }
    }
    if let Some(r) = symbolic_power::solve(&num, x, depth) {
        return Some(r);
    }
    // polynomial in a single generator g(x)
    if let Some(r) = solve_by_generator(&num, x, depth) {
        return Some(r);
    }
    // isolate the only term depending on x and invert
    if let Some(r) = invert(&num, x, depth) {
        return Some(r);
    }
    let flags = checksol::solve_flags();
    // The reference `_tsolve` skips its Lambert step for `bivariate=False`;
    // its closing posify retry (only with `force`) restores it.
    if flags.bivariate || flags.force {
        if let Some(r) = lambert::solve(&num, x) {
            return Some(r.iter().map(lambert::simpler_rational).collect());
        }
    }
    // several logarithms: combine them (log(a) - log(b) = log(a/b)), with
    // the reference `logcombine(lhs, force=flags.get('force', True))`
    let lc = crate::simp::logcombine(&num, flags.force);
    if lc != num {
        if let Some(r) = solve_rec(&lc, x, depth + 1) {
            return Some(r);
        }
    }
    let ex = expand(&num);
    if ex != num {
        return solve_rec(&ex, x, depth + 1);
    }
    // SymPy's as_numer_denom moves exp(-u) to the denominator, so its
    // numerator of C1 + (x + y)*exp(-y) is C1*exp(y) + x + y.
    if let Some(cleared) = exp_denominators_cleared(&num, x) {
        return solve_rec(&cleared, x, depth + 1);
    }
    if let Some(combined) = combined_powers(&num, x) {
        return solve_rec(&combined, x, depth + 1);
    }
    // the `rewrite(exp)` fallback of the reference `_tsolve`
    if let Some(rewritten) = trig::exp_rewrite(&num, x) {
        return solve_rec(&rewritten, x, depth + 1);
    }
    None
}

fn contains_node(e: &Expr, p: impl Fn(&Expr) -> bool) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found |= p(s));
    found
}

fn has_free_parameter(e: &Expr) -> bool {
    !free_symbols(e).is_empty()
}

/// `2*atan(r)` from the tan half-angle rewrite: the reference `simplify`
/// rationalizes `r` (`1/(1 + sqrt(2))` to `-1 + sqrt(2)`), which evaluates
/// `atan` at the tabulated angles. Only that evaluation is adopted.
fn tabulated_atan(e: &Expr) -> Expr {
    let mut map = FxHashMap::default();
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(Func::Atan, a) = s.kind() {
            let v = atan(&crate::simplify::radsimp(&a[0]));
            if !contains_node(&v, |t| matches!(t.kind(), Kind::Func(Func::Atan, _))) {
                map.insert(s.clone(), v);
            }
        }
    });
    if map.is_empty() { e.clone() } else { crate::subs::xreplace(e, &map) }
}

fn has_infinity(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN) {
            found = true;
        }
    });
    found
}

/// `exp(u)` or `b**u` with `u` depending on the unknown (the reference
/// `_ispow` on a dependent part).
fn is_dependent_power(e: &Expr, x: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Exp, _) => true,
        Kind::Pow(_, exponent) => has(exponent, x),
        _ => false,
    }
}

/// The coefficient and the dependent part of a term (`as_independent`).
fn split_term(t: &Expr, x: &Expr) -> (Expr, Expr) {
    match t.kind() {
        Kind::Mul(..) => {
            let (cs, ds): (Vec<Expr>, Vec<Expr>) = t.args().into_iter().partition(|a| !has(a, x));
            (mul(&cs), mul(&ds))
        }
        _ => (Expr::one(), t.clone()),
    }
}

/// Powers of different bases combined like the reference `_invert`.
///
/// A sum of two dependent terms `ai*ad + bi*bd` (no independent part) whose
/// parts are powers of different bases is written `powsimp(ad/bd) =
/// -bi/ai` (`3**x - 2**x` gives `(3/2)**x = 1`, so `[0]`), and a single
/// dependent product of powers is combined by `powsimp(powdenest(...))`
/// (`2**x*3**x - 6` gives `6**x = 6`, so `[1]`). `None` when the combination
/// leaves the parts unchanged.
fn combined_powers(f: &Expr, x: &Expr) -> Option<Expr> {
    let terms = if f.is_add() { f.args() } else { vec![f.clone()] };
    let (dependent, independent): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(|t| has(t, x));
    let options = crate::simptools::powsimp::PowsimpOptions::default();
    match dependent.as_slice() {
        [first, second] if independent.is_empty() => {
            let mut pair = vec![first.clone(), second.clone()];
            crate::logic::ordered(&mut pair);
            let (ai, ad) = split_term(&pair[0], x);
            let (bi, bd) = split_term(&pair[1], x);
            if !(is_dependent_power(&ad, x) || is_dependent_power(&bd, x)) {
                return None;
            }
            if crate::ops::as_base_exp(&ad).0 == crate::ops::as_base_exp(&bd).0 {
                return None;
            }
            let ratio = div(&ad, &bd);
            let lhs = crate::simptools::powsimp::powsimp(&ratio, &options);
            if lhs == ratio {
                return None;
            }
            Some(sub(&lhs, &neg(&div(&bi, &ai))))
        }
        [term] => {
            let (coefficient, part) = split_term(term, x);
            if !matches!(part.kind(), Kind::Mul(..)) || !part.args().iter().any(|a| is_dependent_power(a, x)) {
                return None;
            }
            let lhs = crate::simptools::powsimp::powsimp(&crate::simptools::powsimp::powdenest(&part, false, false), &options);
            if lhs == part {
                return None;
            }
            Some(sub(&lhs, &div(&neg(&add(&independent)), &coefficient)))
        }
        _ => None,
    }
}

/// Whether a quadratic factors in `x` over its coefficient domain; an
/// irreducible quadratic keeps the closed form of `poly_roots`, which is the
/// simplified reference result.
fn splits_quadratic(num: &Expr, x: &Expr) -> bool {
    crate::poly::api::factor_list(num).is_ok_and(|(_, fs)| {
        fs.iter()
            .filter(|(g, _)| has(g, x))
            .map(|(_, k)| k.as_i64().unwrap_or(1))
            .sum::<i64>()
            > 1
    })
}

/// The sum multiplied by the product of its distinct exp denominators
/// (exp(u) for each factor exp(-u) depending on x), or None.
fn exp_denominators_cleared(e: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Add(..) = e.kind() else { return None; };
    let terms = e.args();
    let mut denominators: Vec<Expr> = Vec::new();
    for term in &terms {
        let factors = match term.kind() {
            Kind::Mul(_, facs) => facs.iter().filter(|(_, k)| k.is_one()).map(|(b, _)| b.clone()).collect(),
            _ => vec![term.clone()],
        };
        for factor in factors {
            if let Kind::Func(Func::Exp, args) = factor.kind() {
                if has(&args[0], x) && crate::decompose::could_extract_minus_sign(&args[0]) {
                    let d = neg(&args[0]);
                    if !denominators.contains(&d) {
                        denominators.push(d);
                    }
                }
            }
        }
    }
    if denominators.is_empty() {
        return None;
    }
    let multiplier = exp(&add(&denominators));
    let cleared = add(&terms.iter().map(|t| crate::simp::powsimp(&mul2(t, &multiplier))).collect::<Vec<_>>());
    (cleared != *e).then_some(cleared)
}

/// A sum divided by the gcd of its integer term coefficients.
fn without_integer_content(e: &Expr) -> Expr {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let Kind::Add(..) = e.kind() else { return e.clone(); };
    let terms = e.args();
    let mut g = BigInt::from(0);
    for t in &terms {
        let c = match t.kind() {
            Kind::Num(n) => n.clone(),
            Kind::Mul(c, _) => c.clone(),
            _ => return e.clone(),
        };
        match c.to_bigint() {
            Some(n) if c.is_integer() => g = g.gcd(&n),
            _ => return e.clone(),
        }
    }
    if g <= BigInt::from(1) {
        return e.clone();
    }
    let Some(inverse) = Number::from_bigint(g).inv() else { return e.clone(); };
    let inverse = Expr::num(inverse);
    add(&terms.iter().map(|t| mul2(&inverse, t)).collect::<Vec<_>>())
}

/// The full `roots` machinery (quartics, decomposition such as x**8 in
/// x**2, factorisation), accepted only when it finds every root: SymPy's
/// `roots` may return an incomplete set, which `solve` must not report.
fn complete_roots(num: &Expr, x: &Expr, c: &[Expr]) -> Option<Vec<Expr>> {
    let deg = c.iter().rposition(|t| !t.is_zero())? as u32;
    let rs = crate::poly::roots::roots(num, std::slice::from_ref(x)).ok()?;
    (rs.iter().map(|r| r.1).sum::<u32>() == deg).then(|| rs.into_iter().map(|r| r.0).collect())
}

/// The reference `_solve` with `incomplete=True` (its default): when `roots`
/// misses some roots and `Poly.all_roots` raises `NotImplementedError`
/// (symbolic coefficients), the roots that were found are the solutions,
/// e.g. none for `x**5 + a*x + 1`. This is a partial answer, not a proof
/// that the other roots do not exist.
fn partial_roots(num: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    if !crate::poly::roots::radical_options().incomplete || crate::evalf_mp::max_precision(num).is_some() {
        return None;
    }
    // `all_roots` supports the rational and Gaussian rational domains (and
    // floats); only the other coefficient domains give a partial answer.
    let gaussian = |c: &Expr| {
        let c = expand(c);
        c.is_rational_number() || (re(&c).is_rational_number() && im(&c).is_rational_number())
    };
    if poly_coeffs(num, x)?.iter().all(gaussian) {
        return None;
    }
    let found = crate::poly::roots::roots(num, std::slice::from_ref(x)).ok()?;
    Some(found.into_iter().map(|r| r.0).collect())
}

/// Roots of a polynomial given by its dense coefficients (low degree
/// first). Handles linear, quadratic, binomial and rational roots with
/// deflation; delegates nothing else (returns `None` if a factor of degree
/// >= 3 remains).
pub fn poly_roots(c: &[Expr]) -> Option<Vec<Expr>> {
    let mut c: Vec<Expr> = c.to_vec();
    while c.len() > 1 && c.last().unwrap().is_zero() {
        c.pop();
    }
    let mut out = Vec::new();
    // x = 0 roots
    let mut k = 0;
    while k < c.len() - 1 && c[k].is_zero() {
        k += 1;
    }
    if k > 0 {
        out.push(Expr::zero());
        c.drain(..k);
    }
    let deg = c.len() - 1;
    match deg {
        0 => return Some(out),
        1 => {
            out.push(neg(&div(&c[0], &c[1])));
            return Some(out.into_iter().map(|e| expand(&e)).collect());
        }
        2 => {
            out.extend(quadratic(&c[2], &c[1], &c[0]));
            return Some(out);
        }
        _ => {}
    }
    // binomial a x^n + b
    if c[1..deg].iter().all(|t| t.is_zero()) {
        out.extend(binomial_roots(&c[deg], &c[0], deg as u32));
        return Some(out);
    }
    // rational roots for rational coefficients
    if c.iter().all(|t| t.is_rational_number()) {
        let mut rem = c.clone();
        let mut found = true;
        while found && rem.len() > 3 {
            found = false;
            for r in rational_root_candidates(&rem) {
                if let Some(q) = deflate(&rem, &r) {
                    out.push(Expr::num(r));
                    rem = q;
                    found = true;
                    break;
                }
            }
        }
        if rem.len() - 1 <= 2 {
            out.extend(poly_roots(&rem)?);
            return Some(out);
        }
        // biquadratic a x^4 + b x^2 + c
        if rem.len() == 5 && rem[1].is_zero() && rem[3].is_zero() {
            for r in quadratic(&rem[4], &rem[2], &rem[0]) {
                let s = sqrt(&r);
                out.push(neg(&s));
                out.push(s);
            }
            return Some(out);
        }
        // closed forms for cubics
        if rem.len() == 4 && crate::poly::roots::radical_options().cubics {
            out.extend(cubic(&rem)?);
            return Some(out);
        }
        return None;
    }
    if deg == 3 && crate::poly::roots::radical_options().cubics {
        out.extend(cubic(&c)?);
        return Some(out);
    }
    None
}

fn quadratic(a: &Expr, b: &Expr, c: &Expr) -> Vec<Expr> {
    if c.is_zero() {
        return vec![Expr::zero(), neg(&div(b, a))];
    }
    if b.is_zero() {
        let r = neg(&div(c, a));
        let s = sqrt_squares_out(&composite_simplify(r));
        return vec![neg(&s), s];
    }
    let d = expand(&sub(&b.powi(2), &mul(&[Expr::int(4), a.clone(), c.clone()])));
    let two_a = mul2(&Expr::int(2), a);
    if crate::functions::is_numeric_value(&d) && crate::functions::is_numeric_value(a) && crate::functions::is_numeric_value(b) {
        let bb = neg(&div(b, &two_a));
        let dd = div(&sqrt(&d), &two_a);
        if dd.is_zero() {
            return vec![bb];
        }
        return vec![sub(&bb, &dd), add2(&bb, &dd)];
    }
    let s = sqrt_squares_out(&composite_simplify(d));
    vec![div(&sub(&neg(b), &s), &two_a), div(&add2(&neg(b), &s), &two_a)]
}

/// The reference `roots_quadratic._simplify` for a composite coefficient
/// domain (`ZZ[a]`, `QQ(a)`): the radicand is factored before its squares
/// are taken out. Numerical and `EX` radicands are kept as they are.
fn composite_simplify(e: Expr) -> Expr {
    if crate::functions::is_numeric_value(&e) || !is_rational_function(&e) {
        return e;
    }
    crate::poly::api::factor(&e).unwrap_or(e)
}

/// Whether `e` is built from rationals and symbols with integer powers.
fn is_rational_function(e: &Expr) -> bool {
    let mut ok = true;
    crate::traverse::preorder(e, &mut |x| {
        ok &= match x.kind() {
            Kind::Num(n) => n.is_rational(),
            Kind::Pow(_, k) => k.is_integer_number(),
            Kind::Add(..) | Kind::Mul(..) => true,
            _ => x.is_symbol(),
        };
    });
    ok
}

/// The reference `roots_quadratic._sqrt`: even powers leave the square root,
/// both signs being among the roots (`sqrt(-b**2)` gives `I*b`).
fn sqrt_squares_out(d: &Expr) -> Expr {
    let factors = match d.kind() {
        Kind::Mul(..) => d.args(),
        _ => vec![d.clone()],
    };
    let (mut co, mut other) = (Vec::new(), Vec::new());
    for t in factors {
        match t.kind() {
            Kind::Pow(base, e) if e.as_i64().is_some_and(|n| n % 2 == 0) => {
                co.push(pow(base, &div(e, &Expr::int(2))))
            }
            _ => other.push(t),
        }
    }
    if co.is_empty() {
        return sqrt(d);
    }
    mul2(&mul(&co), &sqrt(&mul(&other)))
}

/// Roots of a x^n + b.
fn binomial_roots(a: &Expr, b: &Expr, n: u32) -> Vec<Expr> {
    let r = neg(&div(b, a));
    let base = pow(&r, &Expr::rational(1, n as i64));
    let mut out = Vec::with_capacity(n as usize);
    for k in 0..n {
        // exp(2 pi I k / n)
        let ang = mul(&[Expr::int(2), Expr::pi(), Expr::rational(k as i64, n as i64)]);
        let w = add2(&cos(&ang), &mul2(&Expr::i(), &sin(&ang)));
        out.push(expand(&mul2(&base, &w)));
    }
    out
}

/// Candidate rational roots p/q (p | a0, q | an).
fn rational_root_candidates(c: &[Expr]) -> Vec<Number> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    // make integer coefficients
    let mut den = BigInt::from(1);
    for t in c {
        if let Some((_, d)) = t.as_number().and_then(|n| n.numer_denom()) {
            den = den.lcm(&d);
        }
    }
    let ints: Vec<BigInt> = c
        .iter()
        .map(|t| {
            let (n, d) = t.as_number().unwrap().numer_denom().unwrap();
            n * (&den / d)
        })
        .collect();
    let a0 = ints[0].clone();
    let an = ints.last().unwrap().clone();
    if a0 == BigInt::from(0) {
        return vec![Number::ZERO];
    }
    let small = |n: &BigInt| num_traits::ToPrimitive::to_u64(&num_traits::Signed::abs(n)).is_some_and(|v| v < 1u64 << 40);
    if !small(&a0) || !small(&an) {
        return vec![];
    }
    let ps = crate::ntheory::divisors(&a0);
    let qs = crate::ntheory::divisors(&an);
    let mut out = Vec::new();
    for p in &ps {
        for q in &qs {
            for s in [1, -1] {
                if let Some(r) = Number::rational(p * BigInt::from(s), q.clone()) {
                    if !out.contains(&r) {
                        out.push(r);
                    }
                }
            }
        }
    }
    out.sort_by(|a, b| a.abs().cmp_value(&b.abs()).unwrap().then(a.cmp_value(b).unwrap()));
    out
}

/// Divide the polynomial by (x - r) if r is a root.
fn deflate(c: &[Expr], r: &Number) -> Option<Vec<Expr>> {
    let n = c.len() - 1;
    let mut q = vec![Number::ZERO; n];
    let mut acc = Number::ZERO;
    for i in (0..=n).rev() {
        let ci = c[i].as_number()?;
        acc = acc.mul(r).add(ci);
        if i > 0 {
            q[i - 1] = acc.clone();
        }
    }
    if !acc.is_zero() {
        return None;
    }
    Some(q.into_iter().map(Expr::num).collect())
}

/// Cardano's formula for a cubic (coefficients low degree first).
fn cubic(c: &[Expr]) -> Option<Vec<Expr>> {
    let (d, cc, b, a) = (&c[0], &c[1], &c[2], &c[3]);
    // depressed cubic t^3 + p t + q with x = t - b/(3a)
    let p = expand(&div(&sub(&mul(&[Expr::int(3), a.clone(), cc.clone()]), &b.powi(2)), &mul(&[Expr::int(3), a.powi(2)])));
    let q = expand(&div(
        &add(&[
            mul(&[Expr::int(2), b.powi(3)]),
            mul(&[Expr::int(-9), a.clone(), b.clone(), cc.clone()]),
            mul(&[Expr::int(27), a.powi(2), d.clone()]),
        ]),
        &mul(&[Expr::int(27), a.powi(3)]),
    ));
    let shift = neg(&div(b, &mul2(&Expr::int(3), a)));
    if p.is_zero() {
        return Some(binomial_roots(&Expr::one(), &q, 3).into_iter().map(|t| expand(&add2(&t, &shift))).collect());
    }
    // C = cbrt(-q/2 + sqrt(q^2/4 + p^3/27))
    let disc = expand(&add2(&div(&q.powi(2), &Expr::int(4)), &div(&p.powi(3), &Expr::int(27))));
    let cc_ = pow(&add2(&neg(&div(&q, &Expr::int(2))), &sqrt(&disc)), &Expr::rational(1, 3));
    let w = add2(&Expr::rational(-1, 2), &mul(&[Expr::i(), sqrt(&Expr::int(3)), Expr::rational(1, 2)]));
    let mut out = Vec::new();
    for k in 0..3 {
        let ck = mul2(&cc_, &pow(&w, &Expr::int(k)));
        let t = sub(&ck, &div(&p, &mul2(&Expr::int(3), &ck)));
        out.push(add2(&t, &shift));
    }
    Some(out)
}

/// Generators: maximal sub-expressions depending on x that are not
/// polynomial in x.
fn generators(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    if !has(e, x) || e == x {
        return;
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) => {
            for a in e.args() {
                generators(&a, x, out);
            }
        }
        Kind::Pow(b, n) if n.is_integer_number() => generators(b, x, out),
        _ => {
            if !out.contains(e) {
                out.push(e.clone());
            }
        }
    }
}

fn solve_by_generator(f: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    let u = Expr::dummy("_u", crate::Assumptions::NONE);
    // exp(k*x + c): write as exp(c) * U**(k*q) with U = exp(x/q)
    let mut gens = Vec::new();
    generators(f, x, &mut gens);
    if gens.is_empty() {
        return None;
    }
    if let Some(r) = power_generator_tsolve(f, x, &gens, depth) {
        return Some(r);
    }
    let rewritten = perfect_power_bases(f, x, &gens);
    let f = rewritten.as_ref().unwrap_or(f);
    if rewritten.is_some() {
        gens.clear();
        generators(f, x, &mut gens);
    }
    // `exp(a)` or `b**a` with one base `b` free of x for every generator:
    // the reference's `decompose_power` takes them all as powers of one
    // generator (`2**(2*x) - 3*2**x + 2` is a quadratic in `2**x`).
    let exponent_of = |g: &Expr| -> Option<(Option<Expr>, Expr)> {
        match g.kind() {
            Kind::Func(Func::Exp, a) => Some((None, a[0].clone())),
            Kind::Pow(b, e) if !has(b, x) => Some((Some(b.clone()), e.clone())),
            _ => None,
        }
    };
    let decomposed: Option<Vec<(Option<Expr>, Expr)>> = gens.iter().map(exponent_of).collect();
    let common_base = decomposed.as_ref().filter(|_| gens.len() > 1).and_then(|d| {
        let base = d.first()?.0.clone();
        d.iter().all(|(b, _)| *b == base).then_some(base)
    });
    if let (Some(decomposed), Some(base)) = (decomposed.as_ref(), common_base) {
        let raise = |e: &Expr| match &base {
            Some(base) => pow(base, e),
            None => exp(e),
        };
        // rational multiples of x in the exponents; for `exp`, rational
        // multiples of `I` (`exp(I*x/2)`, `exp(3*I*x/2)`), the generator
        // being `exp(I*x/q)`
        let mut ks = Vec::new();
        let mut unit: Option<Expr> = None;
        let mut slopes = Vec::new();
        for (_, a) in decomposed {
            let k = diff(a, x);
            if has(&k, x) {
                return None;
            }
            let (r, rest) = if k.is_rational_number() {
                (k.clone(), Expr::one())
            } else if base.is_none() && free_symbols(&k).is_empty() {
                // only the `I*x` multiples of the `rewrite(exp)` of trig
                // functions; other constants keep the earlier paths
                let (r, rest) = crate::ops::as_coeff_mul(&k);
                if rest != Expr::i() {
                    return None;
                }
                (Expr::num(r), rest)
            } else {
                return None;
            };
            if !r.is_rational_number() || unit.as_ref().is_some_and(|u| *u != rest) {
                return None;
            }
            unit = Some(rest);
            slopes.push(k);
            ks.push(r);
        }
        let unit = unit.unwrap_or_else(Expr::one);
        let mut q = num_bigint::BigInt::from(1);
        for k in &ks {
            q = num_integer::Integer::lcm(&q, &k.as_number().unwrap().numer_denom().unwrap().1);
        }
        let qn = Number::from_bigint(q);
        let mut g2 = f.clone();
        for (((g, k), slope), (_, a)) in gens.iter().zip(&ks).zip(&slopes).zip(decomposed) {
            let c0 = expand(&sub(a, &mul2(slope, x)));
            let power = Expr::num(k.as_number().unwrap().mul(&qn));
            g2 = subs1(&g2, g, &mul2(&raise(&c0), &pow(&u, &power)));
        }
        if has(&g2, x) {
            return None;
        }
        let roots = solve_rec(&g2, &u, depth + 1)?;
        if base.is_some() {
            let roots = roots.into_iter().filter(|r| !r.is_zero()).collect();
            return invert_generator(&raise(&div(x, &Expr::num(qn))), roots, x, depth);
        }
        let mut out = Vec::new();
        for r in roots {
            if r.is_zero() {
                continue;
            }
            out.push(div(&mul2(&Expr::num(qn.clone()), &log(&r)), &unit));
        }
        // The reference simplifies the logarithms of the roots on the unit
        // circle along a path not reproduced here: only the lattices whose
        // logarithms all evaluate (`log(-1)`, `log(I)`) are adopted.
        if !unit.is_one() && out.iter().any(|v| contains_node(v, |t| matches!(t.kind(), Kind::Func(Func::Log, _)))) {
            return None;
        }
        return Some(out);
    }
    if gens.len() != 1 {
        return None;
    }
    let g = &gens[0];
    // The reference solves the numerator of `as_numer_denom`, which moves
    // `exp(a)` to the denominator as `1/exp(-a)` when `a` is negative or,
    // when `-a` is not known to be negative, when `a` could extract a minus
    // sign (`ExpBase.as_numer_denom`). Its generator is then `exp(-a)`, so
    // `exp(-x) + 1` gives `I*pi`, not `-I*pi`.
    if let Kind::Func(Func::Exp, arguments) = g.kind() {
        let a = &arguments[0];
        let negated = neg(a);
        let negative = crate::assumptions::is_negative(a) == Some(true)
            || (crate::assumptions::is_negative(&negated) != Some(true)
                && crate::decompose::could_extract_minus_sign(a));
        if negative {
            let generator = exp(&negated);
            // `1/exp(-a)` evaluates back to `exp(a)`: clear it through a Dummy.
            let numerator = numer_denom(&together(&subs1(f, g, &div(&Expr::one(), &u)))).0;
            let numerator = subs1(&numerator, &u, &generator);
            let mut moved = Vec::new();
            generators(&numerator, x, &mut moved);
            if moved == [generator] {
                return solve_rec(&numerator, x, depth + 1);
            }
        }
    }
    if let Some(r) = solve_exp_power(f, g, x, depth) {
        return Some(r);
    }
    let fu = subs1(f, g, &u);
    if has(&fu, x) {
        return None;
    }
    let Some(coefficients) = poly_coeffs(&fu, &u) else {
        return None;
    };
    // A generator that is not a power of a fixed base is inverted by solving
    // `g - r` again, which for a linear `f` is `f` itself: the inversion
    // steps below handle it directly instead of exhausting the depth.
    let exponential = match g.kind() {
        Kind::Func(Func::Exp, _) => true,
        Kind::Pow(base, _) => !has(base, x),
        _ => false,
    };
    if !exponential && coefficients.len() <= 2 {
        return None;
    }
    let roots = solve_rec(&fu, &u, depth + 1)?;
    if exponential && coefficients.len() == 2 && is_active_inversion(g, &roots) {
        return None;
    }
    invert_generator(g, roots, x, depth)
}

/// `b**g` with `b` free of the unknown and other than `E`.
fn is_power_generator(g: &Expr, x: &Expr) -> bool {
    matches!(g.kind(), Kind::Pow(b, e) if !has(b, x) && has(e, x) && *b != Expr::e())
}

/// The first steps of the reference `_tsolve` for an equation whose
/// generators are powers `b**g` of bases other than `E`.
///
/// The reference cannot invert such a generator for a symbol (`_tsolve`
/// raises for `2**x - t` with a Dummy `t`), so `_solve` discards the roots
/// of the polynomial in the generator and calls `_tsolve` on the whole
/// equation, which inverts the single dependent term first: `b**g = r` with
/// `r = b**e` (as `as_base_exp` gives it) is solved as `g = e`, one branch
/// only (`2**(2*x) - 2` gives `[1/2]`). With several dependent terms it
/// solves `factor(powdenest(f))` factor by factor when that is a product,
/// so `4**(2*x) - 3*4**x + 2`, which factors as
/// `(2**x - 1)*(2**x + 1)*(2**(2*x) - 2)`, has the three solutions
/// `0, I*pi/log(2), 1/2`. `None` when neither step applies; the native
/// polynomial and inversion steps then emulate the remaining `_tsolve`
/// steps (the rewrite through `exp`).
fn power_generator_tsolve(f: &Expr, x: &Expr, gens: &[Expr], depth: u32) -> Option<Vec<Expr>> {
    if !gens.iter().all(|g| is_power_generator(g, x)) {
        return None;
    }
    let terms = if f.is_add() { f.args() } else { vec![f.clone()] };
    let (dependent, independent): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(|t| has(t, x));
    if let [term] = dependent.as_slice() {
        let (coefficients, powers): (Vec<Expr>, Vec<Expr>) = match term.kind() {
            Kind::Mul(..) => term.args().into_iter().partition(|a| !has(a, x)),
            _ => (vec![], vec![term.clone()]),
        };
        let [power] = powers.as_slice() else { return None };
        if !is_power_generator(power, x) {
            return None;
        }
        // The base is taken as it is: `(1/4)**x - 4` has two solutions.
        let (base, exponent) = crate::ops::as_base_exp(power);
        let rhs = div(&neg(&add(&independent)), &mul(&coefficients));
        let (rhs_base, rhs_exponent) = crate::ops::as_base_exp(&rhs);
        if rhs_base == base {
            return solve_rec(&sub(&exponent, &rhs_exponent), x, depth + 1);
        }
        // For a negative base, `_tsolve` first guesses the rational exponent
        // `e_rat = nsimplify(log(|rhs|)/log(|b|))` and keeps it when
        // `b**e_rat == rhs`, then adds the principal logarithm solutions
        // (`(-2)**x - 4` gives `[2, log(4)/(log(2) + I*pi)]`).
        if let Some(k) = negative_base_ratio(&rhs, &base) {
            let mut out = solve_rec(&sub(&exponent, &k), x, depth + 1)?;
            for r in invert(f, x, depth)? {
                if !out.contains(&r) {
                    out.push(r);
                }
            }
            return Some(out);
        }
        return None;
    }
    let factored = crate::poly::api::factor(&crate::simptools::powsimp::powdenest(f, false, false)).ok()?;
    let Kind::Mul(..) = factored.kind() else { return None };
    let factors: Vec<Expr> = factored.args().into_iter().filter(|a| has(a, x)).collect();
    if factors.len() < 2 {
        return None;
    }
    let mut out = Vec::new();
    for factor in factors {
        let (b, e) = crate::ops::as_base_exp(&factor);
        let factor = if e.is_integer_number() && e.as_number().is_some_and(|n| n.is_positive()) { b } else { factor };
        // `b**g = 0` has no solution for a nonzero base free of the unknown
        // (the reference solves `b = 0`).
        if is_power_generator(&factor, x) && matches!(factor.kind(), Kind::Pow(b, _) if !b.is_zero()) {
            continue;
        }
        for s in solve_rec(&factor, x, depth + 1)? {
            if !out.contains(&s) {
                out.push(s);
            }
        }
    }
    Some(out)
}

/// `b = p**k` for a positive rational base `b`: an integer perfect power
/// (`4 = 2**2`) or the reciprocal of an integer (`1/4 = 2**(-2)`,
/// `1/3 = 3**(-1)`), with `p` the smallest integer root. `None` otherwise.
fn power_of_integer_base(b: &Expr) -> Option<(Expr, Expr)> {
    let (numerator, denominator) = b.as_number().filter(|n| n.is_rational())?.numer_denom()?;
    let one = num_bigint::BigInt::from(1);
    let (n, sign) = if denominator == one && numerator > one {
        (numerator, 1i64)
    } else if numerator == one && denominator > one {
        (denominator, -1i64)
    } else {
        return None;
    };
    let (p, k) = crate::ntheory::factor::perfect_power(&n, None, true)
        .map(|(p, k)| (p, k as i64))
        .unwrap_or((n, 1));
    if sign == 1 && k == 1 {
        return None;
    }
    Some((Expr::num(Number::from_bigint(p)), Expr::int(sign * k)))
}

/// The reference solves `b**g = r` (with `r` not a power of `b` itself)
/// through `exp(g*log(b)) - r`, and its polynomial conversion expands
/// `log(b)`: a perfect power base is taken as a power of its root
/// (`4**x - 8` is a quadratic in `2**x`, with the solutions `3/2` and
/// `3/2 + I*pi/log(2)`), and `(1/4)**x` is the reference `4**(-x)`.
/// Generators with such bases are rewritten as powers of the root, which
/// also gives a common base to `4**x + 2**x - 6`. When the equation is
/// linear in a single such generator and `r` has the base `b`, the
/// reference solves the exponents directly (`4**x - 4` gives `[1]`), and the
/// generator is left unchanged.
fn perfect_power_bases(f: &Expr, x: &Expr, gens: &[Expr]) -> Option<Expr> {
    let mut replacements = Vec::new();
    for g in gens {
        let Kind::Pow(base, exponent) = g.kind() else { continue };
        if has(base, x) || !has(exponent, x) {
            continue;
        }
        if let Some((root, power)) = power_of_integer_base(base) {
            replacements.push((g.clone(), pow(&root, &expand(&mul2(&power, exponent)))));
        }
    }
    if replacements.is_empty() {
        return None;
    }
    if let ([g], Kind::Pow(base, _)) = (gens, gens[0].kind()) {
        let u = Expr::dummy("_u", crate::Assumptions::NONE);
        let fu = subs1(f, g, &u);
        if let Some(c) = poly_coeffs(&fu, &u).filter(|c| c.len() == 2 && !has(&fu, x)) {
            let rhs = neg(&div(&c[0], &c[1]));
            // `(1/n)**g = r` is the reference `n**(-g) = r`, which it
            // inverts to `n**g = 1/r`.
            let (reference_base, rhs) = match base.as_number().and_then(|n| n.numer_denom()) {
                Some((p, q)) if p == num_bigint::BigInt::from(1) => {
                    (Expr::num(Number::from_bigint(q)), div(&Expr::one(), &rhs))
                }
                _ => (base.clone(), rhs),
            };
            let rhs_base = match rhs.kind() {
                Kind::Pow(b, _) => b.clone(),
                _ => rhs.clone(),
            };
            if rhs_base == reference_base {
                return None;
            }
        }
    }
    Some(subs(f, &replacements))
}

/// The reference inverts a generator `g` once for a symbol `t`
/// (`_vsolve(g - t, x)`) and substitutes each root of the polynomial in `g`
/// into that inverse. Root by root inversion is the fallback when the
/// symbolic inversion fails.
///
/// The reference then simplifies the substituted inverses when there are at
/// most two of them, which takes `log(8)/log(2)` to `3`. The native simplify
/// does not reach these forms, so for a positive real root of a power with
/// a base other than `E` the root is also inverted on its own, and the
/// exact rational exponents found that way (`log_ratio`) are kept.
///
/// Only exponential generators (`exp(a)`, `b**a` with `b` free of x) are
/// inverted this way. Other single generators of the native solver, such as
/// `log((y - 1)/(y + 1))`, are not polynomial generators of the reference,
/// which solves those equations by `_tsolve`; they are inverted root by root,
/// whose forms are the reference's (`(-x**2 - exp(2*C1))/(x**2 - exp(2*C1))`).
thread_local! {
    /// The exponential generators `g` whose symbolic inverse `g - t` is being
    /// solved, with their dummy `t`.
    static ACTIVE_INVERSIONS: std::cell::RefCell<Vec<(Expr, Expr)>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// `g - t` while the inverse of `g` is being computed for that dummy `t`.
///
/// Solving it through its generator asks for the inverse of `g` again (with
/// a fresh dummy), and so on until the depth bound; every level then fails
/// and is retried through the fallback before the direct inversion of the
/// last level succeeds, which costs twice the work per level and per
/// remaining depth. The direct inversion (`exp(a) = t` gives `a = log(t)`,
/// as the reference `_invert`) is the step that answers in the end.
fn is_active_inversion(g: &Expr, roots: &[Expr]) -> bool {
    let [root] = roots else { return false };
    ACTIVE_INVERSIONS.with(|active| active.borrow().iter().any(|(h, t)| h == g && t == root))
}

fn invert_generator(g: &Expr, roots: Vec<Expr>, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    let exponential = match g.kind() {
        Kind::Func(Func::Exp, _) => true,
        Kind::Pow(base, _) => !has(base, x),
        _ => false,
    };
    if !exponential {
        let mut out = Vec::new();
        for r in roots {
            out.extend(solve_rec(&sub(g, &r), x, depth + 1)?);
        }
        return Some(out);
    }
    let t = Expr::dummy("_t", crate::Assumptions::NONE);
    let inverse = {
        ACTIVE_INVERSIONS.with(|active| active.borrow_mut().push((g.clone(), t.clone())));
        struct Pop;
        impl Drop for Pop {
            fn drop(&mut self) { ACTIVE_INVERSIONS.with(|active| { active.borrow_mut().pop(); }); }
        }
        let _pop = Pop;
        solve_rec(&sub(g, &t), x, depth + 1)
    }
    .filter(|inverse| !inverse.is_empty() && inverse.iter().all(|i| has(i, &t)));
    let rational_base = matches!(g.kind(), Kind::Pow(..));
    // The reference simplifies the substituted inverses of `exp` when there
    // are at most two of them; a power of another base is solved through
    // `_tsolve`, whose results the top-level `_solve` always simplifies.
    let simplify_logs = rational_base || roots.len() <= 2;
    let mut out = Vec::new();
    for r in roots {
        if rational_base && crate::evalf::eval_f64(&r).is_some_and(|v| v > 0.0) {
            if let Some(values) = solve_rec(&sub(g, &r), x, depth + 1) {
                if !values.is_empty() && values.iter().all(|v| v.is_rational_number()) {
                    out.extend(values);
                    continue;
                }
            }
        }
        match &inverse {
            Some(inverse) => out.extend(inverse.iter().map(|i| {
                let value = subs1(i, &t, &r);
                if simplify_logs { expanded_log_root(&value) } else { value }
            })),
            None => out.extend(solve_rec(&sub(g, &r), x, depth + 1)?),
        }
    }
    Some(out)
}

/// The part of the reference simplify that the native simplify lacks for
/// the inverse of an exponential generator: `expand_log` takes the
/// logarithm of a radical root apart (`log(2*sqrt(2))` is `3*log(2)/2`), and
/// the quotient is distributed when that cancels a term to a rational
/// (`(log(sqrt(3)) + I*pi)/log(3)` is `1/2 + I*pi/log(3)`, while
/// `(log(3)/2 + I*pi)/log(2)` is kept).
fn expanded_log_root(value: &Expr) -> Expr {
    if !has_log(value) {
        return value.clone();
    }
    let separated = crate::simp::expand_log(value, false);
    // `log._eval_simplify` keeps a logarithm whose expansion is larger by
    // `count_ops` (`log(I*pi)`, not `log(pi) + I*pi/2`). The reference
    // reaches `log(-1 + log(4))` as `log(log(4*exp(-1)))`, whose inner
    // simplification gives `-1 + 2*log(2)`: the logarithms of positive
    // rationals in the argument are expanded. Only numeric arguments take
    // this path; symbolic roots keep the generic expansion below.
    if let Kind::Func(Func::Log, a) = value.kind() {
        if a.len() == 1 && crate::traverse::free_symbols(&a[0]).is_empty() && crate::simptools::methods::count_ops(&separated) > crate::simptools::methods::count_ops(value) {
            return crate::functions::func(Func::Log, vec![expand_rational_logs(&a[0])]);
        }
    }
    let distributed = expand(&separated);
    let cancels = distributed.is_rational_number()
        || distributed.is_add() && distributed.args().iter().any(|t| t.is_rational_number());
    if cancels { distributed } else { separated }
}

fn expand_rational_logs(e: &Expr) -> Expr {
    if let Kind::Func(Func::Log, a) = e.kind() {
        if a.len() == 1 && a[0].is_rational_number() && crate::assumptions::is_positive(&a[0]) == Some(true) {
            return crate::simp::expand_log(e, false);
        }
    }
    if e.is_atom() {
        return e.clone();
    }
    crate::traverse::map_children(e, &mut expand_rational_logs)
}

fn has_log(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found |= matches!(s.kind(), Kind::Func(Func::Log, _)));
    found
}
/// The reference's polynomial generator of `exp(r*t + c)`, for a rational
/// `r = p/q` with `|p| > 1` and `c` free of the unknown, is `exp(t/q)`: the
/// equation is a polynomial of degree `|p|` in it, so `exp(2*x) - 2` has the
/// two solutions `log(2)/2` and `log(2)/2 + I*pi`. Each root is inverted
/// through `exp(t/q)`. `None` when `g` is not such a power.
fn solve_exp_power(f: &Expr, g: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    // `b**(r*t + c)` with `b` free of the unknown is decomposed the same way
    // (`2**(3*x) - 8` has three solutions).
    let (base, argument) = match g.kind() {
        Kind::Func(Func::Exp, arguments) => (None, &arguments[0]),
        Kind::Pow(base, exponent) if !has(base, x) => (Some(base.clone()), exponent),
        _ => return None,
    };
    let raise = |e: &Expr| match &base {
        Some(base) => pow(base, e),
        None => exp(e),
    };
    let argument = expand(argument);
    let terms = if argument.is_add() { argument.args() } else { vec![argument.clone()] };
    let (dependent, constant): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(|t| has(t, x));
    let [term] = dependent.as_slice() else { return None; };
    let Kind::Mul(coefficient, _) = term.kind() else { return None; };
    if !coefficient.is_rational() {
        return None;
    }
    let (p, q) = coefficient.numer_denom()?;
    if num_traits::Signed::abs(&p) <= num_bigint::BigInt::from(1) {
        return None;
    }
    let tail = div(term, &Expr::num(coefficient.clone()));
    let generator = raise(&div(&tail, &Expr::num(Number::from_bigint(q))));
    let u = Expr::dummy("_u", crate::Assumptions::NONE);
    let power = pow(&u, &Expr::num(Number::from_bigint(p)));
    let fu = subs1(f, g, &mul2(&raise(&add(&constant)), &power));
    if has(&fu, x) || poly_coeffs(&numer_denom(&together(&fu)).0, &u).is_none() {
        return None;
    }
    let roots = solve_rec(&fu, &u, depth + 1)?.into_iter().filter(|r| !r.is_zero()).collect();
    invert_generator(&generator, roots, x, depth)
}

/// `f = c*h(x) + rest`, solve `h(x) = -rest/c` by inverting `h`.
fn invert(f: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    let terms = match f.kind() {
        Kind::Add(..) => f.args(),
        _ => vec![f.clone()],
    };
    let (dep, indep): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(|t| has(t, x));
    if dep.len() != 1 {
        return None;
    }
    let t = &dep[0];
    let (c, h) = match t.kind() {
        Kind::Mul(..) => {
            let (cs, hs): (Vec<Expr>, Vec<Expr>) = t.args().into_iter().partition(|a| !has(a, x));
            (mul(&cs), mul(&hs))
        }
        _ => (Expr::one(), t.clone()),
    };
    let rhs = div(&neg(&add(&indep)), &c);
    if let Kind::Pow(b, e) = h.kind() {
        if has(b, x) && has(e, x) {
            return self_power(f, b, e, &rhs, x, depth);
        }
    }
    let targets: Vec<(Expr, Expr)> = match h.kind() {
        Kind::Pow(b, e) if !has(e, x) => {
            let n = e.as_number()?;
            if n.is_integer() {
                return None; // polynomial-like: handled elsewhere
            }
            // b**(p/q) = rhs  =>  b = rhs**(q/p)
            vec![(b.clone(), pow(&rhs, &Expr::num(Number::ONE.div(n)?)))]
        }
        Kind::Pow(b, e) if !has(b, x) => {
            // b**g = rhs  =>  g = log(rhs)/log(b)
            vec![(e.clone(), log_ratio(&rhs, b))]
        }
        Kind::Func(func, a) if a.len() == 1 => {
            let g = a[0].clone();
            if let (Func::Exp, Some(value)) = (func, split_exp_inverse(&g, &rhs, x)) {
                return Some(vec![value]);
            }
            match func {
                Func::Exp => vec![(g, log(&rhs))],
                Func::Log => vec![(g, exp(&rhs))],
                Func::Sin => vec![(g.clone(), asin(&rhs)), (g, sub(&Expr::pi(), &asin(&rhs)))],
                Func::Cos => vec![(g.clone(), acos(&rhs)), (g, sub(&mul2(&Expr::int(2), &Expr::pi()), &acos(&rhs)))],
                Func::Tan => vec![(g, atan(&rhs))],
                // The reference reaches `sec` and `csc` through
                // `rewrite(exp)` and a simplified logarithm. Keep the short
                // inverse forms for the unit argument when they evaluate.
                Func::Sec | Func::Csc => {
                    // Reciprocal sine/cosine never has a finite zero.
                    if rhs.is_zero() || crate::assumptions::exact_queries(|| crate::assumptions::is_zero(&rhs)) == Some(true) {
                        return Some(vec![]);
                    }
                    // `_tsolve` rewrites sec/csc to exponentials before
                    // solving the polynomial in exp(I*x). Even when acos
                    // or asin evaluates exactly, division of only its
                    // principal branches by k loses the other kth roots
                    // (sec(2*x) + 1 needs both -pi/2 and pi/2).
                    if g != *x && trig::rational_affine_argument(&g, x).is_some() {
                        return sec_csc_multiple(matches!(func, Func::Sec), &g, &rhs, x, depth);
                    }
                    let r = div(&Expr::one(), &rhs);
                    let v = if matches!(func, Func::Sec) { acos(&r) } else { asin(&r) };
                    if contains_node(&v, |e| matches!(e.kind(), Kind::Func(Func::Acos | Func::Asin, _)))
                        || has_free_parameter(&v) {
                        return sec_csc_exp(matches!(func, Func::Sec), &g, &rhs, x, depth);
                    }
                    let w = if matches!(func, Func::Sec) { neg(&v) } else { sub(&Expr::pi(), &v) };
                    // a double root of the exp polynomial (`sec(x) + 1`
                    // gives `[pi]`) is one value
                    let turns = div(&sub(&w, &v), &mul2(&Expr::int(2), &Expr::pi()));
                    if turns.as_number().is_some_and(|n| n.is_integer()) {
                        vec![(g, v)]
                    } else if matches!(func, Func::Sec) {
                        vec![(g.clone(), w), (g, v)]
                    } else {
                        vec![(g.clone(), v), (g, w)]
                    }
                }
                Func::Cot => vec![(g, acot(&rhs))],
                Func::Sinh => vec![(g, asinh(&rhs))],
                Func::Cosh => vec![(g.clone(), acosh(&rhs)), (g, neg(&acosh(&rhs)))],
                Func::Tanh => vec![(g, atanh(&rhs))],
                Func::Asin => vec![(g, sin(&rhs))],
                Func::Acos => vec![(g, cos(&rhs))],
                Func::Atan => vec![(g, tan(&rhs))],
                Func::Asinh => vec![(g, sinh(&rhs))],
                Func::Atanh => vec![(g, tanh(&rhs))],
                // SymPy's _invert: LambertW(g) = r gives g = r*exp(r)
                Func::Special(crate::special::SpecialFn::LambertW) => vec![(g, mul2(&rhs, &exp(&rhs)))],
                _ => return None,
            }
        }
        _ => return None,
    };
    let mut out = Vec::new();
    for (g, v) in targets {
        if matches!(v.kind(), Kind::Const(Constant::NaN | Constant::ComplexInfinity)) {
            continue;
        }
        if g == *x {
            out.push(v);
        } else {
            out.extend(solve_rec(&sub(&g, &v), x, depth + 1)?);
        }
    }
    Some(out)
}

/// `sec(x) = r` or `csc(x) = r` through the reference `rewrite(exp)`: with
/// `t = exp(I*x)` the equation is `r*t**2 - 2*t + r = 0` (`sec`) or
/// `r*t**2 - 2*I*t - r = 0` (`csc`), and `x = -I*log(t)`. At the top level
/// the reference simplify takes a positive rational content `1/m` out of the
/// logarithm: `-I*log(1/3 - 2*sqrt(2)*I/3)` is `I*(log(3) - log(1 -
/// 2*sqrt(2)*I))`. Integer multiples use `sec_csc_multiple` below.
fn sec_csc_exp(secant: bool, g: &Expr, rhs: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    if g != x {
        return sec_csc_multiple(secant, g, rhs, x, depth);
    }
    if has(rhs, x) || rhs.is_zero() {
        return None;
    }
    let t = Expr::dummy("t", crate::Assumptions::NONE);
    let t2 = pow(&t, &Expr::int(2));
    let middle = if secant { mul2(&Expr::int(-2), &t) } else { mul(&[Expr::int(-2), Expr::i(), t.clone()]) };
    let last = if secant { rhs.clone() } else { neg(rhs) };
    let p = add(&[mul2(rhs, &t2), middle, last]);
    let roots = solve_rec(&p, &t, depth + 1)?;
    if roots.len() != 2 && !(roots.len() == 1 && rhs.is_number()) {
        return None;
    }
    let minus_i = neg(&Expr::i());
    let mut out = Vec::new();
    for r in roots {
        if r.is_zero() {
            return None;
        }
        let plain = mul2(&minus_i, &log(&r));
        if depth > 0 {
            out.push(plain);
            continue;
        }
        let (c, p) = crate::simptools::primitive(&r);
        let split = (|| {
            if !c.is_rational() || c.is_negative() || c.is_one() {
                return None;
            }
            let inverse = Number::ONE.div(&c)?;
            if !inverse.is_integer() || p.is_number() && p.as_number().is_some() {
                return None;
            }
            Some(mul2(&Expr::i(), &sub(&log(&Expr::num(inverse)), &log(&p))))
        })();
        out.push(split.unwrap_or(plain));
    }
    Some(out)
}

/// `sec(p*x/q+b) = rhs` or `csc(p*x/q+b) = rhs`, with positive coprime p,q.
/// Put t=exp(I*x/q): every quadratic root u of exp(I*(p*x/q+b)) gives
/// all p roots of t**p=u*exp(-I*b), then x=-I*q*log(t). These are the finite
/// exponential representatives in -pi*q < Re(x) <= pi*q, not a claim to
/// enumerate the infinite periodic solution families (solveset owns those).
/// Negative slopes normalize by sec's evenness and csc's oddness. The forms
/// are `-I*log(u)/k` for the principal root and `-I*log(w*u**(1/k))` for the
/// other roots of unity `w` (`-1`, or `-(1 -+ sqrt(3)*I)/2`).
fn sec_csc_multiple(secant: bool, g: &Expr, rhs: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    use num_traits::ToPrimitive;
    if has(rhs, x) || rhs.is_zero() {
        return None;
    }
    let (mut slope, mut offset) = trig::rational_affine_argument(g, x)?;
    let mut target = rhs.clone();
    if slope.is_negative() {
        slope = slope.neg(); offset = neg(&offset);
        if !secant { target = neg(&target); }
    }
    let (p, q) = slope.numer_denom()?;
    let k = p.to_i64()?;
    if k < 1 || offset.has_nonfinite() || target.has_nonfinite() {
        return None;
    }
    let scale = Expr::num(Number::from_bigint(q));
    // The quadratic expressions have the intrinsic condition target != 0:
    // sec/csc never vanish. Specializing the equation to zero was handled
    // above; unlike a symbolic slope, it cannot introduce a constant identity.
    let u = Expr::dummy("u", crate::Assumptions::NONE);
    let u2 = pow(&u, &Expr::int(2));
    let middle = if secant { mul2(&Expr::int(-2), &u) } else { mul(&[Expr::int(-2), Expr::i(), u.clone()]) };
    let last = if secant { target.clone() } else { neg(&target) };
    let p = add(&[mul2(&target, &u2), middle, last]);
    let roots = solve_rec(&p, &u, depth + 1)?;
    // rhs = +-1 gives a double root. It is one u value, but still
    // has all k distinct exponential preimages (u is never zero).
    if roots.len() != 2 && !(roots.len() == 1 && target.is_number()) {
        return None;
    }
    let minus_i = neg(&Expr::i());
    let phase = exp(&mul2(&minus_i, &offset));
    let mut out = Vec::new();
    for root in roots {
        let r = mul2(&root, &phase);
        if r.is_zero() {
            return None;
        }
        out.push(div(&mul(&[scale.clone(), minus_i.clone(), log(&r)]), &Expr::int(k)));
        if k == 2 {
            out.push(mul(&[scale.clone(), minus_i.clone(), log(&neg(&sqrt(&r)))]));
        } else if k == 3 {
            let cube = pow(&r, &Expr::num(Number::ONE.div(&Number::int(3))?));
            let half = Expr::num(Number::int(-1).div(&Number::int(2))?);
            let s3i = mul2(&sqrt(&Expr::int(3)), &Expr::i());
            for w in [sub(&Expr::one(), &s3i), add2(&Expr::one(), &s3i)] {
                out.push(mul(&[scale.clone(), minus_i.clone(), log(&mul(&[half.clone(), w, cube.clone()]))]));
            }
        } else {
            // All roots of t**k = r: the principal root times every kth
            // root of unity. Keep principal radicals compact, rather than
            // converting their nested algebraic coefficients to polar form.
            // No search budget may return only a prefix of these roots.
            let root = pow(&r, &Expr::rational(1, k));
            for j in 1..k {
                let angle = mul(&[Expr::int(2), Expr::pi(), Expr::rational(j, k)]);
                let unity = add2(&cos(&angle), &mul2(&Expr::i(), &sin(&angle)));
                out.push(mul(&[scale.clone(), minus_i.clone(), log(&mul2(&unity, &root))]));
            }
        }
    }
    Some(out)
}

/// `exp(c + k*x) = rhs` with a symbolic `c` free of the unknown and `k` one
/// of `1`, `-1` or a symbol. The reference polynomial generators split
/// `exp(c + k*x)` into `exp(c)*exp(k*x)` (expand with power_exp), so it solves
/// `exp(k*x) = R` with `R = rhs*exp(-c)`, `exp(-x)` through its numerator
/// `1 - R*exp(x)`: `exp(a - x) + 1` gives `log(-exp(a))` (not `a - I*pi`,
/// a different value for real `a`), `exp(a*x + b) - c` gives
/// `log(c*exp(-b))/a`.
fn split_exp_inverse(argument: &Expr, rhs: &Expr, x: &Expr) -> Option<Expr> {
    let argument = expand(argument);
    if !argument.is_add() {
        return None;
    }
    let (dependent, constant): (Vec<Expr>, Vec<Expr>) = argument.args().into_iter().partition(|t| has(t, x));
    let c = add(&constant);
    if free_symbols(&c).is_empty() {
        return None;
    }
    let [term] = dependent.as_slice() else { return None };
    let k = div(term, x);
    if has(&k, x) {
        return None;
    }
    let r = mul2(rhs, &exp(&neg(&c)));
    if k.is_one() {
        Some(log(&r))
    } else if k == Expr::minus_one() {
        Some(log(&div(&Expr::one(), &r)))
    } else if k.is_symbol() {
        Some(div(&log(&r), &k))
    } else {
        None
    }
}

/// `b**e = rhs` with both `b` and `e` depending on the unknown, as the
/// reference `_tsolve` handles `f(x)**g(x) == c`: the solutions of the log
/// form `e*log(b) - log(rhs)`, plus the checked candidates `b = +-r`,
/// `e = k` for each exact `rhs = r**k` (`rhs = 1` gives `e = 0` and
/// `b = +-1`). With `rhs = 0` the solutions are the zeros of `b` where `e`
/// does not vanish.
fn self_power(f: &Expr, b: &Expr, e: &Expr, rhs: &Expr, x: &Expr, depth: u32) -> Option<Vec<Expr>> {
    if rhs.is_zero() {
        let roots = solve_rec(b, x, depth + 1)?;
        return Some(roots.into_iter().filter(|r| !subs1(e, x, r).is_zero()).collect());
    }
    let logarithmic = sub(&mul2(e, &log(b)), &log(rhs));
    let mut log_route_failed = false;
    let mut solutions = solve_rec(&logarithmic, x, depth + 1).or_else(|| {
        // The reference `_tsolve` ends with a posify/force retry, under
        // which log((7**x)**x) = x**2*log(7): its candidates are kept only
        // when they satisfy the original equation.
        if !checksol::solve_flags().force { return None; }
        let forced = crate::simp::expand_log(&logarithmic, true);
        if forced == logarithmic { return None; }
        Some(solve_rec(&forced, x, depth + 1)?.into_iter()
            .filter(|candidate| verify_equation(f, &[(x.clone(), candidate.clone())])).collect())
    }).unwrap_or_else(|| { log_route_failed = true; Vec::new() });
    let mut targets = Vec::new();
    if rhs.is_one() {
        targets.extend([e.clone(), sub(b, &Expr::one()), add(&[b.clone(), Expr::one()])]);
    } else if let Some((p, q)) = rhs.as_number().filter(|n| n.is_rational()).and_then(|n| n.numer_denom()) {
        let zero = num_bigint::BigInt::from(0);
        if p > zero {
            for d in crate::ntheory::divisors(&p) {
                if d <= num_bigint::BigInt::from(1) {
                    continue;
                }
                let (mut k, mut power) = (0u32, num_bigint::BigInt::from(1));
                while &power * &d <= p {
                    power *= &d;
                    k += 1;
                }
                if power != p {
                    continue;
                }
                for s in crate::ntheory::divisors(&num_traits::Signed::abs(&q)) {
                    if num_traits::pow(s.clone(), k as usize) == q {
                        let r = div(&Expr::num(Number::from_bigint(d.clone())), &Expr::num(Number::from_bigint(s)));
                        targets.extend([sub(b, &r), add(&[b.clone(), r]), sub(e, &Expr::int(k as i64))]);
                    }
                }
            }
        }
    }
    for target in targets {
        for candidate in solve_rec(&target, x, depth + 1).unwrap_or_default() {
            if verify_equation(f, &[(x.clone(), candidate.clone())]) {
                solutions.push(candidate);
            }
        }
    }
    let mut unique: Vec<Expr> = Vec::new();
    for solution in solutions {
        if !unique.contains(&solution) {
            unique.push(solution);
        }
    }
    // Without the logarithmic route nothing proves the absence of
    // solutions: decline instead of reporting an empty list.
    if log_route_failed && unique.is_empty() {
        return None;
    }
    Some(unique)
}

/// The rational `k` with `b**k == r` for a negative base `b`, found from
/// `log(|r|)/log(|b|)` (the reference's `e_rat` guess in `_tsolve`); only a
/// `k` for which the principal power is exactly `r` is returned.
fn negative_base_ratio(r: &Expr, b: &Expr) -> Option<Expr> {
    let bv = crate::evalf::eval_f64(b)?;
    let rv = crate::evalf::eval_f64(r)?;
    if !(bv < 0.0) || rv == 0.0 || !rv.is_finite() || bv == -1.0 {
        return None;
    }
    let k = rv.abs().ln() / (-bv).ln();
    for q in 1..=12i64 {
        let p = (k * q as f64).round();
        if (k * q as f64 - p).abs() < 1e-9 * q as f64 && p.abs() < 1e6 {
            let cand = Expr::rational(p as i64, q);
            if crate::expand::expand(&sub(&pow(b, &cand), r)).is_zero() {
                return Some(cand);
            }
        }
    }
    None
}

/// log(r)/log(b), simplified to an exact rational when r is a rational power
/// of b (`log(2*sqrt(2))/log(2)` is `3/2`, as the reference simplifies it).
/// The candidate exponent is proven by exact evaluation of `b**k`.
fn log_ratio(r: &Expr, b: &Expr) -> Expr {
    if let (Some(rv), Some(bv)) = (crate::evalf::eval_f64(r), crate::evalf::eval_f64(b)) {
        if rv > 0.0 && bv > 0.0 && bv != 1.0 {
            let k = rv.ln() / bv.ln();
            for q in 1..=12i64 {
                let p = (k * q as f64).round();
                if (k * q as f64 - p).abs() < 1e-9 * q as f64 && p.abs() < 1e6 {
                    let cand = Expr::rational(p as i64, q);
                    if pow(b, &cand) == *r {
                        return cand;
                    }
                }
            }
        }
    }
    div(&log(r), &log(b))
}

// ----------------------------------------------------------------------
// systems
// ----------------------------------------------------------------------

/// Solution of a system: either a single mapping (linear systems, possibly
/// parametric) or a list of tuples (nonlinear systems).
#[derive(Debug, Clone)]
pub enum SystemSolution {
    Map(Vec<(Expr, Expr)>),
    Tuples(Vec<Vec<Expr>>),
    NoSolution,
}

fn is_linear_in(e: &Expr, vars: &[Expr]) -> bool {
    for v in vars {
        match poly_coeffs(e, v) {
            Some(c) if c.len() <= 2 => {
                if c.len() == 2 && vars.iter().any(|w| has(&c[1], w)) {
                    return false;
                }
            }
            _ => return false,
        }
    }
    true
}

/// Solve a system of equations for `vars`.
pub fn solve_system(eqs: &[Expr], vars: &[Expr]) -> Option<SystemSolution> {
    solve_system_with_check(eqs, vars, true)
}

pub fn solve_system_with_check(eqs: &[Expr], vars: &[Expr], check: bool) -> Option<SystemSolution> {
    // SymPy's solve drops an equation that has free symbols but none of the
    // unknowns: it cannot change the solution of the other equations, which
    // may then be given in terms of its symbols (solve((x - y, y - 3), x) is
    // {x: y}). Equations without any free symbol stay to be decided, and so
    // does one proved nonzero by the assumptions of its symbols (with p
    // nonzero, [x - 1, p] has no solution where SymPy returns {x: 1}).
    let original: Vec<Expr> = eqs.iter().map(as_expr)
        .filter(|e| free_symbols(e).is_empty() || vars.iter().any(|v| has(e, v))
            || crate::assumptions::exact_queries(|| crate::assumptions::is_zero(e)) == Some(false)).collect();
    let inverted: Vec<(Expr, Expr)> = original.iter().map(|e| system_invert(e, vars)).collect();
    let eqs: Vec<Expr> = inverted.iter().map(|(lhs, rhs)| expand(&sub(lhs, rhs))).filter(|e| !e.is_zero()).collect();
    if eqs.iter().all(|e| is_linear_in(e, vars)) {
        let result = linear_system(&eqs, vars)?;
        let result = if check { check_system_solution(result, &original, vars) } else { result };
        // The reference reports the system as nonlinear (tuples rather than
        // a dictionary) when an equation inverts to a single unknown in which
        // it is not linear (`[exp(x) - 2]` gives `[(log(2),)]`).
        let nonlinear = original.iter().zip(&inverted).any(|(equation, (lhs, _))|
            vars.contains(lhs) && !is_linear_in(&expand(equation), std::slice::from_ref(lhs)));
        return Some(match result {
            SystemSolution::Map(assignments) if nonlinear => SystemSolution::Tuples(vec![vars.iter().map(|variable|
                assignments.iter().find(|(key, _)| key == variable).map(|(_, value)| value.clone())
                    .unwrap_or_else(|| variable.clone())).collect()]),
            result => result,
        });
    }
    // A product vanishes on the UNION of its factor zero sets. Solving only
    // for one pivot divides by the other factors and loses whole components
    // when those factors contain another unknown. Preserve supplied factors
    // before expansion, branching on every positive integral factor. Unknown
    // parameter-only factors cannot be discarded as nonzero.
    for (index,equation) in original.iter().enumerate(){
        // Recover expanded components only when a pivot could divide by
        // another unknown. Ordinary equations with a proved nonzero linear
        // pivot need no factorization: that pivot already parametrizes all
        // solutions without a singular denominator.
        let recovered=if equation.is_add(){
            let coefficients:Vec<_>=vars.iter().filter_map(|v|poly_coeffs(equation,v)).collect();
            let regular_linear=coefficients.iter().any(|c|c.len()==2&&
                crate::assumptions::exact_queries(||crate::assumptions::is_zero(&c[1]))==Some(false));
            let singular_pivot=coefficients.iter().any(|c|c.len()>1&&vars.iter().any(|v|has(c.last().unwrap(),v)));
            if !regular_linear&&singular_pivot{crate::poly::api::factor(equation).ok()}else{None}
        }else{None};
        let equation=recovered.as_ref().unwrap_or(equation);
        let Kind::Mul(coefficient,factors)=equation.kind() else{continue;};
        if !coefficient.is_rational()||coefficient.is_zero()
            ||factors.iter().any(|(_,exponent)|!exponent.as_i64().is_some_and(|n|n>0)){continue;}
        let mut dependent=Vec::new();let mut supported=true;
        for (base,_) in factors{
            // A factor proved nonzero everywhere (exp(x)) has an empty zero set.
            if crate::assumptions::exact_queries(||crate::assumptions::is_zero(base))==Some(false){continue;}
            if vars.iter().any(|variable|has(base,variable)){dependent.push(base.clone());}
            else if crate::assumptions::exact_queries(||crate::assumptions::is_zero(base))!=Some(false){supported=false;break;}
        }
        if !supported||dependent.len()<2{continue;}
        let mut rows=Vec::new();
        for factor in dependent{
            let mut branch=original.clone();branch[index]=factor;
            match solve_system_with_check(&branch,vars,check)?{
                SystemSolution::NoSolution=>{},
                SystemSolution::Map(assignments)=>rows.push(vars.iter().map(|variable|
                    assignments.iter().find(|(key,_)|key==variable).map(|(_,value)|value.clone()).unwrap_or_else(||variable.clone())).collect()),
                SystemSolution::Tuples(mut branch_rows)=>rows.append(&mut branch_rows),
            }
        }
        rows.sort_by_cached_key(|row:&Vec<Expr>|row.iter().map(crate::sortkey::sort_key).collect::<Vec<_>>());
        rows.dedup();
        let result=if rows.is_empty(){SystemSolution::NoSolution}else{SystemSolution::Tuples(rows)};
        return Some(if check{check_system_solution(result,&original,vars)}else{result});
    }
    // nonlinear: successive substitution. A system polynomial in the
    // unknowns has its parameters in the coefficient domain, as SymPy's
    // solve_poly_system does.
    let generic = eqs.iter().all(|e| {
        let numerator = numer_denom(&together(e)).0;
        vars.iter().all(|v| poly_coeffs(&numerator, v).is_some())
    });
    let sols = substitute_solve(&eqs, vars, 0, check, generic)?;
    if sols.is_empty() {
        return Some(SystemSolution::NoSolution);
    }
    let mut tuples: Vec<Vec<Expr>> = sols
        .into_iter()
        .map(|m| vars.iter().map(|v| m.get(v).cloned().unwrap_or_else(|| v.clone())).collect())
        .collect();
    tuples.sort_by_cached_key(|t| t.iter().map(crate::sortkey::sort_key).collect::<Vec<_>>());
    tuples.dedup();
    let result = SystemSolution::Tuples(tuples);
    Some(if check { check_system_solution(result, &original, vars) } else { result })
}

/// The reference `_solve_system` replaces each equation `g` by `d - i` for
/// `(i, d) = _invert(g, *symbols)`: independent terms and factors move to the
/// right-hand side, a two-term sum of powers is combined by powsimp, equal
/// functions of the two terms are matched, a function with an inverse is
/// inverted (principal branch, as the reference; `sin` and `cos` are
/// kept), and a radical of a base
/// with several unknowns is raised to the reciprocal power
/// (`exp(x + y) - 1` gives `x + y`, making `[exp(x + y) - 1, x - y]` linear).
/// Only non-polynomial equations change: a polynomial is its own inversion
/// up to the order of its terms.
fn system_invert(eq: &Expr, vars: &[Expr]) -> (Expr, Expr) {
    let polynomial = |e: &Expr| vars.iter().all(|v| poly_coeffs(e, v).is_some());
    if polynomial(eq) || !vars.iter().any(|v| has(eq, v)) {
        return (eq.clone(), Expr::zero());
    }
    let dependent_on = |e: &Expr| vars.iter().any(|v| has(e, v));
    let ispow = |e: &Expr| matches!(e.kind(), Kind::Pow(..) | Kind::Func(Func::Exp, _));
    let options = crate::simptools::powsimp::PowsimpOptions::default();
    let mut lhs = eq.clone();
    let mut rhs = Expr::zero();
    for _ in 0..64 {
        let was = lhs.clone();
        loop {
            let (independent, dependent) = crate::decompose::as_independent(&lhs, vars, None);
            if lhs.is_add() {
                if independent.is_zero() { break; }
                lhs = dependent;
                rhs = sub(&rhs, &independent);
            } else {
                if independent.is_one() { break; }
                lhs = dependent;
                rhs = div(&rhs, &independent);
            }
        }
        if lhs.is_add() && rhs.is_zero() && lhs.args().len() == 2 && !polynomial(&lhs) {
            let mut pair = lhs.args();
            crate::logic::ordered(&mut pair);
            let (ai, ad) = crate::decompose::as_independent(&pair[0], vars, Some(false));
            let (bi, bd) = crate::decompose::as_independent(&pair[1], vars, Some(false));
            if ispow(&ad) || ispow(&bd) {
                let (a_base, a_exp) = crate::ops::as_base_exp(&ad);
                let (b_base, b_exp) = crate::ops::as_base_exp(&bd);
                if a_base == b_base && crate::simptools::extract_additively(&a_exp, &b_exp).is_none() {
                    lhs = crate::simptools::powsimp::powsimp(
                        &crate::simptools::powsimp::powdenest(&div(&ad, &bd), false, false), &options);
                    rhs = neg(&div(&bi, &ai));
                } else {
                    let ratio = div(&ad, &bd);
                    let combined = crate::simptools::powsimp::powsimp(&ratio, &options);
                    if combined != ratio {
                        lhs = combined;
                        rhs = neg(&div(&bi, &ai));
                    }
                }
            } else if ai == neg(&bi) {
                if let (Kind::Func(f, a), Kind::Func(g, b)) = (ad.kind(), bd.kind()) {
                    if f == g && a.len() == 1 && b.len() == 1 {
                        lhs = sub(&a[0], &b[0]);
                    }
                }
            }
        } else if lhs.is_mul() && lhs.args().iter().any(|a| ispow(a)) {
            lhs = crate::simptools::powsimp::powsimp(
                &crate::simptools::powsimp::powdenest(&lhs, false, false), &options);
        }
        if let Kind::Func(f, args) = lhs.kind() {
            let inverse = match f {
                Func::Exp => Some(Func::Log),
                Func::Log => Some(Func::Exp),
                Func::Tan => Some(Func::Atan),
                Func::Asin => Some(Func::Sin),
                Func::Acos => Some(Func::Cos),
                Func::Atan => Some(Func::Tan),
                _ => None,
            };
            // exp(u) = 0 has no solution; log(0) would make the equation zoo.
            let inverse = inverse.filter(|_| !(rhs.is_zero() && *f == Func::Exp));
            if let (Some(inverse), [argument]) = (inverse, args.as_slice()) {
                rhs = func(inverse, vec![rhs.clone()]);
                lhs = argument.clone();
            }
        }
        if let Kind::Pow(base, exponent) = lhs.kind() {
            let integer = exponent.as_i64();
            if !rhs.is_zero() && integer.is_some_and(|n| n < 0) {
                lhs = div(&Expr::one(), &lhs);
                rhs = div(&Expr::one(), &rhs);
            } else if integer.is_none() && exponent.is_number()
                && vars.iter().filter(|v| has(base, v)).count() > 1 {
                rhs = pow(&rhs, &div(&Expr::one(), exponent));
                lhs = base.clone();
            }
        }
        if lhs == was { break; }
    }
    if !dependent_on(&lhs) {
        return (eq.clone(), Expr::zero());
    }
    (lhs, rhs)
}

fn check_system_solution(result: SystemSolution, equations: &[Expr], variables: &[Expr]) -> SystemSolution {
    let valid = |assignments: &[(Expr, Expr)]| {
        if assignments.iter().any(|(variable, value)| !matches_assumptions(variable, value)) { return false; }
        equations.iter().all(|equation| {
            let (_, denominator) = numer_denom(&together(equation));
            let denominator = expand(&subs(&denominator, assignments));
            if denominator.is_zero() || matches!(denominator.kind(), Kind::Const(Constant::ComplexInfinity | Constant::NaN)) { return false; }
            let residual = subs(equation, assignments);
            if matches!(residual.kind(), Kind::Const(Constant::ComplexInfinity | Constant::NaN)) { return false; }
            verify_equation(equation, assignments)
        })
    };
    match result {
        SystemSolution::Map(assignments) => if valid(&assignments) { SystemSolution::Map(assignments) } else { SystemSolution::NoSolution },
        SystemSolution::Tuples(mut rows) => {
            rows.retain(|row| valid(&variables.iter().cloned().zip(row.iter().cloned()).collect::<Vec<_>>()));
            if rows.is_empty() { SystemSolution::NoSolution } else { SystemSolution::Tuples(rows) }
        }
        SystemSolution::NoSolution => SystemSolution::NoSolution,
    }
}

/// Whether a residual equation free of the unknowns rules a candidate out.
///
/// SymPy solves a system polynomial in the unknowns over the coefficient
/// domain of its parameters, so a residual that is a nonzero element of that
/// domain is generically inconsistent: solve([x - y, x - 1], [x]) is [].
/// For the other systems SymPy's checksol rejects a residual that depends on
/// parameters only when it has no function and no radical, the case where the
/// expanded numerator decides the zero test. A number must be proved nonzero.
fn residual_nonzero(e: &Expr, generic: bool) -> bool {
    if crate::assumptions::exact_queries(|| crate::assumptions::is_zero(e)) == Some(false) {
        return true;
    }
    if free_symbols(e).is_empty() {
        return false;
    }
    let numerator = expand(&numer_denom(&together(e)).0);
    if numerator.is_zero() {
        return false;
    }
    if generic {
        return true;
    }
    let mut transcendental = false;
    crate::traverse::preorder(&numerator, &mut |node| match node.kind() {
        Kind::Func(..) | Kind::Derivative(..) | Kind::Piecewise(..) => transcendental = true,
        Kind::Pow(_, exponent) if !exponent.is_integer_number() => transcendental = true,
        Kind::Mul(_, factors) if factors.iter().any(|(_, exponent)| !exponent.is_integer_number()) => transcendental = true,
        _ => {}
    });
    !transcendental
}

fn substitute_solve(eqs: &[Expr], vars: &[Expr], depth: u32, check: bool, generic: bool) -> Option<Vec<FxHashMap<Expr, Expr>>> {
    if depth > 6 {
        return None;
    }
    let eqs: Vec<Expr> = eqs.iter().map(|e| expand(e)).filter(|e| !e.is_zero()).collect();
    if eqs.is_empty() {
        return Some(vec![FxHashMap::default()]);
    }
    if eqs.iter().any(|e| vars.iter().all(|v| !has(e, v)) && residual_nonzero(e, generic)) {
        // A residual independent of the remaining unknowns that is nonzero
        // in the sense of residual_nonzero rules the candidate out. An
        // undecided one is retained, as SymPy keeps a candidate checksol
        // cannot decide.
        return Some(vec![]);
    }
    let eqs:Vec<_>=eqs.into_iter().filter(|e|vars.iter().any(|v|has(e,v))).collect();
    if eqs.is_empty(){return Some(vec![FxHashMap::default()]);}
    // choose the (equation, variable) pair that is simplest to solve:
    // Prefer degree 1, then a coefficient proved nonzero. This avoids
    // introducing a spurious exceptional parameter slice when a regular
    // pivot (for example z in x*y-z) already covers the entire graph.
    let mut best: Option<(usize, usize, (usize,bool))> = None;
    for (i, e) in eqs.iter().enumerate() {
        for (j, v) in vars.iter().enumerate() {
            if !has(e, v) {
                continue;
            }
            let coefficients=poly_coeffs(e,v);
            let deg=coefficients.as_ref().map(|c|c.len()-1).unwrap_or(100);
            let singular=coefficients.as_ref().and_then(|c|c.last()).is_none_or(|leading|
                crate::assumptions::exact_queries(||crate::assumptions::is_zero(leading))!=Some(false));
            let score=(deg,singular);
            if best.is_none_or(|(_, _, d)| score < d) {
                best = Some((i, j, score));
            }
        }
    }
    let (i, j, _) = best?;
    let v = &vars[j];
    let roots = solve_with_check(&eqs[i], v, check)?;
    let mut out = Vec::new();
    let rest_vars: Vec<Expr> = vars.iter().filter(|w| *w != v).cloned().collect();
    for r in roots {
        let rest: Vec<Expr> = eqs.iter().enumerate().filter(|(k, _)| *k != i).map(|(_, e)| subs1(e, v, &r)).collect();
        for mut m in substitute_solve(&rest, &rest_vars, depth + 1, check, generic)? {
            let val = expand(&subs(&r, &m.iter().map(|(a, b)| (a.clone(), b.clone())).collect::<Vec<_>>()));
            m.insert(v.clone(), val);
            out.push(m);
        }
    }
    Some(out)
}

/// Gauss–Jordan elimination over expressions.
fn linear_system(eqs: &[Expr], vars: &[Expr]) -> Option<SystemSolution> {
    let n = vars.len();
    // rows: coefficients | rhs
    let mut rows: Vec<Vec<Expr>> = Vec::new();
    for e in eqs {
        let mut row = Vec::with_capacity(n + 1);
        let mut rest = e.clone();
        for v in vars {
            let c = poly_coeffs(e, v)?;
            let a = c.get(1).cloned().unwrap_or_else(Expr::zero);
            row.push(a.clone());
            rest = expand(&sub(&rest, &mul2(&a, v)));
        }
        row.push(neg(&rest));
        rows.push(row);
    }
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut r = 0;
    for col in 0..n {
        let Some(p) = (r..rows.len()).find(|&i| !expand(&rows[i][col]).is_zero()) else { continue };
        rows.swap(r, p);
        let piv = rows[r][col].clone();
        for k in 0..=n {
            rows[r][k] = reduce_entry(&div(&rows[r][k], &piv));
        }
        for i in 0..rows.len() {
            if i == r {
                continue;
            }
            let f = rows[i][col].clone();
            if f.is_zero() {
                continue;
            }
            for k in 0..=n {
                let v = sub(&rows[i][k], &mul2(&f, &rows[r][k]));
                rows[i][k] = reduce_entry(&v);
            }
        }
        pivots.push((r, col));
        r += 1;
        if r == rows.len() {
            break;
        }
    }
    // inconsistent rows 0 = c
    for row in rows.iter().skip(r) {
        // The coefficients are those of SymPy's linear system domain: a
        // nonzero parametric right-hand side is generically inconsistent.
        if row[..n].iter().all(|c| c.is_zero()) && residual_nonzero(&row[n], true) {
            return Some(SystemSolution::NoSolution);
        }
    }
    let mut out = Vec::new();
    for (i, col) in pivots {
        let mut val = rows[i][n].clone();
        for k in 0..n {
            if k != col && !rows[i][k].is_zero() {
                val = sub(&val, &mul2(&rows[i][k], &vars[k]));
            }
        }
        // Like SymPy, keep a parametric value as one cancelled fraction:
        // expanding it would distribute the shared denominator over every
        // numerator term.
        let numeric = rows[i].iter().all(|c| crate::matrix::ratfunc::as_numer_denom(c).1.is_number());
        out.push((vars[col].clone(), if numeric { expand(&val) } else { val }));
    }
    Some(SystemSolution::Map(out))
}

/// Normalise an elimination entry, cancelling common polynomial factors so
/// fraction-free growth does not pile up (a**9/a**8 ...).
fn reduce_entry(v: &Expr) -> Expr {
    let c = crate::matrix::ratfunc::cancel_opts(v, true);
    if &c == v { together(&expand(v)) } else { c }
}

// ----------------------------------------------------------------------
// nsolve
// ----------------------------------------------------------------------

/// Newton's method for a system `f_i(vars) = 0` from `x0` (double
/// precision, complex arithmetic).
pub fn nsolve(eqs: &[Expr], vars: &[Expr], x0: &[f64]) -> Option<Vec<f64>> {
    nsolve_with_options(eqs, vars, x0, 1e-14, 100, true)
}

/// Newton iteration with an explicit residual check and iteration budget.
/// This real-valued interface never discards an imaginary component.
pub fn nsolve_with_options(eqs: &[Expr], vars: &[Expr], x0: &[f64], tolerance: f64, maxsteps: usize, verify: bool) -> Option<Vec<f64>> {
    use num_complex::Complex64 as C;
    let n = vars.len();
    if n == 0 || eqs.len() != n || x0.len() != n || !tolerance.is_finite() || tolerance <= 0.0 || maxsteps == 0 || x0.iter().any(|x| !x.is_finite()) {
        return None;
    }
    let fs: Vec<Expr> = eqs.iter().map(as_expr).collect();
    let jac: Vec<Vec<Expr>> = fs.iter().map(|f| vars.iter().map(|v| diff(f, v)).collect()).collect();
    let mut x: Vec<C> = x0.iter().map(|v| C::new(*v, 0.0)).collect();
    // The f64 evaluator addresses symbols by name: distinct assumed symbols
    // sharing a name must not silently alias.
    let names: Vec<String> = vars.iter().map(|v| v.to_string()).collect();
    if names.iter().enumerate().any(|(i, name)| names[..i].contains(name)) {
        return None;
    }
    for _ in 0..maxsteps {
        let env = |s: &str| names.iter().position(|nm| nm == s).map(|i| x[i]);
        let fx: Vec<C> = fs.iter().map(|f| crate::evalf::eval_with(f, &env)).collect::<Option<_>>()?;
        if fx.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return None;
        }
        let norm: f64 = fx.iter().map(|z| z.norm()).fold(0.0, f64::max);
        if norm <= tolerance {
            break;
        }
        let mut a: Vec<Vec<C>> = jac
            .iter()
            .map(|row| row.iter().map(|d| crate::evalf::eval_with(d, &env)).collect::<Option<Vec<C>>>())
            .collect::<Option<_>>()?;
        // solve J dx = -f by Gaussian elimination with partial pivoting
        let mut b: Vec<C> = fx.iter().map(|z| -z).collect();
        for col in 0..n {
            if a.iter().flatten().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
                return None;
            }
            let p = (col..n).max_by(|&i, &j| a[i][col].norm().total_cmp(&a[j][col].norm()))?;
            if a[p][col].norm() < 1e-300 {
                return None;
            }
            a.swap(col, p);
            b.swap(col, p);
            for i in col + 1..n {
                let f = a[i][col] / a[col][col];
                for k in col..n {
                    let t = a[col][k];
                    a[i][k] -= f * t;
                }
                let t = b[col];
                b[i] -= f * t;
            }
        }
        let mut dx = vec![C::new(0.0, 0.0); n];
        for i in (0..n).rev() {
            let mut s = b[i];
            for k in i + 1..n {
                s -= a[i][k] * dx[k];
            }
            dx[i] = s / a[i][i];
        }
        let step: f64 = dx.iter().map(|z| z.norm()).fold(0.0, f64::max);
        for i in 0..n {
            x[i] += dx[i];
        }
        if step <= tolerance * x.iter().map(|z| z.norm()).fold(1.0, f64::max) {
            break;
        }
    }
    if x.iter().any(|z| !z.re.is_finite() || !z.im.is_finite() || z.im.abs() > tolerance * z.re.abs().max(1.0)) {
        return None;
    }
    if verify {
        let env = |s: &str| names.iter().position(|nm| nm == s).map(|i| x[i]);
        for f in &fs {
            let residual = crate::evalf::eval_with(f, &env)?;
            if !residual.norm().is_finite() || residual.norm_sqr() > tolerance {
                return None;
            }
        }
    }
    Some(x.into_iter().map(|z| z.re).collect())
}

/// Bracketed real bisection. Residual verification rejects discontinuity poles.
pub fn nsolve_bisect(equation: &Expr, variable: &Expr, mut a: f64, mut b: f64, tolerance: f64, maxsteps: usize, verify: bool) -> Option<f64> {
    if !a.is_finite() || !b.is_finite() || a >= b || !tolerance.is_finite() || tolerance <= 0.0 || maxsteps == 0 {
        return None;
    }
    let expression = as_expr(equation);
    let name = variable.to_string();
    let evaluate = |x: f64| {
        let value = crate::evalf::eval_with(&expression, &|s| (s == name).then(|| num_complex::Complex64::new(x, 0.0)))?;
        (value.re.is_finite() && value.im == 0.0).then_some(value.re)
    };
    let mut fa = evaluate(a)?;
    let fb = evaluate(b)?;
    if fa == 0.0 { return Some(a); }
    if fb == 0.0 { return Some(b); }
    if fa.signum() == fb.signum() { return None; }
    let mut midpoint = a;
    for _ in 0..maxsteps {
        midpoint = a / 2.0 + b / 2.0;
        let fm = evaluate(midpoint)?;
        if fm == 0.0 || (b-a).abs() <= tolerance * midpoint.abs().max(1.0) { break; }
        if fa.signum() == fm.signum() { a = midpoint; fa = fm; } else { b = midpoint; }
    }
    let residual = evaluate(midpoint)?;
    (!verify || residual * residual <= tolerance).then_some(midpoint)
}

/// Convenience: free symbols of a set of equations.
pub fn unknowns(eqs: &[Expr]) -> Vec<Expr> {
    let mut v: Vec<Expr> = Vec::new();
    for e in eqs {
        for s in free_symbols(e) {
            if !v.contains(&s) {
                v.push(s);
            }
        }
    }
    v.sort_by_cached_key(crate::sortkey::sort_key);
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn s(e: &str) -> String {
        let x = Expr::symbol("x");
        let r = solve(&parse_str(e).unwrap(), &x).unwrap();
        format!("[{}]", r.iter().map(|e| e.to_string()).collect::<Vec<_>>().join(", "))
    }

    #[test]
    fn univariate() {
        assert_eq!(s("x-2"), "[2]");
        assert_eq!(s("2*x+3"), "[-3/2]");
        assert_eq!(s("x**2-4"), "[-2, 2]");
        assert_eq!(s("x**2+1"), "[-I, I]");
        assert_eq!(s("x**2-2"), "[-sqrt(2), sqrt(2)]");
        assert_eq!(s("x**2+x-1"), "[-1/2 + sqrt(5)/2, -sqrt(5)/2 - 1/2]");
        assert_eq!(s("x**3-6*x**2+11*x-6"), "[1, 2, 3]");
        assert_eq!(s("x**3-1"), "[1, -1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2]");
        assert_eq!(s("x**4-1"), "[-1, 1, -I, I]");
        assert_eq!(s("exp(x)-2"), "[log(2)]");
        assert_eq!(s("log(x)-1"), "[E]");
        assert_eq!(s("sin(x)"), "[0, pi]");
        assert_eq!(s("cos(x)-1"), "[0, 2*pi]");
        assert_eq!(s("sin(x)-1/2"), "[pi/6, 5*pi/6]");
        assert_eq!(s("sqrt(x)-3"), "[9]");
        assert_eq!(s("1/x-2"), "[1/2]");
        assert_eq!(s("(x-1)/(x-2)"), "[1]");
        assert_eq!(s("exp(2*x)-3*exp(x)+2"), "[0, log(2)]");
        assert_eq!(s("2**x-8"), "[3]");
        assert_eq!(s("x**5-x"), "[-1, 0, 1, -I, I]");
        assert_eq!(s("x**2-2*x+1"), "[1]");
    }

    #[test]
    fn optional_candidate_checks_and_exact_assumptions() {
        use crate::symbol::Assumptions as A;
        let positive = Expr::symbol_with("positive", A::POSITIVE);
        let equation = sub(&positive.powi(2), &Expr::int(4));
        assert_eq!(solve(&equation, &positive).unwrap(), vec![Expr::int(2)]);
        assert_eq!(solve_with_check(&equation, &positive, false).unwrap(), vec![Expr::int(-2), Expr::int(2)]);
        let real = Expr::symbol_with("real", A::REAL);
        assert!(solve(&add2(&real.powi(2), &Expr::one()), &real).unwrap().is_empty());
        let x = Expr::symbol("x");
        let sinc = div(&sin(&x), &x);
        assert_eq!(solve(&sinc, &x).unwrap(), vec![Expr::pi()]);
        assert_eq!(solve_with_check(&sinc, &x, false).unwrap(), vec![Expr::zero(), Expr::pi()]);
        let candidate = sub(&sqrt(&Expr::int(2)), &Expr::int(2));
        assert_eq!(solve(&sub(&positive, &candidate), &positive).unwrap(), Vec::<Expr>::new());
    }

    #[test]
    fn sec_csc_integer_multiples_keep_every_exponential_preimage() {
        let x = Expr::symbol("x");
        for (equation, count) in [("sec(2*x) + 1", 2), ("csc(2*x) - 2", 4),
                                  ("sec(3*x) + 1", 3), ("sec(4*x) - 3", 8),
                                  ("csc(5*x) - 2", 10)] {
            let equation = parse_str(equation).unwrap();
            let roots = solve(&equation, &x).unwrap();
            assert_eq!(roots.len(), count, "{equation}");
            for (i, root) in roots.iter().enumerate() {
                assert!(crate::evalf_mp::max_precision(root).is_none(), "{root}");
                let value = crate::evalf::eval_c64(root).unwrap();
                assert!(crate::evalf::eval_c64(&subs1(&equation, &x, root)).unwrap().norm() < 1e-10,
                        "{equation} at {root}");
                for other in &roots[..i] {
                    assert!((value - crate::evalf::eval_c64(other).unwrap()).norm() > 1e-10);
                }
            }
        }
        assert_eq!(s("sec(2*x) + 1"), "[-pi/2, pi/2]");
    }

    #[test]
    fn sec_csc_shifted_rational_slopes_keep_all_preimages() {
        let x = Expr::symbol("x");
        for (source, count) in [
            ("sec(2*x+1)+1", 2), ("csc(2*x+1)-2", 4),
            ("sec(-2*x+1)+1", 2), ("csc(-2*x+1)-2", 4),
            ("sec(x/2+1)-2", 2), ("csc(3*x/2+1)-2", 6),
            ("sec(-2*x/3+I)-3", 4), ("csc(2*x/3+I)-1", 2),
        ] {
            let equation = parse_str(source).unwrap();
            let roots = solve(&equation, &x).unwrap();
            assert_eq!(roots.len(), count, "{source}");
            for (index, root) in roots.iter().enumerate() {
                assert!(crate::evalf_mp::max_precision(root).is_none());
                assert!(crate::evalf::eval_c64(&subs1(&equation, &x, root)).unwrap().norm() < 1e-9,
                        "{source} at {root}");
                let value = crate::evalf::eval_c64(root).unwrap();
                for previous in &roots[..index] {
                    assert!((value-crate::evalf::eval_c64(previous).unwrap()).norm() > 1e-9);
                }
            }
        }
    }

    #[test]
    fn trig_affine_decomposition_rejects_symbolic_or_nonlinear_slopes() {
        let x = Expr::symbol("x");
        for source in ["a*x+b", "x**2+1", "sin(x)+1", "x**100000000+1", "b"] {
            assert!(trig::rational_affine_argument(&parse_str(source).unwrap(), &x).is_none(), "{source}");
        }
        let (slope, offset) = trig::rational_affine_argument(&parse_str("-3*x/2+b").unwrap(), &x).unwrap();
        assert_eq!(slope, Number::int(-3).div(&Number::int(2)).unwrap());
        assert_eq!(offset, Expr::symbol("b"));
    }

    #[test]
    fn logarithmic_candidate_assumptions_are_certified() {
        use crate::symbol::Assumptions as A;
        let real = Expr::symbol_with("real", A::REAL);
        for head in [sec as fn(&Expr) -> Expr, csc] {
            let equation = sub(&head(&mul2(&Expr::int(4), &real)), &Expr::rational(1, 2));
            let roots = solve(&equation, &real).unwrap();
            assert!(roots.is_empty(), "{equation}: {roots:?}");
        }
        let positive = Expr::symbol_with("positive", A::POSITIVE);
        let negative = Expr::symbol_with("negative", A::NEGATIVE);
        let small = log(&sub(&Expr::int(2), &sqrt(&Expr::int(3))));
        assert!(solve(&sub(&positive, &small), &positive).unwrap().is_empty());
        assert_eq!(solve(&sub(&negative, &small), &negative).unwrap(), vec![small]);
        let t = Expr::symbol("t");
        let undecided = neg(&mul2(&Expr::i(), &log(&abs(&t))));
        assert!(matches_assumptions(&real, &undecided));
    }

    #[test]
    fn systems() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let e1 = parse_str("x+y-3").unwrap();
        let e2 = parse_str("x-y-1").unwrap();
        match solve_system(&[e1, e2], &[x.clone(), y.clone()]).unwrap() {
            SystemSolution::Map(m) => assert_eq!(format!("{:?}", m.iter().map(|(a, b)| format!("{a}: {b}")).collect::<Vec<_>>()), "[\"x: 2\", \"y: 1\"]"),
            other => panic!("{other:?}"),
        }
        let e1 = parse_str("x**2+y**2-1").unwrap();
        let e2 = parse_str("x-y").unwrap();
        match solve_system(&[e1, e2], &[x.clone(), y.clone()]).unwrap() {
            SystemSolution::Tuples(t) => assert_eq!(t.len(), 2),
            other => panic!("{other:?}"),
        }
        let r = nsolve(&[parse_str("cos(x) - x").unwrap()], &[x], &[1.0]).unwrap();
        assert!((r[0] - 0.739_085_133_215_160_6).abs() < 1e-14);
    }

    #[test]
    fn product_system_preserves_every_factor_component(){
        let x=Expr::symbol("x");let y=Expr::symbol("y");
        let equation=parse_str("(x-1)*(y-2)").unwrap();
        let SystemSolution::Tuples(rows)=solve_system(&[equation.clone()],&[x.clone(),y.clone()]).unwrap() else{panic!("component union expected");};
        assert_eq!(rows.len(),2);
        assert!(rows.contains(&vec![Expr::one(),y.clone()]));
        assert!(rows.contains(&vec![x.clone(),Expr::int(2)]));
        let SystemSolution::Tuples(rows)=solve_system(&[equation,parse_str("x+y-4").unwrap()],&[x,y]).unwrap() else{panic!("both branches expected");};
        assert_eq!(rows.len(),2);
        assert!(rows.contains(&vec![Expr::one(),Expr::int(3)]));
        assert!(rows.contains(&vec![Expr::int(2),Expr::int(2)]));
    }

    #[test]
    fn boolean_system_constraints_are_equations_not_discarded_atoms() {
        let x = Expr::symbol("x");
        let equation = sub(&x, &Expr::int(2));
        assert!(matches!(solve_system(&[equation.clone(), Expr::boolean(false)], &[x.clone()]),
            Some(SystemSolution::NoSolution)));
        let Some(SystemSolution::Map(assignments)) = solve_system(
            &[equation, Expr::boolean(true)], &[x.clone()]) else { panic!("x = 2 expected"); };
        assert_eq!(assignments, vec![(x, Expr::int(2))]);
        assert!(as_expr(&Expr::boolean(true)).is_zero());
        assert!(as_expr(&Expr::boolean(false)).is_one());
    }
}
