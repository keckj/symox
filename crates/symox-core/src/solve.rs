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
pub use recast::{solve_atom,solve_atom_system};
pub(crate) mod lambert;

/// `lhs - rhs` for an equation (`Eq`) or the expression itself.
pub fn as_expr(eq: &Expr) -> Expr {
    match eq.kind() {
        Kind::Rel(RelOp::Eq, a, b) => sub(a, b),
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
    let f = as_expr(eq);
    if let Kind::Bool(b) = eq.kind() {
        return if *b { None } else { Some(vec![]) };
    }
    let mut sols = solve_rec(&f, x, 0)?;
    // Infinities reached by log(0), for example, are limiting values rather
    // than finite solutions. This restriction also applies with check=False.
    sols.retain(|s| !matches!(s.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN));
    if check {
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
        sols.retain(|s| verify(&f, x, s));
    }
    sort_solutions(&mut sols);
    Some(sols)
}

fn matches_assumptions(variable: &Expr, value: &Expr) -> bool {
    use crate::assumptions as query;
    let Some(symbol) = variable.as_sym() else { return true; };
    let assumptions = symbol.facts();
    if assumptions.len()==1 { return true; } // default commutativity only
    query::exact_queries(|| {
        let algebraic=query::quadratic_facts(value);
        let integer = if value.as_number().is_some_and(Number::is_float) { None } else { query::is_integer(value).or_else(||algebraic.and_then(|f|f.integer)) };
        let rational = if value.as_number().is_some_and(Number::is_float) { None } else { query::is_rational(value).or_else(||algebraic.and_then(|f|f.rational)) };
        let real = query::is_real(value).or_else(||algebraic.and_then(|f|f.real)).or_else(|| {
            let imaginary_part = im(value);
            if matches!(imaginary_part.kind(), Kind::Func(Func::Im, _)) { None } else { query::is_zero(&imaginary_part) }
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
        let positive=query::is_positive(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s>0));
        let negative=query::is_negative(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s<0));
        let nonnegative=query::is_nonnegative(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s>=0));
        let nonpositive=query::is_nonpositive(value).or_else(||algebraic.and_then(|f|f.sign).map(|s|s<=0));
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
    let precision = assignments.iter().filter_map(|(_, value)| crate::evalf_mp::max_precision(value))
        .chain(crate::evalf_mp::max_precision(f)).max();
    let Some(precision) = precision else {
        if let Some(zero) = crate::assumptions::exact_queries(|| crate::assumptions::is_zero(&v)) {
            return zero;
        }
        if !crate::functions::is_numeric_value(&v) { return true; }
        let mut previous = None;
        for bits in [160, 320] {
            let Some(value) = crate::evalf_mp::eval_mp(&v, bits) else { return true; };
            let scale = residual_scale(f, assignments, bits);
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
    // polynomial in x
    if let Some(c) = poly_coeffs(&num, x) {
        return poly_roots(&c).or_else(|| complete_roots(&num, x, &c))
            .or_else(|| algebraic::roots(&num, x));
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
            out.extend(solve_rec(b, x, depth + 1)?);
        }
        return Some(out);
    }
    // polynomial in a single generator g(x)
    if let Some(r) = solve_by_generator(&num, x, depth) {
        return Some(r);
    }
    // isolate the only term depending on x and invert
    if let Some(r) = invert(&num, x, depth) {
        return Some(r);
    }
    if let Some(r) = lambert::solve(&num, x) {
        return Some(r);
    }
    // several logarithms: combine them (log(a) - log(b) = log(a/b))
    let lc = crate::simp::logcombine(&num, true);
    if lc != num {
        if let Some(r) = solve_rec(&lc, x, depth + 1) {
            return Some(r);
        }
    }
    let ex = expand(&num);
    if ex != num {
        return solve_rec(&ex, x, depth + 1);
    }
    None
}

/// The full `roots` machinery (quartics, decomposition such as x**8 in
/// x**2, factorisation), accepted only when it finds every root: SymPy's
/// `roots` may return an incomplete set, which `solve` must not report.
fn complete_roots(num: &Expr, x: &Expr, c: &[Expr]) -> Option<Vec<Expr>> {
    let deg = c.iter().rposition(|t| !t.is_zero())? as u32;
    let rs = crate::poly::roots::roots(num, std::slice::from_ref(x)).ok()?;
    (rs.iter().map(|r| r.1).sum::<u32>() == deg).then(|| rs.into_iter().map(|r| r.0).collect())
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
        let s = sqrt(&r);
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
    let s = sqrt(&d);
    vec![div(&sub(&neg(b), &s), &two_a), div(&add2(&neg(b), &s), &two_a)]
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
    let all_exp = gens.iter().all(|g| matches!(g.kind(), Kind::Func(Func::Exp, _)));
    if all_exp && gens.len() > 1 {
        // rational multiples of x in the exponents
        let mut ks = Vec::new();
        for g in &gens {
            let Kind::Func(_, a) = g.kind() else { unreachable!() };
            let k = diff(&a[0], x);
            if has(&k, x) || !k.is_rational_number() {
                return None;
            }
            ks.push(k);
        }
        let mut q = num_bigint::BigInt::from(1);
        for k in &ks {
            q = num_integer::Integer::lcm(&q, &k.as_number().unwrap().numer_denom().unwrap().1);
        }
        let qn = Number::from_bigint(q);
        let mut g2 = f.clone();
        for (g, k) in gens.iter().zip(&ks) {
            let Kind::Func(_, a) = g.kind() else { unreachable!() };
            let c0 = expand(&sub(&a[0], &mul2(k, x)));
            let power = Expr::num(k.as_number().unwrap().mul(&qn));
            g2 = subs1(&g2, g, &mul2(&exp(&c0), &pow(&u, &power)));
        }
        if has(&g2, x) {
            return None;
        }
        let roots = solve_rec(&g2, &u, depth + 1)?;
        let mut out = Vec::new();
        for r in roots {
            if r.is_zero() {
                continue;
            }
            out.push(expand(&mul2(&Expr::num(qn.clone()), &log(&r))));
        }
        return Some(out);
    }
    if gens.len() != 1 {
        return None;
    }
    let g = &gens[0];
    let fu = subs1(f, g, &u);
    if has(&fu, x) {
        return None;
    }
    if poly_coeffs(&fu, &u).is_none() {
        return None;
    }
    let roots = solve_rec(&fu, &u, depth + 1)?;
    let mut out = Vec::new();
    for r in roots {
        out.extend(solve_rec(&sub(g, &r), x, depth + 1)?);
    }
    Some(out)
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
            match func {
                Func::Exp => vec![(g, log(&rhs))],
                Func::Log => vec![(g, exp(&rhs))],
                Func::Sin => vec![(g.clone(), asin(&rhs)), (g, sub(&Expr::pi(), &asin(&rhs)))],
                Func::Cos => vec![(g.clone(), acos(&rhs)), (g, sub(&mul2(&Expr::int(2), &Expr::pi()), &acos(&rhs)))],
                Func::Tan => vec![(g, atan(&rhs))],
                Func::Cot => vec![(g, acot(&rhs))],
                Func::Sinh => vec![(g, asinh(&rhs))],
                Func::Cosh => vec![(g.clone(), acosh(&rhs)), (g, neg(&acosh(&rhs)))],
                Func::Tanh => vec![(g, atanh(&rhs))],
                Func::Asin => vec![(g, sin(&rhs))],
                Func::Acos => vec![(g, cos(&rhs))],
                Func::Atan => vec![(g, tan(&rhs))],
                Func::Asinh => vec![(g, sinh(&rhs))],
                Func::Atanh => vec![(g, tanh(&rhs))],
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

/// log(r)/log(b), simplified to an exact rational when r is a power of b.
fn log_ratio(r: &Expr, b: &Expr) -> Expr {
    if let (Some(rv), Some(bv)) = (crate::evalf::eval_f64(r), crate::evalf::eval_f64(b)) {
        if rv > 0.0 && bv > 0.0 {
            let k = rv.ln() / bv.ln();
            let kr = k.round();
            if (k - kr).abs() < 1e-12 && kr.abs() < 1e6 {
                let cand = Expr::int(kr as i64);
                if pow(b, &cand) == *r {
                    return cand;
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
    let original: Vec<Expr> = eqs.iter().map(as_expr).collect();
    let eqs: Vec<Expr> = eqs.iter().map(|e| expand(&as_expr(e))).filter(|e| !e.is_zero()).collect();
    if eqs.iter().all(|e| is_linear_in(e, vars)) {
        let result = linear_system(&eqs, vars)?;
        return Some(if check { check_system_solution(result, &original, vars) } else { result });
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
    // nonlinear: successive substitution
    let sols = substitute_solve(&eqs, vars, 0, check)?;
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

fn substitute_solve(eqs: &[Expr], vars: &[Expr], depth: u32, check: bool) -> Option<Vec<FxHashMap<Expr, Expr>>> {
    if depth > 6 {
        return None;
    }
    let eqs: Vec<Expr> = eqs.iter().map(|e| expand(e)).filter(|e| !e.is_zero()).collect();
    if eqs.is_empty() {
        return Some(vec![FxHashMap::default()]);
    }
    if eqs.iter().any(|e| vars.iter().all(|v| !has(e, v))
        &&crate::assumptions::exact_queries(||crate::assumptions::is_zero(e))==Some(false)) {
        // Only a PROVED nonzero equation independent of the remaining
        // unknowns proves inconsistency. An unresolved parameter constraint
        // is retained by the original-equation verification, as in solve's
        // generic conditional-candidate convention.
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
        for mut m in substitute_solve(&rest, &rest_vars, depth + 1, check)? {
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
        if row[..n].iter().all(|c| c.is_zero())
            &&crate::assumptions::exact_queries(||crate::assumptions::is_zero(&row[n]))==Some(false) {
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
}
