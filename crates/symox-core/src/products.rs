//! Symbolic products: `product(f, (k, a, b))` and held `Product` nodes.
//!
//! Follows SymPy's `Product.doit`/`_eval_product`: constant terms become
//! powers, polynomial terms split over their roots into rising factorials,
//! sums are factored over a common denominator, products and powers are
//! split factor by factor, `b**g(k)` turns into `b**Sum(g(k))`, and short
//! definite ranges are multiplied out. Backwards ranges use Karr's
//! convention `Product(f, (k, a, b)) = 1/Product(f, (k, b + 1, a - 1))`.

use num_traits::{Signed, Zero};
use rustc_hash::FxHashMap;

use crate::assumptions::{is_finite, Signs};
use crate::binders::{depends_on as has, held, limits};
use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, add2, inv, mul, pow, sub};
use crate::special::rf;
use crate::subs::{subs1, xreplace};

/// Ranges of at most this many factors (`b - a < DIRECT`) are multiplied out.
const DIRECT: i64 = 100;

/// `Product(f, (k, a, b))` placeholder.
pub fn product_unevaluated(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Expr {
    held(Func::Product, f, vec![vec![k.clone(), a.clone(), b.clone()]]).expect("valid product limits")
}

/// Sign facts of a held `Product` (SymPy's `Product._eval_is_*`): an empty
/// range gives one; otherwise a zero term gives zero, and over finite limits
/// a positive, nonnegative or real term gives a product of the same kind.
pub fn signs(args: &[Expr]) -> Option<Signs> {
    let all = limits(args);
    if all.iter().any(|l| l.len() == 3 && sub(&l[1], &l[2]).is_one()) {
        return Some(Signs::POS);
    }
    if crate::assumptions::is_zero(&args[0]) == Some(true) {
        return Some(Signs::ZERO);
    }
    let finite = |l: &Vec<Expr>| {
        l.len() == 3
            && l[1..].iter().all(|b| match is_finite(b) {
                Some(f) => f,
                None => is_finite(&l[0]) == Some(true),
            })
    };
    if !all.iter().all(finite) {
        return None;
    }
    Some(match crate::assumptions::real_signs(&args[0])? {
        s @ (Signs::POS | Signs::ZERO | Signs::NONNEG) => s,
        _ => Signs::ANY,
    })
}

/// `term**count` with SymPy's `Pow._eval_power` rule for a base that is a
/// power with a real rational exponent of absolute value below one:
/// `sqrt(a)**n` is `a**(n/2)` for every `n`.
fn power(term: &Expr, count: &Expr) -> Expr {
    if let Kind::Pow(b, e) = term.kind() {
        if let Some(r) = e.as_number().filter(|r| r.is_rational()) {
            if r.abs().cmp_value(&Number::int(1)) == Some(std::cmp::Ordering::Less) {
                return pow(b, &mul(&[e.clone(), count.clone()]));
            }
        }
    }
    pow(term, count)
}

/// `product(f, *limits)`: SymPy's `Product(f, *limits).doit(deep=False)`.
pub fn product(f: &Expr, limits: Vec<Vec<Expr>>) -> Result<Expr, String> {
    let e = held(Func::Product, f, limits)?;
    let Kind::Func(_, args) = e.kind() else { unreachable!() };
    Ok(doit(args, false))
}

/// `Product.doit`: the arguments of a held Product node.
pub fn doit(args: &[Expr], deep: bool) -> Expr {
    // SymPy first renames each index to a Dummy carrying the assumptions
    // shared by both endpoints (integer, positive, ...), then renames back.
    let all = limits(args);
    let mut forward = FxHashMap::default();
    let mut undo = FxHashMap::default();
    for limit in &all {
        if let Some(dummy) = inherited_dummy(&limit[0], &limit[1], &limit[2]) {
            forward.insert(limit[0].clone(), dummy.clone());
            undo.insert(dummy, limit[0].clone());
        }
    }
    if !forward.is_empty() {
        let renamed = xreplace(&Expr::from_kind(Kind::Func(Func::Product, args.to_vec())), &forward);
        let Kind::Func(_, renamed_args) = renamed.kind() else { unreachable!() };
        return xreplace(&doit_renamed(renamed_args, deep), &undo);
    }
    doit_renamed(args, deep)
}

fn doit_renamed(args: &[Expr], deep: bool) -> Expr {
    let all = limits(args);
    let mut f = args[0].clone();
    for (index, limit) in all.iter().enumerate() {
        let (i, mut a, mut b) = (&limit[0], limit[1].clone(), limit[2].clone());
        let dif = sub(&b, &a);
        if crate::assumptions::is_integer(&dif) == Some(true) && crate::assumptions::is_negative(&dif) == Some(true) {
            (a, b) = (add2(&b, &Expr::one()), sub(&a, &Expr::one()));
            f = inv(&f);
        }
        match eval_product(&f, i, &a, &b) {
            Some(g) if !g.is_constant(crate::expr::Constant::NaN) => f = g,
            _ => {
                return held(Func::Product, &powsimp(&f), all[index..].to_vec())
                    .expect("valid product limits");
            }
        }
    }
    if deep { crate::functions::doit(&f, true) } else { powsimp(&f) }
}

/// A Dummy index with the assumptions both endpoints share, or None when
/// the index already has all of them.
fn inherited_dummy(x: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    use crate::assumptions as q;
    use crate::symbol::Assumptions as A;
    let symbol = x.as_sym()?;
    let checks: [(A, fn(&Expr) -> Option<bool>); 9] = [
        (A::NONNEGATIVE, q::is_nonnegative), (A::NONPOSITIVE, q::is_nonpositive),
        (A::POSITIVE, q::is_positive), (A::NEGATIVE, q::is_negative),
        (A::INTEGER, q::is_integer), (A::RATIONAL, q::is_rational),
        (A::FINITE, q::is_finite), (A::ZERO, q::is_zero), (A::REAL, q::is_real),
    ];
    let mut keep = A::NONE;
    let mut added = false;
    for (flag, query) in checks {
        if query(x) == Some(true) {
            keep = keep.with(flag);
        } else if query(a) == Some(true) && query(b) == Some(true) {
            keep = keep.with(flag);
            added = true;
        }
    }
    if !added {
        return None;
    }
    let facts = crate::assumption_facts::FactSet::from_legacy(keep);
    let dummy = crate::symbol::Symbol::dummy_with_facts(symbol.name(), facts, None).ok()?;
    Some(Expr::from_kind(Kind::Sym(dummy)))
}

fn direct(term: &Expr, k: &Expr, a: &Expr, count: i64) -> Expr {
    mul(&(0..count.max(0)).map(|i| subs1(term, k, &add2(a, &Expr::int(i)))).collect::<Vec<_>>())
}

/// `Product._eval_product(term, (k, a, n))`.
pub fn eval_product(term: &Expr, k: &Expr, a: &Expr, n: &Expr) -> Option<Expr> {
    if !has(term, k) {
        if crate::assumptions::is_zero(&sub(term, &Expr::one())) == Some(true) {
            return Some(Expr::one());
        }
        let count = add(&[n.clone(), crate::ops::neg(a), Expr::one()]);
        let infinite = |c: crate::expr::Constant| count.is_constant(c);
        let (up, down) = (infinite(crate::expr::Constant::Infinity), infinite(crate::expr::Constant::NegInfinity));
        if let (true, Some(c)) = (up || down, term.as_number()) {
            // a positive rational to the power oo, as SymPy's Rational does;
            // a -oo count comes from Karr's reversal, which inverts the term
            if c.is_rational() && !c.is_negative() && !c.is_zero() && !c.is_one() {
                let above = c.cmp_value(&Number::int(1)) == Some(std::cmp::Ordering::Greater);
                return Some(if above == up { Expr::constant(crate::expr::Constant::Infinity) } else { Expr::zero() });
            }
        }
        return Some(power(term, &count));
    }
    if a == n {
        return Some(subs1(term, k, a));
    }
    let count = add(&[n.clone(), crate::ops::neg(a), Expr::one()]);
    let dif = sub(n, a);
    let definite = dif.as_number().is_some_and(Number::is_integer);
    let small = dif.as_i64().filter(|d| *d < DIRECT);
    if definite && small.is_some() {
        return Some(direct(term, k, a, small.unwrap() + 1));
    }
    let fallback = |_: ()| if definite { dif.as_i64().map(|d| direct(term, k, a, d + 1)) } else { None };
    if crate::inspection::classify(term, std::slice::from_ref(k), crate::inspection::PolynomialKind::Polynomial) == Some(true) {
        let poly = crate::poly::polyobj::Poly::from_expr(term, std::slice::from_ref(k)).ok()?;
        let degree = poly.degree(0)? as u64;
        let roots = crate::poly::roots::roots_poly(&poly).ok()?;
        let mut factors = Vec::new();
        let mut found = 0u64;
        for (root, multiplicity) in roots {
            found += multiplicity as u64;
            factors.push(pow(&rf(&sub(a, &root), &count), &Expr::int(multiplicity as i64)));
        }
        if found < degree {
            // SymPy would recurse on the quotient by a k-free factor here
            // (and fail); keep the product unevaluated instead.
            return fallback(());
        }
        factors.push(pow(&poly.lc(), &count));
        return Some(mul(&factors));
    }
    match term.kind() {
        Kind::Add(..) => {
            let factored = factor_terms(term, true);
            if factored.is_mul() {
                return eval_product(&factored, k, a, n);
            }
        }
        Kind::Mul(..) => {
            let (with_k, without_k): (Vec<Expr>, Vec<Expr>) = term.args().into_iter().partition(|f| has(f, k));
            let without_k = power(&mul(&without_k), &count);
            if with_k.len() >= 2 {
                let (mut exclude, mut include) = (Vec::new(), Vec::new());
                for t in with_k {
                    match eval_product(&t, k, a, n) {
                        Some(p) => exclude.push(p),
                        None => include.push(t),
                    }
                }
                if exclude.is_empty() {
                    return None;
                }
                let rest = deep_doit(&mul(&include), k, a, n);
                return Some(mul(&[without_k, mul(&exclude), rest]));
            }
            let single = &with_k[0];
            let p = eval_product(single, k, a, n).unwrap_or_else(|| deep_doit(single, k, a, n));
            return Some(mul(&[without_k, p]));
        }
        Kind::Pow(base, exponent) => {
            if !has(base, k) {
                return Some(pow(base, &crate::summation::summation(exponent, k, a, n)));
            }
            if !has(exponent, k) {
                if let Some(p) = eval_product(base, k, a, n) {
                    return Some(pow(&p, exponent));
                }
            }
        }
        Kind::Func(Func::Product, args) => {
            let evaluated = doit(args, true);
            return Some(eval_product(&evaluated, k, a, n).unwrap_or_else(|| product_unevaluated(&evaluated, k, a, n)));
        }
        _ => {}
    }
    fallback(())
}

/// `Product(f, (k, a, n)).doit()`.
fn deep_doit(f: &Expr, k: &Expr, a: &Expr, n: &Expr) -> Expr {
    let e = product_unevaluated(f, k, a, n);
    let Kind::Func(_, args) = e.kind() else { unreachable!() };
    doit(args, true)
}

/// A term as a numeric coefficient and `(base, exponent)` factors.
fn term_factors(t: &Expr) -> (Number, Vec<(Expr, Expr)>) {
    let (c, rest) = crate::ops::as_coeff_mul(t);
    let factors = match rest.kind() {
        Kind::Mul(..) => rest.args(),
        _ if rest.is_one() => vec![],
        _ => vec![rest.clone()],
    };
    (c, factors.iter().map(crate::ops::as_base_exp).collect())
}

/// The smaller of two exponents of a common base, if comparable.
fn min_exponent(x: &Expr, y: &Expr) -> Option<Expr> {
    if x == y {
        return Some(x.clone());
    }
    let n = sub(x, y).as_number()?.clone();
    Some(if n.is_negative() { x.clone() } else { y.clone() })
}

type Factors = Vec<(Expr, Expr)>;

fn product_of(factors: &[(Expr, Expr)]) -> Expr {
    mul(&factors.iter().map(|(b, x)| pow(b, x)).collect::<Vec<_>>())
}

/// Common factors of several factor lists (minimal comparable exponents).
fn common_factors(lists: &[&Factors]) -> Factors {
    let mut out = lists[0].clone();
    for other in &lists[1..] {
        out = out.into_iter().filter_map(|(b, x)| {
            let (_, y) = other.iter().find(|(c, _)| *c == b)?;
            min_exponent(&x, y).map(|m| (b, m))
        }).collect();
    }
    out
}

/// `factors / by` for factor lists where `by` divides `factors`.
fn divide_factors(factors: &[(Expr, Expr)], by: &[(Expr, Expr)]) -> Factors {
    factors.iter().filter_map(|(b, x)| {
        let rest = match by.iter().find(|(c, _)| c == b) { Some((_, y)) => sub(x, y), None => x.clone() };
        (!rest.is_zero()).then(|| (b.clone(), rest))
    }).collect()
}

/// SymPy's `factor_terms(e, fraction=...)` on a single Add: the rational
/// content and the power factors common to every term come out (a common
/// minus sign too when every term is negative); with `fraction` the rest
/// is put over the least common denominator of its terms.
fn factor_terms(e: &Expr, fraction: bool) -> Expr {
    if !matches!(e.kind(), Kind::Add(..)) {
        return e.clone();
    }
    let terms: Vec<(Number, Factors)> = e.args().iter().map(term_factors).collect();
    let mut content = Number::int(1);
    if terms.iter().all(|(c, _)| c.is_rational()) {
        content = terms[0].0.abs();
        for (c, _) in &terms[1..] {
            content = content.gcd(c);
        }
        if content.is_zero() {
            content = Number::int(1);
        }
    }
    if terms.iter().all(|(c, _)| c.is_negative()) {
        content = content.neg();
    }
    // numerator (positive exponent) and denominator factors of each term
    let parts: Vec<(Number, Factors, Factors)> = terms.into_iter().map(|(c, factors)| {
        let (mut numer, mut denom) = (Vec::new(), Vec::new());
        for (b, x) in factors {
            if crate::ops::could_extract_minus_sign(&x) { denom.push((b, crate::ops::neg(&x))); } else { numer.push((b, x)); }
        }
        (c, numer, denom)
    }).collect();
    let gcd_numer = common_factors(&parts.iter().map(|p| &p.1).collect::<Vec<_>>());
    let gcd_denom = common_factors(&parts.iter().map(|p| &p.2).collect::<Vec<_>>());
    let reduced: Vec<(Number, Factors, Factors)> = parts.iter()
        .map(|(c, num, den)| (c.clone(), divide_factors(num, &gcd_numer), divide_factors(den, &gcd_denom)))
        .collect();
    // least common denominator: the largest exponent of every base
    let mut lcd: Factors = Vec::new();
    if fraction {
        for (_, _, den) in &reduced {
            for (b, x) in den {
                match lcd.iter_mut().find(|(c, _)| c == b) {
                    Some(entry) => match min_exponent(&entry.1, x) {
                        Some(m) => if m == entry.1 { entry.1 = x.clone() },
                        None => return e.clone(),
                    },
                    None => lcd.push((b.clone(), x.clone())),
                }
            }
        }
    }
    let inner = add(&reduced.iter().map(|(c, num, den)| {
        let scale = Expr::num(c.div(&content).unwrap_or_else(|| c.clone()));
        let rest = if fraction { product_of(&divide_factors(&lcd, den)) } else { inv(&product_of(den)) };
        mul(&[scale, product_of(num), rest])
    }).collect::<Vec<_>>());
    let outside = mul(&[product_of(&gcd_numer), inv(&product_of(&gcd_denom)), inv(&product_of(&lcd))]);
    crate::poly::api::keep_coeff(&Expr::num(content), &mul(&[outside, inner]))
}

fn coeff_is_negative(e: &Expr) -> bool {
    crate::ops::as_coeff_mul(e).0.is_negative()
}

/// SymPy's `as_base_exp`: `exp(a)` is `E**a` and `(1/q)**x` is `q**(-x)`.
fn sympy_base_exp(t: &Expr) -> (Expr, Expr) {
    match t.kind() {
        Kind::Func(Func::Exp, a) => (Expr::constant(crate::expr::Constant::E), a[0].clone()),
        Kind::Pow(b, x) => match b.as_number().and_then(Number::numer_denom) {
            Some((p, q)) if p == 1.into() && q != 1.into() => (Expr::num(Number::from_bigint(q)), crate::ops::neg(x)),
            _ => (b.clone(), x.clone()),
        },
        _ => (t.clone(), Expr::one()),
    }
}

/// `multiplicity(b, c)` for an integer `b > 1` and a rational `c`: the
/// exponent of `b` in `c` (negative when `b` divides the denominator).
fn rational_multiplicity(b: &num_bigint::BigInt, c: &Number) -> i64 {
    let count = |mut n: num_bigint::BigInt| {
        let mut m = 0;
        while !n.is_zero() && (&n % b).is_zero() {
            n /= b;
            m += 1;
        }
        m
    };
    match c.numer_denom() {
        Some((p, q)) => count(p.abs()) - count(q),
        None => 0,
    }
}

/// SymPy's `powsimp(e)` with its defaults (`deep=False`, `combine='all'`)
/// for the shapes product evaluation produces: powers of a common base
/// combine (`2**k/4 -> 2**(k - 2)`, `exp(a)*exp(b) -> exp(a + b)`), sum
/// exponents of positive bases are factored (`2**(n*(n + 1)/2)`), and
/// bases sharing an exponent are joined when that is valid (`6**x`).
fn powsimp(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Add(..) => add(&e.args().iter().map(powsimp).collect::<Vec<_>>()),
        Kind::Mul(..) => {
            let factors: Vec<Expr> = e.args().iter().map(|f| if f.is_number() { f.clone() } else { powsimp(f) }).collect();
            combine_bases(&combine_exponents(&mul(&factors)))
        }
        Kind::Pow(..) => combine_bases(&combine_exponents(e)),
        _ => e.clone(),
    }
}

/// The `combine='exp'` half of powsimp.
fn combine_exponents(e: &Expr) -> Expr {
    let (mut coeff, rest) = crate::ops::as_coeff_mul(e);
    let factors = if rest.is_mul() { rest.args() } else if rest.is_one() { vec![] } else { vec![rest] };
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    for t in factors {
        let (mut b, mut x) = sympy_base_exp(&t);
        if let Kind::Pow(base, _) = b.kind() {
            if let Some(n) = base.as_number().and_then(Number::to_bigint) {
                // an integer perfect power base: 4**x -> 2**(2*x)
                if let Some((root, m)) = crate::ntheory::perfect_power(&n) {
                    b = Expr::num(Number::from_bigint(root));
                    x = mul(&[Expr::int(m as i64), x]);
                }
            }
        }
        if let Some(n) = b.as_number().and_then(Number::to_bigint) {
            if let Some((root, m)) = crate::ntheory::perfect_power(&n) {
                b = Expr::num(Number::from_bigint(root));
                x = mul(&[Expr::int(m as i64), x]);
            }
        }
        if matches!(b.kind(), Kind::Pow(..) | Kind::Func(Func::Exp, _)) {
            (b, x) = (pow(&b, &x), Expr::one());
        }
        match groups.iter_mut().find(|(g, _)| *g == b) {
            Some((_, xs)) => xs.push(x),
            None => groups.push((b, vec![x])),
        }
    }
    for (b, xs) in &mut groups {
        // 2**k/4 -> 2**(k - 2)
        let Some(base) = b.as_number().and_then(Number::to_bigint) else { continue };
        if coeff.is_one() || base.abs() <= 1.into() || xs.iter().all(Expr::is_number) || !coeff.is_rational() {
            continue;
        }
        let m = rational_multiplicity(&base.abs(), &coeff);
        if m != 0 {
            xs.push(Expr::int(m));
            if let Some(scaled) = b.as_number().unwrap().pow_i64(m).and_then(|p| coeff.div(&p)) {
                coeff = scaled;
            }
        }
    }
    let mut out = vec![Expr::num(coeff)];
    out.extend(groups.into_iter().map(|(b, xs)| pow(&b, &add(&xs))));
    mul(&out)
}

/// The `combine='base'` half of powsimp.
fn combine_bases(e: &Expr) -> Expr {
    use crate::assumptions as q;
    let (coeff, rest) = crate::ops::as_coeff_mul(e);
    let factors = if rest.is_mul() { rest.args() } else if rest.is_one() { vec![] } else { vec![rest] };
    let mut by_exponent: Vec<(Expr, Vec<Expr>)> = Vec::new();
    for t in factors {
        let (mut b, mut x) = sympy_base_exp(&t);
        if b.is_number() && b.as_number().is_some_and(Number::is_rational) {
            let nonnegative = !b.as_number().unwrap().is_negative();
            if nonnegative || q::is_integer(&x) == Some(true) {
                // 2**(2*x) -> 4**x
                let (c, term) = crate::ops::as_coeff_mul(&x);
                if !c.is_one() && !term.is_one() && c.is_rational() {
                    (b, x) = (pow(&b, &Expr::num(c)), term);
                }
            }
        }
        if matches!(x.kind(), Kind::Add(..)) && (q::is_positive(&b) == Some(true) || q::is_integer(&x) == Some(true)) {
            let factored = factor_terms(&x, false);
            // a lone sum times a number stays expanded
            let lone_sum = matches!(crate::ops::as_coeff_mul(&factored).1.kind(), Kind::Add(..));
            if !lone_sum {
                x = factored;
                if coeff_is_negative(&x) {
                    (b, x) = (inv(&b), crate::ops::neg(&x));
                }
            }
        }
        match by_exponent.iter_mut().find(|(y, _)| *y == x) {
            Some((_, bases)) => bases.push(b),
            None => by_exponent.push((x, vec![b])),
        }
    }
    let mut out = vec![Expr::num(coeff)];
    for (x, bases) in by_exponent {
        let joinable = bases.len() == 1 || q::is_integer(&x) == Some(true)
            || bases.iter().filter(|b| q::is_nonnegative(b) != Some(true)).count() <= 1;
        if joinable {
            out.push(pow(&mul(&bases), &x));
        } else {
            let (nonneg, other): (Vec<Expr>, Vec<Expr>) = bases.into_iter().partition(|b| q::is_nonnegative(b) == Some(true));
            out.push(pow(&mul(&nonneg), &x));
            out.extend(other.iter().map(|b| pow(b, &x)));
        }
    }
    // SymPy's Mul writes (1/q)**x as q**(-x)
    let out: Vec<Expr> = out.into_iter().map(|f| match f.kind() {
        Kind::Pow(b, x) if !x.is_number() => match b.as_number().and_then(Number::numer_denom) {
            Some((p, q)) if p == 1.into() && q != 1.into() => pow(&Expr::num(Number::from_bigint(q)), &crate::ops::neg(x)),
            _ => f,
        },
        _ => f,
    }).collect();
    mul(&out)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(f: &str, a: &str, b: &str) -> String {
        let k = Expr::symbol("k");
        product(&parse_str(f).unwrap(), vec![vec![k, parse_str(a).unwrap(), parse_str(b).unwrap()]]).unwrap().to_string()
    }

    #[test]
    fn products_match_sympy() {
        for (f, a, b, expected) in [
            ("k", "1", "n", "factorial(n)"),
            ("k**2", "1", "n", "factorial(n)**2"),
            ("2", "1", "n", "2**n"),
            ("k", "1", "5", "120"),
            ("k+1", "1", "n", "RisingFactorial(2, n)"),
            ("2*k", "1", "n", "2**n*factorial(n)"),
            ("1/k", "1", "n", "1/factorial(n)"),
            ("k", "a", "b", "RisingFactorial(a, -a + b + 1)"),
            ("1+1/k", "1", "n", "RisingFactorial(2, n)/factorial(n)"),
            ("k", "3", "1", "1/2"),
            ("k", "1", "0", "1"),
            ("k", "n", "n", "n"),
            ("cos(k)", "1", "n", "Product(cos(k), (k, 1, n))"),
            ("exp(k)", "1", "n", "Product(exp(k), (k, 1, n))"),
            ("-k", "1", "n", "(-1)**n*factorial(n)"),
            ("x**k", "1", "n", "x**(n**2/2 + n/2)"),
        ] {
            assert_eq!(p(f, a, b), expected, "{f} ({a}, {b})");
        }
    }
}
