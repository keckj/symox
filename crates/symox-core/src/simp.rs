//! Targeted simplification and expansion routines: `trigsimp`,
//! `expand_trig`, `powsimp`, `expand_log`, `logcombine`, `expand_mul`,
//! `expand_power_base`, `expand_power_exp`, `nsimplify`.

use rustc_hash::FxHashMap;

use crate::assumptions::{is_positive, is_real};
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::*;
use crate::number::Number;
use crate::ops::{add, add2, div, mul, mul2, mul_factor, neg, pow, sqrt, sub};
use crate::parser::parse_str;
use crate::pattern::replace_with;
use crate::subs::xreplace;
use crate::traverse::{map_children, node_count};

// ----------------------------------------------------------------------
// expand_trig
// ----------------------------------------------------------------------

/// Chebyshev polynomial T_n(t) (first kind).
fn cheb_t(n: u32, t: &Expr) -> Expr {
    let (mut a, mut b) = (Expr::one(), t.clone());
    if n == 0 {
        return a;
    }
    for _ in 1..n {
        let c = expand(&sub(&mul(&[Expr::int(2), t.clone(), b.clone()]), &a));
        a = b;
        b = c;
    }
    b
}

/// Chebyshev polynomial U_n(t) (second kind).
fn cheb_u(n: u32, t: &Expr) -> Expr {
    let (mut a, mut b) = (Expr::one(), mul2(&Expr::int(2), t));
    if n == 0 {
        return a;
    }
    for _ in 1..n {
        let c = expand(&sub(&mul(&[Expr::int(2), t.clone(), b.clone()]), &a));
        a = b;
        b = c;
    }
    b
}

/// Expand trigonometric functions of sums and integer multiples.
pub fn expand_trig(e: &Expr) -> Expr { expand_trig_deep(e,true) }

pub fn expand_trig_deep(e: &Expr,deep:bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut expand_trig) };
    let Kind::Func(f, args) = e2.kind() else { return e2 };
    if !matches!(f, Func::Sin | Func::Cos | Func::Tan | Func::Sinh | Func::Cosh) {
        return e2;
    }
    let u = &args[0];
    match u.kind() {
        Kind::Add(..) => {
            let terms = crate::sortkey::ordered_terms(u);
            let (a, rest) = (terms[0].clone(), add(&terms[1..]));
            let (sa, ca, sb, cb) = (
                expand_trig(&sin(&a)),
                expand_trig(&cos(&a)),
                expand_trig(&sin(&rest)),
                expand_trig(&cos(&rest)),
            );
            let r = match f {
                Func::Sin => add2(&mul2(&sa, &cb), &mul2(&ca, &sb)),
                Func::Cos => sub(&mul2(&ca, &cb), &mul2(&sa, &sb)),
                Func::Tan => {
                    let (ta, tb) = (expand_trig(&tan(&a)), expand_trig(&tan(&rest)));
                    div(&add2(&ta, &tb), &sub(&Expr::one(), &mul2(&ta, &tb)))
                }
                Func::Sinh => {
                    let (sa, ca, sb, cb) = (sinh(&a), cosh(&a), expand_trig(&sinh(&rest)), expand_trig(&cosh(&rest)));
                    add2(&mul2(&sa, &cb), &mul2(&ca, &sb))
                }
                _ => {
                    let (sa, ca, sb, cb) = (sinh(&a), cosh(&a), expand_trig(&sinh(&rest)), expand_trig(&cosh(&rest)));
                    add2(&mul2(&ca, &cb), &mul2(&sa, &sb))
                }
            };
            if matches!(f,Func::Sinh|Func::Cosh){expand(&r)}else{r}
        }
        Kind::Mul(c, facs) if c.is_integer() && c.is_positive() && !c.is_one() => {
            let n = c.to_i64().unwrap_or(1);
            if n > 64 {
                return e2;
            }
            let n = n as u32;
            let x = crate::ops::add_term(&mul_from(facs), &Number::ONE);
            match f {
                Func::Sin => {
                    let s = sin(&x);
                    if n % 2 == 1 {
                        let sign = if ((n - 1) / 2) % 2 == 0 { 1 } else { -1 };
                        expand(&mul2(&Expr::int(sign), &cheb_t(n, &s)))
                    } else {
                        let sign = if (n / 2 - 1) % 2 == 0 { 1 } else { -1 };
                        expand(&mul(&[Expr::int(sign), cos(&x), cheb_u(n - 1, &s)]))
                    }
                }
                Func::Cos => expand(&cheb_t(n, &cos(&x))),
                Func::Tan => {
                    // tan(n x) through tan(a + b) recursion
                    let t = tan(&x);
                    let mut r = t.clone();
                    for _ in 1..n {
                        r = div(&add2(&r, &t), &sub(&Expr::one(), &mul2(&r, &t)));
                    }
                    crate::simplify::together(&r)
                }
                _ => e2,
            }
        }
        _ => e2,
    }
}

fn mul_from(facs: &[(Expr, Expr)]) -> Expr {
    mul(&facs.iter().map(|(b, x)| mul_factor(b, x)).collect::<Vec<_>>())
}

// ----------------------------------------------------------------------
// trigsimp
// ----------------------------------------------------------------------

struct Rule {
    pat: Expr,
    rep: Expr,
}

fn wilds() -> FxHashMap<&'static str, Expr> {
    let mut m = FxHashMap::default();
    for n in ["a", "b", "c", "d", "n"] {
        m.insert(n, Expr::wild(n, vec![]));
    }
    m
}

fn rule(p: &str, r: &str, w: &FxHashMap<&'static str, Expr>) -> Rule {
    let sub_w = |s: &str| -> Expr {
        let e = parse_str(s).expect("rule parses");
        let map: FxHashMap<Expr, Expr> = w.iter().map(|(k, v)| (Expr::symbol(k), v.clone())).collect();
        xreplace(&e, &map)
    };
    Rule { pat: sub_w(p), rep: sub_w(r) }
}

fn trig_rules() -> &'static Vec<Rule> {
    use std::sync::OnceLock;
    static RULES: OnceLock<Vec<Rule>> = OnceLock::new();
    RULES.get_or_init(|| {
        let w = wilds();
        vec![
            // Pythagorean identities
            rule("a*sin(b)**2 + a*cos(b)**2 + c", "a + c", &w),
            rule("a*cosh(b)**2 - a*sinh(b)**2 + c", "a + c", &w),
            rule("a*tan(b)**2 + a + c", "a/cos(b)**2 + c", &w),
            // double angle
            rule("2*a*sin(b)*cos(b) + c", "a*sin(2*b) + c", &w),
            rule("a*cos(b)**2 - a*sin(b)**2 + c", "a*cos(2*b) + c", &w),
            rule("2*a*cos(b)**2 - a + c", "a*cos(2*b) + c", &w),
            rule("a - 2*a*sin(b)**2 + c", "a*cos(2*b) + c", &w),
            // complements
            rule("a - a*sin(b)**2 + c", "a*cos(b)**2 + c", &w),
            rule("a - a*cos(b)**2 + c", "a*sin(b)**2 + c", &w),
            rule("a*cosh(b)**2 - a + c", "a*sinh(b)**2 + c", &w),
            // sum formulas
            rule("a*sin(b)*cos(d) + a*sin(d)*cos(b) + c", "a*sin(b + d) + c", &w),
            rule("a*cos(b)*cos(d) - a*sin(b)*sin(d) + c", "a*cos(b + d) + c", &w),
            rule("a*cos(b)*cos(d) + a*sin(b)*sin(d) + c", "a*cos(b - d) + c", &w),
            // quotients
            rule("a*sin(b)**n*cos(b)**(-n)", "a*tan(b)**n", &w),
            rule("a*sinh(b)**n*cosh(b)**(-n)", "a*tanh(b)**n", &w),
        ]
    })
}

fn apply_rules(e: &Expr) -> Expr {
    let mut cur = e.clone();
    for _ in 0..16 {
        let mut changed = false;
        for r in trig_rules() {
            let next = replace_with(&cur, &r.pat, &mut |bind, orig| {
                let v = xreplace(&r.rep, bind);
                if node_count(&v) < node_count(orig) || (node_count(&v) == node_count(orig) && v != *orig) {
                    Some(v)
                } else {
                    None
                }
            });
            if next != cur {
                cur = next;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    cur
}

/// sin(u)**2 -> 1 - cos(u)**2 (or the reverse) everywhere.
fn rewrite_squares(e: &Expr, to_cos: bool) -> Expr {
    let e2 = if e.is_atom() { e.clone() } else { map_children(e, &mut |c| rewrite_squares(c, to_cos)) };
    if let Kind::Pow(b, x) = e2.kind() {
        if let (Kind::Func(f, a), Some(k)) = (b.kind(), x.as_i64()) {
            if k >= 2 && k % 2 == 0 {
                let target = if to_cos { Func::Sin } else { Func::Cos };
                if *f == target {
                    let other = if to_cos { cos(&a[0]) } else { sin(&a[0]) };
                    return pow(&sub(&Expr::one(), &other.powi(2)), &Expr::int(k / 2));
                }
            }
        }
    }
    e2
}

fn has_trig(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc | Func::Sinh | Func::Cosh | Func::Tanh, _) => true,
        _ => e.args().iter().any(has_trig),
    }
}

/// Simplify trigonometric expressions using identities.
pub fn trigsimp(e: &Expr) -> Expr {
    if !has_trig(e) {
        return e.clone();
    }
    // simplify inside first
    let e1 = if e.is_atom() { e.clone() } else { map_children(e, &mut trigsimp) };
    let mut best = apply_rules(&e1);
    let mut cands = vec![];
    let ex = expand(&e1);
    if ex != e1 {
        cands.push(apply_rules(&ex));
    }
    for to_cos in [true, false] {
        let r = expand(&rewrite_squares(&e1, to_cos));
        cands.push(apply_rules(&r));
    }
    for c in cands {
        if node_count(&c) < node_count(&best) {
            best = c;
        }
    }
    best
}

// ----------------------------------------------------------------------
// powers and logarithms
// ----------------------------------------------------------------------

/// Combine powers with equal bases (and exponentials), SymPy's `powsimp`.
pub fn powsimp(e: &Expr) -> Expr {
    powsimp_deep(e, true)
}

pub fn powsimp_deep(e: &Expr, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| powsimp_deep(c, deep)) };
    let Kind::Mul(c, facs) = e2.kind() else { return e2 };
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    let mut exps: Vec<Expr> = Vec::new();
    for (b, x) in facs {
        if let Kind::Func(Func::Exp, a) = b.kind() {
            exps.push(mul2(&a[0], x));
            continue;
        }
        if b.is_constant(Constant::E) {
            exps.push(x.clone());
            continue;
        }
        match groups.iter_mut().find(|(g, _)| g == b) {
            Some((_, v)) => v.push(x.clone()),
            None => groups.push((b.clone(), vec![x.clone()])),
        }
    }
    let mut out = vec![Expr::num(c.clone())];
    for (b, xs) in groups {
        out.push(Expr::from_kind(Kind::Pow(b.clone(), add(&xs))).normalize_pow());
    }
    if !exps.is_empty() {
        out.push(exp(&add(&exps)));
    }
    mul(&out)
}

impl Expr {
    /// Re-canonicalise a raw `Pow` node.
    fn normalize_pow(self) -> Expr {
        match self.kind() {
            Kind::Pow(b, x) => pow(b, x),
            _ => self,
        }
    }
}

/// `log(x*y) -> log(x) + log(y)`, `log(x**n) -> n*log(x)` for positive
/// arguments (or everything with `force`).
pub fn expand_log(e: &Expr, force: bool) -> Expr {
    expand_log_deep(e, force, true)
}

pub fn expand_log_deep(e: &Expr, force: bool, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| expand_log_deep(c, force, deep)) };
    let Kind::Func(Func::Log, a) = e2.kind() else { return e2 };
    let ok = |x: &Expr| force || is_positive(x) == Some(true);
    match a[0].kind() {
        Kind::Mul(c, facs) => {
            let mut parts = Vec::new();
            let mut rest = vec![Expr::num(c.clone())];
            for (b, x) in facs {
                let f = mul_factor(b, x);
                if ok(b) && is_real(x) == Some(true) || force {
                    parts.push(expand_log(&log(&f), force));
                } else {
                    rest.push(f);
                }
            }
            if parts.is_empty() {
                return e2;
            }
            let r = mul(&rest);
            if !r.is_one() {
                parts.push(log(&r));
            }
            add(&parts)
        }
        Kind::Pow(b, x) if (ok(b) && is_real(x) == Some(true)) || force => mul2(x, &expand_log(&log(b), force)),
        _ => e2,
    }
}

/// Combine `a*log(x) + b*log(y)` into `log(x**a * y**b)` for positive
/// arguments and real coefficients (or everything with `force`).
pub fn logcombine(e: &Expr, force: bool) -> Expr {
    let e2 = if e.is_atom() { e.clone() } else { map_children(e, &mut |c| logcombine(c, force)) };
    let ok = |x: &Expr| force || is_positive(x) == Some(true);
    let as_coeff_log = |t: &Expr| -> Option<(Expr, Expr)> {
        match t.kind() {
            Kind::Func(Func::Log, a) => Some((Expr::one(), a[0].clone())),
            Kind::Mul(..) => {
                let fs = t.args();
                let i = fs.iter().position(|f| f.is_func(&Func::Log))?;
                let Kind::Func(_, a) = fs[i].kind() else { return None };
                let coeff = mul(&fs.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, f)| f.clone()).collect::<Vec<_>>());
                Some((coeff, a[0].clone()))
            }
            _ => None,
        }
    };
    match e2.kind() {
        Kind::Add(..) => {
            let mut logs = Vec::new();
            let mut rest = Vec::new();
            for t in e2.args() {
                match as_coeff_log(&t) {
                    Some((c, a)) if ok(&a) && (force || is_real(&c) == Some(true)) => logs.push(pow(&a, &c)),
                    _ => rest.push(t),
                }
            }
            if logs.len() < 2 && !(logs.len() == 1 && rest.is_empty()) {
                return e2;
            }
            rest.push(log(&mul(&logs)));
            add(&rest)
        }
        Kind::Mul(..) => match as_coeff_log(&e2) {
            Some((c, a)) if !c.is_one() && ok(&a) && (force || is_real(&c) == Some(true)) => log(&pow(&a, &c)),
            _ => e2,
        },
        _ => e2,
    }
}

/// `(x*y)**a -> x**a * y**a` (for positive factors, or with `force`).
pub fn expand_power_base(e: &Expr, force: bool) -> Expr {
    expand_power_base_deep(e, force, true)
}

pub fn expand_power_base_deep(e: &Expr, force: bool, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| expand_power_base_deep(c, force, deep)) };
    if let Kind::Pow(b, x) = e2.kind() {
        if let Kind::Mul(..) = b.kind() {
            let fs = b.args();
            if force || fs.iter().all(|f| is_positive(f) == Some(true)) {
                return mul(&fs.iter().map(|f| pow(f, x)).collect::<Vec<_>>());
            }
        }
    }
    e2
}

/// `x**(a+b) -> x**a * x**b` when the base is known to be non-zero.
pub fn expand_power_exp(e: &Expr) -> Expr {
    expand_power_exp_deep(e, true)
}

pub fn expand_power_exp_deep(e: &Expr, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else if let Kind::Add(c, terms) = e.kind() {
        // Scaling an expanded child with mul2 would merge its factors again.
        let transformed: Vec<_> = terms.iter().map(|(t,k)| (expand_power_exp_deep(t, true), k.clone())).collect();
        if transformed.iter().zip(terms).all(|((a,_),(b,_))| a == b) { e.clone() }
        else { Expr::from_unevaluated_kind(Kind::Add(c.clone(), transformed)) }
    } else { map_children(e, &mut |c| expand_power_exp_deep(c, deep)) };
    match e2.kind() {
        Kind::Pow(b, x) if x.is_add() && crate::assumptions::is_zero(b) == Some(false) => {
            expanded_product(&x.args().iter().map(|t| pow(b, t)).collect::<Vec<_>>())
        }
        Kind::Func(Func::Exp, a) if a[0].is_add() => expanded_product(&a[0].args().iter().map(exp).collect::<Vec<_>>()),
        Kind::Mul(..) if deep => {
            let args = e2.args();
            let changed: Vec<_> = args.iter().map(|f| expand_power_exp_deep(f, false)).collect();
            if args.iter().zip(&changed).any(|(a,b)| a != b) { expanded_product(&changed) } else { e2 }
        }
        _ => e2,
    }
}

/// The ordinary product constructor merges equal bases, undoing this
/// expansion. Preserve its explicit factors just as for evaluate=False.
pub(crate) fn expanded_product(args: &[Expr]) -> Expr {
    let mut coefficient = Number::ONE;
    let mut factors = Vec::new();
    for arg in args {
        match arg.kind() {
            Kind::Num(n) => coefficient = coefficient.mul(n),
            Kind::Mul(c, fs) => { coefficient = coefficient.mul(c); factors.extend(fs.iter().cloned()); },
            Kind::Pow(base, exponent) => factors.push((base.clone(), exponent.clone())),
            _ => factors.push((arg.clone(), Expr::one())),
        }
    }
    if factors.is_empty() { return Expr::num(coefficient); }
    Expr::from_unevaluated_kind(Kind::Mul(coefficient, factors))
}

/// Distribute products over sums (but do not expand powers).
pub fn expand_mul(e: &Expr) -> Expr {
    expand_mul_deep(e, true)
}

pub fn expand_mul_deep(e: &Expr, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| expand_mul_deep(c, deep)) };
    if let Kind::Mul(..) = e2.kind() {
        let fs = e2.args();
        if fs.iter().any(|f| f.is_add()) {
            let mut acc = vec![Expr::one()];
            for f in fs {
                let terms = if f.is_add() { f.args() } else { vec![f] };
                let mut next = Vec::with_capacity(acc.len() * terms.len());
                for a in &acc {
                    for t in &terms {
                        next.push(mul2(a, t));
                    }
                }
                acc = next;
            }
            return add(&acc);
        }
    }
    e2
}

// ----------------------------------------------------------------------
// nsimplify
// ----------------------------------------------------------------------

/// Best rational approximation of `x` with denominator <= `max_den`.
fn rational_approx(x: f64, max_den: i64, tol: f64) -> Option<(i64, i64)> {
    let (mut h0, mut h1, mut k0, mut k1) = (0i64, 1i64, 1i64, 0i64);
    let mut v = x;
    for _ in 0..64 {
        let a = v.floor();
        if a.abs() > 1e15 {
            return None;
        }
        let a = a as i64;
        let h2 = a.checked_mul(h1)?.checked_add(h0)?;
        let k2 = a.checked_mul(k1)?.checked_add(k0)?;
        if k2 > max_den {
            return None;
        }
        if ((h2 as f64) / (k2 as f64) - x).abs() <= tol * x.abs().max(1.0) {
            return Some((h2, k2));
        }
        (h0, h1, k0, k1) = (h1, h2, k1, k2);
        let frac = v - a as f64;
        if frac.abs() < 1e-300 {
            return None;
        }
        v = 1.0 / frac;
    }
    None
}

/// Find a simple exact expression close to a float (rationals, square
/// roots of rationals, rational multiples of the given constants).
pub fn nsimplify_float(x: f64, constants: &[Expr]) -> Option<Expr> {
    let tol = 1e-12;
    if x == 0.0 {
        return Some(Expr::zero());
    }
    // simplest candidates first: small rationals, multiples of the given
    // constants, square roots, then rationals with larger denominators
    if let Some((p, q)) = rational_approx(x, 10_000, tol) {
        return Some(Expr::rational(p, q));
    }
    for c in constants {
        let cv = crate::evalf::eval_f64(c)?;
        if let Some((p, q)) = rational_approx(x / cv, 10_000, tol) {
            return Some(mul2(&Expr::rational(p, q), c));
        }
    }
    let sq = x * x;
    if let Some((p, q)) = rational_approx(sq, 10_000, 1e-11) {
        let r = sqrt(&Expr::rational(p, q));
        return Some(if x < 0.0 { neg(&r) } else { r });
    }
    if let Some((p, q)) = rational_approx(x, 1_000_000_000, tol) {
        return Some(Expr::rational(p, q));
    }
    None
}

/// Replace floats in `e` by simple exact numbers.
pub fn nsimplify(e: &Expr, constants: &[Expr]) -> Expr {
    match e.kind() {
        Kind::Num(Number::Float(f)) => nsimplify_float(*f, constants).unwrap_or_else(|| e.clone()),
        _ if e.is_atom() => e.clone(),
        _ => {
            if crate::functions::is_numeric_value(e) {
                if let Some(v) = crate::evalf::eval_f64(e) {
                    if let Some(r) = nsimplify_float(v, constants) {
                        return r;
                    }
                }
            }
            match e.kind() {
                // numeric coefficients are not children: go through args()
                Kind::Add(..) => add(&e.args().iter().map(|c| nsimplify(c, constants)).collect::<Vec<_>>()),
                Kind::Mul(..) => mul(&e.args().iter().map(|c| nsimplify(c, constants)).collect::<Vec<_>>()),
                _ => map_children(e, &mut |c| nsimplify(c, constants)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(f: fn(&Expr) -> Expr, s: &str) -> String {
        f(&parse_str(s).unwrap()).to_string()
    }

    #[test]
    fn trig() {
        assert_eq!(t(expand_trig, "sin(2*x)"), "2*sin(x)*cos(x)");
        assert_eq!(t(expand_trig, "cos(2*x)"), "2*cos(x)**2 - 1");
        assert_eq!(t(expand_trig, "sin(x+y)"), "sin(x)*cos(y) + sin(y)*cos(x)");
        assert_eq!(t(expand_trig, "cos(x+y)"), "-sin(x)*sin(y) + cos(x)*cos(y)");
        assert_eq!(t(expand_trig, "sin(3*x)"), "-4*sin(x)**3 + 3*sin(x)");
        assert_eq!(t(expand_trig, "cos(3*x)"), "4*cos(x)**3 - 3*cos(x)");
        assert_eq!(t(trigsimp, "sin(x)**2+cos(x)**2"), "1");
        assert_eq!(t(trigsimp, "2*sin(x)**2+2*cos(x)**2"), "2");
        assert_eq!(t(trigsimp, "sin(x)/cos(x)"), "tan(x)");
        assert_eq!(t(trigsimp, "1-sin(x)**2"), "cos(x)**2");
        assert_eq!(t(trigsimp, "2*sin(x)*cos(x)"), "sin(2*x)");
        assert_eq!(t(trigsimp, "cos(x)**2-sin(x)**2"), "cos(2*x)");
        assert_eq!(t(trigsimp, "sin(x)**4 - cos(x)**4"), "-cos(2*x)");
        assert_eq!(t(trigsimp, "1+tan(x)**2"), "cos(x)**(-2)");
        assert_eq!(t(trigsimp, "sin(x)*cos(y)+cos(x)*sin(y)"), "sin(x + y)");
        assert_eq!(t(trigsimp, "cosh(x)**2-sinh(x)**2"), "1");
        assert_eq!(t(trigsimp, "x*sin(y)**2 + x*cos(y)**2"), "x");
    }

    #[test]
    fn powers_logs() {
        assert_eq!(t(powsimp, "x**a*x**b"), "x**(a + b)");
        assert_eq!(t(powsimp, "exp(x)*exp(y)"), "exp(x + y)");
        assert_eq!(t(powsimp, "2**x*2**y"), "2**(x + y)");
        assert_eq!(t(powsimp, "x**2*x**a"), "x**(a + 2)");
        assert_eq!(t(expand_mul, "x*(y+z)"), "x*y + x*z");
        assert_eq!(nsimplify(&Expr::float(0.5), &[]).to_string(), "1/2");
        assert_eq!(nsimplify(&Expr::float(1.0 / 3.0), &[]).to_string(), "1/3");
        assert_eq!(nsimplify(&Expr::float(2f64.sqrt()), &[]).to_string(), "sqrt(2)");
        assert_eq!(nsimplify(&Expr::float(std::f64::consts::PI), &[Expr::pi()]).to_string(), "pi");
    }
}

