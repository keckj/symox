//! The table of `_create_lookup_table`: rewrites of single factors as sums
//! of `C*x**s*G(eta*x**b)`.
//!
//! Each entry is recognised structurally the way `f.match(formula)` binds
//! the wild symbols `p`, `q` (with `t = p*x**q`), `a`, `b` and `n`; entries
//! are tried in the reference order of their bucket (`_mytype`), and the
//! bindings are polarified exactly like the reference before the Meijer G
//! parameters, the coefficient and the condition are built.
//!
//! Not ported (never reached by the corpus; such factors give `None`): the
//! entries `(t**a - b**a)/(t - b)` and the eight radical forms of sections
//! 12 and 13, and the recursive Mellin-transform rewrite.

use super::{get_coeff_exp, has_infinity, has_x, is_false, G};
use crate::ops::{add, div, mul, mul2, neg, pow, sub};
use crate::special::polar::{polarify, unpolarify};
use crate::special::SpecialFn;
use crate::{Expr, Func, Kind};

/// The bucket of `f` (`_mytype`): function heads, sorted by class key.
fn mytype(f: &Expr, x: &Expr, out: &mut Vec<Func>) {
    if !has_x(f, x) {
        return;
    }
    if let Kind::Func(h, _) = f.kind() {
        if super::is_sympy_function(f) {
            out.push(h.clone());
            return;
        }
    }
    if f.is_atom() {
        return;
    }
    for a in f.args() {
        mytype(&a, x, out);
    }
}

/// `p*x**q` with `p`, `q` free of `x`.
fn match_t(e: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    if e == x {
        return Some((Expr::one(), Expr::one()));
    }
    if let Kind::Pow(b, q) = e.kind() {
        return (b == x && !has_x(q, x)).then(|| (Expr::one(), q.clone()));
    }
    if !e.is_mul() {
        return None;
    }
    let mut coeff = Vec::new();
    let mut found = None;
    for f in crate::simptools::mul_args(e) {
        if !has_x(&f, x) {
            coeff.push(f);
            continue;
        }
        if found.is_some() {
            return None;
        }
        found = Some(match_t(&f, x).filter(|(p, _)| p.is_one())?.1);
    }
    Some((mul(&coeff), found?))
}

/// `p*x**q + c` with exactly one term depending on `x` and `c` nonzero.
fn match_add_t(e: &Expr, x: &Expr) -> Option<(Expr, Expr, Expr)> {
    if !e.is_add() {
        return None;
    }
    let mut consts = Vec::new();
    let mut dep = None;
    for t in crate::simptools::add_args(e) {
        if !has_x(&t, x) {
            consts.push(t);
        } else if dep.replace(t).is_some() {
            return None;
        }
    }
    let (p, q) = match_t(&dep?, x)?;
    let c = add(&consts);
    (!c.is_zero()).then_some((p, q, c))
}

/// `base**k`: the base and the exponent (one for a non-power).
fn base_exp(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Pow(b, k) => (b.clone(), k.clone()),
        _ => (e.clone(), Expr::one()),
    }
}

/// `unpolarify(polarify(v, lift=True), exponents_only=True)`.
fn pol(v: &Expr) -> Expr {
    unpolarify(&polarify(v, true), true)
}

fn g(an: Vec<Expr>, ao: Vec<Expr>, bm: Vec<Expr>, bo: Vec<Expr>, z: Expr) -> RawG {
    RawG { an, ao, bm, bo, z }
}

/// A table G-function before the constructor checks.
struct RawG {
    an: Vec<Expr>,
    ao: Vec<Expr>,
    bm: Vec<Expr>,
    bo: Vec<Expr>,
    z: Expr,
}

/// One instantiated table entry: `(fac, G)` terms, condition and hint.
struct Entry {
    terms: Vec<(Expr, RawG)>,
    cond: Expr,
    hint_false: bool,
}

fn h() -> Expr {
    Expr::rational(1, 2)
}
fn int(i: i64) -> Expr {
    Expr::int(i)
}
fn sqrt_pi() -> Expr {
    crate::ops::sqrt(&Expr::pi())
}
fn polar_minus_one() -> Expr {
    crate::special::apply(SpecialFn::PolarLift, vec![Expr::minus_one()])
}

fn constant(a: &Expr) -> Vec<(Expr, RawG)> {
    vec![
        (a.clone(), g(vec![int(1)], vec![], vec![], vec![int(0)], Expr::zero())),
        (a.clone(), g(vec![], vec![int(1)], vec![int(0)], vec![], Expr::zero())),
    ]
}

/// `constant(a)` with the variable argument filled in.
fn constant_x(a: &Expr, x: &Expr) -> Vec<(Expr, RawG)> {
    let mut v = constant(a);
    for (_, g) in &mut v {
        g.z = x.clone();
    }
    v
}

fn simple(fac: Expr, g: RawG) -> Entry {
    Entry { terms: vec![(fac, g)], cond: Expr::boolean(true), hint_false: false }
}

/// The candidate entries of the bucket of `f`, in table order, already
/// instantiated with their (polarified) bindings.
fn candidates(f: &Expr, x: &Expr) -> Vec<Entry> {
    let mut key = Vec::new();
    mytype(f, x, &mut key);
    let mut out = Vec::new();
    let tq = |p: &Expr, q: &Expr| mul2(p, &pow(x, q));
    match key.as_slice() {
        [] => {
            if !has_x(f, x) {
                out.push(Entry { terms: constant_x(&pol(f), x), cond: Expr::boolean(true), hint_false: false });
                return out;
            }
            // (b + t)**(-a)
            let (base, k) = base_exp(f);
            if !has_x(&k, x) {
                if let Some((p, q, b)) = match_add_t(&base, x) {
                    let (p, q, b, a) = (pol(&p), pol(&q), pol(&b), pol(&neg(&k)));
                    let ua = unpolarify(&a, false);
                    let hint_false = ua.is_integer_number() && ua.as_number().is_some_and(|n| !n.is_positive());
                    let t = tq(&p, &q);
                    out.push(Entry {
                        terms: vec![(
                            div(&pow(&b, &neg(&a)), &crate::functions::gamma(&a)),
                            g(vec![sub(&Expr::one(), &a)], vec![], vec![int(0)], vec![], div(&t, &b)),
                        )],
                        cond: Expr::boolean(true),
                        hint_false,
                    });
                }
            }
        }
        [Func::Heaviside] => heaviside_entries(f, x, &mut out),
        [Func::Abs] => {
            // Abs(t - b)**(-a)
            let (base, k) = base_exp(f);
            if let Kind::Func(Func::Abs, args) = base.kind() {
                if let (Some((p, q, c)), false) = (match_add_t(&args[0], x), has_x(&k, x)) {
                    let (p, q, b, a) = (pol(&p), pol(&q), pol(&neg(&c)), pol(&neg(&k)));
                    let t = tq(&p, &q);
                    let one_a = sub(&Expr::one(), &a);
                    let half = div(&one_a, &int(2));
                    let fac = mul(&[
                        int(2),
                        crate::functions::sin(&div(&mul2(&Expr::pi(), &a), &int(2))),
                        crate::functions::gamma(&one_a),
                        pow(&crate::functions::abs(&b), &neg(&a)),
                    ]);
                    out.push(Entry {
                        terms: vec![(fac, g(vec![one_a.clone()], vec![half.clone()], vec![int(0)], vec![half], div(&t, &b)))],
                        cond: crate::relational::lt(&crate::functions::re(&a), &Expr::one()),
                        hint_false: false,
                    });
                }
            }
        }
        [Func::Exp] => {
            // exp(polar_lift(-1)*t)
            if let Kind::Func(Func::Exp, args) = f.kind() {
                if let Some((p, q)) = match_t(&args[0], x) {
                    let p = mul2(&p, &pow(&polar_minus_one(), &Expr::minus_one()));
                    let (p, q) = (pol(&p), pol(&q));
                    out.push(simple(Expr::one(), g(vec![], vec![], vec![int(0)], vec![], tq(&p, &q))));
                }
            }
        }
        [h1] => {
            if let Some(e) = unary_entry(h1, f, x) {
                out.push(e);
            }
            if *h1 == Func::Log {
                log_entries(f, x, &mut out);
            }
        }
        [Func::Heaviside, Func::Log] | [Func::Log, Func::Heaviside] => log_heaviside_entries(f, x, &mut out),
        _ => {}
    }
    out
}

/// Entries for `F(t)` and `F(a, t)` with a single function head.
fn unary_entry(head: &Func, f: &Expr, x: &Expr) -> Option<Entry> {
    let Kind::Func(fh, args) = f.kind() else { return None };
    if fh != head {
        return None;
    }
    let one_arg = |n: usize| args.len() == n;
    let (order, arg) = match args.len() {
        1 => (None, &args[0]),
        2 if matches!(head, Func::Special(SpecialFn::Expint | SpecialFn::BesselJ | SpecialFn::BesselY | SpecialFn::BesselI | SpecialFn::BesselK)) => {
            if has_x(&args[0], x) {
                return None;
            }
            (Some(pol(&args[0])), &args[1])
        }
        _ => return None,
    };
    let _ = one_arg;
    let (p, q) = match_t(arg, x)?;
    let (p, q) = (pol(&p), pol(&q));
    let t = mul2(&p, &pow(x, &q));
    let t2 = div(&pow(&t, &int(2)), &int(4));
    let pi32 = pow(&Expr::pi(), &Expr::rational(3, 2));
    let e = match (head, order) {
        (Func::Sinh, None) => simple(pi32, g(vec![], vec![int(1)], vec![h()], vec![int(1), int(0)], t2)),
        (Func::Cosh, None) => simple(pi32, g(vec![], vec![h()], vec![int(0)], vec![h(), h()], t2)),
        (Func::Sin, None) => simple(sqrt_pi(), g(vec![], vec![], vec![h()], vec![int(0)], t2)),
        (Func::Cos, None) => simple(sqrt_pi(), g(vec![], vec![], vec![int(0)], vec![h()], t2)),
        (Func::Special(SpecialFn::Sinc), None) => simple(div(&sqrt_pi(), &int(2)), g(vec![], vec![], vec![int(0)], vec![Expr::rational(-1, 2)], t2)),
        (Func::Special(SpecialFn::Ei), None) => {
            let mut terms = constant_x(&mul2(&neg(&Expr::i()), &Expr::pi()), x);
            terms.push((Expr::minus_one(), g(vec![], vec![int(1)], vec![int(0), int(0)], vec![], mul2(&t, &polar_minus_one()))));
            Entry { terms, cond: Expr::boolean(true), hint_false: false }
        }
        (Func::Special(SpecialFn::Si), None) => simple(div(&sqrt_pi(), &int(2)), g(vec![int(1)], vec![], vec![h()], vec![int(0), int(0)], t2)),
        (Func::Special(SpecialFn::Ci), None) => simple(neg(&div(&sqrt_pi(), &int(2))), g(vec![], vec![int(1)], vec![int(0), int(0)], vec![h()], t2)),
        (Func::Special(SpecialFn::Shi), None) => simple(
            div(&mul2(&t, &sqrt_pi()), &int(4)),
            g(vec![h()], vec![], vec![int(0)], vec![Expr::rational(-1, 2), Expr::rational(-1, 2)], mul2(&polar_minus_one(), &t2)),
        ),
        (Func::Special(SpecialFn::Chi), None) => simple(
            neg(&div(&pi32, &int(2))),
            g(vec![], vec![h(), int(1)], vec![int(0), int(0)], vec![h(), h()], t2),
        ),
        (Func::Special(SpecialFn::Expint), Some(a)) => simple(Expr::one(), g(vec![], vec![a.clone()], vec![sub(&a, &Expr::one()), int(0)], vec![], t)),
        (Func::Erf, None) => simple(div(&Expr::one(), &sqrt_pi()), g(vec![int(1)], vec![], vec![h()], vec![int(0)], pow(&t, &int(2)))),
        (Func::Erfc, None) => simple(div(&Expr::one(), &sqrt_pi()), g(vec![], vec![int(1)], vec![int(0), h()], vec![], pow(&t, &int(2)))),
        (Func::Special(SpecialFn::Erfi), None) => simple(
            div(&t, &sqrt_pi()),
            g(vec![h()], vec![], vec![int(0)], vec![Expr::rational(-1, 2)], neg(&pow(&t, &int(2)))),
        ),
        (Func::Special(SpecialFn::FresnelS | SpecialFn::FresnelC), None) => {
            let (b0, b1) = if *head == Func::Special(SpecialFn::FresnelS) {
                (Expr::rational(3, 4), Expr::rational(1, 4))
            } else {
                (Expr::rational(1, 4), Expr::rational(3, 4))
            };
            let z = div(&mul2(&pow(&Expr::pi(), &int(2)), &pow(&t, &int(4))), &int(16));
            simple(h(), g(vec![int(1)], vec![], vec![b0], vec![int(0), b1], z))
        }
        (Func::Special(SpecialFn::BesselJ), Some(a)) => {
            simple(Expr::one(), g(vec![], vec![], vec![div(&a, &int(2))], vec![neg(&div(&a, &int(2)))], t2))
        }
        (Func::Special(SpecialFn::BesselY), Some(a)) => {
            let m = neg(&div(&add(&[a.clone(), Expr::one()]), &int(2)));
            simple(Expr::one(), g(vec![], vec![m.clone()], vec![div(&a, &int(2)), neg(&div(&a, &int(2)))], vec![m], t2))
        }
        (Func::Special(SpecialFn::BesselI), Some(a)) => {
            let m = div(&add(&[Expr::one(), a.clone()]), &int(2));
            simple(Expr::pi(), g(vec![], vec![m.clone()], vec![div(&a, &int(2))], vec![neg(&div(&a, &int(2))), m], t2))
        }
        (Func::Special(SpecialFn::BesselK), Some(a)) => {
            simple(h(), g(vec![], vec![], vec![div(&a, &int(2)), neg(&div(&a, &int(2)))], vec![], t2))
        }
        (Func::Special(SpecialFn::EllipticK), None) => simple(h(), g(vec![h(), h()], vec![], vec![int(0)], vec![int(0)], neg(&t))),
        (Func::Special(SpecialFn::EllipticE), None) => {
            simple(Expr::rational(-1, 4), g(vec![h(), Expr::rational(3, 2)], vec![], vec![int(0)], vec![int(0)], neg(&t)))
        }
        _ => return None,
    };
    Some(e)
}

/// `log(t)**n` (`n` a positive integer): the terms of `make_log1`/`make_log2`.
fn log_power(f: &Expr, x: &Expr) -> Option<(Expr, Expr, i64)> {
    let (base, k) = base_exp(f);
    let Kind::Func(Func::Log, args) = base.kind() else { return None };
    let n = k.as_i64().filter(|n| k.is_integer_number() && *n > 0)?;
    let (p, q) = match_t(&args[0], x)?;
    Some((pol(&p), pol(&q), n))
}

fn factorial(n: i64) -> Expr {
    Expr::num(crate::number::Number::from_bigint(crate::special::integers::factorial(n as u64)))
}

fn make_log1(t: &Expr, n: i64) -> (Expr, RawG) {
    let c = mul2(&pow(&Expr::minus_one(), &int(n)), &factorial(n));
    (c, g(vec![], vec![int(1); n as usize + 1], vec![int(0); n as usize + 1], vec![], t.clone()))
}
fn make_log2(t: &Expr, n: i64) -> (Expr, RawG) {
    (factorial(n), g(vec![int(1); n as usize + 1], vec![], vec![], vec![int(0); n as usize + 1], t.clone()))
}

fn log_entries(f: &Expr, x: &Expr, out: &mut Vec<Entry>) {
    // log(t)**n
    if let Some((p, q, n)) = log_power(f, x) {
        let t = mul2(&p, &pow(x, &q));
        out.push(Entry { terms: vec![make_log1(&t, n), make_log2(&t, n)], cond: Expr::boolean(true), hint_false: false });
    }
    let Kind::Func(Func::Log, args) = f.kind() else { return };
    // log(t + a)
    if let Some((p, q, a)) = match_add_t(&args[0], x) {
        let (p, q, a) = (pol(&p), pol(&q), pol(&a));
        let t = mul2(&p, &pow(x, &q));
        let mut terms = constant_x(&crate::functions::log(&a), x);
        terms.push((Expr::one(), g(vec![int(1), int(1)], vec![], vec![int(1)], vec![int(0)], div(&t, &a))));
        out.push(Entry { terms, cond: Expr::boolean(true), hint_false: false });
    }
    // log(Abs(t - a))
    if let Kind::Func(Func::Abs, inner) = args[0].kind() {
        if let Some((p, q, c)) = match_add_t(&inner[0], x) {
            let (p, q, a) = (pol(&p), pol(&q), pol(&neg(&c)));
            let t = mul2(&p, &pow(x, &q));
            let mut terms = constant_x(&crate::functions::log(&crate::functions::abs(&a)), x);
            terms.push((Expr::pi(), g(vec![int(1), int(1)], vec![h()], vec![int(1)], vec![int(0), h()], div(&t, &a))));
            out.push(Entry { terms, cond: Expr::boolean(true), hint_false: false });
        }
    }
}

/// `Heaviside(arg)` with the default value at zero.
fn heaviside_arg(e: &Expr) -> Option<&Expr> {
    match e.kind() {
        Kind::Func(Func::Heaviside, a) if a.len() == 1 || (a.len() == 2 && a[1] == Expr::rational(1, 2)) => Some(&a[0]),
        _ => None,
    }
}

fn log_heaviside_entries(f: &Expr, x: &Expr, out: &mut Vec<Entry>) {
    let fs = crate::simptools::mul_args(f);
    if fs.len() != 2 {
        return;
    }
    let (hv, lg) = if heaviside_arg(&fs[0]).is_some() { (&fs[0], &fs[1]) } else { (&fs[1], &fs[0]) };
    let (Some(harg), Some((p, q, n))) = (heaviside_arg(hv), log_power(lg, x)) else { return };
    let Some((hp, hq, c)) = match_add_t(harg, x) else { return };
    let (hp, hq, c) = (pol(&hp), pol(&hq), pol(&c));
    if hq != q {
        return;
    }
    let t = mul2(&p, &pow(x, &q));
    // log(t)**n*Heaviside(1 - t)
    if c.is_one() && hp == neg(&p) {
        out.push(Entry { terms: vec![make_log1(&t, n)], cond: Expr::boolean(true), hint_false: false });
    }
    // log(t)**n*Heaviside(t - 1)
    if c == Expr::minus_one() && hp == p {
        out.push(Entry { terms: vec![make_log2(&t, n)], cond: Expr::boolean(true), hint_false: false });
    }
}

/// The four Heaviside entries.
fn heaviside_entries(f: &Expr, x: &Expr, out: &mut Vec<Entry>) {
    let fs = crate::simptools::mul_args(f);
    let (harg, other) = match fs.as_slice() {
        [a] => (heaviside_arg(a), None),
        [a, b] if heaviside_arg(a).is_some() => (heaviside_arg(a), Some(b)),
        [a, b] => (heaviside_arg(b), Some(a)),
        _ => return,
    };
    let Some(harg) = harg else { return };
    let (pbase, pexp) = match other {
        Some(o) => {
            let (b, k) = base_exp(o);
            if has_x(&k, x) {
                return;
            }
            (Some(b), k)
        }
        None => (None, Expr::zero()),
    };
    let a_raw = add(&[pexp.clone(), Expr::one()]);
    let make = |p: Expr, q: Expr, b: Expr, upper: bool| -> Entry {
        let (p, q, b, a) = (pol(&p), pol(&q), pol(&b), pol(&a_raw));
        let t = mul2(&p, &pow(x, &q));
        let fac = mul2(&crate::functions::gamma(&a), &pow(&b, &sub(&a, &Expr::one())));
        let gg = if upper {
            g(vec![a.clone()], vec![], vec![], vec![int(0)], div(&t, &b))
        } else {
            g(vec![], vec![a.clone()], vec![int(0)], vec![], div(&t, &b))
        };
        Entry { terms: vec![(fac, gg)], cond: crate::relational::gt(&b, &Expr::zero()), hint_false: false }
    };
    let same_base = |b: &Expr| pbase.as_ref().is_none_or(|pb| pb == b);
    // Heaviside(t - b)*(t - b)**(a - 1)
    if let Some((p, q, c)) = match_add_t(harg, x) {
        if same_base(harg) {
            out.push(make(p.clone(), q.clone(), neg(&c), true));
            // Heaviside(b - t)*(b - t)**(a - 1)
            out.push(make(neg(&p), q, c, false));
        }
    }
    // Heaviside(x - (b/p)**(1/q))*(t - b)**(a - 1) and the mirrored form.
    let Some(pb) = pbase else { return };
    let Some((hp, hq, hc)) = match_add_t(harg, x) else { return };
    if hq != Expr::one() {
        return;
    }
    let Some((p, q, c)) = match_add_t(&pb, x) else { return };
    if hp.is_one() {
        let b = neg(&c);
        if pow(&div(&b, &p), &div(&Expr::one(), &q)) == neg(&hc) {
            out.push(make(p.clone(), q.clone(), b, true));
        }
    }
    if hp == Expr::minus_one() {
        let (p, b) = (neg(&p), c);
        if pow(&div(&b, &p), &div(&Expr::one(), &q)) == hc {
            out.push(make(p, q, b, false));
        }
    }
}

/// `_rewrite_single(f, x, recursive=False)`.
pub(crate) fn rewrite_single(f: &Expr, x: &Expr) -> Option<(Vec<(Expr, Expr, G)>, Expr)> {
    if let Some(gf) = G::from_expr(f) {
        let factored = crate::poly::api::factor(&gf.z).ok()?;
        let (coeff, m) = crate::decompose::as_coeff_mul(&factored, std::slice::from_ref(x), false);
        let [m] = m.as_slice() else { return None };
        match m.kind() {
            Kind::Pow(b, e) if b == x && e.as_number().is_some_and(|n| n.is_rational()) => {}
            _ if m == x => {}
            _ => return None,
        }
        let g = G::new(gf.an, gf.ao, gf.bm, gf.bo, mul2(&coeff, m))?;
        return Some((vec![(Expr::one(), Expr::zero(), g)], Expr::boolean(true)));
    }
    for entry in candidates(f, x) {
        if entry.hint_false {
            continue;
        }
        let cond = unpolarify(&entry.cond, false);
        if is_false(&super::conds::eval_cond(&cond)) {
            continue;
        }
        let mut res = Vec::new();
        for (fac, rg) in entry.terms {
            let (c, s) = get_coeff_exp(&unpolarify(&fac, true), x)?;
            let Some(g) = G::new(rg.an, rg.ao, rg.bm, rg.bo, rg.z) else { continue };
            if has_infinity(&c) || has_infinity(&s) || has_infinity(&g.z) {
                continue;
            }
            let z = unpolarify(&g.z, true);
            let Some(g) = G::new(g.an, g.ao, g.bm, g.bo, z) else { continue };
            res.push((c, s, g));
        }
        if !res.is_empty() {
            return Some((res, cond));
        }
    }
    None
}
