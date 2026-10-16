//! Printing support for `Order` terms and series (SymPy compatible).
//!
//! SymPy prints an `Add` containing an `Order` with its regular terms in
//! *ascending* lexicographic order (`x - x**3/6 + O(x**5)`), followed by the
//! order terms. The ordering is computed on `as_terms`-style monomials:
//! numeric factors are folded into a complex coefficient, other factors are
//! decomposed into `(generator, integer exponent)`.

use std::cmp::Ordering;

use crate::expr::{Expr, Func, Kind};
use crate::ops::{add_term, as_base_exp, as_coeff_mul, mul_factor, pow};
use crate::printer::to_str;

/// `str` of an `Order` node given its arguments.
pub fn order_str(args: &[Expr]) -> String {
    let mut s = String::from("O(");
    s.push_str(&to_str(&args[0]));
    let pairs: Vec<&[Expr]> = args[1..].chunks(2).collect();
    let all_zero = pairs.iter().all(|p| p.len() == 2 && p[1].is_zero());
    if all_zero {
        if pairs.len() > 1 {
            for p in &pairs {
                s.push_str(", ");
                s.push_str(&to_str(&p[0]));
            }
        }
    } else {
        for p in &pairs {
            s.push_str(", (");
            s.push_str(&to_str(&p[0]));
            if p.len() == 2 {
                s.push_str(", ");
                s.push_str(&to_str(&p[1]));
            }
            s.push(')');
        }
    }
    s.push(')');
    s
}

#[derive(Clone, Debug)]
struct Gen {
    class: (u8, u32),
    name: String,
    exp: f64,
    expr: Expr,
}

fn class_key(b: &Expr) -> (u8, u32) {
    match b.kind() {
        Kind::Num(_) | Kind::Const(_) => (1, 0),
        Kind::Sym(_) => (2, 0),
        Kind::Mul(..) => (3, 0),
        Kind::Add(..) => (3, 1),
        Kind::Pow(..) => (3, 2),
        Kind::Func(Func::Order, _) => (5, 0),
        Kind::Func(f, _) => (
            4,
            match f {
                Func::Exp => 10,
                Func::Log => 11,
                Func::Sin => 20,
                Func::Cos => 21,
                Func::Tan => 22,
                Func::Cot => 23,
                Func::Sinh => 30,
                Func::Cosh => 31,
                Func::Tanh => 32,
                Func::Coth => 33,
                Func::Conjugate => 40,
                Func::Re => 41,
                Func::Im => 42,
                Func::Arg => 43,
                Func::Undefined(_) => 0,
                _ => 10000,
            },
        ),
        _ => (5, 0),
    }
}

fn gen_of(g: &Expr) -> Gen {
    let (b, e) = as_base_exp(g);
    let exp = e.as_number().map(|n| n.to_f64());
    let (b, exp) = match exp {
        Some(v) => (b, v),
        None => (g.clone(), 1.0),
    };
    Gen { class: class_key(&b), name: to_str(&b), exp, expr: g.clone() }
}

fn cmp_gen(a: &Gen, b: &Gen) -> Ordering {
    a.class
        .cmp(&b.class)
        .then_with(|| a.name.cmp(&b.name))
        .then_with(|| a.exp.partial_cmp(&b.exp).unwrap_or(Ordering::Equal))
}

/// SymPy's `decompose_power`.
fn decompose_power(b: &Expr, e: &Expr) -> (Expr, i64) {
    if let Some(n) = e.as_number() {
        if n.is_rational() {
            if let (Some(p), Some(q)) = (n.numer().to_i64(), n.denom().to_i64()) {
                if q == 1 {
                    return (b.clone(), p);
                }
                return (pow(b, &Expr::rational(1, q)), p);
            }
        }
        return (mul_factor(b, e), 1);
    }
    let (k, tail) = as_coeff_mul(e);
    if k.is_minus_one() {
        return (pow(b, &tail), -1);
    }
    if k.is_rational() && !k.is_one() {
        if let (Some(p), Some(q)) = (k.numer().to_i64(), k.denom().to_i64()) {
            let t = crate::ops::mul2(&Expr::rational(1, q), &tail);
            return (pow(b, &t), p);
        }
    }
    (mul_factor(b, e), 1)
}

struct TermKey {
    monom: Vec<(Gen, i64)>,
    re: f64,
    im: f64,
}

fn term_key(t: &Expr) -> TermKey {
    let (c, rest) = as_coeff_mul(t);
    let mut re = c.to_f64();
    let mut im = 0.0;
    let mut monom: Vec<(Gen, i64)> = Vec::new();
    let facs: Vec<(Expr, Expr)> = match rest.kind() {
        Kind::Mul(_, fs) => fs.clone(),
        Kind::Num(_) => vec![],
        _ => vec![as_base_exp(&rest)],
    };
    for (b, e) in facs {
        let f = mul_factor(&b, &e);
        if crate::functions::is_numeric_value(&f) {
            if let Some(z) = crate::evalf::eval_c64(&f) {
                let (r2, i2) = (re * z.re - im * z.im, re * z.im + im * z.re);
                re = r2;
                im = i2;
                continue;
            }
        }
        let (g, k) = decompose_power(&b, &e);
        let gk = gen_of(&g);
        match monom.iter_mut().find(|(h, _)| cmp_gen(h, &gk) == Ordering::Equal && h.expr == gk.expr) {
            Some((_, d)) => *d += k,
            None => monom.push((gk, k)),
        }
    }
    TermKey { monom, re, im }
}

/// Ascending lex comparison of monomials (SymPy `reverse=True` order).
fn cmp_terms(a: &TermKey, b: &TermKey) -> Ordering {
    let mut gens: Vec<&Gen> = a.monom.iter().map(|x| &x.0).chain(b.monom.iter().map(|x| &x.0)).collect();
    gens.sort_by(|x, y| cmp_gen(x, y));
    gens.dedup_by(|x, y| cmp_gen(x, y) == Ordering::Equal && x.expr == y.expr);
    let deg = |m: &[(Gen, i64)], g: &Gen| -> i64 {
        m.iter()
            .filter(|(k, _)| cmp_gen(k, g) == Ordering::Equal && k.expr == g.expr)
            .map(|(_, d)| *d)
            .sum()
    };
    for g in gens {
        let (da, db) = (deg(&a.monom, g), deg(&b.monom, g));
        if da != db {
            return da.cmp(&db);
        }
    }
    // equal monomials: larger coefficient first
    let ka = (a.im != 0.0, a.im, a.re, a.im);
    let kb = (b.im != 0.0, b.im, b.re, b.im);
    kb.partial_cmp(&ka).unwrap_or(Ordering::Equal)
}

fn sort_series_terms(v: &mut Vec<Expr>) {
    let mut keyed: Vec<(TermKey, String, Expr)> = v.drain(..).map(|t| (term_key(&t), to_str(&t), t)).collect();
    keyed.sort_by(|a, b| cmp_terms(&a.0, &b.0).then_with(|| b.1.cmp(&a.1)));
    v.extend(keyed.into_iter().map(|k| k.2));
}

/// Terms of an `Add` containing order terms, in printing order.
pub fn series_ordered_terms(e: &Expr) -> Vec<Expr> {
    let Kind::Add(c, terms) = e.kind() else {
        return vec![e.clone()];
    };
    let mut normal: Vec<Expr> = Vec::with_capacity(terms.len() + 1);
    let mut orders: Vec<Expr> = Vec::new();
    for (t, k) in terms {
        if t.is_func(&Func::Order) {
            orders.push(t.clone());
        } else {
            normal.push(add_term(t, k));
        }
    }
    if !c.is_zero() {
        normal.push(Expr::num(c.clone()));
    }
    sort_series_terms(&mut normal);
    sort_series_terms(&mut orders);
    normal.extend(orders);
    normal
}
