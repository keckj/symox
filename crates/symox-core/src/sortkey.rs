//! SymPy's `default_sort_key` and the printing order of terms and factors
//! (`Add.as_ordered_terms`, `Mul.as_ordered_factors`).

use std::cmp::Ordering;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add_term, as_coeff_mul, mul, mul_factor, pow};

/// A sort key (mirrors the nested tuples SymPy builds).
#[derive(Clone, Debug)]
pub enum SortKey {
    Int(i64),
    Str(String),
    Num(Number),
    Tup(Vec<SortKey>),
}

impl SortKey {
    pub fn compare(&self, o: &SortKey) -> Ordering {
        match (self, o) {
            (SortKey::Int(a), SortKey::Int(b)) => a.cmp(b),
            (SortKey::Str(a), SortKey::Str(b)) => a.cmp(b),
            (SortKey::Num(a), SortKey::Num(b)) => a.cmp_value(b).unwrap_or(Ordering::Equal),
            (SortKey::Tup(a), SortKey::Tup(b)) => {
                for (x, y) in a.iter().zip(b) {
                    let c = x.compare(y);
                    if c != Ordering::Equal {
                        return c;
                    }
                }
                a.len().cmp(&b.len())
            }
            (SortKey::Int(a), SortKey::Num(b)) => Number::Int(*a).cmp_value(b).unwrap_or(Ordering::Equal),
            (SortKey::Num(a), SortKey::Int(b)) => a.cmp_value(&Number::Int(*b)).unwrap_or(Ordering::Equal),
            _ => rank(self).cmp(&rank(o)),
        }
    }
}

fn rank(k: &SortKey) -> u8 {
    match k {
        SortKey::Int(_) | SortKey::Num(_) => 0,
        SortKey::Str(_) => 1,
        SortKey::Tup(_) => 2,
    }
}

fn tup(v: Vec<SortKey>) -> SortKey {
    SortKey::Tup(v)
}

fn class_key3(a: i64, b: i64, name: &str) -> SortKey {
    tup(vec![SortKey::Int(a), SortKey::Int(b), SortKey::Str(name.to_string())])
}

fn func_class_index(f: &Func) -> i64 {
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
        Func::Undefined(_) | Func::Max | Func::Min => 0,
        _ => 10000,
    }
}

fn is_sympy_number(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => matches!(c, Constant::Infinity | Constant::NegInfinity | Constant::NaN),
        _ => false,
    }
}

/// SymPy's `class_key()`.
pub fn class_key(e: &Expr) -> SortKey {
    match e.kind() {
        _ if is_sympy_number(e) => class_key3(1, 0, "Number"),
        Kind::Sym(_) => class_key3(2, 0, "Symbol"),
        Kind::Const(_) | Kind::Bool(_) => class_key3(2, 0, e.class_name()),
        Kind::Mul(..) => class_key3(3, 0, "Mul"),
        Kind::Add(..) => class_key3(3, 1, "Add"),
        Kind::Pow(..) => class_key3(3, 2, "Pow"),
        Kind::Func(f, _) => class_key3(4, func_class_index(f), f.name()),
        _ => class_key3(5, 0, e.class_name()),
    }
}

fn number_key(e: &Expr) -> SortKey {
    let n = match e.kind() {
        Kind::Num(n) => SortKey::Num(n.clone()),
        Kind::Const(Constant::Infinity) => SortKey::Num(Number::float(f64::INFINITY)),
        Kind::Const(Constant::NegInfinity) => SortKey::Num(Number::float(f64::NEG_INFINITY)),
        _ => SortKey::Num(Number::float(f64::NAN)),
    };
    tup(vec![class_key(e), tup(vec![SortKey::Int(0), tup(vec![])]), tup(vec![]), n])
}

fn one_key() -> SortKey {
    number_key(&Expr::one())
}

/// SymPy's `default_sort_key(e)`.
pub fn sort_key(e: &Expr) -> SortKey {
    if is_sympy_number(e) {
        return number_key(e);
    }
    match e.kind() {
        Kind::Sym(s) => tup(vec![
            class_key(e),
            tup(vec![SortKey::Int(1), tup(vec![SortKey::Str(s.to_string())])]),
            one_key(),
            SortKey::Num(Number::ONE),
        ]),
        Kind::Const(_) | Kind::Bool(_) => tup(vec![
            class_key(e),
            tup(vec![SortKey::Int(1), tup(vec![SortKey::Str(crate::printer::to_str(e))])]),
            one_key(),
            SortKey::Num(Number::ONE),
        ]),
        Kind::Rel(..) | Kind::Piecewise(..) => {
            let args: Vec<SortKey> = e.args().iter().map(sort_key).collect();
            tup(vec![
                class_key(e),
                tup(vec![SortKey::Int(args.len() as i64), tup(args)]),
                one_key(),
                SortKey::Num(Number::ONE),
            ])
        }
        _ => {
            let (coeff, rest) = as_coeff_mul(e);
            let (expr, exp) = match rest.kind() {
                Kind::Pow(b, x) => {
                    if b.is_constant(Constant::E) {
                        (crate::functions::exp(x), Expr::one())
                    } else {
                        (b.clone(), x.clone())
                    }
                }
                _ => (rest.clone(), Expr::one()),
            };
            let args: Vec<SortKey> = if expr.is_atom() {
                vec![SortKey::Str(crate::printer::to_str(&expr))]
            } else {
                let parts = match expr.kind() {
                    Kind::Add(..) => ordered_terms(&expr),
                    Kind::Mul(..) => ordered_factors(&expr),
                    _ => expr.args(),
                };
                parts.iter().map(sort_key).collect()
            };
            tup(vec![
                class_key(&expr),
                tup(vec![SortKey::Int(args.len() as i64), tup(args)]),
                sort_key(&exp),
                SortKey::Num(coeff),
            ])
        }
    }
}

pub fn cmp_default(a: &Expr, b: &Expr) -> Ordering {
    sort_key(a).compare(&sort_key(b))
}

/// `Mul.as_ordered_factors()`: the factors (including the numeric
/// coefficient) sorted by `default_sort_key`.
pub fn ordered_factors(e: &Expr) -> Vec<Expr> {
    let mut v = e.args();
    if !e.is_mul() {
        return vec![e.clone()];
    }
    let mut keyed: Vec<(SortKey, Expr)> = v.drain(..).map(|f| (sort_key(&f), f)).collect();
    keyed.sort_by(|a, b| a.0.compare(&b.0));
    keyed.into_iter().map(|x| x.1).collect()
}

/// SymPy's `decompose_power`.
fn decompose_power(f: &Expr) -> (Expr, Expr) {
    let (b, x) = match f.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        Kind::Func(Func::Exp, a) => (Expr::e(), a[0].clone()),
        _ => return (f.clone(), Expr::one()),
    };
    let rebuild = |b: &Expr, x: &Expr| {
        if b.is_constant(Constant::E) {
            crate::functions::exp(x)
        } else {
            pow(b, x)
        }
    };
    if let Some(n) = x.as_number() {
        if n.is_rational() {
            let (p, q) = n.numer_denom().unwrap();
            let base = if q == num_bigint::BigInt::from(1) {
                b
            } else {
                rebuild(&b, &Expr::num(Number::rational(num_bigint::BigInt::from(1), q).unwrap()))
            };
            return (base, Expr::num(Number::from_bigint(p)));
        }
        return (f.clone(), Expr::one());
    }
    let (c, tail) = as_coeff_mul(&x);
    if !c.is_rational() || c.is_one() {
        return (f.clone(), Expr::one());
    }
    if c.is_minus_one() {
        return (rebuild(&b, &tail), Expr::minus_one());
    }
    let (p, q) = c.numer_denom().unwrap();
    let t = mul(&[Expr::num(Number::rational(num_bigint::BigInt::from(1), q).unwrap()), tail]);
    (rebuild(&b, &t), Expr::num(Number::from_bigint(p)))
}

fn is_number_like(e: &Expr) -> bool {
    crate::functions::is_numeric_value(e)
}

fn coeff_key(re: f64, im: f64) -> (bool, f64, f64, f64) {
    (im != 0.0, im, re, im)
}

fn cmp_f(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

/// `Add.as_ordered_terms()` (SymPy's printing order of the terms).
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    let Kind::Add(c, terms) = e.kind() else {
        return vec![e.clone()];
    };
    let mut args: Vec<Expr> = terms.iter().map(|(t, k)| add_term(t, k)).collect();
    if !c.is_zero() {
        args.push(Expr::num(c.clone()));
    }
    // special case Add(Number, Mul(-Number, expr))
    if args.len() == 2 {
        let is_num = |x: &Expr| is_sympy_number(x) || matches!(x.kind(), Kind::Const(Constant::Pi | Constant::E | Constant::EulerGamma));
        let (n, m) = if is_num(&args[0]) { (&args[0], &args[1]) } else { (&args[1], &args[0]) };
        if is_num(n) {
            if let Kind::Mul(k, facs) = m.kind() {
                let positive = match n.kind() {
                    Kind::Num(x) => x.is_positive(),
                    Kind::Const(Constant::Infinity) => true,
                    Kind::Const(Constant::Pi | Constant::E | Constant::EulerGamma) => true,
                    _ => false,
                };
                if facs.len() == 1 && !k.is_one() && k.is_negative() && positive {
                    return vec![n.clone(), m.clone()];
                }
            }
        }
    }
    struct T {
        expr: Expr,
        re: f64,
        im: f64,
        cpart: Vec<(Expr, Expr)>,
        s: String,
    }
    let mut gens: Vec<Expr> = Vec::new();
    let mut ts: Vec<T> = Vec::with_capacity(args.len());
    for t in &args {
        let (k, rest) = as_coeff_mul(t);
        let mut re = k.to_f64();
        let mut im = 0.0;
        let mut cpart: Vec<(Expr, Expr)> = Vec::new();
        if !rest.is_one() {
            let factors: Vec<Expr> = match rest.kind() {
                Kind::Mul(_, facs) => facs.iter().map(|(b, x)| mul_factor(b, x)).collect(),
                _ => vec![rest.clone()],
            };
            for f in factors {
                if is_number_like(&f) {
                    if let Some(z) = crate::evalf::eval_c64(&f) {
                        if z.re.is_finite() && z.im.is_finite() {
                            let (a, b) = (re * z.re - im * z.im, re * z.im + im * z.re);
                            re = a;
                            im = b;
                            continue;
                        }
                    }
                }
                let (base, x) = decompose_power(&f);
                if !gens.contains(&base) {
                    gens.push(base.clone());
                }
                cpart.push((base, x));
            }
        }
        ts.push(T { expr: t.clone(), re, im, cpart, s: crate::printer::to_str(t) });
    }
    let mut keyed: Vec<(SortKey, Expr)> = gens.iter().map(|g| (sort_key(g), g.clone())).collect();
    keyed.sort_by(|a, b| a.0.compare(&b.0));
    let gens: Vec<Expr> = keyed.into_iter().map(|x| x.1).collect();
    let monom = |t: &T| -> Vec<Number> {
        gens.iter()
            .map(|g| {
                t.cpart
                    .iter()
                    .find(|(b, _)| b == g)
                    .and_then(|(_, x)| x.as_number().cloned())
                    .unwrap_or(Number::ZERO)
            })
            .collect()
    };
    let mut items: Vec<(Vec<Number>, (bool, f64, f64, f64), &T)> =
        ts.iter().map(|t| (monom(t), coeff_key(t.re, t.im), t)).collect();
    items.sort_by(|a, b| {
        // negated monomial ascending == monomial descending (lex)
        for (x, y) in a.0.iter().zip(&b.0) {
            let c = y.cmp_value(x).unwrap_or(Ordering::Equal);
            if c != Ordering::Equal {
                return c;
            }
        }
        let (ka, kb) = (&a.1, &b.1);
        ka.0.cmp(&kb.0)
            .then_with(|| cmp_f(ka.1, kb.1))
            .then_with(|| cmp_f(ka.2, kb.2))
            .then_with(|| cmp_f(ka.3, kb.3))
            .then_with(|| a.2.s.cmp(&b.2.s))
    });
    items.into_iter().map(|x| x.2.expr.clone()).collect()
}
