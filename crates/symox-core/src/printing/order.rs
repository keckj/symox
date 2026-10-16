//! SymPy's printing order: `default_sort_key`, `Mul.as_ordered_factors`
//! and `Add.as_ordered_terms` (lexicographic order of the monomials over
//! the generators of the sum).
//!
//! Used by every printer (including `str`) so that terms and factors come
//! out in the same order as in SymPy.

use std::borrow::Cow;
use std::cmp::Ordering;

use rustc_hash::FxHashMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add_term, mul, mul_factor, pow};

// ----------------------------------------------------------------------
// default_sort_key
// ----------------------------------------------------------------------

#[derive(Clone, Debug)]
enum KArgs {
    /// `(0, ())` of numbers
    Empty,
    /// `(1, (str(atom),))`
    Atom(String),
    /// `(len(args), tuple(keys))`
    List(Vec<SortKey>),
}

/// Emulation of SymPy's `default_sort_key`.
#[derive(Clone, Debug)]
pub struct SortKey {
    class: (u8, u32, Cow<'static, str>),
    args: KArgs,
    /// `None` for numbers (their key has an empty exponent part).
    exp: Option<Box<SortKey>>,
    coeff: Number,
}

fn cmp_args(a: &KArgs, b: &KArgs) -> Ordering {
    let len = |k: &KArgs| match k {
        KArgs::Empty => 0,
        KArgs::Atom(_) => 1,
        KArgs::List(v) => v.len(),
    };
    len(a).cmp(&len(b)).then_with(|| match (a, b) {
        (KArgs::Empty, KArgs::Empty) => Ordering::Equal,
        (KArgs::Atom(x), KArgs::Atom(y)) => x.cmp(y),
        (KArgs::List(x), KArgs::List(y)) => {
            for (p, q) in x.iter().zip(y) {
                let r = p.cmp(q);
                if r != Ordering::Equal {
                    return r;
                }
            }
            Ordering::Equal
        }
        (KArgs::Atom(_), _) => Ordering::Less,
        (_, KArgs::Atom(_)) => Ordering::Greater,
        _ => Ordering::Equal,
    })
}

impl PartialEq for SortKey {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for SortKey {}
impl PartialOrd for SortKey {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for SortKey {
    fn cmp(&self, o: &Self) -> Ordering {
        self.class
            .cmp(&o.class)
            .then_with(|| cmp_args(&self.args, &o.args))
            .then_with(|| match (&self.exp, &o.exp) {
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
                (Some(a), Some(b)) => a.cmp(b),
            })
            .then_with(|| self.coeff.cmp_value(&o.coeff).unwrap_or(Ordering::Equal))
    }
}

/// `Function.class_key` ordinals of the known functions.
fn func_ordinal(f: &Func) -> u32 {
    match f {
        Func::Undefined(_) => 0,
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
        _ => 10000,
    }
}

fn class_key(e: &Expr) -> (u8, u32, Cow<'static, str>) {
    match e.kind() {
        Kind::Num(_) => (1, 0, Cow::Borrowed("Number")),
        Kind::Const(c) => match c {
            Constant::Infinity | Constant::NegInfinity | Constant::NaN => (1, 0, Cow::Borrowed("Number")),
            Constant::Pi => (2, 0, Cow::Borrowed("Pi")),
            Constant::E => (2, 0, Cow::Borrowed("Exp1")),
            Constant::I => (2, 0, Cow::Borrowed("ImaginaryUnit")),
            Constant::EulerGamma => (2, 0, Cow::Borrowed("EulerGamma")),
            Constant::ComplexInfinity => (2, 0, Cow::Borrowed("ComplexInfinity")),
        },
        Kind::Sym(_) => (2, 0, Cow::Borrowed("Symbol")),
        Kind::Mul(..) => (3, 0, Cow::Borrowed("Mul")),
        Kind::Add(..) => (3, 1, Cow::Borrowed("Add")),
        Kind::Pow(..) => (3, 2, Cow::Borrowed("Pow")),
        Kind::Func(Func::Max, _) => (5, 0, Cow::Borrowed("Max")),
        Kind::Func(Func::Min, _) => (5, 0, Cow::Borrowed("Min")),
        Kind::Func(f, _) => (4, func_ordinal(f), Cow::Owned(f.name().to_string())),
        Kind::Piecewise(_) => (4, 0, Cow::Borrowed("Piecewise")),
        Kind::Derivative(..) => (5, 0, Cow::Borrowed("Derivative")),
        Kind::Rel(..) => (5, 0, Cow::Owned(e.class_name().to_string())),
        Kind::Bool(_) => (5, 0, Cow::Owned(e.class_name().to_string())),
    }
}

fn one_key() -> SortKey {
    SortKey { class: (1, 0, Cow::Borrowed("Number")), args: KArgs::Empty, exp: None, coeff: Number::ONE }
}

/// SymPy's `default_sort_key(e)`.
pub fn sort_key(e: &Expr) -> SortKey {
    match e.kind() {
        Kind::Num(n) => SortKey { class: class_key(e), args: KArgs::Empty, exp: None, coeff: n.clone() },
        Kind::Const(c @ (Constant::Infinity | Constant::NegInfinity | Constant::NaN)) => {
            let v = match c {
                Constant::Infinity => f64::INFINITY,
                Constant::NegInfinity => f64::NEG_INFINITY,
                _ => f64::NAN,
            };
            SortKey { class: class_key(e), args: KArgs::Empty, exp: None, coeff: Number::Float(v) }
        }
        Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => SortKey {
            class: class_key(e),
            args: KArgs::Atom(crate::printer::to_str(e)),
            exp: Some(Box::new(one_key())),
            coeff: Number::ONE,
        },
        _ => {
            // coeff, expr = as_coeff_Mul()
            let (coeff, rest) = match e.kind() {
                Kind::Mul(c, facs) if !c.is_one() => {
                    let rest = if facs.len() == 1 {
                        mul_factor(&facs[0].0, &facs[0].1)
                    } else {
                        Expr::from_kind(Kind::Mul(Number::ONE, facs.clone()))
                    };
                    (c.clone(), rest)
                }
                _ => (Number::ONE, e.clone()),
            };
            let (expr, exp) = match rest.kind() {
                Kind::Pow(b, x) => (b.clone(), Some(x.clone())),
                _ => (rest, None),
            };
            let args = match expr.kind() {
                Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => {
                    KArgs::Atom(crate::printer::to_str(&expr))
                }
                Kind::Add(..) => KArgs::List(ordered_terms(&expr).iter().map(sort_key).collect()),
                Kind::Mul(..) => KArgs::List(ordered_factors(&expr).iter().map(sort_key).collect()),
                _ => KArgs::List(expr.args().iter().map(sort_key).collect()),
            };
            SortKey {
                class: class_key(&expr),
                args,
                exp: Some(Box::new(match exp {
                    Some(x) => sort_key(&x),
                    None => one_key(),
                })),
                coeff,
            }
        }
    }
}

// ----------------------------------------------------------------------
// factors
// ----------------------------------------------------------------------

/// Sort factors with `default_sort_key` (`Mul.as_ordered_factors`).
pub fn sort_factors(v: &mut [Expr]) {
    if v.len() < 2 {
        return;
    }
    v.sort_by_cached_key(sort_key);
}

/// `Mul.as_ordered_factors()`.
pub fn ordered_factors(e: &Expr) -> Vec<Expr> {
    match e.kind() {
        Kind::Mul(c, facs) => {
            let mut v: Vec<Expr> = Vec::with_capacity(facs.len() + 1);
            if !c.is_one() {
                v.push(Expr::num(c.clone()));
            }
            v.extend(facs.iter().map(|(b, x)| mul_factor(b, x)));
            sort_factors(&mut v);
            v
        }
        _ => vec![e.clone()],
    }
}

// ----------------------------------------------------------------------
// terms
// ----------------------------------------------------------------------

/// `decompose_power`: symbolic base and integer exponent.
fn decompose_power(f: &Expr) -> (Expr, i64) {
    let (b, x) = match f.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        Kind::Func(Func::Exp, a) => (Expr::e(), a[0].clone()),
        _ => return (f.clone(), 1),
    };
    let rebuild = |b: &Expr, t: &Expr| -> Expr {
        if b.is_constant(Constant::E) {
            crate::functions::exp(t)
        } else {
            pow(b, t)
        }
    };
    match x.kind() {
        Kind::Num(Number::Int(p)) => (b, *p),
        Kind::Num(Number::Rational(r)) => {
            use num_traits::ToPrimitive;
            match (r.numer().to_i64(), r.denom().to_i64()) {
                (Some(p), Some(q)) => (rebuild(&b, &Expr::rational(1, q)), p),
                _ => (f.clone(), 1),
            }
        }
        Kind::Num(_) => (f.clone(), 1),
        _ => {
            let (c, tail) = match x.kind() {
                Kind::Mul(c, facs) if c.is_rational() && !c.is_one() => {
                    let tail = if facs.len() == 1 {
                        mul_factor(&facs[0].0, &facs[0].1)
                    } else {
                        Expr::from_kind(Kind::Mul(Number::ONE, facs.clone()))
                    };
                    (c.clone(), tail)
                }
                _ => return (f.clone(), 1),
            };
            if c.is_minus_one() {
                return (rebuild(&b, &tail), -1);
            }
            let (p, q) = match c.numer_denom() {
                Some((p, q)) => {
                    use num_traits::ToPrimitive;
                    match (p.to_i64(), q.to_i64()) {
                        (Some(p), Some(q)) => (p, q),
                        _ => return (f.clone(), 1),
                    }
                }
                None => return (f.clone(), 1),
            };
            let t = if q == 1 { tail } else { mul(&[Expr::rational(1, q), tail]) };
            (rebuild(&b, &t), p)
        }
    }
}

fn is_number_like(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(_))
        || matches!(e.kind(), Kind::Const(c) if matches!(c, Constant::Pi | Constant::E | Constant::EulerGamma))
}

/// Lexicographic comparison of two sparse monomials (sorted by generator
/// index), as dense vectors.
fn cmp_lex(a: &[(usize, i64)], b: &[(usize, i64)]) -> Ordering {
    let (mut i, mut j) = (0, 0);
    loop {
        let ga = a.get(i).map(|x| x.0).unwrap_or(usize::MAX);
        let gb = b.get(j).map(|x| x.0).unwrap_or(usize::MAX);
        if ga == usize::MAX && gb == usize::MAX {
            return Ordering::Equal;
        }
        let (da, db) = if ga == gb {
            let r = (a[i].1, b[j].1);
            i += 1;
            j += 1;
            r
        } else if ga < gb {
            let r = (a[i].1, 0);
            i += 1;
            r
        } else {
            let r = (0, b[j].1);
            j += 1;
            r
        };
        if da != db {
            return da.cmp(&db);
        }
    }
}

fn cmp_f64(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

/// `Add.as_ordered_terms()`: the terms of a sum (with their coefficients,
/// the numeric term included) in SymPy's printing order.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    let Kind::Add(c, terms) = e.kind() else {
        return vec![e.clone()];
    };
    let mut all: Vec<Expr> = terms.iter().map(|(t, k)| add_term(t, k)).collect();
    if !c.is_zero() {
        all.push(Expr::num(c.clone()));
    }
    // special case: a positive number and a negative multiple of one
    // factor print number first ("1 - x", "pi - 2*x")
    if all.len() == 2 {
        let (n, m) = if is_number_like(&all[0]) { (0, 1) } else { (1, 0) };
        if is_number_like(&all[n]) && !is_number_like(&all[m]) {
            let positive = match all[n].kind() {
                Kind::Num(x) => x.is_positive(),
                _ => true,
            };
            if positive {
                if let Kind::Mul(k, facs) = all[m].kind() {
                    if facs.len() == 1 && k.is_negative() && !k.is_one() {
                        return vec![all[n].clone(), all[m].clone()];
                    }
                }
            }
        }
    }
    // decompose every term into coefficient and generator powers
    struct Term {
        expr: Expr,
        re: f64,
        im: f64,
        cpart: Vec<(Expr, i64)>,
    }
    let mut gens: FxHashMap<Expr, usize> = FxHashMap::default();
    let mut gen_list: Vec<Expr> = Vec::new();
    let mut decomposed: Vec<Term> = Vec::with_capacity(all.len());
    for t in all {
        let (coeff, factors): (Number, Vec<Expr>) = match t.kind() {
            Kind::Num(n) => (n.clone(), vec![]),
            Kind::Mul(c, facs) => (c.clone(), facs.iter().map(|(b, x)| mul_factor(b, x)).collect()),
            _ => (Number::ONE, vec![t.clone()]),
        };
        let mut z = num_complex::Complex64::new(coeff.to_f64(), 0.0);
        let mut cpart: Vec<(Expr, i64)> = Vec::new();
        for f in factors {
            if crate::functions::is_numeric_value(&f) {
                if let Some(v) = crate::evalf::eval_c64(&f) {
                    if v.re.is_finite() && v.im.is_finite() {
                        z *= v;
                        continue;
                    }
                }
            }
            let (base, d) = decompose_power(&f);
            if !gens.contains_key(&base) {
                gens.insert(base.clone(), gen_list.len());
                gen_list.push(base.clone());
            }
            // a generator may appear twice (x and sqrt(x) do not merge,
            // but exp(x) and exp(2*x) cannot both be factors): last wins
            cpart.retain(|(g, _)| *g != base);
            cpart.push((base, d));
        }
        decomposed.push(Term { expr: t, re: z.re, im: z.im, cpart });
    }
    // sort the generators with default_sort_key
    let mut order: Vec<usize> = (0..gen_list.len()).collect();
    if order.len() > 1 {
        let keys: Vec<SortKey> = gen_list.iter().map(sort_key).collect();
        order.sort_by(|a, b| keys[*a].cmp(&keys[*b]));
    }
    let mut rank = vec![0usize; gen_list.len()];
    for (r, g) in order.iter().enumerate() {
        rank[*g] = r;
    }
    let mut keyed: Vec<(Vec<(usize, i64)>, Term)> = decomposed
        .into_iter()
        .map(|t| {
            let mut m: Vec<(usize, i64)> = t.cpart.iter().map(|(g, d)| (rank[gens[g]], *d)).collect();
            m.sort_unstable();
            (m, t)
        })
        .collect();
    keyed.sort_by(|(ma, ta), (mb, tb)| {
        cmp_lex(mb, ma)
            .then_with(|| (ta.im != 0.0).cmp(&(tb.im != 0.0)))
            .then_with(|| cmp_f64(ta.im, tb.im))
            .then_with(|| cmp_f64(ta.re, tb.re))
            .then_with(|| sort_key(&ta.expr).cmp(&sort_key(&tb.expr)))
    });
    keyed.into_iter().map(|(_, t)| t.expr).collect()
}

#[cfg(test)]
mod tests {
    use crate::parser::parse_str;

    fn s(x: &str) -> String {
        parse_str(x).unwrap().to_string()
    }

    #[test]
    fn sympy_order() {
        assert_eq!(s("sqrt(x+1) + sin(x)"), "sqrt(x + 1) + sin(x)");
        assert_eq!(s("(x+1)*sin(x)"), "(x + 1)*sin(x)");
        assert_eq!(s("(x+y)*(x-y)*exp(x)"), "(x - y)*(x + y)*exp(x)");
        assert_eq!(s("1/(x+1) + 1/(x-1)"), "1/(x + 1) + 1/(x - 1)");
        assert_eq!(s("2**x + x**2"), "2**x + x**2");
        assert_eq!(s("x + 1 + 1/x"), "x + 1 + 1/x");
        assert_eq!(s("3 - 2*x"), "3 - 2*x");
        assert_eq!(s("sin(x) + cos(y)"), "sin(x) + cos(y)");
        assert_eq!(s("exp(x) + exp(-x)"), "exp(x) + exp(-x)");
        assert_eq!(s("sqrt(x) + x"), "sqrt(x) + x");
        assert_eq!(s("x**2 + 2*x*y + y**2"), "x**2 + 2*x*y + y**2");
    }
}
