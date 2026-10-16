//! SymPy's `default_sort_key`, `as_ordered_terms` and `as_ordered_factors`.
//!
//! These orderings only matter for printing: internally sums and products
//! are kept sorted by hash, which is much cheaper. Reproducing SymPy's
//! orderings exactly is what makes the printed output identical.

use std::cmp::Ordering;

use num_complex::Complex64;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::is_numeric_value;
use crate::number::Number;
use crate::ops::{as_coeff_mul, mul_factor, pow};

/// A Python-tuple like sort key.
#[derive(Clone, Debug)]
pub enum Key {
    Int(i64),
    Str(String),
    Num(Number),
    Tup(Vec<Key>),
}

impl Key {
    fn rank(&self) -> u8 {
        match self {
            Key::Int(_) => 0,
            Key::Num(_) => 0,
            Key::Str(_) => 1,
            Key::Tup(_) => 2,
        }
    }
}

fn cmp_num(a: &Number, b: &Number) -> Ordering {
    a.cmp_value(b).unwrap_or_else(|| a.cmp_struct(b))
}

impl Ord for Key {
    fn cmp(&self, o: &Key) -> Ordering {
        match (self, o) {
            (Key::Int(a), Key::Int(b)) => a.cmp(b),
            (Key::Num(a), Key::Num(b)) => cmp_num(a, b),
            (Key::Int(a), Key::Num(b)) => cmp_num(&Number::Int(*a), b),
            (Key::Num(a), Key::Int(b)) => cmp_num(a, &Number::Int(*b)),
            (Key::Str(a), Key::Str(b)) => a.cmp(b),
            (Key::Tup(a), Key::Tup(b)) => {
                for (x, y) in a.iter().zip(b) {
                    let r = x.cmp(y);
                    if r != Ordering::Equal {
                        return r;
                    }
                }
                a.len().cmp(&b.len())
            }
            _ => self.rank().cmp(&o.rank()),
        }
    }
}
impl PartialOrd for Key {
    fn partial_cmp(&self, o: &Key) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl PartialEq for Key {
    fn eq(&self, o: &Key) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Key {}

fn tup(v: Vec<Key>) -> Key {
    Key::Tup(v)
}

fn class_key(e: &Expr) -> Key {
    let (a, b, n): (i64, i64, String) = match e.kind() {
        Kind::Num(_) => (1, 0, "Number".into()),
        Kind::Const(c) => match c {
            Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN => {
                (1, 0, "Number".into())
            }
            _ => (2, 0, e.class_name().into()),
        },
        Kind::Sym(_) | Kind::Bool(_) => (2, 0, e.class_name().into()),
        Kind::Pow(..) => (3, 2, "Pow".into()),
        Kind::Mul(..) => (3, 0, "Mul".into()),
        Kind::Add(..) => (3, 1, "Add".into()),
        Kind::Func(f @ (Func::Logic(_) | Func::Set(_)), _) => (5, 0, f.name().into()),
        Kind::Func(f, _) => {
            let i = match f {
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
            };
            (4, i, f.name().into())
        }
        Kind::Piecewise(_) => (4, 0, "Piecewise".into()),
        _ => (5, 0, e.class_name().into()),
    };
    tup(vec![Key::Int(a), Key::Int(b), Key::Str(n)])
}

fn number_key(n: Number) -> Key {
    tup(vec![
        tup(vec![Key::Int(1), Key::Int(0), Key::Str("Number".into())]),
        tup(vec![Key::Int(0), tup(vec![])]),
        tup(vec![]),
        Key::Num(n),
    ])
}

/// SymPy's `default_sort_key`.
pub fn sort_key(e: &Expr) -> Key {
    match e.kind() {
        Kind::Num(n) => return number_key(n.clone()),
        Kind::Const(Constant::Infinity) => return number_key(Number::Float(f64::INFINITY)),
        Kind::Const(Constant::NegInfinity) => return number_key(Number::Float(f64::NEG_INFINITY)),
        Kind::Const(Constant::ComplexInfinity | Constant::NaN) => {
            return number_key(Number::Float(f64::MAX))
        }
        _ => {}
    }
    let (coeff, rest) = as_coeff_mul(e);
    let (base, exp) = match rest.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        _ => (rest.clone(), Expr::one()),
    };
    let args: Vec<Key> = if base.is_atom() {
        vec![Key::Str(crate::printer::to_str(&base))]
    } else {
        match base.kind() {
            Kind::Add(..) => ordered_terms(&base).iter().map(sort_key).collect(),
            Kind::Mul(..) => ordered_factors(&base).iter().map(sort_key).collect(),
            _ => base.args().iter().map(sort_key).collect(),
        }
    };
    tup(vec![
        class_key(&base),
        tup(vec![Key::Int(args.len() as i64), tup(args)]),
        sort_key(&exp),
        Key::Num(coeff),
    ])
}

/// SymPy's `as_ordered_factors` for a product (numbers first).
pub fn ordered_factors(e: &Expr) -> Vec<Expr> {
    let mut v = match e.kind() {
        Kind::Mul(..) => e.args(),
        _ => return vec![e.clone()],
    };
    let mut keyed: Vec<(Key, Expr)> = v.drain(..).map(|f| (sort_key(&f), f)).collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0));
    keyed.into_iter().map(|x| x.1).collect()
}

/// SymPy's `decompose_power`: `(base, integer exponent)`.
fn decompose_power(f: &Expr) -> (Expr, i64) {
    let (b, x) = match f.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        Kind::Func(Func::Exp, a) => (Expr::e(), a[0].clone()),
        _ => return (f.clone(), 1),
    };
    if let Some(n) = x.as_number() {
        return match n {
            Number::Int(k) => (b, *k),
            Number::Rational(r) => {
                let q = num_traits::ToPrimitive::to_i64(r.denom()).unwrap_or(1);
                let p = num_traits::ToPrimitive::to_i64(r.numer()).unwrap_or(1);
                (pow(&b, &Expr::rational(1, q)), p)
            }
            _ => (f.clone(), 1),
        };
    }
    let (c, tail) = as_coeff_mul(&x);
    if !c.is_rational() || c.is_one() {
        return (f.clone(), 1);
    }
    if c.is_minus_one() {
        return (pow(&b, &tail), -1);
    }
    let p = c.numer().to_i64().unwrap_or(1);
    let q = c.denom().to_i64().unwrap_or(1);
    let t = crate::ops::mul2(&Expr::rational(1, q), &tail);
    (pow(&b, &t), p)
}

/// SymPy's `as_ordered_terms` (default `lex` order) for a sum.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    let args = match e.kind() {
        Kind::Add(..) => e.args(),
        _ => return vec![e.clone()],
    };
    // special case: Number + Mul(negative Number, expr) prints number first
    if args.len() == 2 {
        let (num, other): (Vec<&Expr>, Vec<&Expr>) =
            args.iter().partition(|a| a.is_number() || is_number_symbol(a));
        if num.len() == 1 && other.len() == 1 {
            if let Kind::Mul(c, facs) = other[0].kind() {
                let positive = match num[0].kind() {
                    Kind::Num(n) => n.is_positive(),
                    _ => true,
                };
                if facs.len() == 1 && c.is_negative() && positive {
                    return vec![num[0].clone(), other[0].clone()];
                }
            }
        }
    }
    struct T {
        term: Expr,
        coeff: Complex64,
        cpart: Vec<(Expr, i64)>,
    }
    let mut terms: Vec<T> = Vec::with_capacity(args.len());
    let mut gens: Vec<Expr> = Vec::new();
    for t in args {
        let (c, rest) = as_coeff_mul(&t);
        let mut coeff = Complex64::new(c.to_f64(), 0.0);
        let mut cpart: Vec<(Expr, i64)> = Vec::new();
        if !rest.is_one() {
            let factors = match rest.kind() {
                Kind::Mul(_, fs) => fs.iter().map(|(b, x)| mul_factor(b, x)).collect(),
                _ => vec![rest.clone()],
            };
            for f in factors {
                if is_numeric_value(&f) {
                    if let Some(z) = crate::evalf::eval_c64(&f) {
                        coeff *= z;
                        continue;
                    }
                }
                let (b, k) = decompose_power(&f);
                if !gens.contains(&b) {
                    gens.push(b.clone());
                }
                cpart.push((b, k));
            }
        }
        terms.push(T { term: t, coeff, cpart });
    }
    let mut gkeys: Vec<(Key, Expr)> = gens.into_iter().map(|g| (sort_key(&g), g)).collect();
    gkeys.sort_by(|a, b| a.0.cmp(&b.0));
    let gens: Vec<Expr> = gkeys.into_iter().map(|x| x.1).collect();
    let mut keyed: Vec<(Vec<i64>, (bool, f64, f64), Expr)> = terms
        .into_iter()
        .map(|t| {
            let mut monom = vec![0i64; gens.len()];
            for (b, k) in &t.cpart {
                let i = gens.iter().position(|g| g == b).unwrap();
                monom[i] = *k;
            }
            (monom, (t.coeff.im != 0.0, t.coeff.im, t.coeff.re), t.term)
        })
        .collect();
    keyed.sort_by(|a, b| {
        // descending lex on monomials, then ascending coefficient
        b.0.cmp(&a.0).then_with(|| {
            let (ai, aim, are) = a.1;
            let (bi, bim, bre) = b.1;
            ai.cmp(&bi)
                .then(aim.partial_cmp(&bim).unwrap_or(Ordering::Equal))
                .then(are.partial_cmp(&bre).unwrap_or(Ordering::Equal))
                .then(aim.partial_cmp(&bim).unwrap_or(Ordering::Equal))
        })
    });
    keyed.into_iter().map(|x| x.2).collect()
}

fn is_number_symbol(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(Constant::Pi | Constant::E | Constant::EulerGamma))
}
