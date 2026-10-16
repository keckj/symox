//! SymPy's canonical printing order: `sort_key` / `default_sort_key` and
//! `Add.as_ordered_terms`, used by the printer to order the factors of a
//! product and the terms of a sum exactly like SymPy.

use std::cmp::Ordering;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add_term, mul_factor, pow};

/// A sort key value, compared like the corresponding Python objects.
#[derive(Clone, Debug)]
pub enum SK {
    Int(i64),
    Str(String),
    Num(Number),
    /// numerical key for infinities / nan (`None` for nan)
    Real(f64),
    Tup(Vec<SK>),
}

impl SK {
    fn rank(&self) -> u8 {
        match self {
            SK::Int(_) | SK::Num(_) | SK::Real(_) => 0,
            SK::Str(_) => 1,
            SK::Tup(_) => 2,
        }
    }
}

fn num_f64(n: &Number) -> f64 {
    n.to_f64()
}

impl PartialEq for SK {
    fn eq(&self, o: &SK) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for SK {}
impl PartialOrd for SK {
    fn partial_cmp(&self, o: &SK) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

impl Ord for SK {
    fn cmp(&self, o: &SK) -> Ordering {
        match (self, o) {
            (SK::Int(a), SK::Int(b)) => a.cmp(b),
            (SK::Str(a), SK::Str(b)) => a.cmp(b),
            (SK::Num(a), SK::Num(b)) => a.cmp_value(b).unwrap_or(Ordering::Equal),
            (SK::Int(a), SK::Num(b)) => Number::Int(*a).cmp_value(b).unwrap_or(Ordering::Equal),
            (SK::Num(a), SK::Int(b)) => a.cmp_value(&Number::Int(*b)).unwrap_or(Ordering::Equal),
            (SK::Real(a), SK::Real(b)) => a.partial_cmp(b).unwrap_or(Ordering::Equal),
            (SK::Real(a), SK::Num(b)) => a.partial_cmp(&num_f64(b)).unwrap_or(Ordering::Equal),
            (SK::Num(a), SK::Real(b)) => num_f64(a).partial_cmp(b).unwrap_or(Ordering::Equal),
            (SK::Real(a), SK::Int(b)) => a.partial_cmp(&(*b as f64)).unwrap_or(Ordering::Equal),
            (SK::Int(a), SK::Real(b)) => (*a as f64).partial_cmp(b).unwrap_or(Ordering::Equal),
            (SK::Tup(a), SK::Tup(b)) => {
                for (x, y) in a.iter().zip(b) {
                    let c = x.cmp(y);
                    if c != Ordering::Equal {
                        return c;
                    }
                }
                a.len().cmp(&b.len())
            }
            _ => self.rank().cmp(&o.rank()),
        }
    }
}

fn tup(v: Vec<SK>) -> SK {
    SK::Tup(v)
}

fn class_key3(a: i64, b: i64, name: &str) -> SK {
    tup(vec![SK::Int(a), SK::Int(b), SK::Str(name.to_string())])
}

fn func_index(f: &Func) -> i64 {
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
    }
}

fn is_number_class(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => matches!(c, Constant::Infinity | Constant::NegInfinity | Constant::NaN),
        _ => false,
    }
}

fn class_key(e: &Expr) -> SK {
    match e.kind() {
        Kind::Num(_) => class_key3(1, 0, "Number"),
        Kind::Const(c) => match c {
            Constant::Infinity | Constant::NegInfinity | Constant::NaN => class_key3(1, 0, "Number"),
            Constant::Pi => class_key3(2, 0, "Pi"),
            Constant::E => class_key3(2, 0, "Exp1"),
            Constant::I => class_key3(2, 0, "ImaginaryUnit"),
            Constant::EulerGamma => class_key3(2, 0, "EulerGamma"),
            Constant::ComplexInfinity => class_key3(2, 0, "ComplexInfinity"),
        },
        Kind::Sym(_) => class_key3(2, 0, "Symbol"),
        Kind::Bool(b) => class_key3(2, 0, if *b { "BooleanTrue" } else { "BooleanFalse" }),
        Kind::Mul(..) => class_key3(3, 0, "Mul"),
        Kind::Add(..) => class_key3(3, 1, "Add"),
        Kind::Pow(..) => class_key3(3, 2, "Pow"),
        Kind::Func(Func::Max, _) | Kind::Func(Func::Min, _) => class_key3(5, 0, e.class_name()),
        Kind::Func(f, _) => class_key3(4, func_index(f), f.name()),
        Kind::Piecewise(_) => class_key3(4, 0, "Piecewise"),
        _ => class_key3(5, 0, e.class_name()),
    }
}

fn number_value_key(e: &Expr) -> SK {
    match e.kind() {
        Kind::Num(n) => SK::Num(n.clone()),
        Kind::Const(Constant::Infinity) => SK::Real(f64::INFINITY),
        Kind::Const(Constant::NegInfinity) => SK::Real(f64::NEG_INFINITY),
        _ => SK::Real(f64::NAN),
    }
}

fn one_key() -> SK {
    tup(vec![class_key3(1, 0, "Number"), tup(vec![SK::Int(0), tup(vec![])]), tup(vec![]), SK::Num(Number::ONE)])
}

fn atom_key(e: &Expr) -> SK {
    tup(vec![
        class_key(e),
        tup(vec![SK::Int(1), tup(vec![SK::Str(crate::printer::to_str(e))])]),
        one_key(),
        SK::Num(Number::ONE),
    ])
}

/// `(coeff, rest)` like `as_coeff_Mul()`.
fn as_coeff_mul(e: &Expr) -> (Number, Expr) {
    match e.kind() {
        Kind::Mul(c, facs) if !c.is_one() => {
            let rest = if facs.len() == 1 {
                mul_factor(&facs[0].0, &facs[0].1)
            } else {
                Expr::from_kind(Kind::Mul(Number::ONE, facs.clone()))
            };
            (c.clone(), rest)
        }
        _ => (Number::ONE, e.clone()),
    }
}

/// SymPy's `default_sort_key`.
pub fn sort_key(e: &Expr) -> SK {
    if is_number_class(e) {
        return tup(vec![class_key(e), tup(vec![SK::Int(0), tup(vec![])]), tup(vec![]), number_value_key(e)]);
    }
    if e.is_atom() {
        return atom_key(e);
    }
    let (coeff, mut ex) = as_coeff_mul(e);
    let mut exp = Expr::one();
    if let Kind::Pow(b, x) = ex.kind() {
        // as_base_exp: (1/q)**e -> q**(-e)
        let (b, x) = match b.as_number() {
            Some(Number::Rational(r)) if r.numer() == &num_bigint::BigInt::from(1) => {
                (Expr::num(Number::from_bigint(r.denom().clone())), crate::ops::neg(x))
            }
            _ => (b.clone(), x.clone()),
        };
        ex = b;
        exp = x;
    }
    let args = if ex.is_atom() {
        tup(vec![SK::Int(1), tup(vec![SK::Str(crate::printer::to_str(&ex))])])
    } else {
        let keys: Vec<SK> = match ex.kind() {
            Kind::Add(..) => ordered_terms(&ex).iter().map(sort_key).collect(),
            Kind::Mul(..) => {
                let mut ks: Vec<SK> = ex.args().iter().map(sort_key).collect();
                ks.sort();
                ks
            }
            _ => ex.args().iter().map(sort_key).collect(),
        };
        tup(vec![SK::Int(keys.len() as i64), tup(keys)])
    };
    tup(vec![class_key(&ex), args, sort_key(&exp), SK::Num(coeff)])
}

// ----------------------------------------------------------------------
// Add.as_ordered_terms
// ----------------------------------------------------------------------

/// `exprtools.decompose_power`.
fn decompose_power(f: &Expr) -> (Expr, i64) {
    let is_exp = matches!(f.kind(), Kind::Func(Func::Exp, _));
    let (b, x) = match f.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        Kind::Func(Func::Exp, a) => (Expr::e(), a[0].clone()),
        _ => return (f.clone(), 1),
    };
    let mk = |t: &Expr| -> Expr {
        if is_exp {
            crate::functions::exp(t)
        } else {
            pow(&b, t)
        }
    };
    if let Some(n) = x.as_number() {
        return match n {
            Number::Int(k) => (b.clone(), *k),
            Number::Rational(r) => {
                let p = num_traits::ToPrimitive::to_i64(r.numer());
                let q = num_traits::ToPrimitive::to_i64(r.denom());
                match (p, q) {
                    (Some(p), Some(q)) => (mk(&Expr::rational(1, q)), p),
                    _ => (f.clone(), 1),
                }
            }
            _ => (f.clone(), 1),
        };
    }
    let (c, tail) = crate::ops::as_coeff_mul(&x);
    if !c.is_rational() || c.is_one() {
        return (f.clone(), 1);
    }
    if c.is_minus_one() {
        return (mk(&tail), -1);
    }
    let (p, q) = c.numer_denom().unwrap();
    match (num_traits::ToPrimitive::to_i64(&p), num_traits::ToPrimitive::to_i64(&q)) {
        (Some(p), Some(q)) => {
            let t = crate::ops::mul(&[Expr::rational(1, q), tail]);
            (mk(&t), p)
        }
        _ => (f.clone(), 1),
    }
}

fn is_number_expr(e: &Expr) -> bool {
    if e.mask() != 0 {
        return false;
    }
    let mut ok = true;
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Func(Func::Undefined(_), _) | Kind::Sym(_) | Kind::Derivative(..)) {
            ok = false;
        }
    });
    ok
}

fn is_number_or_numbersymbol(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => !matches!(c, Constant::I | Constant::ComplexInfinity),
        _ => false,
    }
}

fn is_positive_atom(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_positive(),
        Kind::Const(c) => matches!(c, Constant::Pi | Constant::E | Constant::EulerGamma | Constant::Infinity),
        _ => false,
    }
}

/// Terms of a sum in SymPy's printing order.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    let Kind::Add(c, terms) = e.kind() else {
        return vec![e.clone()];
    };
    let mut all: Vec<Expr> = terms.iter().map(|(t, k)| add_term(t, k)).collect();
    if !c.is_zero() {
        all.push(Expr::num(c.clone()));
    }
    // special case: Add(positive number, Mul(negative number, expr))
    if all.len() == 2 {
        let (num, other): (Vec<&Expr>, Vec<&Expr>) = all.iter().partition(|t| is_number_or_numbersymbol(t));
        if num.len() == 1 && other.len() == 1 {
            if let Kind::Mul(mc, facs) = other[0].kind() {
                let nargs = facs.len() + usize::from(!mc.is_one());
                if nargs == 2 && !mc.is_one() && mc.is_negative() && is_positive_atom(num[0]) {
                    return vec![num[0].clone(), other[0].clone()];
                }
            }
        }
    }
    // as_terms
    struct T {
        term: Expr,
        re: f64,
        im: f64,
        cpart: Vec<(Expr, i64)>,
    }
    let mut ts: Vec<T> = Vec::with_capacity(all.len());
    let mut gens: Vec<Expr> = Vec::new();
    for t in all {
        let (k, rest) = as_coeff_mul(&t);
        let mut z = num_complex::Complex64::new(k.to_f64(), 0.0);
        let mut cpart: Vec<(Expr, i64)> = Vec::new();
        let factors: Vec<Expr> = match rest.kind() {
            Kind::Num(n) if n.is_one() => vec![],
            Kind::Num(n) => {
                z *= n.to_f64();
                vec![]
            }
            Kind::Mul(_, facs) => facs.iter().map(|(b, x)| mul_factor(b, x)).collect(),
            _ => vec![rest.clone()],
        };
        for f in factors {
            if is_number_class(&f) {
                if let SK::Real(v) = number_value_key(&f) {
                    z = num_complex::Complex64::new(v * z.re.signum(), 0.0);
                    continue;
                }
            }
            if is_number_expr(&f) {
                if let Some(v) = crate::evalf::eval_c64(&f) {
                    if v.re.is_finite() && v.im.is_finite() {
                        z *= v;
                        continue;
                    }
                }
            }
            let (b, kk) = decompose_power(&f);
            if !gens.contains(&b) {
                gens.push(b.clone());
            }
            cpart.push((b, kk));
        }
        ts.push(T { term: t, re: z.re, im: z.im, cpart });
    }
    let gkeys: Vec<SK> = gens.iter().map(sort_key).collect();
    let mut order: Vec<usize> = (0..gens.len()).collect();
    order.sort_by(|&a, &b| gkeys[a].cmp(&gkeys[b]));
    let mut pos = vec![0usize; gens.len()];
    for (i, &g) in order.iter().enumerate() {
        pos[g] = i;
    }
    let k = gens.len();
    let mut keyed: Vec<(Vec<i64>, (bool, f64), (f64, f64), Expr)> = ts
        .into_iter()
        .map(|t| {
            let mut m = vec![0i64; k];
            for (b, e) in &t.cpart {
                let i = gens.iter().position(|g| g == b).unwrap();
                m[pos[i]] = -e;
            }
            (m, (t.im != 0.0, t.im), (t.re, t.im), t.term)
        })
        .collect();
    keyed.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1 .0.cmp(&b.1 .0))
            .then_with(|| a.1 .1.partial_cmp(&b.1 .1).unwrap_or(Ordering::Equal))
            .then_with(|| a.2 .0.partial_cmp(&b.2 .0).unwrap_or(Ordering::Equal))
            .then_with(|| a.2 .1.partial_cmp(&b.2 .1).unwrap_or(Ordering::Equal))
    });
    keyed.into_iter().map(|x| x.3).collect()
}

/// Sort factors of a product like `Mul.as_ordered_factors`.
pub fn sort_factors(v: &mut [Expr]) {
    v.sort_by_cached_key(sort_key);
}
