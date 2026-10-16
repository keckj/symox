//! SymPy's `default_sort_key`, `as_ordered_terms` and `as_ordered_factors`.
//!
//! These orderings only matter for printing: internally sums and products
//! are kept sorted by hash, which is much cheaper. Reproducing SymPy's
//! orderings exactly is what makes the printed output identical.

use std::cmp::Ordering;
use std::sync::Arc;
use rustc_hash::FxHashMap;

use num_complex::Complex64;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::is_numeric_value;
use crate::number::Number;
use crate::ops::{as_coeff_mul, mul_factor, pow};
use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand};

/// Presentation arguments deliberately differ from scalar `Expr::args()`:
/// matrix parents live in typed payloads, not in commutative scalar children.
pub(crate) enum PresentationArg {
    Scalar(Expr),
    Matrix(MatrixExpr),
    Text(String),
    Node(&'static str, Vec<PresentationArg>),
}

pub(crate) struct MatrixPresentation {
    pub class: (u8, u32, &'static str),
    pub args: Vec<PresentationArg>,
}

pub(crate) fn matrix_presentation(matrix: &MatrixExpr) -> MatrixPresentation {
    use PresentationArg::{Matrix as M, Scalar as S};
    let mut class = (5, 0, matrix.kind_name());
    let args = match matrix.kind() {
        MatrixKind::Symbol(name) => vec![PresentationArg::Text(name.clone()), S(matrix.rows().clone()), S(matrix.cols().clone())],
        MatrixKind::Explicit(m) => {
            class.2 = "ImmutableDenseMatrix";
            vec![S(matrix.rows().clone()), S(matrix.cols().clone()), PresentationArg::Node("Tuple", m.data().iter().cloned().map(S).collect())]
        }
        MatrixKind::Block(grid) => vec![PresentationArg::Node("ImmutableDenseMatrix", vec![
            S(Expr::int(grid.shape().0 as i64)),
            S(Expr::int(grid.shape().1 as i64)),
            PresentationArg::Node("Tuple", grid.data().iter().cloned().map(M).collect()),
        ])],
        MatrixKind::Identity => vec![S(matrix.rows().clone())],
        MatrixKind::Zero | MatrixKind::One => vec![S(matrix.rows().clone()), S(matrix.cols().clone())],
        MatrixKind::HadamardProduct(factors) | MatrixKind::Kronecker(factors) | MatrixKind::BlockDiag(factors) => factors.iter().cloned().map(M).collect(),
        MatrixKind::Slice(parent, rows, cols) => vec![M(parent.clone()),
            PresentationArg::Node("Tuple", rows.iter().cloned().map(S).collect()),
            PresentationArg::Node("Tuple", cols.iter().cloned().map(S).collect())],
        MatrixKind::HadamardPower(base, exponent) => [base, exponent].into_iter().map(|operand| match operand {
            MatrixOperand::Scalar(e) => S(e.clone()), MatrixOperand::Matrix(m) => M(m.clone()),
        }).collect(),
        MatrixKind::Add(terms) => {
            class = (3, 1, "MatAdd");
            terms.iter().cloned().map(M).collect()
        }
        MatrixKind::Mul(c, factors) => {
            class = (3, 0, "MatMul");
            let mut args = Vec::new();
            if !c.is_one() { args.push(S(c.clone())); }
            args.extend(factors.iter().cloned().map(M)); args
        }
        MatrixKind::Pow(base, exponent) => vec![M(base.clone()), S(exponent.clone())],
        MatrixKind::Inverse(base) => vec![M(base.clone()), S(Expr::minus_one())],
        MatrixKind::Transpose(base) | MatrixKind::Adjoint(base) => vec![M(base.clone())],
        MatrixKind::Conjugate(base) => {
            class.2 = "Adjoint";
            vec![PresentationArg::Node("Transpose", vec![M(base.clone())])]
        }
    };
    MatrixPresentation { class, args }
}

/// Matrix sums order scalar monomials first, then ordered noncommutative
/// factors. A scalar coefficient must not make `2*B` sort before `A`.
pub(crate) fn matrix_ordered_terms<K: Ord>(terms: &[MatrixExpr], mut matrix_key: impl FnMut(&MatrixExpr) -> K) -> Vec<MatrixExpr> {
    let mut gens = Vec::new();
    let mut entries = Vec::new();
    for term in terms {
        let (scalar, factors) = match term.kind() {
            MatrixKind::Mul(c, factors) => (c.clone(), factors.clone()),
            _ => (Expr::one(), vec![term.clone()]),
        };
        let numeric = if is_numeric_value(&scalar) { crate::evalf::eval_c64(&scalar) } else { None };
        let coefficient = numeric.unwrap_or(Complex64::new(1.0, 0.0));
        let mut monomial = Vec::new();
        // MatMul stores its entire scalar coefficient as one argument. For
        // example x*y remains one generator here, matching SymPy's as_terms.
        if numeric.is_none() && !scalar.is_one() {
            let (base, exponent) = decompose_power(&scalar);
            if !gens.contains(&base) { gens.push(base.clone()); }
            monomial.push((base, exponent));
        }
        entries.push((term.clone(), monomial, factors.iter().map(&mut matrix_key).collect::<Vec<_>>(), coefficient));
    }
    gens.sort_by_cached_key(sort_key);
    let mut keyed: Vec<_> = entries.into_iter().map(|(term, monomial, ncpart, coeff)| {
        let mut powers = vec![0i64; gens.len()];
        for (base, exponent) in monomial { powers[gens.iter().position(|g| g == &base).unwrap()] = exponent; }
        (powers, ncpart, coeff, term)
    }).collect();
    keyed.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1))
        .then_with(|| (a.2.im != 0.0).cmp(&(b.2.im != 0.0)))
        .then_with(|| a.2.im.total_cmp(&b.2.im)).then_with(|| a.2.re.total_cmp(&b.2.re)));
    keyed.into_iter().map(|(_, _, _, term)| term).collect()
}

fn composite_key(class: (u8, u32, &str), args: Vec<Key>) -> Key {
    tup(vec![
        tup(vec![Key::Int(class.0.into()), Key::Int(class.1.into()), Key::Str(class.2.into())]),
        tup(vec![Key::Int(args.len() as i64), tup(args)]),
        number_key(Number::ONE), Key::Num(Number::ONE),
    ])
}

pub fn matrix_sort_key(matrix: &MatrixExpr) -> Key {
    SortKeyCache { memo_scalars: true, ..Default::default() }.matrix_key(matrix)
}

/// A Python-tuple like sort key.
#[derive(Clone, Debug)]
pub enum Key {
    Int(i64),
    Str(String),
    Num(Number),
    Tup(Vec<Key>),
    Shared(Arc<[Key]>),
}

impl Key {
    fn rank(&self) -> u8 {
        match self {
            Key::Int(_) => 0,
            Key::Num(_) => 0,
            Key::Str(_) => 1,
            Key::Tup(_) | Key::Shared(_) => 2,
        }
    }

    fn tuple(&self) -> Option<&[Key]> {
        match self { Key::Tup(v) => Some(v), Key::Shared(v) => Some(v), _ => None }
    }
}

fn cmp_num(a: &Number, b: &Number) -> Ordering {
    a.cmp_value(b).unwrap_or_else(|| a.cmp_struct(b))
}

impl Ord for Key {
    fn cmp(&self, o: &Key) -> Ordering {
        if let (Some(a), Some(b)) = (self.tuple(), o.tuple()) {
            if a.len() == b.len() && std::ptr::eq(a.as_ptr(), b.as_ptr()) { return Ordering::Equal; }
            for (x, y) in a.iter().zip(b.iter()) {
                let comparison = x.cmp(y);
                if comparison != Ordering::Equal { return comparison; }
            }
            return a.len().cmp(&b.len());
        }
        match (self, o) {
            (Key::Int(a), Key::Int(b)) => a.cmp(b),
            (Key::Num(a), Key::Num(b)) => cmp_num(a, b),
            (Key::Int(a), Key::Num(b)) => cmp_num(&Number::Int(*a), b),
            (Key::Num(a), Key::Int(b)) => cmp_num(a, &Number::Int(*b)),
            (Key::Str(a), Key::Str(b)) => a.cmp(b),
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

fn share_key(key: Key) -> Key {
    match key {
        Key::Tup(children) => Key::Shared(children.into_iter().map(share_key).collect::<Vec<_>>().into()),
        other => other,
    }
}

pub(crate) fn undefined_ordinal(name: &str) -> u32 {
    match name {
        "exp" => 10, "log" => 11,
        "sin" => 20, "cos" => 21, "tan" => 22, "cot" => 23,
        "sinh" => 30, "cosh" => 31, "tanh" => 32, "coth" => 33,
        "conjugate" => 40, "re" => 41, "im" => 42, "arg" => 43,
        _ => 0,
    }
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
        Kind::Func(f @ (Func::Logic(_) | Func::Set(_) | Func::Subs | Func::Limit), _) => (5, 0, f.name().into()),
        Kind::Func(Func::Undefined(head), _) if head.internal => (5, 0, head.name.to_string()),
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
                Func::Undefined(head) => undefined_ordinal(&head.name).into(),
                Func::Max | Func::Min => 0,
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
    SortKeyCache { memo_scalars: e.has_matrix(), ..Default::default() }.key(e)
}

/// Each scalar and matrix node is keyed once; tuple children retain sharing.
#[derive(Default)]
struct SortKeyCache {
    memo_scalars: bool,
    // Presentation depends on the stored bound names, even when two Subs
    // nodes have equal alpha identities. Retain owners to prevent ID reuse.
    scalars: FxHashMap<usize, (Expr, Key)>,
    matrices: FxHashMap<usize, (MatrixExpr, Key)>,
}

impl SortKeyCache {
    fn key(&mut self, e: &Expr) -> Key {
        if !self.memo_scalars { return self.compute(e); }
        if let Some((_, key)) = self.scalars.get(&e.dag_id()) { return key.clone(); }
        let key = share_key(self.compute(e));
        self.scalars.insert(e.dag_id(), (e.clone(), key.clone()));
        key
    }

    fn presentation_key(&mut self, arg: PresentationArg) -> Key {
        match arg {
            PresentationArg::Scalar(e) => self.key(&e),
            PresentationArg::Matrix(m) => self.matrix_key(&m),
            PresentationArg::Text(s) => composite_key((2, 0, "Str"), vec![Key::Str(s)]),
            PresentationArg::Node(name, args) => {
                let args = args.into_iter().map(|arg| self.presentation_key(arg)).collect();
                composite_key((5, 0, name), args)
            }
        }
    }

    fn matrix_key(&mut self, matrix: &MatrixExpr) -> Key {
        if let Some((_, key)) = self.matrices.get(&matrix.dag_id()) { return key.clone(); }
        let mut presentation = matrix_presentation(matrix);
        if let MatrixKind::Add(terms) = matrix.kind() {
            presentation.args = matrix_ordered_terms(terms, |m| self.matrix_key(m)).into_iter().map(PresentationArg::Matrix).collect();
        }
        let args = presentation.args.into_iter().map(|arg| self.presentation_key(arg)).collect();
        let key = share_key(composite_key(presentation.class, args));
        self.matrices.insert(matrix.dag_id(), (matrix.clone(), key.clone()));
        key
    }

    fn compute(&mut self, e: &Expr) -> Key {
    match e.kind() {
        Kind::Derivative(expression, variables) => {
            let mut args = vec![self.key(expression)];
            for (variable, count) in variables {
                let args_variable = vec![self.key(variable), number_key(Number::Int(*count as i64))];
                args.push(composite_key((5, 0, "Tuple"), args_variable));
            }
            return composite_key((5, 0, "Derivative"), args);
        }
        Kind::Func(Func::MatrixProxy(matrix), _) => return self.matrix_key(matrix),
        Kind::Func(Func::MatrixElement(matrix), args) => {
            let mut keys = vec![self.matrix_key(matrix)];
            keys.extend(args.iter().map(|arg| self.key(arg)));
            return composite_key((5, 0, "MatrixElement"), keys);
        }
        Kind::Func(f @ (Func::MatrixTrace(matrix) | Func::MatrixDeterminant(matrix)), _) => {
            return composite_key((5, 0, f.name()), vec![self.matrix_key(matrix)]);
        }
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
    if matches!(base.kind(), Kind::Derivative(..) | Kind::Func(Func::MatrixProxy(_) | Func::MatrixElement(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_), _)) {
        let mut parts = match self.key(&base) {
            Key::Tup(parts) => parts,
            Key::Shared(parts) => parts.to_vec(),
            _ => unreachable!("expression sort keys are tuples"),
        };
        parts[2] = self.key(&exp);
        parts[3] = Key::Num(coeff);
        return tup(parts);
    }
    let args: Vec<Key> = if base.is_atom() {
        vec![Key::Str(crate::printer::to_str(&base))]
    } else {
        match base.kind() {
            Kind::Add(..) => ordered_terms(&base).iter().map(|arg| self.key(arg)).collect(),
            Kind::Mul(..) => ordered_factors(&base).iter().map(|arg| self.key(arg)).collect(),
            _ => base.args().iter().map(|arg| self.key(arg)).collect(),
        }
    };
    tup(vec![
        class_key(&base),
        tup(vec![Key::Int(args.len() as i64), tup(args)]),
        self.key(&exp),
        Key::Num(coeff),
    ])
    }
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

/// Compare expressions using SymPy default ordering.
pub fn cmp_default(a: &Expr, b: &Expr) -> Ordering {
    sort_key(a).cmp(&sort_key(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_and_owned_sort_keys_have_identical_ordering() {
        let keys: Vec<_> = ["0", "-1", "1/2", "x", "y", "x**2", "sin(x)", "x+y", "x*y"]
            .iter().map(|text| sort_key(&crate::parser::parse_str(text).unwrap())).collect();
        for a in &keys {
            for b in &keys {
                let expected = a.cmp(b);
                assert_eq!(a.cmp(&share_key(b.clone())), expected);
                assert_eq!(share_key(a.clone()).cmp(b), expected);
                assert_eq!(share_key(a.clone()).cmp(&share_key(b.clone())), expected);
            }
        }
    }
}
