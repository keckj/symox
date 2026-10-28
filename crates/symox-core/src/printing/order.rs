//! SymPy's printing order: `default_sort_key`, `Mul.as_ordered_factors`
//! and `Add.as_ordered_terms` (lexicographic order of the monomials over
//! the generators of the sum).
//!
//! Used by every printer (including `str`) so that terms and factors come
//! out in the same order as in SymPy.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::rc::Rc;

use rustc_hash::FxHashMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add_term, mul, mul_factor, pow};
use crate::matrixexpr::{MatrixExpr, MatrixKind};
use crate::sortkey::{matrix_ordered_terms, matrix_presentation, Key, PresentationArg};

// ----------------------------------------------------------------------
// default_sort_key
// ----------------------------------------------------------------------

#[derive(Clone, Debug)]
enum KArgs {
    /// `(0, ())` of numbers
    Empty,
    /// `(1, (str(atom),))`
    Atom(String),
    /// Dummy atoms include their identity index after the unprefixed name.
    Dummy(String, num_bigint::BigInt),
    /// `(len(args), tuple(keys))`
    List(Vec<Rc<SortKey>>),
}

/// Emulation of SymPy's `default_sort_key`.
#[derive(Clone, Debug)]
pub struct SortKey {
    class: (u8, u32, Cow<'static, str>),
    args: KArgs,
    /// `None` for numbers (their key has an empty exponent part).
    exp: Option<Rc<SortKey>>,
    coeff: Key,
}

fn cmp_args(a: &KArgs, b: &KArgs) -> Ordering {
    let len = |k: &KArgs| match k {
        KArgs::Empty => 0,
        KArgs::Atom(_) => 1,
        KArgs::Dummy(..) => 2,
        KArgs::List(v) => v.len(),
    };
    len(a).cmp(&len(b)).then_with(|| match (a, b) {
        (KArgs::Empty, KArgs::Empty) => Ordering::Equal,
        (KArgs::Atom(x), KArgs::Atom(y)) => x.cmp(y),
        (KArgs::Dummy(x, i), KArgs::Dummy(y, j)) => x.cmp(y).then_with(|| i.cmp(j)),
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
        (KArgs::Dummy(..), _) => Ordering::Less,
        (_, KArgs::Dummy(..)) => Ordering::Greater,
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
        if std::ptr::eq(self, o) { return Ordering::Equal; }
        self.class
            .cmp(&o.class)
            .then_with(|| cmp_args(&self.args, &o.args))
            .then_with(|| match (&self.exp, &o.exp) {
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
                (Some(a), Some(b)) => a.cmp(b),
            })
            .then_with(|| self.coeff.cmp(&o.coeff))
    }
}

/// `Function.class_key` ordinals of the known functions.
fn func_ordinal(f: &Func) -> u32 {
    match f {
        Func::Undefined(head) => crate::sortkey::undefined_ordinal(head),
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
        // Classes without `nargs` take any number of arguments (`Naturals0`).
        Func::ExpPolar => 0,
        Func::Special(s) if s.is_transform() || matches!(s, crate::special::SpecialFn::MeijerG
            | crate::special::SpecialFn::Lerchphi | crate::special::SpecialFn::LeviCivita
            | crate::special::SpecialFn::Carmichael) => 0,
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
            Constant::Catalan => (2, 0, Cow::Borrowed("Catalan")),
            Constant::GoldenRatio => (2, 0, Cow::Borrowed("GoldenRatio")),
            Constant::TribonacciConstant => (2, 0, Cow::Borrowed("TribonacciConstant")),
            Constant::ComplexInfinity => (2, 0, Cow::Borrowed("ComplexInfinity")),
        },
        Kind::Sym(_) => (2, 0, Cow::Owned(crate::host::class_name(e).unwrap_or(e.class_name()).into())),
        Kind::Mul(..) => (3, 0, Cow::Borrowed("Mul")),
        Kind::Add(..) => (3, 1, Cow::Borrowed("Add")),
        Kind::Pow(..) => (3, 2, Cow::Borrowed("Pow")),
        Kind::Func(Func::Max, _) => (5, 0, Cow::Borrowed("Max")),
        Kind::Func(Func::Min, _) => (5, 0, Cow::Borrowed("Min")),
        Kind::Func(f @ (Func::SymbolicDerivative | Func::Logic(_) | Func::Set(_) | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::Order | Func::MatrixElement(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..)), _) => (5, 0, Cow::Owned(f.name().into())),
        Kind::Func(Func::Undefined(head), args) if crate::poly::rootof::is_rootof_head(head, args.len()) => (5, 0, Cow::Borrowed("ComplexRootOf")),
        Kind::Func(Func::Undefined(head), _) if head.internal => (5, 0, Cow::Owned(head.name.to_string())),
        Kind::Func(Func::Undefined(head), _) if head.host().and_then(|host| host.class_key()).is_some() => {
            let (a, b) = head.host().and_then(|host| host.class_key()).unwrap();
            (a, b, Cow::Owned(head.name.to_string()))
        }
        Kind::Func(f, _) => (4, func_ordinal(f), Cow::Owned(f.name().to_string())),
        Kind::Piecewise(_) => (4, 0, Cow::Borrowed("Piecewise")),
        Kind::Derivative(..) => (5, 0, Cow::Borrowed("Derivative")),
        Kind::Rel(..) => (5, 0, Cow::Owned(e.class_name().to_string())),
        Kind::Bool(_) => (5, 0, Cow::Owned(e.class_name().to_string())),
    }
}

/// `as_coeff_Mul()` with the coefficient as a sort-key component (a real
/// infinity or nan factor is the coefficient, as SymPy's Number).
fn coeff_mul(e: &Expr) -> (Key, Expr) {
    if let Kind::Mul(c, facs) = e.kind() {
        // A held product: as_coeff_Mul of its own first argument.
        if !e.is_canonical() || c.is_one() && facs.iter().any(|(b, x)| x.is_one() && matches!(b.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN))) {
            let (c, rest) = crate::decompose::as_coeff_mul_head(e, false);
            return match c.kind() {
                Kind::Num(n) => (Key::Num(n.clone()), rest),
                Kind::Const(k) => (Key::SpecialNumber(*k), rest),
                _ => (Key::Num(Number::ONE), e.clone()),
            };
        }
    }
    let (c, rest) = crate::ops::as_coeff_mul(e);
    (Key::Num(c), rest)
}

fn one_key() -> SortKey {
    SortKey { class: (1, 0, Cow::Borrowed("Number")), args: KArgs::Empty, exp: None, coeff: Key::Num(Number::ONE) }
}

fn atom_args(e: &Expr) -> KArgs {
    if let Kind::Sym(symbol) = e.kind() {
        let name = crate::host::atom_name(e, symbol);
        if let Some(index) = symbol.dummy_index() {
            return KArgs::Dummy(name.into(), index.clone());
        }
        return KArgs::Atom(name.into());
    }
    KArgs::Atom(crate::printer::to_str(e))
}

/// Memo table of sort keys (keys of shared sub-expressions are computed
/// once).
#[derive(Default)]
pub struct KeyCache {
    // Sort keys describe original args, not alpha-equivalence classes.
    // Owners keep temporary args alive while their physical IDs are cached.
    map: FxHashMap<usize, (Expr, Rc<SortKey>)>,
    matrices: FxHashMap<usize, (MatrixExpr, Rc<SortKey>)>,
    one: Option<Rc<SortKey>>,
}

impl KeyCache {
    fn one(&mut self) -> Rc<SortKey> {
        self.one.get_or_insert_with(|| Rc::new(one_key())).clone()
    }

    fn composite(&mut self, class: (u8, u32, &'static str), args: KArgs) -> SortKey {
        SortKey { class: (class.0, class.1, Cow::Borrowed(class.2)), args, exp: Some(self.one()), coeff: Key::Num(Number::ONE) }
    }

    fn presentation_arg(&mut self, arg: PresentationArg) -> Rc<SortKey> {
        match arg {
            PresentationArg::Scalar(e) => self.key(&e),
            PresentationArg::Matrix(m) => self.matrix_key(&m),
            PresentationArg::Text(s) => Rc::new(self.composite((2, 0, "Str"), KArgs::Atom(s))),
            PresentationArg::Permutation(s) => Rc::new(self.composite((2, 0, "Permutation"), KArgs::Atom(s))),
            PresentationArg::None => Rc::new(self.composite((0, 0, "NoneType"), KArgs::Atom("None".into()))),
            PresentationArg::Node(name, args) => {
                let args = args.into_iter().map(|arg| self.presentation_arg(arg)).collect();
                Rc::new(self.composite((5, 0, name), KArgs::List(args)))
            }
        }
    }

    fn matrix_key(&mut self, matrix: &MatrixExpr) -> Rc<SortKey> {
        if let Some((_, key)) = self.matrices.get(&matrix.dag_id()) { return key.clone(); }
        let mut presentation = matrix_presentation(matrix);
        if let MatrixKind::Add(terms) = matrix.kind() {
            presentation.args = matrix_ordered_terms(terms, |m| self.matrix_key(m)).into_iter().map(PresentationArg::Matrix).collect();
        }
        let args: Vec<_> = presentation.args.into_iter().map(|arg| self.presentation_arg(arg)).collect();
        let key = Rc::new(self.composite(presentation.class, KArgs::List(args)));
        self.matrices.insert(matrix.dag_id(), (matrix.clone(), key.clone()));
        key
    }

    /// `default_sort_key(e)`.
    pub fn key(&mut self, e: &Expr) -> Rc<SortKey> {
        if let Some((_, k)) = self.map.get(&e.dag_id()) {
            return k.clone();
        }
        let k = Rc::new(self.compute(e));
        self.map.insert(e.dag_id(), (e.clone(), k.clone()));
        k
    }

    fn compute(&mut self, e: &Expr) -> SortKey {
        match e.kind() {
            // Heaviside(x).args is (x, 1/2) in SymPy
            Kind::Func(Func::Heaviside,args) if args.len()==1 => {
                return (*self.key(&Expr::from_kind(Kind::Func(Func::Heaviside,vec![args[0].clone(),Expr::rational(1,2)])))).clone();
            }
            Kind::Func(Func::Set(crate::sets::SetOp::Lambda),args) if args.len()==2 => {
                let signature=crate::sets::lambda_signature(e).unwrap();
                let variables=self.key(&signature);
                let body=self.key(&args[1]);
                return self.composite((5,0,"Lambda"),KArgs::List(vec![variables,body]));
            }
            Kind::Func(Func::Set(crate::sets::SetOp::IdentityFunction),_)=>return self.composite((5,0,"IdentityFunction"),KArgs::List(vec![])),
            Kind::Func(Func::SymbolicDerivative,args)=>{let keys=args.iter().map(|arg|self.key(arg)).collect();return self.composite((5,0,"Derivative"),KArgs::List(keys));},
            Kind::Derivative(expression, variables) => {
                let mut args = vec![self.key(expression)];
                for (variable, count) in variables {
                    args.push(self.presentation_arg(PresentationArg::Node("Tuple", vec![PresentationArg::Scalar(variable.clone()), PresentationArg::Scalar(Expr::int(*count as i64))])));
                }
                return self.composite((5, 0, "Derivative"), KArgs::List(args));
            }
            // Order(expr, (x, 0), ...): the point pairs are Tuples.
            Kind::Func(Func::Order, arguments) if !arguments.is_empty() => {
                let mut args = vec![self.key(&arguments[0])];
                for pair in arguments[1..].chunks(2) {
                    args.push(self.presentation_arg(PresentationArg::Node("Tuple", pair.iter().cloned().map(PresentationArg::Scalar).collect())));
                }
                return self.composite((5, 0, "Order"), KArgs::List(args));
            }
            Kind::Piecewise(pieces) => {
                let args = pieces.iter().map(|(expression, condition)| {
                    self.presentation_arg(PresentationArg::Node("ExprCondPair", vec![PresentationArg::Scalar(expression.clone()), PresentationArg::Scalar(condition.clone())]))
                }).collect();
                return self.composite((4, 0, "Piecewise"), KArgs::List(args));
            }
            // BooleanTrue/BooleanFalse are Basic atoms without arguments.
            Kind::Bool(_) => return SortKey { class: class_key(e), args: KArgs::List(vec![]), exp: Some(self.one()), coeff: Key::Num(Number::ONE) },
            Kind::Func(Func::MatrixProxy(matrix), _) => return (*self.matrix_key(matrix)).clone(),
            Kind::Func(Func::MatrixDotProduct(left,right),_) => {
                let args=vec![self.matrix_key(left),self.matrix_key(right)];
                return self.composite((5,0,"DotProduct"),KArgs::List(args));
            }
            Kind::Func(Func::MatrixElement(matrix), args) => {
                let mut keys = vec![self.matrix_key(matrix)];
                keys.extend(args.iter().map(|arg| self.key(arg)));
                return self.composite((5, 0, "MatrixElement"), KArgs::List(keys));
            }
            Kind::Func(Func::MatrixTrace(matrix), _) => {
                let args = vec![self.matrix_key(matrix)];
                return self.composite((5, 0, "Trace"), KArgs::List(args));
            }
            Kind::Func(Func::MatrixDeterminant(matrix), _) => {
                let args = vec![self.matrix_key(matrix)];
                return self.composite((5, 0, "Determinant"), KArgs::List(args));
            }
            // Number.sort_key: (class_key, (0, ()), (), self)
            Kind::Num(n) => return SortKey { class: class_key(e), args: KArgs::Empty, exp: None, coeff: Key::Num(n.clone()) },
            Kind::Const(c @ (Constant::Infinity | Constant::NegInfinity | Constant::NaN)) => {
                return SortKey { class: class_key(e), args: KArgs::Empty, exp: None, coeff: Key::SpecialNumber(*c) };
            }
            // Atom.sort_key: (class_key, (1, (str(self),)), One.sort_key(), One)
            Kind::Sym(_) | Kind::Const(_) => return SortKey {
                class: class_key(e),
                args: atom_args(e),
                exp: Some(self.one()),
                coeff: Key::Num(Number::ONE),
            },
            _ => {}
        }
        // Expr.sort_key: coeff, expr = as_coeff_Mul(); base, exp = as_base_exp()
        let (coeff, rest) = coeff_mul(e);
        let (expr, exp) = match rest.kind() {
            Kind::Pow(..) => {
                let (b, x) = crate::decompose::as_base_exp(&rest);
                // E**x sorts like the function exp(x)
                if b.is_constant(Constant::E) { (crate::functions::exp(&x), None) } else { (b, Some(x)) }
            }
            _ => (rest, None),
        };
        let heaviside = matches!(expr.kind(), Kind::Func(Func::Heaviside, args) if args.len() == 1);
        if heaviside || matches!(expr.kind(), Kind::Derivative(..) | Kind::Func(Func::SymbolicDerivative,_) | Kind::Piecewise(_) | Kind::Func(Func::Order | Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction),_) | Kind::Func(Func::MatrixProxy(_) | Func::MatrixElement(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..), _)) {
            let mut key = (*self.key(&expr)).clone();
            key.exp = Some(match exp { Some(x) => self.key(&x), None => self.one() });
            key.coeff = coeff;
            return key;
        }
        let args = match expr.kind() {
            Kind::Sym(symbol) if symbol.is_dummy() => KArgs::List(vec![self.key(&expr)]),
            Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) => atom_args(&expr),
            Kind::Bool(_) => KArgs::Atom(crate::printer::to_str(&expr)),
            _ if crate::host::expression_head(&expr).is_some_and(|host| host.atomic()) => {
                KArgs::Atom(crate::printer::to_str(&expr))
            }
            Kind::Add(..) => {
                let t = if expr.is_canonical() { self.ordered_terms(&expr) } else { crate::sortkey::held_ordered_terms(&expr) };
                KArgs::List(t.iter().map(|a| self.key(a)).collect())
            }
            Kind::Mul(..) => {
                let f = self.ordered_factors(&expr);
                KArgs::List(f.iter().map(|a| self.key(a)).collect())
            }
            _ => KArgs::List(expr.args().iter().map(|a| self.key(a)).collect()),
        };
        let exp = match exp {
            Some(x) => self.key(&x),
            None => self.one(),
        };
        SortKey { class: class_key(&expr), args, exp: Some(exp), coeff }
    }

    /// Sort factors with `default_sort_key` (`Mul.as_ordered_factors`).
    pub fn sort_factors(&mut self, v: &mut [Expr]) {
        if v.len() < 2 {
            return;
        }
        if v.iter().all(|e| simple_factor(e).is_some()) {
            v.sort_by(|a, b| cmp_simple(&simple_factor(a).unwrap(), &simple_factor(b).unwrap()));
            return;
        }
        let mut keyed: Vec<(Rc<SortKey>, Expr)> = v.iter().map(|e| (self.key(e), e.clone())).collect();
        keyed.sort_by(|a, b| a.0.cmp(&b.0));
        for (slot, (_, e)) in v.iter_mut().zip(keyed) {
            *slot = e;
        }
    }

    /// `Mul.as_ordered_factors()`.
    pub fn ordered_factors(&mut self, e: &Expr) -> Vec<Expr> {
        if !e.is_mul(){return vec![e.clone()];}
        let(mut commutative,noncommutative)=crate::decompose::args_cnc(e,true);
        self.sort_factors(&mut commutative);
        commutative.extend(noncommutative);commutative
    }
}

/// Cheap view of the sort key of a number or of `atom**number`.
struct Simple<'a> {
    class: (u8, &'static str),
    name: &'a str,
    exp: &'a Number,
}

fn atom_class(e: &Expr) -> Option<(u8, &'static str, &str)> {
    match e.kind() {
        Kind::Sym(s) if !s.is_dummy() && s.host().is_none() => Some((2, if s.is_wild() { "Wild" } else { "Symbol" }, s)),
        Kind::Const(Constant::Pi) => Some((2, "Pi", "pi")),
        Kind::Const(Constant::E) => Some((2, "Exp1", "E")),
        Kind::Const(Constant::I) => Some((2, "ImaginaryUnit", "I")),
        Kind::Const(Constant::EulerGamma) => Some((2, "EulerGamma", "EulerGamma")),
        Kind::Const(Constant::Catalan) => Some((2, "Catalan", "Catalan")),
        Kind::Const(Constant::GoldenRatio) => Some((2, "GoldenRatio", "GoldenRatio")),
        Kind::Const(Constant::TribonacciConstant) => Some((2, "TribonacciConstant", "TribonacciConstant")),
        _ => None,
    }
}

fn simple_factor(e: &Expr) -> Option<Simple<'_>> {
    match e.kind() {
        Kind::Num(n) => Some(Simple { class: (1, "Number"), name: "", exp: n }),
        Kind::Pow(b, x) => {
            let (c, k, name) = atom_class(b)?;
            Some(Simple { class: (c, k), name, exp: x.as_number()? })
        }
        _ => {
            let (c, k, name) = atom_class(e)?;
            Some(Simple { class: (c, k), name, exp: &Number::ONE })
        }
    }
}

fn cmp_simple(a: &Simple, b: &Simple) -> Ordering {
    a.class
        .cmp(&b.class)
        .then_with(|| a.name.cmp(b.name))
        .then_with(|| a.exp.cmp_value(b.exp).unwrap_or(Ordering::Equal))
}

thread_local! {
    static SESSION: std::cell::RefCell<Option<KeyCache>> = const { std::cell::RefCell::new(None) };
}

/// Keeps a sort-key cache alive for the current thread while it exists:
/// printing a whole expression computes the key of every sub-expression
/// once instead of once per enclosing sum (linear instead of quadratic in
/// the depth).
pub struct Session {
    owner: bool,
}

/// Start a printing session (nested sessions share the outer cache).
pub fn session() -> Session {
    SESSION.with(|c| match c.try_borrow_mut() {
        Ok(mut slot) if slot.is_none() => {
            *slot = Some(KeyCache::default());
            Session { owner: true }
        }
        _ => Session { owner: false },
    })
}

impl Drop for Session {
    fn drop(&mut self) {
        if self.owner {
            SESSION.with(|c| {
                if let Ok(mut slot) = c.try_borrow_mut() {
                    *slot = None;
                }
            });
        }
    }
}

/// Run `f` with the session cache when there is one (and it is not in
/// use), with a fresh cache otherwise.
fn with_cache<R>(f: impl FnOnce(&mut KeyCache) -> R) -> R {
    let mut f = Some(f);
    let r = SESSION.with(|c| match c.try_borrow_mut() {
        Ok(mut slot) => slot.as_mut().map(|k| (f.take().unwrap())(k)),
        Err(_) => None,
    });
    match r {
        Some(r) => r,
        None => (f.take().unwrap())(&mut KeyCache::default()),
    }
}

/// SymPy's `default_sort_key(e)`.
pub fn sort_key(e: &Expr) -> Rc<SortKey> {
    with_cache(|k| k.key(e))
}

// ----------------------------------------------------------------------
// factors
// ----------------------------------------------------------------------

/// Sort factors with `default_sort_key` (`Mul.as_ordered_factors`).
/// `sort_factors` of the commutative factors only: factors that may not
/// commute keep their relative order after them, as in
/// `Mul.as_ordered_factors`. Without such a factor this is `sort_factors`.
pub fn sort_commutative_factors(v: &mut Vec<Expr>) {
    if !v.iter().any(Expr::has_noncommutative_host) { return sort_factors(v); }
    let (mut commutative, noncommutative): (Vec<Expr>, Vec<Expr>) =
        v.drain(..).partition(crate::decompose::commutes);
    sort_factors(&mut commutative);
    commutative.extend(noncommutative);
    *v = commutative;
}

pub fn sort_factors(v: &mut [Expr]) {
    if v.len() < 2 {
        return;
    }
    with_cache(|k| k.sort_factors(v))
}

/// `Mul.as_ordered_factors()`.
pub fn ordered_factors(e: &Expr) -> Vec<Expr> {
    with_cache(|k| k.ordered_factors(e))
}

// ----------------------------------------------------------------------
// terms
// ----------------------------------------------------------------------

/// `decompose_power`: symbolic base and integer exponent.
pub(crate) fn decompose_power(f: &Expr) -> (Expr, i64) {
    // base, exp = expr.as_base_exp()
    let (b, x) = match f.kind() {
        Kind::Pow(..) | Kind::Func(Func::Exp, _) => crate::decompose::as_base_exp(f),
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
        || matches!(e.kind(), Kind::Const(c) if c.is_number_symbol())
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

/// `MatAdd.as_ordered_terms()`: the terms of a matrix sum in SymPy's
/// printing order (used by the LaTeX printer).
pub fn matrix_ordered(terms: &[MatrixExpr]) -> Vec<MatrixExpr> {
    with_cache(|k| matrix_ordered_terms(terms, |m| k.matrix_key(m)))
}

/// `Add.as_ordered_terms()`: the terms of a sum (with their coefficients,
/// the numeric term included) in SymPy's printing order.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    if !e.is_add() {
        return vec![e.clone()];
    }
    if !e.is_canonical() {
        return crate::sortkey::held_ordered_terms(e);
    }
    with_cache(|k| k.ordered_terms(e))
}

/// `Add.as_ordered_terms()` special case: a positive Number or NumberSymbol
/// and a negative Number times one factor keep that order ("1 - x").
fn positive_constant_first(all: &[Expr]) -> bool {
    if all.len() != 2 { return false; }
    let numeric = |a: &Expr| matches!(a.kind(), Kind::Num(_) | Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN))
        || matches!(a.kind(), Kind::Const(c) if c.is_number_symbol());
    let (n, m) = if numeric(&all[0]) { (0, 1) } else { (1, 0) };
    if !numeric(&all[n]) || numeric(&all[m]) { return false; }
    // `oo.is_positive` is False (it is not finite).
    let positive = match all[n].kind() { Kind::Num(x) => x.is_positive(), Kind::Const(c) => c.is_number_symbol(), _ => false };
    positive && matches!(all[m].kind(), Kind::Mul(k, facs) if facs.len() == 1 && k.is_negative())
}

/// One term of `Expr.as_terms()`: complex coefficient, commutative
/// generator powers and the factors that may not commute.
struct Term {
    expr: Expr,
    re: f64,
    im: f64,
    cpart: Vec<(usize, i64)>,
    ncpart: Vec<Rc<SortKey>>,
    order: bool,
}

/// `((bool(im), im), (re, im))` with Python tuple semantics: the first
/// unequal component decides, and an unordered (nan) pair is "not less".
fn cmp_coeff(a: &Term, b: &Term) -> Ordering {
    let (ai, bi) = (a.im != 0.0, b.im != 0.0);
    if ai != bi { return ai.cmp(&bi); }
    for (x, y) in [(a.im, b.im), (a.re, b.re), (a.im, b.im)] {
        if x != y { return x.partial_cmp(&y).unwrap_or(Ordering::Equal); }
    }
    Ordering::Equal
}

fn gen_index(base: Expr, gens: &mut FxHashMap<Expr, usize>, gen_list: &mut Vec<Expr>) -> usize {
    if let Some(&i) = gens.get(&base) {
        return i;
    }
    gens.insert(base.clone(), gen_list.len());
    gen_list.push(base);
    gen_list.len() - 1
}

impl KeyCache {
    /// One factor of a term in `Expr.as_terms()`.
    #[allow(clippy::too_many_arguments)]
    fn term_factor(&mut self, factor: Expr, noncommutative: bool, z: &mut num_complex::Complex64, cpart: &mut Vec<(usize, i64)>,
        ncpart: &mut Vec<Rc<SortKey>>, gens: &mut FxHashMap<Expr, usize>, gen_list: &mut Vec<Expr>) {
        // if factor.is_number: coeff *= complex(factor); complex(zoo) is nan+nanj
        if factor.is_constant(Constant::ComplexInfinity) {
            *z *= num_complex::Complex64::new(f64::NAN, f64::NAN);
            return;
        }
        // A host NumberSymbol (hbar): complex(factor) through its evalf.
        if let Kind::Func(Func::Undefined(head), args) = factor.kind() {
            if let Some(host) = head.host().filter(|host| args.is_empty() && host.atomic() && host.number() && host.evaluates_numerically()) {
                if let Some(v) = host.evalf(&factor, 53).as_ref().and_then(crate::evalf::eval_c64) {
                    *z *= v;
                    return;
                }
            }
        }
        // `factor.is_number`: a closed Integral/Sum/Product binds its
        // variables, so the free-symbol mask alone does not decide it.
        if (super::is_numeric(&factor) || matches!(factor.kind(), Kind::Func(Func::Integral | Func::Sum | Func::Product, _)))
            && crate::symbolic_number::is_symbolic_number(&factor) {
            if let Some(v) = crate::sortkey::numeric_coefficient(&factor) {
                *z *= v;
                return;
            }
        }
        // `factor.is_commutative` is None for heads such as Limit, whose
        // direction argument is a Basic Symbol: they go in ncpart even when the
        // sum has no noncommutative symbol.
        if (noncommutative || matches!(factor.kind(), Kind::Func(..))) && crate::decompose::commutativity(&factor) != Some(true) {
            ncpart.push(self.key(&factor));
            return;
        }
        let (base, d) = decompose_power(&factor);
        let g = gen_index(base, gens, gen_list);
        // cpart is a dict: a repeated base keeps the last exponent
        cpart.retain(|(h, _)| *h != g);
        cpart.push((g, d));
    }

    /// `Add.as_ordered_terms()`.
    pub fn ordered_terms(&mut self, e: &Expr) -> Vec<Expr> {
        let Kind::Add(c, terms) = e.kind() else {
            return vec![e.clone()];
        };
        // Add.make_args order: the Number first
        let mut all: Vec<Expr> = Vec::with_capacity(terms.len() + 1);
        if !c.is_zero() {
            all.push(Expr::num(c.clone()));
        }
        all.extend(terms.iter().map(|(t, k)| add_term(t, k)));
        if positive_constant_first(&all) {
            if !is_number_like(&all[0]) { all.swap(0, 1); }
            return all;
        }
        let noncommutative = e.has_noncommutative_host();
        let mut gens: FxHashMap<Expr, usize> = FxHashMap::default();
        let mut gen_list: Vec<Expr> = Vec::new();
        let mut decomposed: Vec<Term> = Vec::with_capacity(all.len());
        for t in all {
            let order = t.is_func(&Func::Order);
            // coeff, _term = term.as_coeff_Mul(); coeff = complex(coeff)
            let (coeff, facs): (Number, &[(Expr, Expr)]) = match t.kind() {
                Kind::Num(n) => (n.clone(), &[]),
                Kind::Mul(c, facs) => (c.clone(), facs.as_slice()),
                _ => (Number::ONE, &[]),
            };
            let mut z = num_complex::Complex64::new(coeff.to_f64(), 0.0);
            let mut cpart: Vec<(usize, i64)> = Vec::with_capacity(facs.len().max(1));
            let mut ncpart = Vec::new();
            // A real infinity (or nan) factor is SymPy's Number coefficient.
            let special = facs.iter().position(|(b, x)| x.is_one() && matches!(b.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN)));
            if let Some(i) = special {
                let value = match facs[i].0.kind() {
                    Kind::Const(Constant::Infinity) => f64::INFINITY,
                    Kind::Const(Constant::NegInfinity) => f64::NEG_INFINITY,
                    _ => f64::NAN,
                };
                z = num_complex::Complex64::new(z.re * value, 0.0);
            }
            let factors: Vec<Expr> = match t.kind() {
                Kind::Num(_) => Vec::new(),
                Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN) => {
                    z = num_complex::Complex64::new(if t.is_constant(Constant::Infinity) { f64::INFINITY } else if t.is_constant(Constant::NegInfinity) { f64::NEG_INFINITY } else { f64::NAN }, 0.0);
                    Vec::new()
                }
                Kind::Mul(..) => Vec::new(),
                _ => vec![t.clone()],
            };
            for factor in factors {
                self.term_factor(factor, noncommutative, &mut z, &mut cpart, &mut ncpart, &mut gens, &mut gen_list);
            }
            for (i, (b, x)) in facs.iter().enumerate() {
                if Some(i) == special { continue; }
                // fast path: symbol ** integer
                if let (false, Kind::Sym(_), Kind::Num(Number::Int(d))) = (noncommutative, b.kind(), x.kind()) {
                    let g = gen_index(b.clone(), &mut gens, &mut gen_list);
                    cpart.retain(|(h, _)| *h != g);
                    cpart.push((g, *d));
                    continue;
                }
                self.term_factor(mul_factor(b, x), noncommutative, &mut z, &mut cpart, &mut ncpart, &mut gens, &mut gen_list);
            }
            decomposed.push(Term { expr: t, re: z.re, im: z.im, cpart, ncpart, order });
        }
        // gens = sorted(gens, key=default_sort_key)
        let mut rank: Vec<usize> = (0..gen_list.len()).collect();
        if gen_list.len() > 1 {
            let keys: Vec<Rc<SortKey>> = gen_list.iter().map(|g| self.key(g)).collect();
            let mut order: Vec<usize> = (0..gen_list.len()).collect();
            order.sort_by(|a, b| keys[*a].cmp(&keys[*b]));
            for (r, g) in order.iter().enumerate() {
                rank[*g] = r;
            }
        }
        let mut keyed: Vec<(Vec<(usize, i64)>, Term)> = decomposed
            .into_iter()
            .map(|t| {
                let mut m: Vec<(usize, i64)> = t.cpart.iter().map(|(g, d)| (rank[*g], *d)).collect();
                m.sort_unstable();
                (m, t)
            })
            .collect();
        // key(term) = (neg(monom), ncpart, coeff): descending lex monomials
        let compare = |(ma, ta): &(Vec<(usize, i64)>, Term), (mb, tb): &(Vec<(usize, i64)>, Term)| {
            cmp_lex(mb, ma).then_with(|| ta.ncpart.cmp(&tb.ncpart)).then_with(|| cmp_coeff(ta, tb))
        };
        let has_order = keyed.iter().any(|(_, t)| t.order);
        let has_nan = keyed.iter().any(|(_, t)| t.re.is_nan() || t.im.is_nan());
        // Python's stable sort; with an Order term both groups are sorted
        // with reverse=True (equal keys keep their order) and the Order terms
        // come last.
        let sort = |v: &mut Vec<(Vec<(usize, i64)>, Term)>, reverse: bool| {
            let cmp = |a: &(Vec<(usize, i64)>, Term), b: &(Vec<(usize, i64)>, Term)| if reverse { compare(b, a) } else { compare(a, b) };
            if has_nan {
                // nan makes the order partial: a stable insertion pass keeps
                // the Python semantics without Rust's total-order contract.
                for i in 1..v.len() { let mut j = i; while j > 0 && cmp(&v[j], &v[j - 1]).is_lt() { v.swap(j, j - 1); j -= 1; } }
            } else {
                v.sort_by(cmp);
            }
        };
        if has_order {
            let (mut normal, mut orders): (Vec<_>, Vec<_>) = keyed.into_iter().partition(|(_, t)| !t.order);
            sort(&mut normal, true);
            sort(&mut orders, true);
            normal.extend(orders);
            keyed = normal;
        } else {
            sort(&mut keyed, false);
        }
        keyed.into_iter().map(|(_, t)| t.expr).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::parser::parse_str;

    fn s(x: &str) -> String {
        parse_str(x).unwrap().to_string()
    }

    #[test]
    fn dummy_coefficient_and_power_keys_wrap_the_atomic_key() {
        use super::*;
        let d=Expr::dummy_with_index("x",crate::Assumptions::NONE,37.into());
        let atom=sort_key(&d);
        for expression in [crate::ops::neg(&d),crate::ops::mul2(&Expr::int(2),&d),d.powi(2)] {
            let key=sort_key(&expression);
            let KArgs::List(children)=&key.args else {panic!("a Dummy expression key wraps its atomic key")};
            assert_eq!(children.len(),1);assert_eq!(*children[0],*atom);
        }
        assert!(matches!(atom.args,KArgs::Dummy(..)));
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
