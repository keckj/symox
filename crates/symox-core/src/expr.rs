//! Core expression tree.
//!
//! An [`Expr`] is a cheap-to-clone reference counted handle on an immutable
//! [`Node`]. Every node caches its structural hash and a 64-bit "symbol mask"
//! (a one-word Bloom filter over the symbols it contains) which makes
//! equality tests, hashing and `has`-style pruning O(1) in the common case.
//!
//! Nodes are always built in canonical form by the constructors in
//! [`crate::ops`]; the raw [`Expr::from_kind`] constructor is reserved for
//! code that already guarantees canonicity.

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, OnceLock};

use rustc_hash::FxHasher;

use crate::number::Number;

/// An undefined function head, distinct from builtin names and held operators.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct UndefinedFunction {
    pub name: Arc<str>,
    pub assumptions: crate::symbol::Assumptions,
    pub nargs: Option<Vec<usize>>,
    pub internal: bool,
}
impl UndefinedFunction {
    pub fn new(name:&str, assumptions:crate::symbol::Assumptions, mut nargs:Option<Vec<usize>>)->Self {
        if let Some(n)=&mut nargs {n.sort_unstable();n.dedup();}
        Self{name:Arc::from(name),assumptions:assumptions.closure().unwrap_or(assumptions),nargs,internal:false}
    }
    pub fn is_internal(&self,name:&str)->bool {self.internal && self.name.as_ref()==name}
    pub fn accepts(&self,n:usize)->bool {self.nargs.as_ref().is_none_or(|v|v.contains(&n))}
}
impl From<&str> for UndefinedFunction {
    fn from(name:&str)->Self {
        let mut function=Self::new(name,crate::symbol::Assumptions::NONE,None);
        function.internal=matches!(name,"Tuple"|"UnevaluatedExpr"|"Integral"|"Sum");function
    }
}
impl From<String> for UndefinedFunction {fn from(name:String)->Self {Self::from(name.as_str())}}
impl std::ops::Deref for UndefinedFunction {type Target=str;fn deref(&self)->&str{&self.name}}
impl AsRef<str> for UndefinedFunction {fn as_ref(&self)->&str{&self.name}}
impl std::fmt::Display for UndefinedFunction {fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{f.write_str(&self.name)}}

/// Named mathematical constants and special values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Constant {
    Pi,
    E,
    I,
    EulerGamma,
    /// `oo`
    Infinity,
    /// `-oo`
    NegInfinity,
    /// `zoo`
    ComplexInfinity,
    /// `nan`
    NaN,
}

impl Constant {
    pub fn name(self) -> &'static str {
        match self {
            Constant::Pi => "pi",
            Constant::E => "E",
            Constant::I => "I",
            Constant::EulerGamma => "EulerGamma",
            Constant::Infinity => "oo",
            Constant::NegInfinity => "-oo",
            Constant::ComplexInfinity => "zoo",
            Constant::NaN => "nan",
        }
    }

    pub fn is_infinite(self) -> bool {
        matches!(
            self,
            Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity
        )
    }
}

/// Known functions. `Undefined` is an applied undefined function `f(x, ...)`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Func {
    Exp,
    Log,
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Acot,
    Atan2,
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Asinh,
    Acosh,
    Atanh,
    Abs,
    Sign,
    Floor,
    Ceiling,
    Gamma,
    LogGamma,
    Erf,
    Erfc,
    Re,
    Im,
    Conjugate,
    Arg,
    Heaviside,
    DiracDelta,
    Max,
    Min,
    /// Special functions (factorial, zeta, Bessel, ...), see [`crate::special`].
    Special(crate::special::SpecialFn),
    /// Landau order term `O(expr, (x, x0), ...)`; args are
    /// `[expr, x1, p1, x2, p2, ...]`. Built only through
    /// [`crate::series::order::order`], absorbing rules live in `ops`.
    Order,
    Undefined(UndefinedFunction),
    /// Held substitution: [expression, Tuple(variables), Tuple(points)].
    Subs,
    /// Boolean connectives (`And`, `Or`, `Not`, ...), see [`crate::logic`].
    Logic(crate::logic::LogicOp),
    /// Sets (`Interval`, `FiniteSet`, `Union`, ...), see [`crate::sets`].
    Set(crate::sets::SetOp),
    /// A scalar entry of an ordered, typed matrix expression.
    MatrixElement(Arc<crate::matrixexpr::MatrixExpr>),
    /// Internal matrix atom used by scalar traversal/substitution APIs.
    MatrixProxy(Arc<crate::matrixexpr::MatrixExpr>),
    MatrixTrace(Arc<crate::matrixexpr::MatrixExpr>),
    MatrixDeterminant(Arc<crate::matrixexpr::MatrixExpr>),
}

impl Func {
    pub fn name(&self) -> &str {
        match self {
            Func::Exp => "exp",
            Func::Log => "log",
            Func::Sin => "sin",
            Func::Cos => "cos",
            Func::Tan => "tan",
            Func::Cot => "cot",
            Func::Sec => "sec",
            Func::Csc => "csc",
            Func::Asin => "asin",
            Func::Acos => "acos",
            Func::Atan => "atan",
            Func::Acot => "acot",
            Func::Atan2 => "atan2",
            Func::Sinh => "sinh",
            Func::Cosh => "cosh",
            Func::Tanh => "tanh",
            Func::Coth => "coth",
            Func::Asinh => "asinh",
            Func::Acosh => "acosh",
            Func::Atanh => "atanh",
            Func::Abs => "Abs",
            Func::Sign => "sign",
            Func::Floor => "floor",
            Func::Ceiling => "ceiling",
            Func::Gamma => "gamma",
            Func::LogGamma => "loggamma",
            Func::Erf => "erf",
            Func::Erfc => "erfc",
            Func::Re => "re",
            Func::Im => "im",
            Func::Conjugate => "conjugate",
            Func::Arg => "arg",
            Func::Heaviside => "Heaviside",
            Func::DiracDelta => "DiracDelta",
            Func::Max => "Max",
            Func::Min => "Min",
            Func::Special(s) => s.name(),
            Func::Order => "Order",
            Func::Undefined(n) => n,
            Func::Subs=>"Subs",
            Func::Logic(op) => op.name(),
            Func::Set(op) => op.name(),
            Func::MatrixElement(_) => "MatrixElement",
            Func::MatrixProxy(_) => "MatrixSymbol",
            Func::MatrixTrace(_) => "Trace",
            Func::MatrixDeterminant(_) => "Determinant",
        }
    }

    /// Look up a known function by its SymPy name.
    pub fn from_name(name: &str) -> Option<Func> {
        Some(match name {
            "exp" => Func::Exp,
            "log" | "ln" => Func::Log,
            "sin" => Func::Sin,
            "cos" => Func::Cos,
            "tan" => Func::Tan,
            "cot" => Func::Cot,
            "sec" => Func::Sec,
            "csc" => Func::Csc,
            "asin" => Func::Asin,
            "acos" => Func::Acos,
            "atan" => Func::Atan,
            "acot" => Func::Acot,
            "atan2" => Func::Atan2,
            "sinh" => Func::Sinh,
            "cosh" => Func::Cosh,
            "tanh" => Func::Tanh,
            "coth" => Func::Coth,
            "asinh" => Func::Asinh,
            "acosh" => Func::Acosh,
            "atanh" => Func::Atanh,
            "Abs" | "abs" => Func::Abs,
            "sign" => Func::Sign,
            "floor" => Func::Floor,
            "ceiling" => Func::Ceiling,
            "gamma" => Func::Gamma,
            "loggamma" => Func::LogGamma,
            "erf" => Func::Erf,
            "erfc" => Func::Erfc,
            "re" => Func::Re,
            "im" => Func::Im,
            "conjugate" => Func::Conjugate,
            "arg" => Func::Arg,
            "Heaviside" => Func::Heaviside,
            "DiracDelta" => Func::DiracDelta,
            "Max" => Func::Max,
            "Min" => Func::Min,
            "Order" => Func::Order,
            "Subs"=>Func::Subs,
            _ => {
                return crate::logic::LogicOp::from_name(name)
                    .map(Func::Logic)
                    .or_else(|| crate::sets::SetOp::from_name(name).map(Func::Set))
                    .or_else(|| crate::special::SpecialFn::from_name(name).map(Func::Special))
            }
        })
    }

    /// Number of arguments, `None` when variadic/unknown.
    pub fn arity(&self) -> Option<usize> {
        match self {
            Func::Atan2 => Some(2),
            Func::MatrixElement(_) => Some(2),
            Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) => Some(0),
            Func::Max | Func::Min | Func::Order | Func::Undefined(_) | Func::Subs | Func::Logic(_) | Func::Set(_) => None,
            Func::Special(s) => s.arity(),
            _ => Some(1),
        }
    }
}

/// Relational operators.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum RelOp {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

impl RelOp {
    pub fn symbol(self) -> &'static str {
        match self {
            RelOp::Eq => "==",
            RelOp::Ne => "!=",
            RelOp::Lt => "<",
            RelOp::Le => "<=",
            RelOp::Gt => ">",
            RelOp::Ge => ">=",
        }
    }
}

/// The payload of a node.
#[derive(Clone, Debug)]
pub enum Kind {
    Num(Number),
    Sym(crate::symbol::Symbol),
    Const(Constant),
    /// `coeff + Σ c_i * t_i`, terms sorted, `t_i` neither number nor carrying
    /// a numeric coefficient, all `c_i != 0`, at least one term.
    Add(Number, Vec<(Expr, Number)>),
    /// `coeff * Π b_i ** e_i`, bases sorted, all exponents non-zero.
    Mul(Number, Vec<(Expr, Expr)>),
    Pow(Expr, Expr),
    Func(Func, Vec<Expr>),
    /// Unevaluated derivative `d^n expr / dv_1^n_1 ...`.
    Derivative(Expr, Vec<(Expr, u32)>),
    /// Relational `lhs op rhs`.
    Rel(RelOp, Expr, Expr),
    /// `Piecewise((e1, c1), (e2, c2), ...)`.
    Piecewise(Vec<(Expr, Expr)>),
    /// Boolean literal (`true`/`false`) for conditions.
    Bool(bool),
}

impl Kind {
    fn rank(&self) -> u8 {
        match self {
            Kind::Num(_) => 0,
            Kind::Const(_) => 1,
            Kind::Sym(_) => 2,
            Kind::Pow(..) => 3,
            Kind::Mul(..) => 4,
            Kind::Add(..) => 5,
            Kind::Func(..) => 6,
            Kind::Derivative(..) => 7,
            Kind::Rel(..) => 8,
            Kind::Piecewise(..) => 9,
            Kind::Bool(_) => 10,
        }
    }
}

pub struct Node {
    hash: u64,
    mask: u64,
    kind: Kind,
    canonical: bool,
    has_matrix: bool,
    free_symbols: OnceLock<Vec<Expr>>,
    ordered_args: OnceLock<Vec<Expr>>,
}

/// A symbolic expression.
#[derive(Clone)]
pub struct Expr(Arc<Node>);

/// Extra bit set in the mask of `Order` nodes (and so of every expression
/// containing one): lets `add`/`mul` detect order terms in O(1).
pub const ORDER_MASK: u64 = 1 << 63;

#[inline]
fn sym_bit(name: &str) -> u64 {
    let mut h = FxHasher::default();
    name.hash(&mut h);
    1u64 << (h.finish() % 64)
}

impl Expr {
    /// Build a node from an already canonical payload.
    pub fn from_kind(kind: Kind) -> Expr {
        Self::from_payload(kind, true)
    }

    /// Preserve the argument order and grouping requested with evaluate=False.
    /// Algorithms may canonicalize these nodes when their operands change.
    pub fn from_unevaluated_kind(kind: Kind) -> Expr {
        Self::from_payload(kind, false)
    }

    fn from_payload(kind: Kind, canonical: bool) -> Expr {
        let mut h = FxHasher::default();
        let mut mask = 0u64;
        let mut has_matrix = false;
        kind.rank().hash(&mut h);
        match &kind {
            Kind::Num(n) => n.hash(&mut h),
            Kind::Sym(s) => {
                s.hash(&mut h);
                mask = sym_bit(s.name());
            }
            Kind::Const(c) => c.hash(&mut h),
            Kind::Add(c, terms) => {
                c.hash(&mut h);
                for (t, k) in terms {
                    h.write_u64(t.hash64());
                    k.hash(&mut h);
                    mask |= t.mask();
                    has_matrix |= t.has_matrix();
                }
            }
            Kind::Mul(c, facs) => {
                c.hash(&mut h);
                for (b, e) in facs {
                    h.write_u64(b.hash64());
                    h.write_u64(e.hash64());
                    mask |= b.mask() | e.mask();
                    has_matrix |= b.has_matrix() || e.has_matrix();
                }
            }
            Kind::Pow(b, e) => {
                h.write_u64(b.hash64());
                h.write_u64(e.hash64());
                mask = b.mask() | e.mask();
                has_matrix = b.has_matrix() || e.has_matrix();
            }
            Kind::Func(f, args) => {
                f.hash(&mut h);
                for a in args {
                    h.write_u64(a.hash64());
                    mask |= a.mask();
                    has_matrix |= a.has_matrix();
                }
                if *f == Func::Order {
                    mask |= ORDER_MASK;
                }
                if let Func::MatrixElement(parent) | Func::MatrixProxy(parent) | Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent) = f {
                    has_matrix = true;
                    for matrix in parent.matrix_symbols() {
                        // Include the whole typed atom, rather than dimensions,
                        // in the dependency mask just as free_symbols does.
                        let mut atom_hash = FxHasher::default();
                        matrix.hash(&mut atom_hash);
                        mask |= 1u64 << (atom_hash.finish() % 64);
                    }
                    for scalar in parent.scalar_free_symbols() {
                        mask |= scalar.mask();
                    }
                    let _ = parent.scalar_map(&mut |scalar| {
                        mask |= scalar.mask();
                        scalar.clone()
                    });
                }
            }
            Kind::Derivative(e, vars) => {
                h.write_u64(e.hash64());
                mask = e.mask();
                has_matrix = e.has_matrix();
                for (v, n) in vars {
                    h.write_u64(v.hash64());
                    n.hash(&mut h);
                    mask |= v.mask();
                    has_matrix |= v.has_matrix();
                }
            }
            Kind::Rel(op, a, b) => {
                op.hash(&mut h);
                h.write_u64(a.hash64());
                h.write_u64(b.hash64());
                mask = a.mask() | b.mask();
                has_matrix = a.has_matrix() || b.has_matrix();
            }
            Kind::Piecewise(pieces) => {
                for (e, c) in pieces {
                    h.write_u64(e.hash64());
                    h.write_u64(c.hash64());
                    mask |= e.mask() | c.mask();
                    has_matrix |= e.has_matrix() || c.has_matrix();
                }
            }
            Kind::Bool(b) => b.hash(&mut h),
        }
        Expr(Arc::new(Node {
            hash: h.finish(),
            mask,
            kind,
            canonical,
            has_matrix,
            free_symbols: OnceLock::new(),
            ordered_args: OnceLock::new(),
        }))
    }

    #[inline]
    pub fn kind(&self) -> &Kind {
        &self.0.kind
    }

    /// Whether this scalar expression contains a typed matrix entry.
    #[inline]
    pub fn has_matrix(&self) -> bool { self.0.has_matrix }

    /// Whether this node was made by a canonical constructor.
    #[inline]
    pub fn is_canonical(&self) -> bool {
        self.0.canonical
    }

    #[inline]
    pub fn hash64(&self) -> u64 {
        self.0.hash
    }

    /// Bloom mask of contained symbols.
    #[inline]
    pub fn mask(&self) -> u64 {
        self.0.mask
    }

    /// Cache symbol discovery on immutable non-atomic nodes. Symbols themselves
    /// must never cache a clone of themselves: that would create an Arc cycle.
    pub(crate) fn cached_free_symbols(&self, init: impl FnOnce() -> Vec<Expr>) -> Vec<Expr> {
        if self.is_atom() {
            return init();
        }
        self.0.free_symbols.get_or_init(init).clone()
    }

    /// The Python facade's canonical argument order differs from the storage
    /// order. Cache it lazily on the DAG node, shared by all Python wrappers.
    pub fn cached_ordered_args(&self, init: impl FnOnce() -> Vec<Expr>) -> Vec<Expr> {
        if self.is_atom() {
            return Vec::new();
        }
        self.0.ordered_args.get_or_init(init).clone()
    }

    #[inline]
    pub fn ptr_eq(&self, o: &Expr) -> bool {
        Arc::ptr_eq(&self.0, &o.0)
    }

    // ---------------------------------------------------------------
    // atoms
    // ---------------------------------------------------------------

    pub fn symbol(name: &str) -> Expr {
        Expr::from_kind(Kind::Sym(crate::symbol::Symbol::new(name)))
    }

    /// A symbol with assumptions (`Symbol('x', positive=True)`).
    pub fn symbol_with(name: &str, a: crate::symbol::Assumptions) -> Expr {
        Expr::from_kind(Kind::Sym(crate::symbol::Symbol::with_assumptions(name, a)))
    }

    /// A pattern wildcard (`Wild('a', exclude=[...])`).
    pub fn wild(name: &str, exclude: Vec<Expr>) -> Expr {
        Expr::from_kind(Kind::Sym(crate::symbol::Symbol::wild(name, exclude)))
    }

    /// The symbol payload, if this is a symbol.
    pub fn as_sym(&self) -> Option<&crate::symbol::Symbol> {
        match self.kind() {
            Kind::Sym(s) => Some(s),
            _ => None,
        }
    }

    pub fn num(n: Number) -> Expr {
        Expr::from_kind(Kind::Num(n))
    }

    pub fn int(i: i64) -> Expr {
        Expr::num(Number::Int(i))
    }

    pub fn rational(n: i64, d: i64) -> Expr {
        match Number::rational_i64(n, d) {
            Some(r) => Expr::num(r),
            None => Expr::constant(if n == 0 {
                Constant::NaN
            } else {
                Constant::ComplexInfinity
            }),
        }
    }

    pub fn float(f: f64) -> Expr {
        if f.is_nan() {
            return Expr::constant(Constant::NaN);
        }
        if f.is_infinite() {
            return Expr::constant(if f > 0.0 {
                Constant::Infinity
            } else {
                Constant::NegInfinity
            });
        }
        Expr::num(Number::float(f))
    }

    pub fn constant(c: Constant) -> Expr {
        Expr::from_kind(Kind::Const(c))
    }

    pub fn boolean(b: bool) -> Expr {
        Expr::from_kind(Kind::Bool(b))
    }

    pub fn zero() -> Expr {
        Expr::int(0)
    }
    pub fn one() -> Expr {
        Expr::int(1)
    }
    pub fn minus_one() -> Expr {
        Expr::int(-1)
    }
    pub fn pi() -> Expr {
        Expr::constant(Constant::Pi)
    }
    pub fn e() -> Expr {
        Expr::constant(Constant::E)
    }
    pub fn i() -> Expr {
        Expr::constant(Constant::I)
    }
    pub fn oo() -> Expr {
        Expr::constant(Constant::Infinity)
    }
    pub fn neg_oo() -> Expr {
        Expr::constant(Constant::NegInfinity)
    }
    pub fn zoo() -> Expr {
        Expr::constant(Constant::ComplexInfinity)
    }
    pub fn nan() -> Expr {
        Expr::constant(Constant::NaN)
    }

    // ---------------------------------------------------------------
    // predicates / accessors
    // ---------------------------------------------------------------

    #[inline]
    pub fn as_number(&self) -> Option<&Number> {
        match self.kind() {
            Kind::Num(n) => Some(n),
            _ => None,
        }
    }

    #[inline]
    pub fn is_number(&self) -> bool {
        matches!(self.kind(), Kind::Num(_))
    }

    #[inline]
    pub fn is_symbol(&self) -> bool {
        matches!(self.kind(), Kind::Sym(_))
    }

    pub fn as_symbol(&self) -> Option<&str> {
        match self.kind() {
            Kind::Sym(s) => Some(s.name()),
            _ => None,
        }
    }

    pub fn is_zero(&self) -> bool {
        matches!(self.kind(), Kind::Num(n) if n.is_zero())
    }

    pub fn is_one(&self) -> bool {
        matches!(self.kind(), Kind::Num(n) if n.is_one())
    }

    pub fn is_minus_one(&self) -> bool {
        matches!(self.kind(), Kind::Num(n) if n.is_minus_one())
    }

    pub fn is_integer_number(&self) -> bool {
        matches!(self.kind(), Kind::Num(n) if n.is_integer())
    }

    pub fn is_rational_number(&self) -> bool {
        matches!(self.kind(), Kind::Num(n) if n.is_rational())
    }

    pub fn as_i64(&self) -> Option<i64> {
        self.as_number().and_then(|n| n.to_i64())
    }

    pub fn is_constant(&self, c: Constant) -> bool {
        matches!(self.kind(), Kind::Const(k) if *k == c)
    }

    pub fn is_add(&self) -> bool {
        matches!(self.kind(), Kind::Add(..))
    }
    pub fn is_mul(&self) -> bool {
        matches!(self.kind(), Kind::Mul(..))
    }
    pub fn is_pow(&self) -> bool {
        matches!(self.kind(), Kind::Pow(..))
    }

    pub fn is_func(&self, f: &Func) -> bool {
        matches!(self.kind(), Kind::Func(g, _) if g == f)
    }

    /// Is this an atom (no sub-expressions)?
    pub fn is_atom(&self) -> bool {
        matches!(
            self.kind(),
            Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_)
        )
    }

    /// Fast negative test: `false` means the expression certainly contains
    /// no `Order` term.
    #[inline]
    pub fn may_have_order(&self) -> bool {
        self.mask() & ORDER_MASK != 0
    }

    /// Fast negative test: `false` means the expression certainly does not
    /// contain symbol `s`.
    #[inline]
    pub fn may_contain(&self, s: &Expr) -> bool {
        self.mask() & s.mask() != 0
    }

    /// SymPy-style class name (`Add`, `Mul`, `Symbol`, `sin`...).
    pub fn class_name(&self) -> &str {
        match self.kind() {
            Kind::Num(Number::Int(_)) | Kind::Num(Number::BigInt(_)) => {
                match self.as_i64() {
                    Some(0) => "Zero",
                    Some(1) => "One",
                    Some(-1) => "NegativeOne",
                    _ => "Integer",
                }
            }
            Kind::Num(Number::Rational(r)) => {
                if r.numer() == &num_bigint::BigInt::from(1) && r.denom() == &num_bigint::BigInt::from(2) {
                    "Half"
                } else {
                    "Rational"
                }
            }
            Kind::Num(Number::Float(_) | Number::MpFloat(..)) => "Float",
            // SymPy's `Dummy` symbols print as `_name`; symox represents them
            // as symbols whose name starts with an underscore.
            Kind::Sym(s) if s.name().starts_with('_') && !s.is_wild() => "Dummy",
            Kind::Sym(_) => "Symbol",
            Kind::Const(c) => match c {
                Constant::Pi => "Pi",
                Constant::E => "Exp1",
                Constant::I => "ImaginaryUnit",
                Constant::EulerGamma => "EulerGamma",
                Constant::Infinity => "Infinity",
                Constant::NegInfinity => "NegativeInfinity",
                Constant::ComplexInfinity => "ComplexInfinity",
                Constant::NaN => "NaN",
            },
            Kind::Add(..) => "Add",
            Kind::Mul(..) => "Mul",
            Kind::Pow(..) => "Pow",
            Kind::Func(f, _) => f.name(),
            Kind::Derivative(..) => "Derivative",
            Kind::Rel(op, ..) => match op {
                RelOp::Eq => "Equality",
                RelOp::Ne => "Unequality",
                RelOp::Lt => "StrictLessThan",
                RelOp::Le => "LessThan",
                RelOp::Gt => "StrictGreaterThan",
                RelOp::Ge => "GreaterThan",
            },
            Kind::Piecewise(_) => "Piecewise",
            Kind::Bool(true) => "BooleanTrue",
            Kind::Bool(false) => "BooleanFalse",
        }
    }

    /// Structural total order (independent of hashes). Used for printing
    /// and as a tie-breaker.
    pub fn cmp_struct(&self, o: &Expr) -> Ordering {
        if self.ptr_eq(o) {
            return Ordering::Equal;
        }
        let (a, b) = (self.kind(), o.kind());
        let r = a.rank().cmp(&b.rank());
        if r != Ordering::Equal {
            return r;
        }
        match (a, b) {
            (Kind::Num(x), Kind::Num(y)) => x.cmp_struct(y),
            (Kind::Sym(x), Kind::Sym(y)) => x.cmp(y),
            (Kind::Const(x), Kind::Const(y)) => x.cmp(y),
            (Kind::Bool(x), Kind::Bool(y)) => x.cmp(y),
            (Kind::Add(c1, t1), Kind::Add(c2, t2)) => t1
                .len()
                .cmp(&t2.len())
                .then_with(|| {
                    for ((a, ka), (b, kb)) in t1.iter().zip(t2) {
                        let r = a.cmp_struct(b).then_with(|| ka.cmp_struct(kb));
                        if r != Ordering::Equal {
                            return r;
                        }
                    }
                    Ordering::Equal
                })
                .then_with(|| c1.cmp_struct(c2)),
            (Kind::Mul(c1, f1), Kind::Mul(c2, f2)) => f1
                .len()
                .cmp(&f2.len())
                .then_with(|| {
                    for ((a, ea), (b, eb)) in f1.iter().zip(f2) {
                        let r = a.cmp_struct(b).then_with(|| ea.cmp_struct(eb));
                        if r != Ordering::Equal {
                            return r;
                        }
                    }
                    Ordering::Equal
                })
                .then_with(|| c1.cmp_struct(c2)),
            (Kind::Pow(b1, e1), Kind::Pow(b2, e2)) => {
                b1.cmp_struct(b2).then_with(|| e1.cmp_struct(e2))
            }
            (Kind::Func(f1, a1), Kind::Func(f2, a2)) => {
                f1.cmp(f2).then_with(|| cmp_slices(a1, a2))
            }
            (Kind::Derivative(e1, v1), Kind::Derivative(e2, v2)) => {
                e1.cmp_struct(e2).then_with(|| {
                    v1.len().cmp(&v2.len()).then_with(|| {
                        for ((a, na), (b, nb)) in v1.iter().zip(v2) {
                            let r = a.cmp_struct(b).then_with(|| na.cmp(nb));
                            if r != Ordering::Equal {
                                return r;
                            }
                        }
                        Ordering::Equal
                    })
                })
            }
            (Kind::Rel(o1, a1, b1), Kind::Rel(o2, a2, b2)) => o1
                .cmp(o2)
                .then_with(|| a1.cmp_struct(a2))
                .then_with(|| b1.cmp_struct(b2)),
            (Kind::Piecewise(p1), Kind::Piecewise(p2)) => {
                p1.len().cmp(&p2.len()).then_with(|| {
                    for ((a, ca), (b, cb)) in p1.iter().zip(p2) {
                        let r = a.cmp_struct(b).then_with(|| ca.cmp_struct(cb));
                        if r != Ordering::Equal {
                            return r;
                        }
                    }
                    Ordering::Equal
                })
            }
            _ => unreachable!(),
        }
    }

    /// Direct children, in a SymPy-compatible way where cheap.
    pub fn args(&self) -> Vec<Expr> {
        use crate::ops::{add_term, mul_factor};
        match self.kind() {
            Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => vec![],
            Kind::Add(c, terms) => {
                let mut v: Vec<Expr> = terms.iter().map(|(t, k)| add_term(t, k)).collect();
                if !c.is_zero() {
                    v.push(Expr::num(c.clone()));
                }
                v
            }
            Kind::Mul(c, facs) => {
                let mut v = Vec::with_capacity(facs.len() + 1);
                if !c.is_one() {
                    v.push(Expr::num(c.clone()));
                }
                v.extend(facs.iter().map(|(b, e)| mul_factor(b, e)));
                v
            }
            Kind::Pow(b, e) => vec![b.clone(), e.clone()],
            Kind::Func(_, a) => a.clone(),
            Kind::Derivative(e, vars) => {
                let mut v = vec![e.clone()];
                for (x, n) in vars {
                    for _ in 0..*n {
                        v.push(x.clone());
                    }
                }
                v
            }
            Kind::Rel(_, a, b) => vec![a.clone(), b.clone()],
            Kind::Piecewise(p) => p
                .iter()
                .flat_map(|(e, c)| [e.clone(), c.clone()])
                .collect(),
        }
    }
}

fn cmp_slices(a: &[Expr], b: &[Expr]) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| {
        for (x, y) in a.iter().zip(b) {
            let r = x.cmp_struct(y);
            if r != Ordering::Equal {
                return r;
            }
        }
        Ordering::Equal
    })
}

fn kind_eq(a: &Kind, b: &Kind) -> bool {
    match (a, b) {
        (Kind::Num(x), Kind::Num(y)) => x == y,
        (Kind::Sym(x), Kind::Sym(y)) => x == y,
        (Kind::Const(x), Kind::Const(y)) => x == y,
        (Kind::Bool(x), Kind::Bool(y)) => x == y,
        (Kind::Add(c1, t1), Kind::Add(c2, t2)) => c1 == c2 && t1 == t2,
        (Kind::Mul(c1, f1), Kind::Mul(c2, f2)) => c1 == c2 && f1 == f2,
        (Kind::Pow(b1, e1), Kind::Pow(b2, e2)) => b1 == b2 && e1 == e2,
        (Kind::Func(f1, a1), Kind::Func(f2, a2)) => f1 == f2 && a1 == a2,
        (Kind::Derivative(e1, v1), Kind::Derivative(e2, v2)) => e1 == e2 && v1 == v2,
        (Kind::Rel(o1, a1, b1), Kind::Rel(o2, a2, b2)) => o1 == o2 && a1 == a2 && b1 == b2,
        (Kind::Piecewise(p1), Kind::Piecewise(p2)) => p1 == p2,
        _ => false,
    }
}

impl PartialEq for Expr {
    #[inline]
    fn eq(&self, o: &Expr) -> bool {
        self.ptr_eq(o) || (self.0.hash == o.0.hash && kind_eq(self.kind(), o.kind()))
    }
}
impl Eq for Expr {}

impl Hash for Expr {
    #[inline]
    fn hash<H: Hasher>(&self, h: &mut H) {
        h.write_u64(self.0.hash)
    }
}

/// Canonical (fast) ordering: by hash, then structurally.
impl Ord for Expr {
    #[inline]
    fn cmp(&self, o: &Expr) -> Ordering {
        if self.ptr_eq(o) {
            return Ordering::Equal;
        }
        self.0.hash.cmp(&o.0.hash).then_with(|| self.cmp_struct(o))
    }
}
impl PartialOrd for Expr {
    fn partial_cmp(&self, o: &Expr) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

impl fmt::Debug for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", crate::printer::srepr(self))
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", crate::printer::to_str(self))
    }
}
