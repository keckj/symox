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
use std::sync::atomic::{AtomicU64,Ordering as AtomicOrdering};

use rustc_hash::{FxHasher, FxHashMap};
use std::cell::RefCell;

use crate::number::Number;
#[path="binding_identity.rs"]
mod binding_identity;

// A comparison may cross scalar/matrix metadata several times. Keep its
// physical-pair memo for the whole call, rather than resetting it at each
// matrix/scalar boundary. RefCell borrows are always released before recursion.
struct ComparisonMemo<T> {
    depth: usize,
    values: Option<FxHashMap<(usize, usize), T>>,
}
impl<T> ComparisonMemo<T> {
    const fn new() -> Self { Self { depth: 0, values: None } }
    fn leave(&mut self) {
        self.depth -= 1;
        if self.depth == 0 { self.values = None; }
    }
}
thread_local! {
    static EQ_MEMO: RefCell<ComparisonMemo<bool>> = const { RefCell::new(ComparisonMemo::new()) };
    static ORDER_MEMO: RefCell<ComparisonMemo<Ordering>> = const { RefCell::new(ComparisonMemo::new()) };
}
struct ComparisonGuard { equality: bool }
impl ComparisonGuard {
    fn equality() -> Self {
        EQ_MEMO.with(|memo| memo.borrow_mut().depth += 1);
        Self { equality: true }
    }
    fn order() -> Self {
        ORDER_MEMO.with(|memo| memo.borrow_mut().depth += 1);
        Self { equality: false }
    }
}
impl Drop for ComparisonGuard {
    fn drop(&mut self) {
        if self.equality { EQ_MEMO.with(|memo| memo.borrow_mut().leave()); }
        else { ORDER_MEMO.with(|memo| memo.borrow_mut().leave()); }
    }
}

/// Immutable facts shared by all applications of one undefined function head.
/// Head identity preserves whether commutative=True was explicitly declared.
#[derive(Debug)]
struct FunctionFacts {
    head: crate::assumption_facts::FactSet,
    applied: crate::assumption_facts::FactSet,
    hash_extra: bool,
}

/// An undefined function head, distinct from builtin names and held operators.
#[derive(Clone, Debug)]
pub struct UndefinedFunction {
    pub name: Arc<str>,
    /// Compatibility projection; full declarations are available through facts().
    pub assumptions: crate::symbol::Assumptions,
    pub nargs: Option<Vec<usize>>,
    pub internal: bool,
    facts: Arc<FunctionFacts>,
    /// Host-language class defining this head, part of its identity.
    host: Option<Arc<dyn crate::host::HostClass>>,
}
impl UndefinedFunction {
    pub fn new(name:&str, assumptions:crate::symbol::Assumptions,nargs:Option<Vec<usize>>)->Self {
        Self::with_facts(name,crate::assumption_facts::FactSet::from_legacy(assumptions),nargs)
            .expect("consistent legacy function assumptions and arity")
    }
    pub fn with_facts(name:&str,given:crate::assumption_facts::FactSet,mut nargs:Option<Vec<usize>>)->Result<Self,String>{
        use crate::assumption_facts::FactSet;
        if let Some(counts)=&mut nargs {
            if counts.is_empty(){return Err("nargs must not be empty".into());}
            counts.sort_unstable();counts.dedup();
        }
        static DEFAULT:OnceLock<Arc<FunctionFacts>>=OnceLock::new();
        if given.is_empty(){
            let facts=DEFAULT.get_or_init(||Arc::new(FunctionFacts{head:FactSet::EMPTY,applied:FactSet::EMPTY,hash_extra:false})).clone();
            return Ok(Self{name:Arc::from(name),assumptions:crate::Assumptions::NONE,nargs,internal:false,facts,host:None});
        }
        // SymPy closes through Symbol, then removes implicit commutativity from
        // the head's _kwargs. An explicit True remains part of head identity.
        let closed=given.closure()?;
        let head=if given.get("commutative").is_none(){closed.without("commutative")}else{closed};
        let applied=head.closure_without_default_commutative()?;
        let assumptions=applied.legacy();
        let old=FactSet::from_legacy(assumptions).closure()?.without("commutative");
        let hash_extra=head!=old;
        let facts=Arc::new(FunctionFacts{head,applied,hash_extra});
        Ok(Self{name:Arc::from(name),assumptions,nargs,internal:false,facts,host:None})
    }
    pub fn head_facts(&self)->crate::assumption_facts::FactSet{self.facts.head}
    pub fn facts(&self)->crate::assumption_facts::FactSet{self.facts.applied}
    pub fn is_internal(&self,name:&str)->bool {self.internal && self.name.as_ref()==name}
    pub fn accepts(&self,n:usize)->bool {self.nargs.as_ref().is_none_or(|v|v.contains(&n))}
    /// The same head defined by a host-language class.
    pub fn with_host(mut self,host:Arc<dyn crate::host::HostClass>)->Self{self.host=Some(host);self}
    pub fn host(&self)->Option<&Arc<dyn crate::host::HostClass>>{self.host.as_ref()}
    fn host_key(&self)->Option<(&str,u64)>{self.host.as_ref().map(|h|(h.class_name(),h.identity()))}
}
impl PartialEq for UndefinedFunction {
    fn eq(&self,other:&Self)->bool{
        self.name==other.name && self.assumptions==other.assumptions && self.nargs==other.nargs && self.internal==other.internal
            && (Arc::ptr_eq(&self.facts,&other.facts)||self.head_facts()==other.head_facts())
            && crate::host::same(self.host.as_ref(),other.host.as_ref())
    }
}
impl Eq for UndefinedFunction {}
impl Hash for UndefinedFunction {
    fn hash<H:Hasher>(&self,state:&mut H){
        // Preserve the old head hash for facts representable by the legacy flags.
        self.name.hash(state);self.assumptions.hash(state);self.nargs.hash(state);self.internal.hash(state);
        if self.facts.hash_extra{self.head_facts().hash(state);}
        if let Some(host)=&self.host{0x4057u16.hash(state);host.identity().hash(state);}
    }
}
impl Ord for UndefinedFunction {
    fn cmp(&self,other:&Self)->Ordering{
        self.name.cmp(&other.name).then_with(||self.assumptions.cmp(&other.assumptions))
            .then_with(||self.head_facts().cmp(&other.head_facts()))
            .then_with(||self.nargs.cmp(&other.nargs)).then_with(||self.internal.cmp(&other.internal))
            .then_with(||self.host_key().cmp(&other.host_key()))
            .then_with(||crate::host::order(self.host.as_ref(),other.host.as_ref()))
    }
}
impl PartialOrd for UndefinedFunction {fn partial_cmp(&self,other:&Self)->Option<Ordering>{Some(self.cmp(other))}}
impl From<&str> for UndefinedFunction {
    fn from(name:&str)->Self {
        let mut function=Self::new(name,crate::symbol::Assumptions::NONE,None);
        function.internal=matches!(name,"Tuple"|"UnevaluatedExpr"|"Integral"|"Sum"|"Product");function
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
    /// Catalan's constant `G`.
    Catalan,
    /// `(1 + sqrt(5))/2`
    GoldenRatio,
    /// Real root of `x**3 - x**2 - x - 1`.
    TribonacciConstant,
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
            Constant::Catalan => "Catalan",
            Constant::GoldenRatio => "GoldenRatio",
            Constant::TribonacciConstant => "TribonacciConstant",
        }
    }

    /// Double precision value of the real named constants.
    pub fn approx_f64(self) -> Option<f64> {
        Some(match self {
            Constant::Pi => std::f64::consts::PI,
            Constant::E => std::f64::consts::E,
            Constant::EulerGamma => 0.577_215_664_901_532_9,
            Constant::Catalan => 0.915_965_594_177_219,
            Constant::GoldenRatio => 1.618_033_988_749_895,
            Constant::TribonacciConstant => 1.839_286_755_214_161_1,
            Constant::Infinity => f64::INFINITY,
            Constant::NegInfinity => f64::NEG_INFINITY,
            _ => return None,
        })
    }

    /// SymPy `NumberSymbol`s: real positive named constants.
    pub fn is_number_symbol(self) -> bool {
        matches!(
            self,
            Constant::Pi
                | Constant::E
                | Constant::EulerGamma
                | Constant::Catalan
                | Constant::GoldenRatio
                | Constant::TribonacciConstant
        )
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
    /// Held limit: [expression, variable, point, Symbol(direction)].
    Limit,
    /// Held calculus binders: [body, Tuple(variable, bounds...), ...].
    Integral,
    Sum,
    Product,
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
    MatrixDotProduct(Arc<crate::matrixexpr::MatrixExpr>, Arc<crate::matrixexpr::MatrixExpr>),
    /// A predicate of the new assumptions (`Q.positive`), see [`crate::ask`].
    Predicate(crate::ask::Pred),
    /// `AppliedPredicate`/`AppliedBinaryRelation`: `[predicate, arguments...]`.
    AppliedPredicate,
    /// Polar exponential retains winding information instead of identifying
    /// arguments modulo 2*pi*I. Append new variants to preserve old node hashes.
    ExpPolar,
    /// Generic Expr constructor used by Basic.fromiter, with symbolic leaves.
    GenericExpr,
}

impl Func {
    pub fn name(&self) -> &str {
        match self {
            Func::Exp => "exp",
            Func::ExpPolar => "exp_polar",
            Func::GenericExpr => "Expr",
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
            Func::Limit=>"Limit",
            Func::Integral=>"Integral",
            Func::Sum=>"Sum",
            Func::Product=>"Product",
            Func::Logic(op) => op.name(),
            Func::Set(op) => op.name(),
            Func::MatrixElement(_) => "MatrixElement",
            Func::MatrixProxy(_) => "MatrixSymbol",
            Func::MatrixTrace(_) => "Trace",
            Func::MatrixDeterminant(_) => "Determinant",
            Func::MatrixDotProduct(_, _) => "DotProduct",
            Func::Predicate(p) => p.class_name(),
            Func::AppliedPredicate => "AppliedPredicate",
        }
    }

    /// Look up a known function by its SymPy name.
    pub fn from_name(name: &str) -> Option<Func> {
        Some(match name {
            "exp" => Func::Exp,
            "exp_polar" => Func::ExpPolar,
            "Expr" => Func::GenericExpr,
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
            "Limit"=>Func::Limit,
            "Integral"=>Func::Integral,
            "Sum"=>Func::Sum,
            "AppliedPredicate" | "AppliedBinaryRelation" => Func::AppliedPredicate,
            "Product"=>Func::Product,
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
            Func::GenericExpr => None,
            Func::Atan2 => Some(2),
            Func::Limit => Some(4),
            Func::MatrixElement(_) => Some(2),
            Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) | Func::Predicate(_) => Some(0),
            Func::AppliedPredicate => None,
            Func::Max | Func::Min | Func::Order | Func::Undefined(_) | Func::Subs | Func::Integral | Func::Sum | Func::Product | Func::Logic(_) | Func::Set(_) => None,
            Func::Special(s) => s.arity(),
            Func::DiracDelta => None,
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
    features: u16,
    binding_identity: Option<Arc<binding_identity::Key>>,
    free_symbols: OnceLock<Vec<Expr>>,
    ordered_args: OnceLock<Vec<Expr>>,
    float_info: OnceLock<FloatInfo>,
    // Three-valued facts have independent exact/numeric banks. This metadata
    // never contributes to structural identity, symbol masks or serialization.
    assumption_cache: AtomicU64,
}

/// Floating-point metadata of an immutable mixed scalar/matrix expression DAG.
#[derive(Clone,Copy,Debug,Default)]
pub struct FloatInfo {
    pub max_precision: Option<u32>,
    pub has_mp: bool,
}
impl FloatInfo {
    fn include(&mut self,other:Self){
        self.max_precision=self.max_precision.max(other.max_precision);
        self.has_mp |= other.has_mp;
    }
    fn number(n:&Number)->Self {Self{max_precision:n.float_precision(),has_mp:matches!(n,Number::MpFloat(..))}}
}

/// A symbolic expression.
#[derive(Clone)]
pub struct Expr(Arc<Node>);

/// A non-owning reference to a node, which identifies it while it is alive.
#[derive(Clone, Debug)]
pub struct WeakExpr(std::sync::Weak<Node>);
impl WeakExpr {
    /// Whether this reference denotes the live node `e`.
    pub fn is(&self, e: &Expr) -> bool { std::ptr::eq(self.0.as_ptr(), Arc::as_ptr(&e.0)) && self.0.strong_count() > 0 }
    pub fn upgrade(&self) -> Option<Expr> { self.0.upgrade().map(Expr) }
}

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
    pub fn from_unevaluated_kind(mut kind: Kind) -> Expr {
        // The first numerical factor has the same public position whether it
        // lives in the compact coefficient slot or in the held factor vector.
        // Normalize that storage choice, while retaining explicit unit factors.
        if let Kind::Mul(Number::Int(1),factors)=&mut kind {
            if let Some((first,power))=factors.first(){
                if power.is_one(){if let Some(number)=first.as_number().filter(|n|!n.is_one()){
                    let coefficient=number.clone();factors.remove(0);
                    kind=Kind::Mul(coefficient,std::mem::take(factors));
                }}
            }
        }
        Self::from_payload(kind, false)
    }

    fn from_payload(kind: Kind, canonical: bool) -> Expr {
        let mut h = FxHasher::default();
        let mut mask = 0u64;
        let mut has_matrix = false;
        let binding_identity=match &kind {Kind::Func(Func::Subs,args)=>binding_identity::build(args),_=>None};
        let mut features:u16=if binding_identity.is_some(){1}else{0};
        if matches!(&kind,Kind::Sym(symbol) if symbol.host().is_some())
            || matches!(&kind,Kind::Func(Func::Undefined(head),_) if head.host().is_some()){features|=256;}
        if matches!(&kind,Kind::Sym(symbol) if symbol.host().is_some_and(|host|host.substitutes()))
            || matches!(&kind,Kind::Func(Func::Undefined(head),_) if head.host().is_some_and(|host|host.substitutes())){features|=512;}
        if matches!(&kind,Kind::Sym(symbol) if symbol.host().is_some_and(|host|host.may_not_commute()))
            || matches!(&kind,Kind::Func(Func::Undefined(head),_) if head.host().is_some_and(|host|host.may_not_commute())){features|=1024;}
        if !canonical && matches!(kind,Kind::Add(..)|Kind::Mul(..)|Kind::Pow(..)){features|=16;}
        if matches!(&kind,Kind::Const(Constant::Infinity|Constant::NegInfinity|Constant::ComplexInfinity|Constant::NaN))
            || matches!(&kind,Kind::Num(Number::Float(value)) if !value.is_finite())
            || matches!(&kind,Kind::Sym(symbol) if symbol.facts().get("finite")==Some(false))
            || matches!(&kind,Kind::Func(Func::Undefined(head),_) if !head.internal && head.facts().get("finite")==Some(false)){features|=32;}
        if match &kind{
            Kind::Num(n)|Kind::Mul(n,_)=>n.has_wide_exponent(),
            Kind::Add(n,terms)=>n.has_wide_exponent()||terms.iter().any(|(_,n)|n.has_wide_exponent()),
            _=>false,
        }{features|=64;}
        kind.rank().hash(&mut h);
        match &kind {
            Kind::Num(n) => n.hash(&mut h),
            Kind::Sym(s) => {
                if s.wild_properties().is_some()||s.exclude().iter().any(Expr::has_wild_properties){features|=128;}
                s.hash(&mut h);
                mask = if s.is_dummy() { 1u64 << (h.finish() % 64) } else { sym_bit(s.name()) };
            }
            Kind::Const(c) => c.hash(&mut h),
            Kind::Add(c, terms) => {
                c.hash(&mut h);
                for (t, k) in terms {
                    h.write_u64(t.hash64());
                    k.hash(&mut h);
                    mask |= t.mask();
                    features |= t.0.features;
                    has_matrix |= t.has_matrix();
                }
            }
            Kind::Mul(c, facs) => {
                c.hash(&mut h);
                for (b, e) in facs {
                    h.write_u64(b.hash64());
                    h.write_u64(e.hash64());
                    mask |= b.mask() | e.mask();
                    features |= b.0.features | e.0.features;
                    has_matrix |= b.has_matrix() || e.has_matrix();
                }
            }
            Kind::Pow(b, e) => {
                h.write_u64(b.hash64());
                h.write_u64(e.hash64());
                mask = b.mask() | e.mask();
                features |= b.0.features | e.0.features;
                has_matrix = b.has_matrix() || e.has_matrix();
            }
            Kind::Func(f, args) => {
                f.hash(&mut h);
                if let Some(key)=&binding_identity {key.hash(&mut h);}
                if *f==Func::Limit {features|=2;}
                for a in args {
                    if binding_identity.is_none(){h.write_u64(a.hash64());}
                    features |= a.0.features;
                    mask |= a.mask();
                    has_matrix |= a.has_matrix();
                }
                if matches!(f,Func::Integral|Func::Sum|Func::Product) {features|=4;}
                if matches!(f,Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction)){features|=8;}
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
                        features |= scalar.0.features;
                    }
                    let _ = parent.scalar_map(&mut |scalar| {
                        mask |= scalar.mask();
                        features |= scalar.0.features;
                        scalar.clone()
                    });
                }
                if let Func::MatrixDotProduct(a,b)=f {
                    has_matrix=true;
                    for parent in [a,b] {
                        for matrix in parent.matrix_symbols(){let mut atom_hash=FxHasher::default();matrix.hash(&mut atom_hash);mask|=1u64<<(atom_hash.finish()%64);}
                        let _=parent.scalar_map(&mut |scalar|{mask|=scalar.mask();features|=scalar.0.features;scalar.clone()});
                    }
                }
            }
            Kind::Derivative(e, vars) => {
                h.write_u64(e.hash64());
                mask = e.mask();
                features |= e.0.features;
                has_matrix = e.has_matrix();
                for (v, n) in vars {
                    h.write_u64(v.hash64());
                    n.hash(&mut h);
                    mask |= v.mask();
                    features |= v.0.features;
                    has_matrix |= v.has_matrix();
                }
            }
            Kind::Rel(op, a, b) => {
                op.hash(&mut h);
                h.write_u64(a.hash64());
                h.write_u64(b.hash64());
                mask = a.mask() | b.mask();
                features |= a.0.features | b.0.features;
                has_matrix = a.has_matrix() || b.has_matrix();
            }
            Kind::Piecewise(pieces) => {
                for (e, c) in pieces {
                    h.write_u64(e.hash64());
                    h.write_u64(c.hash64());
                    mask |= e.mask() | c.mask();
                    features |= e.0.features | c.0.features;
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
            features,
            binding_identity,
            free_symbols: OnceLock::new(),
            ordered_args: OnceLock::new(),
            float_info: OnceLock::new(),
            assumption_cache: AtomicU64::new(0),
        }))
    }

    #[inline]
    pub fn kind(&self) -> &Kind {
        &self.0.kind
    }

    /// Identity of this live immutable allocation, for graph codecs/memoization.
    pub fn dag_id(&self) -> usize { Arc::as_ptr(&self.0) as usize }
    /// A non-owning reference to this node.
    pub fn downgrade(&self) -> WeakExpr { WeakExpr(Arc::downgrade(&self.0)) }

    /// Each shared node is inspected once. Compact coefficients and matrix
    /// metadata are visited directly, without synthesizing public arguments.
    pub fn float_info(&self)->FloatInfo {
        *self.0.float_info.get_or_init(||{
            let mut info=FloatInfo::default();
            match self.kind(){
                Kind::Num(n)=>info=FloatInfo::number(n),
                Kind::Add(c,terms)=>{info.include(FloatInfo::number(c));for (e,c) in terms{info.include(e.float_info());info.include(FloatInfo::number(c));}},
                Kind::Mul(c,factors)=>{info.include(FloatInfo::number(c));for (b,e) in factors{info.include(b.float_info());info.include(e.float_info());}},
                Kind::Pow(b,e)|Kind::Rel(_,b,e)=>{info.include(b.float_info());info.include(e.float_info());},
                Kind::Func(f,args)=>{
                    for arg in args{info.include(arg.float_info());}
                    if let Func::MatrixElement(matrix)|Func::MatrixProxy(matrix)|Func::MatrixTrace(matrix)|Func::MatrixDeterminant(matrix)=f {
                        let mut fields=Vec::new();matrix.scalar_fields(&mut fields);
                        for field in fields{info.include(field.float_info());}
                    }
                    if let Func::MatrixDotProduct(a,b)=f {for matrix in [a,b]{let mut fields=Vec::new();matrix.scalar_fields(&mut fields);for field in fields{info.include(field.float_info());}}}
                },
                Kind::Derivative(e,variables)=>{info.include(e.float_info());for (v,_) in variables{info.include(v.float_info());}},
                Kind::Piecewise(parts)=>{for (e,c) in parts{info.include(e.float_info());info.include(c.float_info());}},
                _=>{},
            }
            info
        })
    }

    #[inline]
    pub(crate) fn cached_assumption(&self,shift:u32,width:u32)->u8 {
        ((self.0.assumption_cache.load(AtomicOrdering::Relaxed)>>shift)&((1u64<<width)-1)) as u8
    }

    /// Publish only a complete fact, without disturbing concurrently written
    /// fields. In particular, an unknown answer is distinct from an empty slot.
    pub(crate) fn remember_assumption(&self,shift:u32,width:u32,value:u8) {
        let mask=((1u64<<width)-1)<<shift;
        debug_assert!(shift+width<=64 && value!=0 && (value as u64)<(1u64<<width));
        let _=self.0.assumption_cache.fetch_update(AtomicOrdering::Relaxed,AtomicOrdering::Relaxed,|old|{
            if old&mask==0{Some(old|((value as u64)<<shift))}else{None}
        });
    }

    /// Whether this scalar expression contains a typed matrix entry.
    #[inline]
    pub fn has_matrix(&self) -> bool { self.0.has_matrix }

    /// A semantic pattern with alpha-bound names cannot use a structural mask.
    pub fn has_binding_identity(&self)->bool {self.0.features & 1 != 0}
    pub fn has_calculus_binding(&self)->bool {self.0.features & 4 != 0}
    pub fn has_lambda_binding(&self)->bool {self.0.features & 8 != 0}
    /// Whether a scalar descendant preserves explicitly held arithmetic.
    pub fn has_held_arithmetic(&self)->bool {self.0.features & 16 != 0}
    pub fn has_nonfinite(&self)->bool {self.0.features & 32 != 0}
    pub(crate) fn has_wide_exponent(&self)->bool {self.0.features & 64 != 0}
    pub fn has_wild_properties(&self)->bool {self.0.features & 128 != 0}
    /// Whether an atom or function head of a host-language subclass occurs.
    pub fn has_host(&self)->bool {self.0.features & 256 != 0}
    /// Whether this node itself is an atom or application of a host class.
    pub fn is_host(&self)->bool {
        self.has_host() && match &self.0.kind {
            Kind::Sym(s)=>s.host().is_some(),
            Kind::Func(Func::Undefined(head),_)=>head.host().is_some(),
            _=>false,
        }
    }
    pub fn has_limit(&self)->bool {self.0.features & 2 != 0}
    /// The host handle of a host atom or application.
    pub fn host_handle(&self)->Option<&Arc<dyn crate::host::HostClass>> {
        if !self.has_host() {return None;}
        match &self.0.kind {
            Kind::Sym(s)=>s.host(),
            Kind::Func(Func::Undefined(head),_)=>head.host(),
            _=>None,
        }
    }
    /// Whether a host node below intercepts substitution (`_eval_subs`), so
    /// that `subs` must visit it even where the replaced pattern is absent.
    pub fn has_host_subs(&self)->bool {self.0.features & 512 != 0}
    /// Whether a host node below may fail to commute (SymPy
    /// `is_commutative` not True), so products must keep its order.
    #[inline]
    pub fn has_noncommutative_host(&self)->bool {self.0.features & 1024 != 0}
    /// The number of strong references to this node.
    pub fn strong_count(&self)->usize {Arc::strong_count(&self.0)}
    /// Visit the nodes this node references strongly through its scalar
    /// payload and its derived caches, once per reference. References held
    /// by other payload (matrix parents, binding keys, atom metadata) are not
    /// reported, so a caller counting owners sees an over-approximation of
    /// the external references of a node, never an under-approximation.
    pub fn for_each_owned<'a>(&'a self, f:&mut dyn FnMut(&'a Expr)) {
        match &self.0.kind {
            Kind::Num(_)|Kind::Sym(_)|Kind::Const(_)|Kind::Bool(_)=>{},
            Kind::Add(_,terms)=>for (t,_) in terms {f(t)},
            Kind::Mul(_,factors)=>for (b,e) in factors {f(b);f(e)},
            Kind::Pow(b,e)|Kind::Rel(_,b,e)=>{f(b);f(e)},
            Kind::Func(_,args)=>for a in args {f(a)},
            Kind::Derivative(e,vars)=>{f(e);for (v,_) in vars {f(v)}},
            Kind::Piecewise(pieces)=>for (e,c) in pieces {f(e);f(c)},
        }
        if let Some(symbols)=self.0.free_symbols.get() {for s in symbols {f(s)}}
        if let Some(args)=self.0.ordered_args.get() {for a in args {f(a)}}
    }

    /// Conservative unchanged-child check for transforms preserving public args.
    pub fn same_presentation(&self,other:&Expr)->bool {
        self.ptr_eq(other) || (!self.has_binding_identity() && !other.has_binding_identity()
            && !self.has_held_arithmetic() && !other.has_held_arithmetic() && self==other)
    }

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

    /// A fresh Dummy atom. Its display name does not determine identity.
    pub fn dummy(name: &str, a: crate::symbol::Assumptions) -> Expr {
        Expr::from_kind(Kind::Sym(crate::symbol::Symbol::dummy(name, a)))
    }

    pub fn dummy_with_index(name: &str, a: crate::symbol::Assumptions, index: num_bigint::BigInt) -> Expr {
        Expr::from_kind(Kind::Sym(crate::symbol::Symbol::dummy_with_index(name, a, index)))
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
            Kind::Sym(s) if s.is_dummy() => "Dummy",
            Kind::Sym(s) if s.is_wild() => "Wild",
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
                Constant::Catalan => "Catalan",
                Constant::GoldenRatio => "GoldenRatio",
                Constant::TribonacciConstant => "TribonacciConstant",
            },
            Kind::Add(..) => "Add",
            Kind::Mul(..) => "Mul",
            Kind::Pow(..) => "Pow",
            Kind::Func(Func::AppliedPredicate, a) if a.first().and_then(crate::ask::as_predicate).is_some_and(|p| p.is_binrel()) => "AppliedBinaryRelation",
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
        if self.ptr_eq(o) { return Ordering::Equal; }
        let rank=self.kind().rank().cmp(&o.kind().rank());
        if rank!=Ordering::Equal { return rank; }
        if self.is_atom() { return self.cmp_struct_uncached(o); }
        let _guard=ComparisonGuard::order();
        let pair=(self.dag_id(),o.dag_id());
        if let Some(result)=ORDER_MEMO.with(|memo|memo.borrow().values.as_ref().and_then(|values|values.get(&pair)).copied()) { return result; }
        let result=self.cmp_struct_uncached(o);
        ORDER_MEMO.with(|memo|memo.borrow_mut().values.get_or_insert_with(FxHashMap::default).insert(pair,result));
        result
    }

    fn cmp_struct_uncached(&self, o: &Expr) -> Ordering {
        let (a,b)=(self.kind(),o.kind());
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
                if *f1==Func::Subs && *f2==Func::Subs && (self.0.binding_identity.is_some() || o.0.binding_identity.is_some()) {
                    self.0.binding_identity.cmp(&o.0.binding_identity)
                }else{f1.cmp(f2).then_with(|| cmp_slices(a1, a2))}
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
        if self.ptr_eq(o) { return true; }
        if self.0.hash!=o.0.hash { return false; }
        if self.is_atom() { return kind_eq(self.kind(),o.kind()); }
        let _guard=ComparisonGuard::equality();
        let pair=(self.dag_id(),o.dag_id());
        if let Some(result)=EQ_MEMO.with(|memo|memo.borrow().values.as_ref().and_then(|values|values.get(&pair)).copied()) { return result; }
        let result=match (&self.0.binding_identity,&o.0.binding_identity) {
            (Some(a),Some(b))=>a==b,
            (None,None)=>kind_eq(self.kind(),o.kind()),
            _=>false,
        };
        EQ_MEMO.with(|memo|memo.borrow_mut().values.get_or_insert_with(FxHashMap::default).insert(pair,result));
        result
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

#[cfg(test)]
mod dag_comparison_tests {
    use super::*;

    fn pair_graph(depth: usize) -> Expr {
        let mut value=Expr::symbol("leaf");
        for _ in 0..depth {
            value=Expr::from_kind(Kind::Func(Func::Undefined("pair".into()),vec![value.clone(),value]));
        }
        value
    }

    #[test]
    fn scalar_dag_eq_and_both_orders_visit_shared_pairs_once() {
        let first=pair_graph(100);let restored=pair_graph(100);
        assert!(first==restored);
        assert_eq!(first.cmp_struct(&restored),Ordering::Equal);
        assert_eq!(first.cmp(&restored),Ordering::Equal);
        // Multiple outer calls must start with independent, freshly owned memos.
        assert!(restored==pair_graph(100));
        assert!(first!=pair_graph(99));
    }

    #[test]
    fn comparison_context_survives_scalar_matrix_boundaries() {
        fn mixed_graph(depth:usize)->Expr {
            let mut value=Expr::symbol("leaf");
            for _ in 0..depth {
                let matrix=crate::matrix::Matrix::new(2,2,vec![value.clone(),value.clone(),value.clone(),value]).unwrap();
                let parent=crate::matrixexpr::MatrixExpr::explicit(matrix);
                value=Expr::from_kind(Kind::Func(Func::MatrixTrace(Arc::new(parent)),vec![]));
            }
            value
        }
        let first=mixed_graph(24);let restored=mixed_graph(24);
        assert!(first==restored);
        assert_eq!(first.cmp_struct(&restored),Ordering::Equal);
        assert_eq!(first.cmp(&restored),Ordering::Equal);
    }
}
