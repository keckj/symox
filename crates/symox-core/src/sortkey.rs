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
use crate::printing::order::decompose_power;
use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand};

/// Presentation arguments deliberately differ from scalar `Expr::args()`:
/// matrix parents live in typed payloads, not in commutative scalar children.
pub(crate) enum PresentationArg {
    Scalar(Expr),
    Matrix(MatrixExpr),
    Text(String),
    Permutation(String),
    None,
    Node(&'static str, Vec<PresentationArg>),
}

pub(crate) struct MatrixPresentation {
    pub class: (u8, u32, &'static str),
    pub args: Vec<PresentationArg>,
}

pub fn matrix_class_key(matrix:&MatrixExpr)->Key{
    let(a,b,name)=match matrix.kind(){
        MatrixKind::Add(_)=>(3,1,"MatAdd"),MatrixKind::Mul(..)=>(3,0,"MatMul"),
        MatrixKind::Explicit(_)=>(5,0,"ImmutableDenseMatrix"),MatrixKind::Conjugate(_)=>(5,0,"Adjoint"),
        _=>(5,0,matrix.kind_name()),
    };
    tup(vec![Key::Int(a),Key::Int(b),Key::Str(name.into())])
}

pub(crate) fn matrix_presentation(matrix: &MatrixExpr) -> MatrixPresentation {
    use PresentationArg::{Matrix as M, Scalar as S};
    let mut class = (5, 0, matrix.kind_name());
    let args = match matrix.kind() {
        MatrixKind::Permutation(p) => vec![PresentationArg::Permutation(crate::matrixexpr::permutation_label(p))],
        MatrixKind::Permute(parent,p,axis) => vec![M(parent.clone()),
            PresentationArg::Permutation(crate::matrixexpr::permutation_label(p)),S(Expr::int(*axis as i64))],
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
        MatrixKind::ApplyFunction(function,parent) => vec![S(function.clone()),M(parent.clone())],
        MatrixKind::Function(lambda) => vec![S(matrix.rows().clone()), S(matrix.cols().clone()), S(lambda.clone())],
        MatrixKind::Identity | MatrixKind::Dft(_) => vec![S(matrix.rows().clone())],
        MatrixKind::GenericIdentity | MatrixKind::GenericZero => vec![],
        MatrixKind::Zero | MatrixKind::One => vec![S(matrix.rows().clone()), S(matrix.cols().clone())],
        MatrixKind::HadamardProduct(factors) | MatrixKind::Kronecker(factors) | MatrixKind::BlockDiag(factors) => factors.iter().cloned().map(M).collect(),
        MatrixKind::Slice(parent, rows, cols) => vec![M(parent.clone()),
            PresentationArg::Node("Tuple", rows.args().into_iter().map(|e| e.map(S).unwrap_or(PresentationArg::None)).collect()),
            PresentationArg::Node("Tuple", cols.args().into_iter().map(|e| e.map(S).unwrap_or(PresentationArg::None)).collect())],
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
        MatrixKind::Transpose(base) | MatrixKind::Adjoint(base)
        | MatrixKind::Diag(base) | MatrixKind::DiagonalOf(base) | MatrixKind::Diagonal(base) => vec![M(base.clone())],
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
    SortKeyCache::default().matrix_key(matrix)
}

/// A Python-tuple like sort key.
#[derive(Clone, Debug)]
pub enum Key {
    Int(i64),
    Index(num_bigint::BigInt),
    Str(String),
    Num(Number),
    SpecialNumber(Constant),
    Tup(Vec<Key>),
    Shared(Arc<[Key]>),
}

impl Key {
    fn rank(&self) -> u8 {
        match self {
            Key::Int(_) | Key::Index(_) => 0,
            Key::Num(_) | Key::SpecialNumber(_) => 0,
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

fn special_rank(c:Constant)->i8{match c{Constant::NegInfinity=>-1,Constant::Infinity=>1,_=>2}}
fn cmp_special(c:Constant,n:&Number)->Ordering{
    let rank=special_rank(c);
    if let Number::Float(f)=n{
        let numeric=if f.is_nan(){Some(2)}else if *f==f64::NEG_INFINITY{Some(-1)}else if *f==f64::INFINITY{Some(1)}else{None};
        if let Some(other)=numeric{return rank.cmp(&other);}
    }
    rank.cmp(&0)
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
            (Key::SpecialNumber(a),Key::SpecialNumber(b))=>special_rank(*a).cmp(&special_rank(*b)),
            (Key::SpecialNumber(a),Key::Num(b))=>cmp_special(*a,b),
            (Key::Num(a),Key::SpecialNumber(b))=>cmp_special(*b,a).reverse(),
            (Key::SpecialNumber(a),Key::Int(b))=>cmp_special(*a,&Number::Int(*b)),
            (Key::Int(a),Key::SpecialNumber(b))=>cmp_special(*b,&Number::Int(*a)).reverse(),
            (Key::SpecialNumber(a),Key::Index(b))=>cmp_special(*a,&Number::from_bigint(b.clone())),
            (Key::Index(a),Key::SpecialNumber(b))=>cmp_special(*b,&Number::from_bigint(a.clone())).reverse(),
            (Key::Int(a), Key::Int(b)) => a.cmp(b),
            (Key::Index(a), Key::Index(b)) => a.cmp(b),
            (Key::Index(a), Key::Int(b)) => a.cmp(&(*b).into()),
            (Key::Int(a), Key::Index(b)) => num_bigint::BigInt::from(*a).cmp(b),
            (Key::Index(a), Key::Num(b)) => cmp_num(&Number::from_bigint(a.clone()), b),
            (Key::Num(a), Key::Index(b)) => cmp_num(a, &Number::from_bigint(b.clone())),
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

/// SymPy `Function.class_key` ordinal of an undefined head: the known names
/// first, then 0 for any argument count (`Naturals0`) and 10000 otherwise.
/// SymPy 1.14 function classes whose `nargs` is finite (generated from
/// `_function_nargs.NARGS`): their held instances order with 10000.
const FIXED_NARGS_CLASSES: &[&str] = &[
    "Abs", "BesselBase", "Chi", "Ci", "Coeff", "DiracDelta", "Ei", "FallingFactorial",
    "FresnelIntegral", "Heaviside", "KroneckerDelta", "LambertW", "Li", "Mod", "Not",
    "ReciprocalHyperbolicFunction", "ReciprocalTrigonometricFunction", "Rem", "RisingFactorial",
    "RoundFunction", "Shi", "Si", "SingularityFunction", "SphericalBesselBase",
    "SphericalHankelBase", "TrigonometricIntegral", "Ynm", "Znm", "_erfs", "_rf", "acos", "acosh",
    "acot", "acoth", "acsc", "acsch", "airyai", "airyaiprime", "airybi", "airybiprime", "andre",
    "appellf1", "arg", "asec", "asech", "asin", "asinh", "assoc_laguerre", "assoc_legendre", "atan",
    "atan2", "atanh", "bell", "bernoulli", "besseli", "besselj", "besselk", "bessely", "beta",
    "betainc", "betainc_regularized", "binomial", "catalan", "ceiling", "chebyshevt",
    "chebyshevt_root", "chebyshevu", "chebyshevu_root", "conjugate", "cos", "cosh", "cot", "coth",
    "csc", "csch", "digamma", "dirichlet_eta", "divisor_sigma", "elliptic_e", "elliptic_f",
    "elliptic_k", "elliptic_pi", "erf", "erf2", "erf2inv", "erfc", "erfcinv", "erfi", "erfinv",
    "euler", "exp", "expint", "factorial", "factorial2", "fibonacci", "floor", "frac", "fresnelc",
    "fresnels", "gamma", "gegenbauer", "genocchi", "hankel1", "hankel2", "harmonic", "hermite",
    "hermite_prob", "hn1", "hn2", "hyper", "im", "jacobi", "jacobi_symbol", "jn",
    "kronecker_symbol", "laguerre", "legendre", "legendre_symbol", "li", "log", "loggamma",
    "lowergamma", "lucas", "marcumq", "mathieuc", "mathieucprime", "mathieus", "mathieusprime",
    "mobius", "motzkin", "multigamma", "partition", "periodic_argument", "polar_lift", "polygamma",
    "polylog", "primenu", "primeomega", "primepi", "principal_branch", "re", "reduced_totient",
    "riemann_xi", "sec", "sech", "sign", "sin", "sinc", "sinh", "stieltjes", "subfactorial", "tan",
    "tanh", "totient", "transpose", "tribonacci", "trigamma", "udivisor_sigma", "uppergamma", "yn",
    "zeta",
];

pub(crate) fn undefined_ordinal(head: &crate::expr::UndefinedFunction) -> u32 {
    match &*head.name {
        "exp" => 10, "log" => 11,
        "sin" => 20, "cos" => 21, "tan" => 22, "cot" => 23,
        "sinh" => 30, "cosh" => 31, "tanh" => 32, "coth" => 33,
        "conjugate" => 40, "re" => 41, "im" => 42, "arg" => 43,
        // A held instance of a fixed-arity SymPy class carries no nargs here.
        name if head.nargs.is_some() || FIXED_NARGS_CLASSES.binary_search(&name).is_ok() => 10000,
        _ => 0,
    }
}

pub fn class_key(e: &Expr) -> Key {
    if let Kind::Func(Func::MatrixProxy(matrix),_)=e.kind(){return matrix_class_key(matrix);}
    let (a, b, n): (i64, i64, String) = match e.kind() {
        Kind::Num(_) => (1, 0, "Number".into()),
        Kind::Const(c) => match c {
            Constant::Infinity | Constant::NegInfinity | Constant::NaN => {
                (1, 0, "Number".into())
            }
            _ => (2, 0, e.class_name().into()),
        },
        // Atom.class_key: the class name, which a host subclass replaces.
        Kind::Sym(_) => (2, 0, crate::host::class_name(e).unwrap_or(e.class_name()).into()),
        Kind::Bool(_) => (5, 0, e.class_name().into()),
        Kind::Pow(..) => (3, 2, "Pow".into()),
        Kind::Mul(..) => (3, 0, "Mul".into()),
        Kind::Add(..) => (3, 1, "Add".into()),
        Kind::Func(f @ (Func::MatrixElement(_)|Func::MatrixTrace(_)|Func::MatrixDeterminant(_)|Func::MatrixDotProduct(..)),_)=>(5,0,f.name().into()),
        // Max and Min are not Function subclasses: Basic.class_key.
        Kind::Func(f @ (Func::SymbolicDerivative | Func::Logic(_) | Func::Set(_) | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::Order | Func::Max | Func::Min), _) => (5, 0, f.name().into()),
        Kind::Func(Func::Undefined(head), args) if crate::poly::rootof::is_rootof_head(head, args.len()) => (5, 0, "ComplexRootOf".into()),
        // RootSum is an Expr subclass, not a Function: Basic.class_key.
        Kind::Func(Func::Undefined(_), _) if crate::printer::root_sum_parts(e).is_some() => (5, 0, "RootSum".into()),
        Kind::Func(Func::Undefined(head), _) if head.internal => (5, 0, head.name.to_string()),
        // An Expr or Atom subclass that is not a Function keeps its class_key.
        Kind::Func(Func::Undefined(head), _) if head.host().and_then(|host| host.class_key()).is_some() => {
            let (a, b) = head.host().and_then(|host| host.class_key()).unwrap();
            (a.into(), b.into(), head.name.to_string())
        }
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
                Func::Undefined(head) => undefined_ordinal(head).into(),
                // Classes with neither `nargs` nor `eval` take any number
                // of arguments (`Naturals0`).
                Func::Max | Func::Min | Func::ExpPolar => 0,
                Func::Special(s) if s.is_transform() || matches!(s, crate::special::SpecialFn::MeijerG
                    | crate::special::SpecialFn::Lerchphi | crate::special::SpecialFn::LeviCivita
                    | crate::special::SpecialFn::Carmichael) => 0,
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
    SortKeyCache::default().key(e)
}

/// Each scalar and matrix node is keyed once; tuple children retain sharing.
#[derive(Default)]
struct SortKeyCache {
    // Presentation depends on the stored bound names, even when two Subs
    // nodes have equal alpha identities. Retain owners to prevent ID reuse.
    scalars: FxHashMap<usize, (Expr, Key)>,
    matrices: FxHashMap<usize, (MatrixExpr, Key)>,
}

impl SortKeyCache {
    fn key(&mut self, e: &Expr) -> Key {
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
            PresentationArg::Permutation(s) => composite_key((2, 0, "Permutation"), vec![Key::Str(s)]),
            PresentationArg::None => composite_key((0, 0, "NoneType"), vec![Key::Str("None".into())]),
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
        Kind::Bool(_)=>return composite_key((5,0,e.class_name()),vec![]),
        // Heaviside(x).args is (x, 1/2) in SymPy, so it sorts as the
        // explicit two-argument node
        Kind::Func(Func::Heaviside,args) if args.len()==1 => {
            return self.key(&Expr::from_kind(Kind::Func(Func::Heaviside,vec![args[0].clone(),Expr::rational(1,2)])));
        }
        Kind::Func(Func::Set(crate::sets::SetOp::Lambda),args) if args.len()==2 => {
            let signature=crate::sets::lambda_signature(e).unwrap();
            return composite_key((5,0,"Lambda"),vec![self.key(&signature),self.key(&args[1])]);
        }
        Kind::Func(Func::Set(crate::sets::SetOp::IdentityFunction),_)=>return composite_key((5,0,"IdentityFunction"),vec![]),
        Kind::Func(Func::SymbolicDerivative,args)=>return composite_key((5,0,"Derivative"),args.iter().map(|arg|self.key(arg)).collect()),
        Kind::Derivative(expression, variables) => {
            let mut args = vec![self.key(expression)];
            for (variable, count) in variables {
                let args_variable = vec![self.key(variable), number_key(Number::Int(*count as i64))];
                args.push(composite_key((5, 0, "Tuple"), args_variable));
            }
            return composite_key((5, 0, "Derivative"), args);
        }
        Kind::Func(Func::Order,arguments) if !arguments.is_empty()=>{
            let mut args=vec![self.key(&arguments[0])];
            for pair in arguments[1..].chunks(2) {
                let keys=pair.iter().map(|arg|self.key(arg)).collect();
                args.push(composite_key((5,0,"Tuple"),keys));
            }
            return composite_key((5,0,"Order"),args);
        }
        Kind::Piecewise(pieces)=>{
            let args=pieces.iter().map(|(expression,condition)|{
                let keys=vec![self.key(expression),self.key(condition)];
                composite_key((5,0,"ExprCondPair"),keys)
            }).collect();
            return composite_key((4,0,"Piecewise"),args);
        }
        Kind::Func(Func::MatrixProxy(matrix), _) => return self.matrix_key(matrix),
        Kind::Func(Func::MatrixDotProduct(left,right),_) => {
            let args=vec![self.matrix_key(left),self.matrix_key(right)];
            return composite_key((5,0,"DotProduct"),args);
        }
        Kind::Func(Func::MatrixElement(matrix), args) => {
            let mut keys = vec![self.matrix_key(matrix)];
            keys.extend(args.iter().map(|arg| self.key(arg)));
            return composite_key((5, 0, "MatrixElement"), keys);
        }
        Kind::Func(f @ (Func::MatrixTrace(matrix) | Func::MatrixDeterminant(matrix)), _) => {
            return composite_key((5, 0, f.name()), vec![self.matrix_key(matrix)]);
        }
        Kind::Num(n) => return number_key(n.clone()),
        Kind::Const(c @ (Constant::Infinity|Constant::NegInfinity|Constant::NaN))=>{
            let mut fields=number_key(Number::ONE).tuple().unwrap().to_vec();fields[3]=Key::SpecialNumber(*c);return tup(fields);
        }
        _ => {}
    }
    let (coeff,rest)=crate::decompose::as_coeff_mul_head(e,false);
    let coefficient=match coeff.kind(){Kind::Num(n)=>Key::Num(n.clone()),Kind::Const(c)=>Key::SpecialNumber(*c),_=>unreachable!("a numeric coefficient")};
    let (base, exp) = match rest.kind() {
        Kind::Pow(..)=>{
            let(b,x)=crate::decompose::as_base_exp(&rest);
            if b.is_constant(Constant::E){(crate::functions::undefined("exp",vec![x]),Expr::one())}else{(b,x)}
        },
        _ => (rest.clone(), Expr::one()),
    };
    if matches!(base.kind(), Kind::Derivative(..) | Kind::Func(Func::SymbolicDerivative,_) | Kind::Piecewise(_) | Kind::Func(Func::Order | Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction),_) | Kind::Func(Func::MatrixProxy(_) | Func::MatrixElement(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..), _)) {
        let mut parts = match self.key(&base) {
            Key::Tup(parts) => parts,
            Key::Shared(parts) => parts.to_vec(),
            _ => unreachable!("expression sort keys are tuples"),
        };
        parts[2] = self.key(&exp);
        parts[3] = coefficient;
        return tup(parts);
    }
    let args: Vec<Key> = if let Kind::Sym(symbol) = base.kind() {
        match symbol.dummy_index() {
            // Dummy.sort_key has its own atomic representation, whereas
            // Expr.sort_key wraps a Dummy base in its complete atomic key.
            Some(_) if !e.is_symbol() => vec![self.key(&base)],
            Some(index) => vec![Key::Str(crate::host::atom_name(&base, symbol).into_owned()), Key::Index(index.clone())],
            None => vec![Key::Str(crate::host::atom_name(&base, symbol).into_owned())],
        }
    } else if base.is_atom() || crate::host::expression_head(&base).is_some_and(|host| host.atomic()) {
        vec![Key::Str(crate::printer::to_str(&base))]
    } else {
        match base.kind() {
            Kind::Add(..) if !base.is_canonical() => self.held_ordered_terms(&base).iter().map(|arg| self.key(arg)).collect(),
            Kind::Add(..) => crate::printing::order::ordered_terms(&base).iter().map(|arg| self.key(arg)).collect(),
            Kind::Mul(..) => self.ordered_factors(&base).iter().map(|arg| self.key(arg)).collect(),
            // RootSum.args is the reference `(expr, fun, poly.gen)`.
            Kind::Func(Func::Undefined(_), _) if crate::printer::root_sum_parts(&base).is_some() => {
                let (gen, poly, fun) = crate::printer::root_sum_parts(&base).unwrap();
                [poly, fun, gen].iter().map(|arg| self.key(arg)).collect()
            }
            _ => base.args().iter().map(|arg| self.key(arg)).collect(),
        }
    };
    tup(vec![
        class_key(&base),
        tup(vec![Key::Int(args.len() as i64), tup(args)]),
        self.key(&exp),
        coefficient,
    ])
    }
}

/// SymPy's `as_ordered_factors` for a product (numbers first).
pub fn ordered_factors(e: &Expr) -> Vec<Expr> {
    SortKeyCache::default().ordered_factors(e)
}

impl SortKeyCache {
    fn ordered_factors(&mut self, e: &Expr) -> Vec<Expr> {
        if !matches!(e.kind(), Kind::Mul(..)) { return vec![e.clone()]; }
        let (commutative, noncommutative) = crate::decompose::args_cnc(e, true);
        let mut keyed: Vec<_> = commutative.into_iter().map(|f| (self.key(&f), f)).collect();
        keyed.sort_by(|a,b| a.0.cmp(&b.0));
        keyed.into_iter().map(|(_, f)| f).chain(noncommutative).collect()
    }
}

/// SymPy's `as_ordered_terms` (default `lex` order) for a sum.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    crate::printing::order::ordered_terms(e)
}

impl SortKeyCache {
/// `as_ordered_terms` of a held sum: `Add.make_args` keeps the held
/// arguments, numeric subtrees evaluate to complex coefficients.
pub(crate) fn held_ordered_terms(&mut self, e: &Expr) -> Vec<Expr> {
    let args = match e.kind() {
        Kind::Add(..) => e.args(),
        _ => return vec![e.clone()],
    };
    // The positive numeric term precedes a negative single-factor term.
    if args.len() == 2 {
        let numeric = |a:&Expr| matches!(a.kind(),Kind::Num(_)|Kind::Const(Constant::Infinity|Constant::NegInfinity|Constant::NaN)) || is_number_symbol(a);
        let (numbers, other):(Vec<_>,Vec<_>)=args.iter().partition(|a|numeric(a));
        if numbers.len()==1 && other.len()==1 && other[0].is_mul() {
            let mut factors=other[0].args();factors.sort_by_key(|a|!numeric(a));
            let positive=match numbers[0].kind(){Kind::Num(n)=>n.is_positive(),_=>is_number_symbol(numbers[0])};
            if factors.len()==2 && positive && factors[0].as_number().is_some_and(Number::is_negative){return vec![numbers[0].clone(),other[0].clone()];}
        }
    }
    struct T {
        term: Expr,
        coeff: Complex64,
        cpart: Vec<(Expr, i64)>,
        ncpart: Vec<Key>,
    }
    let mut terms: Vec<T> = Vec::with_capacity(args.len());
    let mut gens: Vec<Expr> = Vec::new();
    for t in args {
        let (c, rest) = crate::decompose::as_coeff_mul_head(&t, false);
        let mut coeff = numeric_coefficient(&c).unwrap_or(Complex64::new(1.0,0.0));
        let mut cpart: Vec<(Expr, i64)> = Vec::new();
        let mut ncpart=Vec::new();
        if !rest.is_one() {
            let factors = match rest.kind() {
                Kind::Mul(..) => rest.args(),
                _ => vec![rest.clone()],
            };
            for f in factors {
                if crate::symbolic_number::is_symbolic_number(&f) {
                    if let Some(z) = numeric_coefficient(&f) {
                        coeff *= z;
                        continue;
                    }
                }
                if crate::decompose::commutativity(&f)!=Some(true){ncpart.push(self.key(&f));continue;}
                let (b, k) = decompose_power(&f);
                if !gens.contains(&b) {
                    gens.push(b.clone());
                }
                cpart.push((b, k));
            }
        }
        terms.push(T { term: t, coeff, cpart, ncpart });
    }
    let mut gkeys: Vec<(Key, Expr)> = gens.into_iter().map(|g| (self.key(&g), g)).collect();
    gkeys.sort_by(|a, b| a.0.cmp(&b.0));
    let gens: Vec<Expr> = gkeys.into_iter().map(|x| x.1).collect();
    type TermKey=(Vec<i64>,Vec<Key>,(bool,f64,f64),Expr);
    let mut keyed: Vec<TermKey> = terms
        .into_iter()
        .map(|t| {
            let mut monom = vec![0i64; gens.len()];
            for (b, k) in &t.cpart {
                let i = gens.iter().position(|g| g == b).unwrap();
                monom[i] = *k;
            }
            (monom, t.ncpart, (t.coeff.im != 0.0, t.coeff.im, t.coeff.re), t.term)
        })
        .collect();
    let compare=|a:&TermKey,b:&TermKey| {
        // descending lex on monomials, then ascending coefficient
        b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)).then_with(|| {
            let (ai, aim, are) = a.2;
            let (bi, bim, bre) = b.2;
            ai.cmp(&bi).then_with(||{
                // Python tuple comparison stops at unequal NaNs; it does not
                // continue with the real coefficient when neither is less.
                for(x,y)in[(aim,bim),(are,bre),(aim,bim)]{
                    if x!=y{return x.partial_cmp(&y).unwrap_or(Ordering::Equal);}
                }
                Ordering::Equal
            })
        })
    };
    if keyed.iter().any(|entry|entry.2.1.is_nan()||entry.2.2.is_nan()){
        // NaN comparison is partial. A stable insertion pass preserves its
        // Python semantics without imposing Rust sort's total-order contract.
        for i in 1..keyed.len(){let mut j=i;while j>0 && compare(&keyed[j],&keyed[j-1]).is_lt(){keyed.swap(j,j-1);j-=1;}}
    }else{keyed.sort_by(compare);}
    keyed.into_iter().map(|x| x.3).collect()
}
}

fn is_number_symbol(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(c) if c.is_number_symbol())
}

/// `as_ordered_terms` of a held sum (see `SortKeyCache::held_ordered_terms`).
pub(crate) fn held_ordered_terms(e: &Expr) -> Vec<Expr> {
    SortKeyCache::default().held_ordered_terms(e)
}

pub(crate) fn numeric_coefficient(e:&Expr)->Option<Complex64>{
    if matches!(e.kind(),Kind::Const(Constant::NaN|Constant::ComplexInfinity)){return Some(Complex64::new(f64::NAN,f64::NAN));}
    if let Some(value)=crate::evalf::eval_c64(e){return Some(value);}
    // Closed held arithmetic can contain zoo/nan, which the finite evaluator
    // deliberately refuses. Normalize those arithmetic nodes for this numeric
    // presentation only; never feed this fallback into assumption inference.
    if e.has_nonfinite(){
        let value=match e.kind(){
            Kind::Add(..)=>crate::ops::add(&e.args()),
            Kind::Mul(..)=>crate::ops::mul(&e.args()),
            Kind::Pow(base,exponent)=>crate::ops::pow(base,exponent),
            _=>return None,
        };
        if value!=*e{return numeric_coefficient(&value);}
    }
    None
}

thread_local!{
    // `complex(factor)` by value: the same radicals and denominators recur
    // in every term of a printed sum (as fresh nodes, hence the structural
    // key); the table is dropped when it grows past the bound.
    static COMPLEX_MEMO:std::cell::RefCell<rustc_hash::FxHashMap<Expr,Option<Complex64>>>=std::cell::RefCell::new(rustc_hash::FxHashMap::default());
}
const COMPLEX_MEMO_BOUND:usize=1<<12;

/// `complex(factor)` as `Expr.__complex__` computes it, an evalf at 53 bits:
/// a radical such as `31**(2/3)` gets its correctly rounded double, where
/// `numeric_coefficient` (the double evaluator) can be one unit off. Falls
/// back to `numeric_coefficient` outside the multiprecision evaluator.
pub(crate) fn correctly_rounded_coefficient(e:&Expr)->Option<Complex64>{
    if matches!(e.kind(),Kind::Const(Constant::NaN|Constant::ComplexInfinity)){return Some(Complex64::new(f64::NAN,f64::NAN));}
    correctly_rounded(e).or_else(||numeric_coefficient(e))
}

/// The 53-bit evalf of a numeric node as a complex double (`None` when the
/// multiprecision evaluator does not cover it).
///
/// SymPy's evalf keeps the imaginary rounding noise of a real value reached
/// through complex arithmetic (`complex(1/sqrt(2/(3*c) + 2*c))` has an
/// imaginary part of -1e-28 for a complex `c`), and its sign orders the
/// printed terms; the native multiprecision evaluator returns an exact zero
/// there. The double evaluator keeps such noise, so its imaginary part is
/// taken in that case (a proxy for the reference noise, see
/// docs/algorithmic-parity.md, ratintex).
fn correctly_rounded(e:&Expr)->Option<Complex64>{
    if let Some(value)=COMPLEX_MEMO.with(|memo|memo.borrow().get(e).copied()){return value;}
    let value=crate::evalf_mp::eval_mp(e,53).map(|value|{
        let (re,im)=(value.re.to_f64(),value.im.to_f64());
        if im==0.0{
            if let Some(noise)=crate::evalf::eval_c64(e).map(|d|d.im).filter(|n|*n!=0.0&&n.is_finite()){return Complex64::new(re,noise);}
        }
        Complex64::new(re,im)
    });
    COMPLEX_MEMO.with(|memo|{
        let mut memo=memo.borrow_mut();
        if memo.len()>=COMPLEX_MEMO_BOUND{memo.clear();}
        memo.insert(e.clone(),value);
    });
    value
}

/// Compare expressions using SymPy default ordering.
pub fn cmp_default(a: &Expr, b: &Expr) -> Ordering {
    let mut cache = SortKeyCache::default();
    cache.key(a).cmp(&cache.key(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_and_piecewise_keys_retain_constructor_pairs() {
        let parse=|source|crate::parser::parse_str(source).unwrap();
        let x=Expr::symbol("x");
        let order=parse("Order(x**2,x)");
        let expected=composite_key((5,0,"Order"),vec![sort_key(&x.powi(2)),
            composite_key((5,0,"Tuple"),vec![sort_key(&x),sort_key(&Expr::zero())])]);
        assert_eq!(class_key(&order),tup(vec![Key::Int(5),Key::Int(0),Key::Str("Order".into())]));
        assert_eq!(sort_key(&order),expected);
        let piecewise=parse("Piecewise((x,x>0),(0,True))");
        let expected=composite_key((4,0,"Piecewise"),vec![
            composite_key((5,0,"ExprCondPair"),vec![sort_key(&x),sort_key(&parse("x>0"))]),
            composite_key((5,0,"ExprCondPair"),vec![sort_key(&Expr::zero()),sort_key(&Expr::boolean(true))])]);
        assert_eq!(sort_key(&piecewise),expected);
    }

    #[test]
    fn dummy_base_is_wrapped_only_for_composite_expression_keys() {
        let d=Expr::dummy_with_index("x",crate::Assumptions::NONE,37.into());
        let atom=sort_key(&d);
        for expression in [crate::ops::neg(&d),crate::ops::mul2(&Expr::int(2),&d),d.powi(2),crate::ops::mul2(&Expr::int(2),&d.powi(2))] {
            let key=sort_key(&expression);
            let args=key.tuple().unwrap()[1].tuple().unwrap();
            assert_eq!(args[0],Key::Int(1));
            assert_eq!(args[1].tuple().unwrap(),&[atom.clone()]);
        }
        assert_eq!(atom.tuple().unwrap()[1].tuple().unwrap()[0],Key::Int(2));
    }

    #[test]
    fn scalar_keys_retain_shared_children_at_depth_eighty() {
        let mut node=Expr::symbol("x");
        for _ in 0..80 { node=Expr::from_unevaluated_kind(Kind::Func(Func::Undefined("pair".into()),vec![node.clone(),node])); }
        let mut cache=SortKeyCache::default();let first=cache.key(&node);let second=cache.key(&node);
        assert_eq!(first,second);assert!(cache.scalars.len()<250);
        let huge=Number::from_bigint(num_bigint::BigInt::from(1)<<10000);
        assert!(Key::SpecialNumber(Constant::Infinity)>Key::Num(huge.clone()));
        assert!(Key::SpecialNumber(Constant::NegInfinity)<Key::Num(huge.neg()));
    }

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
