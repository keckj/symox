//! Elementary and special functions with SymPy-like automatic evaluation.

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{
    add, could_extract_minus_sign, div, mul, mul2, neg, pow, sub,
};

/// Apply function `f` to `args`, performing automatic simplifications.
pub fn func(f: Func, args: Vec<Expr>) -> Expr {
    use crate::special::trig_eval as te;
    // A host-defined head evaluates only through its own class, as SymPy's
    // Function.__new__ calls cls.eval for each rebuilt application.
    if let Func::Undefined(head) = &f {
        // A held two-argument log(x, b, evaluate=False) evaluates again when
        // it is rebuilt (subs, doit, expand, func(*args)): log.eval(x, b).
        if args.len() == 2 && is_held_log_base(head) {
            return log_base(&args[0], &args[1]);
        }
        if let Some(host) = head.host() {
            if let Some(value) = host.evaluate(&args) { return value; }
            // A new application gets its own handle, which binds the Python
            // instance of this node only.
            if let Some(handle) = host.applied() {
                return Expr::from_kind(Kind::Func(Func::Undefined(head.clone().with_host(handle)), args));
            }
            return Expr::from_kind(Kind::Func(f, args));
        }
    }
    match f {
        Func::Logic(op) => return crate::logic::eval(op, args),
        Func::Set(op) => return crate::sets::eval(op, args),
        Func::Order => return crate::series::order::order_from_args(&args),
        Func::Integral | Func::Sum | Func::Product=>{
            if crate::binders::validate(&f,&args).is_ok(){return crate::binders::held(f,&args[0],crate::binders::limits(&args)).unwrap();}
            return Expr::from_kind(Kind::Func(f,args));
        },
        Func::GenericExpr | Func::Predicate(_) | Func::AppliedPredicate | Func::Subs | Func::Limit=>return Expr::from_kind(Kind::Func(f,args)),
        Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) => return Expr::from_kind(Kind::Func(f, args)),
        _ => {}
    }
    // Complex log(0) is zoo for exact and floating zeros. Handle it before
    // numerical evaluation, which cannot represent complex infinity.
    if matches!(f,Func::Log)&&args.len()==1&&args[0].as_number().is_some_and(Number::is_zero){return Expr::zoo();}
    if args.iter().any(|a| a.is_constant(Constant::NaN)) && !matches!(f, Func::Undefined(_) | Func::Heaviside) {
        return Expr::nan();
    }
    if matches!(f, Func::Sign) && args.len()==1 && crate::assumptions::is_real(&args[0])==Some(true) {
        if let Some(value)=crate::complex_parts::evaluate::sign(&args[0]) { return value; }
    }
    if let Func::Special(s) = f {
        return crate::special::apply(s, args);
    }
    // exact rules SymPy applies before evaluating a float argument
    if args.len() == 1 && has_float(&args[0]) {
        if let Some(v) = te::float_exact(&f, &args[0]) {
            return v;
        }
    }
    // gamma, loggamma, erf and erfc have no float rule in their eval: like
    // SymPy's Function they evaluate only float or float complex arguments,
    // at the largest argument precision.
    let auto_evalf = matches!(f, Func::Gamma | Func::LogGamma | Func::Erf | Func::Erfc);
    if auto_evalf && crate::special::should_evalf(&args) {
        let prec = args.iter().filter_map(crate::special::float_precision).max().unwrap_or(53);
        let node = Expr::from_kind(Kind::Func(f.clone(), args.clone()));
        if let Some(value) = crate::evalf_mp::eval_number(&node, prec) {
            return value;
        }
        if crate::evalf_mp::has_special_host() || prec > 53 {
            return node;
        }
    }
    let r = match &f {
        Func::Exp => te::exp_eval(&args[0]),
        Func::ExpPolar | Func::GenericExpr => None,
        Func::Log => te::log_eval(&args[0]),
        Func::Sin => te::sin_eval(&args[0]),
        Func::Cos => te::cos_eval(&args[0]),
        Func::Tan => te::tan_eval(&args[0]),
        Func::Cot => te::cot_eval(&args[0]),
        Func::Sec => te::sec_eval(&args[0]),
        Func::Csc => te::csc_eval(&args[0]),
        Func::Asin => te::asin_eval(&args[0]),
        Func::Acos => te::acos_eval(&args[0]),
        Func::Atan => te::atan_eval(&args[0]),
        Func::Acot => te::acot_eval(&args[0]),
        Func::Atan2 => te::atan2_eval(&args[0], &args[1]),
        Func::Sinh => te::sinh_eval(&args[0]),
        Func::Cosh => te::cosh_eval(&args[0]),
        Func::Tanh => te::tanh_eval(&args[0]),
        Func::Coth => te::coth_eval(&args[0]),
        Func::Asinh => te::asinh_eval(&args[0]),
        Func::Acosh => te::acosh_eval(&args[0]),
        Func::Atanh => te::atanh_eval(&args[0]),
        Func::Abs if args[0].has_matrix() => eval_abs(&args[0]),
        Func::Abs => crate::complex_parts::evaluate::abs(&args[0]),
        Func::Sign if args[0].has_matrix() => eval_sign(&args[0]),
        Func::Sign => crate::complex_parts::evaluate::sign(&args[0]),
        Func::Floor | Func::Ceiling => eval_floor_ceil(&f, &args[0]),
        Func::Gamma => crate::special::eval_gamma(&args[0]),
        Func::LogGamma => crate::special::eval_loggamma(&args[0]),
        Func::Erf => crate::special::eval_erf(&args[0]),
        Func::Erfc => crate::special::eval_erfc(&args[0]),
        Func::Re | Func::Im if args[0].has_matrix() => eval_re_im(&args[0], matches!(f, Func::Re)),
        Func::Re => crate::complex_parts::evaluate::re_im(&args[0], true),
        Func::Im => crate::complex_parts::evaluate::re_im(&args[0], false),
        Func::Conjugate if args[0].has_matrix() => eval_conjugate(&args[0]),
        Func::Conjugate => crate::complex_parts::evaluate::conjugate(&args[0]),
        Func::Arg if !args[0].has_matrix() => crate::complex_parts::evaluate::arg(&args[0]),
        Func::Arg => match args[0].kind() {
            Kind::Num(n) if n.is_positive() => Some(Expr::zero()),
            Kind::Num(n) if n.is_negative() => Some(Expr::pi()),
            Kind::Const(Constant::I) => Some(div(&Expr::pi(), &Expr::int(2))),
            _ => None,
        },
        // SymPy's default H0 = 1/2 is the one-argument node.
        Func::Heaviside if args.len() == 2 && args[1] == Expr::rational(1, 2) => {
            Some(func(Func::Heaviside, vec![args[0].clone()]))
        }
        Func::Heaviside if args.len() == 2 && !args[0].is_number()
            && crate::assumptions::is_zero(&args[0]) == Some(true) => Some(args[1].clone()),
        Func::Heaviside => match args[0].kind() {
            Kind::Const(Constant::Infinity) => Some(Expr::one()),
            Kind::Const(Constant::NegInfinity) => Some(Expr::zero()),
            Kind::Const(Constant::NaN) => Some(Expr::nan()),
            Kind::Num(n) if n.is_positive() => Some(Expr::one()),
            Kind::Num(n) if n.is_negative() => Some(Expr::zero()),
            Kind::Num(_) => Some(args.get(1).cloned().unwrap_or_else(|| Expr::rational(1, 2))),
            _ if crate::assumptions::fact(&args[0], "extended_positive") == Some(true) => Some(Expr::one()),
            _ if crate::assumptions::fact(&args[0], "extended_negative") == Some(true) => Some(Expr::zero()),
            _ => None,
        },
        Func::DiracDelta => match args[0].kind() {
            Kind::Num(n) if !n.is_zero() => Some(Expr::zero()),
            _ if crate::assumptions::is_real(&args[0]) == Some(true)
                && crate::assumptions::is_nonzero(&args[0]) == Some(true) => Some(Expr::zero()),
            _ if args.len() == 2 && args[1].is_zero() => Some(func(Func::DiracDelta, vec![args[0].clone()])),
            Kind::Mul(coefficient, _) if coefficient.is_negative()
                && (args.len() == 1 || args[1].as_number().is_some_and(|n| n.is_integer() && !n.is_negative())) => {
                let mut positive = args.to_vec();
                positive[0] = neg(&positive[0]);
                let value = func(Func::DiracDelta, positive);
                Some(if args.len() == 2 && crate::assumptions::is_odd(&args[1]) == Some(true) {
                    neg(&value)
                } else { value })
            },
            _ => None,
        },
        Func::Max | Func::Min => return eval_minmax(&f, args),
        Func::Special(_) | Func::Predicate(_) | Func::AppliedPredicate | Func::Order | Func::Undefined(_) | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::Logic(_) | Func::Set(_)
        | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) => None,
    };
    if let Some(v) = r { return v; }
    // Like SymPy's Function.__new__, numerical evaluation happens only when
    // the function's own eval rule returned nothing.
    // SymPy's Function._should_evalf: a call evaluates numerically when every
    // argument is a Float or a sum a + b*I of numbers with a Float part, at
    // the largest argument precision. Other numeric arguments containing a
    // Float (2.5*I, 2.5*pi) go through the function's exact rules instead.
    // The value is computed in multiprecision and rounded once, as mpmath
    // does, never taken from a double-precision evaluation.
    // floor/ceiling of a Float is an Integer in SymPy: use the exact rule
    // arg of a real float is exact (0, pi or nan) as in SymPy
    if !matches!(f, Func::Undefined(_) | Func::ExpPolar | Func::Max | Func::Min | Func::Floor | Func::Ceiling | Func::Heaviside | Func::DiracDelta)
        && !(matches!(f, Func::Arg) && args.first().is_some_and(Expr::is_number))
        && crate::special::should_evalf(&args)
    {
        let precision = args.iter().filter_map(crate::special::float_precision).max().unwrap_or(53);
        let raw = Expr::from_kind(Kind::Func(f.clone(), args.clone()));
        if let Some(value) = crate::evalf_mp::eval_mp(&raw, precision) {
            return crate::evalf_mp::to_expr(value, precision);
        }
    }
    Expr::from_kind(Kind::Func(f, args))
}

/// `im(arg).is_zero` for the real-value check of DiracDelta and Heaviside,
/// given `imaginary = im(arg)`.
///
/// SymPy's assumption rules make every noncommutative expression non-real
/// (`extended_real -> commutative`), so `im(X)` and `re(X)` of a
/// noncommutative `X` are nonzero (`im._eval_is_zero` is
/// `X.is_extended_real`, `re._eval_is_zero` is `X.is_imaginary or
/// X.is_zero`). Products of such parts follow `Mul._eval_is_zero`. Other
/// nodes keep the exact native query.
pub fn imaginary_part_is_zero(imaginary: &Expr) -> Option<bool> {
    let noncommutative = |e: &Expr| !crate::decompose::commutes(e);
    match imaginary.kind() {
        Kind::Func(Func::Re | Func::Im, a) if a.len() == 1 && noncommutative(&a[0]) => Some(false),
        Kind::Mul(c, factors) if imaginary.has_noncommutative_host() => {
            // SymPy's Mul._eval_is_zero over the coefficient and the factors.
            let mut items: Vec<Expr> = Vec::with_capacity(factors.len() + 1);
            if !c.is_one() {
                items.push(Expr::num(c.clone()));
            }
            items.extend(factors.iter().map(|(b, x)| crate::ops::mul_factor(b, x)));
            let (mut zero, mut infinite) = (Some(false), false);
            for item in &items {
                let z = imaginary_part_is_zero(item);
                if z == Some(true) {
                    if infinite {
                        return None;
                    }
                    zero = Some(true);
                } else {
                    if crate::assumptions::fact(item, "finite") != Some(true) {
                        if zero == Some(true) {
                            return None;
                        }
                        infinite = true;
                    }
                    if zero == Some(false) && z.is_none() {
                        zero = None;
                    }
                }
            }
            zero
        }
        _ => crate::assumptions::is_zero_exact(imaginary),
    }
}

/// `Add._eval_is_zero` leaves a sum of nonzero terms undecided; SymPy then
/// settles `zero` through `zero -> integer` and `zero -> algebraic` with
/// Add's `_fuzzy_group(..., quick_exit=True)`: `im(A) + 1` with a
/// noncommutative `A` is not an integer, hence nonzero. Some(false) or None,
/// for a sum whose parts hold a noncommutative argument.
pub(crate) fn noncommutative_sum_nonzero(sum: &Expr) -> Option<bool> {
    let terms = sum.args();
    let algebraic: Vec<Option<bool>> = terms.iter().map(noncommutative_part_algebraic).collect();
    let integer = terms.iter().zip(&algebraic).map(|(t, a)| if *a == Some(false) { Some(false) } else { crate::assumptions::fact(t, "integer") });
    if fuzzy_group(algebraic.iter().copied(), true) == Some(false) || fuzzy_group(integer, true) == Some(false) {
        return Some(false);
    }
    None
}

/// SymPy's `_fuzzy_group`: True when every value is True, None as soon as
/// one is None, and False for exactly one False (any number without
/// `quick_exit`).
fn fuzzy_group(values: impl Iterator<Item = Option<bool>>, quick_exit: bool) -> Option<bool> {
    let mut saw_other = false;
    for value in values {
        match value {
            Some(true) => continue,
            None => return None,
            Some(false) => {
                if quick_exit && saw_other {
                    return None;
                }
                saw_other = true;
            }
        }
    }
    Some(!saw_other)
}

/// `is_algebraic` of a part of `im(z)` or `re(z)`: `re`/`im` defer to their
/// argument (False when it does not commute), Mul uses `_fuzzy_group` and
/// Pow `Pow._eval_is_algebraic` for a rational exponent.
fn noncommutative_part_algebraic(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Func(Func::Re | Func::Im, a) if a.len() == 1 => crate::assumptions::fact(&a[0], "algebraic"),
        Kind::Mul(c, factors) if e.has_noncommutative_host() => {
            let coefficient = std::iter::once(if c.is_float() { None } else { Some(true) });
            fuzzy_group(coefficient.chain(factors.iter().map(|(b, x)| noncommutative_part_algebraic(&crate::ops::mul_factor(b, x)))), true)
        }
        Kind::Pow(b, x) if e.has_noncommutative_host() && x.as_number().is_some_and(|n| n.is_rational() && !n.is_zero()) => {
            if crate::assumptions::is_zero_exact(b) == Some(true) || b.is_one() { return Some(true); }
            if noncommutative_part_algebraic(b) == Some(false) { Some(false) } else { None }
        }
        _ => crate::assumptions::fact(e, "algebraic"),
    }
}

/// Evaluate held operations. Matrix scalar nodes retain their typed parents;
/// abstract symbols stay abstract instead of being materialized implicitly.
pub fn doit(e: &Expr, deep: bool) -> Expr {
    let child = |value: &Expr| if deep { doit(value, true) } else { value.clone() };
    match e.kind() {
        Kind::Func(Func::MatrixDotProduct(a,b), _) => {
            let parents=if deep {a.doit(true).and_then(|a|b.doit(true).map(|b|(a,b)))}else{Ok(((**a).clone(),(**b).clone()))};
            parents.and_then(|(a,b)|crate::matrixcontraction::dot_product(&a,&b,true)).unwrap_or_else(|_|e.clone())
        }
        Kind::Func(Func::MatrixTrace(parent), _) => {
            let parent = if deep { parent.doit(true).ok() } else { Some((**parent).clone()) };
            parent.and_then(|p| p.trace().ok()).unwrap_or_else(|| e.clone())
        }
        Kind::Func(Func::MatrixDeterminant(parent), _) => {
            let parent = if deep { parent.doit(true).ok() } else { Some((**parent).clone()) };
            parent.and_then(|p| p.determinant().ok()).unwrap_or_else(|| e.clone())
        }
        Kind::Func(Func::MatrixElement(parent), args) if args.len() == 2 => {
            // SymPy's MatrixElement.doit: args[0][args[1], args[2]] (the
            // constructor keeps product entries held, `_entry` expands them).
            let parent = if deep { parent.doit(true).ok() } else { Some((**parent).clone()) };
            parent.and_then(|p| p.element(&child(&args[0]), &child(&args[1])).ok()).unwrap_or_else(|| e.clone())
        }
        Kind::Func(f @ (Func::Integral|Func::Sum|Func::Product),args)=>crate::binders::doit(f,args,deep),
        Kind::Func(Func::Special(s), args) if s.is_transform() => {
            let args: Vec<Expr> = args.iter().map(child).collect();
            crate::transforms::doit_default(*s, &args).unwrap_or_else(|| crate::special::apply(*s, args))
        }
        Kind::Derivative(inner, variables) => {
            // Derivative.doit: a changed result that still holds a derivative
            // is done again, so the integrals it creates are done too.
            let rv = crate::diff::diff_evaluated(&child(inner), variables);
            let mut held = false;
            crate::traverse::preorder(&rv, &mut |s| if matches!(s.kind(), Kind::Derivative(..)) { held = true; });
            if rv != *e && held { doit(&rv, deep) } else { rv }
        }
        Kind::Func(Func::Subs,args)=>crate::subs::subs_doit(&doit(&args[0],deep),&args[1].args(),&args[2].args()),
        Kind::Func(Func::Limit,args)=>{
            let body=child(&args[0]);
            let direction=crate::limits::Dir::parse(args[3].as_sym().map(|s|s.name()).unwrap_or("+"));
            direction.and_then(|dir|crate::limits::limit(&body,&args[1],&args[2],dir).ok())
                .unwrap_or_else(||Expr::from_kind(Kind::Func(Func::Limit,vec![body,args[1].clone(),args[2].clone(),args[3].clone()])))
        },
        Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") => child(&args[0]),
        Kind::Func(Func::Undefined(head), args) if args.len() == 2 && is_held_log_base(head) => log_base(&child(&args[0]), &child(&args[1])),
        Kind::Add(..) if !e.is_canonical() => add(&e.args().iter().map(child).collect::<Vec<_>>()),
        Kind::Mul(..) if !e.is_canonical() => mul(&e.args().iter().map(child).collect::<Vec<_>>()),
        Kind::Pow(base, exponent) if !e.is_canonical() => pow(&child(base), &child(exponent)),
        _ if deep => crate::traverse::map_children(e, &mut |value| doit(value, true)),
        _ => e.clone(),
    }
}

/// Lower trace and determinant nodes to finite scalar expressions for code
/// generation. This is explicit, so ordinary simplify/doit preserve abstract
/// matrix symbols as SymPy does.
pub fn materialize_matrix_scalars(e: &Expr) -> crate::matrix::MResult<Expr> {
    let mut failure = None;
    let mapped = crate::traverse::map_children(e, &mut |value| {
        match materialize_matrix_scalars(value) {
            Ok(value) => value,
            Err(error) => { failure = Some(error); value.clone() }
        }
    });
    if let Some(error) = failure { return Err(error); }
    match mapped.kind() {
        Kind::Func(Func::MatrixDotProduct(a,b), _) => crate::matrixcontraction::dot_product(a,b,true),
        Kind::Func(Func::MatrixTrace(parent), _) => parent.trace_explicit(),
        Kind::Func(Func::MatrixDeterminant(parent), _) => parent.determinant_explicit(),
        _ => Ok(mapped),
    }
}

fn has_float(e: &Expr) -> bool {
    e.float_info().max_precision.is_some()
}

/// Free of symbols and undefined functions (so it can be evaluated).
pub fn is_numeric_value(e:&Expr)->bool{
    // Classification is independent of numeric-refinement policy and uses the
    // two bits above the two assumption-cache banks (bits 60 and 61).
    if e.mask()!=0 || e.is_atom(){return is_numeric_value_uncached(e);}
    match e.cached_assumption(60,2){2=>return false,3=>return true,_=>{}}
    let value=is_numeric_value_uncached(e);
    e.remember_assumption(60,2,if value{3}else{2});value
}
fn is_numeric_value_uncached(e: &Expr) -> bool {
    // any symbol (or Order term) inside sets a mask bit: O(1) rejection
    if e.mask() != 0 {
        return false;
    }
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => !c.is_infinite() && *c != Constant::NaN,
        Kind::Sym(_) | Kind::Derivative(..) | Kind::Rel(..) | Kind::Piecewise(_) | Kind::Bool(_) => {
            false
        }
        Kind::Func(Func::Logic(_) | Func::Set(_), _) => false,
        Kind::Func(Func::Undefined(n), a) if crate::poly::rootof::is_rootof_head(n, a.len()) => true,
        Kind::Func(Func::Undefined(n), a) => n.is_internal("Tuple") && a.iter().all(is_numeric_value),
        // Inspect compact payloads instead of args(): reconstructing public
        // terms and powers may canonicalize a shared held arithmetic DAG.
        Kind::Add(_, terms) => terms.iter().all(|(term, coefficient)| {
            coefficient.is_zero() || is_numeric_value(term)
        }),
        Kind::Mul(_, factors) => factors.iter().all(|(base, exponent)| {
            if exponent.is_zero() || base.is_one() {
                true
            } else if base.is_atom() && exponent.is_atom() && !exponent.is_one() {
                // Scalar power simplifications, e.g. oo**-1 -> 0, are bounded.
                is_numeric_value(&crate::ops::pow(base, exponent))
            } else {
                is_numeric_value(base) && is_numeric_value(exponent)
            }
        }),
        Kind::Pow(base, exponent) => is_numeric_value(base) && is_numeric_value(exponent),
        _ => e.args().iter().all(is_numeric_value),
    }
}

// --------------------------------------------------------------------
// constructors
// --------------------------------------------------------------------

macro_rules! unary {
    ($($name:ident => $f:expr),* $(,)?) => {
        $(
            pub fn $name(x: &Expr) -> Expr { func($f, vec![x.clone()]) }
        )*
    };
}

unary! {
    exp => Func::Exp, log => Func::Log,
    sin => Func::Sin, cos => Func::Cos, tan => Func::Tan,
    cot => Func::Cot, sec => Func::Sec, csc => Func::Csc,
    asin => Func::Asin, acos => Func::Acos, atan => Func::Atan, acot => Func::Acot,
    sinh => Func::Sinh, cosh => Func::Cosh, tanh => Func::Tanh, coth => Func::Coth,
    asinh => Func::Asinh, acosh => Func::Acosh, atanh => Func::Atanh,
    abs => Func::Abs, sign => Func::Sign, floor => Func::Floor, ceiling => Func::Ceiling,
    gamma => Func::Gamma, loggamma => Func::LogGamma, erf => Func::Erf, erfc => Func::Erfc,
    re => Func::Re, im => Func::Im, conjugate => Func::Conjugate, arg => Func::Arg,
    heaviside => Func::Heaviside, dirac_delta => Func::DiracDelta,
}

pub fn atan2(y: &Expr, x: &Expr) -> Expr {
    func(Func::Atan2, vec![y.clone(), x.clone()])
}

pub fn max(args: Vec<Expr>) -> Expr {
    func(Func::Max, args)
}

pub fn min(args: Vec<Expr>) -> Expr {
    func(Func::Min, args)
}

/// Logarithm in base `b` (SymPy's `log(x, b)`, with the integer part of
/// the result split off for rational arguments).
pub fn log_base(x: &Expr, b: &Expr) -> Expr {
    crate::special::trig_eval::log_base(x, b)
}

/// The head of SymPy's held `log(x, b, evaluate=False)`: the native log is
/// unary, so the held two-argument application keeps a `log` head with
/// `nargs=(1, 2)`, which evaluates again through `log.eval(x, b)` when it is
/// rebuilt.
pub fn held_log_base_head() -> crate::expr::UndefinedFunction {
    crate::expr::UndefinedFunction::new("log", crate::symbol::Assumptions::NONE, Some(vec![1, 2]))
}

pub fn is_held_log_base(head: &crate::expr::UndefinedFunction) -> bool {
    !head.internal && head.host().is_none() && head.name.as_ref() == "log" && head.nargs.as_deref() == Some(&[1, 2][..])
}

/// Applied undefined function `name(args...)`.
pub fn undefined(name: &str, args: Vec<Expr>) -> Expr {
    func(Func::Undefined(name.into()), args)
}

#[cfg(test)]
#[path="functions_api_tests.rs"]
mod functions_api_tests;

// --------------------------------------------------------------------
// helpers
// --------------------------------------------------------------------

/// `e / I` if `e` is syntactically `I * something`.
fn extract_i(e: &Expr) -> Option<Expr> {
    match e.kind() {
        Kind::Const(Constant::I) => Some(Expr::one()),
        Kind::Mul(c, facs) => {
            let idx = facs
                .iter()
                .position(|(b, x)| b.is_constant(Constant::I) && x.is_one())?;
            let mut fs = facs.clone();
            fs.remove(idx);
            let mut v: Vec<Expr> = fs.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).collect();
            v.push(Expr::num(c.clone()));
            Some(mul(&v))
        }
        _ => None,
    }
}

// --------------------------------------------------------------------
// misc
// --------------------------------------------------------------------

fn eval_abs(a: &Expr) -> Option<Expr> {
    if !a.is_number() {
        if crate::assumptions::is_nonnegative(a) == Some(true) {
            return Some(a.clone());
        }
        if crate::assumptions::is_nonpositive(a) == Some(true) {
            return Some(neg(a));
        }
    }
    match a.kind() {
        Kind::Num(n) => Some(Expr::num(n.abs())),
        Kind::Const(Constant::I) => Some(Expr::one()),
        Kind::Const(c) if c.is_number_symbol() => Some(a.clone()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity) => {
            Some(Expr::oo())
        }
        Kind::Func(Func::Abs, _) => Some(a.clone()),
        Kind::Func(Func::Exp, args) => Some(exp(&re(&args[0]))),
        Kind::Mul(c, facs) => {
            // Abs(c * x) = |c| * Abs(x),  Abs(I*x) = Abs(x)
            let mut rest: Vec<Expr> = Vec::new();
            let mut changed = !c.is_one();
            for (position, (b, e)) in facs.iter().enumerate() {
                if b.is_constant(Constant::I) && e.is_one() {
                    changed = true;
                    continue;
                }
                if matches!(b.kind(), Kind::Const(Constant::Pi | Constant::E))
                    || (b.is_number() && b.as_number().unwrap().is_positive() && e.is_rational_number())
                    || crate::assumptions::is_positive(&crate::ops::mul_factor(b, e)) == Some(true)
                {
                    // positive factor
                    return Some(mul2(
                        &crate::ops::mul_factor(b, e),
                        &abs(&div(a, &crate::ops::mul_factor(b, e))),
                    ));
                }
                // SymPy's Abs.eval splits off factors whose Abs evaluates,
                // e.g. Abs((-1)**n*x) = Abs(x) for an integer n
                let factor = crate::ops::mul_factor(b, e);
                if matches!(factor.kind(), Kind::Pow(..)) {
                    if let Some(known) = eval_abs_pow(&factor) {
                        // rebuild the cofactor: dividing would leave
                        // (-1)**(n + 1)*(-1)**(-n - 1) uncancelled
                        let mut cofactor = vec![Expr::num(c.clone())];
                        cofactor.extend(facs.iter().enumerate().filter(|(j, _)| *j != position)
                            .map(|(_, (b, e))| crate::ops::mul_factor(b, e)));
                        return Some(mul2(&known, &abs(&mul(&cofactor))));
                    }
                }
                rest.push(factor);
            }
            if changed {
                Some(mul2(&Expr::num(c.abs()), &abs(&mul(&rest))))
            } else {
                None
            }
        }
        Kind::Pow(..) => eval_abs_pow(a),
        _ => {
            if could_extract_minus_sign(a) {
                Some(abs(&neg(a)))
            } else {
                None
            }
        }
    }
}

/// SymPy's `Abs.eval` for `base**n` with a real base and an integer `n`:
/// even powers are their own absolute value, `Abs((-1)**n)` is one and
/// otherwise `Abs(base)**n`.
fn eval_abs_pow(a: &Expr) -> Option<Expr> {
    let Kind::Pow(base, exponent) = a.kind() else { return None };
    if exponent.is_number() || crate::assumptions::is_real(base) != Some(true)
        || crate::assumptions::is_integer(exponent) != Some(true)
    {
        return None;
    }
    if crate::assumptions::is_even(exponent) == Some(true) {
        return Some(a.clone());
    }
    if *base == Expr::minus_one() {
        return Some(Expr::one());
    }
    let magnitude = abs(base);
    (magnitude != *base).then(|| pow(&magnitude, exponent))
}

fn eval_sign(a: &Expr) -> Option<Expr> {
    if !a.is_number() {
        if crate::assumptions::is_positive(a) == Some(true) {
            return Some(Expr::one());
        }
        if crate::assumptions::is_negative(a) == Some(true) {
            return Some(Expr::minus_one());
        }
        if crate::assumptions::is_zero(a) == Some(true) {
            return Some(Expr::zero());
        }
    }
    match a.kind() {
        Kind::Num(n) => Some(Expr::int(if n.is_positive() {
            1
        } else if n.is_negative() {
            -1
        } else {
            0
        })),
        Kind::Const(Constant::Infinity) => Some(Expr::one()),
        Kind::Const(c) if c.is_number_symbol() => Some(Expr::one()),
        Kind::Const(Constant::NegInfinity) => Some(Expr::minus_one()),
        Kind::Const(Constant::I) => Some(Expr::i()),
        _ => {
            if could_extract_minus_sign(a) {
                Some(neg(&sign(&neg(a))))
            } else {
                None
            }
        }
    }
}

fn eval_floor_ceil(f: &Func, a: &Expr) -> Option<Expr> {
    let ceil = matches!(f, Func::Ceiling);
    if !a.is_number() && crate::assumptions::is_integer(a) == Some(true) {
        return Some(a.clone());
    }
    match a.kind() {
        // The integer part of a Float with a huge binary exponent cannot be
        // materialized; SymPy raises OverflowError (the Python entry does
        // too) and the function stays held here.
        Kind::Num(Number::MpFloat(value, _)) if !value.integer_part_fits() => None,
        Kind::Num(n) => {
            let fl = n.floor();
            if ceil && fl.cmp_value(n) != Some(std::cmp::Ordering::Equal) {
                Some(Expr::num(fl.add(&Number::ONE)))
            } else {
                Some(Expr::num(fl))
            }
        }
        _ => {
            if matches!(a.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity)) {
                return Some(a.clone());
            }
            if a.is_add() {
                let (integers, remainder): (Vec<_>, Vec<_>) = a.args().into_iter()
                    .partition(|term| crate::assumptions::is_integer(term) == Some(true));
                if !integers.is_empty() {
                    return Some(add(&[add(&integers), func(f.clone(), vec![add(&remainder)])]));
                }
            }
            // A noncommutative argument is neither imaginary nor `I` times a
            // real value (SymPy's assumption rules), so RoundFunction.eval
            // keeps it whole: `floor(x - I*A*B)`.
            if crate::traverse::has(a, &Expr::i()) && crate::decompose::commutes(a) {
                if let Ok(Some((real, imaginary))) = crate::complex_parts::as_real_imag_bounded(a, false, None) {
                    if !imaginary.is_zero() {
                        return Some(add(&[func(f.clone(), vec![real]), mul2(&Expr::i(), &func(f.clone(), vec![imaginary]))]));
                    }
                }
            }
            if let Some(integer) = crate::special::lambert::integer_round(a, if ceil { 1 } else { -1 }) {
                return Some(Expr::num(Number::from_bigint(integer)));
            }
            if crate::functions::is_numeric_value(a) {
                // evaluate numerically with a safety margin
                let v = crate::evalf::eval_f64(a)?;
                let r = if ceil { v.ceil() } else { v.floor() };
                if (v - r).abs() > 1e-9 && r.abs() < 9.0e15 {
                    return Some(Expr::int(r as i64));
                }
            }
            None
        }
    }
}

/// Certified integer truncation for exact rational, dyadic and supported
/// rationally enclosed real expressions. No f64-to-i64 saturation is involved.
pub fn integer_trunc(a: &Expr) -> Option<num_bigint::BigInt> {
    if let Some(n) = a.as_number() {
        return match n {
            Number::Float(f) if f.is_finite() => Some(crate::bigfloat::BigFloat::from_f64(*f).to_int_trunc()),
            Number::MpFloat(value, _) => value.checked_int_trunc(),
            _ => { let (p,q) = n.numer_denom()?; Some(p/q) },
        };
    }
    crate::special::lambert::integer_round(a, 0)
}

#[cfg(test)]
mod certified_integer_tests {
    use super::*;
    fn expression(source: &str) -> Expr { crate::parser::parse_str(source).unwrap() }

    #[test]
    fn integer_parts_are_not_saturated_at_machine_word_size() {
        let huge = num_bigint::BigInt::from(10).pow(100);
        let x = expression("10**100 + sqrt(2)");
        assert_eq!(integer_trunc(&x), Some(&huge + 1));
        assert_eq!(integer_trunc(&neg(&x)), Some(-(&huge + 1_u32)));
        assert_eq!(func(Func::Floor, vec![x.clone()]), Expr::num(Number::from_bigint(&huge + 1)));
        assert_eq!(func(Func::Ceiling, vec![x]), Expr::num(Number::from_bigint(&huge + 2)));
        let point = Expr::float(1e100);
        assert!(func(Func::Floor, vec![point.clone()]).is_integer_number());
        assert_eq!(integer_trunc(&point), Some(crate::bigfloat::BigFloat::from_f64(1e100).to_int_trunc()));
    }

    #[test]
    fn floor_and_ceiling_prove_very_close_integer_boundaries() {
        let below = expression("1 - 10**(-1000)");
        assert_eq!(integer_trunc(&below), Some(0.into()));
        let huge = num_bigint::BigInt::from(10).pow(100);
        let near = expression("sqrt(10**200 + 1)");
        assert_eq!(func(Func::Floor, vec![near.clone()]), Expr::num(Number::from_bigint(huge.clone())));
        assert_eq!(func(Func::Ceiling, vec![near]), Expr::num(Number::from_bigint(huge + 1)));
        assert_eq!(func(Func::Floor, vec![expression("exp(-10**(-1000))")]), Expr::zero());
        assert_eq!(func(Func::Ceiling, vec![expression("exp(10**(-1000))")]), Expr::int(2));
        assert_eq!(integer_trunc(&Expr::pi()), Some(3.into()));
        assert_eq!(func(Func::Floor, vec![Expr::pi()]), Expr::int(3));
        assert_eq!(func(Func::Ceiling, vec![Expr::pi()]), Expr::int(4));
        assert_eq!(func(Func::Floor, vec![Expr::i()]), Expr::i());
        assert_eq!(func(Func::Floor, vec![expression("3/2 + 7*I/3")]), expression("1 + 2*I"));
    }
}

/// Split a numeric value into exact real and imaginary parts when trivially
/// possible.
fn eval_re_im(a: &Expr, want_re: bool) -> Option<Expr> {
    if a.as_sym().is_some_and(|s|s.assumptions().has(crate::Assumptions::IMAGINARY)){
        return Some(if want_re{Expr::zero()}else{mul2(&neg(&Expr::i()),a)});
    }
    if !a.is_number() && crate::assumptions::is_real(a) == Some(true) {
        return Some(if want_re { a.clone() } else { Expr::zero() });
    }
    match a.kind() {
        Kind::Num(_) => {
            Some(if want_re { a.clone() } else { Expr::zero() })
        }
        Kind::Const(c) if c.is_number_symbol() => {
            Some(if want_re { a.clone() } else { Expr::zero() })
        }
        Kind::Const(Constant::I) => Some(if want_re { Expr::zero() } else { Expr::one() }),
        Kind::Add(c, terms) => {
            let mut v = vec![if want_re { Expr::num(c.clone()) } else { Expr::zero() }];
            for (t, k) in terms {
                let tt = crate::ops::add_term(t, k);
                v.push(func(if want_re { Func::Re } else { Func::Im }, vec![tt]));
            }
            Some(add(&v))
        }
        Kind::Mul(c, _) if c.is_float() || !c.is_one() => {
            let rest = div(a, &Expr::num(c.clone()));
            let inner = func(if want_re { Func::Re } else { Func::Im }, vec![rest]);
            Some(mul2(&Expr::num(c.clone()), &inner))
        }
        Kind::Mul(_, _) => {
            if let Some(r) = extract_i(a) {
                // re(I*x) = -im(x), im(I*x) = re(x)
                return Some(if want_re { neg(&im(&r)) } else { re(&r) });
            }
            None
        }
        Kind::Func(Func::Re | Func::Im | Func::Abs, _) => {
            Some(if want_re { a.clone() } else { Expr::zero() })
        }
        _ => None,
    }
}

fn eval_conjugate(a: &Expr) -> Option<Expr> {
    if a.as_sym().is_some_and(|s|s.assumptions().has(crate::Assumptions::IMAGINARY)){
        return Some(neg(a));
    }
    if !a.is_number() && crate::assumptions::is_real(a) == Some(true) {
        return Some(a.clone());
    }
    match a.kind() {
        Kind::Num(_) => Some(a.clone()),
        Kind::Const(c) if c.is_number_symbol() => Some(a.clone()),
        Kind::Const(Constant::I) => Some(neg(a)),
        Kind::Func(Func::Conjugate, x) => Some(x[0].clone()),
        Kind::Pow(base, exponent) if crate::assumptions::is_integer(exponent) == Some(true) => {
            Some(pow(&conjugate(base), exponent))
        }
        Kind::Add(..) | Kind::Mul(..) => {
            let args: Vec<Expr> = a.args().iter().map(conjugate).collect();
            Some(if a.is_add() { add(&args) } else { mul(&args) })
        }
        _ => None,
    }
}

pub(crate) fn contains_i(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::I) => true,
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(contains_i),
    }
}

fn eval_minmax(f: &Func, args: Vec<Expr>) -> Expr {
    let is_max = matches!(f, Func::Max);
    let mut flat: Vec<Expr> = Vec::new();
    for a in args {
        match a.kind() {
            Kind::Func(g, inner) if g == f => flat.extend(inner.iter().cloned()),
            _ => flat.push(a),
        }
    }
    // the identity element drops out and the absorbing element wins
    let (identity, absorbing) = if is_max { (Expr::neg_oo(), Expr::oo()) } else { (Expr::oo(), Expr::neg_oo()) };
    if flat.contains(&absorbing) {
        return absorbing;
    }
    flat.retain(|a| *a != identity);
    // Compare numeric atoms at their stored exact values. Converting them
    // to f64 loses tiny rationals and adjacent large integers.
    let mut best: Option<Expr> = None;
    let mut rest = Vec::new();
    for a in flat {
        if let Some(value) = a.as_number() {
            let replace = best.as_ref().is_none_or(|previous| {
                let old = previous.as_number().unwrap();
                match value.cmp_value(old) {
                    Some(std::cmp::Ordering::Greater) => is_max,
                    Some(std::cmp::Ordering::Less) => !is_max,
                    Some(std::cmp::Ordering::Equal) => {
                        let new_float = value.float_precision();
                        let old_float = old.float_precision();
                        if new_float.is_some() != old_float.is_some() {
                            if is_max { new_float.is_none() } else { new_float.is_some() }
                        } else { new_float > old_float }
                    }
                    None => false,
                }
            });
            if replace { best = Some(a); }
        } else { rest.push(a); }
    }
    if let Some(best) = best { rest.push(best); }
    rest.sort();
    rest.dedup();
    let dual=if is_max {Func::Min}else{Func::Max};
    let mut members:Vec<Option<rustc_hash::FxHashSet<Expr>>>=rest.iter().map(|value| {
        match value.kind() {
            Kind::Func(head,args) if head==&dual=>Some(args.iter().cloned().collect()),
            _=>None,
        }
    }).collect();
    // drop arguments dominated by another one (known sign of the difference)
    let mut i = 0;
    while i < rest.len() {
        let dominated = (0..rest.len()).any(|j| {
            j != i && {
                let a = &rest[j];
                let b = &rest[i];
                if let Some(large)=&members[i] {
                    if large.contains(a) || members[j].as_ref().is_some_and(|small|small.is_subset(large)) {
                        return true;
                    }
                }
                let difference = sub(a,b);
                let (exact,equal) = crate::assumptions::exact_queries(|| {
                    let forward=if is_max { crate::assumptions::is_nonnegative(&difference) }
                        else { crate::assumptions::is_nonpositive(&difference) };
                    let reverse=if forward==Some(true) {
                        if is_max {crate::assumptions::is_nonpositive(&difference)}
                        else {crate::assumptions::is_nonnegative(&difference)}
                    } else {None};
                    (forward,reverse==Some(true))
                });
                if exact == Some(true) { !equal || j<i }
                else if exact == Some(false) { false }
                else if is_numeric_value(a) && is_numeric_value(b) && !contains_i(a) && !contains_i(b) {
                    // Keep unresolved near-ties symbolic rather than discard
                    // an argument because both approximations round alike.
                    match (crate::evalf::eval_f64(a),crate::evalf::eval_f64(b)) {
                        (Some(x),Some(y)) if x.is_finite() && y.is_finite() => {
                            let gap = x-y;
                            let separation = 1e-12*x.abs().max(y.abs()).max(1.0);
                            if is_max { gap > separation } else { gap < -separation }
                        }
                        // Outside the double range the reference still
                        // compares mpmath values (Max(3, exp(710)) is exp(710)).
                        _ => match (crate::evalf_mp::eval_mp_default(a,96),crate::evalf_mp::eval_mp_default(b,96)) {
                            (Some(x),Some(y)) if x.is_real() && y.is_real() => {
                                let gap=if is_max {x.re.sub(&y.re,96)} else {y.re.sub(&x.re,96)};
                                let scale=if x.re.abs().cmp(&y.re.abs()).is_ge() {x.re.abs()} else {y.re.abs()};
                                let scale=if scale.cmp(&crate::bigfloat::BigFloat::one()).is_ge() {scale} else {crate::bigfloat::BigFloat::one()};
                                gap.sign()>0 && gap.cmp(&scale.ldexp(-40)).is_gt()
                            }
                            _ => false,
                        },
                    }
                } else { false }
            }
        });
        if dominated {
            rest.remove(i);
            members.remove(i);
        } else {
            i += 1;
        }
    }
    if rest.len() == 1 {
        return rest.pop().unwrap();
    }
    if rest.is_empty() {
        return if is_max { Expr::neg_oo() } else { Expr::oo() };
    }
    // SymPy: `Expr.__new__(cls, *ordered(args))`
    crate::logic::ordered(&mut rest);
    Expr::from_kind(Kind::Func(f.clone(), rest))
}
