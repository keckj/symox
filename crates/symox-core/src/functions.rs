//! Elementary and special functions with SymPy-like automatic evaluation.

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{
    add, could_extract_minus_sign, div, mul, mul2, neg, pow, sub,
};

/// Apply function `f` to `args`, performing automatic simplifications.
pub fn func(f: Func, args: Vec<Expr>) -> Expr {
    use crate::special::trig_eval as te;
    match f {
        Func::Logic(op) => return crate::logic::eval(op, args),
        Func::Set(op) => return crate::sets::eval(op, args),
        Func::Order => return crate::series::order::order_from_args(&args),
        Func::Integral | Func::Sum | Func::Product=>{
            if crate::binders::validate(&f,&args).is_ok(){return crate::binders::held(f,&args[0],crate::binders::limits(&args)).unwrap();}
            return Expr::from_kind(Kind::Func(f,args));
        },
        Func::Predicate(_) | Func::AppliedPredicate | Func::Subs | Func::Limit=>return Expr::from_kind(Kind::Func(f,args)),
        Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) => return Expr::from_kind(Kind::Func(f, args)),
        _ => {}
    }
    if args.iter().any(|a| a.is_constant(Constant::NaN)) && !matches!(f, Func::Undefined(_)) {
        return Expr::nan();
    }
    if args.iter().any(crate::evalf_mp::has_mp) {
        let precision = args.iter().filter_map(crate::evalf_mp::max_precision).max().unwrap();
        if args.iter().all(is_numeric_value) {
            let raw = Expr::from_kind(Kind::Func(f, args));
            return crate::evalf_mp::eval_mp(&raw, precision)
                .map(|value| crate::evalf_mp::to_expr(value, precision)).unwrap_or(raw);
        }
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
    // numeric evaluation when an argument is a float
    // floor/ceiling of a Float is an Integer in SymPy: use the exact rule
    if !matches!(f, Func::Undefined(_) | Func::Max | Func::Min | Func::Floor | Func::Ceiling)
        && !args.is_empty()
        && args.iter().all(is_numeric_value)
        && args.iter().any(has_float)
    {
        if let Some(v) = crate::evalf::eval_func_numeric(&f, &args) {
            return v;
        }
    }
    let r = match &f {
        Func::Exp => te::exp_eval(&args[0]),
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
        Func::Abs => eval_abs(&args[0]),
        Func::Sign => eval_sign(&args[0]),
        Func::Floor | Func::Ceiling => eval_floor_ceil(&f, &args[0]),
        Func::Gamma => crate::special::eval_gamma(&args[0]),
        Func::LogGamma => crate::special::eval_loggamma(&args[0]),
        Func::Erf => crate::special::eval_erf(&args[0]),
        Func::Erfc => crate::special::eval_erfc(&args[0]),
        Func::Re => eval_re_im(&args[0], true),
        Func::Im => eval_re_im(&args[0], false),
        Func::Conjugate => eval_conjugate(&args[0]),
        Func::Arg => match args[0].kind() {
            Kind::Num(n) if n.is_positive() => Some(Expr::zero()),
            Kind::Num(n) if n.is_negative() => Some(Expr::pi()),
            Kind::Const(Constant::I) => Some(div(&Expr::pi(), &Expr::int(2))),
            _ => None,
        },
        Func::Heaviside => match args[0].kind() {
            Kind::Num(n) if n.is_positive() => Some(Expr::one()),
            Kind::Num(n) if n.is_negative() => Some(Expr::zero()),
            Kind::Num(_) => Some(Expr::rational(1, 2)),
            _ => None,
        },
        Func::DiracDelta => match args[0].kind() {
            Kind::Num(n) if !n.is_zero() => Some(Expr::zero()),
            _ => None,
        },
        Func::Max | Func::Min => return eval_minmax(&f, args),
        Func::Special(_) | Func::Predicate(_) | Func::AppliedPredicate | Func::Order | Func::Undefined(_) | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::Logic(_) | Func::Set(_)
        | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) => None,
    };
    match r {
        Some(v) => v,
        None => Expr::from_kind(Kind::Func(f, args)),
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
        Kind::Func(f @ (Func::Integral|Func::Sum|Func::Product),args)=>crate::binders::doit(f,args,deep),
        Kind::Derivative(inner, variables) => crate::diff::diff_many(&child(inner), variables),
        Kind::Func(Func::Subs,args)=>crate::subs::subs_doit(&doit(&args[0],deep),&args[1].args(),&args[2].args()),
        Kind::Func(Func::Limit,args)=>{
            let body=child(&args[0]);
            let direction=crate::limits::Dir::parse(args[3].as_sym().map(|s|s.name()).unwrap_or("+"));
            direction.and_then(|dir|crate::limits::limit(&body,&args[1],&args[2],dir).ok())
                .unwrap_or_else(||Expr::from_kind(Kind::Func(Func::Limit,vec![body,args[1].clone(),args[2].clone(),args[3].clone()])))
        },
        Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") => child(&args[0]),
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
        Kind::Num(n) => {
            let fl = match n {
                Number::Float(x) => {
                    let v = if ceil { x.ceil() } else { x.floor() };
                    if v.abs() < 9.0e15 {
                        return Some(Expr::int(v as i64));
                    }
                    return Some(Expr::float(v));
                }
                _ => n.floor(),
            };
            if ceil && fl != *n {
                Some(Expr::num(fl.add(&Number::ONE)))
            } else {
                Some(Expr::num(fl))
            }
        }
        _ => {
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
    // fold numeric values
    let mut best: Option<(f64, Expr)> = None;
    let mut rest: Vec<Expr> = Vec::new();
    for a in flat {
        let v = if is_numeric_value(&a) && !crate::functions::contains_i(&a) {
            crate::evalf::eval_f64(&a)
        } else {
            None
        };
        match v {
            Some(v) => {
                let better = match &best {
                    None => true,
                    Some((b, _)) => {
                        if is_max {
                            v > *b
                        } else {
                            v < *b
                        }
                    }
                };
                if better {
                    best = Some((v, a));
                }
            }
            None => rest.push(a),
        }
    }
    if let Some((_, b)) = best {
        rest.push(b);
    }
    rest.sort();
    rest.dedup();
    // drop arguments dominated by another one (known sign of the difference)
    let mut i = 0;
    while i < rest.len() {
        let dominated = (0..rest.len()).any(|j| {
            j != i && {
                let d = sub(&rest[j], &rest[i]);
                if is_max {
                    crate::assumptions::is_positive(&d) == Some(true)
                } else {
                    crate::assumptions::is_negative(&d) == Some(true)
                }
            }
        });
        if dominated {
            rest.remove(i);
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
