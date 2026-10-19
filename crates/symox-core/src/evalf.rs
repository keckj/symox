//! Numerical evaluation (double precision, real and complex).

use num_complex::Complex64 as C64;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{func, is_numeric_value};
use crate::number::Number;
use crate::ops::{add, mul, pow};

const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

// Memoize by physical identity, retaining owners of temporary nodes from doit.
// A root numeric evaluation owns the scope; nested special-function and power
// evaluations reuse it. Environment-dependent eval_with never uses this table.
struct NumericMemo {
    depth:usize,
    values:Option<rustc_hash::FxHashMap<(usize,bool),(Expr,Option<C64>)>>,
}
thread_local! {
    static NUMERIC_MEMO:std::cell::RefCell<NumericMemo>=const{std::cell::RefCell::new(NumericMemo{depth:0,values:None})};
    #[cfg(test)]
    static EVALUATIONS:std::cell::Cell<usize>=const{std::cell::Cell::new(0)};
}
struct NumericScope;
impl NumericScope {
    fn enter()->Self{NUMERIC_MEMO.with(|memo|memo.borrow_mut().depth+=1);Self}
}
impl Drop for NumericScope {
    fn drop(&mut self){NUMERIC_MEMO.with(|memo|{
        let mut memo=memo.borrow_mut();memo.depth-=1;if memo.depth==0{memo.values=None;}
    });}
}
pub(crate) fn with_numeric_cache<R>(evaluation:impl FnOnce()->R)->R{
    let _scope=NumericScope::enter();evaluation()
}


/// Evaluate a numeric expression (no free symbols) to a complex double.
pub fn eval_c64(e:&Expr)->Option<C64>{
    if e.is_atom(){return eval_c64_uncached(e);}
    let _scope=NumericScope::enter();
    let key=(e.dag_id(),crate::assumptions::numeric_refinement_enabled());
    let cached=NUMERIC_MEMO.with(|memo|memo.borrow().values.as_ref().and_then(|values|values.get(&key).map(|(_,value)|*value)));
    if let Some(value)=cached{return value;}
    #[cfg(test)] EVALUATIONS.with(|count|count.set(count.get()+1));
    let value=eval_c64_uncached(e);
    NUMERIC_MEMO.with(|memo|{memo.borrow_mut().values.get_or_insert_with(rustc_hash::FxHashMap::default).insert(key,(e.clone(),value));});
    value
}
fn eval_c64_uncached(e: &Expr) -> Option<C64> {
    Some(match e.kind() {
        Kind::Num(n) => C64::new(n.to_f64(), 0.0),
        Kind::Const(c) => match c {
            Constant::Pi => C64::new(std::f64::consts::PI, 0.0),
            Constant::E => C64::new(std::f64::consts::E, 0.0),
            Constant::I => C64::new(0.0, 1.0),
            Constant::EulerGamma => C64::new(EULER_GAMMA, 0.0),
            Constant::Catalan | Constant::GoldenRatio | Constant::TribonacciConstant => C64::new(c.approx_f64()?, 0.0),
            Constant::Infinity => C64::new(f64::INFINITY, 0.0),
            Constant::NegInfinity => C64::new(f64::NEG_INFINITY, 0.0),
            _ => return None,
        },
        Kind::Add(c, terms) => {
            let mut s = C64::new(c.to_f64(), 0.0);
            for (t, k) in terms {
                s += eval_c64(t)? * k.to_f64();
            }
            s
        }
        Kind::Mul(c, facs) => {
            let mut p = C64::new(c.to_f64(), 0.0);
            for (b, x) in facs {
                p *= cpow(eval_c64(b)?, x)?;
            }
            p
        }
        Kind::Pow(b, x) => cpow(eval_c64(b)?, x)?,
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>eval_c64(&crate::matrixcontraction::dot_product(a,b,true).ok()?)?,
        Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product,_)=>{let evaluated=crate::functions::doit(e,true);if evaluated==*e{return None;}eval_c64(&evaluated)?},
        Kind::Func(Func::Special(s), args) => crate::special::numeric::eval_c64(*s, args)?,
        Kind::Func(f, args) => {
            let a: Option<Vec<C64>> = args.iter().map(eval_c64).collect();
            apply_c64(f, &a?)?
        }
        _ => return None,
    })
}

fn cpow(b: C64, x: &Expr) -> Option<C64> {
    if let Some(k) = x.as_i64() {
        if k.abs() < 1 << 20 {
            return Some(b.powi(k as i32));
        }
    }
    let xe = eval_c64(x)?;
    if b.im == 0.0 && b.re >= 0.0 && xe.im == 0.0 {
        return Some(C64::new(b.re.powf(xe.re), 0.0));
    }
    if b == C64::new(0.0, 0.0) {
        return Some(C64::new(0.0, 0.0));
    }
    Some(b.powc(xe))
}

/// Real value of a numeric expression, if its imaginary part vanishes.
pub fn eval_f64(e: &Expr) -> Option<f64> {
    let z = eval_c64(e)?;
    if z.im.abs() <= 1e-14 * z.re.abs().max(1.0) {
        Some(z.re)
    } else {
        None
    }
}

fn real(z: C64) -> Option<f64> {
    if z.im == 0.0 {
        Some(z.re)
    } else {
        None
    }
}

pub fn apply_c64(f: &Func, a: &[C64]) -> Option<C64> {
    let x = a.first().copied()?;
    let realonly = |v: f64| C64::new(v, 0.0);
    Some(match f {
        Func::Exp => x.exp(),
        Func::Log => x.ln(),
        Func::Sin => x.sin(),
        Func::Cos => x.cos(),
        Func::Tan => x.tan(),
        Func::Cot => x.tan().inv(),
        Func::Sec => x.cos().inv(),
        Func::Csc => x.sin().inv(),
        Func::Asin => {
            if x.im == 0.0 && x.re.abs() <= 1.0 {
                realonly(x.re.asin())
            } else {
                x.asin()
            }
        }
        Func::Acos => {
            if x.im == 0.0 && x.re.abs() <= 1.0 {
                realonly(x.re.acos())
            } else {
                x.acos()
            }
        }
        Func::Atan => x.atan(),
        Func::Acot => x.inv().atan(),
        Func::Atan2 => realonly(real(a[0])?.atan2(real(a[1])?)),
        Func::Sinh => x.sinh(),
        Func::Cosh => x.cosh(),
        Func::Tanh => x.tanh(),
        Func::Coth => x.tanh().inv(),
        Func::Asinh => x.asinh(),
        Func::Acosh => {
            if x.im == 0.0 && x.re >= 1.0 {
                realonly(x.re.acosh())
            } else {
                x.acosh()
            }
        }
        Func::Atanh => {
            if x.im == 0.0 && x.re.abs() < 1.0 {
                realonly(x.re.atanh())
            } else {
                x.atanh()
            }
        }
        Func::Abs => realonly(x.norm()),
        Func::Sign => {
            if x.im == 0.0 {
                realonly(if x.re > 0.0 {
                    1.0
                } else if x.re < 0.0 {
                    -1.0
                } else {
                    0.0
                })
            } else {
                x / x.norm()
            }
        }
        Func::Floor => realonly(real(x)?.floor()),
        Func::Ceiling => realonly(real(x)?.ceil()),
        Func::Gamma => crate::special::numeric::gamma::gamma_c(x)?,
        Func::LogGamma => crate::special::numeric::gamma::loggamma_c(x)?,
        Func::Erf => crate::special::numeric::erf::erf_c(x),
        Func::Erfc => crate::special::numeric::erf::erfc_c(x),
        Func::Re => realonly(x.re),
        Func::Im => realonly(x.im),
        Func::Conjugate => x.conj(),
        Func::Arg => realonly(x.arg()),
        Func::Heaviside => {
            let r = real(x)?;
            realonly(if r > 0.0 {
                1.0
            } else if r < 0.0 {
                0.0
            } else {
                0.5
            })
        }
        Func::DiracDelta => {
            if x == C64::new(0.0, 0.0) {
                return None;
            }
            realonly(0.0)
        }
        Func::Max | Func::Min => {
            let mut best = real(a[0])?;
            for v in &a[1..] {
                let v = real(*v)?;
                best = if matches!(f, Func::Max) { best.max(v) } else { best.min(v) };
            }
            realonly(best)
        }
        Func::Predicate(_) | Func::AppliedPredicate | Func::Order | Func::Undefined(_) | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::Logic(_) | Func::Set(_)
        | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) => return None,
        Func::ExpPolar => {
            if x.im <= -std::f64::consts::PI || x.im > std::f64::consts::PI { return None; }
            x.exp()
        },
        Func::GenericExpr => return None,
        Func::Special(s) => return crate::special::numeric::apply(*s, a),
    })
}

/// Convert a complex double into an expression (`a + b*I`).
pub fn c64_to_expr(z: C64) -> Expr {
    if z.im == 0.0 {
        return Expr::float(z.re);
    }
    let im = Expr::from_kind(Kind::Mul(Number::float(z.im), vec![(Expr::i(), Expr::one())]));
    if z.re == 0.0 {
        return im;
    }
    add(&[Expr::float(z.re), im])
}

/// Evaluate `f(args)` numerically if possible (used for float arguments).
pub fn eval_func_numeric(f:&Func,args:&[Expr])->Option<Expr>{
    with_numeric_cache(||eval_func_numeric_inner(f,args))
}
fn eval_func_numeric_inner(f: &Func, args: &[Expr]) -> Option<Expr> {
    let a: Option<Vec<C64>> = args.iter().map(eval_c64).collect();
    let z = apply_c64(f, &a?)?;
    if !z.re.is_finite() || !z.im.is_finite() {
        return None;
    }
    Some(c64_to_expr(z))
}

/// SymPy's `evalf`: evaluate numbers and constants to floats, keeping
/// symbols. Integer exponents are kept exact.
pub fn evalf(e: &Expr) -> Expr {
    if is_numeric_value(e) {
        if let Some(z) = eval_c64(e) {
            if z.re.is_finite() && z.im.is_finite() {
                // clean tiny imaginary parts produced by rounding
                let z = if z.im.abs() < 1e-15 * z.re.abs() {
                    C64::new(z.re, 0.0)
                } else {
                    z
                };
                return c64_to_expr(z);
            }
        }
    }
    match e.kind() {
        Kind::Num(n) => Expr::float(n.to_f64()),
        Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => e.clone(),
        Kind::Add(..) => add(&e.args().iter().map(evalf).collect::<Vec<_>>()),
        Kind::Mul(..) => mul(&e.args().iter().map(evalf).collect::<Vec<_>>()),
        Kind::Pow(b, x) => {
            let x2 = if x.is_integer_number() { x.clone() } else { evalf(x) };
            pow(&evalf(b), &x2)
        }
        Kind::Func(f, args) => func(f.clone(), args.iter().map(evalf).collect()),
        Kind::Derivative(..) => e.clone(),
        Kind::Rel(op, a, b) => crate::relational::rel(*op, &evalf(a), &evalf(b)),
        Kind::Piecewise(p) => crate::relational::piecewise(
            p.iter().map(|(x, c)| (evalf(x), evalf(c))).collect(),
        ),
    }
}

/// Substitute numeric values for symbols and evaluate to a complex double.
pub fn eval_with(e: &Expr, env: &dyn Fn(&str) -> Option<C64>) -> Option<C64> {
    Some(match e.kind() {
        Kind::Sym(s) => env(s)?,
        Kind::Num(_) | Kind::Const(_) => eval_c64(e)?,
        Kind::Add(c, terms) => {
            let mut s = C64::new(c.to_f64(), 0.0);
            for (t, k) in terms {
                s += eval_with(t, env)? * k.to_f64();
            }
            s
        }
        Kind::Mul(c, facs) => {
            let mut p = C64::new(c.to_f64(), 0.0);
            for (b, x) in facs {
                let bv = eval_with(b, env)?;
                p *= match x.as_i64() {
                    Some(k) if k.abs() < 1 << 20 => bv.powi(k as i32),
                    _ => {
                        let xv = eval_with(x, env)?;
                        if bv.im == 0.0 && bv.re >= 0.0 && xv.im == 0.0 {
                            C64::new(bv.re.powf(xv.re), 0.0)
                        } else {
                            bv.powc(xv)
                        }
                    }
                };
            }
            p
        }
        Kind::Pow(b, x) => {
            let bv = eval_with(b, env)?;
            match x.as_i64() {
                Some(k) if k.abs() < 1 << 20 => bv.powi(k as i32),
                _ => {
                    let xv = eval_with(x, env)?;
                    if bv.im == 0.0 && bv.re >= 0.0 && xv.im == 0.0 {
                        C64::new(bv.re.powf(xv.re), 0.0)
                    } else {
                        bv.powc(xv)
                    }
                }
            }
        }
        Kind::Func(f, args) => {
            let a: Option<Vec<C64>> = args.iter().map(|a| eval_with(a, env)).collect();
            apply_c64(f, &a?)?
        }
        _ => return None,
    })
}

// ----------------------------------------------------------------------
// special functions (real)
// ----------------------------------------------------------------------

/// Gamma function (Lanczos approximation, g=7, n=9).
pub fn gamma_f64(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x.fract() == 0.0 && x <= 0.0 {
        return f64::NAN;
    }
    if x.fract() == 0.0 && x > 0.0 && x < 171.0 {
        let mut acc = 1.0;
        let mut i = 2.0;
        while i < x {
            acc *= i;
            i += 1.0;
        }
        return acc;
    }
    if x < 0.5 {
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * gamma_f64(1.0 - x))
    } else {
        let x = x - 1.0;
        let mut a = G[0];
        let t = x + 7.5;
        for (i, g) in G.iter().enumerate().skip(1) {
            a += g / (x + i as f64);
        }
        (2.0 * std::f64::consts::PI).sqrt() * t.powf(x + 0.5) * (-t).exp() * a
    }
}

/// Error function, accurate to ~1e-15.
pub fn erf_f64(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return -erf_f64(-x);
    }
    if x < 2.5 {
        // Taylor series: erf(x) = 2/sqrt(pi) * sum (-1)^n x^(2n+1) / (n! (2n+1))
        let mut sum = 0.0;
        let mut term = x;
        let x2 = x * x;
        let mut n = 0.0;
        loop {
            let t = term / (2.0 * n + 1.0);
            sum += t;
            if t.abs() < 1e-17 * sum.abs() {
                break;
            }
            n += 1.0;
            term *= -x2 / n;
            if n > 200.0 {
                break;
            }
        }
        sum * 2.0 / std::f64::consts::PI.sqrt()
    } else {
        // continued fraction for erfc
        if x > 27.0 {
            return 1.0;
        }
        let mut f = 0.0;
        for k in (1..60).rev() {
            f = (k as f64 / 2.0) / (x + f);
        }
        let erfc = (-x * x).exp() / (std::f64::consts::PI.sqrt() * (x + f));
        1.0 - erfc
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn special() {
        assert!((gamma_f64(5.0) - 24.0).abs() < 1e-12);
        assert!((gamma_f64(0.5) - std::f64::consts::PI.sqrt()).abs() < 1e-13);
        assert!((erf_f64(1.0) - 0.842_700_792_949_714_9).abs() < 1e-14);
        assert!((erf_f64(3.0) - 0.999_977_909_503_001_4).abs() < 1e-14);
    }
}

#[cfg(test)]
mod dag_cache_tests {
    use super::*;
    fn shared_mul(base:Expr,depth:usize)->Expr{
        let mut value=base;let one=Expr::one();
        for _ in 0..depth{value=Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,vec![(value.clone(),one.clone()),(value,one.clone())]));}value
    }
    #[test]
    fn closed_numeric_dag_is_evaluated_once_per_node(){
        let e=shared_mul(Expr::one(),80);EVALUATIONS.with(|count|count.set(0));
        assert_eq!(eval_c64(&e),Some(C64::new(1.0,0.0)));
        assert_eq!(EVALUATIONS.with(|count|count.get()),80);
        NUMERIC_MEMO.with(|memo|{let memo=memo.borrow();assert_eq!(memo.depth,0);assert!(memo.values.is_none());});
    }
    #[test]
    fn unknown_values_contexts_and_physical_identity_are_memoized(){
        let a=Expr::from_kind(Kind::Func(Func::Sin,vec![Expr::pi()]));
        let b=Expr::from_kind(Kind::Func(Func::Sin,vec![Expr::pi()]));
        assert_eq!(a,b);assert!(!a.ptr_eq(&b));
        let unknown=Expr::from_kind(Kind::Func(Func::Undefined("f".into()),vec![Expr::one()]));
        EVALUATIONS.with(|count|count.set(0));
        with_numeric_cache(||{
            assert!(eval_c64(&a).is_some());assert!(eval_c64(&a).is_some());assert_eq!(EVALUATIONS.with(|count|count.get()),1);
            assert!(eval_c64(&b).is_some());assert_eq!(EVALUATIONS.with(|count|count.get()),2);
            assert_eq!(eval_c64(&unknown),None);assert_eq!(eval_c64(&unknown),None);assert_eq!(EVALUATIONS.with(|count|count.get()),3);
            crate::assumptions::exact_queries(||{assert_eq!(eval_c64(&unknown),None);assert!(eval_c64(&a).is_some());});
            assert_eq!(EVALUATIONS.with(|count|count.get()),5);
            assert!(eval_c64(&a).is_some());assert_eq!(EVALUATIONS.with(|count|count.get()),5);
            NUMERIC_MEMO.with(|memo|assert_eq!(memo.borrow().values.as_ref().unwrap().len(),5));
        });
        NUMERIC_MEMO.with(|memo|assert!(memo.borrow().values.is_none()));
    }
    #[test]
    fn scoped_owners_and_panic_cleanup_are_correct(){
        let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||with_numeric_cache(||{
            // Only the memo owns each temporary after its statement returns.
            for _ in 0..32{let e=Expr::from_kind(Kind::Func(Func::Sin,vec![Expr::pi()]));assert!(eval_c64(&e).is_some());}
            NUMERIC_MEMO.with(|memo|assert_eq!(memo.borrow().values.as_ref().unwrap().len(),32));
            panic!("expected numeric unwind");
        })));
        assert!(result.is_err());
        NUMERIC_MEMO.with(|memo|{let memo=memo.borrow();assert_eq!(memo.depth,0);assert!(memo.values.is_none());});
        assert_eq!(eval_f64(&shared_mul(Expr::one(),80)),Some(1.0));
    }
    #[test]
    fn environments_are_separate_from_closed_numeric_memo(){
        let x=Expr::symbol("x");let e=Expr::from_kind(Kind::Func(Func::Sin,vec![x]));
        with_numeric_cache(||{
            assert_eq!(eval_c64(&e),None);
            let first=eval_with(&e,&|name|(name=="x").then_some(C64::new(0.0,0.0)));
            let second=eval_with(&e,&|name|(name=="x").then_some(C64::new(std::f64::consts::FRAC_PI_2,0.0)));
            assert_eq!(first,Some(C64::new(0.0,0.0)));assert_eq!(second,Some(C64::new(1.0,0.0)));
            assert_eq!(eval_c64(&e),None);
        });
    }
    #[test]
    fn ambiguous_numeric_assumptions_do_not_use_machine_refinement() {
        let mut value=Expr::one();
        for _ in 0..80 {
            value=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![
                (value.clone(),Number::ONE),(value,Number::MINUS_ONE),
            ]));
        }
        EVALUATIONS.with(|count|count.set(0));
        assert_eq!(crate::assumptions::is_positive(&value),None);
        // Assumptions now require increasing multiprecision evaluations.
        assert_eq!(EVALUATIONS.with(|count|count.get()),0);
        NUMERIC_MEMO.with(|memo|{let memo=memo.borrow();assert_eq!(memo.depth,0);assert!(memo.values.is_none());});
    }

    #[test]
    fn multi_argument_numeric_entries_share_their_scope(){
        let e=shared_mul(Expr::one(),80);let args=vec![e.clone(),e];
        EVALUATIONS.with(|count|count.set(0));
        assert!(eval_func_numeric(&Func::Atan2,&args).is_some());
        assert_eq!(EVALUATIONS.with(|count|count.get()),80);
        EVALUATIONS.with(|count|count.set(0));
        assert!(crate::special::numeric::eval_c64(crate::special::SpecialFn::Beta,&args).is_some());
        assert_eq!(EVALUATIONS.with(|count|count.get()),80);
    }
}
