//! Arbitrary precision numerical evaluation of expressions.

use num_bigint::BigInt;
use num_traits::ToPrimitive;

use crate::bigfloat::{ln2, pi, BigFloat};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
mod memo;

/// A complex number with arbitrary precision parts.
#[derive(Clone, Debug)]
pub struct MpComplex {
    pub re: BigFloat,
    pub im: BigFloat,
}

// Precision is part of the key: cancellation refinement can revisit the same
// physical node at several working precisions. Owners prevent address reuse.
// The options are part of the key as well: a chopped or budgeted evaluation
// inside an enclosing numerical cache must not answer an unchopped query.
struct NumericMemo {
    depth: usize,
    values: Option<rustc_hash::FxHashMap<(usize,u32,bool,EvalfOptions),(Expr,Option<MpComplex>)>>,
}
thread_local! {
    static NUMERIC_MEMO:std::cell::RefCell<NumericMemo> = const {
        std::cell::RefCell::new(NumericMemo{depth:0,values:None})
    };
    #[cfg(test)] static EVALUATIONS:std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static OPTIONS:std::cell::Cell<EvalfOptions> = const { std::cell::Cell::new(EvalfOptions{chop:None,budget:None,partial:false,residue:false,cancellation:None}) };
    // The sum whose exhaustion returns its residue (see `EvalfOptions::residue`)
    // and the accuracies, in bits, of the real and imaginary parts it returned.
    static RESIDUE_ROOT:std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static RESIDUE:std::cell::Cell<Option<(i128,i128)>> = const { std::cell::Cell::new(None) };
    static ACHIEVED:std::cell::Cell<Option<u32>> = const { std::cell::Cell::new(None) };
    static EXHAUSTED:std::cell::RefCell<Option<Expr>> = const { std::cell::RefCell::new(None) };
}

/// Options of one public numerical evaluation, as in SymPy's `evalf`.
///
/// `budget` is the number of bits a cancelling sum may add to its target
/// precision (SymPy's `maxprec = max(prec, int(maxn*LG10))`). Without it, a
/// sum may add up to `DEFAULT_CANCELLATION_BUDGET` bits. `chop` is the binary
/// magnitude below which a numerical component is removed; inside a sum it
/// also accepts a component whose rounding error is below that magnitude.
/// With `partial`, a sum that exhausts its budget with at least one correct
/// bit is returned at the accuracy it reached, as SymPy returns a Float of
/// lower precision for an explicit `maxn`; otherwise it stays unevaluated.
/// With `residue` (the public `evalf` without `strict`), a sum evaluated as
/// a whole that exhausts its budget returns, as SymPy's `evalf_add`, its
/// computed value with the accuracy of each part: a part without a correct
/// bit becomes a precision 1 Float, the rounding residue or, for an exact
/// zero, the absolute error bound (SymPy's scaled zero). The default budget
/// is tried first; its exhaustion is then reported with SymPy's default
/// `maxn=100` schedule, so a cancellation resolved within the larger budget
/// keeps its digits.
#[derive(Clone,Copy,Debug,Default,PartialEq,Eq,Hash)]
pub struct EvalfOptions { pub chop:Option<i128>, pub budget:Option<u32>, pub partial:bool, pub residue:bool,
    /// Internal: replaces `DEFAULT_CANCELLATION_BUDGET` for a sum without an
    /// explicit `budget` (same schedule, smaller ceiling).
    pub cancellation:Option<u32> }

pub const DEFAULT_CANCELLATION_BUDGET:u32=4096;
/// SymPy's last restart may overshoot the budget by the missing accuracy;
/// this bounds that overshoot so a huge cancelling term cannot exhaust memory.
const EXPLICIT_WORKING_LIMIT:u32=1<<20;

/// Restores the previous evaluation options when an evaluation ends.
struct OptionsScope(EvalfOptions,Option<Expr>,Option<u32>);
impl OptionsScope {
    fn enter(options:EvalfOptions)->Self{
        Self(OPTIONS.with(|cell|cell.replace(options)),EXHAUSTED.with(|cell|cell.borrow_mut().take()),ACHIEVED.with(|cell|cell.take()))
    }
}
impl Drop for OptionsScope {
    fn drop(&mut self){
        OPTIONS.with(|cell|cell.set(self.0));
        ACHIEVED.with(|cell|cell.set(self.2));
        let previous=self.1.take();
        EXHAUSTED.with(|cell|*cell.borrow_mut()=previous);
    }
}

/// Evaluate `e` with explicit options. The second component is the first sum
/// whose cancellation could not be resolved within the budget, when that
/// left the result incompletely evaluated (SymPy's strict mode raises
/// `PrecisionExhausted` for it).
pub fn evalf_with_options(e:&Expr,prec:u32,options:EvalfOptions)->(Expr,Option<Expr>){
    let _scope=OptionsScope::enter(options);
    let mut value=eval_mp_root(e,prec);
    let residue=RESIDUE.with(|cell|cell.take());
    let complete=value.is_some();
    let mut prec=prec;
    if let (Some(achieved),Some(value))=(ACHIEVED.with(|cell|cell.get()),value.as_mut()) {
        if achieved<prec {
            prec=achieved;
            *value=MpComplex{re:value.re.round_to(prec),im:value.im.round_to(prec)};
        }
    }
    let result=match (value,options.chop) {
        (Some(mut value),Some(cutoff))=>{
            if !value.re.is_zero() && value.re.mag()<cutoff {value.re=BigFloat::zero();}
            if !value.im.is_zero() && value.im.mag()<cutoff {value.im=BigFloat::zero();}
            match residue {Some(accuracies)=>residue_expr(value,prec,accuracies),None=>to_expr(value,prec)}
        },
        (Some(value),None)=>match residue {Some(accuracies)=>residue_expr(value,prec,accuracies),None=>to_expr(value,prec)},
        (None,_)=>evalf(e,prec),
    };
    let exhausted=if complete {None} else {EXHAUSTED.with(|cell|cell.borrow_mut().take())};
    (result,exhausted)
}

/// The public `evalf` without options: `evalf`, where a cancelling sum that
/// exhausts its budget returns its residue (see `EvalfOptions::residue`).
pub fn evalf_reporting_residue(e:&Expr,prec:u32)->Expr{
    let _scope=OptionsScope::enter(EvalfOptions{residue:true,..EvalfOptions::default()});
    evalf(e,prec)
}

/// Evaluate with the default options, independently of an enclosing
/// evaluation's options. Sign and magnitude decisions use this.
pub fn eval_mp_default(e:&Expr,prec:u32)->Option<MpComplex>{
    let _scope=OptionsScope::enter(EvalfOptions::default());
    let _errors=HostErrorScope::enter(false);
    eval_mp(e,prec)
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

/// Share precision-specific values across one outer inference operation.
pub(crate) fn with_numeric_cache<R>(calculation: impl FnOnce() -> R) -> R {
    let _scope=NumericScope::enter();
    calculation()
}

impl MpComplex {
    fn sqrt(&self,p:u32)->Option<MpComplex>{
        if self.is_real(){
            return Some(if self.re.sign()<0{
                MpComplex{re:BigFloat::zero(),im:self.re.neg().sqrt(p)?}
            }else{MpComplex::real(self.re.sqrt(p)?)});
        }
        let working=p+20;
        // The large component uses a sum; recover the small component from
        // 2*Re(sqrt(z))*Im(sqrt(z))=Im(z), avoiding cancellation near a cut.
        let large=self.abs(working).add(&self.re.abs(),working).ldexp(-1).sqrt(working)?;
        let small=self.im.div(&large.ldexp(1),working)?;
        Some(if self.re.sign()<0{
            MpComplex{re:small.abs().round_to(p),im:if self.im.sign()<0{large.neg()}else{large}.round_to(p)}
        }else{MpComplex{re:large.round_to(p),im:small.round_to(p)}})
    }
    pub fn real(re: BigFloat) -> MpComplex {
        MpComplex { re, im: BigFloat::zero() }
    }

    pub fn is_real(&self) -> bool {
        self.im.is_zero()
    }

    pub(crate) fn add(&self, o: &MpComplex, p: u32) -> MpComplex {
        MpComplex { re: self.re.add(&o.re, p), im: self.im.add(&o.im, p) }
    }

    pub(crate) fn mul(&self, o: &MpComplex, p: u32) -> MpComplex {
        if self.is_real() && o.is_real() {
            return MpComplex::real(self.re.mul(&o.re, p));
        }
        let re = self.re.mul(&o.re, p).sub(&self.im.mul(&o.im, p), p);
        let im = self.re.mul(&o.im, p).add(&self.im.mul(&o.re, p), p);
        MpComplex { re, im }
    }

    pub(crate) fn scale(&self, k: &BigFloat, p: u32) -> MpComplex {
        MpComplex { re: self.re.mul(k, p), im: self.im.mul(k, p) }
    }

    pub(crate) fn inv(&self, p: u32) -> Option<MpComplex> {
        if self.is_real() {
            return Some(MpComplex::real(BigFloat::one().div(&self.re, p)?));
        }
        let d = self.re.mul(&self.re, p + 10).add(&self.im.mul(&self.im, p + 10), p + 10);
        Some(MpComplex { re: self.re.div(&d, p)?, im: self.im.neg().div(&d, p)? })
    }

    pub(crate) fn abs(&self, p: u32) -> BigFloat {
        if self.is_real() {
            return self.re.abs();
        }
        self.re
            .mul(&self.re, p + 10)
            .add(&self.im.mul(&self.im, p + 10), p + 10)
            .sqrt(p)
            .unwrap()
    }

    pub(crate) fn neg(&self) -> MpComplex {
        MpComplex { re: self.re.neg(), im: self.im.neg() }
    }

    pub(crate) fn sub(&self, other: &MpComplex, precision: u32) -> MpComplex {
        self.add(&other.neg(), precision)
    }

    pub(crate) fn div(&self, other: &MpComplex, precision: u32) -> Option<MpComplex> {
        Some(self.mul(&other.inv(precision + 8)?, precision))
    }

    pub(crate) fn conjugate(&self) -> MpComplex {
        MpComplex { re: self.re.clone(), im: self.im.neg() }
    }

    pub(crate) fn round_to(&self, precision: u32) -> MpComplex {
        MpComplex { re: self.re.round_to(precision), im: self.im.round_to(precision) }
    }

    fn arg(&self, p: u32) -> BigFloat {
        BigFloat::atan2(&self.im, &self.re, p)
    }

    pub(crate) fn exp(&self, p: u32) -> Option<MpComplex> {
        // mpmath exponents are unbounded. Binary exponents here are i128
        // values bounded by MAX_EXPONENT, which covers exp(2**63) (Torch
        // value ranges compare exp(9223372036854775807)); beyond it the
        // value stays symbolic.
        let r = self.re.checked_exp(p + 10)?;
        if self.im.is_zero() {
            return Some(MpComplex::real(r.round_to(p)));
        }
        let (s, c) = self.im.sin_cos(p + 10);
        Some(MpComplex { re: r.mul(&c, p), im: r.mul(&s, p) })
    }

    pub(crate) fn ln(&self, p: u32) -> Option<MpComplex> {
        crate::elementary_mp::log(self, p)
    }

    fn pow(&self, e: &MpComplex, p: u32) -> Option<MpComplex> {
        if self.re.is_zero() && self.im.is_zero() {
            return if e.re.sign() > 0 { Some(MpComplex::real(BigFloat::zero())) } else { None };
        }
        if self.is_real() && e.is_real() && self.re.sign() > 0 {
            return Some(MpComplex::real(self.re.powf(&e.re, p)?));
        }
        if e.re.mag() > 1 << 20 || e.im.mag() > 1 << 20 { return None; }
        let wp = p + 20 + e.re.mag().max(0) as u32 + e.im.mag().max(0) as u32;
        let l = self.ln(wp)?;
        l.mul(e, wp).exp(p)
    }

    fn powi(&self, n: i64, p: u32) -> Option<MpComplex> {
        if self.is_real() {
            return Some(MpComplex::real(self.re.powi(n, p)?));
        }
        let wp = p + 10 + 64 - n.unsigned_abs().leading_zeros();
        let mut base = self.clone();
        let mut k = n.unsigned_abs();
        let mut acc = MpComplex::real(BigFloat::one());
        while k > 0 {
            if k & 1 == 1 {
                acc = acc.mul(&base, wp);
                if !supported_magnitude(&acc) { return None; }
            }
            k >>= 1;
            if k > 0 {
                base = base.mul(&base, wp);
                if !supported_magnitude(&base) { return None; }
            }
        }
        if n < 0 {
            acc.inv(p)
        } else {
            Some(acc)
        }
    }
}

fn supported_magnitude(value: &MpComplex) -> bool {
    [&value.re, &value.im].iter().all(|v| v.is_zero() || v.mag().unsigned_abs() <= crate::bigfloat::MAX_EXPONENT as u128)
}

fn number(n: &Number, p: u32) -> BigFloat {
    match n {
        Number::Int(i) => BigFloat::from_i64(*i),
        Number::BigInt(b) => BigFloat::from_int(b),
        Number::Rational(r) => BigFloat::from_ratio(r.numer(), r.denom(), p),
        Number::Float(f) => BigFloat::from_f64(*f),
        Number::MpFloat(value, _) => value.round_to(p),
    }
}

/// Highest declared precision among the floating point leaves.
pub fn max_precision(e: &Expr) -> Option<u32> {
    e.float_info().max_precision
}

pub fn has_mp(e: &Expr) -> bool {
    e.float_info().has_mp
}

/// Insert scalar symbol values without prematurely rounding their arithmetic.
/// Sequential mapping semantics are retained, including substitutions into an
/// earlier replacement value. Binding and matrix expressions use the existing
/// algebraic substitution implementation rather than this scalar traversal.
pub fn substitute_for_evalf(e: &Expr, pairs: &[(Expr, Expr)]) -> Expr {
    use rustc_hash::FxHashMap;
    fn scalar(e: &Expr, seen: &mut FxHashMap<usize, Expr>) -> bool {
        if seen.insert(e.dag_id(), e.clone()).is_some() { return true; }
        if e.has_matrix() { return false; }
        match e.kind() {
            Kind::Derivative(..) | Kind::Rel(..) | Kind::Piecewise(..) | Kind::Bool(..) => false,
            Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product
                | Func::Set(..) | Func::Logic(..), _) => false,
            _ => crate::inspection::scalar_children(e).iter().all(|child| scalar(child, seen)),
        }
    }
    fn replace(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<usize, (Expr, Expr)>) -> Expr {
        if e == old { return new.clone(); }
        if e.is_atom() { return e.clone(); }
        if let Some((_, value)) = cache.get(&e.dag_id()) { return value.clone(); }
        let mut child = |value: &Expr| replace(value, old, new, cache);
        let kind = match e.kind() {
            Kind::Add(c, terms) => Kind::Add(c.clone(), terms.iter().map(|(t, k)| (child(t), k.clone())).collect()),
            Kind::Mul(c, factors) => Kind::Mul(c.clone(), factors.iter().map(|(b, p)| (child(b), child(p))).collect()),
            Kind::Pow(b, p) => Kind::Pow(child(b), child(p)),
            Kind::Func(f, args) => Kind::Func(f.clone(), args.iter().map(&mut child).collect()),
            Kind::Rel(op, a, b) => Kind::Rel(*op, child(a), child(b)),
            Kind::Piecewise(branches) => Kind::Piecewise(branches.iter().map(|(v, c)| (child(v), child(c))).collect()),
            _ => return e.clone(),
        };
        let result = Expr::from_unevaluated_kind(kind);
        let result = if result == *e { e.clone() } else { result };
        // Retain source nodes as well as results: a pointer-key cache must not
        // outlive its original node during this transformation.
        cache.insert(e.dag_id(), (e.clone(), result.clone()));
        result
    }
    if pairs.is_empty() { return e.clone(); }
    let mut seen = FxHashMap::default();
    if pairs.iter().any(|(old, _)| !old.is_symbol()) || !scalar(e, &mut seen)
        || pairs.iter().any(|(_, value)| !scalar(value, &mut seen)) {
        return crate::subs::subs(e, pairs);
    }
    pairs.iter().fold(e.clone(), |current, (old, new)| replace(&current, old, new, &mut FxHashMap::default()))
}

/// Evaluate a numeric expression to `prec` bits. Returns `None` when the
/// expression has free symbols or uses an unsupported function.
pub fn eval_mp(e: &Expr, prec: u32) -> Option<MpComplex> {
    eval_mp_bounded(e, prec, true)
}

/// `eval_mp` for callers that only inspect signs and magnitudes: values
/// beyond the decimal bound of materialized Floats are kept, as mpmath keeps
/// them when SymPy compares `exp(16777216)` with a number.
pub(crate) fn eval_mp_unbounded(e: &Expr, prec: u32) -> Option<MpComplex> {
    eval_mp_bounded(e, prec, false)
}

fn eval_mp_bounded(e: &Expr, prec: u32, bounded: bool) -> Option<MpComplex> {
    if let Some((re,im))=crate::complex_parts::literal_parts(e){
        let component=|e:&Expr|match e.as_number()?{
            Number::MpFloat(value,_)=>value.checked_round_to(prec),
            n=>number(n,prec).checked_round_to(prec),
        };
        return Some(MpComplex{re:component(&re)?,im:component(&im)?});
    }
    if let Kind::Func(Func::MatrixDotProduct(a,b),_)=e.kind(){return eval_mp_bounded(&crate::matrixcontraction::dot_product(a,b,true).ok()?,prec,bounded);}
    if matches!(e.kind(),Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product,_)){
        let evaluated=crate::functions::doit(e,true);
        if evaluated==*e{return None;}
        return eval_mp_bounded(&evaluated,prec,bounded);
    }
    let wp = prec + 20;
    let r = ev(e, wp)?;
    if bounded && !supported_magnitude(&r) { return None; }
    // SymPy's evalf computes at prec + 4 bits and rounds that value to prec.
    Some(MpComplex { re: r.re.round_to(prec + 4).round_to(prec), im: r.im.round_to(prec + 4).round_to(prec) })
}

/// Convert a computed value back to native expression numbers, retaining the
/// declared precision. Exact zero components stay exact zero.
/// The numerical value of `e` at `prec` bits as a Float expression.
pub fn eval_number(e: &Expr, prec: u32) -> Option<Expr> {
    eval_mp(e, prec).map(|value| to_expr(value, prec))
}

/// `eval_mp` of `e` evaluated as a whole: under explicit options, a sum `e`
/// follows SymPy's schedule from its target and, under
/// `EvalfOptions::residue`, may return its residue.
fn eval_mp_root(e:&Expr,prec:u32)->Option<MpComplex>{
    if OPTIONS.with(|cell|cell.get())==EvalfOptions::default() {return eval_mp(e,prec);}
    struct Root(usize);
    impl Drop for Root { fn drop(&mut self){RESIDUE_ROOT.with(|cell|cell.set(self.0));} }
    let _root=Root(RESIDUE_ROOT.with(|cell|cell.replace(e.dag_id())));
    RESIDUE.with(|cell|cell.set(None));
    let value=eval_mp(e,prec);
    if value.is_none() {RESIDUE.with(|cell|cell.set(None));}
    value
}

/// The Floats of an exhausted sum: as SymPy's evalf, each part has precision
/// `max(min(prec, accuracy), 1)`.
fn residue_expr(value:MpComplex,prec:u32,(real,imaginary):(i128,i128))->Expr{
    let precision=|accuracy:i128|u32::try_from(accuracy.clamp(1,i128::from(prec))).unwrap_or(prec);
    let re=if value.re.is_zero(){Expr::zero()}else{Expr::num(Number::mp_float(value.re,precision(real)))};
    if value.im.is_zero(){return re;}
    let im=Expr::from_kind(Kind::Mul(Number::mp_float(value.im,precision(imaginary)),vec![(Expr::i(),Expr::one())]));
    if re.is_zero(){im}else{crate::ops::add2(&re,&im)}
}

pub fn to_expr(value: MpComplex, prec: u32) -> Expr {
    let re = if value.re.is_zero() { Expr::zero() } else { Expr::num(Number::mp_float(value.re, prec)) };
    if value.im.is_zero() { return re; }
    let im = Expr::from_kind(Kind::Mul(Number::mp_float(value.im, prec), vec![(Expr::i(), Expr::one())]));
    if re.is_zero() { im } else { crate::ops::add2(&re, &im) }
}

/// Numerical evaluation with native arbitrary precision numbers. Unsupported
/// numerical functions remain symbolic; this never silently computes them at
/// machine precision. Integer powers and symbolic function arguments remain
/// exact, as in SymPy's partial evalf.
pub fn evalf(e: &Expr, prec: u32) -> Expr {
    if e.is_zero() { return Expr::zero(); }
    let options=OPTIONS.with(|cell|cell.get());
    let memoized=EvalfOptions{residue:false,..options}==EvalfOptions::default();
    if memoized { if let Some(value) = memo::get(e,prec,options.residue) { return value; } }
    if let Some(value) = eval_mp_root(e, prec) {
        let value=match RESIDUE.with(|cell|cell.take()) {
            Some(accuracies)=>residue_expr(value,prec,accuracies),
            None=>to_expr(value,prec),
        };
        if memoized {memo::insert(e,&value,prec,options.residue);}
        return value;
    }
    if lower_gamma_pole(e) { return Expr::oo(); }
    if crate::functions::is_numeric_value(e) { return e.clone(); }
    match e.kind() {
        Kind::Add(..) => crate::ops::add(&e.args().iter().map(|a| evalf(a, prec)).collect::<Vec<_>>()),
        Kind::Mul(..) => crate::ops::mul(&e.args().iter().map(|a| evalf(a, prec)).collect::<Vec<_>>()),
        Kind::Pow(base, exponent) => crate::ops::pow(
            &evalf(base, prec),
            &if exponent.is_integer_number() { exponent.clone() } else { evalf(exponent, prec) },
        ),
        Kind::Rel(op, a, b) => crate::relational::rel(*op, &evalf(a, prec), &evalf(b, prec)),
        Kind::Piecewise(pairs) => crate::relational::piecewise(
            pairs.iter().map(|(value, condition)| (evalf(value, prec), condition.clone())).collect(),
        ),
        Kind::Func(Func::Set(crate::sets::SetOp::Lambda),arguments) if arguments.len()==2 =>
            crate::sets::lambda(&arguments[0],&evalf(&arguments[1],prec)),
        // SymPy evalf returns a non-numerical `_eval_evalf` result (such as
        // oo for a Number class whose `_as_mpf_val` is infinite) as it is.
        Kind::Func(Func::Undefined(head),_) if head.host().is_some_and(|host| host.evaluates_numerically()) =>
            match head.host().and_then(|host| host.evalf(e, prec)) {
                Some(value) if !crate::traverse::has(&value, e) => evalf(&value, prec),
                Some(value) => value,
                None => e.clone(),
            },
        _ => e.clone(),
    }
}

/// SymPy's `lowergamma._eval_evalf` returns mpmath's `gammainc(a, 0, z)` as
/// it is: at a pole (a nonpositive integer `a`, `z != 0` with `re(z) >= 0`)
/// that is `+inf`, so the value is `oo`.
fn lower_gamma_pole(e: &Expr) -> bool {
    let Kind::Func(Func::Special(crate::special::SpecialFn::LowerGamma), args) = e.kind() else { return false };
    let [a, z] = args.as_slice() else { return false };
    if !a.as_number().is_some_and(|n| n.is_integer() && !n.is_positive()) {
        return false;
    }
    let Some(value) = ev(z, 64) else { return false };
    !(value.re.is_zero() && value.im.is_zero()) && value.re.sign() >= 0
}

/// Remove absolutely small real/imaginary numerical components by binary
/// magnitude. Partial symbolic coefficients retain ordinary evalf semantics.
/// The value of `e` at `prec` bits with numerical components below the
/// binary magnitude `cutoff` removed, as SymPy's `n(chop=True)`. A sum whose
/// component cancels exactly is resolved once its rounding error falls below
/// the cutoff, where `eval_mp` keeps refining it until its budget is spent.
/// `eval_mp` where a sum without an explicit budget may add at most
/// `cancellation` bits to its target precision before it stays unevaluated.
/// For callers that compare the value with a tolerance: a cancellation of
/// more than that many bits already places the sum below it.
pub(crate) fn eval_mp_cancelling(e:&Expr,prec:u32,cancellation:u32)->Option<MpComplex>{
    let options=OPTIONS.with(|cell|cell.get());
    let _scope=OptionsScope::enter(EvalfOptions{cancellation:Some(cancellation),..options});
    eval_mp(e,prec)
}

pub fn eval_mp_chopped(e:&Expr,prec:u32,cutoff:i128)->Option<MpComplex>{
    let _scope=OptionsScope::enter(EvalfOptions{chop:Some(cutoff),budget:None,partial:false,residue:false,cancellation:None});
    let _errors=HostErrorScope::enter(false);
    let mut value=eval_mp(e,prec)?;
    if !value.re.is_zero() && value.re.mag()<cutoff {value.re=BigFloat::zero();}
    if !value.im.is_zero() && value.im.mag()<cutoff {value.im=BigFloat::zero();}
    Some(value)
}

pub fn evalf_chop(e:&Expr,prec:u32,cutoff:i128)->Expr {
    evalf_with_options(e,prec,EvalfOptions{chop:Some(cutoff),budget:None,partial:false,residue:false,cancellation:None}).0
}

pub(crate) fn ev(e: &Expr, p: u32) -> Option<MpComplex> {
    if e.is_atom(){return ev_uncached(e,p);}
    let _scope=NumericScope::enter();
    let key=(e.dag_id(),p,crate::assumptions::numeric_refinement_enabled(),OPTIONS.with(|cell|cell.get()));
    if let Some(value)=NUMERIC_MEMO.with(|memo|memo.borrow().values.as_ref()
        .and_then(|values|values.get(&key).map(|(_,value)|value.clone()))) {
        return value;
    }
    #[cfg(test)] EVALUATIONS.with(|count|count.set(count.get()+1));
    let result=ev_uncached(e,p);
    NUMERIC_MEMO.with(|memo|memo.borrow_mut().values
        .get_or_insert_with(rustc_hash::FxHashMap::default).insert(key,(e.clone(),result.clone())));
    result
}

fn ev_uncached(e: &Expr, p: u32) -> Option<MpComplex> {
    if e.has_wide_exponent()&&matches!(e.kind(),Kind::Add(..)|Kind::Mul(..)|Kind::Pow(..)){
        return None;
    }
    Some(match e.kind() {
        Kind::Num(n) => MpComplex::real(number(n, p)),
        Kind::Const(c) => match c {
            Constant::Pi => MpComplex::real(pi(p)),
            Constant::E => MpComplex::real(BigFloat::one().exp(p)),
            Constant::I => MpComplex { re: BigFloat::zero(), im: BigFloat::one() },
            Constant::EulerGamma => MpComplex::real(euler_gamma(p)),
            Constant::Catalan => MpComplex::real(catalan(p)),
            Constant::GoldenRatio => MpComplex::real(golden_ratio(p)),
            Constant::TribonacciConstant => MpComplex::real(tribonacci_constant(p)),
            _ => return None,
        },
        Kind::Add(c, terms) => eval_sum(e,c, terms, p)?,
        Kind::Mul(c, facs) => {
            let mut acc = MpComplex::real(number(c, p));
            for (b, x) in facs {
                let v = pow_of(b, x, p)?;
                acc = acc.mul(&v, p);
            }
            acc
        }
        Kind::Pow(b, x) => pow_of(b, x, p)?,
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>ev(&crate::matrixcontraction::dot_product(a,b,true).ok()?,p)?,
        Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product,_)=>{let evaluated=crate::functions::doit(e,true);if evaluated==*e{return None;}ev(&evaluated,p)?},
        Kind::Func(Func::Undefined(head),_) if head.host().is_some_and(|host| host.evaluates_numerically()) => {
            // SymPy evalf falls back to the class's _eval_evalf.
            let value=head.host()?.evalf(e,p)?;
            ev(&value,p)?
        },
        Kind::Func(f, args) => func(f, args, p)?,
        _ => return None,
    })
}

/// Increase working precision when large terms nearly cancel. A fixed guard
/// budget would incorrectly turn sqrt(10**100 + 1) - 10**50 into exact zero.
fn eval_sum(expression:&Expr,constant: &Number, terms: &[(Expr, Number)], precision: u32) -> Option<MpComplex> {
    // Evaluation may collect identical finite operands even when construction
    // was deliberately held. This proves a-a without chasing a numerical zero
    // through ever increasing precisions on each shared-graph prefix.
    if !matches!(constant,Number::Float(value) if !value.is_finite()) {
        let cancelled=crate::assumptions::exact_queries(|| {
            let mut coefficients=rustc_hash::FxHashMap::<usize,Number>::default();
            for (term,coefficient) in terms {
                if coefficient.is_float()
                    || crate::assumptions::is_finite(term)!=Some(true) { return false; }
                let sum=coefficients.entry(term.dag_id()).or_insert(Number::ZERO);
                *sum=sum.add(coefficient);
            }
            coefficients.values().all(Number::is_zero)
        });
        if cancelled {return Some(MpComplex::real(number(constant,precision)));}
    }
    let known_real=crate::assumptions::exact_queries(||crate::assumptions::is_real(expression))==Some(true);
    if known_real&&crate::assumptions::exact_queries(||crate::assumptions::is_zero(expression))==Some(true){
        return Some(MpComplex::real(BigFloat::zero()));
    }
    let options=OPTIONS.with(|cell|cell.get());
    let budget=options.budget.unwrap_or(options.cancellation.unwrap_or(DEFAULT_CANCELLATION_BUDGET));
    let ceiling=precision.saturating_add(budget);
    // Absolute rounding error of the computed sum, as a binary exponent.
    let error_bound=|largest:i128,working:u32|largest.saturating_add(3+(usize::BITS-terms.len().leading_zeros()) as i128)-working as i128;
    // An explicit budget follows SymPy's schedule: terms start at ten guard
    // bits, each restart adds what the accuracy estimate lacks, and the loop
    // stops once the working precision exceeds the target by
    // min(budget, 2*working) bits, the working precision of the terms
    // excluding their ten guard bits. The expression evaluated as a whole is
    // computed at prec+20 bits (see eval_mp) where SymPy's evalf_add targets
    // prec+4 bits; its schedule starts from that target.
    let explicit=options.budget.is_some();
    let root=RESIDUE_ROOT.with(|cell|cell.get())==expression.dag_id();
    let precision=if explicit&&root {precision.saturating_sub(16).max(1)} else {precision};
    let mut working = if explicit {precision+10} else {precision + 32};
    let mut restarts=0u32;
    loop {
        let mut sum = MpComplex::real(number(constant, working));
        let mut largest_real = sum.re.mag();
        let mut largest_imaginary=sum.im.mag();
        for (term, coefficient) in terms {
            let value = ev(term, working)?;
            let value = if coefficient.is_one() { value } else { value.scale(&number(coefficient, working), working) };
            largest_real=largest_real.max(value.re.mag());
            largest_imaginary=largest_imaginary.max(value.im.mag());
            sum = sum.add(&value, working);
        }
        if known_real{sum.im=BigFloat::zero();}
        // Under chop, a component known to lie below the cutoff is zero, as
        // SymPy's chop_parts removes it at each level of evaluation.
        let mut real_chopped=false;let mut imaginary_chopped=false;
        if let Some(cutoff)=options.chop {
            if largest_real>i128::MIN/8 && (sum.re.is_zero()||sum.re.mag()<cutoff) && error_bound(largest_real,working)<cutoff {
                sum.re=BigFloat::zero();real_chopped=true;
            }
            if largest_imaginary>i128::MIN/8 && (sum.im.is_zero()||sum.im.mag()<cutoff) && error_bound(largest_imaginary,working)<cutoff {
                sum.im=BigFloat::zero();imaginary_chopped=true;
            }
        }
        // A large real component must not hide cancellation of a tiny
        // imaginary component (or vice versa), especially near a branch cut.
        if let Some(budget)=options.budget {
            let real_part=!real_chopped&&largest_real>i128::MIN/8;
            let imaginary_part=!known_real&&!imaginary_chopped&&largest_imaginary>i128::MIN/8;
            // As SymPy's complex_accuracy, the error of either part is
            // measured against the larger part; an unresolved zero alone
            // has accuracy -1.
            let resolved_re=real_part&&!sum.re.is_zero();
            let resolved_im=imaginary_part&&!sum.im.is_zero();
            let accuracy=if !real_part&&!imaginary_part {i128::MAX}
            else if !resolved_re&&!resolved_im {-1}
            else {
                let size=|value:&BigFloat,part:bool|if part&&!value.is_zero(){value.mag()}else{i128::MIN/4};
                let error=|part:bool,largest:i128|if part{error_bound(largest,working)}else{i128::MIN/4};
                let largest_size=size(&sum.re,real_part).max(size(&sum.im,imaginary_part));
                largest_size.saturating_sub(error(real_part,largest_real).max(error(imaginary_part,largest_imaginary)))
            };
            if accuracy>=precision as i128 {return Some(sum);}
            // SymPy's evalf_add result once the budget is exhausted: each
            // part with its own accuracy, an exact zero as its error bound.
            let residue=|sum:MpComplex|{
                let part=|value:BigFloat,present:bool,largest:i128|{
                    if !present {return (value,i128::MAX);}
                    let error=error_bound(largest,working);
                    // Without a correct bit, the part is reported as SymPy's
                    // scaled zero: the error bound itself, with accuracy -1.
                    let accuracy=if value.is_zero(){-1}else{value.mag().saturating_sub(error)};
                    if accuracy<1 {(BigFloat::one().ldexp(error),-1)} else {(value,accuracy)}
                };
                let (re,real_accuracy)=part(sum.re,real_part,largest_real);
                let (im,imaginary_accuracy)=part(sum.im,imaginary_part,largest_imaginary);
                RESIDUE.with(|cell|cell.set(Some((real_accuracy,imaginary_accuracy))));
                MpComplex{re,im}
            };
            let root=options.residue&&root;
            let terms_precision=working-10;
            if terms_precision.saturating_sub(precision)>budget.min(terms_precision.saturating_mul(2)) {
                if root {return Some(residue(sum));}
                if options.partial&&accuracy>=1 {
                    let achieved=u32::try_from(accuracy).unwrap_or(u32::MAX);
                    ACHIEVED.with(|cell|cell.set(Some(cell.get().map_or(achieved,|current|current.min(achieved)))));
                    return Some(sum);
                }
                EXHAUSTED.with(|cell|{let mut cell=cell.borrow_mut();if cell.is_none(){*cell=Some(expression.clone());}});
                return None;
            }
            let step=(precision as i128).saturating_sub(accuracy).max(10+(1i128<<restarts.min(60)));
            let next=i128::from(working).saturating_add(step);
            if next>i128::from(EXPLICIT_WORKING_LIMIT) {
                if root {return Some(residue(sum));}
                EXHAUSTED.with(|cell|{let mut cell=cell.borrow_mut();if cell.is_none(){*cell=Some(expression.clone());}});
                return None;
            }
            working=next as u32;
            restarts+=1;
            continue;
        }
        let real_loss=if real_chopped{0}else{largest_real.saturating_sub(sum.re.mag()).max(0)};
        let imaginary_loss=if known_real||imaginary_chopped{0}else{largest_imaginary.saturating_sub(sum.im.mag()).max(0)};
        let cancellation = u64::try_from(real_loss.max(imaginary_loss)).unwrap_or(u64::MAX/4).min(u64::MAX/4);
        let required = precision as u64 + cancellation + 24;
        if required <= working as u64 { return Some(sum); }
        if working >= ceiling {
            // The default budget is exhausted: report the residue with
            // SymPy's default schedule (maxprec = max(prec, 332) bits).
            if options.residue&&root {
                let target=precision.saturating_sub(20).max(1);
                let _scope=OptionsScope::enter(EvalfOptions{budget:Some(target.max(332)),..options});
                return eval_sum(expression,constant,terms,precision);
            }
            EXHAUSTED.with(|cell|{let mut cell=cell.borrow_mut();if cell.is_none(){*cell=Some(expression.clone());}});
            return None;
        }
        working = if (sum.re.is_zero()&&!real_chopped&&largest_real>i128::MIN/8)
            ||(!known_real&&!imaginary_chopped&&sum.im.is_zero()&&largest_imaginary>i128::MIN/8) {
            working.saturating_add((working - precision).max(32)).min(ceiling)
        } else {
            (required.min(ceiling as u64)) as u32
        };
    }
}

fn pow_of(b: &Expr, x: &Expr, p: u32) -> Option<MpComplex> {
    // The relative error of b**x is |x| times that of b, and |x ln b| times
    // the relative error of x: both are evaluated with those extra bits.
    let base = |q: u32| -> Option<MpComplex> {
        let mut bv = ev(b, q)?;
        // A structurally certified real expression must not acquire a different
        // square-root branch from rounding noise in a cancelling imaginary part.
        // Genuine tiny imaginary inputs are retained: this is not numeric chop.
        if !bv.is_real()&&crate::assumptions::exact_queries(||crate::assumptions::is_real(b))==Some(true){
            bv.im=BigFloat::zero();
        }
        Some(bv)
    };
    let bits = |k: i64| 64 - k.unsigned_abs().leading_zeros();
    if let Some(k) = x.as_i64() {
        return base(p + bits(k))?.powi(k, p);
    }
    if let Some((numerator,denominator))=x.as_number().and_then(Number::numer_denom){
        if denominator==2.into(){
            if let Some(n)=numerator.to_i64(){return base(p + bits(n))?.sqrt(p+16)?.powi(n,p);}
        }
    }
    // mpmath's mpc_pow_mpf: an integer or half-integer float exponent is an
    // integer power, of the square root for a half-integer.
    if let Some(n) = x.as_number().filter(|n| n.is_float() && !n.is_zero()) {
        let twice = n.to_bigfloat(n.float_precision()?).ldexp(1);
        let t = twice.to_int_trunc();
        if BigFloat::from_int(&t) == twice {
            if let Some(t) = t.to_i64() {
                let bv = base(p + bits(t))?;
                return if t % 2 == 0 { bv.powi(t / 2, p) } else { bv.sqrt(p + 16)?.powi(t, p) };
            }
        }
    }
    let mut xv = ev(x, p)?;
    let size = |v: &MpComplex| [&v.re, &v.im].into_iter().filter(|c| !c.is_zero()).map(BigFloat::mag).max();
    let extra = |m: Option<i128>| m.unwrap_or(0).clamp(0, 1 << 20) as u32;
    let bv = base(p + extra(size(&xv)))?;
    if let Some(l) = bv.ln(32) {
        let bits = extra(size(&xv).zip(size(&l)).map(|(a, b)| a + b));
        if bits > 0 { xv = ev(x, p + bits)?; }
    }
    bv.pow(&xv, p)
}

/// One argument of a host special-function evaluation: a number, or a
/// (possibly nested) parameter tuple.
#[derive(Clone, Debug)]
pub enum HostArg {
    Number(MpComplex),
    Tuple(Vec<HostArg>),
}

/// An exception raised by the host evaluator, kept for the caller that
/// requested errors (SymPy propagates mpmath's errors from `evalf`).
pub struct HostError(pub Box<dyn std::any::Any + Send>);

/// Host evaluator of a special function head at a binary precision. It
/// returns `Ok(None)` when the value is unavailable; it must never return a
/// value with fewer correct bits than requested.
pub type SpecialHost = fn(&str, &[HostArg], u32) -> Result<Option<MpComplex>, HostError>;

static SPECIAL_HOST: std::sync::OnceLock<SpecialHost> = std::sync::OnceLock::new();

thread_local! {
    static HOST_ERRORS: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static HOST_ERROR: std::cell::RefCell<Option<HostError>> = const { std::cell::RefCell::new(None) };
}

/// Install the multiprecision evaluator used for special functions without a
/// native kernel (the Python bindings install mpmath, as SymPy uses it).
pub fn register_special_host(host: SpecialHost) {
    let _ = SPECIAL_HOST.set(host);
}

/// Is a host evaluator installed?
pub fn has_special_host() -> bool {
    SPECIAL_HOST.get().is_some()
}

/// While alive, the first host error is kept for [`take_host_error`];
/// a nested disabled scope (sign and magnitude decisions) discards them,
/// as SymPy's assumption queries catch evaluation errors.
pub struct HostErrorScope(bool);
impl HostErrorScope {
    pub fn enter(enabled: bool) -> Self {
        let previous = HOST_ERRORS.with(|cell| cell.replace(enabled));
        if enabled && !previous {
            HOST_ERROR.with(|cell| cell.borrow_mut().take());
        }
        Self(previous)
    }
}
impl Drop for HostErrorScope {
    fn drop(&mut self) {
        HOST_ERRORS.with(|cell| cell.set(self.0));
    }
}

/// The first host error raised inside the current enabled scope.
pub fn take_host_error() -> Option<HostError> {
    HOST_ERROR.with(|cell| cell.borrow_mut().take())
}

fn call_host(host: SpecialHost, name: &str, args: &[HostArg], p: u32, report: bool) -> Option<MpComplex> {
    match host(name, args, p) {
        Ok(value) => value,
        Err(error) => {
            if report && HOST_ERRORS.with(|cell| cell.get()) {
                HOST_ERROR.with(|cell| { cell.borrow_mut().get_or_insert(error); });
            }
            None
        }
    }
}

// Argument guard bits; SymPy converts arguments at prec + 5, and polygamma
// works at prec + 12.
const HOST_ARGUMENT_GUARD: u32 = 12;

fn host_arg(e: &Expr, p: u32) -> Option<HostArg> {
    if let Some(items) = crate::special::tuple_items(e) {
        return items.iter().map(|item| host_arg(item, p)).collect::<Option<Vec<_>>>().map(HostArg::Tuple);
    }
    ev(e, p).map(HostArg::Number)
}

fn host_eval(name: &str, args: &[Expr], p: u32) -> Option<MpComplex> {
    let host = *SPECIAL_HOST.get()?;
    refined(p, |q, first| {
        let values = args.iter().map(|arg| host_arg(arg, q + HOST_ARGUMENT_GUARD)).collect::<Option<Vec<_>>>()?;
        call_host(host, name, &values, q, first)
    })
}

/// Do two values agree to `p` bits, relative to the larger part of `b`?
fn agrees(a: &MpComplex, b: &MpComplex, p: u32) -> bool {
    let size = [&b.re, &b.im].into_iter().filter(|part| !part.is_zero()).map(BigFloat::mag).max();
    let Some(size) = size else { return a.re.is_zero() && a.im.is_zero() };
    let wp = p + 16;
    [a.re.sub(&b.re, wp), a.im.sub(&b.im, wp)].iter().all(|d| d.is_zero() || d.mag() < size - p as i128 - 1)
}

// SymPy evaluates a special function once at the target precision; the
// result then carries the conditioning error of the rounded arguments and
// any loss inside mpmath. Evaluating again with more argument and working
// bits until two values agree certifies the returned bits. Values that keep
// disagreeing are not returned; a failure above the first precision returns
// the first value, as SymPy does.
fn refined(p: u32, evaluate: impl Fn(u32, bool) -> Option<MpComplex>) -> Option<MpComplex> {
    let limit = p.saturating_mul(2).saturating_add(512);
    let mut extra = 10u32;
    let mut previous = evaluate(p + extra, true)?;
    loop {
        extra = extra * 2 + 10;
        if extra > limit {
            return None;
        }
        let Some(current) = evaluate(p + extra, false) else {
            return Some(previous.round_to(p));
        };
        if agrees(&previous, &current, p) {
            return Some(current.round_to(p));
        }
        previous = current;
    }
}

/// Whether `exp_polar(w)` is on the principal sheet, `im(w)` in `(-pi, pi]`,
/// where SymPy evaluates it as `exp(w)`. `None` when undecided.
fn exp_polar_principal(w: &Expr, p: u32) -> Option<bool> {
    let ratio = crate::ops::div(w, &crate::ops::mul2(&Expr::i(), &Expr::pi()));
    if let Some(n) = ratio.as_number().filter(|n| n.is_rational()) {
        return Some(n.add(&Number::int(1)).is_positive() && !n.sub(&Number::int(1)).is_positive());
    }
    let wp = p + 20;
    let value = ev(w, wp)?;
    let half_turn = pi(wp);
    let upper = value.im.sub(&half_turn, wp);
    let lower = value.im.add(&half_turn, wp);
    let unresolved = |d: &BigFloat| d.is_zero() || d.mag() < half_turn.mag() - wp as i128 + 8;
    if unresolved(&upper) || unresolved(&lower) {
        return None;
    }
    Some(upper.sign() < 0 && lower.sign() > 0)
}

/// SymPy's `meijerg._eval_evalf`: one non-principal `exp_polar` factor of
/// the argument is passed separately so that the host can evaluate on its
/// sheet through mpmath's `r` parameter.
fn meijerg_host(args: &[Expr], p: u32) -> Option<MpComplex> {
    let host = *SPECIAL_HOST.get()?;
    let [a, b, z] = args else { return None };
    let polar_exponent = |e: &Expr| match e.kind() {
        Kind::Func(Func::ExpPolar, w) if w.len() == 1 && exp_polar_principal(&w[0], p) == Some(false) => Some(Some(w[0].clone())),
        Kind::Func(Func::ExpPolar, w) if w.len() == 1 && exp_polar_principal(&w[0], p).is_none() => None,
        _ => Some(None),
    };
    let (rest, branch) = match z.kind() {
        Kind::Mul(coefficient, factors) => {
            let mut branch = None;
            let mut rest = vec![Expr::num(coefficient.clone())];
            for (base, exponent) in factors {
                let factor = crate::ops::pow(base, exponent);
                match polar_exponent(&factor)? {
                    Some(w) if branch.is_none() => branch = Some(w),
                    Some(_) => return None,
                    None => rest.push(factor),
                }
            }
            (crate::ops::mul(&rest), branch)
        }
        _ => match polar_exponent(z)? {
            Some(w) => (Expr::one(), Some(w)),
            None => (z.clone(), None),
        },
    };
    refined(p, |q, first| {
        let wp = q + HOST_ARGUMENT_GUARD;
        let polar = match &branch {
            Some(w) => ev(w, wp)?,
            None => MpComplex::real(BigFloat::zero()),
        };
        let values = vec![host_arg(a, wp)?, host_arg(b, wp)?, HostArg::Number(ev(&rest, wp)?), HostArg::Number(polar)];
        call_host(host, "meijerg", &values, q, first)
    })
}

fn func(f: &Func, args: &[Expr], p: u32) -> Option<MpComplex> {
    let wp = p + 10;
    if let Some(elementary) = match f {
        Func::Special(s) => crate::elementary_mp::Elementary::of_special(*s),
        _ => crate::elementary_mp::Elementary::of_func(f),
    } {
        if args.len() != 1 { return None; }
        return crate::elementary_mp::evaluate(elementary, &args[0], p);
    }
    if let Func::Special(s) = f {
        if *s == crate::special::SpecialFn::MeijerG {
            return meijerg_host(args, p);
        }
        let value = crate::special::mp::eval(*s, args, p).or_else(|| host_eval(s.name(), args, p))?;
        // SymPy's Ei._eval_evalf: an argument on the ray polar_lift(-1)*(0, oo)
        // lies on the upper side of the cut, Ei(-x) + I*pi.
        if *s == crate::special::SpecialFn::Ei && args.len() == 1 && crate::special::polar::has_polar_func(&args[0]) {
            let lifted = crate::functions::func(Func::ExpPolar, vec![crate::ops::mul2(&Expr::i(), &Expr::pi())]);
            if crate::assumptions::is_positive(&crate::ops::div(&args[0], &lifted)) == Some(true) {
                return Some(MpComplex { re: value.re, im: value.im.add(&pi(wp), wp) }.round_to(p));
            }
        }
        return Some(value);
    }
    if let Func::ExpPolar = f {
        let w = args.first()?;
        if !exp_polar_principal(w, p)? {
            return None;
        }
        let value = ev(w, wp)?.exp(p)?;
        // On the boundary im(w) = pi SymPy keeps the real part.
        let ratio = crate::ops::div(w, &crate::ops::mul2(&Expr::i(), &Expr::pi()));
        return Some(if ratio.as_number().is_some_and(|n| n.is_one()) { MpComplex::real(value.re) } else { value });
    }
    if matches!(f, Func::Undefined(name) if name.is_internal("UnevaluatedExpr")) {
        return ev(args.first()?, p);
    }
    if let Func::Atan2 = f {
        let (y, x) = (ev(&args[0], wp)?, ev(&args[1], wp)?);
        if !(y.is_real() && x.is_real()) {
            return None;
        }
        return Some(MpComplex::real(BigFloat::atan2(&y.re, &x.re, p)));
    }
    if matches!(f, Func::Heaviside | Func::DiracDelta) && !args.is_empty() {
        // An argument substituted by evalf(subs=...) is held: rebuild the
        // node so an exact zero takes the value H0 (SymPy evaluates
        // Heaviside(0, H0) to H0), then decide a numeric argument by sign.
        fn evaluated(e: &Expr) -> Expr {
            if !matches!(e.kind(), Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Func(..)) {
                return e.clone();
            }
            let new: Vec<Expr> = crate::simptools::args(e).iter().map(evaluated).collect();
            crate::simptools::rebuild(e, &new)
        }
        let rebuilt = crate::functions::func(f.clone(), args.iter().map(evaluated).collect());
        if !matches!(rebuilt.kind(), Kind::Func(g, _) if g == f) {
            return ev(&rebuilt, p);
        }
        let a = ev(&args[0], wp)?;
        if !a.is_real() || a.re.is_zero() {
            return None;
        }
        let one = *f == Func::Heaviside && a.re.sign() > 0;
        return Some(MpComplex::real(if one { BigFloat::from_i64(1) } else { BigFloat::zero() }));
    }
    let a = ev(args.first()?, wp)?;
    Some(match f {
        Func::Abs => MpComplex::real(a.abs(p)),
        Func::Re => MpComplex::real(a.re),
        Func::Im => MpComplex::real(a.im),
        Func::Conjugate => MpComplex { re: a.re, im: a.im.neg() },
        Func::Arg => MpComplex::real(a.arg(p)),
        Func::Sign if a.is_real() => MpComplex::real(BigFloat::from_i64(a.re.sign() as i64)),
        // sign(z) = z/|z| for a nonzero complex value, as SymPy's
        // sign.eval returns arg/Abs(arg) for a nonreal number.
        Func::Sign if !a.re.is_zero() || !a.im.is_zero() => a.div(&MpComplex::real(a.abs(wp)), p)?,
        Func::Floor | Func::Ceiling if a.is_real() => {
            let t = a.re.to_int_trunc();
            let tf = BigFloat::from_int(&t);
            let c = a.re.cmp(&tf);
            let v = match (f, c) {
                (Func::Floor, std::cmp::Ordering::Less) => t - 1,
                (Func::Ceiling, std::cmp::Ordering::Greater) => t + 1,
                _ => t,
            };
            MpComplex::real(BigFloat::from_int(&v))
        }
        Func::Erf if a.is_real() => match erf_mp(&a.re, p) {
            Some(value) => MpComplex::real(value),
            None => return host_eval("erf", args, p),
        },
        Func::Erfc if a.is_real() => match erfc_mp(&a.re, p) {
            Some(value) => MpComplex::real(value),
            None => return host_eval("erfc", args, p),
        },
        Func::Gamma if a.is_real() => match gamma_mp(&a.re, p) {
            Some(value) => MpComplex::real(value),
            None => return host_eval("gamma", args, p),
        },
        Func::Erf => return host_eval("erf", args, p),
        Func::Erfc => return host_eval("erfc", args, p),
        Func::Gamma => return host_eval("gamma", args, p),
        Func::LogGamma => return host_eval("loggamma", args, p),
        _ => return None,
    })
}

/// Euler–Mascheroni constant (Brent–McMillan would be faster; this uses
/// the simple exponential integral formula with fixed point arithmetic).
pub fn euler_gamma(prec: u32) -> BigFloat {
    // gamma = H - ln(n) - corrections; use Brent-McMillan B1:
    // with n ~ prec*ln2/4, A = sum_{k} (n^k/k!)^2 (H_k - ln n), B = sum (n^k/k!)^2
    let wp = prec + 40;
    let n = ((wp as f64) * std::f64::consts::LN_2 / 4.0).ceil() as u64 + 1;
    let one = BigInt::from(1) << (wp as usize);
    let ln_n = BigFloat::from_i64(n as i64).ln(wp).unwrap();
    let ln_n_fx = {
        let s = ln_n.exp + wp as i128;
        if s >= 0 {
            &ln_n.man << (s as usize)
        } else {
            &ln_n.man >> ((-s) as usize)
        }
    };
    let n2 = BigInt::from(n * n);
    let mut a = -&ln_n_fx;
    let mut b = one.clone();
    let mut u = a.clone();
    let mut v = one.clone();
    let mut k = 1u64;
    loop {
        let k2 = BigInt::from(k * k);
        u = (&u * &n2) / &k2;
        v = (&v * &n2) / &k2;
        u = &u + (&v / BigInt::from(k));
        if v.bits() == 0 && u.bits() == 0 {
            break;
        }
        a += &u;
        b += &v;
        k += 1;
        if k > 40 * n + 100 {
            break;
        }
    }
    let r = BigFloat::from_int(&a).div(&BigFloat::from_int(&b), prec).unwrap();
    let _ = ln2;
    r
}

/// Catalan's constant with Ramanujan's series
/// `G = pi/8 log(2 + sqrt(3)) + 3/8 sum_n 1/((2n+1)**2 binomial(2n, n))`.
pub fn catalan(prec: u32) -> BigFloat {
    let wp = prec + 30;
    let one = BigInt::from(1) << (wp as usize);
    // r = 1/binomial(2n, n) in fixed point
    let mut r = one.clone();
    let mut sum = BigInt::from(0);
    let mut n = 0u64;
    while r.bits() > 0 {
        let d = BigInt::from(2 * n + 1);
        sum += &r / (&d * &d);
        r = (&r * BigInt::from(n + 1)) / BigInt::from(2 * (2 * n + 1));
        n += 1;
    }
    let series = BigFloat::round(sum, -(wp as i128), wp).mul_i64(3, wp).ldexp(-3);
    let three = BigFloat::from_i64(3).sqrt(wp).unwrap();
    let lg = BigFloat::from_i64(2).add(&three, wp).ln(wp).unwrap();
    pi(wp).mul(&lg, wp).ldexp(-3).add(&series, wp).round_to(prec)
}

/// `GoldenRatio = (1 + sqrt(5))/2`.
pub fn golden_ratio(prec: u32) -> BigFloat {
    let wp = prec + 10;
    BigFloat::from_i64(5).sqrt(wp).unwrap().add(&BigFloat::one(), wp).ldexp(-1).round_to(prec)
}

/// Real root of `x**3 - x**2 - x - 1` (Newton iteration).
pub fn tribonacci_constant(prec: u32) -> BigFloat {
    let wp = prec + 20;
    let mut x = BigFloat::from_f64(1.839_286_755_214_161);
    for _ in 0..64 {
        let x2 = x.mul(&x, wp);
        let x3 = x2.mul(&x, wp);
        let f = x3.sub(&x2, wp).sub(&x, wp).sub(&BigFloat::one(), wp);
        let df = x2.mul_i64(3, wp).sub(&x.mul_i64(2, wp), wp).sub(&BigFloat::one(), wp);
        let step = f.div(&df, wp).unwrap();
        x = x.sub(&step, wp);
        if step.is_zero() || step.mag() < -(wp as i128) {
            break;
        }
    }
    x.round_to(prec)
}

// The power series of erf cancels about x**2*log2(e) bits; beyond this
// square the native kernel defers to the host evaluator.
const ERF_SERIES_LIMIT: f64 = 4096.0;

fn erf_mp(x: &BigFloat, prec: u32) -> Option<BigFloat> {
    // erf x = 2/sqrt(pi) * sum (-1)^n x^(2n+1) / (n! (2n+1)); fine for
    // moderate |x| with enough guard bits (cancellation ~ x^2 / ln 2 bits)
    let x2f = x.to_f64() * x.to_f64();
    // 1 - |erf x| < exp(-x**2): the value rounds to +-1 at this precision.
    if x2f * std::f64::consts::LOG2_E > prec as f64 + 8.0 {
        return Some(BigFloat::from_i64(x.sign() as i64));
    }
    if x2f > ERF_SERIES_LIMIT {
        return None;
    }
    let wp = prec + 20 + (x2f * 1.5) as u32;
    let x2 = x.mul(x, wp);
    let mut term = x.clone();
    let mut sum = x.clone();
    let mut n = 1i64;
    loop {
        term = term.mul(&x2, wp).neg().div_i64(n, wp);
        let t = term.div_i64(2 * n + 1, wp);
        if t.is_zero() || t.mag() < sum.mag() - wp as i128 - 4 {
            break;
        }
        sum = sum.add(&t, wp);
        n += 1;
    }
    let c = BigFloat::from_i64(2).div(&pi(wp).sqrt(wp)?, wp)?;
    Some(sum.mul(&c, prec))
}

/// erfc x, accurate to `prec` bits: `1 - erf x` loses about x**2*log2(e)
/// bits for positive x, which are added to the working precision.
fn erfc_mp(x: &BigFloat, prec: u32) -> Option<BigFloat> {
    let x2f = x.to_f64() * x.to_f64();
    if x.sign() <= 0 || x2f <= 1.0 {
        return Some(BigFloat::one().sub(&erf_mp(x, prec + 10)?, prec));
    }
    if x2f > ERF_SERIES_LIMIT {
        return None;
    }
    let loss = (x2f * std::f64::consts::LOG2_E).ceil() as u32;
    let wp = prec + loss + 16;
    let erf = erf_mp(x, wp)?;
    let value = BigFloat::one().sub(&erf, wp);
    if value.is_zero() {
        return None;
    }
    Some(value.round_to(prec))
}

fn gamma_mp(x: &BigFloat, prec: u32) -> Option<BigFloat> {
    // integers: factorial; otherwise Spouge's approximation with a = prec*ln2/ln(2pi)
    let t = x.to_int_trunc();
    if BigFloat::from_int(&t) == *x {
        let n = num_traits::ToPrimitive::to_i64(&t)?;
        if n <= 0 {
            return None;
        }
        let mut acc = BigInt::from(1);
        for i in 2..n {
            acc *= i;
        }
        return Some(BigFloat::from_int(&acc).round_to(prec));
    }
    let wp = prec.saturating_mul(2) + 80;
    if x.sign() < 0 || x.cmp(&BigFloat { man: BigInt::from(1), exp: -1 }) == std::cmp::Ordering::Less {
        // reflection: gamma(x) = pi / (sin(pi x) gamma(1-x))
        let one_minus = BigFloat::one().sub(x, wp);
        let g = gamma_mp(&one_minus, wp)?;
        let p = pi(wp);
        let s = p.mul(x, wp).sin(wp);
        return p.div(&s.mul(&g, wp), prec);
    }
    let a = (((prec + 40) as f64) * std::f64::consts::LN_2 / (2.0 * std::f64::consts::PI).ln()).ceil() as i64;
    let z = x.sub(&BigFloat::one(), wp);
    let af = BigFloat::from_i64(a);
    // c0 = sqrt(2 pi)
    let mut sum = pi(wp).ldexp(1).sqrt(wp)?;
    let mut fact = BigFloat::one(); // (k-1)!
    for k in 1..a {
        let kf = BigFloat::from_i64(k);
        // c_k = (-1)^(k-1) (a-k)^(k-1/2) e^(a-k) / (k-1)!
        let amk = BigFloat::from_i64(a - k);
        let e1 = amk.powf(&BigFloat::from_i64(2 * k - 1).ldexp(-1), wp)?;
        let e2 = amk.checked_exp(wp)?;
        let mut ck = e1.mul(&e2, wp).div(&fact, wp)?;
        if k % 2 == 0 {
            ck = ck.neg();
        }
        sum = sum.add(&ck.div(&z.add(&kf, wp), wp)?, wp);
        fact = fact.mul(&kf, wp);
    }
    // gamma(z+1) = (z+a)^(z+1/2) e^(-(z+a)) * sum
    let za = z.add(&af, wp);
    let half = BigFloat { man: BigInt::from(1), exp: -1 };
    let pw = za.powf(&z.add(&half, wp), wp)?;
    let ex = za.neg().checked_exp(wp)?;
    Some(pw.mul(&ex, wp).mul(&sum, prec))
}

#[cfg(test)]
mod tests {
    #[test]
    fn shared_dag_evaluation_is_linear_and_precision_specific() {
        use super::*;
        let mut expression=Expr::one();
        for _ in 0..80 {
            expression=Expr::from_unevaluated_kind(Kind::Mul(Number::Int(1),vec![
                (expression.clone(),Expr::one()),(expression,Expr::one())]));
        }
        EVALUATIONS.with(|count|count.set(0));
        assert_eq!(eval_mp(&expression,100).unwrap().re,BigFloat::one());
        assert!(EVALUATIONS.with(|count|count.get())<=80);

        let _scope=NumericScope::enter();
        let sine=crate::functions::func(Func::Sin,vec![Expr::one()]);
        let low=ev(&sine,53).unwrap();
        let high=ev(&sine,200).unwrap();
        assert_ne!(low.re,high.re);
        assert_eq!(high.re,ev_uncached(&sine,200).unwrap().re);
    }

    #[test]
    fn assumption_refinements_share_precision_specific_descendants() {
        let mut value=Expr::one();
        for _ in 0..80 {
            value=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![
                (value.clone(),Number::ONE),(value,Number::MINUS_ONE)]));
        }
        EVALUATIONS.with(|count|count.set(0));
        // The reference also decides that this cancelling sum is not positive.
        assert_eq!(crate::assumptions::is_positive(&value),Some(false));
        let visits=EVALUATIONS.with(|count|count.get());
        assert!(visits>0 && visits<=4*80,"{visits} evaluations");
        NUMERIC_MEMO.with(|memo|{let memo=memo.borrow();assert_eq!(memo.depth,0);assert!(memo.values.is_none());});
    }

    #[test]
    fn held_sum_does_not_collect_inexact_coefficients_at_machine_precision() {
        let term=Expr::pi();
        let expression=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![
            (term.clone(),Number::float(1e16)),(term.clone(),Number::ONE),
            (term,Number::float(-1e16))]));
        let actual=eval_mp(&expression,200).unwrap();
        assert!(actual.re.sub(&pi(220),220).abs().mag() < -190);
        let infinity=Expr::from_unevaluated_kind(Kind::Func(Func::Log,vec![Expr::zero()]));
        let undefined=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![
            (infinity.clone(),Number::ONE),(infinity,Number::MINUS_ONE)]));
        assert!(eval_mp(&undefined,200).is_none());
    }
    use super::*;
    use crate::bigfloat::dps_to_prec;
    use crate::prelude::*;

    #[test]
    fn numerical_symbol_substitution_preserves_cancellation_and_sequential_values() {
        let x=Expr::symbol("x"); let y=Expr::symbol("y"); let z=Expr::symbol("z");
        let expression=add(&[x.clone(),y.clone(),mul(&[Expr::int(-1),z.clone()])]);
        for large in [Expr::num(Number::float(1e16)),
            Expr::num(Number::mp_float(BigFloat::from_decimal("1e100",70).unwrap(),70)),
            Expr::num(Number::mp_float(BigFloat::from_decimal("1e1000",70).unwrap(),70))] {
            let pairs=vec![(x.clone(),large.clone()),(y.clone(),Expr::one()),(z.clone(),large)];
            assert!(crate::subs::subs(&expression,&pairs).is_zero());
            let held=substitute_for_evalf(&expression,&pairs);
            assert_eq!(eval_mp(&held,200).unwrap().re,BigFloat::one(),"held {held}; kind {:?}",held.kind());
            assert_eq!(super::evalf(&held,200),to_expr(MpComplex::real(BigFloat::one()),200));
        }
        let sequential=substitute_for_evalf(&add(&[x.clone(),y.clone()]),&[(x,y.clone()),(y,Expr::int(2))]);
        assert_eq!(eval_mp(&sequential,200).unwrap().re,BigFloat::from_i64(4));
    }

    #[test]
    fn numerical_substitution_preserves_binding_and_algebraic_fallbacks() {
        let parse=|s| crate::parser::parse_str(s).unwrap();
        for expression in [parse("Integral(x,(x,0,1))"),parse("Derivative(f(x),x)"),parse("x**4-2"),parse("Piecewise((x,x>0),(0,True))")] {
            let pairs=if expression.is_pow() || expression.is_add() {vec![(parse("x**2"),Expr::int(3))]}
                else {vec![(parse("x"),Expr::int(2))]};
            assert_eq!(substitute_for_evalf(&expression,&pairs),crate::subs::subs(&expression,&pairs));
        }
    }

    #[test]
    fn complex_square_root_keeps_tiny_genuine_branch_components(){
        let precision=200;
        let tiny=BigFloat{man:BigInt::from(1),exp:-4000};
        for sign in [-1,1]{
            let z=MpComplex{re:BigFloat::from_i64(-1),im:if sign<0{tiny.neg()}else{tiny.clone()}};
            let root=z.sqrt(precision).unwrap();
            assert!(root.re.sign()>0);
            assert_eq!(root.im.sign(),sign);
            let recovered=root.mul(&root,precision);
            assert_eq!(recovered.re,BigFloat::from_i64(-1));
            assert!(recovered.im.sub(&z.im,precision).abs().mag()<tiny.mag()-precision as i128+4);
        }
    }

    #[test]
    fn chopping_uses_native_component_magnitudes_without_machine_underflow() {
        let precision=dps_to_prec(80);
        let tiny=Expr::num(Number::mp_float(BigFloat::from_decimal("1e-1000",precision).unwrap(),precision));
        assert!(evalf_chop(&tiny,precision,-(precision as i128)).is_zero());
        let expression=add(&[Expr::one(),mul(&[Expr::i(),tiny.clone()])]);
        assert_eq!(evalf_chop(&expression,precision,-(precision as i128)),super::evalf(&Expr::one(),precision));
        let symbolic=add(&[Expr::symbol("x"),tiny]);
        assert_eq!(evalf_chop(&symbolic,precision,-(precision as i128)),super::evalf(&symbolic,precision));
    }

    fn s(e: &Expr, dps: u32) -> String {
        let v = eval_mp(e, dps_to_prec(dps)).unwrap();
        v.re.to_string_dps(dps, false)
    }

    #[test]
    fn eval() {
        let e = add(&[Expr::pi(), Expr::e()]);
        assert_eq!(s(&e, 30), "5.85987448204883847382293085463");
        assert_eq!(s(&sqrt(&Expr::int(2)), 40), "1.414213562373095048801688724209698078570");
        assert_eq!(s(&exp(&Expr::pi()), 25), "23.14069263277926900572909");
        assert_eq!(s(&Expr::constant(Constant::EulerGamma), 30), "0.577215664901532860606512090082");
        assert_eq!(s(&gamma(&Expr::rational(1, 3)), 25), "2.678938534707747633655693");
        assert_eq!(s(&erf(&Expr::one()), 25), "0.8427007929497148693412206");
    }

    #[test]
    fn precision_survives_native_expression_arithmetic_and_serialization() {
        let p = dps_to_prec(70);
        let value = crate::bigfloat::BigFloat::from_decimal("0.123456789012345678901234567890123456789", p).unwrap();
        let original = Expr::num(Number::mp_float(value.clone(), p));
        assert_eq!(max_precision(&original), Some(p));
        assert_eq!(crate::parser::parse_str(&crate::printer::srepr(&original)).unwrap(), original);
        let increment = crate::parser::parse_str("Float('1e-60', 70)").unwrap();
        let changed = add(&[original.clone(), increment]);
        assert_ne!(changed, original);
        let delta = sub(&changed, &original);
        let result = eval_mp(&delta, p).unwrap();
        assert!((result.re.to_f64() / 1e-60 - 1.0).abs() < 1e-10);
        let neg = Number::mp_float(value.neg(), p);
        assert_eq!(neg.add(original.as_number().unwrap()), Number::ZERO);
    }

    #[test]
    fn evalf_retains_partial_symbolic_precision() {
        let x = Expr::symbol("x");
        let p = dps_to_prec(80);
        let expression = add(&[x.clone(), Expr::pi()]);
        let evaluated = super::evalf(&expression, p);
        assert_eq!(crate::traverse::free_symbols(&evaluated), vec![x]);
        assert_eq!(max_precision(&evaluated), Some(p));
        assert!(crate::printer::to_str(&evaluated).contains("3.141592653589793238462643383279502884197"));
        let raw = Expr::from_kind(Kind::Func(Func::Undefined("unknown_numeric".into()), vec![Expr::pi()]));
        assert_eq!(super::evalf(&raw, p), raw);
    }

    #[test]
    fn lambda_evalf_changes_its_body_and_preserves_bound_signature() {
        let x=Expr::symbol("x");let y=Expr::symbol("y");
        let signature=crate::special::tuple(vec![x.clone(),y.clone()]);
        let body=add(&[crate::ops::mul2(&Expr::rational(1,3),&x),crate::ops::mul2(&Expr::rational(1,7),&y)]);
        let original=crate::sets::lambda(&signature,&body);
        let precision=dps_to_prec(60);
        let result=super::evalf(&original,precision);
        let (variables,converted)=crate::sets::lambda_parts(&result).unwrap();
        assert_eq!(variables,vec![x,y]);
        assert_eq!(max_precision(&converted),Some(precision));
        assert_eq!(crate::traverse::free_symbols(&result),Vec::<Expr>::new());
    }

    #[test]
    fn cancellation_increases_working_precision() {
        let expression = crate::parser::parse_str("sqrt(10**100 + 1) - 10**50").unwrap();
        let result = eval_mp(&expression, dps_to_prec(60)).unwrap();
        assert!(result.re.sign() > 0);
        assert!(result.re.to_string_dps(50, false).starts_with("5.0000000000000000000000000000000000000000000000000e-51"));
    }

    #[test]
    fn extreme_exponential_inputs_remain_symbolic() {
        let p = dps_to_prec(60);
        let huge = crate::parser::parse_str("Float('1e100', 60)").unwrap();
        for function in [Func::Exp, Func::Sinh, Func::Cosh] {
            for input in [huge.clone(), crate::ops::neg(&huge)] {
                let expression = crate::functions::func(function.clone(), vec![input]);
                // sinh(-x) is -sinh(x) by its eval rule, which runs first.
                assert!(expression.as_number().is_none());
                assert_eq!(super::evalf(&expression, p), expression);
            }
        }
        // tanh saturates to +-1 at any precision, as mpmath returns.
        for (input, sign) in [(huge.clone(), 1), (crate::ops::neg(&huge), -1)] {
            let value = eval_mp(&crate::functions::func(Func::Tanh, vec![input]), p).unwrap();
            assert_eq!(value.re.cmp(&BigFloat::from_i64(sign)), std::cmp::Ordering::Equal);
        }
        // mpmath exponents are unbounded; binary exponents up to
        // MAX_EXPONENT are materialized, so exp(2**63) evaluates as in the
        // reference, while exp(2**110) stays symbolic.
        let exp = |n: Expr| Expr::from_kind(Kind::Func(Func::Exp, vec![n]));
        let power = |k: u32| Expr::num(Number::from_bigint(num_bigint::BigInt::from(1) << k as usize));
        assert!(eval_mp(&exp(Expr::int(2_303_000)), p).is_some());
        let value = eval_mp(&exp(Expr::int(i64::MAX)), p).unwrap();
        assert_eq!(value.re.to_string_dps(19, false), "7.375865886842630054e+4005659580146785380");
        assert!(eval_mp(&exp(crate::ops::neg(&power(63))), p).is_some_and(|v| v.re.sign() > 0));
        assert!(eval_mp(&exp(power(110)), p).is_none());
        assert_eq!(super::evalf(&exp(power(110)), p), exp(power(110)));
    }
}
