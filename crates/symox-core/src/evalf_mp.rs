//! Arbitrary precision numerical evaluation of expressions.

use num_bigint::BigInt;
use num_traits::ToPrimitive;

use crate::bigfloat::{ln2, pi, BigFloat};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;

/// A complex number with arbitrary precision parts.
#[derive(Clone, Debug)]
pub struct MpComplex {
    pub re: BigFloat,
    pub im: BigFloat,
}

// Precision is part of the key: cancellation refinement can revisit the same
// physical node at several working precisions. Owners prevent address reuse.
struct NumericMemo {
    depth: usize,
    values: Option<rustc_hash::FxHashMap<(usize,u32,bool),(Expr,Option<MpComplex>)>>,
}
thread_local! {
    static NUMERIC_MEMO:std::cell::RefCell<NumericMemo> = const {
        std::cell::RefCell::new(NumericMemo{depth:0,values:None})
    };
    #[cfg(test)] static EVALUATIONS:std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
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
        // Keep exponential results inside a decimal exponent magnitude of
        // one million, matching the native decimal constructor's bound.
        if self.re.abs().cmp(&BigFloat::from_i64(2_000_000)) == std::cmp::Ordering::Greater { return None; }
        let r = self.re.checked_exp(p + 10)?;
        if self.im.is_zero() {
            return Some(MpComplex::real(r.round_to(p)));
        }
        let (s, c) = self.im.sin_cos(p + 10);
        Some(MpComplex { re: r.mul(&c, p), im: r.mul(&s, p) })
    }

    pub(crate) fn ln(&self, p: u32) -> Option<MpComplex> {
        if self.is_real() && self.re.sign() > 0 {
            return Some(MpComplex::real(self.re.ln(p)?));
        }
        let a = self.abs(p + 10);
        Some(MpComplex { re: a.ln(p)?, im: self.arg(p) })
    }

    fn pow(&self, e: &MpComplex, p: u32) -> Option<MpComplex> {
        if self.re.is_zero() && self.im.is_zero() {
            return if e.re.sign() > 0 { Some(MpComplex::real(BigFloat::zero())) } else { None };
        }
        if self.is_real() && e.is_real() && self.re.sign() > 0 {
            return Some(MpComplex::real(self.re.powf(&e.re, p)?));
        }
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
    [&value.re, &value.im].iter().all(|v| v.is_zero() || v.mag().unsigned_abs() <= 3_321_930)
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
    if let Some((re,im))=crate::complex_parts::literal_parts(e){
        let component=|e:&Expr|match e.as_number()?{
            Number::MpFloat(value,_)=>value.checked_round_to(prec),
            n=>number(n,prec).checked_round_to(prec),
        };
        return Some(MpComplex{re:component(&re)?,im:component(&im)?});
    }
    if let Kind::Func(Func::MatrixDotProduct(a,b),_)=e.kind(){return eval_mp(&crate::matrixcontraction::dot_product(a,b,true).ok()?,prec);}
    if matches!(e.kind(),Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product,_)){
        let evaluated=crate::functions::doit(e,true);
        if evaluated==*e{return None;}
        return eval_mp(&evaluated,prec);
    }
    let wp = prec + 20;
    let r = ev(e, wp)?;
    if !supported_magnitude(&r) { return None; }
    Some(MpComplex { re: r.re.round_to(prec), im: r.im.round_to(prec) })
}

/// Convert a computed value back to native expression numbers, retaining the
/// declared precision. Exact zero components stay exact zero.
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
    if let Some(value) = eval_mp(e, prec) { return to_expr(value, prec); }
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
        _ => e.clone(),
    }
}

/// Remove absolutely small real/imaginary numerical components by binary
/// magnitude. Partial symbolic coefficients retain ordinary evalf semantics.
pub fn evalf_chop(e:&Expr,prec:u32,cutoff:i64)->Expr {
    let Some(mut value)=eval_mp(e,prec) else {return evalf(e,prec)};
    if !value.re.is_zero() && value.re.mag()<cutoff {value.re=BigFloat::zero();}
    if !value.im.is_zero() && value.im.mag()<cutoff {value.im=BigFloat::zero();}
    to_expr(value,prec)
}

pub(crate) fn ev(e: &Expr, p: u32) -> Option<MpComplex> {
    if e.is_atom(){return ev_uncached(e,p);}
    let _scope=NumericScope::enter();
    let key=(e.dag_id(),p,crate::assumptions::numeric_refinement_enabled());
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
        Kind::Func(f, args) => func(f, args, p)?,
        _ => return None,
    })
}

/// Increase working precision when large terms nearly cancel. A fixed guard
/// budget would incorrectly turn sqrt(10**100 + 1) - 10**50 into exact zero.
fn eval_sum(expression:&Expr,constant: &Number, terms: &[(Expr, Number)], precision: u32) -> Option<MpComplex> {
    let known_real=crate::assumptions::exact_queries(||crate::assumptions::is_real(expression))==Some(true);
    if known_real&&crate::assumptions::exact_queries(||crate::assumptions::is_zero(expression))==Some(true){
        return Some(MpComplex::real(BigFloat::zero()));
    }
    let mut working = precision + 32;
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
        // A large real component must not hide cancellation of a tiny
        // imaginary component (or vice versa), especially near a branch cut.
        let real_loss=largest_real.saturating_sub(sum.re.mag()).max(0);
        let imaginary_loss=if known_real{0}else{largest_imaginary.saturating_sub(sum.im.mag()).max(0)};
        let cancellation = real_loss.max(imaginary_loss) as u64;
        let required = precision as u64 + cancellation + 24;
        if required <= working as u64 { return Some(sum); }
        if working >= precision + 4096 { return None; }
        working = if (sum.re.is_zero()&&largest_real>i64::MIN/8)
            ||(!known_real&&sum.im.is_zero()&&largest_imaginary>i64::MIN/8) {
            (working + (working - precision).max(32)).min(precision + 4096)
        } else {
            (required.min((precision + 4096) as u64)) as u32
        };
    }
}

fn pow_of(b: &Expr, x: &Expr, p: u32) -> Option<MpComplex> {
    let mut bv = ev(b, p)?;
    // A structurally certified real expression must not acquire a different
    // square-root branch from rounding noise in a cancelling imaginary part.
    // Genuine tiny imaginary inputs are retained: this is not numeric chop.
    if !bv.is_real()&&crate::assumptions::exact_queries(||crate::assumptions::is_real(b))==Some(true){
        bv.im=BigFloat::zero();
    }
    if let Some(k) = x.as_i64() {
        return bv.powi(k, p);
    }
    if let Some((numerator,denominator))=x.as_number().and_then(Number::numer_denom){
        if denominator==2.into(){
            if let Some(n)=numerator.to_i64(){return bv.sqrt(p+16)?.powi(n,p);}
        }
    }
    let xv = ev(x, p)?;
    bv.pow(&xv, p)
}

fn func(f: &Func, args: &[Expr], p: u32) -> Option<MpComplex> {
    let wp = p + 10;
    if let Func::Special(s) = f {
        return crate::special::mp::eval(*s, args, p);
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
    let a = ev(args.first()?, wp)?;
    let real_only = |v: Option<BigFloat>| v.map(MpComplex::real);
    Some(match f {
        Func::Exp => a.exp(p)?,
        Func::Log => a.ln(p)?,
        Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc if a.is_real() => {
            let (s, c) = a.re.sin_cos(wp);
            MpComplex::real(match f {
                Func::Sin => s,
                Func::Cos => c,
                Func::Tan => s.div(&c, p)?,
                Func::Cot => c.div(&s, p)?,
                Func::Sec => BigFloat::one().div(&c, p)?,
                _ => BigFloat::one().div(&s, p)?,
            })
        }
        Func::Sin | Func::Cos | Func::Tan => {
            // complex: via exponentials
            let iz = MpComplex { re: a.im.neg(), im: a.re.clone() };
            let e1 = iz.exp(wp)?;
            let e2 = e1.inv(wp)?;
            let half = BigFloat { man: BigInt::from(1), exp: -1 };
            let cosv = e1.add(&e2, wp).scale(&half, wp);
            let diff = e1.add(&e2.scale(&BigFloat::from_i64(-1), wp), wp);
            // sin = (e1 - e2) / (2i)
            let sinv = MpComplex { re: diff.im.mul(&half, wp), im: diff.re.neg().mul(&half, wp) };
            match f {
                Func::Sin => sinv,
                Func::Cos => cosv,
                _ => sinv.mul(&cosv.inv(wp)?, p),
            }
        }
        Func::Asin if a.is_real() => real_only(a.re.asin(p))?,
        Func::Acos if a.is_real() => real_only(a.re.acos(p))?,
        Func::Atan if a.is_real() => MpComplex::real(a.re.atan(p)),
        Func::Acot if a.is_real() => MpComplex::real(BigFloat::one().div(&a.re, wp)?.atan(p)),
        Func::Sinh | Func::Cosh | Func::Tanh if a.is_real() => {
            if a.re.abs().cmp(&BigFloat::from_i64(2_000_000)) == std::cmp::Ordering::Greater { return None; }
            let e1 = a.re.checked_exp(wp)?;
            let e2 = BigFloat::one().div(&e1, wp)?;
            let sh = e1.sub(&e2, wp).ldexp(-1);
            let ch = e1.add(&e2, wp).ldexp(-1);
            MpComplex::real(match f {
                Func::Sinh => sh,
                Func::Cosh => ch,
                _ => sh.div(&ch, p)?,
            })
        }
        Func::Asinh if a.is_real() => {
            // ln(x + sqrt(x^2+1))
            let x = &a.re;
            let s = x.mul(x, wp).add(&BigFloat::one(), wp).sqrt(wp)?;
            MpComplex::real(x.add(&s, wp).ln(p)?)
        }
        Func::Acosh if a.is_real() => {
            let x = &a.re;
            let s = x.mul(x, wp).sub(&BigFloat::one(), wp).sqrt(wp)?;
            MpComplex::real(x.add(&s, wp).ln(p)?)
        }
        Func::Atanh if a.is_real() => {
            let x = &a.re;
            let one = BigFloat::one();
            let q = one.add(x, wp).div(&one.sub(x, wp), wp)?;
            MpComplex::real(q.ln(wp)?.ldexp(-1))
        }
        Func::Abs => MpComplex::real(a.abs(p)),
        Func::Re => MpComplex::real(a.re),
        Func::Im => MpComplex::real(a.im),
        Func::Conjugate => MpComplex { re: a.re, im: a.im.neg() },
        Func::Arg => MpComplex::real(a.arg(p)),
        Func::Sign if a.is_real() => MpComplex::real(BigFloat::from_i64(a.re.sign() as i64)),
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
        Func::Erf if a.is_real() => MpComplex::real(erf_mp(&a.re, p)?),
        Func::Erfc if a.is_real() => MpComplex::real(BigFloat::one().sub(&erf_mp(&a.re, wp)?, p)),
        Func::Gamma if a.is_real() => MpComplex::real(gamma_mp(&a.re, p)?),
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
        let s = ln_n.exp + wp as i64;
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
    let series = BigFloat::round(sum, -(wp as i64), wp).mul_i64(3, wp).ldexp(-3);
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
        if step.is_zero() || step.mag() < -(wp as i64) {
            break;
        }
    }
    x.round_to(prec)
}

fn erf_mp(x: &BigFloat, prec: u32) -> Option<BigFloat> {
    // erf x = 2/sqrt(pi) * sum (-1)^n x^(2n+1) / (n! (2n+1)); fine for
    // moderate |x| with enough guard bits (cancellation ~ x^2 / ln 2 bits)
    let x2f = x.to_f64() * x.to_f64();
    if x2f > 1e6 {
        return Some(BigFloat::from_i64(x.sign() as i64));
    }
    let wp = prec + 20 + (x2f * 1.5) as u32;
    let x2 = x.mul(x, wp);
    let mut term = x.clone();
    let mut sum = x.clone();
    let mut n = 1i64;
    loop {
        term = term.mul(&x2, wp).neg().div_i64(n, wp);
        let t = term.div_i64(2 * n + 1, wp);
        if t.is_zero() || t.mag() < sum.mag() - wp as i64 - 4 {
            break;
        }
        sum = sum.add(&t, wp);
        n += 1;
    }
    let c = BigFloat::from_i64(2).div(&pi(wp).sqrt(wp)?, wp)?;
    Some(sum.mul(&c, prec))
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
            assert!(recovered.im.sub(&z.im,precision).abs().mag()<tiny.mag()-precision as i64+4);
        }
    }

    #[test]
    fn chopping_uses_native_component_magnitudes_without_machine_underflow() {
        let precision=dps_to_prec(80);
        let tiny=Expr::num(Number::mp_float(BigFloat::from_decimal("1e-1000",precision).unwrap(),precision));
        assert!(evalf_chop(&tiny,precision,-(precision as i64)).is_zero());
        let expression=add(&[Expr::one(),mul(&[Expr::i(),tiny.clone()])]);
        assert_eq!(evalf_chop(&expression,precision,-(precision as i64)),super::evalf(&Expr::one(),precision));
        let symbolic=add(&[Expr::symbol("x"),tiny]);
        assert_eq!(evalf_chop(&symbolic,precision,-(precision as i64)),super::evalf(&symbolic,precision));
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
        for function in [Func::Exp, Func::Sinh, Func::Cosh, Func::Tanh] {
            for input in [huge.clone(), crate::ops::neg(&huge)] {
                let expression = crate::functions::func(function.clone(), vec![input]);
                assert!(matches!(expression.kind(), Kind::Func(..)));
                assert_eq!(super::evalf(&expression, p), expression);
            }
        }
        let exponent = crate::parser::parse_str("Float('1000001', 60)").unwrap();
        let expression = crate::ops::pow(&Expr::int(10), &exponent);
        assert!(matches!(expression.kind(), Kind::Pow(..)));
        assert_eq!(super::evalf(&expression, p), expression);
        let boundary = Expr::from_kind(Kind::Func(Func::Exp, vec![Expr::int(2_000_000)]));
        assert!(eval_mp(&boundary, p).is_some());
        let outside = Expr::from_kind(Kind::Func(Func::Exp, vec![Expr::int(2_000_001)]));
        assert!(eval_mp(&outside, p).is_none());
    }
}
