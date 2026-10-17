//! Three-valued assumption queries (`is_positive`, `is_real`, ...).
//!
//! The inference propagates facts from symbol assumptions through sums,
//! products, powers and elementary functions, like SymPy's old assumption
//! system (`expr.is_positive`). `None` means "unknown".

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::symbol::Assumptions as A;
#[path = "assumptions_exact.rs"]
mod exact;
pub(crate) use exact::quadratic_facts;

thread_local! {
    static NUMERIC_REFINEMENT: std::cell::Cell<bool> = const { std::cell::Cell::new(true) };
}

/// Numeric query policy is part of every memoization context.
#[inline]
pub(crate) fn numeric_refinement_enabled()->bool{NUMERIC_REFINEMENT.with(|enabled|enabled.get())}

// Thirty bits per policy: thirteen Boolean slots and a four-bit Signs field.
// Even/odd use the single parity slot; two Boolean slots remain reserved.
const BANK_BITS:u32=30;
const SIGNS_OFFSET:u32=26;
#[derive(Clone,Copy)]
enum Query{Real,Positive,Negative,Nonnegative,Nonpositive,Zero,Nonzero,Integer,Rational,Parity,Finite}

fn algebraic_proof(e:&Expr)->Option<bool>{
    fn power(base:&Expr,exponent:&Expr,cache:&mut rustc_hash::FxHashMap<usize,(Expr,Option<bool>)>)->bool{
        exponent.is_rational_number()&&walk(base,cache)==Some(true)
            &&(exponent.as_number().is_some_and(|n|!n.is_negative())||exact_queries(||is_zero(base))==Some(false))
    }
    fn walk(e:&Expr,cache:&mut rustc_hash::FxHashMap<usize,(Expr,Option<bool>)>)->Option<bool>{
        if let Some((_,answer))=cache.get(&e.dag_id()){return *answer;}
        let answer=match e.kind(){
            Kind::Sym(symbol)=>symbol.facts().get("algebraic"),
            Kind::Func(Func::Undefined(head),_) if !head.internal=>head.facts().get("algebraic"),
            Kind::Num(number) if !number.is_float()||number.is_zero()=>Some(true),
            Kind::Const(Constant::I)=>Some(true),
            Kind::Const(Constant::Pi|Constant::E|Constant::Infinity|Constant::NegInfinity|Constant::ComplexInfinity)=>Some(false),
            Kind::Pow(base,exponent) if power(base,exponent,cache)=>Some(true),
            Kind::Add(c,terms) if !c.is_float()&&terms.iter().all(|(term,c)|!c.is_float()&&walk(term,cache)==Some(true))=>Some(true),
            Kind::Mul(c,factors) if !c.is_float()&&factors.iter().all(|(base,exponent)|power(base,exponent,cache))=>Some(true),
            _=>None,
        };
        cache.insert(e.dag_id(),(e.clone(),answer));answer
    }
    walk(e,&mut rustc_hash::FxHashMap::default())
}

/// An exact declared fact, or a conservative native proof for a scalar value.
/// Additional Symbol facts never fall back to numerical refinement.
pub fn fact(e:&Expr,name:&str)->Option<bool>{
    if let Some(symbol)=e.as_sym(){return symbol.facts().get(name);}
    // Float zero has a known integer value, but SymPy deliberately keeps its
    // rationality/irrationality queries unknown, as for other approximate values.
    if name=="irrational"&&e.as_number().is_some_and(|n|n.is_float()&&n.is_zero()){return None;}
    if let Kind::Func(Func::Undefined(head),_)=e.kind(){
        if !head.internal{return if name=="commutative"{crate::decompose::commutativity(e)}else{head.facts().get(name)};}
    }
    if matches!(name,"algebraic"|"transcendental"){
        if let Some(algebraic)=algebraic_proof(e){
            if name=="algebraic"{return Some(algebraic);}
            if algebraic{return Some(false);}
        }
    }
    match name {
        "real"=>return is_real(e),"positive"=>return is_positive(e),"negative"=>return is_negative(e),
        "nonnegative"=>return is_nonnegative(e),"nonpositive"=>return is_nonpositive(e),
        "zero"=>return is_zero(e),"nonzero"=>return is_nonzero(e),"integer"=>return is_integer(e),
        "rational"=>return is_rational(e),"finite"=>return is_finite(e),"even"=>return is_even(e),"odd"=>return is_odd(e),
        "commutative"=>return crate::decompose::commutativity(e),_=>{},
    }
    let mut seed=Vec::new();
    match e.kind(){
        Kind::Num(n)=>{
            seed.extend([("real",true),("finite",true),("zero",n.is_zero()),("positive",n.is_positive()),("negative",n.is_negative())]);
            if !n.is_float(){seed.push(("rational",true));}
            if let Some(integer)=number_integer(n){seed.push(("integer",integer));}
            if let Some(integer)=n.to_bigint(){
                seed.push(("integer",true));seed.push(("even",num_integer::Integer::is_even(&integer)));
                // Native primality is deterministic on the machine range.
                // Larger probable primes must not become exact assumptions.
                if integer<=num_bigint::BigInt::from(u64::MAX){seed.push(("prime",crate::ntheory::isprime(&integer)));}
                if integer<=num_bigint::BigInt::from(1){seed.push(("composite",false));}
                if integer>num_bigint::BigInt::from(1)&&seed.iter().any(|(k,v)|*k=="prime"&&!*v){seed.push(("composite",true));}
            }
        },
        Kind::Const(Constant::Pi|Constant::E)=>seed.extend([("real",true),("positive",true),("transcendental",true)]),
        Kind::Const(Constant::I)=>seed.extend([("imaginary",true),("algebraic",true)]),
        Kind::Const(Constant::Infinity|Constant::NegInfinity)=>{
            seed.extend([("finite",false),("extended_real",true),("extended_positive",e.is_constant(Constant::Infinity))]);
        },
        Kind::Const(Constant::ComplexInfinity)=>seed.extend([("finite",false),("extended_real",false)]),
        Kind::Const(Constant::EulerGamma)=>seed.extend([("real",true),("positive",true)]),
        Kind::Const(Constant::NaN)=>{},
        _=>{
            if let Some(real)=exact_queries(||is_real(e)){seed.push(("real",real));}
            if let Some(finite)=exact_queries(||is_finite(e)){seed.push(("finite",finite));}
            if let Some(zero)=exact_queries(||is_zero(e)){seed.push(("zero",zero));}
            if let Some(rational)=exact_queries(||is_rational(e)){seed.push(("rational",rational));}
            let positive=exact_queries(||is_positive(e));
            let negative=exact_queries(||is_negative(e));
            let need_imaginary=matches!(name,"imaginary"|"complex");
            let quadratic=if need_imaginary||(positive!=Some(true)&&negative!=Some(true)){
                quadratic_facts(e)
            }else{None};
            if positive==Some(true)||quadratic.and_then(|f|f.sign)==Some(1){seed.push(("positive",true));}
            if negative==Some(true)||quadratic.and_then(|f|f.sign)==Some(-1){seed.push(("negative",true));}
            if need_imaginary{
                if let Some(imaginary)=quadratic.and_then(|f|f.imaginary){seed.push(("imaginary",imaginary));}
            }
        },
    }
    crate::assumption_facts::FactSet::from_kwargs(seed).ok()?.closure().ok()?.get(name)
}

#[cfg(test)]
thread_local!{static COMPUTATIONS:std::cell::Cell<usize>=const{std::cell::Cell::new(0)};}
#[inline]
fn calculation(){#[cfg(test)] COMPUTATIONS.with(|count|count.set(count.get()+1));}

#[inline]
fn cached_bool(e:&Expr,query:Query,compute:fn(&Expr)->Option<bool>)->Option<bool>{
    if let Kind::Sym(symbol)=e.kind(){
        if !matches!(query,Query::Parity){return symbol.facts().get(match query{
            Query::Real=>"real",Query::Positive=>"positive",Query::Negative=>"negative",
            Query::Nonnegative=>"nonnegative",Query::Nonpositive=>"nonpositive",Query::Zero=>"zero",
            Query::Nonzero=>"nonzero",Query::Integer=>"integer",Query::Rational=>"rational",
            Query::Finite=>"finite",Query::Parity=>unreachable!(),
        });}
    }
    if let Kind::Func(Func::Undefined(head),_)=e.kind(){
        if !head.internal && !matches!(query,Query::Parity){return head.facts().get(match query{
            Query::Real=>"real",Query::Positive=>"positive",Query::Negative=>"negative",
            Query::Nonnegative=>"nonnegative",Query::Nonpositive=>"nonpositive",Query::Zero=>"zero",
            Query::Nonzero=>"nonzero",Query::Integer=>"integer",Query::Rational=>"rational",
            Query::Finite=>"finite",Query::Parity=>unreachable!(),
        });}
    }
    // Symbols and literal atoms remain as cheap as the declared-fact lookup.
    if e.is_atom(){return compute(e);}
    let mode=numeric_refinement_enabled();
    let shift=query as u32*2+if mode{BANK_BITS}else{0};
    match e.cached_assumption(shift,2){1=>return None,2=>return Some(false),3=>return Some(true),_=>{}}
    calculation();
    let answer=if mode{crate::evalf::with_numeric_cache(||compute(e))}else{compute(e)};
    let encoded=match answer{None=>1,Some(false)=>2,Some(true)=>3};
    // No marker for an in-progress recursion is ever published.
    e.remember_assumption(shift,2,encoded);
    answer
}

/// Use structural proofs only while checking solver candidate assumptions.
/// This is thread local and restores its previous state even during unwinding.
pub(crate) fn exact_queries<R>(query: impl FnOnce() -> R) -> R {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) { NUMERIC_REFINEMENT.with(|enabled| enabled.set(self.0)); }
    }
    let previous = NUMERIC_REFINEMENT.with(|enabled| enabled.replace(false));
    let _restore = Restore(previous);
    query()
}

/// Possible signs of a real quantity, as a set: bit 0 = negative,
/// bit 1 = zero, bit 2 = positive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Signs(u8);

impl Signs {
    pub const NEG: Signs = Signs(1);
    pub const ZERO: Signs = Signs(2);
    pub const POS: Signs = Signs(4);
    pub const ANY: Signs = Signs(7);
    pub const NONNEG: Signs = Signs(6);
    pub const NONPOS: Signs = Signs(3);
    pub const NONZERO: Signs = Signs(5);

    fn union(self, o: Signs) -> Signs {
        Signs(self.0 | o.0)
    }
    fn subset(self, o: Signs) -> bool {
        self.0 & !o.0 == 0
    }
    fn has(self, o: Signs) -> bool {
        self.0 & o.0 != 0
    }

    fn mul(self, o: Signs) -> Signs {
        let mut r = 0u8;
        for a in [Signs::NEG, Signs::ZERO, Signs::POS] {
            if !self.has(a) {
                continue;
            }
            for b in [Signs::NEG, Signs::ZERO, Signs::POS] {
                if !o.has(b) {
                    continue;
                }
                r |= match (a.0, b.0) {
                    (2, _) | (_, 2) => 2,
                    (1, 1) | (4, 4) => 4,
                    _ => 1,
                };
            }
        }
        Signs(r)
    }

    fn add(self, o: Signs) -> Signs {
        // a + b
        if self == Signs::ZERO {
            return o;
        }
        if o == Signs::ZERO {
            return self;
        }
        let nonneg = self.subset(Signs::NONNEG) && o.subset(Signs::NONNEG);
        let nonpos = self.subset(Signs::NONPOS) && o.subset(Signs::NONPOS);
        if nonneg {
            if self == Signs::POS || o == Signs::POS || !(self.has(Signs::ZERO) && o.has(Signs::ZERO)) {
                return Signs::POS;
            }
            return Signs::NONNEG;
        }
        if nonpos {
            if self == Signs::NEG || o == Signs::NEG || !(self.has(Signs::ZERO) && o.has(Signs::ZERO)) {
                return Signs::NEG;
            }
            return Signs::NONPOS;
        }
        Signs::ANY
    }
}

fn number_sign(n: &Number) -> Signs {
    if n.is_positive() {
        Signs::POS
    } else if n.is_negative() {
        Signs::NEG
    } else {
        Signs::ZERO
    }
}

/// Approximate nonzero integral values do not prove membership in the exact
/// integer set. A fractional binary value does prove non-membership, while
/// Float zero is an explicit exception in the public assumptions protocol.
/// Inspect dyadic exponents directly: a huge positive exponent must never
/// materialize its corresponding arbitrary-precision integer here.
fn number_integer(n:&Number)->Option<bool>{
    if !n.is_float(){return Some(n.is_integer());}
    if n.is_zero(){return Some(true);}
    let integral=match n{
        Number::Float(value) if value.is_finite()=>value.fract()==0.0,
        Number::Float(_)=>return None,
        Number::MpFloat(value,_)=>value.exp>=0||value.man.trailing_zeros().is_some_and(|zeros|zeros>=value.exp.unsigned_abs()),
        _=>unreachable!(),
    };
    if integral{None}else{Some(false)}
}

fn declared(e:&Expr)->A {
    match e.kind(){Kind::Sym(s)=>s.assumptions(),Kind::Func(Func::Undefined(f),_)=>f.assumptions,_=>A::NONE}
}
fn declared_signs(a:A)->Option<Signs>{
    if !a.has(A::REAL){return None;}
    Some(if a.has(A::ZERO){Signs::ZERO}else if a.has(A::POSITIVE){Signs::POS}else if a.has(A::NEGATIVE){Signs::NEG}else if a.has(A::NONNEGATIVE){Signs::NONNEG}else if a.has(A::NONPOSITIVE){Signs::NONPOS}else if a.has(A::NONZERO){Signs::NONZERO}else{Signs::ANY})
}

/// Sign set of `e` if `e` is known to be real (finite), else `None`.
#[inline]
pub fn real_signs(e:&Expr)->Option<Signs>{
    if e.is_atom(){return real_signs_uncached(e);}
    let mode=numeric_refinement_enabled();let shift=SIGNS_OFFSET+if mode{BANK_BITS}else{0};
    let cached=e.cached_assumption(shift,4);
    if cached==8{return None;}else if cached!=0{return Some(Signs(cached));}
    calculation();let answer=if mode{crate::evalf::with_numeric_cache(||real_signs_uncached(e))}else{real_signs_uncached(e)};
    e.remember_assumption(shift,4,answer.map_or(8,|signs|signs.0));answer
}
fn real_signs_uncached(e: &Expr) -> Option<Signs> {
    if !declared(e).is_empty(){return declared_signs(declared(e));}
    match e.kind() {
        Kind::Num(n) => Some(number_sign(n)),
        Kind::Const(c) => match c {
            Constant::Pi | Constant::E | Constant::EulerGamma => Some(Signs::POS),
            _ => None,
        },
        Kind::Sym(s) => {
            let a = s.assumptions();
            if !a.has(A::REAL) {
                return None;
            }
            Some(if a.has(A::ZERO) {
                Signs::ZERO
            } else if a.has(A::POSITIVE) {
                Signs::POS
            } else if a.has(A::NEGATIVE) {
                Signs::NEG
            } else if a.has(A::NONNEGATIVE) {
                Signs::NONNEG
            } else if a.has(A::NONPOSITIVE) {
                Signs::NONPOS
            } else if a.has(A::NONZERO) {
                Signs::NONZERO
            } else {
                Signs::ANY
            })
        }
        Kind::Add(c, terms) => {
            let mut s = number_sign(c);
            for (t, k) in terms {
                let ts = real_signs(t)?.mul(number_sign(k));
                s = s.add(ts);
                if s == Signs::ANY {
                    // keep checking realness of the remaining terms
                    for (t2, _) in terms {
                        real_signs(t2)?;
                    }
                    return numeric_refine(e, Signs::ANY);
                }
            }
            Some(s)
        }
        Kind::Mul(c, facs) => {
            let mut s = number_sign(c);
            for (b, x) in facs {
                s = s.mul(pow_signs(b, x)?);
            }
            Some(s)
        }
        Kind::Pow(b, x) => pow_signs(b, x),
        Kind::Piecewise(pieces) => {
            let mut signs=None;for(value,_)in pieces{let next=real_signs(value)?;signs=Some(signs.map_or(next,|old:Signs|old.union(next)));}signs
        }
        Kind::Func(f, args) => func_signs(f, args),
        _ => None,
    }
}

/// Refine an unknown sign of a numeric expression by evaluating it.
fn numeric_refine(e: &Expr, s: Signs) -> Option<Signs> {
    if s != Signs::ANY || !numeric_refinement_enabled() || !crate::functions::is_numeric_value(e) {
        return Some(s);
    }
    match crate::evalf::eval_f64(e) {
        Some(v) if v > 1e-12 * (1.0 + v.abs()) => Some(Signs::POS),
        Some(v) if v < -1e-12 * (1.0 + v.abs()) => Some(Signs::NEG),
        _ => Some(s),
    }
}

fn pow_signs(b: &Expr, x: &Expr) -> Option<Signs> {
    let bs = real_signs(b);
    if let Some(k) = x.as_i64() {
        let bs = bs?;
        if k == 0 {
            return Some(Signs::POS);
        }
        if k < 0 && bs.has(Signs::ZERO) {
            // possibly division by zero
            return None;
        }
        return Some(if k % 2 == 0 {
            if bs.has(Signs::ZERO) {
                if bs == Signs::ZERO {
                    Signs::ZERO
                } else {
                    Signs::NONNEG
                }
            } else {
                Signs::POS
            }
        } else {
            bs
        });
    }
    // non-integer exponent: need a nonnegative base and a real exponent
    let bs = bs?;
    let xs = real_signs(x)?;
    if bs == Signs::POS {
        return Some(Signs::POS);
    }
    if bs.subset(Signs::NONNEG) {
        if xs == Signs::POS {
            return Some(bs);
        }
        if bs == Signs::ZERO {
            return None;
        }
    }
    None
}

fn func_signs(f: &Func, args: &[Expr]) -> Option<Signs> {
    let a = args.first()?;
    match f {
        Func::Product => crate::products::signs(args),
        Func::Exp => {
            real_signs(a)?;
            Some(Signs::POS)
        }
        Func::Abs => {
            let s = if is_zero(a) == Some(false) { Signs::POS } else { Signs::NONNEG };
            Some(s)
        }
        Func::Log => {
            let s = real_signs(a)?;
            if s != Signs::POS {
                return None;
            }
            if crate::functions::is_numeric_value(a) {
                return numeric_refine(&crate::functions::log(a), Signs::ANY);
            }
            Some(Signs::ANY)
        }
        Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc => {
            real_signs(a)?;
            numeric_refine(&Expr::from_kind(Kind::Func(f.clone(), args.to_vec())), Signs::ANY)
        }
        Func::Atan | Func::Sinh | Func::Tanh | Func::Asinh | Func::Erf | Func::Sign => real_signs(a),
        Func::Cosh => {
            real_signs(a)?;
            Some(Signs::POS)
        }
        Func::Floor | Func::Ceiling => real_signs(a).map(|_| Signs::ANY),
        Func::Re | Func::Im => Some(Signs::ANY),
        Func::Gamma => {
            if real_signs(a)? == Signs::POS {
                Some(Signs::POS)
            } else {
                None
            }
        }
        Func::Max | Func::Min => {
            let mut result:Option<Signs>=None;for argument in args{let next=real_signs(argument)?;result=Some(match result{None=>next,Some(previous)=>{
                let mut bits=0;for i in 0..3{for j in 0..3{if previous.0&(1<<i)!=0&&next.0&(1<<j)!=0{bits|=1<<if matches!(f,Func::Max){i.max(j)}else{i.min(j)};}}}Signs(bits)
            }});}result
        }
        _ => None,
    }
}

// ----------------------------------------------------------------------
// public queries
// ----------------------------------------------------------------------

#[inline]
pub fn is_real(e:&Expr)->Option<bool>{cached_bool(e,Query::Real,is_real_uncached)}
fn is_real_uncached(e: &Expr) -> Option<bool> {
    if real_signs(e).is_some() {
        return Some(true);
    }
    if declared(e).has(A::NOT_REAL){return Some(false);}
    match e.kind() {
        Kind::Pow(base,exponent) if exact::rational_power(base,exponent).and_then(|f| f.real).is_some() => exact::rational_power(base,exponent).and_then(|f| f.real),
        Kind::Const(Constant::I) => Some(false),
        Kind::Const(Constant::NaN)=>None,
        Kind::Const(c) if c.is_infinite()=>Some(false),
        Kind::Sym(s) if s.assumptions().has(A::NOT_REAL) => Some(false),
        Kind::Mul(c, facs) if facs.len() == 1 && facs[0].0.is_constant(Constant::I) && facs[0].1.is_one() && !c.is_zero() => {
            Some(false)
        }
        Kind::Add(_,terms) => {
            let mut nonreal=0;
            for (term,_) in terms { match is_real(term) {Some(true)=>{},Some(false)=>nonreal+=1,None=>return None} }
            if nonreal==1 {Some(false)} else if nonreal==0 {Some(true)} else {None}
        }
        Kind::Mul(_,factors) => {
            let mut nonreal=0;
            for (base,exponent) in factors {
                let factor=crate::ops::mul_factor(base,exponent);
                match is_real(&factor) {
                    Some(true) if exact_queries(|| is_zero(&factor))==Some(false)=>{},
                    Some(false)=>nonreal+=1,
                    _=>return None,
                }
            }
            if nonreal==1 {Some(false)} else if nonreal==0 {Some(true)} else {None}
        }
        _ => {
            if numeric_refinement_enabled() && crate::functions::is_numeric_value(e) {
                if let Some(z) = crate::evalf::eval_c64(e) {
                    if z.im.abs() > 1e-10 * (1.0 + z.re.abs()) {
                        return Some(false);
                    }
                }
            }
            None
        }
    }
}

#[inline]
pub fn is_positive(e:&Expr)->Option<bool>{cached_bool(e,Query::Positive,is_positive_uncached)}
fn is_positive_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NOT_POSITIVE) || exact_queries(||is_real(e))==Some(false){return Some(false);}
    let s = real_signs(e)?;
    if s == Signs::POS {
        Some(true)
    } else if !s.has(Signs::POS) {
        Some(false)
    } else {
        None
    }
}

#[inline]
pub fn is_negative(e:&Expr)->Option<bool>{cached_bool(e,Query::Negative,is_negative_uncached)}
fn is_negative_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NOT_NEGATIVE) || exact_queries(||is_real(e))==Some(false){return Some(false);}
    let s = real_signs(e)?;
    if s == Signs::NEG {
        Some(true)
    } else if !s.has(Signs::NEG) {
        Some(false)
    } else {
        None
    }
}

#[inline]
pub fn is_nonnegative(e:&Expr)->Option<bool>{cached_bool(e,Query::Nonnegative,is_nonnegative_uncached)}
fn is_nonnegative_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NOT_NONNEGATIVE){return Some(false);}
    if exact_queries(||is_real(e))==Some(false){return Some(false);}
    let s = real_signs(e)?;
    if s.subset(Signs::NONNEG) {
        Some(true)
    } else if s == Signs::NEG {
        Some(false)
    } else {
        None
    }
}

#[inline]
pub fn is_nonpositive(e:&Expr)->Option<bool>{cached_bool(e,Query::Nonpositive,is_nonpositive_uncached)}
fn is_nonpositive_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NOT_NONPOSITIVE){return Some(false);}
    if exact_queries(||is_real(e))==Some(false){return Some(false);}
    let s = real_signs(e)?;
    if s.subset(Signs::NONPOS) {
        Some(true)
    } else if s == Signs::POS {
        Some(false)
    } else {
        None
    }
}

#[inline]
pub fn is_zero(e:&Expr)->Option<bool>{cached_bool(e,Query::Zero,is_zero_uncached)}
fn is_zero_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NOT_ZERO){return Some(false);}
    if declared(e).has(A::ZERO){return Some(true);}
    match e.kind() {
        Kind::Num(n) => return Some(n.is_zero()),
        Kind::Sym(s) if s.assumptions().has(A::NOT_ZERO) => return Some(false),
        Kind::Sym(s) if s.assumptions().has(A::ZERO) => return Some(true),
        Kind::Const(Constant::I | Constant::Pi | Constant::E | Constant::EulerGamma) => return Some(false),
        Kind::Const(c) if c.is_infinite()=>return Some(false),
        Kind::Mul(c,factors) if !c.is_zero() && factors.iter().all(|(b,x)| exact_queries(||is_zero(&crate::ops::mul_factor(b,x)))==Some(false)) => return Some(false),
        Kind::Pow(base,exponent) => {
            if let Some(facts)=exact::rational_power(base,exponent) {if facts.zero.is_some(){return facts.zero;}}
            if exact_queries(||is_zero(base))==Some(false) && is_finite(exponent)==Some(true) {return Some(false);}
        },
        _ => {}
    }
    if let Some(s) = real_signs(e) {
        if s == Signs::ZERO {
            return Some(true);
        }
        if !s.has(Signs::ZERO) {
            return Some(false);
        }
    }
    if numeric_refinement_enabled() && crate::functions::is_numeric_value(e) {
        if let Some(z) = crate::evalf::eval_c64(e) {
            if z.norm() > 1e-10 {
                return Some(false);
            }
        }
    }
    None
}

#[inline]
pub fn is_nonzero(e:&Expr)->Option<bool>{cached_bool(e,Query::Nonzero,is_nonzero_uncached)}
fn is_nonzero_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NONZERO){return Some(true);}
    if declared(e).has(A::NOT_NONZERO){return Some(false);}
    if let Kind::Sym(s)=e.kind() {
        if s.assumptions().has(A::NOT_NONZERO){return Some(false);}
        if s.assumptions().has(A::NONZERO){return Some(true);}
    }
    match is_real(e) {
        Some(false)=>Some(false),
        Some(true)=>is_zero(e).map(|zero|!zero),
        None=>is_zero(e).filter(|zero|*zero).map(|_|false),
    }
}

#[inline]
pub fn is_integer(e:&Expr)->Option<bool>{cached_bool(e,Query::Integer,is_integer_uncached)}
fn is_integer_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::INTEGER){return Some(true);}
    if declared(e).has(A::NOT_INTEGER){return Some(false);}
    match e.kind() {
        Kind::Num(n) => number_integer(n),
        Kind::Sym(s) => {
            let a = s.assumptions();
            if a.has(A::INTEGER) {
                Some(true)
            } else if a.has(A::NOT_INTEGER) || a.has(A::IRRATIONAL) || a.has(A::NOT_REAL) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Const(Constant::NaN|Constant::EulerGamma)=>None,
        Kind::Const(_) => Some(false),
        Kind::Add(c, terms) => {
            if is_rational(e)==Some(false){return Some(false);}
            if !c.is_integer() {
                // integer + non-integer rational is not an integer
                if terms.iter().all(|(t, k)| k.is_integer() && is_integer(t) == Some(true)) {
                    return Some(false);
                }
                return None;
            }
            for (t, k) in terms {
                if !(k.is_integer() && is_integer(t) == Some(true)) {
                    return None;
                }
            }
            Some(true)
        }
        Kind::Mul(c, facs) => {
            if is_rational(e)==Some(false){return Some(false);}
            if !c.is_integer() {
                return None;
            }
            for (b, x) in facs {
                let ok = is_integer(b) == Some(true) && x.as_i64().is_some_and(|k| k >= 0);
                if !ok {
                    return None;
                }
            }
            Some(true)
        }
        Kind::Pow(b, x) => {
            if let Some(facts)=exact::rational_power(b,x){if facts.integer.is_some(){return facts.integer;}}
            if is_rational(e)==Some(false){return Some(false);}
            if is_integer(b) == Some(true) && x.as_i64().is_some_and(|k| k >= 0) {
                Some(true)
            } else {
                None
            }
        }
        Kind::Piecewise(pieces) if !pieces.is_empty() && pieces.iter().all(|(value,_)|is_integer(value)==Some(true))=>Some(true),
        Kind::Func(Func::Min | Func::Max,args) if !args.is_empty() && args.iter().all(|a|is_integer(a)==Some(true))=>Some(true),
        Kind::Func(Func::Floor | Func::Ceiling, args) => {
            if is_real(&args[0]) == Some(true) {
                Some(true)
            } else {
                None
            }
        }
        Kind::Func(Func::Sign, args) if is_real(&args[0]) == Some(true) => Some(true),
        Kind::Func(Func::Abs, args) if is_integer(&args[0]) == Some(true) => Some(true),
        _ => None,
    }
}

#[inline]
pub fn is_rational(e:&Expr)->Option<bool>{cached_bool(e,Query::Rational,is_rational_uncached)}
fn is_rational_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::RATIONAL){return Some(true);}
    if declared(e).has(A::NOT_RATIONAL)||declared(e).has(A::IRRATIONAL){return Some(false);}
    if exact_queries(||is_real(e))==Some(false){return Some(false);}
    match e.kind() {
        Kind::Num(n) => if n.is_float(){None}else{Some(true)},
        Kind::Sym(s) => {
            let a = s.assumptions();
            if a.has(A::RATIONAL) {
                Some(true)
            } else if a.has(A::IRRATIONAL) || a.has(A::NOT_RATIONAL) || a.has(A::NOT_REAL) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Const(Constant::Pi | Constant::E | Constant::I) => Some(false),
        Kind::Add(c,terms) => {
            if !c.is_rational(){return None;}
            let mut irrational=0;
            for (term,coefficient) in terms {
                if !coefficient.is_rational(){return None;}
                match is_rational(term){Some(true)=>{},Some(false)=>irrational+=1,None=>return None}
            }
            if irrational==0{Some(true)}else if irrational==1{Some(false)}else{None}
        },
        Kind::Mul(c,factors) => {
            if !c.is_rational(){return None;}
            let mut irrational=0;
            for (b,x) in factors {
                let factor=crate::ops::mul_factor(b,x);
                match is_rational(&factor){
                    Some(false)=>irrational+=1,
                    Some(true) if exact_queries(||is_zero(&factor))==Some(false)=>{},
                    _=>return None,
                }
            }
            if irrational==0{Some(true)}else if irrational==1{Some(false)}else{None}
        },
        Kind::Pow(b,x) => {
            if let Some(facts)=exact::rational_power(b,x){if facts.rational.is_some(){return facts.rational;}}
            let k=x.as_i64()?;
            match is_rational(b) {
                Some(true) if k>=0 || exact_queries(||is_zero(b))==Some(false)=>Some(true),
                Some(false) if k==1 || (k==-1 && exact_queries(||is_zero(b))==Some(false))=>Some(false),
                _=>None,
            }
        },
        Kind::Func(Func::Floor | Func::Ceiling | Func::Sign,args) if is_real(&args[0])==Some(true)=>Some(true),
        Kind::Func(Func::Abs,args)=>is_rational(&args[0]).filter(|r|*r),
        _ => None,
    }
}

/// `Some(true)` if even, `Some(false)` if odd, `None` if unknown / not an
/// integer.
#[inline]
pub fn parity_even(e:&Expr)->Option<bool>{cached_bool(e,Query::Parity,parity_even_uncached)}
fn parity_even_uncached(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => {
            if n.is_float()&&n.is_zero(){return Some(true);}
            let b = n.to_bigint()?;
            Some(num_integer::Integer::is_even(&b))
        }
        Kind::Sym(s) => {
            let a = s.assumptions();
            if a.has(A::EVEN) {
                Some(true)
            } else if a.has(A::ODD) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Mul(c, facs) => {
            if is_integer(e) != Some(true) {
                return None;
            }
            if c.to_bigint().is_some_and(|b| num_integer::Integer::is_even(&b)) {
                return Some(true);
            }
            let mut all_odd = c.to_bigint().is_some_and(|b| num_integer::Integer::is_odd(&b));
            for (b, x) in facs {
                match parity_even(b) {
                    Some(true) if x.as_i64().is_some_and(|k| k > 0) => return Some(true),
                    Some(false) => {}
                    _ => all_odd = false,
                }
            }
            if all_odd {
                Some(false)
            } else {
                None
            }
        }
        Kind::Add(c, terms) => {
            let mut even = parity_even(&Expr::num(c.clone()))?;
            for (t, k) in terms {
                let te = parity_even(&crate::ops::add_term(t, k))?;
                even = even == te;
            }
            Some(even)
        }
        Kind::Pow(b, x) => {
            let k = x.as_i64()?;
            if k <= 0 {
                return None;
            }
            parity_even(b)
        }
        _ => None,
    }
}

pub fn is_even(e: &Expr) -> Option<bool> {
    if let Some(symbol)=e.as_sym(){return symbol.facts().get("even");}
    if let Kind::Func(Func::Undefined(head),_)=e.kind(){if !head.internal{return head.facts().get("even");}}
    if is_integer(e)==Some(false){return Some(false);}
    parity_even(e)
}

pub fn is_odd(e: &Expr) -> Option<bool> {
    if let Some(symbol)=e.as_sym(){return symbol.facts().get("odd");}
    if let Kind::Func(Func::Undefined(head),_)=e.kind(){if !head.internal{return head.facts().get("odd");}}
    if is_integer(e)==Some(false){return Some(false);}
    parity_even(e).map(|b| !b)
}

#[inline]
pub fn is_finite(e:&Expr)->Option<bool>{cached_bool(e,Query::Finite,is_finite_uncached)}
fn is_finite_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::FINITE){return Some(true);}
    match e.kind() {
        Kind::Num(_) => Some(true),
        Kind::Const(Constant::NaN)=>None,
        Kind::Const(c) => Some(!c.is_infinite()),
        Kind::Sym(s) => {
            if s.assumptions().has(A::FINITE) {
                Some(true)
            } else {
                None
            }
        }
        _ => {
            if real_signs(e).is_some() {
                Some(true)
            } else {
                None
            }
        }
    }
}

#[cfg(test)]
mod numeric_literal_tests {
    use super::*;
    use crate::bigfloat::BigFloat;
    use num_bigint::BigInt;

    #[test]
    fn dyadic_integer_queries_never_materialize_extreme_exponents(){
        for exp in [1_000_000_000,i64::MAX]{
            let value=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(3),exp},200));
            assert_eq!(is_integer(&value),None);
            assert_eq!(is_even(&value),None);
            assert_eq!(is_odd(&value),None);
            assert_eq!(fact(&value,"noninteger"),None);
            assert_eq!(is_rational(&value),None);
        }
        for exp in [-1_000_000_000,i64::MIN]{
            let value=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(3),exp},200));
            assert_eq!(is_integer(&value),Some(false));
            assert_eq!(is_even(&value),Some(false));
            assert_eq!(is_odd(&value),Some(false));
            assert_eq!(fact(&value,"noninteger"),Some(true));
        }
        let value=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(4),exp:-2},200));
        assert_eq!(is_integer(&value),None);
        let zero=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(0),exp:i64::MIN},200));
        assert_eq!(is_integer(&zero),Some(true));
        assert_eq!(is_even(&zero),Some(true));
        assert_eq!(is_odd(&zero),Some(false));
    }

    #[test]
    fn noninteger_predicates_do_not_turn_even_false_into_odd_true(){
        let value=Expr::rational(1,2);
        assert_eq!(is_even(&value),Some(false));
        assert_eq!(is_odd(&value),Some(false));
        assert_eq!(parity_even(&value),None);
        let nan=Expr::nan();
        assert_eq!(is_real(&nan),None);
        assert_eq!(is_integer(&nan),None);
        assert_eq!(is_rational(&nan),None);
        assert_eq!(is_finite(&nan),None);
        assert_eq!(is_even(&nan),None);
        assert_eq!(is_odd(&nan),None);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::symbol::Assumptions;

    #[test]
    fn exact_rational_powers_and_complex_nonzero() {
        exact_queries(|| {
            for source in ["sqrt(2)","2**(1/3)","2**(2/3)","sqrt(2)/2","1+sqrt(2)","sqrt(2)-1","I","sqrt(2)*I"] {
                let e=crate::parser::parse_str(source).unwrap();
                assert_eq!(is_rational(&e),Some(false),"{source}");
                assert_eq!(is_integer(&e),Some(false),"{source}");
            }
            let i=Expr::i();
            assert_eq!(is_zero(&i),Some(false));assert_eq!(is_nonzero(&i),Some(false));
            let e=crate::parser::parse_str("(-2)**(1/3)").unwrap();
            assert_eq!(is_real(&e),Some(false));
            let e=crate::parser::parse_str("sqrt(2)+sqrt(3)").unwrap();
            assert_eq!(is_rational(&e),None);
        });
    }

    #[test]
    fn bounded_quadratic_solver_facts() {
        for (source,integer,sign) in [("(1+sqrt(2))*(sqrt(2)-1)",true,1),("sqrt(2)-2",false,-1),("-100+sqrt(2)",false,-1),("-sqrt(2)",false,-1),("2-sqrt(8)",false,-1),("-sqrt(2)-1",false,-1)] {
            let e=crate::parser::parse_str(source).unwrap();
            let facts=quadratic_facts(&e).unwrap();
            assert_eq!(facts.integer,Some(integer));assert_eq!(facts.sign,Some(sign));
        }
        let e=crate::parser::parse_str("1+sqrt(2)*I").unwrap();
        let facts=quadratic_facts(&e).unwrap();
        assert_eq!(facts.real,Some(false));assert_eq!(facts.rational,Some(false));assert_eq!(facts.zero,Some(false));
        let held=Expr::from_unevaluated_kind(Kind::Pow(Expr::int(2),Expr::rational(1_000_001,2)));
        assert!(quadratic_facts(&held).is_none());
    }

    #[test]
    fn inference() {
        let p = Expr::symbol_with("p", Assumptions::POSITIVE);
        let n = Expr::symbol_with("n", Assumptions::NEGATIVE);
        let r = Expr::symbol_with("r", Assumptions::REAL);
        let k = Expr::symbol_with("k", Assumptions::INTEGER);
        let x = Expr::symbol("x");
        assert_eq!(is_positive(&(&p + 1)), Some(true));
        assert_eq!(is_positive(&(&p * &n)), Some(false));
        assert_eq!(is_negative(&(&p * &n)), Some(true));
        assert_eq!(is_nonnegative(&r.powi(2)), Some(true));
        assert_eq!(is_positive(&r.powi(2)), None);
        assert_eq!(is_positive(&(r.powi(2) + 1)), Some(true));
        assert_eq!(is_positive(&exp(&r)), Some(true));
        assert_eq!(is_real(&x), None);
        assert_eq!(is_real(&(&r + &p)), Some(true));
        assert_eq!(is_integer(&(2 * &k + 1)), Some(true));
        assert_eq!(is_odd(&(2 * &k + 1)), Some(true));
        assert_eq!(is_even(&(2 * &k)), Some(true));
        assert_eq!(is_positive(&sqrt(&p)), Some(true));
        assert_eq!(is_positive(&sub(&Expr::pi(), &Expr::int(3))), Some(true));
        assert_eq!(is_real(&Expr::i()), Some(false));
    }
}

#[cfg(test)]
mod dag_cache_tests {
    use super::*;
    use std::hash::{Hash,Hasher};
    use std::collections::hash_map::DefaultHasher;
    fn shared_mul(base:Expr,depth:usize)->Expr{
        let mut value=base;let one=Expr::one();
        for _ in 0..depth{value=Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,vec![(value.clone(),one.clone()),(value,one.clone())]));}
        value
    }
    fn difference()->Expr{Expr::from_kind(Kind::Add(Number::ZERO,vec![(Expr::pi(),Number::ONE),(Expr::e(),Number::Int(-1))]))}
    fn facts(e:&Expr)->[Option<bool>;12]{[
        is_real(e),is_positive(e),is_negative(e),is_nonnegative(e),is_nonpositive(e),is_zero(e),is_nonzero(e),
        is_integer(e),is_rational(e),is_even(e),is_odd(e),is_finite(e),
    ]}
    #[test]
    fn deep_shared_arithmetic_queries_cache_known_and_unknown_answers(){
        let value=shared_mul(Expr::symbol_with("r",A::REAL),80);
        COMPUTATIONS.with(|count|count.set(0));
        let expected=exact_queries(||facts(&value));
        assert_eq!(expected,[Some(true),None,None,None,None,None,None,None,None,None,None,Some(true)]);
        let cold=COMPUTATIONS.with(|count|count.get());assert!(cold<=80*28,"{cold} computations");
        for _ in 0..8{assert_eq!(exact_queries(||facts(&value)),expected);}
        assert_eq!(COMPUTATIONS.with(|count|count.get()),cold);
        assert_eq!(value.cached_assumption(Query::Positive as u32*2,2),1,"unknown is not an empty slot");
        let numeric=facts(&value);assert_eq!(numeric,expected);
        let warm=COMPUTATIONS.with(|count|count.get());for _ in 0..8{assert_eq!(facts(&value),numeric);}
        assert_eq!(COMPUTATIONS.with(|count|count.get()),warm);
        let positive=shared_mul(Expr::symbol_with("p",A::POSITIVE),80);
        assert_eq!(is_positive(&positive),Some(true));assert_eq!(is_zero(&positive),Some(false));
        let integer=shared_mul(Expr::symbol_with("k",A::ODD),80);
        assert_eq!(is_integer(&integer),Some(true));assert_eq!(is_even(&integer),Some(false));assert_eq!(is_odd(&integer),Some(true));
    }
    #[test]
    fn query_modes_never_share_numeric_refinements_in_either_order(){
        for numeric_first in [false,true]{
            let e=difference();
            if numeric_first{assert_eq!(is_positive(&e),Some(true));}
            assert_eq!(exact_queries(||is_positive(&e)),None);
            assert_eq!(exact_queries(||is_negative(&e)),None);
            assert_eq!(exact_queries(||real_signs(&e)),Some(Signs::ANY));
            assert_eq!(is_positive(&e),Some(true));assert_eq!(is_zero(&e),Some(false));
            assert_eq!(real_signs(&e),Some(Signs::POS));
            for _ in 0..8{assert_eq!(exact_queries(||is_positive(&e)),None);assert_eq!(is_positive(&e),Some(true));}
        }
    }
    #[test]
    fn concurrent_modes_publish_only_complete_independent_fields(){
        let e=difference();let product=shared_mul(Expr::symbol_with("p",A::POSITIVE),80);
        std::thread::scope(|scope|{
            for thread in 0..8{let e=e.clone();let product=product.clone();scope.spawn(move||{
                for _ in 0..20{
                    if thread%2==0{assert_eq!(exact_queries(||is_positive(&e)),None);}else{assert_eq!(is_positive(&e),Some(true));}
                    assert_eq!(exact_queries(||is_positive(&product)),Some(true));assert_eq!(is_finite(&product),Some(true));
                }
            });}
        });
        assert_eq!(exact_queries(||is_positive(&e)),None);assert_eq!(is_positive(&e),Some(true));
    }
    #[test]
    fn exact_mode_panic_restores_policy_and_caches_do_not_change_identity(){
        let e=difference();let duplicate=difference();let hash=|e:&Expr|{let mut h=DefaultHasher::new();e.hash(&mut h);h.finish()};
        let before=(hash(&e),e.mask(),crate::printer::srepr(&e));
        let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||exact_queries(||{
            assert_eq!(is_positive(&e),None);panic!("expected query unwind");
        })));
        assert!(result.is_err());assert!(numeric_refinement_enabled());assert_eq!(is_positive(&e),Some(true));
        assert_eq!(e,duplicate);assert_eq!(e.cmp_struct(&duplicate),std::cmp::Ordering::Equal);
        assert_eq!(before,(hash(&e),e.mask(),crate::printer::srepr(&e)));
    }
    #[test]
    fn closed_numeric_shared_add_keeps_classification_and_refinement_linear(){
        let value=shared_mul(Expr::one(),80);
        let e=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![(value.clone(),Number::ONE),(value,Number::Int(-1))]));
        assert_eq!(real_signs(&e),Some(Signs::ANY));
        assert!(crate::functions::is_numeric_value(&e));assert_eq!(e.cached_assumption(60,2),3);
        assert_eq!(crate::evalf::eval_f64(&e),Some(0.0));
        assert_eq!(exact_queries(||real_signs(&e)),Some(Signs::ANY));
    }
}
