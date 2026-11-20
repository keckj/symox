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
    type Cache=rustc_hash::FxHashMap<usize,(Expr,Option<bool>)>;
    // SymPy's `_fuzzy_group`: all True gives True, exactly one False with
    // the others True gives False, anything else is unknown.
    fn group(answers:impl Iterator<Item=Option<bool>>)->Option<bool>{
        let mut saw_false=false;
        for answer in answers{
            match answer{
                Some(true)=>{},
                None=>return None,
                Some(false) if saw_false=>return None,
                Some(false)=>saw_false=true,
            }
        }
        Some(!saw_false)
    }
    // The reference decides `is_zero` of a number with its certified
    // numerical evaluation, so the active refinement policy applies here.
    fn nonzero(e:&Expr)->bool{is_zero(e)==Some(false)}
    fn rational_value(e:&Expr)->bool{e.as_number().is_some_and(|n|!n.is_float())}
    fn exponential(argument:&Expr,cache:&mut Cache)->Option<bool>{
        // exp(r*I*pi) is a root of unity; exp of a nonzero algebraic number
        // or of a nonzero rational multiple of pi is transcendental.
        let quotient=crate::ops::mul2(argument,&crate::ops::pow(&crate::ops::mul2(&Expr::pi(),&Expr::i()),&Expr::int(-1)));
        if rational_value(&quotient){return Some(true);}
        // The structural conditions come before the numerical zero test:
        // the answer is the same, and a nested argument whose algebraicity
        // is unknown is not evaluated numerically.
        let transcendental=walk(argument,cache)==Some(true)
            || rational_value(&crate::ops::mul2(argument,&crate::ops::pow(&Expr::pi(),&Expr::int(-1))));
        if transcendental&&nonzero(argument){return Some(false);}
        None
    }
    fn power(base:&Expr,exponent:&Expr,cache:&mut Cache)->Option<bool>{
        if base.is_zero()||base.is_one(){return Some(true);}
        if base.is_constant(Constant::E){return exponential(exponent,cache);}
        if let Some(n)=exponent.as_number().filter(|n|!n.is_float()){
            let algebraic=walk(base,cache);
            if algebraic==Some(false){return Some(n.is_zero());}
            // Every remaining answer is `algebraic` or unknown.
            if algebraic.is_none(){return None;}
            if nonzero(base){return if n.is_zero(){algebraic.filter(|a|*a)}else{algebraic};}
            if n.is_positive(){return algebraic;}
            return None;
        }
        if walk(base,cache)==Some(true)&&walk(exponent,cache)==Some(true){
            let one=crate::ops::add2(base,&Expr::int(-1));
            if (nonzero(base)&&nonzero(&one))||exact_queries(||is_integer(base))==Some(false){
                return exact_queries(||is_rational(exponent));
            }
        }
        None
    }
    fn walk(e:&Expr,cache:&mut Cache)->Option<bool>{
        if let Some((_,answer))=cache.get(&e.dag_id()){return *answer;}
        let answer=match e.kind(){
            Kind::Sym(symbol)=>symbol.facts().get("algebraic"),
            Kind::Func(Func::Undefined(head),_) if !head.internal=>head.facts().get("algebraic"),
            Kind::Num(number) if !number.is_float()||number.is_zero()=>Some(true),
            Kind::Const(Constant::I|Constant::GoldenRatio|Constant::TribonacciConstant)=>Some(true),
            Kind::Const(Constant::Pi|Constant::E|Constant::Infinity|Constant::NegInfinity|Constant::ComplexInfinity)=>Some(false),
            Kind::Pow(base,exponent)=>power(base,exponent,cache),
            Kind::Add(c,terms) if !c.is_float()&&terms.iter().all(|(_,c)|!c.is_float())=>{
                let answers:Vec<_>=terms.iter().map(|(term,_)|walk(term,cache)).collect();
                group(answers.into_iter())
            },
            Kind::Mul(c,factors) if !c.is_float()=>{
                let answers:Vec<_>=factors.iter().map(|(base,exponent)|power(base,exponent,cache)).collect();
                match group(answers.into_iter()){
                    // One transcendental factor among nonzero algebraic ones.
                    Some(false) if factors.iter().all(|(base,_)|nonzero(base))=>Some(false),
                    Some(false)=>None,
                    answer=>answer,
                }
            },
            Kind::Func(Func::Exp,args) if args.len()==1=>exponential(&args[0],cache),
            Kind::Func(Func::Sin|Func::Cos|Func::Tan|Func::Cot|Func::Sec|Func::Csc,args) if args.len()==1=>{
                let argument=&args[0];
                if walk(argument,cache)==Some(true)&&nonzero(argument){Some(false)}
                else if crate::special::trig_eval::pi_coeff(argument).is_some_and(|c|rational_value(&c)){Some(true)}
                else{None}
            },
            Kind::Func(Func::Log,args) if args.len()==1=>{
                let shifted=crate::ops::add2(&args[0],&Expr::int(-1));
                if shifted.is_zero(){Some(true)}
                else if walk(&args[0],cache)==Some(true)&&nonzero(&shifted){Some(false)}
                else{None}
            },
            Kind::Func(Func::Re|Func::Im|Func::Abs|Func::Conjugate,args) if args.len()==1=>walk(&args[0],cache),
            Kind::Func(Func::Max|Func::Min,args)=>{
                let answers:Vec<_>=args.iter().map(|a|walk(a,cache)).collect();
                if answers.iter().all(|a|*a==Some(true)){Some(true)}
                else if answers.iter().all(|a|*a==Some(false)){Some(false)}
                else{None}
            },
            _=>None,
        };
        cache.insert(e.dag_id(),(e.clone(),answer));answer
    }
    walk(e,&mut rustc_hash::FxHashMap::default())
}

/// `Add._eval_is_extended_positive` (or `_negative` when `positive` is
/// false) over the terms, for a sum that is not a number: every term is
/// extended nonnegative and one is extended positive, with consistent
/// infinite terms.
fn add_extended_sign(e:&Expr,positive:bool)->Option<bool>{
    let (strict,weak,opposite)=if positive{("extended_positive","extended_nonnegative","extended_nonpositive")}else{("extended_negative","extended_nonpositive","extended_nonnegative")};
    let args:Vec<Expr>=e.args().into_iter().filter(|a|exact_queries(||is_zero(a))!=Some(true)).collect();
    if args.is_empty(){return Some(false);}
    let (mut pos,mut nonneg,mut nonpos,mut unknown)=(false,false,false,false);
    let mut saw_inf:Vec<Option<bool>>=Vec::new();
    for a in &args{
        let ispos=fact(a,strict);
        let infinite=exact_queries(||is_finite(a)).map(|f|!f);
        if infinite==Some(true){
            let v=match (ispos,fact(a,weak)){(Some(true),_)|(_,Some(true))=>Some(true),(Some(false),Some(false))=>Some(false),_=>None};
            if !saw_inf.contains(&v){saw_inf.push(v);}
            if saw_inf.contains(&Some(true))&&saw_inf.contains(&Some(false)){return None;}
        }
        if ispos==Some(true){pos=true;continue;}
        if fact(a,weak)==Some(true){nonneg=true;continue;}
        if fact(a,opposite)==Some(true){nonpos=true;continue;}
        infinite?;
        unknown=true;
    }
    if !saw_inf.is_empty(){
        return if saw_inf.len()>1{None}else{saw_inf[0]};
    }
    if unknown{return None;}
    if !nonpos&&pos{return Some(true);}
    if !pos&&!nonneg{return Some(false);}
    None
}

/// An exact declared fact, or a conservative native proof for a scalar value.
/// Additional Symbol facts never fall back to numerical refinement.
thread_local! {
    /// Answers of `fact` during the outermost call: `structural_fact` asks
    /// the same question of every argument (sums, products, `exp`), so on a
    /// graph that reuses sub-expressions each shared node is queried once
    /// instead of once per occurrence in the expanded tree. Keyed by node,
    /// predicate and numeric policy (`exact_queries` may switch it inside a
    /// query); the node is held so that its id cannot be reused.
    static FACTS: std::cell::RefCell<Option<rustc_hash::FxHashMap<(usize, Box<str>, bool), (Expr, Option<bool>)>>> =
        const { std::cell::RefCell::new(None) };
}

pub fn fact(e:&Expr,name:&str)->Option<bool>{
    if e.is_atom(){return fact_uncached(e,name);}
    struct Scope(bool);
    impl Drop for Scope {
        fn drop(&mut self) {
            if self.0 { FACTS.with(|memo| *memo.borrow_mut() = None); }
        }
    }
    let outermost = FACTS.with(|memo| {
        let mut memo = memo.borrow_mut();
        let outermost = memo.is_none();
        if outermost { *memo = Some(rustc_hash::FxHashMap::default()); }
        outermost
    });
    let _scope = Scope(outermost);
    let key = (e.dag_id(), Box::<str>::from(name), numeric_refinement_enabled());
    if let Some(hit) = FACTS.with(|memo| memo.borrow().as_ref().and_then(|m| m.get(&key).map(|(_, r)| *r))) {
        return hit;
    }
    let answer = fact_uncached(e, name);
    FACTS.with(|memo| {
        if let Some(m) = memo.borrow_mut().as_mut() { m.insert(key, (e.clone(), answer)); }
    });
    answer
}

fn fact_uncached(e:&Expr,name:&str)->Option<bool>{
    if name=="commutative" && e.has_noncommutative_host(){return crate::decompose::commutativity(e);}
    if let Some(answer)=crate::host::fact(e,name){return answer;}
    // SymPy's assumption rules make every one of these facts imply
    // `commutative` (`extended_real -> commutative`, `complex ->
    // commutative`), so an expression that does not commute has them all
    // False: `(x*A).is_extended_real` and `(A + x).is_zero` are False.
    if implies_commutative(name)&&nc_settles(e){return Some(false);}
    if name=="polar" &&crate::special::polar::is_polar(e){return Some(true);}
    // Mul._eval_is_polar always settles (any/all over the factors), and
    // Pow defers to its base.
    if name=="polar"{
        match e.kind(){
            Kind::Mul(..)=>return Some(false),
            Kind::Pow(b,_)=>return fact(b,"polar"),
            _=>{}
        }
    }
    if matches!(e.kind(),Kind::Add(..)) && !e.is_number() {
        let sign=match name{"extended_positive"|"extended_nonnegative"=>Some(true),"extended_negative"|"extended_nonpositive"=>Some(false),_=>None};
        if let Some(positive)=sign{
            if add_extended_sign(e,positive)==Some(true)&&structural_fact(e,"extended_real")!=Some(false){return Some(true);}
        }
    }
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
    let native=match name {
        "real"=>Some(is_real(e)),"positive"=>Some(is_positive(e)),"negative"=>Some(is_negative(e)),
        "nonnegative"=>Some(is_nonnegative(e)),"nonpositive"=>Some(is_nonpositive(e)),
        "zero"=>Some(is_zero(e)),"nonzero"=>Some(is_nonzero(e)),"integer"=>Some(is_integer(e)),
        "rational"=>Some(is_rational(e)),"finite"=>return is_finite(e),"even"=>Some(is_even(e)),"odd"=>Some(is_odd(e)),
        "commutative"=>return crate::decompose::commutativity(e),_=>None,
    };
    if let Some(answer)=native{
        // Every one of these facts implies extended_real: a sum, product or
        // exponential proven non-real settles them (SymPy derives them from
        // ``is_extended_real`` through the assumption rules).
        return settle_nonreal(e,answer);
    }
    if implies_extended_real(name)&&structural_kind(e)&&structural_fact(e,"extended_real")==Some(false){return Some(false);}
    // `Expr._eval_is_extended_positive_negative`: the sign of a finite
    // number comes from its evaluation, settled by the minimal polynomial
    // when the evaluation has no significance.
    if matches!(name,"extended_positive"|"extended_negative"|"extended_nonnegative"|"extended_nonpositive")
        && !e.is_atom() && crate::functions::is_numeric_value(e) && exact_queries(||is_finite(e))==Some(true) {
        let answer=match name{
            "extended_positive"=>is_positive(e),"extended_negative"=>is_negative(e),
            "extended_nonnegative"=>is_nonnegative(e),_=>is_nonpositive(e),
        };
        if answer.is_some(){return answer;}
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
        Kind::Const(Constant::EulerGamma|Constant::Catalan)=>seed.extend([("real",true),("positive",true)]),
        Kind::Const(Constant::GoldenRatio|Constant::TribonacciConstant)=>seed.extend([("real",true),("positive",true),("irrational",true),("algebraic",true)]),
        Kind::Const(Constant::NaN)=>{},
        _=>{
            match e.kind() {
                Kind::Func(Func::Re | Func::Im, _) => seed.push(("extended_real",true)),
                Kind::Func(Func::Abs, _) => seed.push(("extended_nonnegative",true)),
                _ => {},
            }
            if let Kind::Func(Func::Special(h),args)=e.kind(){
                if crate::special::numtheory::not_prime(*h){seed.push(("prime",false));}
                if let Some(integer)=crate::special::numtheory::integer_fact(*h,args){seed.push(("integer",integer));}
            }
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
            let mut numeric=false;
            if matches!(name,"extended_nonnegative"|"extended_nonpositive"|"extended_positive"|"extended_negative")&&positive!=Some(true)&&negative!=Some(true){
                if exact_queries(||is_nonnegative(e))==Some(true){seed.push(("nonnegative",true));}
                else if exact_queries(||is_nonpositive(e))==Some(true){seed.push(("nonpositive",true));}
                else {numeric=true;}
            }
            if (numeric||matches!(name,"extended_real"|"positive"|"negative"))&&positive!=Some(true)&&negative!=Some(true){
                match reference_numeric_sign(e){
                    Some(NumericSign::Positive)=>seed.push(("extended_positive",true)),
                    Some(NumericSign::Negative)=>seed.push(("extended_negative",true)),
                    None=>{},
                }
            }
            for structural in structural_names(name){
                if let Some(value)=structural_fact(e,structural){seed.push((structural,value));}
            }
            if need_imaginary{
                if let Some(imaginary)=quadratic.and_then(|f|f.imaginary){seed.push(("imaginary",imaginary));}
                // A reduced rational with denominator two has odd numerator.
                // Its principal power of a strictly negative finite real is
                // purely imaginary, including negative half-integer powers.
                if let Kind::Pow(base,exponent)=e.kind(){
                    if exponent.as_number().and_then(Number::numer_denom).is_some_and(|(_,q)|q==2.into())
                        && exact_queries(||is_negative(base))==Some(true){
                        seed.push(("imaginary",true));
                    }
                    // Pow._eval_is_imaginary: an imaginary base to an odd
                    // integer power is imaginary (even powers are real).
                    if exponent.as_i64().is_some_and(|k|k%2!=0) && fact(base,"imaginary")==Some(true){
                        seed.push(("imaginary",true));
                    }
                }
            }
        },
    }
    crate::assumption_facts::FactSet::from_kwargs(seed).ok()?.closure().ok()?.get(name)
}

/// The structural handlers that can settle a queried fact through the
/// closure of the assumption rules.
fn structural_names(name:&str)->&'static [&'static str]{
    match name{
        "complex"|"infinite"=>&["complex"],
        "transcendental"=>&["complex","algebraic"],
        "algebraic"=>&["algebraic","complex"],
        "irrational"=>&["algebraic"],
        "extended_real"=>&["extended_real"],
        "imaginary"=>&["imaginary"],
        "hermitian"=>&["hermitian","extended_real"],
        "antihermitian"=>&["antihermitian","imaginary"],
        _=>&[],
    }
}

/// A native answer to a fact that implies extended_real (real, zero,
/// nonzero, integer, rational and the signs), settled to False when the
/// native inference is undecided on a number that SymPy's sum, product or
/// exponential handlers prove non-real. SymPy derives these facts from
/// ``is_extended_real`` through the assumption rules.
pub fn settle_nonreal(e:&Expr,answer:Option<bool>)->Option<bool>{
    if answer.is_some(){return answer;}
    // Each of these facts implies real, hence commutative.
    if nc_settles(e){return Some(false);}
    if is_real(e)==Some(false){return Some(false);}
    if answer.is_none()&&structural_kind(e)&&crate::symbolic_number::is_symbolic_number(e)&&fact(e,"extended_real")==Some(false){return Some(false);}
    answer
}

/// True when `e` provably does not commute.
fn noncommuting(e:&Expr)->bool{e.has_noncommutative_host()&&crate::decompose::commutativity(e)==Some(false)}

/// A noncommutative expression whose facts implying `commutative` are all
/// False. SymPy settles them so for symbols, sums, products and functions,
/// but a noncommutative Pow only gets what its own handlers decide
/// (`(A**2).is_positive` is None, `(A**-1).is_real` False), so powers keep
/// the native answers.
fn nc_settles(e:&Expr)->bool{!matches!(e.kind(),Kind::Pow(..))&&noncommuting(e)}

/// Facts whose closure under SymPy's assumption rules contains
/// `commutative`.
pub fn implies_commutative(name:&str)->bool{
    matches!(name,"algebraic"|"complex"|"composite"|"even"|"extended_negative"|"extended_nonnegative"|"extended_nonpositive"
        |"extended_nonzero"|"extended_positive"|"extended_real"|"imaginary"|"integer"|"irrational"|"negative"|"noninteger"
        |"nonnegative"|"nonpositive"|"nonzero"|"odd"|"positive"|"prime"|"rational"|"real"|"transcendental"|"zero")
}

/// Facts that require an extended real value.
fn implies_extended_real(name:&str)->bool{
    matches!(name,"extended_positive"|"extended_negative"|"extended_nonnegative"|"extended_nonpositive"
        |"extended_nonzero"|"irrational"|"noninteger"|"prime"|"composite")
}

/// The expression kinds that ``structural_fact`` handles.
fn structural_kind(e:&Expr)->bool{
    matches!(e.kind(),Kind::Add(..)|Kind::Mul(..)|Kind::Func(Func::Exp,_))
}

/// SymPy's ``Add._eval_is_imaginary``.
fn add_imaginary(args:&[Expr])->Option<bool>{
    let mut nonzero=Vec::new();
    let mut imaginary=Vec::new();
    for a in args{
        if fact(a,"extended_real")==Some(true){
            match is_zero(a){Some(true)=>{},Some(false)=>nonzero.push(a.clone()),None=>return None}
        }else if fact(a,"imaginary")==Some(true){
            imaginary.push(a.clone());
        }else if let Kind::Mul(..)=a.kind(){
            let factors=a.args();
            let units=factors.iter().filter(|f|f.is_constant(Constant::I)).count();
            if units!=1{return None;}
            let rest:Vec<Expr>=factors.iter().filter(|f|!f.is_constant(Constant::I)).cloned().collect();
            if fact(&crate::ops::mul(&rest),"extended_real")!=Some(true){return None;}
            imaginary.push(a.clone());
        }else{return None;}
    }
    if nonzero.len()==args.len(){return None;}
    match is_zero(&crate::ops::add(&nonzero)){
        // The zeroness of a sum of the ``a*I`` is that of the sum of the ``a``.
        Some(true)=>is_zero(&crate::ops::add(&imaginary)).map(|zero|!zero),
        Some(false)=>Some(false),
        None=>None,
    }
}

/// SymPy's ``_fuzzy_group``.
fn fuzzy_group(values:impl IntoIterator<Item=Option<bool>>,quick_exit:bool)->Option<bool>{
    let mut saw_other=false;
    for value in values{
        match value{
            Some(true)=>{},
            None=>return None,
            Some(false)=>{
                if quick_exit&&saw_other{return None;}
                saw_other=true;
            },
        }
    }
    Some(!saw_other)
}

fn fuzzy_or(values:impl IntoIterator<Item=Option<bool>>)->Option<bool>{
    let mut result=Some(false);
    for value in values{
        match value{Some(true)=>return Some(true),None=>result=None,Some(false)=>{}}
    }
    result
}

/// SymPy's ``Mul._eval_real_imag``.
fn mul_real_imag(args:&[Expr],mut real:bool)->Option<bool>{
    let mut zero=Some(false);
    let mut not_re_im:Option<&Expr>=None;
    for t in args{
        let extended_real=fact(t,"extended_real");
        let complex_or_infinite=match (fact(t,"complex"),fact(t,"infinite")){
            (Some(true),_)|(_,Some(true))=>Some(true),(Some(false),Some(false))=>Some(false),_=>None,
        };
        if complex_or_infinite==Some(false)&&extended_real==Some(false){return Some(false);}
        if fact(t,"imaginary")==Some(true){real=!real;}
        else if extended_real==Some(true){
            if zero!=Some(true){
                let z=is_zero(t);
                if z!=Some(true)&&zero==Some(false){zero=z;}
                else if z==Some(true){
                    if args.iter().all(|a|is_finite(a)==Some(true)){return Some(true);}
                    return None;
                }
            }
        }
        else if extended_real==Some(false)||fact(t,"imaginary")==Some(false){
            if not_re_im.is_some(){return None;}
            not_re_im=Some(t);
        }
        else{return None;}
    }
    if let Some(t)=not_re_im{
        if fact(t,"extended_real")==Some(false)&&real{return zero;}
        if fact(t,"imaginary")==Some(false)&&!real{return zero;}
        None
    }else if zero==Some(false)||real{Some(real)}else{None}
}

/// SymPy's ``Mul._eval_herm_antiherm``.
fn mul_herm_antiherm(e:&Expr,args:&[Expr],mut herm:bool)->Option<bool>{
    for t in args{
        let (h,a)=(fact(t,"hermitian"),fact(t,"antihermitian"));
        if h.is_none()||a.is_none(){return None;}
        if h==Some(true){continue;}
        if a==Some(true){herm=!herm;}else{return None;}
    }
    if herm{return Some(true);}
    match is_zero(e){Some(true)=>Some(true),Some(false)=>Some(false),None=>None}
}

/// The SymPy ``Add``, ``Mul``, ``Pow`` and ``exp`` handlers of the facts
/// that the native sign inference does not derive: complex, algebraic,
/// extended real, imaginary, hermitian and antihermitian.
fn structural_fact(e:&Expr,name:&str)->Option<bool>{
    match e.kind(){
        Kind::Add(..)=>match name{
            "complex"|"algebraic"|"extended_real"|"hermitian"|"antihermitian"=>
                fuzzy_group(e.args().iter().map(|a|fact(a,name)),true),
            "imaginary"=>add_imaginary(&e.args()),
            _=>None,
        },
        Kind::Mul(..)=>{
            let args=e.args();
            match name{
                "complex"=>{
                    let comp=fuzzy_group(args.iter().map(|a|fact(a,"complex")),false);
                    if comp==Some(false)&&args.iter().any(|a|fact(a,"infinite")==Some(true)){
                        if args.iter().any(|a|is_zero(a)!=Some(false)){return None;}
                        return Some(false);
                    }
                    comp
                },
                "algebraic"=>{
                    let r=fuzzy_group(args.iter().map(|a|fact(a,"algebraic")),true);
                    if r==Some(false)&&!args.iter().all(|a|is_zero(a)==Some(false)){return None;}
                    r
                },
                "extended_real"=>mul_real_imag(&args,true),
                "imaginary"=>{
                    if args.iter().all(|a|is_zero(a)==Some(false)&&is_finite(a)==Some(true)){mul_real_imag(&args,false)}else{None}
                },
                "hermitian"|"antihermitian"=>mul_herm_antiherm(e,&args,name=="hermitian"),
                _=>None,
            }
        },
        Kind::Pow(base,exponent) if name=="complex"=>{
            if base.is_constant(Constant::E){
                return fuzzy_or([fact(exponent,"complex"),fact(exponent,"extended_negative")]);
            }
            if fact(base,"complex")==Some(true)&&fact(exponent,"complex")==Some(true)&&is_finite(e)==Some(true){Some(true)}else{None}
        },
        Kind::Func(Func::Exp,args) if args.len()==1=>match name{
            "complex"=>fuzzy_or([fact(&args[0],"complex"),fact(&args[0],"extended_negative")]),
            // ``exp._eval_is_extended_real``.
            "extended_real" if fact(&args[0],"extended_real")==Some(true)=>Some(true),
            "extended_real" if fact(&args[0],"imaginary")==Some(true)=>{
                let turns=crate::ops::div(&crate::ops::mul(&[Expr::int(-2),Expr::i(),args[0].clone()]),&Expr::constant(Constant::Pi));
                is_even(&turns)
            },
            // ``exp._eval_is_algebraic``: a nonzero algebraic exponent gives
            // a transcendental value (Lindemann).
            "algebraic" if is_zero(&args[0])==Some(false)&&fact(&args[0],"algebraic")==Some(true)=>Some(false),
            _=>None,
        },
        _=>None,
    }
}

#[cfg(test)]
thread_local!{static COMPUTATIONS:std::cell::Cell<usize>=const{std::cell::Cell::new(0)};}
#[inline]
fn calculation(){#[cfg(test)] COMPUTATIONS.with(|count|count.set(count.get()+1));}

fn query_name(query:Query)->&'static str {
    match query{
        Query::Real=>"real",Query::Positive=>"positive",Query::Negative=>"negative",
        Query::Nonnegative=>"nonnegative",Query::Nonpositive=>"nonpositive",Query::Zero=>"zero",
        Query::Nonzero=>"nonzero",Query::Integer=>"integer",Query::Rational=>"rational",
        Query::Finite=>"finite",Query::Parity=>"even",
    }
}

#[inline]
fn cached_bool(e:&Expr,query:Query,compute:fn(&Expr)->Option<bool>)->Option<bool>{
    if e.is_host() && !matches!(query,Query::Parity) {
        if let Some(answer)=crate::host::fact(e,query_name(query)){return answer;}
    }
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
    let answer=if mode{crate::evalf::with_numeric_cache(||crate::evalf_mp::with_numeric_cache(||compute(e)))}else{compute(e)};
    // An answer derived while a host hook failed is not remembered.
    if e.has_host() && crate::host::failure_pending(){return answer;}
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

/// Decide zero by structural proofs, preserving None when unresolved.
pub fn is_zero_exact(expression: &Expr) -> Option<bool> {
    exact_queries(|| is_zero(expression))
}

/// Query a named assumption using structural proofs only.
pub fn fact_exact(expression: &Expr, name: &str) -> Option<bool> {
    exact_queries(|| fact(expression, name))
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
        Number::MpFloat(value,_)=>value.exp>=0||value.man.trailing_zeros().is_some_and(|zeros|u128::from(zeros)>=value.exp.unsigned_abs()),
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
    calculation();let answer=if mode{crate::evalf::with_numeric_cache(||crate::evalf_mp::with_numeric_cache(||real_signs_uncached(e)))}else{real_signs_uncached(e)};
    if e.has_host() && crate::host::failure_pending(){return answer;}
    e.remember_assumption(shift,4,answer.map_or(8,|signs|signs.0));answer
}
/// The sign set of a host atom or application, from the host's answers to
/// its sign queries; `None` when the host keeps the native facts.
fn host_signs(e:&Expr)->Option<Option<Signs>>{
    if crate::host::fact(e,"real")?!=Some(true){return Some(None);}
    let ask=|name:&str|crate::host::fact(e,name).flatten();
    Some(Some(if ask("zero")==Some(true){Signs::ZERO}else if ask("positive")==Some(true){Signs::POS}
        else if ask("negative")==Some(true){Signs::NEG}else if ask("nonnegative")==Some(true){
            if ask("nonzero")==Some(true){Signs::POS}else{Signs::NONNEG}
        }else if ask("nonpositive")==Some(true){
            if ask("nonzero")==Some(true){Signs::NEG}else{Signs::NONPOS}
        }else if ask("nonzero")==Some(true){Signs::NONZERO}else{Signs::ANY}))
}
fn real_signs_uncached(e: &Expr) -> Option<Signs> {
    if e.is_host(){if let Some(signs)=host_signs(e){return signs;}}
    if !declared(e).is_empty(){return declared_signs(declared(e));}
    if let Some(sign)=exact::cardano_sign(e){
        return Some(match sign{-1=>Signs::NEG,0=>Signs::ZERO,_=>Signs::POS});
    }
    match e.kind() {
        Kind::Num(n) => Some(number_sign(n)),
        Kind::Const(c) => match c {
            c if c.is_number_symbol() => Some(Signs::POS),
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
                match pow_signs(b, x) {
                    Some(f) => s = s.mul(f),
                    None if facs.len() > 1 && !noncommuting(e) => return imaginary_factors_signs(facs).map(|f| number_sign(c).mul(f)),
                    None => return None,
                }
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

/// Numerical refinement is checked at two increasing multiprecisions. It
/// never declares an exact zero from numerical cancellation. The exact-query
/// policy remains separate and does not enter this path.
fn numeric_samples(e: &Expr) -> Option<(crate::evalf_mp::MpComplex, crate::evalf_mp::MpComplex, crate::bigfloat::BigFloat)> {
    exact_queries(|| {
        let first = crate::evalf_mp::eval_mp_unbounded(e,160)?;
        let last = crate::evalf_mp::eval_mp_unbounded(e,320)?;
        let mut scale = first.abs(320);
        let magnitude = last.abs(320);
        if magnitude.cmp(&scale).is_gt() { scale=magnitude; }
        if e.is_add() {
            for term in e.args() {
                let magnitude=crate::evalf_mp::eval_mp_unbounded(&term,320)?.abs(320);
                if magnitude.cmp(&scale).is_gt() { scale=magnitude; }
            }
        }
        Some((first,last,scale))
    })
}

fn stable_component(first:&crate::bigfloat::BigFloat,last:&crate::bigfloat::BigFloat,scale:&crate::bigfloat::BigFloat)->bool {
    !last.is_zero() && first.sign()==last.sign()
        && last.abs().cmp(&scale.ldexp(-96)).is_gt()
        && first.sub(last,320).abs().cmp(&last.abs().ldexp(-64)).is_le()
}

/// Refine an unknown sign only when its real component is stable.
fn numeric_refine(e: &Expr, s: Signs) -> Option<Signs> {
    if s != Signs::ANY || !numeric_refinement_enabled() || !crate::functions::is_numeric_value(e) {
        return Some(s);
    }
    // An evaluation without significance (a cancellation exhausting the
    // precision budget, or a null value) of an algebraic number is settled
    // by its minimal polynomial, as in the reference.
    let null=match numeric_samples(e) {
        Some((first,last,scale)) => {
            if first.im.is_zero() && last.im.is_zero() && stable_component(&first.re,&last.re,&scale) {
                return Some(if last.re.sign()>0 {Signs::POS} else {Signs::NEG});
            }
            numerically_null(&first,&last,&scale)
        }
        None => true,
    };
    if null && algebraic_zero(e) {
        return Some(Signs::ZERO);
    }
    Some(s)
}

enum NumericSign{Positive,Negative}

/// `Expr._eval_is_extended_positive_negative` for the numbers whose class
/// keeps it (every function but `Abs`, `exp`, `log` and `Piecewise`, and a
/// numeric `Add`): the sign of the numerical value when it has
/// significance, and no answer when the evaluation fails, has no
/// significance or is complex. It does
/// not need the value to be finite (`polylog(3, -1/2)` is extended
/// negative while its finiteness stays unknown).
fn reference_numeric_sign(e:&Expr)->Option<NumericSign>{
    let inherits=match e.kind(){
        Kind::Func(Func::Abs|Func::Exp|Func::Log|Func::Undefined(_),_)=>false,
        Kind::Func(..)=>true,
        // A sum of radicals keeps the exact sign machinery (its evaluation
        // can be costly on solver roots); a sum with function values takes
        // the numerical sign.
        // Cheap cached test first; the function search visits each shared
        // node once (a tree walk is exponential on a deep shared graph).
        Kind::Add(..)=>crate::symbolic_number::is_symbolic_number(e)
            &&crate::traverse::any_node(e,&mut |s| matches!(s.kind(),Kind::Func(..))),
        _=>false,
    };
    if !inherits||!numeric_refinement_enabled()||!crate::symbolic_number::is_symbolic_number(e)||!crate::functions::is_numeric_value(e)||crate::poly::rootof::is_rootof(e){return None;}
    let (first,last,scale)=numeric_samples(e)?;
    // `_eval_evalf(2)` of a complex value is a sum without `_prec`: no
    // answer (`polylog(2, 3)`).
    if stable_component(&first.im,&last.im,&scale){return None;}
    if first.im.is_zero()&&last.im.is_zero()&&stable_component(&first.re,&last.re,&scale){
        return Some(if last.re.sign()>0{NumericSign::Positive}else{NumericSign::Negative});
    }
    None
}

fn numerically_null(first:&crate::evalf_mp::MpComplex,last:&crate::evalf_mp::MpComplex,scale:&crate::bigfloat::BigFloat)->bool {
    let tiny=|v:&crate::bigfloat::BigFloat|v.is_zero()||v.abs().cmp(&scale.ldexp(-96)).is_le();
    tiny(&first.re)&&tiny(&first.im)&&tiny(&last.re)&&tiny(&last.im)
}

/// The reference `Expr._eval_is_extended_positive_negative`: a number
/// whose evaluation has no significance, which is algebraic and free of
/// functions, is exactly zero when its minimal polynomial is the generator.
pub(crate) fn algebraic_zero(e:&Expr)->bool {
    // A root object is a number, not a `Function`, for the reference.
    fn has_function(e:&Expr)->bool {
        if crate::poly::rootof::is_rootof(e) {return false;}
        matches!(e.kind(),Kind::Func(..))||e.args().iter().any(has_function)
    }
    // A sum whose identical terms cancel has the generator as minimal
    // polynomial. The check compares shared terms by identity, so it stays
    // linear on a shared graph.
    if let Kind::Add(c,terms)=e.kind() {
        if c.is_zero() && !terms.is_empty() && terms.iter().all(|(t,_)|{
            let mut sum=Number::ZERO;
            for (u,k) in terms { if u==t { sum=sum.add(k); } }
            sum.is_zero()
        }) {return true;}
    }
    // The tree walks below are not shared-DAG aware: a graph whose
    // expanded size is huge is left undecided.
    if crate::traverse::node_count(e)>4096 || has_function(e) {return false;}
    exact_queries(||{
        if algebraic_proof(e)!=Some(true) {return false;}
        matches!(crate::poly::minpoly::minimal_polynomial_qq(e),Ok(p) if p.len()==2&&p[0].is_zero()&&p[1].is_one())
    })
}

fn pow_signs(b: &Expr, x: &Expr) -> Option<Signs> {
    let bs = real_signs(b);
    if let Some(k) = x.as_i64() {
        // Pow._eval_is_extended_positive/negative: an imaginary base (finite
        // and nonzero) to an even integer power is positive when the
        // exponent is 0 mod 4 and negative otherwise ((I*t)**2 = -t**2).
        if bs.is_none() && k % 2 == 0 && fact(b, "imaginary") == Some(true) {
            return Some(if k % 4 == 0 { Signs::POS } else { Signs::NEG });
        }
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
    // Pow._eval_is_extended_positive: an imaginary base to an imaginary
    // power is positive when log(base) is imaginary (I**I = exp(-pi/2)).
    // Its first branch tests `base == exp` and leaves the sign of I**I
    // unknown; the value is still real and nonzero.
    if bs.is_none() && !matches!(x.kind(), Kind::Num(_)) && imaginary_power_of_imaginary(b, x) {
        return Some(if b == x { Signs::NONZERO } else { Signs::POS });
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

/// `base.is_imaginary and exp.is_imaginary and log(base).is_imaginary`.
fn imaginary_power_of_imaginary(b: &Expr, x: &Expr) -> bool {
    fact(b, "imaginary") == Some(true)
        && fact(x, "imaginary") == Some(true)
        && fact(&crate::functions::func(Func::Log, vec![b.clone()]), "imaginary") == Some(true)
}

/// The signs of a product of factors without real signs: an even number of
/// imaginary factors (finite and nonzero) gives a nonzero real value, as
/// `Mul._eval_real_imag` does (`I*i` is real for an imaginary `i`).
fn imaginary_factors_signs(facs: &[(Expr, Expr)]) -> Option<Signs> {
    let mut known = Signs::POS;
    let mut imaginary = 0usize;
    for (b, x) in facs {
        match pow_signs(b, x) {
            Some(s) => known = known.mul(s),
            None => {
                let factor = crate::ops::mul_factor(b, x);
                if fact(&factor, "imaginary") != Some(true) {
                    return None;
                }
                imaginary += 1;
            }
        }
    }
    (imaginary % 2 == 0).then(|| known.mul(Signs::NONZERO))
}

fn func_signs(f: &Func, args: &[Expr]) -> Option<Signs> {
    let a = args.first()?;
    match f {
        Func::Product => crate::products::signs(args),
        Func::Exp => {
            real_signs(a)?;
            Some(Signs::POS)
        }
        // exp_polar of a real number evaluates to a positive number; SymPy
        // leaves the sign of a symbolic exp_polar unknown.
        Func::ExpPolar if a.is_number() && is_real(a) == Some(true) => Some(Signs::POS),
        // exp_polar(c + I*pi) evaluates through exp(c + I*pi) = -exp(c),
        // so the numeric sign test of the reference finds it negative.
        Func::ExpPolar
            if crate::functions::is_numeric_value(a)
                && is_real(&crate::ops::sub(a, &crate::ops::mul2(&Expr::i(), &Expr::pi()))) == Some(true) =>
        {
            Some(Signs::NEG)
        }
        Func::Special(h) => {
            if *h == crate::special::SpecialFn::LambertW {
                let branch = args.get(1).map_or(Some(0), Expr::as_i64)?;
                if crate::special::lambert::real(a, branch) != Some(true) { return None; }
                return Some(if branch == -1 { Signs::NEG } else {
                    match real_signs(a)? {
                        Signs::POS => Signs::POS,
                        Signs::NEG => Signs::NEG,
                        Signs::ZERO => Signs::ZERO,
                        _ => Signs::ANY,
                    }
                });
            }
            if let Some(s) = crate::special::elementary::signs(*h, args) {
                return Some(s);
            }
            let (pos, neg, nonneg) = crate::special::numtheory::sign_facts(*h, args)?;
            let mut bits = 7u8;
            if nonneg || neg == Some(false) {
                bits &= !1;
            }
            if pos == Some(false) {
                bits &= !4;
            }
            if pos == Some(true) {
                bits = 4;
            }
            if neg == Some(true) {
                bits = 1;
            }
            Some(Signs(bits))
        }
        Func::Abs => {
            if is_finite(a)!=Some(true) { return None; }
            let s = if is_zero(a) == Some(false) { Signs::POS } else { Signs::NONNEG };
            Some(s)
        }
        Func::Log => {
            let s = real_signs(a)?;
            if s != Signs::POS {
                return None;
            }
            // log._eval_is_extended_positive and _nonnegative: the sign of a
            // positive argument minus one, decided without numerical
            // refinement (log(2) is positive, log(1/3) negative).
            if let Some(shifted) = real_signs(&crate::ops::add2(a, &Expr::int(-1))) {
                if shifted != Signs::ANY {
                    return Some(shifted);
                }
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
        // These range signs include possible infinities. The finite-real
        // queries below require a separate proof for projections.
        // A numeric projection is finite, so the reference decides its sign
        // numerically (im(CRootOf(x**5 + 11*x - 2, 1)) is negative).
        Func::Re | Func::Im if crate::functions::is_numeric_value(a) => {
            numeric_refine(&Expr::from_kind(Kind::Func(f.clone(), args.to_vec())), Signs::ANY)
        }
        Func::Re | Func::Im => Some(Signs::ANY),
        Func::Arg => Some(Signs::ANY),
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
    // Heaviside.is_real = True (a class fact, whatever H0)
    if let Kind::Func(Func::Heaviside, _) = e.kind() {
        return Some(true);
    }
    if let Kind::Func(Func::Re | Func::Im | Func::Abs,args)=e.kind() {
        return is_finite(&args[0]).filter(|finite|*finite);
    }
    // exp_polar._eval_is_extended_real (finite for a finite argument)
    if let Kind::Func(Func::ExpPolar,args)=e.kind() {
        if is_real(&args[0])==Some(true){return Some(true);}
    }
    if real_signs(e).is_some() {
        return Some(true);
    }
    if declared(e).has(A::NOT_REAL){return Some(false);}
    match e.kind() {
        Kind::Pow(base,exponent) if exact::rational_power(base,exponent).and_then(|f| f.real).is_some() => exact::rational_power(base,exponent).and_then(|f| f.real),
        // Pow._eval_is_extended_real: a negative base to a non-integer
        // rational power is not real (sqrt(n), (-sqrt(2) - 1)**(1/3)).
        Kind::Pow(base,exponent) if exponent.as_number().is_some_and(|n|n.is_rational()&&!n.is_integer())
            && is_negative(base)==Some(true) => Some(false),
        Kind::Const(Constant::I) => Some(false),
        Kind::Func(Func::Special(crate::special::SpecialFn::LambertW), args) =>
            crate::special::lambert::real(&args[0], args.get(1).map_or(Some(0), Expr::as_i64)?),
        Kind::Const(Constant::NaN)=>None,
        Kind::Const(c) if c.is_infinite()=>Some(false),
        Kind::Sym(s) if s.assumptions().has(A::NOT_REAL) => Some(false),
        Kind::Mul(c, facs) if facs.len() == 1 && facs[0].0.is_constant(Constant::I) && facs[0].1.is_one() && !c.is_zero() => {
            Some(false)
        }
        Kind::Add(_,terms) => {
            let mut nonreal=0;
            let mut paired:Option<Vec<bool>>=None;
            for (index,(term,_)) in terms.iter().enumerate(){
                if paired.as_ref().is_some_and(|flags|flags[index]){continue;}
                let answer=is_real(term);
                if answer!=Some(true)&&matches!(term.kind(),Kind::Pow(..)|Kind::Mul(..)){
                    if paired.get_or_insert_with(||exact::paired_conjugate_terms(terms))[index]{continue;}
                }
                match answer{Some(true)=>{},Some(false)=>nonreal+=1,None=>return None}
            }
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
                if let Some((first,last,scale)) = numeric_samples(e) {
                    if stable_component(&first.im,&last.im,&scale) {
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
    let signed = real_signs(e).and_then(|s| {
        if s == Signs::POS {
            Some(true)
        } else if !s.has(Signs::POS) {
            Some(false)
        } else {
            None
        }
    });
    // `Expr._eval_is_positive`: False when the numerical sign excludes it,
    // whatever the finiteness (`polylog(3, 1/2).is_negative`).
    signed.or_else(|| match reference_numeric_sign(e) {
        Some(NumericSign::Negative) => Some(false),
        _ => None,
    })
}

#[inline]
pub fn is_negative(e:&Expr)->Option<bool>{cached_bool(e,Query::Negative,is_negative_uncached)}
fn is_negative_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::NOT_NEGATIVE) || exact_queries(||is_real(e))==Some(false){return Some(false);}
    let signed = real_signs(e).and_then(|s| {
        if s == Signs::NEG {
            Some(true)
        } else if !s.has(Signs::NEG) {
            Some(false)
        } else {
            None
        }
    });
    // `Expr._eval_is_negative`: False when the numerical sign excludes it,
    // whatever the finiteness (`polylog(3, 1/2).is_negative`).
    signed.or_else(|| match reference_numeric_sign(e) {
        Some(NumericSign::Positive) => Some(false),
        _ => None,
    })
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
    // `zero -> commutative`: `(A*B - B*A).is_zero` is False in SymPy.
    if nc_settles(e){return Some(false);}
    match e.kind() {
        Kind::Num(n) => return Some(n.is_zero()),
        Kind::Sym(s) if s.assumptions().has(A::NOT_ZERO) => return Some(false),
        Kind::Sym(s) if s.assumptions().has(A::ZERO) => return Some(true),
        Kind::Const(Constant::I) => return Some(false),
        Kind::Const(c) if c.is_number_symbol() => return Some(false),
        Kind::Const(c) if c.is_infinite()=>return Some(false),
        // ExpBase._eval_is_zero
        Kind::Func(Func::ExpPolar, a) => return Some(a[0].is_constant(Constant::NegInfinity)),
        // SymPy's exp._eval_is_zero: only exp(-oo) vanishes.
        Kind::Func(Func::Exp,args)=>return Some(matches!(args.first().map(Expr::kind),Some(Kind::Const(Constant::NegInfinity)))),
        // log._eval_is_zero: (self.args[0] - 1).is_zero
        Kind::Func(Func::Log,args) if args.len()==1=>{
            if let Some(answer)=is_zero(&crate::ops::add2(&args[0],&Expr::int(-1))){return Some(answer);}
        },
        Kind::Mul(c,factors) if !c.is_zero() && factors.iter().all(|(b,x)| exact_queries(||is_zero(&crate::ops::mul_factor(b,x)))==Some(false)) => return Some(false),
        Kind::Pow(base,exponent) => {
            if let Some(facts)=exact::rational_power(base,exponent) {if facts.zero.is_some(){return facts.zero;}}
            if exact_queries(||is_zero(base))==Some(false) && is_finite(exponent)==Some(true) {return Some(false);}
        },
        // re._eval_is_zero (arg imaginary or zero), im._eval_is_zero (arg
        // extended real), Abs/sign (arg zero): all False for an argument
        // that does not commute.
        Kind::Func(Func::Re|Func::Im|Func::Abs|Func::Sign,args) if args.len()==1&&noncommuting(&args[0]) => return Some(false),
        Kind::Add(..) if e.has_noncommutative_host() && crate::functions::noncommutative_sum_nonzero(e)==Some(false) => return Some(false),
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
        if let Some((first,last,scale)) = numeric_samples(e) {
            if stable_component(&first.re,&last.re,&scale) || stable_component(&first.im,&last.im,&scale) {
                return Some(false);
            }
        }
    }
    // Add._eval_is_zero last: the sign and numeric paths above share the
    // DAG cache and settle the real sums first.
    // Sums of real terms are left to real_signs: the generic extended_real
    // fact is not cached and would revisit a shared DAG exponentially.
    if let Kind::Add(_,terms)=e.kind() {
        if integer_offset_nonzero(e) {return Some(false);}
        if terms.iter().all(|(t,_)|is_real(t)==Some(true)) {return None;}
        return add_is_zero(e);
    }
    None
}

/// `zero -> even -> integer` in the reference fact rules: a sum of integer
/// terms and a rational constant is nonzero when the constant is not an
/// integer (`n - 1/3`), or when the sum is odd (`2*n - 1`).
fn integer_offset_nonzero(e:&Expr)->bool{
    let args=e.args();
    let mut constant=None;
    let mut rest=Vec::with_capacity(args.len());
    for a in &args {
        match a.as_number() {
            Some(c) if constant.is_none() => constant=Some(c.clone()),
            Some(_) => return false,
            None => rest.push(a),
        }
    }
    let Some(c)=constant else {return false};
    if !c.is_rational() || rest.is_empty() || !rest.iter().all(|t|is_integer(t)==Some(true)) {return false;}
    !c.is_integer() || is_odd(e)==Some(true)
}

/// SymPy's `Add._eval_is_zero`: the real terms must sum to a nonzero
/// value next to imaginary ones (`1 + I`, `p + I`, `1 + I*r` are nonzero),
/// or all vanish. Only its decided answers are used.
fn add_is_zero(e:&Expr)->Option<bool>{
    let args=e.args();
    let (mut nonzero,mut zeros,mut im_or_z,mut imaginary,mut unknown)=(Vec::new(),0usize,false,0usize,false);
    for a in &args {
        if fact(a,"extended_real")==Some(true) {
            match is_zero(a) {
                Some(true)=>zeros+=1,
                Some(false)=>nonzero.push(a.clone()),
                None=>unknown=true,
            }
        } else if fact(a,"imaginary")==Some(true) {
            imaginary+=1;
        } else if let Kind::Mul(..)=a.kind() {
            // `a.as_coeff_mul(I)`: exactly the factor I times an extended
            // real coefficient
            let factors=a.args();
            if factors.iter().filter(|f|f.is_constant(Constant::I)).count()!=1 {return None;}
            let rest:Vec<Expr>=factors.into_iter().filter(|f|!f.is_constant(Constant::I)).collect();
            if rest.iter().any(|f|crate::traverse::has_structural(f,&Expr::i())) {return None;}
            if fact(&crate::ops::mul(&rest),"extended_real")!=Some(true) {return None;}
            im_or_z=true;
        } else {
            return None;
        }
    }
    // extended real terms next to a single imaginary one (nonzero, as
    // `imaginary` excludes 0): the sum is not extended real, so not zero
    // (`zero -> extended_real`, how the reference settles `(r + I).is_zero`)
    if imaginary==1 && !im_or_z {return Some(false);}
    if unknown {return None;}
    if zeros==args.len() {return Some(true);}
    if nonzero.is_empty()||nonzero.len()==args.len() {return None;}
    match is_zero(&crate::ops::add(&nonzero)) {
        Some(true) if !im_or_z && imaginary==0 => Some(true),
        Some(true) if !im_or_z && imaginary==1 => Some(false),
        Some(false) => Some(false),
        _ => None,
    }
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
        // re/im._eval_is_algebraic is that of the argument, False when it
        // does not commute, and integer implies algebraic.
        Kind::Func(Func::Re|Func::Im,args) if args.len()==1&&noncommuting(&args[0]) => Some(false),
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
        Kind::Const(Constant::NaN|Constant::EulerGamma|Constant::Catalan)=>None,
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
        Kind::Func(Func::Special(special), args) => crate::special::elementary::integer(*special, args)
            .or_else(|| crate::special::numtheory::integer_fact(*special, args)),
        _ => None,
    }
}

#[inline]
pub fn is_rational(e:&Expr)->Option<bool>{cached_bool(e,Query::Rational,is_rational_uncached)}
fn is_rational_uncached(e: &Expr) -> Option<bool> {
    if declared(e).has(A::RATIONAL){return Some(true);}
    if declared(e).has(A::NOT_RATIONAL)||declared(e).has(A::IRRATIONAL){return Some(false);}
    if exact_queries(||is_real(e))==Some(false){return Some(false);}
    if let Kind::Func(Func::Special(h),args)=e.kind(){
        if crate::special::numtheory::integer_fact(*h,args)==Some(true){return Some(true);}
    }
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
        Kind::Const(Constant::Pi | Constant::E | Constant::I | Constant::GoldenRatio | Constant::TribonacciConstant) => Some(false),
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
    if let Some(answer)=crate::host::fact(e,"even"){return answer;}
    if let Some(symbol)=e.as_sym(){return symbol.facts().get("even");}
    if let Kind::Func(Func::Undefined(head),_)=e.kind(){if !head.internal{return head.facts().get("even");}}
    if is_integer(e)==Some(false){return Some(false);}
    parity_even(e)
}

pub fn is_odd(e: &Expr) -> Option<bool> {
    if let Some(answer)=crate::host::fact(e,"odd"){return answer;}
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
        Kind::Func(Func::Special(special), args) if crate::special::elementary::finite(*special, args).is_some() => {
            crate::special::elementary::finite(*special, args)
        }
        // sign takes values in {0} and the unit circle
        Kind::Func(Func::Sign, _) => Some(true),
        // Heaviside.is_finite = True (a class fact)
        Kind::Func(Func::Heaviside, _) => Some(true),
        // ExpBase._eval_is_finite
        Kind::Func(Func::ExpPolar | Func::Exp, a) => {
            if is_finite(&a[0]) == Some(true) || a[0].is_constant(Constant::NegInfinity) { Some(true) }
            else if a[0].is_constant(Constant::Infinity) { Some(false) }
            else { None }
        }
        Kind::Func(Func::Re | Func::Im | Func::Abs, args) => is_finite(&args[0]).filter(|finite|*finite),
        // Sums and products of finite terms are finite (SymPy's
        // ``Add._eval_is_finite`` and ``Mul._eval_is_finite``), whether or
        // not they are real: 1 + I is finite, hence not extended real.
        Kind::Add(_, terms) if terms.iter().all(|(term, _)| is_finite(term) == Some(true)) => Some(true),
        // ``Add._eval_is_finite``: a single infinite term among finite ones.
        Kind::Add(..) if fuzzy_group(e.args().iter().map(is_finite),true)==Some(false) => Some(false),
        // ``Mul._eval_is_infinite``: an infinite factor and no zero factor.
        Kind::Mul(..) if { let args=e.args(); args.iter().any(|a|is_finite(a)==Some(false))&&args.iter().all(|a|is_zero(a)==Some(false)) } => Some(false),
        Kind::Mul(_, factors) if factors.iter().all(|(base, exponent)| finite_power(base, exponent)) => Some(true),
        Kind::Pow(base, exponent) if finite_power(base, exponent) => Some(true),
        _ => {
            if real_signs(e).is_some() {
                Some(true)
            } else {
                None
            }
        }
    }
}

/// A power of a finite base to a finite exponent that is nonnegative, or of
/// a base that is not zero (``fuzzy_not(base.is_zero)``, so complex bases
/// such as I count), is finite (SymPy's ``Pow._eval_is_finite``).
fn finite_power(base: &Expr, exponent: &Expr) -> bool {
    is_finite(base) == Some(true) && is_finite(exponent) == Some(true)
        && (is_nonnegative(exponent) == Some(true) || is_zero(base) == Some(false))
}

#[cfg(test)]
mod numeric_literal_tests {
    use super::*;
    use crate::bigfloat::BigFloat;
    use num_bigint::BigInt;

    #[test]
    fn dyadic_integer_queries_never_materialize_extreme_exponents(){
        for exp in [1_000_000_000,i64::MAX]{
            let value=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(3),exp:exp.into()},200));
            assert_eq!(is_integer(&value),None);
            assert_eq!(is_even(&value),None);
            assert_eq!(is_odd(&value),None);
            assert_eq!(fact(&value,"noninteger"),None);
            assert_eq!(is_rational(&value),None);
        }
        for exp in [-1_000_000_000,i64::MIN]{
            let value=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(3),exp:exp.into()},200));
            assert_eq!(is_integer(&value),Some(false));
            assert_eq!(is_even(&value),Some(false));
            assert_eq!(is_odd(&value),Some(false));
            assert_eq!(fact(&value,"noninteger"),Some(true));
        }
        let value=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(4),exp:(-2i64).into()},200));
        assert_eq!(is_integer(&value),None);
        let zero=Expr::num(Number::MpFloat(BigFloat{man:BigInt::from(0),exp:i64::MIN.into()},200));
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
mod shared_graph_tests {
    use super::*;
    use crate::functions::func;
    use crate::ops::{add2, sub};

    /// ``e_{k+1} = exp(e_k − y) + sin(e_k) + x`` : each level reuses ``e_k``
    /// twice, so the expanded tree doubles at each level.
    fn shared_chain(depth: usize) -> (Expr, Expr) {
        let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));
        let mut e = x.clone();
        for _ in 0..depth {
            let a = func(Func::Exp, vec![sub(&e, &y)]);
            e = add2(&add2(&a, &func(Func::Sin, vec![e.clone()])), &x);
        }
        (e, y)
    }

    #[test]
    fn sign_queries_on_deep_shared_sums_are_not_tree_walks() {
        // Before: reference_numeric_sign walked the expanded tree of every
        // sum (2^60 nodes here) before its cached symbolic-number test.
        let (e, y) = shared_chain(60);
        let d = sub(&e, &y);
        assert_eq!(is_positive(&d), None);
        assert_eq!(is_zero(&d), None);
        let _ = func(Func::Abs, vec![d.clone()]);
        let _ = crate::functions::min(vec![e.clone(), y.clone()]);
        let _ = crate::relational::lt(&e, &y);
    }

    #[test]
    fn facts_of_deep_shared_real_graphs_are_not_tree_walks() {
        // ``e_{k+1} = exp(min(e_k, y) − e_k)/2 + sin(e_k) + x`` over real
        // symbols: `Min` and the Python wrapper query `extended_real`, and
        // `structural_fact` asked it of every argument (×6 per level).
        let real = crate::symbol::Assumptions::REAL;
        let (x, y) = (Expr::symbol_with("x", real), Expr::symbol_with("y", real));
        let mut e = x.clone();
        for _ in 0..60 {
            let m = crate::functions::min(vec![e.clone(), y.clone()]);
            let a = func(Func::Exp, vec![sub(&m, &e)]);
            e = add2(&add2(&crate::ops::mul(&[Expr::rational(1, 2), a]), &func(Func::Sin, vec![e.clone()])), &x);
        }
        assert_eq!(fact(&e, "extended_real"), Some(true));
        assert_eq!(fact(&sub(&e, &y), "extended_real"), Some(true));
    }

    #[test]
    fn numeric_sign_of_sums_with_functions_is_unchanged() {
        // a closed sum with a function value still takes the numerical sign
        let e = crate::parser::parse_str("sin(1) - 1/2").unwrap();
        assert_eq!(is_positive(&e), Some(true));
        let e = crate::parser::parse_str("exp(-3) - 1/20").unwrap();
        assert_eq!(is_negative(&e), Some(true));
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn cancelling_algebraic_numbers_are_zero_by_minimal_polynomial(){
        // Expr._eval_is_extended_positive_negative: no numerical
        // significance and a minimal polynomial equal to the generator.
        for text in ["-3*sqrt(93)/62 - 1/(3*sqrt(93)/2 + 31/2) + 1/2","sqrt(2)+sqrt(3)-sqrt(5+2*sqrt(6))"]{
            let e=crate::parser::parse_str(text).unwrap();
            assert_eq!(is_zero(&e),Some(true),"{text}");
            assert_eq!(is_positive(&e),Some(false),"{text}");
        }
        let e=crate::parser::parse_str("sqrt(2)+sqrt(3)-sqrt(5+2*sqrt(6))+1").unwrap();
        assert_eq!(is_zero(&e),Some(false));
    }

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
        // The identical terms cancel: zero, as the minimal polynomial of the
        // reference proves on a small graph.
        assert_eq!(real_signs(&e),Some(Signs::ZERO));
        assert!(crate::functions::is_numeric_value(&e));assert_eq!(e.cached_assumption(60,2),3);
        assert_eq!(crate::evalf::eval_f64(&e),Some(0.0));
        assert_eq!(exact_queries(||real_signs(&e)),Some(Signs::ANY));
    }

    #[test]
    fn scaled_cancellation_does_not_acquire_false_numeric_facts(){
        let value=crate::parser::parse_str("10**100*(sin(3)**2+cos(3)**2-1)").unwrap();
        assert_eq!(is_zero(&value),None);
        assert_eq!(is_positive(&value),None);
        assert_eq!(is_negative(&value),None);
        assert_eq!(exact_queries(||is_zero(&value)),None);
    }

    #[test]
    fn tiny_numeric_sign_refinement_does_not_underflow(){
        let value=crate::parser::parse_str("sin(1/10**400)").unwrap();
        assert_eq!(is_positive(&value),Some(true));
        assert_eq!(is_zero(&value),Some(false));
        assert_eq!(is_negative(&value),Some(false));
    }
}
