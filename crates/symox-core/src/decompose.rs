//! Structural expression decomposition without polynomial expansion.
use rustc_hash::{FxHashMap,FxHashSet};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One,ToPrimitive,Signed};
use crate::expr::{Constant,Expr,Kind};
use crate::number::Number;
use crate::{ops,traverse};

fn args(e:&Expr)->Vec<Expr>{traverse::ordered_args(e)}

/// Leading Number coefficient, preserving the unevaluated tail of a product.
pub fn as_coeff_mul_head(e:&Expr,rational:bool)->(Expr,Expr){
    fn product(parent:&Expr,values:&[Expr])->Expr{
        if parent.is_canonical(){return ops::mul(values);}
        match values.len(){0=>Expr::one(),1=>values[0].clone(),_=>Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,values.iter().cloned().map(|a|(a,Expr::one())).collect()))}
    }
    if number_atom(e){
        if !rational||e.is_rational_number(){return(e.clone(),Expr::one());}
        return(Expr::one(),e.clone());
    }
    if e.is_mul(){
        let values=args(e);
        if let Some(first)=values.first().filter(|a|number_atom(a)){
            if !rational||first.is_rational_number(){return(first.clone(),product(e,&values[1..]));}
            if first.as_number().is_some_and(Number::is_negative)||first.is_constant(Constant::NegInfinity){
                let mut tail=vec![ops::neg(first)];tail.extend(values[1..].iter().cloned());
                return(Expr::minus_one(),product(e,&tail));
            }
        }
    }
    (Expr::one(),e.clone())
}

/// Leading Number coefficient of a sum, retaining the public raw tail.
pub fn as_coeff_add_head(e:&Expr,rational:bool)->(Expr,Expr){
    if number_atom(e){
        return if !rational||e.is_rational_number(){(e.clone(),Expr::zero())}
            else{(Expr::zero(),e.clone())};
    }
    if let Kind::Add(constant,terms)=e.kind(){
        let values=args(e);
        if let Some(first)=values.first().filter(|a|number_atom(a)&&(!rational||a.is_rational_number())){
            let tail=&values[1..];
            let tail=match tail.len(){
                0=>Expr::zero(),1=>tail[0].clone(),
                _ if e.is_canonical()&&first.as_number()==Some(constant)=>Expr::from_kind(Kind::Add(Number::ZERO,terms.clone())),
                _=>Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,tail.iter().cloned().map(|e|(e,Number::ONE)).collect())),
            };
            return(first.clone(),tail);
        }
    }
    (Expr::zero(),e.clone())
}

/// Public power view. Internal multiplication keeps the literal base so that
/// this presentation rule cannot change automatic evaluation.
pub fn as_base_exp(e:&Expr)->(Expr,Expr){
    match e.kind(){
        Kind::Const(Constant::I)=>(Expr::minus_one(),Expr::rational(1,2)),
        Kind::Func(crate::expr::Func::Exp,a)=>(Expr::e(),a[0].clone()),
        Kind::Pow(base,exponent)=>{
            if let Some((numerator,denominator))=base.as_number().filter(|n|n.is_rational()).and_then(Number::numer_denom){
                if numerator.is_one() && !denominator.is_one(){return(Expr::num(Number::from_bigint(denominator)),ops::neg(exponent));}
            }
            (base.clone(),exponent.clone())
        },
        _=>(e.clone(),Expr::one()),
    }
}

/// Insertion-ordered factor powers, retaining zero exponents of held products.
/// Memoization counts repeated DAG edges without expanding their occurrences.
pub fn as_powers_dict(e:&Expr)->Vec<(Expr,Expr)>{
    fn visit(e:&Expr,cache:&mut FxHashMap<usize,(Expr,Vec<(Expr,Expr)>)>)->Vec<(Expr,Expr)>{
        if let Some((_,result))=cache.get(&e.dag_id()){return result.clone();}
        let result=if e.is_mul(){
            let mut positions=FxHashMap::default();let mut out:Vec<(Expr,Expr)>=Vec::new();
            for factor in args(e){for(base,exponent)in visit(&factor,cache){
                if let Some(&index)=positions.get(&base){let pair:&mut(Expr,Expr)=&mut out[index];pair.1=ops::add2(&pair.1,&exponent);}
                else{positions.insert(base.clone(),out.len());out.push((base,exponent));}
            }}
            out
        }else if e.is_constant(Constant::NegInfinity){vec![(Expr::minus_one(),Expr::one()),(Expr::oo(),Expr::one())]}
        else{vec![as_base_exp(e)]};
        cache.insert(e.dag_id(),(e.clone(),result.clone()));result
    }
    visit(e,&mut FxHashMap::default())
}

/// `is_commutative` is True. A node free of noncommuting symbols, functions
/// and host classes commutes without a traversal.
#[inline]
pub fn commutes(e:&Expr)->bool{!e.has_noncommutative_host()||commutativity(e)==Some(true)}

/// True when `e` contains a noncommutative Symbol or a noncommutative
/// undefined function application.  The host flag keeps commutative
/// expressions on a constant-time path.
pub fn has_noncommutative_atoms(e:&Expr)->bool{
    fn visit(e:&Expr,seen:&mut rustc_hash::FxHashSet<usize>)->bool{
        if !e.has_noncommutative_host()||!seen.insert(e.dag_id()){return false;}
        match e.kind(){
            Kind::Sym(symbol) if symbol.facts().get("commutative")==Some(false)=>true,
            Kind::Func(crate::expr::Func::Undefined(head),_) if !head.internal&&head.facts().get("commutative")==Some(false)=>true,
            _=>e.args().iter().any(|a|visit(a,seen)),
        }
    }
    e.has_noncommutative_host()&&visit(e,&mut Default::default())
}

/// Native commutativity facts used by the Python expression protocol.
pub fn commutativity(e:&Expr)->Option<bool>{
    fn visit(e:&Expr,cache:&mut FxHashMap<usize,(Expr,Option<bool>)>)->Option<bool>{
        if let Some((_,value))=cache.get(&e.dag_id()){return *value;}
        if let Kind::Func(crate::expr::Func::SymbolicDerivative,args)=e.kind(){
            return args.first().and_then(|body|visit(body,cache));
        }
        if e.has_noncommutative_host(){
            let value=host_commutativity(e,cache);
            cache.insert(e.dag_id(),(e.clone(),value));return value;
        }
        let value=match e.kind(){
            Kind::Func(crate::expr::Func::Limit,_)=>None,
            Kind::Func(crate::expr::Func::MatrixDotProduct(..),_)=>None,
            Kind::Func(crate::expr::Func::Undefined(head),_) if head.is_internal("UnevaluatedExpr")||head.is_internal("Tuple")=>None,
            Kind::Func(crate::expr::Func::Undefined(head),arguments) if !head.internal=>{
                head.facts().get("commutative").or_else(||{
                    let mut unknown=false;
                    for argument in arguments {match visit(argument,cache){Some(false)=>return Some(false),None=>unknown=true,_=>{}}}
                    if unknown{None}else{Some(true)}
                })
            },
            Kind::Add(..)|Kind::Mul(..)=>Some(args(e).iter().all(|a|visit(a,cache)==Some(true))),
            Kind::Pow(base,_)=>visit(base,cache),
            // `ExprCondPair.is_commutative` is that of its expression.
            Kind::Piecewise(pieces)=>{
                let mut unknown=false;
                for (piece,_) in pieces {match visit(piece,cache){Some(false)=>return Some(false),None=>unknown=true,_=>{}}}
                if unknown{None}else{Some(true)}
            },
            _=>Some(true),
        };
        cache.insert(e.dag_id(),(e.clone(),value));value
    }
    /// SymPy's `is_commutative` below a host node that may not commute:
    /// the class declaration or handler of a host node, `fuzzy_and` over
    /// the arguments of functions (`Function._eval_is_commutative`) and the
    /// base and exponent of powers, and True for Add and Mul only when all
    /// their arguments commute (`AssocOp`).
    fn host_commutativity(e:&Expr,cache:&mut FxHashMap<usize,(Expr,Option<bool>)>)->Option<bool>{
        let all=|values:&mut dyn Iterator<Item=Option<bool>>|{
            let mut unknown=false;
            for value in values {match value{Some(false)=>return Some(false),None=>unknown=true,_=>{}}}
            if unknown{None}else{Some(true)}
        };
        if let Some(host)=e.host_handle(){
            if let Some(answer)=host.commutative(e){return answer;}
            if let Some(symbol)=e.as_sym(){return symbol.facts().get("commutative").or(Some(true));}
        }
        match e.kind(){
            Kind::Sym(symbol)=>symbol.facts().get("commutative").or(Some(true)),
            Kind::Func(crate::expr::Func::Undefined(head),arguments) if !head.internal=>head.facts().get("commutative")
                .or_else(||all(&mut arguments.iter().map(|argument|visit(argument,cache)))),
            Kind::Add(..)|Kind::Mul(..)=>Some(args(e).iter().all(|a|visit(a,cache)==Some(true))),
            Kind::Pow(base,exponent)=>all(&mut [base,exponent].into_iter().map(|part|visit(part,cache))),
            Kind::Func(crate::expr::Func::Limit|crate::expr::Func::MatrixDotProduct(..),_)=>None,
            Kind::Func(crate::expr::Func::Undefined(head),_) if head.is_internal("UnevaluatedExpr")||head.is_internal("Tuple")=>None,
            Kind::Func(head,_) if commutative_head(head)=>Some(true),
            Kind::Func(_,arguments)=>all(&mut arguments.iter().map(|argument|visit(argument,cache))),
            Kind::Piecewise(pieces)=>all(&mut pieces.iter().map(|(piece,_)|visit(piece,cache))),
            _=>Some(true),
        }
    }
    visit(e,&mut FxHashMap::default())
}

/// Heads whose SymPy class declares a fact that implies `commutative`
/// (`extended_real` for re, im, Abs and arg, `complex` for sign, `real`
/// for DiracDelta and Heaviside, `integer` for the number-theoretic and
/// Kronecker/Levi-Civita symbols, `is_commutative = True` on
/// `TupleParametersBase` for hyper and meijerg). Their applications commute
/// whatever the arguments.
pub fn commutative_head(head:&crate::expr::Func)->bool{
    use crate::expr::Func;
    use crate::special::SpecialFn as S;
    match head{
        Func::Re|Func::Im|Func::Abs|Func::Sign|Func::Arg|Func::DiracDelta|Func::Heaviside|Func::Max|Func::Min=>true,
        Func::Special(h)=>matches!(h,S::Partition|S::Totient|S::ReducedTotient|S::Mobius|S::DivisorSigma|S::UdivisorSigma
            |S::LegendreSymbol|S::JacobiSymbol|S::KroneckerSymbol|S::PrimeNu|S::PrimeOmega|S::PrimePi
            |S::Hyper|S::MeijerG|S::KroneckerDelta|S::LeviCivita),
        _=>false,
    }
}

/// Split at the first factor whose commutativity is not known to be true.
pub fn args_cnc(e:&Expr,split_one:bool)->(Vec<Expr>,Vec<Expr>){
    let factors=if e.is_mul(){args(e)}else{vec![e.clone()]};
    let boundary=factors.iter().position(|factor|commutativity(factor)!=Some(true)).unwrap_or(factors.len());
    let(mut commutative,noncommutative)=(factors[..boundary].to_vec(),factors[boundary..].to_vec());
    if split_one{
        if let Some(first)=commutative.first().filter(|e|!e.is_minus_one() && (e.as_number().is_some_and(Number::is_negative)||e.is_constant(Constant::NegInfinity))){
            let positive=ops::neg(first);commutative[0]=positive;commutative.insert(0,Expr::minus_one());
        }
    }
    (commutative,noncommutative)
}

/// A syntactic sign choice; assumptions about a symbol's sign do not count.
pub fn could_extract_minus_sign(e:&Expr)->bool{
    fn negate(e:&Expr,cache:&mut FxHashMap<usize,(Expr,Expr)>)->Expr{
        if let Some((_,value))=cache.get(&e.dag_id()){return value.clone();}
        let value=ops::neg(e);
        cache.insert(e.dag_id(),(e.clone(),value.clone()));value
    }
    fn visit(e:&Expr,cache:&mut FxHashMap<usize,(Expr,bool)>,negative_cache:&mut FxHashMap<usize,(Expr,Expr)>)->bool{
        if let Some((_,value))=cache.get(&e.dag_id()){return *value;}
        let value=match e.kind(){
            Kind::Num(n)=>n.is_negative(),
            Kind::Const(Constant::NegInfinity)=>true,
            Kind::Mul(..)=>{
                let first=args(e).into_iter().next();
                first.is_some_and(|a|a.as_number().is_some_and(Number::is_negative)||a.is_constant(Constant::NegInfinity)) && *e!=negate(e,negative_cache)
            },
            Kind::Add(..)=>{
                let terms=args(e);let negative=terms.iter().filter(|term|visit(term,cache,negative_cache)).count();
                match negative.cmp(&(terms.len()-negative)){
                    std::cmp::Ordering::Greater=>true,
                    std::cmp::Ordering::Less=>false,
                    std::cmp::Ordering::Equal=>{
                        crate::sortkey::cmp_default(e,&negate(e,negative_cache)).is_lt()
                    },
                }
            },
            _=>false,
        };
        cache.insert(e.dag_id(),(e.clone(),value));value
    }
    visit(e,&mut FxHashMap::default(),&mut FxHashMap::default())
}

/// The sign choice of the odd/even function rewrites `f(-z) -> -f(z)`.
///
/// Negating a noncommutative sum `z` keeps `-z` as `Mul(-1, z)`, and SymPy
/// reports such a product as negative whatever `z` is.  When `z` itself could
/// extract a minus sign, SymPy then rewrites `f(z)` to `-f(-z)` and back
/// forever (a RecursionError).  Here `-z` is negative exactly when `z` is not,
/// so every rewrite terminates; all other cases are SymPy's rule unchanged.
pub fn could_extract_minus_sign_eval(e:&Expr)->bool{
    if let Kind::Mul(c,factors)=e.kind(){
        if c.is_minus_one() && factors.len()==1 && factors[0].1.is_one() && factors[0].0.is_add() && !commutes(&factors[0].0){
            return !could_extract_minus_sign(&factors[0].0);
        }
    }
    could_extract_minus_sign(e)
}

fn number_atom(e:&Expr)->bool{matches!(e.kind(),Kind::Num(_)|Kind::Const(Constant::Infinity|Constant::NegInfinity|Constant::NaN))}
fn has(e:&Expr,pattern:&Expr)->bool {
    if pattern.is_symbol() || e.has_matrix() || pattern.has_matrix(){return traverse::has(e,pattern);}
    fn visit(e:&Expr,p:&Expr,seen:&mut FxHashMap<usize,Expr>)->bool {
        if e==p{return true;}
        if seen.insert(e.dag_id(),e.clone()).is_some(){return false;}
        if (e.is_add() && p.is_add()) || (e.is_mul() && p.is_mul()){
            let values=args(e);
            if args(p).iter().all(|term|values.contains(term)){return true;}
        }
        e.args().iter().any(|a|visit(a,p,seen))
    }
    visit(e,pattern,&mut FxHashMap::default())
}
fn has_any(e:&Expr,deps:&[Expr])->bool {deps.iter().any(|p|has(e,p))}

/// Rebuild a chosen subset of public arguments without distributing numbers
/// over sums. Ordinary canonical trees retain their compact representation.
fn join(values:&[Expr],add:bool)->Expr {
    if values.is_empty(){return if add{Expr::zero()}else{Expr::one()};}
    if values.len()==1{return values[0].clone();}
    let canonical=if add{ops::add(values)}else{ops::mul(values)};
    let actual=args(&canonical);
    if actual.len()==values.len() && actual.iter().zip(values).all(|(a,b)|a.same_presentation(b)){return canonical;}
    if add {Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,values.iter().cloned().map(|e|(e,Number::ONE)).collect()))}
    else {
        let (coefficient,rest)=match values[0].as_number(){Some(n)=>(n.clone(),&values[1..]),None=>(Number::ONE,values)};
        Expr::from_unevaluated_kind(Kind::Mul(coefficient,rest.iter().cloned().map(|e|(e,Expr::one())).collect()))
    }
}

pub fn as_coeff_add(e:&Expr,deps:&[Expr])->(Expr,Vec<Expr>){
    if e.is_rational_number(){return(e.clone(),Vec::new());}
    if number_atom(e){return(Expr::zero(),vec![e.clone()]);}
    if e.is_add(){
        let values=args(e);
        if !deps.is_empty(){let (dependent,independent):(Vec<_>,Vec<_>)=values.into_iter().partition(|v|has_any(v,deps));return(join(&independent,true),dependent);}
        let (coefficient,mut tail)=as_coeff_add(&values[0],&[]);
        if !coefficient.is_zero(){tail.extend(values[1..].iter().cloned());return(coefficient,tail);}
        return(Expr::zero(),values);
    }
    if !deps.is_empty() && !has_any(e,deps){(e.clone(),Vec::new())}else{(Expr::zero(),vec![e.clone()])}
}

pub fn as_coeff_mul(e:&Expr,deps:&[Expr],rational:bool)->(Expr,Vec<Expr>){
    if number_atom(e){
        if e.is_rational_number() || !rational{return(e.clone(),Vec::new());}
        let negative=e.as_number().is_some_and(Number::is_negative);
        return if negative{(Expr::minus_one(),vec![ops::neg(e)])}else{(Expr::one(),vec![e.clone()])};
    }
    if e.is_mul(){
        let values=args(e);
        if !deps.is_empty(){let (dependent,independent):(Vec<_>,Vec<_>)=values.into_iter().partition(|v|has_any(v,deps));return(join(&independent,false),dependent);}
        if number_atom(&values[0]){
            let (coefficient,mut tail)=as_coeff_mul(&values[0],&[],rational);
            tail.extend(values[1..].iter().cloned());return(coefficient,tail);
        }
        return(Expr::one(),values);
    }
    if !deps.is_empty() && !has_any(e,deps){(e.clone(),Vec::new())}else{(Expr::one(),vec![e.clone()])}
}

pub fn as_independent(e:&Expr,deps:&[Expr],as_add:Option<bool>)->(Expr,Expr){
    if e.is_rational_number() && e.is_zero(){return(Expr::zero(),Expr::zero());}
    let additive=as_add.unwrap_or_else(||e.is_add());
    let identity=if additive{Expr::zero()}else{Expr::one()};
    let symbol_deps:FxHashSet<_>=deps.iter().filter(|p|p.is_symbol()).cloned().collect();
    let other:Vec<_>=deps.iter().filter(|p|!p.is_symbol()).cloned().collect();
    let dependent=|term:&Expr|has_any(term,&other)||(!symbol_deps.is_empty() && traverse::free_symbols(term).iter().any(|s|symbol_deps.contains(s)));
    if (additive && !e.is_add()) || (!additive && !e.is_mul()){
        return if dependent(e){(identity,e.clone())}else{(e.clone(),identity)};
    }
    if !additive && e.has_noncommutative_host(){
        // SymPy: once a noncommutative factor depends on deps, so does the
        // remainder of the noncommutative part.
        let(commutative,noncommutative)=args_cnc(e,true);
        let(mut depend,mut independent):(Vec<_>,Vec<_>)=commutative.into_iter().partition(|f|dependent(f));
        if let Some(start)=noncommutative.iter().position(|f|dependent(f)){
            independent.extend(noncommutative[..start].iter().cloned());
            depend.extend(noncommutative[start..].iter().cloned());
        }else{independent.extend(noncommutative);}
        return(ops::mul(&independent),join(&depend,false));
    }
    let (dependent,independent):(Vec<_>,Vec<_>)=args(e).into_iter().partition(dependent);
    // Unlike the coefficient APIs' raw-argument rebuild, SymPy evaluates the
    // independent part with the ordinary Add/Mul constructor.
    let independent=if additive{ops::add(&independent)}else{ops::mul(&independent)};
    (independent,join(&dependent,additive))
}

/// Preserve a coefficient on a sum; unlike general arithmetic this operation
/// is used to keep a decomposition factored, so it must not distribute it.
fn keep(c:&Number,p:&Expr)->Expr {
    if c.is_one(){return p.clone();}
    if p.is_one(){return Expr::num(c.clone());}
    if p.is_add(){return Expr::from_unevaluated_kind(Kind::Mul(c.clone(),vec![(p.clone(),Expr::one())]));}
    ops::mul(&[Expr::num(c.clone()),p.clone()])
}

fn check_power_size(base:&Number,exponent:&Number)->Result<(),String>{
    if base.is_zero() || base.is_one() || base.is_minus_one(){return Ok(());}
    if let (Some((n,d)),Some((p,q)))=(base.numer_denom(),exponent.numer_denom()){
        if p==BigInt::from(0){return Ok(());}
        let magnitude=(p.abs()/q+BigInt::one()).to_u64();
        if magnitude.and_then(|k|n.bits().max(d.bits()).checked_mul(k)).is_none_or(|bits|bits>1_000_000){
            return Err("rational content power exceeds the supported one-million-bit bound".into());
        }
    }
    Ok(())
}

pub fn as_content_primitive(e:&Expr,radical:bool,clear:bool)->Result<(Number,Expr),String>{
    fn visit(e:&Expr,radical:bool,clear:bool,cache:&mut FxHashMap<usize,(Expr,(Number,Expr))>)->Result<(Number,Expr),String>{
        if let Some((_,result))=cache.get(&e.dag_id()){return Ok(result.clone());}
        let result=match e.kind(){
            Kind::Num(n) if n.is_rational() && !n.is_zero()=>(n.abs(),Expr::int(if n.is_negative(){-1}else{1})),
            Kind::Add(..)=>{
                let mut processed=Vec::new();
                for term in args(e){let(c,p)=visit(&term,radical,clear,cache)?;processed.push(keep(&c,&p));}
                let sum=ops::add(&processed);
                let terms=if sum.is_add(){args(&sum)}else{vec![sum.clone()]};
                let mut content=Number::ZERO;
                let mut coefficients=Vec::new();
                for term in &terms {
                    let (c,p)=ops::as_coeff_mul(term);
                    let (c,p)=if c.is_rational(){(c,p)}else{(Number::ONE,term.clone())};
                    content=content.gcd(&c.abs());coefficients.push((c,p));
                }
                if content.is_zero(){(Number::ONE,sum)}else{
                    let mut primitive=ops::add(&coefficients.iter().map(|(c,p)|keep(&c.div(&content).unwrap(),p)).collect::<Vec<_>>());
                    if !clear && primitive.is_add(){
                        let(n,d)=content.numer_denom().unwrap();
                        if d!=BigInt::one(){
                            let divisor=Number::from_bigint(d);
                            if coefficients.iter().any(|(c,_)|c.div(&content).unwrap().div(&divisor).unwrap().is_integer()){
                                primitive=ops::mul(&[Expr::num(Number::ONE.div(&divisor).unwrap()),primitive]);
                                content=Number::from_bigint(n);
                            }
                        }
                    }
                    if radical && primitive.is_add(){primitive=radical_primitive(&primitive)?;}
                    (content,primitive)
                }
            },
            Kind::Mul(..)=>{
                let mut content=Number::ONE;let mut factors=Vec::new();
                for factor in args(e){let(c,p)=visit(&factor,radical,clear,cache)?;content=content.mul(&c);if !p.is_one(){factors.push(p);}}
                (content,ops::mul(&factors))
            },
            Kind::Pow(base,exponent)=>{
                let (base,exponent)=if let Some((n,d))=base.as_number().and_then(Number::numer_denom){
                    if n.is_one() && d!=BigInt::one(){(Expr::num(Number::from_bigint(d)),ops::neg(exponent))}else{(base.clone(),exponent.clone())}
                }else{(base.clone(),exponent.clone())};
                let(base,exponent)=(&base,&exponent);
                let(c,p)=visit(base,radical,clear,cache)?;
                let(ec,ep)=visit(exponent,radical,clear,cache)?;
                let exponent=keep(&ec,&ep);
                if let Some(n)=exponent.as_i64().filter(|_|!base.is_rational_number()){
                    check_power_size(&c,exponent.as_number().unwrap())?;
                    match c.pow_i64(n){Some(content)=>(content,ops::pow(&p,&exponent)),None=>(Number::ONE,ops::pow(&keep(&c,&p),&exponent))}
                }else if base.is_rational_number() && !base.is_zero(){
                    let(h,tail)=ops::as_coeff_add(&ep);
                    let head=h.mul(&ec);
                    if head.is_rational(){
                        check_power_size(base.as_number().unwrap(),&head)?;
                        let power=ops::pow(base,&Expr::num(head.clone()));
                        let (content,remainder)=if let Some(n)=power.as_number().filter(|n|n.is_rational()){(n.clone(),keep(&ec,&tail))}else{
                            let(n,d)=head.numer_denom().unwrap();let integer=n.div_floor(&d);
                            let extracted=ops::pow(base,&Expr::num(Number::from_bigint(integer.clone())));
                            if let Some(value)=extracted.as_number().filter(|n|n.is_rational()){
                                (value.clone(),keep(&ec,&ops::add(&[tail,Expr::num(head.sub(&Number::from_bigint(integer)).div(&ec).unwrap())])))
                            }else{return Err("rational content exponent is too large".into());}
                        };
                        (content,ops::pow(base,&remainder))
                    }else{(Number::ONE,ops::pow(&keep(&c,&p),&exponent))}
                }else if exponent.is_rational_number(){
                    let raised=ops::pow(&Expr::num(c.clone()),&exponent);
                    let (coefficient,rest)=ops::as_coeff_mul(&raised);
                    if coefficient.is_rational(){(coefficient,ops::mul(&[rest,ops::pow(&p,&exponent)]))}
                    else{(Number::ONE,ops::pow(&keep(&c,&p),&exponent))}
                }else{(Number::ONE,ops::pow(&keep(&c,&p),&exponent))}
            },
            _=>(Number::ONE,e.clone()),
        };
        cache.insert(e.dag_id(),(e.clone(),result.clone()));Ok(result)
    }
    visit(e,radical,clear,&mut FxHashMap::default())
}

fn radical_primitive(e:&Expr)->Result<Expr,String>{
    let terms=args(e);let mut decomposed=Vec::new();let mut common:Option<FxHashMap<u32,BigInt>>=None;
    for term in &terms {
        let factors=if term.is_mul(){args(term)}else{vec![term.clone()]};
        let mut groups:FxHashMap<u32,BigInt>=FxHashMap::default();let mut rest=Vec::new();
        for factor in factors{
            let mut taken=false;
            if let Kind::Pow(base,exponent)=factor.kind(){
                if let (Some(base),Some(exp))=(base.as_number().and_then(Number::to_bigint),exponent.as_number()){
                    if base>BigInt::from(0){if let Some((p,q))=exp.numer_denom(){if q>BigInt::one(){
                        let p=p.to_u32().filter(|p|*p<=4096).ok_or_else(||"radical content exponents exceed the supported bound".to_string())?;
                        let q=q.to_u32().ok_or_else(||"radical content denominator exceeds the supported bound".to_string())?;
                        if base.bits().checked_mul(p as u64).is_none_or(|bits|bits>1_000_000){return Err("radical content power exceeds the supported one-million-bit bound".into());}
                        let entry=groups.entry(q).or_insert_with(BigInt::one);*entry *= base.pow(p);taken=true;
                    }}}
                }
            }
            if !taken{rest.push(factor);}
        }
        if groups.is_empty(){return Ok(e.clone());}
        match &mut common{None=>common=Some(groups.clone()),Some(common)=>common.retain(|q,g|if let Some(value)=groups.get(q){*g=g.gcd(value);true}else{false})}
        decomposed.push((groups,rest));
    }
    let mut common:Vec<_>=common.unwrap_or_default().into_iter().filter(|(_,g)|!g.is_one()).collect();common.sort_by_key(|(q,_)|*q);
    if common.is_empty(){return Ok(e.clone());}
    let factor=ops::mul(&common.iter().map(|(q,g)|ops::pow(&Expr::num(Number::from_bigint(g.clone())),&Expr::rational(1,*q as i64))).collect::<Vec<_>>());
    let mut primitive=Vec::new();
    for (mut groups,mut rest) in decomposed{
        for (q,g) in &common{
            let left=groups.remove(q).unwrap()/g;
            if !left.is_one(){rest.push(ops::pow(&Expr::num(Number::from_bigint(left)),&Expr::rational(1,*q as i64)));}
        }
        for(q,g)in groups{rest.push(ops::pow(&Expr::num(Number::from_bigint(g)),&Expr::rational(1,q as i64)));}
        primitive.push(ops::mul(&rest));
    }
    Ok(ops::mul(&[factor,ops::add(&primitive)]))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(source:&str)->Expr{crate::parser::parse_str(source).unwrap()}
    #[test]
    fn public_factor_views_preserve_exact_and_held_powers(){
        assert_eq!(as_base_exp(&p("I")),(Expr::minus_one(),Expr::rational(1,2)));
        assert_eq!(as_base_exp(&p("Pow(1/2,x,evaluate=False)")),(Expr::int(2),p("-x")));
        assert_eq!(as_base_exp(&p("Pow(2/3,x,evaluate=False)")),(Expr::rational(2,3),p("x")));
        assert_eq!(as_base_exp(&p("exp(x)")),(Expr::e(),p("x")));
        assert_eq!(as_powers_dict(&p("Mul(x,1/x,evaluate=False)")),vec![(p("x"),Expr::zero())]);
        assert_eq!(as_powers_dict(&Expr::neg_oo()),vec![(Expr::minus_one(),Expr::one()),(Expr::oo(),Expr::one())]);
        let mut dag=p("x");
        for _ in 0..80{dag=Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,vec![(dag.clone(),Expr::one()),(dag,Expr::one())]));}
        assert_eq!(as_powers_dict(&dag),vec![(p("x"),Expr::num(Number::from_bigint(BigInt::one()<<80)))]);
        assert_eq!(commutativity(&dag),Some(true));
    }
    #[test]
    fn public_factor_partition_and_syntactic_signs(){
        assert_eq!(args_cnc(&p("-2*x*y"),true),(vec![p("-1"),p("2"),p("x"),p("y")],vec![]));
        assert_eq!(args_cnc(&p("-2*x*y"),false),(vec![p("-2"),p("x"),p("y")],vec![]));
        assert_eq!(args_cnc(&p("Mul(x,UnevaluatedExpr(y),z,evaluate=False)"),true),(vec![p("x")],vec![p("UnevaluatedExpr(y)"),p("z")]));
        for(source,expected)in[("x-y",false),("y-x",true),("1-x",true),("x-1",false),("2-x-y",true),("-2+x+y",false),("3-sqrt(2)",false),("sqrt(2)-3",true),("-x-y+z",true),("x*(y-x)",false),("-oo*x",true),("zoo*x",false),("Add(x,-x,evaluate=False)",false)]{
            assert_eq!(could_extract_minus_sign(&p(source)),expected,"{source}");
        }
        let n=Expr::symbol_with("n",crate::Assumptions::NEGATIVE);
        assert!(!could_extract_minus_sign(&n));
    }
    #[test]
    fn coefficients_dependencies_and_bound_independence(){
        let x=p("x");
        assert_eq!(as_coeff_add(&p("3*x+y+2"),&[]),(Expr::int(2),vec![p("y"),p("3*x")]));
        assert_eq!(as_coeff_add(&p("3*x+y+2"),std::slice::from_ref(&x)),(p("y+2"),vec![p("3*x")]));
        assert_eq!(as_coeff_mul(&p("6*x*y"),std::slice::from_ref(&x),true),(p("6*y"),vec![x.clone()]));
        assert_eq!(as_coeff_add(&p("6*x*y"),&[p("x*y")]),(Expr::zero(),vec![p("6*x*y")]));
        assert_eq!(as_coeff_mul(&p("-2.0*x"),&[],true),(Expr::minus_one(),vec![p("2.0"),x.clone()]));
        assert_eq!(as_independent(&p("x*y+z"),std::slice::from_ref(&x),None),(p("z"),p("x*y")));
        assert_eq!(as_independent(&p("x*y+z"),std::slice::from_ref(&x),Some(false)),(Expr::one(),p("x*y+z")));
        for dep in ["x+y","x*y","x+1","sin(x)"]{
            assert_eq!(as_independent(&p("Add(x,2,evaluate=False)"),&[p(dep)],None),(p("x+2"),Expr::zero()));
        }
        let held=p("Subs(f(x),x,0)");
        assert_eq!(as_independent(&held,std::slice::from_ref(&x),None),(held.clone(),Expr::one()));
        assert_eq!(as_coeff_mul(&held,&[x],true),(Expr::one(),vec![held]));
    }
    #[test]
    fn recursive_content_clear_and_radicals(){
        let cases=[
            ("2+2*x+2*y*(3+3*y)",false,true,"2","1+x+3*y*(1+y)"),
            ("(2+6*x)**2",false,true,"4","(3*x+1)**2"),
            ("x/2+y",false,false,"1","x/2+y"),
            ("x/2+y",false,true,"1/2","x+2*y"),
            ("2*x/3+4*y/9",false,false,"2/9","3*x+2*y"),
            ("2*sqrt(2)+4*sqrt(10)",true,true,"2","sqrt(2)*(1+2*sqrt(5))"),
            ("sqrt(4+4*sqrt(2))",false,true,"2","sqrt(1+sqrt(2))"),
            ("sqrt(3+3*sqrt(2))",false,true,"1","sqrt(3)*sqrt(1+sqrt(2))"),
            ("3**(2+2*x)",false,true,"9","3**(2*x)"),
            ("(-2)**(x+3)",false,true,"-8","(-2)**x"),
            ("(1/2)**(x+2)",false,true,"1/4","2**(-x)"),
        ];
        for(input,radical,clear,c,primitive)in cases{
            let got=as_content_primitive(&p(input),radical,clear).unwrap();
            assert_eq!(got,(p(c).as_number().unwrap().clone(),p(primitive)),"{input}");
        }
    }
}
