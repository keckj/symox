//! Expression collection, with native coefficient partitions and fallible callbacks.
//!
//! Collection never expands sums. Composite patterns are replaced by fresh
//! assumption-preserving dummies, collected, then restored structurally.
use rustc_hash::FxHashMap;
use crate::expr::{Expr,Func,Kind};
use crate::ops::{add,as_coeff_mul,div,mul,neg,pow};
use crate::traverse::{map_children,ordered_args};

#[derive(Clone,Copy,Debug)]
pub struct Options { pub exact:bool, pub discover:bool, pub distribute_order_term:bool }
impl Default for Options {fn default()->Self{Self{exact:false,discover:false,distribute_order_term:true}}}

#[derive(Clone,Debug)]
pub enum Collection { Expression(Expr), Dictionary(Vec<(Expr,Expr)>) }
#[derive(Debug)]
pub enum Error<E>{Unsupported(String),Callback(E),Pattern(crate::pattern::MatchError)}
impl<E> From<String> for Error<E>{fn from(message:String)->Self{Self::Unsupported(message)}}
impl<E> From<crate::pattern::MatchError> for Error<E>{fn from(error:crate::pattern::MatchError)->Self{Self::Pattern(error)}}

type Callback<'a,E>=Option<&'a mut dyn FnMut(&Expr)->Result<Expr,E>>;

struct Collector {
    patterns:Vec<Expr>, restore:FxHashMap<Expr,Expr>, options:Options,
    cache:FxHashMap<usize,(Expr,Expr)>, power_cache:crate::expand::ExpandBudget,
}

/// `evaluate=False` returns the native coefficient partition. With evaluation,
/// recurse into Add/Mul/Pow exactly before applying the coefficient callback.
/// Callback failures are returned unchanged to the caller.
pub fn run<E>(expression:&Expr,patterns:&[Expr],options:Options,evaluate:bool,mut callback:Callback<'_,E>)->Result<Collection,Error<E>>{
    let(mut collector,expression)=Collector::prepare(expression,patterns,options)?;
    let output=if evaluate {
        Collection::Expression(collector.evaluate(&expression,&mut callback)?)
    }else{
        Collection::Dictionary(collector.partition_callback(&expression,&mut callback)?)
    };
    let restore=|e:&Expr|crate::subs::xreplace(e,&collector.restore);
    Ok(match output{
        Collection::Expression(e)=>Collection::Expression(restore(&e)),
        Collection::Dictionary(pairs)=>Collection::Dictionary(pairs.into_iter().map(|(key,value)|(restore(&key),restore(&value))).collect()),
    })
}

fn derived_assumptions(e:&Expr)->crate::Assumptions {
    use crate::assumptions as a;
    let properties=a::exact_queries(||[("real",a::is_real(e)),("positive",a::is_positive(e)),("negative",a::is_negative(e)),
        ("nonnegative",a::is_nonnegative(e)),("nonpositive",a::is_nonpositive(e)),("zero",a::is_zero(e)),
        ("nonzero",a::is_nonzero(e)),("integer",a::is_integer(e)),("rational",a::is_rational(e)),
        ("even",a::is_even(e)),("odd",a::is_odd(e)),("finite",a::is_finite(e))]);
    // These are structural facts, never approximate numerical certificates.
    crate::Assumptions::from_kwargs(properties.into_iter().filter_map(|(key,value)|value.map(|v|(key,v))))
        .unwrap_or(crate::Assumptions::NONE)
}

fn check_pattern(e:&Expr,seen:&mut FxHashMap<usize,Expr>)->Result<(),String>{
    if seen.insert(e.dag_id(),e.clone()).is_some(){return Ok(());}
    if matches!(e.kind(),Kind::Func(Func::MatrixProxy(_),_)){return Err("collect requires scalar patterns".into());}
    for arg in e.args(){check_pattern(&arg,seen)?;}Ok(())
}

fn node_count(e:&Expr,memo:&mut FxHashMap<usize,(Expr,u128)>)->u128{
    if let Some((_,n))=memo.get(&e.dag_id()){return *n;}
    let mut n=1u128;for arg in e.args(){n=n.saturating_add(node_count(&arg,memo));}
    memo.insert(e.dag_id(),(e.clone(),n));n
}

fn derivative_parts(e:&Expr)->Result<Option<(Expr,Expr,u32)>,String>{
    let Kind::Derivative(body,variables)=e.kind()else{return Ok(None)};
    let Some((variable,_))=variables.first()else{return Err("collect requires positive derivative orders".into())};
    if variables.iter().any(|(v,_)|v!=variable){return Err("collect of multivariate derivatives is not implemented".into());}
    let mut order=0u32;for(_,n)in variables{order=order.checked_add(*n).ok_or_else(||"collect derivative order overflow".to_string())?;}
    Ok(Some((body.clone(),variable.clone(),order)))
}

/// Extend exact algebraic substitution to exponential integer powers and
/// derivative towers, the two special matches used by SymPy collect.
fn replace_pattern(e:&Expr,old:&Expr,new:&Expr)->Result<Expr,String>{
    let derivative=derivative_parts(old)?;
    if !matches!(old.kind(),Kind::Func(Func::Exp,_)|Kind::Pow(..)) && derivative.is_none() && !old.is_minus_one(){return Ok(crate::subs::subs1(e,old,new));}
    fn walk(e:&Expr,old:&Expr,new:&Expr,derivative:&Option<(Expr,Expr,u32)>,cache:&mut FxHashMap<usize,(Expr,Expr)>)->Result<Expr,String>{
        if e==old{return Ok(new.clone());}
        // SymPy replaces composite patterns with expr.subs, which leaves
        // ComplexRootOf untouched (its _eval_subs returns self).
        if let Kind::Func(Func::Undefined(head),args)=e.kind(){if crate::poly::rootof::is_rootof_head(head,args.len()){return Ok(e.clone());}}
        if let Some((_,value))=cache.get(&e.dag_id()){return Ok(value.clone());}
        let symbol_power=if let Kind::Pow(base,exponent)=old.kind(){
            if e.is_symbol() && e==base{let ratio=div(&Expr::one(),exponent);if ratio.is_integer_number(){Some(pow(new,&ratio))}else{None}}else{None}
        }else{None};
        let special=if old.is_minus_one() && e.is_one(){Some(neg(new))}else if symbol_power.is_some(){symbol_power}else if let(Kind::Func(Func::Exp,args),Kind::Func(Func::Exp,pattern))=(e.kind(),old.kind()){
            let ratio=div(&args[0],&pattern[0]);if ratio.is_integer_number(){Some(pow(new,&ratio))}else{None}
        }else if let Some((body,var,n))=derivative{
            match derivative_parts(e)?{Some((ebody,evar,en))if ebody==*body && evar==*var && en>=*n=>{
                Some(crate::diff::derivative_unevaluated(new,vec![(var.clone(),en-*n)]))
            },_=>None}
        }else{None};
        let result=if let Some(result)=special{result}else{
            let mut error=None;
            let mut child=|value:&Expr|match walk(value,old,new,derivative,cache){Ok(value)=>value,Err(e)=>{error=Some(e);value.clone()}};
            // Public children include whole factor powers and numerical
            // coefficients, unlike map_children's compact Add/Mul payloads.
            let result=if e.is_add() || e.is_mul(){
                let args=ordered_args(e);let mapped=args.iter().map(&mut child).collect::<Vec<_>>();
                if args.iter().zip(&mapped).all(|(a,b)|a.same_presentation(b)){e.clone()}else if e.is_add(){add(&mapped)}else{mul(&mapped)}
            }else{map_children(e,&mut child)};
            if let Some(error)=error{return Err(error);}result
        };
        cache.insert(e.dag_id(),(e.clone(),result.clone()));Ok(result)
    }
    // Replace whole powers before walking symbol bases: sqrt(x) must become
    // the dummy itself before residual x occurrences become dummy**2.
    let prepared=crate::subs::subs1(e,old,new);
    walk(&prepared,old,new,&derivative,&mut FxHashMap::default())
}

impl Collector {
    fn prepare(expression:&Expr,patterns:&[Expr],options:Options)->Result<(Self,Expr),String>{
        let mut restore=FxHashMap::default();let mut replacements:Vec<(Expr,Expr)>=Vec::new();let mut normalized=Vec::new();
        for pattern in patterns{
            check_pattern(pattern,&mut FxHashMap::default())?;
            if pattern.is_symbol()||neg(pattern).is_symbol()||crate::pattern::has_wild(pattern){normalized.push(pattern.clone());}
            else if let Some((_,dummy))=replacements.iter().find(|(old,_)|old==pattern){normalized.push(dummy.clone());}
            else{
                let dummy=Expr::dummy("collect",derived_assumptions(pattern));
                restore.insert(dummy.clone(),pattern.clone());replacements.push((pattern.clone(),dummy.clone()));normalized.push(dummy);
            }
        }
        let mut weights=FxHashMap::default();
        replacements.sort_by(|(a,_),(b,_)|node_count(b,&mut weights).cmp(&node_count(a,&mut weights)).then_with(||{
            if node_count(a,&mut weights)<10_000{crate::sortkey::sort_key(a).cmp(&crate::sortkey::sort_key(b))}else{a.cmp_struct(b)}
        }));
        let mut expression=expression.clone();
        for(old,new)in replacements{expression=replace_pattern(&expression,&old,&new)?;}
        let mut options=options;
        if options.discover {
            options.discover=false;
            // Basic.has_free on a normalized symbol must respect bound
            // variables. Composite non-Wild patterns have already become
            // fresh symbols, with their assumptions and restoration retained.
            let depends=|e:&Expr|{
                let free=crate::traverse::free_symbols(e);
                normalized.iter().any(|p|if p.is_symbol(){free.contains(p)}else{
                    let positive=neg(p);
                    if positive.is_symbol(){free.contains(&positive)&&crate::traverse::has(e,p)}
                    else{crate::traverse::has(e,p)}
                })
            };
            let terms=if expression.is_add(){ordered_args(&expression)}else{vec![expression.clone()]};
            let mut extra=FxHashMap::default();
            for term in terms {
                if normalized.contains(&term)||!depends(&term){continue;}
                let generator=if term.is_mul(){mul(&ordered_args(&term).into_iter().filter(&depends).collect::<Vec<_>>())}else{term};
                if !normalized.contains(&generator){extra.insert(generator,());}
            }
            let mut extra=extra.into_keys().collect::<Vec<_>>();
            extra.sort_by_cached_key(crate::sortkey::sort_key);
            let simple=extra.iter().all(|p|matches!(p.kind(),Kind::Pow(base,_) if normalized.contains(base)));
            normalized.extend(extra);
            if !simple {
                let (mut nested,expression)=Self::prepare(&expression,&normalized,options)?;
                for value in nested.restore.values_mut(){*value=crate::subs::xreplace(value,&restore);}
                nested.restore.extend(restore);
                return Ok((nested,expression));
            }
            // In the simple-power case exact=None compares each discovered
            // exponent exactly, without the integer-power dummy replacement.
            options.exact=true;
        }
        Ok((Self{patterns:normalized,restore,options,cache:FxHashMap::default(),power_cache:crate::expand::ExpandBudget::default()},expression))
    }

    fn power_base(&mut self,e:&Expr)->Result<Expr,String>{
        if !e.is_pow(){return Ok(e.clone());}
        let options=crate::expand::ExpandOptions{deep:false,basic:false,mul:false,multinomial:false,power_base:true,power_exp:false,log:false,trig:false,func:false,force:false,kroneckerproduct:false};
        crate::expand::expand_with_options(e,&options,&mut self.power_cache).map_err(|e|e.to_string())
    }

    fn evaluate<E>(&mut self,e:&Expr,callback:&mut Callback<'_,E>)->Result<Expr,Error<E>>{
        if callback.is_none(){if let Some((_,value))=self.cache.get(&e.dag_id()){return Ok(value.clone());}}
        let result=match e.kind(){
            Kind::Mul(..)=>{
                let factors=ordered_args(e).iter().map(|arg|self.evaluate(arg,callback)).collect::<Result<Vec<_>,_>>()?;mul(&factors)
            },
            Kind::Pow(base,exponent)=>{let base=self.evaluate(base,callback)?;pow(&base,exponent)},
            _=>{
                let expression=if e.is_add(){
                    let order=crate::series::get_o(e);
                    let mut terms=Vec::new();
                    for term in ordered_args(e){if order.as_ref()!=Some(&term){terms.push(self.evaluate(&term,callback)?);}}
                    if let Some(order)=order{terms.push(order);}add(&terms)
                }else{e.clone()};
                let pairs=self.partition_callback(&expression,callback)?;
                add(&pairs.into_iter().map(|(key,value)|mul(&[key,value])).collect::<Vec<_>>())
            }
        };
        if callback.is_none(){self.cache.insert(e.dag_id(),(e.clone(),result.clone()));}Ok(result)
    }

    fn partition_callback<E>(&mut self,e:&Expr,callback:&mut Callback<'_,E>)->Result<Vec<(Expr,Expr)>,Error<E>>{
        let mut pairs=self.partition(e)?;
        if let Some(callback)=callback{for(_,value)in &mut pairs{*value=callback(value).map_err(Error::Callback)?;}}
        Ok(pairs)
    }

    fn partition<E>(&mut self,e:&Expr)->Result<Vec<(Expr,Expr)>,Error<E>>{
        let order=if self.options.distribute_order_term{crate::series::get_o(e).filter(|o|!self.patterns.iter().any(|p|crate::traverse::has(o,p)))}else{None};
        let expression=if order.is_some(){crate::series::remove_o(e)}else{e.clone()};
        let terms=if expression.is_add(){ordered_args(&expression)}else{vec![expression]};
        let mut groups:Vec<(Expr,Vec<Expr>)>=Vec::new();let mut indices:FxHashMap<Expr,usize>=FxHashMap::default();let mut disliked=Vec::new();
        for term in terms{
            let term=self.power_base(&term)?;
            let factors=if term.is_mul(){ordered_args(&term)}else{vec![term.clone()]};
            let parsed=factors.iter().map(parse_factor).collect::<Result<Vec<_>,_>>()?;
            let mut matched=None;
            for pattern in &self.patterns{
                if crate::pattern::has_wild(pattern)||!pattern.is_symbol()&&!neg(pattern).is_symbol(){
                    if let Some((key,coefficient))=match_factors(&factors,&parsed,pattern,self.options.exact)?{matched=Some((key,coefficient));break;}
                    continue;
                }
                let symbol=if pattern.is_symbol(){pattern.clone()}else{neg(pattern)};
                for(index,factor)in parsed.iter().enumerate(){
                    if factor.base==symbol && factor.symbolic.is_none() && (!self.options.exact || factor.exponent.is_one() && factor.derivative.is_none()){
                        let key=factors[index].clone();
                        let coefficient=mul(&factors.iter().enumerate().filter(|(i,_)|*i!=index).map(|(_,e)|e.clone()).collect::<Vec<_>>());
                        matched=Some((key,coefficient));break;
                    }
                }
                if matched.is_some(){break;}
            }
            if let Some((key,coefficient))=matched{
                let (key,coefficient)=(rebuild(&key),rebuild(&coefficient));
                let index=if let Some(index)=indices.get(&key){*index}else{let index=groups.len();indices.insert(key.clone(),index);groups.push((key,Vec::new()));index};
                groups[index].1.push(coefficient);
            }else{disliked.push(term);}
        }
        let mut result=groups.into_iter().map(|(key,values)|(key,add(&values))).collect::<Vec<_>>();
        let disliked=add(&disliked);if !disliked.is_zero(){result.push((Expr::one(),disliked));}
        if let Some(order)=order{for(_,value)in &mut result{*value=add(&[value.clone(),order.clone()]);}}
        Ok(result)
    }
}

/// SymPy's collect rebuilds a matched term from its parsed factors, so a
/// factor exp(c*t) with a rational c != 1 comes back as Pow(exp(t), c):
/// exp(x/3) becomes exp(x)**(1/3) unless the power evaluates.
fn rebuild(e:&Expr)->Expr{
    let factor=|f:&Expr|->Expr{
        let Kind::Func(Func::Exp,args)=f.kind() else{return f.clone()};
        let Kind::Mul(coefficient,_)=args[0].kind() else{return f.clone()};
        if !coefficient.is_rational()||coefficient.is_one(){return f.clone();}
        let tail=div(&args[0],&Expr::num(coefficient.clone()));
        pow(&crate::functions::func(Func::Exp,vec![tail]),&Expr::num(coefficient.clone()))
    };
    if e.is_mul(){mul(&ordered_args(e).iter().map(factor).collect::<Vec<_>>())}else{factor(e)}
}

struct ParsedFactor {base:Expr,exponent:Expr,symbolic:Option<Expr>,derivative:Option<(Expr,u32)>}
fn parse_factor(e:&Expr)->Result<ParsedFactor,String>{
    let(base,exponent)=match e.kind(){Kind::Pow(base,exponent)=>(base.clone(),exponent.clone()),_=>(e.clone(),Expr::one())};
    let derivative=derivative_parts(&base)?;
    let base=derivative.as_ref().map(|(body,_,_)|body.clone()).unwrap_or(base);
    let (exponent,symbolic)=if exponent.as_number().is_none(){let(coefficient,tail)=as_coeff_mul(&exponent);(Expr::num(coefficient),Some(tail))}else{(exponent,None)};
    Ok(ParsedFactor{base,exponent,symbolic,derivative:derivative.map(|(_,variable,order)|(variable,order))})
}

/// Match product factors without substituting Wilds into the expression.
/// Each factor comparison has its own bindings, as in reference collect;
/// the actual matched factors, not the pattern, become the dictionary key.
fn match_factors<E>(factors:&[Expr],parsed:&[ParsedFactor],pattern:&Expr,exact:bool)->Result<Option<(Expr,Expr)>,Error<E>>{
    let patterns=if pattern.is_mul(){ordered_args(pattern)}else{vec![pattern.clone()]};
    if parsed.len()<patterns.len(){return Ok(None);}
    let mut taken=vec![false;parsed.len()];
    for pattern in patterns{
        let pattern=parse_factor(&pattern)?;
        if pattern.base.is_number()&&pattern.exponent.is_one()&&pattern.symbolic.is_none(){continue;}
        let mut found=None;
        for(index,term)in parsed.iter().enumerate(){
            if taken[index]||crate::pattern::pattern_match(&term.base,&pattern.base)?.is_none(){continue;}
            let symbolic=match (&term.symbolic,&pattern.symbolic){
                (None,None)=>true,(Some(a),Some(b))=>crate::pattern::pattern_match(a,b)?.is_some(),_=>false,
            };
            if symbolic&&(!exact||term.exponent==pattern.exponent&&term.derivative==pattern.derivative){found=Some(index);break;}
        }
        let Some(index)=found else{return Ok(None)};
        taken[index]=true;
    }
    Ok(Some((mul(&factors.iter().zip(&taken).filter(|(_,taken)|**taken).map(|(factor,_)|factor.clone()).collect::<Vec<_>>()),
             mul(&factors.iter().zip(&taken).filter(|(_,taken)|!**taken).map(|(factor,_)|factor.clone()).collect::<Vec<_>>()))))
}

#[cfg(test)]
mod tests{
    use super::*;
    use std::convert::Infallible;
    fn p(text:&str)->Expr{crate::parser::parse_str(text).unwrap()}
    fn expression(text:&str,patterns:&[&str],exact:bool)->Expr{
        let patterns=patterns.iter().map(|s|p(s)).collect::<Vec<_>>();
        match run::<Infallible>(&p(text),&patterns,Options{exact,..Options::default()},true,None).unwrap(){Collection::Expression(e)=>e,_=>unreachable!()}
    }
    fn dictionary(text:&str,patterns:&[&str],exact:bool)->Vec<(Expr,Expr)>{
        let patterns=patterns.iter().map(|s|p(s)).collect::<Vec<_>>();
        match run::<Infallible>(&p(text),&patterns,Options{exact,..Options::default()},false,None).unwrap(){Collection::Dictionary(e)=>e,_=>unreachable!()}
    }
    #[test]
    fn pattern_assumptions_never_promote_numerical_refinement(){
        let difference=p("pi-E");
        let facts=derived_assumptions(&difference);
        let signs=crate::Assumptions::POSITIVE.0|crate::Assumptions::NEGATIVE.0|crate::Assumptions::ZERO.0;
        assert_eq!(facts.0&signs,0);
        assert!(facts.has(crate::Assumptions::REAL));
        // Restoring the query context must keep structural proofs enabled.
        assert!(derived_assumptions(&p("exp(2)")).has(crate::Assumptions::POSITIVE));
    }
    #[test]
    fn discover_compound_generators_and_exact_simple_powers(){
        for(source,patterns,expected)in [
            ("x*exp(x)+y*sin(x)+2*sin(x)+3*x",vec!["x"],"x*exp(x)+3*x+(y+2)*sin(x)"),
            ("a*x*y+x*y+b*x+x",vec!["x","y"],"x*y*(a+1)+x*(b+1)"),
            ("a*x**2+b*x**2+c*x**3+d*x**3",vec!["x"],"x**2*(a+b)+x**3*(c+d)"),
            ("a*sqrt(x)+b*sqrt(x)+c*x+d*x",vec!["x"],"sqrt(x)*(a+b)+x*(c+d)"),
        ]{
            let patterns=patterns.iter().map(|text|p(text)).collect::<Vec<_>>();
            let result=run::<Infallible>(&p(source),&patterns,Options{discover:true,..Options::default()},true,None).unwrap();
            let Collection::Expression(result)=result else{unreachable!()};
            assert_eq!(result,p(expected),"{source}");
        }
    }
    #[test]
    fn wild_factor_patterns_keep_concrete_keys_and_exclusions(){
        let wild=Expr::wild("w",Vec::new());
        let pattern=pow(&wild,&p("y"));
        let source=p("a*x**y-b*x**y");
        let result=run::<Infallible>(&source,&[pattern],Options::default(),false,None).unwrap();
        let Collection::Dictionary(result)=result else{unreachable!()};
        assert_eq!(result,vec![(p("x**y"),p("a-b"))]);
        let wild=Expr::wild("w",vec![p("z")]);
        let pattern=pow(&wild,&p("y"));
        let result=run::<Infallible>(&p("a*x**y+b*x**y+c*z**y"),&[pattern],Options::default(),false,None).unwrap();
        let Collection::Dictionary(result)=result else{unreachable!()};
        assert!(result.contains(&(p("x**y"),p("a+b"))));
        assert!(result.contains(&(Expr::one(),p("c*z**y"))));
    }
    #[test]
    fn partitions_reconstruct_and_preserve_exact_flag(){
        let source="a*x**2+b*x**2+a*x-b*x+c";
        assert_eq!(expression(source,&["x"],false),p("x**2*(a+b)+x*(a-b)+c"));
        assert_eq!(expression(source,&["x"],true),p("a*x**2+b*x**2+x*(a-b)+c"));
        let pairs=dictionary(source,&["x"],false);
        assert_eq!(pairs.len(),3);
        assert!(pairs.contains(&(p("x**2"),p("a+b"))));
        assert!(pairs.contains(&(p("x"),p("a-b"))));
        assert!(pairs.contains(&(Expr::one(),p("c"))));
        for source in ["x**2*y+x*y+x+y","a*sin(x)+b*sin(x)+c","a*x**(2*c)+b*x**(2*c)","a*x*y+b*x*y+a*x+b*x"]{
            for patterns in [vec!["x"],vec!["x","y"],vec!["sin(x)"],vec!["x*y"],vec!["x**c"]]{
                for exact in [false,true]{
                    let pairs=dictionary(source,&patterns,exact);
                    let reconstructed=add(&pairs.iter().map(|(k,v)|mul(&[k.clone(),v.clone()])).collect::<Vec<_>>());
                    assert!(crate::expand::expand(&crate::ops::sub(&reconstructed,&p(source))).is_zero(),"{source} by {patterns:?}");
                }
            }
        }
        assert!(dictionary("0",&["x"],false).is_empty());
    }
    #[test]
    fn composite_patterns_exponentials_and_derivative_towers(){
        for(source,pattern,key)in [
            ("a*exp(2*x)+b*exp(2*x)","exp(x)","exp(2*x)"),
            ("a*exp(-x)+b*exp(-x)","exp(x)","exp(-x)"),
            ("a*x+b*x","sqrt(x)","x"),
            ("a*x+b*x","1/x","x"),
            ("a*x**(2*c)+b*x**(2*c)","x**c","x**(2*c)"),
            ("a*x**2*y**2+b*x**2*y**2","x*y","x**2*y**2"),
            ("a*Derivative(f(x),(x,2))+b*Derivative(f(x),(x,2))","f(x)","Derivative(f(x),(x,2))"),
            ("a*Derivative(f(x),(x,2))+b*Derivative(f(x),(x,2))","Derivative(f(x),x)","Derivative(f(x),(x,2))"),
        ]{
            let pairs=dictionary(source,&[pattern],false);assert_eq!(pairs,vec![(p(key),p("a+b"))],"{source} by {pattern}");
        }
        assert_eq!(dictionary("a*x+b*x+c*sqrt(x)",&["sqrt(x)"],false),vec![(p("x"),p("a+b")),(p("sqrt(x)"),p("c"))]);
        assert_eq!(dictionary("x+1",&["-1"],false),vec![(Expr::minus_one(),Expr::minus_one()),(Expr::one(),p("x"))]);
        assert_eq!(dictionary("a*exp(x/2)+b*exp(x/2)",&["exp(x)"],false),vec![(Expr::one(),p("a*exp(x/2)+b*exp(x/2)"))]);
        assert_eq!(expression("a*Derivative(f(x),(x,2))+b*Derivative(f(x),(x,2))",&["Derivative(f(x),x)"],true),p("a*Derivative(f(x),(x,2))+b*Derivative(f(x),(x,2))"));
    }
    #[test]
    fn callback_recursion_order_is_preserved_and_errors_survive(){
        let e=p("a*x+b*x+x**2");let x=p("x");let mut calls=Vec::new();
        let mut callback=|e:&Expr|->Result<Expr,&'static str>{calls.push(e.clone());Ok(pow(e,&Expr::int(2)))};
        let result=run(&e,&[x.clone()],Options::default(),true,Some(&mut callback)).unwrap();
        let Collection::Expression(result)=result else{unreachable!()};
        assert_eq!(result,p("x**2+x*(a**2+b**2)**2"));
        assert_eq!(calls,vec![Expr::one(),p("a"),Expr::one(),p("b"),Expr::one(),Expr::one(),p("a**2+b**2")]);
        let mut callback=|_:&Expr|->Result<Expr,&'static str>{Err("sentinel")};
        assert!(matches!(run(&e,&[x],Options::default(),true,Some(&mut callback)),Err(Error::Callback("sentinel"))));
    }
    #[test]
    fn order_distribution_and_dag_are_bounded(){
        let pairs=dictionary("a*x+b*x+O(y)",&["x"],false);
        assert_eq!(pairs,vec![(p("x"),p("a+b+O(y)"))]);
        let mut e=p("x");for _ in 0..48{e=crate::functions::undefined("pair",vec![e.clone(),e]);}
        let result=run::<Infallible>(&e,&[p("x")],Options::default(),true,None).unwrap();
        let Collection::Expression(result)=result else{unreachable!()};assert!(result.ptr_eq(&e));
        let dummy=Expr::dummy("x",crate::Assumptions::POSITIVE);let other=Expr::dummy("x",crate::Assumptions::POSITIVE);
        let source=add(&[mul(&[p("a"),dummy.clone()]),mul(&[p("b"),dummy.clone()]),other.clone()]);
        let result=run::<Infallible>(&source,&[dummy.clone()],Options::default(),false,None).unwrap();
        let Collection::Dictionary(result)=result else{unreachable!()};assert!(result.contains(&(dummy,p("a+b"))));assert!(result.contains(&(Expr::one(),other)));
    }
}
