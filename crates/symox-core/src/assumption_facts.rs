//! Immutable, exact declarations and closure of the standard Symbol facts.
//! This is separate from structural expression inference and numerical queries.
#[path = "assumption_facts_data.rs"]
mod data;

/// A finite three-valued fact map. Absence means unknown, never false.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct FactSet { known: u32, positive: u32 }
impl FactSet {
    pub const EMPTY: Self = Self { known: 0, positive: 0 };
    pub const COMMUTATIVE:Self=Self{known:1<<2,positive:1<<2};
    pub const NAMES: &'static [&'static str] = &data::NAMES;
    pub fn from_bits(known:u32,positive:u32)->Result<Self,String> {
        if known>>31!=0 || positive & !known!=0 {return Err("invalid Symbol fact masks".into());}
        let value=Self{known,positive};
        if value.get("commutative")==Some(false){return Err("noncommutative symbols are not implemented".into());}
        value.closure()?;Ok(value)
    }
    pub fn bits(self)->(u32,u32){(self.known,self.positive)}
    pub fn from_kwargs<'a>(kwargs:impl IntoIterator<Item=(&'a str,bool)>)->Result<Self,String>{
        let mut facts=Self::EMPTY;
        for (name,value) in kwargs {
            let index=data::NAMES.iter().position(|n|*n==name).ok_or_else(||format!("unknown assumption {name}"))?;
            if name=="commutative"&&!value{return Err("noncommutative symbols are not implemented".into());}
            facts.merge(1<<index,if value{1<<index}else{0})?;
        }
        facts.closure()?;Ok(facts)
    }
    pub fn get(self,name:&str)->Option<bool>{
        let index=data::NAMES.iter().position(|n|*n==name)?;
        let bit=1<<index;
        (self.known&bit!=0).then_some(self.positive&bit!=0)
    }
    pub fn iter(self)->impl Iterator<Item=(&'static str,bool)>{
        data::NAMES.iter().enumerate().filter_map(move |(i,name)|{
            let bit=1<<i;(self.known&bit!=0).then_some((*name,self.positive&bit!=0))
        })
    }
    pub fn len(self)->usize{self.known.count_ones() as usize}
    pub fn is_empty(self)->bool{self.known==0}
    fn merge(&mut self,known:u32,positive:u32)->Result<(),String>{
        if self.known&known&(self.positive^positive)!=0{return Err("inconsistent assumptions".into());}
        self.known|=known;self.positive|=positive;Ok(())
    }
    pub fn closure(self)->Result<Self,String>{
        if self.is_empty()||self==Self::COMMUTATIVE{return Ok(Self::COMMUTATIVE);}
        let mut facts=self;
        let commutative=1<<data::NAMES.iter().position(|n|*n=="commutative").unwrap();
        facts.merge(commutative,commutative)?;
        loop {
            let before=facts;
            for i in 0..31 {
                let bit=1<<i;
                if facts.known&bit!=0 {
                    let (known,positive)=data::ALPHA[2*i+usize::from(facts.positive&bit!=0)];
                    facts.merge(known,positive)?;
                }
            }
            for (known,positive,conclusion,truth) in data::BETA {
                if facts.known&known==known && facts.positive&known==positive {facts.merge(conclusion,truth)?;}
            }
            if before==facts {return Ok(facts);}
        }
    }
    pub fn from_legacy(flags:crate::Assumptions)->Self{
        if flags.is_empty(){return Self::EMPTY;}
        Self::from_kwargs(crate::Assumptions::NAMES.iter().filter(|(_,flag)|flags.has(*flag)).map(|(name,_)|{
            name.strip_prefix('!').map_or((*name,true),|name|(name,false))
        })).expect("consistent legacy Symbol assumptions")
    }
    pub fn legacy(self)->crate::Assumptions{
        if self.is_empty()||self==Self::COMMUTATIVE{return crate::Assumptions::NONE;}
        let mut flags=crate::Assumptions::NONE;
        for (name,flag) in crate::Assumptions::NAMES {
            let (name,value)=name.strip_prefix('!').map_or((*name,true),|name|(name,false));
            if self.get(name)==Some(value){flags=flags.with(*flag);}
        }
        flags
    }
    pub fn keywords(self)->Vec<String>{self.iter().map(|(name,value)|format!("{name}={}",if value{"True"}else{"False"})).collect()}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn closure_preserves_false_and_beta_rules(){
        let p=FactSet::from_kwargs([("even",true),("positive",true),("prime",false)]).unwrap().closure().unwrap();
        assert_eq!(p.get("composite"),Some(true));assert_eq!(p.get("algebraic"),Some(true));
        let p=FactSet::from_kwargs([("prime",false)]).unwrap().closure().unwrap();
        assert_eq!(p.len(),2);assert_eq!(p.get("integer"),None);
        let p=FactSet::from_kwargs([("extended_real",true),("finite",false)]).unwrap().closure().unwrap();
        assert_eq!(p.get("real"),Some(false));assert_eq!(p.get("extended_nonzero"),Some(true));
        assert!(FactSet::from_kwargs([("integer",true),("even",false),("odd",false)]).is_err());
        assert!(FactSet::from_kwargs([("real",true),("finite",false)]).is_err());
    }
    #[test]
    fn symbol_identity_projection_and_bound_facts_are_lossless(){
        use crate::{Symbol,Expr,Kind};
        let original=Symbol::with_facts("x",FactSet::from_kwargs([("prime",false)]).unwrap()).unwrap();
        let same=Symbol::with_facts("x",FactSet::from_kwargs([("commutative",true),("prime",false)]).unwrap()).unwrap();
        assert_eq!(original,same);assert_ne!(original,Symbol::new("x"));
        let dummy=original.fresh_dummy_like("integration");assert_eq!(dummy.facts(),original.facts());
        let extended=Symbol::with_facts("r",FactSet::from_kwargs([("extended_real",true)]).unwrap()).unwrap();
        assert!(!extended.assumptions().has(crate::Assumptions::REAL));
        let e=Expr::from_kind(Kind::Sym(extended));assert_eq!(crate::assumptions::is_real(&e),None);
        let even=Expr::from_kind(Kind::Sym(Symbol::with_facts("e",FactSet::from_kwargs([("even",false)]).unwrap()).unwrap()));
        assert_eq!(crate::assumptions::is_even(&even),Some(false));assert_eq!(crate::assumptions::is_odd(&even),None);
    }
    #[test]
    fn exact_solver_checks_include_negative_declarations_and_algebraic_roots(){
        use crate::{Expr,Kind,Symbol};
        for (kwargs,constant,expected) in [
            (vec![("prime",false)],4,vec![Expr::int(-2)]),
            (vec![("even",false)],4,vec![]),
            (vec![("imaginary",false)],-1,vec![]),
            (vec![("transcendental",true)],2,vec![]),
        ]{
            let x=Expr::from_kind(Kind::Sym(Symbol::with_facts("x",FactSet::from_kwargs(kwargs).unwrap()).unwrap()));
            let f=crate::ops::sub(&crate::ops::powi(&x,2),&Expr::int(constant));
            assert_eq!(crate::solve::solve(&f,&x).unwrap(),expected);
        }
    }
    #[test]
    fn algebraicity_does_not_materialize_compact_numeric_powers(){
        use crate::{Expr,Kind,Number};
        let value=Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,vec![(Expr::int(2),Expr::int(1_000_000_000))]));
        assert_eq!(crate::assumptions::fact(&value,"algebraic"),Some(true));
        assert_eq!(crate::assumptions::fact(&value,"transcendental"),Some(false));
        assert!(!value.is_canonical());
    }
    #[test]
    fn extended_symbol_facts_control_zero_and_infinite_coefficients(){
        use crate::{Expr,Kind};
        let symbol=|kwargs:Vec<(&str,bool)>|Expr::from_kind(Kind::Sym(crate::Symbol::with_facts("x",FactSet::from_kwargs(kwargs).unwrap()).unwrap()));
        let infinite=symbol(vec![("finite",false)]);
        assert!(infinite.has_nonfinite());
        assert_eq!(crate::ops::mul2(&Expr::zero(),&infinite),Expr::nan());
        assert_eq!(crate::ops::mul2(&infinite,&Expr::zero()),Expr::nan());
        assert!(crate::ops::mul2(&Expr::zoo(),&infinite).is_mul());
        let positive=symbol(vec![("extended_positive",true),("finite",false)]);
        assert_eq!(crate::ops::mul2(&Expr::oo(),&positive),Expr::oo());
        assert_eq!(crate::ops::mul2(&Expr::neg_oo(),&positive),Expr::neg_oo());
        assert_eq!(crate::ops::mul2(&Expr::zoo(),&positive),Expr::zoo());
    }
}
