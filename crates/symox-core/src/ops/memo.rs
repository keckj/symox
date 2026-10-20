//! Completed exact scalar operations, independently bounded on each thread.
//! Floats, binders, dummy/wild/noncommuting symbols and held expressions never
//! enter this cache. Keys own their immutable input graphs. Non-numeric inputs
//! use physical identity, preserving even equal symbols' original declarations.
use std::{cell::RefCell,collections::VecDeque};
use std::hash::{Hash,Hasher};
use rustc_hash::FxHashMap;
use crate::Expr;

const MAX_ENTRIES:usize=512;
const MAX_WEIGHT:usize=65_536;

#[derive(Clone,Copy,PartialEq,Eq,Hash)]
pub(super) enum Operation {Add,Multiply,Power}

#[derive(Clone)]
struct Operand(Expr);
impl PartialEq for Operand{
    fn eq(&self,other:&Self)->bool{
        self.0.ptr_eq(&other.0)||(self.0.is_number()&&other.0.is_number()&&self.0==other.0)
    }
}
impl Eq for Operand{}
impl Hash for Operand{
    fn hash<H:Hasher>(&self,state:&mut H){
        self.0.is_number().hash(state);
        if self.0.is_number(){self.0.hash(state)}else{self.0.dag_id().hash(state)}
    }
}
#[derive(Clone,PartialEq,Eq,Hash)]
struct Key {operation:Operation,a:Operand,b:Operand,policy:bool}
#[derive(Clone,Copy)]
enum Projection {Value,Left,Right}
struct Entry {output:Expr,projection:Projection,weight:usize}
#[derive(Default)]
struct Memo {values:FxHashMap<Key,Entry>,order:VecDeque<Key>,weight:usize}
thread_local!{static MEMO:RefCell<Memo>=RefCell::new(Memo::default());}

pub(super) fn apply(operation:Operation,a:&Expr,b:&Expr,compute:impl FnOnce()->Expr)->Expr{
    let aw=usize::from(a.exact_algebra_weight());
    let bw=usize::from(b.exact_algebra_weight());
    if aw==0||bw==0 || (matches!(operation,Operation::Power)
        && !b.as_number().is_some_and(|n|n.is_integer()&&!n.is_negative())) {
        return compute();
    }
    let key=Key{operation,a:Operand(a.clone()),b:Operand(b.clone()),policy:crate::assumptions::numeric_refinement_enabled()};
    let cached=MEMO.with(|memo|memo.borrow().values.get(&key).map(|entry|match entry.projection{
        Projection::Value=>entry.output.clone(),Projection::Left=>a.clone(),Projection::Right=>b.clone(),
    }));
    if let Some(output)=cached{
        return output;
    }
    // Do not retain a RefCell borrow while arithmetic recursively computes.
    let output=compute();
    let ow=usize::from(output.exact_algebra_weight());
    if ow==0{return output;}
    let weight=aw+bw+ow;
    let projection=if output.ptr_eq(a){Projection::Left}else if output.ptr_eq(b){Projection::Right}else{Projection::Value};
    MEMO.with(|memo|{
        let mut memo=memo.borrow_mut();
        // Recursive arithmetic may have completed the same key already.
        if memo.values.contains_key(&key){return;}
        while memo.values.len()>=MAX_ENTRIES||memo.weight+weight>MAX_WEIGHT{
            let old=memo.order.pop_front().unwrap();
            memo.weight-=memo.values.remove(&old).unwrap().weight;
        }
        memo.weight+=weight;
        memo.order.push_back(key.clone());
        memo.values.insert(key,Entry{output:output.clone(),projection,weight});
    });
    output
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn exact_sum_survives_cache_eviction_and_outgrows_cache_eligibility(){
        let x=Expr::symbol("bounded_binary_sum");
        let mut sum=Expr::zero();
        for i in 1..=400{
            let term=crate::ops::mul2(&Expr::int(i),&crate::ops::pow(&x,&Expr::int(i)));
            sum=crate::ops::add2(&sum,&term);
        }
        let value=crate::subs::subs1(&sum,&x,&Expr::int(2));
        // sum(i*2**i, i=1..n) = (n-1)*2**(n+1)+2, exactly.
        let expected=(num_bigint::BigInt::from(399)<<401usize)+2;
        assert_eq!(value,Expr::num(crate::number::Number::from_bigint(expected)));
        MEMO.with(|memo|{
            let memo=memo.borrow();
            assert!(memo.values.len()<=MAX_ENTRIES);
            assert!(memo.weight<=MAX_WEIGHT);
            assert_eq!(memo.values.len(),memo.order.len());
            assert_eq!(memo.weight,memo.values.values().map(|entry|entry.weight).sum::<usize>());
        });
    }
}
