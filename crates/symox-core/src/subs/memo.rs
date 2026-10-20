//! Reuse complete single-substitution results, preserving binding presentation.
use std::{cell::RefCell,collections::VecDeque};
use crate::expr::Expr;

struct Entry { input:Expr, old:Expr, new:Expr, output:Expr, policy:bool, weight:usize }
#[derive(Default)]
struct Memo { entries:VecDeque<Entry>, weight:usize }
thread_local! { static MEMO:RefCell<Memo> = RefCell::new(Memo::default()); }
const MAX_WEIGHT:usize = 32_768;

pub(super) fn get(input:&Expr,old:&Expr,new:&Expr)->Option<Expr>{
    let policy = crate::assumptions::numeric_refinement_enabled();
    MEMO.with(|memo|memo.borrow().entries.iter().rev().find(|entry|
        input.dag_id()==entry.input.dag_id() && entry.policy==policy
            && old.same_presentation(&entry.old) && new.same_presentation(&entry.new))
        .map(|entry|entry.output.clone()))
}

pub(super) fn insert(input:&Expr,old:&Expr,new:&Expr,output:&Expr){
    let Some(weight)=crate::traverse::retention_weight(&[input,old,new,output],MAX_WEIGHT) else{return;};
    let entry = Entry {input:input.clone(),old:old.clone(),new:new.clone(),output:output.clone(),
        policy:crate::assumptions::numeric_refinement_enabled(),weight};
    MEMO.with(|memo|{
        let mut memo=memo.borrow_mut();
        while memo.entries.len()>=16 || memo.weight+weight>MAX_WEIGHT {
            let evicted=memo.entries.pop_front().unwrap();memo.weight-=evicted.weight;
        }
        memo.weight+=weight;memo.entries.push_back(entry);
    });
}
