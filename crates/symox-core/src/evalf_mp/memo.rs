//! Completed numeric evaluations keyed by input identity and bit precision.
use std::{cell::RefCell,collections::VecDeque};
use crate::expr::Expr;

struct Entry { input:Expr, output:Expr, precision:u32, policy:bool, residue:bool, weight:usize }
#[derive(Default)]
struct Memo { entries:VecDeque<Entry>, weight:usize }
thread_local! { static MEMO:RefCell<Memo> = RefCell::new(Memo::default()); }
const MAX_WEIGHT:usize = 16_384;

pub(super) fn get(input:&Expr,precision:u32,residue:bool)->Option<Expr>{
    let policy=crate::assumptions::numeric_refinement_enabled();
    MEMO.with(|memo|memo.borrow().entries.iter().rev().find(|entry|
        entry.input.dag_id()==input.dag_id() && entry.precision==precision && entry.policy==policy && entry.residue==residue)
        .map(|entry|entry.output.clone()))
}

pub(super) fn insert(input:&Expr,output:&Expr,precision:u32,residue:bool){
    // Bound mantissa storage as well as expression-node retention.
    if precision>16_384 {return;}
    let Some(weight)=crate::traverse::retention_weight(&[input,output],MAX_WEIGHT) else{return;};
    let entry=Entry {input:input.clone(),output:output.clone(),precision,
        policy:crate::assumptions::numeric_refinement_enabled(),residue,weight};
    MEMO.with(|memo|{
        let mut memo=memo.borrow_mut();
        while memo.entries.len()>=16 || memo.weight+weight>MAX_WEIGHT {
            let evicted=memo.entries.pop_front().unwrap();memo.weight-=evicted.weight;
        }
        memo.weight+=weight;memo.entries.push_back(entry);
    });
}
