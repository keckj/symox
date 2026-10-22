//! Bounded reuse between public expansion calls. Keys use physical identities:
//! alpha-equivalent binders and equal numbers at different precision must keep
//! their own presentation. Retained owners prevent pointer-address reuse.
use std::{cell::RefCell, collections::VecDeque};
use crate::expr::Expr;
use super::ExpandOptions;

const MAX_ENTRIES: usize = 16;
const MAX_WEIGHT: usize = 32_768;

struct Entry {
    input: Expr,
    output: Expr,
    options: ExpandOptions,
    limits: (usize, usize),
    numeric_refinement: bool,
    work: usize,
    weight: usize,
}
#[derive(Default)]
struct Memo { entries: VecDeque<Entry>, weight: usize }
thread_local! { static MEMO: RefCell<Memo> = RefCell::new(Memo::default()); }

pub(super) fn get(e: &Expr, options: &ExpandOptions, limits: (usize, usize)) -> Option<(Expr,usize)> {
    let policy = crate::assumptions::numeric_refinement_enabled();
    MEMO.with(|memo| memo.borrow().entries.iter().rev().find(|entry|
        entry.input.dag_id() == e.dag_id() && entry.options == *options
            && entry.limits == limits && entry.numeric_refinement == policy
    ).map(|entry| (entry.output.clone(),entry.work)))
}

pub(super) fn insert(e: &Expr, output: &Expr, options: &ExpandOptions, limits: (usize, usize), work: usize) {
    let Some(weight)=crate::traverse::retention_weight(&[e,output],MAX_WEIGHT) else{return;};
    let entry = Entry { input:e.clone(), output:output.clone(), options:*options, limits,
        numeric_refinement:crate::assumptions::numeric_refinement_enabled(), work, weight };
    MEMO.with(|memo| {
        let mut memo = memo.borrow_mut();
        while memo.entries.len() >= MAX_ENTRIES || memo.weight + weight > MAX_WEIGHT {
            let old = memo.entries.pop_front().unwrap();
            memo.weight -= old.weight;
        }
        memo.weight += weight;
        memo.entries.push_back(entry);
    });
}
