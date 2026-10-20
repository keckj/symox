//! Small exact-value caches. Callers cap integer bit lengths before insertion;
//! values are cloned on lookup, so returned factor maps remain independent.
use std::collections::VecDeque;

pub(super) struct Memo<K,V> { entries: VecDeque<(K,V)> }
impl<K: PartialEq,V: Clone> Memo<K,V> {
    pub(super) fn new() -> Self { Self { entries: VecDeque::new() } }
    pub(super) fn get(&self, key: &K) -> Option<V> {
        self.entries.iter().rev().find(|(stored,_)| stored == key).map(|(_,value)|value.clone())
    }
    pub(super) fn insert(&mut self, key: K, value: V) {
        if self.entries.len() == 32 { self.entries.pop_front(); }
        self.entries.push_back((key,value));
    }
}
