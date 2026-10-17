//! Propositional encoding used by `ask`/`satask` (SymPy's `cnf.py`).
//!
//! Formulas are put in negation normal form exactly like SymPy's `to_NNF`
//! (relationals become `Q` binary relations, composite unary predicates
//! expand into disjunctions of their components). Conversion to CNF
//! distributes like SymPy for small products and switches to the
//! equisatisfiable Plaisted-Greenbaum encoding for large ones, which keeps
//! the satisfiability answers identical without the exponential blowup.

use rustc_hash::FxHashMap;

use super::facts_data::COMPOSITES;
use super::predicates::Pred;
use super::{applied_unchecked, as_applied};
use crate::expr::{Expr, Func, Kind};
use crate::logic::LogicOp;

/// Interned propositional variables: atoms (applied predicates or other
/// boolean expressions) and auxiliary definition variables.
#[derive(Default)]
pub(super) struct Atoms {
    map: FxHashMap<Expr, u32>,
    list: Vec<Option<Expr>>,
}

impl Atoms {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn var(&mut self, e: &Expr) -> u32 {
        if let Some(v) = self.map.get(e) {
            return *v;
        }
        let v = self.list.len() as u32;
        self.list.push(Some(e.clone()));
        self.map.insert(e.clone(), v);
        v
    }
    fn aux(&mut self) -> u32 {
        let v = self.list.len() as u32;
        self.list.push(None);
        v
    }
    pub fn get(&self, v: u32) -> Option<&Expr> {
        self.list.get(v as usize).and_then(|e| e.as_ref())
    }
    pub fn len(&self) -> usize {
        self.list.len()
    }
}

/// A formula in negation normal form over interned variables.
#[derive(Clone, Debug)]
pub(super) enum F {
    True,
    False,
    Lit(u32, bool),
    And(Vec<F>),
    Or(Vec<F>),
}

impl F {
    pub fn negate(self) -> F {
        match self {
            F::True => F::False,
            F::False => F::True,
            F::Lit(v, n) => F::Lit(v, !n),
            F::And(c) => F::Or(c.into_iter().map(F::negate).collect()),
            F::Or(c) => F::And(c.into_iter().map(F::negate).collect()),
        }
    }
}

fn composite(p: Pred) -> Option<&'static [Pred]> {
    COMPOSITES
        .iter()
        .find(|(c, _)| *c == p)
        .map(|(_, parts)| *parts)
}

/// SymPy's `to_NNF(expr, get_composite_predicates())`.
pub(super) fn nnf(e: &Expr, atoms: &mut Atoms) -> F {
    match e.kind() {
        Kind::Bool(true) => F::True,
        Kind::Bool(false) => F::False,
        Kind::Rel(op, l, r) => {
            let atom = applied_unchecked(Pred::from_relop(*op), vec![l.clone(), r.clone()]);
            F::Lit(atoms.var(&atom), false)
        }
        Kind::Func(Func::Logic(op), args) => {
            let sub = |atoms: &mut Atoms| args.iter().map(|a| nnf(a, atoms)).collect::<Vec<_>>();
            match op {
                LogicOp::Not => nnf(&args[0], atoms).negate(),
                LogicOp::And => F::And(sub(atoms)),
                LogicOp::Or => F::Or(sub(atoms)),
                LogicOp::Implies => {
                    let mut c = sub(atoms);
                    let r = c.pop().unwrap_or(F::True);
                    let l = c.pop().unwrap_or(F::True);
                    F::Or(vec![l.negate(), r])
                }
                LogicOp::Equivalent => {
                    let c = sub(atoms);
                    let n = c.len();
                    F::And(
                        (0..n)
                            .map(|i| F::Or(vec![c[i].clone().negate(), c[(i + 1) % n].clone()]))
                            .collect(),
                    )
                }
                LogicOp::ITE => {
                    let c = sub(atoms);
                    F::And(vec![
                        F::Or(vec![c[0].clone().negate(), c[1].clone()]),
                        F::Or(vec![c[0].clone(), c[2].clone()]),
                    ])
                }
                LogicOp::Xor => {
                    // Balanced parity avoids exponential truth-table expansion
                    // and preserves XOR semantics for more than 16 operands.
                    fn parity(c: &[F]) -> F {
                        match c {
                            [] => F::False,
                            [x] => x.clone(),
                            _ => {
                                let mid = c.len() / 2;
                                let a = parity(&c[..mid]);
                                let b = parity(&c[mid..]);
                                F::Or(vec![
                                    F::And(vec![a.clone(), b.clone().negate()]),
                                    F::And(vec![a.negate(), b]),
                                ])
                            }
                        }
                    }
                    parity(&sub(atoms))
                }
            }
        }
        _ => {
            if let Some((p, args)) = as_applied(e) {
                if args.len() == 1 {
                    if let Some(parts) = composite(p) {
                        return F::Or(
                            parts
                                .iter()
                                .map(|q| {
                                    F::Lit(atoms.var(&applied_unchecked(*q, args.to_vec())), false)
                                })
                                .collect(),
                        );
                    }
                }
            }
            F::Lit(atoms.var(e), false)
        }
    }
}

pub(super) type Clause = Vec<i32>;

#[inline]
pub(super) fn lit(v: u32, negated: bool) -> i32 {
    let l = v as i32 + 1;
    if negated {
        -l
    } else {
        l
    }
}

const DISTRIBUTE_LIMIT: usize = 512;

fn cnf_rec(f: &F, atoms: &mut Atoms, extra: &mut Vec<Clause>) -> Vec<Clause> {
    match f {
        F::True => vec![],
        F::False => vec![vec![]],
        F::Lit(v, n) => vec![vec![lit(*v, *n)]],
        F::And(c) => c.iter().flat_map(|x| cnf_rec(x, atoms, extra)).collect(),
        F::Or(c) => {
            let mut parts: Vec<Vec<Clause>> = Vec::with_capacity(c.len());
            for x in c {
                let p = cnf_rec(x, atoms, extra);
                if p.is_empty() {
                    return vec![];
                }
                parts.push(p);
            }
            let size = parts.iter().try_fold(1usize, |acc, p| {
                acc.checked_mul(p.len()).filter(|s| *s <= DISTRIBUTE_LIMIT)
            });
            if size.is_none() {
                for p in parts.iter_mut() {
                    if p.len() > 1 {
                        let v = atoms.aux();
                        for cl in p.iter() {
                            let mut c2 = Vec::with_capacity(cl.len() + 1);
                            c2.push(lit(v, true));
                            c2.extend_from_slice(cl);
                            extra.push(c2);
                        }
                        *p = vec![vec![lit(v, false)]];
                    }
                }
            }
            let mut acc: Vec<Clause> = vec![vec![]];
            for p in parts {
                let mut next = Vec::with_capacity(acc.len() * p.len());
                for a in &acc {
                    for b in &p {
                        let mut c = a.clone();
                        c.extend_from_slice(b);
                        next.push(c);
                    }
                }
                acc = next;
            }
            acc
        }
    }
}

fn normalize(mut clauses: Vec<Clause>) -> Vec<Clause> {
    for c in clauses.iter_mut() {
        c.sort_unstable();
        c.dedup();
    }
    clauses.sort_unstable();
    clauses.dedup();
    clauses
}

/// Clauses of a formula (a set, like SymPy's `CNF`).
pub(super) fn cnf(f: &F, atoms: &mut Atoms) -> Vec<Clause> {
    let mut extra = Vec::new();
    let mut c = cnf_rec(f, atoms, &mut extra);
    c.extend(extra);
    normalize(c)
}

/// `CNF.from_prop(e)`.
pub(super) fn cnf_of(e: &Expr, atoms: &mut Atoms) -> Vec<Clause> {
    let f = nnf(e, atoms);
    cnf(&f, atoms)
}

/// `CNF.from_prop(~e)`.
pub(super) fn cnf_of_not(e: &Expr, atoms: &mut Atoms) -> Vec<Clause> {
    let f = nnf(e, atoms).negate();
    cnf(&f, atoms)
}

/// Atoms (non auxiliary variables) used by clauses.
pub(super) fn clause_atoms(clauses: &[Clause], atoms: &Atoms) -> Vec<Expr> {
    let mut seen = rustc_hash::FxHashSet::default();
    let mut out = Vec::new();
    for c in clauses {
        for l in c {
            let v = l.unsigned_abs() - 1;
            if seen.insert(v) {
                if let Some(e) = atoms.get(v) {
                    out.push(e.clone());
                }
            }
        }
    }
    out
}

/// Satisfiability of a clause set over `nvars` variables.
pub(super) fn satisfiable(nvars: usize, clauses: &[Clause]) -> bool {
    if clauses.iter().any(|c| c.is_empty()) {
        return false;
    }
    crate::logic::sat::solve_dimacs(nvars, clauses).is_some()
}
