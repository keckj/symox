//! A CDCL SAT solver (two watched literals with blockers, first-UIP
//! clause learning with minimisation, VSIDS, phase saving, Luby restarts,
//! learnt clause database reduction) and the Tseitin encoding of boolean
//! expressions used by [`satisfiable`].

use rustc_hash::FxHashMap;

use super::{as_logic, atoms, normal::is_literal, LogicOp};
use crate::expr::{Expr, Kind};

/// A literal: `2 * var + negated`.
pub type Lit = u32;

#[inline]
pub fn mk_lit(v: u32, negated: bool) -> Lit {
    v << 1 | negated as u32
}
#[inline]
fn var(l: Lit) -> usize {
    (l >> 1) as usize
}
#[inline]
fn neg(l: Lit) -> Lit {
    l ^ 1
}

const L_FALSE: u8 = 0;
const L_TRUE: u8 = 1;
const L_UNDEF: u8 = 2;
const NO_REASON: u32 = u32::MAX;

struct Clause {
    lits: Vec<Lit>,
    learnt: bool,
    deleted: bool,
    act: f32,
}

#[derive(Clone, Copy)]
struct Watcher {
    cref: u32,
    blocker: Lit,
}

/// Indexed binary max-heap of variables keyed by activity.
struct VarHeap {
    heap: Vec<u32>,
    pos: Vec<i32>,
}

impl VarHeap {
    fn contains(&self, v: usize) -> bool {
        self.pos.get(v).is_some_and(|&p| p >= 0)
    }
    fn up(&mut self, mut i: usize, act: &[f64]) {
        let v = self.heap[i];
        while i > 0 {
            let p = (i - 1) / 2;
            if act[self.heap[p] as usize] >= act[v as usize] {
                break;
            }
            self.heap[i] = self.heap[p];
            self.pos[self.heap[i] as usize] = i as i32;
            i = p;
        }
        self.heap[i] = v;
        self.pos[v as usize] = i as i32;
    }
    fn down(&mut self, mut i: usize, act: &[f64]) {
        let v = self.heap[i];
        let n = self.heap.len();
        loop {
            let l = 2 * i + 1;
            if l >= n {
                break;
            }
            let r = l + 1;
            let c = if r < n && act[self.heap[r] as usize] > act[self.heap[l] as usize] { r } else { l };
            if act[self.heap[c] as usize] <= act[v as usize] {
                break;
            }
            self.heap[i] = self.heap[c];
            self.pos[self.heap[i] as usize] = i as i32;
            i = c;
        }
        self.heap[i] = v;
        self.pos[v as usize] = i as i32;
    }
    fn insert(&mut self, v: usize, act: &[f64]) {
        if self.pos.len() <= v {
            self.pos.resize(v + 1, -1);
        }
        if self.contains(v) {
            return;
        }
        self.heap.push(v as u32);
        let i = self.heap.len() - 1;
        self.pos[v] = i as i32;
        self.up(i, act);
    }
    fn pop(&mut self, act: &[f64]) -> Option<usize> {
        if self.heap.is_empty() {
            return None;
        }
        let top = self.heap[0] as usize;
        let last = self.heap.pop().unwrap();
        self.pos[top] = -1;
        if !self.heap.is_empty() {
            self.heap[0] = last;
            self.pos[last as usize] = 0;
            self.down(0, act);
        }
        Some(top)
    }
    fn bumped(&mut self, v: usize, act: &[f64]) {
        if self.contains(v) {
            let i = self.pos[v] as usize;
            self.up(i, act);
        }
    }
}

/// A CDCL SAT solver.
pub struct Solver {
    clauses: Vec<Clause>,
    learnts: Vec<u32>,
    watches: Vec<Vec<Watcher>>,
    assigns: Vec<u8>,
    level: Vec<u32>,
    reason: Vec<u32>,
    trail: Vec<Lit>,
    trail_lim: Vec<usize>,
    qhead: usize,
    activity: Vec<f64>,
    var_inc: f64,
    cla_inc: f32,
    heap: VarHeap,
    polarity: Vec<bool>,
    seen: Vec<u8>,
    ok: bool,
    max_learnts: f64,
    /// Statistics.
    pub conflicts: u64,
    pub decisions: u64,
    pub propagations: u64,
}

impl Default for Solver {
    fn default() -> Self {
        Self::new()
    }
}

impl Solver {
    pub fn new() -> Solver {
        Solver {
            clauses: Vec::new(),
            learnts: Vec::new(),
            watches: Vec::new(),
            assigns: Vec::new(),
            level: Vec::new(),
            reason: Vec::new(),
            trail: Vec::new(),
            trail_lim: Vec::new(),
            qhead: 0,
            activity: Vec::new(),
            var_inc: 1.0,
            cla_inc: 1.0,
            heap: VarHeap { heap: Vec::new(), pos: Vec::new() },
            polarity: Vec::new(),
            seen: Vec::new(),
            ok: true,
            max_learnts: 0.0,
            conflicts: 0,
            decisions: 0,
            propagations: 0,
        }
    }

    /// Create a fresh variable.
    pub fn new_var(&mut self) -> u32 {
        let v = self.assigns.len();
        self.assigns.push(L_UNDEF);
        self.level.push(0);
        self.reason.push(NO_REASON);
        self.activity.push(0.0);
        self.polarity.push(false);
        self.seen.push(0);
        self.watches.push(Vec::new());
        self.watches.push(Vec::new());
        self.heap.insert(v, &self.activity);
        v as u32
    }

    #[inline]
    fn value(&self, l: Lit) -> u8 {
        let a = self.assigns[var(l)];
        if a == L_UNDEF {
            L_UNDEF
        } else {
            a ^ (l & 1) as u8
        }
    }

    /// Value of a variable in the current assignment (after a successful
    /// [`Solver::solve`], the model).
    pub fn var_value(&self, v: u32) -> Option<bool> {
        match self.assigns[v as usize] {
            L_TRUE => Some(true),
            L_FALSE => Some(false),
            _ => None,
        }
    }

    #[inline]
    fn decision_level(&self) -> u32 {
        self.trail_lim.len() as u32
    }

    #[inline]
    fn enqueue(&mut self, l: Lit, from: u32) {
        let v = var(l);
        self.assigns[v] = if l & 1 == 0 { L_TRUE } else { L_FALSE };
        self.level[v] = self.decision_level();
        self.reason[v] = from;
        self.trail.push(l);
    }

    fn attach(&mut self, cref: u32) {
        let c = &self.clauses[cref as usize];
        let (a, b) = (c.lits[0], c.lits[1]);
        self.watches[neg(a) as usize].push(Watcher { cref, blocker: b });
        self.watches[neg(b) as usize].push(Watcher { cref, blocker: a });
    }

    /// Add a clause (only at decision level 0). Returns `false` if the
    /// formula became trivially unsatisfiable.
    pub fn add_clause(&mut self, lits: &[Lit]) -> bool {
        if !self.ok {
            return false;
        }
        if self.decision_level() > 0 {
            self.cancel_until(0);
        }
        let mut ls: Vec<Lit> = lits.to_vec();
        ls.sort_unstable();
        ls.dedup();
        let mut out = Vec::with_capacity(ls.len());
        for (i, &l) in ls.iter().enumerate() {
            if i + 1 < ls.len() && ls[i + 1] == neg(l) {
                return true; // tautology
            }
            match self.value(l) {
                L_TRUE => return true,
                L_FALSE => {}
                _ => out.push(l),
            }
        }
        match out.len() {
            0 => {
                self.ok = false;
                false
            }
            1 => {
                self.enqueue(out[0], NO_REASON);
                if self.propagate().is_some() {
                    self.ok = false;
                }
                self.ok
            }
            _ => {
                let cref = self.clauses.len() as u32;
                self.clauses.push(Clause { lits: out, learnt: false, deleted: false, act: 0.0 });
                self.attach(cref);
                true
            }
        }
    }

    fn propagate(&mut self) -> Option<u32> {
        while self.qhead < self.trail.len() {
            let p = self.trail[self.qhead];
            self.qhead += 1;
            self.propagations += 1;
            let false_lit = neg(p);
            let mut ws = std::mem::take(&mut self.watches[p as usize]);
            let (mut i, mut j) = (0usize, 0usize);
            let n = ws.len();
            let mut confl = None;
            while i < n {
                let w = ws[i];
                if self.value(w.blocker) == L_TRUE {
                    ws[j] = w;
                    j += 1;
                    i += 1;
                    continue;
                }
                let cref = w.cref as usize;
                if self.clauses[cref].deleted {
                    i += 1;
                    continue;
                }
                {
                    let lits = &mut self.clauses[cref].lits;
                    if lits[0] == false_lit {
                        lits.swap(0, 1);
                    }
                }
                i += 1;
                let first = self.clauses[cref].lits[0];
                let nw = Watcher { cref: w.cref, blocker: first };
                if first != w.blocker && self.value(first) == L_TRUE {
                    ws[j] = nw;
                    j += 1;
                    continue;
                }
                // look for a new literal to watch
                let len = self.clauses[cref].lits.len();
                let mut found = false;
                for k in 2..len {
                    let l = self.clauses[cref].lits[k];
                    if self.value(l) != L_FALSE {
                        self.clauses[cref].lits.swap(1, k);
                        self.watches[neg(l) as usize].push(nw);
                        found = true;
                        break;
                    }
                }
                if found {
                    continue;
                }
                ws[j] = nw;
                j += 1;
                if self.value(first) == L_FALSE {
                    confl = Some(w.cref);
                    self.qhead = self.trail.len();
                    while i < n {
                        ws[j] = ws[i];
                        j += 1;
                        i += 1;
                    }
                } else {
                    self.enqueue(first, w.cref);
                }
            }
            ws.truncate(j);
            self.watches[p as usize] = ws;
            if confl.is_some() {
                return confl;
            }
        }
        None
    }

    fn bump_var(&mut self, v: usize) {
        self.activity[v] += self.var_inc;
        if self.activity[v] > 1e100 {
            for a in &mut self.activity {
                *a *= 1e-100;
            }
            self.var_inc *= 1e-100;
        }
        self.heap.bumped(v, &self.activity);
    }

    fn bump_clause(&mut self, cref: u32) {
        let c = &mut self.clauses[cref as usize];
        c.act += self.cla_inc;
        if c.act > 1e20 {
            for &l in &self.learnts {
                self.clauses[l as usize].act *= 1e-20;
            }
            self.cla_inc *= 1e-20;
        }
    }

    fn analyze(&mut self, mut confl: u32) -> (Vec<Lit>, u32) {
        let mut learnt: Vec<Lit> = vec![0];
        let mut path_c = 0i32;
        let mut p: Option<Lit> = None;
        let mut index = self.trail.len();
        let dl = self.decision_level();
        loop {
            if self.clauses[confl as usize].learnt {
                self.bump_clause(confl);
            }
            let start = if p.is_some() { 1 } else { 0 };
            let len = self.clauses[confl as usize].lits.len();
            for k in start..len {
                let q = self.clauses[confl as usize].lits[k];
                let v = var(q);
                if self.seen[v] == 0 && self.level[v] > 0 {
                    self.bump_var(v);
                    self.seen[v] = 1;
                    if self.level[v] >= dl {
                        path_c += 1;
                    } else {
                        learnt.push(q);
                    }
                }
            }
            loop {
                index -= 1;
                if self.seen[var(self.trail[index])] != 0 {
                    break;
                }
            }
            let pl = self.trail[index];
            p = Some(pl);
            confl = self.reason[var(pl)];
            self.seen[var(pl)] = 0;
            path_c -= 1;
            if path_c <= 0 {
                break;
            }
        }
        learnt[0] = neg(p.unwrap());
        // basic minimisation: drop literals implied by others in the clause
        let orig: Vec<Lit> = learnt.clone();
        let mut j = 1;
        for i in 1..learnt.len() {
            let v = var(learnt[i]);
            let r = self.reason[v];
            let keep = if r == NO_REASON {
                true
            } else {
                let c = &self.clauses[r as usize];
                c.lits[1..].iter().any(|&l| self.seen[var(l)] == 0 && self.level[var(l)] > 0)
            };
            if keep {
                learnt[j] = learnt[i];
                j += 1;
            }
        }
        learnt.truncate(j);
        for l in orig {
            self.seen[var(l)] = 0;
        }
        // backtrack level
        let bt = if learnt.len() == 1 {
            0
        } else {
            let mut mi = 1;
            for i in 2..learnt.len() {
                if self.level[var(learnt[i])] > self.level[var(learnt[mi])] {
                    mi = i;
                }
            }
            learnt.swap(1, mi);
            self.level[var(learnt[1])]
        };
        (learnt, bt)
    }

    fn cancel_until(&mut self, lvl: u32) {
        if self.decision_level() <= lvl {
            return;
        }
        let lim = self.trail_lim[lvl as usize];
        for c in (lim..self.trail.len()).rev() {
            let l = self.trail[c];
            let v = var(l);
            self.assigns[v] = L_UNDEF;
            self.reason[v] = NO_REASON;
            self.polarity[v] = l & 1 == 1;
            self.heap.insert(v, &self.activity);
        }
        self.qhead = lim;
        self.trail.truncate(lim);
        self.trail_lim.truncate(lvl as usize);
    }

    fn pick_branch(&mut self) -> Option<Lit> {
        while let Some(v) = self.heap.pop(&self.activity) {
            if self.assigns[v] == L_UNDEF {
                return Some(mk_lit(v as u32, self.polarity[v]));
            }
        }
        None
    }

    fn locked(&self, cref: u32) -> bool {
        let c = &self.clauses[cref as usize];
        let l = c.lits[0];
        self.reason[var(l)] == cref && self.value(l) == L_TRUE
    }

    fn reduce_db(&mut self) {
        let mut ls = std::mem::take(&mut self.learnts);
        ls.sort_by(|&a, &b| {
            let (ca, cb) = (&self.clauses[a as usize], &self.clauses[b as usize]);
            (ca.lits.len() > 2)
                .cmp(&(cb.lits.len() > 2))
                .reverse()
                .then(ca.act.partial_cmp(&cb.act).unwrap_or(std::cmp::Ordering::Equal))
        });
        let half = ls.len() / 2;
        let mut keep = Vec::with_capacity(ls.len());
        for (i, &c) in ls.iter().enumerate() {
            let cl = &self.clauses[c as usize];
            if i < half && cl.lits.len() > 2 && !self.locked(c) {
                let cl = &mut self.clauses[c as usize];
                cl.deleted = true;
                cl.lits = Vec::new();
            } else {
                keep.push(c);
            }
        }
        self.learnts = keep;
    }

    fn luby(y: f64, mut x: u64) -> f64 {
        let (mut size, mut seq) = (1u64, 0u32);
        while size < x + 1 {
            seq += 1;
            size = 2 * size + 1;
        }
        while size - 1 != x {
            size = (size - 1) >> 1;
            seq -= 1;
            x %= size;
        }
        y.powi(seq as i32)
    }

    /// Solve. Returns `true` if satisfiable; the model is then available
    /// through [`Solver::var_value`].
    pub fn solve(&mut self) -> bool {
        if !self.ok {
            return false;
        }
        if self.max_learnts == 0.0 {
            self.max_learnts = (self.clauses.len() as f64 / 3.0).max(2000.0);
        }
        let mut restarts = 0u64;
        loop {
            let limit = (Self::luby(2.0, restarts) * 100.0) as u64;
            match self.search(limit) {
                Some(r) => return r,
                None => restarts += 1,
            }
        }
    }

    fn search(&mut self, nof_conflicts: u64) -> Option<bool> {
        let mut conflict_c = 0u64;
        loop {
            if let Some(confl) = self.propagate() {
                self.conflicts += 1;
                conflict_c += 1;
                if self.decision_level() == 0 {
                    self.ok = false;
                    return Some(false);
                }
                let (learnt, bt) = self.analyze(confl);
                self.cancel_until(bt);
                if learnt.len() == 1 {
                    self.enqueue(learnt[0], NO_REASON);
                } else {
                    let cref = self.clauses.len() as u32;
                    let first = learnt[0];
                    self.clauses.push(Clause { lits: learnt, learnt: true, deleted: false, act: 0.0 });
                    self.attach(cref);
                    self.learnts.push(cref);
                    self.bump_clause(cref);
                    self.enqueue(first, cref);
                }
                self.var_inc /= 0.95;
                self.cla_inc /= 0.999;
            } else {
                if conflict_c >= nof_conflicts {
                    self.cancel_until(0);
                    return None;
                }
                if self.learnts.len() as f64 - self.trail.len() as f64 >= self.max_learnts {
                    self.reduce_db();
                    self.max_learnts *= 1.1;
                }
                match self.pick_branch() {
                    None => return Some(true),
                    Some(l) => {
                        self.decisions += 1;
                        self.trail_lim.push(self.trail.len());
                        self.enqueue(l, NO_REASON);
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------
// boolean expressions -> CNF
// ----------------------------------------------------------------------

/// A satisfying assignment: atoms (in `ordered()` order) and values.
pub type Model = Vec<(Expr, bool)>;

struct Encoder {
    solver: Solver,
    atom_var: FxHashMap<Expr, u32>,
    memo: FxHashMap<Expr, Lit>,
    true_lit: Option<Lit>,
}

impl Encoder {
    fn new(atoms: &[Expr]) -> Encoder {
        let mut solver = Solver::new();
        let mut atom_var = FxHashMap::default();
        for a in atoms {
            let v = solver.new_var();
            atom_var.insert(a.clone(), v);
        }
        Encoder { solver, atom_var, memo: FxHashMap::default(), true_lit: None }
    }

    fn fresh(&mut self) -> Lit {
        mk_lit(self.solver.new_var(), false)
    }

    fn lit_of(&mut self, e: &Expr) -> Lit {
        if let Some(&v) = self.atom_var.get(e) {
            return mk_lit(v, false);
        }
        if let Some(&l) = self.memo.get(e) {
            return l;
        }
        let l = match e.kind() {
            Kind::Bool(b) => {
                let t = match self.true_lit {
                    Some(t) => t,
                    None => {
                        let t = self.fresh();
                        self.solver.add_clause(&[t]);
                        self.true_lit = Some(t);
                        t
                    }
                };
                if *b {
                    t
                } else {
                    neg(t)
                }
            }
            _ => match as_logic(e) {
                Some((op, args)) => {
                    let args = args.to_vec();
                    match op {
                        LogicOp::Not => neg(self.lit_of(&args[0])),
                        LogicOp::And => {
                            let ls: Vec<Lit> = args.iter().map(|a| self.lit_of(a)).collect();
                            self.and_lits(&ls)
                        }
                        LogicOp::Or => {
                            let ls: Vec<Lit> = args.iter().map(|a| self.lit_of(a)).collect();
                            self.or_lits(&ls)
                        }
                        LogicOp::Xor => {
                            let ls: Vec<Lit> = args.iter().map(|a| self.lit_of(a)).collect();
                            let mut acc = ls[0];
                            for &b in &ls[1..] {
                                acc = self.xor2(acc, b);
                            }
                            acc
                        }
                        LogicOp::Implies => {
                            let a = self.lit_of(&args[0]);
                            let b = self.lit_of(&args[1]);
                            self.or_lits(&[neg(a), b])
                        }
                        LogicOp::Equivalent => {
                            let ls: Vec<Lit> = args.iter().map(|a| self.lit_of(a)).collect();
                            let p = self.and_lits(&ls);
                            let nls: Vec<Lit> = ls.iter().map(|&l| neg(l)).collect();
                            let q = self.and_lits(&nls);
                            self.or_lits(&[p, q])
                        }
                        LogicOp::ITE => {
                            let c = self.lit_of(&args[0]);
                            let a = self.lit_of(&args[1]);
                            let b = self.lit_of(&args[2]);
                            let t = self.fresh();
                            let s = &mut self.solver;
                            s.add_clause(&[neg(t), neg(c), a]);
                            s.add_clause(&[neg(t), c, b]);
                            s.add_clause(&[t, neg(c), neg(a)]);
                            s.add_clause(&[t, c, neg(b)]);
                            s.add_clause(&[neg(t), a, b]);
                            s.add_clause(&[t, neg(a), neg(b)]);
                            t
                        }
                    }
                }
                None => {
                    // unexpected atom (not collected): give it a variable
                    let v = self.solver.new_var();
                    self.atom_var.insert(e.clone(), v);
                    mk_lit(v, false)
                }
            },
        };
        self.memo.insert(e.clone(), l);
        l
    }

    fn and_lits(&mut self, ls: &[Lit]) -> Lit {
        if ls.len() == 1 {
            return ls[0];
        }
        let t = self.fresh();
        let mut big = Vec::with_capacity(ls.len() + 1);
        big.push(t);
        for &l in ls {
            self.solver.add_clause(&[neg(t), l]);
            big.push(neg(l));
        }
        self.solver.add_clause(&big);
        t
    }

    fn or_lits(&mut self, ls: &[Lit]) -> Lit {
        if ls.len() == 1 {
            return ls[0];
        }
        let t = self.fresh();
        let mut big = Vec::with_capacity(ls.len() + 1);
        big.push(neg(t));
        for &l in ls {
            self.solver.add_clause(&[t, neg(l)]);
            big.push(l);
        }
        self.solver.add_clause(&big);
        t
    }

    fn xor2(&mut self, a: Lit, b: Lit) -> Lit {
        let t = self.fresh();
        let s = &mut self.solver;
        s.add_clause(&[neg(t), a, b]);
        s.add_clause(&[neg(t), neg(a), neg(b)]);
        s.add_clause(&[t, neg(a), b]);
        s.add_clause(&[t, a, neg(b)]);
        t
    }

    /// Assert `e` (top-level conjunctions and clauses are added directly).
    fn assert_true(&mut self, e: &Expr) {
        match e.kind() {
            Kind::Bool(true) => {}
            Kind::Bool(false) => {
                self.solver.add_clause(&[]);
            }
            _ => match as_logic(e) {
                Some((LogicOp::And, args)) => {
                    for a in args.to_vec() {
                        self.assert_true(&a);
                    }
                }
                Some((LogicOp::Or, args)) => {
                    let args = args.to_vec();
                    let ls: Vec<Lit> = args.iter().map(|a| self.lit_of(a)).collect();
                    self.solver.add_clause(&ls);
                }
                _ => {
                    let l = self.lit_of(e);
                    self.solver.add_clause(&[l]);
                }
            },
        }
    }
}

/// Iterator over the models of a boolean expression.
pub struct ModelIter {
    enc: Encoder,
    atoms: Vec<Expr>,
    done: bool,
}

impl Iterator for ModelIter {
    type Item = Model;
    fn next(&mut self) -> Option<Model> {
        if self.done {
            return None;
        }
        if !self.enc.solver.solve() {
            self.done = true;
            return None;
        }
        let mut model = Vec::with_capacity(self.atoms.len());
        let mut block = Vec::with_capacity(self.atoms.len());
        for a in &self.atoms {
            let v = self.enc.atom_var[a];
            let b = self.enc.solver.var_value(v).unwrap_or(false);
            model.push((a.clone(), b));
            block.push(mk_lit(v, b));
        }
        if block.is_empty() || !self.enc.solver.add_clause(&block) {
            self.done = true;
        }
        Some(model)
    }
}

fn encoder_for(e: &Expr) -> (Encoder, Vec<Expr>) {
    let ats = atoms(e);
    let mut enc = Encoder::new(&ats);
    if let Some((LogicOp::And, args)) = as_logic(e) {
        // CNF fast path: clauses of literals go straight to the solver
        for a in args {
            match as_logic(a) {
                Some((LogicOp::Or, lits)) if lits.iter().all(is_literal) => {
                    let ls: Vec<Lit> = lits.iter().map(|l| enc.lit_of(l)).collect();
                    enc.solver.add_clause(&ls);
                }
                _ => enc.assert_true(a),
            }
        }
    } else {
        enc.assert_true(e);
    }
    (enc, ats)
}

/// All models of `e` (lazily).
pub fn satisfiable_all(e: &Expr) -> ModelIter {
    let (enc, atoms) = encoder_for(e);
    ModelIter { enc, atoms, done: false }
}

/// A model of `e`, or `None` if unsatisfiable.
pub fn satisfiable(e: &Expr) -> Option<Model> {
    satisfiable_all(e).next()
}

/// Solve a CNF given as DIMACS-style clauses (`±(var+1)`). Returns the
/// model (one value per variable) or `None`.
pub fn solve_dimacs(nvars: usize, clauses: &[Vec<i32>]) -> Option<Vec<bool>> {
    let mut s = Solver::new();
    for _ in 0..nvars {
        s.new_var();
    }
    for c in clauses {
        let ls: Vec<Lit> = c.iter().map(|&x| mk_lit(x.unsigned_abs() - 1, x < 0)).collect();
        if !s.add_clause(&ls) {
            return None;
        }
    }
    if s.solve() {
        Some((0..nvars as u32).map(|v| s.var_value(v).unwrap_or(false)).collect())
    } else {
        None
    }
}
