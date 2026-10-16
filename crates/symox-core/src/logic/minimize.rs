//! Two-level minimisation: truth tables, Quine–McCluskey prime implicants
//! and an exact (branch and bound) minimum cover, `SOPform`, `POSform`
//! and `simplify_logic`.

use rustc_hash::{FxHashMap, FxHashSet};

use super::normal::{is_cnf, is_dnf, is_literal};
use super::{and, as_logic, atoms, not, or, LogicOp};
use crate::expr::{Expr, Func, Kind};

/// Target normal form of [`simplify_logic`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogicForm {
    Cnf,
    Dnf,
}

// ----------------------------------------------------------------------
// bit-parallel truth tables
// ----------------------------------------------------------------------

/// Truth table of `n` variables: bit `r` of the table is the value on row
/// `r`, where variable `i` is bit `n - 1 - i` of `r` (the first variable is
/// the most significant, like SymPy).
#[derive(Clone, Debug)]
pub struct Table {
    pub n: usize,
    pub words: Vec<u64>,
}

impl Table {
    fn nwords(n: usize) -> usize {
        if n <= 6 {
            1
        } else {
            1 << (n - 6)
        }
    }

    fn tail_mask(n: usize) -> u64 {
        if n >= 6 {
            !0
        } else {
            (1u64 << (1u32 << n)) - 1
        }
    }

    fn constant(n: usize, v: bool) -> Table {
        let w = if v { Self::tail_mask(n) } else { 0 };
        Table { n, words: vec![w; Self::nwords(n)] }
    }

    fn var(n: usize, i: usize) -> Table {
        const PAT: [u64; 6] = [
            0xAAAA_AAAA_AAAA_AAAA,
            0xCCCC_CCCC_CCCC_CCCC,
            0xF0F0_F0F0_F0F0_F0F0,
            0xFF00_FF00_FF00_FF00,
            0xFFFF_0000_FFFF_0000,
            0xFFFF_FFFF_0000_0000,
        ];
        let p = n - 1 - i;
        let nw = Self::nwords(n);
        let words = if p < 6 {
            vec![PAT[p] & Self::tail_mask(n); nw]
        } else {
            (0..nw).map(|w| if (w >> (p - 6)) & 1 == 1 { !0 } else { 0 }).collect()
        };
        Table { n, words }
    }

    fn not(mut self) -> Table {
        let m = Self::tail_mask(self.n);
        for w in &mut self.words {
            *w = !*w & m;
        }
        self
    }

    fn zip(mut self, o: &Table, f: impl Fn(u64, u64) -> u64) -> Table {
        for (a, b) in self.words.iter_mut().zip(&o.words) {
            *a = f(*a, *b);
        }
        self
    }

    pub fn get(&self, r: usize) -> bool {
        (self.words[r >> 6] >> (r & 63)) & 1 == 1
    }

    pub fn count(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Rows where the table is true.
    pub fn ones(&self) -> Vec<u32> {
        let mut v = Vec::with_capacity(self.count());
        for (wi, &w) in self.words.iter().enumerate() {
            let mut w = w;
            while w != 0 {
                let b = w.trailing_zeros();
                v.push((wi as u32) << 6 | b);
                w &= w - 1;
            }
        }
        v
    }
}

/// Truth table of `e` over `vars`; `None` if `e` contains other atoms.
pub fn table_of(e: &Expr, vars: &[Expr]) -> Option<Table> {
    let n = vars.len();
    if n > 20 { return None; }
    let idx: FxHashMap<&Expr, usize> = vars.iter().enumerate().map(|(i, v)| (v, i)).collect();
    fn rec(e: &Expr, n: usize, idx: &FxHashMap<&Expr, usize>) -> Option<Table> {
        if let Some(&i) = idx.get(e) {
            return Some(Table::var(n, i));
        }
        Some(match e.kind() {
            Kind::Bool(b) => Table::constant(n, *b),
            Kind::Func(Func::Logic(op), a) => {
                let mut t: Vec<Table> = Vec::with_capacity(a.len());
                for x in a {
                    t.push(rec(x, n, idx)?);
                }
                match op {
                    LogicOp::And => t.into_iter().reduce(|x, y| x.zip(&y, |a, b| a & b))?,
                    LogicOp::Or => t.into_iter().reduce(|x, y| x.zip(&y, |a, b| a | b))?,
                    LogicOp::Xor => t.into_iter().reduce(|x, y| x.zip(&y, |a, b| a ^ b))?,
                    LogicOp::Not => t.pop()?.not(),
                    LogicOp::Implies => {
                        let b = t.pop()?;
                        t.pop()?.not().zip(&b, |a, b| a | b)
                    }
                    LogicOp::Equivalent => {
                        let all = t.iter().cloned().reduce(|x, y| x.zip(&y, |a, b| a & b))?;
                        let none = t.into_iter().reduce(|x, y| x.zip(&y, |a, b| a | b))?.not();
                        all.zip(&none, |a, b| a | b)
                    }
                    LogicOp::ITE => {
                        let c = t.pop()?;
                        let b = t.pop()?;
                        let a = t.pop()?;
                        let na = a.clone().not();
                        a.zip(&b, |x, y| x & y).zip(&na.zip(&c, |x, y| x & y), |x, y| x | y)
                    }
                }
            }
            _ => return None,
        })
    }
    rec(e, n, &idx)
}

// ----------------------------------------------------------------------
// Quine–McCluskey
// ----------------------------------------------------------------------

/// An implicant: the rows `r` with `r & !mask == val` (mask bits are
/// "don't care" positions).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Implicant {
    pub val: u32,
    pub mask: u32,
}

impl Implicant {
    #[inline]
    fn covers(&self, r: u32) -> bool {
        r & !self.mask == self.val
    }
    fn literals(&self, n: usize) -> u32 {
        n as u32 - self.mask.count_ones()
    }
}

/// All prime implicants of the function true on `on` and free on `dc`.
pub fn prime_implicants(n: usize, on: &[u32], dc: &[u32]) -> Vec<Implicant> {
    let mut cur: Vec<Implicant> = on.iter().chain(dc).map(|&v| Implicant { val: v, mask: 0 }).collect();
    cur.sort_unstable();
    cur.dedup();
    let mut primes = Vec::new();
    while !cur.is_empty() {
        let set: FxHashSet<Implicant> = cur.iter().copied().collect();
        let mut used: FxHashSet<Implicant> = FxHashSet::default();
        let mut next: Vec<Implicant> = Vec::new();
        for im in &cur {
            for b in 0..n {
                let bit = 1u32 << b;
                if im.mask & bit != 0 || im.val & bit != 0 {
                    continue;
                }
                let other = Implicant { val: im.val | bit, mask: im.mask };
                if set.contains(&other) {
                    used.insert(*im);
                    used.insert(other);
                    next.push(Implicant { val: im.val, mask: im.mask | bit });
                }
            }
        }
        for im in &cur {
            if !used.contains(im) {
                primes.push(*im);
            }
        }
        next.sort_unstable();
        next.dedup();
        cur = next;
    }
    primes
}

/// Bitset over rows.
#[derive(Clone, PartialEq, Eq)]
struct Bits(Vec<u64>);

impl Bits {
    fn new(n: usize) -> Bits {
        Bits(vec![0; n.div_ceil(64).max(1)])
    }
    fn set(&mut self, i: usize) {
        self.0[i >> 6] |= 1 << (i & 63);
    }
    fn get(&self, i: usize) -> bool {
        (self.0[i >> 6] >> (i & 63)) & 1 == 1
    }
    fn count(&self) -> u32 {
        self.0.iter().map(|w| w.count_ones()).sum()
    }
    fn and_count(&self, o: &Bits) -> u32 {
        self.0.iter().zip(&o.0).map(|(a, b)| (a & b).count_ones()).sum()
    }
    fn is_empty(&self) -> bool {
        self.0.iter().all(|w| *w == 0)
    }
    fn andnot(&mut self, o: &Bits) {
        for (a, b) in self.0.iter_mut().zip(&o.0) {
            *a &= !b;
        }
    }
    fn first(&self) -> Option<usize> {
        for (i, w) in self.0.iter().enumerate() {
            if *w != 0 {
                return Some(i * 64 + w.trailing_zeros() as usize);
            }
        }
        None
    }
    fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter().enumerate().flat_map(|(wi, &w)| {
            let mut w = w;
            std::iter::from_fn(move || {
                if w == 0 {
                    return None;
                }
                let b = w.trailing_zeros() as usize;
                w &= w - 1;
                Some(wi * 64 + b)
            })
        })
    }
}

/// Minimum cover of the rows `on` by `primes` (fewest implicants, then
/// fewest literals). Exact for moderate
/// sizes (branch and bound with a node budget, greedy beyond it).
pub fn min_cover(n: usize, primes: &[Implicant], on: &[u32], pos: bool) -> Vec<Implicant> {
    let nr = on.len();
    if nr == 0 {
        return vec![];
    }
    let cols: Vec<Bits> = primes
        .iter()
        .map(|p| {
            let mut b = Bits::new(nr);
            for (i, &r) in on.iter().enumerate() {
                if p.covers(r) {
                    b.set(i);
                }
            }
            b
        })
        .collect();
    // literal count only: ties are then broken like SymPy, which ignores
    // negations (`pos` only selects the form being minimised)
    let _ = pos;
    let cost: Vec<u64> = primes.iter().map(|p| p.literals(n) as u64).collect();
    let mut alive: Vec<bool> = cols.iter().map(|c| !c.is_empty()).collect();
    let mut rows = Bits::new(nr);
    for i in 0..nr {
        rows.set(i);
    }
    let mut chosen: Vec<usize> = Vec::new();
    // reduction: essential columns and column dominance
    loop {
        let mut changed = false;
        // essential
        let rlist: Vec<usize> = rows.iter().collect();
        for r in rlist {
            if !rows.get(r) {
                continue;
            }
            let mut only = None;
            let mut cnt = 0;
            for (c, col) in cols.iter().enumerate() {
                if alive[c] && col.get(r) {
                    cnt += 1;
                    only = Some(c);
                    if cnt > 1 {
                        break;
                    }
                }
            }
            if cnt == 1 {
                let c = only.unwrap();
                chosen.push(c);
                alive[c] = false;
                rows.andnot(&cols[c]);
                changed = true;
            }
        }
        if rows.is_empty() {
            break;
        }
        // column dominance (restricted to remaining rows)
        let live: Vec<usize> = (0..cols.len()).filter(|&c| alive[c]).collect();
        let restricted: Vec<Bits> = live
            .iter()
            .map(|&c| {
                let mut b = cols[c].clone();
                for (w, r) in b.0.iter_mut().zip(&rows.0) {
                    *w &= r;
                }
                b
            })
            .collect();
        for (i, &c) in live.iter().enumerate() {
            if !alive[c] {
                continue;
            }
            if restricted[i].is_empty() {
                alive[c] = false;
                changed = true;
                continue;
            }
            for (j, &d) in live.iter().enumerate() {
                if i == j || !alive[d] {
                    continue;
                }
                // d dominates c ?
                let sub = restricted[i].0.iter().zip(&restricted[j].0).all(|(a, b)| a & !b == 0);
                if sub && cost[d] <= cost[c] && (restricted[i] != restricted[j] || d < c || cost[d] < cost[c]) {
                    alive[c] = false;
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            break;
        }
    }
    if !rows.is_empty() {
        let live: Vec<usize> = (0..cols.len()).filter(|&c| alive[c]).collect();
        let mut best: Option<(u32, u64, Vec<usize>)> = None;
        let mut budget: u64 = 200_000;
        let mut cur = Vec::new();
        bnb(&cols, &cost, &live, &rows, &mut cur, 0, &mut best, &mut budget);
        match best {
            Some((_, _, v)) => chosen.extend(v),
            None => {
                // greedy fallback
                let mut rem = rows.clone();
                while !rem.is_empty() {
                    let c = *live
                        .iter()
                        .max_by_key(|&&c| (cols[c].and_count(&rem), std::cmp::Reverse(cost[c])))
                        .unwrap();
                    chosen.push(c);
                    rem.andnot(&cols[c]);
                }
            }
        }
    }
    chosen.sort_unstable();
    chosen.dedup();
    chosen.into_iter().map(|c| primes[c]).collect()
}

#[allow(clippy::too_many_arguments)]
fn bnb(
    cols: &[Bits],
    cost: &[u64],
    live: &[usize],
    rows: &Bits,
    cur: &mut Vec<usize>,
    lits: u64,
    best: &mut Option<(u32, u64, Vec<usize>)>,
    budget: &mut u64,
) {
    if *budget == 0 {
        return;
    }
    *budget -= 1;
    if rows.is_empty() {
        let k = cur.len() as u32;
        if best.as_ref().is_none_or(|b| (k, lits) < (b.0, b.1)) {
            *best = Some((k, lits, cur.clone()));
        }
        return;
    }
    // lower bound on the number of further implicants
    let remaining = rows.count();
    let maxc = live.iter().map(|&c| cols[c].and_count(rows)).max().unwrap_or(0);
    if maxc == 0 {
        return;
    }
    let lb = remaining.div_ceil(maxc);
    if let Some(b) = best {
        if cur.len() as u32 + lb > b.0 {
            return;
        }
    }
    // branch on the row with the fewest candidate columns
    let mut pick = rows.first().unwrap();
    let mut fewest = u32::MAX;
    for r in rows.iter() {
        let k = live.iter().filter(|&&c| cols[c].get(r)).count() as u32;
        if k < fewest {
            fewest = k;
            pick = r;
            if k <= 1 {
                break;
            }
        }
    }
    let mut cands: Vec<usize> = live.iter().copied().filter(|&c| cols[c].get(pick)).collect();
    cands.sort_by_key(|&c| (std::cmp::Reverse(cols[c].and_count(rows)), cost[c], c));
    for c in cands {
        let mut nr = rows.clone();
        nr.andnot(&cols[c]);
        cur.push(c);
        bnb(cols, cost, live, &nr, cur, lits + cost[c], best, budget);
        cur.pop();
    }
}

// ----------------------------------------------------------------------
// SOP / POS
// ----------------------------------------------------------------------

fn term_sop(im: &Implicant, vars: &[Expr]) -> Expr {
    let n = vars.len();
    let mut lits = Vec::new();
    for (i, v) in vars.iter().enumerate() {
        let bit = 1u32 << (n - 1 - i);
        if im.mask & bit != 0 {
            continue;
        }
        lits.push(if im.val & bit != 0 { v.clone() } else { not(v) });
    }
    and(lits)
}

fn term_pos(im: &Implicant, vars: &[Expr]) -> Expr {
    let n = vars.len();
    let mut lits = Vec::new();
    for (i, v) in vars.iter().enumerate() {
        let bit = 1u32 << (n - 1 - i);
        if im.mask & bit != 0 {
            continue;
        }
        lits.push(if im.val & bit != 0 { not(v) } else { v.clone() });
    }
    or(lits)
}

/// Minimal sum of products (`SOPform`) of the function true on the rows
/// `minterms` (row bits: first variable most significant) with
/// `dontcares`.
pub fn sop_form(vars: &[Expr], minterms: &[u32], dontcares: &[u32]) -> Expr {
    let n = vars.len();
    let dc: FxHashSet<u32> = dontcares.iter().copied().collect();
    let on: Vec<u32> = minterms.iter().copied().filter(|r| !dc.contains(r)).collect();
    if on.is_empty() {
        return Expr::boolean(false);
    }
    let primes = prime_implicants(n, &on, dontcares);
    let cover = min_cover(n, &primes, &on, false);
    or(cover.iter().map(|im| term_sop(im, vars)).collect())
}

/// Minimal product of sums (`POSform`) of the function true on the rows
/// `minterms` with `dontcares`.
pub fn pos_form(vars: &[Expr], minterms: &[u32], dontcares: &[u32]) -> Expr {
    // SymPy chooses False for an empty on-set, including an entirely
    // unconstrained function (all remaining rows are don't-cares).
    if minterms.is_empty() { return Expr::boolean(false); }
    let n = vars.len();
    let dc: FxHashSet<u32> = dontcares.iter().copied().collect();
    let on: FxHashSet<u32> = minterms.iter().copied().collect();
    let off: Vec<u32> = (0..(1u32 << n)).filter(|r| !on.contains(r) && !dc.contains(r)).collect();
    if off.is_empty() {
        return Expr::boolean(true);
    }
    let primes = prime_implicants(n, &off, dontcares);
    let cover = min_cover(n, &primes, &off, true);
    and(cover.iter().map(|im| term_pos(im, vars)).collect())
}

/// SymPy's `simplify_logic`: a minimal DNF/CNF equivalent of `e`. Without
/// a `form`, the shorter of the two is guessed from the truth table (like
/// SymPy). Expressions with more than 8 variables are returned unchanged
/// unless `force`.
pub fn simplify_logic(e: &Expr, form: Option<LogicForm>, deep: bool, force: bool) -> Expr {
    simplify_logic_dontcare(e, form, deep, force, None)
}

/// Minimise on the care set. Rows where `dontcare` is true can be assigned
/// either value; they are removed from the on-set before minimisation.
pub fn simplify_logic_dontcare(e: &Expr, form: Option<LogicForm>, deep: bool, force: bool, dontcare: Option<&Expr>) -> Expr {
    if let Some(f) = form {
        let ok = match f {
            LogicForm::Cnf => is_cnf(e),
            LogicForm::Dnf => is_dnf(e),
        };
        if ok && e.args().iter().all(is_literal) {
            return e.clone();
        }
    }
    let e = if deep {
        let changes = atoms(e).into_iter()
            .filter(|r| matches!(r.kind(), Kind::Rel(..)))
            .map(|r| { let simplified = crate::simplify::simplify(&r); (r, simplified) })
            .collect();
        crate::subs::xreplace(e, &changes)
    } else { e.clone() };
    if as_logic(&e).is_none() {
        return e;
    }
    // Relationals and their complements must share one propositional atom,
    // including when the complementary predicate occurs only in dontcare.
    let mut all_atoms = atoms(&e);
    if let Some(dc) = dontcare {
        all_atoms.extend(atoms(dc));
        super::ordered(&mut all_atoms);
        all_atoms.dedup();
    }
    let mut repl: FxHashMap<Expr, Expr> = FxHashMap::default();
    let mut undo: FxHashMap<Expr, Expr> = FxHashMap::default();
    let rels: Vec<Expr> = all_atoms.into_iter().filter(|a| matches!(a.kind(), Kind::Rel(..))).collect();
    let mut k = 0;
    for r in &rels {
        if repl.contains_key(r) { continue; }
        let d = loop {
            let candidate = Expr::symbol(&format!("_logic_{k}"));
            k += 1;
            if !crate::traverse::has(&e, &candidate)
                && dontcare.is_none_or(|dc| !crate::traverse::has(dc, &candidate)) { break candidate; }
        };
        undo.insert(d.clone(), r.clone());
        repl.insert(r.clone(), d.clone());
        let nr = super::negated_rel(r);
        if rels.contains(&nr) && !repl.contains_key(&nr) {
            repl.insert(nr, not(&d));
        }
    }
    let ex = crate::subs::xreplace(&e, &repl);
    let mut vars = atoms(&ex);
    if !force && vars.len() > 8 { return e; }
    let mut dc = dontcare.map(|dc| crate::subs::xreplace(dc, &repl));
    if let Some(d) = &dc {
        let mut combined = vars.clone();
        combined.extend(atoms(d));
        super::ordered(&mut combined);
        combined.dedup();
        if !force && combined.len() > 8 {
            dc = None;
        } else {
            vars = combined;
        }
    }
    let Some(tt) = table_of(&ex, &vars) else { return e };
    let dct = dc.as_ref().and_then(|d| table_of(d, &vars));
    let dcrows = dct.as_ref().map(Table::ones).unwrap_or_default();
    let ones: Vec<u32> = tt.ones().into_iter().filter(|&r| dct.as_ref().is_none_or(|d| !d.get(r as usize))).collect();
    let n = vars.len();
    let big = ones.len() >= (1usize << n) / 2 && n > 0 || (n == 0 && !ones.is_empty());
    let r = if form == Some(LogicForm::Dnf) || (form.is_none() && big) {
        sop_form(&vars, &ones, &dcrows)
    } else {
        pos_form(&vars, &ones, &dcrows)
    };
    crate::subs::xreplace(&r, &undo)
}

/// Rows of the truth table of `e` over `vars` (first variable most
/// significant), with the value of `e` on each row.
pub fn truth_table(e: &Expr, vars: &[Expr]) -> Vec<(Vec<bool>, Expr)> {
    let n = vars.len();
    let mut out = Vec::with_capacity(1 << n);
    if let Some(tt) = table_of(e, vars) {
        for r in 0..(1usize << n) {
            let row: Vec<bool> = (0..n).map(|i| (r >> (n - 1 - i)) & 1 == 1).collect();
            out.push((row, Expr::boolean(tt.get(r))));
        }
        return out;
    }
    for r in 0..(1usize << n) {
        let row: Vec<bool> = (0..n).map(|i| (r >> (n - 1 - i)) & 1 == 1).collect();
        let map: FxHashMap<Expr, Expr> = vars.iter().cloned().zip(row.iter().map(|b| Expr::boolean(*b))).collect();
        out.push((row, crate::subs::xreplace(e, &map)));
    }
    out
}

#[cfg(test)]
mod care_tests {
    use super::*;
    use crate::logic::{eval_bool, ite};

    #[test]
    fn care_set_has_priority_and_extra_variables() {
        let (x,y,z) = (Expr::symbol("x"), Expr::symbol("y"), Expr::symbol("z"));
        let e = or(vec![x.clone(), y.clone()]);
        assert_eq!(simplify_logic_dontcare(&e, None, true, false, Some(&y)), x);
        let dc = or(vec![x.clone(),z.clone()]);
        assert_eq!(simplify_logic_dontcare(&e, None, true, false, Some(&dc)), y);
        let e = ite(&x, &y, &z);
        let dc = and(vec![not(&x),y.clone()]);
        for form in [None, Some(LogicForm::Cnf), Some(LogicForm::Dnf)] {
            let result = simplify_logic_dontcare(&e, form, true, false, Some(&dc));
            for row in 0..8 {
                let val = |v: &Expr| [x.clone(),y.clone(),z.clone()].iter().position(|x|x==v).map(|i| row & (1 << i) != 0);
                if eval_bool(&dc,&val)==Some(false) { assert_eq!(eval_bool(&e,&val),eval_bool(&result,&val)); }
            }
        }
    }

    #[test]
    fn all_two_variable_functions_and_care_sets() {
        let vars = [Expr::symbol("a"), Expr::symbol("b")];
        for on in 0u32..16 {
            let rows: Vec<u32> = (0..4).filter(|r| on & (1 << r) != 0).collect();
            let e = sop_form(&vars, &rows, &[]);
            for dc_bits in 0u32..16 {
                let dcrows: Vec<u32> = (0..4).filter(|r| dc_bits & (1 << r) != 0).collect();
                let dc = sop_form(&vars, &dcrows, &[]);
                for form in [None, Some(LogicForm::Cnf), Some(LogicForm::Dnf)] {
                    let result = simplify_logic_dontcare(&e, form, true, false, Some(&dc));
                    let table = table_of(&result, &vars).unwrap();
                    for row in 0..4 {
                        if dc_bits & (1 << row) == 0 {
                            assert_eq!(table.get(row), on & (1 << row) != 0);
                        }
                    }
                }
            }
        }
    }
}
