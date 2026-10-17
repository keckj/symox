//! The `ask` pipeline (SymPy 1.14 `sympy/assumptions/ask.py`) and
//! `satask` with the class facts of `sathandlers.py`.

use std::cell::{Cell, RefCell};

use rustc_hash::{FxHashMap, FxHashSet};

use super::cnf::{self, Atoms, Clause, F};
use super::facts_data::{ALL_FACTS, COMPOSITES, FACTS_DICT, MATRIX_FACTS, NUMBER_FACTS};
use super::predicates::Pred;
use super::{applied_unchecked, as_applied, handlers, lra, q, AskError};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::logic::LogicOp;

pub(super) type R = Result<Option<bool>, AskError>;

/// Optional native evaluator for undefined `Predicate(name)` handlers:
/// `(name, argument, assumptions)`.
pub type ExternalHandler<'a> = dyn Fn(&str, &Expr, &Expr) -> Result<Option<bool>, AskError> + 'a;

/// Options of [`ask_with`].
#[derive(Default)]
pub struct AskOptions<'a> {
    /// Extra assumptions, like SymPy's `global_assumptions` context.
    pub context: Vec<Expr>,
    /// Handler of undefined predicates.
    pub external: Option<&'a ExternalHandler<'a>>,
}

/// `ask(proposition, assumptions)` with an empty context.
pub fn ask(prop: &Expr, assumptions: &Expr) -> Result<Option<bool>, AskError> {
    ask_with(prop, assumptions, &AskOptions::default())
}

/// `ask(proposition, assumptions, context)`.
pub fn ask_with(
    prop: &Expr,
    assumptions: &Expr,
    opts: &AskOptions<'_>,
) -> Result<Option<bool>, AskError> {
    let cx = Cx::new(&opts.context, opts.external);
    cx.ask(prop, assumptions)
}

/// Query state shared by the handlers of one top level `ask`.
pub(crate) struct Cx<'a> {
    pub context: &'a [Expr],
    external: Option<&'a ExternalHandler<'a>>,
    memo: RefCell<FxHashMap<(Expr, Expr), Option<bool>>>,
    active: RefCell<FxHashSet<(Expr, Expr)>>,
    depth: Cell<u32>,
}

const MAX_DEPTH: u32 = 60;

/// SymPy `BooleanKind`.
pub(super) fn is_boolean_kind(e: &Expr) -> bool {
    match e.kind() {
        Kind::Bool(_) | Kind::Rel(..) => true,
        Kind::Func(Func::Logic(_) | Func::AppliedPredicate, _) => true,
        Kind::Func(Func::Set(crate::sets::SetOp::Contains), _) => true,
        _ => false,
    }
}

/// SymPy `expr.kind is MatrixKind`.
pub(super) fn is_matrix(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::MatrixProxy(_), _))
}

impl<'a> Cx<'a> {
    pub(crate) fn new(context: &'a [Expr], external: Option<&'a ExternalHandler<'a>>) -> Self {
        Cx {
            context,
            external,
            memo: RefCell::default(),
            active: RefCell::default(),
            depth: Cell::new(0),
        }
    }

    pub(super) fn external(&self, p: Pred, arg: &Expr, a: &Expr) -> R {
        match self.external {
            Some(h) => h(p.name(), arg, a),
            None => Ok(None),
        }
    }

    /// `ask(Q.p(e), assumptions)`.
    pub(super) fn askp(&self, p: Pred, e: &Expr, a: &Expr) -> R {
        self.ask(&q(p, e), a)
    }

    /// `ask(Q.p(e), assumptions) is True` (Python truthiness).
    pub(super) fn is(&self, p: Pred, e: &Expr, a: &Expr) -> Result<bool, AskError> {
        Ok(self.askp(p, e, a)? == Some(true))
    }

    /// The full `ask` pipeline.
    pub(crate) fn ask(&self, prop: &Expr, a: &Expr) -> R {
        if !is_boolean_kind(prop) {
            return Err(AskError::Type(
                "proposition must be a valid logical expression".into(),
            ));
        }
        if !is_boolean_kind(a) {
            return Err(AskError::Type(
                "assumptions must be a valid logical expression".into(),
            ));
        }
        let key = (prop.clone(), a.clone());
        if let Some(v) = self.memo.borrow().get(&key) {
            return Ok(*v);
        }
        if self.depth.get() >= MAX_DEPTH || self.active.borrow().contains(&key) {
            return Ok(None);
        }
        self.active.borrow_mut().insert(key.clone());
        self.depth.set(self.depth.get() + 1);
        let r = self.ask_uncached(prop, a);
        self.depth.set(self.depth.get() - 1);
        self.active.borrow_mut().remove(&key);
        if let Ok(v) = &r {
            self.memo.borrow_mut().insert(key, *v);
        }
        r
    }

    fn ask_uncached(&self, prop: &Expr, a: &Expr) -> R {
        let (key, args): (Pred, Vec<Expr>) = if let Some((p, args)) = as_applied(prop) {
            (p, args.to_vec())
        } else if let Kind::Rel(op, l, r) = prop.kind() {
            (Pred::from_relop(*op), vec![l.clone(), r.clone()])
        } else {
            (Pred::IsTrue, vec![prop.clone()])
        };
        let mut atoms = Atoms::new();
        let mut clauses = cnf::cnf_of(a, &mut atoms);
        for c in self.context {
            clauses.extend(cnf::cnf_of(c, &mut atoms));
        }
        let local = extract_all_facts(&clauses, &atoms, &args);
        if !local.is_empty() && !facts_consistent(&local) {
            return Err(AskError::Value(format!(
                "inconsistent assumptions {}",
                crate::printer::to_str(a)
            )));
        }
        if let Some(v) = ask_single_fact(key, &local) {
            return Ok(Some(v));
        }
        if let Some(v) = handlers::eval(self, key, &args, a)? {
            return Ok(Some(v));
        }
        if let Some(v) = self.satask(prop, a)? {
            return Ok(Some(v));
        }
        lra::lra_satask(self, prop, a)
    }

    /// `satask(proposition, assumptions, context)`.
    fn satask(&self, prop: &Expr, a: &Expr) -> R {
        let mut atoms = Atoms::new();
        let props = cnf::cnf_of(prop, &mut atoms);
        let nprops = cnf::cnf_of_not(prop, &mut atoms);
        let acl = cnf::cnf_of(a, &mut atoms);
        let mut ccl = Vec::new();
        for c in self.context {
            ccl.extend(cnf::cnf_of(c, &mut atoms));
        }
        // extract_predargs
        let keys = cnf::clause_atoms(&props, &atoms);
        let mut req: FxHashSet<Expr> = FxHashSet::default();
        for k in &keys {
            find_symbols(k, &mut req);
        }
        let mut lkeys = cnf::clause_atoms(&acl, &atoms);
        for k in cnf::clause_atoms(&ccl, &atoms) {
            if !lkeys.contains(&k) {
                lkeys.push(k);
            }
        }
        lkeys.retain(|k| !matches!(k.kind(), Kind::Bool(_)));
        let lsyms: Vec<FxHashSet<Expr>> = lkeys
            .iter()
            .map(|k| {
                let mut s = FxHashSet::default();
                find_symbols(k, &mut s);
                s
            })
            .collect();
        loop {
            let mut grew = false;
            for s in &lsyms {
                if s.iter().any(|x| req.contains(x)) {
                    for x in s {
                        if req.insert(x.clone()) {
                            grew = true;
                        }
                    }
                }
            }
            if !grew {
                break;
            }
        }
        let mut all_keys = keys;
        for (k, s) in lkeys.iter().zip(&lsyms) {
            if s.iter().any(|x| req.contains(x)) && !all_keys.contains(k) {
                all_keys.push(k.clone());
            }
        }
        let mut exprs: Vec<Expr> = Vec::new();
        let mut seen: FxHashSet<Expr> = FxHashSet::default();
        for k in &all_keys {
            match as_applied(k) {
                Some((_, args)) => {
                    for x in args {
                        if seen.insert(x.clone()) {
                            exprs.push(x.clone());
                        }
                    }
                }
                None => {
                    if seen.insert(k.clone()) {
                        exprs.push(k.clone());
                    }
                }
            }
        }
        // get_all_relevant_facts
        let mut sat: Vec<Clause> = Vec::new();
        let mut all_exprs = exprs.clone();
        let mut frontier = exprs;
        let mut rounds = 0;
        while !frontier.is_empty() && rounds < 64 {
            rounds += 1;
            let mut next = Vec::new();
            for e in &frontier {
                let (facts, newexprs) = class_facts(e, &mut atoms);
                for f in facts {
                    sat.extend(cnf::cnf(&f, &mut atoms));
                }
                for x in newexprs {
                    if seen.insert(x.clone()) {
                        next.push(x.clone());
                        all_exprs.push(x);
                    }
                }
            }
            frontier = next;
        }
        let any_matrix = all_exprs.iter().any(is_matrix);
        let any_number = all_exprs.iter().any(|e| {
            !is_matrix(e)
                && !is_boolean_kind(e)
                && !matches!(e.kind(), Kind::Func(Func::Logic(_), _))
        });
        for e in &all_exprs {
            if any_matrix {
                add_known(&mut sat, MATRIX_FACTS, e, &mut atoms);
            }
            if any_number {
                add_known(&mut sat, NUMBER_FACTS, e, &mut atoms);
            }
        }
        sat.extend(acl);
        sat.extend(ccl);
        let n = atoms.len();
        let mut t = sat.clone();
        t.extend(props);
        let mut f = sat;
        f.extend(nprops);
        let can_true = cnf::satisfiable(n, &t);
        let can_false = cnf::satisfiable(n, &f);
        match (can_true, can_false) {
            (true, true) => Ok(None),
            (true, false) => Ok(Some(true)),
            (false, true) => Ok(Some(false)),
            (false, false) => Err(AskError::Value("Inconsistent assumptions".into())),
        }
    }
}

fn add_known(sat: &mut Vec<Clause>, facts: &[&[(Pred, bool)]], e: &Expr, atoms: &mut Atoms) {
    for cl in facts {
        sat.push(
            cl.iter()
                .map(|(p, n)| cnf::lit(atoms.var(&applied_unchecked(*p, vec![e.clone()])), *n))
                .collect(),
        );
    }
}

/// `pred.atoms(Symbol)`: symbols, including those of matrix shapes.
pub(super) fn find_symbols(e: &Expr, out: &mut FxHashSet<Expr>) {
    // Retain visited owners: public args may synthesize temporary nodes whose
    // allocation addresses can be reused during the same traversal.
    fn visit(e: &Expr, out: &mut FxHashSet<Expr>, seen: &mut FxHashSet<Expr>) {
        if !seen.insert(e.clone()) {
            return;
        }
        match e.kind() {
            Kind::Sym(_) => {
                out.insert(e.clone());
            }
            Kind::Func(
                Func::MatrixProxy(m)
                | Func::MatrixElement(m)
                | Func::MatrixTrace(m)
                | Func::MatrixDeterminant(m),
                args,
            ) => {
                for s in m.scalar_free_symbols() {
                    out.insert(s);
                }
                for d in [m.rows(), m.cols()] {
                    for x in crate::traverse::free_symbols(d) {
                        out.insert(x);
                    }
                }
                for x in args {
                    visit(x, out, seen);
                }
            }
            _ => {
                for x in e.args().iter() {
                    visit(x, out, seen);
                }
            }
        }
    }
    visit(e, out, &mut FxHashSet::default());
}

/// `_extract_all_facts(assump, exprs)`: unit predicates about `exprs`.
fn extract_all_facts(clauses: &[Clause], atoms: &Atoms, exprs: &[Expr]) -> Vec<Vec<(Pred, bool)>> {
    let mut out: Vec<Vec<(Pred, bool)>> = Vec::new();
    'clause: for c in clauses {
        let mut lits = Vec::new();
        for l in c {
            let v = l.unsigned_abs() - 1;
            let Some(atom) = atoms.get(v) else {
                continue 'clause;
            };
            match as_applied(atom) {
                Some((p, args)) if args.len() == 1 && exprs.contains(&args[0]) => {
                    lits.push((p, *l < 0))
                }
                _ => continue 'clause,
            }
        }
        if !lits.is_empty() {
            lits.sort_by_key(|(p, n)| (p.name(), *n));
            lits.dedup();
            if !out.contains(&lits) {
                out.push(lits);
            }
        }
    }
    out
}

/// Satisfiability of the known facts together with the local facts, with
/// one variable per predicate.
fn facts_consistent(local: &[Vec<(Pred, bool)>]) -> bool {
    let mut vars: FxHashMap<Pred, u32> = FxHashMap::default();
    let mut var = |p: Pred| {
        let n = vars.len() as u32;
        *vars.entry(p).or_insert(n)
    };
    let mut clauses: Vec<Clause> = Vec::new();
    for cl in ALL_FACTS
        .iter()
        .copied()
        .chain(local.iter().map(|c| c.as_slice()))
    {
        clauses.push(cl.iter().map(|(p, n)| cnf::lit(var(*p), *n)).collect());
    }
    cnf::satisfiable(vars.len(), &clauses)
}

fn facts_entry(p: Pred) -> Option<&'static (Pred, &'static [Pred], &'static [Pred])> {
    FACTS_DICT.iter().find(|(k, _, _)| *k == p)
}

/// `_ask_single_fact(key, local_facts)`.
fn ask_single_fact(key: Pred, local: &[Vec<(Pred, bool)>]) -> Option<bool> {
    if local.is_empty() {
        return None;
    }
    if local.len() == 1 && local[0].len() == 1 {
        let (f, neg) = local[0][0];
        if neg {
            if let Some((_, req, _)) = facts_entry(key) {
                if req.contains(&f) {
                    return Some(false);
                }
            }
        }
    }
    for c in local {
        if c.len() == 1 && !c[0].1 {
            if let Some((_, req, rej)) = facts_entry(c[0].0) {
                if req.contains(&key) {
                    return Some(true);
                }
                if rej.contains(&key) {
                    return Some(false);
                }
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// class facts (sathandlers.py)
// ---------------------------------------------------------------------------

/// The literal `Q.p(arg)` with composite predicates expanded like `to_NNF`.
fn fq(p: Pred, arg: &Expr, atoms: &mut Atoms) -> F {
    if let Some((_, parts)) = COMPOSITES.iter().find(|(c, _)| *c == p) {
        return F::Or(
            parts
                .iter()
                .map(|x| F::Lit(atoms.var(&q(*x, arg)), false))
                .collect(),
        );
    }
    F::Lit(atoms.var(&q(p, arg)), false)
}

fn implies(a: F, b: F) -> F {
    F::Or(vec![a.negate(), b])
}

fn equiv(a: F, b: F) -> F {
    F::And(vec![
        F::Or(vec![a.clone().negate(), b.clone()]),
        F::Or(vec![b.negate(), a]),
    ])
}

fn allargs(args: &[Expr], atoms: &mut Atoms, f: &mut dyn FnMut(&Expr, &mut Atoms) -> F) -> F {
    F::And(args.iter().map(|x| f(x, atoms)).collect())
}

fn anyarg(args: &[Expr], atoms: &mut Atoms, f: &mut dyn FnMut(&Expr, &mut Atoms) -> F) -> F {
    F::Or(args.iter().map(|x| f(x, atoms)).collect())
}

fn exactlyonearg(args: &[Expr], atoms: &mut Atoms, f: &mut dyn FnMut(&Expr, &mut Atoms) -> F) -> F {
    let lits: Vec<F> = args.iter().map(|x| f(x, atoms)).collect();
    F::Or(
        (0..lits.len())
            .map(|i| {
                let mut conj = vec![lits[i].clone()];
                for (j, l) in lits.iter().enumerate() {
                    if j != i {
                        conj.push(l.clone().negate());
                    }
                }
                F::And(conj)
            })
            .collect(),
    )
}

/// SymPy `Number`, `NumberSymbol` or `ImaginaryUnit`.
fn is_number_class(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => !matches!(c, Constant::ComplexInfinity),
        _ => false,
    }
}

fn class_facts(e: &Expr, atoms: &mut Atoms) -> (Vec<F>, Vec<Expr>) {
    use Pred::*;
    let mut out = Vec::new();
    if is_number_class(e) {
        for (p, name) in [
            (Positive, "positive"),
            (Zero, "zero"),
            (Negative, "negative"),
            (Rational, "rational"),
            (Irrational, "irrational"),
            (Even, "even"),
            (Odd, "odd"),
            (Imaginary, "imaginary"),
            (Prime, "prime"),
            (Composite, "composite"),
        ] {
            if let Some(v) = crate::assumptions::fact(e, name) {
                let l = fq(p, e, atoms);
                out.push(if v { l } else { l.negate() });
            }
        }
        return (out, Vec::new());
    }
    let p = |pr: Pred| move |x: &Expr, atoms: &mut Atoms| fq(pr, x, atoms);
    let np = |pr: Pred| move |x: &Expr, atoms: &mut Atoms| fq(pr, x, atoms).negate();
    match e.kind() {
        Kind::Func(Func::Abs, a) => {
            let arg = &a[0];
            out.push(fq(Nonnegative, e, atoms));
            let l = fq(Zero, arg, atoms).negate();
            let r = fq(Zero, e, atoms).negate();
            out.push(equiv(l, r));
            for pr in [Even, Odd, Integer] {
                let l = fq(pr, arg, atoms);
                let r = fq(pr, e, atoms);
                out.push(implies(l, r));
            }
            (out, vec![arg.clone()])
        }
        Kind::Add(..) => {
            let args = crate::traverse::ordered_args(e);
            for pr in [Positive, Negative, Real, Rational, Integer] {
                let l = allargs(&args, atoms, &mut p(pr));
                let r = fq(pr, e, atoms);
                out.push(implies(l, r));
            }
            let l = exactlyonearg(&args, atoms, &mut np(Integer));
            let r = fq(Integer, e, atoms).negate();
            out.push(implies(l, r));
            let ar = allargs(&args, atoms, &mut p(Real));
            let oi = exactlyonearg(&args, atoms, &mut p(Irrational));
            let ir = fq(Irrational, e, atoms);
            out.push(implies(ar, implies(oi, ir)));
            (out, args)
        }
        Kind::Mul(..) => {
            let args = crate::traverse::ordered_args(e);
            let z = fq(Zero, e, atoms);
            let az = anyarg(&args, atoms, &mut p(Zero));
            out.push(equiv(z, az));
            for pr in [Positive, Real, Rational, Integer] {
                let l = allargs(&args, atoms, &mut p(pr));
                let r = fq(pr, e, atoms);
                out.push(implies(l, r));
            }
            let l = exactlyonearg(&args, atoms, &mut np(Rational));
            let r = fq(Integer, e, atoms).negate();
            out.push(implies(l, r));
            let l = allargs(&args, atoms, &mut p(Commutative));
            let r = fq(Commutative, e, atoms);
            out.push(implies(l, r));
            let l = allargs(&args, atoms, &mut p(Prime));
            let r = fq(Prime, e, atoms).negate();
            out.push(implies(l, r));
            let l = allargs(&args, atoms, &mut |x: &Expr, atoms: &mut Atoms| {
                F::Or(vec![fq(Imaginary, x, atoms), fq(Real, x, atoms)])
            });
            let oi = exactlyonearg(&args, atoms, &mut p(Imaginary));
            let r = fq(Imaginary, e, atoms);
            out.push(implies(l, implies(oi, r)));
            let ar = allargs(&args, atoms, &mut p(Real));
            let oi = exactlyonearg(&args, atoms, &mut p(Irrational));
            let ir = fq(Irrational, e, atoms);
            out.push(implies(ar, implies(oi, ir)));
            let ai = allargs(&args, atoms, &mut p(Integer));
            let ae = anyarg(&args, atoms, &mut p(Even));
            let ev = fq(Even, e, atoms);
            out.push(implies(ai, equiv(ae, ev)));
            (out, args)
        }
        Kind::Pow(b, x) => {
            let conj = |atoms: &mut Atoms, pb: Pred, px: Pred| {
                F::And(vec![
                    fq(pb, b, atoms),
                    fq(px, x, atoms),
                    fq(Nonnegative, x, atoms),
                ])
            };
            let l = conj(atoms, Real, Even);
            out.push(implies(l, fq(Nonnegative, e, atoms)));
            let l = conj(atoms, Nonnegative, Odd);
            out.push(implies(l, fq(Nonnegative, e, atoms)));
            let l = conj(atoms, Nonpositive, Odd);
            out.push(implies(l, fq(Nonpositive, e, atoms)));
            let z = fq(Zero, e, atoms);
            let r = F::And(vec![fq(Zero, b, atoms), fq(Positive, x, atoms)]);
            out.push(equiv(z, r));
            (out, vec![b.clone(), x.clone()])
        }
        Kind::Func(Func::MatrixProxy(m), _) => {
            if let crate::matrixexpr::MatrixKind::Mul(coeff, factors) = m.kind() {
                let args = handlers::matmul_args(coeff, factors);
                let sq = allargs(&args, atoms, &mut p(Square));
                let inv = fq(Invertible, e, atoms);
                let ai = allargs(&args, atoms, &mut p(Invertible));
                out.push(implies(sq, equiv(inv, ai)));
                return (out, args);
            }
            (out, Vec::new())
        }
        _ => (out, Vec::new()),
    }
}

#[allow(dead_code)]
pub(super) fn logic_args(e: &Expr) -> Option<(LogicOp, &[Expr])> {
    crate::logic::as_logic(e)
}
