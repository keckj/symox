//! Formal power series: a port of SymPy 1.14 `sympy/series/formal.py`.
//!
//! The algorithms (`rational_algorithm`, `rational_independent`, `simpleDE`,
//! `exp_re`, `hyper_re`, `rsolve_hypergeometric`, `solve_de`,
//! `hyper_algorithm`, `_compute_fps`) and the term arithmetic of
//! `FormalPowerSeries` and its finite subclasses follow the reference step by
//! step. A series is `(function, x, x0, dir, ak, xk, ind)` where `ak` and
//! `xk` are formula sequences [`Seq`].
//!
//! The public operations the reference calls on intermediate expressions
//! (`simplify`, `limit`, `integrate`, `apart`, `factor`, `cancel`, `doit`,
//! `leadterm`, `roots`, `linsolve` and the `Coeff` class) go through
//! [`Hooks`], so that the host's public entry points (with their own
//! admission and fallbacks) are the ones used, as in the reference. The
//! rest (`rsolve`, differentiation, substitution, `collect`, `expand`,
//! `powsimp`, `Order`, `Piecewise`, Bell polynomials, comparisons) is native.
//!
//! Errors raised by the reference on the same inputs are returned as
//! [`Exc`] with the reference exception class.

use crate::decompose::{as_base_exp, as_coeff_add, as_coeff_add_head, as_coeff_mul, as_independent};
use crate::expr::{Expr, Func, Kind, RelOp};
use crate::functions::func;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sub};
use crate::subs::subs;
use crate::symbol::Assumptions;
use crate::traverse::{free_symbols, has_structural, ordered_args};

/// A reference exception.
#[derive(Clone, Debug)]
pub enum Exc {
    Value(String),
    NotImplemented(String),
    Type(String),
    Index(String),
    Key(String),
}

/// A reference exception, or an error of a hook (propagated unchanged).
#[derive(Debug)]
pub enum Error<E> {
    Exc(Exc),
    Hook(E),
}

impl<E> From<Exc> for Error<E> {
    fn from(e: Exc) -> Self {
        Error::Exc(e)
    }
}

pub type R<T, E> = Result<T, Error<E>>;

/// The public operations of the host used by the reference algorithms.
pub trait Hooks {
    type Error;
    /// `e.simplify()`
    fn simplify(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    /// `e.limit(x, 0)`; a held limit is returned as a `Limit` node.
    fn limit0(&mut self, e: &Expr, x: &Expr) -> Result<Expr, Self::Error>;
    /// `integrate(e, x)`
    fn integrate(&mut self, e: &Expr, x: &Expr) -> Result<Expr, Self::Error>;
    /// `apart(e, x, full=full)`
    fn apart(&mut self, e: &Expr, x: &Expr, full: bool) -> Result<Expr, Self::Error>;
    /// `e.doit()`
    fn doit(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    /// `e.factor()`
    fn factor(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    /// `e.cancel()`
    fn cancel(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    /// `e.leadterm(x)[0]`
    fn leadterm_coeff(&mut self, e: &Expr, x: &Expr) -> Result<Expr, Self::Error>;
    /// `roots(e, x)` as `(root, multiplicity)` pairs.
    fn roots(&mut self, e: &Expr, x: &Expr) -> Result<Vec<(Expr, Expr)>, Self::Error>;
    /// `[i for s in linsolve(eqs, syms) for i in s]`
    fn linsolve(&mut self, eqs: &[Expr], syms: &[Expr]) -> Result<Vec<Expr>, Self::Error>;
    /// `Coeff(p, x, n)`
    fn coeff_node(&mut self, p: &Expr, x: &Expr, n: &Expr) -> Result<Expr, Self::Error>;
    /// Called in unbounded loops (the reference iterates without bound).
    fn tick(&mut self) -> Result<(), Self::Error> {
        Ok(())
    }
}

// ---------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------

fn named(name: &str, args: Vec<Expr>) -> Expr {
    func(Func::from_name(name).unwrap_or_else(|| panic!("unknown function {name}")), args)
}

fn factorial(e: &Expr) -> Expr {
    named("factorial", vec![e.clone()])
}

fn rf(a: &Expr, n: &Expr) -> Expr {
    named("RisingFactorial", vec![a.clone(), n.clone()])
}

fn int(i: i64) -> Expr {
    Expr::int(i)
}

fn dummy_k() -> Expr {
    Expr::dummy("k", Assumptions::NONE)
}

/// `Add.make_args(e)` in the reference argument order.
fn make_add(e: &Expr) -> Vec<Expr> {
    if e.is_add() {
        ordered_args(e)
    } else {
        vec![e.clone()]
    }
}

/// Python truthiness of a value: false only for a zero number or `false`.
fn truthy(e: &Expr) -> bool {
    !(e.is_zero() || matches!(e.kind(), Kind::Bool(false)))
}

/// `e is S.Zero`
fn is_exact_zero(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if n.is_zero() && n.is_integer())
}

fn has(e: &Expr, s: &Expr) -> bool {
    has_structural(e, s)
}

fn has_any(e: &Expr, syms: &[Expr]) -> bool {
    syms.iter().any(|s| has(e, s))
}

/// Children of a node, including the parts of derivatives and piecewise.
fn children(e: &Expr) -> Vec<Expr> {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => vec![],
        Kind::Derivative(inner, vars) => {
            let mut v = vec![inner.clone()];
            v.extend(vars.iter().map(|(s, _)| s.clone()));
            v
        }
        Kind::Piecewise(pairs) => pairs.iter().flat_map(|(a, c)| [a.clone(), c.clone()]).collect(),
        _ => e.args(),
    }
}

fn walk(e: &Expr, f: &mut dyn FnMut(&Expr) -> bool) -> bool {
    if f(e) {
        return true;
    }
    for c in children(e) {
        if walk(&c, f) {
            return true;
        }
    }
    false
}

/// `e.has(Symbol)`
fn has_symbol_atom(e: &Expr) -> bool {
    walk(e, &mut |n| n.is_symbol())
}

/// `e.has(Derivative)`
fn has_derivative(e: &Expr) -> bool {
    walk(e, &mut |n| matches!(n.kind(), Kind::Derivative(..)) || n.is_func(&Func::SymbolicDerivative))
}

/// `e.has(g)` for an undefined function head `g`.
fn has_head(e: &Expr, head: &Func) -> bool {
    walk(e, &mut |n| matches!(n.kind(), Kind::Func(h, _) if h == head))
}

fn is_tuple(e: &Expr) -> bool {
    crate::special::tuple_items(e).is_some()
}

/// Is `e` an instance of the reference `Function` class?
fn is_function_app(e: &Expr) -> bool {
    match e.kind() {
        Kind::Piecewise(_) => true,
        Kind::Func(h, _) => {
            !is_tuple(e)
                && !matches!(
                    h,
                    Func::Sum
                        | Func::Product
                        | Func::Integral
                        | Func::Limit
                        | Func::Subs
                        | Func::Order
                        | Func::SymbolicDerivative
                        | Func::Logic(_)
                        | Func::Set(_)
                )
        }
        _ => false,
    }
}

/// `e.atoms(Function)`, in first-seen order.
fn function_atoms(e: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    walk(e, &mut |n| {
        if is_function_app(n) && !out.contains(n) {
            out.push(n.clone());
        }
        false
    });
    out
}

/// `e.atoms(Function).pop()`: an undefined function application when
/// there is one (the reference uses this on differential equations in an
/// undefined function), else any function application.
fn pop_function_atom(e: &Expr) -> Result<Expr, Exc> {
    let atoms = function_atoms(e);
    if let Some(a) = atoms.iter().find(|a| matches!(a.kind(), Kind::Func(Func::Undefined(_), _))) {
        return Ok(a.clone());
    }
    atoms.into_iter().next().ok_or_else(|| Exc::Key("'pop from an empty set'".into()))
}

fn apply(head: &Func, arg: &Expr) -> Expr {
    func(head.clone(), vec![arg.clone()])
}

/// `Derivative(e, x, n)` (`e` itself for `n = 0`).
fn derivative(e: &Expr, x: &Expr, n: i64) -> Expr {
    if n == 0 {
        e.clone()
    } else {
        crate::diff::diff_n(e, x, n as u32)
    }
}

fn derivative_count(e: &Expr) -> Option<i64> {
    match e.kind() {
        Kind::Derivative(_, vars) => Some(vars.iter().map(|(_, n)| *n as i64).sum()),
        _ => None,
    }
}

/// Truth value of a relational, `TypeError` when it is undecided.
fn truth(op: RelOp, a: &Expr, b: &Expr) -> Result<bool, Exc> {
    let r = crate::relational::rel(op, a, b);
    match r.kind() {
        Kind::Bool(v) => Ok(*v),
        _ => Err(Exc::Type(format!("cannot determine truth value of Relational: {r}"))),
    }
}

/// A Python integer from an expression (`range` arguments).
fn index(e: &Expr) -> Result<i64, Exc> {
    if e.is_integer_number() {
        if let Some(i) = e.as_i64() {
            return Ok(i);
        }
    }
    Err(Exc::Type(format!("'{e}' object cannot be interpreted as an integer")))
}

fn sorted_syms(mut v: Vec<Expr>) -> Vec<Expr> {
    v.sort_by(|a, b| a.to_string().cmp(&b.to_string()).then_with(|| a.cmp_struct(b)));
    v
}

fn free_minus(e: &Expr, drop: &[Expr]) -> Vec<Expr> {
    free_symbols(e).into_iter().filter(|s| !drop.contains(s)).collect()
}

/// Sequential substitution of a solution dict (`expr.subs(dict)`): the
/// reference orders symbols of equal size by their sort key (the name).
fn subs_solution(e: &Expr, syms: &[Expr], values: &[Expr]) -> Expr {
    let mut pairs: Vec<(Expr, Expr)> = syms.iter().cloned().zip(values.iter().cloned()).collect();
    pairs.sort_by(|a, b| a.0.to_string().cmp(&b.0.to_string()).then_with(|| a.0.cmp_struct(&b.0)));
    subs(e, &pairs)
}

fn collect(e: &Expr, pattern: &Expr) -> Expr {
    crate::poly::api::collect(e, std::slice::from_ref(pattern))
}

fn expand(e: &Expr) -> Expr {
    crate::expand::expand(e)
}

fn powsimp(e: &Expr) -> Expr {
    crate::simptools::powsimp::powsimp(e, &crate::simptools::powsimp::PowsimpOptions::default())
}

fn coeff1(e: &Expr, x: &Expr) -> Expr {
    crate::simptools::coeff(e, x, 1)
}

fn is_rational_function(e: &Expr, x: &Expr) -> bool {
    crate::inspection::classify(e, std::slice::from_ref(x), crate::inspection::PolynomialKind::Rational) == Some(true)
}

fn is_polynomial(e: &Expr, x: &Expr) -> bool {
    crate::inspection::classify(e, std::slice::from_ref(x), crate::inspection::PolynomialKind::Polynomial) == Some(true)
}

/// `e.as_coeff_exponent(x)`
fn as_coeff_exponent(e: &Expr, x: &Expr) -> (Expr, Expr) {
    let s = collect(e, x);
    let (c, p) = as_coeff_mul(&s, std::slice::from_ref(x), true);
    if p.len() == 1 {
        let (b, ex) = as_base_exp(&p[0]);
        if &b == x {
            return (c, ex);
        }
    }
    (s, Expr::zero())
}

/// `e.as_coeff_mul(*deps)[1][0]` where `has` selects the dependent factors.
fn first_dependent_factor(e: &Expr, has: &dyn Fn(&Expr) -> bool) -> Result<Expr, Exc> {
    let factors = if e.is_mul() { ordered_args(e) } else { vec![e.clone()] };
    factors
        .into_iter()
        .find(|f| has(f))
        .ok_or_else(|| Exc::Index("tuple index out of range".into()))
}

fn is_limit(e: &Expr) -> bool {
    e.is_func(&Func::Limit)
}

fn lift<T, E>(r: Result<T, E>) -> R<T, E> {
    r.map_err(Error::Hook)
}

// ---------------------------------------------------------------------
// sequences
// ---------------------------------------------------------------------

/// A `SeqFormula(formula, (k, start, stop))`.
#[derive(Clone, Debug)]
pub struct Seq {
    pub formula: Expr,
    pub k: Expr,
    pub start: Expr,
    pub stop: Expr,
}

impl Seq {
    pub fn new(formula: Expr, k: Expr, start: Expr) -> Seq {
        Seq { formula, k, start, stop: Expr::oo() }
    }

    /// `seq.coeff(pt)`
    pub fn coeff(&self, pt: &Expr) -> Result<Expr, Exc> {
        if truth(RelOp::Lt, pt, &self.start)? || truth(RelOp::Gt, pt, &self.stop)? {
            return Err(Exc::Index(format!("Index {pt} out of bounds Interval({}, {})", self.start, self.stop)));
        }
        Ok(subs(&self.formula, &[(self.k.clone(), pt.clone())]))
    }

    /// `seq[i]` (the `i`-th point from the start)
    pub fn ith(&self, i: i64) -> Result<Expr, Exc> {
        self.coeff(&add2(&self.start, &int(i)))
    }

    /// `self + other` (`SeqFormula._add`)
    pub fn add(&self, other: &Seq) -> Seq {
        let formula = add2(&self.formula, &subs(&other.formula, &[(other.k.clone(), self.k.clone())]));
        Seq {
            formula,
            k: self.k.clone(),
            start: crate::functions::max(vec![self.start.clone(), other.start.clone()]),
            stop: crate::functions::min(vec![self.stop.clone(), other.stop.clone()]),
        }
    }

    /// `seq.coeff_mul(c)`
    pub fn coeff_mul(&self, c: &Expr) -> Seq {
        Seq { formula: mul2(&self.formula, c), ..self.clone() }
    }

    pub fn subs(&self, old: &Expr, new: &Expr) -> Seq {
        let p = [(old.clone(), new.clone())];
        Seq { formula: subs(&self.formula, &p), k: self.k.clone(), start: subs(&self.start, &p), stop: subs(&self.stop, &p) }
    }

    fn free_symbols(&self) -> Vec<Expr> {
        let mut out = free_minus(&self.formula, std::slice::from_ref(&self.k));
        for e in [&self.start, &self.stop] {
            for s in free_symbols(e) {
                if !out.contains(&s) && s != self.k {
                    out.push(s);
                }
            }
        }
        out
    }

    /// The held `Sum(formula, (k, start, stop))` encoding the sequence.
    pub fn to_sum(&self) -> Result<Expr, Exc> {
        crate::binders::held(Func::Sum, &self.formula, vec![vec![self.k.clone(), self.start.clone(), self.stop.clone()]])
            .map_err(Exc::Value)
    }

    pub fn from_sum(e: &Expr) -> Result<Seq, Exc> {
        if let Kind::Func(Func::Sum, args) = e.kind() {
            if args.len() == 2 {
                let lim = args[1].args();
                if lim.len() == 3 {
                    return Ok(Seq { formula: args[0].clone(), k: lim[0].clone(), start: lim[1].clone(), stop: lim[2].clone() });
                }
            }
        }
        Err(Exc::Type(format!("{e} does not encode a sequence")))
    }
}

/// `SeqFormula` products: the coefficient of `ffps.bell_coeff_seq` at
/// `pt`, i.e. `(ak.formula*factorial(k)).subs(k, pt)`.
fn bell_coeff(ak: &Seq, pt: i64) -> Expr {
    let formula = mul2(&ak.formula, &factorial(&ak.k));
    subs(&formula, &[(ak.k.clone(), int(pt))])
}

// ---------------------------------------------------------------------
// rational_algorithm, rational_independent
// ---------------------------------------------------------------------

/// `binomial(n, k).rewrite(factorial)`
fn binomial_factorial(n: &Expr, k: &Expr) -> Expr {
    let b = named("binomial", vec![n.clone(), k.clone()]);
    if matches!(b.kind(), Kind::Func(Func::Special(crate::special::SpecialFn::Binomial), _)) {
        div(&factorial(n), &mul2(&factorial(k), &factorial(&sub(n, k))))
    } else {
        b
    }
}

/// `rational_algorithm(f, x, k, order, full)`: `(coeff, sep, start)`.
pub fn rational_algorithm<H: Hooks>(
    h: &mut H,
    f: &Expr,
    x: &Expr,
    k: &Expr,
    order: i64,
    full: bool,
) -> R<Option<(Expr, Expr, Expr)>, H::Error> {
    let mut diff = f.clone();
    let mut ds: Vec<Expr> = Vec::new();
    for i in 0..=order {
        if i > 0 {
            diff = crate::diff::diff(&diff, x);
        }
        if !is_rational_function(&diff, x) {
            ds.push(diff.clone());
            continue;
        }
        let (mut coeff, mut sep) = (Expr::zero(), Expr::zero());
        let mut terms = lift(h.apart(&diff, x, full))?;
        if walk(&terms, &mut |n| matches!(n.kind(), Kind::Func(Func::Undefined(u), _) if u.name.as_ref() == "RootSum")) {
            terms = lift(h.doit(&terms))?;
        }
        for t in make_add(&terms) {
            let (mut num, mut den) = crate::simptools::as_numer_denom(&t);
            if !has(&den, x) {
                sep = add2(&sep, &t);
                continue;
            }
            if den.is_mul() {
                let (i0, i1) = as_independent(&den, std::slice::from_ref(x), None);
                den = i1;
                num = div(&num, &i0);
            }
            let (den, j) = as_base_exp(&den);
            let (mut a, xterm) = as_coeff_add(&den, std::slice::from_ref(x));
            if !truthy(&a) {
                sep = add2(&sep, &t);
                continue;
            }
            let xc = coeff1(&xterm[0], x);
            a = div(&a, &neg(&xc));
            num = div(&num, &pow(&xc, &j));
            let mut ak = mul2(&pow(&int(-1), &j), &num);
            ak = mul2(&ak, &binomial_factorial(&add(&[j.clone(), k.clone(), int(-1)]), k));
            ak = div(&ak, &pow(&a, &add2(&j, k)));
            coeff = add2(&coeff, &ak);
        }
        if crate::assumptions::is_zero(&coeff) == Some(true) {
            return Ok(None);
        }
        if has(&coeff, x) || has(&coeff, &Expr::zoo()) || has(&coeff, &Expr::oo()) || has(&coeff, &Expr::nan()) {
            return Ok(None);
        }
        for j in 0..i {
            coeff = div(&coeff, &add2(k, &int(j + 1)));
            sep = lift(h.integrate(&sep, x))?;
            let d = ds.pop().ok_or_else(|| Exc::Index("pop from empty list".into()))?;
            sep = add2(&sep, &lift(h.limit0(&sub(&d, &sep), x))?);
        }
        return Ok(Some((subs(&coeff, &[(k.clone(), sub(k, &int(i)))]), sep, int(i))));
    }
    Ok(None)
}

/// `rational_independent(terms, x)`
pub fn rational_independent<H: Hooks>(h: &mut H, terms: &[Expr], x: &Expr) -> R<Vec<Expr>, H::Error> {
    if terms.is_empty() {
        return Ok(vec![]);
    }
    let mut ind = vec![terms[0].clone()];
    for t in &terms[1..] {
        let n = as_independent(t, std::slice::from_ref(x), None).1;
        let mut placed = false;
        for i in 0..ind.len() {
            let d = as_independent(&ind[i], std::slice::from_ref(x), None).1;
            let q = lift(h.cancel(&div(&n, &d)))?;
            if is_rational_function(&q, x) {
                ind[i] = add2(&ind[i], t);
                placed = true;
                break;
            }
        }
        if !placed {
            ind.push(t.clone());
        }
    }
    Ok(ind)
}

// ---------------------------------------------------------------------
// simpleDE, exp_re, hyper_re
// ---------------------------------------------------------------------

/// The lazy generator `simpleDE(f, x, g, order)`.
pub struct SimpleDe {
    f: Expr,
    x: Expr,
    g: Func,
    order: i64,
    a: Vec<Expr>,
    found: bool,
    k: i64,
}

impl SimpleDe {
    pub fn new(f: &Expr, x: &Expr, g: &Func, order: i64) -> SimpleDe {
        let a = (0..order.max(0)).map(|i| Expr::symbol(&format!("a{i}"))).collect();
        SimpleDe { f: f.clone(), x: x.clone(), g: g.clone(), order, a, found: false, k: 0 }
    }

    /// The next `(DE, k)`, or `None` when the generator is exhausted.
    pub fn next<H: Hooks>(&mut self, h: &mut H) -> R<Option<(Expr, i64)>, H::Error> {
        let x = self.x.clone();
        let gx = apply(&self.g, &x);
        while self.k < self.order {
            self.k += 1;
            let k = self.k;
            let mut eq = derivative(&self.f, &x, k);
            let mut de = derivative(&gx, &x, k);
            let eq_terms: Vec<Expr> = (0..k).map(|i| mul2(&self.a[i as usize], &derivative(&self.f, &x, i))).collect();
            let de_terms: Vec<Expr> = (0..k).map(|i| mul2(&self.a[i as usize], &derivative(&gx, &x, i))).collect();
            eq = add2(&eq, &add(&eq_terms));
            de = add2(&de, &add(&de_terms));
            let eq = expand(&eq);
            let terms = if eq.is_add() { crate::sortkey::ordered_terms(&eq) } else { vec![eq.clone()] };
            let ind = rational_independent(h, &terms, &x)?;
            if self.found || ind.len() as i64 == k {
                let flat = lift(h.linsolve(&ind, &self.a[..k as usize]))?;
                let n = flat.len().min(self.a.len());
                if n > 0 {
                    self.found = true;
                    de = subs_solution(&de, &self.a[..n], &flat[..n]);
                }
                de = crate::simptools::as_numer_denom(&de).0;
                de = lift(h.factor(&de))?;
                de = first_dependent_factor(&de, &has_derivative)?;
                return Ok(Some((collect(&de, &derivative(&gx, &x, 1)), k)));
            }
        }
        Ok(None)
    }
}

/// `exp_re(DE, r, k)`; `r` is the head of the sequence function.
pub fn exp_re(de: &Expr, r: &Func, k: &Expr) -> Result<Expr, Exc> {
    let mut re = Expr::zero();
    let g = pop_function_atom(de)?;
    let mut mini: Option<i64> = None;
    for t in make_add(de) {
        let (coeff, d) = as_independent(&t, std::slice::from_ref(&g), None);
        let j = derivative_count(&d).unwrap_or(0);
        if mini.is_none_or(|m| j < m) {
            mini = Some(j);
        }
        re = add2(&re, &mul2(&coeff, &apply(r, &add2(k, &int(j)))));
    }
    if let Some(m) = mini.filter(|m| *m != 0) {
        re = subs(&re, &[(k.clone(), sub(k, &int(m)))]);
    }
    Ok(re)
}

/// `hyper_re(DE, r, k)`; `r` is the head of the sequence function.
pub fn hyper_re(de: &Expr, r: &Func, k: &Expr) -> Result<Expr, Exc> {
    let mut re = Expr::zero();
    let g = pop_function_atom(de)?;
    let x = walk_first_symbol(&g).ok_or_else(|| Exc::Key("'pop from an empty set'".into()))?;
    let mut mini: Option<Expr> = None;
    for t in make_add(&expand(de)) {
        let (coeff, d) = as_independent(&t, std::slice::from_ref(&g), None);
        let (c, v) = as_independent(&coeff, std::slice::from_ref(&x), None);
        let l = as_coeff_exponent(&v, &x).1;
        let j = int(derivative_count(&d).unwrap_or(0));
        let term = mul(&[c, rf(&sub(&add2(k, &int(1)), &l), &j), apply(r, &sub(&add2(k, &j), &l))]);
        re = add2(&re, &term);
        let jl = sub(&j, &l);
        let smaller = match &mini {
            None => true,
            Some(m) => truth(RelOp::Lt, &jl, m)?,
        };
        if smaller {
            mini = Some(jl);
        }
    }
    let mini = mini.ok_or_else(|| Exc::Type("unsupported operand type(s) for -: 'Dummy' and 'NoneType'".into()))?;
    re = subs(&re, &[(k.clone(), sub(k, &mini))]);
    let m = Expr::wild("m", vec![]);
    Ok(collect(&re, &apply(r, &add2(k, &m))))
}

fn walk_first_symbol(e: &Expr) -> Option<Expr> {
    let mut found = None;
    walk(e, &mut |n| {
        if n.is_symbol() {
            found = Some(n.clone());
            true
        } else {
            false
        }
    });
    found
}

// ---------------------------------------------------------------------
// rsolve_hypergeometric
// ---------------------------------------------------------------------

struct Re {
    f: Expr,
    p: Expr,
    q: Expr,
    m: Expr,
}

fn transformation_a(t: Re, x: &Expr, k: &Expr, shift: &Expr) -> Re {
    let ks = [(k.clone(), add2(k, shift))];
    Re { f: mul2(&t.f, &pow(x, &neg(shift))), p: subs(&t.p, &ks), q: subs(&t.q, &ks), m: t.m }
}

fn transformation_c(t: Re, x: &Expr, k: &Expr, scale: &Expr) -> Re {
    let ks = [(k.clone(), div(k, scale))];
    Re {
        f: subs(&t.f, &[(x.clone(), pow(x, scale))]),
        p: subs(&t.p, &ks),
        q: subs(&t.q, &ks),
        m: mul2(&t.m, scale),
    }
}

fn transformation_e(t: Re, x: &Expr, k: &Expr) -> Re {
    let ks = [(k.clone(), add2(k, &int(1)))];
    let p = mul2(&subs(&t.p, &ks), &add(&[k.clone(), t.m.clone(), int(1)]));
    let q = mul2(&subs(&t.q, &ks), &add2(k, &int(1)));
    Re { f: crate::diff::diff(&t.f, x), p, q, m: t.m }
}

type Sol = Vec<(Expr, Expr)>;

fn apply_shift(sol: Sol, shift: &Expr) -> Sol {
    sol.into_iter().map(|(r, c)| (r, add2(&c, shift))).collect()
}

fn apply_scale(sol: Sol, scale: &Expr) -> Sol {
    sol.into_iter().map(|(r, c)| (r, div(&c, scale))).collect()
}

fn apply_integrate(sol: Sol, k: &Expr) -> Sol {
    sol.into_iter()
        .map(|(r, c)| {
            let rest = as_coeff_add_head(&c, false).1;
            let d = mul2(&add2(&c, &int(1)), &coeff1(&rest, k));
            (div(&r, &d), add2(&c, &int(1)))
        })
        .collect()
}

fn roots_product<H: Hooks>(h: &mut H, p: &Expr, k: &Expr) -> R<Expr, H::Error> {
    let rs = lift(h.roots(p, k))?;
    Ok(mul(&rs.iter().map(|(r, m)| pow(&rf(&neg(r), k), m)).collect::<Vec<_>>()))
}

fn compute_formula<H: Hooks>(h: &mut H, t: &Re, x: &Expr, k: &Expr, k_max: &Expr) -> R<Sol, H::Error> {
    let mut sol = Vec::new();
    let lo = index(&add2(k_max, &int(1)))?;
    let hi = index(&add(&[k_max.clone(), t.m.clone(), int(1)]))?;
    for i in lo..hi {
        if i < 0 {
            continue;
        }
        let r = div(&lift(h.limit0(&derivative(&t.f, x, i), x))?, &factorial(&int(i)));
        if crate::assumptions::is_zero(&r) == Some(true) {
            continue;
        }
        let kterm = add2(&mul2(&t.m, k), &int(i));
        let mut res = r;
        let p = subs(&t.p, &[(k.clone(), kterm.clone())]);
        let q = subs(&t.q, &[(k.clone(), kterm.clone())]);
        let kinv = [(k.clone(), pow(k, &int(-1)))];
        let c1 = lift(h.leadterm_coeff(&subs(&p, &kinv), k))?;
        let c2 = lift(h.leadterm_coeff(&subs(&q, &kinv), k))?;
        res = mul2(&res, &pow(&div(&neg(&c1), &c2), k));
        res = mul2(&res, &roots_product(h, &p, k)?);
        res = div(&res, &roots_product(h, &q, k)?);
        sol.push((res, kterm));
    }
    Ok(sol)
}

fn lcm_denominators(roots: &[Expr]) -> Result<Expr, Exc> {
    let mut acc = num_bigint::BigInt::from(1);
    for r in roots {
        if crate::assumptions::is_rational(r) != Some(true) {
            continue;
        }
        let d = crate::simptools::as_numer_denom(r).1;
        match d.as_number().and_then(|n| n.to_bigint()) {
            Some(d) => {
                let g = num_integer::Integer::gcd(&acc, &d);
                acc = &acc / &g * &d;
                if acc < num_bigint::BigInt::from(0) {
                    acc = -acc;
                }
            }
            None => return Err(Exc::NotImplemented(format!("native fps: lcm of the symbolic denominator {d}"))),
        }
    }
    Ok(Expr::num(crate::number::Number::from_bigint(acc)))
}

fn root_keys<H: Hooks>(h: &mut H, p: &Expr, k: &Expr) -> R<Vec<Expr>, H::Error> {
    Ok(lift(h.roots(p, k))?.into_iter().map(|(r, _)| r).collect())
}

type RsolveResult = Option<(Sol, Expr, Expr)>;

fn rsolve_hypergeometric_raw<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, p: &Expr, q: &Expr, k: &Expr, m: &Expr) -> R<RsolveResult, H::Error> {
    let proots = root_keys(h, p, k)?;
    let qroots = root_keys(h, q, k)?;
    let mut all = proots;
    for r in qroots {
        if !all.contains(&r) {
            all.push(r);
        }
    }
    let scale = lcm_denominators(&all)?;
    let mut t = transformation_c(Re { f: f.clone(), p: p.clone(), q: q.clone(), m: m.clone() }, x, k, &scale);
    let qroots = root_keys(h, &t.q, k)?;
    let k_min = if qroots.is_empty() { Expr::zero() } else { crate::functions::min(qroots) };
    let shift = add2(&k_min, &t.m);
    t = transformation_a(t, x, k, &shift);
    let l = lift(h.limit0(&mul2(x, &t.f), x))?;
    if !is_limit(&l) && !l.is_zero() {
        return Ok(None);
    }
    let qroots = root_keys(h, &t.q, k)?;
    let k_max = if qroots.is_empty() { Expr::zero() } else { crate::functions::max(qroots) };
    let (mut ind, mut mp) = (Expr::zero(), Expr::neg_oo());
    let n = index(&add(&[k_max.clone(), t.m.clone(), int(1)]))?;
    for i in 0..n {
        let r = div(&lift(h.limit0(&derivative(&t.f, x, i), x))?, &factorial(&int(i)));
        if crate::assumptions::is_finite(&r) == Some(false) {
            let old_f = t.f.clone();
            t = transformation_a(t, x, k, &int(i));
            t = transformation_e(t, x, k);
            let Some((sol, ind, mp)) = rsolve_hypergeometric_raw(h, &t.f, x, &t.p, &t.q, k, &t.m)? else {
                return Err(Exc::Type("cannot unpack non-iterable NoneType object".into()).into());
            };
            let sol = apply_shift(apply_integrate(sol, k), &int(i));
            let mut ind = lift(h.integrate(&ind, x))?;
            ind = add2(&ind, &lift(h.limit0(&sub(&old_f, &ind), x))?);
            return Ok(Some((sol, ind, add2(&mp, &int(1)))));
        } else if truthy(&r) {
            let e = add2(&int(i), &shift);
            ind = add2(&ind, &mul2(&r, &pow(x, &e)));
            if !e.is_integer_number() {
                return Err(Exc::Type(format!("invalid input: {e}")).into());
            }
            let pow_x = div(&e, &scale);
            if truth(RelOp::Gt, &pow_x, &mp)? {
                mp = pow_x;
            }
        }
    }
    ind = subs(&ind, &[(x.clone(), pow(x, &div(&int(1), &scale)))]);
    let sol = compute_formula(h, &t, x, k, &k_max)?;
    let sol = apply_scale(apply_shift(sol, &shift), &scale);
    Ok(Some((sol, ind, mp)))
}

/// `rsolve_hypergeometric(f, x, P, Q, k, m)`: `(formula, ind, start)`.
pub fn rsolve_hypergeometric<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, p: &Expr, q: &Expr, k: &Expr, m: &Expr) -> R<Option<(Expr, Expr, Expr)>, H::Error> {
    let Some((sol_list, mut ind, mp)) = rsolve_hypergeometric_raw(h, f, x, p, q, k, m)? else {
        return Ok(None);
    };
    let mut sol_dict: Vec<(Expr, Expr)> = Vec::new();
    for (res, cond) in sol_list {
        let (mut j, mk) = as_coeff_add_head(&cond, false);
        let c = coeff1(&mk, k);
        let mut res = res;
        if crate::assumptions::is_integer(&j) == Some(false) {
            res = mul2(&res, &pow(x, &named("frac", vec![j.clone()])));
            j = named("floor", vec![j.clone()]);
        }
        res = subs(&res, &[(k.clone(), div(&sub(k, &j), &c))]);
        let cond = crate::relational::rel(RelOp::Eq, &named("Mod", vec![k.clone(), c.clone()]), &named("Mod", vec![j.clone(), c.clone()]));
        if let Some(entry) = sol_dict.iter_mut().find(|(cd, _)| *cd == cond) {
            entry.1 = add2(&entry.1, &res);
        } else {
            sol_dict.push((cond, add2(&Expr::zero(), &res)));
        }
    }
    let mut pieces: Vec<(Expr, Expr)> = sol_dict.into_iter().map(|(c, r)| (r, c)).collect();
    pieces.push((Expr::zero(), Expr::boolean(true)));
    let sol = crate::relational::piecewise(pieces);
    let mut s = if mp == Expr::neg_oo() {
        Expr::zero()
    } else if crate::assumptions::is_integer(&mp) == Some(false) {
        named("ceiling", vec![mp.clone()])
    } else {
        add2(&mp, &int(1))
    };
    if truth(RelOp::Lt, &s, &Expr::zero())? {
        let lo = index(&s)?;
        let term = mul2(&sol, &pow(x, k));
        let mut total = Expr::zero();
        for pt in lo..0 {
            total = add2(&total, &subs(&term, &[(k.clone(), int(pt))]));
        }
        ind = add2(&ind, &total);
        s = Expr::zero();
    }
    Ok(Some((sol, ind, s)))
}

// ---------------------------------------------------------------------
// solve_de, hyper_algorithm
// ---------------------------------------------------------------------

type Found = Option<(Expr, Expr, Expr)>;

fn solve_hyper_re<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, re: &Expr, k: &Expr) -> R<Found, H::Error> {
    if make_add(re).len() != 2 {
        return Ok(None);
    }
    let gs = function_atoms(re);
    if gs.len() != 2 {
        let msg = if gs.len() > 2 { "too many values to unpack (expected 2)".to_string() } else { format!("not enough values to unpack (expected 2, got {})", gs.len()) };
        return Err(Exc::Value(msg).into());
    }
    let (mut p, mut q) = (coeff1(re, &gs[0]), coeff1(re, &gs[1]));
    let arg0 = |g: &Expr| g.args().first().cloned().unwrap_or_else(Expr::zero);
    let mut m = sub(&arg0(&gs[1]), &arg0(&gs[0]));
    if truth(RelOp::Lt, &m, &Expr::zero())? {
        std::mem::swap(&mut p, &mut q);
        m = named("Abs", vec![m]);
    }
    rsolve_hypergeometric(h, f, x, &p, &q, k, &m)
}

fn rsolve_native(re: &Expr, gk: &Expr, init: &[(Expr, Expr)]) -> Result<Option<Expr>, Exc> {
    match crate::recurr::rsolve(re, gk, Some(init)) {
        Ok(v) => Ok(v),
        Err(crate::recurr::RsolveError::Value(m)) => Err(Exc::Value(m)),
        Err(crate::recurr::RsolveError::Type(m)) => Err(Exc::Type(m)),
    }
}

fn solve_explike_de<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, de: &Expr, g: &Func, k: &Expr) -> R<Found, H::Error> {
    let gx = apply(g, x);
    for t in make_add(de) {
        let (coeff, _) = as_independent(&t, std::slice::from_ref(&gx), None);
        if !free_symbols(&coeff).is_empty() {
            return Ok(None);
        }
    }
    let re = exp_re(de, g, k)?;
    let mut f = f.clone();
    let mut init = Vec::new();
    for i in 0..make_add(&re).len() {
        if i > 0 {
            f = crate::diff::diff(&f, x);
        }
        init.push((apply(g, &int(i as i64)), lift(h.limit0(&f, x))?));
    }
    match rsolve_native(&re, &apply(g, k), &init)? {
        Some(sol) if truthy(&sol) => Ok(Some((div(&sol, &factorial(k)), Expr::zero(), Expr::zero()))),
        _ => Ok(None),
    }
}

fn solve_simple<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, de: &Expr, g: &Func, k: &Expr) -> R<Found, H::Error> {
    let re = hyper_re(de, g, k)?;
    let mut f = f.clone();
    let mut init = Vec::new();
    for i in 0..make_add(&re).len() {
        if i > 0 {
            f = crate::diff::diff(&f, x);
        }
        let v = div(&lift(h.limit0(&f, x))?, &factorial(&int(i as i64)));
        init.push((apply(g, &int(i as i64)), v));
    }
    match rsolve_native(&re, &apply(g, k), &init)? {
        Some(sol) if truthy(&sol) => Ok(Some((sol, Expr::zero(), Expr::zero()))),
        _ => Ok(None),
    }
}

fn transform_explike_de<H: Hooks>(h: &mut H, de: &Expr, g: &Func, x: &Expr, order: i64, syms: &[Expr]) -> R<Expr, H::Error> {
    let gx = apply(g, x);
    let mut eq: Vec<Expr> = Vec::new();
    let highest = coeff1(de, &derivative(&gx, x, order));
    for i in 0..order {
        let c = coeff1(de, &derivative(&gx, x, i));
        let c = collect(&expand(&div(&c, &highest)), x);
        eq.extend(make_add(&c));
    }
    let mut temp = Vec::new();
    let mut broke = false;
    for e in &eq {
        if has(e, x) {
            broke = true;
            break;
        } else if has_symbol_atom(e) {
            temp.push(e.clone());
        }
    }
    if !broke {
        eq = temp;
    }
    let mut de = de.clone();
    if !eq.is_empty() {
        let flat = lift(h.linsolve(&eq, syms))?;
        let n = flat.len().min(syms.len());
        if n > 0 {
            de = subs_solution(&de, &syms[..n], &flat[..n]);
            de = lift(h.factor(&de))?;
            de = first_dependent_factor(&de, &has_derivative)?;
            de = collect(&de, &derivative(&gx, x, 1));
        }
    }
    Ok(de)
}

fn transform_de_re<H: Hooks>(h: &mut H, de: &Expr, g: &Func, k: &Expr, order: i64, syms: &[Expr]) -> R<Expr, H::Error> {
    let mut re = hyper_re(de, g, k)?;
    let eq: Vec<Expr> = (1..order).map(|i| coeff1(&re, &apply(g, &add2(k, &int(i))))).collect();
    let flat = lift(h.linsolve(&eq, syms))?;
    let n = flat.len().min(syms.len());
    if n > 0 {
        let m = Expr::wild("m", vec![]);
        re = subs_solution(&re, &syms[..n], &flat[..n]);
        re = lift(h.factor(&re))?;
        re = collect(&crate::simptools::as_numer_denom(&re).0, &apply(g, &add2(k, &m)));
        re = first_dependent_factor(&re, &|e: &Expr| has_head(e, g))?;
        for i in 0..order {
            if truthy(&coeff1(&re, &apply(g, &add2(k, &int(i))))) && i > 0 {
                re = subs(&re, &[(k.clone(), sub(k, &int(i)))]);
                break;
            }
        }
    }
    Ok(re)
}

/// `solve_de(f, x, DE, order, g, k)`
pub fn solve_de<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, de: &Expr, order: i64, g: &Func, k: &Expr) -> R<Found, H::Error> {
    let syms = sorted_syms(free_minus(de, std::slice::from_ref(x)));
    let re = if syms.is_empty() { hyper_re(de, g, k)? } else { transform_de_re(h, de, g, k, order, &syms)? };
    let mut sol = None;
    if free_minus(&re, std::slice::from_ref(k)).is_empty() {
        sol = solve_hyper_re(h, f, x, &re, k)?;
    }
    if sol.is_some() {
        return Ok(sol);
    }
    let mut de = de.clone();
    if !syms.is_empty() {
        de = transform_explike_de(h, &de, g, x, order, &syms)?;
    }
    if free_minus(&de, std::slice::from_ref(x)).is_empty() {
        sol = solve_explike_de(h, f, x, &de, g, k)?;
    }
    Ok(sol)
}

/// `hyper_algorithm(f, x, k, order)`
pub fn hyper_algorithm<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, k: &Expr, order: i64) -> R<Found, H::Error> {
    let g = Func::Undefined("g".into());
    let mut des = Vec::new();
    let mut gen = SimpleDe::new(f, x, &g, order);
    while let Some((de, i)) = gen.next(h)? {
        let sol = solve_de(h, f, x, &de, i, &g, k)?;
        if sol.is_some() {
            return Ok(sol);
        }
        if free_minus(&de, std::slice::from_ref(x)).is_empty() {
            des.push(de);
        }
    }
    for de in des {
        let sol = solve_simple(h, f, x, &de, &g, k)?;
        if sol.is_some() {
            return Ok(sol);
        }
    }
    Ok(None)
}

// ---------------------------------------------------------------------
// _compute_fps
// ---------------------------------------------------------------------

/// `(ak, xk, ind)`
pub type Parts = (Seq, Seq, Expr);

/// Options of `fps`/`compute_fps`.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub hyper: bool,
    pub order: i64,
    pub rational: bool,
    pub full: bool,
}

fn compute_fps_raw<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, x0: &Expr, dir: &Expr, o: Options) -> R<Option<Parts>, H::Error> {
    h_tick(h)?;
    if *x0 == Expr::oo() || *x0 == Expr::neg_oo() {
        let dir = if *x0 == Expr::oo() { int(1) } else { int(-1) };
        let xinv = pow(x, &int(-1));
        let temp = subs(f, &[(x.clone(), xinv.clone())]);
        let Some((ak, xk, ind)) = compute_fps_raw(h, &temp, x, &Expr::zero(), &dir, o)? else {
            return Ok(None);
        };
        return Ok(Some((ak, xk.subs(x, &xinv), subs(&ind, &[(x.clone(), xinv)]))));
    }
    if truthy(x0) || *dir == int(-1) {
        let (rep, back) = if *dir == int(-1) { (add2(&neg(x), x0), add2(&neg(x), x0)) } else { (add2(x, x0), sub(x, x0)) };
        let temp = subs(f, &[(x.clone(), rep)]);
        let Some((ak, xk, ind)) = compute_fps_raw(h, &temp, x, &Expr::zero(), &int(1), o)? else {
            return Ok(None);
        };
        return Ok(Some((ak, xk.subs(x, &back), subs(&ind, &[(x.clone(), back)]))));
    }
    if is_polynomial(f, x) {
        let k = dummy_k();
        let c = lift(h.coeff_node(f, x, &k))?;
        let ak = Seq::new(c, k.clone(), int(1));
        let xk = Seq::new(pow(x, &k), k, int(0));
        // Deliberate deviation: the reference reads the constant term off
        // the unexpanded polynomial, so fps((1 + x)**2) lost its constant 1.
        let ind = crate::simptools::coeff(&expand(f), x, 0);
        return Ok(Some((ak, xk, ind)));
    }
    if f.is_add() {
        let mut result = false;
        let mut ak = Seq::new(Expr::zero(), dummy_k(), int(0));
        let mut ind = Expr::zero();
        let mut xk: Option<Seq> = None;
        for t in make_add(f) {
            match compute_fps_raw(h, &t, x, &Expr::zero(), &int(1), o)? {
                Some((rak, rxk, rind)) => {
                    if !result {
                        result = true;
                        xk = Some(rxk);
                    }
                    let xs = xk.as_ref().expect("set above");
                    let (seq, s, e) = if truth(RelOp::Gt, &rak.start, &ak.start)? {
                        (ak.clone(), ak.start.clone(), rak.start.clone())
                    } else {
                        (rak.clone(), rak.start.clone(), ak.start.clone())
                    };
                    let save = save_terms(&seq, xs, &s, &e)?;
                    ak = ak.add(&rak);
                    ind = add2(&ind, &add2(&rind, &save));
                }
                None => ind = add2(&ind, &t),
            }
        }
        return Ok(if result { Some((ak, xk.expect("set with result"), ind)) } else { None });
    }
    let syms = free_minus(f, std::slice::from_ref(x));
    let (f, symb) = as_independent(&expand(f), &syms, None);
    let k = dummy_k();
    let mut result = None;
    if o.rational {
        result = rational_algorithm(h, &f, x, &k, o.order, o.full)?;
    }
    if result.is_none() && o.hyper {
        result = hyper_algorithm(h, &f, x, &k, o.order)?;
    }
    let Some((formula, rind, start)) = result else {
        return Ok(None);
    };
    let symb = if crate::assumptions::is_zero(&symb) == Some(true) { int(1) } else { powsimp(&symb) };
    let ak = Seq::new(formula, k.clone(), start);
    let xk = Seq::new(powsimp(&mul2(&pow(x, &k), &symb)), k, int(0));
    let ind = powsimp(&mul2(&rind, &symb));
    Ok(Some((ak, xk, ind)))
}

fn h_tick<H: Hooks>(h: &mut H) -> R<(), H::Error> {
    lift(h.tick())
}

/// `Add(*[z[0]*z[1] for z in zip(seq[0:e - s], xk[s:e])])`
fn save_terms(seq: &Seq, xk: &Seq, s: &Expr, e: &Expr) -> Result<Expr, Exc> {
    let n = index(&sub(e, s))?;
    let s = index(s)?;
    let mut terms = Vec::new();
    for i in 0..n.max(0) {
        terms.push(mul2(&seq.ith(i)?, &xk.ith(s + i)?));
    }
    Ok(add(&terms))
}

/// `compute_fps(f, x, x0, dir, hyper, order, rational, full)`; `dir` is
/// already `1` or `-1`.
pub fn compute_fps<H: Hooks>(h: &mut H, f: &Expr, x: &Expr, x0: &Expr, dir: &Expr, o: Options) -> R<Option<Parts>, H::Error> {
    if !has(f, x) {
        return Ok(None);
    }
    if *dir != int(1) && *dir != int(-1) {
        return Err(Exc::Value("Dir must be '+' or '-'".into()).into());
    }
    compute_fps_raw(h, f, x, x0, dir, o)
}

// ---------------------------------------------------------------------
// FormalPowerSeries
// ---------------------------------------------------------------------

/// `FormalPowerSeries(function, x, x0, dir, (ak, xk, ind))`
#[derive(Clone, Debug)]
pub struct Fps {
    pub function: Expr,
    pub x: Expr,
    pub x0: Expr,
    pub dir: Expr,
    pub ak: Seq,
    pub xk: Seq,
    pub ind: Expr,
}

impl Fps {
    /// Decode the host node `(function, x, x0, dir, (Sum(ak), Sum(xk), ind))`.
    pub fn from_args(args: &[Expr]) -> Result<Fps, Exc> {
        if args.len() < 5 {
            return Err(Exc::Index("tuple index out of range".into()));
        }
        let parts = crate::special::tuple_items(&args[4]).ok_or_else(|| Exc::Type("malformed FormalPowerSeries".into()))?;
        if parts.len() != 3 {
            return Err(Exc::Type("malformed FormalPowerSeries".into()));
        }
        Ok(Fps {
            function: args[0].clone(),
            x: args[1].clone(),
            x0: args[2].clone(),
            dir: args[3].clone(),
            ak: Seq::from_sum(&parts[0])?,
            xk: Seq::from_sum(&parts[1])?,
            ind: parts[2].clone(),
        })
    }

    pub fn from_node(e: &Expr) -> Result<Fps, Exc> {
        Fps::from_args(&e.args())
    }

    /// The fifth argument `(Sum(ak), Sum(xk), ind)`.
    pub fn parts_tuple(&self) -> Result<Expr, Exc> {
        Ok(crate::special::tuple(vec![self.ak.to_sum()?, self.xk.to_sum()?, self.ind.clone()]))
    }

    fn with(&self, function: Expr, ak: Seq, xk: Seq, ind: Expr) -> Fps {
        Fps { function, x: self.x.clone(), x0: self.x0.clone(), dir: self.dir.clone(), ak, xk, ind }
    }

    /// `self.free_symbols`
    pub fn free_symbols(&self) -> Vec<Expr> {
        let mut out: Vec<Expr> = Vec::new();
        let mut put = |v: Vec<Expr>| {
            for s in v {
                if !out.contains(&s) {
                    out.push(s);
                }
            }
        };
        for e in [&self.function, &self.x, &self.x0, &self.dir] {
            put(free_symbols(e));
        }
        put(self.ak.free_symbols());
        put(self.xk.free_symbols());
        put(free_symbols(&self.ind));
        out
    }

    /// `_get_pow_x(term)`
    fn pow_x(&self, term: &Expr) -> Expr {
        let dep = as_independent(term, std::slice::from_ref(&self.x), None).1;
        let (xterm, pow_x) = as_base_exp(&dep);
        if !has(&xterm, &self.x) {
            return Expr::zero();
        }
        pow_x
    }

    fn pow_x_numeric(&self, term: &Expr, sym: &[Expr]) -> Expr {
        let p = self.pow_x(term);
        if has_any(&p, sym) {
            as_coeff_add(&p, sym).0
        } else {
            p
        }
    }

    /// `_eval_term(pt)`
    pub fn term<H: Hooks>(&self, h: &mut H, pt: &Expr) -> R<Expr, H::Error> {
        let mut term = match self.xk.coeff(pt) {
            Err(Exc::Index(_)) => Expr::zero(),
            Err(e) => return Err(e.into()),
            Ok(pt_xk) => match self.ak.coeff(pt) {
                Err(Exc::Index(_)) => Expr::zero(),
                Err(e) => return Err(e.into()),
                Ok(c) => mul2(&lift(h.simplify(&c))?, &pt_xk),
            },
        };
        if truthy(&self.ind) {
            let mut ind = Expr::zero();
            let sym = self.free_symbols();
            for t in make_add(&self.ind) {
                let p = self.pow_x_numeric(&t, &sym);
                if pt.is_zero() && truth(RelOp::Lt, &p, &int(1))? {
                    ind = add2(&ind, &t);
                } else if truth(RelOp::Ge, &p, pt)? && truth(RelOp::Lt, &p, &add2(pt, &int(1)))? {
                    ind = add2(&ind, &t);
                }
            }
            term = add2(&term, &ind);
        }
        Ok(collect(&term, &self.x))
    }

    /// `polynomial(n)`
    pub fn polynomial<H: Hooks>(&self, h: &mut H, n: &Expr) -> R<Expr, H::Error> {
        let mut terms = Vec::new();
        let sym = self.free_symbols();
        let n1 = add2(n, &int(1));
        let mut i: i64 = 0;
        loop {
            h_tick(h)?;
            let t = self.term(h, &int(i))?;
            let xp = self.pow_x_numeric(&t, &sym);
            if truth(RelOp::Ge, &xp, n)? {
                break;
            } else if crate::assumptions::is_integer(&xp) == Some(true) && int(i) == n1 {
                break;
            } else if !is_exact_zero(&t) {
                terms.push(t);
            }
            i += 1;
        }
        Ok(add(&terms))
    }

    /// `truncate(n)` for an integer `n`.
    pub fn truncate<H: Hooks>(&self, h: &mut H, n: &Expr) -> R<Expr, H::Error> {
        let pt_xk = self.xk.coeff(n)?;
        let x0 = if self.x0 == Expr::neg_oo() { Expr::oo() } else { self.x0.clone() };
        let sym = self.free_symbols();
        let pexp = self.pow_x(&pt_xk);
        // At +-oo the powers decrease and the reference terms are complete.
        let finite = self.x0 != Expr::oo() && self.x0 != Expr::neg_oo();
        let poly = if finite && has_any(&pexp, &sym) {
            self.polynomial_below(h, &pexp, n, &sym)?
        } else {
            self.polynomial(h, n)?
        };
        Ok(add2(&poly, &crate::series::order::order(&pt_xk, &[(self.x.clone(), x0)])))
    }

    /// The terms below the order term `x**pexp` of `truncate` when `pexp`
    /// is symbolic. Deliberate deviation: the reference compares the
    /// numeric part of each power with `n`, which drops terms its order
    /// term claims (`fps(x**(n + 1)*exp(x), x).truncate(3)` loses
    /// `x**(n + 3)/2` but keeps `O(x**(n + 4))`); the power is compared with
    /// `pexp` instead. For a numeric `pexp` this is `polynomial(n)`. Used at
    /// a finite point only.
    fn polynomial_below<H: Hooks>(&self, h: &mut H, pexp: &Expr, n: &Expr, sym: &[Expr]) -> R<Expr, H::Error> {
        let mut terms = Vec::new();
        let n1 = add2(n, &int(1));
        let mut i: i64 = 0;
        loop {
            h_tick(h)?;
            let t = self.term(h, &int(i))?;
            if !is_exact_zero(&t) {
                let gap = add2(&self.pow_x(&t), &mul2(&int(-1), pexp));
                let gap = if has_any(&gap, sym) { as_coeff_add(&gap, sym).0 } else { gap };
                if truth(RelOp::Ge, &gap, &Expr::zero())? {
                    break;
                }
                terms.push(t);
            } else if crate::assumptions::is_integer(&self.pow_x_numeric(&t, sym)) == Some(true) && int(i) == n1 {
                break;
            }
            i += 1;
        }
        Ok(add(&terms))
    }

    /// `_eval_as_leading_term`: the first nonzero term.
    pub fn leading_term<H: Hooks>(&self, h: &mut H) -> R<Expr, H::Error> {
        let mut i = 0;
        loop {
            h_tick(h)?;
            let t = self.term(h, &int(i))?;
            if !is_exact_zero(&t) {
                return Ok(t);
            }
            i += 1;
        }
    }

    /// `infinite`
    pub fn infinite(&self) -> Result<Expr, Exc> {
        // Sum(...) folds a Piecewise summand (AddWithLimits.__new__).
        let body = crate::integrate::piecewise::piecewise_fold(&mul2(&self.ak.formula, &self.xk.formula), Some(true));
        let s = crate::binders::held(Func::Sum, &body, vec![vec![self.ak.k.clone(), self.ak.start.clone(), self.ak.stop.clone()]])
            .map_err(Exc::Value)?;
        Ok(add2(&self.ind, &s))
    }

    /// `self + other` for two series: a series, or the sum of the functions
    /// when it does not depend on `x`.
    pub fn add_fps(&self, other: &Fps) -> Result<Result<Fps, Expr>, Exc> {
        if self.dir != other.dir {
            return Err(Exc::Value("Both series should be calculated from the same direction.".into()));
        }
        if self.x0 != other.x0 {
            return Err(Exc::Value("Both series should be calculated about the same point.".into()));
        }
        let f = add2(&self.function, &subs(&other.function, &[(other.x.clone(), self.x.clone())]));
        if !free_symbols(&f).contains(&self.x) {
            return Ok(Err(f));
        }
        let ak = self.ak.add(&other.ak);
        let (seq, s, e) = if truth(RelOp::Gt, &self.ak.start, &other.ak.start)? {
            (&other.ak, other.ak.start.clone(), self.ak.start.clone())
        } else {
            (&self.ak, self.ak.start.clone(), other.ak.start.clone())
        };
        let save = save_terms(seq, &self.xk, &s, &e)?;
        let ind = add(&[self.ind.clone(), other.ind.clone(), save]);
        Ok(Ok(self.with(f, ak, self.xk.clone(), ind)))
    }

    /// `self + other` for an `other` free of `x`; `None` means `Add(self, other)`.
    pub fn add_expr(&self, other: &Expr) -> Option<Fps> {
        if has(other, &self.x) {
            return None;
        }
        let f = add2(&self.function, other);
        let ind = add2(&self.ind, other);
        Some(self.with(f, self.ak.clone(), self.xk.clone(), ind))
    }

    /// `-self`
    pub fn neg(&self) -> Fps {
        self.with(neg(&self.function), self.ak.coeff_mul(&int(-1)), self.xk.clone(), neg(&self.ind))
    }

    /// `self*other`; `None` means `Mul(self, other)`.
    pub fn mul_expr(&self, other: &Expr) -> Option<Fps> {
        if has(other, &self.x) {
            return None;
        }
        Some(self.with(mul2(&self.function, other), self.ak.coeff_mul(other), self.xk.clone(), mul2(&self.ind, other)))
    }

    fn piecewise_terms(&self, scale: &dyn Fn(&Expr) -> Expr) -> Result<Expr, Exc> {
        let Kind::Piecewise(pairs) = self.ak.formula.kind() else {
            return Err(Exc::Type(format!("cannot unpack non-iterable {} object", self.ak.formula)));
        };
        let mut form = Vec::new();
        for (e, c) in pairs {
            let mut temp = Expr::zero();
            for t in make_add(e) {
                temp = add2(&temp, &scale(&t));
            }
            form.push((temp, c.clone()));
        }
        Ok(crate::relational::piecewise(form))
    }

    /// `_eval_derivative(v)`
    pub fn derivative(&self, v: &Expr) -> Result<Fps, Exc> {
        let f = crate::diff::diff(&self.function, v);
        let ind = crate::diff::diff(&self.ind, v);
        let pow_xk = self.pow_x(&self.xk.formula);
        let k = &self.ak.k;
        let shift = [(k.clone(), add2(k, &int(1)))];
        let formula = if has(&self.ak.formula, v) {
            let form = self.piecewise_terms(&|t: &Expr| mul2(t, &add2(&pow_xk, &self.pow_x(t))))?;
            subs(&form, &shift)
        } else {
            subs(&mul2(&self.ak.formula, &pow_xk), &shift)
        };
        let ak = Seq { formula, k: k.clone(), start: sub(&self.ak.start, &int(1)), stop: self.ak.stop.clone() };
        Ok(self.with(f, ak, self.xk.clone(), ind))
    }

    /// `integrate(v)`
    pub fn integrate<H: Hooks>(&self, h: &mut H, v: &Expr) -> R<Fps, H::Error> {
        let f = lift(h.integrate(&self.function, v))?;
        let mut ind = lift(h.integrate(&self.ind, v))?;
        ind = add2(&ind, &lift(h.limit0(&sub(&f, &ind), v))?);
        let pow_xk = self.pow_x(&self.xk.formula);
        let k = &self.ak.k;
        let shift = [(k.clone(), sub(k, &int(1)))];
        let formula = if has(&self.ak.formula, v) {
            let form = self.piecewise_terms(&|t: &Expr| div(t, &add(&[pow_xk.clone(), self.pow_x(t), int(1)])))?;
            subs(&form, &shift)
        } else {
            subs(&div(&self.ak.formula, &add2(&pow_xk, &int(1))), &shift)
        };
        let ak = Seq { formula, k: k.clone(), start: add2(&self.ak.start, &int(1)), stop: self.ak.stop.clone() };
        Ok(self.with(f, ak, self.xk.clone(), ind))
    }

    /// `zero_coeff()`
    pub fn zero_coeff<H: Hooks>(&self, h: &mut H) -> R<Expr, H::Error> {
        self.term(h, &Expr::zero())
    }

    /// `coeff_bell(n)` as its list of inner coefficients.
    fn coeff_bell(&self, n: i64) -> Vec<Expr> {
        (1..=n)
            .map(|j| {
                let xs: Vec<Expr> = (1..=(n - j + 1)).map(|pt| bell_coeff(&self.ak, pt)).collect();
                crate::special::bell_incomplete(&int(n), &int(j), &xs)
            })
            .collect()
    }
}

/// The checks of `product` and `compose`.
pub fn check_pair(a: &Fps, b: &Fps) -> Result<(), Exc> {
    if a.dir != b.dir {
        return Err(Exc::Value("Both series should be calculated from the same direction.".into()));
    }
    if a.x0 != b.x0 {
        return Err(Exc::Value("Both series should be calculated about the same point.".into()));
    }
    if a.x != b.x {
        return Err(Exc::Value("Both series should have the same symbol.".into()));
    }
    Ok(())
}

/// The checks of `compose(other)`.
pub fn check_compose<H: Hooks>(h: &mut H, a: &Fps, b: &Fps) -> R<(), H::Error> {
    check_pair(a, b)?;
    let t0 = b.term(h, &Expr::zero())?;
    if !is_exact_zero(&as_coeff_mul(&t0, std::slice::from_ref(&b.x), true).0) {
        return Err(Exc::Value("The formal power series of the inner function should not have any constant coefficient term.".into()).into());
    }
    Ok(())
}

/// The check of `inverse()`.
pub fn check_inverse<H: Hooks>(h: &mut H, a: &Fps) -> R<(), H::Error> {
    let t0 = a.term(h, &Expr::zero())?;
    if crate::assumptions::is_zero(&t0) == Some(true) {
        return Err(Exc::Value("Constant coefficient should exist for an inverse of a formal power series to exist.".into()).into());
    }
    Ok(())
}

/// `FormalPowerSeriesProduct(f, g)._eval_terms(n)`. The reference
/// convolution runs an exact FFT that rejects coefficients containing a
/// `Symbol`; the same coefficients are convolved directly here, with the
/// same rejection.
pub fn product_terms<H: Hooks>(h: &mut H, a: &Fps, b: &Fps, n: i64) -> R<Expr, H::Error> {
    let _ = h;
    let c1: Vec<Expr> = (0..n.max(0)).map(|i| subs(&a.ak.formula, &[(a.ak.k.clone(), int(i))])).collect();
    let c2: Vec<Expr> = (0..n.max(0)).map(|i| subs(&b.ak.formula, &[(b.ak.k.clone(), int(i))])).collect();
    if c1.iter().chain(c2.iter()).any(has_symbol_atom) {
        return Err(Exc::Value("Expected non-symbolic coefficients".into()).into());
    }
    let mut terms = Vec::new();
    for i in 0..n.max(0) {
        let mut conv = Vec::new();
        for j in 0..=i {
            conv.push(crate::simptools::expand_mul(&mul2(&c1[j as usize], &c2[(i - j) as usize])));
        }
        let aks = add(&conv);
        terms.push(mul2(&aks, &a.xk.coeff(&int(i))?));
    }
    Ok(add(&terms))
}

/// `FormalPowerSeriesCompose(f, g)._eval_terms(n)`
pub fn compose_terms<H: Hooks>(h: &mut H, a: &Fps, b: &Fps, n: i64) -> R<Expr, H::Error> {
    let mut terms = vec![a.zero_coeff(h)?];
    for i in 1..n {
        h_tick(h)?;
        let inner = b.coeff_bell(i);
        let seq: Vec<Expr> = (1..=i).map(|pt| mul2(&bell_coeff(&a.ak, pt), &inner[(pt - 1) as usize])).collect();
        let t = mul2(&div(&add(&seq), &factorial(&int(i))), &a.xk.coeff(&int(i))?);
        terms.push(t);
    }
    Ok(add(&terms))
}

/// `FormalPowerSeriesInverse(f)._eval_terms(n)`
pub fn inverse_terms<H: Hooks>(h: &mut H, a: &Fps, n: i64) -> R<Expr, H::Error> {
    let inv = a.zero_coeff(h)?;
    let kx = a.xk.k.clone();
    let inv_formula = pow(&inv, &neg(&add2(&kx, &int(1))));
    let ka = a.ak.k.clone();
    let fact_inv = mul2(&factorial(&ka), &subs(&inv_formula, &[(kx.clone(), ka.clone())]));
    let aux = |pt: i64| -> Expr {
        let sign = if pt % 2 == 1 { int(-1) } else { int(1) };
        mul2(&sign, &subs(&fact_inv, &[(ka.clone(), int(pt))]))
    };
    let mut terms = vec![inv.clone()];
    for i in 1..n {
        h_tick(h)?;
        let inner = a.coeff_bell(i);
        let seq: Vec<Expr> = (1..=i).map(|pt| mul2(&aux(pt), &inner[(pt - 1) as usize])).collect();
        let t = mul2(&div(&add(&seq), &factorial(&int(i))), &a.xk.coeff(&int(i))?);
        terms.push(t);
    }
    Ok(add(&terms))
}

/// `FiniteFormalPowerSeries.truncate(n)`: `terms + O(ffps.xk.coeff(n), (x, x0))`.
pub fn finite_truncate(a: &Fps, terms: &Expr, n: &Expr) -> Result<Expr, Exc> {
    let pt_xk = a.xk.coeff(n)?;
    Ok(add2(terms, &crate::series::order::order(&pt_xk, &[(a.x.clone(), a.x0.clone())])))
}

/// `Coeff.eval(p, x, n)`: `expand(p).coeff(x, n)` for a polynomial `p`
/// and an integer `n`, `None` otherwise.
pub fn coeff_eval(p: &Expr, x: &Expr, n: &Expr) -> Result<Option<Expr>, Exc> {
    if !(is_polynomial(p, x) && crate::assumptions::is_integer(n) == Some(true)) {
        return Ok(None);
    }
    match n.as_i64() {
        // Deliberate deviation: the reference takes `p.coeff(x, n)` of the
        // unexpanded polynomial, which misses the coefficients hidden in
        // products and powers ((x + 1)**2 has no x**2 term for it). The
        // coefficient of the expanded polynomial is the Taylor coefficient.
        Some(i) if n.is_integer_number() => Ok(Some(crate::simptools::coeff(&expand(p), x, i))),
        _ => Err(Exc::Value(format!("{n} is not an integer"))),
    }
}
