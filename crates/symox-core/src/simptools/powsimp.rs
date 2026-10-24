//! `sympy.simplify.powsimp`: `powsimp` and `powdenest`.
use super::exprtools::factor_terms;
use super::{add_args, args, keep_coeff, map_args, mul_args, nodes};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::sortkey::cmp_default;
use crate::{assumptions, decompose, ops};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Combine {
    All,
    Exp,
    Base,
}

#[derive(Clone, Copy, Debug)]
pub struct PowsimpOptions {
    pub deep: bool,
    pub combine: Combine,
    pub force: bool,
}

impl Default for PowsimpOptions {
    fn default() -> Self {
        PowsimpOptions { deep: false, combine: Combine::All, force: false }
    }
}

/// Insertion-ordered map keyed by expressions.
#[derive(Clone, Default)]
struct OMap<K: PartialEq + Clone, V> {
    items: Vec<(K, V)>,
}

impl<K: PartialEq + Clone, V> OMap<K, V> {
    fn new() -> Self {
        OMap { items: Vec::new() }
    }
    fn get(&self, k: &K) -> Option<&V> {
        self.items.iter().find(|(a, _)| a == k).map(|(_, v)| v)
    }
    fn get_mut(&mut self, k: &K) -> Option<&mut V> {
        self.items.iter_mut().find(|(a, _)| a == k).map(|(_, v)| v)
    }
    fn contains(&self, k: &K) -> bool {
        self.get(k).is_some()
    }
    fn pop(&mut self, k: &K) -> Option<V> {
        let i = self.items.iter().position(|(a, _)| a == k)?;
        Some(self.items.remove(i).1)
    }
    fn set(&mut self, k: K, v: V) {
        match self.get_mut(&k) {
            Some(slot) => *slot = v,
            None => self.items.push((k, v)),
        }
    }
    fn entry(&mut self, k: K, default: impl FnOnce() -> V) -> &mut V {
        if let Some(i) = self.items.iter().position(|(a, _)| *a == k) {
            &mut self.items[i].1
        } else {
            self.items.push((k, default()));
            &mut self.items.last_mut().unwrap().1
        }
    }
}

/// `as_base_exp` with SymPy's `exp_polar.as_base_exp`: `exp_polar(x)` is
/// `exp_polar(1)**x` (ExpBase), except `exp_polar(0)`.
fn base_exp(e: &Expr) -> (Expr, Expr) {
    if let Kind::Func(Func::ExpPolar, a) = e.kind() {
        if a.len() == 1 && !a[0].is_zero() {
            return (crate::functions::func(Func::ExpPolar, vec![Expr::one()]), a[0].clone());
        }
    }
    decompose::as_base_exp(e)
}

fn is_polar(e: &Expr) -> bool {
    crate::special::polar::is_polar(e)
}

fn is_exp(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Exp, _))
}


fn commutative(e: &Expr) -> bool {
    decompose::commutativity(e) != Some(false)
}

fn int(n: impl Into<BigInt>) -> Expr {
    Expr::num(Number::from_bigint(n.into()))
}

fn big(e: &Expr) -> Option<BigInt> {
    e.as_number().filter(|n| n.is_integer()).and_then(Number::to_bigint)
}

fn multiplicity_int(p: &BigInt, n: &BigInt) -> i64 {
    let n = n.abs();
    if n.is_zero() || *p < BigInt::from(2) {
        return 0;
    }
    let mut m = 0;
    let mut n = n;
    while (&n % p).is_zero() {
        n /= p;
        m += 1;
    }
    m
}

/// `multiplicity(p, n)` for positive rationals `p` and `n`.
fn multiplicity(p: &Number, n: &Number) -> i64 {
    let (pp, pq) = p.numer_denom().unwrap();
    let (np, nq) = n.numer_denom().unwrap();
    if pq.is_one() {
        if np.is_one() {
            return -multiplicity_int(&pp, &nq);
        }
        return multiplicity_int(&pp, &np) - multiplicity_int(&pp, &nq);
    }
    if pp.is_one() {
        return multiplicity_int(&pq, &nq);
    }
    let like = multiplicity_int(&pp, &np).min(multiplicity_int(&pq, &nq));
    let cross = multiplicity_int(&pq, &np).min(multiplicity_int(&pp, &nq));
    like - cross
}

fn cmp_seq(a: &[Expr], b: &[Expr]) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| {
        for (x, y) in a.iter().zip(b) {
            let c = cmp_default(x, y);
            if c != Ordering::Equal {
                return c;
            }
        }
        Ordering::Equal
    })
}

/// `bkey(b, e)` of powsimp: `((base, q), p)`.
fn bkey(b: &Expr, e: &Expr) -> ((Expr, BigInt), Expr) {
    if e.is_integer_number() {
        return ((b.clone(), BigInt::one()), e.clone());
    }
    if let Some((p, q)) = e.as_number().filter(|n| n.is_rational()).and_then(Number::numer_denom) {
        return ((b.clone(), q), int(p));
    }
    let (c, m) = decompose::as_coeff_mul_head(e, true);
    if !c.is_one() {
        let (cp, cq) = c.as_number().and_then(Number::numer_denom).unwrap_or((BigInt::one(), BigInt::one()));
        if assumptions::is_integer(&m) == Some(true) {
            return ((b.clone(), cq), ops::mul(&[m, int(cp)]));
        }
        return ((ops::pow(b, &m), cq), int(cp));
    }
    ((ops::pow(b, e), BigInt::one()), Expr::one())
}

fn bkey1(x: &Expr) -> ((Expr, BigInt), Expr) {
    let (b, e) = base_exp(x);
    bkey(&b, &e)
}

fn ratq(x: &Expr) -> BigInt {
    bkey1(x).0 .1
}

fn num_lt(a: &Expr, b: &Expr) -> Option<bool> {
    let (x, y) = (a.as_number()?, b.as_number()?);
    Some(x.cmp_value(y)? == Ordering::Less)
}

fn floordiv(a: &Expr, b: &Expr) -> Option<BigInt> {
    let (x, y) = (big(a)?, big(b)?);
    if y.is_zero() {
        return None;
    }
    Some(x.div_floor(&y))
}

/// `powsimp(expr, deep, combine, force)` with `measure=count_ops`.
pub fn powsimp(expr: &Expr, opts: &PowsimpOptions) -> Expr {
    let y = Expr::dummy("y", crate::symbol::Assumptions(0));
    powsimp_y(expr, opts, &y)
}

fn powsimp_y(expr: &Expr, opts: &PowsimpOptions, y: &Expr) -> Expr {
    if expr.is_atom() || matches!(expr.kind(), Kind::Func(Func::MatrixProxy(_), _)) {
        return expr.clone();
    }
    // `expr in (exp_polar(0), exp_polar(1))`
    if let Kind::Func(Func::ExpPolar, a) = expr.kind() {
        if a.len() == 1 && (a[0].is_zero() || a[0].is_one()) {
            return expr.clone();
        }
    }
    if !matches!(expr.kind(), Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Func(..) | Kind::Piecewise(_) | Kind::Rel(..) | Kind::Derivative(..)) {
        return expr.clone();
    }
    let recurse = |a: &Expr, o: PowsimpOptions| powsimp_y(a, &o, y);
    let mut expr = expr.clone();
    if opts.deep || expr.is_add() || (expr.is_mul() && !args(&expr).contains(y)) {
        expr = map_args(&expr, &mut |w| recurse(w, *opts));
    }
    if expr.is_pow() {
        let o = PowsimpOptions { deep: false, ..*opts };
        return ops::div(&recurse(&ops::mul(&[expr.clone(), y.clone()]), o), y);
    }
    if !expr.is_mul() {
        return expr;
    }
    match opts.combine {
        Combine::Exp | Combine::All => combine_exp(&expr, opts, y),
        Combine::Base => combine_base(&expr, opts, y),
    }
}

fn combine_exp(expr: &Expr, opts: &PowsimpOptions, y: &Expr) -> Expr {
    let recurse = |a: &Expr, o: PowsimpOptions| powsimp_y(a, &o, y);
    let mut c_powers: OMap<Expr, Vec<Expr>> = OMap::new();
    let mut nc_part: Vec<Expr> = Vec::new();
    let mut coeff = Number::ONE;
    for term in args(expr) {
        if let Some(n) = term.as_number().filter(|n| n.is_rational()) {
            coeff = coeff.mul(n);
            continue;
        }
        let term = if term.is_pow() { denest_pow(&term) } else { term };
        if commutative(&term) {
            let (mut b, mut e) = base_exp(&term);
            if opts.deep {
                b = recurse(&b, *opts);
                e = recurse(&e, *opts);
            }
            if b.is_pow() || is_exp(&b) {
                b = ops::pow(&b, &e);
                e = Expr::one();
            }
            c_powers.entry(b, Vec::new).push(e);
        } else {
            if let Some(last) = nc_part.last() {
                let (b1, e1) = decompose::as_base_exp(last);
                let (b2, e2) = decompose::as_base_exp(&term);
                if b1 == b2 && commutative(&e1) && commutative(&e2) {
                    let n = nc_part.len();
                    nc_part[n - 1] = ops::pow(&b1, &ops::add(&[e1, e2]));
                    continue;
                }
            }
            nc_part.push(term);
        }
    }
    // `ordered(c_powers.items())`
    let mut order: Vec<usize> = (0..c_powers.items.len()).collect();
    let key_nodes = |i: usize| 1 + nodes(&c_powers.items[i].0) + 1 + c_powers.items[i].1.iter().map(nodes).sum::<usize>();
    order.sort_by(|&i, &j| {
        key_nodes(i)
            .cmp(&key_nodes(j))
            .then_with(|| cmp_default(&c_powers.items[i].0, &c_powers.items[j].0))
            .then_with(|| cmp_seq(&c_powers.items[i].1, &c_powers.items[j].1))
    });
    let mut summed: Vec<Option<Expr>> = vec![None; c_powers.items.len()];
    for i in order {
        let (b, e) = &mut c_powers.items[i];
        if let Some(bn) = b.as_number().filter(|n| n.is_rational() && !n.is_zero()) {
            if !e.iter().all(super::is_number_atom) && !coeff.is_one() && !bn.is_one() && !bn.is_minus_one() {
                let m = multiplicity(&bn.abs(), &coeff.abs());
                if m != 0 {
                    e.push(Expr::int(m));
                    coeff = coeff.div(&bn.pow_i64(m).unwrap()).unwrap();
                }
            }
        }
        summed[i] = Some(ops::add(e));
    }
    let mut cp: OMap<Expr, Expr> = OMap::new();
    for ((b, _), s) in c_powers.items.into_iter().zip(summed) {
        cp.items.push((b, s.unwrap()));
    }
    if !coeff.is_one() {
        let c = Expr::num(coeff.clone());
        if let Some(v) = cp.get_mut(&c) {
            *v = ops::add(&[v.clone(), Expr::one()]);
        } else {
            cp.set(c, Expr::one());
        }
    }
    // Combine b**e with (1/b)**e2.
    let be: Vec<(Expr, Expr)> = cp.items.clone();
    let mut skip: Vec<Expr> = Vec::new();
    for (b, e) in &be {
        if skip.contains(b) {
            continue;
        }
        if assumptions::is_positive(b) == Some(true) || is_polar(b) {
            let binv = ops::inv(b);
            if matches!(b.kind(), Kind::Num(Number::Float(_))) && b.as_number().unwrap().cmp_value(&Number::ONE) == Some(Ordering::Equal) {
                cp.set(b.clone(), Expr::one());
                continue;
            }
            if *b != binv && cp.contains(&binv) {
                if crate::simptools::as_numer_denom(b).0.is_one() {
                    cp.pop(b);
                    let v = cp.get_mut(&binv).unwrap();
                    *v = ops::sub(v, e);
                } else {
                    skip.push(binv.clone());
                    let e2 = cp.pop(&binv).unwrap();
                    if let Some(v) = cp.get_mut(b) {
                        *v = ops::sub(v, &e2);
                    }
                }
            }
        }
    }
    // Combine (-b)**e with b**e.
    let be: Vec<(Expr, Expr)> = cp.items.clone();
    let n1 = Expr::minus_one();
    for (b, e) in &be {
        if !(b.is_symbol() || b.is_add()) {
            continue;
        }
        let nb = ops::neg(b);
        if !(cp.contains(&nb) && cp.contains(b)) {
            continue;
        }
        if assumptions::is_positive(b).is_some() || assumptions::is_integer(e) == Some(true) {
            let moved = if assumptions::is_integer(e) == Some(true) || assumptions::is_negative(b) == Some(true) {
                // Python's `c_powers[-b] += c_powers.pop(b)` reads the slot
                // first, so -b == b (a sum of zoo terms) doubles it.
                let old = cp.get_mut(&nb).unwrap().clone();
                let v = cp.pop(b).unwrap();
                cp.set(nb.clone(), ops::add(&[old, v]));
                e.clone()
            } else {
                let v = cp.pop(&nb).unwrap();
                let old = cp.get_mut(b).map(|slot| slot.clone()).unwrap_or_else(|| v.clone());
                cp.set(b.clone(), ops::add(&[old, v.clone()]));
                v
            };
            if let Some(slot) = cp.get_mut(&n1) {
                *slot = ops::add(&[slot.clone(), moved]);
            } else {
                cp.set(n1.clone(), moved);
            }
        }
    }
    let c_list: Vec<(Expr, Expr)> = cp.items.into_iter().filter(|(_, e)| !is_false_zero(e)).collect();

    // Rational-exponent base collection.
    let mut common_b: OMap<(Expr, BigInt), Expr> = OMap::new();
    let mut bases: Vec<(Expr, BigInt)> = Vec::new();
    for (b, e) in &c_list {
        let (k, e) = bkey(b, e);
        if let Some(v) = common_b.get_mut(&k) {
            *v = ops::add(&[v.clone(), e]);
        } else {
            common_b.set(k.clone(), e);
        }
        if !k.1.is_one() && k.0.is_mul() {
            bases.push(k);
        }
    }
    bases.sort_by(|a, b| cmp_default(&a.0, &b.0).then_with(|| a.1.cmp(&b.1)));
    let measure = |k: &(Expr, BigInt)| super::methods::count_ops(&k.0);
    bases.sort_by(|a, b| measure(b).cmp(&measure(a)));
    let mut done: Vec<(Expr, Expr)> = Vec::new();
    let mut base_idx = 0;
    while base_idx < bases.len() {
        let base = bases[base_idx].clone();
        base_idx += 1;
        if !common_b.contains(&base) {
            continue;
        }
        let (mut b, exponent) = base.clone();
        let mut last = false;
        let mut qlcm = BigInt::one();
        loop {
            let bstart = b.clone();
            let qstart = qlcm.clone();
            let mut bb: Vec<(Expr, BigInt)> = Vec::new();
            let mut ee: Vec<(Expr, Expr)> = Vec::new();
            for bi in mul_args(&b) {
                let (bib, bie) = bkey1(&bi);
                let ok = match common_b.get(&bib) {
                    None => false,
                    Some(v) => num_lt(v, &bie) == Some(false),
                };
                if !ok {
                    ee.clear();
                    bb.clear();
                    break;
                }
                ee.push((bie, common_b.get(&bib).unwrap().clone()));
                bb.push(bib);
            }
            if !ee.is_empty() {
                let mut min1 = floordiv(&ee[0].1, &ee[0].0).unwrap_or_default();
                let mut broke = false;
                for item in ee.iter().skip(1) {
                    let rat = floordiv(&item.1, &item.0).unwrap_or_default();
                    if rat < BigInt::one() {
                        broke = true;
                        break;
                    }
                    min1 = min1.min(rat);
                }
                if !broke {
                    for i in 0..bb.len() {
                        let v = common_b.get_mut(&bb[i]).unwrap();
                        *v = ops::sub(v, &ops::mul(&[int(min1.clone()), ee[i].0.clone()]));
                        update(&mut common_b, &mut bases, &bb[i]);
                    }
                    let v = common_b.entry(base.clone(), Expr::zero);
                    *v = ops::add(&[v.clone(), int(&min1 * &qstart * &exponent)]);
                }
            }
            if last || common_b.items.len() == 1 || common_b.items.iter().all(|(k, _)| k.1.is_one()) {
                break;
            }
            qlcm = mul_args(&bstart).iter().map(ratq).fold(BigInt::one(), |a, q| a.lcm(&q));
            if qlcm.is_one() {
                break;
            }
            b = ops::pow(&bstart, &int(qlcm.clone()));
            qlcm *= &qstart;
            if mul_args(&b).iter().all(|bi| ratq(bi).is_one()) {
                last = true;
            }
        }
        let (b, q) = base.clone();
        let v = common_b.pop(&base).unwrap_or_else(Expr::zero);
        done.push((b, ops::mul(&[v, Expr::num(Number::rational(BigInt::one(), q).unwrap())])));
    }
    let mut c_final = done;
    for ((b, q), e) in common_b.items {
        let nb = if (b.is_pow() || is_exp(&b)) && !q.is_one() && !decompose::as_base_exp(&b).1.is_rational_number() {
            let (bb, be) = decompose::as_base_exp(&b);
            ops::pow(&bb, &ops::div(&be, &int(q)))
        } else {
            ops::pow(&b, &Expr::num(Number::rational(BigInt::one(), q).unwrap()))
        };
        c_final.push((nb, e));
    }
    let mut merged: OMap<Expr, Expr> = OMap::new();
    for (b, e) in c_final {
        merged.set(b, e);
    }
    let parts: Vec<Expr> = merged.items.iter().map(|(b, e)| ops::pow(b, e)).collect();
    let newexpr = ops::mul(&parts);
    if opts.combine == Combine::Exp {
        ops::mul(&[newexpr, ops::mul(&nc_part)])
    } else {
        let o = PowsimpOptions { combine: Combine::Base, ..*opts };
        ops::mul(&[powsimp_y(&ops::mul(&nc_part), &o, y), powsimp_y(&newexpr, &o, y)])
    }
}

fn is_false_zero(e: &Expr) -> bool {
    match e.as_number() {
        Some(n) => n.is_zero(),
        None => false,
    }
}

fn update(common_b: &mut OMap<(Expr, BigInt), Expr>, bases: &mut Vec<(Expr, BigInt)>, b: &(Expr, BigInt)) {
    let Some(v) = common_b.get(b).cloned() else { return };
    let Some(vb) = big(&v) else { return };
    let (newe, r) = vb.div_mod_floor(&b.1);
    if !r.is_zero() {
        return;
    }
    common_b.pop(b);
    if !newe.is_zero() {
        for m in mul_args(&ops::pow(&b.0, &int(newe))) {
            let (k, e) = bkey1(&m);
            let slot = common_b.entry(k.clone(), Expr::zero);
            *slot = ops::add(&[slot.clone(), e]);
            if !k.1.is_one() {
                bases.push(k);
            }
        }
    }
}

fn nonneg_all(b: &Expr) -> bool {
    let (n, d) = crate::simptools::as_numer_denom(b);
    assumptions::is_nonnegative(&n) == Some(true) && assumptions::is_nonnegative(&d) == Some(true)
}

fn coeff_isneg(e: &Expr) -> bool {
    let first = if e.is_mul() { args(e)[0].clone() } else { e.clone() };
    first.as_number().is_some_and(|n| n.is_negative()) || first.is_constant(Constant::NegInfinity)
}

fn count_terms(e: &Expr) -> usize {
    if e.is_add() {
        return args(e).iter().map(count_terms).sum();
    }
    if e.is_mul() {
        return args(e).iter().map(count_terms).product();
    }
    1
}

fn combine_base(expr: &Expr, opts: &PowsimpOptions, y: &Expr) -> Expr {
    let recurse = |a: &Expr| powsimp_y(a, opts, y);
    let mut c_powers: Vec<(Expr, Expr)> = Vec::new();
    let mut nc_part: Vec<Expr> = Vec::new();
    for term in args(expr) {
        if commutative(&term) {
            c_powers.push(base_exp(&term));
        } else {
            nc_part.push(term);
        }
    }
    for item in c_powers.iter_mut() {
        let (b, e) = item.clone();
        if !(nonneg_all(&b) || assumptions::is_integer(&e) == Some(true) || opts.force || is_polar(&b)) {
            continue;
        }
        let (exp_c, exp_t) = decompose::as_coeff_mul_head(&e, true);
        if !exp_c.is_one() && !exp_t.is_one() {
            *item = (ops::pow(&b, &exp_c), exp_t);
        }
    }
    let mut c_exp: OMap<Expr, Vec<Expr>> = OMap::new();
    for (b, e) in c_powers {
        let mut b = b;
        let mut e = if opts.deep { recurse(&e) } else { e };
        if e.is_add() && (assumptions::is_positive(&b) == Some(true) || assumptions::is_integer(&e) == Some(true)) {
            e = factor_terms(&e, true, false, true, true);
            if coeff_isneg(&e) {
                e = ops::neg(&e);
                b = ops::inv(&b);
            }
        }
        c_exp.entry(e, Vec::new).push(b);
    }
    let mut out: OMap<Expr, Vec<Expr>> = OMap::new();
    for (e, bases) in c_exp.items {
        let new_base = if bases.len() == 1 {
            bases[0].clone()
        } else if assumptions::is_integer(&e) == Some(true) || opts.force {
            ops::mul(&bases)
        } else {
            let mut unk: Vec<Expr> = Vec::new();
            let mut nonneg: Vec<Expr> = Vec::new();
            let mut neg: Vec<Expr> = Vec::new();
            for bi in bases {
                if assumptions::is_negative(&bi) == Some(true) {
                    neg.push(bi);
                } else if assumptions::is_nonnegative(&bi) == Some(true) || is_polar(&bi) {
                    // polar can be treated like non-negative
                    nonneg.push(bi);
                } else {
                    unk.push(bi);
                }
            }
            if (unk.len() == 1 && neg.is_empty()) || (neg.len() == 1 && unk.is_empty()) {
                nonneg.append(&mut unk);
                nonneg.append(&mut neg);
            } else if !neg.is_empty() {
                let israt = if e.is_rational_number() {
                    true
                } else {
                    let (p, d) = crate::simptools::as_numer_denom(&e);
                    assumptions::is_integer(&p) == Some(true) && assumptions::is_integer(&d) == Some(true)
                };
                if israt {
                    neg = neg.iter().map(ops::neg).collect();
                    for _ in 0..neg.len() {
                        unk.push(Expr::minus_one());
                    }
                } else {
                    unk.append(&mut neg);
                }
            }
            for b in &unk {
                out.entry(b.clone(), Vec::new).push(e.clone());
            }
            let mut all = nonneg;
            all.extend(neg);
            let mut new_base = ops::mul(&all);
            let xnew = super::expand_hints_shallow_mul(&new_base);
            if add_args(&xnew).len() < count_terms(&new_base) {
                new_base = factor_terms(&xnew, true, false, true, true);
            }
            new_base
        };
        out.entry(new_base, Vec::new).push(e);
    }
    let mut c_part: Vec<Expr> = Vec::new();
    for (b, es) in out.items {
        for e in es {
            c_part.push(ops::pow(&b, &e));
        }
    }
    if c_part.len() + nc_part.len() > 1 {
        // SymPy's Mul.flatten gathers each power's exponent as c*t from its
        // as_coeff_Mul, which distributes the coefficient factor_terms
        // pulled out above: 3*exp(3*(2*a + x**2)) becomes 3*exp(6*a + 3*x**2).
        c_part = c_part.into_iter().map(|p| match p.kind() {
            Kind::Func(Func::Exp, a) => match decompose::as_coeff_mul_head(&a[0], false) {
                (c, t) if !c.is_one() && t.is_add() => crate::functions::exp(&ops::mul(&[c, t])),
                _ => p.clone(),
            },
            Kind::Pow(b, e) => match decompose::as_coeff_mul_head(e, false) {
                (c, t) if !c.is_one() && t.is_add() => ops::pow(b, &ops::mul(&[c, t])),
                _ => p.clone(),
            },
            _ => p.clone(),
        }).collect();
    }
    c_part.extend(nc_part);
    ops::mul(&c_part)
}

/// `powdenest(eq, force)`.
pub fn powdenest(eq: &Expr, force: bool, _polar: bool) -> Expr {
    if force {
        fn denest(b: &Expr, e: &Expr) -> (Option<bool>, Expr) {
            match b.kind() {
                Kind::Pow(bb, be) => denest(bb, &ops::mul(&[be.clone(), e.clone()])),
                Kind::Func(Func::Exp, a) => denest(&Expr::constant(Constant::E), &ops::mul(&[a[0].clone(), e.clone()])),
                _ => (assumptions::is_positive(b), ops::pow(b, e)),
            }
        }
        let mut reps: Vec<(Expr, Expr)> = Vec::new();
        super::preorder(eq, &mut |p| {
            if let Kind::Pow(b, e) = p.kind() {
                if b.is_pow() || is_exp(b) {
                    let (ok, dp) = denest(b, e);
                    if ok != Some(false) && !reps.iter().any(|(k, _)| k == p) {
                        reps.push((p.clone(), dp));
                    }
                }
            }
        });
        let eq = if reps.is_empty() { eq.clone() } else { crate::subs::subs(eq, &reps) };
        let (eq, back) = super::misc::posify(&eq);
        let r = powdenest(&eq, false, false);
        return super::misc::xreplace(&r, &back);
    }
    let new = powsimp(eq, &PowsimpOptions::default());
    transform_pows(&new)
}

fn transform_pows(e: &Expr) -> Expr {
    if e.is_pow() || is_exp(e) {
        return denest_pow(e);
    }
    if e.is_atom() {
        return e.clone();
    }
    map_args(e, &mut transform_pows)
}

fn has_log_term(ei: &Expr) -> bool {
    add_args(ei).iter().any(|a| matches!(a.kind(), Kind::Func(Func::Log, _)))
}

/// `e` has a literal 2 as denominator (`_half` in `Pow._eval_power`).
fn half(e: &Expr) -> bool {
    if let Some(n) = e.as_number() {
        return n.is_rational() && n.denom() == Number::Int(2);
    }
    let (n, d) = crate::simptools::as_numer_denom(e);
    assumptions::is_integer(&n) == Some(true) && d == Expr::int(2)
}

/// `Pow._eval_power(b0, expt)` for a power (or exponential) `b0`, limited
/// to the branches that decide the branch factor without numerics.
fn pow_eval_power(b0: &Expr, expt: &Expr) -> Option<Expr> {
    let (mut b, e) = decompose::as_base_exp(b0);
    if b.is_constant(Constant::NaN) {
        return None;
    }
    let ext_real = if e.as_number().is_some() || e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity) {
        Some(true)
    } else {
        match assumptions::is_real(&e) {
            Some(true) => Some(true),
            Some(false) if assumptions::is_finite(&e) == Some(true) => Some(false),
            _ => None,
        }
    };
    let mut ok = assumptions::is_integer(expt) == Some(true);
    if !ok && ext_real == Some(true) {
        if e.is_minus_one() {
            if half(expt) {
                match assumptions::is_negative(&b) {
                    Some(true) => {
                        return Some(ops::mul(&[ops::pow(&Expr::minus_one(), expt), ops::pow(&ops::neg(&b), &ops::mul(&[e.clone(), expt.clone()]))]));
                    }
                    Some(false) => return Some(ops::pow(&b, &ops::neg(expt))),
                    None => {}
                }
            }
        } else if assumptions::is_even(&e) == Some(true) && assumptions::is_real(&b) == Some(true) {
            b = crate::functions::func(Func::Abs, vec![b]);
        }
        let small = e.as_number().and_then(|n| n.abs().cmp_value(&Number::ONE)).is_some_and(|o| o == Ordering::Less);
        if small || e.is_one() || assumptions::is_nonnegative(&b) == Some(true) {
            ok = true;
        }
    }
    ok.then(|| ops::pow(&b, &ops::mul(&[e, expt.clone()])))
}

/// `_denest_pow(eq)`.
pub fn denest_pow(eq: &Expr) -> Expr {
    let (b, e) = decompose::as_base_exp(eq);
    let new = if (b.is_pow() || is_exp(&b)) && !e.is_one() { pow_eval_power(&b, &e) } else { None };
    let eq = new.as_ref().unwrap_or(eq);
    let (b, e) = decompose::as_base_exp(eq);
    if b.is_constant(Constant::E) && e.is_mul() {
        let mut logs = Vec::new();
        let mut other = Vec::new();
        for ei in args(&e) {
            if has_log_term(&ei) {
                logs.push(ei);
            } else {
                other.push(ei);
            }
        }
        let logs = crate::simp::logcombine(&ops::mul(&logs), false);
        return ops::pow(&crate::functions::func(Func::Exp, vec![logs]), &ops::mul(&other));
    }
    let (_, be) = decompose::as_base_exp(&b);
    if be.is_one() && !(b.is_mul() || b.as_number().is_some_and(|n| n.is_rational() && !n.is_integer()) || assumptions::is_positive(&b) == Some(true)) {
        return eq.clone();
    }
    // polar factors denest freely
    let factors = mul_args(&b);
    let (polars, nonpolars): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(crate::special::polar::is_polar);
    if !polars.is_empty() {
        let rest = powdenest(&ops::pow(&ops::mul(&nonpolars), &e), false, false);
        let pbe: Vec<(Expr, Expr)> = polars.iter().map(decompose::as_base_exp).collect();
        if pbe.len() == 1 && !pbe[0].0.is_mul() && !matches!(pbe[0].0.kind(), Kind::Func(Func::ExpPolar, _)) {
            let (pb, pe) = &pbe[0];
            return ops::mul(&[ops::pow(pb, &ops::mul(&[pe.clone(), e.clone()])), rest]);
        }
        let mut out: Vec<Expr> = pbe
            .iter()
            .map(|(pb, pe)| match pb.kind() {
                // ExpBase.as_base_exp: exp_polar(x) = exp_polar(1)**x
                Kind::Func(Func::ExpPolar, a) if pe.is_one() => {
                    crate::functions::func(Func::ExpPolar, vec![ops::mul(&[a[0].clone(), e.clone()])])
                }
                _ => powdenest(&ops::pow(pb, &ops::mul(&[pe.clone(), e.clone()])), false, false),
            })
            .collect();
        out.push(rest);
        return ops::mul(&out);
    }
    if let Some(n) = big(&b) {
        if let Some((base, k)) = crate::ntheory::perfect_power(&n) {
            return ops::pow(&int(base), &ops::mul(&[e.clone(), Expr::int(k as i64)]));
        }
    }
    if !b.is_mul() || mul_args(&b).iter().any(|s| s.is_atom()) {
        return eq.clone();
    }
    let log = |x: &Expr| crate::functions::func(Func::Log, vec![x.clone()]);
    let mut glogb = crate::simp::expand_log(&log(&b), false);
    if glogb.is_add() {
        let terms = args(&glogb);
        let mut g = terms[0].clone();
        for t in &terms[1..] {
            g = nc_gcd(&g, t);
        }
        if !g.is_one() {
            let (cg, rg) = decompose::as_coeff_mul_head(&g, false);
            let inner: Vec<Expr> = terms.iter().map(|a| ops::div(a, &g)).collect();
            glogb = keep_coeff(&cg, &ops::mul(&[rg, ops::add(&inner)]), true, false);
        }
    }
    let is_log = matches!(glogb.kind(), Kind::Func(Func::Log, _));
    if is_log || !glogb.is_mul() {
        if let Some(a0) = args(&glogb).first() {
            if a0.is_pow() || is_exp(a0) {
                let g2 = denest_pow(a0);
                let (gb, ge) = decompose::as_base_exp(&g2);
                let absl = crate::functions::func(Func::Abs, vec![ge.clone()]);
                let lt = crate::relational::rel(crate::expr::RelOp::Lt, &absl, &Expr::one());
                if matches!(lt.kind(), Kind::Bool(true)) {
                    return ops::pow(&gb, &ops::mul(&[ge, e.clone()]));
                }
            }
        }
        return eq.clone();
    }
    let mut add = Vec::new();
    let mut other = Vec::new();
    for a in args(&glogb) {
        if a.is_add() {
            add.push(a);
        } else {
            other.push(a);
        }
    }
    let inner = crate::simp::logcombine(&ops::mul(&add), false);
    ops::pow(&crate::functions::func(Func::Exp, vec![inner]), &ops::mul(&[e.clone(), ops::mul(&other)]))
}

fn nc_gcd(aa: &Expr, bb: &Expr) -> Expr {
    let (a0, a1) = decompose::as_coeff_mul_head(aa, false);
    let (b0, b1) = decompose::as_coeff_mul_head(bb, false);
    let c = match (a0.as_number(), b0.as_number()) {
        (Some(x), Some(y)) if x.is_rational() && y.is_rational() => {
            let g = super::exprtools::rat_gcd(x, y);
            Expr::num(Number::from_bigint(g.numer_denom().map(|(p, _)| p).unwrap_or_else(BigInt::one)))
        }
        _ => crate::poly::api::gcd(&a0, &b0, &[]).map(|g| crate::simptools::as_numer_denom(&g).0).unwrap_or_else(|_| Expr::one()),
    };
    let ca = decompose::args_cnc(&a1, true).0;
    let cb = decompose::args_cnc(&b1, true).0;
    let mut common: Vec<Expr> = Vec::new();
    for x in ca {
        if cb.contains(&x) && !common.contains(&x) {
            common.push(x);
        }
    }
    keep_coeff(&c, &ops::mul(&common), true, false)
}
