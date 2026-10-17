//! `sympy.core.exprtools`: `Factors`, `Term`, `gcd_terms`, `factor_terms`
//! and `factor_nc` for commutative expressions.
use super::{add_args, add_from_args, args, as_coeff_mul, extract_additively, is_number_atom, keep_coeff, mul_args, primitive};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::{assumptions, decompose, ops, traverse};

/// `decompose_power(expr)`: `(base, integer exponent)`.
pub fn decompose_power(e: &Expr) -> (Expr, Number) {
    let (base, exp) = decompose::as_base_exp(e);
    if let Some(n) = exp.as_number() {
        if let Some((p, q)) = n.numer_denom() {
            let base = if q != 1.into() { ops::pow(&base, &Expr::num(Number::from_bigint(1.into()).div(&Number::from_bigint(q)).unwrap())) } else { base };
            return (base, Number::from_bigint(p));
        }
        return (e.clone(), Number::ONE);
    }
    let (c, tail) = as_coeff_mul(&exp, true);
    let c = c.as_number().cloned().unwrap_or(Number::ONE);
    if c.is_minus_one() {
        return (ops::pow(&base, &tail), Number::int(-1));
    }
    if !c.is_one() {
        let (p, q) = c.numer_denom().unwrap();
        let tail = keep_coeff(&Expr::num(Number::from_bigint(1.into()).div(&Number::from_bigint(q)).unwrap()), &tail, true, false);
        return (ops::pow(&base, &tail), Number::from_bigint(p));
    }
    (e.clone(), Number::ONE)
}

fn isnumber(e: &Expr) -> bool {
    is_number_atom(e)
}

/// Insertion-ordered `{base: exponent}` map of SymPy's `Factors`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Factors {
    pub factors: Vec<(Expr, Expr)>,
}

impl Factors {
    fn get(&self, k: &Expr) -> Option<&Expr> {
        self.factors.iter().find(|(f, _)| f == k).map(|(_, e)| e)
    }
    fn set(&mut self, k: Expr, v: Expr) {
        if let Some(slot) = self.factors.iter_mut().find(|(f, _)| *f == k) {
            slot.1 = v;
        } else {
            self.factors.push((k, v));
        }
    }
    fn remove(&mut self, k: &Expr) {
        self.factors.retain(|(f, _)| f != k);
    }
    fn add_to(&mut self, k: Expr, v: &Expr) {
        let cur = self.get(&k).cloned().unwrap_or_else(Expr::zero);
        self.set(k, ops::add2(&cur, v));
    }

    pub fn zero() -> Factors {
        Factors { factors: vec![(Expr::zero(), Expr::one())] }
    }

    /// `Factors(expr)`.
    pub fn from_expr(e: &Expr) -> Factors {
        let mut out = Factors::default();
        if e.is_one() {
            return out;
        }
        if e.is_zero() {
            return Factors::zero();
        }
        if let Some(n) = e.as_number() {
            let mut n = n.clone();
            if n.is_negative() {
                out.factors.push((Expr::minus_one(), Expr::one()));
                n = n.neg();
            }
            if !n.is_one() {
                if n.is_integer() || n.is_float() {
                    out.factors.push((Expr::num(n), Expr::one()));
                } else if let Some((p, q)) = n.numer_denom() {
                    if p != 1.into() {
                        out.factors.push((Expr::num(Number::from_bigint(p)), Expr::one()));
                    }
                    out.factors.push((Expr::num(Number::from_bigint(q)), Expr::minus_one()));
                }
            }
            return out;
        }
        if e.is_constant(Constant::Infinity) || (e.is_atom() && !e.is_constant(Constant::I)) {
            out.factors.push((e.clone(), Expr::one()));
            return out;
        }
        let (c, nc) = decompose::args_cnc(e, true);
        let i = c.iter().filter(|a| a.is_constant(Constant::I)).count();
        for a in c.iter().filter(|a| !a.is_constant(Constant::I)) {
            let (b, x) = decompose::as_base_exp(a);
            out.add_to(b, &x);
        }
        let keys: Vec<Expr> = out.factors.iter().map(|f| f.0.clone()).collect();
        for f in keys {
            if let Some(n) = f.as_number().filter(|n| n.is_rational() && !n.is_integer()) {
                let (p, q) = n.numer_denom().unwrap();
                let x = out.get(&f).unwrap().clone();
                let (p, q) = (Expr::num(Number::from_bigint(p)), Expr::num(Number::from_bigint(q)));
                out.add_to(p, &x);
                out.add_to(q, &ops::neg(&x));
                out.remove(&f);
            }
        }
        if i > 0 {
            out.add_to(Expr::i(), &Expr::int(i as i64));
        }
        if !nc.is_empty() {
            out.factors.push((super::mul_from_args(&nc), Expr::one()));
        }
        out
    }

    /// `Factors(dict)`: normalises the numeric keys I, -1 and 1.
    pub fn from_map(mut factors: Vec<(Expr, Expr)>) -> Factors {
        let handle: Vec<Expr> = factors
            .iter()
            .filter(|(k, v)| (k.is_constant(Constant::I) || k.is_one() || k.is_minus_one()) && isnumber(v))
            .map(|(k, _)| k.clone())
            .collect();
        if handle.is_empty() {
            return Factors { factors };
        }
        let mut i1 = Expr::one();
        for k in &handle {
            let pos = factors.iter().position(|(f, _)| f == k).unwrap();
            let (k, v) = factors.remove(pos);
            i1 = ops::mul(&[i1, ops::pow(&k, &v)]);
        }
        let mut out = Factors { factors };
        if !i1.is_one() {
            for a in mul_args(&i1) {
                if a.is_minus_one() || a.is_constant(Constant::I) {
                    out.set(a, Expr::one());
                } else if let Kind::Pow(b, x) = a.kind() {
                    out.add_to(b.clone(), x);
                } else if a.is_one() {
                    out.set(a, Expr::one());
                } else if let Some(n) = a.as_number().filter(|n| n.is_minus_one()) {
                    let _ = n;
                    out.set(Expr::one(), Expr::one());
                    out.set(Expr::minus_one(), Expr::one());
                } else {
                    out.set(a, Expr::one());
                }
            }
        }
        out
    }

    pub fn is_zero(&self) -> bool {
        self.factors.len() == 1 && self.factors[0].0.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.factors.is_empty()
    }

    pub fn as_expr(&self) -> Expr {
        let mut parts = Vec::new();
        for (factor, exp) in &self.factors {
            if !exp.is_one() {
                if exp.is_integer_number() {
                    let (b, e) = decompose::as_base_exp(factor);
                    let e = keep_coeff(exp, &e, true, false);
                    parts.push(ops::pow(&b, &e));
                } else {
                    parts.push(ops::pow(factor, exp));
                }
            } else {
                parts.push(factor.clone());
            }
        }
        ops::mul(&parts)
    }

    pub fn mul(&self, other: &Factors) -> Factors {
        if self.is_zero() || other.is_zero() {
            return Factors::zero();
        }
        let mut factors = self.factors.clone();
        for (factor, exp) in &other.factors {
            if let Some(pos) = factors.iter().position(|(f, _)| f == factor) {
                let e = ops::add2(&factors[pos].1, exp);
                if e.is_zero() {
                    factors.remove(pos);
                    continue;
                }
                factors[pos].1 = e;
            } else {
                factors.push((factor.clone(), exp.clone()));
            }
        }
        Factors::from_map(factors)
    }

    pub fn normal(&self, other: &Factors) -> (Factors, Factors) {
        let mut sf = self.clone();
        let mut of = other.clone();
        for (factor, self_exp) in &self.factors {
            let Some(other_exp) = other.get(factor).cloned() else { continue };
            let exp = ops::sub(self_exp, &other_exp);
            if exp.is_zero() {
                sf.remove(factor);
                of.remove(factor);
            } else if isnumber(&exp) {
                if exp.as_number().is_some_and(Number::is_positive) || exp.is_constant(Constant::Infinity) {
                    sf.set(factor.clone(), exp);
                    of.remove(factor);
                } else {
                    sf.remove(factor);
                    of.set(factor.clone(), ops::neg(&exp));
                }
            } else if let Some(r) = extract_additively(self_exp, &other_exp) {
                if !r.is_zero() {
                    sf.set(factor.clone(), r);
                    of.remove(factor);
                } else {
                    sf.remove(factor);
                    of.remove(factor);
                }
            } else {
                let mut other_exp = other_exp;
                let (sc, sa) = super::as_coeff_add(self_exp, false);
                if !sc.is_zero() {
                    let (oc, oa) = super::as_coeff_add(&other_exp, false);
                    let diff = ops::sub(&sc, &oc);
                    let cur = sf.get(factor).cloned().unwrap();
                    match diff.as_number().map(|d| d.cmp_value(&Number::ZERO)) {
                        Some(Some(std::cmp::Ordering::Greater)) => {
                            sf.set(factor.clone(), ops::sub(&cur, &oc));
                            other_exp = oa;
                        }
                        Some(Some(std::cmp::Ordering::Less)) => {
                            sf.set(factor.clone(), ops::sub(&cur, &sc));
                            let oe = of.get(factor).cloned().unwrap();
                            of.set(factor.clone(), ops::sub(&oe, &sc));
                            other_exp = ops::sub(&oa, &diff);
                        }
                        _ => {
                            sf.set(factor.clone(), sa);
                            other_exp = oa;
                        }
                    }
                }
                if !other_exp.is_zero() {
                    of.set(factor.clone(), other_exp);
                } else {
                    of.remove(factor);
                }
            }
        }
        (Factors::from_map(sf.factors), Factors::from_map(of.factors))
    }

    pub fn div(&self, other: &Factors) -> (Factors, Factors) {
        let mut quo = self.clone();
        let mut rem = Factors::default();
        for (factor, exp) in &other.factors {
            let mut exp = exp.clone();
            if let Some(q) = quo.get(factor).cloned() {
                let d = ops::sub(&q, &exp);
                if isnumber(&d) {
                    let sign = d.as_number().and_then(|d| d.cmp_value(&Number::ZERO));
                    if matches!(sign, Some(std::cmp::Ordering::Less | std::cmp::Ordering::Equal)) {
                        quo.remove(factor);
                    }
                    if matches!(sign, Some(std::cmp::Ordering::Greater | std::cmp::Ordering::Equal)) {
                        if !d.is_zero() {
                            quo.set(factor.clone(), d);
                        }
                        continue;
                    }
                    exp = ops::neg(&d);
                } else {
                    if let Some(r) = extract_additively(&q, &exp) {
                        if !r.is_zero() { quo.set(factor.clone(), r) } else { quo.remove(factor) }
                    } else {
                        let mut other_exp = exp.clone();
                        let (sc, sa) = super::as_coeff_add(&q, false);
                        if !sc.is_zero() {
                            let (oc, oa) = super::as_coeff_add(&other_exp, false);
                            let diff = ops::sub(&sc, &oc);
                            match diff.as_number().and_then(|d| d.cmp_value(&Number::ZERO)) {
                                Some(std::cmp::Ordering::Greater) => {
                                    quo.set(factor.clone(), ops::sub(&q, &oc));
                                    other_exp = oa;
                                }
                                Some(std::cmp::Ordering::Less) => {
                                    quo.set(factor.clone(), ops::sub(&q, &sc));
                                    other_exp = ops::sub(&oa, &diff);
                                }
                                _ => {
                                    quo.set(factor.clone(), sa);
                                    other_exp = oa;
                                }
                            }
                        }
                        if !other_exp.is_zero() {
                            rem.set(factor.clone(), other_exp);
                        }
                    }
                    continue;
                }
            }
            rem.set(factor.clone(), exp);
        }
        (Factors::from_map(quo.factors), Factors::from_map(rem.factors))
    }

    pub fn quo(&self, other: &Factors) -> Factors {
        self.div(other).0
    }

    pub fn pow(&self, n: i64) -> Factors {
        if n == 0 {
            return Factors::default();
        }
        Factors::from_map(self.factors.iter().map(|(f, e)| (f.clone(), ops::mul(&[e.clone(), Expr::int(n)]))).collect())
    }

    pub fn gcd(&self, other: &Factors) -> Factors {
        if other.is_zero() {
            return self.clone();
        }
        let mut factors = Vec::new();
        for (factor, exp) in &self.factors {
            if let Some(oe) = other.get(factor) {
                match assumptions::is_negative(&ops::sub(exp, oe)) {
                    Some(true) => factors.push((factor.clone(), exp.clone())),
                    Some(false) => factors.push((factor.clone(), oe.clone())),
                    None => {}
                }
            }
        }
        Factors::from_map(factors)
    }

    pub fn lcm(&self, other: &Factors) -> Factors {
        if self.is_zero() || other.is_zero() {
            return Factors::zero();
        }
        let mut factors = self.factors.clone();
        for (factor, exp) in &other.factors {
            if let Some(pos) = factors.iter().position(|(f, _)| f == factor) {
                let cur = &factors[pos].1;
                // Python's max() keeps the first argument unless the second
                // compares greater.
                let bigger = match (cur.as_number(), exp.as_number()) {
                    (Some(a), Some(b)) => a.cmp_value(b) == Some(std::cmp::Ordering::Greater),
                    _ => assumptions::is_positive(&ops::sub(cur, exp)) == Some(true),
                };
                factors[pos].1 = if bigger { cur.clone() } else { exp.clone() };
            } else {
                factors.push((factor.clone(), exp.clone()));
            }
        }
        Factors::from_map(factors)
    }
}

pub(crate) fn rat_gcd(a: &Number, b: &Number) -> Number {
    match (a.numer_denom(), b.numer_denom()) {
        (Some((p1, q1)), Some((p2, q2))) => {
            use num_integer::Integer;
            if b.is_zero() {
                return b.clone();
            }
            Number::rational(p1.gcd(&p2), q1.lcm(&q2)).unwrap()
        }
        _ => Number::ONE,
    }
}

fn rat_lcm(a: &Number, b: &Number) -> Number {
    match (a.numer_denom(), b.numer_denom()) {
        (Some((p1, q1)), Some((p2, q2))) => {
            use num_integer::Integer;
            Number::rational(p1.lcm(&p2), q1.gcd(&q2)).unwrap()
        }
        _ => a.mul(b),
    }
}

/// SymPy's `Term`: `coeff * numer / denom`.
#[derive(Clone, Debug)]
pub struct Term {
    pub coeff: Number,
    pub numer: Factors,
    pub denom: Factors,
}

impl Term {
    pub fn from_expr(term: &Expr) -> Term {
        let (c, rest) = ops::as_coeff_mul(term);
        let mut coeff = c;
        let factors = if rest.is_one() { vec![] } else { mul_args(&rest) };
        let mut numer: Vec<(Expr, Expr)> = Vec::new();
        let mut denom: Vec<(Expr, Expr)> = Vec::new();
        for factor in factors {
            let (mut base, exp) = decompose_power(&factor);
            if base.is_add() {
                let (cont, b) = primitive(&base);
                base = b;
                coeff = coeff.mul(&ops::num_pow(&cont, &exp).as_number().cloned().unwrap_or(Number::ONE));
            }
            let target = if exp.is_positive() { &mut numer } else { &mut denom };
            let exp = if exp.is_positive() { exp } else { exp.neg() };
            if let Some(slot) = target.iter_mut().find(|(b, _)| *b == base) {
                slot.1 = ops::add2(&slot.1, &Expr::num(exp));
            } else {
                target.push((base, Expr::num(exp)));
            }
        }
        Term { coeff, numer: Factors::from_map(numer), denom: Factors::from_map(denom) }
    }

    pub fn as_expr(&self) -> Expr {
        ops::mul(&[Expr::num(self.coeff.clone()), ops::div(&self.numer.as_expr(), &self.denom.as_expr())])
    }

    pub fn mul(&self, other: &Term) -> Term {
        let coeff = self.coeff.mul(&other.coeff);
        let numer = self.numer.mul(&other.numer);
        let denom = self.denom.mul(&other.denom);
        let (numer, denom) = numer.normal(&denom);
        Term { coeff, numer, denom }
    }

    pub fn inv(&self) -> Term {
        Term { coeff: self.coeff.inv().unwrap_or(Number::ZERO), numer: self.denom.clone(), denom: self.numer.clone() }
    }

    pub fn quo(&self, other: &Term) -> Term {
        self.mul(&other.inv())
    }

    pub fn gcd(&self, other: &Term) -> Term {
        Term { coeff: rat_gcd(&self.coeff, &other.coeff), numer: self.numer.gcd(&other.numer), denom: self.denom.gcd(&other.denom) }
    }

    pub fn lcm(&self, other: &Term) -> Term {
        Term { coeff: rat_lcm(&self.coeff, &other.coeff), numer: self.numer.lcm(&other.numer), denom: self.denom.lcm(&other.denom) }
    }
}

/// `_gcd_terms(terms, isprimitive, fraction)`: `(cont, numer, denom)`.
pub fn gcd_terms_raw(terms: &[Expr], isprimitive: bool, fraction: bool) -> (Expr, Expr, Expr) {
    let mut terms: Vec<Term> = terms.iter().filter(|t| !t.is_zero()).map(Term::from_expr).collect();
    if terms.is_empty() {
        return (Expr::zero(), Expr::zero(), Expr::one());
    }
    let (mut cont, mut numer, denom);
    if terms.len() == 1 {
        cont = Expr::num(terms[0].coeff.clone());
        numer = terms[0].numer.as_expr();
        denom = terms[0].denom.as_expr();
    } else {
        let mut c = terms[0].clone();
        for t in &terms[1..] {
            c = c.gcd(t);
        }
        for t in terms.iter_mut() {
            *t = t.quo(&c);
        }
        let numers: Vec<Expr>;
        let d;
        if fraction {
            let mut dd = terms[0].denom.clone();
            for t in &terms[1..] {
                dd = dd.lcm(&t.denom);
            }
            numers = terms
                .iter()
                .map(|t| ops::mul(&[Expr::num(t.coeff.clone()), t.numer.mul(&dd.quo(&t.denom)).as_expr()]))
                .collect();
            d = dd;
        } else {
            numers = terms.iter().map(Term::as_expr).collect();
            d = Factors::default();
        }
        cont = c.as_expr();
        numer = ops::add(&numers);
        denom = d.as_expr();
    }
    if !isprimitive && numer.is_add() {
        let (c2, n2) = primitive(&numer);
        cont = ops::mul(&[cont, Expr::num(c2)]);
        numer = n2;
    }
    (cont, numer, denom)
}

fn gcd_terms_add(terms: &[Expr], isprimitive: bool, clear: bool, fraction: bool) -> Expr {
    let (cont, mut numer, denom) = gcd_terms_raw(terms, isprimitive, fraction);
    let (coeff, factors) = as_coeff_mul(&cont, false);
    let mut coeff = coeff;
    if !clear {
        let (c, rest) = as_coeff_mul(&coeff, false);
        if !c.is_integer_number() && numer.is_add() {
            if let Some((n, d)) = c.as_number().and_then(Number::numer_denom) {
                let n = Expr::num(Number::from_bigint(n));
                let d = Expr::num(Number::from_bigint(d));
                let new = ops::div(&numer, &d);
                if add_args(&new).iter().any(|a| as_coeff_mul(a, false).0.is_integer_number()) {
                    numer = new;
                    coeff = ops::mul(&[n, rest]);
                }
            }
        }
    }
    keep_coeff(&coeff, &ops::div(&ops::mul(&[factors, numer]), &denom), clear, false)
}

/// `gcd_terms(terms, isprimitive, clear, fraction)` for a single expression.
pub fn gcd_terms(e: &Expr, isprimitive: bool, clear: bool, fraction: bool) -> Expr {
    if decompose::commutativity(e) != Some(true) {
        return e.clone();
    }
    match e.kind() {
        Kind::Add(..) => gcd_terms_add(&args(e), isprimitive, clear, fraction),
        _ if e.is_atom() => e.clone(),
        Kind::Mul(..) => {
            let (c, rest) = ops::as_coeff_mul(e);
            let parts: Vec<Expr> = mul_args(&rest).iter().map(|a| gcd_terms(a, isprimitive, clear, fraction)).collect();
            keep_coeff(&Expr::num(c), &ops::mul(&parts), clear, false)
        }
        _ => traverse::map_children(e, &mut |a| gcd_terms(a, isprimitive, clear, fraction)),
    }
}

/// `gcd_terms` applied to a list of terms (the Add-like input form).
pub fn gcd_terms_list(terms: &[Expr], isprimitive: bool, clear: bool, fraction: bool) -> Expr {
    gcd_terms_add(terms, isprimitive, clear, fraction)
}

fn is_negative_coeff(a: &Expr) -> bool {
    let (c, _) = as_coeff_mul(a, false);
    match c.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Const(Constant::NegInfinity) => true,
        _ => false,
    }
}

/// `factor_terms(expr, radical, clear, fraction, sign)`.
pub fn factor_terms(e: &Expr, radical: bool, clear: bool, fraction: bool, sign: bool) -> Expr {
    fn go(e: &Expr, radical: bool, clear: bool, fraction: bool, sign: bool) -> Expr {
        if e.is_atom() {
            return e.clone();
        }
        let rec = |a: &Expr| go(a, radical, clear, fraction, sign);
        match e.kind() {
            Kind::Add(..) | Kind::Mul(..) => {}
            Kind::Func(Func::Integral | Func::Sum, a) if !a.is_empty() => {
                return factor_sum_int(e, a, radical, clear, fraction, sign);
            }
            _ => return traverse::map_children(e, &mut |a| rec(a)),
        }
        let Ok((cont, p)) = decompose::as_content_primitive(e, radical, clear) else {
            return e.clone();
        };
        let mut cont = Expr::num(cont);
        let p = if p.is_add() {
            let mut list: Vec<Expr> = args(&p).iter().map(rec).collect();
            if list.iter().all(is_negative_coeff) {
                cont = ops::neg(&cont);
                list = list.iter().map(ops::neg).collect();
            }
            let p = add_from_args(&list);
            gcd_terms(&p, true, clear, fraction)
        } else if p.is_mul() {
            ops::mul(&args(&p).iter().map(rec).collect::<Vec<_>>())
        } else if p.is_atom() {
            p
        } else {
            traverse::map_children(&p, &mut |a| rec(a))
        };
        keep_coeff(&cont, &p, clear, sign)
    }
    go(e, radical, clear, fraction, sign)
}

fn factor_sum_int(e: &Expr, a: &[Expr], radical: bool, clear: bool, fraction: bool, sign: bool) -> Expr {
    let Kind::Func(head, _) = e.kind() else { return e.clone() };
    let result = &a[0];
    if result.is_zero() {
        return Expr::zero();
    }
    let limits = &a[1..];
    let wrt: Vec<Expr> = limits
        .iter()
        .map(|l| match l.kind() {
            Kind::Func(_, t) if !t.is_empty() && !l.is_symbol() => t[0].clone(),
            _ => l.clone(),
        })
        .collect();
    let f = factor_terms(result, radical, clear, fraction, sign);
    let (i, d) = decompose::as_independent(&f, &wrt, None);
    let build = |body: Expr| {
        let mut v = vec![body];
        v.extend(limits.iter().cloned());
        crate::functions::func(head.clone(), v)
    };
    if f.is_add() {
        ops::add(&[ops::mul(&[i, build(Expr::one())]), build(d)])
    } else {
        ops::mul(&[i, build(d)])
    }
}

/// `factor_nc(expr)` for commutative expressions.
pub fn factor_nc(e: &Expr) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    if decompose::commutativity(e) != Some(true) {
        return e.clone();
    }
    if !e.is_add() {
        return traverse::map_children(e, &mut factor_nc);
    }
    let expr = ops::add(&args(e).iter().map(|a| super::expand_hints(a, false, false, false, true, false, false)).collect::<Vec<_>>());
    let terms = add_args(&expr);
    let mut c = Expr::one();
    for (i, a) in terms.iter().enumerate() {
        if i == 0 {
            c = a.clone();
        } else {
            c = crate::poly::api::gcd(&c, a, &[]).unwrap_or_else(|_| Expr::one());
        }
    }
    let (mut cn, mut g) = (Expr::one(), Expr::one());
    let mut cc: Vec<Expr> = terms.clone();
    let hit = !c.is_one();
    if hit {
        let (k, rest) = as_coeff_mul(&c, false);
        cn = k;
        g = rest;
        if !g.is_one() {
            cc = cc.iter().map(|t| ops::div(t, &g)).collect();
        }
        cc = cc.iter().map(|t| ops::div(t, &cn)).collect();
    }
    let mid = if hit { ops::add(&cc) } else { expr.clone() };
    let factored = crate::poly::api::factor(&mid).unwrap_or_else(|_| mid.clone());
    let new_mid = super::powsimp::powsimp(&factored, &super::powsimp::PowsimpOptions::default());
    if new_mid.is_pow() {
        return keep_coeff(&cn, &ops::mul(&[g, new_mid]), true, false);
    }
    if new_mid.is_mul() {
        let ok = ops::mul(&[g.clone(), new_mid]);
        let pem = |x: &Expr| super::expand_hints(x, true, true, false, true, false, false);
        if pem(&ok) == pem(&ops::div(&expr, &cn)) {
            return keep_coeff(&cn, &ok, true, false);
        }
    }
    keep_coeff(&cn, &ops::mul(&[g, mid]), true, false)
}
