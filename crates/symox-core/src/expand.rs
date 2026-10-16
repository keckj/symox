//! Expansion of products and powers of sums.
//!
//! The heavy lifting is done on a sparse polynomial representation over
//! "generators" (symbols, functions, radicals, ...) with exponent vectors
//! stored inline, which makes `expand((x+y+z+1)**20)` orders of magnitude
//! faster than a naive term-by-term expression multiplication.

use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use crate::expr::{Expr, Func, Kind};
use crate::functions::exp;
use crate::number::Number;
use crate::ops::{add, add_from_terms, mul, mul_factor, pow};

type Mono = SmallVec<[i32; 8]>;

/// Default `expand()` (hints `mul`, `multinomial`, `power_exp`).
pub fn expand(e: &Expr) -> Expr {
    let mut cache = FxHashMap::default();
    ex(e, &mut cache)
}

fn ex(e: &Expr, cache: &mut FxHashMap<Expr, Expr>) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    if let Some(v) = cache.get(e) {
        return v.clone();
    }
    let r = match e.kind() {
        Kind::Add(c, terms) => {
            let mut v: Vec<Expr> = Vec::with_capacity(terms.len() + 1);
            let mut unchanged: Vec<(Expr, Number)> = Vec::new();
            for (t, k) in terms {
                let nt = ex(t, cache);
                if nt == *t {
                    unchanged.push((t.clone(), k.clone()));
                } else {
                    v.push(mul(&[Expr::num(k.clone()), nt]));
                }
            }
            if v.is_empty() {
                e.clone()
            } else {
                v.push(add_from_terms(c.clone(), unchanged));
                add(&v)
            }
        }
        Kind::Mul(..) => {
            let factors: Vec<Expr> = e.args().iter().map(|f| ex(f, cache)).collect();
            if factors.iter().any(|f| f.is_add()) {
                expand_product(&factors)
            } else {
                mul(&factors)
            }
        }
        Kind::Pow(b, x) => {
            let b2 = ex(b, cache);
            let x2 = ex(x, cache);
            expand_pow(&b2, &x2, cache)
        }
        Kind::Func(Func::Exp, args) => {
            let a = ex(&args[0], cache);
            match a.kind() {
                Kind::Add(..) => {
                    let parts: Vec<Expr> = a.args().iter().map(exp).collect();
                    mul(&parts)
                }
                _ => exp(&a),
            }
        }
        _ => crate::traverse::map_children(e, &mut |c| ex(c, cache)),
    };
    cache.insert(e.clone(), r.clone());
    r
}

fn expand_pow(b: &Expr, x: &Expr, cache: &mut FxHashMap<Expr, Expr>) -> Expr {
    // x**(a+b) -> x**a * x**b, only when the base is known to be non-zero
    // (SymPy's power_exp hint without assumptions)
    let nonzero_base = match b.kind() {
        Kind::Num(n) => !n.is_zero(),
        Kind::Const(c) => !c.is_infinite() && *c != crate::expr::Constant::NaN,
        _ => false,
    };
    if let (Kind::Add(..), true) = (x.kind(), nonzero_base) {
        let parts: Vec<Expr> = x.args().iter().map(|t| pow(b, t)).collect();
        let m = mul(&parts);
        return if m.is_mul() {
            let factors: Vec<Expr> = m.args().iter().map(|f| ex(f, cache)).collect();
            if factors.iter().any(|f| f.is_add()) {
                expand_product(&factors)
            } else {
                mul(&factors)
            }
        } else if let Kind::Pow(b3, x3) = m.kind() {
            if b3.is_add() && x3.is_integer_number() {
                expand_pow(b3, x3, cache)
            } else {
                m
            }
        } else {
            m
        };
    }
    if let (Kind::Add(..), Some(n)) = (b.kind(), x.as_number()) {
        if let Some(k) = n.to_i64() {
            if k > 0 {
                return power_poly(b, k as u32);
            } else if k < 0 {
                return pow(&power_poly(b, (-k) as u32), &Expr::minus_one());
            }
        } else if n.is_rational() && !n.is_integer() {
            // (x+1)**(5/2) -> (x+1)**2 * sqrt(x+1), expanded
            let k = n.floor();
            if let Some(k) = k.to_i64() {
                if k >= 1 {
                    let frac = Expr::num(n.sub(&Number::int(k)));
                    let radical = pow(b, &frac);
                    return expand_product(&[power_poly(b, k as u32), radical]);
                }
            }
        }
    }
    let p = pow(b, x);
    // (a*b)**2 may become a product of expandable powers
    if let Kind::Mul(..) = p.kind() {
        if !matches!(b.kind(), Kind::Mul(..)) || !x.is_integer_number() {
            return p;
        }
        let factors: Vec<Expr> = p.args().iter().map(|f| ex(f, cache)).collect();
        if factors.iter().any(|f| f.is_add()) {
            return expand_product(&factors);
        }
        return mul(&factors);
    }
    p
}

// ----------------------------------------------------------------------
// sparse polynomial machinery
// ----------------------------------------------------------------------

struct Gens {
    list: Vec<Expr>,
    index: FxHashMap<Expr, usize>,
}

impl Gens {
    fn new() -> Self {
        Gens {
            list: Vec::new(),
            index: FxHashMap::default(),
        }
    }
    fn id(&mut self, g: Expr) -> usize {
        if let Some(i) = self.index.get(&g) {
            return *i;
        }
        let i = self.list.len();
        self.index.insert(g.clone(), i);
        self.list.push(g);
        i
    }
}

/// A term as (coefficient, [(gen id, exponent)]).
type RawTerm = (Number, SmallVec<[(usize, i32); 8]>);

fn gen_of(b: &Expr, x: &Expr, gens: &mut Gens) -> (usize, i32) {
    if let Some(n) = x.as_number() {
        if let Some(k) = n.to_i64() {
            if let Ok(k) = i32::try_from(k) {
                return (gens.id(b.clone()), k);
            }
        }
        if let Number::Rational(r) = n {
            if let (Ok(p), true) = (
                i32::try_from(num_traits::ToPrimitive::to_i64(r.numer()).unwrap_or(i64::MAX)),
                r.denom() < &num_bigint::BigInt::from(1 << 20),
            ) {
                let q = num_traits::ToPrimitive::to_i64(r.denom()).unwrap();
                let g = Expr::from_kind(Kind::Pow(b.clone(), Expr::rational(1, q)));
                return (gens.id(g), p);
            }
        }
    }
    (gens.id(mul_factor(b, x)), 1)
}

fn raw_terms(e: &Expr, gens: &mut Gens) -> Vec<RawTerm> {
    let one_term = |t: &Expr, k: &Number, gens: &mut Gens| -> RawTerm {
        let mut m: SmallVec<[(usize, i32); 8]> = SmallVec::new();
        let mut c = k.clone();
        match t.kind() {
            Kind::Num(n) => c = c.mul(n),
            Kind::Mul(mc, facs) => {
                c = c.mul(mc);
                for (b, x) in facs {
                    m.push(gen_of(b, x, gens));
                }
            }
            Kind::Pow(b, x) => m.push(gen_of(b, x, gens)),
            _ => m.push((gens.id(t.clone()), 1)),
        }
        (c, m)
    };
    match e.kind() {
        Kind::Add(c, terms) => {
            let mut v: Vec<RawTerm> = terms.iter().map(|(t, k)| one_term(t, k, gens)).collect();
            if !c.is_zero() {
                v.push((c.clone(), SmallVec::new()));
            }
            v
        }
        _ => vec![one_term(e, &Number::ONE, gens)],
    }
}

type Poly = FxHashMap<Mono, Number>;

fn to_poly(raw: Vec<RawTerm>, n: usize) -> Poly {
    let mut p: Poly = FxHashMap::with_capacity_and_hasher(raw.len(), Default::default());
    for (c, m) in raw {
        let mut mono: Mono = SmallVec::from_elem(0, n);
        for (g, k) in m {
            mono[g] += k;
        }
        match p.get_mut(&mono) {
            Some(v) => *v = v.add(&c),
            None => {
                p.insert(mono, c);
            }
        }
    }
    p.retain(|_, c| !c.is_zero());
    p
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    let (a, b) = if a.len() < b.len() { (b, a) } else { (a, b) };
    let mut out: Poly =
        FxHashMap::with_capacity_and_hasher(a.len().saturating_mul(b.len()).min(1 << 22), Default::default());
    let mut mono: Mono = SmallVec::new();
    for (mb, cb) in b {
        for (ma, ca) in a {
            mono.clear();
            mono.extend(ma.iter().zip(mb.iter()).map(|(x, y)| x + y));
            let c = ca.mul(cb);
            match out.get_mut(&mono) {
                Some(v) => *v = v.add(&c),
                None => {
                    out.insert(mono.clone(), c);
                }
            }
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

fn poly_pow(p: &Poly, n: u32) -> Poly {
    if n == 1 {
        return p.clone();
    }
    // repeated multiplication by the (small) base beats squaring for
    // sparse inputs
    let mut r = p.clone();
    for _ in 1..n {
        r = poly_mul(&r, p);
    }
    r
}

/// Generators that can be assembled into a `Mul` without re-canonicalising.
fn is_simple_gen(g: &Expr) -> bool {
    match g.kind() {
        Kind::Sym(_) | Kind::Add(..) => true,
        Kind::Func(f, _) => !matches!(f, Func::Exp),
        _ => false,
    }
}

fn from_poly(p: Poly, gens: &Gens) -> Expr {
    let n = gens.list.len();
    if gens.list.iter().all(is_simple_gen) {
        // gens in canonical order: monomial factors come out sorted
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&i, &j| gens.list[i].cmp(&gens.list[j]));
        let mut coeff = Number::ZERO;
        let mut extra: Vec<Expr> = Vec::new();
        let mut terms: Vec<(Expr, Number)> = Vec::with_capacity(p.len());
        for (m, c) in p {
            let mut facs: Vec<(Expr, Expr)> = Vec::new();
            for &i in &order {
                if m[i] != 0 {
                    facs.push((gens.list[i].clone(), Expr::int(m[i] as i64)));
                }
            }
            let t = match facs.len() {
                0 => {
                    coeff = coeff.add(&c);
                    continue;
                }
                1 if facs[0].0.is_add() && facs[0].1.is_one() => {
                    // a bare sum must be flattened into the result
                    extra.push(mul(&[Expr::num(c), facs[0].0.clone()]));
                    continue;
                }
                1 => mul_factor(&facs[0].0, &facs[0].1),
                _ => Expr::from_kind(Kind::Mul(Number::ONE, facs)),
            };
            terms.push((t, c));
        }
        let s = add_from_terms(coeff, terms);
        if extra.is_empty() {
            return s;
        }
        extra.push(s);
        return add(&extra);
    }
    let mut v: Vec<Expr> = Vec::with_capacity(p.len());
    for (m, c) in p {
        let mut f: Vec<Expr> = Vec::with_capacity(n + 1);
        f.push(Expr::num(c));
        for (i, k) in m.iter().enumerate() {
            if *k != 0 {
                f.push(pow(&gens.list[i], &Expr::int(*k as i64)));
            }
        }
        v.push(mul(&f));
    }
    // products of radicals may produce new sums (e.g. I*I = -1): flatten
    add(&v)
}

/// Expand a product of already expanded factors.
pub fn expand_product(factors: &[Expr]) -> Expr {
    let mut gens = Gens::new();
    let raws: Vec<Vec<RawTerm>> = factors.iter().map(|f| raw_terms(f, &mut gens)).collect();
    let n = gens.list.len();
    let mut acc: Option<Poly> = None;
    for r in raws {
        let p = to_poly(r, n);
        acc = Some(match acc {
            None => p,
            Some(a) => poly_mul(&a, &p),
        });
    }
    let r = from_poly(acc.unwrap_or_default(), &gens);
    if gens.list.iter().all(is_simple_gen) {
        r
    } else {
        // radical products may need a second pass ((1+sqrt(2))*(1-sqrt(2)))
        renormalize(&r)
    }
}

/// Expand `b**n` for a sum `b` and positive integer `n`.
pub fn power_poly(b: &Expr, n: u32) -> Expr {
    let mut gens = Gens::new();
    let raw = raw_terms(b, &mut gens);
    let p = to_poly(raw, gens.list.len());
    let r = from_poly(poly_pow(&p, n), &gens);
    if gens.list.iter().all(is_simple_gen) {
        r
    } else {
        renormalize(&r)
    }
}

fn renormalize(e: &Expr) -> Expr {
    // terms whose product produced a sum (e.g. sqrt(2)*(..)) are distributed
    if let Kind::Add(..) = e.kind() {
        let parts: Vec<Expr> = e
            .args()
            .iter()
            .map(|t| {
                if let Kind::Mul(_, facs) = t.kind() {
                    if facs.iter().any(|(b, x)| b.is_add() && x.is_one()) {
                        return expand(t);
                    }
                }
                t.clone()
            })
            .collect();
        return add(&parts);
    }
    e.clone()
}
