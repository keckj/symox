//! Port of `sympy.simplify.sqrtdenest`: denesting of square roots
//! (Borodin et al., Fateman, and the SymPy `_denester` algorithm).
use super::methods::count_ops;
use super::radsimp::{rad_rationalize, radsimp, split_surds};
use super::{add_args, args, map_args, mexpand, mul_args, mul_from_args, num_le, num_lt};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::sortkey::cmp_default;
use crate::ops;

fn half() -> Number {
    Number::rational(1.into(), 2.into()).unwrap()
}

fn sqrt(e: &Expr) -> Expr {
    ops::sqrt(e)
}

fn root4(e: &Expr) -> Expr {
    ops::pow(e, &Expr::rational(1, 4))
}

fn is_rational(e: &Expr) -> bool {
    e.as_number().is_some_and(|n| n.is_rational())
}

fn is_integer(e: &Expr) -> bool {
    e.as_number().is_some_and(|n| n.is_integer())
}

fn pow_parts(e: &Expr) -> Option<(&Expr, &Expr)> {
    match e.kind() {
        Kind::Pow(b, x) => Some((b, x)),
        _ => None,
    }
}

/// `is_sqrt(expr)`: a power with exponent `1/2` or `-1/2`.
pub fn is_sqrt(e: &Expr) -> bool {
    pow_parts(e).is_some_and(|(_, x)| x.as_number().is_some_and(|n| n.is_rational() && n.abs() == half()))
}

/// `sqrt_depth(p)`.
pub fn sqrt_depth(p: &Expr) -> usize {
    if p.is_constant(Constant::I) {
        return 1;
    }
    if p.is_atom() {
        return 0;
    }
    if p.is_add() || p.is_mul() {
        return args(p).iter().map(sqrt_depth).max().unwrap_or(0);
    }
    if is_sqrt(p) {
        return sqrt_depth(pow_parts(p).unwrap().0) + 1;
    }
    0
}

/// `is_algebraic(p)`.
pub fn is_algebraic(p: &Expr) -> bool {
    if is_rational(p) {
        return true;
    }
    if p.is_atom() {
        return false;
    }
    if is_sqrt(p) || pow_parts(p).is_some_and(|(_, x)| is_integer(x)) {
        return is_algebraic(pow_parts(p).unwrap().0);
    }
    if p.is_add() || p.is_mul() {
        return args(p).iter().all(is_algebraic);
    }
    false
}

fn subsets(n: usize) -> Vec<Vec<u8>> {
    match n {
        1 => vec![vec![1]],
        2 => vec![vec![1, 0], vec![0, 1], vec![1, 1]],
        3 => vec![vec![1, 0, 0], vec![0, 1, 0], vec![1, 1, 0], vec![0, 0, 1], vec![1, 0, 1], vec![0, 1, 1], vec![1, 1, 1]],
        _ => {
            let b = subsets(n - 1);
            let mut a: Vec<Vec<u8>> = b.iter().map(|x| [x.clone(), vec![0]].concat()).collect();
            let mut unit = vec![0; n - 1];
            unit.push(1);
            a.push(unit);
            a.extend(b.iter().map(|x| [x.clone(), vec![1]].concat()));
            a
        }
    }
}

/// `sqrtdenest(expr, max_iter)`.
pub fn sqrtdenest(e: &Expr, max_iter: usize) -> Expr {
    let mut expr = super::expand_mul(e);
    for _ in 0..max_iter {
        let z = sqrtdenest0(&expr);
        if expr == z {
            return expr;
        }
        expr = z;
    }
    expr
}

/// `_sqrt_match(p)`: `[a, b, r]` with `p = a + b*sqrt(r)`, or empty.
fn sqrt_match(p: &Expr) -> Vec<Expr> {
    let p = mexpand(p);
    if super::is_number_atom(&p) {
        return vec![p, Expr::zero(), Expr::zero()];
    }
    if p.is_add() {
        let mut pargs = args(&p);
        pargs.sort_by(cmp_default);
        if pargs.iter().all(|x| {
            let sq = ops::powi(x, 2);
            is_rational(&sq) && sq.as_number().is_some_and(|n| n.is_positive())
        }) {
            let (r, b, a) = split_surds(&p);
            return vec![a, b, r];
        }
        let mut v: Vec<(usize, Expr, usize)> = pargs.iter().enumerate().map(|(i, x)| (sqrt_depth(x), x.clone(), i)).collect();
        let nmax = v
            .iter()
            .max_by(|a, b| a.0.cmp(&b.0).then_with(|| cmp_default(&a.1, &b.1)).then_with(|| a.2.cmp(&b.2)))
            .cloned()
            .unwrap();
        if nmax.0 == 0 {
            return vec![];
        }
        let (depth, _, i) = nmax;
        let mut r = pargs.remove(i);
        v.remove(i);
        let mut b = Expr::one();
        if r.is_mul() {
            let mut bv = Vec::new();
            let mut rv = Vec::new();
            for x in args(&r) {
                if sqrt_depth(&x) < depth {
                    bv.push(x);
                } else {
                    rv.push(x);
                }
            }
            b = mul_from_args(&bv);
            r = mul_from_args(&rv);
        }
        let mut a1 = Vec::new();
        let mut b1 = vec![b];
        for x in &v {
            if x.0 < depth {
                a1.push(x.1.clone());
            } else {
                let x1 = &x.1;
                if *x1 == r {
                    b1.push(Expr::one());
                } else if x1.is_mul() {
                    let mut x1args = args(x1);
                    if let Some(pos) = x1args.iter().position(|t| *t == r) {
                        x1args.remove(pos);
                        b1.push(ops::mul(&x1args));
                    } else {
                        a1.push(x1.clone());
                    }
                } else {
                    a1.push(x1.clone());
                }
            }
        }
        return vec![ops::add(&a1), ops::add(&b1), ops::powi(&r, 2)];
    }
    let (b, r) = ops::as_coeff_mul(&p);
    if is_sqrt(&r) {
        return vec![Expr::zero(), Expr::num(b), ops::powi(&r, 2)];
    }
    vec![]
}

fn sqrtdenest0(expr: &Expr) -> Expr {
    if is_sqrt(expr) {
        let (n, d) = crate::simptools::as_numer_denom(expr);
        if d.is_one() {
            let mut expr = expr.clone();
            if let Some((base, _)) = pow_parts(&n) {
                if base.is_add() {
                    let mut a = args(base);
                    a.sort_by(cmp_default);
                    if a.len() > 2 && a.iter().all(|x| is_integer(&ops::powi(x, 2))) {
                        if let Some(r) = sqrtdenest_rec(&n) {
                            return r;
                        }
                    }
                    let inner: Vec<Expr> = a.iter().map(sqrtdenest0).collect();
                    expr = sqrt(&mexpand(&ops::add(&inner)));
                }
            }
            return sqrtdenest1(&expr, true);
        }
        return ops::div(&sqrtdenest0(&n), &sqrtdenest0(&d));
    }
    if expr.is_add() {
        let mut cs = Vec::new();
        let mut a = Vec::new();
        for arg in args(expr) {
            let (c, t) = ops::as_coeff_mul(&arg);
            cs.push(c);
            a.push(t);
        }
        if cs.iter().all(|c| c.is_rational()) && a.iter().all(is_sqrt) {
            return sqrt_ratcomb(cs.into_iter().map(Expr::num).collect(), a);
        }
    }
    if expr.is_atom() {
        return expr.clone();
    }
    map_args(expr, &mut sqrtdenest0)
}

/// `_sqrtdenest_rec(expr)`; `None` stands for `SqrtdenestStopIteration`.
fn sqrtdenest_rec(expr: &Expr) -> Option<Expr> {
    let Some((base, _)) = pow_parts(expr) else {
        return Some(sqrtdenest(expr, 3));
    };
    if num_lt(base, &Expr::zero()) {
        let inner = sqrtdenest_rec(&sqrt(&ops::neg(base)))?;
        return Some(ops::mul(&[Expr::i(), inner]));
    }
    let (g, a, b) = split_surds(base);
    let mut a = ops::mul(&[a, sqrt(&g)]);
    let mut b = b;
    if num_lt(&a, &b) {
        std::mem::swap(&mut a, &mut b);
    }
    let c2 = mexpand(&ops::sub(&ops::powi(&a, 2), &ops::powi(&b, 2)));
    let c = if args(&c2).len() > 2 && !c2.is_atom() {
        let (g, a1, b1) = split_surds(&c2);
        let mut a1 = ops::mul(&[a1, sqrt(&g)]);
        let mut b1 = b1;
        if num_lt(&a1, &b1) {
            std::mem::swap(&mut a1, &mut b1);
        }
        let c2_1 = mexpand(&ops::sub(&ops::powi(&a1, 2), &ops::powi(&b1, 2)));
        let c_1 = sqrtdenest_rec(&sqrt(&c2_1))?;
        let d_1 = sqrtdenest_rec(&sqrt(&ops::add(&[a1, c_1])))?;
        let (num, den) = rad_rationalize(&b1, &d_1);
        let s2 = sqrt(&Expr::int(2));
        mexpand(&ops::add(&[ops::div(&d_1, &s2), ops::div(&num, &ops::mul(&[den, s2.clone()]))]))
    } else {
        sqrtdenest1(&sqrt(&c2), true)
    };
    if sqrt_depth(&c) > 1 {
        return None;
    }
    let ac = ops::add(&[a, c]);
    let ac_len = if ac.is_atom() { 0 } else { args(&ac).len() };
    if ac_len >= 2 && count_ops(&ac) >= count_ops(base) {
        return None;
    }
    let d = sqrtdenest(&sqrt(&ac), 3);
    if sqrt_depth(&d) > 1 {
        return None;
    }
    let (num, den) = rad_rationalize(&b, &d);
    let s2 = sqrt(&Expr::int(2));
    let r = ops::add(&[ops::div(&d, &s2), ops::div(&num, &ops::mul(&[den, s2.clone()]))]);
    Some(mexpand(&radsimp(&r, true, 4)))
}

fn sqrtdenest1(expr: &Expr, denester: bool) -> Expr {
    if !is_sqrt(expr) {
        return expr.clone();
    }
    let a = pow_parts(expr).unwrap().0.clone();
    if a.is_atom() {
        return expr.clone();
    }
    let val = sqrt_match(&a);
    if val.is_empty() {
        return expr.clone();
    }
    let (a, b, r) = (val[0].clone(), val[1].clone(), val[2].clone());
    let d2 = mexpand(&ops::sub(&ops::powi(&a, 2), &ops::mul(&[ops::powi(&b, 2), r.clone()])));
    if is_rational(&d2) {
        if d2.as_number().unwrap().is_positive() {
            if let Some(z) = sqrt_numeric_denest(&a, &b, &r, &d2) {
                return z;
            }
        } else {
            let dr2 = mexpand(&ops::neg(&ops::mul(&[d2.clone(), r.clone()])));
            let dr = sqrt(&dr2);
            if is_rational(&dr) {
                if let Some(z) = sqrt_numeric_denest(&mexpand(&ops::mul(&[b.clone(), r.clone()])), &a, &r, &dr2) {
                    return ops::div(&z, &root4(&r));
                }
            }
        }
    } else if let Some(z) = sqrt_symbolic_denest(&a, &b, &r) {
        return z;
    }
    if !denester || !is_algebraic(expr) {
        return expr.clone();
    }
    if let Some(res) = sqrt_biquadratic_denest(expr, &a, &b, &r, &d2) {
        return res;
    }
    let mut av0 = vec![Some(a), Some(b), Some(r), Some(d2)];
    let (z, _) = denest_nested(&[radsimp(&ops::powi(expr, 2), true, 4)], &mut av0, 0, sqrt_depth(expr));
    if av0[1].is_none() {
        return expr.clone();
    }
    if let Some(z) = z {
        if sqrt_depth(&z) == sqrt_depth(expr) && count_ops(&z) > count_ops(expr) {
            return expr.clone();
        }
        return z;
    }
    expr.clone()
}

/// Coefficients of `e` as a polynomial in `y` (highest degree first).
fn poly_coeffs(e: &Expr, y: &Expr) -> Option<Vec<Expr>> {
    let e = crate::expand::expand(e);
    let mut by_degree: Vec<Vec<Expr>> = Vec::new();
    for term in add_args(&e) {
        let mut degree = 0usize;
        let mut rest = Vec::new();
        for f in mul_args(&term) {
            if f == *y {
                degree += 1;
            } else if let Some((b, x)) = pow_parts(&f).filter(|(b, _)| *b == y) {
                let _ = b;
                let n = x.as_number().filter(|n| n.is_integer() && n.is_positive())?;
                degree += usize::try_from(n.to_bigint()?).ok()?;
            } else if crate::traverse::has(&f, y) {
                return None;
            } else {
                rest.push(f);
            }
        }
        if by_degree.len() <= degree {
            by_degree.resize(degree + 1, Vec::new());
        }
        by_degree[degree].push(ops::mul(&rest));
    }
    let mut out: Vec<Expr> = by_degree.iter().map(|v| ops::add(v)).collect();
    while out.len() > 1 && out.last().is_some_and(|c| c.is_zero()) {
        out.pop();
    }
    out.reverse();
    Some(out)
}

fn sqrt_symbolic_denest(a: &Expr, b: &Expr, r: &Expr) -> Option<Expr> {
    let rval = sqrt_match(r);
    if rval.is_empty() {
        return None;
    }
    let (ra, rb, rr) = (&rval[0], &rval[1], &rval[2]);
    if rb.is_zero() {
        return None;
    }
    let y = Expr::dummy("y", crate::symbol::Assumptions::from_kwargs([("positive", true)]).ok()?);
    let sub = ops::div(&ops::sub(&ops::powi(&y, 2), ra), rb);
    let newa = crate::subs::subs(a, &[(sqrt(rr), sub)]);
    let coeffs = poly_coeffs(&newa, &y)?;
    if coeffs.len() != 3 {
        return None;
    }
    let (ca, cb, cc) = (&coeffs[0], ops::add(&[coeffs[1].clone(), b.clone()]), &coeffs[2]);
    let disc = mexpand(&ops::sub(&ops::powi(&cb, 2), &ops::mul(&[Expr::int(4), ca.clone(), cc.clone()])));
    if super::methods::equals(&disc, &Expr::zero()) != Some(true) {
        return None;
    }
    let inner = ops::add(&[sqrt(r), ops::div(&cb, &ops::mul(&[Expr::int(2), ca.clone()]))]);
    let mut z = sqrt(&ops::mul(&[ca.clone(), ops::powi(&inner, 2)]));
    if super::is_num(&z) {
        if let Ok((c, p)) = crate::decompose::as_content_primitive(&z, false, true) {
            z = mexpand(&mul_from_args(&[Expr::num(c), p]));
        }
    }
    Some(z)
}

fn sqrt_numeric_denest(a: &Expr, b: &Expr, r: &Expr, d2: &Expr) -> Option<Expr> {
    let d = sqrt(d2);
    let s = ops::add(&[a.clone(), d.clone()]);
    if sqrt_depth(&s) < sqrt_depth(r) + 1 || is_rational(&ops::powi(&s, 2)) {
        let sign = |e: &Expr| crate::functions::func(Func::Sign, vec![e.clone()]);
        let (mut s1, mut s2) = (sign(&s), sign(b));
        if s1.is_minus_one() && s2.is_minus_one() {
            s1 = Expr::one();
            s2 = Expr::one();
        }
        let res = ops::mul(&[
            ops::add(&[ops::mul(&[s1, sqrt(&ops::add(&[a.clone(), d.clone()]))]), ops::mul(&[s2, sqrt(&ops::sub(a, &d))])]),
            sqrt(&Expr::int(2)),
            Expr::rational(1, 2),
        ]);
        return Some(crate::expand::expand(&res));
    }
    None
}

fn sqrt_biquadratic_denest(expr: &Expr, a: &Expr, b: &Expr, r: &Expr, d2: &Expr) -> Option<Expr> {
    let base = pow_parts(expr)?.0;
    if num_le(r, &Expr::zero()) || num_lt(d2, &Expr::zero()) || b.is_zero() || sqrt_depth(base) < 2 {
        return None;
    }
    for x in [a, b, r] {
        if x.is_atom() {
            continue;
        }
        for y in args(x) {
            let y2 = ops::powi(&y, 2);
            if !(is_integer(&y2) && y2.as_number().is_some_and(|n| n.is_positive())) {
                return None;
            }
        }
    }
    let sqd = mexpand(&sqrtdenest(&sqrt(&radsimp(d2, true, 4)), 3));
    if sqrt_depth(&sqd) > 1 {
        return None;
    }
    let h = Expr::rational(1, 2);
    let x1 = ops::add(&[ops::mul(&[a.clone(), h.clone()]), ops::mul(&[sqd.clone(), h.clone()])]);
    let x2 = ops::sub(&ops::mul(&[a.clone(), h.clone()]), &ops::mul(&[sqd, h]));
    for x in [x1, x2] {
        let big_a = sqrtdenest(&sqrt(&x), 3);
        if sqrt_depth(&big_a) > 1 {
            continue;
        }
        let (bn, bd) = rad_rationalize(b, &mexpand(&ops::mul(&[Expr::int(2), big_a.clone()])));
        let big_b = ops::div(&bn, &bd);
        let mut z = ops::add(&[big_a, ops::mul(&[big_b, sqrt(r)])]);
        if num_lt(&z, &Expr::zero()) {
            z = ops::neg(&z);
        }
        return Some(mexpand(&z));
    }
    None
}

type DenestResult = (Option<Expr>, Option<Vec<u8>>);

fn denest_nested(nested: &[Expr], av0: &mut Vec<Option<Expr>>, h: usize, max_depth_level: usize) -> DenestResult {
    if h > max_depth_level || av0[1].is_none() {
        return (None, None);
    }
    let fallback = |nested: &[Expr]| (Some(sqrt(nested.last().unwrap())), Some(vec![0u8; nested.len()]));
    if av0[0].is_none() && nested.iter().all(super::is_number_atom) {
        for f in subsets(nested.len()) {
            let chosen: Vec<Expr> = (0..f.len()).filter(|&i| f[i] == 1).map(|i| nested[i].clone()).collect();
            let mut p = mexpand(&ops::mul(&chosen));
            if f.iter().filter(|&&v| v == 1).count() > 1 && *f.last().unwrap() == 1 {
                p = ops::neg(&p);
            }
            let sqp = sqrt(&p);
            if is_rational(&sqp) {
                return (Some(sqp), Some(f));
            }
        }
        return fallback(nested);
    }
    let r_big: Expr;
    let values: Vec<Vec<Expr>>;
    let nested2: Vec<Expr>;
    if av0[0].is_some() {
        values = vec![vec![av0[0].clone().unwrap(), av0[1].clone().unwrap()]];
        r_big = av0[2].clone().unwrap();
        nested2 = vec![av0[3].clone().unwrap(), r_big.clone()];
        av0[0] = None;
    } else {
        values = nested.iter().map(sqrt_match).filter(|v| !v.is_empty()).collect();
        let mut rr: Option<Expr> = None;
        for v in &values {
            if !v[2].is_zero() {
                match &rr {
                    Some(x) if *x != v[2] => {
                        av0[1] = None;
                        return (None, None);
                    }
                    Some(_) => {}
                    None => rr = Some(v[2].clone()),
                }
            }
        }
        let Some(x) = rr else {
            return fallback(nested);
        };
        r_big = x;
        let mut n2: Vec<Expr> = values
            .iter()
            .map(|v| ops::sub(&mexpand(&ops::powi(&v[0], 2)), &mexpand(&ops::mul(&[r_big.clone(), ops::powi(&v[1], 2)]))))
            .collect();
        n2.push(r_big.clone());
        nested2 = n2;
    }
    let (d, f) = denest_nested(&nested2, av0, h + 1, max_depth_level);
    let Some(f) = f.filter(|f| !f.is_empty()) else {
        return (None, None);
    };
    let Some(d) = d else {
        return (None, None);
    };
    let n = nested.len();
    if f.len() <= n {
        return (None, None);
    }
    if !(0..n).any(|i| f[i] == 1) {
        let v = values.last().unwrap();
        return (Some(sqrt(&ops::add(&[v[0].clone(), mexpand(&ops::mul(&[v[1].clone(), d]))]))), Some(f));
    }
    let chosen: Vec<Expr> = (0..n).filter(|&i| f[i] == 1).map(|i| nested[i].clone()).collect();
    let p = ops::mul(&chosen);
    let mut v = sqrt_match(&p);
    if v.is_empty() {
        return (None, None);
    }
    if f.iter().position(|&x| x == 1).is_some_and(|i| i < n - 1) && f[n - 1] == 1 {
        v[0] = ops::neg(&v[0]);
        v[1] = ops::neg(&v[1]);
    }
    let s2r = sqrt(&Expr::int(2));
    if f[n] == 0 {
        let vad = mexpand(&ops::add(&[v[0].clone(), d]));
        if num_le(&vad, &Expr::zero()) {
            return fallback(nested);
        }
        if !(sqrt_depth(&vad) <= sqrt_depth(&r_big) + 1 || super::is_number_atom(&ops::powi(&vad, 2))) {
            av0[1] = None;
            return (None, None);
        }
        let sqvad = sqrtdenest1(&sqrt(&vad), false);
        if sqrt_depth(&sqvad) > sqrt_depth(&r_big) + 1 {
            av0[1] = None;
            return (None, None);
        }
        let sqvad1 = radsimp(&ops::inv(&sqvad), true, 4);
        let res = mexpand(&ops::add(&[
            ops::div(&sqvad, &s2r),
            ops::div(&ops::mul(&[v[1].clone(), sqrt(&r_big), sqvad1]), &s2r),
        ]));
        return (Some(res), Some(f));
    }
    let s2 = ops::add(&[mexpand(&ops::mul(&[v[1].clone(), r_big.clone()])), d]);
    if num_le(&s2, &Expr::zero()) {
        return fallback(nested);
    }
    let fr = root4(&mexpand(&r_big));
    let s = sqrt(&s2);
    let res = mexpand(&ops::add(&[
        ops::div(&s, &ops::mul(&[s2r.clone(), fr.clone()])),
        ops::div(&ops::mul(&[v[0].clone(), fr]), &ops::mul(&[s2r, s])),
    ]));
    (Some(res), Some(f))
}

fn sqrt_ratcomb(mut cs: Vec<Expr>, mut a: Vec<Expr>) -> Expr {
    loop {
        let mut found = None;
        'outer: for i in 0..a.len().saturating_sub(1) {
            for j in i + 1..a.len() {
                let s1 = pow_parts(&a[i]).unwrap().0;
                let s2 = pow_parts(&a[j]).unwrap().0;
                let p = mexpand(&ops::mul(&[s1.clone(), s2.clone()]));
                let s = sqrtdenest(&sqrt(&p), 3);
                if s != sqrt(&p) {
                    found = Some((s, i, j));
                    break 'outer;
                }
            }
        }
        let Some((s, i1, i2)) = found else {
            let terms: Vec<Expr> = cs.iter().zip(&a).map(|(c, t)| ops::mul(&[c.clone(), t.clone()])).collect();
            return ops::add(&terms);
        };
        let c2 = cs.remove(i2);
        a.remove(i2);
        let base = pow_parts(&a[i1]).unwrap().0.clone();
        let extra = radsimp(&ops::div(&ops::mul(&[c2, s]), &base), true, 4);
        cs[i1] = ops::add(&[cs[i1].clone(), extra]);
    }
}
