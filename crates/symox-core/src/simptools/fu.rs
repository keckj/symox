//! Native port of `sympy.simplify.fu`: the Fu-Harrington-Wang rules
//! `TR0` ... `TR22`, their helpers and the greedy `fu` driver.
//!
//! The rules follow SymPy's code paths one by one (argument orders,
//! dictionary insertion orders and the places where SymPy keeps or
//! rebuilds expressions), so that the results match SymPy's output.
use super::exprtools::{factor_terms, gcd_terms, Factors};
use super::methods::count_ops;
use super::{add_args, args, as_coeff_mul, as_numer_denom, bottom_up, expand_mul, is_number_atom, mul_args, nodes, num_lt, ordered};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::func;
use crate::number::Number;
use crate::{assumptions as asm, ops};
use num_bigint::BigInt;
use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// small helpers

fn head(e: &Expr) -> Option<&Func> {
    match e.kind() {
        Kind::Func(f, _) => Some(f),
        _ => None,
    }
}

fn is_fn(e: &Expr, f: Func) -> bool {
    matches!(e.kind(), Kind::Func(g, a) if *g == f && a.len() == 1)
}

fn is_sc(e: &Expr) -> bool {
    is_fn(e, Func::Sin) || is_fn(e, Func::Cos)
}

/// `isinstance(e, TrigonometricFunction)`.
pub fn is_trig(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc, a) if a.len() == 1)
}

fn is_hyper(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth, a) if a.len() == 1)
}

fn arg0(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Func(_, a) => a[0].clone(),
        _ => e.clone(),
    }
}

fn base_exp(e: &Expr) -> Option<(Expr, Expr)> {
    match e.kind() {
        Kind::Pow(b, x) => Some((b.clone(), x.clone())),
        _ => None,
    }
}

fn f1(f: Func, a: &Expr) -> Expr {
    func(f, vec![a.clone()])
}
pub fn sin(a: &Expr) -> Expr {
    f1(Func::Sin, a)
}
pub fn cos(a: &Expr) -> Expr {
    f1(Func::Cos, a)
}
pub fn tan(a: &Expr) -> Expr {
    f1(Func::Tan, a)
}
pub fn cot(a: &Expr) -> Expr {
    f1(Func::Cot, a)
}
fn sec(a: &Expr) -> Expr {
    f1(Func::Sec, a)
}
fn csc(a: &Expr) -> Expr {
    f1(Func::Csc, a)
}

fn int(n: i64) -> Expr {
    Expr::int(n)
}
fn mul(v: &[Expr]) -> Expr {
    ops::mul(v)
}
fn mul2(a: &Expr, b: &Expr) -> Expr {
    ops::mul2(a, b)
}
fn add2(a: &Expr, b: &Expr) -> Expr {
    ops::add2(a, b)
}
fn div(a: &Expr, b: &Expr) -> Expr {
    ops::div(a, b)
}
fn pow(a: &Expr, b: &Expr) -> Expr {
    ops::pow(a, b)
}
fn half() -> Expr {
    Expr::rational(1, 2)
}
fn sqrt(n: i64) -> Expr {
    pow(&int(n), &half())
}

fn is_int_true(e: &Expr) -> bool {
    asm::is_integer(e) == Some(true)
}
fn is_pos_true(e: &Expr) -> bool {
    asm::is_positive(e) == Some(true)
}

/// Python truthiness of an expression or `None`.
fn truthy(e: &Option<Expr>) -> bool {
    e.as_ref().is_some_and(|e| !e.is_zero())
}

/// SymPy's `factor`: the polynomial generators are taken from
/// `expr.expand()`, so the arguments of functions and non-integer powers
/// come back expanded (`factor(sin(x*(y - z)))` is `sin(x*y - x*z)`).
pub(crate) fn factor(e: &Expr) -> Expr {
    fn gens(e: &Expr) -> Expr {
        match e.kind() {
            Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) => e.clone(),
            Kind::Add(..) => ops::add(&add_args(e).iter().map(gens).collect::<Vec<_>>()),
            Kind::Mul(..) => mul(&mul_args(e).iter().map(gens).collect::<Vec<_>>()),
            Kind::Pow(b, x) if x.is_integer_number() => pow(&gens(b), x),
            Kind::Func(..) | Kind::Pow(..) => crate::expand::expand(e),
            _ => e.clone(),
        }
    }
    let g = gens(e);
    crate::poly::api::factor(&g).unwrap_or(g)
}

/// `rv.has(f, ...)`.
fn has_any(e: &Expr, fs: &[Func]) -> bool {
    let mut hit = false;
    super::preorder(e, &mut |a| {
        if !hit {
            if let Some(f) = head(a) {
                if fs.contains(f) {
                    hit = true;
                }
            }
        }
    });
    hit
}

/// Numeric `<` between two numbers (False when undecided).
fn nlt(a: &Expr, b: &Expr) -> bool {
    match (a.as_number(), b.as_number()) {
        (Some(x), Some(y)) => x.cmp_value(y) == Some(Ordering::Less),
        _ => num_lt(a, b),
    }
}

/// `min(a, b)` for numbers, keeping `a` unless `b < a`.
fn nmin(a: &Expr, b: &Expr) -> Expr {
    if nlt(b, a) { b.clone() } else { a.clone() }
}

fn func_of(f: &Func, a: &Expr) -> Expr {
    func(f.clone(), vec![a.clone()])
}

/// Re-evaluate the numeric products held inside `e` (SymPy rebuilds
/// `n.func(*n.args)` for numeric `Mul` nodes).
fn reeval_numeric_muls(e: &Expr) -> Expr {
    if !e.has_held_arithmetic() {
        return e.clone();
    }
    let r = super::map_args(e, &mut |a| reeval_numeric_muls(a));
    if r.is_mul() && super::is_num(&r) {
        return mul(&args(&r));
    }
    r
}

// ---------------------------------------------------------------------------
// TR rules

/// `TR0`: `rv.normal().factor().expand()`.
pub fn tr0(rv: &Expr) -> Expr {
    crate::expand::expand(&factor(&super::normal(rv)))
}

/// `TR1`: sec, csc -> 1/cos, 1/sin.
pub fn tr1(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        if is_fn(rv, Func::Sec) {
            div(&Expr::one(), &cos(&arg0(rv)))
        } else if is_fn(rv, Func::Csc) {
            div(&Expr::one(), &sin(&arg0(rv)))
        } else {
            rv.clone()
        }
    })
}

/// `TR2`: tan, cot -> sin/cos, cos/sin.
pub fn tr2(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        if is_fn(rv, Func::Tan) {
            let a = arg0(rv);
            div(&sin(&a), &cos(&a))
        } else if is_fn(rv, Func::Cot) {
            let a = arg0(rv);
            div(&cos(&a), &sin(&a))
        } else {
            rv.clone()
        }
    })
}

type PDict = Vec<(Expr, Option<Expr>)>;

fn pd_get<'a>(d: &'a PDict, k: &Expr) -> Option<usize> {
    d.iter().position(|(b, _)| b == k)
}

/// Position of the key `f(x)` in `d` (SymPy builds `f(x, evaluate=False)`).
fn pd_find_fn(d: &PDict, f: Func, x: &Expr) -> Option<usize> {
    d.iter().position(|(b, _)| is_fn(b, f.clone()) && arg0(b) == *x)
}

fn tr2i_ok(k: &Expr, e: &Expr, half: bool) -> bool {
    (is_int_true(e) || is_pos_true(k))
        && (is_sc(k)
            || (half
                && k.is_add()
                && args(k).len() >= 2
                && args(k).iter().any(|a| mul_args(a).iter().any(|ai| is_fn(ai, Func::Cos)))))
}

fn tr2i_factorize(d: &mut PDict, ddone: &mut Vec<(Expr, Expr)>, half: bool) {
    let mut newk: Vec<(Expr, Expr)> = Vec::new();
    for (k, _) in d.iter() {
        if k.is_add() && args(k).len() > 1 {
            let knew = if half { factor(k) } else { factor_terms(k, false, false, false, true) };
            if knew != *k {
                newk.push((k.clone(), knew));
            }
        }
    }
    if newk.is_empty() {
        return;
    }
    for (k, _) in &newk {
        d.retain(|(b, _)| b != k);
    }
    let prod = mul(&newk.iter().map(|(_, k)| k.clone()).collect::<Vec<_>>());
    for (k, v2) in crate::decompose::as_powers_dict(&prod) {
        // defaultdict(int): reading a missing key inserts 0.
        let pos = match pd_get(d, &k) {
            Some(p) => p,
            None => {
                d.push((k.clone(), Some(Expr::zero())));
                d.len() - 1
            }
        };
        let cur = d[pos].1.clone().unwrap_or_else(Expr::zero);
        let v = add2(&cur, &v2);
        if tr2i_ok(&k, &v, half) {
            d[pos].1 = Some(v);
        } else {
            ddone.push((k, v));
        }
    }
}

fn tr2i_f(rv: &Expr, half: bool) -> Expr {
    if !rv.is_mul() {
        return rv.clone();
    }
    let (n, d) = as_numer_denom(rv);
    if n.is_atom() || d.is_atom() {
        return rv.clone();
    }
    let mut n: PDict = crate::decompose::as_powers_dict(&n).into_iter().map(|(k, v)| (k, Some(v))).collect();
    let mut ndone: Vec<(Expr, Expr)> = Vec::new();
    n.retain(|(k, v)| {
        let v = v.clone().unwrap();
        if tr2i_ok(k, &v, half) {
            true
        } else {
            ndone.push((k.clone(), v));
            false
        }
    });
    if n.is_empty() {
        return rv.clone();
    }
    let mut d: PDict = crate::decompose::as_powers_dict(&d).into_iter().map(|(k, v)| (k, Some(v))).collect();
    let mut ddone: Vec<(Expr, Expr)> = Vec::new();
    d.retain(|(k, v)| {
        let v = v.clone().unwrap();
        if tr2i_ok(k, &v, half) {
            true
        } else {
            ddone.push((k.clone(), v));
            false
        }
    });
    if d.is_empty() {
        return rv.clone();
    }
    tr2i_factorize(&mut n, &mut ndone, half);
    tr2i_factorize(&mut d, &mut ddone, half);
    let mut t: Vec<Expr> = Vec::new();
    for idx in 0..n.len() {
        let k = n[idx].0.clone();
        let nk = n[idx].1.clone();
        if is_fn(&k, Func::Sin) {
            let x = arg0(&k);
            if let Some(p) = pd_find_fn(&d, Func::Cos, &x).filter(|p| d[*p].1 == nk && nk.is_some()) {
                t.push(pow(&tan(&x), nk.as_ref().unwrap()));
                n[idx].1 = None;
                d[p].1 = None;
            } else if half {
                let a1 = add2(&Expr::one(), &cos(&x));
                let found = if is_fn(&cos(&x), Func::Cos) && arg0(&cos(&x)) == x { pd_get(&d, &a1) } else { None };
                if let Some(p) = found.filter(|p| d[*p].1 == nk && nk.is_some()) {
                    t.push(pow(&tan(&div(&x, &int(2))), nk.as_ref().unwrap()));
                    n[idx].1 = None;
                    d[p].1 = None;
                }
            }
        } else if is_fn(&k, Func::Cos) {
            let x = arg0(&k);
            if let Some(p) = pd_find_fn(&d, Func::Sin, &x).filter(|p| d[*p].1 == nk && nk.is_some()) {
                t.push(pow(&tan(&x), &ops::neg(nk.as_ref().unwrap())));
                n[idx].1 = None;
                d[p].1 = None;
            }
        } else if half && k.is_add() {
            let ka = args(&k);
            if ka.len() >= 2 && ka[0].is_one() && is_fn(&ka[1], Func::Cos) {
                let x = arg0(&ka[1]);
                let s = sin(&x);
                if let Some(p) = pd_find_fn(&d, Func::Sin, &x) {
                    let dv = d[p].1.clone();
                    if dv == nk && nk.is_some() && (is_int_true(dv.as_ref().unwrap()) || is_pos_true(&s)) {
                        t.push(pow(&tan(&div(&x, &int(2))), &ops::neg(nk.as_ref().unwrap())));
                        n[idx].1 = None;
                        d[p].1 = None;
                    }
                }
            }
        }
    }
    if t.is_empty() {
        return rv.clone();
    }
    let mut num = t;
    num.extend(n.iter().filter(|(_, e)| truthy(e)).map(|(b, e)| pow(b, e.as_ref().unwrap())));
    let den: Vec<Expr> = d.iter().filter(|(_, e)| truthy(e)).map(|(b, e)| pow(b, e.as_ref().unwrap())).collect();
    let r = div(&mul(&num), &mul(&den));
    let nd: Vec<Expr> = ndone.iter().map(|(b, e)| pow(b, e)).collect();
    let dd: Vec<Expr> = ddone.iter().map(|(b, e)| pow(b, e)).collect();
    mul2(&r, &div(&mul(&nd), &mul(&dd)))
}

/// `TR2i(rv, half)`: sin/cos ratios -> tan.
pub fn tr2i(rv: &Expr, half: bool) -> Expr {
    bottom_up(rv, &mut |rv| tr2i_f(rv, half))
}

fn trig_cofunction(f: &Func) -> Func {
    match f {
        Func::Cos => Func::Sin,
        Func::Sin => Func::Cos,
        Func::Tan => Func::Cot,
        Func::Cot => Func::Tan,
        Func::Sec => Func::Csc,
        _ => Func::Sec,
    }
}

/// `TR3`: induced formulas, sin(-a) = -sin(a) and the reduction to
/// arguments in `[0, pi/4]`.
pub fn tr3(rv: &Expr) -> Expr {
    let rv = reeval_numeric_muls(rv);
    bottom_up(&rv, &mut |rv| {
        if !is_trig(rv) {
            return rv.clone();
        }
        let f = head(rv).unwrap().clone();
        let r = func_of(&f, &super::misc::signsimp(&arg0(rv), true));
        if !is_trig(&r) {
            return r;
        }
        let a = arg0(&r);
        let pi = Expr::pi();
        let p1 = asm::is_positive(&ops::sub(&a, &div(&pi, &int(4))));
        let p2 = asm::is_positive(&ops::sub(&div(&pi, &int(2)), &a));
        if p1 == Some(true) && p2 == Some(true) {
            let g = trig_cofunction(head(&r).unwrap());
            return func_of(&g, &ops::sub(&div(&pi, &int(2)), &a));
        }
        r
    })
}

/// `TR4`: trig functions of the special angles are evaluated.
pub fn tr4(rv: &Expr) -> Expr {
    if !rv.has_held_arithmetic() {
        return rv.clone();
    }
    bottom_up(rv, &mut |x| {
        if is_trig(x) {
            let a = arg0(x);
            let r = div(&a, &Expr::pi());
            if let Some(q) = r.as_number().filter(|n| n.is_rational()).map(|n| n.denom()) {
                if [1, 2, 3, 4, 6].iter().any(|d| q == Number::int(*d)) {
                    let a = if a.is_atom() { a } else { super::rebuild(&a, &args(&a)) };
                    return func_of(head(x).unwrap(), &a);
                }
            }
        }
        x.clone()
    })
}

/// `_TR56(rv, f, g, h, max, pow)`.
pub fn tr56(rv: &Expr, f: Func, g: Func, h: &dyn Fn(&Expr) -> Expr, max: &Expr, powflag: bool) -> Expr {
    bottom_up(rv, &mut |rv| tr56_f(rv, &f, &g, h, max, powflag))
}

fn tr56_f(rv: &Expr, f: &Func, g: &Func, h: &dyn Fn(&Expr) -> Expr, max: &Expr, powflag: bool) -> Expr {
    let Some((b, e)) = base_exp(rv) else { return rv.clone() };
    if !is_fn(&b, f.clone()) {
        return rv.clone();
    }
    if asm::is_real(&e) != Some(true) {
        return rv.clone();
    }
    if asm::is_negative(&e) == Some(true) {
        return rv.clone();
    }
    if asm::is_positive(&ops::sub(&e, max)) == Some(true) {
        return rv.clone();
    }
    if e.is_one() {
        return rv.clone();
    }
    let x = arg0(&b);
    let hg2 = || h(&pow(&func_of(g, &x), &int(2)));
    if e == int(2) {
        return hg2();
    }
    if let Some(n) = e.as_number().filter(|n| n.is_integer()).and_then(|n| n.to_bigint()) {
        let two = BigInt::from(2);
        let half_e = n.clone() / &two;
        if n.clone() % &two == BigInt::from(1) {
            return mul2(&func_of(f, &x), &pow(&hg2(), &Expr::num(Number::from_bigint(half_e))));
        }
        let ee = if n == BigInt::from(4) {
            int(2)
        } else if !powflag {
            Expr::num(Number::from_bigint(half_e))
        } else {
            if crate::ntheory::perfect_power(&n).is_none() {
                return rv.clone();
            }
            Expr::num(Number::from_bigint(half_e))
        };
        return pow(&hg2(), &ee);
    }
    if e.as_number().is_some() {
        // non-integral rational exponents: `rv.exp % 2` is truthy
        return rv.clone();
    }
    // symbolic exponents: Mod(e, 2) decides
    let floor_half = || func(Func::Floor, vec![div(&e, &int(2))]);
    if asm::is_odd(&e) == Some(true) {
        return mul2(&func_of(f, &x), &pow(&hg2(), &floor_half()));
    }
    if asm::is_even(&e) == Some(true) && !powflag {
        return pow(&hg2(), &floor_half());
    }
    rv.clone()
}

fn one_minus(x: &Expr) -> Expr {
    ops::sub(&Expr::one(), x)
}
fn one_plus(x: &Expr) -> Expr {
    add2(&Expr::one(), x)
}
fn minus_one(x: &Expr) -> Expr {
    ops::sub(x, &Expr::one())
}

/// `TR5`: sin**2 -> 1 - cos**2.
pub fn tr5(rv: &Expr, max: &Expr, powflag: bool) -> Expr {
    tr56(rv, Func::Sin, Func::Cos, &one_minus, max, powflag)
}

/// `TR6`: cos**2 -> 1 - sin**2.
pub fn tr6(rv: &Expr, max: &Expr, powflag: bool) -> Expr {
    tr56(rv, Func::Cos, Func::Sin, &one_minus, max, powflag)
}

/// `TR7`: cos(x)**2 -> (1 + cos(2*x))/2.
pub fn tr7(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        if let Some((b, e)) = base_exp(rv) {
            if is_fn(&b, Func::Cos) && e == int(2) {
                return div(&one_plus(&cos(&mul2(&int(2), &arg0(&b)))), &int(2));
            }
        }
        rv.clone()
    })
}

/// `TR8(rv, first)`: products of sin and cos -> sums.
pub fn tr8(rv: &Expr, first: bool) -> Expr {
    bottom_up(rv, &mut |rv| tr8_f(rv, first))
}

fn tr8_f(rv: &Expr, first: bool) -> Expr {
    let pow_ok = match base_exp(rv) {
        Some((b, e)) => is_sc(&b) && (is_int_true(&e) || is_pos_true(&b)),
        None => false,
    };
    if !(rv.is_mul() || pow_ok) {
        return rv.clone();
    }
    if first {
        let (n, d) = as_numer_denom(rv);
        let (n, d) = (expand_mul(&n), expand_mul(&d));
        let newn = tr8(&n, false);
        let newd = tr8(&d, false);
        let mut rv = rv.clone();
        if newn != n || newd != d {
            rv = gcd_terms(&div(&newn, &newd), false, true, true);
            if rv.is_mul() {
                let a = args(&rv);
                if a.len() == 2 && a[0].is_rational_number() && a[1].is_add() {
                    let (c, m) = as_coeff_mul(&rv, false);
                    rv = mul2(&c, &m);
                }
            }
        }
        return rv;
    }
    let mut c: Vec<Expr> = Vec::new();
    let mut s: Vec<Expr> = Vec::new();
    let mut other: Vec<Expr> = Vec::new();
    for a in mul_args(rv) {
        if is_fn(&a, Func::Cos) {
            c.push(arg0(&a));
        } else if is_fn(&a, Func::Sin) {
            s.push(arg0(&a));
        } else if let Some((b, e)) = base_exp(&a).filter(|(b, e)| is_sc(b) && e.is_integer_number() && e.as_number().unwrap().is_positive()) {
            let k = e.as_i64().unwrap_or(0).max(0) as usize;
            let target = if is_fn(&b, Func::Cos) { &mut c } else { &mut s };
            for _ in 0..k {
                target.push(arg0(&b));
            }
        } else {
            other.push(a);
        }
    }
    if !((!c.is_empty() && !s.is_empty()) || c.len() > 1 || s.len() > 1) {
        return rv.clone();
    }
    let mut args_ = other;
    let n = c.len().min(s.len());
    for _ in 0..n {
        let a1 = s.pop().unwrap();
        let a2 = c.pop().unwrap();
        args_.push(div(&add2(&sin(&add2(&a1, &a2)), &sin(&ops::sub(&a1, &a2))), &int(2)));
    }
    while c.len() > 1 {
        let a1 = c.pop().unwrap();
        let a2 = c.pop().unwrap();
        args_.push(div(&add2(&cos(&add2(&a1, &a2)), &cos(&ops::sub(&a1, &a2))), &int(2)));
    }
    if let Some(a) = c.pop() {
        args_.push(cos(&a));
    }
    while s.len() > 1 {
        let a1 = s.pop().unwrap();
        let a2 = s.pop().unwrap();
        args_.push(div(&add2(&ops::neg(&cos(&add2(&a1, &a2))), &cos(&ops::sub(&a1, &a2))), &int(2)));
    }
    if let Some(a) = s.pop() {
        args_.push(sin(&a));
    }
    tr8(&expand_mul(&mul(&args_)), true)
}

/// Shared pairwise driver of `TR9` and `TR10i` (`do` in SymPy).
fn pairwise_do(rv: &Expr, combine: &dyn Fn(&Expr, &Expr) -> Option<Expr>) -> Expr {
    if !rv.is_add() {
        return rv.clone();
    }
    let mut a: Vec<Option<Expr>> = ordered(args(rv)).into_iter().map(Some).collect();
    if a.len() != 2 {
        let mut hit = false;
        for i in 0..a.len() {
            let Some(ai) = a[i].clone() else { continue };
            for j in i + 1..a.len() {
                let Some(aj) = a[j].clone() else { continue };
                let was = add2(&ai, &aj);
                let new = pairwise_do(&was, combine);
                if new != was {
                    a[i] = Some(new);
                    a[j] = None;
                    hit = true;
                    break;
                }
            }
        }
        let mut rv = rv.clone();
        if hit {
            let keep: Vec<Expr> = a.into_iter().filter(|x| truthy(x)).map(Option::unwrap).collect();
            rv = ops::add(&keep);
            if rv.is_add() {
                rv = pairwise_do(&rv, combine);
            }
        }
        return rv;
    }
    let (x, y) = (a[0].clone().unwrap(), a[1].clone().unwrap());
    combine(&x, &y).unwrap_or_else(|| rv.clone())
}

fn tr9_combine(x: &Expr, y: &Expr) -> Option<Expr> {
    let (gcd, n1, n2, mut a, mut b, iscos) = trig_split(x, y, false)?;
    let two = int(2);
    let sum = div(&add2(&a, &b), &two);
    let dif = div(&ops::sub(&a, &b), &two);
    if iscos {
        if n1 == n2 {
            return Some(mul(&[gcd, int(n1), two.clone(), cos(&sum), cos(&dif)]));
        }
        if n1 < 0 {
            std::mem::swap(&mut a, &mut b);
        }
        let sum = div(&add2(&a, &b), &two);
        let dif = div(&ops::sub(&a, &b), &two);
        Some(mul(&[int(-2), gcd, sin(&sum), sin(&dif)]))
    } else {
        if n1 == n2 {
            return Some(mul(&[gcd, int(n1), two.clone(), sin(&sum), cos(&dif)]));
        }
        if n1 < 0 {
            std::mem::swap(&mut a, &mut b);
        }
        let sum = div(&add2(&a, &b), &two);
        let dif = div(&ops::sub(&a, &b), &two);
        Some(mul(&[two, gcd, cos(&sum), sin(&dif)]))
    }
}

/// `TR9`: sums of sin or cos -> products.
pub fn tr9(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        if !rv.is_add() {
            return rv.clone();
        }
        process_common_addends(rv, &mut |e| pairwise_do(e, &tr9_combine), None, true)
    })
}

/// `TR10(rv, first)`: sin/cos of sums -> expanded.
pub fn tr10(rv: &Expr, first: bool) -> Expr {
    bottom_up(rv, &mut |rv| {
        if !is_sc(rv) {
            return rv.clone();
        }
        let arg = arg0(rv);
        if !arg.is_add() {
            return rv.clone();
        }
        let mut a_ = if first { ordered(args(&arg)) } else { args(&arg) };
        let a = a_.pop().unwrap();
        let b = super::add_from_args(&a_);
        let is_sin = is_fn(rv, Func::Sin);
        if b.is_add() {
            if is_sin {
                add2(&mul2(&sin(&a), &tr10(&cos(&b), false)), &mul2(&cos(&a), &tr10(&sin(&b), false)))
            } else {
                ops::sub(&mul2(&cos(&a), &tr10(&cos(&b), false)), &mul2(&sin(&a), &tr10(&sin(&b), false)))
            }
        } else if is_sin {
            add2(&mul2(&sin(&a), &cos(&b)), &mul2(&cos(&a), &sin(&b)))
        } else {
            ops::sub(&mul2(&cos(&a), &cos(&b)), &mul2(&sin(&a), &sin(&b)))
        }
    })
}

fn tr10i_combine(x: &Expr, y: &Expr) -> Option<Expr> {
    let (gcd, n1, n2, a, b, same) = trig_split(x, y, true)?;
    let gcd = mul2(&int(n1), &gcd);
    Some(if same {
        if n1 == n2 { mul2(&gcd, &cos(&ops::sub(&a, &b))) } else { mul2(&gcd, &cos(&add2(&a, &b))) }
    } else if n1 == n2 {
        mul2(&gcd, &sin(&add2(&a, &b)))
    } else {
        mul2(&gcd, &sin(&ops::sub(&b, &a)))
    })
}

fn tr10i_do(e: &Expr) -> Expr {
    pairwise_do(e, &tr10i_combine)
}

fn tr10i_f(rv: &Expr) -> Expr {
    if !rv.is_add() {
        return rv.clone();
    }
    let key2 = |x: &Expr| ordered(crate::traverse::free_symbols(x));
    let mut rv = process_common_addends(rv, &mut tr10i_do, Some(&key2), true);
    while rv.is_add() {
        let mut byrad: Vec<(Expr, Vec<Option<Expr>>)> = Vec::new();
        for a in args(&rv) {
            let mut key = Expr::one();
            if a.is_mul() {
                for ai in args(&a) {
                    if let Some((b, e)) = base_exp(&ai) {
                        if e == half() && b.is_integer_number() {
                            key = ai.clone();
                            break;
                        }
                    }
                }
            }
            match byrad.iter_mut().find(|(k, _)| *k == key) {
                Some(slot) => slot.1.push(Some(a)),
                None => byrad.push((key, vec![Some(a)])),
            }
        }
        let mut newargs: Vec<Expr> = Vec::new();
        let keys: Vec<Expr> = byrad.iter().map(|(k, _)| k.clone()).collect();
        for (ia, ka) in keys.iter().enumerate() {
            for kb in [mul2(&sqrt(3), ka), div(&Expr::one(), &sqrt(3))] {
                let Some(ib) = keys.iter().position(|k| *k == kb) else { continue };
                for i in 0..byrad[ia].1.len() {
                    let Some(x) = byrad[ia].1[i].clone() else { continue };
                    for j in 0..byrad[ib].1.len() {
                        let Some(y) = byrad[ib].1[j].clone() else { continue };
                        let was = add2(&x, &y);
                        let new = tr10i_do(&was);
                        if new != was {
                            newargs.push(new);
                            byrad[ia].1[i] = None;
                            byrad[ib].1[j] = None;
                            break;
                        }
                    }
                }
            }
        }
        if !newargs.is_empty() {
            for (_, v) in &byrad {
                let keep: Vec<Expr> = v.iter().filter(|x| truthy(x)).map(|x| x.clone().unwrap()).collect();
                newargs.push(ops::add(&keep));
            }
            rv = ops::add(&newargs);
        } else {
            rv = tr10i_do(&rv);
            break;
        }
    }
    rv
}

/// `TR10i`: sums of products -> sin/cos of sums.
pub fn tr10i(rv: &Expr) -> Expr {
    bottom_up(rv, &mut tr10i_f)
}

/// `TR11(rv, base)`: double angles -> products.
pub fn tr11(rv: &Expr, base: Option<&Expr>) -> Expr {
    bottom_up(rv, &mut |rv| {
        if !is_sc(rv) {
            return rv.clone();
        }
        let is_cos = is_fn(rv, Func::Cos);
        if let Some(base) = base {
            let f = head(rv).unwrap().clone();
            let mut t = func_of(&f, &mul2(base, &int(2)));
            let mut co = Expr::one();
            if t.is_mul() {
                let (c, m) = as_coeff_mul(&t, false);
                co = c;
                t = m;
            }
            if !is_sc(&t) {
                return rv.clone();
            }
            if arg0(rv) == arg0(&t) {
                let c = cos(base);
                let s = sin(base);
                return if is_cos {
                    div(&ops::sub(&pow(&c, &int(2)), &pow(&s, &int(2))), &co)
                } else {
                    div(&mul(&[int(2), c, s]), &co)
                };
            }
            return rv.clone();
        }
        let a = arg0(rv);
        if is_number_atom(&a) {
            return rv.clone();
        }
        let (c, m) = as_coeff_mul(&a, true);
        let Some(cn) = c.as_number() else { return rv.clone() };
        let Some((p, q)) = cn.numer_denom() else { return rv.clone() };
        if p.clone() % BigInt::from(2) == BigInt::from(0) {
            let arg = div(&mul2(&Expr::num(Number::from_bigint(p / BigInt::from(2))), &m), &Expr::num(Number::from_bigint(q)));
            let c = tr11(&cos(&arg), None);
            let s = tr11(&sin(&arg), None);
            return if is_cos {
                ops::sub(&pow(&c, &int(2)), &pow(&s, &int(2)))
            } else {
                mul(&[int(2), s, c])
            };
        }
        rv.clone()
    })
}

/// `_TR11`: TR11 driven by the half angles present in a ratio.
pub fn tr11_ratio(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        if !is_expr_node(rv) {
            return rv.clone();
        }
        let sincos_args = |flat: &Expr| {
            let mut sins: Vec<Expr> = Vec::new();
            let mut coss: Vec<Expr> = Vec::new();
            for fi in mul_args(flat) {
                let (b, e) = crate::decompose::as_base_exp(&fi);
                if e.is_integer_number() && e.as_number().unwrap().is_positive() {
                    if is_fn(&b, Func::Cos) {
                        if !coss.contains(&arg0(&b)) {
                            coss.push(arg0(&b));
                        }
                    } else if is_fn(&b, Func::Sin) && !sins.contains(&arg0(&b)) {
                        sins.push(arg0(&b));
                    }
                }
            }
            (sins, coss)
        };
        let (n, d) = as_numer_denom(rv);
        let (nsin, ncos) = sincos_args(&n);
        let (dsin, dcos) = sincos_args(&d);
        let handle = |mut rv: Expr, nsin: &Vec<Expr>, dsin: &mut Vec<Expr>, dcos: &mut Vec<Expr>| {
            for narg in nsin {
                let h = div(narg, &int(2));
                let target = if dcos.contains(&h) {
                    &mut *dcos
                } else if dsin.contains(&h) {
                    &mut *dsin
                } else {
                    continue;
                };
                rv = tr11(&rv, Some(&h));
                target.retain(|a| *a != h);
            }
            rv
        };
        let (mut dsin, mut dcos) = (dsin, dcos);
        let (mut nsin2, mut ncos2) = (nsin.clone(), ncos);
        let rv = handle(rv.clone(), &nsin, &mut dsin, &mut dcos);
        handle(rv, &dsin.clone(), &mut nsin2, &mut ncos2)
    })
}

/// `TR12(rv, first)`: tan of sums -> expanded.
pub fn tr12(rv: &Expr, first: bool) -> Expr {
    bottom_up(rv, &mut |rv| {
        if !is_fn(rv, Func::Tan) {
            return rv.clone();
        }
        let arg = arg0(rv);
        if !arg.is_add() {
            return rv.clone();
        }
        let mut a_ = if first { ordered(args(&arg)) } else { args(&arg) };
        let a = a_.pop().unwrap();
        let b = super::add_from_args(&a_);
        let tb = if b.is_add() { tr12(&tan(&b), false) } else { tan(&b) };
        let ta = tan(&a);
        div(&add2(&ta, &tb), &ops::sub(&Expr::one(), &mul2(&ta, &tb)))
    })
}

fn tr12i_dok(di: &Expr) -> Option<(Expr, Expr)> {
    let (g, f, s) = as_f_sign_1(di)?;
    if s.is_minus_one() && f.is_mul() {
        let fa = args(&f);
        if fa.len() == 2 && fa.iter().all(|x| is_fn(x, Func::Tan)) {
            return Some((g, f));
        }
    }
    None
}

fn tr12i_nok(ni: &Expr) -> Option<(Expr, Expr)> {
    if ni.is_add() {
        let a = args(ni);
        if a.len() == 2 && is_fn(&a[0], Func::Tan) && is_fn(&a[1], Func::Tan) {
            return Some((a[0].clone(), a[1].clone()));
        }
    }
    None
}

fn tan_args_sum(t: &[Expr]) -> Expr {
    ops::add(&t.iter().map(arg0).collect::<Vec<_>>())
}

fn tr12i_f(rv: &Expr) -> Expr {
    if !(rv.is_add() || rv.is_mul() || rv.is_pow()) {
        return rv.clone();
    }
    let (n, d) = as_numer_denom(rv);
    if args(&d).is_empty() || args(&n).is_empty() {
        return rv.clone();
    }
    let mut dok: Vec<(Expr, Expr)> = Vec::new();
    let set = |dok: &mut Vec<(Expr, Expr)>, k: Expr, v: Expr| match dok.iter_mut().find(|(a, _)| *a == k) {
        Some(slot) => slot.1 = v,
        None => dok.push((k, v)),
    };
    let mut d_args = mul_args(&d);
    let mut i = 0;
    while i < d_args.len() {
        let di = d_args[i].clone();
        if let Some((g, t)) = tr12i_dok(&di) {
            set(&mut dok, tan_args_sum(&args(&t)), Expr::one());
            d_args[i] = g;
            i += 1;
            continue;
        }
        if di.is_add() {
            let f = factor(&di);
            if f.is_mul() {
                d_args.extend(args(&f));
                d_args[i] = Expr::one();
            }
        } else if let Some((b, e)) = base_exp(&di).filter(|(b, e)| is_int_true(e) || is_pos_true(b)) {
            if let Some((g, t)) = tr12i_dok(&b) {
                set(&mut dok, tan_args_sum(&args(&t)), e.clone());
                d_args[i] = pow(&g, &e);
            } else {
                let f = factor(&di);
                if f.is_mul() {
                    d_args.extend(args(&f));
                    d_args[i] = Expr::one();
                }
            }
        }
        i += 1;
    }
    if dok.is_empty() {
        return rv.clone();
    }
    let mut n_args = mul_args(&factor_terms(&n, false, false, false, true));
    let mut hit = false;
    let mut i = 0;
    while i < n_args.len() {
        let ni = n_args[i].clone();
        let m = match tr12i_nok(&ni) {
            Some(m) => {
                n_args[i] = Expr::one();
                m
            }
            None => {
                if let Some(m) = tr12i_nok(&ops::neg(&ni)) {
                    n_args[i] = Expr::minus_one();
                    m
                } else if ni.is_add() {
                    let f = factor(&ni);
                    if f.is_mul() {
                        n_args.extend(args(&f));
                        n_args[i] = Expr::one();
                    }
                    i += 1;
                    continue;
                } else if let Some((b, _)) = base_exp(&ni).filter(|(b, e)| is_int_true(e) || is_pos_true(b)) {
                    if let Some(m) = tr12i_nok(&b) {
                        n_args[i] = Expr::one();
                        m
                    } else {
                        let f = factor(&ni);
                        if f.is_mul() {
                            n_args.extend(args(&f));
                            n_args[i] = Expr::one();
                        }
                        i += 1;
                        continue;
                    }
                } else {
                    i += 1;
                    continue;
                }
            }
        };
        hit = true;
        let s = add2(&arg0(&m.0), &arg0(&m.1));
        let Some(pos) = dok.iter().position(|(k, _)| *k == s) else {
            // SymPy raises KeyError here
            return rv.clone();
        };
        let ed = dok[pos].1.clone();
        if let Some(newed) = super::extract_additively(&ed, &Expr::one()) {
            if !newed.is_zero() {
                dok[pos].1 = newed;
            } else {
                dok.remove(pos);
            }
        }
        n_args[i] = mul2(&n_args[i], &ops::neg(&tan(&s)));
        i += 1;
    }
    if !hit {
        return rv.clone();
    }
    let extra: Vec<Expr> = dok
        .iter()
        .map(|(s, e)| pow(&minus_one(&ops::add(&add_args(s).iter().map(tan).collect::<Vec<_>>())), e))
        .collect();
    div(&div(&mul(&n_args), &mul(&d_args)), &mul(&extra))
}

/// `TR12i`: (tan(a) + tan(b))/(1 - tan(a)*tan(b)) -> tan(a + b).
pub fn tr12i(rv: &Expr) -> Expr {
    bottom_up(rv, &mut tr12i_f)
}

/// `TR13`: products of tan or cot -> sums.
pub fn tr13(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        if !rv.is_mul() {
            return rv.clone();
        }
        let mut t = Vec::new();
        let mut c = Vec::new();
        let mut other = Vec::new();
        for a in mul_args(rv) {
            if is_fn(&a, Func::Tan) {
                t.push(arg0(&a));
            } else if is_fn(&a, Func::Cot) {
                c.push(arg0(&a));
            } else {
                other.push(a);
            }
        }
        if t.len() < 2 && c.len() < 2 {
            return rv.clone();
        }
        while t.len() > 1 {
            let t1 = t.pop().unwrap();
            let t2 = t.pop().unwrap();
            let s = tan(&add2(&t1, &t2));
            other.push(ops::sub(&Expr::one(), &add2(&div(&tan(&t1), &s), &div(&tan(&t2), &s))));
        }
        if let Some(x) = t.pop() {
            other.push(tan(&x));
        }
        while c.len() > 1 {
            let t1 = c.pop().unwrap();
            let t2 = c.pop().unwrap();
            let s = cot(&add2(&t1, &t2));
            other.push(ops::add(&[Expr::one(), mul2(&cot(&t1), &s), mul2(&cot(&t2), &s)]));
        }
        if let Some(x) = c.pop() {
            other.push(cot(&x));
        }
        mul(&other)
    })
}

fn morrie_f(rv: &Expr, first: bool) -> Expr {
    if !rv.is_mul() {
        return rv.clone();
    }
    if first {
        let (n, d) = as_numer_denom(rv);
        return div(&morrie_f(&n, false), &morrie_f(&d, false));
    }
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    let mut coss: Vec<(Expr, Expr)> = Vec::new();
    let mut other: Vec<Expr> = Vec::new();
    for c in args(rv) {
        let (b, e) = crate::decompose::as_base_exp(&c);
        if e.is_integer_number() && is_fn(&b, Func::Cos) {
            let (co, a) = as_coeff_mul(&arg0(&b), false);
            match groups.iter_mut().find(|(k, _)| *k == a) {
                Some(g) => g.1.push(co),
                None => groups.push((a, vec![co])),
            }
            match coss.iter_mut().find(|(k, _)| *k == arg0(&b)) {
                Some(slot) => slot.1 = e,
                None => coss.push((arg0(&b), e)),
            }
        } else {
            other.push(c);
        }
    }
    let get = |coss: &Vec<(Expr, Expr)>, k: &Expr| coss.iter().position(|(a, _)| a == k);
    let mut new: Vec<Expr> = Vec::new();
    for gi in 0..groups.len() {
        let a = groups[gi].0.clone();
        groups[gi].1.sort_by(|x, y| x.as_number().and_then(|x| y.as_number().and_then(|y| x.cmp_value(y))).unwrap_or(Ordering::Equal));
        while !groups[gi].1.is_empty() {
            let c = &groups[gi].1;
            let mut k = 0u32;
            let ci = c[0].clone();
            let mut cc = ci.clone();
            while c.contains(&cc) {
                k += 1;
                cc = mul2(&cc, &int(2));
            }
            if k > 1 {
                let p2k = pow(&int(2), &int(k as i64));
                let newarg = div(&div(&sin(&mul(&[p2k.clone(), ci.clone(), a.clone()])), &p2k), &sin(&mul2(&ci, &a)));
                let mut take: Option<Expr> = None;
                let mut ccs = Vec::new();
                let mut ok = true;
                for _ in 0..k {
                    cc = div(&cc, &int(2));
                    let key = mul2(&a, &cc);
                    ccs.push(cc.clone());
                    let Some(p) = get(&coss, &key) else {
                        ok = false;
                        break;
                    };
                    let v = coss[p].1.clone();
                    take = Some(match take {
                        Some(t) if !t.is_zero() => nmin(&v, &t),
                        _ => nmin(&v, &v),
                    });
                }
                if !ok {
                    return rv.clone();
                }
                let take = take.unwrap();
                for _ in 0..k {
                    let cc = ccs.pop().unwrap();
                    let key = mul2(&a, &cc);
                    let p = get(&coss, &key).unwrap();
                    coss[p].1 = ops::sub(&coss[p].1, &take);
                    if coss[p].1.is_zero() {
                        let c = &mut groups[gi].1;
                        if let Some(q) = c.iter().position(|x| *x == cc) {
                            c.remove(q);
                        }
                    }
                }
                new.push(pow(&newarg, &take));
            } else {
                let c0 = groups[gi].1.remove(0);
                let key = mul2(&c0, &a);
                let b = cos(&key);
                let e = get(&coss, &key).map(|p| coss[p].1.clone()).unwrap_or_else(Expr::one);
                other.push(pow(&b, &e));
            }
        }
    }
    if new.is_empty() {
        return rv.clone();
    }
    let mut all = new;
    all.extend(other);
    for (a, ks) in &groups {
        for k in ks {
            all.push(cos(&mul2(k, a)));
        }
    }
    mul(&all)
}

/// `TRmorrie`: cos(x)*cos(2*x)*...*cos(2**(k-1)*x) -> sin(2**k*x)/(2**k*sin(x)).
pub fn trmorrie(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| morrie_f(rv, true))
}

#[derive(Clone)]
struct P14 {
    g: Expr,
    isnum: bool,
    e: Expr,
    f: Expr,
    si: Expr,
    a: Expr,
}

fn p14_cmp(x: &P14, y: &P14) -> Ordering {
    let nx = 3 + nodes(&x.g) + nodes(&x.e) + nodes(&x.f) + nodes(&x.a);
    let ny = 3 + nodes(&y.g) + nodes(&y.e) + nodes(&y.f) + nodes(&y.a);
    let c = crate::sortkey::cmp_default;
    nx.cmp(&ny)
        .then_with(|| c(&x.g, &y.g))
        .then_with(|| x.isnum.cmp(&y.isnum))
        .then_with(|| c(&x.e, &y.e))
        .then_with(|| c(&x.f, &y.f))
        .then_with(|| c(&x.si, &y.si))
        .then_with(|| c(&x.a, &y.a))
}

/// `TR14(rv, first)`: factored powers of `1 - sin**2` -> cos**2.
pub fn tr14(rv: &Expr, first: bool) -> Expr {
    bottom_up(rv, &mut |rv| tr14_f(rv, first))
}

fn tr14_f(rv: &Expr, first: bool) -> Expr {
    if !rv.is_mul() {
        return rv.clone();
    }
    if first {
        let (n, d) = as_numer_denom(rv);
        if !d.is_one() {
            let newn = tr14(&n, false);
            let newd = tr14(&d, false);
            if newn != n || newd != d {
                return div(&newn, &newd);
            }
            return rv.clone();
        }
    }
    let mut other: Vec<Expr> = Vec::new();
    let mut process: Vec<P14> = Vec::new();
    for a in args(rv) {
        let (a, e) = if a.is_pow() {
            let (b, e) = crate::decompose::as_base_exp(&a);
            if !(is_int_true(&e) || is_pos_true(&b)) {
                other.push(a);
                continue;
            }
            (b, e)
        } else {
            (a, Expr::one())
        };
        let m = as_f_sign_1(&a).filter(|m| is_sc(&m.1));
        let Some((g, f, si)) = m else {
            other.push(if e.is_one() { a } else { pow(&a, &e) });
            continue;
        };
        process.push(P14 { g, isnum: is_number_atom(&e), e, f, si, a });
    }
    process.sort_by(p14_cmp);
    let nother = other.len();
    while !process.is_empty() {
        let a = process.remove(0);
        if let Some(b) = process.first().cloned() {
            let anum = is_number_atom(&a.e);
            let bnum = is_number_atom(&b.e);
            if anum && bnum {
                if a.f == b.f && a.si != b.si {
                    let b = process.remove(0);
                    let take = nmin(&a.e, &b.e);
                    if b.e != take {
                        let mut rem = b.clone();
                        rem.e = ops::sub(&rem.e, &take);
                        process.insert(0, rem);
                    } else if a.e != take {
                        let mut rem = a.clone();
                        rem.e = ops::sub(&rem.e, &take);
                        process.insert(0, rem);
                    }
                    let t = if is_fn(&a.f, Func::Cos) { sin(&arg0(&a.f)) } else { cos(&arg0(&a.f)) };
                    other.push(pow(&mul(&[int(-1), a.g.clone(), b.g.clone(), pow(&t, &int(2))]), &take));
                    continue;
                }
            } else if a.e == b.e && a.f == b.f && a.si != b.si {
                let b = process.remove(0);
                let take = a.e.clone();
                let t = if is_fn(&a.f, Func::Cos) { sin(&arg0(&a.f)) } else { cos(&arg0(&a.f)) };
                other.push(pow(&mul(&[int(-1), a.g.clone(), b.g.clone(), pow(&t, &int(2))]), &take));
                continue;
            }
        }
        other.push(pow(&a.a, &a.e));
    }
    if other.len() != nother {
        return mul(&other);
    }
    rv.clone()
}

fn odd_exp(e: &Expr) -> bool {
    match e.as_number() {
        Some(n) if n.is_integer() => n.to_bigint().is_some_and(|b| ((b % BigInt::from(2)) + BigInt::from(2)) % BigInt::from(2) == BigInt::from(1)),
        Some(_) => false,
        None => asm::is_odd(e) == Some(true),
    }
}

/// `TR15`: negative powers of sin -> cot.
pub fn tr15(rv: &Expr, max: &Expr, powflag: bool) -> Expr {
    bottom_up(rv, &mut |rv| tr1516_f(rv, Func::Sin, Func::Cot, max, powflag))
}

/// `TR16`: negative powers of cos -> tan.
pub fn tr16(rv: &Expr, max: &Expr, powflag: bool) -> Expr {
    bottom_up(rv, &mut |rv| tr1516_f(rv, Func::Cos, Func::Tan, max, powflag))
}

fn tr1516_f(rv: &Expr, f: Func, g: Func, max: &Expr, powflag: bool) -> Expr {
    let Some((b, e)) = base_exp(rv) else { return rv.clone() };
    if !is_fn(&b, f.clone()) {
        return rv.clone();
    }
    if odd_exp(&e) {
        // SymPy calls TR15 here for both rules
        return div(&tr15(&pow(&b, &add2(&e, &Expr::one())), &int(4), false), &b);
    }
    let ia = div(&Expr::one(), rv);
    let a = tr56(&ia, f, g, &one_plus, max, powflag);
    if a != ia {
        return a;
    }
    rv.clone()
}

/// `TR111`: negative powers of tan, sin, cos -> cot, csc, sec.
pub fn tr111(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        let Some((b, e)) = base_exp(rv) else { return rv.clone() };
        if !(is_pos_true(&b) || (is_int_true(&e) && asm::is_negative(&e) == Some(true))) {
            return rv.clone();
        }
        let ne = ops::neg(&e);
        if is_fn(&b, Func::Tan) {
            pow(&cot(&arg0(&b)), &ne)
        } else if is_fn(&b, Func::Sin) {
            pow(&csc(&arg0(&b)), &ne)
        } else if is_fn(&b, Func::Cos) {
            pow(&sec(&arg0(&b)), &ne)
        } else {
            rv.clone()
        }
    })
}

/// `TR22`: tan**2 -> sec**2 - 1, cot**2 -> csc**2 - 1.
pub fn tr22(rv: &Expr, max: &Expr, powflag: bool) -> Expr {
    bottom_up(rv, &mut |rv| {
        let Some((b, _)) = base_exp(rv) else { return rv.clone() };
        if !(is_fn(&b, Func::Cot) || is_fn(&b, Func::Tan)) {
            return rv.clone();
        }
        let r = tr56(rv, Func::Tan, Func::Sec, &minus_one, max, powflag);
        tr56(&r, Func::Cot, Func::Csc, &minus_one, max, powflag)
    })
}

fn binom(n: i64, k: i64) -> Expr {
    Expr::num(Number::from_bigint(crate::ntheory::binomial(&BigInt::from(n), &BigInt::from(k))))
}

/// `TRpower`: powers of sin and cos -> sums of multiple angles.
pub fn trpower(rv: &Expr) -> Expr {
    bottom_up(rv, &mut |rv| {
        let Some((b, n)) = base_exp(rv) else { return rv.clone() };
        if !is_sc(&b) {
            return rv.clone();
        }
        let x = arg0(&b);
        let Some(n) = n.as_i64().filter(|n| *n > 0 && n.is_integer_number_guard()) else { return rv.clone() };
        let is_cos = is_fn(&b, Func::Cos);
        let c = pow(&int(2), &int(1 - n));
        let mut out;
        if n % 2 == 1 {
            let terms: Vec<Expr> = (0..(n + 1) / 2)
                .map(|k| {
                    let m = mul2(&int(n - 2 * k), &x);
                    if is_cos {
                        mul2(&binom(n, k), &cos(&m))
                    } else {
                        mul(&[binom(n, k), int(if k % 2 == 0 { 1 } else { -1 }), sin(&m)])
                    }
                })
                .collect();
            out = if is_cos {
                mul2(&c, &ops::add(&terms))
            } else {
                mul(&[c, int(if ((n - 1) / 2) % 2 == 0 { 1 } else { -1 }), ops::add(&terms)])
            };
        } else {
            let terms: Vec<Expr> = (0..n / 2)
                .map(|k| {
                    let m = cos(&mul2(&int(n - 2 * k), &x));
                    if is_cos { mul2(&binom(n, k), &m) } else { mul(&[binom(n, k), int(if k % 2 == 0 { 1 } else { -1 }), m]) }
                })
                .collect();
            out = if is_cos {
                mul2(&c, &ops::add(&terms))
            } else {
                mul(&[c, int(if (n / 2) % 2 == 0 { 1 } else { -1 }), ops::add(&terms)])
            };
            out = add2(&out, &mul2(&pow(&int(2), &int(-n)), &binom(n, n / 2)));
        }
        out
    })
}

trait IntGuard {
    fn is_integer_number_guard(&self) -> bool;
}
impl IntGuard for i64 {
    fn is_integer_number_guard(&self) -> bool {
        true
    }
}

/// `L(rv)`: the number of trigonometric functions.
pub fn l_count(rv: &Expr) -> usize {
    let mut n = 0;
    super::preorder(rv, &mut |a| {
        if is_trig(a) {
            n += 1;
        }
    });
    n
}

/// `sincos_to_sum(expr)`.
pub fn sincos_to_sum(expr: &Expr) -> Expr {
    if !has_any(expr, &[Func::Cos, Func::Sin]) {
        return expr.clone();
    }
    tr8(&expand_mul(&trpower(expr)), true)
}

// ---------------------------------------------------------------------------
// helpers

/// `process_common_addends(rv, do, key2, key1)`.
pub fn process_common_addends(
    rv: &Expr,
    do_: &mut dyn FnMut(&Expr) -> Expr,
    key2: Option<&dyn Fn(&Expr) -> Vec<Expr>>,
    key1: bool,
) -> Expr {
    let mut absc: Vec<((Expr, Vec<Expr>), Vec<Expr>)> = Vec::new();
    for a in args(rv) {
        let (c, a) = if key1 {
            let (c, a) = as_coeff_mul(&a, false);
            if nlt(&c, &Expr::zero()) { (ops::neg(&c), ops::neg(&a)) } else { (c, a) }
        } else {
            (Expr::one(), a)
        };
        let k = (c, key2.map(|f| f(&a)).unwrap_or_default());
        match absc.iter_mut().find(|(kk, _)| *kk == k) {
            Some(slot) => slot.1.push(a),
            None => absc.push((k, vec![a])),
        }
    }
    let mut out = Vec::new();
    let mut hit = false;
    for ((c, _), v) in &absc {
        if v.len() > 1 {
            let e = super::add_from_args(v);
            let new = do_(&e);
            if new != e {
                hit = true;
                out.push(mul2(c, &new));
            } else {
                out.push(mul2(c, &ops::add(v)));
            }
        } else {
            out.push(mul2(c, &v[0]));
        }
    }
    if hit {
        return ops::add(&out);
    }
    rv.clone()
}

/// `pow_cos_sin(a, two)` of `trig_split`: `(co, c, s)`.
fn pow_cos_sin(a: &Expr, two: bool) -> Option<(Option<Expr>, Option<Expr>, Option<Expr>)> {
    let mut c: Option<Expr> = None;
    let mut s: Option<Expr> = None;
    let mut co = Expr::one();
    let is_sqrt = |x: &Expr| base_exp(x).is_some_and(|(_, e)| e == half());
    if a.is_mul() {
        let (k, a) = as_coeff_mul(a, false);
        co = k;
        if args(&a).len() > 2 || !two {
            return None;
        }
        let mut list = if a.is_mul() { args(&a) } else { vec![a.clone()] };
        let a = list.remove(0);
        if is_fn(&a, Func::Cos) {
            c = Some(a);
        } else if is_fn(&a, Func::Sin) {
            s = Some(a);
        } else if is_sqrt(&a) {
            co = mul2(&co, &a);
        } else {
            return None;
        }
        if let Some(b) = list.first() {
            if is_fn(b, Func::Cos) {
                if c.is_some() { s = Some(b.clone()) } else { c = Some(b.clone()) }
            } else if is_fn(b, Func::Sin) {
                if s.is_some() { c = Some(b.clone()) } else { s = Some(b.clone()) }
            } else if is_sqrt(b) {
                co = mul2(&co, b);
            } else {
                return None;
            }
        }
        return Some((if co.is_one() { None } else { Some(co) }, c, s));
    } else if is_fn(a, Func::Cos) {
        c = Some(a.clone());
    } else if is_fn(a, Func::Sin) {
        s = Some(a.clone());
    }
    if c.is_none() && s.is_none() {
        return None;
    }
    Some((if co.is_one() { None } else { Some(co) }, c, s))
}

/// `trig_split(a, b, two)`: `(gcd, n1, n2, a, b, iscos_or_same)`.
pub fn trig_split(a: &Expr, b: &Expr, two: bool) -> Option<(Expr, i64, i64, Expr, Expr, bool)> {
    let fa = Factors::from_expr(a);
    let fb = Factors::from_expr(b);
    let (mut ua, mut ub) = fa.normal(&fb);
    let mut gcd = fa.gcd(&fb).as_expr();
    let (mut n1, mut n2) = (1i64, 1i64);
    let m1 = Factors::from_expr(&Expr::minus_one());
    if ua.factors.iter().any(|(k, _)| k.is_minus_one()) {
        ua = ua.quo(&m1);
        n1 = -n1;
    } else if ub.factors.iter().any(|(k, _)| k.is_minus_one()) {
        ub = ub.quo(&m1);
        n2 = -n2;
    }
    let (a, b) = (ua.as_expr(), ub.as_expr());
    let (mut coa, mut ca, mut sa) = pow_cos_sin(&a, two)?;
    let (mut cob, mut cb, mut sb) = pow_cos_sin(&b, two)?;
    if (ca.is_none() && cb.is_some()) || ca.as_ref().is_some_and(|x| is_fn(x, Func::Sin)) {
        std::mem::swap(&mut coa, &mut cob);
        std::mem::swap(&mut ca, &mut cb);
        std::mem::swap(&mut sa, &mut sb);
        std::mem::swap(&mut n1, &mut n2);
    }
    if !two {
        let c = ca.or(sa)?;
        let s = cb.or(sb)?;
        if head(&c) != head(&s) {
            return None;
        }
        return Some((gcd, n1, n2, arg0(&c), arg0(&s), is_fn(&c, Func::Cos)));
    }
    if coa.is_none() && cob.is_none() {
        if let (Some(ca), Some(cb), Some(sa), Some(sb)) = (&ca, &cb, &sa, &sb) {
            if (head(ca) == head(sa)) != (head(cb) == head(sb)) {
                return None;
            }
            let set = [arg0(ca), arg0(sa)];
            if !set.contains(&arg0(cb)) || !set.contains(&arg0(sb)) {
                return None;
            }
            return Some((gcd, n1, n2, arg0(ca), arg0(sa), head(ca) == head(sa)));
        }
    }
    if (ca.is_some() && sa.is_some()) || (cb.is_some() && sb.is_some()) || (ca.is_none() && sa.is_none()) || (cb.is_none() && sb.is_none()) {
        return None;
    }
    let c = ca.or(sa).unwrap();
    let s = cb.or(sb).unwrap();
    if arg0(&c) != arg0(&s) {
        return None;
    }
    let coa = coa.unwrap_or_else(Expr::one);
    let cob = cob.unwrap_or_else(Expr::one);
    let pi = Expr::pi();
    if coa == cob {
        gcd = mul2(&gcd, &sqrt(2));
        return Some((gcd, n1, n2, arg0(&c), div(&pi, &int(4)), false));
    }
    let r = div(&coa, &cob);
    if r == sqrt(3) {
        gcd = mul(&[gcd, int(2), cob]);
        return Some((gcd, n1, n2, arg0(&c), div(&pi, &int(3)), false));
    }
    if r == div(&Expr::one(), &sqrt(3)) {
        gcd = mul(&[gcd, int(2), coa]);
        return Some((gcd, n1, n2, arg0(&c), div(&pi, &int(6)), false));
    }
    None
}

/// `as_f_sign_1(e)`: `e = g*(f + s)` with `s = +/-1`.
pub fn as_f_sign_1(e: &Expr) -> Option<(Expr, Expr, Expr)> {
    if !e.is_add() {
        return None;
    }
    let ea = args(e);
    if ea.len() != 2 {
        return None;
    }
    let (a, b) = (ea[0].clone(), ea[1].clone());
    if a.is_one() || a.is_minus_one() {
        let mut g = Expr::one();
        let (mut a, mut b) = (a, b);
        if b.is_mul() {
            let ba = args(&b);
            if is_number_atom(&ba[0]) && nlt(&ba[0], &Expr::zero()) {
                a = ops::neg(&a);
                b = ops::neg(&b);
                g = ops::neg(&g);
            }
        }
        return Some((g, b, a));
    }
    let fa = Factors::from_expr(&a);
    let fb = Factors::from_expr(&b);
    let (mut ua, mut ub) = fa.normal(&fb);
    let mut gcd = fa.gcd(&fb).as_expr();
    let m1 = Factors::from_expr(&Expr::minus_one());
    let (mut n1, mut n2);
    if ua.factors.iter().any(|(k, _)| k.is_minus_one()) {
        ua = ua.quo(&m1);
        n1 = -1;
        n2 = 1;
    } else if ub.factors.iter().any(|(k, _)| k.is_minus_one()) {
        ub = ub.quo(&m1);
        n1 = 1;
        n2 = -1;
    } else {
        n1 = 1;
        n2 = 1;
    }
    let (mut a, mut b) = (ua.as_expr(), ub.as_expr());
    if a.is_one() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut n1, &mut n2);
    }
    if n1 == -1 {
        gcd = ops::neg(&gcd);
        n2 = -n2;
    }
    if b.is_one() {
        return Some((gcd, a, int(n2)));
    }
    None
}

/// `_osborne(e, d)`: hyperbolic -> trigonometric.
pub fn osborne(e: &Expr, d: &Expr) -> Expr {
    bottom_up(e, &mut |rv| {
        if !is_hyper(rv) {
            return rv.clone();
        }
        let a = arg0(rv);
        let a = if !a.is_add() {
            mul2(&a, d)
        } else {
            super::add_from_args(&args(&a).iter().map(|i| mul2(i, d)).collect::<Vec<_>>())
        };
        let i = Expr::i();
        match head(rv).unwrap() {
            Func::Sinh => mul2(&i, &sin(&a)),
            Func::Cosh => cos(&a),
            Func::Tanh => mul2(&i, &tan(&a)),
            _ => div(&cot(&a), &i),
        }
    })
}

/// `_osbornei(e, d)`: trigonometric -> hyperbolic.
pub fn osbornei(e: &Expr, d: &Expr) -> Expr {
    bottom_up(e, &mut |rv| {
        if !is_trig(rv) {
            return rv.clone();
        }
        let (cst, x) = crate::decompose::as_independent(&arg0(rv), std::slice::from_ref(d), Some(true));
        let a = add2(&super::misc::xreplace(&x, &[(d.clone(), Expr::one())]), &mul2(&cst, &Expr::i()));
        let i = Expr::i();
        let h = |f: Func| f1(f, &a);
        match head(rv).unwrap() {
            Func::Sin => div(&h(Func::Sinh), &i),
            Func::Cos => h(Func::Cosh),
            Func::Tan => div(&h(Func::Tanh), &i),
            Func::Cot => mul2(&h(Func::Coth), &i),
            Func::Sec => div(&Expr::one(), &h(Func::Cosh)),
            _ => div(&i, &h(Func::Sinh)),
        }
    })
}

/// `hyper_as_trig(rv)`: the masked trigonometric form, the Osborne dummy
/// and the `(dummy, trig)` pairs needed by [`hyper_as_trig_back`].
pub fn hyper_as_trig(rv: &Expr) -> (Expr, Expr, Vec<(Expr, Expr)>) {
    let trigs = super::atoms_where(rv, &is_trig);
    let reps: Vec<(Expr, Expr)> = trigs.iter().map(|t| (t.clone(), Expr::dummy("_t", crate::symbol::Assumptions::NONE))).collect();
    let masked = super::misc::xreplace(rv, &reps);
    let back: Vec<(Expr, Expr)> = reps.into_iter().map(|(k, v)| (v, k)).collect();
    let d = Expr::dummy("_d", crate::symbol::Assumptions::NONE);
    (osborne(&masked, &d), d, back)
}

/// The inverse function returned by `hyper_as_trig`.
pub fn hyper_as_trig_back(x: &Expr, d: &Expr, back: &[(Expr, Expr)]) -> Expr {
    let e = super::misc::xreplace(&osbornei(x, d), back);
    crate::poly::api::collect(&super::misc::signsimp(&e, true), &[Expr::i()])
}

// ---------------------------------------------------------------------------
// the greedy driver

/// The objective comparison `measure(a) < measure(b)`.
pub type Less<'a> = &'a dyn Fn(&Expr, &Expr) -> bool;

/// SymPy's default measure `(L(x), x.count_ops())`.
pub fn default_less(a: &Expr, b: &Expr) -> bool {
    (l_count(a), count_ops(a)) < (l_count(b), count_ops(b))
}

/// `min(seq, key=measure)`: the first minimal element.
pub fn min_by_measure(seq: Vec<Expr>, less: Less) -> Expr {
    let mut it = seq.into_iter();
    let mut best = it.next().expect("non-empty");
    for c in it {
        if less(&c, &best) {
            best = c;
        }
    }
    best
}

#[derive(Clone, Copy)]
enum Rule {
    Id,
    T0,
    T3,
    T4,
    T5,
    T6,
    T7,
    T8,
    T9,
    T10,
    T10i,
    T11,
    T12,
    T13,
    Morrie,
}

enum Strat {
    R(Rule),
    Chain(Vec<Strat>),
    Min(Vec<Strat>),
}

fn apply_rule(r: Rule, e: &Expr) -> Expr {
    let four = int(4);
    match r {
        Rule::Id => e.clone(),
        Rule::T0 => tr0(e),
        Rule::T3 => tr3(e),
        Rule::T4 => tr4(e),
        Rule::T5 => tr5(e, &four, false),
        Rule::T6 => tr6(e, &four, false),
        Rule::T7 => tr7(e),
        Rule::T8 => tr8(e, true),
        Rule::T9 => tr9(e),
        Rule::T10 => tr10(e, true),
        Rule::T10i => tr10i(e),
        Rule::T11 => tr11(e, None),
        Rule::T12 => tr12(e, true),
        Rule::T13 => tr13(e),
        Rule::Morrie => trmorrie(e),
    }
}

fn run(s: &Strat, e: &Expr, less: Less) -> Expr {
    match s {
        Strat::R(r) => apply_rule(*r, e),
        Strat::Chain(v) => {
            let mut x = e.clone();
            for s in v {
                x = run(s, &x, less);
            }
            x
        }
        Strat::Min(v) => min_by_measure(v.iter().map(|s| run(s, e, less)).collect(), less),
    }
}

fn ch(v: Vec<Strat>) -> Strat {
    Strat::Chain(v)
}
fn r(x: Rule) -> Strat {
    Strat::R(x)
}

fn ctr1() -> Strat {
    use Rule::*;
    Strat::Min(vec![ch(vec![r(T5), r(T0)]), ch(vec![r(T6), r(T0)]), r(Id)])
}
fn ctr2() -> Strat {
    use Rule::*;
    ch(vec![r(T11), Strat::Min(vec![ch(vec![r(T5), r(T0)]), ch(vec![r(T6), r(T0)]), r(T0)])])
}
fn ctr3() -> Strat {
    use Rule::*;
    Strat::Min(vec![ch(vec![r(Morrie), r(T8), r(T0)]), ch(vec![r(Morrie), r(T8), r(T10i), r(T0)]), r(Id)])
}
fn ctr4() -> Strat {
    use Rule::*;
    Strat::Min(vec![ch(vec![r(T4), r(T10i)]), r(Id)])
}
fn rl1() -> Strat {
    use Rule::*;
    ch(vec![r(T4), r(T3), r(T4), r(T12), r(T4), r(T13), r(T4), r(T0)])
}
fn rl2() -> Strat {
    use Rule::*;
    Strat::Min(vec![
        ch(vec![r(T4), r(T3), r(T10), r(T4), r(T3), r(T11)]),
        ch(vec![r(T5), r(T7), r(T11), r(T4)]),
        ch(vec![ctr3(), ctr1(), r(T9), ctr2(), r(T4), r(T9), r(T9), ctr4()]),
        r(Id),
    ])
}

fn is_expr_node(e: &Expr) -> bool {
    match e.kind() {
        Kind::Rel(..) | Kind::Bool(_) => false,
        Kind::Func(Func::Logic(_) | Func::Set(_), _) => false,
        _ => true,
    }
}

/// `fu(rv, measure)`.
pub fn fu(rv: &Expr, less: Less) -> Expr {
    if !is_expr_node(rv) {
        return match rv.kind() {
            Kind::Bool(_) => rv.clone(),
            _ => super::map_args(rv, &mut |a| fu(a, less)),
        };
    }
    let was = rv.clone();
    let mut rv = tr1(rv);
    if has_any(&rv, &[Func::Tan, Func::Cot]) {
        let rv1 = run(&rl1(), &rv, less);
        if less(&rv1, &rv) {
            rv = rv1;
        }
        if has_any(&rv, &[Func::Tan, Func::Cot]) {
            rv = tr2(&rv);
        }
    }
    if has_any(&rv, &[Func::Sin, Func::Cos]) {
        let rv1 = run(&rl2(), &rv, less);
        let rv2 = tr8(&trmorrie(&rv1), true);
        rv = min_by_measure(vec![was, rv, rv1, rv2], less);
    }
    min_by_measure(vec![tr2i(&rv, false), rv], less)
}

/// Constant used by callers that need `pi`.
pub fn pi_expr() -> Expr {
    Expr::constant(Constant::Pi)
}
