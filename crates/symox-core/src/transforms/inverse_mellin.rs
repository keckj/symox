//! The reference inverse Mellin transform: `_rewrite_gamma` and
//! `_inverse_mellin_transform` of `sympy.integrals.transforms`.
//!
//! `F(s)` is rewritten as `fac * prod Gamma(a_i + s) * prod Gamma(b_j - s) /
//! (...) * C**(-s/e)` (the gamma functions placed on the correct side of
//! the fundamental strip), read as a Meijer G function of `C/x**e` and
//! expanded with `hyperexpand`. The convergence condition is that of the
//! Mellin-Barnes integral along the strip ([L] 5.2).
//!
//! Errors keep the reference's classes: an `IntegralTransformError` lets
//! the caller fall back (sum splitting or a held transform), while
//! `MellinTransformStripError`, `NotImplementedError` (gammas across the
//! strip) and `TypeError` (a strip bound `None` where the rewrite needs it,
//! a non integer gamma coefficient) propagate to the user.

use super::util::{factors, has, poly, truth};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::func;
use crate::integrate::meijerint::G;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::special::SpecialFn;

/// The exceptions of the reference inverse Mellin transform.
#[derive(Clone, Debug)]
pub(crate) enum ImtError {
    /// `IntegralTransformError`.
    Transform(String),
    /// `MellinTransformStripError` (a `ValueError`).
    Strip(String),
    NotImplemented(String),
    Type(String),
    Value(String),
}

impl ImtError {
    /// The message passed through `transform_doit`: an
    /// `IntegralTransformError` is a plain message, the other classes are
    /// tagged so the binding raises them unchanged.
    pub(crate) fn message(self) -> String {
        match self {
            ImtError::Transform(m) => m,
            ImtError::Strip(m) => super::fatal("MellinTransformStripError", &m),
            ImtError::NotImplemented(m) => super::fatal("NotImplementedError", &m),
            ImtError::Type(m) => super::fatal("TypeError", &m),
            ImtError::Value(m) => super::fatal("ValueError", &m),
        }
    }
}

type Res<T> = Result<T, ImtError>;

fn transform_error(msg: &str) -> ImtError {
    ImtError::Transform(format!("Inverse Mellin Transform could not be computed: {msg}"))
}

fn int(i: i64) -> Expr {
    Expr::int(i)
}

fn re(e: &Expr) -> Expr {
    func(Func::Re, vec![e.clone()])
}

fn abs(e: &Expr) -> Expr {
    func(Func::Abs, vec![e.clone()])
}

fn gamma(e: &Expr) -> Expr {
    func(Func::Gamma, vec![e.clone()])
}

fn expand(e: &Expr) -> Expr {
    crate::expand::expand(e)
}

fn is_none(e: &Expr) -> bool {
    *e == super::none_plane()
}

/// Every distinct sub-expression `f(...)` with head in `heads` containing
/// `s`, in a deterministic order.
fn atoms(e: &Expr, s: &Expr, heads: &[Func]) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(e, &mut |x| {
        if let Kind::Func(f, _) = x.kind() {
            if heads.contains(f) && has(x, s) && !out.contains(x) {
                out.push(x.clone());
            }
        }
    });
    out.sort_by(crate::sortkey::cmp_default);
    out
}

/// The truth value of `left(c, is_numer)`: `True`, `False`, `None` or an
/// undecided relational (the reference's `c < b` when `a` is `None`).
enum Side {
    Left,
    Right,
    Unknown,
    Rel,
}

struct Strip {
    a: Option<Expr>,
    b: Option<Expr>,
}

impl Strip {
    /// `left(c, is_numer)`: does the pole at `c` lie to the left of the
    /// fundamental strip?
    fn left(&self, c: &Expr, is_numer: bool) -> Res<Side> {
        let c = expand(&re(c));
        let of = |e: Expr| match truth(&e) {
            Some(true) => Side::Left,
            Some(false) => Side::Right,
            None => Side::Rel,
        };
        match (&self.a, &self.b) {
            (None, Some(b)) if b.is_constant(Constant::Infinity) => return Ok(Side::Left),
            (None, Some(b)) => return Ok(of(crate::relational::lt(&c, b))),
            (Some(a), None) => return Ok(of(crate::relational::le(&c, a))),
            (None, None) => return Ok(Side::Left),
            (Some(a), Some(b)) => {
                if truth(&crate::relational::ge(&c, b)) == Some(true) {
                    return Ok(Side::Right);
                }
                if truth(&crate::relational::le(&c, a)) == Some(true) {
                    return Ok(Side::Left);
                }
                if is_numer {
                    return Ok(Side::Unknown);
                }
                let free = |e: &Expr| !crate::traverse::free_symbols(e).is_empty();
                if free(a) || free(b) || free(&c) {
                    return Ok(Side::Unknown);
                }
                Err(ImtError::Strip("Pole inside critical strip?".into()))
            }
        }
    }
}

/// `bool(e)` of a SymPy comparison.
fn bool_of(e: &Expr) -> Res<bool> {
    truth(e).ok_or_else(|| ImtError::Type("cannot determine truth value of Relational".into()))
}

/// `linear_arg`: `arg == a*s + b`, else the "Unrecognised form" error.
fn linear_arg(arg: &Expr, s: &Expr, fact: &Expr) -> Res<(Expr, Expr)> {
    match poly(arg, s) {
        Some(c) if c.len() == 2 => Ok((c[1].clone(), c[0].clone())),
        _ => Err(unrecognised(fact)),
    }
}

fn unrecognised(fact: &Expr) -> ImtError {
    transform_error(&format!("Unrecognised form '{}'.", crate::printer::to_str(fact)))
}

/// `_rewrite_sin((m, n), s, a, b)`: `1/sin(m*s + n) == fac/(g1*g2)` with
/// the gamma functions on the right side of the strip.
fn rewrite_sin(m: &Expr, n: &Expr, s: &Expr, a: &Option<Expr>) -> Res<(Expr, Expr, Expr)> {
    let pi = Expr::pi();
    let m = crate::simptools::expand_mul(&div(m, &pi));
    let n = crate::simptools::expand_mul(&div(n, &pi));
    let Some(a) = a else {
        return Err(ImtError::Type("unsupported operand type(s) for *: 'NegativeOne' and 'NoneType'".into()));
    };
    let re_n = match crate::complex_parts::as_real_imag(&n, true, None) {
        Ok(Some((r, _))) => r,
        _ => re(&n),
    };
    let r = func(Func::Ceiling, vec![sub(&mul(&[neg(&m), a.clone()]), &re_n)]);
    let ms = mul(&[m, s.clone()]);
    Ok((
        gamma(&add(&[ms.clone(), n.clone(), r.clone()])),
        gamma(&sub(&sub(&sub(&Expr::one(), &n), &r), &ms)),
        mul(&[pow(&Expr::minus_one(), &r), pi]),
    ))
}

/// The rewritten function: `((an, ap), (bm, bq), C, e, fac)` with
/// `F(s) == fac * MellinTransform[G(an, ap; bm, bq; C/x**e)]`.
pub(crate) struct Rewritten {
    pub an: Vec<Expr>,
    pub ap: Vec<Expr>,
    pub bm: Vec<Expr>,
    pub bq: Vec<Expr>,
    pub arg: Expr,
    pub exponent: Expr,
    pub fac: Expr,
}

enum Item {
    E(Expr),
    /// An unevaluated `sin(a)` (the reference's `sin(a, evaluate=False)`).
    Sin(Expr),
}

/// `_rewrite_gamma(f, s, a, b)`.
pub(crate) fn rewrite_gamma(f: &Expr, s: &Expr, a: &Option<Expr>, b: &Option<Expr>) -> Res<Rewritten> {
    let pi = Expr::pi();
    // 1) the common multiplier of s in the gamma and trigonometric arguments
    let coeff_of = |g: &Expr| {
        let mut arg = g.args()[0].clone();
        if arg.is_add() {
            arg = crate::decompose::as_independent(&arg, std::slice::from_ref(s), None).1;
        }
        crate::decompose::as_coeff_mul(&arg, std::slice::from_ref(s), false).0
    };
    let mut mults: Vec<Expr> = atoms(f, s, &[Func::Gamma]).iter().map(coeff_of).collect();
    for g in atoms(f, s, &[Func::Sin, Func::Cos, Func::Tan, Func::Cot]) {
        mults.push(div(&coeff_of(&g), &pi));
    }
    let real = |x: &Expr| crate::assumptions::fact(x, "extended_real") == Some(true);
    let mults: Vec<Expr> = mults.into_iter().map(|x| if real(&x) { abs(&x) } else { x }).collect();
    let common = mults.iter().find(|x| !x.is_rational_number()).cloned().unwrap_or_else(Expr::one);
    let mults: Vec<Expr> = mults.iter().map(|x| div(x, &common)).collect();
    if !(mults.iter().all(|x| x.is_rational_number()) && real(&common)) {
        return Err(ImtError::Transform("Gamma Transform could not be computed: Nonrational multiplier".into()));
    }
    let nd: Vec<(num_bigint::BigInt, num_bigint::BigInt)> =
        mults.iter().filter_map(|x| x.as_number().and_then(|n| n.numer_denom())).collect();
    use num_integer::Integer;
    let lcm = nd.iter().fold(num_bigint::BigInt::from(1), |acc, (_, q)| acc.lcm(q));
    let big = |b: num_bigint::BigInt| Expr::num(crate::number::Number::from_bigint(b));
    let mut sm = div(&common, &big(lcm));
    if sm == common && !mults.is_empty() {
        let g = nd.iter().fold(num_bigint::BigInt::from(0), |acc, (p, _)| acc.gcd(p));
        sm = mul(&[common.clone(), big(g)]);
    }
    let f = crate::subs::subs1(f, s, &div(s, &sm));
    let mut fac = div(&Expr::one(), &sm);
    let exponent = div(&Expr::one(), &sm);
    let strip = Strip {
        a: a.as_ref().map(|a| mul(&[a.clone(), sm.clone()])),
        b: b.as_ref().map(|b| mul(&[b.clone(), sm.clone()])),
    };

    // 2) rewrite every factor as gamma functions or exponentials
    let (numer, denom) = crate::simptools::as_numer_denom(&f);
    let mut args: Vec<(Item, bool)> = factors(&numer).into_iter().map(|e| (Item::E(e), true)).collect();
    args.extend(factors(&denom).into_iter().map(|e| (Item::E(e), false)));
    let mut facs: Vec<Expr> = Vec::new();
    let mut dfacs: Vec<Expr> = Vec::new();
    // (a, c) stands for Gamma(a*s + c)
    let mut numer_gammas: Vec<(Expr, Expr)> = Vec::new();
    let mut denom_gammas: Vec<(Expr, Expr)> = Vec::new();
    let mut exponentials: Vec<Expr> = Vec::new();
    let half_pi = div(&pi, &int(2));
    while let Some((item, is_numer)) = args.pop() {
        macro_rules! ug {
            () => {
                if is_numer { &mut numer_gammas } else { &mut denom_gammas }
            };
        }
        macro_rules! lg {
            () => {
                if is_numer { &mut denom_gammas } else { &mut numer_gammas }
            };
        }
        macro_rules! uf {
            () => {
                if is_numer { &mut facs } else { &mut dfacs }
            };
        }
        let fact = match item {
            Item::Sin(a) => {
                let (g1, g2, fac_) = if is_numer {
                    (gamma(&div(&a, &pi)), gamma(&sub(&Expr::one(), &div(&a, &pi))), pi.clone())
                } else {
                    let sin_a = Expr::from_kind(Kind::Func(Func::Sin, vec![a.clone()]));
                    let (m, n) = linear_arg(&a, s, &sin_a)?;
                    rewrite_sin(&m, &n, s, &strip.a)?
                };
                args.push((Item::E(g1), !is_numer));
                args.push((Item::E(g2), !is_numer));
                uf!().push(fac_);
                continue;
            }
            Item::E(e) => e,
        };
        if !has(&fact, s) {
            uf!().push(fact);
            continue;
        }
        let power = match fact.kind() {
            Kind::Pow(b, e) if b.is_constant(Constant::E) => Some((exp_polar_one(), e.clone())),
            Kind::Pow(b, e) => Some((b.clone(), e.clone())),
            Kind::Func(Func::Exp, a) => Some((exp_polar_one(), a[0].clone())),
            _ => None,
        };
        if let Some((base, exp_)) = power {
            if let Some(k) = exp_.as_i64().filter(|_| exp_.is_integer_number()) {
                let cond = if k < 0 { !is_numer } else { is_numer };
                for _ in 0..k.unsigned_abs() {
                    args.push((Item::E(base.clone()), cond));
                }
                continue;
            }
            if !has(&base, s) {
                let (a, b) = linear_arg(&exp_, s, &fact)?;
                let base = if is_numer { base } else { div(&Expr::one(), &base) };
                exponentials.push(pow(&base, &a));
                facs.push(pow(&base, &b));
                continue;
            }
            return Err(unrecognised(&fact));
        }
        if let Some(c) = poly(&fact, s) {
            let degree = c.len() - 1;
            if degree != 1 {
                let coeff = c[degree].clone();
                let rs = crate::poly::roots::roots(&fact, std::slice::from_ref(s)).unwrap_or_default();
                let roots: Vec<Expr> = if rs.len() == degree {
                    rs.into_iter().map(|(r, _)| r).collect()
                } else if rs.iter().map(|(_, k)| *k as usize).sum::<usize>() == degree {
                    rs.into_iter().flat_map(|(r, k)| std::iter::repeat_n(r, k as usize)).collect()
                } else {
                    return Err(transform_error(&format!("Could not factor '{}'.", crate::printer::to_str(&fact))));
                };
                uf!().push(coeff);
                for r in roots {
                    args.push((Item::E(sub(s, &r)), is_numer));
                }
                continue;
            }
            let (a, c) = (c[1].clone(), c[0].clone());
            uf!().push(a.clone());
            let c = neg(&div(&c, &a));
            let left = match strip.left(&c, is_numer)? {
                Side::Left => true,
                Side::Right | Side::Unknown => false,
                Side::Rel => return Err(ImtError::Type("cannot determine truth value of Relational".into())),
            };
            if left {
                ug!().push((Expr::one(), add(&[neg(&c), Expr::one()])));
                lg!().push((Expr::one(), neg(&c)));
            } else {
                uf!().push(Expr::minus_one());
                ug!().push((Expr::minus_one(), add(&[c.clone(), Expr::one()])));
                lg!().push((Expr::minus_one(), c));
            }
            continue;
        }
        let Kind::Func(head, fargs) = fact.kind() else { return Err(unrecognised(&fact)) };
        let arg = fargs[0].clone();
        match head {
            Func::Gamma => {
                let (a, b) = linear_arg(&arg, s, &fact)?;
                if is_numer {
                    let pole = neg(&div(&b, &a));
                    let positive = bool_of(&crate::relational::gt(&a, &Expr::zero()))?;
                    let across = if positive {
                        matches!(strip.left(&pole, is_numer)?, Side::Right)
                    } else {
                        bool_of(&crate::relational::lt(&a, &Expr::zero()))?
                            && matches!(strip.left(&pole, is_numer)?, Side::Left)
                    };
                    if across {
                        return Err(ImtError::NotImplemented("Gammas partially over the strip.".into()));
                    }
                }
                ug!().push((a, b));
            }
            Func::Sin => args.push((Item::Sin(arg), is_numer)),
            Func::Tan => {
                args.push((Item::Sin(arg.clone()), is_numer));
                args.push((Item::Sin(sub(&half_pi, &arg)), !is_numer));
            }
            Func::Cos => args.push((Item::Sin(sub(&half_pi, &arg)), is_numer)),
            Func::Cot => {
                args.push((Item::Sin(sub(&half_pi, &arg)), is_numer));
                args.push((Item::Sin(arg), !is_numer));
            }
            _ => return Err(unrecognised(&fact)),
        }
    }
    fac = mul(&[fac, mul(&facs), pow(&mul(&dfacs), &Expr::minus_one())]);

    // 3) bring every gamma function to the argument c + s or c - s
    let (mut an, mut ap, mut bm, mut bq) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let two_pi = mul(&[int(2), pi.clone()]);
    for is_numer in [true, false] {
        let mut gammas = if is_numer { std::mem::take(&mut numer_gammas) } else { std::mem::take(&mut denom_gammas) };
        while let Some((a, c)) = gammas.pop() {
            if a != Expr::one() && a != Expr::minus_one() {
                // the multiplication theorem
                let p = abs(&a);
                let newa = div(&a, &p);
                let newc = div(&c, &p);
                let Some(pi64) = a.as_i64().filter(|_| a.is_integer_number()).map(i64::abs) else {
                    return Err(ImtError::Type("a is not an integer".into()));
                };
                for k in 0..pi64 {
                    gammas.push((newa.clone(), add(&[newc.clone(), Expr::rational(k, pi64)])));
                }
                let factor = mul(&[
                    pow(&two_pi, &Expr::rational(1 - pi64, 2)),
                    pow(&p, &sub(&c, &Expr::rational(1, 2))),
                ]);
                if is_numer {
                    fac = mul(&[fac, factor]);
                    exponentials.push(pow(&p, &a));
                } else {
                    fac = div(&fac, &factor);
                    exponentials.push(pow(&p, &neg(&a)));
                }
                continue;
            }
            let plus = a == Expr::one();
            match (is_numer, plus) {
                (true, true) => an.push(sub(&Expr::one(), &c)),
                (true, false) => bm.push(c),
                (false, true) => bq.push(sub(&Expr::one(), &c)),
                (false, false) => ap.push(c),
            }
        }
    }

    // 5) combine the exponentials
    let arg = mul(&exponentials);
    for l in [&mut an, &mut ap, &mut bm, &mut bq] {
        l.sort_by(crate::sortkey::cmp_default);
    }
    Ok(Rewritten { an, ap, bm, bq, arg, exponent, fac })
}

fn exp_polar_one() -> Expr {
    func(Func::ExpPolar, vec![Expr::one()])
}

/// The heads `InverseMellinTransform._compute_transform` accepts in a
/// function of `s`.
fn allowed(f: &Func) -> bool {
    match f {
        Func::Exp | Func::Gamma | Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Cosh | Func::Sinh | Func::Tanh
        | Func::Coth => true,
        Func::Special(SpecialFn::Factorial | SpecialFn::RisingFactorial) => true,
        // not SymPy Function instances
        Func::Logic(_) | Func::Set(_) | Func::Integral | Func::Sum | Func::Product | Func::Order | Func::Subs
        | Func::Limit | Func::GenericExpr | Func::Predicate(_) | Func::AppliedPredicate | Func::MatrixElement(_)
        | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _) => true,
        Func::Undefined(_) => f.name() == "Tuple",
        _ => false,
    }
}

/// The component check of `InverseMellinTransform._compute_transform`.
fn check_components(f: &Expr, s: &Expr) -> Res<()> {
    let mut bad: Option<Expr> = None;
    crate::traverse::preorder(f, &mut |x| {
        if bad.is_some() || !has(x, s) {
            return;
        }
        let rejected = match x.kind() {
            Kind::Func(h, _) => !allowed(h),
            Kind::Piecewise(_) => true,
            _ => false,
        };
        if rejected {
            bad = Some(x.clone());
        }
    });
    match bad {
        Some(b) => Err(transform_error(&format!("Component {} not recognised.", crate::printer::to_str(&b)))),
        None => Ok(()),
    }
}

/// `factor(res, gens=res.atoms(Heaviside))` of a sum of inverse
/// transforms. The generators only fix the main variables of the
/// factorisation, so the plain factorisation gives the same product (the
/// reference raises TypeError for a single Heaviside atom).
fn factor_heaviside(e: &Expr) -> Expr {
    crate::poly::api::factor(e).unwrap_or_else(|_| e.clone())
}

/// `hyperexpand` with the reference's error classes.
fn hyperexpand(g: &Expr) -> Res<Expr> {
    crate::special::hyperexpand_options(g, false, Some("default"), None).map_err(|(kind, m)| match kind {
        "NotImplementedError" => transform_error("Could not calculate integral"),
        "TypeError" => ImtError::Type(m),
        "ValueError" => ImtError::Value(m),
        _ => ImtError::NotImplemented(m),
    })
}

/// `F.rewrite(gamma)`, including the combinatorial functions the
/// component check admits (`factorial(n) = gamma(n + 1)` and the piecewise
/// form of `RisingFactorial`).
fn rewrite_gamma_functions(f: &Expr) -> Expr {
    let f = crate::simptools::bottom_up(f, &mut |e| match e.kind() {
        Kind::Func(Func::Special(SpecialFn::Factorial), a) => gamma(&add(&[a[0].clone(), Expr::one()])),
        Kind::Func(Func::Special(SpecialFn::RisingFactorial), a) => {
            let (x, k) = (&a[0], &a[1]);
            crate::relational::piecewise(vec![
                (div(&gamma(&add(&[x.clone(), k.clone()])), &gamma(x)), crate::relational::gt(x, &Expr::zero())),
                (
                    div(
                        &mul(&[pow(&Expr::minus_one(), k), gamma(&sub(&Expr::one(), x))]),
                        &gamma(&sub(&sub(&Expr::one(), k), x)),
                    ),
                    Expr::boolean(true),
                ),
            ])
        }
        _ => e.clone(),
    });
    crate::special::rewrite::rewrite(&f, &[], "gamma", true)
}

/// `_inverse_mellin_transform(F, s, x, strip, as_meijerg)`: `(f, cond)`.
pub(crate) fn inverse_mellin_core(f: &Expr, s: &Expr, x_: &Expr, strip: (&Option<Expr>, &Option<Expr>), as_meijerg: bool) -> Res<(Expr, Expr)> {
    let x = Expr::dummy("t", crate::symbol::Assumptions::POSITIVE);
    let f = rewrite_gamma_functions(f);
    let candidates = [
        crate::poly::api::factor(&f).unwrap_or_else(|_| f.clone()),
        crate::simptools::expand_mul(&f),
        expand(&f),
    ];
    for g in candidates {
        if g.is_add() {
            let mut ress = Vec::new();
            let mut conds = Vec::new();
            for t in g.args() {
                let (r, c) = inverse_mellin_core(&t, s, &x, strip, as_meijerg)?;
                ress.push(r);
                conds.push(c);
            }
            let mut res = add(&ress);
            if !as_meijerg {
                res = factor_heaviside(&res);
            }
            return Ok((crate::subs::subs1(&res, &x, x_), crate::logic::and(conds)));
        }
        let r = match rewrite_gamma(&g, s, strip.0, strip.1) {
            Ok(r) => r,
            Err(ImtError::Transform(_)) => continue,
            Err(e) => return Err(e),
        };
        let z = div(&r.arg, &pow(&x, &r.exponent));
        let Some(gf) = G::new(r.an, r.ap, r.bm, r.bq, z) else { continue };
        let gexpr = gf.to_expr();
        let mut h = if as_meijerg { gexpr.clone() } else { hyperexpand(&gexpr)? };
        if let Kind::Piecewise(pieces) = h.kind() {
            if pieces.len() == 3 {
                let c = abs(&r.arg);
                h = add(&[
                    mul(&[crate::functions::heaviside(&sub(&x, &c)), pieces[0].0.clone()]),
                    mul(&[crate::functions::heaviside(&sub(&c, &x)), pieces[1].0.clone()]),
                ]);
            }
        }
        // convergence of the Mellin-Barnes integral along the strip ([L] 5.2)
        let pi = Expr::pi();
        let absarg = abs(&func(Func::Arg, vec![gf.z.clone()]));
        let dpi = mul(&[gf.delta(), pi]);
        let c1 = crate::relational::lt(&absarg, &dpi);
        let c2 = crate::logic::and(vec![
            crate::logic::or(vec![
                Expr::boolean(gf.p() != gf.q()),
                crate::relational::le(&add(&[re(&gf.nu()), Expr::one()]), &Expr::zero()),
            ]),
            Expr::boolean(absarg == dpi),
        ]);
        let cond = crate::logic::or(vec![c1, c2]);
        if truth(&cond) == Some(false) {
            return Err(transform_error("does not converge"));
        }
        return Ok((crate::subs::subs1(&mul(&[h, r.fac]), &x, x_), cond));
    }
    Err(transform_error(""))
}

/// `InverseMellinTransform(F, s, x, a, b).doit()` without the sum
/// fallback: the component check followed by `_inverse_mellin_transform`.
pub(crate) fn compute(f: &Expr, s: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Res<(Expr, Expr)> {
    check_components(f, s)?;
    let a = if is_none(a) { None } else { Some(a.clone()) };
    let b = if is_none(b) { None } else { Some(b.clone()) };
    inverse_mellin_core(f, s, x, (&a, &b), false)
}

/// `compute` returning the Meijer G form (`as_meijerg=True`).
pub(crate) fn compute_meijerg(f: &Expr, s: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Res<(Expr, Expr)> {
    check_components(f, s)?;
    let a = if is_none(a) { None } else { Some(a.clone()) };
    let b = if is_none(b) { None } else { Some(b.clone()) };
    inverse_mellin_core(f, s, x, (&a, &b), true)
}
