//! Convergence conditions of the Meijer G integrals (`_check_antecedents_1`,
//! `_check_antecedents`) and their simplification (`_condsimp`,
//! `_eval_cond`, `_my_unpolarify`), as in the reference.

use super::{get_coeff_exp, is_false, G};
use crate::expr::RelOp;
use crate::functions::{abs, cos, re, sin, sign};
use crate::logic::LogicOp;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};
use crate::relational::{eq, ge, gt, le, lt, ne};
use crate::special::polar::unpolarify;
use crate::special::SpecialFn;
use crate::{Expr, Func, Kind};

fn b(v: bool) -> Expr {
    Expr::boolean(v)
}
fn int(i: i64) -> Expr {
    Expr::int(i)
}
fn rat(n: i64, d: i64) -> Expr {
    Expr::rational(n, d)
}
fn and(v: Vec<Expr>) -> Expr {
    crate::logic::and(v)
}
fn or(v: Vec<Expr>) -> Expr {
    crate::logic::or(v)
}
/// `unbranched_argument(e)`.
fn ua(e: &Expr) -> Expr {
    crate::special::apply(SpecialFn::PeriodicArgument, vec![e.clone(), Expr::oo()])
}
fn pi() -> Expr {
    Expr::pi()
}
fn is_pos(e: &Expr) -> bool {
    crate::assumptions::is_positive(e) == Some(true)
}
fn is_neg(e: &Expr) -> bool {
    crate::assumptions::is_negative(e) == Some(true)
}

/// `_check_antecedents_1`: when the integral of one G-function over (0, oo)
/// exists.
pub(crate) fn check_antecedents_1(g: &G, x: &Expr, helper: bool) -> Option<Expr> {
    let delta = g.delta();
    let (eta, _) = get_coeff_exp(&g.z, x)?;
    let (m, n, p, q) = (g.bm.len() as i64, g.an.len() as i64, g.p() as i64, g.q() as i64);
    if p > q {
        let tr = |l: &[Expr]| l.iter().map(|a| sub(&Expr::one(), a)).collect::<Vec<_>>();
        let flipped = G::new(tr(&g.bm), tr(&g.bo), tr(&g.an), tr(&g.ao), div(x, &eta))?;
        return check_antecedents_1(&flipped, x, false);
    }
    let mut tmp: Vec<Expr> = g.bm.iter().map(|bb| lt(&neg(&re(bb)), &Expr::one())).collect();
    tmp.extend(g.an.iter().map(|a| gt(&sub(&Expr::one(), &re(a)), &Expr::one())));
    let cond_3 = and(tmp.clone());
    let cond_4 = gt(&add(&[neg(&re(&g.nu())), rat(q + 1 - p, 2)]), &int(q - p));
    let ae = abs(&ua(&eta));
    let mut conds = Vec::new();
    let tmp1 = vec![b(1 <= n), b(p < q), b(1 <= m)];
    let tmp2 = vec![b(1 <= p), b(1 <= m), b(q == p + 1), b(!(n == 0 && m == p + 1))];
    let mut tmp3 = vec![b(1 <= p), b(q == p)];
    // ceiling(delta/2) with delta = (2*(m+n) - p - q)/2
    let twice = 2 * (m + n) - p - q;
    let kmax = (twice as f64 / 4.0).ceil() as i64;
    for k in 0..=kmax {
        tmp3.push(ne(&ae, &mul2(&sub(&delta, &int(2 * k)), &pi())));
    }
    let tmp = vec![gt(&delta, &Expr::zero()), lt(&ae, &mul2(&delta, &pi()))];
    let extra = if helper { vec![] } else { vec![ne(&eta, &Expr::zero()), cond_3.clone()] };
    for t in [tmp1, tmp2, tmp3] {
        let mut v = t;
        v.extend(tmp.iter().cloned());
        v.extend(extra.iter().cloned());
        conds.push(and(v));
    }
    let extra = if helper { vec![] } else { vec![cond_3.clone()] };
    let mut c2 = vec![b(n == 0), b(p + 1 <= m), b(m <= q), lt(&ae, &mul2(&delta, &pi()))];
    c2.extend(extra.iter().cloned());
    conds.push(and(c2));
    let extra = if helper { vec![] } else { vec![cond_3.clone(), cond_4] };
    let mut c3a = vec![b(p < q), b(1 <= m), gt(&delta, &Expr::zero()), eq(&ae, &mul2(&delta, &pi()))];
    c3a.extend(extra.iter().cloned());
    let mut c3b = vec![b(p <= q - 2), eq(&delta, &Expr::zero()), eq(&ae, &Expr::zero())];
    c3b.extend(extra.iter().cloned());
    conds.push(and(c3a));
    conds.push(and(c3b));
    let mut ce = vec![b(p == q), eq(&delta, &Expr::zero()), eq(&ua(&eta), &Expr::zero()), ne(&eta, &Expr::zero())];
    if !helper {
        ce.push(cond_3.clone());
    }
    let s: Vec<Expr> = g.ap().iter().zip(g.bq().iter()).map(|(a, bb)| sub(bb, a)).collect();
    ce.push(lt(&re(&add(&s)), &Expr::zero()));
    conds.push(and(ce));
    let mut ce2 = vec![gt(&delta, &Expr::zero()), lt(&ae, &mul2(&delta, &pi()))];
    if !helper {
        ce2.push(cond_3);
    }
    conds.push(and(ce2));
    Some(or(conds))
}

/// `_check_antecedents`: when the integral of a product of two G-functions
/// over (0, oo) exists.
pub(crate) fn check_antecedents(g1: &G, g2: &G, x: &Expr) -> Option<Expr> {
    let (sigma, _) = get_coeff_exp(&g1.z, x)?;
    let (omega, _) = get_coeff_exp(&g2.z, x)?;
    let (s, t, u, v) = (g1.bm.len() as i64, g1.an.len() as i64, g1.p() as i64, g1.q() as i64);
    let (m, n, p, q) = (g2.bm.len() as i64, g2.an.len() as i64, g2.p() as i64, g2.q() as i64);
    // Twice b* and c*, kept integral.
    let bstar2 = 2 * (s + t) - (u + v);
    let cstar2 = 2 * (m + n) - (p + q);
    let bstar = rat(bstar2, 2);
    let cstar = rat(cstar2, 2);
    let rho = add(&[g1.nu(), rat(u - v, 2), Expr::one()]);
    let mu = add(&[g2.nu(), rat(p - q, 2), Expr::one()]);
    let phi = q - p - (v - u);
    let eta = sub(&sub(&sub(&Expr::one(), &int(v - u)), &mu), &rho);
    let c1 = {
        let mut ok = true;
        for g in [g1, g2] {
            for i in &g.an {
                for j in &g.bm {
                    let d = sub(i, j);
                    if crate::assumptions::is_integer(&d) == Some(true) && is_pos(&d) {
                        ok = false;
                    }
                }
            }
        }
        b(ok)
    };
    let c2 = and(g1.bm.iter().flat_map(|i| g2.bm.iter().map(move |j| gt(&re(&add(&[Expr::one(), i.clone(), j.clone()])), &Expr::zero()))).collect());
    let c3 = and(g1.an.iter().flat_map(|i| g2.an.iter().map(move |j| lt(&re(&add(&[Expr::one(), i.clone(), j.clone()])), &int(2)))).collect());
    let m32 = rat(-3, 2);
    let c4 = and(g1.an.iter().map(|i| gt(&sub(&mul2(&int(p - q), &re(i)), &re(&mu)), &m32)).collect());
    let c5 = and(g1.bm.iter().map(|i| gt(&sub(&mul2(&int(p - q), &re(&add(&[Expr::one(), i.clone()]))), &re(&mu)), &m32)).collect());
    let c6 = and(g2.an.iter().map(|i| gt(&sub(&mul2(&int(u - v), &re(i)), &re(&rho)), &m32)).collect());
    let c7 = and(g2.bm.iter().map(|i| gt(&sub(&mul2(&int(u - v), &re(&add(&[Expr::one(), i.clone()]))), &re(&rho)), &m32)).collect());
    let inner = add(&[
        mul2(&sub(&rho, &Expr::one()), &int(q - p)),
        int((v - u) * (q - p)),
        mul2(&sub(&mu, &Expr::one()), &int(v - u)),
    ]);
    let c8 = gt(&add(&[int(phi.abs()), mul2(&int(2), &re(&inner))]), &Expr::zero());
    let c9 = gt(&sub(&int(phi.abs()), &mul2(&int(2), &re(&inner))), &Expr::zero());
    let asig = abs(&ua(&sigma));
    let aom = abs(&ua(&omega));
    let c10 = lt(&asig, &mul2(&bstar, &pi()));
    let c11 = eq(&asig, &mul2(&bstar, &pi()));
    let c12 = lt(&aom, &mul2(&cstar, &pi()));
    let c13 = eq(&aom, &mul2(&cstar, &pi()));
    let z0 = crate::functions::exp(&mul(&[neg(&rat(bstar2 + cstar2, 2)), pi(), Expr::i()]));
    let zos = unpolarify(&div(&mul2(&z0, &omega), &sigma), false);
    let zso = unpolarify(&div(&mul2(&z0, &sigma), &omega), false);
    let bc = rat(bstar2 + cstar2, 2);
    let c14 = if zos == pow(&zso, &Expr::minus_one()) {
        and(vec![
            b(phi == 0),
            le(&bc, &Expr::one()),
            or(vec![
                ne(&zos, &Expr::one()),
                lt(&re(&add(&[mu.clone(), rho.clone(), int(v - u)])), &Expr::one()),
                lt(&re(&add(&[mu.clone(), rho.clone(), int(q - p)])), &Expr::one()),
            ]),
        ])
    } else {
        let cond = |z: &Expr| if z.is_one() { b(false) } else { lt(&abs(&crate::functions::arg(&sub(&Expr::one(), z))), &pi()) };
        let a = and(vec![
            b(phi == 0),
            le(&sub(&bc, &Expr::one()), &Expr::zero()),
            or(vec![
                and(vec![ne(&zos, &Expr::one()), cond(&zos)]),
                and(vec![lt(&re(&add(&[mu.clone(), rho.clone(), int(v - u)])), &Expr::one()), eq(&zos, &Expr::one())]),
            ]),
        ]);
        let alt = and(vec![
            b(phi == 0),
            le(&sub(&bc, &Expr::one()), &Expr::zero()),
            or(vec![
                and(vec![ne(&zso, &Expr::one()), cond(&zso)]),
                and(vec![lt(&re(&add(&[mu.clone(), rho.clone(), int(q - p)])), &Expr::one()), eq(&zso, &Expr::one())]),
            ]),
        ]);
        or(vec![a, alt])
    };
    let c15 = if q == p || v == u {
        // psi or theta is NaN: the reference catches the TypeError.
        b(false)
    } else {
        let psi = div(&add(&[mul2(&pi(), &int(q - m - n)), aom.clone()]), &int(q - p));
        let theta = div(&add(&[mul2(&pi(), &int(v - s - t)), asig.clone()]), &int(v - u));
        let po = pow(&abs(&omega), &rat(1, q - p));
        let ps = pow(&abs(&sigma), &rat(1, v - u));
        let lambda_c = add(&[mul(&[int(q - p), po.clone(), cos(&psi)]), mul(&[int(v - u), ps.clone(), cos(&theta)])]);
        if !is_false(&eval_cond(&gt(&lambda_c, &Expr::zero()))) {
            gt(&lambda_c, &Expr::zero())
        } else {
            let ls0 = |k1: &Expr, k2: &Expr| {
                add(&[mul(&[k1.clone(), int(q - p), po.clone(), sin(&psi)]), mul(&[k2.clone(), int(v - u), ps.clone(), sin(&theta)])])
            };
            let (one, mone) = (Expr::one(), Expr::minus_one());
            let (ss, so) = (sign(&ua(&sigma)), sign(&ua(&omega)));
            let lambda_s = crate::relational::piecewise(vec![
                (mul2(&ls0(&one, &one), &ls0(&mone, &mone)), and(vec![eq(&ua(&sigma), &Expr::zero()), eq(&ua(&omega), &Expr::zero())])),
                (mul2(&ls0(&so, &one), &ls0(&so, &mone)), and(vec![eq(&ua(&sigma), &Expr::zero()), ne(&ua(&omega), &Expr::zero())])),
                (mul2(&ls0(&one, &ss), &ls0(&mone, &ss)), and(vec![ne(&ua(&sigma), &Expr::zero()), eq(&ua(&omega), &Expr::zero())])),
                (ls0(&so, &ss), b(true)),
            ]);
            or(vec![
                gt(&lambda_c, &Expr::zero()),
                and(vec![eq(&lambda_c, &Expr::zero()), ne(&lambda_s, &Expr::zero()), gt(&re(&eta), &Expr::minus_one())]),
                and(vec![eq(&lambda_c, &Expr::zero()), eq(&lambda_s, &Expr::zero()), gt(&re(&eta), &Expr::zero())]),
            ])
        }
    };
    let bpos = bstar2 > 0;
    let cpos = cstar2 > 0;
    let bnn = bstar2 >= 0;
    let cnn = cstar2 >= 0;
    let (sigpos, ompos) = (is_pos(&sigma), is_pos(&omega));
    let base = |v: Vec<Expr>| and(v);
    let (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15) =
        (&c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10, &c11, &c12, &c13, &c14, &c15);
    let k = |e: &Expr| e.clone();
    let mut conds = vec![
        base(vec![b(m * n * s * t != 0), b(bpos), b(cpos), k(c1), k(c2), k(c3), k(c10), k(c12)]),
        base(vec![b(u == v), b(bstar2 == 0), b(cpos), b(sigpos), lt(&re(&rho), &Expr::one()), k(c1), k(c2), k(c3), k(c12)]),
        base(vec![b(p == q), b(cstar2 == 0), b(bpos), b(ompos), lt(&re(&mu), &Expr::one()), k(c1), k(c2), k(c3), k(c10)]),
        base(vec![
            b(p == q), b(u == v), b(bstar2 == 0), b(cstar2 == 0), b(sigpos), b(ompos),
            lt(&re(&mu), &Expr::one()), lt(&re(&rho), &Expr::one()), ne(&sigma, &omega), k(c1), k(c2), k(c3),
        ]),
        base(vec![
            b(p == q), b(u == v), b(bstar2 == 0), b(cstar2 == 0), b(sigpos), b(ompos),
            lt(&re(&add(&[mu.clone(), rho.clone()])), &Expr::one()), ne(&omega, &sigma), k(c1), k(c2), k(c3),
        ]),
        base(vec![b(p > q), b(s > 0), b(bpos), b(cnn), k(c1), k(c2), k(c3), k(c5), k(c10), k(c13)]),
        base(vec![b(p < q), b(t > 0), b(bpos), b(cnn), k(c1), k(c2), k(c3), k(c4), k(c10), k(c13)]),
        base(vec![b(u > v), b(m > 0), b(cpos), b(bnn), k(c1), k(c2), k(c3), k(c7), k(c11), k(c12)]),
        base(vec![b(u < v), b(n > 0), b(cpos), b(bnn), k(c1), k(c2), k(c3), k(c6), k(c11), k(c12)]),
        base(vec![b(p > q), b(u == v), b(bstar2 == 0), b(cnn), b(sigpos), lt(&re(&rho), &Expr::one()), k(c1), k(c2), k(c3), k(c5), k(c13)]),
        base(vec![b(p < q), b(u == v), b(bstar2 == 0), b(cnn), b(sigpos), lt(&re(&rho), &Expr::one()), k(c1), k(c2), k(c3), k(c4), k(c13)]),
        base(vec![b(p == q), b(u > v), b(bnn), b(cstar2 == 0), b(ompos), lt(&re(&mu), &Expr::one()), k(c1), k(c2), k(c3), k(c7), k(c11)]),
        base(vec![b(p == q), b(u < v), b(bnn), b(cstar2 == 0), b(ompos), lt(&re(&mu), &Expr::one()), k(c1), k(c2), k(c3), k(c6), k(c11)]),
        base(vec![b(p < q), b(u > v), b(bnn), b(cnn), k(c1), k(c2), k(c3), k(c4), k(c7), k(c11), k(c13)]),
        base(vec![b(p > q), b(u < v), b(bnn), b(cnn), k(c1), k(c2), k(c3), k(c5), k(c6), k(c11), k(c13)]),
        base(vec![b(p > q), b(u > v), b(bnn), b(cnn), k(c1), k(c2), k(c3), k(c5), k(c7), k(c8), k(c11), k(c13), k(c14)]),
        base(vec![b(p < q), b(u < v), b(bnn), b(cnn), k(c1), k(c2), k(c3), k(c4), k(c6), k(c9), k(c11), k(c13), k(c14)]),
        base(vec![b(t == 0), b(s > 0), b(bpos), b(phi > 0), k(c1), k(c2), k(c10)]),
        base(vec![b(s == 0), b(t > 0), b(bpos), b(phi < 0), k(c1), k(c3), k(c10)]),
        base(vec![b(n == 0), b(m > 0), b(cpos), b(phi < 0), k(c1), k(c2), k(c12)]),
        base(vec![b(m == 0), b(n > 0), b(cpos), b(phi > 0), k(c1), k(c3), k(c12)]),
        base(vec![b(s * t == 0), b(bpos), b(cpos), k(c1), k(c2), k(c3), k(c10), k(c12)]),
        base(vec![b(m * n == 0), b(bpos), b(cpos), k(c1), k(c2), k(c3), k(c10), k(c12)]),
    ];
    let mt1 = check_antecedents_1(g1, x, true)?;
    let mt2 = check_antecedents_1(g2, x, true)?;
    conds.push(base(vec![mt2.clone(), b(t == 0), b(u < s), b(bpos), k(c10), k(c1), k(c2), k(c3)]));
    conds.push(base(vec![mt2, b(s == 0), b(v < t), b(bpos), k(c10), k(c1), k(c2), k(c3)]));
    conds.push(base(vec![mt1.clone(), b(n == 0), b(p < m), b(cpos), k(c12), k(c1), k(c2), k(c3)]));
    conds.push(base(vec![mt1, b(m == 0), b(q < n), b(cpos), k(c12), k(c1), k(c2), k(c3)]));
    let r = or(conds.clone());
    if !is_false(&eval_cond(&r)) {
        return Some(r);
    }
    let lt_pi = |e: &Expr, k: i64| lt(e, &mul2(&int(k), &pi()));
    let (bneg, cneg) = (bstar2 < 0, cstar2 < 0);
    let cpi_lt = lt(&mul2(&cstar, &pi()), &aom);
    let bpi_lt = lt(&mul2(&bstar, &pi()), &asig);
    conds.extend([
        base(vec![b(m + n > p), b(t == 0), b(phi == 0), b(s > 0), b(bpos), b(cneg), lt_pi(&aom, m + n - p + 1), k(c1), k(c2), k(c10), k(c14), k(c15)]),
        base(vec![b(m + n > q), b(s == 0), b(phi == 0), b(t > 0), b(bpos), b(cneg), lt_pi(&aom, m + n - q + 1), k(c1), k(c3), k(c10), k(c14), k(c15)]),
        base(vec![b(p == q - 1), b(t == 0), b(phi == 0), b(s > 0), b(bpos), b(cnn), cpi_lt.clone(), k(c1), k(c2), k(c10), k(c14), k(c15)]),
        base(vec![b(p == q + 1), b(s == 0), b(phi == 0), b(t > 0), b(bpos), b(cnn), cpi_lt.clone(), k(c1), k(c3), k(c10), k(c14), k(c15)]),
        base(vec![b(p < q - 1), b(t == 0), b(phi == 0), b(s > 0), b(bpos), b(cnn), cpi_lt.clone(), lt_pi(&aom, m + n - p + 1), k(c1), k(c2), k(c10), k(c14), k(c15)]),
        base(vec![b(p > q + 1), b(s == 0), b(phi == 0), b(t > 0), b(bpos), b(cnn), cpi_lt, lt_pi(&aom, m + n - q + 1), k(c1), k(c3), k(c10), k(c14), k(c15)]),
        base(vec![b(n == 0), b(phi == 0), b(s + t > 0), b(m > 0), b(cpos), b(bneg), lt_pi(&asig, s + t - u + 1), k(c1), k(c2), k(c12), k(c14), k(c15)]),
        base(vec![b(m == 0), b(phi == 0), b(s + t > v), b(n > 0), b(cpos), b(bneg), lt_pi(&asig, s + t - v + 1), k(c1), k(c3), k(c12), k(c14), k(c15)]),
        base(vec![b(n == 0), b(phi == 0), b(u == v - 1), b(m > 0), b(cpos), b(bnn), bpi_lt.clone(), lt(&asig, &mul2(&add(&[bstar.clone(), Expr::one()]), &pi())), k(c1), k(c2), k(c12), k(c14), k(c15)]),
        base(vec![b(m == 0), b(phi == 0), b(u == v + 1), b(n > 0), b(cpos), b(bnn), bpi_lt.clone(), lt(&asig, &mul2(&add(&[bstar.clone(), Expr::one()]), &pi())), k(c1), k(c3), k(c12), k(c14), k(c15)]),
        base(vec![b(n == 0), b(phi == 0), b(u < v - 1), b(m > 0), b(cpos), b(bnn), bpi_lt.clone(), lt_pi(&asig, s + t - u + 1), k(c1), k(c2), k(c12), k(c14), k(c15)]),
        base(vec![b(m == 0), b(phi == 0), b(u > v + 1), b(n > 0), b(cpos), b(bnn), bpi_lt, lt_pi(&asig, s + t - v + 1), k(c1), k(c3), k(c12), k(c14), k(c15)]),
    ]);
    let _ = (is_neg as fn(&Expr) -> bool, ge as fn(&Expr, &Expr) -> Expr);
    Some(or(conds))
}

// ---------------------------------------------------------------------
// _condsimp

fn raw_rel(op: RelOp, a: Expr, b: Expr) -> Expr {
    Expr::from_kind(Kind::Rel(op, a, b))
}

struct Rule {
    op: LogicOp,
    args: Vec<Expr>,
    to: Expr,
    /// Rule 10: the second argument is `And(r, p < q)`.
    and_r: bool,
}

fn rules(p: &Expr, q: &Expr) -> Vec<Rule> {
    use RelOp::*;
    let argp = crate::functions::arg(p);
    let pi = Expr::pi();
    let pi2 = div(&pi, &int(2));
    let two_arg = mul2(&int(2), &argp);
    // SymPy's matcher solves `p**2/2 + 1 == X` for any X (p = sqrt(2*X - 2)),
    // so these patterns accept any argument: match X with `q` directly.
    let pp = q.clone();
    let cosabs = crate::functions::cos(&abs(&argp));
    let r = |op, args, to| Rule { op, args, to, and_r: false };
    let eq_arg = |shift: Expr| raw_rel(Eq, add(&[argp.clone(), shift]), Expr::zero());
    vec![
        r(LogicOp::Or, vec![raw_rel(Eq, p.clone(), q.clone()), raw_rel(Lt, p.clone(), q.clone())], raw_rel(Le, p.clone(), q.clone())),
        r(LogicOp::And, vec![raw_rel(Le, abs(&argp), pi.clone()), raw_rel(Le, abs(&sub(&argp, &mul2(&int(2), &pi))), pi.clone())], eq_arg(neg(&pi))),
        r(LogicOp::And, vec![raw_rel(Le, abs(&add(&[two_arg.clone(), pi.clone()])), pi.clone()), raw_rel(Le, abs(&sub(&two_arg, &pi)), pi.clone())], eq_arg(Expr::zero())),
        r(LogicOp::And, vec![raw_rel(Lt, abs(&add(&[two_arg.clone(), pi.clone()])), pi.clone()), raw_rel(Le, abs(&sub(&two_arg, &pi)), pi.clone())], b(false)),
        r(LogicOp::And, vec![raw_rel(Le, abs(&sub(&argp, &pi2)), pi2.clone()), raw_rel(Le, abs(&add(&[argp.clone(), pi2.clone()])), pi2.clone())], eq_arg(Expr::zero())),
        r(LogicOp::And, vec![raw_rel(Le, abs(&sub(&argp, &pi2)), pi2.clone()), raw_rel(Lt, abs(&add(&[argp.clone(), pi2.clone()])), pi2.clone())], b(false)),
        r(LogicOp::And, vec![raw_rel(Lt, abs(&crate::functions::arg(&pp)), pi.clone()), raw_rel(Ne, abs(&crate::functions::arg(&pp)), pi.clone())], b(true)),
        r(LogicOp::Or, vec![raw_rel(Lt, abs(&crate::functions::arg(&pp)), pi.clone()), raw_rel(Ne, pow(&pp, &Expr::minus_one()), Expr::zero())], b(true)),
        r(LogicOp::And, vec![raw_rel(Le, abs(&argp), pi2.clone()), raw_rel(Le, abs(&sub(&argp, &pi)), pi2.clone())], eq_arg(neg(&pi2))),
        Rule { op: LogicOp::Or, args: vec![raw_rel(Le, p.clone(), q.clone()), raw_rel(Lt, p.clone(), q.clone())], to: raw_rel(Le, p.clone(), q.clone()), and_r: true },
        r(LogicOp::And, vec![raw_rel(Gt, pow(p, &int(2)), Expr::one()), raw_rel(Ne, pow(p, &int(2)), Expr::one())], raw_rel(Gt, pow(p, &int(2)), Expr::one())),
        r(LogicOp::And, vec![raw_rel(Ne, pow(p, &Expr::minus_one()), Expr::one()), raw_rel(Gt, mul2(&cosabs, &abs(p)), Expr::one())], raw_rel(Gt, abs(p), Expr::one())),
        r(LogicOp::And, vec![raw_rel(Ne, p.clone(), int(2)), raw_rel(Gt, mul2(&cosabs, &abs(p)), int(2))], raw_rel(Gt, abs(p), int(2))),
        r(LogicOp::And, vec![raw_rel(Lt, abs(&argp), pi2), raw_rel(Gt, mul2(&cosabs, &crate::ops::sqrt(&abs(&pow(p, &int(2))))), Expr::one())], raw_rel(Gt, pow(p, &int(2)), Expr::one())),
    ]
}

fn logic_parts(e: &Expr) -> Option<(LogicOp, &[Expr])> {
    match e.kind() {
        Kind::Func(Func::Logic(op), a) => Some((*op, a.as_slice())),
        _ => None,
    }
}

fn build(op: LogicOp, args: Vec<Expr>) -> Expr {
    crate::functions::func(Func::Logic(op), args)
}

/// `cond.replace(lambda _: _.is_Relational, f)`.
fn replace_rel(e: &Expr, f: &dyn Fn(&Expr) -> Expr) -> Expr {
    match e.kind() {
        Kind::Rel(..) => f(e),
        Kind::Func(Func::Logic(op), args) => build(*op, args.iter().map(|a| replace_rel(a, f)).collect()),
        _ => e.clone(),
    }
}

fn subs_map(e: &Expr, m: &rustc_hash::FxHashMap<Expr, Expr>) -> Expr {
    let pairs: Vec<(Expr, Expr)> = m.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
    crate::subs::subs(e, &pairs)
}

fn rel_touchup(rel: &Expr, p: &Expr, q: &Expr) -> Expr {
    let Kind::Rel(RelOp::Eq, lhs, rhs) = rel.kind() else { return rel.clone() };
    if !rhs.is_zero() {
        return rel.clone();
    }
    if let Kind::Func(Func::Arg, a) = lhs.kind() {
        return gt(&a[0], &Expr::zero());
    }
    if let Kind::Pow(base, _) = lhs.kind() {
        if let Kind::Func(Func::Arg, a) = base.kind() {
            return gt(&a[0], &Expr::zero());
        }
    }
    let _ = (p, q);
    if let Kind::Func(Func::Special(SpecialFn::PeriodicArgument), a) = lhs.kind() {
        if !crate::special::polar::is_polar(&a[0]) && a[1].is_constant(crate::expr::Constant::Infinity) {
            return gt(&a[0], &Expr::zero());
        }
    }
    rel.clone()
}

type Binding = rustc_hash::FxHashMap<Expr, Expr>;

/// `expr.match(pattern)` through the native port of `matches`, falling back
/// to the structural matcher only where the port declines.
fn one_match(e: &Expr, pattern: &Expr) -> Option<Binding> {
    match crate::matches::reference_match(e, pattern, false) {
        Ok(found) => found,
        Err(crate::matches::Stop::Error(_)) => None,
        Err(crate::matches::Stop::Decline(_)) => crate::pattern::pattern_match(e, pattern).ok().flatten(),
    }
}

/// `arg1.match(pattern)` in `_condsimp`: `Basic.matches` pairs the arguments
/// of two relationals of the same class in order, each one with the bindings
/// of the previous ones substituted.  The native port declines relationals
/// as a whole, so the pairing is done here (`X > 1` matches `p**2 > 1` with
/// `p = sqrt(X)`, as in the reference).
fn cond_match(e: &Expr, pattern: &Expr) -> Option<Binding> {
    let (Kind::Rel(op, l, r), Kind::Rel(pop, pl, pr)) = (e.kind(), pattern.kind()) else {
        return one_match(e, pattern);
    };
    if op != pop {
        return None;
    }
    let mut d = Binding::default();
    for (arg, other) in [(pl, l), (pr, r)] {
        if arg == other {
            continue;
        }
        let arg = subs_map(arg, &d);
        if &arg == other {
            continue;
        }
        d.extend(one_match(other, &arg)?);
    }
    Some(d)
}

/// `_condsimp(cond, first)`.
pub(crate) fn condsimp(cond: &Expr, first: bool) -> Expr {
    let cond = if first { replace_rel(cond, &crate::summation::canonical_coeff) } else { cond.clone() };
    let Some((op, args)) = logic_parts(&cond) else { return cond };
    let p = Expr::wild("p", vec![]);
    let q = Expr::wild("q", vec![]);
    let rules = rules(&p, &q);
    let mut cond = build(op, args.iter().map(|a| condsimp(a, false)).collect());
    let mut change = true;
    while change {
        change = false;
        let Some((op, args)) = logic_parts(&cond) else { break };
        let args = args.to_vec();
        'rules: for rule in &rules {
            if rule.op != op {
                continue;
            }
            for (n, arg1) in args.iter().enumerate() {
                let Some(m) = cond_match(arg1, &rule.args[0]) else { continue };
                if m.is_empty() && arg1 != &rule.args[0] {
                    continue;
                }
                let others: Vec<Expr> = rule.args[1..].iter().map(|a| subs_map(a, &m)).collect();
                let mut used = vec![n];
                for arg2 in &others {
                    for (k, arg3) in args.iter().enumerate() {
                        if used.contains(&k) {
                            continue;
                        }
                        let hit = if rule.and_r {
                            matches!(logic_parts(arg3), Some((LogicOp::And, inner)) if inner.contains(arg2))
                        } else {
                            arg2 == arg3
                        };
                        if hit {
                            used.push(k);
                            break;
                        }
                    }
                }
                if used.len() != others.len() + 1 {
                    continue;
                }
                let mut newargs: Vec<Expr> = args.iter().enumerate().filter(|(k, _)| !used.contains(k)).map(|(_, a)| a.clone()).collect();
                newargs.push(subs_map(&rule.to, &m));
                cond = build(op, newargs);
                change = true;
                break 'rules;
            }
        }
    }
    replace_rel(&cond, &|r| rel_touchup(r, &p, &q))
}

/// `_eval_cond`.
pub(crate) fn eval_cond(cond: &Expr) -> Expr {
    if matches!(cond.kind(), Kind::Bool(_)) {
        return cond.clone();
    }
    condsimp(&crate::functions::doit(cond, true), true)
}

/// `_my_unpolarify`.
pub(crate) fn my_unpolarify(f: &Expr) -> Expr {
    eval_cond(&unpolarify(f, false))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `_condsimp(And(Ne(X, 1), X > 1))` is `X > 1` in the reference: `X > 1`
    /// matches `p**2 > 1` with `p = sqrt(X)` and `Ne(p**2, 1)` becomes `Ne(X, 1)`.
    #[test]
    fn condsimp_merges_square_rule_through_relational_matching() {
        let (s, a) = (Expr::symbol("s"), Expr::symbol("a"));
        let x = add(&[re(&s), re(&a)]);
        let gt1 = raw_rel(RelOp::Gt, x.clone(), Expr::one());
        let cond = build(LogicOp::And, vec![raw_rel(RelOp::Ne, x.clone(), Expr::one()), gt1.clone()]);
        assert_eq!(condsimp(&cond, false), gt1);
        assert_eq!(cond_match(&gt1, &raw_rel(RelOp::Lt, x, Expr::one())), None);
    }
}
