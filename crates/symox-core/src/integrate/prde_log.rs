//! Base-case logarithmic derivative kernels of SymPy `prde.py`:
//! `is_log_deriv_k_t_radical_in_field`, `parametric_log_deriv_heu` and
//! `limited_integrate`, with `D = [Poly(1, x)]` and ZZ/QQ coefficients.
//!
//! Every kernel returns `None` when the reference would raise or reach a
//! branch that is not ported (the Python body then runs unchanged), and
//! `Some(None)` for a proven `None` answer of the reference.
use crate::expr::Expr;
use crate::integrate::prde::{param_rischde, ParamRisch, P};
use crate::ops;
use crate::poly::convert::{from_exprs_gens, rat_to_expr, to_expr};
use crate::poly::mpoly::MPoly;
use crate::poly::{dense, Int, Rat};

fn monic(c: &[Rat]) -> Vec<Rat> {
    if c.is_empty() {
        return Vec::new();
    }
    let l = dense::lc(c).inv();
    c.iter().map(|x| x.mul(&l)).collect()
}

/// `p.as_expr()` of a polynomial in `x`.
fn as_expr(c: &[Rat], x: &Expr) -> Expr {
    to_expr(&dense::to_mpoly(c, 0, 1), std::slice::from_ref(x))
}

/// Primitive integral associate of a rational coefficient vector.
fn integral(c: &[Rat]) -> Vec<Int> {
    let mut l = Int::one();
    for v in c {
        l = l.lcm(&v.d);
    }
    let l = Rat::from(l.abs());
    c.iter().map(|v| v.mul(&l).n).collect()
}

/// The rational roots of a squarefree `s`, ascending, when `s` splits into
/// linear factors over QQ (`real_roots` then lists `deg(s)` rationals).
pub(crate) fn rational_split(s: &[Rat]) -> Option<Vec<Rat>> {
    let (_, factors) = crate::poly::factor_uni::factor_zx(&integral(s));
    let mut roots = Vec::new();
    for (f, k) in &factors {
        if f.len() != 2 {
            return None;
        }
        let r = Rat::new(f[0].neg(), f[1].clone());
        for _ in 0..*k {
            roots.push(r.clone());
        }
    }
    if roots.len() + 1 != s.len() {
        return None;
    }
    roots.sort();
    Some(roots)
}

/// The value at `z = rho` of a bivariate `S(x, z)`, as a polynomial in `x`.
fn eval_z(s: &MPoly<Rat>, rho: &Rat) -> Vec<Rat> {
    let mut out: Vec<Rat> = vec![Rat::zero(); (s.degree(0).max(0) + 1) as usize];
    for (e, c) in s.iter() {
        let mut v = c.clone();
        for _ in 0..e[1] {
            v = v.mul(rho);
        }
        out[e[0] as usize] = out[e[0] as usize].add(&v);
    }
    dense::trim(&mut out);
    out
}

/// One `(s, h)` pair of the base-case `residue_reduce(a, d, DE, z)`.
enum Log {
    /// `h = d` (multiplicity `deg d`).
    Denominator,
    /// The subresultant of that degree, normalized by `invert`.
    Subresultant(MPoly<Rat>),
}

/// Base-case `residue_reduce(fa, fd, DE, z)` (always `b == True` there):
/// the monic `s` and the `h` of every pair of `H` with a non-constant `s`,
/// and the number of pairs of `H` (the constant `(1, h)` pair that
/// `splitfactor_sqf` keeps for a non-unit integer content of `Poly(r, z)`
/// when no factor has multiplicity 1). `None` when the reference raises
/// (`NotInvertible`) or the PRS exceeds its budget.
#[allow(clippy::type_complexity)]
fn residue_reduce(fa: &P, fd: &P) -> Option<(Vec<Rat>, Vec<(Vec<Rat>, Log)>, usize)> {
    let (a, d) = fa.cancel(fd)?;
    let l = d.lc().inv();
    let d: Vec<Rat> = d.c.iter().map(|v| v.mul(&l)).collect();
    if a.is_zero() {
        return Some((d, Vec::new(), 0));
    }
    let a = dense::rem(&a.c.iter().map(|v| v.mul(&l)).collect::<Vec<_>>(), &d);
    let n = d.len() as i64 - 1;
    if n < 1 {
        return None;
    }
    let dd = dense::derivative(&d);
    let big_d = MPoly::from_terms(2, d.iter().enumerate().filter(|(_, c)| !c.is_zero())
        .map(|(k, c)| (vec![k as u32, 0], c.clone())).collect());
    let mut q: Vec<(Vec<u32>, Rat)> = a.iter().enumerate().filter(|(_, c)| !c.is_zero())
        .map(|(k, c)| (vec![k as u32, 0], c.clone())).collect();
    q.extend(dd.iter().enumerate().filter(|(_, c)| !c.is_zero()).map(|(k, c)| (vec![k as u32, 1], c.neg())));
    let big_q = MPoly::from_terms(2, q);
    let prs = crate::poly::prs::subresultants(&big_d, &big_q, 0).ok()?;
    let r = dense::from_mpoly(&prs.resultant, 1);
    let mut h = Vec::new();
    if r.is_empty() {
        return Some((d, h, 0));
    }
    let (_, factors) = crate::poly::sqf::sqf_list_zx(&integral(&r));
    // `Poly(r, z).sqf_list_include()` starts with the content `(c, 1)` when
    // no factor has multiplicity 1; over ZZ (integral coefficients) a
    // non-unit `c` is a pair of `H` (h = d for `deg d == 1`, otherwise the
    // degree-1 subresultant when there is one).
    let mut extra = 0;
    if factors.iter().all(|(f, i)| f.len() < 2 || *i != 1) && r.iter().all(|v| v.d.is_one()) {
        let mut c = Int::zero();
        for v in &r {
            c = c.gcd(&v.n);
        }
        if !c.abs().is_one() && (n == 1 || prs.sequence.iter().any(|p| p.degree(0) == 1)) {
            extra = 1;
        }
    }
    for (f, i) in factors {
        if f.len() < 2 {
            continue;
        }
        let s = monic(&f.iter().map(|v| Rat::from(v.clone())).collect::<Vec<_>>());
        if i as i64 == n {
            h.push((s, Log::Denominator));
            continue;
        }
        let Some(sub) = prs.sequence.iter().rev().find(|p| p.degree(0) == i as i64) else { continue };
        let lc = dense::from_mpoly(&sub.lc_in(0), 1);
        if dense::gcd_field(&lc, &s).len() != 1 {
            return None;
        }
        h.push((s, Log::Subresultant(sub.clone())));
    }
    let count = h.len() + extra;
    Some((d, h, count))
}

fn is_sqf(c: &[Rat]) -> bool {
    c.len() <= 1 || dense::gcd_field(c, &dense::derivative(c)).len() == 1
}

/// Base-case `is_log_deriv_k_t_radical_in_field(fa, fd, DE)`: `(n, u)`.
pub fn is_log_deriv_k_t_radical_in_field(fa: &P, fd: &P, x: &Expr) -> Option<Option<(Int, Expr)>> {
    is_log_deriv_traced(fa, fd, x, &mut Vec::new())
}

/// [`is_log_deriv_k_t_radical_in_field`] that also records, in order, the
/// Dummies the stock body creates: `z`, one unnamed Dummy per pair of `H`
/// (`real_roots`) and `i` (`residue_reduce_derivation`) when every `s`
/// splits.
pub fn is_log_deriv_traced(fa: &P, fd: &P, x: &Expr, trace: &mut Vec<&'static str>) -> Option<Option<(Int, Expr)>> {
    if fd.is_zero() {
        return None;
    }
    let (fa, fd) = fa.cancel(fd)?;
    if fd.deg() == 0 {
        // residue_reduce has nothing to do; `p` is the constant `f` and the
        // base branch keeps only `f == 0` (the empty product).
        trace.extend(["z", "i"]);
        return Some(if fa.is_zero() { Some((Int::one(), Expr::one())) } else { None });
    }
    let (_, pairs, count) = residue_reduce(&fa, &fd)?;
    trace.push("z");
    trace.extend(std::iter::repeat_n("", count));
    let mut roots = Vec::with_capacity(pairs.len());
    for (s, _) in &pairs {
        match rational_split(s) {
            Some(r) => roots.push(r),
            None => return Some(None),
        }
    }
    trace.push("i");
    // The `p in k` check and the base-case conditions: `p` is a constant
    // only if every residue is covered, and the answer needs `fd`
    // squarefree with `deg(fa) < deg(fd)`.
    if !is_sqf(&fd.c) || fa.deg() >= fd.deg() {
        return Some(None);
    }
    let dmonic = monic(&fd.c);
    let mut terms: Vec<(Vec<Rat>, Rat)> = Vec::new();
    for ((_, log), rs) in pairs.iter().zip(&roots) {
        for rho in rs {
            let h = match log {
                Log::Denominator => dmonic.clone(),
                Log::Subresultant(sub) => monic(&eval_z(sub, rho)),
            };
            terms.push((h, rho.clone()));
        }
    }
    // p = f - sum(rho*Dh/h), reduced.
    let (mut num, mut den) = (fa.c.clone(), fd.c.clone());
    for (h, rho) in &terms {
        let dh: Vec<Rat> = dense::derivative(h).iter().map(|v| v.mul(rho)).collect();
        num = dense::sub(&dense::mul(&num, h), &dense::mul(&dh, &den));
        den = dense::mul(&den, h);
        let g = dense::gcd_field(&num, &den);
        if g.len() > 1 {
            num = dense::divrem(&num, &g).0;
            den = dense::divrem(&den, &g).0;
        }
    }
    if num.len() > 1 && den.len() > 1 {
        let g = dense::gcd_field(&num, &den);
        num = dense::divrem(&num, &g).0;
        den = dense::divrem(&den, &g).0;
    }
    if den.len() != 1 || num.len() > 1 {
        return Some(None);
    }
    let mut n = Int::one();
    for (_, rho) in &terms {
        n = n.lcm(&rho.d).abs();
    }
    let nr = Rat::from(n.clone());
    let factors: Vec<Expr> = terms.iter()
        .map(|(h, rho)| ops::pow(&as_expr(h, x), &rat_to_expr(&rho.mul(&nr))))
        .collect();
    Some(Some((n, ops::mul(&factors))))
}

/// `solve([p[i] - c*q[i] ...], c)` for rational pairs: the unique solution,
/// `None` when the system has none or is not determined (`[]`).
fn solve_ratio(eqs: &[(Rat, Rat)]) -> Option<Rat> {
    let c = eqs.iter().find(|(_, q)| !q.is_zero()).map(|(p, q)| p.div(q))?;
    eqs.iter().all(|(p, q)| p.sub(&c.mul(q)).is_zero()).then_some(c)
}

/// Rational `lcm` of the polys module for two numbers (`QQ.lcm`).
fn lcm_number(a: &Rat, b: &Rat) -> Rat {
    let p = a.n.lcm(&b.n).abs();
    let q = a.d.gcd(&b.d).abs();
    Rat::new(p, q)
}

/// Result of the base-case `parametric_log_deriv_heu`: whether
/// `is_log_deriv_k_t_radical_in_field` was reached (it creates a Dummy),
/// and the answer `(n, m, v)` or `None`.
pub struct Heu {
    pub reached: bool,
    pub answer: Option<(Int, Int, Expr)>,
    /// The Dummies created by the stock body, in order (`c1`, then those of
    /// `is_log_deriv_k_t_radical_in_field`).
    pub trace: Vec<&'static str>,
}

fn finish(c: Rat, fa: &P, fd: &P, wa: &P, wd: &P, x: &Expr) -> Option<Heu> {
    let (m, n) = (Rat::from(c.n.clone()), Rat::from(c.d.clone()));
    let nfmwa = P::inferred(vec![n]).mul(fa).mul(wd).sub(&P::inferred(vec![m]).mul(wa).mul(fd));
    let nfmwd = fd.mul(wd);
    let mut trace = vec!["c1"];
    let answer = match is_log_deriv_traced(&nfmwa, &nfmwd, x, &mut trace)? {
        None => None,
        Some((q, v)) => {
            // `Q.is_zero or v.is_zero`: Q >= 1; v is a product of powers.
            if crate::assumptions::is_zero(&v) == Some(true) {
                None
            } else {
                Some((q.mul(&c.d), q.mul(&c.n), v))
            }
        }
    };
    Some(Heu { reached: true, answer, trace })
}

/// Base-case `parametric_log_deriv_heu(fa, fd, wa, wd, DE)` (`B = 0`).
pub fn parametric_log_deriv_heu(fa: &P, fd: &P, wa: &P, wd: &P, x: &Expr) -> Option<Heu> {
    let none = Some(Heu { reached: false, answer: None, trace: vec!["c1"] });
    let (p, _) = fa.div(fd)?;
    let (q, _) = wa.div(wd)?;
    if q.deg() > 0 {
        let top = p.deg().max(q.deg());
        let eqs: Vec<(Rat, Rat)> = (1..=top).map(|i| (p.nth(i), q.nth(i))).collect();
        let Some(c) = solve_ratio(&eqs) else { return none };
        return finish(c, fa, fd, wa, wd, x);
    }
    if p.deg() > 0 {
        return none;
    }
    let c = lcm_number(&fd.lc(), &wd.lc());
    let fm = P::new(monic(&fd.c), true);
    let wm = P::new(monic(&wd.c), true);
    let l = fm.lcm(&wm)?.mul(&P::inferred(vec![c]));
    // Base-case splitfactor(l) is (l, 1): z = gcd(l, l').
    let z = l.gcd(&l.diff());
    let (_, r1) = fa.mul(&l.quo(fd)?).div(&z)?;
    let (_, r2) = wa.mul(&l.quo(wd)?).div(&z)?;
    let eqs: Vec<(Rat, Rat)> = (0..z.deg().max(0)).map(|i| (r1.nth(i), r2.nth(i))).collect();
    let Some(c) = solve_ratio(&eqs) else { return none };
    finish(c, fa, fd, wa, wd, x)
}

/// `Poly(expr, x)` for a polynomial expression with rational coefficients.
fn poly_of(e: &Expr, x: &Expr) -> Option<P> {
    let conv = from_exprs_gens(std::slice::from_ref(e), std::slice::from_ref(x)).ok()?;
    if conv.float || conv.gens.len() != 1 || conv.ngens != 1 {
        return None;
    }
    Some(P::inferred(dense::from_mpoly(&conv.polys[0], 0)))
}

/// Base-case `limited_integrate(fa, fd, G, DE)`: `((Ya, Yd), C)` with the
/// entries of `C` (published as `Poly(c, x, domain=QQ)`).
#[allow(clippy::type_complexity)]
pub fn limited_integrate(fa: &P, fd: &P, g: &[(P, P)], x: &Expr) -> Option<Option<((P, P), Vec<Rat>)>> {
    if fd.is_zero() {
        return None;
    }
    let fa = fa.mul(&P::inferred(vec![fd.lc().inv()]));
    let fd = P::new(monic(&fd.c), true);
    let mut gg = Vec::with_capacity(g.len() + 1);
    gg.push((fa, fd));
    gg.extend(g.iter().cloned());
    let (h, a) = match param_rischde(&P::zero(false), &P::one(false), &gg)? {
        ParamRisch::NoSolution(_) => return Some(None),
        ParamRisch::Solution(h, a) => (h, a),
    };
    let Some(v0) = a.nullspace().into_iter().find(|v| !v[0].is_zero()) else { return Some(None) };
    let r = h.len();
    if r == 0 || v0.len() < r + 1 {
        // `-sum([])` is the int 0, which has no `as_numer_denom`.
        return None;
    }
    let c0 = v0[0].neg();
    let v: Vec<Rat> = v0.iter().map(|e| e.div(&c0)).collect();
    let m = v.len() - r - 1;
    let cs = v[1..=m].to_vec();
    let mut sum = Expr::zero();
    for (i, (ha, hd)) in h.iter().enumerate() {
        let vi = &v[m + 1 + i];
        let prod: Vec<Rat> = if vi.is_zero() { Vec::new() } else { ha.c.iter().map(|e| e.mul(vi)).collect() };
        let term = ops::div(&as_expr(&prod, x), &as_expr(&hd.c, x));
        sum = ops::add2(&sum, &term);
    }
    let y = ops::neg(&sum);
    let (yn, yd) = crate::simptools::as_numer_denom(&y);
    let (ya, yd) = (poly_of(&yn, x)?, poly_of(&yd, x)?);
    if yd.is_zero() {
        return None;
    }
    let ya = ya.mul(&P::inferred(vec![yd.lc().inv()]));
    let yd = P::new(monic(&yd.c), true);
    Some(Some(((ya, yd), cs)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::int::Coeff;
    fn z(c: &[i64]) -> P {
        P::new(c.iter().map(|&v| Rat::from_i64(v)).collect(), false)
    }
    fn x() -> Expr {
        Expr::symbol("x")
    }
    fn show(e: &Expr) -> String {
        format!("{}", e)
    }
    #[test]
    fn log_deriv_base_cases() {
        let x = x();
        let (n, u) = is_log_deriv_k_t_radical_in_field(&z(&[5]), &z(&[0, 2]), &x).unwrap().unwrap();
        assert_eq!((n, show(&u)), (Int::from_i64(2), "x**5".to_string()));
        let (n, u) = is_log_deriv_k_t_radical_in_field(&z(&[1, 3]), &z(&[0, 1, 1]), &x).unwrap().unwrap();
        assert_eq!((n, show(&u)), (Int::one(), "x*(x + 1)**2".to_string()));
        assert!(is_log_deriv_k_t_radical_in_field(&z(&[1]), &z(&[-2, 0, 1]), &x).unwrap().is_none());
        assert!(is_log_deriv_k_t_radical_in_field(&z(&[0, 1]), &z(&[1]), &x).unwrap().is_none());
        assert!(is_log_deriv_k_t_radical_in_field(&z(&[1]), &z(&[0, 0, 1]), &x).unwrap().is_none());
        let (n, u) = is_log_deriv_k_t_radical_in_field(&z(&[]), &z(&[1]), &x).unwrap().unwrap();
        assert_eq!((n, u), (Int::one(), Expr::one()));
    }
    #[test]
    fn limited_integrate_reference_examples() {
        let x = x();
        let g = vec![(z(&[0, 1]), z(&[1, 1]))];
        let ((ya, yd), c) = limited_integrate(&z(&[-1, -1, -5, 3]), &z(&[1, -1, -1, 1]), &g, &x).unwrap().unwrap();
        assert_eq!(ya, z(&[2, -1, 1]));
        assert_eq!(yd, P::new(vec![Rat::from_i64(-1), Rat::one()], true));
        assert_eq!(c, vec![Rat::from_i64(2)]);
        let g = vec![(z(&[1]), z(&[0, 1]))];
        let ((ya, yd), c) = limited_integrate(&z(&[0, 0, 5]), &z(&[3]), &g, &x).unwrap().unwrap();
        let five_ninths = Rat::new(Int::from_i64(5), Int::from_i64(9));
        assert_eq!(ya, P::new(vec![Rat::zero(), Rat::zero(), Rat::zero(), five_ninths], true));
        assert_eq!(yd, P::new(vec![Rat::one()], true));
        assert_eq!(c, vec![Rat::zero()]);
        // param_poly_rischDE returns no polynomial (h == []): C = eye(m), None.
        let q = |n: i64, d: i64| Rat::new(Int::from_i64(n), Int::from_i64(d));
        let fa = P::new(vec![q(-2, 3), q(2, 3), q(-1, 3)], true);
        let g = vec![(P::new(vec![q(1, 2)], true), z(&[-1, 0, 4])), (z(&[4]), z(&[1, 1]))];
        assert!(limited_integrate(&fa, &z(&[0, 0, 2]), &g, &x).unwrap().is_none());
    }
    #[test]
    fn heu_constant_case() {
        // f = 5/(2x), w = 1: c = 0 is not determined by the empty system.
        let x = x();
        let h = parametric_log_deriv_heu(&z(&[5]), &z(&[0, 2]), &z(&[1]), &z(&[1]), &x).unwrap();
        assert!(!h.reached && h.answer.is_none());
    }
}
