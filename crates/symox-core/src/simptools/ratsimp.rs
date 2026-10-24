//! `sympy.simplify.ratsimp`: `ratsimp` and `ratsimpmodprime`.
//!
//! `ratsimpmodprime` follows the second algorithm of Monagan and Pearce
//! ("Rational Simplification Modulo a Polynomial Ideal") exactly as SymPy
//! implements it, including the particular-solution heuristic of
//! `minsolve_linear_system` used to pick the coefficients of the simplified
//! numerator and denominator.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};

use super::as_numer_denom;
use crate::expr::Expr;
use crate::ops::{add, div, mul};
use crate::poly::api::{cancel, reduced};
use crate::poly::convert::{from_exprs, from_exprs_gens, to_expr_f, PResult, PolyError};
use crate::poly::groebner::{reduce, CoefficientField, MonomialOrder};
use crate::poly::int::Rat;
use crate::poly::mpoly::MPoly;
use crate::sortkey::cmp_default;
use crate::symbol::Assumptions;

/// `ratsimp(expr)`: put over a common denominator, cancel and reduce the
/// numerator modulo the denominator.
pub fn ratsimp(e: &Expr) -> Expr {
    let c = cancel(e).unwrap_or_else(|_| e.clone());
    let (f, g) = as_numer_denom(&c);
    match reduced(&f, std::slice::from_ref(&g), &[]) {
        Ok((q, r)) => {
            let rest = div(&r, &g);
            add(&[add(&q), cancel(&rest).unwrap_or(rest)])
        }
        Err(_) => div(&f, &g),
    }
}

fn total_degree(p: &MPoly<Rat>) -> i64 {
    p.total_degree().max(0)
}

/// Dense rational vector helpers for the linear systems.
type Row = Vec<Rat>;

fn row_is_zero(r: &[Rat]) -> bool {
    r.iter().all(Rat::is_zero)
}

/// Reduced row echelon form; returns the pivot columns (one per non-zero row).
fn rref(rows: &mut Vec<Row>, ncols: usize) -> Vec<usize> {
    let mut pivots = Vec::new();
    let mut r = 0;
    for col in 0..ncols {
        if r == rows.len() {
            break;
        }
        let Some(p) = (r..rows.len()).find(|&i| !rows[i][col].is_zero()) else {
            continue;
        };
        rows.swap(r, p);
        let inv = rows[r][col].inv();
        for v in rows[r].iter_mut() {
            *v = v.mul(&inv);
        }
        for i in 0..rows.len() {
            if i != r && !rows[i][col].is_zero() {
                let f = rows[i][col].clone();
                for j in 0..ncols {
                    let t = rows[r][j].mul(&f);
                    rows[i][j] = rows[i][j].sub(&t);
                }
            }
        }
        pivots.push(col);
        r += 1;
    }
    rows.truncate(r);
    pivots
}

/// `solve_linear_system` for a homogeneous system: maps each pivot variable
/// to a linear form (dense coefficient vector) in the free variables.
fn solve_homogeneous(rows: &[Row], nvars: usize) -> BTreeMap<usize, Row> {
    let mut m: Vec<Row> = rows.iter().filter(|r| !row_is_zero(r)).cloned().collect();
    let pivots = rref(&mut m, nvars);
    let mut sol = BTreeMap::new();
    for (row, &p) in m.iter().zip(&pivots) {
        let mut form = vec![Rat::zero(); nvars];
        for j in 0..nvars {
            if j != p && !row[j].is_zero() {
                form[j] = row[j].neg();
            }
        }
        sol.insert(p, form);
    }
    sol
}

/// The linear unknowns `c0, c1, ..., d0, d1, ...` of one degree pair.
struct Unknowns {
    nc: usize,
    syms: Vec<Expr>,
}

impl Unknowns {
    fn new(nc: usize, nd: usize) -> Self {
        let mut syms = Vec::with_capacity(nc + nd);
        for i in 0..nc {
            syms.push(Expr::dummy(&format!("c{}", i), Assumptions::NONE));
        }
        for i in 0..nd {
            syms.push(Expr::dummy(&format!("d{}", i), Assumptions::NONE));
        }
        Unknowns { nc, syms }
    }
    fn name(&self, i: usize) -> String {
        if i < self.nc {
            format!("c{}", i)
        } else {
            format!("d{}", i - self.nc)
        }
    }
    fn expr(&self, form: &[Rat]) -> Expr {
        add(&form
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| mul(&[Expr::num(c.to_number()), self.syms[i].clone()]))
            .collect::<Vec<_>>())
    }
}

/// Substitute determined values into a linear form, returning the new form
/// and its constant part.
fn subs_form(form: &[Rat], konst: &Rat, det: &BTreeMap<usize, Rat>) -> (Row, Rat) {
    let mut f = form.to_vec();
    let mut k = konst.clone();
    for (&v, val) in det {
        if !f[v].is_zero() {
            k = k.add(&f[v].mul(val));
            f[v] = Rat::zero();
        }
    }
    (f, k)
}

/// `minsolve_linear_system(system, *symbols, quick=True)` for a homogeneous
/// system.  Returns `None` when there is no non-trivial pivot solution.
fn minsolve_quick(rows: &[Row], u: &Unknowns) -> Option<BTreeMap<usize, Rat>> {
    let nvars = u.syms.len();
    let s0 = solve_homogeneous(rows, nvars);
    if s0.is_empty() || s0.values().all(|f| row_is_zero(f)) {
        return None;
    }
    // pivot -> (linear part, constant)
    let mut s: BTreeMap<usize, (Row, Rat)> = s0.into_iter().map(|(k, f)| (k, (f, Rat::zero()))).collect();
    let mut determined: BTreeMap<usize, Rat> = BTreeMap::new();
    let update = |s: &mut BTreeMap<usize, (Row, Rat)>, determined: &mut BTreeMap<usize, Rat>| {
        let keys: Vec<usize> = s.keys().copied().collect();
        let mut delete = Vec::new();
        for k in keys {
            let (f, c) = s.get(&k).unwrap().clone();
            let (f, c) = subs_form(&f, &c, determined);
            if row_is_zero(&f) {
                delete.push(k);
                determined.insert(k, c.clone());
            }
            s.insert(k, (f, c));
        }
        for k in delete {
            s.remove(&k);
        }
    };
    update(&mut s, &mut determined);
    while !s.is_empty() {
        // k = max(values, key=(len(free_symbols), default_sort_key))
        let mut best: Option<(&Row, &Rat, usize, Expr)> = None;
        for (f, c) in s.values() {
            let nfree = f.iter().filter(|v| !v.is_zero()).count();
            let e = add(&[u.expr(f), Expr::num(c.to_number())]);
            let better = match &best {
                None => true,
                Some((_, _, bn, be)) => match nfree.cmp(bn) {
                    Ordering::Greater => true,
                    Ordering::Less => false,
                    Ordering::Equal => cmp_default(&e, be) == Ordering::Greater,
                },
            };
            if better {
                best = Some((f, c, nfree, e));
            }
        }
        let (f, _, nfree, _) = best.unwrap();
        let x = (0..nvars).filter(|&i| !f[i].is_zero()).max_by(|&a, &b| u.name(a).cmp(&u.name(b))).unwrap();
        let value = if nfree != 1 {
            Rat::zero()
        } else {
            let mut d0 = determined.clone();
            d0.insert(x, Rat::zero());
            let any_nonzero = s.values().any(|(f, c)| {
                let (f2, c2) = subs_form(f, c, &d0);
                !(row_is_zero(&f2) && c2.is_zero())
            });
            if any_nonzero {
                Rat::zero()
            } else {
                Rat::one()
            }
        };
        determined.insert(x, value);
        update(&mut s, &mut determined);
    }
    Some(determined)
}

/// `minsolve_linear_system(system, *symbols, quick=False)`: the exhaustive
/// search for a solution with fewer non-zero unknowns.
fn minsolve_exhaustive(rows: &[Row], u: &Unknowns) -> BTreeMap<usize, Rat> {
    let n = u.syms.len();
    let Some(mut best) = minsolve_quick(rows, u) else {
        return BTreeMap::new();
    };
    let n0 = best.values().filter(|v| !v.is_zero()).count();
    let mut k = n0 as i64 - 1;
    while k > 1 {
        let mut this: Option<BTreeMap<usize, Rat>> = None;
        for nonzeros in combinations(n, k as usize) {
            let sub: Vec<Row> = rows.iter().map(|r| nonzeros.iter().map(|&i| r[i].clone()).collect()).collect();
            let s = solve_homogeneous(&sub, nonzeros.len());
            if s.is_empty() || s.values().all(|f| row_is_zero(f)) {
                continue;
            }
            // free unknowns within `nonzeros` are set to one
            let mut sol = BTreeMap::new();
            for (&p, f) in &s {
                let mut v = Rat::zero();
                for c in f {
                    v = v.add(c);
                }
                sol.insert(nonzeros[p], v);
            }
            for i in 0..n {
                sol.entry(i).or_insert_with(|| if nonzeros.contains(&i) { Rat::one() } else { Rat::zero() });
            }
            this = Some(sol);
            break;
        }
        match this {
            None => break,
            Some(s) => best = s,
        }
        k -= 1;
    }
    best
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::with_capacity(k);
    fn rec(start: usize, n: usize, k: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        for i in start..n {
            cur.push(i);
            rec(i + 1, n, k, cur, out);
            cur.pop();
        }
    }
    rec(0, n, k, &mut cur, &mut out);
    out
}

/// One recorded minimal solution: the monomial supports of `c_hat`/`d_hat`
/// and the linear system for their coefficients.
#[derive(Clone)]
struct Recorded {
    m1: Vec<Vec<u32>>,
    m2: Vec<Vec<u32>>,
    rows: Vec<Row>,
}

struct Ctx<'a> {
    nv: usize,
    basis: &'a [MPoly<Rat>],
    lms: Vec<Vec<u32>>,
    order: MonomialOrder,
    quick: bool,
    field: CoefficientField,
    tested: HashSet<(i64, i64)>,
}

impl Ctx<'_> {
    fn nf(&self, p: &MPoly<Rat>) -> MPoly<Rat> {
        reduce(&self.field.normalize(p), self.basis, self.order, self.field).1
    }

    /// Monomials of degree `n`, `n-1`, ..., 0 not divisible by a leading
    /// monomial of the basis (SymPy `staircase`).
    fn staircase(&self, n: i64) -> Vec<Vec<u32>> {
        let mut out = Vec::new();
        let mut d = n;
        while d >= 0 {
            if d == 0 {
                out.push(vec![0; self.nv]);
                break;
            }
            for combo in combinations_with_replacement(self.nv, d as usize) {
                let mut m = vec![0u32; self.nv];
                for i in combo {
                    m[i] += 1;
                }
                if self.lms.iter().all(|l| !l.iter().zip(&m).all(|(a, b)| a <= b)) {
                    out.push(m);
                }
            }
            d -= 1;
        }
        out
    }

    fn poly_from(&self, monos: &[Vec<u32>], vals: &[Rat]) -> MPoly<Rat> {
        MPoly::from_terms(
            self.nv,
            monos.iter().zip(vals).filter(|(_, v)| !v.is_zero()).map(|(m, v)| (m.clone(), v.clone())).collect(),
        )
    }

    fn run(
        &mut self,
        a: MPoly<Rat>,
        b: MPoly<Rat>,
        mut allsol: Vec<Recorded>,
        n0: i64,
        d0: i64,
    ) -> PResult<(MPoly<Rat>, MPoly<Rat>, Vec<Recorded>)> {
        let (mut n, mut d) = (n0, d0);
        let (mut c, mut dd) = (a.clone(), b.clone());
        let mut steps = 0i64;
        let maxdeg = total_degree(&a) + total_degree(&b);
        let bound = if self.quick { maxdeg - 1 } else { maxdeg };
        while n + d <= bound {
            if !self.tested.insert((n, d)) {
                break;
            }
            let m1 = self.staircase(n);
            let m2 = self.staircase(d);
            let u = Unknowns::new(m1.len(), m2.len());
            let nvars = u.syms.len();
            // NF(a*d_hat - b*c_hat) = sum d_j NF(a*M2_j) - sum c_i NF(b*M1_i)
            let mut cols: Vec<MPoly<Rat>> = Vec::with_capacity(nvars);
            for m in &m1 {
                cols.push(self.nf(&b.mul_term(m, &Rat::one())).neg());
            }
            for m in &m2 {
                cols.push(self.nf(&a.mul_term(m, &Rat::one())));
            }
            let mut eqs: BTreeMap<Vec<u32>, Row> = BTreeMap::new();
            for (j, col) in cols.iter().enumerate() {
                for (e, cf) in col.iter() {
                    eqs.entry(e.to_vec()).or_insert_with(|| vec![Rat::zero(); nvars])[j] = cf.clone();
                }
            }
            let rows: Vec<Row> = eqs.into_values().filter(|r| !row_is_zero(r)).collect();
            let sol = if rows.is_empty() { None } else { minsolve_quick(&rows, &u) };
            if let Some(sol) = sol.filter(|s| !s.values().all(Rat::is_zero)) {
                let val = |i: usize| sol.get(&i).cloned().unwrap_or_else(Rat::one);
                let cv: Vec<Rat> = (0..m1.len()).map(val).collect();
                let dv: Vec<Rat> = (m1.len()..nvars).map(val).collect();
                c = self.poly_from(&m1, &cv);
                dd = self.poly_from(&m2, &dv);
                if dd.is_zero() {
                    return Err(PolyError::Unsupported("Ideal not prime?".into()));
                }
                allsol.push(Recorded { m1, m2, rows });
                if n + d != maxdeg {
                    allsol = vec![allsol.pop().unwrap()];
                }
                break;
            }
            steps += 1;
            n += 1;
            d += 1;
        }
        if steps > 0 {
            let (c1, d1, s1) = self.run(c, dd, allsol, n, d - steps)?;
            let (c2, d2, s2) = self.run(c1, d1, s1, n - steps, d)?;
            return Ok((c2, d2, s2));
        }
        Ok((c, dd, allsol))
    }
}

fn combinations_with_replacement(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::with_capacity(k);
    fn rec(start: usize, n: usize, k: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        for i in start..n {
            cur.push(i);
            rec(i, n, k, cur, out);
            cur.pop();
        }
    }
    rec(0, n, k, &mut cur, &mut out);
    out
}

fn is_integral(p: &MPoly<Rat>) -> bool {
    p.coeffs().iter().all(|c| c.d.is_one())
}

fn denom_lcm(p: &MPoly<Rat>) -> Rat {
    let mut l = crate::poly::int::Int::one();
    for c in p.coeffs() {
        l = l.lcm(&c.d);
    }
    Rat::new(l, crate::poly::int::Int::one())
}

fn nterms(p: &MPoly<Rat>) -> usize {
    p.len().max(1)
}

/// `ratsimpmodprime(expr, G, *gens, quick=True, polynomial=False, order=...)`.
///
/// `modulus` selects the `GF(p)` domain (as SymPy, the preliminary
/// reduction runs over QQ and only the ideal arithmetic of the main loop is
/// done modulo `p`); `force_field` is `domain=QQ`.
#[allow(clippy::too_many_arguments)]
pub fn ratsimpmodprime(
    expr: &Expr,
    g: &[Expr],
    gens: &[Expr],
    quick: bool,
    polynomial: bool,
    order: MonomialOrder,
    modulus: Option<u64>,
    force_field: bool,
) -> PResult<Expr> {
    let canc = cancel(expr).unwrap_or_else(|_| expr.clone());
    let (num, den) = as_numer_denom(&canc);
    let mut all = vec![num.clone(), den.clone()];
    all.extend_from_slice(g);
    let conv = match if gens.is_empty() { from_exprs(&all) } else { from_exprs_gens(&all, gens) } {
        Ok(c) => c,
        Err(e @ PolyError::Generator(_)) => return Err(e),
        Err(_) => return Ok(expr.clone()),
    };
    let nv = conv.ngens;
    if conv.gens.len() > nv {
        // Coefficients from a polynomial ring: only the trivial constant case.
        let constant = conv.polys[..2].iter().all(|p| p.iter().all(|(e, _)| e[..nv].iter().all(|&k| k == 0)));
        if constant {
            let q = div(&num, &den);
            return Ok(if polynomial { cancel(&q).unwrap_or(q) } else { q });
        }
        return Err(PolyError::Unsupported(
            "ratsimpmodprime over a non-rational coefficient domain is not implemented".into(),
        ));
    }
    let field = force_field || modulus.is_some() || conv.float || conv.polys.iter().any(|p| !is_integral(p));
    let to_e = |p: &MPoly<Rat>| crate::poly::convert::to_expr_prec(p, &conv.gens, conv.prec);
    let basis: Vec<MPoly<Rat>> = conv.polys[2..].to_vec();
    let lms = basis
        .iter()
        .filter(|p| !p.is_zero())
        .map(|p| p.exp(order.leading_index(p)).to_vec())
        .collect();
    let qq = CoefficientField { modulus: None };
    let a = reduce(&conv.polys[0], &basis, order, qq).1;
    let b = reduce(&conv.polys[1], &basis, order, qq).1;
    if polynomial {
        let q = div(&to_e(&a), &to_e(&b));
        return Ok(cancel(&q).unwrap_or(q));
    }
    let fld = CoefficientField { modulus };
    let basis: Vec<MPoly<Rat>> = basis.iter().map(|p| fld.normalize(p)).collect();
    let lms = if modulus.is_some() {
        basis.iter().filter(|p| !p.is_zero()).map(|p| p.exp(order.leading_index(p)).to_vec()).collect()
    } else {
        lms
    };
    let (a, b) = (fld.normalize(&a), fld.normalize(&b));
    let mut ctx = Ctx { nv, basis: &basis, lms, order, quick, field: fld, tested: HashSet::new() };
    let (mut c, mut d, allsol) = ctx.run(a, b, Vec::new(), 0, 0)?;
    if !quick && !allsol.is_empty() {
        let mut best: Option<(usize, MPoly<Rat>, MPoly<Rat>)> = None;
        for rec in &allsol {
            let u = Unknowns::new(rec.m1.len(), rec.m2.len());
            let sol = minsolve_exhaustive(&rec.rows, &u);
            let val = |i: usize| sol.get(&i).cloned().unwrap_or_else(Rat::one);
            let cv: Vec<Rat> = (0..rec.m1.len()).map(val).collect();
            let dv: Vec<Rat> = (rec.m1.len()..u.syms.len()).map(val).collect();
            let cc = ctx.poly_from(&rec.m1, &cv);
            let dd = ctx.poly_from(&rec.m2, &dv);
            let size = nterms(&cc) + nterms(&dd);
            if best.as_ref().map_or(true, |(s, _, _)| size < *s) {
                best = Some((size, cc, dd));
            }
        }
        let (_, cc, dd) = best.unwrap();
        c = cc;
        d = dd;
    }
    if !field {
        let cn = denom_lcm(&c);
        let dn = denom_lcm(&d);
        c = c.scale(&cn);
        d = d.scale(&dn);
        let r = cn.div(&dn);
        c = c.scale(&Rat::new(r.d.clone(), crate::poly::int::Int::one()));
        d = d.scale(&Rat::new(r.n.clone(), crate::poly::int::Int::one()));
    }
    Ok(div(&to_e(&c), &to_e(&d)))
}
