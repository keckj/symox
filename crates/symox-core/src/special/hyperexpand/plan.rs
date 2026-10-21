//! Order reduction, origin construction and shift plans (SymPy's
//! `reduce_order`, `try_shifted_sum`, `try_polynomial`, `try_lerchphi`,
//! `build_hypergeometric_formula`, `devise_plan`, `devise_plan_meijer`) and
//! the Slater sums of `_meijergexpand`.

use crate::functions as f;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::special::{factorial, hyper, lerchphi, rf};
use crate::{Expr, Kind};

use super::ops::{apply_operators, Op};
use super::{
    bucket, hyperexpand_func, is_integer_number, num_cmp, polar_lift, sift_mod1, sort_default, sort_num, Buckets, Formula, HResult,
    HxError, Rewrite, GF, HF,
};

fn rel_err() -> HxError { HxError::Type("cannot determine truth value of Relational".into()) }

fn cmp0(e: &Expr) -> HResult<std::cmp::Ordering> { num_cmp(e).ok_or_else(rel_err) }

fn le0(e: &Expr) -> HResult<bool> { Ok(cmp0(e)? != std::cmp::Ordering::Greater) }
fn gt0(e: &Expr) -> HResult<bool> { Ok(cmp0(e)? == std::cmp::Ordering::Greater) }

fn small(e: &Expr) -> HResult<i64> { e.as_i64().filter(|_| is_integer_number(e)).ok_or_else(|| HxError::Type(format!("'{e}' object cannot be interpreted as an integer"))) }

fn remove_first(v: &mut Vec<Expr>, x: &Expr) -> HResult<()> {
    match v.iter().position(|e| e == x) {
        Some(i) => {
            v.remove(i);
            Ok(())
        }
        None => Err(HxError::Value("list.remove(x): x not in list".into())),
    }
}

fn prod(v: impl IntoIterator<Item = Expr>) -> Expr { mul(&v.into_iter().collect::<Vec<_>>()) }

/// `_reduce_order(ap, bq, gen, key)`.
fn reduce_with(
    ap: &[Expr],
    bq: &[Expr],
    gen: &dyn Fn(&Expr, &Expr) -> Option<Op>,
    key: &dyn Fn(&Expr, &Expr) -> std::cmp::Ordering,
) -> (Vec<Expr>, Vec<Expr>, Vec<Op>) {
    let mut ap = ap.to_vec();
    let mut bq = bq.to_vec();
    ap.sort_by(key);
    bq.sort_by(key);
    let mut nap = Vec::new();
    let mut ops = Vec::new();
    for a in ap {
        let mut found = None;
        for i in 0..bq.len() {
            if let Some(op) = gen(&a, &bq[i]) {
                bq.remove(i);
                found = Some(op);
                break;
            }
        }
        match found {
            Some(op) => ops.push(op),
            None => nap.push(a),
        }
    }
    (nap, bq, ops)
}

/// `reduce_order(func)`.
pub(super) fn reduce_order(func: &HF) -> (HF, Vec<Op>) {
    let (nap, nbq, ops) = reduce_with(&func.ap, &func.bq, &Op::reduce_order, &crate::sortkey::cmp_default);
    (HF::new(nap, nbq), ops)
}

/// `reduce_order_meijer(func)`.
pub(super) fn reduce_order_meijer(func: &GF) -> (GF, Vec<Op>) {
    let negkey = |a: &Expr, b: &Expr| crate::sortkey::cmp_default(&neg(a), &neg(b));
    let (nan, nbq, mut ops1) = reduce_with(&func.an, &func.bq, &Op::meijer_plus, &negkey);
    let (nbm, nap, ops2) = reduce_with(&func.bm, &func.ap, &Op::meijer_minus, &crate::sortkey::cmp_default);
    ops1.extend(ops2);
    (GF::new(nan, nap, nbm, nbq), ops1)
}

/// `do_shifts(fro, to, inc, dec)`: shift `fro` onto `to` in place.
fn do_shifts(
    fro: &mut Vec<Expr>,
    to: &[Expr],
    inc: &dyn Fn(&[Expr], usize) -> HResult<Op>,
    dec: &dyn Fn(&[Expr], usize) -> HResult<Op>,
) -> HResult<Vec<Op>> {
    let mut ops = Vec::new();
    for i in 0..fro.len() {
        let up = gt0(&sub(&to[i], &fro[i]))?;
        let mut guard = 0usize;
        while to[i] != fro[i] {
            guard += 1;
            if guard > 100_000 {
                return Err(HxError::Value("shift plan does not terminate".into()));
            }
            ops.push(if up { inc(fro, i)? } else { dec(fro, i)? });
            fro[i] = add(&[fro[i].clone(), Expr::int(if up { 1 } else { -1 })]);
        }
    }
    Ok(ops)
}

fn cat(a: &[Expr], b: &[Expr]) -> Vec<Expr> { a.iter().chain(b).cloned().collect() }

fn get(b: &Buckets, k: &Expr) -> Vec<Expr> { bucket(b, k).cloned().unwrap_or_default() }

fn set(b: &mut Buckets, k: &Expr, v: Vec<Expr>) {
    match b.iter_mut().find(|(kk, _)| kk == k) {
        Some((_, slot)) => *slot = v,
        None => b.push((k.clone(), v)),
    }
}

fn others(b: &Buckets, key: &Expr) -> Vec<Expr> { b.iter().filter(|(k, _)| k != key).flat_map(|(_, v)| v.iter().cloned()).collect() }

/// `devise_plan(target, origin, z)`.
pub(super) fn devise_plan(target: &HF, origin: &HF, z: &Expr) -> HResult<Vec<Op>> {
    let (abuckets, bbuckets) = (sift_mod1(&target.ap), sift_mod1(&target.bq));
    let (mut nabuckets, mut nbbuckets) = (sift_mod1(&origin.ap), sift_mod1(&origin.bq));
    if abuckets.len() != nabuckets.len() || bbuckets.len() != nbbuckets.len() {
        return Err(HxError::Value("not reachable".into()));
    }
    let mut keys: Vec<Expr> = abuckets.iter().chain(&bbuckets).map(|(k, _)| k.clone()).collect();
    sort_default(&mut keys);
    let mut ops = Vec::new();
    for r in &keys {
        let (mut al, mut nal, mut bk, mut nbk) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        if bucket(&abuckets, r).is_some() {
            al = get(&abuckets, r);
            nal = get(&nabuckets, r);
            if bucket(&nabuckets, r).is_none() {
                nabuckets.push((r.clone(), Vec::new()));
            }
        }
        if bucket(&bbuckets, r).is_some() {
            bk = get(&bbuckets, r);
            nbk = get(&nbbuckets, r);
            if bucket(&nbbuckets, r).is_none() {
                nbbuckets.push((r.clone(), Vec::new()));
            }
        }
        if al.len() != nal.len() || bk.len() != nbk.len() {
            return Err(HxError::Value("not reachable".into()));
        }
        for w in [&mut al, &mut nal, &mut bk, &mut nbk] {
            sort_default(w);
        }
        let aother = others(&nabuckets, r);
        let bother = others(&nbbuckets, r);
        let shifts_a = |nal: &mut Vec<Expr>, nbk: &[Expr], al: &[Expr]| -> HResult<Vec<Op>> {
            let nbkb = cat(nbk, &bother);
            do_shifts(nal, al, &|p, i| Op::shift_a(&p[i]), &|p, i| Op::unshift_a(&cat(p, &aother), &nbkb, i, z))
        };
        let shifts_b = |nal: &[Expr], nbk: &mut Vec<Expr>, bk: &[Expr]| -> HResult<Vec<Op>> {
            let nala = cat(nal, &aother);
            do_shifts(nbk, bk, &|p, i| Op::unshift_b(&nala, &cat(p, &bother), i, z), &|p, i| Op::shift_b(&p[i]))
        };
        if al.is_empty() {
            ops.extend(shifts_b(&[], &mut nbk, &bk)?);
        } else if bk.is_empty() {
            ops.extend(shifts_a(&mut nal, &[], &al)?);
        } else {
            let namax = nal.last().unwrap().clone();
            let amax = al.last().unwrap().clone();
            if le0(&sub(&nbk[0], &namax))? || le0(&sub(&bk[0], &amax))? {
                return Err(HxError::Value("Non-suitable parameters.".into()));
            }
            if gt0(&sub(&namax, &amax))? {
                ops.extend(shifts_a(&mut nal, &nbk.clone(), &al)?);
                ops.extend(shifts_b(&al, &mut nbk, &bk)?);
            } else {
                ops.extend(shifts_b(&nal.clone(), &mut nbk, &bk)?);
                ops.extend(shifts_a(&mut nal, &bk, &al)?);
            }
        }
        set(&mut nabuckets, r, al);
        set(&mut nbbuckets, r, bk);
    }
    ops.reverse();
    Ok(ops)
}

/// `try_shifted_sum(func, z)`.
pub(super) fn try_shifted_sum(func: &HF, z: &Expr) -> HResult<Option<(HF, Vec<Op>, Expr)>> {
    let (ab, bb) = (sift_mod1(&func.ap), sift_mod1(&func.bq));
    let a0 = get(&ab, &Expr::zero());
    if a0.len() != 1 {
        return Ok(None);
    }
    let r = a0[0].clone();
    if le0(&r)? {
        return Ok(None);
    }
    let Some(l) = bucket(&bb, &Expr::zero()) else { return Ok(None) };
    let mut l = l.clone();
    sort_num(&mut l)?;
    let k = l[0].clone();
    if le0(&k)? {
        return Ok(None);
    }
    let mut nap = func.ap.clone();
    remove_first(&mut nap, &r)?;
    let mut nbq = func.bq.clone();
    remove_first(&mut nbq, &k)?;
    let k = sub(&k, &Expr::one());
    let nap: Vec<Expr> = nap.iter().map(|x| sub(x, &k)).collect();
    let nbq: Vec<Expr> = nbq.iter().map(|x| sub(x, &k)).collect();
    let mut ops = Vec::new();
    for n in 0..small(&r)? - 1 {
        ops.push(Op::shift_a(&Expr::int(n + 1))?);
    }
    ops.reverse();
    let mut fac = div(&factorial(&k), &pow(z, &k));
    fac = mul(&[fac, prod(nbq.iter().map(|b| rf(b, &k)))]);
    fac = div(&fac, &prod(nap.iter().map(|a| rf(a, &k))));
    ops.push(Op::mult(&fac));
    let mut p = Expr::zero();
    for n in 0..small(&k)? {
        let n = Expr::int(n);
        let mut m = div(&pow(z, &n), &factorial(&n));
        m = mul(&[m, prod(nap.iter().map(|a| rf(a, &n)))]);
        m = div(&m, &prod(nbq.iter().map(|b| rf(b, &n))));
        p = add(&[p, m]);
    }
    Ok(Some((HF::new(nap, nbq), ops, neg(&p))))
}

/// `try_polynomial(func, z)`.
pub(super) fn try_polynomial(func: &HF, z: &Expr) -> HResult<Option<Expr>> {
    let (ab, bb) = (sift_mod1(&func.ap), sift_mod1(&func.bq));
    let mut a0 = get(&ab, &Expr::zero());
    let mut b0 = get(&bb, &Expr::zero());
    sort_num(&mut a0)?;
    sort_num(&mut b0)?;
    let mut al0 = Vec::new();
    for x in a0 {
        if le0(&x)? {
            al0.push(x);
        }
    }
    let mut bl0 = Vec::new();
    for x in b0 {
        if le0(&x)? {
            bl0.push(x);
        }
    }
    if let Some(last) = bl0.last() {
        let mut all = true;
        for a in &al0 {
            if cmp0(&sub(a, last))? != std::cmp::Ordering::Less {
                all = false;
                break;
            }
        }
        if all {
            return Ok(Some(Expr::oo()));
        }
    }
    let Some(a) = al0.last() else { return Ok(None) };
    let mut fac = Expr::one();
    let mut res = Expr::one();
    for n in 0..small(&neg(a))? {
        let n = Expr::int(n);
        fac = mul(&[fac, z.clone()]);
        fac = div(&fac, &add(&[n.clone(), Expr::one()]));
        fac = mul(&[fac, prod(func.ap.iter().map(|a| add(&[a.clone(), n.clone()])))]);
        fac = div(&fac, &prod(func.bq.iter().map(|b| add(&[b.clone(), n.clone()]))));
        res = add(&[res, fac.clone()]);
    }
    Ok(Some(res))
}

/// `try_lerchphi(func)`.
pub(super) fn try_lerchphi(func: &HF) -> HResult<Option<Formula>> {
    let ab = sift_mod1(&func.ap);
    let mut bb = sift_mod1(&func.bq);
    let mut paired: Vec<(Expr, Vec<Expr>, Vec<Expr>)> = Vec::new();
    for (key, value) in &ab {
        if !key.is_zero() && bucket(&bb, key).is_none() {
            return Ok(None);
        }
        let bvalue = get(&bb, key);
        paired.push((key.clone(), value.clone(), bvalue));
        bb.retain(|(k, _)| k != key);
    }
    if !bb.is_empty() {
        return Ok(None);
    }
    let Some(zi) = paired.iter().position(|(k, _, _)| k.is_zero()) else { return Ok(None) };
    paired[zi].2.push(Expr::one());
    let t = Expr::dummy("t", Default::default());
    let mut numer = Expr::one();
    let mut denom = Expr::one();
    for (_, avalue, bvalue) in &paired {
        if avalue.len() != bvalue.len() {
            return Ok(None);
        }
        for (a, b) in avalue.iter().zip(bvalue) {
            let d = sub(a, b);
            if crate::special::is_positive(&d) == Some(true) {
                numer = mul(&[numer, rf(&add(&[b.clone(), t.clone()]), &d)]);
                denom = mul(&[denom, rf(b, &d)]);
            } else {
                let k = sub(b, a);
                numer = mul(&[numer, rf(a, &k)]);
                denom = mul(&[denom, rf(&add(&[a.clone(), t.clone()]), &k)]);
            }
        }
    }
    let part = crate::poly::api::apart(&div(&numer, &denom), Some(&t)).map_err(|e| HxError::Polynomial(format!("{e:?}")))?;
    let mut monomials: Vec<(Expr, i64)> = Vec::new();
    let mut terms: Vec<(Expr, Vec<(Expr, i64)>)> = Vec::new();
    for arg in crate::simptools::add_args(&part) {
        let (numer, denom) = crate::simptools::as_numer_denom(&arg);
        if !crate::traverse::has(&denom, &t) {
            let p = crate::poly::polyobj::Poly::from_exprs(&[numer.clone()], std::slice::from_ref(&t))
                .map_err(|e| HxError::Polynomial(format!("{e:?}")))?;
            let tms = p[0].terms();
            if tms.len() > 1 {
                return Err(HxError::Type("p should be monomial".into()));
            }
            let (deg, a) = tms.into_iter().next().map(|(m, c)| (m[0] as i64, c)).unwrap_or((0, Expr::zero()));
            monomials.push((div(&a, &denom), deg));
            continue;
        }
        if crate::traverse::has(&numer, &t) {
            return Err(HxError::NotImplemented("Need partial fraction decomposition with linear denominators".into()));
        }
        let (mut indep, deps) = crate::decompose::as_coeff_mul(&denom, std::slice::from_ref(&t), false);
        if deps.len() != 1 {
            return Err(HxError::Value("not enough values to unpack".into()));
        }
        let mut dep = deps[0].clone();
        let mut n = 1i64;
        if let Kind::Pow(b, e) = dep.kind() {
            n = small(e)?;
            dep = b.clone();
        }
        let a;
        if dep == t {
            a = Expr::zero();
        } else if dep.is_add() {
            let (a0, tmp) = crate::decompose::as_independent(&dep, std::slice::from_ref(&t), None);
            let mut b = Expr::one();
            if tmp != t {
                b = crate::decompose::as_independent(&tmp, std::slice::from_ref(&t), None).0;
            }
            if dep != add(&[mul(&[b.clone(), t.clone()]), a0.clone()]) {
                return Err(HxError::NotImplemented(format!("unrecognised form {dep}")));
            }
            a = div(&a0, &b);
            indep = mul(&[indep, pow(&b, &Expr::int(n))]);
        } else {
            return Err(HxError::NotImplemented("unrecognised form of partial fraction".into()));
        }
        let entry = (div(&numer, &indep), n);
        match terms.iter_mut().find(|(k, _)| *k == a) {
            Some((_, l)) => l.push(entry),
            None => terms.push((a, vec![entry])),
        }
    }
    let z = Expr::dummy("z", Default::default());
    monomials.sort_by_key(|m| m.1);
    let mut mon = vec![div(&Expr::one(), &sub(&Expr::one(), &z))];
    if let Some(last) = monomials.last() {
        for k in 0..last.1 {
            let next = mul(&[z.clone(), crate::diff::diff(&mon[k as usize], &z)]);
            mon.push(next);
        }
    }
    let mut coeffs: Vec<(Expr, Vec<Expr>)> = Vec::new();
    let push = |coeffs: &mut Vec<(Expr, Vec<Expr>)>, k: Expr, v: Expr| match coeffs.iter_mut().find(|(kk, _)| *kk == k) {
        Some((_, l)) => l.push(v),
        None => coeffs.push((k, vec![v])),
    };
    for (a, n) in &monomials {
        push(&mut coeffs, Expr::one(), mul(&[a.clone(), mon[*n as usize].clone()]));
    }
    let mut deriv: Vec<(Expr, Vec<(Expr, Expr)>)> = Vec::new();
    let set_deriv = |deriv: &mut Vec<(Expr, Vec<(Expr, Expr)>)>, k: Expr, v: Vec<(Expr, Expr)>| match deriv.iter_mut().find(|(kk, _)| *kk == k) {
        Some((_, l)) => *l = v,
        None => deriv.push((k, v)),
    };
    for (a, l) in terms.iter_mut() {
        for (c, k) in l.iter() {
            push(&mut coeffs, lerchphi(&z, &Expr::int(*k), a), c.clone());
        }
        l.sort_by_key(|x| x.1);
        let top = l.last().unwrap().1;
        for k in 2..=top {
            let lk = lerchphi(&z, &Expr::int(k), a);
            set_deriv(&mut deriv, lk.clone(), vec![(neg(a), lk), (Expr::one(), lerchphi(&z, &Expr::int(k - 1), a))]);
        }
        let l1 = lerchphi(&z, &Expr::one(), a);
        set_deriv(&mut deriv, l1.clone(), vec![(neg(a), l1), (div(&Expr::one(), &sub(&Expr::one(), &z)), Expr::one())]);
    }
    let mut trans: Vec<Expr> = vec![Expr::one()];
    for (k, _) in &deriv {
        match trans.iter().position(|x| x == k) {
            Some(_) => {}
            None => trans.push(k.clone()),
        }
    }
    let idx = |b: &Expr| trans.iter().position(|x| x == b);
    let basis: Vec<Expr> = trans.iter().map(crate::special::expand_func).collect();
    let n = basis.len();
    let mut c = vec![Expr::zero(); n];
    for (b, l) in &coeffs {
        let i = idx(b).ok_or_else(|| HxError::Value(format!("KeyError: {b}")))?;
        c[i] = add(l);
    }
    let mut m = vec![vec![Expr::zero(); n]; n];
    for (b, l) in &deriv {
        let i = idx(b).ok_or_else(|| HxError::Value(format!("KeyError: {b}")))?;
        for (cc, b2) in l {
            let j = idx(b2).ok_or_else(|| HxError::Value(format!("KeyError: {b2}")))?;
            m[i][j] = cc.clone();
        }
    }
    Ok(Some(Formula { func: func.clone(), z, symbols: vec![], b: basis, c, m }))
}

/// `build_hypergeometric_formula(func)`.
pub(super) fn build_hypergeometric_formula(func: &HF) -> HResult<Formula> {
    let z = Expr::dummy("z", Default::default());
    if !func.ap.is_empty() {
        // x*prod(x + b - 1) - z*prod(x + a), coefficients low to high.
        let x = Expr::dummy("x", Default::default());
        let expr = sub(
            &mul(&[x.clone(), prod(func.bq.iter().map(|b| add(&[x.clone(), b.clone(), Expr::minus_one()])))]),
            &mul(&[z.clone(), prod(func.ap.iter().map(|a| add(&[x.clone(), a.clone()])))]),
        );
        let p = crate::poly::polyobj::Poly::from_exprs(&[expr], std::slice::from_ref(&x)).map_err(|e| HxError::Polynomial(format!("{e:?}")))?;
        let all: Vec<Expr> =
            p[0].all_coeffs().map_err(|e| HxError::Polynomial(format!("{e:?}")))?.iter().map(super::ops::norm_coeff).collect();
        let lc = all[0].clone();
        let n = all.len() - 1;
        let mut basis = Vec::new();
        let mut m = vec![vec![Expr::zero(); n]; n];
        for k in 0..n {
            let a = add(&[func.ap[0].clone(), Expr::int(k as i64)]);
            let mut ap = vec![a.clone()];
            ap.extend(func.ap[1..].iter().cloned());
            basis.push(hyper(&ap, &func.bq, &z));
            if k + 1 < n {
                m[k][k] = neg(&a);
                m[k][k + 1] = a;
            }
        }
        let mut c = vec![Expr::zero(); n];
        c[0] = Expr::one();
        let eye: Vec<Vec<Expr>> = (0..n).map(|i| (0..n).map(|j| if i == j { Expr::one() } else { Expr::zero() }).collect()).collect();
        let mut derivs = vec![eye];
        for k in 0..n {
            derivs.push(matmul(&m, &derivs[k]));
        }
        let mut l = all.clone();
        l.reverse();
        let mut res = vec![Expr::zero(); n];
        for (k, cf) in l.iter().enumerate() {
            // C*derivs[k] is the first row of derivs[k].
            for (r, d) in derivs[k][0].iter().enumerate() {
                res[r] = add(&[res[r].clone(), mul(&[cf.clone(), d.clone()])]);
            }
        }
        let piv = derivs[n - 1][0][n - 1].clone();
        for (k, cf) in res.iter().enumerate() {
            m[n - 1][k] = div(&div(&neg(cf), &piv), &lc);
        }
        return Ok(Formula { func: func.clone(), z, symbols: vec![], b: basis, c, m });
    }
    let mut basis = Vec::new();
    let mut bq = func.bq.clone();
    for i in 0..bq.len() {
        basis.push(hyper(&[], &bq, &z));
        bq[i] = add(&[bq[i].clone(), Expr::one()]);
    }
    basis.push(hyper(&[], &bq, &z));
    let n = basis.len();
    let mut c = vec![Expr::zero(); n];
    c[0] = Expr::one();
    let mut m = vec![vec![Expr::zero(); n]; n];
    m[0][n - 1] = div(&z, &prod(func.bq.iter().cloned()));
    for k in 1..n {
        m[k][k - 1] = func.bq[k - 1].clone();
        m[k][k] = neg(&func.bq[k - 1]);
    }
    Ok(Formula { func: func.clone(), z, symbols: vec![], b: basis, c, m })
}

fn matmul(a: &[Vec<Expr>], b: &[Vec<Expr>]) -> Vec<Vec<Expr>> {
    let n = a.len();
    let p = b.first().map_or(0, Vec::len);
    (0..n)
        .map(|i| (0..p).map(|j| add(&(0..b.len()).map(|k| mul(&[a[i][k].clone(), b[k][j].clone()])).collect::<Vec<_>>())).collect())
        .collect()
}

/// `hyperexpand_special(ap, bq, z)`.
pub(super) fn hyperexpand_special(ap: &[Expr], bq: &[Expr], z_: &Expr) -> Expr {
    let z = crate::special::polar::unpolarify(z_, false);
    if z.is_zero() {
        return Expr::one();
    }
    if ap.len() == 2 && bq.len() == 1 {
        let (mut a, mut b, c) = (ap[0].clone(), ap[1].clone(), bq[0].clone());
        if z.is_one() {
            return div(
                &div(&mul(&[f::gamma(&sub(&sub(&c, &a), &b)), f::gamma(&c)]), &f::gamma(&sub(&c, &a))),
                &f::gamma(&sub(&c, &b)),
            );
        }
        let simp = |e: &Expr| crate::simplify::simplify(e);
        if z.is_minus_one() && simp(&add(&[sub(&b, &a), c.clone()])).is_one() {
            std::mem::swap(&mut a, &mut b);
        }
        if z.is_minus_one() && simp(&add(&[sub(&a, &b), c.clone()])).is_one() {
            let half = super::half();
            let b2 = mul(&[b.clone(), half.clone()]);
            if crate::special::is_integer(&b) == Some(true) && crate::special::is_negative(&b) == Some(true) {
                return div(
                    &div(
                        &mul(&[Expr::int(2), f::cos(&mul(&[Expr::pi(), b2.clone()])), f::gamma(&neg(&b)), f::gamma(&add(&[sub(&b, &a), Expr::one()]))]),
                        &f::gamma(&neg(&b2)),
                    ),
                    &f::gamma(&add(&[sub(&b2, &a), Expr::one()])),
                );
            }
            return div(
                &div(&mul(&[f::gamma(&add(&[b2.clone(), Expr::one()])), f::gamma(&add(&[sub(&b, &a), Expr::one()]))]), &f::gamma(&add(&[b.clone(), Expr::one()]))),
                &f::gamma(&add(&[sub(&b2, &a), Expr::one()])),
            );
        }
    }
    hyper(ap, bq, z_)
}

/// `devise_plan_meijer(fro, to, z)`.
pub(super) fn devise_plan_meijer(fro: &GF, to: &GF, z: &Expr) -> HResult<Vec<Op>> {
    let (mut fan, mut fap, mut fbm, mut fbq) = (fro.an.clone(), fro.ap.clone(), fro.bm.clone(), fro.bq.clone());
    let mut ops = Vec::new();
    // Which list, direction, and shifter kind each step uses, in order.
    #[derive(Clone, Copy)]
    enum L {
        An,
        Ap,
        Bm,
        Bq,
    }
    let steps: [(L, i64, u8); 8] =
        [(L::An, 1, 0), (L::Ap, 1, 1), (L::Bm, -1, 2), (L::Bq, -1, 3), (L::An, -1, 4), (L::Ap, -1, 5), (L::Bm, 1, 6), (L::Bq, 1, 7)];
    let mut guard = 0usize;
    'outer: loop {
        guard += 1;
        if guard > 100_000 {
            return Err(HxError::NotImplemented("Could not devise plan.".into()));
        }
        for (list, diff, kind) in steps {
            let (f, t): (&Vec<Expr>, &Vec<Expr>) = match list {
                L::An => (&fan, &to.an),
                L::Ap => (&fap, &to.ap),
                L::Bm => (&fbm, &to.bm),
                L::Bq => (&fbq, &to.bq),
            };
            let counter: Vec<Expr> = match kind {
                0 | 1 => cat(&fbm, &fbq),
                2 | 3 => cat(&fan, &fap),
                _ => vec![],
            };
            let mut hit = None;
            for (idx, (a, b)) in f.iter().zip(t).enumerate() {
                let d = sub(a, b);
                if crate::special::is_integer(&d) == Some(true)
                    && crate::special::is_positive(&div(&neg(&d), &Expr::int(diff))) == Some(true)
                    && counter.iter().all(|x| a != x)
                {
                    hit = Some(idx);
                    break;
                }
            }
            if let Some(i) = hit {
                let op = match kind {
                    0 => Op::meijer_unshift_b(&fan, &fap, &fbm, &fbq, i, z)?,
                    1 => Op::meijer_unshift_d(&fan, &fap, &fbm, &fbq, i, z)?,
                    2 => Op::meijer_unshift_a(&fan, &fap, &fbm, &fbq, i, z)?,
                    3 => Op::meijer_unshift_c(&fan, &fap, &fbm, &fbq, i, z)?,
                    4 => Op::meijer_shift_b(&fan[i]),
                    5 => Op::meijer_shift_d(&fap[i]),
                    6 => Op::meijer_shift_a(&fbm[i]),
                    _ => Op::meijer_shift_c(&fbq[i]),
                };
                let target = match list {
                    L::An => &mut fan,
                    L::Ap => &mut fap,
                    L::Bm => &mut fbm,
                    L::Bq => &mut fbq,
                };
                target[i] = add(&[target[i].clone(), Expr::int(diff)]);
                ops.push(op);
                continue 'outer;
            }
        }
        break;
    }
    if fan != to.an || fap != to.ap || fbm != to.bm || fbq != to.bq {
        return Err(HxError::NotImplemented("Could not devise plan.".into()));
    }
    ops.reverse();
    Ok(ops)
}

/// `do_slater(an, bm, ap, bq, z, zfinal)` of `_meijergexpand`; the
/// condition is a Python bool (a `Bool` node) or a relational.
#[allow(clippy::too_many_arguments)]
pub(super) fn do_slater(
    an: &[Expr],
    bm: &[Expr],
    ap: &[Expr],
    bq: &[Expr],
    z: &Expr,
    zfinal: &Expr,
    ops: &[Op],
    t: &Expr,
) -> HResult<(Expr, Expr)> {
    let func = GF::new(an.to_vec(), bm.to_vec(), ap.to_vec(), bq.to_vec());
    let [_, pbm, pap, _] = func.compute_buckets()?;
    for (i, items) in &pbm {
        if items.len() > 1 {
            let l = bucket(&pap, i).map_or(0, Vec::len);
            if l + 1 < items.len() {
                return Ok((Expr::zero(), Expr::boolean(false)));
            }
        }
    }
    let mut cond = Expr::boolean(an.len() + ap.len() < bm.len() + bq.len());
    if an.len() + ap.len() == bm.len() + bq.len() {
        cond = crate::relational::lt(&f::abs(z), &Expr::one());
    }
    if matches!(cond.kind(), Kind::Bool(false)) {
        return Ok((Expr::zero(), cond));
    }
    let mut res = Expr::zero();
    let neg1 = |n: i64| if n.rem_euclid(2) == 0 { Expr::one() } else { Expr::minus_one() };
    for (m, items) in &pbm {
        if items.len() == 1 {
            let bh = items[0].clone();
            let mut fac = Expr::one();
            let mut bo = bm.to_vec();
            remove_first(&mut bo, &bh)?;
            for bj in &bo {
                fac = mul(&[fac, f::gamma(&sub(bj, &bh))]);
            }
            for aj in an {
                fac = mul(&[fac, f::gamma(&sub(&add(&[Expr::one(), bh.clone()]), aj))]);
            }
            for bj in bq {
                fac = div(&fac, &f::gamma(&sub(&add(&[Expr::one(), bh.clone()]), bj)));
            }
            for aj in ap {
                fac = div(&fac, &f::gamma(&sub(aj, &bh)));
            }
            let nap: Vec<Expr> = an.iter().chain(ap).map(|a| sub(&add(&[Expr::one(), bh.clone()]), a)).collect();
            let nbq: Vec<Expr> = bo.iter().chain(bq).map(|b| sub(&add(&[Expr::one(), bh.clone()]), b)).collect();
            let k = polar_lift(&neg1(ap.len() as i64 - bm.len() as i64));
            let harg = mul(&[k.clone(), zfinal.clone()]);
            let premult = pow(&div(t, &k), &bh);
            let hyp = hyperexpand_func(HF::new(nap, nbq), &harg, ops, t, &premult, &bh, Rewrite::None)?;
            res = add(&[res, mul(&[fac, hyp])]);
        } else {
            let b_ = items[0].clone();
            let ki: Vec<Expr> = items[1..].iter().map(|bi| sub(bi, &b_)).collect();
            let u = ki.len();
            let pm = bucket(&pap, m).cloned().unwrap_or_default();
            let li: Vec<Expr> = pm.iter().take(u + 1).map(|ai| sub(ai, &b_)).collect();
            let mut bo = bm.to_vec();
            for b in items {
                remove_first(&mut bo, b)?;
            }
            let mut ao = ap.to_vec();
            for a in pm.iter().take(u) {
                remove_first(&mut ao, a)?;
            }
            let lu = li.last().cloned().ok_or_else(|| HxError::Value("list index out of range".into()))?;
            let di: Vec<Expr> = li.iter().zip(&ki).map(|(l, k)| sub(l, k)).collect();
            let s = Expr::dummy("s", Default::default());
            let mut integrand = pow(z, &s);
            for b in bm {
                let b = match b.as_number() {
                    Some(n) if n.is_float() && crate::functions::floor(b) == *b => Expr::int(small(&crate::functions::floor(b))?),
                    _ => b.clone(),
                };
                integrand = mul(&[integrand, f::gamma(&sub(&b, &s))]);
            }
            for a in an {
                integrand = mul(&[integrand, f::gamma(&add(&[sub(&Expr::one(), a), s.clone()]))]);
            }
            for b in bq {
                integrand = div(&integrand, &f::gamma(&add(&[sub(&Expr::one(), b), s.clone()])));
            }
            for a in ap {
                integrand = div(&integrand, &f::gamma(&sub(a, &s)));
            }
            let integrand = crate::special::expand_func(&integrand);
            for r in 0..small(&lu)? {
                let point = add(&[b_.clone(), Expr::int(r)]);
                let resid = crate::calculus::residue(&integrand, &s, &point).map_err(HxError::NotImplemented)?;
                let resid = apply_operators(resid, ops, &mut |e: &Expr| Ok(mul(&[z.clone(), crate::diff::diff(e, z)])))?;
                res = sub(&res, &resid);
            }
            let au = add(&[b_.clone(), lu.clone()]);
            let k = polar_lift(&neg1((ao.len() + bo.len() + 1) as i64));
            let harg = mul(&[k.clone(), zfinal.clone()]);
            let premult = pow(&div(t, &k), &au);
            let mut nap: Vec<Expr> = an.iter().chain(ap).map(|a| sub(&add(&[Expr::one(), au.clone()]), a)).collect();
            nap.push(Expr::one());
            let nbq: Vec<Expr> = bm.iter().chain(bq).map(|b| sub(&add(&[Expr::one(), au.clone()]), b)).collect();
            let hyp = hyperexpand_func(HF::new(nap, nbq), &harg, ops, t, &premult, &au, Rewrite::None)?;
            let mut c = div(&pow(&Expr::minus_one(), &lu), &factorial(&lu));
            for i in 0..u.min(di.len()) {
                c = mul(&[c, div(&pow(&Expr::minus_one(), &di[i]), &rf(&add(&[sub(&lu, &li[i]), Expr::one()]), &di[i]))]);
            }
            for a in an {
                c = mul(&[c, f::gamma(&add(&[sub(&Expr::one(), a), au.clone()]))]);
            }
            for b in &bo {
                c = mul(&[c, f::gamma(&sub(b, &au))]);
            }
            for a in &ao {
                c = div(&c, &f::gamma(&sub(a, &au)));
            }
            for b in bq {
                c = div(&c, &f::gamma(&add(&[sub(&Expr::one(), b), au.clone()])));
            }
            res = add(&[res, mul(&[c, hyp])]);
        }
    }
    Ok((res, cond))
}

