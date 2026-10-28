//! Dense univariate polynomials as coefficient vectors (lowest degree
//! first, no trailing zeros; the zero polynomial is the empty vector).

use super::int::{Coeff, Field, Int, Rat};
use super::mpoly::MPoly;

pub fn trim<C: Coeff>(p: &mut Vec<C>) {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
}

#[inline]
pub fn deg<C>(p: &[C]) -> isize {
    p.len() as isize - 1
}

pub fn lc<C: Coeff>(p: &[C]) -> C {
    p.last().cloned().unwrap_or_else(C::zero)
}

pub fn add<C: Coeff>(a: &[C], b: &[C]) -> Vec<C> {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut r = long.to_vec();
    for (i, c) in short.iter().enumerate() {
        r[i].add_assign(c);
    }
    trim(&mut r);
    r
}

pub fn sub<C: Coeff>(a: &[C], b: &[C]) -> Vec<C> {
    let mut r = a.to_vec();
    if r.len() < b.len() {
        r.resize(b.len(), C::zero());
    }
    for (i, c) in b.iter().enumerate() {
        r[i] = r[i].sub(c);
    }
    trim(&mut r);
    r
}

pub fn neg<C: Coeff>(a: &[C]) -> Vec<C> {
    a.iter().map(|c| c.neg()).collect()
}

pub fn scale<C: Coeff>(a: &[C], k: &C) -> Vec<C> {
    super::factor_budget::product(a,std::slice::from_ref(k),a.len(),0,1);
    if k.is_zero() {
        return vec![];
    }
    let mut r: Vec<C> = a.iter().map(|c| c.mul(k)).collect();
    trim(&mut r);
    r
}

pub fn mul<C: Coeff>(a: &[C], b: &[C]) -> Vec<C> {
    super::factor_budget::dense(a.len().saturating_add(b.len()).saturating_sub(1));
    super::factor_budget::charge(a.len().saturating_mul(b.len()) as u64);
    super::factor_budget::product(a,b,a.len().saturating_add(b.len()).saturating_sub(1),0,a.len().min(b.len()));
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![C::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            r[i + j].add_mul(x, y);
        }
    }
    trim(&mut r);
    r
}

/// Integer polynomial product with a 128-bit accumulation fast path.
pub fn mul_int(a: &[Int], b: &[Int]) -> Vec<Int> {
    super::factor_budget::dense(a.len().saturating_add(b.len()).saturating_sub(1));
    super::factor_budget::charge(a.len().saturating_mul(b.len()) as u64);
    super::factor_budget::product(a,b,a.len().saturating_add(b.len()).saturating_sub(1),0,a.len().min(b.len()));
    if super::factor_budget::active()&&!a.is_empty()&&!b.is_empty(){
        let bits=a.iter().map(Int::bits).max().unwrap_or(0)+b.iter().map(Int::bits).max().unwrap_or(0)+64-(a.len().min(b.len()).max(1) as u64).leading_zeros() as u64;
        super::factor_budget::bits(bits);
    }
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let small = |p: &[Int]| -> Option<Vec<i64>> {
        p.iter()
            .map(|c| c.as_i64().filter(|v| v.unsigned_abs() < (1u64 << 52)))
            .collect()
    };
    if let (Some(sa), Some(sb)) = (small(a), small(b)) {
        // |a_i b_j| < 2^104, sums of up to 2^20 terms fit in i128
        if a.len().min(b.len()) < (1 << 20) {
            let mut r = vec![0i128; a.len() + b.len() - 1];
            for (i, &x) in sa.iter().enumerate() {
                if x == 0 {
                    continue;
                }
                for (j, &y) in sb.iter().enumerate() {
                    r[i + j] += x as i128 * y as i128;
                }
            }
            let mut out: Vec<Int> = r.into_iter().map(Int::from_i128).collect();
            trim(&mut out);
            return out;
        }
    }
    mul(a, b)
}

pub fn eval<C: Coeff>(p: &[C], x: &C) -> C {
    let mut acc = C::zero();
    for c in p.iter().rev() {
        acc = acc.mul(x).add(c);
    }
    acc
}

pub fn derivative<C: Coeff>(p: &[C]) -> Vec<C> {
    let mut r: Vec<C> = p
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c.mul(&C::from_i64(i as i64)))
        .collect();
    trim(&mut r);
    r
}

/// Division with remainder over a field.
pub fn divrem<C: Field>(a: &[C], b: &[C]) -> (Vec<C>, Vec<C>) {
    super::factor_budget::dense(a.len());
    super::factor_budget::charge(super::factor_budget::division_work(a.len(), b.len()));
    assert!(!b.is_empty(), "division by zero polynomial");
    if a.len() < b.len() {
        return (vec![], a.to_vec());
    }
    let mut r = a.to_vec();
    let db = b.len() - 1;
    let inv = b[db].inv();
    let mut q = vec![C::zero(); a.len() - db];
    for k in (0..q.len()).rev() {
        let c = r[k + db].mul(&inv);
        if !c.is_zero() {
            for (j, bj) in b.iter().enumerate() {
                r[k + j] = r[k + j].sub(&c.mul(bj));
            }
        }
        q[k] = c;
    }
    r.truncate(db);
    trim(&mut r);
    trim(&mut q);
    (q, r)
}

pub fn rem<C: Field>(a: &[C], b: &[C]) -> Vec<C> {
    divrem(a, b).1
}

pub fn monic<C: Field>(a: &[C]) -> Vec<C> {
    if a.is_empty() {
        return vec![];
    }
    let l = a.last().unwrap().inv();
    scale(a, &l)
}

/// Monic gcd over a field.
pub fn gcd_field<C: Field>(a: &[C], b: &[C]) -> Vec<C> {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    while !b.is_empty() {
        let r = rem(&a, &b);
        a = b;
        // Monic remainders keep the coefficient sizes of the sequence down.
        b = if r.is_empty() { r } else { monic(&r) };
    }
    monic(&a)
}

/// Extended Euclid over a field: `(g, s, t)` with `s a + t b = g` monic.
pub fn gcdex_field<C: Field>(a: &[C], b: &[C]) -> (Vec<C>, Vec<C>, Vec<C>) {
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    let (mut s0, mut s1) = (vec![C::one()], vec![]);
    let (mut t0, mut t1) = (vec![], vec![C::one()]);
    while !r1.is_empty() {
        let (q, r) = divrem(&r0, &r1);
        r0 = std::mem::replace(&mut r1, r);
        let s2 = sub(&s0, &mul(&q, &s1));
        s0 = std::mem::replace(&mut s1, s2);
        let t2 = sub(&t0, &mul(&q, &t1));
        t0 = std::mem::replace(&mut t1, t2);
    }
    if r0.is_empty() {
        return (r0, s0, t0);
    }
    let l = r0.last().unwrap().inv();
    (scale(&r0, &l), scale(&s0, &l), scale(&t0, &l))
}

// ----------------------------------------------------------------------
// integer specific
// ----------------------------------------------------------------------

pub fn content(p: &[Int]) -> Int {
    let mut g = Int::zero();
    for c in p.iter().rev() {
        g = g.gcd(c);
        if g.is_one() {
            break;
        }
    }
    g
}

/// `(content, primitive)` with a positive leading coefficient.
pub fn primitive(p: &[Int]) -> (Int, Vec<Int>) {
    if p.is_empty() {
        return (Int::zero(), vec![]);
    }
    let mut g = content(p);
    if p.last().unwrap().is_negative() {
        g = g.neg();
    }
    if g.is_one() {
        return (g, p.to_vec());
    }
    (g.clone(), p.iter().map(|c| c.div_exact_unchecked(&g)).collect())
}

/// Exact division in `Z[x]` (`None` if not exact).
pub fn div_exact_int(a: &[Int], b: &[Int]) -> Option<Vec<Int>> {
    if b.is_empty() {
        return None;
    }
    if a.is_empty() {
        return Some(vec![]);
    }
    if a.len() < b.len() {
        return None;
    }
    let mut r = a.to_vec();
    let db = b.len() - 1;
    let lb = &b[db];
    let mut q = vec![Int::zero(); a.len() - db];
    for k in (0..q.len()).rev() {
        if r[k + db].is_zero() {
            continue;
        }
        let c = r[k + db].div_exact(lb)?;
        for (j, bj) in b.iter().enumerate() {
            r[k + j].sub_mul(&c, bj);
        }
        q[k] = c;
    }
    if r.iter().any(|c| !c.is_zero()) {
        return None;
    }
    trim(&mut q);
    Some(q)
}

/// Pseudo-remainder `lc(b)^(deg a - deg b + 1) a mod b`.
pub fn prem(a: &[Int], b: &[Int]) -> Vec<Int> {
    let mut r = a.to_vec();
    if r.len() < b.len() {
        return r;
    }
    let db = b.len() - 1;
    let lb = b[db].clone();
    let mut steps = (a.len() - b.len() + 1) as u32;
    while r.len() >= b.len() {
        let k = r.len() - 1 - db;
        let lr = r.last().unwrap().clone();
        for c in r.iter_mut() {
            *c = c.mul(&lb);
        }
        for (j, bj) in b.iter().enumerate() {
            r[k + j].sub_mul(&lr, bj);
        }
        trim(&mut r);
        steps -= 1;
    }
    if steps > 0 {
        let f = lb.pow(steps);
        for c in r.iter_mut() {
            *c = c.mul(&f);
        }
    }
    r
}

pub fn max_norm(p: &[Int]) -> Int {
    p.iter().map(|c| c.abs()).max().unwrap_or_else(Int::zero)
}

pub fn to_rat(p: &[Int]) -> Vec<Rat> {
    p.iter().map(|c| Rat::from(c.clone())).collect()
}

/// `(d, d*p)` with integer coefficients.
pub fn clear_denoms(p: &[Rat]) -> (Int, Vec<Int>) {
    let mut d = Int::one();
    for c in p {
        d = d.lcm(&c.d);
    }
    (d.clone(), p.iter().map(|c| c.n.mul(&d.div_exact_unchecked(&c.d))).collect())
}

// ----------------------------------------------------------------------
// conversions with the sparse representation
// ----------------------------------------------------------------------

/// Dense coefficients of a polynomial depending only on variable `v`.
pub fn from_mpoly<C: Coeff>(p: &MPoly<C>, v: usize) -> Vec<C> {
    if p.is_zero() {
        return vec![];
    }
    let d = p.degree(v) as usize;
    super::factor_budget::dense(d.saturating_add(1));
    let mut r = vec![C::zero(); d + 1];
    for (e, c) in p.iter() {
        r[e[v] as usize].add_assign(c);
    }
    trim(&mut r);
    r
}

pub fn to_mpoly<C: Coeff>(p: &[C], v: usize, n: usize) -> MPoly<C> {
    let mut exps = Vec::with_capacity(p.len() * n);
    let mut coeffs = Vec::with_capacity(p.len());
    for (k, c) in p.iter().enumerate().rev() {
        if c.is_zero() {
            continue;
        }
        let s = exps.len();
        exps.resize(s + n, 0);
        exps[s + v] = k as u32;
        coeffs.push(c.clone());
    }
    MPoly::from_sorted(n, exps, coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zi(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&x| Int::S(x)).collect()
    }

    #[test]
    fn basics() {
        let a = zi(&[-1, 0, 1]); // x^2 - 1
        let b = zi(&[1, 1]); // x + 1
        assert_eq!(div_exact_int(&a, &b), Some(zi(&[-1, 1])));
        assert_eq!(div_exact_int(&a, &zi(&[2, 1])), None);
        assert_eq!(mul_int(&b, &b), zi(&[1, 2, 1]));
        let (q, r) = divrem(&to_rat(&a), &to_rat(&zi(&[0, 2])));
        assert_eq!(q, vec![Rat::zero(), Rat::new(Int::S(1), Int::S(2))]);
        assert_eq!(r, to_rat(&zi(&[-1])));
        let g = gcd_field(&to_rat(&a), &to_rat(&zi(&[1, 2, 1])));
        assert_eq!(g, to_rat(&b));
        assert_eq!(prem(&zi(&[1, 0, 1]), &zi(&[1, 2])), zi(&[5]));
    }
}
