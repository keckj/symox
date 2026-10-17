//! Power residues and discrete logarithms, ported from SymPy's
//! `ntheory.residue_ntheory` (same root orders and error messages).

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::factor::{perfect_power, remove, trailing};
use super::mt::PyRandom;
use super::{factorint, isprime, jacobi_big, mod_inverse};

type R<T> = Result<T, String>;

fn b(n: i64) -> BigInt {
    BigInt::from(n)
}

fn range_step(start: BigInt, stop: &BigInt, step: &BigInt) -> Vec<BigInt> {
    let mut out = Vec::new();
    let mut x = start;
    while &x < stop {
        out.push(x.clone());
        x += step;
    }
    out
}

fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    a.lcm(b)
}

/// `pow(x, e, m)` with Python semantics for negative exponents.
fn pow_mod(x: &BigInt, e: &BigInt, m: &BigInt) -> Option<BigInt> {
    if m.is_one() {
        return Some(BigInt::zero());
    }
    if e.is_negative() {
        let inv = mod_inverse(&x.mod_floor(m), m)?;
        return Some(inv.modpow(&-e, m));
    }
    Some(x.mod_floor(m).modpow(e, m))
}

/// Jacobi symbol `(a/n)` for odd positive `n` (SymPy's `jacobi`).
pub fn jacobi(a: &BigInt, n: &BigInt) -> i32 {
    jacobi_big(a, n)
}

/// Kronecker symbol `(a/n)`.
pub fn kronecker(a: &BigInt, n: &BigInt) -> i32 {
    if n.is_zero() {
        return if a.abs().is_one() { 1 } else { 0 };
    }
    let mut n = n.clone();
    let mut sign = 1;
    if n.is_negative() {
        n = -n;
        if a.is_negative() {
            sign = -1;
        }
    }
    let t = trailing(&n);
    if t > 0 {
        if a.is_even() {
            return 0;
        }
        let r = a.mod_floor(&b(8));
        if t % 2 == 1 && (r == b(3) || r == b(5)) {
            sign = -sign;
        }
        n >>= t;
    }
    sign * jacobi_big(a, &n)
}

/// Multiplicative order of `a` modulo `n`.
pub fn n_order(a: &BigInt, n: &BigInt) -> R<BigInt> {
    if n <= &BigInt::one() {
        return Err("n should be an integer greater than 1".into());
    }
    let a = a.mod_floor(n);
    if a.is_one() {
        return Ok(BigInt::one());
    }
    if !a.gcd(n).is_one() {
        return Err("The two numbers should be relatively prime".into());
    }
    let mut a_order = BigInt::one();
    for (p, e) in factorint(n) {
        let pe = p.pow(e);
        let pe_order: BigInt = (&p - 1) * p.pow(e - 1);
        let mut factors: Vec<(BigInt, u32)> = factorint(&(&p - 1)).into_iter().collect();
        if e > 1 {
            factors.push((p.clone(), e - 1));
        }
        let mut order = BigInt::one();
        for (px, ex) in factors {
            let mut x = a.modpow(&(&pe_order / px.pow(ex)), &pe);
            while !x.is_one() {
                x = x.modpow(&px, &pe);
                order *= &px;
            }
        }
        a_order = lcm(&a_order, &order);
    }
    Ok(a_order)
}

/// Smallest primitive root of an odd prime `p` (first value of SymPy's
/// `_primitive_root_prime_iter`).
fn first_root_prime(p: &BigInt) -> BigInt {
    if p == &b(3) {
        return b(2);
    }
    let r8 = (p % 8u32).to_u32().unwrap();
    let g_min = if r8 == 1 || r8 == 7 { 3 } else { 2 };
    if p < &b(41) {
        return if p == &b(23) { b(5) } else { b(g_min) };
    }
    let pm1: BigInt = p - 1;
    let v: Vec<BigInt> = factorint(&pm1).keys().map(|q| &pm1 / q).collect();
    let mut g = b(g_min);
    while &g < p {
        if v.iter().all(|pw| !g.modpow(pw, p).is_one()) {
            return g;
        }
        g += 1;
    }
    g
}

/// First value of `_primitive_root_prime_power_iter(p, e)`.
fn first_root_prime_power(p: &BigInt, e: u32) -> BigInt {
    let g = first_root_prime(p);
    if e == 1 {
        return g;
    }
    let p2 = p * p;
    let t = (&g - pow_mod(&g, &(2 - p), &p2).unwrap()).mod_floor(&p2);
    if !t.is_zero() {
        g
    } else {
        g + p
    }
}

/// Is `a` a primitive root modulo `p`?
pub fn is_primitive_root(a: &BigInt, p: &BigInt) -> R<bool> {
    if p <= &BigInt::one() {
        return Err("p should be an integer greater than 1".into());
    }
    let a = a.mod_floor(p);
    if !a.gcd(p).is_one() {
        return Err("The two numbers should be relatively prime".into());
    }
    if p <= &b(4) {
        return Ok(a == p - 1);
    }
    let q = if p.is_odd() {
        p.clone()
    } else if !(p % 4u32).is_zero() {
        p / 2
    } else {
        return Ok(false);
    };
    let (group_order, factors): (BigInt, Vec<BigInt>) = if isprime(&q) {
        (&q - 1, factorint(&(&q - 1)).into_keys().collect())
    } else {
        let (q, e) = match perfect_power(&q, None, true) {
            Some(m) => m,
            None => return Ok(false),
        };
        if !isprime(&q) {
            return Ok(false);
        }
        let mut fs: Vec<BigInt> = factorint(&(&q - 1)).into_keys().collect();
        if !fs.contains(&q) {
            fs.push(q.clone());
        }
        (q.pow((e - 1) as u32) * (&q - 1), fs)
    };
    Ok(factors.iter().all(|pr| !a.modpow(&(&group_order / pr), p).is_one()))
}

/// Primitive root modulo `p` (smallest one, or SymPy's fast choice).
pub fn primitive_root(p: &BigInt, smallest: bool) -> R<Option<BigInt>> {
    if p <= &BigInt::one() {
        return Err("p should be an integer greater than 1".into());
    }
    if p <= &b(4) {
        return Ok(Some(p - 1));
    }
    let p_even = p.is_even();
    let q = if !p_even {
        p.clone()
    } else if !(p % 4u32).is_zero() {
        p / 2
    } else {
        return Ok(None);
    };
    let (q, e) = if isprime(&q) {
        (q, 1u32)
    } else {
        match perfect_power(&q, None, true) {
            Some((q, e)) if isprime(&q) => (q, e as u32),
            _ => return Ok(None),
        }
    };
    if !smallest {
        let g = first_root_prime_power(&q, e);
        if p_even {
            return Ok(Some(if g.is_odd() { g } else { g + q.pow(e) }));
        }
        return Ok(Some(g));
    }
    if p_even {
        let mut i = b(3);
        while &i < p {
            if !(&i % &q).is_zero() && is_primitive_root(&i, p)? {
                return Ok(Some(i));
            }
            i += 2;
        }
        return Ok(None);
    }
    let g = first_root_prime(&q);
    if e == 1 || !g.modpow(&(&q - 1), &(&q * &q)).is_one() {
        return Ok(Some(g));
    }
    let mut i: BigInt = &g + 1u32;
    while &i < p {
        if !(&i % &q).is_zero() && is_primitive_root(&i, p)? {
            return Ok(Some(i));
        }
        i += 1;
    }
    Ok(None)
}

fn tonelli_shanks(a: &BigInt, p: &BigInt) -> BigInt {
    let pm1: BigInt = p - 1;
    let s = trailing(&pm1);
    let t: BigInt = p >> s;
    let d = if (p % 12u32) == b(5) {
        b(3)
    } else if matches!((p % 5u32).to_u32(), Some(2 | 3)) {
        b(5)
    } else {
        // any quadratic non-residue gives a valid square root; the
        // caller sorts the pair of roots
        let mut d = b(6);
        while jacobi_big(&d, p) != -1 {
            d += 1;
        }
        d
    };
    let big_a = a.modpow(&t, p);
    let big_d = d.modpow(&t, p);
    let mut m = BigInt::zero();
    for i in 0..s {
        let adm = (&big_a * big_d.modpow(&m, p)) % p;
        let adm = adm.modpow(&(BigInt::one() << (s - 1 - i)), p);
        if adm == pm1 {
            m += BigInt::one() << i;
        }
    }
    (a.modpow(&((&t + 1) / 2), p) * big_d.modpow(&(&m / 2), p)) % p
}

fn sqrt_mod_prime_power(a: &BigInt, p: &BigInt, k: u32) -> Option<Vec<BigInt>> {
    let pk = p.pow(k);
    let a = a.mod_floor(&pk);
    if p == &b(2) {
        if (&a % 8u32) != BigInt::one() {
            return None;
        }
        if k <= 3 {
            return Some(range_step(b(1), &pk, &b(2)));
        }
        let mut r = BigInt::one();
        for nx in 3..k {
            if ((&r * &r - &a) >> nx).is_odd() {
                r += BigInt::one() << (nx - 1);
            }
        }
        let h = BigInt::one() << (k - 1);
        let mut v = vec![r.clone(), &pk - &r, (&r + &h).mod_floor(&pk), (-(&r + &h)).mod_floor(&pk)];
        v.sort();
        return Some(v);
    }
    if jacobi_big(&a, p) != 1 {
        return None;
    }
    let mut res = if (p % 4u32) == b(3) {
        a.modpow(&((p + 1) / 4), p)
    } else if (p % 8u32) == b(5) {
        let mut res = a.modpow(&((p + 3) / 8), p);
        if res.modpow(&b(2), p) != (&a % p) {
            res = res * b(2).modpow(&((p - 1) / 4), p) % p;
        }
        res
    } else {
        tonelli_shanks(&a, p)
    };
    if k > 1 {
        let mut px = p.clone();
        let kbits = 32 - k.leading_zeros();
        for _ in 0..kbits - 1 {
            px = &px * &px;
            let frinv = mod_inverse(&(&res * 2), &px)?;
            res = (&res - (&res * &res - &a) * frinv).mod_floor(&px);
        }
        if k & (k - 1) != 0 {
            let frinv = mod_inverse(&(&res * 2), &pk)?;
            res = (&res - (&res * &res - &a) * frinv).mod_floor(&pk);
        }
    }
    let mut v = vec![res.clone(), &pk - &res];
    v.sort();
    Some(v)
}

fn sqrt_mod1(a: &BigInt, p: &BigInt, n: u32) -> Option<Vec<BigInt>> {
    let pn = p.pow(n);
    let a = a.mod_floor(&pn);
    if a.is_zero() {
        let step = p.pow(n.div_ceil(2));
        return Some(range_step(BigInt::zero(), &pn, &step));
    }
    let (a, r) = remove(&a, p)?;
    if r % 2 == 1 {
        return None;
    }
    let res = sqrt_mod_prime_power(&a, p, n - r as u32)?;
    let m = (r / 2) as u32;
    let step = p.pow(n - m);
    let pm = p.pow(m);
    let mut out = Vec::new();
    for rx in res {
        let mut x = rx * &pm;
        while x < pn {
            out.push(x.clone());
            x += &step;
        }
    }
    Some(out)
}

/// SymPy's `_iproduct2` ordering of a finite Cartesian product.
fn iproduct2<T: Clone>(l1: &[T], l2: &[Vec<T>]) -> Vec<Vec<T>> {
    let (n1, n2) = (l1.len(), l2.len());
    let mut out = Vec::new();
    if n1 == 0 || n2 == 0 {
        return out;
    }
    let mut n = 0usize;
    loop {
        let len1 = (n + 1).min(n1);
        let len2 = (n + 1).min(n2);
        if n > len1 + len2 {
            break;
        }
        let lo = (n + 1).saturating_sub(len1);
        for m in lo..len2 {
            let mut t = vec![l1[n - m].clone()];
            t.extend(l2[m].iter().cloned());
            out.push(t);
        }
        n += 1;
        // once both lists are exhausted, len1 + len2 stays fixed
        if n > n1 + n2 + 1 {
            break;
        }
    }
    out
}

/// SymPy's `iproduct` ordering for finite lists.
fn iproduct<T: Clone>(ls: &[Vec<T>]) -> Vec<Vec<T>> {
    match ls.len() {
        0 => vec![vec![]],
        1 => ls[0].iter().map(|e| vec![e.clone()]).collect(),
        _ => {
            let rest = iproduct(&ls[1..]);
            iproduct2(&ls[0], &rest)
        }
    }
}

/// `itertools.product` ordering.
fn product<T: Clone>(ls: &[Vec<T>]) -> Vec<Vec<T>> {
    let mut out: Vec<Vec<T>> = vec![vec![]];
    for l in ls {
        let mut next = Vec::with_capacity(out.len() * l.len());
        for o in &out {
            for e in l {
                let mut t = o.clone();
                t.push(e.clone());
                next.push(t);
            }
        }
        out = next;
    }
    out
}

/// Garner style CRT of SymPy's `gf_crt1`/`gf_crt2`: the unique value in
/// `[0, prod(ms))`.
fn crt_coprime(us: &[BigInt], ms: &[BigInt]) -> BigInt {
    let p: BigInt = ms.iter().product();
    let mut v = BigInt::zero();
    for (u, m) in us.iter().zip(ms) {
        let e = &p / m;
        let s = mod_inverse(&e, m).unwrap_or_else(BigInt::zero);
        v += &e * (u * s).mod_floor(m);
    }
    v.mod_floor(&p)
}

/// Square roots of `a` modulo `|p|` in SymPy's `sqrt_mod_iter` order.
pub fn sqrt_mod_iter(a: &BigInt, p: &BigInt) -> Vec<BigInt> {
    let p = p.abs();
    let mut v: Vec<Vec<BigInt>> = Vec::new();
    let mut pv = Vec::new();
    let mut use_iproduct = false;
    for (px, ex) in factorint(&p) {
        let rx = if !(a % &px).is_zero() {
            sqrt_mod_prime_power(a, &px, ex)
        } else {
            use_iproduct = true;
            sqrt_mod1(a, &px, ex)
        };
        match rx {
            Some(r) if !r.is_empty() => v.push(r),
            _ => return vec![],
        }
        pv.push(px.pow(ex));
    }
    if v.len() == 1 {
        return v.pop().unwrap();
    }
    let combos = if use_iproduct { iproduct(&v) } else { product(&v) };
    combos.iter().map(|c| crt_coprime(c, &pv)).collect()
}

/// SymPy's `sqrt_mod(a, p)` (single root).
pub fn sqrt_mod(a: &BigInt, p: &BigInt) -> Option<BigInt> {
    let p = p.abs();
    let halfp: BigInt = &p / 2;
    let mut x = None;
    for r in sqrt_mod_iter(a, &p) {
        if r < halfp {
            return Some(r);
        } else if r > halfp {
            return Some(&p - r);
        }
        x = Some(r);
    }
    x
}

/// Is `a` a quadratic residue modulo `p`?
pub fn is_quad_residue(a: &BigInt, p: &BigInt) -> R<bool> {
    if p < &BigInt::one() {
        return Err("p must be > 0".into());
    }
    let mut p = p.clone();
    let mut a = a.mod_floor(&p);
    if a < b(2) || p < b(3) {
        return Ok(true);
    }
    let t = trailing(&p);
    if t > 0 {
        let a_ = &a % (BigInt::one() << t);
        if !a_.is_zero() {
            let r = trailing(&a_);
            if r % 2 == 1 || !((&a_ >> r) & b(6)).is_zero() {
                return Ok(false);
            }
        }
        p >>= t;
        a = a.mod_floor(&p);
        if a < b(2) || p < b(3) {
            return Ok(true);
        }
    }
    let j = jacobi_big(&a, &p);
    if j == -1 || isprime(&p) {
        return Ok(j == 1);
    }
    for (px, ex) in factorint(&p) {
        if !(&a % &px).is_zero() {
            if jacobi_big(&a, &px) != 1 {
                return Ok(false);
            }
        } else {
            let a_ = &a % px.pow(ex);
            if a_.is_zero() {
                continue;
            }
            let (a_, r) = remove(&a_, &px).unwrap();
            if r % 2 == 1 || jacobi_big(&a_, &px) != 1 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn is_nthpow_residue_bign_prime_power(a: &BigInt, n: &BigInt, p: &BigInt, k: u32) -> bool {
    let mut a = a.clone();
    let mut k = k as i64;
    while (&a % p).is_zero() {
        a = a.mod_floor(&p.pow(k.max(0) as u32));
        if a.is_zero() {
            return true;
        }
        let (a2, mu) = remove(&a, p).unwrap();
        a = a2;
        if !(BigInt::from(mu) % n).is_zero() {
            return false;
        }
        k -= mu as i64;
    }
    let k = k.max(1) as u32;
    if p != &b(2) {
        let f: BigInt = p.pow(k - 1) * (p - 1);
        return a.modpow(&(&f / f.gcd(n)), &p.pow(k)).is_one();
    }
    if n.is_odd() {
        return true;
    }
    let c = (trailing(n) + 2).min(k as u64);
    (a.mod_floor(&(BigInt::one() << c))).is_one()
}

/// Is `a` an `n`-th power residue modulo `m`?
pub fn is_nthpow_residue(a: &BigInt, n: &BigInt, m: &BigInt) -> R<bool> {
    // the caller has reduced `a` modulo `m` (m != 0)
    if !m.is_positive() {
        return Err("m must be > 0".into());
    }
    if n.is_negative() {
        return Err("n must be >= 0".into());
    }
    if n.is_zero() {
        if m.is_one() {
            return Ok(false);
        }
        return Ok(a.is_one());
    }
    if a.is_zero() || n.is_one() {
        return Ok(true);
    }
    if n == &b(2) {
        return is_quad_residue(a, m);
    }
    Ok(factorint(m).into_iter().all(|(p, e)| is_nthpow_residue_bign_prime_power(a, n, &p, e)))
}

fn nthroot_mod1(s: &BigInt, q: &BigInt, p: &BigInt) -> R<Vec<BigInt>> {
    let g = first_root_prime(p);
    let pm1: BigInt = p - 1;
    let mut r = s.clone();
    for (qx, ex) in factorint(q) {
        let mut f: BigInt = &pm1 / qx.pow(ex);
        while (&f % &qx).is_zero() {
            f /= &qx;
        }
        let z = &f * mod_inverse(&(-&f).mod_floor(&qx), &qx).unwrap();
        let x = (1 + &z) / &qx;
        let mut t = discrete_log(p, &r.modpow(&f, p), &g.modpow(&(&f * &qx), p), None, None)?;
        for _ in 0..ex {
            r = r.modpow(&x, p) * g.modpow(&(-&z * &t).mod_floor(&pm1), p) % p;
            t /= &qx;
        }
    }
    let mut res = vec![r.clone()];
    let h = g.modpow(&(&pm1 / q), p);
    let mut hx = r;
    let qn = q.to_u64().unwrap_or(0);
    for _ in 0..qn.saturating_sub(1) {
        hx = hx * &h % p;
        res.push(hx.clone());
    }
    res.sort();
    Ok(res)
}

fn nthroot_mod_prime_power(a: &BigInt, n: &BigInt, p: &BigInt, k: u32) -> R<Vec<BigInt>> {
    if !is_nthpow_residue_bign_prime_power(a, n, p, k) {
        return Ok(vec![]);
    }
    let mut a_mod_p = a.mod_floor(p);
    let base_roots = if a_mod_p.is_zero() {
        vec![BigInt::zero()]
    } else if ((p - 1u32) % n).is_zero() {
        nthroot_mod1(&a_mod_p, n, p)?
    } else {
        let mut pa = n.clone();
        let mut pb: BigInt = p - 1;
        let mut bb = BigInt::one();
        if pa < pb {
            std::mem::swap(&mut a_mod_p, &mut bb);
            std::mem::swap(&mut pa, &mut pb);
        }
        while !pb.is_zero() {
            let (q, pc) = pa.div_rem(&pb);
            let c = pow_mod(&bb, &-q, p).ok_or("base is not invertible")? * &a_mod_p % p;
            pa = pb;
            pb = pc;
            a_mod_p = bb;
            bb = c;
        }
        if pa.is_one() {
            vec![a_mod_p]
        } else if pa == b(2) {
            let mut v = sqrt_mod_iter(&a_mod_p, p);
            v.sort();
            v
        } else {
            nthroot_mod1(&a_mod_p, &pa, p)?
        }
    };
    if k == 1 {
        return Ok(base_roots);
    }
    let pk = p.pow(k);
    let a = a.mod_floor(&pk);
    let mut tot: std::collections::BTreeSet<BigInt> = Default::default();
    for root in base_roots {
        let diff = root.modpow(&(n - 1), p) * n % p;
        let mut new_base = p.clone();
        if !diff.is_zero() {
            let m_inv = mod_inverse(&diff, p).unwrap();
            let mut root = root;
            for _ in 0..k - 1 {
                new_base *= p;
                let tmp = (root.modpow(n, &new_base) - &a) * &m_inv;
                root = (root - tmp).mod_floor(&new_base);
            }
            tot.insert(root);
        } else {
            let mut roots_in_base: std::collections::BTreeSet<BigInt> = [root].into_iter().collect();
            for _ in 0..k - 1 {
                new_base *= p;
                let mut new_roots = std::collections::BTreeSet::new();
                for k_ in roots_in_base {
                    if k_.modpow(n, &new_base) != a.mod_floor(&new_base) {
                        continue;
                    }
                    let mut k_ = k_;
                    while !new_roots.contains(&k_) {
                        new_roots.insert(k_.clone());
                        k_ = (&k_ + &new_base / p).mod_floor(&new_base);
                    }
                }
                roots_in_base = new_roots;
            }
            tot.extend(roots_in_base);
        }
    }
    Ok(tot.into_iter().collect())
}

/// All `n`-th roots of `a` modulo `p` (`n >= 3`, `p >= 1`), sorted.
pub fn nthroot_mod_all(a: &BigInt, n: &BigInt, p: &BigInt) -> R<Vec<BigInt>> {
    let mut base = Vec::new();
    let mut pp = Vec::new();
    for (q, e) in factorint(p) {
        let mut roots = nthroot_mod_prime_power(a, n, &q, e)?;
        if roots.is_empty() {
            return Ok(vec![]);
        }
        roots.sort();
        pp.push(q.pow(e));
        base.push(roots);
    }
    let mut ret: Vec<BigInt> = product(&base).iter().map(|c| crt_coprime(c, &pp)).collect();
    ret.sort();
    ret.dedup();
    Ok(ret)
}

/// Sorted quadratic residues modulo `p`.
pub fn quadratic_residues(p: &BigInt) -> Vec<BigInt> {
    let mut r: std::collections::BTreeSet<BigInt> = Default::default();
    let half: BigInt = p.div_floor(&b(2));
    let mut i = BigInt::zero();
    while i <= half {
        r.insert(pow_mod(&i, &b(2), p).unwrap_or_default());
        i += 1;
    }
    r.into_iter().collect()
}

/// SymPy's `igcdex` (extended Euclid with its sign conventions).
pub fn igcdex(a: &BigInt, b_: &BigInt) -> (BigInt, BigInt, BigInt) {
    if a.is_zero() && b_.is_zero() {
        return (b(0), b(1), b(0));
    }
    if a.is_zero() {
        return (b(0), b_.signum(), b_.abs());
    }
    if b_.is_zero() {
        return (a.signum(), b(0), a.abs());
    }
    let (mut a, x_sign) = if a.is_negative() { (-a, -1) } else { (a.clone(), 1) };
    let (mut bb, y_sign) = if b_.is_negative() { (-b_, -1) } else { (b_.clone(), 1) };
    let (mut x, mut y, mut r, mut s) = (b(1), b(0), b(0), b(1));
    while !bb.is_zero() {
        let (q, c) = a.div_mod_floor(&bb);
        let nr = &x - &q * &r;
        let ns = &y - &q * &s;
        a = bb;
        bb = c;
        x = r;
        y = s;
        r = nr;
        s = ns;
    }
    (x * x_sign, y * y_sign, a)
}

/// SymPy's `linear_congruence(a, b, m)` (solutions in its order).
pub fn linear_congruence(a: &BigInt, bb: &BigInt, m: &BigInt) -> Vec<BigInt> {
    if (a % m).is_zero() {
        if (bb % m).is_zero() {
            return range_step(b(0), m, &b(1));
        }
        return vec![];
    }
    let (r, _, g) = igcdex(a, m);
    if !(bb.mod_floor(&g)).is_zero() {
        return vec![];
    }
    let base = (&r * bb).div_floor(&g);
    let step = m.div_floor(&g);
    let gn = g.to_u64().unwrap_or(0);
    (0..gn).map(|t| (&base + &step * t).mod_floor(m)).collect()
}

/// Solutions of `a*x**2 + b*x + c = 0 (mod n)`.  The flag is true when
/// SymPy returns Python ints (false: SymPy Integers).
pub fn quadratic_congruence(a: &BigInt, bb: &BigInt, c: &BigInt, n: &BigInt) -> R<(Vec<BigInt>, bool)> {
    if n <= &BigInt::one() {
        return Err("n should be an integer greater than 1".into());
    }
    let a = a.mod_floor(n);
    let mut bb = bb.mod_floor(n);
    let mut c = c.mod_floor(n);
    if a.is_zero() {
        let all = (bb.mod_floor(n)).is_zero() && (c.mod_floor(n)).is_zero();
        return Ok((linear_congruence(&bb, &(-&c), n), all));
    }
    if n == &b(2) {
        let mut roots = vec![];
        if c.is_zero() {
            roots.push(b(0));
        }
        if ((&bb + &c) % 2u32).is_one() {
            roots.push(b(1));
        }
        return Ok((roots, true));
    }
    if (&a * 2u32).gcd(n).is_one() {
        let inv_a = mod_inverse(&a, n).unwrap();
        bb *= &inv_a;
        c *= &inv_a;
        if bb.is_odd() {
            bb += n;
        }
        bb >>= 1;
        let mut v: Vec<BigInt> = sqrt_mod_iter(&(&bb * &bb - &c), n).into_iter().map(|i| (i - &bb).mod_floor(n)).collect();
        v.sort();
        return Ok((v, true));
    }
    let mut res = std::collections::BTreeSet::new();
    let two_a: BigInt = &a * 2;
    for i in sqrt_mod_iter(&(&bb * &bb - 4 * &a * &c), &(4 * &a * n)) {
        let (q, rem) = (i - &bb).div_mod_floor(&two_a);
        if rem.is_zero() {
            res.insert(q.mod_floor(n));
        }
    }
    Ok((res.into_iter().collect(), true))
}

fn discrete_log_trial_mul(n: &BigInt, a: &BigInt, bb: &BigInt, order: &BigInt) -> R<BigInt> {
    let a = a.mod_floor(n);
    let bb = bb.mod_floor(n);
    let mut x = BigInt::one();
    let mut i = BigInt::zero();
    while &i < order {
        if x == a {
            return Ok(i);
        }
        x = x * &bb % n;
        i += 1;
    }
    Err("Log does not exist".into())
}

fn discrete_log_shanks(n: &BigInt, a: &BigInt, bb: &BigInt, order: &BigInt) -> R<BigInt> {
    let a = a.mod_floor(n);
    let bb = bb.mod_floor(n);
    let m: BigInt = order.sqrt() + 1;
    let mu = m.to_u64().ok_or("order too large")?;
    let mut t: std::collections::HashMap<BigInt, u64> = Default::default();
    let mut x = BigInt::one();
    for i in 0..mu {
        t.insert(x.clone(), i);
        x = x * &bb % n;
    }
    let z = pow_mod(&bb, &-&m, n).ok_or("Log does not exist")?;
    let mut x = a;
    for i in 0..mu {
        if let Some(j) = t.get(&x) {
            return Ok(BigInt::from(i) * &m + j);
        }
        x = x * &z % n;
    }
    Err("Log does not exist".into())
}

fn discrete_log_pollard_rho(n: &BigInt, a: &BigInt, bb: &BigInt, order: &BigInt) -> R<BigInt> {
    let a = a.mod_floor(n);
    let bb = bb.mod_floor(n);
    let mut rng = PyRandom::new(&b(1234));
    let three = b(3);
    let step = |x: &mut BigInt, ea: &mut BigInt, eb: &mut BigInt| {
        let c = (&*x % &three).to_u32().unwrap();
        if c == 0 {
            *x = &a * &*x % n;
            *eb = (&*eb + 1u32).mod_floor(order);
        } else if c == 1 {
            *x = &*x * &*x % n;
            *ea = (&*ea * 2u32).mod_floor(order);
            *eb = (&*eb * 2u32).mod_floor(order);
        } else {
            *x = &bb * &*x % n;
            *ea = (&*ea + 1u32).mod_floor(order);
        }
    };
    for _ in 0..10 {
        let mut aa = rng.randint(&b(1), &(order - 1));
        let mut ba = rng.randint(&b(1), &(order - 1));
        let mut xa = bb.modpow(&aa, n) * a.modpow(&ba, n) % n;
        let (mut xb, mut ab, mut b2) = (xa.clone(), aa.clone(), ba.clone());
        step(&mut xb, &mut ab, &mut b2);
        let mut j = BigInt::zero();
        while &j < order {
            j += 1;
            step(&mut xa, &mut aa, &mut ba);
            step(&mut xb, &mut ab, &mut b2);
            step(&mut xb, &mut ab, &mut b2);
            if xa == xb {
                let r = (&ba - &b2).mod_floor(order);
                if let Some(inv) = mod_inverse(&r, order) {
                    let e = (inv * (&ab - &aa)).mod_floor(order);
                    if (bb.modpow(&e, n) - &a).mod_floor(n).is_zero() {
                        return Ok(e);
                    }
                }
                break;
            }
        }
    }
    Err("Pollard's Rho failed to find logarithm".into())
}

fn discrete_log_pohlig_hellman(n: &BigInt, a: &BigInt, bb: &BigInt, order: &BigInt, order_factors: &[(BigInt, u32)]) -> R<BigInt> {
    let a = a.mod_floor(n);
    let bb = bb.mod_floor(n);
    let mut l = vec![BigInt::zero(); order_factors.len()];
    for (i, (pi, ri)) in order_factors.iter().enumerate() {
        for j in 0..*ri {
            let inner = pow_mod(&bb, &-&l[i], n).ok_or("Log does not exist")?;
            let aj = (&a * inner).modpow(&(order / pi.pow(j + 1)), n);
            let bj = bb.modpow(&(order / pi), n);
            let cj = discrete_log(n, &aj, &bj, Some(pi.clone()), Some(true))?;
            l[i] += cj * pi.pow(j);
        }
    }
    let ms: Vec<BigInt> = order_factors.iter().map(|(p, r)| p.pow(*r)).collect();
    Ok(super::crt(&l, &ms).map(|r| r.0).unwrap_or_default())
}

/// Discrete logarithm of `a` to the base `b` modulo `n`.
pub fn discrete_log(n: &BigInt, a: &BigInt, bb: &BigInt, order: Option<BigInt>, prime_order: Option<bool>) -> R<BigInt> {
    if !n.is_positive() {
        return Err("n should be positive".into());
    }
    if n.is_one() {
        return Ok(BigInt::zero());
    }
    let (order, order_factors) = match order {
        Some(o) => {
            let f = factorint(&o).into_iter().collect::<Vec<_>>();
            (o, f)
        }
        None => {
            let mut factors: Vec<(BigInt, u32)> = Vec::new();
            let add = |p: BigInt, k: u32, fs: &mut Vec<(BigInt, u32)>| match fs.iter_mut().find(|(q, _)| *q == p) {
                Some(slot) => slot.1 += k,
                None => fs.push((p, k)),
            };
            for (px, kx) in factorint(n) {
                if kx > 1 {
                    add(px.clone(), kx - 1, &mut factors);
                }
                for (py, ky) in factorint(&(&px - 1)) {
                    add(py, ky, &mut factors);
                }
            }
            let mut order = BigInt::one();
            for (px, kx) in &factors {
                order *= px.pow(*kx);
            }
            let mut of = Vec::new();
            for (p, e) in factors {
                let mut i = 0;
                for _ in 0..e {
                    if bb.modpow(&(&order / &p), n).is_one() {
                        order /= &p;
                        i += 1;
                    } else {
                        break;
                    }
                }
                if i < e {
                    of.push((p, e - i));
                }
            }
            (order, of)
        }
    };
    let prime_order = prime_order.unwrap_or_else(|| isprime(&order));
    if order < b(1000) {
        return discrete_log_trial_mul(n, a, bb, &order);
    }
    if prime_order {
        if order < b(1_000_000_000_000) {
            return discrete_log_shanks(n, a, bb, &order);
        }
        return discrete_log_pollard_rho(n, a, bb, &order);
    }
    discrete_log_pohlig_hellman(n, a, bb, &order, &order_factors)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(xs: &[i64]) -> Vec<BigInt> {
        xs.iter().map(|x| b(*x)).collect()
    }

    #[test]
    fn residues() {
        assert_eq!(sqrt_mod_iter(&b(4), &b(15)), v(&[7, 13, 2, 8]));
        let mut r = sqrt_mod_iter(&b(17), &b(32));
        r.sort();
        assert_eq!(r, v(&[7, 9, 23, 25]));
        assert_eq!(nthroot_mod_all(&b(11), &b(3), &b(19)).unwrap(), v(&[5, 16, 17]));
        assert_eq!(quadratic_congruence(&b(2), &b(5), &b(3), &b(7)).unwrap().0, v(&[2, 6]));
        assert_eq!(discrete_log(&b(41), &b(15), &b(7), None, None).unwrap(), b(3));
        assert_eq!(primitive_root(&b(18), true).unwrap(), Some(b(5)));
        assert_eq!(primitive_root(&b(50), true).unwrap(), Some(b(3)));
        assert_eq!(primitive_root(&b(97), false).unwrap(), Some(b(5)));
        assert_eq!(n_order(&b(3), &b(7)).unwrap(), b(6));
        assert_eq!(kronecker(&b(3), &b(8)), -1);
        assert_eq!(linear_congruence(&b(3), &b(6), &b(15)), v(&[2, 7, 12]));
    }
}
