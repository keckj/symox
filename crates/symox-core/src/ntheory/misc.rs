//! Integer helpers ported from SymPy: `integer_log`, `num_digits`,
//! Egyptian fractions, periodic continued fractions and binomial rows.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use super::factor::trailing;
use super::mod_inverse;

type R<T> = Result<T, String>;

/// SymPy's `integer_log(n, b)`: `(e, exact)` with `e = floor(log_b |n|)`.
pub fn integer_log(n: &BigInt, b: &BigInt) -> R<(u64, bool)> {
    if b.is_negative() {
        let (e, t) = integer_log(&n.abs(), &-b)?;
        return Ok((e, t && ((e % 2 == 1) == n.is_negative())));
    }
    if b <= &BigInt::one() {
        return Err("base must be 2 or more".into());
    }
    if n.is_negative() {
        if b != &BigInt::from(2) {
            return Ok((0, false));
        }
        let (e, _) = integer_log(&-n, b)?;
        return Ok((e, false));
    }
    if n.is_zero() {
        return Err("n cannot be 0".into());
    }
    if n < b {
        return Ok((0, n.is_one()));
    }
    if b == &BigInt::from(2) {
        let e = n.bits() - 1;
        return Ok((e, trailing(n) == e));
    }
    let t = trailing(b);
    if (BigInt::one() << t) == *b {
        let e = (n.bits() - 1) / t;
        return Ok((e, (BigInt::one() << (t * e)) == *n));
    }
    // floor estimate from the bit lengths, then exact correction
    let mut e = ((n.bits() - 1) as f64 / b.bits() as f64).floor() as u64;
    let mut pw = b.pow(e as u32);
    while &pw > n {
        e -= 1;
        pw /= b;
    }
    loop {
        let nx = &pw * b;
        if &nx > n {
            break;
        }
        pw = nx;
        e += 1;
    }
    Ok((e, &pw == n))
}

/// SymPy's `egypt_harmonic(r)` on `x/y`: the longest harmonic prefix
/// `1, 2, ..., d` with `1 + 1/2 + ... + 1/d <= x/y` and the reduced
/// remainder `(num, den)`.
pub fn egypt_harmonic(x: &BigInt, y: &BigInt) -> (Vec<BigInt>, BigInt, BigInt) {
    let mut rv = Vec::new();
    let (mut an, mut ad) = (BigInt::zero(), BigInt::one());
    let mut d = BigInt::one();
    loop {
        let nn = &an * &d + &ad;
        let nd = &ad * &d;
        if &nn * y <= x * &nd {
            let g = nn.gcd(&nd);
            an = nn / &g;
            ad = nd / &g;
            rv.push(d.clone());
            d += 1;
        } else {
            break;
        }
    }
    let rn = x * &ad - &an * y;
    let rd = y * &ad;
    let g = rn.gcd(&rd);
    if g.is_zero() {
        return (rv, rn, rd);
    }
    (rv, rn / &g, rd / &g)
}

/// SymPy's `egyptian_fraction` on a positive rational `x/y` (reduced).
pub fn egyptian_fraction(x: &BigInt, y: &BigInt, algorithm: &str) -> R<Vec<BigInt>> {
    let one = BigInt::one();
    if y.is_one() && x == &BigInt::from(2) {
        return Ok([1, 2, 3, 6].iter().map(|&i| BigInt::from(i)).collect());
    }
    if *x == y + 1u32 {
        return Ok(vec![one, y.clone()]);
    }
    let (mut rv, rx, ry) = egypt_harmonic(x, y);
    if rx.is_zero() {
        return Ok(rv);
    }
    let post = match algorithm {
        "Greedy" => egypt_greedy(&rx, &ry),
        "Graham Jewett" => egypt_graham_jewett(&rx, &ry),
        "Takenouchi" => egypt_takenouchi(&rx, &ry),
        "Golomb" => egypt_golomb(&rx, &ry),
        _ => return Err("Entered invalid algorithm".into()),
    };
    rv.extend(post);
    Ok(rv)
}

pub fn egypt_greedy(x: &BigInt, y: &BigInt) -> Vec<BigInt> {
    let (mut x, mut y) = (x.clone(), y.clone());
    let mut out = Vec::new();
    loop {
        if x.is_one() {
            out.push(y);
            return out;
        }
        let a = (-&y).mod_floor(&x);
        let q = y.div_floor(&x) + 1u32;
        let bb = &y * &q;
        let c = a.gcd(&bb);
        out.push(q);
        if c > BigInt::one() {
            x = a / &c;
            y = bb / &c;
        } else {
            x = a;
            y = bb;
        }
    }
}

fn first_dup(l: &[BigInt]) -> Option<usize> {
    (0..l.len().saturating_sub(1)).find(|&i| l[i] == l[i + 1])
}

pub fn egypt_graham_jewett(x: &BigInt, y: &BigInt) -> Vec<BigInt> {
    let n: usize = x.try_into().unwrap_or(0);
    let mut l = vec![y.clone(); n];
    loop {
        l.sort();
        let Some(i) = first_dup(&l) else { break };
        l[i + 1] = &l[i] + 1u32;
        let v = &l[i] * (&l[i] + 1u32);
        l.push(v);
    }
    l.sort();
    l
}

pub fn egypt_takenouchi(x: &BigInt, y: &BigInt) -> Vec<BigInt> {
    if x == &BigInt::from(3) {
        if y.is_even() {
            return vec![y / 2, y.clone()];
        }
        let i: BigInt = (y - 1u32) / 2;
        let j = &i + 1u32;
        let k = &j + &i;
        let jk = &j * &k;
        return vec![j, k, jk];
    }
    let n: usize = x.try_into().unwrap_or(0);
    let mut l = vec![y.clone(); n];
    loop {
        l.sort();
        let Some(i) = first_dup(&l) else { break };
        let k = l[i].clone();
        if k.is_even() {
            l[i] = &k / 2;
            l.remove(i + 1);
        } else {
            l[i] = (&k + 1u32) / 2;
            l[i + 1] = &k * (&k + 1u32) / 2;
        }
    }
    l.sort();
    l
}

pub fn egypt_golomb(x: &BigInt, y: &BigInt) -> Vec<BigInt> {
    if x.is_one() {
        return vec![y.clone()];
    }
    let xp = mod_inverse(x, y).unwrap_or_default();
    let mut rv = vec![&xp * y];
    rv.extend(egypt_golomb(&((x * &xp - 1u32) / y), &xp));
    rv.sort();
    rv
}

/// Continued fraction of a rational `p/q` (floor based, `q != 0`).
pub fn cf_rational(p: &BigInt, q: &BigInt) -> Vec<BigInt> {
    let (mut p, mut q) = (p.clone(), q.clone());
    let mut out = Vec::new();
    while !q.is_zero() {
        let (a, r) = p.div_mod_floor(&q);
        out.push(a);
        p = q;
        q = r;
    }
    out
}

/// floor((p + sigma*sqrt(dd)) / q) for non-square `dd`, `r = isqrt(dd)`.
fn floor_quad(p: &BigInt, q: &BigInt, r: &BigInt, neg: bool) -> BigInt {
    let num = if neg { p - r - 1u32 } else { p + r };
    if q.is_positive() {
        num.div_floor(q)
    } else {
        -(num.div_floor(&-q) + 1u32)
    }
}

/// SymPy's `continued_fraction_periodic(p, q, d, s)`: the terms and the
/// start index of the period (`None` for a rational number).
pub fn continued_fraction_periodic(p: &BigInt, q: &BigInt, d: &BigInt, s: &BigInt) -> R<(Vec<BigInt>, Option<usize>)> {
    if d.is_negative() {
        return Err(format!("expected non-negative for `d` but got {d}"));
    }
    if q.is_zero() {
        return Err("The denominator cannot be 0.".into());
    }
    let d = if s.is_zero() { BigInt::zero() } else { d.clone() };
    let sd = d.sqrt();
    if &sd * &sd == d {
        return Ok((cf_rational(&(p + s * sd), q), None));
    }
    let (mut p, mut q, mut s) = (p.clone(), q.clone(), s.clone());
    if q.is_negative() {
        p = -p;
        q = -q;
        s = -s;
    }
    let neg = s.is_negative();
    let mut dd = d * &s * &s;
    if !(&dd - &p * &p).mod_floor(&q).is_zero() {
        dd *= &q * &q;
        p *= &q;
        q = &q * &q;
    }
    let r = dd.sqrt();
    let mut terms = Vec::new();
    let mut seen: std::collections::HashMap<(BigInt, BigInt), usize> = Default::default();
    loop {
        if let Some(&i) = seen.get(&(p.clone(), q.clone())) {
            return Ok((terms, Some(i)));
        }
        seen.insert((p.clone(), q.clone()), terms.len());
        let a = floor_quad(&p, &q, &r, neg);
        p = &a * &q - &p;
        q = (&dd - &p * &p) / &q;
        terms.push(a);
    }
}

/// Row `n` of Pascal's triangle.
pub fn binomial_row(n: u64) -> Vec<BigInt> {
    let mut row = vec![BigInt::one(); n as usize + 1];
    for k in 1..=(n as usize) / 2 {
        let v = &row[k - 1] * BigInt::from(n - k as u64 + 1) / BigInt::from(k as u64);
        row[k] = v.clone();
        row[n as usize - k] = v;
    }
    row
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(x: i64) -> BigInt {
        BigInt::from(x)
    }

    #[test]
    fn misc() {
        assert_eq!(integer_log(&b(125), &b(5)).unwrap(), (3, true));
        assert_eq!(integer_log(&b(17), &b(9)).unwrap(), (1, false));
        assert_eq!(integer_log(&b(-125), &b(-5)).unwrap(), (3, true));
        assert_eq!(integer_log(&b(4), &b(-2)).unwrap(), (2, true));
        assert_eq!(egyptian_fraction(&b(3), &b(7), "Greedy").unwrap(), vec![b(3), b(11), b(231)]);
        assert_eq!(egyptian_fraction(&b(3), &b(7), "Graham Jewett").unwrap(), vec![b(7), b(8), b(9), b(56), b(57), b(72), b(3192)]);
        assert_eq!(egyptian_fraction(&b(3), &b(7), "Takenouchi").unwrap(), vec![b(4), b(7), b(28)]);
        assert_eq!(egyptian_fraction(&b(3), &b(7), "Golomb").unwrap(), vec![b(3), b(15), b(35)]);
        assert_eq!(egyptian_fraction(&b(11), &b(5), "Greedy").unwrap(), vec![b(1), b(2), b(3), b(4), b(9), b(180)]);
        let (t, i) = continued_fraction_periodic(&b(3), &b(2), &b(7), &b(1)).unwrap();
        assert_eq!((t, i), (vec![b(2), b(1), b(4), b(1), b(1)], Some(1)));
        let (t, i) = continued_fraction_periodic(&b(0), &b(1), &b(2), &b(-1)).unwrap();
        assert_eq!((t, i), (vec![b(-2), b(1), b(1), b(2)], Some(3)));
        assert_eq!(binomial_row(4), vec![b(1), b(4), b(6), b(4), b(1)]);
    }
}
