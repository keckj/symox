//! Modular resultants of univariate `ZZ`/`QQ` polynomials (Collins).
//!
//! The Sylvester determinant is computed modulo word primes that do not
//! divide either leading coefficient (the reduced Sylvester matrix then has
//! the same shape, so every image is exact), and the images are combined by
//! the Chinese remainder theorem until the modulus exceeds twice the
//! Hadamard bound `||a||_2^deg(b) ||b||_2^deg(a)`. No image is guessed or
//! stabilised: the result is the exact resultant, equal to the one of the
//! subresultant PRS, with intermediates of word size instead of the
//! coefficient growth of the PRS.
use super::convert::{PResult, PolyError};
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{ToPrimitive, Zero};

/// Largest resultant size (bits of the Hadamard bound) handled here: the
/// CRT accumulator and the result then stay below 32 MiB.
const MAX_BOUND_BITS: u64 = 1 << 28;

/// Largest degree converted to dense images (`rootisolation::MAX_DENSE_DEGREE`).
const MAX_DENSE_DEGREE: i64 = 1 << 20;

fn limit(reason: &str) -> PolyError {
    PolyError::Unsupported(format!("pseudo-division resource limit: {reason}"))
}

/// Barrett reduction for a prime `p < 2^31`; inputs below `2^64`.
#[derive(Clone, Copy)]
struct Modulus {
    p: u64,
    m: u64,
}
impl Modulus {
    fn new(p: u64) -> Self {
        Self { p, m: u64::MAX / p }
    }
    #[inline]
    fn reduce(self, x: u64) -> u64 {
        let q = ((x as u128 * self.m as u128) >> 64) as u64;
        let mut r = x - q * self.p;
        while r >= self.p {
            r -= self.p;
        }
        r
    }
    #[inline]
    fn mul(self, a: u64, b: u64) -> u64 {
        self.reduce(a * b)
    }
    fn pow(self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1 % self.p;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    fn inv(self, a: u64) -> u64 {
        self.pow(a, self.p - 2)
    }
}

/// Integer coefficient split into its sign and base-2^32 digits.
struct Digits {
    negative: bool,
    digits: Vec<u32>,
}
impl Digits {
    fn new(c: &Int) -> Self {
        let (sign, digits) = c.to_big().to_u32_digits();
        Self { negative: sign == Sign::Minus, digits }
    }
    fn residue(&self, m: Modulus) -> u64 {
        let mut r = 0u64;
        for &d in self.digits.iter().rev() {
            r = m.reduce((r << 32) | d as u64);
        }
        if self.negative && r != 0 {
            m.p - r
        } else {
            r
        }
    }
}

/// `res(a, b)` modulo `p` by the Euclidean algorithm, for nonzero leading
/// coefficients and degrees at least one.
fn resultant_mod(a: &[Digits], b: &[Digits], m: Modulus) -> u64 {
    let p = m.p;
    let mut f: Vec<u64> = a.iter().map(|c| c.residue(m)).collect();
    let mut g: Vec<u64> = b.iter().map(|c| c.residue(m)).collect();
    let mut acc = 1u64;
    loop {
        let (df, dg) = (f.len() - 1, g.len() - 1);
        let lead = g[dg];
        if dg == 0 {
            return m.mul(acc, m.pow(lead, df as u64));
        }
        // f <- f mod g, in place.
        let inverse = m.inv(lead);
        if df >= dg {
            for i in (dg..=df).rev() {
                let q = m.mul(f[i], inverse);
                if q != 0 {
                    let negated = p - q;
                    let base = i - dg;
                    for (j, &c) in g[..dg].iter().enumerate() {
                        f[base + j] = m.reduce(f[base + j] + negated * c);
                    }
                }
            }
            f.truncate(dg);
        }
        while f.last() == Some(&0) {
            f.pop();
        }
        if f.is_empty() {
            return 0;
        }
        let dr = f.len() - 1;
        // res(f,g) = (-1)^(df dg) res(g,f) = (-1)^(df dg) lc(g)^(df-dr) res(g, f mod g)
        if df % 2 == 1 && dg % 2 == 1 && acc != 0 {
            acc = p - acc;
        }
        acc = m.mul(acc, m.pow(lead, (df - dr) as u64));
        std::mem::swap(&mut f, &mut g);
    }
}

fn norm_square_bits(c: &[Int]) -> u64 {
    let sum = c.iter().fold(BigInt::zero(), |s, c| {
        let b = c.to_big();
        s + &b * &b
    });
    sum.bits()
}

/// Exact resultant of dense ascending integer polynomials of degree at
/// least one with nonzero leading coefficients.
pub fn resultant_zz(a: &[Int], b: &[Int]) -> PResult<Int> {
    if a.len() < 2 || b.len() < 2 || a[a.len() - 1].is_zero() || b[b.len() - 1].is_zero() {
        return Err(PolyError::Failed(
            "modular resultant needs nonconstant polynomials".into(),
        ));
    }
    let (n, m) = ((a.len() - 1) as u64, (b.len() - 1) as u64);
    let bound = m
        .checked_mul(norm_square_bits(a))
        .and_then(|x| n.checked_mul(norm_square_bits(b)).and_then(|y| x.checked_add(y)))
        .map(|x| x.div_ceil(2))
        .ok_or_else(|| limit("resultant size"))?;
    if bound > MAX_BOUND_BITS {
        return Err(limit("resultant size"));
    }
    // |res| <= 2^bound, so a modulus of more than bound + 1 bits determines
    // its symmetric representative.
    let target = bound + 2;
    let a: Vec<Digits> = a.iter().map(Digits::new).collect();
    let b: Vec<Digits> = b.iter().map(Digits::new).collect();
    let mut value = BigUint::zero();
    let mut modulus = BigUint::from(1u32);
    let mut prime = 1u64 << 31;
    while modulus.bits() <= target {
        prime = previous_prime(prime).ok_or_else(|| limit("resultant primes"))?;
        let field = Modulus::new(prime);
        if a[a.len() - 1].residue(field) == 0 || b[b.len() - 1].residue(field) == 0 {
            continue;
        }
        let image = resultant_mod(&a, &b, field);
        // Garner step: value + modulus * t is congruent to image mod prime.
        let current = (&value % prime).to_u64().unwrap();
        let scale = (&modulus % prime).to_u64().unwrap();
        let t = field.mul((image + prime - current) % prime, field.inv(scale));
        value += &modulus * t;
        modulus *= prime;
    }
    let half = &modulus >> 1u32;
    let result = if value > half {
        BigInt::from(value) - BigInt::from(modulus)
    } else {
        BigInt::from(value)
    };
    Ok(Int::from_big(result))
}

fn previous_prime(mut n: u64) -> Option<u64> {
    while n > 3 {
        n -= 1;
        if super::gf::is_prime(n) {
            return Some(n);
        }
    }
    None
}

fn integral(c: &[Rat]) -> (Int, Vec<Int>) {
    let denominator = c.iter().fold(Int::one(), |l, c| l.lcm(&c.d));
    let values = c
        .iter()
        .map(|c| c.n.mul(&denominator.div_exact(&c.d).expect("lcm multiple")))
        .collect();
    (denominator, values)
}

/// Exact resultant of dense ascending rational polynomials of degree at
/// least one: `res(a/da, b/db) = res(a, b) / (da^deg(b) db^deg(a))`.
pub fn resultant_qq(a: &[Rat], b: &[Rat]) -> PResult<Rat> {
    let (da, za) = integral(a);
    let (db, zb) = integral(b);
    let value = resultant_zz(&za, &zb)?;
    let (n, m) = (a.len() as u32 - 1, b.len() as u32 - 1);
    Ok(Rat::new(value, da.pow(m).mul(&db.pow(n))))
}

/// The resultant in the only generator of two univariate rational
/// polynomials of positive degree (`None` for every other shape, which the
/// subresultant PRS handles).
pub fn resultant_univariate(a: &MPoly<Rat>, b: &MPoly<Rat>) -> PResult<Option<MPoly<Rat>>> {
    if a.nvars != 1 || b.nvars != 1 || a.degree(0) < 1 || b.degree(0) < 1 {
        return Ok(None);
    }
    // Dense images only for a bounded dense storage; sparse giant degrees
    // keep the PRS path and its resource limits.
    if a.degree(0) > MAX_DENSE_DEGREE || b.degree(0) > MAX_DENSE_DEGREE {
        return Ok(None);
    }
    let value = resultant_qq(
        &super::dense::from_mpoly(a, 0),
        &super::dense::from_mpoly(b, 0),
    )?;
    Ok(Some(if value.is_zero() {
        MPoly::zero(1)
    } else {
        MPoly::from_terms(1, vec![(vec![0], value)])
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ints(c: &[i64]) -> Vec<Int> {
        c.iter().map(|&c| Int::S(c)).collect()
    }
    fn prs(a: &[Int], b: &[Int]) -> Int {
        let to = |c: &[Int]| {
            MPoly::from_terms(
                1,
                c.iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(i, c)| (vec![i as u32], Rat::new(c.clone(), Int::one())))
                    .collect(),
            )
        };
        let r = super::super::prs::subresultants(&to(a), &to(b), 0)
            .unwrap()
            .resultant;
        r.coeffs().first().map(|c| c.n.clone()).unwrap_or_else(Int::zero)
    }

    #[test]
    fn agrees_with_the_subresultant_prs() {
        let mut state = 12345u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((state >> 33) % 41) as i64 - 20
        };
        for _ in 0..200 {
            let da = 1 + (next().unsigned_abs() % 7) as usize;
            let db = 1 + (next().unsigned_abs() % 7) as usize;
            let mut a: Vec<i64> = (0..=da).map(|_| next()).collect();
            let mut b: Vec<i64> = (0..=db).map(|_| next()).collect();
            if a[da] == 0 {
                a[da] = 3;
            }
            if b[db] == 0 {
                b[db] = -2;
            }
            let (a, b) = (ints(&a), ints(&b));
            assert_eq!(resultant_zz(&a, &b).unwrap(), prs(&a, &b));
            assert_eq!(resultant_zz(&b, &a).unwrap(), prs(&b, &a));
        }
    }

    #[test]
    fn common_roots_and_closed_forms() {
        // (x-1)(x-2) and (x-1)(x+5) share a root.
        assert_eq!(resultant_zz(&ints(&[2, -3, 1]), &ints(&[-5, 4, 1])).unwrap(), Int::zero());
        // res((x+1)^3, x - 2) = (2+1)^3 (sign (-1)^(3*1) for the swap is absorbed).
        assert_eq!(resultant_zz(&ints(&[1, 3, 3, 1]), &ints(&[-2, 1])).unwrap(), Int::S(-27));
        assert_eq!(resultant_zz(&ints(&[-2, 1]), &ints(&[1, 3, 3, 1])).unwrap(), Int::S(27));
        let half = Rat::new(Int::one(), Int::S(2));
        // res(x/2 - 1, x^2 + 1) = (1/2)^2 * (2^2 + 1)
        assert_eq!(
            resultant_qq(&[Rat::new(Int::S(-1), Int::one()), half.clone()], &[Rat::one(), Rat::zero(), Rat::one()]).unwrap(),
            Rat::new(Int::S(5), Int::S(4))
        );
    }

    #[test]
    fn large_powers_have_the_closed_form() {
        // res((x+1)^40, (x-2)^30 + 1) = ((-3)^30 + 1)^40 since (x+1)^40 is monic.
        let binomial = |n: u32, shift: i64| -> Vec<Int> {
            let mut c = vec![Int::one()];
            for _ in 0..n {
                let mut next = vec![Int::zero(); c.len() + 1];
                for (i, v) in c.iter().enumerate() {
                    next[i + 1] = next[i + 1].add(v);
                    next[i] = next[i].add(&v.mul(&Int::S(shift)));
                }
                c = next;
            }
            c
        };
        let a = binomial(40, 1);
        let mut b = binomial(30, -2);
        b[0] = b[0].add(&Int::one());
        let expected = Int::S(3).pow(30).add(&Int::one()).pow(40);
        assert_eq!(resultant_zz(&a, &b).unwrap(), expected);
    }
}
