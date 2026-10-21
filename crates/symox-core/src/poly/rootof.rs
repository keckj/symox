//! Indexed roots of irreducible polynomials over ZZ (the reference's
//! `ComplexRootOf`).
//!
//! A root object is the head `ComplexRootOf` applied to the defining
//! polynomial, as an expression in its generator, and the integer index of
//! the root. Real roots come first in increasing order, then the non-real
//! roots in the reference's isolation order, each lower half plane root
//! followed by its conjugate. The generator of the defining polynomial is
//! bound: the node has no free symbols, so substitution and differentiation
//! treat it as a constant, as in the reference.
//!
//! The assumptions attached to the node are those the reference derives from
//! its isolating intervals: reality from the real root count, the sign of a
//! real root from its interval, and pure imaginarity from refinement of the
//! complex rectangles.
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

use super::int::{Int, Rat};
use super::polyobj::Poly;
use super::rootisolation as iso;
use crate::assumption_facts::FactSet;
use crate::expr::{Expr, Func, Kind, UndefinedFunction};
use crate::functions::func;

/// Name of the head of root objects.
pub const HEAD: &str = "ComplexRootOf";

/// Whether `head` with `nargs` arguments is the root object head.
pub fn is_rootof_head(head: &UndefinedFunction, nargs: usize) -> bool {
    nargs == 2
        && !head.internal
        && head.host().is_none()
        && head.name.as_ref() == HEAD
        && head.facts().get("algebraic") == Some(true)
}

/// Whether `e` is a root object.
pub fn is_rootof(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Undefined(head), args) if is_rootof_head(head, args.len()))
}

/// The defining expression and the index of a root object.
pub fn parts(e: &Expr) -> Option<(Expr, usize)> {
    let Kind::Func(Func::Undefined(head), args) = e.kind() else {
        return None;
    };
    if !is_rootof_head(head, args.len()) {
        return None;
    }
    let index = Int::from_number(args[1].as_number()?)?.as_i64()?;
    Some((args[0].clone(), usize::try_from(index).ok()?))
}

/// The generator and the integer coefficients, leading first, of a
/// univariate polynomial expression over ZZ.
pub fn integer_poly(expr: &Expr) -> Option<(Expr, Vec<Int>)> {
    let symbols = crate::traverse::free_symbols(expr);
    let [gen] = symbols.as_slice() else {
        return None;
    };
    let p = Poly::from_expr(expr, std::slice::from_ref(gen)).ok()?;
    let coeffs = p
        .all_coeffs()
        .ok()?
        .iter()
        .map(|c| Int::from_number(c.as_number()?))
        .collect::<Option<Vec<_>>>()?;
    Some((gen.clone(), coeffs))
}

/// The location of a root of a square-free polynomial over ZZ.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RootKind {
    pub real: bool,
    pub negative: bool,
    pub imaginary: bool,
    /// The root lies in the lower half plane.
    pub lower: bool,
    /// Number of real roots of the polynomial.
    pub nreal: usize,
}

fn rats(f: &[Int]) -> Vec<Rat> {
    f.iter().map(|c| Rat::new(c.clone(), Int::one())).collect()
}

/// `f(i*y)` for a polynomial `f` in even powers only, leading first.
fn imaginary_axis(f: &[Int]) -> Vec<Rat> {
    let n = f.len() - 1;
    f.iter()
        .enumerate()
        .map(|(k, c)| {
            let degree = n - k;
            let c = Rat::new(c.clone(), Int::one());
            if degree % 4 == 2 {
                c.neg()
            } else {
                c
            }
        })
        .collect()
}

/// Indices, among the non-real roots, of the purely imaginary roots of the
/// irreducible `f`, by the reference's refinement of the rectangles.
fn imaginary_indices(f: &[Int]) -> Option<Vec<bool>> {
    let n = f.len() - 1;
    let even = f.iter().enumerate().all(|(k, c)| c.is_zero() || (n - k) % 2 == 0);
    let mut boxes = iso::isolate_complex_roots_sqf(&rats(f), None, None, None).ok()?;
    if !even {
        return Some(vec![false; boxes.len()]);
    }
    let nimag = iso::count_real_roots(&imaginary_axis(f), None, None);
    let mut potential: Vec<usize> = (0..boxes.len()).collect();
    while potential.len() != nimag {
        if potential.len() < nimag {
            return None;
        }
        let mut kept = Vec::with_capacity(potential.len());
        for &i in &potential {
            let b = &boxes[i];
            let (ax, bx) = (b.ax(), b.bx());
            if ax.mul(&bx) > Rat::zero() {
                continue;
            }
            if ax != bx {
                boxes[i] = boxes[i].refine().ok()?;
            }
            kept.push(i);
        }
        potential = kept;
    }
    let mut flags = vec![false; boxes.len()];
    for i in potential {
        flags[i] = true;
    }
    Some(flags)
}

/// Where root `index` of the irreducible `f` over ZZ lies.
pub fn root_kind(f: &[Int], index: usize) -> Option<RootKind> {
    let degree = f.len().checked_sub(1)?;
    if index >= degree {
        return None;
    }
    let (reals, _) = iso::isolate_real_roots_sqf(f, None, None, None, false).ok()?;
    let nreal = reals.len();
    if index < nreal {
        let negative = reals[index].0.is_negative();
        return Some(RootKind { real: true, negative, imaginary: false, lower: false, nreal });
    }
    let k = index - nreal;
    let imaginary = imaginary_flag(f, k)?;
    Some(RootKind { real: false, negative: false, imaginary, lower: k % 2 == 0, nreal })
}

fn imaginary_flag(f: &[Int], k: usize) -> Option<bool> {
    let n = f.len() - 1;
    let even = f.iter().enumerate().all(|(j, c)| c.is_zero() || (n - j) % 2 == 0);
    if !even {
        return Some(false);
    }
    imaginary_indices(f)?.get(k).copied()
}

/// The assumptions of a root object.
pub fn root_facts(kind: &RootKind) -> Vec<(&'static str, bool)> {
    let mut facts = vec![("algebraic", true), ("complex", true), ("finite", true)];
    if kind.real {
        facts.push(("real", true));
        facts.push(("imaginary", false));
        facts.push(if kind.negative { ("negative", true) } else { ("positive", true) });
    } else {
        facts.push(("real", false));
        facts.push(("imaginary", kind.imaginary));
    }
    facts
}

/// The root object with the given parts and assumptions, without checks.
pub fn node(expr: Expr, index: usize, facts: &[(&str, bool)]) -> Option<Expr> {
    let head = UndefinedFunction::with_facts(HEAD, FactSet::from_kwargs(facts.iter().copied()).ok()?, None).ok()?;
    Some(func(Func::Undefined(head), vec![expr, Expr::int(index as i64)]))
}

/// Root `index` of `expr`, an irreducible polynomial over ZZ in one symbol.
pub fn complex_root(expr: &Expr, index: usize) -> Option<Expr> {
    let (_, f) = integer_poly(expr)?;
    let kind = root_kind(&f, index)?;
    node(expr.clone(), index, &root_facts(&kind))
}

/// The complex conjugate of a non-real root object.
pub fn conjugate(e: &Expr) -> Option<Expr> {
    let (expr, index) = parts(e)?;
    let Kind::Func(Func::Undefined(head), _) = e.kind() else {
        return None;
    };
    match head.facts().get("real") {
        Some(true) => return Some(e.clone()),
        Some(false) => {}
        None => return None,
    }
    let (_, f) = integer_poly(&expr)?;
    let (reals, _) = iso::isolate_real_roots_sqf(&f, None, None, None, false).ok()?;
    let k = index.checked_sub(reals.len())?;
    let other = if k % 2 == 0 { index + 1 } else { index - 1 };
    Some(func(Func::Undefined(head.clone()), vec![expr, Expr::int(other as i64)]))
}

/// The reference's `_integer_basis`: the largest integer `d` with
/// `f(d*y) = m*g(y)` for `g` over ZZ with smaller coefficients.
pub fn integer_basis(f: &[Int]) -> Option<Int> {
    let n = f.len().checked_sub(1)?;
    // Nonzero terms, leading first: (monomial, |coefficient|).
    let mut terms: Vec<(usize, BigInt)> = f
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(k, c)| (n - k, c.to_big().abs()))
        .collect();
    if terms.len() < 2 {
        return None;
    }
    if terms[0].1 < terms[terms.len() - 1].1 {
        let top = terms[0].0;
        terms.reverse();
        for t in &mut terms {
            t.0 = top - t.0;
        }
    } else {
        return None;
    }
    terms.pop();
    if terms.len() == 1 {
        let (m, c) = &terms[0];
        return Int::from_big(c.clone()).exact_root(u32::try_from(*m).ok()?);
    }
    let mut g = BigInt::zero();
    for (_, c) in &terms {
        g = num_integer::Integer::gcd(&g, c);
    }
    let mut divs = crate::ntheory::divisors(&g);
    divs.remove(0);
    'next: for div in divs.into_iter().rev() {
        for (m, c) in &terms {
            let power = num_traits::pow(div.clone(), *m);
            if !(c % &power).is_zero() {
                continue 'next;
            }
        }
        return Some(Int::from_big(div));
    }
    None
}

/// The polynomial `f(d*y)/d**n` of `f` scaled by its integer basis `d`.
pub fn scale_by_basis(f: &[Int], basis: &Int) -> Vec<Int> {
    f.iter()
        .enumerate()
        .map(|(k, c)| c.div_exact_unchecked(&basis.pow(k as u32)))
        .collect()
}

/// The expression of the polynomial with coefficients `f`, leading first.
pub fn poly_expr(f: &[Int], gen: &Expr) -> Expr {
    let n = f.len() - 1;
    let terms: Vec<Expr> = f
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(k, c)| {
            let c = Expr::num(c.to_number());
            match n - k {
                0 => c,
                d => crate::ops::mul(&[c, crate::ops::pow(gen, &Expr::int(d as i64))]),
            }
        })
        .collect();
    crate::ops::add(&terms)
}

/// All roots of the irreducible primitive `f` over ZZ of degree at least
/// two, as the reference's `Poly.all_roots` returns them for a polynomial
/// whose roots are not given by radicals: in the generator `x`, after
/// scaling by the integer basis.
pub fn all_roots_irreducible(f: &[Int]) -> Option<Vec<Expr>> {
    let x = Expr::symbol("x");
    let (coeff, g) = match integer_basis(f) {
        Some(basis) => (Some(basis.clone()), scale_by_basis(f, &basis)),
        None => (None, f.to_vec()),
    };
    let expr = poly_expr(&g, &x);
    let degree = g.len() - 1;
    let (reals, _) = iso::isolate_real_roots_sqf(&g, None, None, None, false).ok()?;
    let nreal = reals.len();
    let imaginary = if nreal < degree { Some(imaginary_indices(&g)?) } else { None };
    let mut roots = Vec::with_capacity(degree);
    for index in 0..degree {
        let kind = if index < nreal {
            RootKind { real: true, negative: reals[index].0.is_negative(), imaginary: false, lower: false, nreal }
        } else {
            let k = index - nreal;
            let flag = *imaginary.as_ref()?.get(k)?;
            RootKind { real: false, negative: false, imaginary: flag, lower: k % 2 == 0, nreal }
        };
        let root = node(expr.clone(), index, &root_facts(&kind))?;
        roots.push(match &coeff {
            Some(c) => crate::ops::mul(&[Expr::num(c.to_number()), root]),
            None => root,
        });
    }
    Some(roots)
}

/// Horner evaluation of `f` (leading first) and its derivative at `z`.
fn horner(f: &[f64], z: num_complex::Complex64) -> (num_complex::Complex64, num_complex::Complex64) {
    let zero = num_complex::Complex64::new(0.0, 0.0);
    let (mut p, mut d) = (zero, zero);
    for &c in f {
        d = d * z + p;
        p = p * z + c;
    }
    (p, d)
}

/// A double precision value of a root object, for numeric ordering only.
///
/// The isolating interval is refined exactly until it is small relative to
/// the root; a complex rectangle is then polished by Newton steps, accepted
/// only when they stay inside the rectangle.
pub fn approx_c64(e: &Expr) -> Option<num_complex::Complex64> {
    use num_complex::Complex64 as C64;
    let (expr, index) = parts(e)?;
    let (_, f) = integer_poly(&expr)?;
    if index + 1 >= f.len() {
        return None;
    }
    let (reals, g) = iso::isolate_real_roots_sqf(&f, None, None, None, false).ok()?;
    let tiny = Rat::new(Int::one(), Int::from(1i64 << 55));
    if index < reals.len() {
        let (s, t) = &reals[index];
        if s == t {
            return Some(C64::new(s.to_f64(), 0.0));
        }
        let mut b = iso::RealBox::from_interval(&g, s, t).ok()?;
        for _ in 0..4096 {
            let (a, c) = (b.a().ok()?, b.b().ok()?);
            if a == c {
                return Some(C64::new(a.to_f64(), 0.0));
            }
            if a.mul(&c) > Rat::zero() {
                let scale = if a.abs() < c.abs() { a.abs() } else { c.abs() };
                if c.sub(&a) < scale.mul(&tiny) {
                    return Some(C64::new(a.add(&c).to_f64() / 2.0, 0.0));
                }
            }
            b = b.refine().ok()?;
        }
        return None;
    }
    let k = index - reals.len();
    let boxes = iso::isolate_complex_roots_sqf(&rats(&f), None, None, None).ok()?;
    let mut b = boxes.into_iter().nth(k)?;
    let imaginary = match e.kind() {
        Kind::Func(Func::Undefined(head), _) => head.facts().get("imaginary") == Some(true),
        _ => false,
    };
    let fc: Vec<f64> = f.iter().map(|c| c.to_f64()).collect();
    let coarse = Rat::new(Int::one(), Int::from(1i64 << 12));
    for step in 0..4096 {
        let (ax, bx, ay, by) = (b.ax(), b.bx(), b.ay(), b.by());
        let big = |u: &Rat, v: &Rat| if u.abs() < v.abs() { v.abs() } else { u.abs() };
        let scale = {
            let sx = big(&ax, &bx);
            let sy = big(&ay, &by);
            if sx < sy { sy } else { sx }
        };
        let (dx, dy) = (b.dx(), b.dy());
        let center = C64::new(ax.add(&bx).to_f64() / 2.0, ay.add(&by).to_f64() / 2.0);
        let done = |z: C64| if imaginary { C64::new(0.0, z.im) } else { z };
        if dx < scale.mul(&tiny) && dy < scale.mul(&tiny) {
            return Some(done(center));
        }
        if step % 8 == 0 && dx < scale.mul(&coarse) && dy < scale.mul(&coarse) {
            let mut z = center;
            for _ in 0..60 {
                let (p, d) = horner(&fc, z);
                if d == C64::new(0.0, 0.0) || !p.is_finite() {
                    break;
                }
                let next = z - p / d;
                if (next - z).norm() <= next.norm() * 4.0 * f64::EPSILON {
                    z = next;
                    break;
                }
                z = next;
            }
            let (lx, hx, ly, hy) = (ax.to_f64(), bx.to_f64(), ay.to_f64(), by.to_f64());
            if z.is_finite() && lx <= z.re && z.re <= hx && ly <= z.im && z.im <= hy {
                let (p, d) = horner(&fc, z);
                let step = if d == C64::new(0.0, 0.0) { f64::INFINITY } else { (p / d).norm() };
                if step <= z.norm() * 64.0 * f64::EPSILON {
                    return Some(done(z));
                }
            }
        }
        b = b.refine().ok()?;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ints(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&c| Int::from(c)).collect()
    }

    #[test]
    fn integer_basis_matches_reference() {
        assert_eq!(integer_basis(&ints(&[1, 0, 0, 0, 512, 1024])), Some(Int::from(4)));
        assert_eq!(integer_basis(&ints(&[1, 0, 0, 0, 16])), Some(Int::from(2)));
        assert_eq!(integer_basis(&ints(&[1, 0, 0, 0, -1, -1])), None);
        assert_eq!(scale_by_basis(&ints(&[1, 0, 0, 0, 512, 1024]), &Int::from(4)), ints(&[1, 0, 0, 0, 2, 1]));
    }

    #[test]
    fn root_kinds_follow_isolation_order() {
        let f = ints(&[1, 0, 0, 0, -1, -1]);
        let k = root_kind(&f, 0).unwrap();
        assert!(k.real && !k.negative && k.nreal == 1);
        let k = root_kind(&f, 1).unwrap();
        assert!(!k.real && k.lower && !k.imaginary);
        assert!(!root_kind(&f, 2).unwrap().lower);
        let g = ints(&[1, 0, 3, 0, 1]);
        assert!((0..4).all(|i| root_kind(&g, i).unwrap().imaginary));
        let h = ints(&[1, 0, 0, 0, 1]);
        assert!((0..4).all(|i| !root_kind(&h, i).unwrap().imaginary));
    }
}
