//! Ground domain detection for matrices (SymPy's `construct_domain` as
//! used by `Matrix.to_DM()`), and conversion of matrices to integer /
//! polynomial representations for the exact algorithms.

use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::ratfunc::{PolyCtx, QPoly};
use super::Matrix;
use crate::expr::{Constant, Expr, Kind};
use crate::number::Number;

/// The domain of a matrix' entries.
#[derive(Clone, Debug)]
pub enum Domain {
    /// Integers.
    ZZ,
    /// Rationals.
    QQ,
    /// Polynomials `ZZ[gens]` (`qq == false`) or `QQ[gens]`.
    Poly { ctx: PolyCtx, qq: bool },
    /// Rational functions `ZZ(gens)` / `QQ(gens)`.
    Frac { ctx: PolyCtx, qq: bool },
    /// Floating point (`RR`, or `RR[...]` / `CC`...): name for printing.
    Inexact(String),
    /// Gaussian rationals and other domains handled generically.
    Other(String),
    /// Generic expressions (`EXRAW`).
    EX,
}

impl Domain {
    pub fn is_rational(&self) -> bool {
        matches!(self, Domain::ZZ | Domain::QQ)
    }

    pub fn is_ex(&self) -> bool {
        matches!(self, Domain::EX)
    }

    /// SymPy's printed name of the domain.
    pub fn name(&self) -> String {
        fn gens(ctx: &PolyCtx) -> String {
            ctx.gens.iter().map(|g| g.to_string()).collect::<Vec<_>>().join(",")
        }
        match self {
            Domain::ZZ => "ZZ".into(),
            Domain::QQ => "QQ".into(),
            Domain::Poly { ctx, qq } => format!("{}[{}]", if *qq { "QQ" } else { "ZZ" }, gens(ctx)),
            Domain::Frac { ctx, qq } => format!("{}({})", if *qq { "QQ" } else { "ZZ" }, gens(ctx)),
            Domain::Inexact(s) | Domain::Other(s) => s.clone(),
            Domain::EX => "EX".into(),
        }
    }
}

/// Is `g` (a generator candidate without free symbols) an algebraic number?
fn is_algebraic_number(g: &Expr) -> bool {
    match g.kind() {
        Kind::Num(n) => n.is_rational(),
        Kind::Const(Constant::I) => true,
        Kind::Pow(b, x) => x.is_rational_number() && is_algebraic_number(b),
        Kind::Add(_, terms) => terms.iter().all(|(t, _)| is_algebraic_number(t)),
        Kind::Mul(_, facs) => facs
            .iter()
            .all(|(b, x)| x.is_rational_number() && is_algebraic_number(b)),
        _ => false,
    }
}

fn contains_float(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_float(),
        _ => e.args().iter().any(contains_float),
    }
}

fn contains_i(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::I) => true,
        _ => e.args().iter().any(contains_i),
    }
}

/// Determine the domain of a list of expressions.
pub fn construct_domain(items: &[Expr]) -> Domain {
    // simple: all numbers
    if items.iter().all(|e| e.is_number()) {
        let mut rat = false;
        let mut float = false;
        for e in items {
            match e.as_number().unwrap() {
                Number::Float(_) => float = true,
                Number::Rational(_) => rat = true,
                _ => {}
            }
        }
        return if float {
            Domain::Inexact("RR".into())
        } else if rat {
            Domain::QQ
        } else {
            Domain::ZZ
        };
    }
    if items.iter().any(contains_float) {
        return Domain::Inexact("RR".into());
    }
    if items.iter().any(contains_i) {
        return Domain::Other("ZZ_I".into());
    }
    if items.iter().all(|e| e.mask() == 0) {
        // numeric (algebraic or transcendental) constants only: EX unless
        // all of them are usable as generators
    }
    let mut parts = Vec::with_capacity(2 * items.len());
    for e in items {
        let (n, d) = super::ratfunc::as_numer_denom(e);
        parts.push(n);
        parts.push(d);
    }
    let Some((ctx, polys)) = PolyCtx::from_exprs(&parts, false) else {
        return Domain::EX;
    };
    // generators that are algebraic numbers or share symbols -> EX
    let mut syms: Vec<Expr> = Vec::new();
    for g in &ctx.gens {
        let fs = crate::traverse::free_symbols(g);
        if fs.is_empty() && is_algebraic_number(g) {
            return Domain::EX;
        }
        // Closed function values need not be algebraically independent:
        // sin(3)^2 + cos(3)^2 = 1. A formal polynomial field can invert an
        // identically singular matrix unless expression pivots see this relation.
        if fs.is_empty() && matches!(g.kind(),Kind::Func(..)) {
            return Domain::EX;
        }
        for s in fs {
            if syms.contains(&s) {
                return Domain::EX;
            }
            syms.push(s);
        }
    }
    let fractions = polys.chunks(2).any(|nd| !nd[1].p.is_constant());
    if !fractions {
        // coefficients are numerators divided by the constant denominators
        let mut qq = false;
        for nd in polys.chunks(2) {
            let d = Rat::new(nd[1].p.lc(), nd[1].den.clone());
            let scale = Rat::new(Int::ONE, nd[0].den.clone()).div(&d);
            for (_, c) in &nd[0].p.terms {
                if !Rat::from_int(c.clone()).mul(&scale).is_integer() {
                    qq = true;
                }
            }
        }
        Domain::Poly { ctx, qq }
    } else {
        let qq = polys.iter().any(|q| !q.den.is_one());
        Domain::Frac { ctx, qq }
    }
}

pub fn matrix_domain(m: &Matrix) -> Domain {
    construct_domain(m.data())
}

/// A structural proof that every entry belongs to a rational polynomial
/// ring in ordinary, commutative symbols. Unknown expressions deliberately
/// keep the general domain-construction path: in particular, function values
/// can have algebraic relations requiring EXRAW coefficient simplification.
pub fn has_rational_polynomial_entries(m: &Matrix) -> bool {
    use rustc_hash::FxHashMap;
    fn exponent(e: &Expr) -> bool {
        e.as_number().is_some_and(|n| n.is_integer() && !n.is_negative())
    }
    fn polynomial(e: &Expr, seen: &mut FxHashMap<Expr, bool>) -> bool {
        if let Some(result) = seen.get(e) {
            return *result;
        }
        let result = match e.kind() {
            Kind::Num(n) => n.is_rational(),
            Kind::Sym(s) => s.host().is_none() && s.scalar_identity().is_none() && crate::decompose::commutes(e),
            Kind::Add(c, terms) => c.is_rational() && terms.iter().all(|(t, c)| {
                c.is_rational() && polynomial(t, seen)
            }),
            Kind::Mul(c, factors) => c.is_rational() && factors.iter().all(|(b, p)| {
                exponent(p) && polynomial(b, seen)
            }),
            Kind::Pow(b, p) => exponent(p) && polynomial(b, seen),
            _ => false,
        };
        seen.insert(e.clone(), result);
        result
    }
    let mut seen = FxHashMap::default();
    m.data().iter().all(|e| polynomial(e, &mut seen))
}

/// Exact rational entries need no floating-point or composite-domain adapter.
pub fn has_rational_entries(m: &Matrix) -> bool {
    m.data().iter().all(Expr::is_rational_number)
}

// ----------------------------------------------------------------------
// integer / rational / polynomial matrices
// ----------------------------------------------------------------------

/// Dense matrix of exact rationals.
#[derive(Clone, Debug, PartialEq)]
pub struct RatMat {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Rat>,
}

impl RatMat {
    pub fn from_matrix(m: &Matrix) -> Option<RatMat> {
        let mut data = Vec::with_capacity(m.len());
        for e in m.data() {
            data.push(Rat::from_number(e.as_number()?)?);
        }
        Some(RatMat { rows: m.rows(), cols: m.cols(), data })
    }

    pub fn to_matrix(&self) -> Matrix {
        Matrix::new_unchecked(
            self.rows,
            self.cols,
            self.data.iter().map(|r| Expr::num(r.to_number())).collect(),
        )
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &Rat {
        &self.data[i * self.cols + j]
    }

    pub fn is_integer(&self) -> bool {
        self.data.iter().all(|r| r.den.is_one())
    }

    /// Scale each row to integers: returns `(integer rows, row scales)`
    /// with `A' = diag(scales) * A`.
    pub fn to_int_rows(&self) -> (Vec<Int>, Vec<Int>) {
        let mut out = Vec::with_capacity(self.data.len());
        let mut scales = Vec::with_capacity(self.rows);
        for i in 0..self.rows {
            let mut l = Int::ONE;
            for j in 0..self.cols {
                l = l.lcm(&self.get(i, j).den);
            }
            for j in 0..self.cols {
                let r = self.get(i, j);
                out.push(if l.is_one() { r.num.clone() } else { r.num.mul(&l.divexact(&r.den)) });
            }
            scales.push(l);
        }
        (out, scales)
    }

    pub fn matmul(&self, o: &RatMat) -> RatMat {
        let (n, m, p) = (self.rows, self.cols, o.cols);
        let mut data = Vec::with_capacity(n * p);
        if self.is_integer() && o.is_integer() {
            for i in 0..n {
                for j in 0..p {
                    let mut s = Int::ZERO;
                    for k in 0..m {
                        let a = &self.get(i, k).num;
                        if a.is_zero() {
                            continue;
                        }
                        s = s.add(&a.mul(&o.get(k, j).num));
                    }
                    data.push(Rat::from_int(s));
                }
            }
        } else {
            for i in 0..n {
                for j in 0..p {
                    let mut s = Rat::zero();
                    for k in 0..m {
                        let a = self.get(i, k);
                        if a.is_zero() {
                            continue;
                        }
                        s = s.add(&a.mul(o.get(k, j)));
                    }
                    data.push(s);
                }
            }
        }
        RatMat { rows: n, cols: p, data }
    }

    pub fn identity(n: usize) -> RatMat {
        let mut data = vec![Rat::zero(); n * n];
        for i in 0..n {
            data[i * n + i] = Rat::one();
        }
        RatMat { rows: n, cols: n, data }
    }

    pub fn pow(&self, mut k: u64) -> RatMat {
        let mut r = RatMat::identity(self.rows);
        let mut b = self.clone();
        let mut first = true;
        while k > 0 {
            if k & 1 == 1 {
                r = if first { b.clone() } else { r.matmul(&b) };
                first = false;
            }
            k >>= 1;
            if k > 0 {
                b = b.matmul(&b);
            }
        }
        r
    }
}

/// Matrix of polynomials over a common context.
#[derive(Clone, Debug)]
pub struct PolyMat {
    pub rows: usize,
    pub cols: usize,
    /// Integer polynomial entries.
    pub data: Vec<MPoly>,
    pub ctx: PolyCtx,
}

impl PolyMat {
    /// Convert a matrix whose entries are polynomials (`Domain::Poly`).
    /// Rows are scaled to integer coefficients: returns the scaled matrix
    /// and the row scales (`A' = diag(scales) * A`).
    pub fn from_matrix(m: &Matrix, ctx: &PolyCtx) -> Option<(PolyMat, Vec<Int>)> {
        let mut qs: Vec<QPoly> = Vec::with_capacity(m.len());
        for e in m.data() {
            qs.push(ctx.to_qpoly(e)?);
        }
        let mut data = Vec::with_capacity(m.len());
        let mut scales = Vec::with_capacity(m.rows());
        for i in 0..m.rows() {
            let mut l = Int::ONE;
            for j in 0..m.cols() {
                l = l.lcm(&qs[i * m.cols() + j].den);
            }
            for j in 0..m.cols() {
                let q = &qs[i * m.cols() + j];
                data.push(q.p.mul_int(&l.divexact(&q.den)));
            }
            scales.push(l);
        }
        Some((PolyMat { rows: m.rows(), cols: m.cols(), data, ctx: ctx.clone() }, scales))
    }

    /// Convert a rational function matrix (`Domain::Frac`): every row is
    /// multiplied by the product of the distinct denominators in it.
    /// Returns the polynomial matrix and the row multipliers.
    pub fn from_frac_matrix(m: &Matrix, ctx: &PolyCtx) -> Option<(PolyMat, Vec<MPoly>)> {
        let n = ctx.nvars();
        let mut data = Vec::with_capacity(m.len());
        let mut scales = Vec::with_capacity(m.rows());
        for i in 0..m.rows() {
            let mut nums = Vec::with_capacity(m.cols());
            let mut dens: Vec<MPoly> = Vec::with_capacity(m.cols());
            for j in 0..m.cols() {
                let (a, b) = super::ratfunc::as_numer_denom(m.get(i, j));
                let pa = ctx.to_qpoly(&a)?;
                let pb = ctx.to_qpoly(&b)?;
                // a/b = (pa.p/pa.den)/(pb.p/pb.den) = (pa.p*pb.den)/(pb.p*pa.den)
                let num = pa.p.mul_int(&pb.den);
                let den = pb.p.mul_int(&pa.den);
                nums.push(num);
                dens.push(den);
            }
            // row multiplier: lcm-ish product of distinct denominators
            let mut l = MPoly::one(n);
            for d in &dens {
                if l.divexact(d).is_none() {
                    let g = super::mpoly::gcd(&l, d);
                    l = l.mul(&d.divexact(&g).unwrap());
                }
            }
            for (num, den) in nums.iter().zip(&dens) {
                data.push(num.mul(&l.divexact(den).unwrap()));
            }
            scales.push(l);
        }
        Some((PolyMat { rows: m.rows(), cols: m.cols(), data, ctx: ctx.clone() }, scales))
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &MPoly {
        &self.data[i * self.cols + j]
    }
}
