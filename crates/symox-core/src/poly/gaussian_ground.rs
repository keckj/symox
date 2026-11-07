//! Polynomials over the composite Gaussian domains `ZZ_I[y...]`,
//! `QQ_I[y...]`, `ZZ_I(y...)` and `QQ_I(y...)`.
//!
//! These polynomials are stored as `EX` polynomials (`ExpressionPoly`)
//! carrying a [`GaussianGround`] marker. Every coefficient is kept in the
//! normal form of the reference ground element: a polynomial in the domain
//! symbols with Gaussian coefficients, or for a fraction field the quotient
//! of two such polynomials reduced as `PolyElement.cancel` does. Arithmetic
//! goes through `EX` and is renormalized afterwards; the algorithms whose
//! result depends on the ring structure (contents, gcds, lcms,
//! factorization, square-free decomposition) follow the reference dense
//! algorithms over the composite domain on a flattened `QQ<I>` polynomial
//! in the generators followed by the symbols.
use super::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN};
use super::convert::{PResult, PolyError};
use super::domains::ExactDomain;
use super::expression_domain::ExpressionPoly;
use super::gaussian as gi;
use super::int::{Coeff, Field, Rat};
use super::mpoly::MPoly;
use crate::expr::Expr;
use std::collections::BTreeMap;
use std::sync::Arc;

type P = MPoly<AN>;

/// The composite Gaussian ground of a polynomial: `ZZ_I` (`integral`) or
/// `QQ_I`, as a polynomial ring (`ring`) or fraction field over `symbols`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GaussianGround {
    pub integral: bool,
    pub ring: bool,
    pub symbols: Vec<Expr>,
    /// `RealField(p)`/`ComplexField(p)` in place of the Gaussian ground
    /// (`integral` is then false).
    pub inexact: Option<super::inexact_ground::Inexact>,
    /// An algebraic field `QQ<a>` in place of the Gaussian ground
    /// (`integral` is then false): `QQ<a>[y...]` and `QQ<a>(y...)`.
    pub algebraic: Option<Arc<NumberField>>,
}

fn scope() -> AlgebraicScope {
    AlgebraicScope::new(gi::field())
}

fn unit_power(k: u32, q: Rat) -> AN {
    match k % 4 {
        0 => gi::make(q, Rat::zero()),
        1 => gi::make(Rat::zero(), q),
        2 => gi::make(q.neg(), Rat::zero()),
        _ => gi::make(Rat::zero(), q.neg()),
    }
}

/// `e` as a polynomial with Gaussian rational coefficients in `vars`.
pub fn gaussian_poly(e: &Expr, vars: &[Expr]) -> PResult<P> {
    let _s = scope();
    let n = vars.len();
    let fail = || PolyError::Failed(format!("{e} is not a Gaussian polynomial"));
    let c = super::convert::from_exprs_gens(std::slice::from_ref(e), vars).map_err(|_| fail())?;
    if c.float {
        return Err(fail());
    }
    let extra = &c.gens[n..];
    if extra.len() > 1 || extra.first().is_some_and(|g| *g != Expr::i()) {
        return Err(fail());
    }
    let terms = c.polys[0]
        .iter()
        .map(|(x, q)| (x[..n].to_vec(), unit_power(x.get(n).copied().unwrap_or(0), q.clone())))
        .collect();
    Ok(MPoly::from_terms(n, terms))
}

/// The reference class of an expression, for `CoercionFailed` messages.
pub(crate) fn class_path(e: &Expr) -> &'static str {
    use crate::expr::Kind;
    match e.kind() {
        Kind::Num(n) if n.is_integer() => "sympy.core.numbers.Integer",
        Kind::Num(n) if n.is_float() => "sympy.core.numbers.Float",
        Kind::Num(_) => "sympy.core.numbers.Rational",
        Kind::Sym(_) => "sympy.core.symbol.Symbol",
        Kind::Add(..) => "sympy.core.add.Add",
        Kind::Mul(..) => "sympy.core.mul.Mul",
        Kind::Pow(..) => "sympy.core.power.Pow",
        _ if *e == Expr::i() => "sympy.core.numbers.ImaginaryUnit",
        _ => "sympy.core.expr.Expr",
    }
}

/// Monic gcd of two polynomials over `QQ_I`.
fn field_poly_gcd(a: &P, b: &P) -> PResult<P> {
    super::algebraic_multi::gcd(a, b)
}

fn exact(a: &P, b: &P) -> PResult<P> {
    super::algebraic_multi::exact(a, b)
}

/// `PolyElement.cofactors` over `QQ<a>`: the zero and monomial shortcuts
/// (`QQ<a>.gcd` is one, `is_nonnegative` the sign of the leading
/// coordinate), otherwise the monic `dmp_inner_gcd` over the field.
fn algebraic_cofactors(f: &P, g: &P) -> PResult<(P, P, P)> {
    let n = f.nvars;
    let (zero, one) = (MPoly::zero(n), MPoly::one(n));
    let zero_gcd = |g: &P| if g.lc().is_negative() { (g.neg(), zero.clone(), one.neg()) } else { (g.clone(), zero.clone(), one.clone()) };
    if f.is_zero() && g.is_zero() {
        return Ok((zero.clone(), zero.clone(), zero));
    }
    if f.is_zero() {
        return Ok(zero_gcd(g));
    }
    if g.is_zero() {
        let (h, cfg, cff) = zero_gcd(f);
        return Ok((h, cff, cfg));
    }
    let monom = |f: &P, g: &P| -> PResult<(P, P, P)> {
        let mut mgcd = f.lm().to_vec();
        for (e, _) in g.iter() {
            for (x, &y) in mgcd.iter_mut().zip(e) {
                *x = (*x).min(y);
            }
        }
        let h = MPoly::monomial(&mgcd, AN::one());
        Ok((h.clone(), exact(f, &h)?, exact(g, &h)?))
    };
    if f.len() == 1 {
        return monom(f, g);
    }
    if g.len() == 1 {
        let (h, cfg, cff) = monom(g, f)?;
        return Ok((h, cff, cfg));
    }
    let h = field_poly_gcd(f, g)?;
    Ok((h.clone(), exact(f, &h)?, exact(g, &h)?))
}

/// Drop the first `k` variables of a polynomial that does not involve them.
fn drop_vars(p: &P, k: usize) -> P {
    MPoly::from_terms(p.nvars - k, p.iter().map(|(e, c)| (e[k..].to_vec(), c.clone())).collect())
}

/// Prepend `k` absent variables.
fn lift(p: &P, k: usize) -> P {
    MPoly::from_terms(
        p.nvars + k,
        p.iter().map(|(e, c)| (std::iter::repeat_n(0, k).chain(e.iter().copied()).collect(), c.clone())).collect(),
    )
}

/// Dense coefficients in the first variable, leading first, each still in
/// `p.nvars` variables with the first one absent.
fn dense(p: &P) -> Vec<P> {
    if p.is_zero() {
        return Vec::new();
    }
    let groups = p.coeffs_in(0);
    let top = groups[0].0;
    let mut groups = groups.into_iter().peekable();
    (0..=top)
        .rev()
        .map(|k| {
            if groups.peek().is_some_and(|(d, _)| *d == k) {
                groups.next().unwrap().1
            } else {
                MPoly::zero(p.nvars)
            }
        })
        .collect()
}

/// The leading coefficient in the first `k` variables, over the rest.
fn lc_ground(p: &P, k: usize) -> P {
    if p.is_zero() {
        return MPoly::zero(p.nvars - k);
    }
    let lm = p.lm()[..k].to_vec();
    MPoly::from_terms(p.nvars - k, p.iter().filter(|(e, _)| e[..k] == lm[..]).map(|(e, c)| (e[k..].to_vec(), c.clone())).collect())
}

impl GaussianGround {
    pub fn domain(&self) -> ExactDomain {
        let ground = match self.inexact {
            Some(k) if k.precision == 53 && !k.complex => ExactDomain::RR,
            Some(k) => ExactDomain::Inexact { precision: k.precision, complex: k.complex },
            None => match &self.algebraic {
                Some(k) => ExactDomain::Algebraic(k.clone()),
                None => ExactDomain::Gaussian { integral: self.integral },
            },
        };
        if self.symbols.is_empty() {
            // the plain inexact field `CC` (or `ComplexField(prec)`)
            return ground;
        }
        ExactDomain::Composite {
            ground: Arc::new(ground),
            ring: self.ring,
            parameters: self.symbols.clone(),
            order: "lex".into(),
        }
    }
    pub fn display(&self) -> String {
        self.domain().to_string()
    }
    pub fn to_field(&self) -> Self {
        Self { ring: false, ..self.clone() }
    }
    pub fn to_ring(&self) -> Self {
        Self { ring: true, ..self.clone() }
    }
    /// The plain field `ComplexField(precision)` as a ground without
    /// symbols.
    pub fn complex_field(precision: u32) -> Self {
        Self { integral: false, ring: true, symbols: Vec::new(), inexact: Some(super::inexact_ground::Inexact { precision, complex: true }), algebraic: None }
    }
    /// Whether this is the plain field `CC` (no domain symbols).
    pub fn is_scalar(&self) -> bool {
        self.symbols.is_empty() && self.inexact.is_some()
    }
    pub fn from_domain(d: &ExactDomain) -> Option<Self> {
        let ExactDomain::Composite { ground, ring, parameters, .. } = d else { return None };
        if parameters.is_empty() {
            return None;
        }
        let inexact = match ground.as_ref() {
            ExactDomain::RR => Some(super::inexact_ground::Inexact { precision: 53, complex: false }),
            ExactDomain::Inexact { precision, complex } => Some(super::inexact_ground::Inexact { precision: *precision, complex: *complex }),
            _ => None,
        };
        if let Some(k) = inexact {
            return Some(Self { integral: false, ring: *ring, symbols: parameters.clone(), inexact: Some(k), algebraic: None });
        }
        if let ExactDomain::Algebraic(k) = ground.as_ref() {
            if gi::is_gaussian_field(k) {
                return None;
            }
            return Some(Self { integral: false, ring: *ring, symbols: parameters.clone(), inexact: None, algebraic: Some(k.clone()) });
        }
        let ExactDomain::Gaussian { integral } = ground.as_ref() else { return None };
        Some(Self { integral: *integral, ring: *ring, symbols: parameters.clone(), inexact: None, algebraic: None })
    }
    /// The coefficient field: `QQ<I>` for the Gaussian grounds.
    pub fn field(&self) -> Arc<NumberField> {
        self.algebraic.clone().unwrap_or_else(gi::field)
    }
    pub(crate) fn scope(&self) -> AlgebraicScope {
        AlgebraicScope::new(self.field())
    }
    /// `e` as a polynomial in the symbols over the coefficient field.
    fn parse(&self, e: &Expr) -> PResult<P> {
        match &self.algebraic {
            None => gaussian_poly(e, &self.symbols),
            Some(k) => Ok(AlgebraicPoly::from_expr(e, &self.symbols, k.clone())?.rep),
        }
    }
    /// `K.canonical_unit`: a unit of the Gaussian ground, or the sign of
    /// the leading coordinate (`AlgebraicField.is_negative`).
    fn unit(&self, a: &AN) -> AN {
        if self.algebraic.is_some() {
            if a.is_negative() { AN::one().neg() } else { AN::one() }
        } else {
            gi::canonical_unit(a)
        }
    }
    /// `PolyElement.cancel` of two ground ring elements followed by the
    /// canonical unit of the denominator, as `FracElement.new` does.
    pub fn cancel_pair(&self, f: &P, g: &P) -> PResult<(P, P)> {
        if self.algebraic.is_none() {
            return gi::expr_cancel(f, g);
        }
        if g.is_zero() {
            return Err(PolyError::Reference { kind: "ZeroDivisionError", message: String::new() });
        }
        let _s = self.scope();
        if f.is_zero() {
            return Ok((f.clone(), MPoly::one(f.nvars)));
        }
        // QQ<a> has no associated ring: `_, p, q = f.cofactors(g)`.
        let (_, p, q) = algebraic_cofactors(f, g)?;
        if q.lc().is_negative() {
            Ok((p.neg(), q.neg()))
        } else {
            Ok((p, q))
        }
    }
    /// `FracField.from_expr` over `QQ<a>`: the reference `_rebuild_expr`,
    /// whose `FracElement` operations cancel at each step (and whose
    /// powers do not), so the stored numerator and denominator are those of
    /// the reference element.
    fn rebuild(&self, e: &Expr) -> PResult<(P, P)> {
        let m = self.m();
        let k = self.algebraic.clone().expect("algebraic ground");
        #[derive(Clone)]
        enum R {
            Ground(AN),
            Frac(P, P),
        }
        let frac = |r: R| -> (P, P) {
            match r {
                R::Ground(c) => (MPoly::constant(c, m), MPoly::one(m)),
                R::Frac(p, q) => (p, q),
            }
        };
        fn walk(g: &GaussianGround, e: &Expr, k: &Arc<NumberField>, m: usize) -> PResult<R> {
            if let Some(i) = g.symbols.iter().position(|s| s == e) {
                let mut x = vec![0u32; m];
                x[i] = 1;
                return Ok(R::Frac(MPoly::monomial(&x, AN::one()), MPoly::one(m)));
            }
            let as_frac = |r: R| -> (P, P) {
                match r {
                    R::Ground(c) => (MPoly::constant(c, m), MPoly::one(m)),
                    R::Frac(p, q) => (p, q),
                }
            };
            if e.is_add() || e.is_mul() {
                let add = e.is_add();
                let mut acc: Option<R> = None;
                for a in e.args() {
                    let b = walk(g, &a, k, m)?;
                    acc = Some(match acc {
                        None => b,
                        Some(R::Ground(x)) if matches!(b, R::Ground(_)) => {
                            let R::Ground(y) = b else { unreachable!() };
                            R::Ground(if add { x.add(&y) } else { x.mul(&y) })
                        }
                        Some(R::Frac(p, q)) if matches!(b, R::Frac(..)) => {
                            let R::Frac(c, d) = b else { unreachable!() };
                            if add {
                                // FracElement.__add__: zero operands are
                                // returned as they are.
                                if c.is_zero() {
                                    R::Frac(p, q)
                                } else if p.is_zero() {
                                    R::Frac(c, d)
                                } else if q == d {
                                    let (x, y) = g.cancel_pair(&p.add(&c), &q)?;
                                    R::Frac(x, y)
                                } else {
                                    let (x, y) = g.cancel_pair(&p.mul(&d).add(&q.mul(&c)), &q.mul(&d))?;
                                    R::Frac(x, y)
                                }
                            } else if p.is_zero() || c.is_zero() {
                                R::Frac(MPoly::zero(m), MPoly::one(m))
                            } else {
                                let (x, y) = g.cancel_pair(&p.mul(&c), &q.mul(&d))?;
                                R::Frac(x, y)
                            }
                        }
                        Some(other) => {
                            // a ground element with a fraction, in either
                            // order: `FracElement.__radd__`/`__rmul__`.
                            let (c, (p, q)) = match (other, b) {
                                (R::Ground(c), f) | (f, R::Ground(c)) => (c, as_frac(f)),
                                _ => unreachable!(),
                            };
                            if !add && (p.is_zero() || c.is_zero()) {
                                R::Frac(MPoly::zero(m), MPoly::one(m))
                            } else {
                                let numer = if add { p.add(&q.scale(&c)) } else { p.scale(&c) };
                                let (x, y) = g.cancel_pair(&numer, &q)?;
                                R::Frac(x, y)
                            }
                        }
                    });
                }
                return acc.ok_or_else(|| g.coercion(e));
            }
            if e.is_pow() {
                let args = e.args();
                if let Some(n) = args[1].as_i64().filter(|n| *n != 1) {
                    let n32 = u32::try_from(n.unsigned_abs()).map_err(|_| g.coercion(e))?;
                    return Ok(match walk(g, &args[0], k, m)? {
                        R::Ground(c) => {
                            if c.is_zero() && n < 0 {
                                return Err(g.coercion(e));
                            }
                            R::Ground(c.pow(n))
                        }
                        R::Frac(p, q) => {
                            // FracElement.__pow__: raw_new, no cancellation
                            if n >= 0 {
                                R::Frac(p.pow(n32), q.pow(n32))
                            } else if p.is_zero() {
                                return Err(PolyError::Reference { kind: "ZeroDivisionError", message: String::new() });
                            } else {
                                R::Frac(q.pow(n32), p.pow(n32))
                            }
                        }
                    });
                }
            }
            let inverse = crate::ops::div(&Expr::one(), e);
            if let Some(i) = g.symbols.iter().position(|s| *s == inverse) {
                let mut x = vec![0u32; m];
                x[i] = 1;
                return Ok(R::Frac(MPoly::one(m), MPoly::monomial(&x, AN::one())));
            }
            if !crate::traverse::free_symbols(e).is_empty() {
                return Err(g.coercion(e));
            }
            AN::from_expr(e, k.clone()).map(R::Ground).map_err(|_| g.coercion(e))
        }
        let _s = self.scope();
        let r = walk(self, e, &k, m)?;
        Ok(match r {
            // FracField.field_new of a ground element: `new(c)`
            R::Ground(c) => self.cancel_pair(&MPoly::constant(c, m), &MPoly::one(m))?,
            r => frac(r),
        })
    }
    fn m(&self) -> usize {
        self.symbols.len()
    }
    pub(crate) fn coercion_error(&self, e: &Expr) -> PolyError {
        self.coercion(e)
    }
    fn coercion(&self, e: &Expr) -> PolyError {
        PolyError::Reference {
            kind: "CoercionFailed",
            message: format!("Cannot convert {e} of type <class '{}'> to {}", class_path(e), self.display()),
        }
    }

    /// The normal form `(p, q)` of the ground element `e`; `q` is one for a
    /// polynomial ring.
    pub fn element(&self, e: &Expr) -> PResult<(P, P)> {
        if self.algebraic.is_some() && !self.ring {
            return self.rebuild(e);
        }
        let _s = self.scope();
        let m = self.m();
        let (n, d) = super::api::as_numer_denom(e);
        let num = self.parse(&n).map_err(|_| self.coercion(e))?;
        let den = self.parse(&d).map_err(|_| self.coercion(e))?;
        if den.is_zero() {
            return Err(self.coercion(e));
        }
        if num.is_zero() {
            return Ok((MPoly::zero(m), MPoly::one(m)));
        }
        if self.ring {
            let (p, q) = if den.is_constant() { (num, den) } else { self.cancel_pair(&num, &den)? };
            if !q.is_constant() {
                return Err(self.coercion(e));
            }
            let p = p.scale(&q.lc().inv());
            if self.integral && !gi::poly_is_integral(&p) {
                return Err(self.coercion(e));
            }
            return Ok((p, MPoly::one(m)));
        }
        gi::expr_cancel(&num, &den)
    }
    fn poly_expr(&self, p: &P) -> Expr {
        let gaussian = if self.algebraic.is_some() { None } else { Some(false) };
        AlgebraicPoly { rep: p.clone(), gens: self.symbols.clone(), number_field: self.field(), gaussian }.as_expr()
    }
    pub fn element_expr(&self, p: &P, q: &P) -> Expr {
        let _algebraic_scope = self.scope();
        if q.is_one() {
            self.poly_expr(p)
        } else {
            crate::ops::div(&self.poly_expr(p), &self.poly_expr(q))
        }
    }
    pub(crate) fn reduce(&self, e: &Expr) -> PResult<Expr> {
        let _algebraic_scope = self.scope();
        let (p, q) = self.element(e)?;
        Ok(self.element_expr(&p, &q))
    }
    /// Convert every coefficient of `p` into this ground, failing with the
    /// reference `CoercionFailed`.
    pub fn convert_poly(&self, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::convert_poly(self, p);
        }
        let mut terms = Vec::new();
        for (m, c) in p.terms() {
            let (a, b) = self.element(&c)?;
            if !a.is_zero() {
                terms.push((m, self.element_expr(&a, &b)));
            }
        }
        let mut out = ExpressionPoly::from_terms(terms, &p.gens)?;
        out.ground = Some(self.clone());
        Ok(out)
    }
    pub fn convert_expr(&self, e: &Expr, gens: &[Expr]) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        let p = ExpressionPoly::from_expr(e, gens)?;
        self.convert_poly(&p)
    }

    // ------------------------------------------------------------------
    // flattened representation
    // ------------------------------------------------------------------

    /// `(num, den)`: `p = num/den` with `num` over the generators followed
    /// by the symbols and `den` the lcm of the coefficient denominators.
    pub fn flatten(&self, p: &ExpressionPoly) -> PResult<(P, P)> {
        let _algebraic_scope = self.scope();
        let _s = self.scope();
        let (k, m) = (p.gens.len(), self.m());
        let mut elements = Vec::new();
        let mut den = MPoly::one(m);
        for (mono, c) in p.terms() {
            let (a, b) = self.element(&c)?;
            if !b.is_one() {
                den = if self.algebraic.is_some() { self.to_ring().elcm(&den, &b)? } else { gi::ring_lcm(&den, &b)? };
            }
            elements.push((mono, a, b));
        }
        let mut terms = Vec::new();
        for (mono, a, b) in elements {
            let a = if den.is_one() { a } else { a.mul(&exact(&den, &b)?) };
            for (e, c) in a.iter() {
                let mut x = mono.clone();
                x.extend_from_slice(e);
                terms.push((x, c.clone()));
            }
        }
        Ok((MPoly::from_terms(k + m, terms), den))
    }
    /// Inverse of [`flatten`].
    pub fn eject(&self, num: &P, den: &P, gens: &[Expr]) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        let (k, m) = (gens.len(), self.m());
        let mut groups: BTreeMap<Vec<u32>, Vec<(Vec<u32>, AN)>> = BTreeMap::new();
        for (e, c) in num.iter() {
            groups.entry(e[..k].to_vec()).or_default().push((e[k..].to_vec(), c.clone()));
        }
        let terms = groups
            .into_iter()
            .map(|(mono, ts)| (mono, self.element_expr(&MPoly::from_terms(m, ts), den)))
            .collect();
        self.convert_poly(&ExpressionPoly::from_terms(terms, gens)?)
    }

    // ------------------------------------------------------------------
    // ground elements (polynomials in the symbols)
    // ------------------------------------------------------------------

    /// The ground domain gcd `ZZ_I.gcd`/`QQ_I.gcd`.
    fn scalar_gcd(&self, a: &AN, b: &AN) -> AN {
        if self.algebraic.is_some() {
            // `Field.gcd` without an associated ring
            AN::one()
        } else if self.integral {
            gi::gcd(a, b)
        } else {
            gi::field_gcd(a, b)
        }
    }
    /// `PolyElement.gcd` in `ZZ_I[y...]`/`QQ_I[y...]`. Gaussian elements are
    /// never nonnegative, so `gcd(0, g) = -g`.
    fn egcd(&self, a: &P, b: &P) -> PResult<P> {
        // `_gcd_zero`: `g` when `g.is_nonnegative`, which over `QQ<a>` is
        // the sign of the leading coordinate of the leading coefficient.
        let signed = |g: &P| if self.algebraic.is_some() && !g.is_zero() && !g.lc().is_negative() { g.clone() } else { g.neg() };
        if a.is_zero() {
            return Ok(signed(b));
        }
        if b.is_zero() {
            return Ok(signed(a));
        }
        if a.len() == 1 || b.len() == 1 {
            let (mono, single, other) = if a.len() == 1 { (a.lm().to_vec(), a, b) } else { (b.lm().to_vec(), b, a) };
            let mut mgcd = mono;
            let mut cgcd = single.lc();
            for (e, c) in other.iter() {
                for (g, &x) in mgcd.iter_mut().zip(e) {
                    *g = (*g).min(x);
                }
                cgcd = self.scalar_gcd(&cgcd, c);
            }
            return Ok(MPoly::monomial(&mgcd, cgcd));
        }
        if self.integral {
            Ok(gi::cofactors(a, b)?.0)
        } else {
            field_poly_gcd(a, b)
        }
    }
    /// `PolyElement.lcm`.
    fn elcm(&self, a: &P, b: &P) -> PResult<P> {
        if a.is_zero() || b.is_zero() {
            return Ok(MPoly::zero(a.nvars));
        }
        if self.integral {
            let (fc, f) = gi::ground_primitive(a, true);
            let (gc, g) = gi::ground_primitive(b, true);
            let c = gi::lcm(&fc, &gc);
            let h = exact(&f.mul(&g), &self.egcd(&f, &g)?)?;
            Ok(h.scale(&c))
        } else {
            Ok(exact(&a.mul(b), &self.egcd(a, b)?)?.monic())
        }
    }
    fn is_ground_one(p: &P) -> bool {
        p.is_one()
    }

    // ------------------------------------------------------------------
    // ring grounds: contents and gcds over the flattened polynomial
    // ------------------------------------------------------------------

    /// `dmp_ground_content` over the ring ground of a polynomial in `k`
    /// generators followed by the symbols, folded over the dense form.
    fn ring_content(&self, f: &P, k: usize) -> PResult<P> {
        if k == 0 {
            return Ok(f.clone());
        }
        let m = self.m();
        if f.is_zero() {
            return Ok(MPoly::zero(m));
        }
        let mut cont = MPoly::zero(m);
        for c in dense(f) {
            let inner = self.ring_content(&drop_vars(&c, 1), k - 1)?;
            cont = self.egcd(&cont, &inner)?;
            if Self::is_ground_one(&cont) {
                break;
            }
        }
        Ok(cont)
    }
    fn ring_primitive(&self, f: &P, k: usize) -> PResult<(P, P)> {
        if f.is_zero() {
            return Ok((MPoly::zero(self.m()), f.clone()));
        }
        let c = self.ring_content(f, k)?;
        if c.is_one() {
            return Ok((c, f.clone()));
        }
        let q = exact(f, &lift(&c, k))?;
        Ok((c, q))
    }
    /// A gcd of two flattened polynomials over `ZZ_I`, defined up to a
    /// ground factor.
    fn flat_gcd(&self, f: &P, g: &P) -> PResult<P> {
        if self.algebraic.is_some() {
            field_poly_gcd(f, g)
        } else if self.integral {
            Ok(gi::cofactors(f, g)?.0)
        } else {
            let (df, f) = gi::clear_denoms(f);
            let (dg, g) = gi::clear_denoms(g);
            // Gaussian-rational content survives the passage to ZZ_I:
            // gcd(df*f, dg*g)/lcm(df, dg). Without this factor a rational
            // common content is silently multiplied by its denominator.
            let scale = AN::rational(Rat::new(1.into(), df.lcm(&dg)), gi::field());
            Ok(gi::cofactors(&f, &g)?.0.scale(&scale))
        }
    }
    /// `dup_rr_prs_gcd`/`dmp_rr_prs_gcd` over the ring ground.
    fn ring_gcd(&self, f: &P, g: &P, k: usize) -> PResult<P> {
        if self.algebraic.is_some() {
            // `dmp_inner_gcd` over `QQ<a>[y...]`, which is neither `ZZ` nor
            // a field: the subresultant gcd with its own content scaling.
            return Ok(super::rr_prs_gcd::ground_inner_gcd(f, g, k)?.0);
        }
        if f.is_zero() {
            return Ok(g.neg());
        }
        if g.is_zero() {
            return Ok(f.neg());
        }
        if k > 1 {
            return self.flat_gcd(f, g);
        }
        let (fc, ff) = self.ring_primitive(f, k)?;
        let (gc, gg) = self.ring_primitive(g, k)?;
        let c = self.egcd(&fc, &gc)?;
        let h = self.flat_gcd(&ff, &gg)?;
        let (_, h) = self.ring_primitive(&h, k)?;
        let unit = gi::canonical_unit(&lc_ground(&h, k).lc());
        Ok(h.mul(&lift(&c.scale(&unit), k)))
    }
    /// `dup_rr_lcm`/`dmp_rr_lcm`.
    fn ring_lcm(&self, f: &P, g: &P, k: usize) -> PResult<P> {
        if f.is_zero() || g.is_zero() {
            return Ok(MPoly::zero(f.nvars));
        }
        let (fc, ff) = self.ring_primitive(f, k)?;
        let (gc, gg) = self.ring_primitive(g, k)?;
        let c = self.elcm(&fc, &gc)?;
        let h = exact(&ff.mul(&gg), &self.ring_gcd(&ff, &gg, k)?)?;
        let c = if k == 1 { c.scale(&self.unit(&lc_ground(&h, k).lc())) } else { c };
        Ok(h.mul(&lift(&c, k)))
    }

    // ------------------------------------------------------------------
    // Poly-level operations
    // ------------------------------------------------------------------

    fn marked(&self, p: &ExpressionPoly) -> ExpressionPoly {
        let mut p = p.clone();
        p.ground = Some(self.clone());
        p
    }
    /// `dmp_ground_monic` through `EX`, renormalized in the fraction field.
    pub fn monic(&self, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::monic(&self.marked(p));
        }
        let field = self.to_field();
        if self.algebraic.is_some() {
            return field.algebraic_monic(p);
        }
        field.convert_poly(&p.monic()?)
    }
    /// `dmp_ground_monic` over `QQ<a>(y...)` (after `to_field`): every
    /// coefficient is the `FracElement` quotient `c/lc`, that is
    /// `cancel(c.numer*lc.denom, c.denom*lc.numer)`, kept as that pair
    /// rather than renormalized from its expression.
    fn algebraic_monic(&self, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        if p.is_zero() {
            return Ok(self.marked(p));
        }
        let (ln, ld) = self.element(&p.lc())?;
        if ln.is_one() && ld.is_one() {
            return self.convert_poly(p);
        }
        let mut terms = Vec::new();
        for (m, c) in p.terms() {
            let (cn, cd) = self.element(&c)?;
            let (a, b) = self.cancel_pair(&cn.mul(&ld), &cd.mul(&ln))?;
            terms.push((m, self.element_expr(&a, &b)));
        }
        let mut out = ExpressionPoly::from_terms(terms, &p.gens)?;
        out.ground = Some(self.clone());
        Ok(out)
    }
    pub fn gcd(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::gcd(&self.marked(f), &self.marked(g));
        }
        let k = f.gens.len();
        let (nf, _) = self.flatten(f)?;
        let (ng, _) = self.flatten(g)?;
        let _s = self.scope();
        let one = MPoly::one(self.m());
        if self.ring {
            return self.eject(&self.ring_gcd(&nf, &ng, k)?, &one, &f.gens);
        }
        if nf.is_zero() && ng.is_zero() {
            return Ok(self.marked(f));
        }
        let h = if nf.is_zero() {
            ng
        } else if ng.is_zero() {
            nf
        } else {
            self.flat_gcd(&nf, &ng)?
        };
        self.monic(&self.eject(&h, &one, &f.gens)?)
    }
    pub fn lcm(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::lcm(&self.marked(f), &self.marked(g));
        }
        let k = f.gens.len();
        let (nf, _) = self.flatten(f)?;
        let (ng, _) = self.flatten(g)?;
        let _s = self.scope();
        let one = MPoly::one(self.m());
        if self.ring {
            return self.eject(&self.ring_lcm(&nf, &ng, k)?, &one, &f.gens);
        }
        if nf.is_zero() || ng.is_zero() {
            return self.convert_poly(&ExpressionPoly::zero(&f.gens));
        }
        let h = exact(&nf.mul(&ng), &self.flat_gcd(&nf, &ng)?)?;
        self.monic(&self.eject(&h, &one, &f.gens)?)
    }
    pub fn cofactors(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(ExpressionPoly, ExpressionPoly, ExpressionPoly)> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::cofactors(&self.marked(f), &self.marked(g));
        }
        let h = self.gcd(f, g)?;
        if h.is_zero() {
            return Ok((h, self.marked(f), self.marked(g)));
        }
        let a = self.convert_poly(&f.exquo(&h)?)?;
        let b = self.convert_poly(&g.exquo(&h)?)?;
        Ok((h, a, b))
    }
    /// `dmp_ground_content` over the ground domain, as an expression.
    pub fn content(&self, p: &ExpressionPoly) -> PResult<Expr> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::content(&self.marked(p));
        }
        let k = p.gens.len();
        let m = self.m();
        if self.ring {
            let (n, _) = self.flatten(p)?;
            let _s = self.scope();
            return Ok(self.element_expr(&self.ring_content(&n, k)?, &MPoly::one(m)));
        }
        let _s = self.scope();
        let (c, d) = self.field_content(p)?;
        Ok(self.element_expr(&c, &d))
    }
    /// `Field.gcd` folded over the dense coefficients.
    fn field_content(&self, p: &ExpressionPoly) -> PResult<(P, P)> {
        let m = self.m();
        let ring = self.to_ring();
        let mut cont = (MPoly::zero(m), MPoly::one(m));
        if p.is_zero() {
            return Ok(cont);
        }
        let terms: BTreeMap<Vec<u32>, Expr> = p.terms().into_iter().collect();
        for mono in dense_monomials(p) {
            let (a, b) = match terms.get(&mono) {
                Some(c) => self.element(c)?,
                None => (MPoly::zero(m), MPoly::one(m)),
            };
            let num = ring.egcd(&cont.0, &a)?;
            let den = ring.elcm(&cont.1, &b)?;
            cont = if self.algebraic.is_some() {
                // `self.convert(p, ring)/q`: `FracElement.new(p, q)`
                self.cancel_pair(&num, &den)?
            } else {
                let e = self.element_expr(&num, &den);
                self.element(&e)?
            };
            if cont.0.is_one() && cont.1.is_one() {
                break;
            }
        }
        Ok(cont)
    }
    pub fn primitive(&self, p: &ExpressionPoly) -> PResult<(Expr, ExpressionPoly)> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::primitive(&self.marked(p));
        }
        let c = self.content(p)?;
        if p.is_zero() || c.is_one() {
            return Ok((c, self.marked(p)));
        }
        let q = p.exquo_ground(&c)?;
        Ok((c, self.convert_poly(&q)?))
    }
    /// `Poly.clear_denoms()`.
    pub fn clear_denoms(&self, p: &ExpressionPoly) -> PResult<(Expr, ExpressionPoly)> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::clear_denoms(&self.marked(p));
        }
        if self.ring {
            return Ok((Expr::one(), self.marked(p)));
        }
        let m = self.m();
        let ring = self.to_ring();
        let _s = self.scope();
        let mut common = MPoly::one(m);
        for (_, c) in p.terms() {
            let (_, d) = self.element(&c)?;
            common = ring.elcm(&common, &d)?;
        }
        let ce = self.element_expr(&common, &MPoly::one(m));
        Ok((ce.clone(), self.convert_poly(&p.mul_ground(&ce)?)?))
    }
    /// `dmp_factor_list` over the composite ground.
    pub fn factor_list(&self, p: &ExpressionPoly) -> PResult<(Expr, Vec<(ExpressionPoly, usize)>)> {
        let _algebraic_scope = self.scope();
        if self.is_scalar() {
            // `dup_factor_list` raises DomainError over `CC` (its exact
            // domain `QQ_I` has the ring `ZZ_I`, which is neither `ZZ` nor
            // a polynomial ring) and `Poly.factor_list` returns the input.
            if p.is_zero() || p.terms().iter().all(|(m, _)| m.iter().all(|&e| e == 0)) {
                return Ok((p.as_expr(), Vec::new()));
            }
            return Ok((Expr::one(), vec![(p.clone(), 1)]));
        }
        if self.inexact.is_some() {
            return Err(super::inexact_ground::factor_error(self));
        }
        let k = p.gens.len();
        let m = self.m();
        let (num, den) = self.flatten(p)?;
        let _s = self.scope();
        let one = MPoly::one(m);
        if num.is_zero() {
            return Ok((Expr::zero(), Vec::new()));
        }
        let mut all_gens = p.gens.clone();
        all_gens.extend(self.symbols.iter().cloned());
        let gaussian = if self.algebraic.is_some() { None } else { Some(self.integral) };
        let flat = AlgebraicPoly { rep: num.clone(), gens: all_gens, number_field: self.field(), gaussian };
        let (_, factors) = flat.factor_list()?;
        let mut out: Vec<(P, u32)> = Vec::new();
        for (f, i) in factors {
            if f.rep.iter().any(|(e, _)| e[..k].iter().any(|&x| x > 0)) {
                out.push((f.rep, i));
            }
        }
        // the unique coefficient: LC(p) / prod LC(f_i)**k_i
        let mut lead = one.clone();
        for (f, i) in &out {
            lead = lead.mul(&lc_ground(f, k).pow(*i));
        }
        let lc_num = lc_ground(&num, k);
        let coeff = if self.ring {
            let q = exact(&lc_num, &lead)?;
            self.element_expr(&q, &one)
        } else {
            let e = self.element_expr(&lc_num, &den.mul(&lead));
            self.reduce(&e)?
        };
        sort_nested(&mut out, k, self.algebraic.is_some());
        let factors = out
            .into_iter()
            .map(|(f, i)| Ok((self.eject(&f, &one, &p.gens)?, i as usize)))
            .collect::<PResult<Vec<_>>>()?;
        Ok((coeff, factors))
    }
    /// `dup_sqf_list` over the composite ground.
    pub fn sqf_list(&self, p: &ExpressionPoly, all: bool) -> PResult<(Expr, Vec<(ExpressionPoly, usize)>)> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::sqf_list(&self.marked(p), all);
        }
        if p.gens.len() != 1 {
            return self.sqf_list_multi(p, all);
        }
        let (coeff, f) = if self.ring {
            let (c, f) = self.primitive(p)?;
            (c, f)
        } else {
            (p.lc(), self.monic(p)?)
        };
        if f.degree(0) <= 0 {
            return Ok((coeff, Vec::new()));
        }
        let mut result = Vec::new();
        let mut i = 1;
        let h = self.convert_poly(&f.diff(1, 0)?)?;
        let (_, mut pp, mut q) = self.cofactors(&f, &h)?;
        loop {
            let d = self.convert_poly(&pp.diff(1, 0)?)?;
            let h = self.convert_poly(&q.sub(&d)?)?;
            if h.is_zero() {
                result.push((pp, i));
                break;
            }
            let (g, p2, q2) = self.cofactors(&pp, &h)?;
            if all || g.degree(0) > 0 {
                result.push((g, i));
            }
            pp = p2;
            q = q2;
            i += 1;
        }
        Ok((coeff, result))
    }
    /// `dmp_sqf_list` in several generators: the square-free part of
    /// multiplicity `i` is the product of the irreducible factors of
    /// multiplicity `i`, computed from the factorization in all variables.
    /// The reference runs `dmp_inner_gcd` over the composite ground, whose
    /// cofactors sometimes lose a unit `-1` (its product is then `-f`); the
    /// product kept here is always `f`. Over a field the parts are monic.
    fn sqf_list_multi(&self, p: &ExpressionPoly, all: bool) -> PResult<(Expr, Vec<(ExpressionPoly, usize)>)> {
        let (mut coeff, factors) = self.factor_list(p)?;
        let top = factors.iter().map(|(_, k)| *k).max().unwrap_or(0);
        let one = self.convert_poly(&ExpressionPoly::from_expr(&Expr::one(), &p.gens)?)?;
        let mut result = Vec::new();
        for i in 1..=top {
            let mut g: Option<ExpressionPoly> = None;
            for (f, k) in &factors {
                if *k == i {
                    g = Some(match g {
                        None => f.clone(),
                        Some(h) => self.convert_poly(&h.mul(f)?)?,
                    });
                }
            }
            match g {
                Some(mut g) => {
                    if !self.ring {
                        let lc = g.lc();
                        g = self.monic(&g)?;
                        coeff = crate::ops::mul2(&coeff, &crate::ops::powi(&lc, i as i64));
                    }
                    result.push((g, i));
                }
                None if all => result.push((one.clone(), i)),
                None => {}
            }
        }
        if !self.ring {
            coeff = self.reduce(&coeff)?;
        }
        Ok((coeff, result))
    }
    /// `dmp_cancel(include=False)`: `(cp, cq, p, q)` with `f/g = cp/cq * p/q`.
    pub fn cancel(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(Expr, Expr, ExpressionPoly, ExpressionPoly)> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::cancel(&self.marked(f), &self.marked(g));
        }
        if self.algebraic.is_some() {
            return self.algebraic_cancel(f, g);
        }
        if self.ring {
            let (_, p, q) = self.cofactors(f, g)?;
            return Ok((Expr::one(), Expr::one(), p, q));
        }
        let ring = self.to_ring();
        let (cq, f1) = self.clear_denoms(f)?;
        let (cp, g1) = self.clear_denoms(g)?;
        let (_, p, q) = ring.cofactors(&ring.convert_poly(&f1)?, &ring.convert_poly(&g1)?)?;
        let _s = self.scope();
        let (a, _) = ring.element(&cp)?;
        let (b, _) = ring.element(&cq)?;
        let h = ring.egcd(&a, &b)?;
        let one = MPoly::one(self.m());
        let (cp, cq) = if h.is_zero() {
            (a, b)
        } else {
            (exact(&a, &h)?, exact(&b, &h)?)
        };
        Ok((self.element_expr(&cp, &one), self.element_expr(&cq, &one), self.convert_poly(&p)?, self.convert_poly(&q)?))
    }
    /// `dmp_cancel(include=False)` over `QQ<a>[y...]` or `QQ<a>(y...)`,
    /// with the reference sign normalization of both sides.
    fn algebraic_cancel(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(Expr, Expr, ExpressionPoly, ExpressionPoly)> {
        let m = self.m();
        let _scope = self.scope();
        let one = MPoly::one(m);
        let (cq, f1, cp, g1) = if self.ring {
            (one.clone(), f.clone(), one.clone(), g.clone())
        } else {
            let ring = self.to_ring();
            let (cq, f1) = self.clear_denoms(f)?;
            let (cp, g1) = self.clear_denoms(g)?;
            let _s = self.scope();
            (ring.element(&cq)?.0, ring.convert_poly(&f1)?, ring.element(&cp)?.0, ring.convert_poly(&g1)?)
        };
        let ring = self.to_ring();
        let k = f.gens.len();
        let (nf, _) = ring.flatten(&f1)?;
        let (ng, _) = ring.flatten(&g1)?;
        let _s = self.scope();
        let (_, mut p, mut q) = super::rr_prs_gcd::ground_inner_gcd(&nf, &ng, k)?;
        let (mut cp, cq) = if self.ring {
            (cp, cq)
        } else {
            let (_, a, b) = algebraic_cofactors(&cp, &cq)?;
            (a, b)
        };
        // K.is_negative(dmp_ground_LC(.)): for a fraction field the sign of
        // the numerator's leading coefficient.
        let negative = |h: &P| !h.is_zero() && h.lc().is_negative();
        let (p_neg, q_neg) = (negative(&p), negative(&q));
        if p_neg && q_neg {
            p = p.neg();
            q = q.neg();
        } else if p_neg {
            cp = cp.neg();
            p = p.neg();
        } else if q_neg {
            cp = cp.neg();
            q = q.neg();
        }
        let p = self.eject(&p, &one, &f.gens)?;
        let q = self.eject(&q, &one, &f.gens)?;
        Ok((self.element_expr(&cp, &one), self.element_expr(&cq, &one), p, q))
    }
    /// `dup_sqf_part`/`dmp_sqf_part`.
    pub fn sqf_part(&self, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let _algebraic_scope = self.scope();
        if self.inexact.is_some() {
            return super::inexact_ground::sqf_part(&self.marked(p));
        }
        if p.is_zero() {
            return Ok(self.marked(p));
        }
        let mut g = p.clone();
        g.ground = Some(self.clone());
        for j in 0..p.gens.len() {
            let d = self.convert_poly(&p.diff(1, j)?)?;
            g = self.gcd(&g, &d)?;
        }
        let s = self.convert_poly(&p.exquo(&g)?)?;
        if self.ring {
            Ok(self.primitive(&s)?.1)
        } else {
            self.monic(&s)
        }
    }
}

/// The monomials of the dense form in the first generator of `p`, leading
/// first (the other generators at their dense positions too).
fn dense_monomials(p: &ExpressionPoly) -> Vec<Vec<u32>> {
    let degrees: Vec<i64> = (0..p.gens.len()).map(|j| p.degree(j).max(0)).collect();
    let mut out = vec![Vec::new()];
    for d in degrees {
        out = out
            .into_iter()
            .flat_map(|prefix: Vec<u32>| {
                (0..=d as u32).rev().map(move |k| {
                    let mut m = prefix.clone();
                    m.push(k);
                    m
                })
            })
            .collect();
    }
    out
}

/// `PolyElement.sort_key` of a ground element: the number of terms, then
/// the terms in decreasing monomial order with the Gaussian order `[y, x]`.
fn element_sort_key(p: &P, algebraic: bool) -> (usize, Vec<(Vec<u32>, Vec<Rat>)>) {
    // an `ANP` compares by its coordinate list, leading first
    let key = |c: &AN| if algebraic { c.coefficients.iter().rev().cloned().collect() } else { gi::element_key(c) };
    (p.len(), p.iter().map(|(e, c)| (e.to_vec(), key(c))).collect())
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum NestedKey {
    Leaf((usize, Vec<(Vec<u32>, Vec<Rat>)>)),
    Node(Vec<NestedKey>),
}

fn nested_key(p: &P, k: usize, algebraic: bool) -> NestedKey {
    if k == 0 {
        return NestedKey::Leaf(element_sort_key(p, algebraic));
    }
    NestedKey::Node(dense(p).iter().map(|c| nested_key(&drop_vars(c, 1), k - 1, algebraic)).collect())
}

/// The reference `_sort_factors(multiple=True)` over a composite ground.
fn sort_nested(factors: &mut [(P, u32)], k: usize, algebraic: bool) {
    factors.sort_by_cached_key(|(p, i)| (dense(p).len(), *i, nested_key(p, k, algebraic)));
}

/// Normalize a marked `EX` polynomial in its ground, falling back to the
/// fraction field and then to plain `EX` (the reference `to_field`).
pub fn normalize_poly(p: ExpressionPoly) -> ExpressionPoly {
    let Some(g) = p.ground.clone() else { return p };
    if let Ok(q) = g.convert_poly(&p) {
        return q;
    }
    if g.ring {
        if let Ok(q) = g.to_field().convert_poly(&p) {
            return q;
        }
    }
    ExpressionPoly { ground: None, ..p }
}

/// The unification descriptor of a polynomial domain: a Gaussian or
/// rational ground, integral or not, as a polynomial ring or fraction field
/// over `symbols` (empty for a bare ground).
#[derive(Clone, Debug)]
pub struct Descriptor {
    pub gaussian: bool,
    pub integral: bool,
    pub ring: bool,
    pub symbols: Vec<Expr>,
    /// An inexact ground (`RR`/`CC` at some precision).
    pub inexact: Option<super::inexact_ground::Inexact>,
    /// An algebraic field ground `QQ<a>` (`gaussian` is then false).
    pub algebraic: Option<Arc<NumberField>>,
}

pub fn descriptor(p: &super::fraction::Polynomial) -> Option<Descriptor> {
    if let Some(e) = &p.expression {
        return e.ground.as_ref().map(|g| Descriptor {
            gaussian: g.inexact.is_none() && g.algebraic.is_none(),
            integral: g.integral,
            ring: g.ring,
            symbols: g.symbols.clone(),
            inexact: g.inexact,
            algebraic: g.algebraic.clone(),
        });
    }
    if let Some(a) = &p.algebraic {
        return Some(match a.gaussian {
            Some(integral) => Descriptor { gaussian: true, integral, ring: true, symbols: Vec::new(), inexact: None, algebraic: None },
            None => Descriptor { gaussian: false, integral: false, ring: true, symbols: Vec::new(), inexact: None, algebraic: Some(a.number_field.clone()) },
        });
    }
    if p.modulus.is_some() {
        return None;
    }
    if let Some(f) = &p.field {
        return match &f.coefficient_domain {
            ExactDomain::Functions { integral, ring, parameters, .. } => {
                Some(Descriptor { gaussian: false, integral: *integral, ring: *ring, symbols: parameters.clone(), inexact: None, algebraic: None })
            }
            ExactDomain::ZZ => Some(Descriptor { gaussian: false, integral: true, ring: true, symbols: Vec::new(), inexact: None, algebraic: None }),
            ExactDomain::QQ => Some(Descriptor { gaussian: false, integral: false, ring: true, symbols: Vec::new(), inexact: None, algebraic: None }),
            ExactDomain::RR => Some(Descriptor {
                gaussian: false,
                integral: false,
                ring: true,
                symbols: Vec::new(),
                inexact: Some(super::inexact_ground::Inexact { precision: 53, complex: false }),
                algebraic: None,
            }),
            ExactDomain::Inexact { precision, complex } => Some(Descriptor {
                gaussian: false,
                integral: false,
                ring: true,
                symbols: Vec::new(),
                inexact: Some(super::inexact_ground::Inexact { precision: *precision, complex: *complex }),
                algebraic: None,
            }),
            _ => None,
        };
    }
    let integral = match p.base.domain {
        super::polyobj::Domain::ZZ => true,
        super::polyobj::Domain::QQ => false,
        super::polyobj::Domain::RR(precision) => {
            let symbols = p.base.domain_gens().to_vec();
            if symbols == [Expr::i()] {
                // the legacy float representation of `CC`
                return Some(Descriptor {
                    gaussian: false,
                    integral: false,
                    ring: true,
                    symbols: Vec::new(),
                    inexact: Some(super::inexact_ground::Inexact { precision, complex: true }),
                    algebraic: None,
                });
            }
            if symbols.iter().any(|s| *s == Expr::i() || (s.is_number() && s.is_pow())) {
                return None;
            }
            return Some(Descriptor {
                gaussian: false,
                integral: false,
                ring: true,
                symbols,
                inexact: Some(super::inexact_ground::Inexact { precision, complex: false }),
                algebraic: None,
            });
        }
    };
    let symbols = p.base.domain_gens().to_vec();
    // radicals select `EX` (or an algebraic field) in the reference
    if symbols.iter().any(|s| *s == Expr::i() || (s.is_number() && s.is_pow())) {
        return None;
    }
    Some(Descriptor { gaussian: false, integral, ring: true, symbols, inexact: None, algebraic: None })
}

fn descriptor_display(d: &Descriptor) -> String {
    let name = d.algebraic.as_ref().map(|k| k.display());
    let ground = match (d.gaussian, d.integral) {
        _ if name.is_some() => name.as_deref().unwrap_or_default(),
        _ if d.inexact.is_some_and(|k| k.complex) => "CC",
        _ if d.inexact.is_some() => "RR",
        (true, true) => "ZZ_I",
        (true, false) => "QQ_I",
        (false, true) => "ZZ",
        (false, false) => "QQ",
    };
    if d.symbols.is_empty() {
        return ground.into();
    }
    let s: Vec<String> = d.symbols.iter().map(|e| e.to_string()).collect();
    if d.ring {
        format!("{ground}[{}]", s.join(","))
    } else {
        format!("{ground}({})", s.join(","))
    }
}

/// The reference `Domain.unify_with_symbols` when one side has a Gaussian
/// ground: `None` when neither does.
pub fn unify(a: &Descriptor, b: &Descriptor, gens: &[Expr]) -> PResult<Option<GaussianGround>> {
    let composite_algebraic = |d: &Descriptor| d.algebraic.is_some() && !d.symbols.is_empty();
    if composite_algebraic(a) || composite_algebraic(b) {
        return unify_algebraic(a, b, gens);
    }
    if a.algebraic.is_some() || b.algebraic.is_some() {
        // a bare `QQ<a>` keeps its own unification paths
        return Ok(None);
    }
    if !a.gaussian && !b.gaussian && a.inexact.is_none() && b.inexact.is_none() {
        return Ok(None);
    }
    for d in [a, b] {
        if d.symbols.iter().any(|s| gens.contains(s)) {
            let gs: Vec<String> = gens.iter().map(|g| g.to_string()).collect();
            let tuple = if gs.len() == 1 { format!("({},)", gs[0]) } else { format!("({})", gs.join(", ")) };
            return Err(PolyError::Reference {
                kind: "UnificationFailed",
                message: format!("Cannot unify {} with {}, given {tuple} generators", descriptor_display(a), descriptor_display(b)),
            });
        }
    }
    let mut integral = a.integral && b.integral;
    let symbols = super::polyobj::unify_gens(&a.symbols, &b.symbols);
    let composite = |d: &Descriptor| !d.symbols.is_empty();
    let frac = (composite(a) && !a.ring) || (composite(b) && !b.ring);
    let poly_frac = (composite(a) && !a.ring && composite(b) && b.ring) || (composite(b) && !b.ring && composite(a) && a.ring);
    if poly_frac && (a.integral || b.integral) && !integral {
        integral = true;
    }
    if symbols.is_empty() {
        // `CC` with anything, or `RR` with `ZZ_I`/`QQ_I`, unify to the plain
        // field `CC` at the largest precision.
        let inexact = a.inexact.is_some() || b.inexact.is_some();
        let complex = [a, b].iter().any(|d| d.inexact.is_some_and(|k| k.complex)) || (inexact && (a.gaussian || b.gaussian));
        if complex {
            let precision = [a, b].iter().filter_map(|d| d.inexact.map(|k| k.precision)).max().unwrap_or(53);
            return Ok(Some(GaussianGround::complex_field(precision)));
        }
        return Ok(None);
    }
    if a.inexact.is_some() || b.inexact.is_some() {
        let precision = [a, b].iter().filter_map(|d| d.inexact.map(|k| k.precision)).max().unwrap_or(53);
        let complex = [a, b].iter().any(|d| d.gaussian || d.inexact.is_some_and(|k| k.complex));
        let inexact = Some(super::inexact_ground::Inexact { precision, complex });
        return Ok(Some(GaussianGround { integral: false, ring: !frac, symbols, inexact, algebraic: None }));
    }
    Ok(Some(GaussianGround { integral, ring: !frac, symbols, inexact: None, algebraic: None }))
}

fn generator_clash(a: &Descriptor, b: &Descriptor, gens: &[Expr]) -> PResult<()> {
    for d in [a, b] {
        if d.symbols.iter().any(|s| gens.contains(s)) {
            let gs: Vec<String> = gens.iter().map(|g| g.to_string()).collect();
            let tuple = if gs.len() == 1 { format!("({},)", gs[0]) } else { format!("({})", gs.join(", ")) };
            return Err(PolyError::Reference {
                kind: "UnificationFailed",
                message: format!("Cannot unify {} with {}, given {tuple} generators", descriptor_display(a), descriptor_display(b)),
            });
        }
    }
    Ok(())
}

/// `unify_with_symbols` when one side is `QQ<a>[y...]` or `QQ<a>(y...)`
/// and the other any exact or inexact ground or a composite over one: the
/// symbols are merged, a fraction field on either side gives a fraction
/// field, and the ground is the reference unification of both grounds.
fn unify_algebraic(a: &Descriptor, b: &Descriptor, gens: &[Expr]) -> PResult<Option<GaussianGround>> {
    generator_clash(a, b, gens)?;
    let symbols = super::polyobj::unify_gens(&a.symbols, &b.symbols);
    let composite = |d: &Descriptor| !d.symbols.is_empty();
    let frac = (composite(a) && !a.ring) || (composite(b) && !b.ring);
    if a.inexact.is_some() || b.inexact.is_some() {
        // `RealField.unify`/`ComplexField.unify` keep the inexact field: an
        // algebraic ground contributes nothing, `ZZ_I`/`QQ_I` a complex one.
        let precision = [a, b].iter().filter_map(|d| d.inexact.map(|k| k.precision)).max().unwrap_or(53);
        let complex = [a, b].iter().any(|d| d.gaussian || d.inexact.is_some_and(|k| k.complex));
        let inexact = Some(super::inexact_ground::Inexact { precision, complex });
        return Ok(Some(GaussianGround { integral: false, ring: !frac, symbols, inexact, algebraic: None }));
    }
    // `AlgebraicField.unify` after the reference swap that puts an
    // algebraic second operand first: the field generated by both original
    // extensions (`QQ_I` counts as `QQ<I>`), through its primitive element.
    let i = Expr::i();
    let field = match (&a.algebraic, &b.algebraic) {
        (Some(x), Some(y)) if x.same_embedding(y) => x.clone(),
        (Some(x), Some(y)) => NumberField::from_generators(&[y.generator.clone(), x.generator.clone()])?,
        (Some(x), None) if b.gaussian => NumberField::from_generators(&[x.generator.clone(), i])?,
        (None, Some(y)) if a.gaussian => NumberField::from_generators(&[y.generator.clone(), i])?,
        (Some(x), None) | (None, Some(x)) => x.clone(),
        (None, None) => return Ok(None),
    };
    Ok(Some(GaussianGround { integral: false, ring: !frac, symbols, inexact: None, algebraic: Some(field) }))
}

/// Whether `e` is an algebraic number in the reference sense of
/// `is_number and is_algebraic`, for the expressions built from numbers,
/// `I` and rational powers; anything else counts as transcendental.
pub(crate) fn is_algebraic_number(e: &Expr) -> bool {
    use crate::expr::Kind;
    match e.kind() {
        Kind::Num(n) => !n.is_float(),
        Kind::Const(_) => *e == Expr::i(),
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().all(is_algebraic_number),
        Kind::Pow(b, x) => is_algebraic_number(b) && x.as_number().is_some_and(|n| n.is_rational()),
        _ => false,
    }
}

/// The reference `_construct_composite` for the coefficients of `e` in
/// `gens` when they involve `I`: `ZZ_I`/`QQ_I` as a polynomial ring or
/// (with a nonconstant denominator, or `field`) a fraction field over the
/// other coefficient generators. `None` when the reference selects another
/// domain (no `I`, no other generator, floats, number-like algebraic
/// generators, or generators sharing symbols).
pub fn infer(e: &Expr, gens: &[Expr], field: bool) -> Option<GaussianGround> {
    let coefficients: Vec<Expr> = if gens.is_empty() {
        vec![e.clone()]
    } else {
        super::expression_domain::dict_from_expr(e, gens).ok()?.into_iter().map(|(_, c)| c).collect()
    };
    infer_coefficients(&coefficients, field)
}

pub fn infer_coefficients(coefficients: &[Expr], field: bool) -> Option<GaussianGround> {
    let i = Expr::i();
    let mut parts = Vec::with_capacity(2 * coefficients.len());
    let (mut numers, mut denoms) = (Vec::new(), Vec::new());
    for c in coefficients {
        let (n, d) = super::api::as_numer_denom(c);
        numers.push(n);
        denoms.push(d);
    }
    parts.extend(numers.iter().cloned());
    parts.extend(denoms.iter().cloned());
    let converted = super::convert::from_exprs(&parts).ok()?;
    if converted.float || !converted.gens.contains(&i) {
        return None;
    }
    let symbols: Vec<Expr> = converted.gens.iter().filter(|g| **g != i).cloned().collect();
    if symbols.is_empty() || symbols.iter().any(|g| crate::traverse::free_symbols(g).is_empty() && is_algebraic_number(g)) {
        return None;
    }
    let mut seen = Vec::new();
    for g in &symbols {
        for s in crate::traverse::free_symbols(g) {
            if seen.contains(&s) {
                return None;
            }
            seen.push(s);
        }
    }
    let _s = scope();
    let polys: Vec<P> = parts.iter().map(|p| gaussian_poly(p, &symbols)).collect::<PResult<_>>().ok()?;
    let k = coefficients.len();
    let fractions = field || polys[k..].iter().any(|d| !d.is_constant());
    let integral = if fractions {
        polys.iter().all(gi::poly_is_integral)
    } else {
        polys[..k].iter().zip(&polys[k..]).all(|(n, d)| gi::poly_is_integral(&n.scale(&d.lc().inv())))
    };
    Some(GaussianGround { integral, ring: !fractions, symbols, inexact: None, algebraic: None })
}
