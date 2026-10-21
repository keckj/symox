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
use super::algebraic::{AlgebraicPoly, AlgebraicScope, AN};
use super::convert::{PResult, PolyError};
use super::domains::ExactDomain;
use super::expression_domain::ExpressionPoly;
use super::gaussian as gi;
use super::int::{Field, Rat};
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
fn class_path(e: &Expr) -> &'static str {
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
        ExactDomain::Composite {
            ground: Arc::new(ExactDomain::Gaussian { integral: self.integral }),
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
    pub fn from_domain(d: &ExactDomain) -> Option<Self> {
        let ExactDomain::Composite { ground, ring, parameters, .. } = d else { return None };
        let ExactDomain::Gaussian { integral } = ground.as_ref() else { return None };
        if parameters.is_empty() {
            return None;
        }
        Some(Self { integral: *integral, ring: *ring, symbols: parameters.clone() })
    }
    fn m(&self) -> usize {
        self.symbols.len()
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
        let _s = scope();
        let m = self.m();
        let (n, d) = super::api::as_numer_denom(e);
        let num = gaussian_poly(&n, &self.symbols).map_err(|_| self.coercion(e))?;
        let den = gaussian_poly(&d, &self.symbols).map_err(|_| self.coercion(e))?;
        if den.is_zero() {
            return Err(self.coercion(e));
        }
        if num.is_zero() {
            return Ok((MPoly::zero(m), MPoly::one(m)));
        }
        if self.ring {
            let (p, q) = if den.is_constant() { (num, den) } else { gi::expr_cancel(&num, &den)? };
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
        AlgebraicPoly { rep: p.clone(), gens: self.symbols.clone(), number_field: gi::field(), gaussian: Some(false) }.as_expr()
    }
    pub fn element_expr(&self, p: &P, q: &P) -> Expr {
        if q.is_one() {
            self.poly_expr(p)
        } else {
            crate::ops::div(&self.poly_expr(p), &self.poly_expr(q))
        }
    }
    fn reduce(&self, e: &Expr) -> PResult<Expr> {
        let (p, q) = self.element(e)?;
        Ok(self.element_expr(&p, &q))
    }
    /// Convert every coefficient of `p` into this ground, failing with the
    /// reference `CoercionFailed`.
    pub fn convert_poly(&self, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
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
        let p = ExpressionPoly::from_expr(e, gens)?;
        self.convert_poly(&p)
    }

    // ------------------------------------------------------------------
    // flattened representation
    // ------------------------------------------------------------------

    /// `(num, den)`: `p = num/den` with `num` over the generators followed
    /// by the symbols and `den` the lcm of the coefficient denominators.
    pub fn flatten(&self, p: &ExpressionPoly) -> PResult<(P, P)> {
        let _s = scope();
        let (k, m) = (p.gens.len(), self.m());
        let mut elements = Vec::new();
        let mut den = MPoly::one(m);
        for (mono, c) in p.terms() {
            let (a, b) = self.element(&c)?;
            if !b.is_one() {
                den = gi::ring_lcm(&den, &b)?;
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
        if self.integral {
            gi::gcd(a, b)
        } else {
            gi::field_gcd(a, b)
        }
    }
    /// `PolyElement.gcd` in `ZZ_I[y...]`/`QQ_I[y...]`. Gaussian elements are
    /// never nonnegative, so `gcd(0, g) = -g`.
    fn egcd(&self, a: &P, b: &P) -> PResult<P> {
        if a.is_zero() {
            return Ok(b.neg());
        }
        if b.is_zero() {
            return Ok(a.neg());
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
        if self.integral {
            Ok(gi::cofactors(f, g)?.0)
        } else {
            let (_, f) = gi::clear_denoms(f);
            let (_, g) = gi::clear_denoms(g);
            Ok(gi::cofactors(&f, &g)?.0)
        }
    }
    /// `dup_rr_prs_gcd`/`dmp_rr_prs_gcd` over the ring ground.
    fn ring_gcd(&self, f: &P, g: &P, k: usize) -> PResult<P> {
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
        let c = if k == 1 { c.scale(&gi::canonical_unit(&lc_ground(&h, k).lc())) } else { c };
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
        let field = self.to_field();
        field.convert_poly(&p.monic()?)
    }
    pub fn gcd(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<ExpressionPoly> {
        let k = f.gens.len();
        let (nf, _) = self.flatten(f)?;
        let (ng, _) = self.flatten(g)?;
        let _s = scope();
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
        let k = f.gens.len();
        let (nf, _) = self.flatten(f)?;
        let (ng, _) = self.flatten(g)?;
        let _s = scope();
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
        let k = p.gens.len();
        let m = self.m();
        if self.ring {
            let (n, _) = self.flatten(p)?;
            let _s = scope();
            return Ok(self.element_expr(&self.ring_content(&n, k)?, &MPoly::one(m)));
        }
        let _s = scope();
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
            let e = self.element_expr(&num, &den);
            cont = self.element(&e)?;
            if cont.0.is_one() && cont.1.is_one() {
                break;
            }
        }
        Ok(cont)
    }
    pub fn primitive(&self, p: &ExpressionPoly) -> PResult<(Expr, ExpressionPoly)> {
        let c = self.content(p)?;
        if p.is_zero() || c.is_one() {
            return Ok((c, self.marked(p)));
        }
        let q = p.exquo_ground(&c)?;
        Ok((c, self.convert_poly(&q)?))
    }
    /// `Poly.clear_denoms()`.
    pub fn clear_denoms(&self, p: &ExpressionPoly) -> PResult<(Expr, ExpressionPoly)> {
        if self.ring {
            return Ok((Expr::one(), self.marked(p)));
        }
        let m = self.m();
        let ring = self.to_ring();
        let _s = scope();
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
        let k = p.gens.len();
        let m = self.m();
        let (num, den) = self.flatten(p)?;
        let _s = scope();
        let one = MPoly::one(m);
        if num.is_zero() {
            return Ok((Expr::zero(), Vec::new()));
        }
        let mut all_gens = p.gens.clone();
        all_gens.extend(self.symbols.iter().cloned());
        let flat = AlgebraicPoly { rep: num.clone(), gens: all_gens, number_field: gi::field(), gaussian: Some(self.integral) };
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
        sort_nested(&mut out, k);
        let factors = out
            .into_iter()
            .map(|(f, i)| Ok((self.eject(&f, &one, &p.gens)?, i as usize)))
            .collect::<PResult<Vec<_>>>()?;
        Ok((coeff, factors))
    }
    /// `dup_sqf_list` over the composite ground.
    pub fn sqf_list(&self, p: &ExpressionPoly, all: bool) -> PResult<(Expr, Vec<(ExpressionPoly, usize)>)> {
        if p.gens.len() != 1 {
            return Err(PolyError::Unsupported(format!(
                "multivariate square-free decomposition over {} is not implemented",
                self.display()
            )));
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
    /// `dmp_cancel(include=False)`: `(cp, cq, p, q)` with `f/g = cp/cq * p/q`.
    pub fn cancel(&self, f: &ExpressionPoly, g: &ExpressionPoly) -> PResult<(Expr, Expr, ExpressionPoly, ExpressionPoly)> {
        if self.ring {
            let (_, p, q) = self.cofactors(f, g)?;
            return Ok((Expr::one(), Expr::one(), p, q));
        }
        let ring = self.to_ring();
        let (cq, f1) = self.clear_denoms(f)?;
        let (cp, g1) = self.clear_denoms(g)?;
        let (_, p, q) = ring.cofactors(&ring.convert_poly(&f1)?, &ring.convert_poly(&g1)?)?;
        let _s = scope();
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
    /// `dup_sqf_part`/`dmp_sqf_part`.
    pub fn sqf_part(&self, p: &ExpressionPoly) -> PResult<ExpressionPoly> {
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
fn element_sort_key(p: &P) -> (usize, Vec<(Vec<u32>, Vec<Rat>)>) {
    (p.len(), p.iter().map(|(e, c)| (e.to_vec(), gi::element_key(c))).collect())
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum NestedKey {
    Leaf((usize, Vec<(Vec<u32>, Vec<Rat>)>)),
    Node(Vec<NestedKey>),
}

fn nested_key(p: &P, k: usize) -> NestedKey {
    if k == 0 {
        return NestedKey::Leaf(element_sort_key(p));
    }
    NestedKey::Node(dense(p).iter().map(|c| nested_key(&drop_vars(c, 1), k - 1)).collect())
}

/// The reference `_sort_factors(multiple=True)` over a composite ground.
fn sort_nested(factors: &mut [(P, u32)], k: usize) {
    factors.sort_by_cached_key(|(p, i)| (dense(p).len(), *i, nested_key(p, k)));
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
}

pub fn descriptor(p: &super::fraction::Polynomial) -> Option<Descriptor> {
    if let Some(e) = &p.expression {
        return e.ground.as_ref().map(|g| Descriptor { gaussian: true, integral: g.integral, ring: g.ring, symbols: g.symbols.clone() });
    }
    if let Some(a) = &p.algebraic {
        return a.gaussian.map(|integral| Descriptor { gaussian: true, integral, ring: true, symbols: Vec::new() });
    }
    if p.modulus.is_some() {
        return None;
    }
    if let Some(f) = &p.field {
        return match &f.coefficient_domain {
            ExactDomain::Functions { integral, ring, parameters, .. } => {
                Some(Descriptor { gaussian: false, integral: *integral, ring: *ring, symbols: parameters.clone() })
            }
            ExactDomain::ZZ => Some(Descriptor { gaussian: false, integral: true, ring: true, symbols: Vec::new() }),
            ExactDomain::QQ => Some(Descriptor { gaussian: false, integral: false, ring: true, symbols: Vec::new() }),
            _ => None,
        };
    }
    let integral = match p.base.domain {
        super::polyobj::Domain::ZZ => true,
        super::polyobj::Domain::QQ => false,
        super::polyobj::Domain::RR => return None,
    };
    let symbols = p.base.domain_gens().to_vec();
    // radicals select `EX` (or an algebraic field) in the reference
    if symbols.iter().any(|s| *s == Expr::i() || (s.is_number() && s.is_pow())) {
        return None;
    }
    Some(Descriptor { gaussian: false, integral, ring: true, symbols })
}

fn descriptor_display(d: &Descriptor) -> String {
    let ground = match (d.gaussian, d.integral) {
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
    if !a.gaussian && !b.gaussian {
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
        return Ok(None);
    }
    Ok(Some(GaussianGround { integral, ring: !frac, symbols }))
}

/// Whether `e` is an algebraic number in the reference sense of
/// `is_number and is_algebraic`, for the expressions built from numbers,
/// `I` and rational powers; anything else counts as transcendental.
fn is_algebraic_number(e: &Expr) -> bool {
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
    Some(GaussianGround { integral, ring: !fractions, symbols })
}
