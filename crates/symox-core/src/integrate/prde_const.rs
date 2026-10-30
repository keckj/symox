//! `constant_system(A, u, DE)` of SymPy `prde.py` over rational function
//! fields `K(x, t0, ...)`.
//!
//! The values are reduced rational functions (`poly::ratfunc::RF`) in one
//! variable set `U` that holds the PolyMatrix domain parameters, the tower
//! variables `DE.T` and the symbols of `DE.D`. The reduced row echelon form
//! is the exact fraction-field one (`ddm::field_irref`), and the derivation
//! is the exact one, `D(f) = sum_k d_k * df/dt_k` over the active levels.
//!
//! Besides the values, the stock body exposes expressions (the PolyMatrix
//! items) and PolyMatrix domains: the items are rebuilt with the same
//! operations as the stock body (`a - b` stays uncancelled, products are
//! cancelled, `D(x)/D(a_ij)` is the quotient of two cancelled expressions)
//! and the domains are computed by the [`Host`] callbacks, which reproduce the
//! `PolyMatrix` inference (`parallel_poly_from_expr(field=True)`), `unify`
//! and membership (`Poly(e, *gens, domain=K)`). Any step the kernel cannot
//! reproduce returns `None` so that the stock body runs unchanged.
use crate::expr::Expr;
use crate::ops;
use crate::poly::int::{Coeff, Field};
use crate::poly::mpoly::MPoly;
use crate::poly::ratfunc::RF;
use crate::poly::Int;

/// PolyMatrix domain operations of the host.
pub trait Host {
    type D: Clone;
    /// Domain of `PolyMatrix(items, gens)`, `None` when not reproduced.
    fn infer(&mut self, items: &[Expr]) -> Option<Self::D>;
    /// `PolyMatrix._unified`: equal domains, `unify`, else the receiver.
    fn unify(&mut self, a: &Self::D, b: &Self::D) -> Self::D;
    /// Canonical expression of `v` as an element of `d` (`Poly(e, *gens,
    /// domain=d)`), `None` when `v` is not an element of `d`.
    fn element(&mut self, d: &Self::D, v: &RF) -> Option<Expr>;
    /// The domain is a field (`rref` ground check).
    fn is_field(&mut self, d: &Self::D) -> bool;
    /// The expression of `v` in the variables `U`.
    fn expr(&mut self, v: &RF) -> Expr;
}

/// The derivation of the active levels and the variables of `DE.T`.
pub struct Tower {
    /// `(index in U, D t_k as a rational function)` for the active levels.
    pub active: Vec<(usize, RF)>,
    /// Indices in `U` of every element of `DE.T` (`has(*DE.T)`).
    pub tower: Vec<usize>,
    /// Index in `U` of the matrix generator when it is one of the variables.
    pub gen: Option<usize>,
}

/// A PolyMatrix: values, items and domain.
pub struct Matrix<D> {
    pub rows: usize,
    pub cols: usize,
    pub values: Vec<RF>,
    pub items: Vec<Expr>,
    pub domain: D,
}

fn depends(p: &MPoly<Int>, v: usize) -> bool {
    p.degree(v) > 0
}

fn involves(v: &RF, var: usize) -> bool {
    depends(&v.n, var) || depends(&v.d, var)
}

/// `df/dv` of a reduced fraction.
pub fn partial(f: &RF, v: usize) -> RF {
    if !involves(f, v) {
        return RF { n: MPoly::zero(f.n.nvars), d: MPoly::one(f.n.nvars) };
    }
    let dn = f.n.derivative(v);
    if !depends(&f.d, v) {
        return RF::new(dn, f.d.clone());
    }
    let dd = f.d.derivative(v);
    RF::new(dn.mul(&f.d).sub(&f.n.mul(&dd)), f.d.mul(&f.d))
}

/// `derivation(f, DE, basic=True)` as a rational function.
pub fn derivation(f: &RF, active: &[(usize, RF)]) -> RF {
    let mut r = RF { n: MPoly::zero(f.n.nvars), d: MPoly::one(f.n.nvars) };
    for (v, d) in active {
        if involves(f, *v) {
            r = r.add(&d.mul(&partial(f, *v)));
        }
    }
    r
}

fn cancel<H: Host>(h: &mut H, v: &RF) -> Option<Expr> {
    crate::poly::api::reference_cancel(&h.expr(v)).ok()
}

/// The stock `constant_system` body. `a` and `u` are the input matrices with
/// their values; `None` whenever the stock body would raise or reach a step
/// that is not reproduced.
pub fn constant_system<H: Host>(h: &mut H, a: Matrix<H::D>, u: Matrix<H::D>, tower: &Tower,
    check: &mut impl FnMut() -> bool) -> Option<(Matrix<H::D>, Matrix<H::D>)> {
    let (m, n) = (a.rows, a.cols);
    if m * n == 0 || u.rows != m || u.cols != 1 || a.values.len() != m * n || u.values.len() != m {
        return None;
    }
    let free_of_gen = |v: &RF| tower.gen.is_none_or(|g| !involves(v, g));
    // Au = A.row_join(u); Au.rref() over the unified ground field.
    let dau = h.unify(&a.domain, &u.domain);
    if !h.is_field(&dau) {
        return None;
    }
    let mut au: Vec<Vec<RF>> = Vec::with_capacity(m);
    for i in 0..m {
        let mut row = a.values[i * n..(i + 1) * n].to_vec();
        row.push(u.values[i].clone());
        for v in &row {
            if !free_of_gen(v) || h.element(&dau, v).is_none() {
                return None;
            }
        }
        au.push(row);
    }
    let mut stop = false;
    let pivots = crate::poly::ddm::field_irref(&mut au, true, &mut || {
        if check() { Ok(()) } else { stop = true; Err(()) }
    }, &mut |_| Ok(()));
    if stop || !matches!(pivots, Ok(Some(_))) {
        return None;
    }
    let mut av = Vec::with_capacity(m * n);
    let mut ai = Vec::with_capacity(m * n);
    let mut uv = Vec::with_capacity(m);
    let mut ui = Vec::with_capacity(m);
    for row in au {
        for (k, v) in row.into_iter().enumerate() {
            let e = h.element(&dau, &v)?;
            if k < n { av.push(v); ai.push(e); } else { uv.push(v); ui.push(e); }
        }
    }
    let mut a = Matrix { rows: m, cols: n, values: av, items: ai, domain: dau.clone() };
    let mut u = Matrix { rows: m, cols: 1, values: uv, items: ui, domain: dau };
    for j in 0..n {
        for i in 0..m {
            if !check() {
                return None;
            }
            let aij = a.values[i * n + j].clone();
            // A[i, j] is read as a Poly over A's domain.
            h.element(&a.domain, &aij)?;
            if !tower.tower.iter().any(|&v| involves(&aij, v)) {
                continue;
            }
            let daij = derivation(&aij, &tower.active);
            if daij.is_zero() {
                return None;
            }
            let c = cancel(h, &daij)?;
            // Rm1 = Ri.applyfunc(lambda x: D(x) / DAij)
            let mut rv = Vec::with_capacity(n);
            let mut ri = Vec::with_capacity(n);
            for k in 0..n {
                let dx = derivation(&a.values[i * n + k], &tower.active);
                ri.push(ops::div(&cancel(h, &dx)?, &c));
                rv.push(dx.div(&daij));
            }
            let dr = h.infer(&ri)?;
            // um1 = D(u[i]) / DAij
            h.element(&u.domain, &u.values[i])?;
            let du = derivation(&u.values[i], &tower.active);
            let um1e = ops::div(&cancel(h, &du)?, &c);
            let um1v = du.div(&daij);
            // A = A - Aj * Rm1; u = u - Aj * um1
            let dp = h.unify(&a.domain, &dr);
            let da_new = h.unify(&a.domain, &dp);
            let du_new = h.unify(&u.domain, &a.domain);
            let mut nav = Vec::with_capacity((m + 1) * n);
            let mut nai = Vec::with_capacity((m + 1) * n);
            let rows = a.rows;
            for r in 0..rows {
                let ajr = a.values[r * n + j].clone();
                for k in 0..n {
                    let p = ajr.mul(&rv[k]);
                    let pe = cancel(h, &p)?;
                    nav.push(a.values[r * n + k].sub(&p));
                    nai.push(ops::sub(&a.items[r * n + k], &pe));
                }
            }
            let mut nuv = Vec::with_capacity(rows + 1);
            let mut nui = Vec::with_capacity(rows + 1);
            for r in 0..rows {
                let q = a.values[r * n + j].mul(&um1v);
                let qe = cancel(h, &q)?;
                nuv.push(u.values[r].sub(&q));
                nui.push(ops::sub(&u.items[r], &qe));
            }
            // A = A.col_join(Rm1); u = u.col_join(Matrix([um1], u.gens))
            let da = h.unify(&da_new, &dr);
            let dm = h.infer(std::slice::from_ref(&um1e))?;
            let dum = h.unify(&du_new, &dm);
            nav.extend(rv);
            nai.extend(ri);
            nuv.push(um1v);
            nui.push(um1e);
            if nav.iter().chain(&nuv).any(|v| !free_of_gen(v)) {
                return None;
            }
            a = Matrix { rows: rows + 1, cols: n, values: nav, items: nai, domain: da };
            u = Matrix { rows: rows + 1, cols: 1, values: nuv, items: nui, domain: dum };
        }
    }
    Some((a, u))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::ratfunc::Scope;

    fn x() -> MPoly<Int> {
        MPoly::var(0, 2)
    }
    fn t() -> MPoly<Int> {
        MPoly::var(1, 2)
    }
    fn c(k: i64) -> MPoly<Int> {
        MPoly::constant(Int::S(k), 2)
    }

    #[test]
    fn prde_const_partial_and_derivation() {
        let _s = Scope::new(2);
        // f = t/x, D = [1, t] (t = exp(x)): Df = t/x - t/x**2.
        let f = RF::new(t(), x());
        assert_eq!(partial(&f, 0), RF::new(t().neg(), x().mul(&x())));
        let active = vec![(0, RF::from_poly(c(1))), (1, RF::from_poly(t()))];
        let df = derivation(&f, &active);
        assert_eq!(df, RF::new(t().mul(&x()).sub(&t()), x().mul(&x())));
        // constants have a zero derivative
        assert!(derivation(&RF::new(c(3), c(7)), &active).is_zero());
    }

    struct Plain;
    impl Host for Plain {
        type D = u8;
        fn infer(&mut self, _: &[Expr]) -> Option<u8> { Some(1) }
        fn unify(&mut self, a: &u8, b: &u8) -> u8 { *a.max(b) }
        fn element(&mut self, _: &u8, _: &RF) -> Option<Expr> { Some(Expr::one()) }
        fn is_field(&mut self, _: &u8) -> bool { true }
        fn expr(&mut self, v: &RF) -> Expr {
            let x = Expr::symbol("x");
            let s = [x];
            ops::div(&crate::poly::convert::to_expr_int(&v.n, &s), &crate::poly::convert::to_expr_int(&v.d, &s))
        }
    }

    #[test]
    fn prde_const_base_case_system() {
        // A = [[1, x]], u = [[0]] over QQ(x), D = [1]: the second column
        // involves x, so the row D(A[0, :])/D(x) = [0, 1] is appended.
        let _s = Scope::new(1);
        let one = RF::from_poly(MPoly::one(1));
        let xv = RF::from_poly(MPoly::var(0, 1));
        let zero = RF { n: MPoly::zero(1), d: MPoly::one(1) };
        let a = Matrix { rows: 1, cols: 2, values: vec![one.clone(), xv.clone()], items: vec![Expr::one(); 2], domain: 1u8 };
        let u = Matrix { rows: 1, cols: 1, values: vec![zero.clone()], items: vec![Expr::zero()], domain: 1u8 };
        let tower = Tower { active: vec![(0, one.clone())], tower: vec![0], gen: None };
        let (a, u) = constant_system(&mut Plain, a, u, &tower, &mut || true).unwrap();
        assert_eq!((a.rows, a.cols), (2, 2));
        // A - Aj*Rm1 = [1, x] - x*[0, 1] = [1, 0]; then [0, 1].
        assert_eq!(a.values, vec![one.clone(), zero.clone(), zero.clone(), one.clone()]);
        assert_eq!(u.values, vec![zero.clone(), zero]);
    }
}
