//! Exact real logarithm conversion for rational linear/quadratic residues.
//!
//! Enumerate the same real roots or the positive-imaginary conjugate as the
//! general resultant algorithm. Keep each logarithm and each Euclidean atan
//! term: a derivative certificate alone would not fix their branch constants.
use crate::expr::Expr;
use crate::functions::{atan, log};
use crate::ops::{add, mul2, sqrt};
use std::sync::Arc;
use crate::poly::{algebraic::{AlgebraicScope, NumberField, AN}, convert, dense,
    int::{Coeff, Field, Rat}, mpoly::MPoly, fraction::Polynomial,
    domains::ExactDomain, polyobj::Domain};

trait RealCoefficient: Field {
    fn expression(&self) -> Expr;
    fn rational(q: &Rat) -> Self;
}
impl RealCoefficient for Rat {
    fn expression(&self) -> Expr { Expr::num(self.to_number()) }
    fn rational(q: &Rat) -> Self { q.clone() }
}
impl RealCoefficient for AN {
    fn expression(&self) -> Expr { self.as_expr() }
    fn rational(q: &Rat) -> Self { AN::rational(q.clone(),AN::zero().field) }
}

fn polynomial<C: RealCoefficient>(p: &[C], x: &Expr) -> Expr {
    add(&p.iter().enumerate().filter(|(_,c)|!c.is_zero()).map(|(i,c)|
        mul2(&c.expression(),&x.powi(i as i64))).collect::<Vec<_>>())
}

fn conjugate<C: RealCoefficient,E>(h: &MPoly<Rat>,u: C,v: C,x: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Expr>,E>
{
    let mut powers=vec![(C::one(),C::zero())];
    for j in 0..h.degree(1).max(0) as usize {
        check()?;
        let (a,b)=&powers[j];
        powers.push((a.mul(&u).sub(&b.mul(&v)),a.mul(&v).add(&b.mul(&u))));
    }
    let length=h.degree(0).max(0) as usize+1;
    let (mut a,mut b)=(vec![C::zero();length],vec![C::zero();length]);
    for (e,c) in h.iter() {
        check()?;
        let coefficient=C::rational(c);
        a[e[0] as usize]=a[e[0] as usize].add(&powers[e[1] as usize].0.mul(&coefficient));
        b[e[0] as usize]=b[e[0] as usize].add(&powers[e[1] as usize].1.mul(&coefficient));
    }
    dense::trim(&mut a);dense::trim(&mut b);
    // Preserve the old zero-polynomial/error paths through the full adapter.
    let Some(arguments)=super::log_to_atan::arguments_field(&a,&b,check)? else{return Ok(None)};
    check()?;
    let norm=dense::add(&dense::mul(&a,&a),&dense::mul(&b,&b));
    let terms=arguments.into_iter().map(|p|mul2(&Expr::int(2),&atan(&polynomial(&p,x)))).collect::<Vec<_>>();
    let result=add(&[mul2(&u.expression(),&log(&polynomial(&norm,x))),mul2(&v.expression(),&add(&terms))]);
    check()?;
    Ok(Some(result))
}

/// Read the exact stored polynomial in (x,t), shared by the public binding
/// and the rational integration driver. General fraction/expression domains
/// never pass through their legacy metadata view.
pub fn rational_input<E>(h: &Polynomial,x: &Expr,t: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<MPoly<Rat>>,E>
{
    if h.algebraic.is_some() || h.expression.is_some() || h.base.modulus.is_some()
        || h.base.ngens!=1 || h.base.true_gens()!=[x.clone()] {return Ok(None);}
    let polynomial=if let Some(field)=&h.field {
        if !matches!(&field.coefficient_domain,ExactDomain::Functions{ring:true,parameters,..} if parameters==&[t.clone()])
            || field.parameters!=[t.clone()] || field.gens!=[x.clone()] || field.rep.degree(0)>32 {return Ok(None);}
        let mut terms=Vec::new();
        for (e,c) in field.rep.iter() {
            check()?;
            let Some(d)=c.d.constant_value() else{return Ok(None)};
            for (parameter,n) in c.n.iter() {
                if terms.len()>=1024 || parameter[0]>32{return Ok(None);}
                terms.push((vec![e[0],parameter[0]],Rat::new(n.clone(),d.clone())));
            }
        }
        MPoly::from_terms(2,terms)
    }else{
        if !matches!(h.base.domain,Domain::ZZ|Domain::QQ) || h.base.gens.iter().any(|g|g!=x && g!=t)
            || h.base.rep.len()>1024 || h.base.rep.degrees().iter().any(|&d|d>32){return Ok(None);}
        let parameter=h.base.gens.iter().position(|g|g==t);
        MPoly::from_terms(2,h.base.rep.iter().map(|(e,c)|
            (vec![e[0],parameter.map_or(0,|i|e[i])],c.clone())).collect())
    };
    Ok(Some(polynomial))
}

/// `h` has exact rational coefficients in `(x,t)`, and `q` is ascending
/// in `t`. Degree/bit limits select only this route; refusals retain the
/// complete general algorithm. A signal is always propagated as an error.
pub fn convert<E>(h: &MPoly<Rat>,q: &[Rat],x: &Expr,t: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Expr>,E>
{
    check()?;
    if q.len()>3 && q.len()<=17 && q.last().is_some_and(|c|!c.is_zero()) {
        return general(h,q,x,t,check);
    }
    if h.nvars!=2 || h.is_zero() || !matches!(q.len(),2|3)
        || q.last().is_some_and(Rat::is_zero) || h.degree(0)>32 || h.degree(1)>32
        || h.len()>1024 || h.coeffs().iter().chain(q).any(|c|c.n.bits().max(c.d.bits())>4096) {
        return Ok(None);
    }
    let expression=convert::to_expr(h,&[x.clone(),t.clone()]);
    let real_term=|root:&Expr|mul2(root,&log(&crate::subs::subs1(&expression,t,root)));
    if q.len()==2 {
        let root=q[0].neg().div(&q[1]).expression();
        check()?;
        let result=real_term(&root);check()?;return Ok(Some(result));
    }
    let four=Rat::from_i64(4);let two=Rat::from_i64(2);
    let discriminant=q[1].mul(&q[1]).sub(&four.mul(&q[2]).mul(&q[0]));
    if discriminant.is_zero(){return Ok(None);}
    let u=q[1].neg().div(&two.mul(&q[2]));
    let radius=discriminant.div(&four.mul(&q[2].mul(&q[2])));
    if !discriminant.is_negative() {
        let v=sqrt(&radius.expression());
        let result=add(&[real_term(&add(&[u.expression(),v.clone()])),
            real_term(&add(&[u.expression(),crate::ops::neg(&v)]))]);
        check()?;return Ok(Some(result));
    }
    let radius=radius.neg();
    let (n,d)=(radius.n.isqrt(),radius.d.isqrt());
    check()?;
    if n.mul(&n)==radius.n && d.mul(&d)==radius.d {
        return conjugate(h,u,Rat::new(n,d),x,check);
    }
    let v=sqrt(&radius.expression());
    let Ok(field)=NumberField::quadratic(&v) else{return Ok(None)};
    let _scope=AlgebraicScope::new(field.clone());
    conjugate(h,AN::rational(u,field.clone()),AN::alpha(field),x,check)
}

/// Real and imaginary parts of `p(.., u+I*v)` for the last variable of `p`.
/// The output variables are the leading ones of `p`, then `u`, then `v`.
fn split(p: &MPoly<Rat>) -> (MPoly<Rat>,MPoly<Rat>) {
    let n=p.nvars;let (mut re,mut im)=(Vec::new(),Vec::new());
    for (e,c) in p.iter() {
        let k=e[n-1];let mut binomial=Rat::one();
        for j in 0..=k {
            if j>0 {binomial=binomial.mul(&Rat::from_i64((k-j+1) as i64)).div(&Rat::from_i64(j as i64));}
            let mut exponent=e[..n-1].to_vec();exponent.push(k-j);exponent.push(j);
            let mut value=c.mul(&binomial);
            if j%4>=2 {value=value.neg();}
            if j%2==0 {re.push((exponent,value));} else {im.push((exponent,value));}
        }
    }
    (MPoly::from_terms(n+1,re),MPoly::from_terms(n+1,im))
}

/// `p` in leading variables (optionally `x`), then `u`, `v`, evaluated at
/// exact `(u,v)`: ascending coefficients in the leading variable.
fn evaluate<C: RealCoefficient>(p: &MPoly<Rat>,u: &C,v: &C) -> Vec<C> {
    let n=p.nvars;
    let length=if n==3 {p.degree(0).max(0) as usize+1} else {1};
    let mut result=vec![C::zero();length];
    for (e,c) in p.iter() {
        let mut value=C::rational(c);
        for _ in 0..e[n-2] {value=value.mul(u);}
        for _ in 0..e[n-1] {value=value.mul(v);}
        let i=if n==3 {e[0] as usize} else {0};
        result[i]=result[i].add(&value);
    }
    dense::trim(&mut result);result
}

/// The reference `roots(f, filter='R')` and the `count_roots` completeness
/// test of `_get_real_roots`, for a rational univariate `f`. `None` is an
/// incomplete root list: the caller declines instead of returning a sum.
fn real_roots(f: &crate::poly::polyobj::Poly) -> Option<Vec<Expr>> {
    let roots=super::super::poly::roots::roots_poly(f).ok()?;
    let real=roots.into_iter().map(|(r,_)|r).filter(|r|r.is_zero() || {
        let (n,d)=crate::simptools::as_numer_denom(r);
        crate::assumptions::is_real(&n)==Some(true) && crate::assumptions::is_real(&d)==Some(true)
    }).collect::<Vec<_>>();
    let coefficients=dense::from_mpoly(&f.rep,0).into_iter().rev().collect::<Vec<_>>();
    (real.len()==crate::poly::rootisolation::count_real_roots(&coefficients,None,None)).then_some(real)
}

/// Whether every real root `r_u` of `f` has degree at most 2 over QQ and
/// gives a rational `c(r_u, v)` (each `v` coefficient of `c`, a polynomial
/// in `u`, reduces to a constant modulo the minimal polynomial of `r_u`).
/// Any other root almost always gives an EX `c(r_u, v)`; declining here is
/// always safe and leaves the costly root computation to the general path
/// alone instead of doing it twice.
fn decidable_resultant_roots(f: &crate::poly::polyobj::Poly,c: &MPoly<Rat>) -> bool {
    if f.rep.nvars!=1 || c.nvars!=2 {return false;}
    let mut coefficients:Vec<Vec<Rat>>=vec![Vec::new();c.degree(1).max(0) as usize+1];
    for (e,value) in c.iter() {
        let row=&mut coefficients[e[1] as usize];
        if row.len()<=e[0] as usize {row.resize(e[0] as usize+1,Rat::zero());}
        row[e[0] as usize].add_assign(value);
    }
    for row in &mut coefficients {dense::trim(row);}
    let (_,integral)=f.rep.clear_denoms();
    let (_,factors)=crate::poly::factor_multi::factor_list_mz(&integral);
    factors.iter().all(|(g,_)|{
        let g=dense::from_mpoly(&g.to_rat(),0);
        if g.len()<=2 {return true;}
        let descending=g.iter().rev().cloned().collect::<Vec<_>>();
        if crate::poly::rootisolation::count_real_roots(&descending,None,None)==0 {return true;}
        g.len()==3 && coefficients.iter().all(|row|dense::rem(row,&g).len()<=1)
    })
}

/// Whether no irreducible factor of degree >= 5 of `f` has a real root.
/// Such roots are generally not found by `roots` (an incomplete list, where
/// the reference returns None) and never lie in a quadratic field, so the
/// costly attempt is left to the general path alone (always safe).
fn findable_real_roots(f: &crate::poly::polyobj::Poly) -> bool {
    if f.rep.nvars!=1 {return false;}
    let (_,integral)=f.rep.clear_denoms();
    let (_,factors)=crate::poly::factor_multi::factor_list_mz(&integral);
    factors.iter().all(|(g,_)|g.degree(0)<5 || {
        let descending=dense::from_mpoly(&g.to_rat(),0).into_iter().rev().collect::<Vec<_>>();
        crate::poly::rootisolation::count_real_roots(&descending,None,None)==0
    })
}

fn rational_poly(e: &Expr,v: &Expr) -> Option<crate::poly::polyobj::Poly> {
    let p=crate::poly::polyobj::Poly::from_expr(e,std::slice::from_ref(v)).ok()?;
    (matches!(p.domain,Domain::ZZ|Domain::QQ) && p.modulus.is_none() && !p.has_domain_gens()
        && p.gens==[v.clone()] && p.ngens==1).then_some(p)
}

/// The common rational or pure quadratic field `QQ(sqrt(k))` of `values`.
/// `Ok(None)` is a field outside these; `Some(None)` is QQ.
fn common_field(values: &[&Expr]) -> Option<Option<Arc<NumberField>>> {
    if values.iter().all(|v|v.is_rational_number()) {return Some(None);}
    let field=crate::poly::algebraic::infer(&values.iter().map(|v|(*v).clone()).collect::<Vec<_>>()).ok()??;
    // The reference EX elements are expanded `p + q*sqrt(k)` sums only for
    // this generator; any other primitive element changes the printed form.
    (field.degree()==2 && field.is_pure_quadratic()
        && field.generator==sqrt(&Expr::num(field.radicand.to_number()))).then_some(Some(field))
}

/// `|value| > 1e-10`, decided with a certified relative error.
fn clearly_nonzero_an(value: &AN) -> bool {
    let c=&value.coefficients;
    let p=c.first().cloned().unwrap_or_else(Rat::zero);
    let q=c.get(1).cloned().unwrap_or_else(Rat::zero);
    let k=value.field.radicand.clone();
    // |p+q*sqrt(k)| = |p^2-q^2*k| / |p-q*sqrt(k)| without cancellation.
    let (pf,qf,kf)=(p.to_f64(),q.to_f64(),k.to_f64());
    if !(pf.is_finite()&&qf.is_finite()&&kf.is_finite()) || kf<=0.0 {return false;}
    let magnitude=if (pf>=0.0)==(qf>=0.0) {pf.abs()+qf.abs()*kf.sqrt()} else {
        let norm=p.mul(&p).sub(&q.mul(&q).mul(&k)).to_f64().abs();
        norm/(pf.abs()+qf.abs()*kf.sqrt())
    };
    magnitude.is_finite() && magnitude>1e-10
}

/// The numeric reference test `D.evalf(chop=True) != 0` where `D` is not
/// in a common quadratic field: only a value far above the evaluation error
/// is decided (as nonzero); anything else is undecided.
fn clearly_nonzero_numeric(d: &MPoly<Rat>,u: &Expr,v: &Expr) -> bool {
    let (Some(u),Some(v))=(crate::evalf::eval_f64(u),crate::evalf::eval_f64(v)) else{return false};
    let (mut value,mut bound)=(0.0f64,0.0f64);
    for (e,c) in d.iter() {
        let term=c.to_f64()*u.powi(e[0] as i32)*v.powi(e[1] as i32);
        value+=term;bound+=term.abs();
    }
    value.is_finite() && bound.is_finite() && value.abs()>1e-10_f64.max(1e-8*bound)
}

fn pair_terms<C: RealCoefficient,E>(a: &MPoly<Rat>,b: &MPoly<Rat>,d: &MPoly<Rat>,
    (u,v): (&C,&C),(ru,rv): (&Expr,&Expr),x: &Expr,nonzero: impl Fn(&C)->bool,
    rational: impl Fn(&C)->bool,check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Option<Expr>>,E>
{
    let value=evaluate(d,u,v).pop().unwrap_or_else(C::zero);
    if !value.is_zero() {
        // Certainly nonzero for the reference numeric test: skipped there.
        return Ok(if nonzero(&value) {Some(None)} else {None});
    }
    check()?;
    let (f,g)=(evaluate(a,u,v),evaluate(b,u,v));
    let Some(arguments)=super::log_to_atan::arguments_divisors(&f,&g,&rational,check)? else{return Ok(None)};
    check()?;
    let norm=dense::add(&dense::mul(&f,&f),&dense::mul(&g,&g));
    let terms=arguments.into_iter().map(|p|mul2(&Expr::int(2),&atan(&polynomial(&p,x)))).collect::<Vec<_>>();
    Ok(Some(Some(add(&[mul2(ru,&log(&polynomial(&norm,x))),mul2(rv,&add(&terms))]))))
}

/// For `q` of odd degree `n`, the row of `v^(n-1)` in `c` is
/// `+-(n*q_n*u + q_(n-1))`, never constant modulo a quadratic factor, so the
/// pre-check above declines as soon as `R(u)` has an irrational real root.
/// Every real root of `q` and every real part of a root of `q` is such a
/// root of `R` (`R` lies in the ideal of `c`, `d`). An irreducible factor
/// `g` of `q` of degree >= 2 has an irrational real root when it has any
/// real root; when all of its roots have rational real parts, its root set
/// is invariant under `t -> 2*r - t` (Galois conjugation keeps the sum of a
/// conjugate pair) with `r` the mean of its roots. So the decline is certain
/// when some factor has a real root or `g(2*r - t) != (-1)^m * g(t)`; this
/// costs one univariate factorization instead of the resultant.
fn odd_degree_certain_decline(q: &[Rat]) -> bool {
    if q.len()%2!=0 {return false;}
    let qpoly=MPoly::from_terms(1,q.iter().enumerate().map(|(i,c)|(vec![i as u32],c.clone())).collect());
    let (_,integral)=qpoly.clear_denoms();
    let (_,factors)=crate::poly::factor_multi::factor_list_mz(&integral);
    factors.iter().any(|(g,_)|{
        let g=dense::from_mpoly(&g.to_rat(),0);
        if g.len()<=2 {return false;}
        let descending=g.iter().rev().cloned().collect::<Vec<_>>();
        if crate::poly::rootisolation::count_real_roots(&descending,None,None)!=0 {return true;}
        let m=g.len()-1;
        let r=g[m-1].neg().div(&g[m].mul(&Rat::from_i64(m as i64)));
        let reflection=vec![r.add(&r),Rat::from_i64(-1)];
        let mut composed:Vec<Rat>=Vec::new();
        for coefficient in descending {
            composed=dense::add(&dense::mul(&composed,&reflection),&[coefficient]);
        }
        let expected=if m%2==0 {g.clone()} else {dense::neg(&g)};
        composed!=expected
    })
}

thread_local! {
    /// The last declined `(h, q)` pairs: a decline depends only on them,
    /// so the public fallback after a declined driver attempt does not pay
    /// for the same refusal twice. Errors (signals) are never recorded.
    static DECLINED: std::cell::RefCell<Vec<(MPoly<Rat>,Vec<Rat>)>>=const{std::cell::RefCell::new(Vec::new())};
}

fn general<E>(h: &MPoly<Rat>,q: &[Rat],x: &Expr,t: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Expr>,E>
{
    if DECLINED.with(|d|d.borrow().iter().any(|(dh,dq)|dh==h && dq.as_slice()==q)) {return Ok(None);}
    let result=general_attempt(h,q,x,t,check)?;
    if result.is_none() {
        DECLINED.with(|d|{
            let mut d=d.borrow_mut();
            if d.len()>=8 {d.remove(0);}
            d.push((h.clone(),q.to_vec()));
        });
    }
    Ok(result)
}

/// The reference general `log_to_real` for a rational auxiliary polynomial
/// `q` of degree >= 3: resultant `R(u)`, real roots of `R` and of each
/// `c(r_u, v)` with the `count_roots` completeness test, the `D` test and
/// `log_to_atan`. Only exactly decided rational and pure quadratic steps are
/// computed here; every other step declines to the unchanged general path.
fn general_attempt<E>(h: &MPoly<Rat>,q: &[Rat],x: &Expr,t: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Expr>,E>
{
    if h.nvars!=2 || h.is_zero() || h.degree(0)>32 || h.degree(1)>32 || h.len()>1024
        || h.coeffs().iter().chain(q).any(|c|c.n.bits().max(c.d.bits())>4096) {return Ok(None);}
    if odd_degree_certain_decline(q) {return Ok(None);}
    check()?;
    let (u,v)=(Expr::dummy("u",crate::symbol::Assumptions::NONE),Expr::dummy("v",crate::symbol::Assumptions::NONE));
    let (a,b)=split(h);
    let qpoly=MPoly::from_terms(1,q.iter().enumerate().map(|(i,c)|(vec![i as u32],c.clone())).collect());
    let (c,d)=split(&qpoly);
    check()?;
    let (ce,de)=(convert::to_expr(&c,&[u.clone(),v.clone()]),convert::to_expr(&d,&[u.clone(),v.clone()]));
    let Ok(resultant)=crate::poly::api::resultant(&ce,&de,std::slice::from_ref(&v)) else{return Ok(None)};
    check()?;
    let Some(r)=rational_poly(&resultant,&u) else{return Ok(None)};
    if !decidable_resultant_roots(&r,&c) {return Ok(None);}
    check()?;
    let Some(real_u)=real_roots(&r) else{return Ok(None)};
    let mut result=Vec::new();
    for ru in &real_u {
        check()?;
        let mut map=rustc_hash::FxHashMap::default();map.insert(u.clone(),ru.clone());
        let Some(cv)=rational_poly(&crate::subs::xreplace(&ce,&map),&v) else{return Ok(None)};
        if cv.rep.is_zero() || !findable_real_roots(&cv) {return Ok(None);}
        check()?;
        let Some(real_v)=real_roots(&cv) else{return Ok(None)};
        check()?;
        let mut paired:Vec<Expr>=Vec::new();
        for rv in real_v {
            let negative=crate::ops::neg(&rv);
            if paired.contains(&rv) || paired.contains(&negative) {continue;}
            if crate::assumptions::fact(&rv,"negative")==Some(true) || crate::decompose::could_extract_minus_sign(&rv) {
                paired.push(negative);
            } else if crate::assumptions::fact(&rv,"zero")!=Some(true) {
                paired.push(rv);
            }
        }
        for rv in &paired {
            check()?;
            let term=match common_field(&[ru,rv]) {
                Some(None)=>{
                    let (Some(un),Some(vn))=(ru.as_number().and_then(Rat::from_number),rv.as_number().and_then(Rat::from_number))
                        else{return Ok(None)};
                    pair_terms(&a,&b,&d,(&un,&vn),(ru,rv),x,|c:&Rat|c.to_f64().abs()>1e-10,|_|true,check)?
                },
                Some(Some(field))=>{
                    let _scope=AlgebraicScope::new(field.clone());
                    let (Ok(un),Ok(vn))=(AN::from_expr(ru,field.clone()),AN::from_expr(rv,field.clone())) else{return Ok(None)};
                    pair_terms(&a,&b,&d,(&un,&vn),(ru,rv),x,clearly_nonzero_an,
                        // A rational or pure `q*sqrt(k)` divisor keeps the
                        // EX `cancel`+`expand` quotient in `p + q*sqrt(k)` form.
                        |c:&AN|c.coefficients.len()<=1 || (c.coefficients.len()==2 && c.coefficients[0].is_zero()),check)?
                },
                None=>if clearly_nonzero_numeric(&d,ru,rv) {Some(None)} else {None},
            };
            match term {
                None=>return Ok(None),
                Some(None)=>{},
                Some(Some(term))=>result.push(term),
            }
        }
    }
    check()?;
    let Some(qp)=rational_poly(&convert::to_expr(&qpoly,std::slice::from_ref(t)),t) else{return Ok(None)};
    let Some(real_t)=real_roots(&qp) else{return Ok(None)};
    let expression=convert::to_expr(h,&[x.clone(),t.clone()]);
    for root in real_t {
        check()?;
        result.push(mul2(&root,&log(&crate::subs::subs1(&expression,t,&root))));
    }
    let result=add(&result);
    check()?;
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn rational_quadratic_roots_and_conjugates_keep_representatives() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::one())]);
        for (q,expected) in [(vec![Rat::one(),Rat::one()],"-log(x-1)"),
            (vec![Rat::one(),Rat::zero(),Rat::one()],"2*atan(x)"),
            (vec![Rat::one(),Rat::zero(),Rat::from_i64(3)],"2*sqrt(3)*atan(sqrt(3)*x)/3")] {
            let result=convert(&h,&q,&x,&t,&mut ||Ok::<_,()>(())).unwrap().unwrap();
            assert_eq!(result,parse_str(expected).unwrap());
        }
    }
    #[test]
    fn repeated_roots_decline_and_signals_never_return_partial_sums() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::one())]);
        assert!(convert(&h,&[Rat::one(),Rat::from_i64(2),Rat::one()],&x,&t,&mut ||Ok::<_,()>(())).unwrap().is_none());
        let mut calls=0;
        let result=convert(&h,&[Rat::one(),Rat::zero(),Rat::from_i64(3)],&x,&t,&mut ||{calls+=1;if calls==6{Err("stop")}else{Ok(())}});
        assert!(matches!(result,Err("stop")));
    }
    #[test]
    fn odd_degree_structural_decline_is_exact() {
        let q=|v:&[i64]|v.iter().map(|&c|Rat::from_i64(c)).collect::<Vec<_>>();
        // 823543*t^7+1: irrational real parts; t^3-2: irrational real root.
        assert!(odd_degree_certain_decline(&q(&[1,0,0,0,0,0,0,823543])));
        assert!(odd_degree_certain_decline(&q(&[-2,0,0,1])));
        // (t+1)*(t^2+1)*(t^2+4): every real part rational; even degree.
        assert!(!odd_degree_certain_decline(&q(&[4,4,5,5,1,1])));
        // (t-1)*((t-2)^2+1)*((t-2)^2+3) = roots with real parts 1 and 2.
        let a=[5i64,-4,1];let b=[7i64,-4,1];
        let mut ab=vec![0i64;5];for i in 0..3 {for j in 0..3 {ab[i+j]+=a[i]*b[j];}}
        let mut f=vec![0i64;6];for i in 0..5 {f[i]-=ab[i];f[i+1]+=ab[i];}
        assert!(!odd_degree_certain_decline(&q(&f)));
        assert!(!odd_degree_certain_decline(&q(&[1,0,0,0,256])));
    }

    #[test]
    fn general_degree_four_quadratic_pairs_and_declines() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        // log part of 1/(x**4+1): h = x + 4*t, q = 256*t**4 + 1.
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::from_i64(4))]);
        let q=vec![Rat::one(),Rat::zero(),Rat::zero(),Rat::zero(),Rat::from_i64(256)];
        let result=convert(&h,&q,&x,&t,&mut ||Ok::<_,()>(())).unwrap().unwrap();
        let text=format!("{result}");
        assert!(text.contains("atan") && text.contains("log") && text.contains("sqrt(2)"),"{text}");
        // Cube-root resultant roots (an EX c(r_u, v)) decline.
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::one())]);
        let cubic=vec![Rat::from_i64(-2),Rat::zero(),Rat::zero(),Rat::one()];
        assert!(convert(&h,&cubic,&x,&t,&mut ||Ok::<_,()>(())).unwrap().is_none());
        // Interruptions propagate and never publish a partial sum.
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::from_i64(4))]);
        for stop in 1..12 {
            let mut calls=0;
            let result=convert(&h,&q,&x,&t,&mut ||{calls+=1;if calls==stop{Err("stop")}else{Ok(())}});
            assert!(matches!(result,Err("stop")) || calls<stop);
        }
    }
}
