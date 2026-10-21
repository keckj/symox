//! Python bindings of the polynomial module: the `Poly` class and the
//! functions `factor`, `factor_list`, `sqf_list`, `gcd`, `lcm`, `cancel`,
//! `apart`, `collect`, `div`, `rem`, `quo`, `resultant`, `discriminant`,
//! `roots`, `nroots`, `degree`, ...

#![allow(non_snake_case)]

#[path = "algebraic_number.rs"]
pub(crate) mod algebraic_number;
#[path="an_rational.rs"]
mod an_rational;
#[path="domains.rs"]
pub(crate) mod domains;
#[path="poly_division.rs"]
pub(crate) mod division;
#[path="poly_real_roots.rs"]
mod root_analysis;
#[path="poly_checked_factors.rs"]
mod checked_factor_views;
#[path="poly_structure.rs"]
mod structural;

use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::expr::Expr;
use symox::poly::api;
use symox::poly::convert::PolyError;
use symox::poly::polyobj::{Domain, Poly};
use symox::poly::fraction::{FractionPoly, Polynomial as PolyValue};
use symox::poly::expression_domain::ExpressionPoly;
use symox::poly::algebraic::{AlgebraicPoly,NumberField,AlgebraicScope};
use std::sync::Arc;
use symox::poly::roots as proots;
use symox::poly::Rat;

use crate::{to_expr, PyExpr};

fn w(e: Expr) -> PyExpr {
    PyExpr { e }
}

/// SymPy's `Poly.cancel(g, include=False)` (`dmp_cancel`) for plain
/// `ZZ`/`QQ` polynomials: `(cp/cq, p, q)` where `p/q` is the cancelled
/// fraction scaled to coprime integer polynomials with positive leading
/// coefficients, `cq`/`cp` being the denominators cleared from `f` and `g`.
/// `a/b` is the already cancelled `f/g`.
fn sympy_cancel_parts(f:&PolyValue,g:&PolyValue,a:&PolyValue,b:&PolyValue)->Option<(Rat,PolyValue,PolyValue)>{
    let plain=|p:&PolyValue|p.field.is_none()&&p.algebraic.is_none()&&p.base.modulus.is_none()&&!p.base.has_domain_gens()&&matches!(p.base.domain,Domain::ZZ|Domain::QQ);
    if ![f,g,a,b].into_iter().all(plain)||a.base.gens!=b.base.gens||a.base.domain!=b.base.domain||b.base.rep.is_zero(){return None;}
    let (cq,_)=f.base.rep.clear_denoms();let (cp,_)=g.base.rep.clear_denoms();
    let (da,ia)=a.base.rep.clear_denoms();let (db,ib)=b.base.rep.clear_denoms();
    // a/b = (ia*db)/(ib*da) and p/q = (cq/cp)*(a/b).
    let n=ia.scale(&db.mul(&cq));let m=ib.scale(&da.mul(&cp));
    let g0=n.content().gcd(&m.content());
    let (n,m)=if g0.is_one()||g0.is_zero(){(n,m)}else{(n.div_coeff(&g0)?,m.div_coeff(&g0)?)};
    let mut sign=1i64;
    let n=if !n.is_zero()&&n.lc().is_negative(){sign=-sign;n.neg()}else{n};
    let m=if m.lc().is_negative(){sign=-sign;m.neg()}else{m};
    let h=cp.gcd(&cq);
    let common=Rat::new(cp.div_exact_unchecked(&h).mul_i64(sign),cq.div_exact_unchecked(&h));
    // The cancelled pair keeps the unified input domain (QQ stays QQ).
    let domain=if f.base.domain==Domain::QQ||g.base.domain==Domain::QQ{Domain::QQ}else{Domain::ZZ};
    let mut p=a.base.clone();p.rep=n.to_rat();p.domain=domain;let mut q=b.base.clone();q.rep=m.to_rat();q.domain=domain;
    Some((common,p.into(),q.into()))
}

pub(crate) fn err(e: PolyError) -> PyErr {
    match &e {
        PolyError::ExactQuotientFailed(message)=>domains::ExactQuotientFailed::new_err(message.clone()),
        PolyError::NotInvertible(message)=>an_rational::NotInvertible::new_err(message.clone()),
        PolyError::Unsupported(s) if s.contains("division by zero") || s.contains("modulo by zero") => PyZeroDivisionError::new_err(s.clone()),
        // The reference raises PolynomialError, a ValueError subclass here.
        PolyError::Generator(message) => Python::with_gil(|py| structural::polynomial_error(py, "PolynomialError", message.clone())),
        PolyError::Reference{kind,message}=>Python::with_gil(|py|reference_error(py,kind,message)),
        _ => PyValueError::new_err(e.to_string()),
    }
}

/// Raises the reference exception class named by an expression-domain error.
pub(crate) fn reference_error(py:Python<'_>,kind:&str,message:&str)->PyErr {
    match kind {
        "ZeroDivisionError"=>return PyZeroDivisionError::new_err(message.to_string()),
        "AssertionError"=>return pyo3::exceptions::PyAssertionError::new_err(message.to_string()),
        "ValueError"=>return pyo3::exceptions::PyValueError::new_err(message.to_string()),
        "TypeError"=>return pyo3::exceptions::PyTypeError::new_err(message.to_string()),
        "NotImplementedError"=>return pyo3::exceptions::PyNotImplementedError::new_err(message.to_string()),
        _=>{}
    }
    let build=||->PyResult<PyErr>{
        let errors=py.import("symox.polys.polyerrors")?;
        let class=errors.getattr(kind)?;
        let value=if kind=="PolynomialDivisionFailed" {
            let (f,g)=message.split_once('\u{1f}').unwrap_or((message,""));
            class.call1((f,g,py.import("symox._core")?.getattr("EX")?))?
        } else {class.call1((message,))?};
        Ok(PyErr::from_value(value))
    };
    build().unwrap_or_else(|e|e)
}

fn nroots_error(py: Python<'_>, error: PolyError) -> PyErr {
    let message = error.to_string();
    if message.starts_with("convergence to root failed") {
        return match py.import("mpmath.libmp.libhyper")
            .and_then(|module| module.getattr("NoConvergence"))
            .and_then(|class| class.call1((message,))) {
                Ok(exception) => PyErr::from_value(exception),
                Err(error) => error,
            };
    }
    err(error)
}

fn gens_of(args: &Bound<'_, PyTuple>) -> PyResult<Vec<Expr>> {
    let mut v = Vec::new();
    for a in args.iter() {
        if let Ok(t) = a.downcast::<PyTuple>() {
            for b in t.iter() {
                v.push(to_expr(&b)?);
            }
        } else if let Ok(l) = a.downcast::<PyList>() {
            for b in l.iter() {
                v.push(to_expr(&b)?);
            }
        } else {
            v.push(to_expr(&a)?);
        }
    }
    Ok(v)
}

fn check_mod(p: u64) -> PyResult<u64> {
    if p < 2 || p >= (1 << 32) || !symox::poly::gf::is_prime(p) {
        return Err(PyValueError::new_err(format!("modulus must be a prime below 2**32, got {p}")));
    }
    Ok(p)
}

fn rat_expr(r: &Rat) -> Expr {
    Expr::num(r.to_number())
}

fn algebraic_extension(extension:Option<&Bound<'_,PyAny>>,expressions:&[Expr],gens:&[Expr],domain:Option<&str>)->PyResult<Option<Arc<NumberField>>> {
    if let Some(d)=domain {
        if extension.is_some(){return Err(PyValueError::new_err("specify either domain or extension"));}
        return NumberField::parse(d).map_err(err);
    }
    let Some(extension)=extension else{return Ok(None)};
    if extension.is_instance_of::<pyo3::types::PyBool>() {
        if !extension.extract::<bool>()?{return Err(PyValueError::new_err("extension=False is not a supported extension option"));}
        return symox::poly::algebraic::infer_in(expressions,gens).map_err(err);
    }
    let generators:Vec<Expr>=if extension.downcast::<PyList>().is_ok()||extension.downcast::<PyTuple>().is_ok(){extension.try_iter()?.map(|e|to_expr(&e?)).collect::<PyResult<_>>()?}else{vec![to_expr(extension)?]};
    if generators.is_empty(){return Err(PyValueError::new_err("an algebraic extension needs a generator"));}
    NumberField::from_generators(&generators).map(Some).map_err(extension_error)
}

/// The reference raises `NotAlgebraic` for a transcendental generator.
fn extension_error(e:symox::poly::convert::PolyError)->PyErr {
    let message=e.to_string();
    if message.ends_with("does not seem to be an algebraic element") {
        return Python::with_gil(|py|structural::polynomial_error(py,"NotAlgebraic",message));
    }
    err(e)
}

/// A univariate polynomial given as a coefficient list (SymPy accepts
/// `roots([1, 0, -1])`).
fn expr_or_coeff_list(obj: &Bound<'_, PyAny>) -> PyResult<(Expr, Vec<Expr>)> {
    if let Ok(p) = obj.downcast::<PyPoly>() {
        let p = &p.get().p;
        return Ok((p.as_expr(), p.true_gens().to_vec()));
    }
    if let Ok(l) = obj.downcast::<PyList>() {
        return coeff_list_in(l, &Expr::symbol("x"));
    }
    Ok((to_expr(obj)?, vec![]))
}

/// The expression of a dense descending coefficient list in the symbol `x`.
fn coeff_list_in(l: &Bound<'_, PyList>, x: &Expr) -> PyResult<(Expr, Vec<Expr>)> {
    let n = l.len();
    let mut terms = Vec::new();
    for (i, c) in l.iter().enumerate() {
        let c = if let Ok(value)=c.extract::<PyRef<'_,domains::PyDomainElement>>(){value.value.as_expr()}else{to_expr(&c)?};
        terms.push(symox::ops::mul(&[c, symox::ops::pow(x, &Expr::int((n - 1 - i) as i64))]));
    }
    Ok((symox::ops::add(&terms), vec![x.clone()]))
}

// ----------------------------------------------------------------------
// Poly
// ----------------------------------------------------------------------

/// A deterministic reduced Groebner basis over QQ, ZZ, or a prime GF(p).
#[pyclass(name = "GroebnerBasis", module = "symox._core", frozen)]
#[derive(Clone)]
pub struct PyGroebnerBasis {
    basis: Vec<PolyValue>,
    generators: Vec<Expr>,
    coefficient_domain: Domain,
    prime: Option<u64>,
    monomial_order: symox::poly::groebner::MonomialOrder,
    return_polys: bool,
    fraction_domain: Option<FractionPoly>,
    algebraic_domain: Option<AlgebraicPoly>,
    /// The reference `EX` coefficient domain.
    expression: bool,
    /// The explicitly passed options, which the reference `GroebnerBasis`
    /// equality compares through `_options`.
    options: String,
}

impl PyGroebnerBasis {
    fn build(
        f: &Bound<'_, PyAny>,
        gens: &Bound<'_, PyTuple>,
        order: &str,
        domain: Option<&Bound<'_,PyAny>>,
        modulus: Option<u64>,
        method: Option<&str>,
        polys: Option<bool>,
        extension: Option<&Bound<'_,PyAny>>,
    ) -> PyResult<Self> {
        let mut basis=Self::build_basis(f,gens,order,domain,modulus,method,polys,extension)?;
        let text=|o:Option<&Bound<'_,PyAny>>|->PyResult<Option<String>>{o.filter(|v|!v.is_none()).map(|v|Ok(v.str()?.to_string())).transpose()};
        basis.options=format!("{:?}",(text(domain)?,modulus,method,polys,text(extension)?));
        Ok(basis)
    }

    #[allow(clippy::too_many_arguments)]
    fn build_basis(
        f: &Bound<'_, PyAny>,
        gens: &Bound<'_, PyTuple>,
        order: &str,
        domain: Option<&Bound<'_,PyAny>>,
        modulus: Option<u64>,
        method: Option<&str>,
        polys: Option<bool>,
        extension: Option<&Bound<'_,PyAny>>,
    ) -> PyResult<Self> {
        let selected=domains::Selection::new(domain,modulus,extension)?;
        let domain=selected.legacy.as_deref();
        let modulus=selected.modulus;
        use symox::poly::convert::{from_exprs, from_exprs_gens};
        use symox::poly::groebner::{buchberger, CoefficientField, MonomialOrder};
        let monomial_order = MonomialOrder::parse(order)
            .ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
        let use_f5b=match method {
            None|Some("buchberger")=>false,
            Some("f5b")=>true,
            Some(other)=>return Err(PyValueError::new_err(format!("'{other}' is not a valid Groebner algorithm (valid are 'buchberger' and 'f5b')"))),
        };
        if domain.is_some_and(|d| d != "ZZ" && d != "QQ" && d != "EX" && d != "ZZ_I" && d != "QQ_I" && !d.starts_with("QQ(") && !d.starts_with("ZZ(") && !d.starts_with("QQ[") && !d.starts_with("ZZ[") && !d.starts_with("QQ<")) {
            return Err(PyNotImplementedError::new_err(
                "Groebner bases support only ZZ, QQ, or modulus=<prime>",
            ));
        }
        if modulus.is_some() && domain.is_some() {
            return Err(PyValueError::new_err("specify either domain or modulus"));
        }
        let mut exprs = Vec::new();
        let mut any_poly = false;
        let mut inferred_gens = Vec::new();
        let mut inferred_prime = None;
        let mut poly_domain = Domain::ZZ;
        let mut inferred_algebraic:Option<Arc<NumberField>>=None;
        let mut inferred_gaussian:Option<Option<bool>>=None;
        let mut inferred_fraction:Option<FractionPoly>=None;
        let mut inferred_expression=false;
        for item in f.try_iter()? {
            let item = item?;
            if let Ok(p) = item.downcast::<PyPoly>() {
                any_poly = true;
                let p = &p.get().p;
                exprs.push(p.as_expr());
                inferred_expression|=p.expression.is_some();
                if let Some(a)=&p.algebraic {
                    inferred_gaussian=Some(match (inferred_gaussian,a.gaussian) {(None|Some(Some(true)),Some(g))=>Some(g),(Some(Some(false)),Some(_))=>Some(false),_=>None});
                    inferred_algebraic=Some(if let Some(previous)=inferred_algebraic {if previous.same_embedding(&a.number_field){previous}else{NumberField::from_generators(&[previous.generator.clone(),a.number_field.generator.clone()]).map_err(err)?}}else{a.number_field.clone()});
                }
                if let Some(field)=&p.field {
                    inferred_fraction=Some(if let Some(current)=inferred_fraction {
                        let combined=domains::unify_native(&current.coefficient_domain,&field.coefficient_domain)?;
                        let symox::poly::domains::ExactDomain::Functions{parameters,..}=&combined else{unreachable!()};
                        let mut prototype=FractionPoly::from_expr(&Expr::zero(),&current.gens,parameters,&combined.to_string()).map_err(err)?;prototype.coefficient_domain=combined;prototype
                    }else{field.clone()});
                }
                inferred_gens = symox::poly::polyobj::unify_gens(&inferred_gens, p.true_gens());
                poly_domain = poly_domain.max(p.domain);
                if let Some(m) = p.modulus {
                    if inferred_prime.is_some_and(|q| q != m) {
                        return Err(PyValueError::new_err("incompatible finite-field domains"));
                    }
                    inferred_prime = Some(m);
                }
            } else {
                exprs.push(to_expr(&item)?);
            }
        }
        // As in SymPy's parallel expression conversion, input Poly objects
        // select the return type and generators; options select the field.
        let prime = modulus.map(check_mod).transpose()?;
        let mut generators = gens_of(gens)?;
        if generators.is_empty() && any_poly {
            generators = inferred_gens;
        }
        let inferred_fraction_display=inferred_fraction.as_ref().map(|field|field.display.clone());
        let domain=domain.or(inferred_fraction_display.as_deref());
        // `EX` is requested, carried by an input Poly or, with explicit
        // generators and no options, selected by `construct_domain`.
        if !generators.is_empty()&&!exprs.is_empty()&&modulus.is_none()&&inferred_prime.is_none()&&extension.is_none() {
            use symox::poly::expression_domain::{ExpressionPoly,joint_selection,guarded};
            let explicit=domain==Some("EX")||(domain.is_none()&&inferred_expression);
            let selection=if explicit {Some(exprs.iter().map(|e|ExpressionPoly::from_expr(e,&generators)).collect::<Result<Vec<_>,_>>().map_err(err)?)}
                else if domain.is_none()&&inferred_algebraic.is_none()&&!any_poly {joint_selection(&exprs,&generators)} else {None};
            if let Some(selection)=selection {
                let sparse:Vec<_>=selection.iter().map(ExpressionPoly::to_sparse).collect();
                let reps=guarded(||if use_f5b {symox::poly::f5b::f5b(&sparse,monomial_order)}else{Ok(symox::poly::fieldalg::buchberger(&sparse,monomial_order))}).map_err(err)?.map_err(err)?;
                let basis=reps.iter().map(|rep|ExpressionPoly::from_sparse(rep,&generators).map(Into::into)).collect::<Result<Vec<PolyValue>,_>>().map_err(err)?;
                return Ok(Self {basis,generators,coefficient_domain:Domain::QQ,prime:None,monomial_order,return_polys:polys.unwrap_or(any_poly),fraction_domain:None,algebraic_domain:None,expression:true,options:String::new()});
            }
        }
        // Without options the reference reads `I` as a coefficient and
        // selects ZZ_I or QQ_I when it is the only nonrational atom.
        if domain.is_none()&&extension.is_none()&&modulus.is_none()&&inferred_algebraic.is_none()&&inferred_prime.is_none()&&!exprs.is_empty(){
            if let Some(selection)=symox::poly::gaussian::default_selection(&exprs,&generators){
                generators=selection.gens;
                inferred_algebraic=Some(symox::poly::gaussian::field());
                inferred_gaussian=Some(Some(selection.integral));
            }
        }
        let algebraic=selected.algebraic_in(extension,&exprs,&generators)?.or(if domain.is_none(){inferred_algebraic}else{None});
        if let Some(number_field)=algebraic {
            if modulus.is_some() || inferred_prime.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
            if generators.is_empty(){generators=symox::poly::algebraic::infer_generators(&exprs,&number_field);}
            let mut prototype=AlgebraicPoly::from_expr(&Expr::zero(),&generators,number_field.clone()).map_err(err)?;
            let inputs=exprs.iter().map(|e|prototype.convert(e).map(|p|p.rep)).collect::<Result<Vec<_>,_>>().map_err(err)?;
            prototype.gaussian=selected.gaussian().or(if domain.is_none()&&extension.is_none(){inferred_gaussian.flatten().map(|g|g&&poly_domain==Domain::ZZ)}else{None});
            if prototype.gaussian==Some(true){if let Some(message)=inputs.iter().find_map(symox::poly::gaussian::integrality_error){return Err(domains::CoercionFailed::new_err(message));}}
            let _scope=AlgebraicScope::new(number_field);
            let reps=if use_f5b {symox::poly::f5b::f5b(&inputs,monomial_order).map_err(err)?}else{symox::poly::fieldalg::buchberger(&inputs,monomial_order)};
            let basis=reps.into_iter().map(|rep|prototype.with_rep(gaussian_ring_basis(&prototype,rep)).into()).collect();
            return Ok(Self {basis,generators,coefficient_domain:Domain::QQ,prime:None,monomial_order,return_polys:polys.unwrap_or(any_poly),fraction_domain:None,algebraic_domain:Some(prototype),expression:false,options:String::new()});
        }
        let mut available=symox::poly::fraction::infer_parameters(&exprs,&[]);
        if let Some(field)=&inferred_fraction {for parameter in &field.parameters {if !available.iter().any(|p|p.as_symbol()==parameter.as_symbol()){available.push(parameter.clone());}}}
        let declared=selected.parameters(&available)?.or_else(||if selected.legacy.is_none()&&selected.native.is_none(){inferred_fraction.as_ref().map(|p|p.parameters.clone())}else{None});
        if generators.is_empty() {
            if let Some(parameters)=&declared {generators=symox::poly::fraction::infer_parameters(&exprs,parameters);}
        }
        for (i, g) in generators.iter().enumerate() {
            if g.is_number() || generators[..i].contains(g) {
                return Err(PyValueError::new_err(
                    "generators must be distinct nonnumeric expressions",
                ));
            }
        }
        let converted = if generators.is_empty() {
            from_exprs(&exprs)
        } else {
            from_exprs_gens(&exprs, &generators)
        }
        .map_err(err)?;
        if converted.ngens == 0 {
            return Err(PyValueError::new_err("groebner() requires generators"));
        }
        if declared.is_some() || converted.ngens != converted.gens.len() {
            if modulus.is_some() {return Err(PyNotImplementedError::new_err("finite rational function coefficient fields are not implemented"));}
            // The reference `construct_domain` composite generators, which
            // may be constants such as `pi` or functions such as `exp(a)`.
            let composite=if declared.is_none()&&domain.is_none(){
                let gens=&converted.gens[..converted.ngens];
                let mut coefficients=Vec::new();
                for e in &exprs {coefficients.extend(symox::poly::expression_domain::dict_from_expr(e,gens).map_err(err)?.into_iter().map(|(_,c)|c));}
                symox::poly::expression_domain::composite_parameters(&coefficients).map_err(err)?
            }else{None};
            let parameters=match &composite {Some((p,_))=>p.clone(),None=>declared.unwrap_or_else(||symox::poly::fraction::infer_parameters(&exprs,&generators))};
            if parameters.is_empty() {return Err(PyNotImplementedError::new_err("algebraic coefficient extensions are not implemented"));}
            let has_fractions=match &composite {Some((_,f))=>*f,None=>converted.gens[converted.ngens..].iter().any(|g|matches!(g.kind(),symox::Kind::Pow(_,e) if e.as_number().is_some_and(|n|n.is_negative())))};
            let display=domain.map(str::to_string).unwrap_or_else(|| {
                let base=if has_fractions || converted.polys.iter().all(|p|p.is_integral()) {"ZZ"}else{"QQ"};
                let names=parameters.iter().map(symox::printer::to_str).collect::<Vec<_>>().join(",");
                if has_fractions {format!("{base}({names})")}else{format!("{base}[{names}]")}
            });
            let mut prototype=FractionPoly::from_expr(&Expr::zero(),&generators,&parameters,&display).map_err(err)?;
            if selected.legacy.is_none()&&selected.native.is_none(){if let Some(inferred)=&inferred_fraction {prototype.coefficient_domain=inferred.coefficient_domain.clone();}}
            selected.apply_fraction(&mut prototype);
            let converted_inputs=exprs.iter().map(|e|prototype.convert(e)).collect::<Result<Vec<_>,_>>().map_err(err)?;
            if display.ends_with(']') && converted_inputs.iter().any(|p|!p.belongs_to_ring()){return Err(PyValueError::new_err("coefficient does not belong to the selected polynomial ring"));}
            let inputs=converted_inputs.into_iter().map(|p|p.rep).collect::<Vec<_>>();
            let _scope=symox::poly::ratfunc::Scope::new(parameters.len());
            let reps=if use_f5b {symox::poly::f5b::f5b(&inputs,monomial_order).map_err(err)?}else{symox::poly::fieldalg::buchberger(&inputs,monomial_order)};
            let basis=reps.into_iter().map(|rep| {
                let p=prototype.with_rep(rep);
                if display.ends_with(']') {p.ring_normalize().into()}else{p.into()}
            }).collect();
            return Ok(Self {basis,generators,coefficient_domain:Domain::QQ,prime:None,monomial_order,return_polys:polys.unwrap_or(any_poly),fraction_domain:Some(prototype),algebraic_domain:None,expression:false,options:String::new()});
        }
        if converted.ngens != converted.gens.len() || converted.float || poly_domain == Domain::RR {
            return Err(PyNotImplementedError::new_err("Groebner bases require exact rational coefficients; coefficient extensions are not implemented"));
        }
        generators = converted.gens;
        if let Some(p) = prime {
            if converted
                .polys
                .iter()
                .any(|f| f.coeffs().iter().any(|c| c.d.rem_u64(p) == 0))
            {
                return Err(PyValueError::new_err(
                    "coefficient denominator is zero in the finite field",
                ));
            }
        }
        let integral = converted.polys.iter().all(|p| p.is_integral());
        let coefficient_domain = match domain {
            Some("ZZ") if !integral => {
                return Err(PyValueError::new_err("coefficients are not integers"))
            }
            Some("QQ") => Domain::QQ,
            Some("ZZ") => Domain::ZZ,
            _ if !integral => Domain::QQ,
            _ => Domain::ZZ,
        };
        let field = CoefficientField { modulus: prime };
        let reps = if use_f5b {
            if let Some(prime)=prime {symox::poly::f5b::f5b_modular(&converted.polys,monomial_order,prime).map_err(err)?}
            else {symox::poly::f5b::f5b(&converted.polys,monomial_order).map_err(err)?}
        } else {buchberger(&converted.polys, monomial_order, field)};
        let basis = reps
            .into_iter()
            .map(|rep| {
                let rep = if coefficient_domain == Domain::ZZ && prime.is_none() {
                    let (_, p) = rep.clear_denoms();
                    p.div_coeff(&p.content()).unwrap().to_rat()
                } else {
                    rep
                };
                Poly {
                    rep,
                    gens: generators.clone(),
                    ngens: generators.len(),
                    domain: coefficient_domain,
                    modulus: prime,
                }.into()
            })
            .collect();
        Ok(Self {
            basis,
            generators,
            coefficient_domain,
            prime,
            monomial_order,
            return_polys: polys.unwrap_or(any_poly),
            fraction_domain: None,
            algebraic_domain: None,
            expression: false,
            options: String::new(),
        })
    }

    fn objects(&self, py: Python<'_>) -> PyResult<Vec<PyObject>> {
        self.basis
            .iter()
            .map(|p| {
                if self.return_polys {
                    poly_obj(py, p.clone())
                } else {
                    Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())
                }
            })
            .collect()
    }

    fn expression_reps(&self)->Vec<symox::poly::MPoly<symox::poly::expression_domain::EX>>{
        self.basis.iter().map(|p|p.expression.as_ref().unwrap().to_sparse()).collect()
    }
    fn expression_input(&self,f:&Bound<'_,PyAny>)->PyResult<symox::poly::MPoly<symox::poly::expression_domain::EX>>{
        let e=if let Ok(p)=f.downcast::<PyPoly>(){p.get().p.as_expr()}else{to_expr(f)?};
        Ok(symox::poly::expression_domain::ExpressionPoly::from_expr(&e,&self.generators).map_err(err)?.to_sparse())
    }
    fn expression_value(&self,rep:&symox::poly::MPoly<symox::poly::expression_domain::EX>)->PyResult<PolyValue>{
        Ok(symox::poly::expression_domain::ExpressionPoly::from_sparse(rep,&self.generators).map_err(err)?.into())
    }

    fn plain_repr(&self) -> String {
        let es: Vec<_> = self
            .basis
            .iter()
            .map(|p| symox::printer::to_str(&p.as_expr()))
            .collect();
        let gs: Vec<_> = self.generators.iter().map(symox::printer::to_str).collect();
        format!(
            "GroebnerBasis([{}], {}, domain='{}', order='{}')",
            es.join(", "),
            gs.join(", "),
            self.domain_descriptor(),
            self.order()
        )
    }

    fn converted_input(&self, f: &Bound<'_, PyAny>) -> PyResult<symox::poly::MPoly<Rat>> {
        let e = if let Ok(p) = f.downcast::<PyPoly>() {
            p.get().p.as_expr()
        } else {
            to_expr(f)?
        };
        let c = symox::poly::convert::from_exprs_gens(&[e], &self.generators).map_err(err)?;
        if c.float || c.ngens != c.gens.len() {
            return Err(PyValueError::new_err(
                "reduction requires exact rational coefficients in the basis generators",
            ));
        }
        let rep = c.polys.into_iter().next().unwrap();
        if self.coefficient_domain == Domain::ZZ && self.prime.is_none() && !rep.is_integral() {
            return Err(PyValueError::new_err("coefficients are not integers in domain ZZ"));
        }
        if let Some(p) = self.prime {
            if rep.coeffs().iter().any(|c| c.d.rem_u64(p) == 0) {
                return Err(PyValueError::new_err(
                    "coefficient denominator is zero in the finite field",
                ));
            }
        }
        Ok(rep)
    }
}

#[pymethods]
impl PyGroebnerBasis {
    #[new]
    #[pyo3(signature = (f, *gens, order="lex", domain=None, modulus=None, method=None, polys=None, extension=None))]
    fn py_new(
        f: &Bound<'_, PyAny>,
        gens: &Bound<'_, PyTuple>,
        order: &str,
        domain: Option<&Bound<'_,PyAny>>,
        modulus: Option<u64>,
        method: Option<&str>,
        polys: Option<bool>,
        extension: Option<&Bound<'_,PyAny>>,
    ) -> PyResult<Self> {
        Self::build(f, gens, order, domain, modulus, method, polys, extension)
    }

    #[getter]
    fn polys(&self) -> Vec<PyPoly> {
        self.basis.iter().cloned().map(|p| PyPoly { p: p.into() }).collect()
    }
    #[getter]
    fn exprs(&self) -> Vec<PyExpr> {
        self.basis.iter().map(|p| w(p.as_expr())).collect()
    }
    #[getter]
    fn gens<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.generators.iter().cloned().map(w))
    }
    #[getter]
    fn order(&self) -> &'static str {
        self.monomial_order.name()
    }
    #[getter]
    fn domain(&self,py:Python<'_>)->PyResult<PyObject>{self.domain_descriptor().object(py)}

    fn domain_descriptor(&self) -> domains::PyDomain {
        use symox::poly::domains::ExactDomain;
        domains::PyDomain {inner:if self.expression {ExactDomain::Expression{raw:false}}else if let Some(field)=&self.algebraic_domain {if let Some(integral)=field.gaussian{ExactDomain::Gaussian{integral}}else{ExactDomain::Algebraic(field.number_field.clone())}}else if let Some(field)=&self.fraction_domain {field.coefficient_domain.clone()}else if let Some(p)=self.prime {ExactDomain::GF(p)}else if self.coefficient_domain==Domain::ZZ {ExactDomain::ZZ}else{ExactDomain::QQ}}
    }

    /// Convert orders through the finite-dimensional quotient algebra.
    fn fglm(&self, order: &str) -> PyResult<Self> {
        use symox::poly::groebner::{fglm, CoefficientField, MonomialOrder};
        let target = MonomialOrder::parse(order)
            .ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
        if target == self.monomial_order { return Ok(self.clone()); }
        if self.expression {
            let reps=self.expression_reps();
            let converted=symox::poly::expression_domain::guarded(||symox::poly::fieldalg::fglm(&reps,self.monomial_order,target,self.generators.len())).map_err(err)?
                .map_err(|error| PyNotImplementedError::new_err(error.to_string()))?;
            let mut result=self.clone();
            result.monomial_order=target;
            result.basis=converted.iter().map(|rep|self.expression_value(rep)).collect::<PyResult<_>>()?;
            return Ok(result);
        }
        if let Some(field) = &self.algebraic_domain {
            let _scope = AlgebraicScope::new(field.number_field.clone());
            let reps: Vec<_> = self.basis.iter().map(|p| p.algebraic.as_ref().unwrap().rep.clone()).collect();
            let converted = symox::poly::fieldalg::fglm(&reps,self.monomial_order,target,self.generators.len())
                .map_err(|error| PyNotImplementedError::new_err(error.to_string()))?;
            let mut result = self.clone();
            result.monomial_order = target;
            result.basis = converted.into_iter().map(|rep| field.with_rep(gaussian_ring_basis(field,rep)).into()).collect();
            return Ok(result);
        }
        if let Some(field) = &self.fraction_domain {
            let _scope = symox::poly::ratfunc::Scope::new(field.parameters.len());
            let reps: Vec<_> = self.basis.iter().map(|p| p.field.as_ref().unwrap().rep.clone()).collect();
            let converted = symox::poly::fieldalg::fglm(&reps,self.monomial_order,target,self.generators.len())
                .map_err(|error| PyNotImplementedError::new_err(error.to_string()))?;
            let mut result = self.clone();
            result.monomial_order = target;
            result.basis = converted.into_iter().map(|rep| {
                let poly = field.with_rep(rep);
                if field.display.ends_with(']') { poly.ring_normalize().into() } else { poly.into() }
            }).collect();
            return Ok(result);
        }
        let reps: Vec<_> = self.basis.iter().map(|p| p.rep.clone()).collect();
        let converted = fglm(&reps,self.monomial_order,target,CoefficientField {modulus:self.prime},self.generators.len())
            .map_err(|e| PyNotImplementedError::new_err(e.to_string()))?;
        let mut result = self.clone();
        result.monomial_order = target;
        result.basis = converted.into_iter().map(|rep| {
            let rep = if self.coefficient_domain == Domain::ZZ && self.prime.is_none() {
                let (_,p) = rep.clear_denoms(); p.div_coeff(&p.content()).unwrap().to_rat()
            } else {rep};
            Poly {rep,gens:self.generators.clone(),ngens:self.generators.len(),domain:self.coefficient_domain,modulus:self.prime}.into()
        }).collect();
        Ok(result)
    }
    #[getter]
    fn is_zero_dimensional(&self) -> bool {
        if self.expression {
            return symox::poly::fieldalg::zero_dimensional(&self.expression_reps(),self.monomial_order,self.generators.len());
        }
        if let Some(field)=&self.algebraic_domain {
            let _scope=AlgebraicScope::new(field.number_field.clone());
            let reps:Vec<_>=self.basis.iter().map(|p|p.algebraic.as_ref().unwrap().rep.clone()).collect();
            return symox::poly::fieldalg::zero_dimensional(&reps,self.monomial_order,self.generators.len());
        }
        if let Some(field)=&self.fraction_domain {
            let reps:Vec<_>=self.basis.iter().map(|p|p.field.as_ref().unwrap().rep.clone()).collect();
            let _scope=symox::poly::ratfunc::Scope::new(field.parameters.len());
            return symox::poly::fieldalg::zero_dimensional(&reps,self.monomial_order,self.generators.len());
        }
        let reps: Vec<_> = self.basis.iter().map(|p| p.rep.clone()).collect();
        symox::poly::groebner::is_zero_dimensional(
            &reps,
            self.monomial_order,
            self.generators.len(),
        )
    }
    fn __len__(&self) -> usize {
        self.basis.len()
    }
    /// `GroebnerBasis.__eq__`: equal bases and options, or the expressions
    /// or polynomials of an iterable.
    fn __eq__(&self,py:Python<'_>,o:&Bound<'_,PyAny>)->PyResult<bool>{
        if let Ok(other)=o.downcast::<PyGroebnerBasis>(){
            let other=other.get();
            return Ok(self.options==other.options&&self.generators==other.generators&&self.monomial_order==other.monomial_order
                &&self.domain_descriptor().to_string()==other.domain_descriptor().to_string()
                &&self.basis.iter().map(|p|p.as_expr()).eq(other.basis.iter().map(|p|p.as_expr())));
        }
        if o.is_instance_of::<pyo3::types::PyString>()||o.is_instance_of::<PyDict>()||o.try_iter().is_err(){return Ok(false);}
        let items=PyList::new(py,o.try_iter()?.collect::<PyResult<Vec<_>>>()?)?;
        if PyList::new(py,self.exprs())?.as_any().eq(&items)?{return Ok(true);}
        PyList::new(py,self.polys())?.as_any().eq(&items)
    }
    fn __ne__(&self,py:Python<'_>,o:&Bound<'_,PyAny>)->PyResult<bool>{Ok(!self.__eq__(py,o)?)}
    fn __hash__(&self)->u64{
        use std::hash::{Hash,Hasher};
        let mut h=std::collections::hash_map::DefaultHasher::new();
        self.plain_repr().hash(&mut h);self.options.hash(&mut h);h.finish()
    }
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyList::new(py, self.objects(py)?)?
            .try_iter()
            .map(|it| it.into_any())
    }
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        Ok(PyList::new(py, self.objects(py)?)?
            .as_any()
            .get_item(key)?
            .unbind())
    }
    /// `_print_GroebnerBasis` prints each basis element with the basis
    /// monomial order; the native printer covers the lexicographic order.
    fn __repr__(slf:&Bound<'_,Self>) -> String {
        let this=slf.get();
        if this.monomial_order.name()!="lex" {
            let printed=slf.py().import("symox.printing.str").and_then(|m|m.getattr("sstr")?.call1((slf,))?.extract::<String>());
            if let Ok(text)=printed {return text;}
        }
        this.plain_repr()
    }
    fn __str__(slf:&Bound<'_,Self>) -> String {Self::__repr__(slf)}

    fn __reduce__(&self,py:Python<'_>)->PyResult<(PyObject,PyObject)> {
        let helper=py.import("symox._core")?.getattr("_rebuild_groebner")?.unbind();
        let parameters=self.fraction_domain.as_ref().map(|f|f.parameters.iter().cloned().map(w).collect::<Vec<_>>());
        let arguments=(self.exprs(),self.generators.iter().cloned().map(w).collect::<Vec<_>>(),parameters,self.domain_descriptor(),self.order(),self.return_polys,Option::<u64>::None,self.options.clone()).into_pyobject(py)?.into_any().unbind();
        Ok((helper,arguments))
    }

    #[pyo3(signature = (f, auto=true))]
    fn reduce(
        &self,
        py: Python<'_>,
        f: &Bound<'_, PyAny>,
        auto: bool,
    ) -> PyResult<(Vec<PyObject>, PyObject)> {
        if self.expression {
            let input=self.expression_input(f)?;
            let reps=self.expression_reps();
            let(mut q,r)=symox::poly::expression_domain::guarded(||symox::poly::fieldalg::reduce(&input,&reps,self.monomial_order)).map_err(err)?;
            if input.is_zero(){q.clear();}
            let object=|rep:&symox::poly::MPoly<symox::poly::expression_domain::EX>|->PyResult<PyObject>{let p=self.expression_value(rep)?;if self.return_polys{poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
            return Ok((q.iter().map(object).collect::<PyResult<_>>()?,object(&r)?));
        }
        if let Some(field)=&self.algebraic_domain {
            let input=field.convert(&to_expr(f)?).map_err(err)?;
            let _scope=AlgebraicScope::new(field.number_field.clone());
            let reps:Vec<_>=self.basis.iter().map(|p|p.algebraic.as_ref().unwrap().rep.clone()).collect();
            if !auto&&field.gaussian==Some(true)&&!symox::poly::gaussian::poly_is_integral(&input.rep){return Err(PyValueError::new_err("coefficients are not Gaussian integers"));}
            let(mut q,r)=symox::poly::fieldalg::reduce(&input.rep,&reps,self.monomial_order);
            if input.rep.is_zero(){q.clear();}
            // The reference reduces over the fraction field and retracts
            // when every quotient and the remainder are integral.
            let integral=q.iter().chain(std::iter::once(&r)).all(symox::poly::gaussian::poly_is_integral);
            let target=if integral{field.clone()}else{field.to_field()};
            let object=|rep|->PyResult<PyObject>{let p=target.with_rep(rep);if self.return_polys{poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
            return Ok((q.into_iter().map(object).collect::<PyResult<_>>()?,object(r)?));
        }
        if let Some(field)=&self.fraction_domain {
                        let input=field.convert(&to_expr(f)?).map_err(err)?;
            if field.display.ends_with(']') && !input.belongs_to_ring(){return Err(PyValueError::new_err("input coefficients do not belong to the basis ring"));}
            let reps:Vec<_>=self.basis.iter().map(|p|p.field.as_ref().unwrap().rep.clone()).collect();
            let _scope=symox::poly::ratfunc::Scope::new(field.parameters.len());
            let(mut q,r)=symox::poly::fieldalg::reduce_filtered(&input.rep,&reps,self.monomial_order,|c|auto || !field.display.ends_with(']') || symox::poly::fraction::coefficient_in_ring(c,field.display.starts_with("ZZ")));
            if input.rep.is_zero(){q.clear();}
            let promote=q.iter().chain(std::iter::once(&r)).any(|p|p.coeffs().iter().any(|c|!symox::poly::fraction::coefficient_in_ring(c,field.display.starts_with("ZZ"))));
            let object=|rep|->PyResult<PyObject>{let mut p=field.with_rep(rep);if promote {p.promote_field();}if self.return_polys {poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
            return Ok((q.into_iter().map(object).collect::<PyResult<_>>()?,object(r)?));
        }
        use symox::poly::groebner::{reduce_with_options, CoefficientField};
        let input = self.converted_input(f)?;
        if !auto
            && self.coefficient_domain == Domain::ZZ
            && self.prime.is_none()
            && !input.is_integral()
        {
            return Err(PyValueError::new_err("coefficients are not integers"));
        }
        let reps: Vec<_> = self.basis.iter().map(|p| p.rep.clone()).collect();
        let (mut q, r) = reduce_with_options(
            &input,
            &reps,
            self.monomial_order,
            CoefficientField {
                modulus: self.prime,
            },
            !auto && self.coefficient_domain == Domain::ZZ && self.prime.is_none(),
        );
        if input.is_zero() { q.clear(); }
        let domain = if self.coefficient_domain == Domain::ZZ
            && self.prime.is_none()
            && (!r.is_integral() || q.iter().any(|p| !p.is_integral()))
        {
            Domain::QQ
        } else {
            self.coefficient_domain
        };
        let object = |rep| -> PyResult<PyObject> {
            let p = Poly {
                rep,
                gens: self.generators.clone(),
                ngens: self.generators.len(),
                domain,
                modulus: self.prime,
            };
            if self.return_polys {
                poly_obj(py, p)
            } else {
                Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())
            }
        };
        Ok((
            q.into_iter().map(object).collect::<PyResult<_>>()?,
            object(r)?,
        ))
    }

    fn contains(&self, f: &Bound<'_, PyAny>) -> PyResult<bool> {
        if self.expression {
            let input=self.expression_input(f)?;
            let reps=self.expression_reps();
            return Ok(symox::poly::expression_domain::guarded(||symox::poly::fieldalg::reduce(&input,&reps,self.monomial_order).1.is_zero()).map_err(err)?);
        }
        if let Some(field)=&self.algebraic_domain {
            let input=field.convert(&to_expr(f)?).map_err(err)?;
            let _scope=AlgebraicScope::new(field.number_field.clone());
            let reps:Vec<_>=self.basis.iter().map(|p|p.algebraic.as_ref().unwrap().rep.clone()).collect();
            return Ok(symox::poly::fieldalg::reduce(&input.rep,&reps,self.monomial_order).1.is_zero());
        }
        if let Some(field)=&self.fraction_domain {
            let input=field.convert(&to_expr(f)?).map_err(err)?;
            if field.display.ends_with(']') && !input.belongs_to_ring(){return Err(PyValueError::new_err("input coefficients do not belong to the basis ring"));}
            let reps:Vec<_>=self.basis.iter().map(|p|p.field.as_ref().unwrap().rep.clone()).collect();
            let _scope=symox::poly::ratfunc::Scope::new(field.parameters.len());
            return Ok(symox::poly::fieldalg::reduce(&input.rep,&reps,self.monomial_order).1.is_zero());
        }
        use symox::poly::groebner::{reduce, CoefficientField};
        let input = self.converted_input(f)?;
        let reps: Vec<_> = self.basis.iter().map(|p| p.rep.clone()).collect();
        Ok(reduce(
            &input,
            &reps,
            self.monomial_order,
            CoefficientField {
                modulus: self.prime,
            },
        )
        .1
        .is_zero())
    }
}

/// A reduced basis over `ZZ_I` is the field basis with cleared
/// denominators, as the reference `groebner` returns for a ring domain.
fn gaussian_ring_basis(field:&AlgebraicPoly,rep:symox::poly::mpoly::MPoly<symox::poly::algebraic::AN>)->symox::poly::mpoly::MPoly<symox::poly::algebraic::AN> {
    if field.gaussian==Some(true){symox::poly::gaussian::clear_denoms(&rep).1}else{rep}
}

/// The `gaussian` option of functions with named keyword parameters.
fn named_gaussian_domain<'py>(py:Python<'py>,gaussian:Option<&Bound<'py,PyAny>>,named:&[(&str,bool)])->PyResult<Option<Bound<'py,PyAny>>> {
    let Some(value)=gaussian.filter(|g|!g.is_none()) else{return Ok(None)};
    let kw=PyDict::new(py);
    kw.set_item("gaussian",value)?;
    for (name,present) in named {if *present{kw.set_item(*name,true)?;}}
    division::gaussian_option(&kw)
}

#[pyfunction]
#[pyo3(signature = (f, *gens, order="lex", domain=None, modulus=None, method=None, polys=None, extension=None, gaussian=None))]
fn groebner<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    order: &str,
    domain: Option<&Bound<'py,PyAny>>,
    modulus: Option<u64>,
    method: Option<&str>,
    polys: Option<bool>,
    extension: Option<&Bound<'py,PyAny>>,
    gaussian: Option<&Bound<'py,PyAny>>,
) -> PyResult<PyGroebnerBasis> {
    let selected=named_gaussian_domain(f.py(),gaussian,&[("domain",domain.is_some_and(|d|!d.is_none())),("extension",extension.is_some_and(|e|!e.is_none())),("modulus",modulus.is_some())])?;
    PyGroebnerBasis::build(f, gens, order, selected.as_ref().or(domain), modulus, method, polys, extension)
}

#[pyfunction]
#[pyo3(signature = (F, *gens, order="lex", domain=None, modulus=None, method=None, extension=None))]
fn is_zero_dimensional(
    F: &Bound<'_, PyAny>,
    gens: &Bound<'_, PyTuple>,
    order: &str,
    domain: Option<&Bound<'_,PyAny>>,
    modulus: Option<u64>,
    method: Option<&str>,
    extension: Option<&Bound<'_,PyAny>>,
) -> PyResult<bool> {
    Ok(
        PyGroebnerBasis::build(F, gens, order, domain, modulus, method, None, extension)?
            .is_zero_dimensional(),
    )
}

#[pyfunction]
fn _solve_poly_system(f: &Bound<'_,PyAny>, symbols: &Bound<'_,PyAny>) -> PyResult<Option<Vec<Vec<PyExpr>>>> {
    let equations = f.try_iter()?.map(|e| to_expr(&e?)).collect::<PyResult<Vec<_>>>()?;
    let variables = symbols.try_iter()?.map(|e| to_expr(&e?)).collect::<PyResult<Vec<_>>>()?;
    let result = symox::poly::systems::solve_poly_system(&equations,&variables).map_err(|e| match e {
        PolyError::Unsupported(message) => PyNotImplementedError::new_err(message),
        other => err(other),
    })?;
    Ok(result.map(|solutions| solutions.into_iter().map(|s| s.into_iter().map(w).collect()).collect()))
}

/// A polynomial with explicit generators (SymPy's `Poly`).
#[pyclass(name = "Poly", module = "symox._core", frozen, subclass)]
#[derive(Clone)]
pub struct PyPoly {
    pub p: PolyValue,
}

fn set_domain(mut p: Poly, domain: Option<&str>) -> PyResult<Poly> {
    match domain {
        None => Ok(p),
        Some("ZZ") => {
            if !p.rep.is_integral() {
                return Err(PyValueError::new_err("coefficients are not integers"));
            }
            p.domain = Domain::ZZ;
            Ok(p)
        }
        Some("QQ") => {
            p.domain = Domain::QQ;
            Ok(p)
        }
        Some("RR") => {
            p.domain = Domain::RR;
            Ok(p)
        }
        Some(d) => Err(PyValueError::new_err(format!("unsupported domain {d}"))),
    }
}

impl PyPoly {
    fn build(expr: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, domain: Option<&str>) -> PyResult<Self> {
        use symox::poly::fraction::{parse_domain,infer_parameters,domain_name};
        if let Some(field)=algebraic_extension(None,&[],&[],domain)? {
            let(e,g2)=expr_or_coeff_list(expr)?;let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}
            if generators.is_empty(){generators=symox::poly::algebraic::infer_generators(std::slice::from_ref(&e),&field);}
            return Ok(Self {p:AlgebraicPoly::from_expr(&e,&generators,field).map_err(err)?.into()});
        }
        let mut generators = gens_of(gens)?;
        if let Ok(q) = expr.downcast::<PyPoly>() {
            if generators.is_empty() && domain.is_none() {return Ok(q.get().clone());}
        }
        let (e,g2) = expr_or_coeff_list(expr)?;
        if generators.is_empty() {generators=g2;}
        let mut available: Vec<_> = symox::traverse::free_symbols(&e).into_iter().collect();
        for g in &generators {if !available.contains(g) {available.push(g.clone());}}
        if let Some(d)=domain {
            if let Some(parameters)=parse_domain(d,&available).map_err(err)? {
                if generators.is_empty() {generators=infer_parameters(std::slice::from_ref(&e),&parameters);}
                let p=FractionPoly::from_expr(&e,&generators,&parameters,d).map_err(err)?;
                if d.ends_with(']') && !p.belongs_to_ring(){return Err(PyValueError::new_err("coefficient does not belong to the selected polynomial ring"));}
                return Ok(Self {p:p.into()});
            }
        }
        // The reference reads `I` as a coefficient and selects ZZ_I or QQ_I
        // when it is the only nonrational coefficient atom.
        if domain.is_none(){
            if let Some(selection)=symox::poly::gaussian::default_selection(std::slice::from_ref(&e),&generators){
                return Self::gaussian_poly(&e,selection.gens,selection.integral);
            }
            if generators.is_empty()&&symox::traverse::has(&e,&Expr::i())&&symox::poly::convert::from_exprs(std::slice::from_ref(&e)).is_ok_and(|c|!c.float&&c.gens==[Expr::i()]){
                return Err(structural::polynomial_error(expr.py(),"GeneratorsNeeded","Cannot initialize from 'dict' without generators".into()));
            }
        }
        let legacy=Poly::from_expr(&e,&generators);
        // With explicit generators the reference `construct_domain` selects
        // EX for algebraic-number or related coefficient generators.
        if domain.is_none()&&!generators.is_empty()&&legacy.as_ref().map_or(true,|p|!p.domain_gens().is_empty()){
            if let Some(p)=Self::expression_selection(&e,&generators)?{return Ok(p);}
            // `_construct_composite` builds a fraction field when a
            // coefficient denominator involves a generator, as for `1/pi`.
            if let Ok(terms)=symox::poly::expression_domain::dict_from_expr(&e,&generators){
                let coefficients:Vec<Expr>=terms.into_iter().map(|(_,c)|c).collect();
                if let Ok(Some(parameters))=symox::poly::expression_domain::composite_fraction_parameters(&coefficients){
                    let display=domain_name(&parameters).replacen("QQ","ZZ",1);
                    return Ok(Self {p:FractionPoly::from_expr(&e,&generators,&parameters,&display).map_err(err)?.into()});
                }
            }
        }
        let p=match legacy{
            Ok(p)=>p,
            // The reference's constructor reports a constant without
            // generators as `GeneratorsNeeded` from its dict representation.
            Err(PolyError::Failed(m)) if generators.is_empty() && m.starts_with("can't initialize from")=>{
                let py=expr.py();
                return Err(structural::polynomial_error(py,"GeneratorsNeeded","Cannot initialize from 'dict' without generators".into()));
            }
            Err(e)=>return Err(err(e)),
        };
        let fractions=p.domain_gens().iter().any(|g| matches!(g.kind(),symox::Kind::Pow(_,e) if e.as_number().is_some_and(|n|n.is_negative())) && !symox::traverse::free_symbols(g).is_empty());
        if domain.is_none() && fractions {
            let parameters=infer_parameters(std::slice::from_ref(&e),p.true_gens());
            let display=domain_name(&parameters).replacen("QQ","ZZ",1);
            return Ok(Self {p:FractionPoly::from_expr(&e,p.true_gens(),&parameters,&display).map_err(err)?.into()});
        }
        if domain.is_some()&&!p.domain_gens().is_empty(){return Err(domains::CoercionFailed::new_err("coefficient does not belong to the selected ground domain"));}
        Ok(Self {p:domains::normalize_poly(set_domain(p,domain)?.into())?})
    }

    /// The EX polynomial of `e` when the reference domain construction
    /// selects EX for its coefficients in `generators`.
    fn expression_selection(e:&Expr,generators:&[Expr])->PyResult<Option<Self>>{
        use symox::poly::expression_domain::{ExpressionPoly,construct_selects_expression};
        let Ok(p)=ExpressionPoly::from_expr(e,generators) else{return Ok(None)};
        let coefficients:Vec<Expr>=p.terms().into_iter().map(|(_,c)|c).collect();
        match construct_selects_expression(&coefficients){
            Ok(true)=>Ok(Some(Self{p:p.into()})),
            _=>Ok(None),
        }
    }

    /// The reference `ZZ_I`/`QQ_I` polynomial of `e`.
    fn gaussian_poly(e:&Expr,mut generators:Vec<Expr>,integral:bool)->PyResult<Self>{
        if generators.is_empty(){generators=symox::poly::algebraic::infer_generators(std::slice::from_ref(e),&symox::poly::gaussian::field());}
        let d=domains::PyDomain{inner:symox::poly::domains::ExactDomain::Gaussian{integral}};
        Ok(Self{p:d.poly(e,&generators)?})
    }

    fn other(&self, o: &Bound<'_, PyAny>) -> Option<PolyValue> {
        if let Ok(q)=o.downcast::<PyPoly>() {
            let q=&q.get().p;
            if self.p.expression.is_some()||q.expression.is_some() {
                self.p.checked_expression_pair(q).ok()??;
                return Some(q.clone());
            }
            if let (Some(a),Some(b))=(&self.p.field,&q.field){domains::unify_native(&a.coefficient_domain,&b.coefficient_domain).ok()?;}
            let gens=symox::poly::polyobj::unify_gens(self.p.true_gens(),q.true_gens());
            if let Some(field)=self.p.algebraic.as_ref().or(q.algebraic.as_ref()) {
                if self.p.field.is_some()||q.field.is_some()||self.p.modulus.is_some()||q.modulus.is_some(){return None;}
                let _=field;
                self.p.checked_algebraic_pair(q).ok()??;
                return Some(q.clone());
            }
            for f in [self.p.field.as_ref(),q.field.as_ref()].into_iter().flatten() {
                if gens.iter().any(|g|f.parameters.contains(g)) || q.modulus.is_some() || self.p.modulus.is_some(){return None;}
            }
            if self.p.field.is_some()||q.field.is_some() {
                self.p.checked_fraction_pair(q).ok()??;
            }
            return Some(q.clone());
        }
        if self.p.expression.is_some() {
            return symox::poly::expression_domain::ExpressionPoly::from_expr(&to_expr(o).ok()?,self.p.true_gens()).ok().map(Into::into);
        }
        if let Some(field)=&self.p.algebraic {return field.convert(&to_expr(o).ok()?).ok().map(Into::into);}
        if let Some(field)=&self.p.field {
            return field.convert(&to_expr(o).ok()?).ok().map(Into::into);
        }
        let e=to_expr(o).ok()?;
        let legacy=Poly::from_expr(&e,self.p.true_gens());
        if legacy.as_ref().map_or(true,|q|!q.domain_gens().is_empty()) {
            if let Ok(Some(q))=Self::expression_selection(&e,self.p.true_gens()){return Some(q.p);}
        }
        legacy.ok().map(|mut q| {
            if q.domain < self.p.domain {q.domain=self.p.domain;}
            q.into()
        })
    }

    fn ex(&self)->Option<&ExpressionPoly>{self.p.expression.as_ref()}
    /// The unified `EX` pair, when either side is over `EX`.
    fn ex_pair(&self,o:&Bound<'_,PyAny>)->PyResult<Option<(ExpressionPoly,ExpressionPoly)>>{
        let other_ex=o.downcast::<PyPoly>().is_ok_and(|q|q.get().p.expression.is_some());
        if self.p.expression.is_none()&&!other_ex{return Ok(None);}
        if let Ok(q)=o.downcast::<PyPoly>(){return self.p.checked_expression_pair(&q.get().p).map_err(err);}
        let q=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial"))?;
        self.p.checked_expression_pair(&q).map_err(err)
    }
    fn binop(&self, o: &Bound<'_, PyAny>, f: impl Fn(&PolyValue, &PolyValue) -> Result<PolyValue,PolyError>, g: impl Fn(&Expr, &Expr) -> Expr) -> PyResult<PyObject> {
        let py = o.py();
        match self.other(o) {
            Some(q) => Ok(PyPoly { p: f(&self.p, &q).map_err(err)? }.into_pyobject(py)?.into_any().unbind()),
            None => {
                if let Ok(q)=o.downcast::<PyPoly>() {
                    self.p.checked_expression_pair(&q.get().p).map_err(err)?;
                    if self.p.algebraic.is_some() || q.get().p.algebraic.is_some() || self.p.field.is_some() || q.get().p.field.is_some() {
                        return Err(PyTypeError::new_err("incompatible coefficient fields; construct a supported joint extension explicitly"));
                    }
                }
                let e = to_expr(o)?;
                Ok(w(g(&self.p.as_expr(), &e)).into_pyobject(py)?.into_any().unbind())
            }
        }
    }

    fn gen_index(&self, g: Option<&Bound<'_, PyAny>>) -> PyResult<usize> {
        match g {
            None => Ok(0),
            Some(g) => {
                if let Ok(i) = g.extract::<i64>() {
                    let count = self.p.true_gens().len() as i64;
                    if i < -count || i >= count {
                        return Err(PyValueError::new_err("polynomial generator index is out of range"));
                    }
                    return Ok(if i < 0 { (count+i) as usize } else { i as usize });
                }
                let e = to_expr(g)?;
                self.p
                    .gen_index(&e)
                    .ok_or_else(|| PyValueError::new_err(format!("{} is not a generator", symox::printer::to_str(&e))))
            }
        }
    }
}

fn poly_obj(py: Python<'_>, p: impl Into<PolyValue>) -> PyResult<PyObject> {
    Ok(PyPoly { p: domains::normalize_poly(p.into())? }.into_pyobject(py)?.into_any().unbind())
}

fn degree_obj(py: Python<'_>, d: Option<u32>) -> PyResult<PyObject> {
    match d {
        Some(k) => Ok(k.into_pyobject(py)?.into_any().unbind()),
        None => Ok(w(Expr::neg_oo()).into_pyobject(py)?.into_any().unbind()),
    }
}

#[pymethods]
impl PyPoly {
    #[new]
    #[pyo3(signature = (rep, *gens, domain=None, modulus=None, extension=None, **_kwargs))]
    fn py_new(rep: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, domain: Option<&Bound<'_,PyAny>>, modulus: Option<u64>, extension:Option<&Bound<'_,PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let normalized=rep.py().import("symox._poly_constructor")?.getattr("construct")?
            .call1((rep,gens,domain,modulus,extension,_kwargs))?;
        Ok(normalized.extract::<PyRef<'_,Self>>()?.clone())
    }

    /// Dense exact ground coefficients need no expression round trip.
    #[staticmethod]
    fn _from_dense_rational(coeffs:Vec<PyExpr>,gen:PyExpr,domain:Option<&str>)->PyResult<Self>{
        use symox::poly::mpoly::MPoly;
        if !matches!(domain,None|Some("ZZ"|"QQ")){
            return Err(PyNotImplementedError::new_err("dense rational construction requires ZZ or QQ"));
        }
        let coefficients=coeffs.iter().map(|coefficient|
            coefficient.e.as_number().and_then(Rat::from_number)
                .ok_or_else(||PyTypeError::new_err("expected exact rational coefficients")))
            .collect::<PyResult<Vec<_>>>()?;
        let integral=coefficients.iter().all(|c|c.d.is_one());
        if domain==Some("ZZ")&&!integral{
            return Err(domains::CoercionFailed::new_err("coefficient does not belong to ZZ"));
        }
        let selected=if domain==Some("QQ")||!integral{Domain::QQ}else{Domain::ZZ};
        let length=coefficients.len();
        let terms=coefficients.into_iter().enumerate().map(|(i,c)|{
            let exponent=u32::try_from(length-1-i).map_err(|_|PyValueError::new_err("polynomial degree limit exceeded"))?;
            Ok((vec![exponent],c))
        }).collect::<PyResult<Vec<_>>>()?;
        Ok(Self{p:Poly::new(MPoly::from_terms(1,terms),vec![gen.e],1,selected).into()})
    }

    /// Backend constructor used after public option/input normalization.
    #[staticmethod]
    #[pyo3(signature=(rep,gens,domain=None,modulus=None,extension=None))]
    fn _from_normalized(rep:&Bound<'_,PyAny>,gens:&Bound<'_,PyTuple>,domain:Option<&Bound<'_,PyAny>>,modulus:Option<u64>,extension:Option<&Bound<'_,PyAny>>)->PyResult<Self>{
        let expr=rep;
        if expr.getattr("is_Matrix").is_ok_and(|flag|flag.is_truthy().unwrap_or(false)) {
            return Err(PyTypeError::new_err("a noncommutative matrix cannot be a scalar polynomial"));
        }
        // As in SymPy, a Poly argument keeps its domain when no domain option is
        // given and the generators are omitted or only reordered.
        if domain.is_none()&&modulus.is_none()&&extension.is_none(){if let Ok(src)=expr.downcast::<PyPoly>(){
            let src=src.get();let wanted=gens_of(gens)?;
            let same_set=wanted.len()==src.p.true_gens().len()&&wanted.iter().all(|g|src.p.true_gens().contains(g));
            if wanted.is_empty(){return Ok(Self{p:src.p.clone()});}
            if same_set{
                if let Some(m)=src.p.modulus{return Self::_from_normalized(expr,gens,None,Some(m),None);}
                if let Ok(d)=src.domain_descriptor(){let d=d.into_pyobject(expr.py())?.into_any();return Self::_from_normalized(expr,gens,Some(&d),None,None);}
            }
        }}
        // A coefficient list is a dense polynomial in the single given
        // generator; expression generators are attached after construction.
        if let (Ok(l),1)=(expr.downcast::<PyList>(),gens.len()){
            let gen=to_expr(&gens.get_item(0)?)?;
            if !gen.is_symbol(){
                let slot=Expr::dummy("x",Default::default());
                let slots=PyTuple::new(expr.py(),[w(slot)])?;
                let built=Self::_from_normalized(l.as_any(),&slots,domain,modulus,extension)?;
                match symox::poly::structure::rename_gens(&built.p,&[gen]){
                    Ok(p)=>return Ok(Self{p}),
                    // Structured coefficient domains keep symbol generators;
                    // read the expression in the composite generator instead.
                    Err(PolyError::Unsupported(_))=>{},
                    Err(e)=>return Err(err(e)),
                }
            }
        }
        if let (Ok(l),1)=(expr.downcast::<PyList>(),gens.len()){
            let (e,_)=coeff_list_in(l,&to_expr(&gens.get_item(0)?)?)?;
            let e=w(e).into_pyobject(expr.py())?.into_any();
            return Self::_from_normalized(&e,gens,domain,modulus,extension);
        }
        let selected=domains::Selection::new(domain,modulus,extension)?;
        let(e,g2)=expr_or_coeff_list(expr)?;
        if let Some(d)=&selected.native {let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}if generators.is_empty()&&!matches!(d.inner,symox::poly::domains::ExactDomain::Expression{..}){if let symox::poly::domains::ExactDomain::Algebraic(field)=&d.inner{generators=symox::poly::algebraic::infer_generators(std::slice::from_ref(&e),field);}else{let excluded=if let symox::poly::domains::ExactDomain::Functions{parameters,..}=&d.inner{parameters.as_slice()}else{&[]};generators=symox::poly::fraction::infer_parameters(std::slice::from_ref(&e),excluded);}}
            return Ok(Self{p:d.poly(&e,&generators)?});}
        let domain=selected.legacy.as_deref();
        let modulus=selected.modulus;
        let mut explicit=gens_of(gens)?;if explicit.is_empty(){explicit=g2;}
        // `extension=True` over only rational and Gaussian rational
        // coefficients selects the reference ZZ_I or QQ_I, not QQ<I>.
        if domain.is_none()&&extension.is_some_and(|x|x.is_instance_of::<pyo3::types::PyBool>()&&x.extract::<bool>().unwrap_or(false)){
            if let symox::poly::algebraic::ExtensionSelection::Gaussian{integral}=symox::poly::algebraic::infer_selection_in(std::slice::from_ref(&e),&explicit).map_err(err)?{
                return Self::gaussian_poly(&e,explicit,integral);
            }
        }
        let mut r=if let Some(field)=algebraic_extension(extension,std::slice::from_ref(&e),&explicit,domain)? {
            let mut generators=explicit;if generators.is_empty(){generators=symox::poly::algebraic::infer_generators(std::slice::from_ref(&e),&field);}
            Self {p:AlgebraicPoly::from_expr(&e,&generators,field).map_err(err)?.into()}
        }else{Self::build(expr,gens,domain)?};
        if let Some(m) = modulus {
            if r.p.field.is_some() || r.p.algebraic.is_some() {return Err(PyNotImplementedError::new_err("finite fraction fields are not implemented"));}
            r.p = r.p.set_modulus(check_mod(m)?).into();
        }
        Ok(r)
    }

    #[getter]
    fn modulus(&self) -> Option<u64> {
        self.p.modulus
    }

    fn __repr__(&self) -> String {
        self.p.to_sympy_string()
    }

    fn __str__(&self) -> String {
        self.p.to_sympy_string()
    }

    fn __reduce__(&self,py:Python<'_>)->PyResult<(PyObject,PyObject)> {
        let helper=py.import("symox._core")?.getattr("_rebuild_poly")?.unbind();
        let parameters=self.p.field.as_ref().map(|f|f.parameters.iter().cloned().map(w).collect::<Vec<_>>());
        let domain=match self.domain_descriptor(){Ok(d)=>d.into_pyobject(py)?.into_any().unbind(),Err(_)=>self.p.domain_str().into_pyobject(py)?.into_any().unbind()};
        let arguments=(w(self.p.as_expr()),self.p.true_gens().iter().cloned().map(w).collect::<Vec<_>>(),parameters,domain,self.p.modulus).into_pyobject(py)?.into_any().unbind();
        Ok((helper,arguments))
    }

    fn _symox_(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    fn as_expr(&self) -> PyExpr {
        w(symox::poly::structure::sympy_expr(&self.p))
    }

    #[getter]
    fn expr(&self) -> PyExpr {
        w(symox::poly::structure::sympy_expr(&self.p))
    }

    /// The same terms over other generators (SymPy's ``Poly.new(rep, *gens)``).
    #[pyo3(signature=(*gens))]
    fn _with_gens(&self,gens:&Bound<'_,PyTuple>)->PyResult<Self>{Ok(Self{p:symox::poly::structure::rename_gens(&self.p,&gens_of(gens)?).map_err(err)?})}

    #[getter]
    fn gens<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.p.true_gens().iter().map(|g| w(g.clone())))
    }

    #[getter]
    fn gen(&self) -> PyExpr {
        w(self.p.true_gens()[0].clone())
    }

    #[getter]
    fn domain(&self,py:Python<'_>)->PyResult<PyObject>{self.domain_descriptor()?.object(py)}

    fn domain_descriptor(&self) -> PyResult<domains::PyDomain> { domains::PyDomain::of(&self.p) }

    fn set_domain(&self,domain:&Bound<'_,PyAny>)->PyResult<Self>{
        let selected=domains::Selection::new(Some(domain),None,None)?;
        if let (Some(f),Some(d))=(self.ex(),&selected.native){
            use symox::poly::domains::ExactDomain as D;
            let fail=|message:String|err(PolyError::Reference{kind:"CoercionFailed",message});
            // Coefficients convert one by one from the leading term, through
            // the target's `from_sympy`.
            match &d.inner {
                D::ZZ=>for (_,c) in f.terms(){
                    if !c.as_number().is_some_and(|n|n.is_integer()||n.is_float()){return Err(fail(format!("expected an integer, got {c}")));}
                },
                D::QQ=>for (_,c) in f.terms(){
                    if c.as_number().is_none(){return Err(fail(format!("expected `Rational` object, got {c}")));}
                },
                D::RR=>{
                    let mut terms=Vec::new();
                    for (m,c) in f.terms(){
                        let v=symox::evalf_mp::evalf(&c,53);
                        if v.as_number().is_none(){return Err(fail(format!("expected real number, got {c}")));}
                        terms.push((m,v));
                    }
                    let e=symox::poly::expression_domain::ExpressionPoly::from_terms(terms,&f.gens).map_err(err)?.as_expr();
                    return Ok(Self{p:d.poly(&e,self.p.true_gens())?});
                }
                _=>{}
            }
        }
        if let Some(d)=selected.native{return Ok(Self{p:d.poly(&self.p.as_expr(),self.p.true_gens())?});}
        let py=domain.py();let gens=PyTuple::new(py,self.p.true_gens().iter().cloned().map(w))?;
        let expr=w(self.p.as_expr()).into_pyobject(py)?.into_any();
        Self::_from_normalized(&expr,&gens,Some(domain),None,None)
    }
    fn to_field(&self)->PyResult<Self>{let d=self.domain_descriptor()?.get_field();Ok(Self{p:d.poly(&self.p.as_expr(),self.p.true_gens())?})}
    fn to_ring(&self)->PyResult<Self>{let d=self.domain_descriptor()?.get_ring()?;Ok(Self{p:d.poly(&self.p.as_expr(),self.p.true_gens())?})}

    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PySet>> {
        let fs = symox::traverse::free_symbols(&self.p.as_expr());
        let mut v: Vec<PyExpr> = Vec::new();
        // SymPy: generators count only when some monomial uses them, plus
        // the symbols of a composite domain (free_symbols_in_domain).
        let used = self.p.true_gens().iter().enumerate().filter(|(i, _)| self.p.degree(*i).is_some_and(|d| d > 0)).map(|(_, g)| g);
        for g in used.chain(self.p.domain_gens().iter()) {
            for s in symox::traverse::free_symbols(g) {
                if !v.iter().any(|x| x.e == s) {
                    v.push(w(s));
                }
            }
        }
        for s in fs {
            if !v.iter().any(|x| x.e == s) {
                v.push(w(s));
            }
        }
        pyo3::types::PySet::new(py, v)
    }

    #[pyo3(signature = (gen=None))]
    fn degree(&self, py: Python<'_>, gen: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
        let i = self.gen_index(gen)?;
        degree_obj(py, self.p.degree(i))
    }

    fn degree_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        // The zero polynomial has degree -oo in every generator, as in SymPy.
        let v = (0..self.p.ngens).map(|i| degree_obj(py, self.p.degree(i))).collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, v)
    }

    fn total_degree(&self, py: Python<'_>) -> PyResult<PyObject> {
        degree_obj(py, Some(self.p.total_degree().unwrap_or(0)))
    }

    fn coeffs(&self) -> Vec<PyExpr> {
        self.p.coeffs().into_iter().map(w).collect()
    }

    fn all_coeffs(&self,py:Python<'_>) -> PyResult<Vec<PyExpr>> {
        if self.p.ngens!=1{return Err(structural::polynomial_error(py,"PolynomialError","multivariate polynomials not supported".into()));}
        Ok(self.p.all_coeffs().map_err(err)?.into_iter().map(w).collect())
    }

    fn deflate<'py>(&self,py:Python<'py>)->PyResult<(Bound<'py,PyTuple>,Self)>{let(j,p)=symox::poly::structure::deflate(&self.p).map_err(err)?;Ok((PyTuple::new(py,j)?,Self{p}))}
    fn terms_gcd<'py>(&self,py:Python<'py>)->PyResult<(Bound<'py,PyTuple>,Self)>{let(j,p)=symox::poly::structure::terms_gcd(&self.p).map_err(err)?;Ok((PyTuple::new(py,j)?,Self{p}))}
    fn exclude(&self)->PyResult<Self>{Ok(Self{p:symox::poly::structure::exclude(&self.p).map_err(err)?})}
    #[pyo3(signature=(*gens,**kwargs))]
    fn reorder(&self,gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Self>{Ok(Self{p:structural::reordered(&self.p,gens,kwargs)?})}
    fn homogeneous_order(&self,py:Python<'_>)->PyResult<PyObject>{use symox::poly::structure::HomogeneousOrder;match symox::poly::structure::homogeneous_order(&self.p).map_err(err)?{
        HomogeneousOrder::Zero=>Ok(w(Expr::neg_oo()).into_pyobject(py)?.into_any().unbind()),HomogeneousOrder::Degree(n)=>Ok(n.into_pyobject(py)?.into_any().unbind()),HomogeneousOrder::Inhomogeneous=>Ok(py.None()),
    }}
    #[pyo3(signature=(*gens))]
    fn has_only_gens(&self,gens:&Bound<'_,PyTuple>)->PyResult<bool>{structural::has_only(&self.p,gens)}
    fn homogenize(&self,s:&Bound<'_,PyAny>)->PyResult<Self>{
        let e=to_expr(s)?;
        if e.as_symbol().is_none(){return Err(PyTypeError::new_err(format!("``Symbol`` expected, got {}",s.get_type().str()?)));}
        Ok(Self{p:symox::poly::structure::homogenize(&self.p,&e).map_err(err)?})
    }
    #[pyo3(signature=(x,m,n=None))]
    fn slice(&self,x:&Bound<'_,PyAny>,m:&Bound<'_,PyAny>,n:Option<&Bound<'_,PyAny>>)->PyResult<Self>{
        let(j,m,n)=match n{None=>(0,x.extract::<i64>()?,m.extract::<i64>()?),Some(n)=>(self.gen_index(Some(x))?,m.extract::<i64>()?,n.extract::<i64>()?)};
        let clamp=|k:i64|k.clamp(0,u32::MAX as i64) as u32;
        Ok(Self{p:symox::poly::structure::slice(&self.p,j,clamp(m),clamp(n)).map_err(err)?})
    }
    fn ltrim(&self,py:Python<'_>,gen:&Bound<'_,PyAny>)->PyResult<Self>{
        let j=self.gen_index(Some(gen))?;
        match symox::poly::structure::ltrim(&self.p,j).map_err(err)?{Some(p)=>Ok(Self{p}),None=>Err(structural::polynomial_error(py,"PolynomialError",format!("Cannot left trim {}",self.__repr__())))}
    }
    fn all_monoms<'py>(&self,py:Python<'py>)->PyResult<Vec<Bound<'py,PyTuple>>>{structural::dense(&self.p,py,false)}
    fn all_terms<'py>(&self,py:Python<'py>)->PyResult<Vec<Bound<'py,PyTuple>>>{structural::dense(&self.p,py,true)}

    fn monoms<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyTuple>>> {
        self.p.monoms().into_iter().map(|m| PyTuple::new(py, m)).collect()
    }

    fn terms<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyTuple>>> {
        let t = self.p.terms();
        let t = if t.is_empty() { vec![(vec![0; self.p.ngens], Expr::zero())] } else { t };
        t.into_iter()
            .map(|(m, c)| {
                let mt = PyTuple::new(py, m)?;
                PyTuple::new(py, [mt.into_any(), w(c).into_pyobject(py)?.into_any()])
            })
            .collect()
    }

    /// The reference `Poly.as_dict`: ascending monomials (`dmp_to_dict`),
    /// `zero=True` keeps the zero monomial of a zero polynomial and
    /// `native=True` gives the domain elements.
    #[pyo3(signature=(native=false,zero=false))]
    fn as_dict<'py>(slf: &Bound<'py, Self>, native: bool, zero: bool) -> PyResult<Bound<'py, PyDict>> {
        let py = slf.py();
        let this = slf.get();
        let d = PyDict::new(py);
        let domain = if native { Some(slf.getattr("domain")?) } else { None };
        let mut terms = this.p.terms();
        terms.reverse();
        if terms.is_empty() && zero {
            terms.push((vec![0; this.p.true_gens().len()], Expr::zero()));
        }
        for (m, c) in terms {
            let value = w(c).into_pyobject(py)?.into_any();
            let value = match &domain { Some(k) => k.call_method1("from_sympy", (value,))?, None => value };
            d.set_item(PyTuple::new(py, m)?, value)?;
        }
        Ok(d)
    }

    fn LC(&self) -> PyExpr {
        w(self.p.lc())
    }

    fn TC(&self) -> PyExpr {
        w(self.p.tc())
    }

    fn LM<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let m = self.p.monoms().into_iter().next().unwrap_or_default();
        PyTuple::new(py, m)
    }

    fn coeff_monomial(&self,monom:&Bound<'_,PyAny>)->PyResult<PyExpr> {
        let powers=if monom.downcast::<PyTuple>().is_ok() || monom.downcast::<PyList>().is_ok() {
            let powers=monom.extract::<Vec<u32>>().map_err(|_|PyValueError::new_err("monomial exponents must be nonnegative integers"))?;
            if powers.len()!=self.p.true_gens().len(){return Err(PyValueError::new_err("exponent of each generator must be specified"));}
            powers
        }else {
            let e=to_expr(monom)?;
            let p=Poly::from_expr(&e,self.p.true_gens()).map_err(err)?;
            if p.gens!=self.p.true_gens() || p.rep.len()!=1 || !p.rep.lc().is_one() || p.domain==Domain::RR {return Err(PyValueError::new_err(format!("expected a monomial, got {e}")));}
            let powers=p.rep.iter().next().unwrap().0.to_vec();
            powers
        };
        Ok(w(self.p.coeff_monomial(&powers)))
    }

    #[pyo3(signature = (*n))]
    fn nth(&self, n: Vec<u32>) -> PyExpr {
        w(self.p.coeff_monomial(&n))
    }

    fn length(&self) -> usize {
        self.p.terms().len()
    }

    #[getter]
    fn is_zero(&self) -> bool {
        self.p.is_zero()
    }

    #[getter]
    fn is_one(&self) -> bool {
        self.p.rep.is_one()
    }

    #[getter]
    fn is_ground(&self) -> bool {
        self.p.rep.iter().all(|(e, _)| e[..self.p.ngens].iter().all(|&k| k == 0))
    }

    #[getter]
    fn is_monomial(&self) -> bool {
        self.p.terms().len() <= 1
    }

    #[getter]
    fn is_univariate(&self) -> bool {
        self.p.ngens == 1
    }

    #[getter]
    fn is_multivariate(&self) -> bool {
        self.p.ngens > 1
    }

    #[getter]
    fn is_monic(&self) -> bool {
        self.p.lc().is_one()
    }

    #[getter]
    fn is_irreducible(&self) -> PyResult<bool> {
        if self.ex().is_some(){return Err(err(PolyError::Reference{kind:"DomainError",message:"factorization not supported over EX".into()}));}
        if let Some(f)=&self.p.algebraic {if f.rep.is_constant(){return Ok(true);}let(_,fs)=f.factor_list().map_err(err)?;return Ok(fs.len()==1 && fs[0].1==1);}
        if let Some(f)=&self.p.field {let(_,fs)=f.factor_list();return Ok(fs.len()==1 && fs[0].1==1);}
        let (_, fs) = self.p.factor_list();
        Ok(fs.len() == 1 && fs[0].1 == 1)
    }

    #[getter]
    fn is_sqf(&self) -> PyResult<bool> {
        if let Some(f)=self.ex(){return f.sqf_p().map_err(err);}
        if let Some(f)=&self.p.algebraic {return Ok(f.sqf_list().1.iter().all(|(_,k)|*k==1));}
        if let Some(f)=&self.p.field {return Ok(f.sqf_list().1.iter().all(|(_,k)|*k==1));}
        let (_, fs) = self.p.sqf_list();
        Ok(fs.iter().all(|(_, k)| *k == 1))
    }

    fn __add__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.checked_add(b), |a, b| symox::ops::add2(a, b))
    }
    fn __radd__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.checked_add(a), |a, b| symox::ops::add2(b, a))
    }
    fn __sub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.checked_sub(b), |a, b| symox::ops::sub(a, b))
    }
    fn __rsub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.checked_sub(a), |a, b| symox::ops::sub(b, a))
    }
    fn __mul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.checked_mul(b), |a, b| symox::ops::mul2(a, b))
    }
    fn __rmul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.checked_mul(a), |a, b| symox::ops::mul2(b, a))
    }
    fn __neg__(&self) -> PyPoly {
        PyPoly { p: (self.p.neg()).into() }
    }
    fn __pos__(&self) -> PyPoly {
        self.clone()
    }
    fn __pow__(&self, k: u32, _m: Option<&Bound<'_, PyAny>>) -> PyResult<PyPoly> {
        Ok(PyPoly { p: self.p.checked_pow(k).map_err(err)? })
    }
    fn __truediv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(w(symox::ops::div(&self.p.as_expr(), &to_expr(o)?)))
    }
    fn __floordiv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: (self.p.quo(&q).map_err(err)?).into() })
    }
    fn __mod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: (self.p.rem(&q).map_err(err)?).into() })
    }
    fn __divmod__<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        self.div(py, o, true)
    }

    fn __eq__(&self, o: &Bound<'_, PyAny>) -> bool {
        if let Ok(q) = o.downcast::<PyPoly>() {
            let q = &q.get().p;
            // SymPy's Poly.__eq__ also requires equal coefficient domains.
            return self.p.true_gens() == q.true_gens() && self.p.domain_str() == q.domain_str() && self.p.modulus == q.modulus && self.p.as_expr() == q.as_expr();
        }
        match self.other(o) {
            Some(q) => self.p.as_expr() == q.as_expr(),
            None => false,
        }
    }

    fn __ne__(&self, o: &Bound<'_, PyAny>) -> bool {
        !self.__eq__(o)
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.p.as_expr().hash(&mut h);
        h.finish()
    }

    fn __bool__(&self) -> bool {
        !self.p.is_zero()
    }

    #[pyo3(signature=(o,auto=true))]
    fn div<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>, auto:bool) -> PyResult<Bound<'py, PyTuple>> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        let (a, b) = self.p.div_with_options(&q,auto).map_err(err)?;
        PyTuple::new(py, [PyPoly { p: (a).into() }, PyPoly { p: (b).into() }])
    }

    #[pyo3(signature=(o,auto=true))]
    fn rem(&self,o:&Bound<'_,PyAny>,auto:bool)->PyResult<PyPoly>{let q=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial"))?;Ok(PyPoly{p:self.p.rem_with_options(&q,auto).map_err(err)?})}

    #[pyo3(signature=(o,auto=true))]
    fn quo(&self,o:&Bound<'_,PyAny>,auto:bool)->PyResult<PyPoly>{let q=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial"))?;Ok(PyPoly{p:self.p.quo_with_options(&q,auto).map_err(err)?})}

    #[pyo3(signature=(o,auto=true))]
    fn exquo(&self, o: &Bound<'_, PyAny>, auto:bool) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: (self.p.exquo_with_options(&q,auto).map_err(err)?).into() })
    }

    fn pdiv<'py>(&self,py:Python<'py>,o:&Bound<'_,PyAny>)->PyResult<Bound<'py,PyTuple>>{
        let other=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial with compatible coefficient domain"))?;
        let(q,r)=self.p.pdiv(&other).map_err(err)?;PyTuple::new(py,[PyPoly{p:q},PyPoly{p:r}])
    }
    fn prem(&self,o:&Bound<'_,PyAny>)->PyResult<PyPoly>{let other=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial with compatible coefficient domain"))?;Ok(PyPoly{p:self.p.prem(&other).map_err(err)?})}
    fn pquo(&self,o:&Bound<'_,PyAny>)->PyResult<PyPoly>{let other=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial with compatible coefficient domain"))?;Ok(PyPoly{p:self.p.pquo(&other).map_err(err)?})}
    fn pexquo(&self,o:&Bound<'_,PyAny>)->PyResult<PyPoly>{let other=self.other(o).ok_or_else(||PyTypeError::new_err("expected a polynomial with compatible coefficient domain"))?;Ok(PyPoly{p:self.p.pexquo(&other).map_err(err)?})}

    #[pyo3(signature=(g,auto=true))]
    fn gcdex(&self,g:&Bound<'_,PyAny>,auto:bool)->PyResult<PyObject> {
        if let Some((a,b))=self.ex_pair(g)?{let(s,t,h)=a.gcdex(&b).map_err(err)?;return Ok(PyTuple::new(g.py(),[s,t,h].map(|p|PyPoly{p:p.into()}))?.into_any().unbind());}
        let q=self.other(g).ok_or_else(||PyTypeError::new_err("incompatible polynomial coefficient fields"))?;
        let values=an_rational::poly_bezout(&self.p,&q,false,auto)?;
        Ok(PyTuple::new(g.py(),values.into_iter().map(|p|PyPoly {p}))?.into_any().unbind())
    }
    #[pyo3(signature=(g,auto=true))]
    fn half_gcdex(&self,g:&Bound<'_,PyAny>,auto:bool)->PyResult<PyObject> {
        if let Some((a,b))=self.ex_pair(g)?{let(s,h)=a.half_gcdex(&b).map_err(err)?;return Ok(PyTuple::new(g.py(),[s,h].map(|p|PyPoly{p:p.into()}))?.into_any().unbind());}
        let q=self.other(g).ok_or_else(||PyTypeError::new_err("incompatible polynomial coefficient fields"))?;
        let values=an_rational::poly_bezout(&self.p,&q,true,auto)?;
        Ok(PyTuple::new(g.py(),values.into_iter().map(|p|PyPoly {p}))?.into_any().unbind())
    }
    #[pyo3(signature=(g,auto=true))]
    fn invert(&self,g:&Bound<'_,PyAny>,auto:bool)->PyResult<PyPoly> {
        if let Some((a,b))=self.ex_pair(g)?{return Ok(PyPoly{p:a.invert(&b).map_err(err)?.into()});}
        let q=self.other(g).ok_or_else(||PyTypeError::new_err("incompatible polynomial coefficient fields"))?;
        Ok(PyPoly {p:an_rational::poly_inverse(&self.p,&q,auto)?})
    }
    #[pyo3(signature=(g,include=false))]
    fn cancel(&self,g:&Bound<'_,PyAny>,include:bool)->PyResult<PyObject> {
        if let Some((a,b))=self.ex_pair(g)?{
            let py=g.py();let(cp,cq,p,q)=a.cancel(&b).map_err(err)?;
            if include{
                let p=p.mul_ground(&cp).map_err(err)?;let q=q.mul_ground(&cq).map_err(err)?;
                return Ok(PyTuple::new(py,[PyPoly{p:p.into()},PyPoly{p:q.into()}])?.into_any().unbind());
            }
            let c=w(symox::ops::div(&cp,&cq)).into_pyobject(py)?.into_any().unbind();
            return Ok(PyTuple::new(py,[c,PyPoly{p:p.into()}.into_pyobject(py)?.into_any().unbind(),PyPoly{p:q.into()}.into_pyobject(py)?.into_any().unbind()])?.into_any().unbind());
        }
        let q=self.other(g).ok_or_else(||PyTypeError::new_err("incompatible polynomial coefficient fields"))?;
        let py=g.py();
        if let Some((cp,cq,p,q))=an_rational::gaussian_poly_cancel(&self.p,&q)?{
            let mut values=vec![];
            if include{
                values.push(PyPoly{p:p.with_rep(p.rep.scale(&cp)).into()}.into_pyobject(py)?.into_any().unbind());
                values.push(PyPoly{p:q.with_rep(q.rep.scale(&cq)).into()}.into_pyobject(py)?.into_any().unbind());
            }else{
                values.push(w(symox::ops::div(&cp.as_expr(),&cq.as_expr())).into_pyobject(py)?.into_any().unbind());
                values.push(PyPoly{p:p.into()}.into_pyobject(py)?.into_any().unbind());
                values.push(PyPoly{p:q.into()}.into_pyobject(py)?.into_any().unbind());
            }
            return Ok(PyTuple::new(py,values)?.into_any().unbind());
        }
        if let Some((ratio,p,q))=an_rational::fraction_poly_cancel(&self.p,&q,include)?{
            let mut values=vec![];
            if let Some(c)=ratio{values.push(w(c).into_pyobject(py)?.into_any().unbind());}
            values.push(PyPoly{p}.into_pyobject(py)?.into_any().unbind());
            values.push(PyPoly{p:q}.into_pyobject(py)?.into_any().unbind());
            return Ok(PyTuple::new(py,values)?.into_any().unbind());
        }
        let(a,b)=an_rational::poly_cancel(&self.p,&q)?;
        if !include{if let Some((c,a,b))=sympy_cancel_parts(&self.p,&q,&a,&b){
            let c=w(Expr::num(c.to_number())).into_pyobject(py)?.into_any().unbind();
            return Ok(PyTuple::new(py,[c,PyPoly{p:a}.into_pyobject(py)?.into_any().unbind(),PyPoly{p:b}.into_pyobject(py)?.into_any().unbind()])?.into_any().unbind());
        }}
        let mut values=vec![];
        if !include{values.push(w(Expr::one()).into_pyobject(py)?.into_any().unbind());}
        values.push(PyPoly {p:a}.into_pyobject(py)?.into_any().unbind());
        values.push(PyPoly {p:b}.into_pyobject(py)?.into_any().unbind());
        Ok(PyTuple::new(py,values)?.into_any().unbind())
    }

    fn gcd(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.checked_gcd(&q).map_err(err)? })
    }

    fn lcm(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.checked_lcm(&q).map_err(err)? })
    }

    fn cofactors<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        let (h, a, b) = self.p.checked_cofactors(&q).map_err(err)?;
        PyTuple::new(py, [PyPoly { p: (h).into() }, PyPoly { p: (a).into() }, PyPoly { p: (b).into() }])
    }

    fn factor_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if let Some(f)=self.ex(){
            // dmp_factor_list raises DomainError over EX; Poly.factor_list
            // then returns the constant or the polynomial itself.
            if f.degree(0)==0 {
                return PyTuple::new(py,[w(f.as_expr()).into_pyobject(py)?.into_any(),PyList::empty(py).into_any()]);
            }
            let item=PyTuple::new(py,[self.clone().into_pyobject(py)?.into_any(),1u32.into_pyobject(py)?.into_any()])?;
            return PyTuple::new(py,[w(Expr::one()).into_pyobject(py)?.into_any(),PyList::new(py,[item])?.into_any()]);
        }
        checked_factor_views::list(&self.p,py,false,false,true)
    }

    fn ground_roots(&self,py:Python<'_>)->PyResult<PyObject>{
        if let Some(f)=self.ex(){
            if f.gens.len()!=1{return Err(err(PolyError::Reference{kind:"MultivariatePolynomialError",message:f.as_expr().to_string()}));}
            let result=PyDict::new(py);
            // The only factor of a nonconstant EX polynomial is itself.
            if f.degree(0)==1{
                let c=f.all_coeffs();
                result.set_item(w(symox::ops::div(&symox::ops::neg(&c[1]),&c[0])),1u32)?;
            }
            return Ok(result.into_any().unbind());
        }
        checked_factor_views::ground(&self.p,py)
    }
    fn gff_list(&self)->PyResult<Vec<(PyPoly,u32)>>{
        if let Some(f)=self.ex(){if f.is_zero(){return Err(PyValueError::new_err("greatest factorial factorization doesn't exist for a zero polynomial"));}return Ok(f.gff_list().map_err(err)?.into_iter().map(|(p,k)|(PyPoly{p:p.into()},k as u32)).collect());}
        checked_factor_views::factorial(&self.p)
    }
    fn factor_list_include(&self)->PyResult<Vec<(PyPoly,u32)>>{
        if self.ex().is_some(){return Ok(vec![(self.clone(),1)]);}
        checked_factor_views::include(&self.p,false,false)
    }
    #[pyo3(signature=(all=false))]
    fn sqf_list_include(&self,all:bool)->PyResult<Vec<(PyPoly,u32)>>{
        if let Some(f)=self.ex(){return Ok(f.sqf_list_include(all).map_err(err)?.into_iter().map(|(p,k)|(PyPoly{p:p.into()},k as u32)).collect());}
        checked_factor_views::include(&self.p,true,all)
    }

    #[pyo3(signature=(all=false))]
    fn sqf_list<'py>(&self, py: Python<'py>, all:bool) -> PyResult<Bound<'py, PyTuple>> {
        if let Some(f)=self.ex(){
            let(c,fs)=f.sqf_list(all).map_err(err)?;
            let items=fs.into_iter().map(|(p,k)|PyTuple::new(py,[PyPoly{p:p.into()}.into_pyobject(py)?.into_any(),(k as u32).into_pyobject(py)?.into_any()])).collect::<PyResult<Vec<_>>>()?;
            return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),PyList::new(py,items)?.into_any()]);
        }
        checked_factor_views::list(&self.p,py,true,all,true)
    }

    fn primitive<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if let Some(f)=self.ex(){let(c,q)=f.primitive().map_err(err)?;return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),PyPoly{p:q.into()}.into_pyobject(py)?.into_any()]);}
        if let Some(a)=&self.p.algebraic{
            if a.gaussian.is_some(){let(c,q)=a.primitive();return PyTuple::new(py,[w(c.as_expr()).into_pyobject(py)?.into_any(),PyPoly{p:q.into()}.into_pyobject(py)?.into_any()]);}
            return PyTuple::new(py,[w(if self.p.is_zero(){Expr::zero()}else{Expr::one()}).into_pyobject(py)?.into_any(),self.clone().into_pyobject(py)?.into_any()]);
        }
        if let Some(f)=&self.p.field {
            let(c,p)=if f.coefficient_domain.is_field(){f.reference_primitive()}else{f.primitive()};let c=symox::poly::fraction::coefficient_expr(&c,&f.parameters);
            return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),poly_obj(py,p)?.into_bound(py)]);
        }
        let (c, p) = self.p.primitive();
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), PyPoly { p: p.into() }.into_pyobject(py)?.into_any()])
    }

    fn content(&self) -> PyExpr {
        if let Some(f)=self.ex(){return w(f.content());}
        if let Some(a)=&self.p.algebraic{if a.gaussian.is_some(){return w(a.primitive().0.as_expr());}return w(if self.p.is_zero(){Expr::zero()}else{Expr::one()});}
        if let Some(f)=&self.p.field {return w(symox::poly::fraction::coefficient_expr(&if f.coefficient_domain.is_field(){f.reference_content()}else{f.primitive().0},&f.parameters));}
        w(rat_expr(&self.p.primitive().0))
    }

    #[pyo3(signature=(auto=true))]
    fn monic(&self, auto:bool) -> PyResult<PyPoly> {
        if let Some(f)=self.ex(){return Ok(PyPoly{p:f.monic().map_err(err)?.into()});}
        let domain = self.domain_descriptor()?;
        let domain = if auto {domain.get_field()} else {domain};
        // Rebuild implicit coefficient rings with their true generators before
        // selecting the leading coefficient. This also preserves typed orders.
        let p = domain.poly(&self.p.as_expr(), self.p.true_gens())?;
        if p.is_zero() {return Ok(PyPoly{p});}
        let divisor = domain.poly(&p.lc(), p.true_gens())?;
        Ok(PyPoly{p:p.exquo_with_options(&divisor,false).map_err(err)?})
    }

    #[pyo3(signature = (*gens))]
    fn diff(&self, gens: &Bound<'_, PyTuple>) -> PyResult<PyPoly> {
        let mut p = self.p.clone();
        if gens.is_empty() {
            return Ok(PyPoly { p: (p.diff(0)).into() });
        }
        for spec in gens.iter() {
            let (generator, order) = if let Ok(spec) = spec.downcast::<PyTuple>() {
                if spec.len()!=2 {
                    return Err(PyValueError::new_err("derivative specification must be (generator, order)"));
                }
                (spec.get_item(0)?, spec.get_item(1)?.extract::<i64>()?)
            } else {(spec,1)};
            let i = self.gen_index(Some(&generator))?;
            for _ in 0..order.max(0) {
                if p.is_zero() {break;}
                p = p.diff(i);
            }
        }
        Ok(PyPoly { p: p.into() })
    }

    #[pyo3(signature = (x, a=None))]
    fn eval(&self, py: Python<'_>, x: &Bound<'_, PyAny>, a: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
        let (i, val) = match a {
            Some(a) => (self.gen_index(Some(x))?, to_expr(a)?),
            None => (0, to_expr(x)?),
        };
        let g = self.p.true_gens()[i].clone();
        let r = symox::subs::subs(&self.p.as_expr(), &[(g, val)]);
        let r = symox::expand::expand(&r);
        if let Some(field)=&self.p.algebraic {
            if self.p.ngens>1 {let mut gens=field.gens.clone();gens.remove(i);return poly_obj(py,AlgebraicPoly::from_expr(&r,&gens,field.number_field.clone()).map_err(err)?);}
            let value=field.convert(&r).map_err(err)?;
            if !value.rep.is_constant(){return Err(PyValueError::new_err("evaluation value must belong to the coefficient field"));}
            return Ok(w(value.as_expr()).into_pyobject(py)?.into_any().unbind());
        }
        if let Some(field)=&self.p.field {
            if self.p.ngens>1 {let mut gens=field.gens.clone();gens.remove(i);return poly_obj(py,FractionPoly::from_expr(&r,&gens,&field.parameters,&field.display).map_err(err)?);}
            let value=field.convert(&r).map_err(err)?;
            if !value.rep.is_constant(){return Err(PyValueError::new_err("evaluation value must belong to the coefficient field"));}
            return Ok(w(value.as_expr()).into_pyobject(py)?.into_any().unbind());
        }
        if self.p.ngens > 1 {
            let mut gens = self.p.true_gens().to_vec();
            gens.remove(i);
                if let Ok(p) = Poly::from_expr(&r, &gens) {
                return poly_obj(py, p);
            }
        }
        Ok(w(r).into_pyobject(py)?.into_any().unbind())
    }

    #[pyo3(signature = (*args))]
    fn __call__(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
        if args.len() > self.p.true_gens().len() {
            return Err(PyValueError::new_err("too many values provided"));
        }
        let mut pairs = Vec::new();
        for (i, a) in args.iter().enumerate() {
            pairs.push((self.p.true_gens()[i].clone(), to_expr(&a)?));
        }
        Ok(w(symox::expand::expand(&symox::subs::subs(&self.p.as_expr(), &pairs))))
    }

    #[pyo3(signature = (n=15, maxsteps=50, cleanup=true))]
    fn nroots(&self, py: Python<'_>, n: u32, maxsteps: u32, cleanup: bool) -> PyResult<Vec<PyExpr>> {
        let generators = self.ex().map_or(self.p.true_gens().len(), |f| f.gens.len());
        if generators > 1 {
            let message = format!("Cannot compute numerical roots of {}", self.__repr__());
            return Err(err(PolyError::Reference { kind: "MultivariatePolynomialError", message }));
        }
        let expression = match self.ex() {
            Some(f) => f.clone(),
            None => {
                let plain = self.p.field.is_none() && self.p.algebraic.is_none();
                if plain {
                    match proots::nroots_poly_options(&self.p, n, maxsteps, cleanup) {
                        Ok(v) => return Ok(v.into_iter().map(w).collect()),
                        // Symbolic coefficient domains are evaluated like EX.
                        Err(PolyError::Unsupported(_)) if self.p.modulus.is_none() && generators == 1 => {}
                        Err(error) => return Err(nroots_error(py, error)),
                    }
                }
                ExpressionPoly::from_expr(&self.p.as_expr(), self.p.true_gens()).map_err(err)?
            }
        };
        let domain = if self.ex().is_some() { "EX".to_string() } else { self.domain_descriptor().map_or_else(|_| "EX".to_string(), |d| d.inner.to_string()) };
        Ok(proots::expression_nroots(&expression, &domain, n, maxsteps, cleanup)
            .map_err(|error| nroots_error(py, error))?.into_iter().map(w).collect())
    }

    fn all_roots(&self) -> PyResult<Vec<PyExpr>> {
        let rs = proots::roots_poly(&self.p).map_err(err)?;
        let mut v = Vec::new();
        for (r, k) in rs {
            for _ in 0..k {
                v.push(w(r.clone()));
            }
        }
        Ok(v)
    }

    fn compose(&self, g: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        if let Some((a,b))=self.ex_pair(g)?{return Ok(PyPoly{p:a.compose(&b).map_err(err)?.into()});}
        if let Some(f)=&self.p.algebraic {if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("algebraic compose requires one generator"));}let q=f.convert(&to_expr(g)?).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        if let Some(f)=&self.p.field {let q=f.convert(&to_expr(g)?).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        let q = self.other(g).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: (self.p.compose(&q)).into() })
    }

    fn decompose(&self) -> PyResult<Vec<PyPoly>> {
        Ok(self.p.decompose().map_err(|error|match error {
            PolyError::Unsupported(message)=>PyNotImplementedError::new_err(message),
            other=>err(other),
        })?.into_iter().map(|p| PyPoly { p }).collect())
    }

    fn shift(&self, a: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        if let Some(f)=self.ex(){return Ok(PyPoly{p:f.shift(&to_expr(a)?).map_err(err)?.into()});}
        if let Some(f)=&self.p.algebraic {if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("algebraic shift requires one generator"));}let argument=symox::ops::add2(&f.gens[0],&to_expr(a)?);let q=f.convert(&argument).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        if let Some(f)=&self.p.field {let argument=symox::ops::add2(&f.gens[0],&to_expr(a)?);let q=f.convert(&argument).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        let e = to_expr(a)?;
        let r = e
            .as_number()
            .and_then(Rat::from_number)
            .ok_or_else(|| PyValueError::new_err("shift expects a rational number"))?;
        Ok(PyPoly { p: (self.p.shift(&r)).into() })
    }

    #[pyo3(signature = (*gens, auto=true))]
    fn integrate(&self, gens: &Bound<'_, PyTuple>, auto:bool) -> PyResult<PyPoly> {
        let mut p=if self.ex().is_some(){self.p.clone()}else{self.p.integration_domain(auto)};
        if gens.is_empty(){return Ok(PyPoly{p:p.integrate_checked(0,false).map_err(err)?});}
        for spec in gens.iter(){
            let (generator,count)=if let Ok(spec)=spec.downcast::<PyTuple>(){
                if spec.len()!=2{return Err(PyValueError::new_err("integration specification must be (generator, order)"));}
                (spec.get_item(0)?,spec.get_item(1)?.extract::<i64>()?)
            }else{(spec,1)};
            let i=if let Ok(index)=generator.extract::<i64>(){let index=if index<0{self.p.ngens as i64+index}else{index};if index<0||index>=self.p.ngens as i64{return Err(PyValueError::new_err("integration generator is out of range"));}index as usize}else{self.gen_index(Some(&generator))?};
            for _ in 0..count.max(0){p=p.integrate_checked(i,false).map_err(err)?;}
        }
        Ok(PyPoly {p})
    }

    #[pyo3(signature=(convert=false))]
    fn clear_denoms<'py>(&self, py: Python<'py>, convert:bool) -> PyResult<Bound<'py, PyTuple>> {
        if let Some(f)=self.ex(){let _=convert;let(c,q)=f.clear_denoms().map_err(err)?;return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),PyPoly{p:q.into()}.into_pyobject(py)?.into_any()]);}
        let domain = self.domain_descriptor()?;
        // Rebuilding from the expression expands composite generators
        // ((C1 + C2*x)*exp(-x) with the generator C1 + C2*x); a ground
        // domain representation is already canonical.
        // Rebuilding would also merge dependent generators ((x, x - 1)
        // from terms_gcd(..., expand=False)), so ZZ and QQ keep theirs.
        let ground = self.p.field.is_none() && self.p.algebraic.is_none()
            && self.p.base.modulus.is_none() && matches!(self.p.base.domain, Domain::ZZ | Domain::QQ);
        let canonical = if ground { self.p.clone() } else { match domain.poly(&self.p.as_expr(),self.p.true_gens()) {
            Ok(p) => p,
            Err(_) if self.p.field.is_none() && self.p.algebraic.is_none() => self.p.clone(),
            Err(e) => return Err(e),
        } };
        let (denominator, numerator) = if let Some(f) = &canonical.field {
            let ring = domain.inner.ring().map_err(err)?;
            let denominator = if matches!(&domain.inner,symox::poly::domains::ExactDomain::Functions{ring:true,..}) {
                ring.convert(&Expr::one()).map_err(err)?
            } else {
                let expression = symox::poly::convert::to_expr_int(&f.denominator(),&f.parameters);
                let value = ring.convert(&expression).map_err(err)?;
                // LCM with one fixes the unit using the associated ring's
                // ground coefficients and monomial order (QQ[a] vs ZZ[a]).
                ring.convert(&Expr::one()).map_err(err)?.lcm(&value).map_err(err)?
            };
            let expression = denominator.as_expr();
            let numerator:PolyValue = f.mul(&f.convert(&expression).map_err(err)?).into();
            let numerator = if convert {
                domains::PyDomain{inner:ring}.poly(&numerator.as_expr(),numerator.true_gens())?
            } else {numerator};
            (expression,numerator)
        } else if let Some(a) = canonical.algebraic.as_ref().filter(|a|a.gaussian.is_some()) {
            // dmp_clear_denoms over QQ_I: the integer lcm of the coefficient
            // denominators, into ZZ_I with convert=True.
            let (d,rep)=symox::poly::gaussian::clear_denoms(&a.rep);
            let mut q=a.with_rep(rep);
            if convert{q.gaussian=Some(true);}
            (Expr::num(Rat::from(d).to_number()),q.into())
        } else if canonical.algebraic.is_some() {
            (Expr::one(),canonical)
        } else {
            let (d,mut p) = canonical.base.clear_denoms();
            if convert && p.domain==Domain::QQ && p.modulus.is_none() {p.domain=Domain::ZZ;}
            (Expr::num(d.to_number()),p.into())
        };
        PyTuple::new(py,[w(denominator).into_pyobject(py)?.into_any(),
                         PyPoly{p:numerator}.into_pyobject(py)?.into_any()])
    }

    fn max_norm(&self) -> PyResult<PyExpr> {
        if let Some(f)=self.ex(){return Ok(w(f.max_norm()));}
        if self.p.field.is_some()||self.p.algebraic.is_some(){return Err(PyTypeError::new_err("rational function coefficients have no ordered absolute norm"));}
        Ok(w(rat_expr(&self.p.max_norm())))
    }

    fn l1_norm(&self) -> PyResult<PyExpr> {
        if let Some(f)=self.ex(){return Ok(w(f.l1_norm().map_err(err)?));}
        if self.p.field.is_some()||self.p.algebraic.is_some(){return Err(PyTypeError::new_err("rational function coefficients have no ordered absolute norm"));}
        Ok(w(rat_expr(&self.p.l1_norm())))
    }

    fn sqf_part(&self) -> PyResult<PyPoly> {
        if let Some(f)=self.ex(){return Ok(PyPoly{p:f.sqf_part().map_err(err)?.into()});}
        checked_factor_views::part(&self.p)
    }

    #[getter]
    fn is_linear(&self) -> bool {
        self.p.is_linear()
    }

    #[getter]
    fn is_quadratic(&self) -> bool {
        self.p.is_quadratic()
    }

    #[getter]
    fn is_homogeneous(&self) -> bool {
        self.p.is_homogeneous()
    }

    fn as_list(&self,py:Python<'_>) -> PyResult<Vec<PyExpr>> {
        self.all_coeffs(py)
    }

    #[pyo3(signature = (*args, **kwargs))]
    fn subs(&self, args: &Bound<'_, PyTuple>,kwargs:Option<&Bound<'_,PyDict>>) -> PyResult<PyObject> {
        let (args,options)=crate::substitution_arguments(args,kwargs)?;
        let py = args.py();
        let e = self.p.as_expr();
        let pairs=crate::subs_pairs(&args)?;
        let r = symox::expand::expand(&crate::checked_subs(py,&e,&pairs,crate::substitution_options(options.as_ref())?)?);
        if let Some(field)=&self.p.algebraic {return poly_obj(py,field.convert(&r).map_err(err)?);}
        if let Some(field)=&self.p.field {return poly_obj(py,field.convert(&r).map_err(err)?);}
        if let Ok(p) = Poly::from_expr(&r, self.p.true_gens()) {
            if !r.is_number() {
                return poly_obj(py, p);
            }
        }
        Ok(w(r).into_pyobject(py)?.into_any().unbind())
    }

    fn termwise(&self, func: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let mut terms = Vec::new();
        for (m, c) in self.p.terms() {
            let mt = PyTuple::new(func.py(), m.clone())?;
            let r = func.call1((mt, w(c)))?;
            let e = to_expr(&r)?;
            let mut f = vec![e];
            for (g, k) in self.p.true_gens().iter().zip(&m) {
                f.push(symox::ops::pow(g, &Expr::int(*k as i64)));
            }
            terms.push(symox::ops::mul(&f));
        }
        let e = symox::ops::add(&terms);
        if let Some(field)=&self.p.algebraic {return Ok(PyPoly {p:field.convert(&e).map_err(err)?.into()});}
        if let Some(field)=&self.p.field {return Ok(PyPoly {p:field.convert(&e).map_err(err)?.into()});}
        let mut p = Poly::from_expr(&e, self.p.true_gens()).map_err(err)?;
        if p.domain < self.p.domain {
            p.domain = self.p.domain;
        }
        if let Some(m) = self.p.modulus {
            p = p.set_modulus(m);
        }
        Ok(PyPoly { p: p.into() })
    }

    fn as_poly(&self) -> PyPoly {
        self.clone()
    }

    fn norm(&self)->PyResult<PyPoly> {
        if self.ex().is_some(){return Err(err(PolyError::Reference{kind:"DomainError",message:"ground domain must be algebraic".into()}));}
        let Some(p)=&self.p.algebraic else{return Err(PyTypeError::new_err("polynomial norm requires an algebraic coefficient field"))};
        Ok(PyPoly {p:p.norm().map_err(err)?.into()})
    }

    #[pyo3(signature=(g,includePRS=false))]
    fn resultant(&self, g: &Bound<'_, PyAny>,includePRS:bool) -> PyResult<PyObject> {
        let q = self.other(g).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        let(result,sequence)=self.p.resultant_prs(&q).map_err(err)?;
        let py=g.py();let result=division::render_eliminated(py,result,true)?;
        if !includePRS{return Ok(result);}
        let sequence=PyList::new(py,sequence.into_iter().map(|p|PyPoly{p}).collect::<Vec<_>>())?.into_any().unbind();
        Ok(PyTuple::new(py,[result,sequence])?.into_any().unbind())
    }

    fn discriminant(&self,py:Python<'_>) -> PyResult<PyObject> {
        division::render_eliminated(py,self.p.discriminant_poly().map_err(err)?,true)
    }

    fn subresultants(&self,g:&Bound<'_,PyAny>)->PyResult<Vec<PyPoly>> {
        let q=self.other(g).ok_or_else(||PyTypeError::new_err("expected a polynomial with a compatible coefficient domain"))?;
        Ok(self.p.subresultants(&q).map_err(err)?.into_iter().map(|p|PyPoly{p}).collect())
    }
    #[pyo3(signature=(auto=true))]
    fn sturm(&self,py:Python<'_>,auto:bool)->PyResult<Vec<PyPoly>>{root_analysis::sturm_method(py,&self.p,auto)}
    #[pyo3(signature=(inf=None,sup=None))]
    fn count_roots(&self,py:Python<'_>,inf:Option<&Bound<'_,PyAny>>,sup:Option<&Bound<'_,PyAny>>)->PyResult<PyExpr>{root_analysis::count_method(py,&self.p,inf,sup)}
    #[pyo3(signature=(all=false,eps=None,inf=None,sup=None,fast=false,sqf=false))]
    fn intervals(&self,py:Python<'_>,all:bool,eps:Option<&Bound<'_,PyAny>>,inf:Option<&Bound<'_,PyAny>>,sup:Option<&Bound<'_,PyAny>>,fast:bool,sqf:bool)->PyResult<PyObject>{root_analysis::intervals_method(py,&self.p,all,eps,inf,sup,fast,sqf)}
    #[pyo3(signature=(s,t,eps=None,steps=None,fast=false,check_sqf=false))]
    fn refine_root(&self,s:&Bound<'_,PyAny>,t:&Bound<'_,PyAny>,eps:Option<&Bound<'_,PyAny>>,steps:Option<i64>,fast:bool,check_sqf:bool)->PyResult<(PyExpr,PyExpr)>{root_analysis::refine_method(&self.p,s,t,eps,steps,fast,check_sqf)}
}

// ----------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------

fn expression_vector(values:&Bound<'_,PyAny>)->PyResult<Vec<Expr>> {
    values.try_iter()?.map(|v|to_expr(&v?)).collect()
}

#[pyfunction]
#[pyo3(signature=(expr,generators,parameters,domain,modulus))]
fn _rebuild_poly(expr:&Bound<'_,PyAny>,generators:&Bound<'_,PyAny>,parameters:Option<&Bound<'_,PyAny>>,domain:&Bound<'_,PyAny>,modulus:Option<u64>)->PyResult<PyPoly> {
    let e=to_expr(expr)?;let generators=expression_vector(generators)?;
    if let Ok(d)=domain.extract::<PyRef<'_,domains::PyDomain>>(){return Ok(PyPoly{p:d.poly(&e,&generators)?});}
    let domain_string=domain.extract::<String>()?;let domain=domain_string.as_str();
    if let Some(parameters)=parameters {
        let parameters=expression_vector(parameters)?;
        return Ok(PyPoly {p:FractionPoly::from_expr(&e,&generators,&parameters,domain).map_err(err)?.into()});
    }
    if let Some(field)=NumberField::parse(domain).map_err(err)?{return Ok(PyPoly {p:AlgebraicPoly::from_expr(&e,&generators,field).map_err(err)?.into()});}
    let mut p=Poly::from_expr(&e,&generators).map_err(err)?;
    let base=domain.split(['[','(']).next().unwrap_or(domain);
    p=set_domain(p,Some(base))?;
    if let Some(m)=modulus {p=p.set_modulus(check_mod(m)?);}
    Ok(PyPoly {p:domains::normalize_poly(p.into())?})
}

#[pyfunction]
#[pyo3(signature=(expressions,generators,parameters,domain,order,polys,modulus,options=None))]
fn _rebuild_groebner(py:Python<'_>,expressions:&Bound<'_,PyAny>,generators:&Bound<'_,PyAny>,parameters:Option<&Bound<'_,PyAny>>,domain:&Bound<'_,PyAny>,order:&str,polys:bool,modulus:Option<u64>,options:Option<String>)->PyResult<PyGroebnerBasis> {
    let generators=expression_vector(generators)?;
    let args=PyTuple::new(py,generators.iter().cloned().map(w))?;
    let mut result=PyGroebnerBasis::build(expressions,&args,order,if modulus.is_some(){None}else{Some(domain)},modulus,None,Some(polys),None)?;
    if let Some(options)=options{result.options=options;}
    if let Some(parameters)=parameters.filter(|_|domain.extract::<String>().is_ok()) {
        let parameters=expression_vector(parameters)?;
        let legacy=domain.extract::<String>()?;
        let prototype=FractionPoly::from_expr(&Expr::zero(),&generators,&parameters,&legacy).map_err(err)?;
        result.basis=result.basis.iter().map(|p|prototype.convert(&p.as_expr()).map(Into::into)).collect::<Result<_,_>>().map_err(err)?;
        result.fraction_domain=Some(prototype);
    }
    Ok(result)
}

fn algebraic_option_poly(f:&Bound<'_,PyAny>,gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Option<AlgebraicPoly>> {
    algebraic_option_polys(f,&[],gens,kwargs)
}

/// Like `algebraic_option_poly`, but the field and the generators are
/// inferred from `f` together with `others`, as the reference does for
/// `parallel_poly_from_expr` in binary operations such as `gcd`.
fn algebraic_option_polys(f:&Bound<'_,PyAny>,others:&[Expr],gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Option<AlgebraicPoly>> {
    let (extension,domain)=match kwargs{
        Some(kwargs)=>(kwargs.get_item("extension")?,match division::gaussian_option(kwargs)? {Some(d)=>Some(d),None=>kwargs.get_item("domain")?}),
        None=>(None,None),
    };
    let extension=extension.filter(|x|!x.is_none());
    let domain=domain.filter(|x|!x.is_none());
    if extension.is_none()&&domain.is_none(){
        return default_gaussian_polys(f,others,gens,kwargs);
    }
    let selected=domains::Selection::new(domain.as_ref(),None,extension.as_ref())?;
    if selected.native.as_ref().is_some_and(|d|!matches!(d.inner,symox::poly::domains::ExactDomain::ZZ|symox::poly::domains::ExactDomain::QQ|symox::poly::domains::ExactDomain::Algebraic(_)|symox::poly::domains::ExactDomain::Gaussian{..})){return Err(PyNotImplementedError::new_err("this expression operation does not support the selected coefficient domain; use its supported Poly method or modulus option"));}
    let(e,g2)=expr_or_coeff_list(f)?;
    let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}
    let mut all=vec![e.clone()];all.extend(others.iter().cloned());
    let Some(field)=selected.algebraic_in(extension.as_ref(),&all,&generators)?else{return Ok(None)};
    if generators.is_empty(){generators=symox::poly::algebraic::infer_generators(&all,&field);}
    if let Some(d)=selected.native.as_ref().filter(|_|selected.gaussian().is_some()) {
        return Ok(d.poly(&e,&generators)?.algebraic);
    }
    AlgebraicPoly::from_expr(&e,&generators,field).map(Some).map_err(err)
}


#[pyfunction]
#[pyo3(signature = (f, *gens, deep=false, modulus=None, **_kwargs))]
fn factor(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, deep: bool, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let value=factor_expression(f,gens,deep,modulus,_kwargs)?;
    rebuild_poly(f,value)
}

/// The reference `_symbolic_factor` maps over the arguments of a `Poly`,
/// so factoring one returns `Poly(factored, *gens)` with its domain
/// inferred again from the factored expression.
fn rebuild_poly(f: &Bound<'_, PyAny>, value: PyExpr) -> PyResult<PyObject> {
    let py=f.py();
    if f.downcast::<PyPoly>().is_err(){return Ok(value.into_pyobject(py)?.into_any().unbind());}
    let mut args=vec![value.into_pyobject(py)?.into_any()];
    for g in f.getattr("gens")?.try_iter()?{args.push(g?);}
    Ok(f.get_type().call1(PyTuple::new(py,args)?)?.unbind())
}

fn factor_expression(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, deep: bool, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    // Over EX every multiplicative component is converted on its own.
    if modulus.is_none()&&f.downcast::<PyPoly>().is_err()&&domains::selects_expression(_kwargs.map(|k|k.get_item("domain")).transpose()?.flatten().as_ref()){
        if !deep{return Ok(w(checked_factor_views::domain_expression(f,gens,_kwargs,false)?));}
        // The reference `deep=True`: factor every sum and product and
        // keep the results that became products or powers.
        let e=to_expr(f)?;
        let mut parts:Vec<Expr>=Vec::new();
        symox::traverse::preorder(&e,&mut |x|{if matches!(x.kind(),symox::expr::Kind::Add(..)|symox::expr::Kind::Mul(..))&&!parts.contains(x){parts.push(x.clone());}});
        let mut partials=Vec::new();
        for p in parts{
            let value=w(p.clone()).into_pyobject(f.py())?.into_any();
            let factored=checked_factor_views::domain_expression(&value,gens,_kwargs,false)?;
            if matches!(factored.kind(),symox::expr::Kind::Mul(..)|symox::expr::Kind::Pow(..))&&factored!=p{partials.push((p,factored));}
        }
        return Ok(w(symox::simptools::misc::xreplace(&e,&partials)));
    }
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if f.downcast::<PyPoly>().is_err()&&!deep{return Ok(w(checked_factor_views::domain_expression(f,gens,_kwargs,false)?));}
        let(c,fs)=p.factor_list().map_err(err)?;let mut factors=vec![c.as_expr()];for(q,k)in fs{factors.push(q.as_expr().powi(k as i64));}return Ok(w(symox::ops::mul(&factors)));
    }
    let _ = gens;
    let e = to_expr(f)?;
    if let Some(p) = modulus {
        return Ok(w(api::factor_mod(&e, check_mod(p)?).map_err(err)?));
    }
    if deep {
        return Ok(w(api::factor_deep(&e).map_err(err)?));
    }
    Ok(w(api::factor(&e).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, modulus=None, **_kwargs))]
fn sqf(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let value=sqf_expression(f,gens,modulus,_kwargs)?;
    rebuild_poly(f,value)
}

fn sqf_expression(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    // Over EX every multiplicative component is converted on its own.
    if modulus.is_none()&&f.downcast::<PyPoly>().is_err()&&domains::selects_expression(_kwargs.map(|k|k.get_item("domain")).transpose()?.flatten().as_ref()){
        return Ok(w(checked_factor_views::domain_expression(f,gens,_kwargs,true)?));
    }
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if f.downcast::<PyPoly>().is_err(){return Ok(w(checked_factor_views::domain_expression(f,gens,_kwargs,true)?));}
        let(c,fs)=p.sqf_list();let mut factors=vec![c.as_expr()];for(q,k)in fs{factors.push(q.as_expr().powi(k as i64));}return Ok(w(symox::ops::mul(&factors)));
    }
    let _ = gens;
    if let Some(p) = modulus {
        return Ok(w(api::sqf_mod(&to_expr(f)?, check_mod(p)?).map_err(err)?));
    }
    Ok(w(api::sqf(&to_expr(f)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, modulus=None, **_kwargs))]
fn factor_list<'py>(f: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'py, PyDict>>) -> PyResult<PyObject> {
    checked_factor_views::global_list(f,gens,modulus,_kwargs,false)
}

#[pyfunction]
#[pyo3(signature = (f, *gens, modulus=None, **_kwargs))]
fn sqf_list<'py>(f: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'py, PyDict>>) -> PyResult<PyObject> {
    checked_factor_views::global_list(f,gens,modulus,_kwargs,true)
}

/// `gcd_list`/`lcm_list` of expressions over an algebraic or Gaussian
/// domain selected by the options; a second argument is a generator.
fn algebraic_option_fold(py:Python<'_>,f:&Bound<'_,PyAny>,g:Option<&Bound<'_,PyAny>>,gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>,modulus:Option<u64>,op:fn(&AlgebraicPoly,&AlgebraicPoly)->AlgebraicPoly)->PyResult<Option<PyObject>>{
    if f.downcast::<PyPoly>().is_ok(){return Ok(None);}
    let Some(items)=list_items(f)? else{return Ok(None)};
    let Some((first,rest))=items.split_first() else{return Ok(None)};
    let mut all_gens:Vec<PyObject>=Vec::new();
    if let Some(g)=g{all_gens.push(g.clone().unbind());}
    for x in gens.iter(){all_gens.push(x.unbind());}
    let all_gens=PyTuple::new(py,all_gens)?;
    let head=w(first.clone()).into_pyobject(py)?.into_any();
    let Some(p)=algebraic_option_polys(&head,rest,&all_gens,kwargs)? else{return Ok(None)};
    if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
    let mut acc=p.clone();
    for x in rest{acc=op(&acc,&p.convert(x).map_err(err)?);}
    Ok(Some(w(acc.as_expr()).into_pyobject(py)?.into_any().unbind()))
}

/// The reference default domain of expressions whose only nonrational
/// coefficient atom is `I` is ZZ_I or QQ_I; other inputs keep their path.
fn default_gaussian_polys(f:&Bound<'_,PyAny>,others:&[Expr],gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Option<AlgebraicPoly>> {
    if f.downcast::<PyPoly>().is_ok()||kwargs.is_some_and(|k|k.get_item("modulus").ok().flatten().is_some_and(|m|!m.is_none())){return Ok(None);}
    let(e,g2)=expr_or_coeff_list(f)?;
    let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}
    let mut all=vec![e.clone()];all.extend(others.iter().cloned());
    let Some(selection)=symox::poly::gaussian::default_selection(&all,&generators) else{return Ok(None)};
    let mut p=AlgebraicPoly::from_expr(&e,&selection.gens,symox::poly::gaussian::field()).map_err(err)?;
    p.gaussian=Some(selection.integral);
    Ok(Some(p))
}

fn list_items(f: &Bound<'_, PyAny>) -> PyResult<Option<Vec<Expr>>> {
    if f.is_instance_of::<PyList>() || f.is_instance_of::<PyTuple>()
        || f.is_instance_of::<pyo3::types::PySet>() || f.is_instance_of::<pyo3::types::PyFrozenSet>() {
        return f.try_iter()?.map(|x| to_expr(&x?)).collect::<PyResult<Vec<_>>>().map(Some);
    }
    if let Ok(value)=f.downcast::<crate::PyExprNode>() {
        return Ok(match value.get().e.kind() {
            symox::Kind::Func(symox::Func::Undefined(head),items) if head.is_internal("Tuple")=>Some(items.to_vec()),
            _=>None,
        });
    }
    // Read the already loaded reference container through the conversion
    // protocol; this does not import or execute reference polynomial code.
    let ty=f.get_type();
    let module=crate::defining_module(ty.as_any());let name=ty.name()?.to_string();
    if (module.starts_with("sympy.") && matches!(name.as_str(),"Tuple"|"FiniteSet"|"EmptySet"))
        || (module=="symox" && matches!(name.as_str(),"FiniteSet"|"EmptySet")) {
        return f.getattr("args")?.try_iter()?.map(|x|to_expr(&x?)).collect::<PyResult<Vec<_>>>().map(Some);
    }
    Ok(None)
}

fn poly_result(py: Python<'_>, e: Expr, as_poly: bool, gens: &[Expr]) -> PyResult<PyObject> {
    if as_poly {
        let p = Poly::from_expr(&e, gens).map_err(err)?;
        return poly_obj(py, p);
    }
    Ok(w(e).into_pyobject(py)?.into_any().unbind())
}

/// Explicit generators and no options: the reference `sring` and
/// `parallel_poly_from_expr` select `EX` when `construct_domain` does.
fn expression_pair(f:&Bound<'_,PyAny>,g:&Bound<'_,PyAny>,gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Option<(ExpressionPoly,ExpressionPoly)>>{
    if gens.is_empty()||kwargs.is_some_and(|k|!k.is_empty())||f.downcast::<PyPoly>().is_ok()||g.downcast::<PyPoly>().is_ok()||list_items(f)?.is_some(){return Ok(None);}
    let gv=gens_of(gens)?;
    let exprs=[to_expr(f)?,to_expr(g)?];
    if exprs.iter().any(rational_function_blocker){return Ok(None);}
    Ok(symox::poly::expression_domain::joint_selection(&exprs,&gv).map(|mut v|{let b=v.pop().unwrap();(v.pop().unwrap(),b)}))
}

/// The reference `parallel_poly_from_expr` of expressions over explicit
/// generators when `_construct_composite` selects a fraction field, as for
/// the coefficients `1/a` or `1/pi`; `None` keeps the other constructors.
fn fraction_polys(items:&[&Bound<'_,PyAny>],gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Option<Vec<FractionPoly>>>{
    if gens.is_empty()||kwargs.is_some_and(|k|!k.is_empty()){return Ok(None);}
    for x in items{if x.downcast::<PyPoly>().is_ok()||list_items(x)?.is_some(){return Ok(None);}}
    let gv=gens_of(gens)?;
    let exprs=items.iter().map(|x|to_expr(x)).collect::<PyResult<Vec<_>>>()?;
    if exprs.iter().any(rational_function_blocker){return Ok(None);}
    Ok(fraction_exprs(&exprs,&gv))
}
fn fraction_exprs(exprs:&[Expr],gens:&[Expr])->Option<Vec<FractionPoly>>{
    let mut coefficients=Vec::new();
    for e in exprs{coefficients.extend(symox::poly::expression_domain::dict_from_expr(e,gens).ok()?.into_iter().map(|(_,c)|c));}
    let parameters=symox::poly::expression_domain::composite_fraction_parameters(&coefficients).ok()??;
    let display=symox::poly::fraction::domain_name(&parameters).replacen("QQ","ZZ",1);
    exprs.iter().map(|e|FractionPoly::from_expr(e,gens,&parameters,&display).ok()).collect()
}

/// The polynomials of an expression-level operation whose only keyword
/// option is `domain=EX`, over the given generators, as the reference
/// `parallel_poly_from_expr`; the flag is whether a Poly was given.
fn explicit_expression_polys(items:&[&Bound<'_,PyAny>],gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>,allowed:&[&str])->PyResult<Option<(Vec<ExpressionPoly>,bool)>>{
    let Some(kwargs)=kwargs else{return Ok(None)};
    let domain=kwargs.get_item("domain")?;
    if !domains::selects_expression(domain.as_ref()){return Ok(None);}
    for (key,_) in kwargs.iter(){let key=key.extract::<String>()?;if key!="domain"&&!allowed.contains(&key.as_str()){return Ok(None);}}
    for x in items{if x.downcast::<PyPoly>().is_err()&&list_items(x)?.is_some(){return Ok(None);}}
    domains::require_expression_generators(domain.as_ref(),!gens.is_empty())?;
    let gv=gens_of(gens)?;
    let polys=items.iter().any(|x|x.downcast::<PyPoly>().is_ok());
    let v=items.iter().map(|x|ExpressionPoly::from_expr(&expr_or_coeff_list(x)?.0,&gv).map_err(err)).collect::<PyResult<Vec<_>>>()?;
    Ok(Some((v,polys)))
}

fn expression_output(py:Python<'_>,p:ExpressionPoly,polys:bool)->PyResult<PyObject>{
    if polys{return poly_obj(py,p);}
    Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())
}

/// Piecewise or boolean content leaves the polynomial construction paths.
fn rational_function_blocker(e:&Expr)->bool{
    use symox::expr::Kind;
    let mut found=false;
    symox::traverse::preorder(e,&mut |x|{if matches!(x.kind(),Kind::Piecewise(_)|Kind::Rel(..)|Kind::Bool(_)){found=true;}});
    found
}

/// The reference `construct_domain(..., extension=True)` proceeds exactly
/// as without the option when a coefficient in the generators is not an
/// algebraic number (`-sqrt(2)*y` in `x`). The options without it are
/// returned then; `None` keeps the given options.
fn inapplicable_extension<'py>(py:Python<'py>,f:&Bound<'py,PyAny>,g:Option<&Bound<'py,PyAny>>,gens:&Bound<'py,PyTuple>,kwargs:Option<&Bound<'py,PyDict>>)->PyResult<Option<Bound<'py,PyDict>>>{
    let Some(k)=kwargs else{return Ok(None)};
    let Some(e)=k.get_item("extension")? else{return Ok(None)};
    if !e.is_instance_of::<pyo3::types::PyBool>()||!e.extract::<bool>()?{return Ok(None);}
    if k.get_item("domain")?.is_some_and(|d|!d.is_none())||division::gaussian_option(k)?.is_some()
        ||k.get_item("modulus")?.is_some_and(|m|!m.is_none()){return Ok(None);}
    let mut generators=gens_of(gens)?;
    let mut exprs=Vec::new();
    if let Some(items)=if f.downcast::<PyPoly>().is_err(){list_items(f)?}else{None} {
        if let Some(g)=g{generators.insert(0,to_expr(g)?);}
        exprs=items;
    } else {
        for x in [Some(f),g].into_iter().flatten() {
            let(e,g2)=expr_or_coeff_list(x)?;
            if gens.is_empty()&&x.downcast::<PyPoly>().is_ok(){for v in g2{if !generators.contains(&v){generators.push(v);}}}
            let(n,d)=api::as_numer_denom(&e);
            exprs.push(n);exprs.push(d);
        }
    }
    if generators.is_empty(){return Ok(None);}
    if !matches!(symox::poly::algebraic::infer_selection_in(&exprs,&generators).map_err(err)?,symox::poly::algebraic::ExtensionSelection::Unchanged){return Ok(None);}
    let copy=PyDict::new(py);
    for(key,value)in k.iter(){if key.extract::<String>().ok().as_deref()!=Some("extension"){copy.set_item(key,value)?;}}
    Ok(Some(copy))
}

/// Compute the greatest common divisor of expressions or polynomials.
#[pyfunction]
#[pyo3(signature = (f, g=None, *gens, modulus=None, **_kwargs))]
fn gcd(py: Python<'_>, f: &Bound<'_, PyAny>, g: Option<&Bound<'_, PyAny>>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let stripped=inapplicable_extension(py,f,g,gens,_kwargs)?;let _kwargs=stripped.as_ref().or(_kwargs);
    if modulus.is_none() {if let (Ok(p),Some(q))=(f.downcast::<PyPoly>(),g){
        let special=|v:&PolyValue|v.algebraic.is_some()||v.field.is_some();
        if special(&p.get().p)||q.downcast::<PyPoly>().is_ok_and(|q|special(&q.get().p)){return poly_obj(py,p.get().gcd(q)?.p);}
    }}
    if let (Some(q),None)=(g,modulus) {if let Some((v,polys))=explicit_expression_polys(&[f,q],gens,_kwargs,&[])? {
        return expression_output(py,v[0].gcd(&v[1]).map_err(err)?,polys);
    }
    if let Some(v)=fraction_polys(&[f,q],gens,_kwargs)? {
        return Ok(w(v[0].gcd(&v[1]).as_expr()).into_pyobject(py)?.into_any().unbind());
    }}
    if let Some(result)=algebraic_option_fold(py,f,g,gens,_kwargs,modulus,AlgebraicPoly::gcd)? {return Ok(result);}
    let sequence=f.downcast::<PyPoly>().is_err()&&list_items(f)?.is_some();
    let others=match g{Some(q) if f.downcast::<PyPoly>().is_err()&&!sequence=>vec![to_expr(q)?],_=>Vec::new()};
    if let Some(p)=if sequence{None}else{algebraic_option_polys(f,&others,gens,_kwargs)?} {
        let q=g.ok_or_else(||PyTypeError::new_err("expected a second polynomial"))?;
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        let q=p.convert(&to_expr(q)?).map_err(err)?;
        let result=p.gcd(&q);
        return Ok(w(result.as_expr()).into_pyobject(py)?.into_any().unbind());
    }
    if let (Some(m), Some(g)) = (modulus, g) {
        let r = api::gcd_mod(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?, check_mod(m)?).map_err(err)?;
        return Ok(w(r).into_pyobject(py)?.into_any().unbind());
    }
    if let Some(q)=g {if let Some((a,b))=expression_pair(f,q,gens,_kwargs)? {
        return Ok(w(a.gcd(&b).map_err(err)?.as_expr()).into_pyobject(py)?.into_any().unbind());
    }}
    let as_poly = f.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if let Some(items) = list_items(f)? {
        if let Some(g) = g {
            gv.insert(0, to_expr(g)?);
        }
        let mut items=items.into_iter();
        let mut acc = items.next().unwrap_or_else(Expr::zero);
        for x in items {
            acc = api::gcd(&acc, &x, &gv).map_err(err)?;
        }
        return Ok(w(acc).into_pyobject(py)?.into_any().unbind());
    }
    let g = g.ok_or_else(|| PyTypeError::new_err("gcd() takes 2 arguments or a sequence"))?;
    let a = to_expr(f)?;
    let b = to_expr(g)?;
    let r = api::gcd(&a, &b, &gv).map_err(err)?;
    let gens = if as_poly { f.downcast::<PyPoly>().unwrap().get().p.true_gens().to_vec() } else { gv };
    poly_result(py, r, as_poly, &gens)
}

#[pyfunction]
#[pyo3(signature = (f, g=None, *gens, modulus=None, **_kwargs))]
fn lcm(py: Python<'_>, f: &Bound<'_, PyAny>, g: Option<&Bound<'_, PyAny>>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let stripped=inapplicable_extension(py,f,g,gens,_kwargs)?;let _kwargs=stripped.as_ref().or(_kwargs);
    if modulus.is_none() {if let (Ok(p),Some(q))=(f.downcast::<PyPoly>(),g){
        let special=|v:&PolyValue|v.algebraic.is_some()||v.field.is_some();
        if special(&p.get().p)||q.downcast::<PyPoly>().is_ok_and(|q|special(&q.get().p)){return poly_obj(py,p.get().lcm(q)?.p);}
    }}
    if let (Some(q),None)=(g,modulus) {if let Some((v,polys))=explicit_expression_polys(&[f,q],gens,_kwargs,&[])? {
        return expression_output(py,v[0].lcm(&v[1]).map_err(err)?,polys);
    }
    if let Some(v)=fraction_polys(&[f,q],gens,_kwargs)? {
        return Ok(w(v[0].lcm(&v[1]).as_expr()).into_pyobject(py)?.into_any().unbind());
    }}
    if let Some(result)=algebraic_option_fold(py,f,g,gens,_kwargs,modulus,AlgebraicPoly::lcm)? {return Ok(result);}
    let sequence=f.downcast::<PyPoly>().is_err()&&list_items(f)?.is_some();
    let others=match g{Some(q) if f.downcast::<PyPoly>().is_err()&&!sequence=>vec![to_expr(q)?],_=>Vec::new()};
    if let Some(p)=if sequence{None}else{algebraic_option_polys(f,&others,gens,_kwargs)?} {
        let q=g.ok_or_else(||PyTypeError::new_err("expected a second polynomial"))?;
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        let q=p.convert(&to_expr(q)?).map_err(err)?;
        let result=p.lcm(&q);
        return Ok(w(result.as_expr()).into_pyobject(py)?.into_any().unbind());
    }
    if let (Some(m), Some(g)) = (modulus, g) {
        let r = api::lcm_mod(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?, check_mod(m)?).map_err(err)?;
        return Ok(w(r).into_pyobject(py)?.into_any().unbind());
    }
    if let Some(q)=g {if let Some((a,b))=expression_pair(f,q,gens,_kwargs)? {
        return Ok(w(a.lcm(&b).map_err(err)?.as_expr()).into_pyobject(py)?.into_any().unbind());
    }}
    let as_poly = f.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if let Some(items) = list_items(f)? {
        if let Some(g) = g {
            gv.insert(0, to_expr(g)?);
        }
        let mut items=items.into_iter();
        let mut acc = items.next().unwrap_or_else(Expr::one);
        for x in items {
            acc = api::lcm(&acc, &x, &gv).map_err(err)?;
        }
        return Ok(w(acc).into_pyobject(py)?.into_any().unbind());
    }
    let g = g.ok_or_else(|| PyTypeError::new_err("lcm() takes 2 arguments or a sequence"))?;
    let r = api::lcm(&to_expr(f)?, &to_expr(g)?, &gv).map_err(err)?;
    let gens = if as_poly { f.downcast::<PyPoly>().unwrap().get().p.true_gens().to_vec() } else { gv };
    poly_result(py, r, as_poly, &gens)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn cofactors<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if f.downcast::<PyPoly>().is_err(){if let Some(p)=algebraic_option_polys(f,&[to_expr(g)?],gens,_kwargs)? {
        let q=p.convert(&to_expr(g)?).map_err(err)?;
        let(h,a,b)=p.cofactors(&q);
        return PyTuple::new(py,[w(h.as_expr()),w(a.as_expr()),w(b.as_expr())]);
    }}
    let (h, a, b) = api::cofactors(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?;
    PyTuple::new(py, [w(h), w(a), w(b)])
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn cancel(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let py=f.py();
    let stripped=if f.downcast::<PyTuple>().is_err(){inapplicable_extension(py,f,None,gens,_kwargs)?}else{None};let _kwargs=stripped.as_ref().or(_kwargs);
    if f.downcast::<PyPoly>().is_err()&&f.downcast::<PyTuple>().is_err()&&domains::selects_expression(_kwargs.map(|k|k.get_item("domain")).transpose()?.flatten().as_ref())&&_kwargs.is_some_and(|k|k.len()==1){
        domains::require_expression_generators(Some(&_kwargs.unwrap().get_item("domain")?.unwrap()),!gens.is_empty())?;
        let e=symox::simptools::misc::signsimp(&to_expr(f)?,true);
        if matches!(e.kind(),symox::expr::Kind::Num(_)){return Ok(w(e).into_pyobject(py)?.into_any().unbind());}
        let e=symox::simptools::exprtools::factor_terms(&e,true,false,false,true);
        let(n,d)=api::as_numer_denom(&e);
        let gv=gens_of(gens)?;
        let p=ExpressionPoly::from_expr(&n,&gv).map_err(err)?;let q=ExpressionPoly::from_expr(&d,&gv).map_err(err)?;
        let(p,q)=p.ring_cancel(&q).map_err(err)?;
        return Ok(w(symox::ops::div(&p.as_expr(),&q.as_expr())).into_pyobject(py)?.into_any().unbind());
    }
    if let Some(result)=an_rational::cancel_input(py,f,gens,_kwargs)?{return Ok(result);}
    if let Ok(pair)=f.downcast::<PyTuple>() {
        if pair.len()!=2{return Err(PyValueError::new_err("cancel expects a numerator/denominator pair"));}
        let result=api::cancel(&symox::ops::div(&to_expr(&pair.get_item(0)?)?,&to_expr(&pair.get_item(1)?)?)).map_err(err)?;
        let(n,d)=api::as_numer_denom(&result);
        return Ok(PyTuple::new(py,[w(Expr::one()),w(n),w(d)])?.into_any().unbind());
    }
    if !gens.is_empty()&&_kwargs.is_none_or(|k|k.is_empty())&&f.downcast::<PyPoly>().is_err(){
        let e=symox::simptools::misc::signsimp(&to_expr(f)?,true);
        if !matches!(e.kind(),symox::expr::Kind::Num(_))&&!rational_function_blocker(&e){
            let e=symox::simptools::exprtools::factor_terms(&e,true,false,false,true);
            let(n,d)=api::as_numer_denom(&e);
            if let Some(mut v)=symox::poly::expression_domain::joint_selection(&[n,d],&gens_of(gens)?){
                let q=v.pop().unwrap();let p=v.pop().unwrap();
                let(p,q)=p.ring_cancel(&q).map_err(err)?;
                return Ok(w(symox::ops::div(&p.as_expr(),&q.as_expr())).into_pyobject(py)?.into_any().unbind());
            }
        }
    }
    Ok(w(api::cancel(&to_expr(f)?).map_err(err)?).into_pyobject(py)?.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, x=None, full=false, **_kwargs))]
fn apart(f: &Bound<'_, PyAny>, x: Option<&Bound<'_, PyAny>>, full: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    if let Some(result)=an_rational::apart_input(f,x,full,_kwargs)?{return Ok(result);}
    let x = match x {
        Some(x) => Some(to_expr(x)?),
        None => None,
    };
    Ok(w(api::apart(&to_expr(f)?, x.as_ref()).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature=(expr,syms,func=None,evaluate=None,exact=Some(false),distribute_order_term=true,**_kwargs))]
fn collect(expr:&Bound<'_,PyAny>,syms:&Bound<'_,PyAny>,func:Option<&Bound<'_,PyAny>>,evaluate:Option<bool>,exact:Option<bool>,distribute_order_term:bool,_kwargs:Option<&Bound<'_,PyDict>>)->PyResult<PyObject>{
    use symox::poly::collect::{self,Collection,Error,Options};
    if let Some(kwargs)=_kwargs{if let Some((key,_))=kwargs.iter().next(){return Err(PyTypeError::new_err(format!("collect() got an unexpected keyword argument '{}'",key.extract::<String>()?)));}}
    let e=to_expr(expr)?;
    let patterns=if let Ok(native)=syms.downcast::<crate::PyExprNode>(){
        if let Some(items)=symox::special::tuple_items(&native.get().e){items.to_vec()}else{vec![native.get().e.clone()]}
    }else if !syms.is_instance_of::<pyo3::types::PyString>(){
        if let Ok(items)=syms.try_iter(){items.map(|value|to_expr(&value?)).collect::<PyResult<Vec<_>>>()?}else{vec![to_expr(syms)?]}
    }else{vec![to_expr(syms)?]};
    let func=func.filter(|value|!value.is_none());
    if func.is_some_and(|value|!value.is_callable()){return Err(PyTypeError::new_err("collect func must be callable"));}
    let mut callback=|value:&Expr|->PyResult<Expr>{to_expr(&func.unwrap().call1((w(value.clone()),))?)};
    let callback=if func.is_some(){Some(&mut callback as &mut dyn FnMut(&Expr)->PyResult<Expr>)}else{None};
    let result=collect::run(&e,&patterns,Options{exact:exact.unwrap_or(false),discover:exact.is_none(),distribute_order_term},evaluate.unwrap_or(true),callback).map_err(|error|match error{Error::Unsupported(message)=>PyNotImplementedError::new_err(message),Error::Callback(error)=>error,Error::Pattern(error)=>crate::wildcards::to_python_error(error)})?;
    let py=expr.py();
    match result{
        Collection::Expression(value)=>Ok(w(value).into_pyobject(py)?.into_any().unbind()),
        Collection::Dictionary(pairs)=>{
            let dictionary=PyDict::new(py);for(key,value)in pairs{dictionary.set_item(w(key),w(value))?;}Ok(dictionary.into_any().unbind())
        }
    }
}

#[pyfunction]
#[pyo3(name = "div", signature = (f, g, *gens, **_kwargs))]
fn py_div<'py>(f: &Bound<'py, PyAny>, g: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, _kwargs: Option<&Bound<'py, PyDict>>) -> PyResult<PyObject> {
    division::run(f,g,gens,_kwargs,division::Operation::Div)
}

#[pyfunction]
#[pyo3(signature = (f, G, *gens, order="lex", domain=None, modulus=None, polys=None, auto=None, extension=None, **_kwargs))]
fn reduced<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, G: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, order: &str, domain: Option<&Bound<'_,PyAny>>, modulus: Option<u64>, polys: Option<bool>, auto: Option<bool>, extension: Option<&Bound<'_,PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    let g=G;
    let gaussian=_kwargs.map(|kw|kw.get_item("gaussian")).transpose()?.flatten();
    let gaussian_domain=named_gaussian_domain(py,gaussian.as_ref(),&[("domain",domain.is_some_and(|d|!d.is_none())),("extension",extension.is_some_and(|e|!e.is_none())),("modulus",modulus.is_some())])?;
    let domain=gaussian_domain.as_ref().or(domain);
    let explicit_domain=domain.is_some_and(|d|!d.is_none());
    if domains::selects_expression(domain)&&modulus.is_none()&&extension.is_none_or(|e|e.is_none())&&_kwargs.is_none_or(|k|k.is_empty()) {
        domains::require_expression_generators(domain,!gens.is_empty())?;
        let gv=gens_of(gens)?;
        let order=symox::poly::groebner::MonomialOrder::parse(order).ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
        let mut any_poly=f.downcast::<PyPoly>().is_ok();
        let input=ExpressionPoly::from_expr(&expr_or_coeff_list(f)?.0,&gv).map_err(err)?.to_sparse();
        let mut reps=Vec::new();
        for item in g.try_iter()? {let item=item?;any_poly|=item.downcast::<PyPoly>().is_ok();reps.push(ExpressionPoly::from_expr(&expr_or_coeff_list(&item)?.0,&gv).map_err(err)?.to_sparse());}
        let(mut q,r)=symox::poly::expression_domain::guarded(||symox::poly::fieldalg::reduce(&input,&reps,order)).map_err(err)?;
        if input.is_zero(){q.clear();}
        let polys=polys.unwrap_or(any_poly);
        let object=|rep:&symox::poly::MPoly<symox::poly::expression_domain::EX>|->PyResult<PyObject>{let p=ExpressionPoly::from_sparse(rep,&gv).map_err(err)?;expression_output(py,p,polys)};
        let q=pyo3::types::PyList::new(py,q.iter().map(object).collect::<PyResult<Vec<_>>>()?)?;
        return PyTuple::new(py,[q.into_any().unbind(),object(&r)?]);
    }
    let mut selected=domains::Selection::new(domain,modulus,extension)?;
    let domain=selected.legacy.as_deref();
    let modulus=selected.modulus;
    use symox::poly::convert::{from_exprs, from_exprs_gens, to_expr_f};
    use symox::poly::groebner::{reduce_with_options, CoefficientField, MonomialOrder};
    let monomial_order = MonomialOrder::parse(order)
        .ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
    // SymPy disables automatic field extension for explicit domains unless
    // auto=True is supplied explicitly.
    let auto = auto.unwrap_or(!explicit_domain);
    if let Some(kwargs) = _kwargs {
        if let Some((key,_)) = kwargs.iter().find(|(key,_)|key.extract::<String>().map_or(true,|k|k!="gaussian")) {
            return Err(PyTypeError::new_err(format!("reduced() got an unsupported keyword argument '{}'",key.extract::<String>()?)));
        }
    }
    if domain.is_some_and(|d| d != "ZZ" && d != "QQ" && d != "RR" && d != "ZZ_I" && d != "QQ_I" && !d.starts_with("QQ(") && !d.starts_with("ZZ(") && !d.starts_with("QQ[") && !d.starts_with("ZZ[") && !d.starts_with("QQ<")) {
        return Err(PyNotImplementedError::new_err("reduced supports ZZ, QQ, RR, and prime modulus domains"));
    }
    if domain.is_some() && modulus.is_some() {
        return Err(PyValueError::new_err("specify either domain or modulus"));
    }
    let mut expressions = Vec::new();
    let mut inferred_gens = Vec::new();
    let mut inferred_domain = Domain::ZZ;
    let mut inferred_descriptor:Option<symox::poly::domains::ExactDomain>=None;
    let mut inferred_algebraic:Option<Arc<NumberField>>=None;
    let mut inferred_prime = None;
    let mut any_poly = false;
    let mut append = |item: &Bound<'_,PyAny>| -> PyResult<()> {
        if let Ok(p) = item.downcast::<PyPoly>() {
            any_poly = true;
            let p = &p.get().p;
            let actual=domains::PyDomain::of(p)?.inner;
            inferred_descriptor=Some(if let Some(previous)=inferred_descriptor.take(){domains::unify_native(&previous,&actual)?}else{actual});
            expressions.push(p.as_expr());
            if let Some(a)=&p.algebraic{inferred_algebraic=Some(if let Some(previous)=inferred_algebraic.take(){if previous.same_embedding(&a.number_field){previous}else{NumberField::from_generators(&[previous.generator.clone(),a.number_field.generator.clone()]).map_err(err)?}}else{a.number_field.clone()});}
            inferred_gens = symox::poly::polyobj::unify_gens(&inferred_gens,p.true_gens());
            inferred_domain = inferred_domain.max(p.domain);
            if let Some(m) = p.modulus {
                if inferred_prime.is_some_and(|q| q != m) {return Err(PyValueError::new_err("incompatible finite-field domains"));}
                inferred_prime = Some(m);
            }
        } else { expressions.push(to_expr(item)?); }
        Ok(())
    };
    append(f)?;
    for divisor in g.try_iter()? { append(&divisor?)?; }
    if !explicit_domain&&modulus.is_none()&&extension.is_none(){if let Some(d)=inferred_descriptor {if matches!(d,symox::poly::domains::ExactDomain::Functions{..}|symox::poly::domains::ExactDomain::Algebraic(_)|symox::poly::domains::ExactDomain::Gaussian{..}){selected=domains::Selection::from_native(d);}}}
    let domain=selected.legacy.as_deref();
    let mut gv = gens_of(gens)?;
    if gv.is_empty() && any_poly {gv = inferred_gens;}
    let algebraic=selected.algebraic_in(extension,&expressions,&gv)?.or(if domain.is_none(){inferred_algebraic}else{None});
    if let Some(number_field)=algebraic {
        if modulus.is_some() || inferred_prime.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if gv.is_empty(){gv=symox::poly::algebraic::infer_generators(&expressions,&number_field);}
        let mut prototype=AlgebraicPoly::from_expr(&Expr::zero(),&gv,number_field.clone()).map_err(err)?;
        let inputs=expressions.iter().map(|e|prototype.convert(e)).collect::<Result<Vec<_>,_>>().map_err(err)?;
        let divisors:Vec<_>=inputs.iter().skip(1).map(|p|p.rep.clone()).collect();
        if divisors.iter().any(|p|p.is_zero()){return Err(PyZeroDivisionError::new_err("polynomial division by zero"));}
        prototype.gaussian=selected.gaussian();
        if prototype.gaussian==Some(true){if let Some(message)=inputs.iter().find_map(|p|symox::poly::gaussian::integrality_error(&p.rep)){return Err(domains::CoercionFailed::new_err(message));}}
        let _scope=AlgebraicScope::new(number_field);
        let(mut q,r)=symox::poly::fieldalg::reduce(&inputs[0].rep,&divisors,monomial_order);
        if inputs[0].rep.is_zero(){q.clear();}
        if !q.iter().chain(std::iter::once(&r)).all(symox::poly::gaussian::poly_is_integral){prototype=prototype.to_field();}
        let object=|rep|->PyResult<PyObject>{let p=prototype.with_rep(rep);if polys.unwrap_or(any_poly){poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
        let qs=PyList::new(py,q.into_iter().map(object).collect::<PyResult<Vec<_>>>()?)?;
        return PyTuple::new(py,[qs.into_any(),object(r)?.into_bound(py)]);
    }
    let available=symox::poly::fraction::infer_parameters(&expressions,&[]);
    let declared=selected.parameters(&available)?;
    if gv.is_empty() {if let Some(parameters)=&declared {gv=symox::poly::fraction::infer_parameters(&expressions,parameters);}}
    for (i,g) in gv.iter().enumerate() {
        if g.is_number() || gv[..i].contains(g) {return Err(PyValueError::new_err("generators must be distinct nonnumeric expressions"));}
    }
    let converted = if gv.is_empty() {from_exprs(&expressions)} else {from_exprs_gens(&expressions,&gv)}.map_err(err)?;
    if converted.ngens == 0 {return Err(PyValueError::new_err("reduced() requires generators"));}
    if converted.polys.iter().skip(1).any(|p| p.is_zero()) {return Err(PyZeroDivisionError::new_err("polynomial division by zero"));}
    let prime = modulus.or(if domain.is_none() {inferred_prime} else {None}).map(check_mod).transpose()?;
    let return_polys = polys.unwrap_or(any_poly);
    if declared.is_some() || converted.ngens != converted.gens.len() {
        let parameters=declared.unwrap_or_else(||symox::poly::fraction::infer_parameters(&expressions,&gv));
        if parameters.is_empty() && domain.is_none() && prime.is_none() && !return_polys && auto && order == "lex" {
            let (q,r)=api::reduced(&expressions[0],&expressions[1..],&gv).map_err(err)?;
            return PyTuple::new(py,[PyList::new(py,q.into_iter().map(w))?.into_any(),w(r).into_pyobject(py)?.into_any()]);
        }
        if prime.is_some() || matches!(domain,Some("ZZ"|"QQ"|"RR")) {
            return Err(PyValueError::new_err("coefficients do not belong to the selected domain"));
        }
        if parameters.is_empty() {return Err(PyNotImplementedError::new_err("algebraic coefficient domains are not yet implemented"));}
        let has_fractions=converted.gens[converted.ngens..].iter().any(|g|matches!(g.kind(),symox::Kind::Pow(_,e) if e.as_number().is_some_and(|n|n.is_negative())));
        let display=domain.map(str::to_string).unwrap_or_else(|| {
            let base=if has_fractions || converted.polys.iter().all(|p|p.is_integral()) {"ZZ"}else{"QQ"};
            let names=parameters.iter().map(symox::printer::to_str).collect::<Vec<_>>().join(",");
            if has_fractions {format!("{base}({names})")}else{format!("{base}[{names}]")}
        });
        let mut prototype=FractionPoly::from_expr(&Expr::zero(),&gv,&parameters,&display).map_err(err)?;
        selected.apply_fraction(&mut prototype);
        let inputs=expressions.iter().map(|e|prototype.convert(e)).collect::<Result<Vec<_>,_>>().map_err(err)?;
        if display.ends_with(']') && inputs.iter().any(|p|!p.belongs_to_ring()){return Err(PyValueError::new_err("coefficient does not belong to the selected polynomial ring"));}
        let divisors:Vec<_>=inputs.iter().skip(1).map(|p|p.rep.clone()).collect();
        if divisors.iter().any(|p|p.is_zero()){return Err(PyZeroDivisionError::new_err("polynomial division by zero"));}
        let _scope=symox::poly::ratfunc::Scope::new(parameters.len());
        let(mut q,r)=symox::poly::fieldalg::reduce_filtered(&inputs[0].rep,&divisors,monomial_order,|c|auto || !display.ends_with(']') || symox::poly::fraction::coefficient_in_ring(c,display.starts_with("ZZ")));
        if inputs[0].rep.is_zero(){q.clear();}
        let promote=q.iter().chain(std::iter::once(&r)).any(|p|p.coeffs().iter().any(|c|!symox::poly::fraction::coefficient_in_ring(c,display.starts_with("ZZ"))));
        let object=|rep|->PyResult<PyObject>{let mut p=prototype.with_rep(rep);if promote {p.promote_field();}if return_polys {poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
        let qs=PyList::new(py,q.into_iter().map(object).collect::<PyResult<Vec<_>>>()?)?;
        return PyTuple::new(py,[qs.into_any(),object(r)?.into_bound(py)]);
    }
    let integral = converted.polys.iter().all(|p| p.is_integral());
    if (domain == Some("ZZ") || prime.is_some()) && (!integral || converted.float) {
        return Err(PyValueError::new_err("coefficients must be integers in the selected domain"));
    }
    let mut output_domain = match domain {
        Some("ZZ") => Domain::ZZ,
        Some("QQ") => Domain::QQ,
        Some("RR") => Domain::RR,
        _ if converted.float => Domain::RR,
        _ if !integral => Domain::QQ,
        _ => inferred_domain,
    };
    let divisors: Vec<_> = converted.polys[1..].iter().map(|p| CoefficientField {modulus:prime}.normalize(p)).collect();
    if divisors.iter().any(|p| p.is_zero()) {return Err(PyZeroDivisionError::new_err("polynomial division by zero in selected domain"));}
    let (mut q,r) = reduce_with_options(&converted.polys[0],&divisors,monomial_order,CoefficientField {modulus:prime},prime.is_none() && !auto && output_domain == Domain::ZZ);
    if converted.polys[0].is_zero() {q.clear();}
    if output_domain == Domain::ZZ && prime.is_none() && (!r.is_integral() || q.iter().any(|p| !p.is_integral())) {output_domain = Domain::QQ;}
    let object = |rep: symox::poly::MPoly<Rat>| -> PyResult<PyObject> {
        if return_polys {
            poly_obj(py,Poly {rep,gens:converted.gens.clone(),ngens:converted.ngens,domain:output_domain,modulus:prime})
        } else {Ok(w(to_expr_f(&rep,&converted.gens,output_domain == Domain::RR)).into_pyobject(py)?.into_any().unbind())}
    };
    let qs = PyList::new(py,q.into_iter().map(object).collect::<PyResult<Vec<_>>>()?)?;
    PyTuple::new(py,[qs.into_any(),object(r)?.into_bound(py)])
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn rem<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Rem)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn quo<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Quo)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn pdiv<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Pdiv)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn prem<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Prem)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn pquo<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Pquo)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn pexquo<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Pexquo)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn exquo<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{division::run(f,g,gens,_kwargs,division::Operation::Exquo)}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, includePRS=false, **_kwargs))]
fn resultant<'py>(f: &Bound<'py, PyAny>, g: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, includePRS:bool, _kwargs: Option<&Bound<'py, PyDict>>) -> PyResult<PyObject> {
    division::run(f,g,gens,_kwargs,division::Operation::Resultant(includePRS))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn discriminant<'py>(f: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, _kwargs: Option<&Bound<'py, PyDict>>) -> PyResult<PyObject> {
    division::run(f,f,gens,_kwargs,division::Operation::Discriminant)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn subresultants<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,_kwargs:Option<&Bound<'py,PyDict>>)->PyResult<PyObject>{
    division::run(f,g,gens,_kwargs,division::Operation::Subresultants)
}

#[pyfunction]
#[pyo3(signature = (f, gen=None))]
fn degree(py: Python<'_>, f: &Bound<'_, PyAny>, gen: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
    if let Ok(p) = f.downcast::<PyPoly>() {
        let p = p.get();
        let i = p.gen_index(gen)?;
        return match p.p.degree(i) {
            Some(k) => Ok(w(Expr::int(k as i64)).into_pyobject(py)?.into_any().unbind()),
            None => Ok(w(Expr::neg_oo()).into_pyobject(py)?.into_any().unbind()),
        };
    }
    let e = to_expr(f)?;
    let g = match gen {
        None => None,
        Some(g) => {
            if g.extract::<i64>().is_ok() {
                None
            } else {
                Some(to_expr(g)?)
            }
        }
    };
    let d = if let Some(g) = &g {
        if !symox::traverse::has(&e, g) {
            Some(if e.is_zero() { None } else { Some(0) })
        } else {
            None
        }
    } else {
        None
    };
    let d = match d {
        Some(d) => d,
        None => api::degree(&e, g.as_ref()).map_err(err)?,
    };
    Ok(w(match d {
        Some(k) => Expr::int(k as i64),
        None => Expr::neg_oo(),
    })
    .into_pyobject(py)?
    .into_any()
    .unbind())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, multiple=false, filter=None, **_kwargs))]
fn roots(py: Python<'_>, f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, multiple: bool, filter: Option<&str>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (e, g2) = expr_or_coeff_list(f)?;
    // SymPy's Poly reads an equation as lhs - rhs.
    let e = match e.kind() {
        symox::expr::Kind::Rel(symox::expr::RelOp::Eq, lhs, rhs) => symox::ops::sub(lhs, rhs),
        _ => e,
    };
    let gv = gens_of(gens)?;
    if f.downcast::<PyList>().is_ok()&&!gv.is_empty(){return Err(PyValueError::new_err("redundant generators given"));}
    let gv = if gv.is_empty() { g2 } else { gv };
    let mut options=proots::RadicalOptions::default();let mut strict=false;let mut predicate=None;
    if let Some(kwargs)=_kwargs{
        for (key,value) in kwargs{
            let key=key.extract::<String>()?;
            match key.as_str(){
                "cubics"=>options.cubics=value.is_truthy()?,
                "quartics"=>options.quartics=value.is_truthy()?,
                "strict"=>strict=value.is_truthy()?,
                "predicate" if !value.is_none()=>predicate=Some(value),
                _=>{},
            }
        }
    }
    if let Some(flt)=filter{if !matches!(flt,"Z"|"Q"|"R"|"I"|"C"){return Err(PyValueError::new_err(format!("Invalid filter: {flt}")));}}
    // Structured coefficient backends are solved from their expression.
    let prepared=f.downcast::<PyPoly>().ok().filter(|p|p.get().p.true_gens()==gv.as_slice()&&p.get().p.algebraic.is_none()&&p.get().p.field.is_none());
    let mut rs = proots::with_radical_options(options,||match prepared{
        Some(p)=>proots::roots_poly(&p.get().p.base),
        None=>proots::roots(&e,&gv),
    }).map_err(err)?;
    if let Some(flt) = filter {
        rs.retain(|(r, _)| r.is_zero()||match flt {
            "Z" => r.is_integer_number(),
            "Q" => r.is_rational_number(),
            "R" => r.mask() == 0 && symox::evalf::eval_c64(r).is_some_and(|z| z.im == 0.0) && !symox::traverse::has(r, &Expr::i()),
            "I" => r.mask() == 0 && symox::evalf::eval_c64(r).is_some_and(|z| z.re.abs() < 1e-300 && z.im != 0.0),
            _ => true,
        });
    }
    if let Some(predicate)=predicate{
        let mut selected=Vec::new();
        for (root,multiplicity) in rs{
            // Reference adds zero roots after its predicate/filter stage.
            if root.is_zero()||predicate.call1((w(root.clone()),))?.is_truthy()?{
                selected.push((root,multiplicity));
            }
        }
        rs=selected;
    }
    if strict{
        let degree=prepared.map(|p|p.get().p.degree(0).unwrap_or(0))
            .or_else(||Poly::from_expr(&e,&gv).ok().map(|p|p.degree(0).unwrap_or(0)));
        if let Some(degree)=degree{
            if rs.iter().map(|(_,multiplicity)|*multiplicity).sum::<u32>()<degree{
                return Err(structural::polynomial_error(py,"UnsolvableFactorError",
                    "Strict mode: some factors cannot be solved in radicals, so a complete list of solutions cannot be returned. Call roots with strict=False to get solutions expressible in radicals (if there are any).".to_string()));
            }
        }
    }
    if multiple {
        // SymPy: ordered(result)
        rs.sort_by_cached_key(|(r, _)| symox::sortkey::sort_key(r));
        let l = PyList::empty(py);
        for (r, k) in rs {
            for _ in 0..k {
                l.append(w(r.clone()))?;
            }
        }
        return Ok(l.into_any().unbind());
    }
    let d = PyDict::new(py);
    for (r, k) in rs {
        d.set_item(w(r), k)?;
    }
    Ok(d.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, n=15, maxsteps=50, cleanup=true))]
fn nroots(py: Python<'_>, f: &Bound<'_, PyAny>, n: u32, maxsteps: u32, cleanup: bool) -> PyResult<Vec<PyExpr>> {
    if let Ok(p) = f.downcast::<PyPoly>() {
        return p.get().nroots(py, n, maxsteps, cleanup);
    }
    if f.downcast::<PyList>().is_err() {
        // The reference builds `Poly(f, greedy=False)`: algebraic and
        // transcendental constants stay in the coefficient domain.
        let p = root_analysis::numeric_root_poly(f)?;
        return p.downcast::<PyPoly>()?.get().nroots(py, n, maxsteps, cleanup);
    }
    let (e, g2) = expr_or_coeff_list(f)?;
    Ok(proots::nroots_options(&e, &g2, n, maxsteps, cleanup)
        .map_err(|error| nroots_error(py,error))?.into_iter().map(w).collect())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn real_roots(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
    let (e, g2) = expr_or_coeff_list(f)?;
    let gv = gens_of(gens)?;
    let gv = if gv.is_empty() { g2 } else { gv };
    let rs = proots::roots(&e, &gv).map_err(err)?;
    let mut out: Vec<(f64, Expr, u32)> = Vec::new();
    for (r, k) in rs {
        if r.mask() == 0 && !symox::traverse::has(&r, &Expr::i()) {
            if let Some(v) = symox::evalf::eval_f64(&r) {
                out.push((v, r, k));
            }
        }
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut v = Vec::new();
    for (_, r, k) in out {
        for _ in 0..k {
            v.push(w(r.clone()));
        }
    }
    Ok(v)
}

/// The generators of the reference `_construct_composite` domain of these
/// coefficients and whether it is a fraction field, or `None` when the
/// reference selects a simple domain or `EX`.
#[pyfunction]
#[pyo3(name = "_composite_generators")]
fn composite_generators(values: Vec<Bound<'_, PyAny>>) -> PyResult<Option<(Vec<PyExpr>, bool)>> {
    let coefficients=values.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    Ok(symox::poly::expression_domain::composite_parameters(&coefficients).map_err(err)?.map(|(gens,fractions)|(gens.into_iter().map(w).collect(),fractions)))
}

#[pyfunction]
#[pyo3(name = "_poly_together", signature = (f, deep=false, **_kwargs))]
fn py_together(f: &Bound<'_, PyAny>, deep: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(w(api::together(&to_expr(f)?)))
}

/// The reference functional `gcdex`, `half_gcdex` and `invert` over Poly
/// objects, or over explicit generators with symbolic coefficients or
/// domain options: `parallel_poly_from_expr`, then the Poly method with
/// `opt.auto` (false when a domain is given), returning Poly objects when
/// one was given unless `polys=False`.
fn functional_euclid<'py>(f: &Bound<'py, PyAny>, g: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, kwargs: Option<&Bound<'py, PyDict>>, method: &str) -> PyResult<Option<PyObject>> {
    let py=f.py();
    let option=|key:&str|->PyResult<Option<Bound<'py,PyAny>>>{Ok(match kwargs{Some(k)=>k.get_item(key)?,None=>None})};
    if let Some(k)=kwargs{for (key,_) in k.iter(){if !matches!(key.extract::<String>()?.as_str(),"auto"|"polys"|"domain"){return Ok(None);}}}
    let items=[f,g];
    let given=items.iter().filter(|x|x.downcast::<PyPoly>().is_ok()).count();
    let domain=option("domain")?;
    let mut inferred:Vec<Expr>=Vec::new();
    if given==0&&gens.is_empty(){
        let exprs=[to_expr(f)?,to_expr(g)?];
        let Ok(converted)=symox::poly::convert::from_exprs(&exprs) else{return Ok(None)};
        if converted.gens.is_empty(){
            // PolificationFailed: the reference computes in the constructed
            // ground domain of the two numbers.
            let constructed=py.import("symox.polys.constructor")?.getattr("construct_domain")?.call1((PyList::new(py,[f,g])?,))?;
            let (domain,values)=(constructed.get_item(0)?,constructed.get_item(1)?);
            let result=match domain.call_method1(method,(values.get_item(0)?,values.get_item(1)?)){
                Err(e) if e.is_instance_of::<pyo3::exceptions::PyNotImplementedError>(py)=>{
                    let failed=py.import("symox")?.getattr("ComputationFailed")?.call1((method,2,"PolificationFailed"))?;
                    return Err(PyErr::from_value(failed));
                }
                r=>r?,
            };
            if method=="invert"{return Ok(Some(domain.call_method1("to_sympy",(result,))?.unbind()));}
            let values=result.try_iter()?.map(|v|domain.call_method1("to_sympy",(v?,))).collect::<PyResult<Vec<_>>>()?;
            return Ok(Some(PyTuple::new(py,values)?.into_any().unbind()));
        }
        if converted.gens.iter().all(|x|x.is_symbol()){return Ok(None);}
        inferred=converted.gens;
    }
    if given<2&&inferred.is_empty(){
        if gens.is_empty(){return Ok(None);}
        let gv=gens_of(gens)?;
        let mut symbolic=false;
        for x in items{let e=to_expr(x)?;symbolic|=symox::traverse::free_symbols(&e).iter().any(|s|!gv.contains(s));}
        if given==0&&!symbolic&&domain.is_none()&&option("auto")?.is_none()&&option("polys")?.is_none(){return Ok(None);}
    }
    let class=py.import("symox")?.getattr("Poly")?;
    let build=|x:&Bound<'py,PyAny>|->PyResult<Bound<'py,PyAny>>{
        if x.downcast::<PyPoly>().is_ok()&&given==2&&domain.is_none(){return Ok(x.clone());}
        let source=if x.downcast::<PyPoly>().is_ok(){x.call_method0("as_expr")?}else{x.clone()};
        let mut args=vec![source];
        if !inferred.is_empty(){for y in &inferred{args.push(w(y.clone()).into_pyobject(py)?.into_any());}}
        else if gens.is_empty(){for y in x.getattr("gens")?.try_iter()?{args.push(y?);}}else{for y in gens.iter(){args.push(y);}}
        let options=PyDict::new(py);
        if let Some(d)=&domain{options.set_item("domain",d)?;}
        class.call(PyTuple::new(py,args)?,Some(&options))
    };
    let (a,b)=(build(f)?,build(g)?);
    let auto=match option("auto")?{Some(v)=>v.is_truthy()?,None=>domain.is_none()};
    let polys=match option("polys")?{Some(v)=>v.is_truthy()?,None=>given>0};
    let options=PyDict::new(py);options.set_item("auto",auto)?;
    let result=a.call_method(method,(b,),Some(&options))?;
    if polys{return Ok(Some(result.unbind()));}
    if method=="invert"{return Ok(Some(result.call_method0("as_expr")?.unbind()));}
    let values=result.try_iter()?.map(|p|p?.call_method0("as_expr")).collect::<PyResult<Vec<_>>>()?;
    Ok(Some(PyTuple::new(py,values)?.into_any().unbind()))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn gcdex(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if let Some((v,polys))=explicit_expression_polys(&[f,g],gens,_kwargs,&[])? {
        let (a,b,h)=v[0].gcdex(&v[1]).map_err(err)?;
        return Ok(PyTuple::new(py,[expression_output(py,a,polys)?,expression_output(py,b,polys)?,expression_output(py,h,polys)?])?.into_any().unbind());
    }
    if let Some(v)=fraction_polys(&[f,g],gens,_kwargs)? {
        let values=an_rational::poly_bezout(&v[0].clone().into(),&v[1].clone().into(),false,true)?;
        return Ok(PyTuple::new(py,values.into_iter().map(|p|w(p.as_expr())))?.into_any().unbind());
    }
    if let Some(selected)=an_rational::select(f,g,gens,_kwargs)?{return an_rational::bezout(py,selected,false);}
    if let Some(v)=functional_euclid(f,g,gens,_kwargs,"gcdex")?{return Ok(v);}
    let (s,t,h)=api::gcdex(&to_expr(f)?,&to_expr(g)?,&gens_of(gens)?).map_err(err)?;
    Ok(PyTuple::new(py,[w(s),w(t),w(h)])?.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn half_gcdex(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if let Some((v,polys))=explicit_expression_polys(&[f,g],gens,_kwargs,&[])? {
        let (a,h)=v[0].half_gcdex(&v[1]).map_err(err)?;
        return Ok(PyTuple::new(py,[expression_output(py,a,polys)?,expression_output(py,h,polys)?])?.into_any().unbind());
    }
    if let Some(v)=fraction_polys(&[f,g],gens,_kwargs)? {
        let values=an_rational::poly_bezout(&v[0].clone().into(),&v[1].clone().into(),true,true)?;
        return Ok(PyTuple::new(py,values.into_iter().map(|p|w(p.as_expr())))?.into_any().unbind());
    }
    if let Some(selected)=an_rational::select(f,g,gens,_kwargs)?{return an_rational::bezout(py,selected,true);}
    if let Some(v)=functional_euclid(f,g,gens,_kwargs,"half_gcdex")?{return Ok(v);}
    let (s,h)=api::half_gcdex(&to_expr(f)?,&to_expr(g)?,&gens_of(gens)?).map_err(err)?;
    Ok(PyTuple::new(py,[w(s),w(h)])?.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn invert(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let py=f.py();
    if let Some((v,polys))=explicit_expression_polys(&[f,g],gens,_kwargs,&[])? {
        return expression_output(py,v[0].invert(&v[1]).map_err(err)?,polys);
    }
    if let Some(v)=fraction_polys(&[f,g],gens,_kwargs)? {
        return Ok(w(an_rational::poly_inverse(&v[0].clone().into(),&v[1].clone().into(),true)?.as_expr()).into_pyobject(py)?.into_any().unbind());
    }
    if let Some(selected)=an_rational::select(f,g,gens,_kwargs)?{return an_rational::inverse(py,selected);}
    if let Some(v)=functional_euclid(f,g,gens,_kwargs,"invert")?{return Ok(v);}
    Ok(w(api::invert(&to_expr(f)?,&to_expr(g)?,&gens_of(gens)?).map_err(err)?).into_pyobject(py)?.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn terms_gcd(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let orig=to_expr(f)?;
    let atomic=orig.args().is_empty()||matches!(orig.kind(),symox::expr::Kind::Rel(..));
    let domain=_kwargs.map(|k|k.get_item("domain")).transpose()?.flatten();
    if domains::selects_expression(domain.as_ref())&&_kwargs.is_some_and(|k|k.len()==1)&&f.downcast::<PyPoly>().is_err(){
        if atomic{return Ok(w(orig));}
        domains::require_expression_generators(domain.as_ref(),!gens.is_empty())?;
        let p=ExpressionPoly::from_expr(&orig,&gens_of(gens)?).map_err(err)?;
        let v=[p];
        // The reference over EX (a ring and a field): clear the
        // denominators, then take the primitive part.
        let (j,p)=v[0].terms_gcd();
        let (denominator,p)=p.clear_denoms().map_err(err)?;
        let (content,p)=p.primitive().map_err(err)?;
        let coeff=symox::ops::div(&content,&denominator);
        let term=symox::ops::mul(&p.gens.iter().zip(&j).map(|(x,k)|x.powi(*k as i64)).collect::<Vec<_>>());
        let coeff=if coeff.is_one(){if term.is_one(){return Ok(w(orig));}Expr::one()}else{coeff};
        return Ok(w(api::keep_coeff(&coeff,&symox::ops::mul(&[term,p.as_expr()]))));
    }
    if !atomic {if let Some(v)=fraction_polys(&[f],gens,_kwargs)? {
        // Over `ZZ(...)` the reference also clears the denominators and
        // takes the primitive part, keeping their quotient as the coefficient.
        let p=&v[0];
        if p.rep.is_zero(){return Ok(w(orig));}
        let n=p.gens.len();
        let j:Vec<u32>=(0..n).map(|i|p.rep.iter().map(|(e,_)|e[i]).min().unwrap_or(0)).collect();
        let p=p.with_rep(symox::poly::mpoly::MPoly::from_terms(n,p.rep.iter().map(|(e,c)|(e.iter().zip(&j).map(|(a,b)|a-b).collect(),c.clone())).collect()));
        let (content,p)=p.primitive();
        let coeff=symox::poly::fraction::coefficient_expr(&content,&p.parameters);
        let term=symox::ops::mul(&p.gens.iter().zip(&j).map(|(x,k)|x.powi(*k as i64)).collect::<Vec<_>>());
        if coeff.is_one()&&term.is_one(){return Ok(w(orig));}
        return Ok(w(api::keep_coeff(&coeff,&symox::ops::mul(&[term,p.as_expr()]))));
    }}
    Ok(w(api::terms_gcd(&to_expr(f)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn compose(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    if let Some(v)=fraction_polys(&[f,g],gens,_kwargs)? {return Ok(w(v[0].compose(&v[1]).as_expr()));}
    Ok(w(api::compose(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn decompose(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
    Ok(api::decompose(&to_expr(f)?, &gens_of(gens)?).map_err(err)?.into_iter().map(w).collect())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn sqf_part<'py>(f: &Bound<'py, PyAny>, gens: &Bound<'py, PyTuple>, _kwargs: Option<&Bound<'py, PyDict>>) -> PyResult<PyObject> {
    checked_factor_views::global_part(f,gens,_kwargs)
}

fn minimal_polynomial_value(model:&symox::poly::algebraic_model::PrimitiveModel,generator:Expr)->Poly {
    let(_,coefficients)=symox::poly::dense::clear_denoms(&model.minimal);
    Poly {rep:symox::poly::dense::to_mpoly(&coefficients,0,1).to_rat(),gens:vec![generator],ngens:1,domain:Domain::QQ,modulus:None}
}
fn minpoly_error(py:Python<'_>,e:symox::poly::minpoly::MinpolyError)->PyErr {
    use symox::poly::minpoly::MinpolyError as M;
    match &e {
        M::NotAlgebraic(_)=>structural::polynomial_error(py,"NotAlgebraic",e.to_string()),
        M::ZeroDivision(_)=>PyZeroDivisionError::new_err(e.to_string()),
        M::Multiple(_)=>PyNotImplementedError::new_err(e.to_string()),
    }
}
/// The ground field of `minimal_polynomial`: `QQ`, `QQ(syms)` or an
/// algebraic field, with the label used in `GeneratorsError`.
enum MinpolyGround {Rational{symbols:Vec<Expr>,label:String},Algebraic(Arc<NumberField>)}
fn minpoly_ground(domain:Option<&Bound<'_,PyAny>>,ex:&Expr)->PyResult<MinpolyGround> {
    let selected=domains::Selection::new(domain,None,None)?;
    if selected.native.is_none()&&selected.legacy.is_none() {
        let mut symbols=symox::traverse::free_symbols(ex);
        symbols.sort_by_key(|s|s.to_string());
        if symbols.is_empty(){return Ok(MinpolyGround::Rational{symbols,label:"QQ".into()});}
        let label=format!("QQ({})",symbols.iter().map(|s|s.to_string()).collect::<Vec<_>>().join(","));
        return Ok(MinpolyGround::Rational{symbols,label});
    }
    if let Some(d)=&selected.native {
        use symox::poly::domains::ExactDomain as D;
        match &d.inner {
            D::ZZ|D::QQ=>return Ok(MinpolyGround::Rational{symbols:vec![],label:d.inner.to_string()}),
            D::Algebraic(k)=>return Ok(MinpolyGround::Algebraic(k.clone())),
            D::Functions{parameters,..}=>return Ok(MinpolyGround::Rational{symbols:parameters.clone(),label:d.inner.to_string()}),
            _=>return Err(PyNotImplementedError::new_err(format!("minimal polynomials over {} are not implemented",d.inner))),
        }
    }
    let text=selected.legacy.clone().unwrap_or_default();
    if text=="QQ"||text=="ZZ"{return Ok(MinpolyGround::Rational{symbols:vec![],label:text});}
    let available=symox::traverse::free_symbols(ex);
    if let Some(parameters)=selected.parameters(&available)? {return Ok(MinpolyGround::Rational{symbols:parameters,label:text});}
    match selected.algebraic(None,&[])? {
        Some(k)=>Ok(MinpolyGround::Algebraic(k)),
        None=>Err(PyNotImplementedError::new_err(format!("minimal polynomials over {text} are not implemented"))),
    }
}
#[pyfunction]
#[pyo3(signature=(ex,x=None,compose=true,polys=false,domain=None))]
fn minimal_polynomial(py:Python<'_>,ex:&Bound<'_,PyAny>,x:Option<&Bound<'_,PyAny>>,compose:bool,polys:bool,domain:Option<&Bound<'_,PyAny>>)->PyResult<PyObject> {
    let value=to_expr(ex)?;
    let pure=x.is_none_or(|x|x.is_none());
    let generator=x.filter(|x|!x.is_none()).map(to_expr).transpose()?.unwrap_or_else(||Expr::symbol("_x"));
    if !generator.is_symbol(){return Err(PyValueError::new_err("minimal polynomial generator must be a symbol"));}
    let ground=minpoly_ground(domain,&value)?;
    let expression=match &ground {
        MinpolyGround::Rational{symbols,label}=>{
            if symbols.contains(&generator){return Err(structural::polynomial_error(py,"GeneratorsError",format!("the variable {generator} is an element of the ground domain {label}")));}
            if !compose&&!symbols.is_empty(){return Err(PyNotImplementedError::new_err("groebner method only works for QQ"));}
            let p=py.allow_threads(||symox::poly::minpoly::minimal_polynomial(&value,symbols)).map_err(|e|minpoly_error(py,e))?;
            symox::poly::minpoly::mpoly_expr(&p,&generator,symbols)
        }
        MinpolyGround::Algebraic(field)=>{
            if !compose{return Err(PyNotImplementedError::new_err("groebner method only works for QQ"));}
            let e=py.allow_threads(||symox::poly::embedding::minimal_polynomial_over(&value,field,&generator)).map_err(|e|match e{
                symox::poly::embedding::OverError::Minpoly(e)=>minpoly_error(py,e),
                symox::poly::embedding::OverError::Poly(e)=>err(e),
            })?;
            e
        }
    };
    if polys {
        let module=py.import("symox")?;
        let cls=module.getattr(if pure{"PurePoly"}else{"Poly"})?;
        let kwargs=PyDict::new(py);kwargs.set_item("field",true)?;
        return Ok(cls.call((w(expression),w(generator)),Some(&kwargs))?.unbind());
    }
    Ok(w(expression).into_pyobject(py)?.into_any().unbind())
}
#[pyfunction]
#[pyo3(signature=(ex,x=None,compose=true,polys=false,domain=None))]
fn minpoly(py:Python<'_>,ex:&Bound<'_,PyAny>,x:Option<&Bound<'_,PyAny>>,compose:bool,polys:bool,domain:Option<&Bound<'_,PyAny>>)->PyResult<PyObject> {minimal_polynomial(py,ex,x,compose,polys,domain)}
#[pyfunction]
#[pyo3(signature=(extension,x=None,ex=false,polys=false))]
fn primitive_element(py:Python<'_>,extension:&Bound<'_,PyAny>,x:Option<&Bound<'_,PyAny>>,ex:bool,polys:bool)->PyResult<PyObject> {
    let generators=expression_vector(extension)?;
    let(generator,model,weights)=symox::poly::algebraic_model::PrimitiveModel::build_many_with_weights(&generators).map_err(err)?;
    let variable=x.map(to_expr).transpose()?.unwrap_or_else(||Expr::symbol("_x"));
    if !variable.is_symbol(){return Err(PyValueError::new_err("primitive polynomial generator must be a symbol"));}
    let polynomial=minimal_polynomial_value(&model,variable);
    let polynomial=if polys{poly_obj(py,polynomial)?}else{w(polynomial.as_expr()).into_pyobject(py)?.into_any().unbind()};
    let weights=weights.into_pyobject(py)?.into_any().unbind();
    let mut values=vec![polynomial,weights];
    if ex {
        let representations=PyList::empty(py);
        for expression in generators {let mut coefficients=model.coefficient(&expression).map_err(err)?;if coefficients.is_empty(){coefficients.push(Rat::from_i64(0));}coefficients.reverse();representations.append(PyList::new(py,coefficients.iter().map(|c|w(rat_expr(c))))?)?;}
        values.push(representations.into_any().unbind());
    }
    // The generator is certified during construction; the returned linear
    // combination and embeddings refer to this same primitive element.
    let _=generator;
    Ok(PyTuple::new(py,values)?.into_any().unbind())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    domains::register(m)?;
    division::register(m)?;
    root_analysis::register(m)?;
    checked_factor_views::register(m)?;
    algebraic_number::register(m)?;
    an_rational::register(m)?;
    m.add_function(wrap_pyfunction!(minimal_polynomial,m)?)?;
    m.add_function(wrap_pyfunction!(minpoly,m)?)?;
    m.add_function(wrap_pyfunction!(primitive_element,m)?)?;
    m.add_class::<PyPoly>()?;
    m.add_function(wrap_pyfunction!(_rebuild_poly,m)?)?;
    m.add_function(wrap_pyfunction!(_rebuild_groebner,m)?)?;
    m.add_class::<PyGroebnerBasis>()?;
    m.add_function(wrap_pyfunction!(groebner, m)?)?;
    m.add_function(wrap_pyfunction!(is_zero_dimensional, m)?)?;
    m.add_function(wrap_pyfunction!(_solve_poly_system, m)?)?;
    m.add_function(wrap_pyfunction!(factor, m)?)?;
    m.add_function(wrap_pyfunction!(sqf, m)?)?;
    m.add_function(wrap_pyfunction!(factor_list, m)?)?;
    m.add_function(wrap_pyfunction!(sqf_list, m)?)?;
    m.add_function(wrap_pyfunction!(gcd, m)?)?;
    m.add_function(wrap_pyfunction!(lcm, m)?)?;
    m.add_function(wrap_pyfunction!(cofactors, m)?)?;
    m.add_function(wrap_pyfunction!(cancel, m)?)?;
    m.add_function(wrap_pyfunction!(apart, m)?)?;
    m.add_function(wrap_pyfunction!(collect, m)?)?;
    m.add_function(wrap_pyfunction!(py_div, m)?)?;
    m.add_function(wrap_pyfunction!(rem, m)?)?;
    m.add_function(wrap_pyfunction!(reduced, m)?)?;
    m.add_function(wrap_pyfunction!(quo, m)?)?;
    m.add_function(wrap_pyfunction!(exquo, m)?)?;
    m.add_function(wrap_pyfunction!(pdiv, m)?)?;
    m.add_function(wrap_pyfunction!(prem, m)?)?;
    m.add_function(wrap_pyfunction!(pquo, m)?)?;
    m.add_function(wrap_pyfunction!(pexquo, m)?)?;
    m.add_function(wrap_pyfunction!(resultant, m)?)?;
    m.add_function(wrap_pyfunction!(discriminant, m)?)?;
    m.add_function(wrap_pyfunction!(subresultants, m)?)?;
    m.add_function(wrap_pyfunction!(degree, m)?)?;
    m.add_function(wrap_pyfunction!(roots, m)?)?;
    m.add_function(wrap_pyfunction!(nroots, m)?)?;
    m.add_function(wrap_pyfunction!(real_roots, m)?)?;
    m.add_function(wrap_pyfunction!(gcdex, m)?)?;
    m.add_function(wrap_pyfunction!(composite_generators, m)?)?;
    m.add_function(wrap_pyfunction!(half_gcdex, m)?)?;
    m.add_function(wrap_pyfunction!(invert, m)?)?;
    m.add_function(wrap_pyfunction!(terms_gcd, m)?)?;
    m.add_function(wrap_pyfunction!(compose, m)?)?;
    m.add_function(wrap_pyfunction!(decompose, m)?)?;
    m.add_function(wrap_pyfunction!(sqf_part, m)?)?;
    m.add_function(wrap_pyfunction!(py_together, m)?)?;
    m.add_function(wrap_pyfunction!(structural::_cyclotomic_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(structural::_orthopoly_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(structural::_apart_full_decomposition, m)?)?;
    m.add_function(wrap_pyfunction!(structural::_dispersion_set, m)?)?;
    m.add_function(wrap_pyfunction!(structural::_is_cyclotomic, m)?)?;
    Ok(())
}
