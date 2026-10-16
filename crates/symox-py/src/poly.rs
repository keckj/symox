//! Python bindings of the polynomial module: the `Poly` class and the
//! functions `factor`, `factor_list`, `sqf_list`, `gcd`, `lcm`, `cancel`,
//! `apart`, `collect`, `div`, `rem`, `quo`, `resultant`, `discriminant`,
//! `roots`, `nroots`, `degree`, ...

#![allow(non_snake_case)]

use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::expr::Expr;
use symox::poly::api;
use symox::poly::convert::PolyError;
use symox::poly::polyobj::{Domain, Poly};
use symox::poly::fraction::{FractionPoly, Polynomial as PolyValue};
use symox::poly::algebraic::{AlgebraicPoly,NumberField,AlgebraicScope};
use std::sync::Arc;
use symox::poly::roots as proots;
use symox::poly::Rat;

use crate::{to_expr, PyExpr};

fn w(e: Expr) -> PyExpr {
    PyExpr { e }
}

fn err(e: PolyError) -> PyErr {
    match &e {
        PolyError::Unsupported(s) if s.contains("division by zero") => PyZeroDivisionError::new_err(s.clone()),
        _ => PyValueError::new_err(e.to_string()),
    }
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

fn algebraic_extension(extension:Option<&Bound<'_,PyAny>>,expressions:&[Expr],domain:Option<&str>)->PyResult<Option<Arc<NumberField>>> {
    if let Some(d)=domain {
        if extension.is_some(){return Err(PyValueError::new_err("specify either domain or extension"));}
        return NumberField::parse(d).map_err(err);
    }
    let Some(extension)=extension else{return Ok(None)};
    if extension.is_instance_of::<pyo3::types::PyBool>() {
        if !extension.extract::<bool>()?{return Err(PyValueError::new_err("extension=False is not a supported extension option"));}
        return symox::poly::algebraic::infer(expressions).map_err(err);
    }
    let generators:Vec<Expr>=if extension.downcast::<PyList>().is_ok()||extension.downcast::<PyTuple>().is_ok(){extension.try_iter()?.map(|e|to_expr(&e?)).collect::<PyResult<_>>()?}else{vec![to_expr(extension)?]};
    if generators.is_empty(){return Err(PyValueError::new_err("an algebraic extension needs a generator"));}
    NumberField::from_generators(&generators).map(Some).map_err(err)
}

/// A univariate polynomial given as a coefficient list (SymPy accepts
/// `roots([1, 0, -1])`).
fn expr_or_coeff_list(obj: &Bound<'_, PyAny>) -> PyResult<(Expr, Vec<Expr>)> {
    if let Ok(p) = obj.downcast::<PyPoly>() {
        let p = &p.get().p;
        return Ok((p.as_expr(), p.true_gens().to_vec()));
    }
    if let Ok(l) = obj.downcast::<PyList>() {
        let x = Expr::symbol("x");
        let n = l.len();
        let mut terms = Vec::new();
        for (i, c) in l.iter().enumerate() {
            let c = to_expr(&c)?;
            terms.push(symox::ops::mul(&[c, symox::ops::pow(&x, &Expr::int((n - 1 - i) as i64))]));
        }
        return Ok((symox::ops::add(&terms), vec![x]));
    }
    Ok((to_expr(obj)?, vec![]))
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
}

impl PyGroebnerBasis {
    fn build(
        f: &Bound<'_, PyAny>,
        gens: &Bound<'_, PyTuple>,
        order: &str,
        domain: Option<&str>,
        modulus: Option<u64>,
        method: Option<&str>,
        polys: Option<bool>,
        extension: Option<&Bound<'_,PyAny>>,
    ) -> PyResult<Self> {
        use symox::poly::convert::{from_exprs, from_exprs_gens};
        use symox::poly::groebner::{buchberger, CoefficientField, MonomialOrder};
        let monomial_order = MonomialOrder::parse(order)
            .ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
        let use_f5b=match method {
            None|Some("buchberger")=>false,
            Some("f5b")=>true,
            Some(other)=>return Err(PyValueError::new_err(format!("'{other}' is not a valid Groebner algorithm (valid are 'buchberger' and 'f5b')"))),
        };
        if domain.is_some_and(|d| d != "ZZ" && d != "QQ" && !d.starts_with("QQ(") && !d.starts_with("ZZ(") && !d.starts_with("QQ[") && !d.starts_with("ZZ[") && !d.starts_with("QQ<")) {
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
        let mut inferred_fraction:Option<FractionPoly>=None;
        for item in f.try_iter()? {
            let item = item?;
            if let Ok(p) = item.downcast::<PyPoly>() {
                any_poly = true;
                let p = &p.get().p;
                exprs.push(p.as_expr());
                if let Some(a)=&p.algebraic {
                    inferred_algebraic=Some(if let Some(previous)=inferred_algebraic {NumberField::from_generators(&[previous.generator.clone(),a.number_field.generator.clone()]).map_err(err)?}else{a.number_field.clone()});
                }
                if let Some(field)=&p.field {
                    if let Some(current)=&mut inferred_fraction {
                        for parameter in &field.parameters {
                            if !current.parameters.iter().any(|p|p.as_symbol()==parameter.as_symbol()){current.parameters.push(parameter.clone());}
                        }
                        symox::poly::convert::sort_gens(&mut current.parameters);
                        let base=if current.display.starts_with("QQ")||field.display.starts_with("QQ"){ "QQ" }else{ "ZZ" };
                        let names=current.parameters.iter().map(symox::printer::to_str).collect::<Vec<_>>().join(",");
                        current.display=if current.display.ends_with(']')&&field.display.ends_with(']'){format!("{base}[{names}]")}else{format!("{base}({names})")};
                    } else {inferred_fraction=Some(field.clone());}
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
        let algebraic=algebraic_extension(extension,&exprs,domain)?.or(if domain.is_none(){inferred_algebraic}else{None});
        if let Some(number_field)=algebraic {
            if modulus.is_some() || inferred_prime.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
            if generators.is_empty(){generators=symox::poly::fraction::infer_parameters(&exprs,&[]);}
            let prototype=AlgebraicPoly::from_expr(&Expr::zero(),&generators,number_field.clone()).map_err(err)?;
            let inputs=exprs.iter().map(|e|prototype.convert(e).map(|p|p.rep)).collect::<Result<Vec<_>,_>>().map_err(err)?;
            let _scope=AlgebraicScope::new(number_field);
            let reps=if use_f5b {symox::poly::f5b::f5b(&inputs,monomial_order).map_err(err)?}else{symox::poly::fieldalg::buchberger(&inputs,monomial_order)};
            let basis=reps.into_iter().map(|rep|prototype.with_rep(rep).into()).collect();
            return Ok(Self {basis,generators,coefficient_domain:Domain::QQ,prime:None,monomial_order,return_polys:polys.unwrap_or(any_poly),fraction_domain:None,algebraic_domain:Some(prototype)});
        }
        let mut available=symox::poly::fraction::infer_parameters(&exprs,&[]);
        if let Some(field)=&inferred_fraction {for parameter in &field.parameters {if !available.iter().any(|p|p.as_symbol()==parameter.as_symbol()){available.push(parameter.clone());}}}
        let declared=domain.map(|d|symox::poly::fraction::parse_domain(d,&available)).transpose().map_err(err)?.flatten();
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
            let parameters=declared.unwrap_or_else(||symox::poly::fraction::infer_parameters(&exprs,&generators));
            if parameters.is_empty() {return Err(PyNotImplementedError::new_err("algebraic coefficient extensions are not implemented"));}
            let has_fractions=converted.gens[converted.ngens..].iter().any(|g|matches!(g.kind(),symox::Kind::Pow(_,e) if e.as_number().is_some_and(|n|n.is_negative())));
            let display=domain.map(str::to_string).unwrap_or_else(|| {
                let base=if has_fractions || converted.polys.iter().all(|p|p.is_integral()) {"ZZ"}else{"QQ"};
                let names=parameters.iter().map(symox::printer::to_str).collect::<Vec<_>>().join(",");
                if has_fractions {format!("{base}({names})")}else{format!("{base}[{names}]")}
            });
            let prototype=FractionPoly::from_expr(&Expr::zero(),&generators,&parameters,&display).map_err(err)?;
            let converted_inputs=exprs.iter().map(|e|prototype.convert(e)).collect::<Result<Vec<_>,_>>().map_err(err)?;
            if display.ends_with(']') && converted_inputs.iter().any(|p|!p.belongs_to_ring()){return Err(PyValueError::new_err("coefficient does not belong to the selected polynomial ring"));}
            let inputs=converted_inputs.into_iter().map(|p|p.rep).collect::<Vec<_>>();
            let _scope=symox::poly::ratfunc::Scope::new(parameters.len());
            let reps=if use_f5b {symox::poly::f5b::f5b(&inputs,monomial_order).map_err(err)?}else{symox::poly::fieldalg::buchberger(&inputs,monomial_order)};
            let basis=reps.into_iter().map(|rep| {
                let p=prototype.with_rep(rep);
                if display.ends_with(']') {p.ring_normalize().into()}else{p.into()}
            }).collect();
            return Ok(Self {basis,generators,coefficient_domain:Domain::QQ,prime:None,monomial_order,return_polys:polys.unwrap_or(any_poly),fraction_domain:Some(prototype),algebraic_domain:None});
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
        domain: Option<&str>,
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
    fn domain(&self) -> String {
        if let Some(field)=&self.algebraic_domain {return field.number_field.display();}
        if let Some(field)=&self.fraction_domain {return field.display.clone();}
        if let Some(p) = self.prime {
            format!("GF({p})")
        } else if self.coefficient_domain == Domain::ZZ {
            "ZZ".into()
        } else {
            "QQ".into()
        }
    }

    /// Convert orders through the finite-dimensional quotient algebra.
    fn fglm(&self, order: &str) -> PyResult<Self> {
        use symox::poly::groebner::{fglm, CoefficientField, MonomialOrder};
        let target = MonomialOrder::parse(order)
            .ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
        if target == self.monomial_order { return Ok(self.clone()); }
        if let Some(field) = &self.algebraic_domain {
            let _scope = AlgebraicScope::new(field.number_field.clone());
            let reps: Vec<_> = self.basis.iter().map(|p| p.algebraic.as_ref().unwrap().rep.clone()).collect();
            let converted = symox::poly::fieldalg::fglm(&reps,self.monomial_order,target,self.generators.len())
                .map_err(|error| PyNotImplementedError::new_err(error.to_string()))?;
            let mut result = self.clone();
            result.monomial_order = target;
            result.basis = converted.into_iter().map(|rep| field.with_rep(rep).into()).collect();
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
    fn __repr__(&self) -> String {
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
            self.domain(),
            self.order()
        )
    }

    fn __reduce__(&self,py:Python<'_>)->PyResult<(PyObject,PyObject)> {
        let helper=py.import("symox._core")?.getattr("_rebuild_groebner")?.unbind();
        let parameters=self.fraction_domain.as_ref().map(|f|f.parameters.iter().cloned().map(w).collect::<Vec<_>>());
        let arguments=(self.exprs(),self.generators.iter().cloned().map(w).collect::<Vec<_>>(),parameters,self.domain(),self.order(),self.return_polys,self.prime).into_pyobject(py)?.into_any().unbind();
        Ok((helper,arguments))
    }

    #[pyo3(signature = (f, auto=true))]
    fn reduce(
        &self,
        py: Python<'_>,
        f: &Bound<'_, PyAny>,
        auto: bool,
    ) -> PyResult<(Vec<PyObject>, PyObject)> {
        if let Some(field)=&self.algebraic_domain {
            let input=field.convert(&to_expr(f)?).map_err(err)?;
            let _scope=AlgebraicScope::new(field.number_field.clone());
            let reps:Vec<_>=self.basis.iter().map(|p|p.algebraic.as_ref().unwrap().rep.clone()).collect();
            let(mut q,r)=symox::poly::fieldalg::reduce(&input.rep,&reps,self.monomial_order);
            if input.rep.is_zero(){q.clear();}
            let object=|rep|->PyResult<PyObject>{let p=field.with_rep(rep);if self.return_polys{poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
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
            let object=|rep|->PyResult<PyObject>{let mut p=field.with_rep(rep);if promote {p.display=p.display.replace('[',"(").replace(']',")");}if self.return_polys {poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
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

#[pyfunction]
#[pyo3(signature = (f, *gens, order="lex", domain=None, modulus=None, method=None, polys=None, extension=None))]
fn groebner(
    f: &Bound<'_, PyAny>,
    gens: &Bound<'_, PyTuple>,
    order: &str,
    domain: Option<&str>,
    modulus: Option<u64>,
    method: Option<&str>,
    polys: Option<bool>,
    extension: Option<&Bound<'_,PyAny>>,
) -> PyResult<PyGroebnerBasis> {
    PyGroebnerBasis::build(f, gens, order, domain, modulus, method, polys, extension)
}

#[pyfunction]
#[pyo3(signature = (f, *gens, order="lex", domain=None, modulus=None, method=None, extension=None))]
fn is_zero_dimensional(
    f: &Bound<'_, PyAny>,
    gens: &Bound<'_, PyTuple>,
    order: &str,
    domain: Option<&str>,
    modulus: Option<u64>,
    method: Option<&str>,
    extension: Option<&Bound<'_,PyAny>>,
) -> PyResult<bool> {
    Ok(
        PyGroebnerBasis::build(f, gens, order, domain, modulus, method, None, extension)?
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
#[pyclass(name = "Poly", module = "symox._core", frozen)]
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
        if let Some(field)=algebraic_extension(None,&[],domain)? {
            let(e,g2)=expr_or_coeff_list(expr)?;let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}
            if generators.is_empty(){generators=symox::poly::fraction::infer_parameters(std::slice::from_ref(&e),&[]);}
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
        let p=Poly::from_expr(&e,&generators).map_err(err)?;
        let fractions=p.domain_gens().iter().any(|g| matches!(g.kind(),symox::Kind::Pow(_,e) if e.as_number().is_some_and(|n|n.is_negative())) && !symox::traverse::free_symbols(g).is_empty());
        if domain.is_none() && fractions {
            let parameters=infer_parameters(std::slice::from_ref(&e),p.true_gens());
            let display=domain_name(&parameters).replacen("QQ","ZZ",1);
            return Ok(Self {p:FractionPoly::from_expr(&e,p.true_gens(),&parameters,&display).map_err(err)?.into()});
        }
        Ok(Self {p:set_domain(p,domain)?.into()})
    }

    fn other(&self, o: &Bound<'_, PyAny>) -> Option<PolyValue> {
        if let Ok(q)=o.downcast::<PyPoly>() {
            let q=&q.get().p;
            let gens=symox::poly::polyobj::unify_gens(self.p.true_gens(),q.true_gens());
            if let Some(field)=self.p.algebraic.as_ref().or(q.algebraic.as_ref()) {
                if self.p.field.is_some()||q.field.is_some()||self.p.modulus.is_some()||q.modulus.is_some(){return None;}
                let number_field=if let(Some(a),Some(b))=(&self.p.algebraic,&q.algebraic){NumberField::from_generators(&[a.number_field.generator.clone(),b.number_field.generator.clone()]).ok()?}else{field.number_field.clone()};
                AlgebraicPoly::from_expr(&self.p.as_expr(),&gens,number_field.clone()).ok()?;
                AlgebraicPoly::from_expr(&q.as_expr(),&gens,number_field).ok()?;
                return Some(q.clone());
            }
            for f in [self.p.field.as_ref(),q.field.as_ref()].into_iter().flatten() {
                if gens.iter().any(|g|f.parameters.contains(g)) || q.modulus.is_some() || self.p.modulus.is_some(){return None;}
            }
            if let Some(prototype)=self.p.field.as_ref().or(q.field.as_ref()) {
                let expressions=[self.p.as_expr(),q.as_expr()];
                let mut parameters=prototype.parameters.clone();
                for p in symox::poly::fraction::infer_parameters(&expressions,&gens){if !parameters.contains(&p){parameters.push(p);}}
                for e in &expressions {FractionPoly::from_expr(e,&gens,&parameters,&prototype.display).ok()?;}
            }
            return Some(q.clone());
        }
        if let Some(field)=&self.p.algebraic {return field.convert(&to_expr(o).ok()?).ok().map(Into::into);}
        if let Some(field)=&self.p.field {
            return field.convert(&to_expr(o).ok()?).ok().map(Into::into);
        }
        let e=to_expr(o).ok()?;
        Poly::from_expr(&e,self.p.true_gens()).ok().map(|mut q| {
            if q.domain < self.p.domain {q.domain=self.p.domain;}
            q.into()
        })
    }

    fn binop(&self, o: &Bound<'_, PyAny>, f: impl Fn(&PolyValue, &PolyValue) -> PolyValue, g: impl Fn(&Expr, &Expr) -> Expr) -> PyResult<PyObject> {
        let py = o.py();
        match self.other(o) {
            Some(q) => Ok(PyPoly { p: (f(&self.p, &q)).into() }.into_pyobject(py)?.into_any().unbind()),
            None => {
                if let Ok(q)=o.downcast::<PyPoly>() {
                    if self.p.algebraic.is_some() || q.get().p.algebraic.is_some() {
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
                if let Ok(i) = g.extract::<usize>() {
                    return Ok(i);
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
    Ok(PyPoly { p: p.into() }.into_pyobject(py)?.into_any().unbind())
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
    #[pyo3(signature = (expr, *gens, domain=None, modulus=None, extension=None, **_kwargs))]
    fn py_new(expr: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, domain: Option<&str>, modulus: Option<u64>, extension:Option<&Bound<'_,PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let(e,g2)=expr_or_coeff_list(expr)?;
        let mut r=if let Some(field)=algebraic_extension(extension,std::slice::from_ref(&e),domain)? {
            let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}if generators.is_empty(){generators=symox::poly::fraction::infer_parameters(std::slice::from_ref(&e),&[]);}
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
        let arguments=(w(self.p.as_expr()),self.p.true_gens().iter().cloned().map(w).collect::<Vec<_>>(),parameters,self.p.domain_str(),self.p.modulus).into_pyobject(py)?.into_any().unbind();
        Ok((helper,arguments))
    }

    fn _symox_(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    fn as_expr(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    #[getter]
    fn expr(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    #[getter]
    fn gens<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.p.true_gens().iter().map(|g| w(g.clone())))
    }

    #[getter]
    fn gen(&self) -> PyExpr {
        w(self.p.true_gens()[0].clone())
    }

    #[getter]
    fn domain(&self) -> String {
        match self.p.modulus {
            Some(p) => format!("GF({p})"),
            None => self.p.domain_str(),
        }
    }

    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PySet>> {
        let fs = symox::traverse::free_symbols(&self.p.as_expr());
        let mut v: Vec<PyExpr> = Vec::new();
        for g in self.p.true_gens() {
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
        let v: Vec<i64> = (0..self.p.ngens).map(|i| self.p.rep.degree(i)).collect();
        PyTuple::new(py, v)
    }

    fn total_degree(&self, py: Python<'_>) -> PyResult<PyObject> {
        degree_obj(py, self.p.total_degree())
    }

    fn coeffs(&self) -> Vec<PyExpr> {
        self.p.coeffs().into_iter().map(w).collect()
    }

    fn all_coeffs(&self) -> PyResult<Vec<PyExpr>> {
        Ok(self.p.all_coeffs().map_err(err)?.into_iter().map(w).collect())
    }

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

    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (m, c) in self.p.terms() {
            d.set_item(PyTuple::new(py, m)?, w(c))?;
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
        if let Some(f)=&self.p.algebraic {if f.rep.is_constant(){return Ok(true);}let(_,fs)=f.factor_list().map_err(err)?;return Ok(fs.len()==1 && fs[0].1==1);}
        if let Some(f)=&self.p.field {let(_,fs)=f.factor_list();return Ok(fs.len()==1 && fs[0].1==1);}
        let (_, fs) = self.p.factor_list();
        Ok(fs.len() == 1 && fs[0].1 == 1)
    }

    #[getter]
    fn is_sqf(&self) -> PyResult<bool> {
        if let Some(f)=&self.p.algebraic {if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate squarefree factorization over algebraic fields is not implemented"));}return Ok(f.sqf_list().1.iter().all(|(_,k)|*k==1));}
        if let Some(f)=&self.p.field {return Ok(f.sqf_list().1.iter().all(|(_,k)|*k==1));}
        let (_, fs) = self.p.sqf_list();
        Ok(fs.iter().all(|(_, k)| *k == 1))
    }

    fn __add__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.add(b), |a, b| symox::ops::add2(a, b))
    }
    fn __radd__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.add(a), |a, b| symox::ops::add2(b, a))
    }
    fn __sub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.sub(b), |a, b| symox::ops::sub(a, b))
    }
    fn __rsub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.sub(a), |a, b| symox::ops::sub(b, a))
    }
    fn __mul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.mul(b), |a, b| symox::ops::mul2(a, b))
    }
    fn __rmul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.mul(a), |a, b| symox::ops::mul2(b, a))
    }
    fn __neg__(&self) -> PyPoly {
        PyPoly { p: (self.p.neg()).into() }
    }
    fn __pos__(&self) -> PyPoly {
        self.clone()
    }
    fn __pow__(&self, k: u32, _m: Option<&Bound<'_, PyAny>>) -> PyPoly {
        PyPoly { p: (self.p.pow(k)).into() }
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
        self.div(py, o)
    }

    fn __eq__(&self, o: &Bound<'_, PyAny>) -> bool {
        if let Ok(q) = o.downcast::<PyPoly>() {
            let q = &q.get().p;
            return self.p.true_gens() == q.true_gens() && self.p.as_expr() == q.as_expr();
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

    fn div<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        let (a, b) = self.p.div(&q).map_err(err)?;
        PyTuple::new(py, [PyPoly { p: (a).into() }, PyPoly { p: (b).into() }])
    }

    fn rem(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        self.__mod__(o)
    }

    fn quo(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        self.__floordiv__(o)
    }

    fn exquo(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: (self.p.exquo(&q).map_err(err)?).into() })
    }

    fn gcd(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        if (self.p.algebraic.is_some()||q.algebraic.is_some()) && symox::poly::polyobj::unify_gens(self.p.true_gens(),q.true_gens()).len()!=1{return Err(PyNotImplementedError::new_err("multivariate gcd over number fields is not implemented"));}
        Ok(PyPoly { p: (self.p.gcd(&q)).into() })
    }

    fn lcm(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        if (self.p.algebraic.is_some()||q.algebraic.is_some()) && symox::poly::polyobj::unify_gens(self.p.true_gens(),q.true_gens()).len()!=1{return Err(PyNotImplementedError::new_err("multivariate gcd over number fields is not implemented"));}
        Ok(PyPoly { p: (self.p.lcm(&q)).into() })
    }

    fn cofactors<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        if (self.p.algebraic.is_some()||q.algebraic.is_some()) && symox::poly::polyobj::unify_gens(self.p.true_gens(),q.true_gens()).len()!=1{return Err(PyNotImplementedError::new_err("multivariate gcd over number fields is not implemented"));}
        let (h, a, b) = self.p.cofactors(&q);
        PyTuple::new(py, [PyPoly { p: (h).into() }, PyPoly { p: (a).into() }, PyPoly { p: (b).into() }])
    }

    fn factor_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if let Some(f)=&self.p.algebraic {
            if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic factorization is not implemented"));}
            let(c,factors)=f.factor_list().map_err(err)?;let list=PyList::empty(py);
            for(p,k)in factors {list.append(PyTuple::new(py,[poly_obj(py,p)?,k.into_pyobject(py)?.into_any().unbind()])?)?;}
            return PyTuple::new(py,[w(c.as_expr()).into_pyobject(py)?.into_any(),list.into_any()]);
        }
        if let Some(f)=&self.p.field {
            let(c,factors)=f.factor_list();let list=PyList::empty(py);
            for(p,k)in factors {list.append(PyTuple::new(py,[poly_obj(py,p)?,k.into_pyobject(py)?.into_any().unbind()])?)?;}
            let c=symox::poly::fraction::coefficient_expr(&c,&f.parameters);
            return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),list.into_any()]);
        }
        let (c, fs) = self.p.factor_list();
        let l = PyList::empty(py);
        for (f, k) in fs {
            l.append(PyTuple::new(py, [poly_obj(py, f)?, k.into_pyobject(py)?.into_any().unbind()])?)?;
        }
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), l.into_any()])
    }

    fn sqf_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if let Some(f)=&self.p.algebraic {
            if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic factorization is not implemented"));}
            let(c,factors)=f.sqf_list();let list=PyList::empty(py);
            for(p,k)in factors {list.append(PyTuple::new(py,[poly_obj(py,p)?,k.into_pyobject(py)?.into_any().unbind()])?)?;}
            return PyTuple::new(py,[w(c.as_expr()).into_pyobject(py)?.into_any(),list.into_any()]);
        }
        if let Some(f)=&self.p.field {
            let(c,factors)=f.sqf_list();let list=PyList::empty(py);
            for(p,k)in factors {list.append(PyTuple::new(py,[poly_obj(py,p)?,k.into_pyobject(py)?.into_any().unbind()])?)?;}
            let c=symox::poly::fraction::coefficient_expr(&c,&f.parameters);
            return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),list.into_any()]);
        }
        let (c, fs) = self.p.sqf_list();
        let l = PyList::empty(py);
        for (f, k) in fs {
            l.append(PyTuple::new(py, [poly_obj(py, f)?, k.into_pyobject(py)?.into_any().unbind()])?)?;
        }
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), l.into_any()])
    }

    fn primitive<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if self.p.algebraic.is_some(){return PyTuple::new(py,[w(if self.p.is_zero(){Expr::zero()}else{Expr::one()}).into_pyobject(py)?.into_any(),self.clone().into_pyobject(py)?.into_any()]);}
        if let Some(f)=&self.p.field {
            let(c,p)=f.primitive();let c=symox::poly::fraction::coefficient_expr(&c,&f.parameters);
            return PyTuple::new(py,[w(c).into_pyobject(py)?.into_any(),poly_obj(py,p)?.into_bound(py)]);
        }
        let (c, p) = self.p.primitive();
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), PyPoly { p: p.into() }.into_pyobject(py)?.into_any()])
    }

    fn content(&self) -> PyExpr {
        if self.p.algebraic.is_some(){return w(if self.p.is_zero(){Expr::zero()}else{Expr::one()});}
        if let Some(f)=&self.p.field {return w(symox::poly::fraction::coefficient_expr(&f.primitive().0,&f.parameters));}
        w(rat_expr(&self.p.primitive().0))
    }

    fn monic(&self) -> PyResult<PyPoly> {
        if let Some(f)=&self.p.algebraic {return Ok(PyPoly {p:f.monic().into()});}
        if let Some(f)=&self.p.field {return Ok(PyPoly {p:f.monic().into()});}
        if self.p.is_zero() {
            return Ok(self.clone());
        }
        let lc = self.p.rep.lc();
        let mut p = self.p.mul_ground(&lc.inv());
        if p.domain == Domain::ZZ && p.modulus.is_none() {
            p.domain = Domain::QQ;
        }
        Ok(PyPoly { p: p.into() })
    }

    #[pyo3(signature = (*gens))]
    fn diff(&self, gens: &Bound<'_, PyTuple>) -> PyResult<PyPoly> {
        let mut p = self.p.clone();
        if gens.is_empty() {
            return Ok(PyPoly { p: (p.diff(0)).into() });
        }
        for g in gens.iter() {
            let i = self.gen_index(Some(&g))?;
            p = p.diff(i);
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
        let mut pairs = Vec::new();
        for (i, a) in args.iter().enumerate() {
            pairs.push((self.p.true_gens()[i].clone(), to_expr(&a)?));
        }
        Ok(w(symox::expand::expand(&symox::subs::subs(&self.p.as_expr(), &pairs))))
    }

    #[pyo3(signature = (n=15, **_kwargs))]
    fn nroots(&self, n: u32, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
        Ok(proots::nroots_poly(&self.p, n).map_err(err)?.into_iter().map(w).collect())
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
        if let Some(f)=&self.p.algebraic {if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("algebraic compose requires one generator"));}let q=f.convert(&to_expr(g)?).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        if let Some(f)=&self.p.field {let q=f.convert(&to_expr(g)?).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        let q = self.other(g).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: (self.p.compose(&q)).into() })
    }

    fn decompose(&self) -> PyResult<Vec<PyPoly>> {
        if self.p.field.is_some()||self.p.algebraic.is_some(){return Err(PyNotImplementedError::new_err("functional decomposition over fraction fields is not yet implemented"));}
        Ok(self.p.decompose().map_err(err)?.into_iter().map(|p| PyPoly { p: p.into() }).collect())
    }

    fn shift(&self, a: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        if let Some(f)=&self.p.algebraic {if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("algebraic shift requires one generator"));}let argument=symox::ops::add2(&f.gens[0],&to_expr(a)?);let q=f.convert(&argument).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        if let Some(f)=&self.p.field {let argument=symox::ops::add2(&f.gens[0],&to_expr(a)?);let q=f.convert(&argument).map_err(err)?;return Ok(PyPoly {p:f.compose(&q).into()});}
        let e = to_expr(a)?;
        let r = e
            .as_number()
            .and_then(Rat::from_number)
            .ok_or_else(|| PyValueError::new_err("shift expects a rational number"))?;
        Ok(PyPoly { p: (self.p.shift(&r)).into() })
    }

    #[pyo3(signature = (*gens))]
    fn integrate(&self, gens: &Bound<'_, PyTuple>) -> PyResult<PyPoly> {
        let mut p = self.p.clone();
        if gens.is_empty() {
            return Ok(PyPoly { p: (p.integrate(0)).into() });
        }
        for g in gens.iter() {
            let i = self.gen_index(Some(&g))?;
            p = p.integrate(i);
        }
        Ok(PyPoly { p: p.into() })
    }

    fn clear_denoms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if self.p.algebraic.is_some(){return PyTuple::new(py,[w(Expr::one()).into_pyobject(py)?.into_any(),self.clone().into_pyobject(py)?.into_any()]);}
        if let Some(f)=&self.p.field {
            let _scope=symox::poly::ratfunc::Scope::new(f.parameters.len());
            let d=f.denominator();let e=symox::poly::convert::to_expr_int(&d,&f.parameters);
            let p=f.with_rep(f.rep.scale(&symox::poly::ratfunc::RF::from_poly(d)));
            return PyTuple::new(py,[w(e).into_pyobject(py)?.into_any(),poly_obj(py,p)?.into_bound(py)]);
        }
        let (d, p) = self.p.clear_denoms();
        PyTuple::new(py, [w(Expr::num(d.to_number())).into_pyobject(py)?.into_any(), PyPoly { p: p.into() }.into_pyobject(py)?.into_any()])
    }

    fn max_norm(&self) -> PyResult<PyExpr> {
        if self.p.field.is_some()||self.p.algebraic.is_some(){return Err(PyTypeError::new_err("rational function coefficients have no ordered absolute norm"));}
        Ok(w(rat_expr(&self.p.max_norm())))
    }

    fn l1_norm(&self) -> PyResult<PyExpr> {
        if self.p.field.is_some()||self.p.algebraic.is_some(){return Err(PyTypeError::new_err("rational function coefficients have no ordered absolute norm"));}
        Ok(w(rat_expr(&self.p.l1_norm())))
    }

    fn sqf_part(&self) -> PyResult<PyPoly> {
        if let Some(f)=&self.p.algebraic {if f.rep.is_zero(){return Ok(self.clone());}if f.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic squarefree factorization is not implemented"));}let mut p=f.convert(&Expr::one()).map_err(err)?;for(q,_)in f.sqf_list().1{p=p.mul(&q);}return Ok(PyPoly {p:p.monic().into()});}
        if let Some(f)=&self.p.field {let mut p=f.convert(&Expr::one()).expect("unit polynomial");for(q,_)in f.factor_list().1{p=p.mul(&q);}return Ok(PyPoly {p:p.monic().into()});}
        Ok(PyPoly { p: (self.p.sqf_part()).into() })
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

    fn as_list(&self) -> PyResult<Vec<PyExpr>> {
        self.all_coeffs()
    }

    #[pyo3(signature = (*args))]
    fn subs(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
        let py = args.py();
        let e = self.p.as_expr();
        let pairs: Vec<(Expr, Expr)> = if args.len() == 2 {
            vec![(to_expr(&args.get_item(0)?)?, to_expr(&args.get_item(1)?)?)]
        } else if args.len() == 1 {
            let d = args.get_item(0)?;
            let d = d.downcast::<PyDict>()?;
            d.iter().map(|(k, v)| Ok((to_expr(&k)?, to_expr(&v)?))).collect::<PyResult<_>>()?
        } else {
            return Err(PyTypeError::new_err("subs expects (old, new) or a dict"));
        };
        let r = symox::expand::expand(&symox::subs::subs(&e, &pairs));
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
        let Some(p)=&self.p.algebraic else{return Err(PyTypeError::new_err("polynomial norm requires an algebraic coefficient field"))};
        Ok(PyPoly {p:p.norm().map_err(err)?.into()})
    }

    fn resultant(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        if let Some(field)=self.p.algebraic.as_ref().or(q.algebraic.as_ref()) {let a=field.convert(&self.p.as_expr()).map_err(err)?;let b=field.convert(&q.as_expr()).map_err(err)?;return Ok(w(a.resultant(&b).map_err(err)?.as_expr()));}
        Ok(w(api::resultant(&self.p.as_expr(), &q.as_expr(), self.p.true_gens()).map_err(err)?))
    }

    fn discriminant(&self) -> PyResult<PyExpr> {
        if let Some(p)=&self.p.algebraic{return Ok(w(p.discriminant().map_err(err)?.as_expr()));}
        Ok(w(api::discriminant(&self.p.as_expr(), self.p.true_gens()).map_err(err)?))
    }
}

// ----------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------

fn expression_vector(values:&Bound<'_,PyAny>)->PyResult<Vec<Expr>> {
    values.try_iter()?.map(|v|to_expr(&v?)).collect()
}

#[pyfunction]
#[pyo3(signature=(expr,generators,parameters,domain,modulus))]
fn _rebuild_poly(expr:&Bound<'_,PyAny>,generators:&Bound<'_,PyAny>,parameters:Option<&Bound<'_,PyAny>>,domain:&str,modulus:Option<u64>)->PyResult<PyPoly> {
    let e=to_expr(expr)?;let generators=expression_vector(generators)?;
    if let Some(parameters)=parameters {
        let parameters=expression_vector(parameters)?;
        return Ok(PyPoly {p:FractionPoly::from_expr(&e,&generators,&parameters,domain).map_err(err)?.into()});
    }
    if let Some(field)=NumberField::parse(domain).map_err(err)?{return Ok(PyPoly {p:AlgebraicPoly::from_expr(&e,&generators,field).map_err(err)?.into()});}
    let mut p=Poly::from_expr(&e,&generators).map_err(err)?;
    let base=domain.split(['[','(']).next().unwrap_or(domain);
    p=set_domain(p,Some(base))?;
    if let Some(m)=modulus {p=p.set_modulus(check_mod(m)?);}
    Ok(PyPoly {p:p.into()})
}

#[pyfunction]
#[pyo3(signature=(expressions,generators,parameters,domain,order,polys,modulus))]
fn _rebuild_groebner(py:Python<'_>,expressions:&Bound<'_,PyAny>,generators:&Bound<'_,PyAny>,parameters:Option<&Bound<'_,PyAny>>,domain:&str,order:&str,polys:bool,modulus:Option<u64>)->PyResult<PyGroebnerBasis> {
    let generators=expression_vector(generators)?;
    let args=PyTuple::new(py,generators.iter().cloned().map(w))?;
    let mut result=PyGroebnerBasis::build(expressions,&args,order,if modulus.is_some(){None}else{Some(domain)},modulus,None,Some(polys),None)?;
    if let Some(parameters)=parameters {
        let parameters=expression_vector(parameters)?;
        let prototype=FractionPoly::from_expr(&Expr::zero(),&generators,&parameters,domain).map_err(err)?;
        result.basis=result.basis.iter().map(|p|prototype.convert(&p.as_expr()).map(Into::into)).collect::<Result<_,_>>().map_err(err)?;
        result.fraction_domain=Some(prototype);
    }
    Ok(result)
}

fn algebraic_option_poly(f:&Bound<'_,PyAny>,gens:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Option<AlgebraicPoly>> {
    let Some(kwargs)=kwargs else{return Ok(None)};
    let extension=kwargs.get_item("extension")?;
    let domain=kwargs.get_item("domain")?.map(|d|d.extract::<String>()).transpose()?;
    let(e,g2)=expr_or_coeff_list(f)?;
    let Some(field)=algebraic_extension(extension.as_ref(),std::slice::from_ref(&e),domain.as_deref())?else{return Ok(None)};
    let mut generators=gens_of(gens)?;if generators.is_empty(){generators=g2;}if generators.is_empty(){generators=symox::poly::fraction::infer_parameters(std::slice::from_ref(&e),&[]);}
    AlgebraicPoly::from_expr(&e,&generators,field).map(Some).map_err(err)
}

fn pair<'py>(py: Python<'py>, c: Expr, fs: Vec<(Expr, Expr)>) -> PyResult<Bound<'py, PyTuple>> {
    let l = PyList::empty(py);
    for (f, k) in fs {
        let kk: PyObject = match k.as_i64() {
            Some(i) => i.into_pyobject(py)?.into_any().unbind(),
            None => w(k).into_pyobject(py)?.into_any().unbind(),
        };
        l.append(PyTuple::new(py, [w(f).into_pyobject(py)?.into_any().unbind(), kk])?)?;
    }
    PyTuple::new(py, [w(c).into_pyobject(py)?.into_any(), l.into_any()])
}

#[pyfunction]
#[pyo3(signature = (f, *gens, deep=false, modulus=None, **_kwargs))]
fn factor(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, deep: bool, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if p.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic factorization is not implemented"));}
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
fn sqf(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if p.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic factorization is not implemented"));}
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
fn factor_list<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if p.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic factorization is not implemented"));}
        let(c,fs)=p.factor_list().map_err(err)?;return pair(py,c.as_expr(),fs.into_iter().map(|(p,k)|(p.as_expr(),Expr::int(k as i64))).collect());
    }
    if let Ok(p) = f.downcast::<PyPoly>() {
        return p.get().factor_list(py);
    }
    if let Some(m) = modulus {
        let (c, fs) = api::factor_list_mod(&to_expr(f)?, check_mod(m)?).map_err(err)?;
        return pair(py, c, fs);
    }
    let gens = gens_of(gens)?;
    let e = to_expr(f)?;
    if !gens.is_empty() {
        let (c, fs) = api::poly_factor_list(&e, &gens).map_err(err)?;
        return pair(py, c, fs.into_iter().map(|(f, k)| (f, Expr::int(k as i64))).collect());
    }
    let (c, fs) = api::factor_list(&e).map_err(err)?;
    pair(py, c, fs)
}

#[pyfunction]
#[pyo3(signature = (f, *gens, modulus=None, **_kwargs))]
fn sqf_list<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if p.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate algebraic factorization is not implemented"));}
        let(c,fs)=p.sqf_list();return pair(py,c.as_expr(),fs.into_iter().map(|(p,k)|(p.as_expr(),Expr::int(k as i64))).collect());
    }
    if let Ok(p) = f.downcast::<PyPoly>() {
        return p.get().sqf_list(py);
    }
    if let Some(m) = modulus {
        let (c, fs) = api::sqf_list_mod(&to_expr(f)?, check_mod(m)?).map_err(err)?;
        return pair(py, c, fs);
    }
    let _ = gens;
    let (c, fs) = api::sqf_list(&to_expr(f)?).map_err(err)?;
    pair(py, c, fs)
}

fn list_items(f: &Bound<'_, PyAny>) -> Option<Vec<Expr>> {
    if f.downcast::<PyList>().is_ok() || f.downcast::<PyTuple>().is_ok() {
        let mut v = Vec::new();
        for x in f.try_iter().ok()? {
            v.push(to_expr(&x.ok()?).ok()?);
        }
        return Some(v);
    }
    None
}

fn poly_result(py: Python<'_>, e: Expr, as_poly: bool, gens: &[Expr]) -> PyResult<PyObject> {
    if as_poly {
        let p = Poly::from_expr(&e, gens).map_err(err)?;
        return poly_obj(py, p);
    }
    Ok(w(e).into_pyobject(py)?.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, g=None, *gens, modulus=None, **_kwargs))]
fn gcd(py: Python<'_>, f: &Bound<'_, PyAny>, g: Option<&Bound<'_, PyAny>>, gens: &Bound<'_, PyTuple>, modulus: Option<u64>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if modulus.is_none() {if let (Ok(p),Some(q))=(f.downcast::<PyPoly>(),g){if p.get().p.algebraic.is_some()||p.get().p.field.is_some(){return poly_obj(py,p.get().gcd(q)?.p);}}}
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        let q=g.ok_or_else(||PyTypeError::new_err("expected a second polynomial"))?;
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if p.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate gcd over number fields is not implemented"));}
        let q=p.convert(&to_expr(q)?).map_err(err)?;
        let result=p.gcd(&q);
        return Ok(w(result.as_expr()).into_pyobject(py)?.into_any().unbind());
    }
    if let (Some(m), Some(g)) = (modulus, g) {
        let r = api::gcd_mod(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?, check_mod(m)?).map_err(err)?;
        return Ok(w(r).into_pyobject(py)?.into_any().unbind());
    }
    let as_poly = f.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if let Some(items) = list_items(f) {
        if let Some(g) = g {
            gv.insert(0, to_expr(g)?);
        }
        let mut acc = Expr::zero();
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
    if modulus.is_none() {if let (Ok(p),Some(q))=(f.downcast::<PyPoly>(),g){if p.get().p.algebraic.is_some()||p.get().p.field.is_some(){return poly_obj(py,p.get().lcm(q)?.p);}}}
    if let Some(p)=algebraic_option_poly(f,gens,_kwargs)? {
        let q=g.ok_or_else(||PyTypeError::new_err("expected a second polynomial"))?;
        if modulus.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if p.gens.len()!=1{return Err(PyNotImplementedError::new_err("multivariate gcd over number fields is not implemented"));}
        let q=p.convert(&to_expr(q)?).map_err(err)?;
        let result=if p.rep.is_zero()||q.rep.is_zero(){p.convert(&Expr::zero()).map_err(err)?}else{p.mul(&q).div(&p.gcd(&q)).map_err(err)?.0.monic()};
        return Ok(w(result.as_expr()).into_pyobject(py)?.into_any().unbind());
    }
    if let (Some(m), Some(g)) = (modulus, g) {
        let r = api::lcm_mod(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?, check_mod(m)?).map_err(err)?;
        return Ok(w(r).into_pyobject(py)?.into_any().unbind());
    }
    let as_poly = f.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if let Some(items) = list_items(f) {
        if let Some(g) = g {
            gv.insert(0, to_expr(g)?);
        }
        let mut acc = Expr::one();
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
    let (h, a, b) = api::cofactors(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?;
    PyTuple::new(py, [w(h), w(a), w(b)])
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn cancel(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = gens;
    Ok(w(api::cancel(&to_expr(f)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, x=None, full=false, **_kwargs))]
fn apart(f: &Bound<'_, PyAny>, x: Option<&Bound<'_, PyAny>>, full: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = full;
    let x = match x {
        Some(x) => Some(to_expr(x)?),
        None => None,
    };
    Ok(w(api::apart(&to_expr(f)?, x.as_ref()).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (expr, syms, func=None, evaluate=true, **_kwargs))]
fn collect(expr: &Bound<'_, PyAny>, syms: &Bound<'_, PyAny>, func: Option<&Bound<'_, PyAny>>, evaluate: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let _ = evaluate;
    let py = expr.py();
    let e = to_expr(expr)?;
    let syms = match list_items(syms) {
        Some(v) => v,
        None => vec![to_expr(syms)?],
    };
    let r = api::collect(&e, &syms);
    match func {
        Some(f) => {
            // apply func to each collected coefficient
            let parts: Vec<Expr> = match r.kind() {
                symox::Kind::Add(..) => r.args(),
                _ => vec![r.clone()],
            };
            let mut out = Vec::new();
            for t in parts {
                out.push(to_expr(&f.call1((w(t),))?)?);
            }
            Ok(w(symox::ops::add(&out)).into_pyobject(py)?.into_any().unbind())
        }
        None => Ok(w(r).into_pyobject(py)?.into_any().unbind()),
    }
}

fn kw_modulus(kw: Option<&Bound<'_, PyDict>>) -> PyResult<Option<u64>> {
    match kw {
        Some(d) => match d.get_item("modulus")? {
            Some(v) if !v.is_none() => Ok(Some(check_mod(v.extract::<u64>()?)?)),
            _ => Ok(None),
        },
        None => Ok(None),
    }
}

fn div_like(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, kw: Option<&Bound<'_, PyDict>>) -> PyResult<(Expr, Expr, Vec<Expr>, bool)> {
    let as_poly = f.downcast::<PyPoly>().is_ok() || g.downcast::<PyPoly>().is_ok();
    if let Some(m) = kw_modulus(kw)? {
        let (q, r) = api::div_mod(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?, m).map_err(err)?;
        return Ok((q, r, vec![], false));
    }
    let mut gv = gens_of(gens)?;
    if gv.is_empty() {
        if let Ok(p) = f.downcast::<PyPoly>() {
            gv = p.get().p.true_gens().to_vec();
        }
    }
    let (q, r) = api::div_expr(&to_expr(f)?, &to_expr(g)?, &gv).map_err(err)?;
    Ok((q, r, gv, as_poly))
}

#[pyfunction]
#[pyo3(name = "div", signature = (f, g, *gens, **_kwargs))]
fn py_div<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if let (Ok(a), Ok(_)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
        return a.get().div(py, g);
    }
    let (q, r, gv, as_poly) = div_like(f, g, gens, _kwargs)?;
    PyTuple::new(py, [poly_result(py, q, as_poly, &gv)?, poly_result(py, r, as_poly, &gv)?])
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, order="lex", domain=None, modulus=None, polys=None, auto=None, extension=None, **_kwargs))]
fn reduced<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, order: &str, domain: Option<&str>, modulus: Option<u64>, polys: Option<bool>, auto: Option<bool>, extension: Option<&Bound<'_,PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    use symox::poly::convert::{from_exprs, from_exprs_gens, to_expr_f};
    use symox::poly::groebner::{reduce_with_options, CoefficientField, MonomialOrder};
    let monomial_order = MonomialOrder::parse(order)
        .ok_or_else(|| PyValueError::new_err("order must be 'lex', 'grlex', or 'grevlex'"))?;
    // SymPy disables automatic field extension for explicit domains unless
    // auto=True is supplied explicitly.
    let auto = auto.unwrap_or(domain.is_none());
    if let Some(kwargs) = _kwargs {
        if let Some((key,_)) = kwargs.iter().next() {
            return Err(PyTypeError::new_err(format!("reduced() got an unsupported keyword argument '{}'",key.extract::<String>()?)));
        }
    }
    if domain.is_some_and(|d| d != "ZZ" && d != "QQ" && d != "RR" && !d.starts_with("QQ(") && !d.starts_with("ZZ(") && !d.starts_with("QQ[") && !d.starts_with("ZZ[") && !d.starts_with("QQ<")) {
        return Err(PyNotImplementedError::new_err("reduced supports ZZ, QQ, RR, and prime modulus domains"));
    }
    if domain.is_some() && modulus.is_some() {
        return Err(PyValueError::new_err("specify either domain or modulus"));
    }
    let mut expressions = Vec::new();
    let mut inferred_gens = Vec::new();
    let mut inferred_domain = Domain::ZZ;
    let mut inferred_algebraic:Option<Arc<NumberField>>=None;
    let mut inferred_prime = None;
    let mut any_poly = false;
    let mut append = |item: &Bound<'_,PyAny>| -> PyResult<()> {
        if let Ok(p) = item.downcast::<PyPoly>() {
            any_poly = true;
            let p = &p.get().p;
            expressions.push(p.as_expr());
            if let Some(a)=&p.algebraic{inferred_algebraic=Some(if let Some(previous)=inferred_algebraic.take(){NumberField::from_generators(&[previous.generator.clone(),a.number_field.generator.clone()]).map_err(err)?}else{a.number_field.clone()});}
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
    let mut gv = gens_of(gens)?;
    if gv.is_empty() && any_poly {gv = inferred_gens;}
    let algebraic=algebraic_extension(extension,&expressions,domain)?.or(if domain.is_none(){inferred_algebraic}else{None});
    if let Some(number_field)=algebraic {
        if modulus.is_some() || inferred_prime.is_some(){return Err(PyValueError::new_err("algebraic extensions cannot use a finite modulus"));}
        if gv.is_empty(){gv=symox::poly::fraction::infer_parameters(&expressions,&[]);}
        let prototype=AlgebraicPoly::from_expr(&Expr::zero(),&gv,number_field.clone()).map_err(err)?;
        let inputs=expressions.iter().map(|e|prototype.convert(e)).collect::<Result<Vec<_>,_>>().map_err(err)?;
        let divisors:Vec<_>=inputs.iter().skip(1).map(|p|p.rep.clone()).collect();
        if divisors.iter().any(|p|p.is_zero()){return Err(PyZeroDivisionError::new_err("polynomial division by zero"));}
        let _scope=AlgebraicScope::new(number_field);
        let(mut q,r)=symox::poly::fieldalg::reduce(&inputs[0].rep,&divisors,monomial_order);
        if inputs[0].rep.is_zero(){q.clear();}
        let object=|rep|->PyResult<PyObject>{let p=prototype.with_rep(rep);if polys.unwrap_or(any_poly){poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
        let qs=PyList::new(py,q.into_iter().map(object).collect::<PyResult<Vec<_>>>()?)?;
        return PyTuple::new(py,[qs.into_any(),object(r)?.into_bound(py)]);
    }
    let available=symox::poly::fraction::infer_parameters(&expressions,&[]);
    let declared=domain.map(|d|symox::poly::fraction::parse_domain(d,&available)).transpose().map_err(err)?.flatten();
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
        let prototype=FractionPoly::from_expr(&Expr::zero(),&gv,&parameters,&display).map_err(err)?;
        let inputs=expressions.iter().map(|e|prototype.convert(e)).collect::<Result<Vec<_>,_>>().map_err(err)?;
        if display.ends_with(']') && inputs.iter().any(|p|!p.belongs_to_ring()){return Err(PyValueError::new_err("coefficient does not belong to the selected polynomial ring"));}
        let divisors:Vec<_>=inputs.iter().skip(1).map(|p|p.rep.clone()).collect();
        if divisors.iter().any(|p|p.is_zero()){return Err(PyZeroDivisionError::new_err("polynomial division by zero"));}
        let _scope=symox::poly::ratfunc::Scope::new(parameters.len());
        let(mut q,r)=symox::poly::fieldalg::reduce_filtered(&inputs[0].rep,&divisors,monomial_order,|c|auto || !display.ends_with(']') || symox::poly::fraction::coefficient_in_ring(c,display.starts_with("ZZ")));
        if inputs[0].rep.is_zero(){q.clear();}
        let promote=q.iter().chain(std::iter::once(&r)).any(|p|p.coeffs().iter().any(|c|!symox::poly::fraction::coefficient_in_ring(c,display.starts_with("ZZ"))));
        let object=|rep|->PyResult<PyObject>{let mut p=prototype.with_rep(rep);if promote {p.display=p.display.replace('[',"(").replace(']',")");}if return_polys {poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}};
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
fn rem(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if let (Ok(a), Ok(_)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
        return Ok(a.get().rem(g)?.into_pyobject(py)?.into_any().unbind());
    }
    let (_, r, gv, as_poly) = div_like(f, g, gens, _kwargs)?;
    poly_result(py, r, as_poly, &gv)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn quo(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if let (Ok(a), Ok(_)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
        return Ok(a.get().quo(g)?.into_pyobject(py)?.into_any().unbind());
    }
    let (q, _, gv, as_poly) = div_like(f, g, gens, _kwargs)?;
    poly_result(py, q, as_poly, &gv)
}

fn pdiv_like(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>) -> PyResult<(Expr, Expr, Vec<Expr>, bool)> {
    let as_poly = f.downcast::<PyPoly>().is_ok() || g.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if gv.is_empty() {
        if let Ok(p) = f.downcast::<PyPoly>() {
            gv = p.get().p.true_gens().to_vec();
        }
    }
    let (q, r) = api::pdiv_expr(&to_expr(f)?, &to_expr(g)?, &gv).map_err(err)?;
    Ok((q, r, gv, as_poly))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn pdiv<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    let (q, r, gv, as_poly) = pdiv_like(f, g, gens)?;
    PyTuple::new(py, [poly_result(py, q, as_poly, &gv)?, poly_result(py, r, as_poly, &gv)?])
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn prem(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (_, r, gv, as_poly) = pdiv_like(f, g, gens)?;
    poly_result(py, r, as_poly, &gv)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn pquo(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (q, _, gv, as_poly) = pdiv_like(f, g, gens)?;
    poly_result(py, q, as_poly, &gv)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn pexquo(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::pexquo_expr(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn exquo(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::exquo(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn resultant(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::resultant(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn discriminant(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::discriminant(&to_expr(f)?, &gens_of(gens)?).map_err(err)?))
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
    let gv = gens_of(gens)?;
    let gv = if gv.is_empty() { g2 } else { gv };
    let mut rs = proots::roots(&e, &gv).map_err(err)?;
    if let Some(flt) = filter {
        rs.retain(|(r, _)| match flt {
            "Z" => r.is_integer_number(),
            "Q" => r.is_rational_number(),
            "R" => r.mask() == 0 && symox::evalf::eval_c64(r).is_some_and(|z| z.im == 0.0) && !symox::traverse::has(r, &Expr::i()),
            "I" => r.mask() == 0 && symox::evalf::eval_c64(r).is_some_and(|z| z.re.abs() < 1e-300 && z.im != 0.0),
            _ => true,
        });
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
#[pyo3(signature = (f, n=15, **_kwargs))]
fn nroots(f: &Bound<'_, PyAny>, n: u32, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
    let (e, g2) = expr_or_coeff_list(f)?;
    Ok(proots::nroots(&e, &g2, n).map_err(err)?.into_iter().map(w).collect())
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

#[pyfunction]
#[pyo3(name = "_poly_together", signature = (f, deep=false, **_kwargs))]
fn py_together(f: &Bound<'_, PyAny>, deep: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(w(api::together(&to_expr(f)?)))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn gcdex<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    let (s, t, h) = api::gcdex(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?;
    PyTuple::new(py, [w(s), w(t), w(h)])
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn half_gcdex<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    let (s, h) = api::half_gcdex(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?;
    PyTuple::new(py, [w(s), w(h)])
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn invert(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::invert(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn terms_gcd(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = gens;
    Ok(w(api::terms_gcd(&to_expr(f)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn compose(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::compose(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn decompose(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
    Ok(api::decompose(&to_expr(f)?, &gens_of(gens)?).map_err(err)?.into_iter().map(w).collect())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn sqf_part(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::sqf_part(&to_expr(f)?, &gens_of(gens)?).map_err(err)?))
}

fn minimal_polynomial_value(model:&symox::poly::algebraic_model::PrimitiveModel,generator:Expr)->Poly {
    let(_,coefficients)=symox::poly::dense::clear_denoms(&model.minimal);
    Poly {rep:symox::poly::dense::to_mpoly(&coefficients,0,1).to_rat(),gens:vec![generator],ngens:1,domain:Domain::QQ,modulus:None}
}
#[pyfunction]
#[pyo3(signature=(ex,x=None,compose=true,polys=false,domain=None))]
fn minimal_polynomial(py:Python<'_>,ex:&Bound<'_,PyAny>,x:Option<&Bound<'_,PyAny>>,compose:bool,polys:bool,domain:Option<&str>)->PyResult<PyObject> {
    if !compose{return Err(PyNotImplementedError::new_err("compose=False minimal polynomials are not implemented"));}
    if domain.is_some_and(|d|d!="QQ"){return Err(PyNotImplementedError::new_err("minimal polynomials currently support ground domain QQ"));}
    let generator=x.map(to_expr).transpose()?.unwrap_or_else(||Expr::symbol("_x"));
    if !generator.is_symbol(){return Err(PyValueError::new_err("minimal polynomial generator must be a symbol"));}
    let model=symox::poly::algebraic_model::PrimitiveModel::build(&to_expr(ex)?).map_err(err)?;
    let p=minimal_polynomial_value(&model,generator);
    if polys{poly_obj(py,p)}else{Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())}
}
#[pyfunction]
#[pyo3(signature=(ex,x=None,compose=true,polys=false,domain=None))]
fn minpoly(py:Python<'_>,ex:&Bound<'_,PyAny>,x:Option<&Bound<'_,PyAny>>,compose:bool,polys:bool,domain:Option<&str>)->PyResult<PyObject> {minimal_polynomial(py,ex,x,compose,polys,domain)}
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
    m.add_function(wrap_pyfunction!(degree, m)?)?;
    m.add_function(wrap_pyfunction!(roots, m)?)?;
    m.add_function(wrap_pyfunction!(nroots, m)?)?;
    m.add_function(wrap_pyfunction!(real_roots, m)?)?;
    m.add_function(wrap_pyfunction!(gcdex, m)?)?;
    m.add_function(wrap_pyfunction!(half_gcdex, m)?)?;
    m.add_function(wrap_pyfunction!(invert, m)?)?;
    m.add_function(wrap_pyfunction!(terms_gcd, m)?)?;
    m.add_function(wrap_pyfunction!(compose, m)?)?;
    m.add_function(wrap_pyfunction!(decompose, m)?)?;
    m.add_function(wrap_pyfunction!(sqf_part, m)?)?;
    m.add_function(wrap_pyfunction!(py_together, m)?)?;
    Ok(())
}
