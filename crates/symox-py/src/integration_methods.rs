//! The Risch, manual and heurisch methods of `integrate`, provided to the
//! native integrator by the Python ports in `symox.integrals`.

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use symox::expr::{Expr, Func, Kind};
use symox::integrate::host_methods::{self, HostMethod, HostOutcome};
use symox::integrate::Conds;

use crate::{to_expr, wrap, PyExprNode};

static LOG_ATAN_POLY_METHODS: std::sync::OnceLock<Vec<(&'static str, PyObject)>> = std::sync::OnceLock::new();
static LOG_ATAN_DOMAIN_FLAG: std::sync::OnceLock<PyObject> = std::sync::OnceLock::new();
static LOG_REAL_POLY_METHODS: std::sync::OnceLock<Vec<(&'static str, PyObject)>> = std::sync::OnceLock::new();

static LRT_POLY_METHODS: std::sync::OnceLock<Vec<(&'static str, PyObject)>> = std::sync::OnceLock::new();

static RATIONAL_POLY_METHODS: std::sync::OnceLock<Vec<(&'static str, PyObject)>> = std::sync::OnceLock::new();
static RATINT_POLY_METHODS: std::sync::OnceLock<Vec<(&'static str, PyObject)>> = std::sync::OnceLock::new();

pub(crate) fn register_rational_reduction(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let class = m.py().get_type::<crate::poly::PyPoly>();
    let methods = ["cofactors", "diff", "degree", "quo"]
        .into_iter().map(|name| Ok((name, class.getattr(name)?.unbind())))
        .collect::<PyResult<Vec<_>>>()?;
    let _ = RATIONAL_POLY_METHODS.set(methods);
    let methods = ["degree", "to_field", "div", "gcdex", "quo", "__neg__", "__mul__", "__add__", "domain", "gens", "is_zero", "as_expr"]
        .into_iter().map(|name| Ok((name, class.getattr(name)?.unbind())))
        .collect::<PyResult<Vec<_>>>()?;
    let _ = LOG_ATAN_POLY_METHODS.set(methods);
    let _ = LOG_ATAN_DOMAIN_FLAG.set(m.py().get_type::<crate::poly::domains::PyDomain>().getattr("is_EX")?.unbind());
    let methods = ["__new__", "_from_normalized", "to_field", "set_domain", "diff", "degree", "LC", "terms", "sqf_list", "primitive", "quo", "invert", "coeffs", "monoms", "gcd", "__bool__", "__sub__", "__mul__", "__pow__", "gens", "domain", "as_expr"]
        .into_iter().map(|name| Ok((name, class.getattr(name)?.unbind())))
        .collect::<PyResult<Vec<_>>>()?;
    let _ = LRT_POLY_METHODS.set(methods);
    let methods = ["__new__", "__init__", "as_expr", "count_roots", "__bool__", "domain", "gens", "degree", "to_field", "div", "gcdex", "quo", "__neg__", "__mul__", "__add__", "__pow__"]
        .into_iter().map(|name| Ok((name, class.getattr(name)?.unbind())))
        .collect::<PyResult<Vec<_>>>()?;
    let _ = LOG_REAL_POLY_METHODS.set(methods);
    let methods = ["cancel", "integrate", "is_zero"]
        .into_iter().map(|name| Ok((name, class.getattr(name)?.unbind())))
        .collect::<PyResult<Vec<_>>>()?;
    let _ = RATINT_POLY_METHODS.set(methods);
    m.add_function(wrap_pyfunction!(_ratint_driver, m)?)?;
    m.add_function(wrap_pyfunction!(_log_to_real, m)?)?;
    m.add_function(wrap_pyfunction!(_ratint_logpart, m)?)?;
    m.add_function(wrap_pyfunction!(_ratint_ratpart, m)?)?;
    m.add_function(wrap_pyfunction!(_log_to_atan, m)?)
}

/// Full QQ orchestration; all other inputs are refused without consuming a
/// user iterator, sympification hook, tuple subclass or mathematical protocol.
#[pyfunction]
pub(crate) fn _ratint_driver(f:&Bound<'_,PyAny>,x:&Bound<'_,PyAny>,flags:&Bound<'_,pyo3::types::PyDict>)
    -> PyResult<Option<(crate::PyExpr,bool)>>
{
    use pyo3::types::{PyBool,PyInt,PyString};
    use crate::poly::PyPoly;
    use symox::integrate::ratint_driver::Input;
    let py=f.py();py.check_signals()?;
    let scalar=|value:&Bound<'_,PyAny>,integer:bool|->PyResult<Option<Expr>>{
        if integer && value.is_exact_instance_of::<PyInt>(){return Ok(Some(to_expr(value)?));}
        if !value.is_instance_of::<PyExprNode>(){return Ok(None);}
        let e=to_expr(value)?;
        Ok(crate::node_type(py,&e).filter(|c|value.get_type().is(c)).map(|_|e))
    };
    let Some(x)=scalar(x,false)? else{return Ok(None)};
    let Some(variable)=x.as_sym() else{return Ok(None)};
    if variable.is_wild() || variable.host().is_some() || variable.scalar_identity().is_some()
        || symox::decompose::commutativity(&x)!=Some(true){return Ok(None);}
    let input=if f.is_exact_instance_of::<PyTuple>() {
        let pair=f.downcast::<PyTuple>()?;
        if pair.len()!=2{return Ok(None);}
        let (Some(p),Some(q))=(scalar(&pair.get_item(0)?,true)?,scalar(&pair.get_item(1)?,true)?) else{return Ok(None)};
        Input::Fraction(p,q)
    }else{
        let Some(f)=scalar(f,false)? else{return Ok(None)};
        Input::Expression(f)
    };
    // Ordinary kwargs keys cannot invoke equality/conversion callbacks.
    if flags.len()>64 || flags.iter().any(|(key,_)|!key.is_exact_instance_of::<PyString>()){return Ok(None);}
    if let Some(real)=flags.get_item("real")? {
        if !real.is_none() && (!real.is_exact_instance_of::<PyBool>() || !real.is_truthy()?){return Ok(None);}
    }
    let symbol=match flags.get_item("symbol")? {
        None=>"t".to_owned(),
        Some(value) if value.is_exact_instance_of::<PyString>()=>{
            let Ok(name)=value.extract::<String>() else{return Ok(None)};name
        },
        Some(value)=>{
            let Some(e)=scalar(&value,false)? else{return Ok(None)};
            let Some(s)=e.as_sym().filter(|s|!s.is_wild() && s.host().is_none() && s.scalar_identity().is_none()) else{return Ok(None)};
            s.name().to_owned()
        },
    };
    let class=py.get_type::<PyPoly>();let own=class.getattr("__dict__")?;
    if own.contains("__getattribute__")? || own.contains("__getattr__")? {return Ok(None);}
    let mut names=std::collections::HashSet::new();
    for registry in [RATINT_POLY_METHODS.get(),RATIONAL_POLY_METHODS.get(),LRT_POLY_METHODS.get(),LOG_REAL_POLY_METHODS.get(),LOG_ATAN_POLY_METHODS.get()] {
        let Some(methods)=registry else{return Ok(None)};
        for (name,original) in methods {
            if !names.insert(*name){continue;}
            match *name {
                "as_expr"=>{
                    let module=py.import("symox.polys._evaluation")?;
                    if !class.getattr(*name)?.is(&module.getattr(*name)?)
                        || !module.getattr("_native_as_expr")?.is(original.bind(py)){return Ok(None);}
                },
                "LC"|"terms"|"coeffs"|"monoms"=>{
                    let module=py.import("symox.polys._poly_methods")?;
                    if !class.getattr(*name)?.is(&module.getattr(*name)?)
                        || !module.getattr("_native")?.get_item(*name)?.is(original.bind(py)){return Ok(None);}
                },
                _=>if !class.getattr(*name)?.is(original.bind(py)){return Ok(None);},
            }
        }
    }
    let Some(flag)=LOG_ATAN_DOMAIN_FLAG.get() else{return Ok(None)};
    if !py.get_type::<crate::poly::domains::PyDomain>().getattr("is_EX")?.is(flag.bind(py)){return Ok(None);}
    let value=symox::integrate::ratint_driver::integrate(&input,&x,&symbol,&mut ||py.check_signals())?;
    py.check_signals()?;Ok(value.map(|(value,logarithmic)|(wrap(value),logarithmic)))
}

/// Select stored rational polynomials only. No sympification, generator
/// protocol or general coefficient-domain conversion is attempted here.
#[pyfunction]
pub(crate) fn _log_to_real(h: &Bound<'_,PyAny>,q: &Bound<'_,PyAny>,x: &Bound<'_,PyAny>,t: &Bound<'_,PyAny>)
    -> PyResult<Option<crate::PyExpr>>
{
    use crate::poly::PyPoly;
    use symox::poly::{dense,polyobj::Domain};
    let py=h.py();py.check_signals()?;
    if !h.is_exact_instance_of::<PyPoly>() || !q.is_exact_instance_of::<PyPoly>() {return Ok(None);}
    let class=h.get_type();let own=class.getattr("__dict__")?;
    if own.contains("__getattribute__")? || own.contains("__getattr__")? {return Ok(None);}
    let Some(methods)=LOG_REAL_POLY_METHODS.get() else{return Ok(None)};
    for (name,original) in methods {
        if *name=="as_expr" {
            let module=py.import("symox.polys._evaluation")?;
            if !class.getattr(*name)?.is(&module.getattr("as_expr")?)
                || !module.getattr("_native_as_expr")?.is(original.bind(py)) {return Ok(None);}
        } else if !class.getattr(*name)?.is(original.bind(py)) {return Ok(None);}
    }
    let Some(flag)=LOG_ATAN_DOMAIN_FLAG.get() else{return Ok(None)};
    if !py.get_type::<crate::poly::domains::PyDomain>().getattr("is_EX")?.is(flag.bind(py)){return Ok(None);}
    let canonical=|value:&Bound<'_,PyAny>|->PyResult<Option<Expr>>{
        if !value.is_instance_of::<PyExprNode>(){return Ok(None);}
        let e=to_expr(value)?;
        if !matches!(e.kind(),Kind::Sym(_)) || !crate::node_type(py,&e).is_some_and(|c|value.get_type().is(&c))
            || symox::decompose::commutativity(&e)!=Some(true) {return Ok(None);}
        Ok(Some(e))
    };
    let (Some(x),Some(t))=(canonical(x)?,canonical(t)?) else{return Ok(None)};
    if x==t{return Ok(None);}
    let (h,q)=(h.downcast::<PyPoly>()?.get(),q.downcast::<PyPoly>()?.get());
    let q=&q.p;
    if q.field.is_some() || q.algebraic.is_some() || q.expression.is_some()
        || q.base.modulus.is_some() || !matches!(q.base.domain,Domain::ZZ|Domain::QQ)
        || q.base.gens!=[t.clone()] || q.base.ngens!=1
        || !matches!(q.base.rep.degree(0),1|2) {return Ok(None);}
    let Some(polynomial)=symox::integrate::log_to_real::rational_input(&h.p,&x,&t,&mut ||py.check_signals())?
        else{return Ok(None)};
    let q=dense::from_mpoly(&q.base.rep,0);
    let result=symox::integrate::log_to_real::convert(&polynomial,&q,&x,&t,&mut ||py.check_signals())?;
    py.check_signals()?;Ok(result.map(wrap))
}

/// The host constructs its input Polys and its optional Dummy before selection.
/// No user protocol is called while testing the stored coefficient domains.
#[pyfunction]
pub(crate) fn _ratint_logpart(f: &Bound<'_,PyAny>,g: &Bound<'_,PyAny>,x: &Bound<'_,PyAny>,t: &Bound<'_,PyAny>)
    -> PyResult<Option<Vec<(crate::poly::PyPoly,crate::poly::PyPoly)>>>
{
    use crate::poly::PyPoly;
    use symox::poly::polyobj::Domain;
    let py=f.py();py.check_signals()?;
    if !f.is_exact_instance_of::<PyPoly>() || !g.is_exact_instance_of::<PyPoly>() {return Ok(None);}
    let class=f.get_type();
    let own=class.getattr("__dict__")?;
    if own.contains("__getattribute__")? || own.contains("__getattr__")? {return Ok(None);}
    let Some(methods)=LRT_POLY_METHODS.get() else{return Ok(None)};
    for (name,original) in methods {
        match *name {
            "as_expr" => {
                let module=py.import("symox.polys._evaluation")?;
                if !class.getattr(*name)?.is(&module.getattr("as_expr")?)
                    || !module.getattr("_native_as_expr")?.is(original.bind(py)) {return Ok(None);}
            },
            "LC" | "coeffs" | "monoms" | "terms" => {
                let module=py.import("symox.polys._poly_methods")?;
                if !class.getattr(*name)?.is(&module.getattr(*name)?)
                    || !module.getattr("_native")?.get_item(*name)?.is(original.bind(py)) {return Ok(None);}
            },
            _ => if !class.getattr(*name)?.is(original.bind(py)) {return Ok(None);},
        }
    }
    let (f,g)=(f.downcast::<PyPoly>()?.get(),g.downcast::<PyPoly>()?.get());
    let plain=|p:&symox::poly::fraction::Polynomial|p.field.is_none() && p.algebraic.is_none()
        && p.expression.is_none() && !p.base.has_domain_gens() && p.base.modulus.is_none()
        && matches!(p.base.domain,Domain::ZZ|Domain::QQ) && p.base.gens.len()==1;
    if !plain(&f.p) || !plain(&g.p) || f.p.base.gens!=g.p.base.gens {return Ok(None);}
    let canonical=|value:&Bound<'_,PyAny>|->PyResult<Option<Expr>>{
        if !value.is_instance_of::<PyExprNode>(){return Ok(None);}
        let e=to_expr(value)?;
        if !matches!(e.kind(),Kind::Sym(_)) || !crate::node_type(py,&e).is_some_and(|c|value.get_type().is(&c)) {return Ok(None);}
        Ok(Some(e))
    };
    let (Some(x),Some(t))=(canonical(x)?,canonical(t)?) else{return Ok(None)};
    if x!=f.p.base.gens[0] || x==t
        || symox::decompose::commutativity(&x)!=Some(true)
        || symox::decompose::commutativity(&t)!=Some(true) {return Ok(None);}
    let Some(pairs)=symox::integrate::logarithmic_reduction::reduce(&f.p.base,&g.p.base,&t,&mut ||py.check_signals())?
        else{return Ok(None)};
    py.check_signals()?;
    Ok(Some(pairs.into_iter().map(|(h,q)|(PyPoly{p:h},PyPoly{p:q})).collect()))
}

/// Exact stock-polynomial selection only; custom domains keep rationaltools'
/// complete host algorithm. Degree bounds select a route, not a capability.
#[pyfunction]
pub(crate) fn _ratint_ratpart(
    numerator: &Bound<'_, PyAny>, denominator: &Bound<'_, PyAny>,
    variable: &Bound<'_, PyAny>,
) -> PyResult<Option<(crate::PyExpr, crate::PyExpr)>> {
    use crate::poly::PyPoly;
    use symox::poly::{dense, convert, polyobj::Domain};
    if !numerator.is_exact_instance_of::<PyPoly>()
        || !denominator.is_exact_instance_of::<PyPoly>() {
        return Ok(None);
    }
    let class = numerator.get_type();
    let Some(methods) = RATIONAL_POLY_METHODS.get() else { return Ok(None) };
    for (name, original) in methods {
        if !class.getattr(*name)?.is(original.bind(numerator.py())) { return Ok(None); }
    }
    // These two methods acquire stock presentation adapters during bootstrap,
    // after the native module has registered its descriptors.
    for (name, module) in [("as_expr", "symox.polys._evaluation"), ("LC", "symox.polys._poly_methods")] {
        if !class.getattr(name)?.is(&numerator.py().import(module)?.getattr(name)?) {
            return Ok(None);
        }
    }
    let f = numerator.downcast::<PyPoly>()?.get();
    let g = denominator.downcast::<PyPoly>()?.get();
    let plain = |p: &symox::poly::fraction::Polynomial| {
        p.field.is_none() && p.algebraic.is_none() && p.expression.is_none()
            && !p.base.has_domain_gens() && p.base.modulus.is_none()
            && matches!(p.base.domain, Domain::ZZ | Domain::QQ)
            && p.base.gens.len() == 1
    };
    if !plain(&f.p) || !plain(&g.p) || f.p.base.gens != g.p.base.gens
        || g.p.base.rep.degree(0) > 128 || f.p.base.rep.degree(0) > 128 {
        return Ok(None);
    }
    let x = &f.p.base.gens[0];
    // No conversion of a user-defined generator or scalar callback.
    if !matches!(x.kind(), Kind::Sym(_))
        || !variable.is_instance_of::<PyExprNode>()
        || !crate::node_type(variable.py(), x).is_some_and(|t| variable.get_type().is(&t)) {
        return Ok(None);
    }
    if to_expr(variable)? != *x { return Ok(None); }
    let py = numerator.py();
    py.check_signals()?;
    let f = dense::from_mpoly(&f.p.base.rep, 0);
    let g = dense::from_mpoly(&g.p.base.rep, 0);
    let Some(result) = symox::integrate::rational_reduction::reduce(
        &f, &g, &mut || py.check_signals(),
    )? else { return Ok(None) };
    let fraction = |a: &[symox::poly::Rat], b: &[symox::poly::Rat]| {
        let expression = |p| convert::to_expr(&dense::to_mpoly(p, 0, 1), std::slice::from_ref(x));
        symox::poly::api::cancel(&symox::ops::div(&expression(a), &expression(b)))
            .map(crate::wrap).map_err(crate::poly::err)
    };
    let rational = fraction(&result.rational_numerator, &result.rational_denominator)?;
    let logarithmic = fraction(&result.logarithmic_numerator, &result.logarithmic_denominator)?;
    py.check_signals()?;
    Ok(Some((rational, logarithmic)))
}

/// Select only exact native univariate QQ/ZZ polynomials. The host adapter
/// separately guards its atan constructor, evaluation context and recursion.
#[pyfunction]
pub(crate) fn _log_to_atan(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>)
    -> PyResult<Option<crate::PyExpr>>
{
    use crate::poly::PyPoly;
    use symox::poly::{dense, convert, polyobj::Domain};
    let py=f.py();
    py.check_signals()?;
    if !f.is_exact_instance_of::<PyPoly>() || !g.is_exact_instance_of::<PyPoly>() { return Ok(None); }
    let class=f.get_type();
    if class.hasattr("__getattribute__")? {
        let own=class.getattr("__dict__")?;
        if own.contains("__getattribute__")? || own.contains("__getattr__")? { return Ok(None); }
    }
    let Some(methods)=LOG_ATAN_POLY_METHODS.get() else { return Ok(None) };
    for (name,original) in methods {
        if *name=="as_expr" {
            let evaluation=py.import("symox.polys._evaluation")?;
            if !class.getattr(*name)?.is(&evaluation.getattr("as_expr")?)
                || !evaluation.getattr("_native_as_expr")?.is(original.bind(py)) { return Ok(None); }
        } else if !class.getattr(*name)?.is(original.bind(py)) { return Ok(None); }
    }
    let Some(original)=LOG_ATAN_DOMAIN_FLAG.get() else { return Ok(None) };
    if !py.get_type::<crate::poly::domains::PyDomain>().getattr("is_EX")?.is(original.bind(py)) { return Ok(None); }
    let (f,g)=(f.downcast::<PyPoly>()?.get(),g.downcast::<PyPoly>()?.get());
    let plain=|p: &symox::poly::fraction::Polynomial| p.field.is_none() && p.algebraic.is_none()
        && p.expression.is_none() && !p.base.has_domain_gens() && p.base.modulus.is_none()
        && matches!(p.base.domain,Domain::ZZ|Domain::QQ) && p.base.gens.len()==1;
    if !plain(&f.p) || !plain(&g.p) || f.p.base.gens!=g.p.base.gens { return Ok(None); }
    let x=&f.p.base.gens[0];
    if !matches!(x.kind(),Kind::Sym(_)) || f.p.base.rep.degree(0)>128 || g.p.base.rep.degree(0)>128 { return Ok(None); }
    let (f,g)=(dense::from_mpoly(&f.p.base.rep,0),dense::from_mpoly(&g.p.base.rep,0));
    let Some(arguments)=symox::integrate::log_to_atan::arguments(&f,&g,&mut ||py.check_signals())? else { return Ok(None) };
    // Match A + recurse(s,t), including native automatic evaluation of each
    // atan argument. No addition formula changes its principal branch.
    let mut terms=Vec::with_capacity(arguments.len());
    for p in arguments {
        py.check_signals()?;
        let arg=convert::to_expr(&dense::to_mpoly(&p,0,1),std::slice::from_ref(x));
        terms.push(symox::ops::mul2(&Expr::int(2),&symox::functions::atan(&arg)));
    }
    let result=symox::ops::add(&terms);
    py.check_signals()?;
    Ok(Some(wrap(result)))
}

/// Expose exact Cartesian numeric polynomial coefficients to heurisch.
#[pyfunction]
pub(crate) fn _heurisch_cartesian_numeric_coefficients(
    polynomial: &Bound<'_, PyAny>, generators: Vec<Bound<'_, PyAny>>,
) -> PyResult<PyObject> {
    let generators = generators.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    let duplicate = generators.iter().enumerate().any(|(index, generator)| generators[..index].contains(generator));
    let noncommutative = generators.iter().any(|generator| symox::decompose::commutativity(generator) == Some(false));
    if duplicate || noncommutative {
        let kind = if duplicate { "duplicated" } else { "non-commutative" };
        let message = format!("{kind} generators: {}", symox::special::tuple(generators));
        return Err(crate::poly::reference_error(polynomial.py(), "GeneratorsError", &message));
    }
    let (value, processed) = symox::integrate::heurisch_cartesian_numeric_coefficients_with_status(
        &to_expr(polynomial)?, &generators,
    ).map_err(|error| pyo3::exceptions::PyNotImplementedError::new_err(error.to_string()))?;
    if processed {
        Ok(wrap(value).into_pyobject(polynomial.py())?.into_any().unbind())
    } else {
        Ok(polynomial.clone().unbind())
    }
}

fn conds_name(conds: Conds) -> &'static str {
    match conds {
        Conds::Piecewise => "piecewise",
        Conds::Separate => "separate",
        Conds::None => "none",
    }
}

/// The host callback of the native integrator. The ports catch what the
/// reference catches (NotImplementedError from `risch_integrate`,
/// PolynomialError from `heurisch`, ValueError and PolynomialError from
/// `manualintegrate`). When a `True` option selects the method, the port is
/// the whole reference computation and any other exception is raised when
/// the computation returns to Python, as the reference raises it out of
/// `integrate` (`integrate(a**(t - s)/(a**t + 1), t, heurisch=True)` raises
/// the TypeError of `factor_list`). In the default chain, where the native
/// rules precede the ported stages, an ordinary exception only leaves the
/// stage to the next method; an interrupt (or another BaseException) is
/// always raised.
fn integration_host(method: HostMethod, f: &Expr, x: &Expr, conds: Conds) -> HostOutcome {
    if symox::host::failure_pending() {
        return HostOutcome::Unavailable;
    }
    Python::with_gil(|py| {
        let run = || -> PyResult<HostOutcome> {
            let name = match method {
                HostMethod::Risch => "risch",
                HostMethod::Manual => "manual",
                HostMethod::Term => "term",
                HostMethod::ManualOnly => "manual_only",
                HostMethod::HeurischOnly => "heurisch_only",
            };
            let module = py.import("symox.integrals._methods")?;
            let result = module.getattr("_native_host")?
                .call1((name, wrap(f.clone()), wrap(x.clone()), conds_name(conds)))?;
            if result.is_none() {
                return Ok(HostOutcome::Unavailable);
            }
            if let Ok(pair) = result.downcast::<PyTuple>() {
                return Ok(HostOutcome::NonElementary(to_expr(&pair.get_item(0)?)?, to_expr(&pair.get_item(1)?)?));
            }
            Ok(HostOutcome::Value(to_expr(&result)?))
        };
        match run() {
            Ok(outcome) => outcome,
            Err(error) => {
                let selected = match method {
                    HostMethod::ManualOnly | HostMethod::HeurischOnly => true,
                    HostMethod::Risch => host_methods::methods().risch == Some(true),
                    HostMethod::Manual | HostMethod::Term => false,
                };
                if selected || !error.is_instance_of::<pyo3::exceptions::PyException>(py) {
                    crate::host::report(py, error, None);
                }
                HostOutcome::Unavailable
            }
        }
    })
}

pub(crate) fn register() {
    host_methods::register_integration_host(integration_host);
}

/// The `(risch, manual, heurisch)` options in force for the current
/// integration (None: default).
#[pyfunction]
pub(crate) fn _integration_methods() -> (Option<bool>, Option<bool>, Option<bool>) {
    let m = host_methods::methods();
    (m.risch, m.manual, m.heurisch)
}

/// Record a held integral node as a `NonElementaryIntegral`.
#[pyfunction]
pub(crate) fn _mark_nonelementary(e: &Bound<'_, PyAny>) -> PyResult<()> {
    host_methods::mark_nonelementary(&to_expr(e)?);
    Ok(())
}

/// Whether `e` is a node recorded as a `NonElementaryIntegral`.
#[pyfunction]
pub(crate) fn _is_nonelementary(e: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(host_methods::is_nonelementary(&to_expr(e)?))
}

/// The `NonElementaryIntegral` object of a recorded node.
pub(crate) fn nonelementary_object(py: Python<'_>, e: &Expr) -> PyResult<Option<PyObject>> {
    if !matches!(e.kind(), Kind::Func(Func::Integral, _)) || !host_methods::is_nonelementary(e) {
        return Ok(None);
    }
    let class = py.import("symox.integrals.risch")?.getattr("NonElementaryIntegral")?;
    let value = Bound::new(py, PyExprNode { e: e.clone() })?;
    Ok(Some(class.call_method1("_from_native", (value,))?.unbind()))
}
