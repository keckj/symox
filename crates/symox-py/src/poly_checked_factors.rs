//! Checked typed factor views; no expression-printer round trips.
use super::{division, err, w, PolyValue, PyPoly};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use symox::poly::{checked_factors as native, factor_budget::Limits};
pub(super) fn list<'py>(
    p: &PolyValue,
    py: Python<'py>,
    squarefree: bool,
    all: bool,
    polys: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    let factors = native::compatible_list(p, squarefree, all, Limits::default()).map_err(err)?;
    render_list(factors, symox::Expr::one(), py, polys)
}
fn render_list<'py>(
    factors: native::Factorization,
    denominator: symox::Expr,
    py: Python<'py>,
    polys: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    let unit = w(symox::ops::div(&factors.unit.as_expr(), &denominator))
        .into_pyobject(py)?
        .into_any();
    let result = pyo3::types::PyList::empty(py);
    for (p, k) in factors.factors {
        let value = if polys {
            PyPoly { p }.into_pyobject(py)?.into_any().unbind()
        } else {
            w(p.as_expr()).into_pyobject(py)?.into_any().unbind()
        };
        result.append(PyTuple::new(
            py,
            [value, k.into_pyobject(py)?.into_any().unbind()],
        )?)?;
    }
    PyTuple::new(py, [unit, result.into_any()])
}
/// Convert one multiplicative component with the caller's options, as the
/// reference `_poly_from_expr(base, opt)`; `None` is a polification failure.
fn component_builder<'a, 'py>(
    py: Python<'py>,
    gens: &'a Bound<'py, PyTuple>,
    kwargs: Option<&'a Bound<'py, PyDict>>,
    failure: &'a std::cell::RefCell<Option<PyErr>>,
) -> impl Fn(&symox::Expr) -> symox::poly::convert::PResult<Option<PolyValue>> + 'a
where
    'py: 'a,
{
    move |base: &symox::Expr| {
        let attempt = || -> PyResult<Option<PolyValue>> {
            let value = w(base.clone()).into_pyobject(py)?.into_any();
            match prepare(&value, gens, None, kwargs, false) {
                Ok((p, _)) => Ok(Some(p)),
                Err(error) => {
                    let failed = py
                        .import("symox.polys.polyerrors")?
                        .getattr("PolificationFailed")?;
                    if error.is_instance(py, &failed) {
                        Ok(None)
                    } else {
                        Err(error)
                    }
                }
            }
        };
        attempt().map_err(|error| {
            let message = error.to_string();
            failure.borrow_mut().get_or_insert(error);
            symox::poly::convert::PolyError::Failed(message)
        })
    }
}
fn with_failure<T>(
    failure: std::cell::RefCell<Option<PyErr>>,
    result: symox::poly::convert::PResult<T>,
) -> PyResult<T> {
    match (failure.into_inner(), result) {
        (Some(error), Err(_)) => Err(error),
        (_, result) => result.map_err(err),
    }
}
/// The reference `_generic_factor_list` for expressions whose options select
/// an algebraic coefficient domain: every multiplicative component of the
/// expression is converted and factored on its own.
fn domain_list<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    squarefree: bool,
    polys: bool,
) -> PyResult<PyObject> {
    let py = f.py();
    let expression = crate::to_expr(f)?;
    let failure = std::cell::RefCell::new(None);
    let result = {
        let build = component_builder(py, gens, kwargs, &failure);
        symox::poly::api::domain_factor_list(&expression, squarefree, &build, Limits::default())
    };
    let (coefficient, factors) = with_failure(failure, result)?;
    let list = pyo3::types::PyList::empty(py);
    for (p, k) in factors {
        let value = if polys {
            PyPoly { p }.into_pyobject(py)?.into_any().unbind()
        } else {
            w(p.as_expr()).into_pyobject(py)?.into_any().unbind()
        };
        list.append(PyTuple::new(py, [value, k.into_pyobject(py)?.into_any().unbind()])?)?;
    }
    Ok(PyTuple::new(py, [w(coefficient).into_pyobject(py)?.into_any(), list.into_any()])?
        .into_any()
        .unbind())
}
/// The reference `_symbolic_factor` for expressions whose options select an
/// algebraic coefficient domain.
pub(super) fn domain_expression<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    squarefree: bool,
) -> PyResult<symox::Expr> {
    let py = f.py();
    let expression = crate::to_expr(f)?;
    let failure = std::cell::RefCell::new(None);
    let result = {
        let build = component_builder(py, gens, kwargs, &failure);
        symox::poly::api::domain_factor(&expression, squarefree, &build, Limits::default())
    };
    with_failure(failure, result)
}
pub(super) fn part(p: &PolyValue) -> PyResult<PyPoly> {
    Ok(PyPoly {
        p: native::squarefree_part(p, Limits::default()).map_err(err)?,
    })
}
fn prepare<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    modulus: Option<u64>,
    kwargs: Option<&Bound<'py, PyDict>>,
    part: bool,
) -> PyResult<(PolyValue, bool)> {
    let kw = kwargs.map_or_else(|| Ok(PyDict::new(f.py())), |kw| kw.copy())?;
    if let Some(value) = kw.get_item("frac")? {
        if part {
            return Err(division::FlagError::new_err(
                "'frac' flag is not allowed in sqf_part",
            ));
        }
        division::boolean(&value, "frac")?;
        kw.del_item("frac")?;
    }
    for flag in ["all", "include"] {
        if kw.contains(flag)? {
            return Err(division::FlagError::new_err(format!(
                "'{flag}' flag is not allowed in this context"
            )));
        }
    }
    if let Some(modulus) = modulus {
        kw.set_item("modulus", modulus)?;
    }
    division::validate_unary(gens, Some(&kw))?;
    let expression = crate::to_expr(f)?;
    symox::poly::factor_preflight::check(&expression, Limits::default()).map_err(err)?;
    // List factorization accepts numeric constants without explicit generators.
    // This temporary, fresh generator cannot escape: constants have no factors.
    if !part
        && expression.is_number()
        && gens.is_empty()
        && kw.get_item("gens")?.is_none_or(|g| g.is_none())
    {
        let generator = symox::expr::Expr::dummy("x", Default::default());
        let args = PyTuple::new(f.py(), [w(generator)])?;
        return division::prepare_unary(f, &args, Some(&kw));
    }
    match division::prepare_unary(f, gens, Some(&kw)) {
        Err(error) if !part && !expression.is_number()
            && f.downcast::<PyPoly>().is_err()
            && division::option_generators(gens, Some(&kw))?.is_empty()
            && error.is_instance_of::<pyo3::exceptions::PyValueError>(f.py()) => {
            // Expansion can cancel every inferred generator. This is a
            // failed polynomial conversion, not a factorization failure;
            // callers use PolificationFailed to select a scalar fallback.
            let expanded = symox::poly::factor_budget::Session::new(Limits::default())
                .run(|| Ok(symox::expand::expand(&expression))).map_err(err)?;
            if !expanded.is_number() { return Err(error); }
            let py = f.py();
            let options = py.import("symox.polys.polyoptions")?.getattr("Options")?
                .call1((gens, &kw))?;
            let value = w(expanded).into_pyobject(py)?;
            let exception = py.import("symox.polys.polyerrors")?
                .getattr("PolificationFailed")?.call1((options, &value, &value))?;
            Err(PyErr::from_value(exception))
        }
        result => result,
    }
}
pub(super) fn global_list<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    modulus: Option<u64>,
    kwargs: Option<&Bound<'py, PyDict>>,
    squarefree: bool,
) -> PyResult<PyObject> {
    if let Some(value) = kwargs.map(|kw| kw.get_item("frac")).transpose()?.flatten() {
        if division::boolean(&value, "frac")? == Some(true) {
            return global_fraction_list(f, gens, modulus, kwargs, squarefree);
        }
    }
    let (p, polys) = prepare(f, gens, modulus, kwargs, false)?;
    let mut symbolic = f.downcast::<PyPoly>().is_err() && modulus.is_none();
    if let Some(kw) = kwargs {
        for name in ["domain", "extension", "modulus"] {
            if kw.get_item(name)?.is_some_and(|v| !v.is_none()) {
                symbolic = false;
            }
        }
        if kw
            .get_item("field")?
            .is_some_and(|v| division::boolean(&v, "field").is_ok_and(|v| v == Some(true)))
        {
            symbolic = false;
        }
    }
    if symbolic {
        let expression = crate::to_expr(f)?;
        let scalar = match expression.kind() {
            symox::expr::Kind::Num(_) | symox::expr::Kind::Const(symox::expr::Constant::I) => true,
            symox::expr::Kind::Add(_, terms) => terms
                .iter()
                .all(|(e, _)| e.is_constant(symox::expr::Constant::I)),
            symox::expr::Kind::Pow(b, e) => b.is_number() && e.is_number(),
            symox::expr::Kind::Mul(_, factors) => factors.iter().all(|(b, e)| {
                (b.is_number() && e.is_number())
                    || (b.is_constant(symox::expr::Constant::I) && e.is_one())
            }),
            _ => false,
        };
        if scalar {
            return Ok(PyTuple::new(
                f.py(),
                [
                    w(expression).into_pyobject(f.py())?.into_any(),
                    pyo3::types::PyList::empty(f.py()).into_any(),
                ],
            )?
            .into_any()
            .unbind());
        }
    }
    let negative = |e: &symox::Expr| match e.kind() {
        symox::expr::Kind::Num(n) | symox::expr::Kind::Mul(n, _) => n.is_negative(),
        _ => false,
    };
    if p.true_gens().iter().any(|g| match g.kind() {
        symox::expr::Kind::Pow(_, exponent) => negative(exponent),
        symox::expr::Kind::Func(symox::expr::Func::Exp, args) => args.first().is_some_and(negative),
        _ => false,
    }) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "a polynomial expected; denominator factors require frac=True",
        ));
    }
    if (p.algebraic.is_some() || p.expression.is_some()) && f.downcast::<PyPoly>().is_err() && modulus.is_none() {
        return domain_list(f, gens, kwargs, squarefree, polys);
    }
    if symbolic {
        if squarefree {
            let (factors, denominator) =
                native::symbolic_list(&p, true, Limits::default()).map_err(err)?;
            return Ok(render_list(factors, denominator, f.py(), polys)?
                .into_any()
                .unbind());
        }
        let mut component_gens = super::gens_of(gens)?;
        if component_gens.is_empty() {
            if let Some(g) = kwargs
                .map(|kw| kw.get_item("gens"))
                .transpose()?
                .flatten()
                .filter(|g| !g.is_none())
            {
                component_gens = if let Ok(g) = g.downcast::<PyTuple>() {
                    super::gens_of(g)?
                } else if let Ok(g) = g.downcast::<pyo3::types::PyList>() {
                    g.iter()
                        .map(|g| crate::to_expr(&g))
                        .collect::<PyResult<_>>()?
                } else {
                    vec![crate::to_expr(&g)?]
                };
            }
        }
        let (coefficient, factors) = symox::poly::api::checked_natural_factor_list(
            &crate::to_expr(f)?,
            &component_gens,
            squarefree,
            &p,
            Limits::default(),
        )
        .map_err(err)?;
        let result = pyo3::types::PyList::empty(f.py());
        for (p, k) in factors {
            // The source component fixes both main generators and coefficient
            // parameters. Attach their exact descriptor without re-inferring
            // generators from the final factor expression.
            let p = super::domains::normalize_poly(p)?;
            let value = if polys {
                PyPoly { p }.into_pyobject(f.py())?.into_any().unbind()
            } else {
                w(p.as_expr()).into_pyobject(f.py())?.into_any().unbind()
            };
            result.append(PyTuple::new(
                f.py(),
                [value, k.into_pyobject(f.py())?.into_any().unbind()],
            )?)?;
        }
        return Ok(PyTuple::new(
            f.py(),
            [
                w(coefficient).into_pyobject(f.py())?.into_any(),
                result.into_any(),
            ],
        )?
        .into_any()
        .unbind());
    }
    Ok(list(&p, f.py(), squarefree, false, polys)?
        .into_any()
        .unbind())
}
pub(super) fn global_part<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyObject> {
    let (p, polys) = prepare(f, gens, None, kwargs, true)?;
    let value = part(&p)?;
    if polys {
        Ok(value.into_pyobject(f.py())?.into_any().unbind())
    } else {
        Ok(w(value.p.as_expr())
            .into_pyobject(f.py())?
            .into_any()
            .unbind())
    }
}

fn is_natural(modulus: Option<u64>, kwargs: &Bound<'_, PyDict>) -> PyResult<bool> {
    if modulus.is_some() {
        return Ok(false);
    }
    for name in ["domain", "extension", "modulus"] {
        if kwargs.get_item(name)?.is_some_and(|v| !v.is_none()) {
            return Ok(false);
        }
    }
    Ok(kwargs
        .get_item("field")?
        .is_none_or(|v| division::boolean(&v, "field").ok() != Some(Some(true))))
}
fn prepare_fraction_side(
    session: &mut symox::poly::factor_budget::Session,
    p: PolyValue,
    expression: &symox::Expr,
    generators: &[symox::Expr],
    natural: bool,
    squarefree: bool,
) -> PyResult<PolyValue> {
    if natural
        && p.true_gens()
            .iter()
            .any(|g| matches!(g.kind(),symox::expr::Kind::Pow(_,e) if !e.is_integer_number()))
    {
        return Err(pyo3::exceptions::PyNotImplementedError::new_err("natural fractional factor lists with noninteger or symbolic power multiplicities are not implemented; use a Poly with explicit formal generators"));
    }
    if natural && squarefree && generators.is_empty() && p.true_gens().len() > 1 {
        let generators = session
            .run(|| {
                symox::poly::factor_fraction::component_generators(expression, Limits::default())
            })
            .map_err(err)?;
        return super::domains::PyDomain::of(&p)?.poly(expression, &generators);
    }
    Ok(p)
}
fn fraction_side(
    p: &PolyValue,
    expression: &symox::Expr,
    generators: &[symox::Expr],
    natural: bool,
    squarefree: bool,
) -> symox::poly::convert::PResult<(symox::Expr, Vec<(PolyValue, u32)>)> {
    let limits = Limits::default();
    if natural {
        if squarefree {
            let (factors, denominator) = native::symbolic_list(p, true, limits)?;
            return Ok((
                symox::ops::div(&factors.unit.as_expr(), &denominator),
                factors.factors,
            ));
        }
        let (c, factors) = symox::poly::api::checked_natural_factor_list(
            expression, generators, false, p, limits,
        )?;
        let factors = factors
            .into_iter()
            .map(|(p, k)| Ok((native::normalize_natural_factor(p)?, k)))
            .collect::<symox::poly::convert::PResult<_>>()?;
        return Ok((c, factors));
    }
    let factors = native::compatible_list(p, squarefree, false, limits)?;
    Ok((factors.unit.as_expr(), factors.factors))
}
fn render_factors(
    py: Python<'_>,
    factors: Vec<(PolyValue, u32)>,
    polys: bool,
) -> PyResult<PyObject> {
    let result = pyo3::types::PyList::empty(py);
    for (p, k) in factors {
        let value = if polys {
            PyPoly { p }.into_pyobject(py)?.into_any().unbind()
        } else {
            w(p.as_expr()).into_pyobject(py)?.into_any().unbind()
        };
        result.append(PyTuple::new(
            py,
            [value, k.into_pyobject(py)?.into_any().unbind()],
        )?)?;
    }
    Ok(result.into_any().unbind())
}

fn global_fraction_list<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    modulus: Option<u64>,
    kwargs: Option<&Bound<'py, PyDict>>,
    squarefree: bool,
) -> PyResult<PyObject> {
    // User conversions, Python domain selection and rendering run with no
    // ambient cancellation. Only native stages activate the same transferable
    // allowance; no private unwind can cross a Python callback or trampoline.
    let py = f.py();
    let limits = Limits::default();
    let mut session = symox::poly::factor_budget::Session::new(limits);
    let kw = kwargs.map_or_else(|| Ok(PyDict::new(py)), |kw| kw.copy())?;
    kw.del_item("frac")?;
    for flag in ["all", "include"] {
        if kw.contains(flag)? {
            return Err(division::FlagError::new_err(format!(
                "'{flag}' flag is not allowed in this context"
            )));
        }
    }
    if let Some(modulus) = modulus {
        kw.set_item("modulus", modulus)?;
    }
    division::validate_unary(gens, Some(&kw))?;
    // A Poly is already a formal polynomial. Its generator may itself be
    // 1/x; interpreting its expression as a rational function loses that
    // contract and its carried descriptor/alias.
    if f.downcast::<PyPoly>().is_ok() {
        let (p, polys) = prepare(f, gens, modulus, Some(&kw), false)?;
        let factors = session
            .run(|| native::compatible_list(&p, squarefree, false, limits))
            .map_err(err)?;
        let c = w(factors.unit.as_expr())
            .into_pyobject(py)?
            .into_any()
            .unbind();
        return Ok(PyTuple::new(
            py,
            [
                c,
                render_factors(py, factors.factors, polys)?,
                pyo3::types::PyList::empty(py).into_any().unbind(),
            ],
        )?
        .into_any()
        .unbind());
    }
    let generators = division::option_generators(gens, Some(&kw))?;
    // Normalize positional/keyword containers once; gens=() must
    // still allow constants without a public artificial generator.
    let parsed_gens = PyTuple::new(py, generators.iter().cloned().map(w))?;
    if kw.contains("gens")? {
        kw.del_item("gens")?;
    }
    let gens = &parsed_gens;
    let expression = crate::to_expr(f)?;
    let parts = session
        .run(|| symox::poly::factor_fraction::split(&expression, &generators, limits))
        .map_err(err)?;
    let ((un, n), (ud, d)) = session
        .run(|| {
            Ok((
                symox::poly::factor_fraction::scalar_content(&parts.numerator, limits)?,
                symox::poly::factor_fraction::scalar_content(&parts.denominator, limits)?,
            ))
        })
        .map_err(err)?;
    let natural = is_natural(modulus, &kw)?;
    let nobj = w(n.clone()).into_pyobject(py)?.into_any();
    let dobj = w(d.clone()).into_pyobject(py)?.into_any();
    let (polyn, polyd, polys) = if natural {
        // In natural mode each side (and each source factor) retains its
        // own inferred generators, as for ordinary symbolic factor lists.
        let (a, polys) = prepare(&nobj, gens, modulus, Some(&kw), false)?;
        let (b, _) = prepare(&dobj, gens, modulus, Some(&kw), false)?;
        (a, b, polys)
    } else {
        let temporary = PyTuple::new(py, [w(symox::Expr::dummy("x", Default::default()))])?;
        let args = if generators.is_empty() && n.is_number() && d.is_number() {
            &temporary
        } else {
            gens
        };
        let (mut a, mut b, polys) = division::prepare_pair(&nobj, &dobj, args, Some(&kw))?;
        if generators.is_empty() && !(n.is_number() && d.is_number()) {
            // Select extensions and coefficient parameters jointly, then
            // infer each side's main generators against that same domain.
            // No domain string is used as identity or to recover symbols.
            let domain = super::domains::PyDomain::of(&a)?;
            let own = |e: &symox::Expr| -> PyResult<PolyValue> {
                let mut gs = if e.is_number() {
                    vec![]
                } else {
                    super::Poly::from_expr(e, &[])
                        .map_err(err)?
                        .true_gens()
                        .to_vec()
                };
                gs.retain(|g| domain.inner.convert(g).is_err());
                if gs.is_empty() {
                    gs.push(symox::Expr::dummy("x", Default::default()));
                }
                domain.poly(e, &gs)
            };
            a = own(&n)?;
            b = own(&d)?;
        }
        (a, b, polys)
    };
    // Coercion can turn a nonzero literal denominator into zero in GF(p).
    // Test this before either factor engine or coefficient division.
    if polyd.is_zero() {
        return Err(pyo3::exceptions::PyZeroDivisionError::new_err(
            "fraction denominator is zero in the selected coefficient domain",
        ));
    }
    let polyn = prepare_fraction_side(&mut session, polyn, &n, &generators, natural, squarefree)?;
    let polyd = prepare_fraction_side(&mut session, polyd, &d, &generators, natural, squarefree)?;
    let ((cn, nf), (cd, df)) = session
        .run(|| {
            Ok((
                fraction_side(&polyn, &n, &generators, natural, squarefree)?,
                fraction_side(&polyd, &d, &generators, natural, squarefree)?,
            ))
        })
        .map_err(err)?;
    let domain = if natural {
        None
    } else {
        Some(super::domains::PyDomain::of(&polyn)?.inner.field())
    };
    let coefficient=session.run(||{
        if let Some(domain)=domain {
            let(nscalar,nexternal)=symox::ops::as_coeff_mul(&un);
            let(dscalar,dexternal)=symox::ops::as_coeff_mul(&ud);
            let numerator=domain.convert(&symox::ops::mul(&[symox::Expr::num(nscalar),cn]))?;
            let denominator=domain.convert(&symox::ops::mul(&[symox::Expr::num(dscalar),cd]))?;
            if denominator.is_zero(){return Err(symox::poly::convert::PolyError::Unsupported("division by zero: fraction denominator is zero in the selected coefficient domain".into()));}
            let quotient=numerator.exquo(&denominator)?.as_expr();
            Ok(symox::ops::mul(&[quotient,symox::ops::div(&nexternal,&dexternal)]))
        }else{Ok(symox::ops::div(&symox::ops::mul(&[un,cn]),&symox::ops::mul(&[ud,cd])))}
    }).map_err(err)?;
    Ok(PyTuple::new(
        py,
        [
            w(coefficient).into_pyobject(py)?.into_any().unbind(),
            render_factors(py, nf, polys)?,
            render_factors(py, df, polys)?,
        ],
    )?
    .into_any()
    .unbind())
}

pub(super) fn ground(p: &PolyValue, py: Python<'_>) -> PyResult<PyObject> {
    let result = PyDict::new(py);
    for (root, multiplicity) in native::ground_roots(p, Limits::default()).map_err(err)? {
        result.set_item(w(root), multiplicity)?;
    }
    Ok(result.into_any().unbind())
}
pub(super) fn include(p: &PolyValue, squarefree: bool, all: bool) -> PyResult<Vec<(PyPoly, u32)>> {
    Ok(native::include(p, squarefree, all, Limits::default())
        .map_err(err)?
        .into_iter()
        .map(|(p, k)| (PyPoly { p }, k))
        .collect())
}
pub(super) fn factorial(p: &PolyValue) -> PyResult<Vec<(PyPoly, u32)>> {
    Ok(symox::poly::factorial::gff_list(p)
        .map_err(err)?
        .into_iter()
        .map(|(p, k)| (PyPoly { p }, k))
        .collect())
}
#[pyfunction]
#[pyo3(signature=(f,*gens,**kwargs))]
fn gff_list<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyObject> {
    let (p, polys) = division::prepare_unary(f, gens, kwargs)?;
    let result = pyo3::types::PyList::empty(f.py());
    for (p, k) in symox::poly::factorial::gff_list(&p).map_err(err)? {
        let value = if polys {
            PyPoly { p }.into_pyobject(f.py())?.into_any().unbind()
        } else {
            w(p.as_expr()).into_pyobject(f.py())?.into_any().unbind()
        };
        result.append(PyTuple::new(
            f.py(),
            [value, k.into_pyobject(f.py())?.into_any().unbind()],
        )?)?;
    }
    Ok(result.into_any().unbind())
}
#[pyfunction]
#[pyo3(signature=(f,*gens,**kwargs))]
fn ground_roots<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyObject> {
    if kwargs.is_some_and(|kw| kw.contains("polys").unwrap_or(false)) {
        return Err(division::FlagError::new_err(
            "'polys' flag is not allowed in ground_roots",
        ));
    }
    let (p, _) = division::prepare_unary(f, gens, kwargs)?;
    ground(&p, f.py())
}
pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ground_roots, m)?)?;
    m.add_function(wrap_pyfunction!(gff_list, m)?)?;
    Ok(())
}
