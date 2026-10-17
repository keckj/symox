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
        let frac = division::boolean(&value, "frac")?.unwrap_or(false);
        kw.del_item("frac")?;
        if frac {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err("frac=True factor lists are not implemented; factor numerator and denominator explicitly"));
        }
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
    division::prepare_unary(f, gens, Some(&kw))
}
pub(super) fn global_list<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    modulus: Option<u64>,
    kwargs: Option<&Bound<'py, PyDict>>,
    squarefree: bool,
) -> PyResult<PyObject> {
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
        return Err(pyo3::exceptions::PyValueError::new_err("a polynomial expected; denominator factors require frac=True, which is not implemented"));
    }
    if symbolic {
        if squarefree {
            let(factors,denominator)=native::symbolic_list(&p,true,Limits::default()).map_err(err)?;
            return Ok(render_list(factors,denominator,f.py(),polys)?.into_any().unbind());
        }
        let mut component_gens=super::gens_of(gens)?;
        if component_gens.is_empty(){if let Some(g)=kwargs.map(|kw|kw.get_item("gens")).transpose()?.flatten().filter(|g|!g.is_none()){
            component_gens=if let Ok(g)=g.downcast::<PyTuple>(){super::gens_of(g)?}else if let Ok(g)=g.downcast::<pyo3::types::PyList>(){g.iter().map(|g|crate::to_expr(&g)).collect::<PyResult<_>>()?}else{vec![crate::to_expr(&g)?]};
        }}
        let(coefficient,factors)=symox::poly::api::checked_natural_factor_list(&crate::to_expr(f)?,&component_gens,squarefree,&p,Limits::default()).map_err(err)?;
        let result=pyo3::types::PyList::empty(f.py());
        for(p,k)in factors{
            // The source component fixes both main generators and coefficient
            // parameters. Attach their exact descriptor without re-inferring
            // generators from the final factor expression.
            let p=super::domains::normalize_poly(p)?;
            let value=if polys{PyPoly{p}.into_pyobject(f.py())?.into_any().unbind()}else{w(p.as_expr()).into_pyobject(f.py())?.into_any().unbind()};
            result.append(PyTuple::new(f.py(),[value,k.into_pyobject(f.py())?.into_any().unbind()])?)?;
        }
        return Ok(PyTuple::new(f.py(),[w(coefficient).into_pyobject(f.py())?.into_any(),result.into_any()])?.into_any().unbind());
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
