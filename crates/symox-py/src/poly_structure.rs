//! Python option handling for sparse structural polynomial views.
use super::{division, w, PolyValue};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyString, PyTuple, PyType};
use symox::Expr;

pub(super) fn polynomial_error(py: Python<'_>, kind: &str, message: String) -> PyErr {
    let class = py
        .import("symox")
        .and_then(|m| m.getattr(kind))
        .and_then(|e| e.downcast_into::<PyType>().map_err(Into::into));
    match class {
        Ok(class) => PyErr::from_type(class, (message,)),
        Err(e) => e,
    }
}
fn label(g: &Expr) -> PyResult<String> {
    if let Some(name) = g.as_symbol() {
        if name.len() > 100_000 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "polynomial structure resource limit: generator label",
            ));
        }
    } else {
        symox::poly::factor_preflight::check(g, symox::poly::factor_budget::Limits::default())
            .map_err(super::err)?;
        let mut stack = vec![g.clone()];
        let mut visits = 0usize;
        let mut text = 0usize;
        while let Some(g) = stack.pop() {
            visits += 1;
            if visits > 4096 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "polynomial structure resource limit: unfolded generator label",
                ));
            }
            text=text.saturating_add(32+match g.kind(){
                symox::Kind::Sym(s)=>s.name().len(),
                symox::Kind::Num(n)=>n.to_string().len(),
                symox::Kind::Func(f,_)=>{
                    if matches!(f,symox::expr::Func::MatrixElement(_)|symox::expr::Func::MatrixProxy(_)|symox::expr::Func::MatrixTrace(_)|symox::expr::Func::MatrixDeterminant(_)|symox::expr::Func::MatrixDotProduct(..)){
                        return Err(pyo3::exceptions::PyNotImplementedError::new_err("default generator sorting of matrix-backed atoms is not implemented; supply an explicit permutation"));
                    }
                    f.name().len()
                },_=>0,
            });
            if text > 1_000_000 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "polynomial structure resource limit: unfolded generator label",
                ));
            }
            stack.extend(g.args());
        }
    }
    Ok(symox::printer::to_str(g))
}
enum Priority {
    Identity(Expr),
    Name(String),
}
fn priority(value: &Bound<'_, PyAny>) -> PyResult<Priority> {
    if let Ok(name) = value.extract::<String>() {
        return Ok(Priority::Name(name));
    }
    Ok(Priority::Identity(crate::to_expr(value)?))
}
fn values<'py>(value: &Bound<'py, PyAny>) -> PyResult<Option<Vec<Bound<'py, PyAny>>>> {
    if let Ok(xs) = value.downcast::<PyTuple>() {
        if xs.len() > 4096 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "polynomial structure resource limit: ordering priorities",
            ));
        }
        return Ok(Some(xs.iter().collect()));
    }
    if let Ok(xs) = value.downcast::<PyList>() {
        if xs.len() > 4096 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "polynomial structure resource limit: ordering priorities",
            ));
        }
        return Ok(Some(xs.iter().collect()));
    }
    Ok(None)
}
// Python's regex \d in the SymPy oracle means Unicode decimal digits. Use the
// standard-library Unicode database only for non-ASCII suffix characters;
// all exponent transformations and sorting remain native and SymPy-free.
fn numeric_suffix(py: Python<'_>, name: &str) -> PyResult<(String, String)> {
    let mut start = name.len();
    let mut reverse = Vec::new();
    for (i, c) in name.char_indices().rev() {
        let digit = if c.is_ascii_digit() {
            Some(c)
        } else if c.is_ascii() {
            None
        } else {
            let result = py
                .import("unicodedata")?
                .call_method1("decimal", (c.to_string(), py.None()))?;
            if result.is_none() {
                None
            } else {
                Some((b'0' + result.extract::<u8>()?) as char)
            }
        };
        if let Some(digit) = digit {
            start = i;
            reverse.push(digit);
        } else {
            break;
        }
    }
    let digits: String = reverse.into_iter().rev().collect();
    Ok((
        name[..start].to_owned(),
        digits.trim_start_matches('0').to_owned(),
    ))
}
pub(super) fn reordered(
    p: &PolyValue,
    gens: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PolyValue> {
    let mut sort = Vec::<String>::new();
    let mut wrt = Vec::<Priority>::new();
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs {
            let key = key.extract::<String>()?;
            if let Ok(s) = value.downcast::<PyString>() {
                if s.to_cow()?.len() > 1_000_000 {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "polynomial structure resource limit: ordering text",
                    ));
                }
            }
            match key.as_str(){
                "sort" if !value.is_none()=>{
                    sort=if let Ok(s)=value.extract::<String>(){s.split('>').map(|s|s.trim().to_owned()).collect()}
                    else if let Some(xs)=values(&value)?{xs.into_iter().map(|v|match priority(&v)?{Priority::Name(n)=>Ok(n),Priority::Identity(g)=>label(&g)}).collect::<PyResult<_>>()?}
                    else{return Err(division::OptionError::new_err("invalid argument for 'sort' option"));};
                },
                "wrt" if !value.is_none()=>{
                    wrt=if let Ok(s)=value.extract::<String>(){
                        if s.trim_end().ends_with(','){return Err(division::OptionError::new_err("Bad input: missing parameter."));}
                        s.split(|c:char|c==','||c.is_whitespace()).filter(|s|!s.is_empty()).map(|s|Priority::Name(s.to_owned())).collect()
                    }else if let Some(xs)=values(&value)?{xs.into_iter().map(|v|priority(&v)).collect::<PyResult<_>>()?}
                    else if value.downcast::<crate::PyExprNode>().is_ok(){vec![priority(&value)?]}
                    else{return Err(division::OptionError::new_err("invalid argument for 'wrt' option"));};
                },
                "order" if !value.is_none()=>{
                    let order=value.extract::<String>()?;
                    if symox::poly::groebner::MonomialOrder::parse(&order).is_none(){return Err(division::OptionError::new_err("order must be lex, grlex or grevlex"));}
                    // As in SymPy this validates a polynomial option; reorder
                    // permutes generators, not the sparse storage invariant.
                },
                "sort"|"wrt"|"order"=>{},
                "domain"|"field"|"extension"|"modulus"|"gens"|"greedy"|"polys"=>return Err(pyo3::exceptions::PyNotImplementedError::new_err(format!("'{key}' does not change a polynomial's generator order; supply an explicit generator permutation"))),
                _=>return Err(division::OptionError::new_err(format!("'{key}' is not a valid option"))),
            }
        }
    }
    let mut requested_count = 0usize;
    for value in gens {
        let n = if let Ok(t) = value.downcast::<PyTuple>() {
            t.len()
        } else if let Ok(t) = value.downcast::<PyList>() {
            t.len()
        } else {
            1
        };
        requested_count = requested_count.saturating_add(n);
    }
    if requested_count != 0 && requested_count != p.ngens {
        return Err(polynomial_error(
            gens.py(),
            "PolynomialError",
            "generators list can differ only up to order of elements".into(),
        ));
    }
    let requested = super::gens_of(gens)?;
    let indices = if !requested.is_empty() {
        if requested.len() != p.ngens
            || requested
                .iter()
                .enumerate()
                .any(|(j, g)| !p.true_gens().contains(g) || requested[..j].contains(g))
        {
            return Err(polynomial_error(
                gens.py(),
                "PolynomialError",
                "generators list can differ only up to order of elements".into(),
            ));
        }
        requested
            .iter()
            .map(|g| p.true_gens().iter().position(|x| x == g).unwrap())
            .collect::<Vec<_>>()
    } else {
        if p.ngens > 4096 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "polynomial structure resource limit: generator sort keys",
            ));
        }
        let mut labels = Vec::with_capacity(p.ngens);
        let mut bytes = 0usize;
        for g in p.true_gens() {
            let name = label(g)?;
            bytes = bytes.saturating_add(name.len());
            if bytes > 1_000_000 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "polynomial structure resource limit: generator labels",
                ));
            }
            labels.push(name);
        }
        let mut keys = Vec::with_capacity(p.ngens);
        for (i, (g, name)) in p.true_gens().iter().zip(&labels).enumerate() {
            let (base, index) = numeric_suffix(gens.py(), name)?;
            let rank = if let Some(j) = wrt.iter().position(|w| match w {
                Priority::Identity(e) => e == g,
                Priority::Name(n) => n == name,
            }) {
                j as i64 - wrt.len() as i64
            } else {
                sort.iter().position(|n| n == &base).map_or_else(
                    || symox::poly::convert::gens_order(&base) as i64,
                    |j| j as i64 + 1,
                )
            };
            keys.push(((rank, base.to_owned(), index.len(), index.to_owned()), i));
        }
        keys.sort_by(|a, b| a.0.cmp(&b.0));
        keys.into_iter().map(|(_, i)| i).collect()
    };
    symox::poly::structure::reorder(p, &indices).map_err(super::err)
}

pub(super) fn has_only(p: &PolyValue, gens: &Bound<'_, PyTuple>) -> PyResult<bool> {
    // Unlike reorder these are separate generator arguments, not a positional
    // list option. Duplicates are harmless; unknown identities are errors.
    let mut indices = vec![];
    for gen in gens {
        let g = crate::to_expr(&gen).map_err(|_| {
            polynomial_error(
                gens.py(),
                "GeneratorsError",
                "argument is not a polynomial generator".into(),
            )
        })?;
        let i = p.true_gens().iter().position(|x| x == &g).ok_or_else(|| {
            polynomial_error(
                gens.py(),
                "GeneratorsError",
                format!("polynomial doesn't have {g} as generator"),
            )
        })?;
        if !indices.contains(&i) {
            indices.push(i);
        }
    }
    symox::poly::structure::has_only_indices(p, &indices).map_err(super::err)
}

pub(super) fn dense<'py>(
    p: &PolyValue,
    py: Python<'py>,
    coefficients: bool,
) -> PyResult<Vec<Bound<'py, PyTuple>>> {
    if p.ngens != 1 {
        return Err(polynomial_error(
            py,
            "PolynomialError",
            "multivariate polynomials not supported".into(),
        ));
    }
    symox::poly::structure::dense_terms(p)
        .map_err(super::err)?
        .into_iter()
        .map(|(m, c)| {
            let monom = PyTuple::new(py, m)?;
            if coefficients {
                PyTuple::new(py, [monom.into_any(), w(c).into_pyobject(py)?.into_any()])
            } else {
                Ok(monom)
            }
        })
        .collect()
}

/// Dense integer coefficients (constant term first) of the `n`-th cyclotomic
/// polynomial, computed natively by the product formula
/// `Phi_{pm}(x) = Phi_m(x^p)/Phi_m(x)` on sparse divisors.
#[pyfunction]
pub(super) fn _cyclotomic_coefficients(py: Python<'_>, n: u64) -> PyResult<Vec<PyObject>> {
    // Memory bound only: the dense image of x^n stays below 2^22 slots.
    if n == 0 || n > 1 << 22 {
        return Err(pyo3::exceptions::PyValueError::new_err(if n == 0 {
            "Cannot generate cyclotomic polynomial of order 0".to_string()
        } else {
            "polynomial structure resource limit: cyclotomic order".to_string()
        }));
    }
    let mut memo = Default::default();
    let int = py.import("builtins")?.getattr("int")?;
    symox::poly::factor_uni::cyclotomic(n, &mut memo)
        .into_iter()
        .map(|c| Ok(int.call1((c.to_string(),))?.unbind()))
        .collect()
}

/// Descending coefficients and SymPy domain name of an orthogonal polynomial
/// (`symox::poly::orthopolys`); the domain is `None` when it must be inferred.
#[pyfunction]
pub(super) fn _orthopoly_coefficients(
    family: &str,
    n: usize,
    params: Vec<Bound<'_, PyAny>>,
) -> PyResult<(Vec<super::PyExpr>, Option<String>)> {
    use symox::poly::orthopolys::{named_poly, Family};
    let family = Family::parse(family)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("unknown polynomial family {family:?}")))?;
    let params = params.iter().map(super::to_expr).collect::<PyResult<Vec<_>>>()?;
    let named = named_poly(family, n, &params).map_err(super::err)?;
    Ok((named.coeffs.into_iter().map(w).collect(), named.domain))
}

/// Dense ascending rational coefficients of a univariate `ZZ`/`QQ` polynomial
/// without parameters (`None` for every other representation).
fn plain_dense(p: &PolyValue) -> Option<(Vec<symox::poly::Rat>, symox::poly::polyobj::Domain)> {
    use symox::poly::polyobj::Domain;
    let base = &p.base;
    if p.field.is_some() || p.algebraic.is_some() || base.modulus.is_some() || base.ngens != 1
        || base.has_domain_gens() || !matches!(base.domain, Domain::ZZ | Domain::QQ)
    {
        return None;
    }
    Some((symox::poly::dense::from_mpoly(&base.rep, 0), base.domain))
}

/// Native `apart_list_full_decomposition` for `ZZ`/`QQ` numerators and the
/// factors of `Q.sqf_list_include(all=True)`: a list of `(D, numer, k)`
/// polynomials in the generator of `P`, or `None` when the inputs need the
/// generic algorithm.
#[pyfunction]
pub(super) fn _apart_full_decomposition(
    p: PyRef<'_, super::PyPoly>,
    factors: Vec<(PyRef<'_, super::PyPoly>, usize)>,
) -> PyResult<Option<Vec<(super::PyPoly, super::PyPoly, usize)>>> {
    use symox::poly::polyobj::{Domain, Poly};
    let Some((num, domain)) = plain_dense(&p.p) else { return Ok(None) };
    let gens = p.p.base.gens.clone();
    let mut dense = Vec::with_capacity(factors.len());
    for (f, n) in &factors {
        match plain_dense(&f.p) {
            Some((d, dom)) if dom == domain && f.p.base.gens == gens => dense.push((d, *n)),
            _ => return Ok(None),
        }
    }
    let terms = match symox::poly::partfrac::full_decomposition(&num, &dense, domain == Domain::ZZ) {
        Ok(t) => t,
        Err(symox::poly::convert::PolyError::Unsupported(_)) => return Ok(None),
        Err(e) => return Err(super::err(e)),
    };
    let make = |c: &[symox::poly::Rat], dom: Domain| super::PyPoly {
        p: Poly::new(symox::poly::dense::to_mpoly(c, 0, 1), gens.clone(), 1, dom).into(),
    };
    Ok(Some(
        terms
            .iter()
            .map(|t| (make(&t.denom, domain), make(&t.numer, Domain::QQ), t.exponent))
            .collect(),
    ))
}

/// Native `dispersionset` for univariate `ZZ`/`QQ` polynomials of positive
/// degree (`None` for other representations).
#[pyfunction]
#[pyo3(signature=(p, q=None))]
pub(super) fn _dispersion_set(
    p: PyRef<'_, super::PyPoly>,
    q: Option<PyRef<'_, super::PyPoly>>,
) -> PyResult<Option<Vec<super::PyExpr>>> {
    let ok = |f: &PolyValue| plain_dense(f).is_some_and(|(d, _)| d.len() > 1);
    if !ok(&p.p) || q.as_ref().is_some_and(|q| !ok(&q.p) || q.p.base.gens != p.p.base.gens) {
        return Ok(None);
    }
    let set = symox::poly::dispersion::dispersion_set(&p.p.base, q.as_ref().map(|q| &q.p.base));
    Ok(Some(set.into_iter().map(|a| w(Expr::num(a.to_number()))).collect()))
}

/// SymPy's `Poly.is_cyclotomic` for univariate `ZZ`/`QQ` polynomials
/// (`False` for every other representation and for non-integral coefficients).
#[pyfunction]
pub(super) fn _is_cyclotomic(p: PyRef<'_, super::PyPoly>) -> bool {
    let Some((asc, _)) = plain_dense(&p.p) else { return false };
    if asc.iter().any(|c| !c.is_integer()) {
        return false;
    }
    let desc: Vec<symox::poly::Int> = asc.iter().rev().map(|c| c.n.clone()).collect();
    symox::poly::structure::cyclotomic_p(&desc, false)
}
