//! Exact rational boundary for real-root analysis. No runtime SymPy calls.
use super::{division, domains, poly_obj, w, PolyValue, PyExpr, PyPoly};
use crate::to_expr;
use num_bigint::BigInt;
use pyo3::exceptions::{PyIndexError, PyNotImplementedError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use symox::expr::{Constant, Expr};
use symox::poly::int::Int;
use symox::poly::polyobj::{Domain, Poly};
use symox::poly::rootisolation::{self as iso, ComplexBox, EdgeRoot, Exclude, IsoError, Qd, RealBox, Rect};
use symox::poly::{dense, Rat};
pyo3::create_exception!(symox, RefinementFailed, PyValueError);

/// `QQ.convert` of a bound or precision, as the reference converts them:
/// a Float becomes the nearest rational under the default real field's
/// denominator limit, infinities and nan belong to no domain.
fn rational(value: &Bound<'_, PyAny>) -> PyResult<Rat> {
    let expression = to_expr(value)?;
    let text = || -> PyResult<String> { Ok(w(expression.clone()).into_pyobject(value.py())?.str()?.to_string()) };
    let illegal = [Constant::Infinity, Constant::NegInfinity, Constant::ComplexInfinity, Constant::NaN];
    if illegal.iter().any(|c| expression.is_constant(*c)) {
        return Err(domains::CoercionFailed::new_err(format!("{} is not in any domain", text()?)));
    }
    match expression.as_number() {
        Some(symox::Number::Float(f)) if !f.is_finite() => {
            Err(domains::CoercionFailed::new_err(format!("{} is not in any domain", text()?)))
        }
        Some(symox::Number::Float(_) | symox::Number::MpFloat(..)) => {
            symox::poly::domains::rationalize_real(&expression, 53, true)
                .map_err(|e| domains::CoercionFailed::new_err(e.to_string()))
        }
        Some(number) => Ok(symox::poly::convert::number_to_rat(number)),
        None => Err(domains::CoercionFailed::new_err(format!(
            "expected `Rational` object, got {}",
            text()?
        ))),
    }
}
fn optional_rational(value: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Rat>> {
    value.filter(|v| !v.is_none()).map(rational).transpose()
}
/// `Poly(f)` as the reference's top-level `count_roots` (not greedy) and
/// `refine_root` (greedy) build it: an expression needs a symbol generator.
fn unary_symbol(f: &Bound<'_, PyAny>, greedy: bool, action: &str) -> PyResult<PolyValue> {
    if let Ok(p) = f.downcast::<PyPoly>() {
        return Ok(p.get().p.clone());
    }
    let py = f.py();
    let module = py.import("symox")?;
    let options = PyDict::new(py);
    if !greedy {
        options.set_item("greedy", false)?;
    }
    let poly = match module.getattr("Poly")?.call((f,), Some(&options)) {
        Ok(poly) => poly,
        Err(e) if e.is_instance(py, &py.import("symox.polys.polyerrors")?.getattr("GeneratorsNeeded")?) => {
            let message = format!("{action} {}, not a polynomial", f.str()?);
            return Err(super::structural::polynomial_error(py, "PolynomialError", message));
        }
        Err(e) => return Err(e),
    };
    if !poly.getattr("gen")?.getattr("is_Symbol")?.is_truthy()? {
        return Err(super::structural::polynomial_error(
            py,
            "PolynomialError",
            "generator must be a Symbol".into(),
        ));
    }
    Ok(poly.downcast::<PyPoly>()?.get().p.clone())
}
/// The reference module-level `nroots` polynomial: `Poly(f, greedy=False)`
/// with a Symbol generator.
pub(crate) fn numeric_root_poly<'py>(f: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    let py = f.py();
    let options = PyDict::new(py);
    options.set_item("greedy", false)?;
    let poly = match py.import("symox")?.getattr("Poly")?.call((f,), Some(&options)) {
        Ok(poly) => poly,
        Err(e) if e.is_instance(py, &py.import("symox.polys.polyerrors")?.getattr("GeneratorsNeeded")?) => {
            let message = format!("Cannot compute numerical roots of {}, not a polynomial", f.str()?);
            return Err(super::structural::polynomial_error(py, "PolynomialError", message));
        }
        Err(e) => return Err(e),
    };
    if !poly.getattr("gen")?.getattr("is_Symbol")?.is_truthy()? {
        return Err(super::structural::polynomial_error(py, "PolynomialError", "generator must be a Symbol".into()));
    }
    Ok(poly)
}
fn unary<'py>(f: &Bound<'py, PyAny>) -> PyResult<PolyValue> {
    let gens = PyTuple::empty(f.py());
    Ok(division::prepare_unary(f, &gens, None)?.0)
}
pub(super) fn sturm_method(py: Python<'_>, p: &PolyValue, auto: bool) -> PyResult<Vec<PyPoly>> {
    Ok(sturm_values(py, p, auto)?.into_iter().map(|p| PyPoly { p }).collect())
}
/// The reference's `Poly.sturm`: exact over ZZ and QQ here, the other
/// domains follow the reference's dense field arithmetic in Python.
fn sturm_values(py: Python<'_>, p: &PolyValue, auto: bool) -> PyResult<Vec<PolyValue>> {
    let plain = p.field.is_none() && p.algebraic.is_none() && p.modulus.is_none() && !p.domain.is_rr();
    if plain && p.ngens == 1 && p.base.gens.len() == 1 {
        if p.domain == Domain::ZZ && !auto {
            return Err(domains::DomainError::new_err("Cannot compute Sturm sequence over ZZ"));
        }
        if p.rep.degree(0) > iso::MAX_DENSE_DEGREE as i64 {
            return Err(degree_limit());
        }
        let f: Vec<Rat> = dense::from_mpoly(&p.base.rep, 0).into_iter().rev().collect();
        return Ok(iso::sturm_sequence(&f)
            .into_iter()
            .map(|s| {
                let asc: Vec<Rat> = s.into_iter().rev().collect();
                let rep = dense::to_mpoly(&asc, 0, 1);
                Poly::new(rep, p.base.gens.clone(), 1, Domain::QQ).into()
            })
            .collect());
    }
    let values = py
        .import("symox.polys._root_counting")?
        .getattr("sturm")?
        .call1((PyPoly { p: p.clone() }, auto))?;
    values
        .try_iter()?
        .map(|v| Ok(v?.downcast::<PyPoly>()?.get().p.clone()))
        .collect()
}
fn iso_error(error: IsoError) -> PyErr {
    match error {
        IsoError::Domain(message) => domains::DomainError::new_err(message),
        IsoError::Refinement(message) => RefinementFailed::new_err(message),
        IsoError::Value(message) => PyValueError::new_err(message),
        IsoError::NotImplemented(message) => PyNotImplementedError::new_err(message),
        IsoError::ZeroDivision => PyZeroDivisionError::new_err("division by zero"),
        IsoError::Index => PyIndexError::new_err("list index out of range"),
    }
}
fn number(r: &Rat) -> PyExpr {
    w(Expr::num(r.to_number()))
}
fn gaussian(re: &Rat, im: &Rat) -> PyExpr {
    w(symox::ops::add2(
        &Expr::num(re.to_number()),
        &symox::ops::mul2(&Expr::i(), &Expr::num(im.to_number())),
    ))
}
/// `QQ.convert(eps)` with the reference's positivity check.
fn positive_eps(eps: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Rat>> {
    let eps = optional_rational(eps)?;
    if eps.as_ref().is_some_and(|e| !e.n.is_positive()) {
        return Err(PyValueError::new_err("'eps' must be a positive rational"));
    }
    Ok(eps)
}
/// The dense coefficients, leading first, of a univariate polynomial over
/// ZZ or QQ; `action` completes the reference's domain error message.
fn dense_rational(py: Python<'_>, p: &PolyValue, multivariate: &str, action: &str) -> PyResult<Vec<Rat>> {
    if p.ngens != 1 {
        return Err(super::structural::polynomial_error(py, "PolynomialError", multivariate.into()));
    }
    let composite = p.base.gens.len() != 1;
    if composite || p.field.is_some() || p.algebraic.is_some() || p.modulus.is_some() || p.domain.is_rr() {
        let domain = match p.modulus {
            Some(prime) if p.base.gens.len() == p.ngens => format!("GF({prime})"),
            _ => p.domain_str(),
        };
        return Err(domains::DomainError::new_err(format!("{action} {domain}")));
    }
    if p.rep.degree(0) > iso::MAX_DENSE_DEGREE as i64 {
        return Err(degree_limit());
    }
    Ok(dense::from_mpoly(&p.base.rep, 0).into_iter().rev().collect())
}
/// The reference's `Poly.count_roots`: real bounds count distinct real
/// roots in `[inf, sup]` with Sturm sequences, a complex bound counts roots
/// in the closed rectangle with those corners.
pub(super) fn count_method(
    py: Python<'_>,
    p: &PolyValue,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    // `(re, im, real)` of a bound; `None` for the matching infinity.
    let parts = |value: Option<&Bound<'_, PyAny>>, infinity: Constant| -> PyResult<Option<(Rat, Rat, bool)>> {
        let Some(value) = value.filter(|v| !v.is_none()) else { return Ok(None) };
        let expression = to_expr(value)?;
        if expression.is_constant(infinity) {
            return Ok(None);
        }
        let object = w(expression).into_pyobject(py)?;
        let pair = object.call_method0("as_real_imag")?;
        let (re, im): (Bound<'_, PyAny>, Bound<'_, PyAny>) = pair.extract()?;
        if !im.is_truthy()? {
            return Ok(Some((rational(&object)?, Rat::zero(), true)));
        }
        Ok(Some((rational(&re)?, rational(&im)?, false)))
    };
    let lower = parts(inf, Constant::NegInfinity)?;
    let upper = parts(sup, Constant::Infinity)?;
    if lower.as_ref().is_none_or(|b| b.2) && upper.as_ref().is_none_or(|b| b.2) {
        let lower = lower.map(|b| b.0);
        let upper = upper.map(|b| b.0);
        return count_real(py, p, lower, upper);
    }
    let f = dense_rational(
        py,
        p,
        "Cannot count roots of a multivariate polynomial",
        "complex root counting is not supported over",
    )?;
    let lower = lower.map(|(a, b, _)| (a, b));
    let upper = upper.map(|(a, b, _)| (a, b));
    let count = iso::count_complex_roots(&f, lower, upper, &Exclude::None).map_err(iso_error)?;
    Ok(w(Expr::int(count)))
}
/// `dup_count_real_roots`: exact over ZZ and QQ; other coefficient domains
/// follow the reference's generic Sturm sequence over their field.
fn count_real(py: Python<'_>, p: &PolyValue, inf: Option<Rat>, sup: Option<Rat>) -> PyResult<PyExpr> {
    let plain = p.field.is_none() && p.algebraic.is_none() && p.modulus.is_none() && !p.domain.is_rr();
    if plain && p.ngens == 1 && p.base.gens.len() == 1 {
        if p.rep.degree(0) > iso::MAX_DENSE_DEGREE as i64 {
            return Err(degree_limit());
        }
        let f: Vec<Rat> = dense::from_mpoly(&p.base.rep, 0).into_iter().rev().collect();
        return Ok(w(Expr::int(iso::count_real_roots(&f, inf.as_ref(), sup.as_ref()) as i64)));
    }
    if plain && p.ngens > 1 && p.base.gens.len() == p.ngens && p.rep.degree(0) > 0 {
        // The reference converts the multivariate representation as if it
        // were univariate and fails with these interpreter errors.
        let message = if p.domain == Domain::ZZ {
            "PythonMPQ() requires numeric or string argument"
        } else {
            "'<' not supported between instances of 'list' and 'int'"
        };
        return Err(pyo3::exceptions::PyTypeError::new_err(message));
    }
    let bound = |b: Option<Rat>| -> Option<PyExpr> { b.map(|r| number(&r)) };
    let count = py
        .import("symox.polys._root_counting")?
        .getattr("count_real_roots")?
        .call1((PyPoly { p: p.clone() }, bound(inf), bound(sup)))?;
    Ok(w(to_expr(&count)?))
}
fn degree_limit() -> PyErr {
    PyValueError::new_err(format!(
        "real-root resource limit: degree exceeds {}",
        iso::MAX_DENSE_DEGREE
    ))
}
fn interval_tuple(py: Python<'_>, s: &Rat, t: &Rat) -> PyResult<PyObject> {
    Ok(PyTuple::new(py, [number(s), number(t)])?.into_any().unbind())
}
fn rectangle_tuple(py: Python<'_>, b: &ComplexBox) -> PyResult<PyObject> {
    Ok(PyTuple::new(py, [gaussian(&b.ax(), &b.ay()), gaussian(&b.bx(), &b.by())])?.into_any().unbind())
}
fn with_multiplicity(py: Python<'_>, value: PyObject, k: u32) -> PyResult<PyObject> {
    Ok(PyTuple::new(py, [value, k.into_pyobject(py)?.into_any().unbind()])?.into_any().unbind())
}
#[allow(clippy::too_many_arguments)]
pub(super) fn intervals_method(
    py: Python<'_>,
    p: &PolyValue,
    all: bool,
    eps: Option<&Bound<'_, PyAny>>,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
    fast: bool,
    sqf: bool,
) -> PyResult<PyObject> {
    let eps = positive_eps(eps)?;
    let inf = optional_rational(inf)?;
    let sup = optional_rational(sup)?;
    let action = if all {
        "isolation of real and complex roots is not supported over"
    } else {
        "isolation of real roots not supported over"
    };
    let f = dense_rational(py, p, "Cannot isolate roots of a multivariate polynomial", action)?;
    let (eps, inf, sup) = (eps.as_ref(), inf.as_ref(), sup.as_ref());
    if !all {
        let result = PyList::empty(py);
        let g = iso::clear_denoms(&f);
        if sqf {
            for (s, t) in iso::isolate_real_roots_sqf(&g, eps, inf, sup, fast).map_err(iso_error)?.0 {
                result.append(interval_tuple(py, &s, &t)?)?;
            }
        } else {
            for ((s, t), k) in iso::isolate_real_roots(&g, eps, inf, sup, fast).map_err(iso_error)? {
                result.append(with_multiplicity(py, interval_tuple(py, &s, &t)?, k)?)?;
            }
        }
        return Ok(result.into_any().unbind());
    }
    let (real, complex) = (PyList::empty(py), PyList::empty(py));
    if sqf {
        let (r, c) = iso::isolate_all_roots_sqf(&f, eps, inf, sup, fast).map_err(iso_error)?;
        for (s, t) in r {
            real.append(interval_tuple(py, &s, &t)?)?;
        }
        for b in c {
            complex.append(rectangle_tuple(py, &b)?)?;
        }
    } else {
        let (r, c) = iso::isolate_all_roots(&f, eps, inf, sup, fast).map_err(iso_error)?;
        for ((s, t), k) in r {
            real.append(with_multiplicity(py, interval_tuple(py, &s, &t)?, k)?)?;
        }
        for (b, k) in c {
            complex.append(with_multiplicity(py, rectangle_tuple(py, &b)?, k)?)?;
        }
    }
    Ok(PyTuple::new(py, [real, complex])?.into_any().unbind())
}
pub(super) fn refine_method(
    p: &PolyValue,
    s: &Bound<'_, PyAny>,
    t: &Bound<'_, PyAny>,
    eps: Option<&Bound<'_, PyAny>>,
    steps: Option<i64>,
    fast: bool,
    check_sqf: bool,
) -> PyResult<(PyExpr, PyExpr)> {
    let py = s.py();
    if check_sqf && !p.is_sqf_checked()? {
        return Err(super::structural::polynomial_error(
            py,
            "PolynomialError",
            "only square-free polynomials supported".into(),
        ));
    }
    let (s, t) = (rational(s)?, rational(t)?);
    let eps = positive_eps(eps)?;
    let steps = match steps {
        Some(n) => Some(n.max(0) as u64),
        None if eps.is_none() => Some(1),
        None => None,
    };
    let f = dense_rational(
        py,
        p,
        "Cannot refine a root of a multivariate polynomial",
        "real root refinement not supported over",
    )?;
    let (a, b) = iso::refine_real_root(&iso::clear_denoms(&f), &s, &t, eps.as_ref(), steps, None, fast)
        .map_err(iso_error)?;
    Ok((number(&a), number(&b)))
}
trait SquareFree {
    fn is_sqf_checked(&self) -> PyResult<bool>;
}
impl SquareFree for PolyValue {
    fn is_sqf_checked(&self) -> PyResult<bool> {
        if self.ngens != 1 || self.field.is_some() || self.algebraic.is_some() || self.modulus.is_some() {
            return Ok(true);
        }
        let f: Vec<Rat> = dense::from_mpoly(&self.base.rep, 0).into_iter().rev().collect();
        Ok(iso::is_squarefree(&f))
    }
}

fn big_rat(value: &(BigInt, BigInt)) -> PyResult<Rat> {
    if value.1 == BigInt::from(0) {
        return Err(PyZeroDivisionError::new_err("division by zero"));
    }
    Ok(Rat::new(Int::from_big(value.0.clone()), Int::from_big(value.1.clone())))
}
fn rat_big(r: &Rat) -> (BigInt, BigInt) {
    (r.n.to_big(), r.d.to_big())
}
fn ints(values: &[BigInt]) -> Vec<Int> {
    values.iter().map(|c| Int::from_big(c.clone())).collect()
}

/// The refinement state of a real isolating interval: the polynomial over
/// ZZ transformed by a Mobius map, as in the reference's `RealInterval`.
#[pyclass(name = "RealIsolation", module = "symox._core", frozen)]
pub struct PyRealIsolation {
    state: RealBox,
}
#[pymethods]
impl PyRealIsolation {
    /// From an isolating interval `(s, t)` of a root of `f` over ZZ.
    #[staticmethod]
    fn from_interval(f: Vec<BigInt>, s: (BigInt, BigInt), t: (BigInt, BigInt)) -> PyResult<Self> {
        let state = RealBox::from_interval(&ints(&f), &big_rat(&s)?, &big_rat(&t)?).map_err(iso_error)?;
        Ok(Self { state })
    }
    /// From the stored Mobius transform and the transformed polynomial.
    #[staticmethod]
    fn from_mobius(f: Vec<BigInt>, mobius: (BigInt, BigInt, BigInt, BigInt), neg: bool) -> Self {
        let m = (
            Int::from_big(mobius.0),
            Int::from_big(mobius.1),
            Int::from_big(mobius.2),
            Int::from_big(mobius.3),
        );
        Self { state: RealBox { f: ints(&f), m, neg } }
    }
    fn a(&self) -> PyResult<(BigInt, BigInt)> {
        Ok(rat_big(&self.state.a().map_err(iso_error)?))
    }
    fn b(&self) -> PyResult<(BigInt, BigInt)> {
        Ok(rat_big(&self.state.b().map_err(iso_error)?))
    }
    #[getter]
    fn neg(&self) -> bool {
        self.state.neg
    }
    #[getter]
    fn mobius(&self) -> (BigInt, BigInt, BigInt, BigInt) {
        let (a, b, c, d) = &self.state.m;
        (a.to_big(), b.to_big(), c.to_big(), d.to_big())
    }
    #[getter]
    fn f(&self) -> Vec<BigInt> {
        self.state.f.iter().map(Int::to_big).collect()
    }
    fn refine(&self) -> PyResult<Self> {
        Ok(Self { state: self.state.refine().map_err(iso_error)? })
    }
    fn refine_size(&self, dx: (BigInt, BigInt)) -> PyResult<Self> {
        Ok(Self { state: self.state.refine_size(&big_rat(&dx)?).map_err(iso_error)? })
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other.downcast::<Self>().is_ok_and(|o| o.get().state == self.state)
    }
}

/// The refinement state of a complex isolating rectangle, as in the
/// reference's `ComplexInterval`.
#[pyclass(name = "ComplexIsolation", module = "symox._core", frozen)]
pub struct PyComplexIsolation {
    state: ComplexBox,
}
#[pymethods]
impl PyComplexIsolation {
    fn ax(&self) -> (BigInt, BigInt) {
        rat_big(&self.state.ax())
    }
    fn ay(&self) -> (BigInt, BigInt) {
        rat_big(&self.state.ay())
    }
    fn bx(&self) -> (BigInt, BigInt) {
        rat_big(&self.state.bx())
    }
    fn by(&self) -> (BigInt, BigInt) {
        rat_big(&self.state.by())
    }
    /// The stored corners `((u, v), (s, t))` in the upper half plane.
    fn corners(&self) -> ((BigInt, BigInt), (BigInt, BigInt), (BigInt, BigInt), (BigInt, BigInt)) {
        let r = &self.state.rect;
        (rat_big(&r.a.0), rat_big(&r.a.1), rat_big(&r.b.0), rat_big(&r.b.1))
    }
    #[getter]
    fn conj(&self) -> bool {
        self.state.conj
    }
    /// The quadrant sequences of the four edges.
    fn quadrants(&self) -> Vec<Vec<&'static str>> {
        self.state.rect.q.iter().map(|q| q.iter().map(|x| x.name()).collect()).collect()
    }
    /// The polynomial over QQ whose root is isolated, leading coefficient first.
    fn poly(&self) -> Vec<(BigInt, BigInt)> {
        self.state.f.iter().map(rat_big).collect()
    }
    /// The isolated roots on the four edges: `((s, t), {index: k}, h)`.
    fn edges<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyList>>> {
        let mut out = Vec::new();
        for edge in &self.state.rect.i {
            let list = PyList::empty(py);
            for r in edge {
                let indices = PyDict::new(py);
                for (k, m) in &r.indices {
                    indices.set_item(*k, *m)?;
                }
                let h: Vec<BigInt> = r.h.iter().map(Int::to_big).collect();
                list.append(((rat_big(&r.s), rat_big(&r.t)), indices, h))?;
            }
            out.push(list);
        }
        Ok(out)
    }
    /// The restrictions of the real and imaginary parts to the four edges.
    fn edge_polys(&self) -> (Vec<Vec<(BigInt, BigInt)>>, Vec<Vec<(BigInt, BigInt)>>) {
        let dense = |fs: &[Vec<Rat>; 4]| fs.iter().map(|f| f.iter().map(rat_big).collect()).collect();
        (dense(&self.state.rect.f1), dense(&self.state.rect.f2))
    }
    /// The real and imaginary parts of `f(x + i*y)`, dense in `x` over
    /// dense polynomials in `y`.
    #[allow(clippy::type_complexity)]
    fn real_imag(&self) -> (Vec<Vec<(BigInt, BigInt)>>, Vec<Vec<(BigInt, BigInt)>>) {
        let (f1, f2) = iso::real_imag(&self.state.f);
        let dense = |fs: Vec<Vec<Rat>>| fs.iter().map(|f| f.iter().map(rat_big).collect()).collect();
        (dense(f1), dense(f2))
    }
    /// Rebuild the state from the reference's constructor arguments.
    #[staticmethod]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn from_parts(
        a: ((BigInt, BigInt), (BigInt, BigInt)),
        b: ((BigInt, BigInt), (BigInt, BigInt)),
        edges: Vec<Vec<(((BigInt, BigInt), (BigInt, BigInt)), Bound<'_, PyDict>, Vec<BigInt>)>>,
        quadrants: Vec<Vec<String>>,
        f1: Vec<Vec<(BigInt, BigInt)>>,
        f2: Vec<Vec<(BigInt, BigInt)>>,
        f: Vec<(BigInt, BigInt)>,
        conj: bool,
    ) -> PyResult<Self> {
        let four = |what: &str, n: usize| -> PyResult<()> {
            if n == 4 {
                Ok(())
            } else {
                Err(PyValueError::new_err(format!("expected four {what}, one per edge")))
            }
        };
        four("lists of edge roots", edges.len())?;
        four("quadrant sequences", quadrants.len())?;
        four("real parts", f1.len())?;
        four("imaginary parts", f2.len())?;
        let mut i: [Vec<EdgeRoot>; 4] = Default::default();
        for (k, edge) in edges.into_iter().enumerate() {
            for ((s, t), indices, h) in edge {
                let mut pairs = Vec::new();
                for (key, value) in indices.iter() {
                    pairs.push((key.extract::<usize>()?, value.extract::<u32>()?));
                }
                i[k].push(EdgeRoot { s: big_rat(&s)?, t: big_rat(&t)?, indices: pairs, h: ints(&h) });
            }
        }
        let mut q: [Vec<Qd>; 4] = Default::default();
        for (k, names) in quadrants.iter().enumerate() {
            for name in names {
                q[k].push(Qd::from_name(name).ok_or_else(|| PyValueError::new_err(format!("unknown quadrant {name}")))?);
            }
        }
        let rats = |f: &Vec<(BigInt, BigInt)>| f.iter().map(big_rat).collect::<PyResult<Vec<Rat>>>();
        let mut g1: [Vec<Rat>; 4] = Default::default();
        let mut g2: [Vec<Rat>; 4] = Default::default();
        for k in 0..4 {
            g1[k] = rats(&f1[k])?;
            g2[k] = rats(&f2[k])?;
        }
        let rect = Rect {
            n: 1,
            a: (big_rat(&a.0)?, big_rat(&a.1)?),
            b: (big_rat(&b.0)?, big_rat(&b.1)?),
            i,
            q,
            f1: g1,
            f2: g2,
        };
        Ok(Self { state: ComplexBox { rect, f: std::sync::Arc::new(rats(&f)?), conj } })
    }
    fn conjugate(&self) -> Self {
        Self { state: self.state.conjugate() }
    }
    fn refine(&self) -> PyResult<Self> {
        Ok(Self { state: self.state.refine().map_err(iso_error)? })
    }
    fn refine_size(&self, dx: (BigInt, BigInt), dy: (BigInt, BigInt)) -> PyResult<Self> {
        Ok(Self { state: self.state.refine_size(&big_rat(&dx)?, &big_rat(&dy)?).map_err(iso_error)? })
    }
    fn is_disjoint(&self, other: PyRef<'_, Self>) -> bool {
        self.state.is_disjoint(&other.state)
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other.downcast::<Self>().is_ok_and(|o| o.get().state == self.state)
    }
}

/// Blackbox real isolating intervals of a square-free `f` over ZZ.
#[pyfunction]
#[pyo3(signature=(f, eps=None, inf=None, sup=None, fast=false))]
fn isolate_real_roots_sqf_state(
    f: Vec<BigInt>,
    eps: Option<(BigInt, BigInt)>,
    inf: Option<(BigInt, BigInt)>,
    sup: Option<(BigInt, BigInt)>,
    fast: bool,
) -> PyResult<Vec<PyRealIsolation>> {
    let eps = eps.as_ref().map(big_rat).transpose()?;
    let inf = inf.as_ref().map(big_rat).transpose()?;
    let sup = sup.as_ref().map(big_rat).transpose()?;
    let (intervals, g) =
        iso::isolate_real_roots_sqf(&ints(&f), eps.as_ref(), inf.as_ref(), sup.as_ref(), fast).map_err(iso_error)?;
    intervals
        .iter()
        .map(|(s, t)| RealBox::from_interval(&g, s, t).map(|state| PyRealIsolation { state }).map_err(iso_error))
        .collect()
}
/// Blackbox complex isolating rectangles of a square-free `f` over QQ.
#[pyfunction]
#[pyo3(signature=(f, eps=None, inf=None, sup=None))]
fn isolate_complex_roots_sqf_state(
    f: Vec<(BigInt, BigInt)>,
    eps: Option<(BigInt, BigInt)>,
    inf: Option<(BigInt, BigInt)>,
    sup: Option<(BigInt, BigInt)>,
) -> PyResult<Vec<PyComplexIsolation>> {
    let f = f.iter().map(big_rat).collect::<PyResult<Vec<_>>>()?;
    let eps = eps.as_ref().map(big_rat).transpose()?;
    let inf = inf.as_ref().map(big_rat).transpose()?;
    let sup = sup.as_ref().map(big_rat).transpose()?;
    Ok(iso::isolate_complex_roots_sqf(&f, eps.as_ref(), inf.as_ref(), sup.as_ref())
        .map_err(iso_error)?
        .into_iter()
        .map(|state| PyComplexIsolation { state })
        .collect())
}
/// Number of roots of `f` over QQ in a closed rectangle, with the
/// reference's edge and corner exclusions.
#[pyfunction]
#[pyo3(signature=(f, inf=None, sup=None, exclude=None))]
fn count_complex_roots_state(
    f: Vec<(BigInt, BigInt)>,
    inf: Option<((BigInt, BigInt), (BigInt, BigInt))>,
    sup: Option<((BigInt, BigInt), (BigInt, BigInt))>,
    exclude: Option<&Bound<'_, PyAny>>,
) -> PyResult<i64> {
    let f = f.iter().map(big_rat).collect::<PyResult<Vec<_>>>()?;
    let point = |p: Option<((BigInt, BigInt), (BigInt, BigInt))>| -> PyResult<Option<(Rat, Rat)>> {
        p.map(|(a, b)| Ok((big_rat(&a)?, big_rat(&b)?))).transpose()
    };
    let exclude = match exclude {
        None => Exclude::None,
        Some(v) if v.is_none() => Exclude::None,
        Some(v) if v.is_instance_of::<pyo3::types::PyBool>() && v.extract::<bool>()? => Exclude::All,
        Some(v) if v.is_instance_of::<pyo3::types::PyBool>() => Exclude::None,
        Some(v) => {
            let names: Vec<String> = v.try_iter()?.map(|x| x?.extract::<String>()).collect::<PyResult<_>>()?;
            Exclude::from_names(&names)
        }
    };
    iso::count_complex_roots(&f, point(inf)?, point(sup)?, &exclude).map_err(iso_error)
}
/// The reference's `dup_isolate_real_roots_list`, returning
/// `((s, t), {index: multiplicity})` pairs.
fn isolate_list(
    py: Python<'_>,
    polys: &[Vec<Rat>],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    strict: bool,
    fast: bool,
) -> PyResult<PyObject> {
    let polys: Vec<_> = polys.iter().map(|f| iso::clear_denoms(f)).collect();
    let result = PyList::empty(py);
    for root in iso::isolate_real_roots_list(&polys, eps, inf, sup, strict, fast).map_err(iso_error)? {
        let indices = PyDict::new(py);
        for (k, m) in &root.indices {
            indices.set_item(w(Expr::int(*k as i64)), w(Expr::int(*m as i64)))?;
        }
        result.append(PyTuple::new(py, [interval_tuple(py, &root.s, &root.t)?, indices.into_any().unbind()])?)?;
    }
    Ok(result.into_any().unbind())
}
#[pyfunction]
#[pyo3(signature=(f,*gens,**kwargs))]
fn sturm<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyObject> {
    let options = PyDict::new(f.py());
    let mut auto = true;
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            if key.extract::<String>()? == "auto" {
                auto = division::boolean(&value, "auto")?.unwrap_or(true);
            } else {
                options.set_item(key, value)?;
            }
        }
    }
    let (p, polys) = division::prepare_unary(f, gens, Some(&options))?;
    let values = sturm_values(f.py(), &p, auto)?;
    Ok(PyList::new(
        f.py(),
        values
            .into_iter()
            .map(|p| {
                if polys {
                    poly_obj(f.py(), p)
                } else {
                    Ok(w(p.as_expr()).into_pyobject(f.py())?.into_any().unbind())
                }
            })
            .collect::<PyResult<Vec<_>>>()?,
    )?
    .into_any()
    .unbind())
}
#[pyfunction]
#[pyo3(signature=(f,inf=None,sup=None))]
fn count_roots(
    f: &Bound<'_, PyAny>,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    count_method(f.py(), &unary_symbol(f, false, "Cannot count roots of")?, inf, sup)
}
#[pyfunction]
#[pyo3(signature=(F,all=false,eps=None,inf=None,sup=None,strict=false,fast=false,sqf=false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn intervals(
    F: &Bound<'_, PyAny>,
    all: bool,
    eps: Option<&Bound<'_, PyAny>>,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
    strict: bool,
    fast: bool,
    sqf: bool,
) -> PyResult<PyObject> {
    let py = F.py();
    if F.downcast::<PyPoly>().is_err() && !F.is_instance_of::<pyo3::types::PyString>() && F.hasattr("__iter__")? {
        let items: Vec<Bound<'_, PyAny>> = F.try_iter()?.collect::<PyResult<_>>()?;
        let module = py.import("symox")?;
        let options = PyDict::new(py);
        options.set_item("domain", "QQ")?;
        let parallel = module.getattr("parallel_poly_from_expr")?.call((PyList::new(py, items)?,), Some(&options))?;
        let (polys, opt): (Vec<Bound<'_, PyAny>>, Bound<'_, PyAny>) = parallel.extract()?;
        if opt.getattr("gens")?.len()? > 1 {
            let error = module.getattr("MultivariatePolynomialError")?;
            return Err(PyErr::from_value(error.call0()?));
        }
        let mut dense_polys = Vec::new();
        for p in &polys {
            let p = p.downcast::<PyPoly>()?.get().p.clone();
            dense_polys.push(dense::from_mpoly(&p.base.rep, 0).into_iter().rev().collect::<Vec<Rat>>());
        }
        let eps = positive_eps(eps)?;
        let inf = optional_rational(inf)?;
        let sup = optional_rational(sup)?;
        return isolate_list(py, &dense_polys, eps.as_ref(), inf.as_ref(), sup.as_ref(), strict, fast);
    }
    let p = if F.downcast::<PyPoly>().is_err() {
        let expression = to_expr(F)?;
        if symox::traverse::free_symbols(&expression).is_empty() {
            return Ok(PyList::empty(py).into_any().unbind());
        }
        unary(F)?
    } else {
        unary(F)?
    };
    intervals_method(py, &p, all, eps, inf, sup, fast, sqf)
}
#[pyfunction]
#[pyo3(signature=(f,s,t,eps=None,steps=None,fast=false,check_sqf=false))]
fn refine_root(
    f: &Bound<'_, PyAny>,
    s: &Bound<'_, PyAny>,
    t: &Bound<'_, PyAny>,
    eps: Option<&Bound<'_, PyAny>>,
    steps: Option<i64>,
    fast: bool,
    check_sqf: bool,
) -> PyResult<(PyExpr, PyExpr)> {
    refine_method(&unary_symbol(f, true, "Cannot refine a root of")?, s, t, eps, steps, fast, check_sqf)
}
/// The root object of index `index` of `expr`, an irreducible polynomial
/// over ZZ in one symbol. `kind`, when known to the caller from its
/// isolating intervals, is `(real, negative, imaginary)`.
#[pyfunction]
#[pyo3(signature=(expr, index, kind=None))]
fn complex_root_node(expr: &Bound<'_, PyAny>, index: usize, kind: Option<(bool, bool, bool)>) -> PyResult<PyExpr> {
    use symox::poly::rootof;
    let expression = to_expr(expr)?;
    let root = match kind {
        Some((real, negative, imaginary)) => {
            let kind = rootof::RootKind { real, negative, imaginary, lower: false, nreal: 0 };
            rootof::node(expression, index, &rootof::root_facts(&kind))
        }
        None => rootof::complex_root(&expression, index),
    };
    root.map(w).ok_or_else(|| PyValueError::new_err("not an indexed root of an irreducible integer polynomial"))
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("RefinementFailed", m.py().get_type::<RefinementFailed>())?;
    m.add_function(wrap_pyfunction!(sturm, m)?)?;
    m.add_function(wrap_pyfunction!(count_roots, m)?)?;
    m.add_function(wrap_pyfunction!(intervals, m)?)?;
    m.add_function(wrap_pyfunction!(refine_root, m)?)?;
    m.add_class::<PyRealIsolation>()?;
    m.add_class::<PyComplexIsolation>()?;
    m.add_function(wrap_pyfunction!(isolate_real_roots_sqf_state, m)?)?;
    m.add_function(wrap_pyfunction!(isolate_complex_roots_sqf_state, m)?)?;
    m.add_function(wrap_pyfunction!(count_complex_roots_state, m)?)?;
    m.add_function(wrap_pyfunction!(complex_root_node, m)?)?;
    Ok(())
}
