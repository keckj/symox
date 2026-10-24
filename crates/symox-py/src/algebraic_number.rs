//! Exact field-element wrappers. Expressions enter the existing scalar engine
//! through `_symox_`; their defining field and power-basis coordinates remain
//! available on this immutable wrapper.
#![allow(non_snake_case)]
use super::{err, minimal_polynomial_value, rat_expr, w, PyPoly};
use crate::{to_expr, PyExpr};
use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySet, PyTuple};
use symox::expr::Expr;
use symox::ops;
use symox::poly::algebraic::{NumberField, AN};
use symox::poly::algebraic_model::PrimitiveModel;
use symox::poly::int::{Coeff, Rat};
use symox::poly::polyobj::{Domain, Poly};

pyo3::create_exception!(symox, IsomorphismFailed, PyValueError);

#[pyclass(name = "AlgebraicNumber", module = "symox._core", frozen)]
#[derive(Clone)]
pub struct PyAlgebraicNumber {
    value: AN,
    alias: Option<Expr>,
    gen: Expr,
}
fn variable(x: Option<&Bound<'_, PyAny>>) -> PyResult<Expr> {
    let x = x
        .map(to_expr)
        .transpose()?
        .unwrap_or_else(|| Expr::symbol("_x"));
    if !x.is_symbol() {
        return Err(PyValueError::new_err(
            "the polynomial generator must be a symbol",
        ));
    }
    Ok(x)
}
fn alias_expr(alias: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Expr>> {
    alias
        .map(|a| {
            if a.is_none() {
                return Ok(None);
            }
            let e = if let Ok(name) = a.extract::<String>() {
                Expr::symbol(&name)
            } else {
                to_expr(a)?
            };
            if !e.is_symbol() {
                return Err(PyValueError::new_err("alias must be a symbol or a string"));
            }
            Ok(Some(e))
        })
        .transpose()
        .map(Option::flatten)
}
fn coefficients(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Rat>> {
    let mut out = obj
        .try_iter()?
        .map(|e| {
            let e = to_expr(&e?)?;
            e.as_number().and_then(Rat::from_number).ok_or_else(|| {
                PyTypeError::new_err("algebraic-number coefficients must be exact rational numbers")
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    out.reverse();
    Ok(out)
}
fn polynomial_value(value: &AN, coefficients: &[Rat]) -> AN {
    let mut out = AN::rational(Rat::zero(), value.field.clone());
    for c in coefficients.iter().rev() {
        out = out
            .mul(value)
            .add(&AN::rational(c.clone(), value.field.clone()));
    }
    out
}
fn polynomial_expression(coefficients: &[Rat], generator: &Expr) -> Expr {
    ops::add(
        &coefficients
            .iter()
            .enumerate()
            .map(|(i, c)| ops::mul(&[rat_expr(c), generator.powi(i as i64)]))
            .collect::<Vec<_>>(),
    )
}
fn polynomial(coefficients: &[Rat], gen: Expr) -> Poly {
    Poly {
        rep: symox::poly::dense::to_mpoly(coefficients, 0, 1),
        gens: vec![gen],
        ngens: 1,
        domain: Domain::QQ,
        modulus: None,
    }
}
fn supplied_minimal(obj: &Bound<'_, PyAny>, root: &Expr) -> PyResult<()> {
    let p = if let Ok(p) = obj.downcast::<PyPoly>() {
        let p = p.get();
        if p.p.true_gens().len() != 1
            || p.p.algebraic.is_some()
            || p.p.field.is_some()
            || p.p.modulus.is_some()
        {
            return Err(PyValueError::new_err(
                "a defining polynomial must be univariate over QQ",
            ));
        }
        Poly::from_expr(&p.p.as_expr(), p.p.true_gens()).map_err(err)?
    } else {
        let expr = to_expr(obj)?;
        let gens = symox::traverse::free_symbols(&expr)
            .into_iter()
            .collect::<Vec<_>>();
        if gens.len() != 1 {
            return Err(PyValueError::new_err(
                "a defining polynomial must have one generator",
            ));
        }
        Poly::from_expr(&expr, &gens).map_err(err)?
    };
    if p.domain.is_rr() {
        return Err(PyValueError::new_err(
            "the defining polynomial must be exact",
        ));
    }
    let mut coefficients = symox::poly::dense::from_mpoly(&p.rep, 0);
    let Some(leading) = coefficients.last().cloned() else {
        return Err(PyValueError::new_err(
            "a defining polynomial cannot be zero",
        ));
    };
    for c in &mut coefficients {
        *c = c.div(&leading);
    }
    let actual = PrimitiveModel::build(root).map_err(err)?;
    if coefficients != actual.minimal {
        return Err(PyValueError::new_err(
            "the supplied polynomial is not the irreducible minimal polynomial of this root",
        ));
    }
    Ok(())
}
impl PyAlgebraicNumber {
    fn build(
        expr: &Bound<'_, PyAny>,
        coeffs: Option<&Bound<'_, PyAny>>,
        alias: Option<&Bound<'_, PyAny>>,
        gen: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let previous = expr.downcast::<PyAlgebraicNumber>().ok();
        let alias = if alias.is_none() {
            previous.and_then(|p| p.get().alias.clone())
        } else {
            alias_expr(alias)?
        };
        let gen = if gen.is_none() {
            previous
                .map(|p| p.get().gen.clone())
                .unwrap_or_else(|| Expr::symbol("_x"))
        } else {
            variable(gen)?
        };
        let mut value = if let Some(previous) = previous {
            previous.get().value.clone()
        } else {
            let root = if let Ok(pair) = expr.downcast::<PyTuple>() {
                if pair.len() != 2 {
                    return Err(PyValueError::new_err(
                        "a defining pair must contain (minimal_polynomial, root)",
                    ));
                }
                let root = to_expr(&pair.get_item(1)?)?;
                supplied_minimal(&pair.get_item(0)?, &root)?;
                root
            } else {
                to_expr(expr)?
            };
            let field = NumberField::for_element(&root).map_err(err)?;
            AN::alpha(field)
        };
        if let Some(coeffs) = coeffs.filter(|c| !c.is_none()) {
            value = polynomial_value(&value, &coefficients(coeffs)?);
        }
        Ok(Self { value, alias, gen })
    }
    fn coords(&self) -> Vec<Rat> {
        let mut c = self.value.coefficients.clone();
        if c.is_empty() {
            c.push(Rat::zero());
        }
        c.reverse();
        c
    }
    fn minpoly_value(&self) -> PyResult<Poly> {
        let model = PrimitiveModel::build(&self.value.field.generator).map_err(err)?;
        Ok(minimal_polynomial_value(&model, self.gen.clone()))
    }
    fn actual_minpoly(&self) -> PyResult<Poly> {
        let model = PrimitiveModel::build(&self.value.as_expr()).map_err(err)?;
        Ok(minimal_polynomial_value(&model, self.gen.clone()))
    }
    /// The reference `as_expr`: the coordinate polynomial at `x` (the
    /// root by default), expanded.
    fn represented(&self, x: Option<&Bound<'_, PyAny>>) -> PyResult<Expr> {
        let x = x
            .map(to_expr)
            .transpose()?
            .unwrap_or_else(|| self.value.field.generator.clone());
        Ok(symox::expand::expand(&polynomial_expression(
            &self.value.coefficients,
            &x,
        )))
    }
    /// A minimal polynomial with the reference class: `PurePoly` unless a
    /// generator was given.
    fn minimal_object(&self, py: Python<'_>, p: Poly) -> PyResult<PyObject> {
        let poly = PyPoly { p: p.into() };
        if self.gen != Expr::symbol("_x") {
            return Ok(Py::new(py, poly)?.into_any());
        }
        Ok(py.import("symox")?.getattr("PurePoly")?.call1((poly,))?.unbind())
    }
}
#[pymethods]
impl PyAlgebraicNumber {
    #[new]
    #[pyo3(signature=(expr,coeffs=None,alias=None,*,gen=None))]
    fn new(
        expr: &Bound<'_, PyAny>,
        coeffs: Option<&Bound<'_, PyAny>>,
        alias: Option<&Bound<'_, PyAny>>,
        gen: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(expr, coeffs, alias, gen)
    }
    #[getter]
    fn root(&self) -> PyExpr {
        w(self.value.field.generator.clone())
    }
    #[getter]
    fn alias(&self) -> Option<PyExpr> {
        self.alias.clone().map(w)
    }
    #[getter]
    fn minpoly(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.minimal_object(py, self.minpoly_value()?)
    }
    fn coeffs(&self) -> Vec<PyExpr> {
        self.coords().iter().map(|c| w(rat_expr(c))).collect()
    }
    fn native_coeffs(&self, py: Python<'_>) -> PyResult<PyObject> {
        let fraction = py.import("fractions")?.getattr("Fraction")?;
        let out = PyList::empty(py);
        for c in self.coords() {
            out.append(fraction.call1((c.n.to_big(), c.d.to_big()))?)?;
        }
        Ok(out.into_any().unbind())
    }
    #[pyo3(signature=(x=None))]
    fn as_expr(&self, x: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        Ok(w(self.represented(x)?))
    }
    #[pyo3(signature=(x=None))]
    fn as_poly(&self, x: Option<&Bound<'_, PyAny>>) -> PyResult<PyPoly> {
        let x = if x.is_some() {
            variable(x)?
        } else {
            self.alias.clone().unwrap_or_else(|| self.gen.clone())
        };
        Ok(PyPoly {
            p: polynomial(&self.value.coefficients, x).into(),
        })
    }
    fn field_element(&self, coeffs: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            value: AN::new(coefficients(coeffs)?, self.value.field.clone()),
            alias: self.alias.clone(),
            gen: self.gen.clone(),
        })
    }
    fn primitive_element(&self) -> Self {
        Self {
            value: AN::alpha(self.value.field.clone()),
            alias: self.alias.clone(),
            gen: self.gen.clone(),
        }
    }
    #[getter]
    fn is_primitive_element(&self) -> bool {
        self.value == AN::alpha(self.value.field.clone())
    }
    #[pyo3(signature=(radicals=true))]
    fn to_primitive_element(&self, radicals: bool) -> PyResult<Self> {
        if !radicals {
            return Err(PyNotImplementedError::new_err(
                "root-isolation representations are not implemented; use radicals=True",
            ));
        }
        if self.is_primitive_element() {
            return Ok(self.clone());
        }
        let field = NumberField::for_element(&self.value.as_expr()).map_err(err)?;
        Ok(Self {
            value: AN::alpha(field),
            alias: None,
            gen: Expr::symbol("_x"),
        })
    }
    fn minpoly_of_element(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.minimal_object(py, self.actual_minpoly()?)
    }
    #[pyo3(signature=(radicals=true,minpoly=None))]
    fn to_root(&self, radicals: bool, minpoly: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        if !radicals {
            return Err(PyNotImplementedError::new_err(
                "root-isolation representations are not implemented; use radicals=True",
            ));
        }
        if let Some(p) = minpoly.filter(|p| !p.is_none()) {
            supplied_minimal(p, &self.value.as_expr())?;
        }
        Ok(w(self.value.as_expr()))
    }
    fn _symox_(&self) -> PyExpr {
        w(self.value.as_expr())
    }
    #[getter]
    fn is_AlgebraicNumber(&self) -> bool {
        true
    }
    #[getter]
    fn is_Number(&self) -> bool {
        false
    }
    #[getter]
    fn is_Atom(&self) -> bool {
        false
    }
    #[getter]
    fn is_number(&self) -> bool {
        true
    }
    #[getter]
    fn is_algebraic(&self) -> bool {
        true
    }
    #[getter]
    fn is_commutative(&self) -> bool {
        true
    }
    #[getter]
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
    #[getter]
    fn is_rational(&self) -> bool {
        self.value.coefficients.len() <= 1
    }
    #[getter]
    fn is_integer(&self) -> bool {
        self.value.coefficients.is_empty()
            || (self.value.coefficients.len() == 1 && self.value.coefficients[0].d.is_one())
    }
    #[getter]
    fn is_nonzero(&self) -> Option<bool> {
        // SymPy's nonzero means real and nonzero, not merely different
        // from zero. Complex field elements must follow the scalar view.
        symox::assumptions::is_nonzero(&self.value.as_expr())
    }
    #[getter]
    fn is_finite(&self) -> bool {
        true
    }
    #[getter]
    fn is_real(&self) -> Option<bool> {
        symox::assumptions::is_real(&self.value.as_expr())
    }
    #[getter]
    fn is_positive(&self) -> Option<bool> {
        symox::assumptions::is_positive(&self.value.as_expr())
    }
    #[getter]
    fn is_negative(&self) -> Option<bool> {
        symox::assumptions::is_negative(&self.value.as_expr())
    }
    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PySet>> {
        PySet::empty(py)
    }
    #[getter]
    fn func(&self, py: Python<'_>) -> PyObject {
        py.get_type::<Self>().into_any().unbind()
    }
    #[getter]
    fn args(&self, py: Python<'_>) -> PyResult<PyObject> {
        let mut args = vec![
            w(self.value.field.generator.clone())
                .into_pyobject(py)?
                .into_any()
                .unbind(),
            PyTuple::new(py, self.coeffs())?.into_any().unbind(),
        ];
        if let Some(alias) = &self.alias {
            args.push(w(alias.clone()).into_pyobject(py)?.into_any().unbind());
        }
        Ok(PyTuple::new(py, args)?.into_any().unbind())
    }
    fn __str__(&self) -> String {
        symox::printer::to_str(&match &self.alias {
            Some(alias) => polynomial_expression(&self.value.coefficients, alias),
            None => symox::expand::expand(&polynomial_expression(
                &self.value.coefficients,
                &self.value.field.generator,
            )),
        })
    }
    fn __repr__(&self) -> String {
        self.__str__()
    }
    fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
        self.args(py)?.bind(py).hash()
    }
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .downcast::<Self>()
            .is_ok_and(|o| self.value == o.get().value && self.alias == o.get().alias)
    }
    fn __neg__(&self) -> PyExpr {
        w(ops::neg(&self.value.as_expr()))
    }
    fn __pos__(&self) -> Self {
        self.clone()
    }
    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(w(ops::add2(&self.value.as_expr(), &to_expr(other)?)))
    }
    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        self.__add__(other)
    }
    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(w(ops::sub(&self.value.as_expr(), &to_expr(other)?)))
    }
    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(w(ops::sub(&to_expr(other)?, &self.value.as_expr())))
    }
    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(w(ops::mul(&[self.value.as_expr(), to_expr(other)?])))
    }
    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        self.__mul__(other)
    }
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        let other = to_expr(other)?;
        if other.is_zero() {
            return Err(PyZeroDivisionError::new_err("division by zero"));
        }
        Ok(w(ops::mul(&[self.value.as_expr(), other.powi(-1)])))
    }
    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        if self.value.is_zero() {
            return Err(PyZeroDivisionError::new_err("division by zero"));
        }
        Ok(w(ops::mul(&[
            to_expr(other)?,
            self.value.as_expr().powi(-1),
        ])))
    }
    fn __pow__(
        &self,
        exponent: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyExpr> {
        if modulo.is_some() {
            return Err(PyTypeError::new_err(
                "modular exponentiation is not supported",
            ));
        }
        Ok(w(ops::pow(&self.value.as_expr(), &to_expr(exponent)?)))
    }
    #[pyo3(signature=(n=15,**kwargs))]
    fn evalf(&self, n: i64, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
        Ok(w(crate::evalf_expr(
            &self.value.as_expr(),
            n,
            None,
            kwargs,
        )?))
    }
    #[pyo3(signature=(n=15,**kwargs))]
    fn n(&self, n: i64, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
        self.evalf(n, kwargs)
    }
    fn __complex__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyComplex>> {
        let value = symox::evalf::eval_c64(&self.value.as_expr())
            .ok_or_else(|| PyTypeError::new_err("cannot evaluate this algebraic number"))?;
        Ok(pyo3::types::PyComplex::from_doubles(py, value.re, value.im))
    }
    fn __float__(&self) -> PyResult<f64> {
        symox::evalf::eval_f64(&self.value.as_expr()).ok_or_else(|| {
            PyTypeError::new_err("cannot convert a nonreal algebraic number to float")
        })
    }
    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let helper = py
            .import("symox._core")?
            .getattr("_rebuild_algebraic_number")?
            .unbind();
        let args = (
            w(self.value.field.generator.clone()),
            self.coeffs(),
            self.alias.clone().map(w),
            w(self.gen.clone()),
        )
            .into_pyobject(py)?
            .into_any()
            .unbind();
        Ok((helper, args))
    }
}
#[pyfunction]
fn _rebuild_algebraic_number(
    root: &Bound<'_, PyAny>,
    coeffs: &Bound<'_, PyAny>,
    alias: Option<&Bound<'_, PyAny>>,
    gen: &Bound<'_, PyAny>,
) -> PyResult<PyAlgebraicNumber> {
    PyAlgebraicNumber::build(root, Some(coeffs), alias, Some(gen))
}

#[pyfunction]
#[pyo3(signature=(extension,theta=None,*,gen=None,alias=None))]
fn to_number_field(
    extension: &Bound<'_, PyAny>,
    theta: Option<&Bound<'_, PyAny>>,
    gen: Option<&Bound<'_, PyAny>>,
    alias: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyAlgebraicNumber> {
    let expressions =
        if extension.downcast::<PyList>().is_ok() || extension.downcast::<PyTuple>().is_ok() {
            extension
                .try_iter()?
                .map(|e| to_expr(&e?))
                .collect::<PyResult<Vec<_>>>()?
        } else {
            vec![to_expr(extension)?]
        };
    if expressions.is_empty() {
        return Err(PyValueError::new_err(
            "cannot compute a primitive element for an empty extension",
        ));
    }
    let (root, _) = PrimitiveModel::build_many(&expressions).map_err(err)?;
    let mut number = if let Some(theta) = theta.filter(|t| !t.is_none()) {
        PyAlgebraicNumber::build(theta, None, alias, gen)?
    } else {
        PyAlgebraicNumber {
            value: AN::alpha(NumberField::for_element(&root).map_err(err)?),
            alias: alias_expr(alias)?,
            gen: variable(gen)?,
        }
    };
    number.value = AN::from_expr(&root, number.value.field.clone()).map_err(|e| {
        if symox::poly::embedding::is_absence(&e) {
            IsomorphismFailed::new_err(format!(
                "{} is not in a subfield of {}",
                symox::printer::to_str(&root),
                symox::printer::to_str(&number.value.field.generator)
            ))
        } else {
            err(e)
        }
    })?;
    Ok(number)
}
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlgebraicNumber>()?;
    m.add("IsomorphismFailed", m.py().get_type::<IsomorphismFailed>())?;
    m.add_function(wrap_pyfunction!(to_number_field, m)?)?;
    m.add_function(wrap_pyfunction!(_rebuild_algebraic_number, m)?)?;
    Ok(())
}
