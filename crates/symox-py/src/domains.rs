//! Typed coefficient domains: no display string participates in identity.
#![allow(non_snake_case)]
use super::{err, w, AlgebraicPoly, Domain, FractionPoly, NumberField, Poly, PolyValue};
use crate::{to_expr, PyExpr};
use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PySlice, PyTuple,PyDict};
use std::hash::{Hash, Hasher};
use std::sync::Arc;
use symox::poly::domains::{DomainElement, ElementValue, ExactDomain};
use symox::poly::int::Coeff;
use symox::Expr;

pyo3::create_exception!(symox, DomainError, PyValueError);
pyo3::create_exception!(symox, CoercionFailed, DomainError);
pyo3::create_exception!(symox, ExactQuotientFailed, DomainError);
fn element_error(e: symox::poly::convert::PolyError) -> PyErr {
    if matches!(&e, symox::poly::convert::PolyError::Failed(_)) {
        coercion(e)
    } else {
        err(e)
    }
}
fn coercion(e: impl std::fmt::Display) -> PyErr {
    let message = e.to_string();
    // A descriptor without an element backend (a fraction field over an
    // algebraic ground) is unsupported, not a failed coercion.
    if message.starts_with("fraction fields over ") && message.ends_with(" are not implemented") {
        return PyNotImplementedError::new_err(message);
    }
    CoercionFailed::new_err(message)
}
fn domain_error(e: impl std::fmt::Display) -> PyErr {
    DomainError::new_err(e.to_string())
}
/// Errors from building a composite domain: shared generators are the
/// reference `GeneratorsError`, unimplemented grounds stay unsupported.
fn composite_error(e: symox::poly::convert::PolyError) -> PyErr {
    match e {
        symox::poly::convert::PolyError::Generator(message) => Python::with_gil(|py| crate::poly::reference_error(py, "GeneratorsError", &message)),
        symox::poly::convert::PolyError::Unsupported(message) => PyNotImplementedError::new_err(message),
        other => domain_error(other),
    }
}
fn conversion_argument<'py>(args:&Bound<'py,PyTuple>,kwargs:Option<&Bound<'py,PyDict>>,alias:&str)->PyResult<Bound<'py,PyAny>> {
    if args.len()>1{return Err(PyTypeError::new_err("domain conversion takes exactly one argument"));}
    let mut result=args.get_item(0).ok();
    if let Some(kwargs)=kwargs{for (key,value) in kwargs.iter(){
        let key:String=key.extract()?;
        if key!="a" && key!=alias && key!="element"{return Err(PyTypeError::new_err(format!("domain conversion got an unexpected keyword argument '{key}'")));}
        if result.is_some(){return Err(PyTypeError::new_err("domain conversion got multiple values for its element"));}
        result=Some(value);
    }}
    result.ok_or_else(||PyTypeError::new_err("domain conversion requires an element"))
}
fn hashed<T: Hash>(v: &T) -> isize {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    v.hash(&mut h);
    let x = h.finish() as isize;
    if x == -1 {
        -2
    } else {
        x
    }
}
fn exact_numeric_parts(e:&Expr)->Option<(symox::number::Number,symox::number::Number)> {
    fn number(e:&Expr)->Option<symox::number::Number> {
        match e.as_number()? {
            symox::number::Number::Float(value) if !value.is_finite()=>None,
            number=>Some(number.clone()),
        }
    }
    let (re,im)=symox::complex_parts::as_real_imag_bounded(e,true,None).ok()??;
    Some((number(&re)?,number(&im)?))
}
fn numeric_hash(number:&symox::number::Number)->Option<isize> {
    match number {
        symox::number::Number::Float(value) if value.is_finite()=>crate::python_hash::exact_hash(
            &symox::number::Number::MpFloat(symox::bigfloat::BigFloat::from_f64(*value),53)),
        number=>crate::python_hash::exact_hash(number),
    }
}

#[pyclass(name = "Domain", module = "symox._core", frozen)]
#[derive(Clone, Debug)]
pub struct PyDomain {
    pub inner: ExactDomain,
}
impl std::fmt::Display for PyDomain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.inner)
    }
}
impl PyDomain {
    pub fn object(&self, py: Python<'_>) -> PyResult<PyObject> {
        match self.inner {
            ExactDomain::ZZ => Ok(py.import("symox._core")?.getattr("ZZ")?.unbind()),
            ExactDomain::QQ => Ok(py.import("symox._core")?.getattr("QQ")?.unbind()),
            ExactDomain::Expression {raw} => Ok(py.import("symox._core")?.getattr(if raw {"EXRAW"} else {"EX"})?.unbind()),
            ExactDomain::RR => Ok(py.import("symox._core")?.getattr("RR")?.unbind()),
            ExactDomain::Inexact {precision:53,complex:true} => Ok(py.import("symox._core")?.getattr("CC")?.unbind()),
            ExactDomain::Gaussian {integral} => Ok(py.import("symox._core")?.getattr(if integral {"ZZ_I"} else {"QQ_I"})?.unbind()),
            _ => Ok(self.clone().into_pyobject(py)?.into_any().unbind()),
        }
    }
    pub fn of(p: &PolyValue) -> PyResult<Self> {
        if let Some(e) = &p.expression {
            if let Some(g) = &e.ground {
                return Ok(Self{inner:g.domain()});
            }
            return Ok(Self{inner:ExactDomain::Expression{raw:false}});
        }
        if p.field.is_none() && p.algebraic.is_none() && p.modulus.is_none() && !p.base.domain_gens().is_empty() && p.domain != Domain::RR {
            // The legacy backend stores `ZZ[pi]` coefficients as extra
            // polynomial generators, which is the reference polynomial ring.
            return Ok(Self{inner:ExactDomain::Functions{integral:p.domain==Domain::ZZ,ring:true,parameters:p.base.domain_gens().to_vec(),order:"lex".into()}});
        }
        if p.field.is_none() && p.algebraic.is_none() && p.modulus.is_none() && p.domain == Domain::RR && !p.base.domain_gens().is_empty() {
            // Legacy floats over coefficient generators: `RR[y]`, or `CC`,
            // `CC[y]` when the imaginary unit is one of them.
            let i = Expr::i();
            let complex = p.base.domain_gens().contains(&i);
            let symbols: Vec<Expr> = p.base.domain_gens().iter().filter(|g| **g != i).cloned().collect();
            if symbols.is_empty() {
                return Ok(Self { inner: ExactDomain::Inexact { precision: 53, complex: true } });
            }
            let ground = symox::poly::gaussian_ground::GaussianGround { integral: false, ring: true, symbols, inexact: Some(symox::poly::inexact_ground::Inexact { precision: 53, complex }) };
            return Ok(Self { inner: ground.domain() });
        }
        if p.field.is_none() && p.algebraic.is_none() && !p.base.domain_gens().is_empty() {
            return Err(PyNotImplementedError::new_err(
                "this legacy coefficient backend has no implemented typed Domain descriptor",
            ));
        }
        Ok(Self {
            inner: if let Some(a) = &p.algebraic {
                match a.gaussian {
                    Some(integral) => ExactDomain::Gaussian { integral },
                    None => ExactDomain::Algebraic(a.number_field.clone()),
                }
            } else if let Some(f) = &p.field {
                f.coefficient_domain.clone()
            } else if let Some(p) = p.modulus {
                ExactDomain::GF(p)
            } else {
                match p.domain {
                    Domain::ZZ => ExactDomain::ZZ,
                    Domain::QQ => ExactDomain::QQ,
                    Domain::RR => ExactDomain::RR,
                }
            },
        })
    }
    pub fn poly(&self, e: &Expr, gens: &[Expr]) -> PyResult<PolyValue> {
        match &self.inner {
            ExactDomain::Expression { raw:false } => {
                if gens.is_empty() {
                    return Err(Python::with_gil(|py|crate::poly::reference_error(py,"GeneratorsError","you have to provide generators because EX domain was requested")));
                }
                symox::poly::expression_domain::ExpressionPoly::from_expr(e,gens).map(Into::into).map_err(crate::poly::err)
            }
            ExactDomain::Expression { .. } | ExactDomain::Inexact {..} => Err(PyNotImplementedError::new_err("this domain is not yet supported as a polynomial coefficient backend")),
            ExactDomain::Gaussian { integral } => {
                let mut p = AlgebraicPoly::from_expr(e, gens, symox::poly::gaussian::field())
                    .map_err(coercion)?;
                if *integral {
                    if let Some(message) = symox::poly::gaussian::integrality_error(&p.rep) {
                        return Err(coercion(message));
                    }
                }
                p.gaussian = Some(*integral);
                Ok(p.into())
            }
            ExactDomain::Algebraic(field) => AlgebraicPoly::from_expr(e, gens, field.clone())
                .map(Into::into)
                .map_err(coercion),
            ExactDomain::Composite { .. } if symox::poly::gaussian_ground::GaussianGround::from_domain(&self.inner).is_some() => {
                let ground = symox::poly::gaussian_ground::GaussianGround::from_domain(&self.inner).unwrap();
                if gens.is_empty() {
                    return Err(PyNotImplementedError::new_err(format!("polynomials over {} require generators", self.inner)));
                }
                ground.convert_expr(e, gens).map(Into::into).map_err(crate::poly::err)
            }
            ExactDomain::Composite { .. } => {
                // Nested rings over ZZ or QQ are stored over their flattened
                // generators (outermost first); the descriptor keeps the
                // nesting. Other composites have no polynomial backend yet.
                let Some(flat) = self.inner.flat_domain() else {
                    return Err(PyNotImplementedError::new_err(format!("polynomials over {} are not implemented", self.inner)));
                };
                let ExactDomain::Functions { parameters, integral, ring, .. } = &flat else { unreachable!() };
                let mut p = FractionPoly::from_expr(e, gens, parameters, &flat.to_string()).map_err(coercion)?;
                if *ring && p.rep.coeffs().iter().any(|c| !symox::poly::fraction::coefficient_in_ring(c, *integral)) {
                    return Err(coercion("coefficient outside selected ring"));
                }
                p.display = self.inner.to_string();
                p.coefficient_domain = self.inner.clone();
                Ok(p.into())
            }
            ExactDomain::Functions {
                parameters,
                ring,
                integral,
                ..
            } => {
                let mut p = FractionPoly::from_expr(e, gens, parameters, &self.inner.to_string())
                    .map_err(coercion)?;
                if *ring
                    && p.rep
                        .coeffs()
                        .iter()
                        .any(|c| !symox::poly::fraction::coefficient_in_ring(c, *integral))
                {
                    return Err(coercion("coefficient outside selected ring"));
                }
                p.coefficient_domain = self.inner.clone();
                Ok(p.into())
            }
            _ => {
                let mut p = Poly::from_expr(e, gens).map_err(coercion)?;
                if p.gens.len() != p.ngens {
                    return Err(coercion("coefficient contains undeclared generators"));
                }
                if matches!(self.inner, ExactDomain::ZZ) && !p.rep.is_integral() {
                    return Err(coercion("expected integer coefficients"));
                }
                p.domain = match self.inner {
                    ExactDomain::ZZ => Domain::ZZ,
                    ExactDomain::RR => Domain::RR,
                    _ => Domain::QQ,
                };
                if let ExactDomain::GF(prime) = self.inner {
                    if !p.rep.is_integral() {
                        return Err(coercion(
                            "finite-field conversion requires integral coefficients",
                        ));
                    }
                    p = p.set_modulus(prime);
                }
                Ok(p.into())
            }
        }
    }
    fn element(
        &self,
        obj: &Bound<'_, PyAny>,
        source: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyDomainElement> {
        if let Ok(value) = obj.extract::<PyRef<'_, PyDomainElement>>() {
            if let Some(source) = source {
                let source = source.extract::<PyRef<'_, Self>>()?;
                if source.inner != value.value.domain {
                    return Err(coercion("source domain does not match element"));
                }
            }
            if value.value.domain == self.inner {
                return Ok(value.clone());
            }
            if let Some(precision)=value.value.domain.inexact_precision() {
                if matches!(self.inner,ExactDomain::ZZ|ExactDomain::QQ) {
                    let rational=symox::poly::domains::rationalize_real(&value.value.as_expr(),precision,true).map_err(coercion)?;
                    return self.inner.convert(&Expr::num(rational.to_number())).map(|value|PyDomainElement{value}).map_err(coercion);
                }
                if matches!(self.inner,ExactDomain::Gaussian{..}) {
                    let numeric=symox::evalf_mp::eval_mp(&value.value.as_expr(),precision).ok_or_else(||coercion("expected a finite numerical value"))?;
                    let re=symox::evalf_mp::to_expr(symox::evalf_mp::MpComplex::real(numeric.re),precision);
                    let im=symox::evalf_mp::to_expr(symox::evalf_mp::MpComplex::real(numeric.im),precision);
                    let real=symox::poly::domains::rationalize_real(&re,precision,true).map_err(coercion)?;
                    let imag=symox::poly::domains::rationalize_real(&im,precision,true).map_err(coercion)?;
                    let expression=symox::ops::add2(&Expr::num(real.to_number()),&symox::ops::mul2(&Expr::i(),&Expr::num(imag.to_number())));
                    return self.inner.convert(&expression).map(|value|PyDomainElement{value}).map_err(coercion);
                }
            }
            // Finite residues are not canonically integers/rationals.
            if let ExactDomain::GF(prime) = value.value.domain {
                if matches!(self.inner, ExactDomain::ZZ) {
                    return self
                        .inner
                        .convert(&value.value.as_expr())
                        .map(|value| PyDomainElement { value })
                        .map_err(coercion);
                }
                if matches!(self.inner, ExactDomain::GF(_)) {
                    let ElementValue::Finite(a) = &value.value.value else {
                        unreachable!()
                    };
                    let n = a.to_rat().n.modulo(&symox::poly::Int::from_u64(prime));
                    return self
                        .inner
                        .convert(&Expr::num(n.to_number()))
                        .map(|value| PyDomainElement { value })
                        .map_err(coercion);
                }
                return Err(coercion("no canonical embedding from a finite field"));
            }
            return self
                .inner
                .convert(&value.value.as_expr())
                .map(|value| PyDomainElement { value })
                .map_err(coercion);
        }
        let e = if let Ok(text) = obj.extract::<String>() {
            if let Some(precision)=self.inner.inexact_precision() {
                match symox::bigfloat::BigFloat::from_decimal(&text,precision) {
                    Ok(value)=>Expr::num(symox::number::Number::mp_float(value,precision)),
                    Err(_)=>symox::parser::parse_str_with_precision(&text,precision).map_err(coercion)?,
                }
            } else {symox::parser::parse_str(&text).map_err(coercion)?}
        } else {
            to_expr(obj)?
        };
        self.inner
            .convert(&e)
            .map(|value| PyDomainElement { value })
            .map_err(|error| self.conversion_failure(&e, error))
    }
    /// The reference coercion message for an element proven outside an
    /// algebraic field; undecided memberships keep their own message.
    fn conversion_failure(&self, e: &Expr, error: symox::poly::convert::PolyError) -> PyErr {
        if matches!(error, symox::poly::convert::PolyError::Reference { .. }) {
            return err(error);
        }
        if matches!(self.inner, ExactDomain::Algebraic(_))
            && symox::poly::embedding::is_absence(&error)
        {
            return coercion(format!(
                "{} is not a valid algebraic number in {}",
                symox::printer::to_str(e),
                self.inner
            ));
        }
        coercion(error)
    }
    fn factory(&self, gens: &Bound<'_, PyTuple>, order: &str, ring: bool) -> PyResult<Self> {
        let parameters = super::gens_of(gens)?;
        if !matches!(self.inner,ExactDomain::ZZ|ExactDomain::QQ|ExactDomain::Functions{..}|ExactDomain::Composite{..}|ExactDomain::Gaussian{..}|ExactDomain::Algebraic(_)|ExactDomain::RR|ExactDomain::Inexact{..}) {
            return Err(PyNotImplementedError::new_err(format!("polynomial rings and fraction fields over {} are not implemented", self.inner)));
        }
        ExactDomain::composite(&self.inner, ring, parameters, order)
            .map(|inner| Self { inner })
            .map_err(composite_error)
    }
}
#[pymethods]
impl PyDomain {
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        self.inner.to_string()
    }
    fn __hash__(&self) -> isize {
        hashed(&self.inner)
    }
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp, py: Python<'_>) -> PyObject {
        let equal = other
            .extract::<PyRef<'_, Self>>()
            .is_ok_and(|o| o.inner == self.inner);
        match op {
            CompareOp::Eq => equal
                .into_pyobject(py)
                .unwrap()
                .to_owned()
                .into_any()
                .unbind(),
            CompareOp::Ne => (!equal)
                .into_pyobject(py)
                .unwrap()
                .to_owned()
                .into_any()
                .unbind(),
            _ => py.NotImplemented(),
        }
    }
    #[getter]
    fn is_ZZ(&self) -> bool {
        matches!(self.inner, ExactDomain::ZZ)
    }
    #[getter]
    fn is_QQ(&self) -> bool {
        matches!(self.inner, ExactDomain::QQ)
    }
    #[getter]
    fn is_RR(&self) -> bool {
        self.inner.inexact_precision().is_some()&&!self.inner.is_complex()
    }
    #[getter]
    fn is_CC(&self)->bool {self.inner.is_complex()}
    #[getter]
    fn is_ZZ_I(&self)->bool {matches!(self.inner,ExactDomain::Gaussian{integral:true})}
    #[getter]
    fn is_QQ_I(&self)->bool {matches!(self.inner,ExactDomain::Gaussian{integral:false})}
    #[getter]
    fn is_GaussianRing(&self)->bool {self.is_ZZ_I()}
    #[getter]
    fn is_GaussianField(&self)->bool {self.is_QQ_I()}
    #[getter]
    fn units(&self)->PyResult<Vec<PyDomainElement>> {
        if !self.is_ZZ_I()&&!self.is_QQ_I(){return Err(domain_error("the finite list of Gaussian units requires a Gaussian domain"));}
        [Expr::one(),Expr::i(),Expr::int(-1),symox::ops::neg(&Expr::i())].into_iter().map(|e|self.inner.convert(&e).map(|value|PyDomainElement{value}).map_err(coercion)).collect()
    }
    #[getter]
    fn is_RealField(&self)->bool {self.is_RR()}
    #[getter]
    fn is_ComplexField(&self)->bool {self.is_CC()}
    #[getter]
    fn precision(&self)->Option<u32> {self.inner.inexact_precision()}
    #[getter]
    fn dps(&self)->Option<u32> {self.precision().map(symox::bigfloat::prec_to_dps)}
    #[getter]
    fn is_FF(&self) -> bool {
        matches!(self.inner, ExactDomain::GF(_))
    }
    #[getter]
    fn is_FiniteField(&self) -> bool {
        self.is_FF()
    }
    #[getter]
    fn is_AlgebraicField(&self) -> bool {
        matches!(self.inner, ExactDomain::Algebraic(_))
    }
    #[getter]
    fn is_FractionField(&self) -> bool {
        matches!(self.inner, ExactDomain::Functions { ring: false, .. } | ExactDomain::Composite { ring: false, .. })
    }
    #[getter]
    fn is_PolynomialRing(&self) -> bool {
        matches!(self.inner, ExactDomain::Functions { ring: true, .. } | ExactDomain::Composite { ring: true, .. })
    }
    // The reference's alias flags.
    #[getter]
    fn is_Poly(&self) -> bool {
        self.is_PolynomialRing()
    }
    #[getter]
    fn is_Frac(&self) -> bool {
        self.is_FractionField()
    }
    #[getter]
    fn is_IntegerRing(&self) -> bool {
        self.is_ZZ()
    }
    #[getter]
    fn is_RationalField(&self) -> bool {
        self.is_QQ()
    }
    #[getter]
    fn is_Algebraic(&self) -> bool {
        self.is_AlgebraicField()
    }
    #[getter]
    fn is_SymbolicDomain(&self) -> bool {
        self.is_EX()
    }
    #[getter]
    fn is_SymbolicRawDomain(&self) -> bool {
        self.is_EXRAW()
    }
    #[getter]
    fn is_Field(&self) -> bool {
        self.inner.is_field()
    }
    #[getter]
    fn is_Ring(&self) -> bool {
        true
    }
    #[getter]
    fn is_Exact(&self) -> bool {
        self.inner.inexact_precision().is_none()
    }
    #[getter]
    fn is_Numerical(&self) -> bool {
        !matches!(self.inner, ExactDomain::Functions { .. } | ExactDomain::Composite { .. } | ExactDomain::Expression { .. })
    }
    #[getter]
    fn is_EX(&self)->bool {matches!(self.inner,ExactDomain::Expression{raw:false})}
    #[getter]
    fn is_EXRAW(&self)->bool {matches!(self.inner,ExactDomain::Expression{raw:true})}
    #[getter]
    fn has_CharacteristicZero(&self) -> bool {
        self.inner.characteristic() == 0
    }
    #[getter]
    fn is_Composite(&self) -> bool {
        matches!(self.inner, ExactDomain::Functions { .. } | ExactDomain::Composite { .. })
    }
    #[getter]
    fn is_PID(&self) -> bool {
        self.is_Field()
            || self.is_ZZ()
            || self.is_ZZ_I()
            || matches!(&self.inner,ExactDomain::Functions{integral:false,ring:true,parameters,..} if parameters.len()==1)
    }
    fn characteristic(&self) -> u64 {
        self.inner.characteristic()
    }
    #[getter]
    fn zero(&self) -> PyResult<PyDomainElement> {
        Ok(PyDomainElement {
            value: self.inner.convert(&Expr::zero()).map_err(err)?,
        })
    }
    #[getter]
    fn one(&self) -> PyResult<PyDomainElement> {
        Ok(PyDomainElement {
            value: self.inner.convert(&Expr::one()).map_err(err)?,
        })
    }
    #[getter]
    fn dtype(&self) -> Self {
        self.clone()
    }
    fn __instancecheck__(&self, obj: &Bound<'_, PyAny>) -> bool {
        obj.extract::<PyRef<'_, PyDomainElement>>()
            .is_ok_and(|e| e.value.domain == self.inner)
    }
    #[pyo3(signature=(value,denominator=None))]
    fn __call__(
        &self,
        value: &Bound<'_, PyAny>,
        denominator: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyDomainElement> {
        let a = self.element(value, None)?;
        if let Some(d) = denominator {
            if self.is_CC() || self.is_ZZ_I() || self.is_QQ_I() {
                let b=self.element(d,None)?;
                let precision=self.precision().unwrap_or(53);
                if !symox::evalf_mp::eval_mp(&a.value.as_expr(),precision).is_some_and(|v|v.is_real()) ||
                   !symox::evalf_mp::eval_mp(&b.value.as_expr(),precision).is_some_and(|v|v.is_real()) {
                    return Err(coercion("complex real and imaginary components must be real"));
                }
                let e=symox::ops::add2(&a.value.as_expr(),&symox::ops::mul2(&Expr::i(),&b.value.as_expr()));
                self.inner.convert(&e).map(|value|PyDomainElement{value}).map_err(coercion)
            } else {a.binary(d, '/', false)}
        } else {
            Ok(a)
        }
    }
    #[pyo3(signature=(element,base=None))]
    fn convert(
        &self,
        element: &Bound<'_, PyAny>,
        base: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyDomainElement> {
        self.element(element, base)
    }
    #[pyo3(signature=(*args,**kwargs))]
    fn from_sympy(&self, args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>) -> PyResult<PyDomainElement> {
        self.element(&conversion_argument(args,kwargs,"expr")?, None)
    }
    #[pyo3(signature=(*args,**kwargs))]
    fn to_sympy(&self, args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>) -> PyResult<PyExpr> {
        Ok(w(self.element(&conversion_argument(args,kwargs,"element")?, None)?.value.as_expr()))
    }
    fn of_type(&self, element: &Bound<'_, PyAny>) -> bool {
        self.__instancecheck__(element)
    }
    #[pyo3(signature=(element,limit=true))]
    fn to_rational(&self,element:&Bound<'_,PyAny>,limit:bool)->PyResult<(num_bigint::BigInt,num_bigint::BigInt)> {
        let precision=self.precision().ok_or_else(||domain_error("to_rational requires a real or complex numerical field"))?;
        let a=self.element(element,None)?;
        let rational=symox::poly::domains::rationalize_real(&a.value.as_expr(),precision,limit).map_err(coercion)?;
        Ok((rational.n.to_big(),rational.d.to_big()))
    }
    fn __contains__(&self, element: &Bound<'_, PyAny>) -> PyResult<bool> {
        match self.element(element, None) {
            Ok(_) => Ok(true),
            // An uncertified membership is not a proven absence.
            Err(e) if e.to_string().contains("could not certify") => Err(PyNotImplementedError::new_err(e.to_string())),
            Err(_) => Ok(false),
        }
    }
    pub fn get_field(&self) -> Self {
        Self {
            inner: self.inner.field(),
        }
    }
    pub fn get_ring(&self) -> PyResult<Self> {
        self.inner
            .ring()
            .map(|inner| Self { inner })
            .map_err(domain_error)
    }
    #[getter]
    fn dom(&self) -> PyResult<Self> {
        match self.inner {
            ExactDomain::Functions { integral, .. } => Ok(Self {
                inner: if integral {
                    ExactDomain::ZZ
                } else {
                    ExactDomain::QQ
                },
            }),
            ExactDomain::Algebraic(_) => Ok(Self {
                inner: ExactDomain::QQ,
            }),
            ExactDomain::Gaussian {integral} => Ok(Self{inner:if integral {ExactDomain::ZZ} else {ExactDomain::QQ}}),
            ExactDomain::Composite { ref ground, .. } => Ok(Self { inner: (**ground).clone() }),
            _ => Err(domain_error(
                "ground domains do not have a coefficient domain",
            )),
        }
    }
    #[getter]
    fn domain(&self) -> PyResult<Self> {
        self.dom()
    }
    #[getter]
    fn order(&self) -> Option<String> {
        match &self.inner {
            ExactDomain::Functions { order, .. } | ExactDomain::Composite { order, .. } => Some(order.clone()),
            _ => None,
        }
    }
    #[getter]
    fn symbols(&self, py: Python<'_>) -> PyResult<PyObject> {
        match &self.inner {
            ExactDomain::Functions { parameters, .. } | ExactDomain::Composite { parameters, .. } => {
                Ok(PyTuple::new(py, parameters.iter().cloned().map(w))?
                    .into_any()
                    .unbind())
            }
            ExactDomain::Algebraic(_) => Ok(PyTuple::new(py, [self.ext(py)?])?.into_any().unbind()),
            _ => Ok(PyTuple::empty(py).into_any().unbind()),
        }
    }
    #[getter]
    fn gens(&self, py: Python<'_>) -> PyResult<PyObject> {
        match &self.inner {
            ExactDomain::Functions { parameters, .. } | ExactDomain::Composite { parameters, .. } => Ok(PyTuple::new(
                py,
                parameters
                    .iter()
                    .map(|p| self.inner.convert(p).map(|value| PyDomainElement { value }))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(coercion)?,
            )?
            .into_any()
            .unbind()),
            _ => self.symbols(py),
        }
    }
    #[getter]
    fn ngens(&self) -> usize {
        match &self.inner {
            ExactDomain::Functions { parameters, .. } | ExactDomain::Composite { parameters, .. } => parameters.len(),
            ExactDomain::Algebraic(_) => 1,
            _ => 0,
        }
    }
    #[getter]
    fn ext(&self, py: Python<'_>) -> PyResult<PyObject> {
        let ExactDomain::Algebraic(k) = &self.inner else {
            return Err(domain_error("this domain is not an algebraic field"));
        };
        let args = PyTuple::new(py, [w(k.generator.clone())])?;
        let kwargs = pyo3::types::PyDict::new(py);
        if let Some(alias) = &k.alias {
            kwargs.set_item("alias", w(alias.clone()))?;
        }
        Ok(py
            .import("symox._core")?
            .getattr("AlgebraicNumber")?
            .call(&args, Some(&kwargs))?
            .unbind())
    }
    /// The defining polynomial of an algebraic field, the representation
    /// of the minimal polynomial of its primitive element.
    #[getter]
    fn r#mod(&self, py: Python<'_>) -> PyResult<PyObject> {
        Ok(self.ext(py)?.bind(py).getattr("minpoly")?.getattr("rep")?.unbind())
    }
    #[getter]
    fn unit(&self) -> PyResult<PyDomainElement> {
        let ExactDomain::Algebraic(k) = &self.inner else {
            return Err(domain_error("unit is the primitive algebraic generator"));
        };
        self.inner
            .convert(&k.generator)
            .map(|value| PyDomainElement { value })
            .map_err(coercion)
    }
    #[pyo3(signature=(*extension,alias=None))]
    fn algebraic_field(
        &self,
        extension: &Bound<'_, PyTuple>,
        alias: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut generators = super::gens_of(extension)?;
        match &self.inner {
            ExactDomain::QQ | ExactDomain::ZZ => {}
            // The reference extends an algebraic field by its generator.
            ExactDomain::Algebraic(k) => generators.insert(0, k.generator.clone()),
            _ => {
                return Err(PyNotImplementedError::new_err(
                    "algebraic fields currently require QQ as their ground field",
                ))
            }
        }
        let field = NumberField::from_generators(&generators).map_err(|e| {
            for g in &generators {
                if let Err(m) = symox::poly::embedding::minimal_rational(g) {
                    return super::minpoly_error(extension.py(), m);
                }
            }
            err(e)
        })?;
        let mut k = (*field).clone();
        if let Some(alias) = alias {
            if !alias.is_none() {
                let a = if let Ok(name) = alias.extract::<String>() {
                    Expr::symbol(&name)
                } else {
                    to_expr(alias)?
                };
                if !a.is_symbol() {
                    return Err(PyValueError::new_err("alias must be a symbol"));
                }
                k.alias = Some(a);
            }
        }
        Ok(Self {
            inner: ExactDomain::Algebraic(Arc::new(k)),
        })
    }
    #[pyo3(signature=(*symbols,order="lex"))]
    fn frac_field(&self, symbols: &Bound<'_, PyTuple>, order: &str) -> PyResult<Self> {
        self.factory(symbols, order, false)
    }
    #[pyo3(signature=(*symbols,order="lex"))]
    fn poly_ring(&self, symbols: &Bound<'_, PyTuple>, order: &str) -> PyResult<Self> {
        self.factory(symbols, order, true)
    }
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if key.downcast::<PySlice>().is_ok() || key.extract::<isize>().is_ok() {
            return Ok(self.__str__().into_pyobject(py)?.get_item(key)?.unbind());
        }
        let args = if let Ok(t) = key.downcast::<PyTuple>() {
            t.clone()
        } else {
            PyTuple::new(py, [key])?
        };
        Ok(self
            .factory(&args, "lex", true)?
            .into_pyobject(py)?
            .into_any()
            .unbind())
    }
    fn startswith(&self, prefix: &str) -> bool {
        self.__str__().starts_with(prefix)
    }
    fn endswith(&self, suffix: &str) -> bool {
        self.__str__().ends_with(suffix)
    }
    fn exquo(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        let a = self.element(a, None)?;
        let b = self.element(b, None)?;
        a.value
            .exquo(&b.value)
            .map(|value| PyDomainElement { value })
            .map_err(err)
    }
    fn quo(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        Ok(self.div(a, b)?.0)
    }
    fn rem(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        Ok(self.div(a, b)?.1)
    }
    fn div(
        &self,
        a: &Bound<'_, PyAny>,
        b: &Bound<'_, PyAny>,
    ) -> PyResult<(PyDomainElement, PyDomainElement)> {
        let a = self.element(a, None)?;
        let b = self.element(b, None)?;
        let (q, r) = a.value.div_rem(&b.value).map_err(err)?;
        Ok((PyDomainElement { value: q }, PyDomainElement { value: r }))
    }
    fn gcd(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        let a = self.element(a, None)?;
        let b = self.element(b, None)?;
        a.value
            .gcd(&b.value)
            .map(|value| PyDomainElement { value })
            .map_err(err)
    }
    fn lcm(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        let a = self.element(a, None)?;
        let b = self.element(b, None)?;
        a.value
            .lcm(&b.value)
            .map(|value| PyDomainElement { value })
            .map_err(err)
    }
    fn numer(&self, a: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        self.element(a, None)?.fraction_part(false)
    }
    fn denom(&self, a: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        self.element(a, None)?.fraction_part(true)
    }
    fn add(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        self.element(a, None)?.binary(b, '+', false)
    }
    fn sub(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        self.element(a, None)?.binary(b, '-', false)
    }
    fn mul(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        self.element(a, None)?.binary(b, '*', false)
    }
    fn neg(&self, a: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        self.element(a, None)?.__neg__()
    }
    fn pow(&self, a: &Bound<'_, PyAny>, n: i64) -> PyResult<PyDomainElement> {
        self.element(a, None)?
            .value
            .pow(n)
            .map(|value| PyDomainElement { value })
            .map_err(err)
    }
    fn invert(&self, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyDomainElement> {
        let a = self.element(a, None)?;
        let b = self.element(b, None)?;
        if let (ExactDomain::ZZ, ElementValue::Rational(a), ElementValue::Rational(b)) =
            (&self.inner, &a.value.value, &b.value.value)
        {
            if b.is_zero() {
                return Err(PyZeroDivisionError::new_err("zero modulus"));
            }
            let (g, s, _) = a.n.gcdex(&b.n);
            if !g.is_one() {
                return Err(super::an_rational::NotInvertible::new_err("zero divisor"));
            }
            return self
                .inner
                .convert(&Expr::num(s.modulo(&b.n).to_number()))
                .map(|value| PyDomainElement { value })
                .map_err(err);
        }
        if matches!(self.inner,ExactDomain::Functions{integral:true,ring:true,..}) {return Err(DomainError::new_err("inversion over a polynomial ring requires a field as its ground domain"));}
        if let ExactDomain::Functions {
            ring: true,
            parameters,
            ..
        } = &self.inner
        {
            let e = symox::poly::api::invert(&a.value.as_expr(), &b.value.as_expr(), parameters)
                .map_err(err)?;
            return self
                .inner
                .convert(&e)
                .map(|value| PyDomainElement { value })
                .map_err(coercion);
        }
        Err(PyNotImplementedError::new_err("modular inversion requires ZZ or a supported univariate polynomial ring; use a**-1 for field inversion"))
    }
    fn gcdex(
        &self,
        a: &Bound<'_, PyAny>,
        b: &Bound<'_, PyAny>,
    ) -> PyResult<(PyDomainElement, PyDomainElement, PyDomainElement)> {
        let a = self.element(a, None)?;
        let b = self.element(b, None)?;
        if let ExactDomain::Functions {
            integral,
            ring: true,
            parameters,
            ..
        } = &self.inner
        {
            if *integral {
                return Err(domain_error("extended gcd requires a field ground domain"));
            }
            let (s, t, g) =
                symox::poly::api::gcdex(&a.value.as_expr(), &b.value.as_expr(), parameters)
                    .map_err(err)?;
            let wrap = |e: &Expr| {
                self.inner
                    .convert(e)
                    .map(|value| PyDomainElement { value })
                    .map_err(coercion)
            };
            return Ok((wrap(&s)?, wrap(&t)?, wrap(&g)?));
        }
        // The reference integer gcdex (gmpy `gcdext`) gives (0, 0, 0) at (0, 0).
        if matches!(self.inner, ExactDomain::ZZ) && a.value.is_zero() && b.value.is_zero() {
            return Ok((a.clone(), a.clone(), a));
        }
        let (s, t, g) = a.value.gcdex(&b.value).map_err(err)?;
        Ok((
            PyDomainElement { value: s },
            PyDomainElement { value: t },
            PyDomainElement { value: g },
        ))
    }
    #[pyo3(signature=(K1,symbols=None))]
    fn unify(&self, K1: PyRef<'_, Self>,symbols:Option<&Bound<'_,PyAny>>) -> PyResult<Self> {
        if let Some(symbols)=symbols.filter(|v|!v.is_none()){
            let symbols=symbols.try_iter()?.map(|v|to_expr(&v?)).collect::<PyResult<Vec<_>>>()?;
            if [&self.inner,&K1.inner].iter().any(|domain|matches!(domain,ExactDomain::Functions{parameters,..} if parameters.iter().any(|p|symbols.contains(p)))){
                let error=K1.py().import("symox")?.getattr("UnificationFailed")?.call1(("domain generators overlap polynomial generators",))?;
                return Err(PyErr::from_value(error));
            }
        }
        unify_native(&self.inner, &K1.inner).map(|inner| Self { inner })
    }
    fn is_zero(&self, a: &Bound<'_, PyAny>) -> PyResult<bool> {
        Ok(self.element(a, None)?.value.is_zero())
    }
    fn is_one(&self, a: &Bound<'_, PyAny>) -> PyResult<bool> {
        Ok(self.element(a, None)?.value == self.one()?.value)
    }
    fn is_unit(&self, a: &Bound<'_, PyAny>) -> PyResult<bool> {
        let polynomial_ring = matches!(self.inner, ExactDomain::Composite { ring: true, .. });
        // PolynomialRing.is_unit reads `a.is_ground` before converting.
        if polynomial_ring && a.extract::<PyRef<'_, PyDomainElement>>().is_err() && !a.hasattr("is_ground")? {
            return Err(pyo3::exceptions::PyAttributeError::new_err(format!(
                "'{}' object has no attribute 'is_ground'",
                a.get_type().name()?
            )));
        }
        let a = self.element(a, None)?;
        if self.is_Field() {
            return Ok(!a.value.is_zero());
        }
        if matches!(self.inner, ExactDomain::Composite { .. } | ExactDomain::Gaussian { .. }) {
            return Ok(ring_unit(&self.inner, &a.value.as_expr()));
        }
        Ok(a.value.pow(-1).is_ok())
    }
    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        if let ExactDomain::Composite { ground, ring, parameters, order } = &self.inner {
            // Rebuilt from its ground domain's factory, like the reference.
            let ground = Self { inner: (**ground).clone() }.object(py)?;
            let factory = ground.bind(py).getattr(if *ring { "poly_ring" } else { "frac_field" })?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("order", order)?;
            let partial = py.import("functools")?.getattr("partial")?.call((factory,), Some(&kwargs))?;
            return Ok((partial.unbind(), PyTuple::new(py, parameters.iter().cloned().map(w))?.into_any().unbind()));
        }
        let (kind, prime, parameters, generator, alias, integral, ring, order) = match &self.inner {
            ExactDomain::Expression {raw} => (if *raw {"EXRAW"} else {"EX"},None,vec![],None,None,false,false,"lex".to_string()),
            ExactDomain::Inexact {precision,complex} => (if *complex {"Complex"} else {"Real"},Some(*precision as u64),vec![],None,None,false,false,"lex".to_string()),
            ExactDomain::Gaussian {integral} => (if *integral {"ZZ_I"} else {"QQ_I"},None,vec![],None,None,false,false,"lex".to_string()),
            ExactDomain::RR => (
                "RR",
                None,
                vec![],
                None,
                None,
                false,
                false,
                "lex".to_string(),
            ),
            ExactDomain::ZZ => (
                "ZZ",
                None,
                vec![],
                None,
                None,
                false,
                false,
                "lex".to_string(),
            ),
            ExactDomain::QQ => (
                "QQ",
                None,
                vec![],
                None,
                None,
                false,
                false,
                "lex".to_string(),
            ),
            ExactDomain::GF(p) => (
                "GF",
                Some(*p),
                vec![],
                None,
                None,
                false,
                false,
                "lex".to_string(),
            ),
            ExactDomain::Algebraic(k) => (
                "AN",
                None,
                vec![],
                Some(w(k.generator.clone())),
                k.alias.clone().map(w),
                false,
                false,
                "lex".to_string(),
            ),
            ExactDomain::Composite { .. } => unreachable!("composite domains are reduced above"),
            ExactDomain::Functions {
                integral,
                ring,
                parameters,
                order,
            } => (
                "RF",
                None,
                parameters.iter().cloned().map(w).collect(),
                None,
                None,
                *integral,
                *ring,
                order.clone(),
            ),
        };
        Ok((
            py.import("symox._core")?
                .getattr("_restore_domain")?
                .unbind(),
            (
                kind, prime, parameters, generator, alias, integral, ring, order,
            )
                .into_pyobject(py)?
                .into_any()
                .unbind(),
        ))
    }
}

#[pyclass(name = "DomainElement", module = "symox._core", frozen)]
#[derive(Clone)]
pub struct PyDomainElement {
    pub value: DomainElement,
}
impl PyDomainElement {
    fn gaussian_coordinate(&self,imaginary:bool)->PyResult<Self> {
        let (x,y)=self.value.gaussian_parts().ok_or_else(||PyTypeError::new_err("coordinates require a Gaussian element"))?;
        let domain=if matches!(self.value.domain,ExactDomain::Gaussian{integral:true}){ExactDomain::ZZ}else{ExactDomain::QQ};
        domain.convert(&Expr::num(if imaginary {y}else{x}.to_number())).map(|value|Self{value}).map_err(coercion)
    }
    fn numerical_coordinate(&self,imaginary:bool)->PyResult<Self> {
        if let Some(precision)=self.value.domain.inexact_precision() {
            let value=symox::evalf_mp::eval_mp(&self.value.as_expr(),precision).ok_or_else(||coercion("cannot evaluate a numerical coordinate"))?;
            let component=if imaginary{value.im}else{value.re};
            let domain=if precision==53 {ExactDomain::RR}else{ExactDomain::Inexact{precision,complex:false}};
            return domain.convert(&symox::evalf_mp::to_expr(symox::evalf_mp::MpComplex::real(component),precision)).map(|value|Self{value}).map_err(coercion);
        }
        if self.value.gaussian_parts().is_some(){return self.gaussian_coordinate(imaginary);}
        if imaginary {self.value.domain.convert(&Expr::zero()).map(|value|Self{value}).map_err(coercion)} else {Ok(self.clone())}
    }
    fn fraction_part(&self, denominator: bool) -> PyResult<Self> {
        if let ElementValue::GaussianFraction(p, q) = &self.value.value {
            let domain = self.value.domain.ring().map_err(err)?;
            let side = if denominator { q } else { p };
            return Ok(Self { value: DomainElement { domain, value: ElementValue::AlgebraicPoly(side.clone()) } });
        }
        if matches!(self.value.domain,ExactDomain::Gaussian{..}) {
            let (n,d)=self.value.gaussian_fraction().map_err(err)?;
            return Ok(Self{value:if denominator {d}else{n}});
        }
        let (domain, e) = match (&self.value.domain, &self.value.value) {
            (ExactDomain::Functions { ring: true, .. } | ExactDomain::Composite { ring: true, .. }, _) => {
                return if denominator {
                    self.value
                        .domain
                        .convert(&Expr::one())
                        .map(|value| Self { value })
                        .map_err(err)
                } else {
                    Ok(self.clone())
                };
            }
            (ExactDomain::ZZ | ExactDomain::QQ, ElementValue::Rational(r)) => (
                ExactDomain::ZZ,
                Expr::num(if denominator {
                    r.d.to_number()
                } else {
                    r.n.to_number()
                }),
            ),
            (
                ExactDomain::Functions {
                    integral,
                    parameters,
                    order,
                    ..
                },
                ElementValue::Function(r),
            ) => (
                ExactDomain::functions(*integral, true, parameters.clone(), order).map_err(err)?,
                symox::poly::convert::to_expr_int(
                    if denominator { &r.d } else { &r.n },
                    parameters,
                ),
            ),
            (domain @ ExactDomain::Composite { ring: false, .. }, ElementValue::Function(r)) => (
                domain.ring().map_err(err)?,
                symox::poly::convert::to_expr_int(if denominator { &r.d } else { &r.n }, &domain.flat_parameters()),
            ),
            (_, _) => (
                self.value.domain.clone(),
                if denominator {
                    Expr::one()
                } else {
                    self.value.as_expr()
                },
            ),
        };
        domain
            .convert(&e)
            .map(|value| Self { value })
            .map_err(coercion)
    }
    fn divide(&self, other: &Bound<'_, PyAny>, reverse: bool) -> PyResult<(Self, Self)> {
        if let Ok(b) = other.extract::<PyRef<'_, Self>>() {
            if !self.value.compatible_embedding(&b.value) {
                return Err(coercion("incompatible coefficient domains"));
            }
        }
        let b = self.parent().element(other, None)?;
        let (q, r) = if reverse {
            b.value.div_rem(&self.value)
        } else {
            self.value.div_rem(&b.value)
        }
        .map_err(err)?;
        Ok((Self { value: q }, Self { value: r }))
    }
    fn ring_operator(
        &self,
        other: &Bound<'_, PyAny>,
        operation: char,
        reverse: bool,
    ) -> PyResult<PyObject> {
        if !matches!(
            self.value.domain,
            ExactDomain::ZZ | ExactDomain::Functions { ring: true, .. } | ExactDomain::Gaussian {integral:true} | ExactDomain::Composite { ring: true, .. }
        ) {
            return Ok(other.py().NotImplemented());
        }
        let (q, r) = self.divide(other, reverse)?;
        if operation == 'q' {
            Ok(q.into_pyobject(other.py())?.into_any().unbind())
        } else if operation == 'r' {
            Ok(r.into_pyobject(other.py())?.into_any().unbind())
        } else {
            Ok((q, r).into_pyobject(other.py())?.into_any().unbind())
        }
    }
    fn operation(
        &self,
        other: &Bound<'_, PyAny>,
        operation: char,
        reverse: bool,
    ) -> PyResult<PyObject> {
        // Let matrix protocols apply a typed scalar on the left, rather
        // than trying to convert the whole matrix into one coefficient.
        if other.getattr("_symox_domainmatrix").is_ok_and(|v|v.is_truthy().unwrap_or(false))
            || other.getattr("_symox_quotient_element").is_ok_and(|v|v.is_truthy().unwrap_or(false))
            || other.getattr("is_Matrix").is_ok_and(|v|v.is_truthy().unwrap_or(false))
            || other.getattr("_symox_sparse").is_ok_and(|v|v.is_truthy().unwrap_or(false))
            // A ring or field element over this domain applies the scalar itself.
            || other.getattr("_symox_ring_element").is_ok_and(|v|v.is_truthy().unwrap_or(false))
            || other.downcast::<crate::matrices::PyMatrix>().is_ok() {
            return Ok(other.py().NotImplemented());
        }
        if other
            .extract::<PyRef<'_, crate::PyExprNode>>()
            .is_ok_and(|e| !e.e.is_number())
            && matches!(
                self.value.domain,
                ExactDomain::ZZ | ExactDomain::QQ | ExactDomain::GF(_) | ExactDomain::RR | ExactDomain::Inexact {..} | ExactDomain::Gaussian {..}
            )
        {
            return Ok(other.py().NotImplemented());
        }
        match self.binary(other, operation, reverse) {
            Ok(value)=>Ok(value.into_pyobject(other.py())?.into_any().unbind()),
            Err(error) if error.is_instance_of::<PyTypeError>(other.py())=>Ok(other.py().NotImplemented()),
            Err(error)=>Err(error),
        }
    }
    fn binary(&self, other: &Bound<'_, PyAny>, operation: char, reverse: bool) -> PyResult<Self> {
        if let Ok(b) = other.extract::<PyRef<'_, Self>>() {
            if !self.value.compatible_embedding(&b.value)
                && !(matches!(self.value.value, ElementValue::Rational(_))
                    && matches!(b.value.value, ElementValue::Rational(_)))
                && !(matches!(self.value.domain,ExactDomain::Gaussian{..})&&matches!(b.value.domain,ExactDomain::Gaussian{..}))
            {
                return Err(coercion(
                    "incompatible coefficient domains; use convert explicitly",
                ));
            }
        }
        let d = PyDomain {
            inner: self.value.domain.clone(),
        };
        let b = d.element(other, None)?;
        let value = if reverse {
            b.value.binary(&self.value, operation)
        } else {
            self.value.binary(&b.value, operation)
        }
        .map_err(element_error)?;
        Ok(Self { value })
    }
}
#[pymethods]
impl PyDomainElement {
    #[getter]
    fn ex(&self)->PyExpr {self.as_expr()}
    /// The reference's ``PolyElement.is_term``: at most one term. Only
    /// polynomial ring elements define it.
    #[getter]
    fn is_term(&self)->PyResult<bool> {
        match (&self.value.domain,&self.value.value) {
            (ExactDomain::Functions{ring:true,..},ElementValue::Function(f))=>Ok(f.n.len()<=1),
            _=>Err(pyo3::exceptions::PyAttributeError::new_err(format!("'{}' object has no attribute 'is_term'",self.value.domain))),
        }
    }
    #[getter]
    fn x(&self)->PyResult<Self> {self.gaussian_coordinate(false)}
    #[getter]
    fn y(&self)->PyResult<Self> {self.gaussian_coordinate(true)}
    #[getter]
    fn real(&self)->PyResult<Self> {self.numerical_coordinate(false)}
    #[getter]
    fn imag(&self)->PyResult<Self> {self.numerical_coordinate(true)}
    fn __abs__(&self)->PyResult<Self> {
        let domain=if self.value.domain.is_complex() {
            let precision=self.value.domain.inexact_precision().unwrap();
            if precision==53 {ExactDomain::RR}else{ExactDomain::Inexact{precision,complex:false}}
        }else{self.value.domain.clone()};
        domain.convert(&symox::functions::abs(&self.value.as_expr())).map(|value|Self{value}).map_err(coercion)
    }
    fn quadrant(&self)->PyResult<usize> {
        let (x,y)=self.value.gaussian_parts().ok_or_else(||PyTypeError::new_err("quadrants require a Gaussian element"))?;
        Ok(if !y.is_zero()&&!y.is_negative(){if !x.is_zero()&&!x.is_negative(){0}else{1}}
            else if y.is_negative(){if x.is_negative(){2}else{3}}else if x.is_negative(){2}else{0})
    }
    #[getter]
    fn numerator(&self) -> PyResult<Self> {
        self.fraction_part(false)
    }
    #[getter]
    fn denominator(&self) -> PyResult<Self> {
        self.fraction_part(true)
    }
    #[getter]
    fn numer(&self) -> PyResult<Self> {
        self.fraction_part(false)
    }
    #[getter]
    fn denom(&self) -> PyResult<Self> {
        self.fraction_part(true)
    }
    fn to_list(&self) -> PyResult<Vec<Self>> {
        let ElementValue::Algebraic(a) = &self.value.value else {
            return Err(PyTypeError::new_err(
                "power-basis coefficients are only defined for algebraic elements",
            ));
        };
        Ok(a.coefficients
            .iter()
            .rev()
            .map(|r| Self {
                value: ExactDomain::QQ.convert(&Expr::num(r.to_number())).unwrap(),
            })
            .collect())
    }
    fn parent(&self) -> PyDomain {
        PyDomain {
            inner: self.value.domain.clone(),
        }
    }
    fn as_expr(&self) -> PyExpr {
        w(self.value.as_expr())
    }
    fn __str__(&self) -> String {
        if let Some(text) = symox::poly::domain_print::composite_str(&self.value) {
            return text;
        }
        self.value.as_expr().to_string()
    }
    fn __repr__(&self) -> String {
        // Composite elements print like the reference `PolyElement`.
        if let Some(text) = symox::poly::domain_print::composite_str(&self.value) {
            return text;
        }
        // The reference GaussianElement repr: `ZZ_I(x, y)`.
        if matches!(self.value.domain,ExactDomain::Gaussian{..}) {
            if let Some((x,y))=exact_numeric_parts(&self.value.as_expr()) {
                return format!("{}({}, {})", self.value.domain, Expr::num(x), Expr::num(y));
            }
        }
        format!("{}({})", self.value.domain, self.value.as_expr())
    }
    fn __bool__(&self) -> bool {
        !self.value.is_zero()
    }
    fn __hash__(&self) -> isize {
        if matches!(self.value.value,ElementValue::Inexact(_)) {
            if let Some((real,imag))=exact_numeric_parts(&self.value.as_expr()) {
                let re=numeric_hash(&real).unwrap();
                let im=numeric_hash(&imag).unwrap();
                let hash=re.wrapping_add(im.wrapping_mul(1_000_003));
                return if hash== -1 {-2} else {hash};
            }
        }
        if let ElementValue::Rational(r) = &self.value.value {
            crate::python_hash::exact_hash(&r.to_number()).unwrap()
        } else {
            hashed(&self.value)
        }
    }
    fn __int__(&self) -> PyResult<num_bigint::BigInt> {
        match &self.value.value {
            ElementValue::Inexact(e) => {
                let value=symox::evalf_mp::eval_mp(e,self.value.domain.inexact_precision().unwrap()).ok_or_else(||coercion("cannot evaluate this numerical element"))?;
                if !value.is_real(){return Err(PyTypeError::new_err("cannot convert a nonreal coefficient to int"));}
                value.re.checked_int_trunc().ok_or_else(||pyo3::exceptions::PyOverflowError::new_err("too many digits in integer"))
            },
            ElementValue::Rational(r) if r.is_integer() => Ok(r.n.to_big()),
            ElementValue::Finite(r) => {
                let r = r.to_rat();
                let ExactDomain::GF(p) = self.value.domain else {
                    unreachable!()
                };
                Ok(r.n.modulo(&symox::poly::Int::from_u64(p)).to_big())
            }
            _ => Err(PyTypeError::new_err(
                "only integral coefficients convert to int",
            )),
        }
    }
    fn __index__(&self) -> PyResult<num_bigint::BigInt> {
        if matches!(self.value.value,ElementValue::Inexact(_)){return Err(PyTypeError::new_err("inexact coefficients cannot be used as integer indices"));}
        self.__int__()
    }
    fn __float__(&self) -> PyResult<f64> {
        match &self.value.value {
            ElementValue::Rational(r) => Ok(r.to_f64()),
            ElementValue::Inexact(e) | ElementValue::Expression(e) => symox::evalf::eval_f64(e).ok_or_else(||PyTypeError::new_err("coefficient has no real numerical value")),
            _ => Err(PyTypeError::new_err(
                "this exact coefficient has no real float conversion",
            )),
        }
    }
    fn __complex__<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,pyo3::types::PyComplex>> {
        let value=symox::evalf::eval_c64(&self.value.as_expr()).ok_or_else(||PyTypeError::new_err("coefficient has no complex numerical value"))?;
        Ok(pyo3::types::PyComplex::from_doubles(py,value.re,value.im))
    }
    fn __add__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '+', false)
    }
    fn __radd__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '+', true)
    }
    fn __sub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '-', false)
    }
    fn __rsub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '-', true)
    }
    fn __mul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '*', false)
    }
    fn __rmul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '*', true)
    }
    fn __truediv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '/', false)
    }
    fn __rtruediv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.operation(o, '/', true)
    }
    fn div(&self, o: &Bound<'_, PyAny>) -> PyResult<(Self, Self)> {
        self.divide(o, false)
    }
    fn quo(&self, o: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(self.divide(o, false)?.0)
    }
    fn rem(&self, o: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(self.divide(o, false)?.1)
    }
    fn exquo(&self, o: &Bound<'_, PyAny>) -> PyResult<Self> {
        let (q, r) = self.divide(o, false)?;
        if !r.value.is_zero() {
            return Err(ExactQuotientFailed::new_err(
                "division is not exact in the selected ring",
            ));
        }
        Ok(q)
    }
    fn __floordiv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.ring_operator(o, 'q', false)
    }
    fn __rfloordiv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.ring_operator(o, 'q', true)
    }
    fn __mod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.ring_operator(o, 'r', false)
    }
    fn __rmod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.ring_operator(o, 'r', true)
    }
    fn __divmod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.ring_operator(o, 'd', false)
    }
    fn __rdivmod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.ring_operator(o, 'd', true)
    }
    fn __neg__(&self) -> PyResult<Self> {
        let zero = self.value.domain.convert(&Expr::zero()).map_err(err)?;
        zero.binary(&self.value, '-')
            .map(|value| Self { value })
            .map_err(err)
    }
    fn __pos__(&self) -> Self {
        self.clone()
    }
    fn __pow__(&self, exponent: i64, modulo: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        if modulo.is_some() {
            return Err(PyTypeError::new_err(
                "modular power is not defined for domain elements",
            ));
        }
        self.value
            .pow(exponent)
            .map(|value| Self { value })
            .map_err(element_error)
    }
    pub(crate) fn __richcmp__(
        &self,
        o: &Bound<'_, PyAny>,
        op: CompareOp,
        py: Python<'_>,
    ) -> PyResult<PyObject> {
        if matches!(self.value.value,ElementValue::Inexact(_)) {
            if o.downcast::<super::algebraic_number::PyAlgebraicNumber>().is_ok() {
                return Ok(match op {CompareOp::Eq=>false.into_pyobject(py)?.to_owned().into_any().unbind(),CompareOp::Ne=>true.into_pyobject(py)?.to_owned().into_any().unbind(),_=>py.NotImplemented()});
            }
            let other=if let Ok(value)=o.extract::<PyRef<'_,Self>>() {
                if !matches!(value.value.value,ElementValue::Inexact(_)|ElementValue::Rational(_)) {None} else {Some(value.value.as_expr())}
            } else {to_expr(o).ok()};
            if let Some((a,b))=exact_numeric_parts(&self.value.as_expr()).zip(other.as_ref().and_then(exact_numeric_parts)) {
                let real=a.0.cmp_value(&b.0);
                let equal=real==Some(std::cmp::Ordering::Equal)&&a.1.cmp_value(&b.1)==Some(std::cmp::Ordering::Equal);
                let answer=match op {
                    CompareOp::Eq=>equal, CompareOp::Ne=>!equal,
                    _=>{
                        if self.value.domain.is_complex()||!a.1.is_zero()||!b.1.is_zero() {return Ok(py.NotImplemented());}
                        let Some(ordering)=real else{return Ok(py.NotImplemented());};
                        match op {CompareOp::Lt=>ordering.is_lt(),CompareOp::Le=>ordering.is_le(),CompareOp::Gt=>ordering.is_gt(),CompareOp::Ge=>ordering.is_ge(),_=>unreachable!()}
                    }
                };
                return Ok(answer.into_pyobject(py)?.to_owned().into_any().unbind());
            }
            return Ok(match op {CompareOp::Eq=>false.into_pyobject(py)?.to_owned().into_any().unbind(),CompareOp::Ne=>true.into_pyobject(py)?.to_owned().into_any().unbind(),_=>py.NotImplemented()});
        }
        if let Ok(other)=o.extract::<PyRef<'_,Self>>() {
            if matches!(self.value.value,ElementValue::Rational(_))&&matches!(other.value.value,ElementValue::Inexact(_)) {
                return other.__richcmp__(&self.clone().into_pyobject(py)?.into_any(),match op {CompareOp::Lt=>CompareOp::Gt,CompareOp::Le=>CompareOp::Ge,CompareOp::Gt=>CompareOp::Lt,CompareOp::Ge=>CompareOp::Le,other=>other},py);
            }
        }
        let d = PyDomain {
            inner: self.value.domain.clone(),
        };
        let b = if let Ok(b) = o.extract::<PyRef<'_, Self>>() {
            if !self.value.compatible_embedding(&b.value)
                && !(matches!(b.value.value, ElementValue::Rational(_))
                    && matches!(self.value.value, ElementValue::Rational(_)))
            {
                return Ok(match op {
                    CompareOp::Eq => false.into_pyobject(py)?.to_owned().into_any().unbind(),
                    CompareOp::Ne => true.into_pyobject(py)?.to_owned().into_any().unbind(),
                    _ => py.NotImplemented(),
                });
            }
            b.clone()
        } else {
            if !matches!(self.value.value, ElementValue::Rational(_)) {
                return Ok(py.NotImplemented());
            }
            match d.element(o, None) {
                Ok(b) => b,
                Err(_) => return Ok(py.NotImplemented()),
            }
        };
        let equal = match (&self.value.value, &b.value.value) {
            (ElementValue::Rational(a), ElementValue::Rational(b)) => a == b,
            _ => self.value == b.value,
        };
        let comparison = match op {
            CompareOp::Eq => equal,
            CompareOp::Ne => !equal,
            _ => {
                let (ElementValue::Rational(a), ElementValue::Rational(b)) =
                    (&self.value.value, &b.value.value)
                else {
                    return Ok(py.NotImplemented());
                };
                let delta = a.sub(b);
                match op {
                    CompareOp::Lt => delta.is_negative(),
                    CompareOp::Le => delta.is_negative() || delta.is_zero(),
                    CompareOp::Gt => !delta.is_negative() && !delta.is_zero(),
                    CompareOp::Ge => !delta.is_negative(),
                    _ => unreachable!(),
                }
            }
        };
        Ok(comparison.into_pyobject(py)?.to_owned().into_any().unbind())
    }
    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        Ok((
            py.import("symox._core")?
                .getattr("_restore_domain_element")?
                .unbind(),
            (self.parent(), self.as_expr())
                .into_pyobject(py)?
                .into_any()
                .unbind(),
        ))
    }
}
#[pyfunction]
#[pyo3(signature=(modulus,symmetric=true))]
fn GF(modulus: u64, symmetric: bool) -> PyResult<PyDomain> {
    if !symmetric {
        return Err(PyNotImplementedError::new_err(
            "nonsymmetric finite-field displays are not implemented",
        ));
    }
    ExactDomain::finite(modulus)
        .map(|inner| PyDomain { inner })
        .map_err(domain_error)
}
#[pyfunction]
#[pyo3(signature=(prec=None,dps=None,complex=false))]
fn _numeric_field(prec:Option<u32>,dps:Option<i64>,complex:bool)->PyResult<PyDomain> {
    if prec.is_some()&&dps.is_some(){return Err(PyTypeError::new_err("specify precision in bits or decimal digits, not both"));}
    let precision=if let Some(dps)=dps {crate::decimal_precision(dps)?} else {prec.unwrap_or(53)};
    let precision=crate::bit_precision(precision)?;
    Ok(PyDomain{inner:if precision==53&&!complex {ExactDomain::RR} else {ExactDomain::Inexact{precision,complex}}})
}
#[pyfunction]
fn _restore_domain(
    py: Python<'_>,
    kind: &str,
    prime: Option<u64>,
    parameters: Vec<PyRef<'_, crate::PyExprNode>>,
    generator: Option<PyRef<'_, crate::PyExprNode>>,
    alias: Option<PyRef<'_, crate::PyExprNode>>,
    integral: bool,
    ring: bool,
    order: &str,
) -> PyResult<PyObject> {
    let inner = match kind {
        "RR" => ExactDomain::RR,
        "ZZ" => ExactDomain::ZZ,
        "QQ" => ExactDomain::QQ,
        "EX" => ExactDomain::Expression {raw:false},
        "EXRAW" => ExactDomain::Expression {raw:true},
        "ZZ_I" => ExactDomain::Gaussian{integral:true},
        "QQ_I" => ExactDomain::Gaussian{integral:false},
        "Real"|"Complex" => ExactDomain::Inexact{precision:prime.and_then(|p|u32::try_from(p).ok()).filter(|p|(1..=crate::MAX_REQUESTED_PRECISION).contains(p)).ok_or_else(||domain_error("invalid numerical precision"))?,complex:kind=="Complex"},
        "GF" => ExactDomain::finite(prime.ok_or_else(|| domain_error("missing characteristic"))?)
            .map_err(domain_error)?,
        "RF" => ExactDomain::functions(
            integral,
            ring,
            parameters.iter().map(|p| p.e.clone()).collect(),
            order,
        )
        .map_err(domain_error)?,
        "AN" => {
            let mut field = (*NumberField::from_generator(
                &generator
                    .ok_or_else(|| domain_error("missing algebraic generator"))?
                    .e,
            )
            .map_err(err)?)
            .clone();
            field.alias = alias.map(|a| a.e.clone());
            ExactDomain::Algebraic(Arc::new(field))
        }
        _ => return Err(domain_error("unknown domain descriptor")),
    };
    PyDomain { inner }.object(py)
}
#[pyfunction]
fn _restore_domain_element(
    domain: PyRef<'_, PyDomain>,
    value: &Bound<'_, PyAny>,
) -> PyResult<PyDomainElement> {
    domain.element(value, None)
}
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDomain>()?;
    m.add_class::<PyDomainElement>()?;
    m.add(
        "ZZ",
        PyDomain {
            inner: ExactDomain::ZZ,
        },
    )?;
    m.add(
        "QQ",
        PyDomain {
            inner: ExactDomain::QQ,
        },
    )?;
    m.add_function(wrap_pyfunction!(GF, m)?)?;
    m.add("EX",PyDomain{inner:ExactDomain::Expression{raw:false}})?;
    m.add("EXRAW",PyDomain{inner:ExactDomain::Expression{raw:true}})?;
    m.add("RR",PyDomain{inner:ExactDomain::RR})?;
    m.add("CC",PyDomain{inner:ExactDomain::Inexact{precision:53,complex:true}})?;
    m.add("ZZ_I",PyDomain{inner:ExactDomain::Gaussian{integral:true}})?;
    m.add("QQ_I",PyDomain{inner:ExactDomain::Gaussian{integral:false}})?;
    m.add_function(wrap_pyfunction!(_numeric_field,m)?)?;
    m.add_function(wrap_pyfunction!(_restore_domain, m)?)?;
    m.add_function(wrap_pyfunction!(_restore_domain_element, m)?)?;
    m.add(
        "ExactQuotientFailed",
        m.py().get_type::<ExactQuotientFailed>(),
    )?;
    m.add("DomainError", m.py().get_type::<DomainError>())?;
    m.add("CoercionFailed", m.py().get_type::<CoercionFailed>())?;
    Ok(())
}

/// Whether a `domain` option is the reference EX (the domain object or the
/// string "EX").
pub fn selects_expression(domain: Option<&Bound<'_, PyAny>>) -> bool {
    let Some(domain) = domain.filter(|d| !d.is_none()) else { return false };
    if let Ok(d) = domain.extract::<PyRef<'_, PyDomain>>() {
        return matches!(d.inner, ExactDomain::Expression { raw: false });
    }
    domain.extract::<String>().is_ok_and(|s| s == "EX")
}

/// The reference option postprocessing: an EX domain needs generators.
pub fn require_expression_generators(domain: Option<&Bound<'_, PyAny>>, generators: bool) -> PyResult<()> {
    if !generators && selects_expression(domain) {
        return Err(Python::with_gil(|py| crate::poly::reference_error(py, "GeneratorsError", "you have to provide generators because EX domain was requested")));
    }
    Ok(())
}

/// The ground, generators, ring flag and order of a polynomial ring or
/// fraction field descriptor.
fn composite_parts(d: &ExactDomain) -> Option<(ExactDomain, Vec<Expr>, bool, String)> {
    match d {
        ExactDomain::Functions { integral, ring, parameters, order } => Some((if *integral { ExactDomain::ZZ } else { ExactDomain::QQ }, parameters.clone(), *ring, order.clone())),
        ExactDomain::Composite { ground, ring, parameters, order } => Some(((**ground).clone(), parameters.clone(), *ring, order.clone())),
        _ => None,
    }
}

/// The reference `Domain.unify` for composite domains: the grounds are
/// unified, the generators merged (receiver's order first), and a fraction
/// field met with a ring over a non-field ground keeps the ground's ring.
fn unify_composite(a: &ExactDomain, b: &ExactDomain) -> PyResult<ExactDomain> {
    let pa = composite_parts(a);
    let pb = composite_parts(b);
    let (ga, sa) = pa.as_ref().map_or((a.clone(), vec![]), |(g, s, _, _)| (g.clone(), s.clone()));
    let (gb, sb) = pb.as_ref().map_or((b.clone(), vec![]), |(g, s, _, _)| (g.clone(), s.clone()));
    let mut domain = unify_native(&ga, &gb)?;
    let symbols = symox::poly::polyobj::unify_gens(&sa, &sb);
    let order = pa.as_ref().or(pb.as_ref()).map(|(_, _, _, o)| o.clone()).unwrap_or_else(|| "lex".into());
    let fa = pa.as_ref().map(|p| !p.2);
    let fb = pb.as_ref().map(|p| !p.2);
    let mixed = (fa == Some(true) && fb == Some(false)) || (fb == Some(true) && fa == Some(false));
    if mixed && (!ga.is_field() || !gb.is_field()) && domain.is_field() && !matches!(domain, ExactDomain::Algebraic(_)) {
        if let Ok(ring) = domain.ring() {
            domain = ring;
        }
    }
    let ring = match (&pa, &pb) {
        (Some((_, _, ra, _)), Some((_, _, rb, _))) => if !ra || *rb { *ra } else { *rb },
        (Some((_, _, ra, _)), None) => *ra,
        (None, Some((_, _, rb, _))) => *rb,
        (None, None) => unreachable!("composite unification needs a composite operand"),
    };
    ExactDomain::composite(&domain, ring, symbols, &order).map_err(composite_error)
}

pub fn unify_native(a: &ExactDomain, b: &ExactDomain) -> PyResult<ExactDomain> {
    if a == b {
        return Ok(a.clone());
    }
    let composite_ground = |d: &ExactDomain| matches!(d, ExactDomain::Composite { .. } | ExactDomain::Gaussian { .. } | ExactDomain::Algebraic(_));
    if (matches!(a, ExactDomain::Composite { .. }) || matches!(b, ExactDomain::Composite { .. }))
        || (matches!(a, ExactDomain::Functions { .. }) && composite_ground(b))
        || (matches!(b, ExactDomain::Functions { .. }) && composite_ground(a))
    {
        if !matches!(a, ExactDomain::Expression { .. }) && !matches!(b, ExactDomain::Expression { .. }) && a.inexact_precision().is_none() && b.inexact_precision().is_none() {
            return unify_composite(a, b);
        }
    }
    let inner = match (a, b) {
        (ExactDomain::Expression {raw:a},ExactDomain::Expression {raw:b}) => ExactDomain::Expression {raw:*a || *b},
        (ExactDomain::Expression {raw}, other) | (other, ExactDomain::Expression {raw}) if other.characteristic()==0 => ExactDomain::Expression {raw:*raw},
        (a,b) if a.inexact_precision().is_some()||b.inexact_precision().is_some() => {
            if a.characteristic()!=0||b.characteristic()!=0||matches!(a,ExactDomain::Functions{..})||matches!(b,ExactDomain::Functions{..}) {
                return Err(domain_error("cannot unify numerical fields with finite or symbolic coefficient domains"));
            }
            let precision=a.inexact_precision().unwrap_or(0).max(b.inexact_precision().unwrap_or(0));
            let complex=a.is_complex()||b.is_complex()||matches!(a,ExactDomain::Gaussian{..})||matches!(b,ExactDomain::Gaussian{..})||matches!(a,ExactDomain::Algebraic(k) if symox::assumptions::is_real(&k.generator)!=Some(true))||matches!(b,ExactDomain::Algebraic(k) if symox::assumptions::is_real(&k.generator)!=Some(true));
            if precision==53&&!complex {ExactDomain::RR} else {ExactDomain::Inexact{precision,complex}}
        },
        (ExactDomain::Gaussian{integral:a},ExactDomain::Gaussian{integral:b})=>ExactDomain::Gaussian{integral:*a && *b},
        (ExactDomain::Gaussian{integral},ExactDomain::ZZ|ExactDomain::QQ)|(ExactDomain::ZZ|ExactDomain::QQ,ExactDomain::Gaussian{integral}) => ExactDomain::Gaussian{integral:*integral&&!matches!(a,ExactDomain::QQ)&&!matches!(b,ExactDomain::QQ)},
        (ExactDomain::ZZ, ExactDomain::QQ) | (ExactDomain::QQ, ExactDomain::ZZ) => ExactDomain::QQ,
        (ExactDomain::Algebraic(k), ExactDomain::ZZ | ExactDomain::QQ)
        | (ExactDomain::ZZ | ExactDomain::QQ, ExactDomain::Algebraic(k)) => {
            ExactDomain::Algebraic(k.clone())
        }
        (ExactDomain::Algebraic(k), ExactDomain::Gaussian { .. }) | (ExactDomain::Gaussian { .. }, ExactDomain::Algebraic(k)) => {
            let generators = if matches!(a, ExactDomain::Algebraic(_)) { [k.generator.clone(), Expr::i()] } else { [Expr::i(), k.generator.clone()] };
            ExactDomain::Algebraic(NumberField::from_generators(&generators).map_err(err)?)
        }
        (ExactDomain::Algebraic(a), ExactDomain::Algebraic(b)) => {
            ExactDomain::Algebraic(if a.same_embedding(b) {
                a.clone()
            } else {
                NumberField::from_generators(&[a.generator.clone(), b.generator.clone()])
                    .map_err(err)?
            })
        }
        (
            ExactDomain::Functions {
                integral: a,
                ring: ra,
                parameters: pa,
                order: oa,
            },
            ExactDomain::Functions {
                integral: b,
                ring: rb,
                parameters: pb,
                order: ob,
            },
        ) => {
            if oa != ob {
                return Err(domain_error(
                    "cannot unify coefficient rings with different monomial orders",
                ));
            }
            // SymPy's `Domain.unify`: the receiver's generator order is kept,
            // and a fraction field unified with a ring over a non-field ground
            // takes the ring of the unified ground (`QQ(a)` and `ZZ[b]` give
            // `ZZ(a,b)`).
            let parameters = symox::poly::polyobj::unify_gens(pa, pb);
            let integral = (*a && *b) || (ra != rb && (*a || *b));
            ExactDomain::functions(integral, *ra && *rb, parameters, oa).map_err(err)?
        }
        (
            ExactDomain::Functions {
                integral,
                ring,
                parameters,
                order,
            },
            ExactDomain::ZZ | ExactDomain::QQ,
        )
        | (
            ExactDomain::ZZ | ExactDomain::QQ,
            ExactDomain::Functions {
                integral,
                ring,
                parameters,
                order,
            },
        ) => ExactDomain::functions(
            *integral && !matches!(a, ExactDomain::QQ) && !matches!(b, ExactDomain::QQ),
            *ring,
            parameters.clone(),
            order,
        )
        .map_err(err)?,
        _ => {
            return Err(domain_error(
                "coefficient domains have no implemented common exact domain",
            ))
        }
    };
    Ok(inner)
}

/// Move exact implicit symbolic coefficients into the actual coefficient-ring
/// backend. A descriptor never labels an unrelated QQ/ZZ polynomial as a ring.
pub fn normalize_poly(p: PolyValue) -> PyResult<PolyValue> {
    if p.field.is_some() || p.algebraic.is_some() || p.expression.is_some() || p.base.domain_gens().is_empty() {
        return Ok(p);
    }
    if p.modulus.is_some() || p.domain == Domain::RR {
        return Ok(p);
    }
    let parameters = p.base.domain_gens().to_vec();
    if parameters.iter().any(|g| !g.is_symbol()) {
        return Ok(p);
    }
    let descriptor =
        ExactDomain::functions(p.domain == Domain::ZZ, true, parameters.clone(), "lex")
            .map_err(err)?;
    let mut field = FractionPoly::from_expr(
        &p.as_expr(),
        p.true_gens(),
        &parameters,
        &descriptor.to_string(),
    )
    .map_err(err)?;
    field.coefficient_domain = descriptor;
    Ok(field.into())
}

/// Legacy strings are parsed only at the public boundary. Native domains carry
/// their generators and embeddings directly into polynomial constructors.
pub struct Selection {
    pub native: Option<PyDomain>,
    pub legacy: Option<String>,
    pub modulus: Option<u64>,
}
impl Selection {
    pub fn from_native(inner: ExactDomain) -> Self {
        let modulus = if let ExactDomain::GF(p) = inner {
            Some(p)
        } else {
            None
        };
        let legacy = if modulus.is_some() {
            None
        } else {
            Some(inner.to_string())
        };
        Self {
            native: Some(PyDomain { inner }),
            legacy,
            modulus,
        }
    }
    pub fn new(
        domain: Option<&Bound<'_, PyAny>>,
        modulus: Option<u64>,
        extension: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let domain = domain.filter(|d| !d.is_none());
        if domain.is_some() && modulus.is_some() {
            return Err(PyValueError::new_err("specify either domain or modulus"));
        }
        if domain.is_some() && extension.is_some_and(|e| !e.is_none()) {
            return Err(PyValueError::new_err("specify either domain or extension"));
        }
        let Some(domain) = domain else {
            return Ok(Self {
                native: None,
                legacy: None,
                modulus,
            });
        };
        if let Ok(d) = domain.extract::<PyRef<'_, PyDomain>>() {
            return Ok(Self {
                native: Some(d.clone()),
                legacy: if matches!(d.inner, ExactDomain::GF(_)) {
                    None
                } else {
                    Some(d.inner.to_string())
                },
                modulus: if let ExactDomain::GF(p) = d.inner {
                    Some(p)
                } else {
                    None
                },
            });
        }
        let text = domain.extract::<String>().map_err(|_| {
            PyTypeError::new_err("domain must be a native Domain or a domain string")
        })?;
        // The reference reads "EX" as the domain object; "EXRAW" is not a
        // domain specification there.
        if text == "EX" {
            return Ok(Self::from_native(ExactDomain::Expression { raw: false }));
        }
        if text == "ZZ_I" || text == "QQ_I" {
            return Ok(Self::from_native(ExactDomain::Gaussian {
                integral: text == "ZZ_I",
            }));
        }
        // The reference `ZZ_I[x,y]` / `QQ_I[x,y]` domain strings.
        for (prefix, integral) in [("ZZ_I[", true), ("QQ_I[", false)] {
            if let Some(body) = text.strip_prefix(prefix).and_then(|t| t.strip_suffix(']')) {
                let gens = body.split(',').map(|g| symox::parser::parse_str(g.trim()).map_err(|e| PyValueError::new_err(e.to_string()))).collect::<PyResult<Vec<_>>>()?;
                let inner = ExactDomain::composite(&ExactDomain::Gaussian { integral }, true, gens, "lex").map_err(composite_error)?;
                return Ok(Self::from_native(inner));
            }
        }
        if text.starts_with("GF(") && text.ends_with(')') {
            let prime = text[3..text.len() - 1]
                .parse::<u64>()
                .map_err(domain_error)?;
            let inner = ExactDomain::finite(prime).map_err(domain_error)?;
            return Ok(Self {
                native: Some(PyDomain { inner }),
                legacy: None,
                modulus: Some(prime),
            });
        }
        Ok(Self {
            native: None,
            legacy: Some(text),
            modulus: None,
        })
    }
    pub fn algebraic(
        &self,
        extension: Option<&Bound<'_, PyAny>>,
        expressions: &[Expr],
    ) -> PyResult<Option<Arc<NumberField>>> {
        self.algebraic_in(extension, expressions, &[])
    }
    /// The selected or inferred field; `extension=True` reads the
    /// coefficients in `gens` when they are known.
    pub fn algebraic_in(
        &self,
        extension: Option<&Bound<'_, PyAny>>,
        expressions: &[Expr],
        gens: &[Expr],
    ) -> PyResult<Option<Arc<NumberField>>> {
        if let Some(PyDomain {
            inner: ExactDomain::Algebraic(field),
        }) = &self.native
        {
            return Ok(Some(field.clone()));
        }
        if self.gaussian().is_some() {
            return Ok(Some(symox::poly::gaussian::field()));
        }
        super::algebraic_extension(extension, expressions, gens, self.legacy.as_deref())
    }
    /// `Some(true)` for a selected `ZZ_I`, `Some(false)` for `QQ_I`.
    pub fn gaussian(&self) -> Option<bool> {
        match &self.native {
            Some(PyDomain {
                inner: ExactDomain::Gaussian { integral },
            }) => Some(*integral),
            _ => None,
        }
    }
    pub fn parameters(&self, available: &[Expr]) -> PyResult<Option<Vec<Expr>>> {
        if let Some(PyDomain {
            inner: ExactDomain::Functions { parameters, .. },
        }) = &self.native
        {
            return Ok(Some(parameters.clone()));
        }
        self.legacy
            .as_deref()
            .map(|d| symox::poly::fraction::parse_domain(d, available))
            .transpose()
            .map_err(err)
            .map(Option::flatten)
    }
    pub fn apply_fraction(&self, p: &mut FractionPoly) {
        if let Some(PyDomain {
            inner: domain @ ExactDomain::Functions { .. },
        }) = &self.native
        {
            p.coefficient_domain = domain.clone();
        }
    }
}

/// `Ring.is_unit` through nested polynomial rings: a ground element whose
/// innermost ring value is one or minus one (`Ring.revert`), or a nonzero
/// element of an innermost field.
fn ring_unit(domain: &ExactDomain, e: &Expr) -> bool {
    let ground = |parameters: &[Expr]| !symox::traverse::free_symbols(e).iter().any(|s| parameters.contains(s));
    match domain {
        ExactDomain::Composite { ground: inner, ring: true, parameters, .. } => ground(parameters) && ring_unit(inner, e),
        ExactDomain::Functions { ring: true, integral, parameters, .. } => {
            ground(parameters) && (!*integral && !e.is_zero() || e.is_one() || symox::ops::neg(e).is_one())
        }
        d if d.is_field() => !e.is_zero(),
        _ => e.is_one() || symox::ops::neg(e).is_one(),
    }
}
