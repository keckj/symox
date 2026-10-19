//! Shared, typed front end for ordinary division and pseudo-division.
use super::*;
use pyo3::create_exception;
use symox::poly::domains::ExactDomain;
use symox::poly::groebner::MonomialOrder;

create_exception!(symox, OptionError, PyValueError);
create_exception!(symox, UnificationFailed, PyValueError);
create_exception!(symox, FlagError, OptionError);

pub enum Operation {
    Div,
    Rem,
    Quo,
    Exquo,
    Pdiv,
    Prem,
    Pquo,
    Pexquo,
    Resultant(bool),
    Discriminant,
    Subresultants,
}
impl Operation {
    fn allow_auto(&self) -> bool {
        matches!(self, Self::Div | Self::Rem | Self::Quo | Self::Exquo)
    }
}
struct Options<'py> {
    auto: bool,
    polys: Option<bool>,
    field: Option<bool>,
    domain: Option<Bound<'py, PyAny>>,
    extension: Option<Bound<'py, PyAny>>,
    modulus: Option<u64>,
    generators: Vec<Expr>,
}
pub(super) fn boolean(value: &Bound<'_, PyAny>, name: &str) -> PyResult<Option<bool>> {
    if value.is_none() {
        return Ok(None);
    }
    if let Ok(v) = value.extract::<bool>() {
        return Ok(Some(v));
    }
    if let Ok(v) = value.extract::<i64>() {
        if v == 0 || v == 1 {
            return Ok(Some(v != 0));
        }
    }
    Err(OptionError::new_err(format!(
        "'{name}' must have a boolean value"
    )))
}
impl<'py> Options<'py> {
    fn parse(
        gens: &Bound<'py, PyTuple>,
        kw: Option<&Bound<'py, PyDict>>,
        allow_auto: bool,
    ) -> PyResult<Self> {
        let mut out = Self {
            auto: true,
            polys: None,
            field: None,
            domain: None,
            extension: None,
            modulus: None,
            generators: gens_of(gens)?,
        };
        let mut explicit_auto = false;
        if let Some(kw) = kw {
            for (key, value) in kw.iter() {
                let key = key.extract::<String>()?;
                match key.as_str() {
                    "auto" => {
                        if !allow_auto {
                            return Err(FlagError::new_err(
                                "'auto' flag is not allowed in this context",
                            ));
                        }
                        if let Some(auto) = boolean(&value, "auto")? {
                            out.auto = auto;
                            explicit_auto = true;
                        }
                    }
                    "polys" => out.polys = boolean(&value, "polys")?,
                    "field" => out.field = boolean(&value, "field")?,
                    "domain" if !value.is_none() => out.domain = Some(value),
                    "extension" if !value.is_none() => {
                        if value.is_instance_of::<pyo3::types::PyBool>()
                            && !value.extract::<bool>()?
                        {
                            return Err(OptionError::new_err(
                                "extension=False is not a valid extension option",
                            ));
                        }
                        out.extension = Some(value);
                    }
                    "modulus" if !value.is_none() => {
                        out.modulus = Some(check_mod(value.extract::<u64>()?)?)
                    }
                    "gens" if !value.is_none() => {
                        if !gens.is_empty() {
                            return Err(OptionError::new_err(
                                "both '*gens' and keyword argument 'gens' supplied",
                            ));
                        }
                        out.generators = if let Ok(t) = value.downcast::<PyTuple>() {
                            gens_of(t)?
                        } else if let Ok(l) = value.downcast::<PyList>() {
                            l.iter().map(|v| to_expr(&v)).collect::<PyResult<_>>()?
                        } else {
                            vec![to_expr(&value)?]
                        };
                    }
                    "order" if !value.is_none() => {
                        let order=value.extract::<String>().map_err(|_|PyNotImplementedError::new_err("division supports the named monomial orders lex, grlex and grevlex"))?;
                        if MonomialOrder::parse(&order).is_none() {
                            return Err(OptionError::new_err(
                                "order must be lex, grlex or grevlex",
                            ));
                        }
                        // Like SymPy's dense division, this is a validated
                        // polynomial option. It does not change recursive
                        // division into sparse monomial reduction.
                    }
                    "domain" | "extension" | "modulus" | "gens" | "order" => {}
                    "greedy" | "sort" | "wrt" | "symmetric" | "composite" | "gaussian"
                    | "split" | "strict" => {
                        return Err(PyNotImplementedError::new_err(format!(
                            "the '{key}' polynomial option is not implemented for global division"
                        )))
                    }
                    _ => {
                        return Err(OptionError::new_err(format!(
                            "'{key}' is not a valid option"
                        )))
                    }
                }
            }
        }
        // An explicit coefficient-domain or field choice disables automatic
        // extension unless auto=True is actually supplied. None is an absent
        // option, matching the reference Options preprocessing.
        if !explicit_auto && (out.domain.is_some() || out.field.is_some()
            || out.modulus.is_some()
            || out.extension.as_ref().is_some_and(|value|
                !(value.is_instance_of::<pyo3::types::PyBool>() && value.extract::<bool>().unwrap_or(false)))) {
            out.auto = false;
        }
        if out.domain.is_some() && out.field.is_some() {
            return Err(OptionError::new_err(
                "'domain' option is not allowed together with 'field'",
            ));
        }
        if out.domain.is_some() && out.modulus.is_some() {
            return Err(OptionError::new_err(
                "'domain' option is not allowed together with 'modulus'",
            ));
        }
        if out.domain.is_some() && out.extension.is_some() {
            return Err(OptionError::new_err(
                "'domain' option is not allowed together with 'extension'",
            ));
        }
        if out.modulus.is_some() && out.extension.is_some() {
            return Err(OptionError::new_err(
                "'modulus' option is not allowed together with 'extension'",
            ));
        }
        for (i, g) in out.generators.iter().enumerate() {
            if g.is_number() || out.generators[..i].contains(g) {
                return Err(OptionError::new_err(
                    "polynomial generators must be distinct nonnumeric expressions",
                ));
            }
        }
        Ok(out)
    }
}
fn compatible(a: &PolyValue, b: &PolyValue) -> PyResult<()> {
    if [a, b]
        .iter()
        .any(|p| p.modulus.is_some() && !p.base.domain_gens().is_empty())
    {
        return Err(UnificationFailed::new_err("finite-field division requires scalar coefficients; finite symbolic coefficient domains are not implemented"));
    }
    if a.modulus != b.modulus && (a.modulus.is_some() || b.modulus.is_some()) {
        return Err(UnificationFailed::new_err(
            "polynomial coefficient fields have incompatible characteristics",
        ));
    }
    if let (Ok(a), Ok(b)) = (domains::PyDomain::of(a), domains::PyDomain::of(b)) {
        // Retain the existing inexact scalar backend. It has no new exact
        // DomainElement implementation, but its polynomial arithmetic works.
        if !matches!(a.inner, ExactDomain::RR) && !matches!(b.inner, ExactDomain::RR) {
            domains::unify_native(&a.inner, &b.inner)
                .map_err(|e| UnificationFailed::new_err(e.to_string()))?;
        } else if !matches!(a.inner, ExactDomain::ZZ | ExactDomain::QQ | ExactDomain::RR)
            || !matches!(b.inner, ExactDomain::ZZ | ExactDomain::QQ | ExactDomain::RR)
        {
            return Err(UnificationFailed::new_err(
                "no implemented common domain for inexact and composite coefficient fields",
            ));
        }
    }
    if a.algebraic.is_some() || b.algebraic.is_some() {
        a.checked_algebraic_pair(b)
            .map_err(|e| UnificationFailed::new_err(e.to_string()))?;
    }
    if let Some(prototype) = a.field.as_ref().or(b.field.as_ref()) {
        let gens = symox::poly::polyobj::unify_gens(a.true_gens(), b.true_gens());
        let expressions = [a.as_expr(), b.as_expr()];
        let mut parameters = prototype.parameters.clone();
        for p in symox::poly::fraction::infer_parameters(&expressions, &gens) {
            if !parameters.contains(&p) {
                parameters.push(p);
            }
        }
        for e in &expressions {
            FractionPoly::from_expr(e, &gens, &parameters, &prototype.display)
                .map_err(|e| UnificationFailed::new_err(e.to_string()))?;
        }
    }
    Ok(())
}
fn common_domain(a: &domains::PyDomain, b: &domains::PyDomain) -> PyResult<ExactDomain> {
    if matches!(a.inner, ExactDomain::ZZ | ExactDomain::QQ | ExactDomain::RR)
        && matches!(b.inner, ExactDomain::ZZ | ExactDomain::QQ | ExactDomain::RR)
        && (matches!(a.inner, ExactDomain::RR) || matches!(b.inner, ExactDomain::RR))
    {
        return Ok(ExactDomain::RR);
    }
    domains::unify_native(&a.inner, &b.inner).map_err(|e| UnificationFailed::new_err(e.to_string()))
}
// extension=True infers algebraic constants; it does not erase a coefficient
// domain already carried by two Poly objects. Reframing different generators
// still follows joint expression inference, just as the ordinary constructor.
fn inherited_extension(
    f: &Bound<'_, PyAny>,
    g: &Bound<'_, PyAny>,
    opts: &Options<'_>,
) -> PyResult<Option<domains::PyDomain>> {
    if !opts.extension.as_ref().is_some_and(|e| {
        e.is_instance_of::<pyo3::types::PyBool>() && e.extract::<bool>().unwrap_or(false)
    }) {
        return Ok(None);
    }
    let (Ok(a), Ok(b)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) else {
        return Ok(None);
    };
    if !opts.generators.is_empty()
        && (!same_generator_set(&opts.generators, a.get().p.true_gens())
            || !same_generator_set(&opts.generators, b.get().p.true_gens()))
    {
        return Ok(None);
    }
    let (Ok(ka), Ok(kb)) = (
        domains::PyDomain::of(&a.get().p),
        domains::PyDomain::of(&b.get().p),
    ) else {
        return Ok(None);
    };
    compatible(&a.get().p, &b.get().p)?;
    Ok(Some(domains::PyDomain {
        inner: common_domain(&ka, &kb)?,
    }))
}
fn promote(p: PolyValue) -> PyResult<PolyValue> {
    let domain=domains::PyDomain::of(&p).map_err(|_|PyNotImplementedError::new_err("field=True for opaque symbolic coefficient domains is not implemented; their default polynomial arithmetic remains supported"))?;
    let field = domain.inner.field();
    if field == domain.inner {
        return Ok(p);
    }
    domains::PyDomain { inner: field }.poly(&p.as_expr(), p.true_gens())
}
// A finite coefficient retains its original residue when the target prime
// changes. Its signed expression is only the associated ZZ representative.
fn coefficient_expression(
    obj: &Bound<'_, PyAny>,
    expression: &Expr,
    target: &domains::PyDomain,
) -> PyResult<Expr> {
    let Ok(poly) = obj.downcast::<PyPoly>() else {
        return Ok(expression.clone());
    };
    let p = &poly.get().p;
    let Some(prime) = p.modulus else {
        return Ok(expression.clone());
    };
    match target.inner {
        ExactDomain::ZZ=>Ok(expression.clone()),
        ExactDomain::GF(_) => {
            let mut terms=Vec::new();
            for (powers,coefficient) in p.terms(){
                let value=coefficient.as_number().and_then(Rat::from_number).filter(|r|r.is_integer()).ok_or_else(||domains::CoercionFailed::new_err("expected scalar finite-field coefficients"))?;
                let representative=value.n.modulo(&symox::poly::Int::from_u64(prime));
                let mut factors=vec![target.inner.convert(&Expr::num(representative.to_number())).map_err(|e|domains::CoercionFailed::new_err(e.to_string()))?.as_expr()];
                for (generator,power) in p.true_gens().iter().zip(powers){if power!=0{factors.push(generator.powi(power as i64));}}
                terms.push(symox::ops::mul(&factors));
            }
            Ok(symox::ops::add(&terms))
        },
        _=>Err(domains::CoercionFailed::new_err("no canonical conversion from a finite coefficient field to this domain; convert to ZZ explicitly first")),
    }
}
fn same_generator_set(a: &[Expr], b: &[Expr]) -> bool {
    a.len() == b.len() && a.iter().all(|g| b.contains(g))
}
fn inferred_generators(
    expressions: &[Expr],
    domain: Option<&domains::PyDomain>,
) -> PyResult<Vec<Expr>> {
    let polynomials = Poly::from_exprs(expressions, &[]).map_err(err)?;
    let mut generators = polynomials[0].true_gens().to_vec();
    if let Some(domain) = domain {
        generators.retain(|g| domain.inner.convert(g).is_err());
    }
    if generators.is_empty() {
        return Err(PyValueError::new_err(
            "can't initialize polynomial division without generators",
        ));
    }
    Ok(generators)
}
struct Inputs {
    a: PolyValue,
    b: PolyValue,
    auto: bool,
    polys: bool,
}
impl Inputs {
    fn build<'py>(
        f: &Bound<'py, PyAny>,
        g: &Bound<'py, PyAny>,
        gens: &Bound<'py, PyTuple>,
        kw: Option<&Bound<'py, PyDict>>,
        operation: &Operation,
    ) -> PyResult<Self> {
        let opts = Options::parse(gens, kw, operation.allow_auto())?;
        let inherited = inherited_extension(f, g, &opts)?;
        let fp = f.downcast::<PyPoly>().ok();
        let gp = g.downcast::<PyPoly>().ok();
        let polys = opts.polys.unwrap_or(fp.is_some() || gp.is_some());
        let (mut a, mut b) = if let (Some(fp), Some(gp)) = (fp, gp) {
            if (opts.generators.is_empty()
                || (opts.generators == fp.get().p.true_gens()
                    && opts.generators == gp.get().p.true_gens()))
                && opts.domain.is_none()
                && (opts.extension.is_none() || inherited.is_some())
                && opts.modulus.is_none()
            {
                compatible(&fp.get().p, &gp.get().p)?;
                let b = fp.get().other(g).ok_or_else(|| {
                    UnificationFailed::new_err("polynomial coefficient domains cannot be unified")
                })?;
                (fp.get().p.clone(), b)
            } else {
                Self::construct(f, g, &opts)?
            }
        } else {
            Self::construct(f, g, &opts)?
        };
        compatible(&a, &b)?;
        if opts.field == Some(true) {
            a = promote(a)?;
            b = promote(b)?;
        }
        Ok(Self {
            a,
            b,
            auto: opts.auto,
            polys,
        })
    }
    fn construct<'py>(
        f: &Bound<'py, PyAny>,
        g: &Bound<'py, PyAny>,
        opts: &Options<'py>,
    ) -> PyResult<(PolyValue, PolyValue)> {
        let expressions = [to_expr(f)?, to_expr(g)?];
        let selection =
            domains::Selection::new(opts.domain.as_ref(), opts.modulus, opts.extension.as_ref())?;
        let inherited = inherited_extension(f, g, opts)?;
        let mut selected = selection.native.clone().or_else(|| inherited.clone());
        if selected.is_none() {
            selected = match selection.legacy.as_deref() {
                Some("ZZ") => Some(domains::PyDomain {
                    inner: ExactDomain::ZZ,
                }),
                Some("QQ") => Some(domains::PyDomain {
                    inner: ExactDomain::QQ,
                }),
                Some("RR") => Some(domains::PyDomain {
                    inner: ExactDomain::RR,
                }),
                _ => None,
            };
        }
        if inherited.is_none() {
            if let Some(field) = selection.algebraic(opts.extension.as_ref(), &expressions)? {
                selected = Some(domains::PyDomain {
                    inner: ExactDomain::Algebraic(field),
                });
            }
        }
        if selected.is_none() {
            if let Some(prime) = selection.modulus {
                selected = Some(domains::PyDomain {
                    inner: ExactDomain::GF(prime),
                });
            }
        }
        // A string is a legacy boundary only; resolve its actual parameters
        // against both inputs once, then use one typed descriptor thereafter.
        if selected.is_none() {
            let mut available = symox::poly::fraction::infer_parameters(&expressions, &[]);
            for v in &opts.generators {
                if v.is_symbol() && !available.contains(v) {
                    available.push(v.clone());
                }
            }
            if let Some(parameters) = selection.parameters(&available)? {
                let text = selection.legacy.as_deref().unwrap();
                selected = Some(domains::PyDomain {
                    inner: ExactDomain::functions(
                        text.starts_with("ZZ"),
                        text.ends_with(']'),
                        parameters,
                        "lex",
                    )
                    .map_err(err)?,
                });
            }
        }
        if selected.is_none()
            && opts.domain.is_none()
            && opts.extension.is_none()
            && opts.modulus.is_none()
            && !opts.generators.is_empty()
        {
            if let (Ok(a), Ok(b)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
                if same_generator_set(&opts.generators, a.get().p.true_gens())
                    && same_generator_set(&opts.generators, b.get().p.true_gens())
                {
                    compatible(&a.get().p, &b.get().p)?;
                    if let (Ok(a), Ok(b)) = (
                        domains::PyDomain::of(&a.get().p),
                        domains::PyDomain::of(&b.get().p),
                    ) {
                        selected = Some(domains::PyDomain {
                            inner: common_domain(&a, &b)?,
                        });
                    }
                }
            }
        }
        let mut generators = opts.generators.clone();
        if generators.is_empty() {
            // With two Poly arguments their generators remain explicit even
            // when a construction option requests another coefficient domain.
            if let (Ok(a), Ok(b)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
                generators =
                    symox::poly::polyobj::unify_gens(a.get().p.true_gens(), b.get().p.true_gens());
            } else {
                generators = inferred_generators(&expressions, selected.as_ref())?;
            }
        }
        if let Some(domain) = selected {
            return Ok((
                domain.poly(
                    &coefficient_expression(f, &expressions[0], &domain)?,
                    &generators,
                )?,
                domain.poly(
                    &coefficient_expression(g, &expressions[1], &domain)?,
                    &generators,
                )?,
            ));
        }
        let py = f.py();
        let args = PyTuple::new(py, generators.iter().cloned().map(w))?;
        // These constructors also retain the established opaque and inexact
        // coefficient backends. Results themselves are never re-polified.
        // Poly inputs are re-polified from their expressions, as in SymPy's
        // parallel_poly_from_expr; the Poly constructor would keep their domain.
        let fe = w(expressions[0].clone()).into_pyobject(py)?.into_any();
        let ge = w(expressions[1].clone()).into_pyobject(py)?.into_any();
        let a = PyPoly::py_new(&fe, &args, opts.domain.as_ref(), None, None, None)?.p;
        let b = PyPoly::py_new(&ge, &args, opts.domain.as_ref(), None, None, None)?.p;
        Ok((a, b))
    }
    fn render(&self, py: Python<'_>, p: PolyValue) -> PyResult<PyObject> {
        if self.polys {
            Ok(PyPoly { p }.into_pyobject(py)?.into_any().unbind())
        } else {
            Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())
        }
    }
    fn calculate(&self, py: Python<'_>, op: Operation) -> PyResult<PyObject> {
        match op {
            Operation::Div => {
                let (q, r) = self.a.div_with_options(&self.b, self.auto).map_err(err)?;
                Ok(
                    PyTuple::new(py, [self.render(py, q)?, self.render(py, r)?])?
                        .into_any()
                        .unbind(),
                )
            }
            Operation::Rem => self.render(
                py,
                self.a.rem_with_options(&self.b, self.auto).map_err(err)?,
            ),
            Operation::Quo => self.render(
                py,
                self.a.quo_with_options(&self.b, self.auto).map_err(err)?,
            ),
            Operation::Exquo => self.render(
                py,
                self.a.exquo_with_options(&self.b, self.auto).map_err(err)?,
            ),
            Operation::Pdiv | Operation::Prem | Operation::Pquo | Operation::Pexquo => {
                let (q, r) = self.a.pdiv(&self.b).map_err(err)?;
                match op {
                    Operation::Pdiv => Ok(PyTuple::new(
                        py,
                        [self.render(py, q)?, self.render(py, r)?],
                    )?
                    .into_any()
                    .unbind()),
                    Operation::Prem => self.render(py, r),
                    Operation::Pexquo if !r.is_zero() => {
                        Err(domains::ExactQuotientFailed::new_err(
                            "pseudo-division has a nonzero remainder",
                        ))
                    }
                    _ => self.render(py, q),
                }
            }
            Operation::Subresultants => {
                let values=self.a.subresultants(&self.b).map_err(err)?;
                Ok(PyList::new(py,values.into_iter().map(|p|self.render(py,p)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind())
            }
            Operation::Discriminant => render_eliminated(py,self.a.discriminant_poly().map_err(err)?,self.polys),
            Operation::Resultant(include_prs) => {
                let(result,sequence)=self.a.resultant_prs(&self.b).map_err(err)?;
                let result=render_eliminated(py,result,self.polys)?;
                if !include_prs {return Ok(result);}
                let sequence=PyList::new(py,sequence.into_iter().map(|p|self.render(py,p)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind();
                Ok(PyTuple::new(py,[result,sequence])?.into_any().unbind())
            }
        }
    }
}
pub fn render_eliminated(py:Python<'_>,p:PolyValue,polys:bool)->PyResult<PyObject> {
    if polys && p.true_gens().len()>1 {
        Ok(PyPoly{p:p.eliminate_first().map_err(err)?}.into_pyobject(py)?.into_any().unbind())
    } else {
        Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())
    }
}
pub(super) fn prepare_unary<'py>(f:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,kw:Option<&Bound<'py,PyDict>>)->PyResult<(PolyValue,bool)> {
    let input=Inputs::build(f,f,gens,kw,&Operation::Discriminant)?;
    Ok((input.a,input.polys))
}
/// Joint coefficient-domain inference for fractional factor lists. The list
/// renderer may then restrict main generators to each side, while retaining
/// this shared descriptor (including unused parameters and embeddings).
pub(super) fn prepare_pair<'py>(f:&Bound<'py,PyAny>,g:&Bound<'py,PyAny>,gens:&Bound<'py,PyTuple>,kw:Option<&Bound<'py,PyDict>>)->PyResult<(PolyValue,PolyValue,bool)> {
    let input=Inputs::build(f,g,gens,kw,&Operation::Discriminant)?;
    Ok((input.a,input.b,input.polys))
}
pub(super) fn option_generators(gens:&Bound<'_,PyTuple>,kw:Option<&Bound<'_,PyDict>>)->PyResult<Vec<Expr>> {
    Ok(Options::parse(gens,kw,false)?.generators)
}
pub(super) fn validate_unary(gens:&Bound<'_,PyTuple>,kw:Option<&Bound<'_,PyDict>>)->PyResult<()> {
    Options::parse(gens,kw,false).map(|_|())
}
pub fn run<'py>(
    f: &Bound<'py, PyAny>,
    g: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kw: Option<&Bound<'py, PyDict>>,
    op: Operation,
) -> PyResult<PyObject> {
    Inputs::build(f, g, gens, kw, &op)?.calculate(f.py(), op)
}
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("OptionError", m.py().get_type::<OptionError>())?;
    m.add("FlagError", m.py().get_type::<FlagError>())?;
    m.add("UnificationFailed", m.py().get_type::<UnificationFailed>())?;
    Ok(())
}
