//! Python bindings: a SymPy compatible API on top of the Rust core.

// SymPy's API uses `is_Add`, `as_coeff_Mul`, ...
#![allow(non_snake_case)]

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use num_bigint::BigInt;
use num_traits::{Signed,Zero,ToPrimitive};
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1};
use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict, PyFloat, PyInt, PyList, PySet, PyString, PyTuple};
use rustc_hash::{FxHashMap, FxHashSet};

mod series;
mod python_hash;
mod expression_inspection;

use symox::expr::{Constant, Expr, Func, Kind, RelOp};
use symox::lambdify::Compiled;
use symox::number::Number;
mod matrices;
mod matrixexpr;
mod matrix_dag;
mod linear_solvers;
mod printing;

use symox::{diff, evalf, expand, functions as fx, ops, parser, printer, relational, simplify, subs, traverse};

mod poly;
mod ntheory;
mod simp;
mod simptools;
mod ode;
mod solve;
mod logic;
mod special;
mod calculus;
mod ask;
mod recurrences;
mod transforms;

/// A symbolic expression.
#[pyclass(name = "Expr", module = "symox._core", frozen, weakref, subclass)]
#[derive(Clone)]
pub struct PyExprNode {
    pub e: Expr,
}

#[derive(Clone)]
pub struct PyExpr {
    pub e: Expr,
}

impl<'py> IntoPyObject<'py> for PyExpr {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        scalar_object(py, self.e).map(|value| value.into_bound(py))
    }
}

impl<'py> IntoPyObject<'py> for &PyExpr {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        scalar_object(py, self.e.clone()).map(|value| value.into_bound(py))
    }
}

impl<'py> FromPyObject<'py> for PyExpr {
    fn extract_bound(value: &Bound<'py, PyAny>) -> PyResult<Self> {
        Ok(wrap(value.downcast::<PyExprNode>()?.get().e.clone()))
    }
}

pub(crate) fn wrap(e: Expr) -> PyExpr {
    PyExpr { e }
}

pub(crate) fn expr_object(py: Python<'_>, e: Expr) -> PyResult<PyObject> {
    if let Kind::Func(Func::MatrixProxy(parent), _) = e.kind() {
        return Ok(matrixexpr::wrap_matrixexpr((**parent).clone()).into_pyobject(py)?.into_any().unbind());
    }
    scalar_object(py, e)
}

fn scalar_object(py: Python<'_>, e: Expr) -> PyResult<PyObject> {
    if matches!(e.kind(), Kind::Func(Func::Undefined(head), _) if head.is_internal("Tuple")) {
        if let Ok(class) = py.import("symox").and_then(|m| m.getattr("Tuple")) {
            if class.downcast::<pyo3::types::PyType>().is_ok() {
                let value = Bound::new(py, PyExprNode { e })?;
                return Ok(class.call_method1("_from_native", (value,))?.unbind());
            }
        }
    }
    let singleton = if e.is_rational_number() {
        if e.is_zero() { Some("Zero") } else if e.is_one() { Some("One") }
        else if e.is_minus_one() { Some("NegativeOne") }
        else if e == Expr::rational(1, 2) { Some("Half") } else { None }
    } else { None };
    if let Some(name) = singleton {
        // The registry is unavailable during its own bootstrap. Once
        // initialized it owns the canonical Python objects for these atoms.
        if let Ok(value) = py.import("symox").and_then(|m| m.getattr("S"))
            .and_then(|s| s.getattr(name)) {
            return Ok(value.unbind());
        }
    }
    Ok(Bound::new(py, PyExprNode { e })?.into_any().unbind())
}

fn arithmetic_object(slf: &PyRef<'_, crate::PyExprNode>, other: &Bound<'_, PyAny>, e: Expr) -> PyResult<PyObject> {
    if e.ptr_eq(&slf.e) {
        return Ok(slf.into_pyobject(slf.py())?.to_owned().into_any().unbind());
    }
    if other.downcast::<crate::PyExprNode>().is_ok_and(|other| e.ptr_eq(&other.get().e)) {
        return Ok(other.clone().unbind());
    }
    expr_object(slf.py(), e)
}

fn matrix_proxy(parent: symox::matrixexpr::MatrixExpr) -> Expr {
    Expr::from_kind(Kind::Func(Func::MatrixProxy(std::sync::Arc::new(parent)), vec![]))
}

fn scalar_expr(obj: &Bound<'_, PyAny>) -> PyResult<Expr> {
    let e = to_expr(obj)?;
    if matches!(e.kind(), Kind::Func(Func::MatrixProxy(_), _)) {
        return Err(PyTypeError::new_err("a matrix expression cannot be used as a scalar"));
    }
    Ok(e)
}

fn replacement_expr(obj: &Bound<'_, PyAny>) -> PyResult<Expr> {
    if obj.downcast::<matrices::PyMatrix>().is_ok() || obj.getattr("_symox_sparse").is_ok_and(|value| value.is_truthy().unwrap_or(false)) {
        return Ok(matrix_proxy(matrixexpr::to_matrixexpr(obj)?));
    }
    to_expr(obj)
}

/// Matrix reconstruction is fallible (dimensions, indices, scalar/matrix
/// types). Validate it before entering the scalar core's infallible traversal.
fn validate_matrix_replacements(py: Python<'_>, e: &Expr, pairs: &[(Expr, Expr)], exact: bool) -> PyResult<()> {
    if !e.has_matrix() { return Ok(()); }
    let mut matrix_pairs = Vec::new();
    let mut scalars = FxHashMap::default();
    for (old, new) in pairs {
        if let Kind::Func(Func::MatrixProxy(old), _) = old.kind() {
            if let Kind::Func(Func::MatrixProxy(new), _) = new.kind() {
                matrix_pairs.push(((**old).clone(), (**new).clone()));
            } else { return Err(PyTypeError::new_err("Matrix replacements must be matrices")); }
        } else { scalars.insert(old.clone(), new.clone()); }
    }
    let mut transform = |value: &Expr| if exact { subs::xreplace(value, &scalars) } else { subs::subs(value, pairs) };
    fn validate(py: Python<'_>, node: &Expr, pairs: &[(Expr, Expr)], matrix_pairs: &[(symox::matrixexpr::MatrixExpr, symox::matrixexpr::MatrixExpr)], transform: &mut dyn FnMut(&Expr)->Expr, visited: &mut FxHashMap<usize,Expr>) -> PyResult<()> {
        if !node.has_matrix() || pairs.iter().any(|(old, _)| old == node) { return Ok(()); }
        if visited.insert(node.dag_id(),node.clone()).is_some(){return Ok(());}
        let parents=match node.kind(){Kind::Func(Func::MatrixDotProduct(a,b),_)=>vec![a,b],Kind::Func(Func::MatrixElement(a)|Func::MatrixTrace(a)|Func::MatrixDeterminant(a)|Func::MatrixProxy(a),_)=>vec![a],_=>Vec::new()};
        for parent in parents{
            let mut failure=None;
            parent.scalar_map(&mut |field|{
                if failure.is_none(){if let Err(error)=validate(py,field,pairs,matrix_pairs,transform,visited){failure=Some(error);}}
                field.clone()
            }).map_err(|e|matrices::err(py,e))?;
            if let Some(error)=failure{return Err(error);}
        }
        if let Kind::Func(Func::MatrixDotProduct(a,b),_)=node.kind(){
            let a=a.map_all(matrix_pairs,transform).map_err(|e|matrices::err(py,e))?;
            let b=b.map_all(matrix_pairs,transform).map_err(|e|matrices::err(py,e))?;
            symox::matrixcontraction::validate(&a,&b).map_err(|e|matrices::err(py,e))?;
        }
        if let Kind::Func(f @ (Func::MatrixElement(parent) | Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent)), indices) = node.kind() {
            parent.map_all(matrix_pairs, transform).and_then(|matrix| {
                if matches!(f, Func::MatrixTrace(_)) { matrix.trace() }
                else if matches!(f, Func::MatrixDeterminant(_)) { matrix.determinant() }
                else { matrix.element(&transform(&indices[0]), &transform(&indices[1])) }
            }).map_err(|error| matrices::err(py, error))?;
        }
        for child in node.args() { validate(py, &child, pairs, matrix_pairs, transform, visited)?; }
        Ok(())
    }
    validate(py, e, pairs, &matrix_pairs, &mut transform, &mut FxHashMap::default())
}

pub(crate) fn complex_parts_options(hints:Option<&Bound<'_,PyDict>>,scalar:bool)->PyResult<Option<Expr>>{
    let mut ignore=None;
    if let Some(hints)=hints{
        for(key,value)in hints.iter(){
            let key:String=key.extract()?;
            if key!="ignore"{return Err(PyNotImplementedError::new_err(format!("as_real_imag hint {key} is not implemented")));}
            if !value.is_none(){
                if !scalar{return Err(PyNotImplementedError::new_err("matrix as_real_imag(ignore=...) is not implemented"));}
                ignore=Some(to_expr(&value)?);
            }
        }
    }
    Ok(ignore)
}

#[pyfunction]
#[pyo3(signature=(values,deep=true))]
fn _as_real_imag_many(values:Vec<PyExpr>,deep:bool)->PyResult<(Vec<PyExpr>,Vec<PyExpr>)>{
    let entries=values.into_iter().map(|entry|entry.e).collect::<Vec<_>>();
    let parts=symox::complex_parts::as_real_imag_many(&entries,deep)
        .map_err(|error|PyNotImplementedError::new_err(error.to_string()))?;
    Ok(parts.into_iter().map(|(real,imag)|(wrap(real),wrap(imag))).unzip())
}

fn function_spec(name:&str,assumptions:&Bound<'_,PyDict>,nargs:Option<Vec<usize>>)->PyResult<symox::expr::UndefinedFunction>{
    let mut kwargs=Vec::new();
    for (key,value) in assumptions.iter(){
        let key:String=key.extract()?;
        if value.is_none(){
            if key=="commutative"{return Err(PyValueError::new_err("commutativity must be True or False"));}
            if !symox::assumption_facts::FactSet::NAMES.contains(&key.as_str()){
                return Err(PyValueError::new_err(format!("unknown assumption {key}")));
            }
        }else{
            if key=="commutative"&&!value.eq(true)?&&!value.eq(false)?{
                return Err(PyValueError::new_err("commutativity must be True or False"));
            }
            kwargs.push((key,value.is_truthy()?));
        }
    }
    let given=symox::assumption_facts::FactSet::from_kwargs(kwargs.iter().map(|(k,v)|(k.as_str(),*v))).map_err(PyValueError::new_err)?;
    symox::expr::UndefinedFunction::with_facts(name,given,nargs).map_err(PyValueError::new_err)
}

fn function_kwargs<'py>(py:Python<'py>,function:&symox::expr::UndefinedFunction)->PyResult<Bound<'py,PyDict>>{
    let kwargs=PyDict::new(py);
    for (key,value) in function.head_facts().iter(){kwargs.set_item(key,value)?;}
    if let Some(counts)=&function.nargs{kwargs.set_item("nargs",PyTuple::new(py,counts)?)?;}
    Ok(kwargs)
}

/// Canonical identity metadata, including all standard true and false facts.
#[pyfunction]
fn _function_facts<'py>(py:Python<'py>,assumptions:&Bound<'py,PyDict>)->PyResult<Bound<'py,PyDict>>{
    function_kwargs(py,&function_spec("",assumptions,None)?)
}

/// Read a SymPy head without importing SymPy, or converting metadata by text.
fn sympy_function_spec(function:&Bound<'_,PyAny>)->PyResult<symox::expr::UndefinedFunction>{
    let name:String=function.getattr("__name__")?.extract()?;
    let options=function.getattr("_kwargs")?;
    let options=options.downcast::<PyDict>()?;
    let selected=PyDict::new(function.py());
    for (key,value) in options.iter(){
        if key.extract::<String>()?!="nargs"{selected.set_item(key,value)?;}
    }
    let counts=function.getattr("nargs")?;
    let nargs=if counts.getattr("is_FiniteSet").is_ok_and(|value|value.is_truthy().unwrap_or(false)){
        Some(counts.try_iter()?.map(|value|to_expr(&value?).and_then(|e|e.as_i64().and_then(|n|usize::try_from(n).ok()).ok_or_else(||PyValueError::new_err("invalid function nargs")))).collect::<PyResult<Vec<_>>>()?)
    }else{None};
    function_spec(&name,&selected,nargs)
}

pub(crate) fn function_head(value:&Bound<'_,PyAny>)->PyResult<Option<symox::expr::UndefinedFunction>>{
    if let (Ok(name),Ok(assumptions),Ok(nargs))=(value.getattr("name"),value.getattr("_assumptions"),value.getattr("_nargs")){
        if let (Ok(name),Ok(assumptions))=(name.extract::<String>(),assumptions.downcast::<PyDict>()){
            return function_spec(&name,assumptions,nargs.extract()?).map(Some);
        }
    }
    let sympy_head=value.get_type().getattr("__mro__").is_ok_and(|mro|mro.try_iter().is_ok_and(|bases|bases.filter_map(Result::ok).any(|base|{
        base.getattr("__name__").is_ok_and(|name|name.extract::<String>().is_ok_and(|name|name=="UndefinedFunction"))
            &&base.getattr("__module__").is_ok_and(|module|module.extract::<String>().is_ok_and(|module|module.starts_with("sympy.")))
    })));
    if sympy_head{sympy_function_spec(value).map(Some)}else{Ok(None)}
}

fn undefined_head<'py>(py:Python<'py>,function:&symox::expr::UndefinedFunction)->PyResult<Bound<'py,PyAny>>{
    py.import("symox")?.getattr("Function")?.call((function.name.as_ref(),),Some(&function_kwargs(py,function)?))
}

thread_local! {
    static ACTIVE_CONVERSIONS: std::cell::RefCell<Vec<usize>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn direct_conversion_hook<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Option<Bound<'py, PyAny>>> {
    if obj.downcast::<crate::PyExprNode>().is_ok() { return Ok(None); }
    if let Ok(hook) = obj.getattr("_symox_") { return Ok(Some(hook)); }
    let ty = obj.get_type();
    let module = ty.module().map(|value| value.to_string()).unwrap_or_default();
    // Native interoperability exporters must never be implicit conversion
    // hooks. Reference objects retain the established structural conversion.
    if module == "symox" || module.starts_with("symox.") { return Ok(None); }
    for base in ty.getattr("__mro__")?.try_iter()? {
        let module = base?.getattr("__module__")?.extract::<String>().unwrap_or_default();
        if module == "sympy" || module.starts_with("sympy.") { return Ok(None); }
    }
    Ok(obj.getattr("_sympy_").ok())
}

fn conversion_hook(obj: &Bound<'_, PyAny>, hook: &Bound<'_, PyAny>) -> PyResult<Expr> {
    struct RecursionGuard;
    impl Drop for RecursionGuard {
        fn drop(&mut self) {
            // The Bound input keeps this conversion inside an attached Python
            // call, including when a user hook raises an exception.
            unsafe { pyo3::ffi::Py_LeaveRecursiveCall() };
            ACTIVE_CONVERSIONS.with(|active| { active.borrow_mut().pop(); });
        }
    }
    let py = obj.py();
    let identifier = obj.as_ptr() as usize;
    if ACTIVE_CONVERSIONS.with(|active| active.borrow().contains(&identifier)) {
        return Err(pyo3::exceptions::PyRecursionError::new_err("cyclic symbolic conversion hook"));
    }
    // Hook calls return before recursive conversion, so Python's ordinary
    // frame counter alone cannot protect cycles between converted objects.
    if unsafe { pyo3::ffi::Py_EnterRecursiveCall(b" while converting a symbolic object\0".as_ptr().cast()) } != 0 {
        return Err(PyErr::fetch(py));
    }
    ACTIVE_CONVERSIONS.with(|active| active.borrow_mut().push(identifier));
    let _guard = RecursionGuard;
    let limit: usize = py.import("sys")?.call_method0("getrecursionlimit")?.extract()?;
    let mut seen = FxHashSet::default();
    seen.insert(identifier);
    let mut owners = vec![obj.clone()];
    let mut converted = hook.call0()?;
    loop {
        if converted.is(owners.last().unwrap()) {
            return Err(PyTypeError::new_err("symbolic conversion hook returned its input"));
        }
        if !seen.insert(converted.as_ptr() as usize) || owners.len() >= limit {
            return Err(pyo3::exceptions::PyRecursionError::new_err("cyclic or unbounded symbolic conversion hooks"));
        }
        // Retain owners so addresses cannot be recycled during cycle checks.
        owners.push(converted.clone());
        if let Some(hook) = direct_conversion_hook(&converted)? {
            converted = hook.call0()?;
        } else {
            return to_expr(&converted);
        }
    }
}

/// Convert any supported Python object into an expression.
pub(crate) fn to_expr(obj: &Bound<'_, PyAny>) -> PyResult<Expr> {
    if let Ok(p) = obj.downcast::<crate::PyExprNode>() {
        return Ok(p.get().e.clone());
    }
    if let Ok(element) = obj.extract::<PyRef<'_, poly::domains::PyDomainElement>>() {
        use symox::poly::domains::{ElementValue, ExactDomain};
        use symox::poly::int::Coeff;
        return match (&element.value.value, &element.value.domain) {
            (ElementValue::Rational(_), _) => Ok(element.value.as_expr()),
            (ElementValue::Inexact(e)|ElementValue::Expression(e),_) => Ok(e.clone()),
            (ElementValue::Algebraic(_),ExactDomain::Gaussian{..}) => Ok(element.value.as_expr()),
            (ElementValue::Finite(value), ExactDomain::GF(prime)) => {
                let rational = value.to_rat();
                let residue = rational.n.modulo(&symox::poly::Int::from_u64(*prime));
                Ok(Expr::num(residue.to_number()))
            }
            _ => Err(PyTypeError::new_err("convert this coefficient explicitly with its domain.to_sympy method")),
        };
    }
    if let Ok(p) = obj.downcast::<matrixexpr::PyMatrixExpr>() {
        return Ok(matrix_proxy(p.get().m.clone()));
    }
    if obj.is_instance_of::<pyo3::types::PyBool>() {
        return Ok(Expr::boolean(obj.extract::<bool>()?));
    }
    if obj.is_instance_of::<PyInt>() {
        if let Ok(i) = obj.extract::<i64>() {
            return Ok(Expr::int(i));
        }
        let b: BigInt = obj.extract()?;
        return Ok(Expr::num(Number::from_bigint(b)));
    }
    if obj.is_instance_of::<PyFloat>() {
        return Ok(Expr::float(obj.extract::<f64>()?));
    }
    if let Ok(c) = obj.downcast::<PyComplex>() {
        let (re, im) = (c.real(), c.imag());
        return Ok(evalf::c64_to_expr(num_complex_c64(re, im)));
    }
    if obj.is_instance_of::<PyString>() {
        let s: String = obj.extract()?;
        return parser::parse_str(&s).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    // fractions.Fraction
    let ty = obj.get_type();
    let tname = ty.name()?.to_string();
    let module = ty.module().map(|m| m.to_string()).unwrap_or_default();
    if tname == "Fraction" {
        let n: BigInt = obj.getattr("numerator")?.extract()?;
        let d: BigInt = obj.getattr("denominator")?.extract()?;
        return Number::rational(n, d)
            .map(Expr::num)
            .ok_or_else(|| PyZeroDivisionError::new_err("zero denominator"));
    }
    // numpy scalars
    if module == "numpy" {
        if let Ok(i) = obj.extract::<i64>() {
            if tname.contains("int") {
                return Ok(Expr::int(i));
            }
        }
        if let Ok(f) = obj.extract::<f64>() {
            return Ok(Expr::float(f));
        }
    }
    // SymPy objects (and anything with a _symox_ hook)
    if let Ok(hook) = obj.getattr("_symox_") {
        return conversion_hook(obj, &hook);
    }
    // Applied undefined SymPy functions can have __module__ = None (or a
    // user module). Identify their SymPy ancestry without importing SymPy
    // during ordinary native expression construction.
    let inherited_sympy = if module.starts_with("sympy") {
        false
    } else if let Ok(mro) = ty.getattr("__mro__") {
        mro.try_iter()?.filter_map(Result::ok).any(|base| {
            base.getattr("__module__").and_then(|m| m.extract::<String>())
                .is_ok_and(|module| module.starts_with("sympy"))
        })
    } else { false };
    if module.starts_with("sympy") || inherited_sympy {
        let py = obj.py();
        let head: String = obj.getattr("func")?.getattr("__name__")?.extract()?;
        if head == "ComplexRootOf" && obj.hasattr("poly")? {
            let polynomial = obj.getattr("poly")?;
            let expression = to_expr(&polynomial.call_method0("as_expr")?)?;
            let generator = to_expr(&polynomial.getattr("gen")?)?;
            let index: i64 = obj.getattr("index")?.extract()?;
            let root = py.import("symox")?.getattr("CRootOf")?.call1((wrap(expression), wrap(generator), index))?;
            return to_expr(&root);
        }
        let applied_undefined=obj.get_type().getattr("__mro__").is_ok_and(|mro|mro.try_iter().is_ok_and(|bases|bases.filter_map(Result::ok).any(|base|base.getattr("__name__").is_ok_and(|name|name.extract::<String>().is_ok_and(|name|name=="AppliedUndef")))));
        if applied_undefined {
            let spec=sympy_function_spec(&obj.getattr("func")?)?;
            let args=obj.getattr("args")?.try_iter()?.map(|arg|to_expr(&arg?)).collect::<PyResult<Vec<_>>>()?;
            return Ok(fx::func(Func::Undefined(spec),args));
        }
        if head.ends_with("Predicate") && head != "AppliedPredicate" {
            let name = obj.getattr("name")?.str()?.to_string();
            let p = if head == "UndefinedPredicate" {symox::ask::Pred::undefined(&name)}
                else {symox::ask::Pred::from_name(&name).ok_or_else(||PyValueError::new_err("unknown predicate"))?};
            return Ok(symox::ask::predicate(p));
        }
        if head=="AppliedPredicate" || head=="AppliedBinaryRelation" {
            let args=obj.getattr("args")?.try_iter()?.map(|a|to_expr(&a?)).collect::<PyResult<Vec<_>>>()?;
            return symox::ask::applied_from_args(args).map_err(ask::error);
        }
        if head=="Tuple" {
            return Ok(fx::undefined("Tuple",obj.getattr("args")?.try_iter()?.map(|arg|to_expr(&arg?)).collect::<PyResult<Vec<_>>>()?));
        }
        if head=="Integral" || head=="Sum" || head=="Product" {
            let body=to_expr(&obj.getattr("function")?)?;
            let limits=obj.getattr("limits")?.try_iter()?.map(|limit|limit?.try_iter()?.map(|value|to_expr(&value?)).collect::<PyResult<Vec<_>>>()).collect::<PyResult<Vec<_>>>()?;
            return symox::binders::held(match head.as_str() {"Integral"=>Func::Integral,"Sum"=>Func::Sum,_=>Func::Product},&body,limits).map_err(PyValueError::new_err);
        }
        if head=="IdentityFunction" {return Ok(symox::sets::identity_function());}
        if head=="Lambda" {
            return symox::sets::lambda_checked(&to_expr(&obj.getattr("signature")?)?,&to_expr(&obj.getattr("expr")?)?).map_err(logic::bad_signature);
        }
        if head=="Subs" {
            let args=obj.getattr("args")?;
            let expression=to_expr(&args.get_item(0)?)?;
            let variables=to_expr(&args.get_item(1)?)?.args();
            let points=to_expr(&args.get_item(2)?)?.args();
            if variables.iter().any(|v|!v.is_symbol()){return Err(PyNotImplementedError::new_err("Subs currently supports symbol variables"));}
            return Ok(symox::subs::held_subs(&expression,variables,points));
        }
        if head=="Derivative" {
            let expression=to_expr(&obj.getattr("expr")?)?;
            let mut variables=Vec::new();
            for item in obj.getattr("variable_count")?.try_iter()?{let item=item?;variables.push((to_expr(&item.get_item(0)?)?,item.get_item(1)?.extract()?));}
            return Ok(diff::derivative_unevaluated(&expression,variables));
        }
        if head=="Limit" {
            let args=obj.getattr("args")?;
            let e=to_expr(&args.get_item(0)?)?;let x=to_expr(&args.get_item(1)?)?;let point=to_expr(&args.get_item(2)?)?;
            let direction=args.get_item(3)?.str()?.to_string();
            let direction=symox::limits::Dir::parse(&direction).ok_or_else(||PyValueError::new_err("invalid limit direction"))?;
            let held=symox::limits::held_limit(&e,&x,&point,direction);
            symox::limits::validate_held(&held).map_err(PyNotImplementedError::new_err)?;
            return Ok(held);
        }
        if head=="DotProduct" {
            let args=obj.getattr("args")?;
            return symox::matrixcontraction::dot_product(&matrixexpr::to_matrixexpr(&args.get_item(0)?)?,&matrixexpr::to_matrixexpr(&args.get_item(1)?)?,false).map_err(|error|matrices::err(py,error));
        }
        if head == "MatrixElement" {
            let args = obj.getattr("args")?;
            let parent = matrixexpr::to_matrixexpr(&args.get_item(0)?)?;
            return parent.element(&to_expr(&args.get_item(1)?)?, &to_expr(&args.get_item(2)?)?)
                .map_err(|error| matrices::err(py, error));
        }
        if matches!(head.as_str(), "Trace" | "Determinant") {
            let args = obj.getattr("args")?;
            let parent = std::sync::Arc::new(matrixexpr::to_matrixexpr(&args.get_item(0)?)?);
            let f = if head == "Trace" { Func::MatrixTrace(parent) } else { Func::MatrixDeterminant(parent) };
            return Ok(Expr::from_kind(Kind::Func(f, vec![])));
        }
        if obj.getattr("is_Matrix").is_ok_and(|value| value.is_truthy().unwrap_or(false)) {
            return Ok(matrix_proxy(matrixexpr::to_matrixexpr(obj)?));
        }
        if head == "Symbol" || head == "Dummy" || head == "Wild" {
            let name: String = obj.getattr("name")?.extract()?;
            let assumptions = obj.getattr("_assumptions_orig")?;
            let assumptions = Some(assumptions.downcast::<PyDict>()?);
            return if head == "Dummy" {
                let index: BigInt = obj.getattr("dummy_index")?.extract()?;
                Ok(py_dummy(Some(&name), Some(index), assumptions)?.e)
            } else if head == "Wild" {
                Ok(py_wild(&name,Some(&obj.getattr("exclude")?),Some(&obj.getattr("properties")?),assumptions)?.e)
            } else { Ok(py_symbol(&name, assumptions)?.e) };
        }
        if matches!(head.as_str(), "Add" | "Mul" | "Pow" | "UnevaluatedExpr") {
            let args = obj.getattr("args")?;
            let args = args.downcast::<PyTuple>()?;
            let evaluate = obj.eq(obj.getattr("func")?.call1(args)?)?;
            let v = args.iter().map(|a| to_expr(&a)).collect::<PyResult<Vec<_>>>()?;
            return construct_node(&head, v, evaluate);
        }
        if !matches!(head.as_str(), "Piecewise" | "Order") && (Func::from_name(&head).is_some() || obj.getattr("is_Function").is_ok_and(|v| v.is_truthy().unwrap_or(false))) {
            let args = obj.getattr("args")?;
            let mut v = args.try_iter()?.map(|a| to_expr(&a?)).collect::<PyResult<Vec<_>>>()?;
            if head == "Xor" && v.iter().any(|value| matches!(value.kind(), Kind::Bool(true))) {
                symox::logic::ordered(&mut v);
                return Ok(Expr::from_kind(Kind::Func(Func::Logic(symox::logic::LogicOp::Xor), v)));
            }
            return Ok(fx::func(Func::from_name(&head).unwrap_or_else(|| Func::Undefined(head.into())), v));
        }
        let srepr = py.import("sympy")?.getattr("srepr")?.call1((obj,))?;
        let s: String = srepr.extract()?;
        return parser::parse_str(&s).map_err(|e| {
            PyValueError::new_err(format!("cannot convert SymPy object {s}: {e}"))
        });
    }
    // External scalar libraries expose the established SymPy conversion
    // protocol. Builtin reference objects above retain their structural path.
    if module != "symox" && !module.starts_with("symox.") {
        if let Ok(hook) = obj.getattr("_sympy_") {
            return conversion_hook(obj, &hook);
        }
    }
    Err(PyTypeError::new_err(format!(
        "cannot convert object of type {tname} to a symox expression"
    )))
}

fn num_complex_c64(re: f64, im: f64) -> num_complex::Complex64 {
    num_complex::Complex64::new(re, im)
}

fn opt_expr(obj: &Bound<'_, PyAny>) -> Option<Expr> {
    to_expr(obj).ok()
}

// Native working precisions use u32. Leave room for doubled-precision
// refinement and guard bits, instead of imposing a decimal-digit policy cap.
pub(crate) const MAX_REQUESTED_PRECISION: u32 = (u32::MAX - 65_536) / 2;

pub(crate) fn bit_precision(precision: u32) -> PyResult<u32> {
    if precision == 0 {
        return Err(PyValueError::new_err("bit precision must be positive"));
    }
    if precision > MAX_REQUESTED_PRECISION {
        return Err(PyValueError::new_err("precision exceeds the native 32-bit working-precision range"));
    }
    Ok(precision)
}

pub(crate) fn decimal_precision(n: impl Into<i64>) -> PyResult<u32> {
    let n=n.into();
    if n<0 {return Ok(1);}
    let n=u32::try_from(n).map_err(|_|PyValueError::new_err("precision exceeds the native 32-bit working-precision range"))?;
    bit_precision(symox::bigfloat::dps_to_prec(n))
}

pub(crate) fn evalf_options(kwargs: Option<&Bound<'_, PyDict>>,precision:u32) -> PyResult<Option<i64>> {
    let mut chop=None;
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs {
            let key: String = key.extract()?;
            match key.as_str() {
                "chop"=>{
                    let flag=if value.is_instance_of::<pyo3::types::PyBool>() {Some(value.extract::<bool>()?)}
                        else if let Ok(expression)=value.extract::<PyRef<'_, crate::PyExprNode>>() {
                            if let Kind::Bool(value)=expression.e.kind(){Some(*value)}else{None}
                        } else {None};
                    if let Some(flag)=flag {chop=flag.then_some(-(precision as i64));continue;}
                    if !value.is_truthy()? {continue;}
                    let tolerance=to_expr(&value)?;
                    let error=||PyValueError::new_err("chop tolerance must be finite, real and positive");
                    let tolerance=symox::evalf_mp::eval_mp(&tolerance,80).filter(|value|value.im.is_zero()&&value.re.sign()>0).ok_or_else(error)?;
                    let log10=tolerance.re.ln(80).ok_or_else(error)?.to_f64()*std::f64::consts::LOG10_E;
                    // Preserve SymPy's historical threshold in the machine
                    // range. Outside it, convert the logarithm itself instead
                    // of overflowing/underflowing the tolerance to f64.
                    let factor=if log10.abs()<=300.0 {3.321}else{std::f64::consts::LOG2_10};
                    let mut bits=(-factor*log10+2.5).round_ties_even() as i64;
                    if bits==3 {bits=2;}
                    chop=Some(4-bits);
                },
                "strict" | "verbose" if !value.is_truthy()? => {},
                "strict" | "verbose" | "maxn" | "quad" => return Err(PyNotImplementedError::new_err(format!("evalf option {key} is not implemented"))),
                _ => return Err(PyTypeError::new_err(format!("unexpected evalf option {key}"))),
            }
        }
    }
    Ok(chop)
}

pub(crate) fn evalf_expr(e: &Expr, n: i64, replacements: Option<&Bound<'_, PyDict>>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Expr> {
    let precision = decimal_precision(n)?;
    let chop=evalf_options(kwargs,precision)?;
    evalf_prepared(e,precision,replacements,chop)
}

pub(crate) fn evalf_prepared(e:&Expr,precision:u32,replacements:Option<&Bound<'_,PyDict>>,chop:Option<i64>)->PyResult<Expr> {
    let e = match replacements {
        Some(mapping) => {
            let pairs = mapping.iter().map(|(old, new)| Ok((to_expr(&old)?, to_expr(&new)?))).collect::<PyResult<Vec<_>>>()?;
            symox::evalf_mp::substitute_for_evalf(e, &pairs)
        }
        None => e.clone(),
    };
    Ok(match chop {Some(cutoff)=>symox::evalf_mp::evalf_chop(&e,precision,cutoff),None=>symox::evalf_mp::evalf(&e, precision)})
}

#[derive(Clone,Copy,PartialEq,Eq)]
enum IntegrationConditions{Piecewise,None,Separate}

fn integrate_options(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<IntegrationConditions> {
    let mut conditions=IntegrationConditions::Piecewise;
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs {
            let key: String = key.extract()?;
            match key.as_str() {
                "meijerg" | "risch" | "heurisch" | "manual" if value.is_none() => {},
                "conds"=>conditions=match value.extract::<String>().ok().as_deref(){
                    Some("piecewise")=>IntegrationConditions::Piecewise,
                    Some("none")=>IntegrationConditions::None,
                    Some("separate")=>IntegrationConditions::Separate,
                    _=>return Err(PyValueError::new_err("conds must be one of separate, piecewise, none")),
                },
                "meijerg" | "risch" | "heurisch" | "manual" => return Err(PyNotImplementedError::new_err(format!("integration option {key} is not implemented"))),
                _ => return Err(PyTypeError::new_err(format!("unexpected integration option {key}"))),
            }
        }
    }
    Ok(conditions)
}

fn expr_to_sympy<'py>(py: Python<'py>, e: &Expr) -> PyResult<Bound<'py, PyAny>> {
    let sm = py.import("sympy")?;
    let convert_args = |args: &[Expr]| -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, args.iter().map(|a| expr_to_sympy(py, a)).collect::<PyResult<Vec<_>>>()?)
    };
    match e.kind() {
        Kind::Num(number) => match number {
            Number::Int(value) => sm.getattr("Integer")?.call1((*value,)),
            Number::BigInt(value) => sm.getattr("Integer")?.call1((value.clone(),)),
            Number::Rational(value) => sm.getattr("Rational")?.call1((value.numer().clone(),value.denom().clone())),
            Number::Float(value) => sm.getattr("Float")?.call1((*value,)),
            Number::MpFloat(value,precision) => {
                // Binary transfer avoids decimal round trips and Python's
                // configurable limit on decimal integer-string conversion.
                let mpf=(u8::from(value.man.is_negative()),value.man.abs(),value.exp,value.man.bits());
                let kwargs=PyDict::new(py);
                kwargs.set_item("precision",*precision)?;
                // SymPy's mpf-tuple constructor canonicalizes zero to Integer(0).
                // Its ordinary numeric constructor preserves Float and precision.
                if value.man.is_zero(){sm.getattr("Float")?.call((0,),Some(&kwargs))}
                else{sm.getattr("Float")?.call((mpf,),Some(&kwargs))}
            },
        },
        Kind::Func(Func::Special(s),args) if s.name().starts_with("HyperRep_")=>py.import("sympy.functions.special.hyper")?.getattr(s.name())?.call1(convert_args(args)?),
        Kind::Func(Func::Predicate(p),_)=>{
            if p.is_undefined(){sm.getattr("Predicate")?.call1((p.name(),))}
            else{sm.getattr("Q")?.getattr(p.name())}
        },
        Kind::Func(Func::AppliedPredicate,args)=>expr_to_sympy(py,&args[0])?.call1(convert_args(&args[1..])?),
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>py.import("sympy")?.getattr("DotProduct")?.call1((matrixexpr::matrixexpr_to_sympy(py,a)?,matrixexpr::matrixexpr_to_sympy(py,b)?)),
        Kind::Func(Func::MatrixProxy(parent), _) => matrixexpr::matrixexpr_to_sympy(py, parent),
        Kind::Func(f @ (Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent)), _) => {
            sm.getattr(f.name())?.call1((matrixexpr::matrixexpr_to_sympy(py, parent)?,))
        }
        Kind::Func(Func::MatrixElement(parent), args) => {
            let parent = matrixexpr::matrixexpr_to_sympy(py, parent)?;
            py.import("sympy.matrices.expressions.matexpr")?.getattr("MatrixElement")?
                .call1((parent, expr_to_sympy(py, &args[0])?, expr_to_sympy(py, &args[1])?))
        }
        Kind::Sym(s) => {
            let kwargs = PyDict::new(py);
            for keyword in s.given_facts().keywords() {
                if let Some((key, value)) = keyword.split_once('=') { kwargs.set_item(key, value == "True")?; }
            }
            if let Some(index) = s.dummy_index() { kwargs.set_item("dummy_index", index)?; }
            if s.is_wild() { kwargs.set_item("exclude", convert_args(s.exclude())?)?; }
            sm.getattr(if s.is_dummy() { "Dummy" } else if s.is_wild() { "Wild" } else { "Symbol" })?.call((s.name(),), Some(&kwargs))
        }
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => {
            let args = convert_args(&e.args())?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("evaluate", e.is_canonical())?;
            sm.getattr(e.class_name())?.call(args, Some(&kwargs))
        }
        Kind::Func(Func::Undefined(name), args) if name.is_internal("Tuple") => sm.getattr("Tuple")?.call1(convert_args(args)?),
        Kind::Func(Func::Undefined(name), args)
            if name.name.as_ref() == "ComplexRootOf" && args.len() == 2
                && name.facts().get("algebraic") == Some(true) => {
            let coefficients = convert_args(&args[0].args())?;
            let generator = sm.getattr("Dummy")?.call1(("x",))?;
            let polynomial = sm.getattr("Poly")?.getattr("from_list")?.call1((coefficients, generator))?;
            sm.getattr("CRootOf")?.call1((polynomial, expr_to_sympy(py, &args[1])?))
        }
        Kind::Func(Func::Undefined(name), args) if !name.is_internal("UnevaluatedExpr") => {
            let kwargs=function_kwargs(py,name)?;
            sm.getattr("Function")?.call((name.as_ref(),),Some(&kwargs))?.call1(convert_args(args)?)
        }
        Kind::Func(Func::Order, _) => sm.getattr("sympify")?.call1((printer::srepr(e),)),
        Kind::Func(Func::Logic(symox::logic::LogicOp::Xor), args) => {
            let kwargs = PyDict::new(py);
            kwargs.set_item("remove_true", false)?;
            sm.getattr("Xor")?.call(convert_args(args)?, Some(&kwargs))
        }
        Kind::Func(Func::Set(symox::sets::SetOp::IdentityFunction),_)=>sm.getattr("S")?.getattr("IdentityFunction"),
        Kind::Func(Func::Set(op), args) if op.is_named() && args.is_empty() => sm.getattr(op.name()),
        Kind::Func(f, args) => sm.getattr(f.name())?.call1(convert_args(args)?),
        Kind::Rel(..) => {
            let kwargs=PyDict::new(py);
            kwargs.set_item("evaluate",e.is_canonical())?;
            sm.getattr(e.class_name())?.call(convert_args(&e.args())?,Some(&kwargs))
        },
        Kind::Piecewise(pairs) => {
            let args = pairs.iter().map(|(value, condition)| {
                PyTuple::new(py, [expr_to_sympy(py, value)?, expr_to_sympy(py, condition)?])
            }).collect::<PyResult<Vec<_>>>()?;
            sm.getattr("Piecewise")?.call1(PyTuple::new(py, args)?)
        }
        Kind::Derivative(value, variables) => {
            let mut args = vec![expr_to_sympy(py, value)?];
            for (variable, order) in variables {
                args.push(PyTuple::new(py, [expr_to_sympy(py, variable)?, order.into_pyobject(py)?.into_any()])?.into_any());
            }
            sm.getattr("Derivative")?.call1(PyTuple::new(py, args)?)
        }
        _ => sm.getattr("sympify")?.call1((printer::srepr(e),)),
    }
}

fn python_args(e: &Expr) -> Vec<Expr> { traverse::ordered_args(e) }

fn operation_count(e: &Expr) -> BigInt {
    fn matrix_count(matrix:&symox::matrixexpr::MatrixExpr,scalars:&mut FxHashMap<usize,(Expr,BigInt)>,matrices:&mut FxHashMap<usize,(symox::matrixexpr::MatrixExpr,BigInt)>)->BigInt{
        use symox::matrixexpr::{MatrixKind as M,MatrixOperand};
        if let Some((_,n))=matrices.get(&matrix.dag_id()){return n.clone();}
        let mut n:BigInt=matrix.children().iter().map(|m|matrix_count(m,scalars,matrices)).sum();
        match matrix.kind(){
            M::Symbol(_)|M::Identity|M::Zero|M::One=>{n+=count(matrix.rows(),scalars)+count(matrix.cols(),scalars);},
            M::Function(lambda)=>{n+=count(matrix.rows(),scalars)+count(matrix.cols(),scalars)+count(lambda,scalars);},
            M::ApplyFunction(function,_)=>n+=count(function,scalars),
            M::Explicit(m)=>n+=m.data().iter().map(|e|count(e,scalars)).sum::<BigInt>(),
            M::Add(args)=>n+=args.len().saturating_sub(1),
            M::Mul(coefficient,args)=>{n+=count(coefficient,scalars)+BigInt::from((args.len()+usize::from(!coefficient.is_one())).saturating_sub(1));},
            M::Pow(_,e)=>n+=count(e,scalars),M::Inverse(_)=>n+=1,
            M::HadamardPower(a,b)=>{for value in [a,b]{if let MatrixOperand::Scalar(e)=value{n+=count(e,scalars);}}},
            M::Slice(_,a,b)=>{for axis in [a,b]{for arg in axis.args().into_iter().flatten(){n+=count(&arg,scalars);}}},
            _=>{},
        }
        matrices.insert(matrix.dag_id(),(matrix.clone(),n.clone()));n
    }
    fn count(e: &Expr, memo: &mut FxHashMap<usize, (Expr, BigInt)>) -> BigInt {
        if let Some((_, value)) = memo.get(&e.dag_id()) { return value.clone(); }
        let value = match e.kind() {
        Kind::Num(n) => BigInt::from(usize::from(n.is_rational() && n.is_negative()) + usize::from(matches!(n, Number::Rational(_)))),
        Kind::Const(Constant::E) => BigInt::from(1),
        Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => BigInt::from(0),
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>{let mut matrices=FxHashMap::default();matrix_count(a,memo,&mut matrices)+matrix_count(b,memo,&mut matrices)},
        Kind::Func(Func::MatrixProxy(_), _) => BigInt::from(0),
        Kind::Func(Func::MatrixElement(_), args) => args.iter().map(|a| count(a, memo)).sum(),
        Kind::Add(c, terms) => {
            let mut result = BigInt::from((terms.len() + usize::from(!c.is_zero())).saturating_sub(1));
            let all_negative = (c.is_zero() || c.is_negative()) && terms.iter().all(|(_, k)| k.is_negative());
            result += usize::from(all_negative);
            if !c.is_zero() { result += usize::from(matches!(c, Number::Rational(_))); }
            for (term, coefficient) in terms {
                result += count(&ops::add_term(term, &coefficient.abs()), memo);
            }
            result
        }
        Kind::Mul(c, factors) => {
            let abs = c.abs();
            let mut numerator = Vec::new();
            let mut denominator = Vec::new();
            for (base, exponent) in factors {
                if exponent.as_number().is_some_and(|n| n.is_negative()) {
                    denominator.push(ops::pow(base, &ops::neg(exponent)));
                } else {
                    numerator.push(ops::mul_factor(base, exponent));
                }
            }
            if abs.is_rational() && !abs.is_integer() {
                let p=abs.numer();
                if !p.is_one(){numerator.push(Expr::num(p));}
                denominator.push(Expr::num(abs.denom()));
            }else if !abs.is_one(){numerator.push(Expr::num(abs));}
            BigInt::from(usize::from(c.is_negative()))
                + numerator.iter().map(|a| count(a, memo)).sum::<BigInt>()
                + denominator.iter().map(|a| count(a, memo)).sum::<BigInt>()
                + numerator.len().saturating_sub(1)
                + denominator.len().saturating_sub(1)
                + usize::from(!denominator.is_empty())
        }
        Kind::Pow(base, exponent) => {
            if base.is_constant(Constant::E){1+count(exponent,memo)}
            else{count(base, memo) + if exponent.is_minus_one() { BigInt::from(1) } else { 1 + count(exponent, memo) }}
        }
        Kind::Func(Func::Subs|Func::Limit|Func::Order|Func::Set(symox::sets::SetOp::IdentityFunction),args) => args.iter().map(|a|count(a,memo)).sum(),
        Kind::Func(Func::Undefined(n), args) if n.is_internal("UnevaluatedExpr") || n.is_internal("Tuple") => args.iter().map(|a| count(a, memo)).sum(),
        _ => 1 + e.args().iter().map(|a| count(a, memo)).sum::<BigInt>(),
        };
        // Synthetic coefficient/power nodes also enter this cache. Retain an
        // owner so their allocation addresses cannot be recycled mid-walk.
        memo.insert(e.dag_id(), (e.clone(), value.clone()));
        value
    }
    count(e, &mut FxHashMap::default())
}

fn unwrap_unevaluated(e: &Expr) -> Expr {
    // Memoized per DAG node: shared subtrees (Lagrangian equations, LU
    // solutions) would otherwise be walked once per path.
    fn walk(e: &Expr, memo: &mut FxHashMap<usize, (Expr, Expr)>) -> Expr {
        if e.is_atom() { return e.clone(); }
        if let Some((_, r)) = memo.get(&e.dag_id()) { return r.clone(); }
        let r = match e.kind() {
            Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") => walk(&args[0], memo),
            _ => traverse::map_children(e, &mut |c| walk(c, memo)),
        };
        memo.insert(e.dag_id(), (e.clone(), r.clone()));
        r
    }
    walk(e, &mut FxHashMap::default())
}

/// Parse the variable list of `diff(expr, x, y, (z, 2), x, 3)`.
fn valid_diff_variable(variable: &Expr) -> bool { diff::is_valid_variable(variable) }

pub(crate) fn diff_vars(args: &Bound<'_, PyTuple>) -> PyResult<Vec<(Expr, u32)>> {
    let order = |value: &Bound<'_, PyAny>| -> PyResult<u32> {
        if value.is_instance_of::<pyo3::types::PyBool>() {
            return Err(PyValueError::new_err("derivative order must be a nonnegative integer"));
        }
        to_expr(value)?.as_i64().and_then(|n| u32::try_from(n).ok())
            .ok_or_else(|| PyValueError::new_err("derivative order must be a nonnegative integer"))
    };
    let mut vars = Vec::new();
    let items: Vec<_> = args.iter().collect();
    let mut i = 0;
    while i < items.len() {
        let item = &items[i];
        let pair = if let Ok(tuple) = item.downcast::<PyTuple>() {
            if tuple.len() != 2 { return Err(PyValueError::new_err("expected (variable, order) pairs")); }
            Some((to_expr(&tuple.get_item(0)?)?, order(&tuple.get_item(1)?)?))
        } else if let Ok(list) = item.downcast::<PyList>() {
            if list.len() != 2 { return Err(PyValueError::new_err("expected (variable, order) pairs")); }
            Some((to_expr(&list.get_item(0)?)?, order(&list.get_item(1)?)?))
        } else {
            let expression = to_expr(item)?;
            if let Kind::Func(Func::Undefined(head), values) = expression.kind() {
                if head.is_internal("Tuple") {
                    if values.len() != 2 { return Err(PyValueError::new_err("expected (variable, order) pairs")); }
                    let n = values[1].as_i64().and_then(|n|u32::try_from(n).ok())
                        .ok_or_else(||PyValueError::new_err("derivative order must be a nonnegative integer"))?;
                    Some((values[0].clone(), n))
                } else { None }
            } else { None }
        };
        let (variable, n) = if let Some(pair) = pair { pair } else {
            let variable = to_expr(item)?;
            let mut n = 1;
            if let Some(next) = items.get(i+1) {
                if !next.is_instance_of::<pyo3::types::PyBool>() &&
                    to_expr(next).ok().is_some_and(|e|e.as_number().is_some_and(|n|n.is_integer())) {
                    n = order(next)?;
                    i += 1;
                }
            }
            (variable, n)
        };
        if matches!(variable.kind(),Kind::Func(Func::MatrixProxy(_),_)) {
            return Err(PyNotImplementedError::new_err("matrix-variable derivatives require tensor expressions"));
        }
        if !valid_diff_variable(&variable) {
            return Err(PyValueError::new_err("cannot differentiate with respect to this expression"));
        }
        vars.push((variable,n));
        i += 1;
    }
    Ok(vars)
}

fn do_diff(e: &Expr, args: &Bound<'_, PyTuple>) -> PyResult<Expr> {
    let vars = diff_vars(args)?;
    if vars.is_empty() {
        let fs = traverse::free_symbols(e);
        return match fs.len() {
            0 => Ok(Expr::zero()),
            1 if matches!(fs[0].kind(),Kind::Func(Func::MatrixProxy(_),_)) => Err(PyNotImplementedError::new_err("matrix-variable derivatives require tensor expressions")),
            1 => Ok(diff::diff(e, &fs[0])),
            _ => Err(PyValueError::new_err(
                "Since there is more than one variable in the expression, the variable(s) of differentiation must be supplied",
            )),
        };
    }
    Ok(diff::diff_many(e, &vars))
}

/// Bind Basic.subs' public argument names before dispatching to each native type.
/// In the reference contract arg2=None selects the mapping/pair-sequence form.
pub(crate) fn substitution_arguments<'py>(args:&Bound<'py,PyTuple>,kwargs:Option<&Bound<'py,PyDict>>)->PyResult<(Bound<'py,PyTuple>,Option<Bound<'py,PyDict>>)> {
    if args.len()>2{return Err(PyTypeError::new_err("subs expects 1 or 2 arguments"));}
    let options=kwargs.map(|d|d.copy()).transpose()?;
    let mut values=args.iter().collect::<Vec<_>>();
    for (index,name) in ["arg1","arg2"].iter().enumerate(){
        let named=if let Some(options)=&options{options.get_item(name)?}else{None};
        if let Some(value)=named{
            if args.len()>index{return Err(PyTypeError::new_err(format!("subs() got multiple values for argument '{name}'")));}
            if values.len()!=index{return Err(PyTypeError::new_err("subs() missing required argument 'arg1'"));}
            values.push(value);
            options.as_ref().unwrap().del_item(name)?;
        }
    }
    if values.is_empty(){return Err(PyTypeError::new_err("subs() missing required argument 'arg1'"));}
    if values.len()==2 && values[1].is_none(){values.pop();}
    Ok((PyTuple::new(args.py(),values)?,options))
}

pub(crate) fn subs_pairs(args: &Bound<'_, PyTuple>) -> PyResult<Vec<(Expr, Expr)>> {
    fn pair(item:&Bound<'_,PyAny>)->PyResult<(Expr,Expr)> {
        let mut values=item.try_iter().map_err(|_|PyTypeError::new_err("subs expects a mapping, an iterable of pairs or (old, new)"))?;
        let old=values.next().transpose()?.ok_or_else(||PyValueError::new_err("substitution pairs must contain exactly two values"))?;
        let new=values.next().transpose()?.ok_or_else(||PyValueError::new_err("substitution pairs must contain exactly two values"))?;
        if values.next().transpose()?.is_some() {return Err(PyValueError::new_err("substitution pairs must contain exactly two values"));}
        Ok((replacement_expr(&old)?,replacement_expr(&new)?))
    }
    let mut pairs = Vec::new();
    match args.len() {
        1 => {
            let a = args.get_item(0)?;
            let mapping=if let Ok(d) = a.downcast::<PyDict>() {
                for (k, v) in d.iter() {
                    pairs.push((replacement_expr(&k)?, replacement_expr(&v)?));
                }
                true
            } else {
                match a.getattr("items") {
                    Ok(items)=>{
                        for item in items.call0()?.try_iter()? {pairs.push(pair(&item?)?);}
                        true
                    },
                    Err(error) if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(a.py())=>{
                        for item in a.try_iter()? {pairs.push(pair(&item?)?);}
                        false
                    },
                    Err(error)=>return Err(error),
                }
            };
            if mapping {
                // A mapping is unordered in SymPy's substitution API: replace
                // larger patterns before symbols they contain. Explicit pair
                // sequences retain the caller's requested sequential order.
                pairs.sort_by_cached_key(|(old, _)| (
                    std::cmp::Reverse(traverse::node_count(old)),
                    symox::sortkey::sort_key(old),
                ));
            }
        }
        2 => pairs.push((replacement_expr(&args.get_item(0)?)?, replacement_expr(&args.get_item(1)?)?)),
        _ => return Err(PyTypeError::new_err("subs expects 1 or 2 arguments")),
    }
    Ok(pairs)
}

pub(crate) fn substitution_options(kwargs:Option<&Bound<'_,PyDict>>)->PyResult<bool> {
    let mut simultaneous=false;
    if let Some(kwargs)=kwargs {
        for (key,value) in kwargs.iter() {
            let key:String=key.extract()?;
            match key.as_str() {
                "simultaneous"=>simultaneous=value.is_truthy()?,
                _=>return Err(PyNotImplementedError::new_err(format!("subs option '{key}' is not implemented"))),
            }
        }
    }
    Ok(simultaneous)
}

/// Share scalar binder/dimension validation with entrywise matrix substitution.
pub(crate) fn checked_subs(py:Python<'_>,expression:&Expr,pairs:&[(Expr,Expr)],simultaneous:bool)->PyResult<Expr> {
    let result=if simultaneous {
        for (old,new) in pairs {symox::binders::check_substitution(expression,old,new).map_err(PyValueError::new_err)?;}
        validate_matrix_replacements(py,expression,pairs,true)?;
        subs::simultaneous(expression,pairs)
    } else {
        let mut result=expression.clone();
        for pair in pairs {
            symox::binders::check_substitution(&result,&pair.0,&pair.1).map_err(PyValueError::new_err)?;
            validate_matrix_replacements(py,&result,std::slice::from_ref(pair),false)?;
            result=subs::subs1(&result,&pair.0,&pair.1);
        }
        result
    };
    symox::limits::validate_held(&result).map_err(PyNotImplementedError::new_err)?;
    symox::binders::validate_tree(&result).map_err(PyValueError::new_err)?;
    symox::sets::validate_lambda_tree(&result).map_err(logic::bad_signature)?;
    Ok(result)
}

#[pymethods]
impl PyExprNode {
    /// Live node identity for Python graph adapters preserving public binders.
    #[getter]
    fn _dag_id(&self) -> usize { self.e.dag_id() }

    #[new]
    fn py_new(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(PyExprNode { e: to_expr(obj)? })
    }

    fn __str__(&self) -> String {
        printer::to_str(&self.e)
    }

    fn __repr__(&self) -> String {
        printer::to_str(&self.e)
    }

    fn _repr_latex_(&self) -> String {
        printing::repr_latex(self)
    }

    fn srepr(&self) -> String {
        printer::srepr(&self.e)
    }

    fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
        match self.e.kind() {
            Kind::Num(Number::Float(value)) => return value.into_pyobject(py)?.hash(),
            Kind::Num(number) => if let Some(hash) = python_hash::exact_hash(number) { return Ok(hash); },
            Kind::Bool(value) => return Ok(isize::from(*value)),
            Kind::Const(Constant::Infinity) => return Ok(314159),
            Kind::Const(Constant::NegInfinity) => return Ok(-314159),
            _ => {},
        }
        let mut h = DefaultHasher::new();
        self.e.hash(&mut h);
        let hash = h.finish() as isize;
        Ok(if hash == -1 { -2 } else { hash })
    }

    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: pyo3::basic::CompareOp) -> PyResult<PyObject> {
        let py = other.py();
        use pyo3::basic::CompareOp::*;
        // AlgebraicNumber retains its defining field and alias as identity.
        // Its _symox_ conversion is for arithmetic, not structural equality.
        if other.downcast::<poly::algebraic_number::PyAlgebraicNumber>().is_ok() {
            return Ok(match op {
                Eq => false.into_pyobject(py)?.to_owned().into_any().unbind(),
                Ne => true.into_pyobject(py)?.to_owned().into_any().unbind(),
                _ => py.NotImplemented(),
            });
        }
        let converted = if let Ok(element) = other.extract::<PyRef<'_, poly::domains::PyDomainElement>>() {
            if matches!(element.value.value,symox::poly::domains::ElementValue::Inexact(_)) {
                let reverse=match op {Lt=>Gt,Le=>Ge,Gt=>Lt,Ge=>Le,other=>other};
                return element.__richcmp__(&self.clone().into_pyobject(py)?.into_any(),reverse,py);
            }
            if matches!(element.value.value, symox::poly::domains::ElementValue::Rational(_)) {
                Some(element.value.as_expr())
            } else { None }
        } else { opt_expr(other) };
        let Some(o) = converted else {
            return Ok(match op {
                Eq => false.into_pyobject(py)?.to_owned().into_any().unbind(),
                Ne => true.into_pyobject(py)?.to_owned().into_any().unbind(),
                _ => py.NotImplemented(),
            });
        };
        Ok(match op {
            Eq => (self.e == o).into_pyobject(py)?.to_owned().into_any().unbind(),
            Ne => (self.e != o).into_pyobject(py)?.to_owned().into_any().unbind(),
            Lt => wrap(relational::lt(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
            Le => wrap(relational::le(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
            Gt => wrap(relational::gt(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
            Ge => wrap(relational::ge(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
        })
    }

    fn __bool__(&self) -> PyResult<bool> {
        match self.e.kind() {
            Kind::Bool(b) => Ok(*b),
            Kind::Num(n) => Ok(!n.is_zero()),
            Kind::Rel(..) => Err(PyTypeError::new_err("cannot determine truth value of Relational")),
            _ => Ok(true),
        }
    }

    // arithmetic -------------------------------------------------------

    fn __add__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if o.downcast::<matrixexpr::PyMatrixExpr>().is_ok() { return Ok(o.py().NotImplemented()); }
        arithmetic_object(&slf, o, ops::add2(&slf.e, &scalar_expr(o)?))
    }
    fn __radd__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if o.downcast::<matrixexpr::PyMatrixExpr>().is_ok() { return Ok(o.py().NotImplemented()); }
        arithmetic_object(&slf, o, ops::add2(&scalar_expr(o)?, &slf.e))
    }
    fn __sub__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if o.downcast::<matrixexpr::PyMatrixExpr>().is_ok() { return Ok(o.py().NotImplemented()); }
        arithmetic_object(&slf, o, ops::sub(&slf.e, &scalar_expr(o)?))
    }
    fn __rsub__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if o.downcast::<matrixexpr::PyMatrixExpr>().is_ok() { return Ok(o.py().NotImplemented()); }
        arithmetic_object(&slf, o, ops::sub(&scalar_expr(o)?, &slf.e))
    }
    fn __mul__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        let py = o.py();
        if slf.e.is_integer_number() && (o.is_instance_of::<PyTuple>() || o.is_instance_of::<PyList>()
            || o.is_instance_of::<PyString>() || o.is_instance_of::<pyo3::types::PyBytes>()) {
            return Ok(py.NotImplemented());
        }
        if o.downcast::<matrices::PyMatrix>().is_ok() || o.downcast::<matrixexpr::PyMatrixExpr>().is_ok() || o.getattr("_symox_sparse").is_ok_and(|marker| marker.is_truthy().unwrap_or(false)) {
            // let Matrix.__rmul__ handle scalar * matrix
            return Ok(py.NotImplemented());
        }
        arithmetic_object(&slf, o, ops::mul2(&slf.e, &scalar_expr(o)?))
    }
    fn __rmul__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if slf.e.is_integer_number() && (o.is_instance_of::<PyTuple>() || o.is_instance_of::<PyList>()
            || o.is_instance_of::<PyString>() || o.is_instance_of::<pyo3::types::PyBytes>()) {
            return Ok(o.py().NotImplemented());
        }
        if o.downcast::<matrices::PyMatrix>().is_ok() || o.downcast::<matrixexpr::PyMatrixExpr>().is_ok() || o.getattr("_symox_sparse").is_ok_and(|marker| marker.is_truthy().unwrap_or(false)) {
            return Ok(o.py().NotImplemented());
        }
        arithmetic_object(&slf, o, ops::mul2(&scalar_expr(o)?, &slf.e))
    }
    fn __truediv__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        arithmetic_object(&slf, o, ops::div(&slf.e, &scalar_expr(o)?))
    }
    fn __rtruediv__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        arithmetic_object(&slf, o, ops::div(&scalar_expr(o)?, &slf.e))
    }
    fn __mod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(checked_mod("Mod", &self.e, &scalar_expr(o)?)?))
    }
    fn __rmod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(checked_mod("Mod", &scalar_expr(o)?, &self.e)?))
    }
    fn __floordiv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        floor_div(&self.e, &scalar_expr(o)?).map(wrap)
    }
    fn __rfloordiv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        floor_div(&scalar_expr(o)?, &self.e).map(wrap)
    }
    fn __divmod__(&self, o: &Bound<'_, PyAny>) -> PyResult<(PyExpr, PyExpr)> {
        let q = scalar_expr(o)?;
        Ok((wrap(floor_div(&self.e, &q)?), wrap(checked_mod("Mod", &self.e, &q)?)))
    }
    fn __rdivmod__(&self, o: &Bound<'_, PyAny>) -> PyResult<(PyExpr, PyExpr)> {
        let p = scalar_expr(o)?;
        Ok((wrap(floor_div(&p, &self.e)?), wrap(checked_mod("Mod", &p, &self.e)?)))
    }
    fn __pow__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>, m: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
        if m.is_some_and(|m| !m.is_none()) {
            return Err(PyTypeError::new_err("modular pow is not supported"));
        }
        arithmetic_object(&slf, o, ops::pow(&slf.e, &scalar_expr(o)?))
    }
    fn __rpow__(slf: PyRef<'_, Self>, o: &Bound<'_, PyAny>, m: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
        if m.is_some_and(|m| !m.is_none()) {
            return Err(PyTypeError::new_err("modular pow is not supported"));
        }
        arithmetic_object(&slf, o, ops::pow(&scalar_expr(o)?, &slf.e))
    }
    fn __neg__(&self, py: Python<'_>) -> PyResult<PyObject> {
        expr_object(py, ops::neg(&self.e))
    }
    fn __pos__(slf: PyRef<'_, Self>) -> PyResult<PyObject> {
        let py = slf.py();
        Ok(slf.into_pyobject(py)?.into_any().unbind())
    }
    fn __abs__(&self, py: Python<'_>) -> PyResult<PyObject> {
        expr_object(py, fx::abs(&self.e))
    }

    fn __float__(&self) -> PyResult<f64> {
        evalf::eval_f64(&self.e).ok_or_else(|| PyTypeError::new_err(format!("cannot convert {} to float", self.e)))
    }

    fn __complex__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyComplex>> {
        let z = evalf::eval_c64(&self.e)
            .ok_or_else(|| PyTypeError::new_err(format!("cannot convert {} to complex", self.e)))?;
        Ok(PyComplex::from_doubles(py, z.re, z.im))
    }

    fn __int__(&self) -> PyResult<BigInt> {
        if let Some(integer) = fx::integer_trunc(&self.e) { return Ok(integer); }
        // Keep the established numerical fallback for functions outside the
        // rational enclosure engine, at native MP precision rather than f64.
        // Stability is numerical evidence; it is not an interval certificate.
        let mut precision = symox::evalf_mp::max_precision(&self.e).unwrap_or(0).saturating_add(64).max(128);
        while precision <= 80_000 {
            let value = symox::evalf_mp::eval_mp(&self.e, precision)
                .ok_or_else(|| PyTypeError::new_err("cannot convert this expression to an integer"))?;
            if !value.is_real() { return Err(PyTypeError::new_err("cannot convert a complex value to an integer")); }
            if !value.re.is_zero() && value.re.mag() > precision as i64 - 64 {
                let required = value.re.mag().saturating_add(128);
                if required > 80_000 { break; }
                precision = required as u32;
                continue;
            }
            let refined = symox::evalf_mp::eval_mp(&self.e, precision * 2)
                .ok_or_else(|| PyTypeError::new_err("cannot refine the expression's integer part"))?;
            if refined.is_real() && value.re.to_int_trunc() == refined.re.to_int_trunc() {
                return Ok(refined.re.to_int_trunc());
            }
            precision *= 2;
        }
        Err(PyValueError::new_err("integer part could not be resolved within the numerical precision budget"))
    }

    fn __index__(&self) -> PyResult<BigInt> {
        self.e
            .as_number()
            .and_then(|n| n.to_bigint())
            .ok_or_else(|| PyTypeError::new_err("only integers can be used as an index"))
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyTuple>)> {
        // Encoding physical nodes once preserves scalar sharing just as it
        // does mixed matrix expressions. Recursive child pickles duplicate a
        // shared binary DAG exponentially and can lose native metadata.
        matrix_dag::reduce_expr(py, &self.e)
    }

    // structure --------------------------------------------------------

    #[getter]
    fn args<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if let Kind::Piecewise(branches)=self.e.kind(){
            let pair=py.import("symox._matrixexpr_methods")?.getattr("_ExprCondPairView")?;
            let arguments=branches.iter().map(|(value,condition)|pair.call1((
                expr_object(py,value.clone())?,expr_object(py,condition.clone())?))).collect::<PyResult<Vec<_>>>()?;
            return PyTuple::new(py,arguments);
        }
        if !self.e.has_matrix() { return PyTuple::new(py,python_args(&self.e).into_iter().map(wrap)); }
        if let Kind::Func(Func::MatrixDotProduct(a,b),_)=self.e.kind(){return PyTuple::new(py,[matrixexpr::matrix_object(py,a)?,matrixexpr::matrix_object(py,b)?]);}
        if let Kind::Func(Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent), _) = self.e.kind() {
            return PyTuple::new(py, [expr_object(py, matrix_proxy((**parent).clone()))?]);
        }
        if let Kind::Func(Func::MatrixElement(parent), args) = self.e.kind() {
            return PyTuple::new(py, [
                expr_object(py, matrix_proxy((**parent).clone()))?,
                expr_object(py, args[0].clone())?, expr_object(py, args[1].clone())?,
            ]);
        }
        PyTuple::new(py, python_args(&self.e).into_iter().map(|e|expr_object(py,e)).collect::<PyResult<Vec<_>>>()?)
    }

    /// SymPy class name of the node (`Add`, `Mul`, `Symbol`, `sin`, ...).
    #[getter]
    fn class_name(&self) -> String {
        self.e.class_name().to_string()
    }

    /// The head of the expression as a callable rebuilding it from args.
    #[getter]
    fn func(&self, py: Python<'_>) -> PyResult<PyObject> {
        if matches!(self.e.kind(),Kind::Piecewise(_)){
            return Ok(py.import("symox")?.getattr("Piecewise")?.unbind());
        }
        if matches!(self.e.kind(),Kind::Func(Func::AppliedPredicate,_)){
            return Ok(py.import("symox.assumptions")?.getattr("AppliedPredicate")?.unbind());
        }
        if let Kind::Func(Func::Undefined(function),_)=self.e.kind(){
            if !function.internal{return Ok(undefined_head(py,function)?.unbind());}
            if function.is_internal("Tuple"){return Ok(py.import("symox._compat")?.getattr("Tuple")?.unbind());}
        }
        if let Kind::Func(Func::Set(op @ (symox::sets::SetOp::Lambda|symox::sets::SetOp::IdentityFunction)),_)=self.e.kind(){return Ok(py.import("symox")?.getattr(op.name())?.unbind());}
        if let Kind::Func(f @ (Func::Integral|Func::Sum|Func::Product),_)=self.e.kind(){return Ok(py.import("symox")?.getattr(f.name())?.unbind());}
        if matches!(self.e.kind(),Kind::Func(Func::Subs,_)){return Ok(py.import("symox")?.getattr("Subs")?.unbind());}
        if matches!(self.e.kind(),Kind::Func(Func::Limit,_)){return Ok(py.import("symox")?.getattr("Limit")?.unbind());}
        if let Kind::Func(f @ (Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _)), _) = self.e.kind() {
            return Ok(py.import("symox._matrixexpr")?.getattr(f.name())?.unbind());
        }
        if matches!(self.e.kind(), Kind::Func(Func::MatrixElement(_), _)) {
            return Ok(py.import("symox._matrixexpr")?.getattr("MatrixElement")?.unbind());
        }
        let name = match self.e.kind() {
            Kind::Add(..) => Some("Add"),
            Kind::Mul(..) => Some("Mul"),
            Kind::Pow(..) => Some("Pow"),
            Kind::Sym(s) => Some(if s.is_dummy() { "Dummy" } else if s.is_wild() { "Wild" } else { "Symbol" }),
            Kind::Num(n) if n.is_integer() => Some("Integer"),
            Kind::Num(n) if n.is_rational() => Some("Rational"),
            Kind::Num(_) => Some("Float"),
            Kind::Func(Func::Undefined(n), _) if n.is_internal("UnevaluatedExpr") => Some("UnevaluatedExpr"),
            _ => None,
        };
        if let Some(name) = name {
            return Ok(py.import("symox._compat")?.getattr(name)?.unbind());
        }
        if let Kind::Func(function, _) = self.e.kind() {
            if !matches!(function, Func::Undefined(_)) {
                if let Ok(public) = py.import("symox")?.getattr(function.name()) {
                    // Native builtin function wrappers carry this marker.
                    // Return their actual identity for expr.func comparisons,
                    // preserving dedicated heads for internal/custom nodes.
                    if public.getattr("_symox_head_name").is_ok_and(|name|
                        name.extract::<String>().is_ok_and(|name| name == function.name())) {
                        return Ok(public.unbind());
                    }
                }
            }
        }
        if let Kind::Rel(op, _, _) = self.e.kind() {
            let name = match op {
                RelOp::Eq => "Eq", RelOp::Ne => "Ne", RelOp::Lt => "Lt",
                RelOp::Le => "Le", RelOp::Gt => "Gt", RelOp::Ge => "Ge",
            };
            return Ok(py.import("symox")?.getattr(name)?.unbind());
        }
        let head = Head::of(&self.e);
        Ok(Py::new(py, head)?.into_any())
    }

    #[getter]
    fn name(&self) -> PyResult<String> {
        match self.e.kind() {
            Kind::Sym(s) => Ok(s.to_string()),
            Kind::Func(Func::Predicate(p), _) => Ok(p.name().to_string()),
            Kind::Func(f, _) => Ok(f.name().to_string()),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("expression has no name")),
        }
    }

    #[getter]
    fn p(&self) -> PyResult<BigInt> {
        self.e
            .as_number()
            .and_then(|n| n.numer_denom())
            .map(|x| x.0)
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("not a rational number"))
    }

    #[getter]
    fn q(&self) -> PyResult<BigInt> {
        self.e
            .as_number()
            .and_then(|n| n.numer_denom())
            .map(|x| x.1)
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("not a rational number"))
    }

    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PySet>> {
        if !self.e.has_matrix() { return PySet::new(py,traverse::free_symbols(&self.e).into_iter().map(wrap)); }
        PySet::new(py, traverse::free_symbols(&self.e).into_iter().map(|e| expr_object(py, e)).collect::<PyResult<Vec<_>>>()?)
    }

    #[getter]
    fn is_Atom(&self) -> bool {
        self.e.is_atom()
    }
    #[getter]
    fn is_Symbol(&self) -> bool {
        self.e.is_symbol()
    }
    #[getter]
    fn is_Dummy(&self) -> bool {
        self.e.as_sym().is_some_and(|s| s.is_dummy())
    }
    #[getter]
    fn is_Wild(&self) -> bool { self.e.as_sym().is_some_and(|s|s.is_wild()) }
    #[getter]
    fn exclude<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let symbol=self.e.as_sym().filter(|s|s.is_wild())
            .ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("only wildcards have exclusions"))?;
        PyTuple::new(py,symbol.exclude().iter().cloned().map(wrap))
    }
    #[getter]
    fn properties<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        if !self.e.as_sym().is_some_and(|s|s.is_wild()) {
            return Err(pyo3::exceptions::PyAttributeError::new_err("only wildcards have properties"));
        }
        Ok(PyTuple::empty(py))
    }
    #[getter]
    fn dummy_index(&self) -> PyResult<BigInt> {
        self.e.as_sym().and_then(|s| s.dummy_index()).cloned()
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("expression is not a Dummy"))
    }
    #[getter]
    fn is_Number(&self) -> bool {
        self.e.is_number() || matches!(self.e.kind(),Kind::Const(Constant::Infinity|Constant::NegInfinity|Constant::NaN))
    }
    #[getter]
    fn is_Integer(&self) -> bool {
        self.e.as_number().is_some_and(|n|n.is_integer())
    }
    #[getter]
    fn is_Rational(&self) -> bool {
        self.e.is_rational_number()
    }
    #[getter]
    fn __class__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        // ABC protocols consult __class__. Expose the existing public numeric
        // facades only for actual rational nodes; registering Expr itself
        // would incorrectly classify symbols and transcendental expressions.
        if self.is_Tuple() {
            let class = py.import("symox")?.getattr("Tuple")?;
            if class.downcast::<pyo3::types::PyType>().is_ok() {
                return Ok(class);
            }
        }
        if self.e.is_rational_number() {
            let name = if self.is_Integer() { "Integer" } else { "Rational" };
            let class = py.import("symox")?.getattr(name)?;
            if class.downcast::<pyo3::types::PyType>().is_ok() {
                return Ok(class);
            }
        }
        Ok(py.get_type::<PyExprNode>().into_any())
    }
    #[getter]
    fn numerator(&self) -> PyResult<BigInt> {
        self.p()
    }
    #[getter]
    fn denominator(&self) -> PyResult<BigInt> {
        self.q()
    }
    #[getter]
    fn is_Float(&self) -> bool {
        self.e.as_number().is_some_and(Number::is_float)
    }

    #[getter]
    fn _prec(&self) -> PyResult<u32> {
        self.e.as_number().and_then(Number::float_precision)
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("only Float numbers have a bit precision"))
    }
    #[getter]
    fn is_Add(&self) -> bool {
        self.e.is_add()
    }
    #[getter]
    fn is_Mul(&self) -> bool {
        self.e.is_mul()
    }
    #[getter]
    fn is_Pow(&self) -> bool {
        self.e.is_pow()
    }
    #[getter]
    fn is_Function(&self) -> bool {
        matches!(self.e.kind(), Kind::Func(f, _) if !matches!(f, Func::GenericExpr | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _)) && !matches!(f,Func::Undefined(head) if head.is_internal("Tuple")))
    }
    #[getter]
    fn is_UndefinedFunction(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Undefined(f),_) if !f.internal)}
    #[getter]
    fn is_Derivative(&self) -> bool {
        matches!(self.e.kind(), Kind::Derivative(..))
    }
    #[getter]
    fn is_Subs(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Subs,_))}
    #[getter]
    fn is_Limit(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Limit,_))}
    #[getter]
    fn is_Tuple(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Undefined(head),_) if head.is_internal("Tuple"))}
    fn _tuple_len(&self) -> PyResult<usize> {
        match self.e.kind() {
            Kind::Func(Func::Undefined(head), arguments) if head.is_internal("Tuple") => Ok(arguments.len()),
            _ => Err(PyTypeError::new_err("only Tuple expressions have sequence length")),
        }
    }
    fn _tuple_item(&self, py: Python<'_>, index: BigInt) -> PyResult<PyObject> {
        let Kind::Func(Func::Undefined(head), arguments) = self.e.kind() else {
            return Err(PyTypeError::new_err("only Tuple expressions support indexing"));
        };
        if !head.is_internal("Tuple") {
            return Err(PyTypeError::new_err("only Tuple expressions support indexing"));
        }
        let index = if index.is_negative() { index+BigInt::from(arguments.len()) } else { index };
        let argument = index.to_usize().and_then(|index| arguments.get(index))
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("tuple index out of range"))?;
        expr_object(py, argument.clone())
    }
    #[getter]
    fn is_Lambda(&self)->bool{symox::sets::lambda_signature(&self.e).is_some()}
    #[getter]
    fn is_IdentityFunction(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Set(symox::sets::SetOp::IdentityFunction),_))}
    #[getter]
    fn signature(&self)->PyResult<PyExpr>{symox::sets::lambda_signature(&self.e).map(wrap).ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("only Lambda has signature"))}
    #[getter]
    fn nargs(&self)->PyResult<PyExpr>{symox::sets::lambda_arity(&self.e).map(|n|wrap(symox::sets::finite_set(vec![Expr::int(n as i64)]))).ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("only Lambda has nargs"))}
    #[getter]
    fn is_identity(&self)->PyResult<bool>{
        let signature=symox::sets::lambda_signature(&self.e).ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("only Lambda has is_identity"))?;
        Ok(symox::sets::lambda_parts(&self.e).is_some_and(|(_,body)|signature==body))
    }
    #[getter]
    fn _lambda_variables<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        let (variables,_)=symox::sets::lambda_parts(&self.e).ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("not a Lambda"))?;
        PyTuple::new(py,variables.into_iter().map(wrap))
    }
    #[getter]
    fn _lambda_body(&self)->PyResult<PyExpr>{symox::sets::lambda_parts(&self.e).map(|(_,body)|wrap(body)).ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("not a Lambda"))}
    #[getter]
    fn _e(&self)->PyResult<PyExpr>{if self.is_Lambda(){Ok(wrap(self.e.clone()))}else{Err(pyo3::exceptions::PyAttributeError::new_err("not a Lambda"))}}
    #[pyo3(signature=(*args))]
    fn __call__(&self,args:&Bound<'_,PyTuple>)->PyResult<PyExpr>{
        if let Some(p)=symox::ask::as_predicate(&self.e){
            return symox::ask::applied(p,args.iter().map(|a|replacement_expr(&a)).collect::<PyResult<_>>()?).map(wrap).map_err(ask::error);
        }
        logic::apply_lambda_expr(&self.e,args)
    }
    #[getter]
    fn is_Integral(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Integral,_))}
    #[getter]
    fn is_Sum(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Sum,_))}
    #[getter]
    fn is_Product(&self)->bool{matches!(self.e.kind(),Kind::Func(Func::Product,_))}
    #[getter]
    fn function(&self)->PyResult<PyExpr>{
        if let Some(function)=transforms::argument(&self.e,0){return Ok(function);}
        match self.e.kind(){Kind::Func(Func::Integral|Func::Sum|Func::Product|Func::AppliedPredicate,args)=>Ok(wrap(args[0].clone())),_=>Err(pyo3::exceptions::PyAttributeError::new_err("only calculus binders have function"))}
    }
    // integral transforms (transforms.rs)
    #[getter]
    fn function_variable(&self)->PyResult<PyExpr>{transforms::argument(&self.e,1).ok_or_else(||transforms::not_a_transform("function_variable"))}
    #[getter]
    fn transform_variable(&self)->PyResult<PyExpr>{transforms::argument(&self.e,2).ok_or_else(||transforms::not_a_transform("transform_variable"))}
    #[getter]
    fn fundamental_plane(&self,py:Python<'_>)->PyResult<PyObject>{transforms::python_attribute(py,self,"fundamental_plane")}
    #[getter]
    fn fundamental_strip(&self,py:Python<'_>)->PyResult<PyObject>{transforms::python_attribute(py,self,"fundamental_strip")}
    #[getter]
    fn as_integral(&self,py:Python<'_>)->PyResult<PyObject>{transforms::python_attribute(py,self,"as_integral")}
    #[getter]
    fn arguments<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        match symox::ask::as_applied(&self.e){
            Some((_,args))=>PyTuple::new(py,args.iter().map(|e|expr_object(py,e.clone())).collect::<PyResult<Vec<_>>>()?),
            None=>Err(pyo3::exceptions::PyAttributeError::new_err("only applied predicates have arguments")),
        }
    }
    #[getter]
    fn arg(&self)->PyResult<PyExpr>{
        match symox::ask::as_applied(&self.e){
            Some((_,[arg]))=>Ok(wrap(arg.clone())),
            Some(_)=>Err(PyTypeError::new_err("arg is allowed only for unary predicates")),
            None=>Err(pyo3::exceptions::PyAttributeError::new_err("only applied predicates have arg")),
        }
    }
    #[getter]
    fn limits<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        match self.e.kind(){Kind::Func(Func::Integral|Func::Sum|Func::Product,args)=>PyTuple::new(py,args[1..].iter().map(|limit|PyTuple::new(py,limit.args().into_iter().map(wrap))).collect::<PyResult<Vec<_>>>()?),_=>Err(pyo3::exceptions::PyAttributeError::new_err("only calculus binders have limits"))}
    }
    #[getter]
    fn bound_symbols<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{
        if self.is_Lambda(){return Ok(self._lambda_variables(py)?.into_any());}
        match self.e.kind(){Kind::Func(f @ (Func::Integral|Func::Sum|Func::Product),args)=>Ok(PyList::new(py,symox::binders::limits(args).into_iter().filter(|limit|*f!=Func::Integral||limit.len()>1).map(|limit|wrap(limit[0].clone())).collect::<Vec<_>>())?.into_any()),_=>Err(pyo3::exceptions::PyAttributeError::new_err("only binders have bound_symbols"))}
    }
    #[getter]
    fn variable_count<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        match self.e.kind(){Kind::Derivative(_,vars)=>PyTuple::new(py,vars.iter().map(|(v,n)|(wrap(v.clone()),*n))),_=>Err(pyo3::exceptions::PyAttributeError::new_err("only derivatives have variable_count"))}
    }
    #[getter]
    fn derivative_count(&self)->PyResult<u64>{
        match self.e.kind(){Kind::Derivative(_,vars)=>Ok(vars.iter().map(|(_,n)|*n as u64).sum()),_=>Err(pyo3::exceptions::PyAttributeError::new_err("only derivatives have derivative_count"))}
    }
    #[getter]
    fn is_Relational(&self) -> bool {
        matches!(self.e.kind(), Kind::Rel(..))
    }
    #[getter]
    fn is_Piecewise(&self) -> bool {
        matches!(self.e.kind(), Kind::Piecewise(..))
    }
    /// Structural SymPy number classification, including nonfinite constants.
    #[getter]
    fn is_number(&self) -> bool {
        symox::symbolic_number::is_symbolic_number(&self.e)
    }
    #[getter]
    fn is_zero(&self) -> Option<bool> {
        symox::assumptions::is_zero(&self.e)
    }
    #[getter]
    fn is_nonzero(&self) -> Option<bool> {
        symox::assumptions::is_nonzero(&self.e)
    }
    #[getter]
    fn is_integer(&self) -> Option<bool> {
        symox::assumptions::is_integer(&self.e)
    }
    #[getter]
    fn is_rational(&self) -> Option<bool> {
        symox::assumptions::is_rational(&self.e)
    }
    #[getter]
    fn is_real(&self) -> Option<bool> {
        symox::assumptions::is_real(&self.e)
    }
    #[getter]
    fn is_positive(&self) -> Option<bool> {
        symox::assumptions::is_positive(&self.e)
    }
    #[getter]
    fn is_negative(&self) -> Option<bool> {
        symox::assumptions::is_negative(&self.e)
    }
    #[getter]
    fn is_nonnegative(&self) -> Option<bool> {
        symox::assumptions::is_nonnegative(&self.e)
    }
    #[getter]
    fn is_nonpositive(&self) -> Option<bool> {
        symox::assumptions::is_nonpositive(&self.e)
    }
    #[getter]
    fn is_even(&self) -> Option<bool> {
        symox::assumptions::is_even(&self.e)
    }
    #[getter]
    fn is_odd(&self) -> Option<bool> {
        symox::assumptions::is_odd(&self.e)
    }
    #[getter]
    fn is_finite(&self) -> Option<bool> {
        symox::assumptions::is_finite(&self.e)
    }
    #[getter]
    fn is_complex(&self) -> Option<bool> {
        symox::assumptions::fact(&self.e,"complex")
    }
    #[getter]
    fn is_algebraic(&self)->Option<bool>{symox::assumptions::fact(&self.e,"algebraic")}
    #[getter]
    fn is_antihermitian(&self)->Option<bool>{symox::assumptions::fact(&self.e,"antihermitian")}
    #[getter]
    fn is_composite(&self)->Option<bool>{symox::assumptions::fact(&self.e,"composite")}
    #[getter]
    fn is_extended_negative(&self)->Option<bool>{symox::assumptions::fact(&self.e,"extended_negative")}
    #[getter]
    fn is_extended_nonnegative(&self)->Option<bool>{symox::assumptions::fact(&self.e,"extended_nonnegative")}
    #[getter]
    fn is_extended_nonpositive(&self)->Option<bool>{symox::assumptions::fact(&self.e,"extended_nonpositive")}
    #[getter]
    fn is_extended_nonzero(&self)->Option<bool>{symox::assumptions::fact(&self.e,"extended_nonzero")}
    #[getter]
    fn is_extended_positive(&self)->Option<bool>{symox::assumptions::fact(&self.e,"extended_positive")}
    #[getter]
    fn is_extended_real(&self)->Option<bool>{symox::assumptions::fact(&self.e,"extended_real")}
    #[getter]
    fn is_hermitian(&self)->Option<bool>{symox::assumptions::fact(&self.e,"hermitian")}
    #[getter]
    fn is_imaginary(&self)->Option<bool>{symox::assumptions::fact(&self.e,"imaginary")}
    #[getter]
    fn is_infinite(&self)->Option<bool>{symox::assumptions::fact(&self.e,"infinite")}
    #[getter]
    fn is_irrational(&self)->Option<bool>{symox::assumptions::fact(&self.e,"irrational")}
    #[getter]
    fn is_noninteger(&self)->Option<bool>{symox::assumptions::fact(&self.e,"noninteger")}
    #[getter]
    fn is_polar(&self)->Option<bool>{symox::assumptions::fact(&self.e,"polar")}
    #[getter]
    fn is_prime(&self)->Option<bool>{symox::assumptions::fact(&self.e,"prime")}
    #[getter]
    fn is_transcendental(&self)->Option<bool>{symox::assumptions::fact(&self.e,"transcendental")}
    /// The assumptions of a symbol as a dict (SymPy's `assumptions0`).
    #[getter]
    fn assumptions0<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        if let Some(s) = self.e.as_sym() {
            for (name,value) in s.facts().iter(){d.set_item(name,value)?;}
        }
        Ok(d)
    }
    /// Left/right hand sides of a relational.
    #[getter]
    fn lhs(&self) -> PyResult<PyExpr> {
        match self.e.kind() {
            Kind::Rel(_, a, _) => Ok(wrap(a.clone())),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("not a relational")),
        }
    }
    #[getter]
    fn rhs(&self) -> PyResult<PyExpr> {
        match self.e.kind() {
            Kind::Rel(_, _, b) => Ok(wrap(b.clone())),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("not a relational")),
        }
    }

    #[getter]
    fn rel_op(&self)->PyResult<&'static str> {
        match self.e.kind() {
            Kind::Rel(operation,_,_)=>Ok(operation.symbol()),
            _=>Err(pyo3::exceptions::PyAttributeError::new_err("not a relational")),
        }
    }

    // algorithms -------------------------------------------------------

    #[pyo3(signature = (*args))]
    fn diff(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
        Ok(wrap(do_diff(&self.e, args)?))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn subs(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
        let (args,options)=substitution_arguments(args,_kwargs)?;
        let pairs = subs_pairs(&args)?;
        Ok(wrap(checked_subs(args.py(),&self.e,&pairs,substitution_options(options.as_ref())?)?))
    }

    /// `expr.match(pattern)`: dict of wildcard bindings or None.
    #[pyo3(name = "match")]
    fn py_match<'py>(&self, py: Python<'py>, pattern: &Bound<'py, PyAny>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let p = to_expr(pattern)?;
        match symox::pattern::pattern_match(&self.e, &p) {
            None => Ok(None),
            Some(m) => {
                let d = PyDict::new(py);
                for (k, v) in m {
                    d.set_item(wrap(k), wrap(v))?;
                }
                Ok(Some(d))
            }
        }
    }

    /// `expr.replace(pattern, value)`: pattern may contain Wilds, value may
    /// be an expression or a callable receiving the bindings by name; a
    /// function head (e.g. `sin`) may be replaced by another callable.
    fn replace(&self, py: Python<'_>, pattern: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        let pat = to_expr(pattern);
        match pat {
            Ok(p) => {
                if let Ok(v) = to_expr(value) {
                    return Ok(wrap(symox::pattern::replace(&self.e, &p, &v)));
                }
                if !value.is_callable() {
                    return Err(PyTypeError::new_err("replacement must be an expression or a callable"));
                }
                let mut err: Option<PyErr> = None;
                let r = symox::pattern::replace_with(&self.e, &p, &mut |bind, _| {
                    let kw = PyDict::new(py);
                    for (k, v) in bind {
                        let name = k.as_sym().map(|s| s.name().to_string()).unwrap_or_default();
                        if kw.set_item(name, wrap(v.clone())).is_err() {
                            return None;
                        }
                    }
                    match value.call((), Some(&kw)).and_then(|r| to_expr(&r)) {
                        Ok(e) => Some(e),
                        Err(e) => {
                            err = Some(e);
                            None
                        }
                    }
                });
                if let Some(e) = err {
                    return Err(e);
                }
                Ok(wrap(r))
            }
            Err(_) => {
                // function head replacement: replace(sin, cos)
                let name: String = pattern.getattr("__name__")?.extract()?;
                let mut err: Option<PyErr> = None;
                let r = replace_heads(&self.e, &name, &mut |args: Vec<Expr>| {
                    let t = match PyTuple::new(py, args.into_iter().map(wrap)) {
                        Ok(t) => t,
                        Err(e) => {
                            err = Some(e);
                            return None;
                        }
                    };
                    match value.call1(t).and_then(|r| to_expr(&r)) {
                        Ok(e) => Some(e),
                        Err(e) => {
                            err = Some(e);
                            None
                        }
                    }
                });
                if let Some(e) = err {
                    return Err(e);
                }
                Ok(wrap(r))
            }
        }
    }

    /// Sub-expressions matching a pattern.
    fn find<'py>(&self, py: Python<'py>, pattern: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PySet>> {
        let p = to_expr(pattern)?;
        PySet::new(py, symox::pattern::find(&self.e, &p).into_iter().map(wrap))
    }

    fn xreplace(&self, rule: &Bound<'_, PyDict>) -> PyResult<PyExpr> {
        let mut m = FxHashMap::default();
        for (k, v) in rule.iter() {
            m.insert(replacement_expr(&k)?, replacement_expr(&v)?);
        }
        validate_matrix_replacements(rule.py(), &self.e, &m.iter().map(|(k,v)|(k.clone(),v.clone())).collect::<Vec<_>>(), true)?;
        let result=subs::xreplace(&self.e,&m);
        symox::limits::validate_held(&result).map_err(PyNotImplementedError::new_err)?;
        symox::binders::validate_tree(&result).map_err(PyValueError::new_err)?;
        symox::sets::validate_lambda_tree(&result).map_err(logic::bad_signature)?;
        Ok(wrap(result))
    }

    #[pyo3(signature = (deep=true, **hints))]
    fn doit(&self, py: Python<'_>, deep: bool, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
        if symox::transforms::as_held(&self.e).is_some() {
            // integral transforms: SymPy's hints and (result, plane, cond) tuples
            return Ok(py.import("symox._transforms")?.getattr("_doit")?.call((wrap(self.e.clone()), deep), hints)?.unbind());
        }
        if let Some(hints)=hints {for (name,value) in hints.iter(){
            if matches!(self.e.kind(),Kind::Func(Func::MatrixDotProduct(_,_),_))&&name.extract::<String>()?=="expand" {
                // SymPy DotProduct accepts this compatibility argument without
                // changing contraction evaluation (it does not expand factors).
                let _:bool=value.extract()?;
            }else{return Err(PyNotImplementedError::new_err("custom doit hints are not implemented"));}
        }}
        symox::limits::check_doit(&self.e,deep).map_err(|error|PyValueError::new_err(error.to_string()))?;
        Ok(wrap(fx::doit(&self.e, deep)).into_pyobject(py)?.into_any().unbind())
    }

    #[pyo3(signature = (deep=true, **hints))]
    fn expand(&self, py: Python<'_>, deep: bool, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
        expr_object(py, expand_hints(&self.e, deep, hints)?)
    }

    fn trigsimp(&self) -> PyExpr {
        wrap(symox::simp::trigsimp(&self.e))
    }

    fn powsimp(&self) -> PyExpr {
        wrap(symox::simp::powsimp(&self.e))
    }

    fn nsimplify(&self) -> PyExpr {
        wrap(symox::simp::nsimplify(&self.e, &[]))
    }

    #[pyo3(signature = (n=15, subs=None, **_kwargs))]
    fn evalf(
        &self,
        n: i64,
        subs: Option<&Bound<'_, PyDict>>,
        _kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyExpr> {
        Ok(wrap(evalf_expr(&self.e, n, subs, _kwargs)?))
    }

    #[pyo3(signature = (n=15, subs=None, **kwargs))]
    fn n(
        &self,
        n: i64,
        subs: Option<&Bound<'_, PyDict>>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyExpr> {
        self.evalf(n, subs, kwargs)
    }

    fn simplify(&self) -> PyExpr {
        wrap(simplify::simplify(&self.e))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn integrate(&self, py:Python<'_>, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
        integrate_result(py,&self.e,args,_kwargs)
    }

    fn together(&self) -> PyExpr {
        wrap(simplify::together(&self.e))
    }

    fn radsimp(&self) -> PyExpr {
        wrap(simplify::radsimp(&self.e))
    }

    fn copy(&self) -> PyExpr { wrap(self.e.clone()) }

    #[pyo3(signature=(*_gens))]
    fn as_expr(&self, _gens:&Bound<'_,PyTuple>) -> PyExpr { wrap(self.e.clone()) }

    #[pyo3(signature = (deep=true, **hints))]
    fn as_real_imag(&self,deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<Option<(PyExpr,PyExpr)>>{
        let ignore=complex_parts_options(hints,true)?;
        symox::complex_parts::as_real_imag(&self.e,deep,ignore.as_ref())
            .map(|pair|pair.map(|(real,imag)|(wrap(real),wrap(imag))))
            .map_err(|error|PyNotImplementedError::new_err(error.to_string()))
    }

    fn conjugate(&self) -> PyExpr { wrap(fx::conjugate(&self.e)) }
    fn adjoint(&self) -> PyExpr { wrap(fx::conjugate(&self.e)) }
    fn transpose(&self) -> PyExpr { wrap(self.e.clone()) }

    #[pyo3(signature=(*types))]
    fn atoms(&self,py:Python<'_>,types:&Bound<'_,PyTuple>)->PyResult<PyObject>{
        expression_inspection::atoms(py,&self.e,types)
    }

    #[pyo3(signature=(*syms))]
    fn is_polynomial(&self,syms:&Bound<'_,PyTuple>)->PyResult<Option<bool>>{
        expression_inspection::classify(&self.e,syms,symox::inspection::PolynomialKind::Polynomial)
    }

    #[pyo3(signature=(*syms))]
    fn is_rational_function(&self,syms:&Bound<'_,PyTuple>)->PyResult<Option<bool>>{
        expression_inspection::classify(&self.e,syms,symox::inspection::PolynomialKind::Rational)
    }

    #[pyo3(signature=(*syms))]
    fn is_algebraic_expr(&self,syms:&Bound<'_,PyTuple>)->PyResult<Option<bool>>{
        expression_inspection::classify(&self.e,syms,symox::inspection::PolynomialKind::Algebraic)
    }

    #[pyo3(signature = (*patterns))]
    fn has(&self, patterns: &Bound<'_, PyTuple>) -> PyResult<bool> {
        for p in patterns.iter() {
            if let Some(target)=function_head(&p)?{
                fn contains(e:&Expr,target:&symox::expr::UndefinedFunction)->bool{
                    let mut pending=vec![e.clone()];let mut seen=FxHashMap::default();
                    while let Some(node)=pending.pop() {
                        if seen.insert(node.dag_id(),node.clone()).is_some(){continue;}
                        if matches!(node.kind(),Kind::Func(Func::Undefined(function),_) if function==target) {return true;}
                        pending.extend(symox::inspection::scalar_children(&node));
                    }
                    false
                }
                if contains(&self.e,&target){return Ok(true);}continue;
            }
            if traverse::has(&self.e, &to_expr(&p)?) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn as_numer_denom<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (n, d) = simplify::numer_denom(&self.e);
        PyTuple::new(py, [wrap(n), wrap(d)])
    }

    fn as_base_exp<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (b,x)=symox::decompose::as_base_exp(&self.e);
        PyTuple::new(py, [wrap(b), wrap(x)])
    }

    fn class_key(&self,py:Python<'_>)->PyResult<PyObject>{
        presentation_key_object(py,&symox::sortkey::class_key(&self.e))
    }

    #[pyo3(signature=(order=None))]
    fn sort_key(&self,py:Python<'_>,order:Option<&str>)->PyResult<PyObject>{
        validate_sort_order(order)?;
        presentation_key_object(py,&symox::sortkey::sort_key(&self.e))
    }

    #[pyo3(signature=(rational=false))]
    fn as_coeff_Add<'py>(&self,py:Python<'py>,rational:bool)->PyResult<Bound<'py,PyTuple>>{
        let(coefficient,tail)=symox::decompose::as_coeff_add_head(&self.e,rational);
        PyTuple::new(py,[wrap(coefficient),wrap(tail)])
    }

    fn as_powers_dict(&self,py:Python<'_>)->PyResult<PyObject>{
        let result=if self.e.is_constant(Constant::NegInfinity){PyDict::new(py).into_any()}
        else{py.import("collections")?.getattr("defaultdict")?.call1((py.get_type::<pyo3::types::PyInt>(),))?};
        for(base,exponent)in symox::decompose::as_powers_dict(&self.e){result.set_item(expr_object(py,base)?,expr_object(py,exponent)?)?;}
        Ok(result.unbind())
    }

    #[pyo3(signature=(cset=false,warn=true,split_1=true))]
    fn args_cnc(&self,py:Python<'_>,cset:bool,warn:bool,split_1:bool)->PyResult<PyObject>{
        let(commutative,noncommutative)=symox::decompose::args_cnc(&self.e,split_1);
        let nc=PyList::new(py,noncommutative.into_iter().map(|e|expr_object(py,e)).collect::<PyResult<Vec<_>>>()?)?;
        let c=if cset{
            let original=commutative.len();
            let c=PySet::new(py,commutative.into_iter().map(|e|expr_object(py,e)).collect::<PyResult<Vec<_>>>()?)?;
            if warn && c.len()!=original{return Err(PyValueError::new_err("repeated commutative arguments"));}
            c.into_any()
        }else{PyList::new(py,commutative.into_iter().map(|e|expr_object(py,e)).collect::<PyResult<Vec<_>>>()?)?.into_any()};
        Ok(PyList::new(py,[c,nc.into_any()])?.into_any().unbind())
    }

    fn could_extract_minus_sign(&self)->bool{symox::decompose::could_extract_minus_sign(&self.e)}

    #[pyo3(signature=(rational=false))]
    fn as_coeff_Mul<'py>(&self, py: Python<'py>,rational:bool) -> PyResult<Bound<'py, PyTuple>> {
        let(c,r)=symox::decompose::as_coeff_mul_head(&self.e,rational);
        PyTuple::new(py,[wrap(c),wrap(r)])
    }

    #[pyo3(signature=(*deps))]
    fn as_coeff_add<'py>(&self,py:Python<'py>,deps:&Bound<'_,PyTuple>)->PyResult<Bound<'py,PyTuple>>{
        let deps=deps.iter().map(|d|to_expr(&d)).collect::<PyResult<Vec<_>>>()?;
        let(c,tail)=symox::decompose::as_coeff_add(&self.e,&deps);
        PyTuple::new(py,[wrap(c).into_pyobject(py)?.into_any(),PyTuple::new(py,tail.into_iter().map(wrap))?.into_any()])
    }

    #[pyo3(signature=(*deps,rational=true,**kwargs))]
    fn as_coeff_mul<'py>(&self,py:Python<'py>,deps:&Bound<'_,PyTuple>,rational:bool,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Bound<'py,PyTuple>>{
        if kwargs.is_some_and(|kw|!kw.is_empty()){return Err(PyNotImplementedError::new_err("unsupported as_coeff_mul options"));}
        let deps=deps.iter().map(|d|to_expr(&d)).collect::<PyResult<Vec<_>>>()?;
        let(c,tail)=symox::decompose::as_coeff_mul(&self.e,&deps,rational);
        PyTuple::new(py,[wrap(c).into_pyobject(py)?.into_any(),PyTuple::new(py,tail.into_iter().map(wrap))?.into_any()])
    }

    #[pyo3(signature=(*deps,**hints))]
    fn as_independent<'py>(&self,py:Python<'py>,deps:&Bound<'_,PyTuple>,hints:Option<&Bound<'_,PyDict>>)->PyResult<Bound<'py,PyTuple>>{
        let mut as_add=None;
        if let Some(hints)=hints{for(key,value)in hints.iter(){let key:String=key.extract()?;if key=="as_Add"{as_add=Some(value.is_truthy()?);}else{return Err(PyNotImplementedError::new_err(format!("as_independent hint '{key}' is not implemented")));}}}
        let deps=deps.iter().map(|d|to_expr(&d)).collect::<PyResult<Vec<_>>>()?;
        let(i,d)=symox::decompose::as_independent(&self.e,&deps,as_add);
        PyTuple::new(py,[wrap(i),wrap(d)])
    }

    #[pyo3(signature=(radical=false,clear=true))]
    fn as_content_primitive<'py>(&self,py:Python<'py>,radical:bool,clear:bool)->PyResult<Bound<'py,PyTuple>>{
        let(c,p)=symox::decompose::as_content_primitive(&self.e,radical,clear).map_err(PyNotImplementedError::new_err)?;
        PyTuple::new(py,[wrap(Expr::num(c)),wrap(p)])
    }

    fn as_coeff_exponent<'py>(&self, py: Python<'py>, x: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        let x = to_expr(x)?;
        let (coefficient, exponent) = if self.e == x {
            (Expr::one(), Expr::one())
        } else {
            match self.e.kind() {
                Kind::Pow(base, exponent) if *base == x => (Expr::one(), exponent.clone()),
                Kind::Mul(c, factors) => {
                    let mut exponent = Expr::zero();
                    let mut rest = vec![Expr::num(c.clone())];
                    for (base, power) in factors {
                        if *base == x {
                            exponent = ops::add2(&exponent, power);
                        } else {
                            rest.push(ops::mul_factor(base, power));
                        }
                    }
                    (ops::mul(&rest), exponent)
                }
                _ => (self.e.clone(), Expr::zero()),
            }
        };
        PyTuple::new(py, [wrap(coefficient), wrap(exponent)])
    }

    fn as_two_terms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let args = python_args(&self.e);
        if args.len() < 2 || !(self.e.is_add() || self.e.is_mul()) {
            return Err(PyValueError::new_err("expected an Add or Mul with at least two arguments"));
        }
        let rest = if self.e.is_add() { ops::add(&args[1..]) } else { ops::mul(&args[1..]) };
        PyTuple::new(py, [wrap(args[0].clone()), wrap(rest)])
    }

    #[getter]
    fn base(&self) -> PyResult<PyExpr> {
        match self.e.kind() {
            Kind::Pow(base, _) => Ok(wrap(base.clone())),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("expression has no base")),
        }
    }

    #[getter]
    fn exp(&self) -> PyResult<PyExpr> {
        match self.e.kind() {
            Kind::Pow(_, exponent) => Ok(wrap(exponent.clone())),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("expression has no exp")),
        }
    }

    #[getter]
    fn is_commutative(&self) -> Option<bool> {
        symox::decompose::commutativity(&self.e)
    }

    fn count_ops(&self) -> BigInt { operation_count(&self.e) }

    fn as_ordered_factors(&self) -> Vec<PyExpr> {
        symox::sortkey::ordered_factors(&self.e).into_iter().map(wrap).collect()
    }

    fn as_ordered_terms(&self) -> Vec<PyExpr> {
        let terms=if self.e.may_have_order(){printer::ordered_terms(&self.e)}else{symox::sortkey::ordered_terms(&self.e)};
        terms.into_iter().map(wrap).collect()
    }

    fn as_coefficients_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        let push = |t: Expr, k: Number| -> PyResult<()> {
            let key = wrap(t);
            let prev = d.get_item(key.clone())?;
            let val = match prev {
                Some(p) => ops::add2(&p.downcast::<crate::PyExprNode>()?.get().e, &Expr::num(k)),
                None => Expr::num(k),
            };
            d.set_item(key, wrap(val))
        };
        match self.e.kind() {
            Kind::Add(c, terms) => {
                for (t, k) in terms {
                    push(t.clone(), k.clone())?;
                }
                if !c.is_zero() {
                    push(Expr::one(), c.clone())?;
                }
            }
            _ => {
                let (c, r) = ops::as_coeff_mul(&self.e);
                push(r, c)?;
            }
        }
        Ok(d)
    }

    /// Coefficient of `x**n` in the expanded expression.
    #[pyo3(signature = (x, n=None))]
    fn coeff(&self, x: &Bound<'_, PyAny>, n: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        let x = to_expr(x)?;
        let n = match n {
            Some(n) => to_expr(n)?,
            None => Expr::one(),
        };
        Ok(wrap(coeff_of(&self.e, &x, &n)))
    }

    /// Convert to a SymPy expression.
    fn to_sympy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        expr_to_sympy(py, &self.e)
    }

    fn _sympy_<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.to_sympy(py)
    }

    fn __getnewargs__(&self) -> (String,) {
        (printer::srepr(&self.e),)
    }
}

/// Replace applications of the function named `name` bottom-up.
fn replace_heads(e: &Expr, name: &str, f: &mut dyn FnMut(Vec<Expr>) -> Option<Expr>) -> Expr {
    let rebuilt = if e.is_atom() {
        e.clone()
    } else {
        traverse::map_children(e, &mut |c| replace_heads(c, name, f))
    };
    if let Kind::Func(g, args) = rebuilt.kind() {
        if g.name() == name {
            if let Some(r) = f(args.clone()) {
                return r;
            }
        }
    }
    rebuilt
}

/// Coefficient of `x**n` in `e` (SymPy `Expr.coeff` semantics: the sum
/// of the cofactors of terms containing exactly `x**n` and nothing else
/// depending on `x`; `n = 0` gives the terms independent of `x`).
fn coeff_of(e: &Expr, x: &Expr, n: &Expr) -> Expr {
    let terms: Vec<Expr> = match e.kind() {
        Kind::Add(..) => e.args(),
        _ => vec![e.clone()],
    };
    let mut out = Vec::new();
    if n.is_zero() {
        for t in &terms {
            if !traverse::has(t, x) {
                out.push(t.clone());
            }
        }
        return ops::add(&out);
    }
    let target = ops::pow(x, n);
    for t in &terms {
        if *t == target {
            out.push(Expr::one());
            continue;
        }
        let factors = match t.kind() {
            Kind::Mul(..) => t.args(),
            _ => continue,
        };
        if let Some(i) = factors.iter().position(|f| *f == target) {
            let rest: Vec<Expr> = factors
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, f)| f.clone())
                .collect();
            let rest = ops::mul(&rest);
            if !traverse::has(&rest, x) {
                out.push(rest);
            }
        }
    }
    ops::add(&out)
}

/// The head (`func`) of an expression: calling it with arguments rebuilds a
/// node of the same kind.
#[pyclass(name = "Head", module = "symox._core", frozen)]
pub struct Head {
    kind: HeadKind,
}

#[derive(PartialEq, Eq, Hash)]
enum HeadKind {
    Add,
    Mul,
    Pow,
    Func(Func),
    Rel(RelOp),
    Atom(Expr),
    Derivative,
    Piecewise,
}

impl Head {
    fn of(e: &Expr) -> Head {
        let kind = match e.kind() {
            Kind::Add(..) => HeadKind::Add,
            Kind::Mul(..) => HeadKind::Mul,
            Kind::Pow(..) => HeadKind::Pow,
            Kind::Func(f, _) => HeadKind::Func(f.clone()),
            Kind::Rel(op, ..) => HeadKind::Rel(*op),
            Kind::Derivative(..) => HeadKind::Derivative,
            Kind::Piecewise(..) => HeadKind::Piecewise,
            _ => HeadKind::Atom(e.clone()),
        };
        Head { kind }
    }
}

#[pymethods]
impl Head {
    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: pyo3::basic::CompareOp) -> PyResult<PyObject> {
        let py = other.py();
        let equal = other.extract::<PyRef<'_, Self>>().is_ok_and(|other| self.kind == other.kind);
        Ok(match op {
            pyo3::basic::CompareOp::Eq => equal.into_pyobject(py)?.to_owned().into_any().unbind(),
            pyo3::basic::CompareOp::Ne => (!equal).into_pyobject(py)?.to_owned().into_any().unbind(),
            _ => py.NotImplemented(),
        })
    }

    fn __hash__(&self) -> isize {
        let mut hasher = DefaultHasher::new();
        self.kind.hash(&mut hasher);
        let hash = hasher.finish() as isize;
        if hash == -1 { -2 } else { hash }
    }

    #[pyo3(signature = (*args))]
    fn __call__(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
        let v: Vec<Expr> = args.iter().map(|a| to_expr(&a)).collect::<PyResult<_>>()?;
        Ok(wrap(match &self.kind {
            HeadKind::Add => ops::add(&v),
            HeadKind::Mul => ops::mul(&v),
            HeadKind::Pow => ops::pow(&v[0], &v[1]),
            HeadKind::Func(f) => fx::func(f.clone(), v),
            HeadKind::Rel(op) => relational::rel(*op, &v[0], &v[1]),
            HeadKind::Atom(e) => e.clone(),
            HeadKind::Derivative => {
                let arguments=PyTuple::new(args.py(),args.iter().skip(1))?;
                let vars=diff_vars(&arguments)?;
                diff::derivative_unevaluated(&v[0], vars)
            }
            HeadKind::Piecewise => {
                let pairs = v.chunks(2).map(|c| (c[0].clone(), c[1].clone())).collect();
                relational::piecewise(pairs)
            }
        }))
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            HeadKind::Add => "Add".into(),
            HeadKind::Mul => "Mul".into(),
            HeadKind::Pow => "Pow".into(),
            HeadKind::Func(f) => f.name().into(),
            HeadKind::Rel(op) => format!("Rel({})", op.symbol()),
            HeadKind::Atom(e) => e.class_name().to_string(),
            HeadKind::Derivative => "Derivative".into(),
            HeadKind::Piecewise => "Piecewise".into(),
        }
    }

    #[getter]
    fn __name__(&self) -> String {
        self.__repr__()
    }
}

// ----------------------------------------------------------------------
// module level functions
// ----------------------------------------------------------------------

/// Build ordinary or explicitly held arithmetic without a Python expression
/// hierarchy: every returned object still uses the Rust expression engine.
#[pyfunction]
#[pyo3(name = "_node", signature = (name, args, evaluate=true))]
fn py_node(name: &str, args: &Bound<'_, PyTuple>, evaluate: bool) -> PyResult<PyExpr> {
    let v: Vec<Expr> = args.iter().map(|a| to_expr(&a)).collect::<PyResult<_>>()?;
    Ok(wrap(construct_node(name, v, evaluate)?))
}

fn construct_node(name: &str, v: Vec<Expr>, evaluate: bool) -> PyResult<Expr> {
    let evaluate = evaluate || (matches!(name, "Add" | "Mul") && v.len() < 2);
    if !evaluate && matches!(name, "Add" | "Mul") && !v.iter().any(|a|a.has_held_arithmetic()||a.has_matrix()) {
        // Holding evaluation is observable through the public argument tree,
        // not through Rust's choice of coefficient/factor storage. Reuse the
        // canonical representation when it already has the requested tree.
        let canonical = match name {
            "Add" => Some(ops::add(&v)),
            "Mul" => Some(ops::mul(&v)),
            _ => None,
        };
        if let Some(canonical) = canonical {
            if python_args(&canonical) == v { return Ok(canonical); }
        }
    }
    Ok(match name {
        "Add" if evaluate => ops::add(&v),
        "Mul" if evaluate => ops::mul(&v),
        "Add" => Expr::from_unevaluated_kind(Kind::Add(Number::Int(0), v.into_iter().map(|a| (a, Number::Int(1))).collect())),
        "Mul" => Expr::from_unevaluated_kind(Kind::Mul(Number::Int(1), v.into_iter().map(|a| (a, Expr::one())).collect())),
        "Pow" => {
            if v.len() != 2 { return Err(PyTypeError::new_err("Pow expects exactly two arguments")); }
            // evaluate=False must never probe the evaluated candidate: even
            // a tiny input such as Pow(2, 10**9) can allocate enormous integers.
            if evaluate { ops::pow(&v[0], &v[1]) } else { Expr::from_unevaluated_kind(Kind::Pow(v[0].clone(), v[1].clone())) }
        }
        "UnevaluatedExpr" => {
            if v.len() != 1 { return Err(PyTypeError::new_err("UnevaluatedExpr expects exactly one argument")); }
            Expr::from_kind(Kind::Func(Func::Undefined(name.into()), v))
        }
        _ => return Err(PyValueError::new_err(format!("unknown arithmetic node {name}"))),
    })
}

#[pyfunction]
fn _restore_arithmetic(name: &str, coefficient: &Bound<'_, PyAny>, children: Vec<(PyExpr, PyExpr)>, canonical: bool) -> PyResult<PyExpr> {
    let coefficient = to_expr(coefficient)?.as_number().cloned()
        .ok_or_else(|| PyTypeError::new_err("arithmetic coefficient must be a number"))?;
    let children: Vec<_> = children.into_iter().map(|(a, b)| (a.e, b.e)).collect();
    let kind = match name {
        "Add" => Kind::Add(coefficient, children.into_iter().map(|(term, coefficient)| {
            coefficient.as_number().cloned().map(|c| (term, c))
                .ok_or_else(|| PyTypeError::new_err("additive coefficient must be a number"))
        }).collect::<PyResult<_>>()?),
        "Mul" => Kind::Mul(coefficient, children),
        "Pow" if children.len() == 1 => Kind::Pow(children[0].0.clone(), children[0].1.clone()),
        _ => return Err(PyValueError::new_err("invalid arithmetic pickle payload")),
    };
    Ok(wrap(if canonical { Expr::from_kind(kind) } else { Expr::from_unevaluated_kind(kind) }))
}

#[pyfunction]
fn _restore_logic(name: &str, children: Vec<PyExpr>) -> PyResult<PyExpr> {
    let op = symox::logic::LogicOp::from_name(name)
        .ok_or_else(|| PyValueError::new_err("Invalid Boolean connective"))?;
    Ok(wrap(Expr::from_kind(Kind::Func(Func::Logic(op), children.into_iter().map(|child| child.e).collect()))))
}

pub(crate) fn validate_sort_order(order:Option<&str>)->PyResult<()>{
    if order.is_some_and(|order|order!="lex"){
        return Err(PyNotImplementedError::new_err("Sort keys currently support the default lexicographic order only"));
    }
    Ok(())
}

pub(crate) fn presentation_key_object(py:Python<'_>,key:&symox::sortkey::Key)->PyResult<PyObject>{
    fn convert(py: Python<'_>, key: &symox::sortkey::Key, tuples: &mut Option<FxHashMap<usize, PyObject>>) -> PyResult<PyObject> {
        use symox::sortkey::Key;
        Ok(match key {
            Key::Int(n) => (*n).into_pyobject(py)?.into_any().unbind(),
            Key::Index(n) => n.clone().into_pyobject(py)?.into_any().unbind(),
            Key::Str(s) => s.as_str().into_pyobject(py)?.into_any().unbind(),
            Key::Num(n) => wrap(Expr::num(n.clone())).into_pyobject(py)?.into_any().unbind(),
            Key::SpecialNumber(c) => wrap(Expr::constant(*c)).into_pyobject(py)?.into_any().unbind(),
            Key::Tup(keys) => PyTuple::new(py, keys.iter().map(|key| convert(py, key, tuples)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind(),
            Key::Shared(keys) => {
                let id = keys.as_ptr() as usize;
                if let Some(value) = tuples.as_ref().and_then(|cache| cache.get(&id)) { return Ok(value.clone_ref(py)); }
                let children = keys.iter().map(|key| convert(py, key, tuples)).collect::<PyResult<Vec<_>>>()?;
                let value = PyTuple::new(py, children)?.into_any().unbind();
                if let Some(cache) = tuples { cache.insert(id, value.clone_ref(py)); }
                value
            },
        })
    }
    let mut tuples = Some(FxHashMap::default());
    convert(py,key,&mut tuples)
}

#[pyfunction]
#[pyo3(name = "_default_sort_key")]
fn py_default_sort_key(py: Python<'_>, obj: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let expression = to_expr(obj)?;
    let key = symox::sortkey::sort_key(&expression);
    presentation_key_object(py, &key)
}

#[pyfunction]
#[pyo3(name = "Symbol", signature = (name, **assumptions))]
fn py_symbol(name: &str, assumptions: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let Some(d) = assumptions else {
        return Ok(wrap(Expr::symbol(name)));
    };
    let mut kw: Vec<(String, bool)> = Vec::new();
    for (k, v) in d.iter() {
        let k: String = k.extract()?;
        if v.is_none() {
            if k=="commutative" {return Err(PyValueError::new_err("Symbol commutativity must be True or False"));}
            continue;
        }
        if k=="commutative"&&!v.eq(true)?&&!v.eq(false)? {return Err(PyValueError::new_err("Symbol commutativity must be True or False"));}
        kw.push((k, v.is_truthy()?));
    }
    let a = symox::assumption_facts::FactSet::from_kwargs(kw.iter().map(|(k, v)| (k.as_str(), *v)))
        .map_err(PyValueError::new_err)?;
    Ok(wrap(Expr::from_kind(Kind::Sym(symox::Symbol::with_facts(name,a).map_err(PyValueError::new_err)?))))
}

#[pyfunction]
#[pyo3(name = "Dummy", signature = (name=None, dummy_index=None, **assumptions))]
fn py_dummy(name: Option<&str>, dummy_index: Option<BigInt>, assumptions: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let given = py_symbol("", assumptions)?.e.as_sym().unwrap().given_facts();
    let symbol = match (name, dummy_index) {
        (Some(name), Some(index)) => symox::Symbol::dummy_with_facts(name, given, Some(index)),
        (Some(name), None) => symox::Symbol::dummy_with_facts(name, given, None),
        (None, None) => symox::Symbol::anonymous_dummy_with_facts(given),
        (None, Some(_)) => return Err(PyValueError::new_err("an explicit dummy_index requires a name")),
    }.map_err(PyValueError::new_err)?;
    Ok(wrap(Expr::from_kind(Kind::Sym(symbol))))
}

#[pyfunction]
#[pyo3(name = "Wild", signature = (name, exclude=None, properties=None, **assumptions))]
fn py_wild(name: &str, exclude: Option<&Bound<'_, PyAny>>, properties: Option<&Bound<'_, PyAny>>, assumptions:Option<&Bound<'_,PyDict>>) -> PyResult<PyExpr> {
    if let Some(properties)=properties {
        if properties.try_iter()?.next().transpose()?.is_some() {
            return Err(PyNotImplementedError::new_err("Wildcard predicate properties are not implemented"));
        }
    }
    let mut ex = Vec::new();
    if let Some(it) = exclude {
        for x in it.try_iter()? {
            ex.push(to_expr(&x?)?);
        }
    }
    let given=py_symbol("",assumptions)?.e.as_sym().unwrap().given_facts();
    Ok(wrap(Expr::from_kind(Kind::Sym(symox::Symbol::wild_with_facts(name,ex,given).map_err(PyValueError::new_err)?))))
}

#[pyfunction]
fn sympify(obj: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    expr_object(obj.py(), to_expr(obj)?)
}

#[pyfunction]
#[pyo3(name = "Integer")]
fn py_integer(obj: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let e = to_expr(obj)?;
    if !e.as_number().is_some_and(|n|n.is_integer()) {
        return Err(PyValueError::new_err(format!("{e} is not an integer")));
    }
    expr_object(obj.py(), e)
}

#[pyfunction]
#[pyo3(name = "Rational", signature = (p, q=None))]
fn py_rational(p: &Bound<'_, PyAny>, q: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
    let py = p.py();
    let p = to_expr(p)?;
    let r = match q {
        Some(q) => {
            let q = to_expr(q)?;
            if q.is_zero() {
                return Err(PyZeroDivisionError::new_err("Rational with zero denominator"));
            }
            ops::div(&p, &q)
        }
        None => p,
    };
    expr_object(py, r)
}

fn float_mpf_parts(payload:&Bound<'_,PyAny>,precision:u32)->PyResult<Expr>{
    let parts=payload.downcast::<PyTuple>().map_err(|_|PyValueError::new_err("malformed mpf tuple"))?;
    if parts.len()!=3&&parts.len()!=4{return Err(PyValueError::new_err("mpf tuples need three or four entries"));}
    let sign=parts.get_item(0)?.extract::<u8>()?;
    if sign>1{return Err(PyValueError::new_err("mpf sign must be zero or one"));}
    let source=parts.get_item(1)?;
    let hexadecimal=source.is_instance_of::<PyString>();
    let mut mantissa=if hexadecimal{
        let text=source.extract::<String>()?;
        BigInt::parse_bytes(text.trim_start_matches("0x").trim_end_matches('L').as_bytes(),16)
            .ok_or_else(||PyValueError::new_err("invalid hexadecimal mpf mantissa"))?
    }else{source.extract::<BigInt>()?};
    if mantissa.is_negative(){return Err(PyValueError::new_err("mpf mantissa must be nonnegative"));}
    let exponent=parts.get_item(2)?.extract::<i64>()?;
    if parts.len()==4&&mantissa.is_zero(){
        let count=parts.get_item(3)?.extract::<i64>()?;
        match (sign,exponent,count){
            (0,-456,-2)=>return Ok(Expr::oo()),
            (1,-789,-3)=>return Ok(Expr::constant(Constant::NegInfinity)),
            (0,-123,-1)=>return Ok(Expr::nan()),
            _=>{},
        }
    }
    if mantissa.is_zero()&&!hexadecimal{return Ok(Expr::zero());}
    if sign==1{mantissa=-mantissa;}
    let value=symox::bigfloat::BigFloat{man:mantissa,exp:exponent}.checked_round_to(precision)
        .ok_or_else(||PyValueError::new_err("mpf value exceeds the native binary exponent counter"))?;
    Ok(Expr::num(Number::mp_float(value,precision)))
}

#[pyfunction]
#[pyo3(name = "Float", signature = (obj, dps=None, *, precision=None))]
fn py_float(obj: &Bound<'_, PyAny>, dps: Option<i64>, precision: Option<u32>) -> PyResult<PyExpr> {
    if dps.is_some() && precision.is_some() { return Err(PyValueError::new_err("specify decimal digits or bit precision, not both")); }
    let mut p = match (dps, precision) {
        (Some(n), _) => decimal_precision(n)?,
        (_, Some(p)) => bit_precision(p)?,
        _ => 53,
    };
    let ty = obj.get_type();
    let decimal = ty.name().is_ok_and(|name| name == "Decimal")
        && ty.module().is_ok_and(|module| module == "decimal");
    if obj.is_instance_of::<PyString>() || decimal {
        let source = if decimal { obj.str()?.to_string() } else { obj.extract::<String>()? };
        let special = match source.trim() {
            "inf" | "+inf" | "Infinity" | "+Infinity" => Some(Expr::oo()),
            "-inf" | "-Infinity" => Some(Expr::constant(Constant::NegInfinity)),
            "nan" | "NaN" => Some(Expr::nan()),
            _ => None,
        };
        if let Some(e) = special { return Ok(wrap(e)); }
        let value = symox::bigfloat::BigFloat::from_decimal(&source, p).map_err(PyValueError::new_err)?;
        return Ok(wrap(Expr::num(Number::mp_float(value, p))));
    }
    if obj.downcast::<PyTuple>().is_ok(){return Ok(wrap(float_mpf_parts(obj,p)?));}
    if obj.downcast::<crate::PyExprNode>().is_err()&&!obj.is_instance_of::<PyInt>()&&!obj.is_instance_of::<PyFloat>(){
        match obj.getattr("_mpf_"){
            Ok(parts)=>{
                if dps.is_none()&&precision.is_none(){
                    let mpmath=obj.py().import("sys")?.getattr("modules")?.call_method1("get",("mpmath",))?;
                    if !mpmath.is_none()&&obj.is_instance(&mpmath.getattr("mpf")?)?{
                        p=bit_precision(obj.getattr("context")?.getattr("prec")?.extract()?)?;
                    }
                }
                return Ok(wrap(float_mpf_parts(&parts,p)?));
            },
            Err(error) if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(obj.py())=>{},
            Err(error)=>return Err(error),
        }
    }
    let e = to_expr(obj)?;
    if matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN)) { return Ok(wrap(e)); }
    if let Some(number) = e.as_number() {
        if dps.is_none() && precision.is_none() && number.is_float() { return Ok(wrap(e)); }
        return Ok(wrap(Expr::num(Number::mp_float(number.to_bigfloat(p), p))));
    }
    match symox::evalf_mp::eval_mp(&e, p) {
        Some(value) if value.is_real() => Ok(wrap(Expr::num(Number::mp_float(value.re, p)))),
        _ => Err(PyTypeError::new_err("Float expects a finite real number")),
    }
}

#[pyfunction]
#[pyo3(name = "srepr")]
fn py_srepr(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    Ok(printer::srepr(&to_expr(obj)?))
}

#[pyfunction]
fn _materialize_matrix_scalars(obj: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let e = to_expr(obj)?;
    if !e.has_matrix() { return Ok(wrap(e)); }
    fx::materialize_matrix_scalars(&e).map(wrap).map_err(|error| matrices::err(obj.py(),error))
}

/// Differentiate an expression with respect to the requested variables.
#[pyfunction]
#[pyo3(name = "diff", signature = (f, *args))]
fn py_diff(f: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    if let Ok(matrix) = f.downcast::<matrixexpr::PyMatrixExpr>() {
        return Ok(matrix.get().differentiate(args)?.into_pyobject(f.py())?.into_any().unbind());
    }
    expr_object(f.py(), do_diff(&to_expr(f)?, args)?)
}

#[pyfunction]
#[pyo3(name = "Derivative", signature = (f, *args, evaluate=false))]
fn py_derivative(f: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>,evaluate:bool) -> PyResult<PyExpr> {
    if f.downcast::<matrixexpr::PyMatrixExpr>().is_ok() {
        return Err(PyNotImplementedError::new_err("held matrix derivatives are not implemented; use diff with a scalar variable"));
    }
    let expression=to_expr(f)?;
    let mut vars = diff_vars(args)?;
    if vars.is_empty(){let symbols=traverse::free_symbols(&expression);match symbols.len(){0=>return Ok(wrap(Expr::zero())),1=>vars.push((symbols[0].clone(),1)),_=>return Err(PyValueError::new_err("specify the differentiation variables"))}}
    if vars.iter().any(|(v,_)|!valid_diff_variable(v)){return Err(PyValueError::new_err("cannot differentiate with respect to this expression"));}
    Ok(wrap(if evaluate{diff::diff_many(&expression,&vars)}else{diff::derivative_unevaluated(&expression,vars)}))
}

#[pyfunction]
fn _undefined_function(name:&str,args:&Bound<'_,PyTuple>,assumptions:&Bound<'_,PyDict>,nargs:Option<Vec<usize>>)->PyResult<PyExpr>{
    let spec=function_spec(name,assumptions,nargs)?;
    if !spec.accepts(args.len()){return Err(PyTypeError::new_err(format!("{name} accepts argument counts {:?}, got {}",spec.nargs.as_ref().unwrap(),args.len())));}
    Ok(wrap(fx::func(Func::Undefined(spec),args.iter().map(|arg|to_expr(&arg)).collect::<PyResult<Vec<_>>>()?)))
}

#[pyfunction]
fn _subs(expr:&Bound<'_,PyAny>,variables:&Bound<'_,PyAny>,points:&Bound<'_,PyAny>)->PyResult<PyExpr>{
    fn values(value:&Bound<'_,PyAny>)->PyResult<Vec<Expr>>{
        if let Ok(e)=to_expr(value){return Ok(match e.kind(){Kind::Func(Func::Undefined(f),args) if f.is_internal("Tuple")=>args.clone(),_=>vec![e]});}
        value.try_iter()?.map(|v|to_expr(&v?)).collect()
    }
    let variables=values(variables)?;let points=values(points)?;
    if variables.len()!=points.len(){return Err(PyValueError::new_err("Subs variables and points must have the same length"));}
    if variables.iter().any(|v|!v.is_symbol()){return Err(PyNotImplementedError::new_err("Subs currently supports symbol variables"));}
    if variables.iter().enumerate().any(|(i,v)|variables[..i].contains(v)){return Err(PyValueError::new_err("Subs variables must be distinct"));}
    Ok(wrap(subs::held_subs(&to_expr(expr)?,variables,points)))
}

/// `integrate(f, x)`, `integrate(f, (x, a, b))`, several variables in turn.
fn integration_limit(value:&Bound<'_,PyAny>)->PyResult<Option<Vec<Expr>>>{
    if value.is_instance_of::<PyTuple>()||value.is_instance_of::<PyList>(){
        return Ok(Some(value.try_iter()?.map(|item|to_expr(&item?)).collect::<PyResult<Vec<_>>>()?));
    }
    if let Ok(value)=value.extract::<PyRef<'_, crate::PyExprNode>>(){
        if let Kind::Func(Func::Undefined(name),arguments)=value.e.kind(){
            if name.is_internal("Tuple"){return Ok(Some(arguments.clone()));}
        }
    }
    Ok(None)
}

fn gamma_integral_conditions(f:&Expr,x:&Expr,lo:&Expr,hi:&Expr)->Option<(Expr,Expr)>{
    symox::integrate::gamma_moment_conditions(f,x,lo,hi).or_else(||{
        if *lo==Expr::oo()&&(hi.is_zero()||*hi==Expr::neg_oo()){
            symox::integrate::gamma_moment_conditions(f,x,hi,lo).map(|(value,condition)|(ops::neg(&value),condition))
        }else{None}
    })
}

fn do_integrate(f: &Expr, args: &Bound<'_, PyTuple>, conditions:IntegrationConditions) -> PyResult<Expr> {
    let mut r = f.clone();
    if args.is_empty() {
        let fs = traverse::free_symbols(f);
        if fs.len() != 1 {
            return Err(PyValueError::new_err("specify the integration variable"));
        }
        return Ok(symox::integrate::integrate(f, &fs[0]));
    }
    for a in args.iter() {
        if let Some(t) = integration_limit(&a)? {
            if !(1..=3).contains(&t.len()){return Err(PyValueError::new_err("expected (x,), (x,b), or (x,a,b)"));}
            let x = t[0].clone();
            if !x.is_symbol(){return Err(PyValueError::new_err("integration variable must be a symbol"));}
            match t.len() {
                1 => r = symox::integrate::integrate(&r, &x),
                2 => {let upper=t[1].clone();
                    let held=symox::binders::held(Func::Integral,&r,vec![vec![x,upper]]).map_err(PyValueError::new_err)?;
                    r=fx::doit(&held,true);
                },
                3 => {
                    let lo = &t[1];let hi=&t[2];
                    r = if conditions==IntegrationConditions::None{
                        gamma_integral_conditions(&r,&x,lo,hi).map(|(value,_)|value)
                            .unwrap_or_else(||symox::integrate::integrate_definite(&r,&x,lo,hi))
                    }else{symox::integrate::integrate_definite(&r,&x,lo,hi)};
                }
                _ => return Err(PyValueError::new_err("expected (x,) or (x, a, b)")),
            }
        } else {
            let variable=to_expr(&a)?;
            if !variable.is_symbol(){return Err(PyValueError::new_err("integration variable must be a symbol"));}
            r = symox::integrate::integrate(&r, &variable);
        }
    }
    Ok(r)
}

fn integrate_result(py:Python<'_>,f:&Expr,args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<PyObject>{
    let conditions=integrate_options(kwargs)?;
    if conditions==IntegrationConditions::Separate&&args.len()>1{
        return Err(PyValueError::new_err("conds=separate not supported in multiple integrals"));
    }
    let mut infinite=false;
    if conditions==IntegrationConditions::Separate&&args.len()==1{
        if let Some(limit)=integration_limit(&args.get_item(0)?)?{
            if limit.len()==3{
                if !limit[0].is_symbol(){return Err(PyValueError::new_err("integration variable must be a symbol"));}
                if let Some((value,condition))=gamma_integral_conditions(f,&limit[0],&limit[1],&limit[2]){
                    return Ok((wrap(value),wrap(condition)).into_pyobject(py)?.into_any().unbind());
                }
                infinite=limit[1..].iter().any(|bound|matches!(bound.kind(),Kind::Const(c) if c.is_infinite()));
            }
        }
    }
    let value=do_integrate(f,args,conditions)?;
    if infinite&&!f.is_zero()&&symox::assumptions::is_finite(&value)==Some(true){
        let mut held=false;
        traverse::preorder(&value,&mut |node|{if matches!(node.kind(),Kind::Func(Func::Integral,_)){held=true;}});
        if !held{return Ok((wrap(value),wrap(Expr::boolean(true))).into_pyobject(py)?.into_any().unbind());}
    }
    Ok(wrap(value).into_pyobject(py)?.into_any().unbind())
}

#[pyfunction]
#[pyo3(name = "summation", signature = (f, *limits))]
fn py_summation(f: &Bound<'_, PyAny>, limits: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    if limits.is_empty(){return Err(PyValueError::new_err("summation requires (k,a,b) limits"));}
    let mut r = to_expr(f)?;
    for l in limits.iter() {
        let t = l
            .downcast::<PyTuple>()
            .map_err(|_| PyTypeError::new_err("summation expects (k, a, b) tuples"))?;
        if t.len() != 3 {
            return Err(PyValueError::new_err("summation expects (k, a, b) tuples"));
        }
        let (k, a, b) = (to_expr(&t.get_item(0)?)?, to_expr(&t.get_item(1)?)?, to_expr(&t.get_item(2)?)?);
        if !k.is_symbol(){return Err(PyValueError::new_err("summation variable must be a symbol"));}
        r = symox::summation::summation(&r, &k, &a, &b);
    }
    Ok(wrap(r))
}

#[pyfunction]
#[pyo3(name = "bernoulli")]
fn py_bernoulli(n: usize) -> PyExpr {
    wrap(symox::summation::bernoulli(n))
}

/// Compute an indefinite or definite symbolic integral using native algorithms.
#[pyfunction]
#[pyo3(name = "integrate", signature = (f, *args, **_kwargs))]
fn py_integrate(py:Python<'_>,f: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    integrate_result(py,&to_expr(f)?,args,_kwargs)
}

/// Parse the common scalar and matrix expansion options once.
pub(crate) fn parse_expand_options(deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<(expand::ExpandOptions,expand::ExpandBudget)> {
    let mut options=expand::ExpandOptions{deep,..Default::default()};
    let mut max_terms=100_000usize;
    if let Some(hints)=hints {
        for (name,value) in hints.iter(){
            let name:String=name.extract()?;
            if name=="max_terms" {
                if value.is_instance_of::<pyo3::types::PyBool>(){return Err(PyTypeError::new_err("max_terms must be a positive integer"));}
                max_terms=value.extract::<usize>().map_err(|_|PyValueError::new_err("max_terms must be a positive integer"))?;
                if max_terms==0{return Err(PyValueError::new_err("max_terms must be a positive integer"));}
                continue;
            }
            let slot=match name.as_str(){
                "basic"=>&mut options.basic,"mul"=>&mut options.mul,"multinomial"=>&mut options.multinomial,
                "power_base"=>&mut options.power_base,"power_exp"=>&mut options.power_exp,"log"=>&mut options.log,
                "trig"=>&mut options.trig,"func"=>&mut options.func,"force"=>&mut options.force,
                _=>return Err(PyNotImplementedError::new_err(format!("expand hint '{name}' is not implemented"))),
            };
            *slot=value.is_truthy()?;
        }
    }
    Ok((options,expand::ExpandBudget::new(max_terms)))
}
/// Like `parse_expand_options`, also accepting SymPy's `complex` hint.
pub(crate) fn parse_scalar_expand_options(deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<(expand::ExpandOptions,expand::ExpandBudget,bool)>{
    let Some(hints)=hints else {let (options,budget)=parse_expand_options(deep,None)?;return Ok((options,budget,false));};
    let Some(complex)=hints.get_item("complex")? else {let (options,budget)=parse_expand_options(deep,Some(hints))?;return Ok((options,budget,false));};
    let rest=hints.copy()?;
    rest.del_item("complex")?;
    let (options,budget)=parse_expand_options(deep,Some(&rest))?;
    Ok((options,budget,complex.is_truthy()?))
}
pub(crate) fn expand_scalar(e:&Expr,options:&expand::ExpandOptions,budget:&mut expand::ExpandBudget,complex:bool)->PyResult<Expr>{
    if complex{
        return symox::complex_parts::expand_complex(e,options).map_err(|error|PyNotImplementedError::new_err(error.to_string()));
    }
    expand::expand_with_options(e,options,budget).map_err(|e|PyValueError::new_err(e.to_string()))
}
fn expand_hints(e:&Expr,deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<Expr>{
    let (options,mut budget,complex)=parse_scalar_expand_options(deep,hints)?;
    expand_scalar(e,&options,&mut budget,complex)
}

#[pyfunction]
#[pyo3(name = "expand", signature = (e, deep=true, **hints))]
fn py_expand(e:&Bound<'_,PyAny>,deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<PyObject>{
    if e.is_instance_of::<matrices::PyMatrix>() {return Ok(e.call_method("expand",(deep,),hints)?.unbind());}
    if let Ok(matrix)=e.downcast::<matrixexpr::PyMatrixExpr>(){
        let (options,mut budget)=parse_expand_options(deep,hints)?;
        let result=matrix.get().m.expand_with_options(&options,&mut budget).map_err(|error|matrices::err(e.py(),error))?;
        return expr_object(e.py(),matrix_proxy(result));
    }
    expr_object(e.py(),expand_hints(&to_expr(e)?,deep,hints)?)
}

#[pyfunction]
#[pyo3(name = "simplify", signature = (e, **_kwargs))]
fn py_simplify(e: &Bound<'_, PyAny>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if _kwargs.is_some_and(|options| !options.is_empty()) {
        return Err(PyNotImplementedError::new_err("custom simplify options are not implemented"));
    }
    if let Ok(matrix) = e.downcast::<matrixexpr::PyMatrixExpr>() {
        let result = matrix.get().m.scalar_map(&mut simplify::simplify).and_then(|matrix| matrix.doit(true))
            .map_err(|error| matrices::err(e.py(), error))?;
        return expr_object(e.py(), matrix_proxy(result));
    }
    expr_object(e.py(), simplify::simplify(&to_expr(e)?))
}

#[pyfunction]
#[pyo3(name = "together")]
fn py_together(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simplify::together(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "radsimp")]
fn py_radsimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simplify::radsimp(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "N", signature = (e, n=15, subs=None, **kwargs))]
fn py_n(e: &Bound<'_, PyAny>, n: i64, subs: Option<&Bound<'_, PyDict>>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(wrap(evalf_expr(&to_expr(e)?, n, subs, kwargs)?))
}

#[pyfunction]
#[pyo3(name = "sqrt")]
fn py_sqrt(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(ops::sqrt(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "cbrt")]
fn py_cbrt(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(ops::pow(&to_expr(e)?, &Expr::rational(1, 3))))
}

#[pyfunction]
#[pyo3(name = "root")]
fn py_root(e: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(ops::pow(&to_expr(e)?, &ops::inv(&to_expr(n)?))))
}

#[pyfunction]
#[pyo3(name = "log", signature = (x, base=None))]
fn py_log(x: &Bound<'_, PyAny>, base: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let x = to_expr(x)?;
    Ok(wrap(match base {
        Some(b) => fx::log_base(&x, &to_expr(b)?),
        None => fx::log(&x),
    }))
}

/// `Mod(p, q)` / `Rem(p, q)` raising SymPy's `ZeroDivisionError` for a zero
/// divisor (the native constructor leaves such nodes unevaluated).
fn checked_mod(name: &str, p: &Expr, q: &Expr) -> PyResult<Expr> {
    if symox::assumptions::is_zero(q) == Some(true) {
        let msg = if name == "Mod" { "Modulo by zero" } else { "Division by zero" };
        return Err(PyZeroDivisionError::new_err(msg));
    }
    let f = if name == "Mod" { symox::special::modulo } else { symox::special::rem };
    Ok(f(p, q))
}

/// `floor(p/q)`; a zero numeric divisor of a number raises like SymPy.
fn floor_div(p: &Expr, q: &Expr) -> PyResult<Expr> {
    if p.is_number() && q.is_zero() {
        return Err(PyZeroDivisionError::new_err("integer division or modulo by zero"));
    }
    Ok(fx::floor(&ops::div(p, q)))
}

/// Apply a known function by name (used by the Python layer).
#[pyfunction]
#[pyo3(name = "_apply", signature = (name, *args, evaluate=true))]
fn py_apply(name: &str, args: &Bound<'_, PyTuple>, evaluate: bool) -> PyResult<PyExpr> {
    if name=="Lambda"{
        if args.len()!=2{return Err(PyTypeError::new_err("Lambda requires a signature and an expression"));}
        return logic::py_lambda(&args.get_item(0)?,&args.get_item(1)?);
    }
    if name=="IdentityFunction" && !args.is_empty(){return Err(PyTypeError::new_err("IdentityFunction takes no arguments"));}
    if name=="Integral" || name=="Sum" || name=="Product" {
        if args.is_empty(){return Err(PyValueError::new_err("calculus constructor requires an integrand"));}
        let limits=PyTuple::new(args.py(),args.iter().skip(1))?;
        return calculus::_calculus(name,&args.get_item(0)?,&limits);
    }
    let v: Vec<Expr> = args.iter().map(|a| to_expr(&a)).collect::<PyResult<_>>()?;
    if evaluate && (name == "Mod" || name == "Rem") && v.len() == 2 {
        return checked_mod(name, &v[0], &v[1]).map(wrap);
    }
    match Func::from_name(name) {
        Some(f) => {
            if let Func::Special(s) = &f {
                if !s.nargs().contains(&v.len()) {
                    return Err(PyTypeError::new_err(format!(
                        "{name} accepts argument counts {:?}, got {}", s.nargs(), v.len()
                    )));
                }
                if let Err((type_error, msg)) = symox::special::numtheory::check_args(*s, &v) {
                    return Err(if type_error { PyTypeError::new_err(msg) } else { PyValueError::new_err(msg) });
                }
            }
            if let Some(n) = f.arity() {
                if v.len() != n {
                    return Err(PyTypeError::new_err(format!("{name} takes exactly {n} argument(s)")));
                }
            }
            Ok(wrap(if evaluate { fx::func(f, v) } else { Expr::from_unevaluated_kind(Kind::Func(f, v)) }))
        }
        None => Ok(wrap(fx::undefined(name, v))),
    }
}

#[pyfunction]
#[pyo3(name = "Rel", signature = (op, a, b, evaluate=true))]
fn py_rel(op: &str, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, evaluate: bool) -> PyResult<PyExpr> {
    let op = match op {
        "==" => RelOp::Eq,
        "!=" => RelOp::Ne,
        "<" => RelOp::Lt,
        "<=" => RelOp::Le,
        ">" => RelOp::Gt,
        ">=" => RelOp::Ge,
        _ => return Err(PyValueError::new_err(format!("unknown relational {op}"))),
    };
    let (a,b)=(to_expr(a)?,to_expr(b)?);
    Ok(wrap(if evaluate { relational::rel(op,&a,&b) } else { Expr::from_unevaluated_kind(Kind::Rel(op,a,b)) }))
}

#[pyfunction]
#[pyo3(name = "Piecewise", signature = (*pieces, evaluate=true))]
fn py_piecewise(pieces: &Bound<'_, PyTuple>, evaluate: bool) -> PyResult<PyExpr> {
    let mut v = Vec::new();
    for p in pieces.iter() {
        let t = integration_limit(&p)?.ok_or_else(||PyTypeError::new_err("Piecewise expects (expr, cond) tuples"))?;
        if t.len() != 2 {
            return Err(PyTypeError::new_err("Piecewise expects pairs of expression and condition"));
        }
        v.push((t[0].clone(),t[1].clone()));
    }
    Ok(wrap(if evaluate { relational::piecewise(v) } else { Expr::from_unevaluated_kind(Kind::Piecewise(v)) }))
}

#[pyfunction]
#[pyo3(name = "free_symbols")]
fn py_free_symbols<'py>(py: Python<'py>, e: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PySet>> {
    PySet::new(py, traverse::free_symbols(&to_expr(e)?).into_iter().map(wrap))
}

#[pyfunction]
#[pyo3(name = "constant")]
fn py_constant(name: &str) -> PyResult<PyExpr> {
    Ok(wrap(match name {
        "pi" => Expr::pi(),
        "E" => Expr::e(),
        "I" => Expr::i(),
        "oo" => Expr::oo(),
        "zoo" => Expr::zoo(),
        "nan" => Expr::nan(),
        "EulerGamma" => Expr::constant(Constant::EulerGamma),
        "Catalan" => Expr::constant(Constant::Catalan),
        "GoldenRatio" => Expr::constant(Constant::GoldenRatio),
        "TribonacciConstant" => Expr::constant(Constant::TribonacciConstant),
        "true" => Expr::boolean(true),
        "false" => Expr::boolean(false),
        _ => return Err(PyValueError::new_err(format!("unknown constant {name}"))),
    }))
}

// ----------------------------------------------------------------------
// lambdify
// ----------------------------------------------------------------------

/// A compiled numerical function.
#[pyclass(name = "CompiledFunction", module = "symox._core", frozen)]
pub struct PyCompiled {
    c: Compiled,
    scalar_out: bool,
    // Printed on demand: rendering large expressions costs more than compiling them.
    exprs: Vec<Expr>,
}

#[pymethods]
impl PyCompiled {
    #[pyo3(signature = (*args))]
    fn __call__(&self, py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
        if args.len() != self.c.n_args() {
            return Err(PyTypeError::new_err(format!(
                "expected {} argument(s), got {}",
                self.c.n_args(),
                args.len()
            )));
        }
        // scalar fast path
        let scalars: Option<Vec<f64>> = args
            .iter()
            .map(|a| if a.is_instance_of::<PyFloat>() || a.is_instance_of::<PyInt>() { a.extract::<f64>().ok() } else { None })
            .collect();
        if let Some(s) = scalars {
            let mut out = vec![0.0; self.c.n_outputs()];
            self.c.call(&s, &mut out);
            return if self.scalar_out {
                Ok(out[0].into_pyobject(py)?.into_any().unbind())
            } else {
                Ok(PyList::new(py, out)?.into_any().unbind())
            };
        }
        // vectorised path: broadcast scalars against 1-d arrays
        let np = py.import("numpy")?;
        let mut arrays: Vec<Vec<f64>> = Vec::with_capacity(args.len());
        let mut n: Option<usize> = None;
        let mut shape: Option<PyObject> = None;
        let mut scalar_flags = Vec::new();
        for a in args.iter() {
            let arr = np.call_method1("asarray", (a, "float64"))?;
            let ndim: usize = arr.getattr("ndim")?.extract()?;
            if ndim == 0 {
                scalar_flags.push(true);
                arrays.push(vec![arr.extract::<f64>()?]);
                continue;
            }
            scalar_flags.push(false);
            let size: usize = arr.getattr("size")?.extract()?;
            if let Some(m) = n {
                if m != size {
                    return Err(PyValueError::new_err("array arguments must have the same size"));
                }
            }
            n = Some(size);
            if shape.is_none() {
                shape = Some(arr.getattr("shape")?.unbind());
            }
            let flat = arr.call_method1("reshape", (-1,))?;
            let flat = np.call_method1("ascontiguousarray", (flat,))?;
            let ro: PyReadonlyArray1<f64> = flat.extract()?;
            arrays.push(ro.as_slice()?.to_vec());
        }
        let n = n.unwrap_or(1);
        for (i, s) in scalar_flags.iter().enumerate() {
            if *s {
                let v = arrays[i][0];
                arrays[i] = vec![v; n];
            }
        }
        let mut outs: Vec<Vec<f64>> = vec![vec![0.0; n]; self.c.n_outputs()];
        {
            let arg_slices: Vec<&[f64]> = arrays.iter().map(|v| v.as_slice()).collect();
            let mut out_slices: Vec<&mut [f64]> = outs.iter_mut().map(|v| v.as_mut_slice()).collect();
            py.allow_threads(|| self.c.call_many(&arg_slices, &mut out_slices));
        }
        let shape = shape.unwrap();
        let mut results = Vec::new();
        for o in outs {
            let a = PyArray1::from_vec(py, o);
            results.push(a.reshape_with_order(
                shape.extract::<Vec<usize>>(py)?,
                numpy::npyffi::NPY_ORDER::NPY_CORDER,
            )?.into_any());
        }
        if self.scalar_out {
            Ok(results.pop().unwrap().unbind())
        } else {
            Ok(PyList::new(py, results)?.into_any().unbind())
        }
    }

    fn __repr__(&self) -> String {
        let text = self.exprs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join(", ");
        format!("<compiled {}>", text)
    }

    #[getter]
    fn nargs(&self) -> usize {
        self.c.n_args()
    }
}

#[pyfunction]
#[pyo3(name = "lambdify", signature = (args, expr, modules=None, **_kwargs))]
fn py_lambdify(
    args: &Bound<'_, PyAny>,
    expr: &Bound<'_, PyAny>,
    modules: Option<&Bound<'_, PyAny>>,
    _kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyCompiled> {
    let _ = modules;
    let arg_list: Vec<Expr> = if let Ok(e) = to_expr(args) {
        vec![e]
    } else {
        args.try_iter()?.map(|a| to_expr(&a?)).collect::<PyResult<_>>()?
    };
    let (exprs, scalar_out): (Vec<Expr>, bool) = if let Ok(e) = to_expr(expr) {
        (vec![e], true)
    } else {
        (expr.try_iter()?.map(|a| to_expr(&a?)).collect::<PyResult<_>>()?, false)
    };
    let exprs: Vec<Expr> = exprs.iter().map(unwrap_unevaluated).collect();
    let c = Compiled::new(&arg_list, &exprs).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyCompiled { c, scalar_out, exprs })
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyExprNode>()?;
    m.add_class::<Head>()?;
    m.add_class::<PyCompiled>()?;
    m.add_function(wrap_pyfunction!(py_node, m)?)?;
    ask::register(m)?;
    m.add_function(wrap_pyfunction!(_restore_arithmetic, m)?)?;
    m.add_function(wrap_pyfunction!(_restore_logic, m)?)?;
    m.add_function(wrap_pyfunction!(py_default_sort_key, m)?)?;
    m.add_function(wrap_pyfunction!(py_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(py_dummy, m)?)?;
    m.add_function(wrap_pyfunction!(sympify, m)?)?;
    m.add_function(wrap_pyfunction!(py_wild, m)?)?;
    m.add_function(wrap_pyfunction!(py_integer, m)?)?;
    m.add_function(wrap_pyfunction!(py_rational, m)?)?;
    m.add_function(wrap_pyfunction!(py_float, m)?)?;
    m.add_function(wrap_pyfunction!(py_srepr, m)?)?;
    m.add_function(wrap_pyfunction!(_materialize_matrix_scalars, m)?)?;
    m.add_function(wrap_pyfunction!(py_diff, m)?)?;
    m.add_function(wrap_pyfunction!(py_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(_undefined_function, m)?)?;
    m.add_function(wrap_pyfunction!(_function_facts, m)?)?;
    m.add_function(wrap_pyfunction!(_as_real_imag_many, m)?)?;
    m.add_function(wrap_pyfunction!(_subs, m)?)?;
    m.add_function(wrap_pyfunction!(py_expand, m)?)?;
    m.add_function(wrap_pyfunction!(py_integrate, m)?)?;
    m.add_function(wrap_pyfunction!(py_summation, m)?)?;
    m.add_function(wrap_pyfunction!(py_bernoulli, m)?)?;
    m.add_function(wrap_pyfunction!(py_simplify, m)?)?;
    m.add_function(wrap_pyfunction!(py_together, m)?)?;
    m.add_function(wrap_pyfunction!(py_radsimp, m)?)?;
    m.add_function(wrap_pyfunction!(py_n, m)?)?;
    m.add_function(wrap_pyfunction!(py_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(py_cbrt, m)?)?;
    m.add_function(wrap_pyfunction!(py_root, m)?)?;
    m.add_function(wrap_pyfunction!(py_log, m)?)?;
    m.add_function(wrap_pyfunction!(py_apply, m)?)?;
    m.add_function(wrap_pyfunction!(py_rel, m)?)?;
    m.add_function(wrap_pyfunction!(py_piecewise, m)?)?;
    m.add_function(wrap_pyfunction!(py_free_symbols, m)?)?;
    m.add_function(wrap_pyfunction!(py_constant, m)?)?;
    m.add_function(wrap_pyfunction!(py_lambdify, m)?)?;
    matrices::register(m)?;
    matrixexpr::register(m)?;
    matrix_dag::register(m)?;
    linear_solvers::register(m)?;
    poly::register(m)?;
    ntheory::register(m)?;
    simp::register(m)?;
    simptools::register(m)?;
    solve::register(m)?;
    ode::register(m)?;
    logic::register(m)?;
    special::register(m)?;
    calculus::register(m)?;
    recurrences::register(m)?;
    transforms::register(m)?;
    series::register(m)?;
    printing::register(m)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
