//! Python bindings for boolean logic, sets and `solveset`.

use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use num_bigint::BigInt;

use symox::expr::{Expr, Kind};
use symox::logic::{self as lg, LogicForm};
use symox::sets::{self as st};

use crate::{to_expr, wrap, PyExpr};

pyo3::create_exception!(_core,BadSignatureError,PyTypeError);
pyo3::create_exception!(_core,BadArgumentsError,PyTypeError);
pub(crate) fn bad_signature(message:String)->PyErr{BadSignatureError::new_err(message)}
pub(crate) fn lambda_argument_error(message:String)->PyErr{
    if message.starts_with("Duplicate symbol") || message.starts_with("Lambda signature"){bad_signature(message)}
    else if message.starts_with("calculus ") || message.starts_with("Integral ") || message.starts_with("Integral/") || message.starts_with("Limit "){PyValueError::new_err(message)}
    else{BadArgumentsError::new_err(message)}
}
/// Python sequences in signatures/arguments are symbolic Tuple nodes.
pub(crate) fn lambda_value(value:&Bound<'_,PyAny>)->PyResult<Expr>{
    if value.is_instance_of::<PyTuple>() || value.is_instance_of::<PyList>(){
        return Ok(symox::special::tuple(value.try_iter()?.map(|v|lambda_value(&v?)).collect::<PyResult<_>>()?));
    }
    to_expr(value)
}
pub(crate) fn apply_lambda_expr(e:&Expr,args:&Bound<'_,PyTuple>)->PyResult<PyExpr>{
    let args=args.iter().map(|v|lambda_value(&v)).collect::<PyResult<Vec<_>>>()?;
    st::lambda_apply(e,&args).map(wrap).map_err(lambda_argument_error)
}

fn exprs(args: &Bound<'_, PyTuple>) -> PyResult<Vec<Expr>> {
    args.iter().map(|a| to_expr(&a)).collect()
}

fn boolean_arg(o: &Bound<'_, PyAny>) -> PyResult<Expr> {
    // Python ints 0/1 behave like False/True (SymPy's `as_Boolean`)
    if let Ok(i) = o.extract::<i64>() {
        if !o.is_instance_of::<PyExpr>() && (i == 0 || i == 1) {
            return Ok(Expr::boolean(i == 1));
        }
    }
    let e = to_expr(o)?;
    if e.is_integer_number() && (e.is_zero() || e.is_one()) {
        return Ok(Expr::boolean(e.is_one()));
    }
    if !lg::is_boolean(&e) {
        return Err(PyTypeError::new_err(format!("expecting bool or Boolean, not `{e}`.")));
    }
    Ok(e)
}

fn bool_args(args: &Bound<'_, PyTuple>) -> PyResult<Vec<Expr>> {
    args.iter().map(|a| boolean_arg(&a)).collect()
}

fn truth_arg(o: &Bound<'_, PyAny>) -> PyResult<Expr> {
    let e = to_expr(o)?;
    match e.kind() {
        Kind::Num(n) => Ok(Expr::boolean(!n.is_zero())),
        Kind::Const(symox::expr::Constant::Infinity | symox::expr::Constant::NegInfinity
            | symox::expr::Constant::NaN | symox::expr::Constant::ComplexInfinity) => Ok(Expr::boolean(true)),
        _ => Ok(e),
    }
}

fn truth_args(args: &Bound<'_, PyTuple>) -> PyResult<Vec<Expr>> {
    args.iter().map(|a| truth_arg(&a)).collect()
}

fn exact(name: &str, v: &[Expr], n: usize) -> PyResult<()> {
    if v.len() != n {
        return Err(PyTypeError::new_err(format!("{name} expects {n} arguments, got {}", v.len())));
    }
    Ok(())
}

macro_rules! variadic {
    ($fname:ident, $py:literal, $f:path, $arguments:path) => {
        #[pyfunction]
        #[pyo3(name = $py, signature = (*args))]
        fn $fname(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
            Ok(wrap($f($arguments(args)?)))
        }
    };
}

variadic!(py_and, "And", lg::and, bool_args);
variadic!(py_or, "Or", lg::or, bool_args);
variadic!(py_xor, "Xor", lg::xor, truth_args);
variadic!(py_nand, "Nand", lg::nand, bool_args);
variadic!(py_nor, "Nor", lg::nor, bool_args);
variadic!(py_xnor, "Xnor", lg::xnor, truth_args);
variadic!(py_equivalent, "Equivalent", lg::equivalent, truth_args);

#[pyfunction]
#[pyo3(name = "Not")]
fn py_not(a: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(lg::not(&truth_arg(a)?)))
}

#[pyfunction]
#[pyo3(name = "Implies", signature = (*args))]
fn py_implies(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let v = truth_args(args)?;
    exact("Implies", &v, 2)?;
    Ok(wrap(lg::implies(&v[0], &v[1])))
}

#[pyfunction]
#[pyo3(name = "ITE", signature = (*args))]
fn py_ite(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let v = bool_args(args)?;
    exact("ITE", &v, 3)?;
    Ok(wrap(lg::ite(&v[0], &v[1], &v[2])))
}

#[pyfunction]
#[pyo3(signature = (expr, simplify = false, force = false))]
fn to_cnf(expr: &Bound<'_, PyAny>, simplify: bool, force: bool) -> PyResult<PyExpr> {
    let e = to_expr(expr)?;
    if simplify && force { check_table_size(lg::atoms(&e).len())?; }
    lg::to_cnf(&e, simplify, force).map(wrap).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (expr, simplify = false, force = false))]
fn to_dnf(expr: &Bound<'_, PyAny>, simplify: bool, force: bool) -> PyResult<PyExpr> {
    let e = to_expr(expr)?;
    if simplify && force { check_table_size(lg::atoms(&e).len())?; }
    lg::to_dnf(&e, simplify, force).map(wrap).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (expr, simplify = true))]
fn to_nnf(expr: &Bound<'_, PyAny>, simplify: bool) -> PyResult<PyExpr> {
    Ok(wrap(lg::to_nnf(&to_expr(expr)?, simplify)))
}

#[pyfunction]
fn eliminate_implications(expr: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(lg::eliminate_implications(&to_expr(expr)?)))
}

#[pyfunction]
fn distribute_and_over_or(expr: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(lg::distribute_and_over_or(&to_expr(expr)?)))
}

#[pyfunction]
fn distribute_or_over_and(expr: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(lg::distribute_or_over_and(&to_expr(expr)?)))
}

#[pyfunction]
fn is_cnf(expr: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(lg::is_cnf(&to_expr(expr)?))
}

#[pyfunction]
fn is_dnf(expr: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(lg::is_dnf(&to_expr(expr)?))
}

#[pyfunction]
#[pyo3(signature = (expr, simplified = true))]
fn is_nnf(expr: &Bound<'_, PyAny>, simplified: bool) -> PyResult<bool> {
    Ok(lg::is_nnf(&to_expr(expr)?, simplified))
}

#[pyfunction]
fn is_literal(expr: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(lg::is_literal(&to_expr(expr)?))
}

fn parse_form(form: Option<&str>) -> PyResult<Option<LogicForm>> {
    match form {
        None => Ok(None),
        Some("cnf") => Ok(Some(LogicForm::Cnf)),
        Some("dnf") => Ok(Some(LogicForm::Dnf)),
        Some(_) => Err(PyValueError::new_err("form can be cnf or dnf only")),
    }
}

/// Row indices from SymPy style term specifications (bit lists, ints or
/// dicts mapping variables to 0/1).
fn terms_to_rows(vars: &[Expr], terms: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<u32>> {
    let n = vars.len();
    check_table_size(n)?;
    let mut rows = Vec::new();
    let Some(t) = terms else { return Ok(rows) };
    if t.is_none() {
        return Ok(rows);
    }
    for item in t.try_iter()? {
        let item = item?;
        if let Ok(i) = item.extract::<u32>() {
            if i >= (1u32 << n) {
                return Err(PyValueError::new_err("minterm index exceeds the number of variables"));
            }
            rows.push(i);
        } else if let Ok(d) = item.downcast::<PyDict>() {
            // partial assignment: expand the unspecified variables
            let mut fixed: Vec<Option<bool>> = vec![None; n];
            for (k, v) in d.iter() {
                let k = to_expr(&k)?;
                let pos = vars.iter().position(|x| *x == k).ok_or_else(|| PyValueError::new_err("unknown variable in term"))?;
                fixed[pos] = Some(v.is_truthy()?);
            }
            for r in 0..(1u32 << n) {
                if (0..n).all(|i| fixed[i].is_none_or(|b| ((r >> (n - 1 - i)) & 1 == 1) == b)) {
                    rows.push(r);
                }
            }
        } else {
            let bits: Vec<i64> = item.extract()?;
            if bits.len() != n {
                return Err(PyValueError::new_err("term length does not match the number of variables"));
            }
            let mut r = 0u32;
            for b in bits {
                r = r << 1 | (b != 0) as u32;
            }
            rows.push(r);
        }
    }
    Ok(rows)
}

fn expr_seq(o: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
    o.try_iter()?.map(|i| to_expr(&i?)).collect()
}

// Dense truth tables and the u32 implicant representation require a bound.
fn check_table_size(n: usize) -> PyResult<()> {
    if n > 20 {
        return Err(PyNotImplementedError::new_err("truth-table operations support at most 20 variables"));
    }
    Ok(())
}

fn validate_terms(vars: &[Expr], on: &[u32], dc: &[u32]) -> PyResult<()> {
    for (i, v) in vars.iter().enumerate() {
        if !lg::is_boolean(v) || vars[..i].contains(v) {
            return Err(PyValueError::new_err("variables must be distinct boolean expressions"));
        }
    }
    if on.iter().any(|r| dc.contains(r)) {
        return Err(PyValueError::new_err("minterms must not also occur in dontcares"));
    }
    Ok(())
}

#[pyfunction]
#[pyo3(name = "SOPform", signature = (variables, minterms, dontcares = None))]
fn sop_form(variables: &Bound<'_, PyAny>, minterms: &Bound<'_, PyAny>, dontcares: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let vars = expr_seq(variables)?;
    let on = terms_to_rows(&vars, Some(minterms))?;
    let dc = terms_to_rows(&vars, dontcares)?;
    validate_terms(&vars, &on, &dc)?;
    Ok(wrap(lg::sop_form(&vars, &on, &dc)))
}

#[pyfunction]
#[pyo3(name = "POSform", signature = (variables, minterms, dontcares = None))]
fn pos_form(variables: &Bound<'_, PyAny>, minterms: &Bound<'_, PyAny>, dontcares: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let vars = expr_seq(variables)?;
    let on = terms_to_rows(&vars, Some(minterms))?;
    let dc = terms_to_rows(&vars, dontcares)?;
    validate_terms(&vars, &on, &dc)?;
    Ok(wrap(lg::pos_form(&vars, &on, &dc)))
}

#[pyfunction]
#[pyo3(signature = (expr, form = None, deep = true, force = false, dontcare = None))]
fn simplify_logic(
    expr: &Bound<'_, PyAny>,
    form: Option<&str>,
    deep: bool,
    force: bool,
    dontcare: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    let e = to_expr(expr)?;
    let form = parse_form(form)?;
    let dc = dontcare.filter(|d| !d.is_none()).map(to_expr).transpose()?;
    if force {
        let mut vars = lg::atoms(&e);
        if let Some(d) = &dc { vars.extend(lg::atoms(d)); }
        lg::ordered(&mut vars);
        vars.dedup();
        check_table_size(vars.len())?;
    }
    Ok(wrap(lg::simplify_logic_dontcare(&e, form, deep, force, dc.as_ref())))
}

#[pyfunction]
#[pyo3(signature = (expr, deep = true))]
fn to_anf(expr: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    Ok(wrap(lg::to_anf(&to_expr(expr)?, deep)))
}

#[pyfunction]
fn is_anf(expr: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(lg::is_anf(&to_expr(expr)?))
}

fn truthvalues(o: &Bound<'_, PyAny>) -> PyResult<Vec<BigInt>> {
    // Accept bools and arbitrary exact integers, as SymPy's bitwise XOR
    // transform does; a float is not an integer truth-table entry.
    o.try_iter()?.map(|v| v?.extract::<BigInt>()).collect()
}

#[pyfunction]
fn anf_coeffs(truthvalues: &Bound<'_, PyAny>) -> PyResult<Vec<BigInt>> {
    let values = self::truthvalues(truthvalues)?;
    if !values.len().is_power_of_two() {
        return Err(PyValueError::new_err(format!("The number of truth values must be a power of two, got {}", values.len())));
    }
    check_table_size(values.len().trailing_zeros() as usize)?;
    Ok(lg::anf_coeffs(&values))
}

#[pyfunction]
#[pyo3(name = "ANFform")]
fn anf_form(variables: &Bound<'_, PyAny>, truthvalues: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let vars = expr_seq(variables)?;
    check_table_size(vars.len())?;
    let values = self::truthvalues(truthvalues)?;
    if values.len() != 1usize << vars.len() {
        return Err(PyValueError::new_err(format!("The number of truth values must be equal to 2^{}, got {}", vars.len(), values.len())));
    }
    for v in &vars {
        if !lg::is_boolean(v) { return Err(PyTypeError::new_err("ANF variables must be Boolean expressions")); }
    }
    let coefficients: Vec<bool> = lg::anf_coeffs(&values).into_iter().map(|c| c == BigInt::from(1)).collect();
    Ok(wrap(lg::anf_from_coeffs(&vars, &coefficients)))
}

fn model_dict<'py>(py: Python<'py>, m: &lg::Model) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    for (k, v) in m {
        d.set_item(wrap(k.clone()), *v)?;
    }
    Ok(d)
}

/// Iterator over the models of a formula (`satisfiable(..., all_models=True)`).
#[pyclass(module = "symox._core", name = "ModelIterator", unsendable)]
pub struct ModelIterator {
    it: lg::sat::ModelIter,
    started: bool,
    finished: bool,
    trivially_true: bool,
}

#[pymethods]
impl ModelIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PyObject>> {
        if self.finished {
            return Ok(None);
        }
        if self.trivially_true {
            self.finished = true;
            let d = PyDict::new(py);
            d.set_item(wrap(Expr::boolean(true)), true)?;
            return Ok(Some(d.into_any().unbind()));
        }
        match self.it.next() {
            Some(m) => {
                self.started = true;
                Ok(Some(model_dict(py, &m)?.into_any().unbind()))
            }
            None => {
                self.finished = true;
                if self.started {
                    Ok(None)
                } else {
                    Ok(Some(false.into_pyobject(py)?.to_owned().into_any().unbind()))
                }
            }
        }
    }
}

#[pyfunction]
#[pyo3(signature = (expr, algorithm = None, all_models = false, minimal = false, use_lra_theory = false))]
fn satisfiable(
    py: Python<'_>,
    expr: &Bound<'_, PyAny>,
    algorithm: Option<&str>,
    all_models: bool,
    minimal: bool,
    use_lra_theory: bool,
) -> PyResult<PyObject> {
    if use_lra_theory || minimal {
        return Err(PyNotImplementedError::new_err("LRA theory and minimal SAT models are not supported"));
    }
    match algorithm {
        None | Some("dpll") | Some("dpll2") => {}
        Some("z3" | "pycosat" | "minisat22") => return Err(PyNotImplementedError::new_err("external SAT solvers are not supported")),
        Some(_) => return Err(PyValueError::new_err("unknown SAT algorithm")),
    }
    let e = boolean_arg(expr)?;
    let trivially_true = matches!(e.kind(), Kind::Bool(true));
    if all_models {
        let it = ModelIterator { it: lg::satisfiable_all(&e), started: false, finished: false, trivially_true };
        return Ok(Py::new(py, it)?.into_any());
    }
    if trivially_true {
        let d = PyDict::new(py);
        d.set_item(wrap(Expr::boolean(true)), true)?;
        return Ok(d.into_any().unbind());
    }
    let r = py.allow_threads(|| lg::satisfiable(&e));
    match r {
        Some(m) => Ok(model_dict(py, &m)?.into_any().unbind()),
        None => Ok(false.into_pyobject(py)?.to_owned().into_any().unbind()),
    }
}

/// Solve a DIMACS style CNF directly (list of clauses of non-zero ints).
#[pyfunction]
fn _sat_dimacs(py: Python<'_>, nvars: usize, clauses: Vec<Vec<i32>>) -> PyResult<PyObject> {
    if nvars > (i32::MAX as usize) || clauses.iter().flatten().any(|&v| v == 0 || v.unsigned_abs() as usize > nvars) {
        return Err(PyValueError::new_err("DIMACS literals must be nonzero and refer to an existing variable"));
    }
    match py.allow_threads(|| lg::sat::solve_dimacs(nvars, &clauses)) {
        Some(m) => Ok(PyList::new(py, m)?.into_any().unbind()),
        None => Ok(false.into_pyobject(py)?.to_owned().into_any().unbind()),
    }
}

#[pyfunction]
#[pyo3(signature = (expr, variables, input = true))]
fn truth_table(py: Python<'_>, expr: &Bound<'_, PyAny>, variables: &Bound<'_, PyAny>, input: bool) -> PyResult<PyObject> {
    let e = to_expr(expr)?;
    let vars = expr_seq(variables)?;
    check_table_size(vars.len())?;
    let rows = lg::truth_table(&e, &vars);
    let out = PyList::empty(py);
    for (bits, v) in rows {
        let val: PyObject = match v.kind() {
            Kind::Bool(_) => wrap(v).into_pyobject(py)?.into_any().unbind(),
            _ => wrap(v).into_pyobject(py)?.into_any().unbind(),
        };
        if input {
            let b: Vec<i64> = bits.iter().map(|&b| b as i64).collect();
            out.append(PyTuple::new(py, [PyList::new(py, b)?.into_any().unbind(), val])?)?;
        } else {
            out.append(val)?;
        }
    }
    Ok(out.into_any().unbind())
}

#[pyfunction]
fn atoms_logic(py: Python<'_>, expr: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    Ok(PyList::new(py, lg::atoms(&to_expr(expr)?).into_iter().map(wrap))?.into_any().unbind())
}

// ----------------------------------------------------------------------
// sets
// ----------------------------------------------------------------------

#[pyfunction]
#[pyo3(name = "Interval", signature = (start, end, left_open = false, right_open = false))]
fn py_interval(start: &Bound<'_, PyAny>, end: &Bound<'_, PyAny>, left_open: bool, right_open: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::interval(&to_expr(start)?, &to_expr(end)?, left_open, right_open)))
}

#[pyfunction]
#[pyo3(name = "FiniteSet", signature = (*args))]
fn py_finite_set(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    Ok(wrap(st::finite_set(exprs(args)?)))
}

fn require_set(e: &Expr) -> PyResult<()> {
    if !st::is_set(e) {
        return Err(PyTypeError::new_err(format!("expected a Set, got {e}")));
    }
    Ok(())
}

fn set_args(args: &Bound<'_, PyTuple>) -> PyResult<Vec<Expr>> {
    // accept a single iterable of sets too
    if args.len() == 1 {
        let a = args.get_item(0)?;
        if to_expr(&a).is_err() {
            let v = expr_seq(&a)?;
            for s in &v { require_set(s)?; }
            return Ok(v);
        }
    }
    let v = exprs(args)?;
    for s in &v {
        if !st::is_set(s) {
            return Err(PyTypeError::new_err(format!("Input args to Union must be Sets, got {s}")));
        }
    }
    Ok(v)
}

#[pyfunction]
#[pyo3(name = "Union", signature = (*args))]
fn py_union(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    Ok(wrap(st::union(set_args(args)?)))
}

#[pyfunction]
#[pyo3(name = "Intersection", signature = (*args))]
fn py_intersection(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    Ok(wrap(st::intersection(set_args(args)?)))
}

#[pyfunction]
#[pyo3(name = "Complement")]
fn py_complement(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let (a, b) = (to_expr(a)?, to_expr(b)?);
    require_set(&a)?;
    require_set(&b)?;
    Ok(wrap(st::complement(&a, &b)))
}

#[pyfunction]
#[pyo3(name = "ProductSet", signature = (*args))]
fn py_product_set(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    Ok(wrap(st::product_set(set_args(args)?)))
}

#[pyfunction]
#[pyo3(name = "Lambda")]
pub(crate) fn py_lambda(signature:&Bound<'_,PyAny>,expr:&Bound<'_,PyAny>)->PyResult<PyExpr>{
    st::lambda_checked(&lambda_value(signature)?,&lambda_value(expr)?).map(wrap).map_err(bad_signature)
}

#[pyfunction]
#[pyo3(signature=(lam,*args))]
fn _lambda_call(lam:&Bound<'_,PyAny>,args:&Bound<'_,PyTuple>)->PyResult<PyExpr>{
    apply_lambda_expr(&to_expr(lam)?,args)
}

#[pyfunction]
#[pyo3(name = "ImageSet")]
fn py_image_set(lam: &Bound<'_, PyAny>, base: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::image_set(&to_expr(lam)?, &to_expr(base)?)))
}

#[pyfunction]
#[pyo3(name = "ConditionSet", signature = (sym, condition, base_set = None))]
fn py_condition_set(sym: &Bound<'_, PyAny>, condition: &Bound<'_, PyAny>, base_set: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let base = match base_set {
        Some(b) => to_expr(b)?,
        None => st::universal_set(),
    };
    Ok(wrap(st::condition_set(&to_expr(sym)?, &to_expr(condition)?, &base)))
}

#[pyfunction]
#[pyo3(name = "Range", signature = (*args))]
fn py_range(args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let v = exprs(args)?;
    for (i, e) in v.iter().enumerate() {
        if e.as_number().is_some_and(|n| !n.is_integer()) || (i == 2 && e.is_zero()) {
            return Err(PyValueError::new_err("Range requires integer arguments and a nonzero step"));
        }
    }
    Ok(wrap(match v.len() {
        1 => st::range(&Expr::zero(), &v[0], &Expr::one()),
        2 => st::range(&v[0], &v[1], &Expr::one()),
        3 => st::range(&v[0], &v[1], &v[2]),
        _ => return Err(PyTypeError::new_err("Range expects 1 to 3 arguments")),
    }))
}

#[pyfunction]
#[pyo3(name = "Contains")]
fn py_contains(x: &Bound<'_, PyAny>, s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::contains_expr(&to_expr(x)?, &to_expr(s)?)))
}

#[pyfunction]
fn _named_set(name: &str) -> PyResult<PyExpr> {
    st::named_set(name).map(wrap).ok_or_else(|| PyValueError::new_err(format!("unknown set {name}")))
}

#[pyfunction]
fn _is_set(e: &Bound<'_, PyAny>) -> bool {
    to_expr(e).is_ok_and(|e| st::is_set(&e))
}

#[pyfunction]
fn _set_contains(s: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::contains(&to_expr(s)?, &to_expr(x)?)))
}

#[pyfunction]
fn _set_as_relational(s: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::as_relational(&to_expr(s)?, &to_expr(x)?)))
}

#[pyfunction]
fn _set_is_subset(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    Ok(st::is_subset(&to_expr(a)?, &to_expr(b)?))
}

fn prop(name: &str, s: &Bound<'_, PyAny>, f: fn(&Expr) -> Option<Expr>) -> PyResult<PyExpr> {
    let e = to_expr(s)?;
    f(&e).map(wrap).ok_or_else(|| PyNotImplementedError::new_err(format!("({e})._{name}")))
}

#[pyfunction]
fn _set_inf(s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    prop("inf", s, st::inf)
}
#[pyfunction]
fn _set_sup(s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    prop("sup", s, st::sup)
}
#[pyfunction]
fn _set_measure(s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    prop("measure", s, st::measure)
}
#[pyfunction]
fn _set_boundary(s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    prop("boundary", s, st::boundary)
}
#[pyfunction]
fn _set_interior(s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    prop("interior", s, st::interior)
}
#[pyfunction]
fn _set_closure(s: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    prop("closure", s, st::closure)
}

#[pyfunction]
#[pyo3(signature = (f, symbol = None, domain = None))]
fn solveset(f: &Bound<'_, PyAny>, symbol: Option<&Bound<'_, PyAny>>, domain: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let f = to_expr(f)?;
    let x = match symbol {
        Some(s) => to_expr(s)?,
        None => {
            let fs = symox::traverse::free_symbols(&f);
            match fs.len() {
                1 => fs[0].clone(),
                0 => Expr::symbol("x"),
                _ => {
                    return Err(PyValueError::new_err(
                        "The independent variable must be specified for a multivariate equation.",
                    ))
                }
            }
        }
    };
    let dom = match domain {
        Some(d) => to_expr(d)?,
        None => st::complexes(),
    };
    require_set(&dom)?;
    if !x.is_symbol() { return Err(PyValueError::new_err("the independent variable must be a symbol")); }
    st::solveset(&f, &x, &dom).map(wrap).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (expr, gen, relational = true, domain = None, continuous = false))]
fn solve_univariate_inequality(
    expr: &Bound<'_, PyAny>,
    gen: &Bound<'_, PyAny>,
    relational: bool,
    domain: Option<&Bound<'_, PyAny>>,
    continuous: bool,
) -> PyResult<PyExpr> {
    let _ = continuous;
    let dom = match domain {
        Some(d) => to_expr(d)?,
        None => st::reals(),
    };
    st::solve_univariate_inequality(&to_expr(expr)?, &to_expr(gen)?, relational, &dom)
        .map(wrap)
        .map_err(PyNotImplementedError::new_err)
}

#[pyfunction]
#[pyo3(signature = (inequalities, symbols = None))]
fn reduce_inequalities(inequalities: &Bound<'_, PyAny>, symbols: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let ineqs = match to_expr(inequalities) {
        Ok(e) => vec![e],
        Err(_) => expr_seq(inequalities)?,
    };
    let syms: Vec<Expr> = match symbols {
        None => vec![],
        Some(s) => match to_expr(s) {
            Ok(e) => vec![e],
            Err(_) => expr_seq(s)?,
        },
    };
    if syms.len() > 1 {
        return Err(PyNotImplementedError::new_err("only univariate inequalities are supported"));
    }
    st::reduce_inequalities(&ineqs, syms.first())
        .map(wrap)
        .map_err(PyNotImplementedError::new_err)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("BadSignatureError",m.py().get_type::<BadSignatureError>())?;
    m.add("BadArgumentsError",m.py().get_type::<BadArgumentsError>())?;
    m.add_class::<ModelIterator>()?;
    for f in [
        wrap_pyfunction!(py_and, m)?,
        wrap_pyfunction!(py_or, m)?,
        wrap_pyfunction!(py_not, m)?,
        wrap_pyfunction!(py_xor, m)?,
        wrap_pyfunction!(py_nand, m)?,
        wrap_pyfunction!(py_nor, m)?,
        wrap_pyfunction!(py_xnor, m)?,
        wrap_pyfunction!(py_implies, m)?,
        wrap_pyfunction!(py_equivalent, m)?,
        wrap_pyfunction!(py_ite, m)?,
        wrap_pyfunction!(to_cnf, m)?,
        wrap_pyfunction!(to_dnf, m)?,
        wrap_pyfunction!(to_nnf, m)?,
        wrap_pyfunction!(to_anf, m)?,
        wrap_pyfunction!(is_anf, m)?,
        wrap_pyfunction!(anf_form, m)?,
        wrap_pyfunction!(anf_coeffs, m)?,
        wrap_pyfunction!(eliminate_implications, m)?,
        wrap_pyfunction!(distribute_and_over_or, m)?,
        wrap_pyfunction!(distribute_or_over_and, m)?,
        wrap_pyfunction!(is_cnf, m)?,
        wrap_pyfunction!(is_dnf, m)?,
        wrap_pyfunction!(is_nnf, m)?,
        wrap_pyfunction!(is_literal, m)?,
        wrap_pyfunction!(sop_form, m)?,
        wrap_pyfunction!(pos_form, m)?,
        wrap_pyfunction!(simplify_logic, m)?,
        wrap_pyfunction!(satisfiable, m)?,
        wrap_pyfunction!(_sat_dimacs, m)?,
        wrap_pyfunction!(truth_table, m)?,
        wrap_pyfunction!(atoms_logic, m)?,
        wrap_pyfunction!(py_interval, m)?,
        wrap_pyfunction!(py_finite_set, m)?,
        wrap_pyfunction!(py_union, m)?,
        wrap_pyfunction!(py_intersection, m)?,
        wrap_pyfunction!(py_complement, m)?,
        wrap_pyfunction!(py_product_set, m)?,
        wrap_pyfunction!(py_lambda, m)?,
        wrap_pyfunction!(_lambda_call, m)?,
        wrap_pyfunction!(py_image_set, m)?,
        wrap_pyfunction!(py_condition_set, m)?,
        wrap_pyfunction!(py_range, m)?,
        wrap_pyfunction!(py_contains, m)?,
        wrap_pyfunction!(_named_set, m)?,
        wrap_pyfunction!(_is_set, m)?,
        wrap_pyfunction!(_set_contains, m)?,
        wrap_pyfunction!(_set_as_relational, m)?,
        wrap_pyfunction!(_set_is_subset, m)?,
        wrap_pyfunction!(_set_inf, m)?,
        wrap_pyfunction!(_set_sup, m)?,
        wrap_pyfunction!(_set_measure, m)?,
        wrap_pyfunction!(_set_boundary, m)?,
        wrap_pyfunction!(_set_interior, m)?,
        wrap_pyfunction!(_set_closure, m)?,
        wrap_pyfunction!(solveset, m)?,
        wrap_pyfunction!(solve_univariate_inequality, m)?,
        wrap_pyfunction!(reduce_inequalities, m)?,
    ] {
        m.add_function(f)?;
    }
    Ok(())
}
