//! Python bindings for `solve` and `nsolve`.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyFrozenSet, PyList, PySet, PyTuple};

use symox::solve::{self as sv, SystemSolution};
use symox::Expr;

use crate::{to_expr, wrap};

fn set_result(py: Python<'_>, variables: &[Expr], solutions: Vec<Vec<Expr>>) -> PyResult<PyObject> {
    let symbols = PyList::new(py, variables.iter().cloned().map(wrap))?;
    let rows = solutions.into_iter().map(|row| PyTuple::new(py, row.into_iter().map(wrap))).collect::<PyResult<Vec<_>>>()?;
    Ok((symbols, PySet::new(py, rows)?).into_pyobject(py)?.into_any().unbind())
}

fn solve_flags(flags: Option<&Bound<'_, PyDict>>) -> PyResult<(bool, bool, bool,symox::poly::roots::RadicalOptions)> {
    let (mut dict, mut set, mut check) = (false, false, true);
    let mut cubics=true;let mut quartics=None;
    if let Some(flags) = flags {
        for (key, value) in flags.iter() {
            let key: String = key.extract()?;
            match key.as_str() {
                "dict" => dict = value.is_truthy()?,
                "set" => set = value.is_truthy()?,
                "check" => check = value.is_truthy()?,
                "cubics"=>cubics=value.is_truthy()?,
                "quartics"=>quartics=Some(value.is_truthy()?),
                "quintics"=>{},
                "simplify" if value.is_truthy()? => {},
                "manual" | "warn" | "minimal" if !value.is_truthy()? => {},
                _ => return Err(PyNotImplementedError::new_err(format!("solve option '{key}' with this value is not implemented"))),
            }
        }
    }
    Ok((dict, set && !dict, check,symox::poly::roots::RadicalOptions{cubics,quartics:quartics.unwrap_or(cubics)}))
}

fn expr_list(o: &Bound<'_, PyAny>) -> PyResult<(Vec<Expr>, bool)> {
    if let Ok(e) = to_expr(o) {
        if let symox::Kind::Func(symox::Func::Undefined(name),arguments)=e.kind(){
            if name.is_internal("Tuple"){return Ok((arguments.clone(),true));}
        }
        return Ok((vec![e], false));
    }
    let mut v = Vec::new();
    for item in o.try_iter()? {
        v.push(to_expr(&item?)?);
    }
    Ok((v, true))
}

/// The scalar solve interface tries a unique undetermined-coefficient result
/// before ordinary algebraic solving. Restrict the automatic attempt to
/// rational functions of independent Symbols: their monomial basis is exact,
/// unlike potentially dependent functional bases. The existing native Python
/// coefficient collector handles sparse/high-degree data and calls the native
/// SYSTEM solver, so this does not recurse into the scalar heuristic.
fn undetermined_solution(py:Python<'_>,equation:&Expr,variables:&[Expr],check:bool,
    options:symox::poly::roots::RadicalOptions)->PyResult<Option<SystemSolution>>{
    if variables.len()<2||variables.iter().any(|v|!v.is_symbol()){return Ok(None);}
    let expression=sv::as_expr(equation);
    let independent:Vec<_>=sv::unknowns(&[expression.clone()]).into_iter().filter(|v|!variables.contains(v)).collect();
    if independent.is_empty()||symox::inspection::classify(&expression,&independent,
        symox::inspection::PolynomialKind::Rational)!=Some(true){return Ok(None);}
    let flags=PyDict::new(py);
    flags.set_item("dict",true)?;flags.set_item("check",check)?;
    flags.set_item("cubics",options.cubics)?;flags.set_item("quartics",options.quartics)?;
    // Coefficient matching chooses canonical dependent coefficients, while
    // the outer interface still preserves the user's requested tuple columns.
    let mut coefficients=variables.to_vec();
    coefficients.sort_by_cached_key(symox::sortkey::sort_key);
    let collector=py.import("symox._solver_checks")?.getattr("solve_undetermined_coeffs")?;
    let result=match collector.call((wrap(expression),PyList::new(py,coefficients.into_iter().map(wrap))?),Some(&flags)){
        Ok(value)=>value,
        Err(error) if error.is_instance_of::<PyNotImplementedError>(py)=>return Ok(None),
        Err(error)=>return Err(error),
    };
    let Ok(solutions)=result.downcast::<PyList>() else{return Ok(None);};
    if solutions.len()!=1{return Ok(None);}
    let solution=solutions.get_item(0)?;
    let Ok(mapping)=solution.downcast::<PyDict>() else{return Ok(None);};
    if mapping.is_empty(){return Ok(None);}
    let assignments=mapping.iter().map(|(key,value)|Ok((to_expr(&key)?,to_expr(&value)?))).collect::<PyResult<Vec<_>>>()?;
    if assignments.iter().any(|(key,_)|!variables.contains(key)){return Ok(None);}
    Ok(Some(SystemSolution::Map(assignments)))
}

#[pyfunction]
fn _nonlinsolve_parametric(f: &Bound<'_, PyAny>, symbols: &Bound<'_, PyAny>) -> PyResult<Vec<(Vec<crate::PyExpr>, Vec<Vec<crate::PyExpr>>)>> {
    let (equations, _) = expr_list(f)?;
    let (variables, _) = expr_list(symbols)?;
    let result = symox::poly::parametric::nonlinsolve(&equations, &variables).map_err(|e| match e {
        symox::poly::convert::PolyError::Unsupported(message) => PyNotImplementedError::new_err(message),
        other => PyValueError::new_err(other.to_string()),
    })?;
    Ok(result.into_iter().map(|b| (b.values.into_iter().map(wrap).collect(), b.excluded.into_iter().map(|r| r.into_iter().map(wrap).collect()).collect())).collect())
}

#[pyfunction]
#[pyo3(signature = (f, *symbols, **flags))]
fn solve(py: Python<'_>, f: &Bound<'_, PyAny>, symbols: &Bound<'_, PyTuple>, flags: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (as_dict, as_set, check,root_options) = solve_flags(flags)?;
    let (eqs, is_system) = expr_list(f)?;
    let unordered_symbols=symbols.len()==1&&symbols.get_item(0).is_ok_and(|s|
        s.is_instance_of::<PySet>()||s.is_instance_of::<PyFrozenSet>());
    let mut vars: Vec<Expr> = Vec::new();
    for s in symbols.iter() {
        let (v, _) = expr_list(&s)?;
        vars.extend(v);
    }
    let inferred_symbols=vars.is_empty();
    if inferred_symbols {
        vars = sv::unknowns(&eqs);
    }else if unordered_symbols{
        vars.sort_by_cached_key(symox::sortkey::sort_key);
    }
    if vars.is_empty() {
        if as_set { return set_result(py, &vars, vec![]); }
        return Ok(PyList::empty(py).into_any().unbind());
    }
    if !is_system && eqs.len() == 1 && vars.len() == 1 {
        let x = &vars[0];
        let sols = symox::poly::roots::with_radical_options(root_options,||
            sv::solve_atom(&eqs[0],x,check).or_else(||sv::solve_with_check(&eqs[0],x,check))).ok_or_else(|| {
            PyNotImplementedError::new_err(format!("could not solve {} for {}", eqs[0], x))
        })?;
        if as_set { return set_result(py, &vars, sols.into_iter().map(|s| vec![s]).collect()); }
        if as_dict {
            let l = PyList::empty(py);
            for s in sols {
                let d = PyDict::new(py);
                d.set_item(wrap(x.clone()), wrap(s))?;
                l.append(d)?;
            }
            return Ok(l.into_any().unbind());
        }
        return Ok(PyList::new(py, sols.into_iter().map(wrap))?.into_any().unbind());
    }
    // Python coefficient collection is outside the lexical native radical
    // scope; nested native system calls establish and restore their own scope.
    let coefficients=if !is_system&&eqs.len()==1{
        undetermined_solution(py,&eqs[0],&vars,check,root_options)?
    }else{None};
    let res=if let Some(solution)=coefficients{solution}else{
        symox::poly::roots::with_radical_options(root_options,||{
            let recast=if vars.len()==1{sv::solve_atom_system(&eqs,&vars[0],check)}else{None};
            recast.or_else(||sv::solve_system_with_check(&eqs,&vars,check))
        })
            .ok_or_else(||PyNotImplementedError::new_err("could not solve the system"))?
    };
    match res {
        SystemSolution::NoSolution if as_set => set_result(py, &vars, vec![]),
        SystemSolution::NoSolution => Ok(PyList::empty(py).into_any().unbind()),
        SystemSolution::Map(m) => {
            if as_set {
                if inferred_symbols||unordered_symbols{
                    return set_result(py,&m.iter().map(|(key,_)|key.clone()).collect::<Vec<_>>(),
                        vec![m.iter().map(|(_,value)|value.clone()).collect()]);
                }
                let row = vars.iter().map(|v| m.iter().find(|(k, _)| k == v).map(|(_, value)| value.clone()).unwrap_or_else(|| v.clone())).collect();
                return set_result(py, &vars, vec![row]);
            }
            let complete=m.len()==vars.len();
            if !is_system&&!as_dict&&!complete&&!inferred_symbols&&!unordered_symbols{
                let row=vars.iter().map(|v|m.iter().find(|(k,_)|k==v)
                    .map(|(_,value)|value.clone()).unwrap_or_else(||v.clone()));
                return Ok(PyList::new(py,[PyTuple::new(py,row.map(wrap))?])?.into_any().unbind());
            }
            let d = PyDict::new(py);
            for (k, v) in m {
                d.set_item(wrap(k), wrap(v))?;
            }
            if as_dict||!is_system&&!complete&&(inferred_symbols||unordered_symbols) {
                return Ok(PyList::new(py, [d])?.into_any().unbind());
            }
            Ok(d.into_any().unbind())
        }
        SystemSolution::Tuples(mut ts) => {
            if as_set {
                if inferred_symbols||unordered_symbols{
                    let selected:Vec<usize>=(0..vars.len()).filter(|&i|ts.iter().any(|row|row[i]!=vars[i])).collect();
                    let labels:Vec<_>=selected.iter().map(|&i|vars[i].clone()).collect();
                    let rows=ts.into_iter().map(|row|selected.iter().map(|&i|row[i].clone()).collect()).collect();
                    return set_result(py,&labels,rows);
                }
                return set_result(py, &vars, ts);
            }
            if !is_system&&!as_dict&&ts.len()==1&&vars.iter().zip(&ts[0]).all(|(key,value)|key!=value){
                let d=PyDict::new(py);
                for (key,value) in vars.iter().zip(&ts[0]){d.set_item(wrap(key.clone()),wrap(value.clone()))?;}
                return Ok(d.into_any().unbind());
            }
            if as_dict||inferred_symbols||unordered_symbols{
                ts.sort_by_cached_key(|row|vars.iter().zip(row).filter(|(key,value)|*key!=*value)
                    .map(|(key,value)|(symox::sortkey::sort_key(key),symox::sortkey::sort_key(value))).collect::<Vec<_>>());
            }
            let l = PyList::empty(py);
            for t in ts {
                if as_dict||inferred_symbols||unordered_symbols {
                    let d = PyDict::new(py);
                    for (k, v) in vars.iter().zip(t) {
                        if *k==v{continue;}
                        d.set_item(wrap(k.clone()), wrap(v))?;
                    }
                    l.append(d)?;
                } else {
                    l.append(PyTuple::new(py, t.into_iter().map(wrap))?)?;
                }
            }
            Ok(l.into_any().unbind())
        }
    }
}

#[pyfunction]
#[pyo3(signature = (f, x, x0=None, **_kwargs))]
fn nsolve(py: Python<'_>, f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, x0: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (mut digits, mut maxsteps, mut verify, mut as_dict) = (15_i64, 100_usize, true, false);
    let mut tolerance = None;
    let mut solver = None;
    if let Some(kwargs) = _kwargs {
        for (key, value) in kwargs.iter() {
            let key: String = key.extract()?;
            match key.as_str() {
                "tol" if value.is_none() => {},
                "tol" => tolerance = Some(to_expr(&value)?),
                "maxsteps" => maxsteps = value.extract()?,
                "verify" => verify = value.extract()?,
                "dict" => as_dict = value.is_truthy()?,
                "solver" => solver = Some(value.extract::<String>()?),
                "prec" => digits = value.extract()?,
                "verbose" if !value.is_truthy()? => {},
                _ => return Err(PyNotImplementedError::new_err(format!("nsolve option '{key}' with this value is not implemented"))),
            }
        }
    }
    let precision = crate::decimal_precision(digits)?;
    if maxsteps == 0 { return Err(PyValueError::new_err("maxsteps must be positive")); }
    let (eqs, sys) = expr_list(f)?;
    let (vars, initial) = if let Some(start) = x0 {
        (expr_list(x)?.0, expr_list(start)?.0)
    } else {
        let vars = sv::unknowns(&eqs);
        if vars.len() != 1 {
            return Err(PyValueError::new_err("nsolve needs an explicit variable list unless there is exactly one unknown"));
        }
        (vars, expr_list(x)?.0)
    };
    let solver = solver.unwrap_or_else(|| if sys { "newton" } else { "secant" }.to_owned());
    if !matches!(solver.as_str(), "newton" | "secant" | "bisect") {
        return Err(PyNotImplementedError::new_err(format!("nsolve solver '{solver}' is not implemented")));
    }
    if vars.is_empty() || vars.iter().any(|v| !v.is_symbol()) || vars.iter().enumerate().any(|(i, v)| vars[..i].contains(v)) {
        return Err(PyValueError::new_err("nsolve unknowns must be distinct symbols"));
    }
    if solver == "bisect" {
        if vars.len() != 1 || eqs.len() != 1 || initial.len() != 2 {
            return Err(PyValueError::new_err("bisection requires one equation, one symbol and a two-endpoint bracket"));
        }
    } else if solver == "secant" {
        if vars.len() != 1 || eqs.len() != 1 || initial.is_empty() || initial.len() > 2 {
            return Err(PyValueError::new_err("secant requires one equation, one unknown and one or two initial points"));
        }
    } else if initial.len() != vars.len() || eqs.len() < vars.len() {
        return Err(PyValueError::new_err("nsolve needs one starting point per unknown and at least as many equations as unknowns"));
    }
    let working = initial.iter().chain(&eqs).filter_map(symox::evalf_mp::max_precision)
        .max().unwrap_or(precision).max(precision) + 64;
    let starts: Vec<_> = initial.iter().map(|point| symox::evalf_mp::eval_mp(point, working)
        .ok_or_else(|| PyValueError::new_err("starting points must be finite numeric values"))).collect::<PyResult<_>>()?;
    let numerical_tolerance = tolerance.as_ref().map(|expression| {
        let value = symox::evalf_mp::eval_mp(expression, working)
            .ok_or_else(|| PyValueError::new_err("tol must be finite, real and positive"))?;
        if !value.is_real() || value.re.sign() <= 0 {
            return Err(PyValueError::new_err("tol must be finite, real and positive"));
        }
        Ok(value.re)
    }).transpose()?;
    let use_mp = solver == "secant" || digits != 15 || eqs.len() > vars.len() || starts.iter().any(|point| !point.is_real())
        || initial.iter().chain(&eqs).any(symox::evalf_mp::has_mp)
        || eqs.iter().any(|equation| symox::traverse::has(equation, &Expr::i()));
    let results: Vec<Expr> = if use_mp {
        let options = symox::nsolve_mp::Options {
            precision,
            tolerance: numerical_tolerance.unwrap_or_else(|| symox::bigfloat::BigFloat::one().ldexp(-(precision as i64) - 8)),
            maxsteps,
            verify,
        };
        let values = if solver == "bisect" {
            symox::nsolve_mp::bisect(&eqs[0], &vars[0], &starts[0], &starts[1], &options).map(|root| vec![root])
        } else if solver == "secant" {
            symox::nsolve_mp::secant(&eqs[0], &vars[0], &starts, &options).map(|root| vec![root])
        } else {
            symox::nsolve_mp::newton(&eqs, &vars, &starts, &options)
        }.map_err(PyValueError::new_err)?;
        values.into_iter().map(|value| symox::evalf_mp::to_expr(value, precision)).collect()
    } else {
        let tolerance = numerical_tolerance.map(|value| value.to_f64()).unwrap_or(1e-14);
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(PyValueError::new_err("tol is not representable at the default precision; request a higher prec"));
        }
        let starts: Vec<f64> = initial.iter().map(|point| symox::evalf::eval_f64(point)
            .ok_or_else(|| PyValueError::new_err("starting points must be finite real values"))).collect::<PyResult<_>>()?;
        let values = if solver == "bisect" {
            sv::nsolve_bisect(&eqs[0], &vars[0], starts[0], starts[1], tolerance, maxsteps, verify).map(|root| vec![root])
        } else {
            sv::nsolve_with_options(&eqs, &vars, &starts, tolerance, maxsteps, verify)
        }.ok_or_else(|| PyValueError::new_err("Could not find root within given tolerance"))?;
        values.into_iter().map(Expr::float).collect()
    };
    if as_dict {
        let mapping = PyDict::new(py);
        for (symbol, value) in vars.iter().zip(results) { mapping.set_item(wrap(symbol.clone()), wrap(value))?; }
        return Ok(PyList::new(py, [mapping])?.into_any().unbind());
    }
    if !sys && results.len() == 1 {
        return Ok(wrap(results[0].clone()).into_pyobject(py)?.into_any().unbind());
    }
    let values = PyList::new(py, results.into_iter().map(wrap))?;
    Ok(py.import("symox")?.getattr("Matrix")?.call1((vars.len(), 1, values))?.unbind())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_nonlinsolve_parametric, m)?)?;
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add_function(wrap_pyfunction!(nsolve, m)?)?;
    Ok(())
}
