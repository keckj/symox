//! Python bindings of the printers, `cse` and code generation.

use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::expr::Expr;
use symox::printing::{
    self as pr, CodeError, CodeSettings, CodegenOptions, CseOptions, CudaSettings, FloatType,
    InvTrigStyle, LatexMode, LatexSettings, Lang, PrettySettings,
};

use crate::{to_expr, wrap};

fn code_err(e: CodeError) -> PyErr {
    match e {
        CodeError::Unsupported(_) => PyNotImplementedError::new_err(e.to_string()),
        CodeError::PiecewiseNoDefault => PyValueError::new_err(e.to_string()),
    }
}

/// One expression or a sequence of expressions; the flag tells whether a
/// single expression was given.
fn expr_list(obj: &Bound<'_, PyAny>) -> PyResult<(Vec<Expr>, bool)> {
    if let Ok(e) = to_expr(obj) {
        return Ok((vec![e], true));
    }
    let v = obj
        .try_iter()
        .map_err(|_| PyTypeError::new_err("expected an expression or a sequence of expressions"))?
        .map(|a| to_expr(&a?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok((v, false))
}

fn float_type(dtype: &str) -> PyResult<FloatType> {
    match dtype {
        "double" | "float64" | "f64" => Ok(FloatType::F64),
        "float" | "float32" | "f32" => Ok(FloatType::F32),
        _ => Err(PyValueError::new_err(format!("unknown dtype {dtype:?} (use 'double' or 'float')"))),
    }
}

// ----------------------------------------------------------------------
// latex / pretty
// ----------------------------------------------------------------------

#[pyfunction(name = "_latex_translate")]
fn py_latex_translate(text: &str) -> String {
    pr::latex::translate(text)
}

#[pyfunction]
#[pyo3(signature = (expr, **settings))]
fn latex(expr: &Bound<'_, PyAny>, settings: Option<&Bound<'_, PyDict>>) -> PyResult<String> {
    let e = to_expr(expr)?;
    let mut s = LatexSettings::default();
    if let Some(kw) = settings {
        for (k, v) in kw.iter() {
            let key: String = k.extract()?;
            match key.as_str() {
                "mode" => {
                    let m: String = v.extract()?;
                    s.mode = match m.as_str() {
                        "plain" => LatexMode::Plain,
                        "inline" => LatexMode::Inline,
                        "equation" => LatexMode::Equation,
                        "equation*" => LatexMode::EquationStar,
                        _ => {
                            return Err(PyValueError::new_err(
                                "'mode' must be one of 'inline', 'plain', 'equation' or 'equation*'",
                            ))
                        }
                    }
                }
                "mul_symbol" => s.mul_symbol = v.extract()?,
                "fold_short_frac" => s.fold_short_frac = v.extract()?,
                "long_frac_ratio" => s.long_frac_ratio = v.extract()?,
                "fold_frac_powers" => s.fold_frac_powers = v.extract()?,
                "fold_func_brackets" => s.fold_func_brackets = v.extract()?,
                "inv_trig_style" => {
                    let m: String = v.extract()?;
                    s.inv_trig_style = match m.as_str() {
                        "abbreviated" => InvTrigStyle::Abbreviated,
                        "full" => InvTrigStyle::Full,
                        "power" => InvTrigStyle::Power,
                        _ => return Err(PyValueError::new_err("unknown inv_trig_style")),
                    }
                }
                "ln_notation" => s.ln_notation = v.extract()?,
                "root_notation" => s.root_notation = v.extract()?,
                "itex" => s.itex = v.extract()?,
                "imaginary_unit" => s.imaginary_unit = v.extract()?,
                "diff_operator" => s.diff_operator = v.extract()?,
                "parenthesize_super" => s.parenthesize_super = v.extract()?,
                "symbol_names" => {
                    let d = v.downcast::<PyDict>()?;
                    for (sym, name) in d.iter() {
                        let se = to_expr(&sym)?;
                        let n = se
                            .as_symbol()
                            .ok_or_else(|| PyTypeError::new_err("symbol_names keys must be symbols"))?
                            .to_string();
                        s.symbol_names.insert(n, name.extract()?);
                    }
                }
                "mat_symbol_style" => s.mat_symbol_bold = v.extract::<Option<String>>()?.as_deref() == Some("bold"),
                // An unknown style falls back to the dagger, as in SymPy.
                "adjoint_style" => s.adjoint_style = v.extract()?,
                "mat_delim" => s.mat_delim = v.extract::<Option<String>>()?.unwrap_or_default(),
                "mat_str" => s.mat_str = v.extract()?,
                "order" | "full_prec" | "decimal_separator" | "perm_cyclic" | "min" | "max"
                | "disable_split_super_sub" => {}
                _ => return Err(PyTypeError::new_err(format!("Unknown setting '{key}'."))),
            }
        }
    }
    crate::host::checked(pr::latex_with(&e, &s))
}

#[pyfunction]
#[pyo3(signature = (expr, use_unicode=None, wrap_line=true, num_columns=None, root_notation=true, use_unicode_sqrt_char=true, **_settings))]
#[allow(clippy::too_many_arguments)]
fn pretty(
    expr: &Bound<'_, PyAny>,
    use_unicode: Option<bool>,
    wrap_line: bool,
    num_columns: Option<usize>,
    root_notation: bool,
    use_unicode_sqrt_char: bool,
    _settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let s = PrettySettings {
        use_unicode: use_unicode.unwrap_or(true),
        wrap_line,
        num_columns: num_columns.unwrap_or(80),
        root_notation,
        use_unicode_sqrt_char,
    };
    if let Ok(matrix) = expr.extract::<PyRef<'_, crate::matrices::PyMatrix>>() {
        return Ok(pr::pretty::pretty_matrix_with(&matrix.m, &s));
    }
    let e = to_expr(expr)?;
    crate::host::checked(pr::pretty_with(&e, &s))
}

// ----------------------------------------------------------------------
// cse
// ----------------------------------------------------------------------

#[pyfunction]
#[pyo3(name = "_cse", signature = (exprs, optimizations=true, ignore=None, prefix="x", order="canonical"))]
fn py_cse<'py>(
    py: Python<'py>,
    exprs: &Bound<'py, PyAny>,
    optimizations: bool,
    ignore: Option<&Bound<'py, PyAny>>,
    prefix: &str,
    order: &str,
) -> PyResult<(Bound<'py, PyList>, Bound<'py, PyList>)> {
    let (v, _) = expr_list(exprs)?;
    let mut opts = CseOptions {
        symbol_prefix: prefix.to_string(),
        optimizations,
        canonical_order: order != "none",
        ..Default::default()
    };
    if let Some(ig) = ignore {
        for s in ig.try_iter()? {
            let e = to_expr(&s?)?;
            if let Some(n) = e.as_symbol() {
                opts.ignore.push(n.to_string());
            }
        }
    }
    let (reps, red) = py.allow_threads(|| pr::cse_with(&v, &opts));
    let reps_py = PyList::empty(py);
    for (s, e) in reps {
        reps_py.append(PyTuple::new(py, [wrap(s), wrap(e)])?)?;
    }
    let red_py = PyList::new(py, red.into_iter().map(wrap))?;
    Ok((reps_py, red_py))
}

// ----------------------------------------------------------------------
// code printers
// ----------------------------------------------------------------------

fn apply_code_kwargs(s: &mut CodeSettings, kw: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let Some(kw) = kw else { return Ok(()) };
    for (k, v) in kw.iter() {
        let key: String = k.extract()?;
        match key.as_str() {
            "user_functions" => {
                let d = v.downcast::<PyDict>()?;
                for (f, n) in d.iter() {
                    let name: String = match f.extract::<String>() {
                        Ok(n) => n,
                        Err(_) => f.getattr("__name__")?.extract()?,
                    };
                    if let Ok(target) = n.extract::<String>() {
                        s.user_functions.insert(name, target);
                    }
                }
            }
            "allow_unknown_functions" => s.allow_unknown_functions = v.extract()?,
            "strict" => s.allow_unknown_functions = !v.extract::<bool>()?,
            "real" | "dtype" | "type" => s.real = float_type(&v.extract::<String>()?)?,
            "math_macros" => s.math_macros = v.extract::<bool>().unwrap_or(true),
            "fast_math" => s.fast_math = v.extract()?,
            "use_rsqrt" => s.use_rsqrt = v.extract()?,
            "inline_pow" => s.inline_pow = v.extract()?,
            "precision" | "contract" | "human" | "full_prec" | "order" | "standard" | "fully_qualified_modules"
            | "inline" | "dereference" | "error_on_reserved" | "reserved_word_suffix" => {}
            _ => return Err(PyTypeError::new_err(format!("Unknown setting '{key}'."))),
        }
    }
    Ok(())
}

fn run_code(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    mut s: CodeSettings,
    kw: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    apply_code_kwargs(&mut s, kw)?;
    let e = to_expr(expr)?;
    let p = pr::CodePrinter::new(&s);
    let body = p.print_top(&e).map_err(code_err)?;
    let body = match assign_to {
        None => body,
        Some(a) => {
            let lhs = match a.extract::<String>() {
                Ok(name) => name,
                Err(_) => p.print(&to_expr(a)?).map_err(code_err)?,
            };
            match s.lang {
                Lang::Python | Lang::NumPy | Lang::Fortran => format!("{lhs} = {body}"),
                _ => format!("{lhs} = {body};"),
            }
        }
    };
    Ok(p.finish(body))
}

#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, standard="c99", **settings))]
fn ccode(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    standard: &str,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    if !matches!(standard.to_lowercase().as_str(), "c99" | "c11" | "c17" | "c89" | "ansi") {
        return Err(PyValueError::new_err(format!("unsupported C standard {standard:?}")));
    }
    run_code(expr, assign_to, CodeSettings::new(Lang::C), settings)
}

#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, standard="c++11", **settings))]
fn cxxcode(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    standard: &str,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let _ = standard;
    run_code(expr, assign_to, CodeSettings::new(Lang::Cxx), settings)
}

#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, **settings))]
fn pycode(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    run_code(expr, assign_to, CodeSettings::new(Lang::Python), settings)
}

fn array_code(
    expr: &Bound<'_, PyAny>,
    symbol_names: Option<&Bound<'_, PyDict>>,
    user_functions: Option<&Bound<'_, PyDict>>,
    reference: Option<String>,
    mask: Option<String>,
    expression_overrides: Option<&Bound<'_, PyDict>>,
    jax: bool,
) -> PyResult<String> {
    let expression = to_expr(expr)?;
    let mut settings = pr::TorchSettings { reference, mask, ..Default::default() };
    if let Some(names) = symbol_names {
        for (symbol, name) in names.iter() {
            let target: String = name.extract()?;
            if let Ok(symbol_name) = symbol.extract::<String>() {
                settings.symbol_names.insert(symbol_name, target);
            } else {
                let symbol = to_expr(&symbol)?;
                if !symbol.is_symbol() && !matches!(symbol.kind(), symox::expr::Kind::Func(symox::expr::Func::MatrixProxy(_), _)) {
                    return Err(PyTypeError::new_err("symbol_names keys must be symbols or strings"));
                }
                settings.exact_symbol_names.insert(symbol, target);
            }
        }
    }
    if let Some(functions) = user_functions {
        for (function, target) in functions.iter() {
            let name = function.extract::<String>().or_else(|_| function.getattr("__name__")?.extract())?;
            settings.user_functions.insert(name, target.extract()?);
        }
    }
    if let Some(overrides) = expression_overrides {
        for (expression, code) in overrides.iter() {
            settings.expression_overrides.insert(to_expr(&expression)?, code.extract()?);
        }
    }
    // Symbols with different assumptions may have the same textual name.
    // Reject those rather than silently aliasing two independent inputs.
    let mut emitted = std::collections::BTreeMap::new();
    for symbol in symox::traverse::free_symbols(&expression) {
        if settings.expression_overrides.contains_key(&symbol) { continue; }
        let matrix_name;
        let name = if let Some(name) = symbol.as_symbol() { name } else if let symox::expr::Kind::Func(symox::expr::Func::MatrixProxy(parent), _) = symbol.kind() {
            if let symox::matrixexpr::MatrixKind::Symbol(name) = parent.kind() { matrix_name = name.clone(); &matrix_name }
            else { return Err(PyTypeError::new_err("expected a matrix symbol")); }
        } else { return Err(PyTypeError::new_err("expected an input symbol")); };
        let default_name = match symbol.kind() {
            symox::expr::Kind::Sym(s) if s.is_dummy() => format!("{}_{index}", s.name(), index=s.dummy_index().unwrap()),
            _ => name.to_owned(),
        };
        let target = settings.exact_symbol_names.get(&symbol)
            .or_else(|| settings.symbol_names.get(name)).map(String::as_str).unwrap_or(&default_name);
        if let Some(previous) = emitted.insert(target.to_string(), symbol.clone()) {
            if previous != symbol {
                return Err(PyValueError::new_err(format!("ambiguous symbol name {name:?}; dummify the inputs before printing")));
            }
        }
    }
    if jax { pr::jaxcode_with(&expression, &settings) } else { pr::torchcode_with(&expression, &settings) }.map_err(code_err)
}


#[pyfunction]
#[pyo3(name = "_torchcode", signature = (expr, symbol_names=None, user_functions=None, reference=None, mask=None, expression_overrides=None))]
fn py_torchcode(
    expr: &Bound<'_, PyAny>,
    symbol_names: Option<&Bound<'_, PyDict>>,
    user_functions: Option<&Bound<'_, PyDict>>,
    reference: Option<String>,
    mask: Option<String>,
    expression_overrides: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    array_code(expr, symbol_names, user_functions, reference, mask, expression_overrides, false)
}

#[pyfunction]
#[pyo3(name = "_jaxcode", signature = (expr, symbol_names=None, user_functions=None, reference=None, mask=None, expression_overrides=None))]
fn py_jaxcode(
    expr: &Bound<'_, PyAny>,
    symbol_names: Option<&Bound<'_, PyDict>>,
    user_functions: Option<&Bound<'_, PyDict>>,
    reference: Option<String>,
    mask: Option<String>,
    expression_overrides: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    array_code(expr, symbol_names, user_functions, reference, mask, expression_overrides, true)
}


#[pyfunction]
#[pyo3(name = "_slangcode", signature = (expr, dtype="float32", symbol_names=None, user_functions=None, expression_overrides=None))]
fn py_slangcode(
    expr: &Bound<'_, PyAny>,
    dtype: &str,
    symbol_names: Option<&Bound<'_, PyDict>>,
    user_functions: Option<&Bound<'_, PyDict>>,
    expression_overrides: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let dtype = match dtype {
        "float32" | "float" => pr::SlangDtype::Float32,
        "float64" | "double" => pr::SlangDtype::Float64,
        _ => return Err(PyValueError::new_err("Slang dtype must be float32 or float64")),
    };
    let mut settings = pr::SlangSettings { dtype, ..Default::default() };
    if let Some(names) = symbol_names {
        for (symbol, name) in names.iter() {
            let target: String = name.extract()?;
            if let Ok(name) = symbol.extract::<String>() {
                settings.symbol_names.insert(name, target);
            } else {
                let symbol = to_expr(&symbol)?;
                if !symbol.is_symbol() { return Err(PyTypeError::new_err("symbol_names keys must be symbols or strings")); }
                settings.exact_symbol_names.insert(symbol, target);
            }
        }
    }
    if let Some(functions) = user_functions {
        for (function, target) in functions.iter() {
            let name = function.extract::<String>().or_else(|_| function.getattr("__name__")?.extract())?;
            settings.user_functions.insert(name, target.extract()?);
        }
    }
    if let Some(overrides) = expression_overrides {
        for (expression, code) in overrides.iter() { settings.expression_overrides.insert(to_expr(&expression)?, code.extract()?); }
    }
    pr::slangcode_with(&to_expr(expr)?, &settings).map_err(code_err)
}

#[pyfunction]
#[pyo3(name = "_slang_type")]
fn py_slang_type(expr: &Bound<'_, PyAny>) -> PyResult<&'static str> {
    pr::slang_type(&to_expr(expr)?).map_err(code_err)
}

#[pyfunction]
#[pyo3(name = "_is_slang_identifier")]
fn py_is_slang_identifier(name: &str) -> bool {
    pr::slang::is_slang_identifier(name)
}


#[pyfunction]
#[pyo3(name = "_syclcode", signature = (expr, dtype="float64", symbol_names=None, user_functions=None, expression_overrides=None))]
fn py_syclcode(
    expr: &Bound<'_, PyAny>,
    dtype: &str,
    symbol_names: Option<&Bound<'_, PyDict>>,
    user_functions: Option<&Bound<'_, PyDict>>,
    expression_overrides: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let dtype = match dtype {
        "float32" | "float" => pr::SyclDtype::Float32,
        "float64" | "double" => pr::SyclDtype::Float64,
        _ => return Err(PyValueError::new_err("SYCL dtype must be float32 or float64")),
    };
    let mut settings = pr::SyclSettings { dtype, ..Default::default() };
    if let Some(names) = symbol_names {
        for (symbol, name) in names.iter() {
            let target: String = name.extract()?;
            if let Ok(name) = symbol.extract::<String>() {
                settings.symbol_names.insert(name, target);
            } else {
                let symbol = to_expr(&symbol)?;
                if !symbol.is_symbol() { return Err(PyTypeError::new_err("symbol_names keys must be symbols or strings")); }
                settings.exact_symbol_names.insert(symbol, target);
            }
        }
    }
    if let Some(functions) = user_functions {
        for (function, target) in functions.iter() {
            let name = function.extract::<String>().or_else(|_| function.getattr("__name__")?.extract())?;
            settings.user_functions.insert(name, target.extract()?);
        }
    }
    if let Some(overrides) = expression_overrides {
        for (expression, code) in overrides.iter() { settings.expression_overrides.insert(to_expr(&expression)?, code.extract()?); }
    }
    pr::syclcode_with(&to_expr(expr)?, &settings).map_err(code_err)
}

#[pyfunction]
#[pyo3(name = "_sycl_type")]
fn py_sycl_type(expr: &Bound<'_, PyAny>) -> PyResult<&'static str> {
    pr::sycl_type(&to_expr(expr)?).map_err(code_err)
}

#[pyfunction]
#[pyo3(name = "_is_sycl_identifier")]
fn py_is_sycl_identifier(name: &str) -> bool {
    pr::sycl::is_sycl_identifier(name)
}


#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, **settings))]
fn numpycode(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    run_code(expr, assign_to, CodeSettings::new(Lang::NumPy), settings)
}

#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, **settings))]
fn rust_code(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    run_code(expr, assign_to, CodeSettings::new(Lang::Rust), settings)
}

#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, standard=77, source_format="fixed", **settings))]
fn fcode(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    standard: u32,
    source_format: &str,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let mut s = CodeSettings::new(Lang::Fortran);
    if ![66, 77, 90, 95, 2003, 2008].contains(&standard) {
        return Err(PyValueError::new_err(format!("Unknown Fortran standard: {standard}")));
    }
    s.fortran_standard = standard;
    s.fortran_free_form = match source_format {
        "fixed" => false,
        "free" => true,
        _ => return Err(PyValueError::new_err(format!("Unknown source format: {source_format}"))),
    };
    run_code(expr, assign_to, s, settings)
}

#[pyfunction]
#[pyo3(signature = (expr, assign_to=None, dtype="double", fast_math=false, use_rsqrt=true, inline_pow=4, **settings))]
#[allow(clippy::too_many_arguments)]
fn cudacode(
    expr: &Bound<'_, PyAny>,
    assign_to: Option<&Bound<'_, PyAny>>,
    dtype: &str,
    fast_math: bool,
    use_rsqrt: bool,
    inline_pow: u32,
    settings: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let cs = CudaSettings { dtype: float_type(dtype)?, fast_math, use_rsqrt, inline_pow };
    run_code(expr, assign_to, cs.to_code_settings(), settings)
}

// ----------------------------------------------------------------------
// codegen
// ----------------------------------------------------------------------

fn symbols_of(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
    let (v, _) = expr_list(obj)?;
    if let Some(bad) = v.iter().find(|a| !a.is_symbol()) {
        return Err(PyTypeError::new_err(format!("arguments must be symbols, got {bad}")));
    }
    Ok(v)
}

#[pyfunction]
#[pyo3(signature = (name, args, exprs, cse=true, dtype="double", vectorized=true, header=true))]
fn ccode_function(
    py: Python<'_>,
    name: &str,
    args: &Bound<'_, PyAny>,
    exprs: &Bound<'_, PyAny>,
    cse: bool,
    dtype: &str,
    vectorized: bool,
    header: bool,
) -> PyResult<String> {
    let a = symbols_of(args)?;
    let (e, _) = expr_list(exprs)?;
    let opts = CodegenOptions {
        cse,
        dtype: float_type(dtype)?,
        vectorized,
        include_header: header,
        ..Default::default()
    };
    py.allow_threads(|| pr::ccode_function(name, &a, &e, &opts)).map_err(code_err)
}

#[pyfunction]
#[pyo3(signature = (name, args, exprs, cse=true, dtype="double", fast_math=false, use_rsqrt=true, inline_pow=4))]
#[allow(clippy::too_many_arguments)]
fn cuda_kernel(
    py: Python<'_>,
    name: &str,
    args: &Bound<'_, PyAny>,
    exprs: &Bound<'_, PyAny>,
    cse: bool,
    dtype: &str,
    fast_math: bool,
    use_rsqrt: bool,
    inline_pow: u32,
) -> PyResult<String> {
    let a = symbols_of(args)?;
    let (e, _) = expr_list(exprs)?;
    let opts = CodegenOptions { cse, dtype: float_type(dtype)?, fast_math, use_rsqrt, inline_pow, ..Default::default() };
    py.allow_threads(|| pr::cuda_kernel(name, &a, &e, &opts)).map_err(code_err)
}

#[pyfunction]
#[pyo3(signature = (name, args, exprs, cse=true, dtype="double", fast_math=false, use_rsqrt=true, inline_pow=4))]
#[allow(clippy::too_many_arguments)]
fn cuda_device_function(
    py: Python<'_>,
    name: &str,
    args: &Bound<'_, PyAny>,
    exprs: &Bound<'_, PyAny>,
    cse: bool,
    dtype: &str,
    fast_math: bool,
    use_rsqrt: bool,
    inline_pow: u32,
) -> PyResult<String> {
    let a = symbols_of(args)?;
    let (e, _) = expr_list(exprs)?;
    let opts = CodegenOptions { cse, dtype: float_type(dtype)?, fast_math, use_rsqrt, inline_pow, ..Default::default() };
    py.allow_threads(|| pr::cuda_device_function(name, &a, &e, &opts)).map_err(code_err)
}

/// `Expr._repr_latex_`: `$\displaystyle ...$` like SymPy.
pub(crate) fn repr_latex(e: &crate::PyExprNode) -> String {
    format!("$\\displaystyle {}$", pr::latex(&e.e))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(latex, m)?)?;
    m.add_function(wrap_pyfunction!(py_latex_translate, m)?)?;
    m.add_function(wrap_pyfunction!(pretty, m)?)?;
    m.add_function(wrap_pyfunction!(py_cse, m)?)?;
    m.add_function(wrap_pyfunction!(ccode, m)?)?;
    m.add_function(wrap_pyfunction!(cxxcode, m)?)?;
    m.add_function(wrap_pyfunction!(pycode, m)?)?;
    m.add_function(wrap_pyfunction!(numpycode, m)?)?;
    m.add_function(wrap_pyfunction!(py_torchcode, m)?)?;
    m.add_function(wrap_pyfunction!(py_jaxcode, m)?)?;
    m.add_function(wrap_pyfunction!(py_slangcode, m)?)?;
    m.add_function(wrap_pyfunction!(py_slang_type, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_slang_identifier, m)?)?;
    m.add_function(wrap_pyfunction!(py_syclcode, m)?)?;
    m.add_function(wrap_pyfunction!(py_sycl_type, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_sycl_identifier, m)?)?;
    m.add_function(wrap_pyfunction!(rust_code, m)?)?;
    m.add_function(wrap_pyfunction!(fcode, m)?)?;
    m.add_function(wrap_pyfunction!(cudacode, m)?)?;
    m.add_function(wrap_pyfunction!(ccode_function, m)?)?;
    m.add_function(wrap_pyfunction!(cuda_kernel, m)?)?;
    m.add_function(wrap_pyfunction!(cuda_device_function, m)?)?;
    Ok(())
}
