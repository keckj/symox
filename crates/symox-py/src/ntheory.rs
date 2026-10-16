//! Python bindings for number theory functions.

use num_bigint::BigInt;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use symox::ntheory as nt;

use crate::to_expr;

fn int_arg(o: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    if let Ok(b) = o.extract::<BigInt>() {
        return Ok(b);
    }
    let e = to_expr(o)?;
    e.as_number()
        .and_then(|n| n.to_bigint())
        .ok_or_else(|| PyTypeError::new_err(format!("expected an integer, got {e}")))
}

#[pyfunction]
fn isprime(n: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(nt::isprime(&int_arg(n)?))
}

#[pyfunction]
fn factorint<'py>(py: Python<'py>, n: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyDict>> {
    let n = int_arg(n)?;
    let d = PyDict::new(py);
    if n.sign() == num_bigint::Sign::Minus {
        d.set_item(-1, 1)?;
    }
    if n == BigInt::from(0) {
        d.set_item(0, 1)?;
        return Ok(d);
    }
    for (p, e) in nt::factorint(&n) {
        d.set_item(p, e)?;
    }
    Ok(d)
}

#[pyfunction]
#[pyo3(signature = (a, b=None))]
fn primerange(a: u64, b: Option<u64>) -> Vec<u64> {
    match b {
        Some(b) => nt::primerange(a, b),
        None => nt::primerange(2, a),
    }
}

#[pyfunction]
fn nextprime(n: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    Ok(nt::nextprime(&int_arg(n)?))
}

#[pyfunction]
fn prevprime(n: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    nt::prevprime(&int_arg(n)?).ok_or_else(|| PyValueError::new_err("no preceding primes"))
}

#[pyfunction]
fn prime(n: u64) -> PyResult<u64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be >= 1"));
    }
    Ok(nt::prime(n))
}

#[pyfunction]
fn primepi(n: u64) -> u64 {
    nt::primepi(n)
}

#[pyfunction]
fn totient(n: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    Ok(nt::totient(&int_arg(n)?))
}

#[pyfunction]
fn divisors(n: &Bound<'_, PyAny>) -> PyResult<Vec<BigInt>> {
    Ok(nt::divisors(&int_arg(n)?))
}

#[pyfunction]
fn divisor_count(n: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    Ok(nt::divisor_count(&int_arg(n)?))
}

#[pyfunction]
#[pyo3(signature = (n, k=1))]
fn divisor_sigma(n: &Bound<'_, PyAny>, k: u32) -> PyResult<BigInt> {
    Ok(nt::divisor_sigma(&int_arg(n)?, k))
}

#[pyfunction]
fn mod_inverse(a: &Bound<'_, PyAny>, m: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    nt::mod_inverse(&int_arg(a)?, &int_arg(m)?).ok_or_else(|| PyValueError::new_err("inverse does not exist"))
}

#[pyfunction]
fn jacobi_symbol(a: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<i32> {
    nt::jacobi_symbol(&int_arg(a)?, &int_arg(n)?).ok_or_else(|| PyValueError::new_err("n should be an odd positive integer"))
}

#[pyfunction]
fn legendre_symbol(a: &Bound<'_, PyAny>, p: &Bound<'_, PyAny>) -> PyResult<i32> {
    nt::legendre_symbol(&int_arg(a)?, &int_arg(p)?).ok_or_else(|| PyValueError::new_err("p should be an odd prime"))
}

#[pyfunction]
fn mobius(n: &Bound<'_, PyAny>) -> PyResult<i32> {
    Ok(nt::mobius(&int_arg(n)?))
}

#[pyfunction]
fn integer_nthroot(n: &Bound<'_, PyAny>, k: u32) -> PyResult<(BigInt, bool)> {
    Ok(nt::integer_nthroot(&int_arg(n)?, k))
}

#[pyfunction]
fn perfect_power(n: &Bound<'_, PyAny>) -> PyResult<Option<(BigInt, u32)>> {
    Ok(nt::perfect_power(&int_arg(n)?))
}

#[pyfunction]
fn n_order(a: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    nt::n_order(&int_arg(a)?, &int_arg(n)?).ok_or_else(|| PyValueError::new_err("a and n must be coprime"))
}

#[pyfunction]
fn primitive_root(p: &Bound<'_, PyAny>) -> PyResult<Option<BigInt>> {
    Ok(nt::primitive_root(&int_arg(p)?))
}

#[pyfunction]
fn crt(m: Vec<BigInt>, v: Vec<BigInt>) -> Option<(BigInt, BigInt)> {
    nt::crt(&v, &m)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(isprime, m)?)?;
    m.add_function(wrap_pyfunction!(factorint, m)?)?;
    m.add_function(wrap_pyfunction!(primerange, m)?)?;
    m.add_function(wrap_pyfunction!(nextprime, m)?)?;
    m.add_function(wrap_pyfunction!(prevprime, m)?)?;
    m.add_function(wrap_pyfunction!(prime, m)?)?;
    m.add_function(wrap_pyfunction!(primepi, m)?)?;
    m.add_function(wrap_pyfunction!(totient, m)?)?;
    m.add_function(wrap_pyfunction!(divisors, m)?)?;
    m.add_function(wrap_pyfunction!(divisor_count, m)?)?;
    m.add_function(wrap_pyfunction!(divisor_sigma, m)?)?;
    m.add_function(wrap_pyfunction!(mod_inverse, m)?)?;
    m.add_function(wrap_pyfunction!(jacobi_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(legendre_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(mobius, m)?)?;
    m.add_function(wrap_pyfunction!(integer_nthroot, m)?)?;
    m.add_function(wrap_pyfunction!(perfect_power, m)?)?;
    m.add_function(wrap_pyfunction!(n_order, m)?)?;
    m.add_function(wrap_pyfunction!(primitive_root, m)?)?;
    m.add_function(wrap_pyfunction!(crt, m)?)?;
    Ok(())
}
