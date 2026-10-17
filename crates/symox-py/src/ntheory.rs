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
    if n == BigInt::from(0) {
        d.set_item(0, 1)?;
        return Ok(d);
    }
    for (p, e) in nt::factorint(&n) {
        d.set_item(p, e)?;
    }
    if n.sign() == num_bigint::Sign::Minus {
        d.set_item(-1, 1)?;
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

fn ve(e: String) -> PyErr {
    PyValueError::new_err(e)
}

#[pyfunction]
fn nt_primerange(a: u64, b: u64) -> Vec<u64> {
    if a >= b {
        return Vec::new();
    }
    nt::primerange(a, b)
}

#[pyfunction]
fn nt_nextprime(n: BigInt, ith: u64) -> BigInt {
    nt::generate::nextprime_ith(&n, ith)
}

#[pyfunction]
fn nt_prevprime(n: BigInt) -> Option<BigInt> {
    nt::generate::prevprime(&n)
}

#[pyfunction]
fn nt_prime(n: u64) -> u64 {
    nt::prime(n)
}

#[pyfunction]
fn nt_factor_list(n: BigInt) -> Vec<(BigInt, u32)> {
    nt::factor::factor_list(&n)
}

#[pyfunction]
fn nt_multiplicity(p: BigInt, n: BigInt) -> PyResult<u64> {
    nt::factor::multiplicity(&p, &n).map_err(ve)
}

#[pyfunction]
fn nt_multiplicity_in_factorial(p: BigInt, n: BigInt) -> PyResult<BigInt> {
    nt::factor::multiplicity_in_factorial(&p, &n).map_err(ve)
}

#[pyfunction]
#[pyo3(signature = (n, s, a, retries, seed, max_steps, f))]
fn nt_pollard_rho(
    py: Python<'_>,
    n: BigInt,
    s: BigInt,
    a: BigInt,
    retries: i64,
    seed: BigInt,
    max_steps: Option<BigInt>,
    f: Option<PyObject>,
) -> PyResult<Option<BigInt>> {
    match f {
        None => nt::factor::pollard_rho::<PyErr>(&n, &s, &a, retries, &seed, max_steps.as_ref(), None),
        Some(f) => {
            let mut g = |x: &BigInt| -> PyResult<BigInt> {
                let r = f.call1(py, (x.clone(),))?;
                let v: BigInt = py.import("builtins")?.getattr("int")?.call1((r,))?.extract()?;
                Ok(v)
            };
            nt::factor::pollard_rho::<PyErr>(&n, &s, &a, retries, &seed, max_steps.as_ref(), Some(&mut g))
        }
    }
}

#[pyfunction]
fn nt_pollard_pm1(n: BigInt, b: u64, a: BigInt, retries: i64, seed: BigInt) -> Option<BigInt> {
    nt::factor::pollard_pm1(&n, b, &a, retries, &seed)
}

#[pyfunction]
fn nt_primefactors(n: BigInt) -> Vec<BigInt> {
    nt::factor::primefactors(&n)
}

#[pyfunction]
fn nt_divisors(n: BigInt, proper: bool, sorted: bool) -> Vec<BigInt> {
    let mut v = nt::factor::divisors_unsorted(&n, proper);
    if sorted {
        v.sort();
    }
    v
}

#[pyfunction]
fn nt_divisor_sigma(n: BigInt, k: u32) -> BigInt {
    nt::factor::divisor_sigma_raw(&n, k)
}

#[pyfunction]
fn nt_abundance(n: BigInt) -> BigInt {
    nt::factor::abundance(&n)
}

#[pyfunction]
fn nt_is_mersenne_prime(n: BigInt) -> bool {
    nt::factor::is_mersenne_prime(&n)
}

#[pyfunction]
fn nt_is_perfect(n: BigInt) -> PyResult<bool> {
    nt::factor::is_perfect(&n).map_err(ve)
}

/// `True`/`False`, or `None` for SymPy's integer `0` answer.
#[pyfunction]
fn nt_is_carmichael(n: BigInt) -> Option<bool> {
    nt::factor::is_carmichael(&n)
}

#[pyfunction]
fn nt_perfect_power(n: BigInt, candidates: Option<Vec<BigInt>>, big: bool) -> Option<(BigInt, u64)> {
    nt::factor::perfect_power(&n, candidates.as_deref(), big)
}

#[pyfunction]
fn nt_iroot(n: BigInt, k: u32) -> (BigInt, bool) {
    nt::factor::iroot(&n, k)
}

#[pyfunction]
fn nt_mersenne_exponents() -> Vec<u64> {
    nt::factor::MERSENNE_PRIME_EXPONENTS.to_vec()
}

#[pyfunction]
fn nt_reduced_totient(n: BigInt) -> BigInt {
    symox::special::numtheory::reduced_totient(&n)
}

#[pyfunction]
fn nt_jacobi(a: BigInt, n: BigInt) -> i32 {
    nt::residue::jacobi(&a, &n)
}

#[pyfunction]
fn nt_kronecker(a: BigInt, n: BigInt) -> i32 {
    nt::residue::kronecker(&a, &n)
}

#[pyfunction]
fn nt_n_order(a: BigInt, n: BigInt) -> PyResult<BigInt> {
    nt::residue::n_order(&a, &n).map_err(ve)
}

#[pyfunction]
fn nt_is_primitive_root(a: BigInt, p: BigInt) -> PyResult<bool> {
    nt::residue::is_primitive_root(&a, &p).map_err(ve)
}

#[pyfunction]
fn nt_primitive_root(p: BigInt, smallest: bool) -> PyResult<Option<BigInt>> {
    nt::residue::primitive_root(&p, smallest).map_err(ve)
}

#[pyfunction]
fn nt_sqrt_mod_iter(a: BigInt, p: BigInt) -> Vec<BigInt> {
    nt::residue::sqrt_mod_iter(&a, &p)
}

#[pyfunction]
fn nt_sqrt_mod(a: BigInt, p: BigInt) -> Option<BigInt> {
    nt::residue::sqrt_mod(&a, &p)
}

#[pyfunction]
fn nt_is_quad_residue(a: BigInt, p: BigInt) -> PyResult<bool> {
    nt::residue::is_quad_residue(&a, &p).map_err(ve)
}

#[pyfunction]
fn nt_is_nthpow_residue(a: BigInt, n: BigInt, m: BigInt) -> PyResult<bool> {
    nt::residue::is_nthpow_residue(&a, &n, &m).map_err(ve)
}

#[pyfunction]
fn nt_nthroot_mod_all(a: BigInt, n: BigInt, p: BigInt) -> PyResult<Vec<BigInt>> {
    nt::residue::nthroot_mod_all(&a, &n, &p).map_err(ve)
}

#[pyfunction]
fn nt_quadratic_residues(p: BigInt) -> Vec<BigInt> {
    nt::residue::quadratic_residues(&p)
}

#[pyfunction]
fn nt_quadratic_congruence(a: BigInt, b: BigInt, c: BigInt, n: BigInt) -> PyResult<(Vec<BigInt>, bool)> {
    nt::residue::quadratic_congruence(&a, &b, &c, &n).map_err(ve)
}

#[pyfunction]
#[pyo3(signature = (n, a, b, order=None, prime_order=None))]
fn nt_discrete_log(n: BigInt, a: BigInt, b: BigInt, order: Option<BigInt>, prime_order: Option<bool>) -> PyResult<BigInt> {
    nt::residue::discrete_log(&n, &a, &b, order, prime_order).map_err(ve)
}

#[pyfunction]
fn nt_igcdex(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    nt::residue::igcdex(&a, &b)
}

#[pyfunction]
fn nt_primorial(n: u64, nth: bool) -> BigInt {
    nt::generate::primorial(n, nth)
}

#[pyfunction]
fn nt_compositepi(n: u64) -> u64 {
    nt::generate::compositepi(n)
}

#[pyfunction]
fn nt_composite(n: u64) -> u64 {
    nt::generate::composite(n)
}

#[pyfunction]
fn nt_totient_range(lo: u64, hi: u64) -> Vec<u64> {
    nt::generate::totient_range(lo, hi)
}

#[pyfunction]
fn nt_mobius_range(lo: u64, hi: u64) -> Vec<i64> {
    nt::generate::mobius_range(lo, hi)
}

#[pyfunction]
fn nt_integer_log(n: BigInt, b: BigInt) -> PyResult<(u64, bool)> {
    nt::misc::integer_log(&n, &b).map_err(ve)
}

#[pyfunction]
fn nt_egyptian_fraction(x: BigInt, y: BigInt, algorithm: &str) -> PyResult<Vec<BigInt>> {
    nt::misc::egyptian_fraction(&x, &y, algorithm).map_err(ve)
}

#[pyfunction]
fn nt_cf_rational(p: BigInt, q: BigInt) -> Vec<BigInt> {
    nt::misc::cf_rational(&p, &q)
}

#[pyfunction]
fn nt_cf_periodic(p: BigInt, q: BigInt, d: BigInt, s: BigInt) -> PyResult<(Vec<BigInt>, Option<usize>)> {
    nt::misc::continued_fraction_periodic(&p, &q, &d, &s).map_err(ve)
}

#[pyfunction]
fn nt_binomial_row(n: u64) -> Vec<BigInt> {
    nt::misc::binomial_row(n)
}

#[pyfunction]
fn nt_npartitions(n: u64) -> BigInt {
    symox::special::integers::partition_number(n)
}

fn ve_or_te(e: String) -> PyErr {
    if e.starts_with("'float'") {
        PyTypeError::new_err(e)
    } else {
        PyValueError::new_err(e)
    }
}

#[pyfunction]
fn nt_prp(kind: &str, n: BigInt, a: BigInt) -> PyResult<bool> {
    use nt::primetest as pt;
    match kind {
        "fermat" => pt::is_fermat_pseudoprime(&n, &a),
        "euler" => pt::is_euler_pseudoprime(&n, &a),
        "euler_jacobi" => pt::is_euler_jacobi_pseudoprime(&n, &a),
        _ => Err(format!("unknown test {kind}")),
    }
    .map_err(ve)
}

#[pyfunction]
fn nt_mr(n: BigInt, bases: Vec<BigInt>) -> bool {
    nt::primetest::mr(&n, &bases)
}

#[pyfunction]
fn nt_lucas_prp(kind: &str, n: BigInt) -> PyResult<bool> {
    use nt::primetest as pt;
    match kind {
        "lucas" => pt::is_lucas_prp(&n),
        "strong" => pt::is_strong_lucas_prp(&n),
        "extra_strong" => Ok(pt::is_extra_strong_lucas_prp(&n)),
        _ => Err(format!("unknown test {kind}")),
    }
    .map_err(ve)
}

#[pyfunction]
fn nt_proth_test(n: BigInt) -> PyResult<bool> {
    nt::primetest::proth_test(&n).map_err(ve)
}

#[pyfunction]
fn nt_csolve(coeffs: Vec<BigInt>, m: BigInt) -> PyResult<Vec<BigInt>> {
    nt::congruence::csolve(&coeffs, &m).map_err(ve)
}

#[pyfunction]
fn nt_binomial_mod(n: BigInt, m: BigInt, k: BigInt) -> PyResult<BigInt> {
    nt::congruence::binomial_mod(&n, &m, &k).map_err(ve)
}

#[pyfunction]
fn nt_ecm(n: BigInt, b1: BigInt, b2: BigInt, max_curve: u64, seed: BigInt) -> PyResult<Vec<BigInt>> {
    nt::ecm::ecm(&n, &b1, &b2, max_curve, &seed).map_err(ve)
}

type PointT = (BigInt, BigInt, BigInt, BigInt);

fn to_point(p: PointT) -> nt::ecm::Point {
    nt::ecm::Point::new(p.0, p.1, p.2, p.3)
}

fn from_point(p: nt::ecm::Point) -> PointT {
    (p.x, p.z, p.a24, p.m)
}

#[pyfunction]
fn nt_point_add(p: PointT, q: PointT, diff: PointT) -> PointT {
    from_point(to_point(p).add(&to_point(q), &to_point(diff)))
}

#[pyfunction]
fn nt_point_double(p: PointT) -> PointT {
    from_point(to_point(p).double())
}

#[pyfunction]
fn nt_point_ladder(p: PointT, k: BigInt) -> PointT {
    from_point(to_point(p).mont_ladder(&k))
}

#[pyfunction]
fn nt_qs_factor(n: BigInt, prime_bound: u64, m: u64, error_term: f64, seed: BigInt) -> PyResult<Vec<(BigInt, u32)>> {
    nt::qs::qs_factor(&n, prime_bound, m, error_term, &seed).map_err(ve_or_te)
}

#[pyfunction]
fn nt_pi_hex_digits(n: u64, prec: u64) -> String {
    nt::bbp::pi_hex_digits(n, prec)
}

#[pyfunction]
fn nt_egypt(algorithm: &str, x: BigInt, y: BigInt) -> PyResult<Vec<BigInt>> {
    use nt::misc as ms;
    Ok(match algorithm {
        "Greedy" => ms::egypt_greedy(&x, &y),
        "Graham Jewett" => ms::egypt_graham_jewett(&x, &y),
        "Takenouchi" => ms::egypt_takenouchi(&x, &y),
        "Golomb" => ms::egypt_golomb(&x, &y),
        _ => return Err(PyValueError::new_err("Entered invalid algorithm")),
    })
}

#[pyfunction]
fn nt_egypt_harmonic(x: BigInt, y: BigInt) -> (Vec<BigInt>, BigInt, BigInt) {
    nt::misc::egypt_harmonic(&x, &y)
}

#[pyfunction]
fn nt_ntt(seq: Vec<BigInt>, p: BigInt, inverse: bool) -> PyResult<Vec<BigInt>> {
    nt::discrete::ntt(&seq, &p, inverse).map_err(ve)
}

#[pyfunction]
fn nt_fwht(seq: Vec<BigInt>) -> Vec<BigInt> {
    nt::discrete::fwht(&seq)
}

#[pyfunction]
fn nt_mobius_transform(seq: Vec<BigInt>, sgn: i32, subset: bool) -> Vec<BigInt> {
    nt::discrete::mobius(&seq, sgn, subset)
}

#[pyfunction]
fn nt_linrec(coeffs: Vec<BigInt>, init: Vec<BigInt>, n: BigInt) -> BigInt {
    nt::discrete::linrec(&coeffs, &init, &n)
}

#[pyfunction]
fn nt_linrec_coeffs(coeffs: Vec<BigInt>, n: BigInt) -> Vec<BigInt> {
    nt::discrete::linrec_coeffs(&coeffs, &n)
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
    m.add_function(wrap_pyfunction!(nt_primerange, m)?)?;
    m.add_function(wrap_pyfunction!(nt_nextprime, m)?)?;
    m.add_function(wrap_pyfunction!(nt_prevprime, m)?)?;
    m.add_function(wrap_pyfunction!(nt_prime, m)?)?;
    m.add_function(wrap_pyfunction!(nt_factor_list, m)?)?;
    m.add_function(wrap_pyfunction!(nt_multiplicity, m)?)?;
    m.add_function(wrap_pyfunction!(nt_multiplicity_in_factorial, m)?)?;
    m.add_function(wrap_pyfunction!(nt_pollard_rho, m)?)?;
    m.add_function(wrap_pyfunction!(nt_pollard_pm1, m)?)?;
    m.add_function(wrap_pyfunction!(nt_primefactors, m)?)?;
    m.add_function(wrap_pyfunction!(nt_divisors, m)?)?;
    m.add_function(wrap_pyfunction!(nt_divisor_sigma, m)?)?;
    m.add_function(wrap_pyfunction!(nt_abundance, m)?)?;
    m.add_function(wrap_pyfunction!(nt_is_mersenne_prime, m)?)?;
    m.add_function(wrap_pyfunction!(nt_is_perfect, m)?)?;
    m.add_function(wrap_pyfunction!(nt_is_carmichael, m)?)?;
    m.add_function(wrap_pyfunction!(nt_perfect_power, m)?)?;
    m.add_function(wrap_pyfunction!(nt_iroot, m)?)?;
    m.add_function(wrap_pyfunction!(nt_mersenne_exponents, m)?)?;
    m.add_function(wrap_pyfunction!(nt_reduced_totient, m)?)?;
    m.add_function(wrap_pyfunction!(nt_jacobi, m)?)?;
    m.add_function(wrap_pyfunction!(nt_kronecker, m)?)?;
    m.add_function(wrap_pyfunction!(nt_n_order, m)?)?;
    m.add_function(wrap_pyfunction!(nt_is_primitive_root, m)?)?;
    m.add_function(wrap_pyfunction!(nt_primitive_root, m)?)?;
    m.add_function(wrap_pyfunction!(nt_sqrt_mod_iter, m)?)?;
    m.add_function(wrap_pyfunction!(nt_sqrt_mod, m)?)?;
    m.add_function(wrap_pyfunction!(nt_is_quad_residue, m)?)?;
    m.add_function(wrap_pyfunction!(nt_is_nthpow_residue, m)?)?;
    m.add_function(wrap_pyfunction!(nt_nthroot_mod_all, m)?)?;
    m.add_function(wrap_pyfunction!(nt_quadratic_residues, m)?)?;
    m.add_function(wrap_pyfunction!(nt_quadratic_congruence, m)?)?;
    m.add_function(wrap_pyfunction!(nt_discrete_log, m)?)?;
    m.add_function(wrap_pyfunction!(nt_igcdex, m)?)?;
    m.add_function(wrap_pyfunction!(nt_primorial, m)?)?;
    m.add_function(wrap_pyfunction!(nt_compositepi, m)?)?;
    m.add_function(wrap_pyfunction!(nt_composite, m)?)?;
    m.add_function(wrap_pyfunction!(nt_totient_range, m)?)?;
    m.add_function(wrap_pyfunction!(nt_mobius_range, m)?)?;
    m.add_function(wrap_pyfunction!(nt_integer_log, m)?)?;
    m.add_function(wrap_pyfunction!(nt_egyptian_fraction, m)?)?;
    m.add_function(wrap_pyfunction!(nt_cf_rational, m)?)?;
    m.add_function(wrap_pyfunction!(nt_cf_periodic, m)?)?;
    m.add_function(wrap_pyfunction!(nt_binomial_row, m)?)?;
    m.add_function(wrap_pyfunction!(nt_npartitions, m)?)?;
    m.add_function(wrap_pyfunction!(nt_prp, m)?)?;
    m.add_function(wrap_pyfunction!(nt_mr, m)?)?;
    m.add_function(wrap_pyfunction!(nt_lucas_prp, m)?)?;
    m.add_function(wrap_pyfunction!(nt_proth_test, m)?)?;
    m.add_function(wrap_pyfunction!(nt_csolve, m)?)?;
    m.add_function(wrap_pyfunction!(nt_binomial_mod, m)?)?;
    m.add_function(wrap_pyfunction!(nt_ecm, m)?)?;
    m.add_function(wrap_pyfunction!(nt_point_add, m)?)?;
    m.add_function(wrap_pyfunction!(nt_point_double, m)?)?;
    m.add_function(wrap_pyfunction!(nt_point_ladder, m)?)?;
    m.add_function(wrap_pyfunction!(nt_qs_factor, m)?)?;
    m.add_function(wrap_pyfunction!(nt_pi_hex_digits, m)?)?;
    m.add_function(wrap_pyfunction!(nt_egypt, m)?)?;
    m.add_function(wrap_pyfunction!(nt_egypt_harmonic, m)?)?;
    m.add_function(wrap_pyfunction!(nt_ntt, m)?)?;
    m.add_function(wrap_pyfunction!(nt_fwht, m)?)?;
    m.add_function(wrap_pyfunction!(nt_mobius_transform, m)?)?;
    m.add_function(wrap_pyfunction!(nt_linrec, m)?)?;
    m.add_function(wrap_pyfunction!(nt_linrec_coeffs, m)?)?;
    Ok(())
}
