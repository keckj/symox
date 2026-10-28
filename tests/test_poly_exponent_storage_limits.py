"""Sparse exponent limits must preserve zeros, cancellations and precision."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

N = 1 << 31
ERRORS = (NotImplementedError,)


def evaluate(module, source):
    namespace = vars(module).copy()
    namespace.update(zip(('x', 'y', 'z'), module.symbols('x y z')))
    return eval(source, namespace)


@pytest.mark.parametrize('domain', ['ZZ', 'QQ', 'ZZ_I', 'QQ_I'])
@pytest.mark.parametrize('power', [False, True])
def test_sparse_polynomial_products_and_powers_refuse_actual_overflow(domain, power):
    x = sx.Symbol('x')
    p = sx.Poly.from_dict({(N,): 1}, x, domain=evaluate(sx, domain))
    assert p.degree() == N
    with pytest.raises(ERRORS, match='exponent storage limit'):
        p**2 if power else p*p
    assert (p*0).is_zero
    assert p**0 == sx.Poly(1, x, domain=evaluate(sx, domain))
    assert p**1 == p


def test_sparse_legal_boundary_and_different_variable_exponents_remain_supported():
    x, y = sx.symbols('x y')
    p = sx.Poly.from_dict({((1 << 32)-1,): 1}, x)
    assert (p*sx.Poly(1, x)).monoms() == [((1 << 32)-1,)]
    a = sx.Poly.from_dict({(N, 0): 1}, x, y)
    b = sx.Poly.from_dict({(0, N): 1}, x, y)
    assert (a*b).monoms() == [(N, N)]


@pytest.mark.parametrize('domain', ['QQ[y]', 'QQ.frac_field(y)', 'QQ_I[y]', 'QQ_I.frac_field(y)'])
def test_polynomial_coefficient_powers_are_guarded_too(domain):
    x, y = sx.symbols('x y')
    K = evaluate(sx, domain)
    p = sx.Poly(y**N*x, x, domain=K)
    if domain.startswith('QQ_I'):
        for result in (p*p, p**2):
            assert result.as_expr() == x**2*y**(2*N)
            assert str(result.domain) == 'EX'
    else:
        with pytest.raises(ERRORS, match='exponent storage limit'):
            p*p
        with pytest.raises(ERRORS, match='exponent storage limit'):
            p**2


@pytest.mark.parametrize('domain', ['QQ[y]', 'QQ.frac_field(y)', 'QQ_I[y]', 'QQ_I.frac_field(y)',
                                     'RealField(140).poly_ring(y)', 'RealField(140).frac_field(y)',
                                     'ComplexField(140).poly_ring(y)', 'ComplexField(140).frac_field(y)'])
def test_domain_and_ring_element_operations_report_representation_errors(domain):
    y = sx.Symbol('y')
    K = evaluate(sx, domain)
    element = K.from_sympy(y**N)
    with pytest.raises(ERRORS, match='exponent storage limit'):
        element*element
    with pytest.raises(ERRORS, match='exponent storage limit'):
        element**2


@pytest.mark.parametrize('complex_', [False, True])
def test_inexact_parameter_guards_do_not_change_regular_140_bit_arithmetic(complex_):
    x, y = sx.symbols('x y')
    base = sx.ComplexField(140) if complex_ else sx.RealField(140)
    K = base.frac_field(y)
    value = sx.Float('1.2345678901234567890123456789', precision=140)*(sx.I if complex_ else 1)
    a = sx.Poly(value*y*x, x, domain=K)
    b = sx.Poly(2*y*x, x, domain=K)
    result = a*b
    assert result.domain.dom.precision == 140
    assert {f._prec for c in result.coeffs() for f in c.atoms(sx.Float)} == {140}
    exact_value = sx.to_sympy(value).xreplace({f: sp.Rational(f) for f in sx.to_sympy(value).atoms(sp.Float)})
    expression = sx.to_sympy(result.as_expr())
    expression = expression.xreplace({f: sp.Rational(f) for f in expression.atoms(sp.Float)})
    assert sp.expand(expression - 2*exact_value*sp.Symbol('y')**2*sp.Symbol('x')**2) == 0


def test_matrix_charpoly_and_products_preserve_cancelled_large_terms():
    x = sx.Symbol('x')
    t = x**N
    matrix = sx.Matrix([[t, t], [t, t]])
    assert matrix.charpoly().all_coeffs() == [1, -2*t, 0]
    assert str(matrix.charpoly().domain) == 'ZZ[x]'
    assert sx.Matrix([[t, t]])*sx.Matrix([[t], [-t]]) == sx.zeros(1)
    nilpotent = sx.Matrix([[t, t], [-t, -t]])
    assert nilpotent**2 == sx.zeros(2)
    assert sx.ones(4).applyfunc(lambda value: value*t).det() == 0


@pytest.mark.parametrize('exponent', [N, (1 << 32)+7])
def test_matrix_charpoly_can_retain_exact_expression_coefficients_beyond_sparse_storage(exponent):
    x = sx.Symbol('x')
    matrix = sx.diag(x**exponent, x**exponent)
    polynomial = matrix.charpoly()
    assert polynomial.all_coeffs() == [1, -2*x**exponent, x**(2*exponent)]
    assert str(polynomial.domain) == 'EX'
    assert matrix.inv() == sx.diag(x**(-exponent), x**(-exponent))


def test_private_storage_cancellation_does_not_emit_a_panic_or_escape_pyo3():
    source = '''import symox as s
x=s.Symbol('x'); p=s.Poly.from_dict({(2**31,):1},x)
for operation in (lambda:p*p,lambda:p**2):
    try: operation()
    except (NotImplementedError,ValueError) as error:
        assert 'exponent storage limit' in str(error)
    else: raise AssertionError('overflow returned a polynomial')
'''
    result = subprocess.run([sys.executable, '-B', '-c', source], capture_output=True, text=True, timeout=15)
    assert result.returncode == 0, result.stderr
    assert 'panicked' not in result.stderr
    assert 'PanicException' not in result.stderr


