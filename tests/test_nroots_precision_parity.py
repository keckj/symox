"""Clustered numerical roots and EX roots, compared with SymPy.

The roots of a polynomial with Float coefficients depend on the last bits of
the coefficients when roots are clustered: (x - 1)*(x - 1.0001)*(x - 1.0002)
has its roots moved by about 4e-8 per unit in the last place. The expansion
of Float products and powers follows the reference association, so the
polynomial is the reference one, and its roots are computed to the requested
precision (errors are measured against high precision roots, not digits).
"""
import mpmath
import pytest
import sympy as sp
import symox as sm

x_sp = sp.Symbol('x')
x_sm = sm.Symbol('x')


def _coefficients(module, expr):
    x = module.Symbol('x')
    return [str(module.Rational(c)) for c in module.Poly(expr, x).all_coeffs()]


FLOAT_PRODUCTS = [
    "(x-1)*(x-Float('1.0001'))*(x-Float('1.0002'))",
    "expand((x-1)*(x-Float('1.0001'))*(x-Float('1.0002')))",
    "(x-1)*(x-Float('1.0001'))*(x-Float('1.0002'))*(x-Float('1.0003'))",
    "(x-1)*(x-Float('1.000001'))*(x-Float('1.000002'))",
    "(x-Float('1.0001'))**3",
    "(x-Float('1.0001'))**5",
    "(x-Float('1.0001'))**2*(x-Float('1.0003'))**3",
    "(Float('1.1')*x-Float('0.3'))**4",
    "(x-Float('1.0001'))**7*(x+Float('3.3'))",
    # storage order and the reference Float ordering (positive first, by mpf tuple)
    "(x-1)*(x-1-Float('0.01'))*(x-1-2*Float('0.01'))",
    "(x-1)*(x-1-Float('1.0e-8'))*(x-1-2*Float('1.0e-8'))*(x-1-3*Float('1.0e-8'))",
    "(x-1)*(x-1-Float('0.0001'))*(x-1-2*Float('0.0001'))*(x-1-3*Float('0.0001'))",
]


@pytest.mark.parametrize("source", FLOAT_PRODUCTS)
def test_float_expansion_matches_reference_bits(source):
    expected = _coefficients(sp, eval(source, {'x': x_sp, 'Float': sp.Float, 'expand': sp.expand}))
    actual = _coefficients(sm, eval(source, {'x': x_sm, 'Float': sm.Float, 'expand': sm.expand}))
    assert actual == expected


def test_multivariate_float_power_matches_reference_bits():
    x, y = sp.symbols('x y')
    expected = [str(sp.Rational(c)) for c in sp.Poly((x + sp.Float('0.1')*y + sp.Float('0.7'))**3, x, y).coeffs()]
    X, Y = sm.symbols('x y')
    actual = [str(sm.Rational(c)) for c in sm.Poly((X + sm.Float('0.1')*Y + sm.Float('0.7'))**3, X, Y).coeffs()]
    assert actual == expected


def test_rational_of_float_is_the_exact_binary_value():
    assert sm.Rational(sm.Float(1.5)) == sm.Rational(3, 2)
    assert str(sm.Rational(sm.Float('1.0001'))) == str(sp.Rational(sp.Float('1.0001')))
    assert str(sm.Rational(sm.Float('1.0001', 30))) == str(sp.Rational(sp.Float('1.0001', 30)))
    assert sm.Rational(sm.Float(1.5), 2) == sm.Rational(3, 4)
    assert sm.Rational(1, sm.Float(2.0)) == sm.Rational(1, 2)


def _error(values, exact):
    """Largest relative error after matching each value to a distinct root."""
    pool = list(exact)
    worst = mpmath.mpf(0)
    for v in values:
        re, im = v.as_real_imag()
        z = mpmath.mpc(mpmath.mpf(str(re)), mpmath.mpf(str(im)))
        i = min(range(len(pool)), key=lambda k: abs(pool[k] - z))
        r = pool.pop(i)
        worst = max(worst, abs(r - z)/abs(r))
    return worst


def _exact_roots(coefficients, digits):
    """High precision roots of the polynomial with these exact coefficients."""
    with mpmath.workdps(4*digits + 200):
        cs = [mpmath.mpf(c.p)/c.q for c in coefficients]
        return mpmath.polyroots(cs, maxsteps=5000, extraprec=2000)


@pytest.mark.parametrize("n", [15, 30, 50])
@pytest.mark.parametrize("source", [
    "(x-1)*(x-1.0001)*(x-1.0002)",
    "expand((x-1)*(x-1.0001)*(x-1.0002))",
    "(x-1)*(x-Float('1.0001'))*(x-Float('1.0002'))*(x-Float('1.0003'))",
    "(x-1)*(x-Float('1.000001'))*(x-Float('1.000002'))",
    "(x-1)*(x-1-Float('1.0e-8'))*(x-1-2*Float('1.0e-8'))*(x-1-3*Float('1.0e-8'))",
    "(x-1)*(x-1-Float('0.000001'))*(x-1-2*Float('0.000001'))",
])
def test_float_cluster_is_as_accurate_as_reference(source, n):
    e_sp = eval(source, {'x': x_sp, 'Float': sp.Float, 'expand': sp.expand})
    e_sm = eval(source, {'x': x_sm, 'Float': sm.Float, 'expand': sm.expand})
    coefficients = [sp.Rational(c) for c in sp.Poly(e_sp, x_sp).all_coeffs()]
    exact = _exact_roots(coefficients, n)
    with mpmath.workdps(4*n + 200):
        reference = _error(sp.nroots(e_sp, n=n), exact)
        actual = _error(sm.nroots(e_sm, n=n), exact)
        assert actual <= max(3*reference, 5*mpmath.mpf(10)**-n)


def test_task_cluster_error_against_true_values():
    roots = sm.nroots((x_sm - 1)*(x_sm - 1.0001)*(x_sm - 1.0002))
    assert [str(r) for r in roots] == [str(r) for r in sp.nroots((x_sp - 1)*(x_sp - 1.0001)*(x_sp - 1.0002))]
    assert abs(float(roots[1]) - 1.0001) < 2e-12
    assert abs(float(roots[2]) - 1.0002) < 2e-12


@pytest.mark.parametrize("n", [15, 30, 50, 100])
@pytest.mark.parametrize("k,d", [(2, 12), (3, 6), (3, 10), (4, 8), (5, 4), (5, 6)])
def test_exact_real_cluster_has_requested_accuracy(k, d, n):
    step = sm.Rational(1, 10**d)
    expr = sm.expand(sm.Mul(*[x_sm - 1 - j*step for j in range(k)]))
    exact = [1 + mpmath.mpf(j)/10**d for j in range(k)]
    with mpmath.workdps(n + 50):
        exact = [1 + mpmath.mpf(j)/10**d for j in range(k)]
        assert _error(sm.nroots(expr, n=n), exact) < 5*mpmath.mpf(10)**-n


@pytest.mark.parametrize("n", [15, 30, 50])
@pytest.mark.parametrize("d", [2, 6, 8])
def test_exact_complex_cluster_has_requested_accuracy(d, n):
    step = sm.Rational(1, 10**d)*(1 + 2*sm.I)
    expr = sm.expand(sm.Mul(*[x_sm - 1 - sm.I - j*step for j in range(3)]))
    with mpmath.workdps(n + 50):
        exact = [mpmath.mpc(1, 1) + j*mpmath.mpc(1, 2)/10**d for j in range(3)]
        assert _error(sm.nroots(expr, n=n), exact) < 5*mpmath.mpf(10)**-n


@pytest.mark.parametrize("n", [15, 30])
def test_high_degree_cluster_has_requested_accuracy(n):
    step = sm.Rational(1, 10**8)
    others = [x_sm - sm.Rational(j, 7) for j in range(1, 28)]
    expr = sm.expand(sm.Mul(*others)*(x_sm - 2 - step)*(x_sm - 2 - 2*step)*(x_sm - 2 - 3*step))
    with mpmath.workdps(n + 50):
        exact = [mpmath.mpf(j)/7 for j in range(1, 28)] + [2 + j*mpmath.mpf(10)**-8 for j in (1, 2, 3)]
        assert _error(sm.nroots(expr, n=n), exact) < 5*mpmath.mpf(10)**-n


def test_exhausted_steps_still_raise_like_reference():
    from mpmath.libmp.libhyper import NoConvergence
    expr = (x_sm - 1)*(x_sm - 1.0001)*(x_sm - 1.0002)
    with pytest.raises(NoConvergence):
        sm.nroots(expr, n=30, maxsteps=5)


EX_ROOTS = [
    "x**2 - (1 + sqrt(2))*x + sqrt(2)",
    "x**2 + sqrt(2)*x + 1",
    "x**2 - sqrt(2)*x + 1",
    "x**2 + sqrt(2)*x + 1 - I",
    "x**2 + sqrt(2)*x + 1 + I",
    "x**2 - sqrt(3)*x - 6",
    "x**2 - 2*sqrt(2)*x + 2",
    "x**2 - (pi + 1)*x + pi",
    "x**2 + I*sqrt(2)*x - 1",
    "x**2 - sqrt(2)*I*x - 1",
    "x**3 - (1 + sqrt(2) + sqrt(3))*x**2 + (sqrt(2) + sqrt(3) + sqrt(6))*x - sqrt(6)",
    "x**4 - sqrt(2)*x**2 + 1/2",
    "2*x**2 - 2*sqrt(3)*x + 1",
    "x**2 - (sqrt(2) + sqrt(3))*x + sqrt(6)",
    "x**2 - sqrt(2)*x - 1",
    "x**2 + x + sqrt(2)",
    "x**2 - 2*x*sqrt(5) + 1",
    "expand((x - sqrt(2))**3)",
    "expand((x - sqrt(2))*(x - 1)*(x + 3))",
]


@pytest.mark.parametrize("source", EX_ROOTS)
def test_ex_roots_match_reference(source):
    expected = sp.roots(sp.sympify(source), x_sp)
    actual = sm.roots(sm.sympify(source), x_sm)
    assert {str(k): v for k, v in actual.items()} == {str(k): v for k, v in expected.items()}
