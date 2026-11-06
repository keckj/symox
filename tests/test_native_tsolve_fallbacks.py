"""The transcendental fallbacks of SymPy 1.14's ``_solve`` are native.

These equations used to leave ``symox._core.solve`` with
NotImplementedError, after which the Python layer ran the ``_tsolve``
steps (``solvers._tsolve_steps.tsolve_step`` through the PyHost), the
piecewise rewrite or the unrad/inversion/masking fallbacks of
``_solve_relations``. The native solver now runs the unrad fallback of
``_solve`` for radicals mixed with other functions of the unknown, then
the ``_tsolve`` steps (``_invert``, the ``Pow``/``sin``/``cos``/
``LambertW`` branches, ``rebuild(rewrite(exp))``, the bivariate block)
with its nested ``solve`` calls kept in Rust. The results are SymPy
1.14's (``audits/performance/2026-11-06-tsolve``).
"""
import pytest

from symox import (Eq, E, I, LambertW, Rational, S, Symbol, exp, log, pi,
                   solve, sqrt, symbols)
import symox
import symox._solve_relations as relations
import symox.solvers._algebraic as algebraic
import symox.solvers._tsolve_steps as steps

x, y, a = symbols('x y a')


@pytest.fixture
def no_python_fallback(monkeypatch):
    """Fail if any Python fallback of ``solve`` is reached."""
    def forbidden(*args, **kwargs):
        raise AssertionError('Python fallback reached')
    monkeypatch.setattr(relations, '_unrad_solution', forbidden)
    monkeypatch.setattr(relations, '_inverted_solution', forbidden)
    monkeypatch.setattr(relations, '_masked_solution', forbidden)
    monkeypatch.setattr(algebraic, '_unrad_resultant', forbidden)
    monkeypatch.setattr(steps, 'tsolve_step', forbidden)
    monkeypatch.setattr(symox, '_solve_piecewise', forbidden)


@pytest.mark.parametrize('eq, expected', [
    # unrad fallback of `_solve` before `_tsolve`
    (sqrt(x)*exp(x) - 1, [LambertW(2)/2]),
    (sqrt(x) - exp(x), [-LambertW(-2)/2]),
    (x**Rational(1, 3)*exp(x) - 1, [LambertW(3)/3]),
    (exp(sqrt(x)) - 3, [log(3)**2]),
    (sqrt(log(x)) - 2, [exp(4)]),
    (sqrt(x)*2**x - 1, [S.Half]),
    # Pow branch of `_tsolve`
    (x**a - 2, [2**(1/a)]),
    (x**(2*y) - 13*x**y + 36, [2**(2/y), 3**(2/y)]),
    (x**x - x, [-1, 1]),
    (x**(2/x) - 2, [2, 4]),
    (sqrt(2)**x - sqrt(sqrt(2)), [S.Half]),
    # invert, logcombine and Lambert W
    (exp(log(5)*x) - 2**x, [0]),
    (x**x - 2, [exp(LambertW(log(2)))]),
    (3*log(x) - x*log(3), [3, -3*LambertW(-log(3)/3)/log(3)]),
    (LambertW(2*x) - y, [y*exp(y)/2]),
    (a/x + exp(x/2), [2*LambertW(-a/2)]),
    (-x**2 + 2**x, [2, 4, -2*LambertW(log(2)/2)/log(2)]),
    (x*exp(x) - 1, [LambertW(1)]),
])
def test_native_tsolve(no_python_fallback, eq, expected):
    assert solve(eq, x) == expected


def test_f_pow_g_rational_candidates(no_python_fallback):
    # f(x)**g(x) == c: the rational candidates d/s of the reference
    assert solve(x**(x + 4) - 4, x) == [-2]
    assert solve((x**2 + 1)**x - 25, x) == [2]


def test_sqrt_lambert_forms(no_python_fallback):
    assert solve(x*exp(sqrt(x)) - 1, x) == [exp(-2*LambertW(S.Half))]
    assert solve(sqrt(x) + log(x), x) == [exp(-2*LambertW(S.Half))]


def test_log_constant_keeps_branch(no_python_fallback):
    # the logcombine step only combines the terms in x: the constant
    # log(-pi/2) = log(pi/2) + I*pi keeps its imaginary part
    assert solve((x**3)**(x/2) + pi/2, x) == [
        exp(LambertW(-2*log(2)/3 + 2*log(pi)/3 + 2*I*pi/3))]


def test_flags(no_python_fallback):
    eq = sqrt(x)*exp(x) - 1
    assert solve(eq, x, check=False) == [LambertW(2)/2]
    assert solve(eq, x, dict=True) == [{x: LambertW(2)/2}]
    assert solve(Eq(sqrt(x)*exp(x), 1), x) == [LambertW(2)/2]
    assert solve(x*exp(x) - 1, x, dict=True) == [{x: LambertW(1)}]


def test_bivariate_type_dummy_is_not_filtered(no_python_fallback):
    # `bivariate_type` solves for a positive Dummy `u`; the nested `_vsolve`
    # does not filter its roots by assumptions, so `u = -sqrt(2)` is kept
    eq = (x**2 - 2*x + 1).subs(x, (log(x) + 3*x)**2 - 1)
    assert solve(eq, x) == [LambertW(3*exp(-sqrt(2)))/3,
                            LambertW(3*exp(sqrt(2)))/3]


def test_take_flag_accepted():
    # `_take` is only read by the public `unrad`; `solve` ignores it
    assert solve(sqrt(x)*exp(x) - 1, x, _take=lambda d: True) == [
        LambertW(2)/2]


@pytest.mark.parametrize('eq', [
    x*(x + 2)**y - 3,
    sqrt(x) + sqrt(x + 1) - exp(x),
])
def test_not_implemented_is_kept(eq):
    # the reference raises; no empty list is returned
    with pytest.raises(NotImplementedError):
        solve(eq, x)


def test_quadratic_primitive_part(no_python_fallback):
    # `roots` takes the primitive part of the polynomial (`preprocess_roots`)
    # before `roots_quadratic`: the content `x` leaves the square root
    C1 = Symbol('C1')
    expected = [(1 - sqrt(x**2*exp(2*C1) + 1))/x,
                (sqrt(x**2*exp(2*C1) + 1) + 1)/x]
    assert solve(-x**2*exp(2*C1) + x*y*(x*y - 2), y) == expected
    # the logcombine branch of `_tsolve` reaches the same polynomial
    eq = -log(x) + log(x*y)/2 + log(x*y - 2)/2 - C1
    assert solve(eq, y, force=True) == expected
