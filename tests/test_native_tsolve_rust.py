"""The native ``_tsolve`` steps (``symox::solve::tsolve``) and the
``sympy.solvers.bivariate`` helpers, against values computed with SymPy
1.14 (hard-coded below)."""
from symox import *  # noqa: F403
from symox import _core
from symox.solvers.bivariate import (
    _filtered_gens, _lambert, _linab, _mexpand, _mostfunc, _solve_lambert,
    bivariate_type)
from symox.solvers._tsolve_steps import _invert
import pytest

del test  # noqa: F821

x, y, a, z = symbols('x y a z')


def _set(v):
    return set(v)


@pytest.mark.parametrize('eq, expected', [
    # f(x)**e = rhs with e free of x
    (lambda: x**y - 2, lambda: [2**(1/y)]),
    ((lambda: (x**2 - 1)**y - 4),
     lambda: [-sqrt(2**(2/y) + 1), sqrt(2**(2/y) + 1)]),
    (lambda: (x + 1)**a - 2, lambda: [2**(1/a) - 1]),
    # Lambert W forms through the bivariate block
    (lambda: x*exp(x**2) - 1, lambda: [exp(-LambertW(2)/2)]),
    (lambda: sqrt(x)*exp(sqrt(x)) - 1, lambda: [exp(-2*LambertW(1))]),
    (lambda: x - exp(-x**2), lambda: [exp(-LambertW(2)/2)]),
    (lambda: x**a*exp(x) - 1, lambda: [a*LambertW(1/a)]),
    (lambda: x**a - exp(x), lambda: [-a*LambertW(-1/a)]),
    (lambda: sqrt(x) + log(x), lambda: [exp(-2*LambertW(S(1)/2))]),
])
def test_solve_through_tsolve(eq, expected):
    assert _set(solve(eq(), x)) == _set(expected())


def test_checksol_rejects_tsolve_candidates():
    # SymPy test_issue_4463: x = 0 makes log(x) infinite
    assert solve(-a*x + 2*x*log(x), x) == [exp(a/2)]
    assert solve(x**x - 2) == [exp(LambertW(log(2)))]
    # check=False keeps the candidate, like the reference
    assert _set(solve(-a*x + 2*x*log(x), x, check=False)) == {0, exp(a/2)}


def test_uselogcombine_force():
    eq = z - log(x) + log(y/(x*(-1 + y**2/x**2)))
    assert _set(solve(eq, x, force=True)) == {
        -sqrt(y*(y - exp(z))), sqrt(y*(y - exp(z)))}


def test_bivariate_helpers():
    assert _filtered_gens(Poly(x + 1/x + exp(x)), x) == {x, exp(x)}
    assert _mostfunc(exp(x) + exp(exp(x) + 2), exp) == exp(exp(x) + 2)
    assert _mostfunc(exp(x) + exp(exp(y) + 2), exp, x) == exp(x)
    assert _linab(x/y - 3/y, x) == (1/y, -3/y, x)
    assert _linab(2*y*exp(x) + 3, x) == (2*y, 3, exp(x))
    assert _lambert(x + log(x) - 2, x) == [LambertW(exp(2))]
    assert _solve_lambert(x*exp(x) - 1, x, [x, exp(x)]) == [LambertW(1)]
    with pytest.raises(NotImplementedError):
        _solve_lambert(x**2 - 2**x, x, [x, 2**x])
    g, p, u = bivariate_type(x*y + 2*x*y + 1, x, y)
    assert (g, p) == (x*y, 3*u + 1)
    g, p, u = bivariate_type((x + y)**2 + 3*(x + y), x, y)
    assert (g, p) == (x + y, u**2 + 3*u)
    assert u.is_positive
    assert _mexpand((x + 1)**2*(x - 1)) == x**3 + x**2 - x - 1


def test_invert():
    assert _invert(exp(x) - 3, x) == (log(3), x)
    assert _invert(log(x + 1) - y, x) == (exp(y) - 1, x)
    rhs, lhs = _invert(atan2(x, 2) - pi/3, x)
    assert rhs == sqrt(3)/3
    assert lhs == x/(sqrt(x**2 + 4) + 2)


def test_native_entry_points():
    for name in ('_tsolve_step', '_tsolve_invert', '_bivariate_solve_lambert',
                 '_bivariate_type', '_bivariate_lambert'):
        assert hasattr(_core, name)


def test_invert_free_of_symbol_returns_early():
    # SymPy's ``_invert`` returns ``(eq, 0)`` for an equation free of the
    # symbol; the native loop used to spin on ``0`` (test_issue_24368)
    from symox import S, Symbol, _core
    x = Symbol('x')
    assert _core._tsolve_invert(S.Zero, x) == (S.Zero, S.Zero)
    assert _core._tsolve_invert(S(3), x) == (S(3), S.Zero)


def test_issue_24368_declines():
    from symox import Mod, S, Symbol, floor, solve
    s2 = Symbol('s2', integer=True, positive=True)
    f = floor(s2/2 - S(1)/2)
    g = f**2/(f + 1) + 2*f/(f + 1) + 1/(f + 1)
    with pytest.raises(NotImplementedError):
        solve(Mod(g, 1)*f + Mod(g, 1), s2)
