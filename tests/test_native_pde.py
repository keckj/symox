"""Native port of SymPy 1.14 ``sympy/solvers/pde.py`` (``symox::ode::pde``).

The first part transcribes ``sympy/solvers/tests/test_pde.py``; the second
compares the results live with SymPy and checks that the public callables
reach the native port without the reference Python body.
"""
import pytest
import sympy as sp

import symox as sx
from symox import _core
from symox import Derivative as D, Eq, Function, S, Symbol, cos, exp, log, sin, symbols
from symox.solvers import pde
from symox.solvers.pde import (checkpdesol, classify_pde, pde_separate,
    pde_separate_add, pde_separate_mul, pdsolve)

raises = pytest.raises

a, b, c, x, y = symbols('a b c x y')


# --- sympy/solvers/tests/test_pde.py ---------------------------------------

def test_pde_separate_add():
    x, y, z, t = symbols("x,y,z,t")
    F, T, X, Y, Z, u = map(Function, 'FTXYZu')

    eq = Eq(D(u(x, t), x), D(u(x, t), t)*exp(u(x, t)))
    res = pde_separate_add(eq, u(x, t), [X(x), T(t)])
    assert res == [D(X(x), x)*exp(-X(x)), D(T(t), t)*exp(T(t))]


def test_pde_separate():
    x, y, z, t = symbols("x,y,z,t")
    F, T, X, Y, Z, u = map(Function, 'FTXYZu')

    eq = Eq(D(u(x, t), x), D(u(x, t), t)*exp(u(x, t)))
    with raises(ValueError):
        pde_separate(eq, u(x, t), [X(x), T(t)], 'div')


def test_pde_separate_mul():
    x, y, z, t = symbols("x,y,z,t")
    c = Symbol("C", real=True)
    Phi = Function('Phi')
    F, R, T, X, Y, Z, u = map(Function, 'FRTXYZu')
    r, theta, z = symbols('r,theta,z')

    eq = Eq(D(F(x, y, z), x) + D(F(x, y, z), y) + D(F(x, y, z), z), 0)

    with raises(ValueError):
        pde_separate_mul(eq, F(x, y, z), [X(x), u(z, z)])
    with raises(ValueError):
        pde_separate_mul(eq, F(x, y, z), [X(x), Y(y)])
    with raises(ValueError):
        pde_separate_mul(eq, F(x, y, z), [X(t), Y(x, y)])

    assert pde_separate_mul(eq, F(x, y, z), [Y(y), u(x, z)]) == \
        [D(Y(y), y)/Y(y), -D(u(x, z), x)/u(x, z) - D(u(x, z), z)/u(x, z)]
    assert pde_separate_mul(eq, F(x, y, z), [X(x), Y(y), Z(z)]) == \
        [D(X(x), x)/X(x), -D(Z(z), z)/Z(z) - D(Y(y), y)/Y(y)]

    wave = Eq(D(u(x, t), t, t), c**2*D(u(x, t), x, x))
    res = pde_separate_mul(wave, u(x, t), [X(x), T(t)])
    assert res == [D(X(x), x, x)/X(x), D(T(t), t, t)/(c**2*T(t))]

    eq = Eq(1/r * D(Phi(r, theta, z), r) + D(Phi(r, theta, z), r, 2) +
            1/r**2 * D(Phi(r, theta, z), theta, 2) + D(Phi(r, theta, z), z, 2), 0)
    res = pde_separate_mul(eq, Phi(r, theta, z), [Z(z), u(theta, r)])
    assert res == [D(Z(z), z, z)/Z(z),
            -D(u(theta, r), r, r)/u(theta, r) -
        D(u(theta, r), r)/(r*u(theta, r)) -
        D(u(theta, r), theta, theta)/(r**2*u(theta, r))]
    eq = Eq(res[1], c)
    res = pde_separate_mul(eq, u(theta, r), [T(theta), R(r)])
    assert res == [D(T(theta), theta, theta)/T(theta),
            -r*D(R(r), r)/R(r) - r**2*D(R(r), r, r)/R(r) - c*r**2]
    res = pde_separate_mul(eq, u(theta, r), [R(r), T(theta)])
    assert res == [r*D(R(r), r)/R(r) + r**2*D(R(r), r, r)/R(r) + c*r**2,
            -D(T(theta), theta, theta)/T(theta)]


def test_issue_11726():
    x, t = symbols("x t")
    f = symbols("f", cls=Function)
    X, T = symbols("X T", cls=Function)

    u = f(x, t)
    eq = u.diff(x, 2) - u.diff(t, 2)
    res = pde_separate(eq, u, [T(x), X(t)])
    assert res == [D(T(x), x, x)/T(x), D(X(t), t, t)/X(t)]


def test_pde_classify():
    f = Function('f')
    eq1 = a*f(x,y) + b*f(x,y).diff(x) + c*f(x,y).diff(y)
    eq2 = 3*f(x,y) + 2*f(x,y).diff(x) + f(x,y).diff(y)
    eq3 = a*f(x,y) + b*f(x,y).diff(x) + 2*f(x,y).diff(y)
    eq4 = x*f(x,y) + f(x,y).diff(x) + 3*f(x,y).diff(y)
    eq5 = x**2*f(x,y) + x*f(x,y).diff(x) + x*y*f(x,y).diff(y)
    eq6 = y*x**2*f(x,y) + y*f(x,y).diff(x) + f(x,y).diff(y)
    for eq in [eq1, eq2, eq3]:
        assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    for eq in [eq4, eq5, eq6]:
        assert classify_pde(eq) == ('1st_linear_variable_coeff',)


def test_checkpdesol():
    f, F = map(Function, ['f', 'F'])
    eq1 = a*f(x,y) + b*f(x,y).diff(x) + c*f(x,y).diff(y)
    eq2 = 3*f(x,y) + 2*f(x,y).diff(x) + f(x,y).diff(y)
    eq3 = a*f(x,y) + b*f(x,y).diff(x) + 2*f(x,y).diff(y)
    for eq in [eq1, eq2, eq3]:
        assert checkpdesol(eq, pdsolve(eq))[0]
    eq4 = x*f(x,y) + f(x,y).diff(x) + 3*f(x,y).diff(y)
    eq5 = 2*f(x,y) + 1*f(x,y).diff(x) + 3*f(x,y).diff(y)
    eq6 = f(x,y) + 1*f(x,y).diff(x) + 3*f(x,y).diff(y)
    assert checkpdesol(eq4, [pdsolve(eq5), pdsolve(eq6)]) == [
        (False, (x - 2)*F(3*x - y)*exp(-x/S(5) - 3*y/S(5))),
         (False, (x - 1)*F(3*x - y)*exp(-x/S(10) - 3*y/S(10)))]
    for eq in [eq4, eq5, eq6]:
        assert checkpdesol(eq, pdsolve(eq))[0]
    sol = pdsolve(eq4)
    sol4 = Eq(sol.lhs - sol.rhs, 0)
    with raises(NotImplementedError):
        checkpdesol(eq4, sol4, solve_for_func=False)


def test_solvefun():
    f, F, G, H = map(Function, ['f', 'F', 'G', 'H'])
    eq1 = f(x,y) + f(x,y).diff(x) + f(x,y).diff(y)
    assert pdsolve(eq1) == Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))
    assert pdsolve(eq1, solvefun=G) == Eq(f(x, y), G(x - y)*exp(-x/2 - y/2))
    assert pdsolve(eq1, solvefun=H) == Eq(f(x, y), H(x - y)*exp(-x/2 - y/2))


def test_pde_1st_linear_constant_coeff_homogeneous():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = 2*u + u.diff(x) + u.diff(y)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(x - y)*exp(-x - y))
    assert checkpdesol(eq, sol)[0]

    eq = 4 + (3*u.diff(x)/u) + (2*u.diff(y)/u)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(2*x - 3*y)*exp(-S(12)*x/13 - S(8)*y/13))
    assert checkpdesol(eq, sol)[0]

    eq = u + (6*u.diff(x)) + (7*u.diff(y))
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(7*x - 6*y)*exp(-6*x/S(85) - 7*y/S(85)))
    assert checkpdesol(eq, sol)[0]

    eq = a*u + b*u.diff(x) + c*u.diff(y)
    sol = pdsolve(eq)
    assert checkpdesol(eq, sol)[0]


def test_pde_1st_linear_constant_coeff():
    f, F = map(Function, ['f', 'F'])
    u = f(x,y)
    eq = -2*u.diff(x) + 4*u.diff(y) + 5*u - exp(x + 3*y)
    sol = pdsolve(eq)
    # symox's simplify distributes exp(-y) (the Python reference body gave
    # the same form); the solutions are equal.
    expected = Eq(f(x,y), (F(4*x + 2*y)*exp(x/2) + exp(x + 4*y)/15)*exp(-y))
    assert sol == expected or (sol.rhs - expected.rhs).expand() == 0
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = (u.diff(x)/u) + (u.diff(y)/u) + 1 - (exp(x + y)/u)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y), F(x - y)*exp(-x/2 - y/2) + exp(x + y)/3)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = 2*u + -u.diff(x) + 3*u.diff(y) + sin(x)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y),
         F(3*x + y)*exp(x/5 - 3*y/5) - 2*sin(x)/5 - cos(x)/5)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = u + u.diff(x) + u.diff(y) + x*y
    sol = pdsolve(eq)
    assert sol.expand() == Eq(f(x, y),
        x + y + (x - y)**2/4 - (x + y)**2/4 + F(x - y)*exp(-x/2 - y/2) - 2).expand()
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]
    eq = u + u.diff(x) + u.diff(y) + log(x)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')


def test_pdsolve_all():
    f, F = map(Function, ['f', 'F'])
    u = f(x,y)
    eq = u + u.diff(x) + u.diff(y) + x**2*y
    sol = pdsolve(eq, hint = 'all')
    keys = ['1st_linear_constant_coeff',
        '1st_linear_constant_coeff_Integral', 'default', 'order']
    assert sorted(sol.keys()) == keys
    assert sol['order'] == 1
    assert sol['default'] == '1st_linear_constant_coeff'
    assert sol['1st_linear_constant_coeff'].expand() == Eq(f(x, y),
        -x**2*y + x**2 + 2*x*y - 4*x - 2*y + F(x - y)*exp(-x/2 - y/2) + 6).expand()


def test_pdsolve_variable_coeff():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = x*(u.diff(x)) - y*(u.diff(y)) + y**2*u - y**2
    sol = pdsolve(eq, hint="1st_linear_variable_coeff")
    assert sol == Eq(u, F(x*y)*exp(y**2/2) + 1)
    assert checkpdesol(eq, sol)[0]

    eq = x**2*u + x*u.diff(x) + x*y*u.diff(y)
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(y*exp(-x))*exp(-x**2/2))
    assert checkpdesol(eq, sol)[0]

    eq = y*x**2*u + y*u.diff(x) + u.diff(y)
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(-2*x + y**2)*exp(-x**3/3))
    assert checkpdesol(eq, sol)[0]

    eq = exp(x)**2*(u.diff(x)) + y
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, y*exp(-2*x)/2 + F(y))
    assert checkpdesol(eq, sol)[0]

    eq = exp(2*x)*(u.diff(y)) + y*u - u
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(x)*exp(-y*(y - 2)*exp(-2*x)/2))


# --- live parity with SymPy -------------------------------------------------

EQUATIONS = [
    'u.diff(x) + u.diff(y)',
    '2*u + 3*u.diff(x) + 4*u.diff(y)',
    'a*u + b*u.diff(x) + c*u.diff(y)',
    'u.diff(x) + u.diff(y) - x',
    '2*u.diff(x) - u.diff(y) + u - x*y',
    '-2*u.diff(x) + 4*u.diff(y) + 5*u - exp(x + 3*y)',
    'u + u.diff(x) + u.diff(y) + x**2*y',
    'u + u.diff(x) + u.diff(y) + log(x)',
    'Eq(u.diff(x) + u.diff(y), x)',
    'Eq(u.diff(x), u.diff(y))',
    'x*u.diff(x) - y*u.diff(y) + y**2*u - y**2',
    'x**2*u + x*u.diff(x) + x*y*u.diff(y)',
    'y*x**2*u + y*u.diff(x) + u.diff(y)',
    'exp(x)**2*u.diff(x) + y',
    'exp(2*x)*u.diff(y) + y*u - u',
    'x*u.diff(x) + u - y',
    'u.diff(x)**2 + u.diff(y)',
    'u.diff(x, 2) + u.diff(y, 2)',
    'u.diff(x)*u + u.diff(y)',
    'u',
]


def build(module, source):
    a, b, c, x, y = module.symbols('a b c x y')
    u = module.Function('u')(x, y)
    names = vars(module) | {'a': a, 'b': b, 'c': c, 'x': x, 'y': y, 'u': u}
    return eval(source, names), u


def outcome(call):
    try:
        return 'ok', call()
    except Exception as error:  # noqa: BLE001 - the exception is compared
        return 'err', type(error).__name__


def convert(value):
    if isinstance(value, dict):
        return [(convert(k), convert(v)) for k, v in value.items()]
    if isinstance(value, (list, tuple)):
        return type(value)(convert(v) for v in value)
    if isinstance(value, Exception):
        return type(value).__name__
    if hasattr(value, 'to_sympy'):
        return sp.srepr(value.to_sympy())
    if isinstance(value, sp.Basic):
        return sp.srepr(value)
    return value


@pytest.mark.parametrize('source', EQUATIONS)
def test_classify_pde_matches_sympy(source):
    eq, u = build(sx, source)
    ref, ru = build(sp, source)
    mine = outcome(lambda: classify_pde(eq))
    theirs = outcome(lambda: sp.classify_pde(ref))
    assert mine == theirs
    mine = outcome(lambda: classify_pde(eq, u, dict=True))
    theirs = outcome(lambda: sp.classify_pde(ref, ru, dict=True))
    assert mine[0] == theirs[0]
    if mine[0] == 'ok':
        assert list(mine[1]) == list(theirs[1])
        assert convert(mine[1]) == convert(theirs[1])


@pytest.mark.parametrize('source', EQUATIONS)
@pytest.mark.parametrize('hint', ['default', '1st_linear_constant_coeff_Integral',
                                  '1st_linear_variable_coeff', 'foo'])
def test_pdsolve_error_or_solution_is_checked_by_sympy(source, hint):
    eq, u = build(sx, source)
    ref, ru = build(sp, source)
    mine = outcome(lambda: pdsolve(eq, hint=hint))
    theirs = outcome(lambda: sp.pdsolve(ref, hint=hint))
    if theirs[0] == 'err':
        if mine[0] == 'ok':
            # symox extensions solve what SymPy cannot; the answer must hold
            check = sp.checkpdesol(ref, mine[1].to_sympy(), ru)
            assert check[0], (source, hint, check)
        else:
            assert mine == theirs
        return
    assert mine[0] == 'ok', (source, hint, mine)
    got = mine[1].to_sympy()
    if got != theirs[1]:
        assert sp.checkpdesol(ref, got.doit(), ru)[0], (source, hint, got, theirs[1])


# SymPy's Integral.doit raises on the one-sided limit (xi, x + y) of the
# unevaluated log integral ('<' between None and Zero); symox keeps the
# integral held. The difference is in doit, outside the pde lot.
KNOWN_ALL_DIFF = {'u + u.diff(x) + u.diff(y) + log(x)': 'TypeError'}


@pytest.mark.parametrize('source', EQUATIONS[:16])
def test_pdsolve_all_keys_match_sympy(source):
    eq, u = build(sx, source)
    ref, ru = build(sp, source)
    mine = outcome(lambda: pdsolve(eq, hint='all'))
    theirs = outcome(lambda: sp.pdsolve(ref, hint='all'))
    if source in KNOWN_ALL_DIFF:
        assert mine[0] == 'ok' and theirs == ('err', KNOWN_ALL_DIFF[source])
        return
    assert mine[0] == theirs[0]
    if mine[0] == 'ok':
        # SymPy iterates over a set of hint names, so its key order follows
        # the string hash seed.
        assert set(mine[1]) == set(theirs[1])
        assert mine[1]['order'] == theirs[1]['order']
        assert mine[1]['default'] == theirs[1]['default']


def test_checkpdesol_raises_like_sympy_on_an_unsolved_form():
    # symox's reference body returned None here (a dead raise after a
    # return); the native port raises as SymPy does.
    f = Function('f')
    eq = x*f(x,y) + f(x,y).diff(x) + 3*f(x,y).diff(y)
    sol = pdsolve(eq)
    with raises(NotImplementedError, match='Unable to test'):
        checkpdesol(eq, Eq(sol.lhs - sol.rhs, 0), solve_for_func=False)


ERRORS = [
    ('pdsolve(f(x, y), f(x, y))', KeyError, 'ordered_hints'),
    ('pdsolve(f(x, y))', ValueError, 'cannot be automatically detected'),
    ("pdsolve(f(x, y).diff(x) + f(x, y), hint='foo')", ValueError, 'Hint not recognized: foo'),
    ("pdsolve(f(x, y).diff(x) + f(x, y), hint='1st_linear_variable_coeff')",
     ValueError, 'does not match hint'),
    ('pdsolve(f(x, y).diff(x)**2 + f(x, y).diff(y)**2 + x)', NotImplementedError, 'psolve: Cannot solve'),
    ('classify_pde(f(x).diff(x))', IndexError, 'tuple index out of range'),
    ('classify_pde(f(x, y, a).diff(x), f(x, y, a))', NotImplementedError, 'two variables'),
]


@pytest.mark.parametrize('source,error,message', ERRORS)
def test_errors_match_sympy(source, error, message):
    for module in (sx, sp):
        names = vars(module) | {'x': module.Symbol('x'), 'y': module.Symbol('y'),
                                'a': module.Symbol('a'), 'f': module.Function('f')}
        with raises(error, match=message):
            eval(source, names)


# --- native dispatch ----------------------------------------------------------

def test_native_entry_points_are_registered():
    for name in ('pde_classify', 'pde_pdsolve', 'pde_solve_hint',
                 'pde_checkpdesol', 'pde_separate', 'pde_ode_order'):
        assert callable(getattr(_core, name))


def test_public_callables_do_not_reach_the_reference_body(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('reference body reached')
    monkeypatch.setattr(pde, '_desolve', forbidden)
    monkeypatch.setattr(pde, '_preprocess', forbidden)
    monkeypatch.setattr(pde, '_separate', forbidden)
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = u + u.diff(x) + u.diff(y) + x**2*y
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
                                '1st_linear_constant_coeff_Integral')
    sol = pdsolve(eq)
    assert checkpdesol(eq, sol) == (True, 0)
    assert set(pdsolve(eq, hint='all')) == {'order', 'default', '1st_linear_constant_coeff',
                                            '1st_linear_constant_coeff_Integral'}
    eq = x*u.diff(x) - y*u.diff(y) + y**2*u - y**2
    assert pdsolve(eq, hint='1st_linear_variable_coeff') == Eq(u, F(x*y)*exp(y**2/2) + 1)
    t = Symbol('t')
    X, T = map(Function, 'XT')
    w = f(x, t)
    assert pde_separate(w.diff(x, 2) - w.diff(t, 2), w, [X(x), T(t)]) == \
        [D(X(x), x, x)/X(x), D(T(t), t, t)/T(t)]


def test_lie_group_hooks_no_longer_call_python_pdsolve():
    from symox.solvers.ode.lie_group import _hooks
    hooks = _hooks()
    assert 'pdsolve' not in hooks
    assert {'doit', 'dsolve'} <= set(hooks)
