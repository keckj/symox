"""Native 2nd_hypergeometric (2F1) steps, the upper-limit Integral.doit and
the solve recast of an undefined-function target, against SymPy 1.14 and
against the reference Python bodies."""
import pytest

import symox as S
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode import hypergeometric as hg

x, u, n, C1 = S.symbols('x u n C1')
f = S.Function('f')
g = S.Function('g')
fx = f(x)
d1 = fx.diff(x)
d2 = fx.diff(x, 2)
R = S.Rational
ENV = {'x': x, 'fx': fx, 'd1': d1, 'd2': d2, 'R': R}


def _eq(src):
    return eval(src, dict(ENV))


# SymPy 1.14: '2nd_hypergeometric' in classify_ode(eq, f(x)) (or the error it
# raises) and dsolve(eq, f(x), hint='2nd_hypergeometric').
HYPER = [
    ('x*(x-1)*d2 + (R(3,2) - 2*x)*d1 + 2*fx', True,
     'Eq(f(x), C1*x**(5/2)*hyper((1/2, 3/2), (7/2,), x) + C2*hyper((-2, -1), (-3/2,), x))'),
    ('x*(x-1)*d2 + (R(7,2)*x)*d1 + fx', True,
     'Eq(f(x), (C1*(1 - x)**(5/2)*hyper((1/2, 2), (7/2,), 1 - x) + C2*hyper((-2, -1/2), (-3/2,), 1 - x))/(x - 1)**(5/2))'),
    ('x*(x-1)*d2 + (3 + R(7,2)*x)*d1 + fx', True,
     'Eq(f(x), (C1*(1 - x)**(11/2)*hyper((1/2, 2), (13/2,), 1 - x) + C2*hyper((-5, -7/2), (-9/2,), 1 - x))/(x - 1)**(11/2))'),
    ('-x**R(5,7)*(-416*x**R(9,7)/9 - 2385*x**R(5,7)/49 + 298*x/3)*fx/(196*(-x**R(6,7) + x)**2*(x**R(6,7) + x)**2) + d2', True,
     'Eq(f(x), x**(45/98)*(C1*x**(4/49)*hyper((-1/2, 1/3), (9/7,), x**(2/7)) + C2*hyper((-11/14, 1/21), (5/7,), x**(2/7)))/(x**(2/7) - 1)**(19/84))'),
    ('(x*x - x)*d2 + (5*x - 1)*d1 + 4*fx', True,
     'Eq(f(x), (C1 + C2*(log(x) + 4/(x + 1)))*hyper((-1, -1), (1,), x)/(x - 1)**3)'),
    # Gap outside the ODE solver (native and reference bodies agree):
    # integrate(x/(x**2 - 1), x) splits into log(x - 1)/2 + log(x + 1)/2
    # where SymPy's ratint gives log(x**2 - 1)/2.
    pytest.param('(1 - x**2)*d2 - 2*x*d1 + 6*fx', True,
     'Eq(f(x), (C1*sqrt(x**2)*hyper((-1/2, 2), (3/2,), x**2) + C2*hyper((-1, 3/2), (1/2,), x**2))*(x**2)**(1/4)/sqrt(x))',
     marks=pytest.mark.xfail(strict=True, reason='integrate(x/(x**2 - 1)) form')),
    ('x*(x-2)*d2 + (3*x - 1)*d1 + fx/4', True,
     'Eq(f(x), (x - 2)**(sqrt(3)/2)*(C1*sqrt(x/(x - 2))*hyper((-sqrt(3)/2, 3/2 - sqrt(3)/2), (3/2,), x/(x - 2))'
     ' + C2*hyper((1 - sqrt(3)/2, -sqrt(3)/2 - 1/2), (1/2,), x/(x - 2)))/(x - 2))'),
    ('x*(1 - x)*d2 + (1 - 3*x)*d1 - fx', False, None),
    ('x*(x-1)*d2 + (3*x - 2)*d1 + fx', IndexError, None),
    ('(x**2 - 1)*(x - 3)*d2 + (x**2+1)*d1 + fx', False, None),
    ('d2 + d1/x + fx/x**2', False, None),
    ('x**2*d2 + x*d1 + (x**2 - 4)*fx', False, None),
]


@pytest.mark.parametrize('src,classified,solution', HYPER)
def test_2f1_matches_sympy(src, classified, solution):
    eq = _eq(src)
    if classified is IndexError:
        with pytest.raises(IndexError):
            S.classify_ode(eq, fx)
        return
    assert ('2nd_hypergeometric' in S.classify_ode(eq, fx)) is classified
    if solution is not None:
        assert str(S.dsolve(eq, fx, hint='2nd_hypergeometric')) == solution


# Without the xfail marks: the native and reference bodies agree on every row.
_PLAIN = [tuple(h.values) if hasattr(h, 'values') else h for h in HYPER]
_PLAIN = [h for h in _PLAIN if h[1] is not IndexError]


def _reference_match(A, B):
    d = hg.equivalence_hypergeometric(A, B, fx)
    if d and d['type'] == '2F1':
        return hg.match_2nd_2F1_hypergeometric(d['I0'], d['k'], d['sing_point'], fx)
    return None


@pytest.mark.parametrize('src,classified,solution', _PLAIN)
def test_2f1_native_is_decided_and_equals_reference(src, classified, solution, monkeypatch):
    eq = _eq(src)
    A, B = hg.match_2nd_hypergeometric(eq, fx)
    native = cn.hyper_2f1_match(A, B, fx)
    assert native is not cn._MISSING
    reference = _reference_match(A, B)
    if reference is None:
        assert native is None
        return
    assert {k: native[k] for k in 'a b c k mobius'.split()} == {k: reference[k] for k in 'a b c k mobius'.split()}
    match_object = dict(native, A=A, B=B)
    sol = cn.hyper_2f1_solution(eq, fx, match_object)
    assert sol is not cn._MISSING
    monkeypatch.setattr(cn, 'hyper_2f1_solution', lambda *args: cn._MISSING)
    assert sol == hg.get_sol_2F1_hypergeometric(eq, fx, dict(match_object))


def test_2f1_index_error_is_left_to_reference():
    # SymPy indexes sing_point[1] of a single singular point and raises.
    A, B = hg.match_2nd_hypergeometric(_eq('x*(x-1)*d2 + (3*x - 2)*d1 + fx'), fx)
    assert cn.hyper_2f1_match(A, B, fx) is cn._MISSING


def test_2f1_parameters_without_basis_return_none():
    # c and c - a - b integers, c != 1: the reference returns None.
    match_object = {'a': S.Integer(1), 'b': S.Integer(1), 'c': S.Integer(3), 'k': S.Integer(1),
                    'mobius': x, 'type': '2F1', 'A': (5*x - 3)/(x**2 - x), 'B': 1/(x**2 - x)}
    assert cn.hyper_2f1_solution(d2, fx, match_object) is None


# SymPy 1.14: Integral(..., (u, b)).doit() keeps the bound variable of an
# unevaluated antiderivative (ExprWithLimits._eval_interval).
UPPER = [
    (S.Integral(1/S.log(u**2), (u, fx)), 'Integral(1/log(u**2), (u, f(x)))'),
    (S.Integral(1/S.log(u**2) + u, (u, fx)), 'Integral((u*log(u**2) + 1)/log(u**2), (u, f(x)))'),
    (S.Integral(2/S.log(u**2), (u, g(x))), '2*Integral(1/log(u**2), (u, g(x)))'),
    (S.Integral(S.exp(u**2) + u, (u, fx)), 'f(x)**2/2 + sqrt(pi)*erfi(f(x))/2'),
]


@pytest.mark.parametrize('integral,expected', UPPER)
def test_upper_limit_doit(integral, expected):
    assert str(integral.doit()) == expected


# SymPy 1.14: solve recasts an undefined application to a Dummy.
SOLVE = [
    # the unrad branch of `_solve`, native (solve/radical.rs)
    (S.Eq(-2*fx/S.sqrt(fx**3), C1 - x), '[4/(C1**2 - 2*C1*x + x**2)]'),
    (S.sqrt(fx) - x, '[x**2]'),
    # the roots of the generator sqrt(f(x)), then its inversion
    (S.Eq(fx**R(3, 2), x), '[x**(2/3), (-x**(1/3)/2 - sqrt(3)*I*x**(1/3)/2)**2, (-x**(1/3)/2 + sqrt(3)*I*x**(1/3)/2)**2]'),
    (S.exp(fx) - x, '[log(x)]'),
]


@pytest.mark.parametrize('equation,expected', SOLVE)
def test_solve_undefined_application(equation, expected):
    assert str(S.solve(equation, fx)) == expected


def test_solve_undefined_application_as_symbol():
    # Native recast of f(x) beyond rational equations: the same values as
    # for a plain symbol (whose form differs from SymPy's, a separate gap).
    y = S.Symbol('y')
    for equation in (S.sqrt(fx) + fx - x, S.log(fx) + fx - x):
        values = S.solve(equation, fx)
        assert values and values == [v.subs(y, fx) for v in S.solve(equation.subs(fx, y), y)]


# SymPy 1.14 dsolve on the corpus rows fact_12, fact_13 and liouville_08.
DSOLVE = [
    (S.exp(d1) - fx**2, 'factorable', 'Eq(Integral(1/log(_y**2), (_y, f(x))), C1 + x)'),
    (d1**2 - fx**3, 'factorable', 'Eq(f(x), 4/(C1**2 - 2*C1*x + x**2))'),
    (x**2*d1 + (fx**2 + n*fx)*d1**2 + d2, 'Liouville',
     'Eq(C1 + C2*lowergamma(1/3, x**3/3) + Integral(exp(_y**3/3)*exp(_y**2*n/2), (_y, f(x))), 0)'),
]


@pytest.mark.parametrize('eq,hint,expected', DSOLVE)
def test_corpus_rows(eq, hint, expected):
    assert str(S.dsolve(eq, fx, hint=hint)) == expected
