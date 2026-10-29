"""Symbolic derivative orders (SymPy issue 22210), expected values from SymPy 1.14.

The closed forms are the reference ``_eval_derivative_n_times`` rules:
an atom gives a Piecewise in the order, a product the Leibniz Sum
(with ``piecewise_fold`` of the summand by the Sum constructor), and any
other body stays a held Derivative.
"""
import pytest
from symox import *
from symox import Dummy
del test

x, y, z = symbols('x y z')
n = Symbol('n', integer=True)
m = Symbol('m', integer=True)
mp = Symbol('mp', integer=True, positive=True)
f, g = Function('f'), Function('g')
nc = Symbol('z', commutative=False)

TERM = 'factorial(n)/(factorial(_k1)*factorial(-_k1 + n))'
SIN_TERM = 'factorial(n)*Derivative(sin(x), (x, Max(0, -_k1 + n)))/(factorial(_k1)*factorial(-_k1 + n))'

CASES = [
    (lambda: diff(x, (x, n)), 'Piecewise((x, Eq(n, 0)), (1, Eq(n, 1)), (0, True))'),
    (lambda: diff(x, (x, mp)), 'Piecewise((1, Eq(mp, 1)), (0, True))'),
    (lambda: diff(x, (x, n + 1)), 'Piecewise((x, Eq(n, -1)), (1, Eq(n, 0)), (0, True))'),
    (lambda: diff(x, (x, 2*n)), 'Piecewise((x, Eq(n, 0)), (0, True))'),
    (lambda: diff(x, (y, n)), 'Piecewise((x, Eq(n, 0)), (0, True))'),
    (lambda: diff(oo, (x, n)), 'Piecewise((oo, Eq(n, 0)), (0, True))'),
    (lambda: diff(nan, (x, n)), 'Piecewise((nan, Eq(n, 0)), (0, True))'),
    (lambda: diff(pi, (x, n)), 'Piecewise((pi, Eq(n, 0)), (0, True))'),
    (lambda: diff(I, (x, n)), 'Piecewise((I, Eq(n, 0)), (0, True))'),
    (lambda: diff(S.Half, (x, n)), 'Piecewise((1/2, Eq(n, 0)), (0, True))'),
    (lambda: diff(Float(2.5), (x, n)), 'Piecewise((2.5, Eq(n, 0)), (0, True))'),
    (lambda: diff(zoo, (x, n)), 'Piecewise((zoo, Eq(n, 0)), (0, True))'),
    (lambda: diff(0, (x, n)), '0'),
    (lambda: diff(nc, (x, n)), 'Piecewise((z, Eq(n, 0)), (0, True))'),
    (lambda: diff(x, (x, n), y), '0'),
    (lambda: diff(x*y, (x, n), (y, 2)), '0'),
    (lambda: diff(x*y, (y, 2), (x, n)), '0'),
    (lambda: diff(x*y, (f(x), n)), '0'),
    (lambda: diff(x, (x, n)).subs(n, 1), '1'),
    (lambda: diff(f(x), (f(x), n)), 'Piecewise((f(x), Eq(n, 0)), (1, Eq(n, 1)), (0, True))'),
    (lambda: diff(Derivative(x, (x, n)), x), 'Piecewise((x, Eq(n, -1)), (1, Eq(n, 0)), (0, True))'),
    (lambda: Derivative(x, (x, n)).doit(), 'Piecewise((x, Eq(n, 0)), (1, Eq(n, 1)), (0, True))'),
    (lambda: (2*Derivative(x, (x, n))).diff(x), '2*Piecewise((x, Eq(n, -1)), (1, Eq(n, 0)), (0, True))'),
    (lambda: diff(2*x, (x, n)),
     'Sum(Piecewise((2*x*%s, Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 0)), (2*%s, Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 1)), (0, True)), (_k1, 0, n))' % (TERM, TERM)),
    (lambda: diff(-x, (x, n)),
     'Sum(Piecewise((-x*%s, Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 0)), (-%s, Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 1)), (0, True)), (_k1, 0, n))' % (TERM, TERM)),
    (lambda: diff(x*sin(x), (x, n)),
     'Sum(Piecewise((x*%s, Eq(_k1, 0)), (%s, Eq(_k1, 1)), (0, True)), (_k1, 0, n))' % (SIN_TERM, SIN_TERM)),
    (lambda: Derivative(x*sin(x), (x, n)).doit(),
     'Sum(Piecewise((x*%s, Eq(_k1, 0)), (%s, Eq(_k1, 1)), (0, True)), (_k1, 0, n))' % (SIN_TERM, SIN_TERM)),
    (lambda: diff(3*sin(x), (x, n)), 'Sum(Piecewise((3*%s, Eq(_k1, 0)), (0, True)), (_k1, 0, n))' % SIN_TERM),
    (lambda: diff(f(x)*g(x), (x, n)),
     'Sum(factorial(n)*Derivative(f(x), (x, _k1))*Derivative(g(x), (x, Max(0, -_k1 + n)))/(factorial(_k1)*factorial(-_k1 + n)), (_k1, 0, n))'),
    (lambda: diff(exp(x)*sin(x), (x, n)),
     'Sum(factorial(n)*Derivative(exp(x), (x, _k1))*Derivative(sin(x), (x, Max(0, -_k1 + n)))/(factorial(_k1)*factorial(-_k1 + n)), (_k1, 0, n))'),
    (lambda: diff(f(x)*g(x)*sin(x), (x, n)),
     'Sum(factorial(n)*Derivative(f(x), (x, _k1))*Derivative(g(x), (x, _k2))*Derivative(sin(x), (x, Max(0, -_k1 - _k2 + n)))/(factorial(_k1)*factorial(_k2)*factorial(-_k1 - _k2 + n)), (_k1, 0, n), (_k2, 0, n))'),
    (lambda: diff(f(x)*x, (f(x), n)),
     'Sum(Piecewise((x*f(x)*%s, Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 0)), (x*%s, Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 1)), (0, True)), (_k1, 0, n))' % (TERM, TERM)),
    (lambda: diff(x*f(x), (x, 1), (x, n)),
     'Sum(Piecewise((x*factorial(n + 1)*Derivative(f(x), (x, Max(0, -_k1 + n + 1)))/(factorial(_k1)*factorial(-_k1 + n + 1)), Eq(_k1, 0)), (factorial(n + 1)*Derivative(f(x), (x, Max(0, -_k1 + n + 1)))/(factorial(_k1)*factorial(-_k1 + n + 1)), Eq(_k1, 1)), (0, True)), (_k1, 0, n + 1))'),
    (lambda: diff(x*nc, (x, n)),
     'Sum(Piecewise((x*factorial(n)*z/(factorial(_k1)*factorial(-_k1 + n)), Eq(_k1, 0) & Eq(Max(0, -_k1 + n), 0)), (0, Eq(_k1, 0)), (factorial(n)*z/(factorial(_k1)*factorial(-_k1 + n)), Eq(_k1, 1) & Eq(Max(0, -_k1 + n), 0)), (0, True)), (_k1, 0, n))'),
    # Bodies without a reference closed form stay held.
    (lambda: diff(sin(x), (x, n)), 'Derivative(sin(x), (x, n))'),
    (lambda: diff(x**2, (x, n)), 'Derivative(x**2, (x, n))'),
    (lambda: diff(x*sin(x) + x, (x, n)), 'Derivative(x*sin(x) + x, (x, n))'),
    (lambda: Derivative(f(x), (x, n)).doit(), 'Derivative(f(x), (x, n))'),
    (lambda: Derivative(f(x), (x, n)).diff(f(x)), 'Derivative(f(x), (x, n), f(x))'),
    (lambda: Derivative(Derivative(f(x), (x, m + 5)), (x, n + 3)), 'Derivative(f(x), (x, m + n + 8))'),
    (lambda: Derivative(sin(x), (x, n)).rewrite(exp), 'Derivative(-I*(exp(I*x) - exp(-I*x))/2, (x, n))'),
    (lambda: Derivative(f(x), (x, n))._eval_derivative(x), 'Derivative(f(x), (x, n + 1))'),
]


@pytest.mark.parametrize('index', range(len(CASES)))
def test_reference_values(index):
    build, expected = CASES[index]
    assert str(build()) == expected


def test_reference_dummy_eq_form():
    # SymPy test_diff_nth_derivative.
    expected = Sum(Piecewise((2*x*factorial(n)/(factorial(y)*factorial(-y + n)),
        Eq(y, 0) & Eq(Max(0, -y + n), 0)),
        (2*factorial(n)/(factorial(y)*factorial(-y + n)), Eq(y, 0) & Eq(Max(0, -y + n), 1)),
        (0, True)), (y, 0, n))
    assert diff(2*x, (x, n)).dummy_eq(expected)
    assert diff(x, (x, n)) == Piecewise((x, Eq(n, 0)), (1, Eq(n, 1)), (0, True))


def test_leibniz_sum_specializes_to_integer_orders():
    exprm = x*sin(x)
    mul_diff = diff(exprm, (x, n))
    assert isinstance(mul_diff, Sum)
    for i in range(5):
        assert mul_diff.subs(n, i).doit() == exprm.diff((x, i)).expand()
    exprm2 = 2*y*x*sin(x)*cos(x)*log(x)*exp(x)
    dex = exprm2.diff((x, n))
    assert isinstance(dex, Sum)
    for i in range(4):
        assert dex.subs(n, i).doit().expand() == exprm2.diff((x, i)).expand()
    assert str(Derivative(x*sin(x), (x, n)).subs(n, 3).doit()) == '-(x*cos(x) + 3*sin(x))'


def test_sum_doit_finishes_held_derivatives():
    k = Symbol('k', integer=True)
    assert Sum(Derivative(sin(x), (x, k)), (k, 0, 2)).doit() == cos(x)
    assert Sum(Piecewise((Derivative(sin(x), (x, k)), Eq(k, 1)), (0, True)), (k, 0, 2)).doit() == cos(x)


def test_order_could_be_zero():
    # SymPy test_order_could_be_zero.
    n0 = Symbol('n', integer=True, nonnegative=True)
    assert diff(y, (x, n0)) == Piecewise((y, Eq(n0, 0)), (0, True))
    assert diff(y, (x, n0 + 1)) is S.Zero
    assert diff(y, (x, mp)) is S.Zero


def test_reference_subs_protocol():
    expr = diff(sin(x), (x, n))
    assert expr.subs(sin(x), cos(-x)) == Derivative(cos(-x), (x, n))
    assert expr.subs(n, 1).doit() == cos(x)
    assert expr.subs(n, 2).doit() == -sin(x)
    d = Derivative(f(x), (x, n))
    # A Float order stays inexact.
    count = d.subs(n, 2.0).args[1][1]
    assert isinstance(count, Float) and count.is_Float and not count.is_Integer
    assert d.subs(n, 0) == f(x)
    # issue 22210
    d2 = 2*Derivative(sin(x), (x, Symbol('d', integer=True)))
    assert d2.simplify() == d2


def test_parity_facts_decide_equalities():
    assert Eq(2*n, 1) is S.false
    assert Eq(2*n + 1, 0) is S.false
    assert Eq(n, Rational(1, 3)) is S.false
    assert (n - Rational(1, 3)).is_zero is False
    assert (2*n - 1).is_zero is False
    assert str(Eq(3*n, 2)) == 'Eq(3*n, 2)'


MESSAGE = 'Cannot give expansion for symbolic count'


@pytest.mark.parametrize('call', [
    lambda d: d.as_finite_difference(),
    lambda d: d.as_finite_difference(wrt=x),
    lambda d: Derivative(f(x), (x, 2.0)).as_finite_difference(),
    lambda d: d.series(x, 0, 2),
    lambda d: d.as_leading_term(x),
    lambda d: d._eval_derivative(y),
    lambda d: d.variables,
])
def test_symbolic_count_errors(call):
    with pytest.raises(TypeError, match=MESSAGE):
        call(Derivative(f(x), (x, n)))


def test_eval_derivative_of_held_derivatives():
    assert hasattr(Derivative(f(x), (x, n)), '_eval_derivative')
    held = Derivative(f(x), x)
    assert held._eval_derivative(x) == Derivative(f(x), (x, 2))
    assert Derivative(f(x, y), x)._eval_derivative(y) == Derivative(f(x, y), x, y)
