"""Differential tests of the w23rust native kernels against SymPy 1.14:
the public ``piecewise_fold`` and the constructor fold of Sum/Product/
Integral (``_core._piecewise_fold``, ``_core._fold_piecewise_function``) and
the variable/order plan of symbolic-order matrix and array derivatives
(``_core._derivative_orders``, ``_core._derivative_order_plan``). Expected
values are SymPy 1.14 output, except the marked intentional corrections."""
import builtins

import pytest

import symox
from symox import *  # noqa: F401,F403
from symox import _core
from symox import _public_calculus

del test

x, y, k = symbols('x y k')
n = Symbol('n', integer=True)
m = Symbol('m')
p = Symbol('p', positive=True)
A, B = symbols('A B', commutative=False)
X = MatrixSymbol('X', 2, 2)
NAMES = dict(vars(symox), x=x, y=y, k=k, n=n, m=m, p=p, A=A, B=B, X=X)

# (call, str(result) or exception type, type name or exception message)
PIECEWISE_FOLD = [
    ('piecewise_fold(Piecewise((x, x > 0), (0, True))*n)', 'Piecewise((n*x, x > 0), (0, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x < 0), (1, True)) + Piecewise((x**2, x < 0), (2, True)))', 'Piecewise((x**2 + x, x < 0), (3, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((Piecewise((x, y > 0), (1, True)), x < 1), (0, True)))', 'Piecewise((x, (y > 0) & (x < 1)), (1, x < 1), (0, True))', 'Piecewise'),
    ('piecewise_fold(sin(Piecewise((x, x > 0), (0, True))))', 'Piecewise((sin(x), x > 0), (0, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x > 0), (0, True))*Piecewise((1, x > 0), (2, True)), evaluate=False)', 'Piecewise((x, x > 0), (0, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x > 0), (0, True))*Piecewise((1, x > 0), (2, True)), evaluate=None)', 'Piecewise((x, x > 0), (0, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x > 0), (y, x > 0), (1, True), evaluate=False)*2, evaluate=None)', 'Piecewise((2*x, x > 0), (2, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x > 0), (y, x > 0), (1, True), evaluate=False)*2, evaluate=False)', 'Piecewise((2*x, x > 0), (2, True))', 'Piecewise'),
    ('piecewise_fold(A*Piecewise((B, x > 0), (1, True))*Piecewise((A, x > 0), (2, True)))', 'Piecewise((A*B*A, x > 0), (2*A, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x > 0), (2, True))**2)', 'Piecewise((x**2, x > 0), (4, True))', 'Piecewise'),
    ('piecewise_fold(2**Piecewise((x, x > 0), (y, True)))', 'Piecewise((2**x, x > 0), (2**y, True))', 'Piecewise'),
    ('piecewise_fold(Max(Piecewise((x, x > 0), (0, True)), Piecewise((1, y > 0), (2, True))))', 'Piecewise((Max(1, x), (x > 0) & (y > 0)), (1, y > 0), (Max(2, x), x > 0), (2, True))', 'Piecewise'),
    ('piecewise_fold(Eq(Piecewise((x, x > 0), (0, True)), 1))', 'Piecewise((Eq(x, 1), x > 0), (False, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x < 0), (1, x < 1), (2, True)) + Piecewise((y, x < 0), (3, True)))', 'Piecewise((x + y, x < 0), (4, x < 1), (5, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x < 0), (1, x < 1), (2, True)) + Piecewise((y, x < 0), (3, x < 1), (4, True)) + x)', 'Piecewise((2*x + y, x < 0), (x + 4, x < 1), (x + 6, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x < 0), (1.5, True)) + Piecewise((y, y < 0), (3, True)) + Piecewise((2*y, y < 0), (5, True)) + Piecewise((1, x < 0), (2, True)))', 'Piecewise((x + 3*y + 1, (x < 0) & (y < 0)), (x + 9, x < 0), (3*y + 3.5, y < 0), (11.5, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, ITE(x > 0, y > 0, y < 0)), (1, True))*2)', 'Piecewise((2*x, ((x <= 0) | (y > 0)) & ((x > 0) | (y < 0))), (2, True))', 'Piecewise'),
    ('piecewise_fold(Sum(Piecewise((n, n > 2), (0, True))*x, (n, 1, k)))', 'Sum(Piecewise((n*x, n > 2), (0, True)), (n, 1, k))', 'Sum'),  # intentional correction (SymPy differs)
    ('piecewise_fold(Sum(Piecewise((n, n > 2), (0, True))*Piecewise((1, x > 0), (2, True)), (n, 1, k)), evaluate=None)', 'Sum(Piecewise((n, (n > 2) & (x > 0)), (2*n, n > 2), (0, True)), (n, 1, k))', 'Sum'),  # intentional correction (SymPy differs)
    ('piecewise_fold(x + y)', 'x + y', 'Add'),
    ('piecewise_fold(Piecewise((x, x > 0), (0, True)))', 'Piecewise((x, x > 0), (0, True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((x, x > 0), (0, True)), evaluate=None)', 'Piecewise((x, x > 0), (0, True))', 'Piecewise'),
    ('piecewise_fold(Abs(Piecewise((x, x > 0), (-x, True))) + Piecewise((y, y > 0), (0, True))*Piecewise((x, x > 1), (1, True)))', 'Piecewise((x*y + Abs(x), (x > 1) & (y > 0)), (Abs(x), x > 1), (y + Abs(x), y > 0), (Abs(x), True))', 'Piecewise'),
    ('piecewise_fold(Piecewise((Piecewise((x, y > 0), (1, True)), x < 1), (Piecewise((2, y > 0), (3, True)), True)) + 1)', 'Piecewise((x + 1, (y > 0) & (x < 1)), (2, x < 1), (3, y > 0), (4, True))', 'Piecewise'),
]
CONSTRUCTORS = [
    ('Sum(Piecewise((x, x > 0), (0, True))*n, (n, 1, 3))', 'Sum(n*Piecewise((x, x > 0), (0, True)), (n, 1, 3))', 'Sum'),
    ('Sum(Piecewise((n, n > 2), (0, True))*x, (n, 1, k))', 'Sum(Piecewise((n*x, n > 2), (0, True)), (n, 1, k))', 'Sum'),
    ('Integral(Piecewise((x, x < 1), (x**2, True))*y + 1, (x, 0, 2))', 'Integral(Piecewise((x*y + 1, x < 1), (x**2*y + 1, True)), (x, 0, 2))', 'Integral'),
    ('Product(Piecewise((n, n > 1), (1, True))*Piecewise((2, n < 3), (n, True)), (n, 1, 4))', 'Product(Piecewise((2*n, (n > 1) & (n < 3)), (n**2, n > 1), (2, n < 3), (n, True)), (n, 1, 4))', 'Product'),
    ('Integral(Piecewise((1, x > y), (0, True))*Piecewise((x, y > 0), (0, True)), (x, 0, 1))', 'Integral(Piecewise((x, (x > y) & (y > 0)), (0, True)), (x, 0, 1))', 'Integral'),
    ('Integral(Piecewise((x, y > 0), (0, True))*Piecewise((1, x > 0), (2, True)))', 'ValueError', 'specify dummy variables for Piecewise((1, x > 0), (2, True))*Piecewise((x, y > 0), (0, True))'),
    ('Integral(Piecewise((x, y > 0), (1, True))*Piecewise((2, x > 0), (3, True)), x)', 'Integral(Piecewise((2*x, (x > 0) & (y > 0)), (2, x > 0), (3*x, y > 0), (3, True)), x)', 'Integral'),
    ('Sum(x*Piecewise((n, x > 0), (1, True)), (n, 0, k))', 'Sum(Piecewise((n*x, x > 0), (x, True)), (n, 0, k))', 'Sum'),
    ('Integral(Piecewise((x, x > 0), (1, True))*Piecewise((y, y > 0), (2, True)), (x, 0, 1), (y, 0, 1))', 'Integral(Piecewise((x*y, (x > 0) & (y > 0)), (2*x, x > 0), (y, y > 0), (2, True)), (x, 0, 1), (y, 0, 1))', 'Integral'),
    ('Product(Piecewise((n, n > 1), (1, True))*Piecewise((2, x < 3), (x, True)), (n, 1, 4))', 'Product(Piecewise((n*Piecewise((2, x < 3), (x, True)), n > 1), (Piecewise((2, x < 3), (x, True)), True)), (n, 1, 4))', 'Product'),
    ('Sum(sin(Piecewise((n, n > 1), (x, True))), (n, 1, k))', 'Sum(Piecewise((sin(n), n > 1), (sin(x), True)), (n, 1, k))', 'Sum'),
]
DERIVATIVES = [
    ('X.diff((x, n))', 'Derivative(X, (x, n))', 'ArrayDerivative'),
    ('ImmutableMatrix([x, x**2]).diff((x, n), y)', 'Matrix([[0], [0]])', 'ImmutableDenseMatrix'),
    ('Array([x, x**2]).diff((x, n))', 'Derivative([x, x**2], (x, n))', 'ArrayDerivative'),
    ('Derivative(Array([x*y, x**2]), (x, 2), (y, m), x)', 'Derivative([x*y, x**2], (x, 2), (y, m), x)', 'Derivative'),
    ('Derivative(Array([x*y, x**2]), (x, 2), (y, m), x).doit()', 'Derivative([0, 0], (y, m))', 'Derivative'),  # intentional correction (SymPy differs)
    ('Array([x*y**2, x**2]).diff(y, (x, n), y)', 'Derivative([x*y**2, x**2], (x, n), (y, 2))', 'ArrayDerivative'),
    ('X.diff((x, -1))', 'ValueError', 'order of differentiation must be nonnegative'),
    ('Matrix([x]).diff((x, n))', 'AttributeError', "'ArrayDerivative' object has no attribute 'as_mutable'"),
    ('X.diff(x, (x, n))', 'Derivative(X, (x, n + 1))', 'ArrayDerivative'),
    ('ImmutableMatrix([x**3]).diff(x, 2, (x, n))', 'Derivative(Matrix([[x**3]]), (x, n + 2))', 'ArrayDerivative'),
    ('Array([x**3, y]).diff((x, n), (x, 2))', 'Derivative([x**3, y], (x, n + 2))', 'ArrayDerivative'),
    ('Array([x**3, y]).diff((y, 2), (x, n))', 'Derivative([x**3, y], (x, n), (y, 2))', 'ArrayDerivative'),
    ('Derivative(Array([x, x**2]), x).doit()', '[1, 2*x]', 'ImmutableDenseNDimArray'),  # intentional correction (SymPy differs)
    ('Derivative(ImmutableMatrix([x**2*y]), (x, 2), y)', 'Derivative(Matrix([[x**2*y]]), (x, 2), y)', 'Derivative'),
    ('Derivative(ImmutableMatrix([x**2*y]), (x, 2), y).doit()', 'Matrix([[2]])', 'ImmutableDenseMatrix'),
    ('Derivative(Array([x, x**2]), (x, n)).diff(x)', 'Derivative([x, x**2], (x, n + 1))', 'ArrayDerivative'),
    ('Derivative(Array([x, x**2]), (x, n)).diff(y)', '0', 'Zero'),
    ('Array([x]).diff((x, p))', 'Derivative([x], (x, p))', 'ArrayDerivative'),
    ('Array([x]).diff((y, p))', '[0]', 'ImmutableDenseNDimArray'),
    ('X.diff((x, n), (y, n))', 'Derivative(X, (x, n), (y, n))', 'ArrayDerivative'),
    ('X.diff((y, n), (x, 3))', '0', 'ZeroMatrix'),
]
EXCEPTIONS = ('ValueError', 'AttributeError', 'TypeError', 'IndexError', 'NotImplementedError')


def _check(call, expected, detail):
    if expected in EXCEPTIONS:
        with pytest.raises(getattr(builtins, expected)) as info:
            eval(call, NAMES)
        assert str(info.value) == detail
        return
    result = eval(call, NAMES)
    assert str(result) == expected
    assert type(result).__name__ == detail


@pytest.mark.parametrize('call,expected,detail', PIECEWISE_FOLD + CONSTRUCTORS + DERIVATIVES)
def test_against_sympy(call, expected, detail):
    _check(call, expected, detail)


@pytest.mark.parametrize('call,expected,detail', PIECEWISE_FOLD + CONSTRUCTORS)
def test_python_fallback_agrees(monkeypatch, call, expected, detail):
    # a replaced constructor sends every fold to the Python algorithm
    monkeypatch.setattr(_core, '_piecewise_fold', lambda *a: pytest.fail('native fold used'))
    monkeypatch.setattr(_core, '_fold_piecewise_function', lambda *a: pytest.fail('native fold used'))
    monkeypatch.setattr(_public_calculus, '_native_fold_admitted', lambda: False)
    _check(call, expected, detail)


def test_native_fold_is_used(monkeypatch):
    calls = []
    native = _core._piecewise_fold
    monkeypatch.setattr(_core, '_piecewise_fold', lambda *a: calls.append(a) or native(*a))
    assert str(piecewise_fold(Piecewise((x, x > 0), (0, True))*n)) == 'Piecewise((n*x, x > 0), (0, True))'
    assert calls
    calls.clear()
    native_function = _core._fold_piecewise_function
    monkeypatch.setattr(_core, '_fold_piecewise_function', lambda *a: calls.append(a) or native_function(*a))
    assert str(Sum(Piecewise((n, n > 2), (0, True))*x, (n, 1, k))) == 'Sum(Piecewise((n*x, n > 2), (0, True)), (n, 1, k))'
    assert calls


def test_admission():
    assert _public_calculus._native_fold_admitted()
    expr = Piecewise((x, x > 0), (0, True))*n
    plain = x + y
    assert piecewise_fold(plain) is plain
    assert _core._has_piecewise(expr) is True
    assert _core._has_piecewise(plain) is False
    assert _core._has_piecewise(Integral(MatrixSymbol('Y', 2, 2)[0, 0], x)) is None
    with evaluate(False):
        assert not _public_calculus._native_fold_admitted()
    original = symox.And
    try:
        symox.And = lambda *a: original(*a)
        assert not _public_calculus._native_fold_admitted()
    finally:
        symox.And = original
    assert _public_calculus._native_fold_admitted()
    # a one-argument non-Piecewise result of evaluate=None is left to Python
    assert _core._piecewise_fold(Piecewise((x, x > 0), (0, True)), None) is not None


def test_derivative_order_kernels():
    orders = _core._derivative_orders
    assert orders((x, 2, y)) == [(x, Integer(2)), (y, Integer(1))]
    assert orders(((x, n), x)) == [(x, n + 1)]
    assert orders(((x, n), (x, -n))) == []
    assert orders((x, 0)) == []
    assert orders((2,)) is None
    assert orders(((x, Rational(1, 2)),)) is None
    assert orders(((x, 1, 2),)) is None
    assert orders(('x',)) is NotImplemented
    assert orders(((x, 2.5),)) is NotImplemented
    with pytest.raises(ValueError, match='order of differentiation must be nonnegative'):
        orders(((x, -1),))
    plan = _core._derivative_order_plan
    assert plan([(y, Integer(1)), (x, n), (y, Integer(1))], {x, y}) == ([(x, n), (y, Integer(2))], 0)
    assert plan([(x, Integer(2)), (y, m)], {x, y}) == ([(x, Integer(2)), (y, m)], 1)
    assert plan([(x, Integer(2))], {x}) == ([(x, Integer(2))], 1)
    assert plan([(y, p)], {x}) is True
    assert plan([(y, n)], {x}) == ([(y, n)], 0)
    assert plan([(x, Integer(1))], ['x']) is NotImplemented
