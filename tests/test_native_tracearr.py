"""Symbolic-order matrix/array derivatives, held Derivative of explicit
matrices and arrays, constructor piecewise_fold. Expected values are
hard-coded from SymPy 1.14."""
import pytest

from symox import (Array, Derivative, Identity, ImmutableMatrix, ImmutableSparseMatrix, Integral,
                   Matrix, MatrixSymbol, MutableDenseNDimArray, MutableSparseNDimArray, Piecewise,
                   Product, SparseMatrix, Sum, Symbol, Trace, ZeroMatrix, diff, latex, sin, srepr,
                   symbols)

x, y, k, m = symbols('x y k m')
n = Symbol('n', integer=True)
X = MatrixSymbol('X', 2, 2)
A_ = Array([x, x**2])


def kind(value):
    return type(value).__name__


def test_matrix_symbol_symbolic_order_is_held():
    d = X.diff((x, n))
    assert str(d) == 'Derivative(X, (x, n))' and kind(d) == 'ArrayDerivative'
    assert d.shape == (2, 2)
    assert srepr(d).startswith('ArrayDerivative(MatrixSymbol(')
    assert srepr(d).endswith("Tuple(Symbol('x'), Symbol('n', integer=True)))")
    assert str(d.doit()) == 'Derivative(X, (x, n))' and kind(d.doit()) == 'ArrayDerivative'
    assert str(d.subs(n, 2)) == 'Derivative(X, (x, 2))' and kind(d.subs(n, 2)) == 'ArrayDerivative'
    assert d.subs(n, 2).doit() == ZeroMatrix(2, 2)
    assert str(d.xreplace({n: 3})) == 'Derivative(X, (x, 3))'
    assert sorted(map(str, d.free_symbols)) == ['X', 'n']
    assert d.expr == X and str(d.variable_count) == '((x, n),)'
    assert d == X.diff((x, n)) and d.is_commutative is False
    assert latex(d) == r'\frac{d^{n}}{d x^{n}} X'
    assert str(diff(X, (x, n))) == 'Derivative(X, (x, n))'
    assert str(X.diff((x, n)).diff(x)) == '0'


@pytest.mark.parametrize('make,expected', [
    (lambda: (x*X).diff((x, n)), 'Derivative(x*X, (x, n))'),
    (lambda: (x*X).diff(x, (x, n)), 'Derivative(x*X, (x, n + 1))'),
    (lambda: (x*X).diff(x, (y, n)), 'Derivative(X, (y, n))'),
    (lambda: (x*y*X).diff((x, n), y), 'Derivative((x*y)*X, (x, n), y)'),
    (lambda: (x*X).diff((x, m)), 'Derivative(x*X, (x, m))'),
    (lambda: (X + x*Identity(2)).diff((x, n)), 'Derivative(x*I + X, (x, n))'),
    (lambda: X.diff((x, n), (x, 2)), 'Derivative(X, (x, n + 2))'),
    (lambda: X.diff((x, n + 1)), 'Derivative(X, (x, n + 1))'),
])
def test_matrix_expression_symbolic_orders(make, expected):
    d = make()
    assert str(d) == expected and kind(d) == 'ArrayDerivative' and d.shape == (2, 2)


def test_matrix_expression_zero_exit():
    assert X.diff((x, n), y) == ZeroMatrix(2, 2)
    assert X.diff(y, (x, n)) == ZeroMatrix(2, 2)
    with pytest.raises(ValueError, match='order of differentiation must be nonnegative'):
        X.diff((x, n), (y, -1))


@pytest.mark.parametrize('make', [
    lambda: X.diff((x, -1)), lambda: Matrix([x]).diff((x, -1)),
    lambda: ImmutableMatrix([x]).diff((x, -1)), lambda: Array([x]).diff((x, -1))])
def test_negative_integer_order_message(make):
    with pytest.raises(ValueError, match='^order of differentiation must be nonnegative$'):
        make()


def test_negated_symbol_order_is_held_and_integer_orders_evaluate():
    p = Symbol('p')
    d = X.diff((x, -p))
    assert str(d) == 'Derivative(X, (x, -p))' and kind(d) == 'ArrayDerivative'
    assert X.diff((x, 2)) == ZeroMatrix(2, 2) and X.diff((x, 0)) == X


def test_explicit_matrix_symbolic_orders():
    d = ImmutableMatrix([x]).diff((x, n))
    assert str(d) == 'Derivative(Matrix([[x]]), (x, n))' and kind(d) == 'ArrayDerivative'
    assert d.shape == (1, 1)
    assert srepr(d) == ("ArrayDerivative(ImmutableDenseMatrix([[Symbol('x')]]), "
                        "Tuple(Symbol('x'), Symbol('n', integer=True)))")
    assert str(d.doit()) == 'Derivative(Matrix([[x]]), (x, n))'
    assert str(d.subs(n, 2)) == 'Derivative(Matrix([[x]]), (x, 2))'
    assert d.subs(n, 2).doit() == ImmutableMatrix([[0]]) and kind(d.subs(n, 2).doit()) == 'ImmutableDenseMatrix'
    held = ImmutableMatrix([[x, y]]).diff((x, n), y)
    assert str(held) == 'Derivative(Matrix([[x, y]]), (x, n), y)' and held.shape == (1, 2)
    assert str(ImmutableMatrix([[x, y]]).diff(y, (x, n))) == 'Derivative(Matrix([[x, y]]), (x, n), y)'
    assert str(ImmutableSparseMatrix([x]).diff((x, n))) == 'Derivative(Matrix([[x]]), (x, n))'
    message = "'ArrayDerivative' object has no attribute 'as_mutable'"
    for make in (lambda: Matrix([x]).diff((x, n)), lambda: diff(Matrix([x]), (x, n)),
                 lambda: Matrix([[x, y]]).diff((x, n), y), lambda: SparseMatrix([x]).diff((x, n))):
        with pytest.raises(AttributeError) as error:
            make()
        assert str(error.value) == message


def test_array_symbolic_orders():
    for make in (lambda: A_.diff((x, n)), lambda: MutableDenseNDimArray([x, x**2]).diff((x, n)),
                 lambda: diff(A_, (x, n)), lambda: MutableSparseNDimArray([x, x**2]).diff((x, n))):
        d = make()
        assert str(d) == 'Derivative([x, x**2], (x, n))' and kind(d) == 'ArrayDerivative'
        assert d.shape == (2,)
    d = A_.diff((x, n))
    assert srepr(d) == ("ArrayDerivative(ImmutableDenseNDimArray(Tuple(Symbol('x'), Pow(Symbol('x'), "
                        "Integer(2))), Tuple(Integer(2))), Tuple(Symbol('x'), Symbol('n', integer=True)))")
    assert str(d.doit()) == 'Derivative([x, x**2], (x, n))'
    assert str(d.subs(n, 1)) == 'Derivative([x, x**2], x)' and kind(d.subs(n, 1)) == 'ArrayDerivative'
    assert str(d.subs(n, 1).doit()) == '[1, 2*x]' and kind(d.subs(n, 1).doit()) == 'ImmutableDenseNDimArray'
    assert str(A_.diff((x, n), y)) == '[0, 0]' and str(A_.diff(y, (x, n))) == '[0, 0]'
    assert str(Array([[x, y]]).diff((x, n))) == 'Derivative([[x, y]], (x, n))'
    assert str(Array([x]).diff((x, n)).diff(x)) == 'Derivative([x], (x, n + 1))'
    assert latex(d) == r'\frac{d^{n}}{d x^{n}} \left[\begin{matrix}x & x^{2}\end{matrix}\right]'


def test_derivative_of_explicit_matrix_is_held():
    d = Derivative(Matrix([x]), (x, n))
    assert str(d) == 'Derivative(Matrix([[x]]), (x, n))' and kind(d) == 'Derivative'
    assert srepr(d) == ("Derivative(ImmutableDenseMatrix([[Symbol('x')]]), "
                        "Tuple(Symbol('x'), Symbol('n', integer=True)))")
    assert str(d.doit()) == 'Derivative(Matrix([[x]]), (x, n))' and kind(d.doit()) == 'Derivative'
    assert not hasattr(d, 'shape')
    assert str(d.subs(n, 2)) == 'Derivative(Matrix([[x]]), (x, 2))' and kind(d.subs(n, 2)) == 'ArrayDerivative'
    assert d.free_symbols == {n, x} and d.expr == ImmutableMatrix([x]) and kind(d.expr) == 'ImmutableDenseMatrix'
    assert d == Derivative(ImmutableMatrix([x]), (x, n))
    one = Derivative(Matrix([x]), x)
    assert str(one) == 'Derivative(Matrix([[x]]), x)' and kind(one) == 'Derivative'
    assert srepr(one) == "Derivative(ImmutableDenseMatrix([[Symbol('x')]]), Tuple(Symbol('x'), Integer(1)))"
    assert one.doit() == ImmutableMatrix([[1]]) and kind(one.doit()) == 'ImmutableDenseMatrix'
    assert Derivative(Matrix([[x, x**2]]), x, 2).doit() == ImmutableMatrix([[0, 2]])
    assert str(Derivative(ImmutableMatrix([x]), (x, n), evaluate=True)) == 'Derivative(Matrix([[x]]), (x, n))'
    assert Derivative(ImmutableMatrix([x]), x, evaluate=True) == ImmutableMatrix([[1]])
    assert Derivative(Matrix([x]), x, evaluate=True) == ImmutableMatrix([[1]])
    assert str(Derivative(Matrix([x]), y).doit()) == '0' and kind(Derivative(Matrix([x]), y).doit()) == 'Zero'
    assert str(latex(d)) == r'\frac{d^{n}}{d x^{n}} \left[\begin{matrix}x\end{matrix}\right]'


def test_derivative_of_array_and_matrix_symbol_is_held():
    d = Derivative(A_, (x, n))
    assert str(d) == 'Derivative([x, x**2], (x, n))' and kind(d) == 'Derivative'
    assert srepr(d) == ("Derivative(ImmutableDenseNDimArray(Tuple(Symbol('x'), Pow(Symbol('x'), "
                        "Integer(2))), Tuple(Integer(2))), Tuple(Symbol('x'), Symbol('n', integer=True)))")
    assert str(d.doit()) == 'Derivative([x, x**2], (x, n))'
    assert str(Derivative(A_, x)) == 'Derivative([x, x**2], x)' and not hasattr(Derivative(A_, x), 'shape')
    assert str(Derivative(A_, y).doit()) == '0'
    assert str(Derivative(Array([x]), x).doit()) == '[1]'
    assert srepr(Derivative(MutableDenseNDimArray([x]), x)) == (
        "Derivative(ImmutableDenseNDimArray(Tuple(Symbol('x')), Tuple(Integer(1))), Tuple(Symbol('x'), Integer(1)))")
    assert str(Derivative(Array([x]), (x, n)).diff(x)) == 'Derivative([x], (x, n + 1))'
    held = Derivative(X, (x, n))
    assert str(held) == 'Derivative(X, (x, n))' and kind(held) == 'Derivative'
    assert str(held.doit()) == 'Derivative(X, (x, n))' and kind(held.doit()) == 'Derivative'
    assert str(Derivative(X, (x, n), y).doit()) == '0' and kind(Derivative(X, (x, n), y).doit()) == 'Zero'


def test_intentional_correction_derivative_of_array_doit():
    # SymPy 1.14 raises "Can't calculate derivative wrt x**2." here (the
    # array entries are taken for variables); the entrywise derivative is
    # [1, 2*x], as Array([x]) gives [1] in SymPy.
    assert str(Derivative(A_, x).doit()) == '[1, 2*x]'
    assert str(Derivative(A_, x, evaluate=True)) == '[1, 2*x]'


@pytest.mark.parametrize('make,expected,function', [
    (lambda: Sum(x*Piecewise((n, x > 0), (1, True)), (n, 0, k)),
     'Sum(Piecewise((n*x, x > 0), (x, True)), (n, 0, k))', 'Piecewise((n*x, x > 0), (x, True))'),
    (lambda: Sum(x*Piecewise((x, k > 0), (1, True)), (n, 0, k)),
     'Sum(x*Piecewise((x, k > 0), (1, True)), (n, 0, k))', 'x*Piecewise((x, k > 0), (1, True))'),
    (lambda: Sum(Piecewise((n, x > 0), (1, True)) + Piecewise((n**2, y > 0), (2, True)), (n, 0, k)),
     'Sum(Piecewise((n**2 + n, (x > 0) & (y > 0)), (n + 2, x > 0), (n**2 + 1, y > 0), (3, True)), (n, 0, k))',
     'Piecewise((n**2 + n, (x > 0) & (y > 0)), (n + 2, x > 0), (n**2 + 1, y > 0), (3, True))'),
    (lambda: Sum(sin(Piecewise((n, n > 2), (1, True))), (n, 0, k)),
     'Sum(Piecewise((sin(n), n > 2), (sin(1), True)), (n, 0, k))', 'Piecewise((sin(n), n > 2), (sin(1), True))'),
    (lambda: Product(x*Piecewise((n, x > 0), (1, True)), (n, 1, k)),
     'Product(Piecewise((n*x, x > 0), (x, True)), (n, 1, k))', 'Piecewise((n*x, x > 0), (x, True))'),
    (lambda: Integral(y*Piecewise((x, x > 0), (1, True)), x),
     'Integral(Piecewise((x*y, x > 0), (y, True)), x)', 'Piecewise((x*y, x > 0), (y, True))'),
    (lambda: Integral(y*Piecewise((x, x > 0), (1, True)) + Piecewise((1, y > 0), (2, True)), (x, 0, 1)),
     'Integral(Piecewise((x*y + Piecewise((1, y > 0), (2, True)), x > 0), '
     '(y + Piecewise((1, y > 0), (2, True)), True)), (x, 0, 1))',
     'Piecewise((x*y + Piecewise((1, y > 0), (2, True)), x > 0), (y + Piecewise((1, y > 0), (2, True)), True))'),
    (lambda: Integral(Piecewise((x, x > 0), (1, True))**2, x),
     'Integral(Piecewise((x**2, x > 0), (1, True)), x)', 'Piecewise((x**2, x > 0), (1, True))'),
    (lambda: Integral(Piecewise((x, k > 0), (1, True)), x),
     'Integral(Piecewise((x, k > 0), (1, True)), x)', 'Piecewise((x, k > 0), (1, True))'),
    (lambda: Sum(x*Piecewise((n, x > 0), (1, True)), (n, 0, k), (x, 0, 2)),
     'Sum(Piecewise((n*x, x > 0), (x, True)), (n, 0, k), (x, 0, 2))', 'Piecewise((n*x, x > 0), (x, True))'),
])
def test_user_built_limits_fold_piecewise(make, expected, function):
    value = make()
    assert str(value) == expected and str(value.function) == function


def test_nested_trace_derivatives():
    A, B = MatrixSymbol('A', 2, 2), MatrixSymbol('B', 2, 2)
    assert str(Trace(A*Trace(B*X)*X).diff(X)) == 'Trace(A*X)*B.T + Trace(B*X)*A.T'
    assert str(Trace(X**3).diff(X)) == '3*X.T**2'
    assert str((Trace(X)*Trace(X)).diff(X)) == '(2*Trace(X))*I'
