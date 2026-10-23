"""Output forms of SymPy 1.14 for codegen, ntheory, lambdify, holonomic and
hyperexpand records of the b2 audit batch; expected strings come from SymPy."""
import re

import symox
from symox import (DeferredVector, Dummy, Matrix, Rational, Symbol, factorint,
                   fcode, hyper, lambdify, sequence, sqrt, symbols)
from symox.abc import n, x, y, z


def test_find_linear_recurrence_coefficients():
    assert sequence(3*5**n + 12).find_linear_recurrence(20, gfvar=x)[0] == [6, -5]
    assert sequence(n**2).find_linear_recurrence(10) == [3, -3, 1]
    assert str(sequence(symox.fibonacci(n)).find_linear_recurrence(10, gfvar=x)) == \
        '([1, 1], -x/(x**2 + x - 1))'


def test_fcode_size_and_sum_keep_noncommutative_order():
    from symox.codegen.ast import FunctionDefinition, Return, real
    from symox.codegen.fnodes import array, size, sum_
    a = Symbol('a', real=True)
    fd = FunctionDefinition(real, 'rms', [array(a, dim=[':'], intent='in')],
                            [Return((sum_(a**2)/size(a))**.5)])
    assert fcode(fd, source_format='free', standard=2003) == (
        'real*8 function rms(a)\n'
        'real*8, dimension(:), intent(in) :: a\n'
        'rms = sqrt(sum(a**2)*1d0/size(a))\n'
        'end function')
    assert fcode(size(a)*2, source_format='free') == '2*size(a)'
    assert fcode(sum_(a) + 1, source_format='free') == '1 + sum(a)'


def test_factorint_limit_keeps_unfactored_cofactors():
    a = 1407633717262338957430697921446883
    assert factorint(a, limit=10000) == {7: 1, 991: 1, 202916782076162456022877024859: 1}
    assert factorint(2**64 + 1, limit=1000) == {18446744073709551617: 1}
    assert factorint(10**20 + 39, limit=50) == {100000000000000000039: 1}
    assert factorint(1001*2**61 - 1, limit=100, use_rho=False, use_pm1=False) == \
        {3: 1, 769382950740969215317: 1}


def test_lambdastr_forms():
    from symox.utilities.lambdify import lambdastr
    assert lambdastr((x, (y, z)), x + y) == 'lambda _0,_1: (lambda x,y,z: (x + y))(_0,_1[0],_1[1])'
    assert lambdastr(x, Rational(1, 3)*x**Rational(-1, 2) + 1/x) == 'lambda x: (1/x + (1/3)/sqrt(x))'
    assert lambdastr([x, y], [sqrt(x), x**y, Matrix([[x, 1]])]) == \
        'lambda x,y: ([sqrt(x), x**y, ImmutableDenseMatrix([[x, 1]])])'
    assert re.fullmatch(r'lambda x: \(q_\d+ \+ x\)|lambda x: \(x \+ q_\d+\)', lambdastr(x, Dummy('q') + x))
    f = symox.Function('f')
    assert re.sub(r'_Dummy_\d+', 'D', lambdastr(f(x), f(x)**2)) == 'lambda D: (D**2)'


def test_lambdify_integer_results():
    from symox.utilities.lambdify import implemented_function
    f = implemented_function('f', lambda t: t + 1)
    assert repr(lambdify(x, f(x))(4)) == '5'
    assert repr(lambdify(x, f(x)**2 + 1)(3)) == '17'
    assert repr(lambdify(x, [f(x), x])(3)) == '[4, 3]'
    assert lambdify(x, f(x) + 1)(2**60) == 2**60 + 2
    assert repr(lambdify(x, f(x)/2)(3)) == '2.0'
    X = DeferredVector('X')
    assert repr(lambdify(X, (X[0] + 2, X[2] + X[1]))([1, 2, 3])) == '(3, 5)'
    assert repr(lambdify(X, (X[0] + 2, X[2]*X[1]))([1.5, 2, 3])) == '(3.5, 6)'
    assert repr(lambdify((x, y), [x + y, x*y])(2, 3)) == '[5, 6]'


def test_from_hyper_initial_conditions():
    from symox.holonomic import from_hyper
    assert str(from_hyper(hyper([], [Rational(3, 2)], x**2/4))) == \
        'HolonomicFunction((-x) + (2)*Dx + (x)*Dx**2, x, 1, [sinh(1), -sinh(1) + cosh(1)])'
    assert str(from_hyper(hyper([1, 1], [Rational(3, 2)], x**2))) == (
        'HolonomicFunction((4*x) + (5*x**2 - 2)*Dx + (x**3 - x)*Dx**2, x, 2, '
        '[-sqrt(3)*I*asin(2)/6, -1/6 + 7*sqrt(3)*I*asin(2)/36])')


def test_hyperexpand_buckets_reduce_numeric_part_only():
    from symox.simplify.hyperexpand import G_Function, Hyper_Function
    a, b = symbols('a b')
    h = Rational(1, 2)
    assert str(G_Function([1, 3, 2, 3*h], [1 + y, y, x, x + 2], [2], [y]).compute_buckets()) == \
        '({0: [3, 2, 1], 1/2: [3/2]}, {y: [y, y + 1], x: [x, x + 2]}, {0: [2]}, {y: [y]})'
    assert str(G_Function([-h, a - 1, a + h], [y - 3, y + 2, 3*y], [], [x, 7]).compute_buckets()) == \
        '({1/2: [-1/2], a: [a - 1], a + 1/2: [a + 1/2]}, {y: [y - 3, y + 2], 3*y: [3*y]}, {}, {x: [x], 0: [7]})'
    assert str(Hyper_Function([h, a + 1, a - h, 2], [b, b + h, 3]).build_invariants()) == \
        '(0, ((0, 1), (1/2, 1), (a, 1), (a + 1/2, 1)), ((0, 1), (b, 1), (b + 1/2, 1)))'
