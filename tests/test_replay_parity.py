"""Regressions found by the wave-8 replay of the SymPy 1.14 docstring audit.

Expected values are hard-coded from SymPy 1.14.
"""
import time

from symox import (FunctionMatrix, Integer, KroneckerDelta, Lambda, MatrixSymbol, Number, Q, Rational, Symbol,
                   ask, cbrt, factorial, sqrt, symbols, sympify)
from symox.assumptions import BinaryRelation
from symox.physics.wigner import gaunt

x = Symbol('x')
i, j = symbols('i j')


def test_square_of_lazy_matrix_entry_is_one_row_times_one_column():
    A = MatrixSymbol('A', 2, 2)
    assert str((A**2)[0, 1]) == 'A[0, 0]*A[0, 1] + A[0, 1]*A[1, 1]'
    Z = FunctionMatrix(3, 3, Lambda((i, j), i*x + j))
    assert (Z**2)[1, 2] == sympify('3*x**2 + 11*x + 6')
    assert (Z**2)[0, 0] == 5*x
    assert (Z**2)[2, 1] == sympify('6*x**2 + 11*x + 3')
    D = FunctionMatrix(3, 3, Lambda((i, j), KroneckerDelta(i, j)*(i + x)))
    assert (D**2)[1, 1] == (x + 1)**2
    assert (D**2)[0, 1] == 0


def test_square_of_large_function_matrix_entry_is_fast():
    Y = FunctionMatrix(1000, 1000, Lambda((i, j), i + j))
    start = time.perf_counter()
    assert (Y**2)[10, 10] == 342923500
    assert time.perf_counter() - start < 10


def test_binary_relation_subclass_prints_and_asks():
    class MyRel(BinaryRelation):
        name = 'R'
        is_reflexive = False

    Q.R = MyRel()
    try:
        @Q.R.register(Number, Number)
        def _(n1, n2, assumptions):
            return ask(Q.zero(n1 - n2 - 1), assumptions)

        assert str(Q.R(2, 1)) == 'Q.R(2, 1)'
        assert repr(Q.R(2, 1)) == 'Q.R(2, 1)'
        assert ask(Q.R(2, 1)) is True
        assert ask(Q.R(1, 1)) is False
    finally:
        del Q.R


def test_gaunt_large_arguments_value_and_speed():
    start = time.perf_counter()
    value = gaunt(1000, 1000, 1200, 9, 3, -12)
    assert str(value.n(64)) == '0.006895004219221134484332976156744208248842039317638217822322799675'
    assert time.perf_counter() - start < 30
    assert str(gaunt(1, 0, 1, 1, 0, -1)) == '-1/(2*sqrt(pi))'
    assert str(gaunt(2, 2, 2, 0, 0, 0)) == 'sqrt(5)/(7*sqrt(pi))'


def test_integer_powers_keep_perfect_power_extraction():
    assert str(sqrt(Rational(10**41, 3))) == '100000000000000000000*sqrt(30)/3'
    assert str(cbrt(72)) == '2*3**(2/3)'
    assert str(Integer(1922)**Rational(2, 3)) == '31*124**(1/3)'
    assert str(sqrt(Integer(2**200*3**4*7))) == '11408855402054064613470328848384*sqrt(7)'
    assert Integer(3**60*5**30)**Rational(1, 30) == 45
    assert str(sqrt(factorial(30))) == '15692092416000*sqrt(1077205)'
    # Large perfect powers with and without small prime factors.
    assert sqrt(Integer(2**130*3**2)) == Integer(2**65*3)
    assert (Integer((2**89 - 1)**3))**Rational(1, 3) == 2**89 - 1
    assert (Integer(7**40*11**20))**Rational(1, 20) == 539
    assert str(sqrt(Integer(2**71*9))) == '103079215104*sqrt(2)'
    assert str(Integer(2**71*3**5)**Rational(2, 5)) == '2415919104*2**(2/5)'
    assert Integer(6**100)**Rational(1, 4) == Integer(28430288029929701376)


def test_expand_keeps_held_reciprocal_of_unchanged_sum():
    from symox import Float, I, Pow, expand
    from symox.physics.control.lti import TransferFunction
    s, a, p = symbols('s a p')
    held = Pow(Float(0.99) + Float(0.2)*I, -1, evaluate=False)
    assert str(held.expand()) == '1/(0.99 + 0.2*I)'
    assert str(TransferFunction(1, s**2 + 2*s + 1, s).eval_frequency(I*0.1)) == '1/(0.99 + 0.2*I)'
    assert str(TransferFunction(s**2, a*s + p, s).eval_frequency(I*2)) == '-4/(2*I*a + p)'
    assert str(expand(1/(1 + 2*I))) == '1/5 - 2*I/5'
    assert str(expand(1/(x*(x + 1)))) == '1/(x**2 + x)'
