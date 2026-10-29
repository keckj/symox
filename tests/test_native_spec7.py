"""Scalar functions applied directly to matrix expressions (wave 21, spec7).

Expected values are hard-coded from SymPy 1.14: ``sin(MatrixSymbol('A', n, n))``
is an unevaluated, noncommutative application equal to ``sin(x).subs(x, A)``,
with the head's own eval rules for a matrix argument.
"""

import pytest

import symox as sx
from symox.matrices.exceptions import NonInvertibleMatrixError
from symox import (MatrixSymbol, Symbol, Identity, ZeroMatrix, Function, Eq, Rational,
                   sin, cos, tan, cot, sec, csc, asin, acos, atan, acot, sinh, cosh, tanh,
                   coth, asinh, acosh, atanh, exp, log, sqrt, cbrt, gamma, erf, erfc, sign,
                   arg, floor, ceiling, loggamma, latex, pretty, ccode)

n = Symbol('n')
x = Symbol('x')
A = MatrixSymbol('A', n, n)
B = MatrixSymbol('B', 2, 2)
C = MatrixSymbol('C', n, 1)

HEADS = [exp, log, sin, cos, tan, cot, sec, csc, asin, acos, atan, acot, sinh, cosh,
         tanh, coth, asinh, acosh, atanh, gamma, erf, erfc]


@pytest.mark.parametrize('head', HEADS, ids=lambda h: h.__name__)
@pytest.mark.parametrize('M', [A, B, C, A.T, A*A, x*A, A + A.T], ids=str)
def test_direct_application_matches_substitution(head, M):
    value = head(M)
    assert value == head(x).subs(x, M)
    assert hash(value) == hash(head(x).subs(x, M))
    assert value.func is head
    assert value.args == (M,)
    assert value.is_commutative is False
    assert str(value) == '%s(%s)' % (head.__name__, M)


def test_sin_matrixsymbol_object():
    value = sin(A)
    assert type(value).__name__ == 'sin'
    assert isinstance(value, sx.Function)
    assert not isinstance(value, sx.MatrixExpr)
    assert value.free_symbols == {A}
    assert value.func(B) == sin(B)
    assert value.diff(x) == 0
    assert value.subs(n, 2) == sin(MatrixSymbol('A', 2, 2))
    assert value.subs(n, 2) != sin(A)
    assert value.subs(A, B) == sin(B)
    assert sin(B).doit() == sin(B)
    assert sx.simplify(value) == value
    assert value.evalf() == value
    assert value.is_zero is False
    assert value.is_real is False
    with pytest.raises(AttributeError):
        value.T
    assert str(sin(sin(A))) == 'sin(sin(A))'
    assert str(sin(2*A)) == 'sin(2*A)'
    assert str(sin(B + B)) == 'sin(2*B)'
    assert str(sin(B.T)) == 'sin(B.T)'
    assert str(sin(Identity(2))) == 'sin(I)'
    assert str(sin(ZeroMatrix(2, 2))) == 'sin(0)'
    assert str(exp(ZeroMatrix(2, 2))) == 'exp(0)'
    assert str(exp(A).expand()) == 'exp(A)'
    assert str(sx.expand_trig(sin(2*A))) == 'sin(2*A)'


def test_printing():
    assert latex(sin(A)) == r'\sin{\left(A \right)}'
    assert latex(exp(A)) == 'e^{A}'
    assert pretty(cos(A)) == 'cos(A)'
    assert ccode(sin(A)) == 'sin(A)'
    assert str(cos(A).rewrite(exp)) == '(exp(-I*A) + exp(I*A))/2'


def test_noncommutative_arithmetic():
    assert str(2*sin(A)) == '2*sin(A)'
    assert str(sin(A) + 1) == '1 + sin(A)'
    assert str(sin(A)*sin(B)) == 'sin(A)*sin(B)'
    assert str(sin(B)*sin(A)) == 'sin(B)*sin(A)'
    assert sin(A)*sin(B) != sin(B)*sin(A)
    assert str(sx.expand((sin(A) + sin(B))**2)) == (
        'sin(A)*sin(B) + sin(A)**2 + sin(B)*sin(A) + sin(B)**2')
    assert Function('f')(A).is_commutative is False
    assert str(Eq(sin(A), sin(B))) == 'Eq(sin(A), sin(B))'


def test_noncommutative_scalar_in_matmul():
    msg = 'noncommutative scalars in MatMul are not supported.'
    with pytest.raises(NotImplementedError, match=msg):
        sin(A)*A
    with pytest.raises(NotImplementedError, match=msg):
        Function('f')(A)*A
    with pytest.raises(NotImplementedError, match=msg):
        Symbol('z', commutative=False)*A
    with pytest.raises(TypeError):
        sin(A) + A


@pytest.mark.parametrize('head', [sin, tan, cot, csc, asin, atan, acot, sinh, tanh, coth,
                                  asinh, atanh, erf, sign], ids=lambda h: h.__name__)
def test_odd_heads(head):
    assert head(-A) == -head(A)
    assert str(head(-A)) == '-%s(A)' % head.__name__
    if head is sign:
        assert str(head(-2*A)) == '-sign(A)'
    else:
        assert str(head(-2*A)) == '-%s(2*A)' % head.__name__


@pytest.mark.parametrize('head', [cos, sec, cosh], ids=lambda h: h.__name__)
def test_even_heads(head):
    assert head(-A) == head(A)
    assert str(head(-2*A)) == '%s(2*A)' % head.__name__


def test_sign_rules():
    assert str(erfc(-A)) == '2 - erfc(A)'
    assert str(erfc(-2*A)) == '2 - erfc(2*A)'
    for head in (acos, acosh, gamma, exp, log):
        assert str(head(-A)) == '%s(-A)' % head.__name__
    assert sign(2*A) == sign(A)
    assert str(sign(-2*A)) == '-sign(A)'
    assert arg(2*A) == arg(A)
    assert str(arg(-A)) == 'arg(-A)'
    assert str(arg(-2*A)) == 'arg(-A)'
    # sign and arg are real: their applications commute.
    assert sign(A).is_commutative is True
    assert arg(A).is_commutative is True
    # (-x)*A does not extract a sign (its first factor is not a number).
    assert str(sin(-x*A)) == 'sin((-x)*A)'
    # A - A.T = MatAdd(-A.T, A) does not extract a sign either.
    assert sin(A - A.T).args == (A - A.T,)
    assert sin(A - A.T) == sin(x).subs(x, A - A.T)


def test_exp_log_and_powers():
    value = exp(log(A))
    assert value == A
    assert type(value).__name__ == 'MatrixSymbol'
    assert str(log(exp(A))) == 'log(exp(A))'
    assert str(sin(log(A))) == 'sin(log(A))'
    assert sin(log(A)).is_commutative is False
    assert type(sqrt(A)).__name__ == 'MatPow'
    assert sqrt(A) == A**Rational(1, 2)
    assert str(sqrt(A)) == 'A**(1/2)'
    assert str(sqrt(-A)) == '(-A)**(1/2)'
    assert str(cbrt(B)) == 'B**(1/3)'
    assert sqrt(Identity(2)) == Identity(2)
    assert sqrt(ZeroMatrix(2, 2)) == ZeroMatrix(2, 2)
    assert cbrt(ZeroMatrix(2, 2)) == ZeroMatrix(2, 2)
    assert ZeroMatrix(n, n)**Rational(1, 3) == ZeroMatrix(n, n)
    with pytest.raises(NonInvertibleMatrixError, match='Matrix det == 0; not invertible'):
        ZeroMatrix(2, 2)**Rational(-1, 2)


def test_reference_errors():
    for head in (floor, ceiling):
        for M in (A, -A, B, Identity(2), ZeroMatrix(2, 2)):
            with pytest.raises(TypeError, match='Mix of Matrix and Scalar symbols'):
                head(M)
    with pytest.raises(NotImplementedError):
        loggamma(A)
    for head in (sec, csc):
        with pytest.raises(NonInvertibleMatrixError, match='Matrix det == 0; not invertible.'):
            head(ZeroMatrix(2, 2))
    assert str(sec(Identity(2))) == 'sec(I)'


def test_unevaluated():
    value = sin(-A, evaluate=False)
    assert str(value) == 'sin(-A)'
    assert value.is_commutative is False


# DiscreteMarkovChain with Str states (SymPy 1.14 values).

def _str_chain():
    from symox.core.symbol import Str
    from symox.stats import DiscreteMarkovChain
    T = sx.Matrix([[Rational(1, 2), Rational(1, 2)], [Rational(1, 3), Rational(2, 3)]])
    return Str, DiscreteMarkovChain('Y', [Str('a'), Str('b')], T), T


def test_str_atom():
    from symox.core.symbol import Str
    a = Str('a')
    assert isinstance(a, sx.Atom) and isinstance(a, sx.Basic)
    assert a == Str('a') and hash(a) == hash(Str('a')) and a != Str('b')
    assert a != Symbol('a')
    assert a.name == 'a' and a.args == () and a.is_Atom
    assert a.free_symbols == set()
    assert str(a) == 'a' and repr(a) == 'a'
    assert sx.srepr(a) == "Str('a')"
    assert str(sx.FiniteSet(a, Str('b'))) == '{a, b}'
    assert str(sx.Tuple(a, Str('b'))) == '(a, b)'
    assert str(Eq(Symbol('y'), a)) == 'Eq(y, a)'
    assert Eq(Symbol('y'), a).rhs == a
    assert type(sx.Tuple(a)[0]) is Str
    with pytest.raises(TypeError):
        Str(1)


def test_str_srepr_in_containers():
    from symox.core.symbol import Str
    assert sx.srepr(sx.Tuple(Str('a'), Str('b'))) == "Tuple(Str('a'), Str('b'))"
    assert sx.srepr(Eq(Symbol('y'), Str('a'))) == "Equality(Symbol('y'), Str('a'))"


def test_markov_chain_str_states():
    Str, Y, T = _str_chain()
    a, b = Str('a'), Str('b')
    from symox.stats import P, E
    assert str(Y.state_space) == '{a, b}'
    assert Y.number_of_states == 2
    assert Y.transition_probabilities == T
    assert P(Eq(Y[1], b), Eq(Y[0], a)) == Rational(1, 2)
    assert P(Eq(Y[2], a), Eq(Y[0], a)) == Rational(5, 12)
    assert P(sx.Ne(Y[1], a), Eq(Y[0], a)) == Rational(1, 2)
    assert P(Eq(Y[1], b), Eq(Y[0], a) & Eq(Y[0], a)) == Rational(1, 2)
    assert str(P(Eq(Y[3], b))) == 'Probability(Eq(Y[3], b))'
    assert str(Eq(Y[1], a)) == 'Eq(Y[1], a)'
    assert Y.stationary_distribution() == sx.Matrix([[Rational(2, 5), Rational(3, 5)]])
    assert Y.fixed_row_vector() == sx.Matrix([[Rational(2, 5), Rational(3, 5)]])
    assert Y.communication_classes() == [([a, b], True, 1)]
    assert Y.canonical_form() == ([a, b], T)
    assert Y.index_of == {a: 0, b: 1}
    states, A, B, C = Y.decompose()
    assert states == [a, b] and A == T and B.shape == (0, 2) and C.shape == (0, 0)
    with pytest.raises(TypeError):
        E(Y[2], Eq(Y[0], a))


def test_markov_chain_str_states_without_transition():
    from symox.core.symbol import Str
    from symox.stats import DiscreteMarkovChain
    Z = DiscreteMarkovChain('Z', [Str('x'), Str('y'), Str('z')])
    assert str(Z.state_space) == '{x, y, z}'
    assert str(Z.transition_probabilities) == '_T'
    assert Z.transition_probabilities.shape == (3, 3)


def test_user_atom_sympyrepr_inside_containers():
    class Label(sx.Atom):
        def _sympyrepr(self, printer):
            return 'Label<>'
    assert sx.srepr(sx.Tuple(Label(), 1)) == 'Tuple(Label<>, Integer(1))'
    assert sx.srepr(Label()) == 'Label<>'


def test_matrix_substituted_into_matmul_coefficient():
    y = Symbol('y')
    A2 = MatrixSymbol('A', 2, 2)
    value = (x*y).subs(x, A2).subs(y, B)
    assert type(value).__name__ == 'MatMul' and str(value) == 'B*A'
    assert str(sx.MatMul(y, A2).subs(y, B)) == 'B*A'
    assert str((x*y).subs({x: A2, y: B})) == 'B*A'
    with pytest.raises(NotImplementedError, match='noncommutative scalars in MatMul'):
        (x*y).subs(x, A2).subs(y, sin(B))
