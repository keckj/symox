"""Derived relations, Boolean atom operands, kinds and Float srepr.

Every case runs on SymPy and Symox and compares the result, or the
exception type and message.
"""
import warnings

import pytest
import sympy
import symox

MODULES = {'sympy': sympy, 'symox': symox}


def outcome(expression, name):
    m = MODULES[name]
    x, y = m.symbols('x y')
    namespace = {'m': m, 'x': x, 'y': y, 't': m.true, 'f': m.false, 'e': m.Eq(x, 1), 'a': m.And(x > 0, x < 1)}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            result = eval(expression, namespace)
        except Exception as error:  # noqa: BLE001 - the error is compared
            return ('error', type(error).__name__, str(error))
    return ('value', str(result))


RELATIONS = ['m.Lt(x, y)', 'm.Gt(x, 1)', 'm.Le(1, x)', 'm.Ge(-x, y)', 'm.Eq(-x, 1)', 'm.Ne(2, x)',
             'm.Gt(-x, -y)', 'm.Lt(-y, x)', 'm.Gt(-y, x)', 'm.Lt(-y, -x)', 'm.Eq(x, -y)', 'm.Eq(-x, y)',
             'm.Lt(x, 2, evaluate=False).reversed', 'm.Eq(x, t)', 'm.Eq(x + 1, y + 2)']
ATTRIBUTES = ['canonical', 'reversed', 'negated', 'reversedsign', 'lhs', 'rhs', 'rel_op', 'lts', 'gts',
              'strict', 'weak', 'func.__name__']

DERIVED = [relation + '.' + attribute for relation in RELATIONS for attribute in ATTRIBUTES]

BOOLEAN_ATOMS = [atom + '.' + attribute for atom in ('t', 'f')
                 for attribute in ('negated', 'canonical', 'as_set()', '_op_priority', 'is_Boolean', 'kind')]

OPERANDS = [
    't+1', '1+t', 'x+t', 't+x', 'x-t', 't-1', '1-t', 't*x', 'x*t', 'x**t', 't**2', '2**t', 't/2', '2/t',
    'x/t', 't%2', '2%t', 'x%t', 't//2', 'x//t', '-t', '+t', 'abs(t)', '~t', 'divmod(x, t)', 'divmod(t, 2)',
    'x+e', 'e+1', '1+e', 'x*e', '-e', 'x+a', 'a*2', '-a',
    'x<t', 't<x', '1<t', 't<1', 'x<e', 'e<1', 'x<a', 'a<1', 'x>=f', 't==1', 't!=1',
    'x+m.Tuple(1)', 'x<m.Tuple(1)', 'x<m.Interval(0, 1)', 'm.Interval(0, 1)<x',
    'x+m.Q.positive(x)', 'm.Contains(x, m.Interval(0, 1))+1',
]

INEQUALITIES = [
    'x<m.I', 'm.I<x', 'x<m.nan', 'x<m.zoo', 'x>=m.I+1', 'x>=2*m.I', 'x<m.Symbol("r", real=True)+m.I',
    'm.Lt(x, t)', 'm.Lt(t, 1)', 'm.Le(x, e)', 'm.Lt(x, m.I)', 'm.Gt(x, m.nan)', 'm.Lt(x, t, evaluate=False)',
    'm.Lt(t, m.I)', 'm.Ne(t, 1)', 'm.Lt(x, m.Tuple(1))', 'm.Lt(m.Interval(0, 1), 1)', 'm.Min(x, 2) < 1',
    'm.Lt(m.Symbol("y", real=False), 1)', 'm.Ge(m.Symbol("y", real=False), 1)',
    'm.Gt(m.Symbol("y", real=False), m.Symbol("z", real=False))', 'm.Lt(m.Q.positive(x), 1)',
    'm.Lt(m.Symbol("p", positive=True), 0)', 'm.Ge(m.Symbol("p", positive=True), 0)',
    'm.Float(0.36) <= m.exp(16777216)', 'm.exp(16777216) <= m.exp(4294967296)', 'm.exp(2**37) > m.exp(2**36)',
    'm.exp(-16777216) < m.Float("1e-1000")', 'm.exp(2**40) - m.exp(2**40 - 1) > 0',
]

FINITENESS = [expression + '.' + attribute
              for expression in ('(1+m.I)', '(2*m.I)', '(m.Symbol("r", real=True)+m.I)', '(m.I/2)',
                                 '(m.sqrt(2)*m.I)', '(m.Symbol("f", finite=True)*m.I)',
                                 '(1/m.Symbol("f", finite=True))', '(m.Symbol("f", finite=True)**2)',
                                 '(m.Symbol("n", nonzero=True, finite=True)**-2)', '(m.pi*m.I)')
              for attribute in ('is_finite', 'is_extended_real', 'is_real')]

KINDS = ['x.kind', '(x+1).kind', '(2*x).kind', '(x**2).kind', 'm.sin(x).kind', 'm.exp(x).kind', 'm.Integer(2).kind',
         'm.pi.kind', 'm.I.kind', 'm.exp(x).exp', 'm.exp(x).base', 'm.exp_polar(x).exp', 'm.exp_polar(x).kind', 'm.zoo.kind', 'e.kind', 't.kind', 'a.kind', 'm.Tuple(1, 2).kind',
         'm.Mod(x, 2).kind', 'm.Derivative(m.Function("g")(x), x).kind', 'm.Integral(x, x).kind',
'm.Matrix([[x]]).kind', 'm.Matrix([[1, 2]]).kind',
         'm.Matrix([[1.5]]).kind', 'm.ImmutableMatrix([[x, 1]]).kind', '(2*m.Matrix([[x]])).kind']

SREPR = ['m.srepr(m.Float(3.0))', 'm.srepr(m.Float("1.5", 30))', 'm.srepr(x + 1.5)', 'm.srepr(m.Float(0))',
         'm.srepr(m.Float(1e20))', 'm.srepr(m.Float(-2.5e-7))', 'm.srepr(m.Float("0.1", precision=50))',
         'm.srepr(m.Float("0.1", precision=52))', 'm.srepr(m.Float("0.1", precision=54))', 'm.srepr(m.Float(1)/3)']

CASES = DERIVED + BOOLEAN_ATOMS + OPERANDS + INEQUALITIES + FINITENESS + KINDS + SREPR


@pytest.mark.parametrize('expression', CASES)
def test_protocol_matches_sympy(expression):
    assert outcome(expression, 'symox') == outcome(expression, 'sympy')


def test_derived_relations_are_unevaluated():
    x = symox.Symbol('x')
    assert str(symox.Lt(x, 1).reversed) == '1 > x'
    assert str(symox.Lt(1, 2, evaluate=False).negated) == '1 >= 2'
    assert str(symox.Lt(1, 2, evaluate=False).reversedsign) == '-1 > -2'
    assert symox.Lt(1, 2, evaluate=False).weak is symox.true
