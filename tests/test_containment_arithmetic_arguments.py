"""Containment follows explicit expression arguments, not compact storage."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('formula', ['3*x', 'x+3', '3*x+y', 'x+3*y', 'x*y',
    'x+y', 'x+1', 'x**3', 'x/3', '3*x*exp(-x)*y', '-x', '-x-y',
    '2*x**3', 'x**(y+1)', 'Rational(1, 3)*x'])
@pytest.mark.parametrize('pattern', [None, 0, 1, -1, 2, 3, 4])
def test_literal_containment_matches_reference(formula, pattern):
    native, reference = sx.sympify(formula), sp.sympify(formula)
    assert native.has(pattern) == reference.has(pattern)
    assert native.has(None, pattern) == reference.has(None, pattern)


def test_none_does_not_prevent_matching_another_pattern_or_class():
    x = sx.Symbol('x')
    assert not x.has(None)
    assert not x.has()
    assert x.has(None, x)
    assert x.has(None, sx.Symbol)
    assert (3*x).has(3)
