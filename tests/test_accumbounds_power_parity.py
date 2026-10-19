"""Accumulation bounds cover fractional and infinite exponent domains."""
import pytest
import sympy as sp
import symox as sx


def convert(value):
    return value.to_sympy() if hasattr(value,'to_sympy') else sp.sympify(value)


@pytest.mark.parametrize('bounds',[(0,1),(1,2),(2,4),(-3,-2),(-2,-1),(-1,0),(-1,1),(-sx.Rational(1,2),sx.Rational(1,2))])
@pytest.mark.parametrize('exponent',[sx.oo,-sx.oo])
def test_infinite_exponent_limits_match_reference(bounds,exponent):
    actual=sx.AccumBounds(*bounds)**exponent
    expected=sp.AccumBounds(*(convert(v) for v in bounds))**convert(exponent)
    assert convert(actual)==expected


@pytest.mark.parametrize('bounds',[(0,1),(1,2),(2,4),(0,4)])
@pytest.mark.parametrize('exponent',[sx.Rational(1,2),sx.Rational(2,3),sx.Rational(-1,2),sx.Rational(-2,3),sx.Float('0.5')])
def test_nonnegative_base_fractional_powers(bounds,exponent):
    actual=sx.AccumBounds(*bounds)**exponent
    expected=sp.AccumBounds(*bounds)**convert(exponent)
    assert convert(actual)==expected


def test_interval_exponents_and_metadata_cache_clear():
    actual=sx.AccumBounds(1,2)**sx.AccumBounds(-1,2)
    expected=sp.AccumBounds(1,2)**sp.AccumBounds(-1,2)
    assert convert(actual)==expected
    x=sx.Symbol('x');root=sx.CRootOf(x**5-x-1,0)
    before=root.evalf(40)
    sx.CRootOf.clear_cache()
    assert root.evalf(40)==before
