"""Reused exact Add construction keeps presentation and foreign protocols."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize("size", [10,100,1000])
def test_repeated_monomial_sums_keep_exact_coefficients(size):
    x,rx=sx.Symbol("x"),sp.Symbol("x")
    terms=[i*x**i for i in range(size)]
    expected=sp.Add(*(i*rx**i for i in range(size)))
    for _ in range(3):
        result=sx.Add(*terms)
        assert result.to_sympy()==expected
    assert sx.Add(x,-x)==0
    assert sx.Add(x,-x)==0


def test_held_evaluation_keeps_input_order_after_a_warmed_sum():
    x,y,z=sx.symbols("x y z")
    terms=(x*x,2*y,3*z*z*z)
    sx.Add(*terms);sx.Add(*terms)
    value=sx.Add(*reversed(terms),evaluate=False)
    assert value.args==tuple(reversed(terms))
    with sx.evaluate(False):
        assert sx.Add(*terms).args==terms


def test_float_precision_and_distinct_symbol_assumptions_are_preserved():
    x=sx.Symbol("x")
    for digits in [30,100,30]:
        value=sx.Float("1.25",digits)
        result=sx.Add(x,value)
        assert next(iter(result.atoms(sx.Float)))._prec==value._prec
    for positive in [True,False,True]:
        symbol=sx.Symbol("v",positive=positive)
        result=sx.Add(symbol,symbol)
        assert next(iter(result.free_symbols)).is_positive is positive


def test_foreign_conversion_is_not_hashed_or_repeated_by_constructor_reuse():
    x=sx.Symbol("x");calls=[]
    class Foreign:
        def _symox_(self):
            calls.append(True)
            return x
        def __hash__(self):
            raise AssertionError("foreign hash")
    value=Foreign()
    assert sx.Add(value,1)==x+1
    assert calls==[True]


def test_alpha_equal_binders_keep_their_own_presentations():
    x,y,t,z=sx.symbols("x y t z");f=sx.Function("f")
    a,b=sx.Subs(f(x),x,t),sx.Subs(f(y),y,t)
    assert a==b and a.args!=b.args
    for value in [a,b,a,b]:
        result=sx.Add(value,z)
        assert next(arg for arg in result.args if arg.class_name=="Subs").args==value.args
