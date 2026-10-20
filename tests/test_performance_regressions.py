"""Performance fixes retain exact answers, limits and presentation semantics."""
import random

import pytest
import symox as sx
import sympy as sp
from symox.combinatorics import Permutation
from sympy.combinatorics import Permutation as Reference


@pytest.mark.parametrize("degree", [20,50,100])
def test_linear_factor_expansion_counts_distinct_univariate_monomials(degree):
    x,rx = sx.Symbol("x"),sp.Symbol("x")
    value = sx.prod(x+i for i in range(1,degree+1))
    expected = sp.expand(sp.prod(rx+i for i in range(1,degree+1)))
    actual = value.expand(max_terms=degree+1)
    assert actual.to_sympy() == expected
    assert value.expand().to_sympy() == expected
    if degree == 20:
        assert sx.factor(actual).to_sympy() == sp.factor(expected)


def test_univariate_power_support_and_warm_expansion_limits():
    x = sx.Symbol("x")
    value = (x*x+x+1)**50
    expected = sp.expand((sp.Symbol("x")**2+sp.Symbol("x")+1)**50)
    assert value.expand(max_terms=101).to_sympy() == expected
    assert value.expand(max_terms=101).to_sympy() == expected
    with pytest.raises(ValueError,match="limit"):
        value.expand(max_terms=100)
    y = sx.Symbol("y")
    with pytest.raises(ValueError,match="limit|budget"):
        ((x+y)**10000).expand()


def test_repeated_expansion_keeps_options_and_alpha_equal_binder_names():
    x,y,t = sx.symbols("x y t")
    f,g = sx.Function("f"),sx.Function("g")
    left,right = sx.Subs(f(x),x,t),sx.Subs(f(y),y,t)
    assert left == right and left.args != right.args
    a,b = g(left),g(right)
    for _ in range(3):
        assert a.expand().args[0].args == left.args
        assert b.expand().args[0].args == right.args
    value = sx.sin((x+y)**3)
    assert value.expand(deep=False) == value
    assert value.expand().to_sympy() == sp.sin((sp.Symbol("x")+sp.Symbol("y"))**3).expand()
    assert value.expand(deep=False) == value


@pytest.mark.parametrize("degree", [0,1,2,7,100,1000])
def test_permutation_fast_operations_keep_images_and_arbitrary_integer_powers(degree):
    rng = random.Random(183)
    a,b = list(range(degree)),list(range(degree))
    rng.shuffle(a);rng.shuffle(b)
    p,q = Permutation(a),Permutation(b)
    P,Q = Reference(a),Reference(b)
    assert (p*q).array_form == (P*Q).array_form
    assert (~p).array_form == (~P).array_form
    for n in [-10**30-1,-65,-64,-17,-2,-1,0,1,2,3,4,17,63,64,65,10**30+1]:
        result = p**n
        assert result.array_form == (P**n).array_form
        assert result.size == degree
        assert result.rank() == Reference(result.array_form).rank()


def test_permutation_results_clear_cached_metadata_and_preserve_subclass_construction():
    p = Permutation([1,2,0,3])
    p.rank();p.parity();p.cycle_structure
    result = p*p
    assert result.rank() == Reference(result.array_form).rank()
    assert result.parity() == Reference(result.array_form).parity()
    seen = []
    class Custom(Permutation):
        def __init__(self,*args,**kwargs):
            seen.append(True)
            super().__init__(*args,**kwargs)
    custom = Custom([1,0])
    for result in [custom*custom,~custom,custom**17]:
        assert isinstance(result,Custom)
    assert len(seen) == 4
    for invalid in [[0,0],[1,2],[-1,0]]:
        with pytest.raises(ValueError):Permutation(invalid)


def test_repeated_substitution_keeps_replacement_precision_and_bound_names():
    x,y,t = sx.symbols("x y t")
    expression = sx.Function("f")(x)
    for digits in [30,100,30]:
        result = expression.subs(x,sx.Float("1.25",digits))
        assert result.args[0]._prec == sx.Float("1.25",digits)._prec
    a,b = sx.Subs(expression,x,t),sx.Subs(sx.Function("f")(y),y,t)
    assert a == b and a.args != b.args
    pair = sx.Function("g")(x)
    for _ in range(2):
        assert pair.subs(x,a).args[0].args == a.args
        assert pair.subs(x,b).args[0].args == b.args


def test_repeated_numeric_evaluation_keeps_precision_substitutions_and_chop():
    expression = sx.exp(sx.pi)*sx.sqrt(2)+sx.log(3)
    reference = sp.exp(sp.pi)*sp.sqrt(2)+sp.log(3)
    for digits in [30,100,300,100,30]:
        result = expression.evalf(digits).to_sympy()
        assert result._prec == sp.Float(1,digits)._prec
        assert abs(sp.N((result-reference.evalf(digits+20))/reference,digits+10)) < sp.Rational(1,10**(digits-2))
    x = sx.Symbol("x")
    value = sx.exp(x)
    for integer in [1,2,1]:
        assert abs(float(value.evalf(50,subs={x:integer}))-float(sp.exp(integer))) < 1e-12
    tiny = sx.Rational(1,10**1000)
    value = 1+sx.I*tiny
    assert value.evalf(60,chop=True).as_real_imag()[1] == 0
    assert value.evalf(60).as_real_imag()[1] != 0
    assert value.evalf(60,chop=True).as_real_imag()[1] == 0


def test_repeated_number_theory_keeps_exact_answers_and_independent_factor_dicts():
    integer = 2**64+1
    expected = sp.factorint(integer)
    first = sx.factorint(integer)
    assert first == expected
    first.clear()
    assert sx.factorint(integer) == expected
    assert sx.factorint(-integer) == expected | {-1:1}
    for n in [10**5,10**6,10**5]:
        assert sx.primepi(n).to_sympy() == sp.primepi(n)
    value = (2**127-1)**12
    for big in [True,False,True]:
        assert sx.perfect_power(value,big=big) == sp.perfect_power(value,big=big)
    for candidates in [[2],[3],[4],[5],None]:
        assert sx.perfect_power(value,candidates=candidates) == sp.perfect_power(value,candidates=candidates)


@pytest.mark.parametrize("digits", [2,5,15,30,50,100,300])
def test_faster_logarithm_constant_preserves_requested_accuracy(digits):
    actual = sx.log(2).evalf(digits).to_sympy()
    expected = sp.log(2).evalf(digits+30)
    assert abs(sp.N(actual-expected,digits+20)) < sp.Rational(1,10**digits)


def test_characteristic_polynomial_fallback_keeps_context_independent_conversion():
    x,t = sx.symbols("x t")
    matrix = sx.Matrix([[sx.sin(x),1],[0,sx.log(x)]])
    expected = matrix.charpoly(t).as_expr().to_sympy()
    with sx.evaluate(False):
        actual = matrix.charpoly(t).as_expr().to_sympy()
    assert actual == expected
    coefficients = [sx.Integer(1),-sx.sin(x),sx.Integer(0)]
    baseline = sx.PurePoly._from_coeffs(coefficients,t,'EX').as_expr().to_sympy()
    with sx.evaluate(False):
        result = sx.PurePoly._from_coeffs(coefficients,t,'EX').as_expr().to_sympy()
    assert result == baseline
