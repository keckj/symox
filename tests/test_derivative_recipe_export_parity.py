"""Recipe export must not change ordinary canonical or explicit-held exports."""
from collections import Counter
import pickle
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source',[
    'arg(sqrt(-2+I))',
    'cos(1)/sin(1)',
    'atan(cos(atan(Rational(1,2))/2)/sin(atan(Rational(1,2))/2))',
    'sin(1)**-2+cos(1)**-3',
    '1/(log(2)+sqrt(3))',
])
@pytest.mark.parametrize('method',['method','function'])
def test_canonical_closed_reciprocals_keep_reference_identity_and_hash(source,method):
    actual=eval(source,vars(sx));expected=eval(source,vars(sp))
    converted=actual.to_sympy() if method=='method' else sx.to_sympy(actual)
    # The failing R4 expression already had the same str and srepr. Require
    # actual SymPy identity semantics instead of accepting printed agreement.
    assert converted==expected
    assert hash(converted)==hash(expected)
    assert converted.is_commutative==expected.is_commutative


@pytest.mark.parametrize('context',['root','canonical_sum','held_product'])
def test_explicit_evaluate_false_exports_follow_the_r3b_convention(context):
    x=sx.Symbol('x');X=sp.Symbol('x');raw=sx.Pow(2,7,evaluate=False)
    if context=='root':
        actual=raw;expected=sp.Pow(2,7,evaluate=False)
    elif context=='canonical_sum':
        actual=raw+x;expected=X+128
    else:
        actual=sx.Mul(raw,x,evaluate=False)
        expected=sp.Mul(sp.Pow(2,7,evaluate=False),X,evaluate=False)
    assert actual.to_sympy()==expected
    assert hash(actual.to_sympy())==hash(expected)


@pytest.mark.parametrize('order',[2**32+1,2**128+1])
@pytest.mark.parametrize('slope',['2','Rational(2,3)','2+I'])
@pytest.mark.parametrize('context',['coefficient','product','sum'])
def test_actual_large_derivative_recipes_remain_exact_and_unmaterialized(order,slope,context):
    x=sx.Symbol('x');X=sp.Symbol('x')
    a=eval(slope,vars(sx));A=eval(slope,vars(sp))
    result=sx.diff(sx.exp(a*x),x,order)
    coefficient,phase=result.as_independent(x)
    power=sp.Pow(A,order,evaluate=False)
    if context=='coefficient':actual=coefficient;expected=power
    elif context=='product':actual=result;expected=sp.Mul(power,sp.exp(A*X),evaluate=False)
    else:actual=coefficient+x;expected=sp.Add(power,X,evaluate=False)
    assert actual._has_deferred_numeric_power
    for candidate in (actual,pickle.loads(pickle.dumps(actual))):
        converted=candidate.to_sympy()
        if context=='coefficient':
            assert converted==expected and hash(converted)==hash(expected)
        else:
            # These NEW recipes intentionally use held Add/Mul. Their stored
            # argument order may differ; compare every exact commutative term
            # or factor, including the literal huge power. No sampling.
            assert converted.func is expected.func
            assert all(value.is_commutative is True for value in converted.args)
            assert Counter(converted.args)==Counter(expected.args)
