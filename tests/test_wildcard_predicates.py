"""Predicate constraints participate in native matching and backtracking."""
import copy
import pickle

import pytest
import sympy as sp
import symox as sx


def integer_value(value):
    return value.is_integer


def symbol_free(value):
    return not value.free_symbols


@pytest.mark.parametrize('value',[-2,0,3,sp.Rational(1,2),sp.Symbol('x')])
def test_predicate_binding_matches_reference(value):
    pattern=sx.Wild('a',properties=[integer_value])
    reference=sp.Wild('a',properties=[integer_value])
    actual=sx.sympify(value).match(pattern)
    expected=value.match(reference) if isinstance(value,sp.Basic) else sp.sympify(value).match(reference)
    assert (actual is None)==(expected is None)
    if actual is not None:assert actual=={pattern:sx.sympify(value)}


def test_predicate_identity_survives_expression_copy_and_pickle():
    x=sx.Symbol('x')
    a=sx.Wild('a',exclude=[x],properties=[integer_value,symbol_free],positive=True)
    assert a.properties==(integer_value,symbol_free)
    assert a==sx.Wild('a',exclude=[x],properties=[integer_value,symbol_free],positive=True)
    assert a!=sx.Wild('a',exclude=[x],properties=[symbol_free,integer_value],positive=True)
    assert a!=sx.Wild('a',exclude=[x],properties=[integer_value],positive=True)
    assert a!=sx.Wild('a',exclude=[x],positive=True)
    assert sx.sympify(a.to_sympy())==a
    for expression in [a,2*a+sx.sin(a)]:
        for clone in [copy.copy(expression),copy.deepcopy(expression),pickle.loads(pickle.dumps(expression))]:
            assert clone==expression
            assert all(w.properties==(integer_value,symbol_free) for w in clone.atoms(sx.Wild))
    assert sx.Integer(3).match(a)=={a:sx.Integer(3)}
    assert x.match(a) is None


def test_predicate_rejection_backtracks_to_another_function_argument():
    f=sx.Function('f');a=sx.Wild('a',properties=[lambda value:value.is_positive])
    b=sx.Wild('b',properties=[lambda value:value.is_negative])
    assert (f(-2)+f(3)).match(f(a)+f(b))=={a:sx.Integer(3),b:sx.Integer(-2)}


def test_commutative_groups_and_neutral_values_remain_available():
    x,y=sx.symbols('x y')
    a=sx.Wild('a',properties=[lambda value:value.has(x) and value.has(y)])
    b=sx.Wild('b',properties=[integer_value])
    assert (x+y+1).match(a+b)=={a:x+y,b:sx.Integer(1)}
    z=sx.Wild('zero',properties=[lambda value:value==0])
    rest=sx.Wild('rest')
    assert (x+y).match(z+rest)=={z:sx.Integer(0),rest:x+y}
    assert (3*x*y).match(a*b)=={a:x*y,b:sx.Integer(3)}
    unit=sx.Wild('unit',properties=[lambda value:value==1])
    assert (x*y).match(unit*rest)=={unit:sx.Integer(1),rest:x*y}


def test_predicates_short_circuit_and_exclusions_run_first():
    calls=[]
    def reject(value):calls.append(value);return False
    def forbidden(value):raise AssertionError('predicate should not be called')
    a=sx.Wild('a',properties=[reject,forbidden])
    assert sx.Integer(3).match(a) is None and calls==[sx.Integer(3)]
    x=sx.Symbol('x');a=sx.Wild('a',exclude=[x],properties=[forbidden])
    assert (x+1).match(a) is None


@pytest.mark.parametrize('operation',['match','matches','replace','find','collect'])
def test_predicate_errors_propagate_unchanged(operation):
    error=RuntimeError('predicate sentinel')
    def failure(value):raise error
    a=sx.Wild('a',properties=[failure]);x=sx.Symbol('x')
    with pytest.raises(RuntimeError) as caught:
        if operation=='match':x.match(a)
        elif operation=='matches':a.matches(x)
        elif operation=='replace':x.replace(a,sx.Integer(2))
        elif operation=='find':x.find(a)
        else:sx.collect(x,a)
    assert caught.value is error


def test_reciprocal_constraints_cover_canonical_numbers():
    x=sx.Symbol('x');a=sx.Wild('a',properties=[symbol_free])
    for coefficient in [sx.Integer(2),sx.Rational(2,3),sx.Integer(-3)]:
        assert (coefficient*x).match(x/a)=={a:1/coefficient}
    assert sx.Integer(0).match(1/a) is None
    assert (x*sx.Symbol('y')).match(x/a) is None


def test_replacement_and_find_keep_constraint_metadata():
    x=sx.Symbol('x');f=sx.Function('f');a=sx.Wild('a',properties=[integer_value])
    expression=f(2)+f(x)
    assert expression.find(f(a))=={f(2)}
    assert expression.replace(f(a),lambda a:a+1)==3+f(x)


def test_noncallable_predicates_are_rejected_when_evaluated():
    a=sx.Wild('a',properties=[3])
    assert a.properties==(3,)
    with pytest.raises(TypeError):sx.Integer(1).match(a)


def test_unhashable_predicates_keep_reference_hash_errors():
    class Predicate:
        __hash__=None
        def __init__(self,result):self.result=result
        def __call__(self,value):return self.result
    reject=sx.Wild('a',properties=[Predicate(False)])
    assert sx.Integer(1).match(reject) is None
    accept=sx.Wild('a',properties=[Predicate(True)])
    with pytest.raises(TypeError):hash(accept)
    with pytest.raises(TypeError):sx.Integer(1).match(accept)
    with pytest.raises(TypeError):sx.Integer(1).replace(accept,lambda a:a+1)
    with pytest.raises(TypeError):sx.Integer(1).find(accept)


def test_matrix_pickle_preserves_predicate_bindings():
    a=sx.Wild('a',properties=[integer_value])
    matrix=sx.ImmutableMatrix([[a,sx.sin(a)]])
    restored=pickle.loads(pickle.dumps(matrix))
    assert restored==matrix
    assert restored[0,0].properties==(integer_value,)
    assert sx.Integer(2).match(restored[0,0])=={a:sx.Integer(2)}
    assert sx.Rational(1,2).match(restored[0,0]) is None


def test_legacy_linear_fallback_enforces_predicates():
    x=sx.Symbol('x')
    accept=sx.Wild('a',properties=[lambda value:value.has(x)])
    reject=sx.Wild('a',properties=[symbol_free])
    assert (2/x).match(accept*x,old=True)=={accept:2/x**2}
    assert (2/x).match(reject*x,old=True) is None
