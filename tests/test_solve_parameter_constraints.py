"""An unknown parameter constraint is not a proof of no solutions."""
import pytest
import symox as sx
import sympy as sp


def foreign(value):
    if isinstance(value, dict):return {foreign(k):foreign(v) for k,v in value.items()}
    if isinstance(value, (tuple,list,set)):return type(value)(foreign(v) for v in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


@pytest.mark.parametrize('source', [
    '[x-y,y-3]', '[x-1,y-3]', '[x*x-2,y-3]', '[x*x-2,sin(y)]',
])
@pytest.mark.parametrize('flags', [{}, {'dict':True}, {'set':True}])
def test_generic_candidates_retain_undecided_parameter_constraints(source, flags):
    x,y=sx.symbols('x y')
    equations=eval(source,dict(x=x,y=y,sin=sx.sin))
    assert foreign(sx.solve(equations,x,**flags))==sp.solve(foreign(equations),x.to_sympy(),**flags)


@pytest.mark.parametrize('source', ['[x-y,1]', '[x*x-2,1]'])
def test_proved_constant_inconsistency_still_gives_no_solutions(source):
    x,y=sx.symbols('x y')
    assert sx.solve(eval(source,dict(x=x,y=y)),x)==[]


def test_assumption_proved_nonzero_parameter_is_a_real_inconsistency():
    x=sx.Symbol('x');p=sx.Symbol('p',nonzero=True)
    assert sx.solve([x-1,p],x)==[]
    assert sx.solve([x*x-2,p],x)==[]
