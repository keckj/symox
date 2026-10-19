"""Undefined applications have specific native classes across reconstruction."""
import copy
import pickle

import pytest
import symox as sx
from symox.core.function import UndefinedFunction


@pytest.mark.parametrize('assumptions',[{}, {'real':True}, {'positive':True}])
def test_application_class_preserves_native_head_identity(assumptions):
    x=sx.Symbol('x');f=sx.Function('f',nargs=1,**assumptions);value=f(x)
    assert isinstance(type(value),UndefinedFunction)
    assert issubclass(type(type(value)),UndefinedFunction)
    assert type(value).__name__=='f'
    assert type(value)(x)==value
    assert value.func==f and value.is_UndefinedFunction
    assert not callable(value)
    assert type(sx.Function('g')(x)) is not type(value)
    assert type(f(x+1)) is type(value)
    for restored in [copy.copy(value),copy.deepcopy(value),pickle.loads(pickle.dumps(value)),
                     next(iter((value+1).atoms(sx.Function))),
                     value.subs(x,x+1).subs(x+1,x)]:
        assert restored==value
        assert type(restored) is type(value)
    assert isinstance(f,UndefinedFunction)
    assert isinstance(type(value),sx.FunctionClass)
    assert not isinstance(type(x),UndefinedFunction)


def test_application_pattern_remains_distinct_from_its_head_class():
    x=sx.Symbol('x');f=sx.Function('f')
    integer=sx.Wild('integer',properties=[lambda value:value.is_Integer])
    expression=f(2)+f(x)
    assert expression.find(f(integer))=={f(2)}
    assert expression.find(type(f(x)))=={f(2),f(x)}
    assert expression.replace(f(integer),integer+1)==3+f(x)
