"""CSE recipe reconstruction must retain observable attribute protocols."""
import importlib
import pytest
import symox as s


@pytest.mark.parametrize('customization',['root_getattribute','base_getattribute','base_getattr'])
def test_recipe_cse_observes_the_same_attribute_callbacks_as_the_existing_route(monkeypatch,customization):
    module=importlib.import_module('symox.simplify.cse_main')
    x,y=s.symbols('x y')
    power=s.diff(s.exp(2*x),x,524289).as_independent(x)[0]
    expressions=[power*(x+1),power*(y+1)]
    calls=[]
    if customization=='root_getattribute':
        kind=type(power)
        def getattribute(self,name):
            if name=='func':calls.append(name)
            return object.__getattribute__(self,name)
        monkeypatch.setattr(kind,'__getattribute__',getattribute)
    elif customization=='base_getattribute':
        kind=type(s.Integer(2))
        def getattribute(self,name):
            if name in ('_symox_matrixexpr','_symox_sparse'):calls.append(name)
            return object.__getattribute__(self,name)
        monkeypatch.setattr(kind,'__getattribute__',getattribute)
    else:
        kind=type(s.Integer(2))
        def getattr_(self,name):
            if name in ('_symox_matrixexpr','_symox_sparse'):calls.append(name)
            raise AttributeError(name)
        monkeypatch.setattr(kind,'__getattr__',getattr_,raising=False)
    assert not module._stock_recipe_reconstruction(power,list(power.args))
    actual=s.cse(expressions);observed=list(calls);calls.clear()
    monkeypatch.setattr(module,'_stock_recipe_reconstruction',lambda *args:False)
    expected=s.cse(expressions)
    assert actual==expected
    assert observed==calls and observed
