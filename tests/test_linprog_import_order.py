"""Optional stock Tensor dispatch keeps the complete QQ preparation native."""
import subprocess
import sys

import pytest


@pytest.mark.parametrize('order',['linprog-first','optional-first'])
@pytest.mark.parametrize('modules',[('diffgeom','vector'),('vector','diffgeom')])
def test_exact_linprog_dispatch_and_real_callbacks_after_diffgeom_import(order,modules):
    source=r'''
import functools,importlib,sys
import symox as S
if ORDER=='linprog-first':
    lp=importlib.import_module('symox.solvers.simplex')
    from symox import stats
for module in MODULES:importlib.import_module('symox.'+module)
if ORDER!='linprog-first':
    lp=importlib.import_module('symox.solvers.simplex')
    from symox import stats
import math,symox.stats._joint as joint
cache=joint.lru_cache

def no_recurrence(*args,**kwargs):raise AssertionError('stock wrapper disabled native Gaussian recurrence')
joint.lru_cache=no_recurrence
assert int(joint._normal_raw((1200,),[S.Integer(0)],S.ImmutableMatrix([[1]])))==math.prod(range(1,1200,2))
joint.lru_cache=cache
codes={lp._matrix_rows.__code__,lp._bounds.__code__,lp._simplex.__code__,S.sympify.__code__}
calls=[]
def profile(frame,event,arg):
    if event=='call' and frame.f_code in codes:calls.append(frame.f_code.co_name)
sys.setprofile(profile)
try:
    assert lp.linprog([1,2],[[-1,-1]],[-3],[[1,-1]],[1],bounds={0:(None,None),1:(0,4)})==(4,[2,1])
finally:sys.setprofile(None)
assert calls==[],calls

# A callback can copy every public attribute of the stock wrapper. Its code
# still differs, so neither plain nor @wraps overrides may bypass preparation.
original=S.simplify
native=S._core._native_linprog
for decorated in (False,True):
    calls=[]
    def override(value,*args,**kwargs):
        calls.append(value)
        raise KeyboardInterrupt
    if decorated:override=functools.wraps(original)(override)
    S.simplify=override
    def forbidden(*args):raise AssertionError('real callback reached native preparation')
    S._core._native_linprog=forbidden
    try:lp.linprog([1],[[1]],[2])
    except KeyboardInterrupt:assert len(calls)==1
    else:raise AssertionError('real callback was skipped')
    S._core._native_linprog=native
    S.simplify=original
'''.replace('ORDER',repr(order)).replace('MODULES',repr(modules))
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)


def test_unrelated_converter_keeps_native_preparation(monkeypatch):
    import importlib
    import symox as S
    from symox.core.sympify import converter
    lp=importlib.import_module('symox.solvers.simplex')
    class Unrelated:pass
    def forbidden(value):raise AssertionError('unrelated converter invoked')
    monkeypatch.setitem(converter,Unrelated,forbidden)
    calls=[];original=S._core._native_linprog
    def native(*args):
        result=original(*args);calls.append(result);return result
    monkeypatch.setattr(S._core,'_native_linprog',native)
    assert lp.linprog([1],bounds=(2,3))==(2,[2])
    assert len(calls)==1 and calls[0][0]=='ok'


@pytest.mark.parametrize('key',[int,object])
def test_applicable_integer_converter_stays_observable(key,monkeypatch):
    import importlib
    import symox as S
    from symox.core.sympify import converter
    lp=importlib.import_module('symox.solvers.simplex')
    calls=[]
    def interrupted(value):calls.append(value);raise KeyboardInterrupt
    monkeypatch.setitem(converter,key,interrupted)
    with pytest.raises(KeyboardInterrupt):lp.linprog([1])
    assert calls==[1]


def test_converter_eligibility_never_hashes_or_compares_custom_keys():
    from symox._native_dispatch import stock_integer_converters
    calls=[]
    class Collision:
        def __hash__(self):calls.append('hash');return hash(int)
        def __eq__(self,other):calls.append('eq');raise AssertionError('custom equality')
    key=Collision();registry={key:lambda value:value};calls.clear()
    assert not stock_integer_converters(registry)
    assert calls==[]


def test_wrapper_chain_rejects_cycles_and_hidden_user_callbacks():
    source=r'''
from types import FunctionType
import symox as S
import symox.diffgeom,symox.vector
from symox._general_simplify import simplify
from symox._native_dispatch import stock_scalar_dispatch
wrapper=S.simplify
index=wrapper.__code__.co_freevars.index('original')

def clone(delegate):
    cells=list(wrapper.__closure__)
    cell=(lambda value:lambda:value)(delegate).__closure__[0]
    cells[index]=cell
    value=FunctionType(wrapper.__code__,wrapper.__globals__,closure=tuple(cells))
    value.__dict__.update(wrapper.__dict__)
    value.__wrapped__=delegate
    return value,cell

def callback(value):raise KeyboardInterrupt
hidden,_=clone(callback)
assert not stock_scalar_dispatch(hidden,simplify,'simplify')
cycle,cell=clone(None)
cell.cell_contents=cycle
cycle.__wrapped__=cycle
assert not stock_scalar_dispatch(cycle,simplify,'simplify')
assert stock_scalar_dispatch(wrapper,simplify,'simplify')
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
