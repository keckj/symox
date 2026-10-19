"""Function-head metadata must correspond to executable keyword bindings."""
import inspect
import subprocess
import sys
import pytest
import sympy as sp
import symox as sx
from symox._function_signatures import _PARAMETERS


@pytest.mark.parametrize('name',sorted(_PARAMETERS))
def test_native_heads_bind_required_keywords_and_reject_invalid_shapes(name):
    function=sx._make(name)
    contract=inspect.signature(getattr(sp,name))
    args=[];kwargs={}
    for p in contract.parameters.values():
        if p.kind==p.VAR_POSITIONAL:
            args.extend([2,3,1] if name=='lerchphi' else [2,3]);break
        if p.default is p.empty:
            value=3 if p.name in ('n','p') else 2
            args.append(value);kwargs[p.name]=value
    # Held construction isolates the call contract from simplification rules.
    positional=function(*args,evaluate=False)
    if kwargs:
        keyword=function(**kwargs,evaluate=False)
        assert keyword.srepr()==positional.srepr()
        first=next(iter(kwargs))
        with pytest.raises(TypeError):function(*args,**{first:kwargs[first]},evaluate=False)
    with pytest.raises(TypeError):function(*args,_unknown_api_option=True,evaluate=False)
    if any(p.default is p.empty and p.kind!=p.VAR_POSITIONAL for p in contract.parameters.values()):
        with pytest.raises(TypeError):function(evaluate=False)


@pytest.mark.parametrize('name,kwargs,reference_args',[
    ('sin',{'arg':2},(2,)),('atan2',{'y':2,'x':3},(2,3)),
    ('LambertW',{'x':2,'k':None},(2,)),('DiracDelta',{'arg':2,'k':0},(2,0)),
    ('binomial',{'n':5,'k':2},(5,2)),('besselj',{'nu':2,'z':3},(2,3)),
    ('harmonic',{'n':3,'m':None},(3,)),('beta',{'x':2,'y':3},(2,3)),
    ('zeta',{'s':2,'a':None},(2,)),
])
def test_public_keyword_calls_match_reference(name,kwargs,reference_args):
    actual=getattr(sx,name)(**kwargs)
    # SymPy's function classes display parameter names but often reject their
    # keyword spelling. Native keyword calls extend the positional operation.
    expected=getattr(sp,name)(*reference_args)
    assert actual.to_sympy()==expected


def test_function_keyword_adapter_has_no_reference_runtime_dependency():
    code='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('reference import attempted')
sys.meta_path.insert(0,Block())
import symox as s
x=s.Symbol('x')
assert s.sin(arg=x)==s.sin(x)
assert s.binomial(n=5,k=2)==10
assert s.LambertW(x=2,k=None)==s.LambertW(2)
'''
    subprocess.run([sys.executable,'-c',code],check=True)
