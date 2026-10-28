"""Exact QQ/ZZ rounding; preserve the old numeric protocol elsewhere."""
import math
import subprocess
import sys
import warnings

import pytest
import symox as sx
from symox.polys.domains import QQ, ZZ, RR, CC, EX, GF


PAIRS=[(0,1),(1,2),(-1,2),(3,2),(-3,2),(6,3),(-6,3),
       (((1<<100)+1),3),(-((1<<100)+1),3),
       ((1<<160)+1,1<<80),((1<<160)-1,1<<80),
       (-((1<<160)+1),1<<80),(-((1<<160)-1),1<<80),
       ((1<<5000)+1,(1<<5000)+2),(-((1<<5000)+1),(1<<5000)+2)]


@pytest.mark.parametrize('n,d',PAIRS)
def test_rational_rounding_has_exact_integer_certificates(n,d):
    value=QQ(n,d)
    floor=math.floor(value);ceil=math.ceil(value)
    assert type(floor) is int and type(ceil) is int
    quotient,remainder=divmod(n,d)
    assert floor==quotient and ceil==quotient+(remainder!=0)
    assert floor*d<=n<(floor+1)*d
    assert (ceil-1)*d<n<=ceil*d


@pytest.mark.parametrize('n',[0,1,-1,(1<<100)+1,-((1<<100)+1),(1<<5000)+1,-((1<<5000)+1)])
def test_integer_rounding_never_converts_to_float(n):
    value=ZZ(n)
    assert math.floor(value)==math.ceil(value)==n
    assert type(math.floor(value)) is int


def outcome(operation,value):
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter('always')
        try:result=('value',operation(value))
        except BaseException as error:result=('error',type(error).__name__,str(error))
    return result,[(type(w.message).__name__,str(w.message)) for w in observed]


def compare_old(monkeypatch,operation,value,calls=None):
    new=outcome(operation,value)
    trace=list(calls) if calls is not None else None
    if calls is not None:calls.clear()
    with monkeypatch.context() as patch:
        patch.delattr(sx._core.DomainElement,'__'+operation.__name__+'__')
        old=outcome(operation,value)
    assert new==old
    if calls is not None:assert calls==trace
    return new


@pytest.mark.parametrize('factory',[
    lambda:RR(1.25),lambda:RR(-1.25),lambda:sx.RealField(150)(sx.Rational((1<<100)+1,3)),
    lambda:CC(1,0),lambda:CC(1,2),lambda:EX(sx.Rational(1,3)),
    lambda:EX(sx.Symbol('rounding_x')),lambda:GF(7)(2),
    lambda:QQ.algebraic_field(sx.sqrt(2))(sx.sqrt(2)),
    lambda:ZZ.poly_ring(sx.Symbol('rounding_x'))(sx.Symbol('rounding_x')),
])
@pytest.mark.parametrize('operation',[math.floor,math.ceil])
def test_other_domains_keep_old_results_or_errors(monkeypatch,factory,operation):
    value=factory()
    if value.parent().is_RealField and value.parent().precision>53:
        # This followup corrects the formerly preserved RR precision gap.
        # Keep the same witness and certify its actual stored dyadic instead.
        from fractions import Fraction
        sign,mantissa,exponent,_=value.parent().to_sympy(value)._mpf_
        literal=Fraction((-1 if sign else 1)*mantissa)*Fraction(2)**exponent
        result=operation(value)
        if operation is math.floor:assert result<=literal<result+1
        else:assert result-1<literal<=result
    else:
        compare_old(monkeypatch,operation,value)


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
@pytest.mark.parametrize('kind',['count','wrongint','wrongstr','wrongnone','interrupt','memory','value','nan','posinf','neginf','subfloat'])
def test_current_float_override_is_observed_exactly(monkeypatch,operation,kind):
    value=QQ((1<<100)+1,3);calls=[]
    class F(float):pass
    def convert(self):
        calls.append('float')
        if kind=='count':return -2.75
        if kind=='wrongint':return 2
        if kind=='wrongstr':return '2.75'
        if kind=='wrongnone':return None
        if kind=='interrupt':raise KeyboardInterrupt('float hook')
        if kind=='memory':raise MemoryError('float hook')
        if kind=='value':raise ValueError('float hook')
        if kind=='subfloat':return F(2.75)
        return {'nan':float('nan'),'posinf':float('inf'),'neginf':-float('inf')}[kind]
    monkeypatch.setattr(sx._core.DomainElement,'__float__',convert)
    compare_old(monkeypatch,operation,value,calls)


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
def test_float_descriptor_and_host_getattribute_keep_special_lookup(monkeypatch,operation):
    value=QQ(7,3);calls=[]
    def descriptor(self):
        calls.append('descriptor')
        def convert():calls.append('convert');return -2.75
        return convert
    def hostile(self,name):raise AssertionError('ordinary getattr: '+name)
    monkeypatch.setattr(sx._core.DomainElement,'__float__',property(descriptor))
    monkeypatch.setattr(sx._core.DomainElement,'__getattribute__',hostile)
    result=compare_old(monkeypatch,operation,value,calls)
    assert result[0][0]=='value'
    assert calls==['descriptor','convert']


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
def test_deleted_float_slot_keeps_index_fallback(monkeypatch,operation):
    value=QQ(7,3);calls=[]
    monkeypatch.delattr(sx._core.DomainElement,'__float__')
    def index(self):calls.append('index');return 17
    monkeypatch.setattr(sx._core.DomainElement,'__index__',index)
    assert compare_old(monkeypatch,operation,value,calls)[0]==('value',17)
    assert calls==['index']


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
def test_explicit_rounding_override_has_precedence(monkeypatch,operation):
    value=QQ(7,3);marker=object();calls=[]
    def rounding(self):calls.append('round');return marker
    def forbidden(self):raise AssertionError('float should not run')
    monkeypatch.setattr(sx._core.DomainElement,'__float__',forbidden)
    monkeypatch.setattr(sx._core.DomainElement,'__'+operation.__name__+'__',rounding)
    assert operation(value) is marker and calls==['round']


def test_float_slot_restoration_recovers_exact_path(monkeypatch):
    value=QQ((1<<100)+1,3)
    with monkeypatch.context() as patch:
        patch.setattr(sx._core.DomainElement,'__float__',lambda self:2.75)
        assert math.floor(value)==2 and math.ceil(value)==3
    assert math.floor(value)==((1<<100)+1)//3
    assert math.ceil(value)==-(-((1<<100)+1)//3)


def test_domain_element_subclass_policy_is_unchanged():
    with pytest.raises(TypeError,match='acceptable base type'):
        type('ForeignElement',(sx._core.DomainElement,),{})


def test_no_reference_runtime_and_live_slot_proof_in_fresh_process():
    script=r'''
import sys,importlib.abc,math
class NoReference(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('runtime SymPy')
sys.meta_path.insert(0,NoReference())
import symox as s
from symox.polys.domains import QQ
T=s._core.DomainElement
native=T.__dict__['__float__']
# A change installed before creating the value must not become a trust anchor.
T.__float__=lambda self:-2.25
v=QQ((1<<100)+1,3)
assert math.floor(v)==-3 and math.ceil(v)==-2
T.__float__=native
assert math.floor(v)==((1<<100)+1)//3
assert math.ceil(v)==-(-((1<<100)+1)//3)
assert not any(k=='sympy' or k.startswith('sympy.') for k in sys.modules)
'''
    subprocess.run([sys.executable,'-c',script],check=True)


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
def test_missing_both_numeric_slots_keeps_original_error(monkeypatch,operation):
    value=QQ(7,3)
    monkeypatch.delattr(sx._core.DomainElement,'__float__')
    monkeypatch.delattr(sx._core.DomainElement,'__index__')
    result=compare_old(monkeypatch,operation,value)
    assert result[0]==('error','TypeError','must be real number, not symox._core.DomainElement')
