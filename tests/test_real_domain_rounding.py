"""Directed rounding certifies the stored dyadic, not the input's decimal text."""
from fractions import Fraction
import math
import subprocess
import sys

import pytest
import symox as sx
from test_domain_element_rounding import compare_old


def field_value(mantissa,exponent,precision):
    value=sx.Float((int(mantissa<0),abs(mantissa),exponent,abs(mantissa).bit_length()),precision=precision)
    return sx.RealField(precision).convert(value)


def stored(value):
    expression=value.parent().to_sympy(value)
    if expression.is_Float:
        sign,mantissa,exponent,_=expression._mpf_
        return (-1 if sign else 1)*mantissa,exponent
    assert expression.is_Integer
    return int(expression),0


def certify(value):
    mantissa,exponent=stored(value)
    literal=Fraction(mantissa)*Fraction(2)**exponent
    lower,upper=math.floor(value),math.ceil(value)
    assert type(lower) is int and type(upper) is int
    assert lower<=literal<lower+1
    assert upper-1<literal<=upper
    assert lower==literal.numerator//literal.denominator
    assert upper==-(-literal.numerator//literal.denominator)
    return lower,upper


@pytest.mark.parametrize('precision',[53,80,200,500])
@pytest.mark.parametrize('mantissa,exponent',[(0,0),(1,-1),(-1,-1),(3,-1),(-3,-1),(12,-2),(-12,-2),
    ((1<<100)+1,-3),(-((1<<100)+1),-3),((1<<160)+1,-80),((1<<160)-1,-80),
    (-((1<<160)+1),-80),(-((1<<160)-1),-80),(1,-1074),(-1,-1074),(1,-1075),(-1,-1075)])
def test_rounding_certifies_the_literal_at_each_declared_precision(precision,mantissa,exponent):
    certify(field_value(mantissa,exponent,precision))


def test_the_concrete_200_bit_witness_uses_its_stored_dyadic():
    value=sx.RealField(200).convert(sx.QQ((1<<100)+1,3),sx.QQ)
    mantissa,exponent=stored(value)
    assert (mantissa,exponent)==(1071292029505993517027974728228286835414954148789526359004501,-101)
    assert certify(value)==(422550200076076467165567735125,422550200076076467165567735126)


@pytest.mark.parametrize('precision',[53,80,200,500])
@pytest.mark.parametrize('sign',[-1,1])
def test_feasible_large_integer_exponents_are_materialized_exactly(precision,sign):
    value=field_value(sign*3,100_000,precision)
    assert math.floor(value)==math.ceil(value)==sign*(3<<100_000)


@pytest.mark.parametrize('precision',[53,200])
@pytest.mark.parametrize('sign',[-1,1])
@pytest.mark.parametrize('exponent',[-100_000,-(1<<190)])
def test_sparse_tiny_values_need_no_giant_denominator(precision,sign,exponent):
    value=field_value(sign*3,exponent,precision)
    mantissa,actual_exponent=stored(value)
    # 0<abs(m)<2**bit_length(m), and bit_length(m)+e<=0, prove
    # -1<value<0 or 0<value<1 without materializing 2**(-e).
    assert mantissa and abs(mantissa).bit_length()+actual_exponent<=0
    assert math.floor(value)==(-1 if mantissa<0 else 0)
    assert math.ceil(value)==(0 if mantissa<0 else 1)


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
@pytest.mark.parametrize('sign',[-1,1])
def test_positive_exponents_beyond_existing_capacity_report_overflow(operation,sign):
    value=field_value(sign, (1<<31)+1,200)
    with pytest.raises(OverflowError,match='too many digits in integer'):operation(value)


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
@pytest.mark.parametrize('kind',['count','wrong','interrupt','nan','inf','neginf'])
def test_real_field_float_override_keeps_original_callbacks(monkeypatch,operation,kind):
    value=field_value(3,100_000,200);calls=[]
    def convert(self):
        calls.append('float')
        if kind=='count':return -2.75
        if kind=='wrong':return 3
        if kind=='interrupt':raise KeyboardInterrupt('RR float hook')
        return {'nan':float('nan'),'inf':float('inf'),'neginf':-float('inf')}[kind]
    monkeypatch.setattr(sx._core.DomainElement,'__float__',convert)
    compare_old(monkeypatch,operation,value,calls)


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
def test_custom_float_precedes_native_materialization_capacity(monkeypatch,operation):
    value=field_value(1,(1<<190),200);calls=[]
    def convert(self):calls.append('float');return 2.75
    monkeypatch.setattr(sx._core.DomainElement,'__float__',convert)
    assert compare_old(monkeypatch,operation,value,calls)[0]==('value',2 if operation is math.floor else 3)
    assert calls==['float']


@pytest.mark.parametrize('operation',[math.floor,math.ceil])
@pytest.mark.parametrize('domain',['complex','expression'])
def test_other_numeric_domains_keep_the_complete_old_path(monkeypatch,operation,domain):
    literal=sx.Float((0,(1<<180)+1,-80,181),precision=200)
    value=(sx.ComplexField(200) if domain=='complex' else sx.EX).convert(literal)
    compare_old(monkeypatch,operation,value)


def test_no_reference_runtime_for_wide_sparse_and_large_integer_exponents():
    script=r'''
import sys,importlib.abc,math
class NoReference(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('runtime reference')
sys.meta_path.insert(0,NoReference())
import symox as s
K=s.RealField(200)
v=K.convert(s.QQ((1<<100)+1,3),s.QQ)
assert math.floor(v)==422550200076076467165567735125
assert math.ceil(v)==422550200076076467165567735126
for negative in [0,1]:
 tiny=K.convert(s.Float((negative,1,-(1<<190),1),precision=200))
 assert math.floor(tiny)==(-1 if negative else 0)
 assert math.ceil(tiny)==(0 if negative else 1)
large=K.convert(s.Float((0,3,100000,2),precision=200))
assert math.floor(large)==math.ceil(large)==3<<100000
assert not any(k=='sympy' or k.startswith('sympy.') for k in sys.modules)
'''
    subprocess.run([sys.executable,'-c',script],check=True)
