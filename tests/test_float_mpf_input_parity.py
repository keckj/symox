"""Float accepts exact mpf payloads and mpmath objects without decimal detours."""
import subprocess
import sys

import mpmath
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('precision',[1,3,24,53,167])
@pytest.mark.parametrize('payload',[
    (0,3,-2),(1,3,-2),(0,3,-2,2),(1,3,-2,2),
    (1,'0xff',-3,8),(0,'ffL',-3,8),
    (0,0,0),(0,0,0,0),(0,'0x0',0,0),
    (0,0,-456,-2),(1,0,-789,-3),(0,0,-123,-1),
    (0,1,10**12,1),(1,1,-10**12,1),
])
def test_mpf_tuple_shapes_precision_and_special_values(payload,precision):
    actual=sx.Float(payload,precision=precision).to_sympy()
    expected=sp.Float(payload,precision=precision)
    assert type(actual) is type(expected)
    if expected is sp.nan:assert actual is sp.nan
    else:assert actual==expected


@pytest.mark.parametrize('digits',[15,80])
@pytest.mark.parametrize('options',[{},{'precision':167},{'dps':20}])
def test_mpmath_inputs_follow_reference_precision(digits,options):
    with mpmath.workdps(digits):
        for literal in ['.1','-.125','0','+inf','-inf','nan']:
            value=mpmath.mpf(literal)
            actual=sx.Float(value,**options).to_sympy();expected=sp.Float(value,**options)
            assert type(actual) is type(expected)
            if expected is sp.nan:assert actual is sp.nan
            else:assert actual==expected
            if expected.is_Float:assert actual._prec==expected._prec


def test_large_binary_mantissa_is_rounded_without_integer_strings():
    integer=10**6000+123
    payload=(0,integer,-20000,integer.bit_length())
    actual=sx.Float(payload,precision=167)
    assert actual.to_sympy()==sp.Float(payload,precision=167)


def test_mpf_input_handling_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
import mpmath
assert Float((0,3,-2),30)==Float('.75',30)
assert Float(mpmath.mpf('.75'),30)==Float('.75',30)
assert Float((0,0,-456,-2))==oo
assert Float((0,1,10**12,1),precision=53).is_Float
'''
    subprocess.run([sys.executable,'-c',source],check=True)
