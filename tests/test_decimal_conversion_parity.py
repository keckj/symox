"""Decimal parsing and printing keep certified rounding at large exponents."""
import random
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('precision',[1,2,3,24,53,167])
@pytest.mark.parametrize('exponent',[-10**12,-1000001,-1000000,-5000,-130,0,130,5000,1000000,1000001,10**12])
def test_large_exponent_float_values_and_printing_match_reference(precision,exponent):
    for coefficient in ['1','-1.234567890123456789','9.999999999999999999']:
        literal=f'{coefficient}e{exponent}'
        actual=sx.Float(literal,precision=precision)
        expected=sp.Float(literal,precision=precision)
        assert actual.to_sympy()==expected
        if precision>=24:assert str(actual)==str(expected)


def test_random_decimal_rounding_and_binary_roundtrip():
    rng=random.Random(917)
    for _ in range(250):
        precision=rng.choice([1,2,3,24,53,107,257,1000])
        digits=''.join(str(rng.randrange(10)) for _ in range(rng.randrange(1,80)))
        exponent=rng.choice([-10000,-500,-3,0,3,500,10000])
        literal=f"{'-' if rng.randrange(2) else ''}0.{digits}e{exponent}"
        actual=sx.Float(literal,precision=precision);expected=sp.Float(literal,precision=precision)
        assert actual.to_sympy()==expected
        if precision>=24 and not actual.is_zero:assert str(actual)==str(expected)


@pytest.mark.parametrize('literal',['1.015625','-1.015625','1.234375','-1.234375'])
def test_decimal_display_uses_reference_half_away_convention(literal):
    actual=sx.Float(literal,precision=24);expected=sp.Float(literal,precision=24)
    assert str(actual)==str(expected)


@pytest.mark.parametrize('literal',['1e2000000000000000000','-1e-2000000000000000000'])
def test_binary_representable_exponents_beyond_decimal_validation_range(literal):
    for options in [{},{'dps':20},{'precision':167}]:
        actual=sx.Float(literal,**options);expected=sp.Float(literal,**options)
        assert actual.to_sympy()==expected
        assert actual._prec==expected._prec
    for module in [sx,sp]:
        with pytest.raises(ValueError):module.Float(literal,'')


@pytest.mark.parametrize('sign',[-1,1])
def test_literal_evaluation_and_numerical_fields_do_not_apply_transcendental_magnitude_caps(sign):
    source=f'1e{sign*10**12}'
    value=sx.Float(source,30);reference=sp.Float(source,30)
    for digits in [5,50]:
        assert sx.N(value,digits).to_sympy()==sp.N(reference,digits)
        complex_value=value+sx.Float('.125',30)*sx.I
        expected=reference+sp.Float('.125',30)*sp.I
        assert sx.N(complex_value,digits).to_sympy()==sp.N(expected,digits)
    for factory in [sx.RealField,sx.ComplexField]:
        field=factory(prec=167)
        element=field.convert(value)
        assert element==element
        assert field.to_sympy(element).to_sympy()==sp.Float(reference,precision=167)


def test_compact_literals_have_explicit_materialized_ratio_limits():
    value=sx.Float('1e1000000000000',20)
    with pytest.raises(NotImplementedError,match='resource budget'):
        sx.nsimplify(value,rational=True)
    with pytest.raises(sx.CoercionFailed,match='resource budget'):
        sx.nsimplify(value,rational=True,rational_conversion='exact')


@pytest.mark.parametrize('source',[lambda:10**6000+123,lambda:sx.Integer(10**6000+123)])
def test_float_integer_inputs_do_not_pass_through_decimal_integer_strings(source):
    value=source()
    for precision in [53,257]:
        actual=sx.Float(value,precision=precision)
        integer=10**6000+123
        expected=sp.Float((0,integer,0,integer.bit_length()),precision=precision)
        assert actual.to_sympy()==expected
    automatic=sx.Float(value)
    assert automatic==sx.Float(10**6000+123,6001)


def test_large_decimal_exponent_conversion_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
for sign in [-1,1]:
    value=Float('1e'+str(sign*10**12),30)
    assert str(value).startswith('1.00000000000000000000000000000e')
    assert N(value,50)._prec==Float(1,50)._prec
    field=RealField(prec=167)
    assert field.convert(value)==field.convert(value)
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=15)
