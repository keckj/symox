"""Decimal and exact binary rationalization retain arbitrary precision."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('literal',['0','0.1','-0.125','0.333333333333333',
    '1.4142135623730951','3.141592653589793','1.2345678901234567',
    '-1234.56789','0.0001234567'])
@pytest.mark.parametrize('conversion',['base10','exact'])
def test_float_rationalization_modes_match_reference(literal,conversion):
    value=sx.Float(literal)
    actual=sx.nsimplify(value,rational=True,rational_conversion=conversion)
    expected=sp.nsimplify(value.to_sympy(),rational=True,rational_conversion=conversion)
    assert actual.is_Rational
    assert actual.to_sympy()==expected


@pytest.mark.parametrize('literal',['0.1','0.333333333333333','1.4142135623730951',
    '3.141592653589793','-1234.56789'])
@pytest.mark.parametrize('tolerance',[sx.Rational(1,10),sx.Rational(1,1000),1,10])
@pytest.mark.parametrize('conversion',['base10','exact'])
def test_explicit_rational_tolerance_matches_reference(literal,tolerance,conversion):
    value=sx.Float(literal);T=sx.sympify(tolerance).to_sympy()
    actual=sx.nsimplify(value,rational=True,tolerance=tolerance,rational_conversion=conversion)
    expected=sp.nsimplify(value.to_sympy(),rational=True,tolerance=T,rational_conversion=conversion)
    assert actual.to_sympy()==expected


@pytest.mark.parametrize('conversion',['base10','exact'])
def test_symbolic_coefficients_exponents_and_complex_parts(conversion):
    x=sx.Symbol('x')
    expression=sx.Float('.333333333333333')*x**sx.Float('.5')+sx.Float('.1')*sx.I+sx.Float('1.4142135623730951')
    result=sx.nsimplify(expression,rational=True,rational_conversion=conversion)
    expected=sp.nsimplify(expression.to_sympy(),rational=True,rational_conversion=conversion)
    assert result.to_sympy()==expected
    assert result.atoms(sx.Float)==set()
    assert result.free_symbols=={x}


def test_high_precision_and_extreme_exponents_do_not_underflow_to_zero():
    tiny=sx.Float('1e-1000',100);large=sx.Float('1e1000',100)
    assert sx.nsimplify(tiny,rational=True)==sx.Rational(1,10**1000)
    assert sx.nsimplify(large,rational=True)==sx.Integer(10**1000)
    third=sx.Float('0.'+'3'*90,100)
    assert sx.nsimplify(third,rational=True)==sx.Rational(1,3)
    for value in [tiny,large,sx.Float('1.234567890123456789123456789',80)]:
        result=sx.nsimplify(value,rational=True,rational_conversion='exact')
        assert result.to_sympy()==sp.Rational(value.to_sympy())
        assert result!=0


def test_rational_flag_selects_conversion_instead_of_irrational_recognition():
    value=sx.Float('1.4142135623730951')
    assert sx.nsimplify(value,rational=False)==sx.sqrt(2)
    assert sx.nsimplify(value,rational=True).is_Rational
    assert sx.nsimplify(sx.Float('.1'),rational=True,full=True)==sx.Rational(1,10)
    assert sx.nsimplify(sx.Float('.1'),rational=True,rational_conversion='exact')!=sx.Rational(1,10)
    for exact in [sx.pi,sx.sqrt(2),sx.Rational(1,3),sx.Symbol('x')]:
        assert sx.nsimplify(exact,rational=True)==exact


def test_invalid_options_and_large_conversion_exponents():
    value=sx.Float('.5')
    for tolerance in [0,-1,float('nan'),float('inf'),sx.I,sx.Symbol('epsilon')]:
        with pytest.raises(ValueError):sx.nsimplify(value,rational=True,tolerance=tolerance)
    with pytest.raises(ValueError):sx.nsimplify(value,rational=True,rational_conversion='unknown')
    with pytest.raises(TypeError):sx.nsimplify(value,rational='yes')
    with pytest.raises(NotImplementedError):sx.nsimplify(value,rational=False,tolerance=.01)
    with pytest.raises(NotImplementedError):sx.nsimplify(value,full=True)
    huge_exponent=sx.Float('1e-1000000',20)
    assert sx.nsimplify(huge_exponent,rational=True)==sx.Rational(1,10**1000000)
    exact=sx.nsimplify(huge_exponent,rational=True,rational_conversion='exact')
    assert exact.to_sympy()==sp.Rational(huge_exponent.to_sympy())


@pytest.mark.parametrize('precision',[100001,200003])
@pytest.mark.parametrize('conversion',['base10','exact'])
def test_rational_conversion_preserves_precision_beyond_former_cap(precision,conversion):
    value=sx.Float('0.125',precision=precision)
    assert sx.nsimplify(value,rational=True,rational_conversion=conversion)==sx.Rational(1,8)
    value=sx.Float('1.234567890123456789012345678901234567890123456789',precision=precision)
    result=sx.nsimplify(value,rational=True,rational_conversion=conversion)
    if conversion=='exact':assert result.to_sympy()==sp.Rational(value.to_sympy())
    else:assert sx.Float(result,precision=precision)==value


def test_rationalization_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
value=Float('.1')
assert nsimplify(value,rational=True)==Rational(1,10)
assert nsimplify(value,rational=True,rational_conversion='exact')==Rational(3602879701896397,36028797018963968)
assert nsimplify(Float('3.141592653589793'),rational=True,tolerance=Rational(1,1000))==Rational(355,113)
assert nsimplify(Float('1e-1000',100),rational=True)==Rational(1,10**1000)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
