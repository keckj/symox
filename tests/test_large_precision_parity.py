"""Precision is determined by the request rather than a 10,000-digit policy cap."""
import pickle

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('precision',[1,40001,65536])
@pytest.mark.parametrize('literal',['1.234567890123456789','-0.1','1e-100'])
def test_explicit_binary_precision_matches_reference_rounding(precision,literal):
    actual=sx.Float(literal,precision=precision)
    expected=sp.Float(literal,precision=precision)
    assert actual._prec==expected._prec
    assert actual.to_sympy()._mpf_==expected._mpf_


@pytest.mark.parametrize('digits',[-3,-1,0,12001,20000])
def test_requested_decimal_precision_matches_reference_bits(digits):
    actual=sx.Float('1.234567890123456789',digits)
    expected=sp.Float('1.234567890123456789',digits)
    assert actual._prec==expected._prec
    assert actual.to_sympy()._mpf_==expected._mpf_
    value=sx.N(sx.Rational(1,3),n=digits)
    reference=sp.N(sp.Rational(1,3),n=digits)
    assert value._prec==reference._prec
    assert value.to_sympy()._mpf_==reference._mpf_


@pytest.mark.parametrize('precision',[1,40001,65536])
def test_numerical_fields_preserve_high_precision_and_pickle(precision):
    domain,reference=sx.RealField(precision),sp.RealField(precision)
    actual,expected=domain('0.1'),reference('0.1')
    assert domain.to_rational(actual,limit=False)==tuple(map(int,reference.to_rational(expected,limit=False)))
    assert pickle.loads(pickle.dumps(domain)).precision==precision
    assert pickle.loads(pickle.dumps(actual))==actual
    complex_field=sx.ComplexField(precision)
    z=complex_field('0.1','-0.2')
    assert z.parent().precision==precision
    assert pickle.loads(pickle.dumps(z))==z


def test_unrepresentable_working_precision_raises_without_wrapping_the_counter():
    for construct in [lambda:sx.Float('1.25',precision=2**32-1),
                      lambda:sx.Float('1.25',2**32-1),
                      lambda:sx.RealField(2**32-1),
                      lambda:sx.ComplexField(2**32-1)]:
        with pytest.raises(ValueError,match='working-precision'):
            construct()


def test_large_exact_numbers_cross_interoperability_without_decimal_integer_parsing():
    numerator=10**6000+1;denominator=10**6000+3
    assert sx.Integer(numerator).to_sympy()==sp.Integer(numerator)
    assert sx.Rational(numerator,denominator).to_sympy()==sp.Rational(numerator,denominator)
