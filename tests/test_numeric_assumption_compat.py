"""Exact facts and unknown answers for numerical scalar values."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
from sympy.core.assumptions import _assume_defined
import symox as sx


SOURCES = [
    'Integer(0)', 'Integer(1)', 'Integer(-1)', 'Integer(2)', 'Integer(4)',
    'Integer(97)', 'Integer(18446744073709551557)', '-Integer(2)**100',
    'Rational(1,2)', 'Rational(-17,8)',
    'Float(0)', 'Float(2)', 'Float(-2)', 'Float(2.5)', 'Float(-2.5)',
    'Float("0",80)', 'Float("2",80)', 'Float("-2",80)',
    'Float("2.5",80)', 'Float("-2.5",80)',
    'Float("1e100",80)', 'Float("-1e100",80)',
    'Float("1e-100",80)', 'Float("-1e-100",80)',
    'I', 'oo', '-oo', 'zoo', 'nan', 'sqrt(2)', 'pi', 'E', 'EulerGamma',
    'Float(2)+Float(2.5)', 'Rational(1,2)+Rational(1,3)',
    'oo+pi', '1-oo', 'nan*sqrt(2)', 'sqrt(2)+1', 'sqrt(2)-1',
]


@pytest.mark.parametrize('source', SOURCES)
def test_all_standard_numeric_facts_match_sympy(source):
    actual, expected = sx.sympify(source), sp.sympify(source)
    for name in sorted(_assume_defined):
        assert getattr(actual, 'is_'+name) is getattr(expected, 'is_'+name), (source, name)


@pytest.mark.parametrize('source', ['Float(0)', 'Float(2)', 'Float(2.5)',
                                   'Float("2",80)', 'Float("2.5",80)'])
def test_approximate_integer_queries_preserve_type_and_roundtrips(source):
    value = sx.sympify(source)
    before = sx.srepr(value)
    assert value.is_Float and not value.is_Integer and not value.is_Rational
    for _ in range(4):
        assert value.is_integer is sp.sympify(source).is_integer
    assert sx.srepr(value) == before
    assert sx.from_sympy(sx.to_sympy(value)) == value
    for protocol in range(pickle.HIGHEST_PROTOCOL+1):
        restored = pickle.loads(pickle.dumps(value, protocol=protocol))
        assert restored.is_Float and not restored.is_Integer
        assert restored.is_integer is value.is_integer
        assert restored.is_even is value.is_even and restored.is_odd is value.is_odd


@pytest.mark.parametrize('source', ['Rational(1,2)', 'I', 'sqrt(2)', 'pi', 'oo', '-oo', 'zoo'])
def test_proven_nonintegers_are_neither_even_nor_odd(source):
    value = sx.sympify(source)
    assert value.is_integer is False
    assert value.is_even is False and value.is_odd is False


def test_nan_unknown_facts_do_not_change_exceptional_arithmetic():
    for name in _assume_defined - {'commutative'}:
        assert getattr(sx.nan, 'is_'+name) is None
    assert sx.nan.is_commutative is True
    assert sx.nan+1 == sx.nan
    assert sx.nan*0 == sx.nan
    assert 0*sx.oo == sx.nan and sx.oo-sx.oo == sx.nan
    assert sx.sign(sx.oo) == 1 and sx.sign(-sx.oo) == -1
    assert sx.exp(-sx.oo) == 0 and sx.exp(sx.oo) == sx.oo


def test_numeric_queries_without_sympy_runtime():
    code = r'''
import sys
class NoSymPy:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise ModuleNotFoundError('SymPy disabled')
sys.meta_path.insert(0, NoSymPy())
import symox as s
assert s.Float(0).is_integer is True
assert s.Float(0).is_even is True and s.Float(0).is_odd is False
assert s.Float(0).is_rational is None and s.Float(0).is_irrational is None
assert s.Float(2).is_integer is None
assert s.Float('2.5',80).is_integer is False
assert s.Float('2.5',80).is_prime is False
assert s.Float('2.5',80).is_composite is False
assert s.Float('2.5',80).is_noninteger is True
assert s.sqrt(2).is_even is False and s.sqrt(2).is_odd is False
assert (s.sqrt(2)-1).is_extended_positive is True
assert s.oo.is_positive is False and s.oo.is_extended_positive is True
assert (-s.oo).is_negative is False and (-s.oo).is_extended_negative is True
assert s.nan.is_finite is None and s.nan.is_real is None
assert s.EulerGamma.is_integer is None
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
