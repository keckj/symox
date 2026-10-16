"""Public factor views agree with SymPy without changing held expressions."""
from collections import defaultdict
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def expression(module, source):
    names = ('Mul', 'Add', 'Pow', 'Rational', 'Float', 'Symbol', 'exp', 'sin',
             'sqrt', 'UnevaluatedExpr', 'Limit', 'Function', 'oo', 'zoo', 'nan',
             'I', 'E', 'pi')
    namespace = {name: getattr(module, name) for name in names}
    namespace.update(zip(('x', 'y', 'z', 'a', 'b'), module.symbols('x y z a b')))
    return eval(source, {'__builtins__': {}}, namespace)


def converted(value):
    if isinstance(value, (list, tuple, set)):
        return type(value)(converted(item) for item in value)
    return sx.to_sympy(value)


@pytest.mark.parametrize('source', [
    'x', '0', '1', '-2', 'Rational(1,2)', 'Rational(-1,2)', 'I', '-I',
    'E', 'pi', 'oo', '-oo', 'zoo', 'nan', 'sin(x)', 'exp(x)',
    'x**a', 'Pow(Rational(1,2), x, evaluate=False)',
    'Pow(Rational(1,123456789123456789123456789), x, evaluate=False)',
    'Pow(Rational(2,3), x, evaluate=False)',
    'Pow(Float(0.5), x, evaluate=False)',
    'Pow(exp(x), 2, evaluate=False)', 'Function("exp")(x)',
])
def test_public_base_exponent_view_matches_oracle_without_mutating_args(source):
    value, reference = expression(sx, source), expression(sp, source)
    if isinstance(value, int):
        value, reference = sx.sympify(value), sp.sympify(reference)
    original = sx.srepr(value), value.args
    assert converted(value.as_base_exp()) == reference.as_base_exp()
    assert (sx.srepr(value), value.args) == original
    assert pickle.loads(pickle.dumps(value)).args == value.args


@pytest.mark.parametrize('source', [
    'x', '0', '1', '-2', 'Rational(-2,3)', 'Float(-2.5)',
    'I', '-I', 'oo', '-oo', 'zoo', 'nan', 'E', 'pi', 'sin(x)',
    'exp(x)', 'x**a', '2*x**a*y**b', '-2*I*x', 'x+y',
    'Mul(x,x,evaluate=False)', 'Mul(1,x,x,evaluate=False)',
    'Mul(x,1/x,evaluate=False)', 'Mul(x**a,x**b,evaluate=False)',
    'Mul(Mul(x,y,evaluate=False),x,evaluate=False)',
    'Mul(Pow(Rational(1,2),x,evaluate=False),2**y,evaluate=False)',
    'Mul(exp(x),exp(y),evaluate=False)',
    'Mul(-oo,x,evaluate=False)',
    'Mul(0,x,evaluate=False)', 'Function("exp")(x)',
])
def test_power_dictionary_keeps_held_factors_and_default_factory(source):
    value, reference = expression(sx, source), expression(sp, source)
    value, reference = sx.sympify(value), sp.sympify(reference)
    before = sx.srepr(value)
    actual, expected = value.as_powers_dict(), reference.as_powers_dict()
    assert isinstance(actual, type(expected))
    if isinstance(expected, defaultdict):
        assert isinstance(actual, defaultdict) and actual.default_factory is int
    else:
        assert type(actual) is dict
    assert {sx.to_sympy(base): sx.to_sympy(exponent) for base, exponent in actual.items()} == expected
    assert sx.srepr(value) == before
    if isinstance(expected, defaultdict):
        assert actual[sx.Symbol('absent')] == expected[sp.Symbol('absent')] == 0
    else:
        with pytest.raises(KeyError):
            actual[sx.Symbol('absent')]


@pytest.mark.parametrize('source', [
    'x', '-2', 'Rational(-2,3)', 'Float(-2.5)', '-1', 'oo', '-oo',
    'zoo', 'nan', 'I', 'E', 'x+y', '-2*x*y',
    'Mul(-2,x,y,evaluate=False)',
    'UnevaluatedExpr(x)', 'Limit(x,x,1)',
    'Mul(x,UnevaluatedExpr(y),z,evaluate=False)',
    'Mul(x,Limit(y,y,1),z,evaluate=False)',
])
@pytest.mark.parametrize('cset', [False, True])
@pytest.mark.parametrize('split', [False, True])
def test_commutative_partition_preserves_order_unknown_factors_and_split_sign(source, cset, split):
    value, reference = expression(sx, source), expression(sp, source)
    value, reference = sx.sympify(value), sp.sympify(reference)
    actual = value.args_cnc(cset=cset, split_1=split)
    assert isinstance(actual, list) and isinstance(actual[1], list)
    assert isinstance(actual[0], set if cset else list)
    assert converted(actual) == reference.args_cnc(cset=cset, split_1=split)


@pytest.mark.parametrize('warn', [True, False])
def test_commutative_sets_validate_repeated_held_factors(warn):
    x = sx.Symbol('x')
    held = sx.Mul(x, x, evaluate=False)
    assert held.args_cnc() == [[x, x], []]
    if warn:
        with pytest.raises(ValueError, match='repeated'):
            held.args_cnc(cset=True, warn=warn)
    else:
        assert held.args_cnc(cset=True, warn=warn) == [{x}, []]


@pytest.mark.parametrize('source', [
    'x', '-2', 'Float(-2.5)', 'Rational(-2,3)', 'I', 'E', 'oo', '-oo',
    'zoo', 'nan', 'x-y', 'y-x', '1-x', 'x-1', '2-x-y', '-2+x+y',
    '3-sqrt(2)', 'sqrt(2)-3', '-x-y+z', 'x-y+z', 'x*(y-x)',
    '-x*(y-x)', '-oo*x', 'zoo*x', '-2.5*x',
    'Add(x,-x,evaluate=False)', 'Add(-1,x,y,evaluate=False)',
    'Mul(-1,Add(x,-y,evaluate=False),evaluate=False)',
    'Mul(x,Add(y,-x,evaluate=False),evaluate=False)',
    'Symbol("negative", negative=True)',
])
def test_minus_sign_is_syntactic_and_uses_canonical_sum_tie_breaking(source):
    value, reference = expression(sx, source), expression(sp, source)
    value, reference = sx.sympify(value), sp.sympify(reference)
    assert value.could_extract_minus_sign() == reference.could_extract_minus_sign()


@pytest.mark.parametrize('source', ['oo', '-oo', 'nan', 'zoo', 'I', 'pi', 'Float(0.5)', 'Rational(1,2)'])
def test_number_class_includes_numeric_infinities_and_nan(source):
    value, reference = expression(sx, source), expression(sp, source)
    assert value.is_Number == reference.is_Number
    assert isinstance(value, sx.Number) == isinstance(reference, sp.Number)


@pytest.mark.xfail(strict=True, reason='Held arithmetic still differs in normalization and default_sort_key')
@pytest.mark.parametrize('source', [
    'Add(Mul(-1,2,zoo,evaluate=False),Add(-oo,oo,x,evaluate=False),evaluate=False)',
    'Add(oo,Mul(-1,sin(x),-1,evaluate=False),evaluate=False)',
    'Add(y,Add(y,-2,-1,evaluate=False),evaluate=False)',
    'Add(Mul(-1,x,Rational(1,2),evaluate=False),y,evaluate=False)',
])
def test_held_sign_ties_track_remaining_sympy_ordering_gaps(source):
    assert expression(sx, source).could_extract_minus_sign() == expression(sp, source).could_extract_minus_sign()


def test_assumptions_dummy_identity_and_held_zero_exponents_are_retained():
    a, b = sx.Dummy('x', positive=True), sx.Dummy('x', positive=True)
    value = sx.Mul(a**2, b, a**-2, evaluate=False)
    powers = value.as_powers_dict()
    assert powers == {a: 0, b: 1}
    assert len(powers) == 2 and a != b
    assert converted(value.args_cnc()) == value.to_sympy().args_cnc()


def test_factor_inspection_does_not_import_sympy():
    code = '''
import sys
class DenySympy:
    def find_spec(self, fullname, *args):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('unexpected SymPy dependency')
sys.meta_path.insert(0, DenySympy())
import symox as s
x=s.Symbol('x')
assert s.Pow(s.S.Half,x,evaluate=False).as_base_exp() == (2,-x)
assert s.I.as_base_exp() == (-1,s.S.Half)
assert s.Mul(x,1/x,evaluate=False).as_powers_dict() == {x:0}
assert (-2*x).args_cnc() == [[-1,2,x],[]]
assert (1-x).could_extract_minus_sign()
assert s.oo.is_Number and s.nan.is_Number and not s.zoo.is_Number
assert not any(name == 'sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable, '-c', code], check=True, capture_output=True, text=True)
