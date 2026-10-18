"""Algebraic field metadata survives structural expression protocols."""
from collections import defaultdict
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def examples():
    return [sx.AlgebraicNumber(sx.sqrt(2)),
            sx.AlgebraicNumber(-sx.sqrt(2)),
            sx.AlgebraicNumber(sx.sqrt(2), [2, 3], alias='a'),
            sx.AlgebraicNumber(sx.I),
            sx.AlgebraicNumber(sx.sqrt(2)/2),
            sx.AlgebraicNumber(sx.sqrt(2), [sx.Rational(1, 2), 0]),
            sx.AlgebraicNumber(2**sx.Rational(1, 3)+sx.I)]


def convert(value):
    if hasattr(value, 'to_sympy'):
        return value.to_sympy()
    if isinstance(value, dict):
        return {convert(k): convert(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return type(value)(convert(v) for v in value)
    return value


@pytest.mark.parametrize('index', range(7))
@pytest.mark.parametrize('method', ['as_base_exp', 'as_coeff_Add', 'as_coeff_Mul',
    'as_coeff_add', 'as_coeff_mul', 'args_cnc', 'as_coefficients_dict',
    'as_content_primitive', 'as_dummy', 'as_independent', 'as_numer_denom',
    'as_ordered_factors', 'as_ordered_terms', 'as_powers_dict', 'copy',
    'could_extract_minus_sign', 'getO', 'getn', 'removeO', 'doit',
    'transpose', 'primitive', 'class_key', 'sort_key', 'count_ops'])
def test_structural_views_preserve_algebraic_number(index, method):
    a = examples()[index]
    actual, expected = getattr(a, method)(), getattr(a.to_sympy(), method)()
    assert convert(actual) == expected


@pytest.mark.parametrize('index', range(7))
def test_atoms_queries_and_identity(index):
    a = examples()[index]
    reference = a.to_sympy()
    assert convert(a.atoms()) == reference.atoms()
    for actual_type, expected_type in [(sx.Symbol, sp.Symbol), (sx.Integer, sp.Integer),
                                        (sx.AlgebraicNumber, sp.AlgebraicNumber), (sx.Pow, sp.Pow)]:
        assert convert(a.atoms(actual_type)) == reference.atoms(expected_type)
        assert a.count(actual_type) == reference.count(expected_type)
        assert convert(a.find(actual_type)) == reference.find(expected_type)
    needles = [sx.Integer(0), sx.Integer(1), sx.Integer(2), sx.sqrt(2), sx.Symbol('a')]
    for needle in needles:
        assert a.has(needle) == reference.has(needle.to_sympy())
        assert a.has_free(needle) == reference.has_free(needle.to_sympy())
        assert a.has_xfree({needle}) == reference.has_xfree({needle.to_sympy()})
    assert convert(a.find(lambda node: node.is_Number)) == reference.find(lambda node: node.is_Number)
    assert a.is_same(a) and a.dummy_eq(a)
    assert not a.is_same(a.as_expr())
    assert a.equals(a.as_expr()) is True
    assert pickle.loads(pickle.dumps(a)) == a
    assert sx.AlgebraicNumber.fromiter(a.args) == a
    assert a.fromiter(a.args) == a


@pytest.mark.parametrize('index', range(7))
def test_scalar_constant_calculus_and_simplifiers(index):
    a = examples()[index]
    x = sx.Symbol('x')
    for method in ['series', 'nseries', 'aseries', 'as_leading_term', 'fps', 'fourier_series']:
        assert getattr(a, method)(x) is a
    assert list(a.lseries(x)) == [a]
    assert a.leadterm(x) == (a, sx.Integer(0))
    assert a.taylor_term(0, x) is a
    assert a.taylor_term(3, x) == 0
    assert a.diff(x) == 0
    assert sp.simplify(a.integrate(x).to_sympy()-a.to_sympy().as_expr()*x.to_sympy()) == 0
    assert a.limit(x, 0) is a
    for method in ['apart', 'cancel', 'combsimp', 'factor', 'gammasimp',
                   'normal', 'nsimplify', 'powsimp', 'radsimp', 'ratsimp',
                   'rcall', 'refine', 'simplify', 'together', 'trigsimp', 'expand']:
        actual = getattr(a, method)()
        assert isinstance(actual, sx.AlgebraicNumber)
        assert actual == a
    assert a.is_constant(x) and a.is_polynomial(x)
    assert a.is_rational_function(x) and a.is_algebraic_expr(x)
    assert a.is_hypergeometric(x) and a.is_meromorphic(x, 0)


@pytest.mark.parametrize('index', range(7))
def test_conjugation_and_real_imaginary_parts(index):
    a = examples()[index]
    result = a.conjugate()
    assert isinstance(result, sx.AlgebraicNumber)
    assert sp.simplify(result.as_expr().to_sympy()-sp.conjugate(a.as_expr().to_sympy())) == 0
    assert a.adjoint() == result
    re, im = a.as_real_imag()
    def scalar(v):
        return v.as_expr() if isinstance(v, sx.AlgebraicNumber) else v
    assert (scalar(re)+sx.I*scalar(im)-a.as_expr()).expand().simplify() == 0


@pytest.mark.parametrize('index', range(7))
def test_structural_substitution_and_replacement(index):
    a = examples()[index]
    b = sx.Symbol('b')
    reference = a.to_sympy()
    replacements = [(sx.sqrt(2), sx.sqrt(3)), (sx.Integer(2), sx.Integer(3)),
                    (sx.Symbol('a'), b)]
    for old, new in replacements:
        actual = a.subs(old, new)
        assert convert(actual) == reference.subs(old.to_sympy(), new.to_sympy())
        actual = a.xreplace({old: new})
        assert convert(actual) == reference.xreplace({old.to_sympy(): new.to_sympy()})
    assert a.xreplace({a: b}) == b
    assert a.subs(a, b) == b
    assert a.replace(a, b) == b
    assert a.replace(lambda node: node == a, lambda node: b) == b
    changed = a.replace(sx.Integer(2), sx.Integer(3))
    assert convert(changed) == reference.replace(sp.Integer(2), sp.Integer(3))


def test_unordered_and_simultaneous_substitution():
    a = sx.AlgebraicNumber(sx.sqrt(2), [2, 3], alias='a')
    replacements = {sx.Integer(2): sx.Integer(3), sx.sqrt(2): sx.sqrt(5)}
    expected = a.to_sympy().subs(convert(replacements))
    assert convert(a.subs(replacements)) == expected
    # SymPy 1.14 tries to store temporary symbolic masks in rational field
    # coordinates here and raises. Check the simultaneous definition directly.
    result = a.subs(replacements, simultaneous=True)
    assert result.root == sx.sqrt(5)
    assert result.coeffs() == [sx.Integer(3), sx.Integer(3)]
    assert result.alias == a.alias
    x, y = sx.symbols('x y')
    a = sx.AlgebraicNumber(sx.sqrt(2), alias='x')
    assert a.subs({x: y, y: x}, simultaneous=True).alias == y


@pytest.mark.parametrize('index', range(7))
def test_integral_generator_conversion(index):
    a = examples()[index]
    result = a.to_algebraic_integer()
    assert isinstance(result, sx.AlgebraicNumber)
    assert result.minpoly.LC() == 1
    expected = a.to_sympy().to_algebraic_integer()
    assert convert(result) == expected
    assert result.coeffs() == a.coeffs()


def test_matching_comparison_and_alias_dependency():
    values = examples()
    for a in values:
        assert a.matches(a) == {}
        assert a.match(a) == {}
        wild = sx.Wild('w')
        assert a.match(wild) == {wild: a}
        assert a.as_coeff_exponent(a) == (sx.Integer(1), sx.Integer(1))
        for b in values:
            assert a.compare(b) == a.to_sympy().compare(b.to_sympy())
    wild = sx.Wild('w')
    pattern = sx.AlgebraicNumber(sx.sqrt(2), alias=wild)
    assert pattern.matches(sx.AlgebraicNumber(sx.sqrt(2), alias='a')) == {wild: sx.Symbol('a')}
    a = values[2]
    alias = a.alias
    for name in ['as_coeff_add', 'as_coeff_mul', 'as_independent',
                 'is_polynomial', 'is_rational_function', 'is_algebraic_expr']:
        assert convert(getattr(a, name)(alias)) == getattr(a.to_sympy(), name)(alias.to_sympy())


def test_numeric_rounding_and_native_runtime():
    a = sx.AlgebraicNumber(sx.sqrt(2))
    assert a.round(5).to_sympy() == a.to_sympy().round(5)
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x'); a=s.AlgebraicNumber(s.sqrt(2),[2,3],alias='a')
assert a.diff(x)==0
assert a.simplify()==a
assert a.subs(s.Symbol('a'),x).alias==x
assert a.count_ops()==2
assert s.AlgebraicNumber(s.sqrt(2)/2).to_algebraic_integer().minpoly.LC()==1
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
