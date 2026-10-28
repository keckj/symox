"""Native domain proofs avoid adapters without changing matrix semantics."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('kind', [
    'Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix',
])
@pytest.mark.parametrize('source', [
    '[[Rational(1, 3), 2], [3, Rational(7, 5)]]',
    '[[x**2 + y/3, x*y], [y + 1, (x + y)**3]]',
    '[[x, 1, y], [0, x + y, 2], [1, y**2, x/7]]',
])
def test_polynomial_charpoly_stays_native(kind, source, monkeypatch):
    def matrix(module):
        x, y = module.symbols('x y')
        return getattr(module, kind)(eval(source, {
            'x': x, 'y': y, 'Rational': module.Rational,
        }))

    original = matrix(sx)
    expected = matrix(sp).charpoly()

    def unexpected(*args, **kwargs):
        raise AssertionError('polynomial entries unnecessarily entered domain inference')

    from symox.polys import constructor
    monkeypatch.setattr(constructor, 'construct_domain', unexpected)
    actual = original.charpoly()
    assert str(actual.domain) == str(expected.domain)
    assert [sp.expand(c.to_sympy() - r) for c, r in
            zip(actual.all_coeffs(), expected.all_coeffs())] == [0] * len(expected.all_coeffs())
    assert original == matrix(sx)


@pytest.mark.parametrize('kind', [
    'Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix',
])
@pytest.mark.parametrize('method', [None, 'DM'])
def test_rational_inverse_bypasses_numeric_adapters(kind, method, monkeypatch):
    from symox import _matrix_inverse as adapter
    data = [[sx.Rational(1, 3), 2], [3, sx.Rational(7, 5)]]
    matrix = getattr(sx, kind)(data)
    expected = sp.Matrix([[sp.Rational(1, 3), 2], [3, sp.Rational(7, 5)]]).inv(method=method)

    def unexpected(*args, **kwargs):
        raise AssertionError('exact rational inverse entered an inexact adapter')

    monkeypatch.setattr(adapter, '_inexact_domain', unexpected)
    monkeypatch.setattr(adapter, '_inexact_composite_DM', unexpected)
    actual = matrix.inv(method=method)
    assert type(actual) is type(matrix)
    assert matrix * actual == sx.eye(2)
    assert sp.Matrix([[v.to_sympy() for v in row] for row in actual.tolist()]) == expected


def test_rational_inverse_still_reports_shape_and_singularity():
    with pytest.raises(sx.NonSquareMatrixError):
        sx.zeros(2, 3).inv()
    with pytest.raises(sx.NonInvertibleMatrixError):
        sx.Matrix([[1, 2], [2, 4]]).inv()
    assert sx.zeros(0).inv() == sx.zeros(0)


def test_exraw_charpoly_still_simplifies_coefficients():
    x = sx.Symbol('x')
    matrix = sx.Matrix([[sx.sin(x)**2 + sx.cos(x)**2, 1], [1, 1]])
    assert not matrix._has_rational_polynomial_entries()
    assert matrix.charpoly().all_coeffs() == [1, -2, 0]


def test_noncommutative_symbols_do_not_prove_a_commutative_domain():
    a = sx.Symbol('a', commutative=False)
    assert not sx.Matrix([[a, 0], [0, 1]])._has_rational_polynomial_entries()


def test_inexact_composite_domain_and_precision_are_preserved():
    x = sx.Symbol('x')
    matrix = sx.Matrix([[sx.Float('1.234567890123456789', 45)*x, 1], [0, x + 2]])
    assert not matrix._has_rational_polynomial_entries()
    actual = matrix.charpoly()
    X = sp.Symbol('x')
    expected = sp.Matrix([[sp.Float('1.234567890123456789', 45)*X, 1], [0, X + 2]]).charpoly()
    assert str(actual.domain) == str(expected.domain)
    for value, ref in zip(actual.all_coeffs(), expected.all_coeffs()):
        assert value.to_sympy() == ref


def test_generator_identity_and_custom_coefficient_callback():
    x, y = sx.symbols('x y', real=True)
    matrix = sx.Matrix([[x, 1], [0, 2]])
    assert matrix.charpoly(y).gen == y
    assert matrix.charpoly(x).gen != x
    seen = []
    actual = matrix.charpoly(y, simplify=lambda c: seen.append(c) or c.subs(x, 0))
    assert len(seen) == 3
    assert actual.all_coeffs() == [1, -2, 0]


def test_polynomial_charpoly_respects_subclass_domain_override():
    calls = []

    class CustomMatrix(sx.MutableDenseMatrix):
        def to_DM(self, *args, **kwargs):
            calls.append(self)
            return super().to_DM(*args, **kwargs)

    x = sx.Symbol('x')
    result = CustomMatrix([[x, 1], [0, 2]]).charpoly()
    assert len(calls) == 1
    assert result.all_coeffs() == [1, -x - 2, 2*x]


def test_exact_matrix_kernels_do_not_load_mpmath():
    source = '''
import sys
import symox as s
assert 'mpmath' not in sys.modules
class Block:
    def find_spec(self, name, *args):
        if name == 'mpmath' or name.startswith('mpmath.'):
            raise AssertionError('exact native kernel imported mpmath')
sys.meta_path.insert(0, Block())
x,y = s.symbols('x y')
assert s.Matrix([[x+y, 1], [0, x]]).charpoly().all_coeffs() == [1, -2*x-y, x*x+x*y]
assert s.Matrix([[1,2],[3,4]]).inv() == s.Matrix([[-2,1],[s.Rational(3,2),-s.Rational(1,2)]])
assert 'mpmath' not in sys.modules
'''
    subprocess.run([sys.executable, '-c', source], check=True, timeout=30)
