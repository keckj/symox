"""Differential checks of graph factors and exact domain conversions."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.exceptions import DMBadInputError, DMDomainError, DMShapeError


KINDS = [sx.Matrix, sx.ImmutableMatrix, sx.SparseMatrix, sx.ImmutableSparseMatrix]


@pytest.mark.parametrize('cls', KINDS)
@pytest.mark.parametrize('method,options', [
    ('connected_components_decomposition', {}),
    ('strongly_connected_components_decomposition', {}),
    ('strongly_connected_components_decomposition', {'lower': False}),
])
@pytest.mark.parametrize('rows', [
    [[1, 0, 2], [0, 3, 0], [4, 0, 5]],
    [[1, 0, 2, 0], [3, 4, 0, 0], [0, 0, 5, 0], [0, 6, 7, 8]],
    [[0, 0], [0, 0]], [[1]], [],
])
def test_graph_decomposition(cls, method, options, rows):
    matrix = cls(rows)
    if not rows and method == 'strongly_connected_components_decomposition':
        with pytest.raises(ValueError):
            getattr(matrix, method)(**options)
        with pytest.raises(ValueError):
            getattr(sp.Matrix(rows), method)(**options)
        return
    p, blocks = getattr(matrix, method)(**options)
    expected_p, expected_blocks = getattr(sp.Matrix(rows), method)(**options)
    assert p.as_explicit().to_sympy() == expected_p.as_explicit()
    assert blocks.as_explicit().to_sympy() == expected_blocks.as_explicit()
    explicit_p = p.as_explicit()
    assert explicit_p.T * blocks.as_explicit() * explicit_p == matrix
    assert (p.T * blocks * p).as_explicit() == matrix
    assert explicit_p.T * explicit_p == sx.eye(matrix.rows)
    if method == 'connected_components_decomposition':
        assert isinstance(blocks, sx.BlockDiagMatrix)
    else:
        assert isinstance(blocks, sx.BlockMatrix)


@pytest.mark.parametrize('cls', KINDS)
def test_graph_non_square_and_symbolic(cls):
    for method in ('connected_components_decomposition', 'strongly_connected_components_decomposition'):
        with pytest.raises(sx.NonSquareMatrixError):
            getattr(cls([[1, 2]]), method)()
        x = sx.Symbol('x')
        matrix = cls([[0, 0, x], [0, 1, 0], [x, 0, 0]])
        p, blocks = getattr(matrix, method)()
        assert p.as_explicit().T * blocks.as_explicit() * p.as_explicit() == matrix


@pytest.mark.parametrize('cls', KINDS)
@pytest.mark.parametrize('case', range(8))
@pytest.mark.parametrize('field', [False, True])
def test_domain_inference_and_roundtrip(cls, case, field):
    x, y = sx.symbols('x y')
    examples = [[], [[0, 0]], [[1, 2], [3, 4]], [[sx.Rational(1, 2), 1]],
                [[x, y]], [[x/2, y]], [[1/x, y/(x+1)]], [[1/(2*x), x]]]
    matrix = cls(examples[case])
    dm = matrix.to_DM(field=field)
    expected = sp.Matrix(matrix.to_sympy()).to_DM(field=field)
    assert str(dm.domain) == str(expected.domain)
    assert dm.shape == expected.shape
    assert dm.to_Matrix() == matrix
    assert all(v.parent() == dm.domain for v in dm.to_list_flat())
    assert pickle.loads(pickle.dumps(dm)) == dm


@pytest.mark.parametrize('cls', KINDS)
def test_explicit_domains_and_conversion_failures(cls):
    x = sx.Symbol('x')
    matrix = cls([[1, 2], [3, 4]])
    for domain in (sx.QQ, sx.ZZ.poly_ring(x), sx.QQ.algebraic_field(sx.sqrt(2))):
        dm = matrix.to_DM(domain=domain)
        assert dm.domain == domain
        assert dm.to_Matrix() == matrix
    assert matrix.to_DM(domain=sx.GF(5)).to_Matrix() == sx.Matrix([[1, 2], [-2, -1]])
    assert matrix.to_DM(domain=sx.GF(3)).to_Matrix() == sx.Matrix([[1, -1], [0, 1]])
    with pytest.raises(TypeError):
        matrix.to_DM(domain=sx.QQ, field=True)
    with pytest.raises(sx.CoercionFailed):
        cls([[sx.Rational(1, 2)]]).to_DM(domain=sx.ZZ)
    formal = cls([[sx.sqrt(2)*x]]).to_DM()
    assert formal.domain == sx.EXRAW and formal.to_Matrix() == cls([[sx.sqrt(2)*x]])
    assert cls([[sx.sin(x)]]).to_DM().domain == sx.ZZ[sx.sin(x)]
    numerical = cls([[sx.Float('1.2')]]).to_DM()
    assert numerical.domain == sx.RR and numerical.to_Matrix() == cls([[sx.Float('1.2')]])
    with pytest.raises(sx.OptionError):
        matrix.to_DM(unsupported=True)
    algebraic = cls([[sx.sqrt(2), 1]])
    dm = algebraic.to_DM(extension=True)
    assert dm.domain.is_AlgebraicField and dm.to_Matrix() == algebraic


def test_domain_matrix_owns_entries_and_converts():
    matrix = sx.Matrix([[1, 2], [3, 4]])
    dm = matrix.to_DM()
    matrix[0, 0] = 99
    data = dm.to_dod()
    data[0][0] = sx.ZZ(99)
    assert dm.to_Matrix()[0, 0] == 1
    assert dm.copy() == dm
    assert dm.convert_to(sx.QQ).to_Matrix() == dm.to_Matrix()
    assert dm.to_field().domain == sx.QQ
    assert dm.to_field().choose_domain().domain == sx.ZZ
    assert dm.to_dense().to_list() == dm.to_sparse().to_list()
    assert dm.matmul(dm).to_Matrix() == dm.to_Matrix()**2
    assert dm.add(dm).to_Matrix() == 2*dm.to_Matrix()
    assert dm.det() == sx.ZZ(-2)
    assert DomainMatrix(dm.to_list(), dm.shape, dm.domain).to_sparse() == dm
    with pytest.raises(DMDomainError):
        dm.add(dm.to_field())
    with pytest.raises(DMShapeError):
        dm.matmul(sx.zeros(1, 2).to_DM())
    with pytest.raises(DMBadInputError):
        DomainMatrix([[1]], (1, 2), sx.ZZ)
    with pytest.raises(DMBadInputError):
        DomainMatrix({0: {2: 1}}, (1, 2), sx.ZZ)


def test_without_sympy():
    subprocess.run([sys.executable, '-c', '''
import sys
class BlockSymPy:
    def find_spec(self, name, *args):
        if name == 'sympy' or name.startswith('sympy.'):
            raise RuntimeError('SymPy imported at runtime')
sys.meta_path.insert(0, BlockSymPy())
import symox as s
x = s.Symbol('x')
m = s.Matrix([[x, 0, 2], [0, 1, 0], [3, 0, 4]])
for method in ('connected_components_decomposition', 'strongly_connected_components_decomposition'):
    p, b = getattr(m, method)()
    assert p.as_explicit().T * b.as_explicit() * p.as_explicit() == m
assert m.to_DM().to_Matrix() == m
assert m.to_DM(field=True).domain == s.ZZ.frac_field(x)
assert s.factorial(x).is_hypergeometric(x) is True
assert s.sin(x).is_hypergeometric(x) is False
a = s.Wild('a')
assert (a*x).matches(2*x) == {a: s.Integer(2)}
assert (~(x & s.Symbol('y'))).to_nnf() == s.to_nnf(~(x & s.Symbol('y')))
assert s.sin(x).dir(x, -1) == -1
assert s.Limit(s.sin(x)/x, x, 0).dir == '+'
'''], check=True)
