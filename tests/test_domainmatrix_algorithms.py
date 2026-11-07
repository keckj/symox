"""Domain-sensitive linear algebra, checked against exact identities and SymPy."""
import pickle
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.polys.matrices import DM, DomainMatrix, DomainScalar
from symox.polys.matrices.exceptions import (DMBadInputError, DMDomainError,
    DMShapeError, DMNotAField, DMNonInvertibleMatrixError, DMNonSquareMatrixError)


@pytest.mark.parametrize('kind', [sx.Matrix, sx.ImmutableMatrix, sx.SparseMatrix,
                                  sx.ImmutableSparseMatrix])
@pytest.mark.parametrize('case', range(9))
def test_new_domain_inference(kind, case):
    x = sx.Symbol('x')
    precise = sx.Float('0.12345678901234567890123456789', 60)
    values = [sx.I, sx.I/2, precise, precise+sx.I, sx.sqrt(2),
              sx.sin(x), sx.exp(x), sx.sqrt(2)+x, sx.Float('0.0')]
    domains = [sx.ZZ_I, sx.QQ_I, sx.RealField(precise._prec),
               sx.ComplexField(precise._prec), sx.EXRAW, sx.ZZ[sx.sin(x)],
               sx.ZZ[sx.exp(x)], sx.EXRAW, sx.RR]
    matrix = kind([[values[case], 0], [0, 1]])
    actual = matrix.to_DM()
    assert actual.domain == domains[case]
    for left, right in zip(actual.to_Matrix(), matrix):
        assert sp.simplify(left.to_sympy()-right.to_sympy()) == 0
    assert pickle.loads(pickle.dumps(actual)) == actual


DOMAINS = [sx.ZZ, sx.QQ, sx.GF(2), sx.GF(5), sx.ZZ_I, sx.QQ_I,
           sx.QQ.algebraic_field(sx.sqrt(2)), sx.EX, sx.EXRAW]


def same(a, b):
    # The reference ``DomainMatrix.__eq__`` also compares the storage format.
    return a.to_sparse() == b.to_sparse()


def reference(matrix):
    from sympy.polys.matrices import DomainMatrix as Reference
    reference_domains = [sp.ZZ, sp.QQ, sp.GF(2), sp.GF(5), sp.ZZ_I, sp.QQ_I,
                         sp.QQ.algebraic_field(sp.sqrt(2)), sp.EX, sp.EXRAW]
    return Reference.from_Matrix(matrix.to_Matrix().to_sympy()).convert_to(
        reference_domains[DOMAINS.index(matrix.domain)])


@pytest.mark.parametrize('domain', DOMAINS)
@pytest.mark.parametrize('fmt', ['dense', 'sparse'])
def test_charpoly_determinant_adjugate(domain, fmt):
    generator = random.Random(426)
    for size in range(6):
        rows = [[generator.randrange(-5, 6) for _ in range(size)] for _ in range(size)]
        actual = DomainMatrix(rows, (size, size), domain, fmt=fmt)
        expected = reference(actual)
        assert expected.domain.from_sympy(domain.to_sympy(actual.det()).to_sympy()) == expected.det()
        assert [expected.domain.from_sympy(domain.to_sympy(v).to_sympy())
                for v in actual.charpoly()] == expected.charpoly()
        adjugate, determinant = actual.adj_det()
        assert same(actual*adjugate, DomainMatrix.eye(size, domain)*determinant)
        assert same(adjugate*actual, actual*adjugate)
        assert determinant == actual.det()
        # Cayley-Hamilton also tests small-characteristic arithmetic.
        polynomial = DomainMatrix.zeros((size, size), domain)
        identity = DomainMatrix.eye(size, domain)
        for coefficient in actual.charpoly():
            polynomial = polynomial*actual + identity*coefficient
        assert same(polynomial, DomainMatrix.zeros((size, size), domain))


@pytest.mark.parametrize('domain', DOMAINS)
@pytest.mark.parametrize('method', ['auto', 'GJ', 'FF'])
def test_rectangular_rref_and_kernel(domain, method):
    generator = random.Random(119)
    for m, n in [(0, 0), (0, 3), (3, 0), (1, 4), (4, 1), (3, 5), (5, 3)]:
        rows = [[generator.randrange(-2, 3) for _ in range(n)] for _ in range(m)]
        if m:
            rows[0] = [0]*n
        actual = DomainMatrix(rows, (m, n), domain)
        expected, expected_pivots = reference(actual).rref()
        reduced, pivots = actual.rref(method=method)
        assert pivots == expected_pivots
        assert reduced.to_Matrix().to_sympy() == expected.to_Matrix()
        numerator, denominator, den_pivots = actual.rref_den(method=method)
        assert den_pivots == pivots
        assert same(numerator.to_field()/denominator, reduced)
        assert numerator.domain == domain
        kernel = actual.nullspace()
        assert kernel.shape == (n-len(pivots), n)
        assert same(actual*kernel.transpose(), DomainMatrix.zeros((m, n-len(pivots)), domain))
        assert kernel.rank() == n-len(pivots)
        if domain.is_Field:
            normalized = actual.nullspace(divide_last=True)
            assert same(actual*normalized.transpose(), DomainMatrix.zeros((m, normalized.shape[0]), domain))
            rowspace, colspace = actual.rowspace(), actual.columnspace()
            assert rowspace.rank() == colspace.rank() == len(pivots)


@pytest.mark.parametrize('domain', DOMAINS)
@pytest.mark.parametrize('method', [None, 'rref', 'charpoly'])
def test_solve_and_inverse_denominator(domain, method):
    actual = DM([[0, 1, 0], [1, 2, 3], [0, 0, 1]], domain)
    rhs = DM([[1, 2], [3, 4], [5, 6]], domain)
    numerator, denominator = actual.solve_den(rhs, method=method)
    assert same(actual*numerator, rhs*denominator)
    inverse, divisor = actual.inv_den(method=method)
    assert same(actual*inverse, DomainMatrix.eye(3, domain)*divisor)
    if domain.is_Field:
        assert same(actual.inv()*actual, DomainMatrix.eye(3, domain))
        assert same(actual**-2, actual.inv()*actual.inv())
    else:
        with pytest.raises(DMNotAField):
            actual.inv()


@pytest.mark.parametrize('domain', [sx.ZZ.poly_ring(sx.Symbol('x')), sx.QQ_I,
    sx.QQ.algebraic_field(sx.sqrt(2)), sx.EX])
def test_nonrational_coefficients(domain):
    x = sx.Symbol('x')
    value = x if domain.is_Composite else sx.I if domain == sx.QQ_I else (
        sx.sqrt(2) if domain.is_AlgebraicField else sx.sin(x))
    matrix = DM([[value, 1], [0, value+1]], domain)
    # SymPy 1.14 keeps a converted EX element as given and normalizes the
    # arithmetic, so the determinant equals the product of the elements.
    assert matrix.det() == domain(value)*domain(value+1)
    adjugate, determinant = matrix.adj_det()
    assert same(matrix*adjugate, DomainMatrix.eye(2, domain)*determinant)
    inverse, denominator = matrix.inv_den()
    # Formal EX equality need not prove all equivalent rational trees.
    difference = (matrix*inverse - DomainMatrix.eye(2, domain)*denominator).to_Matrix()
    assert all(sx.simplify(entry) == 0 for entry in difference)
    assert matrix.rref()[1] == (0, 1)


@pytest.mark.parametrize('bits', [53, 120, 240])
@pytest.mark.parametrize('complex_', [False, True])
def test_numeric_arithmetic_precision(bits, complex_):
    domain = sx.ComplexField(bits) if complex_ else sx.RealField(bits)
    value = domain('0.123456789012345678901234567890123456789')
    if complex_:
        value += domain(sx.I)
    actual = DM([[value, 1], [2, 3]], domain)
    inverse = actual.inv()
    assert inverse.domain == domain
    identity = (actual*inverse).to_list()
    for i in range(2):
        for j in range(2):
            assert domain.almosteq(identity[i][j], int(i == j))
    assert domain.almosteq(actual.det(), 3*value-2)
    assert all(v.parent() == domain for v in (actual+actual).to_list_flat())


def test_views_unification_and_scalar_results():
    a = DM([[1, 2, 0], [3, 4, 5]], sx.ZZ)
    assert a[-1, -1].element == sx.ZZ(5)
    assert isinstance(a[0, 0], DomainScalar)
    assert a[::-1, ::-1].to_Matrix() == sx.Matrix([[5, 4, 3], [0, 2, 1]])
    assert a.extract([1, 1, 0], [2, 0]).to_Matrix() == sx.Matrix([[5, 3], [5, 3], [0, 1]])
    for index in [(2, 0), (0, 3), (-3, 0)]:
        with pytest.raises(IndexError):
            a[index]
    assert a.hstack(a).shape == (2, 6)
    assert a.vstack(a).shape == (4, 3)
    assert same(a.transpose().transpose(), a)
    assert same(sx.ZZ(2)*a, a.mul(2))
    rational = a.convert_to(sx.QQ)
    assert (a+rational).domain == sx.QQ
    assert (a*rational.transpose()).domain == sx.QQ
    assert same(DM([[1, 2], [3, 4]], sx.ZZ)**0, DomainMatrix.eye(2, sx.ZZ))
    q = DM([[sx.QQ(1, 2), sx.QQ(1, 3)]], sx.QQ)
    denominator, numerator = q.clear_denoms(convert=True)
    assert denominator.domain == numerator.domain == sx.ZZ
    assert denominator.element == sx.ZZ(6)
    assert denominator*q == numerator.convert_to(sx.QQ)
    denominator2, numerator2 = q.clear_denoms()
    assert numerator2.domain == sx.QQ and denominator2 == denominator
    assert numerator2 == denominator*q
    assert a[0, 0] + DomainScalar(sx.QQ(1, 2), sx.QQ) == DomainScalar(sx.QQ(3, 2), sx.QQ)


@pytest.mark.parametrize('kind',[sx.Matrix,sx.ImmutableMatrix,sx.SparseMatrix,sx.ImmutableSparseMatrix])
@pytest.mark.parametrize('domain',[sx.ZZ,sx.QQ,sx.QQ_I,sx.RR,sx.CC,sx.EX])
def test_typed_coefficient_on_the_left_of_explicit_matrices(kind,domain):
    scalar=domain.convert(sx.I if domain in (sx.CC,sx.QQ_I) else 2)
    matrix=kind([[1,2],[0,3]])
    assert scalar*matrix == matrix*scalar.as_expr()


def test_singular_inconsistent_and_invalid_inputs():
    matrix = DM([[1, 2], [2, 4]], sx.QQ)
    for method in [None, 'rref', 'charpoly']:
        with pytest.raises(DMNonInvertibleMatrixError):
            matrix.inv_den(method=method)
    with pytest.raises(DMNonInvertibleMatrixError):
        matrix.inv()
    with pytest.raises(DMNonInvertibleMatrixError):
        DM([[1], [2]], sx.QQ).solve_den(DM([[1], [3]], sx.QQ))
    rectangular = DM([[1, 2]], sx.QQ)
    for method in [rectangular.det, rectangular.charpoly, rectangular.inv]:
        with pytest.raises(DMNonSquareMatrixError):
            method()
    # The reference rejects an underdetermined system by its shape.
    with pytest.raises(DMShapeError):
        rectangular.solve_den(DM([[1]], sx.QQ))
    # The reference rejects an unknown rref method with ValueError.
    with pytest.raises(ValueError):
        matrix.rref(method='unimplemented')
    with pytest.raises(DMDomainError):
        matrix.add(matrix.convert_to(sx.ZZ))
    # The reference DDM.hstack checks the row counts with an assert.
    with pytest.raises(AssertionError):
        matrix.hstack(rectangular)
    for rows in [[], [[0, 0]], [[0, 0], [0, 0]]]:
        zero = DomainMatrix.from_list(rows, sx.ZZ)
        assert zero.rank() == 0


def test_empty_and_overdetermined_systems():
    for domain in [sx.ZZ, sx.QQ, sx.GF(5)]:
        empty = DomainMatrix.zeros((0, 0), domain)
        numerator, denominator = empty.inv_den()
        assert numerator == empty and denominator == domain.one
        a = DM([[0, 0], [1, 2], [0, 1], [2, 4]], domain)
        solution = DM([[3], [4]], domain)
        numerator, denominator = a.solve_den(a*solution)
        assert numerator == solution*denominator


def test_no_sympy_import():
    source = '''
import sys
class Block:
    def find_spec(self, fullname, *args):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0, Block())
from symox import ZZ, QQ, GF, I, Matrix, Float
from symox.polys.matrices import DM, DomainMatrix
for K in [ZZ, QQ, GF(5)]:
    A = DM([[1,2],[3,5]], K)
    N, d = A.inv_den()
    assert A*N == (DomainMatrix.eye(2,K)*d).to_dense()
    assert A.rref()[1] == (0,1)
assert Matrix([[I]]).to_DM().domain.is_ZZ_I
assert Matrix([[Float('0.1',60)]]).to_DM().domain.precision > 190
'''
    subprocess.run([sys.executable, '-c', source], check=True)
