"""Actual subclass reconstruction and result-class precedence."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('operation', ['add', 'sub', 'product', 'matmul', 'scale',
    'rscale', 'divide', 'power', 'negate', 'transpose', 'T', 'H', 'C', 'slice',
    'extract', 'reshape', 'row', 'col', 'applyfunc', 'subs', 'xreplace',
    'copy', 'inverse', 'rref', 'cholesky'])
def test_matrix_operations_preserve_the_real_python_subclass(operation):
    def run(lib):
        class Custom(lib.Matrix):
            def label(self):
                return 'custom matrix'
        x = lib.Symbol('x', real=True)
        matrix = Custom([[x+3, 1], [1, 3]])
        functions = dict(add=lambda:matrix+matrix, sub=lambda:matrix-matrix,
            product=lambda:matrix*matrix, matmul=lambda:matrix@matrix,
            scale=lambda:matrix*2, rscale=lambda:2*matrix, divide=lambda:matrix/2,
            power=lambda:matrix**2, negate=lambda:-matrix,
            transpose=lambda:matrix.transpose(), T=lambda:matrix.T,
            H=lambda:matrix.H, C=lambda:matrix.C,
            slice=lambda:matrix[:,:1], extract=lambda:matrix.extract([1],[0,1]),
            reshape=lambda:matrix.reshape(1,4), row=lambda:matrix.row(0),
            col=lambda:matrix.col(0), applyfunc=lambda:matrix.applyfunc(lambda v:v+1),
            subs=lambda:matrix.subs(x,2), xreplace=lambda:matrix.xreplace({x:lib.Integer(2)}),
            diff=lambda:matrix.diff(x), copy=lambda:matrix.copy(),
            inverse=lambda:matrix.inv(), rref=lambda:matrix.rref()[0],
            cholesky=lambda:matrix.subs(x,2).cholesky())
        result = functions[operation]()
        assert type(result) is Custom
        assert result.label() == 'custom matrix'
        return result
    actual, expected = run(sx), run(sp)
    assert actual.shape == expected.shape
    assert all(sp.simplify(a-b)==0 for a,b in zip(actual.to_sympy(), expected))


def test_reconstruction_bypasses_a_public_constructor_with_different_arguments():
    for lib in (sx,sp):
        class Custom(lib.Matrix):
            def __new__(cls, value):
                return super().__new__(cls, [[value]])
        matrix = Custom(2)
        for result in [matrix+matrix, matrix.T, matrix.inv(), matrix._new(1,1,[7])]:
            assert type(result) is Custom
        assert (matrix+matrix)[0,0] == 4


@pytest.mark.parametrize('operation', ['add','sub','product'])
@pytest.mark.parametrize('priorities', [(4,4),(8,4),(4,8)])
def test_arithmetic_uses_the_reference_result_class_precedence(operation, priorities):
    for lib in (sx,sp):
        class Left(lib.Matrix):
            _class_priority=priorities[0]
        class Right(lib.Matrix):
            _class_priority=priorities[1]
        left,right=Left([[2]]),Right([[3]])
        if lib is sx:
            from symox.matrices.matrixbase import classof
        else:
            from sympy.matrices.matrixbase import classof
        assert classof(left,right) is (Left if priorities[0]>priorities[1] else Right)
        result=dict(add=lambda:left+right,sub=lambda:left-right,
                    product=lambda:left*right)[operation]()
        left_wins = priorities[0]>=priorities[1] if operation in {'add','sub'} else priorities[0]>priorities[1]
        assert type(result) is (Left if left_wins else Right)


def test_immutable_peer_and_explicit_conversions_select_their_target_category():
    for lib in (sx,sp):
        class Custom(lib.Matrix):pass
        matrix=Custom([[2]])
        immutable=lib.ImmutableMatrix([[3]])
        assert isinstance(matrix+immutable,lib.ImmutableMatrix)
        assert isinstance(matrix.as_immutable(),lib.ImmutableMatrix)
        assert type(matrix.as_mutable()) is lib.Matrix
        with pytest.raises(TypeError):
            Custom._new([[2]],copy=False)


def _class_label(lib, cls, other, value):
    if isinstance(value, (list, tuple)):
        return tuple(_class_label(lib, cls, other, item) for item in value)
    if type(value) is cls:
        return 'subclass'
    if type(value) is other:
        return 'peer'
    if isinstance(value, lib.ImmutableMatrix):
        return 'ImmutableMatrix'
    if isinstance(value, lib.MatrixBase):
        return 'Matrix'
    return type(value).__name__


_RESULT_CASES = {
    'col_join subclass base': lambda m: m.a.col_join(m.B),
    'col_join base subclass': lambda m: m.B.col_join(m.a),
    'col_join subclass immutable': lambda m: m.a.col_join(m.I),
    'col_join immutable subclass': lambda m: m.I.col_join(m.a),
    'col_join priority': lambda m: m.a.col_join(m.m2),
    'col_join priority reversed': lambda m: m.m2.col_join(m.a),
    'row_join subclass base': lambda m: m.a.row_join(m.B),
    'row_join base subclass': lambda m: m.B.row_join(m.a),
    'hstack base class': lambda m: m.lib.Matrix.hstack(m.a, m.B),
    'hstack subclass on base': lambda m: m.M.hstack(m.B, m.a),
    'hstack immutable first': lambda m: m.M.hstack(m.I, m.a),
    'hstack priority': lambda m: m.lib.Matrix.hstack(m.a, m.m2),
    'hstack empty': lambda m: m.M.hstack(),
    'vstack subclass on base': lambda m: m.M.vstack(m.B, m.a),
    'instance zeros': lambda m: m.a.zeros(2),
    'diag cls': lambda m: m.M.diag(1, 2, cls=m.lib.Matrix),
    'diag immutable cls': lambda m: m.M.diag(1, 2, cls=m.lib.ImmutableMatrix),
    'elementwise subclass base': lambda m: m.a.multiply_elementwise(m.B),
    'elementwise base subclass': lambda m: m.B.multiply_elementwise(m.a),
    'elementwise subclass immutable': lambda m: m.a.multiply_elementwise(m.I),
    'elementwise immutable base': lambda m: m.I.multiply_elementwise(m.B),
    'elementwise immutable subclass': lambda m: m.I.multiply_elementwise(m.a),
    'cross subclass base': lambda m: m.v3.cross(m.b3),
    'cross base subclass': lambda m: m.b3.cross(m.v3),
    'cross immutable subclass': lambda m: m.lib.ImmutableMatrix([0, 1, 0]).cross(m.v3),
    'cross subclass immutable': lambda m: m.v3.cross(m.lib.ImmutableMatrix([0, 1, 0])),
    'LUsolve subclass base': lambda m: m.a.LUsolve(m.v),
    'LUsolve base subclass': lambda m: m.B.LUsolve(m.M([1, 2])),
    'LUsolve subclass immutable': lambda m: m.a.LUsolve(m.lib.ImmutableMatrix([1, 2])),
    'LUsolve immutable subclass': lambda m: m.I.LUsolve(m.M([1, 2])),
    'LUsolve immutable base': lambda m: m.I.LUsolve(m.v),
    'solve base subclass': lambda m: m.B.solve(m.M([1, 2])),
    'solve subclass base': lambda m: m.a.solve(m.v),
    'solve subclass immutable': lambda m: m.a.solve(m.lib.ImmutableMatrix([1, 2])),
    'solve immutable subclass': lambda m: m.I.solve(m.M([1, 2])),
    'gauss_jordan_solve base subclass': lambda m: m.B.gauss_jordan_solve(m.M([1, 2])),
    'col_insert base subclass': lambda m: m.B.col_insert(0, m.M([1, 2])),
    'col_insert subclass base': lambda m: m.a.col_insert(0, m.v),
    'col_insert subclass immutable': lambda m: m.a.col_insert(0, m.lib.ImmutableMatrix([1, 2])),
    'row_insert priority': lambda m: m.a.row_insert(0, m.M2([[1, 2]])),
    'jacobian': lambda m: m.M([m.x, m.x**2]).jacobian(m.lib.Matrix([m.x])),
    'LDLsolve base subclass': lambda m: m.B.LDLsolve(m.M([1, 2])),
    'lower_triangular_solve base subclass': lambda m: m.B.lower_triangular_solve(m.M([1, 2])),
    'cholesky_solve base subclass': lambda m: m.B.cholesky_solve(m.M([1, 2])),
    'add priority': lambda m: m.a + m.m2,
    'add priority reversed': lambda m: m.m2 + m.a,
    'mul priority': lambda m: m.a * m.m2,
    'mul priority reversed': lambda m: m.m2 * m.a,
    'diag function': lambda m: m.lib.diag(m.a, m.B),
    'jordan_block': lambda m: m.M.jordan_block(2, 3),
    'companion': lambda m: m.M.companion(m.lib.Poly(m.x**2+1, m.x)),
    'wilkinson': lambda m: m.M.wilkinson(1),
    'orthogonalize': lambda m: m.M.orthogonalize(m.M([1, 2]), m.M([0, 1])),
    'irregular': lambda m: m.M.irregular(2, m.a, m.B),
    'evalf': lambda m: m.a.evalf(),
    'adjugate': lambda m: m.a.adjugate(),
    'cofactor_matrix': lambda m: m.a.cofactor_matrix(),
    'exp': lambda m: m.M([[1, 0], [0, 4]]).exp(),
    'diagonal': lambda m: m.a.diagonal(),
    'upper_triangular': lambda m: m.a.upper_triangular(),
    'minor_submatrix': lambda m: m.a.minor_submatrix(0, 0),
    'normalized': lambda m: m.M([1, 2]).normalized(),
    'pinv': lambda m: m.M([[1, 0], [0, 4]]).pinv(),
    'rot90': lambda m: m.a.rot90(),
    'permute': lambda m: m.a.permute([[0, 1]]),
    'vec': lambda m: m.a.vec(),
    'vech': lambda m: m.M([[1, 0], [0, 4]]).vech(),
    'replace': lambda m: m.a.replace(m.x, 2),
    'integrate': lambda m: m.a.integrate(m.x),
    'limit': lambda m: m.a.limit(m.x, 0),
    'as_real_imag': lambda m: m.a.as_real_imag(),
    'nullspace': lambda m: m.M([[1, 1], [1, 1]]).nullspace(),
    'columnspace': lambda m: m.a.columnspace(),
    'diagonalize': lambda m: m.M([[1, 0], [0, 4]]).diagonalize(),
    'LUdecomposition': lambda m: m.a.LUdecomposition()[:2],
    'QRdecomposition': lambda m: m.M([[1, 0], [0, 4]]).QRdecomposition(),
    'jordan_form': lambda m: m.M([[1, 0], [0, 4]]).jordan_form(),
    'inverse_LU': lambda m: m.a.inverse_LU(),
    'echelon_form': lambda m: m.a.echelon_form(),
    'LDLdecomposition': lambda m: m.M([[1, 0], [0, 4]]).LDLdecomposition(),
    'cholesky': lambda m: m.M([[1, 0], [0, 4]]).cholesky(),
    'analytic_func': lambda m: m.M([[1, 0], [0, 4]]).analytic_func(m.lib.exp(m.x), m.x),
    'kronecker_product': lambda m: m.lib.kronecker_product(m.a, m.B),
    'elementary_row_op': lambda m: m.a.elementary_row_op('n->kn', row=0, k=2),
    'singular_value_decomposition': lambda m: m.M([[1, 0], [0, 4]]).singular_value_decomposition(),
    'rank_decomposition': lambda m: m.M([[1, 0], [0, 4]]).rank_decomposition(),
    'project': lambda m: m.M([1, 2]).project(m.M([1, 1])),
    'abs': lambda m: abs(m.a),
    'mod': lambda m: m.M([[3, 5]]) % 2,
    'evalf subs': lambda m: m.a.evalf(subs={m.x: 1}),
    'rmul integer': lambda m: 3 * m.a,
    'diff function': lambda m: m.lib.diff(m.a, m.x),
    'simplify function': lambda m: m.lib.simplify(m.a),
    'getitem lists': lambda m: m.a[[0, 1], [0]],
    'extract': lambda m: m.a.extract([0], [1]),
}


class _Fixture:
    def __init__(self, lib):
        class M(lib.Matrix):
            pass

        class M2(lib.Matrix):
            _class_priority = 5

        self.lib, self.M, self.M2, self.x = lib, M, M2, lib.Symbol('x')
        self.a = M([[2, self.x], [1, 3]])
        self.B = lib.Matrix([[1, 0], [0, 4]])
        self.I = lib.ImmutableMatrix([[1, 0], [0, 4]])
        self.m2 = M2([[1, 0], [0, 4]])
        self.v = lib.Matrix([1, 2])
        self.v3 = M([1, 0, 2])
        self.b3 = lib.Matrix([0, 1, 0])


@pytest.mark.parametrize('name', list(_RESULT_CASES))
def test_subclass_result_classes_match_sympy(name):
    labels = []
    for lib in (sx, sp):
        fixture = _Fixture(lib)
        result = _RESULT_CASES[name](fixture)
        labels.append(_class_label(lib, fixture.M, fixture.M2, result))
    assert labels[0] == labels[1]


# SymPy warns about the rejected nested-list form before raising.
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_matrix_constructor_accepts_reference_options():
    x = sx.Symbol('x')
    for lib in (sx, sp):
        source = [[1, 2], [3, 4]]
        assert lib.Matrix(2, 2, [1, 2, 3, 4], copy=False) == lib.Matrix(source)
        with pytest.raises(TypeError):
            lib.Matrix(source, copy=False)
        assert lib.Matrix(2, 2, [1, 2, 3, 4], evaluate=True) == lib.Matrix(source)
        assert lib.ImmutableMatrix(source, copy=True) == lib.ImmutableMatrix(source)
        assert lib.SparseMatrix(source, copy=False) == lib.SparseMatrix(source)
    assert sx.Matrix([[x + x]], evaluate=False)[0, 0] == 2*x


def test_subclass_iteration_and_numeric_operators():
    for lib in (sx, sp):
        class M(lib.Matrix):
            pass
        matrix = M([[-3, 5], [7, -1]])
        assert [int(value) for value in matrix] == [-3, 5, 7, -1]
        assert list(abs(matrix)) == [3, 5, 7, 1]
        assert list(matrix % 4) == [1, 1, 3, 3]
        assert type(abs(matrix)) is M and type(matrix % 4) is M
