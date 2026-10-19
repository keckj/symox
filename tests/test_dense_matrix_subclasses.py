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
