"""Public constructor classes must expose operations on their native instances."""
import pytest
import symox as sx


@pytest.mark.parametrize('constructor',[sx.ImmutableMatrix,sx.ImmutableDenseMatrix,sx.MatrixBase])
def test_explicit_constructor_unbound_operations(constructor):
    matrix=sx.ImmutableMatrix([[1,2],[3,4]])
    assert constructor.det(matrix)==-2
    assert constructor.inv(matrix)==matrix.inv()
    assert constructor.applyfunc(matrix,f=lambda value:value+1)==matrix.applyfunc(lambda value:value+1)
    assert isinstance(constructor.applyfunc(matrix,f=lambda value:value+1),sx.ImmutableMatrix)


@pytest.mark.parametrize('constructor',[sx.MatrixExpr,sx.MatrixSymbol,sx.MatAdd,sx.MatMul])
def test_symbolic_constructor_unbound_operations(constructor):
    matrix=sx.MatrixSymbol('A',2,2)
    replacement=sx.MatrixSymbol('B',2,2)
    assert constructor.subs(matrix,matrix,replacement)==replacement
    assert constructor.det(matrix)==matrix.det()
    assert constructor.has(matrix,matrix)


def test_delegation_preserves_factory_mutability():
    assert isinstance(sx.ImmutableMatrix.zeros(2),sx.ImmutableMatrix)
    assert isinstance(sx.ImmutableMatrix.eye(2),sx.ImmutableMatrix)
    assert isinstance(sx.ImmutableMatrix.diag(1,2),sx.ImmutableMatrix)
