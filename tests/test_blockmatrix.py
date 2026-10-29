"""True matrix-valued block grids, partition validation and safe block rewrites."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def matrices(module, prefix=''):
    A=module.MatrixSymbol(prefix+'A',2,2);B=module.MatrixSymbol(prefix+'B',2,1)
    C=module.MatrixSymbol(prefix+'C',1,2);D=module.MatrixSymbol(prefix+'D',1,1)
    return A,B,C,D


def block(module, prefix=''):
    A,B,C,D=matrices(module,prefix)
    return module.BlockMatrix([[A,B],[C,D]])


def same(actual, expected):
    actual=actual.to_sympy()
    if not isinstance(actual,sp.MatrixBase):actual=actual.as_explicit()
    if not isinstance(expected,sp.MatrixBase):expected=expected.as_explicit()
    assert actual.shape==expected.shape
    assert all(sp.simplify(a-b)==0 for a,b in zip(actual,expected))


def test_constructor_interface_interchange_and_immutability():
    A,B,C,D=matrices(sx);value=sx.BlockMatrix([[A,B],[C,D]])
    assert value.shape==(3,3) and value.blockshape==(2,2)
    assert value.rowblocksizes==[2,1] and value.colblocksizes==[2,1]
    assert isinstance(value,sx.BlockMatrix) and isinstance(value,sx.MatrixExpr)
    assert isinstance(sx.BlockDiagMatrix(A,D),sx.BlockMatrix)
    assert value.args==(value.blocks,)
    assert value.func(*value.args)==value
    assert value.blocks[0,1]==B
    assert value.is_structurally_symmetric
    assert value.structurally_equal(block(sx,'other'))
    assert not value.structurally_equal(sx.BlockDiagMatrix(sx.MatrixSymbol('E',1,1),sx.MatrixSymbol('F',2,2)))
    assert sx.sympify(sx.srepr(value))==value
    assert sx.MatrixExpr(value.to_sympy())==value
    same(value,block(sp))
    assert value.free_symbols=={A,B,C,D}
    cell=sx.Matrix([[1,2]])
    explicit=sx.BlockMatrix([[cell,sx.ImmutableMatrix([[3]])]])
    cell[0,0]=7
    same(explicit,sp.ImmutableMatrix([[1,2,3]]))
    row=sx.BlockMatrix([A,B]);assert row.shape==(2,3) and row.blockshape==(1,2)
    assert sx.BlockMatrix(sx.BlockDiagMatrix(A,D).blocks).blocks[1,1]==D


@pytest.mark.parametrize('factory',[
    lambda:sx.BlockMatrix([]), lambda:sx.BlockMatrix([[]]),
    lambda:sx.BlockMatrix([[sx.MatrixSymbol('A',2,2)],[sx.MatrixSymbol('B',1,1),sx.MatrixSymbol('C',1,1)]]),
    lambda:sx.BlockMatrix([[sx.MatrixSymbol('A',2,2),sx.MatrixSymbol('B',1,1)]]),
    lambda:sx.BlockMatrix([[sx.MatrixSymbol('A',2,2),sx.MatrixSymbol('B',2,1)],[sx.MatrixSymbol('C',1,1),sx.MatrixSymbol('D',1,2)]]),
    lambda:sx.BlockMatrix([[sx.MatrixSymbol('A',2,2),0]]),
])
def test_invalid_block_partitions_and_scalar_cells_are_rejected(factory):
    with pytest.raises((ValueError,TypeError)):
        factory()


def test_symbolic_partition_equality_must_be_provable():
    n,m=sx.symbols('n m',integer=True,positive=True)
    A=sx.MatrixSymbol('A',n,n);B=sx.MatrixSymbol('B',m,n)
    with pytest.raises(ValueError):sx.BlockMatrix([[A,B]])
    value=sx.BlockMatrix([[A,sx.MatrixSymbol('C',n,m)]])
    assert value.shape==(n,n+m)
    assert value.free_symbols=={A,sx.MatrixSymbol('C',n,m)}
    assert value.xreplace({n:2,m:3}).shape==(2,5)


@pytest.mark.parametrize('operation',[
    lambda s,X,Y:s.block_collapse(X+Y),
    lambda s,X,Y:s.block_collapse(X*Y),
    lambda s,X,Y:s.block_collapse(Y*X),
    lambda s,X,Y:s.block_collapse(3*X),
    lambda s,X,Y:s.block_collapse(X**2),
    lambda s,X,Y:s.block_collapse(X**3),
    lambda s,X,Y:X.T,lambda s,X,Y:X.adjoint(),lambda s,X,Y:X.conjugate(),
])
def test_block_arithmetic_matches_sympy(operation):
    value=operation(sx,block(sx),block(sx,'E'))
    expected=operation(sp,block(sp),block(sp,'E'))
    same(value,expected)
    assert value.func(*value.args)==value
    # srepr lists MatAdd terms in printing order (SymPy 1.14): compare printed forms.
    assert sx.srepr(sx.sympify(sx.srepr(value)))==sx.srepr(value)


def test_ordered_product_cells_and_unaligned_partitions():
    X,Y=block(sx),block(sx,'other')
    A,B,C,D=matrices(sx);E,F,G,H=matrices(sx,'other')
    result=sx.block_collapse(X*Y)
    assert result.blocks[0,0]==A*E+B*G
    assert result.blocks[1,1]==C*F+D*H
    assert sx.block_collapse(X*Y)!=sx.block_collapse(Y*X)
    unaligned=sx.BlockDiagMatrix(sx.MatrixSymbol('U',1,1),sx.MatrixSymbol('V',2,2))
    assert sx.block_collapse(X*unaligned)==X*unaligned
    assert sx.block_collapse(X+unaligned)==X+unaligned
    assert sx.block_collapse(sx.BlockMatrix([[A]]))==A
    assert sx.block_collapse(X.inverse())==X.inverse()  # No unproved Schur pivot.
    singular=sx.BlockDiagMatrix(sx.MatrixSymbol('R',2,3),sx.MatrixSymbol('S',3,2))
    assert sx.block_collapse(sx.Inverse(singular))==sx.Inverse(singular)
    with pytest.raises(Exception) as error:singular.inverse()
    assert type(error.value).__name__=='NonInvertibleMatrixError'


def test_elements_substitutions_trace_scalar_collapse_and_derivative():
    X=block(sx);A,B,C,D=matrices(sx)
    assert X[0,2]==B[0,0] and X[2,0]==C[0,0]
    assert X[-1,-1]==D[0,0]  # Normalize global negative indices before selecting block.
    assert X[-1,0]==C[0,0]
    assert sx.trace(X)==sx.trace(A)+sx.trace(D)
    assert sx.block_collapse(sx.Trace(X))==sx.trace(A)+sx.trace(D)
    assert isinstance(sx.det(X),sx.Determinant)
    x=sx.Symbol('x',real=True)
    value=sx.BlockMatrix([[x*A,B],[C,x*x*D]])
    same(value.diff(x),sp.BlockMatrix([[sp.MatrixSymbol('A',2,2),sp.ZeroMatrix(2,1)],[sp.ZeroMatrix(1,2),2*sp.Symbol('x',real=True)*sp.MatrixSymbol('D',1,1)]]))
    replacement=sx.BlockMatrix([[A,2*B],[3*C,D]])
    assert X.xreplace({B:2*B,C:3*C})==replacement
    assert X.subs(B,2*B).blocks[0,1]==2*B
    for key in [(3,0),(0,3),(-4,0)]:
        with pytest.raises(IndexError):X[key]


@pytest.mark.parametrize('protocol',range(pickle.HIGHEST_PROTOCOL+1))
def test_pickle_retains_typed_blocks_and_shared_cells(protocol):
    X=block(sx)
    restored=pickle.loads(pickle.dumps(X,protocol=protocol))
    assert restored==X and hash(restored)==hash(X)
    assert restored.args==(restored.blocks,)
    assert isinstance(restored,sx.BlockMatrix)
    assert sx.sympify(sx.srepr(restored))==X


def test_large_block_power_collapse_retains_shared_dag():
    source='''
import pickle
import symox as s
A,B,C,D=[s.MatrixSymbol(n,1,1) for n in 'ABCD']
value=s.block_collapse(s.BlockMatrix([[A,B],[C,D]])**1000000000)
assert value.shape==(2,2)
restore,args=value.__reduce__()
assert len(args[1])<10000
assert len(pickle.dumps(value))<1000000
assert pickle.loads(pickle.dumps(value))==value
'''
    result=subprocess.run([sys.executable,'-c',source],capture_output=True,text=True,timeout=25)
    assert result.returncode==0,result.stderr


def test_import_path_and_native_execution_without_sympy():
    from symox.matrices.expressions.blockmatrix import BlockMatrix,block_collapse
    assert BlockMatrix is sx.BlockMatrix and block_collapse is sx.block_collapse
    source='''
import sys
sys.modules['sympy']=None
import symox as s
A,B,C,D=[s.MatrixSymbol(n,1,1) for n in 'ABCD']
X=s.BlockMatrix([[A,B],[C,D]])
assert s.block_collapse(X*X).blocks[0,1]==A*B+B*D
assert X.as_explicit().shape==(2,2)
'''
    result=subprocess.run([sys.executable,'-c',source],capture_output=True,text=True,timeout=10)
    assert result.returncode==0,result.stderr


def test_empty_block_diagonal_power_and_inverse_collapse():
    empty=sx.BlockDiagMatrix()
    assert sx.block_collapse(sx.MatPow(empty,2)).shape==(0,0)
    assert sx.block_collapse(sx.Inverse(empty)).shape==(0,0)
