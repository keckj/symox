"""Native permutation matrix semantics, structural protocols and codecs."""
import itertools
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from sympy.combinatorics import Permutation as SP
from symox.combinatorics import Permutation
from symox._matrix_permutations import MatrixPermute, PermutationMatrix


CASES = [list(images) for size in range(5)
         for images in itertools.permutations(range(size))]
CASES += [[1,0,3,2,4],[0,2,1,4,5,3]]


def converted(value):
    if isinstance(value,(tuple,list)):
        return type(value)(converted(v) for v in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


@pytest.mark.parametrize('images',CASES)
def test_permutation_matrices_and_structural_roundtrips(images):
    permutation = Permutation(images)
    matrix = PermutationMatrix(permutation)
    expected = sp.PermutationMatrix(SP(images))
    assert isinstance(matrix,sx.MatrixExpr) and isinstance(matrix,PermutationMatrix)
    assert matrix.args == (permutation,)
    assert matrix.func(*matrix.args) == matrix
    assert matrix.shape == (len(images),len(images))
    assert matrix.is_Identity == permutation.is_Identity
    assert matrix.free_symbols == set()
    assert matrix.as_explicit().to_sympy() == expected.as_explicit()
    assert matrix.to_sympy() == expected
    assert sx.MatrixExpr(expected) == matrix
    assert sx.sympify(sx.srepr(matrix)) == matrix
    assert sx.sympify(sp.srepr(expected)) == matrix
    assert converted(matrix.sort_key()) == expected.sort_key()
    assert matrix.atoms() == {permutation}
    assert matrix.has(permutation) and matrix.has(Permutation)
    assert matrix.has_free(permutation)
    assert not matrix.has(sx.Integer(len(images)))
    assert matrix.diff(sx.Symbol('x')).as_explicit() == sx.zeros(len(images))
    assert matrix.det() == permutation.signature()
    assert sx.trace(matrix) == sum(i==j for i,j in enumerate(images))
    assert matrix.T.to_sympy() == expected.T
    assert matrix.adjoint().to_sympy() == expected.adjoint()
    assert matrix.conjugate() == matrix
    assert matrix.inverse().to_sympy() == expected.inverse()
    if permutation.is_Identity:
        assert matrix.doit() == sx.Identity(len(images))
    for exponent in [-7,-1,0,1,2,7,10**18]:
        assert (matrix**exponent).as_explicit().to_sympy() == (expected**exponent).as_explicit()
    block_form = matrix.rewrite(sx.BlockDiagMatrix)
    assert block_form.as_explicit() == matrix.as_explicit()
    assert matrix.rewrite(sx.sin,sx.BlockDiagMatrix) == matrix


@pytest.mark.parametrize('images',CASES)
def test_composition_preserves_permutation_convention(images):
    q_images = list(reversed(range(len(images))))
    p,q = Permutation(images),Permutation(q_images)
    first,second = PermutationMatrix(p),PermutationMatrix(q)
    assert (first*second).as_explicit() == first.as_explicit()*second.as_explicit()
    assert (first*second).as_explicit() == PermutationMatrix(p*q).as_explicit()
    assert first*first.T == sx.Identity(len(images))
    assert first.T*first == sx.Identity(len(images))
    for predicate in [sx.Q.square,sx.Q.invertible,sx.Q.orthogonal,sx.Q.unitary,
                      sx.Q.normal,sx.Q.fullrank,sx.Q.integer_elements,
                      sx.Q.real_elements,sx.Q.complex_elements]:
        assert sx.ask(predicate(first)) is True
    assert sx.ask(sx.Q.singular(first)) is False
    assert sx.ask(sx.Q.symmetric(first)) == (p*p).is_Identity
    for predicate in [sx.Q.diagonal,sx.Q.upper_triangular,sx.Q.lower_triangular]:
        assert sx.ask(predicate(first)) == p.is_Identity


@pytest.mark.parametrize('axis',[0,1])
@pytest.mark.parametrize('images',[[],[0],[1,0],[1,2,0],[1,0,2,3],[1,0,3,2]])
@pytest.mark.parametrize('kind',['explicit','symbol','function'])
def test_row_column_permutation_and_resizing(axis,images,kind):
    A=sx.MatrixSymbol('A',3,4)
    expected_A=sp.MatrixSymbol('A',3,4)
    if kind=='explicit':
        A=sx.Matrix(3,4,list(range(12)))
        expected_A=sp.Matrix(3,4,list(range(12)))
    elif kind=='function':
        i,j=sx.symbols('i j')
        I,J=sp.symbols('i j')
        A=sx.FunctionMatrix(3,4,sx.Lambda((i,j),i+2*j))
        expected_A=sp.FunctionMatrix(3,4,sp.Lambda((I,J),I+2*J))
    p=Permutation(images)
    try:
        expected=sp.MatrixPermute(expected_A,SP(images),axis)
    except ValueError:
        with pytest.raises((ValueError,sx.ShapeError)):
            MatrixPermute(A,p,axis)
        return
    value=MatrixPermute(A,p,axis)
    assert value.to_sympy()==expected
    assert value.as_explicit().to_sympy()==expected.as_explicit()
    assert value.func(*value.args)==value
    assert converted(value.sort_key())==expected.sort_key()
    assert sx.MatrixExpr(expected)==value
    assert sx.sympify(sx.srepr(value))==value
    assert sx.sympify(sp.srepr(expected))==value
    assert value.rewrite(sx.MatMul).as_explicit()==value.as_explicit()
    assert value.rewrite('MatMul').as_explicit()==value.as_explicit()
    assert value.T.as_explicit()==value.as_explicit().T
    assert value.adjoint().as_explicit()==value.as_explicit().adjoint()
    assert value.conjugate().as_explicit()==value.as_explicit().conjugate()
    assert value.evalf(30).shape==value.shape
    restored=pickle.loads(pickle.dumps(value))
    assert restored==value and hash(restored)==hash(value)


@pytest.mark.parametrize('axis',[0,1])
def test_permutation_calculus_substitutions_and_composition(axis):
    x=sx.Symbol('x',real=True)
    A=sx.MatrixSymbol('A',3,3)
    B=sx.MatrixSymbol('B',3,3)
    p=Permutation([1,2,0])
    value=MatrixPermute(x*x*A,p,axis)
    assert value.diff(x).as_explicit()==2*x*MatrixPermute(A,p,axis).as_explicit()
    assert value.subs(x,3).as_explicit()==9*MatrixPermute(A,p,axis).as_explicit()
    assert value.subs(A,B)==MatrixPermute(x*x*B,p,axis)
    assert value.xreplace({A:B})==MatrixPermute(x*x*B,p,axis)
    assert MatrixPermute(A,p,0).subs(sx.Integer(0),sx.Integer(1))==MatrixPermute(A,p,1)
    assert MatrixPermute(MatrixPermute(A,p,axis),~p,axis).doit()==A
    assert MatrixPermute(A,PermutationMatrix(p),axis)==MatrixPermute(A,p,axis)
    assert MatrixPermute(sx.Identity(3),p,axis).doit().as_explicit()==(
        PermutationMatrix(p if axis==0 else ~p).as_explicit())
    assert value.det()==x**6*sx.det(A)
    assert value.free_symbols=={x,A}
    assert value.has(p) and value.has_free(p)


def test_symbolic_indices_are_substitutable():
    i,j=sx.symbols('i j',integer=True)
    P=PermutationMatrix(Permutation([2,0,1]))
    value=P[i,j]
    A=sx.MatrixSymbol('A',3,2)
    row=MatrixPermute(A,Permutation([2,0,1]))
    entry=row[i,j]
    for r in range(-3,3):
        for c in range(-3,3):
            assert value.subs({i:r,j:c})==P[r,c]
        for c in range(2):
            assert entry.subs({i:r,j:c})==A[[2,0,1][r],c]
    assert P[-1,-2]==1
    assert P[-3,-1]==1


@pytest.mark.parametrize('protocol',range(pickle.HIGHEST_PROTOCOL+1))
def test_pickle_shared_matrix_parents(protocol):
    A=sx.MatrixSymbol('A',3,3)
    P=PermutationMatrix(Permutation([1,2,0]))
    permuted=MatrixPermute(A,P,0)
    values=[P,permuted,permuted.T,P+P,sx.Trace(permuted),sx.Determinant(permuted),
            sx.MatrixElement(permuted,0,1)]
    for value in values:
        restored=pickle.loads(pickle.dumps(value,protocol=protocol))
        assert restored==value and hash(restored)==hash(value)
        assert restored.free_symbols==value.free_symbols


def test_invalid_constructor_and_pickle_metadata():
    A=sx.MatrixSymbol('A',3,3)
    for invalid in [[1,0],1,sx.Symbol('p')]:
        with pytest.raises(ValueError):
            PermutationMatrix(invalid)
        with pytest.raises(ValueError):
            MatrixPermute(A,invalid)
    for axis in [-1,2,1.5,sx.Symbol('axis')]:
        with pytest.raises(ValueError):
            MatrixPermute(A,Permutation([1,0]),axis)
    with pytest.raises(ValueError):
        MatrixPermute([1,2],Permutation([1,0]))
    value=MatrixPermute(A,Permutation([1,2,0]))
    restore,(version,nodes,root)=value.__reduce__()
    tag,canonical,fields=nodes[root]
    assert tag=='MatrixPermute'
    for invalid in [(fields[0],[1,1,0],0),(fields[0],[1,0],0),
                    (fields[0],[1,2,0],2)]:
        modified=list(nodes)
        modified[root]=(tag,canonical,invalid)
        with pytest.raises(ValueError):
            restore(version,modified,root)


def test_no_sympy_import_and_graph_factors():
    source='''
import sys,pickle
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import Matrix,PermutationMatrix,MatrixPermute,MatrixSymbol,Identity
from symox.combinatorics import Permutation
p=Permutation([1,2,0]); P=PermutationMatrix(p)
assert P*P.T==Identity(3)
A=MatrixSymbol('A',3,3)
M=MatrixPermute(A,p)
assert pickle.loads(pickle.dumps(M))==M
matrix=Matrix([[1,0,2],[0,3,0],[4,0,5]])
P,B=matrix.connected_components_decomposition()
assert isinstance(P,PermutationMatrix)
assert (P.T*B*P).as_explicit()==matrix
'''
    subprocess.run([sys.executable,'-c',source],check=True)
