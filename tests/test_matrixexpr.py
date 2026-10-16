"""Typed matrix expressions: noncommutative algebra, elements and interop."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def equal(actual,expected):
    actual=actual.to_sympy()
    assert actual.shape==expected.shape
    difference=(actual-expected).doit().expand().doit()
    if all(sp.sympify(dimension).is_Integer for dimension in expected.shape):
        difference=difference.as_explicit() if not isinstance(difference,sp.MatrixBase) else difference
    if isinstance(difference,sp.MatrixBase):
        assert all(sp.simplify(entry)==0 for entry in difference)
    else:
        assert difference == sp.ZeroMatrix(*expected.shape)


@pytest.mark.parametrize('operation',[
    lambda A,B:A+B,lambda A,B:A+A,lambda A,B:2*A+3*A,
    lambda A,B:2*(A+B),lambda A,B:(A+B)*A,
    lambda A,B:A*B,lambda A,B:B*A,lambda A,B:A**3,
    lambda A,B:A**0,lambda A,B:A*A**-1,
    lambda A,B:A*B*B**-1,lambda A,B:(A**2)**3,
    lambda A,B:(A+B).T,lambda A,B:(A*B).T,
    lambda A,B:(A*B).adjoint(),lambda A,B:A.T.T,
    lambda A,B:A.adjoint().adjoint(),lambda A,B:(2*A).T,
    lambda A,B:(A*B).inverse(),lambda A,B:A/2,
    lambda A,B:-A+B,
])
def test_matrix_algebra_against_sympy(operation):
    a,b=sx.MatrixSymbol('A',2,2),sx.MatrixSymbol('B',2,2)
    A,B=sp.MatrixSymbol('A',2,2),sp.MatrixSymbol('B',2,2)
    equal(operation(a,b),operation(A,B))


def test_noncommutative_products_and_expansion():
    A,B,C=(sx.MatrixSymbol(name,2,2) for name in 'ABC')
    assert A*B != B*A
    assert (A*B).args == (A,B)
    assert (A*B).T == B.T*A.T
    assert (A*B).H == B.H*A.H
    equal(((A+B)*C).expand(),(sp.MatrixSymbol('A',2,2)+sp.MatrixSymbol('B',2,2))*sp.MatrixSymbol('C',2,2))
    assert ((A+B)*C).expand() == A*C+B*C


def test_classes_args_and_reconstruction():
    A=sx.MatrixSymbol('A',2,2)
    B=sx.MatrixSymbol('B',2,2)
    values=[A,A+B,A*B,A**3,A.T,A.H,A.inverse(),sx.Identity(2),sx.ZeroMatrix(2,3)]
    classes=[sx.MatrixSymbol,sx.MatAdd,sx.MatMul,sx.MatPow,sx.Transpose,sx.Adjoint,sx.Inverse,sx.Identity,sx.ZeroMatrix]
    for value,kind in zip(values,classes):
        assert isinstance(value,sx.MatrixExpr)
        assert isinstance(value,kind)
        assert value.func(*value.args)==value
        assert pickle.loads(pickle.dumps(value))==value
        assert hash(pickle.loads(pickle.dumps(value)))==hash(value)
    assert isinstance(A.inverse(),sx.MatPow)
    assert sx.sympify(A) is A or sx.sympify(A)==A
    assert sx.Add(A,A)==2*A
    assert sx.Mul(A,B)==A*B
    assert sx.Pow(A,2)==A**2
    assert not A.is_commutative


def test_raw_constructors_and_doit():
    A=sx.MatrixSymbol('A',2,2)
    raw=sx.MatAdd(A,A)
    assert raw.args==(A,A)
    assert raw.doit()==2*A
    raw=sx.MatMul(A,sx.Identity(2))
    assert raw.args==(A,sx.Identity(2))
    assert raw.doit()==A
    assert raw.xreplace({})==raw
    assert sx.MatPow(A,0).doit()==sx.Identity(2)
    assert sx.Transpose(A+A).doit()==2*A.T


def test_symbolic_dimensions_substitution_and_shape_errors():
    n,m,p=sx.symbols('n m p',integer=True,positive=True)
    A=sx.MatrixSymbol('A',n,m)
    B=sx.MatrixSymbol('B',m,p)
    assert (A*B).shape==(n,p)
    assert (A+B.T).shape==(n,m)
    assert A.free_symbols=={A}
    assert sx.Identity(n).free_symbols=={n}
    assert sx.ZeroMatrix(n,m).free_symbols=={n,m}
    assert A.xreplace({n:m,m:p}).shape==(m,p)
    assert A.subs({n:2,m:3}).shape==(2,3)
    with pytest.raises(NotImplementedError): A.as_explicit()
    with pytest.raises(sx.ShapeError): sx.MatrixSymbol('C',2,3)*sx.MatrixSymbol('D',4,2)
    with pytest.raises(sx.ShapeError): sx.MatrixSymbol('C',2,3)+sx.MatrixSymbol('D',2,4)
    with pytest.raises(sx.ShapeError): sx.MatrixSymbol('C',2,3).inverse()
    for dim in [-1,sx.Rational(1,2),sx.pi]:
        with pytest.raises(ValueError): sx.MatrixSymbol('bad',dim,2)


def test_matrix_replacements_and_free_symbols():
    x=sx.Symbol('x',positive=True)
    A,B,C=(sx.MatrixSymbol(name,2,2) for name in 'ABC')
    expression=x*A*B+C
    assert expression.free_symbols=={x,A,B,C}
    assert expression.xreplace({A:C,B:A,x:2})==2*C*A+C
    assert expression.subs(A,C)==x*C*B+C
    assert expression.subs({A:sx.Identity(2),B:sx.ZeroMatrix(2,2)})==C
    restored=pickle.loads(pickle.dumps(expression))
    assert restored==expression and x in restored.free_symbols
    assert any(getattr(v,'is_positive',None) is True for v in restored.free_symbols if not isinstance(v,sx.MatrixExpr))


def test_scalar_elements_and_finite_contraction():
    A=sx.MatrixSymbol('A',2,3);B=sx.MatrixSymbol('B',3,2)
    element=A[0,1]
    assert isinstance(element,sx.MatrixElement)
    assert element.args==(A,0,1)
    assert element.func(*element.args)==element
    assert element.free_symbols=={A}
    assert (A*B)[0,1]==sum((A[0,k]*B[k,1] for k in range(3)),sx.Integer(0))
    assert (A.T)[1,0]==element
    assert element.xreplace({A:sx.Matrix([[1,2,3],[4,5,6]])})==2
    assert element.subs(A,sx.Matrix([[1,2,3],[4,5,6]]))==2
    assert pickle.loads(pickle.dumps(element))==element
    assert (element+3).free_symbols=={A}
    assert sx.Identity(3)[0,0]==1 and sx.Identity(3)[0,1]==0
    assert sx.ZeroMatrix(2,3)[1,2]==0
    with pytest.raises(IndexError): A[2,0]
    with pytest.raises(IndexError): A[0,3]


def test_explicit_conversion_and_dense_interop():
    A=sx.MatrixSymbol('A',2,2)
    B=sx.MatrixSymbol('B',2,2)
    dense=sx.Matrix([[1,2],[3,4]])
    for expr,expected in [(A+dense,sp.MatrixSymbol('A',2,2)+sp.Matrix([[1,2],[3,4]])),(dense+A,sp.Matrix([[1,2],[3,4]])+sp.MatrixSymbol('A',2,2)),(A*dense,sp.MatrixSymbol('A',2,2)*sp.Matrix([[1,2],[3,4]])),(dense*A,sp.Matrix([[1,2],[3,4]])*sp.MatrixSymbol('A',2,2))]:
        equal(expr,expected)
        result=expr.as_explicit()
        assert isinstance(result,sx.ImmutableMatrix)
        assert result.to_sympy()==expected.as_explicit()
    assert (A*B).as_explicit()==A.as_explicit()*B.as_explicit()
    assert (A**2).as_explicit()==A.as_explicit()**2
    assert A.xreplace({A:dense}).as_explicit()==dense
    assert (sx.Identity(2)*A)==A
    assert (A*sx.ZeroMatrix(2,3))==sx.ZeroMatrix(2,3)


def test_sympy_interchange_and_import_paths():
    from symox.matrices.expressions.matexpr import MatrixSymbol
    from symox.matrices.expressions.matmul import MatMul
    from symox.matrices import MatrixExpr
    assert MatrixSymbol is sx.MatrixSymbol and MatMul is sx.MatMul and MatrixExpr is sx.MatrixExpr
    n=sp.Symbol('n',integer=True,positive=True)
    A=sp.MatrixSymbol('A',n,n);B=sp.MatrixSymbol('B',n,n)
    for reference in [A,A+B,A*B,A.T,A.adjoint(),A**3,A.inverse(),sp.Identity(n),sp.ZeroMatrix(n,n)]:
        value=sx.MatrixExpr(reference)
        equal(value,reference)
        assert sx.MatrixExpr(value.to_sympy())==value


def test_native_runtime_needs_no_sympy():
    source='''
import sys
import symox as s
A=s.MatrixSymbol('A',2,2); B=s.MatrixSymbol('B',2,2)
assert A*B != B*A
assert (A*B)[0,0] == A[0,0]*B[0,0]+A[0,1]*B[1,0]
assert (A*B).as_explicit() == A.as_explicit()*B.as_explicit()
assert A[0,0].free_symbols == {A}
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,capture_output=True,text=True)


def test_srepr_parser_roundtrips_native_and_sympy():
    n=sx.Symbol('n',integer=True,positive=True)
    A=sx.MatrixSymbol('A',n,n);B=sx.MatrixSymbol('B',n,n)
    values=[A,A+B,A*B,2*A,A**3,A.T,A.adjoint(),A.inverse(),sx.Identity(n),sx.ZeroMatrix(n,n),A[0,1],(A*B)[0,1],A+sx.MatrixSymbol('C',n,n)]
    for value in values:
        assert sx.sympify(sx.srepr(value))==value
        reference=value.to_sympy()
        assert sx.sympify(sp.srepr(reference))==value
    finite=sx.MatrixSymbol('F',2,2)
    for matrix in [sx.Matrix([[1,2],[3,4]]),sx.Matrix([[1]]),sx.Matrix(0,0,[])]:
        value=finite+matrix if matrix.shape==(2,2) else sx.MatrixExpr(matrix)
        assert sx.sympify(sx.srepr(value))==value
        assert sx.sympify(sp.srepr(value.to_sympy()))==value


def test_xreplace_matrix_replacements_are_simultaneous():
    n=sx.Symbol('n',integer=True,positive=True)
    A=sx.MatrixSymbol('A',n,n);B=sx.MatrixSymbol('B',n,n)
    assert A.xreplace({A:B,n:3})==B
    expression=2*A
    assert expression.xreplace({A:B,n:3})==2*B
    raw=sx.MatMul(sx.Identity(n),A)
    assert raw.xreplace({})==raw


def test_explicit_wrapper_equality_and_hash_match_immutable_matrices():
    dense=sx.ImmutableMatrix([[1,0],[2,3]])
    wrapper=sx.MatrixExpr(dense)
    assert wrapper==dense and dense==wrapper
    assert hash(wrapper)==hash(dense)
    sparse=sx.ImmutableSparseMatrix(dense)
    assert wrapper==sparse and sparse==wrapper
    assert hash(wrapper)==hash(sparse)
    assert len({wrapper,dense,sparse})==1
