"""Native scalar calculus and scalar invariants of typed matrix expressions."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def assert_entries_equal(actual, expected):
    actual=actual.as_explicit().to_sympy()
    expected=expected.as_explicit() if isinstance(expected,sp.MatrixExpr) else expected
    assert actual.shape==expected.shape
    for a,b in zip(actual,expected):
        assert sp.cancel(sp.expand(a-b))==0


@pytest.mark.parametrize('operation',[
    lambda x,y,A,B:x*A+y*B,
    lambda x,y,A,B:(x*A+B)*(A+y*B),
    lambda x,y,A,B:(x*A+B)**2,
    lambda x,y,A,B:(x*A+B)**3,
    lambda x,y,A,B:(x*A+B)**-1,
    lambda x,y,A,B:(x*A+B)**-2,
    lambda x,y,A,B:(x*A+B).T,
    lambda x,y,A,B:(x*A+B).adjoint(),
    lambda x,y,A,B:(x*A+B).conjugate(),
    lambda x,y,A,B:x**3*A,
    lambda x,y,A,B:A*B,
])
def test_scalar_derivatives_match_explicit_sympy(operation):
    x,y=sx.symbols('x y',real=True);A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    X,Y=sp.symbols('x y',real=True);a=sp.MatrixSymbol('A',2,2);b=sp.MatrixSymbol('B',2,2)
    expression=operation(x,y,A,B);reference=operation(X,Y,a,b)
    result=expression.diff(x)
    assert result.shape==expression.shape
    expected=reference.as_explicit().diff(X)
    # Distinct, noncommuting matrices keep the scalar-parameter oracle exact
    # without spending the test on giant rational normalization in eight
    # unrelated matrix-entry symbols. The ordered identities below retain
    # abstract A and B independently of this evaluation certificate.
    replacements={A:sx.Matrix([[1,2],[0,1]]),B:sx.Matrix([[3,0],[1,4]])}
    oracle_replacements={a:sp.ImmutableMatrix([[1,2],[0,1]]),b:sp.ImmutableMatrix([[3,0],[1,4]])}
    assert_entries_equal(result.xreplace(replacements),expected.xreplace(oracle_replacements))
    assert sx.diff(expression,x)==result
    assert expression.diff(x,0)==expression


def test_ordered_product_and_inverse_formulas():
    x=sx.Symbol('x');A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    F=x*A+B
    assert ((x*A+B)*(A+x*B)).diff(x)==A*(A+x*B)+(x*A+B)*B
    assert (F**-1).diff(x)==-(F**-1)*A*(F**-1)
    assert ((F**2).diff(x))==A*F+F*A
    assert A.diff(x)==sx.ZeroMatrix(2,2)
    assert (A**sx.Symbol('n',integer=True)).diff(x)==sx.ZeroMatrix(2,2)
    with pytest.raises(NotImplementedError): (F**x).diff(x)
    with pytest.raises(NotImplementedError): (F**sx.Symbol('n',integer=True)).diff(x)


def test_mixed_higher_derivatives_and_explicit_entries():
    x,y=sx.symbols('x y',real=True);A=sx.MatrixSymbol('A',2,2)
    expression=x**3*y**2*A
    assert expression.diff((x,2),y)==12*x*y*A
    assert sx.diff(expression,x,2,y)==12*x*y*A
    assert expression.diff((x,4))==sx.ZeroMatrix(2,2)
    explicit=sx.MatrixExpr(sx.Matrix([[x**3,sx.sin(x)],[y,x*y]]))
    expected=sp.Matrix([[3*sp.Symbol('x',real=True)**2,sp.cos(sp.Symbol('x',real=True))],[0,sp.Symbol('y',real=True)]])
    assert_entries_equal(explicit.diff(x),expected)
    assert sx.MatrixExpr(sx.Matrix([[x]])).diff().as_explicit()==sx.ImmutableMatrix([[1]])
    n=sx.Symbol('n',integer=True,positive=True)
    assert sx.MatrixSymbol('N',n,n).diff(n)==sx.ZeroMatrix(n,n)
    assert sx.Identity(n).diff(n)==sx.ZeroMatrix(n,n)
    assert sx.ZeroMatrix(0,3).diff(x).shape==(0,3)


def test_complex_scalar_guards_and_invalid_variables():
    z=sx.Symbol('z');x=sx.Symbol('x',real=True);A=sx.MatrixSymbol('A',2,2)
    assert (z**2*A).diff(z)==2*z*A
    assert sx.Adjoint(A).diff(z)==sx.ZeroMatrix(2,2)
    assert sx.Adjoint(x*A).diff(x)==A.adjoint()
    assert sx.MatrixConjugate(x*A).diff(x)==A.conjugate()
    for expression in [sx.Adjoint(z*A),sx.MatrixConjugate(z*A),(z*A).adjoint(),(z*A).conjugate()]:
        with pytest.raises(NotImplementedError,match='real'): expression.diff(z)
    jacobian=(z*A).diff(A)
    assert jacobian.shape==(2,2,2,2)
    assert jacobian[0,1,0,1]==z and jacobian[0,1,1,0]==0
    with pytest.raises(NotImplementedError): (z*A).diff(A[0,0])
    with pytest.raises(ValueError): (z*A).diff(2)
    with pytest.raises(NotImplementedError): A.diff()
    with pytest.raises(ValueError): (x*A).diff()
    with pytest.raises(NotImplementedError): (x*A).diff(x,evaluate=False)
    with pytest.raises(NotImplementedError): (x*A).diff(x,unknown=True)


def test_scalar_matrix_invariants_and_exact_rules():
    n=sx.Symbol('n',integer=True,positive=True);x=sx.Symbol('x',real=True)
    A=sx.MatrixSymbol('A',n,n);B=sx.MatrixSymbol('B',n,n)
    assert isinstance(sx.Trace(A),sx.Trace)
    assert isinstance(sx.Determinant(A),sx.Determinant)
    assert sx.Trace(A).args==(A,)
    assert sx.Determinant(A).args==(A,)
    assert sx.Trace(A).free_symbols=={A}
    assert sx.Determinant(A).free_symbols=={A}
    for value in [sx.Trace(A),sx.Determinant(A),sx.Trace(x*A+B),sx.Determinant(x*A+B)]:
        assert value.func(*value.args)==value
        assert pickle.loads(pickle.dumps(value))==value
        assert sx.sympify(sx.srepr(value))==value
        assert sx.sympify(sp.srepr(value.to_sympy()))==value
    assert sx.trace(A)==sx.Trace(A)
    assert sx.det(A)==sx.Determinant(A)
    assert sx.trace(x*A+B)==x*sx.Trace(A)+sx.Trace(B)
    assert sx.det(x*A*B)==x**n*sx.Determinant(A)*sx.Determinant(B)
    assert (A**3).det()==sx.Determinant(A)**3
    assert A.inverse().det()==1/sx.Determinant(A)
    assert sx.trace(A.T)==sx.Trace(A)
    assert sx.trace(A.adjoint())==sx.conjugate(sx.Trace(A))
    assert sx.det(A.adjoint())==sx.conjugate(sx.Determinant(A))
    assert sx.Trace(A+B).doit()==sx.Trace(A)+sx.Trace(B)
    assert sx.Determinant(A*B).doit()==sx.Determinant(A)*sx.Determinant(B)
    assert sx.trace(sx.Identity(n))==n
    assert sx.det(sx.Identity(n))==1
    assert sx.trace(sx.ZeroMatrix(n,n))==0
    assert sx.det(sx.ZeroMatrix(n,n))==0
    assert sx.det(sx.ZeroMatrix(0,0))==1
    for constructor in [sx.Trace,sx.Determinant]:
        with pytest.raises(sx.NonSquareMatrixError): constructor(sx.MatrixSymbol('Rect',2,3))


def test_trace_of_rectangular_products_keeps_sympy_order_and_explicit_lowering():
    A=sx.MatrixSymbol('A',2,3);B=sx.MatrixSymbol('B',3,2)
    # SymPy 1.14 MatMul._eval_trace does not canonicalize cycles or transposes.
    assert sx.trace(A*B)!=sx.trace(B*A)
    assert str(sx.trace(A*B))=='Trace(A*B)' and str(sx.trace(B*A))=='Trace(B*A)'
    assert str(sx.trace((A*B).T))=='Trace(B.T*A.T)'
    assert str(sx.trace(2*A*B))=='2*Trace(A*B)'
    square=sx.MatrixSymbol('Square',2,2)
    assert square.trace_explicit()==square[0,0]+square[1,1]
    assert square.determinant_explicit()==square[0,0]*square[1,1]-square[0,1]*square[1,0]
    assert sx.det(A*B)==sx.Determinant(A*B)
    n=sx.Symbol('n',integer=True,positive=True)
    with pytest.raises(NotImplementedError): sx.MatrixSymbol('N',n,n).trace_explicit()
    with pytest.raises(NotImplementedError): sx.MatrixSymbol('N',n,n).determinant_explicit()
    explicit=sx.MatrixExpr(sx.Matrix([[1,2],[3,4]]))
    assert explicit.trace()==5 and explicit.det()==-2


def test_trace_derivative_and_substitution_to_explicit_matrix():
    x=sx.Symbol('x',real=True);A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    assert sx.Trace(x**3*A).diff(x)==3*x**2*sx.Trace(A)
    assert sx.Trace(A).diff(x)==0
    assert sx.Trace(A).subs(A,sx.Matrix([[1,2],[3,4]])).doit()==5
    assert sx.Determinant(A).xreplace({A:sx.Matrix([[1,2],[3,4]])}).doit()==-2
    assert sx.Trace(A).xreplace({A:B})==sx.Trace(B)
    held=sx.Trace(x*A+B)
    assert held.xreplace({})==held
    assert held.subs(sx.Symbol('unrelated'),2)==held
    held=sx.Determinant(A*B)
    assert held.xreplace({})==held
    assert held.subs(sx.Symbol('unrelated'),2)==held


def test_matrix_calculus_runs_without_sympy_import():
    source='''
import sys
import symox as s
x=s.Symbol('x',real=True); A=s.MatrixSymbol('A',2,2); B=s.MatrixSymbol('B',2,2)
F=x*A+B
assert (F**-1).diff(x)==-(F**-1)*A*(F**-1)
assert s.trace(F)==x*s.Trace(A)+s.Trace(B)
assert s.det(A*B)==s.Determinant(A)*s.Determinant(B)
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,capture_output=True,text=True)


def test_determinant_derivative_valid_at_a_singular_matrix():
    x=sx.Symbol('x',real=True);A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    derivative=sx.Determinant(x*A+B).diff(x)
    # det(diag(x,1))=x, so its derivative at the singular point x=0 is 1.
    # det(F)*trace(F^-1*F') introduces a spurious 0/0 at this point.
    values={A:sx.Matrix([[1,0],[0,0]]),B:sx.Matrix([[0,0],[0,1]]),x:0}
    assert derivative.xreplace(values)==1
    assert derivative.subs(A,sx.Matrix([[1,0],[0,0]])).subs(B,sx.Matrix([[0,0],[0,1]])).subs(x,0)==1


def test_nonreal_conjugate_trace_derivative_remains_symbolic():
    z=sx.Symbol('z');A=sx.MatrixSymbol('A',2,2)
    derivative=sx.Trace((z*A).adjoint()).diff(z)
    Z=sp.Symbol('z');a=sp.MatrixSymbol('A',2,2)
    expected=sp.diff(sp.Trace((Z*a).adjoint()).doit(),Z)
    assert derivative.to_sympy()==expected


def test_held_zero_and_identity_derivatives_short_circuit():
    z=sx.Symbol('z');A=sx.MatrixSymbol('A',2,2)
    dependent=sx.Adjoint(z*A)
    assert sx.MatPow(dependent,0).diff(z)==sx.ZeroMatrix(2,2)
    assert sx.MatMul(0,dependent).diff(z)==sx.ZeroMatrix(2,2)
