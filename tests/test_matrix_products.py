"""Typed elementwise and tensor products, evaluated without a SymPy runtime."""
import ast
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx


def same_entries(actual, expected):
    actual=actual.as_explicit().to_sympy()
    expected=expected.as_explicit() if isinstance(expected,sp.MatrixExpr) else expected
    assert actual.shape==expected.shape
    for a,b in zip(actual,expected):
        assert sp.simplify(a-b)==0


def test_one_matrix_is_entrywise_identity_and_rectangular():
    one=sx.OneMatrix(2,3);A=sx.MatrixSymbol('A',2,3)
    assert isinstance(one,sx.OneMatrix) and one.is_OneMatrix
    assert one.shape==(2,3) and one[0,2]==1
    assert one.as_explicit()==sx.ImmutableMatrix.ones(2,3)
    assert one.T==sx.OneMatrix(3,2)
    assert one.conjugate()==one
    assert sx.hadamard_product(A,one)==A
    assert sx.hadamard_product(A,sx.ZeroMatrix(2,3))==sx.ZeroMatrix(2,3)
    assert sx.HadamardPower(A,0).doit()==one
    assert sx.HadamardPower(A,1).doit()==A
    assert sx.OneMatrix(2,2).trace()==2 and sx.OneMatrix(2,2).det()==0
    assert sx.OneMatrix(1,1).det()==1 and sx.OneMatrix(0,0).det()==1
    with pytest.raises(sx.NonInvertibleMatrixError): sx.OneMatrix(2,2).inverse()
    assert sx.OneMatrix(0,3).as_explicit().shape==(0,3)
    n,m=sx.symbols('n m',integer=True,positive=True)
    assert sx.OneMatrix(n,m).free_symbols=={n,m}
    assert sx.OneMatrix(n,m).diff(n)==sx.ZeroMatrix(n,m)


def test_hadamard_canonicalization_raw_constructors_and_methods():
    A=sx.MatrixSymbol('A',2,3);B=sx.MatrixSymbol('B',2,3)
    raw=sx.HadamardProduct(B,A,A)
    assert raw.args==(B,A,A) and isinstance(raw,sx.HadamardProduct)
    assert raw.xreplace({})==raw
    canonical=raw.doit()
    assert canonical==sx.hadamard_product(A,A,B)==sx.hadamard_product(B,A,A)
    assert set(canonical.args)=={sx.HadamardPower(A,2),B}
    assert raw[1,2]==A[1,2]**2*B[1,2]
    assert sx.hadamard_product(A)==A
    assert sx.HadamardProduct(A).doit()==A
    assert sx.hadamard_power(A,1)==A
    assert sx.hadamard_power(2,3)==8
    with pytest.raises(ValueError): sx.hadamard_power(A,B)
    assert sx.HadamardPower(2,3)==8
    assert sx.HadamardPower(A,B).base==A and sx.HadamardPower(A,B).exp==B
    assert sx.HadamardPower(2,A).base==2 and sx.HadamardPower(2,A).exp==A
    with pytest.raises(ValueError): sx.HadamardProduct()
    with pytest.raises(TypeError): sx.hadamard_product()
    with pytest.raises(TypeError): sx.HadamardProduct(A,2)
    with pytest.raises(sx.ShapeError): sx.HadamardProduct(A,sx.MatrixSymbol('C',3,2))
    with pytest.raises(sx.ShapeError): sx.HadamardPower(A,sx.MatrixSymbol('C',3,2))
    with pytest.raises(NotImplementedError): sx.HadamardProduct(A,B,check=False)


@pytest.mark.parametrize('operation',[
    lambda s,A,B:s.HadamardProduct(A,B),
    lambda s,A,B:s.HadamardPower(A,3),
    lambda s,A,B:s.HadamardPower(A,s.Rational(1,2)),
    lambda s,A,B:s.HadamardPower(A,B),
    lambda s,A,B:s.HadamardPower(2,A),
    lambda s,A,B:s.HadamardProduct(A,B).T,
    lambda s,A,B:s.HadamardProduct(A,B).adjoint(),
    lambda s,A,B:s.HadamardPower(A,3).adjoint(),
])
def test_hadamard_entries_against_sympy(operation):
    A=sx.MatrixSymbol('A',2,3);B=sx.MatrixSymbol('B',2,3)
    a=sp.MatrixSymbol('A',2,3);b=sp.MatrixSymbol('B',2,3)
    same_entries(operation(sx,A,B),operation(sp,a,b))


def test_matrix_exponent_transpose_has_correct_rectangular_shape():
    # SymPy 1.14's direct .T rule forgets to transpose a matrix exponent.
    # Certify the mathematical entries independently of that oracle rule.
    A=sx.MatrixSymbol('A',2,3);B=sx.MatrixSymbol('B',2,3)
    actual=sx.HadamardPower(A,B).T
    assert actual.shape==(3,2)
    assert actual==sx.HadamardPower(A.T,B.T)
    expected=sp.Matrix(2,3,lambda i,j:sp.MatrixSymbol('A',2,3)[i,j]**sp.MatrixSymbol('B',2,3)[i,j]).T
    same_entries(actual,expected)


def test_kronecker_order_rectangularity_and_identities():
    A=sx.MatrixSymbol('A',2,3);B=sx.MatrixSymbol('B',3,2)
    product=sx.KroneckerProduct(A,B)
    assert product.shape==(6,6) and product.args==(A,B)
    assert product!=sx.KroneckerProduct(B,A)
    assert product[0,3]==A[0,1]*B[0,1]
    assert sx.KroneckerProduct(A).doit()==A
    assert sx.kronecker_product(A)==A
    assert sx.KroneckerProduct(sx.Identity(2),sx.Identity(3))==sx.Identity(6)
    assert sx.KroneckerProduct()==sx.ImmutableMatrix([[1]])
    assert sx.kronecker_product(sx.ZeroMatrix(0,2),B).shape==(0,4)
    assert sx.kronecker_product(2*A,B)==2*sx.kronecker_product(A,B)
    assert sx.KroneckerProduct(A,sx.KroneckerProduct(B,A)).doit()==sx.KroneckerProduct(A,B,A).doit()
    assert product.T==sx.kronecker_product(A.T,B.T)
    assert product.adjoint()==sx.kronecker_product(A.adjoint(),B.adjoint())
    assert product.conjugate()==sx.kronecker_product(A.conjugate(),B.conjugate())
    a=sp.MatrixSymbol('A',2,3);b=sp.MatrixSymbol('B',3,2)
    same_entries(product,sp.KroneckerProduct(a,b))
    # Factors are rectangular: splitting traces/determinants is invalid.
    assert product.trace()==sx.Trace(product)
    assert product.det()==sx.Determinant(product)
    with pytest.raises(TypeError): sx.kronecker_product()
    with pytest.raises(TypeError): sx.KroneckerProduct(A,2)
    with pytest.raises(NotImplementedError): sx.KroneckerProduct(A,B,check=False)


def test_square_kronecker_scalar_invariants_and_inverse():
    A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',3,3)
    product=sx.kronecker_product(A,B)
    assert product.trace()==sx.Trace(A)*sx.Trace(B)
    assert product.det()==sx.Determinant(A)**3*sx.Determinant(B)**2
    assert product.inverse()==sx.kronecker_product(A.inverse(),B.inverse())


def test_symbolic_shapes_indices_substitutions_and_interchange():
    n,m=sx.symbols('n m',integer=True,positive=True);x=sx.Symbol('x')
    A=sx.MatrixSymbol('A',n,m);B=sx.MatrixSymbol('B',n,m)
    values=[sx.OneMatrix(n,m),sx.HadamardProduct(A,B),sx.HadamardPower(A,x),sx.HadamardPower(x,A),sx.HadamardPower(A,B),sx.KroneckerProduct(A,B)]
    for value in values:
        assert value.func(*value.args)==value
        assert pickle.loads(pickle.dumps(value))==value
        assert sx.sympify(sx.srepr(value))==value
        assert sx.sympify(sp.srepr(value.to_sympy()))==value
        assert sx.MatrixExpr(value.to_sympy())==value
        assert value.xreplace({})==value
        assert value.subs(n,2).subs(m,3).shape==((4,9) if isinstance(value,sx.KroneckerProduct) else (2,3))
    assert sx.HadamardPower(x,A).free_symbols=={x,A}
    i,j=sx.symbols('i j',integer=True,nonnegative=True)
    product=sx.KroneckerProduct(A,B)
    element=product[i,j]
    aa=sx.MatrixSymbol('A',2,3);bb=sx.MatrixSymbol('B',2,3)
    finite=element.xreplace({n:2,m:3,i:3,j:5})
    assert finite==sx.KroneckerProduct(aa,bb)[3,5]
    raw=sx.KroneckerProduct(sx.MatrixSymbol('Unknown',n,n),sx.MatrixSymbol('MaybeZero',sx.Symbol('k'),sx.Symbol('k')))
    assert isinstance(raw[i,j],sx.MatrixElement)


def test_explicit_products_preserve_immutable_values_and_nonfinite_zero_terms():
    A=sx.Matrix([[1,2],[3,4]]);B=sx.Matrix([[2,3],[4,5]])
    result=sx.hadamard_product(A,B)
    assert isinstance(result,sx.ImmutableMatrix) and result==sx.ImmutableMatrix([[2,6],[12,20]])
    result=sx.kronecker_product(A,B)
    assert isinstance(result,sx.ImmutableMatrix)
    assert np.array(result.tolist(),dtype=int).tolist()==np.kron([[1,2],[3,4]],[[2,3],[4,5]]).tolist()
    for nonfinite in [sx.oo,-sx.oo,sx.zoo,sx.nan]:
        for operands in [(sx.ZeroMatrix(1,1),sx.Matrix([[nonfinite]])),(sx.Matrix([[nonfinite]]),sx.ZeroMatrix(1,1))]:
            for function in [sx.hadamard_product,sx.kronecker_product]:
                value=function(*operands)
                assert value.as_explicit()[0,0] is sx.nan or value.as_explicit()[0,0]==sx.nan


@pytest.mark.parametrize('operation',[
    lambda s,x,A,B:s.HadamardProduct(x*A,B),
    lambda s,x,A,B:s.HadamardPower(x*A+B,3),
    lambda s,x,A,B:s.HadamardPower(x*A+B,-2),
    lambda s,x,A,B:s.HadamardPower(A,x),
    lambda s,x,A,B:s.HadamardPower(2,x*A),
    lambda s,x,A,B:s.HadamardPower(x+2,A),
    lambda s,x,A,B:s.HadamardPower(x*A+B,x*B),
    lambda s,x,A,B:s.KroneckerProduct(x*A,B+x*A),
])
def test_scalar_derivatives_match_entrywise_sympy(operation):
    x=sx.Symbol('x',real=True);A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    X=sp.Symbol('x',real=True);a=sp.MatrixSymbol('A',2,2);b=sp.MatrixSymbol('B',2,2)
    expression=operation(sx,x,A,B)
    derivative=expression.diff(x)
    native={A:sx.Matrix([[1,2],[3,4]]),B:sx.Matrix([[2,3],[4,5]])}
    oracle={a:sp.ImmutableMatrix([[1,2],[3,4]]),b:sp.ImmutableMatrix([[2,3],[4,5]])}
    expected=operation(sp,X,a,b).as_explicit().diff(X).xreplace(oracle)
    same_entries(derivative.xreplace(native),expected)


def test_symbolic_variable_exponent_requires_explicit_entries():
    n=sx.Symbol('n',integer=True,positive=True);x=sx.Symbol('x',real=True);A=sx.MatrixSymbol('A',n,n)
    # SymPy 1.14: HadamardProduct(Lambda(_d, log(_d)).(A), HadamardPower(A, x))
    assert str(sx.HadamardPower(A,x).diff(x))=="HadamardProduct(Lambda(_d, log(_d)).(A), HadamardPower(A, x))"
    assert sx.HadamardPower(A,3).diff(x)==sx.ZeroMatrix(n,n)
    assert sx.HadamardPower(2,x*A).diff(x).shape==(n,n)


@pytest.mark.parametrize('mode',['hadamard','matrix_power','scalar_power','kronecker'])
def test_numpy_expression_source_and_tensor_source_generation(mode):
    A=sx.MatrixSymbol('A',2,2);B=sx.MatrixSymbol('B',2,2)
    values={'A':np.array([[1.,2.],[3.,4.]]),'B':np.array([[2.,3.],[4.,5.]])}
    if mode=='hadamard': expression=sx.hadamard_product(A,B); expected=values['A']*values['B']
    elif mode=='matrix_power': expression=sx.HadamardPower(A,B); expected=values['A']**values['B']
    elif mode=='scalar_power': expression=sx.HadamardPower(2,A); expected=2**values['A']
    else: expression=sx.kronecker_product(A,B); expected=np.kron(values['A'],values['B'])
    explicit=expression.as_explicit()
    rows,cols=explicit.shape
    source='numpy.array(['+','.join('['+','.join(sx.numpycode(explicit[i,j]) for j in range(cols))+']' for i in range(rows))+'])'
    actual=eval(source,{'numpy':np,**values})
    np.testing.assert_allclose(actual,expected)
    for generate in [sx.extra.torchcode_function,sx.extra.jaxcode_function]:
        source=generate('matrix_product',(A,B),expression)
        ast.parse(source)
        assert 'Hadamard' not in source and 'Kronecker' not in source
        assert 'symox' not in source and 'sympy' not in source


def test_import_paths_and_native_runtime_without_sympy():
    from symox.matrices.expressions.hadamard import HadamardProduct,HadamardPower
    from symox.matrices.expressions.kronecker import KroneckerProduct
    from symox.matrices.expressions.special import OneMatrix
    assert (HadamardProduct,HadamardPower,KroneckerProduct,OneMatrix)==(sx.HadamardProduct,sx.HadamardPower,sx.KroneckerProduct,sx.OneMatrix)
    source='''
import sys
import symox as s
A=s.MatrixSymbol('A',2,3); B=s.MatrixSymbol('B',2,3)
assert s.hadamard_product(A,s.OneMatrix(2,3))==A
assert s.HadamardPower(2,A)[0,0]==2**A[0,0]
assert s.kronecker_product(A,B).shape==(4,9)
assert s.kronecker_product(A,B)[0,4]==A[0,1]*B[0,1]
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,capture_output=True,text=True)


def test_explicit_entries_keep_matrix_metadata_substitutions():
    A=sx.MatrixSymbol('A',2,2)
    explicit=sx.MatrixExpr(A.as_explicit())
    values=sx.ImmutableMatrix([[1,2],[3,4]])
    assert explicit.xreplace({A:values}).as_explicit()==values
    assert explicit.subs(A,values).as_explicit()==values
    n=sx.Symbol('n',integer=True,positive=True)
    X=sx.MatrixSymbol('X',n,n);Y=sx.MatrixSymbol('Y',n,n)
    entry_matrix=sx.MatrixExpr(sx.Matrix([[X[0,0]]]))
    replaced=entry_matrix.xreplace({X:Y,n:3})
    assert replaced.as_explicit()[0,0]==Y[0,0]
