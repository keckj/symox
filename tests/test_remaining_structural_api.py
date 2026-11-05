"""Behavioral checks for the remaining structural array, matrix and set APIs."""
import itertools
import math
import subprocess
import sys
import numpy as np
import pytest
import symox as S
import sympy as R


def test_geometry_hierarchy_and_affine_transformations():
    from symox.geometry.entity import GeometrySet,rotate,scale,translate
    from symox.geometry.line import LinearEntity2D,LinearEntity3D
    from symox.geometry.util import are_coplanar,find
    for entity in [S.Line((0,0),(1,1)),S.Segment((0,0),(1,2)),S.Circle((0,0),2),S.Triangle((0,0),(1,0),(0,1))]:
        assert isinstance(entity,(GeometrySet,S.Set))
    assert isinstance(S.Line((0,0),(1,1)),LinearEntity2D)
    assert isinstance(S.Line3D((0,0,0),(1,2,3)),LinearEntity3D)
    for a in [0,S.pi/2,S.pi/3]:
        assert S.Point(2,1).transform(rotate(a))==S.Point(2,1).rotate(a)
    assert S.Point(2,3).transform(scale(2,4,(1,1)))==S.Point(3,9)
    assert S.Point(2,3).transform(translate(1,-2))==S.Point(3,1)
    assert are_coplanar(S.Plane((0,0,0),normal_vector=(0,0,1)),S.Line((0,0),(1,1)))
    assert not are_coplanar(S.Line3D((0,0,0),(1,0,0)),S.Line3D((0,1,1),(0,1,2)))
    x=S.Symbol('x');dummy=S.Dummy('x')
    assert find('x',dummy)==dummy
    with pytest.raises(ValueError):find('x',x+dummy)


def test_matrix_domain_companion_and_block_inverse():
    from symox.matrices.expressions import CompanionMatrix,MatrixSet
    from symox.matrices.kind import MatrixKind
    from symox.matrices.expressions.blockmatrix import blockinverse_2x2,reblock_2x2,deblock,bounds
    from symox.matrices.immutable import ImmutableRepMatrix,sympify_mpmath_matrix
    from symox.matrices.sparse import SparseRepMatrix
    from symox.matrices.matrixbase import a2idx
    x=S.Symbol('x');m=CompanionMatrix(S.Poly(x**3+2*x+3,x))
    assert S.expand(m.charpoly(x).as_expr()-(x**3+2*x+3))==0
    assert MatrixKind() is MatrixKind()
    assert MatrixSet(2,2,S.Reals).contains(S.eye(2))==S.true
    assert MatrixSet(2,2,S.Reals).contains(S.Matrix([[1,S.I],[0,1]]))==S.false
    assert MatrixSet(2,2,S.Reals).contains(S.ones(3,3))==S.false
    assert isinstance(ImmutableRepMatrix([[1]]),ImmutableRepMatrix)
    assert isinstance(SparseRepMatrix([[1]]),SparseRepMatrix)
    assert a2idx(-1,4)==3
    with pytest.raises(IndexError):a2idx(4,4)
    concrete=S.Matrix([[2,1],[3,4]])
    block=S.BlockMatrix([[S.Matrix([[2]]),S.Matrix([[1]])],[S.Matrix([[3]]),S.Matrix([[4]])]])
    assert blockinverse_2x2(S.Inverse(block)).as_explicit()==concrete.inv()
    blocks=S.BlockMatrix([[S.Matrix([[3*i+j]]) for j in range(3)] for i in range(3)])
    assert deblock(reblock_2x2(blocks)).as_explicit()==blocks.as_explicit()
    assert bounds([1,10,50])==[(0,1),(1,11),(11,61)]


@pytest.mark.parametrize('lo,hi,left,right',[(-1,1,False,False),(0,2,True,True),(3,5,False,True),(-4,3,True,False),(1,2,False,False)])
def test_complex_angle_normalization(lo,hi,left,right):
    from symox.sets.fancysets import normalize_theta_set,PolarComplexRegion,CartesianComplexRegion
    native=normalize_theta_set(S.Interval(lo*S.pi,hi*S.pi,left,right))
    reference=R.sets.fancysets.normalize_theta_set(R.Interval(lo*R.pi,hi*R.pi,left,right))
    for k in range(17):assert bool(native.contains(S.Rational(k,8)*S.pi))==bool(reference.contains(R.Rational(k,8)*R.pi))
    region=PolarComplexRegion(S.ProductSet(S.Interval(1,2),S.Interval(-S.pi/2,S.pi/2)))
    assert region.contains(S.I)==S.true and region.contains(-1)==S.false
    assert CartesianComplexRegion(S.ProductSet(S.Interval(0,1),S.Interval(0,2))).contains(1+S.I)==S.true


def test_elementwise_set_arithmetic():
    from symox.sets.sets import set_add,set_sub,set_mul,set_div,set_pow
    a=S.FiniteSet(1,2);b=S.FiniteSet(3,5)
    for fn,op in [(set_add,lambda x,y:x+y),(set_sub,lambda x,y:x-y),(set_mul,lambda x,y:x*y),(set_div,lambda x,y:x/y),(set_pow,lambda x,y:x**y)]:
        assert fn(a,b)==S.FiniteSet(*(op(x,y) for x in a for y in b))
    assert set_add(S.Interval.open(0,1),S.Interval(2,3))==S.Interval.open(2,4)
    assert set_mul(S.Interval.Lopen(1,3),-2)==S.Interval.Ropen(-6,-2)


def test_symbolic_array_operations_against_concrete_engine():
    from symox.tensor.array.expressions import ArrayTensorProduct,ArrayContraction,ArrayDiagonal,PermuteDims,Reshape,ArrayAdd,ArrayElementwiseApplyFunc,ArraySymbol,array_derive,convert_array_to_matrix,convert_matrix_to_array,convert_array_to_indexed,convert_indexed_to_array
    a=S.Array([[1,2],[3,4]]);b=S.Array([[5,6],[7,8]])
    tensor=ArrayTensorProduct(a,b)
    assert tensor.as_explicit()==S.tensorproduct(a,b)
    assert ArrayContraction(tensor,(1,2)).as_explicit()==S.Array(a.tomatrix()*b.tomatrix())
    assert ArrayDiagonal(tensor,(1,2)).as_explicit()==S.tensordiagonal(S.tensorproduct(a,b),(1,2))
    assert PermuteDims(tensor,(3,1,2,0)).as_explicit()==S.permutedims(S.tensorproduct(a,b),(3,1,2,0))
    assert Reshape(a,(4,)).as_explicit()==a.reshape(4)
    assert ArrayAdd(a,b).as_explicit()==a+b
    assert ArrayElementwiseApplyFunc(S.sin,a).as_explicit()==a.applyfunc(S.sin)
    A,B=ArraySymbol('A',(2,2)),ArraySymbol('B',(2,2))
    product=ArrayContraction(ArrayTensorProduct(A,B),(1,2))
    # SymPy 1.14: array symbols are not matrices, the contraction is kept
    assert convert_array_to_matrix(product)==product
    assert convert_array_to_matrix(convert_matrix_to_array(S.MatrixSymbol('A',2,2).T))==S.MatrixSymbol('A',2,2).T
    assert array_derive(A,A).as_explicit()==S.Array([S.Integer(int(i==k and j==l)) for i,j,k,l in itertools.product(range(2),repeat=4)],(2,2,2,2))
    i,j=S.symbols('i j',integer=True)
    assert convert_indexed_to_array(A[i,j])==A
    assert i in convert_array_to_indexed(product,[i,j]).free_symbols
    n=S.Symbol('n',integer=True,positive=True);X=ArraySymbol('X',(n,n))
    expression=ArrayContraction(X,(0,1))[()]
    assert expression.is_Sum and expression.function.has(expression.limits[0][0])


def test_tensor_partial_derivative_components():
    from symox.tensor.tensor import TensorIndexType,TensorHead,tensor_indices
    from symox.tensor.toperators import PartialDerivative
    L=TensorIndexType('L',dim=2);i,j=tensor_indices('i j',L);A=TensorHead('A',[L]);B=TensorHead('B',[L])
    x,y=S.symbols('x y');derivative=PartialDerivative(A(i),B(j))
    assert derivative.get_free_indices()==[i,-j]
    values=derivative.replace_with_arrays({A(i):[S.sin(x),S.log(x)*y**3],B(i):[x,y]})
    assert values==S.Array([[S.cos(x),0],[y**3/x,3*y**2*S.log(x)]])
    assert derivative.replace_with_arrays({A(i):[S.sin(x),S.log(x)*y**3],B(i):[x,y]},[-j,i])==S.permutedims(values,(1,0))
    divergence=PartialDerivative(A(i),B(i))
    assert divergence.get_free_indices()==[]
    assert divergence.replace_with_arrays({A(i):[x*x,x*y],B(i):[x,y]})==3*x


def test_compiled_numpy_ufunc_multiple_outputs(tmp_path):
    from symox.utilities.autowrap import ufuncify
    x,y=S.symbols('x y');function=ufuncify([x,y],[x*x+y,S.sin(x)-y],tempdir=tmp_path)
    assert isinstance(function,np.ufunc) and function.nin==2 and function.nout==2
    inputs=np.arange(6.).reshape(2,3);first,second=function(inputs,np.arange(3.))
    np.testing.assert_allclose(first,inputs**2+np.arange(3.))
    np.testing.assert_allclose(second,np.sin(inputs)-np.arange(3.))
    out=np.empty_like(inputs);assert function(inputs,1,out=(out,None))[0] is out


def test_f2py_autowrap(tmp_path):
    from symox.utilities.autowrap import autowrap,binary_function
    x=S.Symbol('x');function=autowrap(x*x+S.sin(x),tempdir=tmp_path)
    assert function(.25)==pytest.approx(.25**2+math.sin(.25))


def test_probability_domains_and_typed_spaces():
    from symox.stats import Normal,Die,E,P,pspace
    from symox.stats.rv import SinglePSpace,IndependentProductPSpace,RandomDomain,Density
    from symox.stats.crv import SingleContinuousDomain,ProductContinuousDomain,ConditionalContinuousDomain,ContinuousPSpace,SingleContinuousPSpace
    from symox.stats.frv import SingleFiniteDomain,ProductFiniteDomain,ConditionalFiniteDomain,FinitePSpace,SingleFinitePSpace,ProductFinitePSpace
    x,y=S.symbols('x y');X=Normal('remaining_normal',0,1);D=Die('remaining_die',6)
    assert isinstance(X.pspace,(SinglePSpace,SingleContinuousPSpace))
    assert isinstance(D.pspace,SingleFinitePSpace)
    assert isinstance(pspace(D+Die('remaining_die2',6)),ProductFinitePSpace)
    assert isinstance(X.pspace.domain,RandomDomain)
    assert X.pspace.compute_expectation(X**2)==1
    assert D.pspace.probability(D>3)==S.Rational(1,2)
    assert Density(D).doit().dict==D.distribution.dict
    a=SingleFiniteDomain(x,[1,2]);b=SingleFiniteDomain(y,[0,1]);product=ProductFiniteDomain(a,b)
    assert len(list(product))==4
    assert len(list(ConditionalFiniteDomain(product,x+y>2)))==1
    finite=FinitePSpace(a,{1:S.Rational(1,4),2:S.Rational(3,4)})
    assert finite.compute_expectation(x)==S.Rational(7,4)
    assert finite.probability(x>1)==S.Rational(3,4)
    domains=ProductContinuousDomain(SingleContinuousDomain(x,S.Interval(0,1)),SingleContinuousDomain(y,S.Interval(0,1)))
    continuous=ContinuousPSpace(domains,S.Lambda((x,y),4*x*y))
    assert continuous.compute_expectation(x*y,evaluate=True)==S.Rational(4,9)
    assert ContinuousPSpace(SingleContinuousDomain(x,S.Interval(0,1)),S.Lambda(x,2*x)).probability(x>S.Rational(1,2))==S.Rational(3,4)
    correlated=FinitePSpace(product,{frozenset([(x,1),(y,0)]):S.Rational(1,2),frozenset([(x,2),(y,1)]):S.Rational(1,2)})
    U,V=correlated.variables
    assert E(U*V)==1
    samples=correlated.sample(size=(20,),seed=42)
    assert all(float(v)==float(u)-1 for u,v in zip(samples[U],samples[V]))
