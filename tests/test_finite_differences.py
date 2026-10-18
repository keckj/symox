"""Fornberg coefficients, exact polynomial moments and shifted derivatives."""
import itertools
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox import finite_diff_weights,apply_finite_diff,differentiate_finite


def equivalent(actual,expected):
    assert sp.cancel(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('points', [[0],[0,1],[-2,-1,0,1,3],
    [sx.Rational(1,3),sx.Rational(-2,5),sx.Rational(3,7)], [0,sx.I,1+sx.I]])
@pytest.mark.parametrize('order',[0,1,2,4])
@pytest.mark.parametrize('center',[0,1,sx.Rational(2,3)])
def test_weight_tables_match_reference(points,order,center):
    native=finite_diff_weights(order,points,center)
    reference=sp.finite_diff_weights(order,[sx.sympify(p).to_sympy() for p in points],
                                     sx.sympify(center).to_sympy())
    for a,b in zip(native,reference):
        for aa,bb in zip(a,b):
            for value,expected in zip(aa,bb):
                equivalent(value,expected)


def test_every_prefix_reproduces_exact_polynomial_derivatives():
    rng=random.Random(173)
    grid=rng.sample(range(-20,21),8)
    center=sx.Rational(3,7)
    weights=finite_diff_weights(5,grid,center)
    for derivative in range(6):
        for prefix in range(1,len(grid)+1):
            for power in range(prefix):
                actual=sum(weights[derivative][prefix-1][j]*grid[j]**power for j in range(prefix))
                expected=(sx.factorial(power)/sx.factorial(power-derivative)*center**(power-derivative)
                          if derivative<=power else sx.S.Zero)
                assert sx.cancel(actual-expected)==0


def test_symbolic_irregular_grid_and_interpolation():
    x,h=sx.symbols('x h');X,H=sp.symbols('x h')
    grid=[x-h,x,x+2*h,x+3*h]
    actual=finite_diff_weights(3,grid,x+h/3)
    expected=sp.finite_diff_weights(3,[X-H,X,X+2*H,X+3*H],X+H/3)
    for order,prefix,position in itertools.product(range(4),repeat=3):
        equivalent(actual[order][prefix][position],expected[order][prefix][position])
    assert sx.cancel(apply_finite_diff(0,grid,[p**3 for p in grid],x+h/3)-(x+h/3)**3)==0
    assert sx.cancel(apply_finite_diff(2,grid,(p**3 for p in grid),x)-6*x)==0


@pytest.mark.parametrize('order',[1,2,3,4])
@pytest.mark.parametrize('kind',['unit','step','grid','offset','variable_step'])
def test_derivative_stencils_match_reference(order,kind):
    x,h=sx.symbols('x h');X,H=sp.symbols('x h');f=sx.Function('f');F=sp.Function('f')
    dx=sx.Function('dx');DX=sp.Function('dx')
    parameters={'unit':{},'step':dict(points=h),'grid':dict(points=[x+i*h for i in range(order+2)]),
        'offset':dict(points=h,x0=x+3*h),'variable_step':dict(points=dx(x),x0=x+h)}[kind]
    reference={'unit':{},'step':dict(points=H),'grid':dict(points=[X+i*H for i in range(order+2)]),
        'offset':dict(points=H,x0=X+3*H),'variable_step':dict(points=DX(X),x0=X+H)}[kind]
    equivalent(sx.Derivative(f(x),x,order).as_finite_difference(**parameters),
               sp.Derivative(F(X),X,order).as_finite_difference(**reference))


@pytest.mark.parametrize('kind',['product','mixed','nested','variable_step','evaluated'])
def test_differentiate_finite_matches_reference(kind):
    x,y,h=sx.symbols('x y h');X,Y,H=sp.symbols('x y h')
    f,g,dx=[sx.Function(n) for n in ['f','g','dx']]
    F,G,DX=[sp.Function(n) for n in ['f','g','dx']]
    native={'product':(f(x)*g(x),(x,),dict(points=[x-h,x+h])),
        'mixed':(f(x,y),(x,y),{}),
        'nested':(f(x)*g(x).diff(x),(x,),{}),
        'variable_step':(f(x)*g(x).diff(x),(),dict(points=dx(x))),
        'evaluated':(f(x)*g(x),(x,),dict(evaluate=True))}[kind]
    reference={'product':(F(X)*G(X),(X,),dict(points=[X-H,X+H])),
        'mixed':(F(X,Y),(X,Y),{}),
        'nested':(F(X)*G(X).diff(X),(X,),{}),
        'variable_step':(F(X)*G(X).diff(X),(),dict(points=DX(X))),
        'evaluated':(F(X)*G(X),(X,),dict(evaluate=True))}[kind]
    with sp.utilities.exceptions.ignore_warnings(sp.utilities.exceptions.SymPyDeprecationWarning):
        expected=sp.differentiate_finite(reference[0],*reference[1],**reference[2])
    equivalent(differentiate_finite(native[0],*native[1],**native[2]),expected)


def test_partial_derivative_preserves_other_variables():
    x,y,h=sx.symbols('x y h');X,Y,H=sp.symbols('x y h')
    f=sx.Function('f');F=sp.Function('f')
    value=sx.Derivative(f(x,y),x,2,y)
    reference=sp.Derivative(F(X,Y),X,2,Y)
    equivalent(value.as_finite_difference(h,wrt=x),reference.as_finite_difference(H,wrt=X))
    equivalent(value.as_finite_difference(h),reference.as_finite_difference(H))


@pytest.mark.parametrize('power',[1,2,3])
def test_nested_mixed_derivative_stencils(power):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y');f=sx.Function('f');F=sp.Function('f')
    expr=sx.Derivative(f(x,y),x,y)**power
    actual=differentiate_finite(expr,x,y)
    # SymPy 1.14 incorrectly returns zero for powers > 1 here. Compare
    # against the definition: apply the two centered differences to the
    # power of the inner centered mixed difference.
    half=sp.Rational(1,2)
    def inner(a,b):
        return F(a-half,b-half)-F(a-half,b+half)-F(a+half,b-half)+F(a+half,b+half)
    expected=sum(i*j*inner(X+i*half,Y+j*half)**power for i,j in itertools.product([-1,1],repeat=2))
    equivalent(actual,expected)
    pending=[actual]
    while pending:
        node=pending.pop()
        assert not node.is_Derivative and not node.is_Subs
        pending.extend(child for child in node.args if isinstance(child,sx.Expr))


def test_invalid_grids_orders_and_resource_limits():
    x=sx.Symbol('x');f=sx.Function('f')
    for order in [-1,sx.Rational(1,2),x,sx.oo,1000]:
        with pytest.raises((ValueError,TypeError)):
            finite_diff_weights(order,[0,1])
    for grid in [[],[0,0],[x,x],[0,sx.oo],range(1000)]:
        with pytest.raises(ValueError):
            finite_diff_weights(1,grid)
    with pytest.raises(ValueError):
        finite_diff_weights(50,range(100))
    with pytest.raises(ValueError):
        apply_finite_diff(1,[0,1],[1])
    with pytest.raises(ValueError):
        sx.Derivative(f(x),x).as_finite_difference(0)
    with pytest.raises(ValueError):
        sx.Derivative(f(x),x,3).as_finite_difference([0,1])


def test_integral_floating_order_matches_reference():
    assert finite_diff_weights(2.0,[-1,0,1],0)==finite_diff_weights(2,[-1,0,1],0)


def test_without_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols,Function
from symox.calculus.finite_diff import finite_diff_weights,apply_finite_diff,differentiate_finite
x=symbols('x');f=Function('f')
assert finite_diff_weights(2,[-1,0,1],0)[2][-1]==[1,-2,1]
assert apply_finite_diff(1,[-1,0,1],[1,0,1])==0
assert differentiate_finite(f(x),x,2)==f(x-1)-2*f(x)+f(x+1)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
