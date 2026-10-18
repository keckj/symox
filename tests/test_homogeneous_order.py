"""Homogeneous degrees follow positive scaling and principal branches."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def expression(module,source):
    x,y,n,a=module.symbols('x y n a');f=module.Function('f')
    return module.sympify(eval(source,vars(module)|{'x':x,'y':y,'n':n,'a':a,'f':f}))


@pytest.mark.parametrize('source',[
    '0','3','pi','I','a','a*x*x','x+y','x*x+y*y','x*x+y',
    'x*y/(x+y)','sqrt(x*x+y*y)','x*x/sqrt(x*x+y*y)',
    'sqrt(x*y)','x**n','sqrt(x**I)',
    'sin(y/x)','cos(x/y)','exp(y/x)','log(y/x)',
    'sin(x)','exp(x+y)','x*x+sin(y)',
    'f(x/y)','f(x)','Derivative(f(x),x,evaluate=False)',
    'Order(x*x,x)','(x+y)**3/(x*x+y*y)',
    '1/(x*x+y*y)','x**2*f(x)/sqrt(x*x+f(x)**2)',
])
def test_homogeneous_scaling_degrees_match_reference(source):
    x,y=sx.symbols('x y');f=sx.Function('f')
    X,Y=sp.symbols('x y');F=sp.Function('f')
    actual=sx.homogeneous_order(expression(sx,source),x,y,f(x))
    expected=sp.homogeneous_order(expression(sp,source),X,Y,F(X))
    assert (actual.to_sympy() if actual is not None else None)==expected


def test_function_coordinates_require_their_arguments_and_preserve_order():
    x,y=sx.symbols('x y');f=sx.Function('f')
    assert sx.homogeneous_order(f(x),f(x)) is None
    assert sx.homogeneous_order(f(x),f(x),x)==1
    assert sx.homogeneous_order(f(x,y),f(y,x),x,y) is None
    assert sx.homogeneous_order(x*x+f(x),x,f(x)) is None
    assert sx.homogeneous_order(x*x*f(x)/sx.sqrt(x*x+f(x)**2),x,f(x))==2
    assert sx.homogeneous_order(x*x*f(y),x,f(x))==2


@pytest.mark.parametrize('source,degree',[
    ('Abs(x)',1),('Abs(x*x+y*y)',2),('Abs(x**I)',0),
    ('conjugate(x*x+y*y)',2),('re(x)',1),('im(x)',1),
    ('sign(x)',0),('arg(x)',0),
])
def test_positive_scaling_proves_absolute_values_and_complex_projections(source,degree):
    x,y=sx.symbols('x y');value=expression(sx,source)
    assert sx.homogeneous_order(value,x,y)==degree
    # These direct scaling witnesses supplement the structural proof and
    # cover complex coordinates where a real-only shortcut would be wrong.
    for scale in [sx.Rational(1,3),sx.Integer(2),sx.Integer(5)]:
        point={x:2+sx.I,y:1-2*sx.I}
        original=value.subs(point)
        scaled=value.subs({x:scale*point[x],y:scale*point[y]})
        assert abs(complex((scaled-scale**degree*original).evalf(40)))<1e-12


def test_complex_inner_degrees_do_not_cross_principal_power_branches():
    x,y=sx.symbols('x y');n=sx.Symbol('n')
    assert sx.homogeneous_order(x**n,x)==n
    assert sx.homogeneous_order(x**(n+1)*y**(-n),x,y)==1
    assert sx.homogeneous_order(sx.sqrt(x**sx.I),x) is None
    assert sx.homogeneous_order((x**sx.I)**2,x)==2*sx.I


def test_large_degrees_and_repeated_graphs_do_not_expand_polynomials():
    x,y=sx.symbols('x y');pair=sx.Function('pair')
    degree=sx.Integer(10**100)
    assert sx.homogeneous_order(x**degree,x)==degree
    value=x/y
    for _ in range(80):value=pair(value,value)
    assert sx.homogeneous_order(value,x,y)==0


def test_missing_coordinates_and_nonscalar_inputs_raise():
    x=sx.Symbol('x')
    with pytest.raises(ValueError):sx.homogeneous_order(x)
    with pytest.raises(TypeError):sx.homogeneous_order(sx.Matrix([[x]]),x)
    with pytest.raises(TypeError):sx.homogeneous_order(sx.Eq(x,0),x)


def test_homogeneous_degrees_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y,n=symbols('x y n');f=Function('f')
assert homogeneous_order(x*x/sqrt(x*x+y*y),x,y)==1
assert homogeneous_order(x**n,x)==n
assert homogeneous_order(f(x),f(x),x)==1
assert homogeneous_order(sqrt(x**I),x) is None
'''
    subprocess.run([sys.executable,'-c',source],check=True)
