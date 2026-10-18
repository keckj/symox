"""Equation validation and scalar linear isolation, including excluded poles."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox import checksol,solve_linear,solve_linear_system_LU,solve_undetermined_coeffs


def test_standard_solver_imports():
    from symox.solvers.solvers import checksol as check
    from symox.solvers.solveset import solveset
    assert check is checksol and solveset is sx.solveset


@pytest.mark.parametrize('source,point,expected',[
    ('x*x-2','sqrt(2)',True),('x*x-2','1',False),
    ('sin(x)','pi',True),('sin(x)','pi/2',False),
    ('exp(x)-1','0',True),('log(x)','1',True),('log(x)','0',False),
    ('1/(x-1)','1',False),('sqrt(x)-2','4',True),('sqrt(x)+2','4',False),
    ('x*x+1','I',True),('x*x-1','-1',True),('x*x-1','2',False),
])
@pytest.mark.parametrize('numerical',[False,True])
def test_checksol_exact_roots_poles_and_principal_branches(source,point,expected,numerical):
    x=sx.Symbol('x');X=sp.Symbol('x')
    assert checksol(sx.sympify(source),x,sx.sympify(point),numerical=numerical) is expected
    assert sp.checksol(sp.sympify(source),X,sp.sympify(point),numerical=numerical) is expected


def test_checksol_equation_collections_and_mappings_are_not_mutated():
    x,y=sx.symbols('x y');solution={x:1,y:2}
    equations=[sx.Eq(x+y,3),sx.Eq(2*x-y,0)]
    assert checksol(equations,solution) is True
    assert solution=={x:1,y:2}
    assert checksol(sx.Matrix([x+y-3,2*x-y]),solution) is True
    assert checksol([x+y-3,x-2],solution) is False
    assert checksol([],solution) is True
    assert checksol(sx.true,{}) is True and checksol(sx.false,{}) is False
    assert checksol(x+y,{x:y,y:0}) is True
    assert checksol(sx.sin(x)**2+sx.cos(x)**2-1,{}) is True


def test_checksol_numerical_tolerance_is_explicit_and_optional():
    x=sx.Symbol('x')
    approximate=sx.Float(2**.5)
    assert checksol(x*x-2,x,approximate) is True
    assert checksol(x*x-2,x,approximate,numerical=False) is False
    assert checksol(x,x,sx.Rational(1,10**20)) is False
    assert checksol(x,x,sx.Float('1e-12')) is True
    assert checksol(x,x,sx.Float('1e-12'),numerical=False) is False
    assert checksol(x,x,sx.oo) is False
    assert checksol(x,x,sx.nan) is False


def test_checksol_unknown_results_options_and_warnings():
    x,y=sx.symbols('x y');f=sx.Function('f')
    assert checksol(f(y),{x:1},minimal=True,numerical=False) is None
    with pytest.warns(RuntimeWarning,match='Could not verify'):
        assert checksol(f(y),{x:1},warn=True) is None
    assert checksol(f(x),{f(x):0}) is True
    with pytest.raises(ValueError):checksol(x,x)
    with pytest.raises(TypeError):checksol(x,{x:0},unknown=True)


@pytest.mark.parametrize('source',[
    'x+1','x+y+z','2*x+y-1','x*x-1','x*y+1','x*x+y',
    '1/(1/x-2)','1/x-y','1/x+1/y-1','x/(x-1)',
    'sin(x)+y','x*x-x*(x+1)','x**2*(1/x-y**2/x)',
])
def test_scalar_linear_isolation_matches_reference_or_rational_fallback(source):
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    native=solve_linear(sx.sympify(source));reference=sp.solve_linear(sp.sympify(source))
    a,b=(value.to_sympy() for value in native)
    if reference[0] in (0,X,Y,Z):
        assert a==reference[0]
        assert sp.simplify(b-reference[1])==0
    else:
        # Different uncancelled common denominators represent the same
        # rational equation. All original poles are checked separately.
        assert sp.cancel(a/b-reference[0]/reference[1])==0
    if a in (X,Y,Z):
        assert sp.checksol(sp.sympify(source),a,b) is True


def test_scalar_linear_options_equalities_and_non_symbol_rejection():
    x,y,z=sx.symbols('x y z')
    assert solve_linear(x+y+z,symbols=[z])==(z,-x-y)
    assert solve_linear(x+y+z,exclude=[x])==(y,-x-z)
    assert solve_linear(sx.Eq(2*x+1,3))==(x,1)
    assert solve_linear(x+1,rhs=3)==(x,2)
    assert solve_linear(x-x)==(0,1)
    assert solve_linear(x+1,symbols=[y])==(0,1)
    assert solve_linear(1/(1/x-2))==(0,0)
    with pytest.raises(ValueError):solve_linear(sx.Eq(x,1),rhs=2)
    with pytest.raises(ValueError):solve_linear(x+1,symbols=[sx.sin(x)])


def test_linear_isolation_evaluates_relevant_held_derivatives():
    x,y=sx.symbols('x y')
    expression=sx.Derivative(x*x,x)+y
    assert solve_linear(expression,symbols=[x])==(x,-y/2)


@pytest.mark.parametrize('kind',[sx.Matrix,sx.ImmutableMatrix,sx.SparseMatrix])
def test_lu_system_solution_reconstructs_every_equation(kind):
    x,y,z=sx.symbols('x y z')
    augmented=kind([[1,2,0,1],[3,2,2,1],[2,0,0,1]])
    original=augmented.copy()
    solution=solve_linear_system_LU(augmented,[x,y,z])
    assert solution=={x:sx.Rational(1,2),y:sx.Rational(1,4),z:-sx.Rational(1,2)}
    assert augmented==original
    assert augmented[:,:-1]*sx.Matrix([solution[x],solution[y],solution[z]])==sx.Matrix([1,1,1])


def test_lu_system_rejects_invalid_shapes_and_nonunique_systems():
    x,y=sx.symbols('x y')
    with pytest.raises(ValueError):solve_linear_system_LU(sx.eye(2),[x,y])
    with pytest.raises(ValueError):solve_linear_system_LU(sx.Matrix([[1,2,3]]),[x,y])
    with pytest.raises(ValueError):solve_linear_system_LU(sx.Matrix([[1,2]]),[])
    with pytest.raises(sx.NonInvertibleMatrixError) as error:
        solve_linear_system_LU(sx.Matrix([[1,1,2],[2,2,4]]),[x,y])
    assert 'invertible' in str(error.value).lower() or 'singular' in str(error.value).lower()


def test_native_solver_checks_without_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols,sqrt,Matrix
from symox import checksol,solve_linear,solve_linear_system_LU
x,y=symbols('x y')
assert checksol(x*x-2,x,sqrt(2)) is True
assert solve_linear(1/(1/x-2))==(0,0)
assert solve_linear_system_LU(Matrix([[1,1,3],[1,-1,1]]),[x,y])=={x:2,y:1}
'''
    subprocess.run([sys.executable,'-c',source],check=True)


@pytest.mark.parametrize('source',[
    'a*x+b-2*x-3','a*x*x+b*x+c-(x*x+2*x+3)',
    'a*a*x+b-x','a*exp(x)+b-exp(x)-y',
    'a*sin(x)+b*cos(x)-2*sin(x)-3*cos(x)',
    '(a*x+b)/(x+1)-2','a-2','(a+b)*x+b-c',
])
@pytest.mark.parametrize('options',[{},dict(dict=True),dict(set=True)])
def test_undetermined_coefficients_match_reference(source,options):
    x,y,a,b,c=sx.symbols('x y a b c');X,Y,A,B,C=sp.symbols('x y a b c')
    actual=solve_undetermined_coeffs(sx.sympify(source),[a,b,c],**options)
    expected=sp.solve_undetermined_coeffs(sp.sympify(source),[A,B,C],**options)
    def convert(value):
        if isinstance(value,dict):return {convert(k):convert(v) for k,v in value.items()}
        if isinstance(value,(list,tuple,set)):return type(value)(convert(v) for v in value)
        return value.to_sympy() if isinstance(value,sx.Expr) else value
    assert convert(actual)==expected


def test_undetermined_coefficients_explicit_variables_parameters_and_invalid_generators():
    x,y,a,b,c=sx.symbols('x y a b c')
    expression=a*x*x+b*x+c-(x*x+2*x+3)
    assert solve_undetermined_coeffs(expression,[a,b,c],x)=={a:1,b:2,c:3}
    assert solve_undetermined_coeffs(a*x+b-y,[a,b],[x])=={a:0,b:y}
    assert solve_undetermined_coeffs(a**x-2**x,[a],x) is None
    assert solve_undetermined_coeffs(x+1,[a],x)==[]
    with pytest.raises(ValueError):solve_undetermined_coeffs(x+1,[],x)
    with pytest.raises(ValueError):solve_undetermined_coeffs(x+1,[a*a],x)
