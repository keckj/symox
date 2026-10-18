"""PDE method entry points perform actual characteristic integrations."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.solvers import pde


CASES = [
    ('diff(u,x)+diff(u,y)','1st_linear_constant_coeff_homogeneous'),
    ('2*diff(u,x)+3*diff(u,y)+4*u','1st_linear_constant_coeff_homogeneous'),
    ('diff(u,x)+diff(u,y)-x','1st_linear_constant_coeff'),
    ('2*diff(u,x)-diff(u,y)+u-x*y','1st_linear_constant_coeff'),
    ('x*diff(u,x)-y*diff(u,y)+y*y*u-y*y','1st_linear_variable_coeff'),
    ('x*diff(u,x)+y*diff(u,y)+x*y*u-x*y','1st_linear_variable_coeff'),
    ('(x+y)*diff(u,x)+diff(u,y)','1st_linear_variable_coeff'),
    ('(y+x)*diff(u,y)+diff(u,x)','1st_linear_variable_coeff'),
    ('x*diff(u,x)+u-y','1st_linear_variable_coeff'),
    ('y*diff(u,y)+u-x','1st_linear_variable_coeff'),
]


def build(module,source):
    x,y = module.symbols('x y')
    u = module.Function('u')(x,y)
    return eval(source,vars(module)|{'x':x,'y':y,'u':u}),u


@pytest.mark.parametrize('source,hint',CASES)
@pytest.mark.parametrize('protocol',['native','wild'])
def test_specialized_methods_use_native_and_wild_match_protocols(source,hint,protocol):
    equation,function = build(sx,source)
    data = sx.classify_pde(equation,function,dict=True)[hint]
    if protocol=='wild':
        b,c,d,e = [sx.Wild(name) for name in ['b','c','d','e']]
        data = {'b':b,'c':c,'d':d,'e':e,b:data['coefficients'][0],
                c:data['coefficients'][1],d:data['coefficients'][2],e:-data['forcing']}
    method = getattr(pde,'pde_'+hint)
    solution = method(equation,function,1,data,sx.Function('H'))
    assert solution.has(sx.Function('H'))
    assert sx.checkpdesol(equation,solution,function)==(True,0)
    assert sp.checkpdesol(equation.to_sympy(),solution.to_sympy().doit(),function.to_sympy())[0]
    public = sx.pdsolve(equation,function,solvefun=sx.Function('H'))
    assert sx.checkpdesol(equation,public,function)==(True,0)


def test_constant_pde_method_retains_integral_until_evaluation():
    equation,function = build(sx,'diff(u,x)+diff(u,y)-x')
    data = sx.classify_pde(equation,function,dict=True)['1st_linear_constant_coeff']
    solution = pde.pde_1st_linear_constant_coeff(equation,function,1,data,sx.Function('H'))
    assert solution.has(sx.Integral)
    evaluated = solution.doit()
    assert not evaluated.has(sx.Integral)
    assert sx.checkpdesol(equation,evaluated,function)==(True,0)


@pytest.mark.parametrize('source',[
    'diff(u,x)+I*diff(u,y)','diff(u,x)+I*diff(u,y)+2*u-x',
    'I*diff(u,x)+diff(u,y)+u-exp(y)',
])
def test_isotropic_complex_characteristics_have_nondegenerate_coordinates(source):
    equation,function = build(sx,source)
    solution = sx.pdsolve(equation,function)
    assert sx.checkpdesol(equation,solution,function)==(True,0)
    assert sp.checkpdesol(equation.to_sympy(),solution.to_sympy().doit(),function.to_sympy())[0]


def test_unknown_forcing_remains_outside_characteristic_dummy_scope():
    x,y = sx.symbols('x y')
    u,G = sx.Function('u')(x,y),sx.Function('G')
    equation = x*u.diff(x)-y*u.diff(y)+y*y*u-G(x,y)
    solution = sx.pdsolve(equation,u)
    assert solution.has(sx.Integral)
    assert sx.checkpdesol(equation,solution,u)==(True,0)


def test_characteristic_constants_cannot_capture_parameters_in_forcing_only():
    x,y,C1,C2 = sx.symbols('x y C1 C2')
    u = sx.Function('u')(x,y)
    equation = x*u.diff(x)-y*u.diff(y)+C1*y*y*u-C2*y*y
    solution = sx.pdsolve(equation,u)
    assert {C1,C2} <= solution.free_symbols
    assert sx.checkpdesol(equation,solution,u)==(True,0)


def test_invalid_method_metadata_cannot_solve_a_different_equation():
    equation,function = build(sx,'diff(u,x)+diff(u,y)-x')
    data = sx.classify_pde(equation,function,dict=True)['1st_linear_constant_coeff']
    with pytest.raises(ValueError):
        pde.pde_1st_linear_constant_coeff(equation,function,2,data,sx.Function('F'))
    with pytest.raises(ValueError):
        pde.pde_1st_linear_constant_coeff(equation,function,1,data|{'forcing':sx.Integer(1)},sx.Function('F'))
    with pytest.raises(ValueError):
        pde.pde_1st_linear_constant_coeff_homogeneous(equation,function,1,data,sx.Function('F'))


def test_specialized_pde_methods_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.solvers import pde
x,y=symbols('x y');u=Function('u')(x,y)
eq=x*u.diff(x)-y*u.diff(y)+y*y*u-y*y
match=classify_pde(eq,u,dict=True)['1st_linear_variable_coeff']
solution=pde.pde_1st_linear_variable_coeff(eq,u,1,match,Function('H'))
assert checkpdesol(eq,solution,u)==(True,0)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
