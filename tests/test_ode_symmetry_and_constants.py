"""Arbitrary constants and point symmetries have independent certificates."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from sympy.solvers.ode import constantsimp as reference_constantsimp
from sympy.solvers.ode import checkinfsol as reference_checkinfsol
from symox.solvers.ode import constantsimp, infinitesimals, checkinfsol, dsolve


def namespace(module):
    x,y = module.symbols('x y')
    f = module.Function('f')(x)
    return vars(module) | {'x':x,'y':y,'f':f,'F':module.Function('F'),
                           **dict(zip(['C1','C2','C3','a','b','m','n'],module.symbols('C1 C2 C3 a b m n')))}


def convert(value):
    if hasattr(value,'to_sympy'):
        return value.to_sympy()
    if isinstance(value,(list,tuple,set,frozenset)):
        return type(value)(convert(item) for item in value)
    if isinstance(value,dict):
        return {convert(key):convert(item) for key,item in value.items()}
    return value


@pytest.mark.parametrize('source',[
    '0','x','C1','C1+C2','C1*C2','2*C1*x','C1+2+x','C1*C2+2+C2+C3*x',
    '(C1*cos(x)+C2*cos(x))*exp(x)','C1*cos(x)+C2*cos(x)+C3*sin(x)',
    'exp(C1+x)','exp(C1)+2+(exp(C1)+2)*x',
    'C1+x+Integral(x,(x,1,2))','sin(C1)+x','log(C1)+x',
    '(C1+C2)*x+C3','C1*x+C1*x*x','C1*x+2*C1*x*x',
    'C1*x+C2*x+C2*x*x','(2*C1+x)/(C1+x)','C1*x+exp(C1)',
    'exp(C1)*x+exp(C1)*sin(x)','2*C1*sin(x)+C2*cos(x)',
    'Eq(C1+x,C2)','Eq(2*C1*x,C2)','Eq(C1*x,C2*x+C3)',
    'Eq(C1,x+2)','Eq(x+C1,0)','Eq(C1+C2*x,C3)',
    '[2*C1*x,C1+2+x]','(2*C1*x,C1+2+x)',
])
def test_arbitrary_constant_absorption_matches_reference(source):
    native,reference = namespace(sx),namespace(sp)
    actual = constantsimp(eval(source,native),{native['C1'],native['C2'],native['C3']})
    expected_input = eval(source,reference)
    if isinstance(expected_input,int):
        expected_input = sp.sympify(expected_input)
    expected = reference_constantsimp(expected_input,{reference['C1'],reference['C2'],reference['C3']})
    assert convert(actual)==expected


def test_unlisted_parameters_and_repeated_constants_keep_their_correlations():
    x,C1,a = sx.symbols('x C1 a')
    value = (2*C1+x)/(C1+x)
    assert constantsimp(value,{C1})==value
    assert a in constantsimp(C1*a+sx.exp(a),{C1}).free_symbols
    assert constantsimp(C1*x+sx.exp(C1),{C1})==C1*x+sx.exp(C1)
    # SymPy leaves an expression unchanged when a constant is not a symbol
    assert constantsimp(C1*x,{C1*x})==C1*x


@pytest.mark.parametrize('source,hint',[
    ('diff(f,x)-x*x*f','default'),
    ('diff(f,x)-x*x*f','abaco1_simple'),
    ('diff(f,x)+a*f-b*exp(x)','abaco1_simple'),
    ('diff(f,x)+2*x*f-x*exp(-x*x)','abaco1_simple'),
    ('(1+2*x)*diff(f,x)+2-4*exp(-f)','abaco1_simple'),
    ('x*x*diff(f,x)-f+x*x*exp(x-1/x)','abaco1_simple'),
    ('diff(f,x)-(f/x)*(x*log(x*x/f)+2)','abaco1_product'),
    ('diff(f,x)-F(a*x+b*f)','abaco2_similar'),
    ('diff(f,x)-f*f/(sin(f-x)-x*x+2*x*f)','abaco2_similar'),
    ('diff(f,x)-x**(a-1)*f**(1-b)*F(x**a/a+f**b/b)','abaco2_unique_unknown'),
    ('diff(f,x)+tan(F(x*x+f*f)+atan(x/f))','abaco2_unique_unknown'),
    ('diff(f,x)-exp(x+f)','abaco2_unique_general'),
    ('x**(n*(m+1)-m)*diff(f,x)-a*f**n-b*x**(n*(m+1))','linear'),
    ('diff(f,x)-(3*(1+x*x/f**2)*atan(f/x)+(1-2*f)/x+(1-3*f)*x/f**2)','function_sum'),
    ('x*x*diff(f,x)+x*f+f*f+x*x','bivariate'),
    ('x*x*(-f*f+diff(f,x))-a*x*x*f+2-a*x','bivariate'),
    ('x*diff(f,x)+1-f*f','chi'),
    ('diff(f,x)-f','all'),
])
def test_every_named_symmetry_search_returns_proved_generators(source,hint):
    native = namespace(sx)
    equation = eval(source,native)
    function,x = native['f'],native['x']
    found = infinitesimals(equation,function,hint=hint)
    assert found
    assert all(valid and residual==0 for valid,residual in checkinfsol(equation,found,function))
    reference = reference_checkinfsol(equation.to_sympy(),convert(found),function.to_sympy())
    # SymPy's checker can leave cancellable symbolic powers in its residual.
    # Demand an independent exact zero certificate, including that case.
    assert all(sp.simplify(sp.powsimp(sp.expand(residual)))==0
               for _,residual in reference)
    # A perturbation must be rejected by its actual determining residual.
    xi,eta = sx.Function('xi')(x,function),sx.Function('eta')(x,function)
    perturbed = dict(found[0])
    perturbed[eta] += x*function+x*x
    checked,residual = checkinfsol(equation,[perturbed],function)[0]
    assert not checked and residual!=0


def test_simple_symmetry_values_match_the_reference_basis():
    from sympy.solvers.ode import infinitesimals as reference
    native = namespace(sx)
    equation = native['f'].diff(native['x'])-native['x']**2*native['f']
    assert convert(infinitesimals(equation,native['f'],hint='abaco1_simple'))==reference(
        equation.to_sympy(),native['f'].to_sympy(),hint='abaco1_simple')


@pytest.mark.parametrize('source',[
    'diff(f,x)-x*x*f','diff(f,x)-x*f*f','x*diff(f,x)+1-f*f',
    'diff(f,x)+x*f-x','diff(f,x)-1/(x+f)',
])
def test_lie_solver_constructs_explicit_or_implicit_first_integrals(source):
    native = namespace(sx)
    equation,function = eval(source,native),native['f']
    solution = dsolve(equation,function,hint='lie_group')
    solutions = solution if isinstance(solution,list) else [solution]
    for item in solutions:
        assert sx.checkodesol(equation,item,function)==(True,0)
        assert sp.checkodesol(equation.to_sympy(),item.to_sympy(),function.to_sympy())[0]


def test_custom_symmetry_is_verified_and_integrated():
    native = namespace(sx)
    x,f = native['x'],native['f']
    equation = x*f.diff(x)+1-f*f
    solution = dsolve(equation,f,hint='lie_group',xi=sx.sqrt(f-1)/sx.sqrt(f+1),eta=0)
    assert sx.checkodesol(equation,solution,f)==(True,0)
    # SymPy does not verify a supplied symmetry; these still give the
    # reference solution
    expected = sx.Eq(f,(-native['C1']-x**2)/(-native['C1']+x**2))
    assert dsolve(equation,f,hint='lie_group',xi=0,eta=1)==expected
    assert dsolve(equation,f,hint='lie_group',xi=1)==expected


@pytest.mark.parametrize('initial',[False,True])
def test_unknown_slope_preserves_held_integral_binding_and_initial_data(initial):
    native = namespace(sx)
    x,f = native['x'],native['f']
    equation = f.diff(x)-sx.Function('G')(f)
    conditions = {sx.Function('f')(0):sx.Integer(1)} if initial else None
    # SymPy's solve finds no explicit solution of the held quadrature, so the
    # Lie group method fails like the reference
    with pytest.raises(NotImplementedError,match='cannot be solved by the lie group method'):
        dsolve(equation,f,hint='lie_group',xi=1,eta=0,ics=conditions)


def test_lie_initial_conditions_preserve_external_constant_names():
    native = namespace(sx)
    x,f,C1 = native['x'],native['f'],native['C1']
    solution = dsolve(f.diff(x)-f,f,hint='lie_group',ics={sx.Function('f')(0):C1})
    assert solution==sx.Eq(f,C1*sx.exp(x))
    assert sx.checkodesol(f.diff(x)-f,solution,f)==(True,0)


def test_symmetry_argument_validation_and_held_order_limits():
    native = namespace(sx)
    x,f = native['x'],native['f']
    with pytest.raises(ValueError):
        infinitesimals(f.diff(x)-f,f,hint='unknown')
    with pytest.raises(NotImplementedError):
        infinitesimals(f.diff(x,2)+f,f)
    with pytest.raises(NotImplementedError):
        checkinfsol(sx.Derivative(f,x,10**9,evaluate=False),[],f)
    # SymPy indexes the generator by xi(x, f(x))
    with pytest.raises(KeyError):
        checkinfsol(f.diff(x)-f,[{}],f)


def test_ode_completion_has_no_sympy_runtime_dependency():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.solvers.ode import constantsimp,infinitesimals,checkinfsol
x,C1,C2=symbols('x C1 C2');f=Function('f')(x)
assert constantsimp((C1+C2)*x,{C1,C2})==C1*x
eq=f.diff(x)-x*x*f
assert all(valid and residual==0 for valid,residual in checkinfsol(eq,infinitesimals(eq,f),f))
assert checkodesol(eq,dsolve(eq,f,hint='lie_group'),f)==(True,0)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
