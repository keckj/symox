"""Non-Symbol targets use solve's formal algebraic-atom convention."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


def context(module):
    x=module.Symbol('x');f=module.Function('f')
    return dict(x=x,f=f,exp=module.exp,sin=module.sin,log=module.log,Derivative=module.Derivative)


def foreign(value):
    if isinstance(value,dict):return {foreign(k):foreign(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,set)):return type(value)(foreign(v) for v in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


def equivalent(actual,expected):
    if isinstance(expected,dict):
        assert actual.keys()==expected.keys()
        for key in expected:equivalent(actual[key],expected[key])
    elif isinstance(expected,(list,tuple)):
        assert isinstance(actual,type(expected)) and len(actual)==len(expected)
        for a,b in zip(actual,expected):equivalent(a,b)
    elif isinstance(expected,set):
        assert len(actual)==len(expected)
        remaining=list(actual)
        for wanted in expected:
            for i,got in enumerate(remaining):
                if all(sp.simplify(a-b)==0 for a,b in zip(got,wanted)):
                    remaining.pop(i);break
            else:raise AssertionError((actual,expected))
    else:assert sp.simplify(actual-expected)==0


@pytest.mark.parametrize('source,target',[
    ('f(x)-3','f(x)'),
    ('Derivative(f(x),x)-f(x)-x','f(x)'),
    ('Derivative(f(x),(x,2))+f(x)**2-1','f(x)'),
    ('x+exp(2*x)','exp(x)'),
    ('exp(2*x+3)-2','exp(x)'),
    ('exp(x)+exp(-x)-2','exp(x)'),
    ('x**4-1','x**2'),
    ('sin(x)**2-2*sin(x)+1','sin(x)'),
    ('log(x)**2-1','log(x)'),
    ('f(1)**2+1','f(1)'),
    ('[f(x)-3,x-2]','f(x)'),
])
@pytest.mark.parametrize('flags',[{},{'dict':True},{'set':True},{'check':False}])
def test_reference_atom_roots_labels_and_output_shapes(source,target,flags):
    native,reference=context(sx),context(sp)
    actual=sx.solve(eval(source,native),eval(target,native),**flags)
    expected=sp.solve(eval(source,reference),eval(target,reference),**flags)
    equivalent(foreign(actual),expected)


def test_reference_atom_convention_does_not_copy_analytic_target_ranges():
    x=sx.Symbol('x');f=sx.Function('r',real=True)
    assert sx.solve(sx.exp(x),sx.exp(x))==[0]
    assert sx.solve(f(x)**2+1,f(x))==[-sx.I,sx.I]


@pytest.mark.parametrize('flags',[{},{'dict':True},{'set':True},{'check':False}])
def test_affine_exponential_atom_keeps_both_values_when_reference_misses_them(flags):
    x=sx.Symbol('x');target=sx.exp(2*x+1)
    values=[-sx.sqrt(5),sx.sqrt(5)]
    expected=([target],{(value,) for value in values}) if flags.get('set') else (
        [{target:value} for value in values] if flags.get('dict') else values)
    assert sx.solve(sx.exp(4*x+2)-5,target,**flags)==expected
    # exp(4*x+2)=target² for every complex x. Factoring the quadratic proves
    # these are the ONLY possible values, and exact preimages witness both;
    # this is a complete value-set proof, not sampled equation verification.
    t=sp.Symbol('t')
    assert sp.expand((t-sp.sqrt(5))*(t+sp.sqrt(5)))==t*t-5
    first=sp.log(5)/4-sp.Rational(1,2)
    second=first+sp.I*sp.pi/2
    assert sp.simplify(sp.exp(2*first+1)-sp.sqrt(5))==0
    assert sp.simplify(sp.exp(2*second+1)+sp.sqrt(5))==0


@pytest.mark.parametrize('source,options',[
    ('f(x)**5-f(x)-1',{}),('f(x)**3-f(x)+1',{'cubics':False}),
])
def test_complete_indexed_roots_survive_atom_restoration_without_dummy_leaks(source,options):
    native,reference=context(sx),context(sp)
    target=native['f'](native['x']);other=reference['f'](reference['x'])
    actual=sx.solve(eval(source,native),target,**options)
    expected=sp.solve(eval(source,reference),other,**options)
    assert len(actual)==len(expected)
    assert all(not value.free_symbols for value in actual)
    assert set(foreign(actual))==set(expected)


def test_private_dummy_names_cannot_capture_user_symbols():
    x=sx.Symbol('x');user=sx.Symbol('_solve_target');f=sx.Function('f')(x)
    assert sx.solve(f-user,f)==[user]
    derivative=sx.Derivative(f,x)
    assert sx.solve(f-derivative-user,f)==[derivative+user]


def test_native_atom_solver_never_imports_reference():
    script='''
import sys
class NoReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as sx
x=sx.Symbol('x');f=sx.Function('f')(x)
assert sx.solve(sx.Derivative(f,x)-f-x,f)==[sx.Derivative(f,x)-x]
assert sx.solve(x**4-1,x**2,set=True)==([x**2],{(-1,),(1,)})
assert sx.solve(x+sx.exp(2*x),sx.exp(x))==[-sx.sqrt(-x),sx.sqrt(-x)]
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=30)
