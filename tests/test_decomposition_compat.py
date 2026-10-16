"""Differential coverage for structural coefficient/content decompositions."""
import pickle
import subprocess
import sys
import pytest
import sympy as sp
import symox as sm


def expressions(module):
    x,y,z=module.symbols('x y z')
    f=module.Function('f')
    return x,y,z,f


def convert(value):
    if isinstance(value,tuple):return tuple(convert(v) for v in value)
    return sm.to_sympy(value)


@pytest.mark.parametrize('source',[
    '0','Rational(-2,3)','Float(2)','Float(-2)','oo','-oo','zoo','nan','I','pi',
    'sqrt(2)','x','3*x+y+2','6*x*y','x*(y+1)','x+y','x*y+z','sin(x)',
    'Subs(f(x),x,0)','Limit(f(x),x,0)',
    'Mul(2,Float(3),x,evaluate=False)','Add(x,2,evaluate=False)',
])
@pytest.mark.parametrize('dep_source',['none','x','y','x*y'])
def test_coefficient_tuples_match_oracle(source,dep_source):
    scopes=[]
    for m in [sm,sp]:
        x,y,z,f=expressions(m)
        scope=dict(vars(m),x=x,y=y,z=z,f=f)
        e=eval(source,scope)
        if isinstance(e,int):e=m.sympify(e)
        deps=() if dep_source=='none' else (eval(dep_source,scope),)
        scopes.append((e,deps))
    (actual,deps),(expected,spdeps)=scopes
    assert convert(actual.as_coeff_add(*deps))==expected.as_coeff_add(*spdeps)
    for rational in [False,True]:
        assert convert(actual.as_coeff_mul(*deps,rational=rational))==expected.as_coeff_mul(*spdeps,rational=rational)


@pytest.mark.parametrize('source',['0','Float(0)','x','x+y','x*y','x*y+z','x*(y+1)','sin(x)','Subs(f(x),x,0)','Limit(f(x)+y,x,0)'])
@pytest.mark.parametrize('hint',[None,False,True])
def test_independent_split_uses_free_symbol_scope(source,hint):
    scopes=[]
    for m in [sm,sp]:
        x,y,z,f=expressions(m);scope=dict(vars(m),x=x,y=y,z=z,f=f)
        e=eval(source,scope);e=m.sympify(e)
        scopes.append((e,x))
    (actual,x),(expected,sx)=scopes
    kwargs={} if hint is None else {'as_Add':hint}
    assert convert(actual.as_independent(x,**kwargs))==expected.as_independent(sx,**kwargs)


@pytest.mark.parametrize('dep_source',['x+y','x*y','x+1','sin(x)'])
def test_independent_part_of_held_add_is_reconstructed_canonically(dep_source):
    values=[]
    for m in [sm,sp]:
        x,y,z,f=expressions(m)
        dep=eval(dep_source,dict(vars(m),x=x,y=y,z=z,f=f))
        values.append(m.Add(x,2,evaluate=False).as_independent(dep))
    assert convert(values[0])==values[1]
    assert values[0][0].args==(2,sm.Symbol('x'))


@pytest.mark.parametrize('source',[
    'Rational(-2,3)','Float(2)','Float(-2)','0','I','oo','sin(2*x+2)',
    '2+2*x+2*y*(3+3*y)','(2+6*x)**2','(2+6*x)**(2*y)',
    '5*x*(1+y)+2*x*(3+3*y)','3*z*(1+y)+Float(2)*x*(3+3*y)',
    'x/2+y','2*x/3+4*y/9','2*x/3+4*y/3',
    'sqrt(4+4*sqrt(2))','sqrt(3+3*sqrt(2))','2*sqrt(2)+4*sqrt(10)',
    'Pow(3,(5+y)/2,evaluate=False)','Pow(4,(1+y)/2,evaluate=False)',
    'Pow(3,2+2*x,evaluate=False)','Pow(-2,x+3,evaluate=False)',
    'Pow(Rational(1,2),x+2,evaluate=False)','Pow(2,x-Rational(1,2),evaluate=False)',
])
@pytest.mark.parametrize('radical,clear',[(False,True),(False,False),(True,True),(True,False)])
def test_recursive_rational_content_matches_oracle(source,radical,clear):
    values=[]
    for m in [sm,sp]:
        x,y,z,f=expressions(m);e=m.sympify(eval(source,dict(vars(m),x=x,y=y,z=z,f=f)))
        values.append(e.as_content_primitive(radical=radical,clear=clear))
    assert convert(values[0])==values[1]
    content,primitive=values[0]
    assert pickle.loads(pickle.dumps(primitive))==primitive
    assert sm.from_sympy(sm.to_sympy(primitive))==primitive


def test_options_are_explicit_and_preserved_precision_is_not_numeric_equality():
    x=sm.Symbol('x');low=sm.Float('1.25');high=sm.Float('1.25',80)
    assert low._prec<high._prec and float(low)==float(high)
    f=sm.Function('pair')
    for a,b in [(low,high),(high,low)]:
        held=f(a,b)
        assert held.args[0]._prec==a._prec and held.args[1]._prec==b._prec
    coefficient,tail=(-high*x).as_coeff_mul()
    assert coefficient==-1 and tail[0]._prec==high._prec
    with pytest.raises(NotImplementedError):x.as_coeff_mul(unknown=True)
    with pytest.raises(NotImplementedError):x.as_independent(x,unknown=True)
    with pytest.raises(NotImplementedError):sm.Pow(2,600_000,evaluate=False).as_content_primitive()


def test_dag_function_construction_and_decomposition_have_no_sympy_dependency():
    subprocess.run([sys.executable,'-c',"""
import sys
sys.modules['sympy']=None
import symox as s
x,y=s.symbols('x y');pair=s.Function('pair')
dag=x
for _ in range(40):dag=pair(dag,dag)
assert dag.func==pair
assert (3*x+y+2).as_coeff_add(x)==(y+2,(3*x,))
assert (6*x*y).as_coeff_mul(x)==(6*y,(x,))
assert (x*y+y).as_independent(x)==(y,x*y)
assert (x/2+y).as_content_primitive(clear=False)==(1,x/2+y)
"""],check=True,timeout=20)
