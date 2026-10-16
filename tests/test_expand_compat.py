"""Independent expansion hints, bound-variable presentation and resource limits."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def expression(module,case):
    x,y,z,a,b=module.symbols('x y z a b')
    p,q=module.symbols('p q',positive=True)
    n=module.Symbol('n',nonzero=True)
    return [
        (x+y)*(x+z),z*(x+y)**2,module.sin((x+y)**2),
        module.exp(x*(y+z)),2*n**(a+b),2*module.exp(a+b),
        (x+y)**module.Rational(5,2),module.log(p*q),
        module.sin(x+y+z),(p*q)**a,
    ][case]


@pytest.mark.parametrize('case',range(10))
@pytest.mark.parametrize('flags',[
    {},{'deep':False},{'mul':False},{'multinomial':False},
    {'mul':False,'multinomial':False},{'power_exp':False},
    {'power_base':False},{'log':False},{'basic':False},
    {'trig':True},{'trig':True,'deep':False},{'force':True},
    {'trig':True,'mul':False},{'trig':True,'mul':False,'multinomial':False},
])
def test_expand_hints_match_sympy(case,flags):
    value=expression(sx,case)
    expected=expression(sp,case).expand(**flags)
    actual=value.expand(**flags)
    assert actual.to_sympy()==expected
    assert str(actual)==str(expected)
    assert sx.expand(value,**flags)==actual
    assert pickle.loads(pickle.dumps(actual))==actual


def test_expand_preserves_distinct_alpha_equal_subs_arguments():
    x,y,t=sx.symbols('x y t')
    f=sx.Function('f')
    left,right=sx.Subs(f(x),x,t),sx.Subs(f(y),y,t)
    assert left==right and left.args!=right.args
    original=sx.Function('pair')(left,right)
    for flags in ({},{'deep':False},{'mul':False},{'power_exp':False}):
        actual=original.expand(**flags)
        assert actual.args[0].args==left.args
        assert actual.args[1].args==right.args


def test_shared_dag_and_huge_power_expansion_are_bounded():
    source='''
import symox as s
x,y=s.symbols('x y')
pair=s.Function('pair')
value=x
for _ in range(80): value=pair(value,value)
expanded=value.expand()
assert expanded==value
assert len(expanded.__reduce__()[1][1])==81
for n in [10**9,2**64,-2**63]:
    original=(x+y)**n
    before=hash(original)
    try: original.expand(max_terms=1000)
    except ValueError as error: assert 'limit' in str(error) or 'budget' in str(error)
    else: raise AssertionError('huge expansion accepted')
    assert hash(original)==before
'''
    run=subprocess.run([sys.executable,'-c',source],capture_output=True,text=True,timeout=15)
    assert run.returncode==0,run.stderr


def test_mul_and_power_work_limits_are_checked_before_expansion():
    x,y,z=sx.symbols('x y z')
    with pytest.raises(ValueError,match='limit|budget'):
        ((x+y)*(x+z)*(y+z)).expand(max_terms=3)
    with pytest.raises(ValueError,match='budget'):
        ((x+y)**10000).expand()
    assert ((x+y)**10000).expand(multinomial=False)==(x+y)**10000
    assert ((x**2147483647+y)**2).expand()==x**4294967294+2*x**2147483647*y+y**2


@pytest.mark.parametrize('limit',[0,-1,1.5,True,None])
def test_invalid_term_limits_fail_explicitly(limit):
    with pytest.raises((TypeError,ValueError)):
        sx.expand(sx.Symbol('x'),max_terms=limit)


@pytest.mark.parametrize('mutable',[False,True])
@pytest.mark.parametrize('deep',[False,True])
def test_explicit_matrix_entries_follow_the_same_options(mutable,deep):
    x,y=sx.symbols('x y')
    cls=sx.Matrix if mutable else sx.ImmutableMatrix
    value=cls([[sx.sin((x+y)**2),(x+y)**2]])
    expected=sp.Matrix(value.to_sympy()).applyfunc(lambda e:e.expand(deep=deep,mul=False))
    assert value.expand(deep=deep,mul=False).to_sympy()==expected
    assert sx.expand(value,deep=deep,mul=False)==value.expand(deep=deep,mul=False)
