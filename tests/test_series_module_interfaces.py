"""Series module imports exercise their native algorithms, not just names."""
import importlib
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.series import approximants,gruntz,limit_seq
from sympy.series.approximants import approximants as reference_approximants


@pytest.mark.parametrize('source',[
    '[0]*5','[1]*12','[fibonacci(k) for k in range(16)]',
    '[lucas(k) for k in range(16)]','[lucas(k)**2 for k in range(16)]',
    '[lucas(k)+fibonacci(k+2) for k in range(16)]',
    '[Rational(1,factorial(k)) for k in range(10)]',
    '[1,0,1,0,1,0,1,0]','[0,0,1,2,3,4,5,6]',
    '[(a+1)**k for k in range(10)]',
])
def test_pade_approximants_match_reference_and_reproduce_coefficients(source):
    a,x=sx.symbols('a x');A,X=sp.symbols('a x')
    native=eval(source,vars(sx)|{'a':a})
    reference=eval(source,vars(sp)|{'a':A})
    actual=list(approximants(native,x,simplify=True))
    expected=list(reference_approximants(reference,X,simplify=True))
    expected=[value for value in expected if value is not sp.nan]
    if source=='[(a+1)**k for k in range(10)]':
        # SymPy 1.14 can emit a spurious final zero after exact termination.
        assert len(actual)==2
        assert all(sp.simplify(value.subs(X,0)-reference[0])!=0
                   for value in expected[len(actual):])
        expected=expected[:len(actual)]
    assert len(actual)==len(expected)
    for left,right in zip(actual,expected):
        assert sp.cancel(left.to_sympy()-right)==0
    if actual:
        expanded=actual[-1].to_sympy().series(X,0,len(native)).removeO().expand()
        assert all(sp.simplify(expanded.coeff(X,k)-reference[k])==0 for k in range(len(native)))


@pytest.mark.parametrize('source,point,direction',[
    ('sin(x)/x',0,'+'),('sin(x)/x',0,'-'),('1/x',0,'-'),
    ('exp(x)/x**5','oo','+'),('log(x)/x','oo','+'),
    ('(1+1/x)**x','oo','+'),('(x*x+3)/(2*x*x-1)','oo','+'),
    ('exp(-exp(x))*exp(x)','oo','+'),
])
def test_gruntz_uses_exact_native_mrv_limits(source,point,direction):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=eval(source,vars(sx)|{'x':x})
    reference=eval(source,vars(sp)|{'x':X})
    p=sx.sympify(point)
    assert gruntz(expression,x,p,direction).to_sympy()==sp.limit(reference,X,p.to_sympy(),dir=direction)


@pytest.mark.parametrize('source',[
    '(5*n**3+3*n**2+4)/(3*n**3+4*n-5)',
    '(1+1/n)**n','log(n)/n','factorial(n)/factorial(n+1)',
    'gamma(n+2)/gamma(n+1)/n','2**n/factorial(n)','factorial(n)/2**n',
    '(-1)**n','(-1)**n/n','(-1)**n*(1+1/n)','cos(pi*n)',
    'fibonacci(n)/fibonacci(n+1)',
    'binomial(2*n,n)/Sum(binomial(2*k,k),(k,1,n))',
    'Sum(k*k,(k,1,n))/n**3','7',
])
def test_sequence_limits_match_reference_with_exact_discrete_certificates(source):
    n,k=sx.symbols('n k',integer=True,positive=True)
    N,K=sp.symbols('n k',integer=True,positive=True)
    actual=limit_seq(eval(source,vars(sx)|{'n':n,'k':k}),n)
    expected=sp.limit_seq(sp.sympify(eval(source,vars(sp)|{'n':N,'k':K})),N)
    assert actual is not None
    converted=actual.to_sympy()
    if isinstance(expected,sp.AccumBounds):assert converted==expected
    elif converted in (sp.oo,-sp.oo):assert converted==expected
    else:assert sp.simplify(converted-expected)==0


def test_sequence_limit_validation_and_unsupported_unknown_slope():
    n,a=sx.symbols('n a')
    assert limit_seq(3)==3
    assert limit_seq((n+1)/n)==1
    assert limit_seq(sx.Function('F')(n),n) is None
    with pytest.raises(ValueError):limit_seq(n+a)
    # As in the reference: an index that is not a free symbol leaves the
    # expression unchanged, and a large trial count is only an upper bound.
    assert limit_seq(n,n+1)==n
    assert limit_seq(n,n,trials=100000)==sx.oo
    with pytest.raises(ValueError):list(approximants(range(257)))


def test_lucas_fibonacci_ratio_uses_the_exact_binet_limit():
    # SymPy 1.14's sequence-limit rewrite handles Fibonacci but overlooks Lucas.
    # Their exact Binet formulas give 2*GoldenRatio - 1 = sqrt(5).
    n=sx.Symbol('n',integer=True,positive=True)
    assert limit_seq(sx.lucas(n)/sx.fibonacci(n),n)==2*sx.S.GoldenRatio-1


def test_module_paths_use_the_existing_native_limit_order_and_residue_implementations():
    series=importlib.import_module('symox.series')
    x=sx.Symbol('x')
    for name in ['Limit','Order','O','limit','residue','difference_delta']:
        assert getattr(series,name) is getattr(sx,name)
    assert series.limit(sx.sin(x)/x,x,0)==1
    assert series.Limit(sx.sin(x)/x,x,0).doit()==1
    assert series.O(x*x,x)==sx.O(x*x,x)
    assert series.residue(sx.exp(x)/x**3,x,0)==sx.Rational(1,2)
    assert sx.expand(series.difference_delta(x*x,x))==2*x+1
    for module,name in [('limits','limit'),('order','Order'),('residues','residue'),
                        ('approximants','approximants'),('gruntz','gruntz'),('limitseq','limit_seq')]:
        assert getattr(importlib.import_module('symox.series.'+module),name) is getattr(series,name)


def test_accumulation_interval_arithmetic_matches_reference():
    native=sx.AccumBounds(-2,3);reference=sp.AccumBounds(-2,3)
    for actual,expected in [(native+2,reference+2),(native*2,reference*2),
                            (native**2,reference**2),(native**3,reference**3),
                            (native**-1,reference**-1),(-native,-reference)]:
        assert actual.to_sympy()==expected
    assert sx.AccumBounds(2,2)==2
    assert sx.AccumBounds(-sx.oo,sx.oo)**2==sx.AccumBounds(0,sx.oo)
    assert sx.AccumBounds(0,sx.oo)**-1==sx.AccumBounds(0,sx.oo)
    with pytest.raises(ValueError):sx.AccumBounds(2,1)
    with pytest.raises(ValueError):sx.AccumBounds(sx.I,2)


def test_series_completion_has_no_sympy_runtime_dependency():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.series import approximants,gruntz,limit_seq,Order,Limit
n=Symbol('n',integer=True,positive=True);x=Symbol('x')
assert cancel(list(approximants([1]*8,x))[-1]-1/(1-x))==0
assert limit_seq((-1)**n/n,n)==0
assert gruntz(sin(x)/x,x,0)==1
assert Limit(sin(x)/x,x,0).doit()==1
'''
    subprocess.run([sys.executable,'-c',source],check=True)
