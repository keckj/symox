"""Integer rounding must preserve arbitrary magnitudes and exact boundaries."""
from decimal import Decimal
import subprocess
import sys
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source',[
    '10**100 + sqrt(2)', '-10**100 - sqrt(2)',
    'sqrt(10**200 + 1)', '-sqrt(10**200 + 1)',
    '1 - 10**(-1000)', '-1 + 10**(-1000)',
    'exp(10**(-1000))', 'exp(-10**(-1000))',
    'pi', '-pi', 'E', '-E', 'exp(5)', 'sqrt(2) + sqrt(3)',
])
def test_exact_real_integer_parts(source):
    a=sx.sympify(source);b=sp.sympify(source)
    # SymPy 1.14's int conversion is off by one on these radicals and is
    # undecided on the tiny exponentials. Use exact inequalities as the oracle.
    n=10**100
    proved={
        'sqrt(10**200 + 1)':(n,n,n+1),
        '-sqrt(10**200 + 1)':(-n,-n-1,-n),
        'exp(10**(-1000))':(1,1,2),
        'exp(-10**(-1000))':(0,0,1),
    }
    truncation,floor,ceiling=proved.get(source,(None,None,None))
    if truncation is None:
        assert int(a)==int(b)
        assert sx.Integer(a).to_sympy()==sp.Integer(b)
    else:
        assert int(a)==truncation
        assert sx.Integer(a)==truncation
    for name in ['floor','ceiling']:
        result=getattr(sx,name)(a)
        assert result.is_Integer
        expected=(floor if name=='floor' else ceiling) if source in proved else getattr(sp,name)(b)
        assert result.to_sympy()==expected


@pytest.mark.parametrize('value',[1e100,-1e100,1.999,-1.999,2.0,-2.0,0.0])
@pytest.mark.parametrize('digits',[15,50,100])
def test_float_integer_parts_are_exact_for_the_stored_value(value,digits):
    a=sx.Float(value,digits);b=sp.Float(value,digits)
    assert int(a)==int(b)
    for name in ['floor','ceiling']:
        result=getattr(sx,name)(a)
        assert result.is_Integer
        assert result.to_sympy()==getattr(sp,name)(b)


def test_floor_complex_parts_and_integer_translations():
    n=sx.Symbol('n',integer=True);N=sp.Symbol('n',integer=True)
    x=sx.Symbol('x',real=True);X=sp.Symbol('x',real=True)
    for a,b in [(sx.I,sp.I),(sx.Rational(3,2)+sx.Rational(7,3)*sx.I,sp.Rational(3,2)+sp.Rational(7,3)*sp.I),
                (n+x,N+X),(n+sx.sqrt(2),N+sp.sqrt(2)),
                (sx.Float('2.5',70)+sx.I*sx.Float('-3.5',70),sp.Float('2.5',70)+sp.I*sp.Float('-3.5',70))]:
        for name in ['floor','ceiling']:
            assert getattr(sx,name)(a).to_sympy()==getattr(sp,name)(b)
    assert sx.Integer(Decimal('2.9'))==2
    assert sx.Integer(True)==1


def test_exact_function_results_precede_mp_numeric_dispatch():
    assert sx.arg(sx.Float(-2,100))==sx.pi
    assert sx.Max(sx.Float('1.5',100),2).is_Integer
    assert sx.DiracDelta(sx.Float('1e-1000',100))==0
    assert sx.Heaviside(sx.Float('-1e-1000',100))==0


def test_certified_rounding_has_no_sympy_runtime_dependency():
    code='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('reference import attempted')
sys.meta_path.insert(0,Block())
import symox as s
assert int(10**100+s.sqrt(2))==10**100+1
assert s.ceiling(s.sqrt(10**200+1))==10**100+1
assert s.floor(s.exp(-s.Rational(1,10**1000)))==0
assert s.ceiling(s.exp(s.Rational(1,10**1000)))==2
'''
    subprocess.run([sys.executable,'-c',code],check=True)
