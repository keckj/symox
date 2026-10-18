"""Numerical chopping compares native binary magnitudes at declared precision."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('precision',[5,15,40,80])
@pytest.mark.parametrize('source',[
    'Float("1e-100",120)',
    'Float("-1e-100",120)',
    '1+I*Float("1e-100",120)',
    'Float("1e-100",120)+2*I',
    'Float("1e-100",120)*(1+I)',
    'Float("1e-10",120)+I*Float("1e-100",120)',
    'pi+I*Float("1e-100",120)',
])
def test_precision_based_chopping_matches_reference(precision,source):
    expression=eval(source,vars(sx));foreign=eval(source,vars(sp))
    actual=expression.evalf(precision,chop=True)
    expected=foreign.evalf(precision,chop=True)
    assert actual.to_sympy()==expected
    assert sx.N(expression,precision,chop=True)==actual
    assert expression.n(precision,chop=True)==actual


@pytest.mark.parametrize('tolerance',[sx.Rational(1,10**4),sx.Float('1e-10'),sx.Rational(1,10**30),1])
@pytest.mark.parametrize('literal',['1e-4','2e-4','1e-10','2e-10','1e-30','2e-30','1e-50','-1e-50'])
def test_explicit_chop_tolerances_match_reference(tolerance,literal):
    expression=sx.Float(literal,80)+sx.I*sx.Float('1e-60',80)
    actual=expression.evalf(50,chop=tolerance)
    foreign=expression.to_sympy().evalf(50,chop=sx.sympify(tolerance).to_sympy())
    assert actual.to_sympy()==foreign


@pytest.mark.parametrize('literal',['1e-15','1e-16','5e-17','2e-17','1e-17'])
def test_boolean_chopping_uses_binary_threshold_boundaries(literal):
    expression=sx.Float(literal,80)
    assert expression.evalf(15,chop=True).to_sympy()==expression.to_sympy().evalf(15,chop=True)


@pytest.mark.parametrize('source',[
    'x+Float("1e-30")','Float("1e-30")*x','sin(x+Float("1e-30"))',
    'Function("f")(Float("1e-30"))',
])
def test_symbolic_partial_evaluation_preserves_small_coefficients(source):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=eval(source,dict(vars(sx),x=x));foreign=eval(source,dict(vars(sp),x=X))
    assert expression.evalf(chop=True).to_sympy()==foreign.evalf(chop=True)


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
def test_matrix_chopping_and_substitutions(kind):
    x=sx.Symbol('x');X=sp.Symbol('x')
    matrix=getattr(sx,kind)([[x,sx.Float('1e-50',80)],[1+sx.I*sx.Float('1e-50',80),0]])
    reference=getattr(sp,kind)([[X,sp.Float('1e-50',80)],[1+sp.I*sp.Float('1e-50',80),0]])
    original=matrix.copy()
    actual=matrix.evalf(30,subs={x:sx.Float('1e-50',80)},chop=True)
    expected=reference.evalf(30,subs={X:sp.Float('1e-50',80)},chop=True)
    assert actual.to_sympy()==expected and matrix==original
    assert type(actual) is type(matrix)
    assert matrix.n(30,subs={x:sx.Float('1e-50',80)},chop=True)==actual


def test_extreme_decimal_tolerances_do_not_underflow_in_threshold_conversion():
    for exponent in [1000,2000]:
        tolerance=sx.Float('1e-'+str(exponent),100)
        assert sx.Float('1e-'+str(exponent),100).evalf(80,chop=tolerance)==0
        assert sx.Float('2e-'+str(exponent),100).evalf(80,chop=tolerance)!=0
        assert sx.Float('1e-'+str(exponent+1000),100).evalf(80,chop=tolerance)==0


def test_zero_and_false_disable_chopping_and_native_booleans_are_flags():
    expression=sx.Float('1e-50',80)
    assert expression.evalf(chop=False)==expression.evalf(chop=0)==expression.evalf(chop=sx.false)
    assert expression.evalf(chop=sx.true)==0
    assert sx.pi.evalf(40,chop=True)==sx.pi.evalf(40)
    for tolerance in [-1,sx.I,sx.oo,sx.nan,sx.Symbol('epsilon')]:
        with pytest.raises(ValueError):expression.evalf(chop=tolerance)


def test_chop_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
tiny=Float('1e-1000',100)
assert tiny.evalf(chop=True)==0
assert N(1+I*tiny,50,chop=True)==Float(1,50)
assert Matrix([[tiny,1]]).evalf(chop=True)==Matrix([[0,Float(1)]])
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
