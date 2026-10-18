"""Exact comparisons must precede rounded arithmetic in extrema and relations."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def pairs():
    yield sx.S.Zero,sx.Rational(1,10**1000)
    yield sx.Integer(10**100),sx.Integer(10**100+1)
    yield sx.Rational(10**150+1,10**150),sx.S.One
    yield sx.Float('.1'),sx.Rational(1,10)
    yield sx.Float(2**53),sx.Integer(2**53+1)
    yield sx.Float(1),sx.Rational(10**100+1,10**100)
    yield sx.Float('1e-1000',100),sx.Rational(1,10**1000)
    yield sx.Float('1.234567890123456789',80),sx.Rational(1234567890123456789,10**18)
    yield sx.Float('1e500',80),sx.Integer(10**500)
    yield sx.Float('.5'),sx.Rational(1,2)
    yield sx.Float(0),sx.S.Zero
    yield sx.oo,sx.Integer(10**1000)
    yield -sx.oo,sx.Rational(1,10**1000)
    yield -sx.oo,sx.oo


@pytest.mark.parametrize('first,second',list(pairs()))
@pytest.mark.parametrize('reverse',[False,True])
def test_stored_numeric_values_and_extrema_match_exact_reference(first,second,reverse):
    if reverse:first,second=second,first
    a,b=first.to_sympy(),second.to_sympy()
    exact_a=sp.Rational(a) if a.is_Float else a
    exact_b=sp.Rational(b) if b.is_Float else b
    for native,oracle in [(sx.Eq,sp.Eq),(sx.Ne,sp.Ne),(sx.Lt,sp.Lt),
                          (sx.Le,sp.Le),(sx.Gt,sp.Gt),(sx.Ge,sp.Ge)]:
        actual=native(first,second).to_sympy();expected=oracle(exact_a,exact_b)
        assert actual==expected
    for native,oracle in [(sx.Min,sp.Min),(sx.Max,sp.Max)]:
        actual=native(first,second).to_sympy();expected=oracle(a,b)
        assert actual==expected


def test_numeric_atom_extrema_preserve_exact_selected_values():
    tiny=sx.Rational(1,10**1000)
    assert sx.Max(0,tiny)==tiny
    assert sx.Min(0,-tiny)==-tiny
    assert sx.Max(10**100,10**100+1)==10**100+1
    assert sx.Min(10**100+1,10**100)==10**100
    values=[sx.Rational(k,10**1000) for k in range(-20,21)]
    assert sx.Max(*values)==sx.Rational(20,10**1000)
    assert sx.Min(*values)==-sx.Rational(20,10**1000)


def test_infinity_relations_distinguish_finite_and_extended_real_symbols():
    finite=sx.Symbol('finite',real=True);extended=sx.Symbol('extended',extended_real=True)
    assert sx.Gt(sx.oo,finite)==sx.true
    assert sx.Lt(-sx.oo,finite)==sx.true
    assert sx.Eq(sx.oo,finite)==sx.false
    assert sx.Ge(sx.oo,extended)==sx.true
    assert sx.Le(-sx.oo,extended)==sx.true
    assert sx.Gt(sx.oo,extended) not in (sx.true,sx.false)
    assert sx.Eq(sx.nan,sx.nan)==sx.false
    assert sx.Ne(sx.nan,sx.nan)==sx.true


def test_nearby_constant_extrema_keep_unresolved_arguments():
    lower=sx.sqrt(2);upper=sx.Rational(14142135623730951,10**16)
    for name in ['Min','Max']:
        result=getattr(sx,name)(lower,upper)
        assert result.has(lower) and result.has(upper)
        assert result.to_sympy()==getattr(sp,name)(sp.sqrt(2),upper.to_sympy())
    assert sx.Min(3,sx.pi)==3
    assert sx.Max(3,sx.pi)==sx.pi


def test_symbolic_extrema_still_use_signs_and_nested_flattening():
    x,y=sx.symbols('x y');p=sx.Symbol('p',positive=True)
    assert sx.Max(x,x+1)==x+1
    assert sx.Min(x,x+p)==x
    assert sx.Max(x,sx.Max(x,y))==sx.Max(x,y)
    assert sx.Min()==sx.oo and sx.Max()==-sx.oo
    assert sx.Min(x,-sx.oo)==-sx.oo
    assert sx.Max(x,sx.oo)==sx.oo


def test_nonreal_extrema_are_rejected_and_conversion_hooks_run_once():
    for function in [sx.Min,sx.Max]:
        for value in [sx.I,sx.nan,sx.zoo,sx.sqrt(-2),sx.Matrix([[1]])]:
            with pytest.raises((ValueError,TypeError)):function(value,1)
    calls=[]
    class Converted:
        def _symox_(self):calls.append(1);return sx.Integer(3)
    assert sx.Max(Converted(),2)==3
    assert calls==[1]


def test_exact_comparisons_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
tiny=Rational(1,10**1000)
assert Max(0,tiny)==tiny
assert Min(10**100+1,10**100)==10**100
assert Gt(Float('.1'),Rational(1,10))==true
assert Lt(Float(2**53),Integer(2**53+1))==true
'''
    subprocess.run([sys.executable,'-c',source],check=True)
