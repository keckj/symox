"""Wildcards retain their type, exclusion metadata and assumption facts."""
import copy
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def test_wildcard_type_flags_heads_and_queries():
    a=sx.Wild('a');x=sx.Symbol('x')
    assert isinstance(a,sx.Wild) and isinstance(a,sx.Symbol)
    assert a.is_Wild and not x.is_Wild
    assert a.func is sx.Wild and a.class_name=='Wild'
    assert a.class_key()==sp.Wild('a').class_key()
    assert a.sort_key()==sp.Wild('a').sort_key()
    expression=a*x+sx.sin(a)
    assert expression.atoms(sx.Wild)=={a}
    assert expression.atoms(sp.Wild)=={a}
    assert expression.atoms(a)=={a}
    assert expression.find(sx.Wild)=={a}
    assert expression.has(sx.Wild)
    assert expression.find(sx.Symbol)=={a,x}
    assert expression.replace(sx.Wild,lambda:x)==x*x+sx.sin(x)


@pytest.mark.parametrize('assumptions',[{}, {'positive':True},{'integer':True},{'real':False}])
def test_exclusions_and_assumptions_roundtrip_reference(assumptions):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    a=sx.Wild('a',exclude=[y,x,x],**assumptions)
    expected=sp.Wild('a',exclude=[Y,X,X],**assumptions)
    assert a.exclude==(y,x,x) and a.properties==()
    assert a.to_sympy()==expected
    assert sx.sympify(expected)==a
    assert a==sx.Wild('a',exclude=[y,x,x],**assumptions)
    assert a!=sx.Wild('a',exclude=[x,y],**assumptions)
    assert sx.Wild('a',exclude=[x])!=sx.Wild('a',exclude=[y])
    assert sx.Wild('a',**assumptions)!=a
    assert a.srepr()==sp.srepr(expected)
    assert a.is_positive==expected.is_positive
    for copied in [pickle.loads(pickle.dumps(a)),copy.copy(a),copy.deepcopy(a)]:
        assert copied==a and copied.exclude==(y,x,x)
        assert copied.to_sympy()==expected
        assert copied.is_Wild
    value=sp.sympify('2')
    assert {k.to_sympy():v.to_sympy() for k,v in sx.Integer(2).match(a).items()}==value.match(expected)
    assert (x+1).match(a) is None


def test_exclusion_input_is_copied_and_nested_parser_metadata_is_retained():
    x=sx.Symbol('x');excluded=[x]
    a=sx.Wild('a',exclude=excluded,positive=True)
    excluded.clear()
    assert a.exclude==(x,)
    assert sx.sympify("Wild('a',exclude=[Symbol('x')],positive=True)")==a
    nested=sx.Wild('outer',exclude=[a])
    parsed=sx.sympify("Wild('outer',exclude=[Wild('a',exclude=[Symbol('x')],positive=True)])")
    assert parsed==nested and parsed.exclude==(a,)
    ordinary=sx.sympify("Wild('ordinary',real=True)")
    assert ordinary.is_Wild and ordinary.is_real


def test_foreign_wildcards_can_be_used_as_native_patterns():
    x=sx.Symbol('x');X=sp.Symbol('x');A=sp.Wild('a',exclude=[X])
    assert (3*x).match(A*X)=={sx.sympify(A):sx.Integer(3)}
    assert (sx.sin(x)).match(A) is None
    result=(sx.sin(x)).replace(sp.sin(sp.Wild('argument')),lambda argument:sx.cos(argument))
    assert result==sx.cos(x)


def test_symbols_honors_wild_and_custom_constructors():
    x=sx.Symbol('x')
    a,b=sx.symbols('a b',cls=sx.Wild,exclude=[x])
    assert a.is_Wild and b.is_Wild and a.exclude==(x,)
    assert sx.symbols('w:3',cls=sx.Wild)==tuple(sx.Wild('w%d'%i) for i in range(3))
    calls=[]
    def constructor(name,**facts):calls.append((name,facts));return 'created:'+name
    assert sx.symbols('x y',cls=constructor,real=True)==('created:x','created:y')
    assert calls==[('x',{'real':True}),('y',{'real':True})]
    with pytest.raises(TypeError):sx.symbols('x',cls=object())


def test_predicate_properties_are_enforced_in_native_and_foreign_patterns():
    for a in [sx.Wild('a',properties=[lambda value:value.is_positive]),
              sx.sympify(sp.Wild('a',properties=[lambda value:value.is_positive]))]:
        assert sx.Integer(2).match(a)=={a:sx.Integer(2)}
        assert sx.Integer(-2).match(a) is None
    assert sx.Wild('a',properties=[]).properties==()
    with pytest.raises(AttributeError):_ = sx.Symbol('x').exclude


def test_wildcard_protocols_do_not_import_sympy():
    source='''
import sys,pickle
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');a=Wild('a',exclude=[x],positive=True)
assert a.is_Wild and a.exclude==(x,) and a.is_positive
assert (3*x).match(a*x)=={a:Integer(3)}
assert (a+x).find(Wild)=={a}
assert pickle.loads(pickle.dumps(a))==a
assert sympify("Wild('a',exclude=[Symbol('x')],positive=True)")==a
def integer_value(value):return value.is_integer
a=Wild('constrained',properties=[integer_value])
assert Integer(3).match(a)=={a:Integer(3)}
assert Rational(1,2).match(a) is None
assert pickle.loads(pickle.dumps(a)).properties==(integer_value,)
assert gcd(Tuple(6,9))==3 and lcm(Tuple(6,9))==18
'''
    subprocess.run([sys.executable,'-c',source],check=True)
