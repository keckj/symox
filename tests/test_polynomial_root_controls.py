"""Radical controls must keep complete solve results and strict roots errors."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('source,options',[
    ('x**3-x+1',{'cubics':False}),
    ('x**4-x+1',{'quartics':False}),
    ('x**4-x+1',{'cubics':False}),
    ('x**6-x**2+1',{'cubics':False}),
    ('(x-2)*(x**3-x+1)**2',{'cubics':False}),
    ('(x**3-x+1)/(x-1)',{'cubics':False}),
])
@pytest.mark.parametrize('output',[{}, {'dict':True},{'set':True}])
def test_disabled_formulas_keep_all_exact_indexed_roots(source,options,output):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    actual=sx.solve(sx.sympify(source),x,**options,**output)
    expected=sp.solve(sp.sympify(source),rx,**options,**output)
    if output.get('dict'):
        actual=[{k.to_sympy():v.to_sympy() for k,v in row.items()} for row in actual]
        assert actual==expected
    elif output.get('set'):
        actual=([v.to_sympy() for v in actual[0]],{tuple(v.to_sympy() for v in row) for row in actual[1]})
        assert actual==expected
    else:assert {v.to_sympy() for v in actual}==set(expected)


@pytest.mark.parametrize('assumptions',[{'real':True},{'positive':True},{'negative':True}])
def test_indexed_roots_respect_unknown_assumptions(assumptions):
    x,rx=sx.Symbol('x',**assumptions),sp.Symbol('x',**assumptions)
    actual=sx.solve(x**3-x+1,x,cubics=False)
    expected=sp.solve(rx**3-rx+1,rx,cubics=False)
    assert {v.to_sympy() for v in actual}==set(expected)


@pytest.mark.parametrize('source,options',[
    ('x**3-x+1',{'cubics':False}),
    ('x**4-x+1',{'quartics':False}),
    ('(x-2)*(x**3-x+1)',{'cubics':False}),
    ('x**4+1',{'quartics':False}),
    ('x**4-5*x**2+4',{'cubics':False,'quartics':False}),
])
def test_roots_respects_formula_flags_and_retains_special_factorizations(source,options):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    actual=sx.roots(sx.sympify(source),x,**options)
    expected=sp.roots(sp.sympify(source),rx,**options)
    assert {root.to_sympy():m for root,m in actual.items()}==expected


def test_quartic_formula_can_be_enabled_independently_of_cubics():
    x=sx.Symbol('x');f=x**4-x+1
    roots=sx.solve(f,x,cubics=False,quartics=True)
    assert len(roots)==4 and all(not isinstance(root,sx.CRootOf) for root in roots)
    assert sx.roots(f,x,cubics=False,quartics=True)==sx.roots(f,x)


@pytest.mark.parametrize('source,options',[
    ('x**5-x-1',{}),('(x-1)*(x**5-x-1)',{}),
    ('x**3-x+1',{'cubics':False}),
    ('x**2+1',{'filter':'R'}),
    ('x**2-1',{'predicate':lambda root:False}),
])
def test_strict_roots_raises_instead_of_silently_returning_an_incomplete_dictionary(source,options):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    from symox.polys.polyerrors import UnsolvableFactorError
    from sympy.polys.polyerrors import UnsolvableFactorError as ReferenceError
    with pytest.raises(ReferenceError):sp.roots(sp.sympify(source),rx,strict=True,**options)
    with pytest.raises(UnsolvableFactorError):sx.roots(sx.sympify(source),x,strict=True,**options)


@pytest.mark.parametrize('multiple',[False,True])
def test_predicate_filtering_preserves_reference_zero_root_and_multiplicities(multiple):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    for predicate in [lambda root:False,lambda root:root>0]:
        actual=sx.roots(x**3*(x-1)**2,x,predicate=predicate,multiple=multiple)
        expected=sp.roots(rx**3*(rx-1)**2,rx,predicate=predicate,multiple=multiple)
        converted=([v.to_sympy() for v in actual] if multiple else
                   {v.to_sympy():m for v,m in actual.items()})
        assert converted==expected
    actual=sx.roots(x*(x*x+1),x,filter='I')
    assert {v.to_sympy():m for v,m in actual.items()}==sp.roots(rx*(rx*rx+1),rx,filter='I')


@pytest.mark.parametrize('invalid',['','X','real'])
def test_invalid_filter_is_rejected(invalid):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    with pytest.raises(ValueError):sp.roots(rx*rx-1,rx,filter=invalid)
    with pytest.raises(ValueError):sx.roots(x*x-1,x,filter=invalid)


def test_coefficient_list_rejects_redundant_generator_like_reference():
    with pytest.raises(ValueError,match='redundant generators'):sp.roots([1,0,-1],sp.Symbol('t'))
    with pytest.raises(ValueError,match='redundant generators'):sx.roots([1,0,-1],sx.Symbol('t'))


@pytest.mark.parametrize('generator',['y+1','sin(y)','2'])
def test_polynomial_input_keeps_formal_generator_exponent_data(generator):
    p=sx.Poly([1,0,-1],sx.sympify(generator))
    reference=sp.Poly([1,0,-1],sp.sympify(generator))
    actual=sx.roots(p,strict=True)
    assert {r.to_sympy():m for r,m in actual.items()}==sp.roots(reference,strict=True)


def test_root_controls_leave_scope_before_python_predicates():
    x=sx.Symbol('x');seen=[]
    def predicate(root):
        seen.append(root)
        assert len(sx.roots(x**3-x+1,x))==3
        return True
    assert sx.roots(x-1,x,cubics=False,quartics=False,predicate=predicate)=={sx.Integer(1):1}
    assert seen==[sx.Integer(1)]
    assert len(sx.solve(x**3-x+1,x))==3


def test_root_controls_and_strict_errors_never_import_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
from symox.polys.polyerrors import UnsolvableFactorError
x=Symbol('x')
solutions=solve(x**3-x+1,x,cubics=False)
assert len(solutions)==3 and all(isinstance(root,CRootOf) for root in solutions)
try:roots(x**5-x-1,x,strict=True)
except UnsolvableFactorError:pass
else:raise AssertionError('missing strict error')
assert roots(x*x-1,x,predicate=lambda root:root>0)=={Integer(1):1}
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=20)
