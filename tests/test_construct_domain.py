"""Native inference preserves coefficient values, precision and containers."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.polys.constructor import construct_domain


@pytest.mark.parametrize('case', range(12))
@pytest.mark.parametrize('field', [False,True])
@pytest.mark.parametrize('shape', ['list','tuple','dict','generator'])
def test_exact_and_numeric_inference(case,field,shape):
    x,y = sx.symbols('x y')
    X,Y = sp.symbols('x y')
    examples = [[],[1,2],[1,sx.Rational(1,3)],[x,y],[x/2,y],
        [1/x,y/(x+1)],[sx.I],[sx.I/2],[sx.Float('0.123456789123456789',60)],
        [sx.Float('0.123456789123456789',60)+sx.I],[sx.sqrt(2)],[sx.pi]]
    reference = [[],[1,2],[1,sp.Rational(1,3)],[X,Y],[X/2,Y],
        [1/X,Y/(X+1)],[sp.I],[sp.I/2],[sp.Float('0.123456789123456789',60)],
        [sp.Float('0.123456789123456789',60)+sp.I],[sp.sqrt(2)],[sp.pi]]
    obj, expected_obj = examples[case],reference[case]
    if shape == 'dict':
        obj = {(i,0):value for i,value in enumerate(obj)}
        expected_obj = {(i,0):value for i,value in enumerate(expected_obj)}
    elif shape == 'tuple':
        obj,expected_obj = tuple(obj),tuple(expected_obj)
    elif shape == 'generator':
        obj,expected_obj = iter(obj),iter(expected_obj)
    domain,elements = construct_domain(obj,field=field)
    expected_domain,expected_elements = sp.construct_domain(expected_obj,field=field)
    if case == 11:
        # Native composite domains require Symbol generators. Constants such
        # as pi stay formal rather than becoming an independent ring variable.
        assert domain == sx.EX
    else:
        assert str(domain) == str(expected_domain)
    if case in (8,9):
        assert domain.precision == expected_domain.precision
    if shape == 'dict':
        assert elements.keys() == expected_elements.keys()
        elements,expected_elements = list(elements.values()),list(expected_elements.values())
    assert len(elements) == len(expected_elements)
    for element,expected in zip(elements,expected_elements):
        actual = domain.to_sympy(element).to_sympy()
        expected = expected_domain.to_sympy(expected)
        assert sp.simplify(actual-expected) == 0
        assert element.parent() == domain


def test_scalar_input_and_expression_fallback():
    x = sx.Symbol('x')
    for obj in [3,sx.Rational(1,3),sx.I,sx.sqrt(2),x,sx.sin(x)]:
        domain,element = construct_domain(obj)
        assert element.parent() == domain
        assert sx.simplify(domain.to_sympy(element)-obj) == 0
    assert construct_domain(sx.sin(x))[0] == sx.EX
    assert construct_domain(x,composite=False)[0] == sx.EX
    assert construct_domain([x,1],composite=True)[0] == sx.ZZ.poly_ring(x)
    assert construct_domain([]) == (sx.ZZ,[])
    assert construct_domain({}) == (sx.ZZ,{})


@pytest.mark.parametrize('values', [[2],[sx.sqrt(2)],[sx.sqrt(2),sx.sqrt(3)],[sx.I]])
def test_algebraic_extension_inference(values):
    domain,elements = construct_domain(values,extension=True)
    expected,_ = sp.construct_domain([sx.sympify(v).to_sympy() for v in values],extension=True)
    assert str(domain) == str(expected)
    for value,element in zip(values,elements):
        assert sx.simplify(domain.to_sympy(element)-value) == 0


def test_invalid_options_and_no_silent_generator_independence():
    from symox.polys.polyerrors import OptionError
    x = sx.Symbol('x')
    for options in [{'unknown':True},{'field':'yes'},{'composite':3},{'extension':False}]:
        with pytest.raises(OptionError):
            construct_domain([1],**options)
    with pytest.raises(NotImplementedError):
        construct_domain([sx.sin(x)],composite=True)
    with pytest.raises(NotImplementedError):
        construct_domain([1],extension=sx.sqrt(2))


def test_no_sympy_import():
    source = '''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols, QQ, Float
from symox.polys.constructor import construct_domain
x=symbols('x')
K,v=construct_domain([x/2,1])
assert K==QQ.poly_ring(x) and K.to_sympy(v[0])==x/2
R,v=construct_domain([Float('0.1',60)])
assert R.precision>190 and v[0].parent()==R
'''
    subprocess.run([sys.executable,'-c',source],check=True)
