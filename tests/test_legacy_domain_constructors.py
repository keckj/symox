"""Legacy backend spellings construct native rings and exact coefficients."""
from decimal import Decimal
from fractions import Fraction
import importlib
import pickle
import subprocess
import sys

import pytest
import symox as sx
from symox.polys import domains
from symox.polys.domains import PythonRational
from sympy.polys.domains import PythonRational as ReferenceRational


@pytest.mark.parametrize('name',['PythonIntegerRing','GMPYIntegerRing','ZZ_python','ZZ_gmpy'])
def test_integer_backend_names_use_native_exact_arithmetic(name):
    constructor=getattr(domains,name);ring=constructor()
    assert ring==sx.ZZ and isinstance(ring,constructor)
    assert ring(3)*ring(4)==ring(12)
    assert ring.div(ring(17),ring(5))==(ring(3),ring(2))
    assert ring.gcd(ring(30),ring(18))==ring(6)
    assert pickle.loads(pickle.dumps(ring))==ring


@pytest.mark.parametrize('name',['PythonRationalField','GMPYRationalField','QQ_python','QQ_gmpy'])
def test_rational_backend_names_convert_the_python_rational_protocol(name):
    constructor=getattr(domains,name);field=constructor()
    assert field==sx.QQ and isinstance(field,constructor)
    value=field(PythonRational(2,3))
    assert field.to_sympy(value)==sx.Rational(2,3)
    assert value/field(2)==field(sx.Rational(1,3))


@pytest.mark.parametrize('name',['PythonFiniteField','GMPYFiniteField','FF_python','FF_gmpy'])
@pytest.mark.parametrize('symmetric',[False,True])
def test_finite_field_backend_names_validate_the_native_display_option(name,symmetric):
    if not symmetric:
        with pytest.raises(NotImplementedError):getattr(domains,name)(7,symmetric)
        return
    field=getattr(domains,name)(7,symmetric)
    assert field==sx.GF(7,symmetric=symmetric)
    assert field(3)*field(5)==field(1)
    assert field.exquo(field(2),field(3))==field(3)
    assert field.to_sympy(field(6))==(-1 if symmetric else 6)
    with pytest.raises((ValueError,TypeError)):getattr(domains,name)(0)


@pytest.mark.parametrize('value',[0,1,-2,'2/3','-1.125',0.125,Decimal('1.25'),Fraction(3,7)])
def test_python_rational_normalization_and_coercion_match_reference(value):
    actual=PythonRational(value);expected=ReferenceRational(value)
    assert (actual.numerator,actual.denominator)==(expected.numerator,expected.denominator)
    assert hash(actual)==hash(expected)
    assert int(actual)==int(expected) and float(actual)==float(expected)
    assert bool(actual)==bool(expected)
    assert sx.sympify(actual)==sx.Rational(actual.numerator,actual.denominator)
    assert pickle.loads(pickle.dumps(actual))==actual


@pytest.mark.parametrize('p,q',[(2,3),(-6,8),(0,5),(7,-3)])
@pytest.mark.parametrize('operation',['add','sub','mul','truediv','pow'])
def test_python_rational_arithmetic_is_exact_and_closed(p,q,operation):
    import operator
    actual=PythonRational(p,q);expected=ReferenceRational(p,q)
    rhs=2 if operation=='pow' else PythonRational(4,5)
    reference_rhs=2 if operation=='pow' else ReferenceRational(4,5)
    result=getattr(operator,operation)(actual,rhs)
    reference=getattr(operator,operation)(expected,reference_rhs)
    assert isinstance(result,PythonRational)
    assert (result.numerator,result.denominator)==(reference.numerator,reference.denominator)
    assert -actual==PythonRational(-p,q) and abs(actual)>=0


def test_domain_module_paths_and_rational_rejections():
    for module,name in [('pythonrational','PythonRational'),('pythonintegerring','PythonIntegerRing'),
                        ('pythonrationalfield','PythonRationalField'),('pythonfinitefield','PythonFiniteField'),
                        ('gmpyintegerring','GMPYIntegerRing'),('gmpyrationalfield','GMPYRationalField'),
                        ('gmpyfinitefield','GMPYFiniteField')]:
        assert getattr(importlib.import_module('symox.polys.domains.'+module),name) is getattr(domains,name)
    q=PythonRational(2,3)
    assert 1-q==PythonRational(1,3) and 2/q==3 and q**-1==PythonRational(3,2)
    with pytest.raises(ZeroDivisionError):PythonRational(1,0)
    with pytest.raises(ZeroDivisionError):q/0
    with pytest.raises(TypeError):PythonRational(Fraction(1,2),3)
    with pytest.raises(TypeError):q//2
    with pytest.raises(TypeError):q%2


def test_domain_aliases_do_not_import_sympy_or_gmpy():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name.split('.')[0] in ('sympy','gmpy','gmpy2'):
            raise AssertionError('Foreign arithmetic import attempted')
sys.meta_path.insert(0,Block())
from symox.polys.domains import *
assert QQ_python()(PythonRational(2,3))*QQ_gmpy()(3)==QQ(2)
assert GMPYFiniteField(7)(3)*PythonFiniteField(7)(5)==GF(7)(1)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
