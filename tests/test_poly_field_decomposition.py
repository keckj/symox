"""Functional decomposition stays in the selected exact coefficient domain."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def compare(expression,domain):
    x=sx.Symbol('x');X=sp.Symbol('x')
    polynomial=sx.Poly(expression,x,domain=domain)
    reference=sp.Poly(expression.to_sympy(),X,domain=domain)
    actual=polynomial.decompose();expected=reference.decompose()
    assert len(actual)==len(expected)
    for part,foreign in zip(actual,expected):
        assert sx.cancel(part.as_expr()-sx.sympify(str(foreign.as_expr())))==0
        assert str(part.domain)==str(foreign.domain)
        assert part.gens==polynomial.gens
    rebuilt=actual[0]
    for inner in actual[1:]:rebuilt=rebuilt.compose(inner)
    assert rebuilt==polynomial
    assert sx.decompose(expression,x,domain=domain)==[part.as_expr() for part in actual]
    assert sx.decompose(expression,x,domain=domain,polys=True)==actual
    assert sx.decompose(polynomial)==actual
    return actual


@pytest.mark.parametrize('domain',['QQ(t)','QQ(t,u)','QQ[t]','ZZ[t]'])
@pytest.mark.parametrize('source',[
    '0','1','x','x**3+x+1',
    '(x**2+t*x)**2',
    '(x**2+t*x)**3+2*(x**2+t*x)+1',
    '((x**2+t*x)**2+3*(x**2+t*x))**2+1',
])
def test_symbolic_coefficient_domains(source,domain):
    x,t,u=sx.symbols('x t u')
    compare(sx.sympify(eval(source,{'x':x,'t':t,'u':u})),domain)


@pytest.mark.parametrize('source',[
    '(x**2+x/t)**2+t',
    '(x**2+(t+u)/(t-u)*x)**3+u',
    '(t*x**2+x)**2',
    '((x**2+x/(t+1))**2+x**2+x/(t+1))**2',
])
def test_fractional_coefficients_and_nested_compositions(source):
    x,t,u=sx.symbols('x t u')
    compare(eval(source,{'x':x,'t':t,'u':u}),'QQ(t,u)')


@pytest.mark.parametrize('domain,root',[
    ('QQ<sqrt(2)>',sx.sqrt(2)),
    ('QQ<I>',sx.I),
    ('QQ<sqrt(2) + sqrt(3)>',sx.sqrt(2)+sx.sqrt(3)),
])
@pytest.mark.parametrize('power',[2,3,4])
def test_algebraic_number_fields(domain,root,power):
    x=sx.Symbol('x');inner=x*x+root*x
    compare(inner**power+2*inner+root,domain)


@pytest.mark.parametrize('domain',['ZZ','QQ','ZZ[t]','QQ[t]','QQ(t)'])
def test_decomposition_never_leaves_a_coefficient_ring(domain):
    x,t=sx.symbols('x t')
    compare(4*x**4+4*x**3+x*x,domain)
    if '[' in domain or '(' in domain:compare((t*x*x+x)**2,domain)


@pytest.mark.parametrize('prime',[3,5,7,11])
@pytest.mark.parametrize('source',['x**4+2*x**2+1','(x**2+2*x)**2+3','x**3+x+1'])
def test_finite_field_arithmetic_is_not_rational_arithmetic(prime,source):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=eval(source,{'x':x});polynomial=sx.Poly(expression,x,modulus=prime)
    reference=sp.Poly(expression.to_sympy(),X,modulus=prime)
    actual=polynomial.decompose();expected=reference.decompose()
    assert [part.as_expr().to_sympy() for part in actual]==[part.as_expr() for part in expected]
    assert all(part.modulus==prime for part in actual)
    reconstructed=actual[0]
    for inner in actual[1:]:reconstructed=reconstructed.compose(inner)
    assert reconstructed==polynomial


def test_characteristic_obstruction_raises_instead_of_inverting_zero():
    x=sx.Symbol('x')
    with pytest.raises(sx.NotInvertible):sx.Poly(x**4+1,x,modulus=2).decompose()


def test_multivariate_and_sparse_huge_degrees_fail_before_dense_allocation():
    x,y,t=sx.symbols('x y t')
    for domain in ['QQ(t)','QQ<sqrt(2)>']:
        with pytest.raises(NotImplementedError):sx.Poly(x*x+y,x,y,domain=domain).decompose()
        with pytest.raises(NotImplementedError):sx.Poly(x**1000000000+1,x,domain=domain).decompose()
    with pytest.raises(NotImplementedError):sx.Poly(x**1000000000+1,x).decompose()


def test_field_decomposition_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,t=symbols('x t')
for domain,inner in [('QQ(t)',x*x+x/t),('QQ<sqrt(2)>',x*x+sqrt(2)*x)]:
    p=Poly(inner**3+inner,x,domain=domain)
    parts=p.decompose()
    assert len(parts)==2 and parts[0].compose(parts[1])==p
    assert decompose(p)==parts
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
