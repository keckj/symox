"""ODE parameters must remain distinct from generated solver symbols."""
import subprocess
import sys
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('name',['C1','C2','C3','_y0','_y1','_y2','_r','_u','_v'])
@pytest.mark.parametrize('initial',[False,True])
def test_parameters_are_not_captured_by_constants_or_placeholders(name,initial):
    x=sx.Symbol('x');X=sp.Symbol('x');a=sx.Symbol(name);A=sp.Symbol(name)
    f=sx.Function('f');F=sp.Function('f')
    equation=f(x).diff(x)-a*f(x)
    foreign=F(X).diff(X)-A*F(X)
    options={'ics':{f(0):3}} if initial else {}
    reference_options={'ics':{F(0):3}} if initial else {}
    actual=sx.dsolve(equation,f(x),**options)
    expected=sp.dsolve(foreign,F(X),**reference_options)
    assert actual.to_sympy()==expected
    assert sx.checkodesol(equation,actual,f(x))==(True,0)
    if initial:assert actual.rhs.subs(x,0)==3


@pytest.mark.parametrize('name',['C1','_y0','_y1','_r'])
def test_the_independent_variable_can_use_a_reserved_name(name):
    x=sx.Symbol(name);X=sp.Symbol(name);f=sx.Function('f');F=sp.Function('f')
    equation=f(x).diff(x)-x*f(x)
    actual=sx.dsolve(equation,f(x))
    expected=sp.dsolve(F(X).diff(X)-X*F(X),F(X))
    assert actual.to_sympy()==expected
    assert sx.checkodesol(equation,actual,f(x))==(True,0)


def test_initial_values_can_contain_names_of_integration_constants():
    x=sx.Symbol('x');C1,C2=sx.symbols('C1 C2');f=sx.Function('f')
    equation=f(x).diff(x)-f(x)
    for value in [C1+1,C1+C2]:
        solution=sx.dsolve(equation,f(x),ics={f(0):value})
        assert solution.lhs==f(x)
        assert (solution.rhs-value*sx.exp(x)).expand()==0
        assert solution.rhs.subs(x,0)==value
        assert sx.checkodesol(equation,solution,f(x))==(True,0)


def test_multiple_constants_skip_all_parameter_names():
    x=sx.Symbol('x');C1,C2,C3=sx.symbols('C1 C2 C3');f=sx.Function('f')
    equation=f(x).diff(x,2)-C1*f(x).diff(x)
    solution=sx.dsolve(equation,f(x))
    assert solution==sx.Eq(f(x),C2+C3*sx.exp(C1*x))
    assert sx.checkodesol(equation,solution,f(x))==(True,0)


def test_parameter_assumptions_are_preserved_during_internal_renaming():
    x=sx.Symbol('x');a=sx.Symbol('_y0',positive=True);f=sx.Function('f')
    equation=f(x).diff(x)-a*f(x)
    solution=sx.dsolve(equation,f(x),ics={f(0):2})
    assert solution==sx.Eq(f(x),2*sx.exp(a*x))
    assert a in solution.free_symbols and a.is_positive is True


def test_initial_conditions_for_other_functions_or_excessive_orders_are_rejected():
    x=sx.Symbol('x');f,g=sx.Function('f'),sx.Function('g')
    equation=f(x).diff(x)-f(x)
    # SymPy rejects these with ValueError. A derivative order such as
    # 2**32-1 exhausts memory in SymPy too and is not exercised.
    for condition in [f(),g(0)]:
        with pytest.raises(ValueError,match="Invalid boundary conditions"):
            sx.dsolve(equation,f(x),ics={condition:1})


def test_parameter_isolation_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');f=Function('f')
for name in ['C1','_y0','_r']:
    parameter=Symbol(name)
    equation=f(x).diff(x)-parameter*f(x)
    solution=dsolve(equation,f(x),ics={f(0):3})
    assert solution==Eq(f(x),3*exp(parameter*x))
    assert checkodesol(equation,solution,f(x))==(True,0)
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=15)
