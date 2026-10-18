"""Exact equations, closure identities and numeric continuation, with an oracle."""
import math
import subprocess
import sys

import pytest
import symox as s
from symox.holonomic import (DifferentialOperator, DifferentialOperators,
    HolonomicFunction, HolonomicSequence, RecurrenceOperator,
    RecurrenceOperators, expr_to_holonomic, from_hyper, from_meijerg)

x = s.Symbol('x')


def assert_annihilates(h, expression):
    import sympy as sp
    variable = x.to_sympy()
    reference = expression.to_sympy()
    residual = sum(c.to_sympy()*sp.diff(reference,variable,i)
                   for i,c in enumerate(h.annihilator.coefficients))
    assert sp.simplify(sp.expand_func(residual)) == 0


@pytest.mark.parametrize('power',[1,2,3,5])
def test_differential_ore_relation(power):
    R,D = DifferentialOperators(s.QQ.old_poly_ring(x),'Dx')
    assert D*x == x*D+1
    product = D**power*x**3
    f = s.exp(x)*s.sin(x)
    applied = sum(c*f.diff(x,i) for i,c in enumerate(product.coefficients))
    assert s.simplify(applied-(x**3*f).diff(x,power)) == 0
    assert R.derivative_operator == D
    assert (D+1)**0 == R(1)


def test_shift_ore_relation_and_fibonacci():
    n = s.Symbol('n',integer=True)
    R,S = RecurrenceOperators(s.ZZ.old_poly_ring(n),'Sn')
    assert S*n == (n+1)*S
    assert S**3*n**2 == (n+3)**2*S**3
    assert n*S*n+1-S**2*n == RecurrenceOperator([1,n**2+n,-n-2],R)
    fib = HolonomicSequence(S**2-S-1,[0,1])
    assert [fib.eval(i) for i in range(12)] == [0,1,1,2,3,5,8,13,21,34,55,89]
    factorial = HolonomicSequence(S-(n+1),[1])
    assert factorial.eval(10) == math.factorial(10)
    with pytest.raises(ValueError):
        HolonomicSequence(n*S,[1]).eval(1)


@pytest.mark.parametrize('expression',[
    s.exp(x),s.exp(x**2),s.sin(x),s.cos(2*x+1),s.sinh(x),s.cosh(x),
    s.log(1+x),s.sqrt(1+x),1/(1+x**2),s.erf(x),s.erfc(x),
    s.exp(x)*s.sin(x),s.sin(x)+s.exp(x),s.sin(x)**2,
    s.besselj(1,x),s.besseli(2,x),s.airyai(x),s.airybi(x),
    s.hyper([1,2],[3],x),s.meijerg([],[],[0],[],x),
])
def test_conversions_satisfy_reference_equation(expression):
    h = expr_to_holonomic(expression,x,initcond=False)
    assert_annihilates(h,expression)
    assert h.to_expr() == expression
    assert h.y0 is None


@pytest.mark.parametrize('expression',[s.exp(x),s.sin(x),s.cos(x),s.log(1+x),
    s.sqrt(1+x),s.exp(x)*s.sin(x),s.hyper([1,2],[3],x)])
def test_taylor_coefficients_and_recurrence(expression):
    import sympy as sp
    h = expr_to_holonomic(expression,x)
    values = h.series(10,coefficient=True)
    reference = expression.to_sympy().series(x.to_sympy(),0,10).removeO()
    assert all(sp.simplify(value.to_sympy()-reference.coeff(x.to_sympy(),i)) == 0 for i,value in enumerate(values))
    sequence,start = h.to_sequence()[0]
    assert all(sequence.eval(i) == value for i,value in enumerate(values))
    assert h.series(6).getO() == s.Order(x**6,x)


def test_manual_closure_and_initial_conditions():
    R,D = DifferentialOperators(s.QQ.old_poly_ring(x),'Dx')
    sine = HolonomicFunction(D**2+1,x,0,[0,1])
    exponential = HolonomicFunction(D-1,x,0,[1])
    product = sine*exponential
    total = sine+exponential
    for h,f in [(product,s.sin(x)*s.exp(x)),(total,s.sin(x)+s.exp(x)),
                (sine.diff(),s.cos(x)),(sine.composition(x**2),s.sin(x**2))]:
        assert_annihilates(h,f)
        assert s.expand(h.series(9,order=False)-s.series(f,x,0,9).removeO()) == 0
    assert sine.to_expr() == s.sin(x)
    assert sine.integrate(x).y0 == [0,0,1]
    assert s.expand(sine.integrate(x).series(8,order=False)-(1-s.cos(x)).series(x,0,8).removeO()) == 0
    assert (s.Integer(2)*sine).y0 == [0,2]
    assert sine.diff(x,0) is sine
    import sympy as sp
    assert sp.simplify(sp.hyperexpand(sine.to_hyper().to_sympy())-s.sin(x).to_sympy()) == 0


def test_composition_shift_and_new_variable():
    h = expr_to_holonomic(s.exp(x),x)
    composed = h.composition(x**2+1)
    assert_annihilates(composed,s.exp(x**2+1))
    shifted = h.shift_x(2)
    assert shifted.x0 == -2 and shifted.y0 == [1]
    t = s.Symbol('t')
    renamed = shifted.change_x(t)
    assert renamed.x == t and renamed.to_expr() == s.exp(t+2)
    assert h.change_ics(1).y0 == [s.E]
    assert s.simplify(h.integrate((x,0,x)).to_expr()-(s.exp(x)-1)) == 0


def test_hypergeometric_and_meijerg_converters():
    hyper = s.hyper([s.Rational(1,2)],[s.Rational(3,2)],x**2)
    h = from_hyper(hyper)
    assert_annihilates(h,hyper)
    assert h.series(8,coefficient=True) == [1,0,s.Rational(1,3),0,s.Rational(1,10),0,s.Rational(1,42),0]
    g = s.meijerg([],[],[0],[],x)
    h = from_meijerg(g,initcond=False)
    assert_annihilates(h,g)
    assert h.to_meijerg() == g and h.y0 is None
    assert from_meijerg(g).series(4,coefficient=True) == [1,-1,s.Rational(1,2),-s.Rational(1,6)]
    with pytest.raises(TypeError):
        from_hyper(s.sin(x))


def test_parameters_and_overridden_initial_values():
    a = s.Symbol('a')
    h = expr_to_holonomic(s.sin(a*x),x)
    assert_annihilates(h,s.sin(a*x))
    assert h.series(4,order=False) == a*x-a**3*x**3/6
    assert h.to_sequence()[0][0].eval(5) == a**5/120
    hyper = from_hyper(s.hyper([a],[2],x))
    assert hyper.series(3,coefficient=True) == [1,a/2,(a*a+a)/12]
    assert expr_to_holonomic(s.exp(x),x,y0=[2]).to_expr() == 2*s.exp(x)
    original = expr_to_holonomic(s.exp(x),x)
    assert original.composition(2*x,0,[3]).to_expr() == 3*s.exp(2*x)
    assert expr_to_holonomic(1,x).diff().series(3,order=False) == 0
    R,D = DifferentialOperators(s.QQ.old_poly_ring(x),'D')
    sine = HolonomicFunction(D**2+1,x,0,[0,1])
    product = sine*x
    assert product.y0 == [0,0,2]
    assert product.degree() == 2
    assert product.series(6,order=False) == x**2-x**4/6
    assert product.to_sequence()[0][0].eval(4) == -s.Rational(1,6)


@pytest.mark.parametrize('target',[1,-1,s.I/2])
def test_numeric_continuation(target):
    R,D = DifferentialOperators(s.QQ.old_poly_ring(x),'D')
    h = HolonomicFunction(D**2+1,x,0,[0,1])
    actual = h.evalf(target,h=0.01,derivatives=True)
    expected = [s.sin(target),s.cos(target)]
    for value,reference in zip(actual,expected):
        assert float(s.Abs(value-reference).evalf()) < 1e-9
    assert len(h.evalf([0,1])) == 2


def test_failures_are_explicit():
    R,D = DifferentialOperators(s.QQ.old_poly_ring(x),'D')
    with pytest.raises(ValueError):
        DifferentialOperators(s.QQ.poly_ring(x,s.Symbol('y')),'D')
    with pytest.raises(ValueError):
        D**-1
    with pytest.raises(ValueError):
        HolonomicFunction(D-1,x).series()
    with pytest.raises(NotImplementedError):
        expr_to_holonomic(s.tan(x),x)
    with pytest.raises(NotImplementedError):
        expr_to_holonomic(s.exp(x),x).composition(s.sin(x))
    with pytest.raises(ValueError):
        HolonomicFunction(x*D-1,x,0,[1]).series()
    with pytest.raises(NotImplementedError):
        HolonomicFunction(x*D-1,x,0,[0]).series()
    with pytest.raises(ValueError):
        HolonomicFunction(D-1,x,0,[1]).evalf(1,h=0)


def test_no_runtime_reference_dependency():
    subprocess.run([sys.executable,'-c',r'''
import sys
class NoSympy:
 def find_spec(self,name,*args):
  if name == 'sympy' or name.startswith('sympy.'):
   raise RuntimeError('runtime reference import: '+name)
sys.meta_path.insert(0,NoSympy())
import symox as s
from symox.holonomic import *
x=s.Symbol('x')
h=expr_to_holonomic(s.exp(x)*s.sin(x),x)
assert h.series(4,order=False)==x+x**2+x**3/3
assert h.to_sequence()[0][0].eval(4)==0
assert 'sympy' not in sys.modules
'''],check=True)
