"""Native core protocol behavior, arbitrary precision and independent probes."""
from fractions import Fraction
import math
import mpmath
import pytest
import symox as S
import sympy as R


@pytest.mark.parametrize('name,expression',[('Zero','0'),('One','1'),('NegativeOne','-1'),('Half','Rational(1,2)'),('Infinity','oo'),('NegativeInfinity','-oo'),('ComplexInfinity','zoo'),('NaN','nan'),('ImaginaryUnit','I'),('Pi','pi'),('Exp1','E')])
def test_singleton_number_constructors_and_type_protocol(name,expression):
    from symox.core import numbers
    constructor=getattr(numbers,name);value=S.sympify(expression)
    assert constructor()==value and constructor() is constructor()
    assert isinstance(value,constructor)
    assert not isinstance(S.Symbol('z'),constructor)


@pytest.mark.parametrize('value',[0,1,-1,1.0,0.5,1.25,float('inf')])
def test_numeric_integer_valued_and_exact_float_comparison(value):
    from symox.core.numbers import int_valued,equal_valued
    from sympy.core.numbers import int_valued as reference
    assert int_valued(value)==reference(value)
    if math.isfinite(value):assert int_valued(S.Float(value))==reference(R.Float(value))
    assert equal_valued(S.Float(.5),S.Rational(1,2))
    assert not equal_valued(S.Float(.1),S.Rational(1,10))
    assert not equal_valued(S.Float('1.000000000000000000000000000000000000001',60),1)


def test_numeric_conversions_preserve_high_precision():
    from symox.core.numbers import sympify_mpmath,sympify_fractions,all_close
    context=mpmath.mp.clone();context.dps=70;value=context.sqrt(2)
    native=sympify_mpmath(value)
    assert abs(context.mpf(str(native.evalf(65)))-value)<context.mpf('1e-64')
    assert sympify_fractions(Fraction(11,17))==S.Rational(11,17)
    x=S.Symbol('x')
    assert all_close(2*x+S.Float('0.5'),2*x+S.Rational(1,2))
    assert not all_close(2*x,3*x)


def test_core_fact_database_and_managed_properties():
    from symox.core.assumptions import StdFactKB,ManagedProperties
    facts=StdFactKB({'integer':True,'positive':True})
    assert facts['rational'] and facts['real'] and facts['nonzero'] and not facts['negative']
    with pytest.raises(ValueError):facts.deduce_all_facts({'integer':False})
    # Reference: the deprecated metaclass prepares the class assumptions; the
    # instance storage ``_assumptions`` is what Basic.__new__ provides.
    with pytest.warns(S.utilities.exceptions.SymPyDeprecationWarning):
        class Example(metaclass=ManagedProperties):
            def _eval_is_integer(self):return True
        class Stored(Example):
            def __init__(self):self._assumptions=type(self).default_assumptions
    with pytest.raises(AttributeError):Example().is_rational
    value=Stored()
    assert value.is_rational and value.is_real


def test_core_builders_containers_and_unique_names():
    from symox.core.expr import ExprBuilder,unchanged
    from symox.core.containers import OrderedSet,TupleKind
    from symox.core.symbol import Str,uniquely_named_symbol,disambiguate
    from symox.core.function import Application,AppliedUndef,DefinedFunction,UndefinedFunction
    x,y=S.symbols('x y');builder=ExprBuilder(S.Add,[x,ExprBuilder(S.Mul,[2,y])])
    assert builder.build()==x+2*y
    assert unchanged(S.sin,x) and not unchanged(S.cos,0)
    values=OrderedSet([3,1,3,2]);assert list(values)==[3,1,2]
    assert values.pop(last=False)==3
    assert TupleKind(S.NumberKind,S.BooleanKind) is TupleKind(S.NumberKind,S.BooleanKind)
    assert Str('x').free_symbols==frozenset() and Str('x')!=x
    assert str(uniquely_named_symbol('x',[x,S.Symbol('x0')]))=='x1'
    a,b=S.Symbol('x'),S.Symbol('x',integer=True);separated=disambiguate(a,b)
    assert len(separated.free_symbols)==2
    f=UndefinedFunction('custom');assert isinstance(f,UndefinedFunction)
    assert isinstance(f(x),AppliedUndef) and isinstance(S.sin(x),DefinedFunction) and isinstance(S.sin(x),Application)


def test_core_factor_algebra_and_traversal():
    from symox.core.exprtools import Factors,Term,decompose_power,decompose_power_rat
    from symox.core.traversal import iterargs,walk
    from symox.core.multidimensional import apply_on_element
    x,y=S.symbols('x y');a,b=Factors(x**3*y),Factors(x*y**2)
    assert a.gcd(b).as_expr()==x*y
    assert a.lcm(b).as_expr()==x**3*y**2
    assert a.quo(b).as_expr()==x**2  # SymPy: Factors({x: 2}); rem is Factors({y: 1})
    assert (Term(2*x/(x+1))*Term((x+1)/y)).as_expr()==2*x/y
    for expression in [x,x**2,x**S.Rational(3,2),S.exp(2*y/3),x**(3*y/2)]:
        base,power=decompose_power(expression)
        assert S.simplify(base**power-expression)==0
    assert list(iterargs(S.sin(x+y)))==[S.sin(x+y),x+y,x,y]
    assert apply_on_element(lambda a,b:a+b,[[1,[2,3]],10],{},0)==[11,[12,13]]


def test_numeric_protocol_bits_and_complex_roundtrips():
    from symox.core.evalf import evalf,quad_to_mpmath,complex_accuracy,add_terms,scaled_zero,iszero,check_target,check_convergence,as_mpmath
    context=mpmath.mp.clone();context.prec=170
    for native,reference in [(S.pi,context.pi),(S.Rational(1,7),context.mpf(1)/7),(S.exp(S.I),context.exp(context.j))]:
        result=evalf(native,150,{})
        assert abs(quad_to_mpmath(result,context)-reference)<context.mpf(2)**-145
        assert complex_accuracy(result)>=150
    assert evalf(S.Integer(0),100,{})==(None,None,None,None)
    zero,accuracy=add_terms([(context.mpf(1)._mpf_,100),(context.mpf(-1)._mpf_,100)],100,80)
    assert iszero(zero,scaled=True) and accuracy==-1
    z,p=scaled_zero(100,-1);assert scaled_zero(z)==(1,1,100,1)
    with pytest.raises(S.PrecisionExhausted):check_target(S.pi,(context.pi._mpf_,None,2,None),100)
    n=S.Symbol('n');assert check_convergence(n,n+2,n)==(0,1,2)
    assert abs(as_mpmath(S.sqrt(2),100,{})**2-2)<context.mpf(2)**-95


def test_reproducible_numeric_probes_and_shared_singletons():
    from symox.core.random import seed,random_complex_number,verify_numerically,test_derivative_numerically
    from symox.core.singleton import Singleton,SingletonRegistry
    seed(123);a=random_complex_number();seed(123);assert random_complex_number()==a
    x=S.Symbol('x');assert verify_numerically(S.sin(x)**2+S.cos(x)**2,1,x)
    assert not verify_numerically(S.sin(x),S.cos(x),x)
    assert test_derivative_numerically(S.exp(x),x)
    class LocalSingleton(metaclass=Singleton):pass
    assert LocalSingleton() is LocalSingleton() and S.S.LocalSingleton is LocalSingleton()
    assert isinstance(S.S,SingletonRegistry)
