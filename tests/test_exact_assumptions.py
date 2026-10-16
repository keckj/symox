"""Exact algebraic assumptions and real nonzero versus complex not-zero."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


@pytest.mark.parametrize("source",[
    "sqrt(2)","2**Rational(1,3)","2**Rational(2,3)",
    "Rational(8,27)**Rational(1,3)","Rational(4,9)**Rational(-1,2)",
    "Rational(1,9)**Rational(-1,2)","(-2)**Rational(1,3)",
    "(-8)**Rational(2,3)","sqrt(-2)","sqrt(2)/2","1+sqrt(2)",
    "sqrt(2)-1","sqrt(2)*sqrt(3)","I","sqrt(2)*I","1+sqrt(2)*I",
])
def test_rational_power_predicates_match_sympy(source):
    actual,expected=sm.sympify(source),sp.sympify(source)
    for name in ["is_integer","is_rational","is_real","is_zero","is_nonzero"]:
        assert getattr(actual,name)==getattr(expected,name),(source,name)


@pytest.mark.parametrize("flags",[
    {"nonzero":True},{"nonzero":False},{"zero":False},{"imaginary":True},
    {"nonzero":False,"zero":False},{"real":True,"nonzero":False},
    {"real":True,"zero":False},{"rational":False},{"real":False},
])
def test_symbol_nonzero_flags_and_roundtrips(flags):
    x,X=sm.Symbol("x",**flags),sp.Symbol("x",**flags)
    for name in ["is_real","is_zero","is_nonzero","is_integer","is_rational"]:
        assert getattr(x,name)==getattr(X,name),(flags,name)
    assert pickle.loads(pickle.dumps(x))==x
    assert sm.sympify(x.srepr())==x
    converted=sm.to_sympy(x)
    for name in ["is_real","is_zero","is_nonzero"]:
        assert getattr(converted,name)==getattr(X,name)
    assert sm.from_sympy(converted)==x


def test_inconsistent_nonzero_and_imaginary_assumptions_reject():
    for flags in [{"nonzero":True,"imaginary":True},{"zero":False,"real":True,"nonzero":False},
                  {"rational":False,"integer":True}]:
        # The second case is real, not zero, and not nonzero simultaneously.
        with pytest.raises(ValueError):sm.Symbol("x",**flags)


@pytest.mark.parametrize("flags",[
    {"integer":True},{"rational":True},{"nonzero":True},{"nonzero":False},
    {"zero":False},{"imaginary":True},{"rational":False},{"real":True},
])
@pytest.mark.parametrize("source",["x**2-2","x**3-2","x**2+2","x**2-x","(x-1)*(x**2-2)"])
def test_solver_exact_algebraic_candidate_filtering(flags,source):
    x,X=sm.Symbol("x",**flags),sp.Symbol("x",**flags)
    equation=sm.sympify(source).xreplace({sm.Symbol("x"):x})
    oracle=sp.sympify(source).xreplace({sp.Symbol("x"):X})
    actual=sm.solve(equation,x)
    expected=sp.solve(oracle,X)
    remaining=list(map(sm.to_sympy,actual))
    for root in expected:
        found=next((i for i,value in enumerate(remaining) if sp.simplify(value-root)==0),None)
        assert found is not None,(source,flags,remaining,expected)
        remaining.pop(found)
    assert not remaining,(source,flags,remaining,expected)


@pytest.mark.parametrize("source",["-sqrt(2)","2-sqrt(8)","-sqrt(2)-1","-100+sqrt(2)"])
def test_quadratic_sign_filter_respects_generator_embedding(source):
    candidate=sm.sympify(source)
    negative=sm.Symbol("n",negative=True)
    positive=sm.Symbol("p",positive=True)
    assert sm.solve(negative-candidate,negative)==[candidate]
    assert sm.solve(positive-candidate,positive)==[]


def test_near_integer_and_near_zero_candidates_use_exact_proofs():
    integer=sm.Symbol("n",integer=True)
    huge=sm.sqrt(sm.Integer(10)**60+1)
    assert sm.solve(integer-huge,integer)==[]
    near_zero=huge-sm.Integer(10)**30
    positive=sm.Symbol("p",positive=True)
    negative=sm.Symbol("n",negative=True)
    assert sm.solve(positive-near_zero,positive)==[near_zero]
    assert sm.solve(negative-near_zero,negative)==[]
    unrestricted=sm.Symbol("z",zero=False)
    assert sm.solve(unrestricted*unrestricted+1,unrestricted)==[-sm.I,sm.I]


def test_rational_factors_require_nonzero_proof_before_irrational_inference():
    r=sm.Symbol("r",rational=True)
    nonzero=sm.Symbol("r",rational=True,nonzero=True)
    assert (r*sm.sqrt(2)).is_rational is None
    assert (nonzero*sm.sqrt(2)).is_rational is False
    assert (sm.sqrt(2)+sm.sqrt(3)).is_rational is None


def test_complex_not_zero_internal_consumers_remain_valid_without_sympy():
    x,y=sm.symbols("x y")
    assert sm.expand_power_exp(sm.I**(x+y))==sm.I**x*sm.I**y
    assert sm.I.is_nonzero is False and sm.I.is_zero is False
    image=sm.ImageSet(sm.Lambda(x,sm.I*x),sm.Reals)
    assert image.contains(2*sm.I) == sm.true
    assert image.contains(2) == sm.false
    subprocess.run([sys.executable,"-c","""
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x',integer=True)
assert s.solve(x*x-2,x)==[]
assert s.sqrt(2).is_rational is False
assert s.I.is_nonzero is False
z=s.Symbol('z',zero=False)
assert z.is_real is None and z.is_zero is False
"""],check=True,timeout=30)
