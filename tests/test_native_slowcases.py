"""Nested exponential solve: same solution sets as SymPy 1.14, and the
algebraicity answers the reordered structural tests must keep."""
import time

import pytest
import symox as sx

x=sx.Symbol('x')


def tower(n,shift=0):
    e=x
    for _ in range(n):e=sx.exp(e)+shift
    return e


# Expected values printed by SymPy 1.14.0.
@pytest.mark.parametrize('equation,expected',[
    (lambda:tower(3,1)-5,'[log(1 - log(-1 + log(4))) + I*pi]'),
    (lambda:tower(4)-2,'[log(log(-log(log(2))) + I*pi)]'),
    (lambda:tower(3)-3,'[log(log(log(3)))]'),
    (lambda:sx.exp(sx.exp(x))-sx.exp(2),'[log(2)]'),
])
def test_nested_exp_matches_sympy(equation,expected):
    assert str(sx.solve(equation(),x))==expected


# SymPy 1.14 prints these with log(I*pi) and 2*log(2); symox keeps the
# equal forms log(pi) + I*pi/2 and log(4) (unchanged by this wave).
# Certificate: exp(log(z)) == z makes the substituted equation exactly 0.
@pytest.mark.parametrize('equation',[
    lambda:sx.exp(sx.exp(sx.exp(x)))+1,
    lambda:tower(5)+1,
    lambda:tower(4,1)-5,
    lambda:tower(6,1)-5,
    lambda:tower(7)-2,
])
def test_nested_exp_single_exact_root(equation):
    e=equation()
    roots=sx.solve(e,x)
    assert len(roots)==1
    assert e.subs(x,roots[0])==0


def test_nested_exp_depth_is_not_exponential():
    start=time.perf_counter()
    roots=sx.solve(tower(7)-2,x)
    assert len(roots)==1
    assert time.perf_counter()-start<20


@pytest.mark.parametrize('expression,expected',[
    (lambda:sx.exp(sx.sqrt(2)),False),
    (lambda:sx.exp(sx.I*sx.pi/3),True),
    (lambda:sx.exp(sx.E),None),
    (lambda:sx.exp(2*sx.I),False),
    (lambda:sx.sin(sx.sqrt(3)),False),
    (lambda:sx.cos(sx.log(2)),None),
    (lambda:sx.log(1+sx.sqrt(2)),False),
    (lambda:sx.log(1+sx.E),None),
    (lambda:sx.sqrt(sx.E),False),
    (lambda:sx.Pow(sx.pi,sx.Rational(1,3)),False),
    (lambda:sx.sqrt(1+sx.sqrt(2)),True),
    (lambda:sx.exp(1+sx.I*sx.pi),False),
])
def test_algebraic_answers_match_sympy(expression,expected):
    assert expression().is_algebraic is expected


def test_reducible_quintic_log_part_is_split_by_factor():
    # SymPy 1.14 does not finish this integral in 1800 s (memory grows in the
    # general log_to_real substitution); the estimate gate splits q(t).
    f=1/(x**5+x+1)
    r=sx.integrate(f,x)
    assert not r.has(sx.RootSum)
    d=sx.diff(r,x)-f
    for point in (sx.Rational(1,3),sx.Rational(-5,7),sx.Integer(2)):
        assert abs(sx.N(d.subs(x,point),30))<sx.Float('1e-25')


def test_substitution_estimate_below_threshold_keeps_general_path():
    from symox.integrals import rationaltools as rt
    u,v=sx.symbols('u v')
    a=u**4+v**4*x+u*v
    assert rt._lrt_substitution_log10(a,u,v,sx.Rational(1,2),sx.sqrt(2)/2)<rt._LRT_SPLIT_LOG10
    big=sx.sqrt(sum(sx.sqrt(k)*sx.Rational(k,7) for k in range(2,200)))
    assert rt._lrt_substitution_log10(a,u,v,sx.Rational(1,2),big)>rt._LRT_SPLIT_LOG10
