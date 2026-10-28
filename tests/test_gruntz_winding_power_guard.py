"""Principal roots cannot treat a winding complex scale as positive."""
import pytest
import sympy as sp
import symox as sx


def data():
    n=sx.Symbol('n',positive=True)
    delta=(n+1)**(1+sx.I)-n**(1+sx.I)
    return n,delta


@pytest.mark.parametrize('form',['ratio','inverse','root_square'])
def test_winding_roots_keep_the_original_undecided_limit(form):
    n,delta=data()
    f={'ratio':sx.sqrt(delta)/sx.sqrt(-delta),
       'inverse':sx.sqrt(-delta)/sx.sqrt(delta),
       'root_square':sx.sqrt(delta**2)/delta}[form]
    result=sx.limit(f,n,sx.oo)
    assert isinstance(result,sx.Limit)
    assert result==sx.Limit(f,n,sx.oo)


def test_two_exact_subsequences_have_distinct_principal_root_limits():
    # delta=n**I*h(1/n), h(t)=((1+t)**(1+I)-1)/t -> 1+I.
    # The derivative certifies h's removable limit on positive t, where log
    # is analytic. For n_j=exp(2*pi*j+offset), n_j**I is exactly exp(I*offset).
    # The two limiting deltas are ±sqrt(2)*I, both away from either root cut.
    t=sp.Symbol('t',positive=True);j=sp.Symbol('j',integer=True,positive=True)
    coefficient=sp.diff((1+t)**(1+sp.I),t).subs(t,0)
    assert coefficient==1+sp.I
    limits=[]
    for offset in [sp.pi/4,5*sp.pi/4]:
        carrier=sp.simplify(sp.expand_complex(sp.exp(sp.I*(2*sp.pi*j+offset))))
        value=sp.simplify(carrier*coefficient)
        limits.append(sp.simplify(sp.expand_complex(sp.sqrt(value)/sp.sqrt(-value))))
    assert limits==[sp.I,-sp.I]


@pytest.mark.parametrize('k',[sx.I,1+sx.I])
def test_complex_exponents_keep_certified_rational_scale_cancellation(k):
    n=sx.Symbol('n',positive=True)
    f=(n+1)**k/((n+1)**(k+1)-n**(k+1))
    result=sx.limit(f,n,sx.oo)
    assert not result.has(sx.Limit)
    assert sx.simplify(result-1/(k+1))==0


@pytest.mark.parametrize('power',[2,3])
def test_integer_powers_of_two_differences_have_no_root_branch_loss(power):
    n,delta=data();other=(n+2)**(1+sx.I)-n**(1+sx.I)
    result=sx.limit(delta**power/other**power,n,sx.oo)
    assert not result.has(sx.Limit)
    assert sx.simplify(result-sx.Rational(1,2)**power)==0


def test_scale_independent_root_does_not_disable_complex_rational_limits():
    n,delta=data();coefficient=sx.sqrt(1+sx.I)
    result=sx.limit(coefficient*(n+1)**sx.I/delta,n,sx.oo)
    assert not result.has(sx.Limit)
    assert sx.simplify(result-coefficient/(1+sx.I))==0


def test_proven_real_scale_keeps_principal_root_normalization():
    n=sx.Symbol('n',positive=True);k=sx.Symbol('k',positive=True)
    value=n**k
    result=sx.limit(sx.sqrt(value)/sx.sqrt(-value),n,sx.oo)
    assert not result.has(sx.Limit)
    assert sx.simplify(result+sx.I)==0


def test_noninteger_power_guard_also_preserves_the_third_root_branches():
    n,delta=data();q=sx.Rational(1,3)
    f=delta**q/(-delta)**q
    assert sx.limit(f,n,sx.oo)==sx.Limit(f,n,sx.oo)
