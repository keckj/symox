"""Correlated Gaussian moments: independent MGF coefficients and protocols."""
import json
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp
import symox as sx
from symox import stats
from symox.stats._joint import _normal_raw


def mgf_coefficient(powers,mean,cov):
    """Coefficient in exp(mu*t) times the separate quadratic exponentials."""
    dimension=len(powers)
    zero=(0,)*dimension
    polynomial={zero:sp.S.One}
    factors=[]
    for i,mu in enumerate(mean):
        factor={}
        for n in range(powers[i]+1):
            exponent=list(zero);exponent[i]=n
            factor[tuple(exponent)]=mu**n/sp.factorial(n)
        factors.append(factor)
    for i in range(dimension):
        for j in range(i,dimension):
            coefficient=cov[i][j]/(2 if i==j else 1)
            limit=powers[i]//2 if i==j else min(powers[i],powers[j])
            factor={}
            for n in range(limit+1):
                exponent=list(zero);exponent[i]+=n;exponent[j]+=n
                factor[tuple(exponent)]=coefficient**n/sp.factorial(n)
            factors.append(factor)
    for factor in factors:
        result={}
        for left,a in polynomial.items():
            for right,b in factor.items():
                exponent=tuple(x+y for x,y in zip(left,right))
                if all(n<=maximum for n,maximum in zip(exponent,powers)):
                    result[exponent]=result.get(exponent,0)+a*b
        polynomial=result
    return sp.expand(polynomial.get(tuple(powers),0)*sp.prod(sp.factorial(n) for n in powers))


def equal(actual,expected):
    assert sp.expand(sx.to_sympy(actual)-expected)==0


@pytest.mark.parametrize('powers',[(0,0),(1,1),(2,2),(4,2),(3,3),(6,4)])
def test_bivariate_raw_moments_match_generating_function(powers):
    mean=sx.ImmutableMatrix([1,-2])
    cov=sx.ImmutableMatrix([[2,sx.Rational(1,3)],[sx.Rational(1,3),3]])
    expected=mgf_coefficient(powers,[sp.Integer(1),sp.Integer(-2)],[[sp.Integer(2),sp.Rational(1,3)],[sp.Rational(1,3),sp.Integer(3)]])
    equal(_normal_raw(powers,mean,cov),expected)


@pytest.mark.parametrize('powers',[(2,2,2),(4,3,2)])
def test_three_correlated_coordinates_match_generating_function(powers):
    mean=sx.ImmutableMatrix([1,-2,3]);cov=sx.ImmutableMatrix([[2,1,0],[1,3,1],[0,1,4]])
    expected=mgf_coefficient(powers,list(map(sp.Integer,[1,-2,3])),[[sp.Integer(x) for x in row] for row in [[2,1,0],[1,3,1],[0,1,4]]])
    equal(_normal_raw(powers,mean,cov),expected)


def test_symbolic_correlations_match_wick_pairings():
    a,b,c=sx.symbols('joint_a joint_b joint_c')
    value=_normal_raw((2,2),[sx.Integer(0)]*2,sx.ImmutableMatrix([[a,c],[c,b]]))
    assert sx.expand(value-(a*b+2*c*c))==0


def test_public_joint_expectation_keeps_correlations():
    X=stats.MultivariateNormal('correlated_native',[1,-2],[[2,sx.Rational(1,3)],[sx.Rational(1,3),3]])
    expected=mgf_coefficient((4,2),[sp.Integer(1),sp.Integer(-2)],[[sp.Integer(2),sp.Rational(1,3)],[sp.Rational(1,3),sp.Integer(3)]])
    equal(stats.E(X[0]**4*X[1]**2),expected)
    assert stats.covariance(X[0],X[1])==sx.Rational(1,3)


def test_dense_kernel_avoids_python_recurrence_cache(monkeypatch):
    import symox.stats._joint as implementation
    def unavailable(*args,**kwargs):
        raise AssertionError('native recurrence entered Python lru_cache')
    monkeypatch.setattr(implementation,'lru_cache',unavailable)
    assert _normal_raw((2,2),[sx.Integer(0)]*2,sx.ImmutableMatrix([[2,1],[1,3]]))==8


def test_high_order_centered_moment_has_no_recursive_call_chain():
    value=_normal_raw((1200,),[sx.Integer(0)],sx.ImmutableMatrix([[1]]))
    equal(value,sp.factorial2(1199))


def test_float_rounding_and_precision_match_frozen_release():
    fixture=json.loads((Path(__file__).parent/'fixtures/normal_raw_float_c461.json').read_text())
    for case in fixture['cases']:
        digits=case['digits']
        mean=sx.ImmutableMatrix([sx.Float('1.125',digits),sx.Float('-0.75',digits)])
        cov=sx.ImmutableMatrix([[sx.Float('2.125',digits),sx.Float('0.2',digits)],[sx.Float('0.2',digits),sx.Float('3.5',digits)]])
        value=_normal_raw(tuple(case['powers']),mean,cov)
        assert str(value)==case['value']
        assert getattr(value,'_prec',None)==case['precision']


def test_float_accuracy_is_checked_independently_at_higher_precision():
    digits=80
    mean=sx.ImmutableMatrix([sx.Float('1.125',digits),sx.Float('-0.75',digits)])
    cov=sx.ImmutableMatrix([[sx.Float('2.125',digits),sx.Float('0.2',digits)],[sx.Float('0.2',digits),sx.Float('3.5',digits)]])
    powers=(6,4);value=_normal_raw(powers,mean,cov)
    exact_mean=[sp.Rational(sx.to_sympy(v)) for v in mean]
    exact_cov=[[sp.Rational(sx.to_sympy(cov[i,j])) for j in range(2)] for i in range(2)]
    reference=mgf_coefficient(powers,exact_mean,exact_cov)
    with mp.workdps(digits+50):
        expected=mp.mpf(str(reference.evalf(digits+50)))
        actual=mp.mpf(str(value))
        assert abs(actual-expected)<=abs(expected)*mp.mpf(10)**(3-digits)


def test_matrix_subclass_indexing_remains_observable():
    class CustomCovariance(sx.MutableDenseMatrix):
        def __getitem__(self,key):
            raise RuntimeError('custom covariance indexing')
    with pytest.raises(RuntimeError,match='custom covariance indexing'):
        _normal_raw((2,),[sx.Integer(0)],CustomCovariance([[1]]))


def test_modified_standard_matrix_getter_uses_original_protocol(monkeypatch):
    cov=sx.ImmutableMatrix([[1]])
    def unavailable(self,key):
        raise RuntimeError('patched matrix getter')
    monkeypatch.setattr(type(cov),'__getitem__',unavailable)
    with pytest.raises(RuntimeError,match='patched matrix getter'):
        _normal_raw((2,),[sx.Integer(0)],cov)


def test_host_and_function_parameters_retain_fallback():
    class Mean(sx.Symbol):
        pass
    host=Mean('host_joint_mean',real=True)
    cov=sx.ImmutableMatrix([[1]])
    assert sx._core._stats_normal_raw([2],[host],cov) is None
    x=sx.Symbol('joint_parameter_x',real=True)
    assert sx._core._stats_normal_raw([2],[sx.sin(x)],cov) is None
    equal(_normal_raw((2,),[sx.sin(x)],cov),1+sp.sin(sx.to_sympy(x))**2)


def test_modified_final_simplifier_keeps_interruptions(monkeypatch):
    def interrupted(value):
        raise KeyboardInterrupt
    monkeypatch.setattr(sx,'simplify',interrupted)
    with pytest.raises(KeyboardInterrupt):
        _normal_raw((2,),[sx.Integer(0)],sx.ImmutableMatrix([[1]]))


def test_invalid_and_lazy_shapes_keep_existing_results_and_errors():
    cov=sx.ImmutableMatrix([[1]])
    assert _normal_raw((1,0),[sx.Integer(4)],cov)==4
    assert sx._core._stats_normal_raw([1,0],[4],cov) is None
    with pytest.raises(IndexError):
        _normal_raw((0,1),[sx.Integer(4)],cov)
    assert _normal_raw((0,0),object(),object())==1


def test_unsupported_orders_and_state_budget_are_not_zero_moments():
    cov=sx.ImmutableMatrix([[1]])
    for powers in [[-1],[sx.Symbol('n')],[sx.Float(2.0)],[2**40]]:
        assert sx._core._stats_normal_raw(powers,[0],cov) is None
    assert sx._core._stats_normal_raw([1000,1000],[0,0],sx.eye(2)) is None


@pytest.mark.parametrize('order,value,precision',[(sx.Integer(2),'7',None),(sx.Float(2.0),'7.00000000000000',53),(sx.Float(0.0),'1',None)])
def test_non_builtin_concrete_orders_preserve_existing_rounding(order,value,precision):
    result=_normal_raw((order,),[sx.Integer(2)],sx.ImmutableMatrix([[3]]))
    assert str(result)==value
    assert getattr(result,'_prec',None)==precision


def test_abstract_covariance_and_held_parameters_keep_fallback():
    cov=sx.MatrixSymbol('abstract_joint_covariance',2,2)
    assert sx._core._stats_normal_raw([2,2],[0,0],cov) is None
    assert sx.expand(_normal_raw((2,2),[0,0],cov)-cov[0,0]*cov[1,1]-2*cov[0,1]**2)==0
    x=sx.Symbol('held_joint_mean')
    mean=sx.Add(x,1,evaluate=False)
    dense=sx.ImmutableMatrix([[1]])
    assert sx._core._stats_normal_raw([2],[mean],dense) is None
    equal(_normal_raw((2,),[mean],dense),(sx.to_sympy(x)+1)**2+1)
