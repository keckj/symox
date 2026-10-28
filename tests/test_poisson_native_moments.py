"""Touchard moments and independent Poisson polynomials use native arithmetic."""
import json
from pathlib import Path

import pytest
import sympy as sp
import symox as sx
from symox import stats


def exact_equal(actual,expected):
    assert sp.expand(sx.to_sympy(actual)-expected)==0


@pytest.mark.parametrize('order',[0,1,2,5,12,24])
@pytest.mark.parametrize('symbolic',[False,True])
def test_raw_poisson_moment_matches_touchard_polynomial(order,symbolic):
    rate=sx.Symbol('touchard_rate',nonnegative=True) if symbolic else sx.Rational(3,2)
    model=stats.density(stats.Poisson('touchard',rate))
    exact_equal(model.moment(sx.Integer(order)),sp.bell(order,sx.to_sympy(rate)))


@pytest.mark.parametrize('order',[0,1,2,8])
def test_zero_rate_keeps_degenerate_probability_mass(order):
    X=stats.Poisson('zero_rate',0)
    expected=1 if order==0 else 0
    assert stats.density(X).moment(order)==expected
    assert sx._core._stats_poisson_moment(sx.Integer(0),sx.Integer(order))==expected
    assert stats.E(X**order)==expected
    assert stats.variance(X)==0


def test_raw_moment_does_not_call_python_stirling(monkeypatch):
    rate=sx.Symbol('native_touchard_rate',positive=True)
    model=stats.density(stats.Poisson('native_touchard',rate))
    def unavailable(*args,**kwargs):
        raise AssertionError('Poisson moments must generate one native Stirling row')
    monkeypatch.setattr(sx,'stirling',unavailable)
    exact_equal(model.moment(20),sp.bell(20,sx.to_sympy(rate)))


@pytest.mark.parametrize('order',[2,4,6])
def test_independent_poisson_sum_matches_convolution(order):
    a,b=sx.symbols('poisson_rate_a poisson_rate_b',positive=True)
    X,Y=stats.Poisson('independent_poisson_X',a),stats.Poisson('independent_poisson_Y',b)
    from symox.stats._queries import _native_polynomial_query,random_symbols
    expression=(X+Y)**order
    assert _native_polynomial_query(expression,random_symbols(expression)) is not None
    exact_equal(stats.E(expression),sp.bell(order,sx.to_sympy(a+b)))
    assert stats.variance(X+Y)==a+b
    assert stats.variance(2*X)==4*a


@pytest.mark.parametrize('order',[2,4,6])
def test_poisson_normal_polynomials_match_independent_exact_moments(order):
    X,N=stats.Poisson('mixed_poisson',2),stats.Normal('mixed_normal',0,1)
    # Normal even moments are (j-1)!!; odd moments vanish.
    expected=sum(sp.binomial(order,j)*sp.bell(order-j,2)*sp.factorial2(j-1)
                 for j in range(0,order+1,2))
    exact_equal(stats.E((X+N)**order),expected)


def test_custom_moment_callback_keeps_general_route():
    X=stats.Poisson('custom_poisson_callback',2)
    model=stats.density(X)
    calls=[]
    def custom(order):
        calls.append(order)
        return sx.Integer(100)+order
    model._moment=custom
    from symox.stats._queries import _native_polynomial_query,random_symbols
    assert _native_polynomial_query(X,random_symbols(X)) is None
    assert stats.E(X)==101
    assert calls


def test_overridden_distribution_method_keeps_general_route(monkeypatch):
    X=stats.Poisson('overridden_poisson_method',2)
    model=stats.density(X)
    monkeypatch.setattr(type(model),'moment',lambda self,order:sx.Integer(200)+order)
    from symox.stats._queries import _native_polynomial_query,random_symbols
    assert _native_polynomial_query(X,random_symbols(X)) is None
    assert stats.E(X)==201


def test_random_rate_is_not_declared_independent():
    Y=stats.Poisson('hierarchical_rate',2)
    X=stats.Poisson('hierarchical_poisson',Y)
    from symox.stats._queries import _native_polynomial_query,random_symbols
    assert _native_polynomial_query(X+Y,random_symbols(X+Y)) is None
    atoms=[sx.Symbol('generator_x'),sx.Symbol('generator_y')]
    models=[('Poisson',(atoms[1],),(),()),('Poisson',(sx.Integer(2),),(),())]
    assert sx._core._stats_polynomial_query(atoms[0]+atoms[1],atoms,models,False) is None


def test_host_rate_retains_python_formula_and_interrupts(monkeypatch):
    class Rate(sx.Symbol):
        pass
    rate=Rate('host_poisson_rate',positive=True)
    assert sx._core._stats_poisson_moment(rate,sx.Integer(3)) is None
    model=stats.density(stats.Poisson('host_poisson',rate))
    def interrupted(*args,**kwargs):
        raise KeyboardInterrupt
    monkeypatch.setattr(sx,'stirling',interrupted)
    with pytest.raises(KeyboardInterrupt):
        model.moment(3)


@pytest.mark.parametrize('source',['Symbol("n",integer=True,nonnegative=True)','Float(2.0)','Float(0.0)','Integer(-1)','Integer(4097)'])
def test_ineligible_orders_select_existing_fallback(source):
    order=eval(source,vars(sx))
    assert sx._core._stats_poisson_moment(sx.Integer(2),order) is None
    if not order.is_Integer:
        from symox.stats._univariate import _poisson_moment
        assert _poisson_moment(order,sx.Integer(2)) is None


def test_floating_rate_rounding_matches_previous_release():
    fixture=json.loads((Path(__file__).parent/'fixtures/poisson_float_moments_c461.json').read_text())
    for case in fixture['cases']:
        rate=sx.Float(case['rate'],case['digits'])
        model=stats.density(stats.Poisson('floating_poisson',rate))
        value=model.moment(case['order'])
        assert str(value)==case['value']
        assert getattr(value,'_prec',None)==case['precision']


def test_invalid_rate_validation_is_unchanged():
    with pytest.raises(ValueError,match='nonnegative'):
        stats.Poisson('invalid_poisson_rate',-1)



def test_host_parameter_free_symbols_is_not_inspected_speculatively():
    class Rate(sx.Symbol):
        @property
        def free_symbols(self):
            raise AssertionError('host parameter protocol inspected during native eligibility')
    rate=Rate('opaque_poisson_rate',positive=True)
    X=stats.Poisson('opaque_poisson',rate)
    from symox.stats._queries import _native_polynomial_query,random_symbols
    assert _native_polynomial_query(X,random_symbols(X)) is None
