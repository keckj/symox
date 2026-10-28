"""Native independent polynomial moments retain complete general fallbacks."""
import pytest
import sympy as sp
import symox as sx
from symox import stats


def equal(a,b):
    if isinstance(a,sx.Basic):a=sx.to_sympy(a)
    assert sp.simplify(sp.expand(a-b)) == 0


@pytest.mark.parametrize('n',[0,1,4,12])
@pytest.mark.parametrize('order',[0,1,2,3,5,7])
def test_binomial_raw_moments_equal_the_whole_finite_sum(n,order):
    p=sx.Symbol('p',positive=True)
    variable=stats.Binomial('native_binomial_%s_%s'%(n,order),n,p,2,-3)
    P=sp.Symbol('p',positive=True)
    expected=sum((2*j-3*(n-j))**order*sp.binomial(n,j)*P**j*(1-P)**(n-j) for j in range(n+1))
    equal(variable.distribution.moment(order),expected)
    equal(stats.E(variable**order),expected)


@pytest.mark.parametrize('order',range(9))
def test_normal_raw_moments_match_gaussian_generating_function(order):
    mu=sx.Symbol('native_mu',real=True);sigma=sx.Symbol('native_sigma',positive=True)
    x=stats.Normal('native_normal',mu,sigma)
    M=sp.Symbol('native_mu',real=True);S=sp.Symbol('native_sigma',positive=True)
    # Independent combinatorial coefficients of exp(mu*t+sigma²*t²/2).
    expected=sum(sp.factorial(order)*M**(order-2*j)*S**(2*j)/
                 (sp.factorial(order-2*j)*2**j*sp.factorial(j)) for j in range(order//2+1))
    equal(x.distribution.moment(order),expected)
    equal(stats.E(x**order),expected)


def test_symbolic_binomial_variance_is_native_and_exact(monkeypatch):
    import symox.stats._queries as queries
    n=sx.Symbol('native_n',integer=True,nonnegative=True)
    p=sx.Symbol('native_p',nonnegative=True)
    x=stats.Binomial('native_symbolic_binomial',n,p,3,-2)
    def forbidden(*args,**kwargs):raise AssertionError('native polynomial query entered Python enumeration')
    monkeypatch.setattr(queries,'_worlds',forbidden)
    N,P=sp.symbols('native_n native_p',integer=True)
    # The outcome is -2*n+5*Binomial(n,p).
    actual=stats.variance(x).subs({n:sx.Symbol('N'),p:sx.Symbol('P')})
    equal(actual,25*sp.Symbol('N')*sp.Symbol('P')*(1-sp.Symbol('P')))
    assert stats.E(x).has(x._atom) is False


def test_independent_polynomial_moment_products_match_direct_joint_sum(monkeypatch):
    import symox.stats._queries as queries
    x=stats.Binomial('native_joint_x',3,sx.Rational(1,3))
    y=stats.Die('native_joint_y',4)
    def forbidden(*args,**kwargs):raise AssertionError('native polynomial query entered Python enumeration')
    monkeypatch.setattr(queries,'_worlds',forbidden)
    expected=sum(sp.binomial(3,i)*sp.Rational(1,3)**i*sp.Rational(2,3)**(3-i)*
                 sp.Rational(1,4)*(i*j+i**2+2)**3 for i in range(4) for j in range(1,5))
    equal(stats.E((x*y+x**2+2)**3),expected)


def test_mutated_finite_table_keeps_its_actual_mass_and_coefficients():
    x=stats.Binomial('native_mutated',3,sx.Rational(1,2))
    original=x.distribution._table
    try:
        x.distribution._table={sx.Integer(0):sx.Integer(1),sx.Integer(1):sx.Integer(1)}
        assert stats.E(x+1)==3
        assert stats.E(x)==1
    finally:x.distribution._table=original


def test_custom_moment_callback_is_not_bypassed():
    x=stats.Normal('native_override',0,1)
    original=x.distribution._moment
    try:
        x.distribution._moment=lambda n:sx.Integer(7)
        assert stats.E(x)==7
        assert stats.variance(x)==-42
    finally:x.distribution._moment=original


def test_conditions_and_held_queries_retain_their_existing_paths():
    x=stats.Binomial('native_conditioned',3,sx.Rational(1,2))
    assert stats.E(x,x>1)==sx.Rational(9,4)
    assert stats.variance(x,x>1)==sx.Rational(3,16)
    assert stats.variance(x,evaluate=False).doit()==sx.Rational(3,4)


def test_nonpolynomial_expectation_uses_the_existing_complete_finite_route():
    x=stats.Die('native_nonpolynomial',4)
    equal(stats.E(sx.sin(x)),sum(sp.sin(j) for j in range(1,5))/4)


@pytest.mark.parametrize('precision',[53,140])
def test_float_moments_keep_precision_and_the_general_route(precision):
    p=sx.Float('0.1234567890123456789',precision=precision)
    x=stats.Binomial('native_float_%s'%precision,4,p)
    actual=stats.variance(x)
    table={sp.Integer(k):sp.Rational(sx.to_sympy(v)) for k,v in x.distribution.dict.items()}
    exact=sum(k*k*v for k,v in table.items())-sum(k*v for k,v in table.items())**2
    assert abs(sp.Rational(sx.to_sympy(actual))-exact) < sp.Rational(2)**(-(precision-10))


def test_unsupported_native_query_is_distinct_from_a_zero_expectation():
    x=sx.Symbol('native_scalar_x')
    assert sx._core._stats_polynomial_query(sx.sin(x),[x],[('Finite',(),[1,2],[sx.Rational(1,2)]*2)],False) is None
    assert sx._core._stats_polynomial_query(x,[x],[('Finite',(),[-1,1],[sx.Rational(1,2)]*2)],False)==0


def test_large_degree_finite_query_keeps_the_complete_nonexpanding_fallback():
    x=stats.Bernoulli('native_large_degree_x',sx.Rational(1,2))
    y=stats.Bernoulli('native_large_degree_y',sx.Rational(1,2))
    assert stats.E((x+y)**1000)==sx.Rational(2+2**1000,4)


def test_known_nonfinite_parameter_retains_the_general_model_route():
    x=sx.Symbol('native_query_atom')
    value=sx.Symbol('native_nonfinite_value',finite=False)
    assert sx._core._stats_polynomial_query(x,[x],[('Binomial',[1,sx.Rational(1,2),value,0],(),())],False) is None


def test_lazy_finite_table_is_not_consumed_by_a_speculative_query():
    x=stats.Binomial('native_lazy_table',3,sx.Rational(1,2))
    original=x.distribution._table;calls=[]
    def table():
        calls.append(1)
        return {sx.Integer(1):sx.Integer(1)}
    try:
        x.distribution._table=table
        assert stats.E(sx.sin(x))==sx.sin(1)
        assert calls==[1]
    finally:x.distribution._table=original


def test_overridden_table_property_keeps_its_probability_space_protocol():
    from symox.stats._base import Distribution
    class Custom(Distribution):
        @property
        def dict(self):return {sx.Integer(1):sx.Integer(1)}
    x=stats.Binomial('native_custom_table',3,sx.Rational(1,2))
    original=type(x.distribution)
    try:
        x.distribution.__class__=Custom
        assert stats.E(x)==1
    finally:x.distribution.__class__=original


def test_held_repeated_random_factors_keep_every_exponent():
    x=stats.Die('native_held_repeated',6)
    assert stats.E(sx.Mul(x,x,evaluate=False))==sx.Rational(91,6)
    assert stats.E(sx.Mul(x,x,x,evaluate=False))==sx.Rational(147,2)
