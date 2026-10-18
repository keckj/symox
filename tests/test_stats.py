"""Probability checks and independent references for the public statistics API."""
import importlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest
import symox as S
import symox.stats as st
import sympy as R
import sympy.stats as ref

Q = S.Rational

C={
'Arcsin':(0,2),'Benini':(2,3,1),'Beta':(2,3),'BetaNoncentral':(2,3,2),'BetaPrime':(2,5),'BoundedPareto':(2,1,4),'Cauchy':(0,2),'Chi':(3,),'ChiNoncentral':(3,2),'ChiSquared':(4,),'Dagum':(2,5,1),'Davis':(1,5,1),'Erlang':(3,2),'ExGaussian':(1,2,3),'Exponential':(2,),'ExponentialPower':(0,2,3),'FDistribution':(5,8),'FisherZ':(5,8),'Frechet':(5,2,0),'Gamma':(3,2),'GammaInverse':(5,2),'Gompertz':(2,3),'Gumbel':(2,1),'Kumaraswamy':(2,3),'Laplace':(1,2),'Levy':(0,2),'Logistic':(1,2),'LogCauchy':(0,2),'LogLogistic':(2,5),'LogitNormal':(0,2),'LogNormal':(1,2),'Lomax':(5,2),'Moyal':(1,2),'Maxwell':(2,),'Nakagami':(2,3),'Normal':(1,2),'GaussianInverse':(2,3),'Pareto':(1,5),'PowerFunction':(2,0,2),'QuadraticU':(0,2),'RaisedCosine':(0,2),'Rayleigh':(2,),'Reciprocal':(1,4),'StudentT':(5,),'ShiftedGompertz':(2,3),'Trapezoidal':(0,1,2,3),'Triangular':(0,2,1),'Uniform':(0,2),'UniformSum':(3,),'VonMises':(0,2),'Wald':(2,3),'Weibull':(2,3),'WignerSemicircle':(2,),
'FlorySchulz':(S.Rational(1,3),),'Geometric':(S.Rational(1,3),),'Hermite':(2,3),'Logarithmic':(S.Rational(1,3),),'NegativeBinomial':(3,S.Rational(1,3)),'Poisson':(2,),'Skellam':(2,3),'YuleSimon':(5,),'Zeta':(5,)}


def same(actual,expected):
    assert S.simplify(actual-S.sympify(expected)) == 0


def reference_args(args):return tuple(R.sympify(str(a)) for a in args)


@pytest.mark.parametrize('name,args',C.items())
def test_density_matches_reference(name,args):
    X = getattr(st,name)('density_'+name,*args)
    Y = getattr(ref,name)('Y',*reference_args(args))
    point = Q(1,2) if name in ('Beta','BetaNoncentral','Kumaraswamy','LogitNormal') else S.Integer(2) if name in st._univariate.DISCRETE else Q(3,2)
    ours,theirs = st.density(X)(point),ref.density(Y)(R.sympify(str(point)))
    if name == 'BetaNoncentral':
        # Compare the actual Poisson mixture at the first 35 components.
        j = next(v for v in ours.free_symbols if str(v).startswith('_j')) if ours.free_symbols else ours.args[1].args[0]
        ours = S.Add(*(ours.args[0].subs(j,k) for k in range(35)))
    assert float(ours.evalf()) == pytest.approx(float(theirs.evalf()),rel=1e-11,abs=1e-12)
    assert isinstance(ours,S.Expr)
    assert st.pspace(X).distribution is X.distribution


@pytest.mark.parametrize('name,args',C.items())
def test_distribution_support_and_cdf_endpoints(name,args):
    X = getattr(st,name)('support_'+name,*args);model = X.distribution
    assert model.cdf(-S.oo) == 0
    assert model.cdf(S.oo) == 1
    if model.lower != -S.oo:assert model.bounded_pdf(model.lower-1) == 0
    if model.upper != S.oo:assert model.bounded_pdf(model.upper+1) == 0
    assert st.E(X**0) == 1


@pytest.mark.parametrize('name,args',[(n,a) for n,a in C.items() if n not in ('BetaNoncentral','ChiNoncentral','Davis','Trapezoidal','VonMises','UniformSum') and n not in st._univariate.DISCRETE])
def test_cdf_derivative_is_density(name,args):
    model = getattr(st,name)('cdf_'+name,*args).distribution
    x = S.Symbol('x',real=True)
    point = Q(2,5) if name in ('Beta','Kumaraswamy','LogitNormal') else Q(7,5)
    value = model._cdf(x).diff(x).subs(x,point)
    assert float(value.evalf()) == pytest.approx(float(model.pdf(point).evalf()),rel=2e-9,abs=2e-10)


FINITE = {
    'DiscreteUniform':([1,2,5],), 'Die':(6,), 'Bernoulli':(Q(1,3),), 'Coin':(Q(1,3),),
    'Binomial':(5,Q(1,3)), 'BetaBinomial':(5,2,3), 'Hypergeometric':(10,4,3),
    'Rademacher':(), 'IdealSoliton':(6,), 'RobustSoliton':(10,Q(1,10),Q(1,10)),
}


@pytest.mark.parametrize('name,args',FINITE.items())
def test_finite_normalization_and_enumeration(name,args):
    X = getattr(st,name)('finite_'+name,*args);table = st.density(X).dict
    same(sum(table.values()),1)
    same(st.E(X),sum(value*weight for value,weight in table.items()))
    same(st.variance(X),sum(value**2*weight for value,weight in table.items())-st.E(X)**2)
    for value,weight in table.items():same(st.P(S.Eq(X,value)),weight)
    assert st.sample(X,size=(2,3),seed=0).shape == (2,3)


@pytest.mark.parametrize('name,args',[(n,a) for n,a in FINITE.items() if n not in ('Coin','RobustSoliton')])
def test_finite_tables_match_reference(name,args):
    X = getattr(st,name)('finite_ref_'+name,*args)
    converted = [list(a) if isinstance(a,list) else R.sympify(str(a)) for a in args]
    Y = getattr(ref,name)('Y',*converted)
    table = ref.density(Y).dict
    assert {str(k):float(v) for k,v in st.density(X).dict.items() if v != 0} == pytest.approx({str(k):float(v) for k,v in table.items() if v != 0})


MOMENTS = [
 ('Normal',(2,3),2,9),('Uniform',(1,4),Q(5,2),Q(3,4)),('Gamma',(2,3),6,18),
 ('Exponential',(3,),Q(1,3),Q(1,9)),('Beta',(2,3),Q(2,5),Q(1,25)),
 ('Binomial',(6,Q(1,3)),2,Q(4,3)),('Poisson',(3,),3,3),('Geometric',(Q(1,3),),3,6),
 ('NegativeBinomial',(3,Q(1,3)),6,18),('ChiSquared',(3,),3,6),
 ('Laplace',(2,3),2,18),('Logistic',(2,3),2,3*S.pi**2),
 ('StudentT',(5,),0,Q(5,3)),('Pareto',(1,5),Q(5,4),Q(5,48)),
 ('FDistribution',(5,8),Q(4,3),Q(88,45)),('GammaInverse',(5,2),Q(1,2),Q(1,12)),
 ('GaussianInverse',(2,3),2,Q(8,3)),('Gumbel',(2,1),1+2*S.EulerGamma,2*S.pi**2/3),
 ('Moyal',(1,2),1+2*(S.EulerGamma+S.log(2)),2*S.pi**2),
 ('UniformSum',(3,),Q(3,2),Q(1,4)),('Triangular',(0,2,1),1,Q(1,6)),
 ('Trapezoidal',(0,1,2,3),Q(3,2),Q(5,12)),('QuadraticU',(0,2),1,Q(3,5)),
 ('WignerSemicircle',(2,),0,1),('Arcsin',(0,2),1,Q(1,2)),
 ('FlorySchulz',(Q(1,3),),5,12),('Hermite',(2,3),8,14),('Skellam',(2,3),-1,5),
 ('YuleSimon',(5,),Q(5,4),Q(25,48)),('Lomax',(5,2),Q(1,2),Q(5,12)),
]


@pytest.mark.parametrize('name,args,mean,var',MOMENTS)
def test_exact_moments(name,args,mean,var):
    X = getattr(st,name)('moment_'+name,*args)
    same(st.E(X),mean);same(st.variance(X),var)
    same(st.moment(X,2),var+scalar_square(mean))


def scalar_square(value):return S.sympify(value)**2


def test_queries_two_dice_and_transforms():
    X,Y = st.Die('X',6),st.Die('Y',6)
    same(st.P(X+Y>10),Q(1,12));same(st.P(X>3,X+Y>10),1)
    same(st.E(X,X>3),5);same(st.E(X+Y),7)
    same(st.covariance(X,X+Y),Q(35,12));same(st.covariance(X,Y),0)
    assert st.independent(X,Y) and st.dependent(X,X+Y)
    assert st.random_symbols(X+Y) == {X,Y}
    assert st.where(X>3).set == S.FiniteSet(4,5,6)
    assert st.median(X) == S.FiniteSet(3,4)
    assert st.cdf(X)[S.Integer(3)] == Q(1,2)
    assert st.quantile(X)(Q(1,3)) == 2
    assert st.E(st.given(X**2,X>3)) == Q(77,3)
    assert st.density(X+Y)[S.Integer(7)] == Q(1,6)
    same(st.factorial_moment(X,3),35)
    same(st.H(X),S.log(6));same(st.entropy(st.Bernoulli('coin',Q(1,2)),base=2),1)
    same(st.std(X)**2,Q(35,12));same(st.correlation(X,X),1)
    same(st.coskewness(X,X,X),0)


def test_normal_queries_and_truncation():
    X = st.Normal('normal_query',0,1);Y = st.Normal('normal_query_y',0,2)
    same(st.P(X>0),Q(1,2));same(st.P(X+Y>0),Q(1,2))
    same(st.E(X,X>0),S.sqrt(2/S.pi))
    same(st.variance(X,X>0),1-2/S.pi)
    same(st.E(X,S.Eq(X,2)),2)
    same(st.E(st.given(2*X,X>0)),2*S.sqrt(2/S.pi))
    same(st.density(st.given(X,X>0))(1),2*st.density(X)(1))
    same(st.cmoment(X,4),3);same(st.smoment(X,4),3);same(st.skewness(X),0);same(st.kurtosis(X),3)
    same(st.median(X).args[0],0)
    same(st.moment_generating_function(X)(1),S.exp(Q(1,2)))
    same(st.characteristic_function(X)(1),S.exp(-Q(1,2)))
    same(st.entropy(X),(1+S.log(2*S.pi))/2)
    same(st.moment_generating_function(st.Poisson('pois_mgf',2))(1),S.exp(2*(S.E-1)))
    assert all(st.sample(X,X>0,size=20,seed=3)>0)


def test_handmade_laws():
    x = S.Symbol('x');k = S.Symbol('k',integer=True)
    X = st.ContinuousRV(x,2*x,S.Interval(0,1));D = st.DiscreteRV(k,2**(-k),S.Range(1,S.oo))
    F = st.FiniteRV('F',{0:Q(1,3),1:Q(2,3)})
    assert st.density(X)(Q(1,2)) == 1
    assert st.density(D)(2) == Q(1,4)
    same(st.E(F),Q(2,3))
    assert st.density(X).free_symbols == set()
    assert st.density(D).free_symbols == set()
    assert st.sample(X,seed=1) == pytest.approx(np.sqrt(np.random.default_rng(1).random()),abs=1e-8)
    assert st.FiniteDistributionHandmade({0:1}).pdf(0) == 1
    assert st.ContinuousDistributionHandmade(S.Lambda(x,2*x),S.Interval(0,1)).pdf(Q(1,2)) == 1
    assert st.DiscreteDistributionHandmade(S.Lambda(k,2**(-k)),S.Naturals).pdf(2) == Q(1,4)


@pytest.mark.parametrize('name,args',[(n,a) for n,a in C.items() if n not in ('BetaNoncentral','ChiNoncentral','Davis','Gompertz','ShiftedGompertz','VonMises')])
def test_numeric_sampling_shape_seed_and_support(name,args):
    X = getattr(st,name)('sampling_'+name,*args)
    first,second = st.sample(X,size=(2,3),seed=42),st.sample(X,size=(2,3),seed=42)
    assert first.shape == (2,3)
    np.testing.assert_array_equal(first,second)
    assert np.all(np.isfinite(first))
    assert np.all(first >= float(X.distribution.lower))
    assert np.all(first <= float(X.distribution.upper))
    if X.distribution.kind == 'discrete':assert np.all(first == first.astype(int))


def test_sampling_queries_iterators_and_empirical_density():
    X,Y = st.Die('sample_x',6),st.Die('sample_y',6)
    result = st.sample(X+Y,size=30,seed=1)
    assert np.all((result>=2)&(result<=12))
    assert len(list(st.sample_iter(X,numsamples=4,seed=3))) == 4
    table = st.sampling_density(X,numsamples=50,seed=4)
    same(sum(table.values()),1)
    same(st.P(X>3,numsamples=6000,seed=3),st.P(X>3,numsamples=6000,seed=3))
    assert abs(float(st.E(X,numsamples=6000,seed=3))-3.5)<0.08
    kde = st.sampling_density(st.Normal('kde',0,1),numsamples=30,seed=1)
    assert float(kde(0).evalf())>0


def test_joint_gaussian_and_dependence():
    X = st.MultivariateNormal('joint_normal',[1,2],[[2,1],[1,3]])
    same(st.E(X[0]),1);same(st.variance(X[1]),3);same(st.covariance(X[0],X[1]),1)
    same(st.E((X[0]-1)**2*(X[1]-2)**2),8)
    same(st.P(X[0]+X[1]>3),Q(1,2))
    assert st.dependent(X[0],X[1])
    Y = st.MultivariateNormal('independent_normal',[0,0],S.eye(2))
    assert st.independent(Y[0],Y[1])
    assert st.Normal('normal_vector',[0,0],S.eye(2)).dimension == 2
    same(st.marginal_distribution(X,0)(0),st.density(X[0])(0))
    np.testing.assert_allclose(st.sample(X,size=30000,seed=2).mean(axis=0),[1,2],atol=0.03)
    same(st.density(Y)(0,0),1/(2*S.pi))


def test_dirichlet_multinomial_and_ewens():
    D = st.Dirichlet('dirichlet',2,3,4)
    same(st.E(D[0]),Q(2,9));same(st.covariance(D[0],D[1]),-Q(1,135))
    assert st.MultivariateBeta('mb',2,3).dimension == 2
    np.testing.assert_allclose(st.sample(D,size=15,seed=1).sum(axis=1),np.ones(15))
    M = st.Multinomial('multinomial',4,Q(1,3),Q(2,3))
    same(st.E(M[0]),Q(4,3));same(st.covariance(M[0],M[1]),-Q(8,9))
    same(st.P(S.Eq(M[0]+M[1],4)),1)
    N = st.NegativeMultinomial('negative_multinomial',3,Q(1,5),Q(3,10))
    same(st.E(N[0]),Q(6,5));same(st.covariance(N[0],N[1]),Q(18,25))
    E = st.MultivariateEwens('ewens',4,2)
    same(st.P(S.Eq(sum((i+1)*x for i,x in enumerate(E)),4)),1)
    np.testing.assert_array_equal(st.sample(E,size=20,seed=2)@np.arange(1,5),np.full(20,4))


def test_other_joint_laws_and_marginals():
    T = st.MultivariateT('joint_t',[0,0],[[1,0],[0,2]],5)
    same(st.variance(T[0]),Q(5,3));same(st.E(T[0]**2*T[1]**2),Q(50,3))
    assert st.sample(T,size=3,seed=4).shape == (3,2)
    L = st.MultivariateLaplace('joint_laplace',[1,2],[[2,1],[1,3]])
    same(st.E(L[0]),1);same(st.variance(L[0]),3);same(st.covariance(L[0],L[1]),3)
    NG = st.NormalGamma('normal_gamma',0,1,3,2)
    same(st.variance(NG[0]),1);same(st.E(NG[1]),Q(3,2));same(st.covariance(NG[0]**2,NG[1]),-Q(1,2))
    G = st.GeneralizedMultivariateLogGamma('log_gamma',1,2,[1,2],[1,1])
    O = st.GeneralizedMultivariateLogGammaOmega('log_gamma_o',S.eye(2),2,[1,2],[1,1])
    same(st.E(G[0]),st.E(O[0]));same(st.E(G[0]),1-S.EulerGamma)
    assert st.sample(G,size=3,seed=2).shape == (3,2)
    x,y = S.symbols('x y',real=True)
    H = st.JointRV('hand_joint',S.Lambda((x,y),S.exp(-(x*x+y*y)/2)/(2*S.pi)))
    same(st.density(H)(0,0),1/(2*S.pi))
    assert 'Integral' in str(st.marginal_distribution(H,0)(x))


def test_matrix_laws_and_covariances():
    W = st.Wishart('wishart',5,S.eye(2))
    assert st.E(W) == 5*S.eye(2)
    same(st.variance(W[0,0]),10);same(st.variance(W[0,1]),5)
    same(st.covariance(W[0,0],W[1,1]),0)
    M = st.MatrixGamma('matrix_gamma',Q(5,2),2,S.eye(2))
    assert st.E(M) == st.E(W)
    same(st.density(M)(S.eye(2)),st.density(W)(S.eye(2)))
    N = st.MatrixNormal('matrix_normal',S.zeros(2),[[2,1],[1,3]],[[3,1],[1,2]])
    same(st.covariance(N[0,0],N[1,1]),1)
    same(st.E(N[0,0]**2*N[1,1]**2),38)
    assert st.ExpectationMatrix(N).doit() == S.zeros(2)
    assert st.VarianceMatrix(N).doit() == 5*S.Matrix([[3,1],[1,2]])
    assert st.CrossCovarianceMatrix(N,N).doit() == st.VarianceMatrix(N).doit()
    T = st.MatrixStudentT('matrix_t',5,S.zeros(2),S.eye(2),S.eye(2))
    same(st.variance(T[0,0]),Q(1,3))
    for X in (W,M,N,T):
        draws = st.sample(X,size=(2,3),seed=1)
        assert draws.shape == (2,3,2,2)
        assert np.all(np.isfinite(draws))
    draws = st.sample(W,size=20000,seed=3)
    np.testing.assert_allclose(draws.mean(axis=0),5*np.eye(2),atol=0.08)
    assert np.all(np.linalg.eigvalsh(draws)>0)
    V = st.MatrixNormal('column_normal',[[0],[0]],[[2,1],[1,3]],[[2]])
    assert st.VarianceMatrix(V).shape == (2,2)
    assert st.VarianceMatrix(V).doit() == S.Matrix([[4,2],[2,6]])
    A,B = st.Die('column_a',6),st.Die('column_b',6)
    left,right = S.Matrix([A,B]),S.Matrix([A,B,A+B])
    assert st.VarianceMatrix(left).doit() == Q(35,12)*S.eye(2)
    cross = st.CrossCovarianceMatrix(left,right)
    assert cross.shape == (2,3)
    assert cross.doit() == Q(35,12)*S.Matrix([[1,0,1],[0,1,1]])


def test_symbolic_operators():
    X,Y = st.Die('held_x',6),st.Die('held_y',6)
    assert isinstance(st.P(X>3,evaluate=False),st.Probability)
    same(st.Probability(X>3).doit(),Q(1,2))
    same(st.Expectation(X).doit(),Q(7,2))
    same(st.Variance(X).doit(),Q(35,12))
    same(st.Covariance(X,X+Y).doit(),Q(35,12))
    same(st.Moment(X,2).doit(),Q(91,6))
    same(st.CentralMoment(X,2).doit(),Q(35,12))
    same((2*st.Expectation(X)+st.Expectation(Y)).doit(),Q(21,2))
    same(st.Variance(X).rewrite(st.Expectation).doit(),Q(35,12))
    same(st.Covariance(X,Y).rewrite(st.Expectation).doit(),0)
    same(st.Expectation(X+Y).expand().doit(),7)
    a = S.Symbol('a');same(st.Expectation(a*X).subs(a,2).doit(),7)


def test_markov_chains_contexts_and_stationary_law():
    T = S.Matrix([[Q(3,4),Q(1,4)],[Q(1,2),Q(1,2)]])
    C = st.DiscreteMarkovChain('chain',[0,1],T)
    same(st.P(S.Eq(C[2],1),S.Eq(C[0],0)),Q(5,16))
    same(st.E(C[2],S.Eq(C[0],0)),Q(5,16))
    assert C.stationary_distribution() == S.Matrix([[Q(2,3),Q(1,3)]])
    override = st.TransitionMatrixOf(C,S.Matrix([[0,1],[1,0]]))
    same(st.P(S.Eq(C[1],1),S.And(S.Eq(C[0],0),override)),1)
    same(st.P(S.Eq(C[1],1),S.Eq(C[0],0)),Q(1,4))
    states = st.StochasticStateSpaceOf(C,[0,1]);assert states is not None
    D = st.ContinuousMarkovChain('ct_chain',[0,1],[[-1,1],[2,-2]])
    same(st.P(S.Eq(D(1),1),S.Eq(D(0),0)),(1-S.exp(-3))/3)
    override_q = st.GeneratorMatrixOf(D,S.Matrix([[-2,2],[1,-1]]))
    same(st.P(S.Eq(D(1),1),S.And(S.Eq(D(0),0),override_q)),2*(1-S.exp(-3))/3)
    assert st.sample(C[2],S.Eq(C[0],0),size=15,seed=1).shape == (15,)
    assert next(st.sample_stochastic_process(C)) in C.state_space
    assert st.StochasticProcess('abstract').state_space == S.Reals
    assert st.DiscreteTimeStochasticProcess('discrete')[2].time == 2


def test_increment_processes_and_bridges():
    B = st.BernoulliProcess('bernoulli',Q(1,3))
    same(st.P(S.Eq(B[1]+B[2],2)),Q(1,9));same(st.covariance(B[1],B[2]),0)
    assert st.independent(B[1],B[2])
    W = st.WienerProcess('wiener')
    same(st.covariance(W(1),W(3)),1);same(st.E(W(1)**2*W(3)**2),5)
    same(st.E(W(1),S.Eq(W(2),2)),1);same(st.variance(W(1),S.Eq(W(2),2)),Q(1,2))
    same(st.P(W(1)>1,S.Eq(W(2),2)),Q(1,2))
    same(st.variance(W(3)-W(1)),2)
    P = st.PoissonProcess('poisson_process',2)
    same(st.E(P(3),S.Eq(P(1),4)),8);same(st.E(P(1),S.Eq(P(3),6)),2)
    same(st.covariance(P(1),P(3)),2)
    same(st.E(P(1)**2*P(2)),34)
    G = st.GammaProcess('gamma_process',2,3)
    same(st.E(G(3),S.Eq(G(1),4)),7);same(st.E(G(1),S.Eq(G(3),6)),2)
    same(st.covariance(G(1),G(3)),Q(3,4))
    for process in (W,P,G):
        assert st.E(process(0)) == 0
        values = st.sample(process(3)-process(1),size=30,seed=2)
        assert np.all(np.isfinite(values))
        if process is not W:assert np.all(values>=0)
        assert next(st.sample_stochastic_process(process)) == 0


@pytest.mark.parametrize('family',['GaussianOrthogonalEnsemble','GaussianUnitaryEnsemble','GaussianSymplecticEnsemble','CircularOrthogonalEnsemble','CircularUnitaryEnsemble','CircularSymplecticEnsemble'])
def test_ensembles_matrix_structure_and_spectral_normalization(family):
    model = getattr(st,family)('ensemble_'+family,2)
    draws = st.sample(model,size=3,seed=2)
    assert draws.shape == (3,)+model.shape
    for draw in draws:
        if family.startswith('Gaussian'):np.testing.assert_allclose(draw,draw.conj().T,atol=1e-12)
        else:np.testing.assert_allclose(draw@draw.conj().T,np.eye(model.rows),atol=1e-12)
    if family == 'GaussianSymplecticEnsemble':
        eigenvalues = np.linalg.eigvalsh(draws[0]);np.testing.assert_allclose(eigenvalues[::2],eigenvalues[1::2],atol=1e-12)
    spacing = st.level_spacing_distribution(model)
    f = R.lambdify([R.Symbol(str(v)) for v in spacing.variables],R.sympify(str(spacing.expr)),'mpmath')
    assert float(mp.quad(lambda x:f(x),[0,1,mp.inf])) == pytest.approx(1,abs=1e-9)
    assert float(mp.quad(lambda x:x*f(x),[0,1,mp.inf])) == pytest.approx(1,abs=1e-9)
    singleton = getattr(st,family)('single_'+family,1)
    pdf = st.joint_eigen_distribution(singleton)
    f = R.lambdify([R.Symbol(str(v)) for v in pdf.variables],R.sympify(str(pdf.expr)),'mpmath')
    domain = [0,2*mp.pi] if family.startswith('Circular') else [-mp.inf,0,mp.inf]
    assert float(mp.quad(lambda x:f(x),domain)) == pytest.approx(1,abs=1e-9)


def test_abstract_ensembles_and_eigenvalue_pushforward():
    for factory in (st.GaussianEnsemble,st.CircularEnsemble):
        G = factory('abstract_ensemble_'+factory.__name__,2)
        assert 'beta' in str(st.joint_eigen_distribution(G))
    X = st.Normal('eigen_x',0,1)
    D = st.JointEigenDistribution(S.Matrix([[X,0],[0,2*X]]))
    assert D.dimension == 2
    assert D.sample(np.random.default_rng(1)).shape == (2,)
    assert 'Integral' in str(D.pdf(0,0))
    with pytest.raises(ValueError):st.joint_eigen_distribution(S.eye(2))


@pytest.mark.parametrize('name,args',[
 ('Normal',(0,0)),('Beta',(-1,2)),('Gamma',(0,1)),('Uniform',(2,1)),('Exponential',(-1,)),
 ('Binomial',(-1,Q(1,2))),('Die',(0,)),('Bernoulli',(2,)),('Geometric',(0,)),
 ('Hypergeometric',(3,4,2)),('Poisson',(-1,)),('Dirichlet',(0,1)),
 ('MultivariateNormal',([0,0],[[1,2],[2,1]])),('Multinomial',(2,Q(1,2),Q(1,3))),
 ('DiscreteMarkovChain',([0,1],[[1,1],[0,1]])),('Wishart',(0,S.eye(2))),
 ('MatrixNormal',(S.zeros(2),S.eye(3),S.eye(2))),
])
def test_invalid_parameters_raise(name,args):
    with pytest.raises((ValueError,TypeError)):getattr(st,name)('invalid_'+name,*args)


def test_export_inventory_and_no_reference_runtime():
    assert set(st.__all__) == set(ref.__all__)
    for name in st.__all__:assert callable(getattr(st,name))
    code = """
import sys
class Block:
 def find_spec(self,fullname,path=None,target=None):
  if fullname == 'sympy' or fullname.startswith('sympy.'):
   raise AssertionError('runtime SymPy import: '+fullname)
sys.meta_path.insert(0,Block())
import symox as S
from symox.stats import *
X=Normal('X',0,1)
assert E(X)==0 and variance(X)==1 and P(X>0)==S.Rational(1,2)
assert E(Die('D',6))==S.Rational(7,2)
assert E(MultivariateNormal('M',[0,0],S.eye(2))[0])==0
assert covariance(WienerProcess('W')(1),WienerProcess('V')(2))==0
"""
    result = subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
    assert result.returncode == 0,result.stderr


def test_conditioning_preserves_joint_dependence_and_point_samples():
    X = st.MultivariateNormal('conditional_joint',[0,0],[[1,Q(1,2)],[Q(1,2),1]])
    condition = S.Eq(X[1],2)
    same(st.E(X[0],condition),1);same(st.variance(X[0],condition),Q(3,4))
    same(st.P(X[0]>1,condition),Q(1,2))
    same(st.E(st.given(X[0],condition)),1)
    same(st.density(X[0],condition)(1),2/S.sqrt(6*S.pi))
    draws = st.sample(X[0],condition,size=20000,seed=4)
    assert draws.mean() == pytest.approx(1,abs=0.025)
    assert draws.var() == pytest.approx(0.75,abs=0.03)
    NG = st.NormalGamma('conditional_normal_gamma',1,2,3,4)
    same(st.variance(NG[0],S.Eq(NG[1],2)),Q(1,4))
    same(st.E(NG[1],S.Eq(NG[0],2)),Q(7,10))
    Y = st.Normal('point_normal',0,1)
    np.testing.assert_array_equal(st.sample(Y,S.Eq(Y,3),size=4,seed=1),np.full(4,3))
    with pytest.raises(ValueError):st.sample(st.Die('impossible_die',6),S.Eq(st.Die('impossible_die',6),9))
    assert isinstance(st.variance(Y,evaluate=False),st.Variance)
    assert isinstance(st.covariance(Y,Y,evaluate=False),st.Covariance)


@pytest.mark.parametrize('module,names',[
 ('crv_types',list(st._univariate.CONTINUOUS)+['ContinuousRV']),
 ('frv_types',list(st._univariate.FINITE)+['FiniteRV']),
 ('drv_types',list(st._univariate.DISCRETE)+['DiscreteRV']),
 ('joint_rv_types',st._joint.__all__),('matrix_distributions',st._matrix.__all__),
 ('stochastic_process_types',st._processes.__all__),('random_matrix_models',st._ensembles.__all__),
 ('symbolic_probability',st._symbolic.__all__),('rv_interface',st._queries.__all__),
])
def test_compatible_submodule_imports(module,names):
    imported = importlib.import_module('symox.stats.'+module)
    for name in names:assert getattr(imported,name) is getattr(st,name)


def test_probability_space_helpers():
    from symox.stats import rv
    X,Y = st.Die('space_x',6),st.Die('space_y',6)
    assert rv.is_random(X+Y)
    assert rv.pspace_independent(X,Y)
    assert rv.rv_subs(X+Y) == X.symbol+Y.symbol
    Z = st.Die('space_x',4)
    assert rv.rs_swap([X],[Z]) == {X:Z}
    assert rv.sumsets([S.FiniteSet(1),S.FiniteSet(2)]) == frozenset([S.Integer(1),S.Integer(2)])
    assert rv.expectation(X) == Q(7,2)


def test_sampling_a_probability_space_preserves_component_dependence():
    X = st.MultivariateNormal('space_dependent',[0,0],[[1,Q(3,4)],[Q(3,4),1]])
    values = st.pspace(X[0]+X[1]).sample(size=20000,seed=1)
    assert np.cov(values[X[0]],values[X[1]])[0,1] == pytest.approx(0.75,abs=0.03)
    Y = st.Binomial('scaled_binomial',2,Q(1,2),2,0)
    same(st.density(Y)(2),Q(1,2));same(st.P(S.Eq(Y,2)),Q(1,2))
    assert st.E(st.NormalGamma('heavy_ng',0,1,Q(1,4),2)[0]) == S.nan
    n = S.Symbol('n',integer=True,positive=True)
    D = st.Die('symbolic_die',n)
    same(st.E(D**3),n*(n+1)**2/4)
    assert D is st.Die('symbolic_die',n)
    same(st.density(st.ChiNoncentral('central_chi',3,0))(1),st.density(st.Chi('central_chi_ref',3))(1))
    same(st.density(st.BetaNoncentral('central_beta',2,3,0))(Q(1,2)),Q(3,2))
    Z = st.Skellam('zero_skellam',0,0)
    assert st.density(Z)(0) == 1 and st.density(Z)(-1) == 0
    U = st.Uniform('conditioned_uniform',0,1)
    with pytest.raises(ValueError):st.E(U,S.Eq(U,2))


@pytest.mark.parametrize('beta,family',[(1,'Orthogonal'),(2,'Unitary'),(4,'Symplectic')])
def test_two_eigenvalue_joint_densities_are_normalized(beta,family):
    G = getattr(st,'Gaussian'+family+'Ensemble')('gaussian_norm_'+family,2)
    distribution = st.joint_eigen_distribution(G)
    variables = [R.Symbol(str(v)) for v in distribution.variables]
    f = R.lambdify(variables,R.sympify(str(distribution.expr)),'mpmath')
    gaussian_u_mass = mp.sqrt(2*mp.pi/beta)
    mass = gaussian_u_mass*mp.quad(lambda v:f(v/mp.sqrt(2),-v/mp.sqrt(2)),[-mp.inf,0,mp.inf])
    assert float(mass) == pytest.approx(1,abs=1e-10)
    C = getattr(st,'Circular'+family+'Ensemble')('circular_norm_'+family,2)
    distribution = st.joint_eigen_distribution(C)
    variables = [R.Symbol(str(v)) for v in distribution.variables]
    f = R.lambdify(variables,R.sympify(str(distribution.expr)),'mpmath')
    mass = 2*mp.pi*mp.quad(lambda theta:mp.re(f(theta,0)),[0,mp.pi,2*mp.pi])
    assert float(mass) == pytest.approx(1,abs=1e-10)


def test_matrix_scalar_reductions_and_ensemble_expectations():
    x = S.Matrix([[Q(3,2)]])
    W = st.Wishart('scalar_wishart',5,[[2]])
    same(st.density(W)(x),st.density(st.Gamma('scalar_gamma',Q(5,2),4))(Q(3,2)))
    N = st.MatrixNormal('scalar_mn',[[1]],[[2]],[[3]])
    same(st.density(N)(x),st.density(st.Normal('scalar_normal',1,S.sqrt(6)))(Q(3,2)))
    T = st.MatrixStudentT('scalar_mt',5,[[1]],[[2]],[[3]])
    scale = S.sqrt(Q(6,5))
    assert float(st.density(T)(x).evalf()) == pytest.approx(float((st.density(st.StudentT('scalar_t',5))((Q(3,2)-1)/scale)/scale).evalf()),rel=1e-12)
    for cls in [st.GaussianOrthogonalEnsemble,st.GaussianUnitaryEnsemble,st.GaussianSymplecticEnsemble,st.CircularUnitaryEnsemble]:
        A = cls('ensemble_mean_'+cls.__name__,2)
        assert st.E(A) == S.zeros(*A.shape)
