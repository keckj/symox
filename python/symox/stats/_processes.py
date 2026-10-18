"""Finite-state Markov models and processes with independent increments."""
from itertools import product
import operator
import numpy as np
import symox as S
from ._base import Distribution,variable,scalar,check,positive,HALF
from ._joint import matrix,_normal_raw
from ._univariate import continuous_model,discrete_model,finite_model

_CONTEXT = {}


def _parts(condition):
    return condition.args if getattr(condition,'class_name','') == 'And' else (() if condition is None else (condition,))


class StochasticProcess:
    _stats_process = True
    def __init__(self,sym,state_space=S.Reals,**kwargs):
        if kwargs:raise TypeError('unknown stochastic-process options')
        self.symbol = S.Symbol(str(sym));self.name = str(sym);self.state_space = state_space
        self.index_set = S.Interval(0,S.oo);self._observations = {}
    def __repr__(self):return self.name
    def __call__(self,time):return self._at(time)
    def _model(self,time):
        x = S.Dummy('state',real=True)
        model = Distribution(time);model.family = 'AbstractProcess'
        model._pdf = lambda value:S.Function(self.name+'_density')(time,value)
        model._support = self.state_space
        return model
    def _at(self,time):
        time = scalar(time);check(time >= 0,'time must be nonnegative')
        if time not in self._observations:
            rv = variable('%s(%s)' % (self.name,time),self._model(time),group=self)
            rv.process,rv.time = self,time;self._observations[time] = rv
        return self._observations[time]
    def finite_worlds(self,variables):return None
    def query(self,operation,expr,condition):return None
    def sample(self,seed=None):
        rng = np.random.default_rng(seed)
        for time in range(1000000000000):
            observations = [self._at(time)]
            yield self.sample_mapping(rng,observations)[observations[0]._atom]


class DiscreteTimeStochasticProcess(StochasticProcess):
    def __init__(self,sym,state_space=S.Reals,**kwargs):
        super().__init__(sym,state_space,**kwargs);self.index_set = S.Range(0,S.oo)
    def __getitem__(self,time):return self._at(time)
    def _at(self,time):
        time = scalar(time);check(time.is_integer is not False,'discrete time must be an integer')
        return super()._at(time)


class _Markov:
    def _initialize_chain(self,sym,state_space,transition,continuous):
        StochasticProcess.__init__(self,sym)
        self.continuous = continuous
        self.index_set = S.Interval(0,S.oo) if continuous else S.Range(0,S.oo)
        self._transition = None if transition is None else matrix(transition)
        if self._transition is not None:check(self._transition.rows == self._transition.cols,'transition matrix must be square')
        if state_space is None:
            state_space = list(range(self._transition.rows)) if self._transition is not None else []
        self._states = tuple(scalar(x) for x in (state_space.args if isinstance(state_space,S.FiniteSet) else state_space))
        if len(set(self._states)) != len(self._states):raise ValueError('states must be distinct')
        self.state_space = S.FiniteSet(*self._states)
        if self._transition is not None:self._validate(self._transition,self._states)
    def _validate(self,T,states):
        check(T.rows == T.cols == len(states),'matrix and state-space dimensions differ')
        if isinstance(T,S.Matrix):
            for i in range(T.rows):
                for j in range(T.cols):
                    if not self.continuous or i != j:check(T[i,j] >= 0,'transition rates/probabilities must be nonnegative')
                check(S.simplify(sum(T[i,j] for j in range(T.cols))) == (0 if self.continuous else 1),'transition rows must sum to zero (generator) or one')
    @property
    def transition_probabilities(self):return self._transition if not self.continuous else None
    @property
    def generator_matrix(self):return self._transition if self.continuous else None
    def _model(self,time):
        model = Distribution(time);model.kind = 'finite';model.family = 'MarkovState';model._support = self.state_space
        # An initial law is deliberately left unspecified: transition data
        # determine conditional, rather than unconditional, probabilities.
        model._pdf = lambda x:S.Function(self.name+'_initial_probability')(time,x)
        return model
    def _context(self,condition):
        states,T = self._states,self._transition
        ordinary = []
        for item in _parts(condition):
            if item in _CONTEXT:
                process,kind,value = _CONTEXT[item]
                if process is not self:ordinary.append(item);continue
                if kind == 'states':states = value
                else:T = value
            else:ordinary.append(item)
        if T is None or not states:raise ValueError('a transition/generator matrix and a finite state space are required')
        self._validate(T,states)
        return states,T,S.And(*ordinary) if ordinary else None
    def _transition_at(self,T,elapsed):
        return (elapsed*T).exp() if self.continuous else T**elapsed
    def query(self,operation,expr,condition):
        from ._queries import random_symbols,_bindings,_indicator,_truth
        states,T,gate = self._context(condition)
        variables = random_symbols(expr)|random_symbols(gate)
        if any(v.process is not self for v in variables):raise NotImplementedError('Markov queries require observations from one process')
        variables = sorted(variables,key=lambda v:float(v.time) if v.time.is_number else str(v.time))
        bindings = _bindings(gate)
        pinned = [(v,bindings[v._atom]) for v in variables if v._atom in bindings]
        if not pinned:
            from ._symbolic import Probability,Expectation
            return Probability(expr,condition) if operation == 'P' else Expectation(expr,condition)
        times = sorted({v.time for v in variables},key=lambda t:float(t) if t.is_number else str(t))
        earliest = times[0]
        anchors = [value for v,value in pinned if v.time == earliest]
        if not anchors:raise ValueError('conditioning must specify a state at the earliest observation time')
        initial = anchors[0];check(initial in states,'conditioned state is outside the state space')
        if len(times)>1 and len(states)**(len(times)-1)>1000000:raise RuntimeError('Markov path enumeration exceeds one million states')
        kernels = [self._transition_at(T,b-a) for a,b in zip(times,times[1:])]
        numerator,denominator = S.Integer(0),S.Integer(0)
        for tail in product(states,repeat=len(times)-1):
            path = (initial,)+tail;weight = S.Integer(1)
            for kernel,a,b in zip(kernels,path,path[1:]):weight *= kernel[states.index(a),states.index(b)]
            mapping = {v._atom:path[times.index(v.time)] for v in variables}
            indicator = S.Integer(1) if gate is None else _indicator(gate.subs(mapping,simultaneous=True))
            denominator += weight*indicator
            outcome = scalar(expr).subs(mapping,simultaneous=True)
            numerator += weight*indicator*(_indicator(outcome) if operation == 'P' else outcome)
        check(denominator != 0,'conditioning event has zero probability')
        return S.simplify(numerator/denominator)
    def stationary_distribution(self,condition_set=False):
        if self._transition is None:raise ValueError('transition data are required')
        n = len(self._states);symbols = S.symbols('_stationary0:'+str(n))
        equations = list((S.Matrix(1,n,symbols)*self._transition-(S.zeros(1,n) if self.continuous else S.Matrix(1,n,symbols))))+[sum(symbols)-1]
        solutions = S.solve(equations,symbols,dict=True)
        if len(solutions) != 1:raise ValueError('stationary distribution is not unique')
        return S.ImmutableMatrix([[solutions[0].get(x,x) for x in symbols]])
    def sample_mapping(self,rng,variables,condition=None):
        from ._queries import _bindings
        states,T,_ = self._context(condition);bindings = _bindings(condition)
        times = sorted({v.time for v in variables},key=float)
        pinned = {v.time:bindings[v._atom] for v in variables if v._atom in bindings}
        initial = pinned.get(times[0],states[0])
        if any(t != times[0] for t in pinned):raise NotImplementedError('sampling a Markov bridge requires a single earliest state condition')
        path = {times[0]:initial};previous = times[0];state = initial
        for time in times[1:]:
            kernel = self._transition_at(T,time-previous);row = np.array([float(kernel[states.index(state),j].evalf()) for j in range(len(states))])
            row = np.maximum(row,0);row /= row.sum();state = states[int(rng.choice(len(states),p=row))]
            path[time] = state;previous = time
        return {v._atom:path[v.time] for v in variables}
    def sample(self,seed=None):
        rng = np.random.default_rng(seed);states,T,_ = self._context(None);state = states[0]
        kernel = self._transition_at(T,S.Integer(1))
        while True:
            yield state
            row = np.array([float(kernel[states.index(state),j].evalf()) for j in range(len(states))]);row = np.maximum(row,0);row /= row.sum()
            state = states[int(rng.choice(len(states),p=row))]


class DiscreteMarkovChain(_Markov,DiscreteTimeStochasticProcess):
    def __init__(self,sym,state_space=None,trans_probs=None):self._initialize_chain(sym,state_space,trans_probs,False)


class ContinuousMarkovChain(_Markov,StochasticProcess):
    def __init__(self,sym,state_space=None,gen_mat=None):self._initialize_chain(sym,state_space,gen_mat,True)


def _context_predicate(process,kind,value):
    if not isinstance(process,_Markov):raise TypeError('Markov context requires a Markov chain')
    if kind == 'states':
        value = tuple(scalar(x) for x in (value.args if isinstance(value,S.FiniteSet) else value));argument = S.Tuple(*value)
    else:value = matrix(value);argument = value
    predicate = S.Predicate({'transition':'TransitionMatrixOf','generator':'GeneratorMatrixOf','states':'StochasticStateSpaceOf'}[kind])(process.symbol,argument)
    _CONTEXT[predicate] = process,kind,value
    return predicate


def TransitionMatrixOf(process,matrix):return _context_predicate(process,'transition',matrix)
def GeneratorMatrixOf(process,matrix):return _context_predicate(process,'generator',matrix)
def StochasticStateSpaceOf(process,state_space):return _context_predicate(process,'states',state_space)


class BernoulliProcess(DiscreteTimeStochasticProcess):
    def __init__(self,sym,p,success=1,failure=0):
        self.p,self.success,self.failure = map(scalar,(p,success,failure))
        finite_model('Bernoulli',(self.p,self.success,self.failure))
        super().__init__(sym,S.FiniteSet(self.success,self.failure))
    def _model(self,time):return finite_model('Bernoulli',(self.p,self.success,self.failure))
    def finite_worlds(self,variables):
        choices = [list(v.distribution.dict.items()) for v in variables]
        for entries in product(*choices):
            weight = S.Integer(1);mapping = {}
            for v,(value,p) in zip(variables,entries):mapping[v._atom] = value;weight *= p
            yield mapping,weight
    def joint_moment(self,factors):return S.Mul(*(rv.distribution.moment(power) for rv,power in factors.items()))
    def independent_components(self,first,second):return first != second
    def sample_mapping(self,rng,variables,condition=None):return {v._atom:scalar(v.distribution.sample(rng,())) for v in variables}


class _IncrementProcess(StochasticProcess):
    def _model(self,time):
        if time == 0:return finite_model('Bernoulli',(S.Integer(0),S.Integer(1),S.Integer(0)))
        if self.kind == 'Wiener':return continuous_model('Normal',(S.Integer(0),S.sqrt(time)))
        if self.kind == 'Poisson':return discrete_model('Poisson',(self.lamda*time,))
        return continuous_model('Gamma',(self.gamma*time,1/self.lamda))
    def mean_at(self,t):return S.Integer(0) if self.kind == 'Wiener' else (self.lamda*t if self.kind == 'Poisson' else self.gamma*t/self.lamda)
    def covariance_at(self,s,t):return S.Min(s,t)*(1 if self.kind == 'Wiener' else self.lamda if self.kind == 'Poisson' else self.gamma/self.lamda**2)
    def normal_linear(self,expr,variables):
        if self.kind != 'Wiener':return None
        try:p = S.Poly(expr,*(v._atom for v in variables))
        except (ValueError,TypeError):return None
        if p.total_degree()>1:return None
        mean = expr.subs({v._atom:0 for v in variables},simultaneous=True)
        var = sum(expr.diff(a._atom)*expr.diff(b._atom)*S.Min(a.time,b.time) for a in variables for b in variables)
        return mean,S.simplify(var)
    def joint_moment(self,factors):
        variables = sorted(factors,key=lambda v:float(v.time) if v.time.is_number else str(v.time));powers = [factors[v] for v in variables]
        if self.kind == 'Wiener':return _normal_raw(powers,[S.Integer(0)]*len(variables),S.Matrix([[S.Min(a.time,b.time) for b in variables] for a in variables]))
        if sum(powers)==1:return self.mean_at(variables[0].time)
        if sum(powers)==2:
            a,b = (variables*2)[:2]
            return self.mean_at(a.time)*self.mean_at(b.time)+self.covariance_at(a.time,b.time)
        if any(not v.time.is_number for v in variables):raise NotImplementedError('higher process moments require ordered numeric observation times')
        increments = [S.Dummy('increment'+str(i)) for i in range(len(variables))]
        polynomial = S.Mul(*(sum(increments[:i+1])**factors[v] for i,v in enumerate(variables)))
        times = [S.Integer(0)]+[v.time for v in variables]
        models = [self._model(b-a) for a,b in zip(times,times[1:])]
        return S.simplify(sum(coefficient*S.Mul(*(model.moment(power) for model,power in zip(models,monomial))) for monomial,coefficient in S.Poly(polynomial,*increments).terms()))
    def _gaussian_condition(self,variables,bindings):
        from ._queries import random_symbols
        observed = [v for v in random_symbols(S.Add(*bindings)) if v.process is self]
        observed = [v for v in observed if v._atom in bindings]
        if not observed:return [S.Integer(0)]*len(variables),S.Matrix([[S.Min(a.time,b.time) for b in variables] for a in variables])
        observed = [v for v in observed if v.time != 0]
        if not observed:return [S.Integer(0)]*len(variables),S.Matrix([[S.Min(a.time,b.time) for b in variables] for a in variables])
        K = S.Matrix([[S.Min(a.time,b.time) for b in observed] for a in observed]);cross = S.Matrix([[S.Min(a.time,b.time) for b in observed] for a in variables])
        weights = cross*K.inv();mean = weights*S.Matrix([bindings[v._atom] for v in observed])
        covariance = S.Matrix([[S.Min(a.time,b.time) for b in variables] for a in variables])-weights*cross.T
        return [mean[i,0] for i in range(len(variables))],covariance
    def query(self,operation,expr,condition):
        if condition is None:return None
        from ._queries import random_symbols,_bindings,_truth,_indicator,P,E
        bindings = _bindings(condition)
        if not bindings or _truth(scalar(condition).subs(bindings,simultaneous=True)) is not True:return None
        all_variables = random_symbols(expr)|random_symbols(condition)
        if any(v.process is not self for v in all_variables):return None
        for v in all_variables:
            if v.time == 0 and v._atom in bindings:check(bindings[v._atom] == 0,'increment processes start at zero')
        expression = scalar(expr).subs(bindings,simultaneous=True)
        variables = sorted(random_symbols(expression),key=lambda v:float(v.time) if v.time.is_number else str(v.time))
        if not variables:return _indicator(expression) if operation == 'P' else expression
        if self.kind == 'Wiener':
            mean,cov = self._gaussian_condition(variables,bindings)
            if operation == 'E':
                try:poly = S.Poly(expression,*(v._atom for v in variables))
                except (ValueError,TypeError):raise NotImplementedError('conditional Wiener expectation requires a polynomial')
                return S.simplify(sum(coefficient*_normal_raw(powers,mean,cov) for powers,coefficient in poly.terms()))
            if expression.is_Relational:
                delta = expression.lhs-expression.rhs
                if S.Poly(delta,*(v._atom for v in variables)).total_degree()<=1:
                    coefficients = [delta.diff(v._atom) for v in variables]
                    m = delta.subs({v._atom:a for v,a in zip(variables,mean)},simultaneous=True)
                    var = S.simplify(sum(coefficients[i]*coefficients[j]*cov[i,j] for i in range(len(variables)) for j in range(len(variables))))
                    if var == 0:return _indicator(S.Rel(m,0,expression.rel_op))
                    from ._univariate import Normal
                    Z = Normal('_wiener_conditional',m,S.sqrt(var));return P(S.Rel(Z,0,expression.rel_op))
            raise NotImplementedError('conditional Wiener probabilities require a linear relational event')
        observed = [(v,bindings[v._atom]) for v in all_variables if v._atom in bindings]
        if len(observed)!=1 or len(variables)!=1:raise NotImplementedError('conditioned Poisson/Gamma queries require one observation and one target time')
        anchor,value = observed[0];target = variables[0];dt = S.simplify(target.time-anchor.time)
        if dt.is_nonnegative:
            model = self._model(dt);Z = variable('_conditioned_increment',model)
            replacement = value+Z
        elif dt.is_negative:
            if self.kind == 'Poisson':
                check(value.is_integer is not False and value.is_nonnegative is not False,'a Poisson state must be a nonnegative integer')
                model = finite_model('Binomial',(value,target.time/anchor.time,S.Integer(1),S.Integer(0)));replacement = variable('_poisson_bridge',model)
            else:
                model = continuous_model('Beta',(self.gamma*target.time,self.gamma*(anchor.time-target.time)));replacement = value*variable('_gamma_bridge',model)
        else:raise NotImplementedError('conditional increment query requires ordered times')
        expression = expression.subs(target._atom,scalar(replacement))
        return P(expression) if operation == 'P' else E(expression)
    def sample_mapping(self,rng,variables,condition=None):
        from ._queries import _bindings,_truth
        bindings = _bindings(condition)
        times = sorted({v.time for v in variables},key=float)
        if bindings and self.kind == 'Wiener':
            mean,cov = self._gaussian_condition(variables,bindings)
            values = rng.multivariate_normal(np.array([float(a) for a in mean]),np.array([[float(cov[i,j]) for j in range(len(variables))] for i in range(len(variables))]))
            return {v._atom:scalar(value) for v,value in zip(variables,values)}
        anchors = [(v.time,bindings[v._atom]) for v in variables if v._atom in bindings]
        if anchors:
            if len(anchors)!=1 or any(t<anchors[0][0] for t in times):raise NotImplementedError('conditioned path sampling currently supports a single initial state')
            previous,state = anchors[0]
        else:previous,state = S.Integer(0),S.Integer(0)
        path = {}
        for t in times:
            if t == previous:path[t] = state;continue
            model = self._model(t-previous);state += scalar(model.sample(rng,()));path[t] = state;previous = t
        return {v._atom:path[v.time] for v in variables}
    def sample(self,seed=None):
        rng = np.random.default_rng(seed);state = S.Integer(0);model = self._model(S.Integer(1))
        while True:
            yield state
            state += scalar(model.sample(rng,()))


class WienerProcess(_IncrementProcess):
    kind = 'Wiener'
    def __init__(self,sym):super().__init__(sym,S.Reals)


class PoissonProcess(_IncrementProcess):
    kind = 'Poisson'
    def __init__(self,sym,lamda):
        self.lamda = scalar(lamda);positive(self.lamda);super().__init__(sym,S.Range(0,S.oo))


class GammaProcess(_IncrementProcess):
    kind = 'Gamma'
    def __init__(self,sym,lamda,gamma):
        self.lamda,self.gamma = scalar(lamda),scalar(gamma);positive(self.lamda,self.gamma);super().__init__(sym,S.Interval(0,S.oo))


def sample_stochastic_process(process):
    if not isinstance(process,StochasticProcess):raise ValueError('process must be a StochasticProcess')
    return process.sample()


__all__ = ['StochasticProcess','DiscreteTimeStochasticProcess','DiscreteMarkovChain','TransitionMatrixOf','StochasticStateSpaceOf','GeneratorMatrixOf','ContinuousMarkovChain','BernoulliProcess','PoissonProcess','WienerProcess','GammaProcess','sample_stochastic_process']
