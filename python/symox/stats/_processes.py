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


class JointDistributionHandmade:
    """SymPy ``JointDistributionHandmade``: a joint density given as a
    ``Lambda`` over the process values."""
    def __init__(self,pdf):self.args = (pdf,)
    @property
    def pdf(self):return self.args[0]
    def __call__(self,*args):return self.args[0](*args)
    @property
    def free_symbols(self):return set()
    def __repr__(self):return 'JointDistributionHandmade(%s)' % (self.args[0],)
    __str__ = __repr__
    def __eq__(self,other):return isinstance(other,JointDistributionHandmade) and self.args == other.args
    def __hash__(self):return hash(self.args)


class StochasticProcess:
    _stats_process = True
    _value_format = '%s(%s)'
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
            rv = variable(self._value_format % (self.name,time),self._model(time),group=self)
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
    _value_format = '%s[%s]'
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
        self._symbolic = self._transition is not None and not isinstance(self._transition,S.DenseMatrix)
        if self._transition is not None:check(self._transition.rows == self._transition.cols,'transition matrix must be square')
        if self._symbolic and not self._transition.rows.is_Integer:
            # A transition MatrixSymbol of symbolic size: the states are
            # Range(0, n) and only matrix-power queries are available.
            if state_space is not None:raise NotImplementedError('named states of a symbolic-size chain are not supported')
            self._states = ();self.state_space = S.Range(0,self._transition.rows);return
        # SymPy (MarkovProcess._sanity_checks): a chain without transition
        # data has the transition MatrixSymbol ``_T`` (of size a Dummy ``n``
        # when the states are not given either). Queries keep using the
        # transition data of their conditions.
        self._default_transition = None
        if self._transition is None and state_space is None:
            size = S.Dummy('n',integer=True,nonnegative=True)
            self._default_transition = S.MatrixSymbol('_T',size,size)
            self._states = ();self.state_space = S.Range(0,size);return
        if state_space is None:
            state_space = list(range(self._transition.rows)) if self._transition is not None else []
        self._states = tuple(scalar(x) for x in (state_space.args if isinstance(state_space,S.FiniteSet) else state_space))
        if len(set(self._states)) != len(self._states):raise ValueError('states must be distinct')
        self.state_space = S.FiniteSet(*self._states)
        if self._transition is not None:self._validate(self._transition,self._states)
        else:self._default_transition = S.MatrixSymbol('_T',len(self._states),len(self._states))
    def _validate(self,T,states):
        check(T.rows == T.cols == len(states),'matrix and state-space dimensions differ')
        if isinstance(T,S.DenseMatrix):
            for i in range(T.rows):
                for j in range(T.cols):
                    if not self.continuous or i != j:check(T[i,j] >= 0,'transition rates/probabilities must be nonnegative')
                gap = sum(T[i,j] for j in range(T.cols))-(0 if self.continuous else 1)
                check(gap == 0 or S.simplify(gap) == 0,'transition rows must sum to zero (generator) or one')
    @property
    def transition_probabilities(self):
        if self.continuous:return None
        return self._transition if self._transition is not None else getattr(self,'_default_transition',None)
    @property
    def generator_matrix(self):
        if not self.continuous:return None
        return self._transition if self._transition is not None else getattr(self,'_default_transition',None)
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
        if T is not None and not isinstance(T,S.DenseMatrix) and T.rows.is_Integer:
            # An explicit-size MatrixSymbol is used entrywise: T[i, j].
            if not states:states = tuple(S.Integer(i) for i in range(T.rows))
            check(T.rows == len(states),'matrix and state-space dimensions differ')
            T = T.as_explicit();self._stochastic_symbol = True
        elif T is not None:self._validate(T,states);self._stochastic_symbol = False
        if T is None or not states:raise ValueError('a transition/generator matrix and a finite state space are required')
        return states,T,S.And(*ordinary) if ordinary else None
    def _transition_at(self,T,elapsed):
        if elapsed.is_Rational or not isinstance(T,S.DenseMatrix):
            return (elapsed*T).exp() if self.continuous else T**elapsed
        if not self.continuous:
            # SymPy: trans_probs**(rv.key - given.key)
            try:return T**elapsed
            except (NotImplementedError,S.MatrixError):pass
        return self._symbolic_kernel(T,elapsed)
    def _symbolic_kernel(self,T,elapsed):
        # The reference transition function at a float or symbolic time:
        # Q*exp(t*D)*Q**-1 (CTMC) or Q*D**t*Q**-1 (DTMC) when T is
        # diagonalizable, else exp(t*G); evaluated at the elapsed time.
        cache = self.__dict__.setdefault('_kernels',{})
        key = (S.ImmutableMatrix(T),self.continuous)
        if key not in cache:
            t = S.Dummy('t');n = T.rows
            try:diagonalizable = T.is_diagonalizable()
            except Exception:diagonalizable = False
            if diagonalizable:
                Q,D = T.diagonalize()
                powers = [S.exp(t*D[i,i]) if self.continuous else D[i,i]**t for i in range(n)]
                kernel = Q*S.diag(*powers)*Q.inv()
            elif self.continuous:kernel = (t*T).exp()
            else:raise NotImplementedError('symbolic powers of a non-diagonalizable transition matrix are not implemented')
            cache[key] = (t,kernel)
        t,kernel = cache[key]
        return kernel.subs(t,elapsed)
    def _symbolic_time_query(self,operation,expr,states,T,gate):
        # The reference _symbolic_probability: one Relational on one random
        # indexed symbol, given one Eq at another (possibly symbolic) time.
        from ._queries import random_symbols
        from ._base import REGISTRY
        from ._symbolic import Probability
        relations = {'Equality','Unequality','StrictGreaterThan','GreaterThan','StrictLessThan','LessThan'}
        name = getattr(expr,'class_name','')
        targets = random_symbols(expr);given = random_symbols(gate)
        if operation != 'P' or name not in relations or len(targets) != 1 or len(given) != 1 or getattr(gate,'class_name','') != 'Equality':
            raise NotImplementedError('Markov queries at symbolic times support only P(Rel(X(a), b), Eq(X(c), d))')
        rv, = targets;start, = given
        curr = gate.rhs if gate.lhs == start._atom else gate.lhs
        if expr.lhs == rv._atom:nxt = expr.rhs
        else:
            nxt = expr.lhs
            name = {'StrictGreaterThan':'StrictLessThan','StrictLessThan':'StrictGreaterThan','GreaterThan':'LessThan','LessThan':'GreaterThan'}.get(name,name)
        kernel = self._transition_at(T,rv.time-start.time)
        index = range(len(states))
        def element(j):
            # Row ``curr`` of the kernel at column j (the state is a symbol).
            return S.Piecewise(*[(kernel[i,j],S.Eq(curr,states[i])) for i in index])
        if name in ('Equality','Unequality'):
            prob = S.Piecewise(*[(element(j),S.Eq(nxt,states[j])) for j in index],(S.Integer(0),True))
            if name == 'Unequality':prob = 1-prob
        else:
            upper = 0 if name in ('GreaterThan','StrictLessThan') else 1
            total = S.Add(*[S.Piecewise((element(j),S.Ge(states[j],nxt+upper)),(S.Integer(0),True)) for j in index])
            prob = total if name in ('GreaterThan','StrictGreaterThan') else 1-total
        return S.Piecewise((prob,rv.time > start.time),(Probability(expr),True))
    def _symbolic_size_query(self,operation,expr,condition):
        from ._queries import _bindings
        bindings = _bindings(condition)
        target = _bindings(expr) if operation == 'P' else None
        if not bindings or not target or len(bindings) != 1 or len(target) != 1 or self.continuous:
            raise NotImplementedError('a chain of symbolic size supports only P(Eq(Y[a], b), Eq(Y[c], d))')
        (start,i), = bindings.items();(stop,j), = target.items()
        from ._base import REGISTRY
        elapsed = REGISTRY[stop].time-REGISTRY[start].time
        check(elapsed >= 0,'the conditioning time must not follow the queried time')
        return (self._transition**elapsed)[i,j]
    def query(self,operation,expr,condition):
        from ._queries import random_symbols,_bindings,_indicator,_truth
        if self._symbolic and not self._transition.rows.is_Integer and not any(item in _CONTEXT for item in _parts(condition)):
            return self._symbolic_size_query(operation,expr,condition)
        states,T,gate = self._context(condition)
        variables = random_symbols(expr)|random_symbols(gate)
        if any(v.process is not self for v in variables):raise NotImplementedError('Markov queries require observations from one process')
        variables = sorted(variables,key=lambda v:float(v.time) if v.time.is_number else str(v.time))
        bindings = _bindings(gate)
        pinned = [(v,bindings[v._atom]) for v in variables if v._atom in bindings]
        if not pinned:
            from ._symbolic import Probability,Expectation
            return Probability(expr,condition) if operation == 'P' else Expectation(expr,condition)
        if any(not v.time.is_number for v in variables):return self._symbolic_time_query(operation,expr,states,T,gate)
        times = sorted({v.time for v in variables},key=lambda t:float(t) if t.is_number else str(t))
        earliest = times[0]
        anchors = [value for v,value in pinned if v.time == earliest]
        if not anchors:raise ValueError('conditioning must specify a state at the earliest observation time')
        initial = anchors[0];check(initial in states,'conditioned state is outside the state space')
        if len(times)>1 and len(states)**(len(times)-1)>1000000:raise RuntimeError('Markov path enumeration exceeds one million states')
        kernels = [self._transition_at(T,b-a) for a,b in zip(times,times[1:])]
        # Symbolic state labels are distinct by construction: events are
        # decided on their indices, while expectations keep the labels.
        codes = None if any(x.is_number for x in states) else {x:S.Integer(i) for i,x in enumerate(states)}
        event = scalar(expr) if codes is None else scalar(expr).xreplace(codes)
        coded_gate = gate if codes is None or gate is None else gate.xreplace(codes)
        numerator,denominator = S.Integer(0),S.Integer(0)
        for tail in product(states,repeat=len(times)-1):
            path = (initial,)+tail;weight = S.Integer(1)
            for kernel,a,b in zip(kernels,path,path[1:]):weight *= kernel[states.index(a),states.index(b)]
            mapping = {v._atom:path[times.index(v.time)] for v in variables}
            coded = mapping if codes is None else {k:codes[v] for k,v in mapping.items()}
            indicator = S.Integer(1) if gate is None else _indicator(coded_gate.subs(coded,simultaneous=True))
            denominator += weight*indicator
            if operation == 'P':numerator += weight*indicator*_indicator(event.subs(coded,simultaneous=True))
            else:
                if any(type(value).__name__ == 'Str' for value in mapping.values()):
                    # SymPy: a Str state is not an expression of a weighted sum.
                    raise TypeError("unsupported operand type(s) for *: 'Rational' and 'Str'")
                numerator += weight*indicator*scalar(expr).subs(mapping,simultaneous=True)
        if getattr(self,'_stochastic_symbol',False) and gate is not None:
            # Rows of a transition MatrixSymbol sum to one, so conditioning on
            # the initial state alone has probability one (as in SymPy).
            rest = [item for item in _parts(gate) if not (getattr(item,'class_name','') in ('Eq','Equality') and any(v._atom in item.free_symbols and v.time == earliest for v in variables) and len(random_symbols(item)) == 1)]
            if not rest:denominator = S.Integer(1)
        check(denominator != 0,'conditioning event has zero probability')
        return S.simplify(numerator/denominator)
    def stationary_distribution(self,condition_set=False):
        T = self._transition if self._transition is not None else getattr(self,'_default_transition',None)
        if T is None:raise ValueError('transition data are required')
        if not isinstance(T,S.DenseMatrix):
            # The reference symbolic version.
            n = T.rows
            if n == 0:return S.ImmutableMatrix(S.Matrix([[]]))
            wm = S.MatrixSymbol('wm',1,n)
            if condition_set:return S.ConditionSet(wm,S.Eq(wm*T,wm))
            return MatrixLambda((wm,T),S.Eq(wm*T,wm))
        n = len(self._states)
        if n == 0:return S.ImmutableMatrix(S.Matrix([[]]))
        if not self.continuous:
            # The reference linear system: (P - I)^T with its first row
            # replaced by ones, solved by linsolve (free parameters tau0...).
            a = S.Matrix(T-S.eye(n)).T
            a[0,0:n] = S.ones(1,n)
            b = S.zeros(n,1);b[0,0] = 1
            return S.ImmutableMatrix([list(list(S.linsolve((a,b)))[0])])
        symbols = S.symbols('_stationary0:'+str(n))
        equations = list(S.Matrix(1,n,symbols)*T)+[sum(symbols)-1]
        solutions = S.solve(equations,symbols,dict=True)
        if len(solutions) != 1:raise ValueError('stationary distribution is not unique')
        return S.ImmutableMatrix([[solutions[0].get(x,x) for x in symbols]])
    def limiting_distribution(self):
        """The solution of w Q = 0 with unit sum (a chain with a generator)."""
        if self._transition is None:return None
        if self._symbolic:raise NotImplementedError('a Lambda over matrix symbols is not supported')
        n = self._transition.rows;w = S.IndexedBase('w');wi = [w[i] for i in range(n)]
        equations = list(S.Matrix([wi])*self._transition)+[sum(wi)-1]
        return S.ImmutableMatrix([[S.simplify(v) for v in list(S.linsolve(equations,wi))[0]]])
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


class MatrixLambda:
    """``Lambda((wm, T), expr)`` over matrix symbols, which the native
    ``Lambda`` (scalar signatures) does not hold: the reference value of a
    symbolic ``stationary_distribution``. Calls substitute the arguments."""
    is_Lambda = True
    def __init__(self,signature,expr):
        self.signature = S.Tuple(*signature);self.expr = expr
    variables = property(lambda self:self.signature)
    args = property(lambda self:(self.signature,self.expr))
    nargs = property(lambda self:S.FiniteSet(len(self.signature)))
    free_symbols = property(lambda self:self.expr.free_symbols-set(self.signature))
    def __call__(self,*values):
        if len(values) != len(self.signature):
            raise TypeError('%s takes %d arguments (%d given)' % (self,len(self.signature),len(values)))
        return self.expr.xreplace(dict(zip(self.signature,values)))
    def __eq__(self,other):return isinstance(other,MatrixLambda) and self.args == other.args
    def __hash__(self):return hash(('Lambda',self.signature,self.expr))
    def __repr__(self):return 'Lambda(%s, %s)' % (self.signature,self.expr)
    __str__ = __repr__


class DiscreteMarkovChain(_Markov,DiscreteTimeStochasticProcess):
    def __init__(self,sym,state_space=None,trans_probs=None):self._initialize_chain(sym,state_space,trans_probs,False)
    def fixed_row_vector(self):return self.stationary_distribution()
    @property
    def limiting_distribution(self):return self.fixed_row_vector()
    @property
    def number_of_states(self):return S.Integer(len(self._states)) if self._states else self.state_space.size
    @property
    def index_of(self):return {state:S.Integer(i) for i,state in enumerate(self._states)}

    def communication_classes(self):
        """The reference ``communication_classes``: Tarjan's strongly
        connected components, their recurrence and their periods."""
        from ..utilities.iterables import strongly_connected_components
        T = self.transition_probabilities
        if T is None or not isinstance(T,S.DenseMatrix):
            raise NotImplementedError("Cannot perform the operation with a symbolic matrix.")
        n = len(self._states)
        V = list(range(n))
        E = [(i,j) for i in V for j in V if T[i,j] != 0]
        classes = strongly_connected_components((V,E))
        recurrence = [];periods = []
        for class_ in classes:
            submatrix = T.extract(class_,class_)
            is_recurrent = S.true
            for row in submatrix.tolist():
                if (sum(row)-1) != 0:
                    is_recurrent = S.false
                    break
            recurrence.append(is_recurrent)
            non_tree_edge_values = set()
            visited = {class_[0]};newly_visited = {class_[0]};level = {class_[0]:0}
            current_level = 0;done = False
            while not done:
                done = len(visited) == len(class_)
                current_level += 1
                for i in newly_visited:
                    newly_visited = {j for j in class_ if T[i,j] != 0}
                    new_tree_edges = newly_visited.difference(visited)
                    for j in new_tree_edges:level[j] = current_level
                    new_non_tree_edges = newly_visited.intersection(visited)
                    non_tree_edge_values |= {level[i]-level[j]+1 for j in new_non_tree_edges}
                    visited = visited.union(new_tree_edges)
            positive = {v for v in non_tree_edge_values if v > 0}
            if len(positive) == 0:periods.append(len(class_))
            elif len(positive) == 1:periods.append(positive.pop())
            else:periods.append(S.igcd(*positive))
        classes = [[scalar(self._states[i]) for i in class_] for class_ in classes]
        return list(zip(classes,recurrence,map(S.Integer,periods)))

    def decompose(self):
        T = self.transition_probabilities
        r_states = [];t_states = []
        for states,recurrent,period in self.communication_classes():
            if recurrent:r_states += states
            else:t_states += states
        states = r_states+t_states
        index = self.index_of
        indexes = [int(index[state]) for state in states]
        nr,nt = len(r_states),len(t_states)
        A = S.Matrix(nr,nr,lambda i,j:T[indexes[i],indexes[j]])
        B = S.Matrix(nt,nr,lambda i,j:T[indexes[nr+i],indexes[j]])
        C = S.Matrix(nt,nt,lambda i,j:T[indexes[nr+i],indexes[nr+j]])
        return states,A.as_immutable(),B.as_immutable(),C.as_immutable()

    def canonical_form(self):
        states,A,B,C = self.decompose()
        O = S.zeros(A.shape[0],C.shape[1])
        return states,S.BlockMatrix([[A,O],[B,C]]).as_explicit()


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
    def _pmf(self,x):
        return S.Piecewise((self.p,S.Eq(x,self.success)),(1-self.p,S.Eq(x,self.failure)),(0,True))
    def density(self,x):
        return self._pmf(scalar(x))
    def joint_distribution(self,*args):
        # SymPy StochasticProcess.joint_distribution: the product of the
        # marginal densities, as a Lambda over the process values.
        values = []
        for arg in args:
            if getattr(scalar(arg),'is_Number',False):arg = self[arg]
            elif getattr(arg,'process',None) is not self:
                raise ValueError('Expected a RandomIndexedSymbol or key not  %s' % (type(arg),))
            values.append(arg)
        names = tuple(S.Symbol(v.name) for v in values)
        return JointDistributionHandmade(S.Lambda(names,S.Mul(*(self._pmf(n) for n in names))))
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
        if self.kind != 'Poisson' and _has_contains(condition):
            return self._contains_query(operation,scalar(expr),scalar(condition))
        if self.kind == 'Poisson':
            result = self._poisson_query(operation,expr,condition)
            if result is not None:return result
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
    def _contains_query(self,operation,expr,condition):
        """The reference convention for conditions ``Contains(t, interval)``:
        the value at the symbolic time t is the increment over that interval,
        distributed as the process at time length (a Poisson count of mean
        lamda*length, a Normal(0, sqrt(length)), a Gamma(gamma*length,
        1/lamda)). Counts over pairwise
        disjoint intervals are independent; overlapping intervals give an
        unevaluated query."""
        from itertools import combinations
        from ._queries import random_symbols,P,E
        from ._univariate import Poisson
        args = condition.args if getattr(condition,'class_name','') in ('And','Or') else (condition,)
        if not all(_has_contains(arg) for arg in args):
            raise ValueError('If given condition is passed with `Contains`, then please pass the evaluated condition with its corresponding information in terms of intervals of each time stamp to be passed in given condition.')
        lengths,intervals = {},[]
        for v in sorted(random_symbols(expr),key=lambda v:str(v.time)):
            for arg in args:
                if arg.has(v.time) and v.time.is_Symbol:
                    interval = arg.args[1]
                    if isinstance(interval,(list,tuple,set)):interval = S.FiniteSet(*interval)
                    elif not hasattr(interval,'sup'):
                        # Boolean arguments hold structural copies of the
                        # set; rebuild the set object to reach its bounds.
                        interval = getattr(S,type(interval).__name__)(*interval.args)
                    length = interval.sup-interval.inf
                    if length == S.oo:raise ValueError('%s should have finite bounds' % str(v.name))
                    if length == 0:length = interval.sup
                    lengths[v] = length;intervals.append(interval)
        if len(intervals) > 1 and not all(S.Intersection(*pair) == S.EmptySet for pair in combinations(intervals,2)):
            from ._symbolic import Expectation,Probability
            return Probability(expr,condition) if operation == 'P' else Expectation(expr,condition)
        if any(v not in lengths for v in random_symbols(expr)):
            expr = expr.subs({v._atom:self(length)._atom for v,length in lengths.items()},simultaneous=True)
            return P(expr) if operation == 'P' else E(expr)
        ordered = sorted(lengths,key=lambda v:str(v.time))
        if operation == 'P' and self.kind == 'Poisson':
            totals = _PoissonEvents([v._atom for v in ordered],[self.lamda*lengths[v] for v in ordered],False).probabilities([expr])
            if totals is not None:return S.simplify(totals[0])
        fresh = {}
        for v in ordered:
            PoissonProcess._increments += 1
            name = '_%s_count%d' % (self.name,PoissonProcess._increments)
            fresh[v._atom] = scalar(Poisson(name,self.lamda*lengths[v]) if self.kind == 'Poisson' else variable(name,self._model(lengths[v])))
        expr = expr.subs(fresh,simultaneous=True)
        return P(expr) if operation == 'P' else E(expr)
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


_BOOLEAN_HEADS = ('And','Or','Not')
_RELATIONAL_HEADS = ('Equality','Unequality','StrictLessThan','LessThan','StrictGreaterThan','GreaterThan')


def _has_contains(expr):
    return expr is not None and any(getattr(node,'class_name','') == 'Contains' for node in S.preorder_traversal(expr))


def _event_constants(event,atoms):
    """The numeric constants c of an event built with And/Or/Not from
    relationals ``X op c`` (one process value against a real number), or
    None when the event has another shape."""
    head = getattr(event,'class_name','')
    if event in (S.true,S.false):return []
    if head in _BOOLEAN_HEADS:
        constants = []
        for arg in event.args:
            found = _event_constants(arg,atoms)
            if found is None:return None
            constants += found
        return constants
    if head not in _RELATIONAL_HEADS:return None
    lhs,rhs = event.args
    free = event.free_symbols
    if len(free) != 1 or not free <= atoms:return None
    if lhs in atoms and rhs.is_number and rhs.is_extended_real:return [rhs]
    if rhs in atoms and lhs.is_number and lhs.is_extended_real:return [lhs]
    return None


class _PoissonEvents:
    """Exact probabilities of events on the values of a Poisson process.

    The values are either the process at increasing numeric times (dependent,
    with independent increments) or independent Poisson variables. Every
    relational compares one value with a number, so all values at or above
    K = floor(max constant) + 1 behave alike and are lumped in one state;
    the probability is then a finite sum over the joint states."""
    _LIMIT = 20000

    def __init__(self,atoms,means,cumulative):
        self.atoms,self.means,self.cumulative = atoms,means,cumulative

    @staticmethod
    def _pmf(k,mean):
        if mean == 0:return S.Integer(1) if k == 0 else S.Integer(0)
        return S.exp(-mean)*mean**k/S.factorial(k)

    def probabilities(self,events):
        constants = []
        for event in events:
            found = _event_constants(event,set(self.atoms))
            if found is None:return None
            constants += found
        K = max([int(S.floor(c))+1 for c in constants]+[1])
        if K < 1 or (K+1)**len(self.atoms) > self._LIMIT:return None
        totals = [S.Integer(0) for _ in events]
        def visit(index,previous,weight,states):
            if weight == 0:return
            if index == len(self.atoms):
                rule = {a:(K if s is None else s) for a,s in zip(self.atoms,states)}
                for position,event in enumerate(events):
                    if S.sympify(event).subs(rule,simultaneous=True) == S.true:totals[position] += weight
                return
            mean = self.means[index]
            if self.cumulative and previous is None:
                visit(index+1,None,weight,states+[None]);return
            start = previous if self.cumulative else 0
            remaining = S.Integer(1)
            for value in range(start,K):
                p = self._pmf(value-start,mean);remaining -= p
                visit(index+1,value,weight*p,states+[value])
            visit(index+1,None,weight*remaining,states+[None])
        visit(0,0,S.Integer(1),[])
        return totals


class WienerProcess(_IncrementProcess):
    kind = 'Wiener'
    def __init__(self,sym):super().__init__(sym,S.Reals)


class PoissonProcess(_IncrementProcess):
    kind = 'Poisson'
    _increments = 0
    def __init__(self,sym,lamda):
        self.lamda = scalar(lamda);positive(self.lamda);super().__init__(sym,S.Naturals0)

    def __add__(self,other):
        if not isinstance(other,PoissonProcess):raise ValueError('Only instances of Poisson Process can be merged')
        return PoissonProcess(S.Dummy(self.symbol.name+other.symbol.name),self.lamda+other.lamda)

    def split(self,l1,l2):
        if S.sympify(l1+l2) != self.lamda:raise ValueError('Sum of l1 and l2 should be %s' % str(self.lamda))
        return (PoissonProcess(S.Dummy('l1'),l1),PoissonProcess(S.Dummy('l2'),l2))

    def _poisson_query(self,operation,expr,condition):
        from ._queries import random_symbols
        expr = scalar(expr)
        if condition is not None:condition = scalar(condition)
        if _has_contains(condition):return self._contains_query(operation,expr,condition)
        if operation != 'P':return None
        variables = random_symbols(expr)|random_symbols(condition)
        if not variables or any(getattr(v,'process',None) is not self for v in variables):return None
        if any(not v.time.is_number for v in variables):return None
        if condition is None and len({v.time for v in variables}) < 2:return None
        variables = sorted(variables,key=lambda v:float(v.time))
        times = [S.Integer(0)]+[v.time for v in variables]
        means = [self.lamda*(b-a) for a,b in zip(times,times[1:])]
        events = [expr] if condition is None else [S.And(expr,condition),condition]
        totals = _PoissonEvents([v._atom for v in variables],means,True).probabilities(events)
        if totals is None:return None
        if condition is None:return S.simplify(totals[0])
        if totals[1] == 0:return None
        return S.simplify(totals[0]/totals[1])

class GammaProcess(_IncrementProcess):
    kind = 'Gamma'
    def __init__(self,sym,lamda,gamma):
        self.lamda,self.gamma = scalar(lamda),scalar(gamma);positive(self.lamda,self.gamma);super().__init__(sym,S.Interval(0,S.oo))


def sample_stochastic_process(process):
    if not isinstance(process,StochasticProcess):raise ValueError('process must be a StochasticProcess')
    return process.sample()


__all__ = ['StochasticProcess','DiscreteTimeStochasticProcess','DiscreteMarkovChain','TransitionMatrixOf','StochasticStateSpaceOf','GeneratorMatrixOf','ContinuousMarkovChain','BernoulliProcess','PoissonProcess','WienerProcess','GammaProcess','sample_stochastic_process']
