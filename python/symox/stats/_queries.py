"""Exact probability, expectation, transforms and conditional sampling."""
from collections import defaultdict
from itertools import product
import operator
import numpy as np
import symox as S
from ._base import REGISTRY,RandomSymbol,Distribution,Domain,PSpace,HALF,scalar,check,variable,integral,sum_expr


def random_symbols(expr):
    if expr is None:return set()
    if isinstance(expr,RandomSymbol):return {expr}
    if isinstance(expr,(list,tuple,set)):return set().union(*(random_symbols(v) for v in expr))
    if getattr(expr,'_stats_joint',False):return set(expr.components)
    if getattr(expr,'_stats_matrix',False):return {expr}
    try:free = scalar(expr).free_symbols
    except (ValueError,TypeError):free = getattr(expr,'free_symbols',set())
    return {REGISTRY[v] for v in free if v in REGISTRY}


def pspace(expr):
    if getattr(expr,'_stats_joint',False) or getattr(expr,'_stats_matrix',False):return expr.pspace
    variables = random_symbols(expr)
    if not variables:raise ValueError('expression contains no random variables')
    from ._spaces import make_space
    return make_space(variables)


def _truth(condition):
    if condition is True or condition == S.true:return True
    if condition is False or condition == S.false:return False
    return None


def _indicator(condition):
    truth = _truth(condition)
    return S.Integer(int(truth)) if truth is not None else S.Piecewise((1,condition),(0,True))


def _finite_indicator(condition,mapping,variables):
    result = condition.subs(mapping,simultaneous=True)
    replacements = {}
    for rv in variables:
        model = rv.distribution
        if model.kind == 'finite' and model._table is not None:
            labels = [v for v in model.dict if v.is_Symbol]
            for a in labels:
                for b in labels:
                    if a != b:replacements[S.Eq(a,b)] = S.false
    return _indicator(result.xreplace(replacements) if replacements else result)


def _groups(variables):
    groups = defaultdict(list)
    for rv in sorted(variables,key=lambda rv:rv.name):groups[rv.group].append(rv)
    return groups


def _worlds(variables):
    choices = []
    count = 1
    for group,values in _groups(variables).items():
        if hasattr(group,'finite_worlds'):
            available = group.finite_worlds(values)
            if available is None:return None
            states = list(available)
        elif len(values) == 1 and values[0].distribution.kind == 'finite' and values[0].distribution._table is not None:
            rv = values[0]
            states = [({rv._atom:value},probability) for value,probability in rv.distribution.dict.items()]
        else:return None
        count *= len(states)
        if count > 1000000:raise RuntimeError('finite probability enumeration exceeds one million joint outcomes')
        choices.append(states)
    def states():
        for combination in product(*choices):
            values,weight = {},S.Integer(1)
            for mapping,probability in combination:values.update(mapping);weight *= probability
            yield values,weight
    return states()


def _bindings(condition):
    if condition is None:return {}
    parts = condition.args if getattr(condition,'class_name','') == 'And' else (condition,)
    bindings = {}
    for part in parts:
        if getattr(part,'class_name','') == 'Equality':
            a,b = part.lhs,part.rhs
            if a in REGISTRY and not random_symbols(b):bindings[a] = b
            elif b in REGISTRY and not random_symbols(a):bindings[b] = a
    for atom,value in bindings.items():
        rv = REGISTRY[atom];model = rv.distribution
        if rv.process is not None or hasattr(rv.group,'joint_moment'):continue
        if model.kind == 'continuous':check(S.And(value>=model.lower,value<=model.upper),'conditioned value is outside the distribution support')
        elif model.kind == 'finite' and model._table is not None and value in model.dict:
            check(model.dict[value] != 0,'conditioning event has zero probability')
        else:check(model.bounded_pdf(value).is_zero is not True,'conditioned value is outside the distribution support')
    return bindings


def _joint_conditioning(expression,condition):
    bindings = _bindings(condition)
    if not bindings:return None
    replacements = {}
    for group,values in _groups(random_symbols(expression)|random_symbols(condition)).items():
        if hasattr(group,'conditioned_replacements'):
            mapping = group.conditioned_replacements(bindings)
            if mapping:replacements.update(mapping)
    if not replacements:return None
    return scalar(expression).subs(replacements,simultaneous=True),condition.subs(replacements,simultaneous=True)


def _condition_set(condition,rv):
    model,x = rv.distribution,rv._atom
    truth = _truth(condition)
    if truth is not None:return model.set if truth else S.EmptySet
    name = getattr(condition,'class_name','')
    if name == 'And':return S.Intersection(model.set,*(_condition_set(part,rv) for part in condition.args))
    if name == 'Or':return S.Union(*(_condition_set(part,rv) for part in condition.args))
    if name == 'Not':return S.Complement(model.set,_condition_set(condition.args[0],rv))
    if condition.is_Relational:
        if condition.rel_op in ('==','!='):
            solution = S.solveset(condition.lhs-condition.rhs,x,domain=S.Reals)
            return S.Intersection(model.set,solution) if condition.rel_op == '==' else S.Complement(model.set,solution)
        return S.Intersection(model.set,S.solve_univariate_inequality(condition,x,relational=False))
    raise NotImplementedError('this event cannot be reduced to a univariate support set')


def _eval_sum_direct(expr,i,a,b):
    """SymPy's eval_sum_direct: a finite sum term by term, after pulling the
    factors and terms free of the index out of it."""
    dif = b-a
    if expr.is_Mul:
        without_i,with_i = expr.as_independent(i)
        if without_i != 1:
            s = _eval_sum_direct(with_i,i,a,b)
            if s:
                r = without_i*s
                if r is not S.nan:return r
        else:
            L,R = expr.as_two_terms()
            if not L.has(i):
                sR = _eval_sum_direct(R,i,a,b)
                if sR:return L*sR
            if not R.has(i):
                sL = _eval_sum_direct(L,i,a,b)
                if sL:return sL*R
    try:expr = S.apart(expr,i)
    except (S.PolynomialError,NotImplementedError,ValueError,TypeError):pass
    if expr.is_Add:
        without_i,with_i = expr.as_independent(i)
        if without_i != 0:
            s = _eval_sum_direct(with_i,i,a,b)
            if s:
                r = without_i*(dif+1)+s
                if r is not S.nan:return r
        else:
            L,R = expr.as_two_terms()
            lsum = _eval_sum_direct(L,i,a,b)
            rsum = _eval_sum_direct(R,i,a,b)
            if None not in (lsum,rsum):
                r = lsum+rsum
                if r is not S.nan:return r
    return S.Add(*[expr.subs(i,a+j) for j in range(int(dif)+1)])


def _direct_range_probability(model,region):
    """SymPy's SingleDiscretePSpace.eval_prob on a bounded Range:
    summation(pmf(n*step), (n, inf/step, sup/step - 1)).doit(), which for
    fewer than 100 terms is eval_sum_direct and is not simplified. Returns
    None outside that case."""
    if model.kind != 'discrete' or not isinstance(region,S.Range) or not region.size.is_Integer:return None
    start,stop,step = region.args
    a,b = start/step,stop/step-1
    if not (a.is_Integer and b.is_Integer):return None
    dif = b-a
    if dif < 0 or dif >= 100:return None
    n = S.Symbol('n',integer=True)
    summand = model.pdf(n*step)
    if n not in summand.free_symbols:return summand*(dif+1)
    if dif == 0:return summand.subs(n,a)
    if isinstance(summand,S.Piecewise):return None
    return _eval_sum_direct(summand,n,a,b)


def _measure(model,region):
    if region == S.EmptySet:return S.Integer(0)
    if region == model.set:return S.Integer(1)
    if isinstance(region,S.Interval):
        a,b = region.start,region.end
        if model.kind == 'continuous':return S.simplify(model.cdf(b)-model.cdf(a))
        lower = S.floor(a)+1 if region.left_open else S.ceiling(a)
        upper = S.ceiling(b)-1 if region.right_open else S.floor(b)
        return S.simplify(model.cdf(upper)-model.cdf(lower-1))
    if isinstance(region,S.FiniteSet):
        return S.Integer(0) if model.kind == 'continuous' else S.Add(*(model.bounded_pdf(v) for v in region.args))
    if isinstance(region,S.Union):return S.Add(*(_measure(model,part) for part in region.args))
    if isinstance(region,S.Complement):return _measure(model,region.args[0])-_measure(model,S.Intersection(*region.args))
    if isinstance(region,S.Range):return model.cdf(region.stop-1)-model.cdf(region.start-1)
    raise NotImplementedError('probability of this support-set form is not implemented')


def _normal_linear(expr):
    variables = random_symbols(expr)
    if not variables:return None
    if all(v.distribution.family == 'Normal' and v.process is None and not hasattr(v.group,'joint_moment') for v in variables):
        atoms = [v._atom for v in variables]
        try:poly = S.Poly(scalar(expr),*atoms)
        except (TypeError,ValueError):return None
        if poly.total_degree() > 1:return None
        mean = scalar(expr).subs({v._atom:v.distribution.args[0] for v in variables},simultaneous=True)
        var = S.Add(*(scalar(expr).diff(v._atom)**2*v.distribution.args[1]**2 for v in variables))
        return mean,S.simplify(var)
    groups = _groups(variables)
    if len(groups) == 1:
        group,values = next(iter(groups.items()))
        if hasattr(group,'normal_linear'):return group.normal_linear(scalar(expr),values)
    return None


def _process_query(operation,expr,condition):
    variables = random_symbols(expr)|random_symbols(condition) if condition is not None else random_symbols(expr)
    processes = {v.process for v in variables if getattr(v,'process',None) is not None}
    if len(processes) == 1:
        return next(iter(processes)).query(operation,expr,condition)
    return None


def _independent_parts(parts):
    """Group Boolean parts into blocks whose random variables share no
    probability space; parts free of random variables form their own block."""
    blocks = []
    for part in parts:
        groups = {v.group for v in random_symbols(part)}
        joined = [block for block in blocks if block[0]&groups]
        rest = [block for block in blocks if not block[0]&groups]
        merged = (set(groups).union(*(b[0] for b in joined)),[p for b in joined for p in b[1]]+[part])
        blocks = rest+[merged]
    return [parts for _,parts in blocks]


def _factored_probability(condition,variables):
    """Exact probability of an And, Or or Not of events over independent
    probability spaces: the product rule P(A & B) = P(A) P(B) and the
    complement rule P(A | B) = 1 - (1 - P(A)) (1 - P(B))."""
    if any(getattr(v,'process',None) is not None for v in variables):return None
    name = getattr(condition,'class_name','')
    if name == 'Not':return S.simplify(1-P(condition.args[0]))
    if name not in ('And','Or'):return None
    blocks = _independent_parts(condition.args)
    if len(blocks) < 2:return None
    def value(parts):
        event = S.And(*parts) if name == 'And' else S.Or(*parts)
        return P(event) if random_symbols(event) else _indicator(event)
    if name == 'And':return S.simplify(S.Mul(*(value(parts) for parts in blocks)))
    return S.simplify(S.expand(1-S.Mul(*(1-value(parts) for parts in blocks))))


def _independent_condition(expr,condition,variables):
    """Drop the conjuncts of a condition that are independent of ``expr``:
    E(X | A & B) = E(X | A) when B lives on probability spaces disjoint from
    those of X and A. None when nothing can be dropped."""
    if any(getattr(v,'process',None) is not None for v in variables):return None
    if not random_symbols(condition):return None
    parts = condition.args if getattr(condition,'class_name','') == 'And' else (condition,)
    target = {v.group for v in random_symbols(expr)}
    kept,dropped = [],[]
    for block in _independent_parts(parts):
        groups = {v.group for part in block for v in random_symbols(part)}
        (dropped if groups and not groups&target else kept).append(block)
    if not dropped:return None
    check(S.And(*(P(S.And(*block)) != 0 for block in dropped)),'conditioning event has zero probability')
    kept = [part for block in kept for part in block]
    return S.And(*kept) if kept else S.true


def _observed(condition):
    """A random symbol used as a condition means that its value is observed:
    SymPy reads ``given(expr, Y)`` as ``given(expr, Eq(Y, Y.symbol))``. Return
    the substitution of the variable by its symbol, or None."""
    if condition is None:return None
    if isinstance(condition,RandomSymbol):rv = condition
    else:
        try:atom = scalar(condition)
        except (ValueError,TypeError):return None
        rv = REGISTRY.get(atom) if getattr(atom,'is_Symbol',False) or getattr(atom,'is_Dummy',False) else None
        if rv is None:return None
    if getattr(rv,'process',None) is not None or getattr(rv,'_stats_joint',False) or getattr(rv,'_stats_matrix',False):return None
    if rv.group is not rv._atom:raise NotImplementedError('conditioning on a component of a joint distribution is not supported')
    return {rv._atom:rv.symbol}


def P(condition,given_condition=None,numsamples=None,evaluate=True,**kwargs):
    if not evaluate:
        from ._symbolic import Probability
        return Probability(condition,given_condition)
    if numsamples is not None:
        values = sample(condition,given_condition,size=(operator.index(numsamples),),**kwargs)
        return S.Float(float(np.mean(values)))
    observed = _observed(given_condition)
    if observed is not None:
        event = scalar(condition).xreplace(observed)
        return P(event,None,**kwargs) if random_symbols(event) else _indicator(event)
    condition = scalar(condition)
    if given_condition is not None:given_condition = scalar(given_condition)
    special = _process_query('P',condition,given_condition)
    if special is not None:return special
    truth = _truth(condition)
    if truth is not None:return S.Integer(int(truth))
    variables = random_symbols(condition)|random_symbols(given_condition)
    if not variables:
        if condition.free_symbols:
            from ._symbolic import Probability
            return Probability(condition,given_condition)
        raise ValueError('probability requires a Boolean event')
    worlds = _worlds(variables)
    if worlds is not None:
        numerator,denominator = S.Integer(0),S.Integer(0)
        for mapping,weight in worlds:
            gate = S.Integer(1) if given_condition is None else _finite_indicator(given_condition,mapping,variables)
            denominator += weight*gate
            numerator += weight*gate*_finite_indicator(condition,mapping,variables)
        check(denominator != 0,'conditioning event has zero probability')
        return S.simplify(numerator/denominator)
    conditioned = _joint_conditioning(condition,given_condition)
    if conditioned is not None:
        event,gate = conditioned
        return P(event,None if _truth(gate) is True else gate)
    bindings = _bindings(given_condition)
    if bindings:
        condition = condition.subs(bindings,simultaneous=True)
        given_condition = given_condition.subs(bindings,simultaneous=True)
        if _truth(given_condition) is True:return P(condition)
    if given_condition is not None:
        denominator = P(given_condition)
        check(denominator != 0,'conditioning event has zero probability')
        return S.simplify(P(S.And(condition,given_condition))/denominator)
    factored = _factored_probability(condition,variables)
    if factored is not None:return factored
    # Linear combinations of independent or jointly normal variables reduce
    # to a single exact normal distribution before event integration.
    if condition.is_Relational and len(variables)>1:
        expression = condition.lhs-condition.rhs
        normal = _normal_linear(expression)
        if normal is not None:
            from ._univariate import Normal
            if normal[1] == 0:return _indicator(S.Rel(normal[0],0,condition.rel_op))
            Z = Normal('_linear_event',normal[0],S.sqrt(normal[1]))
            return P(S.Rel(Z,0,condition.rel_op))
    if len(variables) == 1:
        rv = next(iter(variables))
        try:
            region = _condition_set(condition,rv)
            # SymPy sums the pmf over a bounded range without simplifying
            # the sum ((9*t**3/2 + 9*t**2/2 + 3*t + 1)*exp(-3*t)).
            direct = _direct_range_probability(rv.distribution,region)
            if direct is not None:return direct
            return S.simplify(_measure(rv.distribution,region))
        except NotImplementedError:pass
    # A general event is a well-defined indicator expectation, not an
    # unverified numerical probability or a guessed independence product.
    return _expect(S.Piecewise((1,condition),(0,True)),variables)


def _expect(expr,variables=None):
    expr = scalar(expr)
    variables = random_symbols(expr) if variables is None else set(variables)
    if not variables:return expr
    worlds = _worlds(variables)
    if worlds is not None:return S.simplify(S.Add(*(weight*expr.subs(values,simultaneous=True) for values,weight in worlds)))
    atoms = [v._atom for v in sorted(variables,key=lambda v:v.name)]
    try:polynomial = S.Poly(expr,*atoms)
    except (ValueError,TypeError):polynomial = None
    if polynomial is not None:
        terms = []
        for powers,coefficient in polynomial.terms():
            grouped = defaultdict(dict)
            for atom,power in zip(atoms,powers):
                if power:
                    rv = REGISTRY[atom];grouped[rv.group][rv] = power
            value = coefficient
            for group,factors in grouped.items():
                if hasattr(group,'joint_moment'):value *= group.joint_moment(factors)
                elif len(factors) == 1:
                    rv,power = next(iter(factors.items()));value *= rv.distribution.moment(S.Integer(power))
                else:raise ValueError('dependent random variables need a joint probability model')
            terms.append(value)
        return S.simplify(S.Add(*terms))
    if expr.is_Add:return S.Add(*(_expect(part) for part in expr.args))
    if expr.is_Mul and not any(getattr(v,'process',None) is not None for v in variables):
        # E(f(X) g(Y)) = E(f(X)) E(g(Y)) over independent probability spaces
        blocks = _independent_parts(expr.args)
        if len(blocks) > 1:return S.simplify(S.Mul(*(_expect(S.Mul(*block)) for block in blocks)))
    groups = _groups(variables)
    if len(groups) == 1:
        group,values = next(iter(groups.items()))
        if hasattr(group,'expectation'):return group.expectation(expr,values)
    value = expr
    for group,rvs in groups.items():
        if len(rvs) != 1:raise NotImplementedError('nonpolynomial expectation over this dependent space is not implemented')
        rv = rvs[0];model = rv.distribution
        if model.kind == 'finite':value = S.Add(*(p*value.subs(rv._atom,k) for k,p in model.dict.items()))
        else:
            x = S.Dummy(rv.name,real=True)
            value = value.subs(rv._atom,x)*model.pdf(x)
            value = sum_expr(value,x,model.lower,model.upper) if model.kind == 'discrete' else integral(value,x,model.lower,model.upper)
    return value


def _truncated_moment(model,n,a,b):
    family,args = model.family,model.args
    if family == 'Normal':
        mu,sigma = args;alpha,beta = (a-mu)/sigma,(b-mu)/sigma
        phi = lambda x:S.Integer(0) if x in (S.oo,-S.oo) else S.exp(-x*x/2)/S.sqrt(2*S.pi)
        sequence = [model.cdf(b)-model.cdf(a),phi(alpha)-phi(beta)]
        for j in range(2,n+1):
            boundary = (S.Integer(0) if alpha == -S.oo else alpha**(j-1)*phi(alpha))-(S.Integer(0) if beta == S.oo else beta**(j-1)*phi(beta))
            sequence.append(boundary+(j-1)*sequence[j-2])
        return S.Add(*(S.binomial(n,j)*mu**(n-j)*sigma**j*sequence[j] for j in range(n+1)))
    if family == 'Uniform':return (b**(n+1)-a**(n+1))/((n+1)*(args[1]-args[0]))
    if family == 'Beta':
        alpha,beta = args
        return S.beta(alpha+n,beta)/S.beta(alpha,beta)*S.betainc_regularized(alpha+n,beta,a,b)
    if family in ('Gamma','Exponential','ChiSquared'):
        if family == 'Gamma':k,theta = args
        elif family == 'Exponential':k,theta = S.Integer(1),1/args[0]
        else:k,theta = args[0]/2,S.Integer(2)
        gamma_at = lambda x:S.gamma(k+n) if x == S.oo else S.lowergamma(k+n,x/theta)
        return theta**n*(gamma_at(b)-gamma_at(a))/S.gamma(k)
    x = S.Dummy('truncated_x',real=True)
    return integral(x**n*model.pdf(x),x,a,b)


def E(expr,condition=None,numsamples=None,evaluate=True,**kwargs):
    if not evaluate:
        from ._symbolic import Expectation
        return Expectation(expr,condition)
    if numsamples is not None:return S.Float(float(np.mean(sample(expr,condition,size=(operator.index(numsamples),),**kwargs))))
    if getattr(expr,'_stats_matrix',False):return expr.mean()
    if getattr(expr,'_stats_joint',False):return S.Matrix([E(v,condition) for v in expr.components])
    if isinstance(expr,S.DenseMatrix) or getattr(expr,'_symox_sparse',False):return expr.applyfunc(lambda value:E(value,condition,**kwargs))
    observed = _observed(condition)
    if observed is not None:return E(scalar(expr).xreplace(observed),None,**kwargs)
    expr = scalar(expr)
    if condition is not None:condition = scalar(condition)
    special = _process_query('E',expr,condition)
    if special is not None:return special
    if condition is None:return _expect(expr)
    variables = random_symbols(expr)|random_symbols(condition)
    reduced = _independent_condition(expr,condition,variables)
    if reduced is not None:return E(expr,None if reduced == S.true else reduced)
    worlds = _worlds(variables)
    if worlds is not None:
        numerator,denominator = S.Integer(0),S.Integer(0)
        for values,weight in worlds:
            gate = _finite_indicator(condition,values,variables)
            denominator += weight*gate
            numerator += weight*gate*expr.subs(values,simultaneous=True)
        check(denominator != 0,'conditioning event has zero probability')
        return S.simplify(numerator/denominator)
    conditioned = _joint_conditioning(expr,condition)
    if conditioned is not None:
        value,gate = conditioned
        return E(value,None if _truth(gate) is True else gate)
    bindings = _bindings(condition)
    if bindings:
        reduced = condition.subs(bindings,simultaneous=True)
        if _truth(reduced) is True:return E(expr.subs(bindings,simultaneous=True))
    denominator = P(condition)
    check(denominator != 0,'conditioning event has zero probability')
    if len(variables) == 1:
        rv = next(iter(variables));region = _condition_set(condition,rv)
        intervals = region.args if isinstance(region,S.Union) else (region,)
        try:poly = S.Poly(expr,rv._atom)
        except ValueError:poly = None
        total = S.Integer(0)
        for interval in intervals:
            if not isinstance(interval,S.Interval):break
            if poly is not None and rv.distribution.kind == 'continuous':
                total += S.Add(*(coefficient*_truncated_moment(rv.distribution,int(powers[0]),interval.start,interval.end) for powers,coefficient in poly.terms()))
            else:
                x = S.Dummy('conditional_x',real=True)
                total += integral(expr.subs(rv._atom,x)*rv.distribution.pdf(x),x,interval.start,interval.end)
        else:return S.simplify(total/denominator)
    return S.simplify(_expect(expr*S.Piecewise((1,condition),(0,True)),variables)/denominator)


def variance(X,condition=None,**kwargs):
    if kwargs.pop('evaluate',True) is False:
        from ._symbolic import Variance
        return Variance(X,condition)
    if getattr(X,'_stats_matrix',False):return X.variance()
    mean = E(X,condition,**kwargs)
    return S.simplify(E(scalar(X)**2,condition,**kwargs)-mean**2)


def std(X,condition=None,**kwargs):return S.sqrt(variance(X,condition,**kwargs))
def covariance(X,Y,condition=None,**kwargs):
    if kwargs.pop('evaluate',True) is False:
        from ._symbolic import Covariance
        return Covariance(X,Y,condition)
    return S.simplify(E(scalar(X)*scalar(Y),condition,**kwargs)-E(X,condition,**kwargs)*E(Y,condition,**kwargs))
def correlation(X,Y,condition=None,**kwargs):return S.simplify(covariance(X,Y,condition,**kwargs)/(std(X,condition,**kwargs)*std(Y,condition,**kwargs)))


def moment(X,n,c=0,condition=None,*,evaluate=True,**kwargs):
    if not evaluate:
        from ._symbolic import Moment
        return Moment(X,n,c,condition)
    n = scalar(n);check(S.And(n>=0,S.Contains(n,S.Integers)),'moment order must be a nonnegative integer')
    return E((scalar(X)-scalar(c))**n,condition,**kwargs)


def cmoment(X,n,condition=None,*,evaluate=True,**kwargs):
    if not evaluate:
        from ._symbolic import CentralMoment
        return CentralMoment(X,n,condition)
    return moment(X,n,E(X,condition),condition,**kwargs)


def smoment(X,n,condition=None,**kwargs):return S.simplify(cmoment(X,n,condition,**kwargs)/std(X,condition,**kwargs)**scalar(n))
def skewness(X,condition=None,**kwargs):return smoment(X,3,condition,**kwargs)
def kurtosis(X,condition=None,**kwargs):return smoment(X,4,condition,**kwargs)


def factorial_moment(X,n,condition=None,**kwargs):
    n = operator.index(n)
    if n < 0:raise ValueError('factorial moment order must be nonnegative')
    return E(S.Mul(*(scalar(X)-i for i in range(n))),condition,**kwargs)


def coskewness(X,Y,Z,condition=None,**kwargs):
    numerator = E((scalar(X)-E(X,condition))*(scalar(Y)-E(Y,condition))*(scalar(Z)-E(Z,condition)),condition)
    return S.simplify(numerator/(std(X,condition)*std(Y,condition)*std(Z,condition)))


def independent(a,b):
    left,right = random_symbols(a),random_symbols(b)
    if not left or not right:return True
    if not {v.group for v in left}&{v.group for v in right}:return True
    if _normal_linear(a) is not None and _normal_linear(b) is not None:return covariance(a,b) == 0
    if len(left) == len(right) == 1:
        x,y = next(iter(left)),next(iter(right))
        if x.group is y.group and hasattr(x.group,'independent_components'):return x.group.independent_components(x,y)
    return False


def dependent(a,b):return not independent(a,b)


def density(expr,condition=None,evaluate=True,numsamples=None,**kwargs):
    if numsamples is not None:return sampling_density(expr,condition,numsamples=numsamples,**kwargs)
    if condition is not None:expr = given(expr,condition)
    if isinstance(expr,RandomSymbol):return expr.distribution
    if getattr(expr,'_stats_joint',False) or getattr(expr,'_stats_matrix',False):return expr.distribution
    expression = scalar(expr);variables = random_symbols(expression)
    worlds = _worlds(variables)
    if worlds is not None:
        table = {}
        for values,weight in worlds:
            value = expression.subs(values,simultaneous=True)
            table[value] = table.get(value,S.Integer(0))+weight
        from ._spaces import FiniteDensity
        return FiniteDensity({value:S.simplify(weight) for value,weight in table.items()})
    normal = _normal_linear(expression)
    if normal is not None:
        from ._univariate import continuous_model
        if normal[1] == 0:return {normal[0]:S.Integer(1)}
        return continuous_model('Normal',(normal[0],S.sqrt(normal[1])))
    if len(variables) == 1:
        rv = next(iter(variables));x,y = rv._atom,S.Dummy('density_value',real=True)
        roots = S.solve(S.Eq(expression,y),x)
        if roots:
            terms = []
            for root in roots:
                jacobian = S.Abs(expression.diff(x).subs(x,root))
                terms.append(S.Piecewise((rv.distribution.pdf(root)/jacobian,S.And(root>=rv.distribution.lower,root<=rv.distribution.upper)),(0,True)))
            return S.Lambda(y,S.Add(*terms))
    y = S.Dummy('density_value',real=True)
    return S.Lambda(y,E(S.DiracDelta(y-expression)))


def cdf(expr,condition=None,evaluate=True,**kwargs):
    model = density(expr,condition,evaluate=evaluate)
    if isinstance(model,dict):
        total,table = S.Integer(0),{}
        for value,p in sorted(model.items(),key=lambda pair:float(pair[0])):total += p;table[value] = S.simplify(total)
        return table
    if isinstance(model,Distribution) and model.kind == 'finite' and model._table is not None:
        total,table = S.Integer(0),{}
        for value,p in sorted(model.dict.items(),key=lambda pair:float(pair[0])):total += p;table[value] = S.simplify(total)
        return table
    mass = getattr(model,'_symbolic_mass',None) if isinstance(model,Distribution) else None
    if mass is not None and model._table is None:
        # SymPy SingleFinitePSpace.compute_cdf for a symbolic size.
        k,ki = S.Dummy('k'),S.Dummy('ki')
        term = S.Piecewise((mass(ki),S.And(ki >= model.lower,ki <= model.upper)),(0,True))
        return S.Lambda(k,S.Sum(term,(ki,model.lower,k)))
    x = S.Dummy('z',real=True)
    if isinstance(model,Distribution):return S.Lambda(x,model.cdf(x,**kwargs))
    t = S.Dummy('t',real=True)
    return S.Lambda(x,integral(model(t),t,-S.oo,x))


def quantile(expr,evaluate=True,**kwargs):
    model = density(expr)
    if isinstance(model,dict):
        from ._univariate import FiniteDistributionHandmade
        model = FiniteDistributionHandmade(model)
    p = S.Dummy('p',real=True)
    return S.Lambda(p,model.quantile(p))


def median(X,evaluate=True,**kwargs):
    model = density(X)
    if isinstance(model,dict) or getattr(model,'kind',None) == 'finite':
        table = model if isinstance(model,dict) else model.dict
        total,values = S.Integer(0),[]
        for value,probability in sorted(table.items(),key=lambda pair:float(pair[0])):
            after = total+probability
            if _truth(S.And(total<=HALF,after>=HALF)) is not False:values.append(value)
            total = after
        return S.FiniteSet(*values)
    return S.FiniteSet(model.quantile(HALF))


def moment_generating_function(expr,condition=None,evaluate=True,**kwargs):
    model = density(expr,condition,evaluate=evaluate);t = S.Dummy('t',real=True)
    if isinstance(model,dict):return S.Lambda(t,S.Add(*(p*S.exp(t*v) for v,p in model.items())))
    if not isinstance(model,Distribution):
        x = S.Dummy('mgf_x',real=True)
        return S.Lambda(t,integral(S.exp(t*x)*model(x),x,-S.oo,S.oo))
    return S.Lambda(t,model.mgf(t))


def characteristic_function(expr,condition=None,evaluate=True,**kwargs):
    model = density(expr,condition,evaluate=evaluate);t = S.Dummy('t',real=True)
    if isinstance(model,dict):return S.Lambda(t,S.Add(*(p*S.exp(S.I*t*v) for v,p in model.items())))
    if not isinstance(model,Distribution):
        x = S.Dummy('cf_x',real=True)
        return S.Lambda(t,integral(S.exp(S.I*t*x)*model(x),x,-S.oo,S.oo))
    if model.family == 'Cauchy':mu,scale = model.args;value = S.exp(S.I*mu*t-scale*S.Abs(t))
    elif model.family == 'Levy':mu,c = model.args;value = S.exp(S.I*mu*t-S.sqrt(-2*S.I*c*t))
    elif model.family == 'StudentT':
        nu, = model.args
        value = S.Piecewise((1,S.Eq(t,0)),(2**(1-nu/2)*(S.sqrt(nu)*S.Abs(t))**(nu/2)*S.besselk(nu/2,S.sqrt(nu)*S.Abs(t))/S.gamma(nu/2),True))
    else:value = model.mgf(S.I*t)
    return S.Lambda(t,value)


def entropy(expr,condition=None,**kwargs):
    base = kwargs.pop('base',kwargs.pop('b',None))
    if kwargs:raise TypeError('unknown entropy options')
    model = density(expr,condition)
    if isinstance(model,dict) or getattr(model,'kind',None) == 'finite':
        table = model if isinstance(model,dict) else model.dict
        result = -S.Add(*(p*S.log(p) for p in table.values() if p != 0))
    elif not isinstance(model,Distribution):
        x = S.Dummy('entropy_x',real=True);result = integral(-model(x)*S.log(model(x)),x,-S.oo,S.oo)
    elif model._entropy is not None:result = model._entropy
    else:
        x = S.Dummy('entropy_x',real=True);body = -model.pdf(x)*S.log(model.pdf(x))
        result = sum_expr(body,x,model.lower,model.upper) if model.kind == 'discrete' else integral(body,x,model.lower,model.upper)
    return S.simplify(result if base is None else result/S.log(scalar(base)))


H = entropy


def where(condition,given_condition=None,**kwargs):
    condition = scalar(condition)
    if given_condition is not None:condition = S.And(condition,scalar(given_condition))
    variables = sorted(random_symbols(condition),key=lambda rv:rv.name)
    if len(variables) == 1:
        rv = variables[0]
        if rv.distribution.kind == 'finite' and rv.distribution._table is not None:
            support = S.FiniteSet(*(value for value in rv.distribution.dict if _truth(condition.subs(rv._atom,value)) is True))
        else:support = _condition_set(condition,rv)
        return Domain(variables,support,condition)
    worlds = _worlds(variables)
    if worlds is not None:
        support = S.FiniteSet(*(S.Tuple(*(values[v._atom] for v in variables)) for values,weight in worlds if _truth(condition.subs(values,simultaneous=True)) is True))
        return Domain(variables,support,condition)
    return Domain(variables,condition=condition)


def given(expr,condition=None,**kwargs):
    if condition is None:return expr
    if isinstance(expr,(tuple,list)) or getattr(expr,'class_name','') == 'Tuple':return _given_tuple(expr,condition,**kwargs)
    observed = _observed(condition)
    if observed is not None:return scalar(expr).xreplace(observed)
    condition = scalar(condition)
    expression = scalar(expr)
    if not random_symbols(condition):return expr
    if not {v.group for v in random_symbols(expression)}&{v.group for v in random_symbols(condition)}:return expr
    conditioned = _joint_conditioning(expression,condition)
    if conditioned is not None:
        value,gate = conditioned
        return value if _truth(gate) is True else given(value,gate,**kwargs)
    bindings = _bindings(condition)
    if bindings and _truth(condition.subs(bindings,simultaneous=True)) is True:return expression.subs(bindings,simultaneous=True)
    variables = random_symbols(expression)|random_symbols(condition)
    if len(variables) != 1:raise NotImplementedError('conditional random-variable construction requires one scalar source')
    source = next(iter(variables));model0 = source.distribution
    if model0.kind == 'finite':
        from ._univariate import FiniteRV
        normalizer = P(condition);check(normalizer != 0,'conditioning event has zero probability')
        table = {}
        for value,weight in model0.dict.items():
            outcome = expression.subs(source._atom,value)
            probability = S.simplify(weight*_indicator(condition.subs(source._atom,value))/normalizer)
            if probability != 0:table[outcome] = table.get(outcome,S.Integer(0))+probability
        return FiniteRV(source.name+'|given',table)
    normalizer = P(condition);check(normalizer != 0,'conditioning event has zero probability')
    model = Distribution(source._atom,condition);model.family = 'Conditional';model.kind = model0.kind
    model._bound_symbols = {source._atom}
    model.lower,model.upper = model0.lower,model0.upper
    model._pdf = lambda x:model0.pdf(x)*_indicator(condition.subs(source._atom,x))/normalizer
    model._moment = lambda n:E(source._atom**n,condition)
    model._sampler = lambda rng,size:_sample_values(source._atom,condition,size,rng)
    model._support = _condition_set(condition,source)
    result = variable(source.name+'|given',model)
    return result if expression == source._atom else expression.subs(source._atom,scalar(result))


def _given_tuple(items,condition,**kwargs):
    """Condition every entry of a tuple on the same event. Each entry is
    conditioned on its own only when this preserves the joint law: at most
    one entry may depend on the event, and it may share no probability space
    with the other entries."""
    items = [scalar(v) for v in items]
    observed = _observed(condition)
    if observed is not None:return S.Tuple(*(v.xreplace(observed) for v in items))
    condition = scalar(condition)
    event = {v.group for v in random_symbols(condition)}
    groups = [{v.group for v in random_symbols(v)} for v in items]
    touched = [i for i,g in enumerate(groups) if g&event]
    if len(touched)>1 or any(groups[i]&groups[j] for i in touched for j in range(len(items)) if j != i):
        # SymPy: the random symbols of the tuple and the event move to the
        # conditional space of their joint space, so the entries stay
        # jointly distributed (and dependent).
        swap = _conditional_space(items,condition)
        return S.Tuple(*(v.xreplace(swap) for v in items))
    return S.Tuple(*(given(v,condition,**kwargs) if i in touched else v for i,v in enumerate(items)))


class _ConditionalFiniteSpace:
    """The joint law of finite random symbols given an event (SymPy
    ``FinitePSpace.conditional_space``): explicit weighted joint outcomes."""
    def __init__(self,outcomes):
        self.outcomes = outcomes
    def finite_worlds(self,values):
        table = {}
        for outcome,weight in self.outcomes:
            key = tuple(outcome[rv] for rv in values)
            table[key] = table.get(key,S.Integer(0))+weight
        return [({rv._atom:k for rv,k in zip(values,key)},weight) for key,weight in table.items()]


class _ConditionalContinuousSpace:
    """The joint law of continuous random symbols given an equation
    ``g = 0`` linear in one of them: SymPy's ``DiracDelta(g)`` conditional
    space, ``f(x, y(x))/|dg/dy|/Z`` after integrating out ``y``."""
    def __init__(self,variables,equation,normalizer):
        self.variables,self.equation,self.normalizer = variables,equation,normalizer
        # The unconditional laws of the variables.
        self.laws = {v:v.distribution for v in variables}
    def eliminate(self,keep=None):
        """(free, eliminated, solution, density) with ``keep`` free."""
        for eliminated in reversed(self.variables):
            if eliminated is keep:continue
            x = eliminated._atom;slope = S.diff(self.equation,x)
            if slope == 0 or slope.has(x) or random_symbols(slope)-set(self.variables):continue
            if any(slope.has(v._atom) for v in self.variables):continue
            solution = S.solve(self.equation,x)
            if len(solution) != 1:continue
            free = [v for v in self.variables if v is not eliminated]
            density = S.Mul(*(self.laws[v].pdf(v._atom) for v in free))
            density *= self.laws[eliminated].bounded_pdf(solution[0])/abs(slope)
            return free,eliminated,solution[0],density
        return None
    def expectation(self,expr,values=None):
        free,eliminated,solution,density = self.eliminate()
        body = scalar(expr).xreplace({eliminated._atom:solution})
        coordinates = [(rv,S.Dummy(rv.name,real=True)) for rv in free]
        body = (body*density/self.normalizer).xreplace({rv._atom:x for rv,x in coordinates})
        limits = [(x,rv.distribution.lower,rv.distribution.upper) for rv,x in coordinates]
        return _erf_simplify(S.integrate(body,*limits))
    def joint_moment(self,factors):
        return self.expectation(S.Mul(*(rv._atom**power for rv,power in factors.items())))
    def marginal(self,rv):
        free,eliminated,solution,density = self.eliminate(keep=rv)
        rest = [(w,S.Dummy(w.name,real=True)) for w in free if w is not rv]
        value = density.xreplace({w._atom:d for w,d in rest})/self.normalizer
        if rest:value = S.integrate(value,*[(d,w.distribution.lower,w.distribution.upper) for w,d in rest])
        return _erf_simplify(value)


def _erf_simplify(value):
    """simplify, with ``erfc`` written through ``erf`` so that the Gaussian
    integrals ``erf(a) + erfc(a)`` collapse."""
    if value.has(S.erfc):value = value.rewrite(S.erf)
    return S.simplify(value)


def _conditional_space(items,condition):
    """The map from the random symbols of ``items`` and ``condition`` (and
    their joint spaces) to their copies in the conditional space."""
    variables = set(random_symbols(condition))
    for item in items:variables |= random_symbols(item)
    variables = sorted(variables,key=lambda v:v.name)
    if any(v.process is not None for v in variables):
        raise NotImplementedError('conditioning process values of a tuple jointly is not supported')
    worlds = _worlds(variables)
    if worlds is not None:
        outcomes,total = [],S.Integer(0)
        for values,weight in worlds:
            if _truth(condition.subs(values,simultaneous=True)) is True:
                outcomes.append(({rv:values[rv._atom] for rv in variables},weight));total += weight
        check(total != 0,'conditioning event has zero probability')
        outcomes = [(outcome,weight/total) for outcome,weight in outcomes]
        space = _ConditionalFiniteSpace(outcomes)
        swap,copies = {},{}
        for rv in variables:
            table = {}
            for outcome,weight in outcomes:table[outcome[rv]] = table.get(outcome[rv],S.Integer(0))+weight
            model = Distribution();model.kind = 'finite';model.family = 'FiniteConditional';model._table = table
            copies[rv] = variable(rv.name,model,group=space);swap[rv._atom] = copies[rv]._atom
        space.outcomes = [({copies[rv]:value for rv,value in outcome.items()},weight) for outcome,weight in outcomes]
        return swap
    if (len(variables) > 1 and all(v.distribution.kind == 'continuous' and v.group is v._atom for v in variables)
            and getattr(condition,'class_name','') == 'Equality'):
        space = _ConditionalContinuousSpace(variables,condition.lhs-condition.rhs,S.Integer(1))
        parts = space.eliminate()
        if parts is not None:
            free,eliminated,solution,density = parts
            coordinates = [(v,S.Dummy(v.name,real=True)) for v in free]
            limits = [(c,v.distribution.lower,v.distribution.upper) for v,c in coordinates]
            normalizer = _erf_simplify(S.integrate(density.xreplace({v._atom:c for v,c in coordinates}),*limits))
            if not normalizer.has(S.Integral):
                check(normalizer != 0,'conditioning event has zero probability')
                space.normalizer = normalizer
                marginals = {v:space.marginal(v) for v in variables}
                swap = {}
                for v in variables:
                    model = Distribution();model.family = 'ContinuousConditional'
                    model.lower,model.upper = v.distribution.lower,v.distribution.upper
                    model._pdf = lambda z,m=marginals[v],a=v._atom:m.subs(a,z)
                    swap[v._atom] = variable(v.name,model,group=space)._atom
                space.variables = [REGISTRY[swap[v._atom]] for v in variables]
                space.laws = {REGISTRY[swap[v._atom]]:v.distribution for v in variables}
                space.equation = space.equation.xreplace(swap)
                return swap
    if all(v.distribution.kind == 'continuous' for v in variables) and not getattr(condition,'class_name','') == 'Equality':
        raise NotImplementedError('Multivariate Inequalities not yet implemented')
    raise NotImplementedError('conditioning several dependent entries of a tuple jointly is not supported')


def _sample_values(expr,condition,size,rng):
    if condition is not None:
        conditioned = _joint_conditioning(expr,condition)
        if conditioned is not None:
            value,gate = conditioned
            return _sample_values(value,None if _truth(gate) is True else gate,size,rng)
    variables = random_symbols(expr)|random_symbols(condition)
    groups = _groups(variables)
    bindings = _bindings(condition)
    for atom,value in bindings.items():
        rv = REGISTRY[atom]
        if rv.process is None and not hasattr(rv.group,'joint_moment'):
            check(rv.distribution.bounded_pdf(value).is_zero is not True,'conditioned value is outside the distribution support')
    samples = []
    count = int(np.prod(size)) if size else 1
    for i in range(count):
        for trial in range(100000):
            mapping = {}
            for group,rvs in groups.items():
                if getattr(group,'_stats_process',False):mapping.update(group.sample_mapping(rng,rvs,condition=condition))
                elif hasattr(group,'sample_mapping'):mapping.update(group.sample_mapping(rng,rvs))
                else:
                    for rv in rvs:mapping[rv._atom] = bindings.get(rv._atom) if rv._atom in bindings else scalar(rv.distribution.sample(rng,()))
            if condition is not None and _finite_indicator(condition,mapping,variables) != 1:continue
            result = scalar(expr).subs(mapping,simultaneous=True)
            truth = _truth(result)
            if truth is not None:samples.append(truth)
            else:
                try:samples.append(int(result) if result.is_Integer else float(result.evalf()))
                except (ValueError,TypeError):samples.append(result)
            break
        else:raise RuntimeError('conditional sampler exceeded its rejection limit')
    result = np.asarray(samples).reshape(size)
    return result.item() if not size else result


def sample_iter(expr,condition=None,size=(),library='scipy',numsamples=S.oo,seed=None,**kwargs):
    if library not in ('numpy','scipy'):raise NotImplementedError('sampling supports the NumPy/scipy compatibility selectors')
    if kwargs:raise TypeError('unknown sampling options')
    if isinstance(size,int):size = (size,)
    size = tuple(operator.index(n) for n in size)
    if any(n < 0 for n in size):raise ValueError('sample dimensions must be nonnegative')
    rng = seed if isinstance(seed,np.random.Generator) else np.random.default_rng(seed)
    if condition is not None:condition = scalar(condition)
    count = 0
    while numsamples == S.oo or count < operator.index(numsamples):
        if (getattr(expr,'_stats_joint',False) or getattr(expr,'_stats_matrix',False)) and condition is None:yield expr.sample(rng,size)
        elif isinstance(expr,RandomSymbol) and condition is None and expr.process is None:yield expr.distribution.sample(rng,size)
        else:yield _sample_values(expr,condition,size,rng)
        count += 1


def sample(expr,condition=None,size=(),library='scipy',numsamples=1,seed=None,**kwargs):
    iterator = sample_iter(expr,condition,size,library,numsamples,seed,**kwargs)
    return next(iterator) if numsamples == 1 else list(iterator)


def sampling_density(expr,given_condition=None,library='scipy',numsamples=1,seed=None,**kwargs):
    count = operator.index(numsamples)
    if count < 1:raise ValueError('sampling density needs a positive sample count')
    values = np.asarray(sample(expr,given_condition,size=(count,),library=library,seed=seed,**kwargs))
    variables = random_symbols(expr)
    if all(rv.distribution.kind != 'continuous' for rv in variables):
        unique,counts = np.unique(values,return_counts=True)
        return {scalar(v.item() if hasattr(v,'item') else v):S.Rational(int(n),count) for v,n in zip(unique,counts)}
    x = S.Dummy('kde_x',real=True)
    bandwidth = 1.06*float(np.std(values,ddof=1))*count**(-.2) if count>1 else 0
    if bandwidth == 0:return S.Lambda(x,S.DiracDelta(x-scalar(float(values[0]))))
    width = S.Float(bandwidth)
    return S.Lambda(x,S.Add(*(S.exp(-(x-S.Float(float(v)))**2/(2*width**2)) for v in values))/(count*width*S.sqrt(2*S.pi)))


__all__ = ['P','E','H','density','where','given','sample','cdf','median','characteristic_function',
 'pspace','sample_iter','variance','std','skewness','kurtosis','covariance','dependent','entropy',
 'independent','random_symbols','correlation','factorial_moment','moment','cmoment','sampling_density',
 'moment_generating_function','smoment','quantile','coskewness']
