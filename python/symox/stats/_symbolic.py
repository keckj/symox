"""Held statistical operators with executable exact counterparts."""
from functools import wraps
import symox as S
from ._base import scalar

HELD = {}


def _native(arg):
    if isinstance(arg,StatisticalExpression):return arg._expr
    return arg if isinstance(arg,S.Matrix) else scalar(arg)


def _held_nodes(expr):
    """The held statistical nodes inside a native expression."""
    found,seen,stack = [],set(),[expr]
    while stack:
        node = stack.pop()
        if not isinstance(node,S.Basic) or node in seen:continue
        seen.add(node)
        if node in HELD:found.append(node)
        else:stack.extend(getattr(node,'args',()))
    return found


def _map_held(expr,function):
    """Apply ``function`` to every held node of ``expr`` (inputs first),
    keeping the held registry consistent."""
    if isinstance(expr,StatisticalExpression):expr = expr._expr
    if not isinstance(expr,S.Basic):return expr
    nodes = _held_nodes(expr)
    if not nodes:return expr
    return expr.xreplace({node:_native(function(HELD[node])) for node in nodes})


def _rebuild(held,inputs):
    return type(held)(*inputs,**held.options)


def _held_replace(expr,mapping):
    """xreplace that also rebuilds the held statistical nodes it reaches."""
    expr = _native(expr) if expr is not None else None
    if expr is None:return None
    def visit(held):return _rebuild(held,[_held_replace(i,mapping) if i is not None else None for i in held.inputs])
    return _map_held(expr,visit).xreplace(mapping)


def _is_random(expr):
    from ._queries import random_symbols
    return bool(random_symbols(_native(expr)))


_TARGETS = ('Expectation','Probability','Integral','Sum')


def _target_name(target):
    name = target if isinstance(target,str) else getattr(target,'__name__',None) or getattr(target,'class_name',None)
    return name if name in _TARGETS else None


def rewrite_held(expr,target):
    """Rewrite the held statistical nodes of ``expr`` (inputs first), as the
    reference rewrite does bottom-up."""
    name = _target_name(target)
    if name is None:return _native(expr)
    def visit(held):
        rebuilt = _rebuild(held,[rewrite_held(i,name) if i is not None else None for i in held.inputs])
        return rebuilt._rewrite(name)
    return _map_held(expr,visit)


def _variables(expr):
    from ._queries import random_symbols
    return sorted(random_symbols(_native(expr)),key=lambda v:v.name)


def _binder(rv,body,lower,upper):
    return S.Sum(body,(rv._atom if hasattr(rv,'_atom') else rv,lower,upper)) if rv.distribution.kind == 'discrete' else S.Integral(body,(rv._atom,lower,upper))


def _independent(variables):
    return all(v.group is v._atom and v.process is None for v in variables)


def _unevaluated_expectation(expr,condition=None):
    """``Expectation(expr, condition).doit(evaluate=False)`` of the reference:
    the exact integral or sum over the density, unevaluated."""
    expr = _native(expr)
    variables = _variables(expr)
    if not variables:return expr
    if condition is None:
        if len(variables) == 1 and _independent(variables) and variables[0].distribution.kind != 'finite':
            rv = variables[0];model = rv.distribution
            return _binder(rv,expr*model.pdf(rv._atom),model.lower,model.upper)
        return Expectation(expr)
    from ._queries import random_symbols,_condition_set
    condition = _native(condition)
    found = random_symbols(condition)|set(variables)
    if len(found) == 1 and _independent(found) and variables[0].distribution.kind == 'continuous':
        rv = variables[0];model = rv.distribution
        try:region = _condition_set(condition,rv)
        except NotImplementedError:region = None
        if isinstance(region,S.Interval):
            y = S.Dummy(rv.name,real=True)
            norm = S.Integral(model.pdf(y),(y,region.start,region.end))
            return S.Integral(expr*model.pdf(rv._atom)/norm,(rv._atom,region.start,region.end))
    return Expectation(expr,condition)


class StatisticalExpression:
    operation = None
    def __init__(self,*args,**kwargs):
        self.inputs = args;self.options = kwargs
        native = [S.MatrixExpr(arg) if isinstance(arg,S.Matrix) else _native(arg) for arg in args if arg is not None]
        self._expr = S._core._apply(type(self).__name__,*native,evaluate=False)
        HELD[self._expr] = self
    @property
    def args(self):return tuple(_native(arg) for arg in self.inputs if arg is not None)
    @property
    def func(self):return type(self)
    @property
    def free_symbols(self):return self._expr.free_symbols
    def _symox_(self):return self._expr
    def __repr__(self):
        import re
        from ._queries import random_symbols
        text = str(self._expr)
        # random symbols print by their names, not as the native dummies
        for v in random_symbols(self._expr):
            text = re.sub(r'(?<![\w])_%s(?![\w])' % re.escape(v.name),v.name.replace('\\','\\\\'),text)
        return text
    __str__ = __repr__
    def __hash__(self):return hash(self._expr)
    def __eq__(self,other):
        if isinstance(other,StatisticalExpression):return self._expr == other._expr
        return isinstance(other,S.Basic) and self._expr == other
    @property
    def _condition(self):return self.inputs[-1] if self.inputs else None
    def _inputs_doit(self,**hints):
        return tuple(arg.doit(**hints) if isinstance(arg,S.Basic) and _held_nodes(arg) else arg for arg in self.inputs)
    def doit(self,**hints):
        from . import _queries
        inputs = self._inputs_doit(**hints) if hints.get('deep',True) else self.inputs
        return getattr(_queries,self.operation)(*inputs,**self.options)
    def evaluate_integral(self):return self.rewrite('Integral').doit()
    evaluate_sum = evaluate_integral
    def subs(self,*args,**kwargs):
        inputs = tuple((arg if isinstance(arg,S.Matrix) else scalar(arg)).subs(*args,**kwargs) if arg is not None else None for arg in self.inputs)
        return type(self)(*inputs,**self.options)
    def xreplace(self,mapping):return _held_replace(self._expr,mapping)
    def has(self,*patterns):return self._expr.has(*(_native(p) if isinstance(p,StatisticalExpression) else p for p in patterns))
    def expand(self,**kwargs):return self
    def rewrite(self,*targets,**kwargs):
        result = self
        for target in targets:result = rewrite_held(result,target)
        return HELD.get(result,result) if isinstance(result,S.Basic) else result
    def _rewrite(self,name):return self
    def simplify(self,**kwargs):return self
    def __add__(self,other):return self._expr+_native_other(other)
    def __radd__(self,other):return _native_other(other)+self._expr
    def __sub__(self,other):return self._expr-_native_other(other)
    def __rsub__(self,other):return _native_other(other)-self._expr
    def __mul__(self,other):return self._expr*_native_other(other)
    def __rmul__(self,other):return _native_other(other)*self._expr
    def __truediv__(self,other):return self._expr/_native_other(other)
    def __rtruediv__(self,other):return _native_other(other)/self._expr
    def __pow__(self,other):return self._expr**_native_other(other)
    def __neg__(self):return -self._expr


def _native_other(other):
    return other._expr if isinstance(other,StatisticalExpression) else other


class Probability(StatisticalExpression):
    operation = 'P'
    def __init__(self,prob,condition=None,**kwargs):super().__init__(prob,condition,**kwargs)
    def _rewrite(self,name):
        if name not in ('Integral','Sum'):return self
        condition,given_condition = self.inputs
        from ._queries import random_symbols,_condition_set,P
        event = _native(condition)
        found = random_symbols(event)
        if given_condition is None and len(found) == 1 and _independent(found):
            rv = next(iter(found));model = rv.distribution
            if model.kind == 'continuous':
                try:region = _condition_set(event,rv)
                except NotImplementedError:region = None
                if isinstance(region,S.Interval):
                    z = S.Dummy('z',real=True)
                    return S.Integral(model.pdf(z),(z,region.start,region.end))
        return P(condition,given_condition)


class Expectation(StatisticalExpression):
    operation = 'E'
    def __new__(cls,expr,condition=None,**kwargs):
        if isinstance(expr,S.Matrix):return ExpectationMatrix(expr,condition)
        if condition is None and not _is_random(expr):return _native(expr)
        return super().__new__(cls)
    def __init__(self,expr,condition=None,**kwargs):super().__init__(expr,condition,**kwargs)
    def expand(self,**kwargs):
        expr,condition = self.inputs;expr = _native(expr)
        if not _is_random(expr):return expr
        if expr.is_Add:return S.Add(*(_native(_expand_held(Expectation(part,condition))) for part in expr.args))
        expanded = S.expand(expr)
        if expanded.is_Add:return S.Add(*(_native(_expand_held(Expectation(part,condition))) for part in expanded.args))
        if expr.is_Mul:
            rv = [a for a in expr.args if _is_random(a)];nonrv = [a for a in expr.args if not _is_random(a)]
            return S.Mul(*nonrv)*_native(Expectation(S.Mul(*rv),condition))
        return self
    def _rewrite(self,name):
        expr,condition = self.inputs;expr = _native(expr)
        if name in ('Integral','Sum'):return _unevaluated_expectation(expr,condition)
        if name != 'Probability':return self
        variables = _variables(expr)
        if len(variables) > 1:raise NotImplementedError('rewriting an expectation of several random variables as a probability')
        if not variables:return expr
        rv = variables[0];model = rv.distribution
        if model.kind == 'finite':raise NotImplementedError('rewriting a finite expectation as a probability')
        symbol = S.Symbol(rv.name.lower()) if rv.name[0].isupper() else S.Symbol(rv.name+'_1')
        body = _held_replace(expr,{rv._atom:symbol})*_native(Probability(S.Eq(rv._atom,symbol),condition))
        limits = (symbol,model.set.inf,model.set.sup)
        return S.Integral(body,limits) if model.kind == 'continuous' else S.Sum(body,limits)


def _expand_held(value):
    return value.expand() if isinstance(value,StatisticalExpression) else value


class Variance(StatisticalExpression):
    operation = 'variance'
    def __init__(self,arg,condition=None,**kwargs):super().__init__(arg,condition,**kwargs)
    def expand(self,**kwargs):
        from itertools import combinations
        from ._queries import random_symbols
        arg,condition = self.inputs;arg = _native(arg)
        if not _is_random(arg):return S.Integer(0)
        if arg in {v._atom for v in random_symbols(arg)}:return self
        if arg.is_Add:
            rv = [a for a in arg.args if _is_random(a)]
            variances = S.Add(*(_native(_expand_held(Variance(x,condition))) for x in rv))
            covariances = S.Add(*(2*_native(_expand_held(Covariance(x,y,condition))) for x,y in combinations(rv,2)))
            return variances+covariances
        if arg.is_Mul:
            rv = [a for a in arg.args if _is_random(a)];nonrv = [a**2 for a in arg.args if not _is_random(a)]
            if not rv:return S.Integer(0)
            return S.Mul(*nonrv)*_native(Variance(S.Mul(*rv),condition))
        return self
    def _rewrite(self,name):
        X,condition = self.inputs;X = _native(X)
        if name == 'Expectation':return Expectation(X**2,condition)-Expectation(X,condition)**2
        if name == 'Probability':return rewrite_held(rewrite_held(self,'Expectation'),'Probability')
        if name in ('Integral','Sum'):
            # variance(X, condition, evaluate=False) of the reference
            return rewrite_held(Expectation((X-_native(Expectation(X,condition)))**2,condition),'Integral')
        return self


def _sort_key(expr):
    from ._queries import random_symbols
    value = _native(expr)
    return S.default_sort_key(value.xreplace({v._atom:v.symbol for v in random_symbols(value)}))


def _coefficient_rv_pairs(expr):
    def split(m):
        return S.Mul(*(a for a in m.args if not _is_random(a))),S.Mul(*(a for a in m.args if _is_random(a)))
    if expr.is_Add:
        return [split(a) if a.is_Mul else (S.Integer(1),a) for a in expr.args if _is_random(a)]
    if expr.is_Mul:return [split(expr)]
    return [(S.Integer(1),expr)] if _is_random(expr) else []


class Covariance(StatisticalExpression):
    operation = 'covariance'
    def __init__(self,arg1,arg2,condition=None,**kwargs):
        if kwargs.pop('evaluate',True) and not isinstance(arg1,S.Matrix) and not isinstance(arg2,S.Matrix):
            arg1,arg2 = sorted([arg1,arg2],key=_sort_key)
        super().__init__(arg1,arg2,condition,**kwargs)
    def expand(self,**kwargs):
        from ._queries import random_symbols
        arg1,arg2,condition = self.inputs;arg1,arg2 = _native(arg1),_native(arg2)
        if arg1 == arg2:return _expand_held(Variance(arg1,condition))
        if not _is_random(arg1) or not _is_random(arg2):return S.Integer(0)
        atoms = {v._atom for v in random_symbols(arg1)|random_symbols(arg2)}
        if arg1 in atoms and arg2 in atoms:return self
        left,right = _coefficient_rv_pairs(S.expand(arg1)),_coefficient_rv_pairs(S.expand(arg2))
        return S.Add(*(a*b*_native(Covariance(r1,r2,condition)) for a,r1 in left for b,r2 in right))
    def _rewrite(self,name):
        X,Y,condition = self.inputs;X,Y = _native(X),_native(Y)
        if name == 'Expectation':return Expectation(X*Y,condition)-Expectation(X,condition)*Expectation(Y,condition)
        if name == 'Probability':return rewrite_held(rewrite_held(self,'Expectation'),'Probability')
        if name in ('Integral','Sum'):
            # covariance(X, Y, condition, evaluate=False) of the reference
            return _unevaluated_or_held((X-_native(Expectation(X,condition)))*(Y-_native(Expectation(Y,condition))),condition)
        return self


def _unevaluated_or_held(expr,condition):
    variables = _variables(expr)
    if len(variables) == 1:return rewrite_held(Expectation(expr,condition),'Integral')
    return Expectation(expr,condition)


class Moment(StatisticalExpression):
    operation = 'moment'
    def __init__(self,X,n,c=0,condition=None,**kwargs):super().__init__(X,n,c,condition,**kwargs)
    def doit(self,**hints):return rewrite_held(self,'Expectation').doit(**hints)
    def _rewrite(self,name):
        X,n,c,condition = self.inputs
        value = Expectation((_native(X)-_native(c))**_native(n),condition)
        return value if name == 'Expectation' else rewrite_held(value,name)


class CentralMoment(StatisticalExpression):
    operation = 'cmoment'
    def __init__(self,X,n,condition=None,**kwargs):super().__init__(X,n,condition,**kwargs)
    def doit(self,**hints):return rewrite_held(self,'Expectation').doit(**hints)
    def _rewrite(self,name):
        X,n,condition = self.inputs
        mu = Expectation(X,condition)
        value = Moment(X,n,_native(mu),condition)._rewrite('Expectation')
        return value if name == 'Expectation' else rewrite_held(value,name)


def _matrix_value(value):
    from ._matrix import MATRIX_REGISTRY
    try:return MATRIX_REGISTRY.get(value,value)
    except TypeError:return value


class ExpectationMatrix(StatisticalExpression):
    def __init__(self,expr,condition=None):super().__init__(expr,condition)
    def doit(self,**hints):
        from ._queries import E
        expr,condition = self.inputs
        return E(_matrix_value(expr),condition)


class VarianceMatrix(StatisticalExpression):
    def __init__(self,arg,condition=None):super().__init__(arg,condition)
    @property
    def shape(self):
        value = _matrix_value(self.inputs[0]);dimension = value.rows if value.cols == 1 else value.cols
        return dimension,dimension
    def doit(self,**hints):
        expr,condition = self.inputs;expr = _matrix_value(expr)
        if getattr(expr,'_stats_matrix',False):return expr.variance()
        from ._queries import covariance
        if not isinstance(expr,S.Matrix):raise TypeError('matrix variance requires a random matrix or an explicit matrix of random scalars')
        if expr.cols == 1:
            return S.Matrix([[covariance(expr[i,0],expr[j,0],condition) for j in range(expr.rows)] for i in range(expr.rows)])
        return S.Matrix([[S.Add(*(covariance(expr[i,j],expr[i,k],condition) for i in range(expr.rows))) for k in range(expr.cols)] for j in range(expr.cols)])


class CrossCovarianceMatrix(StatisticalExpression):
    def __init__(self,arg1,arg2,condition=None):super().__init__(arg1,arg2,condition)
    @property
    def shape(self):
        left,right = map(_matrix_value,self.inputs[:2])
        return (left.rows,right.rows) if left.cols == right.cols == 1 else (left.cols,right.cols)
    def doit(self,**hints):
        left,right,condition = self.inputs;left,right = _matrix_value(left),_matrix_value(right)
        column_vectors = left.cols == right.cols == 1
        if not column_vectors and left.rows != right.rows:raise ValueError('cross covariance requires equal row counts')
        if left is right:return VarianceMatrix(left,condition).doit()
        if getattr(left,'_stats_matrix',False) and getattr(right,'_stats_matrix',False):return S.zeros(*self.shape)
        from ._queries import covariance
        if column_vectors:return S.Matrix([[covariance(left[i,0],right[j,0],condition) for j in range(right.rows)] for i in range(left.rows)])
        return S.Matrix([[S.Add(*(covariance(left[i,j],right[i,k],condition) for i in range(left.rows))) for k in range(right.cols)] for j in range(left.cols)])


def install():
    if getattr(S.Expr,'_stats_doit_installed',False):return
    original = S.Expr.doit
    @wraps(original)
    def doit(self,*args,**kwargs):
        if self in HELD:return HELD[self].doit(**kwargs)
        if kwargs.get('deep',True):
            replacements,seen,stack = {},set(),[self]
            while stack:
                node = stack.pop()
                if not isinstance(node,S.Expr) or node in seen:continue
                seen.add(node)
                if node in HELD:replacements[node] = HELD[node].doit(**kwargs)
                else:stack.extend(getattr(node,'args',()))
            if replacements:return original(self.xreplace(replacements),*args,**kwargs)
        return original(self,*args,**kwargs)
    S.Expr.doit = doit;S.Expr._stats_doit_installed = True
    original_rewrite = S.Expr.rewrite
    @wraps(original_rewrite)
    def rewrite(self,*args,**kwargs):
        targets = [a for a in args if _target_name(a) is not None]
        if targets and len(targets) == len(args) and _held_nodes(self):
            result = self
            for target in targets:result = rewrite_held(result,target)
            return result
        return original_rewrite(self,*args,**kwargs)
    S.Expr.rewrite = rewrite


__all__ = ['Probability','Expectation','Variance','Covariance','Moment','CentralMoment',
 'ExpectationMatrix','VarianceMatrix','CrossCovarianceMatrix']
