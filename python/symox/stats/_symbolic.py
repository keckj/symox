"""Held statistical operators with executable exact counterparts."""
from functools import wraps
import symox as S
from ._base import scalar

HELD = {}


def _native(arg):
    if isinstance(arg,StatisticalExpression):return arg._expr
    return arg if isinstance(arg,S.DenseMatrix) else scalar(arg)


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
        native = [S.MatrixExpr(arg) if isinstance(arg,S.DenseMatrix) else _native(arg) for arg in args if arg is not None]
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
        return tuple(_input_doit(arg,**hints) for arg in self.inputs)
    def doit(self,**hints):
        from . import _queries
        inputs = self._inputs_doit(**hints) if hints.get('deep',True) else self.inputs
        return getattr(_queries,self.operation)(*inputs,**self.options)
    def evaluate_integral(self):return self.rewrite('Integral').doit()
    evaluate_sum = evaluate_integral
    def subs(self,*args,**kwargs):
        inputs = tuple((arg if isinstance(arg,S.DenseMatrix) else scalar(arg)).subs(*args,**kwargs) if arg is not None else None for arg in self.inputs)
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


def _input_doit(arg,**hints):
    # Basic.doit(deep=True) of one argument: nested held operators evaluate.
    if isinstance(arg,StatisticalExpression):return arg.doit(**hints)
    if isinstance(arg,S.Basic) and _held_nodes(arg):return arg.doit(**hints)
    return arg


def _held_expectation(value):
    held = HELD.get(value) if isinstance(value,S.Basic) else value
    return isinstance(held,Expectation)


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
        if isinstance(expr,S.DenseMatrix):return ExpectationMatrix(expr,condition)
        if condition is None and not _is_random(expr):return _native(expr)
        return super().__new__(cls)
    def __init__(self,expr,condition=None,**kwargs):super().__init__(expr,condition,**kwargs)
    def doit(self,**hints):
        # Reference Expectation.doit: nested expectations are kept with
        # deep=False, and linearity skips the held expectation terms.
        expr,condition = self.inputs
        if hints.get('numsamples') or not hints.get('evaluate',True) or self.options:
            return super().doit(**hints)
        if hints.get('deep',True):expr = _input_doit(expr,**hints)
        native = _native(expr)
        if not _is_random(native) or _held_expectation(native):return native
        if not isinstance(native,S.Basic):
            from ._queries import E
            return E(expr,condition)
        if condition is None and native.is_Add:
            return S.Add(*(_native(Expectation(arg)) if _held_expectation(arg) else _native(_input_doit(Expectation(arg),**hints)) for arg in native.args))
        if condition is None and native.is_Mul and any(_held_expectation(node) for node in _held_nodes(native)):return native
        from ._queries import E
        return E(expr,condition)
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
        if kwargs.pop('evaluate',True) and not isinstance(arg1,S.DenseMatrix) and not isinstance(arg2,S.DenseMatrix):
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


def _mnative(value):
    """The native form of a matrix operand: held matrix statistics and
    random matrix symbols become their native matrix symbols."""
    if isinstance(value,_MatrixStatistic):return value._symox_()
    if getattr(value,'_stats_matrix',False):return value.symbol
    if isinstance(value,S.DenseMatrix):return value
    return S.sympify(value)


def _mheld(value):
    """The held matrix statistic or random matrix symbol a native operand stands for."""
    from ._matrix import MATRIX_REGISTRY
    try:
        if value in _MATRIX_HELD:return _MATRIX_HELD[value]
        return MATRIX_REGISTRY.get(value,value)
    except TypeError:return value


def _is_random_matrix(value):
    value = _mheld(value)
    if getattr(value,'_stats_matrix',False):return True
    if isinstance(value,_MatrixStatistic):return _is_random_matrix(value.inputs[0])
    from ._queries import random_symbols
    try:return bool(random_symbols(_mnative(value)))
    except Exception:return False


def _is_random_matrix_symbol(value):
    return getattr(_mheld(value),'_stats_matrix',False)


def _shape(value):
    value = _mheld(value)
    return tuple(value.shape)


def _product(values):
    result = S.Integer(1)
    for value in values:result = result*_mnative(value)
    return result


def _sum(values):
    values = [_mnative(v) for v in values]
    if not values:return S.Integer(0)
    result = values[0]
    for value in values[1:]:result = result+value
    return result


def _transpose(value):
    value = _mnative(value)
    return value.T if hasattr(value,'T') else value


def _is_head(value,name):
    return getattr(value,'is_'+name,False) is True


_MATRIX_HELD = {}


class _MatrixStatistic(StatisticalExpression):
    """Held matrix statistic: a matrix expression of known shape whose
    native stand-in is a matrix symbol printing as the held operator, so
    products keep their order and shape inside native matrix expressions."""
    _symox_matrixexpr = True
    is_Matrix = False
    def _native_symbol(self):
        symbol = getattr(self,'_symbol',None)
        if symbol is None:
            arguments = ', '.join(str(_mheld(a)) for a in self.inputs if a is not None)
            symbol = S.MatrixSymbol('%s(%s)'%(type(self).__name__,arguments),*self.shape)
            self._symbol = symbol;_MATRIX_HELD[symbol] = self
        return symbol
    def _symox_(self):return self._native_symbol()
    @property
    def rows(self):return self.shape[0]
    @property
    def cols(self):return self.shape[1]
    @property
    def T(self):return self._symox_().T
    def transpose(self):return self.T
    def __add__(self,other):return self._symox_()+_mnative(other)
    def __radd__(self,other):return _mnative(other)+self._symox_()
    def __sub__(self,other):return self._symox_()-_mnative(other)
    def __rsub__(self,other):return _mnative(other)-self._symox_()
    def __mul__(self,other):return self._symox_()*_mnative(other)
    def __rmul__(self,other):return _mnative(other)*self._symox_()
    def __matmul__(self,other):return self._symox_()*_mnative(other)
    def __rmatmul__(self,other):return _mnative(other)*self._symox_()
    def __neg__(self):return -self._symox_()


def _expectation_matrix(expr,condition=None):
    return ExpectationMatrix(expr,condition) if _is_random_matrix(expr) else _mnative(expr)


class ExpectationMatrix(_MatrixStatistic):
    def __new__(cls,expr,condition=None):
        if condition is None and not isinstance(expr,S.DenseMatrix) and not _is_random_matrix(expr):return _mnative(expr)
        return super().__new__(cls)
    def __init__(self,expr,condition=None):super().__init__(expr,condition)
    @property
    def shape(self):return _shape(self.inputs[0])
    def doit(self,**hints):
        from ._queries import E
        expr,condition = self.inputs
        return E(_matrix_value(expr),condition)
    def expand(self,**hints):
        # Reference ExpectationMatrix.expand
        expr,condition = self.inputs;expr = _mnative(expr)
        if not _is_random_matrix(expr):return expr
        if _is_head(expr,'MatAdd'):
            return _sum(_mnative(_expand_held(_expectation_matrix(a,condition))) for a in expr.args)
        expanded = expr.expand() if hasattr(expr,'expand') else expr
        if _is_head(expanded,'MatAdd'):
            return _sum(_mnative(_expand_held(_expectation_matrix(a,condition))) for a in expanded.args)
        if _is_head(expr,'MatMul'):
            rv,nonrv,postnon = [],[],[]
            for a in expr.args:
                if _is_random_matrix(a):
                    if rv:rv.extend(postnon)
                    else:nonrv.extend(postnon)
                    postnon = [];rv.append(a)
                elif getattr(a,'is_Matrix',False) or getattr(a,'_symox_matrixexpr',False):postnon.append(a)
                else:nonrv.append(a)
            if not nonrv:return self
            return _product(nonrv)*_mnative(_expectation_matrix(_product(rv),condition))*_product(postnon)
        return self


def _variance_matrix_expand(arg,condition):
    held = VarianceMatrix(arg,condition)
    return held.expand()


class VarianceMatrix(_MatrixStatistic):
    def __init__(self,arg,condition=None):
        if not isinstance(arg,S.DenseMatrix) and not getattr(_mheld(arg),'_stats_matrix',False):
            shape = _shape(arg)
            if 1 not in shape:
                from ..matrices import ShapeError
                raise ShapeError('Expression is not a vector')
        super().__init__(arg,condition)
    @property
    def shape(self):
        value = _mheld(self.inputs[0]);rows,cols = value.shape
        dimension = rows if cols == 1 else cols
        return dimension,dimension
    def doit(self,**hints):
        expr,condition = self.inputs;expr = _matrix_value(expr)
        if getattr(expr,'_stats_matrix',False):return expr.variance()
        from ._queries import covariance
        if not isinstance(expr,S.DenseMatrix):raise TypeError('matrix variance requires a random matrix or an explicit matrix of random scalars')
        if expr.cols == 1:
            return S.Matrix([[covariance(expr[i,0],expr[j,0],condition) for j in range(expr.rows)] for i in range(expr.rows)])
        return S.Matrix([[S.Add(*(covariance(expr[i,j],expr[i,k],condition) for i in range(expr.rows))) for k in range(expr.cols)] for j in range(expr.cols)])
    def expand(self,**hints):
        # Reference VarianceMatrix.expand
        from itertools import combinations
        arg,condition = self.inputs
        if isinstance(arg,S.DenseMatrix):return self
        arg = _mnative(arg)
        if not _is_random_matrix(arg):return S.ZeroMatrix(*self.shape)
        if _is_random_matrix_symbol(arg):return self
        if _is_head(arg,'MatAdd'):
            rv = [a for a in arg.args if _is_random_matrix(a)]
            variances = _sum(_mnative(_expand_held(VarianceMatrix(x,condition))) for x in rv)
            covariances = [2*_mnative(_expand_held(CrossCovarianceMatrix(x,y,condition))) for x,y in combinations(rv,2)]
            return _sum([variances]+covariances)
        if _is_head(arg,'MatMul'):
            nonrv = [a for a in arg.args if not _is_random_matrix(a)]
            rv = [a for a in arg.args if _is_random_matrix(a)]
            if not rv:return S.ZeroMatrix(*self.shape)
            if not nonrv or len(rv) > 1:return self
            coefficient = _product(nonrv)
            return coefficient*_mnative(VarianceMatrix(_product(rv),condition))*_transpose(coefficient)
        return self


def _mul_nonrv_rv(m):
    rv = [a for a in m.args if _is_random_matrix(a)]
    nonrv = [a for a in m.args if not _is_random_matrix(a)]
    return _product(nonrv),_product(rv)


def _expand_single_argument(expr):
    # CrossCovarianceMatrix._expand_single_argument of the reference
    if _is_random_matrix_symbol(expr):return [(S.Integer(1),expr)]
    if _is_head(expr,'MatAdd'):
        out = []
        for a in expr.args:
            if _is_head(a,'MatMul'):out.append(_mul_nonrv_rv(a))
            elif _is_random_matrix(a):out.append((S.Integer(1),a))
        return out
    if _is_head(expr,'MatMul'):return [_mul_nonrv_rv(expr)]
    if _is_random_matrix(expr):return [(S.Integer(1),expr)]
    return []


class CrossCovarianceMatrix(_MatrixStatistic):
    def __init__(self,arg1,arg2,condition=None):
        matrix_laws = all(getattr(_mheld(v),'_stats_matrix',False) for v in (arg1,arg2))
        if not isinstance(arg1,S.DenseMatrix) and not isinstance(arg2,S.DenseMatrix) and not matrix_laws:
            shape1,shape2 = _shape(arg1),_shape(arg2)
            if 1 not in shape1 or 1 not in shape2 or shape1[1] != shape2[1]:
                from ..matrices import ShapeError
                raise ShapeError('Expression is not a vector')
        super().__init__(arg1,arg2,condition)
    @property
    def shape(self):
        left,right = map(_mheld,self.inputs[:2])
        if isinstance(left,S.DenseMatrix) or isinstance(right,S.DenseMatrix):
            return (left.rows,right.rows) if left.cols == right.cols == 1 else (left.cols,right.cols)
        return (left.shape[0],right.shape[0]) if left.shape[1] == 1 and right.shape[1] == 1 else (1,1)
    def doit(self,**hints):
        left,right,condition = self.inputs;left,right = _matrix_value(left),_matrix_value(right)
        column_vectors = left.cols == right.cols == 1
        if not column_vectors and left.rows != right.rows:raise ValueError('cross covariance requires equal row counts')
        if left is right:return VarianceMatrix(left,condition).doit()
        if getattr(left,'_stats_matrix',False) and getattr(right,'_stats_matrix',False):return S.zeros(*self.shape)
        from ._queries import covariance
        if column_vectors:return S.Matrix([[covariance(left[i,0],right[j,0],condition) for j in range(right.rows)] for i in range(left.rows)])
        return S.Matrix([[S.Add(*(covariance(left[i,j],right[i,k],condition) for i in range(left.rows))) for k in range(right.cols)] for j in range(left.cols)])
    def expand(self,**hints):
        # Reference CrossCovarianceMatrix.expand
        arg1,arg2,condition = self.inputs
        if isinstance(arg1,S.DenseMatrix) or isinstance(arg2,S.DenseMatrix):return self
        arg1,arg2 = _mnative(arg1),_mnative(arg2)
        if arg1 == arg2:return _expand_held(VarianceMatrix(arg1,condition))
        if not _is_random_matrix(arg1) or not _is_random_matrix(arg2):return S.ZeroMatrix(*self.shape)
        if _is_random_matrix_symbol(arg1) and _is_random_matrix_symbol(arg2):return self
        left = _expand_single_argument(arg1.expand())
        right = _expand_single_argument(arg2.expand())
        return _sum(a*_mnative(CrossCovarianceMatrix(r1,r2,condition))*_transpose(b) for a,r1 in left for b,r2 in right)


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
