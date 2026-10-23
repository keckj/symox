"""Probability models over native symbolic scalars and exact domains."""
from collections.abc import Mapping
from itertools import product
from math import prod
import operator
import numpy as np
import symox as S

REGISTRY = {}
_RV_CACHE = {}
HALF = S.Rational(1,2)


def scalar(value):
    if isinstance(value,np.generic):value = value.item()
    return S.sympify(value)


def check(condition,message):
    if condition is False or condition == S.false:raise ValueError(message)


def positive(*values):
    for value in values:
        check(value > 0,'distribution parameters must be positive')
        check(value.is_real is not False,'distribution parameters must be real')


def integral(expr,var,lower,upper):
    return S.Integral(expr,(var,lower,upper))


def sum_expr(expr,var,lower,upper):
    if lower.is_Integer and upper.is_Integer and 0 <= int(upper-lower) <= 4096:
        return S.Add(*(expr.subs(var,j) for j in range(int(lower),int(upper)+1)))
    return S.Sum(expr,(var,lower,upper))


class Distribution:
    kind = 'continuous'
    family = 'Handmade'
    _argnames = ()
    def __init__(self,*args):
        self.args = tuple(scalar(a) for a in args)
        self.lower,self.upper = -S.oo,S.oo
        self._support = None
        self._table = None
        self._pdf = None
        self._cdf = None
        # The complete reference-form CDF (support Piecewise included), when
        # SymPy writes one; ``_cdf`` stays the in-support formula.
        self._cdf_form = None
        self._cdf_integrated = None
        self._moment = None
        self._mgf = None
        self._quantile = None
        self._sampler = None
        self._entropy = None
    def __repr__(self):return '%sDistribution(%s)' % (self.family,', '.join(map(str,self.args)))
    def __getattr__(self,name):
        if name in self._argnames:return self.args[self._argnames.index(name)]
        raise AttributeError(name)
    @property
    def is_Finite(self):return self.kind == 'finite'
    @property
    def is_Discrete(self):return self.kind in ('finite','discrete')
    @property
    def is_Continuous(self):return self.kind == 'continuous'
    @property
    def free_symbols(self):
        def collect(value):
            if isinstance(value,(list,tuple)):return set().union(*(collect(a) for a in value))
            return getattr(value,'free_symbols',set())
        if self.kind == 'finite' and self._table is not None:
            return set().union(*(collect(p) for p in self.dict.values()))
        return set().union(*(collect(a) for a in self.args))-getattr(self,'_bound_symbols',set())
    @property
    def set(self):
        if self._support is not None:return self._support
        if getattr(self,'_held_density',False) and self._table is None:
            if self.upper.is_Integer:return S.FiniteSet(*S.Range(self.lower,self.upper+1))
            return S.Range(self.lower,self.upper+1)
        if self.kind == 'finite' and self._table is not None:return S.FiniteSet(*self.dict)
        if self.kind == 'discrete':return S.Range(self.lower,self.upper+1)
        return S.Interval(self.lower,self.upper)
    @property
    def dict(self):
        if self._table is None and getattr(self,'_held_density',False):return HeldDensity(self)
        if self._table is None:raise ValueError('the support is not an explicit finite table')
        return self._table() if callable(self._table) else dict(self._table)
    def __call__(self,*args):return self.pdf(*args)
    def remake(self,function):
        """The same family with ``function`` applied to every parameter, or
        None when the parameters do not change or cannot be rebuilt."""
        params = getattr(self,'_params',None)
        if params is None:return None
        def apply(value):
            if isinstance(value,(list,tuple)):return type(value)(apply(v) for v in value)
            return function(value) if hasattr(value,'subs') else value
        values = tuple(apply(v) for v in params)
        if values == params:return None
        return self._remake(values)
    def pdf(self,x):
        x = scalar(x)
        if self._pdf is not None:return self._pdf(x)
        if self.kind == 'finite':return S.Piecewise(*((p,S.Eq(x,value)) for value,p in self.dict.items()),(0,True))
        raise ValueError('a density must be specified for a distribution')
    pmf = pdf
    def bounded_pdf(self,x):
        value = self.pdf(x)
        if self.kind == 'finite':return value
        condition = S.And(x >= self.lower,x <= self.upper)
        if self.kind == 'discrete':condition = S.And(condition,S.Contains(x,S.Integers))
        return S.Piecewise((value,condition),(0,True))
    def cdf(self,x,**kwargs):
        x = scalar(x)
        if x == -S.oo:return S.Integer(0)
        if x == S.oo:return S.Integer(1)
        if self.kind == 'finite':
            return S.Add(*(S.Piecewise((p,x >= value),(0,True)) for value,p in self.dict.items()))
        # SymPy integrates the density whenever integration options are
        # given, and otherwise returns the family's own CDF expression.
        if kwargs and self._cdf_integrated is not None:return self._cdf_integrated(x)
        if self._cdf_form is not None:return self._cdf_form(x)
        t = S.Dummy('t',real=True)
        if self._cdf is not None:value = self._cdf(x)
        elif self.kind == 'discrete':value = sum_expr(self.pdf(t),t,self.lower,S.floor(x))
        else:value = integral(self.pdf(t),t,self.lower,x)
        return S.Piecewise((0,x < self.lower),(1,x >= self.upper),(value,True))
    def moment(self,n):
        n = scalar(n)
        if n == 0:return S.Integer(1)
        if self._moment is not None:
            result = self._moment(n)
            if result is not None:return scalar(result)
        if self.kind == 'finite':return S.Add(*(v**n*p for v,p in self.dict.items()))
        x = S.Dummy('x',integer=True) if self.kind == 'discrete' else S.Dummy('x',real=True)
        expression = x**n*self.pdf(x)
        return sum_expr(expression,x,self.lower,self.upper) if self.kind == 'discrete' else integral(expression,x,self.lower,self.upper)
    def mgf(self,t):
        t = scalar(t)
        if t == 0:return S.Integer(1)
        if self._mgf is not None:return self._mgf(t)
        if self.kind == 'finite':return S.Add(*(p*S.exp(t*v) for v,p in self.dict.items()))
        x = S.Dummy('x',real=True)
        integrand = S.exp(t*x)*self.pdf(x)
        return sum_expr(integrand,x,self.lower,self.upper) if self.kind == 'discrete' else integral(integrand,x,self.lower,self.upper)
    def quantile(self,p):
        p = scalar(p)
        if p.is_number:check(S.And(p >= 0,p <= 1),'quantile probability must lie in [0,1]')
        if self._quantile is not None:return self._quantile(p)
        if self.kind == 'finite':
            total,pieces = S.Integer(0),[]
            for value,probability in sorted(self.dict.items(),key=lambda pair:float(pair[0])):
                total += probability
                pieces.append((value,p <= total))
            return S.Piecewise((S.nan,S.Or(p < 0,p > 1)),*pieces)
        x = S.Dummy('q',real=True)
        roots = S.solve(S.Eq(self.cdf(x),p),x)
        if len(roots) == 1:return roots[0]
        raise NotImplementedError('this distribution has no implemented symbolic quantile')
    def sample(self,rng,size=()):
        if self.free_symbols:raise ValueError('sampling requires numeric distribution parameters')
        if self._sampler is not None:return self._sampler(rng,size)
        if self.kind == 'finite':
            table = self.dict
            values = list(table)
            indices = rng.choice(len(values),size=size,p=np.array([float(table[v]) for v in values]))
            def convert(value):
                try:return int(value) if value.is_Integer else float(value)
                except (ValueError,TypeError):return value
            result = np.asarray([convert(values[int(i)]) for i in np.asarray(indices).flat]).reshape(size)
            return result.item() if not size else result
        return self._inverse_samples(rng,size)
    def _inverse_samples(self,rng,size):
        # Numeric inverse transform. A closed CDF/quantile is preferred; the
        # continuous fallback integrates the actual density over its support.
        import mpmath
        mp = mpmath.mp.clone();mp.dps = 25
        x = S.Dummy('sample_x',real=True)
        expression = self._cdf(x) if self._cdf is not None else None
        density = self.pdf(x)
        def evaluate(expr,value):
            result = expr.subs(x,S.Float(str(value),25)).evalf(25)
            return mp.mpf(str(result))
        def cumulative(value):
            if expression is not None:return evaluate(expression,value)
            return mp.quad(lambda t:evaluate(density,t),[float(self.lower),value])
        def draw(p):
            if self._quantile is not None:return float(self._quantile(S.Float(float(p))).evalf())
            if self.kind == 'discrete':
                if self.lower == -S.oo:raise NotImplementedError('two-sided discrete inverse sampling needs a specialized sampler')
                total = 0.0;index = int(self.lower)
                while total < p:
                    total += float(self.pdf(S.Integer(index)).evalf());index += 1
                    if index > 1000000:raise RuntimeError('inverse sampling exceeded its work limit')
                return index-1
            lo,hi = float(self.lower),float(self.upper)
            if not np.isfinite(lo):
                lo = -1.0
                for i in range(1024):
                    if cumulative(lo) <= p:break
                    lo *= 2
                else:raise RuntimeError('failed to bracket a quantile')
            if not np.isfinite(hi):
                hi = max(1.0,lo+1)
                for i in range(1024):
                    if cumulative(hi) >= p:break
                    hi = max(2*hi,hi+1)
                else:raise RuntimeError('failed to bracket a quantile')
            for i in range(80):
                midpoint = (lo+hi)/2
                if cumulative(midpoint) < p:lo = midpoint
                else:hi = midpoint
            return (lo+hi)/2
        uniforms = np.asarray(rng.random(size))
        result = np.array([draw(float(p)) for p in uniforms.flat]).reshape(size)
        return result.item() if not size else result


class HeldDensity:
    """The unevaluated density table of a distribution with symbolic size,
    printed as ``Density(...)``; ``subs`` then ``doit`` evaluates it."""
    def __init__(self,distribution):self.distribution = distribution
    def __repr__(self):return 'Density(%r)' % (self.distribution,)
    __str__ = __repr__
    @property
    def args(self):return (self.distribution,)
    @property
    def free_symbols(self):return self.distribution.free_symbols
    def subs(self,*args,**kwargs):
        model = self.distribution.remake(lambda v:v.subs(*args,**kwargs))
        return self if model is None else HeldDensity(model)
    def doit(self,**hints):
        return self if self.distribution._table is None else self.distribution.dict
    def _table(self,*args):raise NotImplementedError('a distribution with symbolic size has no explicit density table')
    items = keys = values = __iter__ = __len__ = __getitem__ = _table
    def __eq__(self,other):return isinstance(other,HeldDensity) and repr(self) == repr(other)
    def __hash__(self):return hash(repr(self))


class RandomAtom(S.Dummy):
    """The native atom of a random symbol in arithmetic: an independent
    Dummy that prints as the plain name, as SymPy prints a RandomSymbol."""
    _symox_user_hooks = True
    _symox_nominal = 'Dummy'
    def _sympystr(self,printer):return self.name
    def _latex(self,printer):return printer._print(S.Symbol(self.name))
    def _pretty(self,printer):return printer._print(S.Symbol(self.name))
    # The random symbol of this atom (a Tuple entry or an arithmetic
    # operand) keeps its probability space, as SymPy's RandomSymbol.
    pspace = property(lambda self:REGISTRY[self].pspace)
    distribution = property(lambda self:REGISTRY[self].distribution)


class RandomSymbol:
    """A probability-space scalar coercible to an independent native Dummy."""
    is_random = True
    is_Symbol = True
    is_Atom = True
    is_commutative = True
    def __init__(self,name,distribution,group=None):
        self.name = str(name)
        self.distribution = distribution
        self.symbol = S.Symbol(self.name,real=True) if distribution.kind != 'finite' else S.Symbol(self.name)
        assumptions = {'real':True} if distribution.kind != 'finite' or distribution._table is None or all(v.is_real is True for v in distribution.dict) else {}
        if distribution.kind == 'discrete':assumptions['integer'] = True
        self._atom = RandomAtom(self.name,**assumptions)
        self.group = self._atom if group is None else group
        self.process = None
        REGISTRY[self._atom] = self
    def _symox_(self):return self._atom
    def __repr__(self):return self.name
    def __str__(self):return self.name
    def __hash__(self):return hash(self._atom)
    def __eq__(self,other):return isinstance(other,RandomSymbol) and self._atom == other._atom
    @property
    def args(self):return self.symbol,self.pspace
    @property
    def free_symbols(self):return {self}
    @property
    def pspace(self):
        from ._spaces import make_space
        return make_space((self,))
    @property
    def is_Finite(self):return self.distribution.is_Finite
    @property
    def is_Continuous(self):return self.distribution.is_Continuous
    @property
    def is_Discrete(self):return self.distribution.is_Discrete
    def __getattr__(self,name):return getattr(self._atom,name)
    def subs(self,*args,**kwargs):
        """Substituting into the parameters gives the variable of the same
        name with the new distribution, as in SymPy ``X.subs(mu, 2)``."""
        if self.group is self._atom and self.process is None:
            model = self.distribution.remake(lambda v:v.subs(*args,**kwargs))
            if model is not None:return variable(self.name,model)
        return self._atom.subs(*args,**kwargs)
    def __add__(self,other):return self._atom+other
    def __radd__(self,other):return other+self._atom
    def __sub__(self,other):return self._atom-other
    def __rsub__(self,other):return other-self._atom
    def __mul__(self,other):return self._atom*other
    def __rmul__(self,other):return other*self._atom
    def __truediv__(self,other):return self._atom/other
    def __rtruediv__(self,other):return other/self._atom
    def __pow__(self,other):return self._atom**other
    def __rpow__(self,other):return other**self._atom
    def __neg__(self):return -self._atom
    def __abs__(self):return S.Abs(self._atom)
    def __lt__(self,other):return self._atom < other
    def __le__(self,other):return self._atom <= other
    def __gt__(self,other):return self._atom > other
    def __ge__(self,other):return self._atom >= other


def variable(name,distribution,group=None):
    if group is None:
        try:key = str(name),distribution.family,distribution.args,tuple(distribution.dict.items()) if distribution.kind == 'finite' and distribution._table is not None else None
        except (ValueError,TypeError):key = None
        if key is not None and key in _RV_CACHE:return _RV_CACHE[key]
    else:key = None
    result = RandomSymbol(name,distribution,group)
    if key is not None:_RV_CACHE[key] = result
    return result


from ._domains import RandomDomain


class Domain(RandomDomain):
    def __init__(self,variables,set=None,condition=None):
        self.symbols = frozenset(v.symbol for v in variables)
        self.variables = tuple(variables)
        self.set = set if set is not None else (variables[0].distribution.set if len(variables) == 1 else S.ProductSet(*(v.distribution.set for v in variables)))
        self.condition = condition
    def __repr__(self):return 'Domain(%s, %s)' % (', '.join(map(str,self.symbols)),self.set)
    def as_boolean(self):
        return self.condition if self.condition is not None else S.And(*(S.Contains(v._atom,v.distribution.set) for v in self.variables))


class PSpace:
    def __init__(self,variables):
        self.values = frozenset(variables)
        self.variables = tuple(sorted(variables,key=lambda v:v.name))
        self.domain = Domain(self.variables)
    @property
    def distribution(self):
        if len(self.variables) != 1:raise AttributeError('a product space has several distributions')
        return self.variables[0].distribution
    @property
    def symbol(self):
        if len(self.variables) != 1:raise AttributeError('a product space has several symbols')
        return self.variables[0].symbol
    def probability(self,condition,**kwargs):
        from ._queries import P
        return P(condition,**kwargs)
    def compute_expectation(self,expr,**kwargs):
        from ._queries import E
        return E(expr,**kwargs)
    def sample(self,size=(),library='numpy',seed=None):
        if library not in ('numpy','scipy'):raise NotImplementedError('sampling supports the NumPy/scipy compatibility selectors')
        if isinstance(size,int):size = (size,)
        size = tuple(operator.index(n) for n in size)
        if any(n<0 for n in size):raise ValueError('sample dimensions must be nonnegative')
        rng = np.random.default_rng(seed)
        from ._queries import _groups
        result = {}
        for group,variables in _groups(self.variables).items():
            if getattr(group,'_stats_joint',False):
                draw = group.sample(rng,size)
                result.update({v:draw[...,v.index] for v in variables})
            elif getattr(group,'_stats_matrix',False):
                draw = group.sample(rng,size)
                result.update({v:draw if v is group else draw[...,v.matrix_index[0],v.matrix_index[1]] for v in variables})
            elif hasattr(group,'sample_mapping'):
                count = int(np.prod(size)) if size else 1
                mappings = [group.sample_mapping(rng,variables) for _ in range(count)]
                for v in variables:
                    values = np.asarray([m[v._atom] for m in mappings],dtype=object).reshape(size)
                    result[v] = values.item() if not size else values
            else:
                for v in variables:result[v] = v.distribution.sample(rng,size)
        return result
    def __eq__(self,other):return isinstance(other,PSpace) and self.values == other.values
    def __hash__(self):return hash(self.values)
