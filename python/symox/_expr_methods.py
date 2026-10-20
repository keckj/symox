"""Public expression methods composed from native symbolic operations."""
from . import _core

_native_has = _core.Expr.has
_native_coefficients_dict = _core.Expr.as_coefficients_dict
_native_ordered_terms = _core.Expr.as_ordered_terms
_native_ordered_factors = _core.Expr.as_ordered_factors
_native_coeff = _core.Expr.coeff
_native_match = _core.Expr.match


def coeff(self, x, n=1, right=False, _first=True):
    """Select explicit multiplicative factors without expanding or factoring."""
    from . import S, Add, Mul, sympify
    try:x=sympify(x)
    except (TypeError,ValueError):return S.Zero
    if not isinstance(x,_core.Expr):return S.Zero
    if isinstance(n,bool):raise ValueError('Coefficient power must be an integer')
    power=sympify(n)
    if getattr(power,'is_Integer',False) is not True:raise ValueError('Coefficient power must be an integer')
    n=int(power)
    if x==0:return S.Zero
    if x==self:return S.One if n==1 else S.Zero
    terms=self.args if self.is_Add else (self,)
    if x==1:
        return Add(*(term for term in terms if term.as_coeff_Mul()[0]==1))
    if n==0:
        if x.is_Add and self.is_Add:
            coefficient=coeff(self,x,right=right)
            if coefficient==0:return S.Zero
            pieces=coefficient.args if coefficient.is_Add else (coefficient,)
            return self-Add(*(part*x for part in pieces))
        return self.as_independent(x,as_Add=True)[0]
    if x.is_number is True and abs(n)>100000:
        raise NotImplementedError('Numeric coefficient powers are limited to magnitude 100000')
    target=x**n
    if self.is_commutative is not True or target.is_commutative is not True:
        if self.is_commutative is True:return S.Zero
        if right:raise NotImplementedError('Right scalar coefficients require commutative factors')
        return _native_coeff(self,x,n)
    wanted,_=target.args_cnc(cset=True,warn=False)
    result=[]
    for term in terms:
        factors,ordered=term.args_cnc(cset=True,warn=False)
        if wanted.issubset(factors):
            result.append(Mul(*(list(factors-wanted)+ordered)))
    return Add(*result)


def as_ordered_terms(self, order=None, data=False):
    if order is None:
        ordered=_native_ordered_terms(self)
        if not data:return ordered
        terms,generators=self.as_terms()
        # Equality alone merges alpha-equivalent substitutions whose bound
        # names differ. Presentation keys preserve the constructor data.
        by_presentation={(term,term.sort_key()):information for term,information in terms}
        return [(term,by_presentation[(term,term.sort_key())]) for term in ordered],generators
    from .polys.orderings import monomial_key
    reverse=isinstance(order,str) and order.startswith('rev-')
    selected=order[4:] if reverse else order
    monomial=monomial_key(selected)
    def descending(value):
        return tuple(descending(item) for item in value) if isinstance(value,tuple) else -value
    def key(item):
        _,((real,imaginary),powers,noncommutative)=item
        return (descending(monomial(powers)),
                tuple(factor.sort_key(order=selected) for factor in noncommutative),
                ((bool(imaginary),imaginary),(real,imaginary)))
    terms,generators=self.as_terms()
    if any(getattr(term,'is_Order',False) for term,_ in terms):
        ordinary=[item for item in terms if not getattr(item[0],'is_Order',False)]
        remainders=[item for item in terms if getattr(item[0],'is_Order',False)]
        ordered=sorted(ordinary,key=key,reverse=True)+sorted(remainders,key=key,reverse=True)
    else:ordered=sorted(terms,key=key,reverse=reverse)
    return (ordered,generators) if data else [term for term,_ in ordered]


def as_ordered_factors(self, order=None):
    if order is None:return _native_ordered_factors(self)
    if not self.is_Mul:return [self]
    commutative,noncommutative=self.args_cnc()
    return sorted(commutative,key=lambda factor:factor.sort_key(order=order))+noncommutative


def refine(self, assumption=True, **options):
    from . import refine as refine_expression
    if 'assumptions' in options:
        if assumption is not True:raise TypeError('Supply assumption only once')
        assumption=options.pop('assumptions')
    if options:raise TypeError('Unknown refine options: '+', '.join(sorted(options)))
    return refine_expression(self,assumption)


def collect(self, syms, func=None, evaluate=True, exact=False, distribute_order_term=True):
    from . import collect as collect_expression
    return collect_expression(self,syms,func,evaluate,exact,distribute_order_term)


def as_coefficients_dict(self, *syms):
    """Collect numeric coefficients, or factors independent of the supplied generators."""
    from collections import defaultdict
    from . import S
    if not syms:
        if self.is_Add:
            result=defaultdict(int)
            for term in self.args:
                coefficient,generator=term.as_coeff_Mul()
                result[generator]+=coefficient
            return result
        return defaultdict(int,_native_coefficients_dict(self))
    result=defaultdict(int)
    independent,dependent=self.as_independent(*syms,as_Add=True)
    terms=dependent.args if dependent.is_Add else (dependent,)
    for term in terms:
        if term==0:continue
        coefficient,generator=term.as_independent(*syms,as_Add=False)
        result[generator]+=coefficient
    if independent!=0:result[S.One]+=independent
    return result


def count_ops(self, visual=False):
    from ._operation_counts import count_ops as count
    return count(self,visual=visual)


def nsimplify(self, constants=(), tolerance=None, full=False):
    from . import nsimplify as simplify_number
    return simplify_number(self,constants=constants,tolerance=tolerance,full=full)


def _structural_nodes(expression):
    from ._matrixexpr_methods import _children, _TupleView
    from . import Expr, Set, _wrap_set
    pending=[expression]
    seen={}
    while pending:
        node=pending.pop()
        if type(node) in (tuple,list):
            pending.extend(node)
            continue
        if isinstance(node,int):
            node=_core.sympify(node)
        if isinstance(node,Expr):node=_wrap_set(node)
        try:
            underlying=node._e if isinstance(node,Set) else node
            key=(type(node),id(node) if isinstance(node,_TupleView) else getattr(underlying,'_dag_id',node))
            if key in seen:
                continue
            seen[key]=node
        except TypeError:
            pass
        yield node
        pending.extend(_children(node))


def has(self, *patterns):
    """Native literal containment plus structural class/function queries."""
    from . import Function
    from ._simptools import _head_name
    from ._structural_queries import matcher
    literals,tests=[],[]
    for pattern in patterns:
        # Scalar native expressions cannot contain None. SymPy accepts it as
        # a query and returns False rather than trying to coerce it to Expr.
        if pattern is None:
            continue
        foreign_undefined=isinstance(pattern,type) and any(
            base.__name__=='UndefinedFunction' and base.__module__.startswith('sympy.')
            for base in type(pattern).__mro__)
        if foreign_undefined:
            literals.append(pattern)
        elif isinstance(pattern,_core.Head):
            tests.append(lambda node,head=pattern:getattr(node,'func',None)==head)
        elif pattern is Function:
            tests.append(matcher(pattern,matching=False))
        elif isinstance(pattern,type):
            tests.append(matcher(pattern,matching=False))
        elif not isinstance(pattern,Function) and (name:=_head_name(pattern)):
            relations={'Eq':'==','Ne':'!=','Lt':'<','Le':'<=','Gt':'>','Ge':'>='}
            name='Rel(%s)' % relations[name] if name in relations else name
            tests.append(lambda node,name=name:
                not getattr(node,'is_UndefinedFunction',False)
                and getattr(getattr(node,'func',None),'__name__',None)==name)
        else:
            literals.append(pattern)
    if literals and _native_has(self,*literals):
        return True
    return bool(tests) and any(test(node) for node in _structural_nodes(self) for test in tests)


def match(self, pattern, old=False):
    from . import sympify, expand
    pattern=sympify(pattern)
    result=_native_match(self,pattern)
    if result is not None or not old or self.is_Symbol:
        return result
    wilds={node for node in _structural_nodes(pattern)
           if getattr(node,'is_Symbol',False) and str(node)==node.name+'_'}
    if len(wilds)!=1:return None
    wild=next(iter(wilds))
    residual=pattern-self
    slope=residual.diff(wild)
    if slope==0 or slope.has(wild):return None
    offset=residual.xreplace({wild:_core.sympify(0)})
    candidate=-offset/slope
    if candidate.has(wild):return None
    # Let the native binder check wildcard exclusions on the solved value.
    if _native_match(candidate,wild) is None:return None
    if expand(pattern.xreplace({wild:candidate})-self)!=0:return None
    return {wild:candidate}


def matches(self, expr, repl_dict=None, old=False):
    """Match this pattern against expr, extending a copy of repl_dict."""
    expr = _core.sympify(expr)
    bindings = {} if repl_dict is None else dict(repl_dict)
    # Arithmetic and Wild matches can replace previous bindings in SymPy.
    # Function arguments, in contrast, substitute existing bindings first.
    wild = self.is_Symbol and str(self) == self.name + '_'
    arithmetic = self.is_Add or self.is_Mul or self.is_Pow or wild
    result = expr.match(self if arithmetic else self.xreplace(bindings),old=old)
    if result is None:
        return None
    bindings.update(result)
    return bindings


def is_hypergeometric(self, k):
    """Test whether the consecutive-term ratio is rational in k."""
    from . import hypersimp
    pending, seen = [self], set()
    while pending:
        node = pending.pop()
        if node in seen:
            continue
        seen.add(node)
        if node.is_Piecewise:
            return None
        pending.extend(node.args)
    return hypersimp(self, k) is not None


def to_nnf(self, simplify=True):
    from . import to_nnf as convert
    return convert(self, simplify=simplify)


def _direction(self, x, cdir):
    """First nonzero derivative coefficient at zero, along cdir."""
    from . import S, log
    x, cdir = _core.sympify(x), _core.sympify(cdir)
    if self.is_zero:
        return S.Zero
    current = self
    # Flat functions may have all derivatives zero at the origin. Bound the
    # search rather than returning a fabricated nonzero direction.
    for order in range(1, 257):
        current = current.diff(x)
        coefficient = current.subs(x, 0)
        if coefficient == S.NaN:
            coefficient = current.limit(x, 0)
        if coefficient == S.ComplexInfinity:
            try:
                coefficient, _ = current.leadterm(x)
                if coefficient.has(log(x)):
                    coefficient = current.limit(x, 0)
            except ValueError:
                coefficient = current.limit(x, 0)
        if coefficient != 0 or current == 0:
            return coefficient*cdir**order
    raise NotImplementedError('Direction requires more than 256 derivatives at zero')


class _Direction:
    """Dispatch Limit.dir metadata and Expr.dir(x, cdir) on the shared type."""
    def __get__(self, instance, owner):
        if instance is None:
            return _direction
        if instance.is_Limit:
            return str(instance.args[3])
        return _direction.__get__(instance, owner)


def install():
    from ._ordered_keys import install as install_ordered_keys
    install_ordered_keys()
    # SymPy's sympify returns objects carrying its ``__sympy__`` marker
    # unchanged, so proxies forwarding attribute access to an expression
    # behave as that expression. Exposing ``__sympy__`` itself would make the
    # reference library keep native objects unconverted; this tuple-valued
    # marker survives attribute forwarding (non-callable, non-Expr values are
    # returned as is) and names the forwarded expression.
    _core.Expr._symox_forwarded_ = property(lambda self: (self,))
    _core.Expr.has = has
    _core.Expr.as_coefficients_dict = as_coefficients_dict
    _core.Expr.coeff = coeff
    _core.Expr.count_ops = count_ops
    _core.Expr.nsimplify = nsimplify
    _core.Expr.as_ordered_terms = as_ordered_terms
    _core.Expr.as_ordered_factors = as_ordered_factors
    _core.Expr.refine = refine
    _core.Expr.collect = collect
    from ._expression_replace import replace
    _core.Expr.replace = replace
    from ._expr_leading import leadterm, as_leading_term
    _core.Expr.leadterm = leadterm
    _core.Expr.as_leading_term = as_leading_term
    from ._structural_queries import find, count
    _core.Expr.find = find
    _core.Expr.count = count
    _core.Expr.matches = matches
    _core.Expr.match = match
    _core.Expr.is_hypergeometric = is_hypergeometric
    _core.Expr.to_nnf = to_nnf
    _core.Expr.dir = _Direction()
    from ._expr_analysis import (
        FromIter, as_dummy, aseries, canonical_variables,
        extract_branch_factor, is_meromorphic,
    )
    _core.Expr.fromiter = FromIter()
    _core.Expr.as_dummy = as_dummy
    _core.Expr.canonical_variables = property(canonical_variables)
    _core.Expr.aseries = aseries
    _core.Expr.extract_branch_factor = extract_branch_factor
    _core.Expr.is_meromorphic = is_meromorphic
    from .series import fps, fourier_series
    _core.Expr.fps = fps
    _core.Expr.fourier_series = lambda self, limits=None: fourier_series(self, limits)
