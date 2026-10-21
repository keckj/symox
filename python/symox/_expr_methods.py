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
        return _noncommutative_coeff(self,target,right,_first)
    wanted,_=target.args_cnc(cset=True,warn=False)
    result=[]
    for term in terms:
        factors,ordered=term.args_cnc(cset=True,warn=False)
        if wanted.issubset(factors):
            result.append(Mul(*(list(factors-wanted)+ordered)))
    return Add(*result)


def _noncommutative_coeff(self, x, right, _first):
    """SymPy's ``Expr.coeff`` after ``x = x**n`` when a factor may not commute."""
    from functools import reduce
    from . import S, Add, Mul
    def incommon(l1, l2):
        if not l1 or not l2:return []
        for i in range(min(len(l1),len(l2))):
            if l1[i]!=l2[i]:return l1[:i]
        return l1[:]
    def find(l, sub, first=True):
        if not sub or not l or len(sub)>len(l):return None
        n=len(sub)
        if not first:l.reverse();sub.reverse()
        for i in range(len(l)-n+1):
            if all(l[i+j]==sub[j] for j in range(n)):break
        else:i=None
        if not first:l.reverse();sub.reverse()
        if i is not None and not first:i=len(l)-(i+n)
        return i
    co2=[];co=[]
    args=Add.make_args(self)
    self_c=self.is_commutative
    x_c=x.is_commutative
    if self_c and not x_c:return S.Zero
    if _first and self.is_Add and not self_c and not x_c:
        xargs=Mul.make_args(x)
        d=Add(*[i for i in Add.make_args(self.as_independent(x)[1]) if all(xi in Mul.make_args(i) for xi in xargs)])
        rv=coeff(d,x,right=right,_first=False)
        if not rv.is_Add or not right:return rv
        c_part,nc_part=zip(*[i.args_cnc() for i in rv.args])
        if len(set(map(tuple,c_part)))>1:return rv
        return Add(*[Mul(*i) for i in nc_part])
    one_c=self_c or x_c
    xargs,nx=x.args_cnc(cset=True,warn=bool(not x_c))
    for a in args:
        margs,nc=a.args_cnc(cset=True,warn=bool(not self_c))
        if nc is None:nc=[]
        if len(xargs)>len(margs):continue
        resid=margs.difference(xargs)
        if len(resid)+len(xargs)==len(margs):
            if one_c:co2.append(Mul(*(list(resid)+nc)))
            else:co.append((resid,nc))
    if one_c:
        return Add(*co2) if co2 else S.Zero
    if not co:return S.Zero
    if all(n==co[0][1] for r,n in co):
        ii=find(co[0][1],nx,right)
        if ii is not None:
            if not right:return Mul(Add(*[Mul(*r) for r,c in co]),Mul(*co[0][1][:ii]))
            return Mul(*co[0][1][ii+len(nx):])
    beg=reduce(incommon,(n[1] for n in co))
    if beg:
        ii=find(beg,nx,right)
        if ii is not None:
            if not right:
                gcdc=co[0][0]
                for i in range(1,len(co)):
                    gcdc=gcdc.intersection(co[i][0])
                    if not gcdc:break
                return Mul(*(list(gcdc)+beg[:ii]))
            m=ii+len(nx)
            return Add(*[Mul(*(list(r)+n[m:])) for r,n in co])
    end=list(reversed(reduce(incommon,(list(reversed(n[1])) for n in co))))
    if end:
        ii=find(end,nx,right)
        if ii is not None:
            if not right:return Add(*[Mul(*(list(r)+n[:-len(end)+ii])) for r,n in co])
            return Mul(*end[ii+len(nx):])
    hit=None
    for i,(r,n) in enumerate(co):
        ii=find(n,nx,right)
        if ii is not None:
            if not hit:hit=ii,r,n
            else:break
    else:
        if hit:
            ii,r,n=hit
            if not right:return Mul(*(list(r)+n[:ii]))
            return Mul(*n[ii+len(nx):])
    return S.Zero


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
    from . import Set, _wrap_set
    from ._core import Expr
    pending=[expression]
    seen={}
    while pending:
        node=pending.pop()
        kind=type(node)
        if kind is tuple or kind is list:
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
    """The reference pattern matching (symox._matching)."""
    from . import sympify
    from ._matching import matches as structural_matches
    return structural_matches(sympify(pattern), self, old=old)


def matches(self, expr, repl_dict=None, old=False):
    """Match this pattern against expr, extending a copy of repl_dict."""
    from ._matching import matches as structural_matches
    return structural_matches(self, _core.sympify(expr), repl_dict, old=old)


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
