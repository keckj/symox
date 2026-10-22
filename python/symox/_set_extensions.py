"""Structural sets whose membership is evaluated with native scalar predicates."""
from itertools import combinations, product
from ._core import Expr
from . import (Set, FiniteSet, EmptySet, ProductSet, Union, ImageSet,
               Lambda, Interval, S, sympify, And, Or, Eq, I, re, im, Abs, arg, pi)


class _StructuralSet(Set):
    def __new__(cls,*args,**kwargs):
        return object.__new__(cls)

    @property
    def func(self):
        return type(self)

    @property
    def class_name(self):
        return type(self).__name__

    def __repr__(self):
        return '%s(%s)' % (self.class_name, ', '.join(map(str,self.args)))
    __str__=__repr__

    def __eq__(self,other):
        return type(self) is type(other) and self.args==other.args

    def __hash__(self):
        return hash((type(self),self.args))

    def __reduce__(self):
        return type(self),self.args

    def _symox_(self):
        raise NotImplementedError('this structural set has no native scalar representation')

    def to_sympy(self):
        import sympy
        return getattr(sympy,type(self).__name__)(*(a.to_sympy() if hasattr(a,'to_sympy') else a for a in self.args))
    _sympy_=to_sympy

    @property
    def free_symbols(self):
        return set().union(*(getattr(a,'free_symbols',set()) for a in self.args))

    @property
    def is_Atom(self):
        return False

    def subs(self,*args,**kwargs):
        return self.func(*(a.subs(*args,**kwargs) if hasattr(a,'subs') else a for a in self.args))

    def xreplace(self,rule):
        if self in rule:
            return rule[self]
        return self.func(*(a.xreplace(rule) if hasattr(a,'xreplace') else a for a in self.args))

    def _algebra(self,*args):
        raise NotImplementedError('general algebra of structural sets is not implemented')
    union=intersect=intersection=complement=symmetric_difference=_algebra
    __or__=__ror__=__add__=__radd__=__and__=__rand__=__sub__=__rsub__=__xor__=__mul__=__pow__=_algebra
    as_relational=isdisjoint=is_disjoint=_algebra

    def is_subset(self,other):
        if self==other:
            return True
        return None
    issubset=is_subset

    @property
    def measure(self):
        raise NotImplementedError('measure is undefined for this structural set')
    inf=sup=boundary=interior=closure=measure


class PowerSet(_StructuralSet):
    def __init__(self,arg,evaluate=None):
        if not isinstance(arg,Set):
            raise ValueError('PowerSet requires a set')
        self.arg=arg

    @property
    def args(self):
        return (self.arg,)

    @property
    def is_empty(self):
        return False

    @property
    def is_finite_set(self):
        return self.arg.is_finite_set

    @property
    def is_iterable(self):
        return self.arg.is_finite_set is True

    def __len__(self):
        return 2**len(self.arg)

    def __iter__(self):
        if self.arg.is_finite_set is not True:
            raise NotImplementedError('power set iteration requires a finite base')
        values=[] if self.arg==EmptySet else list(self.arg)
        if len(values)>16:
            raise NotImplementedError('power set enumeration is limited to 65536 subsets')
        yield EmptySet
        for j,value in enumerate(values):
            for mask in range(1<<j):
                yield FiniteSet(value,*(values[i] for i in range(j) if mask>>i&1))

    def contains(self,other):
        if not isinstance(other,Set):
            return S.false
        result=other.is_subset(self.arg)
        if result is None:
            raise NotImplementedError('subset membership could not be decided')
        return sympify(result)

    def is_subset(self,other):
        if isinstance(other,PowerSet):
            return self.arg.is_subset(other.arg)
        return super().is_subset(other)

    def rewrite(self,target):
        if target is FiniteSet and self.is_finite_set:
            return FiniteSet(*self)
        return self


class DisjointUnion(_StructuralSet):
    def __init__(self,*sets):
        if any(not isinstance(s,Set) for s in sets):
            raise TypeError('DisjointUnion requires sets')
        self.sets=tuple(sets)

    @property
    def args(self):
        return self.sets

    @property
    def is_empty(self):
        statuses=[s.is_empty for s in self.sets]
        return True if all(s is True for s in statuses) else False if False in statuses else None

    @property
    def is_finite_set(self):
        statuses=[s.is_finite_set for s in self.sets]
        return False if False in statuses else True if all(s is True for s in statuses) else None

    @property
    def is_iterable(self):
        return all(s.is_iterable or s==EmptySet for s in self.sets)

    def __len__(self):
        return sum(0 if s==EmptySet else len(s) for s in self.sets)

    def __iter__(self):
        # Round-robin also handles infinite iterable component sets fairly.
        active=[(i,iter(s)) for i,s in enumerate(self.sets) if s!=EmptySet]
        while active:
            next_active=[]
            for i,it in active:
                try:
                    yield (next(it),sympify(i))
                    next_active.append((i,it))
                except StopIteration:
                    pass
            active=next_active

    def contains(self,other):
        if not isinstance(other,(tuple,list)) or len(other)!=2:
            return S.false
        value,index=other
        index=sympify(index)
        return Or(*(And(Eq(index,i),s.contains(value)) for i,s in enumerate(self.sets)))

    def rewrite(self,target):
        if target is Union:
            return Union(*(ProductSet(s,FiniteSet(i)) for i,s in enumerate(self.sets)))
        return self


class ComplexRegion(_StructuralSet):
    def __init__(self,sets,polar=False):
        if polar not in (True,False):
            raise ValueError('polar must be boolean')
        pieces=sets.args if isinstance(sets,Union) else (sets,)
        if any(not isinstance(s,ProductSet) or len(s.args)!=2 for s in pieces):
            raise ValueError('complex regions require unions of two-factor Cartesian products')
        self.sets,self.polar=sets,polar
        if polar:
            # Nonnegative radii and principal angles provide an unambiguous
            # membership test. General angular normalization needs set images.
            for piece in pieces:
                radius,theta=piece.args
                if radius.is_subset(Interval(0,S.Infinity)) is not True or theta.is_subset(Interval(0,2*pi)) is not True:
                    raise NotImplementedError('polar regions require r >= 0 and angles in [0, 2*pi]')

    @property
    def args(self):
        return self.sets,self.polar

    @property
    def psets(self):
        return self.sets.args if isinstance(self.sets,Union) else (self.sets,)

    @property
    def a_interval(self):
        return Union(*(s.args[0] for s in self.psets))

    @property
    def b_interval(self):
        return Union(*(s.args[1] for s in self.psets))

    @property
    def is_empty(self):
        return self.sets.is_empty

    @property
    def is_finite_set(self):
        return self.sets.is_finite_set

    @property
    def is_iterable(self):
        return all(all(s.is_iterable or s==EmptySet for s in piece.args) for piece in self.psets)

    def __iter__(self):
        seen=set()
        from . import exp
        for piece in self.psets:
            if any(s.is_finite_set is not True for s in piece.args):
                raise NotImplementedError('complex region enumeration requires finite factors')
            for a,b in product(*(list(s) if s!=EmptySet else [] for s in piece.args)):
                value=a*exp(I*b) if self.polar else a+I*b
                if value not in seen:
                    seen.add(value)
                    yield value

    def contains(self,other):
        other=sympify(other)
        if other.is_number and other.is_finite is False:
            return S.false
        if not self.polar:
            return Or(*(And(a.contains(re(other)),b.contains(im(other))) for a,b in (s.args for s in self.psets)))
        angle=arg(other)
        if angle.is_negative is True:
            angle+=2*pi
        if angle.is_nonnegative is not True and other!=0:
            raise NotImplementedError('the principal angle could not be normalized')
        return Or(*(And(a.contains(Abs(other)),Or(b.contains(angle),b.contains(angle+2*pi))) if other!=0 else a.contains(0) for a,b in (s.args for s in self.psets)))

    @property
    def measure(self):
        from . import Intersection, Integer
        pieces=self.psets
        if len(pieces)>12:
            raise NotImplementedError('region area unions are limited to 12 rectangles')
        area=Integer(0)
        for size in range(1,len(pieces)+1):
            for selected in combinations(pieces,size):
                a=Intersection(*(p.args[0] for p in selected))
                b=Intersection(*(p.args[1] for p in selected))
                if a==EmptySet or b==EmptySet or a.measure==0 or b.measure==0:
                    continue
                if self.polar:
                    if not isinstance(a,Interval):
                        raise NotImplementedError('polar area requires interval radii')
                    term=(a.end**2-a.start**2)/2*b.measure
                else:
                    term=a.measure*b.measure
                area+=(-1)**(size+1)*term
        return area


def _invertible_in_set(func, setv):
    """The reference ``is_function_invertible_in_set``."""
    from . import Dummy, exp, log
    if func in (exp, log):
        return True
    u = Dummy('u')
    fdiff = func(u).diff(u)
    if (fdiff > 0) == True or (fdiff < 0) == True:  # noqa: E712
        return True
    return None


def _image_interval(f, x):
    """The reference ``_set_function(Lambda, Interval)``."""
    from . import (Dummy, Intersection, Complement, Min, Max, diff, limit,
                   true)
    from .calculus.singularities import singularities
    from .solvers.solveset import solveset
    from .solvers._periodic_inequality import _nsort_reals
    expr = f.expr
    if len(expr.free_symbols) > 1 or len(f.variables) != 1:
        return None
    var = f.variables[0]
    if not var.is_real:
        if expr.subs(var, Dummy(real=True)).is_real is False:
            return None
    if expr.is_Piecewise:
        result = S.EmptySet
        domain_set = x
        for (p_expr, p_cond) in expr.args:
            if p_cond is true or p_cond == True:  # noqa: E712
                intrvl = domain_set
            else:
                intrvl = Intersection(domain_set, p_cond.as_set())
            if p_expr.is_Number:
                image = FiniteSet(p_expr)
            else:
                image = imageset(Lambda(var, p_expr), intrvl)
            result = Union(result, image)
            domain_set = Complement(domain_set, intrvl)
            if domain_set == EmptySet:
                break
        return result
    if not x.start.is_comparable or not x.end.is_comparable:
        return None
    try:
        sing = list(singularities(expr, var, x))
        if len(sing) > 1:
            sing = _nsort_reals(sing)
    except NotImplementedError:
        return None
    if x.left_open:
        _start = limit(expr, var, x.start, dir="+")
    elif x.start not in sing:
        _start = f(x.start)
    if x.right_open:
        _end = limit(expr, var, x.end, dir="-")
    elif x.end not in sing:
        _end = f(x.end)
    if len(sing) == 0:
        soln_expr = solveset(diff(expr, var), var)
        if not (isinstance(soln_expr, FiniteSet) or soln_expr == EmptySet):
            return None
        solns = list(soln_expr)
        extr = [_start, _end] + [f(i) for i in solns if i.is_real and i in x]
        start, end = Min(*extr), Max(*extr)
        left_open, right_open = False, False
        if _start <= _end:
            if start == _start and start not in solns:
                left_open = x.left_open
            if end == _end and end not in solns:
                right_open = x.right_open
        else:
            if start == _end and start not in solns:
                left_open = x.right_open
            if end == _start and end not in solns:
                right_open = x.left_open
        return Interval(start, end, left_open, right_open)
    return imageset(f, Interval(x.start, sing[0], x.left_open, True)) + \
        Union(*[imageset(f, Interval(sing[i], sing[i + 1], True, True))
                for i in range(0, len(sing) - 1)]) + \
        imageset(f, Interval(sing[-1], x.end, True, x.right_open))


def _image_range(f, self):
    """The reference ``_set_function(FunctionUnion, Range)``."""
    from . import Range, expand_mul
    if not self:
        return EmptySet
    if not isinstance(f.expr, Expr):
        return None
    if self.size == 1:
        return FiniteSet(f(self[0]))
    x = f.variables[0]
    expr = f.expr
    if x not in expr.free_symbols or x in expr.diff(x).free_symbols:
        return None
    if self.start.is_finite:
        F = f(self.step*x + self.start)
    else:
        F = f(-self.step*x + self[-1])
    F = expand_mul(F)
    if F != expr:
        return imageset(x, F, Range(self.size))
    return None


def _match_real_imag(expr):
    """The reference ``match_real_imag``."""
    r_, i_ = expr.as_independent(I, as_Add=True)
    if i_ == 0 and r_.is_real:
        return (r_, i_)
    i_ = i_.as_coefficient(I)
    if i_ and i_.is_real and r_.is_real:
        return (r_, i_)
    return (None, None)


def _image_integers(f, self):
    """The reference ``_set_function(FunctionUnion, Integers)``."""
    from . import Add, Wild, Float
    expr = f.expr
    if not isinstance(expr, Expr):
        return None
    n = f.variables[0]
    if expr == abs(n):
        return S.Naturals0
    c = f(0)
    fx = f(n) - c
    f_x = f(-n) - c
    neg_count = lambda e: sum(t.could_extract_minus_sign() for t in Add.make_args(e))
    if neg_count(f_x) < neg_count(fx):
        expr = f_x + c
    a = Wild('a', exclude=[n])
    b = Wild('b', exclude=[n])
    match = expr.match(a*n + b)
    if match and match[a] and (not match[a].atoms(Float) and not match[b].atoms(Float)):
        a, b = match[a], match[b]
        if a in [1, -1]:
            b = Add(*[bi for bi in Add.make_args(b) if not bi.is_integer])
        if b.is_number and a.is_real:
            br, bi = _match_real_imag(b)
            if br and br.is_comparable and a.is_comparable:
                br %= a
                b = br + S.ImaginaryUnit*bi
        elif b.is_number and a.is_imaginary:
            br, bi = _match_real_imag(b)
            ai = a/S.ImaginaryUnit
            if bi and bi.is_comparable and ai.is_comparable:
                bi %= ai
                b = br + S.ImaginaryUnit*bi
        expr = a*n + b
    if expr != f.expr:
        return ImageSet(Lambda(n, expr), S.Integers)
    return None


def _image_naturals(f, self):
    """The reference ``_set_function(FunctionUnion, Naturals)``."""
    from . import Range, oo
    expr = f.expr
    if not isinstance(expr, Expr):
        return None
    x = f.variables[0]
    if not expr.free_symbols - {x}:
        if expr == abs(x):
            if self == S.Naturals:
                return self
            return S.Naturals0
        step = expr.coeff(x)
        c = expr.subs(x, 0)
        if c.is_Integer and step.is_Integer and expr == step*x + c:
            if self == S.Naturals:
                c += step
            if step > 0:
                if step == 1:
                    if c == 0:
                        return S.Naturals0
                    elif c == 1:
                        return S.Naturals
                return Range(c, oo, step)
            return Range(c, -oo, step)
    return None


def _set_function(f, x):
    """The reference ``set_function`` dispatch for a ``Lambda``."""
    from . import Symbol, Intersection, Range, oo
    _x = Symbol('x')
    if x == EmptySet:
        return x
    if isinstance(x, FiniteSet):
        return FiniteSet(*map(f, x))
    if x == S.Integers:
        return _image_integers(f, x)
    if x == S.Naturals or x == S.Naturals0:
        return _image_naturals(f, x)
    if isinstance(x, type(S.Reals)):
        if not isinstance(f.expr, Expr):
            return None
        return _image_interval(f, Interval(-oo, oo))
    if isinstance(x, Interval):
        return _image_interval(f, x)
    if isinstance(x, Range):
        return _image_range(f, x)
    if isinstance(x, Union):
        return Union(*(imageset(f, a) for a in x.args))
    if isinstance(x, Intersection):
        if _invertible_in_set(f, x):
            return Intersection(*(imageset(f, a) for a in x.args))
        return ImageSet(Lambda(_x, f(_x)), x)
    return ImageSet(Lambda(_x, f(_x)), x)


def imageset(*args):
    """The reference ``imageset``."""
    from . import Symbol, Tuple
    if len(args) < 2:
        raise ValueError('imageset expects at least 2 args, got: %s' % len(args))
    if isinstance(args[0], (Symbol, tuple, Tuple)) and len(args) > 2:
        f = Lambda(args[0], args[1])
        set_list = args[2:]
    else:
        f = args[0]
        set_list = args[1:]
    if isinstance(f, Expr) and getattr(f, 'is_Lambda', False):
        pass
    elif callable(f):
        import inspect
        from . import Dummy
        from .core.symbol import uniquely_named_symbol
        nargs = getattr(f, 'nargs', {})
        if nargs:
            if len(nargs) != 1:
                raise NotImplementedError('This function can take more than 1 arg')
            N = nargs.args[0]
            s = 'x' if N == 1 else [Symbol('x%i' % i) for i in range(1, N + 1)]
        else:
            s = inspect.signature(f).parameters
        dexpr = sympify(f(*[Dummy() for i in s]))
        var = tuple(uniquely_named_symbol(Symbol(i), dexpr) for i in s)
        f = Lambda(var, f(*var))
    else:
        raise TypeError("expecting lambda, Lambda, or FunctionClass, not '%s'." % type(f).__name__)
    if any(not isinstance(s, Set) for s in set_list):
        raise ValueError('arguments after mapping should be sets, not %s'
                         % [type(s).__name__ for s in set_list])
    signature = getattr(f, 'signature', None) or f.variables
    if len(signature) != len(set_list):
        raise ValueError('Incompatible signature')
    if len(set_list) == 1:
        set_ = set_list[0]
        try:
            if len(f.variables) != 1:
                raise TypeError
            r = _set_function(f, set_)
            if r is None:
                raise TypeError
            if r == EmptySet:
                return r
        except TypeError:
            r = ImageSet(f, set_)
        if isinstance(r, ImageSet):
            f, set_ = r.lamda, r.base_set
        if f.variables[0] == f.expr:
            return set_
        if isinstance(set_, ImageSet):
            if len(set_.lamda.variables) == 1 and len(f.variables) == 1:
                x = set_.lamda.variables[0]
                y = f.variables[0]
                return imageset(Lambda(x, f.expr.subs(y, set_.lamda.expr)), set_.base_set)
        if r is not None:
            return r
    return ImageSet(f, *set_list)
