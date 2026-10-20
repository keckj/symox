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


def imageset(*args):
    if len(args)<2:
        raise TypeError('imageset requires a map and base sets')
    if isinstance(args[0],(Expr,)) and getattr(args[0],'is_Lambda',False):
        function,bases=args[0],args[1:]
    elif callable(args[0]) and not isinstance(args[0],Expr):
        import inspect
        from . import Dummy
        function=args[0]
        n=len(inspect.signature(function).parameters)
        variables=tuple(Dummy('x%d'%i) for i in range(n))
        function,bases=Lambda(variables,function(*variables)),args[1:]
    else:
        if len(args)<3:
            raise TypeError('imageset requires variables, expression, and base sets')
        function,bases=Lambda(args[0],args[1]),args[2:]
    if len(function.variables)!=len(bases):
        raise ValueError('map arity does not match base sets')
    if all(s.is_finite_set is True for s in bases):
        size=1
        for s in bases:
            size*=0 if s==EmptySet else len(s)
        if size>100000:
            raise NotImplementedError('explicit image exceeds 100000 points')
        return FiniteSet(*(function(*values) for values in product(*(list(s) if s!=EmptySet else [] for s in bases))))
    return ImageSet(function,*bases)
