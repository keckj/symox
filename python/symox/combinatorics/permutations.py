"""Finite permutations in array and cycle notation."""
import math
import operator


class Permutation:
    """A bijection of range(size). Multiplication applies the left factor first."""
    args = ()
    is_Atom = True
    is_Number = is_Symbol = is_Add = is_Mul = is_Pow = is_Function = False
    free_symbols = frozenset()

    @property
    def func(self):
        return type(self)

    @classmethod
    def class_key(cls):
        return (2, 0, 'Permutation')

    def sort_key(self, order=None):
        from .. import S
        return (self.class_key(), (1, (str(self),)), S.One.sort_key(), S.One)

    def doit(self, **hints):
        return self

    def to_sympy(self):
        from sympy.combinatorics import Permutation
        return Permutation(self.array_form)

    _sympy_ = to_sympy

    def __init__(self, *args, size=None, **kwargs):
        # Reference construction accepts historical/extension keywords without
        # using them. Validation and the explicit size control still apply.
        self._rank=None
        self._parity_cache=None
        self._cycle_structure_cache=None
        if size is not None:
            size = int(size)
        if not args:
            self._images = tuple(range(size or 0))
            return
        if len(args) > 1:
            self._images = tuple(Cycle(*args).list(size))
            return
        images = args[0]
        if isinstance(images, Cycle):
            self._images = tuple(images.list(size))
            return
        if isinstance(images, Permutation):
            images = images.array_form
        else:
            try:
                images = list(images)
            except TypeError:
                maximum = operator.index(images)
                if size is not None and maximum + 1 > size:
                    raise ValueError('size is too small for the maximum index')
                self._images = tuple(range(maximum + 1))
                return
        nested = [isinstance(i, (list, tuple, range)) for i in images]
        if any(nested):
            if not all(nested):
                raise ValueError('mixed array and cycle notation')
            cycles = [[int(i) for i in cycle] for cycle in images]
            for cycle in cycles:
                if any(i < 0 for i in cycle) or len(set(cycle)) != len(cycle):
                    raise ValueError('invalid cycle indices')
            degree = max((max(cycle, default=-1) for cycle in cycles), default=-1) + 1
            degree = max(degree, size or 0)
            # Maintain inverse positions so overlapping cycles compose in
            # O(degree + number of supplied cycle entries), without rescanning
            # the entire image array for each cycle.
            images = list(range(degree))
            inverse = images.copy()
            for cycle in cycles:
                positions = [inverse[i] for i in cycle]
                for position, target in zip(positions, cycle[1:] + cycle[:1]):
                    images[position] = target
                    inverse[target] = position
            self._images = tuple(images)
            return
        images = tuple(int(i) for i in images)
        if set(images) != set(range(len(images))):
            raise ValueError('Array must contain each index exactly once')
        if size is not None:
            if size < len(images):
                raise ValueError('size cannot shrink a permutation')
            images += tuple(range(len(images), size))
        self._images = images

    @property
    def array_form(self):
        return list(self._images)

    @property
    def size(self):
        return len(self._images)

    def __call__(self, *indices):
        if len(indices) != 1:
            return self * type(self)(Cycle(*indices), size=self.size)
        value = indices[0]
        try:
            index = operator.index(value)
        except TypeError:
            if not hasattr(value, '__iter__'):
                raise TypeError('permutation index must be an integer') from None
            if len(value) != self.size:
                raise TypeError('sequence must have the permutation size')
            return [value[i] for i in self._images]
        if index < 0 or index > self.size:
            raise TypeError('permutation index out of range')
        return self._images[index]

    def __mul__(self, other):
        if not isinstance(other, Permutation):
            return NotImplemented
        size = max(self.size, other.size)
        left = self._images if self.size == size else self._images + tuple(range(self.size, size))
        right = other._images if other.size == size else other._images + tuple(range(other.size, size))
        return type(self)._from_images([right[index] for index in left])

    @classmethod
    def _from_images(cls, images):
        """Internal construction from a bijection established by an algorithm.

        Public input still goes through all constructor checks. Subclasses keep
        their constructor behavior; only this concrete class uses the shortcut.
        """
        if cls is not Permutation:
            return cls(images)
        result = object.__new__(cls)
        result._images = tuple(images)
        result._rank = result._parity_cache = result._cycle_structure_cache = None
        return result

    def __invert__(self):
        inverse = [0]*self.size
        for i, image in enumerate(self._images):
            inverse[image] = i
        return type(self)._from_images(inverse)

    def __pow__(self, n):
        exponent = operator.index(n)
        # Small powers need only a few array-composition passes. Cycle
        # decomposition remains linear in degree for arbitrarily large powers.
        if abs(exponent) <= 64:
            base = list((~self)._images if exponent < 0 else self._images)
            remaining = abs(exponent)
            result = None
            while remaining:
                if remaining & 1:
                    result = base if result is None else [base[i] for i in result]
                remaining >>= 1
                if remaining >= 2 and not remaining & 1:
                    base = [base[base[base[i]]] for i in base]
                    remaining >>= 1
                elif remaining:
                    base = [base[i] for i in base]
            return type(self)._from_images(range(self.size) if result is None else result)
        images = list(range(self.size))
        for cycle in self.cyclic_form:
            shift = exponent % len(cycle)
            for index,value in enumerate(cycle):
                images[value] = cycle[(index+shift)%len(cycle)]
        return type(self)._from_images(images)

    @property
    def cyclic_form(self):
        seen, cycles = set(), []
        for start in range(self.size):
            if start in seen:
                continue
            cycle, i = [], start
            while i not in seen:
                seen.add(i)
                cycle.append(i)
                i = self._images[i]
            if len(cycle) > 1:
                cycles.append(cycle)
        return cycles

    def order(self):
        order = 1
        for cycle in self.cyclic_form:
            length = len(cycle)
            order = order*length // math.gcd(order, length)
        return order

    def signature(self):
        return (-1)**sum(len(cycle)-1 for cycle in self.cyclic_form)

    @property
    def is_even(self):
        return self.signature() == 1

    @property
    def is_Identity(self):
        return all(i == image for i, image in enumerate(self._images))

    def __hash__(self):
        return hash(self._images)

    def __eq__(self, other):
        if not isinstance(other, Permutation):
            return NotImplemented
        return self._images == other._images

    def __repr__(self):
        return 'Permutation(%r)' % self.array_form

    def __str__(self):
        prefix = '(%d)' % (self.size-1) if self.size and self._images[-1] == self.size-1 else ''
        cycles = ''.join('(%s)' % ' '.join(map(str, cycle)) for cycle in self.cyclic_form)
        return prefix+cycles or '()'


class Cycle(dict):
    """A finite cycle mapping; composition applies the left mapping first."""
    def __init__(self,*args):
        super().__init__()
        if len(args)==1 and isinstance(args[0],Permutation):
            self.update(enumerate(args[0].array_form))
        elif len(args)==1 and isinstance(args[0],dict):
            self.update(args[0])
            if set(self)!=set(self.values()):raise ValueError('cycle mapping must be a bijection')
        elif args:
            values=[operator.index(v) for v in args]
            if min(values)<0 or len(set(values))!=len(values):raise ValueError('invalid cycle indices')
            self.update(zip(values,values[1:]+values[:1]))

    @property
    def size(self):
        return max(self,default=-1)+1

    def __missing__(self,key):
        return key

    def __call__(self,*args):
        other=Cycle(*args)
        keys=set(self)|set(other)
        return Cycle({i:other[self[i]] for i in keys})

    def list(self,size=None):
        if size is None:
            if not self:raise ValueError('must give size for empty Cycle')
            size=self.size
        else:
            largest_moved=max((i for i in self if self[i]!=i),default=0)
            size=max(operator.index(size),largest_moved+1)
        return [self[i] for i in range(size)]

    def copy(self):
        return Cycle(self)

    def __repr__(self):
        return str(Permutation(self.list()))

Perm=Permutation


class AppliedPermutation:
    def __new__(cls,perm,x,evaluate=None):
        from symox import sympify,Integer
        if not isinstance(perm,Permutation):raise ValueError('a permutation application requires a permutation')
        x=sympify(x)
        if evaluate is not False and x.is_Integer:
            index=int(x)
            if not 0<=index<perm.size:raise IndexError('permutation argument out of range')
            return Integer(perm(index))
        result=object.__new__(cls);result.perm,result.x=perm,x;return result
    args=property(lambda self:(self.perm,self.x))
    free_symbols=property(lambda self:self.x.free_symbols)
    def _symox_(self):
        from symox import Piecewise,Eq,nan
        return Piecewise(*[(self.perm(i),Eq(self.x,i)) for i in range(self.perm.size)],(nan,True))
    def doit(self,**hints):return type(self)(self.perm,self.x)
    def subs(self,*args,**kwargs):return type(self)(self.perm,self.x.subs(*args,**kwargs))
    def __repr__(self):return 'AppliedPermutation(%r, %s)'%(self.perm,self.x)
    def __eq__(self,other):return isinstance(other,AppliedPermutation) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))
    def __add__(self,other):return self._symox_()+other
    __radd__=__add__
    def __mul__(self,other):return self._symox_()*other
    __rmul__=__mul__


def _apply(self,index):return AppliedPermutation(self,index)
Permutation.apply=_apply

from ._permutation_algorithms import install as _install_algorithms
_install_algorithms(Permutation)
