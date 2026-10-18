"""Finite permutations in array form, used by exact Galois group results."""
import math
import operator


class Permutation:
    """A bijection of range(size), constructed from its array of images.

    Array construction is supported; SymPy's positional cycle constructor is
    not part of this subset. Multiplication applies the left factor first.
    """
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

    def __init__(self, images, *, size=None):
        if isinstance(images, Cycle):
            images = images.list()
        if isinstance(images, Permutation):
            images = images.array_form
        images = tuple(operator.index(i) for i in images)
        if sorted(images) != list(range(len(images))):
            raise ValueError('Array must contain each index exactly once')
        if size is not None:
            size = operator.index(size)
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

    def __call__(self, index):
        return self._images[operator.index(index)]

    def __mul__(self, other):
        if not isinstance(other, Permutation):
            return NotImplemented
        size = max(self.size, other.size)
        left = self._images + tuple(range(self.size, size))
        right = other._images + tuple(range(other.size, size))
        return type(self)([right[left[i]] for i in range(size)])

    def __invert__(self):
        inverse = [0]*self.size
        for i, image in enumerate(self._images):
            inverse[image] = i
        return type(self)(inverse)

    def __pow__(self, exponent):
        exponent = operator.index(exponent)
        base = self if exponent >= 0 else ~self
        exponent = abs(exponent)
        result = type(self)(range(self.size))
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

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
        size=self.size if size is None else operator.index(size)
        if size<self.size:raise ValueError('cycle size cannot be reduced')
        return [self[i] for i in range(size)]

    def copy(self):
        return Cycle(self)

    def __repr__(self):
        return str(Permutation(self.list()))
