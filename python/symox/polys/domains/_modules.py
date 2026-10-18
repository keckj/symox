"""Finite free modules, field subspaces, and exact quotient arithmetic.

Polynomial ideals use the native Groebner engine over a field. Integer ideals
use their positive gcd. Submodule elimination currently requires a field.
"""
import operator


class FreeModule:
    def __init__(self, ring, rank):
        rank = operator.index(rank)
        if rank < 0:
            raise ValueError('module rank must be nonnegative')
        self.ring, self.rank = ring, rank

    def __repr__(self):
        return f'{self.ring}**{self.rank}'

    def __eq__(self, other):
        return isinstance(other, FreeModule) and (self.ring, self.rank) == (other.ring, other.rank)

    def convert(self, value):
        if isinstance(value, ModuleElement):
            if value.module != self:
                raise ValueError('incompatible free modules')
            return value
        coordinates = tuple(self.ring.convert(a) for a in value)
        if len(coordinates) != self.rank:
            raise ValueError('wrong number of module coordinates')
        return ModuleElement(self, coordinates)

    __call__ = convert

    def __contains__(self, value):
        from ... import CoercionFailed
        try:
            self.convert(value)
            return True
        except (CoercionFailed, TypeError, ValueError):
            return False

    @property
    def basis(self):
        k = self.ring
        return tuple(self([k.one if i == j else k.zero for j in range(self.rank)])
                     for i in range(self.rank))

    def submodule(self, *generators):
        return SubModule(self, generators)


class ModuleElement:
    def __init__(self, module, coordinates):
        self.module, self.coordinates = module, coordinates

    def __iter__(self):
        return iter(self.coordinates)

    def __len__(self):
        return len(self.coordinates)

    def __getitem__(self, key):
        return self.coordinates[key]

    def __repr__(self):
        return '[' + ', '.join(map(str, self.coordinates)) + ']'

    def __eq__(self, other):
        try:
            return self.coordinates == self.module.convert(other).coordinates
        except (TypeError, ValueError):
            return False

    def __add__(self, other):
        other = self.module.convert(other)
        return self.module([a+b for a, b in zip(self, other)])

    def __neg__(self):
        return self.module([-a for a in self])

    def __sub__(self, other):
        return self + (-self.module.convert(other))

    def __mul__(self, scalar):
        scalar = self.module.ring.convert(scalar)
        return self.module([scalar*a for a in self])

    __rmul__ = __mul__


class SubModule:
    def __init__(self, container, generators):
        if not container.ring.is_Field:
            raise NotImplementedError('submodule elimination currently requires a field')
        self.container, self.ring = container, container.ring
        self.gens = tuple(container(g) for g in generators)
        rows = [list(g) for g in self.gens]
        pivots, pivot_row = [], 0
        for column in range(container.rank):
            selected = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
            if selected is None:
                continue
            rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
            scale = self.ring.revert(rows[pivot_row][column])
            rows[pivot_row] = [scale*a for a in rows[pivot_row]]
            for i in range(len(rows)):
                if i != pivot_row and rows[i][column]:
                    scale = rows[i][column]
                    rows[i] = [a-scale*b for a, b in zip(rows[i], rows[pivot_row])]
            pivots.append(column)
            pivot_row += 1
        self._rows = tuple(container(row) for row in rows[:pivot_row])
        self._pivots = tuple(pivots)

    def __repr__(self):
        return '<' + ', '.join(map(str, self.gens)) + '>'

    def reduce_element(self, value):
        value = self.container(value)
        for pivot, row in zip(self._pivots, self._rows):
            value = value - row*value[pivot]
        return value

    def contains(self, value):
        return all(not a for a in self.reduce_element(value))

    __contains__ = contains

    def is_full_module(self):
        return len(self._rows) == self.container.rank

    def is_zero(self):
        return not self._rows

    def submodule(self, *generators):
        if not all(self.contains(g) for g in generators):
            raise ValueError('generators must belong to the containing submodule')
        return SubModule(self.container, generators)


class Ideal:
    def __init__(self, ring, generators):
        from ... import DomainError, groebner
        self.ring = ring
        self.gens = tuple(ring.convert(g) for g in generators)
        self._basis = None
        if ring.is_Field:
            self._generator = ring.one if any(self.gens) else ring.zero
        elif ring.is_ZZ or getattr(ring, 'is_ZZ_I', False):
            g = ring.zero
            for a in self.gens:
                g = ring.gcd(g, a)
            self._generator = g*ring.canonical_unit(g) if getattr(ring, 'is_ZZ_I', False) else ring.abs(g)
        elif ring.is_PolynomialRing and ring.dom.is_Field:
            self._basis = groebner([ring.to_sympy(g) for g in self.gens],
                                   *ring.symbols, domain=ring.dom, order=ring.order)
            self._generator = None
        else:
            raise DomainError('ideal reduction requires ZZ, a field, or a polynomial ring over a field')

    def __repr__(self):
        return '<' + ','.join(map(str, self.gens)) + '>'

    def reduce_element(self, element):
        element = self.ring.convert(element)
        if self._basis is not None:
            remainder = self._basis.reduce(self.ring.to_sympy(element))[1]
            return self.ring.convert(remainder)
        if self._generator:
            return self.ring.rem(element, self._generator)
        return element

    def contains(self, element):
        return not self.reduce_element(element)

    __contains__ = contains

    def is_zero(self):
        return all(not g for g in self.gens)

    def is_whole_ring(self):
        return self.contains(self.ring.one)

    def subset(self, other):
        self._check(other)
        return all(other.contains(g) for g in self.gens)

    def _check(self, other):
        if not isinstance(other, Ideal) or other.ring != self.ring:
            raise ValueError('ideals must belong to the same ring')

    def __eq__(self, other):
        return (isinstance(other, Ideal) and other.ring == self.ring and
                self.subset(other) and other.subset(self))

    def union(self, other):
        self._check(other)
        return Ideal(self.ring, self.gens + other.gens)

    __add__ = union

    def product(self, other):
        self._check(other)
        return Ideal(self.ring, [a*b for a in self.gens for b in other.gens])

    __mul__ = product


class QuotientRing:
    def __init__(self, ring, ideal):
        if ideal.ring != ring:
            raise ValueError('the ideal must belong to the base ring')
        self.ring, self.base_ideal = ring, ideal

    def __repr__(self):
        return f'{self.ring}/{self.base_ideal}'

    def __eq__(self, other):
        return (isinstance(other, QuotientRing) and self.ring == other.ring and
                self.base_ideal == other.base_ideal)

    def convert(self, value):
        if isinstance(value, QuotientElement):
            if value.parent() != self:
                raise ValueError('incompatible quotient rings')
            return value
        return QuotientElement(self, self.base_ideal.reduce_element(value))

    __call__ = convert

    @property
    def zero(self):
        return self(self.ring.zero)

    @property
    def one(self):
        return self(self.ring.one)

    @property
    def gens(self):
        return tuple(self(g) for g in self.ring.gens)

    def to_sympy(self, value):
        return self.ring.to_sympy(self.convert(value).data)

    def revert(self, value):
        from ... import DomainError
        from ..polyerrors import NotReversible
        value = self.convert(value)
        if self.base_ideal.is_zero():
            return self(self.ring.revert(value.data))
        if self.base_ideal.is_whole_ring():
            raise NotReversible('the zero ring has no nonzero unit')
        if self.base_ideal._generator is not None:
            return self(self.ring.invert(value.data, self.base_ideal._generator))
        basis = list(self.base_ideal._basis)
        if len(self.ring.symbols) == 1 and len(basis) == 1:
            return self(self.ring.invert(value.data, self.ring.convert(basis[0])))
        raise DomainError('quotient inversion currently requires a principal univariate ideal')


class QuotientElement:
    _symox_quotient_element = True

    def __init__(self, ring, data):
        self._ring, self.data = ring, data

    def parent(self):
        return self._ring

    def __repr__(self):
        return f'{self.data} + {self._ring.base_ideal}'

    def __bool__(self):
        return bool(self.data)

    def __eq__(self, other):
        try:
            return self.data == self._ring.convert(other).data
        except (TypeError, ValueError):
            return False

    def __add__(self, other):
        return self._ring(self.data + self._ring.convert(other).data)

    __radd__ = __add__

    def __neg__(self):
        return self._ring(-self.data)

    def __sub__(self, other):
        return self + (-self._ring.convert(other))

    def __rsub__(self, other):
        return self._ring.convert(other) - self

    def __mul__(self, other):
        return self._ring(self.data * self._ring.convert(other).data)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * self._ring.revert(other)

    def __rtruediv__(self, other):
        return self._ring.convert(other) / self

    def __pow__(self, n):
        n = operator.index(n)
        if n < 0:
            return self._ring.revert(self) ** -n
        result, base = self._ring.one, self
        while n:
            if n & 1:
                result *= base
            base *= base
            n >>= 1
        return result
