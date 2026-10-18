"""Domain matrices backed by native coefficient-domain elements.

This conversion-oriented subset supports implemented exact and numerical
domains. Algebraic field inference is opt-in; general expressions remain
formal and numerical fields retain the source floating-point precision.
"""
import operator

from .exceptions import DMBadInputError, DMDomainError, DMShapeError
from ._algorithms import DomainMatrixAlgorithms


from ..constructor import _infer_domain


class DomainMatrix(DomainMatrixAlgorithms):
    """A matrix of native elements of one coefficient domain.

    Construction accepts a rectangular list of rows or a dictionary of
    dictionaries, with an explicit shape and native Domain. Conversion does
    not retain aliases to mutable matrices or input containers.
    """
    _symox_domainmatrix = True

    def __init__(self, rows, shape, domain, *, fmt=None):
        from ... import Domain
        if not isinstance(domain, Domain):
            raise TypeError('domain must be a native symox Domain')
        if len(shape) != 2:
            raise DMBadInputError('shape must contain two dimensions')
        shape = tuple(operator.index(n) for n in shape)
        if min(shape) < 0:
            raise DMBadInputError('shape must be nonnegative')
        if fmt not in (None, 'dense', 'sparse'):
            raise ValueError("fmt must be 'dense' or 'sparse'")
        self.shape, self.domain = shape, domain
        self._fmt = fmt or ('sparse' if isinstance(rows, dict) else 'dense')
        data = {}
        if isinstance(rows, dict):
            entries = ((i, j, value) for i, row in rows.items() for j, value in row.items())
        else:
            rows = [list(row) for row in rows]
            if len(rows) != shape[0] or any(len(row) != shape[1] for row in rows):
                raise DMBadInputError('Input rows do not match shape')
            entries = ((i, j, value) for i, row in enumerate(rows) for j, value in enumerate(row))
        for i, j, value in entries:
            i, j = operator.index(i), operator.index(j)
            if not (0 <= i < shape[0] and 0 <= j < shape[1]):
                raise DMBadInputError('Matrix index out of range')
            value = domain.convert(value)
            if not domain.is_zero(value):
                data.setdefault(i, {})[j] = value
        self._data = data

    @classmethod
    def from_Matrix(cls, matrix, fmt='sparse', **options):
        from ... import sympify
        values = {key: sympify(v) for key, v in matrix.todok().items()}
        inference_values = (list(values.values()) if getattr(matrix, '_symox_sparse', False)
                            else [v for row in matrix.tolist() for v in row])
        domain = _infer_domain(inference_values, **options)
        rows = {}
        for (i, j), value in values.items():
            rows.setdefault(i, {})[j] = value
        return cls(rows, matrix.shape, domain, fmt=fmt)

    def to_Matrix(self):
        from ... import Matrix
        return Matrix.from_dok(*self.shape, {
            (i, j): self.domain.to_sympy(value)
            for i, row in self._data.items() for j, value in row.items()})

    def to_dod(self):
        return {i: dict(row) for i, row in self._data.items()}

    def to_dok(self):
        return {(i, j): value for i, row in self._data.items() for j, value in row.items()}

    def to_list(self):
        zero = self.domain.zero
        return [[self._data.get(i, {}).get(j, zero) for j in range(self.shape[1])]
                for i in range(self.shape[0])]

    def to_list_flat(self):
        return [v for row in self.to_list() for v in row]

    def convert_to(self, domain):
        return type(self)(self._data, self.shape, domain, fmt=self._fmt)

    def to_field(self):
        return self.convert_to(self.domain.get_field())

    def choose_domain(self, **options):
        return type(self).from_Matrix(self.to_Matrix(), fmt=self._fmt, **options)

    def copy(self):
        return self.convert_to(self.domain)

    def to_dense(self):
        return type(self)(self._data, self.shape, self.domain, fmt='dense')

    def to_sparse(self):
        return type(self)(self._data, self.shape, self.domain, fmt='sparse')

    def __eq__(self, other):
        if not isinstance(other, DomainMatrix):
            return NotImplemented
        return (self.shape == other.shape and self.domain == other.domain
                and self._data == other._data)

    def __repr__(self):
        rows = self.to_list() if self._fmt == 'dense' else self.to_dod()
        return 'DomainMatrix(%r, %r, %s)' % (rows, self.shape, self.domain)

    def _check(self, other, *, product=False, rows_only=False):
        if not isinstance(other, DomainMatrix):
            raise TypeError('Expected a DomainMatrix')
        if self.domain != other.domain:
            raise DMDomainError('Domain mismatch')
        if (self.shape[0] != other.shape[0] if rows_only else
                self.shape[1] != other.shape[0] if product else self.shape != other.shape):
            raise DMShapeError('Shape mismatch')

    def add(self, other):
        self._check(other)
        rows = self.to_dod()
        for i, row in other._data.items():
            output = rows.setdefault(i, {})
            for j, value in row.items():
                output[j] = output.get(j, self.domain.zero) + value
        return type(self)(rows, self.shape, self.domain, fmt=self._fmt)

    def sub(self, other):
        self._check(other)
        return self.add(-other)

    def neg(self):
        return self.mul(-self.domain.one)

    __neg__ = neg

    def mul(self, scalar):
        scalar = self.domain.convert(scalar)
        return type(self)({i: {j: value*scalar for j, value in row.items()}
                           for i, row in self._data.items()},
                          self.shape, self.domain, fmt=self._fmt)

    scalarmul = mul

    def __add__(self, other):
        if not isinstance(other, DomainMatrix):
            return NotImplemented
        first, second = self.unify(other)
        return first.add(second)

    def __sub__(self, other):
        if not isinstance(other, DomainMatrix):
            return NotImplemented
        first, second = self.unify(other)
        return first.sub(second)

    def __mul__(self, other):
        if isinstance(other, DomainMatrix):
            first, second = self.unify(other)
            return first.matmul(second)
        from .domainscalar import DomainScalar
        if isinstance(other, DomainScalar):
            domain = self.domain.unify(other.domain)
            return self.convert_to(domain).mul(domain.convert(other.element))
        try:
            return self.mul(other)
        except (TypeError, ValueError):
            return NotImplemented

    __rmul__ = __mul__
    __matmul__ = __mul__

    def __truediv__(self, scalar):
        from .domainscalar import DomainScalar
        if isinstance(scalar, DomainScalar):
            scalar = scalar.element
        if not self.domain.is_Field:
            from .exceptions import DMNotAField
            raise DMNotAField('Matrix division requires a field')
        scalar = self.domain.convert(scalar)
        return self.mul(self.domain.exquo(self.domain.one, scalar))

    def __pow__(self, exponent):
        exponent = operator.index(exponent)
        n = self._square()
        base = self
        if exponent < 0:
            base, exponent = self.inv(), -exponent
        result = self.eye(n, self.domain)
        while exponent:
            if exponent & 1:
                result = result.matmul(base)
            exponent >>= 1
            if exponent:
                base = base.matmul(base)
        return result

    def matmul(self, other):
        self._check(other, product=True)
        rows = {}
        for i, row in self._data.items():
            output = rows.setdefault(i, {})
            for k, value in row.items():
                for j, factor in other._data.get(k, {}).items():
                    output[j] = output.get(j, self.domain.zero) + value*factor
        return type(self)(rows, (self.shape[0], other.shape[1]), self.domain, fmt=self._fmt)

    def transpose(self):
        rows = {}
        for i, row in self._data.items():
            for j, value in row.items():
                rows.setdefault(j, {})[i] = value
        return type(self)(rows, self.shape[::-1], self.domain, fmt=self._fmt)

    T = property(transpose)

    def extract(self, rowslist, colslist):
        def indices(values, size):
            output = []
            for value in values:
                value = operator.index(value)
                if not -size <= value < size:
                    raise IndexError('Matrix index out of range')
                output.append(value % size)
            return output
        rowslist = indices(rowslist, self.shape[0])
        colslist = indices(colslist, self.shape[1])
        rows = {i: {j: self._data[row][col] for j, col in enumerate(colslist)
                     if col in self._data.get(row, {})}
                for i, row in enumerate(rowslist)}
        return type(self)(rows, (len(rowslist), len(colslist)), self.domain, fmt=self._fmt)

    def __getitem__(self, key):
        from .domainscalar import DomainScalar
        i, j = key
        if isinstance(i, slice) or isinstance(j, slice):
            rows = range(*i.indices(self.shape[0])) if isinstance(i, slice) else [i]
            cols = range(*j.indices(self.shape[1])) if isinstance(j, slice) else [j]
            return self.extract(rows, cols)
        i, j = operator.index(i), operator.index(j)
        if not -self.shape[0] <= i < self.shape[0] or not -self.shape[1] <= j < self.shape[1]:
            raise IndexError('Matrix index out of range')
        return DomainScalar(self._data.get(i % self.shape[0], {}).get(
            j % self.shape[1], self.domain.zero), self.domain)

    @classmethod
    def zeros(cls, shape, domain, *, fmt='sparse'):
        return cls({}, shape, domain, fmt=fmt)

    @classmethod
    def ones(cls, shape, domain):
        return cls([[domain.one]*shape[1] for _ in range(shape[0])], shape, domain)

    @classmethod
    def eye(cls, shape, domain):
        if isinstance(shape, int):
            shape = (shape, shape)
        return cls({i: {i: domain.one} for i in range(min(shape))}, shape, domain)

    @classmethod
    def diag(cls, diagonal, domain, shape=None):
        diagonal = list(diagonal)
        shape = shape or (len(diagonal), len(diagonal))
        return cls({i: {i: value} for i, value in enumerate(diagonal)}, shape, domain)

    @classmethod
    def from_list(cls, rows, domain):
        rows = [list(row) for row in rows]
        return cls(rows, (len(rows), len(rows[0]) if rows else 0), domain)

    def unify(self, *others, fmt=None):
        matrices = (self,) + others
        domain = self.domain
        for matrix in others:
            domain = domain.unify(matrix.domain)
        if fmt not in (None, 'dense', 'sparse'):
            raise ValueError("fmt must be 'dense' or 'sparse'")
        return tuple(type(self)(m._data, m.shape, domain, fmt=fmt or m._fmt)
                     for m in matrices)

    def hstack(self, *others):
        rows, width = self.to_dod(), self.shape[1]
        for matrix in others:
            self._check(matrix, rows_only=True)
            for i, row in matrix._data.items():
                rows.setdefault(i, {}).update({j+width: v for j, v in row.items()})
            width += matrix.shape[1]
        return type(self)(rows, (self.shape[0], width), self.domain, fmt=self._fmt)

    def vstack(self, *others):
        matrices = (self,) + others
        if any(m.domain != self.domain for m in others):
            raise DMDomainError('Domain mismatch')
        if any(m.shape[1] != self.shape[1] for m in others):
            raise DMShapeError('Column count mismatch')
        rows, offset = {}, 0
        for matrix in matrices:
            rows.update({i+offset: dict(row) for i, row in matrix._data.items()})
            offset += matrix.shape[0]
        return type(self)(rows, (offset, self.shape[1]), self.domain, fmt=self._fmt)

    def clear_denoms(self, convert=False):
        from ... import DomainError
        from .domainscalar import DomainScalar
        try:
            ring = self.domain.get_ring()
        except DomainError:
            ring = self.domain
        numerator, denominator = self._denominator_in(ring)
        if not convert:
            numerator = numerator.convert_to(self.domain)
        return DomainScalar(denominator, ring), numerator

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.polys.matrices', 'DM')
