"""Matrices of polynomials, after the reference ``sympy.polys.polymatrix``.

The reference stores a ``DomainMatrix`` over ``K[gens]``. Here the entries
are kept as expressions together with the generators and the ground domain;
items are converted to ``Poly`` in those generators when they are read, as
the reference ``__getitem__`` does. ``rref`` and ``nullspace`` work in the
ground field and require ground entries, like the reference; the reduced row
echelon form is unique, and entries are normalized by ``cancel`` so that zero
tests are exact for rational function fields.
"""
from symox.core.expr import Expr
from symox.core.sympify import _sympify

from symox.polys.polyerrors import CoercionFailed
from symox.polys.polytools import Poly, parallel_poly_from_expr
from symox.polys.domains import QQ


def _cancel(e):
    from symox import cancel
    return cancel(e)


class MutablePolyDenseMatrix:
    """
    A mutable matrix of objects from poly module or to operate with them.

    Examples
    ========

    >>> from sympy.polys.polymatrix import PolyMatrix
    >>> from sympy import Symbol, Poly
    >>> x = Symbol('x')
    >>> pm1 = PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(x**3, x), Poly(-1 + x, x)]])
    >>> v1 = PolyMatrix([[1, 0], [-1, 0]], x)
    >>> pm1*v1
    PolyMatrix([
    [    x**2 + x, 0],
    [x**3 - x + 1, 0]], ring=QQ[x])
    """

    def __new__(cls, *args, ring=None):
        if not args:
            if ring is None:
                raise TypeError("The ring needs to be specified for an empty PolyMatrix")
            rows, cols, items, gens = 0, 0, [], ()
        elif isinstance(args[0], list):
            elements, gens = args[0], args[1:]
            if not elements:
                rows, cols, items = 0, 0, []
            elif isinstance(elements[0], (list, tuple)):
                rows, cols = len(elements), len(elements[0])
                items = [e for row in elements for e in row]
            else:
                rows, cols = len(elements), 1
                items = elements
        elif [type(a) for a in args[:3]] == [int, int, list]:
            rows, cols, items, gens = args[0], args[1], args[2], args[3:]
        elif [type(a) for a in args[:3]] == [int, int, type(lambda: 0)]:
            rows, cols, func, gens = args[0], args[1], args[2], args[3:]
            items = [func(i, j) for i in range(rows) for j in range(cols)]
        else:
            raise TypeError("Invalid arguments")
        if len(gens) == 1 and isinstance(gens[0], tuple):
            gens = gens[0]
        return cls.from_list(rows, cols, items, gens, ring)

    @classmethod
    def from_list(cls, rows, cols, items, gens, ring):
        items = [item if isinstance(item, Poly) else _sympify(item) for item in items]
        if ring is not None:
            if isinstance(ring, str):
                from symox import Dummy
                ring = Poly(0, Dummy(), domain=ring).domain
            domain, gens = ring, tuple(gens)
        elif items and all(isinstance(item, Poly) for item in items):
            p = items[0]
            for p2 in items[1:]:
                p, _ = p.unify(p2)
            domain, gens = p.domain, tuple(p.gens)
        else:
            polys, info = parallel_poly_from_expr(items, gens, field=True)
            domain, gens = info['domain'], tuple(info['gens'])
        exprs = [item.as_expr() if isinstance(item, Poly) else item for item in items]
        return cls._new(rows, cols, exprs, gens, domain)

    @classmethod
    def _new(cls, rows, cols, exprs, gens, domain):
        obj = object.__new__(cls)
        obj._rows, obj._cols = rows, cols
        obj._items = list(exprs)
        obj.gens = tuple(gens)
        obj.domain = domain
        return obj

    @property
    def ring(self):
        return self.domain[self.gens] if self.gens else self.domain

    def _like(self, rows, cols, exprs, domain=None):
        return self._new(rows, cols, exprs, self.gens,
                         self.domain if domain is None else domain)

    def _unified(self, other):
        if self.domain == other.domain:
            return self.domain
        try:
            return self.domain.unify(other.domain)
        except Exception:
            return self.domain

    def to_Matrix(self):
        from symox import Matrix
        if not self._rows or not self._cols:
            from symox import zeros
            return zeros(self._rows, self._cols)
        return Matrix(self._rows, self._cols, [_cancel(e) for e in self._items])

    @classmethod
    def from_Matrix(cls, other, *gens, ring=None):
        return cls(*other.shape, list(other), *gens, ring=ring)

    def set_gens(self, gens):
        return self.from_Matrix(self.to_Matrix(), gens)

    def __repr__(self):
        if self.rows * self.cols:
            return 'Poly' + repr(self.to_Matrix())[:-1] + f', ring={self.ring})'
        else:
            return f'PolyMatrix({self.rows}, {self.cols}, [], ring={self.ring})'

    @property
    def shape(self):
        return (self._rows, self._cols)

    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    def __len__(self):
        return self.rows * self.cols

    def flat(self):
        return list(self._items)

    def _poly(self, e):
        return Poly(e, *self.gens, domain=self.domain)

    def _entry(self, i, j):
        return self._items[i*self._cols + j]

    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self._poly(item) for item in self._items[key]]
        elif isinstance(key, int):
            return self._poly(self._items[key])
        i, j = key
        if isinstance(i, int) and isinstance(j, int):
            if i < 0:
                i += self._rows
            if j < 0:
                j += self._cols
            if not (0 <= i < self._rows and 0 <= j < self._cols):
                raise IndexError("index out of range")
            return self._poly(self._entry(i, j))
        rows = range(self._rows)[i] if isinstance(i, slice) else [i % self._rows]
        cols = range(self._cols)[j] if isinstance(j, slice) else [j % self._cols]
        items = [self._entry(r, c) for r in rows for c in cols]
        return self._like(len(rows), len(cols), items)

    def __iter__(self):
        return iter(self[:])

    def __eq__(self, other):
        if not isinstance(self, type(other)):
            return NotImplemented
        return (self.shape == other.shape and self.gens == other.gens
                and all(_cancel(a - b) == 0 for a, b in zip(self._items, other._items)))

    def __add__(self, other):
        if isinstance(other, type(self)):
            if self.shape != other.shape:
                from symox.polys.matrices.exceptions import DMShapeError
                raise DMShapeError("shape")
            return self._like(self._rows, self._cols,
                [a + b for a, b in zip(self._items, other._items)], self._unified(other))
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, type(self)):
            if self.shape != other.shape:
                from symox.polys.matrices.exceptions import DMShapeError
                raise DMShapeError("shape")
            return self._like(self._rows, self._cols,
                [a - b for a, b in zip(self._items, other._items)], self._unified(other))
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, type(self)):
            if self._cols != other._rows:
                from symox.polys.matrices.exceptions import DMShapeError
                raise DMShapeError("shape")
            n = self._cols
            items = []
            for i in range(self._rows):
                for j in range(other._cols):
                    items.append(sum((self._entry(i, k)*other._entry(k, j) for k in range(n)),
                                     _sympify(0)))
            return self._new(self._rows, other._cols, [_cancel(e) for e in items],
                             self.gens, self._unified(other))
        elif isinstance(other, int):
            other = _sympify(other)
        if isinstance(other, Expr):
            return self._like(self._rows, self._cols, [_cancel(e*other) for e in self._items])
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, int):
            other = _sympify(other)
        if isinstance(other, Expr):
            return self._like(self._rows, self._cols, [_cancel(other*e) for e in self._items])
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Poly):
            other = other.as_expr()
        elif isinstance(other, int):
            other = _sympify(other)
        if not isinstance(other, Expr):
            return NotImplemented
        if other.is_zero:
            raise ZeroDivisionError("division by zero")
        return self._like(self._rows, self._cols, [_cancel(e/other) for e in self._items])

    def __neg__(self):
        return self._like(self._rows, self._cols, [-e for e in self._items])

    def transpose(self):
        items = [self._entry(i, j) for j in range(self._cols) for i in range(self._rows)]
        return self._like(self._cols, self._rows, items)

    def row_join(self, other):
        if self._rows != other._rows:
            from symox.polys.matrices.exceptions import DMShapeError
            raise DMShapeError("Expected the same number of rows")
        items = []
        for i in range(self._rows):
            items.extend(self._items[i*self._cols:(i+1)*self._cols])
            items.extend(other._items[i*other._cols:(i+1)*other._cols])
        return self._like(self._rows, self._cols + other._cols, items, self._unified(other))

    def col_join(self, other):
        if self._cols != other._cols:
            from symox.polys.matrices.exceptions import DMShapeError
            raise DMShapeError("Expected the same number of columns")
        return self._like(self._rows + other._rows, self._cols,
                          self._items + other._items, self._unified(other))

    def applyfunc(self, func):
        M = self.to_Matrix().applyfunc(func)
        return self.from_Matrix(M, self.gens)

    @classmethod
    def eye(cls, n, gens):
        if not isinstance(gens, tuple):
            gens = (gens,)
        items = [_sympify(1 if i == j else 0) for i in range(n) for j in range(n)]
        return cls._new(n, n, items, gens, QQ)

    @classmethod
    def zeros(cls, m, n, gens):
        if not isinstance(gens, tuple):
            gens = (gens,)
        return cls._new(m, n, [_sympify(0)]*(m*n), gens, QQ)

    def _ground_check(self, name):
        if not (self.domain.is_Field and all(p.is_ground for p in self)):
            raise ValueError("PolyMatrix %s is only for ground field elements" % name)

    def _rref_rows(self):
        rows = [[_cancel(self._entry(i, j)) for j in range(self._cols)] for i in range(self._rows)]
        pivots = []
        r = 0
        for j in range(self._cols):
            if r == self._rows:
                break
            pivot = None
            for i in range(r, self._rows):
                if rows[i][j] != 0:
                    pivot = i
                    break
            if pivot is None:
                continue
            rows[r], rows[pivot] = rows[pivot], rows[r]
            p = rows[r][j]
            rows[r] = [_cancel(e/p) for e in rows[r]]
            for i in range(self._rows):
                if i != r and rows[i][j] != 0:
                    c = rows[i][j]
                    rows[i] = [_cancel(a - c*b) for a, b in zip(rows[i], rows[r])]
            pivots.append(j)
            r += 1
        return rows, tuple(pivots)

    def rref(self, simplify='ignore', normalize_last='ignore'):
        # If this is K[x] then computes RREF in ground field K.
        self._ground_check('rref')
        rows, pivots = self._rref_rows()
        return self._like(self._rows, self._cols, [e for row in rows for e in row]), pivots

    def nullspace(self):
        # If this is K[x] then computes nullspace in ground field K.
        self._ground_check('nullspace')
        rows, pivots = self._rref_rows()
        basis = []
        for j in sorted(set(range(self._cols)) - set(pivots)):
            vec = [_sympify(0)]*self._cols
            vec[j] = _sympify(1)
            for i, p in enumerate(pivots):
                if rows[i][j] != 0:
                    vec[p] = -rows[i][j]
            basis.append(self._like(self._cols, 1, vec))
        return basis

    def rank(self):
        return self.cols - len(self.nullspace())


MutablePolyMatrix = PolyMatrix = MutablePolyDenseMatrix
