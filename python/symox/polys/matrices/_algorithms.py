"""Linear algebra over native typed coefficients, without expression coercion."""

from .exceptions import (DMBadInputError, DMNotAField,
                         DMNonSquareMatrixError, DMNonInvertibleMatrixError,
                         DMShapeError)


class DomainMatrixAlgorithms:
    def _square(self):
        if self.shape[0] != self.shape[1]:
            raise DMNonSquareMatrixError('Matrix must be square')
        return self.shape[0]

    def rref(self, *, method='auto'):
        if method not in ('auto', 'GJ', 'FF'):
            raise DMBadInputError("Supported RREF methods: 'auto', 'GJ', 'FF'")
        if method == 'FF' or (method == 'auto' and not self.domain.is_Field):
            numerator, denominator, pivots = self.rref_den(method='FF')
            return numerator.to_field() / denominator, pivots
        matrix = self.to_field()
        domain, rows = matrix.domain, matrix.to_list()
        m, n = self.shape
        pivots = []
        for column in range(n):
            row = len(pivots)
            pivot = next((i for i in range(row, m)
                          if not domain.is_zero(rows[i][column])), None)
            if pivot is None:
                continue
            rows[row], rows[pivot] = rows[pivot], rows[row]
            divisor = rows[row][column]
            rows[row] = [domain.exquo(v, divisor) for v in rows[row]]
            for i in range(m):
                if i == row or domain.is_zero(rows[i][column]):
                    continue
                coefficient = rows[i][column]
                rows[i] = [v - coefficient*w for v, w in zip(rows[i], rows[row])]
                rows[i][column] = domain.zero
            pivots.append(column)
        return type(self)(rows, self.shape, domain, fmt=self._fmt), tuple(pivots)

    def rref_den(self, *, method='auto', keep_domain=True):
        if method not in ('auto', 'GJ', 'FF'):
            raise DMBadInputError("Supported RREF methods: 'auto', 'GJ', 'FF'")
        if method == 'GJ' or (method == 'auto' and self.domain.is_Field):
            matrix, pivots = self.rref(method='GJ')
            if not keep_domain:
                return matrix, matrix.domain.one, pivots
            return matrix._denominator_in(self.domain, pivots)
        domain, rows = self.domain, self.to_list()
        m, n = self.shape
        pivots, denominator = [], domain.one
        # Fraction-free Gauss-Jordan: every pivot column has the common
        # denominator on its diagonal; divisions are exact in the input ring.
        for column in range(n):
            row = len(pivots)
            pivot = next((i for i in range(row, m)
                          if not domain.is_zero(rows[i][column])), None)
            if pivot is None:
                continue
            rows[row], rows[pivot] = rows[pivot], rows[row]
            value = rows[row][column]
            for i in range(m):
                if i == row:
                    continue
                factor = rows[i][column]
                rows[i] = [domain.exquo(value*a - factor*b, denominator)
                           if j != column else domain.zero
                           for j, (a, b) in enumerate(zip(rows[i], rows[row]))]
            denominator = value
            pivots.append(column)
        result = type(self)(rows, self.shape, domain, fmt=self._fmt)
        return result, denominator, tuple(pivots)

    def _denominator_in(self, target, extra=None):
        if self.domain == target:
            result = self.copy(), target.one
        else:
            denominator = target.one
            for row in self._data.values():
                for value in row.values():
                    denominator = target.lcm(
                        denominator, target.convert(self.domain.denom(value)))
            scaled = self.mul(self.domain.convert(denominator)).convert_to(target)
            result = scaled, denominator
        return result if extra is None else (*result, extra)

    def rank(self):
        return len(self.rref()[1])

    def nullspace(self, divide_last=False):
        if divide_last:
            if not self.domain.is_Field:
                raise DMNotAField('Cannot normalize vectors over a non-field')
            reduced, pivots = self.rref()
        else:
            reduced, denominator, pivots = self.rref_den()
            reduced = reduced.mul(self.domain.canonical_unit(denominator))
        return reduced.nullspace_from_rref(pivots)

    def nullspace_from_rref(self, pivots=None):
        domain = self.domain
        if pivots is None:
            pivots = tuple(min(row) for _, row in sorted(self._data.items()) if row)
        pivots = tuple(pivots)
        if len(pivots) > self.shape[0] or any(not 0 <= p < self.shape[1] for p in pivots):
            raise DMBadInputError('Invalid RREF pivot columns')
        free = [j for j in range(self.shape[1]) if j not in pivots]
        scale = self._data.get(0, {}).get(pivots[0], domain.one) if pivots else domain.one
        rows = {}
        for i, column in enumerate(free):
            row = {column: scale}
            for k, pivot in enumerate(pivots):
                value = -self._data.get(k, {}).get(column, domain.zero)
                if not domain.is_zero(value):
                    row[pivot] = value
            rows[i] = row
        return type(self)(rows, (len(free), self.shape[1]), domain, fmt=self._fmt)

    def rowspace(self):
        if not self.domain.is_Field:
            raise DMNotAField('Rowspace requires a field')
        # Select independent original rows, including when leading rows vanish.
        _, pivots = self.transpose().rref()
        return self.extract(pivots, range(self.shape[1]))

    def columnspace(self):
        if not self.domain.is_Field:
            raise DMNotAField('Columnspace requires a field')
        _, pivots = self.rref()
        return self.extract(range(self.shape[0]), pivots)

    def inv(self):
        n = self._square()
        if not self.domain.is_Field:
            raise DMNotAField('Matrix inverse requires a field')
        reduced, pivots = self.hstack(self.eye(n, self.domain)).rref()
        if pivots[:n] != tuple(range(n)):
            raise DMNonInvertibleMatrixError('Matrix is singular')
        return reduced.extract(range(n), range(n, 2*n))

    def inv_den(self, method=None):
        n = self._square()
        return self.solve_den(self.eye(n, self.domain), method=method)

    def solve_den(self, b, method=None):
        if method not in (None, 'rref', 'charpoly'):
            raise DMBadInputError("method must be 'rref' or 'charpoly'")
        self._check(b, product=False, rows_only=True)
        m, n = self.shape
        if method == 'charpoly':
            adjugate, determinant = self.adj_det()
            if self.domain.is_zero(determinant):
                raise DMNonInvertibleMatrixError('Matrix is singular')
            return adjugate.matmul(b), determinant
        if m < n:
            raise DMNonInvertibleMatrixError('Underdetermined system')
        reduced, denominator, pivots = self.hstack(b).rref_den()
        if pivots != tuple(range(n)):
            raise DMNonInvertibleMatrixError('System has no unique solution')
        return reduced.extract(range(n), range(n, n+b.shape[1])), denominator

    def charpoly(self):
        """Division-free Samuelson-Berkowitz coefficients in descending order."""
        n = self._square()
        domain, rows = self.domain, self.to_list()
        coefficients = [domain.one]
        # Grow trailing principal submatrices. The Toeplitz column consists
        # of 1, -a, -R*C, -R*A*C, ...; convolve it with charpoly(A).
        for start in range(n-1, -1, -1):
            size = n-start
            column = [domain.one, -rows[start][start]]
            vector = [rows[i][start] for i in range(start+1, n)]
            for power in range(size-1):
                column.append(-domain.sum(rows[start][j]*vector[j-start-1]
                                          for j in range(start+1, n)))
                if power < size-2:
                    vector = [domain.sum(rows[i][j]*vector[j-start-1]
                                         for j in range(start+1, n))
                              for i in range(start+1, n)]
            coefficients = [domain.sum(column[i-j]*coefficients[j]
                                       for j in range(min(i, size-1)+1))
                            for i in range(size+1)]
        return coefficients

    def adj_det(self):
        n = self._square()
        identity = self.eye(n, self.domain)
        if not n:
            return identity, self.domain.one
        coefficients = self.charpoly()
        adjugate = identity
        for coefficient in coefficients[1:-1]:
            adjugate = self.matmul(adjugate).add(identity.mul(coefficient))
        sign = self.domain(-1 if (n-1) % 2 else 1)
        return adjugate.mul(sign), -sign*coefficients[-1]

    def det(self):
        n = self._square()
        if not n:
            return self.domain.one
        domain, rows = self.domain, self.to_list()
        denominator, sign = domain.one, domain.one
        for k in range(n-1):
            pivot = next((i for i in range(k, n)
                          if not domain.is_zero(rows[i][k])), None)
            if pivot is None:
                return domain.zero
            if pivot != k:
                rows[k], rows[pivot] = rows[pivot], rows[k]
                sign = -sign
            value = rows[k][k]
            for i in range(k+1, n):
                for j in range(k+1, n):
                    rows[i][j] = domain.exquo(
                        value*rows[i][j] - rows[i][k]*rows[k][j], denominator)
                rows[i][k] = domain.zero
            denominator = value
        return sign*rows[-1][-1]
