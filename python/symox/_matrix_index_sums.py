"""Recognize matrix chains in explicitly contracted index expressions."""
from collections import Counter


def from_index_summation(expr, first_index=None, last_index=None, dimensions=None):
    from . import (sympify, MatrixElement, MatrixExpr, MatMul, MatAdd,
                   Trace, HadamardProduct, S)
    expr = sympify(expr)
    dimensions = {} if dimensions is None else dict(dimensions)

    def same(a, b):
        return a == b or (a-b).expand() == 0

    def parse(expression, sums, start, end):
        if expression.is_Sum:
            inner, *limits = expression.args
            sums = dict(sums)
            for limit in limits:
                limit = limit if isinstance(limit, (tuple, list)) else limit.args
                if len(limit) != 3:
                    raise ValueError('matrix contractions require finite summation bounds')
                index, lower, upper = limit
                if index in sums:
                    raise ValueError('a contraction index cannot be bound twice')
                sums[index] = (lower, upper)
            return parse(inner, sums, start, end)
        if expression.is_Add:
            terms = [parse(term, sums, start, end) for term in expression.args]
            if all(isinstance(term, MatrixExpr) for term in terms):
                return MatAdd(*terms).doit()
            if any(isinstance(term, MatrixExpr) for term in terms):
                raise ValueError('an indexed sum cannot mix scalar and matrix results')
            return sum(terms, S.Zero)

        factors = expression.args if expression.is_Mul else (expression,)
        edges, scalar = [], S.One
        sizes = dict(dimensions)
        for factor in factors:
            if isinstance(factor, MatrixElement):
                matrix, row, column = factor.args
                if not row.is_Symbol or not column.is_Symbol:
                    raise NotImplementedError('matrix-chain inference requires symbolic indices')
                for index, size in [(row, matrix.rows), (column, matrix.cols)]:
                    if index in sizes and not same(sizes[index], size):
                        raise ValueError('inconsistent dimensions for a contracted index')
                    sizes[index] = size
                edges.append((matrix, row, column))
            elif factor.is_Pow and isinstance(factor.args[0], MatrixElement):
                base, power = factor.args
                if power.is_Integer and 0 < int(power) <= 64:
                    matrix, row, column = base.args
                    for index, size in [(row, matrix.rows), (column, matrix.cols)]:
                        if index in sizes and not same(sizes[index], size):
                            raise ValueError('inconsistent index dimensions')
                        sizes[index] = size
                    edges.extend([(matrix, row, column)]*int(power))
                else:
                    raise NotImplementedError('indexed powers require a positive integer exponent')
            else:
                scalar *= factor
        if not edges:
            if sums:
                from . import Sum
                return Sum(scalar, *((index, *bounds) for index, bounds in sums.items())).doit()
            return scalar
        indices = {index for _, row, column in edges for index in (row, column)}
        if scalar.has(*indices):
            raise NotImplementedError('index-dependent scalar factors require array inference')
        for index, (lower, upper) in sums.items():
            if index not in sizes:
                # An unused finite index contributes its number of values.
                scalar *= upper-lower+1
            elif not same(lower, S.Zero) or not same(upper, sizes[index]-1):
                raise ValueError('matrix contractions must sum over the entire index dimension')

        free = indices.difference(sums)
        if start is not None and start not in free:
            raise ValueError('first_index must be a free index')
        if end is not None and end not in free:
            raise ValueError('last_index must be a free index')
        if not sums and len(free) == 2:
            # Parallel free edges represent an elementwise product.
            row, column = sorted(free, key=lambda index: index.sort_key())
            if start is not None:
                row, column = start, next(i for i in free if i != start)
            elif end is not None:
                column, row = end, next(i for i in free if i != end)
            matrices = []
            for matrix, i, j in edges:
                if (i, j) == (row, column):
                    matrices.append(matrix)
                elif (j, i) == (row, column):
                    matrices.append(matrix.T)
                else:
                    break
            else:
                return scalar*(matrices[0] if len(matrices) == 1 else HadamardProduct(*matrices).doit())

        degree = Counter(index for _, i, j in edges for index in (i, j))
        if any(degree[index] != 2 for index in sums if index in indices):
            raise NotImplementedError('a matrix-chain contraction requires each bound index twice')
        if any(degree[index] != 1 for index in free) or len(free) not in (0, 2):
            raise NotImplementedError('the indexed expression has an unsupported array rank or diagonal')

        remaining, matrix_result = list(edges), None
        # Process the open chain first; closed connected components become
        # scalar traces and commute with the resulting matrix expression.
        components = []
        if free:
            first = start
            if first is None:
                first = next(i for i in free if i != end) if end is not None else min(free, key=lambda i: i.sort_key())
            components.append((first, False))
        while remaining:
            first, closed = components.pop(0) if components else (remaining[0][1], True)
            current, matrices = first, []
            while True:
                selected = next((position for position, (_, i, j) in enumerate(remaining)
                                 if current in (i, j)), None)
                if selected is None:
                    break
                matrix, row, column = remaining.pop(selected)
                if row == current:
                    matrices.append(matrix)
                    current = column
                else:
                    matrices.append(matrix.T)
                    current = row
                if current == first:
                    break
            product = MatMul(*matrices).doit()
            if closed:
                if current != first:
                    raise ValueError('disconnected open index chain')
                scalar *= Trace(product)
            else:
                if current not in free or current == first or (end is not None and current != end):
                    raise ValueError('index chain does not end at the requested free index')
                matrix_result = product
        return scalar if matrix_result is None else scalar*matrix_result

    return parse(expr, {}, first_index, last_index)
