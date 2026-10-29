"""Array-expression and matrix-derivative timings.

usage: python bench.py  (PYTHONPATH selects the symox tree; set it to the
python/ directory of the tree to measure). Prints one JSON line:
{case: [first_call_s, median_repeated_s]}. First call includes the lazy
import of the converter modules for the first case only.
"""
import json
import statistics
import sys
import time

REPEAT = 10


def cases():
    from symox import MatrixSymbol, Sum, Trace, sqrt, symbols
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    x = MatrixSymbol('x', 3, 1)
    y = MatrixSymbol('y', 3, 1)
    d, i, j = symbols('d i j')
    M = MatrixSymbol('M', d, d)
    N = MatrixSymbol('N', d, d)

    def conv_m2a():
        from symox.tensor.array.expressions import convert_matrix_to_array
        return convert_matrix_to_array(A*B*A.T)

    def conv_a2m():
        from symox.tensor.array.expressions import (convert_array_to_matrix,
            convert_matrix_to_array)
        return convert_array_to_matrix(convert_matrix_to_array(A*B*A.T))

    def conv_i2a():
        from symox.tensor.array.expressions import convert_indexed_to_array
        return convert_indexed_to_array(Sum(M[i, j]*N[j, i], (i, 0, d - 1), (j, 0, d - 1)))

    new_api = []
    try:  # names added by this task (absent from the base tree)
        from symox import Array
        from symox.tensor.array.arrayop import Flatten
        from symox.matrices.expressions.determinant import det
        big = Array(range(600)).reshape(20, 30)
        new_api = [('list(Flatten(20x30 array))', lambda: list(Flatten(big))),
                   ('determinant.det(A)', lambda: det(A))]
    except ImportError:
        pass

    return new_api + [
        ('A.diff(A)', lambda: A.diff(A)),
        ('(A*A).diff(A)', lambda: (A*A).diff(A)),
        ('Trace(A*B).diff(A)', lambda: Trace(A*B).diff(A)),
        ('(x.T*x).diff(x)', lambda: (x.T*x).diff(x)),
        ('sqrt(x.T*y).diff(x)', lambda: sqrt(x.T*y).diff(x)),
        ('convert_matrix_to_array', conv_m2a),
        ('convert_array_to_matrix', conv_a2m),
        ('convert_indexed_to_array', conv_i2a),
    ]


def main():
    import symox  # noqa: F401  (import cost excluded)
    result = {}
    for name, f in cases():
        t = time.perf_counter(); f(); first = time.perf_counter() - t
        times = []
        for _ in range(REPEAT):
            t = time.perf_counter(); f(); times.append(time.perf_counter() - t)
        result[name] = [first, statistics.median(times)]
    json.dump(result, sys.stdout)
    print()


if __name__ == '__main__':
    main()
