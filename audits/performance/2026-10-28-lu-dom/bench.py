"""Before/after timings of public LU / Gauss-Jordan calls beyond QQ.

usage: python bench.py CASE first|repeat  -> prints one JSON line
       python bench.py list
Each case returns a canonical string of its result; the driver checks that
both binaries give identical strings (a wrong result is never a speed-up).
"""
import json, sys, time, random


def build(case):
    from symox import GF, QQ_I, QQ, sqrt, I, Matrix, Rational
    from symox.polys.matrices import DomainMatrix
    rng = random.Random(20261028)
    def dm(K, m, n, gen, sparse=False):
        M = DomainMatrix([[K(gen()) for _ in range(n)] for _ in range(m)], (m, n), K)
        return M.to_sparse() if sparse else M.to_dense()
    F = GF(10007)
    A2 = QQ.algebraic_field(sqrt(2))
    def gf():
        return rng.randrange(10007)
    def qi():
        return QQ_I(rng.randint(-9, 9), rng.randint(-9, 9))
    def a2():
        return A2.from_sympy(rng.randint(-9, 9) + rng.randint(-9, 9) * sqrt(2))
    cases = {
        'gf_lu_12': lambda: (lambda M: lambda: M.lu())(dm(F, 12, 12, gf)),
        'gf_lu_solve_12': lambda: (lambda M, b: lambda: M.lu_solve(b))(dm(F, 12, 12, gf), dm(F, 12, 1, gf)),
        'gf_rref_dense_12': lambda: (lambda M: lambda: M.rref(method='GJ_dense'))(dm(F, 12, 12, gf)),
        'gf_rref_sparse_12': lambda: (lambda M: lambda: M.rref(method='GJ'))(dm(F, 12, 12, gf, True)),
        'alg_lu_6': lambda: (lambda M: lambda: M.lu())(dm(A2, 6, 6, a2)),
        'alg_lu_solve_6': lambda: (lambda M, b: lambda: M.lu_solve(b))(dm(A2, 6, 6, a2), dm(A2, 6, 2, a2)),
        'alg_rref_dense_6': lambda: (lambda M: lambda: M.rref(method='GJ_dense'))(dm(A2, 6, 6, a2)),
        'alg_rref_sparse_6': lambda: (lambda M: lambda: M.rref(method='GJ'))(dm(A2, 6, 6, a2, True)),
        'qqi_rref_dense_8': lambda: (lambda M: lambda: M.rref(method='GJ_dense'))(dm(QQ_I, 8, 8, qi)),
        'qqi_rref_sparse_8': lambda: (lambda M: lambda: M.rref(method='GJ'))(dm(QQ_I, 8, 8, qi, True)),
        'qqi_lu_8_fallback': lambda: (lambda M: lambda: M.lu())(dm(QQ_I, 8, 8, qi)),
        'gf_lu_3_small': lambda: (lambda M: lambda: M.lu())(dm(F, 3, 3, gf)),
        'matrix_LUsolve_sqrt2_4': lambda: (lambda M, b: lambda: M.LUsolve(b))(
            Matrix(4, 4, lambda i, j: rng.randint(-5, 5) + rng.randint(-5, 5) * sqrt(2)), Matrix(4, 1, lambda i, j: rng.randint(-5, 5))),
        'matrix_LUsolve_I_4': lambda: (lambda M, b: lambda: M.LUsolve(b))(
            Matrix(4, 4, lambda i, j: rng.randint(-5, 5) + rng.randint(-5, 5) * I), Matrix(4, 1, lambda i, j: rng.randint(-5, 5))),
    }
    if case == 'list':
        return list(cases)
    return cases[case]()


def canon(result):
    if isinstance(result, tuple):
        return '|'.join(canon(r) for r in result)
    if hasattr(result, 'to_Matrix'):
        return repr(result.domain) + str(result.to_Matrix())
    return str(result)


if __name__ == '__main__':
    case = sys.argv[1]
    if case == 'list':
        print(' '.join(build('list'))); sys.exit()
    mode = sys.argv[2]
    call = build(case)
    if mode == 'first':
        t = time.perf_counter(); r = call(); dt = time.perf_counter() - t
        print(json.dumps({'case': case, 'mode': mode, 'seconds': dt, 'result': canon(r)}))
    else:
        call()
        times = []
        for _ in range(15):
            t = time.perf_counter(); r = call(); times.append(time.perf_counter() - t)
        times.sort()
        print(json.dumps({'case': case, 'mode': mode, 'seconds': times[len(times) // 2], 'result': canon(r)}))
