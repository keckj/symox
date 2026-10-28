"""One fresh-process measurement: python bench_case.py <case> -> JSON line."""
import json, random, statistics, sys, time

case = sys.argv[1]
from symox import ZZ, QQ, GF, symbols
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.normalforms import smith_normal_form, smith_normal_decomp, invariant_factors

x = symbols('x')


def rnd(r, c, lo, hi, seed):
    g = random.Random(seed)
    return [[g.randint(lo, hi) for _ in range(c)] for _ in range(r)]


def dm(K, rows):
    return DomainMatrix([[K(v) for v in r] for r in rows], (len(rows), len(rows[0])), K)


def qq(rows, seed):
    g = random.Random(seed)
    return DomainMatrix([[QQ(v, g.randint(1, 5)) for v in r] for r in rows], (len(rows), len(rows[0])), QQ)


def charmat(n, seed):
    # x*I - A over QQ[x], A with small integer entries and a repeated block.
    R = QQ[x]
    A = rnd(n, n, -2, 2, seed)
    return DomainMatrix([[R(x if i == j else 0) - R(A[i][j]) for j in range(n)] for i in range(n)], (n, n), R)


def qqx(r, c, seed):
    R = QQ[x]
    g = random.Random(seed)
    return DomainMatrix([[R(sum(g.randint(-3, 3) * x**k for k in range(g.randint(0, 2) + 1))) for _ in range(c)]
                         for _ in range(r)], (r, c), R)


DOC = dm(ZZ, [[12, 6, 4], [3, 9, 6], [2, 16, 14]])
JORDAN = (lambda R: DomainMatrix([[R(x if i == j else 0) - R([[2, 1, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]][i][j])
                                    for j in range(4)] for i in range(4)], (4, 4), R))(QQ[x])


def make(f, m):
    return lambda: f(m)


CASES = {
    'zz_snf_doc_3x3': make(smith_normal_form, DOC),
    'zz_invf_doc_3x3': make(invariant_factors, DOC),
    'zz_decomp_doc_3x3': make(smith_normal_decomp, DOC),
    'qq_snf_3x3': make(smith_normal_form, qq(rnd(3, 3, -9, 9, 1), 1)),
    'qq_decomp_6x6': make(smith_normal_decomp, qq(rnd(6, 6, -9, 9, 2), 2)),
    'qq_invf_10x10': make(invariant_factors, qq(rnd(10, 10, -9, 9, 3), 3)),
    'gf7_invf_10x10': make(invariant_factors, dm(GF(7), rnd(10, 10, 0, 6, 4))),
    'gf101_decomp_8x8': make(smith_normal_decomp, dm(GF(101), rnd(8, 8, 0, 100, 5))),
    'gf65537_snf_20x20': make(smith_normal_form, dm(GF(65537), rnd(20, 20, 0, 65536, 6))),
    'qqx_invf_jordan_4x4': make(invariant_factors, JORDAN),
    'qqx_decomp_jordan_4x4': make(smith_normal_decomp, JORDAN),
    'qqx_snf_char_6x6': make(smith_normal_form, charmat(6, 7)),
    'qqx_decomp_3x4': make(smith_normal_decomp, qqx(3, 4, 8)),
}
f = CASES[case]
t0 = time.perf_counter(); result = f(); first = time.perf_counter() - t0
n = 1
while True:
    t0 = time.perf_counter()
    for _ in range(n): f()
    if time.perf_counter() - t0 > 0.05 or n > 100000: break
    n *= 2
reps = []
for _ in range(7):
    t0 = time.perf_counter()
    for _ in range(n): f()
    reps.append((time.perf_counter() - t0) / n)
print(json.dumps({'case': case, 'first_s': first, 'repeat_s': statistics.median(reps), 'result': repr(result)[:4000]}))
