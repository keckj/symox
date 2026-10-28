"""One fresh-process measurement: python bench_case.py <case> -> JSON line."""
import json, random, statistics, sys, time

case = sys.argv[1]
import symox
from symox import ZZ, Matrix, Poly, symbols, cyclotomic_poly
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.normalforms import smith_normal_form, smith_normal_decomp, invariant_factors
from symox.polys.numberfields.basis import round_two
from symox.polys.numberfields.primes import prime_decomp
from symox.polys.numberfields import _normalforms as nfh

x = symbols('x')


def dm(rows):
    return DomainMatrix([[ZZ(v) for v in r] for r in rows], (len(rows), len(rows[0])), ZZ)


def rnd(r, c, lo, hi, seed):
    g = random.Random(seed)
    return [[g.randint(lo, hi) for _ in range(c)] for _ in range(r)]


DOC = dm([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
CASES = {
    'snf_doc_3x3': lambda: smith_normal_form(DOC),
    'invf_doc_3x3': lambda: invariant_factors(DOC),
    'decomp_doc_3x3': lambda: smith_normal_decomp(DOC),
    'decomp_6x6': (lambda m: lambda: smith_normal_decomp(m))(dm(rnd(6, 6, -100, 100, 1))),
    'invf_8x8': (lambda m: lambda: invariant_factors(m))(dm(rnd(8, 8, -100, 100, 2))),
    'decomp_4x6_big': (lambda m: lambda: smith_normal_decomp(m))(dm(rnd(4, 6, -10**20, 10**20, 3))),
    'invf_12x12': (lambda m: lambda: invariant_factors(m))(dm(rnd(12, 12, -9, 9, 4))),
    'matrix_snf_4x4': (lambda M: lambda: __import__("symox.matrices.normalforms", fromlist=["x"]).smith_normal_form(M, domain=ZZ))(
        Matrix(rnd(4, 4, -20, 20, 5))),
    'round_two_x3m2': lambda: round_two(Poly(x**3 - 2)),
    'round_two_cyclo7': lambda: round_two(Poly(cyclotomic_poly(7, x))),
    'prime_decomp_5': lambda: prime_decomp(5, T=Poly(x**3 + x**2 - 2*x + 8)),
    'nf_hnf_5x8': (lambda m: lambda: nfh.hermite_normal_form(m))(dm(rnd(5, 8, -50, 50, 6))),
    'nf_hnf_mod_4x4': (lambda m: lambda: nfh.hermite_normal_form(m, D=int(abs(m.det())) or 1))(
        dm(rnd(4, 4, -30, 30, 7))),
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
print(json.dumps({'case': case, 'first_s': first, 'repeat_s': statistics.median(reps), 'result': repr(result)[:2000]}))
