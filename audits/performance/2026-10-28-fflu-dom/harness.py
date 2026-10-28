"""Before/after timing of the fraction-free DomainMatrix entry points.

Usage: python harness.py CASE first|repeat  (PYTHONPATH selects the tree)
       python harness.py --drive BEFORE_PYTHONPATH AFTER_PYTHONPATH OUT.json
"""
import json, os, statistics, subprocess, sys, time

CASES = ['zzi_fflu', 'zzi_det', 'zzi_charpoly', 'zzi_rref_den',
         'gf_rref_den', 'gf_det', 'gf_fflu', 'qqa_det', 'qqa_rref_den',
         'zzx_det', 'zzx_rref_den', 'zzx_charpoly', 'zzx_fflu',
         'matrix_bareiss_zzi', 'matrix_bareiss_zzx']


def build(case):
    from symox import I, Matrix, sqrt, symbols
    from symox.polys.domains import GF, QQ, ZZ, ZZ_I
    from symox.polys.matrices import DomainMatrix
    x = symbols('x')
    def entries(n, f):
        return [[f(i, j, (7*i*i + 3*j*j + 5*i*j + i + 2*j) % 11 - 5) for j in range(n)]
                for i in range(n)]
    def dm(K, rows):
        return DomainMatrix([[K.from_sympy(e) for e in r] for r in rows],
                            (len(rows), len(rows[0])), K)
    kind = case.split('_')[0]
    if kind == 'zzi':
        M = dm(ZZ_I, entries(7, lambda i, j, v: v + ((i + 2*j) % 5 - 2)*I))
    elif kind == 'gf':
        M = dm(GF(10007), entries(14, lambda i, j, v: v*(i + 1) + j*j))
    elif kind == 'qqa':
        K = QQ.algebraic_field(sqrt(2))
        M = dm(K, entries(7, lambda i, j, v: v + ((i*j) % 3 - 1)*sqrt(2)))
    elif kind == 'zzx':
        M = dm(ZZ[x], entries(5, lambda i, j, v: v + ((i + j) % 3 - 1)*x + int(i == j)*x**2))
    elif case == 'matrix_bareiss_zzi':
        return lambda: Matrix(entries(6, lambda i, j, v: v + ((i + 2*j) % 5 - 2)*I)).det(method='bareiss')
    elif case == 'matrix_bareiss_zzx':
        return lambda: Matrix(entries(5, lambda i, j, v: v + ((i + j) % 3 - 1)*x + int(i == j)*x**2)).det(method='bareiss')
    op = case.split('_', 1)[1]
    if op == 'fflu':
        return lambda: M.fflu()
    if op == 'det':
        return lambda: M.det()
    if op == 'charpoly':
        return lambda: M.charpoly()
    if op == 'rref_den':
        return lambda: M.rref_den()
    raise ValueError(case)


def run(case, mode):
    call = build(case)
    if mode == 'first':
        t = time.perf_counter(); call(); return [time.perf_counter() - t]
    call()
    out = []
    for _ in range(7):
        t = time.perf_counter(); call(); out.append(time.perf_counter() - t)
    return out


def drive(before, after, path):
    results = {}
    for case in CASES:
        for label, tree in (('before', before), ('after', after)):
            env = dict(os.environ, PYTHONPATH=tree)
            first = []
            for _ in range(5):
                p = subprocess.run([sys.executable, __file__, case, 'first'], env=env,
                                   capture_output=True, text=True, check=True)
                first += json.loads(p.stdout)
            p = subprocess.run([sys.executable, __file__, case, 'repeat'], env=env,
                               capture_output=True, text=True, check=True)
            results.setdefault(case, {})[label] = {
                'first_runs': first, 'first_median': statistics.median(first),
                'repeat_runs': json.loads(p.stdout),
                'repeat_median': statistics.median(json.loads(p.stdout))}
        r = results[case]
        print('%-20s first %9.3f ms -> %9.3f ms (x%.2f)  repeat %9.3f ms -> %9.3f ms (x%.2f)' % (
            case, 1e3*r['before']['first_median'], 1e3*r['after']['first_median'],
            r['before']['first_median']/r['after']['first_median'],
            1e3*r['before']['repeat_median'], 1e3*r['after']['repeat_median'],
            r['before']['repeat_median']/r['after']['repeat_median']), flush=True)
    with open(path, 'w') as f:
        json.dump(results, f, indent=1)


if __name__ == '__main__':
    if sys.argv[1] == '--drive':
        drive(*sys.argv[2:5])
    else:
        print(json.dumps(run(sys.argv[1], sys.argv[2])))
