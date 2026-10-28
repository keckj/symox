"""Before/after harness: python bench.py CASE -> JSON line (first call, median of repeats).

Run once per fresh process with PYTHONPATH pointing at the tree under test.
"""
import hashlib, json, statistics, sys, time

CASES = {
    'res_witness_200_150': "resultant(expand((x + 1)**200), expand((x - 2)**150 + 1), x)",
    'res_60_45': "resultant(expand((x + 1)**60), expand((x - 2)**45 + 1), x)",
    'res_dense_40_30': "resultant(sum((3*k**3 - 7*k + 11)*x**k for k in range(41)), sum((5*k**2 + 2*k - 13)*x**k for k in range(31)), x)",
    'res_small_zz': "resultant(3*x**4 - 2*x + 5, 7*x**3 + x**2 - 4, x)",
    'res_small_qq': "resultant(x**3/2 + x - Rational(1, 3), Rational(3, 5)*x**2 - 7, x)",
    'res_prs_small': "resultant(x**3 - 2*x + 1, 2*x**2 + 3, x, includePRS=True)",
    'res_bivariate': "resultant(x**2 + y, x*y - 1, x)",
    'poly_res_25_20': "Poly(expand((x - 3)**25 + x), x).resultant(Poly(expand((2*x + 1)**20 - 5), x))",
    'disc_dense_25': "discriminant(sum((k**2 - 3*k + 1)*x**k for k in range(26)), x)",
    'disc_small': "discriminant(5*x**5 - 3*x**2 + 7)",
    'cyclotomic_20001': "cyclotomic_poly(20001, x)",
    'cyclotomic_15015': "cyclotomic_poly(15015, x)",
    'cyclotomic_19999': "cyclotomic_poly(19999, x)",
    'factor_x1155_m1': "factor_list(x**1155 - 1)",
    'decompogen_4000': "decompogen(x**4000 + x**2, x)",
    'decompogen_4100': "decompogen(x**4100 + x**2, x)",
}


def main():
    name = sys.argv[1]
    repeats = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    t0 = time.perf_counter()
    from symox import (Poly, Rational, cyclotomic_poly, decompogen, discriminant,
                       expand, factor_list, resultant, symbols)
    import_time = time.perf_counter() - t0
    x, y = symbols('x y')
    env = dict(locals())
    code = compile(CASES[name], name, 'eval')
    out = {'case': name, 'import': import_time}
    try:
        t = time.perf_counter(); value = eval(code, env); out['first'] = time.perf_counter() - t
        times = []
        for _ in range(repeats):
            t = time.perf_counter(); eval(code, env); times.append(time.perf_counter() - t)
        out['repeat_median'] = statistics.median(times)
        # Deterministic across processes (str hashes are randomized).
        out['fingerprint'] = hashlib.sha256(repr(value).encode()).hexdigest()[:16]
    except Exception as error:  # a refusal is reported, never timed as a result
        out['error'] = '%s: %s' % (type(error).__name__, str(error)[:120])
        out['first'] = time.perf_counter() - t
    print(json.dumps(out))


if __name__ == '__main__':
    import sys as _s
    _s.set_int_max_str_digits(0)
    main()
