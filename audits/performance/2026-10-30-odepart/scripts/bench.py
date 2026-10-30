"""usage: bench.py sympy|symox|symox-stock OUT.json [runs]

Timings of dsolve(eq, f(x), hint=h, simplify=s) on the cases below. Each
(case, simplify) is measured in `runs` fresh processes: first call (module
imported beforehand) and median of 5 repeated calls; the reported values
are medians over the runs. Also times the helpers directly
(`_undetermined_coefficients_match`, `_solve_undetermined_coefficients`)
on the undetermined-coefficients cases.

`symox-stock` removes the three new kernels from `_core` (reference bodies).
"""
import json, subprocess, sys, statistics

CASES = [
    ('uc_2nd', "f(x).diff(x,2) + 3*f(x).diff(x) + 2*f(x) - x**2*exp(-x) - sin(2*x)", 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('uc_3rd', "f(x).diff(x,3) - f(x).diff(x) - x*exp(x)*cos(x) - 3", 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('uc_res4', "f(x).diff(x,4) + 2*f(x).diff(x,2) + f(x) - x**2*sin(x)", 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('uc_poly', "f(x).diff(x,2) - f(x) - x**3 - 2*x + 1", 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('vop_sec', "f(x).diff(x,2) + f(x) - 1/cos(x)", 'nth_linear_constant_coeff_variation_of_parameters'),
    ('vop_expx', "f(x).diff(x,2) - 2*f(x).diff(x) + f(x) - exp(x)/x", 'nth_linear_constant_coeff_variation_of_parameters'),
    ('vop_3rd', "f(x).diff(x,3) - f(x).diff(x) - exp(2*x)", 'nth_linear_constant_coeff_variation_of_parameters'),
    ('euler_uc', "x**2*f(x).diff(x,2) - 2*x*f(x).diff(x) + 2*f(x) - log(x)", 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients'),
    ('euler_vop', "x**2*f(x).diff(x,2) - 2*x*f(x).diff(x) + 2*f(x) - x**4", 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters'),
]

CHILD = r'''
import json, sys, time, statistics
mod, src, hint, simp = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == 'True'
M = __import__('sympy' if mod == 'sympy' else 'symox')
if mod == 'symox-stock':
    from symox import _core
    for k in ('ode_uc_trialset', 'ode_uc_solve', 'ode_vop_solve'):
        if hasattr(_core, k):
            delattr(_core, k)
ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
x = M.Symbol('x'); f = M.Function('f'); ns.update(x=x, f=f)
eq = eval(src, ns)
call = lambda: M.dsolve(eq, f(x), hint=hint, simplify=simp)
t = time.perf_counter(); r = call(); first = time.perf_counter() - t
rep = []
for _ in range(5):
    t = time.perf_counter(); call(); rep.append(time.perf_counter() - t)
out = {'first': first, 'repeat': statistics.median(rep), 'result': str(r)}
if 'undetermined' in hint and 'euler' not in hint:
    nh = __import__(M.__name__ + '.solvers.ode.nonhomogeneous', fromlist=['x'])
    single = __import__(M.__name__ + '.solvers.ode.single', fromlist=['x'])
    prob = single.SingleODEProblem(eq, f(x), x)
    solver = single.NthLinearConstantCoeffUndeterminedCoefficients(prob)
    assert solver.matches()
    r = dict(solver.r); forcing = r[-1]
    hom = M.Add(prob.eq_high_order_free, -forcing)
    rep = []
    for _ in range(5):
        t = time.perf_counter(); ts = nh._undetermined_coefficients_match(forcing, x, f(x), hom)['trialset']; rep.append(time.perf_counter() - t)
    out['trialset'] = statistics.median(rep)
    roots, _ = nh._get_const_characteristic_eq_sols(r, f(x), prob.order)
    consts = prob.get_numbered_constants(num=len(roots))
    r.update({'list': roots, 'sol': M.Eq(f(x), M.Add(*[c*g for c, g in zip(consts, roots)]))})
    rep = []
    for _ in range(5):
        t = time.perf_counter(); nh._solve_undetermined_coefficients(prob.eq, f(x), prob.order, r, ts); rep.append(time.perf_counter() - t)
    out['solve'] = statistics.median(rep)
print(json.dumps(out))
'''


def main():
    mod, out = sys.argv[1], sys.argv[2]
    runs = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    res = {}
    for name, src, hint in CASES:
        for simp in (True, False):
            samples = []
            for _ in range(runs):
                p = subprocess.run([sys.executable, '-c', CHILD, mod, src, hint, str(simp)], capture_output=True, text=True, timeout=600)
                try:
                    samples.append(json.loads(p.stdout.strip().splitlines()[-1]))
                except Exception:
                    samples.append({'error': p.stderr[-500:]})
            ok = [s for s in samples if 'first' in s]
            entry = {'result': ok[0]['result'] if ok else None}
            for k in ('first', 'repeat', 'trialset', 'solve'):
                vals = [s[k] for s in ok if k in s]
                entry[k] = statistics.median(vals) if vals else None
            if not ok:
                entry['error'] = samples[0].get('error')
            key = '%s/%s' % (name, simp)
            res[key] = entry
            print(key, {k: (round(v * 1e3, 2) if isinstance(v, float) else v) for k, v in entry.items() if k != 'result'}, flush=True)
    json.dump(res, open(out, 'w'), indent=1)


main()
