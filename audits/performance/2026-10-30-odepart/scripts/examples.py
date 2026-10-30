"""Run the SymPy 1.14 test_single.py examples of the ported hints.

usage: examples.py sympy|symox|symox-stock OUT.jsonl [TEST_SINGLE_DIR]

Examples marked dsolve_too_slow are skipped. With symox, TEST_SINGLE_DIR holds a copy of sympy/solvers/ode/tests/
test_single.py with `sympy` replaced by `symox`. For each example of the
lists below: str(dsolve(eq, func, hint=hint)) and
str(dsolve(eq, func, hint=hint, simplify=False)), or 'EXC <type>', with
the wall time of each call (first call in the process, warm caches).
"""
import json
import os
import signal
import sys
import time

# ODE_SKIP: comma separated example names (symox exceeds the memory cap on
# var_of_parameters_15, before and after this change)
SKIP = set(filter(None, os.environ.get('ODE_SKIP', '').split(',')))
LISTS = ('nth_linear_undetermined_coefficients', 'nth_linear_var_of_parameters',
         'euler_undetermined_coeff', 'euler_var_para')
if sys.argv[1] == 'sympy':
    import sympy as L
    from sympy.solvers.ode.tests import test_single as T
else:
    sys.path.insert(0, sys.argv[3])
    import symox as L
    import test_single as T
    if sys.argv[1] == 'symox-stock':
        from symox import _core
        for _k in ('ode_uc_trialset', 'ode_uc_solve', 'ode_vop_solve'):
            if hasattr(_core, _k):
                delattr(_core, _k)


def alarm(*a):
    raise TimeoutError()


signal.signal(signal.SIGALRM, alarm)


def call(fn):
    signal.alarm(60)
    try:
        return str(fn())
    except TimeoutError:
        return 'TIMEOUT'
    except Exception as e:  # noqa: BLE001
        return 'EXC %s' % type(e).__name__
    finally:
        signal.alarm(0)


with open(sys.argv[2], 'w') as out:
    for lst in LISTS:
        for d in getattr(T, '_get_examples_ode_sol_' + lst):
            if d['dsolve_too_slow'] or d['example_name'] in SKIP:
                continue
            eq, func, hint = d['eq'], d['func'], d['hint']
            t0 = time.perf_counter()
            hint_sol = call(lambda: L.dsolve(eq, func, hint=hint))
            t1 = time.perf_counter()
            nosimp = call(lambda: L.dsolve(eq, func, hint=hint, simplify=False))
            t2 = time.perf_counter()
            rec = {'name': d['example_name'], 'hint': hint, 'hint_sol': hint_sol, 'nosimp': nosimp,
                   't_hint_sol': round(t1 - t0, 4), 't_nosimp': round(t2 - t1, 4)}
            out.write(json.dumps(rec) + '\n')
            out.flush()
