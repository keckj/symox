"""Per-example pass/fail of SymPy 1.14 test_single.py example lists for the
nth_algebraic, factorable, 2nd_linear_bessel and lie_group hints.

usage: per_example.py <test_single.py copy dir> <out.json>
Each example runs in a forked child with a 300 s limit."""
import json
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, sys.argv[1])
import test_single as T  # noqa: E402

LISTS = {
    'nth_algebraic': T._get_examples_ode_sol_nth_algebraic,
    'factorable': T._get_examples_ode_sol_factorable,
    '2nd_linear_bessel': T._get_examples_ode_sol_2nd_linear_bessel,
    'lie_group': T._get_examples_ode_sol_lie_group,
}


def run_one(q, hint, name):
    for ex in LISTS[hint]:
        if ex['example_name'] == name:
            try:
                T._test_particular_example(ex['hint'], ex, solver_flag=True)
                q.put(('pass', ''))
            except BaseException as e:  # noqa: BLE001
                q.put(('fail', '%s: %s' % (type(e).__name__, str(e)[:300])))
            return


def main():
    out = {}
    ctx = mp.get_context('fork')
    for hint, fn in LISTS.items():
        for ex in fn:
            name = ex['example_name']
            q = ctx.Queue()
            p = ctx.Process(target=run_one, args=(q, hint, name))
            t0 = time.time()
            p.start()
            p.join(300)
            if p.is_alive():
                p.kill()
                p.join()
                res = ('timeout', '')
            else:
                try:
                    res = q.get(timeout=5)
                except Exception:
                    res = ('crash', str(p.exitcode))
            out['%s/%s' % (hint, name)] = {'status': res[0], 'detail': res[1], 'slow': ex['slow'],
                                            'secs': round(time.time() - t0, 1)}
            print(hint, name, res[0], round(time.time() - t0, 1), res[1][:120], flush=True)
    json.dump(out, open(sys.argv[2], 'w'), indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
