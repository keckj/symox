"""Share of nth_algebraic / factorable / 2nd_linear_bessel decisions taken
by the native kernels on every example of SymPy 1.14 test_single.py.

usage: native_share.py <dir holding the symox copy of test_single.py> <out.json>

Each wrapper of python/symox/solvers/ode/_classify_native.py is counted as
"native" when it returns a decision and "declined" when it returns None
(the reference body then runs). classify_ode runs on every example, then
dsolve with the example hint when that hint is one of the three."""
import json
import multiprocessing as mp
import sys
from collections import Counter

sys.path.insert(0, sys.argv[1])
import test_single as T  # noqa: E402
from symox import classify_ode, dsolve  # noqa: E402
from symox.solvers.ode import _classify_native as cn  # noqa: E402

WRAPPERS = ('nth_algebraic', 'factorable_match', 'factorable_solutions', 'bessel_match', 'bessel_solution')
HINTS = ('nth_algebraic', 'factorable', '2nd_linear_bessel')


def child(q, idx):
    counts = Counter()
    declined = []
    ex = T._get_all_examples()[idx]
    for name in WRAPPERS:
        orig = getattr(cn, name)

        def wrapped(*a, _orig=orig, _name=name, **k):
            v = _orig(*a, **k)
            counts[(_name, 'declined' if v is None else 'native')] += 1
            if v is None:
                declined.append(_name)
            return v
        setattr(cn, name, wrapped)
    try:
        classify_ode(ex['eq'], ex['func'])
        if ex['hint'] in HINTS and not ex['dsolve_too_slow']:
            dsolve(ex['eq'], ex['func'], hint=ex['hint'], simplify=ex['simplify_flag'])
    except Exception:  # noqa: BLE001
        pass
    q.put(({'%s:%s' % k: v for k, v in counts.items()}, declined))


def main():
    ctx = mp.get_context('fork')
    total = Counter()
    declines = {}
    examples = T._get_all_examples()
    for i, ex in enumerate(examples):
        q = ctx.Queue()
        p = ctx.Process(target=child, args=(q, i))
        p.start()
        p.join(120)
        if p.is_alive():
            p.kill()
            p.join()
            total['timeout'] += 1
            continue
        try:
            counts, declined = q.get(timeout=5)
        except Exception:  # noqa: BLE001
            total['crash'] += 1
            continue
        total.update(counts)
        if declined:
            declines[ex['example_name']] = sorted(set(declined))
    out = {'examples': len(examples), 'counts': dict(sorted(total.items())), 'declines': declines}
    json.dump(out, open(sys.argv[2], 'w'), indent=1, sort_keys=True)
    print(json.dumps(out['counts']))
    print('declining examples:', len(declines))


if __name__ == '__main__':
    main()
