"""pytest plugin / runner: counts native RootSum calls and declines (None)."""
import atexit, collections, os
from symox import _core, _public_roots
counts = collections.Counter()
def wrap(name):
    f = getattr(_core, name)
    def g(*a, **k):
        r = f(*a, **k)
        counts[(name, 'decline' if r is None else 'native')] += 1
        return r
    setattr(_core, name, g)
for n in ('_root_sum_transform', '_root_sum_build', '_root_sum_new', '_root_sum_doit', '_root_sum_evalf'):
    wrap(n)
@atexit.register
def report():
    with open(os.environ.get('RSCOUNT_OUT', 'rscount.txt'), 'a') as fh:
        fh.write(os.environ.get('RSCOUNT_LABEL', '') + ' ' + repr(sorted(counts.items())) + '\n')
