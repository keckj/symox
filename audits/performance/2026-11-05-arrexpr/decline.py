"""decline.py <worktree> <out.json>: per-case decline reason of the native a2m
port, for the 'convert' pass of the diff_a2m corpus."""
import sys, json, runpy, collections
wt, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, wt + '/python')
import symox.tensor.array.expressions._a2m as A
import symox.tensor.array.expressions.from_array_to_matrix as F
REASONS = []
_run = A.run
def run(op, exprs, extra=None, removed=False):
    try:
        return _run(op, exprs, extra, removed)
    except A._Decline as e:
        REASONS.append('Decline(enc) %s' % op); raise
    except A.Unsupported as e:
        REASONS.append('%s: %s' % (op, e)); raise
A.run = run
src = open(wt + '/audits/performance/2026-10-30-arr2mat/diff_a2m.py').read()
src = src.split("with open(out, 'w') as fh:")[0]
sys.argv = ['x', 'symox', '/dev/null']
g = {'__name__': 'diffa2m'}
exec(compile(src, 'diff_a2m', 'exec'), g)
res = {}
for name, e in g['cases']:
    if isinstance(e, str):
        continue
    del REASONS[:]
    try:
        g['ca'](e)
    except Exception as ex:
        pass
    if REASONS:
        res[name] = REASONS[0]
json.dump(res, open(out, 'w'), indent=1)
c = collections.Counter(res.values())
for k, v in c.most_common():
    print(v, k)
print('total', len(res))
