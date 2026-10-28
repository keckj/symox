"""Summarize raw.jsonl: median over 3 fresh processes per tree."""
import json, statistics, sys
rows = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else 'raw.jsonl')]
cases = []
for r in rows:
    if r['data']['case'] not in cases:
        cases.append(r['data']['case'])
def med(case, tree, key):
    v = [r['data'].get(key) for r in rows if r['tree'] == tree and r['data']['case'] == case]
    v = [x for x in v if x is not None]
    return statistics.median(v) if v else None
def status(case, tree):
    d = [r['data'] for r in rows if r['tree'] == tree and r['data']['case'] == case]
    errs = {x['error'] for x in d if 'error' in x}
    fps = {str(x.get('fingerprint')) for x in d if 'error' not in x}
    return (sorted(errs)[0] if errs else 'ok'), fps
print('| case | before first (s) | after first (s) | ratio | before repeat (s) | after repeat (s) | ratio | result |')
print('|---|---|---|---|---|---|---|---|')
f = lambda v: '-' if v is None else '%.4g' % v
for c in cases:
    sb, fb = status(c, 'before'); sa, fa = status(c, 'after')
    b1, a1 = med(c, 'before', 'first'), med(c, 'after', 'first')
    br, ar = med(c, 'before', 'repeat_median'), med(c, 'after', 'repeat_median')
    r1 = '-' if sb != 'ok' or not a1 else '%.2fx' % (b1 / a1)
    rr = '-' if not (br and ar) else '%.2fx' % (br / ar)
    if sb != 'ok':
        res = 'before refused (%s); after ok' % sb[:60] if sa == 'ok' else 'both refused'
    else:
        res = 'identical' if fb == fa else 'DIFFERENT %s %s' % (fb, fa)
    print('| %s | %s | %s | %s | %s | %s | %s | %s |' % (c, f(b1), f(a1), r1, f(br), f(ar), rr, res))
