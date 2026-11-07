"""compare.py REF TAG [TAG2]: classify a symox replay against the SymPy one.

SAME     same str of the evaluated side (or same exception) as SymPy
DIFF     both produce a value, different str, the assert holds in symox
         whenever it holds in SymPy (form-only / equivalent value)
DECLINED symox raises (or holds the input) where SymPy computes
WRONG    the assert holds in SymPy and fails in symox with a value
With TAG2, lists the keys that changed category.
"""
import json, os, sys

raw = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'raw')


def load(tag):
    return json.load(open(os.path.join(raw, tag + '.json')))


def classify(ref, got):
    if got is None:
        return 'DECLINED'
    if got[0] == ref[0] and (got[0] is not None or got[1] == ref[1]):
        return 'SAME'
    if isinstance(got[0], str) and got[0].startswith(('EXC', 'TIMEOUT')) and not (ref[0] or '').startswith('EXC'):
        return 'DECLINED'
    if ref[1] is True and got[1] is not True:
        return 'WRONG'
    return 'DIFF'


ref = load(sys.argv[1])
tags = sys.argv[2:]
res = {t: load(t) for t in tags}
cats = {t: {} for t in tags}
for t in tags:
    for k, v in sorted(ref.items()):
        cats[t][k] = classify(v, res[t].get(k))
for t in tags:
    counts = {}
    for c in cats[t].values():
        counts[c] = counts.get(c, 0) + 1
    print(t, len(ref), json.dumps(dict(sorted(counts.items()))))
verbose = '-v' in os.environ.get('CMP', '')
for k in sorted(ref):
    cs = [cats[t][k] for t in tags]
    if len(tags) == 2 and cs[0] != cs[1] or verbose and cs[-1] != 'SAME':
        print(' | '.join([k] + cs + [str(ref[k][0])[:70]] + [str(res[t].get(k, ['-'])[0])[:70] for t in tags]))
