"""Compare a symox replay of the radical corpus with the SymPy 1.14 replay.

usage: compare.py SYMPY.json SYMOX.json [--verbose]

Categories:
  SAME       identical printed result (or the same exception class);
  EQUIV      same solution multiset, values equal at random points;
  SUBSET     symox misses solutions SymPy returns (an empty symox list
             against a non-empty SymPy list is counted as WRONG_EMPTY);
  SUPERSET   symox returns more solutions than SymPy;
  DIFF       other differences;
  DECLINED   symox raises where SymPy returns a result;
  EXTRA      symox returns a result where SymPy raises;
  ERRDIFF    both raise, different exception classes;
  TIMEOUT    either side timed out.
"""
import json
import random
import sys

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import sympy  # noqa: E402
from cases import POSITIVE, SYMBOLS  # noqa: E402

NS = {k: getattr(sympy, k) for k in dir(sympy) if not k.startswith('_')}
for _n in SYMBOLS.split():
    NS[_n] = sympy.Symbol(_n)
for _n in POSITIVE.split():
    NS[_n] = sympy.Symbol(_n, positive=True)


def values(text):
    obj = eval(text, dict(NS))
    if isinstance(obj, tuple) and len(obj) == 2 and isinstance(obj[1], set):
        obj = [t[0] if len(t) == 1 else t for t in obj[1]]
    if isinstance(obj, dict):
        obj = [obj]
    out = []
    for v in obj:
        if isinstance(v, dict):
            v = tuple(v[k] for k in sorted(v, key=str))
            v = v[0] if len(v) == 1 else sympy.Tuple(*v)
        out.append(sympy.sympify(v))
    return out


def close(a, b, points):
    for point in points:
        try:
            va = complex(sympy.N(a.subs(point), 30))
            vb = complex(sympy.N(b.subs(point), 30))
        except (TypeError, ValueError, ZeroDivisionError):
            return sympy.simplify(a - b) == 0
        if abs(va - vb) > 1e-9 * max(1.0, abs(va), abs(vb)):
            return False
    return True


def matching(sa, sb):
    free = set()
    for v in sa + sb:
        free |= v.free_symbols
    rng = random.Random(1)
    points = [{s: sympy.Rational(rng.randint(3, 40), rng.randint(2, 9)) for s in free} for _ in range(3)]
    left = list(sb)
    unmatched = []
    for a in sa:
        for i, b in enumerate(left):
            if a == b or close(a, b, points):
                left.pop(i)
                break
        else:
            unmatched.append(a)
    return unmatched, left


def classify(ref, got):
    if ref.get('error') == 'Timeout' or got.get('error') == 'Timeout':
        return 'TIMEOUT'
    if not ref['ok'] and not got['ok']:
        return 'SAME' if ref['error'] == got['error'] else 'ERRDIFF'
    if not got['ok']:
        return 'DECLINED'
    if not ref['ok']:
        return 'EXTRA'
    if ref['str'] == got['str']:
        return 'SAME'
    try:
        sr, sg = values(ref['str']), values(got['str'])
    except Exception:
        return 'DIFF'
    if sr and not sg:
        return 'WRONG_EMPTY'
    extra, missing = matching(sg, sr)
    if not extra and not missing:
        return 'EQUIV'
    if missing and not extra:
        return 'SUBSET'
    if extra and not missing:
        return 'SUPERSET'
    return 'DIFF'


def main():
    ref = json.load(open(sys.argv[1]))
    got = json.load(open(sys.argv[2]))
    verbose = '--verbose' in sys.argv
    counts = {}
    for cid in ref:
        c = classify(ref[cid], got[cid])
        counts[c] = counts.get(c, 0) + 1
        if verbose or c not in ('SAME', 'EQUIV'):
            print('%-12s %-20s sympy=%s\n%34s symox=%s' % (
                c, cid, ref[cid]['str'][:200] if ref[cid]['ok'] else ref[cid]['error'], '',
                got[cid]['str'][:200] if got[cid]['ok'] else got[cid]['error'] + ': ' + got[cid]['str'][:100]))
    print(json.dumps(dict(sorted(counts.items())), sort_keys=True))


if __name__ == '__main__':
    main()
