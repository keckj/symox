"""Classify a symox run against SymPy: python compare.py raw/sympy.jsonl raw/X.jsonl

SAME      identical printed result, or the same exception class;
DECLINED  symox raises NotImplementedError where SymPy answers;
DIFF      another printed result (or another exception) with the same value,
          or a value that cannot be compared automatically;
WRONG     a different value (the difference of both values does not
          simplify to zero), or an answer where SymPy raises.
A SymPy timeout (CaseTimeout) is no reference answer: an answer of symox
there is DIFF (its correctness is checked by hand in the README).
Prints one line per case and the totals."""
import json, sys
from sympy import sympify, simplify, Tuple


def load(path):
    return {r['case']: r for r in map(json.loads, open(path))}


def same_value(a, b):
    try:
        x, y = sympify(a), sympify(b)
        if isinstance(x, Tuple) or isinstance(y, Tuple):
            return len(x) == len(y) and all(simplify(u - v) == 0 for u, v in zip(x, y))
        return simplify(x - y) == 0
    except Exception:
        return None


ref, run = load(sys.argv[1]), load(sys.argv[2])
counts = {'SAME': 0, 'DIFF': 0, 'DECLINED': 0, 'WRONG': 0}
for case, r in ref.items():
    q = run.get(case)
    if q is None:
        kind = 'DECLINED'
    elif r['error'] is not None and r['error'].startswith('CaseTimeout'):
        kind = 'SAME' if q['error'] is not None and q['error'].startswith('CaseTimeout') else 'DIFF'
    elif r['error'] is not None:
        if q['error'] is None:
            kind = 'WRONG'
        elif q['error'].split(':')[0] == r['error'].split(':')[0]:
            kind = 'SAME'
        else:
            kind = 'DIFF'
    elif q['error'] is not None:
        kind = 'DECLINED' if q['error'].startswith('NotImplementedError') else 'DIFF'
    elif q['result'] == r['result']:
        kind = 'SAME'
    else:
        ok = same_value(r['value'], q['value']) if r['value'] and q['value'] else None
        kind = 'WRONG' if ok is False else 'DIFF'
    counts[kind] += 1
    print(f'{kind:8} | {case}')
print(' '.join(f'{k}={v}' for k, v in counts.items()), f'total={len(ref)}')
