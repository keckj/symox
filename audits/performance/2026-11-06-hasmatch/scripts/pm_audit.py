"""Classify the logged calls of the old ``pattern_match`` by internal caller.

The log (env ``SYMOX_PM_AUDIT`` of an instrumented build) holds one line per
call: ``caller \t expr \t pattern \t old-answer \t reference-answer``, where
the reference answer comes from ``matches::reference_match`` (the native
port of SymPy's ``matches``) or is ``DECLINED``.  Declined calls are
replayed in SymPy 1.14 (``sympify`` of the printed forms, Wilds ``p q r a b
c d``), which is the arbiter for both.

Usage: pm_audit.py LOG OUT.json
"""
import collections
import json
import signal
import sys

import sympy

# Wilds print with a trailing underscore; names of functions that also
# occur as plain symbols in the logged conditions are read as symbols.
WILDS = {name + '_': sympy.Wild(name) for name in 'p q r a b c d n'.split()}
WILDS.update({name: sympy.Symbol(name) for name in ('beta', 'gamma', 'zeta', 'lambda_', 'E1', 'S1')})


def normal(answer):
    if answer in ('None', 'ERROR', 'DECLINED'):
        return answer
    return '; '.join(sorted(part for part in answer.split('; ') if part))


class _Timeout(Exception):
    pass


def _alarm(*_):
    raise _Timeout()


def sympy_answer(expr, pattern, limit=10):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(limit)
    try:
        ns = dict(WILDS)
        e = sympy.sympify(expr, locals=ns)
        p = sympy.sympify(pattern, locals=ns)
        m = e.match(p)
    except _Timeout:
        return 'TIMEOUT'
    except Exception as error:
        return 'UNPARSED %s' % type(error).__name__
    finally:
        signal.alarm(0)
    if m is None:
        return 'None'
    return '; '.join(sorted('%s=%s' % (k, v) for k, v in m.items()))


def main(log, out):
    seen = set()
    rows = []
    for line in open(log):
        parts = line.rstrip('\n').split('\t')
        if len(parts) != 5:
            continue
        key = tuple(parts)
        if key in seen:
            continue
        seen.add(key)
        rows.append(parts)
    counts = collections.Counter()
    details = collections.defaultdict(list)
    for caller, expr, pattern, mine, ref in rows:
        mine, ref = normal(mine), normal(ref)
        if ref == 'DECLINED':
            ref = sympy_answer(expr, pattern)
            source = 'sympy'
        else:
            source = 'reference'
        if ref.startswith('UNPARSED') or ref == 'TIMEOUT':
            verdict = ref.split()[0]
        elif mine == ref:
            verdict = 'SAME'
        elif mine == 'None':
            verdict = 'OLD-NONE'
        elif ref == 'None':
            verdict = 'OLD-SPURIOUS'
        else:
            verdict = 'OLD-OTHER-BINDING'
        counts[caller, source, verdict] += 1
        if verdict != 'SAME':
            details[caller].append(dict(expr=expr, pattern=pattern, old=mine, reference=ref,
                                        source=source, verdict=verdict))
    for key in sorted(counts):
        print('%-12s %-10s %-18s %d' % (key + (counts[key],)))
    json.dump(dict(counts={'|'.join(k): v for k, v in counts.items()}, details=details),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
