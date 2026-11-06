"""Differential replay of the matcher corpus.

usage:
  compare.py sympy CORPUS.json.gz OUT.json      # SymPy 1.14 Expr.match
  compare.py symox CORPUS.json.gz OUT.json      # symox implementations
  compare.py report SYMPY.json SYMOX.json CORPUS.json.gz OUT.json

symox implementations (one process, sympy imports blocked):
  python_port  symox._matching.matches: the Python port of the reference
               protocol, i.e. the public Expr.match before this lot;
  pattern      _core.Expr.match: the native backtracking pattern_match;
  native       _core._reference_match: the native port (symox::matches),
               'declined' when it falls back;
  public       Expr.match after this lot (native, else the Python port).

A result is the srepr of each binding, compared as a dict (and its
insertion order, reported separately).  Errors are compared by class name.
"""
import gzip, json, signal, sys, time


class _Alarm(Exception):
    pass


def _alarm(*_):
    raise _Alarm()


def _render(srepr, value):
    if value is None:
        return None
    return [[srepr(k), srepr(v)] for k, v in value.items()]


def _run(fn, srepr, limit=20):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(limit)
    start = time.perf_counter()
    try:
        value = fn()
        out = {'ok': True, 'value': _render(srepr, value)}
    except _Alarm:
        out = {'ok': False, 'error': 'DEADLINE'}
    except RecursionError:
        out = {'ok': False, 'error': 'RecursionError'}
    except Exception as error:
        out = {'ok': False, 'error': type(error).__name__}
    finally:
        signal.alarm(0)
    out['seconds'] = time.perf_counter() - start
    return out


def replay_sympy(cases):
    import sympy
    ns = vars(sympy)
    rows = []
    for case in cases:
        expr, pattern = eval(case['expr'], ns), eval(case['pattern'], ns)
        rows.append({'id': case['id'], 'sympy': _run(lambda: expr.match(pattern, old=case['old']), sympy.srepr)})
    return rows


def replay_symox(cases):
    class NoReference:
        def find_spec(self, fullname, path=None, target=None):
            if fullname == 'sympy' or fullname.startswith('sympy.'):
                raise ImportError('reference forbidden')
    sys.meta_path.insert(0, NoReference())
    import symox
    from symox import _core, _matching, _expr_methods
    ns = vars(symox)
    rows = []

    def native(expr, pattern, old):
        handled, value = _core._reference_match(expr, pattern, old)
        if not handled:
            raise _Declined(value)
        return value

    class _Declined(Exception):
        pass
    for case in cases:
        expr, pattern = eval(case['expr'], ns), eval(case['pattern'], ns)
        old = case['old']
        row = {'id': case['id']}
        row['python_port'] = _run(lambda: _matching.matches(pattern, expr, old=old), symox.srepr)
        row['pattern'] = _run(lambda: _expr_methods._native_match(expr, pattern), symox.srepr)
        if hasattr(_core, '_reference_match'):
            r = _run(lambda: native(expr, pattern, old), symox.srepr)
            if not r['ok'] and r['error'] == '_Declined':
                r['error'] = 'declined'
            row['native'] = r
        row['public'] = _run(lambda: expr.match(pattern, old=old), symox.srepr)
        rows.append(row)
    return rows


def _key(r):
    # symox prints its WildFunction (a Wild subclass) as Wild('g'): the
    # srepr differs from SymPy's WildFunction('g') for the same binding.
    if not r['ok']:
        return ('error', r['error'])
    if r['value'] is None:
        return ('value', None)
    return ('value', tuple(sorted((k.replace('WildFunction(', 'Wild('), v) for k, v in r['value'])))


def report(sympy_rows, symox_rows, cases):
    from collections import Counter, defaultdict
    by_id = {c['id']: c for c in cases}
    ref = {r['id']: r['sympy'] for r in sympy_rows}
    impls = [k for k in ('python_port', 'pattern', 'native', 'public') if k in symox_rows[0]]
    counts = {k: Counter() for k in impls}
    per_source = {k: defaultdict(Counter) for k in impls}
    order = Counter()
    diffs = defaultdict(list)
    timing = defaultdict(float)
    for row in symox_rows:
        r = ref[row['id']]
        src = by_id[row['id']]['source']
        timing['sympy'] += r['seconds']
        for k in impls:
            mine = row[k]
            timing[k] += mine['seconds']
            if k == 'native' and not mine['ok'] and mine['error'] == 'declined':
                verdict = 'declined'
            elif _key(mine) == _key(r):
                verdict = 'equal'
                if k in ('native', 'public') and r['ok'] and r['value'] and [k2 for k2, _ in mine['value']] != [k2.replace('WildFunction(', 'Wild(') for k2, _ in r['value']]:
                    order[k] += 1
            else:
                verdict = 'different'
                diffs[k].append(row['id'])
            counts[k][verdict] += 1
            per_source[k][src][verdict] += 1
    return {'implementations': impls, 'counts': {k: dict(v) for k, v in counts.items()},
            'per_source': {k: {s: dict(c) for s, c in v.items()} for k, v in per_source.items()},
            'binding_order_differs_when_equal': dict(order),
            'different_ids': dict(diffs), 'seconds': dict(timing)}


def main():
    mode = sys.argv[1]
    if mode in ('sympy', 'symox'):
        cases = json.loads(gzip.open(sys.argv[2]).read())['cases']
        rows = replay_sympy(cases) if mode == 'sympy' else replay_symox(cases)
        with open(sys.argv[3], 'w') as fh:
            json.dump(rows, fh)
        return
    sympy_rows = json.load(open(sys.argv[2]))
    symox_rows = json.load(open(sys.argv[3]))
    cases = json.loads(gzip.open(sys.argv[4]).read())['cases']
    out = report(sympy_rows, symox_rows, cases)
    with open(sys.argv[5], 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in ('counts', 'binding_order_differs_when_equal', 'seconds')}, indent=1))


if __name__ == '__main__':
    main()
