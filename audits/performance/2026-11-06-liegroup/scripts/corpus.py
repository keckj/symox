"""Differential corpus of the native ``lie_group`` solver (wave 44).

usage:
  corpus.py sympy OUT.json            # SymPy 1.14
  corpus.py symox OUT.json            # symox on PYTHONPATH (sympy blocked)
  corpus.py report SYMPY.json SYMOX.json OUT.json   # SymPy validates symox

Each case is a Python expression evaluated in the package namespace:
``infinitesimals`` for every hint (default, all, the nine heuristics),
``checkinfsol`` of the ``all`` infinitesimals, ``'lie_group' in
classify_ode`` and ``dsolve(..., hint='lie_group')``. The equations are
those of SymPy's ``test_lie_group.py`` and the 20 ``lie_group_*`` rows of
``test_single.py`` (also the lie_group rows of
``audits/performance/2026-11-05-ode``).

A result is its ``srepr`` (dicts as ordered pair lists) or the error class
and message. The report classes each case as
  SAME      equal value (dict order ignored), or equal error class;
  DECLINED  symox raises (or times out) where SymPy answers;
  WRONG     symox answers and SymPy rejects the answer (``checkinfsol`` of
            the infinitesimals, ``checkodesol`` of the solutions);
  DIFF      any other difference (``validated_by_sympy`` records whether
            SymPy's check accepted symox's differing answer).
``Dummy`` indices are numbered by first appearance before comparing.
For symox the run also counts the native calls (``_core.ode_lie_*``) that
declined to the Python reference body.
"""
import json, signal, sys, time

HINTS = ['default', 'all', 'abaco1_simple', 'abaco1_product', 'abaco2_similar',
         'abaco2_unique_unknown', 'abaco2_unique_general', 'linear',
         'function_sum', 'bivariate', 'chi']

EQS = {
    # test_lie_group.py
    'h1_eq': 'Eq(f(x).diff(x), x**2*f(x))',
    'h1_eq1': 'f(x).diff(x) + a*f(x) - c*exp(b*x)',
    'h1_eq2': 'f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)',
    'h1_eq3': '(1 + 2*x)*f(x).diff(x) + 2 - 4*exp(-f(x))',
    'h1_eq4': 'f(x).diff(x) - (a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0)**Rational(-1, 2)',
    'h1_eq5': 'x**2*f(x).diff(x) - f(x) + x**2*exp(x - (1/x))',
    'h1_eq6': 'f(x).diff(x) - (f(x)/x)*(x*log(x**2/f(x)) + 2)',
    'h1_eq7': 'x*(f(x).diff(x)) + 1 - f(x)**2',
    'h3_eq1': 'x**2*f(x).diff(x) + x*f(x) + f(x)**2 + x**2',
    'h3_eq2': 'x**2*(-f(x)**2 + f(x).diff(x))- a*x**2*f(x) + 2 - a*x',
    'fsum': 'f(x).diff(x) - (3*(1 + x**2/f(x)**2)*atan(f(x)/x) + (1 - 2*f(x))/x + (1 - 3*f(x))*(x/f(x)**2))',
    'sim1': 'f(x).diff(x) - F(a*x + b*f(x))',
    'sim2': 'f(x).diff(x) - (f(x)**2 / (sin(f(x) - x) - x**2 + 2*x*f(x)))',
    'uu1': 'f(x).diff(x) - x**(a - 1)*(f(x)**(1 - b))*F(x**a/a + f(x)**b/b)',
    'uu2': 'f(x).diff(x) + tan(F(x**2 + f(x)**2) + atan(x/f(x)))',
    'uu3': '(x*f(x).diff(x) + f(x) + 2*x)**2 -4*x*f(x) -4*x**2 -4*a',
    'lin1': 'x**(n*(m + 1) - m)*(f(x).diff(x)) - a*f(x)**n -b*x**(n*(m + 1))',
    'kamke': 'x**2*(a*f(x)**2+(f(x).diff(x))) + b*x**alpha + c',
    'riccati': 'f(x).diff(x) - f(x)**2',
    # test_single.py lie_group_01..20
    'lg01': 'x*f(x).diff(x)*(f(x)+4) + (f(x)**2) -2*f(x)-2*x',
    'lg03': 'Eq(x**7*Derivative(f(x), x) + 5*x**3*f(x)**2 - (2*x**2 + 2)*f(x)**3, 0)',
    'lg04': 'f(x).diff(x) - (f(x) - x*log(x))**2/x**2 + log(x)',
    'lg05': 'f(x).diff(x)**2',
    'lg11': 'x**2*f(x)**2 + x*Derivative(f(x), x)',
    'lg13': 'diff(f(x),x) + f(x)*cos(x) - exp(2*x)',
    'lg14': 'diff(f(x),x) + f(x)*cos(x) - sin(2*x)/2',
    'lg15': 'x*diff(f(x),x) + f(x) - x*sin(x)',
    'lg16': 'x*diff(f(x),x) - f(x) - x/log(x)',
    'lg17': '(f(x).diff(x)-f(x)) * (f(x).diff(x)+f(x))',
    'lg18': 'f(x).diff(x) * (f(x).diff(x) - f(x))',
    'lg20': 'f(x).diff(x)*(f(x).diff(x)+f(x))',
}

# dsolve(hint='lie_group') rows: test_single lie_group_01..20 (02, 06..10, 12
# and 19 repeat an equation above) and test_lie_group's other equations.
DSOLVE = ['lg01', 'lg03', 'lg04', 'lg05', 'h1_eq', 'h1_eq1', 'h1_eq2', 'h1_eq3',
          'h1_eq5', 'lg11', 'lg13', 'lg14', 'lg15', 'lg16', 'lg17', 'lg18', 'lg20',
          'h1_eq6', 'h1_eq7', 'h3_eq1', 'riccati', 'fsum', 'sim2']

SYMBOLS = 'x a b c a0 a1 a2 a3 a4 m n alpha'


def cases():
    out = []
    for name, eq in EQS.items():
        for hint in HINTS:
            out.append({'id': f'inf:{name}:{hint}', 'kind': 'inf', 'eq': name,
                        'code': f'infinitesimals({eq}, hint={hint!r})'})
        out.append({'id': f'check:{name}', 'kind': 'check', 'eq': name,
                    'code': f'checkinfsol({eq}, infinitesimals({eq}, hint="all"))'})
        out.append({'id': f'classify:{name}', 'kind': 'classify', 'eq': name,
                    'code': f'"lie_group" in classify_ode({eq}, f(x))'})
    for name in DSOLVE:
        eq = EQS[name]
        out.append({'id': f'dsolve:{name}', 'kind': 'dsolve', 'eq': name,
                    'code': f'dsolve({eq}, f(x), hint="lie_group")'})
    out.append({'id': 'dsolve:user', 'kind': 'dsolve', 'eq': 'h1_eq7',
                'code': 'dsolve(x*(f(x).diff(x)) + 1 - f(x)**2, f(x), hint="lie_group", '
                        'xi=sqrt(f(x) - 1)/sqrt(f(x) + 1), eta=0)'})
    out.append({'id': 'dsolve:user2', 'kind': 'dsolve', 'eq': 'h1_eq',
                'code': 'dsolve(Eq(f(x).diff(x), x**2*f(x)), f(x), hint="lie_group", xi=0, eta=f(x))'})
    return out


class _Alarm(Exception):
    pass


def _alarm(*_):
    raise _Alarm()


def render(value, srepr):
    if isinstance(value, dict):
        return {'dict': [[render(k, srepr), render(v, srepr)] for k, v in value.items()]}
    if isinstance(value, (list, tuple)):
        return {type(value).__name__: [render(v, srepr) for v in value]}
    if isinstance(value, bool):
        return value
    return srepr(value)


def namespace(mod):
    ns = dict(vars(mod))
    for s in SYMBOLS.split():
        ns[s] = mod.Symbol(s)
    ns['f'] = mod.Function('f')
    ns['F'] = mod.Function('F')
    from importlib import import_module
    ode = import_module(mod.__name__ + '.solvers.ode')
    ns['infinitesimals'] = ode.infinitesimals
    ns['checkinfsol'] = ode.checkinfsol
    return ns


def run(code, ns, srepr, limit):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(limit)
    start = time.perf_counter()
    try:
        out = {'ok': True, 'value': render(eval(code, ns), srepr)}
    except _Alarm:
        out = {'ok': False, 'error': 'DEADLINE', 'message': ''}
    except RecursionError:
        out = {'ok': False, 'error': 'RecursionError', 'message': ''}
    except Exception as error:
        out = {'ok': False, 'error': type(error).__name__, 'message': str(error)[:300]}
    finally:
        signal.alarm(0)
    out['seconds'] = round(time.perf_counter() - start, 3)
    return out


def isolated(code, ns, srepr, limit, counts):
    import os
    import select
    read, write = os.pipe()
    start = time.perf_counter()
    pid = os.fork()
    if pid == 0:
        os.close(read)
        try:
            out = run(code, ns, srepr, limit)
            data = json.dumps([out, dict(counts)]).encode()
        except BaseException as error:  # rendering failure
            data = json.dumps([{'ok': False, 'error': type(error).__name__,
                                'message': str(error)[:300], 'seconds': 0}, {}]).encode()
        with os.fdopen(write, 'wb') as pipe:
            pipe.write(data)
        os._exit(0)
    os.close(write)
    chunks = []
    deadline = start + limit + 10
    with os.fdopen(read, 'rb') as pipe:
        while True:
            left = deadline - time.perf_counter()
            if left <= 0 or not select.select([pipe], [], [], left)[0]:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
                return {'ok': False, 'error': 'DEADLINE', 'message': '',
                        'seconds': round(time.perf_counter() - start, 3)}, {}
            chunk = os.read(pipe.fileno(), 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
    _, status = os.waitpid(pid, 0)
    if not chunks:
        return {'ok': False, 'error': 'CRASH', 'message': str(status),
                'seconds': round(time.perf_counter() - start, 3)}, {}
    out, native = json.loads(b''.join(chunks))
    return out, native


def replay(backend, limit):
    if backend == 'symox':
        class NoReference:
            def find_spec(self, fullname, path=None, target=None):
                if fullname == 'sympy' or fullname.startswith('sympy.'):
                    raise ImportError('reference forbidden')
        sys.meta_path.insert(0, NoReference())
    mod = __import__(backend)
    ns = namespace(mod)
    counts = {}
    if backend == 'symox':
        from symox import _core
        for name in dir(_core):
            if name.startswith('ode_lie_') and name != 'ode_lie_group_remove':
                def wrap(fn, name=name):
                    def call(*args):
                        value = fn(*args)
                        key = name + (':declined' if value is None else ':native')
                        counts[key] = counts.get(key, 0) + 1
                        return value
                    return call
                setattr(_core, name, wrap(getattr(_core, name)))
    rows = []
    for case in cases():
        counts.clear()
        row = {'id': case['id'], 'kind': case['kind'], 'eq': case['eq']}
        if backend == 'symox':
            # SIGALRM does not interrupt native code: each symox case runs in
            # a forked child that is killed at the deadline.
            row['result'], row['native'] = isolated(case['code'], ns, mod.srepr, limit, counts)
        else:
            row['result'] = run(case['code'], ns, mod.srepr, limit)
        rows.append(row)
        print(case['id'], row['result'].get('error', 'ok'), row['result']['seconds'], file=sys.stderr, flush=True)
    return rows


def _dummies(text):
    """Number the ``Dummy`` indices of an ``srepr`` by first appearance."""
    import re
    seen = {}
    return re.sub(r'dummy_index=(\d+)',
                  lambda m: 'dummy_index=%d' % seen.setdefault(m.group(1), len(seen)), text)


def _normal(value):
    if isinstance(value, str):
        return _dummies(value)
    if isinstance(value, dict) and 'dict' in value:
        return {'dict': sorted(_dummies(json.dumps(p)) for p in value['dict'])}
    if isinstance(value, list):
        return [_normal(v) for v in value]
    if isinstance(value, dict):
        return {k: [_normal(v) for v in vs] for k, vs in value.items()}
    return value


def _unrender(value, ns):
    if isinstance(value, dict) and 'dict' in value:
        return {_unrender(k, ns): _unrender(v, ns) for k, v in value['dict']}
    if isinstance(value, dict):
        (kind, items), = value.items()
        seq = [_unrender(v, ns) for v in items]
        return tuple(seq) if kind == 'tuple' else seq
    if isinstance(value, bool):
        return value
    return eval(value, ns)


def report(ref_rows, sx_rows, limit):
    import sympy
    from sympy.solvers.ode import checkinfsol, checkodesol
    ns = namespace(sympy)
    by_id = {r['id']: r for r in ref_rows}
    out = {'counts': {}, 'by_kind': {}, 'cases': {}, 'native': {}}
    for row in sx_rows:
        ref = by_id[row['id']]['result']
        got = row['result']
        if ref['ok'] and got['ok']:
            verdict = 'SAME' if _normal(ref['value']) == _normal(got['value']) else None
        elif not ref['ok'] and not got['ok']:
            verdict = 'SAME' if ref['error'] == got['error'] else 'DIFF'
        elif ref['ok']:
            verdict = 'DECLINED'
        else:
            verdict = None
        if verdict is None:
            # symox answers differently: SymPy checks the answer.
            eq = eval(EQS[row['eq']], ns)
            value = _unrender(got['value'], ns)
            def check():
                if row['kind'] == 'inf':
                    return all(ok for ok, _ in checkinfsol(eq, value))
                if row['kind'] == 'dsolve':
                    sols = value if isinstance(value, list) else [value]
                    return all(checkodesol(eq, s, f=ns['f'](ns['x']))[0] for s in sols)
                return True
            res = run('check()', {'check': check}, sympy.srepr, limit)
            verdict = 'DIFF' if res['ok'] and res['value'] is True else (
                'WRONG' if res['ok'] else 'DIFF')
            checked = res['ok'] and res['value'] is True
        else:
            checked = None
        out['counts'][verdict] = out['counts'].get(verdict, 0) + 1
        kind = out['by_kind'].setdefault(row['kind'], {})
        kind[verdict] = kind.get(verdict, 0) + 1
        out['cases'][row['id']] = {'verdict': verdict, 'sympy': ref.get('error', 'ok'),
                                   'symox': got.get('error', 'ok'),
                                   'sympy_s': ref['seconds'], 'symox_s': got['seconds'],
                                   'native': row.get('native', {}),
                                   'validated_by_sympy': checked}
        for k, v in row.get('native', {}).items():
            out['native'][k] = out['native'].get(k, 0) + v
    return out


def main():
    mode = sys.argv[1]
    if mode in ('sympy', 'symox'):
        rows = replay(mode, 120)
        json.dump(rows, open(sys.argv[2], 'w'), indent=0)
    else:
        ref, sx = json.load(open(sys.argv[2])), json.load(open(sys.argv[3]))
        json.dump(report(ref, sx, 120), open(sys.argv[4], 'w'), indent=1)
        print(json.dumps(json.load(open(sys.argv[4]))['counts']))


if __name__ == '__main__':
    main()
