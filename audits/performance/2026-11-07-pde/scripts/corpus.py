"""Differential corpus of the native ``pde`` module (wave 45).

usage:
  corpus.py sympy OUT.json            # SymPy 1.14
  corpus.py symox OUT.json            # symox on PYTHONPATH (sympy blocked)
  corpus.py report SYMPY.json SYMOX.json OUT.json   # SymPy validates symox

Each case is a Python expression evaluated in the package namespace:
``classify_pde`` (tuple and dict), ``pdsolve`` for the default hint, every
matching hint, ``all``/``all_Integral``/``best`` and wrong hints,
``checkpdesol`` (of ``pdsolve`` results and of hand-written candidates) and
``pde_separate``/``_add``/``_mul``. The equations are those of SymPy's
``sympy/solvers/tests/test_pde.py`` and an in-house set (symbolic
coefficients, power reduction, ``Eq`` with a right-hand side, the error
paths of ``_desolve``, ``checkpdesol`` and ``pde_separate``).

A result is its ``srepr`` (dicts as ordered pair lists, exceptions as their
class and message) or the error class and message. The report classes each
case as
  SAME      equal value (dict order ignored), or equal error class;
  DECLINED  symox raises (or times out) where SymPy answers;
  WRONG     symox answers and SymPy rejects the answer (``checkpdesol`` of
            a ``pdsolve`` solution);
  DIFF      any other difference (``validated_by_sympy`` records whether
            SymPy's check accepted symox's differing answer).
``Dummy`` indices are numbered by first appearance before comparing.
For symox the run also counts the native calls (``_core.pde_*``) that
declined to the Python reference body.
"""
import json, signal, sys, time

EQS = {
    # test_pde.py
    'c1': 'a*u + b*u.diff(x) + c*u.diff(y)',
    'c2': '3*u + 2*u.diff(x) + u.diff(y)',
    'c3': 'a*u + b*u.diff(x) + 2*u.diff(y)',
    'v4': 'x*u + u.diff(x) + 3*u.diff(y)',
    'v5': 'x**2*u + x*u.diff(x) + x*y*u.diff(y)',
    'v6': 'y*x**2*u + y*u.diff(x) + u.diff(y)',
    'c5': '2*u + 1*u.diff(x) + 3*u.diff(y)',
    'c6': 'u + 1*u.diff(x) + 3*u.diff(y)',
    'sf1': 'u + u.diff(x) + u.diff(y)',
    'h1': '2*u + u.diff(x) + u.diff(y)',
    'h2': '4 + (3*u.diff(x)/u) + (2*u.diff(y)/u)',
    'h3': 'u + (6*u.diff(x)) + (7*u.diff(y))',
    'k1': '-2*u.diff(x) + 4*u.diff(y) + 5*u - exp(x + 3*y)',
    'k2': '(u.diff(x)/u) + (u.diff(y)/u) + 1 - (exp(x + y)/u)',
    'k3': '2*u + -u.diff(x) + 3*u.diff(y) + sin(x)',
    'k4': 'u + u.diff(x) + u.diff(y) + x*y',
    'k5': 'u + u.diff(x) + u.diff(y) + log(x)',
    'all1': 'u + u.diff(x) + u.diff(y) + x**2*y',
    'vc1': 'x*(u.diff(x)) - y*(u.diff(y)) + y**2*u - y**2',
    'vc2': 'x**2*u + x*u.diff(x) + x*y*u.diff(y)',
    'vc3': 'y*x**2*u + y*u.diff(x) + u.diff(y)',
    'vc4': 'exp(x)**2*(u.diff(x)) + y',
    'vc5': 'exp(2*x)*(u.diff(y)) + y*u - u',
    # in-house
    'm1': 'u.diff(x) + u.diff(y)',
    'm2': 'u.diff(x) - u.diff(y) + x',
    'm3': 'a*u.diff(x) + b*u.diff(y) - x*y',
    'm4': 'u.diff(x) + 2*u.diff(y) + 3*u - exp(x)*cos(y)',
    'm5': 'Eq(u.diff(x) + u.diff(y), x + y)',
    'm6': 'Eq(u.diff(x), u.diff(y))',
    'm7': 'u**2*u.diff(x) + u**3*u.diff(y) + u**3',
    'm8': 'u**2*u.diff(x) + u**2*u.diff(y) - x*u**2',
    'm9': 'x*u.diff(x) + u - y',
    'm10': 'y*u.diff(y) + u - x',
    'm11': '(x + y)*u.diff(x) + u.diff(y)',
    'm12': '(y + x)*u.diff(y) + u.diff(x)',
    'm13': 'x*u.diff(x) + y*u.diff(y) + x*y*u - x*y',
    'm14': 'u.diff(x) + x*u.diff(y)',
    'm15': 'u.diff(x) + I*u.diff(y)',
    'm16': 'u.diff(x) + I*u.diff(y) + 2*u - x',
    'm17': 'u.diff(x) + u.diff(y) - F(x - y)',
    'm18': 'u.diff(x) + u.diff(y) - xi',
    'm19': 'x*u.diff(x) - y*u.diff(y) + y*y*u - G(x, y)',
    'm20': 'u.diff(x)*u.diff(y)',
    'm21': 'u.diff(x, 2) + u.diff(y, 2)',
    'm22': 'u + x',
    'm23': 'sin(x)*u.diff(x) + cos(y)*u.diff(y)',
    'm24': 'u.diff(x) + y*u',
    'm25': 'x*u.diff(y) + u - exp(y)',
    'm26': '2*u.diff(x) + 3*u.diff(y) + 4*u',
    'm27': '2*u.diff(x) - u.diff(y) + u - x*y',
    'm28': 'u.diff(x) + u.diff(y) + u - exp(x + y)',
    'm29': 'x*u.diff(x) + y*u.diff(y)',
    'm30': 'u.diff(x) + exp(x)*u.diff(y) + u',
}

HINTS = ['1st_linear_constant_coeff_homogeneous', '1st_linear_constant_coeff',
         '1st_linear_constant_coeff_Integral', '1st_linear_variable_coeff']

SYMBOLS = 'a b c x y z t r theta xi eta'

# checkpdesol of hand-written candidates (pde, sol)
CHECKS = [
    ('m1', 'Eq(u, F(x - y))'),
    ('m1', 'F(x - y)'),
    ('m1', 'Eq(F(x - y), u)'),
    ('m1', 'Eq(x - y, u)'),
    ('m1', '[Eq(u, x - y), Eq(u, x + y)]'),
    ('m1', '{Eq(u, x - y), Eq(u, x + y)}'),
    ('m1', '(Eq(u, x - y), Eq(u, x + 2*y))'),
    ('m1', 'Eq(u**2, (x - y)**2)'),
    ('m1', 'Eq(u - x + y, 0)'),
    ('m1', 'Eq(sin(u), x - y)'),
    ('k1', 'x + y'),
    ('m2', 'Eq(u, F(x + y) - x**2/2)'),
    ('m3', 'Eq(u, x)'),
]

SEPARATE = [
    "pde_separate_add(Eq(D(w(x, t), x), D(w(x, t), t)*exp(w(x, t))), w(x, t), [X(x), T(t)])",
    "pde_separate(Eq(D(w(x, t), x), D(w(x, t), t)*exp(w(x, t))), w(x, t), [X(x), T(t)], 'div')",
    "pde_separate_mul(Eq(D(P(x, y, z), x) + D(P(x, y, z), y) + D(P(x, y, z), z), 0), P(x, y, z), [X(x), w(z, z)])",
    "pde_separate_mul(Eq(D(P(x, y, z), x) + D(P(x, y, z), y) + D(P(x, y, z), z), 0), P(x, y, z), [X(x), Y(y)])",
    "pde_separate_mul(Eq(D(P(x, y, z), x) + D(P(x, y, z), y) + D(P(x, y, z), z), 0), P(x, y, z), [X(t), Y(x, y)])",
    "pde_separate_mul(Eq(D(P(x, y, z), x) + D(P(x, y, z), y) + D(P(x, y, z), z), 0), P(x, y, z), [Y(y), w(x, z)])",
    "pde_separate_mul(Eq(D(P(x, y, z), x) + D(P(x, y, z), y) + D(P(x, y, z), z), 0), P(x, y, z), [X(x), Y(y), Z(z)])",
    "pde_separate_mul(Eq(D(w(x, t), t, t), C**2*D(w(x, t), x, x)), w(x, t), [X(x), T(t)])",
    "pde_separate_mul(Eq(1/r * D(Phi(r, theta, z), r) + D(Phi(r, theta, z), r, 2) + 1/r**2 * D(Phi(r, theta, z), theta, 2) + D(Phi(r, theta, z), z, 2), 0), Phi(r, theta, z), [Z(z), w(theta, r)])",
    "pde_separate_mul(Eq(-D(w(theta, r), r, r)/w(theta, r) - D(w(theta, r), r)/(r*w(theta, r)) - D(w(theta, r), theta, theta)/(r**2*w(theta, r)), C), w(theta, r), [T(theta), R(r)])",
    "pde_separate_mul(Eq(-D(w(theta, r), r, r)/w(theta, r) - D(w(theta, r), r)/(r*w(theta, r)) - D(w(theta, r), theta, theta)/(r**2*w(theta, r)), C), w(theta, r), [R(r), T(theta)])",
    "pde_separate(w(x, t).diff(x, 2) - w(x, t).diff(t, 2), w(x, t), [T(x), X(t)])",
    "pde_separate(w(x, t).diff(x, 2) - w(x, t).diff(t, 2), w(x, t), [T(x), X(t)], strategy='add')",
    "pde_separate(w(x, t).diff(x) + w(x, t).diff(t), w(x, t), [X(x), T(t)])",
    "pde_separate(w(x, t).diff(x) + w(x, t).diff(t), w(x, t), [X(x), T(t)], strategy='add')",
    "pde_separate(w(x, t).diff(x)*t + w(x, t).diff(t)*x, w(x, t), [X(x), T(t)])",
    "pde_separate(w(x, t).diff(x) - x*t*w(x, t).diff(t), w(x, t), [X(x), T(t)])",
    "pde_separate(w(x, t).diff(x)*w(x, t).diff(t), w(x, t), [X(x), T(t)])",
    "pde_separate(Eq(w(x, t).diff(x), w(x, t).diff(t) + w(x, t)), w(x, t), [X(x), T(t)])",
    "pde_separate(w(x, t).diff(x), w(x, t), [X(x), T(x, t)])",
    "pde_separate(x*w(x, t), w(x, t), [X(x), T(t)])",
    "pde_separate(w(x, t).diff(x, t), w(x, t), [X(x), T(t)])",
    "pde_separate(w(x, t).diff(x, 2) + w(x, t).diff(t) + w(x, t), w(x, t), [X(x), T(t)])",
]


def cases():
    out = []
    for name, eq in EQS.items():
        out.append({'id': f'classify:{name}', 'kind': 'classify', 'eq': name,
                    'code': f'classify_pde({eq})'})
        out.append({'id': f'classify_dict:{name}', 'kind': 'classify', 'eq': name,
                    'code': f'classify_pde({eq}, u, dict=True)'})
        out.append({'id': f'pdsolve:{name}', 'kind': 'pdsolve', 'eq': name,
                    'code': f'pdsolve({eq})'})
        out.append({'id': f'pdsolve_func:{name}', 'kind': 'pdsolve', 'eq': name,
                    'code': f'pdsolve({eq}, u, solvefun=G)'})
        for hint in HINTS:
            out.append({'id': f'pdsolve:{name}:{hint}', 'kind': 'pdsolve', 'eq': name,
                        'code': f'pdsolve({eq}, u, hint={hint!r})'})
        for hint in ('all', 'all_Integral', 'best', 'foo'):
            out.append({'id': f'pdsolve:{name}:{hint}', 'kind': 'all', 'eq': name,
                        'code': f'pdsolve({eq}, u, hint={hint!r})'})
        out.append({'id': f'check:{name}', 'kind': 'check', 'eq': name,
                    'code': f'checkpdesol({eq}, pdsolve({eq}, u), u)'})
    # direct methods
    for name in ('h1', 'c1', 'k1', 'k4', 'm3', 'vc1', 'vc3', 'm9', 'm11'):
        eq = EQS[name]
        for hint in HINTS:
            out.append({'id': f'method:{name}:{hint}', 'kind': 'method', 'eq': name,
                        'code': f'pde_{hint.removesuffix("_Integral")}({eq}, u, 1, '
                                f'classify_pde({eq}, u, dict=True)[{hint!r}], H)'})
    # argument errors
    extra = [
        ('err:3var', 'classify_pde(P(x, y, z).diff(x), P(x, y, z))'),
        ('err:1var', 'classify_pde(f(x).diff(x) + f(x), f(x))'),
        ('err:nofunc', 'classify_pde(x + y)'),
        ('err:twofunc', 'pdsolve(u.diff(x) + g(x, y).diff(y))'),
        ('err:eq_dict', 'classify_pde(Eq(u.diff(x), x), dict=True)'),
        ('err:eq_zero_dict', 'classify_pde(Eq(u.diff(x) + u.diff(y), 0), dict=True)'),
        ('err:pow', 'pdsolve((u.diff(x) + u.diff(y))**2, u)'),
        ('err:hint_pow', 'pdsolve((u.diff(x) + u.diff(y) - x)**3, u, hint="1st_linear_constant_coeff")'),
        ('err:orderkw', 'pdsolve(u.diff(x) + u.diff(y), u, ics={u: 1})'),
        ('err:derivexpr', 'pdsolve(Derivative(x*u, x) + u.diff(y), u)'),
        ('err:derivexpr_c', 'classify_pde(Derivative(x*u, x) + u.diff(y), u)'),
        ('err:sfun', 'pdsolve(u.diff(x) + u.diff(y) + u, u, solvefun=Function("K"))'),
        ('err:checkfunc', 'checkpdesol(x + y, Eq(u, x))'),
        ('err:checkfunc2', 'checkpdesol(x + y, Eq(u, x), u)'),
        ('err:check_nosolve', 'checkpdesol(x*u + u.diff(x) + 3*u.diff(y), Eq(u - F(3*x - y), 0), u, solve_for_func=False)'),
        ('err:check_sol4', 'checkpdesol(x*u + u.diff(x) + 3*u.diff(y), Eq(pdsolve(x*u + u.diff(x) + 3*u.diff(y)).lhs - pdsolve(x*u + u.diff(x) + 3*u.diff(y)).rhs, 0), solve_for_func=False)'),
        ('check:list45', 'checkpdesol(x*u + u.diff(x) + 3*u.diff(y), [pdsolve(2*u + 1*u.diff(x) + 3*u.diff(y)), pdsolve(u + 1*u.diff(x) + 3*u.diff(y))])'),
    ]
    for cid, code in extra:
        out.append({'id': cid, 'kind': 'extra', 'eq': '', 'code': code})
    for i, (name, sol) in enumerate(CHECKS):
        out.append({'id': f'checkhand:{i}:{name}', 'kind': 'check', 'eq': name,
                    'code': f'checkpdesol({EQS[name]}, {sol}, u)'})
        out.append({'id': f'checkhand_nofunc:{i}:{name}', 'kind': 'check', 'eq': name,
                    'code': f'checkpdesol({EQS[name]}, {sol})'})
    for i, code in enumerate(SEPARATE):
        out.append({'id': f'separate:{i}', 'kind': 'separate', 'eq': '', 'code': code})
    return out


class _Alarm(Exception):
    pass


def _alarm(*_):
    raise _Alarm()


def render(value, srepr):
    if isinstance(value, BaseException):
        return {'exc': [type(value).__name__, str(value)[:300]]}
    if isinstance(value, dict):
        return {'dict': [[render(k, srepr), render(v, srepr)] for k, v in value.items()]}
    if isinstance(value, (list, tuple, set)):
        items = [render(v, srepr) for v in value]
        if isinstance(value, set):
            items = sorted(items, key=json.dumps)
        return {type(value).__name__: items}
    if isinstance(value, (bool, str)) or value is None:
        return value if not isinstance(value, str) else {'str': value}
    if isinstance(value, int):
        return {'int': value}
    return srepr(value)


def namespace(mod):
    ns = dict(vars(mod))
    for s in SYMBOLS.split():
        ns[s] = mod.Symbol(s)
    ns['C'] = mod.Symbol('C', real=True)
    for name in 'f g F G H K P R T X Y Z w Phi'.split():
        ns[name] = mod.Function(name)
    ns['u'] = ns['f'](ns['x'], ns['y'])
    ns['D'] = mod.Derivative
    from importlib import import_module
    pde = import_module(mod.__name__ + '.solvers.pde')
    for name in dir(pde):
        if name.startswith('pde_') or name in ('pdsolve', 'classify_pde', 'checkpdesol'):
            ns[name] = getattr(pde, name)
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
            if name.startswith('pde_'):
                def wrap(fn, name=name):
                    def call(*args, **kwargs):
                        value = fn(*args, **kwargs)
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
        return {k: (_normal(vs) if isinstance(vs, list) else vs) for k, vs in value.items()}
    return value


def _unrender(value, ns):
    if isinstance(value, dict) and 'dict' in value:
        return {_unrender(k, ns): _unrender(v, ns) for k, v in value['dict']}
    if isinstance(value, dict):
        (kind, items), = value.items()
        if kind in ('str', 'int'):
            return items
        if kind == 'exc':
            return NotImplementedError(items[1])
        seq = [_unrender(v, ns) for v in items]
        return tuple(seq) if kind == 'tuple' else seq
    if isinstance(value, bool) or value is None:
        return value
    return eval(value, ns)


def report(ref_rows, sx_rows, limit):
    import sympy
    from sympy.solvers.pde import checkpdesol
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
            # symox answers differently: SymPy checks pdsolve solutions.
            def check():
                if row['kind'] != 'pdsolve' or not row['eq']:
                    return None
                eq = eval(EQS[row['eq']], ns)
                value = _unrender(got['value'], ns)
                if not isinstance(value, sympy.Equality):
                    return False
                return bool(checkpdesol(eq, value.doit(), ns['u'])[0])
            res = run('check()', {'check': check}, sympy.srepr, limit)
            checked = res['ok'] and res['value'] is True
            verdict = 'WRONG' if res['ok'] and res['value'] is False else 'DIFF'
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
        rows = replay(mode, 60)
        json.dump(rows, open(sys.argv[2], 'w'), indent=0)
    else:
        ref, sx = json.load(open(sys.argv[2])), json.load(open(sys.argv[3]))
        json.dump(report(ref, sx, 60), open(sys.argv[4], 'w'), indent=1)
        print(json.dumps(json.load(open(sys.argv[4]))['counts']))


if __name__ == '__main__':
    main()
