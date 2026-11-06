"""Harvest the differential pattern-matching corpus from SymPy 1.14.

usage: harvest.py WORKTREE OUT.json.gz

Three sources, all recorded as (expr, pattern, old) srepr strings that both
SymPy and symox evaluate:

* ``ode``: every ``.match(`` call made from ``sympy/solvers/ode/*.py`` while
  SymPy runs ``classify_ode`` and the archived hint of the 315 ODE corpus
  inputs (audits/sympy-1.14.0/corpus-2026-10-26/ode-cases.json.gz);
* ``test_match``: every ``match``/``matches`` call issued directly by the
  functions of ``sympy/core/tests/test_match.py``;
* ``random``: seeded variations (pattern instances, perturbed instances and
  unrelated expressions) over the ODE pattern shapes.

Wild ``exclude`` lists are kept (SymPy's srepr drops them).  Only the inputs
are stored here; ``compare.py`` evaluates them.
"""
import gzip, inspect, json, multiprocessing as mp, os, random, signal, sys
from pathlib import Path

import sympy
from sympy import Basic, Wild, srepr
from sympy.printing.repr import ReprPrinter

ODE_DIR = os.path.join(os.path.dirname(sympy.__file__), 'solvers', 'ode') + os.sep


class _Printer(ReprPrinter):
    def _print_Wild(self, e):
        args = [repr(e.name)]
        if e.exclude:
            args.append('exclude=(%s,)' % ', '.join(self._print(x) for x in e.exclude))
        if e.properties:
            raise ValueError('Wild properties are not serializable')
        return 'Wild(%s)' % ', '.join(args)


def wsrepr(e):
    return _Printer().doprint(e)


_records = []


def _install_ode_probe():
    original = Basic.match

    def match(self, pattern, old=False):
        caller = sys._getframe(1)
        if caller.f_code.co_filename.startswith(ODE_DIR):
            try:
                _records.append({'source': 'ode', 'caller': '%s:%s' % (
                    os.path.basename(caller.f_code.co_filename), caller.f_code.co_name),
                    'expr': wsrepr(self), 'pattern': wsrepr(pattern), 'old': bool(old)})
            except Exception:
                pass
        return original(self, pattern, old)
    Basic.match = match


class _Alarm(Exception):
    pass


def _alarm(*_):
    raise _Alarm()


def _limit():
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (6 * 1024 ** 3, 6 * 1024 ** 3))


def _ode_case(case):
    _records.clear()
    signal.signal(signal.SIGALRM, _alarm)
    ns = vars(sympy)
    for method in ('classify', 'hint'):
        signal.alarm(30)
        try:
            eq = eval(case['equation'], ns); func = eval(case['function'], ns)
            if method == 'classify':
                sympy.classify_ode(eq, func)
            else:
                sympy.dsolve(eq, func, hint=case['hint'], simplify=case['simplify_flag'])
        except (_Alarm, MemoryError, RecursionError):
            pass
        except Exception:
            pass
        finally:
            signal.alarm(0)
    return [dict(r, case=case['id']) for r in _records]


def harvest_ode(wt):
    cases = json.loads(gzip.open(Path(wt) / 'audits/sympy-1.14.0/corpus-2026-10-26/ode-cases.json.gz').read())['cases']
    _install_ode_probe()
    with mp.get_context('fork').Pool(24, initializer=_limit, maxtasksperchild=4) as pool:
        out = []
        for rows in pool.imap(_ode_case, cases, chunksize=1):
            out.extend(rows)
    return out


def harvest_test_match():
    from sympy.core import mul, add, power, operations, symbol, basic, function
    targets = [basic.Basic, operations.AssocOp, mul.Mul, add.Add, power.Pow, symbol.Wild,
               function.WildFunction, function.Function]
    records = []
    me = os.path.join('core', 'tests', 'test_match.py')

    def wrap(cls, name):
        original = cls.__dict__[name]

        def method(self, *args, **kw):
            caller = sys._getframe(1)
            if caller.f_code.co_filename.endswith(me):
                try:
                    if name == 'match':
                        expr, pattern = self, args[0]
                        old = args[1] if len(args) > 1 else kw.get('old', False)
                    else:
                        if (len(args) > 1 and args[1]) or kw.get('repl_dict'):
                            raise ValueError('repl_dict')
                        pattern, expr = self, args[0]
                        old = args[2] if len(args) > 2 else kw.get('old', False)
                    records.append({'source': 'test_match', 'caller': caller.f_code.co_name,
                                    'kind': name, 'expr': wsrepr(expr), 'pattern': wsrepr(pattern), 'old': bool(old)})
                except Exception:
                    pass
            return original(self, *args, **kw)
        setattr(cls, name, method)
        return original
    saved = []
    for cls in targets:
        for name in ('match', 'matches'):
            if name in cls.__dict__:
                saved.append((cls, name, wrap(cls, name)))
    import sympy.core.tests.test_match as tm
    for name, fn in sorted(vars(tm).items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn()
            except Exception:
                pass
    for cls, name, original in saved:
        setattr(cls, name, original)
    return records


def harvest_random(seed=20261106, count=500):
    from sympy import (Symbol, Function, Rational, Integer, sin, cos, exp, log, sqrt, I, pi,
                       Derivative)
    rng = random.Random(seed)
    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    f = Function('f'); fx = f(x); df = fx.diff(x); d2 = fx.diff(x, 2)
    W = lambda n, ex=(): Wild(n, exclude=list(ex))
    a, b, c, d = (W(n, [x]) for n in 'abcd')
    af, bf, cf = (W(n + 'f', [x, fx, df, d2]) for n in 'abc')
    u, v, w = W('u'), W('v'), W('w')
    atoms = [x, y, z, Integer(2), Integer(-3), Rational(1, 2), pi, I, sqrt(2)]
    funcs = [sin, cos, exp, log]

    def leaf(depth=0):
        r = rng.random()
        if depth > 2 or r < 0.45:
            return rng.choice(atoms)
        if r < 0.6:
            return rng.choice(funcs)(leaf(depth + 1))
        if r < 0.75:
            return leaf(depth + 1) * leaf(depth + 1)
        if r < 0.9:
            return leaf(depth + 1) + leaf(depth + 1)
        return leaf(depth + 1) ** rng.choice([Integer(2), Integer(-1), Rational(1, 2), y])

    patterns = [a * x + b, a * x ** b, a * x ** b + c, a + b * sin(c * x), a * exp(b * x),
                u * v, u + v, u * v + w, u ** v, sin(u) * v, a * x ** 2 + b * x + c,
                c ** 2 * x ** (2 * a), a * (x - b) ** c, a * x, a - b, d ** 2,
                af * d2 + bf * df + cf * fx, af * df + bf * fx ** 2 + cf * fx / x + d / x ** 2,
                a + b * x, af * d2 + bf * df ** 2 + cf * df, u * exp(v), log(u) + v,
                x ** a * exp(b), (a * x + b) ** c, u * sin(v) + w * cos(v)]
    records = []
    seen = set()
    while len(records) < count:
        p = rng.choice(patterns)
        kind = rng.random()
        try:
            wilds = sorted(p.atoms(Wild), key=str)
            if kind < 0.5:
                e = p.xreplace({q: leaf() for q in wilds})
                tag = 'instance'
            elif kind < 0.75:
                e = p.xreplace({q: leaf() for q in wilds}) + rng.choice([x, y, 1, x ** 2, sin(x)])
                tag = 'perturbed'
            else:
                e = leaf() + leaf() * rng.choice([x, fx, df, d2, x ** 2])
                tag = 'unrelated'
            key = (srepr(e), wsrepr(p))
        except Exception:
            continue
        if key in seen:
            continue
        seen.add(key)
        records.append({'source': 'random', 'caller': tag, 'expr': key[0], 'pattern': key[1], 'old': False})
    return records


def main():
    wt, out = sys.argv[1], sys.argv[2]
    rows = harvest_test_match()
    rows += harvest_random()
    rows += harvest_ode(wt)
    seen, unique = {}, []
    for r in rows:
        key = (r['expr'], r['pattern'], r['old'])
        if key in seen:
            seen[key]['count'] = seen[key].get('count', 1) + 1
            continue
        seen[key] = r
        unique.append(r)
    for i, r in enumerate(unique):
        r['id'] = i
    with gzip.open(out, 'wt') as fh:
        json.dump({'sympy': sympy.__version__, 'cases': unique}, fh)
    from collections import Counter
    print(len(rows), 'calls', len(unique), 'unique', Counter(r['source'] for r in unique))


if __name__ == '__main__':
    main()
