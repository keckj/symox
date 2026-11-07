"""Statement-level run of the SymPy 1.14 test modules of the lot.

Usage: python upstream.py sympy|symox OUT.jsonl

Modules: ``polys/tests/test_polytools.py`` (test functions about factor,
sqf, cancel, gcd, roots or extensions), ``polys/tests/test_factortools.py``,
``polys/tests/test_rootoftools.py`` and ``polys/numberfields/tests``
(``test_minpoly``, ``test_subfield``, ``test_numbers``, ``test_utilities``).
Each test function runs statement by statement (a failing statement does
not stop the function, unlike pytest); every assertion, ``raises`` block or
expression statement is recorded with its outcome. For symox the module
text is rewritten (``sympy`` -> ``symox``) and each ``from symox... import
name`` that symox lacks becomes a stub that raises NotImplementedError when
used, so a missing internal API shows as a failed statement.
"""
import ast, importlib, json, re, signal, sys

mode, out = sys.argv[1], sys.argv[2]
import sympy
root = sympy.__file__.rsplit('/', 1)[0]
MODULES = [
    ('polys/tests/test_polytools.py', r'factor|sqf|cancel|gcd|lcm|root|extension|algebraic'),
    ('polys/tests/test_factortools.py', r''),
    ('polys/tests/test_rootoftools.py', r''),
    ('polys/numberfields/tests/test_minpoly.py', r''),
    ('polys/numberfields/tests/test_subfield.py', r''),
    ('polys/numberfields/tests/test_numbers.py', r''),
    ('polys/numberfields/tests/test_utilities.py', r''),
]


class _Missing:
    def __init__(self, n):
        self.n = n

    def __call__(self, *a, **k):
        raise NotImplementedError('missing ' + self.n)

    def __getattr__(self, a):
        raise NotImplementedError('missing ' + self.n)


def _imp(m, n):
    try:
        return getattr(importlib.import_module(m), n)
    except Exception:
        return _Missing(m + '.' + n)


def alarm(*a):
    raise TimeoutError


signal.signal(signal.SIGALRM, alarm)
# Statements on which symox runs native code that the 30 s alarm cannot
# interrupt (more than 4 minutes at full CPU): recorded as failed without
# running them, for symox only.
HANG = {
    'test_polytools.py::test_factor_large#9',
    'test_polytools.py::test_cancel#80',
    'test_minpoly.py::test_minimal_polynomial_sq#9',
}
recs = []
fh = open(out, 'w')
for rel, pattern in MODULES:
    path = root + '/' + rel
    src = open(path).read()
    if mode == 'symox':
        src = re.sub(r'\bsympy\b', 'symox', src)
    tree = ast.parse(src)
    ns = {'__name__': 'upstream'}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name.startswith('test_'):
            continue
        if isinstance(node, ast.ImportFrom) and node.module and node.module.split('.')[0] in ('sympy', 'symox'):
            for a in node.names:
                if a.name == '*':
                    try:
                        mod = importlib.import_module(node.module)
                        ns.update({k: v for k, v in vars(mod).items() if not k.startswith('_')})
                    except Exception:
                        pass
                    continue
                ns[a.asname or a.name] = _imp(node.module, a.name)
            continue
        try:
            signal.alarm(60)
            exec(compile(ast.Module(body=[node], type_ignores=[]), path, 'exec'), ns)
        except Exception:
            pass
        finally:
            signal.alarm(0)
    for f in tree.body:
        if not (isinstance(f, ast.FunctionDef) and f.name.startswith('test_')):
            continue
        if pattern and not re.search(pattern, f.name):
            continue
        if any(isinstance(d, ast.Name) and d.id in ('slow', 'XFAIL', 'skip') or
               isinstance(d, ast.Call) and getattr(d.func, 'id', '') in ('skip', 'slow')
               for d in f.decorator_list):
            continue
        loc = dict(ns)
        for k, st in enumerate(f.body):
            seg = ast.get_source_segment(src, st) or ''
            ok, err = True, None
            case = f'{rel.rsplit("/", 1)[1]}::{f.name}#{k}'
            if mode == 'symox' and case in HANG:
                recs.append({'case': case, 'ok': False, 'error': 'HardTimeout: not interruptible, not run', 'src': seg[:300]})
                fh.write(json.dumps(recs[-1]) + '\n')
                fh.flush()
                continue
            signal.alarm(30)
            try:
                exec(compile(ast.Module(body=[st], type_ignores=[]), path, 'exec'), loc)
            except BaseException as e:
                if isinstance(e, KeyboardInterrupt):
                    raise
                ok, err = False, type(e).__name__ + ': ' + str(e)[:200]
            finally:
                signal.alarm(0)
            if isinstance(st, (ast.Assert, ast.With, ast.Expr, ast.For)):
                recs.append({'case': case, 'ok': ok,
                             'error': err, 'src': seg[:300]})
                fh.write(json.dumps(recs[-1]) + '\n')
                fh.flush()
fh.close()
print(f'{mode}: {sum(r["ok"] for r in recs)}/{len(recs)} statements pass')
