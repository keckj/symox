"""Statement-level run of SymPy 1.14 ``test_polytools.py`` on extensions.

Usage: python upstream_polytools.py sympy|symox OUT.jsonl

The test functions of ``sympy/polys/tests/test_polytools.py`` that mention
``extension`` or ``algebraic_field`` are run statement by statement (a
failing statement does not stop the function, unlike pytest); the
assertions that involve an extension, an algebraic field or a composite
ground are run and recorded (the other assertions are skipped; set-up
statements always run). For symox the module text is rewritten
(``sympy`` -> ``symox``) and each ``from symox... import name`` that symox
lacks becomes a stub that raises NotImplementedError when used.
"""
import ast, importlib, json, re, signal, sys

mode, out = sys.argv[1], sys.argv[2]
import sympy
path = sympy.__file__.rsplit('/', 1)[0] + '/polys/tests/test_polytools.py'
src = open(path).read()
if mode == 'symox':
    src = re.sub(r'\bsympy\b', 'symox', src)


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


tree = ast.parse(src)
ns = {}
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name.startswith('test_'):
        continue
    if isinstance(node, ast.ImportFrom) and node.module and node.module.split('.')[0] in ('sympy', 'symox'):
        for a in node.names:
            ns[a.asname or a.name] = _imp(node.module, a.name)
        continue
    exec(compile(ast.Module(body=[node], type_ignores=[]), path, 'exec'), ns)

KEY = re.compile(r'extension|algebraic_field|frac_field|QQ_I\[|ZZ_I\[|\bsqrt\(')


def alarm(*a):
    raise TimeoutError


signal.signal(signal.SIGALRM, alarm)
recs = []
for f in tree.body:
    if not (isinstance(f, ast.FunctionDef) and f.name.startswith('test_')):
        continue
    if not re.search(r'extension|algebraic_field', ast.get_source_segment(src, f)):
        continue
    loc = dict(ns)
    for k, st in enumerate(f.body):
        seg = ast.get_source_segment(src, st)
        # set-up statements always run; assertions only when they involve
        # an extension (the others, e.g. large factorizations, are skipped)
        if not KEY.search(seg) and not isinstance(st, (ast.Assign, ast.Expr, ast.With)):
            continue
        ok, err = True, None
        signal.alarm(30)
        try:
            exec(compile(ast.Module(body=[st], type_ignores=[]), path, 'exec'), loc)
        except BaseException as e:
            if isinstance(e, KeyboardInterrupt):
                raise
            ok, err = False, type(e).__name__ + ': ' + str(e)[:200]
        finally:
            signal.alarm(0)
        if KEY.search(seg) and isinstance(st, (ast.Assert, ast.With, ast.Expr)):
            recs.append({'case': f'{f.name}#{k}', 'ok': ok, 'error': err, 'src': seg[:300]})
with open(out, 'w') as fh:
    for r in recs:
        fh.write(json.dumps(r) + '\n')
print(f'{mode}: {sum(r["ok"] for r in recs)}/{len(recs)} extension statements pass')
