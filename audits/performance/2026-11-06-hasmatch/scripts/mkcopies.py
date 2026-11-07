"""mkcopies.py OUTDIR: SymPy 1.14 test files with sympy imports rewritten to
symox (a missing name raises NotImplementedError when used)."""
import ast, os, sys, glob
import sympy
OUT = sys.argv[1]
ROOT = os.path.dirname(sympy.__file__)
FILES = ['simplify/tests/test_trigsimp.py', 'simplify/tests/test_simplify.py', 'simplify/tests/test_fu.py',
         'simplify/tests/test_radsimp.py', 'integrals/tests/test_meijerint.py', 'integrals/tests/test_transforms.py',
         'integrals/tests/test_laplace.py', 'core/tests/test_basic.py', 'core/tests/test_match.py']
os.makedirs(OUT, exist_ok=True)
for rel in FILES:
    src = open(os.path.join(ROOT, rel)).read()
    tree = ast.parse(src)
    lines = src.split("\n")
    out = ["def _missing(name):\n    def f(*a, **k):\n        if len(a) == 1 and callable(a[0]) and not k:\n            return f  # a decorator: the decorated test fails when run\n        raise NotImplementedError('missing ' + name)\n    return f\n"]
    last = 0
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("sympy"):
            out.extend(lines[last:node.lineno - 1])
            mod = node.module.replace("sympy", "symox")
            for a in node.names:
                nm = a.asname or a.name
                out.append("try:\n    from %s import %s as %s\nexcept ImportError:\n    try:\n        from symox import %s as %s\n    except ImportError:\n        %s = _missing(%r)" % (mod, a.name, nm, a.name, nm, nm, nm))
            last = node.end_lineno
    out.extend(lines[last:])
    name = 'test_' + rel.split('/')[0] + '_' + os.path.basename(rel)[5:]
    open(os.path.join(OUT, name), "w").write("\n".join(out).replace("sympy", "symox"))
open(os.path.join(OUT, 'conftest.py'), 'w').write('''
import signal, pytest
class _Deadline(Exception):
    pass
def _alarm(*a):
    raise _Deadline()
@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(120)
    try:
        yield
    finally:
        signal.alarm(0)
''')
