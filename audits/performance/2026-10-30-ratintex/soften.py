import ast, sys
src = open(sys.argv[1]).read()
tree = ast.parse(src)
lines = src.splitlines(keepends=True)
out = []; pos = 0
offs = [0]
for l in lines: offs.append(offs[-1] + len(l))
head = ("import importlib as _il\n"
        "def _missing(name):\n"
        "    def f(*a, **k): raise NotImplementedError(name + ' not provided by symox')\n"
        "    return f\n"
        "def _soft(mod, name):\n"
        "    try: m = _il.import_module(mod)\n"
        "    except ImportError: return _missing(name)\n"
        "    return getattr(m, name, None) or _missing(name)\n")
edits = []
for node in tree.body:
    if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith('symox'):
        start = offs[node.lineno - 1]; end = offs[node.end_lineno]
        rep = ''.join(f"{a.asname or a.name} = _soft({node.module!r}, {a.name!r})\n" for a in node.names)
        edits.append((start, end, rep))
res = src
for s, e, r in reversed(edits): res = res[:s] + r + res[e:]
# insert head after docstring
i = res.index('\n', res.index('"""', 3)) + 1 if res.startswith('"""') else 0
open(sys.argv[2], 'w').write(res[:i] + head + res[i:])
