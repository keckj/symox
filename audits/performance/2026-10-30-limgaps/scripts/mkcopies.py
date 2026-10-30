import ast, sys, glob, os
T = sys.argv[1]
for f in ["limits", "gruntz", "series", "nseries", "order", "formal", "aseries", "lseries", "demidovich"]:
    src = open(os.path.join(T, "test_%s.py" % f)).read()
    tree = ast.parse(src)
    lines = src.split("\n")
    out = ["def _missing(name):\n    def f(*a, **k):\n        raise NotImplementedError('missing ' + name)\n    return f\n"]
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
    open("copies/test_%s.py" % f, "w").write("\n".join(out).replace("sympy", "symox"))
