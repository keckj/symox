"""Statement-level replay of sympy/series/tests/test_formal.py.

tl.py LIB OUT TESTFILE      (LIB is sympy or symox)

Every ``assert`` (also inside loops) is rewritten into a recording call:
for ``assert lhs == rhs`` (or ``is`` / ``!=``) the str of ``lhs`` and the
truth of the whole assert are recorded; for any other assert, the truth.
Every ``raises(E, f)`` records the exception class raised by ``f()`` (or
``NO EXC <str>``). Each top-level statement of a test body runs under a 20 s
alarm. For symox, the source has ``sympy`` replaced by ``symox``.
"""
import ast, json, os, signal, sys, time, warnings
warnings.simplefilter("ignore")
sys.setrecursionlimit(10000)

LIB, OUT, FNAME = sys.argv[1:4]


class Alarm(Exception):
    pass


def _alarm(*a):
    raise Alarm()


signal.signal(signal.SIGALRM, _alarm)
REC = {}


def _show(v):
    if isinstance(v, (list, tuple)):
        return type(v).__name__ + "(" + ", ".join(map(_show, v)) + ")"
    return str(v)


def _rec(key, line, lhs_f, test_f):
    t = time.perf_counter()
    try:
        r = _show(lhs_f()) if lhs_f is not None else None
    except Alarm:
        raise
    except BaseException as e:
        r = "EXC " + type(e).__name__
    try:
        ok = bool(test_f())
    except Alarm:
        raise
    except BaseException as e:
        ok = "EXC " + type(e).__name__
    REC[key] = [r, ok, round(time.perf_counter() - t, 4), line]


def _rec_raises(key, line, exc, f):
    t = time.perf_counter()
    try:
        v = f()
        r = "NO EXC " + _show(v)[:80]
    except Alarm:
        raise
    except BaseException as e:
        r = "EXC " + type(e).__name__
    ok = r == "EXC " + exc.__name__
    REC[key] = [r, ok, round(time.perf_counter() - t, 4), line]


class Rewrite(ast.NodeTransformer):
    def __init__(self, fn):
        self.fn, self.idx = fn, 0

    def key(self):
        self.idx += 1
        return "%s:%d" % (self.fn, self.idx)

    def visit_Assert(self, node):
        t = node.test
        lhs = t.left if isinstance(t, ast.Compare) and len(t.ops) == 1 else None
        lam = lambda e: ast.Lambda(args=ast.arguments(posonlyargs=[], args=[], kwonlyargs=[], kw_defaults=[], defaults=[]), body=e)
        call = ast.Call(func=ast.Name("_rec", ast.Load()),
                        args=[ast.Constant(self.key()), ast.Constant(node.lineno),
                              lam(lhs) if lhs is not None else ast.Constant(None), lam(t)], keywords=[])
        return ast.copy_location(ast.Expr(call), node)

    def visit_Expr(self, node):
        v = node.value
        if isinstance(v, ast.Call) and isinstance(v.func, ast.Name) and v.func.id == "raises" and len(v.args) == 2:
            call = ast.Call(func=ast.Name("_rec_raises", ast.Load()),
                            args=[ast.Constant(self.key()), ast.Constant(node.lineno), v.args[0], v.args[1]], keywords=[])
            return ast.copy_location(ast.Expr(call), node)
        return self.generic_visit(node)


def main():
    src = open(FNAME).read()
    if LIB == "symox":
        src = src.replace("sympy", "symox")
    tree = ast.parse(src)
    ns = {"__name__": "tl_formal", "_rec": _rec, "_rec_raises": _rec_raises}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name.startswith("test"):
            continue
        if isinstance(node, ast.ImportFrom):
            for a in node.names:
                one = ast.ImportFrom(module=node.module, names=[a], level=0)
                try:
                    exec(compile(ast.fix_missing_locations(ast.Module([one], [])), FNAME, "exec"), ns)
                except BaseException:
                    nm = a.asname or a.name

                    def missing(*args, _nm=nm, **kw):
                        raise NotImplementedError("missing " + _nm)
                    ns[nm] = missing
            continue
        try:
            exec(compile(ast.Module([node], []), FNAME, "exec"), ns)
        except BaseException:
            pass
    for fn in tree.body:
        if not isinstance(fn, ast.FunctionDef) or not fn.name.startswith("test"):
            continue
        rw = Rewrite(fn.name)
        body = [ast.fix_missing_locations(rw.visit(st)) for st in fn.body]
        loc = dict(ns)
        for st in body:
            signal.alarm(20)
            try:
                exec(compile(ast.Module([st], []), FNAME, "exec"), loc)
            except Alarm:
                REC["%s:stmt%d" % (fn.name, st.lineno)] = ["TIMEOUT", False, 20, st.lineno]
            except BaseException:
                pass
            finally:
                signal.alarm(0)
    json.dump(REC, open(OUT, "w"), indent=0, sort_keys=True)


main()
