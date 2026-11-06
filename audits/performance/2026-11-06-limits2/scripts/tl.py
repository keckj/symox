"""Statement-level evaluation of SymPy series test files.

tl.py LIB OUT FILE...   (LIB is sympy or symox)

For each test function, the body statements are executed in order in a fresh
namespace (the module header imported with LIB substituted for sympy).  For
every ``assert lhs == rhs`` (or ``lhs is rhs``) whose source mentions limit,
series, O( / Order, gruntz or leadterm, the lhs is evaluated under a 20 s
alarm and its str (or the exception class) is recorded under file:line.
"""
import ast, json, signal, sys, time, os, warnings
warnings.simplefilter("ignore")

LIB, OUT = sys.argv[1], sys.argv[2]
FILES = sys.argv[3:]
KEYS = ("limit", "series", "O(", "Order", "gruntz", "leadterm", "Limit", "as_leading_term", "lseries", "heuristics")


class Alarm(Exception):
    pass


def _alarm(*a):
    raise Alarm()


signal.signal(signal.SIGALRM, _alarm)


def header_ns(tree, src, fname):
    ns = {"__name__": "tl_" + os.path.basename(fname)}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name.startswith("test"):
            continue
        code = ast.get_source_segment(src, node)
        if LIB == "symox":
            code = code.replace("sympy", "symox")
        try:
            exec(compile(code, fname, "exec"), ns)
        except BaseException:
            pass
    return ns


def main():
    out = {}
    core = None
    stats = None
    if LIB == "symox":
        from symox import _core as core
        if hasattr(core, "_limit_trace"):
            core._limit_trace(True)
        else:
            core = None
        stats = _count_accum_paths()
    for fname in FILES:
        src = open(fname).read()
        tree = ast.parse(src)
        base = header_ns(tree, src, fname)
        short = os.path.basename(fname)
        for fn in tree.body:
            if not isinstance(fn, ast.FunctionDef) or not fn.name.startswith("test"):
                continue
            ns = dict(base)
            idx = 0
            for st in fn.body:
                seg = ast.get_source_segment(src, st) or ""
                if LIB == "symox":
                    seg_exec = seg.replace("sympy", "symox")
                else:
                    seg_exec = seg
                if isinstance(st, ast.Assert) and isinstance(st.test, ast.Compare) and len(st.test.ops) == 1 \
                        and isinstance(st.test.ops[0], (ast.Eq, ast.Is)) and any(k in seg for k in KEYS):
                    lhs = ast.get_source_segment(src, st.test.left)
                    if LIB == "symox":
                        lhs = lhs.replace("sympy", "symox")
                    t = time.perf_counter()
                    signal.alarm(20)
                    try:
                        v = eval(compile(lhs, fname, "eval"), ns)
                        r = str(v)
                    except Alarm:
                        r = "TIMEOUT"
                    except BaseException as e:
                        r = "EXC " + type(e).__name__
                    finally:
                        signal.alarm(0)
                    idx += 1
                    out["%s:%s:%d" % (short, fn.name, idx)] = [r, round(time.perf_counter() - t, 4), st.lineno]
                    continue
                signal.alarm(20)
                try:
                    exec(compile(seg_exec, fname, "exec"), ns)
                except BaseException:
                    pass
                finally:
                    signal.alarm(0)
    json.dump(out, open(OUT, "w"), indent=0, sort_keys=True)
    if core is not None:
        json.dump(core._limit_trace(False), open(OUT[:-5] + ".trace.json", "w"), indent=0)
    if stats is not None:
        json.dump(stats, open(OUT[:-5] + ".stats.json", "w"), indent=0, sort_keys=True)


def _count_accum_paths():
    """Count how the limits declined by the native ``_limit`` end in
    ``oscillating_limit``: before the port, the native bounded decomposition
    (``_core._limit_bounded_affine``), the Python walker, or held; after it
    (``_core._oscillating_limit`` present), the native bounded decomposition,
    the native walker, or held."""
    import symox.calculus._accum_limit as AL
    from symox import _core
    native = hasattr(_core, "_oscillating_limit")
    walker = "native_walker" if native else "python_walker"
    stats = {"native_bounded_affine": 0, walker: 0, "held": 0}
    hits = [0]
    depth = [0]
    orig_osc = AL.oscillating_limit

    if native:
        def affine(e, z, z0, dir):
            # The native result came from the bounded decomposition when it
            # applies on that side (both sides for '+-').
            try:
                sides = ('+', '-') if str(dir) == '+-' else (str(dir),)
                return all(_core._limit_bounded_affine(e, z, z0, d) is not None for d in sides)
            except Exception:
                return False
    else:
        orig_ba = AL._bounded_affine

        def ba(*a):
            r = orig_ba(*a)
            if r is not None:
                hits[0] += 1
            return r
        AL._bounded_affine = ba

    def osc(*a):
        depth[0] += 1
        before = hits[0]
        try:
            r = orig_osc(*a)
        finally:
            depth[0] -= 1
        if depth[0] == 0:
            if r is None:
                key = "held"
            elif native:
                key = "native_bounded_affine" if affine(*a) else walker
            else:
                key = "native_bounded_affine" if hits[0] > before else walker
            stats[key] += 1
        return r

    AL.oscillating_limit = osc
    return stats


main()
