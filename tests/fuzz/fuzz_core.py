"""Differential fuzzer: symox versus SymPy.

A seeded random expression generator builds the *same* expression tree in
SymPy and in symox through their Python APIs and compares the results of
many operations (printing, srepr round trips, expand, diff, subs, evalf,
free_symbols, equality/hash, pickling, ``func(*args)`` reconstruction) as
well as the higher level algorithms (integrate, solve, series, limit,
summation, simplify/trigsimp, lambdify, the parser).

Every mismatch is classified:

* ``crash``    symox raised (or panicked) where SymPy did not,
* ``died``     the worker process died (abort, stack overflow, segfault),
* ``timeout``  symox did not finish within the per case time limit,
* ``slow``     symox finished but took far longer than expected,
* ``wrong``    numerically different values (a real bug),
* ``form``     same value, different canonical form (printing/structure),
* ``weak``     symox gave up (unevaluated) where SymPy succeeded.

Failing cases are shrunk by greedy delta debugging on the expression tree
and deduplicated; findings go to ``<out>/findings.jsonl``.

Usage::

    python tests/fuzz/fuzz_core.py --mode core --jobs 2 --time 600
    python tests/fuzz/fuzz_core.py --mode all --seeds 0:20000
    python tests/fuzz/fuzz_core.py --replay core 1234       # one case, verbose

Each case runs in a separate long lived worker process (restarted after a
crash or timeout), so a panic or stack overflow in symox never kills the
fuzzer.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import multiprocessing as mp
import os
import pickle
import queue
import random
import sys
import time
import traceback

# ---------------------------------------------------------------------------
# expression trees
# ---------------------------------------------------------------------------
#
# ('S', name)            symbol from the pool
# ('N', int)             Integer
# ('Q', p, q)            Rational
# ('F', float)           Float
# ('C', name)            constant: pi E I oo -oo zoo nan EulerGamma
# ('+'|'-'|'*'|'/'|'**', a, b)
# ('neg', a)
# ('f', name, a, ...)    function application
# ('rel', op, a, b)      relational (top level only)

SYMS = {
    "x": {},
    "y": {},
    "z": {},
    "r": {"real": True},
    "p": {"positive": True},
    "q": {"negative": True},
    "n": {"integer": True},
    "m": {"integer": True, "positive": True},
    "k": {"nonnegative": True},
    "w": {"nonzero": True},
}

UNARY = [
    "exp", "log", "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos",
    "atan", "acot", "sinh", "cosh", "tanh", "coth", "asinh", "acosh", "atanh",
    "Abs", "sign", "floor", "ceiling", "gamma", "loggamma", "erf", "erfc",
    "re", "im", "conjugate", "arg", "sqrt", "Heaviside",
]
COMMON_UNARY = ["exp", "log", "sin", "cos", "tan", "sqrt", "Abs", "atan", "sinh", "cosh"]
BINARY_FUNCS = ["atan2", "Max", "Min", "log2"]
RELOPS = ["Eq", "Ne", "Lt", "Le", "Gt", "Ge"]
CONSTS = ["pi", "E", "I", "oo", "-oo", "zoo", "nan", "EulerGamma"]
FLOATS = [0.5, 1.25, -2.75, 0.1, 3.3, 2.0, -0.5, 1e-3, 1.5e5]


def tree_size(t):
    if t[0] in ("S", "N", "Q", "F", "C"):
        return 1
    if t[0] == "f":
        return 1 + sum(tree_size(a) for a in t[2:])
    if t[0] == "rel":
        return 1 + tree_size(t[2]) + tree_size(t[3])
    return 1 + sum(tree_size(a) for a in t[1:])


def to_src(t):
    """Python source (valid in both namespaces) for a tree."""
    h = t[0]
    if h == "S":
        return t[1]
    if h == "N":
        return f"Integer({t[1]})"
    if h == "Q":
        return f"Rational({t[1]}, {t[2]})"
    if h == "F":
        return f"Float({t[1]!r})"
    if h == "C":
        return "(-oo)" if t[1] == "-oo" else t[1]
    if h == "neg":
        return f"(-{to_src(t[1])})"
    if h in ("+", "-", "*", "/", "**"):
        return f"({to_src(t[1])} {h} {to_src(t[2])})"
    if h == "f":
        if t[1] == "log2":
            return f"log({to_src(t[2])}, {to_src(t[3])})"
        return f"{t[1]}({', '.join(to_src(a) for a in t[2:])})"
    if h == "rel":
        return f"{t[1]}({to_src(t[2])}, {to_src(t[3])})"
    raise ValueError(t)


class Lib:
    """One library (sympy or symox) plus its symbol pool."""

    def __init__(self, mod):
        self.mod = mod
        self.name = mod.__name__
        self.syms = {n: mod.Symbol(n, **a) for n, a in SYMS.items()}
        self.ns = {n: getattr(mod, n) for n in dir(mod) if not n.startswith("_")}
        self.ns.update(self.syms)
        self.consts = {
            "pi": mod.pi, "E": mod.E, "I": mod.I, "oo": mod.oo, "-oo": -mod.oo,
            "zoo": mod.zoo, "nan": mod.nan, "EulerGamma": mod.EulerGamma,
        }

    def build(self, t):
        m = self.mod
        h = t[0]
        if h == "S":
            return self.syms[t[1]]
        if h == "N":
            return m.Integer(t[1])
        if h == "Q":
            return m.Rational(t[1], t[2])
        if h == "F":
            return m.Float(t[1])
        if h == "C":
            return self.consts[t[1]]
        if h == "neg":
            return -self.build(t[1])
        if h in ("+", "-", "*", "/", "**"):
            a, b = self.build(t[1]), self.build(t[2])
            if h == "+":
                return a + b
            if h == "-":
                return a - b
            if h == "*":
                return a * b
            if h == "/":
                return a / b
            return a ** b
        if h == "f":
            if t[1] == "log2":
                return m.log(self.build(t[2]), self.build(t[3]))
            return getattr(m, t[1])(*[self.build(a) for a in t[2:]])
        if h == "rel":
            return getattr(m, t[1])(self.build(t[2]), self.build(t[3]))
        raise ValueError(t)


# ---------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------


class Gen:
    def __init__(self, rng, syms=None, funcs=None, consts=True, floats=True, inf=True):
        self.rng = rng
        self.syms = syms or list(SYMS)
        self.funcs = funcs or UNARY
        self.consts = consts
        self.floats = floats
        self.inf = inf

    def atom(self):
        r = self.rng
        u = r.random()
        if u < 0.45:
            return ("S", r.choice(self.syms))
        if u < 0.70:
            return ("N", r.choice([0, 1, 2, 3, -1, -2, 4, 5, 6, 8, 9, 12, -3, 10, 100, 2**70]) if r.random() < 0.97 else r.randint(-10**30, 10**30))
        if u < 0.82:
            p = r.choice([1, -1, 2, 3, -3, 5, 7])
            q = r.choice([2, 3, 4, 5, 7, 12])
            return ("Q", p, q) if p % q else ("N", p // q)
        if u < 0.88 and self.floats:
            return ("F", r.choice(FLOATS))
        if self.consts:
            c = r.choice(CONSTS[:3] * 6 + (CONSTS[3:] if self.inf else []) + ["EulerGamma"])
            return ("C", c)
        return ("S", r.choice(self.syms))

    def exponent(self, depth):
        r = self.rng
        u = r.random()
        if u < 0.45:
            return ("N", r.choice([2, 3, -1, -2, 2, 3, 4, -3, 0, 5]))
        if u < 0.75:
            p, q = r.choice([(1, 2), (-1, 2), (1, 3), (2, 3), (3, 2), (-3, 2), (-1, 3), (5, 2), (1, 4)])
            return ("Q", p, q)
        if u < 0.85 and self.floats:
            return ("F", r.choice([0.5, 2.0, -1.5, 0.3]))
        return self.expr(min(depth, 1))

    def expr(self, depth):
        r = self.rng
        if depth <= 0 or r.random() < 0.18:
            return self.atom()
        u = r.random()
        if u < 0.24:
            return ("+", self.expr(depth - 1), self.expr(depth - 1))
        if u < 0.32:
            return ("-", self.expr(depth - 1), self.expr(depth - 1))
        if u < 0.52:
            return ("*", self.expr(depth - 1), self.expr(depth - 1))
        if u < 0.58:
            return ("/", self.expr(depth - 1), self.expr(depth - 1))
        if u < 0.72:
            return ("**", self.expr(depth - 1), self.exponent(depth - 1))
        if u < 0.75:
            return ("neg", self.expr(depth - 1))
        if u < 0.95:
            f = r.choice(self.funcs if r.random() < 0.5 else COMMON_UNARY)
            return ("f", f, self.expr(depth - 1))
        f = r.choice(BINARY_FUNCS)
        return ("f", f, self.expr(depth - 1), self.expr(depth - 1))


def gen_core(rng):
    g = Gen(rng)
    t = g.expr(rng.choice([2, 3, 3, 4, 4, 5]))
    if rng.random() < 0.04:
        t = ("rel", rng.choice(RELOPS), t, g.expr(2))
    return t


# ---------------------------------------------------------------------------
# numerics
# ---------------------------------------------------------------------------

_sympy = None
_symox = None
SP = None  # Lib(sympy)
SX = None  # Lib(symox)


def init_libs():
    global _sympy, _symox, SP, SX
    if SP is None:
        import sympy
        import symox
        _sympy, _symox = sympy, symox
        SP, SX = Lib(sympy), Lib(symox)


def rand_point(rng, names):
    """Exact (Rational/Integer, possibly complex) values respecting assumptions."""
    pt = {}
    for n in names:
        a = SYMS.get(n, {})

        def rr(lo, hi):
            return (rng.randint(int(lo * 100), int(hi * 100)) or 37) / 100.0

        if a.get("integer"):
            v = rng.randint(1, 4) if a.get("positive") else rng.randint(-3, 3)
        elif a.get("positive"):
            v = rr(0.1, 2.5)
        elif a.get("negative"):
            v = -rr(0.1, 2.5)
        elif a.get("nonnegative"):
            v = abs(rr(0, 2.5))
        elif a.get("real"):
            v = rr(-2.5, 2.5)
        else:
            v = complex(rr(-2, 2), rr(-2, 2))
        pt[n] = v
    return pt


def sp_num(v):
    S = _sympy
    if isinstance(v, complex):
        return S.Rational(round(v.real * 100), 100) + S.I * S.Rational(round(v.imag * 100), 100)
    if isinstance(v, int):
        return S.Integer(v)
    return S.Rational(round(v * 100), 100)


def sx_num(v):
    X = _symox
    if isinstance(v, complex):
        return X.Rational(round(v.real * 100), 100) + X.I * X.Rational(round(v.imag * 100), 100)
    if isinstance(v, int):
        return X.Integer(v)
    return X.Rational(round(v * 100), 100)


INF = complex(float("inf"), 0)
NAN = complex(float("nan"), 0)


def to_complex(v):
    """sympy number -> python complex (inf/nan aware); None if not a number."""
    S = _sympy
    if v is S.nan:
        return NAN
    if v is S.zoo:
        return NAN  # an unsigned complex infinity has no numerical direction
    if v == -S.oo:
        return -INF
    if v == S.oo:
        return INF
    try:
        return complex(v)
    except (TypeError, ValueError):
        if getattr(v, "is_number", False) and v.has(S.oo, -S.oo, S.zoo):
            return INF
        if getattr(v, "is_number", False) and v.has(S.nan):
            return NAN
        return None
    except OverflowError:
        return INF


def sp_value(e, pt, prec=30):
    """Value of a sympy expression at an exact point (None if not numeric)."""
    S = _sympy
    sub = {SP.syms[n]: sp_num(v) for n, v in pt.items() if n in SP.syms}
    try:
        v = e.xreplace(sub)
        if not getattr(v, "is_number", False):
            return None
        v = v.evalf(prec)
    except (ZeroDivisionError, OverflowError, ValueError, TypeError, RecursionError, AttributeError, KeyError, NotImplementedError):
        return None
    return to_complex(v)


def close(a, b, tol):
    if cmath.isnan(a) or cmath.isnan(b):
        return cmath.isnan(a) and cmath.isnan(b)
    if cmath.isinf(a) or cmath.isinf(b):
        return a == b
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def sx_to_sp(e):
    """symox expression -> sympy (via srepr: structure + assumptions)."""
    return _sympy.sympify(_symox.srepr(e))


def values_differ(s, xs, rng, npts=3, tol=1e-9, names=None):
    """Compare two sympy expressions numerically; returns a description of
    the first differing point or None."""
    if names is None:
        names = sorted({str(a) for a in s.free_symbols | xs.free_symbols})
    if any(n not in SYMS for n in names):
        return None
    hits = 0
    for _ in range(npts):
        pt = rand_point(rng, names)
        a = sp_value(s, pt)
        b = sp_value(xs, pt)
        if a is None or b is None:
            continue
        hits += 1
        if close(a, b, tol):
            continue
        # confirm at higher precision (sympy's evalf may lose digits)
        a2 = sp_value(s, pt, 60)
        b2 = sp_value(xs, pt, 60)
        if a2 is None or b2 is None or close(a2, b2, tol):
            continue
        return f"at {pt}: sympy={a2} symox={b2}"
    return None


def has_float(t):
    if t[0] == "F":
        return True
    return any(isinstance(a, tuple) and has_float(a) for a in t[1:])


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------


class Finding(Exception):
    def __init__(self, cls, msg):
        super().__init__(msg)
        self.cls = cls
        self.msg = msg


class SympySkip(Exception):
    """SymPy itself failed: the case is not usable."""


STATE = None  # shared memory with the parent: what is running right now


def stage(s):
    if STATE is not None:
        b = s.encode("utf8", "replace")[:1000]
        STATE.value = b


def sp_call(f, *a):
    stage("sympy")
    try:
        return f(*a)
    except Exception as e:  # noqa: BLE001
        raise SympySkip(f"{type(e).__name__}: {e}") from None
    except RecursionError as e:
        raise SympySkip(str(e)) from None


SLOW = 2.0  # seconds: a symox call slower than that is reported


def sx_call(what, f, *a):
    stage("symox:" + what)
    t0 = time.perf_counter()
    try:
        r = f(*a)
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as e:  # noqa: BLE001  (PanicException is a BaseException)
        name = type(e).__name__
        raise Finding("crash", f"{what}: {name}: {str(e)[:300]}") from None
    dt = time.perf_counter() - t0
    if dt > SLOW:
        raise Finding("slow", f"{what}: {dt:.1f}s")
    return r


def compare_exprs(what, s, x, rng, tolerant=False, strict_str=True):
    """s sympy result, x symox result."""
    ss = str(s)
    xs = sx_call("str", str, x)
    if ss == xs:
        return
    try:
        xsp = sx_to_sp(x)
    except Exception as e:  # noqa: BLE001
        raise Finding("form", f"{what}: srepr of symox result not sympifiable: {e}; symox={xs} sympy={ss}") from None
    if isinstance(s, _sympy.Basic) and xsp == s:
        if strict_str:
            raise Finding("form", f"{what}: same tree, different str: symox={xs} sympy={ss}")
        return
    if isinstance(s, _sympy.Basic) and not s.has(_sympy.Derivative, _sympy.Integral) and not getattr(s, "is_Relational", False) and not isinstance(s, _sympy.logic.boolalg.Boolean):
        d = values_differ(s, xsp, rng, tol=1e-6 if tolerant else 1e-9)
        if d:
            raise Finding("wrong", f"{what}: symox={xs} sympy={ss} {d}")
    raise Finding("form", f"{what}: symox={xs} sympy={ss}")


def built(t):
    s = sp_call(SP.build, t)
    x = sx_call("build", SX.build, t)
    return s, x


def chk_build(t, rng):
    s, x = built(t)
    compare_exprs("build", s, x, rng, tolerant=has_float(t))


def chk_srepr(t, rng):
    s, x = built(t)
    r = sx_call("srepr", _symox.srepr, x)
    x2 = sx_call("sympify(srepr(x))", _symox.sympify, r)
    if x2 != x:
        difference = values_differ(sx_to_sp(x), sx_to_sp(x2), rng)
        raise Finding("wrong" if difference else "form", f"srepr roundtrip: {r} -> {_symox.srepr(x2)} {difference or ''}")
    if hash(x2) != hash(x):
        raise Finding("wrong", f"hash differs after srepr roundtrip: {r}")
    sr = sp_call(_sympy.srepr, s)
    x3 = sx_call("sympify(sympy.srepr)", _symox.sympify, sr)
    if str(s) == str(x) and x3 != x:
        raise Finding("form", f"sympify(sympy.srepr(e)) != e: {sr} -> {x3} (expected {x})")
    compare_exprs("sympify(sympy.srepr)", s, x3, rng, tolerant=has_float(t))
    # the str printer must parse back (plain symbol pool: assumptions are lost in str)
    xs = str(x)
    try:
        sp_back = _sympy.parse_expr(xs, local_dict=dict(SP.syms))
    except Exception as e:  # noqa: BLE001
        raise Finding("form", f"str(symox) not parseable by sympy: {xs}: {e}") from None
    if not isinstance(sp_back, _sympy.Basic) or getattr(sp_back, "is_Relational", False) or isinstance(sp_back, _sympy.logic.boolalg.Boolean):
        return
    d = values_differ(sx_to_sp(x), sp_back, rng)
    if d:
        raise Finding("wrong", f"str printer changes the value: str={xs} srepr={r} {d}")
    x4 = sx_call("sympify(str)", _symox.sympify, xs)
    if x4 != x and all(not SYMS.get(str(v)) for v in x.free_symbols):
        raise Finding("form", f"sympify(str(e)) != e: {xs}")


def chk_eqhash(t, rng):
    _, x = built(t)
    x2 = sx_call("build", SX.build, t)
    if not (x == x2) or hash(x) != hash(x2):
        raise Finding("wrong", f"rebuild differs: {x} vs {x2}")
    if x != x:
        raise Finding("wrong", "x != x")
    y = sx_call("pickle", lambda: pickle.loads(pickle.dumps(x)))
    if y != x or hash(y) != hash(x) or str(y) != str(x):
        raise Finding("wrong", f"pickle roundtrip: {_symox.srepr(x)} -> {_symox.srepr(y)}")
    # func(*args) reconstruction, recursively
    stack = [x]
    seen = 0
    while stack and seen < 200:
        e = stack.pop()
        seen += 1
        args = sx_call("args", lambda: e.args)
        if not args:
            continue
        if e.is_Relational or e.class_name in ("BooleanTrue", "BooleanFalse"):
            continue
        r = sx_call("func(*args)", lambda: e.func(*args))
        if r != e:
            raise Finding("form", f"func(*args) != e for {e} ({e.class_name}): got {r} args={args}")
        stack.extend(args)


def chk_free(t, rng):
    s, x = built(t)
    a = {str(v) for v in s.free_symbols}
    b = {str(v) for v in sx_call("free_symbols", lambda: x.free_symbols)}
    if a != b:
        sx = str(x)
        if str(s) == sx:
            raise Finding("wrong", f"free_symbols: symox={sorted(b)} sympy={sorted(a)} of {sx}")
        # canonical forms differ, the difference is reported by chk_build
        bb = {str(v) for v in sx_to_sp(x).free_symbols}
        if bb != b:
            raise Finding("wrong", f"free_symbols: symox={sorted(b)} actual={sorted(bb)} of {sx}")


def chk_expand(t, rng):
    s, x = built(t)
    es = sp_call(_sympy.expand, s)
    ex = sx_call("expand", _symox.expand, x)
    compare_exprs("expand", es, ex, rng, tolerant=has_float(t))


def chk_diff(t, rng):
    s, x = built(t)
    if getattr(s, "is_Relational", False):
        return
    if any(arg.args[0].is_real is not True for arg in s.atoms(_sympy.Abs)):
        return  # documented real-variable Abs derivative is outside this oracle
    v = rng.choice(list(SYMS))
    ds = sp_call(_sympy.diff, s, SP.syms[v])
    dx = sx_call(f"diff(., {v})", _symox.diff, x, SX.syms[v])
    if ds.has(_sympy.Abs, _sympy.sign, _sympy.re, _sympy.im, _sympy.arg, _sympy.conjugate) or s.has(_sympy.Abs, _sympy.sign, _sympy.re, _sympy.im, _sympy.arg, _sympy.conjugate, _sympy.floor, _sympy.ceiling, _sympy.Heaviside):
        # real-variable formulas are a documented deviation
        if ds.has(_sympy.Derivative):
            return
    compare_exprs(f"diff(., {v})", ds, dx, rng, tolerant=has_float(t))


def chk_subs(t, rng):
    s, x = built(t)
    v = rng.choice(list(SYMS))
    g = Gen(rng)
    rep = g.expr(rng.choice([0, 0, 1, 2]))
    rs, rx = built(rep)
    ss = sp_call(lambda: s.subs(SP.syms[v], rs))
    xx = sx_call(f"subs({v}, {to_src(rep)})", lambda: x.subs(SX.syms[v], rx))
    compare_exprs(f"subs({v}, {to_src(rep)})", ss, xx, rng, tolerant=has_float(t) or has_float(rep))


def chk_evalf(t, rng):
    s, x = built(t)
    if getattr(s, "is_Relational", False) or isinstance(s, _sympy.logic.boolalg.Boolean):
        return
    names = sorted(str(a) for a in s.free_symbols)
    pt = rand_point(rng, names)
    ref = sp_value(s, pt, 30)
    if ref is None:
        return
    if s.has(_sympy.zoo):
        return  # complex infinity has no direction for a numerical comparison
    sub = {SX.syms[n]: sx_num(v) for n, v in pt.items()}
    xv = sx_call("subs(point)", lambda: x.subs(sub) if sub else x)
    xe = sx_call("evalf", lambda: xv.evalf())
    try:
        val = complex(xe)
    except (TypeError, ValueError):
        sxe = str(xe)
        if sxe in ("oo", "-oo", "zoo") or "oo" in sxe:
            val = INF
        elif sxe == "nan" or "nan" in sxe:
            val = NAN
        else:
            raise Finding("weak", f"evalf did not give a number: {sxe} (sympy {ref}) at {pt}") from None
    except OverflowError:
        val = INF
    # double precision: cancellations are allowed to lose digits
    if not close(ref, val, 1e-6):
        # suppress catastrophic cancellation: evaluate the sympy reference at 15 digits too
        ref15 = sp_value(s, pt, 15)
        if ref15 is not None and close(ref15, val, 1e-6):
            return
        raise Finding("wrong", f"evalf at {pt}: symox={val} sympy={ref} (symox subs gave {xv})")


CORE_CHECKS = {
    "build": chk_build,
    "srepr": chk_srepr,
    "eqhash": chk_eqhash,
    "free": chk_free,
    "expand": chk_expand,
    "diff": chk_diff,
    "subs": chk_subs,
    "evalf": chk_evalf,
}

# ---------------------------------------------------------------------------
# higher level algorithms
# ---------------------------------------------------------------------------

X_ONLY = ["x", "x", "x", "p", "y"]


def gen_univariate(rng, depth=None, funcs=None, syms=("x",), inf=False):
    g = Gen(rng, syms=list(syms), funcs=funcs or ["exp", "log", "sin", "cos", "tan", "sqrt", "atan", "sinh", "cosh", "asin", "tanh"], consts=True, floats=False, inf=inf)
    return g.expr(depth if depth is not None else rng.choice([1, 2, 2, 3, 3]))


def gen_integrand(rng):
    """Integrands from families with elementary antiderivatives."""
    x = ("S", "x")
    u = rng.random()
    g = Gen(rng, syms=["x"], floats=False, inf=False, consts=False)
    if u < 0.25:
        # polynomial times elementary
        poly = g.expr(1)
        f = rng.choice(["exp", "sin", "cos", "sinh", "cosh"])
        a = ("N", rng.choice([1, 2, 3, -1, -2]))
        return ("*", poly, ("f", f, ("*", a, x)))
    if u < 0.45:
        # rational function
        num = g.expr(2)
        den = ("+", ("**", x, ("N", rng.choice([1, 2, 3]))), ("N", rng.choice([1, 2, -1, 4, -4])))
        return ("/", num, den)
    if u < 0.55:
        return ("*", ("**", ("f", rng.choice(["sin", "cos"]), x), ("N", rng.randint(1, 4))), ("**", ("f", rng.choice(["sin", "cos"]), x), ("N", rng.randint(0, 3))))
    if u < 0.65:
        return ("*", ("**", x, ("N", rng.randint(-1, 3))), ("f", "log", x))
    if u < 0.75:
        return ("**", ("+", ("**", x, ("N", 2)), ("N", rng.choice([1, -1, 2]))), ("Q", rng.choice([1, -1, -3]), 2))
    return gen_univariate(rng, funcs=["exp", "log", "sin", "cos", "sqrt", "atan"])


def chk_integrate(t, rng):
    s, x = built(t)
    X = SX.syms["x"]
    F = sx_call("integrate", lambda: _symox.integrate(x, X))
    if "Integral" in str(F):
        # symox gave up; does sympy do it quickly? (weak, low priority)
        return
    dF = sx_call("diff(F)", lambda: _symox.diff(F, X))
    Fs = sx_to_sp(F)
    d = values_differ(s, sx_to_sp(dF), rng, npts=3, tol=1e-7)
    if d:
        raise Finding("wrong", f"diff(integrate(f)) != f: F={F} {d}")
    if Fs.has(_sympy.nan, _sympy.zoo):
        raise Finding("wrong", f"antiderivative contains nan/zoo: {F}")


def gen_solve(rng):
    x = ("S", "x")
    u = rng.random()
    g = Gen(rng, syms=["y", "p"], floats=False, inf=False, consts=False)
    if u < 0.5:
        deg = rng.randint(1, 4)
        t = ("N", rng.randint(-5, 5))
        for i in range(1, deg + 1):
            c = ("N", rng.choice([1, -1, 2, -2, 3, 0, 5])) if rng.random() < 0.8 else g.atom()
            if i == deg and c == ("N", 0):
                c = ("N", 1)
            t = ("+", t, ("*", c, ("**", x, ("N", i))))
        return t
    if u < 0.65:
        # product of linear factors
        t = ("N", 1)
        for _ in range(rng.randint(1, 3)):
            t = ("*", t, ("-", x, g.atom() if rng.random() < 0.3 else ("N", rng.randint(-4, 4))))
        return t
    if u < 0.8:
        f = rng.choice(["exp", "log", "sin", "cos", "tan", "sqrt", "atan", "sinh"])
        return ("-", ("f", f, ("+", ("*", ("N", rng.choice([1, 2, -1])), x), ("N", rng.randint(-2, 2)))), ("Q", rng.choice([1, -1, 3]), rng.choice([1, 2])) if rng.random() < 0.7 else g.atom())
    if u < 0.9:
        return ("-", ("**", x, ("Q", rng.choice([1, 3, 2]), rng.choice([2, 3]))), ("N", rng.randint(-3, 5)))
    return ("+", ("/", ("N", 1), ("-", x, ("N", rng.randint(-2, 2)))), ("/", ("N", rng.randint(1, 3)), ("+", x, ("N", rng.randint(-2, 2)))))


def chk_solve(t, rng):
    s, x = built(t)
    X = SX.syms["x"]
    try:
        sols = sx_call("solve", lambda: _symox.solve(x, X))
    except Finding as f:
        if "NotImplementedError" in f.msg:
            return
        raise
    for sol in sols:
        res = sx_call("subs", lambda: x.subs(X, sol))
        rs = sx_to_sp(res)
        names = sorted(str(a) for a in rs.free_symbols)
        for _ in range(3):
            pt = rand_point(rng, names)
            v = sp_value(rs, pt)
            if v is None:
                continue
            if cmath.isnan(v) or cmath.isinf(v) or abs(v) > 1e-8 * max(1.0, _scale(s, sx_to_sp(sol), pt)):
                raise Finding("wrong", f"solution {sol} does not satisfy the equation: residual {v} at {pt}")
    # compare the number of solutions with sympy (polynomials only: sympy is fast there)
    if s.is_polynomial(SP.syms["x"]) and not s.free_symbols - {SP.syms["x"]}:
        ss = sp_call(_sympy.solve, s, SP.syms["x"])
        if len(ss) != len(sols):
            sset = {complex(_sympy.N(v)) for v in ss}
            xset = {complex(_sympy.N(sx_to_sp(v))) for v in sols}
            if len(sset) != len(xset) or any(min(abs(a - b) for b in xset) > 1e-8 for a in sset):
                raise Finding("weak" if len(sols) < len(ss) else "wrong", f"solutions: symox={sols} sympy={ss}")
        else:
            try:
                compare_exprs("solve", _sympy.FiniteSet(*ss), _symox_finite(sols), rng, strict_str=False)
            except Finding:
                if [str(v) for v in ss] != [str(v) for v in sols]:
                    raise Finding("form", f"solve: symox={[str(v) for v in sols]} sympy={[str(v) for v in ss]}") from None


def _symox_finite(sols):
    class _F:
        def __str__(self):
            return "{" + ", ".join(str(v) for v in sols) + "}"
    return _F()


def _scale(s, sol, pt):
    try:
        v = sp_value(sol, pt)
        return 1.0 + abs(v) ** max(1, int(_sympy.degree(s, SP.syms["x"]))) if v is not None else 1.0
    except Exception:  # noqa: BLE001
        return 1.0


def gen_series(rng):
    return gen_univariate(rng, funcs=["exp", "log", "sin", "cos", "tan", "sqrt", "atan", "sinh", "cosh", "asin", "tanh", "acos", "erf", "gamma", "asinh", "atanh"], syms=("x", "x", "y"))


def chk_series(t, rng):
    s, x = built(t)
    n = rng.randint(2, 6)
    try:
        ss = sp_call(lambda: _sympy.series(s, SP.syms["x"], 0, n))
    except SympySkip:
        return
    xs = sx_call(f"series(n={n})", lambda: _symox.series(x, SX.syms["x"], 0, n))
    if str(ss) == str(xs):
        return
    xsp = sx_to_sp(xs)
    a, b = ss.removeO(), xsp.removeO()
    oa, ob = ss.getO(), xsp.getO()
    d = values_differ(a, b, rng, npts=2, tol=1e-8)
    if d:
        raise Finding("wrong", f"series(n={n}): symox={xs} sympy={ss} {d}")
    if (oa is None) != (ob is None) or (oa is not None and oa != ob):
        raise Finding("wrong" if oa is not None and ob is not None and str(oa) != str(ob) else "form", f"series(n={n}) order term: symox={xs} sympy={ss}")
    raise Finding("form", f"series(n={n}): symox={xs} sympy={ss}")


def gen_limit(rng):
    return gen_univariate(rng, funcs=["exp", "log", "sin", "cos", "tan", "sqrt", "atan", "sinh", "cosh", "tanh", "asin"], syms=("x", "x", "x", "p"))


def chk_limit(t, rng):
    s, x = built(t)
    pt = rng.choice([0, 0, "oo", "-oo", 1, 2, "pi"])
    d = rng.choice(["+", "-"]) if pt not in ("oo", "-oo") else ("-" if pt == "oo" else "+")
    pts = {"oo": SP.mod.oo, "-oo": -SP.mod.oo, "pi": SP.mod.pi}.get(pt, pt)
    ptx = {"oo": SX.mod.oo, "-oo": -SX.mod.oo, "pi": SX.mod.pi}.get(pt, pt)
    ls = sp_call(lambda: _sympy.limit(s, SP.syms["x"], pts, d))
    if isinstance(ls, _sympy.Limit) or ls.has(_sympy.Limit) or ls.has(_sympy.AccumBounds):
        return
    lx = sx_call(f"limit(x->{pt}{d})", lambda: _symox.limit(x, SX.syms["x"], ptx, d))
    if isinstance(lx, _symox.Limit):
        raise Finding("weak", f"limit(x->{pt}{d}) unevaluated, sympy={ls}")
    compare_exprs(f"limit(x->{pt}{d})", ls, lx, rng)


def gen_sum(rng):
    k = ("S", "m")
    u = rng.random()
    g = Gen(rng, syms=["m", "y"], floats=False, inf=False, consts=False, funcs=["exp"])
    if u < 0.4:
        t = ("N", rng.randint(-3, 3))
        for i in range(1, rng.randint(1, 4)):
            t = ("+", t, ("*", ("N", rng.choice([1, -1, 2, 3])), ("**", k, ("N", i))))
        return t
    if u < 0.6:
        return ("*", g.expr(1), ("**", rng.choice([("N", 2), ("Q", 1, 2), ("S", "y"), ("Q", -1, 3)]), k))
    if u < 0.75:
        a = rng.randint(0, 2)
        return ("/", ("N", 1), ("*", ("+", k, ("N", a)), ("+", k, ("N", a + rng.randint(1, 2)))))
    if u < 0.85:
        return ("/", ("N", 1), ("**", k, ("N", rng.choice([2, 3, 4]))))
    return ("*", ("**", k, ("N", rng.randint(0, 2))), ("**", ("N", rng.choice([-1, 2, 3])), k))


def chk_sum(t, rng):
    s, x = built(t)
    lo = rng.choice([0, 1, 1, 2])
    hi_kind = rng.choice(["num", "sym", "oo"])
    K = SX.syms["m"]
    Ks = SP.syms["m"]
    if hi_kind == "num":
        hi = rng.randint(lo, lo + 6)
        if lo == 0 and s.subs(Ks, 0).has(_sympy.zoo, _sympy.nan, _sympy.oo):
            return
        ref = sp_call(lambda: _sympy.Add(*[s.subs(Ks, i) for i in range(lo, hi + 1)]))
        r = sx_call(f"summation(m, {lo}, {hi})", lambda: _symox.summation(x, (K, lo, hi)))
        compare_exprs(f"summation(m={lo}..{hi})", _sympy.simplify(ref), r, rng, strict_str=False)
        return
    if lo == 0 and t[0] == "/":
        lo = 1
    if hi_kind == "oo":
        ref = sp_call(lambda: _sympy.summation(s, (Ks, lo, _sympy.oo)))
        if ref.has(_sympy.Sum):
            return
        r = sx_call(f"summation(m, {lo}, oo)", lambda: _symox.summation(x, (K, lo, _symox.oo)))
        if "Sum" in str(r):
            raise Finding("weak", f"summation(m={lo}..oo) unevaluated, sympy={ref}")
        compare_exprs(f"summation(m={lo}..oo)", ref, r, rng, strict_str=False)
        return
    N = SX.syms["p"]
    r = sx_call(f"summation(m, {lo}, p)", lambda: _symox.summation(x, (K, lo, N)))
    if "Sum" in str(r):
        return
    rs = sx_to_sp(r)
    for hi in range(lo, lo + 5):
        ref = _sympy.Add(*[s.subs(Ks, i) for i in range(lo, hi + 1)])
        got = rs.subs(SP.syms["p"], hi)
        d = values_differ(ref, got, rng, npts=1, names=sorted(str(a) for a in (ref.free_symbols | got.free_symbols)))
        if d:
            raise Finding("wrong", f"summation(m={lo}..p) = {r} wrong at p={hi}: {d}")


def gen_simplify(rng):
    trig = rng.random() < 0.6
    funcs = ["sin", "cos", "tan", "sinh", "cosh", "cot"] if trig else ["exp", "log", "sqrt", "sin", "cos"]
    g = Gen(rng, syms=["x", "y", "p"], funcs=funcs, floats=False, inf=False)
    t = g.expr(rng.choice([2, 3, 3, 4]))
    if trig and rng.random() < 0.5:
        a = ("S", rng.choice(["x", "y"]))
        t = ("+", t, ("+", ("**", ("f", "sin", a), ("N", 2)), ("**", ("f", "cos", a), ("N", 2))))
    return t


def chk_simplify(t, rng):
    s, x = built(t)
    for name in ("simplify", "trigsimp", "together", "radsimp", "powsimp", "expand_trig", "expand_log", "logcombine", "expand_mul", "nsimplify"):
        if name == "nsimplify":
            continue
        f = getattr(_symox, name)
        r = sx_call(name, lambda: f(x))
        d = values_differ(s, sx_to_sp(r), rng, npts=3)
        if d:
            raise Finding("wrong", f"{name}: {x} -> {r} {d}")


def gen_lambdify(rng):
    g = Gen(rng, syms=["x", "y", "p"], funcs=["exp", "log", "sin", "cos", "tan", "sqrt", "atan", "sinh", "cosh", "tanh", "Abs", "asin", "acos", "erf", "sign", "floor", "ceiling", "gamma", "atan2", "Max", "Min"], floats=True, inf=False, consts=True)
    g.consts = True
    t = g.expr(rng.choice([2, 3, 4]))
    return t


def chk_lambdify(t, rng):
    import numpy as np
    s, x = built(t)
    if s.has(_sympy.I) or s.has(_sympy.zoo, _sympy.nan, _sympy.oo):
        return
    args_s = [SP.syms[n] for n in ("x", "y", "p")]
    args_x = [SX.syms[n] for n in ("x", "y", "p")]
    fs = sp_call(lambda: _sympy.lambdify(args_s, s, modules="numpy"))
    fx = sx_call("lambdify", lambda: _symox.lambdify(args_x, x))
    N = 16
    xs = np.array([rng.uniform(-3, 3) for _ in range(N)])
    ys = np.array([rng.uniform(-3, 3) for _ in range(N)])
    ps = np.array([rng.uniform(0.05, 3) for _ in range(N)])
    with np.errstate(all="ignore"):
        try:
            a = np.broadcast_to(np.asarray(fs(xs, ys, ps), dtype=complex), (N,))
        except Exception:  # noqa: BLE001
            return
    b = sx_call("lambdified call", lambda: fx(xs, ys, ps))
    b = np.broadcast_to(np.asarray(b, dtype=complex), (N,))
    for i in range(N):
        va, vb = complex(a[i]), complex(b[i])
        if cmath.isnan(va) or cmath.isinf(va):
            continue  # numpy domain error (e.g. sqrt(-1) is nan for real arrays)
        if abs(va.imag) > 0:
            continue
        if not close(va, vb, 1e-8):
            # sympy's own lambdify may be less accurate: check the exact value
            ref = sp_value(s, {"x": float(xs[i]), "y": float(ys[i]), "p": float(ps[i])})
            if ref is not None and close(ref, va, 1e-8) and not close(ref, vb, 1e-6):
                raise Finding("wrong", f"lambdify at x={xs[i]!r}, y={ys[i]!r}, p={ps[i]!r}: symox={vb} numpy={va} exact={ref}")


PARSER_JUNK = list("()+-*/^,.= xyz0123456789IEeij[]{}<>!&|~%@:;'\"\\") + ["**", "sin(", "Symbol(", "Integer(", "oo", "nan", "pi", "lambda", "_", "1e", "Rational(", "Float(", "Add(", "Mul(", "Pow(", "x=", "\n", "\t", "\x00", "é"]


def gen_parser(rng):
    return gen_core(rng)


def chk_parser(t, rng):
    s, x = built(t)
    # 1. the plain-symbol version of sympy's str output parses to the same thing
    plain = {v: _sympy.Symbol(str(v)) for v in s.free_symbols}
    sp_plain = s.xreplace(plain)
    st = str(sp_plain)
    xp = sx_call(f"sympify({st!r})", _symox.sympify, st)
    compare_exprs(f"sympify({st!r})", sp_plain, xp, rng)
    # 2. random mutations of the string: symox must raise ValueError/TypeError or parse
    for _ in range(3):
        cs = list(st)
        for _ in range(rng.randint(1, 3)):
            op = rng.random()
            i = rng.randrange(len(cs) + 1)
            if op < 0.4 and cs:
                del cs[min(i, len(cs) - 1)]
            elif op < 0.8:
                cs.insert(i, rng.choice(PARSER_JUNK))
            else:
                cs[i:i] = cs[: rng.randint(0, len(cs))]
        junk = "".join(cs)
        stage("symox:parse junk " + repr(junk))
        try:
            _symox.sympify(junk)
        except (ValueError, TypeError, ZeroDivisionError, SyntaxError, OverflowError):
            pass
        except BaseException as e:  # noqa: BLE001
            if isinstance(e, (KeyboardInterrupt, SystemExit)):
                raise
            raise Finding("crash", f"sympify({junk!r}): {type(e).__name__}: {e}") from None


MODES = {
    "core": (gen_core, CORE_CHECKS),
    "integrate": (gen_integrand, {"integrate": chk_integrate}),
    "solve": (gen_solve, {"solve": chk_solve}),
    "series": (gen_series, {"series": chk_series}),
    "limit": (gen_limit, {"limit": chk_limit}),
    "summation": (gen_sum, {"summation": chk_sum}),
    "simplify": (gen_simplify, {"simplify": chk_simplify}),
    "lambdify": (gen_lambdify, {"lambdify": chk_lambdify}),
    "parser": (gen_parser, {"parser": chk_parser}),
}


def check_seed(seed, check):
    return (seed * 1000003 + sum(map(ord, check))) & 0xFFFFFFFF


def run_check(mode, check, t, seed):
    """Run one check; returns (cls, msg) or None. cls 'skip' for sympy failures."""
    init_libs()
    fn = MODES[mode][1][check]
    rng = random.Random(check_seed(seed, check))
    stage(f"{mode}:{check}:{to_src(t)}")
    try:
        fn(t, rng)
    except Finding as f:
        return (f.cls, f.msg)
    except SympySkip:
        return ("skip", "")
    except RecursionError:
        return ("skip", "recursion")
    except Exception as e:  # noqa: BLE001
        # an exception in the harness itself, typically sympy choking on a
        # symox result (e.g. while evaluating it)
        tb = traceback.format_exc(limit=3)
        if "symox" in tb and "sympy" not in tb.split("\n")[-2]:
            return ("crash", f"{check}: {type(e).__name__}: {str(e)[:300]}")
        return ("skip", f"harness: {type(e).__name__}: {e}")
    return None


def run_case(mode, seed, only=None):
    gen, checks = MODES[mode]
    t = gen(random.Random(seed))
    out = []
    for name in checks:
        if only and name not in only:
            continue
        r = run_check(mode, name, t, seed)
        if r is not None and r[0] != "skip":
            out.append((name, r[0], r[1]))
    return t, out


# ---------------------------------------------------------------------------
# worker processes
# ---------------------------------------------------------------------------


def worker_main(conn, state):
    global STATE
    STATE = state
    sys.setrecursionlimit(10000)
    init_libs()
    while True:
        try:
            task = conn.recv()
        except EOFError:
            return
        kind = task[0]
        if kind == "case":
            _, mode, seed = task
            t, res = run_case(mode, seed)
            conn.send(("ok", t, res))
        elif kind == "check":
            _, mode, check, t, seed = task
            conn.send(("ok", run_check(mode, check, t, seed)))
        elif kind == "quit":
            return


class Worker:
    def __init__(self, ctx, target=worker_main, stop=None):
        self.ctx = ctx
        self.target = target
        self.stop = stop
        self.proc = None
        self.restarts = 0
        self.start()

    def start(self):
        self.conn, child = self.ctx.Pipe()
        self.state = self.ctx.Array("c", 1024, lock=False)
        self.proc = self.ctx.Process(target=self.target, args=(child, self.state), daemon=True)
        self.proc.start()
        child.close()

    def kill(self):
        try:
            self.proc.kill()
            self.proc.join(5)
        except Exception:  # noqa: BLE001
            pass
        finally:
            self.conn.close()

    def call(self, task, timeout, stop=None):
        """Returns ('ok', ...), ('timeout', state) or ('died', exitcode, state)."""
        stop = self.stop if stop is None else stop
        try:
            self.conn.send(task)
        except (BrokenPipeError, OSError):
            self.kill()
            self.start()
            self.restarts += 1
            self.conn.send(task)
        deadline = time.monotonic() + timeout
        while True:
            if stop is not None and stop.is_set():
                raise FuzzStopped()
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            if self.conn.poll(min(0.05, remaining)):
                try:
                    return self.conn.recv()
                except (EOFError, ConnectionResetError, OSError):
                    break
        state = self.state.value.decode("utf8", "replace")
        alive = self.proc.is_alive()
        self.kill()
        code = self.proc.exitcode
        self.start()
        self.restarts += 1
        if alive:
            return ("timeout", state)
        return ("died", code, state)


class FuzzStopped(Exception):
    """Cooperative cancellation, including during a shrinking worker call."""


# ---------------------------------------------------------------------------
# shrinking
# ---------------------------------------------------------------------------


def subtrees(t, path=()):
    yield path, t
    if t[0] in ("+", "-", "*", "/", "**"):
        yield from subtrees(t[1], path + (1,))
        yield from subtrees(t[2], path + (2,))
    elif t[0] == "neg":
        yield from subtrees(t[1], path + (1,))
    elif t[0] == "f":
        for i in range(2, len(t)):
            yield from subtrees(t[i], path + (i,))
    elif t[0] == "rel":
        yield from subtrees(t[2], path + (2,))
        yield from subtrees(t[3], path + (3,))


def replace_at(t, path, new):
    if not path:
        return new
    i = path[0]
    lst = list(t)
    lst[i] = replace_at(t[i], path[1:], new)
    return tuple(lst)


def candidates(t):
    """One step simplifications of t, smallest first."""
    out = []
    for path, sub in subtrees(t):
        h = sub[0]
        reps = []
        if h in ("+", "-", "*", "/", "**"):
            reps += [sub[1], sub[2]]
        elif h == "neg":
            reps += [sub[1]]
        elif h == "f":
            reps += list(sub[2:])
            if sub[1] not in ("log2", "atan2", "Max", "Min"):
                pass
        elif h == "rel":
            reps += [sub[2], sub[3]]
        if h == "N" and abs(sub[1]) > 1:
            reps += [("N", 1), ("N", 2), ("N", sub[1] // 2)]
        if h == "N" and sub[1] < 0:
            reps += [("N", -sub[1])]
        if h == "Q":
            reps += [("N", 1), ("Q", 1, 2), ("N", sub[1])]
        if h == "F":
            reps += [("N", 1), ("F", 0.5), ("Q", 1, 2)]
        if h == "C":
            reps += [("N", 1), ("S", "x")]
        if h == "S" and sub[1] != "x":
            reps += [("S", "x"), ("N", 2)]
        if h not in ("N", "S") and path:
            reps += [("S", "x"), ("N", 1), ("N", 2)]
        if h == "f" and len(sub) == 3:
            # same argument, simpler function
            for g in ("exp", "sin", "log"):
                if sub[1] != g:
                    reps.append(("f", g, sub[2]))
        for r in reps:
            if r != sub:
                c = replace_at(t, path, r)
                if c[0] != "rel" or t[0] == "rel":
                    out.append(c)
    out.sort(key=tree_size)
    return out


def same_failure(a, b):
    if a is None or b is None:
        return False
    if a[0] != b[0]:
        return False
    if a[0] in ("crash",):
        # same exception type and message head
        ka = a[1].split(":")[0:3]
        kb = b[1].split(":")[0:3]
        return ka[:2] == kb[:2]
    if a[0] == "died":
        return True
    if a[0] in ("wrong", "form", "weak"):
        # same operation
        return a[1].split(":")[0].split("(")[0] == b[1].split(":")[0].split("(")[0]
    return True


def eval_candidate(worker, mode, check, t, seed, timeout):
    r = worker.call(("check", mode, check, t, seed), timeout)
    if r[0] == "ok":
        return r[1]
    if r[0] == "timeout":
        return ("timeout", r[1])
    return ("died", f"exit {r[1]}: {r[2]}")


def shrink(worker, mode, check, t, seed, failure, timeout, budget=400):
    best = t
    used = 0
    seen = {t}
    progress = True
    while progress and used < budget:
        progress = False
        for c in candidates(best):
            if c in seen or tree_size(c) > tree_size(best):
                continue
            seen.add(c)
            used += 1
            r = eval_candidate(worker, mode, check, c, seed, timeout)
            if r is not None and r[0] == "timeout" and failure[0] != "timeout":
                continue
            if same_failure(failure, r):
                best, failure = c, r
                progress = True
                break
            if used >= budget:
                break
    return best, failure


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def signature(mode, check, failure, t):
    """Deduplication key: operation + class + minimal source with numbers abstracted."""
    import re
    src = to_src(t)
    norm = re.sub(r"Integer\(-?\d+\)", "N", src)
    norm = re.sub(r"Rational\(-?\d+, \d+\)", "Q", norm)
    norm = re.sub(r"Float\([^)]*\)", "F", norm)
    head = failure[1].split(":")[0] if failure[0] in ("wrong", "form", "weak") else failure[1][:60]
    if failure[0] == "crash":
        head = failure[1].split(":")[0] + ":" + failure[1].split(":")[1] if ":" in failure[1] else failure[1]
    return f"{mode}|{check}|{failure[0]}|{head}|{norm}"


def drive(args):
    ctx = mp.get_context("spawn")
    os.makedirs(args.out, exist_ok=True)
    modes = list(MODES) if args.mode == "all" else args.mode.split(",")
    lo, hi = (int(v) for v in args.seeds.split(":"))
    sigs = set()
    fpath = os.path.join(args.out, "findings.jsonl")
    if os.path.exists(fpath):
        with open(fpath) as fh:
            for line in fh:
                try:
                    sigs.add(json.loads(line)["sig"])
                except (ValueError, KeyError):
                    pass
    stats = {"cases": 0, "findings": 0, "new": 0, "by_class": {}, "timeouts": 0,
             "deaths": 0, "skipped_sympy_timeouts": 0, "interrupted": False}
    jobs = min(args.jobs, (hi - lo) * len(modes))
    deadline = time.monotonic() + args.time if args.time else None
    stop = ctx.Event()
    task_q, res_q = ctx.Queue(), ctx.Queue()
    procs = []
    seeds = iter((m, s) for s in range(lo, hi) for m in modes)
    inflight = 0
    last_report = time.monotonic()
    try:
        with open(fpath, "a") as fh:
            for i in range(jobs):
                proc = ctx.Process(target=driver_slave, args=(i, task_q, res_q, args.timeout,
                                                            args.shrink_budget, stop))
                proc.start()
                procs.append(proc)
                task_q.put(next(seeds))
                inflight += 1
            while inflight:
                remaining = deadline - time.monotonic() if deadline else 1.0
                if remaining <= 0:
                    stats["interrupted"] = True
                    break
                try:
                    msg = res_q.get(timeout=min(0.25, remaining))
                except queue.Empty:
                    dead = [p.exitcode for p in procs if not p.is_alive()]
                    if dead:
                        raise RuntimeError(f"fuzz driver exited unexpectedly: {dead}")
                    continue
                if msg[0] == "done":
                    inflight -= 1
                    stats["cases"] += 1
                    for k in ("timeouts", "deaths", "skipped_sympy_timeouts"):
                        stats[k] += msg[2].get(k, 0)
                    nxt = next(seeds, None)
                    if nxt is not None and not (deadline and time.monotonic() >= deadline):
                        task_q.put(nxt)
                        inflight += 1
                elif msg[0] == "driver_error":
                    raise RuntimeError(f"fuzz driver {msg[1]} failed:\n{msg[2]}")
                elif msg[0] == "finding":
                    rec = msg[1]
                    stats["findings"] += 1
                    stats["by_class"][rec["cls"]] = stats["by_class"].get(rec["cls"], 0) + 1
                    if rec["sig"] not in sigs:
                        sigs.add(rec["sig"])
                        stats["new"] += 1
                        fh.write(json.dumps(rec) + "\n")
                        fh.flush()
                        if args.verbose:
                            print(f"[{rec['cls']}] {rec['mode']}/{rec['check']} {rec['min_src']}\n"
                                  f"    {rec['min_msg'][:400]}", flush=True)
                if time.monotonic() - last_report > 30:
                    last_report = time.monotonic()
                    print(f"stats: {json.dumps(stats)}", flush=True)
    except KeyboardInterrupt:
        stats["interrupted"] = True
    finally:
        # Cancellation is checked inside every worker call, including shrinking.
        # Drivers close and reap their own native workers before exiting.
        stop.set()
        for _ in procs:
            task_q.put(None)
        for p in procs:
            p.join(2)
            if p.is_alive():
                p.terminate()
                p.join(2)
        task_q.close()
        task_q.cancel_join_thread()
        res_q.close()
        with open(os.path.join(args.out, "stats.json"), "w") as sf:
            json.dump(stats, sf)
    print(f"final stats: {json.dumps(stats)}", flush=True)
    return stats


def driver_slave(idx, task_q, res_q, timeout, budget, stop):
    """Feeds cases to its own worker, shrinks failures, and always reaps it."""
    ctx = mp.get_context("spawn")
    w = Worker(ctx, stop=stop)
    try:
        while True:
            if stop.is_set():
                return
            task = task_q.get()
            if task is None:
                return
            mode, seed = task
            st = {}
            r = w.call(("case", mode, seed), timeout)
            fails = []
            if r[0] == "ok":
                t, res = r[1], r[2]
                fails = [(c, (cls, msg)) for c, cls, msg in res]
            else:
                state = r[1] if r[0] == "timeout" else r[2]
                parts = state.split(":", 2)
                t = MODES[mode][0](random.Random(seed))
                if state.startswith("sympy"):
                    st["skipped_sympy_timeouts"] = 1
                else:
                    # which check? re-run each check separately to find out
                    st["timeouts" if r[0] == "timeout" else "deaths"] = 1
                    for c in MODES[mode][1]:
                        rr = eval_candidate(w, mode, c, t, seed, timeout)
                        if rr is not None and rr[0] in ("timeout", "died", "crash", "wrong", "form", "weak", "slow"):
                            if rr[0] == "timeout" and "sympy" in rr[1].split(":")[0]:
                                continue
                            fails.append((c, rr))
                del parts
            for check, failure in fails:
                mt, mf = shrink(w, mode, check, t, seed, failure, timeout, budget)
                rec = {
                    "mode": mode, "seed": seed, "check": check, "cls": failure[0],
                    "src": to_src(t), "msg": failure[1][:2000],
                    "min_src": to_src(mt), "min_msg": (mf[1] if mf else failure[1])[:2000],
                    "min_tree": repr(mt),
                }
                rec["sig"] = signature(mode, check, mf or failure, mt)
                res_q.put(("finding", rec))
            res_q.put(("done", seed, st))
    except FuzzStopped:
        pass
    except BaseException:
        res_q.put(("driver_error", idx, traceback.format_exc()))
    finally:
        w.kill()


def replay(mode, seed, check=None, timeout=20):
    # Replay is isolated too: reproducing a crash must not kill the frontend.
    t = MODES[mode][0](random.Random(seed))
    worker = Worker(mp.get_context("spawn"))
    try:
        if check:
            failure = eval_candidate(worker, mode, check, t, seed, timeout)
            res = [(check, *failure)] if failure and failure[0] != "skip" else []
        else:
            result = worker.call(("case", mode, seed), timeout)
            if result[0] == "ok":
                t, res = result[1:]
            else:
                res = [("case", result[0], str(result[1:]))]
    finally:
        worker.kill()
    print("tree:", to_src(t))
    for c, cls, msg in res:
        print(f"  [{cls}] {c}: {msg}")
    if not res:
        print("  no findings")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", default="core", help="core, integrate, solve, series, limit, summation, simplify, lambdify, parser, all, or a comma list")
    ap.add_argument("--seeds", default="0:100000000")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--time", type=float, default=0, help="stop after this many seconds")
    ap.add_argument("--timeout", type=float, default=20.0, help="per case timeout (seconds)")
    ap.add_argument("--shrink-budget", type=int, default=300)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"))
    ap.add_argument("--replay", nargs="+", metavar="ITEM", help="MODE SEED [CHECK], isolated with --timeout")
    ap.add_argument("--fail-on", default="", help="comma-separated finding classes that give exit status 1")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    if args.jobs < 1:
        ap.error("--jobs must be at least 1")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        ap.error("--timeout must be finite and positive")
    if not math.isfinite(args.time) or args.time < 0:
        ap.error("--time must be finite and nonnegative")
    if args.shrink_budget < 0:
        ap.error("--shrink-budget must be nonnegative")
    classes = {"crash", "died", "timeout", "slow", "wrong", "form", "weak"}
    fail_on = set(args.fail_on.split(",")) if args.fail_on else set()
    if not fail_on <= classes:
        ap.error("--fail-on contains an unknown finding class")
    if args.replay:
        if not 2 <= len(args.replay) <= 3 or args.replay[0] not in MODES:
            ap.error("--replay expects MODE SEED [CHECK] with a supported mode")
        try:
            seed = int(args.replay[1])
        except ValueError:
            ap.error("replay SEED must be an integer")
        check = args.replay[2] if len(args.replay) == 3 else None
        if check and check not in MODES[args.replay[0]][1]:
            ap.error("unknown check for replay mode")
        return int(any(cls in fail_on for _, cls, _ in replay(args.replay[0], seed, check, args.timeout)))
    modes = list(MODES) if args.mode == "all" else args.mode.split(",")
    if any(mode not in MODES for mode in modes):
        ap.error("--mode contains an unknown mode")
    try:
        lo, hi = (int(v) for v in args.seeds.split(":"))
        if lo < 0 or hi < lo:
            raise ValueError
    except ValueError:
        ap.error("--seeds must be a nonnegative START:STOP range")
    stats = drive(args)
    return int(any(stats["by_class"].get(cls, 0) for cls in fail_on))


if __name__ == "__main__":
    sys.exit(main())
