"""Differential tests of the printers, cse and code generation against SymPy.

Every case is a Python expression evaluated in a SymPy namespace and in a
symox namespace; the outputs of ``latex``, ``pretty``, ``ccode``,
``cxxcode``, ``pycode``, the NumPy printer and ``rust_code`` are compared
string by string. When the two canonical forms differ already (``str``
differs), a printer mismatch is reported as ``str-diff`` (not a printer
bug).

``cse`` is checked semantically: substituting the replacements back must
give the original expression and the number of replacements must be
comparable with SymPy's.

Generated C functions are compiled with ``cc`` (when available) and
evaluated through ``ctypes``; the CUDA kernels are compiled as C++ with a
small host shim (``__global__``, ``blockIdx``, intrinsics) and executed on
the CPU.

Run directly for a report, or through pytest.
"""

from __future__ import annotations

import ctypes
import os
import random
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import sympy
from sympy.printing.numpy import NumPyPrinter

import symox

CASES = r"""
x
x + y
x - y
-x - y
x*y
-x
2*x
-2*x
x/2
-x/2
x/y
x*y/z
x*y/(2*z)
1/x
2/x
-2/x
1/(2*x)
x**2
x**3
x**-1
x**-2
x**y
x**(-y)
x**(y + 1)
x**(2*y)
(x + 1)**2
(x + 1)/(x - 1)
1/(x + y)
-1/(x + 1)
(x*y)**-1
x*(x + y)**-2
x**2*y**-3*z
sqrt(x)
sqrt(x + 1)
1/sqrt(x)
x**Rational(1, 3)
x**Rational(-1, 3)
x**Rational(2, 3)
x**Rational(3, 2)
x**Rational(-3, 2)
sqrt(2)
sqrt(2)*x
sqrt(2)/2
2**x
3**x
(-2)**x
Rational(1, 2)**x
Rational(1, 2)
Rational(-3, 4)
Rational(3, 4)*x**2
Rational(1, 3)
x + Rational(1, 2)
2*3**x
2*x*3**x
3*x**2*y
x**2 - 2*x + 1
x**2 + 1/x
x**3 + 3*x**2*y + 3*x*y**2 + y**3
pi
pi*x
pi/2
pi*x/2
2*pi*x
E
E*x
I*x
oo
-oo
nan
zoo
EulerGamma
x + 0.1
1.5*x
2.5*x**2 + 0.25
1e-20*x
123456789.123*x
x**0.5
x**2.5
exp(x)
exp(-x)
exp(x + y)
exp(x)*sin(x)
log(x)
log(2)
log(x)**2
sin(x)
cos(x)
tan(x)
sin(x) + cos(y)
sin(x)**2
sin(x)*cos(x)
x*sin(x)
sin(x)/x
sin(x + y)
sin(2*x)
sin(x/y)
sqrt(sin(x))
sin(x)/cos(y)**2
sec(x)
cot(x)
csc(x)
asin(x)
acos(x)
atan(x)
acot(x)
atan2(y, x)
sinh(x)
cosh(x)
tanh(x)
coth(x)
asinh(x)
acosh(x)
atanh(x)
asin(x)**2
Abs(x)
Abs(x + y)
Abs(x)**3
sign(x)
floor(x)
ceiling(x)
floor(x)**2
gamma(x)
loggamma(x)
erf(x)
erfc(x)
Heaviside(x)
DiracDelta(x)
re(x)
im(x)
conjugate(x)
Max(x, y)
Max(x, y, z)
Min(x, 1)
f(x)
f(x, y)
f(x)**2
Derivative(f(x), x)
Derivative(f(x), x, 2)
Derivative(f(x, y), x, y)
Eq(x, y)
Ne(x, y)
x < y
x <= y
x > 1
x >= y
Piecewise((x, x > 0), (0, True))
Piecewise((x, x > 0), (x**2, x < -1), (1, True))
alpha + beta
alpha*Gamma
x1 + y2
theta_i
xdot
alphahat
x**2*exp(x)*sin(y)/(1 + z)
(x + y)**2*sin(x + y) + cos(x + y)/(x + y)
-x**2/(2*y) + 3*z
x*y*z + x**2*y + y**2 + 1
""".strip().splitlines()


def namespace(mod):
    ns = {name: getattr(mod, name) for name in dir(mod) if not name.startswith("_")}
    names = "x y z alpha beta Gamma x1 y2 theta_i xdot alphahat"
    for s in mod.symbols(names):
        ns[s.name] = s
    ns["f"] = mod.Function("f")
    return ns


NS_S = namespace(sympy)
NS_X = namespace(symox)


def _sympy_numpy(e):
    return NumPyPrinter().doprint(e)


PRINTERS = {
    "latex": (sympy.latex, symox.latex),
    "pretty": (lambda e: sympy.pretty(e, use_unicode=True), symox.pretty),
    "ccode": (sympy.ccode, symox.ccode),
    "cxxcode": (sympy.cxxcode, symox.cxxcode),
    "pycode": (sympy.pycode, symox.pycode),
    "numpy": (_sympy_numpy, symox.numpycode),
    "rust_code": (sympy.rust_code, symox.rust_code),
}

# Deliberate deviations (documented in the module docs / report).
KNOWN_DEVIATIONS = {
    # SymPy's rust printer rewrites Max/Min into an if-else chain that mixes
    # integer and float literals; symox emits `x.max(y)` method chains.
    ("rust_code", "Max(x, y)"),
    ("rust_code", "Max(x, y, z)"),
    ("rust_code", "Min(x, 1)"),
    # SymPy's rust printer drops the parentheses around TypeCast'ed sums
    # ("PI*x + y" for pi*(x + y)); symox keeps them (valid code).
    # NumPy conjugate is supported by symox (numpy.conjugate).
    ("numpy", "conjugate(x)"),
    # SymPy raises (TypeError / KeyError) where symox prints valid code.
    ("rust_code", "Rational(1, 2)**x"),
    ("numpy", "zoo"),
}


def run_printer(fn, e):
    try:
        return "ok", fn(e)
    except Exception as ex:  # noqa: BLE001
        return "err", type(ex).__name__


def compare_printers(verbose=True):
    stats = {}
    bad = []
    for case in CASES:
        es = eval(case, NS_S)
        try:
            ex = eval(case, NS_X)
        except Exception as err:  # noqa: BLE001
            stats["eval-error"] = stats.get("eval-error", 0) + 1
            if verbose:
                print(f"[eval-error ] {case}: {err}")
            bad.append(("eval", case, str(err), ""))
            continue
        str_same = str(es) == str(ex)
        for pname, (fs, fx) in PRINTERS.items():
            ks, vs = run_printer(fs, es)
            kx, vx = run_printer(fx, ex)
            if ks == "err" and kx == "err":
                status = "both-error"
            elif ks == kx and vs == vx:
                status = "same"
            elif (pname, case) in KNOWN_DEVIATIONS:
                status = "deviation"
            elif not str_same:
                status = "str-diff"
            else:
                status = "DIFFERENT"
            stats[(pname, status)] = stats.get((pname, status), 0) + 1
            if status in ("DIFFERENT", "str-diff") and verbose:
                print(f"[{pname:9s} {status:9s}] {case}\n    symox: {vx!r}\n    sympy: {vs!r}")
            if status == "DIFFERENT":
                bad.append((pname, case, vx, vs))
    if verbose:
        summary = {}
        for (p, st), n in sorted(stats.items(), key=str):
            summary.setdefault(p, {})[st] = n
        print("\nprinter summary:")
        for p, d in summary.items():
            print(f"  {p:10s} {d}")
    return bad


# ----------------------------------------------------------------------
# cse
# ----------------------------------------------------------------------

CSE_CASES = r"""
(x + y)**2 + sin(x + y)
[x + y + z, sin(x + y + w)]
x**-2 + x**2*sin(y)
exp(-x*y) + cos(x*y)
[sin(x)*cos(y) + cos(x)*sin(y), sin(x)*cos(y) - cos(x)*sin(y)]
(x*y*z + 1)/(x*y*z - 1) + sqrt(x*y*z)
[(x + y)*(x + z), (x + y)*(y + z), (x + z)*(y + z)]
diff(sin(x*y)*exp(x**2 + y), x, 3)
[diff(exp(x*y)/(1 + x**2 + y**2), x), diff(exp(x*y)/(1 + x**2 + y**2), y)]
expand((x + y + z + 1)**4) + expand((x + y + z)**3)
""".strip().splitlines()


def _back(reps, red):
    out = []
    for r in red:
        for s, v in reversed(reps):
            r = r.xreplace({s: v})
        out.append(r)
    return out


def compare_cse(verbose=True):
    bad = []
    nss = dict(NS_S, w=sympy.Symbol("w"))
    nsx = dict(NS_X, w=symox.Symbol("w"))
    for case in CSE_CASES:
        ex = eval(case, nsx)
        ex = list(ex) if isinstance(ex, list) else [ex]
        # run SymPy's cse on the very same expressions (symox's canonical
        # forms, e.g. of derivatives, may differ structurally)
        es = [sympy.sympify(str(e), locals=nss) for e in ex]
        rs, reds = sympy.cse(es)
        rx, redx = symox.cse(ex)
        back = _back(rx, redx)
        ok_back = [str(b) == str(e) or sympy.simplify(sympy.sympify(str(b)) - sympy.sympify(str(e))) == 0
                   for b, e in zip(back, ex)]
        ok_count = abs(len(rx) - len(rs)) <= max(1, len(rs) // 4)
        status = "ok" if all(ok_back) and ok_count else "BAD"
        if verbose:
            print(f"[cse {status}] {case[:60]:60s} symox: {len(rx):3d} reps   sympy: {len(rs):3d} reps")
        if status != "ok":
            bad.append((case, len(rx), len(rs), ok_back))
    return bad


# ----------------------------------------------------------------------
# compiling generated C and CUDA code
# ----------------------------------------------------------------------

CC = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
CXX = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")

CODEGEN_CASES = [
    "sin(x)*cos(y) + x**2/y",
    "exp(-x**2 - y**2)*(x + y)**3 + sqrt(x**2 + y**2)",
    "atan2(y, x) + log(1 + x**2) + Abs(x - y) + floor(x)",
    "Piecewise((x**2, x > y), (y - x, True)) + Max(x, y)",
    "1/(1 + exp(-x*y)) + tanh(x) + asinh(y) + erf(x) + gamma(1 + x**2)",
    "sec(x) + cot(y) + coth(x + 1) + sign(y) + pi*E*x + EulerGamma*y",
    "x**Rational(1, 3) + x**Rational(2, 3) + 1/sqrt(1 + y**2) + x**-2",
    "cosh(x)**4 - sinh(y)**3 + x**5*y**-3",
]

CUDA_SHIM = r"""
#include <cmath>
#include <cstdio>
using namespace std;
struct dim3_ { long x; };
static dim3_ blockIdx = {0}, threadIdx = {0}, blockDim = {1}, gridDim = {1};
#define __global__ extern "C"
#define __device__
#define __forceinline__ inline
#define __restrict__
static inline double rsqrt(double x) { return 1.0 / sqrt(x); }
static inline float rsqrtf(float x) { return 1.0f / sqrtf(x); }
#define __sinf(x) sinf(x)
#define __cosf(x) cosf(x)
#define __tanf(x) tanf(x)
#define __expf(x) expf(x)
#define __logf(x) logf(x)
#define __powf(x, y) powf(x, y)
"""


def _random_inputs(n, rnd):
    x = np.array([rnd.uniform(0.3, 2.0) for _ in range(n)])
    y = np.array([rnd.uniform(0.2, 1.7) for _ in range(n)])
    return x, y


def _reference(case, x, y):
    xs, ys = sympy.symbols("x y")
    e = eval(case, NS_S)
    f = sympy.lambdify((xs, ys), e, modules=["mpmath"])
    return np.array([float(f(a, b)) for a, b in zip(x, y)])


def _compile_and_run(src, lang, nouts, x, y, dtype=np.float64, compiler=None):
    tmp = tempfile.mkdtemp(prefix="symox_codegen_")
    ext = "c" if lang == "c" else "cpp"
    path = os.path.join(tmp, f"gen.{ext}")
    lib = os.path.join(tmp, "gen.so")
    with open(path, "w") as fh:
        fh.write(src)
    cmd = [compiler, "-O2", "-shared", "-fPIC", "-o", lib, path, "-lm"]
    if lang == "c":
        cmd.insert(1, "-std=c99")
        cmd.insert(2, "-D_DEFAULT_SOURCE")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"compilation failed:\n{src}\n{proc.stderr}")
    so = ctypes.CDLL(lib)
    n = len(x)
    xa = np.ascontiguousarray(x, dtype=dtype)
    ya = np.ascontiguousarray(y, dtype=dtype)
    outs = [np.zeros(n, dtype=dtype) for _ in range(nouts)]
    ptr = lambda a: a.ctypes.data_as(ctypes.c_void_p)  # noqa: E731
    so.f.restype = None
    so.f(ctypes.c_long(n), ptr(xa), ptr(ya), *[ptr(o) for o in outs])
    shutil.rmtree(tmp, ignore_errors=True)
    return outs


def compare_codegen(verbose=True):
    if CC is None:
        if verbose:
            print("no C compiler found: skipping compile tests")
        return []
    bad = []
    rnd = random.Random(7)
    x, y = _random_inputs(64, rnd)
    xs, ys = symox.symbols("x y")
    exprs = [eval(c, NS_X) for c in CODEGEN_CASES]
    refs = [_reference(c, x, y) for c in CODEGEN_CASES]
    for use_cse in (True, False):
        src = symox.ccode_function("f", [xs, ys], exprs, cse=use_cse)
        outs = _compile_and_run(src, "c", len(exprs), x, y, compiler=CC)
        for case, o, r in zip(CODEGEN_CASES, outs, refs):
            ok = np.allclose(o, r, rtol=1e-10, atol=1e-12)
            if not ok:
                bad.append(("C", use_cse, case))
            if verbose:
                print(f"[C cse={use_cse!s:5s} {'ok' if ok else 'BAD'}] {case}")
    if CXX is not None:
        for dtype, fast, tol in (("double", False, 1e-10), ("float", False, 2e-4), ("float", True, 2e-3)):
            src = CUDA_SHIM + symox.cuda_kernel("f", [xs, ys], exprs, dtype=dtype, fast_math=fast)
            npd = np.float64 if dtype == "double" else np.float32
            outs = _compile_and_run(src, "cuda", len(exprs), x, y, dtype=npd, compiler=CXX)
            for case, o, r in zip(CODEGEN_CASES, outs, refs):
                ok = np.allclose(o.astype(np.float64), r, rtol=tol, atol=tol)
                if not ok:
                    bad.append(("CUDA", dtype, fast, case))
                if verbose:
                    print(f"[CUDA {dtype:6s} fast={fast!s:5s} {'ok' if ok else 'BAD'}] {case}")
    return bad


# ----------------------------------------------------------------------
# pytest entry points
# ----------------------------------------------------------------------


def test_printers_against_sympy():
    bad = compare_printers(verbose=False)
    assert not bad, "\n".join(f"{p} {c}: symox={x!r} sympy={s!r}" for p, c, x, s in bad)


def test_cse_semantics():
    bad = compare_cse(verbose=False)
    assert not bad, bad


def test_generated_code_compiles_and_runs():
    bad = compare_codegen(verbose=False)
    assert not bad, bad


def test_repr_latex():
    x = symox.Symbol("x")
    assert (x**2 / 2)._repr_latex_() == sympy.Symbol("x")._repr_latex_().replace("x", r"\frac{x^{2}}{2}")


if __name__ == "__main__":
    b1 = compare_printers()
    print()
    b2 = compare_cse()
    print()
    b3 = compare_codegen()
    print(f"\nprinter mismatches: {len(b1)}, cse failures: {len(b2)}, codegen failures: {len(b3)}")
    sys.exit(1 if (b1 or b2 or b3) else 0)
