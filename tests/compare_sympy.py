"""Differential testing against SymPy.

Every case is a Python expression evaluated in a SymPy namespace and in a
symox namespace. We check

* ``str`` equality (the canonical forms should match SymPy's), and
* semantic equality (numerically, at random complex points).

Run directly for a report, or through pytest (``test_against_sympy``).
"""

from __future__ import annotations

import cmath
import random
import sys

import sympy
import symox

CASES = r"""
x + y
x + x
2*x + 3*x
x - x
x*x
x**2*x**3
x*(1/x)
(x+1)/(x+1)
x/y/x
-(x+y)
2*(x+y)
x*(x+y)
Rational(1,2)*(x+1)
(1-x)/2
x/2
-x/2
1/x
2/x
x*y/z
x**-2
-x**2
x**Rational(1,2)
x**Rational(-1,2)
1/(x+1)
(x+1)**2
(-x)**2
(-x)**3
(x*y)**2
(x**y)**2
(x**2)**Rational(1,2)
(x**Rational(1,2))**2
x**y*x**z
x*x**y
x**2*x**y
x**y*x**(-y)
x**y*x**y
x**(2*y)*x**y
2**x*2**y
2**x*2
x**(y+1)*x
Rational(3,2)*x*y
12**Rational(1,3)
8**Rational(1,2)
4**Rational(1,4)
Rational(1,2)**Rational(1,2)
(-8)**Rational(1,3)
sqrt(2)*sqrt(3)
sqrt(2)*sqrt(6)
2**Rational(1,3)*3**Rational(1,3)
2**Rational(1,2)*3**Rational(1,3)
sqrt(-4)
sqrt(2*x)
sqrt(-2*x)
(2*x)**Rational(1,3)
(Rational(1,2)*x)**Rational(1,2)
Rational(4,9)**Rational(1,2)
27**Rational(2,3)
1024**Rational(1,10)
(-1)**Rational(1,2)
(-1)**Rational(3,2)
(-1)**Rational(5,3)
(-1)**Rational(2,3)
(-2)**Rational(3,2)
sqrt(8)/2
I**3
I*I
(1+I)**2
I*x + I*y
exp(x)*exp(y)
exp(x)*exp(-x)
exp(x)*exp(y)*exp(-y)
exp(x)**2
E**x
E*E
E*exp(x)
exp(1)
exp(0)
exp(log(x))
log(exp(x))
exp(2*log(x))
exp(I*pi)
exp(I*pi/2)
log(1)
log(E)
log(-1)
log(-2)
log(Rational(1,2))
sin(-x)
cos(-x)
tan(-x)
sin(0)
sin(pi)
sin(pi/6)
sin(pi/4)
sin(pi/3)
sin(pi/2)
sin(2*pi/3)
sin(5*pi/6)
sin(7*pi/6)
sin(pi/12)
sin(5*pi/12)
cos(pi/4)
cos(pi/3)
cos(2*pi/3)
cos(pi)
cos(3*pi/2)
tan(pi/4)
tan(pi/3)
tan(pi/6)
tan(pi/12)
sin(x + pi)
sin(x + 2*pi)
sin(x + pi/2)
cos(x + pi)
cos(x + pi/2)
sin(I*x)
cos(I*x)
sin(asin(x))
cos(acos(x))
tan(atan(x))
sin(acos(x))
cos(asin(x))
asin(1)
asin(Rational(1,2))
acos(0)
acos(-1)
atan(1)
atan(sqrt(3))
atan(oo)
sinh(0)
cosh(0)
sinh(-x)
cosh(-x)
Abs(-3)
Abs(-x)
Abs(2*x)
Abs(I)
sign(-5)
floor(Rational(7,2))
ceiling(Rational(7,2))
floor(-Rational(7,2))
gamma(5)
gamma(Rational(1,2))
gamma(Rational(3,2))
x+oo
2*oo
oo-oo
0*oo
oo+1
x*oo
-x*oo
2*x*oo
1/S(0)
2.0*x
x + 1.5
x + 0.0
0.0*x
2.0*(x+1)
x**2.0
sin(0.5)
exp(1.0)
sqrt(2.0)
S(2.0)
S(0.5)
S(1e-5)
S(1e20)
2.5*x
-2.5*x
2.5*I
1.0 + 2.0*I
x/(2*y)
2*x/(3*y**2)
-x/(y*z)
x**2/(x+1)
(x+1)/(x-1)
sqrt(x)/x
x/sqrt(y)
1/(2*x)
-1/x
I/2
-I*x/3
x*y**-1*z**-2
exp(-x)
2*exp(-x)/x
x**y**z
(x**y)**z
x + 2*y + 3 + x**2 + y**3
x*y*z + x**2*y + y**2 + 1
-1 - x
(x+y)*(x-y)
sin(x)**2 + cos(x)**2
2**x
pi*x
expand((x+y)**2)
expand((x+y)**3)
expand((x+y+z)**3)
expand((x+1)*(x-1))
expand((x+1)**2*(x-1))
expand(x*(x+y))
expand((x+y)*(x+z)*(y+z))
expand((x+1)**-2)
expand(2*(x+1)**2)
expand((sqrt(2)+1)**2)
expand((sqrt(2)+1)*(sqrt(2)-1))
expand((x+I)**2)
expand((1+I)**3)
expand(exp(x+y))
expand(x**(y+2))
expand(2**(x+1))
expand(pi**(x+y))
expand((x+1)**Rational(5,2))
expand((x+sin(x))**2)
expand((x+1/x)**2)
expand((2*x+3*y)**4)
expand((x*y+1)**2)
diff(x**2, x)
diff(x**3*y, x)
diff(sin(x), x)
diff(cos(x), x)
diff(tan(x), x)
diff(exp(x**2), x)
diff(log(x), x)
diff(sin(x)*cos(x), x)
diff(x**x, x)
diff(2**x, x)
diff(sqrt(x), x)
diff(1/x, x)
diff(atan(x), x)
diff(asin(x), x)
diff(acos(x), x)
diff(sinh(x), x)
diff(cosh(x), x)
diff(tanh(x), x)
diff(exp(sin(x)), x)
diff(x**2*sin(x), x, 2)
diff(x**4, x, 3)
diff(x*y, x, y)
diff(f(x), x)
diff(f(x)*x, x)
diff(f(x), x, 2)
diff(sin(f(x)), x)
diff(erf(x), x)
diff(Abs(x), x)
(x**2+2*x+1).subs(x, 1)
(x**2+2*x+1).subs(x, y)
(x+y).subs({x: 1, y: 2})
(x*y).subs(x, y)
(x + y + z).subs(x + y, 1)
(x**4).subs(x**2, y)
(sin(x)+cos(x)).subs(x, pi)
(x*y*z).subs(x*y, 2)
(exp(x)).subs(x, log(2))
Eq(x, x)
Eq(1, 2)
Lt(1, 2)
Lt(x, 1)
Piecewise((x, x > 0), (0, True))
Piecewise((1, S(1) > 0), (0, True))
Max(1, 2, x)
Min(1, 2)
Max(x, x)
atan2(1, 1)
atan2(1, -1)
atan2(0, -1)
(x+y+z)**0
0**x
1**x
x**0
x**1
(-oo)**2
""".strip().splitlines()


# Deliberate deviations from SymPy (documented in README).
KNOWN_DEVIATIONS = {
    # SymPy returns the complex-aware formula; symox uses the real one.
    "diff(Abs(x), x)",
}


def namespace(mod):
    ns = {name: getattr(mod, name) for name in dir(mod) if not name.startswith("_")}
    x, y, z = mod.symbols("x y z")
    ns.update(x=x, y=y, z=z, f=mod.Function("f"))
    return ns


def numeric_equal(a, b, syms, trials=4):
    """Compare two SymPy expressions numerically at random points."""
    rnd = random.Random(1234)
    a_f = sympy.lambdify(syms, a, modules=["mpmath"])
    b_f = sympy.lambdify(syms, b, modules=["mpmath"])
    for _ in range(trials):
        pt = [complex(rnd.uniform(0.2, 1.5), rnd.uniform(-0.3, 0.3)) for _ in syms]
        try:
            va = complex(a_f(*pt))
            vb = complex(b_f(*pt))
        except (ZeroDivisionError, OverflowError, ValueError, TypeError):
            continue
        if cmath.isnan(va) and cmath.isnan(vb):
            continue
        if abs(va - vb) > 1e-8 * max(1.0, abs(va), abs(vb)):
            return False
    return True


def compare(case):
    ns_s = namespace(sympy)
    ns_x = namespace(symox)
    exp_s = eval(case, ns_s)
    try:
        exp_x = eval(case, ns_x)
    except Exception as e:  # noqa: BLE001
        return "error", f"{type(e).__name__}: {e}", str(exp_s)
    s_s, s_x = str(exp_s), str(exp_x)
    if s_s == s_x:
        return "same", s_x, s_s
    # semantic comparison
    try:
        back = sympy.sympify(s_x, locals={"f": sympy.Function("f")})
    except Exception as e:  # noqa: BLE001
        return "unparseable", s_x, s_s
    if isinstance(exp_s, sympy.Basic) and isinstance(back, sympy.Basic):
        if back == exp_s:
            return "same-tree", s_x, s_s
        syms = sorted(exp_s.free_symbols | back.free_symbols, key=str)
        if not getattr(exp_s, "is_Relational", False) and not exp_s.has(sympy.Derivative):
            try:
                if numeric_equal(back, exp_s, syms):
                    return "equivalent", s_x, s_s
            except Exception:  # noqa: BLE001
                pass
    return "DIFFERENT", s_x, s_s


def run(verbose=True):
    stats = {}
    bad = []
    for case in CASES:
        status, got, want = compare(case)
        stats[status] = stats.get(status, 0) + 1
        if status not in ("same", "same-tree"):
            if verbose:
                print(f"[{status:11s}] {case:40s} symox: {got:35s} sympy: {want}")
            if status in ("DIFFERENT", "error", "unparseable") and case not in KNOWN_DEVIATIONS:
                bad.append((case, got, want))
    if verbose:
        print("\nsummary:", stats, f"(total {len(CASES)})")
    return bad


def test_against_sympy():
    bad = run(verbose=False)
    assert not bad, "\n".join(f"{c}: symox={g} sympy={w}" for c, g, w in bad)


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
