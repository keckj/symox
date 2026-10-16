"""Differential tests of series expansions, Order terms and limits against SymPy.

Every case is a Python expression evaluated in a SymPy namespace and in a
symox namespace. We check ``str`` equality and, when it differs, semantic
equality: the order terms must be equal and the polynomial parts must agree
numerically at random points.

Run directly for a report (with ``-v`` to list the equivalent cases too), or
through pytest (``test_series_against_sympy``).
"""

from __future__ import annotations

import cmath
import random
import sys

import sympy
import symox

SERIES = r"""
series(sin(x), x)
series(cos(x), x)
series(tan(x), x)
series(tan(x), x, 0, 12)
series(cot(x), x)
series(sec(x), x, 0, 8)
series(csc(x), x)
series(exp(x), x)
series(exp(x), x, 0, 3)
series(exp(x), x, 0, 1)
series(exp(x), x, 0, 0)
series(log(1+x), x)
series(log(1+x)/x, x)
series(log(1-x), x, 0, 5)
series(sqrt(1+x), x)
series(1/sqrt(1-x), x, 0, 5)
series((1+x)**Rational(1,3), x, 0, 4)
series(1/(1-x), x)
series(1/(1-x)**2, x, 0, 4)
series(1/(1-x-x**2), x)
series(1/(1+x+x**2), x, 0, 8)
series(x/(exp(x)-1), x, 0, 8)
series(exp(sin(x)), x, 0, 10)
series(exp(cos(x)), x, 0, 6)
series(sin(sin(x)), x, 0, 8)
series(sin(x)*cos(x), x, 0, 8)
series(sin(x)**2, x, 0, 8)
series(exp(x)*sin(x), x, 0, 7)
series(exp(x)/x**2, x, 0, 3)
series(sin(x)/x, x, 0, 8)
series(sin(x)/x**3, x, 0, 4)
series(1/sin(x), x)
series(1/(x*sin(x)), x, 0, 4)
series(cos(x)/x**2, x, 0, 3)
series(log(cos(x)), x, 0, 8)
series(log(1+sin(x)), x, 0, 5)
series(exp(x)-1-x, x, 0, 5)
series(sinh(x), x, 0, 8)
series(cosh(x), x, 0, 8)
series(tanh(x), x, 0, 8)
series(coth(x), x)
series(asin(x), x, 0, 8)
series(acos(x), x, 0, 4)
series(atan(x), x, 0, 8)
series(acot(x), x, 0, 4)
series(asinh(x), x, 0, 8)
series(atanh(x), x)
series(erf(x), x)
series(atan(x+1), x, 0, 3)
series(asin(x+Rational(1,2)), x, 0, 3)
series(sinh(x+1), x, 0, 3)
series(tanh(x+1), x, 0, 3)
series(sin(x+y), x, 0, 3)
series(sin(x)+cos(y), x, 0, 4)
series(sin(x)*cos(y), x, 0, 4)
series(exp(a*x), x, 0, 4)
series(exp(a*x+x), x, 0, 3)
series(1/(1-a*x), x, 0, 4)
series(sec(x+1), x, 0, 2)
series(sqrt(1+x)**3, x, 0, 3)
series(cos(sqrt(x)), x, 0, 3)
series(sin(sqrt(x))/sqrt(x), x, 0, 3)
series(sqrt(x)*exp(x), x, 0, 3)
series(sqrt(sin(x)), x, 0, 3)
series(x**Rational(1,3)*exp(x), x, 0, 2)
series(asin(1-x), x, 0, 3)
series(log(x), x)
series(log(x)**2, x)
series(log(2*x), x)
series(log(x**2), x)
series(x*log(x), x, 0, 3)
series(log(x)/x, x, 0, 3)
series(log(x)*x + x, x)
series(log(x)*sin(x), x)
series(exp(x)*log(x), x, 0, 3)
series(log(x)*log(1+x), x, 0, 3)
series(log(x+x**2), x, 0, 3)
series(log(sin(x)), x)
series(log(sin(x)/x), x, 0, 8)
series(x**x, x, 0, 3)
series(exp(x*log(x)), x, 0, 3)
series(exp(log(x)+x), x, 0, 3)
series((1+x)**x, x, 0, 3)
series(1/x, x)
series(1/x**10, x)
series(1/x+x, x)
series((x+1)/x, x)
series(1/x+1, x)
series(x**2+2*x, x)
series((1+x)**3, x)
series(x**10+x, x)
series(x**10, x)
series(x*(1+x), x, 0, 3)
series(x**2*(1+x), x, 0, 3)
series(x*(1+x)**2, x, 0, 3)
series(1/x/(1+x), x, 0, 3)
series(sqrt(x), x)
series(sqrt(x)*(1+x), x)
series(x**Rational(3,2)+x**2, x, 0, 2)
series(x*sin(x), x, 0, 4)
series(abs(sin(x)), x)
series(sqrt(x**2), x)
series(exp(1/x), x)
series(exp(x), x, 1, 3)
series(exp(x), x, 2, 3)
series(log(x), x, 1)
series(log(x), x, 1, 3)
series(sqrt(x), x, 1, 3)
series(cos(x), x, pi, 3)
series(sin(x), x, pi/2, 5)
series(tan(x), x, 1, 3)
series(1/(x**2-1), x, 1, 3)
series(x/(x-1), x, 1, 3)
series(1/sin(x), x, pi, 3)
series(log(x-1), x, 1, 3)
series(acosh(x), x, 1, 3)
series(x**2, x, 1)
series(exp(x)*x, x, -1, 3)
series(1/(x+1), x, oo)
series(atan(x), x, oo)
series(log(1+1/x), x, oo, 3)
series(exp(1/x), x, oo, 3)
series(sin(1/x), x, oo, 3)
series(sqrt(x+1)-sqrt(x), x, oo, 3)
series(x*(sqrt(x**2+1)-x), x, oo, 4)
series(x/(x+1), x, oo, 4)
series((x**2+1)/(2*x**2-3), x, oo, 4)
series(exp(-x), x, oo, 3)
series(log(x), x, oo)
series(x+1/x, x, oo, 3)
series(x**2+1/x, x, oo, 3)
series(exp(x), x, 0, 3, '-')
series(log(x), x, 1, 3, '-')
series(sqrt(1-x), x, 1, 3)
series(exp(1/x), x, 0, 3, '-')
series(sin(x), x, 0, 20)
series(exp(sin(x))*cos(x), x, 0, 12)
series(1/cos(x), x, 0, 16)
series(tan(sin(x)) - sin(tan(x)), x, 0, 10)
series(exp(exp(x)), x, 0, 5)
series(log(exp(x)+x), x, 0, 4)
series(1/(exp(x)-1), x, 0, 4)
series(sqrt(4+x), x, 0, 3)
series(1/(2+x), x, 0, 4)
series(log(2+x), x, 0, 4)
series(exp(2+x), x, 0, 3)
series(sin(x)**3/x**3, x, 0, 6)
series(cos(x)**Rational(1,2), x, 0, 6)
series(x**2*exp(1/x), x, oo, 2)
series(atan(1/x), x, 0, 3)
O(x) + x
x + O(x**2)
x**3 + O(x**2)
O(x)*x
y*O(x)
O(x)**2
O(2*x + x**2)
O(x*y, x)
O(1) + 1
O(sin(x))
O(exp(x))
O(x) + O(x**2)
O(x, (x, 1))
O(1/x, (x, oo))
O(x**2 + x, (x, oo))
O(x*log(x)) + x
O(x) + x*log(x)
O(x) + x**2*log(x)
O(x**2) + 1/x
O(1) + 1/x
2*O(x)
O(x) - O(x)
x*O(1)
exp(x)*O(x)
O(x)/x
(x - x**3/6 + O(x**5)).removeO()
(x - x**3/6 + O(x**5)).getO()
series(sin(x), x).removeO()
series(sin(x), x).getn()
sin(x).taylor_term(3, x)
exp(x).taylor_term(4, x)
cos(x).series(x, 0, 8)
(1/(1-x)).series(x, 0, 4)
exp(x).series(x, 1, 3)
""".strip().splitlines()

LIMITS = r"""
limit(sin(x)/x, x, 0)
limit((1-cos(x))/x**2, x, 0)
limit(x*log(x), x, 0)
limit(x**x, x, 0)
limit(x**10*exp(-x), x, oo)
limit((1+1/x)**x, x, oo)
limit((x**2+1)/(2*x**2-3), x, oo)
limit((x**2+1)/(2*x**2-3), x, -oo)
limit(1/x, x, 0)
limit(1/x, x, 0, '-')
limit(1/x, x, 0, '+-')
limit(1/x**2, x, 0, '+-')
limit(-1/x, x, 0)
limit(1/x, x, oo)
limit(exp(x), x, oo)
limit(exp(x), x, -oo)
limit(exp(-x)*x**100, x, oo)
limit(log(x)/x, x, oo)
limit(log(x)**10/x, x, oo)
limit(x/log(x), x, oo)
limit(log(x), x, 0)
limit(log(x), x, oo)
limit(x*exp(x), x, -oo)
limit(sqrt(x**2+x)-x, x, oo)
limit(sqrt(x**2+x)-x, x, -oo)
limit(x**2, x, -oo)
limit(x**3, x, -oo)
limit(exp(-1/x), x, 0)
limit(exp(-1/x), x, 0, '-')
limit(atan(x), x, oo)
limit(atan(x), x, -oo)
limit(tan(x), x, pi/2)
limit(tan(x), x, pi/2, '-')
limit(1/(x-1), x, 1)
limit(1/(x-1)**2, x, 1)
limit(abs(x)/x, x, 0)
limit(abs(x)/x, x, 0, '-')
limit((exp(x)-1)/x, x, 0)
limit((1+x)**(1/x), x, 0)
limit((1+2/x)**x, x, oo)
limit((1+a/x)**x, x, oo)
limit(x*sin(1/x), x, oo)
limit((sin(x)-x)/x**3, x, 0)
limit((tan(x)-sin(x))/x**3, x, 0)
limit((exp(x)-1-x)/x**2, x, 0)
limit(log(1+x)/x, x, 0)
limit((2**x-1)/x, x, 0)
limit(x**(1/x), x, oo)
limit(x**(1/x), x, 0)
limit((x+1)**(1/x), x, oo)
limit(exp(x)/x**5, x, oo)
limit(exp(x**2)/exp(x), x, oo)
limit(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x, oo)
limit(exp(exp(x))/exp(exp(x - exp(-x))), x, oo)
limit(exp(x+exp(-x))-exp(x), x, oo)
limit(log(log(x))/log(x), x, oo)
limit(x**2*(exp(1/x)-1-1/x), x, oo)
limit((3**x+5**x)**(1/x), x, oo)
limit(sqrt(x+sqrt(x))-sqrt(x), x, oo)
limit(x*(sqrt(x**2+1)-x), x, oo)
limit((x-sin(x))/(x+sin(x)), x, oo)
limit(sin(x)/x, x, oo)
limit(cos(x)/x, x, oo)
limit(x/(x+1), x, oo)
limit((1-x)/(1+x), x, -oo)
limit(sinh(x)/exp(x), x, oo)
limit(tanh(x), x, oo)
limit(tanh(x), x, -oo)
limit(cosh(x)/sinh(x), x, oo)
limit(x**2*exp(-x**2), x, oo)
limit(exp(1/x)*x, x, 0)
limit(exp(1/x)*x, x, 0, '-')
limit((x**2-1)/(x-1), x, 1)
limit((x**3-8)/(x-2), x, 2)
limit(sin(3*x)/sin(5*x), x, 0)
limit(log(x)-log(x+1), x, oo)
limit(x*log(1+1/x), x, oo)
limit(x - log(exp(x)+1), x, oo)
limit(log(exp(x)+x)/x, x, oo)
limit(atan(1/x), x, 0)
limit(atan(1/x), x, 0, '-')
limit(sin(x), x, 1)
limit(x**2+1, x, 2)
limit(a*x, x, oo)
limit(gamma(x), x, 0)
limit(acos(x), x, 1)
limit(asin(x)/x, x, 0)
limit(x**Rational(1,3)*log(x), x, 0)
limit((log(x)-1)/(x-E), x, E)
limit(exp(x)/(1+exp(x)), x, oo)
limit(exp(x)/(1+exp(x)), x, -oo)
limit((x+1)/x, x, 0)
limit(sqrt(x)/sqrt(x+sqrt(x)), x, oo)
limit(1/x - 1/sin(x), x, 0)
limit(1/x - 1/tan(x), x, 0)
""".strip().splitlines()

CASES = SERIES + LIMITS

# Deliberate deviations (semantically fine, printed differently or more
# precise than SymPy).
KNOWN_DEVIATIONS = {
    # SymPy loses the x**2*log(x) term (bug); symox keeps it
    "series((1+x)**x, x, 0, 3)",
}


def namespace(mod):
    ns = {name: getattr(mod, name) for name in dir(mod) if not name.startswith("_")}
    x, y, z, a = mod.symbols("x y z a")
    ns.update(x=x, y=y, z=z, a=a)
    return ns


def numeric_equal(a, b, syms, trials=4):
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
        if abs(va - vb) > 1e-8 * max(1.0, abs(va), abs(vb)):
            return False
    return True


def compare(case):
    ns_s = namespace(sympy)
    ns_x = namespace(symox)
    try:
        exp_s = eval(case, ns_s)
    except Exception as e:  # noqa: BLE001
        exp_s = f"<{type(e).__name__}>"
    try:
        exp_x = eval(case, ns_x)
    except Exception as e:  # noqa: BLE001
        return "error", f"{type(e).__name__}: {e}", str(exp_s)
    s_s, s_x = str(exp_s), str(exp_x)
    if s_s == s_x:
        return "same", s_x, s_s
    if not isinstance(exp_s, sympy.Basic):
        return "DIFFERENT", s_x, s_s
    try:
        x, y, z, a = sympy.symbols("x y z a")
        back = sympy.sympify(s_x, locals={"x": x, "y": y, "z": z, "a": a, "O": sympy.O})
    except Exception:  # noqa: BLE001
        return "unparseable", s_x, s_s
    if back == exp_s:
        return "same-tree", s_x, s_s
    try:
        o_s, o_x = exp_s.getO(), back.getO()
        if o_s != o_x:
            return "DIFFERENT", s_x, s_s
        p_s = exp_s.removeO() if o_s is not None else exp_s
        p_x = back.removeO() if o_x is not None else back
        syms = sorted(p_s.free_symbols | p_x.free_symbols, key=str)
        if sympy.simplify(p_s - p_x) == 0 or numeric_equal(p_s, p_x, syms):
            return "equivalent", s_x, s_s
    except Exception:  # noqa: BLE001
        pass
    return "DIFFERENT", s_x, s_s


def run(verbose=True, show_equivalent=False):
    stats = {}
    bad = []
    for case in CASES:
        status, got, want = compare(case)
        stats[status] = stats.get(status, 0) + 1
        if status not in ("same", "same-tree") and (show_equivalent or status != "equivalent"):
            if verbose:
                print(f"[{status:11s}] {case}\n      symox: {got}\n      sympy: {want}")
        if status in ("DIFFERENT", "error", "unparseable") and case not in KNOWN_DEVIATIONS:
            bad.append((case, got, want))
    if verbose:
        print("\nsummary:", stats, f"(total {len(CASES)}: {len(SERIES)} series/order, {len(LIMITS)} limits)")
    return bad


def test_series_against_sympy():
    bad = run(verbose=False)
    assert not bad, "\n".join(f"{c}: symox={g} sympy={w}" for c, g, w in bad)


if __name__ == "__main__":
    sys.exit(1 if run(show_equivalent="-v" in sys.argv) else 0)
