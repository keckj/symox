"""Differential tests of the polynomial subsystem against SymPy.

Every case is a Python expression evaluated in a SymPy namespace and in a
symox namespace. Results are compared by ``str`` first; when the strings
differ, a semantic comparison is attempted:

* expressions: structural equality after conversion, then numerical
  equality at random complex points;
* dicts (``roots``): same keys up to numerical equality, same
  multiplicities;
* lists/tuples: element-wise;
* ``Poly``: same printed form required.

Run directly for a report, or through pytest (``test_against_sympy``).
"""

from __future__ import annotations

import cmath
import random
import sys

import sympy
import symox

CASES = r"""
factor(x**2 - 1)
factor(x**4 - 1)
factor(x**6 - 1)
factor(x**8 - 1)
factor(x**12 - 1)
factor(x**20 - 1)
factor(x**15 + 1)
factor(x**30 - 1)
factor(x**10 - 1)
factor(x**7 - 1)
factor(x**12 + 1)
factor(x**9 - x)
factor(x**2 + 2*x + 1)
factor(2*x**2 + 4*x + 2)
factor(x**3 - x)
factor(x**5 - x)
factor(x**4 + 4)
factor(x**4 + x**2 + 1)
factor(x**6 + 1)
factor(6*x**2 + 5*x + 1)
factor(6*x**3 - 11*x**2 + 6*x - 1)
factor(x**4 - 10*x**2 + 1)
factor(expand((x**2 - 2)*(x**2 - 3)))
factor(x**8 - 256)
factor(3*x**3 + 3)
factor(x/2 + 1)
factor(x**2/2 + x + Rational(1, 2))
factor(-x**2 + 1)
factor(-x - 1)
factor(2*x + 2)
factor(-2*x - 2)
factor(x**2 + 1)
factor(x)
factor(5)
factor(0)
factor(expand((x - 1)**3*(x + 2)**2))
factor(expand((x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5)*(x - 6)*(x - 7)*(x - 8)*(x - 9)*(x - 10)))
factor(expand((2*x + 3)**3*(5*x - 7)**2*(x**2 + x + 1)))
factor(expand((x**3 + 2*x + 7)*(x**4 - 3*x + 1)*(2*x**2 - 5)))
factor(x**16 - 1)
factor(x**24 - 1)
factor(x**18 + 1)
factor(x**2 - y**2)
factor(x**3*y - x*y**3)
factor(expand((x + y)**2*(x - y)))
factor(x*y + x*z + y**2 + y*z)
factor(x**2*y**2 - 1)
factor(x**4 - y**4)
factor(x**3 + y**3)
factor(x**3 - y**3)
factor(expand((x + y + z)**3))
factor(expand((x*y + 1)*(x - y)))
factor(6*x**2*y + 4*x*y**2)
factor(x**2 + 2*x*y + y**2 - z**2)
factor(expand((x + 2*y)*(3*x - y)*(x + y + 1)))
factor(x**2*y - y)
factor(x**6 - y**6)
factor(expand((x**2 + y)*(y**2 + x)))
factor(a*x + a*y + b*x + b*y)
factor(expand((x*y*z + 1)*(x + y + z)))
factor(expand((x**2*y - 3*z)*(x*z**2 + y + 2)))
factor(expand((x + y)**5))
factor(x**4*y**2 - x**2*y**4)
factor(expand((a + b)*(a - b)*(a**2 + b**2 + 1)))
factor((x**2 - 1)/(x**2 + 2*x + 1))
factor(1/x + 1/y)
factor((x**2 - 4)/(x**2 + 4*x + 4))
factor(x/(x**2 - x))
factor(sin(x)**2 - 1)
factor(exp(2*x) - 1)
factor(sqrt(2)*x + sqrt(2))
factor(x**2*sin(y) - sin(y))
factor(cos(x)**2 + 2*cos(x) + 1)
factor(x**2 + 4*x*sqrt(3) + 12)
factor_list(2*x**2 - 2)
factor_list(x**5 - x)
factor_list(-x**2 + x)
factor_list(x**2/2 - Rational(1, 2))
factor_list(x**4*y - x**2*y**3)
factor_list(x*y + x)
factor_list(x**2*y - y)
factor_list(expand((x - y)**2*(x + y)))
factor_list(x**20 - 1)
factor_list(expand(3*(x + 1)**2*(x - 2)))
factor_list(5)
sqf_list(x**3 - x**2 - x + 1)
sqf_list(expand((x - 1)**2*(x + 2)**3))
sqf_list(expand(2*(x + 1)**2*(x**2 + 1)))
sqf_list(expand((x + y)**2*(x - y)))
sqf(x**3 + 2*x**2 + x)
sqf(expand((x - 1)**2*(x + 2)**3))
gcd(x**2 - 1, x**2 + 2*x + 1)
gcd(2*x**2 - 2, 4*x + 4)
gcd(-x + 1, x**2 - 1)
gcd(6, 4)
gcd(x*y**2, x**2*y)
gcd(x/2 + 1, x + 2)
gcd(expand((x + y)**3*(x - y)), expand((x + y)**2*(x + 2*y)))
gcd(x**4 - 1, x**6 - 1)
gcd(x**3 + 1, x**2 + 3*x + 2)
gcd(Rational(1, 2), Rational(1, 3))
gcd(x**2 - 1, x**2 - 3*x + 2)
gcd(expand((x*y + z)**2*(x - 1)), expand((x*y + z)*(x + 1)))
gcd(expand((x + 1)**4*(x - 3)**2), expand((x + 1)**2*(x - 3)**5))
gcd([x**2 - 1, x**2 - 2*x + 1, x**3 - 1])
gcd(2*x, 4)
lcm(x**2 - 1, x + 1)
lcm(6, 4)
lcm(x*y, y*z)
lcm(x**2 - y**2, x + y)
lcm(x/2, x + 1)
lcm(2*x + 2, 3*x - 3)
cancel((x**2 - 1)/(x - 1))
cancel((x**2 - 2)/(x - sqrt(2)))
cancel((2*x + 2)/(4*x**2 - 4))
cancel(x/2 + 1/(3*x))
cancel((x*y - y)/(x**2 - 1))
cancel(1/(-x + 1))
cancel((-x - 1)/(x**2 - 1))
cancel((x**3 - y**3)/(x - y))
cancel((x**2 + 2*x*y + y**2)/(x + y))
cancel(1/(1/x + 1))
cancel((exp(2*x) - 1)/(exp(x) + 1))
cancel((x**4 - 1)/(x**2 + 1))
cancel((x**2 - x)/x**2)
cancel((a*x + a)/(b*x + b))
cancel(0)
cancel(x)
cancel((x + 1)**2/(x**2 - 1))
cancel((x**2 + 3*x + 2)/(x**2 + 4*x + 3))
cancel(x**2/(2*x))
cancel((6*x**2 + 3*x)/(9*x))
cancel(1/x + 1/y)
cancel((x**5 - 1)/(x - 1))
cancel((sin(x)**2 - 1)/(sin(x) - 1))
apart(1/(x**2 - 1))
apart((x**2 + 1)/(x**3 - x))
apart(1/(x*(x + 1)**2))
apart(x**3/(x**2 - 1))
apart(1/(x**2 + 1)/(x - 1))
apart(1/(x**3 + 1))
apart((3*x + 5)/(x**2 + x - 2))
apart(1/(x**4 - 1))
apart((x + 1)/(x - 1))
apart(1/(2*x + 2))
apart(x/(x + 1)**3)
apart((x**4 + 1)/(x**2*(x - 1)))
apart(1/(x**2*(x**2 + 1)))
collect(a*x**2 + b*x**2 + a*x - b*x + c, x)
collect(x**2*y + x*y + x + y, [x, y])
collect(expand((x + y + 1)**2), x)
collect(a*sin(x) + b*sin(x) + c, sin(x))
collect(x**3*a + x**3*b + x*a + c, x)
div(x**2 + 1, 2*x)
rem(x**3 + 1, x**2 - x)
quo(x**3 + 1, x**2 - x)
div(x**3 - 2*x + 1, x - 1)
div(3*x**4 + 2*x + 1, x**2 + 1)
rem(x**5, x**2 + x + 1)
quo(x**4 - 1, x - 1)
div(x**2 + y**2, x + y)
resultant(x**2 - 1, x - 2)
resultant(x**2 + y, x*y - 1, x)
resultant(x**3 - 2, x**2 + x + 1)
resultant(x**2 - 2*x + 1, x**2 - 1)
discriminant(x**2 + b*x + c, x)
discriminant(x**3 + x + 1)
discriminant(a*x**2 + b*x + c, x)
discriminant(x**4 - 1)
roots(x**2 - 2)
roots(x**3 - 1)
roots(x**2 + x + 1)
roots((x - 1)**2*(x + 2))
roots(x**4 - 5*x**2 + 4)
roots(x**3 - 2)
roots(x**2 - x - 1)
roots(x**4 - 1)
roots(x**6 - 1)
roots(2*x**2 - 3*x + 1)
roots(x**3 - 6*x**2 + 11*x - 6)
roots(x**4 + 1)
roots(x**3 + x + 1)
roots(x**4 - 10*x**2 + 1)
roots(x**2 + 1)
roots(x**5 - x)
roots(x**2 - a, x)
roots(x**2 + b*x + c, x)
roots(x**8 - 1)
roots(x**3 - 3*x + 1)
roots(4*x**2 + 4*x + 1)
roots(x**2 + 2*x + 5)
roots(x**4 + x**2 + 1)
roots(x**5 + 1)
roots(x)
roots(x**3)
roots(x**2*(x - 1))
nroots(x**3 - 2)
nroots(x**2 + 1)
nroots(x**5 - x - 1)
nroots(x**4 + x + 1)
nroots(x**3 - 6*x**2 + 11*x - 6)
nroots(x**2 - 2)
nroots(3*x**3 - x + 7)
Poly(x**2 + 1, x)
Poly(x*y + 1, x)
Poly(x/2 + 1)
Poly(x*y + y**2*z + 3)
Poly(x**2 + 2*x*y + 1, y)
Poly(a*x**2 + b, x)
Poly(sqrt(2)*x)
Poly(pi*x + sin(y))
Poly(1/x + x)
Poly(x**2*y - 3*x + y/3, x, y)
Poly(x + y + z + a + b + t + c)
Poly(x**2 + 1)*2
Poly(x + 1)**2
Poly(x, x) + Poly(y, y)
Poly(x**2 + 1) - Poly(x)
Poly(x**2 + 1).degree()
Poly(x**3*y + y, x, y).degree(y)
Poly(2*x**3 + x).coeffs()
Poly(2*x**3 + x).all_coeffs()
Poly(x**2 - 1).factor_list()
Poly(x**2 + 2*x + 1).LC()
div(Poly(x**2 + 1), Poly(2*x))
Poly(x**2 - 1).gcd(Poly(x - 1))
Poly(x**2 - 1).as_expr()
Poly(x**2*y + 1).terms()
Poly(x**2 + 1).gens
Poly(x**3 + x).diff(x)
Poly(x**2 - 3*x + 2).eval(1)
Poly(x**2*y + x*y**2).monoms()
Poly(6*x**2 + 3*x).primitive()
Poly(x**2 + 1, x, domain='QQ')
Poly(2*x**2 + 4).monic()
Poly(x**4 - 1).sqf_list()
degree(x**3*y + y**2, y)
degree(x**2 + 1)
degree(5)
degree(x**5 + x**2*y, x)
factor(x**2/4 - 1)
factor(x**2 + I*x)
factor(10**20*x**2 - 10**20)
factor(x**4 + x**3 + x**2 + x + 1)
factor(x + 1/x)
factor(exp(x)**2 - 1)
factor((x + 1)*(x + 2))
factor(x**2 + x*y + y**2)
factor(expand((x - y + 2)**3*(x + y)))
factor(x**3 - 3*x**2*y + 3*x*y**2 - y**3)
factor(4*x**4 - 4*x**2*y**2 + y**4)
factor(x**2 - 2)
factor(25*x**2 - 49)
factor(x**6 - 2*x**3 + 1)
factor(expand((x*y**2 + 1)*(x**2*y + x + 1)))
factor(expand((x**2*y*z + y + 1)*(x*z**2 - y**3 + z)))
factor(x**2*y**2 + x**2 + y**2 + 1)
factor(sqrt(x)**2 - 1)
factor(x**3/3 + x**2/2 + x/6)
factor(-3*x**2 + 3)
factor(-x**2*y + y)
factor(Rational(3, 4)*x**2 - Rational(3, 4))
factor(x**4 - 2*x**2 + 1)
factor(2*x**2*y - 8*y)
factor(x**5 + x + 1)
factor(x**8 + x**4 + 1)
factor(expand((x**2 + 1)**3*(x - 5)))
factor(x*y - x - y + 1)
factor(x**2 - x*y - 2*y**2)
factor(Eq(x**2 - 1, 0))
factor(x**8 - 40*x**6 + 352*x**4 - 960*x**2 + 576)
factor(x**50 - 1)
factor(x**36 - 1)
factor(x**4*y - y)
factor(sin(x)*cos(x)**2 - sin(x))
factor(sin(x)**2 + 2*sin(x)*cos(x) + cos(x)**2, deep=True)
cancel(x/y)
cancel(1/(x*y) + 1/(x*z))
cancel((x**2 - 1)/(x**2 - x))
cancel((x**2 + x)/(x**2 - 1) + 1/(x - 1))
cancel(sqrt(2)*x/sqrt(2))
cancel((2*x**2 - 2)/(x + 1)/3)
gcd(x**2*y + x*y**2, x*y)
gcd(0, x**2 + 1)
gcd(x**2 + 1, 0)
gcd(x, x)
gcd(-x**2 + 1, -x + 1)
gcd(3*x**2 - 3, 6*x + 6)
lcm(x**2 + x, x**2 - x)
lcm(0, x)
apart(1/(x**2 + 2*x + 1))
apart((2*x**2 + 3)/(x**3 + x))
apart(1/(x*(x - 1)*(x - 2)))
apart((x**3 + 2)/((x - 1)**2*(x + 1)))
apart(1/(x**2 - 2))
apart(x**2/(x**2 + 1))
apart(1/((x + 1)*(x**2 + x + 1)))
apart(1/((x + a)*(x + b)), x)
apart(x/((x - a)*(x + a)), x)
apart(1/(x**2 - a**2), x)
apart((a*x + 1)/(x*(x + a)), x)
apart(1/(x**2*(x + a)), x)
apart(x**2/(x + a), x)
apart(1/((x + a)**2*(x + b)), x)
apart((x + c)/((x + a)*(x + b)), x)
apart(1/(a*x**2 + b*x), x)
apart(1/(x*(a*x + 1)), x)
apart(1/(x**2 + a), x)
apart(a/(x**2 - 1), x)
roots(x**3 + 2*x**2 + 3*x + 4)
roots(x**6 - 2)
roots(2*x**3 - 16)
roots(x**3 - x**2 - x - 1)
roots(x**4 - 4*x**2 + 1)
roots(x**4 + x**3 + x**2 + x + 1)
roots(x**2 - 2*x + 1)
roots(x**3 - 3*x**2 + 3*x - 1)
roots(9*x**2 - 4)
roots(x**2/2 - 1)
roots(x**4 - 3*x**2 + 2)
roots(x**3 + x)
roots(x**6 + 1)
roots(x**12 - 1)
roots(x**5 + 1)
roots(x**8 - 1)
roots(x**10 - 1)
roots(x**5 - 1)
roots(x**8 + 1)
roots(x**5 - 32)
roots(x**4 + 16)
roots(x**3 + 8)
roots(x**6 + 64)
roots(x**2 - 3*x + 2, multiple=True)
roots(x**3 - 7*x + 6, filter='Z')
nroots(x**8 - 1)
nroots(x**2 - x - 1)
nroots(x**3 + x**2/3 - 2)
nroots((x - 1)**2*(x + 2))
real_roots(x**3 - 7*x + 6)
degree(x**2*y, x)
Poly(x**2*y).degree(y)
Poly(x**2 + x*y, y)
Poly(3*x**2 + 2*x + 1).all_coeffs()
Poly(x**3 - x, domain='QQ').factor_list()
Poly(x**2 + 1) + 1
Poly(x**2 + 1) * x
Poly(x**2 + y**2, x, y).coeffs()
Poly(0, x)
Poly(5, x)
Poly(x**2 + 1).rem(Poly(x + 1))
Poly(x**4 - 1).quo(Poly(x - 1))
Poly(x**2 - 1).exquo(Poly(x - 1))
Poly(x**2 + 2*x + 1).sqf_list()
Poly(x*y**2 + x**2*y + 1).total_degree()
Poly(x**2 + 1).is_irreducible
Poly(x**2 - 1).is_irreducible
Poly(2*x + 1).LC()
Poly(2*x + 1).TC()
Poly(x**2 + 1).length()
sqf_list(x**5 - x**4 - x + 1)
sqf_list(2*x**2 + 4*x + 2)
sqf_list(x**2*(x + 1)**3)
sqf(x**2 + 2*x + 1)
factor_list(x**2 + 1)
factor_list(x)
factor_list(2*x)
factor_list(-2*x**2 + 2)
factor_list(x**2*y + x*y**2)
factor_list(expand((x + 1)**3*(y - 2)**2))
factor_list(-716*x*z)
factor_list(432*w**2*x**2*y**6)
resultant(x**2 + 1, x**2 - 1)
resultant(2*x**2 + 1, 3*x - 1)
resultant(x*y + 1, x + y, y)
discriminant(x**2 - 2)
discriminant(3*x**3 + 2*x**2 + x + 5)
discriminant(x**5 - x - 1)
collect(x*y + x*z + y, x)
collect(x**2 + x**2*y + x + 1, x)
collect(a*x + b*x + a*y + b*y, [x, y])
collect(exp(x)*a + exp(x)*b, exp(x))
div(x**4 + 1, x**2 + 1)
div(2*x**3 + 4*x, 2*x)
quo(x**2 - y**2, x - y)
rem(x**2 + y, x + y, x)
div(x*y + y**2, x + y, y)
pdiv(x**2 + 1, 2*x - 4)
prem(x**2 + 1, 2*x - 4)
pquo(x**2 + 1, 2*x - 4)
pdiv(x**2 + y, x*y - 1)
prem(x**3*y + x, x*y - 1, x)
pexquo(x**2 - 1, 2*x - 2)
prem(x**2 + y**2, x - y)
pdiv(3*x**3 + 1, 2*x + 1)
LC(3*x**2 + y)
LM(3*x**2 + y)
LT(3*x**2*y + y)
primitive(6*x**2 + 4)
together(1/x + 1/y)
together(x/(x + 1) + 1/x)
together(1/(x + 1) + 1/(x - 1))
together(x + 1/x)
together(1 + exp(-x))
together(1/x**2 + 1/x)
together(sin(x)/2 + cos(x)/3)
together(sqrt(x)/x)
simplify((x**2 - 1)/(x - 1))
simplify((x + 1)**2 - 1)
simplify(x/(x**2 + x))
simplify((x**3 - 1)/(x**2 + x + 1))
simplify(2*x + 2*y)
simplify(x*y + x*z)
simplify((x**2 + 3*x + 2)/(x + 1))
simplify(1/(x + 1) + 1/(x - 1))
""".strip().splitlines()


# Deliberate deviations from SymPy (documented in the report).
KNOWN_DEVIATIONS: set[str] = set()


def namespace(mod):
    ns = {name: getattr(mod, name) for name in dir(mod) if not name.startswith("_")}
    x, y, z, t, w, a, b, c = mod.symbols("x y z t w a b c")
    ns.update(x=x, y=y, z=z, t=t, w=w, a=a, b=b, c=c)
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


def to_sympy(obj):
    """Convert a symox result into the corresponding SymPy object."""
    if isinstance(obj, dict):
        return {to_sympy(k): to_sympy(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return type(obj)(to_sympy(v) for v in obj)
    if isinstance(obj, symox.Poly):
        return sympy.sympify(str(obj.as_expr()))
    if isinstance(obj, symox.Expr):
        return sympy.sympify(str(obj))
    return obj


def sem_equal(s, x):
    """Semantic equality of a SymPy result and a converted symox result."""
    if isinstance(s, sympy.Poly):
        return False
    if isinstance(s, dict) and isinstance(x, dict):
        if len(s) != len(x):
            return False
        used = set()
        for kx, vx in x.items():
            found = False
            for i, (ks, vs) in enumerate(s.items()):
                if i in used or vs != vx:
                    continue
                if sem_equal(ks, kx):
                    used.add(i)
                    found = True
                    break
            if not found:
                return False
        return True
    if isinstance(s, (list, tuple)) and isinstance(x, (list, tuple)):
        return len(s) == len(x) and all(sem_equal(a, b) for a, b in zip(s, x))
    if isinstance(s, sympy.Basic) and isinstance(x, sympy.Basic):
        if s == x:
            return True
        syms = sorted(s.free_symbols | x.free_symbols, key=str)
        if not syms:
            try:
                return abs(complex(s.evalf(30)) - complex(x.evalf(30))) <= 1e-10 * max(1.0, abs(complex(s.evalf())))
            except (TypeError, ValueError):
                return False
        try:
            return numeric_equal(s, x, syms)
        except Exception:  # noqa: BLE001
            return False
    return s == x


def compare(case):
    exp_s = eval(case, namespace(sympy))
    try:
        exp_x = eval(case, namespace(symox))
    except Exception as e:  # noqa: BLE001
        return "error", f"{type(e).__name__}: {e}", str(exp_s)
    s_s, s_x = str(exp_s), str(exp_x)
    if s_s == s_x:
        return "same", s_x, s_s
    try:
        back = to_sympy(exp_x)
    except Exception:  # noqa: BLE001
        return "unparseable", s_x, s_s
    if str(back) == s_s:
        return "same-tree", s_x, s_s
    if sem_equal(exp_s, back):
        return "equivalent", s_x, s_s
    return "DIFFERENT", s_x, s_s


def run(verbose=True):
    stats = {}
    bad = []
    for case in CASES:
        status, got, want = compare(case)
        stats[status] = stats.get(status, 0) + 1
        if status not in ("same", "same-tree"):
            if verbose:
                print(f"[{status:11s}] {case:50s}\n      symox: {got}\n      sympy: {want}")
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
