"""Accuracy of symox' double precision special functions against mpmath.

Every function is evaluated at many points through ``symox._core._numeric``
(or the expression API for the core heads ``gamma``, ``loggamma``, ``erf``,
``erfc``) and compared with mpmath at 30 digits. Run directly for a report,
or through pytest.
"""

from __future__ import annotations

import cmath
import math
import sys

import mpmath as mp

import symox
from symox import _core

# References are computed at 30 digits inside ``run`` only, so importing this
# module during collection does not change the global mpmath precision.
DPS = 30

REL = 1e-13  # tolerance used by the assertion (typical errors are ~1e-15)


def lin(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


X_POS = [0.05, 0.1, 0.3, 0.5, 0.77, 1.0, 1.3, 1.9, 2.5, 3.7, 5.2, 7.5, 9.9, 13.4, 18.0, 24.0, 31.5, 47.0]
X_ALL = [-7.3, -3.6, -2.5, -1.2, -0.5, -0.1] + X_POS
X_UNIT = [-0.95, -0.7, -0.4, -0.13, 0.0, 0.2, 0.45, 0.6, 0.85, 0.99]
Z_C = [complex(0.5, 0.5), complex(1.2, -0.7), complex(-0.8, 1.3), complex(2.5, 1.5), complex(-1.5, -2.2),
       complex(0.1, 3.0), complex(4.0, -0.3)]
ORDERS = [0.0, 1.0, 2.0, 0.5, 1.5, -0.5, 2.5, 3.3, -1.7, 5.0, 7.5, 10.0]


def sx_core(name):
    f = getattr(symox, name)

    def g(*a):
        v = f(*[symox.Float(x) if isinstance(x, float) else symox.sympify(x) for x in a])
        return complex(v)

    return g


def sx(name):
    def g(*a):
        return _core._numeric(name, *a)

    return g


def ref_catalan(x):
    return 4 ** x * mp.gamma(x + 0.5) / (mp.sqrt(mp.pi) * mp.gamma(x + 2))


def ref_lucas(x):
    phi = (1 + mp.sqrt(5)) / 2
    return phi ** x + mp.cospi(x) * phi ** (-x)


def ref_jn(n, z):
    return mp.sqrt(mp.pi / (2 * z)) * mp.besselj(n + 0.5, z)


def ref_yn(n, z):
    return mp.sqrt(mp.pi / (2 * z)) * mp.bessely(n + 0.5, z)


def pts1(xs):
    return [(x,) for x in xs]


def pts2(a, b):
    return [(x, y) for x in a for y in b]


CASES = [
    # (label, symox callable, reference, points)
    ("gamma", sx_core("gamma"), mp.gamma, pts1([x for x in X_ALL + [0.001, 60.5, 120.3, 170.2]])),
    ("gamma_c", sx_core("gamma"), mp.gamma, pts1(Z_C)),
    ("loggamma", sx_core("loggamma"), mp.loggamma, pts1(X_ALL + [150.0, 1e4])),
    ("loggamma_c", sx_core("loggamma"), mp.loggamma, pts1(Z_C)),
    ("erf", sx_core("erf"), mp.erf, pts1(lin(-5, 5, 41))),
    ("erfc", sx_core("erfc"), mp.erfc, pts1(lin(-3, 26, 59))),
    ("erf_c", sx_core("erf"), mp.erf, pts1(Z_C)),
    ("erfi", sx("erfi"), mp.erfi, pts1(lin(-12, 12, 49))),
    ("erfinv", sx("erfinv"), mp.erfinv, pts1([-0.999999, -0.9, -0.5, -0.1, 1e-5, 0.3, 0.7, 0.95, 0.999, 0.9999999])),
    ("erfcinv", sx("erfcinv"), lambda x: mp.erfinv(1 - mp.mpf(x)), pts1([1e-12, 1e-5, 0.01, 0.3, 0.9, 1.3, 1.9])),
    ("factorial", sx("factorial"), mp.factorial, pts1([0.5, 2.5, 7.7, -0.5, -2.3, 30.1])),
    ("binomial", sx("binomial"), mp.binomial, pts2([2.5, 7.0, -3.5, 10.2], [0.5, 2.0, 3.7])),
    ("rf", sx("RisingFactorial"), mp.rf, pts2([0.5, 2.5, -3.5, 10.2], [0.5, 2.0, 3.7, -1.0])),
    ("ff", sx("FallingFactorial"), mp.ff, pts2([0.5, 2.5, 10.2], [0.5, 2.0, 3.7])),
    ("harmonic", sx("harmonic"), mp.harmonic, pts1([0.5, 2.5, 10.3, 100.7])),
    ("catalan", sx("catalan"), ref_catalan, pts1([0.5, 2.5, 7.3])),
    ("fibonacci", sx("fibonacci"), mp.fibonacci, pts1([0.5, 2.5, 7.3, -1.5])),
    ("lucas", sx("lucas"), ref_lucas, pts1([0.5, 2.5, 7.3, -1.5])),
    ("digamma", sx("digamma"), mp.digamma, pts1(X_ALL)),
    ("digamma_c", sx("digamma"), mp.digamma, pts1(Z_C)),
    ("polygamma", sx("polygamma"), mp.polygamma, pts2([1.0, 2.0, 3.0, 5.0], X_POS + [-2.5, -0.3])),
    ("beta", sx("beta"), mp.beta, pts2([0.5, 2.5, 7.0], [0.3, 1.5, 9.5])),
    ("lowergamma", sx("lowergamma"), lambda a, x: mp.gammainc(a, 0, x), pts2([0.5, 1.0, 2.5, 7.3, 20.0], [0.1, 1.0, 3.0, 8.0, 25.0])),
    ("uppergamma", sx("uppergamma"), lambda a, x: mp.gammainc(a, x), pts2([0.5, 1.0, 2.5, 7.3, -0.5, -2.0, 0.0], [0.1, 1.0, 3.0, 8.0, 25.0])),
    ("zeta", sx("zeta"), mp.zeta, pts1([-11.5, -5.0, -2.5, -0.5, 0.3, 0.5, 0.9, 1.1, 1.5, 2.0, 3.0, 7.5, 30.0])),
    ("zeta_c", sx("zeta"), mp.zeta, pts1([complex(0.5, 14.0), complex(1.5, 2.0), complex(-1.5, 3.0), complex(2.0, -1.0)])),
    ("hurwitz", sx("zeta"), mp.zeta, pts2([1.5, 2.0, 3.0, 4.5], [0.3, 1.0, 2.5, 10.0])),
    ("dirichlet_eta", sx("dirichlet_eta"), mp.altzeta, pts1([-3.5, -1.0, 0.0, 0.5, 1.0, 1.0001, 2.0, 2.5, 10.0])),
    ("polylog", sx("polylog"), mp.polylog, pts2([2.0, 3.0, 1.0, 0.0, -1.0, -2.0, 4.0], [-3.0, -0.9, -0.3, 0.2, 0.5, 0.7, 0.95])),
    ("polylog_real_s", sx("polylog"), mp.polylog, pts2([2.5, 0.5, 1.5], [-0.9, -0.3, 0.2, 0.5, 0.7, 0.95])),
    ("polylog_big", sx("polylog"), mp.polylog, pts2([2.0, 3.0], [2.0, 5.0, -3.0, -10.0])),
    ("lerchphi", sx("lerchphi"), mp.lerchphi, [(0.5, 2.0, 1.5), (-0.3, 1.5, 0.7), (0.9, 3.0, 2.0)]),
    ("fresnels", sx("fresnels"), mp.fresnels, pts1(lin(-7, 7, 57))),
    ("fresnelc", sx("fresnelc"), mp.fresnelc, pts1(lin(-7, 7, 57))),
    ("Ei", sx("Ei"), mp.ei, pts1([x for x in lin(-30, 60, 61) if x != 0])),
    ("Ei_c", sx("Ei"), mp.ei, pts1(Z_C)),
    ("E1", sx("E1"), mp.e1, pts1([0.01, 0.3, 1.0, 1.5, 3.0, 10.0, 40.0])),
    ("li", sx("li"), mp.li, pts1([0.1, 0.5, 0.9, 1.5, 2.5, 10.0, 1e5])),
    ("Li", sx("Li"), lambda x: mp.li(x, offset=True), pts1([1.5, 2.5, 10.0, 1e5])),
    ("Si", sx("Si"), mp.si, pts1(lin(-30, 30, 61))),
    ("Ci", sx("Ci"), mp.ci, pts1([x for x in lin(-30, 30, 61) if x != 0])),
    ("Shi", sx("Shi"), mp.shi, pts1(lin(-50, 50, 41))),
    ("Chi", sx("Chi"), mp.chi, pts1([x for x in lin(-50, 50, 40) if x != 0])),
    ("expint", sx("expint"), mp.expint, pts2([0.0, 2.0, 3.0, 5.0, 0.5, 2.5], [0.1, 0.9, 1.5, 4.0, 20.0])),
    ("besselj", sx("besselj"), mp.besselj, pts2(ORDERS, X_POS)),
    ("bessely", sx("bessely"), mp.bessely, pts2(ORDERS, X_POS)),
    ("besseli", sx("besseli"), mp.besseli, pts2(ORDERS, X_POS[:-2])),
    ("besselk", sx("besselk"), mp.besselk, pts2(ORDERS, X_POS)),
    ("besselj_c", sx("besselj"), mp.besselj, pts2([0.0, 1.0, 2.5], Z_C)),
    ("besseli_c", sx("besseli"), mp.besseli, pts2([0.0, 1.0, 2.5], Z_C)),
    ("besselk_c", sx("besselk"), mp.besselk, pts2([0.0, 1.0, 2.5], [z for z in Z_C if z.real > 0])),
    ("hankel1", sx("hankel1"), mp.hankel1, pts2([0.0, 1.5], [0.5, 2.5, 30.0])),
    ("jn", sx("jn"), ref_jn, pts2([0.0, 1.0, 3.0], [0.5, 2.5, 30.0])),
    ("yn", sx("yn"), ref_yn, pts2([0.0, 1.0, 3.0], [0.5, 2.5, 30.0])),
    ("airyai", sx("airyai"), mp.airyai, pts1(lin(-15, 10, 51))),
    ("airybi", sx("airybi"), mp.airybi, pts1(lin(-15, 10, 51))),
    ("airyaiprime", sx("airyaiprime"), lambda x: mp.airyai(x, derivative=1), pts1(lin(-15, 10, 51))),
    ("airybiprime", sx("airybiprime"), lambda x: mp.airybi(x, derivative=1), pts1(lin(-15, 10, 51))),
    ("elliptic_k", sx("elliptic_k"), mp.ellipk, pts1([-5.0, -0.5, 0.0, 0.3, 0.7, 0.99, 2.0])),
    ("elliptic_e", sx("elliptic_e"), mp.ellipe, pts1([-5.0, -0.5, 0.0, 0.3, 0.7, 0.99])),
    ("elliptic_e2", sx("elliptic_e"), mp.ellipe, pts2([0.3, 1.2, 2.5, -4.0], [-0.5, 0.3, 0.8])),
    ("elliptic_f", sx("elliptic_f"), mp.ellipf, pts2([0.3, 1.2, 2.5, -4.0], [-0.5, 0.3, 0.8])),
    ("elliptic_pi", sx("elliptic_pi"), mp.ellippi, pts2([-0.5, 0.3, 0.6], [-0.5, 0.3, 0.8])),
    ("elliptic_pi3", sx("elliptic_pi"), mp.ellippi, [(0.3, 0.5, 0.4), (-0.5, 1.2, 0.3), (0.6, 2.5, 0.2)]),
    ("LambertW", sx("LambertW"), mp.lambertw, pts1([-0.3678, -0.3, -0.1, 0.0, 0.5, 1.0, 10.0, 1e10])),
    ("LambertW_m1", sx("LambertW"), lambda x, k: mp.lambertw(x, int(k)), pts2([-0.3678, -0.3, -0.1, -1e-5], [-1.0])),
    ("LambertW_c", sx("LambertW"), mp.lambertw, pts1(Z_C + [-1.0, -5.0])),
    ("legendre", sx("legendre"), mp.legendre, pts2([2.5, 0.3, 3.0, 7.7], [-0.4, 0.0, 0.5, 0.9])),
    ("chebyshevt", sx("chebyshevt"), mp.chebyt, pts2([2.5, 0.3, 3.0, 7.7], [-0.4, 0.0, 0.5, 0.9, 1.5])),
    ("chebyshevu", sx("chebyshevu"), mp.chebyu, pts2([2.5, 0.3, 3.0, 7.7], [-0.4, 0.0, 0.5, 0.9, 1.5])),
    ("hermite", sx("hermite"), mp.hermite, pts2([2.5, 0.3, 3.0, 7.7], [-1.4, 0.0, 0.5, 1.9])),
    ("ynm", sx("Ynm"), mp.spherharm, [(2.0, 1.0, 0.5, 0.3), (3.0, -2.0, 1.2, 2.0), (5.0, 3.0, 2.0, -1.0)]),
]


def relerr(got, want):
    want = complex(want)
    if got is None:
        return math.inf
    got = complex(got)
    if cmath.isnan(got) or cmath.isinf(got):
        return math.inf
    scale = max(abs(want), 1e-300)
    return abs(got - want) / scale


def run(verbose=True):
    with mp.workdps(DPS):
        return _run(verbose)


def _run(verbose):
    bad = []
    report = []
    for label, f, ref, pts in CASES:
        worst = 0.0
        worst_pt = None
        for p in pts:
            try:
                want = ref(*p)
            except (ValueError, ZeroDivisionError, TypeError):
                continue
            w = complex(want)
            if cmath.isinf(w) or cmath.isnan(w):
                continue
            try:
                got = f(*p)
            except Exception:  # noqa: BLE001
                got = None
            e = relerr(got, w)
            # near zeros of the function compare absolutely
            if abs(w) < 1e-6 and got is not None:
                e = min(e, abs(complex(got) - w))
            if e > worst:
                worst, worst_pt = e, (p, got, w)
            if e > REL:
                bad.append((label, p, got, w, e))
        report.append((label, len(pts), worst, worst_pt))
    if verbose:
        for label, n, worst, wp in report:
            flag = "" if worst <= 1e-14 else "   <-- " + repr(wp)
            print(f"{label:16s} {n:4d} pts  max rel err {worst:9.2e}{flag}")
        print(f"\n{len(bad)} failures (tol {REL})")
    return bad, report


def test_special_numeric():
    bad, _ = run(verbose=False)
    assert not bad, "\n".join(f"{b[0]}{b[1]}: got {b[2]} want {b[3]} (err {b[4]:.2e})" for b in bad[:30])


if __name__ == "__main__":
    b, _ = run()
    sys.exit(1 if b else 0)
