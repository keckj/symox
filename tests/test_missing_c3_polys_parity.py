"""Polys parity for apart(full=True), partfrac lists, Poly domains, Poly.per,
CRootOf.eval_approx, which_all_roots, nested RootSum.doit and sfield.

Expected values are the SymPy 1.14 outputs."""

from symox import (
    CRootOf, Float, I, Lambda, Poly, QQ, RR, Rational, RootSum, Symbol,
    exp, log, sqrt, symbols,
)
from symox.polys.fields import field, sfield
from symox.polys.partfrac import apart, apart_list, assemble_partfrac_list
from symox.polys.polyclasses import DMP
from symox.polys.rings import ring, sring

x, y, t = symbols('x y t')


def test_apart_full_rootsum():
    f = 1/(x**5 + 1)
    assert str(apart(f, full=False)) == (
        "-(x**3 - 2*x**2 + 3*x - 4)/(5*(x**4 - x**3 + x**2 - x + 1)) + 1/(5*(x + 1))")
    full = apart(f, full=True)
    assert str(full) == (
        "-RootSum(_w**4 - _w**3 + _w**2 - _w + 1, Lambda(_a, _a/(-_a + x)))/5 + 1/(5*(x + 1))")
    assert str(full.doit()) == str(apart(f, full=True).doit())
    assert str(apart(1/(x**2 + 1), x, full=True)) == (
        "-RootSum(_w**2 + 1, Lambda(_a, _a/(-_a + x)))/2")


def test_apart_full_rational_factors():
    f = (x**2 + 1)/(x*(x - 1)**2)
    assert apart(f, x, full=True) == apart(f, x)
    assert str(apart(y/(x**2 - 2), x, full=True)) == (
        "y*RootSum(_w**2 - 2, Lambda(_a, _a/(-_a + x)))/4")


def test_apart_list_and_assemble():
    f = (2*x**3 - 2*x)/(x**2 - 2*x + 1)
    pfd = apart_list(f)
    assert str(pfd) == "(1, Poly(2*x + 4, x, domain='ZZ'), [(Poly(_w - 1, _w, domain='ZZ'), Lambda(_a, 4), Lambda(_a, -_a + x), 1)])"
    assert assemble_partfrac_list(pfd) == 2*x + 4 + 4/(x - 1)
    pfd = apart_list(-2/(x**2 - 2))
    assert str(assemble_partfrac_list(pfd)) == "-RootSum(_w**2 - 2, Lambda(_a, _a/(-_a + x)))/2"
    assert str(assemble_partfrac_list(pfd).doit()) == "sqrt(2)/(2*(x + sqrt(2))) - sqrt(2)/(2*(x - sqrt(2)))"


def test_rootsum_doit_inside_expression():
    a = Symbol('a')
    r = RootSum(Poly(x**2 - 2, x), Lambda(a, 1/(y - a)), auto=False)
    assert str((r/2).doit()) == "1/(2*(y + sqrt(2))) + 1/(2*(y - sqrt(2)))"


def test_poly_real_domain_aliases():
    p = Poly(x + sqrt(2), x, domain='R')
    assert str(p) == "Poly(1.0*x + 1.4142135623731, x, domain='RR')"
    assert Poly(x/3, x, domain='RR_53').domain == RR
    assert str(Poly(x/2 + 1, x, domain='Q')) == "Poly(1/2*x + 1, x, domain='QQ')"
    assert Poly(2*x, x, domain='Z').domain.is_ZZ


def test_poly_per_with_dmp():
    p = Poly(x**2 + 1, x)
    rep = DMP([QQ(1), QQ(0), QQ(1)], QQ)
    assert str(p.per(rep)) == "Poly(x**2 + 1, x, domain='QQ')"
    q = Poly(x*y + 1, x, y)
    assert str(q.per(DMP([QQ(1), QQ(2)], QQ), remove=0)) == "Poly(y + 2, y, domain='QQ')"


def test_crootof_eval_approx():
    r = CRootOf(x**3 + x + 3, 0)
    assert str(r.eval_approx(2)) == "-1.2"
    s = CRootOf(x**2 + 1, 0)
    assert str(s.eval_approx(2)) == "-1.0*I"
    assert s.eval_approx(2).as_real_imag()[0] == 0


def test_which_all_roots():
    p = Poly(x**4 - 1)
    assert p.which_all_roots([-1, 1, -I, I, 0]) == [-1, 1, -I, I]
    assert p.which_all_roots([-1, 1, -I, I, 0, 0]) == [-1, 1, -I, I]
    q = Poly(x**3 + x + 3)
    roots = q.all_roots()
    assert q.which_all_roots(roots + [Rational(1, 3)]) == roots


def test_sfield_infers_function_generators():
    K, f = sfield((x*log(x) + 4*x**2)*exp(1/x + log(x)/3)/x**2)
    assert str(K) == ("Rational function field in x, exp(1/x), log(x), x**(1/3) "
                      "over ZZ with lex order")
    assert str(f) == ("(4*x**2*(exp(1/x)) + x*(exp(1/x))*(log(x)))"
                      "/((x**(1/3))**5)")
    K, f = sfield(x/y + 1)
    assert str(f) == "(x + y)/y"
    R, p = sring(x + exp(x))
    assert str(R) == "Polynomial ring in x, exp(x) over ZZ with lex order"
    assert str(p) == "x + (exp(x))"


def test_ring_element_str():
    R, a, b = ring('a,b', QQ)
    assert str(a/2 + b**2) == "1/2*a + b**2"
    assert str(Rational(-3, 2)*a**2) == "-3/2*a**2"
    assert str(-a*b - Rational(1, 3)) == "-a*b - 1/3"
    R, a, b = ring('a,b', RR)
    assert str(a*2.5 - 1/3.) == "2.5*a - 0.333333333333333"
    assert str(R(0)) == "0.0"
    Ry, s = ring('t', QQ[y])
    assert str((y + 1)*s**2 - y*s + 3) == "(y + 1)*t**2 - y*t + 3"
    K, u, v = field('u,v', QQ)
    assert str(Rational(1, 2)*u/v) == "u/(2*v)"
    assert str((u + 1)/(u - v)) == "(u + 1)/(u - v)"
    R, p = sring(x + 2.5*y)
    assert str(p) == "x + 2.5*y"
    R, p = sring(x + sqrt(2), extension=True)
    assert str(R) == "Polynomial ring in x over QQ<sqrt(2)> with lex order"
