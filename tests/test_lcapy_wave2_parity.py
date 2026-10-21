"""Regressions found by running Lcapy against Symox (wave 2).

The expected values are SymPy 1.14 outputs, hardcoded so the tests run
without SymPy or Lcapy.
"""
import symox as sm


def test_gcd_keeps_its_docstring():
    # Lcapy builds its own docstrings from sympy.gcd.__doc__ at import time.
    assert sm.gcd.__doc__


def test_cancel_and_simplify_pull_radical_factors():
    R1, R2, T, k = sm.symbols('R1 R2 T k_B', positive=True)
    e = sm.sqrt(4*R1*T*k + 4*R2*T*k)
    expected = '2*sqrt(T)*sqrt(k_B)*sqrt(R1 + R2)'
    assert str(sm.cancel(e)) == expected
    assert str(sm.factor_terms(e)) == expected
    assert str(sm.simplify(e)) == expected
    assert sm.simplify(e - 2*sm.sqrt(T)*sm.sqrt(k)*sm.sqrt(R1 + R2)) == 0


def test_simplify_keeps_shorter_expansion_after_factor_terms():
    n = sm.Symbol('n')
    assert str(sm.simplify(n**3 - n)) == 'n**3 - n'
    assert str(sm.simplify(n*(n**2 - 1))) == 'n**3 - n'


def test_simplify_applies_exptrigsimp():
    e = sm.exp(2*sm.I)/(1 - sm.exp(4*sm.I))
    assert str(sm.simplify(sm.im(e))) == '1/(2*sin(2))'
    assert sm.simplify(sm.re(e)) == 0
    assert str(sm.simplify(e)) == 'I/(2*sin(2))'
    x = sm.Symbol('x', real=True)
    assert str(sm.simplify(sm.im(sm.exp(sm.I*x)/(1 - sm.exp(2*sm.I*x))))) == '1/(2*sin(x))'


def test_unevaluated_product_numbers_are_coefficients():
    x = sm.Symbol('x')
    p = sm.Poly(sm.Mul(-1, 2, x, evaluate=False), x)
    assert repr(p) == "Poly(-2*x, x, domain='ZZ')"
    p = sm.Poly(sm.Mul(-1, 2, x, evaluate=False))
    assert repr(p) == "Poly(-2*x, x, domain='ZZ')"
    p = sm.Poly(sm.Mul(2, 3, x, evaluate=False), x)
    assert repr(p) == "Poly(6*x, x, domain='ZZ')"


def test_gaussian_fraction_coefficients_build_a_polynomial():
    # SymPy selects ZZ_I(N,pi); Symox has no ZZ_I fraction-field backend yet
    # and keeps the coefficients exact in EX (see docs/algorithmic-parity.md)
    # instead of raising.
    n, N = sm.symbols('n N')
    e = sm.Mul(-1, 2, sm.I, sm.pi, 1/N, n, evaluate=False)
    p = e.as_poly(n)
    assert p.is_linear
    assert p.gens == (n,)
    assert p.as_expr() == -2*sm.I*sm.pi*n/N
    assert p.coeffs() == [-2*sm.I*sm.pi/N]
    x, y = sm.symbols('x y')
    q = sm.Poly(sm.I*sm.pi*x/y, x)
    assert q.as_expr() == sm.I*sm.pi*x/y
    assert (q**2).as_expr() == -sm.pi**2*x**2/y**2
