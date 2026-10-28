"""Native modular resultant, cyclotomic product formula and decomposition.

Expected values are hard-coded from SymPy 1.14 or follow from exact
certificates (closed forms of resultants, prod_{d | n} Phi_d = x**n - 1,
composition of the decomposition chain).
"""
from symox import (Poly, Rational, cyclotomic_poly, decompogen, discriminant,
                   divisors, expand, resultant, symbols)
from symox.solvers.decompogen import compogen

x, y = symbols('x y')
R = Rational


def test_resultant_values_match_the_reference():
    assert resultant(x**3/2 + x - R(1, 3), R(3, 5)*x**2 - 7, x) == R(-58823, 500)
    assert resultant(3*x**4 - 2*x + 5, 7*x**3 + x**2 - 4, x) == 260750
    assert resultant(x**2 - 1, x - 1, x) == 0
    assert resultant(2*x + 1, 4*x**2 - 1, x) == 0
    assert resultant(x**2 + 1, x - 1, x, polys=True) == 2
    assert Poly(x**5 + 3*x + 1, x).resultant(Poly(x**2 - 2, x)) == -97
    assert Poly(x**5/3 + 3*x + 1, x).resultant(Poly(x**2 - R(1, 2), x)) == R(-1081, 288)


def test_resultant_with_prs_is_unchanged():
    assert resultant(x**3 - 2*x + 1, 2*x**2 + 3, x, includePRS=True) == (
        155, [x**3 - 2*x + 1, 2*x**2 + 3, 4 - 14*x, 155])
    value, prs = Poly(x**5 + 3*x + 1, x).resultant(Poly(x**2 - 2, x), includePRS=True)
    assert value == -97
    assert prs == [Poly(x**5 + 3*x + 1, x), Poly(x**2 - 2, x), Poly(7*x + 1, x),
                   Poly(-97, x, domain='ZZ')]


def test_multivariate_resultant_keeps_the_prs_path():
    assert resultant(x**2 + y, x*y - 1, x) == y**3 + 1


def test_resultant_closed_form_with_large_coefficients():
    # res((x+1)**60, g) = g(-1)**60 since (x+1)**60 is monic.
    f, g = expand((x + 1)**60), expand((x - 2)**45 + 1)
    assert resultant(f, g, x) == (3**45 - 1)**60
    assert resultant(g, f, x) == (3**45 - 1)**60
    # Degrees of the reported witness, closed form ((-3)**150 + 1)**200.
    f, g = expand((x + 1)**200), expand((x - 2)**150 + 1)
    assert resultant(f, g, x) == (3**150 + 1)**200


def test_discriminants_match_the_reference():
    assert discriminant(x**2 + 3*x + 1) == 5
    assert discriminant(x**3/2 - x + R(1, 7)) == R(365, 196)
    assert discriminant(5*x**5 - 3*x**2 + 7) == 4684860425
    assert discriminant(4*x**4 - 1, x) == -16384
    assert discriminant(x**2 - 2*x + 1) == 0
    assert Poly(x**3 - 2, x).discriminant() == -108
    assert Poly(x**3/3 - 2, x).discriminant() == -12


def _coefficients(n):
    return Poly(cyclotomic_poly(n, x), x).all_coeffs()


def test_cyclotomic_product_certificate():
    for n in list(range(1, 80)) + [105, 210, 385, 1155]:
        product = Poly(1, x)
        for d in divisors(n):
            product *= Poly(cyclotomic_poly(d, x), x)
        assert product == Poly(x**n - 1, x)
    assert _coefficients(105) == [
        1, 1, 1, 0, 0, -1, -1, -2, -1, -1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1, 0,
        -1, 0, -1, 0, -1, 0, -1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -2, -1,
        -1, 0, 0, 1, 1, 1]


def _fingerprint(cs):
    return sum(int(v)*(i + 1) for i, v in enumerate(cs)) % (2**61 - 1)


def test_cyclotomic_orders_beyond_the_former_bound():
    cs = _coefficients(20001)
    assert len(cs) - 1 == 12992
    assert max(abs(v) for v in cs) == 2
    assert _fingerprint(cs) == 6497
    cs = _coefficients(30030)
    assert (len(cs) - 1, max(abs(v) for v in cs), _fingerprint(cs)) == (5760, 23, 2881)


def test_decompogen_beyond_the_former_degree_bound():
    assert decompogen(x**4100 + x**2, x) == [x**2050 + x, x**2]
    assert decompogen(x**6000 + 2*x**3 + 1, x) == [x**2000 + 2*x + 1, x**3]
    assert decompogen(x**8194 + x**2, x) == [x**4097 + x, x**2]
    chain = decompogen(expand((x**2 + x + 1)**3), x)
    assert chain == [x**3 + 3*x**2 + 3*x + 1, x**2 + x]
    assert expand(compogen(chain, x)) == expand((x**2 + x + 1)**3)
