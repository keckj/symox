"""RR and CC polynomial domains follow the floating-point reference results.

Every expected string below is the output of the reference SymPy 1.14.
"""
import pytest
import symox as sx
from symox import Poly, RR, CC, symbols

x, y = symbols('x y')


def test_frac_field_and_poly_ring_attributes():
    assert str(RR.frac_field(x).field) == 'Rational function field in x over RR with lex order'
    assert str(CC.frac_field(x).field) == 'Rational function field in x over CC with lex order'
    ring = RR.poly_ring(x).ring
    assert str(ring) == 'Polynomial ring in x over RR with lex order'
    assert ring.domain == RR
    assert ring.ngens == 1
    assert ring.symbols == (x,)
    with pytest.raises(AttributeError):
        RR.ring


@pytest.mark.parametrize('call,expected', [
    (lambda: Poly(1.5*x - 1.5*y, x, y).content(), '1.00000000000000'),
    (lambda: Poly(3.0*x*y + 6.0*y, x, y).primitive(),
     "(1.00000000000000, Poly(3.0*x*y + 6.0*y, x, y, domain='RR'))"),
    (lambda: sx.primitive(3.0*x*y + 6.0*y), '(1.00000000000000, 3.0*x*y + 6.0*y)'),
    (lambda: sx.content(3.0*x*y + 6.0*y, x, y), '1.00000000000000'),
])
def test_rr_content_is_one(call, expected):
    assert str(call()) == expected


@pytest.mark.parametrize('call,expected', [
    (lambda: Poly(1.5*x, x).invert(Poly(x**2 + 1, x)),
     "Poly(-0.666666666666667*x, x, domain='RR')"),
    (lambda: Poly(x + 0.5, x).invert(Poly(x**2 + 1.0, x)),
     "Poly(-0.8*x + 0.4, x, domain='RR')"),
    (lambda: Poly(1.5*x + 1, x).gcdex(Poly(x**2 + 1, x)),
     "(Poly(-0.461538461538462*x + 0.307692307692308, x, domain='RR'), "
     "Poly(0.692307692307692, x, domain='RR'), Poly(1.0, x, domain='RR'))"),
    (lambda: Poly(1.5*x + 1, x).half_gcdex(Poly(x**2 + 1, x)),
     "(Poly(-0.461538461538462*x + 0.307692307692308, x, domain='RR'), "
     "Poly(1.0, x, domain='RR'))"),
    (lambda: sx.invert(1.5*x, x**2 + 1), '-0.666666666666667*x'),
    (lambda: sx.invert(x + 0.5, x**2 + 1), '0.4 - 0.8*x'),
    (lambda: sx.invert(2.0*x + 1, x**2 - 2), '0.285714285714286*x - 0.142857142857143'),
    (lambda: sx.gcdex(1.5*x + 1, x**2 + 1),
     '(0.307692307692308 - 0.461538461538462*x, 0.692307692307692, 1.00000000000000)'),
    (lambda: sx.gcdex(x + 0.5, x**2 - 1),
     '(1.33333333333333*x - 0.666666666666667, -1.33333333333333, 1.00000000000000)'),
    (lambda: sx.half_gcdex(1.5*x + 1, x**2 + 1),
     '(0.307692307692308 - 0.461538461538462*x, 1.00000000000000)'),
])
def test_rr_extended_euclid_stays_in_floats(call, expected):
    assert str(call()) == expected


def test_rr_invert_of_zero_divisor_raises():
    with pytest.raises(sx.polys.polyerrors.NotInvertible, match='zero divisor'):
        sx.invert(x**2 - 1.0, x - 1)


@pytest.mark.parametrize('call,expected', [
    (lambda: Poly(x**2 - 1.0, x).cancel(Poly(x - 1, x)),
     "(1.00000000000000, Poly(1.0*x + 1.0, x, domain='RR'), Poly(1.0, x, domain='RR'))"),
    (lambda: Poly(1.5*x**2 - 1.5, x).cancel(Poly(x + 1, x)),
     "(1.00000000000000, Poly(1.5*x - 1.5, x, domain='RR'), Poly(1.0, x, domain='RR'))"),
    (lambda: Poly(2.0*x**2 - 2.0, x).cancel(Poly(2*x + 2, x)),
     "(1.00000000000000, Poly(2.0*x - 2.0, x, domain='RR'), Poly(2.0, x, domain='RR'))"),
    (lambda: Poly(1.5*x**2 - 1.5, x).cancel(Poly(3*x + 3, x), include=True),
     "(Poly(1.5*x - 1.5, x, domain='RR'), Poly(3.0, x, domain='RR'))"),
    (lambda: Poly(2.0*x**2 - 2.0, x).cancel(Poly(0.5*x + 0.5, x)),
     "(1.00000000000000, Poly(2.0*x - 2.0, x, domain='RR'), Poly(0.5, x, domain='RR'))"),
])
def test_rr_poly_cancel_keeps_floats(call, expected):
    assert str(call()) == expected
