"""Polynomials over RR[y], CC[y] and RR(y) (SymPy 1.14 reference outputs)."""
import pytest

import symox
from symox.polys.domains import RR, CC

NAMESPACE = {name: getattr(symox, name) for name in dir(symox) if not name.startswith('_')}
NAMESPACE.update(RR=RR, CC=CC)
NAMESPACE.update(zip(('x', 'y', 'z'), symox.symbols('x y z')))


def _run(source):
    try:
        return repr(eval(source, dict(NAMESPACE)))
    except Exception as error:  # the reference error class and message
        return 'EXC %s %s' % (type(error).__name__, str(error)[:150])


CASES = [
    ('Poly(1.5*x**2*y-1.5*y,x).factor_list()', "EXC TypeError unsupported operand type(s) for *: 'PolyElement' and 'PolyElement'"),
    ('Poly(2.0*y,x).factor_list()', "EXC TypeError unsupported operand type(s) for *: 'PolyElement' and 'PolyElement'"),
    ('Poly(1.5*x*y,x).factor_list()', "EXC TypeError unsupported operand type(s) for *: 'PolyElement' and 'PolyElement'"),
    ("Poly(1.5*x**2-1.5,x,domain='RR[y]').factor_list()", "EXC TypeError unsupported operand type(s) for *: 'PolyElement' and 'PolyElement'"),
    ('Poly(x**2*y/2.0-1/y,x).factor_list()', "EXC TypeError bad operand type for abs(): 'FracElement'"),
    ('Poly(I*x**2*y-1.5*y,x).factor_list()', "EXC TypeError unsupported operand type(s) for *: 'PolyElement' and 'PolyElement'"),
    ('Poly(1.5*x**2*y-3*x*y+1.5*y,x).sqf_list()', "(1.0*y, [(Poly(0.5*x - 0.5, x, domain='RR[y]'), 2)])"),
    ('Poly(1.5*x**2*y-3*x*y+1.5*y,x).sqf_part()', "Poly(0.5*x - 0.5, x, domain='RR[y]')"),
    ('sqf(1.5*x**2*y-3*x*y+1.5*y,x)', '1.5*y*(1.0*x - 1.0)**2'),
    ('Poly(1.5*x**2*y-3*x*y+1.5*y,x).gcd(Poly(x*y-y,x))', "Poly(1.0*y*x - 1.0*y, x, domain='RR[y]')"),
    ('Poly(1.5*x**2*y-3*x*y+1.5*y,x).cancel(Poly(x*y-y,x))', "(1.00000000000000, Poly(1.5*x - 1.5, x, domain='RR[y]'), Poly(1.0, x, domain='RR[y]'))"),
    ('cancel((1.5*x**2*y-3*x*y+1.5*y)/(x*y-y),x)', '1.5*x - 1.5'),
    ('cancel((1.5*x*y+3*y)/(x*y+2*y),x)', '1.50000000000000'),
    ('cancel(1.5*x/(x*y),x)', '1.5/y'),
    ('Poly(1, x).cancel(Poly(51.05*x*y - 1.0, x))', "(1.00000000000000, Poly(1.0, x, domain='RR[y]'), Poly(51.05*y*x - 1.0, x, domain='RR[y]'))"),
    ('Poly(x**2*y-1.0,x).gcd(Poly(x*y-1.0,x))', "Poly(1.0, x, domain='RR[y]')"),
    ('Poly(0.1*x**2*y+0.3*x*y+0.2*y,x).gcd(Poly(0.1*x*y+0.1*y,x))', "Poly(1.0*y, x, domain='RR[y]')"),
    ('Poly(0.1*x**2*y+0.3*x*y+0.2*y,x).cofactors(Poly(0.1*x*y+0.1*y,x))', "(Poly(1.0*y, x, domain='RR[y]'), Poly(0.1*x**2 + 0.3*x + 0.2, x, domain='RR[y]'), Poly(0.1*x + 0.1, x, domain='RR[y]'))"),
    ('Poly(0.1*x**2*y+0.3*x*y+0.2*y,x).div(Poly(0.1*x*y+0.1*y,x))', "(Poly(1.0*x + 2.0, x, domain='RR[y]'), Poly(8.67361737988403e-17*y, x, domain='RR[y]'))"),
    ('Poly(0.1*x**2+0.3*x*y+0.2*y**2,x).div(Poly(0.3*x+0.1*y,x))', "(Poly(0.333333333333333*x + 0.888888888888889*y, x, domain='RR[y]'), Poly(0.111111111111111*y**2, x, domain='RR[y]'))"),
    ('Poly(x**2+y,x,domain=RR[y]).div(Poly(2*x*y+1,x,domain=RR[y]))', "(Poly(0.5/y*x - 0.25/y**2, x, domain='RR(y)'), Poly(0.25*(4.0*y**3 + 1.0)/y**2, x, domain='RR(y)'))"),
    ("Poly(x**2*y/2.0+1/y,x).div(Poly(x/y+1,x,domain='RR(y)'))", 'EXC OptionError expected a valid domain specification, got RR(y)'),
    ('Poly(x**2*y/2.0+1/y,x).gcd(Poly(x*y+1.0,x))', "Poly(1.0, x, domain='RR(y)')"),
    ('Poly(x**2*y/2.0+1/y,x).monic()', "Poly(1.0*x**2 + 2.0/y**2, x, domain='RR(y)')"),
    ('Poly(x**2*y/2.0+1/y,x).cancel(Poly(x*y+1.0,x))', "(1.0/y, Poly(0.5*y**2*x**2 + 1.0, x, domain='RR(y)'), Poly(1.0*y*x + 1.0, x, domain='RR(y)'))"),
    ('Poly(1.5*x*y+2*y,x).lcm(Poly(3*x+4,x))', "Poly(1.5*y*x + 2.0*y, x, domain='RR[y]')"),
    ('Poly(I*x*y+2*y,x).gcd(Poly(y*x-2*I*y,x))', "Poly(y*x - 2*I*y, x, domain='ZZ_I[y]')"),
    ('Poly(I*x*y+2*y,x).cancel(Poly(y*x-2*I*y,x))', "(1, Poly(I, x, domain='ZZ_I[y]'), Poly(1, x, domain='ZZ_I[y]'))"),
    ('Poly(-1.5*x*y-2,x).cancel(Poly(-y*x,x))', "(1.00000000000000, Poly(1.5*y*x + 2.0, x, domain='RR[y]'), Poly(1.0*y*x, x, domain='RR[y]'))"),
    ('Poly(x*y+2.0,x).content()', '1.00000000000000'),
    ('Poly(x*y-2.0*y,x).primitive()', "(1.0*y, Poly(1.0*x - 2.0, x, domain='RR[y]'))"),
    ('Poly(-x*y-2.0*y,x).primitive()', "(1.0*y, Poly(-1.0*x - 2.0, x, domain='RR[y]'))"),
    ('Poly(I*x*y-2.0*y,x).primitive()', "(1.0*y, Poly(1.0*I*x - 2.0, x, domain='CC[y]'))"),
    ('Poly(1.5*x*y+2,x).terms_gcd()', "((0,), Poly(1.5*y*x + 2.0, x, domain='RR[y]'))"),
    ('Poly(1.5*x**2*y+2*x,x).terms_gcd()', "((1,), Poly(1.5*y*x + 2.0, x, domain='RR[y]'))"),
    ('roots(Poly(1.5*x**2*y-1.5*y,x))', '{-1.00000000000000: 1, 1.00000000000000: 1}'),
    ('Poly(1.5*x**2*y-1.5*y,x).nroots()', 'EXC DomainError Numerical domain expected, got RR[y]'),
    ('solve(1.5*x**2*y-1.5*y,x)', '[-1.00000000000000, 1.00000000000000]'),
    ('together(1.5/x+2/y)', '(2*x + 1.5*y)/(x*y)'),
    ('Poly(2.0*x*y**2+4.0*y,x).clear_denoms()', "(1, Poly(2.0*y**2*x + 4.0*y, x, domain='RR[y]'))"),
    ('Poly(x/(2.0*y)+1,x).clear_denoms()', "(1.0*y, Poly(0.5*x + 1.0*y, x, domain='RR(y)'))"),
    ('Poly(1.5*x*y+2,x).as_expr()', '1.5*x*y + 2.0'),
    ('Poly(1.5*x*y+2,x).eval({x:1})', '1.5*y + 2.0'),
    ('Poly(1.5*x*y+2,x).subs(y,2)', '3.0*x + 2.0'),
    ('Poly(1.5*x*y+2,x)==Poly(1.5*x*y+2,x)', 'True'),
    ('Poly(1.5*x*y+2,x).is_monic', 'False'),
    ('Poly(1.5*x*y+2,x).sqf_norm()', 'EXC DomainError ground domain must be algebraic'),
    ('Poly(1.5*x*y+2,x).integrate()', "Poly(0.75*y*x**2 + 2.0*x, x, domain='RR(y)')"),
    ('Poly(1.5*x*y+2,x).compose(Poly(x+y,x))', "Poly(1.5*y*x + 1.5*y**2 + 2.0, x, domain='RR[y]')"),
    ('Poly(1.5*x*y+2,x)**2', "Poly(2.25*y**2*x**2 + 6.0*y*x + 4.0, x, domain='RR[y]')"),
    ('Poly(1.5*x*y+2,x).pdiv(Poly(y*x+1,x))', "(Poly(1.5*y, x, domain='RR[y]'), Poly(0.5*y, x, domain='RR[y]'))"),
    ('Poly(1.5*x*y+2,x).prem(Poly(3*x+1,x))', "Poly(6.0 - 1.5*y, x, domain='RR[y]')"),
    ('Poly(1.5*x*y+2,x).resultant(Poly(3*x+1,x))', '1.5*y - 6.0'),
    ('Poly(1.5*x*y+2,x).discriminant()', '1.00000000000000'),
    ('Poly(1.5*x**2*y+2,x).discriminant()', '-12.0*y'),
    ('Poly(1.5*x*y+2,x).half_gcdex(Poly(3*x+1,x))', "(Poly(-3.0/(1.5*y - 6.0), x, domain='RR(y)'), Poly(1.0, x, domain='RR(y)'))"),
    ('Poly(1.5*x*y+2,x).to_ring()', "Poly(1.5*y*x + 2.0, x, domain='RR[y]')"),
    ('Poly(1.5*x*y+2,x).abs()', "EXC TypeError bad operand type for abs(): 'PolyElement'"),
    ('-Poly(1.5*x*y+2,x)', "Poly(-1.5*y*x - 2.0, x, domain='RR[y]')"),
    ('Poly(1.5*x*y+2,x).max_norm()', "EXC TypeError bad operand type for abs(): 'PolyElement'"),
    ('Poly(1.5*x*y+2,x).l1_norm()', "EXC TypeError bad operand type for abs(): 'PolyElement'"),
]


@pytest.mark.parametrize("source, expected", CASES)
def test_inexact_composite_reference(source, expected):
    assert _run(source) == expected
