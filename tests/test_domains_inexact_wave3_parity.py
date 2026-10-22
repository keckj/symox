"""Inexact and composite domains, wave 3 (SymPy 1.14 reference outputs).

``construct_domain`` over RealField/ComplexField grounds, ``CC.invert``,
element display (``2``, ``MPQ(1,2)``, ``mpf``, ``PolyElement`` over RR/CC),
unification with inexact grounds and ``parallel_poly_from_expr``.
"""
import pytest

import symox
from symox.polys.domains import RR, CC, ZZ, QQ, EX, RealField, ComplexField
from symox.polys.constructor import construct_domain
from symox.polys.polytools import parallel_poly_from_expr

NAMESPACE = {name: getattr(symox, name) for name in dir(symox) if not name.startswith('_')}
NAMESPACE.update(RR=RR, CC=CC, ZZ=ZZ, QQ=QQ, EX=EX, RealField=RealField, ComplexField=ComplexField,
                 construct_domain=construct_domain, parallel_poly_from_expr=parallel_poly_from_expr)
NAMESPACE.update(zip(('x', 'y', 'z', 't'), symox.symbols('x y z t')))


def _run(source):
    try:
        return repr(eval(source, dict(NAMESPACE)))
    except Exception as error:  # the reference error class and message
        return 'EXC %s %s' % (type(error).__name__, str(error)[:150])


CASES = [
    ('Poly(1.5*x*z + y, x, z)', "Poly(1.5*x*z + 1.0*y, x, z, domain='RR[y]')"),
    ('Poly(1.5*x*z + y*x, x, z)', "Poly(1.5*x*z + 1.0*y*x, x, z, domain='RR[y]')"),
    ('Poly(1.5*x + y*z, x, z)', "Poly(1.5*x + 1.0*y*z, x, z, domain='RR[y]')"),
    ('Poly(1.5*x*z + y*t, x, z)', "Poly(1.5*x*z + 1.0*t*y, x, z, domain='RR[y,t]')"),
    ('Poly(x/2.0 + y/z, x)', "Poly(0.5*x + 1.0*y/z, x, domain='RR(y,z)')"),
    ('Poly(1.5*x*z + y, x, z).as_expr()', '1.5*x*z + 1.0*y'),
    ('Poly(1.5*x*z + y, x, z).LC()', '1.50000000000000'),
    ('Poly(1.5*x*z + y, x, z).coeffs()', '[1.50000000000000, 1.0*y]'),
    ('Poly(1.5*x*z + y, x, z).eval(x, 1)', "Poly(1.5*z + 1.0*y, z, domain='RR[y]')"),
    ('Poly(2.5*x**2*z + y*x, x, z).degree(z)', '1'),
    ('Poly(1.5*x*z + y, z, x)', "Poly(1.5*z*x + 1.0*y, z, x, domain='RR[y]')"),
    ('Poly(1.5*x*y*z + 1, x, z)', "Poly(1.5*y*x*z + 1.0, x, z, domain='RR[y]')"),
    ('Poly(1.5*x*z + y, x, z).div(Poly(x, x, z))', "(Poly(1.5*z, x, z, domain='RR[y]'), Poly(1.0*y, x, z, domain='RR[y]'))"),
    ('Poly(1.5*x*z + sqrt(2)*y, x, z)', "Poly(1.5*x*z + sqrt(2)*y, x, z, domain='EX')"),
    ('construct_domain([1.5*x])', '(RR[x], [1.5*x])'),
    ('construct_domain([1.5*x, y])', '(RR[x,y], [1.5*x, y])'),
    ('construct_domain([1.5*x], composite=True)', '(RR[x], [1.5*x])'),
    ('construct_domain([1.5*x/y])', '(RR(x,y), [1.5*x/y])'),
    ('construct_domain([1.5*x/y], field=True)', '(RR(x,y), [1.5*x/y])'),
    ('construct_domain([1.5*x], field=True)', '(RR(x), [1.5*x])'),
    ('construct_domain([1.5*I*x])', '(CC[x], [(0.0 + 1.5j)*x])'),
    ('construct_domain([1.5*x + I])', '(CC[x], [(1.5 + 0.0j)*x + (0.0 + 1.0j)])'),
    ('construct_domain([1.5, x])', '(RR[x], [1.5, x])'),
    ('construct_domain([1.5, x], composite=False)', '(EX, [EX(1.50000000000000), EX(x)])'),
    ('construct_domain([2*x], composite=True)', '(ZZ[x], [2*x])'),
    ('construct_domain([x, 2.0])', '(RR[x], [x, 2.0])'),
    ('construct_domain([x*y, 2.0])', '(RR[x,y], [x*y, 2.0])'),
    ('construct_domain([1.5*x, sqrt(2)])', '(EX, [EX(1.5*x), EX(sqrt(2))])'),
    ('construct_domain([1.5*x, pi])', '(RR[x,pi], [1.5*x, pi])'),
    ('construct_domain([1/(x+1.5)])', '(RR(x), [1.0/(x + 1.5)])'),
    ('construct_domain([RealField(100)(1.5)*x])', '(RR[x], [1.5*x])'),
    ('construct_domain([Float(1.5, 30)*x])', '(RR[x], [1.5*x])'),
    ('construct_domain([Float(1.5, 30)*x, Float(2.5, 20)])', '(RR[x], [1.5*x, 2.5])'),
    ('construct_domain({(0,): 1.5*x, (1,): 2})', '(RR[x], {(0,): 1.5*x, (1,): 2.0})'),
    ('CC.invert(CC(1), CC(2))', "EXC TypeError unsupported operand type(s) for %: 'mpc' and 'mpc'"),
    ('RR.invert(RR(1), RR(2))', "mpf('1.0')"),
    ('CC.invert(CC(1,1), CC(2))', "EXC TypeError unsupported operand type(s) for %: 'mpc' and 'mpc'"),
    ('ZZ(2)', '2'),
    ('QQ(2)', 'MPQ(2,1)'),
    ('QQ(1,2)', 'MPQ(1,2)'),
    ('ZZ[x](2)', '2'),
    ('ZZ[x](x)', 'x'),
    ('QQ[x](2)', '2'),
    ('ZZ.frac_field(x)(2)', '2'),
    ('repr(ZZ(2))', "'2'"),
    ('str(ZZ(2))', "'2'"),
    ('repr(QQ(1,2))', "'MPQ(1,2)'"),
    ('ZZ.convert(2)', '2'),
    ('EX(2)', 'EX(2)'),
    ('RR(2)', "mpf('2.0')"),
    ('CC(2)', "mpc(real='2.0', imag='0.0')"),
    ('CC(1,2)', "mpc(real='1.0', imag='2.0')"),
    ('parallel_poly_from_expr([x*y + sqrt(2), x + y], x)', "([Poly(y*x + sqrt(2), x, domain='EX'), Poly(x + y, x, domain='EX')], {'gens': (x,), 'domain': EX, 'polys': False})"),
    ('parallel_poly_from_expr([x*y + sin(y), x + y], x)', "([Poly(y*x + sin(y), x, domain='EX'), Poly(x + y, x, domain='EX')], {'gens': (x,), 'domain': EX, 'polys': False})"),
    ('parallel_poly_from_expr([x*y, x + 1.5], x)', "([Poly(1.0*y*x, x, domain='RR[y]'), Poly(1.0*x + 1.5, x, domain='RR[y]')], {'gens': (x,), 'domain': RR[y], 'polys': False})"),
    ('parallel_poly_from_expr([x*y + 1.5*z, x + y], x)', "([Poly(1.0*y*x + 1.5*z, x, domain='RR[y,z]'), Poly(1.0*x + 1.0*y, x, domain='RR[y,z]')], {'gens': (x,), 'domain': RR[y,z], 'polys': False})"),
    ('parallel_poly_from_expr([x*y, x + I], x)', "([Poly(y*x, x, domain='ZZ_I[y]'), Poly(x + I, x, domain='ZZ_I[y]')], {'gens': (x,), 'domain': ZZ_I[y], 'polys': False})"),
    ('parallel_poly_from_expr([x + exp(y), x*y], x)', "([Poly(x + exp(y), x, domain='EX'), Poly(y*x, x, domain='EX')], {'gens': (x,), 'domain': EX, 'polys': False})"),
    ('parallel_poly_from_expr([x + 1.5, x*y], x, z)', "([Poly(1.0*x + 1.5, x, z, domain='RR[y]'), Poly(1.0*y*x, x, z, domain='RR[y]')], {'gens': (x, z), 'domain': RR[y], 'polys': False})"),
    ('RR[x]', 'RR[x]'),
    ('RR[x, y]', 'RR[x,y]'),
    ('RR.frac_field(x)', 'RR(x)'),
    ('RR.frac_field(x, y)', 'RR(x,y)'),
    ('RealField(100)[x]', 'RR[x]'),
    ('RealField(100)[x](1.5)', '1.5'),
    ('RealField(100).frac_field(x)', 'RR(x)'),
    ('CC[x]', 'CC[x]'),
    ('CC[x](1.5)', '(1.5 + 0.0j)'),
    ('CC[x](x + 1)', 'x + (1.0 + 0.0j)'),
    ('CC[x](I*x + 1)', '(0.0 + 1.0j)*x + (1.0 + 0.0j)'),
    ('CC[x](x)*CC[x](x+I)', 'x**2 + (0.0 + 1.0j)*x'),
    ('RR[x](x/2)', '0.5*x'),
    ('RR[x](x)**2', 'x**2'),
    ('RR[x](1.5*x + 2)', '1.5*x + 2.0'),
    ('RR[x].gens', '(x,)'),
    ('RR[x].domain', 'RR'),
    ('RR[x].dom', 'RR'),
    ('RR[x].is_Exact', 'False'),
    ('RR[x].get_field()', 'RR(x)'),
    ('RR[x].get_exact()', 'QQ[x]'),
    ('RR.frac_field(x)(1/x)', '1.0/x'),
    ('RR.frac_field(x)(x/2)', '0.5*x'),
    ('RR.frac_field(x).get_ring()', 'RR[x]'),
    ('RR[x] == RR[x]', 'True'),
    ('RR[x] == RealField(100)[x]', 'False'),
    ('CC[x].is_ComplexField', 'False'),
    ('RealField(100)[x].domain.precision', '100'),
    ('RR[x].convert(1.5)', '1.5'),
    ('RR[x].convert(ZZ(2), ZZ)', '2.0'),
    ('CC[x].from_sympy(I*x + 1.5)', '(0.0 + 1.0j)*x + (1.5 + 0.0j)'),
    ('CC[x].to_sympy(CC[x](x + 1))', '1.0*x + 1.0'),
    ('RR[x].unify(QQ)', 'RR[x]'),
    ('RR[x].unify(CC)', 'CC[x]'),
    ('RR[x].unify(QQ[y])', 'RR[x,y]'),
    ('RR.frac_field(x).unify(ZZ[y])', 'RR(x,y)'),
    ('ZZ[x].unify(RR)', 'RR[x]'),
    ('RR[x].characteristic()', '0'),
    ('RR[x].zero', '0.0'),
    ('CC[x].one', '(1.0 + 0.0j)'),
    ('RR[x](1.5) == 1.5', 'True'),
    ('RR[x].is_PolynomialRing', 'True'),
    ('str(RR[x])', "'RR[x]'"),
    ('str(CC.frac_field(x))', "'CC(x)'"),
    ('RR[x].ngens', '1'),
    ('RR[x].symbols', '(x,)'),
]


@pytest.mark.parametrize("source, expected", CASES)
def test_domains_inexact_wave3_reference(source, expected):
    assert _run(source) == expected
