"""SymPy 1.14 parity for Poly methods, polytools helpers, polyfuncs,
special and orthogonal polynomials, partial fraction lists, monomials,
orderings and the polyerrors hierarchy.

Every case is evaluated with ``sympy`` and ``symox`` namespaces and the
printed results (or exception class and message) are compared.
"""
import pickle
import re

import pytest
import sympy
import symox

from symox.polys.monomials import Monomial as XMonomial
from sympy.polys.monomials import Monomial as SMonomial


def _namespace(lib):
    ns = {name: getattr(lib, name) for name in dir(lib) if not name.startswith('_')}
    if lib is sympy:
        from sympy.polys.orthopolys import gegenbauer_poly, spherical_bessel_fn
    else:
        from symox.polys.orthopolys import gegenbauer_poly, spherical_bessel_fn
    ns['gegenbauer_poly'] = gegenbauer_poly
    ns['spherical_bessel_fn'] = spherical_bessel_fn
    ns['M'] = lib
    ns['pickle'] = pickle
    x, y, z, a, b = lib.symbols('x y z a b')
    ns.update(x=x, y=y, z=z, a=a, b=b, p=lib.symbols('p', positive=True))
    return ns


NS = {'s': _namespace(sympy), 'x': _namespace(symox)}
_DUMMY = re.compile(r'_(x|w|a|t)\b')
# SymPy names the internal ground type in coercion messages; symox does not.
_GROUND_TYPE = re.compile(r" of type <class '[^']*'>")


def show(value):
    if isinstance(value, (XMonomial, SMonomial)):
        return 'Monomial%s' % (tuple(value.exponents),)
    if isinstance(value, type):
        return value.__name__
    if isinstance(value, dict):
        return '{%s}' % ', '.join('%s: %s' % (show(k), show(v)) for k, v in value.items())
    if isinstance(value, (set, frozenset)):
        return '{%s}' % ', '.join(sorted(show(v) for v in value))
    if isinstance(value, tuple):
        return '(%s)' % ', '.join(show(v) for v in value)
    if isinstance(value, list):
        return '[%s]' % ', '.join(show(v) for v in value)
    return _DUMMY.sub(r'D\1', str(value))


def evaluate(lib, code):
    try:
        return show(eval(code, dict(NS[lib])))
    except Exception as exc:  # compared as class name and message
        message = _GROUND_TYPE.sub('', _DUMMY.sub(r'D\1', str(exc)))
        return 'EXC %s: %s' % (type(exc).__name__, message)


POLY_METHODS = [
    "Poly(2*x**2*y+3*x*y**3+1,x,y).terms(order='grlex')",
    "Poly(2*x**2*y+3*x*y**3+x**3,x,y).coeffs(order='grevlex')",
    "Poly(2*x**2*y+3*x*y**3+1,x,y).monoms(order='ilex')",
    "Poly(2*x**2*y+3*x*y**3,x,y).LC(order='grlex')",
    "Poly(2*x**2*y+3*x*y**3+5,x,y).EC()",
    "Poly(2*x**2*y+3*x*y**3,x,y).LM(order='grlex')",
    "Poly(2*x**2*y+3*x*y**3+5,x,y).EM()",
    "Poly(2*x**2*y+3*x*y**3,x,y).LT()",
    "Poly(2*x**2*y+3*x*y**3,x,y).ET(order='grevlex')",
    "(Poly(0,x).terms(), Poly(0,x).coeffs(), Poly(0,x).LC(), Poly(0,x).EC(), Poly(0,x).monoms())",
    "Poly(x*y+1,x,y).as_list()",
    "(Poly(0,x).as_list(), Poly(0,x,y).as_list(), Poly(0,x,y,z).as_list(), Poly(y,x,y).as_list())",
    "Poly(-x**2+2*x-3).abs()",
    "abs(Poly(-x**2-1))",
    "Poly(a*x-1, x, domain='ZZ[a]').abs()",
    "Poly(x-2).neg()",
    "Poly(x+1).add(Poly(y))",
    "Poly(x+1).add(x)",
    "Poly(x+1).sub(Poly(x/2))",
    "Poly(x+1).mul(Poly(x-1))",
    "Poly(x+y,x,y).sqr()",
    "Poly(x+1).pow(3)",
    "Poly(x+1).pow(-1)",
    "Poly(x+1).pow(Rational(1,2))",
    "Poly(x+1).add_ground(2)",
    "Poly(x+1).sub_ground(Rational(1,2))",
    "Poly(x+1).mul_ground(3)",
    "Poly(x+1).mul_ground(Rational(1,2))",
    "Poly(0,x).mul_ground(5)",
    "(Poly(Float(1.5)*x+1).mul_ground(2), Poly(Float(1.5)*x+1).add_ground(Rational(1,2)))",
    "Poly(Float(3.0)*x+1).quo_ground(2)",
    "Poly(a*x+1, x).mul_ground(a)",
    "Poly(a*x+a, x).exquo_ground(a)",
    "Poly(a*x+1, x).exquo_ground(a)",
    "Poly(-5*x+7).quo_ground(3)",
    "Poly(x+3, modulus=5).quo_ground(2)",
    "Poly(x+3).quo_ground(0)",
    "Poly(x+3).exquo_ground(0)",
    "Poly(6*x+4).exquo_ground(2)",
    "Poly(6*x+3).exquo_ground(2)",
    "Poly(6*x+3, domain=QQ).exquo_ground(2)",
    "(Poly(x+3,modulus=5).mul_ground(3), Poly(x+3,modulus=5).neg(), Poly(x+3,modulus=5).sqr())",
    "Poly(3*x**2+4,modulus=5).terms()",
    "(Poly(x+1).eq(Poly(x+1,domain=QQ)), Poly(x+1).eq(Poly(x+1,domain=QQ),strict=True), Poly(x).ne(Poly(y)))",
    "Poly.from_dict({(1,0):2,(0,1):3}, x, y)",
    "Poly.from_dict({(1,):2})",
    "Poly.from_list([1,2,3], x)",
    "Poly.from_list([1,2,3], x, domain=QQ)",
    "Poly.from_list([1,2,3], x, y)",
    "Poly.from_poly(Poly(x/2+1), x, y)",
    "Poly.from_poly(Poly(x/2+1))",
    "Poly.from_expr(x*y+1, y)",
    "Poly(x,x).unify(Poly(y/2,y))",
    "Poly(x,x).unify(3)",
    "Poly(x,x).unify(y)",
    "(Poly(x*y,x,y).one, Poly(x,domain=QQ).zero)",
    "Poly(x**2+1).coeff(x)",
    "Poly(2*x+4).as_content_primitive()",
    "(Poly(x/2+1).retract(), Poly(x+1,domain=QQ).retract(), Poly(x+1,domain=QQ).retract(field=True))",
    "Poly(x+1, extension=sqrt(2)).retract()",
    "Poly(x/2+Float(0.25)).to_exact()",
    "Poly(x/2+7).set_modulus(3)",
    "Poly(x**2+7).set_modulus(5)",
    "Poly(x+a*y,x,domain='ZZ[a,y]').inject()",
    "Poly(x+a*y,x,domain='ZZ[a,y]').inject(front=True)",
    "Poly(x/2+a*y,x,domain='QQ[a,y]').inject()",
    "Poly(x*y+a*x+1,x,y,a).eject(x)",
    "Poly(x*y+a*x+1,x,y,a).eject(a)",
    "Poly(x*y/2+a*x+1,x,y,a).eject(a)",
    "Poly(x/2+1).rat_clear_denoms(Poly(x/3))",
    "Poly(x/2+1).rat_clear_denoms(x)",
    "Poly(x/a+1, x).rat_clear_denoms(Poly(x/b, x))",
    "Poly(x+sqrt(2), extension=sqrt(2)).lift()",
    "Poly(3*x+I, x, extension=I).lift()",
    "Poly(x+sqrt(2)*y, x, y, extension=sqrt(2)).lift()",
    "Poly(x-2**Rational(1,3), extension=2**Rational(1,3)).lift()",
    "Poly(x+1).lift()",
    "Poly(x**2-1, extension=sqrt(3)).sqf_norm()",
    "Poly(x**2-y,x,y, extension=sqrt(3)).sqf_norm()",
    "Poly(x**2+1, extension=I).sqf_norm()",
    "Poly(x**3-2, extension=2**Rational(1,3)).sqf_norm()",
    "Poly(1+x+x**2).revert(5)",
    "Poly(1-x+3*x**2-x**5).revert(8)",
    "Poly(2+x/3).revert(4)",
    "Poly(x+1).revert(0)",
    "Poly(x, domain=QQ).revert(3)",
    "Poly(2+x, modulus=5).revert(3)",
    "Poly(x**4-x**2+1).nth_power_roots_poly(2)",
    "Poly(x**2-x/2+3).nth_power_roots_poly(3)",
    "Poly(x**3-2*x+5).nth_power_roots_poly(3)",
    "Poly(x**3-2*x+5).nth_power_roots_poly(0)",
    "Poly(x*y+1,x,y).nth_power_roots_poly(2)",
    "Poly(x*y+x,x,y).shift_list([1,2])",
    "Poly(x*y+x,x,y).shift_list([1])",
    "Poly(x**3-x**2).real_roots()",
    "Poly(x**3-x**2).real_roots(multiple=False)",
    "Poly(x**3+x).all_roots()",
    "Poly((x**2+4)*(x**2+1)*(x-3)).all_roots()",
    "Poly((x**2+x+1)**2*(x+1)).all_roots(multiple=False)",
    "Poly(x**4+1).all_roots()",
    "(Poly(x**2-2).root(0), Poly(x**2-2).root(-1))",
    "Poly(x**2-2).same_root(sqrt(2), Float(1.41421356237))",
    "Poly(x**2-2).same_root(sqrt(2), Float(1.4))",
    "sorted(Poly(x*(x+3)*(x+5)).dispersionset())",
    "sorted(Poly(x).dispersionset(Poly(x+7)))",
    "(Poly(x*(x+3)*(x+5)).dispersion(), Poly(3,x).dispersion())",
    "Poly(x**2+1).dispersion(Poly(x**2+4*x+5))",
    "Poly(x*y).dispersion()",
    "Poly(5*x**2+7*x-4).trunc(3)",
    "Poly(7*x*y+11*y+5,x,y).trunc(3)",
    "Poly(x**2+y+1,x,y).homogenize(z)",
    "Poly(x**4+3*x**3+2*x+1).slice(1,3)",
    "Poly(y**2+y*z,x,y,z).ltrim(y)",
    "(Poly(2*x+4).is_primitive, Poly(x**2+x+1).is_cyclotomic, Poly(x**2+x+2).is_cyclotomic)",
    "Poly(x**2+1).transform(Poly(x+1), Poly(x-1))",
    "(Poly(x*y+1,x,y).args, Poly(x).func.__name__)",
    "Poly(x*y+2,x,y).atoms()",
    "Poly(x*y+2,x).atoms(Symbol)",
    "(Poly(x*y,x).has(y), Poly(x).has(z))",
    "Poly(2*x**2*y-3*x+y/2+1,x,y).count_ops()",
    "(Poly(x).compare(Poly(x+1)), Poly(x+1).compare(Poly(x)), Poly(x).compare(Poly(x)))",
    "Poly(x+1).xreplace({x:y})",
    "Poly(2*x*y,y,domain='ZZ[x]').doit()",
    "Poly(x+1).copy()",
    "Poly(x+1).replace(x, z)",
    "(Poly(x+1).subs(x, a), Poly(x+1).subs(x,2), Poly(x*y,x,y).subs(x,2), Poly(x*y+a,x,y).subs(a,2))",
    "(Poly(x*y,x,y).subs(x,y), Poly(x*y,x,y).subs({x:a}))",
    "Poly(x+1).match(Wild('w'))",
    "Poly(x).sort_key() < Poly(x+1).sort_key()",
    "Poly(x+a*y,x).free_symbols_in_domain",
    "pickle.loads(pickle.dumps(Poly(x**2+1).mul_ground(3)))",
    "Poly(p*x+1, x).mul_ground(p)",
]

TOP_LEVEL = [
    "poly_from_expr(x+1)",
    "poly_from_expr(x+1, domain=QQ)",
    "poly_from_expr(x*y+1, order='grlex')",
    "poly_from_expr(Poly(x+1))",
    "poly_from_expr(x*y+1, y)",
    "poly_from_expr([1,2])",
    "poly_from_expr(3)",
    "poly_from_expr(x+sqrt(2), extension=True)",
    "poly_from_expr(x+sqrt(2), extension=sqrt(2))",
    "(lambda o: (o.order, o.auto, o.expand, o.gens, o.domain, o.polys))(poly_from_expr(x+1)[1])",
    "parallel_poly_from_expr([x+1, y/2])",
    "parallel_poly_from_expr([3, y])",
    "parallel_poly_from_expr([Poly(x+1,x), y])",
    "parallel_poly_from_expr([3, 4])",
    "parallel_poly_from_expr((x+1, x*y), y, x)",
    "parallel_poly_from_expr([x+1, x**2], modulus=3)",
    "parallel_poly_from_expr([x+1, y], domain='QQ')",
    "degree_list(x**2*y+1)",
    "degree_list(3)",
    "degree_list(Poly(x**2*y, x, y, z))",
    "degree_list(x**2*y+1, y)",
    "(total_degree(x**2*y+x**5), total_degree(3), total_degree(x*y**2, x))",
    "total_degree(Poly(x**2*y+x, x, y))",
    "trunc(5*x**2+7*x-4, 3)",
    "trunc(Poly(5*x**2+7*x-4), 3)",
    "monic(2*x*y+3, x)",
    "monic(3)",
    "monic(x/2+Rational(1,3))",
    "monic(2*x+4, polys=True)",
    "(content(6*x+4), content(x/2+Rational(1,3)), content(6*x*y+4*x))",
    "content(3)",
    "gff(x**2-x)",
    "sqf_norm(x**2-2, extension=sqrt(3))",
    "sqf_norm(x**2-2)",
    "all_roots(x**4-1)",
    "all_roots(x**3-x**2, multiple=False)",
    "all_roots(3)",
    "all_roots(x*y)",
    "real_roots(x**3-x**2, multiple=False)",
    "nth_power_roots_poly(x**2-2, 2)",
    "nth_power_roots_poly(x**2-2, 2, polys=True)",
    "nth_power_roots_poly(Poly(x**2-2), 2)",
    "(LC(2*x**2*y+3*x*y**3, order='grlex'), LM(2*x**2*y+3*x*y**3, order='grlex'), LT(2*x**2*y+3*x*y**3, order='grlex'))",
    "LT(3)",
    "(LC(Poly(3*x+1)), LM(Poly(3*x+1)), LT(Poly(3*x*y+1,x,y)))",
    "(LC(x*y**2+x**2*y, y), LM(x*y**2+x**2*y, y), LT(x*y**2+x**2*y, y))",
]

POLYFUNCS = [
    "horner(x**4+2*x**3+x+5)",
    "horner(x+x*y+x**2*y+x*y**2, y)",
    "horner(x+x*y+x**2*y+x*y**2)",
    "horner(3)",
    "horner(x**3/2+x/3+1)",
    "horner(x**2*y**2+x*y+y+1, x, y)",
    "horner(a*x**2+b*x+1, x)",
    "interpolate([1,4,9,16], x)",
    "interpolate([(1,1),(2,4),(3,9)], x)",
    "interpolate({-1:2, 1:7, 3:-5}, x)",
    "interpolate([(1,1),(2,4)], 2)",
    "interpolate([(1,1),(2,4),(3,9)], 5)",
    "interpolate([(1,a),(2,a+1),(3, 2*a)], x)",
    "interpolate([Float(1.5), Float(2.0), Float(3.5)], x)",
    "interpolate([(Rational(1,2), 1), (Rational(3,2), 2), (3, 0)], x)",
    "interpolate([(a, 1), (b, 2)], x)",
    "interpolate([7], x)",
    "rational_interpolate([(1,-210),(2,-35),(3,105),(4,231),(5,350),(6,465)], 2)",
    "rational_interpolate([(1,-210),(2,-35),(3,105),(4,231),(5,350),(6,465)], 3, X=y)",
    "rational_interpolate([(0,1),(1,Rational(1,2)),(2,Rational(1,5)),(3,Rational(1,10))], 0)",
    "rational_interpolate([(1,2)], 2)",
    "rational_interpolate([1, 2, 5], 1)",
    "viete(a*x**2+2*x+3, [y,z], x)",
    "viete(x**3+2*x+3)",
    "viete(x**2/2+x+3)",
    "viete(x**2+1, [y])",
    "viete(x*y+1)",
    "viete(Poly(3,x))",
    "symmetrize(x**2+y**2)",
    "symmetrize(x**2+y**2, formal=True)",
    "symmetrize([x**2+y**2, x-y])",
    "symmetrize([x**2+y**2, x-y], formal=True)",
    "symmetrize(x**2-y)",
    "symmetrize(a*x**2+y**2, x, y)",
    "symmetrize(x**3+y**3+z**3)",
    "symmetrize(x**2/2+y**2/2+x*y*z)",
    "symmetrize(x**2+y**2, formal=True, symbols=[a,b])",
    "symmetrize(3)",
    "symmetrize(x**2*y+y**2*z+z**2*x + x*y**2+y*z**2+z*x**2)",
    "symmetrize(x**3+y**3+z**3-3*x*y*z, formal=True)",
    "symmetrize(x**3*y+y**2, x, y, formal=True)",
    "symmetrize(Rational(1,3)*x*y + Rational(2,5)*(x+y))",
    "symmetrize([x+y, x-y, x*y], x, y)",
    "symmetrize(a*x+a*y+x*y, x, y, formal=True)",
    "symmetrize(p**2+x**2)",
]

SPECIAL = [
    "(cyclotomic_poly(6, x), cyclotomic_poly(12, x, polys=True), cyclotomic_poly(1,x))",
    "cyclotomic_poly(6, polys=True)",
    "cyclotomic_poly(0, x)",
    "cyclotomic_poly(105, x)",
    "[cyclotomic_poly(n, x) for n in [2,3,4,5,8,10,15,30,36]]",
    "cyclotomic_poly(1155, x).coeff(x, 100)",
    "cyclotomic_poly(Integer(9), x)",
    "pickle.loads(pickle.dumps(cyclotomic_poly(9,x,polys=True)))",
    "(swinnerton_dyer_poly(1,x), swinnerton_dyer_poly(3,x), swinnerton_dyer_poly(2,x,polys=True))",
    "swinnerton_dyer_poly(4,x)",
    "swinnerton_dyer_poly(4,x,polys=True)",
    "Poly(swinnerton_dyer_poly(5, x)).coeffs()[:4]",
    "swinnerton_dyer_poly(0,x)",
    "swinnerton_dyer_poly(2, polys=True)",
    "(symmetric_poly(2,x,a,b), symmetric_poly(0,x,y), symmetric_poly(2,[x,a,b],polys=True))",
    "symmetric_poly(3,x,a)",
    "symmetric_poly(3,x)",
    "symmetric_poly(3, [x,y,z,a])",
    "symmetric_poly(1)",
    "interpolating_poly(3, x)",
    "interpolating_poly(2, x, [1,2], [a,b])",
    "interpolating_poly(3, x, [0,1,2], [1,a,b])",
    "interpolating_poly(2, x, [x,2], [a,b])",
    "interpolating_poly(2, y)",
    "Poly(random_poly(x, 5, -3, 3)).degree()",
    "random_poly(x, 4, 2, 2)",
    "random_poly(x, 3, 1, 1, polys=True)",
]

ORTHO = [
    "[chebyshevt_poly(n,x) for n in range(6)]",
    "[chebyshevu_poly(n,x) for n in range(6)]",
    "chebyshevt_poly(4,x,polys=True)",
    "chebyshevt_poly(-1,x)",
    "chebyshevu_poly(-1,x)",
    "chebyshevt_poly(30, x).coeff(x, 10)",
    "chebyshevt_poly(3.7, x)",
    "chebyshevu_poly(3, polys=True)",
    "[hermite_poly(n,x) for n in range(6)]",
    "[hermite_prob_poly(n,x) for n in range(6)]",
    "hermite_prob_poly(-2,x)",
    "hermite_poly(-2,x)",
    "hermite_poly(15, x, polys=True).LC()",
    "hermite_poly(2)",
    "[legendre_poly(n,x) for n in range(6)]",
    "legendre_poly(3,x,polys=True)",
    "legendre_poly(Rational(7,2),x)",
    "legendre_poly(20, x).coeff(x, 10)",
    "legendre_poly(2, polys=True)",
    "legendre_poly(-1, x)",
    "[laguerre_poly(n,x) for n in range(5)]",
    "laguerre_poly(2,x,alpha=a)",
    "laguerre_poly(2,x,alpha=1/a)",
    "laguerre_poly(2,x,alpha=Float(0.5))",
    "laguerre_poly(2,x,alpha=a,polys=True)",
    "laguerre_poly(3, x, alpha=Rational(1,2), polys=True)",
    "laguerre_poly(2, x, alpha=a+b)",
    "laguerre_poly(-1,x)",
    "[jacobi_poly(n,1,2,x) for n in range(4)]",
    "jacobi_poly(2,a,b,x)",
    "jacobi_poly(2,a,b,x,polys=True)",
    "jacobi_poly(3,Rational(1,2),Rational(-1,2),x,polys=True)",
    "jacobi_poly(2, Float(0.5), Float(1.5), x)",
    "jacobi_poly(3, 0, 0, x, polys=True)",
    "jacobi_poly(2, a, Rational(1,2), x, polys=True)",
    "jacobi_poly(-1,a,b,x)",
    "[gegenbauer_poly(n,a,x) for n in range(4)]",
    "gegenbauer_poly(3,Rational(3,2),x,polys=True)",
    "gegenbauer_poly(2, a, x, polys=True)",
    "gegenbauer_poly(4, Rational(1,2), x)",
    "gegenbauer_poly(-1,a,x)",
    'laguerre_poly(2,x,alpha=pi)',
    'laguerre_poly(2,x,alpha=I)',
    'laguerre_poly(2,x,alpha=Float(0.5)*a)',
    'laguerre_poly(2,x,alpha=sin(a))',
    'laguerre_poly(Float(2.0),x)',
    'laguerre_poly(3,x,alpha=1.1)',
    'jacobi_poly(3,0.3,0.7,x)',
    'jacobi_poly(4,0.3,0.7,x,polys=True)',
    'gegenbauer_poly(5,0.25,x,polys=True)',
    '[spherical_bessel_fn(n, x) for n in range(-4, 6)]',
    'spherical_bessel_fn(3, x, polys=True)',
    'spherical_bessel_fn(2, polys=True)',
    'spherical_bessel_fn(2)',
    'spherical_bessel_fn(-1, x, polys=True)',
    'gegenbauer_poly(3, a+b, x)',
    'jacobi_poly(2, a, a, x, polys=True)',
    'jacobi_poly(2, -1, -1, x)',
    'jacobi_poly(3, -2, -3, x)',
    'laguerre_poly(3, x, alpha=-3)',
    'legendre_poly(40, x).coeff(x, 20)',
    'hermite_poly(3, a)',
    'chebyshevt_poly(5, y+1)',
    'chebyshevt_poly(3, polys=True)',
    'jacobi_poly(2, a, b)',
    'laguerre_poly(2, alpha=a, polys=True)',
    'gegenbauer_poly(2, Rational(1,3), polys=True)',
    'laguerre_poly(0, x, alpha=a, polys=True)',
    'jacobi_poly(1, a, b, x, polys=True)',
    'gegenbauer_poly(1, a, x, polys=True)',
    'gegenbauer_poly(0, a, x)',
    'laguerre_poly(2, x, a)',
    'legendre_poly(Integer(3), x)',
]

# Coefficient grouping of as_expr() and dense coefficient-list construction
# (including composite generators such as ``y + 1``).
CONSTRUCTION = [
    'Poly(x*(a+1), x).as_expr()',
    'Poly(x*(a+1)+x**2*(a-1)/2, x).as_expr()',
    'Poly(x*(pi+2), x).as_expr()',
    'Poly(x*(pi+2), x)',
    'Poly([1, a+1, 3], x)',
    'Poly([1, a+1, 3], x).as_expr()',
    'Poly([16, 0, -20, 0, 5, 0], y+1)',
    'Poly([16, 0, -20, 0, 5, 0], y+1).as_expr()',
    "Poly(x*(a+1), x, domain='ZZ(a)').as_expr()",
    'Poly(x*(a+1)/(b+1), x).as_expr()',
    'Poly(x*y*(a+1), x, y).as_expr()',
    'Poly(x*(sqrt(2)+1), x).as_expr()',
    'Poly(x*(sqrt(2)+1), x, extension=True).as_expr()',
    'Poly([1, 2, 3], y)',
    "Poly([1, 2, 3], y, domain='QQ')",
    'Poly([1, a, 3], x)',
    "Poly([1, 2], y+1, domain='QQ')",
    'Poly([a, 2], sin(y))',
]

PARTFRAC = [
    "apart_list((x**2+1)/(x+1)/(x-2)**2, x)",
    "apart_list(2*x/(x**2-1))",
    "apart_list(x)",
    "apart_list(1/(x**2-2), x)",
    "apart_list((x+Rational(1,2))/(x**2+3*x+2), x)",
    "apart_list((x**3+1)/(x-1), x)",
    "apart_list(1/(x**2+1)/(x-1), x)",
    "apart_list(1/(x*y+1))",
    "apart_list(1/(x**2-2), x, extension=sqrt(2))",
    "assemble_partfrac_list(apart_list((x**2+1)/(x+1)/(x-2)**2, x))",
    "assemble_partfrac_list(apart_list((2*x+3)/(x**3-x), x))",
    "assemble_partfrac_list(apart_list(x**4/(x**2-4)/(x+1)**3, x))",
    "assemble_partfrac_list(apart_list((x+Rational(1,2))/(x**2+3*x+2), x))",
    "assemble_partfrac_list(apart_list((x**3+1)/(x-1)**2, x))",
    "assemble_partfrac_list(apart_list(1/(x-1)**4/(x+2)**2, x))",
    "assemble_partfrac_list([1, Poly(0,x), [([1,2], Lambda(a, 3), Lambda(a, x-a), 1)]])",
    'apart_list(1/(x-2)**2, x)',
    'apart_list(1/(x**3-2), x)',
    'apart_list(Rational(1,2)/(x+1)**2, x)',
    'apart_list(x/(2*x**2+3), x)',
    'apart_list(x/(x**2-1)**2, x)',
    'apart_list((x**5+3*x-1)/((x**2+x+1)**3*(x-3)), x)',
    'apart_list(1/((2*x-1)*(x+1)**2), x)',
    'apart_list((3*x+1)/(4*x**2-1)**2, x)',
    'apart_list(1/(x**4+1), x)',
    'apart_list(1/(x*(x-1)**3*(x**2+2)**2), x)',
    'apart_list((x**3-2*x)/(x**2-2)**3, x)',
    'apart_list(-1/(x-1)**2, x)',
    'apart_list(-x/(3-x)**2, x)',
    'apart_list(2/(4*x**2+4*x+1), x)',
    'apart_list((x+1)/x**3, x)',
    'assemble_partfrac_list(apart_list(x/(x**2-1)**2, x))',
    'assemble_partfrac_list(apart_list(1/((2*x-1)*(x+1)**2), x))',
    'assemble_partfrac_list(apart_list(-x/(3-x)**2, x))',
]

# Native recognition of cyclotomic polynomials (dup_cyclotomic_p).
CYCLOTOMIC = [
    '[cyclotomic_poly(n, x, polys=True).is_cyclotomic for n in range(1, 40)]',
    '[Poly(x**n - 1, x).is_cyclotomic for n in range(1, 8)]',
    'Poly(x**4 + x**2 + 1, x).is_cyclotomic',
    "Poly(x**2 + x + 1, x, domain='QQ').is_cyclotomic",
    'Poly(x**2 + x/2 + 1, x).is_cyclotomic',
    'Poly(x**2 + x + 1, x, modulus=7).is_cyclotomic',
    'Poly(x**2 + y + 1, x, y).is_cyclotomic',
    'Poly(x**2 + a*x + 1, x).is_cyclotomic',
    'Poly(x**2 + sqrt(2)*x + 1, x, extension=True).is_cyclotomic',
    'Poly(1, x).is_cyclotomic',
    'Poly(x, x).is_cyclotomic',
    'Poly(x + 1, x).is_cyclotomic',
    'Poly(-x - 1, x).is_cyclotomic',
    'Poly(x**8 - x**7 + x**5 - x**4 + x**3 - x + 1, x).is_cyclotomic',
    'Poly(x**8 + x**7 + x**5 + x**4 + x**3 + x + 1, x).is_cyclotomic',
    'Poly((x**2+x+1).subs(x, x**3), x).is_cyclotomic',
    'Poly(x**2 - x + 1, x).is_cyclotomic',
    'Poly(x**4 + 1, x).is_cyclotomic',
    'Poly(x**4 - x**2 + 1, x).is_cyclotomic',
    'Poly(x**6 + x**5 + x**4 + x**3 + x**2 + x + 1, x).is_cyclotomic',
]

# Native ZZ/QQ dispersion sets and the generic fallback.
DISPERSION = [
    'sorted(Poly(x**4-3*x**2+1).dispersionset(Poly(x**4-3*x**2+1).shift(-3)))',
    'sorted(Poly((x-3)*(x+3)).dispersionset())',
    'sorted(Poly(x/2+Rational(1,3)).dispersionset(Poly(3*x - 4)))',
    'sorted(Poly((2*x+1)**2*(2*x+7)).dispersionset())',
    'sorted(Poly(x**2+1, x).dispersionset(Poly(x**2-2*x+2, x)))',
    'Poly(x**2+1, x).dispersion(Poly(y, x))',
    'Poly(Float(1.5)*x+1).dispersion()',
    'Poly(x+a, x).dispersion(Poly(x+a+3, x))',
    'sorted(Poly(x*(x-1)**3*(x+5)**2, x).dispersionset(Poly((x+2)*(x-8), x)))',
]

# SymPy's dmp_cancel normalization of Poly.cancel (sign and content).
CANCEL = [
    'Poly(-1, x).cancel(Poly((x-1)**2, x))',
    'Poly(-x, x).cancel(Poly((3-x)**2, x))',
    'Poly(-x, x).cancel(Poly(x**2+1, x))',
    'Poly(x, x).cancel(Poly(-x**2+1, x))',
    'Poly(-2*x, x).cancel(Poly(-4*x**2+2, x))',
    'Poly(Rational(3,7)*x, x).cancel(Poly(x**2+Rational(1,3), x))',
    'Poly(Rational(-3,7)*x, x).cancel(Poly(x**2+Rational(1,3), x))',
    "Poly(Rational(3,7)*x, x, domain='QQ').cancel(Poly(x**2+Rational(1,3), x))[1].div(Poly(x**2+Rational(1,3), x).cancel(Poly(1,x))[1], auto=True)",
    "Poly(3*x, x, domain='QQ').cancel(Poly(x**2+1, x, domain='QQ'))",
    "Poly(-3*x, x, domain='QQ').cancel(Poly(x**2+1, x, domain='QQ'))",
    "Poly(-3*x, x, domain='QQ').cancel(Poly(-x**2+1, x, domain='QQ'))",
    'Poly(-1, x).cancel(Poly((x-1)**2, x), include=True)',
    'Poly(x**2-1, x).cancel(Poly(-2*x+2, x))',
    'Poly(x**2-1, x).cancel(Poly(-2*x+2, x), include=True)',
    'Poly(Rational(1,2)*x, x).cancel(Poly(-3*x**2+3, x))',
    'Poly(Rational(1,2)*x, x).cancel(Poly(-3*x**2+3, x), include=True)',
]

BASIC_AND_ROOTS = [
    'Poly(2*x**2-1).make_monic_over_integers_by_scaling_roots()',
    'Poly(x**2/3+x/2+1).make_monic_over_integers_by_scaling_roots()',
    'Poly(x**2+1).make_monic_over_integers_by_scaling_roots()',
    'Poly(x*y+1).make_monic_over_integers_by_scaling_roots()',
    'Poly(x**2+sqrt(2), x).make_monic_over_integers_by_scaling_roots()',
    'Poly(x**3-x).which_real_roots([0, 1, -1, 2, 1])',
    'Poly(x**2-2).which_real_roots([sqrt(2), -sqrt(2), 1, 2])',
    'Poly(x**2+1).which_all_roots([I, -I, 1])',
    'Poly(x**3-x).which_all_roots([5, 6])',
    'Poly(x*y).which_real_roots([0])',
    'Poly(x*y).which_all_roots([0])',
    '(Poly(x**2+2*x+1).count(x), Poly(x**2+2*x+1).count(Symbol))',
    '(Poly(x**2+1).has_free(x), Poly(x**2+1).has_free(y), Poly(x**2+1).has_free())',
    'Poly(x**2+1).has_free([x])',
    '(Poly(x**2+1).has_xfree({x}), Poly(x**2+1).has_xfree({y}))',
    'Poly(x**2+1).has_xfree([x])',
    '(Poly(x**2+1).is_same(Poly(x**2+1)), Poly(x**2+1).is_same(x**2+1), Poly(x+1).is_same(Poly(x+2)))',
    '(Poly(x**2+1).dummy_eq(Poly(x**2+1)), Poly(x**2+1).dummy_eq(Poly(x**2+2)))',
    'Poly(-x).could_extract_minus_sign()',
    '(Poly(x**2+1).matches(Poly(x**2+1)), Poly(x**2+1).matches(x**2+1))',
    'Poly(x**2+2*x+1).rcall(2)',
    'Poly.fromiter([x**2, x])',
    'Poly(sqrt(x**2)*y, y).refine(Q.positive(x))',
    "Poly(x, x, domain='QQ').refine()",
    "Poly(x**2+2*x, x, domain='QQ').simplify()",
    '(Poly(x**2+1).is_hypergeometric(x), Poly(x*y+1, x).is_hypergeometric(y))',
]

MONOMIALS = [
    "sorted(map(str, itermonomials([x,y],2)))",
    "sorted(map(str, itermonomials([x,y],[2,1],[1,0])))",
    "sorted(map(str, itermonomials([x,y,z],3,2)))",
    "list(itermonomials([x],[1],[2]))",
    "list(itermonomials([x,y],-1))",
    "(str(Monomial((1,2))), str(Monomial((1,2),(x,y))), Monomial((1,2))*Monomial((3,4)), Monomial((1,2)).as_expr(x,y))",
    "(Monomial((3,4))/Monomial((1,2)), Monomial((1,2))**3, Monomial((1,5)).gcd(Monomial((2,3))), Monomial((1,5)).lcm((2,3)))",
    "Monomial((1,2))/Monomial((2,1))",
    "Monomial(x**2*y)",
    "Monomial(2*x)",
    "(lex((1,2)), grlex((1,2)), grevlex((1,2,3)), ilex((1,2)), igrlex((1,2)), igrevlex((1,2,3)))",
    "(str(igrevlex), repr(grlex), lex == grlex, ilex == ilex)",
    "sorted([(1,2),(2,0),(0,3),(1,1)], key=grevlex)",
    "pickle.loads(pickle.dumps(igrlex)) == igrlex",
]

ERRORS = [
    "[issubclass(getattr(M, n), M.BasePolynomialError) for n in ['ComputationFailed','EvaluationFailed','ExtraneousFactors','GeneratorsError','GeneratorsNeeded','HeuristicGCDFailed','HomomorphismFailed','MultivariatePolynomialError','NotAlgebraic','NotReversible','OperationNotSupported','PolificationFailed','PolynomialDivisionFailed','UnivariatePolynomialError']]",
    "[issubclass(getattr(M, n), M.PolynomialError) for n in ['GeneratorsError','GeneratorsNeeded','MultivariatePolynomialError','UnivariatePolynomialError','PolificationFailed','NotAlgebraic','ComputationFailed']]",
    "str(M.ComputationFailed('LT', 1, M.PolificationFailed({}, 3, 3)))",
    "str(M.PolynomialDivisionFailed(Poly(x), Poly(0, x), QQ))[:40]",
    "str(M.OperationNotSupported(Poly(x), 'foo')).split(' by ')[0]",
    "pickle.loads(pickle.dumps(M.NotReversible('zero is not reversible'))).args",
]


def _check(code):
    assert evaluate('x', code) == evaluate('s', code)


@pytest.mark.parametrize('code', POLY_METHODS)
def test_poly_methods(code):
    _check(code)


@pytest.mark.parametrize('code', TOP_LEVEL)
def test_top_level_polytools(code):
    _check(code)


@pytest.mark.parametrize('code', POLYFUNCS)
def test_polyfuncs(code):
    _check(code)


@pytest.mark.parametrize('code', SPECIAL)
def test_special_polys(code):
    _check(code)


@pytest.mark.parametrize('code', ORTHO)
def test_orthogonal_polys(code):
    _check(code)


@pytest.mark.parametrize('code', CONSTRUCTION)
def test_construction_and_as_expr(code):
    _check(code)


@pytest.mark.parametrize('code', CYCLOTOMIC)
def test_is_cyclotomic(code):
    _check(code)


@pytest.mark.parametrize('code', DISPERSION)
def test_dispersion(code):
    _check(code)


@pytest.mark.parametrize('code', CANCEL)
def test_poly_cancel_normalization(code):
    _check(code)


@pytest.mark.parametrize('code', PARTFRAC)
def test_partial_fraction_lists(code):
    _check(code)


@pytest.mark.parametrize('code', BASIC_AND_ROOTS)
def test_basic_protocol_and_root_selection(code):
    _check(code)


@pytest.mark.parametrize('code', MONOMIALS)
def test_monomials_and_orderings(code):
    _check(code)


@pytest.mark.parametrize('code', ERRORS)
def test_polyerrors(code):
    _check(code)


def test_exports():
    for name in ['all_roots', 'content', 'degree_list', 'gff', 'monic',
                 'nth_power_roots_poly', 'parallel_poly_from_expr',
                 'poly_from_expr', 'sqf_norm', 'total_degree', 'trunc', 'horner',
                 'interpolate', 'rational_interpolate', 'symmetrize', 'viete',
                 'cyclotomic_poly', 'interpolating_poly', 'random_poly',
                 'swinnerton_dyer_poly', 'symmetric_poly', 'chebyshevt_poly',
                 'chebyshevu_poly', 'hermite_poly', 'hermite_prob_poly',
                 'jacobi_poly', 'laguerre_poly', 'legendre_poly', 'itermonomials',
                 'Monomial', 'lex', 'grlex', 'grevlex', 'ilex', 'igrlex',
                 'igrevlex', 'apart_list', 'assemble_partfrac_list']:
        assert hasattr(symox, name) and hasattr(sympy, name)
        assert hasattr(symox.polys, name)
    assert not hasattr(symox, 'gegenbauer_poly') and not hasattr(sympy, 'gegenbauer_poly')
    import symox.polys.polyerrors as xe
    import sympy.polys.polyerrors as se
    for name in se.__all__:
        if hasattr(xe, name):
            assert getattr(symox, name) is getattr(xe, name)


def test_basic_protocol_on_poly():
    x, y = symox.symbols('x y')
    p = symox.Poly(x*y + 1, x, y)
    assert p.is_Poly and not p.is_Number and p.is_commutative
    assert p.is_real is None and p.is_Atom is False
    assert p.func(*p.args) == p
    assert p.find(x) == {x}
    assert p.as_dummy() == p and p.rewrite(symox.cos) == p
    assert symox.Poly.class_key() == sympy.Poly.class_key()


def test_documented_deviations():
    x, y, a = symox.symbols('x y a')
    # SymPy 1.14 overwrites the x**2 term (x**2 + x*y); symox keeps both terms.
    assert symox.Poly(x**2 + y + 1, x, y).homogenize(x) == symox.Poly(2*x**2 + x*y, x, y)
    # Non-linear factors assemble into a RootSum, as in SymPy 1.14.
    assert str(symox.assemble_partfrac_list(symox.apart_list(1/(x**2 + 1)/(x - 1), x))) == (
        'RootSum(_w**2 + 1, Lambda(_a, (_a**2/4 + _a/4)/(-_a + x))) + 1/(2*(x - 1))')
    # Composite moduli build a polynomial over the ring Z/4Z.
    assert repr(symox.Poly(x**2 + 7).set_modulus(4)) == 'Poly(x**2 - 1, x, modulus=4)'
    # SymPy overflows on candidates whose value is exactly zero but the
    # count is too small (non-real roots); symox keeps the exact roots.
    assert symox.Poly(x**2 + 1).which_real_roots([symox.I, -symox.I]) == [symox.I, -symox.I]
    # Parametric coefficient rings follow the reference partial fractions.
    assert str(symox.assemble_partfrac_list(symox.apart_list(a/(x - a)/(x + 1), x))) == (
        '(-2*a/(a**2 + 2*a + 1) + (-a**2 + a)/(a**2 + 2*a + 1))/(x + 1)'
        ' + (2*a**2/(a**2 + 2*a + 1) + (-a**2 + a)/(a**2 + 2*a + 1))/(-a + x)')


def test_random_poly_bounds():
    x = symox.Symbol('x')
    for _ in range(20):
        p = symox.random_poly(x, 6, -5, 5, polys=True)
        assert p.degree() == 6
        assert all(-5 <= c <= 5 for c in p.all_coeffs())
