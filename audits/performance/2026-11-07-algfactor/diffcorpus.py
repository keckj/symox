"""Differential corpus for the algfactor lot (wave 45): factorization,
square-free decomposition, cancel/gcd and roots with algebraic numbers,
natural (algebraic numbers as generators), with `extension=`/`gaussian=`
and `greedy=False`, CRootOf extensions and exception classes.

Usage: python diffcorpus.py sympy|symox OUT.jsonl [case names...]

Each record keeps the printed result (or the exception class and message)
and a `value`, the result as an expression, used by compare.py to tell a
representation difference from a wrong value. A factor list
`(c, [(f, k), ...])` has the value `c*prod(f**k)`."""
import json, sys, time

mode, out = sys.argv[1], sys.argv[2]
if mode == 'sympy':
    import sympy as M
else:
    import symox as M
names = ('symbols', 'sqrt', 'Poly', 'QQ', 'ZZ', 'QQ_I', 'ZZ_I', 'I', 'pi', 'S', 'Rational',
         'cancel', 'factor', 'gcd', 'lcm', 'sqf_list', 'factor_list', 'sqf', 'cofactors', 'div',
         'quo', 'rem', 'resultant', 'real_roots', 'all_roots', 'nroots', 'ground_roots', 'roots',
         'CRootOf', 'minimal_polynomial', 'primitive_element', 'to_number_field', 'AlgebraicNumber',
         'field_isomorphism', 'expand')
(symbols, sqrt, Poly, QQ, ZZ, QQ_I, ZZ_I, I, pi, S, Rational, cancel, factor, gcd, lcm, sqf_list,
 factor_list, sqf, cofactors, div, quo, rem, resultant, real_roots, all_roots, nroots, ground_roots,
 roots, CRootOf, minimal_polynomial, primitive_element, to_number_field, AlgebraicNumber,
 field_isomorphism, expand) = (getattr(M, n) for n in names)

x, y, z, t = symbols('x y z t')
s2, s3 = sqrt(2), sqrt(3)
c2 = 2**Rational(1, 3)
r = CRootOf(x**3 + x + 3, 0)
rc = CRootOf(x**3 + x + 3, 1)
r2 = CRootOf(x**3 - 2, 0)
q2 = CRootOf(x**2 - 2, 1)

cases = [
    # natural factor lists: algebraic numbers and constants are generators
    ('fl (x-s2)(x+s2)', lambda: factor_list((x - s2)*(x + s2))),
    ('fl (x-s2*t)**2', lambda: factor_list((x - s2*t)**2)),
    ('fl (x-s2*pi)(x+s2*pi)', lambda: factor_list((x - s2*pi)*(x + s2*pi))),
    ('fl (x-s2*pi)(x+s2*pi) x', lambda: factor_list((x - s2*pi)*(x + s2*pi), x)),
    ('fl (x+1)(x-s2)**2', lambda: factor_list((x + 1)*(x - s2)**2)),
    ('fl (x-s2)**2(x+s2)', lambda: factor_list((x - s2)**2*(x + s2))),
    ('fl x**2-2', lambda: factor_list(x**2 - 2)),
    ('fl expanded (x-s2)**2', lambda: factor_list(x**2 - 2*s2*x + 2)),
    ('fl 2x-2s2', lambda: factor_list(2*x - 2*s2)),
    ('fl s2*x+2', lambda: factor_list(s2*x + 2)),
    ('fl (s2x+1)(x-s3)', lambda: factor_list((s2*x + 1)*(x - s3))),
    ('fl (x-s2)(y+s3)', lambda: factor_list((x - s2)*(y + s3))),
    ('fl (x-cbrt2)(x+1)', lambda: factor_list((x - c2)*(x + 1))),
    ('fl (x-I)(x+I)', lambda: factor_list((x - I)*(x + I))),
    ('fl (x-s2*I)(x+1)', lambda: factor_list((x - s2*I)*(x + 1))),
    ('fl (x-s2)**3', lambda: factor_list((x - s2)**3)),
    ('fl x(x-s2)(x**2-3)', lambda: factor_list(x*(x - s2)*(x**2 - 3))),
    ('fl frac (x-s2)/(x+1)', lambda: factor_list((x - s2)/(x + 1), frac=True)),
    ('fl frac (x-s2)/(x+s2)**2', lambda: factor_list((x - s2)/(x + s2)**2, frac=True)),
    ('factor (x-s2)(x+s2)', lambda: factor((x - s2)*(x + s2))),
    ('factor expanded (x-s2)**2', lambda: factor(x**2 - 2*s2*x + 2)),
    ('factor (x-s2*t)**2', lambda: factor((x - s2*t)**2)),
    ('factor (x-s2*pi)(x+s2*pi)', lambda: factor((x - s2*pi)*(x + s2*pi))),
    ('fl (x-s2)(x+s2) x', lambda: factor_list((x - s2)*(x + s2), x)),
    ('fl (x-s2*t)**2 x', lambda: factor_list((x - s2*t)**2, x)),
    ('fl (x-s2*t)**2 x t', lambda: factor_list((x - s2*t)**2, x, t)),
    ('fl (x-pi)**2(x+pi)', lambda: factor_list((x - pi)**2*(x + pi))),
    ('fl (xy-s2)(x+y)', lambda: factor_list((x*y - s2)*(x + y))),
    ('fl polys (x-s2)**2(x+1)', lambda: factor_list((x - s2)**2*(x + 1), polys=True)),
    ('fl 3(x-s2)**2', lambda: factor_list(3*(x - s2)**2)),
    ('fl -(x-s2)(x+s2)', lambda: factor_list(-(x - s2)*(x + s2))),
    ('fl s2*(x-1)**2', lambda: factor_list(s2*(x - 1)**2)),
    ('fl (x-s2)**2*(x-s3)**3', lambda: factor_list((x - s2)**2*(x - s3)**3)),
    ('fl (x-s2*y)**2*(x+y)', lambda: factor_list((x - s2*y)**2*(x + y))),
    ('factor (x-s2*y)**2*(x+y)', lambda: factor((x - s2*y)**2*(x + y))),
    ('fl xy-s2x+y-s2', lambda: factor_list(x*y - s2*x + y - s2)),
    ('factor xy-s2x+y-s2', lambda: factor(x*y - s2*x + y - s2)),
    ('fl (x-s2)(x+s2)(x-pi)', lambda: factor_list((x - s2)*(x + s2)*(x - pi))),
    ('fl (x+s2)**2 x**2', lambda: factor_list(x**2*(x + s2)**2)),
    ('fl (x-1/s2)(x+1)', lambda: factor_list((x - 1/s2)*(x + 1))),
    ('fl (2x-s2)(x+s2)', lambda: factor_list((2*x - s2)*(x + s2))),
    ('fl ((x-s2)**2-3)', lambda: factor_list((x - s2)**2 - 3)),
    # constant products: the generators are the algebraic numbers themselves
    ('fl (1-s2)(1+s2)', lambda: factor_list((1 - s2)*(1 + s2))),
    ('fl (1-s2)(1+s2) polys', lambda: factor_list((1 - s2)*(1 + s2), polys=True)),
    ('fl (1-s2)(1+s2) x', lambda: factor_list((1 - s2)*(1 + s2), x)),
    ('fl s2(1-s2)', lambda: factor_list(s2*(1 - s2))),
    ('fl (1-s3)**2(1+s3)**2', lambda: factor_list((1 - s3)**2*(1 + s3)**2)),
    ('fl (s2-s3)(s2+s3)', lambda: factor_list((s2 - s3)*(s2 + s3))),
    ('fl (1+pi)(1-pi)', lambda: factor_list((1 + pi)*(1 - pi))),
    ('fl (s2+pi)**2', lambda: factor_list((s2 + pi)**2)),
    ('sqfl (1-s2)(1+s2)', lambda: sqf_list((1 - s2)*(1 + s2))),
    ('sqfl (1-s2)**2(1+s2)', lambda: sqf_list((1 - s2)**2*(1 + s2))),
    ('factor (1-s2)(1+s2)', lambda: factor((1 - s2)*(1 + s2))),
    # natural square-free lists
    ('sqfl (x+1)(x-s2)**2', lambda: sqf_list((x + 1)*(x - s2)**2)),
    ('sqfl (x-s2)**2(x+s2)', lambda: sqf_list((x - s2)**2*(x + s2))),
    ('sqfl (x-s2*t)**2', lambda: sqf_list((x - s2*t)**2)),
    ('sqfl (x**2-1)(x+1)', lambda: sqf_list((x**2 - 1)*(x + 1))),
    ('sqfl (x-1)**2(x+1)', lambda: sqf_list((x - 1)**2*(x + 1))),
    ('sqfl (x**2-1)**2(x+1)', lambda: sqf_list((x**2 - 1)**2*(x + 1))),
    ('sqf (x-s2)**2(x+1)', lambda: sqf((x - s2)**2*(x + 1))),
    ('sqf (x**2-1)(x+1)', lambda: sqf((x**2 - 1)*(x + 1))),
    ('sqfl (x-s2)(x+s2)', lambda: sqf_list((x - s2)*(x + s2))),
    ('sqfl expanded', lambda: sqf_list(x**3 - s2*x**2 - 2*x + 2*s2)),
    ('sqfl (x-pi)**2(x+1)', lambda: sqf_list((x - pi)**2*(x + 1))),
    ('sqfl 3(x-s2)**2(x+1)**2', lambda: sqf_list(3*(x - s2)**2*(x + 1)**2)),
    ('sqfl (x-s2*pi)**2', lambda: sqf_list((x - s2*pi)**2)),
    ('sqfl (x-s2)**2 x', lambda: sqf_list((x - s2)**2, x)),
    ('sqfl (x-y)**2(x+y)', lambda: sqf_list((x - y)**2*(x + y))),
    ('sqfl (x**2-y**2)(x+y)', lambda: sqf_list((x**2 - y**2)*(x + y))),
    ('sqfl x**2(x+1)**2', lambda: sqf_list(x**2*(x + 1)**2)),
    ('sqfl (x-s2)**2(x-s3)**3', lambda: sqf_list((x - s2)**2*(x - s3)**3)),
    ('sqf (x-s2*t)**2', lambda: sqf((x - s2*t)**2)),
    # Poly methods (generators x, sqrt(2))
    ('P fl (x-s2)(x+s2)', lambda: Poly((x - s2)*(x + s2)).factor_list()),
    ('P sqfl (x-s2)**2(x+1)', lambda: Poly((x - s2)**2*(x + 1)).sqf_list()),
    ('P fl expanded (x-s2)**2', lambda: Poly(x**2 - 2*s2*x + 2).factor_list()),
    ('P fl (x-s2*t)**2', lambda: Poly((x - s2*t)**2).factor_list()),
    ('P fl x,s2 explicit', lambda: Poly(x**2 - 2*s2*x + 2, x, s2).factor_list()),
    # extension= and gaussian=
    ('factor x**2-2 ext s2', lambda: factor(x**2 - 2, extension=s2)),
    ('fl x**2-2 ext s2', lambda: factor_list(x**2 - 2, extension=s2)),
    ('factor x**2-2t**2 ext s2', lambda: factor(x**2 - 2*t**2, extension=s2)),
    ('factor x**2-2t**2 x ext s2', lambda: factor(x**2 - 2*t**2, x, extension=s2)),
    ('fl x**2-2t**2 ext s2', lambda: factor_list(x**2 - 2*t**2, extension=s2)),
    ('fl x**2-2t**2 x ext s2', lambda: factor_list(x**2 - 2*t**2, x, extension=s2)),
    ('sqfl x**2-2t**2 x ext s2', lambda: sqf_list(x**2 - 2*t**2, x, extension=s2)),
    ('factor x**2-2y**2 ext s2', lambda: factor(x**2 - 2*y**2, extension=s2)),
    ('factor x**4-4 ext s2', lambda: factor(x**4 - 4, extension=s2)),
    ('factor x**4+4 ext I', lambda: factor(x**4 + 4, extension=I)),
    ('factor x**2+1 gaussian', lambda: factor(x**2 + 1, gaussian=True)),
    ('factor x**2+y**2 gaussian', lambda: factor(x**2 + y**2, gaussian=True)),
    ('fl x**2+y**2 gaussian', lambda: factor_list(x**2 + y**2, gaussian=True)),
    ('factor x**4-1 ext I', lambda: factor(x**4 - 1, extension=I)),
    ('factor x**4-1 ext [s2,I]', lambda: factor(x**4 - 1, extension=[s2, I])),
    ('factor x**4-5x**2+6 ext [s2,s3]', lambda: factor(x**4 - 5*x**2 + 6, extension=[s2, s3])),
    ('factor x**2-3 ext s2', lambda: factor(x**2 - 3, extension=s2)),
    ('factor (x-s2y)**2 ext s2', lambda: factor(x**2 - 2*s2*x*y + 2*y**2, extension=s2)),
    ('sqfl expanded ext s2', lambda: sqf_list(x**3 - s2*x**2 - 2*x + 2*s2, extension=s2)),
    ('sqfl (x-s2)**2(x+1) ext s2', lambda: sqf_list((x - s2)**2*(x + 1), extension=s2)),
    ('factor (x-s2)**2(x+1) ext True', lambda: factor((x - s2)**2*(x + 1), extension=True)),
    ('factor expanded ext True', lambda: factor(x**2 - 2*s2*x + 2, extension=True)),
    ('fl expanded ext True', lambda: factor_list(x**2 - 2*s2*x + 2, extension=True)),
    ('factor x**3-2 ext cbrt2', lambda: factor(x**3 - 2, extension=c2)),
    ('factor x**6-2 ext cbrt2', lambda: factor(x**6 - 2, extension=c2)),
    ('factor x**2y**2-2 ext s2', lambda: factor(x**2*y**2 - 2, extension=s2)),
    ('factor x**2-2y**2z**2 ext s2', lambda: factor(x**2 - 2*y**2*z**2, extension=s2)),
    ('fl x**2-2 gaussian', lambda: factor_list(x**2 - 2, gaussian=True)),
    ('P fl x**2-2 ext s2', lambda: Poly(x**2 - 2, x, extension=s2).factor_list()),
    ('P fl x**2-2y**2 ext s2', lambda: Poly(x**2 - 2*y**2, x, y, extension=s2).factor_list()),
    ('P fl x**4+1 ext s2', lambda: Poly(x**4 + 1, x, extension=s2).factor_list()),
    ('factor x**4+1 ext s2', lambda: factor(x**4 + 1, extension=s2)),
    ('factor x**4+1 ext [s2,I]', lambda: factor(x**4 + 1, extension=[s2, I])),
    ('P fl x**2+y**2 gaussian', lambda: Poly(x**2 + y**2, x, y, gaussian=True).factor_list()),
    ('P sqfl (x-s2)**2(x-1) ext s2', lambda: Poly((x - s2)**2*(x - 1), x, extension=s2).sqf_list()),
    ('factor x**2-2 x ext s2', lambda: factor(x**2 - 2, x, extension=s2)),
    ('factor x**2-2pi**2 ext s2', lambda: factor(x**2 - 2*pi**2, extension=s2)),
    ('factor x**2-2pi**2 x ext s2', lambda: factor(x**2 - 2*pi**2, x, extension=s2)),
    ('factor s2x**2+2x+s2/2 ext True', lambda: factor(s2*x**2 + 2*x + s2/2, extension=True)),
    ('fl s2x**2+2x+s2/2 ext True', lambda: factor_list(s2*x**2 + 2*x + s2/2, extension=True)),
    ('factor x**2+2Ix-1 ext True', lambda: factor(x**2 + 2*I*x - 1, extension=True)),
    ('factor x**2-2 ext s3', lambda: factor(x**2 - 2, extension=s3)),
    ('factor x**2y-2y ext s2', lambda: factor(x**2*y - 2*y, extension=s2)),
    ('factor (x**2-2)/(x**2-8) ext s2', lambda: factor((x**2 - 2)/(x**2 - 8), extension=s2)),
    ('sqf expanded ext s2', lambda: sqf(x**3 - s2*x**2 - 2*x + 2*s2, extension=s2)),
    ('factor x**3y-2xy ext s2', lambda: factor(x**3*y - 2*x*y, extension=s2)),
    ('factor x**4-4y**4 ext [s2,I]', lambda: factor(x**4 - 4*y**4, extension=[s2, I])),
    ('factor x**2+2y**2 ext s2*I', lambda: factor(x**2 + 2*y**2, extension=s2*I)),
    ('factor x**2-3y**2t**2 ext s3', lambda: factor(x**2 - 3*y**2*t**2, extension=s3)),
    ('fl x**2+y**2 ext I', lambda: factor_list(x**2 + y**2, extension=I)),
    ('sqfl (x-s2y)**2(x+y) ext s2', lambda: sqf_list((x - s2*y)**2*(x + y), extension=s2)),
    ('factor (x-s2y)**2(x+y) ext s2', lambda: factor((x - s2*y)**2*(x + y), extension=s2)),
    ('factor x**2-2y**2+x-s2y ext s2', lambda: factor(x**2 - 2*y**2 + x - s2*y, extension=s2)),
    ('factor xy-s2x+y-s2 ext s2', lambda: factor(x*y - s2*x + y - s2, extension=s2)),
    ('factor x**2y+s2xy ext True', lambda: factor(x**2*y + s2*x*y, extension=True)),
    ('fl (x**2-2)(y**2-2) ext s2', lambda: factor_list((x**2 - 2)*(y**2 - 2), extension=s2)),
    ('factor x**2-2 ext [s2,s3]', lambda: factor(x**2 - 2, extension=[s2, s3])),
    ('factor x**2-6 ext [s2,s3]', lambda: factor(x**2 - 6, extension=[s2, s3])),
    # cancel, gcd, lcm with greedy=False
    ('cancel greedy=False x+s2', lambda: cancel((x**2 - 2)/(x + s2), greedy=False)),
    ('cancel greedy=False x-s2', lambda: cancel((x**2 - 2)/(x - s2), greedy=False)),
    ('cancel greedy=False t', lambda: cancel((x**2 - 2*t**2)/(x + s2*t), greedy=False)),
    ('cancel greedy=False pi', lambda: cancel((x**2 - pi)/(x + sqrt(pi)), greedy=False)),
    ('cancel default x+s2', lambda: cancel((x**2 - 2)/(x + s2))),
    ('cancel ext True x+s2', lambda: cancel((x**2 - 2)/(x + s2), extension=True)),
    ('cancel greedy=False gens x', lambda: cancel((x**2 - 2)/(x + s2), x, greedy=False)),
    ('cancel greedy=False x**3-2x', lambda: cancel((x**3 - 2*x)/(x**2 + s2*x), greedy=False)),
    ('cancel greedy=False s3', lambda: cancel((x**2 - 3)/(x - s3), greedy=False)),
    ('cancel greedy=False I', lambda: cancel((x**2 + 1)/(x - I), greedy=False)),
    ('cancel greedy=False 2x+2s2', lambda: cancel((x**2 - 2)/(2*x + 2*s2), greedy=False)),
    ('cancel greedy=False +1', lambda: cancel((x**2 - 2)/(x + s2) + 1, greedy=False)),
    ('cancel greedy=False y', lambda: cancel((x**2 - 2*y**2)/(x + s2*y), greedy=False)),
    ('cancel greedy=True x+s2', lambda: cancel((x**2 - 2)/(x + s2), greedy=True)),
    ('cancel greedy=False s2x+2', lambda: cancel((s2*x + 2)/(x + s2), greedy=False)),
    ('cancel gens x x+s2', lambda: cancel((x**2 - 2)/(x + s2), x)),
    ('cancel greedy=False no gens', lambda: cancel(s2/(1 + s2), greedy=False)),
    ('cancel greedy=False tuple', lambda: cancel(((x**2 - 2), (x + s2)), greedy=False)),
    ('cancel greedy=False cbrt', lambda: cancel((x**3 - 2)/(x - c2), greedy=False)),
    ('cancel greedy=False s2 s3', lambda: cancel((x**2 - 5 - 2*sqrt(6))/(x - s2 - s3), greedy=False)),
    ('gcd greedy=False', lambda: gcd(x**2 - 2, x + s2, greedy=False)),
    ('lcm greedy=False', lambda: lcm(x**2 - 2, x + s2, greedy=False)),
    ('gcd default', lambda: gcd(x**2 - 2, x + s2)),
    ('gcd ext True', lambda: gcd(x**2 - 2, x + s2, extension=True)),
    ('lcm ext True', lambda: lcm(x**2 - 2, x + s2, extension=True)),
    ('gcd greedy=False t', lambda: gcd(x**2 - 2*t**2, x + s2*t, greedy=False)),
    ('cofactors greedy=False', lambda: cofactors(x**2 - 2, x + s2, greedy=False)),
    ('div greedy=False', lambda: div(x**2 - 2, x + s2, greedy=False)),
    ('quo greedy=False', lambda: quo(x**2 - 2, x + s2, greedy=False)),
    ('rem greedy=False', lambda: rem(x**2 - 2, x + s2, greedy=False)),
    ('gcd gens x', lambda: gcd(x**2 - 2, x + s2, x)),
    ('lcm gens x', lambda: lcm(x**2 - 2, x + s2, x)),
    ('P cancel ZZ EX', lambda: Poly(x**2 - 2, x).cancel(Poly(x + s2, x))),
    ('resultant greedy=False', lambda: resultant(x**2 - 2, x + s2, greedy=False)),
    ('gcd greedy=False pi', lambda: gcd(x**2 - pi, x + sqrt(pi), greedy=False)),
    ('gcd greedy=False I', lambda: gcd(x**2 + 1, x + I, greedy=False)),
    ('lcm greedy=False I', lambda: lcm(x**2 + 1, x + I, greedy=False)),
    # roots
    ('real_roots x**2-s2 ext', lambda: real_roots(x**2 - s2, extension=True)),
    ('real_roots P(x**2-s2) ext', lambda: real_roots(Poly(x**2 - s2), extension=True)),
    ('real_roots P(x**2-s2,ext) ', lambda: real_roots(Poly(x**2 - s2, extension=True))),
    ('all_roots x**2+s2 ext', lambda: all_roots(x**2 + s2, extension=True)),
    ('all_roots P(x**2+s2) ext', lambda: all_roots(Poly(x**2 + s2), extension=True)),
    ('all_roots P(x**2+s2,x) ext', lambda: all_roots(Poly(x**2 + s2, x), extension=True)),
    ('real_roots x**2-s2', lambda: real_roots(x**2 - s2)),
    ('all_roots x**2+s2', lambda: all_roots(x**2 + s2)),
    ('all_roots P(x**2+s2,x)', lambda: all_roots(Poly(x**2 + s2, x))),
    ('real_roots P(xy-s2) ext', lambda: real_roots(Poly(x*y - s2), extension=True)),
    ('all_roots P(x**2-2y) ext', lambda: all_roots(Poly(x**2 - 2*y), extension=True)),
    ('nroots x**2-s2', lambda: nroots(x**2 - s2)),
    ('nroots P(x**2-s2)', lambda: nroots(Poly(x**2 - s2))),
    ('nroots P(x**2-s2,x)', lambda: nroots(Poly(x**2 - s2, x))),
    ('nroots P(x**2-s2,ext)', lambda: nroots(Poly(x**2 - s2, extension=True))),
    ('nroots x**2-s2*I', lambda: nroots(x**2 - s2*I)),
    ('real_roots x**3-s2x ext', lambda: real_roots(x**3 - s2*x, extension=True)),
    ('all_roots x**2-I ext', lambda: all_roots(x**2 - I, extension=True)),
    ('all_roots x**2-s2*I ext', lambda: all_roots(x**2 - s2*I, extension=True)),
    ('real_roots x**2-2 mult=False ext', lambda: real_roots(x**2 - 2, multiple=False, extension=True)),
    ('ground_roots expanded ext', lambda: ground_roots(x**2 - 2*s2*x + 2, x, extension=True)),
    ('roots expanded', lambda: roots(x**2 - 2*s2*x + 2, x)),
    ('P real_roots ext', lambda: Poly(x**2 - s2, x, extension=True).real_roots()),
    ('P all_roots ext', lambda: Poly(x**2 - s2, x, extension=True).all_roots()),
    ('P nroots ext', lambda: Poly(x**2 - s2, x, extension=True).nroots()),
    ('real_roots P(x**2-s2,x,EX) ext', lambda: real_roots(Poly(x**2 - s2, x, domain='EX'), extension=True)),
    ('all_roots P(x**4-s2*x**2) ext', lambda: all_roots(Poly(x**4 - s2*x**2), extension=True)),
    ('real_roots P(x**2-s3, gaussian) ext', lambda: real_roots(Poly(x**2 - 2, x, gaussian=True), extension=True)),
    ('real_roots P(x**2-pi) ext', lambda: real_roots(Poly(x**2 - pi), extension=True)),
    # CRootOf extensions
    ('K QQ<r>', lambda: QQ.algebraic_field(r)),
    ('K QQ<r> ext', lambda: QQ.algebraic_field(r).ext),
    ('K QQ<r> mod', lambda: QQ.algebraic_field(r).mod),
    ('K QQ<r2,s2>', lambda: QQ.algebraic_field(r2, s2)),
    ('K QQ<s2,r2> ext', lambda: QQ.algebraic_field(s2, r2).ext),
    ('K QQ<r> from r**3', lambda: QQ.algebraic_field(r).from_sympy(r**3)),
    ('K QQ<r> roundtrip', lambda: (lambda K: K.to_sympy(K.from_sympy(r**2 + 1)))(QQ.algebraic_field(r))),
    ('minpoly r', lambda: minimal_polynomial(r, y)),
    ('minpoly rootof5', lambda: minimal_polynomial(CRootOf(x**5 - x - 1, 0), y)),
    ('minpoly r+1', lambda: minimal_polynomial(r + 1, y)),
    ('minpoly r**2', lambda: minimal_polynomial(r**2, y)),
    ('minpoly r+s2', lambda: minimal_polynomial(r + s2, y)),
    ('minpoly rc', lambda: minimal_polynomial(rc, y)),
    ('minpoly r*rc', lambda: minimal_polynomial(r*rc, y)),
    ('minpoly CRootOf(x**4-2,1)', lambda: minimal_polynomial(CRootOf(x**4 - 2, 1), y)),
    ('minpoly r2 polys', lambda: minimal_polynomial(r2, y, polys=True)),
    ('primitive_element [r,s2]', lambda: primitive_element([r, s2], y)),
    ('primitive_element [r2,s3]', lambda: primitive_element([r2, s3], y)),
    ('primitive_element [r,q2]', lambda: primitive_element([r, CRootOf(x**2 - 2, 0)], y)),
    ('to_number_field r**2 r', lambda: to_number_field(r**2, r)),
    ('to_number_field r r', lambda: to_number_field(r, r)),
    ('to_number_field r+1 r', lambda: to_number_field(r + 1, r)),
    ('to_number_field q2 s2', lambda: to_number_field(q2, s2)),
    ('to_number_field s2 q2', lambda: to_number_field(s2, q2)),
    ('to_number_field r2**2 r2', lambda: to_number_field(r2**2, r2)),
    ('field_isomorphism r r', lambda: field_isomorphism(r, r)),
    ('AlgebraicNumber r minpoly', lambda: AlgebraicNumber(r).minpoly),
    ('AlgebraicNumber r coeffs', lambda: AlgebraicNumber(r, [1, 2]).as_expr()),
    ('Poly r*y+1 ext True', lambda: Poly(r*y + 1, y, extension=True)),
    ('Poly r ext True', lambda: Poly(r, y, extension=True)),
    ('factor y**3+y+3 ext r', lambda: factor(y**3 + y + 3, extension=r)),
    ('factor y**3-2 ext r2', lambda: factor(y**3 - 2, extension=r2)),
    ('fl y**3-2 ext r2', lambda: factor_list(y**3 - 2, extension=r2)),
    ('factor y**2-2 ext q2', lambda: factor(y**2 - 2, extension=q2)),
    ('Poly y-r2 ext True fl', lambda: Poly(y - r2, y, extension=True).factor_list()),
    # exception classes
    ('exc factor x ext s2', lambda: factor(x**2 - 2*t**2, x, extension=s2)),
    ('exc ZZ[t](1/t)', lambda: ZZ[t](1/t)),
    ('exc QQ_I[t](1/t)', lambda: QQ_I[t](1/t)),
    ('exc ZZ_I[t](1/t)', lambda: ZZ_I[t](1/t)),
    ('exc QQ[t](1/t)', lambda: QQ[t](1/t)),
    ('exc QQ<s2>[t](1/t)', lambda: QQ.algebraic_field(s2)[t](1/t)),
    ('exc ZZ[t](s2)', lambda: ZZ[t](s2)),
    ('ZZ[t](t/2)', lambda: ZZ[t](t/2)),
    ('exc gcd x ext s2', lambda: gcd(x**2 - 2*t**2, x - s2*t, x, extension=s2)),
    ('exc Poly x-s2t x ext s2', lambda: Poly(x - s2*t, x, extension=s2)),
    ('exc ZZ[t](1/(t+1))', lambda: ZZ[t](1/(t + 1))),
    ('exc QQ_I[t](I/t)', lambda: QQ_I[t](I/t)),
    ('exc factor x**2-2pi**2 x ext s2 list', lambda: factor_list(x**2 - 2*pi**2, x, extension=s2)),
]


def value(r):
    """The result as a SymPy-parsable expression string, or None."""
    try:
        if (isinstance(r, tuple) and len(r) == 2 and isinstance(r[1], list)
                and all(isinstance(p, tuple) and len(p) == 2 for p in r[1])):
            c = value(r[0])
            parts = [f'({value(f)})**({k})' for f, k in r[1]]
            return None if c is None or 'None' in ''.join(parts) else '*'.join([f'({c})'] + parts)
        if isinstance(r, (tuple, list)):
            parts = [value(v) for v in r]
            return None if None in parts else '(' + ', '.join(parts) + ',)'
        if isinstance(r, dict):
            return value(sorted(((str(k), v) for k, v in r.items()), key=str))
        if hasattr(r, 'as_expr') and callable(r.as_expr):
            return str(r.as_expr())
        if hasattr(r, 'ex'):
            return str(r.ex)
        if isinstance(r, (bool, int, str)):
            return str(r)
        if hasattr(r, 'is_Number') or hasattr(r, 'free_symbols'):
            return str(r)
    except Exception:
        return None
    return None


import signal


class CaseTimeout(Exception):
    pass


def _alarm(*_):
    raise CaseTimeout('no result within %d s' % LIMIT)


LIMIT = 60
signal.signal(signal.SIGALRM, _alarm)
only = sys.argv[3:]
assert len({n for n, _ in cases}) == len(cases), 'duplicate case names'

with open(out, 'w') as fh:
    for name, fn in cases:
        if only and name not in only:
            continue
        signal.alarm(LIMIT)
        val = None
        try:
            t0 = time.perf_counter()
            res = fn(); text = str(res); err = None; val = value(res)
            dt = time.perf_counter() - t0
        except (Exception, BaseException) as e:  # a native panic is a BaseException
            if isinstance(e, KeyboardInterrupt):
                raise
            text = None; err = f'{type(e).__name__}: {str(e)[:200]}'; dt = None
        finally:
            signal.alarm(0)
        fh.write(json.dumps({'case': name, 'result': text, 'value': val, 'error': err, 'seconds': dt}) + '\n')
        fh.flush()
        print(name, '|', err or text[:150], flush=True)
