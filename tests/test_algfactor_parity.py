"""Factorization and cancellation with algebraic numbers (wave 45, lot algfactor).

Natural factor lists whose generators include algebraic numbers (the
reconstruction certificate of ``checked_factors::verify_expression_factors``
and the component-wise square-free combination), ``greedy=False`` for
``gcd``/``lcm``/``cofactors``/``cancel`` (``convert::nongreedy_generators``
and the ``EX`` polynomials), the conversion errors of polynomial rings and
algebraic fields, and the ``extension`` option of ``real_roots`` and
``all_roots``. Every case runs live on SymPy 1.14 and on symox and compares
the printed results, or the exception classes.
"""
import pytest

import sympy
import symox


def env(M):
    x, y, t = M.symbols('x y t')
    return dict(M=M, x=x, y=y, t=t, s2=M.sqrt(2), s3=M.sqrt(3), pi=M.pi, I=M.I)


CASES = {
    # natural factor lists (generators x and sqrt(2))
    'fl (x-s2)(x+s2)': lambda M, x, s2, **_: M.factor_list((x - s2)*(x + s2)),
    'fl (x-s2*t)**2': lambda M, x, s2, t, **_: M.factor_list((x - s2*t)**2),
    'fl (x-s2*pi)(x+s2*pi)': lambda M, x, s2, pi, **_: M.factor_list((x - s2*pi)*(x + s2*pi)),
    'fl (x+1)(x-s2)**2': lambda M, x, s2, **_: M.factor_list((x + 1)*(x - s2)**2),
    'fl (x-s2)(y+s3)': lambda M, x, y, s2, s3, **_: M.factor_list((x - s2)*(y + s3)),
    'fl (x-I)(x+I)': lambda M, x, I, **_: M.factor_list((x - I)*(x + I)),
    'fl 3(x-s2)**2': lambda M, x, s2, **_: M.factor_list(3*(x - s2)**2),
    'fl -(x-s2)(x+s2)': lambda M, x, s2, **_: M.factor_list(-(x - s2)*(x + s2)),
    'fl (2x-s2)(x+s2)': lambda M, x, s2, **_: M.factor_list((2*x - s2)*(x + s2)),
    'fl (x-s2*y)**2(x+y)': lambda M, x, y, s2, **_: M.factor_list((x - s2*y)**2*(x + y)),
    'fl frac': lambda M, x, s2, **_: M.factor_list((x - s2)/(x + 1)),
    'fl (x-s2)(x+s2)(x-pi)': lambda M, x, s2, pi, **_: M.factor_list((x - s2)*(x + s2)*(x - pi)),
    'fl (x-s2)(x-pi)': lambda M, x, s2, pi, **_: M.factor_list((x - s2)*(x - pi)),
    'fl frac=True s2': lambda M, x, s2, **_: M.factor_list((x - s2)**2*(x + 1)/(x + s2), frac=True),
    'sqfl frac=True s2': lambda M, x, s2, **_: M.sqf_list((x - s2)**2*(x + 1)/(x + s2), frac=True),
    'sqfl frac=True rational': lambda M, x, **_: M.sqf_list((x**2 - 1)*(x + 1)/((x + 2)*(x**2 - 4)), frac=True),
    'sqfl frac s2 denominator': lambda M, x, s2, **_: M.sqf_list((x**2 - 1)*(x + 1)/(x - s2), frac=True),
    'sqfl (x**2-1)**2(x+1)': lambda M, x, **_: M.sqf_list((x**2 - 1)**2*(x + 1)),
    'fl (1-s2)(1+s2)': lambda M, s2, **_: M.factor_list((1 - s2)*(1 + s2)),
    'sqfl (1-s2)(1+s2)': lambda M, s2, **_: M.sqf_list((1 - s2)*(1 + s2)),
    'fl (1-s2)(1+s2) polys': lambda M, s2, **_: M.factor_list((1 - s2)*(1 + s2), polys=True),
    'fl (1-s2)(1+s2) x': lambda M, x, s2, **_: M.factor_list((1 - s2)*(1 + s2), x),
    'fl s2(1-s2)': lambda M, s2, **_: M.factor_list(s2*(1 - s2)),
    'fl (1-s3)**2(1+s3)**2': lambda M, s3, **_: M.factor_list((1 - s3)**2*(1 + s3)**2),
    'factor (x-s2)(x+s2)': lambda M, x, s2, **_: M.factor((x - s2)*(x + s2)),
    'sqfl (x+1)(x-s2)**2': lambda M, x, s2, **_: M.sqf_list((x + 1)*(x - s2)**2),
    'sqfl (x-s2*t)**2': lambda M, x, s2, t, **_: M.sqf_list((x - s2*t)**2),
    'sqfl (x**2-1)(x+1)': lambda M, x, **_: M.sqf_list((x**2 - 1)*(x + 1)),
    'sqfl (x**2-y**2)(x+y)': lambda M, x, y, **_: M.sqf_list((x**2 - y**2)*(x + y)),
    'sqfl 3(x-s2)**2(x+1)**2': lambda M, x, s2, **_: M.sqf_list(3*(x - s2)**2*(x + 1)**2),
    # extensions
    'factor ext s2': lambda M, x, s2, **_: M.factor(x**2 - 2, extension=s2),
    'factor ext s2 multivariate': lambda M, x, t, s2, **_: M.factor(x**2 - 2*t**2, extension=s2),
    'factor gaussian': lambda M, x, t, **_: M.factor(x**2 + t**2, gaussian=True),
    'factor ext True': lambda M, x, s2, s3, **_: M.factor(x**2 - s2*x - s3*x + s2*s3, extension=True),
    'factor_list ext s2 multivariate': lambda M, x, t, s2, **_: M.factor_list(x**2 - 2*t**2, extension=s2),
    # greedy=False
    'cancel greedy=False': lambda M, x, s2, **_: M.cancel((x**2 - 2)/(x + s2), greedy=False),
    'cancel greedy=False t': lambda M, x, t, s2, **_: M.cancel((x**2 - 2*t**2)/(x + s2*t), greedy=False),
    'cancel greedy=False pi': lambda M, x, pi, **_: M.cancel((x**2 - pi)/(x + M.sqrt(pi)), greedy=False),
    'cancel greedy=False constant': lambda M, s2, **_: M.cancel(s2/(1 + s2), greedy=False),
    'cancel greedy=False tuple': lambda M, x, s2, **_: M.cancel((x**2 - 2, x + s2), greedy=False),
    'cancel greedy=False s2 s3': lambda M, x, s2, s3, **_: M.cancel((x**2 - 5 - 2*M.sqrt(6))/(x - s2 - s3), greedy=False),
    'cancel default': lambda M, x, s2, **_: M.cancel((x**2 - 2)/(x + s2)),
    'gcd greedy=False': lambda M, x, s2, **_: M.gcd(x**2 - 2, x + s2, greedy=False),
    'gcd greedy=False t': lambda M, x, t, s2, **_: M.gcd(x**2 - 2*t**2, x - s2*t, greedy=False),
    'lcm greedy=False': lambda M, x, s2, **_: M.lcm(x**2 - 2, x + s2, greedy=False),
    'cofactors greedy=False': lambda M, x, s2, **_: M.cofactors(x**2 - 2, x + s2, greedy=False),
    'gcd greedy=False rational': lambda M, x, **_: M.gcd(x**2 - 1, x + 1, greedy=False),
    # roots with extension
    'real_roots Poly ext': lambda M, x, s2, **_: M.real_roots(M.Poly(x**2 - s2), extension=True),
    'all_roots Poly ext': lambda M, x, s2, **_: M.all_roots(M.Poly(x**2 + s2), extension=True),
    'real_roots expr ext': lambda M, x, s2, s3, **_: M.real_roots(M.expand((x - s2)*(x - s3)), extension=True),
    'real_roots Poly xy ext': lambda M, x, y, s2, **_: M.real_roots(M.Poly(x*y - s2), extension=True),
    'all_roots Poly x**2-2y ext': lambda M, x, y, **_: M.all_roots(M.Poly(x**2 - 2*y), extension=True),
    # CRootOf extensions
    'QQ<CRootOf>': lambda M, x, **_: M.QQ.algebraic_field(M.CRootOf(x**3 + x + 3, 0)),
    'minpoly CRootOf+s2': lambda M, x, y, s2, **_: M.minimal_polynomial(M.CRootOf(x**3 + x + 3, 0) + s2, y),
    'primitive_element CRootOf s2': lambda M, x, y, s2, **_: M.primitive_element([M.CRootOf(x**3 + x + 3, 0), s2], y),
    'to_number_field CRootOf': lambda M, x, **_: M.to_number_field(M.CRootOf(x**3 + x + 3, 0)**2, M.CRootOf(x**3 + x + 3, 0)),
    # conversion errors
    'factor x ext s2 symbolic': lambda M, x, t, s2, **_: M.factor(x**2 - 2*t**2, x, extension=s2),
    'factor x ext s2 pi': lambda M, x, pi, s2, **_: M.factor(x**2 - 2*pi**2, x, extension=s2),
    'gcd x ext s2': lambda M, x, t, s2, **_: M.gcd(x**2 - 2*t**2, x - s2*t, x, extension=s2),
    'Poly x ext s2 sqrt3': lambda M, x, s3, s2, **_: M.Poly(x - s3, x, extension=s2),
    'ZZ[t](1/t)': lambda M, t, **_: M.ZZ[t](1/t),
    'QQ[t](1/t)': lambda M, t, **_: M.QQ[t](1/t),
    'ZZ[t](sqrt2)': lambda M, t, s2, **_: M.ZZ[t](s2),
    'ZZ[t](x)': lambda M, t, x, **_: M.ZZ[t](x),
    'ZZ[t](t/2)': lambda M, t, **_: M.ZZ[t](t/2),
    'ZZ[t].from_sympy(1/t)': lambda M, t, **_: M.ZZ[t].from_sympy(1/t),
    'ZZ[t].convert(1/t)': lambda M, t, **_: M.ZZ[t].convert(1/t),
    'QQ_I[t].convert(1/t)': lambda M, t, **_: M.QQ_I[t].convert(1/t),
    'QQ_I[t](1/t)': lambda M, t, **_: M.QQ_I[t](1/t),
    'ZZ(t)(sqrt2)': lambda M, t, s2, **_: M.ZZ.frac_field(t)(s2),
    'ZZ[t](2t+1)': lambda M, t, **_: M.ZZ[t](2*t + 1),
    'ZZ[t](1.5)': lambda M, t, **_: M.ZZ[t](1.5),
}


def run(M, name):
    try:
        return str(CASES[name](**env(M)))
    except Exception as e:  # compare the exception class only
        return type(e).__name__


@pytest.mark.parametrize('name', sorted(CASES))
def test_against_sympy(name):
    assert run(symox, name) == run(sympy, name)


def test_constant_cancel_expands_like_reference():
    # sring finds no generator: the reference returns ``f.expand()``.
    def go(M):
        s2 = M.sqrt(2)
        return str(M.cancel(s2*(1 + s2), greedy=False))
    assert go(symox) == go(sympy)
