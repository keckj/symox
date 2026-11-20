"""Ground-domain and list polynomial contracts, compared with SymPy 1.14."""

import pytest
import sympy
import symox


def cases(M):
    x, y = M.symbols('x y')
    return {
        'cofactors integers': lambda: M.cofactors(8, 6),
        'cofactors rational': lambda: M.cofactors(M.Rational(2, 3), M.Rational(4, 9)),
        'cofactors negative': lambda: M.cofactors(-8, 6),
        'cofactors zero': lambda: M.cofactors(0, -6),
        'cofactors both zero error': lambda: M.cofactors(0, 0),
        'cofactors polynomial zero': lambda: M.cofactors(0, 0, x),
        'cofactors floats': lambda: M.cofactors(3., 9.),
        'cofactors float zero': lambda: M.cofactors(0., 0.),
        'cofactors explicit QQ': lambda: M.cofactors(8, 6, x, domain='QQ'),
        'cofactors explicit QQ polys': lambda: M.cofactors(8, 6, x, domain='QQ', polys=True),
        'cofactors polynomial float': lambda: M.cofactors(3., 9., x),
        'gcd floats': lambda: M.gcd(3., 9.),
        'gcd float and integer': lambda: M.gcd(3., 9),
        'gcd float zero': lambda: M.gcd(0., 0.),
        'gcd polynomial float zero': lambda: M.gcd(0., 0., x),
        'gcd explicit QQ': lambda: M.gcd(8, 6, x, domain='QQ'),
        'gcd explicit QQ polys': lambda: M.gcd(8, 6, x, domain='QQ', polys=True),
        'gcd rational polynomial': lambda: M.gcd(M.Rational(1, 2), M.Rational(1, 3), x),
        'lcm floats': lambda: M.lcm(3., 9.),
        'lcm signed floats': lambda: M.lcm(-3., 9.),
        'lcm float zero': lambda: M.lcm(0., 9.),
        'lcm polynomial floats': lambda: M.lcm(3., 9., x),
        'lcm explicit QQ': lambda: M.lcm(8, 6, x, domain='QQ'),
        'lcm rational polynomial': lambda: M.lcm(M.Rational(1, 2), M.Rational(1, 3), x),
        'gcdex auto false inferred': lambda: M.gcdex(x + 1, x - 1, auto=False),
        'gcdex auto false explicit': lambda: M.gcdex(x + 1, x - 1, x, auto=False),
        'gcdex auto false QQ inferred': lambda: M.gcdex(x/2 + 1, x - 1, auto=False),
        'gcdex domain QQ inferred': lambda: M.gcdex(x + 1, x - 1, domain='QQ'),
        'gcdex domain ZZ inferred': lambda: M.gcdex(x + 1, x - 1, domain='ZZ'),
        'gcdex polys inferred': lambda: M.gcdex(x + 1, x - 1, polys=True),
        'half_gcdex auto false inferred': lambda: M.half_gcdex(x + 1, x - 1, auto=False),
        'half_gcdex QQ inferred': lambda: M.half_gcdex(x/2 + 1, x - 1, auto=False),
        'half_gcdex polys inferred': lambda: M.half_gcdex(x + 1, x - 1, polys=True),
        **{
            f'{name} {label}': (lambda name=name, args=args, kw=kw:
                getattr(M, name)(*args, **kw))
            for name in ('gcd_list', 'lcm_list')
            for label, args, kw in [
                ('empty', ([],), {}),
                ('singleton signed', ([-3],), {}),
                ('integers', ([8, 6, 12],), {}),
                ('rational', ([M.Rational(2, 3), M.Rational(4, 9)],), {}),
                ('floats', ([3., 9.],), {}),
                ('float zeros', ([0., 0.],), {}),
                ('polynomial', ([2*x, 3*x, x**2],), {}),
                ('multivariate', ([x*y, x*y**2, x**2*y],), {}),
                ('polys', ([2*x, 3*x],), {'polys': True}),
                ('domain QQ', ([2*x, 3*x],), {'domain': 'QQ'}),
                ('domain QQ object', ([2*x, 3*x],), {'domain': M.QQ, 'polys': True}),
                ('field', ([2*x, 4*x],), {'field': True}),
                ('field polys', ([2*x, 4*x],), {'field': True, 'polys': True}),
                ('unknown option', ([x, x**2],), {'unknown': 3}),
                ('conflicting options', ([x, x**2],), {'field': True, 'domain': 'ZZ'}),
                ('domain QQ polys', ([2*x, 3*x],), {'domain': 'QQ', 'polys': True}),
                ('numbers explicit x', ([8, 6], x), {}),
                ('numbers explicit x QQ', ([8, 6], x), {'domain': 'QQ'}),
                ('empty explicit x', ([], x), {'polys': True}),
                ('empty no gens error', ([],), {'polys': True}),
                ('numbers no gens error', ([8, 6],), {'polys': False}),
                ('bad flag', ([x, x**2],), {'auto': False}),
                ('modulus', ([x**2 - 1, x**2 + x + 1], x), {'modulus': 3}),
                ('modulus polys', ([x**2 - 1, x**2 + x + 1], x), {'modulus': 3, 'polys': True}),
                ('Poly generators QQ same', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x, domain='QQ')], x), {'polys': True}),
                ('Poly generators QQ replacement', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x, domain='QQ')], y), {'polys': True}),
                ('Poly generators QQ addition', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x, domain='QQ')], y, x), {'polys': True}),
                ('Poly generators QQ permutation', ([M.Poly(x*y, x, y, domain='QQ'), M.Poly(x*y**2, x, y, domain='QQ')], y, x), {'polys': True}),
                ('Poly generators QQ domain error', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x, domain='QQ')], y), {'domain': 'QQ', 'polys': True}),
                ('Poly generators GF same', ([M.Poly(x**2 + 1, x, modulus=2), M.Poly(x + 1, x, modulus=2)], x), {'polys': True}),
                ('Poly generators GF replacement', ([M.Poly(x**2 + 1, x, modulus=2), M.Poly(x + 1, x, modulus=2)], y), {'polys': True}),
                ('Poly generators GF addition', ([M.Poly(x**2 + 1, x, modulus=2), M.Poly(x + 1, x, modulus=2)], y, x), {'polys': True}),
                ('Poly generators QQ triple', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x, domain='QQ'), M.Poly(6*x, x, domain='QQ')],), {'polys': True}),
                ('Poly generators GF triple', ([M.Poly(x**2 + 1, x, modulus=2), M.Poly(x + 1, x, modulus=2), M.Poly(x + 1, x, modulus=2)],), {'polys': True}),
                ('Poly GF', ([M.Poly(x**2 + 1, x, modulus=2), M.Poly(x + 1, x, modulus=2)],), {'polys': True}),
                ('Poly mixed GF error', ([M.Poly(x**2 + 1, x, modulus=2), M.Poly(x + 1, x, modulus=3)],), {'polys': True}),
                ('Poly retains QQ', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x)],), {'polys': True}),
                ('Poly expressions', ([M.Poly(2*x, x, domain='QQ'), M.Poly(4*x, x)],), {'polys': False}),
                ('extension polys', ([x**2 - 2, x - M.sqrt(2)], x), {'extension': M.sqrt(2), 'polys': True}),
            ]
        },
    }


def snapshot(value):
    if isinstance(value, tuple):
        return tuple(snapshot(v) for v in value)
    if isinstance(value, (sympy.Poly, symox.Poly)):
        return ('Poly', str(value.as_expr()), tuple(map(str, value.gens)), str(value.domain))
    return (str(value), getattr(value, '_prec', None) if getattr(value, 'is_Float', False) else None)


def run(M, name):
    try:
        return ('value', snapshot(cases(M)[name]()))
    except Exception as exc:
        return ('error', type(exc).__name__)


@pytest.mark.parametrize('name', sorted(cases(sympy)))
def test_reference_contract(name):
    assert run(symox, name) == run(sympy, name)


@pytest.mark.parametrize('name', ['gcd', 'lcm', 'cofactors', 'gcd_list', 'lcm_list'])
def test_float_ground_precision(name):
    def evaluate(M):
        # An input beyond binary64's exponent range and a rational that cannot
        # be represented at 53 bits exercise conversion at the joint precision.
        a = M.Float('1.23456789012345678901234567890123456789e400', 80)
        b = M.Rational(2**270 + 1, 2**269 - 1)
        result = getattr(M, name)([a, b]) if name.endswith('_list') else getattr(M, name)(a, b)
        return result
    actual, expected = evaluate(symox), evaluate(sympy)
    for a, b in zip(actual if isinstance(actual, tuple) else (actual,), expected if isinstance(expected, tuple) else (expected,)):
        assert a._prec == b._prec
        assert str(a) == str(b)


def test_cofactor_reconstruction_and_bezout():
    x = symox.Symbol('x')
    for a, b in [(8, 6), (symox.Rational(2, 3), symox.Rational(4, 9)), (x**3 - x, x**2 - 1)]:
        h, f, g = symox.cofactors(a, b)
        assert symox.expand(h*f - a) == 0
        assert symox.expand(h*g - b) == 0
    a, b = x/2 + 1, x - 1
    s, t, h = symox.gcdex(a, b, auto=False)
    assert symox.expand(s*a + t*b - h) == 0


@pytest.mark.parametrize('name', ['gcd_list', 'lcm_list'])
def test_iterator_admission(name):
    x = symox.Symbol('x')
    assert getattr(symox, name)((p for p in [x, x**2])) == (x if name == 'gcd_list' else x**2)


@pytest.mark.parametrize('precision', [1, 2, 52, 54, 80, 100])
@pytest.mark.parametrize('name', ['gcd', 'lcm', 'cofactors', 'gcd_list', 'lcm_list'])
def test_realfield_output_precision(name, precision):
    def evaluate(M):
        value = M.Float(2, precision=precision)
        return getattr(M, name)([value, 3]) if name.endswith('_list') else getattr(M, name)(value, 3)
    assert snapshot(evaluate(symox)) == snapshot(evaluate(sympy))
