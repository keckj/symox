"""ratintex: EX log_to_atan recursion, RootSum order, binomial preimages.

Every expected value is computed by SymPy 1.14 in the same process and
compared through its printed form.
"""
import importlib

import pytest
import sympy as sp
import symox as sx

rt = importlib.import_module('symox.integrals.rationaltools')
ref = importlib.import_module('sympy.integrals.rationaltools')
x = sx.Symbol('x')
xp = sp.Symbol('x')
a, ap = sx.Symbol('a'), sp.Symbol('a')


def both(text):
    return sx.sympify(text, locals={'x': x, 'a': a}), sp.sympify(text, locals={'x': xp, 'a': ap})


@pytest.mark.parametrize('f,g,fdom,gdom', [
    ('x + 1/2', 'sqrt(3)/2', 'QQ', 'EX'),
    ('x**2 + sqrt(2)*x', 'x + sqrt(2)', 'EX', 'EX'),
    ('x**2 + sqrt(2)', 'sqrt(2)*x + 1', 'EX', 'EX'),
    ('x**3 + sqrt(3)*x + 1', 'sqrt(2)*x**2 - x', 'EX', 'EX'),
    ('2*x + 1', 'sqrt(5)', 'ZZ', 'EX'),
])
def test_ex_log_to_atan_is_the_reference_recursion(f, g, fdom, gdom):
    (fs, fp), (gs, gp) = both(f), both(g)
    F, G = sx.Poly(fs, x, domain=fdom), sx.Poly(gs, x, domain=gdom)
    expected = ref.log_to_atan(sp.Poly(fp, xp, domain=fdom), sp.Poly(gp, xp, domain=gdom))
    native = sx._core._log_to_atan_ex(F, G)
    assert native is not None
    assert str(native) == str(expected)
    assert str(rt.log_to_atan(F, G)) == str(expected)


def test_rational_log_to_atan_route_is_still_taken():
    F, G = sx.Poly(x**2 - 1, x, domain='QQ'), sx.Poly(2*x + 3, x, domain='QQ')
    assert sx._core._log_to_atan(F, G) is not None
    assert sx._core._log_to_atan_ex(F, G) is None


@pytest.mark.parametrize('source', ['1/(x**5 + 1)', 'x/(x**3 - a*x + 1)', '1/(x**2 + 2)'])
def test_real_false_rootsum_prints_in_reference_order(source):
    s, p = both(source)
    assert str(rt.ratint(s, x, real=False)) == str(ref.ratint(p, xp, real=False))


@pytest.mark.parametrize('source', ['x**4 + 2*x**2 + 3', 'x**4 - x**2 + 3', 'x**8 + x**4 + 3',
                                    '(x**4 - 2*x**2 - 3)**2', 'x**4 + x**2 + 1'])
def test_decomposed_binomial_preimages_are_the_reference_ones(source):
    s, p = both(source)
    assert str(sx.roots(s, x)) == str(sp.roots(p, xp))


@pytest.mark.parametrize('source', ['x**8 - 2*x**4 + x**2 - 3', 'x**6 + x**2 + 1'])
def test_nested_radical_binomials_keep_a_complete_zero_set(source):
    # Documented form gap (nested_binomial): the quadratic formula is kept.
    s, p = both(source)
    found = sx.roots(s, x)
    assert sum(found.values()) == sp.degree(p, xp)
    values = sorted((complex(sx.N(r, 30)) for r in found), key=lambda z: (round(z.real, 9), round(z.imag, 9)))
    expected = sorted((complex(r) for r in sp.Poly(p, xp).nroots(n=30)), key=lambda z: (round(z.real, 9), round(z.imag, 9)))
    assert all(abs(u - v) < 1e-20 for u, v in zip(values, expected))


@pytest.mark.parametrize('text', [
    '-14945472*31**(1/3)*(2*31**(1/3) + 2**(1/3)*(3*sqrt(93) + 31)**(2/3))'
    '*(3*sqrt(93) + 31)**(1/3)*(3*sqrt(3)*31**(5/6) + 31*31**(1/3))',
    '2**(1/3)*(31**(1/3) + 1)*(31**(1/3) + 3**(1/2))',
    '(2**(1/3) + 3**(1/3))*(2**(2/3) + 5**(1/2))*x',
    '(sqrt(2) + sqrt(3))*(sqrt(2) - sqrt(3))*(sqrt(6) + 1)',
])
def test_expand_of_mixed_radicals_follows_the_reference_order(text):
    # Mul.flatten combines numeric radicals by association order, so the
    # reference Mul._expandsums splitting decides 1922**(1/3) versus
    # 2**(1/3)*31**(2/3); the native expand reproduces that order.
    s, p = both(text)
    assert str(sx.expand(s)) == str(sp.expand(p))
    flags = dict(deep=False, mul=True, power_exp=False, power_base=False, multinomial=False, log=False)
    assert str(s.expand(**flags)) == str(p.expand(**flags))


@pytest.mark.parametrize('text', [
    '92163364291974*31**(2/3) + 296263766063322*sqrt(3)*31**(1/6)',
    '31**(2/3) + sqrt(3)*31**(1/6)',
    '1922**(1/3) + 2**(1/3)*31**(2/3)',
    '31**(2/3) + sqrt(3)',
    '5**(2/3) + sqrt(3)*7**(1/6)',
    '31**(2/3) + 3**(1/3)*31**(1/6)',
    # keys one unit apart once each factor is correctly rounded
    '395018354751096*sqrt(3)*31**(1/6)/D + 122884485722632*31**(2/3)/D'.replace('D',
        '(4613477382612*62**(1/3)*(3*sqrt(93) + 31)**(2/3)'
        ' + 478395319050*2**(1/3)*sqrt(3)*31**(5/6)*(3*sqrt(93) + 31)**(2/3)'
        ' + 14335343519193*2**(2/3)*sqrt(93)*(3*sqrt(93) + 31)**(1/3)'
        ' + 138245046437961*2**(2/3)*(3*sqrt(93) + 31)**(1/3)'
        ' + 92163364291974*31**(2/3) + 296263766063322*sqrt(3)*31**(1/6))'),
])
def test_print_order_of_numeric_radicals_is_the_reference_one(text):
    # as_ordered_terms folds the numeric factors into a complex float in
    # the canonical Mul.make_args order (the last bit decides) and keeps the
    # canonical Add.make_args order between equal keys.
    s, p = both(text)
    assert str(s) == str(p)
    assert [str(t) for t in s.as_ordered_terms()] == [str(t) for t in p.as_ordered_terms()]
