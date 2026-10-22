"""(-2)**(1/6) and other negative bases with rational exponents, compared with SymPy.

SymPy keeps (-2)**(1/6) as one power and joins (-1)**t with a positive base of
the same exponent inside a product. With extension=True this gives QQ<(-2)**(1/6)>
instead of the degree twelve field of 2**(1/6) + (-1)**(1/6).
"""
import pytest
import sympy as sp
import symox as sm

PROBES = [
    '(-2)**R(1,6)',
    '(-2)**R(5,6)',
    '(-2)**R(7,6)',
    '(-2)**R(1,3)',
    '(-2)**R(2,3)',
    '(-2)**R(-1,6)',
    '(-6)**R(1,6)',
    '(-4)**R(1,6)',
    '(-8)**R(1,5)',
    '(-3)**R(1,4)',
    'R(-2,3)**R(1,6)',
    'R(-1,3)**R(1,3)',
    'R(-2,3)**R(-1,3)',
    '(-2)**R(1,6)*(-2)**R(1,6)',
    '(-2)**R(1,6)*2**R(1,6)',
    '(-2)**R(1,3)*3**R(1,3)',
    '((-2)**R(1,6))**3',
    '((-2)**R(1,6))**6',
    '((-2)**R(1,6))**2',
    'sqrt((-2)**R(1,3))',
    '(-2)**R(1,6)*(-1)**R(1,6)',
    '-(-2)**R(1,6)',
    '1/(-2)**R(1,6)',
    '(-2)**R(1,6)+2**R(1,6)*(-1)**R(1,6)',
    '(-2)**R(1,3)*(-3)**R(1,3)',
    '(-2)**R(1,3)*(-2)**R(2,3)',
    'N((-2)**R(1,6))',
    '((-2)**R(1,6)).is_real',
    '((-2)**R(1,6)).is_positive',
    '((-2)**R(1,6)).is_extended_real',
    '((-2)**R(1,6)).is_algebraic',
    '((-2)**R(1,6)).is_zero',
    '((-2)**R(1,3)).is_real',
    'expand(((-2)**R(1,6)+1)**2)',
    '((-2)**R(1,6)).as_base_exp()',
    '((-2)**R(1,6)).conjugate()',
    'Abs((-2)**R(1,6))',
    'arg((-2)**R(1,6))',
    're((-2)**R(1,6))',
    'simplify((-2)**R(1,6))',
    'powsimp((-1)**R(1,6)*2**R(1,6))',
    'sqrt(-2)',
    '(-2)**R(3,2)',
    '(-2)**R(-3,2)',
    '(-9)**R(1,2)',
    '(-27)**R(1,3)',
    'R(-8,27)**R(1,3)',
    '(-10)**R(2,3)',
    '(-1)**R(1,6)*2**R(1,6)',
    'Mul((-1)**R(1,6),2**R(1,6))',
    '(-2)**R(1,6)*x',
    'srepr((-2)**R(1,6)*x)',
    'Poly(x-a6, x, extension=True)',
    'Poly(x-a6, x, extension=True).domain.ext',
    'Poly(x-a6, x, extension=True).domain.mod',
    'Poly(x**2-a6, x, extension=True)',
    'minimal_polynomial(a6, x)',
    'minimal_polynomial(a6+1, x)',
    'minimal_polynomial((-2)**R(1,3), x)',
    'minimal_polynomial((-3)**R(1,4), x)',
    'minimal_polynomial((-6)**R(1,6), x)',
    'factor(x**6+2, extension=True)',
    'factor(x**3+2, extension=(-2)**R(1,3))',
    'factor(x**2-a6**2, extension=True)',
    'factor(x**2-(-2)**R(1,3), extension=True)',
    'factor((x-a6)*(x+1), extension=True)',
    'gcd(x**6+2, x-a6, extension=True)',
    'gcd(x**6+2, x**2-(-2)**R(1,3), extension=True)',
    'primitive_element([a6], x)',
    'primitive_element([(-2)**R(1,3), sqrt(3)], x)',
    'field_isomorphism(a6, a6)',
    'to_number_field(a6, a6)',
    'sqf_list(x**2 - a6, extension=True)',
    'roots(x**6+2, x)',
    'Poly(x**6+2, x).all_roots()',
    'QQ.algebraic_field(a6)',
    'QQ.algebraic_field(a6).ext',
    'factor(x**4 + 3, extension=(-3)**R(1,4))',
    'cancel((x**2-(-2)**R(1,3))/(x-a6), extension=True)',
]


def _namespace(module):
    ns = {}
    exec('from %s import *' % module.__name__, ns)
    exec('from %s.polys.domains import ZZ, QQ, RR, CC, EX\nR=Rational\n'
         'x,y,z,t=symbols("x y z t")\na6=(-2)**Rational(1,6)' % module.__name__, ns)
    return ns


SP = _namespace(sp)
SM = _namespace(sm)


def _show(ns, probe):
    try:
        value = eval(probe, ns)
    except Exception as error:
        return 'EXC ' + type(error).__name__
    text = str(value)
    if hasattr(value, 'domain'):
        text += ' || ' + repr(value.domain)
    return text


@pytest.mark.parametrize('probe', PROBES)
def test_negative_base_power_matches_reference(probe):
    assert _show(SM, probe) == _show(SP, probe)


@pytest.mark.parametrize('items', [['a6'], ['a6', '1'], ['(-2)**R(1,3)'], ['a6', 'sqrt(3)']])
def test_construct_domain_field_matches_reference(items):
    def run(ns):
        return eval('construct_domain([%s], extension=True)' % ', '.join(items), ns)
    (dm, elems), (dp, elemp) = run(SM), run(SP)
    assert str(dm) == str(dp)
    assert [dm.to_sympy(e) for e in elems] == [sm.sympify(str(dp.to_sympy(e))) for e in elemp]


def test_extension_field_of_minus_two_sixth_root():
    x = sm.Symbol('x')
    a = (-2)**sm.Rational(1, 6)
    assert str(a) == '(-2)**(1/6)' and a.as_base_exp() == (-2, sm.Rational(1, 6))
    p = sm.Poly(x - a, x, extension=True)
    assert str(p.domain) == 'QQ<(-2)**(1/6)>'
    assert sm.minimal_polynomial(a, x) == x**6 + 2
    assert (-1)**sm.Rational(1, 6)*2**sm.Rational(1, 6) == a
    assert sm.sqrt(sm.sqrt(2)/2) == 2**sm.Rational(3, 4)/2


def test_nroots_intentional_corrections_are_more_accurate():
    # The reference raises NoConvergence or loses digits on these inputs.
    x = sm.Symbol('x')
    R = sm.Rational
    assert [str(r) for r in sm.nroots((x - R(1, 3))**3)] == ['0.333333333333333'] * 3
    assert [str(r) for r in sm.nroots(x**2 + 2*sm.I*x - 1)] == ['-1.0*I', '-1.0*I']
    assert [str(r) for r in sm.nroots(x**2 - 1e20*x + 1)] == ['1.00000000000000e-20', '1.00000000000000e+20']
    assert [str(r) for r in sm.nroots(x**2 - 2*x + 1 + R(1, 10**10))] == ['1.0 - 1.0e-5*I', '1.0 + 1.0e-5*I']
    double = sm.nroots(x**4 - sm.sqrt(2)*x**2 + R(1, 2))
    assert [str(r) for r in double] == ['-0.840896415253715'] * 2 + ['0.840896415253715'] * 2
