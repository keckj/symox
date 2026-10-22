"""Core gaps closed against the reference SymPy 1.14.

The expected strings are the outputs of SymPy 1.14 on the same inputs,
except where a test documents an intentional difference.
"""
import pickle
import time

import pytest

import symox as sx
from symox import (DiracDelta, Heaviside, I, Integer, Piecewise, Rational, S, Symbol,
                   ccode, exp, expand, im, integrate, laplace_transform, latex, oo,
                   pretty, re, srepr, symbols)

x, y, t, s, w = symbols('x y t s w')
a, b = symbols('a b', real=True)


# --- expand(complex=True), as_real_imag, re/im: no default size limit ---

def _terms(e):
    return len(sx.Add.make_args(e))


@pytest.mark.parametrize('build, expected', [
    (lambda: expand((a + I*b)**200, complex=True), 201),
    (lambda: expand((a + I*b)**129, complex=True), 130),
    (lambda: expand((a + I*b)**-130, complex=True), 131),
    (lambda: expand((a + b + I)**70, complex=True), 2556),
    (lambda: sx.expand_complex((a + I)**200), 102),
    (lambda: ((a + I*b)**150).as_real_imag()[0], 76),
    (lambda: re((a + I*b)**140), 71),
    (lambda: im((a + I*b)**140), 70),
    (lambda: expand(sx.Add(*[(a + I*b)**k for k in range(1, 60)]), complex=True), 1829),
])
def test_complex_expansion_has_no_default_limit(build, expected):
    start = time.perf_counter()
    assert _terms(build()) == expected
    # SymPy 1.14 takes 0.3 to 2.2 s on each case
    assert time.perf_counter() - start < 10


def test_complex_expansion_of_large_numeric_powers():
    assert expand((1 + I)**200, complex=True) == Integer(2)**100
    value = expand((2 + 3*I)**300, complex=True)
    assert value.is_Add and sx.expand((2 + 3*I)**300) == value


def test_as_real_imag_large_power_of_complex_symbol():
    real, imag = ((x + I)**130).as_real_imag()
    assert _terms(expand(imag)) == 4290


def test_re_of_a_large_sum():
    symbols_ = [Symbol('s%d' % k) for k in range(2000)]
    value = re(sx.Add(*[symbols_[k]*I**k for k in range(2000)]))
    assert _terms(value) == 2000


# --- Heaviside(x, H0): a two-argument node ---

HEAVISIDE_CASES = [
    (lambda: Heaviside(x, 0), 'Heaviside(x, 0)', "Heaviside(Symbol('x'), Integer(0))"),
    (lambda: Heaviside(x, 1), 'Heaviside(x, 1)', "Heaviside(Symbol('x'), Integer(1))"),
    (lambda: Heaviside(x, Rational(1, 3)), 'Heaviside(x, 1/3)', None),
    (lambda: Heaviside(x, S.Half), 'Heaviside(x)', "Heaviside(Symbol('x'))"),
    (lambda: Heaviside(x, 0.5), 'Heaviside(x, 0.5)', None),
    (lambda: Heaviside(0, 0.5), '0.500000000000000', None),
    (lambda: Heaviside(x, y), 'Heaviside(x, y)', "Heaviside(Symbol('x'), Symbol('y'))"),
    (lambda: Heaviside(0, 0), '0', 'Integer(0)'),
    (lambda: Heaviside(0, y), 'y', "Symbol('y')"),
    (lambda: Heaviside(0), '1/2', 'Rational(1, 2)'),
    (lambda: Heaviside(-2, 0), '0', None),
    (lambda: Heaviside(3, 0), '1', None),
    (lambda: Heaviside(Symbol('z', zero=True), 3), '3', None),
    (lambda: Heaviside(Symbol('p', positive=True), 0), '1', None),
    (lambda: Heaviside(-Symbol('p', positive=True), 7), '0', None),
    (lambda: Heaviside(x, 0, evaluate=False), 'Heaviside(x, 0)', None),
    (lambda: Heaviside(0, 0, evaluate=False), 'Heaviside(0, 0)', "Heaviside(Integer(0), Integer(0))"),
    (lambda: Heaviside(x, 1).rewrite(Piecewise), 'Piecewise((0, x < 0), (1, True))', None),
    (lambda: Heaviside(x, 0).rewrite(Piecewise), 'Piecewise((0, x <= 0), (1, True))', None),
    (lambda: Heaviside(x, y).rewrite(Piecewise), 'Piecewise((0, x < 0), (y, Eq(x, 0)), (1, True))', None),
    (lambda: Heaviside(x, 0).diff(x), 'DiracDelta(x)', None),
    (lambda: Heaviside(x, 0).diff(x, 2), 'DiracDelta(x, 1)', None),
    (lambda: Heaviside(x, y).diff(x), 'DiracDelta(x)', None),
    (lambda: Heaviside(x, 0).subs(x, 0), '0', None),
    (lambda: Heaviside(x - 1, 0).subs(x, 1), '0', None),
    (lambda: Heaviside(x, y).subs(y, 1), 'Heaviside(x, 1)', None),
    (lambda: Heaviside(x, y).subs(x, 0), 'y', None),
    (lambda: Heaviside(x - 1, y).subs({x: 1, y: 3}), '3', None),
    (lambda: Heaviside(x, 0).subs(x, y), 'Heaviside(y, 0)', None),
    (lambda: Heaviside(x, 0).xreplace({x: 0}), '0', None),
    (lambda: Heaviside(x, 0).func(*Heaviside(x, 0).args), 'Heaviside(x, 0)', None),
    (lambda: 2*Heaviside(x, 0) + Heaviside(x), '2*Heaviside(x, 0) + Heaviside(x)', None),
    (lambda: Heaviside(x, 0)**2, 'Heaviside(x, 0)**2', None),
    (lambda: expand(Heaviside(x, 0)*(x + 1)), 'x*Heaviside(x, 0) + Heaviside(x, 0)', None),
    (lambda: Heaviside(x, 0).as_real_imag(), '(Heaviside(x, 0), 0)', None),
    (lambda: sx.conjugate(Heaviside(x, 0)), 'Heaviside(x, 0)', None),
    (lambda: integrate(Heaviside(x, 0), (x, -1, 1)), '1', None),
    (lambda: integrate(Heaviside(x - 1, 1)*exp(-x), (x, 0, oo)), 'exp(-1)', None),
    (lambda: integrate(Heaviside(x, 1)*x, (x, -1, 2)), '2', None),
    (lambda: laplace_transform(Heaviside(t, 1), t, s, noconds=True), '1/s', None),
    (lambda: sx.fourier_transform(Heaviside(x, 0)*exp(-x), x, w), '1/(2*I*pi*w + 1)', None),
]


@pytest.mark.parametrize('build, text, tree', HEAVISIDE_CASES)
def test_two_argument_heaviside(build, text, tree):
    value = build()
    assert str(value) == text
    if tree is not None:
        assert srepr(value) == tree


def test_heaviside_args_and_identity():
    assert Heaviside(x, 0).args == (x, Integer(0))
    assert Heaviside(x).args == (x, S.Half)
    assert Heaviside(x, 0).pargs == (x, Integer(0))
    assert Heaviside(x).pargs == (x,)
    assert Heaviside(x, 0) != Heaviside(x)
    assert Heaviside(x, S.Half) == Heaviside(x)
    assert Heaviside(x, y).free_symbols == {x, y}
    assert Heaviside(x, 0).is_real is True
    assert pickle.loads(pickle.dumps(Heaviside(x, 0))) == Heaviside(x, 0)
    import sympy as sp
    native = sx.to_sympy(Heaviside(x, 0))
    assert native == sp.Heaviside(sp.Symbol('x'), 0)
    assert sx.from_sympy(sp.Heaviside(sp.Symbol('x'), 0)) == Heaviside(x, 0)
    assert sx.from_sympy(sp.Heaviside(sp.Symbol('x'))) == Heaviside(x)


def test_heaviside_value_at_zero_does_not_change_transforms():
    # SymPy prints 2*exp(-2*s)/s + exp(-2*s)/s**2; symox keeps the factored
    # form also for Heaviside(t - 2) (a separate, recorded form gap)
    value = laplace_transform(Heaviside(t - 2, 1)*t, t, s, noconds=True)
    assert value == laplace_transform(Heaviside(t - 2)*t, t, s, noconds=True)
    assert sx.expand(value) == 2*exp(-2*s)/s + exp(-2*s)/s**2
    assert str(integrate(Heaviside(x, 0), x)) == 'x*Heaviside(x)'
    assert sx.lambdify(x, Heaviside(x, 0))(0) == 0
    assert sx.lambdify(x, Heaviside(x, 1))(0) == 1


def test_heaviside_printing():
    assert latex(Heaviside(x, 0)) == r'\theta\left(x, 0\right)'
    assert latex(Heaviside(x)) == r'\theta\left(x\right)'
    assert pretty(Heaviside(x, 0), use_unicode=True) == 'θ(x, 0)'
    assert ccode(Heaviside(x, 0)) == '(((x <= 0) ? (\n   0\n)\n: (\n   1\n)))'
    assert ccode(Heaviside(x, y)) == ('(((x < 0) ? (\n   0\n)\n: ((x == 0) ? (\n   y\n)\n'
                                      ': (\n   1\n))))')


def test_heaviside_errors():
    with pytest.raises(AttributeError):
        Heaviside(x, None)
    with pytest.raises(ValueError):
        Heaviside(I, 0)


# --- a number times a noncommutative sum stays a product (Mul.flatten) ---

A, B, C = symbols('A B C', commutative=False)


@pytest.mark.parametrize('build, tree', [
    (lambda: -(A + B), "Mul(Integer(-1), Add(Symbol('A', commutative=False), Symbol('B', commutative=False)))"),
    (lambda: 2*(A + B), "Mul(Integer(2), Add(Symbol('A', commutative=False), Symbol('B', commutative=False)))"),
    (lambda: -(A + 1), "Mul(Integer(-1), Add(Integer(1), Symbol('A', commutative=False)))"),
    (lambda: 2*(x + y), "Add(Mul(Integer(2), Symbol('x')), Mul(Integer(2), Symbol('y')))"),
    (lambda: expand(-(A + B)), "Add(Mul(Integer(-1), Symbol('A', commutative=False)), "
                               "Mul(Integer(-1), Symbol('B', commutative=False)))"),
])
def test_number_times_noncommutative_sum(build, tree):
    assert srepr(build()) == tree


def test_noncommutative_sum_differences():
    assert str(-(A + B) + (A + B)) == 'A + B - (A + B)'
    assert expand(-(A + B) + (A + B)) == 0
    assert (A*(B + C)).equals(A*B + A*C) is True


# --- powers that cannot be materialized ---

def test_held_power_with_huge_integer_part():
    value = sx.Pow(2, Rational(10**20 + 1, 2))
    assert value.is_Pow and value.exp == Rational(10**20 + 1, 2)


def test_infeasible_binomial_expansion_raises_quickly():
    start = time.perf_counter()
    with pytest.raises(NotImplementedError, match='limit'):
        sx.Pow(x + I, 10**9, evaluate=False).as_real_imag()
    assert time.perf_counter() - start < 5
