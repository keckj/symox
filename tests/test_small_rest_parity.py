"""Small remaining SymPy 1.14 parity gaps: printing order, container
pretty printing, division options and orthogonal polynomial errors.

Every expected value below was produced by the reference SymPy 1.14."""

import pytest

from symox import (Derivative, Function, I, Integral, Limit, OptionError,
                   Poly, Product, QQ, Rational, Subs, Sum, Tuple, assoc_laguerre,
                   div, im, jacobi, laguerre, latex, pdiv, pi, prem,
                   pretty, quo, re, sqrt, symbols)
from symox.physics.quantum.constants import hbar

x, y, n, a, b = symbols('x y n a b')
A, B = symbols('A B', commutative=False)
f = Function('f')


@pytest.mark.parametrize('expr, text, tex', [
    (lambda: Sum(x, (x, 1, n)) + Limit(1/x, x, 0),
     "Sum(x, (x, 1, n)) + Limit(1/x, x, 0, dir='+')",
     r'\sum_{x=1}^{n} x + \lim_{x \to 0^+} \frac{1}{x}'),
    (lambda: Sum(x, (x, 1, n)) + Integral(x, (x, 0, 1)),
     'Sum(x, (x, 1, n)) + Integral(x, (x, 0, 1))',
     r'\sum_{x=1}^{n} x + \int\limits_{0}^{1} x\, dx'),
    (lambda: Product(x, (x, 1, n)) + Limit(1/x, x, 0),
     "Product(x, (x, 1, n)) + Limit(1/x, x, 0, dir='+')",
     r'\prod_{x=1}^{n} x + \lim_{x \to 0^+} \frac{1}{x}'),
    (lambda: Subs(f(x), x, y) + Integral(x, (x, 0, 1)),
     'Subs(f(x), x, y) + Integral(x, (x, 0, 1))',
     r'\left. f{\left(x \right)} \right|_{\substack{ x=y }} + \int\limits_{0}^{1} x\, dx'),
    (lambda: Derivative(f(x), x) + Sum(x, (x, 1, n)),
     'Derivative(f(x), x) + Sum(x, (x, 1, n))',
     r'\frac{d}{d x} f{\left(x \right)} + \sum_{x=1}^{n} x'),
    (lambda: Subs(f(x), x, y) + Sum(x, (x, 1, n)),
     'Subs(f(x), x, y) + Sum(x, (x, 1, n))',
     r'\left. f{\left(x \right)} \right|_{\substack{ x=y }} + \sum_{x=1}^{n} x'),
    (lambda: x + re(A) + im(B) + 1, 'x + re(A) + im(B) + 1',
     r'x + \operatorname{re}{\left(A\right)} + \operatorname{im}{\left(B\right)} + 1'),
    (lambda: re(A) + A, 're(A) + A', r'\operatorname{re}{\left(A\right)} + A'),
    (lambda: re(A) + im(B) + im(x), 're(A) + im(B) + im(x)',
     r'\operatorname{re}{\left(A\right)} + \operatorname{im}{\left(B\right)} + \operatorname{im}{\left(x\right)}'),
    (lambda: x + hbar, 'x + hbar', r'x + \hbar'),
    (lambda: x**2 + y + 2*hbar, 'x**2 + y + 2*hbar', r'x^{2} + y + 2 \hbar'),
    (lambda: Sum(x, (x, 1, y)) + hbar, 'Sum(x, (x, 1, y)) + hbar', r'\sum_{x=1}^{y} x + \hbar'),
])
def test_add_order_str_latex(expr, text, tex):
    e = expr()
    assert str(e) == text
    assert latex(e) == tex


def test_re_im_noncommutative_and_hbar_number():
    assert re(A).is_commutative is True
    assert im(B).is_commutative is True
    assert hbar.is_number is True
    assert pretty(x + hbar) == 'x + ℏ'


@pytest.mark.parametrize('obj, uni, asc', [
    (lambda: [x**2, x/2], '⎡ 2  x⎤\n⎢x , ─⎥\n⎣    2⎦', '  2  x \n[x , -]\n     2 '),
    (lambda: (x/2,), '⎛x ⎞\n⎜─,⎟\n⎝2 ⎠', ' x  \n(-,)\n 2  '),
    (lambda: {x/2: y, 1: x**2}, '⎧    2  x   ⎫\n⎨1: x , ─: y⎬\n⎩       2   ⎭',
     '     2  x    \n{1: x , -: y}\n        2    '),
    (lambda: {x/2, y}, '⎧x   ⎫\n⎨─, y⎬\n⎩2   ⎭', ' x    \n{-, y}\n 2    '),
    (lambda: frozenset([x/3, y]), '         ⎛⎧x   ⎫⎞\nfrozenset⎜⎨─, y⎬⎟\n         ⎝⎩3   ⎭⎠',
     '           x     \nfrozenset({-, y})\n           3     '),
    (lambda: set(), 'set()', 'set()'),
    (lambda: Tuple(x/2, 1), '⎛x   ⎞\n⎜─, 1⎟\n⎝2   ⎠', ' x    \n(-, 1)\n 2    '),
    (lambda: [x/2, [x, (1, x/3)]], '⎡x  ⎡   ⎛   x⎞⎤⎤\n⎢─, ⎢x, ⎜1, ─⎟⎥⎥\n⎣2  ⎣   ⎝   3⎠⎦⎦',
     ' x          x   \n[-, [x, (1, -)]]\n 2          3   '),
])
def test_container_pretty(obj, uni, asc):
    o = obj()
    assert pretty(o) == uni
    assert pretty(o, use_unicode=False) == asc


def test_division_greedy_and_symmetric_options():
    assert div(x**3 + 2*x + 1, x - 2, greedy=True) == (x**2 + 2*x + 6, 13)
    assert div(x**2 + 1, x + 1, greedy=None, symmetric=None) == (x - 1, 2)
    assert div(x**2 - 1, x - 1, auto=False) == (x + 1, 0)
    q, r = div(x**2 + 1, x + 2, modulus=5, symmetric=True, polys=True)
    assert str(q) == 'Poly(x - 2, x, modulus=5)' and str(r) == 'Poly(0, x, modulus=5)'
    q, r = div(sqrt(2)*x**2 + 1, sqrt(3)*x + pi, greedy=False)
    assert str(q) == 'sqrt(6)*x/3 - sqrt(2)*pi/3' and str(r) == '1 + sqrt(2)*pi**2/3'
    assert str(quo(x**2 + sqrt(2)*x, sqrt(2)*x, greedy=False)) == 'sqrt(2)*x/2 + 1'
    q, r = div(x**2 + sqrt(2)*x, sqrt(2)*x, x, greedy=False, polys=True)
    assert str(q) == "Poly(sqrt(2)/2*x + 1, x, domain='EX')"
    assert str(r) == "Poly(0, x, domain='EX')"
    assert str(pdiv(x**2 + sqrt(2)*x, sqrt(2)*x + 1, greedy=False)) == '(sqrt(2)*x + 1, -1)'
    assert prem(x**2 + sqrt(2)*x, sqrt(2)*x + 1, greedy=False) == -1


@pytest.mark.parametrize('kwargs, message', [
    ({'symmetric': True}, "'symmetric' option is only allowed together with 'modulus'"),
    ({'domain': QQ, 'greedy': True}, "'domain' option is not allowed together with 'greedy'"),
    ({'greedy': True, 'modulus': 3}, "'greedy' option is not allowed together with 'modulus'"),
    ({'modulus': 3, 'greedy': True}, "'modulus' option is not allowed together with 'greedy'"),
    ({'field': True, 'domain': QQ}, "'field' option is not allowed together with 'domain'"),
    ({'extension': True, 'greedy': False}, "'extension' option is not allowed together with 'greedy'"),
    ({'greedy': 'a'}, "'greedy' must have a boolean value assigned, got a"),
])
def test_division_option_errors(kwargs, message):
    with pytest.raises(OptionError) as info:
        div(x**2 + 1, x, **kwargs)
    assert str(info.value) == message


@pytest.mark.parametrize('call, message', [
    (lambda: jacobi(Rational(-1, 2), a, b, x), 'Cannot generate Jacobi polynomial of index -1/2'),
    (lambda: jacobi(Rational(-3, 2), 1, 2, Rational(1, 3)), 'Cannot generate Jacobi polynomial of index -3/2'),
    (lambda: jacobi(-1.5, 1, 2, 0.3), 'Cannot generate Jacobi polynomial of index -1.50000000000000'),
    (lambda: assoc_laguerre(Rational(-1, 2), 1, x), 'The index n must be nonnegative integer (got -1/2)'),
    (lambda: assoc_laguerre(-1.5, a, x), 'The index n must be nonnegative integer (got -1.50000000000000)'),
    (lambda: assoc_laguerre(Rational(-1, 2), 0, x), 'Error: n should be an integer.'),
    (lambda: laguerre(Rational(1, 2), x), 'Error: n should be an integer.'),
    (lambda: laguerre(I, x), 'Error: n should be an integer.'),
    (lambda: laguerre(-1 + I, x), 'Error: n should be an integer.'),
    (lambda: jacobi(n, a, b, x).subs(n, Rational(-1, 2)), 'Cannot generate Jacobi polynomial of index -1/2'),
    (lambda: assoc_laguerre(n, a, x).subs(n, Rational(-1, 2)), 'The index n must be nonnegative integer (got -1/2)'),
])
def test_orthopoly_negative_noninteger_degree_errors(call, message):
    with pytest.raises(ValueError) as info:
        call()
    assert str(info.value) == message


def test_orthopoly_negative_degree_accepted_forms():
    # The symmetric Jacobi branches are taken before the degree check, so
    # no error is raised. SymPy then sets gegenbauer(-1/2, ...) to 0, which
    # holds only for negative integer degrees; symox keeps it unevaluated.
    assert jacobi(Rational(-1, 2), a, a, x).has(Rational(-1, 2))
    assert jacobi(-1, a, a, x) == 0
    assert str(assoc_laguerre(-2.0, 0, x)) == 'exp(x)*laguerre(1.0, -x)'
    assert str(assoc_laguerre(-2, 0, x)) == '(x + 1)*exp(x)'
    held = jacobi(Rational(-1, 2), 1, 2, Rational(1, 3), evaluate=False)
    assert str(held.evalf()) == '0.734036281085578'
