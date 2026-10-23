"""Printing parity, audit batch print-a3: sets, Order, hyper, meijerg,
Integral, algebraic and expression domain elements.

The expected strings are hard-coded from SymPy 1.14.
"""
import pytest

from symox import (EX, EXRAW, QQ, Complement, ConditionSet, Eq, FiniteSet, ImageSet,
                   Integral, Intersection, Interval, Lambda, O, Range, Rational, S,
                   Union, construct_domain, hyper, latex, meijerg, oo, pi, pretty,
                   sin, solveset, sqrt, stationary_points)
from symox.abc import a, n, x, y

CASES = [
    ('FiniteSet(x/2, 3)',
     '    x \n{3, -}\n    2 ',
     '⎧   x⎫\n⎨3, ─⎬\n⎩   2⎭',
     '\\left\\{3, \\frac{x}{2}\\right\\}'),
    ('Interval.Lopen(-1, x)',
     '(-1, x]',
     '(-1, x]',
     '\\left(-1, x\\right]'),
    ('Union(Interval(0, 1), FiniteSet(3))',
     '[0, 1] U {3}',
     '[0, 1] ∪ {3}',
     '\\left[0, 1\\right] \\cup \\left\\{3\\right\\}'),
    ('Complement(Interval(0, x), FiniteSet(y))',
     '[0, x] \\ {y}',
     '[0, x] \\ {y}',
     '\\left[0, x\\right] \\setminus \\left\\{y\\right\\}'),
    ('Intersection(Interval(0, x), Interval(y, 5))',
     '[0, x] n [y, 5]',
     '[0, x] ∩ [y, 5]',
     '\\left[0, x\\right] \\cap \\left[y, 5\\right]'),
    ('S.Integers',
     'Integers',
     'ℤ',
     '\\mathbb{Z}'),
    ('S.Naturals0',
     'Naturals0',
     'ℕ₀',
     '\\mathbb{N}_0'),
    ('S.Reals',
     '(-oo, oo)',
     'ℝ',
     '\\mathbb{R}'),
    ('S.EmptySet',
     'EmptySet',
     '∅',
     '\\emptyset'),
    ('ImageSet(Lambda(n, n*pi/2), S.Integers)',
     ' pi*n                 \n{---- | n in Integers}\n  2                   ',
     '⎧π⋅n │      ⎫\n⎨─── │ n ∊ ℤ⎬\n⎩ 2  │      ⎭',
     '\\left\\{\\frac{\\pi n}{2}\\; \\middle|\\; n \\in \\mathbb{Z}\\right\\}'),
    ('ImageSet(Lambda((x, y), x + y), S.Integers, S.Naturals)',
     '{x + y | x in Integers, y in Naturals}',
     '{x + y │ x ∊ ℤ, y ∊ ℕ}',
     '\\left\\{x + y\\; \\middle|\\; x \\in \\mathbb{Z}, y \\in \\mathbb{N}\\right\\}'),
    ('ConditionSet(x, Eq(sin(x), x), S.Reals)',
     '{x | x in (-oo, oo) and sin(x) = x}',
     '{x │ x ∊ ℝ ∧ (sin(x) = x)}',
     '\\left\\{x\\; \\middle|\\; x \\in \\mathbb{R} \\wedge \\sin{\\left(x \\right)} = x \\right\\}'),
    ('ConditionSet(x, Eq(x**2 + sin(x), 0))',
     '      2              \n{x | x  + sin(x) = 0}',
     '⎧  │ ⎛ 2             ⎞⎫\n⎨x │ ⎝x  + sin(x) = 0⎠⎬\n⎩  │                  ⎭',
     '\\left\\{x\\; \\middle|\\; x^{2} + \\sin{\\left(x \\right)} = 0 \\right\\}'),
    ('solveset(sin(x), x)',
     '{2*n*pi | n in Integers} U {2*n*pi + pi | n in Integers}',
     '{2⋅n⋅π │ n ∊ ℤ} ∪ {2⋅n⋅π + π │ n ∊ ℤ}',
     '\\left\\{2 n \\pi\\; \\middle|\\; n \\in \\mathbb{Z}\\right\\} \\cup \\left\\{2 n \\pi + \\pi\\; \\middle|\\; n \\in \\mathbb{Z}\\right\\}'),
    ('Range(-4, 5)',
     '{-4, -3, ..., 4}',
     '{-4, -3, …, 4}',
     '\\left\\{-4, -3, \\ldots, 4\\right\\}'),
    ('Range(10, 0, -3)',
     '{10, 7, 4, 1}',
     '{10, 7, 4, 1}',
     '\\left\\{10, 7, \\ldots, 1\\right\\}'),
    ('Range(0, oo, 2)',
     '{0, 2, ...}',
     '{0, 2, …}',
     '\\left\\{0, 2, \\ldots\\right\\}'),
    ('stationary_points(x**3 - 3*x, x, S.Reals)',
     '{-1, 1}',
     '{-1, 1}',
     '\\left\\{-1, 1\\right\\}'),
    ('O(x**6)',
     ' / 6\\\nO\\x /',
     ' ⎛ 6⎞\nO⎝x ⎠',
     'O\\left(x^{6}\\right)'),
    ('O(1/x, (x, oo))',
     ' /1         \\\nO|-; x -> oo|\n \\x         /',
     ' ⎛1       ⎞\nO⎜─; x → ∞⎟\n ⎝x       ⎠',
     'O\\left(\\frac{1}{x}; x\\rightarrow \\infty\\right)'),
    ('O(x*y)',
     'O(x*y; (x, y) -> (0, 0))',
     'O(x⋅y; (x, y) → (0, 0))',
     'O\\left(x y; \\left( x, \\  y\\right)\\rightarrow \\left( 0, \\  0\\right)\\right)'),
    ('1 + x + O(x**2)',
     '         / 2\\\n1 + x + O\\x /',
     '         ⎛ 2⎞\n1 + x + O⎝x ⎠',
     '1 + x + O\\left(x^{2}\\right)'),
    ('hyper([1, 2], [4], x)',
     '  _            \n |_  /1, 2 |  \\\n |   |     | x|\n2  1 \\ 4   |  /',
     ' ┌─  ⎛1, 2 │  ⎞\n ├─  ⎜     │ x⎟\n2╵ 1 ⎝ 4   │  ⎠',
     '{{}_{2}F_{1}\\left(\\begin{matrix} 1, 2 \\\\ 4 \\end{matrix}\\middle| {x} \\right)}'),
    ('hyper((a,), (), x**2)',
     '  _          \n |_  /a |  2\\\n |   |  | x |\n1  0 \\  |   /',
     ' ┌─  ⎛a │  2⎞\n ├─  ⎜  │ x ⎟\n1╵ 0 ⎝  │   ⎠',
     '{{}_{1}F_{0}\\left(\\begin{matrix} a \\\\  \\end{matrix}\\middle| {x^{2}} \\right)}'),
    ('meijerg([[1, 2], [3]], [[4], []], x)',
     ' __1, 2 /1, 2  3 |  \\\n/__     |        | x|\n\\_|3, 1 \\ 4      |  /',
     '╭─╮1, 2 ⎛1, 2  3 │  ⎞\n│╶┐     ⎜        │ x⎟\n╰─╯3, 1 ⎝ 4      │  ⎠',
     '{G_{3, 1}^{1, 2}\\left(\\begin{matrix} 1, 2 & 3 \\\\4 &  \\end{matrix} \\middle| {x} \\right)}'),
    ('meijerg((1, 2), (a, 4), (5,), [], x)',
     ' __1, 2 /1, 2  4, a |  \\\n/__     |           | x|\n\\_|4, 1 \\ 5         |  /',
     '╭─╮1, 2 ⎛1, 2  4, a │  ⎞\n│╶┐     ⎜           │ x⎟\n╰─╯4, 1 ⎝ 5         │  ⎠',
     '{G_{4, 1}^{1, 2}\\left(\\begin{matrix} 1, 2 & 4, a \\\\5 &  \\end{matrix} \\middle| {x} \\right)}'),
    ('Integral(sin(x), x)',
     '  /         \n |          \n | sin(x) dx\n |          \n/           ',
     '⌠          \n⎮ sin(x) dx\n⌡          ',
     '\\int \\sin{\\left(x \\right)}\\, dx'),
    ('Integral(x**2, (x, 0, 1))',
     '  1      \n  /      \n |       \n |   2   \n |  x  dx\n |       \n/        \n0        ',
     '1      \n⌠      \n⎮  2   \n⎮ x  dx\n⌡      \n0      ',
     '\\int\\limits_{0}^{1} x^{2}\\, dx'),
    ('Integral(x*y, x, y)',
     '  /  /         \n |  |          \n |  | x*y dx dy\n |  |          \n/  /           ',
     '⌠ ⌠          \n⎮ ⎮ x⋅y dx dy\n⌡ ⌡          ',
     '\\iint x y\\, dx\\, dy'),
    ('Integral(x + 1, x)',
     '  /          \n |           \n | (x + 1) dx\n |           \n/            ',
     '⌠           \n⎮ (x + 1) dx\n⌡           ',
     '\\int \\left(x + 1\\right)\\, dx'),
    ('2*Integral(x, (x, 0, 1)) + y',
     '        1     \n        /     \n       |      \ny + 2* |  x dx\n       |      \n      /       \n      0       ',
     '      1     \n      ⌠     \ny + 2⋅⎮ x dx\n      ⌡     \n      0     ',
     'y + 2 \\int\\limits_{0}^{1} x\\, dx'),
    ('Integral(1/x, (x, 1, oo))',
     ' oo     \n  /     \n |      \n |  1   \n |  - dx\n |  x   \n |      \n/       \n1       ',
     '∞     \n⌠     \n⎮ 1   \n⎮ ─ dx\n⎮ x   \n⌡     \n1     ',
     '\\int\\limits_{1}^{\\infty} \\frac{1}{x}\\, dx'),
]


@pytest.mark.parametrize('source, ascii_form, unicode_form, tex', CASES)
def test_pretty_and_latex(source, ascii_form, unicode_form, tex):
    value = eval(source)
    assert pretty(value, use_unicode=False) == ascii_form
    assert pretty(value, use_unicode=True) == unicode_form
    assert latex(value) == tex


def test_integral_latex_modes():
    value = Integral(x**2, (x, 0, 1))
    assert latex(value, mode='inline') == '$\\int_{0}^{1} x^{2}\\, dx$'
    assert latex(value, itex=True) == '\\int_{0}^{1} x^{2}\\, dx'
    assert latex(value) == '\\int\\limits_{0}^{1} x^{2}\\, dx'


def test_algebraic_and_expression_elements():
    K = QQ.algebraic_field(sqrt(2))
    assert repr(K.convert(1 + sqrt(2))) == 'ANP([MPQ(1,1), MPQ(1,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)'
    assert str(K.convert(sqrt(2)/3)) == 'ANP([MPQ(1,3), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)'
    assert repr(EX.convert(sin(x))) == 'EX(sin(x))'
    assert str(EX.convert(x + 1)) == 'EX(x + 1)'
    assert repr(EXRAW.convert(x/2)) == 'x/2'


def test_construct_domain_extension_generators():
    assert str(construct_domain([1 + sqrt(2), 1], extension=True)[0]) == 'QQ<sqrt(2)>'
    assert str(construct_domain([sqrt(2)*sqrt(3), 1], extension=True)[0]) == 'QQ<sqrt(6)>'
