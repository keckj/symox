"""Wave-4 core wrong-result parity (WP-G1/WP-G2) against SymPy 1.14.

Expected strings were produced by the reference SymPy 1.14 on the same
expressions: infinity tables of Add/Mul/Pow, comp/all_close, check_assumptions,
unchanged, WildFunction matching, bound function variables, Factors/Term and
the satisfiability-based Boolean.equals.
"""
import pytest

import symox
from symox import Function, Integral, Symbol, WildFunction, symbols
from symox.core.numbers import all_close, comp
from symox.core.expr import unchanged
from symox.core.traversal import iterfreeargs
from symox.core.exprtools import Factors, Term
from symox.abc import x, y, z, n

f, g = Function('f'), Function('g')
A, B, C = symbols('A B C')
F = WildFunction('F')

# Evaluation namespace of the probes; symox.test (the suite runner) must not
# become a module-level name collected by pytest.
NAMESPACE = {name: getattr(symox, name) for name in symox.__all__ if name != 'test'} \
    if hasattr(symox, '__all__') else {}
NAMESPACE.update(
    S=symox.S, all_close=all_close, comp=comp, unchanged=unchanged,
    iterfreeargs=iterfreeargs, Factors=Factors, Term=Term,
    x=x, y=y, z=z, n=n, f=f, g=g, A=A, B=B, C=C, F=F,
    p=Symbol('p', positive=True), r=Symbol('r', real=True), q=Symbol('q', negative=True))


def _render(expression):
    try:
        value = eval(expression, dict(NAMESPACE))
    except Exception as error:
        return 'EXC ' + type(error).__name__
    if isinstance(value, list):
        return '[' + ', '.join(map(str, value)) + ']'
    if isinstance(value, set):
        return 'set:' + str(sorted(map(str, value)))
    return str(value)


CORE_CASES = [
    ('(zoo) + (zoo)', 'nan'),
    ('(zoo) ** (zoo)', 'nan'),
    ('(zoo) - (zoo)', 'nan'),
    ('(oo) ** (zoo)', 'nan'),
    ('(zoo) ** (-oo)', '0'),
    ('(-oo) ** (zoo)', 'nan'),
    ('(I*oo) ** (zoo)', 'nan'),
    ('(-I*oo) ** (zoo)', 'nan'),
    ('(0) ** (zoo)', 'nan'),
    ('(-1) ** (zoo)', 'nan'),
    ('(I) ** (zoo)', 'nan'),
    ('(x) ** (zoo)', 'nan'),
    ('(2) ** (zoo)', 'nan'),
    ('(S.Half) ** (zoo)', 'nan'),
    ('(1.5) ** (zoo)', 'nan'),
    ('(oo) ** (oo)', 'oo'),
    ('(oo) ** (-oo)', '0'),
    ('(-oo) ** (oo)', 'nan'),
    ('(oo) ** (I*oo)', 'nan'),
    ('(oo) ** (-I*oo)', 'nan'),
    ('(-1) ** (oo)', 'nan'),
    ('(oo) ** (I)', 'nan'),
    ('(2) ** (oo)', 'oo'),
    ('(S.Half) ** (oo)', '0'),
    ('(1.5) ** (oo)', 'oo'),
    ('(-oo) ** (-oo)', 'nan'),
    ('(-oo) ** (I*oo)', 'nan'),
    ('(-oo) ** (-I*oo)', 'nan'),
    ('(0) ** (-oo)', 'zoo'),
    ('(-1) ** (-oo)', 'nan'),
    ('(-oo) ** (I)', 'nan'),
    ('(2) ** (-oo)', '0'),
    ('(-oo) ** (S.Half)', 'oo*I'),
    ('(S.Half) ** (-oo)', 'oo'),
    ('(-oo) ** (1.5)', '-oo*I'),
    ('(1.5) ** (-oo)', '0'),
    ('(I*oo) + (I*oo)', 'oo*I'),
    ('(I*oo) - (I*oo)', 'nan'),
    ('(I*oo) + (-I*oo)', 'nan'),
    ('(I*oo) - (-I*oo)', 'oo*I'),
    ('(0) ** (I*oo)', 'nan'),
    ('(1) ** (I*oo)', 'nan'),
    ('(I*oo) + (I)', 'oo*I'),
    ('(I*oo) - (I)', 'oo*I'),
    ('(I*oo) ** (S.Half)', 'oo*sqrt(I)'),
    ('(I*oo) ** (1.5)', 'oo*I**1.5'),
    ('comp(pi.n(4), 3.142)', 'True'),
    ('comp(sqrt(2).n(3), 1.41)', 'True'),
    ('comp(S(1)/3, 0.333)', 'True'),
    ('comp(1+I, 1.0001+I)', 'True'),
    ('comp(1e-10, 0)', 'True'),
    ('all_close(x, 1.0*x)', 'True'),
    ('all_close(x, x + 1e-10)', 'True'),
    ('all_close(x+y, 1.0*x+y)', 'True'),
    ('all_close(exp(x), exp(1.0*x))', 'True'),
    ('all_close(x*y, 1.0*x*y)', 'True'),
    ('all_close(pi, 3.14159265)', 'False'),
    ('all_close(x, y)', 'True'),
    ('all_close(sqrt(2), 1.4142135)', 'False'),
    ('check_assumptions(Symbol("x", positive=True), 3)', 'None'),
    ('check_assumptions(Symbol("p", positive=True), Symbol("q", positive=True))', 'True'),
    ('check_assumptions(S.Half, integer=None)', 'True'),
    ('unchanged(Piecewise, (x, x > 1), (0, True))', 'True'),
    ('unchanged(Integral, x, (x, 0, 1))', 'True'),
    ('F.match(F)', '{F_: F_}'),
    ('F.matches(F)', '{F_: F_}'),
    ('list(iterfreeargs(Integral(f(x), (f(x), 1))))', '[Integral(f(x), (f(x), 1)), 1]'),
    ('Integral(f(x), (f(x), 1, g(y))).free_symbols', "set:['y']"),
    ('Integral(f(x), (f(x), 1, 2)).free_symbols', 'set:[]'),
    ('Factors(x*y**2).gcd(Factors(x*y/z))', 'Factors({x: 1, y: 1})'),
    ('Factors(x*y**2).lcm(Factors(x*y/z))', 'Factors({x: 1, y: 2, z: -1})'),
    ('Factors(x*y**2).quo(Factors(x*y/z))', 'Factors({y: 1})'),
    ('Factors(x*y**2).rem(Factors(x*y/z))', 'Factors({z: -1})'),
    ('Factors(x*y**2).div(Factors(x*y/z))', '(Factors({y: 1}), Factors({z: -1}))'),
    ('Factors(x*y**2).normal(Factors(x*y/z))', '(Factors({y: 1}), Factors({z: -1}))'),
    ('Factors((x*y**2).as_powers_dict()).gcd(Factors((x*y/z).as_powers_dict()))', 'Factors({x: 1, y: 1})'),
    ('Factors(x**y).gcd(Factors(x**(y+1)))', 'Factors({x: y})'),
    ('Factors(x**(y+2)).div(Factors(x**(y+1)))', '(Factors({x: 1}), Factors({}))'),
    ('Factors(x**(y+2)).normal(Factors(x**3))', '(Factors({x: y}), Factors({x: 1}))'),
    ('Factors(-x*I)', 'Factors({-1: 1, I: 1, x: 1})'),
    ('Factors(x**2*y).normal(x*y**2)', '(Factors({x: 1}), Factors({y: 1}))'),
    ('Factors(x**2).div(x**3)', '(Factors({}), Factors({x: 1}))'),
    ('Factors(x**2/y).gcd(Factors(y/x))', 'Factors({x: -1, y: -1})'),
    ('Term(4*x*y/z).lcm(Term(2*x/z**2))', 'Term(4, Factors({x: 1, y: 1}), Factors({z: 2}))'),
    ('(A >> B).equals(~B >> ~A)', 'True'),
    ('(A | B).equals(~(~A & ~B))', 'True'),
    ('Xor(A, B).equals(Or(And(A, ~B), And(~A, B)))', 'True'),
    ('(A & (B | C)).equals((A & B) | (A & C))', 'True'),
    ('Equivalent(A, B).equals((A >> B) & (B >> A))', 'True'),
    ('true.equals(false)', 'False'),
    ('(A & (x > 1)).equals(A & (x > 1))', 'EXC NotImplementedError'),
    ('ITE(A, B, C).equals((A & B) | (~A & C))', 'True'),
    ('Nand(A,B).equals(~A | ~B)', 'True'),
    ('comp(1, 1.0)', 'True'),
    ('Term(4*x*y/z).gcd(Term(2*x/z**2))', 'Term(2, Factors({x: 1}), Factors({z: 1}))'),
]

TABLE_CASES = [
    ('Pow(0, -oo)', 'zoo'),
    ('Pow(0, zoo)', 'nan'),
    ('Pow(0, I*oo)', 'nan'),
    ('Pow(0, -I*oo)', 'nan'),
    ('Pow(0, I)', 'nan'),
    ('Pow(0, p)', '0'),
    ('Pow(0, q)', 'zoo'),
    ('Pow(0, 1+I)', 'nan'),
    ('Pow(0, -1+I)', 'nan'),
    ('Pow(0, 0.5*I)', 'nan'),
    ('Pow(1, I*oo)', 'nan'),
    ('Pow(1, -I*oo)', 'nan'),
    ('Pow(-1, oo)', 'nan'),
    ('Pow(-1, -oo)', 'nan'),
    ('Pow(-1, zoo)', 'nan'),
    ('Pow(2, oo)', 'oo'),
    ('Pow(2, -oo)', '0'),
    ('Pow(2, zoo)', 'nan'),
    ('Pow(-2, oo)', 'zoo'),
    ('Pow(-2, -oo)', '0'),
    ('Pow(-2, zoo)', 'nan'),
    ('Pow(S.Half, oo)', '0'),
    ('Pow(S.Half, -oo)', 'oo'),
    ('Pow(S.Half, zoo)', 'nan'),
    ('Pow(-S.Half, oo)', '0'),
    ('Pow(-S.Half, -oo)', 'nan'),
    ('Pow(-S.Half, zoo)', 'nan'),
    ('Pow(1.5, oo)', 'oo'),
    ('Pow(1.5, -oo)', '0'),
    ('Pow(1.5, zoo)', 'nan'),
    ('Pow(-1.5, oo)', 'zoo'),
    ('Pow(-1.5, -oo)', 'nan'),
    ('Pow(-1.5, zoo)', 'nan'),
    ('Pow(I, zoo)', 'nan'),
    ('Pow(1+I, zoo)', 'nan'),
    ('Pow(oo, oo)', 'oo'),
    ('Pow(oo, -oo)', '0'),
    ('Pow(oo, zoo)', 'nan'),
    ('Pow(oo, I*oo)', 'nan'),
    ('Pow(oo, -I*oo)', 'nan'),
    ('Pow(oo, I)', 'nan'),
    ('Pow(oo, p)', 'oo'),
    ('Pow(oo, q)', '0'),
    ('Pow(oo, 1+I)', 'zoo'),
    ('Pow(oo, -1+I)', '0'),
    ('Pow(oo, 0.5*I)', 'nan'),
    ('Pow(-oo, oo)', 'nan'),
    ('Pow(-oo, -oo)', 'nan'),
    ('Pow(-oo, zoo)', 'nan'),
    ('Pow(-oo, I*oo)', 'nan'),
    ('Pow(-oo, -I*oo)', 'nan'),
    ('Pow(-oo, I)', 'nan'),
    ('Pow(-oo, S.Half)', 'oo*I'),
    ('Pow(-oo, 1.5)', '-oo*I'),
    ('Pow(-oo, -S.Half)', '0'),
    ('Pow(-oo, 1+I)', 'zoo'),
    ('Pow(-oo, -1+I)', '0'),
    ('Pow(-oo, 0.5*I)', 'nan'),
    ('Pow(zoo, -oo)', '0'),
    ('Pow(zoo, zoo)', 'nan'),
    ('Pow(I*oo, zoo)', 'nan'),
    ('Pow(I*oo, S.Half)', 'oo*sqrt(I)'),
    ('Pow(I*oo, 1.5)', 'oo*I**1.5'),
    ('Pow(I*oo, -S.Half)', '0'),
    ('Pow(-I*oo, zoo)', 'nan'),
    ('Pow(-I*oo, S.Half)', 'oo*sqrt(-I)'),
    ('Pow(-I*oo, 1.5)', 'oo*(-I)**1.5'),
    ('Pow(-I*oo, -S.Half)', '0'),
    ('Pow(x, zoo)', 'nan'),
    ('Pow(p, zoo)', 'nan'),
    ('Pow(p+1, oo)', 'oo'),
    ('Pow(p+1, zoo)', 'nan'),
    ('Pow(q-1, oo)', 'zoo'),
    ('Pow(q-1, zoo)', 'nan'),
    ('Pow(2*I, zoo)', 'nan'),
    ('Pow(pi, oo)', 'oo'),
    ('Pow(pi, zoo)', 'nan'),
    ('Pow(-pi, oo)', 'zoo'),
    ('Pow(-pi, zoo)', 'nan'),
    ('Pow(S.Half*I, zoo)', 'nan'),
    ('Add(zoo, zoo)', 'nan'),
    ('Add(zoo, zoo + x)', 'nan'),
    ('Add(I*oo, I*oo)', 'oo*I'),
    ('Add(I*oo, -I*oo)', 'nan'),
    ('Add(I*oo, I)', 'oo*I'),
    ('Add(I*oo, 2*I*oo)', 'oo*I'),
    ('Add(I*oo, oo*I + 1)', '1 + oo*I'),
    ('Add(-I*oo, -I*oo)', '-oo*I'),
    ('Add(-I*oo, I)', '-oo*I'),
    ('Add(-I*oo, 2*I*oo)', 'nan'),
    ('Add(-I*oo, oo*I + 1)', 'nan'),
    ('Add(x, x*oo)', 'oo*x'),
    ('Add(I, 2*I*oo)', 'oo*I'),
    ('Add(I, oo*I + 1)', '1 + oo*I'),
    ('Add(2*I*oo, 2*I*oo)', 'oo*I'),
    ('Add(2*I*oo, oo*I + 1)', '1 + oo*I'),
    ('Add((1+I)*oo, (1+I)*oo)', 'oo*(1 + I)'),
    ('Add(x*oo, x*oo)', 'oo*x'),
    ('Add(oo*I + 1, oo*I + 1)', '2 + oo*I'),
    ('Add(zoo + x, zoo + x)', 'nan'),
]


@pytest.mark.parametrize('expression, expected', CORE_CASES)
def test_core_reference_values(expression, expected):
    assert _render(expression) == expected


@pytest.mark.parametrize('expression, expected', TABLE_CASES)
def test_infinity_tables(expression, expected):
    assert _render(expression) == expected


def test_function_integration_variable_is_bound():
    expr = Integral(f(x)*y, (f(x), 1, g(z)))
    assert expr.free_symbols == {y, z}
    assert Integral(f(x), (f(x), 1, 2)).free_symbols == set()


def test_boolean_equals_with_relationals_is_not_decided():
    with pytest.raises(NotImplementedError):
        (A & (x > 1)).equals(A & (x > 1))
