"""Noncommutative adjoint/transpose, radsimp, nonlinsolve, assumption
errors and srepr, checked against SymPy 1.14 printed and srepr forms."""
import pickle

import pytest

import symox as sp
from symox.core.facts import InconsistentAssumptions

A, B, C = sp.symbols('A B C', commutative=False)
x, y, z = sp.symbols('x y z')
f = sp.Function('f')
g = sp.Function('g', commutative=False)
H = sp.Symbol('H', commutative=False, hermitian=True)
K = sp.Symbol('K', commutative=False, antihermitian=True)
NS = dict(sp=sp, A=A, B=B, C=C, x=x, y=y, z=z, f=f, g=g, H=H, K=K)

# (expression, str in SymPy 1.14, srepr in SymPy 1.14 or "" when not compared)
REFERENCE = [
    ('sp.adjoint(A*B)', 'adjoint(B)*adjoint(A)', "Mul(adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('sp.transpose(A*B)', 'transpose(B)*transpose(A)', "Mul(transpose(Symbol('B', commutative=False)), transpose(Symbol('A', commutative=False)))"),
    ('sp.adjoint(A)', 'adjoint(A)', "adjoint(Symbol('A', commutative=False))"),
    ('sp.transpose(A)', 'transpose(A)', "transpose(Symbol('A', commutative=False))"),
    ('sp.adjoint(A*B*C)', 'adjoint(C)*adjoint(B)*adjoint(A)', "Mul(adjoint(Symbol('C', commutative=False)), adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('sp.transpose(A*B*C)', 'transpose(C)*transpose(B)*transpose(A)', "Mul(transpose(Symbol('C', commutative=False)), transpose(Symbol('B', commutative=False)), transpose(Symbol('A', commutative=False)))"),
    ('sp.adjoint(2*A*B)', '2*adjoint(B)*adjoint(A)', "Mul(Integer(2), adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('sp.transpose(x*A*B)', 'x*transpose(B)*transpose(A)', "Mul(Symbol('x'), transpose(Symbol('B', commutative=False)), transpose(Symbol('A', commutative=False)))"),
    ('sp.adjoint(x*A*B)', 'conjugate(x)*adjoint(B)*adjoint(A)', "Mul(conjugate(Symbol('x')), adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('sp.adjoint(A + B)', 'adjoint(A) + adjoint(B)', "Add(adjoint(Symbol('A', commutative=False)), adjoint(Symbol('B', commutative=False)))"),
    ('sp.transpose(A + B)', 'transpose(A) + transpose(B)', "Add(transpose(Symbol('A', commutative=False)), transpose(Symbol('B', commutative=False)))"),
    ('sp.adjoint(A*B + C)', 'adjoint(B)*adjoint(A) + adjoint(C)', "Add(Mul(adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False))), adjoint(Symbol('C', commutative=False)))"),
    ('sp.adjoint(A**2)', 'adjoint(A)**2', "Pow(adjoint(Symbol('A', commutative=False)), Integer(2))"),
    ('sp.transpose(A**2)', 'transpose(A)**2', "Pow(transpose(Symbol('A', commutative=False)), Integer(2))"),
    ('sp.adjoint(A**-1)', 'adjoint(A)**(-1)', "Pow(adjoint(Symbol('A', commutative=False)), Integer(-1))"),
    ('sp.transpose(A**x)', 'transpose(A**x)', "transpose(Pow(Symbol('A', commutative=False), Symbol('x')))"),
    ('sp.adjoint(sp.I*A*B)', '-I*adjoint(B)*adjoint(A)', "Mul(Integer(-1), I, adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('sp.adjoint(sp.adjoint(A*B))', 'A*B', "Mul(Symbol('A', commutative=False), Symbol('B', commutative=False))"),
    ('sp.transpose(sp.transpose(A*B))', 'A*B', "Mul(Symbol('A', commutative=False), Symbol('B', commutative=False))"),
    ('sp.adjoint(sp.transpose(A*B))', 'conjugate(A)*conjugate(B)', "Mul(conjugate(Symbol('A', commutative=False)), conjugate(Symbol('B', commutative=False)))"),
    ('sp.conjugate(sp.adjoint(A*B))', 'transpose(B)*transpose(A)', "Mul(transpose(Symbol('B', commutative=False)), transpose(Symbol('A', commutative=False)))"),
    ('(A*B).adjoint()', 'adjoint(B)*adjoint(A)', "Mul(adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('(A*B).transpose()', 'transpose(B)*transpose(A)', "Mul(transpose(Symbol('B', commutative=False)), transpose(Symbol('A', commutative=False)))"),
    ('sp.adjoint(f(A))', 'adjoint(f(A))', "adjoint(Function('f')(Symbol('A', commutative=False)))"),
    ('sp.adjoint(A*B).is_commutative', 'False', ''),
    ('sp.adjoint(A).is_commutative', 'False', ''),
    ('sp.adjoint(A*B).args', '(adjoint(B), adjoint(A))', ''),
    ('sp.adjoint(A*B) - sp.adjoint(B)*sp.adjoint(A)', '0', 'Integer(0)'),
    ('sp.adjoint(A*B).subs(A, x)', 'conjugate(x)*adjoint(B)', "Mul(conjugate(Symbol('x')), adjoint(Symbol('B', commutative=False)))"),
    ('sp.transpose(A*B).subs({A: x, B: y})', 'x*y', "Mul(Symbol('x'), Symbol('y'))"),
    ('sp.adjoint(sp.exp(A))', 'exp(adjoint(A))', "exp(adjoint(Symbol('A', commutative=False)))"),
    ('sp.transpose(sp.sin(A))', 'transpose(sin(A))', "transpose(sin(Symbol('A', commutative=False)))"),
    ('sp.radsimp(A/(1 + sp.sqrt(2)))', '-A + sqrt(2)*A', "Add(Mul(Integer(-1), Symbol('A', commutative=False)), Mul(Pow(Integer(2), Rational(1, 2)), Symbol('A', commutative=False)))"),
    ('sp.radsimp(A*B/(1 + sp.sqrt(2)))', '-A*B + sqrt(2)*A*B', "Add(Mul(Integer(-1), Symbol('A', commutative=False), Symbol('B', commutative=False)), Mul(Pow(Integer(2), Rational(1, 2)), Symbol('A', commutative=False), Symbol('B', commutative=False)))"),
    ('sp.radsimp(A/(sp.sqrt(2) + sp.sqrt(3)))', '-sqrt(2)*A + sqrt(3)*A', "Add(Mul(Integer(-1), Pow(Integer(2), Rational(1, 2)), Symbol('A', commutative=False)), Mul(Pow(Integer(3), Rational(1, 2)), Symbol('A', commutative=False)))"),
    ('sp.radsimp((A + B)/(1 - sp.sqrt(3)))', '(-sqrt(3)*A - A - sqrt(3)*B - B)/2', "Mul(Rational(1, 2), Add(Mul(Integer(-1), Pow(Integer(3), Rational(1, 2)), Symbol('A', commutative=False)), Mul(Integer(-1), Symbol('A', commutative=False)), Mul(Integer(-1), Pow(Integer(3), Rational(1, 2)), Symbol('B', commutative=False)), Mul(Integer(-1), Symbol('B', commutative=False))))"),
    ('sp.radsimp(A/(x + sp.sqrt(2)))', '(x*A - sqrt(2)*A)/(x**2 - 2)', "Mul(Pow(Add(Pow(Symbol('x'), Integer(2)), Integer(-2)), Integer(-1)), Add(Mul(Symbol('x'), Symbol('A', commutative=False)), Mul(Integer(-1), Pow(Integer(2), Rational(1, 2)), Symbol('A', commutative=False))))"),
    ('sp.nonlinsolve([A*B - 1], [B])', '{(A**(-1),)}', "FiniteSet(Tuple(Pow(Symbol('A', commutative=False), Integer(-1))))"),
    ('sp.nonlinsolve([A*B - 1], B)', '{(A**(-1),)}', "FiniteSet(Tuple(Pow(Symbol('A', commutative=False), Integer(-1))))"),
    ('sp.nonlinsolve([B*A - 1], [B])', '{(A**(-1),)}', "FiniteSet(Tuple(Pow(Symbol('A', commutative=False), Integer(-1))))"),
    ('sp.nonlinsolve([B*A - C], [B])', '{(C*A**(-1),)}', "FiniteSet(Tuple(Mul(Symbol('C', commutative=False), Pow(Symbol('A', commutative=False), Integer(-1)))))"),
    ('sp.nonlinsolve([x*A - A], [x])', '{(1,)}', 'FiniteSet(Tuple(Integer(1)))'),
    ('sp.nonlinsolve([A - B], [A])', '{(B,)}', "FiniteSet(Tuple(Symbol('B', commutative=False)))"),
    ('sp.srepr(sp.adjoint(A*B))', "Mul(adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))", ''),
    ('sp.Matrix([[A, B], [C, A]]).H', 'Matrix([[adjoint(A), adjoint(C)], [adjoint(B), adjoint(A)]])', "MutableDenseMatrix([[adjoint(Symbol('A', commutative=False)), adjoint(Symbol('C', commutative=False))], [adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False))]])"),
    ('sp.adjoint(H)', 'H', "Symbol('H', commutative=False, hermitian=True)"),
    ('sp.transpose(H)', 'conjugate(H)', "conjugate(Symbol('H', commutative=False, hermitian=True))"),
    ('sp.adjoint(K)', '-K', ''),
    ('sp.transpose(K)', '-conjugate(K)', ''),
    ('sp.adjoint(H*A)', 'adjoint(A)*H', "Mul(adjoint(Symbol('A', commutative=False)), Symbol('H', commutative=False, hermitian=True))"),
    ('sp.conjugate(sp.transpose(A*B))', 'adjoint(B)*adjoint(A)', "Mul(adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False)))"),
    ('sp.transpose(sp.adjoint(A))', 'conjugate(A)', "conjugate(Symbol('A', commutative=False))"),
    ('sp.adjoint(sp.conjugate(A))', 'transpose(A)', "transpose(Symbol('A', commutative=False))"),
    ('sp.adjoint(A).subs(A, B)', 'adjoint(B)', "adjoint(Symbol('B', commutative=False))"),
    ('sp.adjoint(A).xreplace({A: x})', 'conjugate(x)', "conjugate(Symbol('x'))"),
    ('sp.srepr(sp.sin)', 'sin', ''),
    ('sp.srepr(sp.Function("q"))', "Function('q')", ''),
    ('sp.adjoint(sp.Piecewise((A, x > 0), (B, True)))', 'Piecewise((adjoint(A), x > 0), (adjoint(B), True))', "Piecewise(ExprCondPair(adjoint(Symbol('A', commutative=False)), StrictGreaterThan(Symbol('x'), Integer(0))), ExprCondPair(adjoint(Symbol('B', commutative=False)), true))"),
    ('sp.adjoint(x)', 'conjugate(x)', "conjugate(Symbol('x'))"),
    ('sp.transpose(x)', 'x', "Symbol('x')"),
    ('sp.Matrix([[A, B], [C, A]]).adjoint()', 'Matrix([[adjoint(A), adjoint(C)], [adjoint(B), adjoint(A)]])', "MutableDenseMatrix([[adjoint(Symbol('A', commutative=False)), adjoint(Symbol('C', commutative=False))], [adjoint(Symbol('B', commutative=False)), adjoint(Symbol('A', commutative=False))]])"),
    ('sp.Matrix([[x, sp.I]]).H', 'Matrix([[conjugate(x)], [-I]])', "MutableDenseMatrix([[conjugate(Symbol('x'))], [Mul(Integer(-1), I)]])"),
]


@pytest.mark.parametrize('code, text, representation', REFERENCE)
def test_matches_sympy(code, text, representation):
    value = eval(code, NS)
    assert str(value) == text
    if representation:
        assert sp.srepr(value) == representation


@pytest.mark.parametrize('code', [
    "sp.Symbol('A', commutative=False, real=True)",
    "sp.Symbol('A', commutative=False, positive=True)",
    "sp.Symbol('A', commutative=False, integer=True)",
    "sp.Symbol('A', commutative=False, complex=True)",
    "sp.Symbol('A', commutative=False, zero=True)",
    "sp.Symbol('A', commutative=False, nonzero=True)",
    "sp.Dummy('D', commutative=False, real=True)",
    "sp.Function('h', commutative=False, real=True)",
])
def test_inconsistent_assumptions(code):
    with pytest.raises(InconsistentAssumptions):
        eval(code, NS)
    # SymPy's InconsistentAssumptions derives from ValueError
    with pytest.raises(ValueError):
        eval(code, NS)


def test_adjoint_nodes_keep_identity():
    e = sp.adjoint(A)*sp.transpose(B)
    assert pickle.loads(pickle.dumps(e)) == e
    assert sp.adjoint(A*B) == sp.adjoint(B)*sp.adjoint(A)
    assert sp.adjoint(A*B).free_symbols == {A, B}
    assert sp.conjugate(sp.adjoint(A)) == sp.transpose(A)
    assert sp.conjugate(sp.transpose(A)) == sp.adjoint(A)


def test_nonlinsolve_intentional_corrections():
    # SymPy 1.14 returns {(C*A**(-1),)}, which does not satisfy A*B = C.
    sol = sp.nonlinsolve([A*B - C], [B])
    assert sol == sp.FiniteSet(sp.Tuple(A**-1*C))
    (b,), = sol
    assert sp.expand(A*b - C) == 0
    # SymPy 1.14 gives a noncommutative value to the commutative x;
    # consistent with solve(A*x - B, x) == [].
    assert sp.solve(A*x - B, x) == []
    assert sp.nonlinsolve([A*x - B], [x]) == sp.EmptySet
    # SymPy 1.14 claims every complex A commutes with B.
    with pytest.raises(NotImplementedError):
        sp.nonlinsolve([A*B - B*A], [A])
    # SymPy 1.14 crashes with AttributeError.
    with pytest.raises(NotImplementedError):
        sp.nonlinsolve([A*B - 1, x - 1], [B, x])


def test_function_srepr_round_trips_facts():
    # Intentional: SymPy drops the facts of applied undefined functions,
    # symox keeps them so srepr round-trips (pickling, to_sympy).
    for e in (g(x), g(A), sp.Function('h', real=True)(x)):
        text = sp.srepr(e)
        assert text.startswith("Function('")
        assert sp.sympify(text) == e
        assert sp.sympify(text).is_commutative == e.is_commutative
