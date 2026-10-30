"""The native ``unrad`` (``symox::solve::unrad``) and the native radical and
negative-base steps of ``solveset`` (``symox::sets::solveset_radical``),
against values computed with SymPy 1.14 (hard-coded below)."""
from symox import *  # noqa: F403
from symox.solvers.solvers import unrad
import pytest

del test  # noqa: F821

x, y, t, a, b, c, d = symbols('x y t a b c d')
s = Symbol('s')


def _norm(rv):
    """``(eq, cov)`` with every Dummy renamed ``s`` and ``eq`` expanded."""
    if rv is None:
        return 'None'
    eq, cov = rv
    ds = set(eq.atoms(Dummy))
    for z in cov:
        ds |= z.atoms(Dummy)
    rep = {dd: s for dd in ds}
    return '%s ;; %s' % (expand(eq.xreplace(rep)), [z.xreplace(rep) for z in cov])


# SymPy 1.14 test_unrad1 / test_unrad_slow and solve's unrad calls:
# str of the normalized result, or the error class and message.
UNRAD_CASES = [
    ('sqrt(x) + sqrt(x + 1) + sqrt(1 - sqrt(x)) + 3',
     'ERR NotImplementedError: Cannot remove all radicals'),
    ('sqrt(x) + (x + 1)**Rational(1, 3) + 2*sqrt(y)',
     'ERR NotImplementedError: Cannot remove all radicals'),
    ('S(1)',
     'None'),
    ('sqrt(x)',
     'x ;; []'),
    ('sqrt(x) + 1',
     'x - 1 ;; []'),
    ('sqrt(x) + root(x, 3) + 2',
     's**3 + s**2 + 2 ;; [s, s**6 - x]'),
    ('sqrt(x)*root(x, 3) + 2',
     'x**5 - 64 ;; []'),
    ('sqrt(x) + (x + 1)**Rational(1, 3)',
     '-x**3 + x**2 + 2*x + 1 ;; []'),
    ('sqrt(x) + sqrt(x + 1) + sqrt(2*x)',
     '2*x + 2*sqrt(2)*x - 1 ;; []'),
    ('sqrt(x) + sqrt(x + 1) + 2',
     '16*x - 9 ;; []'),
    ('sqrt(x) + sqrt(x + 1) + sqrt(1 - x)',
     '5*x**2 - 4*x ;; []'),
    ('a*sqrt(x) + b*sqrt(x) + c*sqrt(y) + d*sqrt(y)',
     'a**2*x + 2*a*b*x + b**2*x - c**2*y - 2*c*d*y - d**2*y ;; []'),
    ('sqrt(x) + sqrt(1 - x)',
     '2*x - 1 ;; []'),
    ('sqrt(x) + sqrt(1 - x) - 3',
     'x**2 - x + 16 ;; []'),
    ('sqrt(x) + sqrt(1 - x) + sqrt(2 + x)',
     '5*x**2 - 2*x + 1 ;; []'),
    ('sqrt(x) + sqrt(1 - x) + sqrt(2 + x) - 3',
     '25*x**4 + 376*x**3 + 1256*x**2 - 2272*x + 784 ;; []'),
    ('sqrt(x) + sqrt(1 - x) + sqrt(2 + x) - sqrt(1 - 2*x)',
     '41*x**4 + 40*x**3 + 232*x**2 - 160*x + 16 ;; []'),
    ('sqrt(x) + sqrt(x + 1)',
     '1 ;; []'),
    ('sqrt(x) + sqrt(x + 1) + sqrt(1 - sqrt(x))',
     '16*x**2 - 9*x ;; []'),
    ('(sqrt(x) + root(x + 1, 3) + 2*sqrt(y), y)',
     '2*sqrt(x)*(x + 1)**(1/3) + x - 4*y + (x + 1)**(2/3) ;; []'),
    ('sqrt(x/(1 - x)) + (x + 1)**Rational(1, 3)',
     'x**5 - x**4 - x**3 + 2*x**2 + x - 1 ;; []'),
    ('(sqrt(x/(1 - x)) + 2*sqrt(y), y)',
     '4*x*y + x - 4*y ;; []'),
    ('(sqrt(x)*sqrt(1 - x) + 2, x)',
     'x**2 - x + 4 ;; []'),
    ('root(x + 1, 5) - root(x, 3)',
     '-x**5 + x**3 + 3*x**2 + 3*x + 1 ;; []'),
    ('(x + root(x, 3) + root(x, 3)**2 + sqrt(y), x)',
     's**3 + s**2 + s + sqrt(y) ;; [s, s**3 - x]'),
    ('sqrt(x) + root(x, 3) + y',
     's**3 + s**2 + y ;; [s, s**6 - x]'),
    ('sqrt(x) - root(x + 1, 3)*sqrt(x + 2) + 2',
     '-s**10 - 8*s**8 - 24*s**6 + 12*s**5 + 22*s**4 + 160*s**3 + 212*s**2 + 192*s + 56 ;; [s, s**2 - x]'),
    ('root(x + 1, 3) + root(x, 3)',
     '2*x + 1 ;; []'),
    ('sqrt(x) + sqrt(x + 1) + sqrt(1 - x) - 6*sqrt(5)/5',
     '15625*x**4 + 173000*x**3 + 355600*x**2 - 817920*x + 331776 ;; []'),
    ('root(x, 4) + root(x, 4)**3 - 1',
     's**3 + s - 1 ;; [s, s**4 - x]'),
    ('root(x, 2) + root(x, 2)**3 - 1',
     'x**3 + 2*x**2 + x - 1 ;; []'),
    ('x**0.5',
     'None'),
    ('t + root(x + y, 5) + root(x + y, 5)**3',
     's**3 + s + t ;; [s, s**5 - x - y]'),
    ('(x + root(x + y, 5) + root(x + y, 5)**3, y)',
     's**3 + s + x ;; [s, s**5 - x - y]'),
    ('(x + root(x + y, 5) + root(x + y, 5)**3, x)',
     's**5 + s**3 + s - y ;; [s, s**5 - x - y]'),
    ('root(x - 1, 3) + root(x + 1, 5) + root(2, 5)',
     's**5 + 5*2**(1/5)*s**4 + s**3 + 10*2**(2/5)*s**3 + 10*2**(3/5)*s**2 + 5*2**(4/5)*s + 4 ;; [s, s**3 - x + 1]'),
    ('(root(x, 2) + root(x, 3) + root(x, 4)).subs(x, x**5 - x + 1)',
     'ERR NotImplementedError: Cannot remove all radicals'),
    ('sqrt(x + root(x + 1, 3)) - root(x + 1, 3) - 2',
     '-s**3 + s**2 + 3*s + 5 ;; [s, s**3 - x - 1]'),
    ('root(x**2 + 1, 3) - root(x**2 - 1, 5) - 2',
     's**5 - 10*s**4 + 39*s**3 - 80*s**2 + 80*s - 30 ;; [s, s**3 - x**2 - 1]'),
    ('sqrt(x + root(x + 1, 2)) - root(x + 1, 3) - 2',
     's**6 - 2*s**5 - 7*s**4 - 3*s**3 + 26*s**2 + 40*s + 25 ;; [s, s**3 - x - 1]'),
    ('sqrt(x + sqrt(root(x, 3) - 1)) - root(x, 6) - 2',
     's**12 - 2*s**8 - 8*s**7 - 8*s**6 + s**4 + 8*s**3 + 23*s**2 + 32*s + 17 ;; [s, s**6 - x]'),
    ('root(cosh(x), 3)/x*root(x + 1, 5) - 1',
     '-x**15 + x**3*cosh(x)**5 + 3*x**2*cosh(x)**5 + 3*x*cosh(x)**5 + cosh(x)**5 ;; []'),
    ("S('(x+y)**(2*y/3) + (x+y)**(1/3) + 1')",
     'None'),
    ("(S('(x+y)**(2*y/3) + (x+y)**(1/3) + 1'), x)",
     's + s**(2*y) + 1 ;; [s, s**3 - x - y]'),
    ('(x**(S.Half/y) + y, x)',
     'x**(1/y) - y**2 ;; []'),
    ('root(x + 1, 3) - (root(x, 3) + root(x, 5))',
     '3*s**13 + 3*s**11 + s**9 - 1 ;; [s, s**15 - x]'),
    ('root(x + 1, 3) - (root(x, 3) + root(x, 5)) - 2',
     '3*s**13 + 3*s**11 + 6*s**10 + s**9 + 12*s**8 + 6*s**6 + 12*s**5 + 12*s**3 + 7 ;; [s, s**15 - x]'),
    ('root(x, 3) - root(x + 1, 4)/2 + root(x + 2, 3)',
     '4096*s**10 + 960*s**9 + 48*s**8 - s**7 - 1728*s ;; [s, s**4 - x - 1]'),
    ('root(x, 3) + root(x + 1, 4) - root(x + 2, 3)/2',
     '343*s**13 + 2904*s**12 + 1344*s**11 + 512*s**10 - 1323*s**9 - 3024*s**8 - 1728*s**7 + 1701*s**5 + 216*s**4 - 729*s ;; [s, s**4 - x - 1]'),
    ('root(x, 3)/2 - root(x + 1, 4) + root(x + 2, 3)',
     '729*s**13 - 216*s**12 + 1728*s**11 - 512*s**10 + 1701*s**9 - 3024*s**8 + 1344*s**7 + 1323*s**5 - 2904*s**4 + 343*s ;; [s, s**4 - x - 1]'),
    ('root(x, 3)/2 - root(x + 1, 4) + root(x + 2, 3) - 2',
     '-729*s**13 - 1242*s**12 - 18496*s**10 - 129701*s**9 - 388602*s**8 - 453312*s**7 + 612864*s**6 + 3337173*s**5 + 6332418*s**4 + 7134912*s**3 + 5064768*s**2 + 2111913*s + 398034 ;; [s, s**4 - x - 1]'),
    ('root(x, 3) - root(y, 3) + root(x, 5)',
     's**15 + 3*s**13 + 3*s**11 + s**9 - y ;; [s, s**15 - x]'),
    ('root(x, 3) + root(y, 3) + root(x*y, 4)',
     's**13*y + 3*s**12*y**2 + 3*s**11*y**3 + s**10*y**4 + 3*s**9*y**3 - 21*s**8*y**4 + 3*s**7*y**5 + 3*s**5*y**5 + 3*s**4*y**6 + s*y**7 ;; [s, s**4 - x*y]'),
    ('root(x, 3) + root(y, 3) + root(x*y, 5)',
     'ERR NotImplementedError: Cannot remove all radicals'),
    ('Eq(-x**(S(1)/5) + x**(S(1)/3), -3**(S(1)/3) - (-1)**(S(3)/5)*3**(S(1)/5))',
     's**5 - s**3 + 3**(1/3) + (-1)**(3/5)*3**(1/5) ;; [s, s**15 - x]'),
    ('(sqrt(x) - 1)**2 - (sqrt(x) - 1)**3',
     'x**3 - 6*x**2 + 9*x - 4 ;; []'),
    ('((x/(x + 1) + 3)**(-2), x)',
     'None'),
    ('sqrt(x - y)*exp(t*sqrt(x - y)) - exp(t*sqrt(x - y))',
     'None'),
    ('sqrt(1 + sqrt(1 - 4*x**2)) - x*(1 + sqrt(1 + 2*sqrt(1 - 4*x**2)))',
     '256*x**19 - 1280*x**17 + 3168*x**15 - 4752*x**13 + 4737*x**11 - 3220*x**9 + 1476*x**7 - 418*x**5 + 64*x**3 - 4*x ;; []'),
    ('sqrt(x) + sqrt(x + 1) - sqrt(1 - x) - sqrt(2 + x)',
     '5*x**4 + 4*x**3 - 2*x**2 - 4*x + 1 ;; []'),
    ('sqrt(x + 1) + root(x, 3) - 2',
     's**3 - s**2 + 4*s - 3 ;; [s, s**3 - x]'),
    ('sqrt(sqrt(x + 1)) + x**Rational(1, 3) - 2',
     's**4 - 9*s**3 + 24*s**2 - 32*s + 15 ;; [s, s**3 - x]'),
    ('sqrt(x) + root(sqrt(x) + 1, 3) - 2',
     's**3 - 6*s**2 + 13*s - 7 ;; [s, s**2 - x]'),
    ('sqrt(2*x - 1) - sqrt(x - 4) - 2',
     'x**2 - 18*x + 65 ;; []'),
    ('x + 1 - root(x**4 + 4*x**3 - x, 4)',
     '6*x**2 + 5*x + 1 ;; []'),
    ('sqrt(17*x - sqrt(x**2 - 5)) - 7',
     '144*x**2 - 833*x + 1203 ;; []'),
    ('sqrt(x) - sqrt(x - 1) + sqrt(sqrt(x))',
     '16*x**3 - 17*x**2 + 2*x - 1 ;; []'),
    ('root(x, 3) + root(x, 5) - 2',
     's**5 + s**3 - 2 ;; [s, s**15 - x]'),
    ('sqrt(x) + root(x, 3) + root(x + 1, 5) - 2',
     'ERR NotImplementedError: Cannot remove all radicals'),
]


@pytest.mark.parametrize('src, expected', UNRAD_CASES)
def test_unrad_matches_sympy(src, expected):
    arg = eval(src)
    try:
        out = _norm(unrad(*arg) if isinstance(arg, tuple) else unrad(arg))
    except Exception as e:  # noqa: BLE001
        out = 'ERR %s: %s' % (type(e).__name__, str(e)[:60])
    assert out == expected


def test_unrad_equality_and_non_expr():
    assert unrad(Eq(sqrt(x), 2)) == unrad(sqrt(x) - 2)
    assert unrad([sqrt(x)]) is None


def test_unrad_recursion_flags():
    # SymPy 1.14 accepts the recursion flags cov, n and rpt
    c = []
    assert str(unrad(sqrt(x) + 1, x, cov=c)) == '(x - 1, [])' and c == []


def test_unrad_take():
    # SymPy 1.14: a custom _take selects the generators; a taken generator
    # that is not a power has no ``exp``
    assert str(unrad(sqrt(x) + root(y, 3) - 2, _take=lambda d: d == root(y, 3) or d.is_Symbol)) == \
        '(x**(3/2) + 12*sqrt(x) - 6*x + y - 8, [])'
    with pytest.raises(AttributeError, match="'Symbol' object has no attribute 'exp'"):
        unrad(sqrt(x) + root(y, 3) - 2, _take=lambda d: d.has(y))
    with pytest.raises(AttributeError, match="'Symbol' object has no attribute 'exp'"):
        unrad(sqrt(x) + sqrt(x + 1) - y, _take=lambda d: d.has(x))


@pytest.mark.parametrize('src, expected', [
    # SymPy 1.14
    ("solveset(sqrt(x) + sqrt(x + 1) - 3, x, S.Reals)", "{16/9}"),
    ("solveset(sqrt(x) + sqrt(x + 1) - 3, x)", "{16/9}"),
    ("solveset(sqrt(x) - sqrt(x - 1) - 2, x, S.Reals)", "EmptySet"),
    ("solveset(sqrt(2*x - 1) - sqrt(x - 4) - 2, x, S.Reals)", "{5, 13}"),
    ("solveset(sqrt(x + 1) + root(x, 3) - 2, x, S.Reals)",
     "{(-11/(9*(47/54 + sqrt(93)/6)**(1/3)) + 1/3 + (47/54 + sqrt(93)/6)**(1/3))**3}"),
    ("solveset(sqrt(x) + root(x, 3) - 2, x, S.Reals)", "{1}"),
    # issue 17882: unrad gives None, then _solve_as_rational
    ("solveset(-8*x**2/(9*(x**2 - 1)**(S(4)/3)) + 4/(3*(x**2 - 1)**(S(1)/3)), x, S.Complexes)",
     "{-sqrt(3), sqrt(3)}"),
    ("solveset(sqrt(x**2 - 1) - x + 2, x, S.Reals)", "EmptySet"),
    ("solveset(x - sqrt(x + 6), x, S.Reals)", "{3}"),
    ("solveset(sqrt(x - 1) + sqrt(x + 2) - 3, x, S.Complexes)", "{2}"),
    ("solveset(1/sqrt(x) + sqrt(x) - 2, x, S.Reals)", "{1}"),
    ("solveset(sqrt(x) + 2, x, S.Reals)", "EmptySet"),
    ("solveset(sqrt(x)*3 - 6, x, S.Reals)", "{4}"),
    ("solveset((-2)**x - 4, x, S.Reals)", "{2}"),
    ("solveset((-2)**x + 8, x, S.Reals)", "{3}"),
    ("solveset((-2)**x - 2, x, S.Reals)", "EmptySet"),
    ("solveset(3*(-2)**x - 12, x, S.Reals)", "{2}"),
    ("solveset((-2)**(2*x) - 16, x, S.Reals)", "{2}"),
    ("solveset((-2)**x - 4, x, S.Complexes)", "ConditionSet(x, Eq((-2)**x - 4, 0), Complexes)"),
    ("solve(sqrt(x) + sqrt(x + 1) - 3, x)", "[16/9]"),
    ("solve(sqrt(x) + root(x, 3) - 2, x)", "[1]"),
    ("solve(root(x + 1, 3) + root(x, 3) - 1, x)", "[0]"),
    ("solve(sqrt(x - 1) - x + 3, x)", "[5]"),
    # intentional corrections (docs/algorithmic-parity.md): SymPy raises
    # ValueError from integer_log for these bases ...
    ("solveset((-1)**x - 1, x, S.Reals)", "ImageSet(Lambda(_n, 2*_n), Integers)"),
    ("solveset((-1)**x + 1, x, S.Reals)", "ImageSet(Lambda(_n, 2*_n + 1), Integers)"),
    ("solveset((-S(1)/2)**x - 4, x, S.Reals)", "{-2}"),
    # ... and returns EmptySet although x = 0 is a solution
    ("solveset((-3)**x - 2**x, x, S.Reals)", "ConditionSet(x, Eq((-3)**x - 2**x, 0), Reals)"),
])
def test_steps_match_sympy(src, expected):
    assert str(eval(src)) == expected


def test_radical_step_keeps_undecided_in_conditionset():
    from symox.solvers._solveset_steps import radical_step
    # not a radical equation: the native result is kept
    assert radical_step(x**2 - 1, x, S.Reals, solveset) is None
    # a single radical with a negative right side has no real solution
    assert radical_step(sqrt(x) + 2, x, S.Reals, solveset) == S.EmptySet


def test_negative_base_guard():
    from symox.solvers._solveset_steps import negative_base_guard
    p = Symbol('p', positive=True)
    assert negative_base_guard(p**x - 4, x, S.Reals, FiniteSet(1)) == FiniteSet(1)
    assert negative_base_guard(y**x - 4, x, S.Reals, FiniteSet(1)) ==         ConditionSet(x, Eq(y**x - 4, 0), S.Reals)
    assert negative_base_guard(y**x - 4, x, S.Complexes, FiniteSet(1)) == FiniteSet(1)


def test_radical_step_resultant_fallback():
    # SymPy 1.14: unrad raises NotImplementedError('Cannot remove all
    # radicals') and so does solveset. The established resultant elimination
    # keeps the checked root 0.
    eq = sqrt(x) + sqrt(x + 1) + sqrt(x + 4) + sqrt(x + 9) - 6
    with pytest.raises(NotImplementedError, match='Cannot remove all radicals'):
        unrad(eq)
    assert solveset(eq, x) == FiniteSet(0)
    assert solveset(eq, x, S.Reals) == FiniteSet(0)
    assert solve(eq, x) == [0]


def test_radical_step_keeps_parametric_candidate():
    # SymPy 1.14 (test_issue_11534): the candidate -y/sqrt(y**2 + 1) is
    # valid at y = 0 and stays in the ConditionSet.
    xr, yr = symbols('xr yr', real=True)
    r = solveset(-yr + xr/sqrt(-xr**2 + 1), xr, S.Reals)
    assert r == Union(Complement(ConditionSet(xr, Eq(xr/sqrt(1 - xr**2) - yr, 0),
                                              FiniteSet(-yr/sqrt(yr**2 + 1))), FiniteSet(-1, 1)),
                      Complement(FiniteSet(yr/sqrt(yr**2 + 1)), FiniteSet(-1, 1)))
