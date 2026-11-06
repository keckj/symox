"""Native is_log_deriv_k_t_radical_in_field and parametric_log_deriv_heu on
exp/primitive towers (symox::integrate::prde_tower).

The cases are seeded random towers
(audits/performance/2026-11-05-prde/scripts/random_towers.py); the expected
outputs (srepr, exception type and message, number of Dummies created by a
repeated call) are hard-coded from SymPy 1.14 (sympy.integrals.prde).
"""
import pytest

from symox import Symbol, Poly, Dummy, srepr, sympify, symbols, QQ
from symox.integrals.risch import DifferentialExtension
from symox.integrals.prde import (is_log_deriv_k_t_radical_in_field,
                                  parametric_log_deriv_heu)
from symox.integrals import _prde_native

# [D, case (None: parametric_log_deriv_heu), args, expected, dummies]
CASES = [
    [['1', '2*t0'], 'auto', ['3*t0**2 - t0', '-t0**2 + t0*x - 2'], None, 1],
    [['1', '1/x', '1/(x*t0)'], 'tan', ['-t0*t1', '2*t1*x/(x + 1) - 1'], None, 1],
    [['1', '(1 + t0**2)'], 'primitive', ['2*x', '3*t0'], None, 1],
    [['1', '1/x'], 'base', ['t0 + 2*x', 't0**2 + 2*t0 - 2'], None, 1],
    [['1', 't0', '2*x*t1'], 'auto', ['-2*t0*t1 - 2*t0', '2'], None, 2],
    [['1', 't0', '1/x'], 'auto', ['t0*t1*x + t0*x - 1', 't1*x + x'], None, 5],
    [['1', '1/(x + 1)'], 'exp', ['1', '-x'], None, 3],
    [['1', 't0', '1/x'], 'auto', ['2*x + 1', 'x + 1'], None, 5],
    [['1', '1/x'], 'auto', ['-3*t0*x + 3*x - 1', '3*t0*x - 3*x'], None, 5],
    [['1', 't0', 't1/x'], 'exp', ['t0**2 - x + 1', 't0**2 + x'], None, 3],
    [['1', 'x*t0'], 'auto', ['-1', 'x**2 + x'], "(Integer(1), Mul(Pow(Symbol('x'), Integer(-1)), Add(Symbol('x'), Integer(1))))", 6],
    [['1', 'x*t0'], 'auto', ['9*x - 3', '3*x**2 - 2*x'], "(Integer(2), Add(Pow(Symbol('x'), Integer(6)), Mul(Integer(-1), Integer(2), Pow(Symbol('x'), Integer(5))), Mul(Rational(4, 3), Pow(Symbol('x'), Integer(4))), Mul(Integer(-1), Rational(8, 27), Pow(Symbol('x'), Integer(3)))))", 6],
    [['1', '1/(x + 1)'], 'auto', ['1/(x + 1)', '-2'], "(Integer(2), Pow(Add(Symbol('x'), Integer(1)), Integer(-1)))", 5],
    [['1', '1/x'], 'auto', ['t0 + x + 1', 't0*x**2 + t0*x'], "(Integer(1), Mul(Symbol('t0'), Symbol('x'), Pow(Add(Symbol('x'), Integer(1)), Integer(-1))))", 6],
    [['1', '1/(x + 1)'], 'primitive', ['0', '1'], '(Integer(1), Integer(1))', 4],
    [['1', '1/x'], 'auto', ['3', 'x + 1'], "(Integer(1), Add(Pow(Symbol('x'), Integer(3)), Mul(Integer(3), Pow(Symbol('x'), Integer(2))), Mul(Integer(3), Symbol('x')), Integer(1)))", 5],
    [['1', '2*t0'], 'base', ['0', '1'], '(Integer(1), Integer(1))', 2],
    [['1', '1/x', '1/(x*t0)'], 'auto', ['12*t0**3 + 24*t0**2 + 12*t0 + 3', '4*t0**4*x + 8*t0**3*x + 4*t0**2*x - 2*t0*x - 2*x'], "(Integer(2), Mul(Pow(Add(Mul(Integer(8), Pow(Symbol('t0'), Integer(3))), Mul(Integer(24), Pow(Symbol('t0'), Integer(2))), Mul(Integer(24), Symbol('t0')), Integer(8)), Integer(-1)), Add(Mul(Integer(8), Pow(Symbol('t0'), Integer(9))), Mul(Integer(24), Pow(Symbol('t0'), Integer(8))), Mul(Integer(24), Pow(Symbol('t0'), Integer(7))), Mul(Integer(-1), Integer(4), Pow(Symbol('t0'), Integer(6))), Mul(Integer(-1), Integer(24), Pow(Symbol('t0'), Integer(5))), Mul(Integer(-1), Integer(12), Pow(Symbol('t0'), Integer(4))), Mul(Integer(6), Pow(Symbol('t0'), Integer(3))), Mul(Integer(6), Pow(Symbol('t0'), Integer(2))), Integer(-1))))", 8],
    [['1', '1/x', '1/(x*t0)'], 'base', ['0', 't0/(t0 + 2) + 2*t1**2'], '(Integer(1), Integer(1))', 2],
    [['1', '1/x'], 'auto', ['-6*t0**3*x - 5*t0**2*x**2 - 10*t0**2*x + 2*t0*x**2 + 2*t0*x + 8*t0 - 3*x**2 - 6*x', '3*t0**3*x**3 + 6*t0**3*x**2 - 3*t0**2*x**3 - 12*t0**2*x**2 - 12*t0**2*x + 3*t0*x**3 + 6*t0*x**2 - 3*x**3 - 12*x**2 - 12*x'], "(Integer(3), Mul(Add(Pow(Symbol('x'), Integer(3)), Mul(Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(12), Symbol('x')), Integer(8)), Pow(Add(Mul(Pow(Symbol('t0'), Integer(5)), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(3), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(6), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(2))), Mul(Integer(4), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(3))), Mul(Integer(12), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(2))), Mul(Integer(12), Pow(Symbol('t0'), Integer(3)), Symbol('x')), Mul(Integer(-1), Integer(4), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(12), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(12), Pow(Symbol('t0'), Integer(2)), Symbol('x')), Mul(Integer(-1), Integer(8), Pow(Symbol('t0'), Integer(2))), Mul(Integer(3), Symbol('t0'), Pow(Symbol('x'), Integer(3))), Mul(Integer(12), Symbol('t0'), Pow(Symbol('x'), Integer(2))), Mul(Integer(12), Symbol('t0'), Symbol('x')), Mul(Integer(-1), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(6), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(12), Symbol('x')), Integer(-8)), Integer(-1))))", 7],
    [['1', '1/x', '1/(x*t0)'], 'auto', ['0', '1'], '(Integer(1), Integer(1))', 6],
    [['1', 'x*t0'], 'auto', ['6*t0*x**2 + 2*t0 - 3', '12*t0*x - 18*x'], "(Integer(6), Add(Mul(Pow(Symbol('t0'), Integer(3)), Symbol('x')), Mul(Integer(-1), Rational(9, 2), Pow(Symbol('t0'), Integer(2)), Symbol('x')), Mul(Rational(27, 4), Symbol('t0'), Symbol('x')), Mul(Integer(-1), Rational(27, 8), Symbol('x'))))", 7],
    [['1', '1/(x + 1)'], 'auto', ['5*t0**2*x**3 + 10*t0**2*x**2 + 6*t0**2*x + 3*t0**2 + 6*t0*x**3 + 8*t0*x**2 + 10*t0*x - 5*x**3 - 7*x**2 - 3*x', '2*t0**2*x**4 + 4*t0**2*x**3 + 4*t0**2*x**2 + 2*t0**2*x + 2*t0*x**3 + 4*t0*x**2 + 4*t0*x + 2*t0 - 2*x**4 - 4*x**3 - 4*x**2 - 2*x'], "(Integer(2), Mul(Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1)), Integer(-1)), Add(Mul(Pow(Symbol('t0'), Integer(6)), Pow(Symbol('x'), Integer(7))), Mul(Integer(2), Pow(Symbol('t0'), Integer(6)), Pow(Symbol('x'), Integer(6))), Mul(Integer(3), Pow(Symbol('t0'), Integer(6)), Pow(Symbol('x'), Integer(5))), Mul(Integer(2), Pow(Symbol('t0'), Integer(6)), Pow(Symbol('x'), Integer(4))), Mul(Pow(Symbol('t0'), Integer(6)), Pow(Symbol('x'), Integer(3))), Mul(Integer(3), Pow(Symbol('t0'), Integer(5)), Pow(Symbol('x'), Integer(6))), Mul(Integer(6), Pow(Symbol('t0'), Integer(5)), Pow(Symbol('x'), Integer(5))), Mul(Integer(9), Pow(Symbol('t0'), Integer(5)), Pow(Symbol('x'), Integer(4))), Mul(Integer(6), Pow(Symbol('t0'), Integer(5)), Pow(Symbol('x'), Integer(3))), Mul(Integer(3), Pow(Symbol('t0'), Integer(5)), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(3), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(7))), Mul(Integer(-1), Integer(6), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(6))), Mul(Integer(-1), Integer(6), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(5))), Mul(Integer(6), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(3))), Mul(Integer(6), Pow(Symbol('t0'), Integer(4)), Pow(Symbol('x'), Integer(2))), Mul(Integer(3), Pow(Symbol('t0'), Integer(4)), Symbol('x')), Mul(Integer(-1), Integer(6), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(6))), Mul(Integer(-1), Integer(12), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(5))), Mul(Integer(-1), Integer(17), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(4))), Mul(Integer(-1), Integer(10), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(3), Pow(Symbol('t0'), Integer(3)), Pow(Symbol('x'), Integer(2))), Mul(Integer(2), Pow(Symbol('t0'), Integer(3)), Symbol('x')), Pow(Symbol('t0'), Integer(3)), Mul(Integer(3), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(7))), Mul(Integer(6), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(6))), Mul(Integer(6), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(5))), Mul(Integer(-1), Integer(6), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(3))), Mul(Integer(-1), Integer(6), Pow(Symbol('t0'), Integer(2)), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Integer(3), Pow(Symbol('t0'), Integer(2)), Symbol('x')), Mul(Integer(3), Symbol('t0'), Pow(Symbol('x'), Integer(6))), Mul(Integer(6), Symbol('t0'), Pow(Symbol('x'), Integer(5))), Mul(Integer(9), Symbol('t0'), Pow(Symbol('x'), Integer(4))), Mul(Integer(6), Symbol('t0'), Pow(Symbol('x'), Integer(3))), Mul(Integer(3), Symbol('t0'), Pow(Symbol('x'), Integer(2))), Mul(Integer(-1), Pow(Symbol('x'), Integer(7))), Mul(Integer(-1), Integer(2), Pow(Symbol('x'), Integer(6))), Mul(Integer(-1), Integer(3), Pow(Symbol('x'), Integer(5))), Mul(Integer(-1), Integer(2), Pow(Symbol('x'), Integer(4))), Mul(Integer(-1), Pow(Symbol('x'), Integer(3))))))", 7],
    [['1', '2*t0'], 'primitive', ['1', '2*x'], "(Integer(2), Symbol('x'))", 5],
    [['1', '1/x'], 'auto', ['-4*t0*x + 2*t0 + 7*x - 5', '6*t0*x**2 - 6*t0*x - 6*x**2 + 6*x'], "(Integer(6), Mul(Pow(Add(Pow(Symbol('x'), Integer(4)), Mul(Integer(-1), Integer(2), Pow(Symbol('x'), Integer(3))), Pow(Symbol('x'), Integer(2))), Integer(-1)), Add(Pow(Symbol('t0'), Integer(3)), Mul(Integer(-1), Integer(3), Pow(Symbol('t0'), Integer(2))), Mul(Integer(3), Symbol('t0')), Integer(-1))))", 6],
    [['1', '1/x', '1/(x*t0)'], 'primitive', ['-t0*t1 - 1', '3*t0*t1*x'], "(Integer(3), Mul(Pow(Symbol('t1'), Integer(-1)), Pow(Symbol('x'), Integer(-1))))", 8],
    [['1', '1/x'], 'auto', ['3', '2*t0*x - x'], "(Integer(2), Add(Pow(Symbol('t0'), Integer(3)), Mul(Integer(-1), Rational(3, 2), Pow(Symbol('t0'), Integer(2))), Mul(Rational(3, 4), Symbol('t0')), Rational(-1, 8)))", 5],
    [['1', '1/x', '1/(x*t0)'], 'auto', ['t0*t1**2*x**2 + 4*t0*t1**2*x + 4*t0*t1**2 + 2*t0*t1*x**2 + 9*t0*t1*x + 8*t0*t1 - x - 2', '2*t0**2*t1**2*x**3 + 8*t0**2*t1**2*x**2 + 8*t0**2*t1**2*x - 2*t0*t1*x**2 - 4*t0*t1*x'], "(Integer(2), Mul(Pow(Add(Symbol('x'), Integer(2)), Integer(-1)), Add(Mul(Symbol('t0'), Pow(Symbol('t1'), Integer(2)), Symbol('x')), Mul(Integer(2), Symbol('t0'), Pow(Symbol('t1'), Integer(2))), Mul(Integer(-1), Symbol('t1')))))", 8],
    [['1', 'x*t0'], 'auto', ['6*t0*x**3 + 12*t0*x**2 + 6*t0*x - 6*x**2 - 6*x + 3', '2*t0*x**2 + 4*t0*x + 2*t0 - 3*x - 3'], "(Integer(1), Mul(Pow(Add(Mul(Integer(2), Symbol('x')), Integer(2)), Integer(-1)), Add(Mul(Integer(2), Pow(Symbol('t0'), Integer(3)), Symbol('x')), Mul(Integer(2), Pow(Symbol('t0'), Integer(3))), Mul(Integer(-1), Integer(3), Pow(Symbol('t0'), Integer(2))))))", 6],
    [['1', '2*t0'], 'primitive', ['-1', '-x'], "(Integer(1), Symbol('x'))", 5],
    [['1', 'x*t0'], 'auto', ['0', '3*t0'], '(Integer(1), Integer(1))', 5],
    [['1', '1/x', '1/(x*t0)'], 'auto', ['-2', '6*t0*x - 3*x'], "(Integer(3), Mul(Integer(2), Pow(Add(Mul(Integer(2), Symbol('t0')), Integer(-1)), Integer(-1))))", 7],
    [['1', '(1 + t0**2)'], 'primitive', ['0', '1'], '(Integer(1), Integer(1))', 4],
    [['1', '(1 + t0**2)'], 'base', ['0', '-t0'], '(Integer(1), Integer(1))', 2],
    [['1', '1/x', '1/(x*t0)'], 'auto', ['t0*t1 + t0 + 1', '2*t0*t1*x + 2*t0*x'], "(Integer(2), Add(Mul(Symbol('t1'), Symbol('x')), Symbol('x')))", 8],
    [['1', '(1 + t0**2)'], 'primitive', ['1', '3*x**2 + 4*x + 1'], "(Integer(2), Mul(Add(Mul(Integer(3), Symbol('x')), Integer(1)), Pow(Add(Mul(Integer(3), Symbol('x')), Integer(3)), Integer(-1))))", 5],
    [['1', 'x*t0'], 'auto', ['t0*x', '2*t0 + 2'], "(Integer(2), Add(Symbol('t0'), Integer(1)))", 6],
    [['1', '1/(x + 1)'], 'auto', ['1', '2*x'], "(Integer(2), Symbol('x'))", 5],
    [['1', 'x*t0'], 'auto', ['6*t0**2*x**3 + 24*t0**2*x**2 + 24*t0**2*x - 2', '3*t0**2*x**2 + 12*t0**2*x + 12*t0**2 + 2*x + 4'], "(Integer(1), Mul(Pow(Add(Mul(Integer(3), Symbol('x')), Integer(6)), Integer(-1)), Add(Mul(Integer(3), Pow(Symbol('t0'), Integer(2)), Symbol('x')), Mul(Integer(6), Pow(Symbol('t0'), Integer(2))), Integer(2))))", 6],
    [['1', '1/x'], 'auto', ['2*x + 4', 'x**2 + x'], "(Integer(1), Mul(Pow(Symbol('x'), Integer(4)), Pow(Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1)), Integer(-1))))", 5],
    [['1', '1/(x + 1)'], 'tan', ['0', '-t0**2/(x + 2) - t0 + 2'], ['raise', 'NotImplementedError', 'The hypertangent case is not yet implemented for is_log_deriv_k_t_radical_in_field()'], 2],
    [['1', '1/x', '1/(x*t0)'], 'bogus', ['2*t0*t1*x - 4*t1 + 2*x', '-t0*t1**2*x + t0*t1*x**2 - t0*x'], ['raise', 'ValueError', "case must be one of {'primitive', 'exp', 'tan', 'base', 'auto'}, not bogus"], 3],
    [['1', '1/(x + 1)'], 'other_linear', ['0', '3 - 2*t0'], ['raise', 'ValueError', 'The other_linear case is not supported in this function.'], 2],
    [['1', '(1 + t0**2)'], 'auto', ['-3*t0**2', '-2*t0 + 2*x'], ['raise', 'NotImplementedError', 'sorted roots not supported over QQ'], 2],
    [['1', '1/(x + 1)'], 'tan', ['-9*t0**3*x**2 - 15*t0**3*x - 6*t0**3 + 2*t0**2*x**4 + 10*t0**2*x**3 + 16*t0**2*x**2 + 8*t0**2*x - 4*t0*x**4 - 16*t0*x**3 - 16*t0*x**2 + 12*t0*x + 12*t0 - 8*x**3 - 40*x**2 - 76*x - 56', '3*t0**3*x**3 + 9*t0**3*x**2 + 6*t0**3*x + 2*t0**2*x**5 + 10*t0**2*x**4 + 16*t0**2*x**3 + 8*t0**2*x**2 - 6*t0*x**2 - 18*t0*x - 12*t0 - 4*x**4 - 20*x**3 - 32*x**2 - 16*x'], ['raise', 'PolynomialError', '1/(27*x**6 + 54*x**5 + 27*x**4) contains an element of the set of generators.'], 1],
    [['1', '2*t0'], None, ['-2*t0**2*x - 1', '-2*x/(x + 2)', '0', '-2'], None, 1],
    [['1', '1/(x + 1)'], None, ['3', '1', '3', '2*x/(x + 1)'], None, 1],
    [['1', 't0'], None, ['0', '-2*x', '0', '3'], None, 1],
    [['1', 't0', 't1/x'], None, ['-1', '-2*t0*t1 - 1/(x + 1)', 't1', '1'], None, 2],
    [['1', 't0', '1/x'], None, ['x', '-2*t1/(t0 + 1) + 1', 't1 - 1', '2*t1/(t0 + 1)'], None, 1],
    [['1', '1/x', 't1'], None, ['3*t0*t1**2 + 3*t0*t1 + 3*x', '-2', 't0*t1', '3*t1 - 2'], None, 1],
    [['1', '1/x', '1/(x*t0)'], None, ['-2*x', '-2*t0 - t1', '3*t1**2', '-1/(x + 2)'], None, 2],
    [['1', '1/(x + 1)'], None, ['1', '3*t0**2 - 2*t0*x', '2', '2'], None, 1],
    # Regression (random towers): SymPy residue_reduce gives ([], False)
    # through a constant S_i of splitfactor_sqf; the native kernel used to
    # decline on these.
    # rt 90 (exp(x), log(x))
    [['1', 't0', '1/x'], 'primitive', ['2*t1/(x + 2) + 1', '-2*t1**2 - 2*t1 + 3'], None, 1],
    # rt 128 (log(x), log(log(x)))
    [['1', '1/x', '1/(x*t0)'], 'auto', ['3*t1 + 2/(x + 1)', '-2*t0*t1**2 - 2*t0*t1/(t0 + 1) - x'], None, 1],
    # rt 286 (log(x+1))
    [['1', '1/(x + 1)'], None, ['-1', 't0**2*x + 2*t0 - 1', 't0**2 + 2', 't0*x + 3'], None, 2],
    # rt 403 (exp(x))
    [['1', 't0'], 'auto', ['-2*t0', '-2*t0**2 + t0/(x + 2) + x'], None, 1],
    # rt 478 (log(x+1))
    [['1', '1/(x + 1)'], None, ['1 - t0', 't0**2 - t0 + 1', 't0**2 + t0', 't0/(x + 2) - 1'], None, 2],
    # rt 586 (exp(x), log(x))
    [['1', 't0', '1/x'], 'auto', ['-t1**2 - t1/(x + 1)', '2 - t1**2'], None, 1],
    # case='base' above the base level with a residue term: the stock u is
    # the uncancelled Mul(*[Pow(h, rho*n)]).
    [['1', '1/x'], 'base', ['1', 'x*t0'], "(Integer(1), Symbol('t0'))", 3],
]


# Cases the native kernels leave to the stock body (index: reason).
DECLINES = {
    # explicit case 'exp' on a primitive level: the stock exp branch runs
    # on a level whose derivation is not exp-like.
    6: 'explicit exp case on a primitive level',
    # real_roots over QQ(...) raises in the stock body (kept content of
    # gcd(r, kappa(r)) and a non-integral monic s).
    43: 'stock real_roots NotImplementedError',
    59: 'base branch above the base level with residue terms',
}

# Case 44 (the kkinv PolynomialError of SymPy residue_reduce) is reproduced
# by the native kernel since 2026-11-06 (prde3); no recorded gap remains.
GAPS = {}


def _make(D, args):
    syms = [Symbol('x')] + [Symbol('t%d' % i) for i in range(len(D) - 1)]
    ns = {str(s): s for s in syms}
    ex = lambda s: sympify(s, locals=ns)
    DE = DifferentialExtension(extension={
        'D': [Poly(ex(d), syms[i]) for i, d in enumerate(D)], 'T': syms})
    return DE, [Poly(ex(a), syms[-1]) for a in args]


def _run(D, case, args):
    DE, P = _make(D, args)
    try:
        if case is None:
            r = parametric_log_deriv_heu(*P, DE)
        else:
            r = is_log_deriv_k_t_radical_in_field(P[0], P[1], DE, case)
        return None if r is None else srepr(r)
    except Exception as e:
        return ['raise', type(e).__name__, str(e)]


def _params():
    for i, c in enumerate(CASES):
        marks = [pytest.mark.xfail(reason=GAPS[i], strict=True)] if i in GAPS else []
        yield pytest.param(*c, marks=marks, id=str(i))


@pytest.mark.parametrize('D, case, args, expected, dummies', list(_params()))
def test_tower_against_sympy(D, case, args, expected, dummies):
    _run(D, case, args)
    d0 = Dummy().dummy_index
    assert _run(D, case, args) == expected
    assert Dummy().dummy_index - d0 - 1 == dummies


@pytest.mark.parametrize('index', range(len(CASES)))
def test_native_answers(index):
    # The native tower kernel answers (or raises the stock exception)
    # rather than declining to the stock body, except on DECLINES.
    D, case, args, expected, dummies = CASES[index]
    DE, P = _make(D, args)
    try:
        if case is None:
            r = _prde_native.parametric_log_deriv_heu(*P, DE, None)
        else:
            r = _prde_native.is_log_deriv_k_t_radical_in_field(P[0], P[1], DE, case, None)
    except Exception:
        r = ('raised',)
    assert (r is None) == (index in DECLINES)


def test_poly_free_symbols_include_domain_symbols():
    # Poly.free_symbols includes the symbols of a fraction-field ground
    # domain (stock real_roots relies on it).
    x, z = symbols('x z')
    assert Poly(z - sympify(3)/2, z, domain='QQ(x)').free_symbols == {x, z}
    assert Poly(1, z, domain='QQ(x)').free_symbols == {x}
    assert Poly(x*z, z, domain='QQ[x]').free_symbols == {x, z}
    assert Poly(z, z, domain=QQ).free_symbols == {z}
