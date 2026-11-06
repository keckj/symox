"""Native parametric Risch DE above the base level (prde3, 2026-11-06):
prde_no_cancel_b_small above the base level, the b_equal
NotImplementedError, the hypertangent and non-linear levels
(prde_special_denom tan arm with real_imag/recognize_log_derivative,
bound_degree tan/other_nonlinear arm, the non-linear NotImplementedError),
and the ValueError of an other_linear level.

Expected outputs (structural serialization of Polys: expression, gens,
domain; exception type and message; Dummies created by a repeated call) are
hard-coded from SymPy 1.14 (sympy.integrals.prde). Each case also checks
that the native tower kernel answered (STATS counter) instead of declining
to the stock body.
"""
import pytest

from symox import Symbol, Poly, Dummy, srepr, sympify
from symox.integrals import prde
from symox.integrals import _prde_native
from symox.integrals.risch import DifferentialExtension

# [function, D, args, expected, dummies, native STATS counter]
CASES = [['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '2', ['t0 + x', 'x*t0 - 1'], 2],
  [[['Poly', "Add(Mul(Rational(1, 2), Symbol('x')), Rational(-1, 4))", ["Symbol('t0')"], 'ZZ(x)'],
    ['Poly', 'Rational(-1, 2)', ["Symbol('t0')"], 'QQ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
   '[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(-1), '
   'Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(-1), '
   'Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(-1)]])'],
  1, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '0', ['t0**2 + 1', '2*t0'], 2],
  [[['Poly', 'Integer(1)', ["Symbol('t0')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
    ['Poly', "Symbol('t0')", ["Symbol('t0')"], 'ZZ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0)], [Integer(-1), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), '
   'Integer(0)], [Integer(0), Integer(-1), Integer(0), Integer(0), Integer(1), Integer(0), '
   'Integer(0)], [Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(-1), '
   'Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(-1)]])'],
  1, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '2*t0', ['t0 + x'], 1],
  ['raise', 'NotImplementedError', 'prde_no_cancel_b_equal() is not yet implemented.'], 0,
  'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '-3*t0', ['t0 + x'], 1],
  ['raise', 'NotImplementedError',
   'non-linear and hypertangent cases have not yet been implemented'],
  0, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', 't0**2 + 1', ['t0**3', 'x*t0'], 1],
  [[['Poly', "Add(Symbol('t0'), Integer(-1))", ["Symbol('t0')"], 'ZZ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), '
   'Integer(0), Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(1), '
   'Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(-1)]])'],
  0, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['2', '-x', ['x*t0', '1'], 1],
  [[['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
   '[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(-1), '
   'Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), '
   'Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(-1)]])'],
  1, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', 't0**3 + 1'], ['1', 't0', ['t0**2 + x', 't0'], 2],
  [[['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'], ['Poly', 'Integer(1)', ["Symbol('t0')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), '
   'Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(1), '
   'Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(-1)]])'],
  0, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', 't0**2'], ['1', '1', ['t0', 'x'], 1],
  [[['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
    ['Poly', "Add(Symbol('x'), Integer(-1))", ["Symbol('t0')"], 'ZZ(x)'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'],
    ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
   '[Integer(-1), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), '
   'Integer(-1), Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(0), '
   'Integer(0), Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), '
   'Integer(0), Integer(0), Integer(-1)]])'],
  1, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', 't0 + 1'], ['1', '1', ['t0', 'x'], 1],
  ['raise', 'NotImplementedError',
   'non-linear and hypertangent cases have not yet been implemented'],
  0, 'param_poly_rischDE_tower'],
 ['param_poly_rischDE', ['1', '1/x', '1 + t1**2'], ['1', '1', ['t1 + t0', 'x'], 1],
  [[['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'],
    ['Poly', "Symbol('t0')", ["Symbol('t1')"], 'ZZ(t0)'],
    ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'],
    ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'],
    ['Poly', "Add(Symbol('x'), Integer(-1))", ["Symbol('t1')"], 'ZZ(x)'],
    ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'],
    ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'],
    ['Poly', 'Integer(0)', ["Symbol('t1')"], 'QQ'], ['Poly', 'Integer(0)', ["Symbol('t1')"], 'ZZ'],
    ['Poly', 'Integer(0)', ["Symbol('t1')"], 'ZZ']],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
   '[Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(-1), Integer(0), Integer(1), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0)], [Integer(0), Integer(-1), Integer(0), Integer(0), Integer(0), '
   'Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0)], [Integer(0), Integer(-1), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
   '[Integer(0), Integer(-1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
   'Integer(-1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0)], [Integer(0), Integer(-1), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), '
   'Integer(0)], [Integer(0), Integer(-1), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(-1), Integer(0)], '
   '[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(0), Integer(-1)]])'],
  6, 'param_poly_rischDE_tower'],
 ['param_rischDE', ['1', '1 + t0**2'], ['1', 't0', [['1', 't0'], ['t0', '1']]],
  ['raise', 'NotImplementedError',
   'non-linear and hypertangent cases have not yet been implemented'],
  5, 'param_rischDE_tower'],
 ['param_rischDE', ['1', '1 + t0**2'], ['2*t0', 't0**2 + 1', [['1', 't0**2 + 1']]],
  [[[['Poly', "Mul(Integer(-1), Symbol('t0'))", ["Symbol('t0')"], 'QQ'],
     ['Poly', 'Integer(2)', ["Symbol('t0')"], 'QQ']]],
   'MutableDenseMatrix([[Integer(1), Integer(0)], [Integer(0), Integer(1)]])'],
  1, 'param_rischDE_tower'],
 ['param_rischDE', ['1', '1 + t0**2'], ['x', '1', [['t0', '1'], ['1', '1']]],
  [[[['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
     ['Poly', 'Integer(1)', ["Symbol('t0')"], 'QQ']],
    [['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
     ['Poly', 'Integer(1)', ["Symbol('t0')"], 'QQ']],
    [['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'],
     ['Poly', 'Integer(1)', ["Symbol('t0')"], 'ZZ']],
    [['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'],
     ['Poly', 'Integer(1)', ["Symbol('t0')"], 'ZZ']]],
   'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
   '[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
   'Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), '
   'Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(1), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
   'Integer(0), Integer(1)]])'],
  5, 'param_rischDE_tower'],
 ['param_rischDE', ['1', 't0**2'], ['1', 't0', [['1', 't0']]],
  [[[['Poly', 'Integer(1)', ["Symbol('t0')"], 'QQ'],
     ['Poly', 'Integer(1)', ["Symbol('t0')"], 'QQ']]],
   'MutableDenseMatrix([[Integer(-1), Integer(1)]])'],
  1, 'param_rischDE_tower'],
 ['param_rischDE', ['1', 't0 + 1'], ['1', '1', [['t0', '1']]],
  ['raise', 'ValueError',
   "case must be one of {'exp', 'tan', 'primitive', 'base'}, not other_linear."],
  1, 'param_rischDE_tower'],
 ['limited_integrate', ['1', '1 + t0**2'], ['2*t0', 't0**2 + 1', [['1', '1']]],
  ['raise', 'NotImplementedError',
   'non-linear and hypertangent cases have not yet been implemented'],
  1, 'limited_integrate_tower'],
 ['prde_no_cancel_b_small', ['1', '1 + t0**2'], ['3', ['t0**2 + x', 'x*t0'], 2],
  [[[['Poly', "Add(Mul(Rational(1, 3), Symbol('x')), Rational(-4, 9))", ["Symbol('t0')"], 'ZZ(x)'],
     ['Poly', 'Integer(0)', ["Symbol('t0')"], 'QQ'],
     ['Poly', "Symbol('t0')", ["Symbol('t0')"], 'ZZ'],
     ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
    'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
    '[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(-1), '
    'Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(-1), '
    'Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), '
    'Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0), '
    'Integer(0), Integer(-1)]])'],
   [['Poly', "Add(Mul(Integer(-1), Integer(3), Symbol('t0')), Symbol('x'), Integer(-1))",
     ["Symbol('t0')"], 'ZZ[x]'],
    ['Poly', "Mul(Symbol('t0'), Symbol('x'))", ["Symbol('t0')"], 'ZZ[x]']]],
  1, 'prde_no_cancel_b_small_tower'],
 ['prde_no_cancel_b_small', ['1', 't0**3 + 1'], ['t0', ['t0**2 + x', '2'], 1],
  [[[['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ'],
     ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
    'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), '
    'Integer(1), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0)], '
    '[Integer(1), Integer(0), Integer(-1), Integer(0)], [Integer(0), Integer(1), Integer(0), '
    'Integer(-1)]])'],
   [['Poly', "Add(Pow(Symbol('t0'), Integer(2)), Symbol('x'))", ["Symbol('t0')"], 'ZZ[x]'],
    ['Poly', 'Integer(2)', ["Symbol('t0')"], 'ZZ']]],
  0, 'prde_no_cancel_b_small_tower'],
 ['prde_no_cancel_b_large', ['1', '1 + t0**2'], ['t0**2', ['t0**3 + x', 't0'], 1],
  [[[['Poly', "Add(Symbol('t0'), Integer(-1))", ["Symbol('t0')"], 'ZZ'],
     ['Poly', 'Integer(0)', ["Symbol('t0')"], 'ZZ']],
    'MutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0)], [Integer(0), '
    'Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(-1), Integer(0)], '
    '[Integer(0), Integer(1), Integer(0), Integer(-1)]])'],
   [['Poly', "Add(Symbol('x'), Integer(-1))", ["Symbol('t0')"], 'ZZ[x]'],
    ['Poly', "Symbol('t0')", ["Symbol('t0')"], 'ZZ']]],
  0, 'prde_no_cancel_b_large_tower']]

# Random nonlinear-corpus calls 190 and 454: a param_rischDE (nested in
# prde_no_cancel_b_small for 190) with no polynomial solution, v == 0 and W
# empty: SymPy builds Matrix(0, m, [], t) from `wl[:m] + wl[-v:]` and its
# nullspace is the identity.
CASES += [['param_poly_rischDE',
  ['1', '2 + 2*t0**2'],
  ['1', '(-1/(x + 1))*t0**0', ['(1)*t0**0 + (3/(x + 2))*t0**1 + (1*x)*t0**2'], 1],
  [[['Poly', "Mul(Rational(1, 2), Symbol('t0'), Symbol('x'))", ["Symbol('t0')"], 'QQ[x]']],
   'MutableDenseMatrix([[Integer(1), Integer(0)], [Integer(1), Integer(0)], [Integer(1), '
   'Integer(-1)]])'],
  2,
  'param_poly_rischDE_tower'],
 ['param_rischDE',
  ['1', 'x*(1 + t0**2)'],
  ['(1)*t0**0', '(-2*x/(x + 0))*t0**0', [['(0)*t0**0 + (1)*t0**1 + (-2)*t0**2', 'x']]],
  [[[['Poly', "Mul(Integer(-1), Integer(2), Symbol('t0'))", ["Symbol('t0')"], 'ZZ(x)'],
     ['Poly', "Pow(Symbol('x'), Integer(2))", ["Symbol('t0')"], 'ZZ(x)']]],
   'MutableDenseMatrix([[Integer(1), Integer(0)], [Integer(0), Integer(1)]])'],
  10,
  'param_rischDE_tower']]


def _ser(o):
    if isinstance(o, Poly):
        return ['Poly', srepr(o.as_expr()), [srepr(g) for g in o.gens], str(o.domain)]
    if isinstance(o, (list, tuple)):
        return [_ser(x) for x in o]
    if hasattr(o, 'to_Matrix'):
        o = o.to_Matrix()
    return srepr(o)


def _build(kind, D, args):
    T = [Symbol('x')] + [Symbol('t%d' % i) for i in range(len(D) - 1)]
    ns = {str(s): s for s in T}
    ex = lambda s: sympify(s, locals=ns)
    DE = DifferentialExtension(extension={
        'D': [Poly(ex(d), T[i]) for i, d in enumerate(D)], 'T': T})
    P = lambda s: Poly(ex(s), T[-1])
    if kind == 'param_poly_rischDE':
        a, b, q, n = args
        return (P(a), P(b), [P(x) for x in q], n, DE), 2
    if kind in ('param_rischDE', 'limited_integrate'):
        fa, fd, G = args
        return (P(fa), P(fd), [(P(x), P(y)) for x, y in G], DE), 2
    b, q, n = args
    return (P(b), [P(x) for x in q], n, DE), 1


def _run(kind, D, args):
    a, qi = _build(kind, D, args)
    try:
        out = _ser(getattr(prde, kind)(*a))
        if kind.startswith('prde_no_cancel'):
            # Q is updated in place, as in the stock body.
            out = [out, _ser(a[qi])]
        return out
    except Exception as e:
        return ['raise', type(e).__name__, str(e)]


@pytest.mark.parametrize('kind, D, args, expected, dummies, stat', CASES,
                         ids=[str(i) for i in range(len(CASES))])
def test_against_sympy_and_native(kind, D, args, expected, dummies, stat):
    _run(kind, D, args)
    before = dict(_prde_native.STATS)
    d0 = Dummy().dummy_index
    assert _run(kind, D, args) == expected
    assert Dummy().dummy_index - d0 - 1 == dummies
    # The native tower kernel answered (or raised the stock exception).
    assert _prde_native.STATS[stat] == before[stat] + 1


def test_stock_fallback_matches(monkeypatch):
    # With the native kernels disabled the stock body gives the same results.
    monkeypatch.setattr(_prde_native, '_tower', lambda *a, **k: None)
    for kind, D, args, expected, dummies, stat in CASES:
        assert _run(kind, D, args) == expected


def test_polymatrix_generator_order():
    # Random nonlinear-corpus call 4 (tan level over log): SymPy 1.14 infers
    # the PolyMatrix domain with `_sort_gens` (x before t0), QQ(x,t0); the
    # native kernel used to give the alphabetical QQ(t0,x).
    x, t0, t1 = Symbol('x'), Symbol('t0'), Symbol('t1')
    DE = DifferentialExtension(extension={
        'D': [Poly(1, x), Poly(1/x, t0), Poly(1 + t1**2, t1)], 'T': [x, t0, t1]})
    P = lambda e: Poly(e, t1)
    before = _prde_native.STATS['param_poly_rischDE_tower']
    h, M = prde.param_poly_rischDE(P(1), P(-1/(t0 + 1) - 2*t1 + t1**2),
        [P(2), P(1 - 2*t1), P(2*x/(t0 + 2) + 3/(x + 2)*t1)], 1, DE)
    assert _prde_native.STATS['param_poly_rischDE_tower'] == before + 1
    assert [p.is_zero for p in h] == [True, True, True]
    assert str(M.domain) == 'QQ(x,t0)' and M.gens == (t1,) and M.shape == (6, 6)


@pytest.mark.parametrize('values, ring, field', [
    ('[a*x]', 'ZZ[x,a]', 'ZZ(x,a)'),
    ('[a, 1/x]', 'ZZ(x,a)', 'ZZ(x,a)'),
    ('[x/(t0 + 2), 1]', 'ZZ(x,t0)', 'ZZ(x,t0)'),
    ('[t0, x]', 'ZZ[x,t0]', 'ZZ(x,t0)'),
    ('[b*y, a]', 'ZZ[y,a,b]', 'ZZ(y,a,b)'),
])
def test_construct_domain_generator_order(values, ring, field):
    # SymPy 1.14 outputs: the generators of a symbol-only composite domain
    # follow `_sort_gens` (the order of the native `_composite_generators`).
    from symox.polys.constructor import construct_domain
    v = sympify(values)
    assert str(construct_domain(v)[0]) == ring
    assert str(construct_domain(v, field=True)[0]) == field
