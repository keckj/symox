"""Native parametric Risch DE above the base level
(crates/symox-py/src/prde_param.rs): param_rischDE, param_poly_rischDE,
prde_cancel_liouvillian and limited_integrate on exp/primitive towers.

The expected outputs (srepr of the Poly expressions with their generators
and domains, PolyMatrix entries and domain, exception type and message;
a sha256 of that encoding for the long ones) and the number of Dummies
created by the call are hard-coded from SymPy 1.14 (sympy.integrals.prde).
"""
import hashlib

import pytest

from symox import Poly, Dummy, srepr, sympify, symbols
from symox.integrals import prde
from symox.integrals.risch import DifferentialExtension as DEx
from symox.integrals import _prde_native

x, t, t0, t1 = symbols('x t t0 t1')


def enc(v):
    if isinstance(v, Poly):
        return ('P', srepr(v.as_expr()), tuple(srepr(g) for g in v.gens), str(v.domain))
    if isinstance(v, (list, tuple)):
        return tuple(enc(i) for i in v)
    if hasattr(v, 'rows') and hasattr(v, 'domain'):
        return ('M', v.rows, v.cols, tuple(srepr(v[i, j].as_expr())
                for i in range(v.rows) for j in range(v.cols)), str(v.domain))
    if v is None:
        return None
    return srepr(v)


def P(e, g, **kw):
    return Poly(sympify(e), g, **kw)


def D1():
    return DEx(extension={'D': [Poly(1, x), Poly(1/x, t)]})


def Dexp():
    return DEx(extension={'D': [Poly(1, x), Poly(-t, t)]})


def Dt():
    return DEx(extension={'D': [Poly(1, x), Poly(t, t)]})


def D3():
    return DEx(extension={'D': [Poly(1, x), Poly(1/x, t0), Poly(1/(x*t0), t1)]})


def G2():
    return [(P('1', t), P('1', t)), (P('t', t), P('x', t))]


CALLS = {
    'cancel_prim': lambda: prde.prde_cancel_liouvillian(P('-1/(x - 1)', t), [P('-x + 1', t), P('1', t)], 1, D1()),
    'cancel_exp': lambda: prde.prde_cancel_liouvillian(Poly(0, t, domain='QQ[x]'), [Poly(1, t, domain='QQ(x)')], 0, Dexp()),
    'risch_exp': lambda: prde.param_rischDE(P('-t - 1', t, field=True), P('t + x', t, field=True),
        [(Poly(t + x, t, domain='ZZ(x)'), Poly(1, t, domain='QQ')), (Poly(0, t, domain='QQ'), Poly(1, t, domain='QQ'))], Dt()),
    'risch_prim': lambda: prde.param_rischDE(P('1', t), P('x', t), [(P('1', t), P('x*t', t)), (P('t', t), P('1', t))], D1()),
    'poly_prim': lambda: prde.param_poly_rischDE(P('1', t), P('0', t), [P('t', t), P('1/x', t)], 2, D1()),
    'poly_exp_large': lambda: prde.param_poly_rischDE(P('1', t), P('t', t), [P('t**2', t), P('t', t)], 2, Dt()),
    'poly_spde': lambda: prde.param_poly_rischDE(P('t', t), P('1', t), [P('t**2 + x', t), P('t', t)], 2, Dt()),
    'poly_neg': lambda: prde.param_poly_rischDE(P('t', t), P('1', t), [P('t**2 + x', t), P('t', t)], -1, Dt()),
    'li_prim': lambda: prde.limited_integrate(P('1', t), P('x', t), [(P('1', t), P('x', t))], D1()),
    'li_prim2': lambda: prde.limited_integrate(P('t', t), P('x', t), [(P('1', t), P('x', t))], D1()),
    'li_exp': lambda: prde.limited_integrate(P('t', t), P('1', t), [(P('1', t), P('x', t))], Dt()),
    'li_none': lambda: prde.limited_integrate(P('1', t), P('t', t), [(P('t', t), P('1', t))], Dt()),
    'li_tower3': lambda: prde.limited_integrate(P('1', t1), P('x*t0', t1), [(P('1', t1), P('x', t1))], D3()),
    'risch_tower3': lambda: prde.param_rischDE(P('0', t1), P('1', t1), [(P('1', t1), P('x*t0', t1)), (P('t1', t1), P('x', t1))], D3()),
    # weak_normalizer: real_roots of the resultant over ZZ(x) (cleared
    # denominators, primitive part, retracted to ZZ), symbolic roots.
    'wn_root1': lambda: prde.param_rischDE(P('t + 1', t), P('t + x', t), G2(), Dt()),
    'wn_symroot': lambda: prde.param_rischDE(P('x*t', t), P('t + x', t), G2(), Dt()),
    'wn_strip': lambda: prde.param_rischDE(P('x*t + 1', t), P('x*t', t), G2(), Dt()),
    'wn_prim_root': lambda: prde.param_rischDE(P('2/x', t), P('t + x', t), G2(), D1()),
    'li_exp_den': lambda: prde.limited_integrate(P('1', t), P('t + 1', t), [(P('t', t), P('t + 1', t))], Dt()),
    'li_prim_den': lambda: prde.limited_integrate(P('1', t), P('x*t + x', t), [(P('1', t), P('x*t + x', t))], D1()),
    'li_prim_base_den': lambda: prde.limited_integrate(P('t', t), P('x - 1', t), [(P('1', t), P('x - 1', t))], D1()),
    'cancel_exp_den': lambda: prde.prde_cancel_liouvillian(Poly(0, t, domain='QQ(x)'),
        [Poly(1/(x + 1), t, domain='QQ(x)'), Poly(t/(x - 1), t, domain='QQ(x)')], 1, Dexp()),
    # prde_no_cancel_b_small, base case of the coefficients (b == 0).
    'poly_prim_bsmall': lambda: prde.param_poly_rischDE(P('1', t), P('0', t),
        [P('t/(x+1)', t), P('1/(x**2-1)', t), P('t**2', t)], 3, D1()),
    # Primitive bound_degree with deg(b) == deg(a) and alpha = Dz/z in k:
    # beta, derivation(z) and limited_integrate one level down.
    'bd_prim_eq': lambda: prde.param_rischDE(P('-1', t), P('x', t), [(P('1', t), P('1', t)), (P('t', t), P('1', t))], D1()),
    'bd_prim_eq2': lambda: prde.param_rischDE(P('-2', t), P('x', t), [(P('t', t), P('1', t)), (P('x*t**2', t), P('1', t))], D1()),
    'bd_prim_eq3': lambda: prde.param_rischDE(P('-x - 1', t), P('x', t), [(P('t', t), P('1', t)), (P('1', t), P('x', t))], D1()),
    'li_prim_eq': lambda: prde.limited_integrate(P('t', t), P('1', t), [(P('1', t), P('x', t))], D1()),
}

# name: (repr of the encoded SymPy 1.14 result, or its sha256; Dummies)
EXPECTED = {
    'cancel_prim': ('sha256:35c38cb4dacdad2f1c6d0ad1fed1630b0fd643af0985ae5908aa3a828c92ce58', 4),
    'cancel_exp': ('(((\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\'), (\'P\', "Symbol(\'x\')", ("Symbol(\'t\')",), \'ZZ(x)\')), (\'M\', 1, 3, (\'Integer(-1)\', \'Integer(0)\', \'Integer(1)\'), \'QQ\'))', 1),
    'risch_exp': ('sha256:077323056b94890d52797fe206b1a265036d4a461639a9abcc29f2ac96173d69', 5),
    'risch_prim': ("('raise', 'SympifyError', 'SympifyError: None')", 4),
    'poly_prim': ('sha256:035229f5ac31b15afa705f6aefd10c17bffbb8153b24533590366ed7be3a8f80', 3),
    'poly_exp_large': ('(((\'P\', "Add(Symbol(\'t\'), Integer(-1))", ("Symbol(\'t\')",), \'ZZ\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'ZZ\')), (\'M\', 2, 4, (\'Integer(1)\', \'Integer(0)\', \'Integer(-1)\', \'Integer(0)\', \'Integer(0)\', \'Integer(1)\', \'Integer(0)\', \'Integer(-1)\'), \'QQ\'))', 0),
    'poly_spde': ('sha256:951631df5edc445dc5802d2c9bb6af1f76b5e3aee090c7cd94efcf4ef55192f5', 0),
    'poly_neg': ("((), ('M', 3, 2, ('Integer(1)', 'Integer(0)', 'Integer(0)', 'Integer(1)', 'Integer(0)', 'Integer(0)'), 'QQ(x)'))", 0),
    'li_prim': ('(((\'P\', "Symbol(\'t\')", ("Symbol(\'t\')",), \'ZZ\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(0)\', ("Symbol(\'t\')",), \'QQ\'),))', 3),
    'li_prim2': ('(((\'P\', "Mul(Rational(1, 2), Pow(Symbol(\'t\'), Integer(2)))", ("Symbol(\'t\')",), \'QQ\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(0)\', ("Symbol(\'t\')",), \'QQ\'),))', 4),
    'li_exp': ('(((\'P\', "Mul(Rational(1, 2), Symbol(\'t\'), Symbol(\'x\'))", ("Symbol(\'t\')",), \'QQ[x]\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(0)\', ("Symbol(\'t\')",), \'QQ\'),))', 6),
    'li_none': ('(((\'P\', \'Integer(-1)\', ("Symbol(\'t\')",), \'ZZ\'), (\'P\', "Symbol(\'t\')", ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(0)\', ("Symbol(\'t\')",), \'QQ\'),))', 7),
    'li_tower3': ('(((\'P\', "Symbol(\'t1\')", ("Symbol(\'t1\')",), \'ZZ\'), (\'P\', \'Integer(1)\', ("Symbol(\'t1\')",), \'QQ\')), ((\'P\', \'Integer(0)\', ("Symbol(\'t1\')",), \'QQ\'),))', 6),
    'risch_tower3': ('sha256:96d30931120552d7c83c501cd1ff264275d120c18eee719af23e8585c23dcb4c', 37),
    'wn_root1': ('sha256:b208a461eb403ba452ad07a3773c08ce3835b528670378ce190098b4c0e75298', 9),
    'wn_symroot': ("('raise', 'NotImplementedError', 'sorted roots not supported over ZZ[x]')", 2),
    'wn_strip': ('sha256:7104d536a8ef464b2cbcc5614f5a8b921ed39b84f988f11acb28ec9e6dc66938', 2),
    'wn_prim_root': ("('raise', 'NotImplementedError', 'sorted roots not supported over ZZ[x]')", 2),
    'li_exp_den': ('(((\'P\', "Symbol(\'x\')", ("Symbol(\'t\')",), \'ZZ[x]\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(-1)\', ("Symbol(\'t\')",), \'QQ\'),))', 2),
    'li_prim_den': ('(((\'P\', \'Integer(0)\', ("Symbol(\'t\')",), \'ZZ\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ(x)\'),))', 2),
    'li_prim_base_den': ('None', 4),
    'cancel_exp_den': ('sha256:d39569d4a83549e0bac22e194c6cbb3cea63be3533caef1157a43bb344ebc028', 3),
    'poly_prim_bsmall': ('sha256:e550948f91ab3ce29f4d3bbd8026a2e9b32839e3d064c575eae59087bd74a888', 4),
    'bd_prim_eq': ('sha256:e5cd042ceffcadbd4b7aba3c628e378806b7707b74e96f32337ab16471175c3f', 11),
    'bd_prim_eq2': ('sha256:1602119c2ca9bdc68f53eec595821168f256315c1a4f07e150c2d2a81477861e', 13),
    'bd_prim_eq3': ('sha256:35f657f7d1089f6219846293579ed167f32b187e9b1ecb398f1a966f9da6053d', 10),
    'li_prim_eq': ('(((\'P\', "Add(Mul(Symbol(\'t\'), Symbol(\'x\')), Mul(Integer(-1), Symbol(\'x\')))", ("Symbol(\'t\')",), \'ZZ[x]\'), (\'P\', \'Integer(1)\', ("Symbol(\'t\')",), \'QQ\')), ((\'P\', \'Integer(0)\', ("Symbol(\'t\')",), \'QQ\'),))', 11),
}

# The stock derivation(z) adds None (z rational in x): Poly + None raises
# SympifyError as in SymPy (strict sympify of `_polifyit`, since prde3);
# the native kernel declines there and the stock body raises it.
GAPS = {}

# Cases left to the stock body by the native kernel: real_roots whose
# preprocessed polynomial is not over ZZ (the stock NotImplementedError).
DECLINES = {'risch_prim', 'wn_symroot', 'wn_prim_root'}

STAT = {'cancel': 'prde_cancel_liouvillian_tower', 'risch': 'param_rischDE_tower',
        'poly': 'param_poly_rischDE_tower', 'wn': 'param_rischDE_tower',
        'bd': 'param_rischDE_tower', 'li': 'limited_integrate_tower'}


def _run(name):
    try:
        r = enc(CALLS[name]())
    except Exception as e:
        r = ('raise', type(e).__name__, str(e))
    s = repr(r)
    return s if len(s) <= 300 else 'sha256:' + hashlib.sha256(s.encode()).hexdigest()


def _params():
    for name in CALLS:
        marks = [pytest.mark.xfail(reason=GAPS[name], strict=True)] if name in GAPS else []
        yield pytest.param(name, marks=marks, id=name)


@pytest.mark.parametrize('name', list(_params()))
def test_param_tower_against_sympy(name):
    expected, dummies = EXPECTED[name]
    d0 = Dummy().dummy_index
    assert _run(name) == expected
    assert Dummy().dummy_index - d0 - 1 == dummies


@pytest.mark.parametrize('name', list(CALLS))
def test_native_answers(name):
    # The native kernel answers (or raises the stock exception) rather than
    # declining to the stock body, except on DECLINES.
    stat = STAT[name.split('_')[0]]
    before = _prde_native.STATS[stat]
    _run(name)
    assert (_prde_native.STATS[stat] > before) == (name not in DECLINES)
