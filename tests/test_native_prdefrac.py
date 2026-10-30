"""Wave 30 (prdefrac): PolyMatrix exact-rational entries and native admission
helpers on the parametric Risch DE path.

Public results are hard-coded from SymPy 1.14; the fast paths are compared
with the stock ``Poly`` construction they replace (same Poly, domain, gens),
and every changed protocol must decline to the stock code.
"""
import os
import subprocess
import sys

import pytest

import symox as s
from symox import _core, _poly_constructor as constructor
from symox.integrals.risch import NonElementaryIntegral, risch_integrate
from symox.polys import _polymatrix_rational as rational
from symox.polys import polymatrix as pm

x, y, z = s.symbols('x y z')
t = s.Dummy('t')
LOG = s.log
QQ, ZZ = s.QQ, s.ZZ


# SymPy 1.14: risch_integrate(f, x) for each f, as (str, top-level classes).
EXPECTED = [
    (LOG(LOG(x)), 'x*log(log(x)) + Integral(-1/log(x), x)', ['Mul', 'NonElementaryIntegral']),
    (LOG(x)*LOG(x + 1), '(x*log(x) - x)*log(x + 1) + Integral((-x*log(x) + x)/(x + 1), x)',
     ['Mul', 'NonElementaryIntegral']),
    (1/(x*LOG(x)), 'log(log(x))', ['log']),
    (s.exp(x)/(1 + s.exp(x))**2, '-1/(exp(x) + 1)', ['Mul']),
    (x*s.exp(x)*LOG(x), 'Integral(x*exp(x)*log(x), x)', ['NonElementaryIntegral']),
]


@pytest.mark.parametrize('f,text,kinds', EXPECTED)
def test_risch_integrate_matches_sympy(f, text, kinds):
    result = risch_integrate(f, x)
    assert str(result) == text
    assert [type(a).__name__ for a in s.Add.make_args(result)] == kinds


def test_loglog_parts_are_exact():
    result = risch_integrate(LOG(LOG(x)), x)
    elementary, integral = s.Add.make_args(result)
    assert elementary == x*LOG(LOG(x))
    assert isinstance(integral, NonElementaryIntegral) and integral.args == (-1/LOG(x), (x,))


DOMAINS = [QQ, QQ.frac_field(y), ZZ.frac_field(y), QQ.frac_field(y, z), ZZ.frac_field(y, z)]
ENTRIES = [s.S.Zero, s.S.One, s.S.NegativeOne, s.S.Half, s.Integer(7), s.Rational(-5, 3)]


@pytest.mark.parametrize('domain', DOMAINS)
@pytest.mark.parametrize('gen', [x, t])
def test_ground_poly_is_the_stock_poly(domain, gen):
    m = pm.PolyMatrix([[e] for e in ENTRIES], gen, ring=domain)
    for e in ENTRIES:
        fast = rational.ground_poly(m, e)
        stock = s.Poly(e, gen, domain=domain)
        assert fast is not None
        assert fast == stock and fast.domain == stock.domain and fast.gens == stock.gens
        assert fast.as_expr() == e and fast.is_ground
    assert rational.ground_items(m)


@pytest.mark.parametrize('domain', DOMAINS[1:])
def test_fraction_field_rref_nullspace_equal_stock(monkeypatch, domain):
    rows = [[1, s.Rational(1, 2), 0], [2, 1, 3]]
    m = pm.PolyMatrix(rows, t, ring=domain)
    fast = (m.rref(), m.nullspace())
    monkeypatch.setattr(rational, 'ground_items', lambda matrix: False)
    monkeypatch.setattr(rational, 'ground_poly', lambda matrix, e: None)
    stock = (m.rref(), m.nullspace())
    assert fast == stock
    assert all(v.domain == domain and v.gens == (t,) for v in fast[1])


def test_declines_outside_the_admitted_domains():
    cases = [
        (pm.PolyMatrix([[1]], x, ring=QQ.frac_field(x, y)), s.S.One),  # shared generator
        (pm.PolyMatrix([[1]], t, ring=s.RR), s.S.One),
        (pm.PolyMatrix([[1]], t, ring=ZZ), s.S.One),  # not a field
        (pm.PolyMatrix([[1]], t, ring=QQ.algebraic_field(s.sqrt(2)).frac_field(y)), s.S.One),
        (pm.PolyMatrix([[1]], t, ring=QQ), x),  # not an exact rational
        (pm.PolyMatrix([[1]], t, ring=QQ), s.Float(1)),
        (pm.PolyMatrix([[1]], (t, x), ring=QQ), s.S.One),  # two generators
    ]
    for m, e in cases:
        assert rational.ground_poly(m, e) is None
    assert not rational.ground_items(pm.PolyMatrix([[x]], t, ring=QQ.frac_field(y)))
    with pytest.raises(Exception, match='disjoint'):
        pm.PolyMatrix([[1]], x, ring=QQ.frac_field(x, y))[0, 0]


def test_shared_generator_keeps_stock_error():
    m = pm.PolyMatrix([[1, 2]], x, ring=QQ.frac_field(x))
    with pytest.raises(Exception) as fast:
        m.rref()
    m2 = pm.PolyMatrix([[1, 2]], x, ring=QQ.frac_field(x))
    with pytest.raises(Exception) as stock:
        s.Poly(1, x, domain=QQ.frac_field(x))
    assert type(fast.value) is type(stock.value)
    assert m2 is not m


def _matrix():
    return pm.PolyMatrix([[1, 2]], t, ring=QQ.frac_field(y))


@pytest.mark.parametrize('mode', ['construct', 'options', 'code', 'defaults', 'poly', 'sympify',
                                  'method', 'module'])
def test_changed_protocol_declines(monkeypatch, mode):
    m = _matrix()
    one = s.Poly(1, t, domain=QQ)
    assert rational.ground_poly(m, s.S.One) is not None and rational.common_ring([one]) is not None
    if mode == 'construct':
        monkeypatch.setattr(constructor, 'construct', lambda *a, **k: None)
    elif mode == 'options':
        monkeypatch.setattr(constructor, '_options', lambda *a, **k: None)
    elif mode == 'code':
        def changed(gens, domain, modulus, extension, kwargs):
            return gens, {}
        monkeypatch.setattr(constructor._options, '__code__', changed.__code__)
    elif mode == 'defaults':
        monkeypatch.setattr(constructor.construct, '__defaults__', (None, None, None, {}))
    elif mode == 'poly':
        monkeypatch.setattr(pm, 'Poly', lambda *a, **k: None)
    elif mode == 'sympify':
        monkeypatch.setattr(s, 'sympify', lambda *a, **k: None)
    elif mode == 'method':
        monkeypatch.setattr(pm.MutablePolyDenseMatrix, '_poly', lambda self, e: None)
    else:
        monkeypatch.setitem(sys.modules, 'symox._poly_constructor', object())
    assert rational.ground_poly(m, s.S.One) is None
    assert not rational.ground_items(m)
    assert rational.common_ring([one]) is None


def test_common_ring_is_the_unify_chain():
    p, q = s.Poly(t + 1, t, domain=QQ), s.Poly(2*t, t, domain=QQ)
    assert rational.common_ring([p, q]) == (QQ, (t,))
    r = p
    for other in (q,):
        r = r.unify(other)[0]
    assert r.domain == QQ
    assert rational.common_ring([s.Poly(t, t, domain=ZZ), s.Poly(t, t, domain=QQ)]) is None
    assert rational.common_ring([s.Poly(t, t, domain=QQ), s.Poly(x, x, domain=QQ)]) is None
    assert rational.common_ring([]) is None
    m = pm.PolyMatrix([[p], [q]])
    assert m.domain == QQ and m.gens == (t,) and list(m) == [p, q]


class Key(str):
    pass


def test_str_keyed():
    check = _core._str_keyed
    assert check({}) and check({'a': 1, 'b': 2})
    assert not check({1: 2}) and not check({'a': 1, Key('b'): 2})
    assert not check(type('D', (dict,), {})(a=1)) and not check([('a', 1)]) and not check(None)
    assert check(_core.__dict__)


def test_function_states_unchanged():
    check = _core._function_states_unchanged
    if check(()) is None:
        pytest.skip('public CPython function getters unavailable')

    def outer():
        value = object()

        def inner(a, b=1, *, c=2):
            return value, a, b, c
        return inner, value
    f, value = outer()
    record = (f, f.__code__, f.__defaults__, f.__kwdefaults__, (value,))
    assert check((record,)) is True
    assert check((record[:4] + ((object(),),),)) is False
    assert check(((len, None, None, None, ()),)) is False
    assert check([record]) is False and check((record[:4],)) is False
    f.__defaults__ = (2,)
    assert check((record,)) is False
    f.__defaults__ = record[2]
    assert check((record,)) is True
    f.__code__ = (lambda a, b=1, *, c=2: value).__code__
    assert check((record,)) is False


def test_no_audit_event_on_the_fast_paths():
    source = '''
import sys
import symox as s
from symox.polys import polymatrix as p
from symox.integrals.risch import risch_integrate
x,y=s.symbols('x y')
m=p.PolyMatrix([[1,2,3],[3,7,8]],x,ring=s.QQ.frac_field(y))
events=[]
def hook(event,args):
    if event=='object.__getattr__':events.append(args[1])
sys.addaudithook(hook)
from symox.polys import _polymatrix_rational as r
assert r.ground_items(m) and r.ground_poly(m,s.S.One) is not None
assert events==[],events
'''
    subprocess.run([sys.executable, '-c', source], check=True, env=os.environ.copy())


def test_constant_system_domain_is_sympy():
    # SymPy 1.14 test_prde.py::test_constant_system: the derivation step
    # builds ``D(e)/DAij`` = ``2*(x**2/2 + x + 1/2)/(x**2 + 2*x + 1)``, whose
    # ``expand`` is 1, so that row is over QQ and the result over QQ(x)[t].
    from symox.integrals.prde import constant_system
    from symox.integrals.risch import DifferentialExtension
    T = s.Symbol('t')
    A = pm.PolyMatrix([[-(x + 3)/(x - 1), (x + 1)/(x - 1), 1],
                       [-x - 3, x + 1, x - 1],
                       [2*(x + 3)/(x - 1), 0, 0]], T)
    u = pm.PolyMatrix([[(x + 1)/(x - 1)], [x + 1], [0]], T)
    DE = DifferentialExtension(extension={'D': [s.Poly(1, x)]})
    R = QQ.frac_field(x)[T]
    B, v = constant_system(A, u, DE)
    assert (B, v) == (pm.PolyMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 1]], ring=R),
                      pm.PolyMatrix([0, 1, 0, 0], ring=R))
    assert B.ring == R and v.ring == R


# SymPy 1.14 ``expand``/``expand_mul`` of products with one denominator
# factor: ``Mul._eval_expand_mul`` rebuilds ``n/d`` from ``fraction``, and the
# rebuilt numerator distributes a number over a sum.
_B = x**2 + 2*x + 1
SINGLE_DENOMINATOR = [
    (s.Mul(2, x**2/2 + x + s.Rational(1, 2), s.Pow(_B, -1), evaluate=False), s.S.One),
    (s.Mul(2, x/2 + 1, s.Pow(x + 2, -1), evaluate=False), s.S.One),
    (s.Mul(3, x/3 + y/3, s.Pow(x + y, -1), evaluate=False), s.S.One),
    ((2/_B)*s.Pow(2/_B, -1), s.S.One),
]


@pytest.mark.parametrize('e,expected', SINGLE_DENOMINATOR)
def test_expand_rebuilds_single_denominator_fraction(e, expected):
    assert s.srepr(s.expand(e)) == s.srepr(expected)
    assert s.Poly(e, t).domain == ZZ and s.Poly(e, t, field=True).domain == QQ


def test_polymatrix_polynomial_ring_argument_splits_domain_and_gens():
    R = QQ.frac_field(x)[t]
    m = pm.PolyMatrix([[1, 0], [0, 1]], ring=R)
    assert m.gens == (t,) and m.domain == QQ.frac_field(x) and m.ring == R
    assert m != pm.PolyMatrix([[1, 0], [0, 1]], t)  # QQ[t] is another ring
    q = pm.PolyMatrix([[1, 0], [0, (x - 1)/(x + 1)]], t)
    assert q.domain == ZZ.frac_field(x)
