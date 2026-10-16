"""Exact sparse ring division and dense Poly division use distinct contracts."""
import pickle
import subprocess
import sys

import pytest
import symox as sm


def rings(ground, order):
    import sympy as sy
    a = sm.Symbol('a', positive=True)
    b = sm.Dummy('b', real=True)
    params = (a, b)
    k = getattr(sm, ground).poly_ring(*params, order=order)
    oracle = getattr(sy, ground).poly_ring(*(sm.to_sympy(v) for v in params), order=order)
    return k, oracle, a, b


def same_element(k, actual, oracle, expected):
    import sympy as sy
    assert actual.parent() == k
    assert sy.expand(sm.to_sympy(k.to_sympy(actual)) - oracle.to_sympy(expected)) == 0
    assert k.of_type(actual)
    assert pickle.loads(pickle.dumps(actual)) == actual


@pytest.mark.parametrize('ground', ['ZZ', 'QQ'])
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
@pytest.mark.parametrize('case', range(6))
def test_ring_division_oracle_operators_and_independent_certificate(ground, order, case):
    k, oracle, a, b = rings(ground, order)
    f, g = [(a*a+2*a, 2*a), (a*b+1, 2*a+2*b),
            ((a+2*b*b)*(a+b)+1, a+2*b*b),
            ((a+b)*(a-b), a-b), (0, a+b), (a+b, 3)][case]
    u, v = k.convert(f), k.convert(g)
    ou, ov = oracle.convert(sm.to_sympy(sm.sympify(f))), oracle.convert(sm.to_sympy(sm.sympify(g)))
    expected = oracle.div(ou, ov)
    q, r = k.div(u, v)
    for actual, value in zip((q, r), expected):
        same_element(k, actual, oracle, value)
    assert u == q*v+r
    assert k.div(r, v) == (k.zero, r)
    assert k.quo(u, v) == q and k.rem(u, v) == r
    assert u.div(v) == divmod(u, v) == (q, r)
    assert u.quo(v) == u//v == q and u.rem(v) == u%v == r
    if r == k.zero:
        assert k.exquo(u, v) == u.exquo(v) == u/v == q
    else:
        for operation in [lambda: k.exquo(u, v), lambda: u.exquo(v), lambda: u/v]:
            with pytest.raises(sm.ExactQuotientFailed):
                operation()
    assert k.numer(u) == u and k.denom(u) == k.one
    assert u.numerator == u.numer == u and u.denominator == u.denom == k.one
    for operation in [lambda: k.div(u, k.zero), lambda: u//k.zero, lambda: u%0, lambda: u/0]:
        with pytest.raises(ZeroDivisionError):
            operation()


@pytest.mark.parametrize('ground', ['ZZ', 'QQ'])
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
def test_ring_reversed_scalars_and_parameter_order(ground, order):
    k, oracle, a, b = rings(ground, order)
    u = k.convert(2*a+2)
    q, r = k.div(k(7), u)
    assert divmod(7, u) == (q, r)
    assert 7//u == q and 7%u == r
    assert (7//u)*u+7%u == k(7)
    q, r = k.div(u, k(2))
    assert u//2 == q and u%2 == r and u/2 == q
    reversed_ring = getattr(sm, ground).poly_ring(b, a, order=order)
    other = reversed_ring.convert(2*a+2)
    with pytest.raises(sm.CoercionFailed):
        divmod(u, other)
    # Domain methods explicitly convert their operands; native operators
    # reject implicit mixing of two different typed parent rings.
    assert k.div(u, other) == (k.one, k.zero)
    assert k.symbols == (a, b) and q.parent().symbols == (a, b)
    assert q.parent().order == order


@pytest.mark.parametrize('ground', ['ZZ', 'QQ'])
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
def test_ring_gcd_lcm_declared_order_and_divisibility(ground, order):
    k, oracle, a, b = rings(ground, order)
    f = a+2*b*b
    for left, right in [(f, f*(a+b)), (2*a, 4*a), (0, 3), (sm.Rational(2,3)*a, sm.Rational(4,5)*a)] if ground == 'QQ' else [(f, f*(a+b)), (2*a, 4*a), (0, 3)]:
        u, v = k.convert(left), k.convert(right)
        ou, ov = oracle.convert(sm.to_sympy(sm.sympify(left))), oracle.convert(sm.to_sympy(sm.sympify(right)))
        gcd, lcm = k.gcd(u, v), k.lcm(u, v)
        same_element(k, gcd, oracle, oracle.gcd(ou, ov))
        same_element(k, lcm, oracle, oracle.lcm(ou, ov))
        for value in [u, v]:
            if gcd:
                assert k.exquo(value, gcd)*gcd == value
            if value:
                assert k.exquo(lcm, value)*value == lcm


@pytest.mark.parametrize('kind', ['QQ', 'RF_QQ', 'RF_ZZ', 'GF', 'AN'])
@pytest.mark.parametrize('case', range(4))
def test_field_gcdex_returns_gcd_and_exact_bezout(kind, case):
    a = sm.Symbol('a', positive=True)
    k = {'QQ': sm.QQ, 'RF_QQ': sm.QQ.frac_field(a), 'RF_ZZ': sm.ZZ.frac_field(a),
         'GF': sm.GF(7), 'AN': sm.QQ.algebraic_field(sm.sqrt(2))}[kind]
    pairs = [(2, 3), (0, 3), (2, 0), (0, 0)]
    if kind == 'QQ':
        pairs[0] = (sm.Rational(2,3), sm.Rational(4,5))
    if kind.startswith('RF'):
        pairs[0] = (2*a/(a+1), 3*a/(a+2))
    u, v = (k.convert(value) for value in pairs[case])
    s, t, h = k.gcdex(u, v)
    assert s*u+t*v == h
    # GF/AN gcd(0,0) in SymPy 1.14 returns 1, which has no Bezout certificate.
    assert h == (k.zero if case == 3 else k.gcd(u, v))
    assert all(value.parent() == k for value in [s, t, h])


@pytest.mark.parametrize('ground', ['ZZ', 'QQ'])
def test_polynomial_ring_gcdex_and_inverse_supported_ground(ground):
    import sympy as sy
    x = sm.Symbol('x')
    k = getattr(sm, ground).poly_ring(x)
    u, v = k.convert(x+1), k.convert(x*x+1)
    if ground == 'ZZ':
        with pytest.raises(sm.DomainError):
            k.gcdex(u, v)
        with pytest.raises(sm.DomainError):
            k.invert(u, v)
    else:
        s, t, h = k.gcdex(u, v)
        assert s*u+t*v == h == k.one
        inverse = k.invert(u, v)
        assert (inverse*u)%v == k.one
        oracle = sy.QQ.poly_ring(sm.to_sympy(x))
        same_element(k, inverse, oracle, oracle.invert(oracle.convert(sm.to_sympy(x+1)), oracle.convert(sm.to_sympy(x*x+1))))


@pytest.mark.parametrize('ground', ['ZZ', 'QQ', 'ZZring', 'QQring', 'ZZfield', 'QQfield'])
@pytest.mark.parametrize('auto', [False, True])
@pytest.mark.parametrize('case', range(5))
def test_poly_division_auto_normalization_and_recursive_contract(ground, auto, case):
    import sympy as sy
    x, y = sm.symbols('x y')
    a = sm.Symbol('a', positive=True)
    k = {'ZZ': sm.ZZ, 'QQ': sm.QQ, 'ZZring': sm.ZZ.poly_ring(a), 'QQring': sm.QQ.poly_ring(a), 'ZZfield': sm.ZZ.frac_field(a), 'QQfield': sm.QQ.frac_field(a)}[ground]
    ok = {'ZZ': sy.ZZ, 'QQ': sy.QQ, 'ZZring': sy.ZZ.poly_ring(sm.to_sympy(a)), 'QQring': sy.QQ.poly_ring(sm.to_sympy(a)), 'ZZfield': sy.ZZ.frac_field(sm.to_sympy(a)), 'QQfield': sy.QQ.frac_field(sm.to_sympy(a))}[ground]
    f, g = [(x+1, 2*x+2), (x*x+2*x, 2*x), (x*y+2*x, 2*x), (2*x*x+2*x+1, 2*x+2), (x*y+x,y)][case]
    if 'ring' in ground or 'field' in ground:
        f, g = a*f, a*g
    p, q = sm.Poly(f, x, y, domain=k), sm.Poly(g, x, y, domain=k)
    op, oq = sy.Poly(sm.to_sympy(f), sm.to_sympy(x), sm.to_sympy(y), domain=ok), sy.Poly(sm.to_sympy(g), sm.to_sympy(x), sm.to_sympy(y), domain=ok)
    quot, rem = p.div(q, auto=auto)
    assert sm.expand(quot.as_expr()*q.as_expr()+rem.as_expr()-p.as_expr()) == 0
    for method in ['div', 'quo', 'rem', 'exquo']:
        try:
            expected = getattr(op, method)(oq, auto=auto)
        except sy.polys.polyerrors.ExactQuotientFailed:
            with pytest.raises(sm.ExactQuotientFailed):
                getattr(p, method)(q, auto=auto)
            continue
        actual = getattr(p, method)(q, auto=auto)
        for v, ov in zip(actual if method == 'div' else (actual,), expected if method == 'div' else (expected,)):
            assert sy.expand(sm.to_sympy(v.as_expr())-ov.as_expr()) == 0
            assert str(v.domain) == str(ov.domain)
            assert all(c in v.domain for c in v.all_coeffs()) if len(v.gens) == 1 else all(c in v.domain for c in v.coeffs())


@pytest.mark.parametrize('ground', ['ZZ', 'QQ'])
@pytest.mark.parametrize('case', range(7))
def test_poly_ring_lcm_cofactors_and_exact_quotients(ground, case):
    import sympy as sy
    x = sm.Symbol('x')
    a = sm.Symbol('a', positive=True)
    k = getattr(sm, ground).poly_ring(a, order='grevlex')
    ok = getattr(sy, ground).poly_ring(sm.to_sympy(a), order='grevlex')
    left, right = [(a*x, a*x*x), (2*a*x, 4*a*x*x), (2*(a+1)*x, 4*(a+1)*x*x),
                   (2*x+1, 4*x+2), (x+1, 2*x+2), (0, 2*(a+1)*x),
                   (2*(a+1)*(x+1), 4*(a+1)*(x+1)**2)][case]
    p, q = sm.Poly(left, x, domain=k), sm.Poly(right, x, domain=k)
    op, oq = sy.Poly(sm.to_sympy(left), sm.to_sympy(x), domain=ok), sy.Poly(sm.to_sympy(right), sm.to_sympy(x), domain=ok)
    g, u, v = p.cofactors(q)
    lcm = p.lcm(q)
    for actual, expected in zip([g, u, v, lcm], [*op.cofactors(oq), op.lcm(oq)]):
        assert actual.domain == k
        assert sy.expand(sm.to_sympy(actual.as_expr())-expected.as_expr()) == 0
        assert all(c in k for c in actual.coeffs())
    assert sm.expand(g.as_expr()*u.as_expr()-p.as_expr()) == 0
    assert sm.expand(g.as_expr()*v.as_expr()-q.as_expr()) == 0
    for value in [p, q]:
        if value:
            quotient = lcm.exquo(value, auto=False)
            assert quotient.domain == k
            assert sm.expand(quotient.as_expr()*value.as_expr()-lcm.as_expr()) == 0


def test_native_division_without_sympy_runtime():
    code = '''
import importlib.abc,sys,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise ImportError('SymPy blocked')
sys.meta_path.insert(0,Block())
import symox as s
x=s.Symbol('x');a=s.Dummy('a',positive=True);b=s.Symbol('b',real=True)
for ground in [s.ZZ,s.QQ]:
 for order in ['lex','grlex','grevlex']:
  k=ground.poly_ring(a,b,order=order);u=k.convert(a*a+2*a);v=k.convert(2*a)
  q,r=divmod(u,v);assert u==q*v+r and q.parent()==k
  assert pickle.loads(pickle.dumps(q)).parent()==k
  p=s.Poly(2*(a+1)*x,x,domain=k);z=s.Poly(4*(a+1)*x*x,x,domain=k)
  g,c,d=p.cofactors(z);assert g*c==p and g*d==z
  l=p.lcm(z);assert l.exquo(p,auto=False).domain==k
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result = subprocess.run([sys.executable, '-c', code], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr


def test_poly_mixed_ground_promotion_and_incompatible_parent_safety():
    x = sm.Symbol('x')
    a = sm.Symbol('a', positive=True)
    b = sm.Dummy('b', real=True)
    k = sm.ZZ.poly_ring(a, b, order='grlex')
    qk = sm.QQ.poly_ring(a, b, order='grlex')
    p = sm.Poly(a*x+a, x, domain=k)
    q = sm.Poly(a*x, x, domain=qk)
    for result in [p+q, p*q, p.gcd(q), p.lcm(q), *p.cofactors(q), *p.div(q, auto=False)]:
        assert result.domain == qk
        assert result.domain.symbols == (a, b)
        assert result.domain.order == 'grlex'
        assert all(c in qk for c in result.coeffs())
    incompatible = sm.Poly(a*x, x, domain=sm.QQ.poly_ring(a,b,order='grevlex'))
    for operation in [lambda:p+incompatible, lambda:p*incompatible, lambda:p.div(incompatible),
                      lambda:p.lcm(incompatible), lambda:p.cofactors(incompatible)]:
        with pytest.raises(TypeError):
            operation()


@pytest.mark.parametrize('a,b', [(-7,3),(7,-3),(-7,-3),(7,3)])
def test_integer_euclidean_operators_preserve_floor_identity(a,b):
    u, v = sm.ZZ(a), sm.ZZ(b)
    q, r = divmod(u,v)
    assert (int(q), int(r)) == divmod(a,b)
    assert u == q*v+r
    assert u//b == q and a//v == q
    assert u%b == r and a%v == r
    assert q.parent() == r.parent() == sm.ZZ
