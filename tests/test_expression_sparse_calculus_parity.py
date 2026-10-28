"""Sparse calculus keeps domain rounding and exact factorial coefficients."""
import math
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

N = 1 << 28


def exact(value):
    e = sx.to_sympy(value) if isinstance(value, sx.Expr) else sp.sympify(value)
    return e.xreplace({f: sp.Rational(f) for f in e.atoms(sp.Float)})


def domains(precision, kind):
    a, A = sx.Symbol('a'), sp.Symbol('a')
    native = sx.RealField(precision) if kind.startswith('real') else sx.ComplexField(precision)
    oracle = sp.RealField(precision) if kind.startswith('real') else sp.ComplexField(precision)
    if kind.endswith('ring'):
        return native.poly_ring(a), oracle.poly_ring(A)
    if kind.endswith('field'):
        return native.frac_field(a), oracle.frac_field(A)
    return native, oracle


def coefficient(kind):
    c = sx.Float('1.234567890123456789012345678901234567890123456789', precision=211)
    if not kind.startswith('real'):
        c += sx.Float('0.7654321098765432109876543210987654321', precision=211)*sx.I
    a = sx.Symbol('a')
    return c*(a+1) if kind.endswith('ring') else c*a/(a+1) if kind.endswith('field') else c


def equal_dict(actual, wanted):
    assert set(actual.as_dict()) == set(wanted)
    for m, value in actual.as_dict().items():
        assert sp.cancel(exact(value)-exact(wanted[m])) == 0


def reference_coefficient(oracle, value):
    """Load stored dyadics into the reference ground without evalf(p')."""
    base = oracle.dom if oracle.is_Composite else oracle
    def scalar(c):
        real, imag = c.as_real_imag()
        return base.dtype(real, imag) if base.is_CC else base.dtype(real)
    if not oracle.is_Composite:
        return scalar(exact(value))
    numerator, denominator = sp.fraction(sp.cancel(exact(value)))
    n, d = (sp.Poly(e, *oracle.symbols, domain=sp.QQ_I) for e in (numerator, denominator))
    # Canonical monic parameter denominator; the numerical dyadics are then
    # converted directly, without ComplexField.from_sympy's display evalf.
    leading = d.LC()
    n, d = (p.mul_ground(sp.expand_complex(1/leading)) for p in (n, d))
    ring = oracle.field.ring if oracle.is_Field else oracle.ring
    make = lambda p: ring.from_dict({m: scalar(c) for m, c in p.terms()})
    return oracle.dtype(make(n), make(d)) if oracle.is_Field else make(n)


def ground_operation(native, oracle, coefficient, factor, divide=False):
    value = reference_coefficient(oracle, coefficient)
    factor = oracle(factor)
    value = oracle.exquo(value, factor) if divide else oracle.mul(value, factor)
    rendered = oracle.to_sympy(value)
    # Polynomial construction retains the established Dom round trip, including
    # constant fraction denominators rendered at p' then converted at ground p.
    return native.from_sympy(sx.from_sympy(rendered)).as_expr()


@pytest.mark.parametrize('precision', [1, 53, 139, 140, 211])
@pytest.mark.parametrize('kind', ['complex', 'real_ring', 'complex_ring', 'real_field', 'complex_field'])
def test_large_sparse_diff_and_integrate_have_coefficient_and_domain_certificates(precision, kind):
    K, oracle = domains(precision, kind)
    x, y = sx.symbols('x y')
    f = sx.Poly.from_dict({(N+3, 2): coefficient(kind), (1, N+5): 2, (0, 0): 1}, x, y, domain=K)
    for j, variable in enumerate((x, y)):
        for order in (0, 1, 2, 5):
            derivative = f.diff((variable, order))
            expected = {}
            for m, c in f.as_dict().items():
                if m[j] < order:
                    continue
                k = math.prod(range(m[j]-order+1, m[j]+1))
                v = ground_operation(K, oracle, c, k)
                powers = list(m); powers[j] -= order
                expected[tuple(powers)] = v
            equal_dict(derivative, expected)
            assert derivative.domain == K and derivative.gens == f.gens
            primitive = f.integrate((variable, order))
            target = K.get_field()
            expected = {}
            for m, c in f.as_dict().items():
                k = math.prod(range(m[j]+1, m[j]+order+1))
                powers = list(m); powers[j] += order
                expected[tuple(powers)] = ground_operation(target, oracle.get_field(), c, k, divide=True)
            equal_dict(primitive, expected)
            assert primitive.domain == target and primitive.gens == f.gens


@pytest.mark.parametrize('precision,degree,order', [(1, 7, 2), (53, 5, 3), (139, 7, 2), (140, 7, 2), (211, 17, 4)])
def test_a_requested_derivative_order_matches_single_factor_reference_rounding(precision, degree, order):
    x = sx.Symbol('x'); X = sp.Symbol('x')
    value = sx.Float('1.234567890123456789012345678901234567890123456789', precision=211)
    f = sx.Poly(value*x**degree, x, domain=sx.ComplexField(precision))
    reference = sp.Poly(sx.to_sympy(value)*X**degree, X, domain=sp.ComplexField(precision))
    assert exact(f.LC()) == exact(reference.LC())
    assert exact(f.diff((x, order)).LC()) == exact(reference.diff((X, order)).LC())
    # Separate specs deliberately keep the reference's sequential rounding.
    assert exact(f.diff(*([x]*order)).LC()) == exact(reference.diff(*([X]*order)).LC())
    assert exact(f.diff((x, order)).LC()) != exact(f.diff(*([x]*order)).LC())


@pytest.mark.parametrize('precision', [1, 53, 139, 140, 211])
def test_integration_specs_keep_their_individual_grouping(precision):
    x = sx.Symbol('x'); X = sp.Symbol('x')
    value = sx.Float('1.234567890123456789012345678901234567890123456789', precision=211)
    f = sx.Poly(value*x**7, x, domain=sx.ComplexField(precision))
    reference = sp.Poly(sx.to_sympy(value)*X**7, X, domain=sp.ComplexField(precision))
    for specs, oracle_specs in [(((x, 3),), ((X, 3),)), ((x, x, x), (X, X, X)), (((x, 0),), ((X, 0),))]:
        assert exact(f.integrate(*specs).LC()) == exact(reference.integrate(*oracle_specs).LC())


@pytest.mark.parametrize('coefficient_source', ['1', 'sin(a)', 'Float("1.234567890123456789",precision=211)'])
def test_exact_ex_calculus_handles_factors_beyond_i64(coefficient_source):
    x, y, a = sx.symbols('x y a')
    c = sx.sympify(eval(coefficient_source, dict(vars(sx), a=a)))
    f = sx.Poly.from_dict({(1 << 31, 2): c, (0, 0): 3}, x, y, domain=sx.EX)
    factor = (1 << 31)*((1 << 31)-1)*((1 << 31)-2)
    result = f.diff((x, 3))
    reference = sp.EX.from_sympy(sx.to_sympy(c))*sp.EX(factor)
    equal_dict(result, {((1 << 31)-3, 2): sp.EX.to_sympy(reference)})
    q = f.integrate((x, 3))
    denominator = ((1 << 31)+1)*((1 << 31)+2)*((1 << 31)+3)
    expected = sp.EX.exquo(sp.EX.from_sympy(sx.to_sympy(c)), sp.EX(denominator))
    equal_dict(q, {((1 << 31)+3, 2): sp.EX.to_sympy(expected), (3, 0): sp.Rational(1, 2)})
    assert f.diff((x, (1 << 31)+1)).is_zero
    h = sx.Poly.from_dict({(40,): 1}, x, domain=sx.EX)
    equal_dict(h.diff((x, 21)), {(19,): math.prod(range(20, 41))})


@pytest.mark.parametrize('precision', [1, 139, 211])
@pytest.mark.parametrize('complex_', [False, True])
def test_auto_domains_zero_and_true_coordinate_limit(precision, complex_):
    a, x, y = sx.symbols('a x y')
    base = sx.ComplexField(precision) if complex_ else sx.RealField(precision)
    K = base.poly_ring(a)
    p = sx.Poly.from_dict({(7, 1): coefficient('complex_ring' if complex_ else 'real_ring')}, x, y, domain=K)
    assert p.integrate((x, 0)).domain == K.get_field()
    assert p.integrate((x, -1)).domain == K.get_field()
    assert p.integrate((x, 0), auto=False).domain == K
    assert p.integrate((x, 1), auto=False).domain == K  # divisor8 is exact
    z = sx.Poly(0, x, y, domain=K)
    assert z.integrate((x, (1 << 63)-1)).is_zero
    assert z.integrate((x, 0)).domain == K.get_field()
    q = sx.Poly.from_dict({((1 << 32)-1,): 1}, x, domain=K)
    with pytest.raises(NotImplementedError, match='exponent storage limit'):
        q.integrate()
    assert q.diff((x, (1 << 32))).is_zero


@pytest.mark.parametrize('name', ['EX', 'ComplexField(1)', 'ComplexField(139)', 'ComplexField(211)', 'ComplexField(140).frac_field(a)'])
def test_zero_evaluation_preserves_ground_and_scalar_zero_metadata(name):
    x, y, a = sx.symbols('x y a')
    K = eval(name, dict(vars(sx), a=a))
    p = sx.Poly.from_dict({(N, 2): 3, (0, N): 2, (0, 0): 1}, x, y, domain=K)
    result = p.eval(x, sx.Float(0, precision=211))
    equal_dict(result, {(N,): 2, (0,): 1})
    assert result.domain == K and result.gens == (y,)
    assert exact(result.eval(y, 0)) == 1
    q = sx.Poly.from_dict({(N,): 1}, x, domain=K)
    zero = q.eval(0)
    reference_K = eval(name, dict(vars(sp), a=sp.Symbol('a')))
    reference = sp.Poly(sp.Symbol('x')**2, sp.Symbol('x'), domain=reference_K).eval(0)
    assert exact(zero) == exact(reference)
    assert {v._prec for v in zero.atoms(sx.Float)} == {v._prec for v in reference.atoms(sp.Float)}


def test_invalid_generator_is_checked_even_when_the_order_is_zero():
    x, y = sx.symbols('x y')
    p = sx.Poly(x, x, domain=sx.CC)
    for operation in (lambda: p.diff((y, 0)), lambda: p.integrate((y, 0))):
        with pytest.raises((ValueError, sx.PolynomialError)):
            operation()


def test_sparse_calculus_has_no_runtime_sympy_dependency():
    source = '''import sys,importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0,Block())
import symox as s
x,y=s.symbols('x y');N=2**28
for K in (s.EX,s.ComplexField(211)):
    p=s.Poly.from_dict({(N,2):1,(0,N):2},x,y,domain=K)
    assert p.diff((x,3)).degree()==N-3
    assert p.integrate((x,2)).degree()==N+2
    assert p.eval(x,0).gens==(y,)
'''
    result = subprocess.run([sys.executable, '-B', '-c', source], capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    assert result.stderr == ''


@pytest.mark.parametrize('precision', [1, 53, 139, 140, 211])
def test_compact_multivariate_complex_calculus_matches_reference_binary_inputs(precision):
    x,y=sx.symbols('x y'); X,Y=sp.symbols('x y')
    K=sx.ComplexField(precision); O=sp.ComplexField(precision)
    f=sx.Poly.from_dict({(7,3):coefficient('complex'),(2,4):sx.I/3,(0,0):2},x,y,domain=K)
    g=sp.Poly.from_dict({m:sx.to_sympy(c) for m,c in f.as_dict().items()},X,Y,domain=O)
    equal_dict(f,g.as_dict())  # Identical stored binary inputs before arithmetic.
    for j,variable in enumerate((x,y)):
        for order in (0,1,2,3,21):
            equal_dict(f.diff((variable,order)),g.diff((g.gens[j],order)).as_dict())
            equal_dict(f.integrate((variable,order)),g.integrate((g.gens[j],order)).as_dict())


@pytest.mark.parametrize('coefficient_source', ['sin(a)', 'sqrt(2)', 'Float("1.234567890123456789",precision=139)'])
def test_expression_coefficients_and_dependent_generators_keep_reference_calculus(coefficient_source):
    a,x=sx.symbols('a x'); A,X=sp.symbols('a x')
    c=sx.sympify(eval(coefficient_source,dict(vars(sx),a=a)))
    gens=(sx.sin(x),x+1); oracle_gens=(sp.sin(X),X+1)
    f=sx.Poly.from_dict({(5,2):c,(1,3):sx.Rational(2,3),(0,0):1},*gens,domain=sx.EX)
    g=sp.Poly.from_dict({m:sx.to_sympy(v) for m,v in f.as_dict().items()},*oracle_gens,domain=sp.EX)
    for j,variable in enumerate(gens):
        for order in (0,1,3):
            equal_dict(f.diff((variable,order)),g.diff((oracle_gens[j],order)).as_dict())
            equal_dict(f.integrate((variable,order)),g.integrate((oracle_gens[j],order)).as_dict())
    equal_dict(f.eval(gens[0],0),g.eval(oracle_gens[0],0).as_dict())


def test_calculus_input_protocols_and_errors_are_observable_once():
    x=sx.Symbol('x'); events=[]
    class Input:
        def __init__(self,value): self.value=value
        def _sympy_(self): events.append(self.value);return self.value
    p=sx.Poly.from_dict({(N,):Input(sx.Rational(1,3))},x,domain=sx.EX)
    assert events==[sx.Rational(1,3)]
    p.diff((x,3));p.integrate((x,2));p.eval(0)
    assert events==[sx.Rational(1,3)]
    failure=RuntimeError('coefficient callback')
    class Bad:
        def _sympy_(self): raise failure
    with pytest.raises(RuntimeError) as caught:
        sx.Poly.from_dict({(N,):Bad()},x,domain=sx.EX)
    assert caught.value is failure


@pytest.mark.parametrize('precision', [1,53,139,140,211])
def test_big_exact_factor_with_expression_float_uses_reference_rounding(precision):
    x=sx.Symbol('x'); n=1<<31
    c=sx.Float('1.234567890123456789123456789',precision=precision)
    f=sx.Poly.from_dict({(n,):c},x,domain=sx.EX)
    factor=n*(n-1)*(n-2)
    expected=sp.EX(factor)*sp.EX(sx.to_sympy(c))
    equal_dict(f.diff((x,3)),{(n-3,):sp.EX.to_sympy(expected)})




def test_ex_cancellation_finishes_before_the_output_coordinate_bound():
    a,x=sx.symbols('a x')
    c=(a*a-1)/(a-1)-a-1
    p=sx.Poly.from_dict({((1<<32)-1,):c,(0,):3},x,domain=sx.EX)
    assert not p.is_zero and c.is_zero is None
    equal_dict(p.integrate(),{(1,):3})
    q=sx.Poly.from_dict({((1<<32)-1,):a},x,domain=sx.EX)
    with pytest.raises(NotImplementedError,match='exponent storage limit'):
        q.integrate()
