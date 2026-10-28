"""Recursive sparse products preserve DMP grouping and stored binary values."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

SCALE = 1 << 28


def exact(value):
    expression = sx.to_sympy(value) if isinstance(value, sx.Expr) else sp.sympify(value)
    return expression.xreplace({f: sp.Rational(f) for f in expression.atoms(sp.Float)})


def domain(module, precision, kind, a):
    base = module.RealField(precision) if kind.startswith('real') else module.ComplexField(precision)
    if kind.endswith('ring'):
        return base.poly_ring(a)
    if kind.endswith('field'):
        return base.frac_field(a)
    return base


def inputs(module, precision, kind):
    x, y, a = module.symbols('x y a')
    K = domain(module, precision, kind, a)
    first = module.Float('1.23456789012345678901234567890123456789', precision=211)
    other = module.Float('0.7654321098765432109876543210987654321', precision=211)
    if not kind.startswith('real'):
        first += other*module.I
    if kind.endswith('ring'):
        first *= a+1
    elif kind.endswith('field'):
        first *= a/(a+1)
    F = {(2, 1): first, (1, 2): module.Float('-0.3333333333333333333333333333', precision=211),
         (1, 1): module.Float('1.0e18', precision=211), (0, 0): 1}
    G = {(1, 1): module.Float('0.875', precision=211),
         (1, 0): module.Float('0.12345678901234567890123456789', precision=211),
         (0, 1): module.Float('-0.875', precision=211), (0, 0): -1}
    return (x, y, a), K, F, G


@pytest.mark.parametrize('precision', [1, 53, 139, 140, 211])
@pytest.mark.parametrize('kind', ['complex', 'real_ring', 'complex_ring', 'real_field', 'complex_field'])
def test_sparse_products_squares_and_powers_match_dense_binary_certificates(precision, kind):
    (x, y, a), K, F, G = inputs(sx, precision, kind)
    (_, _, _), oracle_K, OF, OG = inputs(sp, precision, kind)
    X, Y = sp.symbols('x y')
    oracle_f = sp.Poly.from_dict(OF, X, Y, domain=oracle_K)
    oracle_g = sp.Poly.from_dict(OG, X, Y, domain=oracle_K)
    scaled = lambda terms: {tuple(SCALE*k for k in m): c for m, c in terms.items()}
    f = sx.Poly.from_dict(scaled(F), x, y, domain=K)
    g = sx.Poly.from_dict(scaled(G), x, y, domain=K)
    # The reference CC(a) path above53 bits has a documented preexisting
    # dtype precision loss. Preserve the native dense result for those three
    # contexts, and certify the reference gap independently below.
    dense_f = sx.Poly.from_dict(F, x, y, domain=K)
    dense_g = sx.Poly.from_dict(G, x, y, domain=K)
    retained_precision = kind == 'complex_field' and precision > 53
    operations = [(f*g, oracle_f*oracle_g, dense_f*dense_g),
                  (f*f, oracle_f*oracle_f, dense_f*dense_f)]
    operations += [(f**n, oracle_f**n, dense_f**n) for n in [0, 1, 2, 3, 4]]
    for actual, oracle, native_dense in operations:
        reference = native_dense if retained_precision else oracle
        assert actual.domain == K
        expected = scaled(reference.as_dict())
        assert set(actual.monoms()) == set(expected)
        for m, coefficient in actual.as_dict().items():
            # Compare exact rational images of every stored Float component.
            assert sp.cancel(exact(coefficient)-exact(expected[m])) == 0, (precision, kind, m)
        floats = {v._prec for c in actual.coeffs() for v in c.atoms(sx.Float)}
        expected_floats = {v._prec for c in oracle.coeffs() for v in c.atoms(sp.Float)}
        # Domain p and to_sympy's displayed precision differ at 1/139/211.
        assert floats == expected_floats
        assert (actual.domain if kind == 'complex' else actual.domain.dom).precision == precision


@pytest.mark.parametrize('coefficient', ['1', 'sin(a)', 'sqrt(2)'])
def test_exact_ex_product_and_power_have_symbolic_reconstruction_certificates(coefficient):
    x, y, a = sx.symbols('x y a')
    X, Y, A = sp.symbols('x y a')
    c = eval(coefficient, dict(vars(sx), a=a))
    C = sp.sympify(coefficient, locals={'a': A})
    f = sx.Poly.from_dict({(1 << 30, 0): 1, (0, 1 << 30): c, (0, 0): 1}, x, y, domain=sx.EX)
    g = sx.Poly.from_dict({(1 << 30, 0): 1, (0, 1 << 30): -c}, x, y, domain=sx.EX)
    F = sp.Poly(X+C*Y+1, X, Y, domain=sp.EX)
    G = sp.Poly(X-C*Y, X, Y, domain=sp.EX)
    for actual, reference in [(f*g, F*G), (f*f, F*F), (f**3, F**3)]:
        expected = {tuple((1 << 30)*k for k in m): v for m, v in reference.as_dict().items()}
        assert set(actual.monoms()) == set(expected)
        for m, value in actual.as_dict().items():
            assert sp.cancel(exact(value)-expected[m]) == 0


def test_last_legal_coordinate_and_true_storage_overflow_are_distinct():
    x = sx.Symbol('x')
    n = (1 << 31)-1
    for K in (sx.EX, sx.ComplexField(140)):
        p = sx.Poly.from_dict({(n,): 1, (0,): 1}, x, domain=K)
        for result in (p*p, p**2):
            assert {m: exact(c) for m, c in result.as_dict().items()} == {(2*n,): 1, (n,): 2, (0,): 1}
        q = sx.Poly.from_dict({(1 << 31,): 1, (0,): 1}, x, domain=K)
        for operation in (lambda: q*q, lambda: q**2):
            with pytest.raises(NotImplementedError, match='exponent storage limit'):
                operation()
        assert (q*0).is_zero and q**0 == sx.Poly(1, x, domain=K)


def test_input_conversion_callbacks_are_not_reentered_by_native_products():
    calls = []
    class Coefficient:
        def _sympy_(self):
            calls.append('convert')
            return sx.Float('1.234567890123456789', precision=211) + sx.I
    x, y = sx.symbols('x y')
    p = sx.Poly.from_dict({(SCALE, 0): Coefficient(), (0, SCALE): 1}, x, y, domain=sx.ComplexField(211))
    assert calls
    initial = tuple(calls)
    for result in (p*p, p**2, p**3):
        assert not result.is_zero
    assert tuple(calls) == initial


def test_sparse_products_do_not_import_sympy_at_runtime():
    source = '''import sys,importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0,Block())
import symox as s
x,y=s.symbols('x y');N=2**28
for K in (s.EX,s.ComplexField(211)):
    p=s.Poly.from_dict({(N,0):1,(0,N):2,(0,0):1},x,y,domain=K)
    assert (p*p).degree_list()==(2*N,2*N)
    assert (p**3).degree_list()==(3*N,3*N)
'''
    result = subprocess.run([sys.executable, '-B', '-c', source], capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    assert result.stderr == ''


@pytest.mark.parametrize('precision', [139, 140, 211])
def test_complex_fraction_domain_precision_retention_has_an_exact_rounding_certificate(precision):
    (x, y, a), K, F, G = inputs(sx, precision, 'complex_field')
    (X, Y, A), oracle_K, OF, OG = inputs(sp, precision, 'complex_field')
    f = sx.Poly.from_dict(F, x, y, domain=K)
    g = sx.Poly.from_dict(G, x, y, domain=K)
    oracle_f = sp.Poly.from_dict(OF, X, Y, domain=oracle_K)
    oracle_g = sp.Poly.from_dict(OG, X, Y, domain=oracle_K)
    for native, reference in [(f, oracle_f), (g, oracle_g)]:
        assert set(native.monoms()) == set(reference.monoms())
        for m, value in native.as_dict().items():
            assert sp.cancel(exact(value)-exact(reference.as_dict()[m])) == 0
    pairs = [(exact(c), exact(d)) for m, c in f.as_dict().items() for n, d in g.as_dict().items()
             if tuple(i+j for i, j in zip(m, n)) == (2, 1)]
    ideal = sum((c*d for c, d in pairs), sp.S.Zero)
    scale = sum((abs(sp.re(c.subs(A, 1)))+abs(sp.im(c.subs(A, 1))))*
                (abs(sp.re(d.subs(A, 1)))+abs(sp.im(d.subs(A, 1)))) for c, d in pairs)
    bound = sp.Rational(2)**(-precision+4)*max(scale, 1)
    def coefficient_error(value):
        # The common parameter denominator is exactly a+1. Two rounded
        # scalar operations contribute to each numerator coefficient;
        # 16*2**(-p) times the input product norm bounds their rounding.
        error = sp.Poly(sp.cancel((exact(value)-ideal)*(A+1)), A)
        return max([abs(sp.re(c))+abs(sp.im(c)) for c in error.all_coeffs()]+[sp.S.Zero])
    native_value = (f*g).as_dict()[(2, 1)]
    reference_value = (oracle_f*oracle_g).as_dict()[(2, 1)]
    assert coefficient_error(native_value) <= bound
    assert coefficient_error(reference_value) > bound
