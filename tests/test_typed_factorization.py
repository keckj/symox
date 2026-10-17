"""Typed factor/square-free lists apply their options and prove reconstruction."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

from test_global_poly_division import selected

OPS = ['factor_list', 'sqf_list', 'sqf_part']
KINDS = ['ZZ', 'QQ', 'GF', 'AN', 'QQRF', 'ZZRF', 'QQring', 'ZZring']


def expression_of(value):
    return value.as_expr() if isinstance(value, sm.Poly) else value


def compare_result(actual, expected):
    import sympy as sy
    if isinstance(expected, (tuple, list)):
        assert isinstance(actual, type(expected)) and len(actual) == len(expected)
        for a, e in zip(actual, expected):
            compare_result(a, e)
    elif isinstance(expected, sy.Poly):
        assert isinstance(actual, sm.Poly)
        assert tuple(sm.to_sympy(g) for g in actual.gens) == expected.gens
        assert str(actual.domain) == str(expected.domain)
        assert sy.expand(sm.to_sympy(actual.as_expr())-expected.as_expr()) == 0
    else:
        assert not isinstance(actual, sm.Poly)
        assert sy.expand(sm.to_sympy(actual)-expected) == 0


def reconstruct(p, result):
    unit, factors = result
    product = sm.Poly(unit, *p.gens, domain=p.domain)
    for factor, multiplicity in factors:
        assert isinstance(multiplicity, int) and multiplicity > 0
        f = factor if isinstance(factor, sm.Poly) else sm.Poly(factor, *p.gens, domain=p.domain)
        assert f.domain == p.domain and f.gens == p.gens
        assert all(c in p.domain for c in f.coeffs())
        restored = pickle.loads(pickle.dumps(f))
        assert restored.domain == p.domain and restored.gens == p.gens
        assert restored.as_expr() == f.as_expr()
        product *= f**multiplicity
    assert sm.Poly(sm.expand(product.as_expr()-p.as_expr()), *p.gens, domain=p.domain).is_zero


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
@pytest.mark.parametrize('case', range(4))
def test_typed_method_and_global_factors_oracle_metadata_and_certificates(kind, order, case):
    import sympy as sy
    x = sm.Symbol('x', real=True)
    k, ok, a, unused = selected(kind, order)
    c = a/(unused+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    expression = [sm.Integer(0), 3*c, (2*x+1)**2*(x-2), c*(c*x+1)**2*(x-3)**3][case]
    p = sm.Poly(expression, x, domain=k)
    op = sy.Poly(sm.to_sympy(expression), sm.to_sympy(x), domain=ok)
    for name in OPS:
        actual = getattr(p, name)()
        if p.is_zero and kind == 'GF' and name in ['sqf_list','sqf_part']:
            # SymPy 1.14 raises AssertionError on the zero GF polynomial.
            if name=='sqf_part':
                assert actual.is_zero and actual.domain == p.domain
                assert sm.sqf_part(p,polys=True).is_zero
            else:
                assert actual==(0,[])
                assert sm.sqf_list(p,polys=True)==(0,[])
            continue
        expected = getattr(op, name)()
        compare_result(actual, expected)
        if name.endswith('list'):
            reconstruct(p, actual)
        else:
            assert actual.domain == k and actual.gens == p.gens
            assert all(c in k for c in actual.coeffs())
            if not p.is_zero:
                assert p.rem(actual).is_zero
                assert actual.gcd(actual.diff()).degree() == 0
        for value, oracle in [(p, op), (expression, sm.to_sympy(expression))]:
            result = getattr(sm, name)(value, x, domain=k, polys=True)
            expected = getattr(sy, name)(oracle, sm.to_sympy(x), domain=ok, polys=True)
            if kind=='GF' and name.endswith('list'):
                # SymPy's symbolic Expr path can keep an integer content
                # outside GF(p): 8 and -2 denote the same GF(5) coefficient.
                assert ok.convert(sm.to_sympy(result[0]))==ok.convert(expected[0])
                compare_result(result[1],expected[1])
                reconstruct(p,result)
            else:
                compare_result(result, expected)
        polys = actual[1] if name.endswith('list') else [(actual, 1)]
        for factor, _ in polys:
            if kind.endswith('RF') or kind.endswith('ring'):
                assert factor.domain.symbols == (a, unused)
                assert factor.domain.order == order


@pytest.mark.parametrize('name', OPS)
@pytest.mark.parametrize('input_kind', ['expression', 'poly'])
@pytest.mark.parametrize('polys', [None, False, True])
def test_global_result_rendering(name, input_kind, polys):
    import sympy as sy
    x, y = sm.symbols('x y')
    expression = sm.Rational(1,3)*(2*x+y)**2*(x-1)
    value = sm.Poly(expression, x, y) if input_kind == 'poly' else expression
    oracle = sy.Poly(sm.to_sympy(expression), sm.to_sympy(x), sm.to_sympy(y)) if input_kind == 'poly' else sm.to_sympy(expression)
    actual = getattr(sm, name)(value, x, y, polys=polys)
    expected = getattr(sy, name)(oracle, sm.to_sympy(x), sm.to_sympy(y), polys=polys)
    compare_result(actual, expected)


@pytest.mark.parametrize('kind', ['ZZ', 'QQ', 'AN', 'QQRF', 'QQring', 'ZZring'])
@pytest.mark.parametrize('coefficient', [1, sm.Rational(1,3)])
def test_factor_list_and_squarefree_list_have_distinct_normalization(kind, coefficient):
    import sympy as sy
    x = sm.Symbol('x')
    k, ok, _, _ = selected(kind)
    if coefficient != 1 and kind in ['ZZ', 'ZZring']:
        return
    expression = coefficient*(2*x+1)**2
    p = sm.Poly(expression, x, domain=k)
    op = sy.Poly(sm.to_sympy(expression), sm.to_sympy(x), domain=ok)
    for name in OPS:
        compare_result(getattr(p, name)(), getattr(op, name)())
    if kind in ['QQ', 'QQRF']:
        assert p.factor_list()[1][0][0].as_expr() == 2*x+1
        assert p.sqf_list()[1][0][0].as_expr() == x+sm.Rational(1,2)
        assert p.sqf_part().as_expr() == x+sm.Rational(1,2)
    elif kind.endswith('ring'):
        assert p.sqf_part().as_expr() == 2*x+1
    for included in [p.factor_list_include(), p.sqf_list_include()]:
        product = sm.Poly(1, x, domain=k)
        for factor, multiplicity in included:
            product *= factor**multiplicity
        assert sm.Poly(sm.expand(product.as_expr()-p.as_expr()), x, domain=k).is_zero


@pytest.mark.parametrize('kind', KINDS)
def test_squarefree_all_slots_are_actual_domain_units(kind):
    import sympy as sy
    x = sm.Symbol('x')
    k, ok, _, _ = selected(kind)
    p = sm.Poly(2*(2*x+1)**3*(x-3)**5, x, domain=k)
    result = p.sqf_list(all=True)
    assert [multiplicity for _, multiplicity in result[1]] == list(range(1, 6))
    reconstruct(p, result)
    assert all(f.domain == k for f, _ in result[1])
    if kind != 'GF':
        op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x), domain=ok)
        compare_result(result, op.sqf_list(all=True))


@pytest.mark.parametrize('prime,exponent', [(2, 8), (3, 9), (5, 10)])
def test_finite_field_inseparability_has_correct_distinct_factors(prime, exponent):
    import sympy as sy
    x = sm.Symbol('x')
    p = sm.Poly((x*x+1)**exponent*(x+1), x, modulus=prime)
    op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x), modulus=prime)
    for name in OPS:
        actual = getattr(sm, name)(p, polys=True)
        compare_result(actual, getattr(op, name)())
        if name.endswith('list'):
            reconstruct(p, actual)
    part = p.sqf_part()
    assert part.gcd(part.diff()).degree() == 0
    assert p.rem(part).is_zero
    if prime == 2:
        assert part.as_expr() == x+1
    actual = sm.sqf_part((x*x+1)**2, x, modulus=2)
    assert actual == x+1


@pytest.mark.parametrize('name', OPS)
@pytest.mark.parametrize('domain', ['QQ', 'GF', 'AN', 'RF', 'ring'])
def test_poly_construction_options_are_really_applied(name, domain):
    import sympy as sy
    x = sm.Symbol('x')
    a = sm.Symbol('a', positive=True)
    p = sm.Poly((x*x-2)**2, x)
    op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x))
    if domain == 'QQ':
        kw, okw = {'domain':sm.QQ}, {'domain':sy.QQ}
    elif domain == 'GF':
        kw, okw = {'modulus':2}, {'modulus':2}
    elif domain == 'AN':
        kw, okw = {'extension':sm.sqrt(2)}, {'extension':sy.sqrt(2)}
    elif domain == 'RF':
        kw, okw = {'domain':sm.QQ.frac_field(a)}, {'domain':sy.QQ.frac_field(sm.to_sympy(a))}
    else:
        kw, okw = {'domain':sm.ZZ.poly_ring(a)}, {'domain':sy.ZZ.poly_ring(sm.to_sympy(a))}
    compare_result(getattr(sm, name)(p, polys=True, **kw), getattr(sy, name)(op, polys=True, **okw))


@pytest.mark.parametrize('name', OPS)
@pytest.mark.parametrize('keyword', ['nonsense', 'auto', 'all', 'include', 'polys'])
def test_kwargs_validation(name, keyword):
    x = sm.Symbol('x')
    value = 2 if keyword == 'polys' else True
    expected = sm.FlagError if keyword in ['auto', 'all', 'include'] else sm.OptionError
    with pytest.raises(expected):
        getattr(sm, name)(x, x, **{keyword:value})


@pytest.mark.parametrize('name', ['factor_list', 'sqf_list'])
@pytest.mark.parametrize('value', [False, None, True, 2])
def test_frac_flag_has_an_explicit_contract(name, value):
    x = sm.Symbol('x')
    if value is True:
        assert getattr(sm, name)(1/(x+1), x, frac=True) == (1, [], [(x+1, 1)])
    elif value == 2:
        with pytest.raises(sm.OptionError):
            getattr(sm, name)(x+1, x, frac=value)
    else:
        assert getattr(sm, name)(x+1, x, frac=value) == (1, [(x+1, 1)])


@pytest.mark.parametrize('name',['factor_list','sqf_list'])
def test_default_frac_false_does_not_treat_denominators_as_generators(name):
    x=sm.Symbol('x')
    for expression in [1/x,(x*x-1)/(x+1),sm.exp(-x)]:
        with pytest.raises(ValueError,match='polynomial expected'):
            getattr(sm,name)(expression)


@pytest.mark.parametrize('name', OPS)
def test_incompatible_options_and_explicit_generators(name):
    x, y = sm.symbols('x y')
    for kwargs in [{'domain':sm.QQ, 'modulus':2}, {'domain':sm.QQ, 'extension':sm.sqrt(2)},
                   {'domain':sm.QQ, 'field':True}, {'modulus':2, 'extension':sm.I}]:
        with pytest.raises(sm.OptionError):
            getattr(sm, name)(x*x+1, x, **kwargs)
    with pytest.raises(sm.OptionError):
        getattr(sm, name)(x, x, gens=(y,))
    result = getattr(sm, name)((x+y)**2, gens=(y, x), polys=True)
    factors = result[1] if name.endswith('list') else [(result, 1)]
    assert all(p.gens == (y, x) for p, _ in factors)


def test_legacy_expression_coefficients_preserve_mathematically_correct_results():
    x, a = sm.symbols('x a')
    expression = sm.sin(a)**2*(x+1)**2*(x-2)
    p = sm.Poly(expression, x)
    for result in [p.factor_list(), p.sqf_list()]:
        value = result[0]
        for factor, power in result[1]:
            value *= factor.as_expr()**power
        assert sm.expand(value-expression) == 0
    assert p.sqf_part().as_expr() == (x+1)*(x-2) or sm.expand(p.sqf_part().as_expr()-(x+1)*(x-2)) == 0
    assert sm.sqf_part(expression, x, polys=True).as_expr() == p.sqf_part().as_expr()


@pytest.mark.parametrize('name', ['factor_list', 'sqf_list'])
@pytest.mark.parametrize('constant', [0, 3, sm.Rational(2,3), -7])
def test_numeric_constants_without_generators(name, constant):
    import sympy as sy
    compare_result(getattr(sm, name)(constant), getattr(sy, name)(sm.to_sympy(constant)))


@pytest.mark.parametrize('name', ['factor_list','sqf_list'])
@pytest.mark.parametrize('case', range(5))
def test_symbolic_expression_content_and_individual_factor_domains(name,case):
    import sympy as sy
    x,a,b=sm.symbols('x a b')
    expression=[x/a,(a*x+1)**2/a**2,sm.Rational(1,3)*(a*x+1)**2,
                (a*x+1)*(b*x+1)/((a+1)*(b+1)),x/3+sm.Rational(1,2)][case]
    actual=getattr(sm,name)(expression,x,polys=True)
    expected=getattr(sy,name)(sm.to_sympy(expression),sm.to_sympy(x),polys=True)
    compare_result(actual,expected)
    value=actual[0]
    for factor,power in actual[1]:value*=factor.as_expr()**power
    assert sm.cancel(value-expression)==0


@pytest.mark.parametrize('name',['factor_list','sqf_list'])
@pytest.mark.parametrize('polys',[False,True])
@pytest.mark.parametrize('constant',[sm.I,2+sm.I,sm.sqrt(2),2*sm.sqrt(2),sm.pi,sm.exp(2)])
def test_symbolic_constant_coefficients_and_constant_generators(name,polys,constant):
    import sympy as sy
    compare_result(getattr(sm,name)(constant,polys=polys),getattr(sy,name)(sm.to_sympy(constant),polys=polys))


@pytest.mark.parametrize('name',['factor_list'])
@pytest.mark.parametrize('polys',[False,True])
@pytest.mark.parametrize('explicit_gens',[False,True])
@pytest.mark.parametrize('case',range(10))
def test_natural_factor_order_uses_source_component_polynomials(name,polys,explicit_gens,case):
    import sympy as sy
    x,y,z,w=sm.symbols('x y z w')
    expressions=[x**4*y-x*x*y**3,x*y+x,x*x*y-y,x*x*y+x*y*y,-716*x*z,
                 432*w*w*x*x*y**6,sm.expand((x*x-1)*(y*y-1)),
                 (x*x-1)*(y*y-1),(x+1)*y+(x+1)*z,
                 (x+1)*sm.expand((x+1)*(y+1))]
    expression=expressions[case]
    gens=(x,y,z,w) if explicit_gens else ()
    actual=getattr(sm,name)(expression,*gens,polys=polys)
    expected=getattr(sy,name)(sm.to_sympy(expression),*(sm.to_sympy(g) for g in gens),polys=polys)
    compare_result(actual,expected)
    value=actual[0]
    for factor,power in actual[1]:value*=expression_of(factor)**power
    assert sm.expand(value-expression)==0


@pytest.mark.parametrize('name', OPS)
def test_preflight_rejects_expansion_before_allocation(name):
    x = sm.Symbol('x')
    for expression in [sm.Pow(x+1, 10**9, evaluate=False), sm.Pow(2, 10**9, evaluate=False),
                       sm.Pow(sum(sm.symbols('a b c d e f g h')), 1000, evaluate=False)]:
        with pytest.raises(ValueError, match='resource limit before expansion'):
            getattr(sm, name)(expression, x)
    result = getattr(sm, name)((x+1)**3, x)
    assert result is not None


def test_no_sympy_runtime_and_typed_pickle():
    program = r'''
import sys, importlib.abc, pickle
class NoSymPy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0, NoSymPy())
import symox as s
x=s.Dummy('x', real=True)
a=s.Symbol('a', positive=True)
unused=s.Dummy('unused')
for domain in [s.QQ,s.GF(2),s.QQ.frac_field(a,unused,order='grevlex'),s.ZZ.poly_ring(a,unused,order='grlex'),s.QQ.algebraic_field(s.I,alias='alpha')]:
    p=s.Poly((x*x+1)**2*(x+1),x,domain=domain)
    for op in [s.factor_list,s.sqf_list]:
        c,fs=op(p,polys=True)
        product=s.Poly(c,x,domain=domain)
        for f,k in fs:
            restored=pickle.loads(pickle.dumps(f))
            assert restored.domain==domain and restored.gens==(x,)
            product*=restored**k
        assert s.Poly(s.expand(product.as_expr()-p.as_expr()),x,domain=domain).is_zero
    part=s.sqf_part(p,polys=True)
    assert part.domain==domain and p.rem(part).is_zero
    assert part.gcd(part.diff()).degree()==0
assert s.sqf_part((x*x+1)**2,x,modulus=2)==x+1
'''
    subprocess.run([sys.executable, '-c', program], check=True, timeout=30)
