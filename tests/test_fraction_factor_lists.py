"""Fractional lists preserve both sides, exact domains and positive powers."""
import os
import pickle
import subprocess
import sys

import pytest
import symox as sm

from test_global_poly_division import selected
from test_typed_factorization import compare_result, expression_of

OPS = ['factor_list', 'sqf_list']


def expression_product(factors):
    value = sm.Integer(1)
    for factor, exponent in factors:
        assert isinstance(exponent, int) and exponent > 0
        value *= expression_of(factor)**exponent
    return value


def certify(expression, result, gens=(), domain=None, known_parts=None):
    coefficient, numerator, denominator = result
    assert isinstance(result, tuple) and len(result) == 3
    n, d = known_parts if known_parts is not None else expression.as_numer_denom()
    difference = sm.expand(coefficient*expression_product(numerator)*d-n*expression_product(denominator))
    if domain is None:
        assert difference == 0
    else:
        # Certify the cross product in the selected ring, rather than clearing
        # finite coefficients or checking only approximate sample evaluations.
        assert sm.Poly(difference, *gens, domain=domain.get_field()).is_zero


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('polys', [None, False, True])
@pytest.mark.parametrize('case', range(12))
def test_natural_fraction_oracle_and_cross_product(op, polys, case):
    import sympy as sy
    x, y, a = sm.symbols('x y a')
    f = [(x*x-1)/(x-1), 1/x+1/x**2, 1/(x+1)+1/(x*x+2*x+1),
         x/y, (x+1)/(y+sm.Rational(1,2)), x/a,
         sm.I*x/(x+1), sm.sqrt(2)*x/y,
         3*(x-1)**2*(x+2)**3/((x-1)**4*(x+3)**2),
         (x*y+x)/(x*x*y-y), sm.pi*x/y, sm.exp(-x)*(x+1)][case]
    opts = {'frac':True, 'polys':polys}
    actual = getattr(sm, op)(f, **opts)
    expected = getattr(sy, op)(sm.to_sympy(f), **opts)
    compare_result(actual, expected)
    certify(f, actual)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring'])
@pytest.mark.parametrize('order', ['lex','grlex','grevlex'])
@pytest.mark.parametrize('polys', [False, True])
def test_selected_domain_metadata_units_and_certificates(op, kind, order, polys):
    import sympy as sy
    x = sm.Dummy('x', real=True)
    k, ok, a, unused = selected(kind, order)
    c = a/(unused+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    # Expand both sides before interoperability conversion: otherwise SymPy's
    # constructor cancels identical literal factors before reaching its API.
    n = sm.expand((c*x+1)**2*(x-3))
    d = sm.expand(3*(c*x+1)*(x+2)**3)
    f = n/d
    opts = {'frac':True,'domain':k,'polys':polys}
    actual = getattr(sm,op)(f,x,**opts)
    expected = getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(x),frac=True,domain=ok,polys=polys)
    if kind == 'GF':
        # SymPy returns the unreduced rational quotient of its two units.
        # Native lists reduce that quotient in GF(p), including zero checks.
        assert actual[0] in k
        assert k(actual[0]) == k(sm.from_sympy(expected[0]))
        compare_result(actual[1:],expected[1:])
    elif kind == 'QQring' and op == 'sqf_list':
        # SymPy 1.14 scales these two square-free factors by 1/18 and
        # 1/12 without updating their contents. Its triple does not reconstruct
        # the original fraction. Retain our independently certified factors.
        sn, sd = sm.to_sympy(f).as_numer_denom()
        def oracle_product(fs):
            return sy.prod((p.as_expr() if isinstance(p,sy.Poly) else p)**m for p,m in fs)
        assert sy.cancel(expected[0]*oracle_product(expected[1])*sd-sn*oracle_product(expected[2])) != 0
        assert actual[0] == sm.Rational(1,3)
        def expected_factor(value):
            value=sm.to_sympy(value)
            return sy.Poly(value,sm.to_sympy(x),domain=ok) if polys else value
        compare_result(actual[1:],
                       ([(expected_factor(x-3),1),(expected_factor(a*x+1),2)],
                        [(expected_factor(a*x+1),1),(expected_factor(x+2),3)]))
    else:
        assert sy.cancel(sm.to_sympy(actual[0])-expected[0]) == 0
        compare_result(actual[1:],expected[1:])
    for factor, _ in actual[1]+actual[2]:
        if polys:
            assert factor.domain == k and factor.gens == (x,)
            restored = pickle.loads(pickle.dumps(factor))
            assert restored.domain == k and restored.gens == (x,)
            assert all(v in k for v in restored.coeffs())
    certify(f,actual,(x,),k)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('polys', [None, False, True])
@pytest.mark.parametrize('constant', [0, sm.Rational(2,3), sm.I, 2+sm.I, sm.sqrt(2), sm.pi])
def test_constants(op, polys, constant):
    import sympy as sy
    actual = getattr(sm,op)(constant,frac=True,polys=polys)
    compare_result(actual,getattr(sy,op)(sm.to_sympy(constant),frac=True,polys=polys))
    compare_result(getattr(sm,op)(constant,gens=(),frac=True,polys=polys),getattr(sy,op)(sm.to_sympy(constant),gens=(),frac=True,polys=polys))


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('polys', [None, False, True])
def test_poly_is_a_formal_numerator_including_reciprocal_generator(op,polys):
    import sympy as sy
    x = sm.Symbol('x')
    p = sm.Poly((1/x+1)**2,1/x,domain=sm.QQ)
    oracle = sy.Poly(sm.to_sympy(p.as_expr()),1/sm.to_sympy(x),domain=sy.QQ)
    actual = getattr(sm,op)(p,frac=True,polys=polys)
    compare_result(actual,getattr(sy,op)(oracle,frac=True,polys=polys))
    assert actual[2] == []


@pytest.mark.parametrize('op',OPS)
def test_common_factors_are_not_cancelled_and_squarefree_groups_are_certified(op):
    x=sm.Symbol('x')
    f=(x*x-1)/(x-1)
    result=getattr(sm,op)(f,x,frac=True,polys=True)
    assert sm.expand(expression_product(result[1])-(x*x-1)) == 0
    assert expression_product(result[2]) == x-1
    certify(f,result)
    # SymPy's component-wise sqf_list multiplies the components of equal
    # multiplicity, so its multiplicity-one group can retain a repeated
    # factor: `x**3 + 4*x**2 + 5*x + 2` is `(x + 1)**2*(x + 2)`. symox
    # follows it (wave 45, lot algfactor), as `sqf` already did.
    import sympy as sy
    f=(x+1)*sm.expand((x+1)*(x+2))/(x-3)
    result=getattr(sm,op)(f,x,frac=True,polys=True)
    oracle=getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(x),frac=True,polys=True)
    assert [(str(a.as_expr()),k) for a,k in result[1]+result[2]] == [(str(a.as_expr()),k) for a,k in oracle[1]+oracle[2]]
    certify(f,result)


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('options',[{'modulus':2},{'domain':sm.GF(2)}])
def test_zero_denominator_after_coercion(op,options):
    x=sm.Symbol('x')
    for f in [1/(2*x+2),1/(2*x*x+4),sm.Rational(1,2)/(x+1)]:
        with pytest.raises(ZeroDivisionError,match='denominator.*zero'):
            getattr(sm,op)(f,x,frac=True,**options)


@pytest.mark.parametrize('op',OPS)
def test_options_and_bounded_split(op):
    x=sm.Symbol('x')
    for opts in [{'unknown':True},{'frac':'invalid'},{'domain':sm.QQ,'modulus':2}]:
        with pytest.raises(sm.OptionError):
            getattr(sm,op)(x,**({'frac':True}|opts))
    with pytest.raises(sm.FlagError):
        getattr(sm,op)(x,frac=True,auto=False)
    for f in [sm.Pow(x+1,-10**9,evaluate=False),sm.Pow(2,-10**9,evaluate=False)]:
        with pytest.raises(ValueError,match='resource limit'):
            getattr(sm,op)(f,x,frac=True)
    with pytest.raises(NotImplementedError,match='noninteger or symbolic'):
        getattr(sm,op)(1/sm.sqrt(x),frac=True)


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('kind',['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring'])
def test_poly_input_preserves_native_descriptor_and_default_result_type(op,kind):
    x=sm.Dummy('x',positive=True)
    k,_,a,unused=selected(kind,'grevlex')
    if kind=='AN':k=sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha')
    p=sm.Poly((x+1)**2*(x-2),x,domain=k)
    result=getattr(sm,op)(p,frac=True)
    assert result[2] == []
    assert all(isinstance(f,sm.Poly) and f.domain==k and f.gens==(x,) for f,_ in result[1])
    certify(p.as_expr(),result,(x,),k)
    expression_result=getattr(sm,op)(p,frac=True,polys=False)
    assert all(not isinstance(f,sm.Poly) for f,_ in expression_result[1])


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('domain_option',['domain','modulus','field','extension'])
def test_selected_options_and_inferred_side_generators(op,domain_option):
    import sympy as sy
    x,y=sm.symbols('x y')
    options={'domain':{'domain':sm.QQ},'modulus':{'modulus':5},'field':{'field':True},'extension':{'extension':sm.sqrt(2)}}[domain_option]
    oracle_options={'domain':{'domain':sy.QQ},'modulus':{'modulus':5},'field':{'field':True},'extension':{'extension':sy.sqrt(2)}}[domain_option]
    f=(x+1)/(y+2)
    result=getattr(sm,op)(f,frac=True,polys=True,**options)
    compare_result(result,getattr(sy,op)(sm.to_sympy(f),frac=True,polys=True,**oracle_options))
    assert result[1][0][0].gens==(x,) and result[2][0][0].gens==(y,)


@pytest.mark.parametrize('op',OPS)
def test_extension_inference_uses_both_sides_and_preserves_primitive_embedding(op):
    x,y=sm.symbols('x y')
    f=(x+sm.sqrt(2))/(y+sm.sqrt(3))
    result=getattr(sm,op)(f,frac=True,extension=True,polys=True)
    pn,pd=result[1][0][0],result[2][0][0]
    assert pn.domain==pd.domain
    assert sm.sqrt(2) in pn.domain and sm.sqrt(3) in pn.domain
    assert pn.gens==(x,) and pd.gens==(y,)
    certify(f,result,(x,y),pn.domain)


@pytest.mark.parametrize('op',OPS)
def test_explicit_generator_order_and_parameters(op):
    import sympy as sy
    x,y=sm.symbols('x y');a=sm.Symbol('a',positive=True)
    f=(a+x)/(a+y)
    result=getattr(sm,op)(f,gens=(y,x),frac=True,polys=True)
    expected=getattr(sy,op)(sm.to_sympy(f),gens=(sm.to_sympy(y),sm.to_sympy(x)),frac=True,polys=True)
    compare_result(result,expected)
    for p,_ in result[1]+result[2]:
        assert p.gens==(y,x)
        assert p.domain.symbols==(a,)


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('case',range(8))
def test_nested_rational_arithmetic_has_exact_cross_product(op,case):
    import sympy as sy
    x,y=sm.symbols('x y')
    f=[1/(x/y+1),(1/x+1/y)/(x+y),1/(1/(x+1)+1/(x-1)),
       (1/x-1/y)**2,(x/(x+1)+y/(y+1))**-2,
       (sm.Rational(1,2)/x+sm.Rational(2,3)/y)/(x-1),
       (1/(x+1)+1/(x+2))/(1/(x+3)+1/(x+4)),
       1/(-x-1)+1/(x+2)][case]
    # Independently written common-denominator identities avoid relying on
    # the legacy Expr.as_numer_denom treatment of nested reciprocal sums.
    known_parts=[(y,x+y),(x+y,x*y*(x+y)),(x*x-1,2*x),
                 ((y-x)**2,x*x*y*y),((x+1)**2*(y+1)**2,(2*x*y+x+y)**2),
                 (3*y+4*x,6*x*y*(x-1)),
                 ((2*x+3)*(x+3)*(x+4),(x+1)*(x+2)*(2*x+7)),
                 (sm.Integer(-1),(x+1)*(x+2))][case]
    result=getattr(sm,op)(f,frac=True,polys=True)
    certify(f,result,known_parts=known_parts)
    n=sm.to_sympy(expression_product(result[1]));d=sm.to_sympy(expression_product(result[2]))
    assert sy.cancel(sm.to_sympy(result[0])*n/d-sm.to_sympy(f)) == 0


@pytest.mark.parametrize('op',OPS)
def test_python_conversion_callbacks_and_errors_stay_outside_native_cancellation(op):
    x=sm.Symbol('x')
    calls=[]
    class Converted:
        def _symox_(self):
            # Reentrant native arithmetic is allowed during user conversion.
            calls.append(sm.Poly((x+1)**3,x).degree())
            return (x*x-1)/(x-1)
    result=getattr(sm,op)(Converted(),x,frac=True,polys=True)
    assert calls == [3]
    certify((x*x-1)/(x-1),result)
    class Broken:
        def _symox_(self):
            raise LookupError('conversion sentinel')
    with pytest.raises(LookupError,match='conversion sentinel'):
        getattr(sm,op)(Broken(),x,frac=True)
    assert getattr(sm,op)(1/(x+1),x,frac=True) == (1,[],[(x+1,1)])


def test_no_sympy_runtime():
    program=r'''
import importlib.abc,sys,pickle
class NoSymPy(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0,NoSymPy())
import symox as s
x=s.Dummy('x',real=True); a=s.Symbol('a',positive=True); b=s.Dummy('unused')
for k in [s.ZZ,s.QQ,s.GF(5),s.QQ.frac_field(a,b,order='grevlex'),s.QQ.poly_ring(a,b,order='grlex'),s.QQ.algebraic_field(s.I,alias='alpha')]:
    for op in [s.factor_list,s.sqf_list]:
        f=(x*x-1)/(x-1)
        c,ns,ds=op(f,x,domain=k,frac=True,polys=True)
        n=s.Integer(1);d=s.Integer(1)
        for fs,isnum in [(ns,True),(ds,False)]:
            for p,m in fs:
                p=pickle.loads(pickle.dumps(p));assert p.domain==k and p.gens==(x,)
                if isnum:n*=p.as_expr()**m
                else:d*=p.as_expr()**m
        assert s.Poly(s.expand(c*n*(x-1)-(x*x-1)*d),x,domain=k.get_field()).is_zero
'''
    subprocess.run([sys.executable,'-c',program],check=True,timeout=30,env=os.environ.copy())
