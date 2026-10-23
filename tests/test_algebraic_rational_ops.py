"""Prospective AN rational-operation tests; activate after the current lot.

Native operations must prove polynomial identities independently of SymPy.
The oracle additionally checks canonical mathematical results and return types.
"""
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm

FIELDS = ['sqrt(2)', 'I', '2**(1/3)', 'sqrt(2)+sqrt(3)',
          '2**(1/3)+sqrt(2)', '2**(1/3)+I*sqrt(3)',
          'sqrt(2)+sqrt(3)+sqrt(5)']


def pair(field):
    x = sm.Symbol('x', real=True)
    sx = sp.Symbol('x', real=True)
    return x, sx, sm.sympify(field), sp.sympify(field)


def as_oracle(expr, sx):
    return sp.sympify(str(expr), locals={'x': sx})


def equal_polynomial(actual, expected, x, sx, alpha):
    assert sp.Poly(as_oracle(actual, sx)-expected, sx, extension=alpha).is_zero


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('polys', [False, True])
def test_bezout_coprime_and_nontrivial_gcd(field, polys):
    x, sx, a, sa = pair(field)
    for common in [sm.Integer(1), x-a]:
        f=common*(x*x+a*x+1)
        g=common*(x+2*a)
        sf=as_oracle(f,sx); sg=as_oracle(g,sx)
        F=sm.Poly(f,x,extension=a) if polys else f
        G=sm.Poly(g,x,extension=a) if polys else g
        s,t,h=sm.gcdex(F,G,x,extension=a)
        os,ot,oh=sp.gcdex(sf,sg,sx,extension=sa)
        if polys:
            assert all(isinstance(p,sm.Poly) and p.domain==F.domain and p.gens==F.gens for p in [s,t,h])
            assert (s*F+t*G)==h
        else:
            assert sm.Poly(s*f+t*g-h,x,extension=a).is_zero
        for actual,expected in zip((s,t,h),(os,ot,oh)):
            equal_polynomial(actual.as_expr() if polys else actual,expected,x,sx,sa)
        ss,hh=sm.half_gcdex(F,G,x,extension=a)
        assert sm.Poly((ss.as_expr() if polys else ss)*f-(hh.as_expr() if polys else hh),x,extension=a).rem(sm.Poly(g,x,extension=a)).is_zero


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('polys', [False, True])
def test_inverse_modulus_certificate_and_oracle(field,polys):
    x,sx,a,sa=pair(field)
    f=x+a;g=x*x-a*a-1
    F=sm.Poly(f,x,extension=a) if polys else f
    G=sm.Poly(g,x,extension=a) if polys else g
    inverse=sm.invert(F,G,x,extension=a)
    expected=sp.invert(sx+sa,sx*sx-sa*sa-1,sx,extension=sa)
    if polys:
        assert isinstance(inverse,sm.Poly) and inverse.domain==F.domain
        assert (F*inverse).rem(G).as_expr()==1
        inverse=inverse.as_expr()
    else:
        assert sm.Poly(f*inverse-1,x,extension=a).rem(sm.Poly(g,x,extension=a)).is_zero
    equal_polynomial(inverse,expected,x,sx,sa)


@pytest.mark.parametrize('field',FIELDS)
def test_cancel_common_factor_and_tuple_scaling(field):
    x,sx,a,sa=pair(field)
    for unit_n,unit_d in [(1,1),(2,2),(a,3),(a+1,a-1)]:
        numerator=unit_n*(x-a)*(x+1)
        denominator=unit_d*(x-a)*(x+2)
        value=sm.cancel(numerator/denominator,x,extension=a)
        n,d=value.as_numer_denom()
        assert sm.Poly(n*denominator-d*numerator,x,extension=a).is_zero
        assert sm.Poly(n,x,extension=a).gcd(sm.Poly(d,x,extension=a)).degree()==0
        expected=sp.cancel(as_oracle(numerator,sx)/as_oracle(denominator,sx),sx,extension=sa)
        en,ed=sp.fraction(expected)
        equal_polynomial(n*sm.sympify(str(ed))-d*sm.sympify(str(en)),0,x,sx,sa)
        coefficient,p,q=sm.cancel((numerator,denominator),x,extension=a)
        assert sm.Poly(coefficient*p*denominator-q*numerator,x,extension=a).is_zero
        assert sm.Poly(p,x,extension=a).gcd(sm.Poly(q,x,extension=a)).degree()==0


@pytest.mark.parametrize('field',FIELDS)
@pytest.mark.parametrize('kind',['split','repeated','improper','irreducible_quadratic'])
def test_apart_crt_identity_and_irreducible_pole_blocks(field,kind):
    x,sx,a,sa=pair(field)
    if kind=='split': f=1/((x-a)*(x+a+1))
    elif kind=='repeated': f=(x+a)/((x-a)**3*(x+1)**2)
    elif kind=='improper': f=(x**5+a*x+1)/((x-a)*(x+2))
    else: f=(x+a)/((x-a)*(x*x+1))
    result=sm.apart(f,x,extension=a)
    # Exact cross-multiplication in the native quotient field: this verifies
    # every CRT inverse and every digit of the factor-adic expansion.
    n,d=f.as_numer_denom();rn,rd=result.as_numer_denom()
    assert sm.Poly(rn*d-rd*n,x,extension=a).is_zero
    expected=sp.apart(as_oracle(f,sx),sx,extension=sa)
    delta=as_oracle(result,sx)-expected
    dn,_=sp.fraction(sp.together(delta))
    assert sp.Poly(dn,sx,extension=sa).is_zero
    # An unchanged rational function satisfies the identity too; require
    # each returned proper-fraction term to contain one irreducible factor.
    for term in sp.Add.make_args(as_oracle(result,sx)):
        _,den=sp.fraction(term)
        if den.has(sx):
            assert len(sp.factor_list(den,sx,extension=sa)[1])==1
    assert isinstance(result,sm.Add)==isinstance(expected,sp.Add)


@pytest.mark.parametrize('field',FIELDS[:4])
def test_seeded_bezout_polynomials_have_exact_certificates(field):
    x,sx,a,sa=pair(field)
    rng=random.Random(63917)
    for _ in range(8):
        fc=[rng.randrange(-4,5)+a*rng.randrange(-2,3) for _ in range(4)]
        gc=[rng.randrange(-3,4)+a*rng.randrange(-2,3) for _ in range(3)]
        f=sum((coefficient*x**i for i,coefficient in enumerate(fc)),sm.Integer(0))
        g=sum((coefficient*x**i for i,coefficient in enumerate(gc)),sm.Integer(0))
        s,t,h=sm.gcdex(f,g,x,extension=a)
        assert sm.Poly(s*f+t*g-h,x,extension=a).is_zero
        H=sm.Poly(h,x,extension=a)
        assert H.is_monic and sm.Poly(f,x,extension=a).rem(H).is_zero and sm.Poly(g,x,extension=a).rem(H).is_zero
        expected=sp.gcd(as_oracle(f,sx),as_oracle(g,sx),sx,extension=sa)
        equal_polynomial(h,expected,x,sx,sa)


def test_noninvertible_zero_and_unsupported_are_explicit():
    x,y=sm.symbols('x y');a=sm.sqrt(2)
    from symox.polys.polyerrors import NotInvertible
    with pytest.raises(NotInvertible):
        sm.invert((x-a)*(x+1),(x-a)*(x+2),x,extension=a)
    with pytest.raises(NotInvertible):
        sm.invert(0,x+a,x,extension=a)
    with pytest.raises(ZeroDivisionError):
        sm.invert(0,0,x,extension=a)
    # Supported since the full=True port; SymPy 1.14 output.
    assert str(sm.apart(1/(x*x-2),x,extension=a,full=True))=='RootSum(_w**2 - 2, Lambda(_a, _a/(-_a + x)))/4'
    for call in [lambda:sm.gcdex(x+y,x-y,x,y,extension=a),
                 lambda:sm.apart(1/(y*x*x-2),x,extension=a),
                 lambda:sm.invert(x+a,x*x-3,x,extension=a,modulus=5)]:
        with pytest.raises((NotImplementedError,ValueError,TypeError)):
            call()


def test_native_only_process_with_bezout_crt_and_cancel():
    script=r'''
import sys,importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy runtime import')
sys.meta_path.insert(0,Block())
import symox as s
x=s.Symbol('x')
for a in [s.sqrt(2),s.I,2**s.Rational(1,3),2**s.Rational(1,3)+s.sqrt(2),s.sqrt(2)+s.sqrt(3)+s.sqrt(5)]:
    f=x+a;g=x*x-a*a-1
    u,v,h=s.gcdex(f,g,x,extension=a)
    assert s.Poly(u*f+v*g-h,x,extension=a).is_zero
    inv=s.invert(f,g,x,extension=a)
    assert s.Poly(f*inv-1,x,extension=a).rem(s.Poly(g,x,extension=a)).is_zero
    result=s.cancel((x*x-a*a)/(x-a),x,extension=a)
    assert s.Poly(result-x-a,x,extension=a).is_zero
    value=1/((x-a)**2*(x+1))
    apart=s.apart(value,x,extension=a)
    n,d=value.as_numer_denom();an,ad=apart.as_numer_denom()
    assert s.Poly(n*ad-an*d,x,extension=a).is_zero
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',script],check=True)


@pytest.mark.parametrize('field',FIELDS)
def test_poly_method_routes_and_return_type_options(field):
    x,sx,a,sa=pair(field)
    p=sm.Poly(x+a,x,extension=a)
    q=sm.Poly(x*x-a*a-1,x,extension=a)
    for auto in [False,True]:
        s,t,h=p.gcdex(q,auto=auto)
        assert s*p+t*q==h
        assert (p*p.invert(q,auto=auto)).rem(q).is_one
        u,v=p.half_gcdex(q,auto=auto)
        assert (u*p-v).rem(q).is_zero
    common=sm.Poly(x-a,x,extension=a)
    c,n,d=(p*common).cancel(q*common)
    assert c==1 and n*(q*common)==d*(p*common)
    assert (p*common).cancel(q*common,include=True)==(n,d)
    s,t,h=sm.gcdex(p,q,polys=False)
    assert isinstance(s,sm.Expr) and isinstance(t,sm.Expr) and isinstance(h,sm.Expr)
    s,t,h=sm.gcdex(p.as_expr(),q.as_expr(),x,extension=a,polys=True)
    assert all(isinstance(v,sm.Poly) for v in [s,t,h])
    assert s*p+t*q==h


@pytest.mark.parametrize('field',FIELDS)
def test_zeros_constants_and_nonprimitive_domain_metadata(field):
    x,sx,a,sa=pair(field)
    zero=sm.Poly(0,x,extension=a)
    one=sm.Poly(1,x,extension=a)
    nonzero=sm.Poly(2*x+a,x,extension=a)
    u,v,h=sm.gcdex(zero,nonzero)
    assert u*zero+v*nonzero==h and h.is_monic
    u,h=sm.half_gcdex(nonzero,zero)
    assert u*nonzero==h
    for call in [lambda:sm.gcdex(zero,zero),lambda:sm.gcdex(nonzero,zero),lambda:zero.invert(zero)]:
        with pytest.raises(ZeroDivisionError):call()
    assert nonzero.invert(one).is_zero
    assert sm.apart(sm.Integer(3),x,extension=a)==3
    assert sm.cancel(sm.Integer(3),x,extension=a)==3
    c,n,d=sm.cancel((zero,nonzero))
    assert c==1 and n.is_zero and d.is_one
    assert n.domain==zero.domain and d.domain==zero.domain


@pytest.mark.parametrize('domain',[None,'QQ','QQ(a)','QQ<sqrt(2)>'])
def test_coeff_monomial_expression_tuple_and_domains(domain):
    x=sm.Symbol('x',positive=True);y=sm.Symbol('y',real=True);a=sm.Symbol('a')
    coefficient=sm.sqrt(2) if domain=='QQ<sqrt(2)>' else a if domain=='QQ(a)' else sm.Rational(2,3) if domain=='QQ' else sm.Integer(2)
    options={} if domain is None else {'domain':domain}
    p=sm.Poly(coefficient*x*x*y+3*x+y+5,x,y,**options)
    assert p.coeff_monomial(x*x*y)==coefficient
    assert p.coeff_monomial((2,1))==coefficient
    assert p.coeff_monomial(x)==3
    assert p.coeff_monomial(1)==5
    assert p.coeff_monomial((0,2))==0
    for monomial in [2*x,x+y,0,1/x,sm.Symbol('z'),(2,),(-1,0),(1,2,3),(sm.Rational(1,2),0),sm.Float(1)]:
        with pytest.raises((ValueError,TypeError,NotImplementedError)):
            p.coeff_monomial(monomial)


def test_independent_fields_unify_and_incompatible_domains_raise():
    x=sm.Symbol('x');a=2**sm.Rational(1,3);b=3**sm.Rational(1,3)
    sx=sp.Symbol('x');sa=2**sp.Rational(1,3);sb=3**sp.Rational(1,3)
    p=sm.Poly(x+a,x,extension=a);q=sm.Poly(x+b,x,extension=b)
    sp_p=sp.Poly(sx+sa,sx,extension=sa);sp_q=sp.Poly(sx+sb,sx,extension=sb)
    # Like the reference, independent number fields unify to the field of
    # a primitive element.
    calls=[lambda m,f,g:f+g,lambda m,f,g:f*g,lambda m,f,g:f.div(g),lambda m,f,g:f.gcdex(g),lambda m,f,g:f.half_gcdex(g),
           lambda m,f,g:f.invert(g),lambda m,f,g:f.cancel(g),lambda m,f,g:m.gcdex(f,g),lambda m,f,g:m.invert(f,g)]
    for call in calls:
        actual=call(sm,p,q);expected=call(sp,sp_p,sp_q)
        if not isinstance(actual,tuple):actual,expected=(actual,),(expected,)
        assert len(actual)==len(expected)
        for u,v in zip(actual,expected):
            if isinstance(v,sp.Poly):
                assert str(u.domain)==str(v.domain)
                assert sp.expand(sp.sympify(str(u.as_expr()))-v.as_expr())==0
            else:
                assert sp.sympify(str(u))==v
    # Gap: the reference unifies these to QQ<2**(1/3)>(a); the native
    # backend has no algebraic rational-function field yet and says so.
    for domain in ['QQ(a)']:
        r=sm.Poly(x+sm.Symbol('a'),x,domain=domain)
        with pytest.raises((ValueError,TypeError,NotImplementedError)):p.gcdex(r)
        with pytest.raises((ValueError,TypeError,NotImplementedError)):sm.gcdex(p,r)
    r=sm.Poly(x+1,x,modulus=5)
    with pytest.raises((ValueError,TypeError,NotImplementedError)):sm.invert(p,r)
    exception=sm.NotInvertible('zero divisor')
    import pickle
    assert type(pickle.loads(pickle.dumps(exception))) is type(exception)
