"""Pure polynomial operations preserve coefficients while ignoring generator names."""
import pickle
import random
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from symox._purepoly import PurePoly,install

install()


def assert_same(actual,expected):
    if isinstance(expected,sp.PurePoly):
        assert isinstance(actual,PurePoly)
        assert isinstance(actual,sx.Poly)
        assert str(actual.domain)==str(expected.domain)
        assert tuple(g.to_sympy() for g in actual.gens)==expected.gens
        assert sp.cancel(actual.as_expr().to_sympy()-expected.as_expr())==0
    elif isinstance(expected,(list,tuple)):
        assert type(actual) is type(expected)
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):assert_same(a,b)
    else:
        assert (actual.to_sympy() if hasattr(actual,'to_sympy') else actual)==expected


@pytest.mark.parametrize('method',['add','sub','mul','div','rem','quo','exquo',
    'pdiv','prem','pquo','pexquo','gcd','lcm','cofactors','gcdex',
    'half_gcdex','resultant','subresultants','compose'])
def test_binary_algorithms_align_generator_slots(method):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    f=PurePoly(x*x-1,x,domain='QQ');g=PurePoly(y-1,y,domain='QQ')
    F=sp.PurePoly(X*X-1,X,domain='QQ');G=sp.PurePoly(Y-1,Y,domain='QQ')
    assert_same(getattr(f,method)(g),getattr(F,method)(G))


@pytest.mark.parametrize('domain',['ZZ','QQ','GF(7)','QQ<sqrt(2)>'])
def test_random_exact_coefficient_arithmetic_and_alpha_equality(domain):
    x,y,z,w=sx.symbols('x y z w');X,Y,Z,W=sp.symbols('x y z w')
    rng=random.Random(72)
    for _ in range(30):
        coefficients=[rng.randrange(-5,6) for _ in range(6)]
        a,b,c,d,e,f=coefficients
        left=a*x*x+b*x*y+c*y*y+d*x+e*y+f
        right=b*z*z+c*z*w+d*w*w+e*z+f*w+a
        L=a*X*X+b*X*Y+c*Y*Y+d*X+e*Y+f
        R=b*Z*Z+c*Z*W+d*W*W+e*Z+f*W+a
        p=PurePoly(left,x,y,domain=domain);q=PurePoly(right,z,w,domain=domain)
        P=sp.PurePoly(L,X,Y,domain=domain);Q=sp.PurePoly(R,Z,W,domain=domain)
        assert_same(p+q,P+Q);assert_same(p-q,P-Q);assert_same(p*q,P*Q)
        twin=PurePoly(left,x,y,domain=domain)._with_gens(z,w)
        assert p==twin and hash(p)==hash(twin)
        assert (p==q)==(P==Q)


def test_constructor_factories_and_result_type():
    x=sx.Symbol('x');X=sp.Symbol('x')
    for method,args in [('from_list',([1,0,-1],x)),
                        ('from_dict',({(2,):1,(0,):-1},x)),
                        ('from_expr',(x*x-1,x)),
                        ('from_poly',(sx.Poly(x*x-1,x),))]:
        ref_args=([1,0,-1],X) if method=='from_list' else (
            ({(2,):1,(0,):-1},X) if method=='from_dict' else (
            (X*X-1,X) if method=='from_expr' else (sp.Poly(X*X-1,X),)))
        assert_same(getattr(PurePoly,method)(*args),getattr(sp.PurePoly,method)(*ref_args))
    assert_same(PurePoly(rep={2:1,0:2},gens=(x,)),sp.PurePoly(rep={2:1,0:2},gens=(X,)))


def test_differentiation_evaluation_containers_and_pickle():
    x,y,a=sx.symbols('x y a');X,Y,A=sp.symbols('x y a')
    p=PurePoly(a*x*x+x*y+y,x,y);P=sp.PurePoly(A*X*X+X*Y+Y,X,Y)
    assert p.free_symbols=={a}
    assert p.func is PurePoly
    for actual,expected in [(p.diff(x),P.diff(X)),(p.integrate(x),P.integrate(X)),
            (p.eval(x,2),P.eval(X,2)),(p.eval({x:2,y:3}),P.eval({X:2,Y:3})),
            (p.factor_list(),P.factor_list()),(p.primitive(),P.primitive()),
            (p.one,P.one),(p.zero,P.zero),(p.as_expr(y,x),P.as_expr(Y,X)),
            (p.as_expr({x:y}),P.as_expr({X:Y}))]:assert_same(actual,expected)
    result=pickle.loads(pickle.dumps(p))
    assert type(result) is PurePoly and result==p
    assert_same(result,P)


def test_different_generator_counts_and_functional_generators():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    p=PurePoly(x+1,x);q=PurePoly(x+y,x,y)
    assert p!=q
    with pytest.raises(sx.UnificationFailed):p+q
    assert PurePoly(x)==PurePoly(sx.exp(x))
    assert sp.PurePoly(X)==sp.PurePoly(sp.exp(X))
    assert p==sx.Poly(y+1,y)
    assert p==PurePoly(y+1,y,domain='QQ')
    assert_same(p+sx.Poly(y+2,y),sp.PurePoly(X+1,X)+sp.Poly(Y+2,Y))


def test_exact_native_arithmetic_does_not_import_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
import symox as s
from symox._purepoly import PurePoly,install
install()
x,y=s.symbols('x y')
p=PurePoly(x*x-1,x);q=PurePoly(y-1,y)
assert p.gcd(q)==q
assert p.div(q)[1].is_zero
assert isinstance(p.diff(),PurePoly)
assert p.all_coeffs()==[1,0,-1]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
