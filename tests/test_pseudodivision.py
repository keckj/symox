"""Typed pseudo-division over exact rings, with independent certificates."""
import pickle
import subprocess
import sys

import pytest
import symox as sm
from test_global_poly_division import compare, selected

OPS = ['pdiv','prem','pquo','pexquo']


def leading_in_first(p):
    if p.is_zero:
        return sm.Integer(0)
    degree=int(p.degree())
    terms=[]
    for powers,coef in p.as_dict().items():
        if powers[0] == degree:
            value=coef
            for generator,power in zip(p.gens[1:],powers[1:]):
                value *= generator**power
            terms.append(value)
    return sm.Add(*terms)


def certificate(p,g,q,r):
    exponent=0 if p.is_zero else max(int(p.degree())-int(g.degree())+1,0)
    residue=sm.expand(leading_in_first(g)**exponent*p.as_expr()-q.as_expr()*g.as_expr()-r.as_expr())
    assert sm.Poly(residue,*p.gens,domain=p.domain).is_zero
    assert r.is_zero or int(r.degree()) < int(g.degree())
    assert q.domain == r.domain == p.domain
    assert all(c in p.domain for v in [q,r] for c in v.coeffs())


@pytest.mark.parametrize('kind',['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring'])
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
@pytest.mark.parametrize('case',range(5))
def test_methods_and_globals_match_oracle_and_stay_in_coefficient_ring(kind,order,case):
    import sympy as sy
    x,y=sm.symbols('x y')
    k,ok,a,b=selected(kind,order)
    c=a/(b+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind=='AN' else sm.Integer(2)
    divisor=c*(y+1)*x+y+2
    f,g=[(x**3+(c+1)*x+y,divisor), ((x+y+1)*divisor,divisor),
         (sm.Integer(0),divisor), (x+y,x*x+c), (x*x+y*x+1,c)][case]
    p,d=sm.Poly(f,x,y,domain=k),sm.Poly(g,x,y,domain=k)
    op,od=sy.Poly(sm.to_sympy(f),sm.to_sympy(x),sm.to_sympy(y),domain=ok),sy.Poly(sm.to_sympy(g),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    q,r=p.pdiv(d)
    compare((q,r),op.pdiv(od))
    certificate(p,d,q,r)
    for result in [q,r]:
        restored=pickle.loads(pickle.dumps(result))
        assert restored.domain==k and restored.as_expr()==result.as_expr()
        if kind.endswith('RF') or kind.endswith('ring'):
            assert result.domain.symbols==(a,b) and result.domain.order==order
    for name in OPS:
        try:
            expected=getattr(op,name)(od)
        except sy.polys.polyerrors.ExactQuotientFailed:
            with pytest.raises(sm.ExactQuotientFailed):getattr(p,name)(d)
            with pytest.raises(sm.ExactQuotientFailed):getattr(sm,name)(p,d)
            with pytest.raises(sm.ExactQuotientFailed):getattr(sm,name)(f,g,x,y,domain=k,polys=True)
        else:
            compare(getattr(p,name)(d),expected)
            compare(getattr(sm,name)(p,d),expected)
            compare(getattr(sm,name)(f,g,x,y,domain=k,polys=True),expected)


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('polys',[None,False,True])
@pytest.mark.parametrize('input_kind',['expr','poly','mixed'])
def test_global_output_flags_and_joint_inference(op,polys,input_kind):
    import sympy as sy
    x=sm.Symbol('x');sx=sm.to_sympy(x)
    f,g=(x+1)*(2*x+1),2*x+1
    sf,sg=sm.to_sympy(f),sm.to_sympy(g)
    if input_kind!='expr':f,sf=sm.Poly(f,x),sy.Poly(sf,sx)
    if input_kind=='poly':g,sg=sm.Poly(g,x),sy.Poly(sg,sx)
    compare(getattr(sm,op)(f,g,polys=polys),getattr(sy,op)(sf,sg,polys=polys))


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('kind',['GF','AN','QQRF','ZZring'])
@pytest.mark.parametrize('option',['field','extension','gens'])
def test_common_options_preserve_or_convert_domains_explicitly(op,kind,option):
    import sympy as sy
    x,y,z=sm.symbols('x y z')
    k,ok,a,b=selected(kind,'grevlex')
    p,d=sm.Poly((x+1)*(2*x+y),x,y,domain=k),sm.Poly(2*x+y,x,y,domain=k)
    sp=sy.Poly(sm.to_sympy(p.as_expr()),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    sd=sy.Poly(sm.to_sympy(d.as_expr()),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    kwargs={'field':True} if option=='field' else {'extension':True} if option=='extension' else {'gens':(z,)}
    oracle_kwargs={'field':True} if option=='field' else {'extension':True} if option=='extension' else {'gens':(sm.to_sympy(z),)}
    compare(getattr(sm,op)(p,d,**kwargs),getattr(sy,op)(sp,sd,**oracle_kwargs))


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('value',[True,False,None,2,'bad'])
def test_auto_is_forbidden_with_real_flag_error(op,value):
    from symox.polys.polyerrors import FlagError,OptionError
    x=sm.Symbol('x')
    with pytest.raises(FlagError) as failure:getattr(sm,op)(x*x,2*x,x,auto=value)
    assert sm.FlagError is FlagError and issubclass(FlagError,OptionError)
    assert type(pickle.loads(pickle.dumps(failure.value))) is FlagError
    with pytest.raises(TypeError):getattr(sm.Poly(x*x,x),op)(sm.Poly(2*x,x),auto=value)


@pytest.mark.parametrize('op',OPS)
def test_unknown_options_incompatible_fields_and_zero_divisors(op):
    x=sm.Symbol('x')
    for kwargs in [{'bogus':None},{'domain':sm.ZZ,'field':True},{'domain':sm.QQ,'extension':True},{'modulus':5,'extension':True}]:
        with pytest.raises(sm.OptionError):getattr(sm,op)(x+1,x,x,**kwargs)
    for k in [sm.ZZ,sm.QQ,sm.GF(5),sm.QQ.algebraic_field(sm.sqrt(2))]:
        p,d=sm.Poly(x*x+1,x,domain=k),sm.Poly(0,x,domain=k)
        with pytest.raises(ZeroDivisionError):getattr(p,op)(d)
        with pytest.raises(ZeroDivisionError):getattr(sm,op)(p,d)
    with pytest.raises(sm.UnificationFailed):getattr(sm,op)(sm.Poly(x,x,domain=sm.GF(5)),sm.Poly(x,x,domain=sm.GF(7)))
    with pytest.raises(ZeroDivisionError):getattr(sm,op)(x*x+1,5*x,x,modulus=5)


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('prime',[2,5,7])
def test_modular_and_typed_finite_coercion(op,prime):
    import sympy as sy
    x=sm.Symbol('x');sx=sm.to_sympy(x)
    p,d=sm.Poly((3*x+1)*(x+1),x,domain=sm.GF(5)),sm.Poly(x+1,x,domain=sm.GF(5))
    sp,sd=sy.Poly((3*sx+1)*(sx+1),sx,domain=sy.GF(5)),sy.Poly(sx+1,sx,domain=sy.GF(5))
    compare(getattr(sm,op)(p,d,modulus=prime),getattr(sy,op)(sp,sd,modulus=prime))


@pytest.mark.parametrize('op',OPS)
def test_algebraic_alias_and_parameter_identity_roundtrip(op):
    x=sm.Symbol('x')
    a,b=sm.Dummy('a',positive=True),sm.Symbol('b',real=True)
    for k in [sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha'),sm.ZZ.poly_ring(a,b,order='grlex'),sm.QQ.frac_field(a,b,order='grevlex')]:
        p,d=sm.Poly(x*x-1,x,domain=k),sm.Poly(x-1,x,domain=k)
        for result in [getattr(p,op)(d),getattr(sm,op)(p,d,extension=True)]:
            for v in result if op=='pdiv' else (result,):
                assert v.domain==k
                assert pickle.loads(pickle.dumps(v)).domain==k
                if k.is_AlgebraicField:assert v.domain.ext.alias==sm.Symbol('alpha')
                else:assert v.domain.symbols==(a,b) and v.domain.order==k.order


@pytest.mark.parametrize('op',OPS)
@pytest.mark.parametrize('coefficient',['radical','function','imaginary','float'])
def test_historical_opaque_coefficients_still_compute(op,coefficient):
    import sympy as sy
    x,a=sm.symbols('x a')
    c={'radical':sm.sqrt(2),'function':sm.sin(a),'imaginary':sm.I,'float':0.5*a}[coefficient]
    f,g=(c*x+1)*(x+1),c*x+1
    compare(getattr(sm,op)(f,g,x),getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x)))


def test_sparse_extreme_degrees_are_safe_and_growth_is_rejected():
    x,y=sm.symbols('x y')
    largest=2**32-1
    p=sm.Poly(x**largest,x)
    q,r=p.pdiv(sm.Poly(1,x))
    assert q==p and r.is_zero and int(q.degree())==largest
    with pytest.raises(ValueError,match='resource limit'):
        sm.Poly(x**1_000_000_000,x).pdiv(sm.Poly(2*x,x))
    with pytest.raises(ValueError,match='exponent overflow'):
        sm.Poly(x*x*y**largest,x,y).pdiv(sm.Poly(x*y**largest+1,x,y))


def test_pseudo_division_without_sympy_runtime():
    code='''
import importlib.abc,sys,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('blocked')
sys.meta_path.insert(0,Block())
import symox as s
x,y=s.symbols('x y');a=s.Dummy('a',positive=True)
for k in [s.ZZ,s.QQ,s.GF(5),s.QQ.algebraic_field(s.sqrt(2)),s.QQ.frac_field(a),s.ZZ.poly_ring(a)]:
 p=s.Poly(x*x+y*x+1,x,y,domain=k);g=s.Poly(2*x+y,x,y,domain=k)
 q,r=s.pdiv(p,g);assert q.domain==r.domain==k
 residue=4*p.as_expr()-q.as_expr()*g.as_expr()-r.as_expr()
 assert s.Poly(residue,x,y,domain=k).is_zero
 assert s.pquo(p,g)==q and s.prem(p,g)==r
 for v in [q,r]:assert pickle.loads(pickle.dumps(v)).domain==k
 try:s.pexquo(p,g)
 except s.ExactQuotientFailed:assert not r.is_zero
 else:assert r.is_zero
try:s.pdiv(x,x,x,auto=None)
except s.FlagError as e:assert type(pickle.loads(pickle.dumps(e))) is s.FlagError
else:raise AssertionError('auto ignored')
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
