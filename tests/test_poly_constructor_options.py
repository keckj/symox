"""Actual constructor inputs/options, including keyword and factory calls."""
import subprocess
import sys
import pickle

import pytest
import symox as sx
import sympy as sp


def compare(native,reference,kind):
    assert isinstance(native,getattr(sx,kind))
    assert tuple(g.to_sympy() for g in native.gens)==reference.gens
    assert {m:c.to_sympy() for m,c in native.as_dict().items()}==reference.as_dict()
    assert str(native.domain)==str(reference.domain)
    convert=lambda value:value.to_sympy() if hasattr(value,'to_sympy') else value
    assert tuple(map(convert,native.degree_list()))==reference.degree_list()
    assert convert(native.total_degree())==reference.total_degree()
    assert native.length()==reference.length()
    for flag in ['is_zero','is_one','is_ground','is_monomial','is_linear','is_quadratic','is_homogeneous']:
        assert getattr(native,flag)==getattr(reference,flag),flag


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('source',[
    'x**2+2*x+1','x*y+y**2+2','Rational(1,3)*x+2',
])
@pytest.mark.parametrize('form',['scalar','tuple','list','nested','positional'])
def test_generator_keyword_forms_are_real_constructor_arguments(kind,source,form):
    actual,reference=getattr(sx,kind),getattr(sp,kind)
    x,y=sx.symbols('x y');rx,ry=sp.symbols('x y')
    if form=='positional':
        native=actual(sx.sympify(source),x,y);expected=reference(sp.sympify(source),rx,ry)
    else:
        ng,rg={'scalar':(x,rx),'tuple':((x,y),(rx,ry)),
               'list':([x,y],[rx,ry]),'nested':(([x,y],),([rx,ry],))}[form]
        native=actual(rep=sx.sympify(source),gens=ng)
        expected=reference(rep=sp.sympify(source),gens=rg)
    compare(native,expected,kind)


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('data',[
    {2:3,0:5},{(3,):2,(1,):sx.Rational(1,3)},
    [1,0,2],(1,0,2),[],{},[sx.Rational(1,3),2,0],
])
@pytest.mark.parametrize('options',[{}, {'field':True},{'domain':'QQ'}])
def test_dense_and_sparse_input_containers(kind,data,options):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    if isinstance(data,dict):reference={m:sp.sympify(str(c)) for m,c in data.items()}
    else:reference=type(data)(sp.sympify(str(c)) for c in data)
    native=getattr(sx,kind)(data,gens=x,**options)
    expected=getattr(sp,kind)(reference,gens=rx,**options)
    compare(native,expected,kind)


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('factory',['from_dict','from_list','from_expr','from_poly'])
def test_factories_honor_field_and_keyword_generators(kind,factory):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    native,reference=getattr(sx,kind),getattr(sp,kind)
    nrep,rrep={
        'from_dict':({(2,):1,(0,):2},{(2,):1,(0,):2}),
        'from_list':([1,0,2],[1,0,2]),
        'from_expr':(x**2+2,rx**2+2),
        'from_poly':(native(x**2+2,x),reference(rx**2+2,rx)),
    }[factory]
    compare(getattr(native,factory)(nrep,gens=(x,),field=True),
            getattr(reference,factory)(rrep,gens=(rx,),field=True),kind)


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('source',[
    '(x+1)**3','(x+y)**2*(y+1)','x**Rational(3,2)+sqrt(x)',
    'exp(2*x)+exp(x)','(x+1)**2+2*(x+1)',
])
def test_unexpanded_generator_inference_uses_exact_factors(kind,source):
    compare(getattr(sx,kind)(sx.sympify(source),expand=False),
            getattr(sp,kind)(sp.sympify(source),expand=False),kind)


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('options',[
    {'sort':'y>x'},{'sort':['y','x']},{'wrt':'y'},
    {'wrt':'y, z'},{'wrt':['z','x']},{'sort':'y>x','wrt':'z'},
])
def test_generator_sorting_and_priority(kind,options):
    compare(getattr(sx,kind)(sx.sympify('x+y+z'),**options),
            getattr(sp,kind)(sp.sympify('x+y+z'),**options),kind)


@pytest.mark.parametrize('options',[
    {'foo':True},{'field':False,'domain':'ZZ'},
    {'greedy':False,'domain':'ZZ'},{'expand':2},
    {'symmetric':False},{'extension':False},{'sort':object()},
    {'wrt':object()},{'method':3},{'symbols':3},
])
def test_invalid_options_use_reference_exception_class(options):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    try:sp.Poly(rx+1,rx,**options)
    except Exception as expected:
        with pytest.raises(Exception) as actual:sx.Poly(x+1,x,**options)
        assert type(actual.value).__name__==type(expected).__name__
    else:pytest.fail('case must be invalid in reference')


@pytest.mark.parametrize('options',[{'order':'lex'},{'split':False}])
def test_unimplemented_reference_options_are_not_silently_ignored(options):
    x,rx=sx.Symbol('x'),sp.Symbol('x')
    with pytest.raises(NotImplementedError):sp.Poly(rx+1,rx,**options)
    with pytest.raises(NotImplementedError):sx.Poly(x+1,x,**options)


@pytest.mark.parametrize('data',[{(2,):1},[1,0,2],(1,0,2)])
def test_data_requires_generators(data):
    with pytest.raises(Exception) as expected:sp.Poly(data)
    with pytest.raises(Exception) as actual:sx.Poly(data)
    assert type(actual.value).__name__==type(expected.value).__name__


def test_input_nonmutation_and_generator_consumption_once():
    x=sx.Symbol('x');data={(2,):1,(0,):sx.Rational(1,3)};saved=data.copy()
    assert sx.Poly(data,x,field=True).domain==sx.QQ
    assert data==saved
    seen=[]
    def coefficients():
        for c in [1,0,2]:seen.append(c);yield c
    assert sx.Poly(coefficients(),x).as_dict()=={(2,):sx.Integer(1),(0,):sx.Integer(2)}
    assert seen==[1,0,2]


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('generator',['sin(y)','y+1'])
@pytest.mark.parametrize('domain',['algebraic','polynomial','fraction'])
def test_structured_domains_and_composite_generators_preserve_exact_coefficients(kind,generator,domain):
    a,y=sx.symbols('a y');ra,ry=sp.symbols('a y')
    k,rk={
        'algebraic':(sx.QQ.algebraic_field(sx.sqrt(2)),sp.QQ.algebraic_field(sp.sqrt(2))),
        'polynomial':(sx.ZZ.poly_ring(a),sp.ZZ.poly_ring(ra)),
        'fraction':(sx.QQ.frac_field(a),sp.QQ.frac_field(ra)),
    }[domain]
    c,rc=(sx.sqrt(2),sp.sqrt(2)) if domain=='algebraic' else (a,ra)
    g,rg=sx.sympify(generator),sp.sympify(generator)
    p=getattr(sx,kind)([k.convert(c),k.one],g,domain=k)
    reference=getattr(sp,kind)([rk.convert(rc),rk.one],rg,domain=rk)
    compare(p,reference,kind)
    compare(p*p,reference*reference,kind)
    compare(p.diff(),reference.diff(),kind)
    quotient,remainder=(p*p).div(p)
    compare(quotient,reference,kind)
    assert remainder.is_zero


@pytest.mark.parametrize('generator',['y+1','sin(y)','sqrt(y)'])
@pytest.mark.parametrize('protocol',range(pickle.HIGHEST_PROTOCOL+1))
def test_purepoly_interchange_and_pickle_preserve_formal_generator_slots(generator,protocol):
    g=sx.sympify(generator);rg=sp.sympify(generator)
    native=sx.PurePoly([1,0,2],g,domain=sx.QQ)
    reference=sp.PurePoly([1,0,2],rg,domain=sp.QQ)
    assert native.to_sympy().as_dict()==reference.as_dict()
    assert native.to_sympy().gens==reference.gens
    restored=pickle.loads(pickle.dumps(native,protocol=protocol))
    compare(restored,reference,'PurePoly')


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('generator',['2','sqrt(2)','I','pi'])
@pytest.mark.parametrize('protocol',[0,pickle.HIGHEST_PROTOCOL])
def test_numeric_formal_generator_slots_survive_arithmetic_and_pickle(kind,generator,protocol):
    g,rg=sx.sympify(generator),sp.sympify(generator)
    native=getattr(sx,kind)([1,0,2],g)
    reference=getattr(sp,kind)([1,0,2],rg)
    compare(native,reference,kind)
    compare(native*native,reference*reference,kind)
    compare(pickle.loads(pickle.dumps(native,protocol=protocol)),reference,kind)


@pytest.mark.parametrize('kind',['Poly','PurePoly'])
@pytest.mark.parametrize('domain',['polynomial','fraction'])
def test_unifying_structured_coefficient_domains_remaps_formal_generators(kind,domain):
    a,b,y,z=sx.symbols('a b y z');ra,rb,ry,rz=sp.symbols('a b y z')
    nk1,nk2=(sx.ZZ.poly_ring(a),sx.ZZ.poly_ring(b)) if domain=='polynomial' else (sx.QQ.frac_field(a),sx.QQ.frac_field(b))
    rk1,rk2=(sp.ZZ.poly_ring(ra),sp.ZZ.poly_ring(rb)) if domain=='polynomial' else (sp.QQ.frac_field(ra),sp.QQ.frac_field(rb))
    p=getattr(sx,kind)([a,1],y+1,domain=nk1)
    q=getattr(sx,kind)([b,2],z+1,domain=nk2)
    rp=getattr(sp,kind)([ra,1],ry+1,domain=rk1)
    rq=getattr(sp,kind)([rb,2],rz+1,domain=rk2)
    # PurePoly aligns generator positions; ordinary Poly unifies their names.
    compare(p+q,rp+rq,kind)
    compare(p*q,rp*rq,kind)


def test_constructor_normalization_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y')
assert Poly({(2,):1,(0,):2},gens=x,field=True).domain==QQ
assert PurePoly((1,0,2),gens=x).degree()==2
assert Poly((x+1)**3,expand=False).gens==(x+1,)
assert Poly(x+y,wrt=y).gens==(y,x)
k=QQ.algebraic_field(sqrt(2))
p=Poly([k.unit,k.one],y+1,domain=k)
assert p.degree()==1 and (p*p).degree()==2
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=20)


def test_sparse_formal_generator_metadata_never_expands_a_million_degree_binomial():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
a,y=symbols('a y')
degree=10**6
p=Poly({(degree,):a,(0,):1},y+1,domain=ZZ[a])
assert p.degree()==degree and p.length()==2
assert p.degree_list()==(degree,) and p.total_degree()==degree
assert p.diff().degree()==degree-1
assert p.as_dict()=={(degree,):a,(0,):Integer(1)}
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=20)
