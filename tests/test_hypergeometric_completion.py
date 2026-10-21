import pickle
import pytest
import mpmath as mp
import symox as s
import sympy as sp
from symox.functions.special.hyper import TupleArg, TupleParametersBase
import importlib
sh=importlib.import_module('symox.functions.special.hyper')
ph=importlib.import_module('sympy.functions.special.hyper')


@pytest.mark.parametrize('parameters',[
    ([],[],[0],[]),([],[],[s.Rational(1,2)],[]),
    ([],[1],[0],[]),([],[s.Rational(3,2)],[s.Rational(1,2)],[]),
    ([],[],[0,s.Rational(1,2)],[]),([0],[],[0],[]),
])
def test_meijer_metadata_derivatives_and_closed_forms(parameters):
    z,t=s.symbols('z t',positive=True)
    an,ao,bm,bo=parameters
    f=s.meijerg(an,ao,bm,bo,z)
    expected=sp.meijerg(*[[v.to_sympy() if hasattr(v,'to_sympy') else v for v in group] for group in parameters],z.to_sympy())
    assert isinstance(f,TupleParametersBase)
    assert f.to_sympy()==expected
    assert f==s.meijerg([an,ao],[bm,bo],z)
    assert f.argument==z and tuple(f.an)==tuple(map(s.sympify,an))
    assert f.nu.to_sympy()==expected.nu and f.delta.to_sympy()==expected.delta
    assert sp.simplify(f.integrand(t).to_sympy()-expected.integrand(t.to_sympy()))==0
    assert sp.simplify(f.diff(z).to_sympy()-expected.diff(z.to_sympy()))==0
    assert f.fdiff()==f.diff(z)
    assert f.get_period().to_sympy()==expected.get_period()
    actual=s.hyperexpand(f)
    wanted=sp.hyperexpand(expected)
    if len(an)+len(ao)<len(bm)+len(bo):
        assert actual!=f
        assert sp.simplify(actual.to_sympy()-wanted)==0
    assert pickle.loads(pickle.dumps(f))==f


@pytest.mark.parametrize('z',[s.Rational(1,2),2,5])
def test_meijer_numeric_and_repeated_pole_guard(z):
    for parameters in [([],[],[s.Rational(1,3)],[]),([],[2],[0],[]),([],[],[0,s.Rational(1,2)],[])]:
        f=s.meijerg(*parameters,z)
        got=complex(f.evalf())
        expected=complex(mp.meijerg([list(parameters[0]),list(parameters[1])],[[float(v) for v in parameters[2]],list(parameters[3])],float(z)))
        assert abs(got-expected)<3e-12*max(1,abs(expected))
    boundary=s.meijerg([],[2],[0],[],1)
    assert boundary.evalf()==boundary
    f=s.meijerg([],[],[0,0],[],z)
    assert s.hyperexpand(f)==f
    with pytest.raises(ValueError):s.meijerg([1],[],[0],[],z)
    with pytest.raises(ValueError):s.meijerg([],[],[s.oo],[],z)
    with pytest.raises(TypeError):s.meijerg([],[],z)


@pytest.mark.parametrize('values',[(1,2,3,4,.2,.3),(s.Rational(1,2),1,2,3,.1,.25),(-3,2,3,4,2.,3.),(2,1,3,5,-.4,.2)])
def test_appell_f1_native_numeric_and_derivatives(values):
    a,b1,b2,c,x,y=values
    actual=s.appellf1(*values)
    expected=mp.appellf1(*[float(v) for v in values])
    assert abs(complex(actual.evalf())-complex(expected))<5e-12*max(1,abs(expected))
    x,y=s.symbols('x y');a,b1,b2,c=map(s.sympify,(a,b1,b2,c))
    f=s.appellf1(a,b1,b2,c,x,y)
    assert f.diff(x)==a*b1/c*s.appellf1(a+1,b1+1,b2,c+1,x,y)
    assert f.diff(y)==a*b2/c*s.appellf1(a+1,b1,b2+1,c+1,x,y)
    assert f.fdiff(5)==f.diff(x)
    oracle=sp.appellf1(*(v.to_sympy() for v in (a,b1,b2,c,x,y)))
    assert f.to_sympy()==oracle
    assert s.appellf1(a,b1,b2,c,0,0)==1


REP_NAMES=['power1','power2','log1','log2','atanh','asin1','asin2','sqrts1','sqrts2','cosasin','sinasin']


@pytest.mark.parametrize('suffix',REP_NAMES)
@pytest.mark.parametrize('z',[s.Rational(1,4),s.Rational(-1,4),2,-2])
def test_hyper_representatives_rewrite_and_branch_formulas(suffix,z):
    cls=getattr(sh,'HyperRep_'+suffix);oracle=getattr(ph,'HyperRep_'+suffix)
    leading=(s.Rational(2,3),) if suffix in ('power1','power2','sqrts1','sqrts2','cosasin','sinasin') else ()
    f=cls(*leading,z);reference=oracle(*(v.to_sympy() if hasattr(v,'to_sympy') else v for v in leading+(z,)))
    assert isinstance(f,s.HyperRep)
    assert f.to_sympy()==reference
    for rule in ('nonrep','nonrepsmall'):
        got=f.rewrite(rule).to_sympy();want=reference.rewrite(rule)
        assert sp.simplify(got-want)==0
    for branch in (-1,0,1,2):
        args=leading+(s.Rational(3,2),s.Integer(branch))
        got=cls._expr_big(*args).to_sympy()
        want=reference._expr_big(*(v.to_sympy() for v in args))
        assert sp.simplify(got-want)==0


def test_tuple_helpers_and_nested_representative_rewrite():
    x=s.Symbol('x')
    values=TupleArg(s.sin(x),x*x+1)
    assert values.as_leading_term(x)==TupleArg(x,1)
    assert values.limit(x,0)==TupleArg(0,1)
    assert isinstance(s.hyper([1],[2],x),TupleParametersBase)
    f=sh.HyperRep_log1(x)
    assert (1+f).rewrite('nonrepsmall')==1+s.log(1-x)
    f=sh.HyperRep_log1(s.Rational(1,2)*s.exp_polar(2*s.pi*s.I))
    want=ph.HyperRep_log1(sp.Rational(1,2)*sp.exp_polar(2*sp.pi*sp.I)).rewrite('nonrep')
    assert sp.simplify(f.rewrite('nonrep').to_sympy()-want)==0


def test_meijer_piecewise_reexpansion_matches_reference():
    # The reference re-expands the G-function kept in the default branch.
    z=s.Symbol('z'); zr=sp.Symbol('z')
    for f, g in [(s.meijerg([],[2],[0],[],z), sp.meijerg([],[2],[0],[],zr)),
                 (s.meijerg([1],[],[],[0],z), sp.meijerg([1],[],[],[0],zr))]:
        expanded=s.hyperexpand(f)
        want=sp.hyperexpand(g)
        assert str(expanded)==str(want)
        assert str(s.hyperexpand(expanded))==str(sp.hyperexpand(want))
        assert s.expand_func(expanded)==expanded
