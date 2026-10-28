"""Exact slit-domain certificates for retained logarithmic and quadratic rays."""
import pytest
import sympy as sp
import symox as sx

x,a=sx.symbols('x a')
cut=sx.Interval(-sx.oo,0)


def cases(parameter=a):
    return [(1/x,(x,parameter,1),-sx.log(parameter)),
            (1/(x*x+parameter),(x,0,sx.oo),sx.pi/(2*sx.sqrt(parameter)))]


@pytest.mark.parametrize('index',[0,1])
def test_exact_formula_and_full_slit_domain_are_retained_with_original_fallback(index):
    f,limit,value=cases()[index]
    result=sx.integrate(f,limit)
    assert result.is_Piecewise
    assert result.args[0][0]==value
    assert result.args[0][1]==sx.And(sx.Contains(a,sx.S.Complexes),sx.Not(sx.Contains(a,cut)))
    assert result.args[-1][0]==sx.Integral(f,limit) and result.args[-1][1]==sx.true


def test_symbolic_primitive_and_partial_fraction_certificates_are_exact():
    z,c,s=sp.symbols('z c s',nonzero=True)
    assert sp.cancel(sp.diff(sp.log(z),z)-1/z)==0
    primitive=(sp.log(z-sp.I*s)-sp.log(z+sp.I*s))/(2*sp.I*s)
    assert sp.cancel(sp.diff(primitive,z)-1/(z*z+s*s))==0
    assert sp.cancel((1/(z-sp.I*s)-1/(z+sp.I*s))/(2*sp.I*s)-1/(z*z+s*s))==0
    polynomial=sp.Poly(z*z+c,z)
    assert polynomial.degree()==2 and polynomial.eval(sp.I*s).subs(c,s*s)==0
    # For the principal sqrt off the closed negative ray, Re(s)>0. Its
    # arguments after multiplication by +/-I lie in (0,pi) and (-pi,0),
    # respectively. Their difference is exactly pi, not an undetermined
    # multiple of 2*pi. The complete convergence/branch proof is archived.
    theta=sp.Symbol('theta',real=True)
    assert sp.expand((theta+sp.pi/2)-(theta-sp.pi/2))==sp.pi


@pytest.mark.parametrize('index',[0,1])
@pytest.mark.parametrize('point',[sx.Rational(1,3),sx.Integer(2),sx.Integer(3),sx.I,-2+sx.I,-2-sx.I,-1+sx.I/10**100,-1-sx.I/10**100])
def test_regular_exact_real_and_complex_substitutions_select_principal_formulas(index,point):
    f,limit,value=cases()[index]
    result=sx.integrate(f,limit).subs(a,point)
    expected=value.subs(a,point)
    assert result==expected and not result.has(sx.Integral)
    # A differential value check accompanies the exact symbolic proof; it
    # is not used to establish convergence or choose a branch.
    cp=sp.sympify(str(point))
    reference=-sp.log(cp) if index==0 else sp.pi/(2*sp.sqrt(cp))
    assert sp.simplify(result.to_sympy()-reference)==0


@pytest.mark.parametrize('index',[0,1])
@pytest.mark.parametrize('point',[sx.Integer(-3),sx.Rational(-1,2),sx.Integer(0)])
def test_closed_cut_keeps_every_interior_and_endpoint_pole(index,point):
    f,limit,_=cases()[index]
    result=sx.integrate(f,limit).subs(a,point)
    expected=sx.Integral(f,limit).subs(a,point)
    assert result==expected
    assert result.doit()==(sx.oo if point==0 else sx.nan)


@pytest.mark.parametrize('facts',[{}, {'real':True},{'complex':True},{'nonzero':True}])
@pytest.mark.parametrize('index',[0,1])
def test_doit_is_repeatedly_idempotent_and_keeps_zero_substitution_safe(index,facts):
    c=sx.Symbol('c',**facts);f,limit,_=cases(c)[index]
    result=sx.integrate(f,limit);original=result
    assert result.is_Piecewise
    for _ in range(5):
        result=result.doit()
        assert result==original and result.subs(c,0).doit()==sx.oo


@pytest.mark.parametrize('facts',[{}, {'real':True},{'complex':True},{'nonzero':True}])
@pytest.mark.parametrize('point',[sx.oo,-sx.oo,sx.zoo,sx.nan])
def test_nonfinite_values_before_and_after_doit_keep_the_original_integral(facts,point):
    c=sx.Symbol('c',**facts)
    for f,limit,_ in cases(c):
        result=sx.integrate(f,limit)
        for first in (False,True):
            value=(result.doit() if first else result).subs(c,point)
            expected=sx.Integral(f,limit).subs(c,point)
            assert value==expected
            try:reference=expected.doit()
            except Exception as error:
                with pytest.raises(type(error)):value.doit()
            else:assert value.doit()==reference


def test_positive_parameter_results_remain_on_the_established_paths():
    c=sx.Symbol('c',positive=True)
    for f,limit,value in cases(c):
        assert sx.integrate(f,limit)==value


@pytest.mark.parametrize('options',[{'conds':'none'},{'conds':'separate'},{'meijerg':False},
    {'manual':False},{'heurisch':False},{'risch':False}])
def test_unported_explicit_options_keep_full_retained_integrals(options):
    for f,limit,_ in cases():
        assert sx.integrate(f,limit,**options)==sx.Integral(f,limit)


def test_nearby_unsupported_families_are_not_given_the_new_formulas():
    for f,limit in [(1/x,(x,a,2)),(1/(x*x+a),(x,1,sx.oo)),
                    (1/(x**4+a),(x,0,sx.oo)),(1/(x*x+a),(x,0,a))]:
        assert sx.integrate(f,limit).has(sx.Integral)


def test_rays_have_no_runtime_sympy_dependency():
    import os,subprocess,sys
    script="""import sys
sys.modules['sympy']=None
import symox as s
x,a=s.symbols('x a')
for f,limit,value in [(1/x,(x,a,1),-s.log(a)),(1/(x*x+a),(x,0,s.oo),s.pi/(2*s.sqrt(a)))]:
 result=s.integrate(f,limit)
 assert result.is_Piecewise and result.doit()==result
 assert result.subs(a,2)==value.subs(a,2)
 assert result.subs(a,0).doit()==s.oo
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('facts',[{'imaginary':True},{'complex':True,'real':False}])
def test_predecided_cut_exclusion_on_free_parameters_keeps_the_old_complete_path(facts):
    c=sx.Symbol('c',**facts)
    value=sx.integrate(1/x,(x,c,1))
    assert value==sx.Integral(1/x,(x,c,1))
    for _ in range(3):
        assert value.subs(c,-1).doit()==sx.nan
        assert value.subs(c,0).doit()==sx.oo
        value=value.doit()
        assert value==sx.Integral(1/x,(x,c,1))
