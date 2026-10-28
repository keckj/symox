"""Certified Piecewise specialization returns to the complete integral safely."""
import pytest
import symox as sx

x,a,b=sx.symbols('x a b')


def families(parameter=a):
    return [(1/(x+parameter),(x,0,1)),(1/(x-parameter),(x,0,1)),
            (1/(1+parameter*x),(x,0,1)),(1/x,(x,parameter,1)),
            (1/(x*x+parameter),(x,0,sx.oo))]


def apply(expression,old,new,method):
    if method=='subs':return expression.subs(old,new)
    if method=='xreplace':return expression.xreplace({old:new})
    return expression.subs({old:new},simultaneous=True)


@pytest.mark.parametrize('method',['subs','xreplace','simultaneous'])
@pytest.mark.parametrize('index',range(5))
@pytest.mark.parametrize('facts',[{}, {'imaginary':True}])
def test_symbolic_parameters_revert_to_the_exact_original_integral(index,method,facts):
    c=sx.Symbol('c',**facts);f,limit=families()[index]
    result=sx.integrate(f,limit);original=sx.Integral(f,limit)
    assert result.is_Piecewise
    actual=apply(result,a,c,method);expected=apply(original,a,c,method)
    assert actual==expected and actual.has(sx.Integral)
    for point in (-2,sx.Rational(-1,2),0,2):
        # Once the complete path is restored, subsequent substitutions and
        # doit have precisely the established Integral behavior.
        assert apply(actual,c,point,method)==apply(expected,c,point,method)
        assert apply(actual,c,point,method).doit()==apply(expected,c,point,method).doit()
    assert actual.doit()==expected.doit()


@pytest.mark.parametrize('method',['subs','xreplace','simultaneous'])
@pytest.mark.parametrize('index',range(5))
def test_symbolic_expression_and_repeated_doit_preserve_the_full_path(index,method):
    c=sx.Symbol('c',imaginary=True);f,limit=families()[index]
    result=sx.integrate(f,limit);original=sx.Integral(f,limit)
    actual=apply(result,a,c+1,method);expected=apply(original,a,c+1,method)
    assert actual==expected
    for _ in range(3):
        assert actual==expected
        actual=actual.doit();expected=expected.doit()
    assert apply(actual,c,-1,method).doit()==apply(expected,c,-1,method).doit()


@pytest.mark.parametrize('method',['subs','xreplace','simultaneous'])
@pytest.mark.parametrize('index',range(5))
@pytest.mark.parametrize('point',[0,2,sx.I,-2+sx.I])
def test_numeric_specializations_keep_normal_formulas_and_conditions(index,method,point):
    f,limit=families()[index];result=sx.integrate(f,limit)
    actual=apply(result,a,point,method)
    # Expected uses the same branch values with individually substituted
    # numeric conditions; no symbolic parameter is introduced here.
    branches=[(apply(pair[0],a,point,method),apply(pair[1],a,point,method)) for pair in result.args]
    assert actual==sx.Piecewise(*branches)


@pytest.mark.parametrize('method',['subs','xreplace','simultaneous'])
def test_partial_numeric_specialization_keeps_the_remaining_parameter_guard(method):
    c=sx.Symbol('c',positive=True);parameter=a+c
    original=sx.Integral(1/(x+parameter),(x,0,1));result=sx.integrate(1/(x+parameter),(x,0,1))
    actual=apply(result,a,1,method);expected=apply(original,a,1,method)
    assert actual==expected
    assert apply(actual,c,-sx.Rational(3,2),method).doit()==sx.nan


@pytest.mark.parametrize('method',['subs','xreplace','simultaneous'])
def test_fully_numeric_composed_parameter_specialization_keeps_its_value(method):
    f=1/(x+a+b);result=sx.integrate(f,(x,0,1));original=sx.Integral(f,(x,0,1))
    if method=='subs':actual=result.subs([(a,2),(b,3)])
    elif method=='xreplace':actual=result.xreplace({a:2,b:3})
    else:actual=result.subs({a:2,b:3},simultaneous=True)
    assert actual.doit()==original.subs({a:2,b:3}).doit()
    if method!='subs':assert not actual.has(sx.Integral)


@pytest.mark.parametrize('method',['subs','xreplace'])
def test_foreign_piecewise_with_same_guard_and_wrong_value_is_not_transformed(method):
    f,limit=families()[0];certified=sx.integrate(f,limit)
    condition=certified.args[0][1];original=sx.Integral(f,limit)
    foreign=sx.Piecewise((a+1,condition),(original,True))
    c=sx.Symbol('c',imaginary=True)
    actual=apply(foreign,a,c,method)
    expected=sx.Piecewise((c+1,apply(condition,a,c,method)),(apply(original,a,c,method),True))
    assert actual==expected and actual.is_Piecewise
    assert actual!=apply(original,a,c,method)


@pytest.mark.parametrize('method',['subs','xreplace'])
def test_foreign_domain_or_original_integral_is_not_authorized_by_a_matching_formula(method):
    f,limit=families()[0];certified=sx.integrate(f,limit);value=certified.args[0][0]
    different=sx.Integral(1/(x*x+a),(x,0,1))
    wrong_domain=sx.And(sx.Contains(a,sx.S.Complexes),sx.Not(sx.Contains(-a,sx.Interval(0,2))))
    for foreign in (sx.Piecewise((value,wrong_domain),(sx.Integral(f,limit),True)),
                    sx.Piecewise((value,a>0),(sx.Integral(f,limit),True)),
                    sx.Piecewise((value,certified.args[0][1]),(different,True))):
        actual=apply(foreign,a,b,method)
        expected=sx.Piecewise(*[(apply(pair[0],a,b,method),apply(pair[1],a,b,method)) for pair in foreign.args])
        assert actual==expected and actual.is_Piecewise


@pytest.mark.parametrize('method',['subs','xreplace','simultaneous'])
def test_bound_and_identity_substitutions_do_not_remove_the_piecewise(method):
    f,limit=families()[0];result=sx.integrate(f,limit)
    assert apply(result,a,a,method)==result
    assert apply(result,b,a,method)==result
    actual=apply(result,x,b,method)
    if method!='xreplace':assert actual==result
    else:
        assert actual.is_Piecewise and actual.args[-1][0]==sx.Integral(1/(b+a),(b,0,1))
    # A parameter replaced by the binding variable is a different operation:
    # preserve the established Integral capture behavior.
    assert apply(result,a,x,method)==apply(sx.Integral(f,limit),a,x,method)


def test_exact_whole_node_and_fallback_replacements_keep_their_priority():
    f,limit=families()[0];original=sx.Integral(f,limit);result=sx.integrate(f,limit)
    assert result.xreplace({result:7,a:b})==7
    assert result.subs(result,7)==7
    expected=sx.Piecewise((result.args[0][0],result.args[0][1]),(b,True))
    assert result.subs(original,b)==expected
    assert result.xreplace({original:b})==expected


def test_masked_simultaneous_numeric_rules_retain_the_established_integral_value():
    f,limit=families()[0];result=sx.integrate(f,limit);original=sx.Integral(f,limit)
    actual=result.subs({1/a:2},simultaneous=True)
    expected=original.subs({1/a:2},simultaneous=True)
    assert actual==expected and actual.doit()==sx.log(3)
    assert result.subs(1/a,2)==sx.log(sx.Rational(3,2))+sx.log(2)


@pytest.mark.parametrize('method',['xreplace','simultaneous'])
def test_maps_are_simultaneous_and_do_not_resubstitute_new_symbols(method):
    c=sx.Symbol('c',imaginary=True);f,limit=families()[0]
    result=sx.integrate(f,limit);original=sx.Integral(f,limit)
    rules={a:c,c:-sx.Rational(1,2)}
    actual=result.xreplace(rules) if method=='xreplace' else result.subs(rules,simultaneous=True)
    expected=original.xreplace(rules) if method=='xreplace' else original.subs(rules,simultaneous=True)
    assert actual==expected and c in actual.free_symbols
    assert apply(actual,c,-sx.Rational(1,2),'subs').doit()==sx.nan


def test_algebraic_parameter_patterns_and_wrapped_piecewise_are_supported():
    f=1/(x+a*a);original=sx.Integral(f,(x,0,1));result=sx.integrate(f,(x,0,1))
    assert result.subs(a*a,b)==original.subs(a*a,b)
    assert result.xreplace({a*a:b})==original.xreplace({a*a:b})
    f,limit=families()[0];original=sx.Integral(f,limit);result=sx.integrate(f,limit)
    assert (sx.sin(result)+1).subs(a,b)==sx.sin(original.subs(a,b))+1


def test_custom_protocol_and_fatal_errors_keep_the_original_piecewise_path():
    calls=[]
    class Hook(sx.Symbol):
        def _eval_subs(self,old,new):calls.append((old,new));raise KeyboardInterrupt('foreign hook')
    h=Hook('h');condition=sx.And(sx.Contains(a,sx.S.Complexes),sx.Not(sx.Contains(-a,sx.Interval(0,1))))
    foreign=sx.Piecewise((h+1,condition),(sx.Integral(1/(x+a),(x,0,1)),True))
    with pytest.raises(KeyboardInterrupt,match='foreign hook'):foreign.subs(a,b)
    assert calls
    # A failed protocol call must not poison later ordinary substitutions.
    result=sx.integrate(1/(x+a),(x,0,1))
    assert result.subs(a,b)==sx.Integral(1/(x+b),(x,0,1))


def test_invalid_exact_bound_replacement_preserves_the_integral_error():
    f,limit=families()[0];original=sx.Integral(f,limit);result=sx.integrate(f,limit)
    rule={x:2,a:b}
    try:expected=original.xreplace(rule)
    except Exception as error:
        with pytest.raises(type(error),match='^'+__import__('re').escape(str(error))+'$'):result.xreplace(rule)
    else:assert result.xreplace(rule)==expected


def test_hardening_has_no_runtime_sympy_dependency():
    import os,subprocess,sys
    script="""import sys
sys.modules['sympy']=None
import symox as s
x,a=s.symbols('x a');b=s.Symbol('b',imaginary=True)
for f,limit in [(1/(x+a),(x,0,1)),(1/x,(x,a,1)),(1/(x*x+a),(x,0,s.oo))]:
 result=s.integrate(f,limit);original=s.Integral(f,limit)
 for method in ('subs','xreplace'):
  actual=result.subs(a,b) if method=='subs' else result.xreplace({a:b})
  expected=original.subs(a,b) if method=='subs' else original.xreplace({a:b})
  assert actual==expected
  assert actual.subs(b,-s.Rational(1,2)).doit()==expected.subs(b,-s.Rational(1,2)).doit()
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('part',['value','guard','fallback'])
def test_exact_maps_that_override_formula_or_domain_keep_ordinary_piecewise_semantics(part):
    f,limit=families()[0];result=sx.integrate(f,limit);original=sx.Integral(f,limit)
    old={'value':result.args[0][0],'guard':result.args[0][1],'fallback':original}[part]
    new=sx.true if part=='guard' else sx.Integer(7)
    rule={a:b,old:new}
    expected=sx.Piecewise(*[(pair[0].xreplace(rule),pair[1].xreplace(rule)) for pair in result.args])
    assert result.xreplace(rule)==expected


def test_irrelevant_rules_and_bound_renaming_do_not_change_parameter_substitution_policy():
    f,limit=families()[0];result=sx.integrate(f,limit);original=sx.Integral(f,limit)
    z=sx.Symbol('z');c=sx.Symbol('c',imaginary=True)
    rule={a:c,z:1,x:b}
    assert result.xreplace(rule)==original.xreplace(rule)
    assert result.subs(rule,simultaneous=True)==original.subs(rule,simultaneous=True)


def test_reconstructed_and_pickled_certified_rules_keep_the_same_recognition():
    import pickle
    d=sx.Dummy('d');f=1/(d+a)
    result=sx.integrate(f,(d,0,1))
    for rebuilt in (result.func(*result.args),pickle.loads(pickle.dumps(result))):
        original=rebuilt.args[-1][0]
        assert rebuilt.subs(a,b)==original.subs(a,b)
        assert rebuilt.xreplace({a:b})==original.xreplace({a:b})
