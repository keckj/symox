"""Targeted simplification flags must affect traversal or fail explicitly."""
import pytest
import sympy as sp
import symox as sx


def oracle(function, expression, **options):
    actual = getattr(sx,function)(expression,**options)
    expected = getattr(sp,function)(expression.to_sympy(),**options)
    assert str(actual) == str(expected)
    return actual


@pytest.mark.parametrize('deep',[False,True])
@pytest.mark.parametrize('function,kind',[
    ('expand_log','log'), ('expand_mul','mul'), ('expand_power_base','base'),
])
def test_deep_controls_nested_function_arguments(function,kind,deep):
    x,y,z=sx.symbols('x y z',positive=True)
    child = {'log':sx.log(x*y),'mul':x*(y+z),'base':(x*y)**z}[kind]
    expression=sx.sin(child)
    oracle(function,expression,deep=deep)
    if not deep: assert getattr(sx,function)(expression,deep=False) == expression


@pytest.mark.parametrize('deep',[False,True])
def test_expand_log_root_is_transformed_even_when_shallow(deep):
    x,y=sx.symbols('x y',positive=True)
    oracle('expand_log',sx.log(x*y),deep=deep)


def test_explicit_force_is_applied_only_by_supported_routines():
    x,y,z=sx.symbols('x y z')
    assert sx.expand_log(sx.log(x*y)) == sx.log(x*y)
    oracle('expand_log',sx.log(x*y),force=True)
    assert sx.expand_power_base((x*y)**z) == (x*y)**z
    oracle('expand_power_base',(x*y)**z,force=True)


@pytest.mark.parametrize('function',['expand_power_exp','expand'])
@pytest.mark.parametrize('kind',['pow','exp','coefficient','nested'])
def test_power_exponent_expansion_preserves_explicit_factors(function,kind):
    x=sx.Symbol('x',nonzero=True)
    a,b,y=sx.symbols('a b y')
    base = sx.exp(a+b) if kind=='exp' else x**(a+b)
    expression=2*base+y if kind=='coefficient' else sx.sin(base) if kind=='nested' else base
    actual=oracle(function,expression)
    assert actual != expression
    assert sx.powsimp(actual) == sx.powsimp(expression)
    if function=='expand':
        unchanged=sx.expand(expression,power_exp=False)
        assert unchanged == expression


@pytest.mark.parametrize('deep',[False,True])
def test_power_expansion_shallow_traversal(deep):
    a,b=sx.symbols('a b')
    expression=sx.sin(sx.exp(a+b))
    oracle('expand_power_exp',expression,deep=deep)
    x=sx.Symbol('x',nonzero=True)
    oracle('expand_power_exp',2*x**(a+b),deep=deep)


def test_disabling_power_exp_still_expands_inside_exponential_arguments():
    x,y,z=sx.symbols('x y z')
    for expression in [sx.exp(x*(y+z)),sx.exp((x+y)**2)]:
        oracle('expand',expression,power_exp=False)


def test_powsimp_deep_and_exp_mode():
    a,b=sx.symbols('a b')
    held=sx.Mul(sx.exp(a),sx.exp(b),evaluate=False)
    expression=sx.sin(held)
    assert sx.powsimp(expression,deep=False,combine='exp') == expression
    result=sx.powsimp(expression,deep=True,combine='exp')
    assert result == sx.sin(sx.exp(a+b))
    assert result != expression
    for options in [{'force':True},{'combine':'base'}]:
        with pytest.raises(NotImplementedError): sx.powsimp(expression,**options)
    with pytest.raises(ValueError): sx.powsimp(expression,combine='invalid')


@pytest.mark.parametrize('options',[{'tolerance':1e-5},{'full':True},{'rational':True},{'rational':False}])
def test_unsupported_nsimplify_options_are_not_silently_ignored(options):
    with pytest.raises(NotImplementedError): sx.nsimplify(1.4142135623730951,**options)


@pytest.mark.parametrize('tolerance',[0,-1,float('inf'),float('nan')])
def test_nsimplify_invalid_tolerances(tolerance):
    with pytest.raises(ValueError): sx.nsimplify(.5,tolerance=tolerance)


def test_nsimplify_existing_constants_and_defaults_preserved():
    assert sx.nsimplify(.5) == sx.Rational(1,2)
    assert sx.nsimplify(1.4142135623730951) == sx.sqrt(2)
    assert sx.nsimplify(3.141592653589793,constants=[sx.pi]) == sx.pi
    assert sx.nsimplify(.5,full=False,tolerance=None,rational=None) == sx.Rational(1,2)
