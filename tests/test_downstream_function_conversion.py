"""Function reconstruction and accurate scalar conversion used downstream."""
import inspect
import itertools
import math
import pickle

import pytest
import sympy as sp
import symox as sx


def function_heads(module):
    return [module.Function('f'), module.Function('g'),
            module.Function('f', real=True), module.Function('f', nargs=1),
            module.Function('f', nargs=(1, 2))]


def test_undefined_function_instance_queries_keep_head_assumptions_and_arity():
    actual, expected = function_heads(sx), function_heads(sp)
    for i,j in itertools.product(range(len(actual)), repeat=2):
        for arguments in [(1,), (1,2)]:
            try:
                reference = expected[i](*arguments)
            except TypeError:
                continue
            assert isinstance(actual[i](*arguments), actual[j]) == isinstance(reference, expected[j])
    for value in [1, sx.Symbol('x'), sx.sin(sx.Symbol('x')), actual[0]]:
        assert not isinstance(value, actual[0])
    restored = pickle.loads(pickle.dumps(actual[0](sx.Symbol('x'))))
    assert isinstance(restored, actual[0])


@pytest.mark.parametrize('items', [(), (1,2), (sx.Symbol('x'),sx.Tuple(1,2))])
def test_tuple_reconstruction_options_do_not_change_container_evaluation(items):
    for options in [{}, {'evaluate':False}, {'evaluate':True}, {'unused':True}]:
        result = sx.Tuple(*items, **options)
        assert result == sx.Tuple(*items)
        assert result.func(*result.args, **options) == result
    assert inspect.signature(sx.Tuple) == inspect.signature(sp.Tuple)
    value = sx.Tuple(sx.Symbol('x'), sx.Integer(1), sympify=False)
    assert value.args == (sx.Symbol('x'),sx.Integer(1))


@pytest.mark.parametrize('power', [20,50,100,300])
def test_float_conversion_keeps_cancellation_at_target_binary_precision(power):
    actual = sx.sqrt(sx.Integer(10)**power+1)-sx.Integer(10)**(power//2)
    # Rationalization is an exact identity and gives an independent stable
    # reference, even beyond the default SymPy float evaluation budget.
    expected = 1/(sp.sqrt(sp.Integer(10)**power+1)+sp.Integer(10)**(power//2))
    assert float(actual) == float(expected)
    assert complex(actual) == complex(expected)
    assert float(actual) > 0


@pytest.mark.parametrize('power', [20,100,1000])
def test_float_rejects_genuine_small_imaginary_components(power):
    actual = 1+sx.I*sx.Rational(1,10**power)
    expected = 1+sp.I*sp.Rational(1,10**power)
    with pytest.raises(TypeError): float(expected)
    with pytest.raises(TypeError): float(actual)
    assert complex(actual) == complex(expected)


@pytest.mark.parametrize('name', ['nan','oo','zoo'])
def test_explicit_nonfinite_scalar_conversions_match_reference(name):
    actual, expected = getattr(sx,name), getattr(sp,name)
    for operation in [float,complex]:
        try:
            reference = operation(expected)
        except TypeError:
            with pytest.raises(TypeError): operation(actual)
            continue
        native = operation(actual)
        components = lambda value: (value.real,value.imag) if isinstance(value,complex) else (value,)
        for result, target in zip(components(native),components(reference)):
            assert math.isnan(result) if math.isnan(target) else result == target


def test_held_derivative_has_public_constructor_identity_and_reconstructs():
    x = sx.Symbol('x')
    f = sx.Function('f')
    expression = sx.Derivative(f(x), (x,2))
    assert expression.func is sx.Derivative
    assert expression.func(*expression.args,evaluate=False) == expression


@pytest.mark.parametrize('names', [r'\Delta\ H \Delta\ S R T', r'a\,b c',
    r'a\:b x:3', r'a\ b', r'a\,b,', r'foo\ bar:2', 'x y z',
    r'alpha\ beta gamma\:delta', '\0 x\\ y'])
def test_symbol_names_preserve_escaped_separators(names):
    actual, expected = sx.symbols(names), sp.symbols(names)
    assert isinstance(actual,tuple) == isinstance(expected,tuple)
    assert [value.name for value in (actual if isinstance(actual,tuple) else (actual,))] == [
        value.name for value in (expected if isinstance(expected,tuple) else (expected,))]


@pytest.mark.parametrize('output', ['scalar','list','empty'])
def test_compiled_numerical_function_accepts_only_zero_dimensional_arguments(output):
    import numpy as np
    x,y=sx.symbols('x y')
    expr=x+2*y if output=='scalar' else [x+2*y,x-y] if output=='list' else []
    function=sx.lambdify((x,y),expr)
    for arguments in [(np.array(2.),np.array(3.)), (sx.Rational(2),sx.Rational(3)),
                      (np.float64(2),np.float64(3))]:
        result=function(*arguments)
        assert result == (8 if output=='scalar' else [8,-1] if output=='list' else [])


def test_callable_slots_are_confined_to_lambdas_and_predicates():
    import symox.diffgeom  # Loading optional geometry must not add scalar call slots.
    x=sx.Symbol('x')
    for value in [x,sx.Integer(1),sx.pi,sx.sin(x),sx.Tuple(1,2),sx.Q.positive(x)]:
        assert not callable(value)
        with pytest.raises(TypeError):value(1)
    for value in [sx.Lambda(x,x+1),sx.Id,sx.Q.positive]:
        assert callable(value)
        assert callable(pickle.loads(pickle.dumps(value)))
    assert sx.Lambda(x,x+1)(2)==3
    assert sx.Id(2)==2
    assert sx.Q.positive(x)==sx.Q.positive(x)
    transformed=sx.Lambda(x,x+sx.Symbol('a')).subs(sx.Symbol('a'),2)
    assert callable(transformed) and transformed(3)==5
    from symox.diffgeom import Differential
    from symox.diffgeom.rn import R2_r
    field=R2_r.base_scalar(0)
    zero=Differential(Differential(field))
    assert callable(zero) and zero(R2_r.base_vector(0))==0
    assert pickle.loads(pickle.dumps(zero))(R2_r.base_vector(0))==0
    assert not callable(sx.S.Zero)


@pytest.mark.parametrize('power',[1,2,3])
def test_differential_field_callability_survives_pickle_and_identity_operations(power):
    from symox.diffgeom.rn import R2_r
    field=R2_r.base_scalar(0)
    point=R2_r.point([2,3])
    value=field**power
    assert callable(value) and value(point)==2**power
    clone=pickle.loads(pickle.dumps(value))
    assert callable(clone) and clone(point)==2**power
    assert (clone+1)(point)==2**power+1
    assert (-clone)(point)==-(2**power)
    assert not callable(sx.Dummy('ordinary_dummy')**2)


def test_native_scalar_reconstruction_does_not_require_caller_imports():
    from symox.diffgeom.rn import R2_r
    x = sx.Symbol('restricted_x')
    namespace = dict(x=x, Lambda=sx.Lambda, Tuple=sx.Tuple,
                     field=R2_r.base_scalar(0), Derivative=sx.Derivative)
    def restricted(source):
        return eval(source, {'__builtins__': {}}, namespace)
    assert restricted('-x') == -x
    assert restricted('x-x') is sx.S.Zero
    assert isinstance(restricted('Tuple(x, 1)'), sx.Tuple)
    function = restricted('Lambda(x, x+1)')
    assert callable(function) and function(2) == 3
    derivative = restricted('Derivative(x**2, x, evaluate=False)')
    assert eval('value.func', {'__builtins__': {}}, dict(value=derivative)) is sx.Derivative
    assert callable(restricted('field**2'))


def test_numerical_functions_preserve_complex_arrays_and_masked_array_protocols():
    import numpy as np
    x=sx.Symbol('x')
    f=sx.lambdify(x,x*x+1)
    for value in [2+3j,np.array([1+2j,3-4j]),
                  np.ma.array([1.,2.,3.],mask=[False,True,False])]:
        result=f(value)
        expected=value*value+1
        assert type(result) is type(expected)
        np.testing.assert_array_equal(result,expected)
        if isinstance(expected,np.ma.MaskedArray):
            np.testing.assert_array_equal(result.mask,expected.mask)


def test_numerical_functions_use_external_arithmetic_instead_of_float_coercion():
    class Tagged:
        def __init__(self,value):self.value=value
        def __float__(self):raise AssertionError('Operand metadata must be retained')
        def __pow__(self,exponent):return type(self)(self.value**exponent)
        def __add__(self,other):return type(self)(self.value+other)
        __radd__=__add__
    x=sx.Symbol('x')
    f=sx.lambdify(x,x*x+1)
    result=f(Tagged(3))
    assert isinstance(result,Tagged) and result.value==10


@pytest.mark.parametrize('name', ['__add__','__radd__','__sub__','__rsub__','__mul__',
    '__rmul__','__truediv__','__rtruediv__','__pow__','__rpow__'])
def test_unsupported_scalar_conversion_allows_foreign_reflected_operations(name):
    class External:pass
    value=External()
    assert getattr(sx.Integer(7),name)(value) is NotImplemented
    if not name.startswith('__r'):
        reflected='__r'+name[2:]
        calls=[]
        def operation(self,other):
            calls.append(other)
            return 'foreign result'
        setattr(External,reflected,operation)
        import operator
        function={'__add__':operator.add,'__sub__':operator.sub,'__mul__':operator.mul,
                  '__truediv__':operator.truediv,'__pow__':operator.pow}[name]
        assert function(sx.Integer(7),value)=='foreign result'
        assert calls==[sx.Integer(7)]


@pytest.mark.parametrize('error',[TypeError,ValueError,RuntimeError])
def test_conversion_hook_errors_propagate_without_using_foreign_arithmetic(error):
    class External:
        def _sympy_(self):raise error('conversion hook failed')
        def __rmul__(self,other):raise AssertionError('A conversion error must propagate')
    with pytest.raises(error,match='conversion hook failed'):
        sx.Integer(7)*External()


@pytest.mark.parametrize('offset,sign',[(1,1),(-1,1),(1,-1),(-1,-1)])
def test_literal_complex_subnormal_components_are_rounded_only_once(offset,sign):
    from fractions import Fraction
    fine=2**200
    numerator=sign*(fine+offset)
    denominator=2**1075*fine
    value=sx.Rational(numerator,denominator)
    assert complex(1+sx.I*value).imag.hex()==float(Fraction(numerator,denominator)).hex()
