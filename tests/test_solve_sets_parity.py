"""Solution sets of solve against hardcoded SymPy 1.14 outputs: equations in
powers of a non-E base, equations with Float coefficients, hyperbolic
equations, exponentials of negative arguments, mixed power bases and the
logarithmic Lambert forms."""

import pytest

import symox as sx


NAMES = 'x a b'


def solve(source, **flags):
    namespace = dict(zip(NAMES.split(), sx.symbols(NAMES)))
    for name in ('exp', 'log', 'sqrt', 'S', 'Rational', 'Eq', 'sinh', 'cosh', 'tanh', 'sech', 'coth', 'asinh'):
        namespace[name] = getattr(sx, name)
    return str(sx.solve(eval(source, namespace), namespace['x'], **flags))


POWER_BASE = [
    ('4**(2*x) - 3*4**x + 2', '[0, 1/2, I*pi/log(2)]'),
    ('2**(4*x) - 3*2**(2*x) + 2', '[0, 1/2, I*pi/log(2)]'),
    ('2**(2*x) - 2', '[1/2]'),
    ('4**(2*x) - 4', '[1/2]'),
    ('(2**x - 1)*(2**(2*x) - 2)', '[0, 1/2]'),
    ('2**(2*x) - 3*2**x + 2', '[0, 1]'),
    ('8**x - 2**x', '[0, I*pi/log(2)]'),
    ('4**x - 2**(x + 1)', '[1]'),
    ('2**(2*x) + 2**x - 2', '[0, 1 + I*pi/log(2)]'),
]


@pytest.mark.parametrize('source, expected', POWER_BASE)
def test_power_base_equations_give_the_reference_solutions(source, expected):
    assert solve(source) == expected


def test_power_base_with_check_false_gives_the_reference_solutions():
    assert solve('4**(2*x) - 3*4**x + 2', check=False) == '[0, 1/2, I*pi/log(2)]'


FLOATS = [
    ('4.0**x - 8', '[1.50000000000000, 1.5 + 4.53236014182719*I]'),
    ('exp(2.0*x) - 3', '[0.549306144334055, 0.549306144334055 + 3.14159265358979*I]'),
    ('2.0**x - 8', '[3.00000000000000]'),
    ('2**x - 8.0', '[3.00000000000000]'),
    ('0.5**x - 4', '[-2.00000000000000]'),
    ('1.5**x - 2', '[1.70951129135145]'),
    ('exp(0.5*x) - 2', '[1.38629436111989]'),
    ('2.5*exp(x) - 1', '[-0.916290731874155]'),
    ('log(2.0*x) - 1', '[1.35914091422952]'),
    ('log(x) - 0.5', '[1.64872127070013]'),
    ('x**0.5 - 2', '[4.00000000000000]'),
    ('x**1.5 - 8', '[4.00000000000000]'),
    ('x**2 - 0.25', '[-0.500000000000000, 0.500000000000000]'),
    ('1.0*x**3 - 1', '[1.00000000000000, -0.5 - 0.866025403784439*I, -0.5 + 0.866025403784439*I]'),
]


@pytest.mark.parametrize('source, expected', FLOATS)
def test_float_equations_are_rationalized_and_floated_back(source, expected):
    assert solve(source) == expected


def test_rational_true_keeps_the_exact_solutions():
    assert solve('x**2 - 2.0', rational=True) == '[-sqrt(2), sqrt(2)]'


def test_rational_false_keeps_the_floats():
    expected = ('[1.25992104989487, -0.629960524947437 - 1.09112363597172*I, '
                '-0.629960524947437 + 1.09112363597172*I]')
    assert solve('x**3 - 2.0', rational=False) == expected


def test_float_rationalization_applies_to_dict_output():
    assert solve('4.0**x - 8', dict=True) == (
        '[{x: 1.50000000000000}, {x: 1.5 + 4.53236014182719*I}]')


HYPERBOLIC = [
    ('sinh(x) - 1', '[log(-1 + sqrt(2)) + I*pi, log(1 + sqrt(2))]'),
    ('cosh(x) - 2', '[log(2 - sqrt(3)), log(sqrt(3) + 2)]'),
    ('tanh(x) - S(1)/2', '[log(sqrt(3)) + I*pi, log(sqrt(3))]'),
    ('sinh(x) - cosh(x) - 1', '[I*pi]'),
    ('cosh(x) - 1', '[0]'),
    ('sinh(2*x) - 1', '[log(sqrt(-1 + sqrt(2))) - I*pi/2, log(sqrt(-1 + sqrt(2))) + I*pi/2, '
                      'log(sqrt(1 + sqrt(2))) + I*pi, log(sqrt(1 + sqrt(2)))]'),
    ('sinh(x)**2 - 1', '[log(-1 + sqrt(2)) + I*pi, log(1 + sqrt(2)) + I*pi, log(-1 + sqrt(2)), log(1 + sqrt(2))]'),
    ('x - sinh(a)', '[sinh(a)]'),
    ('cosh(x)*x', '[0, -I*pi/2, I*pi/2]'),
    ('sech(x) - S(1)/2', '[log(2 - sqrt(3)), log(sqrt(3) + 2)]'),
    ('coth(x) - 2', '[log(sqrt(3)) + I*pi, log(sqrt(3))]'),
    ('Eq(sinh(x), 3)', '[log(-3 + sqrt(10)) + I*pi, log(3 + sqrt(10))]'),
    ('asinh(x) - 1', '[sinh(1)]'),
    ('sinh(x) - a', '[log(a - sqrt(a**2 + 1)), log(a + sqrt(a**2 + 1))]'),
]


@pytest.mark.parametrize('source, expected', HYPERBOLIC)
def test_hyperbolic_equations_are_rewritten_through_exp(source, expected):
    assert solve(source) == expected


def test_hyperbolic_equation_without_closed_form_is_unsupported():
    with pytest.raises(NotImplementedError):
        solve('sinh(x) + x')


NEGATIVE_EXPONENTIALS = [
    ('exp(-x) + 1', '[I*pi]'),
    ('exp(-x) - 2', '[-log(2)]'),
    ('exp(-2*x) + 1', '[-I*pi/2, I*pi/2]'),
    ('exp(-x/2) + 1', '[2*I*pi]'),
    ('exp(1 - x) + 2', '[-log(2) + 1 + I*pi]'),
    ('exp(-x) + exp(-2*x) - 1', '[log(-1/2 + sqrt(5)/2) + I*pi, log(1/2 + sqrt(5)/2)]'),
    ('3*exp(-x) + exp(x) - 4', '[0, log(3)]'),
    ('x*exp(-x) - 1', '[-LambertW(-1)]'),
    ('exp(-x) - x', '[LambertW(1)]'),
]


@pytest.mark.parametrize('source, expected', NEGATIVE_EXPONENTIALS)
def test_exponentials_of_negative_arguments_use_the_reference_generator(source, expected):
    assert solve(source) == expected


def test_log_sum_keeps_only_the_checked_root():
    assert solve('log(x) + log(x - 1)') == '[1/2 + sqrt(5)/2]'


def test_solutions_are_filtered_by_the_sign_assumptions():
    x = sx.Symbol('x', negative=True)
    y = sx.Symbol('y', positive=True)
    assert sx.solve(x - sx.log(2), x) == []
    assert sx.solve(x**2 - 4, x) == [-2]
    assert sx.solve(sx.exp(x) - 2, x) == []
    assert sx.solve(y + sx.log(2), y) == []
    assert sx.solve(y - sx.log(sx.S(1)/3), y) == []
    assert str(sx.solve(y**2 - sx.log(2), y, dict=True)) == '[{y: sqrt(log(2))}]'
    assert sx.solve(sx.exp(-y) + 1, y) == []
    assert sx.solve(sx.exp(-x) - sx.S(1)/2, x) == []


MIXED_BASES = [
    ('3**x - 2**x', '[0]'),
    ('2**x*3**x - 6', '[1]'),
    ('exp(x) - 2**x', '[0]'),
    ('exp(x)*2**x - 3', '[log(3)/(log(2) + 1)]'),
    ('2**x - 3**(x + 1)', '[log(3)/(-log(3) + log(2))]'),
]


@pytest.mark.parametrize('source, expected', MIXED_BASES)
def test_powers_of_different_bases_are_combined(source, expected):
    assert solve(source) == expected


LOGARITHMIC_LAMBERT = [
    ('x**x - 4', '[2]'),
    ('x**x - 3', '[exp(LambertW(log(3)))]'),
    ('x**x - 27', '[3]'),
    ('x**x - a', '[exp(LambertW(log(a)))]'),
    ('x**x - 1', '[0, 1]'),
    ('x**x', '[]'),
    ('x**x - S(1)/4', '[exp(LambertW(-log(4)))]'),
    ('x**(x + 1) - 8', '[2]'),
    ('x**x - 2', '[exp(LambertW(log(2)))]'),
    ('2*x*log(x) - 1', '[exp(LambertW(1/2))]'),
    ('x*log(x) - a', '[exp(LambertW(a))]'),
    ('x*log(x) - 1', '[exp(LambertW(1))]'),
    ('(x + 1)*log(x + 1) - 2', '[-1 + exp(LambertW(2))]'),
    ('x*log(2*x) - 1', '[exp(LambertW(2))/2]'),
    ('x*log(x) + S(1)/4', '[exp(LambertW(-1/4)), exp(LambertW(-1/4, -1))]'),
    ('x**2*log(x) - 1', '[exp(LambertW(2)/2)]'),
    ('x*log(x) - log(4)', '[2]'),
    ('x*log(x) - log(12)', '[exp(LambertW(log(12)))]'),
    ('x*log(x) + log(4)', '[exp(LambertW(-log(4)))]'),
    ('x*log(x) - 3*exp(3)', '[exp(3)]'),
    ('x*exp(x) - 2*exp(2)', '[2]'),
    ('x*exp(x) - exp(1)', '[1]'),
    ('x*exp(x) - sqrt(2)', '[LambertW(sqrt(2))]'),
]


@pytest.mark.parametrize('source, expected', LOGARITHMIC_LAMBERT)
def test_logarithmic_lambert_forms_give_the_reference_solutions(source, expected):
    assert solve(source) == expected
