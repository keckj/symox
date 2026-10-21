"""Solution sets of solve against hardcoded SymPy 1.14 outputs: equations in
powers of a non-E base and equations with Float coefficients."""

import pytest

import symox as sx


NAMES = 'x a b'


def solve(source, **flags):
    namespace = dict(zip(NAMES.split(), sx.symbols(NAMES)))
    for name in ('exp', 'log', 'sqrt', 'S', 'Rational'):
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
