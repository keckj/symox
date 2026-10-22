"""Solve forms against hardcoded SymPy 1.14 outputs: Lambert equations in a
non-E base, exponentials of shifted arguments, the canonical key order of
solution dicts, the _invert preprocessing of systems and the Float components
of systems solved with rational=False. The real W_-1 branches that SymPy
misses are kept as intentional corrections and checked numerically."""

import pytest

import symox as sx


NAMES = 'x y z a b c'


def namespace():
    names = dict(zip(NAMES.split(), sx.symbols(NAMES)))
    for name in ('exp', 'log', 'sqrt', 'S', 'Rational', 'sin', 'solve', 'I', 'pi'):
        names[name] = getattr(sx, name)
    return names


def run(source):
    return str(eval(source, namespace()))


LAMBERT = [
    ('solve(x*2**x + Rational(1, 4), x)', '[-4, LambertW(-log(2)/4)/log(2)]'),
    ('solve(x*2**x - 1, x)', '[LambertW(log(2))/log(2)]'),
    ('solve(x - 2**(-x), x)', '[LambertW(log(2))/log(2)]'),
    ('solve(x**2 - 2**x, x)', '[2, 4, -2*LambertW(log(2)/2)/log(2)]'),
    ('solve(x**(2*x) - 4, x)', '[exp(LambertW(log(2)))]'),
    ('solve(x*log(x) - log(9)/2, x)', '[exp(LambertW(log(3)))]'),
    ('solve(x*log(x) + 2*log(2), x)', '[exp(LambertW(-log(4)))]'),
    ('solve(x**x - sqrt(2), x)', '[exp(LambertW(log(2)/2))]'),
    ('solve(exp(x) - exp(2*x)/3 - 1, x)', '[log(3)/2 - I*pi/6, log(3)/2 + I*pi/6]'),
]

SHIFTED_EXP = [
    ('solve(exp(a - x) + 1, x)', '[log(-exp(a))]'),
    ('solve(exp(a + x) + 1, x)', '[log(-exp(-a))]'),
    ('solve(exp(x - a) + 1, x)', '[log(-exp(a))]'),
    ('solve(exp(a - x) - 1, x)', '[log(exp(a))]'),
    ('solve(exp(a - x) - 2, x)', '[log(exp(a)/2)]'),
    ('solve(exp(a - x) + 2, x)', '[log(-exp(a)/2)]'),
    ('solve(exp(2*a - x) + 1, x)', '[log(-exp(2*a))]'),
    ('solve(exp(-a - x) + 3, x)', '[log(-exp(-a)/3)]'),
    ('solve(exp(a*x + b) - c, x)', '[log(c*exp(-b))/a]'),
]

KEY_ORDER = [
    ('solve([x + y - 2, x - y], [y, x], dict=True)', '[{x: 1, y: 1}]'),
    ('solve([x + y + z - 1, x - y, y - z], [z, y, x], dict=True)', '[{x: 1/3, y: 1/3, z: 1/3}]'),
    ('solve([x*y - 1, x + y - 2], [y, x], dict=True)', '[{x: 1, y: 1}]'),
    ('solve([x**2 - y, y - 4], [y, x], dict=True)', '[{x: -2, y: 4}, {x: 2, y: 4}]'),
    ('solve([y - x**2, x + y - 2], [y, x])', '[(4, -2), (1, 1)]'),
]

INVERTED_SYSTEMS = [
    ('solve([exp(2*x) - 1], x)', '[(0,)]'),
    ('solve([exp(x) - 2], x)', '[(log(2),)]'),
    ('solve([exp(2*x) - 1, y - 1], [x, y])', '[(0, 1)]'),
    ('solve([exp(x + y) - 2, x - y], [x, y])', '{x: log(2)/2, y: log(2)/2}'),
    ('solve([exp(x*y) - 1, x - 1], [x, y])', '[(1, 0)]'),
    ('solve([log(x + y) - 1, x - y], [x, y])', '{x: E/2, y: E/2}'),
    ('solve([exp(x) - exp(y), x + y - 2], [x, y])', '{x: 1, y: 1}'),
    ('solve([log(x) - log(y), x + y - 2], [x, y])', '{x: 1, y: 1}'),
    ('solve([sqrt(x + y) - 2, x - y], [x, y])', '{x: 2, y: 2}'),
    ('solve([exp(2*x + y) - 3, x - y], [x, y])', '{x: log(3)/3, y: log(3)/3}'),
    ('solve([exp(x) + exp(y) - 2, x - y], [x, y])', '[(0, 0)]'),
    ('solve([sin(x + y) - 1, x - y], [x, y])', '[(pi/4, pi/4)]'),
    ('solve([2*exp(x - y) - 4, x + y], [x, y])', '{x: log(2)/2, y: -log(2)/2}'),
    ('solve([exp(x)*exp(y) - 1, x - y], [x, y])', '{x: 0, y: 0}'),
    ('solve([exp(x + y) - 1, x - y])', '{x: 0, y: 0}'),
    ('solve([exp(x + y) - 1, x - y], [x, y])', '{x: 0, y: 0}'),
    ('solve([exp(x + y) - a, x - y], [x, y])', '{x: log(a)/2, y: log(a)/2}'),
    ('solve([exp(x**2) - 1, y - x], [x, y])', '[(0, 0)]'),
    ('solve([log(x*y) - 1, x - 2], [x, y])', '[(2, E/2)]'),
    ('solve([exp(-x - y) - 1, x - 2*y], [x, y])', '{x: 0, y: 0}'),
    ('solve([1/exp(x + y) - 2, x - y], [x, y])', '{x: -log(2)/2, y: -log(2)/2}'),
]

FLOAT_COMPONENTS = [
    ('solve([0.25*x**2 - 1, x - 2], x, rational=False)', '[(2.00000000000000,)]'),
    ('solve([x - 1, 0.5*y - x], [x, y], rational=False)', '{x: 1.00000000000000, y: 2.00000000000000}'),
]


@pytest.mark.parametrize('source, expected', LAMBERT)
def test_lambert_forms_in_a_non_e_base_give_the_reference_solutions(source, expected):
    assert run(source) == expected


@pytest.mark.parametrize('source, expected', SHIFTED_EXP)
def test_exponentials_of_shifted_arguments_give_the_reference_solutions(source, expected):
    assert run(source) == expected


@pytest.mark.parametrize('source, expected', KEY_ORDER)
def test_solution_dicts_use_the_reference_key_order(source, expected):
    assert run(source) == expected


def test_solution_dict_keys_are_in_canonical_order():
    result = eval('solve([x + y + z - 1, x - y, y - z], [z, y, x], dict=True)', namespace())
    assert [str(k) for k in result[0]] == ['x', 'y', 'z']


@pytest.mark.parametrize('source, expected', INVERTED_SYSTEMS)
def test_systems_are_inverted_like_the_reference(source, expected):
    assert run(source) == expected


@pytest.mark.parametrize('source, expected', FLOAT_COMPONENTS)
def test_float_components_stay_float_with_rational_false(source, expected):
    assert run(source) == expected


# SymPy 1.14 returns only the principal Lambert branch for these equations;
# the extra real solutions come from the W_-1 branch and are kept.
EXTRA_REAL_BRANCH = [
    ('2**x - 4*x', '4'),
    ('3**x - 3*x', '1'),
    ('2**x - x - 2', '2'),
]


@pytest.mark.parametrize('source, extra', EXTRA_REAL_BRANCH)
def test_real_lambert_branch_missed_by_the_reference_is_kept(source, extra):
    names = namespace()
    equation = eval(source, names)
    roots = sx.solve(equation, names['x'])
    assert extra in [str(root) for root in roots]
    for root in roots:
        assert abs(complex(sx.N(equation.subs(names['x'], root), 30))) < 1e-20


@pytest.mark.parametrize('source', ['exp(x) - x - 2', 'exp(x) - 3*x'])
def test_both_real_lambert_branches_are_valid_solutions(source):
    names = namespace()
    equation = eval(source, names)
    roots = sx.solve(equation, names['x'])
    assert len(roots) == 2
    for root in roots:
        assert abs(complex(sx.N(equation.subs(names['x'], root), 30))) < 1e-20
