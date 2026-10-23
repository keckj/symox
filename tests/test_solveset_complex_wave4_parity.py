"""Native solveset over the complex domain against hardcoded SymPy 1.14
outputs: exponential polynomials inverted into ImageSets, trigonometric and
hyperbolic equations rewritten in exp, and the systems that inherit them.
Where SymPy prints the same set in another form, the sets are compared through
their points in a bounded window."""

import cmath

import pytest

import symox as sx


x, y, a = sx.symbols('x y a')
p = sx.symbols('p', positive=True)
r = sx.symbols('r', real=True)
DOMAINS = {'Complexes': sx.S.Complexes, 'Reals': sx.S.Reals,
           'Interval(0, oo)': sx.Interval(0, sx.oo), 'Interval(0, 10)': sx.Interval(0, 10)}


def parse(source):
    namespace = {name: getattr(sx, name) for name in dir(sx) if not name.startswith('_')}
    namespace.update(x=x, y=y, a=a, p=p, r=r, _n=sx.Symbol('_n'))
    return eval(source, namespace)


def solveset(source, domain):
    return sx.solveset(parse(source), x, DOMAINS[domain])


def points(s, bound=12.0):
    """The complex points of a set built from ImageSets over the integers."""
    if isinstance(s, sx.Union):
        return set().union(*(points(arg, bound) for arg in s.args))
    if isinstance(s, sx.Complement):
        return points(s.args[0], bound) - points(s.args[1], bound)
    if s == sx.S.EmptySet:
        return set()
    if isinstance(s, sx.FiniteSet):
        values = [complex(sx.N(v, 30)) for v in s.args]
    else:
        assert isinstance(s, sx.ImageSet) and s.base_set == sx.S.Integers, s
        values = [complex(sx.N(s.lamda(n), 30)) for n in range(-40, 41)]
    return {(round(v.real, 8), round(v.imag, 8)) for v in values if abs(v) <= bound}


EXACT = [
    ('exp(x) - exp(-x)', 'Complexes', 'ImageSet(Lambda(_n, _n*I*pi), Integers)'),
    ('exp(2*x) - 3*exp(x) + 2', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers))'),
    ('sin(x) + cos(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi + 3*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + 7*pi/4), Integers))'),
    ('exp(x) - 1', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'),
    ('exp(x) + 1', 'Complexes', 'ImageSet(Lambda(_n, I*(2*_n*pi + pi)), Integers)'),
    ('exp(2*x) - 1', 'Complexes', 'ImageSet(Lambda(_n, _n*I*pi), Integers)'),
    ('exp(x)**2 - 4', 'Complexes', 'ImageSet(Lambda(_n, _n*I*pi + log(4)/2), Integers)'),
    ('exp(3*x) - 1', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi/3), Integers)'),
    ('exp(x) - exp(-x)', 'Reals', '{0}'),
    ('exp(2*x) - 3*exp(x) + 2', 'Reals', '{0, log(2)}'),
    ('2**x - 8', 'Reals', '{3}'),
    ('sin(x) + cos(x)', 'Reals', 'Union(ImageSet(Lambda(_n, 2*_n*pi + 3*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + 7*pi/4), Integers))'),
    ('exp(x) - a', 'Complexes', 'ImageSet(Lambda(_n, I*(2*_n*pi + arg(a)) + log(Abs(a))), Integers)'),
    ('exp(a*x) - 2', 'Complexes', 'ImageSet(Lambda(_n, (2*_n*I*pi + log(2))/a), Integers)'),
    ('exp(x/2) - 3', 'Complexes', 'ImageSet(Lambda(_n, 4*_n*I*pi + 2*log(3)), Integers)'),
    ('exp(x) - y', 'Complexes', 'ImageSet(Lambda(_n, I*(2*_n*pi + arg(y)) + log(Abs(y))), Integers)'),
    ('exp(x)*(x-1)', 'Complexes', '{1}'),
    ('x*(exp(x)-1)', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'),
    ('(exp(x)-1)/x', 'Complexes', 'Complement(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), {0})'),
    ('exp(x)/(exp(x)-2) - 3', 'Complexes', 'Complement(ImageSet(Lambda(_n, 2*_n*I*pi + log(3)), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers))'),
    ('sinh(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + I*pi), Integers))'),
    ('cosh(x) - 1', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'),
    ('sinh(x) - 1', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*I*pi + log(1 + sqrt(2))), Integers), ImageSet(Lambda(_n, 2*_n*I*pi - log(1 + sqrt(2)) + I*pi), Integers))'),
    ('tanh(x) - 2', 'Complexes', 'ImageSet(Lambda(_n, _n*I*pi + atanh(2)), Integers)'),
    ('sinh(x) + cosh(x) - 2', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers)'),
    ('sin(x)**2 - cos(x)**2', 'Complexes', 'Union(ImageSet(Lambda(_n, _n*pi + pi/4), Integers), ImageSet(Lambda(_n, _n*pi + 3*pi/4), Integers))'),
    ('sin(x) + cos(x) - 1', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers))'),
    ('tan(x) - 1', 'Complexes', 'ImageSet(Lambda(_n, _n*pi + pi/4), Integers)'),
    ('sin(x) - I', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi + I*log(1 + sqrt(2))), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi - I*log(1 + sqrt(2))), Integers))'),
    ('sin(x) + 2*cos(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi - atan(2) + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi - atan(2) + 2*pi), Integers))'),
    ('sin(x) + cos(x) - 1', 'Reals', 'Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers))'),
    ('log(x) - 1', 'Complexes', '{E}'),
    ('log(x) - y', 'Complexes', '{exp(y)}'),
    ('log(x - 1) - 2', 'Complexes', '{1 + exp(2)}'),
    ('exp(x) - x', 'Complexes', 'ConditionSet(x, Eq(-x + exp(x), 0), Complexes)'),
    ('exp(x) + x', 'Complexes', 'ConditionSet(x, Eq(x + exp(x), 0), Complexes)'),
    ('x**2*exp(x) - 1', 'Complexes', 'ConditionSet(x, Eq(x**2*exp(x) - 1, 0), Complexes)'),
    ('exp(I*x) - 1', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*pi), Integers)'),
    ('exp(I*x) + exp(-I*x)', 'Complexes', 'ImageSet(Lambda(_n, _n*pi + pi/2), Integers)'),
    ('exp(x) + exp(-x) - 2', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'),
    ('exp(x) + exp(-x)', 'Complexes', 'ImageSet(Lambda(_n, I*(2*_n*pi + pi)/2), Integers)'),
    ('exp(2*x) + exp(x) + 1', 'Complexes', 'Union(ImageSet(Lambda(_n, I*(2*_n*pi - 2*pi/3)), Integers), ImageSet(Lambda(_n, I*(2*_n*pi + 2*pi/3)), Integers))'),
    ('exp(4*x) - 1', 'Complexes', 'ImageSet(Lambda(_n, _n*I*pi/2), Integers)'),
    ('exp(x/3) - 1', 'Complexes', 'ImageSet(Lambda(_n, 6*_n*I*pi), Integers)'),
    ('exp(x+1) - 2', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi - 1 + log(2)), Integers)'),
    ('2*exp(x) - 3', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi + log(3/2)), Integers)'),
    ('exp(x) - exp(1)', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi + 1), Integers)'),
    ('exp(-x) - 2', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi - log(2)), Integers)'),
    ('exp(p*x) - 1', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi/p), Integers)'),
    ('exp(x) - p', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi + log(p)), Integers)'),
    ('exp(x) - r', 'Complexes', 'ImageSet(Lambda(_n, I*(2*_n*pi + arg(r)) + log(Abs(r))), Integers)'),
    ('exp(x) - exp(-x)', 'Interval(0, oo)', '{0}'),
    ('x**x - 1', 'Complexes', 'ConditionSet(x, Eq(x**x - 1, 0), Complexes)'),
    ('2**x - x', 'Complexes', 'ConditionSet(x, Eq(2**x - x, 0), Complexes)'),
    ('exp(x)-1', 'Interval(0, 10)', '{0}'),
    ('sinh(x) - cosh(x)', 'Complexes', 'EmptySet'),
    ('exp(x)*sin(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))'),
    ('(exp(x)-1)*(exp(x)-2)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers))'),
    ('exp(x)**3 - 8', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi/3 + log(8)/3), Integers)'),
    ('E**x - 1', 'Complexes', 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'),
    ('exp(2*x) - 4', 'Reals', '{log(4)/2}'),
]


@pytest.mark.parametrize('source,domain,expected', EXACT)
def test_same_set_as_reference(source, domain, expected):
    assert str(solveset(source, domain)) == expected


SAME_SET = [
    ('sin(x) - cos(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi + 5*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/4), Integers))'),
    ('sin(x) - cos(x)', 'Reals', 'Union(ImageSet(Lambda(_n, 2*_n*pi + 5*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/4), Integers))'),
    ('sin(x)*cos(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))'),
    ('sin(x)*cos(x)', 'Reals', 'Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))'),
    ('sin(2*x) - sin(x)', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/3), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/3), Integers))'),
    ('sin(2*x) - sin(x)', 'Reals', 'Union(ImageSet(Lambda(_n, 2*_n*pi + pi), Integers), ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/3), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/3), Integers))'),
    ('cos(x) - 2', 'Complexes', 'Union(ImageSet(Lambda(_n, 2*_n*pi - acos(2)), Integers), ImageSet(Lambda(_n, 2*_n*pi + acos(2)), Integers))'),
    ('exp(x)/(exp(x) - 2) - 3', 'Complexes', 'Complement(ImageSet(Lambda(_n, 2*_n*I*pi + log(3)), Integers), ImageSet(Lambda(_n, 2*_n*I*pi + log(2)), Integers))'),
]


@pytest.mark.parametrize('source,domain,expected', SAME_SET)
def test_same_set_in_another_form(source, domain, expected):
    got, want = points(solveset(source, domain)), points(parse(expected))
    assert want and got == want


# SymPy cannot invert a rational base and returns a ConditionSet, or an
# incomplete set for mixed bases; symox inverts b**x = exp(x*log(b)) completely
DIVERGENT = [
    ('2**x - 8', 'ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + 3), Integers)'),
    ('3**x - 9', 'ImageSet(Lambda(_n, 2*_n*I*pi/log(3) + 2), Integers)'),
    ('2**(2*x) - 2**x - 2', 'Union(ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + 1), Integers), ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + I*pi/log(2)), Integers))'),
    ('4**x - 2**x - 2', 'Union(ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + 1), Integers), ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + I*pi/log(2)), Integers))'),
    # SymPy adds a spurious nan
    ('exp(2*x) - exp(x)', 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'),
    # SymPy keeps the root 1/2 - sqrt(5)/2 where log(x) + log(x - 1) is 2*I*pi + log(-1)
    ('log(x) + log(x - 1)', '{1/2 + sqrt(5)/2}'),
]


@pytest.mark.parametrize('source,expected', DIVERGENT)
def test_complete_where_reference_is_not(source, expected):
    got = solveset(source, 'Complexes')
    assert str(got) == expected
    f = parse(source)
    for re_, im_ in points(got, 30.0):
        assert abs(complex(sx.N(f.subs(x, complex(re_, im_)), 20))) < 1e-6


@pytest.mark.parametrize('source', ['2**x - 3**x', '2**x*3**(2*x) - 6', 'x**x - 1',
                                    'exp(x) - x', '2**x - x'])
def test_undecided_is_condition_set(source):
    # never an incomplete finite set: 2**x = 3**x also holds at 2*I*pi*n/log(3/2)
    assert isinstance(solveset(source, 'Complexes'), sx.ConditionSet)


def test_combined_numeric_bases_are_one_base():
    # 2**x*3**x is 6**x, as in SymPy; 6**x - 6 is then solved through
    # exp(x*log(6)) (documented divergence: SymPy returns a ConditionSet)
    assert str(solveset('2**x*3**x - 6', 'Complexes')) == \
        'ImageSet(Lambda(_n, 2*_n*I*pi/log(6) + 1), Integers)'


SYSTEMS = [
    ('[exp(x) - sin(y), y**2 - 4]', '{(ImageSet(Lambda(_n, I*(2*_n*pi + pi) + log(sin(2))), Integers), -2), (ImageSet(Lambda(_n, 2*_n*I*pi + log(sin(2))), Integers), 2)}'),
    ('[exp(y) - x, x**2 + 1]', '{(-I, ImageSet(Lambda(_n, I*(2*_n*pi - pi/2)), Integers)), (I, ImageSet(Lambda(_n, I*(2*_n*pi + pi/2)), Integers))}'),
    ('[sin(x)**2 - 1, y - 1]', '{(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), 1), (ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers), 1)}'),
    ('[2**x - 8, y - x]', '{(3, 3), (ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + 3), Integers), ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + 3), Integers))}'),
]


@pytest.mark.parametrize('system,expected', SYSTEMS)
def test_systems(system, expected):
    from symox.solvers.solveset import substitution
    equations = parse(system)
    # the first entry matches SymPy; SymPy returns the unsolved
    # {(exp(y), y)}, {(x, 1)} and {(y, y)} for the others
    assert str(sx.nonlinsolve(equations, [x, y])) == expected
    assert str(substitution(equations, [x, y])) == expected


@pytest.mark.parametrize('source,expected', [
    ('exp(x) - 2', '[log(2)]'),
    ('exp(2*x) - 3*exp(x) + 2', '[0, log(2)]'),
    ('2**x - 8', '[3]'),
    ('exp(x) - exp(-x)', '[0, I*pi]'),
    ('sinh(x) - 1', '[log(-1 + sqrt(2)) + I*pi, log(1 + sqrt(2))]'),
])
def test_solve_unchanged(source, expected):
    assert str(sx.solve(parse(source), x)) == expected


def test_order_subs_takes_the_limit_of_the_first_solution():
    # Order._eval_subs takes the limit of log(Abs(x)) + I*arg(x), the first
    # point of the complete ImageSet; arg(x) -> 0 as x -> 0+
    assert str(sx.O(x**3).subs(x, sx.exp(x))) == 'O(exp(3*x), (x, -oo))'
    assert [sx.limit(sx.arg(x), x, 0), sx.limit(sx.arg(x), x, 0, '-'), sx.limit(sx.arg(x - 1), x, 1, '-')] == [0, sx.pi, sx.pi]
