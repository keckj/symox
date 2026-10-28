"""A nonfinite-parameter fallback must survive partial substitution.

The protected predicate explicitly permits nonfinite parameters. Finite-symbol
facts such as exp(z)!=0 therefore cannot decide its held equation before the
parameter guard is known. The ordinary relation protocol stays unchanged.
"""
import itertools
import pytest
import symox as sx


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals, sx.Interval(-2, 3), sx.FiniteSet(-1, 0, 1)])
@pytest.mark.parametrize('target', [0, 1, 2])
@pytest.mark.parametrize('shift', [-sx.oo, sx.oo, sx.zoo, sx.nan, 0, 1, -1, sx.I])
def test_every_parameter_order_agrees_with_simultaneous_specialization(domain, target, shift):
    x, a, y = sx.symbols('x a y')
    equation = sx.exp(sx.exp(x)+a)-y
    result = sx.solveset(equation, x, domain)
    substitutions = {y: target, a: shift}
    expected = result.subs(substitutions, simultaneous=True)
    for pairs in itertools.permutations(substitutions.items()):
        actual = result
        for parameter, value in pairs:
            actual = actual.subs(parameter, value)
        assert actual == expected
    assert result.subs(y, target, simultaneous=True).subs(a, shift) == expected
    if shift == -sx.oo:
        # Every x in these domains is finite. exp(x) is finite, so the
        # outer exponent is -oo and its exponential is exactly zero.
        assert expected == (domain if target == 0 else sx.EmptySet)
    elif shift in (0, 1, -1, sx.I) and target == 0:
        # A finite exponential has no zero; this is an absence proof.
        assert expected == sx.EmptySet


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
@pytest.mark.parametrize('coefficient', [-sx.oo, sx.oo])
def test_nonfinite_multiplier_preserves_unsupported_half_lines(domain, coefficient):
    x, a, y = sx.symbols('x a y')
    equation = sx.exp(sx.exp(a*x))-y
    result = sx.solveset(equation, x, domain)
    for target in (0, 1, 2):
        expected = result.subs({a: coefficient, y: target}, simultaneous=True)
        assert result.subs(y, target).subs(a, coefficient) == expected
        assert result.subs(a, coefficient).subs(y, target) == expected
        if target == 1:
            # The unsupported original equation includes a real half-line.
            assert expected != sx.EmptySet


def test_multiple_nonfinite_guards_preserve_all_held_equations():
    x, a, b, y = sx.symbols('x a b y')
    guard = ~(sx.Contains(a, sx.Complexes) & sx.Contains(b, sx.Complexes))
    first = sx.Eq(sx.exp(a+sx.exp(x)), y, evaluate=False)
    second = sx.Eq(sx.exp(b+sx.exp(x)), y, evaluate=False)
    result = sx.ConditionSet(x, guard & first & second, sx.Reals)
    partial = result.subs(y, 0)
    assert partial != sx.EmptySet
    assert len(partial.condition.atoms(sx.Equality)) == 2
    assert partial.subs(a, -sx.oo).subs(b, -sx.oo) == sx.Reals
    assert partial.subs(a, -sx.oo).subs(b, 0) == sx.EmptySet
    assert partial.subs(a, 0).subs(b, -sx.oo) == sx.EmptySet


def test_xreplace_and_simultaneous_partial_substitution_share_the_guard_context():
    x, a, y = sx.symbols('x a y')
    result = sx.solveset(sx.exp(sx.exp(x)+a)-y, x)
    for partial in [result.xreplace({y: sx.S.Zero}), result.subs(y, 0, simultaneous=True)]:
        assert partial != sx.EmptySet
        assert partial.subs(a, -sx.oo) == sx.Complexes
        assert partial.xreplace({a: -sx.oo}) == sx.Complexes


def test_membership_and_bound_symbol_substitution_keep_the_predicate():
    x, a, y = sx.symbols('x a y')
    result = sx.solveset(sx.exp(sx.exp(x)+a)-y, x).subs(y, 0)
    assert result.subs(x, 2) == result  # the ConditionSet binder is protected
    membership = result.contains(0)
    assert membership.subs(a, -sx.oo) == sx.true
    assert membership.subs(a, 0) == sx.false


def test_guarded_negation_and_piecewise_condition_do_not_evaluate_the_equation_early():
    x, a, y = sx.symbols('x a y')
    guard = ~sx.Contains(a, sx.Complexes)
    held = sx.Eq(sx.exp(a+sx.exp(x)), y, evaluate=False)
    condition = guard & ~held
    partial = condition.subs(y, 0)
    assert partial.subs(a, -sx.oo) == sx.false
    piece = sx.Piecewise((sx.S.One, guard & held), (sx.S.Zero, True))
    assert piece.subs(y, 0).subs(a, -sx.oo) == 1
    assert piece.subs(y, 0).subs(a, 0) == 0


def test_ordinary_equations_and_finite_domain_predicates_keep_their_protocol():
    x, a, y = sx.symbols('x a y')
    held = sx.Eq(sx.exp(a+sx.exp(x)), y, evaluate=False)
    assert held.subs(y, 0) == sx.false
    assert sx.Eq(sx.nan, 0) == sx.false
    finite = sx.Contains(a, sx.Complexes) & held
    assert finite.subs(y, 0) == sx.false
    guard = ~sx.Contains(a, sx.Complexes)
    condition = guard & held
    assert condition.subs(held, sx.false) == sx.false
    assert condition.subs(held, sx.true) == guard


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals, sx.FiniteSet(-1, 0, 1)])
@pytest.mark.parametrize('values,undefined', [
    ((-sx.oo, 0), False), ((-sx.oo, 1), False),
    ((-sx.oo, sx.I), False), ((-sx.oo, sx.oo), True),
    ((sx.oo, -sx.oo), True), ((sx.nan, 0), True),
])
def test_two_shifts_retain_arithmetic_and_guard_in_every_order(domain, values, undefined):
    x, a, b, y = sx.symbols('x a b y')
    result = sx.solveset(sx.exp(sx.exp(x)+a+b)-y, x, domain)
    assignments = [(a, values[0]), (b, values[1]), (y, 0)]
    expected = result.subs(dict(assignments), simultaneous=True)
    for pairs in itertools.permutations(assignments):
        actual = result
        for parameter, value in pairs:
            actual = actual.subs(parameter, value)
        assert actual == expected
    if undefined:
        assert expected != sx.EmptySet
        assert isinstance(expected, sx.ConditionSet)
        assert expected.condition.has(sx.nan)
    else:
        assert expected == domain


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
@pytest.mark.parametrize('values', [(-sx.oo, 0, 1), (-sx.oo, sx.oo, 0), (0, -sx.oo, sx.oo)])
def test_three_shifts_do_not_erase_a_guard_that_has_become_true(domain, values):
    x, a, b, c, y = sx.symbols('x a b c y')
    result = sx.solveset(sx.exp(sx.exp(x)+a+b+c)-y, x, domain)
    assignments = [(a, values[0]), (b, values[1]), (c, values[2]), (y, 0)]
    expected = result.subs(dict(assignments), simultaneous=True)
    for pairs in itertools.permutations(assignments):
        actual = result
        for parameter, value in pairs:
            actual = actual.subs(parameter, value)
        assert actual == expected
    assert expected == domain if sx.oo not in values else expected != sx.EmptySet


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
def test_zero_multiplier_stays_conditional_until_its_inner_parameters_decide(domain):
    x, a, b, c, y = sx.symbols('x a b c y')
    result = sx.solveset(sx.exp(c*sx.exp(a*x)+b)-y, x, domain)
    partial = result.subs({c: 0, b: 0, y: 1}, simultaneous=True)
    assert partial.subs(a, 1) == domain
    specialized = partial.subs(a, sx.oo)
    assert specialized != sx.EmptySet
    assert specialized != domain
    assert specialized.contains(-1) == sx.true
    # exp(oo) is infinite, so 0*exp(oo) is undefined at x=1.
    assert specialized.contains(1) != sx.true
    assert result.subs({c: 0, b: 0, y: 1, a: sx.oo}, simultaneous=True) == specialized


def test_closed_unknown_function_values_keep_their_nonfinite_fallback():
    x, a, y = sx.symbols('x a y')
    f = sx.Function('f')
    result = sx.solveset(sx.exp(sx.exp(x)+f(a))-y, x)
    partial = result.subs(y, 0).subs(a, 1)
    assert partial != sx.EmptySet
    assert partial.subs(f(1), -sx.oo) == sx.Complexes
    assert partial.subs(f(1), 0) == sx.EmptySet


def test_guarded_substitution_preserves_a_function_host_protocol():
    calls = []
    x, a, y = sx.symbols('x a y')

    class H(sx.Function):
        def _eval_subs(self, old, new):
            calls.append((old, new))
            if old == a:
                return new+1

    result = sx.solveset(sx.exp(sx.exp(x)+H(a))-y, x)
    partial = result.subs(y, 0)
    calls.clear()
    assert partial.subs(a, 1) == sx.EmptySet
    assert (a, sx.S.One) in calls


def test_a_disjunctive_guard_keeps_the_ordinary_relation_protocol():
    a, b, y = sx.symbols('a b y')
    guard = ~(sx.Contains(a, sx.Complexes) & sx.Contains(b, sx.Complexes))
    condition = guard | sx.Eq(sx.exp(a+b), y, evaluate=False)
    assert condition.subs(y, 0).subs(a, sx.nan) == sx.true


@pytest.mark.parametrize('coefficient,member,nonmember', [(sx.oo, 1, -1), (-sx.oo, -1, 1)])
def test_negative_inner_exponential_can_make_the_outer_exponential_zero(coefficient, member, nonmember):
    x, a, y = sx.symbols('x a y')
    result = sx.solveset(sx.exp(-sx.exp(a*x))-y, x, sx.Reals)
    first = result.subs(y, 0).subs(a, coefficient)
    second = result.subs(a, coefficient).subs(y, 0)
    simultaneous = result.subs({a: coefficient, y: 0}, simultaneous=True)
    assert first == second == simultaneous
    assert first != sx.EmptySet
    assert first.contains(member) == sx.true
    assert first.contains(nonmember) == sx.false


def test_opposing_nonfinite_subexpressions_are_not_cancelled_as_finite_arithmetic():
    x, a, b, y = sx.symbols('x a b y')
    result = sx.solveset(sx.exp(a-sx.exp(b*x))-y, x, sx.Reals)
    assignments = [(a, sx.oo), (b, sx.oo), (y, 0)]
    expected = result.subs(dict(assignments), simultaneous=True)
    for pairs in itertools.permutations(assignments):
        actual = result
        for parameter, value in pairs:
            actual = actual.subs(parameter, value)
        assert actual == expected
    assert expected != sx.EmptySet
    assert isinstance(expected, sx.ConditionSet)


def test_a_domain_containing_infinity_keeps_the_nonfinite_bound_equation():
    x, a, y = sx.symbols('x a y')
    result = sx.solveset(sx.exp(sx.exp(x)+a)-y, x, sx.UniversalSet)
    first = result.subs(y, 0).subs(a, -sx.oo)
    second = result.subs(a, -sx.oo).subs(y, 0)
    simultaneous = result.subs({y: 0, a: -sx.oo}, simultaneous=True)
    assert first == second == simultaneous
    assert first != sx.UniversalSet
    assert first != sx.EmptySet
    assert first.contains(0) == sx.true
    assert first.contains(sx.oo) != sx.true


def test_a_mixed_finite_domain_retains_the_undefined_infinite_point():
    x, a, y = sx.symbols('x a y')
    domain = sx.FiniteSet(0, sx.oo)
    result = sx.solveset(sx.exp(sx.exp(x)+a)-y, x, domain)
    specialized = result.subs(y, 0).subs(a, -sx.oo)
    assert specialized == result.subs(a, -sx.oo).subs(y, 0)
    assert specialized == result.subs({a: -sx.oo, y: 0}, simultaneous=True)
    assert specialized.contains(0) == sx.true
    assert specialized != domain
    assert specialized.contains(sx.oo) != sx.true


def test_negative_infinity_remains_an_accepted_bound_solution_in_an_extended_domain():
    x, y = sx.symbols('x y')
    result = sx.solveset(sx.exp(sx.exp(x))-y, x, sx.FiniteSet(-sx.oo, 0))
    assert result.subs(y, 1) == sx.FiniteSet(-sx.oo)
    assert result.subs(y, sx.E) == sx.FiniteSet(0)


@pytest.mark.parametrize('assumptions', [{'extended_real': True}, {'finite': False}])
@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
def test_parameters_explicitly_accepting_infinity_keep_their_partial_substitution(assumptions, domain):
    x, y = sx.symbols('x y')
    a = sx.Symbol('a', **assumptions)
    result = sx.solveset(sx.exp(sx.exp(x)+a)-y, x, domain)
    assert result.subs(y, 0).subs(a, -sx.oo) == domain
    assert result.subs(a, -sx.oo).subs(y, 0) == domain
    assert result.subs({a: -sx.oo, y: 0}, simultaneous=True) == domain
