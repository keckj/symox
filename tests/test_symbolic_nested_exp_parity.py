"""Complete guarded logarithmic inverses of nested affine exponentials.

For a nonzero finite y, exp(exp(x))=y is equivalent to
exp(x)=log(y)+2*pi*I*n != 0, followed by every complex log branch k.
The index domain, not just residuals of returned samples, is part of the answer.
"""
from itertools import product

import pytest
import sympy as sp
import symox as sx


x, y = sx.symbols('x y')


def complex_image(result):
    """Return the integer variables, composed map and guarded index domain."""
    while isinstance(result, sx.ConditionSet):
        result = result.base_set
    assert isinstance(result, sx.ImageSet), result
    variables = result.lamda.variables
    if len(variables) == 1 and isinstance(result.base_set, sx.ImageSet):
        inner_variables, inner_body, base = complex_image(result.base_set)
        return inner_variables, result.lamda.expr.xreplace({variables[0]: inner_body}), base
    return variables, result.lamda.expr, result.base_set


def finite_values(result):
    assert isinstance(result, sx.FiniteSet), result
    return list(result)


def assert_whole_domain(result, domain):
    """Prove full coverage using only parameter conditions on subsets of D."""
    def coverage(subset):
        if subset == domain:
            return sp.true
        if subset == sx.EmptySet:
            return sp.false
        if isinstance(subset, sx.Union):
            return sp.Or(*(coverage(part) for part in subset.args))
        assert isinstance(subset, sx.ConditionSet), subset
        # Its base is recursively a subset of D. With a condition independent
        # of the point being selected, it covers that base exactly when the
        # parameter condition holds.
        assert not subset.condition.has(subset.sym)
        return sp.And(sx.to_sympy(subset.condition), coverage(subset.base_set))

    # In particular, G*D union ~G*D is D for every parameter value, without
    # relying on the native union constructor to absorb conditional subsets.
    assert sp.simplify_logic(coverage(result)) is sp.true


def test_symbolic_family_has_all_two_integer_branches_and_the_exact_nonzero_guard():
    result = sx.solveset(sx.exp(sx.exp(x))-y, x)
    assert result.free_symbols == {y}
    assert result.subs(y, 0) == sx.EmptySet
    variables, body, base = complex_image(result)
    assert len(variables) == 2 and all(v.is_Dummy for v in variables)
    k, n = map(sx.to_sympy, variables)
    Y = sp.Symbol('y')
    first_log = sp.log(sp.Abs(Y))+sp.I*(sp.arg(Y)+2*sp.pi*n)
    expected = sp.log(sp.Abs(first_log))+sp.I*(sp.arg(first_log)+2*sp.pi*k)
    assert sp.simplify(sx.to_sympy(body)-expected) == 0
    assert isinstance(base, sx.ConditionSet)
    assert set(base.sym) == set(variables)
    assert sx.to_sympy(base.condition) == sp.Ne(first_log, 0)
    # The bound indices must not leak into the set's free parameters.
    assert base.free_symbols == {y}


@pytest.mark.parametrize('target', [sx.Integer(1), sx.Integer(2), sx.Integer(-1), sx.Integer(-2), sx.I, 2+sx.I])
def test_substitution_matches_direct_family_and_excludes_exactly_zero_preimages(target):
    generic = sx.solveset(sx.exp(sx.exp(x))-y, x).subs(y, target)
    direct = sx.solveset(sx.exp(sx.exp(x))-target, x)
    variables, body, base = complex_image(generic)
    direct_vars, direct_body, direct_base = complex_image(direct)
    # Alpha conversion compares the complete symbolic maps and index guards.
    rename = dict(zip(direct_vars, variables))
    assert direct_body.xreplace(rename) == body
    renamed_base = direct_base.xreplace(rename)
    if renamed_base != base:
        # Direct inversion can certify that n+1/2 or n+1/4 never vanishes
        # for any integer n, while substitution retains the exact guard.
        assert isinstance(base, sx.ConditionSet)
        assert renamed_base == base.base_set
        assert sp.simplify(sp.expand(sx.to_sympy(base.condition))) is sp.true
    assert not sx.to_sympy(body).has(sp.Float)
    k, n = variables
    target_sp = sx.to_sympy(target)
    for kv, nv in product(range(-2, 3), repeat=2):
        first_log = sp.log(sp.Abs(target_sp))+sp.I*(sp.arg(target_sp)+2*sp.pi*nv)
        admissible = first_log != 0
        assert base.contains(sx.Tuple(kv, nv)) == (sx.true if admissible else sx.false)
        if not admissible:
            continue
        actual = sx.to_sympy(body.xreplace({k: sx.Integer(kv), n: sx.Integer(nv)}))
        expected = sp.log(sp.Abs(first_log))+sp.I*(sp.arg(first_log)+2*sp.pi*kv)
        assert abs((actual-expected).evalf(65)) < sp.Float('1e-50', 65)
        assert abs((sp.exp(sp.exp(actual))-target_sp).evalf(65)) < sp.Float('1e-48', 65)


@pytest.mark.parametrize('target', [0, sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonzero_finite_target_is_required(target):
    assert sx.solveset(sx.exp(sx.exp(x))-y, x).subs(y, target) == sx.EmptySet


@pytest.mark.parametrize('target', [0, 1, sx.Rational(1, 2), -1, sx.I])
def test_real_nested_exponential_range_excludes_one_and_nonreal_values(target):
    generic = sx.solveset(sx.exp(sx.exp(x))-y, x, sx.Reals)
    assert generic.subs(y, target) == sx.EmptySet


@pytest.mark.parametrize('target', [2, 3, sx.E, sx.exp(sx.exp(2))])
def test_real_family_is_the_unique_double_log_on_its_full_range(target):
    result = sx.solveset(sx.exp(sx.exp(x))-y, x, sx.Reals).subs(y, target)
    actual, = finite_values(result)
    expected = sp.log(sp.log(sx.to_sympy(target)))
    assert sp.simplify(sx.to_sympy(actual)-expected) == 0


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
def test_zero_inner_affine_coefficient_retains_the_whole_domain_only_for_the_constant_value(domain):
    a = sx.Symbol('a', real=True)
    result = sx.solveset(sx.exp(sx.exp(a*x))-y, x, domain)
    collapsed = result.subs(a, 0)
    assert collapsed.subs(y, sx.E) == domain
    assert collapsed.subs(y, 2) == sx.EmptySet
    assert result.subs({a: 1, y: 0}) == sx.EmptySet


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
def test_each_affine_level_keeps_its_zero_multiplier_branch(domain):
    a, b, c, d = sx.symbols('a b c d', real=True)
    equation = sx.exp(c*sx.exp(a*x+b)+d)-y
    result = sx.solveset(equation, x, domain)
    assert_whole_domain(result.subs({c: 0, d: 0, y: 1}), domain)
    assert result.subs({c: 0, d: 0, y: 2}) == sx.EmptySet
    assert result.subs({a: 0, b: 0, c: 1, d: 0, y: sx.E}) == domain
    assert result.subs({a: 0, b: 0, c: 1, d: 0, y: 2}) == sx.EmptySet


def test_unknown_constant_equation_is_conditional_instead_of_a_false_empty_set():
    a = sx.Symbol('a')
    result = sx.solveset(a, x)
    assert isinstance(result, sx.ConditionSet)
    assert result.subs(a, 0) == sx.Complexes
    assert result.subs(a, 1) == sx.EmptySet


@pytest.mark.parametrize('source', ['exp(exp(x/a))-y', '(exp(exp(x))-y)/a', 'exp(exp(x)+1/a)-y'])
def test_parameter_poles_are_not_mistaken_for_zero_multiplier_solutions(source):
    a = sx.Symbol('a')
    equation = eval(source, dict(vars(sx), x=x, y=y, a=a))
    result = sx.solveset(equation, x)
    assert result.subs(a, 0).subs(y, 2) == sx.EmptySet
    assert result.subs(a, 0).subs(y, 0) == sx.EmptySet


@pytest.mark.parametrize('source', ['exp(2*exp(3*exp(x)))-y', 'exp(exp(exp(exp(x))))-y'])
def test_affine_maps_compose_before_further_logs_without_losing_integer_dimensions(source):
    equation = eval(source, dict(vars(sx), x=x, y=y))
    result = sx.solveset(equation, x)
    variables, body, base = complex_image(result)
    depth = 3 if source.startswith('exp(2') else 4
    assert len(variables) == depth
    assert result.free_symbols == {y}
    assert isinstance(base, sx.ConditionSet)
    assert result.subs(y, 0) == sx.EmptySet
    assert not body.has(equation+y)


def test_free_parameter_names_cannot_be_captured_by_the_affine_image_placeholder():
    z = sx.Symbol('z')
    parameter = sx.Symbol('x')
    result = sx.solveset(sx.exp(sx.exp(2*z))-parameter, z)
    assert result.free_symbols == {parameter}
    assert result.subs(parameter, 0) == sx.EmptySet


def test_tuple_condition_membership_substitutes_components_simultaneously():
    k, n = sx.symbols('k n', integer=True)
    domain = sx.ConditionSet(sx.Tuple(k, n), sx.Ne(n, 0) & sx.Eq(k+n, 0),
                             sx.ProductSet(sx.Integers, sx.Integers))
    assert domain.free_symbols == set()
    assert domain.contains(sx.Tuple(0, 0)) == sx.false
    assert domain.contains(sx.Tuple(-1, 1)) == sx.true
    assert domain.contains(sx.Tuple(1, 1)) == sx.false
    assert domain.contains(sx.Tuple(1,)) == sx.false


@pytest.mark.parametrize('domain', [sx.Complexes, sx.Reals])
@pytest.mark.parametrize('coefficient', [sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonfinite_multiplier_keeps_the_original_equation_instead_of_dividing_by_infinity(domain, coefficient):
    a = sx.Symbol('a')
    equation = sx.exp(sx.exp(a*x))-y
    result = sx.solveset(equation, x, domain)
    for target in (1, 2):
        substitution = {a: coefficient, y: target}
        specialized_equation = equation.subs(substitution)
        if specialized_equation == sx.nan:
            # An undefined specialization is a retained calculation, not a
            # proof of absence obtained from ordinary Eq(nan,0)->False.
            expected = sx.ConditionSet(x, sx.Eq(sx.nan, 0, evaluate=False), domain)
            assert result.subs(substitution) == expected
        else:
            original_equation = sx.Eq(specialized_equation, 0)
            # Keep the true nonfinite guard as provenance while the extended
            # calculation still depends on x. This exact tautology changes
            # neither the original equation nor any of its solutions.
            assert sx.Contains(coefficient, sx.Complexes) == sx.false
            retained_guard = ~sx._core._apply("Contains", coefficient, sx.Complexes, evaluate=False)
            expected = sx.ConditionSet(x, retained_guard & original_equation, domain)
            assert result.subs(substitution) == expected
            assert sx.solveset(specialized_equation, x, domain) == sx.ConditionSet(x, original_equation, domain)
    if coefficient in (sx.oo, -sx.oo):
        # For a=+oo, x<0 has exp(a*x)=0; a=-oo gives x>0.
        # Its unsupported equation must survive even when finite-a inversion
        # proves no real roots. Returning EmptySet would lose that half-line.
        assert result.subs({a: coefficient, y: 1}) != sx.EmptySet


def test_nonfinite_inner_shift_survives_an_empty_finite_preimage():
    a = sx.Symbol('a')
    equation = sx.exp(sx.exp(x)+a)
    result = sx.solveset(equation, x)
    specialized = equation.subs(a, -sx.oo)
    assert result.subs(a, -sx.oo) == sx.ConditionSet(x, sx.Eq(specialized, 0), sx.Complexes)


def test_clearing_a_parameter_denominator_preserves_its_nonfinite_specialization():
    a = sx.Symbol('a')
    equation = (sx.exp(sx.exp(x))-y)/a
    result = sx.solveset(equation, x)
    assert result.subs({a: sx.oo, y: 2}) == sx.Complexes


def test_partial_target_substitution_preserves_a_later_nonfinite_shift():
    # Also fails on the pre-wave native 8323eddaa68e4c1f2c55cff6933c84fa:
    # its plain original ConditionSet loses the Eq on substituting y=0.
    # Substituting a=-oo first returns Complexes there; at x=1 the original
    # expression specialized in either parameter order is exactly zero.
    a = sx.Symbol('a')
    equation = sx.exp(sx.exp(x)+a)-y
    result = sx.solveset(equation, x)
    assert result.subs(y, 0).subs(a, -sx.oo) == sx.Complexes
    assert result.subs(a, -sx.oo).subs(y, 0) == sx.Complexes
    assert result.subs({y: 0, a: -sx.oo}, simultaneous=True) == sx.Complexes


def test_cartesian_modulus_presentation_is_identical_on_the_whole_real_plane():
    # Integer branch indices are real. Thus the newer sqrt presentation in
    # the legacy display contracts equals the former Abs at every index.
    real, imaginary = sp.symbols('real imaginary', real=True)
    assert sp.simplify(sp.Abs(real+sp.I*imaginary)-sp.sqrt(real**2+imaginary**2)) == 0


def test_triple_exponential_one_excludes_each_zero_preimage_on_the_entire_integer_domain():
    result = sx.solveset(sx.exp(sx.exp(sx.exp(x)))-1, x)
    variables, body, base = complex_image(result)
    k2, k1, n = map(sx.to_sympy, variables)
    first = 2*sp.pi*sp.I*n
    second = sp.log(2*sp.pi*sp.Abs(n))+sp.I*(sp.arg(sp.I*n)+2*sp.pi*k1)
    expected = sp.log(sp.Abs(second))+sp.I*(sp.arg(second)+2*sp.pi*k2)
    assert sp.simplify(sx.to_sympy(body)-expected) == 0
    assert sx.to_sympy(base.condition) == sp.And(sp.Ne(first, 0), sp.Ne(second, 0))
