"""All finite exponential preimages for sec/csc integer multiples.

SymPy's solve rewrites these equations as a polynomial in exp(I*x). Each
nonzero polynomial root gives its principal logarithm; division of only the
inverse trig branches by the multiple loses some of those representatives.
"""
import pytest
import sympy as sp
import symox as sx


def _module_equation(module, name, multiple, rhs, assumptions):
    x = module.Symbol('x', **assumptions)
    return x, getattr(module, name)(multiple*x) - module.sympify(rhs)


def _assert_numeric_bijection(actual, expected, substitutions=None):
    """Match every exact output with a distinct reference at high precision."""
    substitutions = substitutions or {}
    pool = [value.subs(substitutions).evalf(65) for value in expected]
    tolerance = sp.Float('1e-48', 65)
    for value in actual:
        assert not value.has(sp.Float)
        numeric = value.subs(substitutions).evalf(65)
        matches = [i for i, reference in enumerate(pool)
                   if abs(numeric-reference) < tolerance]
        assert len(matches) == 1, (value, pool)
        pool.pop(matches[0])
    assert not pool


@pytest.mark.parametrize('name,multiple,rhs', [
    ('sec', 2, '-1'),
    ('sec', 2, '1'),
    ('sec', 3, '-1'),
    ('sec', 3, '1'),
    ('csc', 2, '2'),
    ('csc', 2, '-2'),
    ('csc', 3, '1'),
    ('csc', 3, '-1'),
    ('sec', 4, '3'),
    ('csc', 4, '2'),
    ('sec', -4, '3'),
    ('csc', -3, '2'),
])
def test_complete_reference_representatives(name, multiple, rhs):
    x, equation = _module_equation(sx, name, multiple, rhs, {})
    reference_x, reference_equation = _module_equation(sp, name, multiple, rhs, {})
    actual = [sx.to_sympy(root) for root in sx.solve(equation, x)]
    expected = sp.solve(reference_equation, reference_x)
    assert len(actual) == len(expected)
    _assert_numeric_bijection(actual, expected)
    for root in actual:
        assert abs(reference_equation.subs(reference_x, root).evalf(65)) < sp.Float('1e-48', 65)


def test_exact_inverse_double_root_keeps_both_signs():
    x = sx.Symbol('x')
    assert sx.solve(sx.sec(2*x) + 1, x) == [-sx.pi/2, sx.pi/2]


@pytest.mark.parametrize('assumptions', [{'real': True}, {'positive': True}, {'negative': True}])
def test_candidate_assumptions_are_applied_after_all_preimages(assumptions):
    x, equation = _module_equation(sx, 'sec', 2, '-1', assumptions)
    reference_x, reference_equation = _module_equation(sp, 'sec', 2, '-1', assumptions)
    actual = [sx.to_sympy(root) for root in sx.solve(equation, x)]
    expected = sp.solve(reference_equation, reference_x)
    _assert_numeric_bijection(actual, expected)


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_reciprocal_trig_has_no_finite_zero(name):
    x = sx.Symbol('x')
    assert sx.solve(getattr(sx, name)(4*x), x) == []


@pytest.mark.parametrize('name,rhs', [('sec', '3'), ('csc', '2')])
@pytest.mark.parametrize('multiple', [5, 6, 7])
def test_general_integer_multiple_has_every_algebraic_preimage(name, rhs, multiple):
    # An independent polynomial computation checks all exp(I*x) images,
    # without asking solve to simplify many equivalent logarithmic radicals.
    x, equation = _module_equation(sx, name, multiple, rhs, {})
    actual = [sx.to_sympy(root) for root in sx.solve(equation, x)]
    assert len(actual) == 2*multiple
    t = sp.Symbol('t')
    rhs = sp.sympify(rhs)
    polynomial = (rhs*t**(2*multiple) - 2*t**multiple + rhs if name == 'sec'
                  else rhs*t**(2*multiple) - 2*sp.I*t**multiple - rhs)
    expected = list(sp.Poly(polynomial, t).nroots(n=65, maxsteps=200))
    images = [sp.exp(sp.I*root).evalf(65) for root in actual]
    tolerance = sp.Float('1e-48', 65)
    for value in images:
        matches = [i for i, reference in enumerate(expected) if abs(value-reference) < tolerance]
        assert len(matches) == 1
        expected.pop(matches[0])
    assert not expected


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_nonreal_roots_are_rejected_for_a_real_unknown(name):
    x = sx.Symbol('x', real=True)
    assert sx.solve(getattr(sx, name)(4*x) - sx.Rational(1, 2), x) == []


def test_symbolic_coefficient_keeps_all_preimages():
    x, a = sx.symbols('x a')
    roots = [sx.to_sympy(root) for root in sx.solve(sx.sec(4*x) - a, x)]
    assert len(roots) == 8
    reference_a = sp.Symbol('a')
    t = sp.Symbol('t')
    for value in (sp.Integer(3), sp.Integer(-3), 2*sp.I):
        expected = list(sp.Poly(value*t**8 - 2*t**4 + value, t).nroots(n=65, maxsteps=200))
        for root in roots:
            assert not root.has(sp.Float)
            image = sp.exp(sp.I*root.subs(reference_a, value)).evalf(65)
            matches = [i for i, reference in enumerate(expected)
                       if abs(image-reference) < sp.Float('1e-48', 65)]
            assert len(matches) == 1
            expected.pop(matches[0])
        assert not expected


@pytest.mark.parametrize('assumptions', [{'positive': True}, {'negative': True}])
@pytest.mark.parametrize('sign', [-1, 1])
def test_real_logarithm_candidate_sign_follows_certified_algebraic_argument(assumptions, sign):
    p = sx.Symbol('p', **assumptions)
    candidate = sx.log(2 + sign*sx.sqrt(3))
    actual = [sx.to_sympy(root) for root in sx.solve(p-candidate, p)]
    reference_p = sp.Symbol('p', **assumptions)
    expected = sp.solve(reference_p-sp.log(2+sign*sp.sqrt(3)), reference_p)
    _assert_numeric_bijection(actual, expected)


@pytest.mark.parametrize('assumptions', [{'positive': True}, {'negative': True}, {'real': True}])
def test_symbolic_logarithm_sign_remains_undecided(assumptions):
    p = sx.Symbol('p', **assumptions)
    t = sx.Symbol('t', positive=True)
    assert sx.solve(p-sx.log(t), p) == [sx.log(t)]


def test_check_false_still_disables_candidate_assumption_filtering():
    x = sx.Symbol('x', real=True)
    assert len(sx.solve(sx.sec(4*x)-sx.Rational(1, 2), x, check=False)) == 8
