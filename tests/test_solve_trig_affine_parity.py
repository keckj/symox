"""Complete finite exponential representatives for shifted sec/csc arguments.

For reduced nonzero p/q the certificate is the degree 2*abs(p)
polynomial in t=exp(I*x/q), with all nonzero roots retained. These finite
principal-log representatives are solve's convention; solveset owns the
infinite periodic families. Nonintegral slopes extend SymPy 1.14's
NotImplementedError rather than claiming differential parity there.
"""
import pytest
import sympy as sp
import symox as sx


def _equation(module, name, slope, shift, rhs, assumptions=None):
    x = module.Symbol('x', **(assumptions or {}))
    parse = lambda value: module.sympify(value)
    return x, getattr(module, name)(parse(slope)*x+parse(shift))-parse(rhs)


def _roots(name, slope, shift, rhs, **options):
    x, equation = _equation(sx, name, slope, shift, rhs)
    return [sx.to_sympy(root) for root in sx.solve(equation, x, **options)]


def _bijection(actual, expected):
    remaining = [value.evalf(65) for value in expected]
    tolerance = sp.Float('1e-48', 65)
    for value in actual:
        numeric = value.evalf(65)
        matches = [i for i, reference in enumerate(remaining)
                   if abs(numeric-reference) < tolerance]
        assert len(matches) == 1, (value, remaining)
        remaining.pop(matches[0])
    assert not remaining


def _polynomial(name, slope, shift, rhs):
    slope, shift, rhs = map(sp.sympify, (slope, shift, rhs))
    if slope < 0:
        slope, shift = -slope, -shift
        if name == 'csc':
            rhs = -rhs
    p, q = map(int, slope.as_numer_denom())
    t = sp.Symbol('t')
    phase = sp.exp(sp.I*shift)
    polynomial = (rhs*phase**2*t**(2*p)-2*phase*t**p+rhs if name == 'sec'
                  else rhs*phase**2*t**(2*p)-2*sp.I*phase*t**p-rhs)
    return t, polynomial, p, q


def _complete_images(roots, name, slope, shift, rhs):
    t, polynomial, p, q = _polynomial(name, slope, shift, rhs)
    # Square-free over the exact coefficient field before numeric isolation:
    # rhs=+-1 gives p distinct images, not 2*p repeated ones.
    exact = sp.Poly(polynomial, t).sqf_part()
    expected = exact.nroots(n=65, maxsteps=300)
    images = [sp.exp(sp.I*root/q) for root in roots]
    assert len(images) == exact.degree()
    _bijection(images, expected)
    for root, image in zip(roots, images):
        assert not root.has(sp.Float)
        assert -sp.pi*q < sp.re(root.evalf(65)) <= sp.pi*q
        assert abs(polynomial.subs(t, image).evalf(65)) < sp.Float('1e-47', 65)


@pytest.mark.parametrize('name,slope,shift,rhs', [
    ('sec', '2', '1', '-1'), ('sec', '-2', '1', '-1'),
    ('csc', '2', '1', '2'), ('csc', '-2', '1', '2'),
    ('sec', '1', '10', '-1'), ('csc', '1', '10', '2'),
    ('sec', '2', 'I', '2'), ('csc', '2', '1+I', '-2'),
    ('sec', '2', 'pi/3', '2'), ('csc', '2', '-pi/4', '-1'),
])
def test_all_shifted_representatives_match_sympy(name, slope, shift, rhs):
    actual = _roots(name, slope, shift, rhs)
    x, equation = _equation(sp, name, slope, shift, rhs)
    expected = sp.solve(equation, x)
    _bijection(actual, expected)
    _complete_images(actual, name, slope, shift, rhs)


@pytest.mark.parametrize('name,slope,shift,rhs', [
    ('sec', '3', '1', '3'), ('csc', '-3', 'I', '2'),
    ('sec', '1/2', '1', '2'), ('csc', '3/2', '1', '2'),
    ('sec', '-2/3', '1+I', '3'), ('csc', '-3/2', '-10', '-1'),
    ('sec', '4/3', 'pi/3', '-1'), ('csc', '2/3', 'I', '1'),
])
def test_all_affine_polynomial_preimages(name, slope, shift, rhs):
    _complete_images(_roots(name, slope, shift, rhs), name, slope, shift, rhs)


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_exact_polynomial_residuals(name):
    roots = _roots(name, '2', '1', '2')
    t, polynomial, _, q = _polynomial(name, '2', '1', '2')
    for root in roots:
        image = sp.exp(sp.I*root/q)
        assert sp.simplify(sp.expand(polynomial.subs(t, image))) == 0


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_symbolic_shift_preserves_every_specialized_branch(name):
    x, b = sx.symbols('x b')
    actual = [sx.to_sympy(root) for root in sx.solve(getattr(sx, name)(2*x+b)-2, x)]
    assert len(actual) == 4
    for shift in (0, 1, 10, -10, sp.I, sp.pi/3):
        specialized = [root.subs(sp.Symbol('b'), shift) for root in actual]
        _complete_images(specialized, name, '2', shift, '2')


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_symbolic_target_specializes_at_double_roots(name):
    x, a = sx.symbols('x a')
    roots = [sx.to_sympy(root) for root in sx.solve(getattr(sx, name)(2*x+1)-a, x)]
    assert len(roots) == 4
    for target in (-1, 1, 2, -2, 2*sp.I):
        # A double quadratic root merges equal representatives; it must
        # retain both distinct square-root preimages.
        specialized = list(set(sp.simplify(root.subs(sp.Symbol('a'), target)) for root in roots))
        _complete_images(specialized, name, '2', '1', target)
    assert sx.solve(getattr(sx, name)(2*x+1), x) == []
    assert all(root.subs(sp.Symbol('a'), 0).has(sp.nan, sp.zoo, sp.oo) for root in roots)


@pytest.mark.parametrize('assumptions', [{'real': True}, {'positive': True}, {'negative': True}])
def test_shifted_assumption_decisions_match_reference(assumptions):
    x, equation = _equation(sx, 'sec', '2', '1', '-1', assumptions)
    reference_x, reference_equation = _equation(sp, 'sec', '2', '1', '-1', assumptions)
    actual = [sx.to_sympy(root) for root in sx.solve(equation, x)]
    _bijection(actual, sp.solve(reference_equation, reference_x))


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_complex_shift_and_real_domain(name):
    x = sx.Symbol('x', real=True)
    assert sx.solve(getattr(sx, name)(2*x+sx.I)-2, x) == []
    assert len(sx.solve(getattr(sx, name)(2*x+sx.I)-2, x, check=False)) == 4


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_zero_slopes_constants_and_symbolic_slope_specialization(name):
    x, a = sx.symbols('x a')
    original = getattr(sx, name)(a*x+sx.pi/2)-1
    # Preserve solve's existing constant-equation convention (no unknown
    # assignment), including the csc identity at a=0.
    assert sx.solve(original.subs(a, 0), x) == []
    assert sx.solve(getattr(sx, name)(sx.pi/4)-2, x) == []
    assert sx.solve(getattr(sx, name)(0*x+sx.pi/4)-2, x) == []


@pytest.mark.parametrize('options', [
    {'simplify': False}, {'check': False}, {'cubics': False}, {'quartics': False},
])
def test_affine_routes_preserve_solver_options(options):
    actual = _roots('csc', '3/2', '1', '2', **options)
    _complete_images(actual, 'csc', '3/2', '1', '2')


def test_dictionary_and_set_output_protocols():
    x = sx.Symbol('x')
    equation = sx.csc(2*x+1)-2
    plain = sx.solve(equation, x)
    assert sx.solve(equation, x, dict=True) == [{x: root} for root in plain]
    symbols, solutions = sx.solve(equation, x, set=True)
    assert symbols == [x]
    assert set(solutions) == {(root,) for root in plain}


@pytest.mark.parametrize('name', ['sec', 'csc'])
@pytest.mark.parametrize('rational', [True, False])
def test_shifted_high_precision_target_keeps_requested_precision(name, rational):
    x = sx.Symbol('x')
    target_text = '2.12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789'
    target = sx.Float(target_text, 90)
    roots = [sx.to_sympy(root) for root in sx.solve(getattr(sx, name)(2*x+1)-target, x, rational=rational)]
    assert len(roots) == 4
    t, polynomial, _, q = _polynomial(name, '2', '1', sp.Float(target_text, 90))
    expected = sp.Poly(polynomial, t).nroots(n=75, maxsteps=300)
    _bijection([sp.exp(sp.I*root/q) for root in roots], expected)
    for root in roots:
        atoms = root.atoms(sp.Float)
        if rational:
            assert not atoms
        else:
            assert atoms and all(value._prec >= 280 for value in atoms)


@pytest.mark.parametrize('name', ['sec', 'csc'])
def test_symbolic_target_proven_zero_has_no_finite_preimage(name):
    x = sx.Symbol('x')
    zero = sx.Symbol('zero', zero=True)
    assert sx.solve(getattr(sx, name)(2*x+1)-zero, x) == []
