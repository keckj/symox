"""Differential classification, characteristics and separation regressions."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def ode(module, source):
    x = module.Symbol('x')
    f = module.Function('f')(x)
    return eval(source, vars(module) | {'x':x, 'f':f}), f


@pytest.mark.parametrize('source,required', [
    ('diff(f,x)', '1st_linear'),
    ('diff(f,x)-f', 'separable'),
    ('diff(f,x)+x*f-x', '1st_linear'),
    ('diff(f,x)-x*f*f', 'Bernoulli'),
    ('diff(f,x)+f-x*f*f', 'Bernoulli'),
    ('diff(f,x,2)+f', 'nth_linear_constant_coeff_homogeneous'),
    ('diff(f,x,3)-f', 'nth_linear_constant_coeff_homogeneous'),
    ('diff(f,x,2)+f-x', 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('diff(f,x,2)+f-sin(x)', 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('diff(f,x,2)+f-exp(2*x)', 'nth_linear_constant_coeff_undetermined_coefficients'),
    ('diff(f,x,2)+f-log(x)', 'nth_linear_constant_coeff_variation_of_parameters'),
    ('x*x*diff(f,x,2)+x*diff(f,x)-f', 'nth_linear_euler_eq_homogeneous'),
    ('(x*x+1)*diff(f,x)+2*x*f', '1st_exact'),
    ('(x+cos(f))*diff(f,x)+f+2*x', '1st_exact'),
    ('Eq(diff(f,x),f)', 'separable'),
])
def test_native_ode_family_labels_are_reference_classes(source, required):
    native, function = ode(sx, source)
    foreign, reference_function = ode(sp, source)
    hints = sx.classify_ode(native, function)
    reference_hints = sp.classify_ode(foreign, reference_function)
    assert required in hints
    assert hints == reference_hints
    assert sx.classify_ode(native) == hints
    data = sx.classify_ode(native, function, dict=True)
    reference_data = sp.classify_ode(foreign, reference_function, dict=True)
    assert data['order'] == sp.ode_order(foreign, reference_function)
    assert set(data) == set(reference_data)
    # The reference's default is not always the first hint (f' = 0 lists
    # nth_algebraic first and defaults to the constant coefficient method).
    assert data['default'] == reference_data['default']


@pytest.mark.parametrize('source', [
    'diff(f,x)-f', 'diff(f,x)-x*f*f', 'diff(f,x)+f-x*f*f',
    'diff(f,x,2)+f', 'diff(f,x,2)+f-x',
    'x*x*diff(f,x,2)+x*diff(f,x)-f',
    '(x*x+1)*diff(f,x)+2*x*f',
])
def test_accepted_ode_hints_produce_valid_solutions(source):
    native, function = ode(sx, source)
    reference, reference_function = ode(sp, source)
    for hint in sx.classify_ode(native, function):
        # Truncated series, held integrals and the conserved-quantity forms
        # (whose check raises in the reference) cannot be checked exactly.
        # 2nd_linear_bessel matches Euler equations, where the reference
        # returns the incomplete Eq(f(x), C1).
        if (hint.endswith('_Integral') or 'series' in hint
                or 'conserved' in hint or hint == '2nd_linear_bessel'):
            continue
        solution = sx.dsolve(native, function, hint=hint)
        solutions = solution if isinstance(solution, list) else [solution]
        assert solutions
        for item in solutions:
            assert sx.checkodesol(native, item, function) == (True, 0)
            assert sp.checkodesol(reference, item.to_sympy(), reference_function)[0]


def test_exact_implicit_equations_are_solved_and_checked():
    x = sx.Symbol('x')
    f = sx.Function('f')(x)
    equation = (x+sx.cos(f))*f.diff(x)+f+2*x
    # The reference default is 'factorable', whose Lie group search raises
    # NotImplementedError after a long computation; name the exact method.
    solution = sx.dsolve(equation, f, hint='1st_exact')
    assert solution.lhs.has(sx.sin(f))
    # The implicit solution cannot be solved for f(x); like the reference,
    # checkodesol raises NotImplementedError unless solve_for_func=False.
    with pytest.raises(NotImplementedError):
        sx.checkodesol(equation, solution, f)
    assert sx.checkodesol(equation, solution, f, solve_for_func=False) == (True, 0)
    wrong = sx.Eq(x*x+x*f+sx.sin(f)+x, sx.Symbol('C1'))
    checked, residual = sx.checkodesol(equation, wrong, f, solve_for_func=False)
    assert not checked
    assert residual != 0


def test_classification_handles_zero_order_and_unimplemented_families():
    x = sx.Symbol('x')
    f = sx.Function('f')(x)
    # As in the reference, an equation without derivatives is algebraic
    # and f'**2 = f is reached through the factorable and Lie group methods.
    assert sx.classify_ode(f-1, f) == ('nth_algebraic', 'nth_algebraic_Integral')
    assert sx.classify_ode(f-1, f, dict=True)['default'] == 'nth_algebraic'
    assert sx.classify_ode(f.diff(x)**2-f, f) == ('factorable', 'lie_group')
    with pytest.raises(ValueError):
        sx.dsolve(f.diff(x,2)+f, f, hint='separable')


def test_exact_constants_cannot_capture_existing_user_symbols():
    x, C1 = sx.symbols('x C1')
    f = sx.Function('f')(x)
    equation = (x+sx.cos(f))*f.diff(x)+f+2*x+C1
    solution = sx.dsolve(equation, f, hint='1st_exact')
    assert C1 in solution.free_symbols
    assert sx.Symbol('C2') in solution.free_symbols
    assert sx.checkodesol(equation, solution, f, solve_for_func=False) == (True, 0)


def pde(module, source):
    x, y = module.symbols('x y')
    u = module.Function('u')(x,y)
    namespace = vars(module) | {'x':x, 'y':y, 'u':u}
    return eval(source, namespace), u


@pytest.mark.parametrize('source', [
    'diff(u,x)+diff(u,y)', '2*diff(u,x)+3*diff(u,y)+4*u',
    'diff(u,x)+diff(u,y)-x', '2*diff(u,x)-diff(u,y)-x*y',
    'diff(u,x)+2*diff(u,y)+u-exp(x+y)',
    '-2*diff(u,x)+4*diff(u,y)+5*u-exp(x+3*y)',
    'diff(u,x)', 'diff(u,y)', 'diff(u,x)+2*u-y',
    'diff(u,y)-u-x',
    'x*diff(u,x)+y*diff(u,y)',
    '2*x*diff(u,x)+3*y*diff(u,y)+4*u',
    'diff(u,x)+x*diff(u,y)',
    'y*diff(u,x)+diff(u,y)',
    'x*diff(u,x)+diff(u,y)',
    'diff(u,x)+y*diff(u,y)',
    'x*diff(u,x)+u-y',
    'y*diff(u,y)+u-x',
])
def test_linear_pde_classification_and_solutions(source):
    native, function = pde(sx, source)
    reference, reference_function = pde(sp, source)
    assert sx.classify_pde(native, function) == sp.classify_pde(reference, reference_function)
    solution = sx.pdsolve(native, function)
    assert sx.checkpdesol(native, solution, function) == (True, 0)
    assert sp.checkpdesol(reference, solution.to_sympy(), reference_function)[0]
    assert sx.pdsolve(native) == solution
    info = sx.classify_pde(native, function, dict=True)
    assert info['order'] == 1
    assert info['default'] == info['ordered_hints'][0]


@pytest.mark.parametrize('source', [
    'diff(u,x)+diff(u,y)', '2*diff(u,x)+3*diff(u,y)+4*u',
    'diff(u,x)+diff(u,y)-x', 'diff(u,x)+2*diff(u,y)+u-exp(x+y)',
    '-2*diff(u,x)+4*diff(u,y)+5*u-exp(x+3*y)',
])
def test_constant_characteristic_form_matches_reference(source):
    native, function = pde(sx, source)
    reference, reference_function = pde(sp, source)
    actual = sx.pdsolve(native, function).rhs.to_sympy()
    expected = sp.pdsolve(reference, reference_function).rhs
    assert sp.simplify(actual-expected) == 0


def test_pde_arbitrary_functions_and_characteristic_coordinates_do_not_capture():
    x, y, xi, eta = sx.symbols('x y xi eta')
    u, F = sx.Function('u')(x,y), sx.Function('F')
    equation = u.diff(x)+u.diff(y)+xi*u-eta
    solution = sx.pdsolve(equation, u, solvefun=sx.Function('H'))
    assert solution.has(sx.Function('H'))
    # The antiderivative of eta*exp(xi*t) keeps the xi = 0 branch, as the
    # reference does, so each branch of the solution is checked separately.
    assert solution.has(sx.Piecewise)
    for value in (sx.Integer(0), sx.Integer(3), sx.Symbol('k', nonzero=True)):
        assert sx.checkpdesol(equation.subs(xi, value), solution.subs(xi, value), u) == (True, 0)
    forcing = u.diff(x)+u.diff(y)-F(x-y)
    solution = sx.pdsolve(forcing, u)
    assert solution.has(sx.Function('F_'))
    assert sx.checkpdesol(forcing, solution, u) == (True, 0)


def test_integral_pde_hint_remains_held_and_satisfies_equation():
    x, y = sx.symbols('x y')
    u, forcing = sx.Function('u')(x,y), sx.Function('G')(x,y)
    equation = 2*u.diff(x)+3*u.diff(y)+4*u-forcing
    # Fresh dummy masks repeatedly hit different free-coordinate masks.
    # Every choice must give the same proof of the characteristic identity.
    for _ in range(80):
        solution = sx.pdsolve(equation, u, hint='1st_linear_constant_coeff_Integral')
        assert solution.has(sx.Integral)
        assert sx.checkpdesol(equation, solution, u) == (True, 0)


def test_pde_solution_checks_lists_reversed_equalities_and_incorrect_values():
    x, y = sx.symbols('x y')
    u = sx.Function('u')(x,y)
    equation = u.diff(x)+u.diff(y)
    assert sx.checkpdesol(equation, sx.Eq(x-y,u), u) == (True, 0)
    assert sx.checkpdesol(equation, [sx.Eq(u,x-y), sx.Eq(u,x+y)], u) == [(True,0),(False,2)]
    # The reference pdsolve accepts dict= but returns the equation.
    assert sx.pdsolve(equation, u, dict=True) == sx.Eq(u, sx.Function('F')(x-y))
    assert sx.pdsolve(equation, u, hint='all')['order'] == 1


@pytest.mark.parametrize('source,strategy,groups', [
    ('diff(u,x,2)-diff(u,y,2)', 'mul', ['X(x)','Y(y)']),
    ('diff(u,x)-diff(u,y,2)', 'mul', ['X(x)','Y(y)']),
    ('diff(u,x,2)+diff(u,y,2)', 'mul', ['X(x)','Y(y)']),
    ('diff(u,x,2)-x*diff(u,y,2)', 'mul', ['X(x)','Y(y)']),
    ('diff(u,x)-exp(u)*diff(u,y)', 'add', ['X(x)','Y(y)']),
    ('diff(u,x)+diff(u,y)', 'add', ['X(x)','Y(y)']),
])
def test_pde_separation_matches_reference(source, strategy, groups):
    results = []
    for module in (sx, sp):
        equation, function = pde(module, source)
        x, y = function.args
        namespace = {'x':x, 'y':y, 'X':module.Function('X'), 'Y':module.Function('Y')}
        separated = [eval(group, namespace) for group in groups]
        results.append(module.pde_separate(equation, function, separated, strategy=strategy))
    native, reference = results
    assert native is not None and reference is not None
    assert all(sp.simplify(left.to_sympy()-right) == 0 for left, right in zip(native, reference))


def test_pde_separation_handles_groups_and_rejects_inconsistent_variables():
    x, y, z = sx.symbols('x y z')
    u = sx.Function('u')(x,y,z)
    X, V = sx.Function('X')(x), sx.Function('V')(y,z)
    result = sx.pde_separate_mul(u.diff(x,2)-u.diff(y,2)-u.diff(z,2), u, [X,V])
    assert result is not None
    assert not result[0].has(y,z)
    assert not result[1].has(x)
    with pytest.raises(ValueError):
        sx.pde_separate(u.diff(x), u, [X,sx.Function('W')(x,z)])
    assert sx.classify_pde(sx.Function('f')(x,y).diff(x,2), sx.Function('f')(x,y)) == ()


def test_differential_interfaces_do_not_import_sympy():
    source = '''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y');f=Function('f')(x);u=Function('u')(x,y)
assert 'separable' in classify_ode(f.diff(x)-f,f)
equation=u.diff(x)+u.diff(y)-x
assert checkpdesol(equation,pdsolve(equation,u),u)==(True,0)
assert pde_separate_mul(u.diff(x,2)-u.diff(y,2),u,[Function('X')(x),Function('Y')(y)])
g=Function('g')(x)
solution=dsolve([f.diff(x)-f-g,g.diff(x)-g],[f,g],ics={f.subs(x,0):1,g.subs(x,0):2})
assert solution[0].rhs.subs(x,0)==1
assert solution[1].rhs.subs(x,0)==2
'''
    subprocess.run([sys.executable,'-c',source],check=True)


def test_held_integral_forms_follow_reference_doit():
    # Integral.doit factors an undone (x, b) limit, and Derivative.doit is done
    # again when its result still holds a derivative.
    for module in (sx, sp):
        x, y = module.symbols('x y')
        f, g = module.Function('f'), module.Function('g')
        held = module.Integral(3*(3*f(x)-2*g(x))*module.exp(x)/13, (x, 2*y)).doit()
        assert held.as_coeff_Mul()[0] == module.Rational(3, 13)
    x, y = sx.symbols('x y')
    u, forcing = sx.Function('u')(x,y), sx.Function('G')(x,y)
    equation = 2*u.diff(x)+3*u.diff(y)+4*u-forcing
    solution = sx.pdsolve(equation, u, hint='1st_linear_constant_coeff_Integral')
    assert sx.simplify(equation.subs(u, solution.rhs).doit()) == 0
