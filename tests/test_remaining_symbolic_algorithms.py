"""Independent identities for the inventory's symbolic algorithm adapters."""
import importlib
import math

import mpmath
import pytest
import symox as S

H = importlib.import_module('symox.simplify.hyperexpand')
O = importlib.import_module('symox.solvers.ode.ode')
C = importlib.import_module('symox.simplify.cse_main')
T = importlib.import_module('symox.simplify.simplify')
V = importlib.import_module('symox.solvers.solveset')
x = S.Symbol('x')
f = S.Function('f')(x)
g = S.Function('g')(x)
h = S.Function('h')(x)


def test_constants_avoid_existing_symbols_and_function_names():
    c1 = S.Symbol('C1')
    expr = c1+x+S.Function('C2')(x)
    assert O.get_numbered_constants(expr, 2) == S.symbols('C3 C4')
    assert O.constant_renumber(S.Symbol('C9')+S.Symbol('C7')*x) == S.Symbol('C2')+S.Symbol('C1')*x
    assert O.constant_renumber(S.Symbol('a')*x, variables=[x]) == S.Symbol('C1')*x


def _residual(eq, solution):
    from symox.solvers._differential import replace_function
    return S.expand(replace_function(eq, f, solution.rhs.removeO()))


@pytest.mark.parametrize('point', [0, 2])
def test_nonlinear_taylor_series_satisfies_ode_and_initial_value(point):
    equation = f.diff(x)-f**2
    result = S.dsolve(equation, f, hint='1st_power_series', x0=point, n=7,
                      ics={f.func(point): 3})
    expected = sum(3**(k+1)*(x-point)**k for k in range(7))
    assert S.expand(result.rhs.removeO()-expected) == 0
    shifted = _residual(equation, result).subs(x, x+point).expand()
    assert all(shifted.coeff(x, k) == 0 for k in range(6))


def test_second_order_series_with_forcing_and_ics():
    equation = f.diff(x, 2)+f-x
    result = S.dsolve(equation, f, hint='2nd_power_series_ordinary', n=8,
                     ics={f.func(0): 2, S.Subs(f.diff(x), x, 0): 3})
    assert result.rhs.removeO().subs(x, 0) == 2
    assert result.rhs.removeO().diff(x).subs(x, 0) == 3
    residual = _residual(equation, result)
    assert all(residual.coeff(x, k) == 0 for k in range(6))


def test_frobenius_two_nonresonant_branches():
    equation = x*x*f.diff(x, 2)+x*f.diff(x)+(x*x-S.Rational(1, 9))*f
    result = S.dsolve(equation, f, hint='2nd_power_series_regular', n=6)
    # Verify each truncated branch separately after removing its local O term.
    c1, c2 = S.symbols('C1 C2')
    for constant, other in [(c1, c2), (c2, c1)]:
        branch = result.rhs.subs({constant: 1, other: 0})
        branch = S.bottom_up(branch, lambda node: S.S.Zero if node.is_Order else node)
        from symox.solvers._differential import replace_function
        residual = S.expand(replace_function(equation, f, branch))
        normalized = S.expand(residual/x**(S.Rational(1, 3) if constant == c2 else -S.Rational(1, 3)))
        assert all(S.simplify(normalized.coeff(x, k)) == 0 for k in range(6))


def test_repeated_frobenius_roots_are_explicitly_rejected():
    with pytest.raises(NotImplementedError, match='Repeated'):
        S.dsolve(x*x*f.diff(x, 2)+x*f.diff(x)+x*x*f, f, hint='2nd_power_series_regular')


def test_system_inspection_and_triangular_nonlinear_solution():
    equations = [f.diff(x)-f, g.diff(x)-f*g]
    match = O.classify_sysode(equations, [f, g])
    assert not match['is_linear']
    assert match['order'] == {f: 1, g: 1}
    solutions = S.dsolve(equations,[f,g])
    assert '2nd_power_series_ordinary' in S.classify_ode(f.diff(x,2)+f,f,n=6)
    from symox.solvers._differential import replace_function
    for equation in equations:
        for solution in reversed(solutions):
            equation = replace_function(equation, solution.lhs, solution.rhs)
        assert S.simplify(equation) == 0
    assert len(set().union(*(sol.free_symbols for sol in solutions))-{x}) == 2


def test_linear_system_adapter_and_solution_scoring():
    match = O.classify_sysode([f.diff(x)-g, g.diff(x)+f], [f, g])
    assert match['is_linear']
    assert len(O.sysode_linear_2eq_order1(match)) == 2
    assert O.ode_sol_simplicity(S.Eq(f, x), f) == -2
    assert O.ode_sol_simplicity(S.Eq(f**2, x), f) == -1


def _op_numeric(operator, function, point, mp):
    coeffs = list(reversed(operator._poly.all_coeffs()))
    value = mp.mpf('0')
    for k, coefficient in enumerate(coeffs):
        coefficient = coefficient.subs(x, S.Rational(str(point)))
        value += mp.mpf(str(coefficient.evalf(35)))*mp.diff(lambda t: function(point*mp.exp(t)), 0, k)
    return value


@pytest.mark.parametrize('kind', ['upper_up', 'upper_down', 'lower_down', 'lower_up'])
def test_hyper_contiguous_operators_against_independent_mpmath(kind):
    mp = mpmath.mp.clone(); mp.dps = 35
    ap, bq = [S.Rational(7, 3)], [S.Rational(8, 3)]
    aa, bb = [mp.mpf(7)/3], [mp.mpf(8)/3]
    if kind == 'upper_up':
        operator = H.ShiftA(ap[0]); aa[0] += 1
    elif kind == 'upper_down':
        operator = H.UnShiftA(ap, bq, 0, x); aa[0] -= 1
    elif kind == 'lower_down':
        operator = H.ShiftB(bq[0]); bb[0] -= 1
    else:
        operator = H.UnShiftB(ap, bq, 0, x); bb[0] += 1
    point = mp.mpf('0.2')
    original = lambda z: mp.hyper([mp.mpf(7)/3], [mp.mpf(8)/3], z)
    assert abs(_op_numeric(operator, original, point, mp)-mp.hyper(aa, bb, point)) < mp.mpf('1e-28')


@pytest.mark.parametrize('group,direction', [(g, d) for g in range(4) for d in (-1, 1)])
def test_meijer_contiguous_operators_against_independent_mpmath(group, direction):
    mp = mpmath.mp.clone(); mp.dps = 35
    groups = [[S.Rational(1, 5)], [S.Rational(7, 5)], [S.Rational(3, 10)], [S.Rational(9, 5)]]
    numeric = [[mp.mpf(str(v.evalf(35))) for v in entries] for entries in groups]
    original = lambda z: mp.meijerg(numeric[:2], numeric[2:], z)
    classes = ([H.MeijerUnShiftB, H.MeijerUnShiftD, H.MeijerShiftA, H.MeijerShiftC] if direction > 0 else
               [H.MeijerShiftB, H.MeijerShiftD, H.MeijerUnShiftA, H.MeijerUnShiftC])
    klass = classes[group]
    operator = klass(*groups, 0, x) if 'UnShift' in klass.__name__ else klass(groups[group][0])
    target = [list(entries) for entries in numeric]; target[group][0] += direction
    point = mp.mpf('0.3')
    assert abs(_op_numeric(operator, original, point, mp)-mp.meijerg(target[:2], target[2:], point)) < mp.mpf('1e-25')


def test_shift_plan_composition_and_reduction():
    target, origin = H.Hyper_Function([4], [5]), H.Hyper_Function([2], [3])
    result = H.apply_operators(origin(x), H.devise_plan(target, origin, x), lambda y: x*y.diff(x))
    assert abs(float((result-target(x)).subs(x, S.Rational(1, 5)).evalf())) < 1e-12
    reduced, operators = H.reduce_order(H.Hyper_Function([3], [1]))
    assert reduced.sizes == (0, 0)
    result = H.apply_operators(S.exp(x), operators, lambda y: x*y.diff(x))
    assert S.simplify(result-S.exp(x)*(1+2*x+x*x/2)) == 0


def test_hyper_formula_basis_and_finite_polynomials():
    formula = H.Formula(H.Hyper_Function([2], []), x, (1-x)**-2, [])
    residual = x*formula.B.diff(x)-formula.M*formula.B
    assert all(S.simplify(value) == 0 for value in residual)
    polynomial = H.try_polynomial(H.Hyper_Function([-3, 2], [4]), x)
    assert polynomial == 1-S.Rational(3, 2)*x+S.Rational(9, 10)*x*x-x**3/5
    formula = H.FormulaCollection().lookup_origin(H.Hyper_Function([], []))
    assert formula.closed_form == S.exp(formula.z)
    meijer = H.MeijerFormulaCollection().lookup_origin(H.G_Function([], [], [2], []))
    assert meijer.closed_form == meijer.z**2*S.exp(-meijer.z)


def test_cse_dependency_sort_and_release_preserve_values():
    a, b, c = S.symbols('a b c')
    replacements = C.reps_toposort([(c, b+1), (b, a+1), (a, x*x)])
    assert [symbol for symbol, value in replacements] == [a, b, c]
    program, outputs = C.cse_release_variables(replacements, [c+b])
    environment = {x: S.Integer(3)}
    for symbol, value in program:
        if value is None:
            environment.pop(symbol, None)
        else:
            environment[symbol] = value.subs(environment)
    assert environment[outputs[0]] == 21
    with pytest.raises(ValueError): C.reps_toposort([(a, b), (b, a)])


def test_cse_argument_tracking_updates_both_indices():
    tracker = C.FuncArgTracker([x+f+g, x+f+h])
    common = tracker.func_to_argset[0]&tracker.func_to_argset[1]
    assert tracker.get_common_arg_candidates(common) == {0: 2, 1: 2}
    tracker.stop_arg_tracking(1)
    assert tracker.get_subset_candidates(common) == {0}


def test_exact_root_and_sum_product_identities():
    assert T.nthroot(90+34*S.sqrt(7), 3) == 3+S.sqrt(7)
    i = S.Symbol('i', integer=True)
    first, second = S.Sum(i, (i, 1, 3)), S.Sum(i, (i, 4, 6))
    assert T.sum_add(first, second, 1).doit() == 21
    assert T.product_mul(S.Product(i, (i, 1, 3)), S.Product(i, (i, 4, 6)), 1).doit() == 720
    expression, rhs = T.clear_coefficients(4*x+6, 10)
    assert S.solve(expression-rhs, x) == [S.Integer(1)]


def test_groebner_trig_identity_and_solver_helpers():
    expression = (S.sin(x)**4-S.cos(x)**4)/(S.sin(x)**2-S.cos(x)**2)
    assert S.trigsimp(expression, method='groebner') == 1
    assert V.linear_coeffs(2*x+3*g+4, x, g) == [2, 3, 4]
    assert V.domain_check(1/(x-1), x, 0)
    assert not V.domain_check(1/(x-1), x, 1)


def test_commensurate_angle_identity_from_wigner_rotation():
    beta=S.Symbol('beta',real=True)
    expression=-6*S.sin(beta/2)**8+12*S.sin(beta/2)**6-6*S.sin(beta/2)**4+3*S.sin(beta)**4/8
    assert S.trigsimp(expression)==0
    assert S.trigsimp(expression,method='groebner')==0
