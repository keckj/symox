"""Solver options are applied or explicitly rejected, never silently ignored."""
import pytest
import sympy as sp
import symox as s


def test_solve_set_and_dict_outputs():
    x, y = s.symbols("x y")
    for eqs, variables in [(x*x-1, [x]), ([x+y-3, x-y-1], [x, y]), ([x+y-1], [x, y]), ([x-1, x-2], [x])]:
        actual_variables, actual_rows = s.solve(eqs, variables, set=True)
        foreign = [e.to_sympy() for e in eqs] if isinstance(eqs, list) else eqs.to_sympy()
        expected_variables, expected_rows = sp.solve(foreign, [v.to_sympy() for v in variables], set=True)
        assert [v.to_sympy() for v in actual_variables] == expected_variables
        assert {tuple(v.to_sympy() for v in row) for row in actual_rows} == expected_rows
    assert s.solve(x-1, x, dict=True, set=True) == [{x: 1}]
    assert s.solve(x-1, x, check=True, simplify=True) == [1]


@pytest.mark.parametrize("flags", [dict(manual=True), dict(mispelled=True)])
def test_solve_flags_accepted_like_reference(flags):
    # SymPy 1.14 implements manual=True and ignores flags it never reads
    x = s.Symbol("x")
    assert s.solve(x-1, x, **flags) == [1]
    assert sp.solve(sp.Symbol("x")-1, sp.Symbol("x"), **flags) == [1]


@pytest.mark.parametrize("flags", [dict(simplify=False), dict(rational=False)])
def test_solve_reference_flags_are_applied(flags):
    x = s.Symbol("x")
    assert s.solve(x-1, x, **flags) == [1]
    assert [r.to_sympy() for r in s.solve(x**2-2*x-1, x, **flags)] == sp.solve(sp.Symbol("x")**2-2*sp.Symbol("x")-1, sp.Symbol("x"), **flags)


def test_nsolve_convergence_controls_and_real_bracketing():
    x, y = s.symbols("x y")
    root = s.nsolve(s.cos(x)-x, x, 1, tol=1e-13, maxsteps=30)
    assert abs(float(root)-0.7390851332151607) < 1e-13
    with pytest.raises(ValueError, match="Could not find root"):
        s.nsolve(x*x-2, x, 100, maxsteps=1)
    # Explicitly disabling verification permits inspecting an unfinished iterate.
    assert float(s.nsolve(x*x-2, x, 100, maxsteps=1, verify=False)) > 40
    bracket = s.nsolve(x*x-2, x, (0, 2), solver="bisect", tol=1e-13)
    assert abs(float(bracket)-2**.5) < 1e-12
    with pytest.raises(ValueError):
        s.nsolve(x*x+1, x, (-1, 1), solver="bisect")
    with pytest.raises(ValueError):
        s.nsolve(1/x, x, (-1, 1), solver="bisect")
    mapping = s.nsolve([x+y-3, x-y-1], [x, y], [0, 0], dict=True)
    assert len(mapping) == 1 and abs(float(mapping[0][x])-2) < 1e-14
    vector = s.nsolve([x+y-3, x-y-1], [x, y], [0, 0])
    assert isinstance(vector, s.Matrix) and vector.shape == (2, 1)


@pytest.mark.parametrize("flags", [dict(verbose=True), dict(solver="muller"), dict(solver="halley")])
def test_nsolve_unsupported_flags_are_explicit(flags):
    x = s.Symbol("x")
    with pytest.raises(NotImplementedError):
        s.nsolve(x*x-2, x, 1, **flags)


@pytest.mark.parametrize("flags,error,message", [
    (dict(solver="foo"), ValueError, "could not recognize solver"),
    (dict(solver="anderson"), ValueError, "expected interval of 2 points, got 1"),
    (dict(solver="mdnewton"), KeyError, "norm"),
    (dict(method="secant"), ValueError, 'keyword to use is "solver"'),
])
def test_nsolve_solver_errors_like_reference(flags, error, message):
    # mpmath's findroot registry: SymPy 1.14 raises these exact errors
    for module in (s, sp):
        x = module.Symbol("x")
        with pytest.raises(error, match=message):
            module.nsolve(x*x-2, x, 1, **flags)


def test_nsolve_ignores_unknown_keywords_like_reference():
    x, y = s.symbols("x y")
    assert str(s.nsolve(x*x-2, x, 1, unknown=True)) == "1.41421356237310"
    assert str(s.nsolve(x*x-2, x, 1, modules=["mpmath"])) == "1.41421356237310"
    # a system always uses MDNewton, whatever valid solver name is given
    root = s.nsolve([x+y-3, x-y-1], [x, y], [0, 0], solver="anderson")
    assert [float(v) for v in root] == [2.0, 1.0]


@pytest.mark.parametrize("flags", [
    dict(force=False), dict(force=True), dict(warn=True), dict(bivariate=False),
    dict(minimal=True), dict(numerical=False), dict(simplify=False),
    dict(force=False, bivariate=False),
])
def test_solve_flags_like_reference_quadratic(flags):
    x = s.Symbol("x")
    assert s.solve(x**2-1, x, **flags) == [-1, 1]


@pytest.mark.parametrize("text,flags", [
    ("x*exp(x) - 1", dict(bivariate=False)),
    ("x*exp(x) - 1", dict(force=False)),
    ("x**x - 2", dict(force=False)),
    ("exp(2*x) - 3*exp(x) + 2", dict(force=False)),
    ("3**x - 9", dict(force=False)),
    ("log(x**2) - 2", dict(bivariate=False, force=False)),
    ("sqrt(x) + y", dict(minimal=True)),
    ("sqrt(x + 1) - sqrt(x) - 1", dict(force=False)),
    ("sqrt(x)*sqrt(x + 1) - 2", dict(check=False)),
    ("x**y - 2", dict(dict=True)),
    ("x**y - 2", dict(set=True)),
])
def test_solve_flags_like_reference(text, flags):
    expected = sp.solve(sp.sympify(text), sp.Symbol("x"), **flags)
    result = s.solve(s.sympify(text), s.Symbol("x"), **flags)
    assert sp.sympify(str(result)) == sp.sympify(str(expected))


def test_solve_without_lambert_or_posify_is_not_proven_empty():
    # SymPy returns [] after its integer guesses fail, although x**x = 2
    # has the real root exp(LambertW(log(2))): symox reports unsupported
    x = s.Symbol("x")
    with pytest.raises(NotImplementedError):
        s.solve(x**x-2, x, bivariate=False, force=False)
    assert s.solve(x**x-4, x, bivariate=False, force=False) == [2]
    assert s.solve(x**x-2, x, force=False) == [s.exp(s.LambertW(s.log(2)))]
    # without both Lambert and posify, neither solves x*exp(x) = 1
    for module in (s, sp):
        y = module.Symbol("x")
        with pytest.raises(NotImplementedError):
            module.solve(y*module.exp(y)-1, y, bivariate=False, force=False)


@pytest.mark.parametrize("flags", [dict(tol=0), dict(tol=float("nan")), dict(maxsteps=0)])
def test_nsolve_invalid_controls(flags):
    x = s.Symbol("x")
    with pytest.raises(ValueError):
        s.nsolve(x*x-2, x, 1, **flags)


def test_ode_unknown_hints_fail_explicitly():
    x = s.Symbol("x")
    f = s.Function("f")(x)
    ode = s.diff(f, x)-f
    assert s.dsolve(ode, f, hint="default") == s.dsolve(ode, f)
    # SymPy's dsolve rejects an unknown hint with ValueError
    with pytest.raises(ValueError, match="Hint not recognized: nonsense"):
        s.dsolve(ode, f, hint="nonsense")
    assert s.dsolve(ode, f, simplify=False) == s.Eq(f, s.Symbol("C1")*s.exp(x))


def test_expand_and_simplify_options_are_applied_or_rejected():
    x, a, b = s.symbols("x a b")
    nonzero = s.Symbol("nonzero", nonzero=True)
    for expression in (nonzero**(a+b), s.exp(a+b)):
        assert s.expand(expression, power_exp=True).is_Mul
        assert s.expand(expression, power_exp=True).to_sympy() == sp.expand(expression.to_sympy(), power_exp=True)
        assert s.expand(expression, power_exp=False) == expression
    collision = s.Symbol("_sx_expand_held_0")
    for expression in (s.exp(x*(a+b)), s.exp(a)*s.exp(b), collision+s.exp(a+b)):
        assert s.expand(expression, power_exp=False).to_sympy() == sp.expand(expression.to_sympy(), power_exp=False)
    for operation in (s.expand, lambda expression, **kw: expression.expand(**kw)):
        with pytest.raises(NotImplementedError):
            operation(x, unknown=True)
        assert operation(x, deep=False) == x
    assert s.simplify(x, ratio=1) == x
    with pytest.raises((TypeError, NotImplementedError)):
        s.integrate(x, x, unknown=True)


@pytest.mark.parametrize("expression,point", [
    ("exp(x)/x**5", "0"), ("1/(x*(x+1))", "0"),
    ("1/(x*(x+1))", "-1"), ("sin(x)/x", "0"),
    ("1/sin(x)", "0"), ("x**2", "3"),
    ("exp(x)/(x-2)**3", "2"), ("1/(x*x+1)", "I"),
])
def test_residues_against_sympy(expression, point):
    x = s.Symbol("x")
    actual = s.residue(s.sympify(expression), x, s.sympify(point))
    expected = sp.residue(sp.sympify(expression), sp.Symbol("x"), sp.sympify(point))
    assert sp.simplify(actual.to_sympy()-expected) == 0


def test_residue_branch_points_are_not_assumed_poles():
    x = s.Symbol("x")
    with pytest.raises(NotImplementedError):
        s.residue(s.sqrt(x)/x, x, 0)
    with pytest.raises(ValueError):
        s.residue(x, x, x)
